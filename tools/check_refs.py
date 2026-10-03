"""Check every BibTeX entry against Crossref: DOI resolves, title, year, volume, first author.

Usage: python tools/check_refs.py paper/refs.bib [paper/main.tex]
"""
import difflib
import re
import sys

import requests


def fields(body):
    out, i = {}, 0
    for m in re.finditer(r'(\w+)\s*=\s*\{', body):
        if m.start() < i:
            continue
        j, depth = m.end(), 1
        while depth:
            depth += {'{': 1, '}': -1}.get(body[j], 0)
            j += 1
        out[m.group(1).lower()] = body[m.end():j - 1].strip()
        i = j
    return out


def entries(path):
    s = open(path).read()
    out = {}
    for m in re.finditer(r'@(\w+)\{([^,\s]+),', s):
        j, depth = s.index('{', m.start()) + 1, 1
        while depth:
            depth += {'{': 1, '}': -1}.get(s[j], 0)
            j += 1
        out[m.group(2)] = (m.group(1), fields(s[m.end():j - 1]))
    return out


def norm(t):
    return re.sub(r'\s+', ' ', re.sub(r'[^a-z0-9 ]', ' ', re.sub(r'[{}\\]|\$[^$]*\$', '', t.lower()))).strip()


def main(bibpath, texpath=None):
    bib = entries(bibpath)
    if texpath:
        tex = open(texpath).read()
        cited = {k.strip() for g in re.findall(r'\\cite[pt]?\{([^}]*)\}', tex) for k in g.split(',')}
        print('cited but missing:', sorted(cited - set(bib)))
        print('never cited:', sorted(set(bib) - cited))
    problems = 0
    for key, (typ, f) in sorted(bib.items()):
        doi = f.get('doi')
        if not doi:
            print(f'NO-DOI   {key} ({typ}): {f.get("title", "")[:70]}')
            continue
        r = requests.get(f'https://api.crossref.org/works/{doi}', timeout=30)
        if r.status_code != 200:
            print(f'BAD-DOI  {key}: {doi}'); problems += 1; continue
        m = r.json()['message']
        issues = []
        ct = (m.get('title') or [''])[0]
        sim = difflib.SequenceMatcher(None, norm(ct), norm(f.get('title', ''))).ratio()
        if sim < 0.9:
            issues.append(f'title {sim:.2f} crossref="{ct[:90]}"')
        years = {p[0] for k in ('issued', 'published-print', 'published-online') for p in m.get(k, {}).get('date-parts', []) if p and p[0]}
        try:
            by = int(f.get('year', '0'))
        except ValueError:
            by = 0
        if by not in years:
            issues.append(f'year {by} not in crossref {sorted(years)}')
        if m.get('volume') and f.get('volume') and m['volume'] != f['volume']:
            issues.append(f'volume {f["volume"]} vs {m["volume"]}')
        fam = [a.get('family', '') for a in m.get('author', [])]
        if fam and norm(fam[0]) and norm(fam[0]).split()[-1] not in norm(f.get('author', '')):
            issues.append(f'first author {fam[0]!r} missing')
        if issues:
            problems += 1
            print(f'CHECK    {key}: ' + '; '.join(issues))
    print(f'{len(bib)} entries, {problems} with issues')


if __name__ == '__main__':
    main(*sys.argv[1:])
