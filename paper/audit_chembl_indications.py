import requests, collections, sys
B="https://www.ebi.ac.uk/chembl/api/data/drug_indication.json"
def pull(dz):
    rows=[];off=0
    while True:
        r=requests.get(B,params={"mesh_heading__icontains":dz,"limit":1000,"offset":off},timeout=60).json()
        rows+=r["drug_indications"]
        if not r["page_meta"]["next"]: break
        off+=1000
    return rows
for dz in sys.argv[1:]:
    rows=pull(dz)
    mols={}
    for x in rows:
        m=x["molecule_chembl_id"]; ph=float(x["max_phase_for_ind"] or 0)
        mols.setdefault(m,[]).append((x["mesh_heading"],ph))
    heads=collections.Counter(x["mesh_heading"] for x in rows)
    best={m:max(p for _,p in v) for m,v in mols.items()}
    ph=collections.Counter(best.values())
    print(f"\n=== {dz}: {len(rows)} indication rows, {len(mols)} unique molecules")
    print(" max phase FOR THIS INDICATION:",dict(sorted(ph.items())))
    print(" approved (phase 4) for indication:",sum(1 for v in best.values() if v==4))
    print(" mesh headings matched:",heads.most_common(12))
