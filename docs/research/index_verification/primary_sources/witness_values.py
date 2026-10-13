# Witness values for competing literature conventions (03b_PRIMARY_SOURCE_DEFINITIONS.md).
# Pure Python, no package import. Run: python3 -I witness_values.py
import math, itertools
from collections import deque
def G(edges):
    n=max(max(e) for e in edges)+1; adj=[set() for _ in range(n)]
    for u,v in edges: adj[u].add(v); adj[v].add(u)
    return n,edges,adj
def bfs(adj,s):
    d=[None]*len(adj); d[s]=0; q=deque([s])
    while q:
        x=q.popleft()
        for y in adj[x]:
            if d[y] is None: d[y]=d[x]+1; q.append(y)
    return d
graphs={'P4':[(0,1),(1,2),(2,3)],'K13':[(0,1),(0,2),(0,3)],'C4':[(0,1),(1,2),(2,3),(3,0)],'C5':[(0,1),(1,2),(2,3),(3,4),(4,0)]}
for nm,E in graphs.items():
    n,E,adj=G(E); deg=[len(a) for a in adj]; D=[bfs(adj,s) for s in range(n)]; m=len(E)
    W=sum(D[u][v] for u,v in itertools.combinations(range(n),2))
    abc=sum(math.sqrt((deg[u]+deg[v]-2)/(deg[u]*deg[v])) for u,v in E)
    # PI edge (min convention), PI strict (both endpoints closer), PI_v, Sz
    PI=PIs=PIv=Sz=0
    for u,v in E:
        nu=sum(1 for x in range(n) if D[u][x]<D[v][x]); nv=sum(1 for x in range(n) if D[v][x]<D[u][x])
        PIv+=nu+nv; Sz+=nu*nv
        mu=mv=su=sv=0
        for x,y in E:
            a=min(D[u][x],D[u][y]); b=min(D[v][x],D[v][y])
            mu+=a<b; mv+=b<a
            su+= (D[u][x]<D[v][x] and D[u][y]<D[v][y]); sv+=(D[v][x]<D[u][x] and D[v][y]<D[u][y])
        PI+=mu+mv; PIs+=su+sv
    # edge wiener variants
    We0=We1=We4=0
    for (a,b),(c,d) in itertools.combinations(E,2):
        dd=[D[a][c],D[a][d],D[b][c],D[b][d]]
        We1+=min(dd); We4+=max(dd)
        shared = len({a,b}&{c,d})>0
        We0+= 1 if shared else min(dd)+1
    Hb=-sum(d/(2*m)*math.log2(d/(2*m)) for d in deg); He=Hb*math.log(2)
    print(nm,'W',W,'Sz',Sz,'ABC',round(abc,4),'ABC*sqrt2',round(abc*math.sqrt(2),4),'PI(min)',PI,'PI(both-endpoints)',PIs,'PI_v',PIv,'We0=W(L)',We0,'We1',We1,'We4',We4,'Ideg bits',round(Hb,4),'nats',round(He,4))
# R_k parentheses witness on K13 (Delta=3), k=2
Delta=3;k=2
def mk(d,paren=True): return (Delta+k-d) if k<=d else ((Delta+k-d)%Delta if paren else Delta-d+(k%Delta))
deg=[3,1,1,1];E=graphs['K13']
print('K13 R2_M1 paren',sum(mk(deg[u])+mk(deg[v]) for u,v in E),'noparen',sum(mk(deg[u],False)+mk(deg[v],False) for u,v in E))
E=graphs['P4'];deg=[1,2,2,1];Delta=2
for k in (1,2):
    def mk2(d): return (Delta+k-d) if k<=d else (Delta+k-d)%Delta
    print('P4 k',k,[mk2(d) for d in deg])
