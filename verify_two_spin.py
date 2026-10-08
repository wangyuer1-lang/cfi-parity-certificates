import json
from pathlib import Path
import networkx as nx
from itertools import combinations
from math import comb


def cfi_graph(B, twist_edges=frozenset()):
    G=nx.Graph()
    # canonical edge labels as sorted tuples
    edges=[tuple(sorted(e)) for e in B.edges()]
    inc={v:[] for v in B.nodes()}
    for e in edges:
        u,v=e; inc[u].append(e); inc[v].append(e)
    # outer nodes ('o',v,e,bit) bit 1=a,0=b
    for v in B.nodes():
        for e in inc[v]:
            for bit in (0,1):
                G.add_node(('o',v,e,bit))
        # middle nodes indexed by even subsets encoded frozenset of incident edges
        lst=inc[v]
        for r in range(len(lst)+1):
            if r%2: continue
            for S in combinations(lst,r):
                S=frozenset(S)
                m=('m',v,S)
                G.add_node(m)
                for e in lst:
                    # adjacent a(bit1) if e in S, else b(bit0)
                    bit=1 if e in S else 0
                    G.add_edge(m,('o',v,e,bit))
    for e in edges:
        u,v=e
        c=1 if e in twist_edges else 0
        if c==0:
            G.add_edge(('o',u,e,0),('o',v,e,0))
            G.add_edge(('o',u,e,1),('o',v,e,1))
        else:
            G.add_edge(('o',u,e,0),('o',v,e,1))
            G.add_edge(('o',u,e,1),('o',v,e,0))
    return G

class Factor:
    def __init__(self, scope, table):
        self.scope=tuple(scope)
        self.table=dict(table)

def multiply(factors):
    if not factors: return Factor((),{0:1})
    scope=[]
    for f in factors:
        for v in f.scope:
            if v not in scope: scope.append(v)
    pos={v:i for i,v in enumerate(scope)}
    out={}
    for mask in range(1<<len(scope)):
        val=1
        for f in factors:
            fm=0
            for j,v in enumerate(f.scope):
                if (mask>>pos[v])&1: fm|=1<<j
            val*=f.table.get(fm,0)
            if not val: break
        if val: out[mask]=val
    return Factor(scope,out)

def sum_out(f,v):
    idx=f.scope.index(v)
    scope=list(f.scope); scope.pop(idx)
    out={}
    for mask,val in f.table.items():
        low=mask & ((1<<idx)-1)
        high=mask>>(idx+1)
        nm=low | (high<<idx)
        out[nm]=out.get(nm,0)+val
    return Factor(scope,out)

def greedy_min_fill_order(G):
    H=G.copy(); order=[]
    while H:
        best=None; bestkey=None
        for v in H.nodes():
            nbr=list(H.neighbors(v))
            fill=0
            for i in range(len(nbr)):
                for j in range(i+1,len(nbr)):
                    if not H.has_edge(nbr[i],nbr[j]): fill+=1
            key=(fill,len(nbr))
            if bestkey is None or key<bestkey:
                bestkey=key; best=v
        nbr=list(H.neighbors(best))
        for i in range(len(nbr)):
            for j in range(i+1,len(nbr)):
                H.add_edge(nbr[i],nbr[j])
        H.remove_node(best); order.append(best)
    return order

def two_spin_Z(G,beta,gamma,lam):
    nodes=list(G.nodes())
    factors=[]
    for v in nodes:
        factors.append(Factor((v,),{0:1,1:lam}))
    for u,v in G.edges():
        factors.append(Factor((u,v),{0:beta,1:1,2:1,3:gamma}))
    order=greedy_min_fill_order(G)
    for v in order:
        hit=[f for f in factors if v in f.scope]
        factors=[f for f in factors if v not in f.scope]
        prod=multiply(hit)
        factors.append(sum_out(prod,v))
    prod=multiply(factors)
    assert prod.scope==()
    return prod.table.get(0,0)

def cubic_formula(n,beta,gamma,lam):
    return (2**(n//2+1) * lam**(4*n) * (beta*gamma-1)**(9*n//2) * (beta**3-lam**2*gamma**3)**n)

if __name__=='__main__':
    B=nx.complete_graph(4)
    e0=tuple(sorted(next(iter(B.edges()))))
    G0=cfi_graph(B,set())
    G1=cfi_graph(B,{e0})
    print('K4 CFI |V|,|E|',G0.number_of_nodes(),G0.number_of_edges())
    beta,gamma,lam=2,3,5
    z0=two_spin_Z(G0,beta,gamma,lam)
    z1=two_spin_Z(G1,beta,gamma,lam)
    diff=z0-z1
    rhs=cubic_formula(4,beta,gamma,lam)
    print('Z0-Z1',diff)
    print('formula',rhs)
    assert diff==rhs
    first={'base':'K4','computed_difference':diff,'formula':rhs}
    print('match',diff==rhs)

    # triangle base: direct full subset enumeration vs general formula at numeric s,t
    B=nx.cycle_graph(3)
    e0=tuple(sorted(next(iter(B.edges()))))
    G0=cfi_graph(B,set()); G1=cfi_graph(B,{e0})
    s,t=2,3
    def Q_brute(G,s,t):
        nd=list(G.nodes()); N=len(nd); ans=0
        edgeidx=[(nd.index(u),nd.index(v)) for u,v in G.edges()]
        for mask in range(1<<N):
            k=mask.bit_count(); ee=0
            for i,j in edgeidx:
                if (mask>>i)&1 and (mask>>j)&1: ee+=1
            ans += s**k*t**ee
        return ans
    qdiff=Q_brute(G0,s,t)-Q_brute(G1,s,t)
    n=B.number_of_nodes(); m=B.number_of_edges()
    def D(d):
        A=1; C=1
        for j in range(d+1):
            fac=(1+s*(t**(d-j)))**comb(d,j)
            if j%2==0:A*=fac
            else:C*=fac
        return A-C
    rhs=(2**(m-n+1))*s**(2*m)*(t-1)**m
    for v in B.nodes(): rhs*=D(B.degree(v))
    print('C3 CFI |V|',G0.number_of_nodes())
    print('Qdiff brute',qdiff)
    print('Qdiff formula',rhs)
    assert qdiff==rhs
    result={'checks':[first,{'base':'C3','computed_difference':qdiff,'formula':rhs}],'all_assertions_passed':True}
    Path(__file__).with_name('two_spin_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print('match',qdiff==rhs)
