"""Independent deletion-contraction and flow-value checks, exact arithmetic."""
from functools import lru_cache
from collections import Counter
from pathlib import Path
import itertools as it
import importlib.util
import json
import math
import sympy as s

ROOT = Path(__file__).resolve().parent
q = s.symbols('q')

def norm(edges):
    vs = sorted({v for e in edges for v in e})
    ix = {v:i for i,v in enumerate(vs)}
    return tuple(sorted(tuple(sorted((ix[u], ix[v]))) for u,v in edges))

def plus(a,b,sgn=1):
    out = [0]*max(len(a),len(b))
    for i,c in enumerate(a):out[i]+=c
    for i,c in enumerate(b):out[i]+=sgn*c
    while len(out)>1 and out[-1]==0:out.pop()
    return tuple(out)

def times(a,b):
    out = [0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    while len(out)>1 and out[-1]==0:out.pop()
    return tuple(out)

def expr(p):return s.expand(sum(c*q**i for i,c in enumerate(p)))
def ev(p,x):return sum(c*x**i for i,c in enumerate(p))

@lru_cache(None)
def flow(es):
    if not es:return (1,)
    loops = sum(u==v for u,v in es)
    if loops:
        p=flow(norm([(u,v) for u,v in es if u!=v]))
        for _ in range(loops):p=times(p,(-1,1))
        return p
    adj={}
    for j,(u,v) in enumerate(es):
        adj.setdefault(u,[]).append((v,j))
        adj.setdefault(v,[]).append((u,j))
    if any(len(ns)==1 for ns in adj.values()):return (0,)
    for v,ns in adj.items():
        if len(ns)==2:
            (a,j),(b,k)=ns
            return flow(norm([e for i,e in enumerate(es) if i not in (j,k)]+[(a,b)]))
    # Split components, and reject a bridge, before deletion-contraction.
    visited=set();components=[];bridge=False;low={};disc={};clock=0
    def dfs(v,parent_edge,comp):
        nonlocal clock,bridge
        visited.add(v);comp.add(v);disc[v]=low[v]=clock;clock+=1
        for w,j in adj[v]:
            if j==parent_edge:continue
            if w not in visited:
                dfs(w,j,comp);low[v]=min(low[v],low[w])
                if low[w]>disc[v]:bridge=True
            else:low[v]=min(low[v],disc[w])
    for v in adj:
        if v not in visited:
            comp=set();dfs(v,-1,comp);components.append(comp)
    if bridge:return (0,)
    if len(components)>1:
        out=(1,)
        for comp in components:out=times(out,flow(norm([e for e in es if e[0] in comp])))
        return out
    a,b=es[0];rest=es[1:]
    contracted=norm([(a if u==b else u,a if v==b else v) for u,v in rest])
    return plus(flow(contracted),flow(norm(rest)),-1)

def partitions(xs):
    if not xs:yield ();return
    x,*tail=xs
    for p in partitions(tail):
        yield ((x,),)+p
        for i in range(len(p)):
            yield tuple(((x,)+b if i==j else b) for j,b in enumerate(p))

@lru_cache(None)
def allowed(d):return tuple(p for p in partitions(tuple(range(d))) if all(len(b)>=2 for b in p))
def mu(p):return (-1)**(len(p)-1)*math.factorial(len(p)-1)

def cumulant_flow(edges):
    vertices=sorted({v for e in edges for v in e})
    incid={v:[j for j,e in enumerate(edges) if v in e] for v in vertices}
    ans=(0,);count=0
    for choice in it.product(*(allowed(len(incid[v])) for v in vertices)):
        image={};coeff=1;new=0
        for v,p in zip(vertices,choice):
            coeff*=mu(p)
            for block in p:
                for k in block:image[v,incid[v][k]]=new
                new+=1
        es=norm([(image[u,j],image[v,j]) for j,(u,v) in enumerate(edges)])
        ans=plus(ans,tuple(coeff*x for x in flow(es)));count+=1
    return ans,count

def weighted_flows(edges,order,klein=False):
    vertices=sorted({v for e in edges for v in e})
    inc={v:[] for v in vertices}
    for j,(u,v) in enumerate(edges):inc[u].append((j,1));inc[v].append((j,-1))
    def total(xs):
        if not klein:return sum(xs)%order
        z=0
        for x in xs:z^=x
        return z
    def neg(x):return x if klein else (-x)%order
    ans=0;count=0
    for xs in it.product(range(1,order),repeat=len(edges)):
        outgoing={v:[xs[j] if sign==1 else neg(xs[j]) for j,sign in inc[v]] for v in vertices}
        if any(total(vals) for vals in outgoing.values()):continue
        count+=1;weight=1
        for vals in outgoing.values():
            weight*=sum(mu(p) for p in allowed(len(vals)) if all(total([vals[k] for k in block])==0 for block in p))
        ans+=weight
    return ans,count

def main():
    assert flow(norm([(0,1),(1,2),(2,0)]))==(-1,1)
    assert expr(flow(norm([(0,1)]*3)))==s.expand((q-1)*(q-2))
    assert expr(flow(norm(list(it.combinations(range(4),2)))))==s.expand((q-1)*(q-2)*(q-3))
    spec=importlib.util.spec_from_file_location('cert',ROOT/'verify_certificates.py')
    cert=importlib.util.module_from_spec(spec);spec.loader.exec_module(cert)
    for edges in [cert.E0,cert.E1]:
        eset=set(edges)
        assert not any({(a,b),(a,c),(b,c)}<=eset for a,b,c in it.combinations(range(16),3))
        assert len(cert.squares(edges))==2
    bases=[expr(flow(norm(e))) for e in [cert.E0,cert.E1]]
    assert bases[0]==bases[1]
    values=[]
    for es in [cert.E0,cert.E1]:
        rows=[]
        for cyc in cert.squares(es):rows.append([expr(flow(norm(e))) for n,e in cert.modified_graphs(es,cyc)])
        values.append(rows)
    dc=s.expand(sum(r[0] for r in values[0])-sum(r[0] for r in values[1]))
    dx=s.expand(sum(r[1] for r in values[0])-sum(r[1] for r in values[1]))
    diff=s.factor((q-4)*((-q*q+4*q-6)*dc+q*dx))
    assert s.expand(diff+(q-1)*(q-2)*(q-3)**2*(q-4)*(q+2))==0
    result={'all_q_independent':{'common_flow':str(s.factor(bases[0])),'contracted_difference':str(s.factor(dc)),
            'crossed_difference':str(s.factor(dx)),'z4_difference':str(diff)}}
    print('Deletion-contraction all-q certificate:',result['all_q_independent'],flush=True)
    star=[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0),(0,5),(5,6),(6,0)]
    chain=[(0,1),(1,2),(2,0),(2,3),(3,4),(4,2),(4,5),(5,6),(6,4)]
    graphs={**{'theta'+str(d):[(0,1)]*d for d in range(2,7)},
            'doubled_triangle':[(0,1)]*2+[(1,2)]*2+[(2,0)]*2,
            'K4_plus_parallel':list(it.combinations(range(4),2))+[(0,1)],
            'three_triangles_star':star,'three_triangles_chain':chain}
    result['cumulant_graphs']={}
    for name,es in graphs.items():
        polynomial,split_count=cumulant_flow(es)
        rows={}
        for order,klein in [(2,False),(3,False),(4,False),(4,True)]:
            value,count=weighted_flows(es,order,klein)
            assert value==ev(polynomial,order),(name,order,klein,value,polynomial)
            assert count==ev(flow(norm(es)),order)
            rows[('V4' if klein else 'Z'+str(order))]={'weighted':value,'ordinary':count}
        result['cumulant_graphs'][name]={'polynomial':str(s.factor(expr(polynomial))),
                                      'vertex_splittings':split_count,'groups':rows}
        print(name,result['cumulant_graphs'][name],flush=True)
    assert result['cumulant_graphs']['theta5']['groups']['Z3']['weighted']==90
    assert result['cumulant_graphs']['doubled_triangle']['groups']['Z3']['weighted']==-10
    assert result['cumulant_graphs']['three_triangles_star']['groups']['Z2']['weighted']==16
    assert result['cumulant_graphs']['three_triangles_chain']['groups']['Z2']['weighted']==4
    result['flow_coefficients']=[{'base':list(flow(norm(es))),
        'modifications':[{'cycle':list(cyc),'contracted':list(flow(norm(cert.modified_graphs(es,cyc)[0][1]))),
                          'crossed':list(flow(norm(cert.modified_graphs(es,cyc)[1][1])))}
                         for cyc in cert.squares(es)]} for es in [cert.E0,cert.E1]]
    (ROOT/'independent_audit_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print('PASS: independent flow recursion, vertex partitions, and 36 finite-group flow checks.',flush=True)

if __name__=='__main__':main()
