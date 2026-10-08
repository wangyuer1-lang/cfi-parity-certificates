"""Exact certificate for the 16-vertex CFI parity witness; Python standard library only.
Run: python verify_certificates.py
No external input files are required. Outputs are written beside the script.
"""
import itertools as it
import json
import math
from pathlib import Path

ADJ=[[6,10,12],[7,8,11],[5,7,9],[9,13,14],[5,11,15],[2,4,12],
     [0,7,14],[1,2,6],[1,9,13],[2,3,8],[0,11,15],[1,4,10],
     [0,5,13],[3,8,12],[3,6,15],[4,10,14]]
E0=sorted((u,v) for u,ns in enumerate(ADJ) for v in ns if u<v)
E1=sorted([e for e in E0 if e not in [(4,11),(10,15)]]+[(4,10),(11,15)])

def rank_counts(n,edges):
    m=len(edges);parent=list(range(n));size=[1]*n
    counts=[[0]*(m+1) for _ in range(n+1)]
    def visit(i,k,s):
        if i==m:
            counts[k][s]+=1
            return
        visit(i+1,k,s)
        a,b=edges[i]
        while parent[a]!=a:a=parent[a]
        while parent[b]!=b:b=parent[b]
        if a==b:
            visit(i+1,k,s+1)
        else:
            if size[a]<size[b]:a,b=b,a
            parent[b]=a;size[a]+=size[b]
            visit(i+1,k-1,s+1)
            size[a]-=size[b];parent[b]=b
    visit(0,n,0)
    assert all(sum(row[s] for row in counts)==math.comb(m,s) for s in range(m+1))
    return counts

def flow_from_counts(n,m,counts,q=5):
    return sum(c*(-1)**(m-s)*q**(s-n+k)
               for k,row in enumerate(counts) for s,c in enumerate(row) if c)

def flow_coefficients(n,m,counts):
    out=[0]*(m+1)
    for k,row in enumerate(counts):
        for s,c in enumerate(row):
            if c:
                assert s-n+k>=0
                out[s-n+k]+=c*(-1)**(m-s)
    while len(out)>1 and out[-1]==0:out.pop()
    return out

def connected_after_deleting(n,edges,removed):
    remaining=set(range(n))-set(removed)
    seen={next(iter(remaining))};todo=list(seen)
    while todo:
        v=todo.pop()
        for a,b in edges:
            w=b if a==v else a if b==v else None
            if w in remaining and w not in seen:seen.add(w);todo.append(w)
    return seen==remaining

def squares(edges):
    edge_set=set(edges);cycles=set()
    for nodes in it.combinations(range(16),4):
        if sum((a,b) in edge_set for a,b in it.combinations(nodes,2))!=4:continue
        for tail in it.permutations(nodes[1:]):
            c=(nodes[0],)+tail
            if all(tuple(sorted((c[i],c[(i+1)%4]))) in edge_set for i in range(4)):
                cycles.add(min(c,(c[0],)+tuple(reversed(c[1:]))))
    return sorted(cycles)

def modified_graphs(edges,cycle):
    selected=set(cycle);others=sorted(set(range(16))-selected)
    outside=[]
    for v in cycle:
        neighbors=[b if a==v else a for a,b in edges if a==v or b==v]
        outside.append(next(x for x in neighbors if x not in selected))
    ix={v:i for i,v in enumerate(others)}
    contracted=[]
    for a,b in edges:
        if a in selected and b in selected:continue
        contracted.append((12 if a in selected else ix[a],12 if b in selected else ix[b]))
    crossed=[(ix[a],ix[b]) for a,b in edges if a not in selected and b not in selected]
    crossed.extend([(ix[outside[0]],ix[outside[2]]),(ix[outside[1]],ix[outside[3]])])
    return (13,contracted),(12,crossed)

def main():
    tables=[];results=[]
    for j,edges in enumerate([E0,E1]):
        assert len(edges)==len(set(edges))==24
        assert all(sum(v in e for e in edges)==3 for v in range(16))
        assert all(connected_after_deleting(16,edges,r)
                   for k in range(3) for r in it.combinations(range(16),k))
        table=rank_counts(16,edges);tables.append(table)
        F=flow_from_counts(16,24,table)
        assert F==3336 and table[1][15]==223080
        rows=[]
        for cycle in squares(edges):
            values=[];coefficients=[]
            for n,es in modified_graphs(edges,cycle):
                modified_counts=rank_counts(n,es)
                values.append(flow_from_counts(n,len(es),modified_counts))
                coefficients.append(flow_coefficients(n,len(es),modified_counts))
            fc,fx=values
            rows.append({'cycle':cycle,'contracted':fc,'crossed':fx,'weight':7*F-11*fc+5*fx,'contracted_coefficients':coefficients[0],'crossed_coefficients':coefficients[1]})
        total=sum(row['weight'] for row in rows)
        assert len(rows)==2 and total==[6480,6816][j]
        results.append({'graph':j,'F5':F,'flow_coefficients':flow_coefficients(16,24,table),'squares':rows,'z4':total})
        print('B'+str(j),results[-1],flush=True)
    assert tables[0]==tables[1]
    out=Path(__file__).resolve().parent
    (out/'rank_counts.json').write_text(json.dumps(tables[0]),encoding='utf-8')
    report={'subsets_each':sum(map(sum,tables[0])),'rank_tables_equal':True,
            'flow_values':{q:flow_from_counts(16,24,tables[0],q) for q in [3,4,5,6]},
            'witnesses':results}
    (out/'certificate_results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    assert (6480-6816)*2**9==-172032
    print('PASS: equal full Tutte polynomials; unequal square coefficient; all arithmetic exact.')

if __name__=='__main__':main()
