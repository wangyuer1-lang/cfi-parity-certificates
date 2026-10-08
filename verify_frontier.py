"""Independent connectivity-frontier certificate. Python standard library only.

Reads explicit witness edge lists. Does not import any other checker.
Polynomials are coefficient lists in ascending powers of q.
"""
from collections import defaultdict
from itertools import combinations, permutations
from math import comb
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent

def canonical(labels):
    rename = {}
    return tuple(rename.setdefault(x, len(rename)) for x in labels)

def frontier_counts(n, edges):
    """Count edge subsets by completed components and number of selected edges."""
    first, last = {}, {}
    for i, edge in enumerate(edges):
        for v in edge:
            first.setdefault(v, i)
            last[v] = i
    active = []
    states = {((), n-len(first), 0): 1}
    peak = 1
    for i, (u, v) in enumerate(edges):
        added = sorted(x for x in {u, v} if first[x] == i)
        boundary = active + added
        iu, iv = boundary.index(u), boundary.index(v)
        keep = [j for j, x in enumerate(boundary) if last[x] > i]
        nxt = defaultdict(int)
        for (part, k, s), multiplicity in states.items():
            p = list(part)
            start = max(p, default=-1)+1
            p += list(range(start, start+len(added)))
            for selected in (0, 1):
                merged = p
                if selected and p[iu] != p[iv]:
                    old, new = p[iv], p[iu]
                    merged = [new if a == old else a for a in p]
                surviving = [merged[j] for j in keep]
                closed = len(set(merged)-set(surviving))
                key = (canonical(surviving), k+closed, s+selected)
                nxt[key] += multiplicity
        active = [boundary[j] for j in keep]
        states = nxt
        peak = max(peak, len(states))
    counts = [[0]*(len(edges)+1) for _ in range(n+1)]
    for (part, k, s), count in states.items():
        assert not part
        counts[k][s] += count
    assert all(sum(row[s] for row in counts) == comb(len(edges), s)
               for s in range(len(edges)+1))
    return counts, peak

def flow_coefficients(n, edges, counts):
    out = [0]*(len(edges)+1)
    for k, row in enumerate(counts):
        for s, value in enumerate(row):
            if value:
                exponent = s-n+k
                assert exponent >= 0
                out[exponent] += (-1)**(len(edges)-s)*value
    while len(out)>1 and out[-1]==0: out.pop()
    return out

def add(a, b, scale=1):
    out = [0]*max(len(a), len(b))
    for i, x in enumerate(a): out[i] += x
    for i, x in enumerate(b): out[i] += scale*x
    while len(out)>1 and out[-1]==0: out.pop()
    return out

def mul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i+j] += x*y
    while len(out)>1 and out[-1]==0: out.pop()
    return out

def properties(n, edges):
    adj = [set() for _ in range(n)]
    assert all(u != v for u,v in edges)
    assert len({tuple(sorted(e)) for e in edges}) == len(edges)
    for u,v in edges: adj[u].add(v); adj[v].add(u)
    assert all(len(ns)==3 for ns in adj)
    assert not any(adj[u] & adj[v] for u,v in edges)
    deletion_checks = 0
    for size in range(3):
        for deleted in combinations(range(n), size):
            remaining = set(range(n))-set(deleted)
            seen = {min(remaining)}
            todo = list(seen)
            for u in todo:
                for v in adj[u] & remaining-seen:
                    seen.add(v); todo.append(v)
            assert seen == remaining
            deletion_checks += 1
    squares = []
    for group in combinations(range(n), 4):
        if not all(len(adj[u] & set(group)) == 2 for u in group): continue
        for tail in permutations(group[1:]):
            cycle = (group[0],)+tail
            if cycle[1]>cycle[-1]: continue
            if all(cycle[(i+1)%4] in adj[cycle[i]] for i in range(4)):
                squares.append(cycle)
    assert len(squares)==2
    return adj, squares, deletion_checks

def modifications(n, edges, cycle, adj):
    chosen=set(cycle)
    outside=[next(iter(adj[v]-chosen)) for v in cycle]
    remain=sorted(set(range(n))-chosen)
    rename={v:i for i,v in enumerate(remain)}
    centre=len(remain)
    contraction=[]
    for u,v in edges:
        if u in chosen and v in chosen: continue
        contraction.append((centre if u in chosen else rename[u],
                            centre if v in chosen else rename[v]))
    crossing=[(rename[u],rename[v]) for u,v in edges if u not in chosen and v not in chosen]
    crossing += [(rename[outside[0]],rename[outside[2]]),
                 (rename[outside[1]],rename[outside[3]])]
    return [(centre+1, contraction),(centre,crossing)]

def main():
    # Independently check loops, parallel edges and isolated vertices.
    assert frontier_counts(2,[(0,0)])[0] == [[0,0],[0,0],[1,1]]
    assert flow_coefficients(2,[(0,1)]*2,frontier_counts(2,[(0,1)]*2)[0]) == [-1,1]
    inputs=json.loads((ROOT/'inputs'/'witness_graphs.json').read_text(encoding='utf-8'))
    results=[]
    A=mul(mul([-1,1],[-2,1]),[-3,1])
    P=[-352,458,-264,83,-14,1]
    R=mul([-3,1],[-37,30,-9,1])
    for j,graph in enumerate(inputs['graphs']):
        n=graph['vertices'];edges=graph['edges']
        adj,squares,num=properties(n,edges)
        counts,peak=frontier_counts(n,edges)
        f=flow_coefficients(n,edges,counts)
        rows=[]
        for index,cycle in enumerate(squares):
            polys=[]
            for nv,es in modifications(n,edges,cycle,adj):
                c,_=frontier_counts(nv,es)
                polys.append(flow_coefficients(nv,es,c))
            exceptional=j==1 and cycle==(3,9,8,13)
            assert polys[0]==mul(A,add(P,[1]) if exceptional else P)
            assert polys[1]==mul(A,add(R,[-5,2]) if exceptional else R)
            weight=add(add(mul([-3,-3,1],f),mul(mul([-4,1],[6,-4,1]),polys[0]),-1),mul([0,-4,1],polys[1]))
            rows.append({'cycle':cycle,'contracted':polys[0],'crossed':polys[1],'square_weight':weight})
        results.append({'graph':graph['name'],'counts':counts,'flow':f,'squares':rows,
                        'deletion_checks':num,'frontier_peak_states':peak})
        print(graph['name'], 'frontier certificate complete; peak states',peak,flush=True)
    assert results[0]['counts']==results[1]['counts']
    assert results[0]['counts'][1][15]==223080
    totals=[add(r['squares'][0]['square_weight'],r['squares'][1]['square_weight']) for r in results]
    difference=add(totals[0],totals[1],-1)
    predicted=[-1]
    for factor in ([-1,1],[-2,1],[-3,1],[-3,1],[-4,1],[2,1]):predicted=mul(predicted,factor)
    assert difference==predicted
    report={'algorithm':'connectivity frontier dynamic programming','witnesses':results,
            'square_difference':difference,'all_assertions_passed':True}
    (ROOT/'frontier_results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('PASS: full rank arrays and symbolic square polynomials, independently recomputed.')

if __name__=='__main__': main()
