"""Formal full cubic tensors and direct exterior traces; exact arithmetic."""
from itertools import combinations, product
from pathlib import Path
import json
import random
import sympy as s

q,v=s.symbols('q v')

def canon(values):
    labels={}
    return tuple(labels.setdefault(x,len(labels)) for x in values)

def c(pair,x): return int(pair[0]==x)-int(pair[1]==x)

def oriented_delta(a,b):
    return int(a==b)-int(a==b[::-1])

def local_check():
    pairs=list(combinations(range(6),2))
    patterns=sorted({canon(sum(pp,())) for pp in product(pairs,repeat=3)})
    for word in patterns:
        pp=[tuple(word[2*i:2*i+2]) for i in range(3)]
        L=[s.Integer(1),s.Integer(1)]
        for bits in product(range(2),repeat=3):
            z=[pp[i][0] if bits[i] else pp[i][1] for i in range(3)]
            f=q+3*v+v*v*sum(z[i]==z[j] for i,j in combinations(range(3),2))+v**3*int(z[0]==z[1]==z[2])
            L[sum(bits)%2]*=f
        colors=range(max(word)+1)
        T0=sum(c(pp[0],x)*c(pp[1],x)*c(pp[2],x) for x in colors)
        T1=0
        for i,j in combinations(range(3),2):
            k=3-i-j
            T1+=oriented_delta(pp[i],pp[j])*sum(c(pp[i],x)**2*c(pp[k],x) for x in colors)
        rhs=-v**3*(q*q+3*q*v-v**3)*((q+3*v)*T0+v*v*T1)
        assert s.expand(L[0]-L[1]-rhs)==0,word
    return len(patterns)

def lift(M):
    d=M.rows;pairs=list(combinations(range(d),2))
    T=s.zeros(len(pairs))
    for i,(a,b) in enumerate(pairs):
        for j,(c,d_) in enumerate(pairs):
            T[i,j]=M[a,c]*int(b==d_)-M[b,c]*int(a==d_)+int(a==c)*M[b,d_]-int(b==c)*M[a,d_]
    return T

def trace_check():
    rng=random.Random(8102026);count=0
    for dim in range(2,8):
        for _ in range(4):
            Ms=[]
            for i in range(4):
                R=s.Matrix(dim,dim,lambda a,b:rng.randrange(-3,4))
                M=R+R.T
                M[dim-1,dim-1]-=s.trace(M)
                assert M==M.T and s.trace(M)==0
                Ms.append(M)
            def tr(word):
                out=s.eye(dim)
                for i in word:out=out*Ms[i]
                return s.trace(out)
            lhs=s.trace(lift(Ms[0])*lift(Ms[1])*lift(Ms[2])*lift(Ms[3]))
            rhs=(dim-4)*tr([0,1,2,3])-2*tr([0,1,3,2])-2*tr([0,2,1,3])
            rhs+=tr([0,1])*tr([2,3])+tr([0,2])*tr([1,3])+tr([0,3])*tr([1,2])
            assert lhs==rhs,(dim,lhs,rhs)
            count+=1
    return count

if __name__=='__main__':
    report={'full_formal_cubic_patterns':local_check(),
            'direct_general_trace_checks':trace_check(),
            'all_assertions_passed':True}
    Path(__file__).with_name('cubic_local_results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(report)
