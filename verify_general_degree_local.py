"""Exact complete local tensors, without projecting to an edge channel."""
from itertools import combinations,product
from fractions import Fraction
from pathlib import Path
import json
from verify_extensions import allowed,mu

def mul(a,b,D):
    out=[0]*(D+1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b[:D+1-i]):out[i+j]+=x*y
    return out

def tensor_series(u,R,h,d):
    q=len(u);pairs=list(combinations(range(q),2));D=d;cache={}
    for colors in product(range(q),repeat=d):
        f=[0]*(D+1)
        for x in range(q):
            p=[h[x]]+[0]*D
            for z in colors:p=mul(p,[u[x]*u[z],R[x][z]],D)
            f=[a+b for a,b in zip(f,p)]
        cache[colors]=f
    for indices in product(range(len(pairs)),repeat=d):
        pp=[pairs[i] for i in indices];L=[[1]+[0]*D,[1]+[0]*D]
        for bits in product(range(2),repeat=d):
            colors=tuple(p[0] if bit else p[1] for p,bit in zip(pp,bits))
            parity=sum(bits)%2
            L[parity]=mul(L[parity],cache[colors],D)
        outer=1
        for a,b in pp:outer*=h[a]*h[b]
        yield pp,[outer*(a-b) for a,b in zip(*L)]

def prod(xs):
    ans=1
    for x in xs:ans*=x
    return ans

def potts_local(d,q):
    # Independent full-tensor identity suggested by the proof repair:
    # [v^d]H_d = (-1)^d q^(2^(d-1)) kappa(Y_1,...,Y_d),
    # Y_i(x)=delta(x,a_i)-delta(x,b_i), x uniform on [q].
    R=[[int(x==y) for y in range(q)] for x in range(q)]
    count=0;nonzero=0
    for pp,coeff in tensor_series([1]*q,R,[1]*q,d):
        assert not any(coeff[:d]),(q,d,pp,coeff)
        Y=[[int(x==a)-int(x==b) for a,b in pp] for x in range(q)]
        kap=Fraction(0)
        for p in allowed(d):
            val=Fraction(mu(p))
            for b in p:val*=Fraction(sum(prod(Y[x][i] for i in b) for x in range(q)),q)
            kap+=val
        expected=(-1)**d*q**(2**(d-1))*kap
        assert coeff[d]==expected,(q,d,pp,coeff,expected)
        count+=1;nonzero+=coeff[d]!=0
    return {'degree':d,'q':q,'complete_tensor_coordinates':count,'nonzero_leading_coordinates':nonzero}

def main():
    results={'potts_full_local':[],'rank_one_full_local':[]}
    for d,q in [(2,2),(3,2),(4,2),(5,2),(6,2),(7,2),(8,2),
                (2,3),(3,3),(4,3),(5,3),(6,3),(2,4),(3,4),(4,4)]:
        r=potts_local(d,q);results['potts_full_local'].append(r);print(r,flush=True)
    R=[[2,-1,3],[-1,4,2],[3,2,-2]]
    for d in [2,3,4,5,6]:
        count=0
        for pp,coeff in tensor_series([1,2,3],R,[1,2,3],d):
            assert not any(coeff[:d]);count+=1
        results['rank_one_full_local'].append({'degree':d,'coordinates':count})
    # Why the minimum degree condition cannot be silently removed.
    pp,coeff=next(tensor_series([1,2],[[1,0],[0,1]],[1,1],1))
    assert coeff[0]!=0
    results['leaf_local_constant']={'pair':pp,'coefficients':coeff}
    Path(__file__).with_name('full_local_audit_results.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    print('PASS: full local tensors and the stronger cumulant coefficient, exact Python integers.',flush=True)

if __name__=='__main__':main()
