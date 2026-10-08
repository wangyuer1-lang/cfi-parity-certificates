import itertools as it
import json
from pathlib import Path
import sympy as s

q=s.symbols('q')
def phi(w):
    if not w: return q-1
    value=(-1)**len(w)*(q-1)
    for mask in range(1,1<<len(w)):
        a=[w[i] for i in range(len(w)) if mask>>i&1]
        value+=(-1)**(len(w)-len(a))*s.prod(q*int(x==y)-1 for x,y in zip(a,a[1:]+a[:1]))
    return s.expand(value)

def exterior(w):
    value=0
    for mask in range(1<<len(w)):
        a=[w[i] for i in range(len(w)) if mask>>i&1]
        b=[w[i] for i in range(len(w)) if not(mask>>i&1)]
        value+=phi(a)*phi(b)-phi(a+b)
    return s.expand(value/2)

def canonical(w):
    labs={};out=[]
    for x in w:
        if x not in labs:labs[x]=len(labs)+1
        out.append(labs[x])
    return tuple(out)

patterns=sorted({canonical(w) for w in it.product(range(4),repeat=4)})
rows=[]
for w in patterns:
    g=lambda i,j:q*int(w[i]==w[j])-1
    d4=q-len(set(w))+sum(s.prod(q*int(a==x)-1 for x in w) for a in set(w))
    adj=g(0,1)*g(2,3)+g(0,3)*g(1,2)
    opposite=g(0,2)*g(1,3)
    rhs=(q-9)*(q-4)*d4/q+(q*q-3*q-3)*adj+q*(q-4)*opposite
    lhs=exterior(w)
    assert s.expand(lhs-rhs)==0,w
    assert s.expand(phi(w)-(q-4)*d4/q-adj)==0,w
    rows.append({'pattern':''.join(map(str,w)),'trace':str(s.factor(lhs))})
representatives=['1111','1112','1122','1212','1123','1213','1234']
lines=[]
for row in rows:
    if row['pattern'] in representatives:
        value=s.sympify(row['trace'])
        lines.append(f"${row['pattern']}$ & ${s.latex(value)}$ \\\\")
out=Path(__file__).parent
(out/'trace_results.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
(out/'trace_table.tex').write_text('\n'.join(lines),encoding='utf-8')
print('All 15 equality patterns verified as identities in Q[q].')
print('\n'.join(lines))
