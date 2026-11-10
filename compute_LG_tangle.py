from sympy import *
from sympy.polys.matrices import DomainMatrix
from math import comb
from itertools import product
import re as re
import os
import pickle
import time

# =========================
# FOLDERS/PKLs
# =========================
cwd = os.getcwd()
pkl_folder = "tangle_pkls"
pkl_dir = os.path.join(cwd, pkl_folder)

def create_fold(fpath):
    if not os.path.exists(fpath):
        os.makedirs(fpath)
        print(f"Created folder: {fpath}")
    else:
        print(f"Folder already exists: {fpath}")

def disk_load(fname):
    if os.path.exists(fname):
        with open(fname, "rb") as f:
            return pickle.load(f)
    return None

def disk_save(fname, obj):
    with open(fname, "wb") as f:
        pickle.dump(obj, f)

# =========================
# SYMBOLS
# =========================
q  = Symbol('q')
t  = Symbol('t')
t1 = Symbol('t1')
t0 = Symbol('t0')
s  = Symbol('s')
s0 = Symbol('s0')
s1 = Symbol('s1')
a  = Symbol('alpha')
b  = Symbol('beta')
g  = Symbol('gamma')
z  = Symbol('z')
lambda1 = Symbol('lambda1')
lambda2 = Symbol('lambda2')

# =========================
# WORKER FUNCTIONS
# =========================
def iszero(x):
    return (simplify(x)==zeros(shape(x)[0],shape(x)[1]))

def sparseM(x):
    return  ImmutableSparseMatrix(x)

def KP(*factors):
    prod = Matrix([[1]])
    for x in factors:
        prod = kronecker_product(prod,x)
    return  sparseM(prod)

def qnum(k):
    return (q**k-q**(-k))*(q-q**(-1))**(-1)

def qfact(k):
    if k == 1 or k == 0:
        return 1
    return qnum(k)*qfact(k-1)

def qexp(x,k):
    if k==0:
        return x**0
    return qexp(x,k-1)+x**k*q**(comb(k,2))/qfact(k)

def subsst(x):
    return simplify(x.subs([(a,log(s)/log(q)),(b,log(t)/log(q))]))

def subst0(x):
    return factor(x.subs([(s,t0**Rational(1,2)),(q,(t0*t1)**Rational(-1,2))]))

def swapsst(x):
    return x.subs([(t,z), (s,t), (z,s)])

#%% Matrix Definitions
# each generator X is specified on the representation V(n,alpha) by "X(n)" or "X(n,a)"

# setup algorithm for determining matrix entries
def F1gen(n,i,j):
    if j%4 == 0 and i==j+4 and (j//4)+1<=n:
        return 1
    if j%4 == 1 and i==j+1:
        return 1
    if j%4 == 2 and i==j+4 and (j//4)+1<=n:
        return qnum(2)
    if j%4 == 2 and i==j+2*4-1 and (j//4)+2<=n:
        return -1
    if j%4 == 3 and i==j+4 and (j//4)+1<=n:
        return 1
    else:
        return 0

def F2gen(n,i,j):
    if j%4 == 0 and i==j+1:
        return 1
    if j%4 == 2 and i==j+1:
        return 1
    else:
        return 0

def E1gen(n,i,j):
    if j%4 == 0 and i==j-4:
        return qnum(j//4)*qnum(n-(j//4)+1)
    if j%4 == 1 and i==j-4:
        return qnum(j//4)*qnum(n-(j//4)+1)
    if j%4 == 2 and i==j-1:
        return qnum(n+1-2*(j//4))
    if j%4 == 2 and i==j-4 and (j//4)-1>=0:
        return qnum(j//4)*qnum(n-(j//4)+1)
    if j%4 == 3 and i==j-4 and (j//4)-1>=0:
        return qnum(j//4)*qnum(n-(j//4)+1)
    else:
        return 0

def E2gen(n,x,i,j):
    if j%4 == 1 and i==j-1:
        return qnum(x+(j//4))
    if j%4 == 2 and i==j+4-2 and (j//4)+1<=n:
        return qnum(x+(j//4))
    if j%4 == 3 and i==j-1:
        return qnum(x+(j//4)+1)
    if j%4 == 3 and i==j+4-2 and (j//4)+1<=n:
        return -qnum(x+(j//4))
    else:
        return 0

def K1gen(n,i,j):
    if i!=j:
        return 0
    if j%4 == 0:
        return q**(n-2*(j//4))
    if j%4 == 1:
        return q**(n-2*(j//4)+1)
    if j%4 == 2:
        return q**(n-2*(j//4+1)+1)
    if j%4 == 3:
        return q**(n-2*(j//4+1)+2)

def K2gen(n,x,i,j):
    if i!=j:
        return 0
    if j%4 == 0:
        return q**(x+(j//4))
    if j%4 == 1:
        return q**(x+(j//4))
    if j%4 == 2:
        return q**(x+((j//4)+1))
    if j%4 == 3:
        return q**(x+((j//4)+1))

def H1gen(n,i,j):
    if i!=j:
        return 0
    if j%4 == 0:
        return (n-2*(j//4))
    if j%4 == 1:
        return (n-2*(j//4)+1)
    if j%4 == 2:
        return (n-2*(j//4+1)+1)
    if j%4 == 3:
        return (n-2*(j//4+1)+2)

def H2gen(n,x,i,j):
    if i!=j:
        return 0
    if j%4 == 0:
        return (x+(j//4))
    if j%4 == 1:
        return (x+(j//4))
    if j%4 == 2:
        return (x+((j//4)+1))
    if j%4 == 3:
        return (x+((j//4)+1))

def bosgen(n,i,j):
    if i!=j:
        return 0
    if j%4 == 0:
        return (-1)**n
    if j%4 == 1:
        return (-1)**(n+1)
    if j%4 == 2:
        return (-1)**(n+1)
    if j%4 == 3:
        return (-1)**n

# function to output the matrix
def F1(n):
    return (ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: F1gen(n,i,j)))

def F2(n):
    return (ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: F2gen(n,i,j)))

def F12(n):
    return F2(n)*F1(n)-q**-1*F1(n)*F2(n)

def E1(n):
    return (ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: E1gen(n,i,j)))

def E2(n,x):
    return simplify(ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: E2gen(n,x,i,j)))

def E12(n,x):
    return E1(n)*E2(n,x)-q*E2(n,x)*E1(n)

def K1(n):
    return (ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: K1gen(n,i,j)))

def K2(n,x):
    return simplify(ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: K2gen(n,x,i,j)))

def H1(n):
    return (ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: H1gen(n,i,j)))

def H2(n,x):
    return (ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: H2gen(n,x,i,j)))

def bos(n):
    return (ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: bosgen(n,i,j)))

def Pgen(n,m,i,j):
    if i//(4*(n+1)) == j%(4*(m+1)) and i%(4*(n+1)) == j//(4*(m+1)):
        return 1
    else:
        return 0

def P(n,m):
    return (Matrix(16*(n+1)*(m+1), 16*(n+1)*(m+1), lambda i, j: Pgen(n,m,i,j)))

def pivot(n):
    return (bos(n)*subsst(K2(0,a)).inv()**2)

def coevgen(n,i):
    if i//(4*(n+1)) == i%(4*(n+1)):
        return 1
    else:
        return 0

def coev(n):
    return (Matrix(16*(n+1)**2, 1, lambda i,j: coevgen(n,i)))

def evgen(n,i):
    if i//(4*(n+1)) == i%(4*(n+1)):
        return 1
    else:
        return 0

def ev(n):
    return (Matrix(1,16*(n+1)**2, lambda i,j: evgen(n,j)))

#%% V(0,a) and dual

def matrix_constructions():
    ### V(0,a)
    id1 = eye(4)
    id11 = eye(16)
    ve1 = subsst(E1(0))
    ve2 = subsst(E2(0,a))
    vf1 = subsst(F1(0))
    vf2 = subsst(F2(0))
    vk1 =  subsst(K1(0))
    vk2 =  subsst(K2(0,a))
    vb = bos(0)
    vh1 = H1(0)
    vh2 = H2(0,a)

    ve12 = ve1*ve2-q*ve2*ve1
    vf12 = vf2*vf1-q**-1*vf1*vf2

    ### V(0,a)*
    se1 = (-ve1*vk1.inv()).transpose()
    se2 = (-vb*ve2*vk2.inv()).transpose()
    sf1 = (-vk1*vf1).transpose()
    sf2 = (-vb*vk2*vf2).transpose()
    sk1 = (vk1.inv()).transpose()
    sk2 = (vk2.inv()).transpose()
    sb = vb
    sh1 = (-vh1).transpose()
    sh2 = (-vh2).transpose()

    se12 = se1*se2-q*se2*se1
    sf12 = sf2*sf1-q**-1*sf1*sf2

    perm = P(0,0)

    ### Rvv
    exph1vv = simplify(-(2*KP(vh2, vh2)+KP(vh2,vh1)+KP(vh1,vh2)).diagonal())
    eHvv = subsst(s**2*diag(*[q**(i-exph1vv[0]) for i in exph1vv]))
    exp1vv   = id11 + (q-q**-1) * KP(ve1, vf1)
    exp12vv  = id11 - (q-q**-1) * KP(ve12*vb, vf12)
    exp2vv   = id11 - (q-q**-1) * KP(ve2*vb, vf2)
    Bosvv = (id11+KP(vb, id1)+KP(id1, vb)-KP(vb, vb))/2
    Rvv = perm*Bosvv*eHvv*exp2vv*exp12vv*exp1vv

    ### Rvs
    exph1vs = simplify(-(2*KP(vh2, sh2)+KP(vh2,sh1)+KP(vh1,sh2)).diagonal())
    eHvs = subsst(s**-2*diag(*[q**(i-exph1vs[0]) for i in exph1vs]))
    exp1vs   = id11 + (q-q**-1) * KP(ve1, sf1)
    exp12vs  = id11 - (q-q**-1) * KP(ve12*vb, sf12)
    exp2vs   = id11 - (q-q**-1) * KP(ve2*vb, sf2)
    Bosvs = (id11+KP(vb, id1)+KP(id1, sb)-KP(vb, sb))/2
    Rvs = perm*Bosvs*eHvs*exp2vs*exp12vs*exp1vs

    ### Rsv
    exph1sv = simplify(-(2*KP(sh2, vh2)+KP(sh2,vh1)+KP(sh1,vh2)).diagonal())
    eHsv = subsst(s**-2*diag(*[q**(i-exph1sv[0]) for i in exph1sv]))
    exp1sv   = id11 + (q-q**-1) * KP(se1, vf1)
    exp12sv  = id11 - (q-q**-1) * KP(se12*sb, vf12)
    exp2sv   = id11 - (q-q**-1) * KP(se2*sb, vf2)
    Bossv = (id11+KP(sb, id1)+KP(id1, vb)-KP(sb, vb))/2
    Rsv = perm*Bossv*eHsv*exp2sv*exp12sv*exp1sv

    ### Rss
    exph1ss = simplify(-(2*KP(sh2, sh2)+KP(sh2,sh1)+KP(sh1,sh2)).diagonal())
    eHss = subsst(s**2*diag(*[q**(i-exph1ss[0]) for i in exph1ss]))
    exp1ss   = id11 + (q-q**-1) * KP(se1, sf1)
    exp12ss  = id11 - (q-q**-1) * KP(se12*sb, sf12)
    exp2ss   = id11 - (q-q**-1) * KP(se2*sb, sf2)
    Bosss = (id11+KP(sb, id1)+KP(id1, sb)-KP(sb, sb))/2
    Rss = perm*Bosss*eHss*exp2ss*exp12ss*exp1ss

    h    = simplify(pivot(0))

    cl = coev(0)
    cr = KP(id1,h**-1)*coev(0)
    el = ev(0)
    er = ev(0)*KP(h, id1)

    return {'er':er, 'el':el, 'cr':cr, 'cl':cl, 'u':Rvv, 'ri':Rvs, 'li':Rsv, 'd':Rss,
          'ui':Rvv.inv(), 'r':Rsv.inv(), 'l':Rvs.inv(), 'di':Rss.inv(), 'id':id1}

def save_LG_RT():
    LG_RT_data = matrix_constructions()  # returns dict with keys Rvv, Rvs, Rsv, Rss, pivot, coev, ev
    RTfname = os.path.join(cwd, 'LG_RT_data.pkl')
    disk_save(RTfname, LG_RT_data)
    print("LG data saved to {LG_RT_data.pkl}")

#%%

# input a comma string drawing the tangle from bottom to top:
# this is the tangle word
# e.g: "0er, 1er, 0u, 0u, 1cl, 0cl"
# the first digit is #strands to left of the etangle
# e = evaluation, left/right
# c = coevaluation, left/right
# u, r, l, d = positive crossings
# ui, ri, li, di = negative crossings


def widener(elementary_tangle):
    if elementary_tangle in {"cl", "cr"}:
        return 2
    if elementary_tangle in {"el", "er"}:
        return - 2
    return 0

def filler(left_strands, bottom_width, elementary_tangle):
    dim = 4
    fname = os.path.join(pkl_dir, f"{elementary_tangle}_at_left_{left_strands}_bottom_{bottom_width}.pkl")
    slice_diagram = disk_load(fname)
    if slice_diagram is not None:
        return
    right_strands = bottom_width - left_strands - 2
    if elementary_tangle in {"cl", "cr"}:
        right_strands = bottom_width - left_strands
    if elementary_tangle in {'id'}:
        right_strands = bottom_width - left_strands - 1
    RTfname = os.path.join(cwd, 'LG_RT_data.pkl')
    RT = disk_load(RTfname)
    slice_diagram = KP(eye(dim**left_strands), RT[elementary_tangle], eye(dim**right_strands))
    disk_save(fname, slice_diagram)
    return

def filler_DMout(left_strands, bottom_width, elementary_tangle):
    dim = 4
    fname = os.path.join(pkl_dir, f"{elementary_tangle}_at_left_{left_strands}_bottom_{bottom_width}.pkl")
    slice_diagram = disk_load(fname)
    if slice_diagram is not None:
        return slice_diagram.to_DM()

def generate_pkls(entries,width):
    for entry in entries:
        left_strands = int(''.join([char for char in entry if char.isdigit()]))
        elementary_tangle = ''.join([char for char in entry if char.isalpha()])
        filler(left_strands, width, elementary_tangle)
        width += widener(elementary_tangle)

def tangle_invariant(tangle, bottom_width):
    dim = 4
    entries = re.split(',|;', tangle)
    width = bottom_width
    generate_pkls(entries,width)
    A = eye(dim**width).to_DM()
    for i in range(len(entries)//2):
        left_strands = int(''.join([char for char in entries[i] if char.isdigit()]))
        elementary_tangle = ''.join([char for char in entries[i] if char.isalpha()])
        M = filler_DMout(left_strands, width, elementary_tangle)
        A = M*A
        width += widener(elementary_tangle)
    B = eye(dim**width).to_DM()
    shift = len(entries)//2
    for j in range((len(entries)+1)//2):
        i = j + shift
        left_strands = int(''.join([char for char in entries[i] if char.isdigit()]))
        elementary_tangle = ''.join([char for char in entries[i] if char.isalpha()])
        M = filler_DMout(left_strands, width, elementary_tangle)
        B = M*B
        width += widener(elementary_tangle)
    return B*A

def compute_LG(name, tangle, bottom_width):
    dim = 4
    if name == '':
        return tangle_invariant(tangle, bottom_width).to_Matrix()
    start_time = time.perf_counter()
    output_matrix = tangle_invariant(tangle, bottom_width).to_Matrix()
    end_time = time.perf_counter()
    cwd = os.getcwd()
    fname = os.path.join(cwd, f"{name}.pkl")
    disk_save(fname, output_matrix)

    info_file_path = f"{name}_computation.txt"
    outputtime = str(end_time - start_time)
    with open(info_file_path, "w") as f:
        f.write(f"{name}: \nInput tangle: {tangle} \nBottom width: {bottom_width}\nElapsed time: {outputtime}")
        if shape(output_matrix) == (dim,dim):
            value = expand(output_matrix[0,0])
            if expand(output_matrix) != eye(dim)*value:
                print("Warning: A scalar matrix expected, but was not computed")
                f.write(f"Warning: A scalar matrix expected, but was not computed\nValue of LG?: {value}")
            else:
                factor_value = factor(value)
                f.write(f"\nValue of LG: {factor_value}")
    return output_matrix

def main():
    create_fold(pkl_dir)
    save_LG_RT()


if __name__ == "__main__":
    main()