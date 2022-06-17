import numpy as np
from sympy.utilities.iterables import multiset_permutations
import scipy.special
import scipy.sparse as spr
import scipy.sparse.linalg as splng
import matplotlib.pyplot as plt
from numba import jit

@jit
def Binom(n,k):
    ret = 1
    for i in range(min(k,n-k)):
        ret *= n-i
        ret /= i+1
    return ret

#state to label
@jit
def S2L(L,l_vec):
    N = L/2
    l_counter = 1
    n_counter = N
    for i in range(L):
        if l_vec[i] == 0:
            l_counter += Binom(L-(i+1),n_counter-1)
        if l_vec[i] == 1:
            n_counter -= 1
        if n_counter == 0:
            return l_counter

#label to state
@jit
def L2S(L,l):
    N = L/2
    n_counter = N
    l_counter = l
    l_vec = np.zeros(L)
    for i in range(L):
        binom = Binom(L-(i+1),n_counter-1)
        if n_counter == 0:
            l_vec[i] = 0
            continue
        if l_counter > binom:
            l_vec[i] = 0
            l_counter -= binom
            continue
        if l_counter <= binom:
            l_vec[i] = 1
            n_counter -= 1
    return l_vec

#state to label for subsystem
@jit
def S2L_Binary(vec,l):
    c = 0
    for i in range(l):
        if vec[l-1-i] == 1:
            c += 2**i
    return c

#perform hopping operations
def Hop(a):
    a = np.ndarray.astype(a,int)
    b_lst = []
    for i in range(L-1):
        if a[i] == 1 and a[i+1] == 0: #hopping to the right
            b = np.copy(a)
            b[i] = 0
            b[i+1] = 1
            b_lst.append(b)
        if a[i] == 0 and a[i+1] == 1: #hopping to the left
            b = np.copy(a)
            b[i] = 1
            b[i+1] = 0
            b_lst.append(b)
    if a[L-1] == 1 and a[0] == 0: #open boundaries right
        b = np.copy(a)
        b[L-1] = 0
        b[0] = 1
        b_lst.append(b)
    if a[L-1] == 0 and a[0] == 1: #open boundaries left
        b = np.copy(a)
        b[L-1] = 1
        b[0] = 0
        b_lst.append(b)
    ns = np.shape(b_lst)[0]
    return b_lst,ns

#diagonal elements
def Diag():
    Dim = int(Binom(L,L/2))
    row = np.empty(Dim*L)
    col = np.empty(Dim*L)
    data = np.empty(Dim*L)
    nnz = 0
    for i in range(Dim):
        A = L2S(L,i+1) #select vector
        v = 0 #nearest-neighbour interaction counter
        g = 0 #field interaction counter
        h = 0 #disorder interaction counter
        for k in range(L):
            g = g-gamma*k*A[k] #count field interactions
            h = h+np.random.normal(0,w)*A[k] #count disorder interactions
        for k in range(L-1):
            v = v+V*A[k]*A[k+1] #count nearest-neighbour interactions
        v = v+V*A[L-1]*A[0] #open boundary interaction
        row[nnz] = i
        col[nnz] = i
        data[nnz] = g+h+v
        nnz += 1
    return spr.csc_matrix((data[:nnz],(row[:nnz],col[:nnz])),shape=(Dim,Dim))

#off-diagonal elements
def Offdiag():
    Dim = int(Binom(L,L/2))
    row = np.empty(Dim*L)
    col = np.empty(Dim*L)
    data = np.empty(Dim*L)
    nnz = 0
    for i in range(Dim):
        A = L2S(L,i+1) #generate Fock state
        B_lst,ns = Hop(A) #perform hopping
        for k in range(ns): #loop through all states generated by hopping
            B = B_lst[k]
            j = int(S2L(L,B)-1) #convert to label
            row[nnz] = i
            col[nnz] = j
            data[nnz] = J
            nnz += 1
    return spr.csc_matrix((data[:nnz],(row[:nnz],col[:nnz])),shape=(Dim,Dim))

#Floquet operator
def Floquet():
    H_0 = Diag()
    F = splng.expm(-1j*H_0*T_0)@splng.expm(-1j*H_1*T_1)
    return F

def Weights():
    Dim = int(Binom(L,L/2))
    f = int(np.ceil(navg/Dim))
    mW = np.zeros((f,Dim))
    I = np.zeros(f)
    U = np.zeros((f,Dim,Dim),dtype=complex)
    for x in range(f):
        H = Floquet() #generate Floquet operator
        U[x] = np.linalg.eig(H.todense())[1] #calculate evecs
    return U

def IPR(q):
    f = int(np.ceil(navg/Dim))
    I = np.zeros(f)
    for x in range(f):
        I[x] = np.mean(np.sum(np.abs(U[x])**(2*q),axis=1))
    return np.mean(I)

#scaling exponent
def Tau():
    tau = np.zeros(nq)
    for i in range(nq):
        q = q_lst[i]
        I = IPR(q) #IPR at q
        tau[i] = -np.log(I)/np.log(Dim) #scaling exponent
    return tau

#extrapolated scaling exponent
def Tau_ext():
    tau_ext = np.zeros((nl-1,nq))
    m = np.zeros((nl-1,nq))
    Dim_2 = int(scipy.special.comb(L_lst[nl-1],L_lst[nl-1]/2)) #larger dimensionality
    for i in range(nl-1):
        a = 1/np.log(Dim_2)-1/np.log(Dim_lst[i])
        for j in range(nq):
            m[i,j] = (tau[nl-1,j]-tau[i,j])/a #gradient of tau vs 1/lnN
            tau_ext[i,j] = tau[nl-1,j]-m[i,j]/np.log(Dim_2) #y-intercept
    return m,tau_ext

#Parameters
L_lst = np.array([6,8,10,12,14,16]) #system size
J = 0.5 #hopping amplitude
V = 0.5 #interaction strength
gamma = 5 #potential strength
w = 0.05 #disorder strength
T_1 = 0.1 #kicked period
f = 2 #fraction
T_0 = 2*np.pi/gamma*f #unkicked period
q_lst = np.linspace(0,3,16) #fractal dimension parameter
navg = 3000 #total averaging

#calculate extrapolated scaling exponent
nq = np.size(q_lst)
nl = np.size(L_lst)
tau = np.zeros((nl,nq))
Dim_lst = np.zeros(nl)
for i in range(nl):
    L = L_lst[i]
    Dim_lst[i] = int(Binom(L,L/2))
    Dim = Dim_lst[i]
    #H_1 = Offdiag()
    #U = Weights()
    #tau[i] = Tau()

tau = np.load('Stark_Files/Data/Tau/Clean/Tau_Clean_2.0.npy')

m,tau_ext = Tau_ext()

#extrapolation
l_1 = L_lst[nl-2]
l_2 = 200
x_lst = np.zeros(2)
x_lst[0] = int(Binom(l_1,l_1/2))
x_lst[1] = int(Binom(l_2,l_2/2))

plt.figure(figsize=(8,6),dpi=80)
plt.subplot(1,2,1)
c = 0
for i in range(1,nq,2):
    plt.plot(1/np.log(Dim_lst),tau[:,i],'C{}-x'.format(c),label=r'$q={}$'.format(np.round(q_lst[i],1)))
    plt.plot(1/np.log(x_lst),m[nl-2,i]/np.log(x_lst)+tau_ext[nl-2,i],'C{}--'.format(c))
    c = c+1
plt.xlabel(r'$1/\ln\mathcal{D}$')
plt.ylabel(r'$\tau_q$')
plt.ylim([-1,2])
#plt.legend(ncol=2,fontsize=10)
ax = plt.gca()
plt.annotate('(i)',(0.92,0.03),xycoords='axes fraction')
plt.annotate(r'$\gamma T_0/2\pi=2.0$',(0.02,0.92),xycoords='axes fraction')
ax.set_aspect(0.109)
plt.subplot(1,2,2)
for i in range(nl-1):
    plt.plot(q_lst,tau_ext[i],'-x',label=r'$L\in[{},{}]$'.format(L_lst[i],L_lst[nl-1]))
plt.axhline(y=0,linestyle='--',xmin=0.35,color='k'.format(i))
plt.annotate('CUE',(0.89,0.86),xycoords='axes fraction')
plt.annotate('Loc.',(0.9,0.28),xycoords='axes fraction')
plt.annotate('(j)',(0.92,0.03),xycoords='axes fraction')
plt.plot(q_lst,q_lst-1,'k--')
plt.xlabel(r'$q$')
plt.ylim([-1,2])
#plt.legend()
plt.tight_layout()
ax = plt.gca()
ax.set_yticklabels([])
ax.set_aspect(1)