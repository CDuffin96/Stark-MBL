import numpy as np
from sympy.utilities.iterables import multiset_permutations
import scipy.special
import scipy.sparse as spr
import scipy.sparse.linalg as splng
import matplotlib.pyplot as plt

#Swapping operators
def Operators(L):
    c = []
    for i in range(L-1):
        op = np.eye(L)
        op[i,i] = 0
        op[i+1,i+1] = 0
        op[i,i+1] = 1
        op[i+1,i] = 1
        c.append(op)
    op = np.eye(L) #open boundaries
    op[L-1,L-1] = 0
    op[0,0] = 0
    op[L-1,0] = 1
    op[0,L-1] = 1
    c.append(op)
    return c

#Fock states
def Hilbert(L):
    a = np.zeros(L)
    a[0:int(L/2)] = 1 #fill left half of sites
    a = np.ndarray.astype(a,int)
    vec = list(multiset_permutations(a)) #find all permutations
    bits = np.packbits(vec,axis=-1) #decimal values of Fock states
    return vec,bits

#locate regions where swap operations are allowed
def Swap(a,n=2):
    ret = np.cumsum(a,axis=0,dtype=float)
    ret[n:] = ret[n:]-ret[:-n]
    swap = np.where(ret[n-1:]==1)
    return swap[0],np.size(swap[0])

def Diag(L):
    Dim = int(scipy.special.comb(L,L/2))
    H = np.zeros((Dim,Dim)) #Hamiltonian
    c = Operators(L)
    vec,bits = Hilbert(L)
    for i in range(Dim):
        A = vec[i] #select vector
        v = 0 #interaction counter
        for k in range(L):
            n = -gamma*k #tilted potential
            H[i,i] = H[i,i]+n*A[k] #apply tilted potential
            H[i,i] = H[i,i]+np.random.normal(0,w)*A[k] #apply disorder
        for k in range(L-1):
            v = v+A[k]*A[k+1] #count interactions
        v = v+A[L-1]*A[0] #open boundary interaction
        H[i,i] = H[i,i]+V*v #apply interactions
    return H

def Offdiag(L):
    Dim = int(scipy.special.comb(L,L/2))
    H = np.zeros((Dim,Dim)) #Hamiltonian
    c = Operators(L)
    vec,bits = Hilbert(L)
    for i in range(Dim):
        A = vec[i]
        swap,ns = Swap(A) #find sites where swap is allowed
        for k in range(ns):
            k = swap[k]
            B = c[k]@A #apply swap operator
            j = np.where((np.packbits(B.astype(int))==bits).all(axis=1))[0][0] #find index of new vector
            H[i,j] = 1 #apply hopping elements to matrix
            H[j,i] = 1
        if A[0]+A[L-1] == 1: #open boundaries
            B = c[L-1]@A
            j = np.where((np.packbits(B.astype(int))==bits).all(axis=1))[0][0]
            H[i,j] = 1
            H[j,i] = 1
    return H

def IPR(L,T_0):
    Dim = int(scipy.special.comb(L,L/2))
    f = int(np.ceil(navg/Dim))
    I = np.zeros(f)
    for x in range(f):
        H_Diag = Diag(L)
        F = splng.expm(-1j*H_Diag*T_0)@splng.expm(-1j*J*H_Offdiag*T_1) #Floquet operator
        U = np.linalg.eig(F)[1] #calculate evecs
        I[x] = np.mean(np.sum(np.abs(U)**4,axis=1)) #average over evecs
    return np.mean(I)

#Parameters
L_lst = np.array([8,10,12,14,16]) #system size
J = 0.5 #hopping amplitude
V = 0.5 #interaction strength
gamma = 5 #potential strength
w = 0.05 #disorder strength
navg = 3000 #total averaging
T_1 = 0.1 #kicked period
T_0_lst = 2*np.pi*np.arange(2.5,25,2.5)*0.1/gamma #unkicked period

ng = np.size(T_0_lst)
nl = np.size(L_lst)
I = np.zeros((nl,ng))
Dim_lst = np.zeros(nl)
for i in range(nl):
    L = L_lst[i]
    Dim_lst[i] = int(scipy.special.comb(L,L/2))
    H_Offdiag = Offdiag(L)
    for j in range(ng):
        I[i,j] = IPR(L,T_0_lst[j])

plt.figure(1)
for i in range(nl):
    plt.plot(T_0_lst*gamma/(2*np.pi),I[i],'-x',label=r'$L={}$'.format(L_lst[i]))
plt.xlabel(r'$\gamma T_0/2\pi$')
plt.ylabel(r'$I(T_0)$')
plt.legend(ncol=2,fontsize=10)
