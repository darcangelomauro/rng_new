import numpy
import scipy
import math
from numpy import matrix
from numpy import linalg
from scipy import linalg
import copy
import itertools


def decomp(p, q):
    if p != 0:
        if p % 2 == 0:
            num20 = p/2
            num10 = 0
        else:
            num20 = (p-1)/2
            num10 = 1
    else:
        num20 = 0
        num10 = 0

    if q != 0:
        if q % 2 == 0:
            num02 = q/2
            num01 = 0
        else:
            num02 = (q-1)/2
            num01 = 1
    else:
        num02 = 0
        num01 = 0

    num11 = 0
    if num10 == 1 and num01 == 1:
        num10 = 0
        num01 = 0
        num11 = 1
    

    return [int(num20), int(num02), int(num11), int(num10), int(num01)]

def build_chiral(gamma, s):
    chiral = gamma[0]
    for i in range(len(gamma)-1):
        chiral = chiral*gamma[i+1]
    factor = pow(1j, s*(s+1.)/2.)
    chiral = chiral*factor

    return chiral

def tensor_modules(gamma1, gamma2, s):
    
    size1 = math.sqrt(gamma1[0].size)
    size2 = math.sqrt(gamma2[0].size)

    chiral1 = build_chiral(gamma1, s)
    id2 = numpy.identity(int(size2))

    res = []

    for item1 in gamma1:
        res.append(numpy.kron(item1, id2))
    for item2 in gamma2:
        res.append(numpy.kron(chiral1, item2))

    return res

def order(gamma):
    herm = []
    antiherm = []
    for i in range(len(gamma)):
        if (gamma[i].H == gamma[i]).all():
            herm.append(gamma[i])
        else:
            antiherm.append(gamma[i])
    gamma = [item for item in herm]
    for item in antiherm:
        gamma.append(item)
    return gamma

def order_tuple(res):
    herm = []
    antiherm = []
    for i in range(len(res)):
        if (res[i][1].H == res[i][1]).all():
            herm.append((res[i][0] + "_H", res[i][1]))
        else:
            antiherm.append((res[i][0] + "_A", 1j*res[i][1]))

    res = [item for item in herm]
    for item in antiherm:
        res.append(item)
    return res

def masterchef(p, q):

    s = 0

    structure = decomp(p, q)
    print("num20: " + str(structure[0]))
    print("num02: " + str(structure[1]))
    print("num11: " + str(structure[2]))
    print("num10: " + str(structure[3]))
    print("num01: " + str(structure[4]))

    gamma_10 = [matrix([[1]])]
    gamma_01 = [matrix([[-1j]])]
    gamma_20 = [matrix([[1, 0], [0, -1]]), matrix([[0, 1], [1, 0]])]
    gamma_02 = [matrix([[1j, 0], [0, -1j]]), matrix([[0, 1], [-1, 0]])]
    gamma_11 = [matrix([[1, 0], [0, -1]]), matrix([[0, 1], [-1, 0]])]

    gamma_pq = []

    if structure[0] != 0:
        gamma_pq = [item for item in gamma_20]
        s = -2 % 8
        for i in range(structure[0]-1):
            temp = tensor_modules(gamma_pq, gamma_20, s)
            s = (s-2) % 8
            gamma_pq = [item for item in temp]
        for i in range(structure[1]):
            temp = tensor_modules(gamma_pq, gamma_02, s)
            s = (s+2) % 8
            gamma_pq = [item for item in temp]
        if structure[2]:
            temp = tensor_modules(gamma_pq, gamma_11, s)
            gamma_pq = [item for item in temp]
        else:
            if structure[3]:
                temp = tensor_modules(gamma_pq, gamma_10, s)
                s = (s-1) % 8
                gamma_pq = [item for item in temp]
            elif structure[4]:
                temp = tensor_modules(gamma_pq, gamma_01, s)
                s = (s+1) % 8
                gamma_pq = [item for item in temp]
    elif structure[1] != 0:
        gamma_pq = [item for item in gamma_02]
        s = 2 % 8
        for i in range(structure[1]-1):
            temp = tensor_modules(gamma_pq, gamma_02, s)
            s = (s+2) % 8
            gamma_pq = [item for item in temp]
        if structure[2]:
            temp = tensor_modules(gamma_pq, gamma_11, s)
            gamma_pq = [item for item in temp]
        else:
            if structure[3]:
                temp = tensor_modules(gamma_pq, gamma_10, s)
                s = (s-1) % 8
                gamma_pq = [item for item in temp]
            elif structure[4]:
                temp = tensor_modules(gamma_pq, gamma_01, s)
                s = (s+1) % 8
                gamma_pq = [item for item in temp]
    elif structure[2] != 0:
        gamma_pq = [item for item in gamma_11]
        s = 0 % 8
    else:
        if structure[3]:
            gamma_pq = [item for item in gamma_10]
            s = -1 % 8
        elif structure[4]:
            gamma_pq = [item for item in gamma_01]
            s = 1 % 8
        else:
            print("module (p, q) not recognized")

    print("s: " + str(s))

    return order(gamma_pq)

def comm(A):
    dim = A.shape[0]
    I = numpy.identity(dim)
    return scipy.linalg.kron(A, I) - scipy.linalg.kron(I, A.T)

def anticomm(A):
    dim = A.shape[0]
    I = numpy.identity(dim)
    return scipy.linalg.kron(A, I) + scipy.linalg.kron(I, A.T)

def innprod(a, b):
    c = a.H*b
    return numpy.trace(c)


def indprod(gamma, p, q):
    newgamma = [item for item in gamma]
    newgamma_idx = [str(i) for i in range(len(gamma))]
    proposed_even = []
    proposed_even_idx = []
    proposed_odd = []
    proposed_odd_idx = []
    previous = []
    previous_idx = []
    s = (q-p) % 8
    n = p+q

    #products of two
    for i in range(len(gamma)):
        for j in range(len(gamma)):
            k = str(i)+str(j)
            v = gamma[i]*gamma[j]
            proposed_even.append(v)
            proposed_even_idx.append(k)
            previous.append(v)
            previous_idx.append(k)

    #other products
    for i in range(n-2):
        temp = []
        temp_idx = []
        for j in range(len(previous)):
            for k in range(len(gamma)):
                newkey = str(previous_idx[j]) + str(k)
                newval = previous[j]*gamma[k]
                if i % 2 == 0:
                    proposed_odd.append(newval)
                    proposed_odd_idx.append(newkey)
                else:
                    proposed_even.append(newval)
                    proposed_even_idx.append(newkey)
                temp.append(newval)
                temp_idx.append(newkey)
        previous = [item for item in temp]
        previous_idx = [item for item in temp_idx]

    for i in range(len(proposed_odd)):
        check = 1
        for j in range(len(newgamma)):
            if innprod(proposed_odd[i], newgamma[j]) != 0:
                check = 0
        if check:
            newgamma.append(proposed_odd[i])
            newgamma_idx.append(proposed_odd_idx[i])

    if s % 2 != 0:
        for i in range(len(proposed_even)):
            check = 1
            for j in range(len(newgamma)):
                if innprod(proposed_even[i], newgamma[j]) != 0:
                    check = 0
            if check:
                newgamma.append(proposed_even[i])
                newgamma_idx.append(proposed_even_idx[i])
            
    res = []
    for i in range(len(newgamma)):
        res.append((newgamma_idx[i], newgamma[i]))

    return order_tuple(res)


def read_mat(nH, nL, dim, f):
    MAT_ls = []
    for i in range(nH+nL):
        temp_str = f.readline()
        MAT_ls.append(temp_str.split(' '))

    MAT_fl = []
    for item in MAT_ls:
        MAT_fl.append([float(item1) for item1 in item if item1 != '\n'])

    MAT = []
    for i in range(nH+nL):
        MAT.append(numpy.zeros(shape=(dim,dim), dtype=complex))


    for i in range(dim):
        for j in range(dim):
            
            idx = 2*(i*dim + j)

            for k in range(nH+nL):
                MAT[k][i, j] = complex(MAT_fl[k][idx], MAT_fl[k][idx+1])

    return MAT


# bootstrap resampling for error estimates
def bootstrap(ls, nresamp):

    N = len(ls)

    vals = []

    # repeat resampling 'nresamp' times
    for i in range(nresamp):
        
        # list of resampled data
        new_ls = []
        
        # build new list of N elements with repetitions 
        for j in range(N):
            # pick random index
            idx = numpy.random.randint(N)
            # append that element in new list
            new_ls.append(ls[idx])

        vals.append(numpy.mean(new_ls))

    return math.sqrt(numpy.var(vals))

# bootstrap resampling for variance
def bootstrap_var(ls, nresamp):

    N = len(ls)

    vals = []

    # repeat resampling 'nresamp' times
    for i in range(nresamp):
        
        # list of resampled data
        new_ls = []
        new_ls2 = []
        
        # build new list of N elements with repetitions 
        for j in range(N):
            # pick random index
            idx = numpy.random.randint(N)
            # append that element in new list
            new_ls.append(ls[idx])
            new_ls2.append(ls[idx]*ls[idx])

        mean = numpy.mean(new_ls)
        mean2 = numpy.mean(new_ls2)
        vals.append(mean2 - mean*mean)

    return numpy.mean(vals), math.sqrt(numpy.var(vals))


def build_dirac(gammas, mats, nH, dim, dimG):
    size = dim*dim*dimG
    D = numpy.zeros((size, size), dtype=complex)

    for i, gamma in enumerate(gammas):
        if i<nH:
            A = scipy.linalg.kron(gamma[1], anticomm(mats[i]))
            D = numpy.add(D, A)
        else:
            A = scipy.linalg.kron(gamma[1], comm(mats[i]))
            D = numpy.add(D, A)

    return D





# input multicode
multicode = input("Code\n")

# read single G codes and store them in args
filenameARGS = multicode + "_varG_args.txt"
with open(filenameARGS) as r_ARGS:
    args = r_ARGS.readline().split(" ")
args = [item for item in args if len(item) >= 5]

filenameEIG_cumul = multicode + "_varG_EIG.txt"

p = int(input('p: '))
q = int(input('q: '))
gamma = masterchef(p, q)
print("************products***********")
gamma_prod = indprod(gamma, p, q)
for item in gamma_prod:
    print("gamma_" + str(item[0]))
    print(item[1])
dimG = gamma_prod[0][1].shape[0]
print('dimG = ' + str(dimG))


with open(filenameEIG_cumul, 'w') as w_EIG_cumul:

    # iterate on args
    for arg in args:
        print("Processing: " + arg)
        
        # read parameters
        filenameDATA = arg + "_data.txt"
        r_DATA = open(filenameDATA, "r")
        params = r_DATA.readline().split(" ")
        r_DATA.close()
        dim = int(params[0])
        nH = int(params[1])
        nL = int(params[2])
        SCALE = float(params[3])
        G = float(params[4])
        Nsw = int(params[5])
        GAP = int(params[6])
        

        filenameHL = arg + "_simHL.txt"
        filenameEIG = arg + "_simEIG.txt"

        EIG = []
        with open(filenameHL) as r_HL, open(filenameEIG, 'w') as w_EIG:

            for h in range(int(math.ceil(Nsw/GAP))):
                
                # write traces
                mats = read_mat(nH, nL, dim, r_HL)

                D = build_dirac(gamma_prod, mats, nH, dim, dimG)

                eigvals = scipy.linalg.eigvalsh(D)
                for val in eigvals:
                    w_EIG.write(str(val) + ' ')
                w_EIG.write('\n')
                EIG.append(eigvals)

        EIG = [list(x) for x in zip(*EIG)]
        w_EIG_cumul.write(str(G) + ' ')
        for item in EIG:
            w_EIG_cumul.write(str(numpy.mean(item)) + ' ' + str(bootstrap(item, 100)) + ' ')
        w_EIG_cumul.write('\n')





