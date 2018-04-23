import numpy


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


# frobenius inner product btween generic matrices
def frobinnprod(M1, M2):
    return numpy.trace(numpy.matmul(M1.getH(), M2))

# frobenius norm of generic matrix
def frobnorm(M):
    return numpy.sqrt(frobinnprod(M, M).real)

# frobenius norm of hermitian matrix
def frobnorm_h(H):
    return numpy.sqrt(numpy.trace(numpy.matmul(H, H)).real)

# frobenius norm of antihermitian matrix
def frobnorm_a(A):
    return numpy.sqrt(-numpy.trace(numpy.matmul(A, A)).real)

# simplified version of inner product for our purposes
def simpinnprod(A, B, sign):
    return sign*(numpy.trace(numpy.matmul(A, B)).imag)

# computes angle between [H_i, H_j] and L_k
def angle_case1(Hi, Hj, Lk):
    
    COMM = numpy.matmul(Hi, Hj) - numpy.matmul(Hj, Hi)
    num = simpinnprod(COMM, Lk, 1)
    den = frobnorm_a(COMM)*frobnorm_h(Lk)

    return numpy.arccos(num/den)

# computes angle between [L_i, L_j] and L_k
def angle_case2(Li, Lj, Lk):
    
    COMM = numpy.matmul(Li, Lj) - numpy.matmul(Lj, Li)
    num = simpinnprod(COMM, Lk, -1)
    den = frobnorm_a(COMM)*frobnorm_h(Lk)

    return numpy.arccos(num/den)

# computes angle between [L_i, H_j] and H_k
def angle_case3(Li, Hj, Hk):
    
    COMM = numpy.matmul(Li, Hj) - numpy.matmul(Hj, Li)
    num = simpinnprod(COMM, Hk, -1)
    den = frobnorm_a(COMM)*frobnorm_h(Hk)

    return numpy.arccos(num/den)



# input multicode
multicode = input("Code\n")

# read single G codes and store them in args
filenameARGS = multicode + "_varG_args.txt"
with open(filenameARGS) as r_ARGS:
    args = r_ARGS.readline().split(" ")
args = [item for item in args if len(item) >= 5]

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
    filenameTR = arg + "_simTR.txt"
    filenameCOMM1 = arg + "_simCOMM1.txt"
    filenameCOMM2 = arg + "_simCOMM2.txt"
    filenameCOMM3 = arg + "_simCOMM3.txt"

    with open(filenameHL) as r_HL, open(filenameTR, 'w') as w_TR, open(filenameCOMM1, 'w') as w_COMM1, open(filenameCOMM2, 'w') as w_COMM2, open(filenameCOMM3, 'w') as w_COMM3:

        for h in range(int(Nsw/GAP)):
            mats = read_mat(nH, nL, dim, r_HL)
            mats2 = [numpy.matmul(mat, mat) for mat in mats]
            mtrless = [ mat - (numpy.trace(mat)/dim)*numpy.identity(dim) for mat in mats]
            mtrless2 = [numpy.matmul(mat, mat) for mat in mtrless]

            for i in range(nH+nL):
                w_TR.write(str(numpy.trace(mats[i]).real) + ' ' + str(numpy.trace(mats2[i]).real) + ' ' + str(numpy.trace(mtrless2[i]).real) + ' ')
            w_TR.write('\n')

            # write commutators
            
            # CASE 1
            for i in range(nH):
                for j in range(i+1, nH):
                    for k in range(nH, nH+nL):
                        w_COMM1.write(str(angle_case1(mats[i], mats[j], mats[k])) + ' ')
            w_COMM1.write('\n')

            # CASE 2
            for i in range(nH, nH+nL):
                for j in range(i+1, nH+nL):
                    for k in range(nH, nH+nL):
                        w_COMM2.write(str(angle_case2(mats[i], mats[j], mats[k])) + ' ')
            w_COMM2.write('\n')

            # CASE 3
            for i in range(nH, nH+nL):
                for j in range(nH):
                    for k in range(nH):
                        w_COMM3.write(str(angle_case3(mats[i], mats[j], mats[k])) + ' ')
            w_COMM3.write('\n')





