import numpy
import math


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



# frobenius norm of generic matrix
def frobnorm2(M):
    return frobinnprod(M, M).real

# frobenius norm of hermitian matrix
def frobnorm2_h(H):
    return numpy.trace(numpy.matmul(H, H)).real

# frobenius norm of antihermitian matrix
def frobnorm2_a(A):
    return -numpy.trace(numpy.matmul(A, A)).real

# computes angle between [H_i, H_j] and L_k
def angle2_case1(Hi, Hj, Lk):
    
    COMM = numpy.matmul(Hi, Hj) - numpy.matmul(Hj, Hi)
    num = math.pow(simpinnprod(COMM, Lk, 1), 2.)
    den = frobnorm2_a(COMM)*frobnorm2_h(Lk)

    return (num/den)

# computes angle between [L_i, L_j] and L_k
def angle2_case2(Li, Lj, Lk):
    
    COMM = numpy.matmul(Li, Lj) - numpy.matmul(Lj, Li)
    num = math.pow(simpinnprod(COMM, Lk, -1), 2.)
    den = frobnorm2_a(COMM)*frobnorm2_h(Lk)

    return (num/den)

# computes angle between [L_i, H_j] and H_k
def angle2_case3(Li, Hj, Hk):
    
    COMM = numpy.matmul(Li, Hj) - numpy.matmul(Hj, Li)
    num = math.pow(simpinnprod(COMM, Hk, -1), 2.)
    den = frobnorm2_a(COMM)*frobnorm2_h(Hk)

    return (num/den)





# input multicode
multicode = input("Code\n")

# read single G codes and store them in args
filenameARGS = multicode + "_varG_args.txt"
with open(filenameARGS) as r_ARGS:
    args = r_ARGS.readline().split(" ")
args = [item for item in args if len(item) >= 5]

filenameS_cumul = multicode + "_varG_S.txt"
filenameTR_cumul = multicode + "_varG_TR.txt"
filenameF_cumul = multicode + "_varG_F.txt"
filenameCOMM1_cumul = multicode + "_varG_COMM1.txt"
filenameCOMM2_cumul = multicode + "_varG_COMM2.txt"
filenameCOMM3_cumul = multicode + "_varG_COMM3.txt"


with open(filenameS_cumul, 'w') as w_S_cumul, open(filenameTR_cumul, 'w') as w_TR_cumul, open(filenameF_cumul, 'w') as w_F_cumul, open(filenameCOMM1_cumul, 'w') as w_COMM1_cumul, open(filenameCOMM2_cumul, 'w') as w_COMM2_cumul, open(filenameCOMM3_cumul, 'w') as w_COMM3_cumul:

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
        

        filenameS = arg + "_simS.txt"
        filenameHL = arg + "_simHL.txt"
        filenameTR = arg + "_simTR.txt"
        filenameCOMM1 = arg + "_simCOMM1.txt"
        filenameCOMM2 = arg + "_simCOMM2.txt"
        filenameCOMM3 = arg + "_simCOMM3.txt"

        S = []
        TR = []
        F = []
        COMM1 = []
        COMM2 = []
        COMM3 = []
        with open(filenameHL) as r_HL, open(filenameS) as r_S, open(filenameTR, 'w') as w_TR, open(filenameCOMM1, 'w') as w_COMM1, open(filenameCOMM2, 'w') as w_COMM2, open(filenameCOMM3, 'w') as w_COMM3:

            for h in range(int(math.ceil(Nsw/GAP))):
                
                # write action
                S.append(float(r_S.readline()))
                
                # write traces
                mats = read_mat(nH, nL, dim, r_HL)
                mats2 = [numpy.matmul(mat, mat) for mat in mats]
                mtrless = [ mat - (numpy.trace(mat)/dim)*numpy.identity(dim) for mat in mats]
                mtrless2 = [numpy.matmul(mat, mat) for mat in mtrless]

                TR_temp = []
                for i in range(nH+nL):
                    temp1 = numpy.trace(mats[i]).real
                    temp11 = temp1*temp1
                    temp2 = numpy.trace(mats2[i]).real
                    temp3 = numpy.trace(mtrless2[i]).real
                    TR_temp.append(temp1)
                    TR_temp.append(temp11)
                    TR_temp.append(temp2)
                    TR_temp.append(temp3)
                    w_TR.write(str(temp1) + ' ' + str(temp11) + ' ' + str(temp2) + ' ' + str(temp3) + ' ')
                w_TR.write('\n')
                TR.append(TR_temp)


                # write F
                num = 0
                den = 0
                for i in range(nH):
                    temp1 = numpy.trace(mats[i]).real
                    temp2 = numpy.trace(mats2[i]).real
                    num += temp1*temp1
                    den += temp2

                try:
                    F.append( num/(dim*den) )
                except:
                    pass


                # write commutators
                
                # CASE 1
                COMM1_temp = []
                for i in range(nH):
                    for j in range(i+1, nH):
                        for k in range(nH, nH+nL):
                            temp = angle2_case1(mats[i], mats[j], mats[k])
                            COMM1_temp.append(temp)
                            w_COMM1.write(str(temp) + ' ')
                w_COMM1.write('\n')
                COMM1.append(COMM1_temp)

                # CASE 2
                COMM2_temp = []
                for i in range(nH, nH+nL):
                    for j in range(i+1, nH+nL):
                        for k in range(nH, nH+nL):
                            temp = angle2_case2(mats[i], mats[j], mats[k])
                            COMM2_temp.append(temp)
                            w_COMM2.write(str(temp) + ' ')
                w_COMM2.write('\n')
                COMM2.append(COMM2_temp)

                # CASE 3
                COMM3_temp = []
                for i in range(nH, nH+nL):
                    for j in range(nH):
                        for k in range(nH):
                            temp = angle2_case3(mats[i], mats[j], mats[k])
                            COMM3_temp.append(temp)
                            w_COMM3.write(str(temp) + ' ')
                w_COMM3.write('\n')
                COMM3.append(COMM3_temp)


        
        w_S_cumul.write(str(G) + ' ' + str(numpy.mean(S)) + ' ' + str(bootstrap(S, 100)) + '\n')
        
        if len(F):
            var, var_err = bootstrap_var(F, 100)
            w_F_cumul.write(str(G) + ' ' + str(numpy.mean(F)) + ' ' + str(bootstrap(F, 100)) + ' ' + str(var) + ' ' + str(var_err) + '\n')


        if len(COMM1):
            COMM1 = [list(x) for x in zip(*COMM1)]
            w_COMM1_cumul.write(str(G) + ' ')
            for item in COMM1:
                w_COMM1_cumul.write(str(numpy.mean(item)) + ' ' + str(bootstrap(item, 100)) + ' ')
            w_COMM1_cumul.write('\n')
        if len(COMM2):
            COMM2 = [list(x) for x in zip(*COMM2)]
            w_COMM2_cumul.write(str(G) + ' ')
            for item in COMM2:
                w_COMM2_cumul.write(str(numpy.mean(item)) + ' ' + str(bootstrap(item, 100)) + ' ')
            w_COMM2_cumul.write('\n')
        if len(COMM3):
            COMM3 = [list(x) for x in zip(*COMM3)]
            w_COMM3_cumul.write(str(G) + ' ')
            for item in COMM3:
                w_COMM3_cumul.write(str(numpy.mean(item)) + ' ' + str(bootstrap(item, 100)) + ' ')
            w_COMM3_cumul.write('\n')

        TR = [list(x) for x in zip(*TR)]
        w_TR_cumul.write(str(G) + ' ')
        for item in TR:
            w_TR_cumul.write(str(numpy.mean(item)) + ' ' + str(bootstrap(item, 100)) + ' ')
        w_TR_cumul.write('\n')







