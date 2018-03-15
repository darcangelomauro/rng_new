import numpy

multicode = input("Code\n")
filenameARGS = multicode + "_varG_args.txt"
filenameCURVE = multicode + "_varG_F.txt"
r_ARGS = open(filenameARGS, "r")
w_CURVE = open(filenameCURVE, "w")
args = r_ARGS.readline().split(" ")
args = [item for item in args if len(item) >= 5]
r_ARGS.close()

for code in args:
    print("Processing: " + code)
    filenameHL = code + "_simHL.txt"
    filenameF = code + "_simF.txt"
    filenameTR = code + "_simTR.txt"

    r_HL = open(filenameHL, "r")
    w_F = open(filenameF, "w")
    w_TR = open(filenameTR, "w")

    # read parameters
    filenameD = code + "_data.txt"
    r_D = open(filenameD, "r")
    params = r_D.readline().split(" ")
    r_D.close()
    dim = int(params[0])
    nH = int(params[1])
    nL = int(params[2])
    SCALE = float(params[3])
    G = float(params[4])
    Nsw = int(params[5])
    GAP = int(params[6])

    obs = []
    for h in range(int(Nsw/GAP)):

        if not h%1000:
            print(h)
        
        H_ls = []
        for i in range(nH+nL):
            temp_str = r_HL.readline()
            temp_ls = temp_str.split(" ")
            if i < nH:
                H_ls.append(temp_ls)

        H_fl = []
        for item in H_ls:
            H_fl.append([float(item1) for item1 in item if item1 != '\n'])


        H = []
        H2 = []
        for i in range(nH):
            H.append(numpy.zeros(shape=(dim,dim), dtype=complex))


        for i in range(dim):
            for j in range(dim):
                
                idx = 2*(i*dim + j)

                for k in range(nH):
                    H[k][i, j] = complex(H_fl[k][idx], H_fl[k][idx+1])

        for mat in H:
            H2.append(numpy.matmul(mat, mat))


        numerator = 0
        denominator = 0
        for i in range(nH):
            temp = numpy.real(numpy.trace(H[i]))
            temp2 = numpy.real(numpy.trace(H2[i]))
            w_TR.write(str(temp) + " ")
            numerator += temp*temp
            denominator += temp2

        w_TR.write("\n")
        w_F.write(str(numerator/(dim*denominator)) + "\n")
        obs.append(numerator/(dim*denominator))

    bin_obs = numpy.array_split(obs, 100)
    bin_obs = [numpy.mean(item) for item in bin_obs]

    mean = numpy.mean(bin_obs)
    var = numpy.var(bin_obs)

    w_CURVE.write(str(G) + " " + str(mean) + " " + str(numpy.sqrt(var/len(bin_obs))) + " " + str(var) + "\n")


        
    r_HL.close()
    w_F.close()
    w_TR.close()

w_CURVE.close()
    





