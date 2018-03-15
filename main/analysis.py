import numpy

code = input("Code\n")
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

for h in range(Nsw):

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
    for i in range(nH):
        H.append(numpy.zeros(shape=(dim,dim), dtype=complex))


    for i in range(dim):
        for j in range(dim):
            
            idx = 2*(i*dim + j)

            for k in range(nH):
                H[k][i, j] = complex(H_fl[k][idx], H_fl[k][idx+1])


    obs = 0
    for i in range(nH):
        temp = numpy.real(numpy.trace(H[i]))
        w_TR.write(str(temp) + " ")
        obs += temp*temp

    w_TR.write("\n")
    w_F.write(str(obs) + "\n")

    
r_HL.close()
w_F.close()
w_TR.close()
    





