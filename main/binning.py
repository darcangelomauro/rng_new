import numpy

code = input("Code\n")
filenameERR = code + "_simERR.txt"
filenameF = code + "_simF.txt"

r_F = open(filenameF, "r")
w_ERR = open(filenameERR, "w")


# read observable
F = []
while True:
    temp_str = r_F.readline()
    if not temp_str:
        break
    F.append(float(temp_str))


for k in range(1, 10000):
    if not k%1000:
        print(k)
    bin_F = numpy.array_split(F, k)
    bin_F = [numpy.mean(item) for item in bin_F]

    #mean = numpy.mean(bin_F)
    var = numpy.var(bin_F)
    
    w_ERR.write(str(k) + " " + str(var) + "\n")

        
r_F.close()
w_ERR.close()
    





