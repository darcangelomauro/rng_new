import numpy

code = input("Code\n")
filenameF = code + "_simF.txt"
filenameA = code + "_autF.txt"

r_F = open(filenameF, "r")
w_A = open(filenameA, "w")

obs = []


while True:
    


    obs_str = r_F.readline()
    if not obs_str:
        break
    obs.append(float(obs_str))

r_F.close()

n = len(obs)
mean = numpy.mean(obs)

gamma0 = 0
for j in range(n):
    gamma0 += (obs[j]-mean)*(obs[j]-mean)
w_A.write(str(0) + " " + str(1) + "\n")

for i in range(1,10000):
    if not i%1000:
        print(i)
    gamma = 0
    for j in range(n-i):
        gamma += (obs[j]-mean)*(obs[j+i]-mean)
    gamma /= gamma0
    #gamma /= float(n-i)

    w_A.write(str(i) + " " + str(gamma) + "\n")


w_A.close()
    
