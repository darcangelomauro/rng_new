
nH = int(input('nH\n'))
nL = int(input('nL\n'))

print('COMM1')
for i in range(nH):
    for j in range(i+1, nH):
        for k in range(nH, nH+nL):
            print('[H' + str(i+1) + ',H' + str(j+1) + ']L' + str(k-nH+1), end=' ')
print('')


print('COMM2')
for i in range(nH, nH+nL):
    for j in range(i+1, nH+nL):
        for k in range(nH, nH+nL):
            print('[L' + str(i-nH+1) + ',L' + str(j-nH+1) + ']L' + str(k-nH+1), end=' ')
print('')


print('COMM3')
for i in range(nH, nH+nL):
    for j in range(nH):
        for k in range(nH):
            print('[L' + str(i-nH+1) + ',H' + str(j+1) + ']H' + str(k+1), end=' ')
print('')
