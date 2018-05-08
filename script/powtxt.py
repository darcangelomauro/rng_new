
nH = int(input('nH\n'))
nL = int(input('nL\n'))

for i in range(nH):
    print('TrH' + str(i+1) + ' TrH' + str(i+1) + '^2 Tr%H' + str(i+1) + '^2', end = ' ')
for i in range(nH, nH+nL):
    print('TrL' + str(i+1-nH) + ' TrL' + str(i+1-nH) + '^2 Tr%L' + str(i+1-nH) + '^2', end = ' ')
print('')
