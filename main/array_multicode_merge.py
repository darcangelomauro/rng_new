import math

# input multicode
multicode = input('Code\n')
array_size = int(input('Array size\n'))
array_first = int(input('First array index\n'))

filenames = []
for idx in range(array_first, array_first+array_size):
    filename_args = multicode + str(idx) + '_varG_args.txt'
    try:
        with open(filename_args) as r_args:
            args_list = r_args.readline().split(' ')
        args_list = [item for item in args_list if len(item) >= 2]
        filenames.append(args_list)
    except:
        print('Missing ' + str(idx) + '. Ignoring.')


N = len(filenames)
print(N)
for i, item in enumerate(filenames):
    print(str(i) + ' ' + str(len(item)))
filenames = [list(x) for x in zip(*filenames)]

with open(multicode + '_varG_args.txt', 'w') as w_args, open(multicode + '_varG_G_args.txt', 'w') as w_G_args:
    for Gval in filenames:
        
        # ugly thing to build name
        name_list = Gval[0].split('_')
        code = name_list[-1]

        # concatenate array simS and simHL
        with open(multicode + '_' + code + '_simS.txt', 'w') as outfile:
            for i, fname in enumerate(Gval):
                with open(fname + '_simS.txt') as infile:
                    for line in infile:
                        outfile.write(line)
        with open(multicode + '_' + code + '_simHL.txt', 'w') as outfile:
            for i, fname in enumerate(Gval):
                with open(fname + '_simHL.txt') as infile:
                    for line in infile:
                        outfile.write(line)

        # create a data file
        with open(multicode + '_' + code + '_data.txt', 'w') as outfile, open(multicode + str(array_first) + '_' + code + '_data.txt') as infile:
            line_list = (infile.readline()).split(' ')
            G = float(line_list[4])
            nsw = int(line_list[5])
            gap = int(line_list[6])
            nmeas = int(math.ceil(nsw/gap))
            nmeas *= N
            nsw1 = nmeas*gap
            line_list[5] = str(nsw1)

            outfile.write(' '.join(line_list))

        # create varG_args and varG_G_args files
        w_args.write(multicode + '_' + code + ' ')
        w_G_args.write(str(G) + ' ' + multicode + '_' + code + '\n')

# create varG_data file
with open(multicode + '_varG_data.txt', 'w') as outfile:
    with open(multicode + str(array_first) + '_varG_data.txt') as infile:
        for line in infile:
            if line.find('Nsw') != -1:
                outfile.write('Nsw, GAP and ADJ refer to a single job in the array' + '\n')
                outfile.write(line)
            elif line.find('ADJ') != -1:
                outfile.write(line)
                outfile.write('Job array size: ' + str(N) + '\n')
                outfile.write('Nmeas: ' + str(nmeas) + '\n')
            else:
                outfile.write(line)


