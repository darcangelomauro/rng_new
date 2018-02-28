
world_size = int(input("Input world size\n"))

first_rank = int(input("Input first rank\n"))

r_init_fileobj = open("master_init.txt", "r")
r_gval_fileobj = open("master_gval.txt", "r")

for i in range(first_rank, first_rank+world_size):
    
    w_init_filename = "init" + str(i) + ".txt"
    w_gval_filename = "gval" + str(i) + ".txt"

    w_init_fileobj = open(w_init_filename, "w")
    w_gval_fileobj = open(w_gval_filename, "w")

    for line in r_init_fileobj:
        w_init_fileobj.write(line + "\n")
    for line in r_gval_fileobj:
        w_gval_fileobj.write(line + "\n")

    w_init_fileobj.close()
    w_gval_fileobj.close()

r_init_fileobj.close()
r_gval_fileobj.close()
