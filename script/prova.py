import itertools

def generate_combinations(a, n):
    A = []
    for i in range(n):
        A.append(a)
    res1 = itertools.product(*A)
    res2 = [''.join(list(item)) for item in res1]
    return res2


a = ['0', '1', '2', '3', '4', '5']

ls = generate_combinations(a, 4)

idx_ls = [5, 2, 3, 5]

idx = 0

for num in range(len(idx_ls)-1):
    idx += idx_ls[num]
    idx *= 6
idx += idx_ls[-1]

print(ls[idx])

print(type(idx_ls))

