import math 
import scipy.special

# convert 'num' to base 'base' and add zeroes to reach length 'length' 
def convert_base(num, base, length):
    ls = []
    while num != 0:
        ls.append(num % base)
        num = int(num/base)
    while len(ls) < length:
        ls.append(0)
    return ls

# return first cyclic permutation among all cyclic permutations 
def first_cyclic(ls):

    n = len(ls)
    b = [[ls[i - j] for i in range(n)] for j in range(n)]
    b.sort()

    return b[0]

def sort_lists(ls1, ls2):
    ls = []
    ls.append("_".join(ls1))
    ls.append("_".join(ls2))
    ls.sort()
    ls1_new = ls[0].split('_')
    
    if ls1_new == ls1:
        return 0
    else:
        return 1


n = input('input n: ')
n = int(n)

for s in range(1, n+1):
    terms = []
    print('s: ' + str(s))
    print('n-s: ' + str(n-s))

    # generate i indices
    i = []
    for num_i in range(n-s):
        i.append('i[' + str(num_i) + ']')
        
    # generate k indices and check that their sum is <= s 
    for num_k in range(int(math.pow(s+1, n-s))):
        k = convert_base(num_k, s+1, n-s)

        sum_k = 0
        for num in k:
            sum_k += num

        if sum_k <= s:
            
            # generate p indices and check that p[j] <= k[j] for all j
            for num_p in range(int(math.pow(s+1, n-s))):
                p = convert_base(num_p, s+1, n-s)

                check_p = 1
                for idx in range(len(p)):
                    if p[idx] > k[idx]:
                        check_p = 0
                        break

                if check_p:

                    # loop over l index
                    for l in range(s - sum_k + 1):

                        # generate q indices
                        for num_q in range(int(math.pow(2, n-s))):
                            q = convert_base(num_q, 2, n-s)

                            
                            # NOW THE INDICES ARE ALL SET

                            # 1. gamma matrices
                            first = []

                            for idx_1 in range(n-s):
                                for idx_k in range(k[idx_1]):
                                    first.append('gamma[uM]')
                                first.append('gamma[' + i[idx_1] + ']')
                            for idx_1 in range(s - sum_k):
                                first.append('gamma[uM]')

                            first = first_cyclic(first)

                            #print(first)


                            # 2. normal matrices
                            second = []
                            
                            for idx_2 in range(n-s):
                                for idx_p in range(p[idx_2]):
                                    second.append('dM')
                                for idx_q in range(q[idx_2]):
                                    second.append('MAT[' + i[idx_2] + ']')
                            for idx_2 in range(l):
                                second.append('dM')

                            if len(second):
                                second = first_cyclic(second)
                            else:
                                second.append('I')

                            #print(second)


                            # 3. transposed matrices
                            third = []
                            
                            for idx_3 in range(s - sum_k - l):
                                third.append('dM')
                            for idx_3 in range(n-s):
                                for idx_1q in range(1-q[idx_3]):
                                    third.append('MAT[' + i[idx_3] + ']')
                                for idx_pk in range(k[idx_3] - p[idx_3]):
                                    third.append('dM')

                            if len(third):
                                third = first_cyclic(third)
                            else:
                                third.append('I')

                            #print(third)

                            # order second and third
                            if sort_lists(second, third):
                                second, third = third, second


                            # 4. binomial coefficients
                            binom = 1

                            for idx_4 in range(n-s):
                                binom *= scipy.special.comb(k[idx_4], p[idx_4], exact=True)
                            binom *= scipy.special.comb(s-sum_k, l, exact=True)

                            #print(binom)


                            # 5. e factors
                            e = []

                            sum_p = 0
                            for num in p:   
                                sum_p += num

                            if (s - sum_p - l) % 2 != 0:
                                e.append('e[uM]')
                            
                            for idx_5 in range(n-s):
                                if not q[idx_5]:
                                    e.append('e[' + i[idx_5] + ']')

                            #print(e)


                            # 6. simplify gammas (gamma[i] squares to e[i]*identity)
                            check_6 = 1
                            while len(first) > 1 and check_6:
                                check_6 = 0
                                for idx_6 in range(len(first)-1):
                                    if first[idx_6] == first[idx_6+1]:
                                        temp = first.pop(idx_6)
                                        del first[idx_6]
                                        e.append(temp.replace('gamma', 'e'))
                                        if not len(first):
                                            first.append('I_g')
                                        if len(first) > 1:
                                            check_6 = 1
                                        break
                            first = first_cyclic(first)


                            # 7. simplify e (e[i] always squares to 1)
                            for item in e:
                                occ = e.count(item)
                                e = [value for value in e if value != item]
                                if occ % 2:
                                    e.append(item)


                            # 8. 





                            # add everything to the full lists
                            terms.append([[binom] + e] + [first] + [second] + [third])
    
    for item in terms:
        print(item)








