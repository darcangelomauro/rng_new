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


def contract_list(ls):
    ls_new = []
    for item in ls:
        check_1 = 1
        for item_new in ls_new:
            if len(item) == len(item_new):
                check_2 = 1
                for idx in range(1, len(item)):
                    if item[idx] != item_new[idx]:
                        check_2 = 0
                if check_2:
                    item_new[0] = item_new[0] + '+' + item[0]
                    check_1 = 0

        if check_1:
            ls_new.append(item)

    for item in ls_new:
        item[0] = [item[0]]

    return ls_new

n = input('input n: ')
n = int(n)


print('gsl_complex pfac;')
print('int jfac;')
print('gsl_complex res = GSL_COMPLEX_ZERO;')

for s in range(1, n+1):
    print('// CASE: s=' + str(s))
    terms = []
    
    '''
    print('s: ' + str(s))
    print('n-s: ' + str(n-s))
    '''

    # generate i indices
    i = []
    if n-s:
        print('int* i' + str(s) + ' = malloc(' + str(n-s) + '*sizeof(int));')
    for num_i in range(n-s):
        i.append('i' + str(s) + '[' + str(num_i) + ']')
        print('for(' + i[num_i] + '=0; ' + i[num_i] + '<nHL; ' + i[num_i] + '++)')
        print('{')

        
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

                            if not len(e):
                                e.append('1')

                            #print(e)

                            # 5.1 i factor
                            temp_51 = []
                            for idx_51 in range(n-s):
                                temp_51.append('e[i[' + str(idx_51) + ']]')
                            if not len(temp_51):
                                temp_51.append('0')
                            j = ['0.5*(' + str(n) + '-' + str(s) + '*e[uM]' + '-' + '-'.join(temp_51) + ')']



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
                                elif not len(e):
                                    e.append('1')


                            # add everything to the full lists
                            terms.append([str(binom) + '*' + '*'.join(e)] + [j] + [first] + [second] + [third])

    
    # 8. group terms together if they are the same
    terms = contract_list(terms)

    # 9. dim factors
    for term in terms:
        I_count = term.count(['I'])
        I_g_count = term.count(['I_g'])

        if I_count:
            term[0][0] = 'dim*(' + term[0][0] + ')'
            I_count -= 1
            while I_count:
                term[0][0] = 'dim*' + term[0][0]
                I_count -= 1

            while I_g_count:
                term[0][0] = 'dimG*' + term[0][0]
                I_g_count -= 1
        else:
            if I_g_count:
                term[0][0] = 'dimG*(' + term[0][0] + ')'
                I_g_count -= 1
                while I_g_count:
                    term[0][0] = 'dimG*' + term[0][0]
                    I_g_count -= 1
        
    terms = [[item for item in term if item != ['I'] and item != ['I_g']] for term in terms]


                


    # 10. find needed matrix products and traces

    products = []

    for term in terms:
        for idx_10 in range(2, len(term)):
            check_10 = 1
            for product in products:
                if product == term[idx_10]:
                    check_10 = 0
            if check_10:
                products.append(term[idx_10])
    
    ''' 
    for item in terms:
        print(item)
    print('')
    for item in products:
        print(item)
    '''

    # 11. write code to declare traces
    for item in products:
        if len(item) > 1:
            name = 'array' + str(s) + str(products.index(item))
            print('gsl_matrix_complex** ' + name + ' = malloc(' + str(len(item)) + '*sizeof(gsl_matrix_complex*));')
            for idx_11 in range(len(item)):
                print(name + '[' + str(idx_11) + '] = ' + item[idx_11] + ';')

            res = ''.join(item)
            if res.find('gamma') != -1:
                dim = 'dimG'
            else:
                dim = 'dim'
            print('gsl_matrix_complex* ' + res + ' = gsl_matrix_complex_alloc(' + dim + ', ' + dim + ');')
            print('matrix_multiprod(' + name + ', ' + str(len(item)) + ', ' + res + ');')
            print('gsl_complex tr' + res + ' = trace(' + res + ');')
            print('gsl_matrix_complex_free(' + res + ');')
            print('free(' + name + ');')
        else:
            print('gsl_complex tr' + ''.join(item) + ' = trace(' + ''.join(item) + ');')

    
    # 12. write code to compute action
    for term in terms:
        
        pre = term[0][0]
        jfac = term[1][0]
        print('jfac = (int)' + jfac + ' % 4;')
        print('if(jfac == 0)')
        print('pfac = gsl_complex_rect(' + pre + ', 0);')
        print('if(jfac == 1)')
        print('pfac = gsl_complex_rect(0, ' + pre + ');')
        print('if(jfac == 2)')
        print('pfac = gsl_complex_rect(-' + pre + ', 0);')
        print('if(jfac == 3)')
        print('pfac = gsl_complex_rect(0, -' + pre + ');')
        compl = []
        for idx_12 in range(2, len(term)):
            compl.append('tr' + ''.join(term[idx_12]))
        
        print('res = gsl_complex_add(res, ', end='') 
        if len(compl) == 3:
            print('gsl_complex_mul(gsl_complex_mul(gsl_complex_mul(' + compl[0] + ', ' + compl[1] + '), ' + compl[2] + '), pfac)', end='')
        elif len(compl) == 2:
            print('gsl_complex_mul(gsl_complex_mul(' + compl[0] + ', ' + compl[1] + '), pfac)', end='')
        elif len(compl) == 1:
            print('gsl_complex_mul(' + compl[0] + ', pfac)', end='')
        print(');')

    for val in range(n-s):
        print('}')
    
    if n-s:
        print('free(i' + str(s) + ');')



        

















