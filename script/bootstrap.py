import math
import numpy



def bootstrap(ls, nresamp):

    N = len(ls)

    vals = []

    # repeat resampling 'nresamp' times
    for i in range(nresamp):
        
        # list of resampled data
        new_ls = []
        
        # build new list of N elements with repetitions 
        for j in range(N):
            # pick random index
            idx = numpy.random.randint(N)
            # append that element in new list
            new_ls.append(ls[idx])

        vals.append(numpy.mean(new_ls))

    return math.sqrt(numpy.var(vals))



