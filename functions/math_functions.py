import numpy as np
import math

#tried to use this initially but failed at k = 40 due to overflow
#tried return log(product) but realised that was a recursive log return so...
def factorial(x:np.int32
              ) -> np.float32:

    if x > 1:
        return np.float32(x * factorial(x-1))
    else:
        return np.float32(1)

#solved overflow issue
def log_factorial(x:np.int32
                  ) -> np.float32:

    fac_sum = 0.0
    for i in range(1, x + 1):
        fac_sum += math.log(i)

    return np.float32(fac_sum)

def poisson_distribution(lam:np.float32|np.int32,
                         k:np.int32
                         ) -> np.int32:

    if k == 0:
        p = ((lam ** k) * np.exp(-lam)) / (factorial(k))
        return np.float32(p)
    else:
        p = k * math.log(lam) - lam - log_factorial(k)
        return np.exp(np.float32(p))

def get_vandermonde_matrix(x,y):
    #pretty confident this is correct
    #going to have to create the vandermonde matrix from the xi as j is treated as an exponent here.
    #Vandermonde matrix has to be square so just creating a square matrix of 0s
    V_m = np.zeros((len(x),len(x)))
    for j in range(V_m.shape[1]):
        for i in range(V_m.shape[0]):
            V_m[i,j] = x[i]**j

    return V_m