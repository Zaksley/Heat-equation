import numpy as np
from math import *

a = np.array( [ [4, -2, -4], [-2, 10, 5], [-4, 5, 6] ])

#print(a)

def facto_cholesky(A):
        #Initialize L
    L = np.zeros(shape=a.shape)

    lines = a.shape[0]
    columns = a.shape[1]
    for i in range(0, columns):
        for j in range(i, lines):

            sum = 0
            value = 0

                # Cas diagonale
            if (i == j):
                
                for k in range(0, i-1):
                    sum += L[i][k]**2
                value = sqrt(A[j][i] - sum)

            else:
                for k in range(0, i-1):
                    sum += L[j][k] * L[i][k]
                value = (A[j][i] - sum) / L[i][i]
                
            L[j][i] = value

    return L

L = facto_cholesky(a)
print(L)
