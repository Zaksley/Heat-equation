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

                # Cas diagonale
            if (i == j):
                sum = 0
                for k in range(0, i):
                    sum += L[i][k]**2
                value = sqrt(A[j][i] - sum)

            else:
                sum = 0
                for k in range(0, i):
                    sum += L[j][k] * L[i][k]
                value = (A[j][i] - sum) / L[i][i]
                
            L[j][i] = value

    return L

def test_facto_cholesky():

    #Test 1
    A = np.array( [ [4, -2, -4], [-2, 10, 5], [-4, 5, 6] ])
    L = facto_cholesky(A)
    Lt = L.T
    test_A = np.dot(L, Lt)

    if ( (A == test_A).all()):
        return True
    return False
    

#L = facto_cholesky(a)
#print(test_facto_cholesky())

"""
1) Complexité n^3
"""

def symetrik_generator():
    print("test")

