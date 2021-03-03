import numpy as np
from math import *
from random import randint

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

 # A = B (généré aléat) + Bt ==> A symétrique
def symetrik_generator(N):
    B = np.zeros(shape=(3, 3))
    Bt = B.T
    lines = B.shape[0]
    columns = B.shape[1]

        # Erreur dans la value N
    if (N < 0): return 0

    nb_zero = (lines * (columns-1)/2) - N

        # Rempli la matrice avec valeurs aléatoires sauf cas particulier éléments diagonaux
    for i in range(0, columns):
        sum = 0
        for j in range(0, lines):
                # On gère les éléments diagonaux à la fin
            if (i == j):
                continue

            value = randint(-100, 100)
            sum += abs(value)
            B[i][j] = value

        B[i][i] = randint(sum+1, sum+100)

    A = B + Bt

        # Annule un maximum d'éléments non diagonaux en fonction du paramètre N
    for i in range(0, columns):
        for j in range(0, lines):
            if (i == j):
                continue

            if (nb_zero > 0):
                A[i][j] = 0
                A[j][i] = 0
                nb_zero -= 1
            else:
                break

    print(A)
    return A

symetrik_generator()


