import numpy as np
import time
from math import *
from random import randint


#a = np.array( [ [4, -2, -4], [-2, 10, 5], [-4, 5, 6] ])
#print(a)

def facto_cholesky(A):
        #Initialize L
    L = np.zeros(shape=A.shape)

    lines = A.shape[0]
    columns = A.shape[1]
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
    B = np.zeros(shape=(4, 4))
    Bt = B.T
    lines = B.shape[0]
    columns = B.shape[1]

        # Erreur dans la value N
    if (N < 0): return 0

    nb_zero = ( (lines-1) * columns / 2 ) - N

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
                if (A[i][j] != 0 and A[j][i] != 0):
                    A[i][j] = 0
                    A[j][i] = 0
                    nb_zero -= 1
            else:
                break
    return A

"""
def test_symetrik_generator():
    for i in range(20):
        print(symetrik_generator(3))
"""

def facto_cholesky_incomplete(A):
        #Initialize L
    L = np.zeros(shape=A.shape)

    lines = A.shape[0]
    columns = A.shape[1]
    for i in range(0, columns):
        for j in range(i, lines):

                # Version incomplete de cholesky
            # => On ne traite plus les éléments nuls 
            if (A[i][j] == 0):
                L[i][j] = 0
                continue


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


def test_facto_cholesky_incomplete():

    nb_tests = 100
    nb_true = 0
    # Batterie de test sur 100 matrices générés
    for z in range(nb_tests):

            # N dépends de la taille des matrices
        # => Attention donc à la shape fourni par symetrik_generator()
        A = symetrik_generator(randint(0, 3))
        L = facto_cholesky_incomplete(A)
        Lt = L.T
        test_A = np.dot(L, Lt)

        if ( np.allclose(A, test_A)):
            nb_true += 1
    
    
    return (nb_true == nb_tests)


#print(test_facto_cholesky_incomplete())


def test_time_cholesky():

    nb_tests = 1000
    tps1 = time.time()
    
    # Batterie de test sur 100 matrices générés - CHOLESKY NORMAL
    for z in range(nb_tests):
        A = symetrik_generator(0)
        L = facto_cholesky(A)

    tps2 = time.time()
    time_cholesky = (tps2 - tps1)

    tps1 = time.time()
    # Batterie de test sur 100 matrices générés - CHOLESKY INCOMPLETE
    for z in range(nb_tests):
        A = symetrik_generator(4)
        L = facto_cholesky_incomplete(A)

    tps2 = time.time()
    time_cholesky_incomplete = (tps2 - tps1)

    print(time_cholesky)
    print(time_cholesky_incomplete)


test_time_cholesky()


