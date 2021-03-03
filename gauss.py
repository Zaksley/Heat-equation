# Modules necessaires
import matplotlib.pyplot as plt
import numpy as np
import math as m

"""
    @brief Cette fonction implémente le Pivot de Gauss une résolution de système linéaire.
    Resolution de Tx = b, avec x supposee inconnue.

    @param T : Matrice numpy triangulaire supérieure.
    @param b : Vecteur numpy de même hauteur que T.
    @return Le vecteur x solution de l'equation Tx = b.
"""
def pivotGauss(T,b):
    n = np.size(b)
    x = np.zeros(n)

    for i in range(n-1, -1, -1):
        tx = 0

        for j in range(i+1, n):
            tx +=  T[i][j]*x[j]

        x[i] = (b[i] - tx) / T[i][i]

    return x



def print_res(msg,res):
    if res :
        print("[SUCCESS] --> ",end='')
    else :
        print("[FAILED] --> ",end='')
    print(msg)

triangulaire_test_1 = np.array([[2., -2., -1.],
                                [0., 2., 5.],
                                [0., 0., 1.]])
triangulaire_test_2 = np.array([[1., 3., -1.],
                                [0., 1., 2.],
                                [0., 0.,  1.]])
vector_test_1 = np.array([1., 2., 4.])
vector_test_2 = np.array([3., 1., 1.])

res1 = pivotGauss(triangulaire_test_1, vector_test_1)
res2 = pivotGauss(triangulaire_test_2, vector_test_2)
print(res1)
print(res2)

print_res("Gaussian elimination results in [-6.5, -9,  4].",
          np.around(res1[0],3) == -6.5 and np.around(res1[1],3) == -9 and np.around(res1[2],3) == 4)
print_res("Gaussian elimination results in [7, -1,  1].",
          np.around(res2[0],3) == 7 and np.around(res2[1],3) == -1 and np.around(res2[2],3) == 1)