# Modules nécessaires
import matplotlib.pyplot as plt
import numpy as np
import math as m

### Question 1 : Par quels aspects cette implémentation ne respecte t’elle pas des standards de codage très sains ?
"""
Ap = A*p
variables n'ont pas des noms faciles
10**6 et 10**-10 en dur ?
"""

### Question 2 : Par quel aspect cette méthode est-elle différente de la méthode décrite mathématiquement juste avant ?
# Quel est le gain de complexité occasionné ? Est-il systématique ?
"""

"""

### Question 3 : Implémentation de la méthode du gradient conjugué



"""
    @brief Cette fonction met en oeuvre la méthode du Gradient Conjugué.
    Résolution de Ax = b, avec x supposée inconnue.

    @param A : Matrice numpy symétrique définie positive.
    @param b : Vecteur numpy de même hauteur que A.
    @param x : Vecteur numpy représentant l'origine de la descente de gradient. 
    @param imax : Le nombre d'itération maximale. 
    @param imax : La précision servant à stopper la boucle d'itération. 
    @return Le vecteur x solution de l'équation Ax = b avec la précision spécifiée.
"""
def conjugateGradient(A, b, x, imax=10**6, precision=1e-10):

    # Initialisation de r0 et p0
    r = b - np.dot(A, x)
    p = r
    rsquare_old = np.dot(r, np.transpose(r))

    for i in range(imax): #Nombre d'itérations arbitraire

        # Calcul de alpha_k+1
        Ap = np.dot(A, p)
        alpha = rsquare_old / np.dot(np.transpose(p), Ap)
        # Mise à jour de x et calcul de r_k+1
        x += alpha * p
        r -= alpha * Ap
        rsquare_new = np.dot(r, np.transpose(r))

        # On a atteint la précision attendue
        if(m.sqrt(rsquare_new) < precision): return x

        # Sinon, calcul de p_k+1 et passage à l'itération suivante
        p = r + p * (rsquare_new/rsquare_old)
        rsquare_old = rsquare_new

    # La précision n'a pas été atteinte, on renvoie quand même le résultat
    return x


def print_res(msg,res):
    if res :
        print("[SUCCESS] --> ",end='')
    else :
        print("[FAILED] --> ",end='')
    print(msg)

# Tests de comportement
matrice_test_1 = np.array([[4., -2., -4.],
                           [-2., 10., 5.],
                           [-4., 5.,  6.]])
matrice_test_2 = np.array([[4., 2., -9.],
                           [2., 10., -5.],
                           [-9., -5.,  6.]])
vector_test_1 = np.array([1., 2., 4.])
vector_test_2 = np.array([1., 4., 8.])
zeros_1 = np.zeros(3)
zeros_2 = np.zeros(3)

res1 = conjugateGradient(matrice_test_1, vector_test_1, zeros_1)
res2 = conjugateGradient(matrice_test_2, vector_test_2, zeros_2)
print(res1)
print(res2)
print_res("Conjugate gradient 1, is result close to [ 3.861, -1.111,  4.167].", 
          np.around(res1[0],3) == 3.861 and np.around(res1[1],3) == -1.111 and np.around(res1[2],3) == 4.167)
print_res("Conjugate gradient 2, is result close to [-1.570,  0.348, -0.732].", 
          np.around(res2[0],3) == -1.570 and np.around(res2[1],3) == 0.348 and np.around(res2[2],3) == -0.732)

### Question 4 : Implémentation de la méthode du gradient conjugué avec préconditionneur

# Renvoie la solution de l'équation Ax = b sous forme de vecteur
def conjugateGradientPrecond(A, b, x):

    # Initialisation de r0, z0 et p0
    r = b - np.dot(A, x)
    z = np.dot(np.invert(M), r) #TODO retrouver ce qu'est M (paramètre d'entrée non ?)
    p = z
    rsquare_old = np.dot(r, np.transpose(z))

    for i in range(10**6): #Nombre d'itérations arbitraire

        #TODO à partir de là
        # Calcul de alpha_k+1
        Ap = np.dot(A, p)
        alpha = rsquare_old / np.dot(np.transpose(p), Ap)
        # Mise à jour de x et calcul de r_k+1
        x += alpha * p
        r -= alpha * Ap
        rsquare_new = np.dot(r, np.transpose(r))

        # On a atteint la précision attendue
        if(m.sqrt(rsquare_new) < 1e-10): return x

        # Sinon, calcul de p_k+1 et passage à l'itération suivante
        p = r + p * (rsquare_new/rsquare_old)
        rsquare_old = rsquare_new

    # La précision n'a pas été atteinte, on renvoie quand même le résultat
    return x


















