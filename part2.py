# Modules necessaires
import matplotlib.pyplot as plt
import numpy as np
import math as m

### Question 1 : Par quels aspects cette implementation ne respecte t’elle pas des standards de codage tres sains ?
"""
Ap = A*p
variables n'ont pas des noms faciles
10**6 et 10**-10 en dur ?
"""

### Question 2 : Par quel aspect cette methode est-elle differente de la methode decrite mathematiquement juste avant ?
# Quel est le gain de complexite occasionne ? Est-il systematique ?
"""

"""

### Question 3 : Implementation de la methode du gradient conjugue


"""
    @brief Cette fonction met en oeuvre la methode du Gradient Conjugue.
    Resolution de Ax = b, avec x supposee inconnue.

    @param A : Matrice numpy symetrique definie positive.
    @param b : Vecteur numpy de même hauteur que A.
    @param x : Vecteur numpy representant l'origine de la descente de gradient.
    @param imax : Le nombre d'iteration maximale.
    @param imax : La precision servant a stopper la boucle d'iteration.
    @return Le vecteur x solution de l'equation Ax = b avec la precision specifiee.
"""
def conjugateGradient(A, b, x, imax=10**6, precision=1e-10):

    # Initialisation de r0 et p0
    r = b - np.dot(A, x)
    p = r
    rsquare_old = np.dot(r, r.T)
    
    for i in range(imax): #Nombre d'iterations arbitraire

        # Calcul de alpha_k+1
        Ap = np.dot(A, p)
        alpha = rsquare_old / np.dot(p.T, Ap)
        # Mise a jour de x et calcul de r_k+1
        x += np.dot(alpha, p)
        r -= np.dot(alpha, Ap)
        
        rsquare_new = np.dot(r.T, r)

        # On a atteint la precision attendue
        if(m.sqrt(rsquare_new) < precision): 
            return x

        # Sinon, calcul de p_k+1 et passage a l'iteration suivante
        p = r + p * (rsquare_new/rsquare_old)
        rsquare_old = rsquare_new

    # La precision n'a pas ete atteinte, on renvoie quand même le resultat
    return x

def conjugateGradient2(A, b, x_0, imax=10**6, precision=1e-10):
    r = b - np.dot(A,x_0)
    p = r
    x = x_0
    for i in range(imax): #Nombre d'iterations arbitraire
        a_k = (np.dot(r.T,r))/(np.dot(p.T,np.dot(A,p)))
        x_next = x + (a_k * p)
        r_next = r - (a_k * np.dot(A,p))
        if np.linalg.norm(r_next) < precision: 
            return x_next
        b_k = (np.dot(r_next.T,r_next))/(np.dot(r.T,r))
        p = r_next + (b_k * p)
        #passage au k+1
        r = r_next
        x = x_next
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

def tests():
    res1 = conjugateGradient(matrice_test_1, vector_test_1, zeros_1)
    res2 = conjugateGradient(matrice_test_2, vector_test_2, zeros_2)
    print(res1)
    print(res2)
    print_res("Conjugate gradient 1, is result close to [ 3.861, -1.111,  4.167].", 
            np.around(res1[0],3) == 3.861 and np.around(res1[1],3) == -1.111 and np.around(res1[2],3) == 4.167)
    print_res("Conjugate gradient 2, is result close to [-1.570,  0.348, -0.732].", 
            np.around(res2[0],3) == -1.570 and np.around(res2[1],3) == 0.348 and np.around(res2[2],3) == -0.732)

    from part1 import symetrik_generator
    import random
    ReelA = symetrik_generator(random.randint(0, 10),10)
    Reelx = np.array([random.uniform(0, 10) for iter in range(10)])
    Reelx = Reelx.reshape((10,1))
    Reelb = np.dot(ReelA,Reelx)
    tmp_x = np.zeros((10,1))
    res_x = conjugateGradient2(ReelA,Reelb,tmp_x)
    res_allclose = np.allclose(np.dot(ReelA, res_x), Reelb)
    print_res("Test matrice size 10, with numpy vectors.", res_allclose)


### Question 4 : Implementation de la methode du gradient conjugue avec preconditionneur


"""
    @brief Cette fonction implemente le Pivot de Gauss une résolution de systeme lineaire.
    Resolution de Tx = b, avec x supposee inconnue.

    @param T : Matrice numpy triangulaire.
    @param b : Vecteur numpy de même hauteur que T.
    @param sup : 1 si T est superieure, 0 si elle est inferieure
    @return Le vecteur x solution de l'equation Tx = b.
"""
def pivotGauss(T,b,sup=1):
    n = np.size(b)
    x = np.zeros((n,1))

    if(sup == 1): #T est triangulaire superieure
        for i in range(n-1, -1, -1):
            tx = 0

            for j in range(i+1, n):
                tx +=  T[i][j]*x[j]

            x[i] = (b[i] - tx) / T[i][i]

    else: #T est triangulaire inferieure
        for i in range(n):
            tx = 0

            for j in range(i):
                tx += T[i][j]*x[j]

            x[i] = (b[i] - tx) / T[i][i]

    return x

"""
    @brief Cette fonction met en oeuvre la methode du Gradient Conjugue avec preconditionneur.
    Resolution de Ax = b, avec x supposee inconnue et en utilisant le preconditionneur M.

    @param A : Matrice numpy symetrique definie positive.
    @param T : Matrice numpy triangulaire superieure (issue d'une decomposition de Cholesky).
    @param b : Vecteur numpy de même hauteur que A.
    @param x : Vecteur numpy representant l'origine de la descente de gradient.
    @param imax : Le nombre d'iteration maximale.
    @param imax : La precision servant a stopper la boucle d'iteration.
    @return Le vecteur x solution de l'equation Ax = b avec la precision specifiee.
"""
def conjugateGradientPrecondIter(A, T, b, x, imax=10**6, precision=1e-10):
    # Initialisation de r0, z0 et p0
    r = b - np.dot(A, x)
    # On utilise le pivot de Gauss pour ne pas avoir a calculer M^-1
    # On passe pour l'instant T en parametre, mais on peut le changer pour M avec les fonctions de la partie 1
    z = pivotGauss(T, r, 1) # Tz = r
    z = pivotGauss(np.transpose(T), z, 0) #T'z = z
    p = z
    res = []
    for i in range(imax): #Nombre d'iterations arbitraire
        a_k = np.dot(r.T,z)/(np.dot(p.T,np.dot(A,p)))
        x_next = x + a_k * p
        r_next = r - a_k * np.dot(A,p)
        # Si on a atteint la precision attendue
        res.append(x_next)
        if(np.linalg.norm(r_next) < precision): 
            return x_next
        # Sinon, calcul de z_k+1, p_k+1 et passage a l'iteration suivante
        z_next = pivotGauss(T, r_next, 1) 
        z_next = pivotGauss(np.transpose(T), z_next, 0)
        b_k = np.dot(r_next.T,z_next)/(np.dot(r_next.T,r_next))
        p = z_next + b_k * p 
        x = x_next
        r = r_next
        z = z_next
    # La precision n'a pas ete atteinte, on renvoie quand même le resultat
    return x

triangulaire_test = np.array([[6., -2., -1.],
                              [0., 7., 2.],
                              [0., 0., 9.]])
#Decomposition de cholesky de array([[ 41., -16.,  -9.],   symetrique definie positive (a diagonale dominante)
#                                    [-16.,  53.,  18.],
#                                    [ -9.,  18.,  81.]])

# res3 = conjugateGradientPrecond(matrice_test_1, triangulaire_test, vector_test_1, zeros_1)
# res4 = conjugateGradientPrecond(matrice_test_2, triangulaire_test, vector_test_2, zeros_2)
# print(res3)
# print(res4)
# print_res("Conjugate gradient (Precond) 1, is close to [ 3.861, -1.111,  4.167].",
#          np.around(res3[0],3) == 3.861 and np.around(res3[1],3) == -1.111 and np.around(res3[2],3) == 4.167)
# print_res("Conjugate gradient (Precond) 2, is close to [-1.570,  0.348, -0.732].",
#          np.around(res4[0],3) == -1.570 and np.around(res4[1],3) == 0.348 and np.around(res4[2],3) == -0.732)













