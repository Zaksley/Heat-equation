import matplotlib.pyplot as plt
import numpy as np 
import random
import math as m

from part2 import conjugateGradient2,pivotGauss
from part1 import symetrik_generator

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
            return res
        else:
            print(np.linalg.norm(r_next))
        # Sinon, calcul de z_k+1, p_k+1 et passage a l'iteration suivante
        z_next = pivotGauss(T, r_next, 1) 
        z_next = pivotGauss(np.transpose(T), z_next, 0)
        b_k = np.dot(r_next.T,z_next)/(np.dot(r_next.T,r_next))
        p = z_next + b_k * p 
        x = x_next
        r = r_next
        z = z_next
    # La precision n'a pas ete atteinte, on renvoie quand même le resultat
    return res



def conjugateGradient2Iter(A, b, x_0, imax=10**6, precision=1e-10):
    r = b - np.dot(A,x_0)
    p = r
    x = x_0
    res = []
    for i in range(imax): #Nombre d'iterations arbitraire
        a_k = (np.dot(r.T,r))/(np.dot(p.T,np.dot(A,p)))
        x_next = x + (a_k * p)
        r_next = r - (a_k * np.dot(A,p))
        res.append(x_next)
        if np.linalg.norm(r_next) < precision: 
            return res
        b_k = (np.dot(r_next.T,r_next))/(np.dot(r.T,r))
        p = r_next + (b_k * p)
        #passage au k+1
        r = r_next
        x = x_next
    return res


N = 100
x = np.arange(0,N, 1)
ys = []
ys2 = []
from part1 import facto_cholesky
j = N
for i in range(1):
    print(j)
    # for k in range(1):
    ReelA = symetrik_generator(random.randint(0, j),j)
    Reelx = np.array([random.uniform(0, 10) for iter in range(j)])
    Reelx = Reelx.reshape((j,1))
    Reelb = np.dot(ReelA,Reelx)
    tmp_x = np.zeros((j,1)) 
    GuessedXs = conjugateGradient2Iter(ReelA,Reelb,tmp_x)
    GuessedXs2 = conjugateGradientPrecondIter(ReelA,facto_cholesky(ReelA),Reelb,tmp_x)
    tmp_diffs = []
    tmp_diffs2 = []
    for guess in GuessedXs:
        diff = np.linalg.norm(np.true_divide(np.abs(Reelx-guess),Reelx))
        tmp_diffs.append(diff)
    for guess2 in GuessedXs2:
        diff2 = np.linalg.norm(np.true_divide(np.abs(Reelx-guess2),Reelx))
        tmp_diffs2.append(diff2)
    ys.append(tmp_diffs)
    ys2.append(tmp_diffs2)

for e in ys:
    print("ys")
    x_tmp = np.arange(0,len(e), 1)
    plt.plot(x_tmp,e, label="Sans Precond", color="blue")

for e in ys2:
    print("ys2")
    x_tmp = np.arange(0,len(e), 1)
    plt.plot(x_tmp,e, label="Avec Precond", color="red")


plt.yscale("log")
plt.show()
