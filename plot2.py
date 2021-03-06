import matplotlib.pyplot as plt
import numpy as np 
import random
import math as m

from part2 import conjugateGradient2
from part1 import symetrik_generator



def conjugateGradientIter(A, b, x, imax=10**6, precision=1e-10):

    # Initialisation de r0 et p0
    r = b - np.dot(A, x)
    p = r
    rsquare_old = np.dot(r, r.T)
    cpt = 0
    for i in range(imax): #Nombre d'iterations arbitraire
        cpt += 1
        # Calcul de alpha_k+1
        Ap = np.dot(A, p)
        alpha = rsquare_old / np.dot(p.T, Ap)
        # Mise a jour de x et calcul de r_k+1
        x += np.dot(alpha, p)
        r -= np.dot(alpha, Ap)
        
        rsquare_new = np.dot(r.T, r)

        # On a atteint la precision attendue
        if(m.sqrt(rsquare_new) < precision): 
            return cpt

        # Sinon, calcul de p_k+1 et passage a l'iteration suivante
        p = r + p * (rsquare_new/rsquare_old)
        rsquare_old = rsquare_new

    # La precision n'a pas ete atteinte, on renvoie quand même le resultat
    return cpt

# y = []

# for i in range(1000):
#     tmp = np.array([random.randint(0, 9) for iter in range(i)])
#     A = symetrik_generator(random.randint(0, i),i)
#     tmp2 = np.zeros(i)
#     res = conjugateGradientIter(A,tmp,tmp2)
#     y.append(res)
#     print(i)


# x = np.arange(0,1000, 1)
# plt.plot(x,y)
# plt.show()

x = np.arange(0,50, 1)
y = []

for j in range(50):
    print(j)
    tmp_y = 0
    for k in range(50):
        ReelA = symetrik_generator(random.randint(0, j),j)
        Reelx = np.array([random.uniform(0, 10) for iter in range(j)])
        Reelx = Reelx.reshape((j,1))
        Reelb = np.dot(ReelA,Reelx)
        tmp_x = np.zeros((j,1)) 
        GuessedX = conjugateGradient2(ReelA,Reelb,tmp_x)
        diff = np.linalg.norm(GuessedX-Reelx)
        tmp_y += diff
    y.append(tmp_y)

plt.plot(x,y)
plt.yscale("log")
plt.show()
