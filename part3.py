import numpy as np 
import numpy.linalg as linalg
import math as math 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def tridiag (N ):
    N2=N*N 
    A=np.diag ( [-4]*N2 ) 

    B=np.eye ( N2 , N2 , 1)
    for i in range (1 , N):
        B[i*(N-1)][i*N]=0  
    C=B.T  
    D=np.eye (N2 , N2 , N)
    E=D.T
    A=A+B+C+D+E 
    A=A/(N+1)*(N+1)

    return (A)


#Le vecteur b en vrai c'est la fonction f ( x , y ) 
#si on suppose que la fonction f(x,y) est representee par la matrice B alors : 
# B=fromfunction(f,dim) 
#
#b = reshape de B 
def f(x,y) :

        if x in [1/2 , 3/4] and y in [1/2 , 3/4] :
            return -25 
        else :
           return 0 
def b(N , f ):
    X=np.linspace ( 0 , 1 , N+1 ) 
    Y=np.linspace ( 0 , 1 , N+1 )
    B=np.zeros ( (N,N) ) ; 
    for i in range ( 0 , N) :
        for j in range ( 0 , N ) :
            B[i][j] = f( X[i] , Y[j]) 

    b=B.reshape ( N*N,1)
    return (b)

b(4 , f)
# resoudre le systeme dans le cas d'un radiateur mis au centre du carre sachant que f(x , y)  l'apport de chaleur exterieur en ce meme point : 
#alors placer un radiateur au centre du carre
# alors "normalement" , f( 1/2 , 1/2 ) = 1 et f( x , y) = 0 Vx, y != 1/2 

# donc b de f est 


from part2 import conjugateGradient 
def temperature ( N  , f ) :
    A=tridiag (N) # construction de A 
    B= b( N , f)
    print ("la matrice de la fonction f est " , B )
    t=np.zeros ( (N,N)) 
    x=t.reshape ( N*N , 1 ) 
    return conjugateGradient( A , B , x, imax=10**6, precision=1e-10) 
     

   
   
def graphe(N , f  ) :
    VX = np.linspace(0, 1,  N)
    VY = np.linspace(0, 1, N)
    X,Y = np.meshgrid(VX, VY)
    temp = temperature ( N , f ) 
    Z=[] 
    for i in range ( 0 , N) : 
        for j in range (0 , N ) : 
            Z=np.append (Z , temp[i*(N-1) + j ])
    fig = plt.figure()
    plt.plot (X, Y, Z)
    plt.show()
    plt.close

graphe ( 4 , f )

