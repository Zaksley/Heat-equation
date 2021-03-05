import numpy as np 
import numpy.linalg as linalg
import math as math 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import PIL
PIL.__version__ 
from PIL import Image

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
    print (A) 
    A=A/((N+1)**2)
    print (A)

    return (A)

tridiag ( 4 )
#Le vecteur b en vrai c'est la fonction f ( x , y ) 
#si on suppose que la fonction f(x,y) est representee par la matrice B alors : 
# B=fromfunction(f,dim) 
#
#b = reshape de B 
def f(x,y) :

        if float(x) in [0.5 , float(4/7)] or  float(y) in [0.5 , float(4/7)] :
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

    
    print ("la matrice de lafocntion" , B ) 
    b=B.reshape ( N*N,1)
    return (b) 
b(16 , f )



# resoudre le systeme dans le cas d'un radiateur mis au centre du carre sachant que f(x , y)  l'apport de chaleur exterieur en ce meme point : 
#alors placer un radiateur au centre du carre
# alors "normalement" , f( 1/2 , 1/2 ) = 1 et f( x , y) = 0 Vx, y != 1/2 

# donc b de f est 


from part2 import conjugateGradient 
def temperature ( N  , f ) :
    A=tridiag (N) # construction de A 
    B= b( N , f)
    t=np.zeros ( (N,N)) 
    x=t.reshape ( N*N , 1 ) 
    t=conjugateGradient( A , B , x, imax=10**6, precision=1e-10) 
    for i in range (0 , N*N) :  
            if t[i] != 0 : 
                print ( "temperature non nulle en " , i  )
    return conjugateGradient( A , B , x, imax=10**6, precision=1e-10) 
     

   
   
def graphe3d(N , f  ) :
    VX = np.linspace(0, 1.0, N)
    VY = np.linspace(0, 1.0, N)
    X,Y = np.meshgrid(VX, VY)
    x=[i for i in range (0 , N )] 
    y=[i for i in range ( 0 , N ) ] 
    
    def t( i , j  ) : 
        temp = temperature ( N , f ) 
        temp=temp.reshape (N , N )
        return ( temp [i][j] )
    Z = t(x , y)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.view_init(40, -30)
    ax.plot_surface(X, Y, Z)
    plt.show()
    input('press <ENTER> to continue') 
    plt.close
graphe3d (  8 , f  )
def imagepix ( f , N ) :
    
    size = (N,N)
    im = Image.new('RGB',size)
    pix = im.load() 
    temp = temperature ( N , f ) 
    temp=temp.reshape (N , N )
    maxi=np.amax ( temp )
    print ("le max" , np.amax ( temp ))
    for i in range(size[0]):
        for j in range(size[1]):
                                    
            pix[i,j] = ( int ( temp[i][j] / maxi  *255 ) , int (temp[i][j] / maxi* 69 ) ,  0)  
          
    im.save('ima4.png')

imagepix ( f , 20)
def image ( f , N) : 
    rouge = [0.807, 0.066, 0.149]
    jaune = [0.802, 0.101, 0.066] 
    # Image de départ, entièrement noire
    image = np.zeros((N, N,3)) 
    temp = temperature ( N , f ) 
    temp=temp.reshape (N , N )
    # Remplissage  
    maxi=np.amax ( temp )
    for i in range (0, N ) : 
        for j in range (0, N) : 
            
            image[i][j] = [ temp[i][j] / maxi , temp[i][j] / maxi* 69/255  ,  0 ] 
    plt.imsave ('image.png' , image)
    
    plt.show()
image ( f , 100)