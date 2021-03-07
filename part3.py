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
    
    A=A/((N+1)**2)
  

    return (A)

tridiag ( 4 )



" Pour un radiateur au centre "

def b_centre(N , T ):
    B=np.zeros( (N*N,1) ) 
    B[ (N//2)*N + (N//2)] = T*(-1) 
    return (B) 



" Resolution dans le cas d'un radiateur au centre "


from part2 import conjugateGradient2 
def temperature_centre( N  , T ) :
    A=tridiag (N) # construction de A 
    B= b_centre( N , T) # construction de b 
    t=np.eye  (N ) 
    x=t.reshape ( N*N , 1 ) 
    return conjugateGradient2( A ,  B , x, imax=10**6, precision=1e-10) 


" Image correspondante "
def imagecentre ( N , T ) : 
    temp = temperature_centre( N , T ) 
    temp=temp.reshape (N , N )
    fig, ax = plt.subplots()
    X=np.linspace ( 0 , 1 , N+1 ) 
    Y=np.linspace ( 0 , 1 , N+1 )
   
    ax.set_title("temperature obtenue pour un radiateur au centre ")
    fig.tight_layout()

    plt.pcolor(temp, cmap=plt.cm.hot, vmin=0, vmax=50000)
    plt.show() 

    

##########################################################################################################

" Pour un radiateur mis au nord "


def b_cote(N , T ):
    B=np.zeros ( (N*N , 1 ) ) 
    for i in range ( N) :
        B[ (N-1)*N + i ] = T*(-1) 
    return (B) 


" temperature correspondante "

def temperature_cote ( N  , T ) :
    A=tridiag (N) # construction de A 
    B= b_cote( N , T) #construction de B 
    t=np.eye  (N )
    x=t.reshape ( N*N , 1 ) 
    
    return conjugateGradient2( A , B , x, imax=10**6, precision=1e-5) 


" Image correspondante "
def imagecote ( N , T ) : 
    temp = temperature_cote ( N , T ) 
    temp=temp.reshape (N , N )
    fig, ax = plt.subplots()
    X=np.linspace ( 0 , 1 , N+1 ) 
    Y=np.linspace ( 0 , 1 , N+1 )
   
    ax.set_title("temperature obtenue pour un radiateur au centre ")
    fig.tight_layout()

    plt.pcolor(temp, cmap=plt.cm.hot, vmin=0, vmax=50000)
    plt.show() 

