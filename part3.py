import numpy as np 
import numpy.linalg as linalg
import math as math 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import PIL
PIL.__version__ 
from PIL import Image
" Pour la premiére question : " 
	
   "t ( i , j ) = 4 *b(i , j ) * h^2 - t(i+1 , j) - t(i-1 , j)-t(i , j+1)-t(i , j-1)) "
   	
   "les indices des coefficients varient entre [ i-1 , i+1 ] et [ j-1 , j+1] "
   "ce qui induit nécessairement lutilisation des blocs tridiagonales . En effet  , "
    "soit j fixe dans [0 ,N ]"
   	"t(j,1) = a*t(j+1 , 1) + b*t(j+1,1) + c*t(j+1,0)"
   	
    "t(j,2) = a*t(j+1 , 2) + b*t(j+1,2) + c*t(j+1,1)"
 
    "t(j,3) = a*t(j+1 , 3) + b*t(j+1,3) + c*t(j+1,2) "
   
   	
   "t(j,n-1) = a*t(j+1 , n-1) + b*t(j+1,n-1) + c*t(j+1,n-2) "


"La construction de la matrice tridiagonale par bloc "
def tridiag (N ):
    N2=N*N 
    A=np.diag ( [-4]*N2 )  #On commence par définir la diagonale principale 

    B=np.eye ( N2 , N2 , 1) #La diagonale supérieure des blocs supérieurs 
    for i in range (1 , N): #Annulation des termes pour conserver une matrice tridiagonale par blocs 
        B[i*(N-1)][i*N]=0  
    C=B.T  #la diagonale inférieure des blocs supérieurs  
   
    D=np.eye (N2 , N2 , N) #La diagonale  supérieurs des blocs inférieurs  
    E=D.T #la diagonale inférieure des blocs inférieurs 
    A=A+B+C+D+E 
    
    A=A/((N+1)**2)
  

    return (A)





" Pour un radiateur au centre " 
#construction de la matrice colonne du second membre 

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
    return conjugateGradient2( A ,  B , x, imax=10**6, precision=1e-10) #appel a conjugateGradient


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

#construction de la matrice colonne du second membre 
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

