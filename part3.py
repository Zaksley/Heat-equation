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
#Le vecteur b en vrai c'est la fonction f ( x , y ) 
#si on suppose que la fonction f(x,y) est representee par la matrice B alors : 
# B=fromfunction(f,dim) 
#
#b = reshape de B 
#def f(x,y, N) :

 #       if x in [ N//2 , 3*N//4 ] and y in [ N//2 , 3*N//4 ] :
  #          return -25  
   #     else : 
    #       return 0 
def b_centre(N , T ):
   # X=np.linspace ( 0 , 1 , N+1 ) 
   # Y=np.linspace ( 0 , 1 , N+1 )
    B=np.zeros( (N*N,1) ) 
    B[ (N//2)*N + (N//2)] = T*(-1) 
    return (B) 




# resoudre le systeme dans le cas d'un radiateur mis au centre du carre sachant que f(x , y)  l'apport de chaleur exterieur en ce meme point : 
#alors placer un radiateur au centre du carre
# alors "normalement" , f( 1/2 , 1/2 ) = 1 et f( x , y) = 0 Vx, y != 1/2 

# donc b de f est 


from part2 import conjugateGradient2 
def temperature_centre( N  , T ) :
    A=tridiag (N) # construction de A 
    B= b_centre( N , T)
    t=np.eye  (N )
    x=t.reshape ( N*N , 1 ) 
   # TE = np.eye ( N )
    #tmp=conjugateGradient( A , B , x, imax=10**6, precision=1e-10)
    #tmp=tmp.reshape( N ,N )
   # for i in range ( 0 , N) :
    #        if ( tmp[i][i]!= 0 ) :
     #           print (" l'indice diagonale non nulle   " , i )
    return conjugateGradient2( A ,  B , x, imax=10**6, precision=1e-10) 



def graphe3d(N , f  ) :
    VX = np.linspace(0, 1.0, N)
    VY = np.linspace(0, 1.0, N)
    X,Y = np.meshgrid(VX, VY)
    x=[i for i in range (0 , N )] 
    y=[i for i in range ( 0 , N ) ] 
    
    def t( i , j  ) : 
        temp = temperature_centre( N , f ) 
        temp=temp.reshape(N , N )
        return temp[i][j]
    Z = t(x , y)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.view_init(40, -30)
    ax.plot_surface(X, Y, Z)
    plt.show()
    input('press <ENTER> to continue') 
    plt.close




def imagepix_centre ( T, N ) :
    size = (N,N)
    im = Image.new('RGB',(N , N ))
    print ("1  ça passe  ?")
    pix = im.load() 
    print (" 2 ça passe ?")
    temp = temperature_centre( N , T ) 
    print (" 3 ça passe ?")
    temp=temp.reshape (N , N )
    temp=temp.T
    maxi=np.amax ( temp )
    print (" ça passe ?")
    for i in range(0 , size[0]):
        for j in range( 0 , size[1]):
            #if ( temp[i][j] < 1 ) :
            #    pix[i , j] = ( int (temp[i][j] *255) , int (temp[i][j]*69) , 0 ) 
           # else : 
                if ( temp[i][j] <0 ) :
                    pix[i , j] = ( 0.04 , 0.05 , 0.06 )
                else  :  
                    # print (" le coefficient du rouge ," , temp[i][j] / (len(str(maxi))) *255/(T*T) ) 
                    pix[i,j] = ( int ( temp[i][j] / (len(str(maxi)) ) *255/(T*T)*10 ) , int (temp[i][j] / (len(str(maxi)) )* 6.9*10/(T*T) )  ,  0)  
                
           # print ( "pixel " , i , j , "est " , pix[i , j ] , "\n" )
    im.save('radiateuraucentre.png')


def image ( T, N) : 
    rouge = [0.807, 0.066, 0.149]
    jaune = [0.802, 0.101, 0.066] 
    # Image de départ, entièrement noire
    image = np.zeros((N, N,3)) 
    temp = temperature_centre ( N , T ) 
    temp=temp.reshape (N , N )
    # Remplissage  
    maxi=np.amax ( temp )
    for i in range (0, N ) : 
        for j in range (0, N) : 
            if ( temp[i][j] < 1 and temp[i][j] > 0 ) :
                image[i][j] = ( temp[i][j]  , temp[i][j]*69/255 , 0 ) 
            else : 
               image[i][j] = ( temp[i][j] / maxi   , temp[i][j] / maxi* 69/255  ,  0)  
            if ( temp[i][j] *255 < 0) :
                image[i][j]= ( temp[i][j] *(-1) , temp[i][j] *69/255 * (-1) , 0 )
             
    plt.imsave ('image.png' , image)
    
    plt.show()




"""
    BGR : Blue Green Red
    temp :       [0, 50, 100, 150, 200]
    temp color : [B, BG, G,   GR,  R]
"""
def tempToColor(t, B=0, G=1500, R=3000):
    if t < B :
        return (0,0,255)
    if t > R : 
        return (255,0,0)
    r = 0
    g = 0
    b = 0
    if t in [B,G]:
        g = 0
    # r = (t - G)*(t -(R+(R-G)))*-255
    # g = (t - B)*(t - R)*-255
    # b = (t - (B-(G-B)))*(t - G)*-255
    # print(str(r)+" "+str(g)+" "+str(b))
    return (int(r),int(g),int(b))
        
    


def imagepix_cote ( T , N ) :
    
    size = (N,N)
    im = Image.new('RGB',(N , N ))
    pix = im.load() 
    temp = temperature_cote ( N , T ) 
    temp=temp.reshape (N , N )  
    maxi=np.amax ( temp )
    
    for i in range(0 , size[0]):
        for j in range( 0 , size[1]):
            print(str(i)+str(j))
            pix[i,j] = tempToColor(temp[i][j])
            # if ( temp[i][j] < 1 and temp[i][j] > 0 ) :
            #     pix[i , j] = ( int (temp[i][j] *255) , int (temp[i][j]*69) , 0 ) 
            # else : 
            #     pix[i,j] = ( int ( temp[i][j] / len(str(maxi))  *255 ) , int (temp[i][j] / len(str(maxi))* 69 ) ,  0)  
            # if ( temp[i][j] *255 < 0 ) :
            #      pix[i , j] = ( int (temp[i][j] *255 *10 )*(-1) , int (temp[i][j]*69 *10) * (-1) , 0 )
           # print ( "pixel " , i , j , "est " , pix[i , j ] , "\n" )
    im.save('radiateurcote.png')



#cote
def b_cote(N , T ):
   # X=np.linspace ( 0 , 1 , N+1 ) 
   # Y=np.linspace ( 0 , 1 , N+1 )
    B=np.zeros ( (N*N , 1 ) ) 
    for i in range ( N) :
        B[ (N-1)*N + i ] = T*(-1) 
    return (B) 




def temperature_cote ( N  , T ) :
    A=tridiag (N) # construction de A 
    B= b_cote( N , T)
    t=np.eye  (N )
    x=t.reshape ( N*N , 1 ) 
    #tmp=conjugateGradient( A , B , x, imax=10**6, precision=1e-10)
    #tmp=tmp.reshape( N ,N )
   # for i in range ( 0 , N) :
    #        if ( tmp[i][i]!= 0 ) :
     #           print (" l'indice diagonale non nulle   " , i )
    return conjugateGradient2( A , B , x, imax=10**6, precision=1e-5) 
   # print ("avec la fct prédifinie" ,  linalg.solve ( A , B ) ) 
    #return linalg.solve ( A , B )



#imagepix_cote ( 20, 8) 
imagepix_centre ( 1000, 20) 
A=  temperature_centre( 4 , 25 ) 
#print (A)



def imagehector ( N , T ) : 
    temp = temperature_centre ( N , T ) 
    temp=temp.reshape (N , N )
    
    fig, ax = plt.subplots()
    # im = ax.imshow(temp, cmap=plt.cm.hot) 
    # We want to show all ticks...
    X=np.linspace ( 0 , 1 , N+1 ) 
    Y=np.linspace ( 0 , 1 , N+1 )
    # ax.set_xticks(np.arange(X))
    # ax.set_yticks(np.arange(Y)) 
    
    ax.set_title("temperature ")
    fig.tight_layout()

    plt.pcolor(temp, cmap=plt.cm.hot, vmin=0, vmax=50000)
    # plt.colorbar()
    plt.show() 


imagehector ( 100 , 25 )
