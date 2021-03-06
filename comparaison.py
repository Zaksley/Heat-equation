import numpy as np 
import numpy.linalg as linalg
import math as math 
import matplotlib.pyplot as plt

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







from part2 import conjugateGradient2  , conjugateGradientPrecond
def temperature_centresolve ( N  , T ) :
    A=tridiag (N) # construction de A 
    B= b_centre( N , T)
    t=np.eye  (N )
    x=t.reshape ( N*N , 1 ) 
    #tmp=conjugateGradient( A , B , x, imax=10**6, precision=1e-10)
    #tmp=tmp.reshape( N ,N )
   # for i in range ( 0 , N) :
    #        if ( tmp[i][i]!= 0 ) :
     #           print (" l'indice diagonale non nulle   " , i )
    # return conjugateGradient2( A , B , x, imax=10**6, precision=1e-10) 
    #print ("avec la fct prédifinie" ,  linalg.solve ( A , B ) ) 
    return linalg.solve ( A , B ) 



def b_centre(N , T ):
   # X=np.linspace ( 0 , 1 , N+1 ) 
   # Y=np.linspace ( 0 , 1 , N+1 )
    B=np.zeros ( (N*N , 1) ) ; 
    B[ (N//2)*N + (N//2)] = T*(-1) 
    return (B)  



def temperature_centregrad ( N  , T ) :
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
    
    

def comparaison ( N , T ) : 
    grad=temperature_centregrad ( N , T ) 
    solve=temperature_centresolve ( N , T ) 
    x=np.zeros ( (N*N ,  1) )
    for i in range (0 , N*N ) : 
        x[i]=solve[i]-grad[i] 
    for i in range ( 0 , N*N ) :
        print ( " la " , i , "eme difference est " , x[i] )  
    

comparaison ( 4 , 25 )