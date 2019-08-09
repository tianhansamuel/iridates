from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from numpy import linalg as LA
from math import exp

LL=70  ## longeur du systeme
NX=100 ## nombre du pas

MM=zeros([4*(LL+1),4*(LL+1)],complex) ## initialisation d'une matrice du taille 4*(LL+1) puisque l'Hamitonian 1D avec kx conserve est un bloc 4*4


RZ=sqrt(3)
dd1=-sqrt(3)*0.5
dd2=sqrt(3)*0.5
SS=4*(LL+1)*NX*2  ## nombre de point trace

EE=[3]*SS   ## un array de SS dimension 
KK=[0]*SS   ## Kx
KX=[0]*(NX+1)
SXB=[0]*(NX+1)
SYB=[0]*(NX+1)
SZB=[0]*(NX+1)
SXB1=[0]*(NX+1)
SYB1=[0]*(NX+1)
SZB1=[0]*(NX+1)
Ct=0        ## compteur
ff=open('shitade_cylinder.dat','w')  ## ouvrir un autre fichier pour tracer avec gnuplot
ff1=open('shitade_cylinder_spin.dat','w')
#ff2=open('shitade_cylinder_spin_t.dat','w')
ff3=open('shitade_cylinder_spin_layer.dat','w')

SX=ones([NX+1,LL+1])
SY=ones([NX+1,LL+1])
SZ=ones([NX+1,LL+1])
SX1=ones([NX+1,LL+1])
SY1=ones([NX+1,LL+1])
SZ1=ones([NX+1,LL+1])
LAYER=[0]*(LL+1)
SXL=[0]*(LL+1)
SYL=[0]*(LL+1)
SZL=[0]*(LL+1)

NBB=30
SXT=[0]*NBB
SYT=[0]*NBB
SZT=[0]*NBB


t1=0.5           ## couplage spin-orbit  t=1
ttt=[0]*NBB
#for ib in range(1,NBB):
# t1=1.5/NBB*ib
# ttt[ib]=t1
for j in range(0,NX+1):
    k=2*pi*j*1.0/(sqrt(3)*NX)  ## attention c'est le PZB pour le systeme 1D.
    KX[j]=k
    for i in range(0,LL):

####################################################################
## Les elements de matrice pour des couplages dans le bulk  ########

## 1er voisins
    	MM[4*i+1][4*i]=MM[4*i][4*i+1]=2*cos(k*dd2)   #cib^dagger cia
    	MM[4*i+3][4*i+2]=MM[4*i+2][4*i+3]=2*cos(k*dd2)   #cia^dagger cib
    	MM[4*(i+1)+1][4*i]=MM[4*i][4*(i+1)+1]=MM[4*(i+1)+3][4*i+2]=MM[4*i+2][4*(i+1)+3]=1
#####################################################################
## 2nd voisins
    	MM[4*(i+1)][4*i]=MM[4*(i+1)+3][4*i+3]=1J*t1*(cos(k*dd1)-1J*sin(k*dd1)) #ci+1aup^dagger ciaup ci+1bdown^dagger cibdown    z
    	MM[4*(i+1)+2][4*i+2]=MM[4*(i+1)+1][4*i+1]=-1J*t1*(cos(k*dd1)-1J*sin(k*dd1)) #ci+1adown^dagger ciadown ci+1bup^dagger cibup     z
    	MM[4*i+1][4*(i+1)+1]=MM[4*i+2][4*(i+1)+2]=1J*t1*(cos(k*dd2)-1J*sin(k*dd2)) #cib^daggerup ci+1bup ciadown^dagger ci+1adown     z
    	MM[4*i][4*(i+1)]=MM[4*i+3][4*(i+1)+3]=-1J*t1*(cos(k*dd2)-1J*sin(k*dd2)) #ciaup^dagger ci+1aup  cibdown^dagger cibdown        z
    	
        MM[4*i][4*i+2]=MM[4*i+2][4*i]=2*t1*sin(k*RZ)  #x
    	MM[4*i+1][4*i+3]=MM[4*i+3][4*i+1]=-2*t1*sin(k*RZ)#x
### A    	
        MM[4*i][4*(i+1)+2]=-1J*1J*t1*(cos(k*dd1)-1J*sin(k*dd1))#y
        MM[4*i+2][4*(i+1)]=1J*1J*t1*(cos(k*dd1)-1J*sin(k*dd1))#y

	MM[4*(i+1)][4*i+2]=1J*1J*t1*(cos(k*dd2)-1J*sin(k*dd2))#y
        MM[4*(i+1)+2][4*i]=-1J*1J*t1*(cos(k*dd2)-1J*sin(k*dd2))#y
### B
        MM[4*i+1][4*(i+1)+3]=1J*1J*t1*(cos(k*dd1)-1J*sin(k*dd1))#y
        MM[4*i+3][4*(i+1)+1]=-1J*1J*t1*(cos(k*dd1)-1J*sin(k*dd1))#y
        
        MM[4*(i+1)+1][4*i+3]=-1J*1J*t1*(cos(k*dd2)-1J*sin(k*dd2))#y
        MM[4*(i+1)+3][4*i+1]=1J*1J*t1*(cos(k*dd2)-1J*sin(k*dd2))#y

#####################################################################
##  Les couplages aux bords
##  1er voisin

    MM[4*LL+1][4*LL]=MM[4*LL+3][4*LL+2]=2*cos(k*dd2)   #cib^dagger cia
    MM[4*LL][4*LL+1]=MM[4*LL+2][4*LL+3]=2*cos(k*dd2)   #cia^dagger cib
    	
#######################################################################
##  2nd voisin        

    MM[4*LL][4*LL+2]=MM[4*LL+2][4*LL]=2*t1*sin(k*RZ)  #x
    MM[4*LL+1][4*LL+3]=MM[4*LL+3][4*LL+1]=-2*t1*sin(k*RZ)#x
    	
###############################################################

    VAL,VEC=LA.eig(MM)   ## diagonalisation Val valeur propre, Vec vecteur propre
    VAL=VAL.real
    VAL_sorted = np.sort(VAL)
    VEC_sorted = VEC[:,VAL.argsort()]
    
    for pos in range (0,LL+1): ## pos eme couche 
        LAYER[pos]=pos 
        SZ[j,pos]=(VEC_sorted[4*pos,2*(LL+1)].conjugate()*VEC_sorted[4*pos,2*(LL+1)]-VEC_sorted[4*pos+3,2*(LL+1)].conjugate()*VEC_sorted[4*pos+3,2*(LL+1)]).real+(VEC_sorted[4*pos+1,2*(LL+1)+1].conjugate()*VEC_sorted[4*pos+1,2*(LL+1)+1]-VEC_sorted[4*pos+3,2*(LL+1)+1].conjugate()*VEC_sorted[4*pos+3,2*(LL+1)+1]).real
        SX[j,pos]=(VEC_sorted[4*pos+1,2*(LL+1)].conjugate()*VEC_sorted[4*pos+3,2*(LL+1)]+VEC_sorted[4*pos+3,2*(LL+1)].conjugate()*VEC_sorted[4*pos+1,2*(LL+1)]).real+(VEC_sorted[4*pos+1,2*(LL+1)+1].conjugate()*VEC_sorted[4*pos+3,2*(LL+1)+1]+VEC_sorted[4*pos+3,2*(LL+1)+1].conjugate()*VEC_sorted[4*pos+1,2*(LL+1)+1]).real
        SY[j,pos]=(-VEC_sorted[4*pos+1,2*(LL+1)].conjugate()*VEC_sorted[4*pos+3,2*(LL+1)]+VEC_sorted[4*pos+3,2*(LL+1)].conjugate()*VEC_sorted[4*pos+1,2*(LL+1)]).imag+(-VEC_sorted[4*pos+1,2*(LL+1)+1].conjugate()*VEC_sorted[4*pos+3,2*(LL+1)+1]+VEC_sorted[4*pos+3,2*(LL+1)+1].conjugate()*VEC_sorted[4*pos+1,2*(LL+1)+1]).imag

    for ii in range(0,4*(LL+1)):
      KK[Ct]=k
      ff.write(str(KK[Ct])+'\t')
      EE[Ct]=VAL[ii]
      ff.write(str(EE[Ct].real)+'\n')
      Ct=Ct+1
    SXB[j]=SX[j,0]
    SYB[j]=SY[j,0]
    SZB[j]=-SY[j,0]    
    ff1.write(str(KX[j])+'\t')
    ff1.write(str(SXB[j])+'\t')
    ff1.write(str(SYB[j])+'\t')
    ff1.write(str(SYB[j])+'\n')
print VAL_sorted[2*(LL+1)]
print VAL_sorted[2*(LL+1)+1]
SXM=SYM=SZM=0
for j in range(0,NX+1):
    if fabs(SXB[j])>SXM:
       SXM=fabs(SXB[j])
       SYM=SYB[j]
       SZM=SZB[j]
 #SXT[ib]=SXM       
 #SYT[ib]=SYM       
 #SZT[ib]=SZM       
for i in range(0,LL+1):
    SXM=SYM=SZM=0
    for j in range(0,NX):
       if fabs(SX[j,i])>SXM:
          SXM=fabs(SX[j,i])
       if fabs(SY[j,i])>SYM:
          SYM=fabs(SX[j,i])
       if fabs(SZ[j,i])>SZM:
          SZM=fabs(SX[j,i])
    ff3.write(str(i)+'\t')
    ff3.write(str(SXM)+'\t')
    ff3.write(str(SYM)+'\t')
    ff3.write(str(SZM)+'\n')
 #ff2.write(str(ttt[ib])+'\t')
 #ff2.write(str(SXT[ib])+'\t')
 #ff2.write(str(SYT[ib])+'\t')
 #ff2.write(str(SZT[ib])+'\n')
#print SXM
#print SYM
#print SZM
#print SY[:,0]
#print SY[:,LL]

#print 'SX'
#print SX[:,0]
#print 'SY'
#print SY[:,0]
#print 'SZ'
#print SZ[:,0]  
#plot(KK,EE, 'r.')
#plot(KX,SXB1,'r')
#plot(KX,SYB1,'g')
#plot(KX,SZB1,'b')
#plot(LAYER,SX[LMAX1,:],'r')
#plot(LAYER,SY[LMAX2,:],'g')
#plot(LAYER,SZ[LMAX3,:],'b')
#plot(ttt,SXT,'r')
#plot(ttt,SYT,'g')
#plot(ttt,SZT,'b')
#plt.xlabel('t1/t')
#plt.xlabel('k_x')
#plt.ylabel('energy')
#plt.ylabel('energy and spin_polarization_upper_band_lower_edge')
#plt.title('x_red,y_green,z_blue')
#plt.xlabel('layer')
#plt.ylabel('spin_polarization_upper_band')
#plt.grid()
plt.show()


