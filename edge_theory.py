from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from numpy import linalg as LA
from math import exp


ff=open('edge_theory.dat','w')

t=1
t1=0.5
NB=100
EE=[0]*(NB-1)
EE2=[0]*(NB-1)
KK=[0]*(NB-1)
LP1=[0]*(NB-1)
LP2=[0]*(NB-1)



for i in range(0,NB-1):
    kx=pi*(i+1)*1.0/(sqrt(3)*NB)+pi/(2*sqrt(3))  ## attention c'est le PZB pour le systeme 1D.
    KK[i]=kx
    CC1=t**4+16*t1**4*(sin(sqrt(3)*kx))**2+8*t1**2*t**2+8*t**4*(cos(sqrt(3)/2*kx))**2*cos(sqrt(3)*kx)
    CC2=(t**2+4*t1**2*(1-cos(sqrt(3)*kx)))
    E=4*t**2*(cos(sqrt(3)*0.5*kx))**2+4*t1**2*(sin(sqrt(3)*kx))**2+4*t**4*(cos(sqrt(3)*0.5*kx))**4/(4*t1**4*cos(sqrt(3)*kx)-t**2-4*t1**2)  
    #E=t**2+4*t1**2+4*t1**2*cos(sqrt(3)*kx)+4*t**2*(cos(sqrt(3)/2*kx))**2+4*t1**2*(sin(sqrt(3)*kx))**2-CC1/CC2
    if E>=0:
      EE[i]=sqrt(E)
      EE2[i]=-EE[i]
    else:
      EE[i]=3
      EE2[i]=-3
    L=(1+4*t1**2)/(2*t1**2*cos(sqrt(3)*kx))
    M=-L/2-sqrt((L/2)**2-1)
    N=cos(sqrt(3)/2*kx)/(t1**2*cos(sqrt(3)*kx)*(1+1.0/M))
    lamb1=abs(N/2-sqrt(abs(N**2-4*M)/2))
    lamb2=abs(N/2+sqrt(abs(N**2-4*M)/2))
    LP1[i]=-log(lamb1)
    LP2[i]=-log(lamb2)
      
    ff.write(str(KK[i])+'\t')
    ff.write(str(EE[i])+'\t')
    ff.write(str(EE2[i])+'\t')
    ff.write(str(LP1[i])+'\t')
    ff.write(str(LP2[i])+'\n')
LX=LY=50

NB1=(LX+1)*(LY+1)
E0=zeros([(LX+1),(LY+1)])
K0=[0]*(LX+1)
EE0=[0]*(LX+1)

#for i1 in range(0,LY+1):
#  ky=4*pi*i1/(3*LY)      
#  ymin=-4*pi/(3*sqrt(3))+1.0/(sqrt(3))*fabs(ky)
#  ymax=4*pi/(3*sqrt(3))-1.0/(sqrt(3))*fabs(ky)
#  E=30
#  for i2 in range(0,LX+1):
#   kx=(ymax-ymin)*i2*1.0/LX+ymin
#   E0[i1,i2]=sqrt(1+4*(cos(sqrt(3)*0.5*kx))**2+4*cos(sqrt(3)*0.5*kx)*cos(1.5*ky)+t1**2*((2*sin(sqrt(3)*kx))**2+(2*sin(sqrt(3)*0.5*kx+1.5*ky))**2+(2*sin(sqrt(3)*0.5*kx-1.5*ky))**2))
#for i2 in range(0,LY+1):
#  E=30
#  pod=0
#  for i1 in range(0,LX+1):
#     if E0[i1,i2]<E:
#        E=E0[i1,i2]
#        pod=i1
#  ky=4*pi*pod/(3*LY)      
#  ymin=1.0/(sqrt(3))*fabs(ky)
#  ymax=8*pi/(3*sqrt(3))-1.0/(sqrt(3))*fabs(ky)   
#  EE0[i2]=E
#  K0[i2]=(ymax-ymin)*i2*1.0/LX+ymin
#print K0
#print EE0
plot(KK,LP1,'r')
plot(KK,LP2,'b')
print kx
#plot(K0,EE0,'g')
#plot(KK,ET,'r')
#plot(KK,ET1,'b')
plt.grid()
plt.show()

