from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt


LX=10
LY=10
N=20
S=zeros([N])
J=zeros([N])
SH=0.48
for i in range(1,N):
  J[i]=SH*i/N
  J2=J[i]
  Sz=0
  for i1 in range(1,2*LX):
      qx=2*pi*(i1-LX)/30         
      ymin=max(-4*pi/(3*sqrt(3))-1/(sqrt(3))*qx,-4*pi/(3*sqrt(3))+1/(sqrt(3))*qx)
      ymax=min(4*pi/(3*sqrt(3))-1/(sqrt(3))*qx,4*pi/(3*sqrt(3))+1/(sqrt(3))*qx)
      for i2 in range(1,2*LY):
            qy=(ymax-ymin)*i2/60+ymin
            gamma1=cos(qx)+2*cos(0.5*qx)*cos(sqrt(3)/2*qy)
            gamma2=cos(sqrt(3)*qy)+cos(3/2*qx+sqrt(3)/2*qy)+cos(3/2*qx-sqrt(3)/2*qy)
            gammaz=3-2*J2+J2/2*(cos(sqrt(3)*qy)+cos(3/2*qx-sqrt(3)/2*qy))
            gammaxy=gamma1-J2*gamma2-J2/2*(cos(sqrt(3)*qy)-cos(3/2*qx-sqrt(3)/2*qy))
            gammaxx=cos(3/2*qx-sqrt(3)/2*qy)
            gammayy=cos(sqrt(3)*qy)
            Sz=Sz+gammaz/(sqrt((3-2*J2+gamma1-J2*gamma2+J2*gammaxx)*(3-2*J2-gamma1+J2*gamma2+J2*gammayy)))-1
  S[i]=Sz/(LX*LY)
xlabel('J2/J1')
ylabel('Sz')
plot(J,S)
plt.show()
