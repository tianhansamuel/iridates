from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt

NN=30
LX=20
LY=20
theta=pi/2
phi=pi/2

for i in range(0,NN):
   J2=0.5*i/NN+0.1
   for i1 in range(0,2*LX):
      qx=2*pi*(i1-LX)/(3*LX)         
      ymin=max(-4*pi/(3*sqrt(3))-1/(sqrt(3))*qx,-4*pi/(3*sqrt(3))+1/(sqrt(3))*qx)
      ymax=min(4*pi/(3*sqrt(3))-1/(sqrt(3))*qx,4*pi/(3*sqrt(3))+1/(sqrt(3))*qx)
      for i2 in range(0,2*LY):
         qy=(ymax-ymin)*i2/(2*LY)+ymin
         gamma1=cos(qx)+2*cos(0.5*qx)*cos(sqrt(3)/2*qy)
         waz=cos(sqrt(3)*qy)
         wazc=sin(sqrt(3)*qy)
         wax=cos(3/2*qx+sqrt(3)/2*qy)
         waxc=sin(3/2*qx+sqrt(3)/2*qy)
         way=cos(3/2*qx-sqrt(3)/2*qy)
         wayc=sin(3/2*qx-sqrt(3)/2*qy)
         cox=sin(theta)**2
         coy1=(cos(theta)*cos(phi))**2-(sin(phi))**2
         coy2=sin(2*phi)*cos(theta)
         coz1=(sin(phi)*cos(theta))**2-(cos(phi))**2
         gammaxy=(waz*cox+coy1*way+coy2*waxc+coz1*way-coy2*wayc)**2+(cox*waz-coy2*wax+coy1*waxc+coy2*way+coz1*wayc)**2         
         gammaz=3+2*J2-2*J2*cos(2*theta)*cos(sqrt(3)*qy)+2*J2*(cos(phi))**2*(cos(2*theta)+(sin(phi))**2)*cos(3/2*qx+sqrt(3)/2*qy)+2*J2*(sin(phi)**2*cos(2*theta)+(cos(phi))**2)*cos(3/2*qx-sqrt(3)/2*qy)
         if (gammaz**2-(sqrt(gamma1)+J2*sqrt(gammaxy))**2)<0:
               print J2



