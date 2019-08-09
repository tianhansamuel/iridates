from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

nx=100
ny=100
EQ=[[0 for i in range(nx)] for j in range(ny)]
J2=0.001
LX=10
LY=10

ff=open('energy_fluctuation.dat','w')
for v1 in range(0,nx+1):
    phi=pi*v1*1.0/nx
    for v2 in range(0,ny+1):
       theta=pi*v2*1.0/ny
       E=0.0     
       for i1 in range(0,2*LX+1):
        qx=2*pi*(i1-LX)*1.0/(3*LX)      
        ymin=-4*pi*1.0/(3*sqrt(3))+1.0/(sqrt(3))*fabs(qx)
        ymax=4*pi*1.0/(3*sqrt(3))-1.0/(sqrt(3))*fabs(qx)
        for i2 in range(0,4*LY+1):
         qy=2*pi*i2*1.0/(3*LX)+ymin
         if qy<ymax:
          gamma1=1+4.0*(cos(sqrt(3)*0.5*qy))**2+4.0*cos(1.5*qx)*cos(sqrt(3)*0.5*qy)
          waz=cos(-sqrt(3)*qy)
          wazc=sin(-sqrt(3)*qy)
          wax=cos(1.5*qx+sqrt(3)*0.5*qy)
          waxc=sin(1.5*qx+sqrt(3)*0.5*qy)
          way=cos(-1.5*qx+sqrt(3)*0.5*qy)
          wayc=sin(-1.5*qx+sqrt(3)*0.5*qy)
          coz=(sin(theta))**2
          cox1=(cos(theta)*cos(phi))**2-(sin(phi))**2
          cox2=sin(2.0*phi)*cos(theta)
          coy1=(sin(phi)*cos(theta))**2-(cos(phi))**2
          gammaxy=(coz*waz+cox1*wax+coy1*way+cox2*waxc-cox2*wayc)**2+(coz*wazc+cox1*waxc+coy1*wayc-cox2*wax+cox2*way)**2
          tt1=-2.0*J2*((sin(theta)*cos(phi))**2)*wax         
          tt2=-2.0*J2*((sin(phi)*sin(theta))**2)*way
          tt3=-2.0*J2*(cos(theta))**2*waz
          gammaz=3+2.0*J2+tt1+tt2+tt3
          E=E+sqrt(gammaz**2-(sqrt(gamma1)+J2*sqrt(gammaxy))**2)
       E=(E-884.4)*100  
       ff.write(str(E)+'\t')
    ff.write('\n')
X=np.arange(0,1,1.0/nx)
Y=np.arange(0,1,1.0/ny)
X,Y=np.meshgrid(X, Y)
plt.ylabel('phi/pi')
plt.xlabel('theta/pi')
surf=ax.plot_surface(X, Y, EQ, rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=False)
plt.show()
