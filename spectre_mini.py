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

J2=0.01
LX=20
LY=20
EQ=[[0 for i in range(2*LX)] for j in range(2*LY)]
theta=0
phi=pi*1.0/6
ff=open('energy_fluctuation.txt','w')
for i1 in range(0,2*LX):
      qx=2*pi*(i1-LX)/(3*LX)         
      ymin=max(-4*pi/(3*sqrt(3))-1.0/(sqrt(3))*qx,-4*pi/(3*sqrt(3))+1.0/(sqrt(3))*qx)
      ymax=min(4*pi/(3*sqrt(3))-1.0/(sqrt(3))*qx,4*pi/(3*sqrt(3))+1.0/(sqrt(3))*qx)
      for i2 in range(0,2*LY):
          qy=(ymax-ymin)*i2/(2*LY)+ymin
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
          EQ[i1][i2]=sqrt(gammaz**2-(sqrt(gamma1)+J2*sqrt(gammaxy))**2)
          ff.write(str(EQ[i1][i2])+'\t')
X=np.arange(-2*pi/3,2*pi/3,4*pi/(6*LX))
Y=np.arange(-2*pi/3,2*pi/3,4*pi/(6*LY))
X,Y=np.meshgrid(X, Y)
plt.xlabel('qy')
plt.ylabel('qx')
surf=ax.plot_surface(X, Y, EQ, rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=False)
plt.show()
