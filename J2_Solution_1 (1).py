from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')

import mpl_toolkits.mplot3d.axes3d as p3
nx=1000
ny=1000
FF=[[0 for i in range(nx)] for j in range(ny)]

SS=nx*ny
t=[0]*SS
s=[0]*SS
l=0
ff=open('theta_phi.dat','w')
fftt=open('sphere_axis.dat','w')
for i in range(0,nx):
  theta=2*pi*i*1.0/nx
  for j in range(0,ny):
    phi=2*pi*j*1.0/ny

    a=(cos(theta))**2+0.25*(sin(theta))**2+0.75*(sin(phi)*cos(theta))**2+0.75*(cos(phi))**2-1-sqrt(3)*0.5*sin(theta)*cos(theta)*sin(phi)

    b=-3.0*cos(theta)*sin(theta)*sin(phi)+0.5*sqrt(3)*(sin(theta))**2-0.5*sqrt(3)*(sin(phi)*cos(theta))**2+sqrt(3)*0.5*(cos(phi))**2

    c=(sin(phi)*sin(theta))**2+0.75*(sin(theta))**2+0.25*(sin(phi)*sin(theta))**2+sqrt(3)*0.5*sin(theta)*cos(theta)*sin(phi)+0.25*(cos(phi))**2-1

    aa=0.25*(cos(theta))**2+(sin(theta))**2+0.75*(sin(phi)*sin(theta))**2-1-sqrt(3)*0.5*sin(theta)*cos(theta)*sin(phi)-1+0.75*(cos(phi))**2  ## a prime

    bb=sqrt(3)*0.5*((sin(phi)*sin(theta))**2-(cos(theta))**2+(cos(phi))**2)+3*sin(theta)*cos(theta)*sin(phi)

    cc=0.75*(cos(theta))**2+0.25*(sin(phi)*sin(theta))**2+(sin(phi)*cos(theta))**2+0.25*(cos(phi))**2-1+sqrt(3)*0.5*sin(theta)*cos(theta)*sin(phi)
    FF[i][j]=(a*cc-aa*c)**2+(b*cc-bb*c)*(aa*b-a*bb)

for i1 in range(1,nx-1):
  for j1 in range(1,ny-1):
    if (FF[i1][j1]*FF[i1+1][j1]<0) or (FF[i1][j1]*FF[i1][j1+1]<0) or (FF[i1][j1]*FF[i1-1][j1]<0) or (FF[i1][j1]*FF[i1][j1-1]<0):
       t[l]=2*pi*i1*1.0/nx
       s[l]=2*pi*j1*1.0/nx
       ff.write(str(2*pi*i1*1.0/nx))
       ff.write('\t')
       ff.write(str(2*pi*j1*1.0/nx))
       ff.write('\n')
       l=l+1
xx=[0]*l
yy=[0]*l
zz=[0]*l
for i2 in range(0,l):
    xx[i2]=sin(t[i2])*cos(s[i2])
    fftt.write(str(xx[i2])+'\t')
    yy[i2]=-cos(t[i2])*cos(s[i2])
    fftt.write(str(yy[i2])+'\t')
    zz[i2]=sin(s[i2])
    fftt.write(str(zz[i2])+'\n')

plot(t,s, 'r_')
plt.xlabel('theta')
plt.ylabel('phi')
plt.show()
fig = plt.figure(1)
ax = p3.Axes3D(fig)
ax.scatter3D(xx,yy,zz)
ax.set_xlabel('X',fontsize=10)
ax.set_ylabel('Y',fontsize=10)
ax.set_zlabel('Z',fontsize=10)
x=arange(-1,1,0.5)
y=arange(-1,1,0.5)
z=arange(-1,1,0.5)
plt.show()


