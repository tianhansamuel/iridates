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
nx=200
ny=200
FF=[[0 for i in range(nx+1)] for j in range(ny+1)]
CC=2.0
SS=(nx+1)*(ny+1)
t=[0]*SS
s=[0]*SS
l1=0
FF1=[[[0 for i in range(nx+1)] for j in range(ny+1)] for k in range(nx+1)]
FF2=[[[0 for i in range(nx+1)] for j in range(ny+1)] for k in range(nx+1)]
ff=open('theta_phi.dat','w')
fftt=open('sphere_axis.dat','w')
for l in range (0,nx+1):
 alpha=2*pi*l*1.0/nx
 for i in range(0,nx+1):
   theta=2*pi*i*1.0/nx
   for j in range(0,ny+1):
     phi=2*pi*j*1.0/ny
     N1=(cos(theta)*cos(alpha)-sin(alpha)*sin(phi)*sin(theta))**2+(cos(alpha+2*pi*1.0/3)*sin(theta)+sin(alpha+2*pi*1.0/3)*sin(phi)*cos(theta))**2+(sin(alpha-2*pi*1.0/3)*cos(phi))**2
     N2=(cos(alpha-2*pi*1.0/3)*cos(theta)-sin(alpha-2*pi*1.0/3)*sin(phi)*sin(theta))**2+(cos(alpha)*sin(theta)+sin(alpha)*sin(phi)*cos(theta))**2+(cos(phi)*sin(alpha+2*pi*1.0/3))**2
     FF1[l][i][j]=N1-1
     FF2[l][i][j]=N2-1


for l in range(1,nx-1): 

 for i1 in range(1,nx-1):
  for j1 in range(1,ny-1):
    if (FF1[l][i1][j1]*FF1[l][i1+1][j1]<0) or (FF1[l][i1][j1]*FF1[l][i1][j1+1]<0) or (FF1[l][i1][j1]*FF1[l][i1-1][j1]<0) or (FF1[l][i1][j1]*FF1[l][i1][j1-1]<0):
     if (FF2[l][i1][j1]*FF2[l][i1+1][j1]<0) or (FF2[l][i1][j1]*FF2[l][i1][j1+1]<0) or (FF2[l][i1][j1]*FF2[l][i1-1][j1]<0) or (FF2[l][i1][j1]*FF2[l][i1][j1-1]<0):
       t[l1]=2*pi*i1*1.0/nx
       s[l1]=2*pi*j1*1.0/nx
       ff.write(str(2*pi*i1*1.0/nx))
       ff.write('\t')
       ff.write(str(2*pi*j1*1.0/nx))
       ff.write('\n')
       l1=l1+1
     
xx=[0]*l1
yy=[0]*l1
zz=[0]*l1
for i2 in range(0,l):
    xx[i2]=sin(t[i2])*cos(s[i2])
    fftt.write(str(xx[i2])+'\t')
    yy[i2]=-cos(t[i2])*cos(s[i2])
    fftt.write(str(yy[i2])+'\t')
    zz[i2]=sin(s[i2])
    fftt.write(str(zz[i2])+'\n')
    

plot(t,s, 'r')
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

