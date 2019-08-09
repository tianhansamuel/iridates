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
CC=0.0001
l=0
ff=open('theta_phi.dat','w')
fftt=open('sphere_axis.dat','w')

phi=0.615479709
theta=5*pi/4

for i in range(0,nx+1):
   alpha=2*pi*i*1.0/nx
   N1=(cos(theta)*cos(alpha)-sin(alpha)*sin(phi)*sin(theta))**2+(cos(alpha+2*pi*1.0/3)*sin(theta)+sin(alpha+2*pi*1.0/3)*sin(phi)*cos(theta))**2+(sin(alpha-2*pi*1.0/3)*cos(phi))**2
   N2=(cos(alpha-2*pi*1.0/3)*cos(theta)-sin(alpha-2*pi*1.0/3)*sin(phi)*sin(theta))**2+(cos(alpha)*sin(theta)+sin(alpha)*sin(phi)*cos(theta))**2+(cos(phi)*sin(alpha+2*pi*1.0/3))**2
   if fabs(N1-1)<CC/nx and fabs(N2-1)<CC/nx:
     print alpha
     l=l+1

print l


