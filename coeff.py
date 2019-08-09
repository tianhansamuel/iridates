from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt

nx=20
ny=20
EQ=[[0 for i in range(nx)] for j in range(ny)]
J2=2
LX=5
LY=5
SS=10*LX*LY
LL=0
t=[0]*SS
s=[0]*SS
S1=0
S2=0
S3=0
n=0
bb=4   
for i1 in range(0,2*LX+1):
        kx=2*pi*(i1-LX)/(3*LX)      
        ymin=-4*pi/(3*sqrt(3))+1/(sqrt(3))*fabs(kx)
        ymax=4*pi/(3*sqrt(3))-1/(sqrt(3))*fabs(kx)
        for i2 in range(0,5*LY+1):
         ky=2*pi*i2/(3*LX)+ymin
         if ky<ymax:
          g=(3+2*cos(sqrt(3)*ky)+4*cos(sqrt(3)/2*ky)*cos(1.5*kx))
          rx=J2*sin(-3.0/2*kx-sqrt(3)/2*ky)
          ry=J2*sin(3.0/2*kx-sqrt(3)/2*ky)
          rz=J2*sin(sqrt(3)*ky)
          phix=(-3.0/2*kx-sqrt(3)/2*ky)/bb
          phiy=(3.0/2*kx-sqrt(3)/2*ky)/bb
          phiz=sqrt(3)*ky/bb
          rr=J2*(sin(-3.0/2*kx-sqrt(3)/2*ky+phix)+sin(-3.0/2*kx-sqrt(3)/2*ky+phiy)+sin(-3.0/2*kx-sqrt(3)/2*ky+phiz))
          E=sqrt(g+rx**2+ry**2+rz**2)
          S1=S1+rx*rr/E
          S2=S2+ry*rr/E
          S3=S3+rz*rr/E
          n=n+1
print S1/n
print S2/n
print S3/n
          
