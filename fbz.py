from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt

N=10
SS=(N+1)*(N+1)*6
t=[0]*SS
s=[0]*SS   
LL=0

for i in range(0,N+1):
 for j in range(0,i):
  px=j*1.0/N-0.5*i/N
  py=sqrt(3)*0.5*i/N
  for j in range(0,6):
   qqx=0.5*px+sqrt(3)*0.5*py
   qqy=-0.5*sqrt(3)*px+0.5*py
   px=qqx
   py=qqy
   t[LL]=px
   s[LL]=py
   LL=LL+1



     
plot(t,s, 'r.')
plt.show()
