from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from numpy import linalg as LA
from math import exp

wid=20
t=1
t1=0.5
N=0
N2=75
N3=75
LX=86
LY=100
SIZE=200
ff=open('mott_transition.dat','w')
ff1=open('mott_transition_1.dat','w')
ff2=open('mott_transition_2.dat','w')
yy=[0]*(N2+1)
t1=0
UC=1.8336923761
hh1=0
hh2=0
ff.write(str(t1)+'\t')
ff.write(str(UC))
ff.write('\n')
for i1 in range(0,N2+N3+1):
   t1=0
   y1=20-(20-UC)*i1*1.0/(N2+N3)
   ff1.write(str(t1)+'\t')
   ff1.write(str(y1)+'\t')       
   zz1=-1
   ff1.write(str(zz1)+'\n')     
   
for ii in range(1,SIZE+1):
   t1=2.0/SIZE*ii
   QQx=0
   QQxx=0
   xixi=10
   summ=0
   N=0
   print t1
   ff.write(str(t1))
   ff.write('\t')
   for i1 in range(0,2*LX):
      kx=2*pi*(i1-LX)/(3*LX)         
      ymin=max(-4*pi/(3*sqrt(3))-1.0/(sqrt(3))*kx,-4*pi/(3*sqrt(3))+1.0/(sqrt(3))*kx)
      ymax=min(4*pi/(3*sqrt(3))-1.0/(sqrt(3))*kx,4*pi/(3*sqrt(3))+1.0/(sqrt(3))*kx)
      for i2 in range(0,2*LY):
          ky=(ymax-ymin)*i2/(2*LY)+ymin
          N=N+1
          gg=t**2*(3+2*cos(sqrt(3)*ky)+4*cos(sqrt(3)*1.0/2*ky)*cos(1.5*kx))
          mrz=(2*t1*sin(1.5*kx+sqrt(3)*1.0/2*ky))**2
          mrxy=((2*t1*sin(1.5*kx+sqrt(3)*1.0/2*ky))**2+(2*t1*sin(1.5*kx-sqrt(3)*1.0/2*ky))**2)
          r=2*t1*(-sin(sqrt(3)*ky)+2*cos(1.5*kx)*sin(sqrt(3)*0.5*ky))
          EE=sqrt(gg+mrz+mrxy)
          QQx=QQx+gg/(EE*t)
          QQxx=QQxx-((2*t1*sin(1.5*kx+sqrt(3)/2*ky))**2/EE+(2*t1*sin(1.5*kx-sqrt(3)/2*ky))**2/EE+(2*t1*sin(sqrt(3)*ky))**2/EE)/(2*t1)
   QQx=QQx/(3*N)
   QQxx=QQxx/(3*N)
   
   for i1 in range(0,2*LX):
      kx=2*pi*(i1-LX)/(3*LX)         
      ymin=max(-4*pi/(3*sqrt(3))-1.0/(sqrt(3))*kx,-4*pi/(3*sqrt(3))+1.0/(sqrt(3))*kx)
      ymax=min(4*pi/(3*sqrt(3))-1.0/(sqrt(3))*kx,4*pi/(3*sqrt(3))+1.0/(sqrt(3))*kx)
      for i2 in range(0,2*LY):
          ky=(ymax-ymin)*i2/(2*LY)+ymin
          g2=2*cos(sqrt(3)*ky)+4*cos(sqrt(3)*1.0/2*ky)*cos(1.5*kx)
          gg=t**2*(3+2*cos(sqrt(3)*ky)+4*cos(sqrt(3)*1.0/2*ky)*cos(1.5*kx))
          xi=-QQx*sqrt(gg)+QQxx*t1*g2
          if xi<xixi:
             xixi=xi

   for i1 in range(0,2*LX):
      kx=2*pi*(i1-LX)/(3*LX)         
      ymin=max(-4*pi/(3*sqrt(3))-1.0/(sqrt(3))*kx,-4*pi/(3*sqrt(3))+1.0/(sqrt(3))*kx)
      ymax=min(4*pi/(3*sqrt(3))-1.0/(sqrt(3))*kx,4*pi/(3*sqrt(3))+1.0/(sqrt(3))*kx)
      for i2 in range(0,2*LY):
          ky=(ymax-ymin)*i2/(2*LY)+ymin
          g2=2*cos(sqrt(3)*ky)+4*cos(sqrt(3)*1.0/2*ky)*cos(1.5*kx)
          gg=t**2*(3+2*cos(sqrt(3)*ky)+4*cos(sqrt(3)/2*ky)*cos(1.5*kx))
          xi=-QQx*sqrt(gg)+QQxx*t1*g2
          if xi>xixi:
              summ=summ+1/sqrt(4*(xi-xixi))
   UC=(N/summ)**2
   
   print UC
   h=2*(17-1.83)
   c=1.83*2-17
   UM=h*1.0/(1+exp(-1.5*t1))+c
   for l2 in range(0,N3):
       ff1.write(str(t1)+'\t')
       yy[l2]=20-(20-UM)*l2*1.0/N3
       ff1.write(str(yy[l2])+'\t')       
       zz1=-1+(t1)*l2*1.0/N3
       zz2=1-(2-t1)*l2*1.0/N3  
       if t1<1+0.5/SIZE:
         ff1.write(str(zz1)+'\n')     
       else:
         ff1.write(str(zz2)+'\n')  
   for ll in range(0,N2+1):
       ff1.write(str(t1)+'\t')
       yy[ll]=UM-(UM-UC)*ll*1.0/N2
       ff1.write(str(yy[ll])+'\t')
       zz=exp(t1**3*(UC-20)*ll*0.01/N2)*(t1-1)#t1-0.2*ll*(t1-1)/N2
       ff1.write(str(zz)+'\n')        
   ff1.write('\n') 
   ff.write(str(UC)+'\n')
ff.close()
