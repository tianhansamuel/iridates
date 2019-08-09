from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from numpy import linalg as LA
from math import exp

t=1
t1=lamb=0.1
LX=50
LY=50
EE=zeros([16*LX*LY,1])
KK=zeros([16*LX*LY,1])
QQ=0
QX=0
QY=0
QZ=0
Nombre=0
for i1 in range(0,2*LX):
      kx=2*pi*(i1-LX)/(3*LX)         
      ymin=max(-4*pi/(3*sqrt(3))-1.0/(sqrt(3))*kx,-4*pi/(3*sqrt(3))+1.0/(sqrt(3))*kx)
      ymax=min(4*pi/(3*sqrt(3))-1.0/(sqrt(3))*kx,4*pi/(3*sqrt(3))+1.0/(sqrt(3))*kx)
      for i2 in range(0,2*LY):
          ky=(ymax-ymin)*i2/(2*LY)+ymin
          Nombre=Nombre+1
          g=t*(cos(-kx/2+sqrt(3)/2*ky)+1j*sin(-kx/2+sqrt(3)/2*ky)+cos(-kx/2-sqrt(3)/2*ky)+1j*sin(-kx/2-sqrt(3)/2*ky)+cos(kx)+1j*sin(kx))
          HAM=zeros([4,4],complex)
          rz=cos(kx)+1j*sin(kx)
          rx=cos(-kx/2+sqrt(3)/2*ky)+1j*sin(-kx/2+sqrt(3)/2*ky)
          ry=cos(-kx/2-sqrt(3)/2*ky)+1j*sin(-kx/2-sqrt(3)/2*ky)
          HAM[0][0]=HAM[0][2]=HAM[1][1]=HAM[1][3]=HAM[2][0]=HAM[2][2]=HAM[3][1]=HAM[3][3]=0
          HAM[0][1]=g+1j*lamb*rz
          HAM[0][3]=lamb*(1j*rx+ry)
          HAM[1][0]=g.conjugate()-1j*lamb*rz.conjugate()
          HAM[1][2]=-1j*lamb*rx.conjugate()-lamb*ry.conjugate()
          HAM[2][1]=lamb*(1j*rx-ry)
          HAM[2][3]=g-1j*rz
          HAM[3][0]=lamb*(-1j*rx.conjugate()+ry.conjugate())
          HAM[3][2]=g.conjugate()+1j*lamb*rz.conjugate()
          VAL,VEC=LA.eig(HAM)
          MM=VEC
          for i in range(0,4):
             N=sqrt((MM[i][0].conjugate()*MM[i][0]+MM[i][1].conjugate()*MM[i][1]+MM[i][2].conjugate()*MM[i][2]+MM[i][3].conjugate()*MM[i][3]).real)
             MM[i][0]=MM[i][0]/N
             MM[i][1]=MM[i][1]/N
             MM[i][2]=MM[i][2]/N
             MM[i][3]=MM[i][3]/N
          MM=MM/det(MM)
          T=linalg.inv(MM)
          for i in range(0,4):
             a=T[i][0]
             b=T[i][1]
             c=T[i][2]
             d=T[i][3]
             T[i][1]=a
             T[i][2]=b
             T[i][3]=c
             T[i][0]=d
          for i in range(0,4):
             N=sqrt((T[0][i].conjugate()*T[0][i]+T[1][i].conjugate()*T[1][i]+T[2][i].conjugate()*T[2][i]+T[3][i].conjugate()*T[3][i]).real)
             T[0][i]=T[0][i]/N
             T[1][i]=T[1][i]/N
             T[2][i]=T[2][i]/N
             T[3][i]=T[3][i]/N
          QX=QX+2*(1j*rx*((T[0][0]).conjugate()*T[3][0]+(T[0][1]).conjugate()*T[3][1]+(T[2][0]).conjugate()*T[1][0]+(T[2][0]).conjugate()*T[1][0])).real
          QY=QY+(1j*(ry*(-1j*T[0][0].conjugate()*T[3][0]-1j*T[0][1].conjugate()*T[3][1]+1j*T[2][0].conjugate()*T[1][0]+1j*T[2][1].conjugate()*T[1][1]))-1j*(ry.conjugate()*(1j*T[0][0]*T[3][0].conjugate()+1j*T[0][1]*T[3][1].conjugate()-1j*T[2][0]*T[1][0].conjugate()-1j*T[2][1]*T[1][1].conjugate()))).real
          QZ=QZ+(((1j*rz*(T[0][0].conjugate()*T[1][0]+T[0][1].conjugate()*T[1][1]-T[2][0].conjugate()*T[3][0]-T[2][1].conjugate()*T[3][1])-1j*rz.conjugate()*(T[0][0]*T[1][0].conjugate()+T[0][1]*T[1][1].conjugate()-T[2][0]*T[3][0].conjugate()-T[2][1]*T[3][1].conjugate())))).real
          QQ=QQ+2*(g*T[0][0].conjugate()*T[1][0]+g*T[0][1].conjugate()*T[1][1]+g*T[2][0].conjugate()*T[3][0]+g*T[2][1].conjugate()*T[3][1]).real
QX=2*QX.real/(Nombre)
QY=2*QY.real/(Nombre)
QZ=2*QZ.real/(Nombre)
QQ=2*QQ/(Nombre)
#print QX
#print QY
#print QZ
#print QQ
xixi=10
for i1 in range(0,2*LX):
      kx=2*pi*(i1-LX)/(3*LX)         
      ymin=max(-4*pi/(3*sqrt(3))-1.0/(sqrt(3))*kx,-4*pi/(3*sqrt(3))+1.0/(sqrt(3))*kx)
      ymax=min(4*pi/(3*sqrt(3))-1.0/(sqrt(3))*kx,4*pi/(3*sqrt(3))+1.0/(sqrt(3))*kx)
      for i2 in range(0,2*LY):
          ky=(ymax-ymin)*i2/(2*LY)+ymin
          gx=2*t1*cos(-0.5*kx+sqrt(3)*0.5*ky)
          gy=2*t1*cos(-0.5*kx-sqrt(3)*0.5*ky)
          gz=2*t1*cos(ky)
          gg=sqrt(t**2*(3+2*cos(sqrt(3)*ky)+4*cos(sqrt(3)*1.0/2*ky)*cos(1.5*kx)))
          xi=(QQ*gg+gx*QX+gy*QY+gz*QZ).real
          if xi<xixi:
             xixi=xi.real
summ=0
for i1 in range(0,2*LX):
      kx=2*pi*(i1-LX)/(3*LX)         
      ymin=max(-4*pi/(3*sqrt(3))-1.0/(sqrt(3))*kx,-4*pi/(3*sqrt(3))+1.0/(sqrt(3))*kx)
      ymax=min(4*pi/(3*sqrt(3))-1.0/(sqrt(3))*kx,4*pi/(3*sqrt(3))+1.0/(sqrt(3))*kx)
      for i2 in range(0,2*LY):
          ky=(ymax-ymin)*i2/(2*LY)+ymin
          gx=2*t1*cos(-0.5*kx+sqrt(3)*0.5*ky)
          gy=2*t1*cos(-0.5*kx-sqrt(3)*0.5*ky)
          gz=2*t1*cos(ky)
          gg=sqrt(t**2*(3+2*cos(sqrt(3)*ky)+4*cos(sqrt(3)*1.0/2*ky)*cos(1.5*kx)))
          xi=(QQ*gg+gx*QX+gy*QY+gz*QZ).real
          if xi>xixi:
              summ=summ+1/sqrt(4*(xi-xixi))
print (Nombre/summ)**2

