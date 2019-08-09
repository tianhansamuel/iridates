from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from numpy import linalg as LA
from math import exp


rx=-sqrt(3)*0.5
ry=-0.5


t=1
LX=10
LY=10
PhiM=2*pi
P1=zeros([4,4],'complex')
P2=zeros([4,4],'complex')
P3=zeros([4,4],'complex')
P4=zeros([4,4],'complex')
HH=zeros([4,4],'complex')
sigmax=zeros([4,4],'complex')
sigmay=zeros([4,4],'complex')
sigmaz=zeros([4,4],'complex')
sigmaz[0][0]=sigmaz[1][1]=1
sigmaz[2][2]=sigmaz[3][3]=-1
sigmax[0][2]=sigmax[1][3]=sigmax[2][0]=sigmax[3][1]=1
sigmay[0][2]=sigmay[1][3]=-1j
sigmay[2][0]=sigmay[3][1]=1j

K0=1
KX=1
KY=1
KZ=1
SZ1=SZ2=SX1=SX2=SY1=SY2=0
NNNN=0
NB=10
SSX=[0]*NB
SSY=[0]*NB
SSZ=[0]*NB
for lll in range(0,NB):
 t1=lamb=3.0/NB*lll
 for i1 in range(0,2*LX+1):
  py=2*pi*(i1-LX)/(3*LX)      
  ymin=-4*pi/(3*sqrt(3))+1/(sqrt(3))*fabs(py)
  ymax=4*pi/(3*sqrt(3))-1/(sqrt(3))*fabs(py)
  for i2 in range(0,4*LY+1):
   px=2*pi*i2/(3*LX)+ymin
   if px<ymax:
    for i3 in range(0,2*LX+1):
     qy=2*pi*(i3-LX)/(3*LX)      
     ymin=-4*pi/(3*sqrt(3))+1/(sqrt(3))*fabs(qy)
     ymax=4*pi/(3*sqrt(3))-1/(sqrt(3))*fabs(qy)
     for i4 in range(0,4*LY+1):
      qx=2*pi*i4/(3*LX)+ymin
      if qx<ymax:
       kx=px-qx
       ky=py-qy
       if kx*ky!=0:
        POX=(px+qx)*0.5
        POY=(py+qy)*0.5
	h1=cos((sqrt(3)*kx))+1j*sin((sqrt(3)*kx))
	h2=cos(sqrt(3)*0.5*kx+1.5*ky)+1j*sin(sqrt(3)*0.5*kx+1.5*ky)
	h3=cos(-sqrt(3)*0.5*kx+1.5*ky)+1j*sin(-sqrt(3)*0.5*kx+1.5*ky)
	h4=cos(-sqrt(3)*kx)+1j*sin(-sqrt(3)*kx)
	h5=cos(-sqrt(3)*0.5*kx-1.5*ky)+1j*sin(-sqrt(3)*0.5*kx-1.5*ky)
	h6=cos(sqrt(3)*0.5*kx-1.5*ky)+1j*sin(sqrt(3)*0.5*kx-1.5*ky)
	h=6-(h1+h2+h3+h4+h5+h6)########

	A1=2*1j*PhiM*(1-h1)/h*(cos(POY)-1j*sin(POY))*(cos((sqrt(3)*kx)/2)-1j*sin((sqrt(3)*kx)/2))
	A2=2*1j*PhiM*(1-h2)/h*(cos(-sqrt(3)*0.5*POX+0.5*POY)-1j*sin(-sqrt(3)*0.5*POX+0.5*POY))*(cos(sqrt(3)*0.25*kx+1.5*ky/2)-1j*sin(sqrt(3)*0.25*kx+1.5*ky/2))
	A3=2*1j*PhiM*(1-h3)/h*(cos(-sqrt(3)*0.5*POX-0.5*POY)-1j*sin(-sqrt(3)*0.5*POX-0.5*POY))*(cos(-sqrt(3)*0.25*kx+1.5*ky/2)-1j*sin(-sqrt(3)*0.25*kx+1.5*ky/2))
	A4=2*1j*PhiM*(1-h4)/h*(cos(sqrt(3)*kx/2)+1j*sin(sqrt(3)*kx/2))*(cos(POY)+1j*sin(POY))
	A5=2*1j*PhiM*(1-h5)/h*(cos(sqrt(3)*0.5*POX-0.5*POY)-1j*sin(sqrt(3)*0.5*POX-0.5*POY))*(cos(sqrt(3)*0.25*kx+1.5*ky/2)+1j*sin(sqrt(3)*0.25*kx+1.5*ky/2))
	A6=2*1j*PhiM*(1-h6)/h*(cos(sqrt(3)*0.5*POX+0.5*POY)-1j*sin(sqrt(3)*0.5*POX+0.5*POY))*(cos(sqrt(3)*0.25*kx-1.5*ky/2)-1j*sin(sqrt(3)*0.25*kx-1.5*ky/2))######

	AA2=-PhiM*(2-h3-h4)/(3*h)*(cos(-sqrt(3)*0.25*kx+0.25*ky)-1j*sin(-sqrt(3)*0.25*kx+0.25*ky))*cos(-sqrt(3)*0.5*POX-1.5*POY)
	AA1=-PhiM*(2-h1-h2)/(3*h)*(cos(sqrt(3)*0.25*kx+0.25*ky)-1j*sin(sqrt(3)*0.25*kx+0.25*ky))*cos(-sqrt(3)*0.5*POX+1.5*POY)
	AA3=-PhiM*(2-h5-h6)/(3*h)*(cos(0.5*ky)+1j*sin(0.5*ky))*cos(sqrt(3)*POX)

	AB2=-PhiM*(2-h1-h6)/(3*h)*(cos(sqrt(3)*0.25*kx-0.25*ky)-1j*sin(sqrt(3)*0.25*kx-0.25*ky))*cos(sqrt(3)*0.5*POX+1.5*POY)
	AB1=-PhiM*(2-h4-h5)/(3*h)*(cos(-sqrt(3)*0.25*kx-0.25*ky)-1j*sin(-sqrt(3)*0.25*kx-0.25*ky))*cos(sqrt(3)*0.5*POX-1.5*POY)
        AB3=-PhiM*(2-h2-h3)/(3*h)*(cos(0.5*ky)-1j*sin(0.5*ky))*cos(sqrt(3)*POX)
	
	c1=K0*(A1+A3+A5)
	c2=K0*(A2+A4+A6)
	cx=t1*KX*AA2
	cy=t1*KY*AA1
	cz=t1*KZ*AA3
	cx1=t1*KX*AB2
	cy1=t1*KY*AB1
	cz1=t1*KZ*AB3

	m1=cos(py)+cos(sqrt(3)*0.5*px-0.5*py)+cos(-sqrt(3)*0.5*px-0.5*py)
	m2=-(sin(py)+sin(sqrt(3)*0.5*px-0.5*py)+sin(-sqrt(3)*0.5*px-0.5*py))
	mx=2*t1*sin(-px*sqrt(3)*0.5-1.5*py)
	my=2*t1*sin(-px*sqrt(3)*0.5+1.5*py)
        mz=2*t1*sin(px*sqrt(3))
	M1=sqrt(m1**2+m2**2+mx**2+my**2+mz**2)
        m1=m1/M1
        m2=m2/M1
        mx=mx/M1
        my=my/M1
        mz=mz/M1

	mm1=cos(qy)+cos(sqrt(3)*0.5*qx-0.5*qy)+cos(-sqrt(3)*0.5*qx-0.5*qy)
	mm2=-(sin(qy)+sin(sqrt(3)*0.5*qx-0.5*qy)+sin(-sqrt(3)*0.5*qx-0.5*qy))
	mmx=2*t1*sin(-qx*sqrt(3)*0.5-1.5*qy)
	mmy=2*t1*sin(-qx*sqrt(3)*0.5+1.5*qy)
        mmz=2*t1*sin(qx*sqrt(3))
	M2=sqrt(mmx**2+mmy**2+mmz**2+mm1**2+mm2**2)
        mm1=mm1/M2
        mm2=mm2/M2
        mmx=mmx/M2
        mmy=mmy/M2
        mmz=mmz/M2

	P1[0][0]=P1[3][3]=1-mmz
	P1[1][1]=P1[2][2]=1+mmz
	P1[0][1]=P1[2][3]=-mm1+1j*mm2
	P1[1][0]=P1[3][2]=-mm1-1j*mm2
	P1[2][0]=-mmx-1j*mmy
	P1[3][1]=mmx+1j*mmy
	P1[0][2]=-mmx+1j*mmy
	P1[1][3]=mmx-1j*mmy

	P2[0][0]=P2[3][3]=1+mz
	P2[1][1]=P2[2][2]=1-mz
	P2[0][1]=P2[2][3]=m1-1j*m2
	P2[1][0]=P2[3][2]=m1+1j*m2
	P2[2][0]=mx+1j*my
	P2[3][1]=-mx-1j*my
	P2[0][2]=mx-1j*my
	P2[1][3]=-mx+1j*my
	
	P3[0][0]=P3[3][3]=1-mz
	P3[1][1]=P3[2][2]=1+mz
	P3[0][1]=P3[2][3]=-m1+1j*m2
	P3[1][0]=P3[3][2]=-m1-1j*m2
	P3[2][0]=-mx-1j*my
	P3[3][1]=mx+1j*my
	P3[0][2]=-mx+1j*my
	P3[1][3]=mx-1j*my

	P4[0][0]=P4[3][3]=1+mmz
	P4[1][1]=P4[2][2]=1-mmz
	P4[0][1]=P4[2][3]=mm1-1j*mm2
	P4[1][0]=P4[3][2]=mm1+1j*mm2
	P4[2][0]=mmx+1j*mmy
	P4[3][1]=-mmx-1j*mmy
	P4[0][2]=mmx-1j*mmy
	P4[1][3]=-mmx+1j*mmy

	HH[0][0]=cz
	HH[2][2]=-cz
	HH[1][1]=-cz1
	HH[3][3]=cz1
	HH[1][0]=HH[3][2]=c1
	HH[0][1]=HH[2][3]=c2
	HH[2][0]=cx+1j*cy
	HH[3][1]=-cx1-1j*cy1
	HH[0][2]=cx-1j*cy	
	HH[1][3]=-cx1+1j*cy1

        HH=HH*(cos(kx*rx+ky*ry)-1j*sin(kx*rx+ky*ry))/(M1+M2)
        REP=P1.dot(HH)
	REP1=REP.dot(P2)
        RR=P3.dot(HH)
        RRR=RR.dot(P4)
        REP1=REP1+RRR
	REP2=REP1.dot(sigmaz)
	SZ1=(SZ1+REP2[0][0]+REP2[2][2])
	SZ2=(SZ2+REP2[1][1]+REP2[3][3])
	REP2=REP1.dot(sigmax)
	SX1=(SX1+REP2[0][0]+REP2[2][2])
	SX2=(SX2+REP2[1][1]+REP2[3][3])
	REP2=REP1.dot(sigmay)
	SY1=(SY1+REP2[0][0]+REP2[2][2])
	SY2=(SY2+REP2[1][1]+REP2[3][3])
	NNNN=NNNN+1
 SSX[lll]=SX1
 SSY[lll]=SY1
 SSZ[lll]=SZ1

plot(SSX,'r')
plot(SSY,'g')
plot(SSZ,'b')
plt.show()
#print rx
#print ry
#print "SZ1" 
#print SZ1/NNNN
#print 'SZ2' 
#print SZ2/NNNN
#print "SX1" 
#print SX1/NNNN
#print 'SX2' 
#print SX2/NNNN
#print "SY1" 
#print SY1/NNNN
#print 'SY2' 
#print SY2/NNNN
