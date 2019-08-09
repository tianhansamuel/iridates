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
ff=open('spin_texture_position.dat','w')
K0=1
KX=1
KY=1
KZ=1

NB=20
NB1=15
ttt=[0]*(2*NB+1)
SSX=[0]*(2*NB+1)
SSY=[0]*(2*NB+1)
SSZ=[0]*(2*NB+1)
rrr=[0]*(2*NB+1)

t1=lamb=1

rx=sqrt(3)*0.5
ry=-0.5
ddx=sqrt(3)
ddy=0
for lll in range(1,2*NB+1):
 t1=lamb=3.5/(2*NB+1)*lll
 print t1
 #ry=2*(lll-NB)+1 
 ttt[lll]=t1
 SZ1=SZ2=SX1=SX2=SY1=SY2=0
 NNNN=0
 for i in range(0,NB1):
  for j in range(0,i):
   px=(j*1.0/NB1-0.5*i/NB1)*4*pi/(3*sqrt(3))+0.001/NB1
   py=sqrt(3)*0.5*i/NB1*4*pi/(3*sqrt(3))+0.001/NB1
   for k in range(0,6):
    qqx=0.5*px+sqrt(3)*0.5*py
    qqy=-0.5*sqrt(3)*px+0.5*py
    px=qqx
    py=qqy
   
    for i1 in range(0,NB1):
     for j1 in range(0,i1):
      qx=(j1*1.0/NB1-0.5*i1/NB1)*4*pi/(3*sqrt(3))+0.001/NB1
      qy=sqrt(3)*0.5*i1/NB1*4*pi/(3*sqrt(3))+0.001/NB1
      for k2 in range(0,6):
       qqx=0.5*qx+sqrt(3)*0.5*qy
       qqy=-0.5*sqrt(3)*qx+0.5*qy
       qx=qqx
       qy=qqy
       NNNN=NNNN+1

       kx=px-qx
       ky=py-qy
       h1=cos((sqrt(3)*kx))+1j*sin((sqrt(3)*kx))
       h2=cos(sqrt(3)*0.5*kx+1.5*ky)+1j*sin(sqrt(3)*0.5*kx+1.5*ky)
       h3=cos(-sqrt(3)*0.5*kx+1.5*ky)+1j*sin(-sqrt(3)*0.5*kx+1.5*ky)
       h4=cos(-sqrt(3)*kx)+1j*sin(-sqrt(3)*kx)
       h5=cos(-sqrt(3)*0.5*kx-1.5*ky)+1j*sin(-sqrt(3)*0.5*kx-1.5*ky)
       h6=cos(sqrt(3)*0.5*kx-1.5*ky)+1j*sin(sqrt(3)*0.5*kx-1.5*ky)
       h=6-(h1+h2+h3+h4+h5+h6)
       if h!=0:
        POX=(px+qx)*0.5
        POY=(py+qy)*0.5
	########

	A1=1j*PhiM*(1-h1)/h*(cos(POY)-1j*sin(POY))*(cos((sqrt(3)*kx)/2)-1j*sin((sqrt(3)*kx)/2))
	A2=1j*PhiM*(1-h2)/h*(cos(-sqrt(3)*0.5*POX+0.5*POY)-1j*sin(-sqrt(3)*0.5*POX+0.5*POY))*(cos(sqrt(3)*0.25*kx+1.5*ky/2)-1j*sin(sqrt(3)*0.25*kx+1.5*ky/2))
	A3=1j*PhiM*(1-h3)/h*(cos(-sqrt(3)*0.5*POX-0.5*POY)-1j*sin(-sqrt(3)*0.5*POX-0.5*POY))*(cos(-sqrt(3)*0.25*kx+1.5*ky/2)-1j*sin(-sqrt(3)*0.25*kx+1.5*ky/2))
	A4=1j*PhiM*(1-h4)/h*(cos(sqrt(3)*kx/2)+1j*sin(sqrt(3)*kx/2))*(cos(POY)+1j*sin(POY))
	A5=1j*PhiM*(1-h5)/h*(cos(sqrt(3)*0.5*POX-0.5*POY)-1j*sin(sqrt(3)*0.5*POX-0.5*POY))*(cos(sqrt(3)*0.25*kx+1.5*ky/2)+1j*sin(sqrt(3)*0.25*kx+1.5*ky/2))
	A6=1j*PhiM*(1-h6)/h*(cos(sqrt(3)*0.5*POX+0.5*POY)-1j*sin(sqrt(3)*0.5*POX+0.5*POY))*(cos(sqrt(3)*0.25*kx-1.5*ky/2)-1j*sin(sqrt(3)*0.25*kx-1.5*ky/2))######

	AA2=-PhiM*(2-h3-h4)/(3*h)*(cos(-sqrt(3)*0.25*kx+0.25*ky)-1j*sin(-sqrt(3)*0.25*kx+0.25*ky))*cos(-sqrt(3)*0.5*POX-1.5*POY)
	AA1=-PhiM*(2-h1-h2)/(3*h)*(cos(sqrt(3)*0.25*kx+0.25*ky)-1j*sin(sqrt(3)*0.25*kx+0.25*ky))*cos(-sqrt(3)*0.5*POX+1.5*POY)
	AA3=-PhiM*(2-h5-h6)/(3*h)*(cos(0.5*ky)+1j*sin(0.5*ky))*cos(sqrt(3)*POX)

	AB2=-PhiM*(2-h1-h6)/(3*h)*(cos(sqrt(3)*0.25*kx-0.25*ky)-1j*sin(sqrt(3)*0.25*kx-0.25*ky))*cos(sqrt(3)*0.5*POX+1.5*POY)
	AB1=-PhiM*(2-h4-h5)/(3*h)*(cos(-sqrt(3)*0.25*kx-0.25*ky)-1j*sin(-sqrt(3)*0.25*kx-0.25*ky))*cos(sqrt(3)*0.5*POX-1.5*POY)
        AB3=-PhiM*(2-h2-h3)/(3*h)*(cos(0.5*ky)-1j*sin(0.5*ky))*cos(sqrt(3)*POX)
	
	c1=K0*(A1+A3+A5)
	c2=K0*(A2+A4+A6)
	cx=2*t1*KX*AA2
	cy=2*t1*KY*AA1
	cz=2*t1*KZ*AA3
	cx1=2*t1*KX*AB2
	cy1=2*t1*KY*AB1
	cz1=2*t1*KZ*AB3

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

        HH=HH/(M1+M2)
        REP=P1.dot(HH)
	REP1=REP.dot(P2)
        RR=P3.dot(HH)
        RRR=RR.dot(P4)
        REP1=REP1*(cos(kx*rx+ky*ry)+1j*sin(kx*rx+ky*ry)+cos(kx*(rx-ddx)+ky*(ry-ddy))+1j*sin(kx*(rx-ddx)+ky*(ry-ddy)))+RRR*(cos(kx*rx+ky*ry)-1j*sin(kx*rx+ky*ry)+cos(kx*(rx-ddx)+ky*(ry-ddy))-1j*sin(kx*(rx-ddx)+ky*(ry-ddy)))
	REP2=REP1.dot(sigmaz)
        SZ1=SZ1+(REP2[0][0]+REP2[2][2]).real
	SZ2=SZ2+(REP2[1][1]+REP2[3][3]).real
	REP2=REP1.dot(sigmax)
	SX1=SX1+(REP2[0][0]+REP2[2][2]).real
	SX2=SX2+(REP2[1][1]+REP2[3][3]).real
	REP2=REP1.dot(sigmay)
	SY1=(SY1+REP2[0][0]+REP2[2][2]).real
        SY2=(SY2+REP2[1][1]+REP2[3][3]).real

 SSX[lll]=SX1/NNNN
 SSY[lll]=SY1/NNNN
 SSZ[lll]=SZ1/NNNN
 rrr[lll]=SX1/SZ1
 ff.write(str(ttt[lll])+'\t')
 ff.write(str(SSX[lll])+'\t')
 ff.write(str(SSY[lll])+'\t')
 ff.write(str(SSZ[lll])+'\t')
 ff.write(str(rrr[lll])+'\n')
     
#rrr[0]=rrr[1]


#plot(ttt,rrr,'r')
#plot(ttt,SSY,'g')
#plot(ttt,SSZ,'b')
#plt.grid()
plt.xlabel('t1/t')
#plt.ylabel('spin_polarization')
#plt.show()

#ff.write(str(rx)+'\t')
#ff.write(str(ry)+'\t')
#ff.write(str(SX2/NNNN)+'\t')
#ff.write(str(SY2/NNNN)+'\t')
#ff.write(str(SZ2/NNNN)+'\n')
print rx
print ry
print "SZ1" 
print SZ1/NNNN
print 'SZ2' 
print SZ2/NNNN
print "SX1" 
print SX1/NNNN
print 'SX2' 
print SX2/NNNN
print "SY1" 
print SY1/NNNN
print 'SY2' 
print SY2/NNNN


