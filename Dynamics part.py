import numpy as np
from numpy import linalg as la
import random as ran
from matplotlib import pyplot as plt

def LJ(r): #r is a variable of distance between 2 particles
    if r==0:
        P=0 #P is potential between 2 particles
    else:P=4*(r**-12-r**-6)
    return P
def Force(r):
    r=la.norm(r)
    if r==0:
        F=0
    else:F=-24*(-2*r**-13+r**-7)
    return F
def Z(r):
    if r==0:
        a=1
    else:a=r
    return a
N=4
Steps=1


R=np.zeros((N,3))
R[0]=([0,0,0])
R[1]=([1,0,0])
R[2]=([0,1,0])
R[3]=([1,1,0])
V = np.zeros((N, 3))
M=1
dt=0.01
S=[]
P=[]
ff=[]



for step in range(Steps): #Verlet Algorithm
    F=np.zeros((N,3))
    for val in range (N): #Force update according to position
        for i in range(N):
            d =la.norm(R[i]-R[val])
            for j in range(3):
                F[val][j] = F[val][j] + (R[val][j] - R[i][j]) * Force(d) / Z(d)

    for val in range(N): #Velocity and position update according to acc
        for j in range (3):
            V[val][j]=V[val][j]+F[val][j]*dt/(2*M)
            R[val][j]=R[val][j]+V[val][j]*dt
    print(V)
    F = np.zeros((N, 3))#Force update 2nd time
    for val in range(N):
        for i in range(N):
            d = la.norm(R[i] - R[val])
            for j in range(3):
                F[val][j] = F[val][j] + (R[val][j] - R[i][j]) * Force(d) / Z(d)

    for val in range(N):#velocity 2nd update
        for j in range(3):
            V[val][j]=V[val][j]+F[val][j]*dt/(2*M)
    S.append(step)
    P.append(R[0][1])
    ff.append(F[0][1])
    print(V)
print(ff,S,P)

fig=plt.figure()
ax1=fig.add_subplot()
ax1.plot(S,P,'red')
ax1.set_ylabel('Positional sample')

ax2=ax1.twinx()
ax2.plot(S,ff,'blue')
ax2.set_ylabel("Force sample")

