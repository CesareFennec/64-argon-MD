import numpy as np
from numpy import linalg as la
from numpy import random as ra
np.set_printoptions(threshold=np.inf)
r=[]
L=4
N=64
points = set()
while len(points) < N:
    p = (ra.randint(0, L), ra.randint(0, L), ra.randint(0, L))
    if p not in points:
        points.add(p)
        r.append(p)
print(r)
R=np.array(r,dtype=float)
print(R)
for val in range(N):
    for j in range(3):
        R[val][j]=R[val][j]+np.random.uniform(0,0.1)-0.05
        while R[val][j]>=4:
            R[val][j]=R[val][j]-4
            continue
        while R[val][j]<0:
            R[val][j]=R[val][j]+4
            continue
print(R)
D=np.zeros((N,N))
for val in range(N):
    for i in range(N):
        D[val][i]=la.norm(R[val]-R[i])
        if D[val][i]<0.81:
            print(val,i)