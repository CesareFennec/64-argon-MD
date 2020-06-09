import numpy as np
from numpy import linalg as la
from numpy import random as ra
np.set_printoptions(threshold=np.inf)
N=16 #systemic particle number
R=np.zeros((N,3))

for a in range(N):
    R[a][0] = ra.uniform(4)
    R[a][1] = ra.uniform(4)
    R[a][2] = ra.uniform(4)
    #give a^t{th} particle a set of position

    b=0
    while b<a:
        if b==0:
            b=b+1
            if la.norm(R[a]-R[b])<1:
                R[a][0] = ra.uniform(4)
                R[a][1] = ra.uniform(4)
                R[a][2] = ra.uniform(4)
                b = 0
        elif la.norm(R[a]-R[b])<1:
            R[a][0] = ra.uniform(4)
            R[a][1] = ra.uniform(4)
            R[a][2] = ra.uniform(4)
            b=0
        else:b=b+1
print(R)

D=np.zeros((N,N))
for c in range(N):
    for d in range(N):
        D[c][d]=la.norm(R[c]-R[d])
print(D)