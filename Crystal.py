import numpy as np
V=np.zeros((4,4,4,3))
R=np.zeros((64,3))
r1=np.array([1,0,0])
r2=np.array([0.5,(3/4)**0.5,0])
r3=np.array([0.5,(1/36)**0.5,(2/3)**0.5])
#fcc with order initial position

for i in range(4):
    for j in range(4):
        for k in range(4):
            V[i][j][k]=[i,j,k]

A=np.reshape(V,(64,3))
print(A)
for val in range(64):
    for a in range(3):
        R[val][a]=A[val][a]*(r1[a]+r2[a]+r3[a])
print(R)
