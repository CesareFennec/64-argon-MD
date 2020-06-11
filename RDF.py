import numpy as np
from Dynamics import R
from Dynamics import N
from numpy import linalg as la
from matplotlib import pyplot as plt
D=np.zeros((N,N))
d=[]
for i in range(N):
    for j in range(i):
       D[i][j]=la.norm(R[i]-R[j])
       if D[i][j]!=0:
           d.append(D[i][j])

plt.hist(d,bins=80,facecolor="blue",edgecolor="black")
plt.xlabel("Radius")
plt.ylabel("counting")
plt.title("RDF")
plt.show()