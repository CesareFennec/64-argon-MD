import numpy as np

N = 108
R = np.zeros(N * 3)

def zero(V):
    V = np.reshape(V,(3,N))
    for i in range(3):
        v = np.sum(V[i])
        V[i] = V[i] - v / N
    return V
R = np.loadtxt('R.txt', delimiter=' ')
V = np.loadtxt('V.txt', delimiter=' ')
V=zero(V)
print(sum(V[0]),sum(V[1]),sum(V[2]))
V=np.reshape(V,3*N)
