import numpy as np

N = 108
R = np.zeros(N * 3)

def zero(V):
    V = np.zeros(N*3)
    V = np.reshape(3,N)
    v = np.sum(V)
    V = V - v / N
    return V

R = np.loadtxt('R.txt', delimiter=' ')
V = np.loadtxt('V.txt', delimiter=' ')

