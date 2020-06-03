import numpy as np
from numpy import linalg as la

np.seterr(divide='ignore', invalid='ignore')
N = 64  # particle numbers
Steps = 10  # simulate steps


def LJ(r):
    if r == 0:
        p = 0
    else:
        p = 4 * (-r ** -6 + r ** -12)
    return p


def force(r):
    if r == 0:
        f = 0
    else:
        f = 24 * (r ** -7 + r ** -13)
    return f


R = np.zeros((N, 3))  # vectors
V = np.zeros((N, 3))  # sqr velocity
F = np.zeros((N, 3))
M = np.zeros(N)
dt = 0



for val in range(N):
    M[val] = 1
    for i in range(N):
        for j in range(3):
            F[val][j] = F[val][j] + (R[i]-R[val])[j]*force(la.norm(R[i] - R[val])) / (la.norm(R[i] - R[val]))
            V[val][j] = V[val][j] + dt * F[val][j] / (2 * M[val])
            R[val][j] = R[val][j] + dt * V[val][j]
            F[val][j] = F[val][j] + LJ(la.norm(R[val] - R[i])) * (R[i][j] - R[val][j]) / (la.norm(R[val] - R[i]))
            V[val][j] = V[val][j] + dt * F[val][j] / (2 * M[val])
