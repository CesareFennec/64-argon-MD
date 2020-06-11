import random as ran
import numpy as np
from numpy import linalg as la

N = 64  # particle numbers
Steps = 20  # simulate steps
x = np.zeros(N)  # coordinates
y = np.zeros(N)
z = np.zeros(N)
vx = np.zeros(N)  # velocity
vy = np.zeros(N)
vz = np.zeros(N)
R = np.zeros((N, 3))  # vectors
V = np.zeros((N, 3))  # sqr velocity
v = np.zeros(N)  # velocity
D = np.zeros((N, N, 3))
d = np.zeros((N, N, 3))
dd = np.zeros((N, N))
F = np.zeros((N, N))
Fd = np.zeros((N, N, 3))
Fg = np.zeros((N, 3))
Acc = np.zeros((N, 3))
P = np.zeros((N, N))  # potential matrix
np.seterr(divide='ignore', invalid='ignore')
m = 1
dt = 10 ** -3


def LJ(r):
    if r == 0:
        P = 0
    else:
        P = 4 * (-r ** -6 + r ** -12)
    return P


def f(r):
    if r == 0:
        f = 0
    else:
        f = -24 * (-2 * r ** -13 + r ** -6)
    return f


# initial positions generator
for val in range(N):
    x[val] = ran.uniform(0, 4)
    y[val] = ran.uniform(0, 4)
    z[val] = ran.uniform(0, 4)
    R[val] = ([x[val], y[val], z[val]])  # directions
    D[val] = R - R[val]  # dist vevtors
    vx[val] = ran.uniform(1.44, 1.48)
    vy[val] = ran.uniform(1.44, 1.48)
    vz[val] = ran.uniform(1.44, 1.48)
    V[val] = ([vx[val], vy[val], vz[val]])  # velocity
    v[val] = (vx[val] ** 2 + vy[val] ** 2 + vz[val] ** 2)  # sum of square of velocity
    for i in range(N):
        d[val][i] = D[val][i] / (la.norm(D[val][i]))  # normalized distance
        dd[val][i] = la.norm(D[val][i])  # Dist metrix
        F[val][i] = f(dd[val][i])  # Force
        P[val][i] = LJ(dd[val][i])
        Fd[val][i] = F[val][i] * d[val][i]
        Fg[val] = Fg[val] + F[val][i]
        Acc[val] = Fg[val] / m

# step-by-step calc
for step in range(Steps):
    p = 0.5 * np.sum(P)
    k = 0.5 * np.sum(v)
    print(k, p, k + p, R)
    for val in range(N):
        V[val] = Acc[val] * dt + V[val]
        R[val] = R[val] + V[val] * dt
        v[val] = (vx[val] ** 2 + vy[val] ** 2 + vz[val] ** 2)
        D[val] = R - R[val]

        for i in range(N):
            d[val][i] = D[val][i] / (la.norm(D[val][i]))
            dd[val][i] = la.norm(D[val][i])
            P[val][i] = LJ(dd[val][i])
            F[val][i] = f(dd[val][i])
            Fd[val][i] = F[val][i] * d[val][i]
            Fg[val] = Fg[val] + F[val][i]
            Acc[val] = Fg[val] / m
