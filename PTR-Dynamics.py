import numpy as np
from numpy import linalg as la
from matplotlib import pyplot as plt
from matplotlib import rc
from math import sqrt

rc('mathtext', default='regular')


def LJ(r):  # r is a variable of distance between 2 particles
    if r == 0:
        P = 0  # P is potential between 2 particles
    else:
        P = 4 * (r ** -12 - r ** -6)
    return P


def Force(r):
    r = la.norm(r)
    if r == 0 :
        F = 0
    else:
        F = -24 * (-2 * r ** -13 + r ** -7)
    return F


def Z(r):
    if r == 0:
        a = 1
    else:
        a = r
    return a


N = 64
Steps = 5000
R=np.zeros(N*3)
from QuasiCluster import B
for val in range(N):
    for j in range(3):
        R[val + j * N]=B[val][j]
print(R)

V = np.zeros(N * 3)
M = 1
dt = 0.0001
S = []  # steps
P = []  # positional
ff = []  # force constants
ke = []  # kinetics
p = []  # potential
e = []  # energy
F = np.zeros(N * 3)
for step in range(Steps):  # Verlet Algorithm
    for val in range(N):  # Velocity and position update according to acc
        for j in range(3):
            V[val + j * N] = V[val + j * N] + F[val + j * N] * dt / (2 * M)
            R[val + j * N] = R[val + j * N] + V[val + j * N] * dt

    print(V,'V',R,'R',step)

    F = np.zeros(N * 3)  # Force update 2nd time
    P = 0
    for val in range(N):
        for i in range(N):
            d = sqrt((R[i] - R[val]) ** 2 + (R[i + N] - R[val + N]) ** 2 + (R[i + 2 * N] - R[val + 2 * N]) ** 2)
            P = P + LJ(d)
            for j in range(3):
                F[val + N * j] = F[val + N * j] + (R[val + N * j] - R[i + N * j]) * Force(d) / Z(d)
    print(F,'F',step)
    A = P / 2
    p.append(A)
    KE=0
    for val in range(N):  # velocity 2nd update
        for j in range(3):
            V[val + N * j] = V[val + N * j] + F[val + N * j] * dt / (2 * M)
        KE = KE + 0.5 * M * ((V[val]) ** 2 + (V[val + N * 1]) ** 2 + (V[val + N * 2]) ** 2)
    ke.append(KE)
    E = KE + P / 2
    print(V,'V',step)
    print(step, KE, P, E)
    e.append(E)
    S.append(step)
    ff.append(F[0])
    
    #this version works on quad-cluster with D4h symmetry
