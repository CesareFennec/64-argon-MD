from __future__ import division
import numpy as np
from numpy import linalg as la
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import rc
from math import sqrt
import time
rc('mathtext', default='regular')

L = 5.173403878  #cubic box side length rou*sigma^3=0.78


def LJ(r):  # r is a variable of distance between 2 particles
    if r == 0:
        P = 0
    elif r > L / 2:  # potential cutoff
        P = 0
    else:
        P = 4 * (r ** -12 - r ** -6)  # P is potential between 2 particles
    return P


def Force(r):
    r = la.norm(r)
    if r == 0:
        F = 0
    elif r > L / 2: #Force cutoff
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

tbegin=time.time()
N = 108
Steps =200000
R = np.zeros(N * 3)
from Config import R #initial config position
from Config import V #velocity
M = 1
dt = 0.0025
S = []  # steps
P = []  # positional
ff = []  # force constants
ke = []  # kinetics
p = []  # potential
e = []  # energy

F = np.zeros(N * 3)
RDF=np.zeros(500)
#EXAM KINETIC ENERGY
KE=0
for i in range(N*3):
    KE=KE+0.5*V[i]*V[i]
print(KE)

# Velocity Verlet Algorithm
for step in range(Steps):
    # in this loop time mooo forward
    for val in range(N):  # Velocity and position update according to acc
        for j in range(3):
            V[val + j * N] = V[val + j * N] + F[val + j * N] * dt / (2 * M)
            R[val + j * N] = R[val + j * N] + V[val + j * N] * dt
            # positional pbc begins
            while R[val + j * N] > L:
                R[val + j * N] = R[val + j * N] - L
            while R[val + j * N] < 0:
                R[val + j * N] = R[val + j * N] + L
            continue
            # PBC of position ends
    #this loop updated velocity and position

    F = np.zeros(N * 3)  # Force update loop
    P = 0
    # reconstruct force loop
    for val in range(N):  # we gonna process force on val^{th} particle
        for i in range(N):  # it is the i^{th} particle's force on the val^{th} particle
            a = np.zeros(3)
            for j in range(3):
                #force pbc: pseudo vector 'a' generation
                if R[val + j * N] - R[i + j * N] > L / 2:
                    a[j] = R[val + j * N] - R[i + j * N] - L
                elif R[val + j * N] - R[i + j * N] < -L / 2:
                    a[j] = R[val + j * N] - R[i + j * N] + L
                else:
                    a[j] = R[val + j * N] - R[i + j * N]
                # a pbc precess of i^{th} particle 's position to val
            d = la.norm(a)
            P = P + LJ(d)
            F[val + 0 * N] = F[val + 0 * N] + a[0] * Force(d) / Z(d)
            F[val + 1 * N] = F[val + 1 * N] + a[1] * Force(d) / Z(d)
            F[val + 2 * N] = F[val + 2 * N] + a[2] * Force(d) / Z(d)

            if d>0 and d!=0 and step%4==1:
			    b=(d/L)*1000//1
			    b=int(b)
			    RDF[b]=RDF[b]+1


    A = P / 2
    p.append(A)

    KE = 0
    for val in range(N):  # velocity 2nd update
        for j in range(3):
            V[val + N * j] = V[val + N * j] + F[val + N * j] * dt / (2 * M)
    for i in range(N*3):
        KE=KE+0.5*V[i]*V[i]
    ke.append(KE)
    E = KE + P / 2
    print(KE)
    e.append(E)
    S.append(step)
    print(step)


ama = max(ke)
ami = min(ke)
sig = np.var(ke)
ava=np.average(ke)
print(ama, ami, sig, ava)

tend=time.time()
print((tend-tbegin)/Steps)

sumRDF=np.sum(RDF)
RDF=RDF/sumRDF
fig=plt.figure()

