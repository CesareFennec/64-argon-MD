from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from numpy.linalg import norm
import time

from math import sqrt
from numba import jit,vectorize,int64
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

def Force(a):
    r = norm(a)
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

def main():
    tb=time.time()
    N = 108
    Steps = 1600

    from Config import R #initial config position
    from Config import V #velocity
    RDF=np.zeros(500)
    norm=np.linalg.norm
    M = 1
    dt = 0.0025
    S = []  # steps
    P = []  # positional
    ff = []  # force constants
    ke = []  # kinetics
    p = []  # potential
    e = []  # energy
    rrr=[]
    vvv=[]

    V_vcc=np.zeros((Steps//4+1,324))
    Coeff=np.zeros(100)
    F = np.zeros(N * 3)
# Velocity Verlet Algorithm
    for step in range(Steps):
        V = V + F * dt / (2 * M)
        R = R + V * dt

    # in this loop time mooo forward
        for val in range(N):  # Velocity and position update according to acc
            for j in range(3):
                # positional pbc begins
                while R[val + j * N] > L:
                    R[val + j * N] = R[val + j * N] - L
                while R[val + j * N] < 0:
                    R[val + j * N] = R[val + j * N] + L
                pass
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
                d = norm(a)
                P = P + LJ(d)
                F[val + 0 * N] = F[val + 0 * N] + a[0] * Force(d) / Z(d)
                F[val + 1 * N] = F[val + 1 * N] + a[1] * Force(d) / Z(d)
                F[val + 2 * N] = F[val + 2 * N] + a[2] * Force(d) / Z(d)
                if step % 4==1 and d<0.5*L and d !=0:
                    b=((d/L)*1000) // 1
                    b=int(b)
                    RDF[b]=RDF[b]+1  #rdf modules
        A = P / 2
        p.append(A)

        V= V+ F * dt / (2 * M)
        # velocity 2nd update
        KE=0.5*sum(V*V)
        ke.append(KE)

        if step%4==0:
            s = step // 4
            V_vcc[s] = V
        else:pass


        E = KE + P / 2
        print(KE)
        e.append(E)
        S.append(step)
        print(step)

    for j in range(100):
        for i in range(Steps//4+1-j):
            Coeff[j]=Coeff[j]+np.sum(V_vcc[i]*V_vcc[j+i])/np.sum(V_vcc[i]*V_vcc[i])




#DATA WRITE MODULE


    ama = max(ke)
    ami = min(ke)
    sig = np.var(ke)
    ava=np.average(ke)

    print(ama, ami, sig, ava)



    for i in range(100):
        Coeff[i]=Coeff[i]/(Steps//4-i)

    c=np.arange(-0.01,1,0.01)
    y=np.zeros(100)
    fig = plt.figure()
    ax0 = fig.add_subplot(111)
    lns0 = ax0.plot(c,Coeff,'red', label='Cvv')
    lns1 =ax0.plot(c,y, 'blue')
    lns = lns1
    ax0.set_xlabel('steps 0.0025 tao per step')
    ax0.set_ylabel("energy emma=118.7K*kb")
    plt.legend()
    plt.show()


    te=time.time()
    print(te-tb,'s')


main()
