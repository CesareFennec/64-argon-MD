import numpy as np
from numpy import linalg as la
from scipy import constants as C
from matplotlib import pyplot as plt
import random as ran

# initial positions and velocity
r1x = ran.uniform(0.9, 5)
r1y = ran.uniform(0.9, 5)
r1z = ran.uniform(0.9, 5)
r2x = ran.uniform(0.9, 5)
r2y = ran.uniform(0.9, 5)
r2z = ran.uniform(0.9, 5)
r3x = ran.uniform(0.9, 5)
r3y = ran.uniform(0.9, 5)
r3z = ran.uniform(0.9, 5)
v1x = ran.uniform(0.9, 5)
v1y = ran.uniform(0.9, 5)
v1z = ran.uniform(0.9, 5)
v2x = ran.uniform(0.9, 5)
v2y = ran.uniform(0.9, 5)
v2z = ran.uniform(0.9, 5)
v3x = ran.uniform(0.9, 5)
v3y = ran.uniform(0.9, 5)
v3z = ran.uniform(0.9, 5)
# position vectors
R1 = np.array([r1x, r1y, r1z])
R2 = np.array([r2x, r2y, r2z])
R3 = np.array([r3x, r3y, r3z])
R = np.array([R1, R2, R3])
# velocity vectors
V1 = np.array([v1x, v1y, v1z])
V2 = np.array([v2x, v2y, v2z])
V3 = np.array([v3x, v3y, v3z])
V = np.array([V1, V2, V3])
# interacting vectors of position
D12 = -R1 + R2
D21 = -R2 + R1
D13 = -R1 + R3
D31 = -R3 + R1
D23 = -R2 + R3
D32 = -R3 + R2
# noramlized interaction vectors of position
d12 = D12 / (la.norm(D12))
d21 = D12 / (la.norm(D21))
d13 = D12 / (la.norm(D13))
d31 = D12 / (la.norm(D31))
d23 = D12 / (la.norm(D23))
d32 = D12 / (la.norm(D32))
d1 = la.norm(D12)
d2 = la.norm(D23)
d3 = la.norm(D13)
p1 = d1 ** -12 - d1 ** -6
p2 = d2 ** -12 - d2 ** -6
p3 = d3 ** -12 - d3 ** -6
P = p1 + p2 + p3
# initial force
f1 = 6 * d1 ** -7 - 12 * d1 ** -13
f2 = 6 * d2 ** -7 - 12 * d2 ** -13
f3 = 6 * d3 ** -7 - 12 * d3 ** -13
F12 = f1 * d12
F21 = f1 * d21
F23 = f2 * d23
F32 = f2 * d32
F31 = f3 * d31
F13 = f3 * d13
F1 = F21 + F31
F2 = F12 + F32
F3 = F13 + F23
F = np.array([F1, F2, F3])
m = 1
dt = 0.001
for i in range(1000):
    A = F / m
    V = V + A * dt
    [V1, V2, V3] = np.array(V)
    V1 = np.array([v1x, v1y, v1z])
    V2 = np.array([v2x, v2y, v2z])
    V3 = np.array([v3x, v3y, v3z])
    v1 = la.norm(V1) ** 2
    v2 = la.norm(V2) ** 2
    v3 = la.norm(V3) ** 2
    KE = 0.5 * m * (v1 + v2 + v3)
    R = R + V * dt
    [R1, R2, R3] = np.array(R)
    # interacting vectors of position
    D12 = -R1 + R2
    D21 = -R2 + R1
    D13 = -R1 + R3
    D31 = -R3 + R1
    D23 = -R2 + R3
    D32 = -R3 + R2
    # noramlized interaction vectors of position
    d12 = D12 / (la.norm(D12))
    d21 = D12 / (la.norm(D21))
    d13 = D12 / (la.norm(D13))
    d31 = D12 / (la.norm(D31))
    d23 = D12 / (la.norm(D23))
    d32 = D12 / (la.norm(D32))
    d1 = la.norm(D12)
    d2 = la.norm(D23)
    d3 = la.norm(D13)
    p1 = d1 ** -12 - d1 ** -6
    p2 = d2 ** -12 - d2 ** -6
    p3 = d3 ** -12 - d3 ** -6
    P = p1 + p2 + p3
    # initial force
    f1 = 6 * d1 ** -7 - 12 * d1 ** -13
    f2 = 6 * d2 ** -7 - 12 * d2 ** -13
    f3 = 6 * d3 ** -7 - 12 * d3 ** -13
    F12 = f1 * d12
    F21 = f1 * d21
    F23 = f2 * d23
    F32 = f2 * d32
    F31 = f3 * d31
    F13 = f3 * d13
    F1 = F21 + F31
    F2 = F12 + F32
    F3 = F13 + F23
    F = np.array([F1, F2, F3])
    E = P + KE
    print(P, KE, E)
