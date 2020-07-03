import numpy as np
from numpy import random as ra

np.set_printoptions(threshold=np.inf)
r = []
L = 5
N = 108
points = set()
while len(points) < N:
    p = (ra.randint(0, L), ra.randint(0, L), ra.randint(0, L))
    if p not in points:
        points.add(p)
        r.append(p)
print(r)
P = np.array(r, dtype=float)
print(P)

for val in range(N):
    for j in range(3):
        P[val][j] = P[val][j] + ra.uniform(0, 0.1) - 0.05
        while P[val][j] >= L:
            P[val][j] = P[val][j] - L
            continue
        while P[val][j] < 0:
            P[val][j] = P[val][j] + L
            continue
R = np.zeros(N * 3)
for val in range(N):
    for j in range(3):
        R[val + j * N] = P[val][j] * 1.034680776  # coeff for number density
print(R)
