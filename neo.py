import numpy as np
from scipy import constants as C
from matplotlib import pyplot as plt
d=188*10**-12
M=6.69719964739627e-26
kb=C.Boltzmann
emma=kb*118.2

sigma=341.9
a=0.94
print()
b=1200/sigma
r=np.arange(a,b,0.001)
P=4*(-r**-6+r**-12)
print(r,P)
plt.title('LJ potential of Argon for ref 1')
plt.plot(r,P)
plt.show()
plt.savefig("1d LJP.PNG")