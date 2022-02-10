# all quantities are dimensionless

import numpy as np

N_max = 10
I_max = 5
tau = 1 / N_max
h = 1 / I_max
mesh = np.zeros((N_max, I_max))

# initial condition
for i in range(I_max):
    mesh[0, i] = (i + 1) * 10.

# mesh filling
for n in range(N_max-1):
    for i in range(I_max-1):
        mesh[n+1, i] = mesh[n, i] - tau / h * (mesh[n, i+1] - mesh[n, i])
    mesh[n+1, I_max-1] = mesh[n+1, I_max-2]

print('Mesh for rho*S.\n')
print('The time axis is directed down, the coordinate axis is directed to the right.\n')
print(mesh)