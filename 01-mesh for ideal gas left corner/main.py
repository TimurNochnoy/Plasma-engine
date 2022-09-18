# all quantities are dimensionless

import matplotlib.pyplot as plt
import numpy as np

# mesh options

N_max = 100     # there will be (N_max + 1) values in the mesh (enter the node number starting from 0)
I_max = 50      # there will be (I_max + 1) values in the mesh (enter the node number starting from 0)
tau = 1 / N_max
h = 1 / I_max

# u_mesh_1 for rho*S, u_mesh_2 for rho*v*S, u_mesh_3 for rho*energy*S
# u = {rho*S, rho*v*S, rho*energy*S}

u_mesh_1 = np.zeros((N_max + 1, I_max + 1))
u_mesh_2 = np.zeros((N_max + 1, I_max + 1))
u_mesh_3 = np.zeros((N_max + 1, I_max + 1))

# f_mesh_1 for rho*v*S, f_mesh_2 for rho*v^2*S, f_mesh_3 for rho*v*energy*S
# f(u) = {rho*v*S, rho*v^2*S, rho*v*energy*S}

f_mesh_1 = np.zeros((N_max + 1, I_max + 1))
f_mesh_2 = np.zeros((N_max + 1, I_max + 1))
f_mesh_3 = np.zeros((N_max + 1, I_max + 1))

# g(u) = {0, -S*dp/dz, -p*d/dz(v*S)}

# initial condition

rho_0 = 1.
v_0 = 0.1
p_0 = 1.
gamma = 1.67

S_mesh = np.zeros(I_max + 1)
v_mesh = np.zeros(I_max + 1)
p_mesh = np.zeros(I_max + 1)

for i in range(I_max + 1):
    S_mesh[i] = 2 * (i * h - 0.5) ** 2 + 0.5

for i in range(I_max + 1):
    v_mesh[i] = v_0

for i in range(I_max + 1):
    p_mesh[i] = p_0

for i in range(I_max + 1):
    u_mesh_1[0, i] = rho_0 * S_mesh[i]
    u_mesh_2[0, i] = rho_0 * v_0 * S_mesh[i]
    u_mesh_3[0, i] = p_0 * S_mesh[i]
    f_mesh_1[0, i] = u_mesh_1[0, i] * v_mesh[i]
    f_mesh_2[0, i] = u_mesh_2[0, i] * v_mesh[i]
    f_mesh_3[0, i] = u_mesh_3[0, i] * v_mesh[i]

# mesh filling

for n in range(N_max):

    for i in range(1, I_max):
        v_mesh[i] = f_mesh_1[n, i] / u_mesh_1[n, i]
    v_mesh[I_max] = f_mesh_1[n, I_max] / u_mesh_1[n, I_max]
    v_mesh[I_max] = v_mesh[I_max - 1]

    for i in range(1, I_max + 1):
        u_mesh_1[n + 1, i] = u_mesh_1[n, i] - tau / h * (f_mesh_1[n, i] - f_mesh_1[n, i - 1])

        f_mesh_1[n + 1, i] = u_mesh_1[n + 1, i] * v_mesh[i]
    u_mesh_1[n + 1, 0] = u_mesh_1[n + 1, 1]
    f_mesh_1[n + 1, 0] = f_mesh_1[n + 1, 1]

    for i in range(1, I_max + 1):
        u_mesh_2[n + 1, i] = u_mesh_2[n, i] - tau / h * \
                             (f_mesh_2[n, i] - f_mesh_2[n, i - 1] + 1 / gamma * S_mesh[i] * (p_mesh[i] - p_mesh[i - 1]))

        f_mesh_2[n + 1, i] = u_mesh_2[n + 1, i] * v_mesh[i]
    u_mesh_2[n + 1, 0] = u_mesh_2[n + 1, 1]
    f_mesh_2[n + 1, 0] = f_mesh_2[n + 1, 1]

    for i in range(1, I_max + 1):
        u_mesh_3[n + 1, i] = u_mesh_3[n, i] - \
                             tau / h * (gamma * (f_mesh_3[n, i] - f_mesh_3[n, i - 1]) -
                                        (gamma - 1) * v_mesh[i] * S_mesh[i] * (p_mesh[i] - p_mesh[i - 1]))

        f_mesh_3[n + 1, i] = u_mesh_3[n + 1, i] * v_mesh[i]
    u_mesh_3[n + 1, 0] = u_mesh_3[n + 1, 1]
    f_mesh_3[n + 1, 0] = f_mesh_3[n + 1, 1]

# print numerical results

'''
print('The time axis is directed down, the coordinate axis is directed to the right.\n')
print('Mesh for rho*S.\n')
print(u_mesh_1)
print('\n')

print('Mesh for rho*v*S.\n')
print(u_mesh_2)
print('\n')

print('Mesh for rho*energy*S.\n')
print(u_mesh_3)
print('\n')
'''

# checkout

x_lst = np.zeros(I_max + 1)
y_lst = np.zeros(I_max + 1)

for i in range(I_max + 1):
    y_lst[i] = u_mesh_1[0, i] / S_mesh[i]

for i in range(1, I_max + 1):
    x_lst[i] = x_lst[i-1] + h

# for rho

'''
plt.title('Zero time layer for rho')
plt.xlabel('channel length coordinate')
plt.ylabel('rho value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()


for i in range(I_max + 1):
    y_lst[i] = u_mesh_1[1, i] / S_mesh[i]

plt.title('First time layer for rho')
plt.xlabel('channel length coordinate')
plt.ylabel('rho value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()

'''

for i in range(I_max + 1):
    y_lst[i] = u_mesh_1[N_max, i] / S_mesh[i]

plt.title('Last time layer for rho')
plt.xlabel('channel length coordinate')
plt.ylabel('rho value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()

# for v

'''
for i in range(I_max + 1):
    y_lst[i] = f_mesh_1[0, i] / u_mesh_1[0, i]

plt.title('Zero time layer for velocity')
plt.xlabel('channel length coordinate')
plt.ylabel('velocity value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()

for i in range(I_max + 1):
    y_lst[i] = f_mesh_1[1, i] / u_mesh_1[1, i]

plt.title('First time layer for velocity')
plt.xlabel('channel length coordinate')
plt.ylabel('velocity value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()
'''

for i in range(I_max + 1):
    y_lst[i] = f_mesh_1[N_max, i] / u_mesh_1[N_max, i]

plt.title('Last time layer for velocity')
plt.xlabel('channel length coordinate')
plt.ylabel('velocity value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()

# for energy

'''
for i in range(I_max + 1):
    y_lst[i] = u_mesh_3[0, i] / (S_mesh[i] * (f_mesh_1[0, i] / u_mesh_1[0, i]))

plt.title('Zero time layer for energy')
plt.xlabel('channel length coordinate')
plt.ylabel('energy value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()

for i in range(I_max + 1):
    y_lst[i] = u_mesh_3[1, i] / (S_mesh[i] * (f_mesh_1[1, i] / u_mesh_1[1, i]))

plt.title('First time layer for energy')
plt.xlabel('channel length coordinate')
plt.ylabel('energy value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()
'''

for i in range(I_max + 1):
    y_lst[i] = u_mesh_3[N_max, i] / (S_mesh[i] * (f_mesh_1[N_max, i] / u_mesh_1[N_max, i]))

plt.title('Last time layer for energy')
plt.xlabel('channel length coordinate')
plt.ylabel('energy value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()
