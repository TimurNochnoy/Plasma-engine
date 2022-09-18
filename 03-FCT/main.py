# all quantities are dimensionless

import matplotlib.pyplot as plt
import numpy as np

# mesh options

time = 0.5
I_max = 50      # there will be (I_max + 1) values in the mesh (enter the node number starting from 0)
tau = 0.01
h = 1 / I_max
N_max = round(time / tau)
print(N_max)

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

v_0 = 0.1
gamma = 1.67
mu_0 = 0.7

S_mesh = np.zeros(I_max + 1)
v_mesh = np.zeros(I_max + 1)
p_mesh = np.zeros(I_max + 1)
rho_mesh = np.zeros(I_max + 1)

for i in range(I_max + 1):
    S_mesh[i] = 2 * (i * h - 0.5) ** 2 + 0.5

for i in range(I_max + 1):
    v_mesh[i] = v_0

for i in range(I_max + 1):
    p_mesh[i] = 1 - 0.8 * i * h

for i in range(I_max + 1):
    rho_mesh[i] = 1 - 0.8 * i * h

for i in range(I_max + 1):
    u_mesh_1[0, i] = rho_mesh[i] * S_mesh[i]
    u_mesh_2[0, i] = rho_mesh[i] * v_0 * S_mesh[i]
    u_mesh_3[0, i] = p_mesh[i] * S_mesh[i]
    f_mesh_1[0, i] = u_mesh_1[0, i] * v_mesh[i]
    f_mesh_2[0, i] = u_mesh_2[0, i] * v_mesh[i]
    f_mesh_3[0, i] = u_mesh_3[0, i] * v_mesh[i]

# mesh filling

for n in range(N_max):

    # filling mesh for u

    for i in range(1, I_max):
        mu_left = mu_0 + (abs(v_mesh[i - 1]) + abs(v_mesh[i])) / 4
        mu_right = mu_0 + (abs(v_mesh[i]) + abs(v_mesh[i + 1])) / 4
        C_left = tau / h * (v_mesh[i - 1] + v_mesh[i]) / 4 + tau / h * mu_left
        C_0 = 1 - tau / h * (v_mesh[i + 1] - v_mesh[i - 1]) / 4 - tau / h * (mu_left + mu_right)
        C_right = - tau / h * (v_mesh[i] + v_mesh[i + 1]) / 4 + tau / h * mu_right

        u_mesh_1[n + 1, i] = C_left * u_mesh_1[n, i - 1] + C_0 * u_mesh_1[n, i] + C_right * u_mesh_1[n, i + 1]
    u_mesh_1[n + 1, 0] = rho_mesh[0] * S_mesh[0]
    u_mesh_1[n + 1, I_max] = u_mesh_1[n + 1, I_max - 1]

    for i in range(1, I_max):
        mu_left = mu_0 + (abs(v_mesh[i - 1]) + abs(v_mesh[i])) / 4
        mu_right = mu_0 + (abs(v_mesh[i]) + abs(v_mesh[i + 1])) / 4
        C_left = tau / h * (v_mesh[i - 1] + v_mesh[i]) / 4 + tau / h * mu_left
        C_0 = 1 - tau / h * (v_mesh[i + 1] + v_mesh[i - 1]) / 4 - tau / h * (mu_left + mu_right)
        C_right = - tau / h * (v_mesh[i] + v_mesh[i + 1]) / 4 + tau / h * mu_right

        u_mesh_2[n + 1, i] = C_left * u_mesh_2[n, i - 1] + C_0 * u_mesh_2[n, i] + C_right * u_mesh_2[n, i + 1] -\
                             tau / (2 * h) * 1 / gamma * S_mesh[i] * (p_mesh[i + 1] - p_mesh[i - 1])
    u_mesh_2[n + 1, 0] = rho_mesh[0] * v_mesh[0] * S_mesh[0]
    u_mesh_2[n + 1, I_max] = u_mesh_2[n + 1, I_max - 1]

    for i in range(1, I_max):
        mu_left = mu_0 + (abs(v_mesh[i - 1]) + abs(v_mesh[i])) / 4
        mu_right = mu_0 + (abs(v_mesh[i]) + abs(v_mesh[i + 1])) / 4
        C_left = tau / h * (v_mesh[i - 1] + v_mesh[i]) / 4 + tau / h * mu_left
        C_0 = 1 - tau / h * (v_mesh[i + 1] + v_mesh[i - 1]) / 4 - tau / h * (mu_left + mu_right)
        C_right = - tau / h * (v_mesh[i] + v_mesh[i + 1]) / 4 + tau / h * mu_right

        u_mesh_3[n + 1, i] = C_left * u_mesh_3[n, i - 1] + C_0 * u_mesh_3[n, i] + C_right * u_mesh_3[n, i + 1] +\
                             tau / (2 * h) * (gamma - 1) * v_mesh[i] * S_mesh[i] * (p_mesh[i + 1] - p_mesh[i - 1])
    u_mesh_3[n + 1, 0] = p_mesh[0] * S_mesh[0]
    u_mesh_3[n + 1, I_max] = u_mesh_3[n + 1, I_max - 1]

    # calculate the velocity as u_2 / u_1

    v_mesh[0] = u_mesh_2[n + 1, 1]  # proven analytically
    for i in range(1, I_max + 1):
        v_mesh[i] = u_mesh_2[n + 1, i] / u_mesh_1[n + 1, i]

    # filling mesh for f(u)

    for i in range(I_max + 1):
        f_mesh_1[n + 1, i] = u_mesh_1[n + 1, i] * v_mesh[i]

    for i in range(I_max + 1):
        f_mesh_2[n + 1, i] = u_mesh_2[n + 1, i] * v_mesh[i]

    for i in range(I_max + 1):
        f_mesh_3[n + 1, i] = u_mesh_3[n + 1, i] * v_mesh[i]

# checkout

x_lst = np.zeros(I_max + 1)
y_lst = np.zeros(I_max + 1)

for i in range(1, I_max + 1):
    x_lst[i] = x_lst[i-1] + h

# enter the time layer for which you want to display the values
n = N_max

# for rho

for i in range(I_max + 1):
    y_lst[i] = u_mesh_1[n, i] / S_mesh[i]

plt.title('Last time layer for rho')
plt.xlabel('channel length coordinate')
plt.ylabel('rho value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()

# for v

for i in range(I_max + 1):
    y_lst[i] = f_mesh_1[n, i] / u_mesh_1[n, i]

plt.title('Last time layer for velocity')
plt.xlabel('channel length coordinate')
plt.ylabel('velocity value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()

# for energy

for i in range(I_max + 1):
    y_lst[i] = u_mesh_3[n, i] / ((gamma - 1) * u_mesh_1[n, i])

plt.title('Last time layer for energy')
plt.xlabel('channel length coordinate')
plt.ylabel('energy value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()
