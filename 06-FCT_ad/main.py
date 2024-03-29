# all quantities are dimensionless

import matplotlib.pyplot as plt
import numpy as np

# mesh options

I_max = 100
N_max = 1000

# u = {rho*S, rho*v*S, rho*energy*S, H_phi*S}
# g(u) = {0, -S*d/dz(p + H_phi^2 / 2), -p*d/dz(v*S), 0}

# u_mesh_1 for rho*S, u_mesh_2 for rho*v*S, u_mesh_3 for rho*energy*S, u_mesh_4 for H_phi*S

u_td_1 = np.zeros((N_max + 1, I_max + 1))
u_td_2 = np.zeros((N_max + 1, I_max + 1))
u_td_3 = np.zeros((N_max + 1, I_max + 1))
u_td_4 = np.zeros((N_max + 1, I_max + 1))

u_ad_1 = np.zeros((N_max + 1, I_max + 1))
u_ad_2 = np.zeros((N_max + 1, I_max + 1))
u_ad_3 = np.zeros((N_max + 1, I_max + 1))
u_ad_4 = np.zeros((N_max + 1, I_max + 1))

# constants

k = 0.5
gamma = 1.67
beta = 1
mu_0 = 0.7

# constants on current layer

mu_left = np.zeros(I_max + 1)
mu_right = np.zeros(I_max + 1)
C_left = np.zeros(I_max + 1)
C_0 = np.zeros(I_max + 1)
C_right = np.zeros(I_max + 1)

mu_ad_left = np.zeros(I_max + 1)
mu_ad_right = np.zeros(I_max + 1)

# create meshes for parameters of problem

rho_mesh = np.zeros(I_max + 1)
v_mesh = np.zeros(I_max + 1)
e_mesh = np.zeros(I_max + 1)
p_mesh = np.zeros(I_max + 1)
H_phi_mesh = np.zeros(I_max + 1)
S_mesh = np.zeros(I_max + 1)

# initial condition

rho_0 = 1
v_0 = 0.1
e_0 = beta / (2 * (gamma - 1))
p_0 = beta / 2 * rho_0 ** gamma

# grid steps

h = 1 / I_max
tau = k * h / (2 * mu_0 + v_0)

# filling meshes for parameters

for i in range(I_max + 1):
    rho_mesh[i] = rho_0
    v_mesh[i] = v_0
    e_mesh[i] = e_0
    p_mesh[i] = p_0
    H_phi_mesh[i] = 1 - 0.9 * i * h
    S_mesh[i] = 2 * (i * h - 0.5) ** 2 + 0.5

# filling meshes for u

for i in range(I_max + 1):
    u_td_1[0, i] = rho_mesh[i] * S_mesh[i]
    u_td_2[0, i] = rho_mesh[i] * v_mesh[i] * S_mesh[i]
    u_td_3[0, i] = rho_mesh[i] * e_mesh[i] * S_mesh[i]
    u_td_4[0, i] = H_phi_mesh[i] * S_mesh[i]

    u_ad_1[0, i] = rho_mesh[i] * S_mesh[i]
    u_ad_2[0, i] = rho_mesh[i] * v_mesh[i] * S_mesh[i]
    u_ad_3[0, i] = rho_mesh[i] * e_mesh[i] * S_mesh[i]
    u_ad_4[0, i] = H_phi_mesh[i] * S_mesh[i]

# evolution in time

for n in range(N_max):

    # left boundary condition for current layer

    rho_mesh[0] = 1
    v_mesh[0] = u_td_2[n, 1] / (rho_mesh[0] * S_mesh[0])
    e_mesh[0] = beta / (2 * (gamma - 1))
    H_phi_mesh[0] = 1

    u_td_1[n + 1, 0] = rho_mesh[0] * S_mesh[0]
    u_td_2[n + 1, 0] = rho_mesh[0] * v_mesh[0] * S_mesh[0]
    u_td_3[n + 1, 0] = rho_mesh[0] * e_mesh[0] * S_mesh[0]
    u_td_4[n + 1, 0] = H_phi_mesh[0] * S_mesh[0]

    u_ad_1[n + 1, 0] = rho_mesh[0] * S_mesh[0]
    u_ad_2[n + 1, 0] = rho_mesh[0] * v_mesh[0] * S_mesh[0]
    u_ad_3[n + 1, 0] = rho_mesh[0] * e_mesh[0] * S_mesh[0]
    u_ad_4[n + 1, 0] = H_phi_mesh[0] * S_mesh[0]

    # constants for current layer

    for i in range(1, I_max):
        mu_left[i] = mu_0 + (abs(v_mesh[i - 1]) + abs(v_mesh[i])) / 4
        mu_right[i] = mu_0 + (abs(v_mesh[i]) + abs(v_mesh[i + 1])) / 4
        C_left[i] = tau / h * (v_mesh[i - 1] + v_mesh[i]) / 4 + tau / h * mu_left[i]
        C_0[i] = 1 - tau / h * (v_mesh[i + 1] - v_mesh[i - 1]) / 4 - tau / h * (mu_left[i] + mu_right[i])
        C_right[i] = - tau / h * (v_mesh[i] + v_mesh[i + 1]) / 4 + tau / h * mu_right[i]

        mu_ad_left[i] = mu_left[i] - abs(v_mesh[i - 1] + v_mesh[i]) / 4
        mu_ad_right[i] = mu_right[i] - abs(v_mesh[i] + v_mesh[i + 1]) / 4

    # filling central points of grid

    for i in range(1, I_max):
        u_td_1[n + 1, i] = C_left[i] * u_td_1[n, i - 1] + C_0[i] * u_td_1[n, i] + C_right[i] * u_td_1[n, i + 1]

        u_td_2[n + 1, i] = C_left[i] * u_td_2[n, i - 1] + C_0[i] * u_td_2[n, i] + \
                             C_right[i] * u_td_2[n, i + 1] - tau / (2 * h) * (p_mesh[i + 1] - p_mesh[i - 1] +
                                                                                H_phi_mesh[i + 1] ** 2 / 2 -
                                                                                H_phi_mesh[i - 1] ** 2 / 2) * S_mesh[i]

        u_td_3[n + 1, i] = C_left[i] * u_td_3[n, i - 1] + C_0[i] * u_td_3[n, i] + \
                             C_right[i] * u_td_3[n, i + 1] - \
                             tau / (2 * h) * p_mesh[i] * (v_mesh[i + 1] * S_mesh[i + 1] - v_mesh[i - 1] * S_mesh[i - 1])

        u_td_4[n + 1, i] = C_left[i] * u_td_4[n, i - 1] + C_0[i] * u_td_4[n, i] + C_right[i] * u_td_4[n, i + 1]

    # right boundary condition for current layer

    u_td_1[n + 1, I_max] = u_td_1[n + 1, I_max - 1]
    u_td_2[n + 1, I_max] = u_td_2[n + 1, I_max - 1]
    u_td_3[n + 1, I_max] = u_td_3[n + 1, I_max - 1]
    u_td_4[n + 1, I_max] = u_td_4[n + 1, I_max - 1]

    # filling meshes for antidiffusion

    for i in range(1, I_max):
        u_ad_1[n + 1, i] = u_td_1[n + 1, i] * (1 + tau / h * (mu_ad_right[i] + mu_ad_left[i])) -\
                           tau / h * (mu_ad_right[i] * u_td_1[n + 1, i + 1] + mu_ad_left[i] * u_td_1[n + 1, i - 1])

        u_ad_2[n + 1, i] = u_td_2[n + 1, i] * (1 + tau / h * (mu_ad_right[i] + mu_ad_left[i])) - \
                           tau / h * (mu_ad_right[i] * u_td_2[n + 1, i + 1] + mu_ad_left[i] * u_td_2[n + 1, i - 1])

        u_ad_3[n + 1, i] = u_td_3[n + 1, i] * (1 + tau / h * (mu_ad_right[i] + mu_ad_left[i])) - \
                           tau / h * (mu_ad_right[i] * u_td_3[n + 1, i + 1] + mu_ad_left[i] * u_td_3[n + 1, i - 1])

        u_ad_4[n + 1, i] = u_td_4[n + 1, i] * (1 + tau / h * (mu_ad_right[i] + mu_ad_left[i])) - \
                           tau / h * (mu_ad_right[i] * u_td_4[n + 1, i + 1] + mu_ad_left[i] * u_td_4[n + 1, i - 1])

    # right boundary condition for current layer for antidiffusion part

    u_ad_1[n + 1, I_max] = u_td_1[n + 1, I_max - 1]
    u_ad_2[n + 1, I_max] = u_td_2[n + 1, I_max - 1]
    u_ad_3[n + 1, I_max] = u_td_3[n + 1, I_max - 1]
    u_ad_4[n + 1, I_max] = u_td_4[n + 1, I_max - 1]

    # update parameters of problem

    for i in range(I_max + 1):
        rho_mesh[i] = u_ad_1[n + 1, i] / S_mesh[i]
        v_mesh[i] = u_ad_2[n + 1, i] / u_ad_1[n + 1, i]
        e_mesh[i] = u_ad_3[n + 1, i] / u_ad_1[n + 1, i]
        p_mesh[i] = beta / 2 * rho_mesh[i] ** gamma
        H_phi_mesh[i] = u_ad_4[n + 1, i] / S_mesh[i]

    # update tau

    tau = k * h / (2 * mu_0 + np.max(v_mesh))

# checkout

x_lst = np.zeros(I_max + 1)
y_lst = np.zeros(I_max + 1)

for i in range(1, I_max + 1):
    x_lst[i] = x_lst[i - 1] + h

# enter the time layer for which you want to display the values
n = N_max

# for rho

for i in range(I_max + 1):
    y_lst[i] = u_ad_1[n, i] / S_mesh[i]

plt.title('Last time layer for rho')
plt.xlabel('channel length coordinate')
plt.ylabel('rho value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()

# for v

for i in range(I_max + 1):
    y_lst[i] = u_ad_2[n, i] / u_ad_1[n, i]

plt.title('Last time layer for velocity')
plt.xlabel('channel length coordinate')
plt.ylabel('velocity value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()

# for energy

for i in range(I_max + 1):
    y_lst[i] = u_ad_3[n, i] / u_ad_1[n, i]

plt.title('Last time layer for energy')
plt.xlabel('channel length coordinate')
plt.ylabel('energy value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()

# for strength

for i in range(I_max + 1):
    y_lst[i] = u_ad_4[n, i] / S_mesh[i]

plt.title('Last time layer for strength')
plt.xlabel('channel length coordinate')
plt.ylabel('strength H_phi value')
plt.grid()
plt.plot(x_lst, y_lst)
plt.show()
