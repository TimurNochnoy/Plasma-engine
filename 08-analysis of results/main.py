# all quantities are dimensionless

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def algorithm(I_max, N_max, k, gamma, beta, mu_0, rho_0, v_0, a, b):
    # initial condition

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
        S_mesh[i] = a * (i * h - 0.5) ** 2 + b

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

        u_cor_1[0, i] = rho_mesh[i] * S_mesh[i]
        u_cor_2[0, i] = rho_mesh[i] * v_mesh[i] * S_mesh[i]
        u_cor_3[0, i] = rho_mesh[i] * e_mesh[i] * S_mesh[i]
        u_cor_4[0, i] = H_phi_mesh[i] * S_mesh[i]

    # evolution in time

    for n in range(N_max):

        # left boundary condition for current layer

        rho_mesh[0] = 1
        v_mesh[0] = u_cor_2[n, 1] / (rho_mesh[0] * S_mesh[0])
        e_mesh[0] = beta / (2 * (gamma - 1))
        H_phi_mesh[0] = 1

        # left boundary condition for current layer for transport part

        u_td_1[n + 1, 0] = rho_mesh[0] * S_mesh[0]
        u_td_2[n + 1, 0] = rho_mesh[0] * v_mesh[0] * S_mesh[0]
        u_td_3[n + 1, 0] = rho_mesh[0] * e_mesh[0] * S_mesh[0]
        u_td_4[n + 1, 0] = H_phi_mesh[0] * S_mesh[0]

        # left boundary condition for current layer for antidiffusion part

        u_ad_1[n + 1, 0] = rho_mesh[0] * S_mesh[0]
        u_ad_2[n + 1, 0] = rho_mesh[0] * v_mesh[0] * S_mesh[0]
        u_ad_3[n + 1, 0] = rho_mesh[0] * e_mesh[0] * S_mesh[0]
        u_ad_4[n + 1, 0] = H_phi_mesh[0] * S_mesh[0]

        # left boundary condition for current layer for correction part

        u_cor_1[n + 1, 0] = rho_mesh[0] * S_mesh[0]
        u_cor_2[n + 1, 0] = rho_mesh[0] * v_mesh[0] * S_mesh[0]
        u_cor_3[n + 1, 0] = rho_mesh[0] * e_mesh[0] * S_mesh[0]
        u_cor_4[n + 1, 0] = H_phi_mesh[0] * S_mesh[0]

        # constants for current layer

        for i in range(1, I_max):
            mu_left[i] = mu_0 + (abs(v_mesh[i - 1]) + abs(v_mesh[i])) / 4
            mu_right[i] = mu_0 + (abs(v_mesh[i]) + abs(v_mesh[i + 1])) / 4
            C_left[i] = tau / h * (v_mesh[i - 1] + v_mesh[i]) / 4 + tau / h * mu_left[i]
            C_0[i] = 1 - tau / h * (v_mesh[i + 1] - v_mesh[i - 1]) / 4 - tau / h * (mu_left[i] + mu_right[i])
            C_right[i] = - tau / h * (v_mesh[i] + v_mesh[i + 1]) / 4 + tau / h * mu_right[i]

            mu_ad_left[i] = mu_left[i] - abs(v_mesh[i - 1] + v_mesh[i]) / 4
            mu_ad_right[i] = mu_right[i] - abs(v_mesh[i] + v_mesh[i + 1]) / 4

        # filling central points of grid for transport part

        for i in range(1, I_max):
            u_td_1[n + 1, i] = C_left[i] * u_td_1[n, i - 1] + C_0[i] * u_td_1[n, i] + C_right[i] * u_td_1[n, i + 1]

            u_td_2[n + 1, i] = C_left[i] * u_td_2[n, i - 1] + C_0[i] * u_td_2[n, i] + \
                               C_right[i] * u_td_2[n, i + 1] - tau / (2 * h) * (p_mesh[i + 1] - p_mesh[i - 1] +
                                                                                H_phi_mesh[i + 1] ** 2 / 2 -
                                                                                H_phi_mesh[i - 1] ** 2 / 2) * S_mesh[i]

            u_td_3[n + 1, i] = C_left[i] * u_td_3[n, i - 1] + C_0[i] * u_td_3[n, i] + \
                               C_right[i] * u_td_3[n, i + 1] - \
                               tau / (2 * h) * p_mesh[i] * (
                                       v_mesh[i + 1] * S_mesh[i + 1] - v_mesh[i - 1] * S_mesh[i - 1])

            u_td_4[n + 1, i] = C_left[i] * u_td_4[n, i - 1] + C_0[i] * u_td_4[n, i] + C_right[i] * u_td_4[n, i + 1]

        # right boundary condition for current layer for transport part

        u_td_1[n + 1, I_max] = u_td_1[n + 1, I_max - 1]
        u_td_2[n + 1, I_max] = u_td_2[n + 1, I_max - 1]
        u_td_3[n + 1, I_max] = u_td_3[n + 1, I_max - 1]
        u_td_4[n + 1, I_max] = u_td_4[n + 1, I_max - 1]

        # filling central points of grid for antidiffusion part

        for i in range(1, I_max):
            u_ad_1[n + 1, i] = u_td_1[n + 1, i] * (1 + tau / h * (mu_ad_right[i] + mu_ad_left[i])) - \
                               tau / h * (mu_ad_right[i] * u_td_1[n + 1, i + 1] + mu_ad_left[i] * u_td_1[n + 1, i - 1])

            u_ad_2[n + 1, i] = u_td_2[n + 1, i] * (1 + tau / h * (mu_ad_right[i] + mu_ad_left[i])) - \
                               tau / h * (mu_ad_right[i] * u_td_2[n + 1, i + 1] + mu_ad_left[i] * u_td_2[n + 1, i - 1])

            u_ad_3[n + 1, i] = u_td_3[n + 1, i] * (1 + tau / h * (mu_ad_right[i] + mu_ad_left[i])) - \
                               tau / h * (mu_ad_right[i] * u_td_3[n + 1, i + 1] + mu_ad_left[i] * u_td_3[n + 1, i - 1])

            u_ad_4[n + 1, i] = u_td_4[n + 1, i] * (1 + tau / h * (mu_ad_right[i] + mu_ad_left[i])) - \
                               tau / h * (mu_ad_right[i] * u_td_4[n + 1, i + 1] + mu_ad_left[i] * u_td_4[n + 1, i - 1])

        # right boundary condition for current layer for antidiffusion part

        u_ad_1[n + 1, I_max] = u_ad_1[n + 1, I_max - 1]
        u_ad_2[n + 1, I_max] = u_ad_2[n + 1, I_max - 1]
        u_ad_3[n + 1, I_max] = u_ad_3[n + 1, I_max - 1]
        u_ad_4[n + 1, I_max] = u_ad_4[n + 1, I_max - 1]

        # filling meshes for correction part parameters

        for i in range(2, I_max - 1):
            s = np.sign(u_td_1[n + 1, i + 1] - u_td_1[n + 1, i])
            f_cor_left_1[i] = s * max(0, min(s * (u_td_1[n + 1, i + 1] - u_td_1[n + 1, i]),
                                             abs(tau / h * mu_ad_left[i] * (u_td_1[n + 1, i] - u_td_1[n + 1, i - 1])),
                                             s * (u_td_1[n + 1, i - 1] - u_td_1[n + 1, i - 2])))
            f_cor_right_1[i] = s * max(0, min(s * (u_td_1[n + 1, i + 2] - u_td_1[n + 1, i + 1]),
                                              abs(tau / h * mu_ad_right[i] * (u_td_1[n + 1, i + 1] - u_td_1[n + 1, i])),
                                              s * (u_td_1[n + 1, i] - u_td_1[n + 1, i - 1])))
            s = np.sign(u_td_2[n + 1, i + 1] - u_td_2[n + 1, i])
            f_cor_left_2[i] = s * max(0, min(s * (u_td_2[n + 1, i + 1] - u_td_2[n + 1, i]),
                                             abs(tau / h * mu_ad_left[i] * (u_td_2[n + 1, i] - u_td_2[n + 1, i - 1])),
                                             s * (u_td_2[n + 1, i - 1] - u_td_2[n + 1, i - 2])))
            f_cor_right_2[i] = s * max(0, min(s * (u_td_2[n + 1, i + 2] - u_td_2[n + 1, i + 1]),
                                              abs(tau / h * mu_ad_right[i] * (u_td_2[n + 1, i + 1] - u_td_2[n + 1, i])),
                                              s * (u_td_2[n + 1, i] - u_td_2[n + 1, i - 1])))
            s = np.sign(u_td_3[n + 1, i + 1] - u_td_3[n + 1, i])
            f_cor_left_3[i] = s * max(0, min(s * (u_td_3[n + 1, i + 1] - u_td_3[n + 1, i]),
                                             abs(tau / h * mu_ad_left[i] * (u_td_3[n + 1, i] - u_td_3[n + 1, i - 1])),
                                             s * (u_td_3[n + 1, i - 1] - u_td_3[n + 1, i - 2])))
            f_cor_right_3[i] = s * max(0, min(s * (u_td_3[n + 1, i + 2] - u_td_3[n + 1, i + 1]),
                                              abs(tau / h * mu_ad_right[i] * (u_td_3[n + 1, i + 1] - u_td_3[n + 1, i])),
                                              s * (u_td_3[n + 1, i] - u_td_3[n + 1, i - 1])))
            s = np.sign(u_td_4[n + 1, i + 1] - u_td_4[n + 1, i])
            f_cor_left_4[i] = s * max(0, min(s * (u_td_4[n + 1, i + 1] - u_td_4[n + 1, i]),
                                             abs(tau / h * mu_ad_left[i] * (u_td_4[n + 1, i] - u_td_4[n + 1, i - 1])),
                                             s * (u_td_4[n + 1, i - 1] - u_td_4[n + 1, i - 2])))
            f_cor_right_4[i] = s * max(0, min(s * (u_td_4[n + 1, i + 2] - u_td_4[n + 1, i + 1]),
                                              abs(tau / h * mu_ad_right[i] * (u_td_4[n + 1, i + 1] - u_td_4[n + 1, i])),
                                              s * (u_td_4[n + 1, i] - u_td_4[n + 1, i - 1])))

        # filling central points of grid for correction part

        for i in range(2, I_max - 1):
            u_cor_1[n + 1, i] = u_td_1[n + 1, i] + (f_cor_left_1[i] - f_cor_right_1[i])
            u_cor_2[n + 1, i] = u_td_2[n + 1, i] + (f_cor_left_2[i] - f_cor_right_2[i])
            u_cor_3[n + 1, i] = u_td_3[n + 1, i] + (f_cor_left_3[i] - f_cor_right_3[i])
            u_cor_4[n + 1, i] = u_td_4[n + 1, i] + (f_cor_left_4[i] - f_cor_right_4[i])

        # boundary condition for current layer for correction part

        u_cor_1[n + 1, 1] = u_ad_1[n + 1, 1]
        u_cor_2[n + 1, 1] = u_ad_2[n + 1, 1]
        u_cor_3[n + 1, 1] = u_ad_3[n + 1, 1]
        u_cor_4[n + 1, 1] = u_ad_4[n + 1, 1]

        u_cor_1[n + 1, I_max - 1] = u_ad_1[n + 1, I_max - 1]
        u_cor_2[n + 1, I_max - 1] = u_ad_2[n + 1, I_max - 1]
        u_cor_3[n + 1, I_max - 1] = u_ad_3[n + 1, I_max - 1]
        u_cor_4[n + 1, I_max - 1] = u_ad_4[n + 1, I_max - 1]

        u_cor_1[n + 1, I_max] = u_cor_1[n + 1, I_max - 1]
        u_cor_2[n + 1, I_max] = u_cor_2[n + 1, I_max - 1]
        u_cor_3[n + 1, I_max] = u_cor_3[n + 1, I_max - 1]
        u_cor_4[n + 1, I_max] = u_cor_4[n + 1, I_max - 1]

        # update parameters of problem

        for i in range(I_max + 1):
            rho_mesh[i] = u_cor_1[n + 1, i] / S_mesh[i]
            v_mesh[i] = u_cor_2[n + 1, i] / u_cor_1[n + 1, i]
            e_mesh[i] = u_cor_3[n + 1, i] / u_cor_1[n + 1, i]
            p_mesh[i] = beta / 2 * rho_mesh[i] ** gamma
            H_phi_mesh[i] = u_cor_4[n + 1, i] / S_mesh[i]

        # update tau

        tau = k * h / (2 * mu_0 + np.max(v_mesh))


# mesh options

I_max_main = 100
N_max_main = 1000
k_main = 0.5
gamma_main = 1.67
beta_main = 1
mu_0_main = 0.7
rho_0_main = 1
v_0_main = 0.1
a_main = 2
b_main = 0.5

# u = {rho*S, rho*v*S, rho*energy*S, H_phi*S}
# g(u) = {0, -S*d/dz(p + H_phi^2 / 2), -p*d/dz(v*S), 0}

# u_mesh_1 for rho*S, u_mesh_2 for rho*v*S, u_mesh_3 for rho*energy*S, u_mesh_4 for H_phi*S

u_td_1 = np.zeros((N_max_main + 1, I_max_main + 1))
u_td_2 = np.zeros((N_max_main + 1, I_max_main + 1))
u_td_3 = np.zeros((N_max_main + 1, I_max_main + 1))
u_td_4 = np.zeros((N_max_main + 1, I_max_main + 1))

u_ad_1 = np.zeros((N_max_main + 1, I_max_main + 1))
u_ad_2 = np.zeros((N_max_main + 1, I_max_main + 1))
u_ad_3 = np.zeros((N_max_main + 1, I_max_main + 1))
u_ad_4 = np.zeros((N_max_main + 1, I_max_main + 1))

u_cor_1 = np.zeros((N_max_main + 1, I_max_main + 1))
u_cor_2 = np.zeros((N_max_main + 1, I_max_main + 1))
u_cor_3 = np.zeros((N_max_main + 1, I_max_main + 1))
u_cor_4 = np.zeros((N_max_main + 1, I_max_main + 1))

# create meshes for correction part parameters

f_cor_left_1 = np.zeros(I_max_main + 1)
f_cor_right_1 = np.zeros(I_max_main + 1)
f_cor_left_2 = np.zeros(I_max_main + 1)
f_cor_right_2 = np.zeros(I_max_main + 1)
f_cor_left_3 = np.zeros(I_max_main + 1)
f_cor_right_3 = np.zeros(I_max_main + 1)
f_cor_left_4 = np.zeros(I_max_main + 1)
f_cor_right_4 = np.zeros(I_max_main + 1)

# constants on current layer

mu_left = np.zeros(I_max_main + 1)
mu_right = np.zeros(I_max_main + 1)
C_left = np.zeros(I_max_main + 1)
C_0 = np.zeros(I_max_main + 1)
C_right = np.zeros(I_max_main + 1)

mu_ad_left = np.zeros(I_max_main + 1)
mu_ad_right = np.zeros(I_max_main + 1)

# create meshes for parameters of problem

rho_mesh = np.zeros(I_max_main + 1)
v_mesh = np.zeros(I_max_main + 1)
e_mesh = np.zeros(I_max_main + 1)
p_mesh = np.zeros(I_max_main + 1)
H_phi_mesh = np.zeros(I_max_main + 1)
S_mesh = np.zeros(I_max_main + 1)

# algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_0_main)

# checkout

x_lst = np.zeros(I_max_main + 1)
y_lst = np.zeros(I_max_main + 1)

for i in range(1, I_max_main + 1):
    x_lst[i] = x_lst[i - 1] + 1 / I_max_main

# enter the time layer for which you want to display the values
N = N_max_main

# animation of the establishment of a stationary regime of plasma flow in the channel

algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_0_main, a_main, b_main)

# for rho

fig = plt.figure()
images = []

for n in range(0, N_max_main + 1, 5):
    plt.title('Plasma density')
    plt.xlabel('channel length coordinate')
    plt.ylabel('rho value')
    for i in range(I_max_main + 1):
        y_lst[i] = u_cor_1[n, i] / S_mesh[i]
    image = plt.plot(x_lst, y_lst, 'g')
    plt.grid()
    images.append(image)

ani = animation.ArtistAnimation(fig, images, interval=1, repeat_delay=100)
ani.save("02 - establishment of rho.gif", writer='pillow')

# for velocity

fig = plt.figure()
images = []

for n in range(0, N_max_main + 1, 5):
    plt.title('Plasma velocity')
    plt.xlabel('channel length coordinate')
    plt.ylabel('velocity value')
    for i in range(I_max_main + 1):
        y_lst[i] = u_cor_2[n, i] / u_cor_1[n, i]
    image = plt.plot(x_lst, y_lst, 'b')
    plt.grid()
    images.append(image)

ani = animation.ArtistAnimation(fig, images, interval=1, repeat_delay=100)
ani.save("02 - establishment of velocity.gif", writer='pillow')

# for energy

fig = plt.figure()
images = []

for n in range(0, N_max_main + 1, 5):
    plt.title('Plasma energy')
    plt.xlabel('channel length coordinate')
    plt.ylabel('energy value')
    for i in range(I_max_main + 1):
        y_lst[i] = u_cor_3[n, i] / u_cor_1[n, i]
    image = plt.plot(x_lst, y_lst, 'r')
    plt.grid()
    images.append(image)

ani = animation.ArtistAnimation(fig, images, interval=1, repeat_delay=100)
ani.save("02 - establishment of energy.gif", writer='pillow')

# for strength

fig = plt.figure()
images = []

for n in range(0, N_max_main + 1, 5):
    plt.title('Azimuthal magnetic field')
    plt.xlabel('channel length coordinate')
    plt.ylabel('strength value')
    for i in range(I_max_main + 1):
        y_lst[i] = u_cor_4[n, i] / S_mesh[i]
    image = plt.plot(x_lst, y_lst, 'y')
    plt.grid()
    images.append(image)

ani = animation.ArtistAnimation(fig, images, interval=1, repeat_delay=100)
ani.save("02 - establishment of strength.gif", writer='pillow')
'''
# checking that the velocities intersect at the center

plt.title('Plasma and fast magnetosonic velocity')
plt.xlabel('channel length coordinate')
plt.ylabel('velocity value')

algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_2[N, i] / u_cor_1[N, i]
plt.plot(x_lst, y_lst, label='Plasma velocity')

for i in range(I_max_main + 1):
    y_lst[i] = np.sqrt(H_phi_mesh[i] ** 2 / rho_mesh[i] + gamma_main * p_mesh[i] / rho_mesh[i])
plt.plot(x_lst, y_lst, label='Fast magnetosonic velocity')

plt.grid()
plt.legend()
plt.show()

# geometry of channel

plt.title('Geometry of channel')
plt.xlabel('z coordinate')
plt.ylabel('r coordinate')

for i in range(I_max_main + 1):
    y_lst[i] = 3.6 * (i * 1 / I_max_main - 0.5) ** 2 + 0.1
plt.plot(x_lst, y_lst, label='S = 3.6(z - 0.5)^2 + 0.1')

for i in range(I_max_main + 1):
    y_lst[i] = 2 * (i * 1 / I_max_main - 0.5) ** 2 + 0.5
plt.plot(x_lst, y_lst, label='S = 2(z - 0.5)^2 + 0.5')

for i in range(I_max_main + 1):
    y_lst[i] = 0.8 * (i * 1 / I_max_main - 0.5) ** 2 + 0.8
plt.plot(x_lst, y_lst, label='S = 0.8(z - 0.5)^2 + 0.8')

plt.grid()
plt.legend()
plt.show()

# FOR VARIATION OF BETA

# for rho

plt.title('Plasma density')
plt.xlabel('channel length coordinate')
plt.ylabel('rho value')
plt.grid()

algorithm(I_max_main, N_max_main, k_main, gamma_main, 0.1, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_1[N, i] / S_mesh[i]
plt.plot(x_lst, y_lst, label='beta = 0.1')

algorithm(I_max_main, N_max_main, k_main, gamma_main, 0.5, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_1[N, i] / S_mesh[i]
plt.plot(x_lst, y_lst, label='beta = 0.5')

algorithm(I_max_main, N_max_main, k_main, gamma_main, 1.0, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_1[N, i] / S_mesh[i]
plt.plot(x_lst, y_lst, label='beta = 1.0')

algorithm(I_max_main, N_max_main, k_main, gamma_main, 10.0, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_1[N, i] / S_mesh[i]
plt.plot(x_lst, y_lst, label='beta = 10.0')

plt.legend()
plt.show()

# for v

plt.title('Plasma velocity')
plt.xlabel('channel length coordinate')
plt.ylabel('velocity value')
plt.grid()

algorithm(I_max_main, N_max_main, k_main, gamma_main, 0.1, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_2[N, i] / u_cor_1[N, i]
plt.plot(x_lst, y_lst, label='beta = 0.1')

algorithm(I_max_main, N_max_main, k_main, gamma_main, 0.5, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_2[N, i] / u_cor_1[N, i]
plt.plot(x_lst, y_lst, label='beta = 0.5')

algorithm(I_max_main, N_max_main, k_main, gamma_main, 1.0, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_2[N, i] / u_cor_1[N, i]
plt.plot(x_lst, y_lst, label='beta = 1.0')

algorithm(I_max_main, N_max_main, k_main, gamma_main, 10.0, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_2[N, i] / u_cor_1[N, i]
plt.plot(x_lst, y_lst, label='beta = 10.0')

plt.legend()
plt.show()

# for energy

plt.title('Plasma energy')
plt.xlabel('channel length coordinate')
plt.ylabel('energy value')
plt.grid()

algorithm(I_max_main, N_max_main, k_main, gamma_main, 0.1, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_3[N, i] / u_cor_1[N, i]
plt.plot(x_lst, y_lst, label='beta = 0.1')

algorithm(I_max_main, N_max_main, k_main, gamma_main, 0.5, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_3[N, i] / u_cor_1[N, i]
plt.plot(x_lst, y_lst, label='beta = 0.5')

algorithm(I_max_main, N_max_main, k_main, gamma_main, 1.0, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_3[N, i] / u_cor_1[N, i]
plt.plot(x_lst, y_lst, label='beta = 1.0')

algorithm(I_max_main, N_max_main, k_main, gamma_main, 10.0, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_3[N, i] / u_cor_1[N, i]
plt.plot(x_lst, y_lst, label='beta = 10.0')

plt.legend()
plt.show()

# for strength

plt.title('Azimuthal magnetic field')
plt.xlabel('channel length coordinate')
plt.ylabel('strength value')
plt.grid()

algorithm(I_max_main, N_max_main, k_main, gamma_main, 0.1, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_4[N, i] / S_mesh[i]
plt.plot(x_lst, y_lst, label='beta = 0.1')

algorithm(I_max_main, N_max_main, k_main, gamma_main, 0.5, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_4[N, i] / S_mesh[i]
plt.plot(x_lst, y_lst, label='beta = 0.5')

algorithm(I_max_main, N_max_main, k_main, gamma_main, 1.0, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_4[N, i] / S_mesh[i]
plt.plot(x_lst, y_lst, label='beta = 1.0')

algorithm(I_max_main, N_max_main, k_main, gamma_main, 10.0, mu_0_main, rho_0_main, v_0_main, a_main, b_main)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_4[N, i] / S_mesh[i]
plt.plot(x_lst, y_lst, label='beta = 10.0')

plt.legend()
plt.show()

# FOR VARIATION OF GEOMETRY OF CHANNEL

# for rho

plt.title('Plasma density')
plt.xlabel('channel length coordinate')
plt.ylabel('rho value')
plt.grid()

algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_0_main, 3.6, 0.1)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_1[N, i] / S_mesh[i]
plt.plot(x_lst, y_lst, label='S = 3.6(z - 0.5)^2 + 0.1')

algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_0_main, 2.0, 0.5)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_1[N, i] / S_mesh[i]
plt.plot(x_lst, y_lst, label='S = 2(z - 0.5)^2 + 0.5')

algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_0_main, 0.8, 0.8)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_1[N, i] / S_mesh[i]
plt.plot(x_lst, y_lst, label='S = 0.8(z - 0.5)^2 + 0.8')

plt.legend()
plt.show()

# for v

plt.title('Plasma velocity')
plt.xlabel('channel length coordinate')
plt.ylabel('velocity value')
plt.grid()

algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_0_main, 3.6, 0.1)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_2[N, i] / u_cor_1[N, i]
plt.plot(x_lst, y_lst, label='S = 3.6(z - 0.5)^2 + 0.1')

algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_0_main, 2.0, 0.5)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_2[N, i] / u_cor_1[N, i]
plt.plot(x_lst, y_lst, label='S = 2(z - 0.5)^2 + 0.5')

algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_0_main, 0.8, 0.8)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_2[N, i] / u_cor_1[N, i]
plt.plot(x_lst, y_lst, label='S = 0.8(z - 0.5)^2 + 0.8')

plt.legend()
plt.show()

# for energy

plt.title('Plasma energy')
plt.xlabel('channel length coordinate')
plt.ylabel('energy value')
plt.grid()

algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_0_main, 3.6, 0.1)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_3[N, i] / u_cor_1[N, i]
plt.plot(x_lst, y_lst, label='S = 3.6(z - 0.5)^2 + 0.1')

algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_0_main, 2.0, 0.5)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_3[N, i] / u_cor_1[N, i]
plt.plot(x_lst, y_lst, label='S = 2(z - 0.5)^2 + 0.5')

algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_0_main, 0.8, 0.8)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_3[N, i] / u_cor_1[N, i]
plt.plot(x_lst, y_lst, label='S = 0.8(z - 0.5)^2 + 0.8')

plt.legend()
plt.show()

# for strength

plt.title('Azimuthal magnetic field')
plt.xlabel('channel length coordinate')
plt.ylabel('strength value')
plt.grid()

algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_0_main, 3.6, 0.1)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_4[N, i] / S_mesh[i]
plt.plot(x_lst, y_lst, label='S = 3.6(z - 0.5)^2 + 0.1')

algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_0_main, 2.0, 0.5)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_4[N, i] / S_mesh[i]
plt.plot(x_lst, y_lst, label='S = 2(z - 0.5)^2 + 0.5')

algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_0_main, 0.8, 0.8)
for i in range(I_max_main + 1):
    y_lst[i] = u_cor_4[N, i] / S_mesh[i]
plt.plot(x_lst, y_lst, label='S = 0.8(z - 0.5)^2 + 0.8')

plt.legend()
plt.show()
'''