# all quantities are dimensionless

# import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def algorithm(i_max, n_max, k, gamma, beta, mu_0, rho_0, v_z_0, v_phi_0, h_z_s):
    # initial condition

    e_0 = beta / (2 * (gamma - 1))
    p_0 = beta / 2 * rho_0 ** gamma

    # grid steps

    h = 1 / i_max
    tau = k * h / (2 * mu_0 + v_z_0)

    # filling meshes for parameters

    for i in range(i_max + 1):
        rho_mesh[i] = rho_0
        v_phi_mesh[i] = v_phi_0
        v_z_mesh[i] = v_z_0
        e_mesh[i] = e_0
        p_mesh[i] = p_0
        H_phi_mesh[i] = 1 - 0.9 * i * h
        H_z_mesh[i] = h_z_s / S_mesh[i]

    # filling meshes for u

    for i in range(i_max + 1):
        u_td_1[0, i] = rho_mesh[i] * S_mesh[i]
        u_td_2[0, i] = rho_mesh[i] * v_z_mesh[i] * S_mesh[i]
        u_td_3[0, i] = rho_mesh[i] * v_phi_mesh[i] * S_mesh[i]
        u_td_4[0, i] = rho_mesh[i] * e_mesh[i] * S_mesh[i]
        u_td_5[0, i] = H_phi_mesh[i] * S_mesh[i]
        u_td_6[0, i] = H_z_mesh[i] * S_mesh[i]

        u_ad_1[0, i] = rho_mesh[i] * S_mesh[i]
        u_ad_2[0, i] = rho_mesh[i] * v_z_mesh[i] * S_mesh[i]
        u_ad_3[0, i] = rho_mesh[i] * v_phi_mesh[i] * S_mesh[i]
        u_ad_4[0, i] = rho_mesh[i] * e_mesh[i] * S_mesh[i]
        u_ad_5[0, i] = H_phi_mesh[i] * S_mesh[i]
        u_ad_6[0, i] = H_z_mesh[i] * S_mesh[i]

        u_cor_1[0, i] = rho_mesh[i] * S_mesh[i]
        u_cor_2[0, i] = rho_mesh[i] * v_z_mesh[i] * S_mesh[i]
        u_cor_3[0, i] = rho_mesh[i] * v_phi_mesh[i] * S_mesh[i]
        u_cor_4[0, i] = rho_mesh[i] * e_mesh[i] * S_mesh[i]
        u_cor_5[0, i] = H_phi_mesh[i] * S_mesh[i]
        u_cor_6[0, i] = H_z_mesh[i] * S_mesh[i]

    # evolution in time

    for n in range(n_max):

        # left boundary condition for current layer

        rho_mesh[0] = 1
        v_phi_mesh[0] = 0
        v_z_mesh[0] = u_cor_2[n, 1] / (rho_mesh[0] * S_mesh[0])
        e_mesh[0] = beta / (2 * (gamma - 1))
        H_phi_mesh[0] = 1

        # left boundary condition for current layer for transport part

        u_td_1[n + 1, 0] = rho_mesh[0] * S_mesh[0]
        u_td_2[n + 1, 0] = rho_mesh[0] * v_z_mesh[0] * S_mesh[0]
        u_td_3[n + 1, 0] = rho_mesh[0] * v_phi_mesh[0] * S_mesh[0]
        u_td_4[n + 1, 0] = rho_mesh[0] * e_mesh[0] * S_mesh[0]
        u_td_5[n + 1, 0] = H_phi_mesh[0] * S_mesh[0]
        u_td_6[n + 1, 0] = H_z_mesh[0] * S_mesh[0]

        # left boundary condition for current layer for anti-diffusion part

        u_ad_1[n + 1, 0] = rho_mesh[0] * S_mesh[0]
        u_ad_2[n + 1, 0] = rho_mesh[0] * v_z_mesh[0] * S_mesh[0]
        u_ad_3[n + 1, 0] = rho_mesh[0] * v_phi_mesh[0] * S_mesh[0]
        u_ad_4[n + 1, 0] = rho_mesh[0] * e_mesh[0] * S_mesh[0]
        u_ad_5[n + 1, 0] = H_phi_mesh[0] * S_mesh[0]
        u_ad_6[n + 1, 0] = H_z_mesh[0] * S_mesh[0]

        # left boundary condition for current layer for correction part

        u_cor_1[n + 1, 0] = rho_mesh[0] * S_mesh[0]
        u_cor_2[n + 1, 0] = rho_mesh[0] * v_z_mesh[0] * S_mesh[0]
        u_cor_3[n + 1, 0] = rho_mesh[0] * v_phi_mesh[0] * S_mesh[0]
        u_cor_4[n + 1, 0] = rho_mesh[0] * e_mesh[0] * S_mesh[0]
        u_cor_5[n + 1, 0] = H_phi_mesh[0] * S_mesh[0]
        u_cor_6[n + 1, 0] = H_z_mesh[0] * S_mesh[0]

        # constants for current layer

        for i in range(1, i_max):
            mu_left[i] = mu_0 + (abs(v_z_mesh[i - 1]) + abs(v_z_mesh[i])) / 4
            mu_right[i] = mu_0 + (abs(v_z_mesh[i]) + abs(v_z_mesh[i + 1])) / 4
            C_left[i] = tau / h * (v_z_mesh[i - 1] + v_z_mesh[i]) / 4 + tau / h * mu_left[i]
            C_0[i] = 1 - tau / h * (v_z_mesh[i + 1] - v_z_mesh[i - 1]) / 4 - tau / h * (mu_left[i] + mu_right[i])
            C_right[i] = - tau / h * (v_z_mesh[i] + v_z_mesh[i + 1]) / 4 + tau / h * mu_right[i]

            mu_ad_left[i] = mu_left[i] - abs(v_z_mesh[i - 1] + v_z_mesh[i]) / 4
            mu_ad_right[i] = mu_right[i] - abs(v_z_mesh[i] + v_z_mesh[i + 1]) / 4

        # filling central points of grid for transport part

        for i in range(1, i_max):
            u_td_1[n + 1, i] = (C_left[i] * u_td_1[n, i - 1] +
                                C_0[i] * u_td_1[n, i] +
                                C_right[i] * u_td_1[n, i + 1])

            u_td_2[n + 1, i] = (C_left[i] * u_td_2[n, i - 1] +
                                C_0[i] * u_td_2[n, i] +
                                C_right[i] * u_td_2[n, i + 1] -
                                tau / (2 * h) * (p_mesh[i + 1] - p_mesh[i - 1] +
                                                 H_phi_mesh[i + 1] ** 2 / 2 - H_phi_mesh[i - 1] ** 2 / 2) * S_mesh[i])

            u_td_3[n + 1, i] = (C_left[i] * u_td_3[n, i - 1] +
                                C_0[i] * u_td_3[n, i] +
                                C_right[i] * u_td_3[n, i + 1] +
                                tau / (2 * h) * (H_phi_mesh[i + 1] * H_z_mesh[i + 1] * S_mesh[i + 1] -
                                                 H_phi_mesh[i - 1] * H_z_mesh[i - 1] * S_mesh[i - 1]))

            u_td_4[n + 1, i] = (C_left[i] * u_td_4[n, i - 1] +
                                C_0[i] * u_td_4[n, i] +
                                C_right[i] * u_td_4[n, i + 1] -
                                tau / (2 * h) * p_mesh[i] * (v_z_mesh[i + 1] * S_mesh[i + 1] -
                                                             v_z_mesh[i - 1] * S_mesh[i - 1]))

            u_td_5[n + 1, i] = (C_left[i] * u_td_5[n, i - 1] +
                                C_0[i] * u_td_5[n, i] +
                                C_right[i] * u_td_5[n, i + 1] +
                                tau / (2 * h) * (H_z_mesh[i + 1] * v_phi_mesh[i + 1] * S_mesh[i + 1] -
                                                 H_z_mesh[i - 1] * v_phi_mesh[i - 1] * S_mesh[i - 1]))

            u_td_6[n + 1, i] = (C_left[i] * u_td_6[n, i - 1] +
                                C_0[i] * u_td_6[n, i] +
                                C_right[i] * u_td_6[n, i + 1])

        # right boundary condition for current layer for transport part

        u_td_1[n + 1, i_max] = u_td_1[n + 1, i_max - 1]
        u_td_2[n + 1, i_max] = u_td_2[n + 1, i_max - 1]
        u_td_3[n + 1, i_max] = u_td_3[n + 1, i_max - 1]
        u_td_4[n + 1, i_max] = u_td_4[n + 1, i_max - 1]
        u_td_5[n + 1, i_max] = u_td_5[n + 1, i_max - 1]
        u_td_6[n + 1, i_max] = u_td_6[n + 1, i_max - 1]

        # filling central points of grid for anti-diffusion part

        for i in range(1, i_max):
            u_ad_1[n + 1, i] = (u_td_1[n + 1, i] * (1 + tau / h * (mu_ad_right[i] + mu_ad_left[i])) -
                                tau / h * (mu_ad_right[i] * u_td_1[n + 1, i + 1] +
                                           mu_ad_left[i] * u_td_1[n + 1, i - 1]))

            u_ad_2[n + 1, i] = (u_td_2[n + 1, i] * (1 + tau / h * (mu_ad_right[i] + mu_ad_left[i])) -
                                tau / h * (mu_ad_right[i] * u_td_2[n + 1, i + 1] +
                                           mu_ad_left[i] * u_td_2[n + 1, i - 1]))

            u_ad_3[n + 1, i] = (u_td_3[n + 1, i] * (1 + tau / h * (mu_ad_right[i] + mu_ad_left[i])) -
                                tau / h * (mu_ad_right[i] * u_td_3[n + 1, i + 1] +
                                           mu_ad_left[i] * u_td_3[n + 1, i - 1]))

            u_ad_4[n + 1, i] = (u_td_4[n + 1, i] * (1 + tau / h * (mu_ad_right[i] + mu_ad_left[i])) -
                                tau / h * (mu_ad_right[i] * u_td_4[n + 1, i + 1] +
                                           mu_ad_left[i] * u_td_4[n + 1, i - 1]))

            u_ad_5[n + 1, i] = (u_td_5[n + 1, i] * (1 + tau / h * (mu_ad_right[i] + mu_ad_left[i])) -
                                tau / h * (mu_ad_right[i] * u_td_5[n + 1, i + 1] +
                                           mu_ad_left[i] * u_td_5[n + 1, i - 1]))

            u_ad_6[n + 1, i] = (u_td_6[n + 1, i] * (1 + tau / h * (mu_ad_right[i] + mu_ad_left[i])) -
                                tau / h * (mu_ad_right[i] * u_td_6[n + 1, i + 1] +
                                           mu_ad_left[i] * u_td_6[n + 1, i - 1]))

        # right boundary condition for current layer for anti-diffusion part

        u_ad_1[n + 1, i_max] = u_ad_1[n + 1, i_max - 1]
        u_ad_2[n + 1, i_max] = u_ad_2[n + 1, i_max - 1]
        u_ad_3[n + 1, i_max] = u_ad_3[n + 1, i_max - 1]
        u_ad_4[n + 1, i_max] = u_ad_4[n + 1, i_max - 1]
        u_ad_5[n + 1, i_max] = u_ad_5[n + 1, i_max - 1]
        u_ad_6[n + 1, i_max] = u_ad_6[n + 1, i_max - 1]

        # filling meshes for correction part parameters

        for i in range(2, i_max - 1):
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

            s = np.sign(u_td_5[n + 1, i + 1] - u_td_5[n + 1, i])
            f_cor_left_5[i] = s * max(0, min(s * (u_td_5[n + 1, i + 1] - u_td_5[n + 1, i]),
                                             abs(tau / h * mu_ad_left[i] * (u_td_5[n + 1, i] - u_td_5[n + 1, i - 1])),
                                             s * (u_td_5[n + 1, i - 1] - u_td_5[n + 1, i - 2])))
            f_cor_right_5[i] = s * max(0, min(s * (u_td_5[n + 1, i + 2] - u_td_5[n + 1, i + 1]),
                                              abs(tau / h * mu_ad_right[i] * (u_td_5[n + 1, i + 1] - u_td_5[n + 1, i])),
                                              s * (u_td_5[n + 1, i] - u_td_5[n + 1, i - 1])))

            s = np.sign(u_td_6[n + 1, i + 1] - u_td_6[n + 1, i])
            f_cor_left_6[i] = s * max(0, min(s * (u_td_6[n + 1, i + 1] - u_td_6[n + 1, i]),
                                             abs(tau / h * mu_ad_left[i] * (u_td_6[n + 1, i] - u_td_6[n + 1, i - 1])),
                                             s * (u_td_6[n + 1, i - 1] - u_td_6[n + 1, i - 2])))
            f_cor_right_6[i] = s * max(0, min(s * (u_td_6[n + 1, i + 2] - u_td_6[n + 1, i + 1]),
                                              abs(tau / h * mu_ad_right[i] * (u_td_6[n + 1, i + 1] - u_td_6[n + 1, i])),
                                              s * (u_td_6[n + 1, i] - u_td_6[n + 1, i - 1])))

        # filling central points of grid for correction part

        for i in range(2, i_max - 1):
            u_cor_1[n + 1, i] = u_td_1[n + 1, i] + (f_cor_left_1[i] - f_cor_right_1[i])
            u_cor_2[n + 1, i] = u_td_2[n + 1, i] + (f_cor_left_2[i] - f_cor_right_2[i])
            u_cor_3[n + 1, i] = u_td_3[n + 1, i] + (f_cor_left_3[i] - f_cor_right_3[i])
            u_cor_4[n + 1, i] = u_td_4[n + 1, i] + (f_cor_left_4[i] - f_cor_right_4[i])
            u_cor_5[n + 1, i] = u_td_5[n + 1, i] + (f_cor_left_5[i] - f_cor_right_5[i])
            u_cor_6[n + 1, i] = u_td_6[n + 1, i] + (f_cor_left_6[i] - f_cor_right_6[i])

        # boundary condition for current layer for correction part

        u_cor_1[n + 1, 1] = u_ad_1[n + 1, 1]
        u_cor_2[n + 1, 1] = u_ad_2[n + 1, 1]
        u_cor_3[n + 1, 1] = u_ad_3[n + 1, 1]
        u_cor_4[n + 1, 1] = u_ad_4[n + 1, 1]
        u_cor_5[n + 1, 1] = u_ad_5[n + 1, 1]
        u_cor_6[n + 1, 1] = u_ad_6[n + 1, 1]

        u_cor_1[n + 1, i_max - 1] = u_ad_1[n + 1, i_max - 1]
        u_cor_2[n + 1, i_max - 1] = u_ad_2[n + 1, i_max - 1]
        u_cor_3[n + 1, i_max - 1] = u_ad_3[n + 1, i_max - 1]
        u_cor_4[n + 1, i_max - 1] = u_ad_4[n + 1, i_max - 1]
        u_cor_5[n + 1, i_max - 1] = u_ad_5[n + 1, i_max - 1]
        u_cor_6[n + 1, i_max - 1] = u_ad_6[n + 1, i_max - 1]

        u_cor_1[n + 1, i_max] = u_cor_1[n + 1, i_max - 1]
        u_cor_2[n + 1, i_max] = u_cor_2[n + 1, i_max - 1]
        u_cor_3[n + 1, i_max] = u_cor_3[n + 1, i_max - 1]
        u_cor_4[n + 1, i_max] = u_cor_4[n + 1, i_max - 1]
        u_cor_5[n + 1, i_max] = u_cor_5[n + 1, i_max - 1]
        u_cor_6[n + 1, i_max] = u_cor_6[n + 1, i_max - 1]

        # update parameters of problem

        for i in range(i_max + 1):
            rho_mesh[i] = u_cor_1[n + 1, i] / S_mesh[i]
            v_phi_mesh[i] = u_cor_3[n + 1, i] / u_cor_1[n + 1, i]
            v_z_mesh[i] = u_cor_2[n + 1, i] / u_cor_1[n + 1, i]
            e_mesh[i] = u_cor_4[n + 1, i] / u_cor_1[n + 1, i]
            p_mesh[i] = beta / 2 * rho_mesh[i] ** gamma
            H_phi_mesh[i] = u_cor_5[n + 1, i] / S_mesh[i]
            H_z_mesh[i] = h_z_s / S_mesh[i]

        # update tau

        tau = k * h / (2 * mu_0 + np.max(v_z_mesh))


# mesh options

I_max_main = 100
N_max_main = 5000
k_main = 0.5
gamma_main = 1.67
beta_main = 1
mu_0_main = 0.7
rho_0_main = 1
v_z_0_main = 0.1
v_phi_0_main = 0
H_z_S_main = 0.49
h_main = 1 / I_max_main

# u = {rho*S, rho*v_z*S, rho*v_phi*S, rho*energy*S, H_phi*S, H_z*S}

u_td_1 = np.zeros((N_max_main + 1, I_max_main + 1))
u_td_2 = np.zeros((N_max_main + 1, I_max_main + 1))
u_td_3 = np.zeros((N_max_main + 1, I_max_main + 1))
u_td_4 = np.zeros((N_max_main + 1, I_max_main + 1))
u_td_5 = np.zeros((N_max_main + 1, I_max_main + 1))
u_td_6 = np.zeros((N_max_main + 1, I_max_main + 1))

u_ad_1 = np.zeros((N_max_main + 1, I_max_main + 1))
u_ad_2 = np.zeros((N_max_main + 1, I_max_main + 1))
u_ad_3 = np.zeros((N_max_main + 1, I_max_main + 1))
u_ad_4 = np.zeros((N_max_main + 1, I_max_main + 1))
u_ad_5 = np.zeros((N_max_main + 1, I_max_main + 1))
u_ad_6 = np.zeros((N_max_main + 1, I_max_main + 1))

u_cor_1 = np.zeros((N_max_main + 1, I_max_main + 1))
u_cor_2 = np.zeros((N_max_main + 1, I_max_main + 1))
u_cor_3 = np.zeros((N_max_main + 1, I_max_main + 1))
u_cor_4 = np.zeros((N_max_main + 1, I_max_main + 1))
u_cor_5 = np.zeros((N_max_main + 1, I_max_main + 1))
u_cor_6 = np.zeros((N_max_main + 1, I_max_main + 1))

# create meshes for correction part parameters

f_cor_left_1 = np.zeros(I_max_main + 1)
f_cor_right_1 = np.zeros(I_max_main + 1)
f_cor_left_2 = np.zeros(I_max_main + 1)
f_cor_right_2 = np.zeros(I_max_main + 1)
f_cor_left_3 = np.zeros(I_max_main + 1)
f_cor_right_3 = np.zeros(I_max_main + 1)
f_cor_left_4 = np.zeros(I_max_main + 1)
f_cor_right_4 = np.zeros(I_max_main + 1)
f_cor_left_5 = np.zeros(I_max_main + 1)
f_cor_right_5 = np.zeros(I_max_main + 1)
f_cor_left_6 = np.zeros(I_max_main + 1)
f_cor_right_6 = np.zeros(I_max_main + 1)

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
v_phi_mesh = np.zeros(I_max_main + 1)
v_z_mesh = np.zeros(I_max_main + 1)
e_mesh = np.zeros(I_max_main + 1)
p_mesh = np.zeros(I_max_main + 1)
H_phi_mesh = np.zeros(I_max_main + 1)
H_z_mesh = np.zeros(I_max_main + 1)
S_mesh = np.zeros(I_max_main + 1)

# checkout

x_lst = np.zeros(I_max_main + 1)
y_lst = np.zeros(I_max_main + 1)

for i_m in range(1, I_max_main + 1):
    x_lst[i_m] = x_lst[i_m - 1] + 1 / I_max_main

# enter the time layer for which you want to display the values
N = N_max_main

# search critical flux for difference parameters for variation of beta

beta_main = 1.0
H_z_S_main = 0.58
H_z_S_cr_1 = 0.54
H_z_S_cr_2 = 0.6

# S = 1.78(z - 0.5)^2 + 0.555 -> H_z_S_main = 0.58
# S = -3z^3 + 6z^2 - 3z + 1 -> HzS_cr_1 = 0.54
# S = 3(z-1)^3 + 6(z-1)^2 + 3(z-1) + 1 -> HzS_cr_2 = 0.6

plt.title('Plasma density')
plt.xlabel('channel length coordinate')
plt.ylabel('rho value')
plt.grid()

for i_m in range(I_max_main + 1):
    S_mesh[i_m] = 3 * (i_m * h_main - 1) ** 3 + 6 * (i_m * h_main - 1) ** 2 + 3 * (i_m * h_main - 1) + 1
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_cr_2)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_1[N, i_m] / S_mesh[i_m]
plt.plot(x_lst, y_lst)

plt.show()

plt.title('Critical longitudinal magnetic field')
plt.xlabel('channel length coordinate')
plt.ylabel('velocities value')
plt.grid()

for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_2[N, i_m] / u_cor_1[N, i_m]
plt.plot(x_lst, y_lst, label='Longitudinal velocity')
for i_m in range(I_max_main + 1):
    y_lst[i_m] = H_z_mesh[i_m] / np.sqrt(rho_mesh[i_m])
plt.plot(x_lst, y_lst, label='Alfven velocity')
for i_m in range(I_max_main + 1):
    y_lst[i_m] = np.sqrt(1 / 2 * (H_phi_mesh[i_m] ** 2 / rho_mesh[i_m] +
                                  H_z_mesh[i_m] ** 2 / rho_mesh[i_m] +
                                  gamma_main * p_mesh[i_m] / rho_mesh[i_m]) +
                         1 / 2 * np.sqrt((H_phi_mesh[i_m] ** 2 / rho_mesh[i_m] +
                                          H_z_mesh[i_m] ** 2 / rho_mesh[i_m] +
                                          gamma_main * p_mesh[i_m] / rho_mesh[i_m]) ** 2 -
                                         4 * gamma_main * p_mesh[i_m] / rho_mesh[i_m] *
                                         H_z_mesh[i_m] ** 2 / rho_mesh[i_m]))
plt.plot(x_lst, y_lst, label='Fast magneto-sonic velocity')
plt.legend()
plt.show()

# geometry of channel
plt.title('Геометрии канала')
plt.xlabel(r'z')
plt.ylabel(r'r')
for i_m in range(I_max_main + 1):
    y_lst[i_m] = 1.78 * (i_m * h_main - 0.5) ** 2 + 0.555
plt.plot(x_lst, y_lst, label=r'$S = 1.78(z - 0.5)^2 + 0.555$')
for i_m in range(I_max_main + 1):
    y_lst[i_m] = - 3 * (i_m * h_main) ** 3 + 6 * (i_m * h_main) ** 2 - 3 * (i_m * h_main) + 1
plt.plot(x_lst, y_lst, label=r'$S = -3z^3 + 6z^2 - 3z + 1$')
for i_m in range(I_max_main + 1):
    y_lst[i_m] = 3 * (i_m * h_main - 1) ** 3 + 6 * (i_m * h_main - 1) ** 2 + 3 * (i_m * h_main - 1) + 1
plt.plot(x_lst, y_lst, label=r'$S = 3(z-1)^3 + 6(z-1)^2 + 3(z-1) + 1$')
plt.grid()
plt.legend()
plt.show()

# FOR VARIATION OF GEOMETRY OF CHANNEL

# for rho

plt.title(r'Плотность плазмы $(\rho)$')
plt.xlabel(r'z')
plt.ylabel(r'$\rho$')
plt.grid()

for i_m in range(I_max_main + 1):
    y_lst[i_m] = 1.78 * (i_m * h_main - 0.5) ** 2 + 0.555
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_main)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_1[N, i_m] / S_mesh[i_m]
plt.plot(x_lst, y_lst, label=r'$S = 1.78(z - 0.5)^2 + 0.555$')

for i_m in range(I_max_main + 1):
    S_mesh[i_m] = - 3 * (i_m * h_main) ** 3 + 6 * (i_m * h_main) ** 2 - 3 * (i_m * h_main) + 1
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_cr_1)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_1[N, i_m] / S_mesh[i_m]
plt.plot(x_lst, y_lst, label=r'$S = -3z^3 + 6z^2 - 3z + 1$')

for i_m in range(I_max_main + 1):
    S_mesh[i_m] = 3 * (i_m * h_main - 1) ** 3 + 6 * (i_m * h_main - 1) ** 2 + 3 * (i_m * h_main - 1) + 1
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_cr_2)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_1[N, i_m] / S_mesh[i_m]
plt.plot(x_lst, y_lst, label=r'$S = 3(z-1)^3 + 6(z-1)^2 + 3(z-1) + 1$')

plt.legend()
plt.show()

# for v_z

plt.title(r'Продольная скорость $(v_z)$')
plt.xlabel(r'z')
plt.ylabel(r'$v_z$')
plt.grid()

for i_m in range(I_max_main + 1):
    y_lst[i_m] = 1.78 * (i_m * h_main - 0.5) ** 2 + 0.555
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_main)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_2[N, i_m] / u_cor_1[N, i_m]
plt.plot(x_lst, y_lst, label=r'$S = 1.78(z - 0.5)^2 + 0.555$')

for i_m in range(I_max_main + 1):
    S_mesh[i_m] = - 3 * (i_m * h_main) ** 3 + 6 * (i_m * h_main) ** 2 - 3 * (i_m * h_main) + 1
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_cr_1)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_2[N, i_m] / u_cor_1[N, i_m]
plt.plot(x_lst, y_lst, label=r'$S = -3z^3 + 6z^2 - 3z + 1$')

for i_m in range(I_max_main + 1):
    S_mesh[i_m] = 3 * (i_m * h_main - 1) ** 3 + 6 * (i_m * h_main - 1) ** 2 + 3 * (i_m * h_main - 1) + 1
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_cr_2)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_2[N, i_m] / u_cor_1[N, i_m]
plt.plot(x_lst, y_lst, label=r'$S = 3(z-1)^3 + 6(z-1)^2 + 3(z-1) + 1$')

plt.legend()
plt.show()

# for v_phi

plt.title(r'Азимутальная скорость $(v_\phi)$')
plt.xlabel(r'z')
plt.ylabel(r'$v_\phi$')
plt.grid()

for i_m in range(I_max_main + 1):
    y_lst[i_m] = 1.78 * (i_m * h_main - 0.5) ** 2 + 0.555
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_main)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_3[N, i_m] / u_cor_1[N, i_m]
plt.plot(x_lst, y_lst, label=r'$S = 1.78(z - 0.5)^2 + 0.555$')

for i_m in range(I_max_main + 1):
    S_mesh[i_m] = - 3 * (i_m * h_main) ** 3 + 6 * (i_m * h_main) ** 2 - 3 * (i_m * h_main) + 1
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_cr_1)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_3[N, i_m] / u_cor_1[N, i_m]
plt.plot(x_lst, y_lst, label=r'$S = -3z^3 + 6z^2 - 3z + 1$')

for i_m in range(I_max_main + 1):
    S_mesh[i_m] = 3 * (i_m * h_main - 1) ** 3 + 6 * (i_m * h_main - 1) ** 2 + 3 * (i_m * h_main - 1) + 1
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_cr_2)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_3[N, i_m] / u_cor_1[N, i_m]
plt.plot(x_lst, y_lst, label=r'$S = 3(z-1)^3 + 6(z-1)^2 + 3(z-1) + 1$')

plt.legend()
plt.show()
'''
# for energy

plt.title('Plasma energy')
plt.xlabel('channel length coordinate')
plt.ylabel('energy value')
plt.grid()

for i_m in range(I_max_main + 1):
    y_lst[i_m] = 1.78 * (i_m * h_main - 0.5) ** 2 + 0.555
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_main)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_4[N, i_m] / u_cor_1[N, i_m]
plt.plot(x_lst, y_lst, label='S = 1.78(z - 0.5)^2 + 0.555')

for i_m in range(I_max_main + 1):
    S_mesh[i_m] = - 3 * (i_m * h_main) ** 3 + 6 * (i_m * h_main) ** 2 - 3 * (i_m * h_main) + 1
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_cr_1)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_4[N, i_m] / u_cor_1[N, i_m]
plt.plot(x_lst, y_lst, label='S = -3z^3 + 6z^2 - 3z + 1')

for i_m in range(I_max_main + 1):
    S_mesh[i_m] = 3 * (i_m * h_main - 1) ** 3 + 6 * (i_m * h_main - 1) ** 2 + 3 * (i_m * h_main - 1) + 1
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_cr_2)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_4[N, i_m] / u_cor_1[N, i_m]
plt.plot(x_lst, y_lst, label='S = 3(z-1)^3 + 6(z-1)^2 + 3(z-1) + 1')

plt.legend()
plt.show()
'''
# for H_phi

plt.title(r'Азимутальное магнитное поле $H_\phi$')
plt.xlabel(r'z')
plt.ylabel(r'$H_\phi$')
plt.grid()

for i_m in range(I_max_main + 1):
    y_lst[i_m] = 1.78 * (i_m * h_main - 0.5) ** 2 + 0.555
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_main)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_5[N, i_m] / u_cor_1[N, i_m]
plt.plot(x_lst, y_lst, label=r'$S = 1.78(z - 0.5)^2 + 0.555$')

for i_m in range(I_max_main + 1):
    S_mesh[i_m] = - 3 * (i_m * h_main) ** 3 + 6 * (i_m * h_main) ** 2 - 3 * (i_m * h_main) + 1
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_cr_1)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_5[N, i_m] / u_cor_1[N, i_m]
plt.plot(x_lst, y_lst, label=r'$S = -3z^3 + 6z^2 - 3z + 1$')

for i_m in range(I_max_main + 1):
    S_mesh[i_m] = 3 * (i_m * h_main - 1) ** 3 + 6 * (i_m * h_main - 1) ** 2 + 3 * (i_m * h_main - 1) + 1
algorithm(I_max_main, N_max_main, k_main, gamma_main, beta_main, mu_0_main, rho_0_main, v_z_0_main, v_phi_0_main,
          H_z_S_cr_2)
for i_m in range(I_max_main + 1):
    y_lst[i_m] = u_cor_5[N, i_m] / u_cor_1[N, i_m]
plt.plot(x_lst, y_lst, label=r'$S = 3(z-1)^3 + 6(z-1)^2 + 3(z-1) + 1$')

plt.legend()
plt.show()

# for H_z

plt.title('Longitudinal magnetic field')
plt.xlabel('channel length coordinate')
plt.ylabel('H_z value')
plt.grid()

for i_m in range(I_max_main + 1):
    y_lst[i_m] = 1.78 * (i_m * h_main - 0.5) ** 2 + 0.555
for i_m in range(I_max_main + 1):
    y_lst[i_m] = H_z_S_main / S_mesh[i_m]
plt.plot(x_lst, y_lst, label='S = 1.78(z - 0.5)^2 + 0.555')

for i_m in range(I_max_main + 1):
    S_mesh[i_m] = - 3 * (i_m * h_main) ** 3 + 6 * (i_m * h_main) ** 2 - 3 * (i_m * h_main) + 1
for i_m in range(I_max_main + 1):
    y_lst[i_m] = H_z_S_cr_1 / S_mesh[i_m]
plt.plot(x_lst, y_lst, label='S = -3z^3 + 6z^2 - 3z + 1')

for i_m in range(I_max_main + 1):
    S_mesh[i_m] = 3 * (i_m * h_main - 1) ** 3 + 6 * (i_m * h_main - 1) ** 2 + 3 * (i_m * h_main - 1) + 1
for i_m in range(I_max_main + 1):
    y_lst[i_m] = H_z_S_cr_2 / S_mesh[i_m]
plt.plot(x_lst, y_lst, label='S = 3(z-1)^3 + 6(z-1)^2 + 3(z-1) + 1')

plt.legend()
plt.show()
