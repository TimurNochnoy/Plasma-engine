#include<stdio.h>
#include<cmath>
#include<iostream>
#include<fstream>

int main() {

    // task data

    double gamma = 1.67;
    double beta = 1.0;

    // discrete solution area

    double T = 10;
    double t = 0;
    double dt = 0.0001;

    int L_max = 1000;

    double dz = double(1.0 / L_max);

    // creating arrays for axes

    double *z = new double [L_max + 1];

    for (int i = 0; i < L_max + 1; i++) {
        z[i] = i * dz;
    }

    // creating array for area

    double* S = new double[L_max + 1];

    for (int i = 0; i < L_max + 1; i++) {
        S[i] = 2 * pow((z[i] - 0.5), 2) + 0.5;
    }

    // creating arrays for values

    double* rho = new double[L_max + 1];
    double* v_phi = new double[L_max + 1];
    double* v_z = new double[L_max + 1];
    double* p = new double[L_max + 1];
    double* e = new double[L_max + 1];
    double* H_phi = new double[L_max + 1];
    double* H_z = new double[L_max + 1];

    for (int i = 0; i < L_max + 1; i++) {
        rho[i] = 0;
        v_phi[i] = 0;
        v_z[i] = 0;
        p[i] = 0;
        e[i] = 0;
        H_phi[i] = 0;
        H_z[i] = 0.1;
    }

    // creating arrays for components of vector

    double* u_1 = new double[L_max + 1];
    double* u_2 = new double[L_max + 1];
    double* u_3 = new double[L_max + 1];
    double* u_4 = new double[L_max + 1];
    double* u_5 = new double[L_max + 1];
    double* u_6 = new double[L_max + 1];

    double* u_1_next = new double[L_max + 1];
    double* u_2_next = new double[L_max + 1];
    double* u_3_next = new double[L_max + 1];
    double* u_4_next = new double[L_max + 1];
    double* u_5_next = new double[L_max + 1];
    double* u_6_next = new double[L_max + 1];

    for (int i = 0; i < L_max + 1; i++) {
        u_1[i] = 0;
        u_2[i] = 0;
        u_3[i] = 0;
        u_4[i] = 0;
        u_5[i] = 0;
        u_6[i] = 0;

        u_1_next[i] = 0;
        u_2_next[i] = 0;
        u_3_next[i] = 0;
        u_4_next[i] = 0;
        u_5_next[i] = 0;
        u_6_next[i] = 0;
    }

    // create stencil variable

    double stencil_1 = 0;
    double stencil_2 = 0;
    double stencil_3 = 0;

    double* f_right = new double[L_max + 1];
    double* f_left = new double[L_max + 1];

    for (int i = 0; i < L_max + 1; i++) {
        f_right[i] = 0;
        f_left[i] = 0;
    }

    // initial conditions

    for (int i = 0; i < L_max + 1; i++) {
        rho[i] = 1 - 0.9 * z[i];
        v_phi[i] = 0;
        v_z[i] = 0.1;
        p[i] = beta / 2 * pow(rho[i], gamma);
        e[i] = beta / (2 * (gamma - 1));
        H_phi[i] = 1 - 0.9 * z[i];
        H_z[i] = 0 / S[i];
    }

    for (int i = 0; i < L_max + 1; i++) {
        u_1[i] = rho[i] * S[i];
        u_2[i] = rho[i] * v_z[i] * S[i];
        u_3[i] = rho[i] * v_phi[i] * S[i];
        u_4[i] = rho[i] * e[i] * S[i];
        u_5[i] = H_phi[i] * S[i];
        u_6[i] = H_z[i] * S[i];
    }

    // start time

    while (t < T) {

        // Lax-Friedrichs method

        for (int i = 0; i < L_max + 1; i++) {
            if (i != 0 && i != L_max) {
                u_1_next[i] = 0.5 * (u_1[i + 1] + u_1[i - 1]) + 
                              dt * (0 - 
                                    (u_1[i + 1] * v_z[i + 1] - u_1[i - 1] * v_z[i - 1]) / (2 * dz));
                u_2_next[i] = 0.5 * (u_2[i + 1] + u_2[i - 1]) + 
                              dt * (-S[i] * ((p[i + 1] + pow(H_phi[i + 1], 2) / 2.0) - (p[i - 1] + pow(H_phi[i - 1], 2) / 2.0)) / (2 * dz) - 
                                    (u_2[i + 1] * v_z[i + 1] - u_2[i - 1] * v_z[i - 1]) / (2 * dz));
                u_3_next[i] = 0.5 * (u_3[i + 1] + u_3[i - 1]) + 
                              dt * ((H_phi[i + 1] * H_z[i + 1] * S[i + 1] - H_phi[i - 1] * H_z[i - 1] * S[i - 1]) / (2 * dz) - 
                                    (u_3[i + 1] * v_z[i + 1] - u_3[i - 1] * v_z[i - 1]) / (2 * dz));
                u_4_next[i] = 0.5 * (u_4[i + 1] + u_4[i - 1]) + 
                              dt * (-p[i] * (v_z[i + 1] * S[i + 1] - v_z[i - 1] * S[i - 1]) / (2 * dz) - 
                                    (u_4[i + 1] * v_z[i + 1] - u_4[i - 1] * v_z[i - 1]) / (2 * dz));
                u_5_next[i] = 0.5 * (u_5[i + 1] + u_5[i - 1]) + 
                              dt * ((H_z[i + 1] * v_phi[i + 1] * S[i + 1] - H_z[i - 1] * v_phi[i - 1] * S[i - 1]) / (2 * dz) - 
                                    (u_5[i + 1] * v_z[i + 1] - u_5[i - 1] * v_z[i - 1]) / (2 * dz));
                u_6_next[i] = 0.5 * (u_6[i + 1] + u_6[i - 1]);
            } 
            else if (i == 0) {
                u_1_next[i] = 1 * S[i];
                u_2_next[i] = 1 * v_z[i] * S[i];
                u_3_next[i] = 1 * 0 * S[i];
                u_4_next[i] = 1 * (beta / (2 * (gamma - 1))) * S[i];
                u_5_next[i] = 1 * S[i];
                u_6_next[i] = H_z[i] * S[i];
            }
            else if (i == L_max) {
                u_1_next[i] = u_1_next[i - 1];
                u_2_next[i] = u_2_next[i - 1];
                u_3_next[i] = u_3_next[i - 1];
                u_4_next[i] = u_4_next[i - 1];
                u_5_next[i] = u_5_next[i - 1];
                u_6_next[i] = u_6_next[i - 1];
            }
        }

        // calculate stencils

        for (int i = 2; i < L_max - 1; i++) {
            stencil_1 = 1.0 / 3.0 * u_1[i - 2] * v_z[i - 2] - 7.0 / 6.0 * u_1[i - 1] * v_z[i - 1] + 11.0 / 6.0 * u_1[i] * v_z[i];
            stencil_2 = -1.0 / 6.0 * u_1[i - 1] * v_z[i - 1] + 5.0 / 6.0 * u_1[i] * v_z[i] + 1.0 / 3.0 * u_1[i + 1] * v_z[i + 1];
            stencil_3 = 1.0 / 3.0 * u_1[i] * v_z[i] + 5.0 / 6.0 * u_1[i + 1] * v_z[i + 1] - 1.0 / 6.0 * u_1[i + 2] * v_z[i + 2];

            f_right[i] = 1.0 / 10.0 * stencil_1 + 3.0 / 5.0 * stencil_2 + 3.0 / 10.0 * stencil_3;
            f_left[L_max - i] = 1.0 / 10.0 * stencil_1 + 3.0 / 5.0 * stencil_2 + 3.0 / 10.0 * stencil_3;
        }

        // data update

        for (int i = 0; i < L_max + 1; i++) {
            rho[i] = u_1_next[i] / S[i];
            v_phi[i] = u_3_next[i] / u_1_next[i];
            v_z[i] = u_2_next[i] / u_1_next[i];
            p[i] = beta / 2 * pow(rho[i], gamma);
            e[i] = u_4_next[i] / u_1_next[i];
            H_phi[i] = u_5_next[i] / S[i];
            H_z[i] = u_6_next[i] / S[i];

            u_1[i] = rho[i] * S[i];
            u_2[i] = rho[i] * v_z[i] * S[i];
            u_3[i] = rho[i] * v_phi[i] * S[i];
            u_4[i] = rho[i] * e[i] * S[i];
            u_5[i] = H_phi[i] * S[i];
            u_6[i] = H_z[i] * S[i];
        }
        // time step

        t += dt;

        // print data

        printf("t=%lf rho=%lf v_phi=%lf v_z=%lf p=%lf e=%lf H_phi=%lf H_z=%lf\n", 
                t, rho[1], v_phi[1], v_z[1], p[1], e[1], H_phi[1], H_z[1]);
    }

    // output results for Tecplot

    std::ofstream outfile;

    outfile.open("OUTPUT.dat", std::ios::out | std::ios::trunc);
    outfile << "TITLE=\"" << "Graphics" << "\"" << "\n";
    outfile << R"(VARIABLES= "x", "rho", "v_phi", "v_z", "p", "e", "H_phi", "H_z")" << "\n";
    outfile << R"(ZONE T="zone")" << ", I=" << L_max + 1 << ", F=POINT" << "\n";
    for (int i = 0; i < L_max + 1; i++) {
        outfile << z[i] << " " << rho[i] << " " << v_phi[i] << " " << v_z[i] << " " << 
                   p[i] << " " << e[i] << " " << H_phi[i] << " " << H_z[i] << "\n";
    }
    outfile.close();

    // free memory

    delete[] z;
    delete[] S;

    delete[] rho;
    delete[] v_phi;
    delete[] v_z;
    delete[] p;
    delete[] e;
    delete[] H_phi;
    delete[] H_z;

    delete[] u_1;
    delete[] u_2;
    delete[] u_3;
    delete[] u_4;
    delete[] u_5;
    delete[] u_6;

    delete[] u_1_next;
    delete[] u_2_next;
    delete[] u_3_next;
    delete[] u_4_next;
    delete[] u_5_next;
    delete[] u_6_next;

	return 0;
}
