#include <iostream>
#include <fstream>
#include <cmath>

void memory_allocation(double** &array, int rows, int columns) {

	array = new double* [rows + 1];
    for (int i = 0; i < rows + 1; i++) {
        array[i] = new double[columns + 1];
        for (int j = 0; j < columns + 1; j++) {
            array[i][j] = 0;
        }
    }
}

void memory_clearing(double **arr, int rows) {
	for (int i = 0; i < rows + 1; i++) {
		delete [] arr[i];
	}
	delete [] arr;
}

int round(double d) {
	return (int)floor(d + 0.5);
}

double max_for_dt(double *mu_l, double *mu_m, double *R, double dz, double dy, double L, double M) {
	double value;
	double maxim = 0.0;

	for (int l = 1; l < L; l++) {
		for (int m = 1; m < M; m++) {
			value = mu_l[l] / dz + mu_m[m] / (dy * R[l]);

			if (maxim < value) {
				maxim = value;
			}
		}
	}

	return maxim;
}

int sign(double n) {
	int res = 0;

	if (n == 0) {
		res = 0;
	} else {
		res = (int)(n / abs(n));
	}

	return res;
}

double max(double a, double b) {
	double res;

	if (a > b) {
		res = a;
	} else {
		res = b;
	}

	return res;
}

double min(double a, double b, double c) {
	double res;

	if ((a <= b) && (a <= c)) {
		res = a;
	} else if ((b <= a) && (b <= c)) {
		res = b;
	} else if ((c <= a) && (c <= b)) {
		res = c;
	}

	return res;
}

double r1(double z) {
	return 0.55 - pow((z - 0.5), 2);
}

double r2(double z) {
	return 0.8;
}

double der_r1(double z) {
	return -2 * (z - 0.5);
}

double der_r2(double z) {
	return 0;
}


int main() {
	// mesh options

	int L_max = 100;
	int M_max = 100;
	int N_max = 5;

	// U = {rho*r*R, rho*u*r*R, rho*v*r*R, rho*w*r*R, rho*energy*r*R, H_phi*R, H_z*r*R, H_y*r*R}

	double **u_td_1_next;
	double **u_td_2_next;
	double **u_td_3_next;
	double **u_td_4_next;
	double **u_td_5_next;
	double **u_td_6_next;
	double **u_td_7_next;
	double **u_td_8_next;
	memory_allocation(u_td_1_next, L_max, M_max);
	memory_allocation(u_td_2_next, L_max, M_max);
	memory_allocation(u_td_3_next, L_max, M_max);
	memory_allocation(u_td_4_next, L_max, M_max);
	memory_allocation(u_td_5_next, L_max, M_max);
	memory_allocation(u_td_6_next, L_max, M_max);
	memory_allocation(u_td_7_next, L_max, M_max);
	memory_allocation(u_td_8_next, L_max, M_max);

	double **u_td_1;
	double **u_td_2;
	double **u_td_3;
	double **u_td_4;
	double **u_td_5;
	double **u_td_6;
	double **u_td_7;
	double **u_td_8;
	memory_allocation(u_td_1, L_max, M_max);
	memory_allocation(u_td_2, L_max, M_max);
	memory_allocation(u_td_3, L_max, M_max);
	memory_allocation(u_td_4, L_max, M_max);
	memory_allocation(u_td_5, L_max, M_max);
	memory_allocation(u_td_6, L_max, M_max);
	memory_allocation(u_td_7, L_max, M_max);
	memory_allocation(u_td_8, L_max, M_max);

	// constants
	
	double k = 0.5;
	double gamma = 1.67;
	double beta = 1.0;
	double mu_0_l = 0.7;
	double mu_0_m = 0.7;
	double r_0 = (r1(0) + r2(0)) / 2.0;

	// constants on current layer

	double **C_0;
	double **C_l_left;
	double **C_l_right;
	double **C_m_left;
	double **C_m_right;
	memory_allocation(C_0, L_max, M_max);
	memory_allocation(C_l_left, L_max, M_max);
	memory_allocation(C_l_right, L_max, M_max);
	memory_allocation(C_m_left, L_max, M_max);
	memory_allocation(C_m_right, L_max, M_max);

	double *mu_l_left = (double *)malloc((L_max + 1) * sizeof(double));
	double *mu_l_right = (double *)malloc((L_max + 1) * sizeof(double));
	double *mu_m_left = (double *)malloc((M_max + 1) * sizeof(double));
	double *mu_m_right = (double *)malloc((M_max + 1) * sizeof(double));

	// create meshes for parameters of problem

	double **rho;
	double **v_r;
	double **v_phi;
	double **v_z;
	double **e;
	double **p;
	double **P;
	double **H_r;
	double **H_phi;
	double **H_z;
	memory_allocation(rho, L_max, M_max);
	memory_allocation(v_r, L_max, M_max);
	memory_allocation(v_phi, L_max, M_max);
	memory_allocation(v_z, L_max, M_max);
	memory_allocation(e, L_max, M_max);
	memory_allocation(p, L_max, M_max);
	memory_allocation(P, L_max, M_max);
	memory_allocation(H_r, L_max, M_max);
	memory_allocation(H_phi, L_max, M_max);
	memory_allocation(H_z, L_max, M_max);

	double **r;
	double **r_z;
	memory_allocation(r, L_max, M_max);
	memory_allocation(r_z, L_max, M_max);
	double *R = (double *)malloc((L_max + 1) * sizeof(double));

	double **v_y;
	double **H_y;
	memory_allocation(v_y, L_max, M_max);
	memory_allocation(H_y, L_max, M_max);

	// filling zeros

	for (int l = 0; l < L_max + 1; l++) {
		for (int m = 0; m < M_max + 1; m++) {
			u_td_1_next[l][m] = 0;
			u_td_2_next[l][m] = 0;
			u_td_3_next[l][m] = 0;
			u_td_4_next[l][m] = 0;
			u_td_5_next[l][m] = 0;
			u_td_6_next[l][m] = 0;
			u_td_7_next[l][m] = 0;
			u_td_8_next[l][m] = 0;

			u_td_1[l][m] = 0;
			u_td_2[l][m] = 0;
			u_td_3[l][m] = 0;
			u_td_4[l][m] = 0;
			u_td_5[l][m] = 0;
			u_td_6[l][m] = 0;
			u_td_7[l][m] = 0;
			u_td_8[l][m] = 0;

			C_0[l][m] = 0;
			C_l_left[l][m] = 0;
			C_l_right[l][m] = 0;
			C_m_left[l][m] = 0;
			C_m_right[l][m] = 0;

			rho[l][m] = 0;
			v_r[l][m] = 0;
			v_phi[l][m] = 0;
			v_z[l][m] = 0;
			e[l][m] = 0;
			p[l][m] = 0;
			P[l][m] = 0;
			H_r[l][m] = 0;
			H_phi[l][m] = 0;
			H_z[l][m] = 0;

			r[l][m] = 0;
			r_z[l][m] = 0;

			v_y[l][m] = 0;
			H_y[l][m] = 0;
		}
	}

	for (int l = 0; l < L_max + 1; l++) {
		mu_l_left[l] = 0;
		mu_l_right[l] = 0;
		R[l] = 0;
	}

	for (int m = 0; m < M_max + 1; m++) {
		mu_m_left[m] = 0;
		mu_m_right[m] = 0;
	}

	// grid steps

	double dz = 1.0 / L_max;
	double dy = 1.0 / M_max;
	double dt = 1.0 / N_max;

	// initial condition

	for (int l = 0; l < L_max + 1; l++) {
		for (int m = 0; m < M_max + 1; m++) {
			rho[l][m] = 1.0;
			v_r[l][m] = 0;
			v_phi[l][m] = 0;
			v_z[l][m] = 0;
			H_r[l][m] = 0;
			H_phi[l][m] = 0;
			H_z[l][m] = 0;
			
			e[l][m] = beta / (2.0 * (gamma - 1.0));
			p[l][m] = (gamma - 1.0) * rho[l][m] * e[l][m];
			P[l][m] = p[l][m] + 1.0 / 2.0 * (pow(H_z[l][m], 2) + pow(H_r[l][m], 2) + pow(H_phi[l][m], 2));
		}
	}

	for (int l = 0; l < L_max + 1; l++) {
		R[l] = r2(l * dz) - r1(l * dz);

		for (int m = 0; m < M_max + 1; m++) {
			r[l][m] = (1 - m * dy) * r1(l * dz) + m * dy * r2(l * dz);
			r_z[l][m] = (1 - m * dy) * der_r1(l * dz) + m * dy * der_r2(l * dz);

			v_y[l][m] = v_r[l][m] - v_z[l][m] * r_z[l][m];
			H_y[l][m] = H_r[l][m] - H_z[l][m] * r_z[l][m];
		}
	}

	// filling meshes for u

	for (int l = 0; l < L_max + 1; l++) {
		for (int m = 0; m < M_max + 1; m++) {
			u_td_1[l][m] = rho[l][m] * r[l][m] * R[l];
			u_td_2[l][m] = rho[l][m] * v_z[l][m] * r[l][m] * R[l];
			u_td_3[l][m] = rho[l][m] * v_r[l][m] * r[l][m] * R[l];
			u_td_4[l][m] = rho[l][m] * v_phi[l][m] * r[l][m] * R[l];
			u_td_5[l][m] = rho[l][m] * e[l][m] * r[l][m] * R[l];
			u_td_6[l][m] = H_phi[l][m] * R[l];
			u_td_7[l][m] = H_z[l][m] * r[l][m] * R[l];
			u_td_8[l][m] = H_y[l][m] * r[l][m];
		}
	}

	// evolution in time

	for (int n = 0; n < N_max; n++) {
		
		// left boundary condition for current layer

		for (int m = 0; m < M_max + 1; m++) {
			rho[0][m] = 1.0;
			v_phi[0][m] = 0;
			v_y[0][m] = 0;
			H_y[0][m] = 0;
			H_phi[0][m] = r_0 / r[0][m];
			H_z[0][m] = 0;
			e[0][m] = beta / (2.0 * (gamma - 1.0));
		}

		// up an down boundary condition for current layer

		for (int l = 0; l < L_max + 1; l++) {
			v_y[l][0] = 0;
			v_y[l][M_max] = 0;
			H_y[l][0] = 0;
			H_y[l][M_max] = 0;
		}
		
		// left boundary condition for current layer for transport part

		for (int m = 0; m < M_max + 1; m++) {
			u_td_1_next[0][m] = rho[0][m] * r[0][m] * R[0];
			u_td_2_next[0][m] = rho[0][m] * v_z[0][m] * r[0][m] * R[0];
			u_td_3_next[0][m] = rho[0][m] * v_r[0][m] * r[0][m] * R[0];
			u_td_4_next[0][m] = rho[0][m] * v_phi[0][m] * r[0][m] * R[0];
			u_td_5_next[0][m] = rho[0][m] * e[0][m] * r[0][m] * R[0];
			u_td_6_next[0][m] = H_phi[0][m] * R[0];
			u_td_7_next[0][m] = H_z[0][m] * r[0][m] * R[0];
			u_td_8_next[0][m] = H_y[0][m] * r[0][m];
		}

		// up and down boundary condition for current layer for transport part

		for (int l = 0; l < L_max + 1; l++) {
			u_td_1_next[l][0] = rho[l][0] * r[l][0] * R[l];
			u_td_2_next[l][0] = rho[l][0] * v_z[l][0] * r[l][0] * R[l];
			u_td_3_next[l][0] = rho[l][0] * v_r[l][0] * r[l][0] * R[l];
			u_td_4_next[l][0] = rho[l][0] * v_phi[l][0] * r[l][0] * R[l];
			u_td_5_next[l][0] = rho[l][0] * e[l][0] * r[l][0] * R[l];
			u_td_6_next[l][0] = H_phi[l][0] * R[l];
			u_td_7_next[l][0] = H_z[l][0] * r[l][0] * R[l];
			u_td_8_next[l][0] = H_y[l][0] * r[l][0];

			u_td_1_next[l][M_max] = rho[l][M_max] * r[l][M_max] * R[l];
			u_td_2_next[l][M_max] = rho[l][M_max] * v_z[l][M_max] * r[l][M_max] * R[l];
			u_td_3_next[l][M_max] = rho[l][M_max] * v_r[l][M_max] * r[l][M_max] * R[l];
			u_td_4_next[l][M_max] = rho[l][M_max] * v_phi[l][M_max] * r[l][M_max] * R[l];
			u_td_5_next[l][M_max] = rho[l][M_max] * e[l][M_max] * r[l][M_max] * R[l];
			u_td_6_next[l][M_max] = H_phi[l][M_max] * R[l];
			u_td_7_next[l][M_max] = H_z[l][M_max] * r[l][M_max] * R[l];
			u_td_8_next[l][M_max] = H_y[l][M_max] * r[l][M_max];
		}

		// constants for current layer

		for (int l = 1; l < L_max; l++) {
			for (int m = 1; m < M_max; m++) {
				mu_l_left[l] = mu_0_l + (abs(v_y[l - 1][m]) + abs(v_y[l][m])) / 4.0;
				mu_l_right[l] = mu_0_l + (abs(v_y[l][m]) + abs(v_y[l + 1][m])) / 4.0;
				mu_m_left[m] = mu_0_m + (abs(v_y[l][m - 1]) + abs(v_y[l][m])) / 4.0;
				mu_m_right[m] = mu_0_m + (abs(v_y[l][m]) + abs(v_y[l][m + 1])) / 4.0;

				C_0[l][m] = 1 - dt / dz * (v_z[l + 1][m] - v_z[l - 1][m]) / 4.0 -
								dt / (dy * R[l]) * (v_z[l][m + 1] - v_z[l][m - 1]) / 4.0 -
								dt / dz * (mu_l_left[l] + mu_l_right[l]) - 
								dt / (dy * R[l]) * (mu_m_left[m] + mu_m_right[m]);
				C_l_left[l][m] = dt / dz * ((v_z[l - 1][m] + v_z[l][m]) / 4.0 + mu_l_left[l]);
				C_l_right[l][m] = dt / dz * (- (v_z[l][m] + v_z[l + 1][m]) / 4.0 + mu_l_right[l]);
				C_m_left[l][m] = dt / (dy * R[l]) * ((v_y[l][m - 1] + v_y[l][m]) / 4.0 + mu_m_left[m]);
				C_m_right[l][m] = dt / (dy * R[l]) * (- (v_y[l][m] + v_y[l][m + 1]) / 4.0 + mu_m_right[m]);
			}
		}

		// filling central points of grid for transport part

		for (int l = 1; l < L_max; l++) {
			for (int m = 1; m < M_max; m++) {
				u_td_1_next[l][m] = C_0[l][m] * u_td_1[l][m] + 
									C_l_left[l][m] * u_td_1[l - 1][m] + 
									C_l_right[l][m] * u_td_1[l + 1][m] +
									C_m_left[l][m] * u_td_1[l][m - 1] +
									C_m_right[l][m] * u_td_1[l][m + 1];

				u_td_2_next[l][m] = C_0[l][m] * u_td_2[l][m] + 
									C_l_left[l][m] * u_td_2[l - 1][m] + 
									C_l_right[l][m] * u_td_2[l + 1][m] +
									C_m_left[l][m] * u_td_2[l][m - 1] +
									C_m_right[l][m] * u_td_2[l][m + 1] -
									dt / (2.0 * dz) * ((P[l + 1][m] - pow(H_z[l + 1][m], 2)) * r[l + 1][m] * R[l + 1] - 
													   (P[l - 1][m] - pow(H_z[l - 1][m], 2)) * r[l - 1][m] * R[l - 1]) +
									dt / (2.0 * dy) * ((r_z[l][m + 1] * P[l][m + 1] + H_z[l][m + 1] * H_y[l][m + 1]) * r[l][m + 1] -
													   (r_z[l][m - 1] * P[l][m - 1] + H_z[l][m - 1] * H_y[l][m - 1]) * r[l][m - 1]);
				
				u_td_3_next[l][m] = C_0[l][m] * u_td_3[l][m] + 
									C_l_left[l][m] * u_td_3[l - 1][m] + 
									C_l_right[l][m] * u_td_3[l + 1][m] +
									C_m_left[l][m] * u_td_3[l][m - 1] +
									C_m_right[l][m] * u_td_3[l][m + 1] + 
									dt / (2.0 * dz) * (H_z[l + 1][m] * H_r[l + 1][m] * r[l + 1][m] * R[l + 1] - 
													   H_z[l - 1][m] * H_r[l - 1][m] * r[l - 1][m] * R[l - 1]) -
									dt / (2.0 * dy) * ((P[l][m + 1] - H_r[l][m + 1] * H_y[l][m + 1]) * r[l][m + 1] -
													   (P[l][m - 1] - H_r[l][m - 1] * H_y[l][m - 1]) * r[l][m - 1]) +
									dt * (rho[l][m] * pow(v_phi[l][m], 2) + P[l][m] - pow(H_phi[l][m], 2)) * R[l];

				u_td_4_next[l][m] = C_0[l][m] * u_td_4[l][m] + 
									C_l_left[l][m] * u_td_4[l - 1][m] + 
									C_l_right[l][m] * u_td_4[l + 1][m] +
									C_m_left[l][m] * u_td_4[l][m - 1] +
									C_m_right[l][m] * u_td_4[l][m + 1] + 
									dt / (2.0 * dz) * (H_phi[l + 1][m] * H_z[l + 1][m] * r[l + 1][m] * R[l + 1] - 
													   H_phi[l - 1][m] * H_z[l - 1][m] * r[l - 1][m] * R[l - 1]) +
									dt / (2.0 * dy) * (H_z[l][m + 1] * H_y[l][m + 1] * r[l][m + 1] -
													   H_z[l][m - 1] * H_y[l][m - 1] * r[l][m - 1]) +
									dt * (- rho[l][m] * v_r[l][m] * v_phi[l][m] + H_phi[l][m] * H_r[l][m]) * R[l];

				u_td_5_next[l][m] = C_0[l][m] * u_td_5[l][m] + 
									C_l_left[l][m] * u_td_5[l - 1][m] + 
									C_l_right[l][m] * u_td_5[l + 1][m] +
									C_m_left[l][m] * u_td_5[l][m - 1] +
									C_m_right[l][m] * u_td_5[l][m + 1] -
									dt * p[l][m] * (1.0 / (2.0 * dz) * (v_z[l + 1][m] * r[l + 1][m] * R[l + 1] -
																		v_z[l - 1][m] * r[l - 1][m] * R[l - 1]) + 
													1.0 / (2.0 * dy) * (v_y[l][m + 1] * r[l][m + 1] -
																		v_y[l][m - 1] * r[l][m - 1]));

				u_td_6_next[l][m] = C_0[l][m] * u_td_6[l][m] + 
									C_l_left[l][m] * u_td_6[l - 1][m] + 
									C_l_right[l][m] * u_td_6[l + 1][m] +
									C_m_left[l][m] * u_td_6[l][m - 1] +
									C_m_right[l][m] * u_td_6[l][m + 1] +
									dt / (2.0 * dz) * (H_phi[l + 1][m] * v_phi[l + 1][m] * R[l + 1] - 
													   H_phi[l - 1][m] * v_phi[l - 1][m] * R[l - 1]) + 
									dt / (2.0 * dy) * (H_y[l][m + 1] * v_phi[l][m + 1] -
													   H_y[l][m - 1] * v_phi[l][m - 1]);
				
				u_td_7_next[l][m] = C_0[l][m] * u_td_7[l][m] + 
									C_l_left[l][m] * u_td_7[l - 1][m] + 
									C_l_right[l][m] * u_td_7[l + 1][m] +
									C_m_left[l][m] * u_td_7[l][m - 1] +
									C_m_right[l][m] * u_td_7[l][m + 1] +
									dt / (2.0 * dy) * (H_y[l][m + 1] * v_z[l][m + 1] * r[l][m + 1] -
													   H_y[l][m - 1] * v_z[l][m - 1] * r[l][m - 1]);

				u_td_8_next[l][m] = C_0[l][m] * u_td_8[l][m] + 
									C_l_left[l][m] * u_td_8[l - 1][m] + 
									C_l_right[l][m] * u_td_8[l + 1][m] +
									C_m_left[l][m] * u_td_8[l][m - 1] +
									C_m_right[l][m] * u_td_8[l][m + 1] + 
									dt / (2.0 * dz) * (H_z[l + 1][m] * v_y[l + 1][m] * r[l + 1][m] - 
													   H_z[l - 1][m] * v_y[l - 1][m] * r[l - 1][m]);
			}
		}

		// right boundary condition for current layer for transport part

		for (int m = 0; m < M_max; m++) {
			u_td_1_next[L_max][m] = u_td_1_next[L_max - 1][m];
			u_td_2_next[L_max][m] = u_td_2_next[L_max - 1][m];
			u_td_3_next[L_max][m] = u_td_3_next[L_max - 1][m];
			u_td_4_next[L_max][m] = u_td_4_next[L_max - 1][m];
			u_td_5_next[L_max][m] = u_td_5_next[L_max - 1][m];
			u_td_6_next[L_max][m] = u_td_6_next[L_max - 1][m];
			u_td_7_next[L_max][m] = u_td_7_next[L_max - 1][m];
			u_td_8_next[L_max][m] = u_td_8_next[L_max - 1][m];
		}
		
		// update parameters of problem

		for (int l = 0; l < L_max + 1; l++) {
			for (int m = 0; m < M_max + 1; m++) {
				rho[l][m] = u_td_1_next[l][m] / (r[l][m] * R[l]);
				v_r[l][m] = u_td_3_next[l][m] / u_td_1_next[l][m];
				v_phi[l][m] = u_td_4_next[l][m] / u_td_1_next[l][m];
				v_z[l][m] = u_td_2_next[l][m] / u_td_1_next[l][m];
				v_y[l][m] = v_r[l][m] - v_z[l][m] * r_z[l][m];

				H_y[l][m] = u_td_8_next[l][m] / r[l][m];
				H_phi[l][m] = u_td_6_next[l][m] / R[l];
				H_z[l][m] = u_td_7_next[l][m] / (r[l][m] * R[l]);
				H_r[l][m] = H_y[l][m] + H_z[l][m] * r_z[l][m];

				e[l][m] = u_td_5_next[l][m] / u_td_1_next[l][m];
				p[l][m] = (gamma - 1) * rho[l][m] * e[l][m];
				P[l][m] = p[l][m] + 1.0 / 2.0 * (pow(H_z[l][m], 2) + pow(H_r[l][m], 2) + pow(H_phi[l][m], 2));
			}
		}

		// update u_td

		for (int l = 0; l < L_max + 1; l++) {
			for (int m = 0; m < M_max + 1; m++) {
				u_td_1[l][m] = u_td_1_next[l][m];
				u_td_2[l][m] = u_td_2_next[l][m];
				u_td_3[l][m] = u_td_3_next[l][m];
				u_td_4[l][m] = u_td_4_next[l][m];
				u_td_5[l][m] = u_td_5_next[l][m];
				u_td_6[l][m] = u_td_6_next[l][m];
				u_td_7[l][m] = u_td_7_next[l][m];
				u_td_8[l][m] = u_td_8_next[l][m];
			}
		}
		
		// update tau
		
		dt = k / (2.0 * max_for_dt(mu_l_right, mu_m_right, R, dz, dy, L_max, M_max));
		
	}

	// checkout

	double *z_lst = (double *)malloc((L_max + 1) * sizeof(double));
	double *y_lst = (double *)malloc((M_max + 1) * sizeof(double));

	double **rho_lst;
	double **v_r_lst;
	double **v_phi_lst;
	double **v_z_lst;
	double **e_lst;
	double **p_lst;
	double **P_lst;
	double **H_r_lst;
	double **H_phi_lst;
	double **H_z_lst;
	memory_allocation(rho_lst, L_max, M_max);
	memory_allocation(v_r_lst, L_max, M_max);
	memory_allocation(v_phi_lst, L_max, M_max);
	memory_allocation(v_z_lst, L_max, M_max);
	memory_allocation(e_lst, L_max, M_max);
	memory_allocation(p_lst, L_max, M_max);
	memory_allocation(P_lst, L_max, M_max);
	memory_allocation(H_r_lst, L_max, M_max);
	memory_allocation(H_phi_lst, L_max, M_max);
	memory_allocation(H_z_lst, L_max, M_max);

	z_lst[0] = 0;
	for (int l = 1; l < L_max + 1; l++) {
		z_lst[l] = z_lst[l - 1] + dz;
	}
	y_lst[0] = 0;
	for (int m = 1; m < M_max + 1; m++) {
		y_lst[m] = y_lst[m - 1] + dy;
	}

	for (int l = 0; l < L_max + 1; l++) {
		for (int m = 0; m < M_max + 1; m++) {
			rho_lst[l][m] = u_td_1_next[l][m] / (r[l][m] * R[l]);
			v_r_lst[l][m] = u_td_3_next[l][m] / u_td_1_next[l][m];
			v_phi_lst[l][m] = u_td_4_next[l][m] / u_td_1_next[l][m];
			v_z_lst[l][m] = u_td_2_next[l][m] / u_td_1_next[l][m];

			H_phi_lst[l][m] = u_td_6_next[l][m] / R[l];
			H_z_lst[l][m] = u_td_7_next[l][m] / (r[l][m] * R[l]);
			H_r_lst[l][m] = H_y[l][m] + H_z[l][m] * r_z[l][m];

			e_lst[l][m] = u_td_5_next[l][m] / u_td_1_next[l][m];
			p_lst[l][m] = (gamma - 1) * rho[l][m] * e[l][m];
			P_lst[l][m] = p[l][m] + 1.0 / 2.0 * (pow(H_z[l][m], 2) + pow(H_r[l][m], 2) + pow(H_phi[l][m], 2));
		}
	}

	// output results in file

	std :: ofstream f("data_MGD_CPP.txt");

	for (int l = 0; l < L_max + 1; l++) {
		for (int m = 0; m < M_max + 1; m++) {
			f << z_lst[l] << " " << y_lst[m] << " " << rho_lst[l][m] << " " << v_r_lst[l][m] << " " <<
				 v_phi_lst[l][m] << " " << v_z_lst[l][m] << " " << H_phi_lst[l][m] << " " << H_z_lst[l][m] << " " <<
				 H_r_lst[l][m] << " " << e_lst[l][m] << " " << p_lst[l][m] << " " << P_lst[l][m] << " " << "\n";
		}
	}

	f.close();

	// free memory
/*
	memory_clearing(u_td_1, L_max);
	memory_clearing(u_td_2, L_max);
	memory_clearing(u_td_3, L_max);
	memory_clearing(u_td_4, L_max);
	memory_clearing(u_td_5, L_max);
	memory_clearing(u_td_6, L_max);
	memory_clearing(u_td_7, L_max);
	memory_clearing(u_td_8, L_max);

	memory_clearing(u_td_1_next, L_max);
	memory_clearing(u_td_2_next, L_max);
	memory_clearing(u_td_3_next, L_max);
	memory_clearing(u_td_4_next, L_max);
	memory_clearing(u_td_5_next, L_max);
	memory_clearing(u_td_6_next, L_max);
	memory_clearing(u_td_7_next, L_max);
	memory_clearing(u_td_8_next, L_max);

	memory_clearing(C_0, L_max);
	memory_clearing(C_l_left, L_max);
	memory_clearing(C_l_right, L_max);
	memory_clearing(C_m_left, L_max);
	memory_clearing(C_m_right, L_max);

	delete [] mu_l_left;
	delete [] mu_l_right;
	delete [] mu_m_left;
	delete [] mu_m_right;

	memory_clearing(rho, L_max);
	memory_clearing(v_r, L_max);
	memory_clearing(v_phi, L_max);
	memory_clearing(e, L_max);
	memory_clearing(p, L_max);
	memory_clearing(P, L_max);
	memory_clearing(H_r, L_max);
	memory_clearing(H_r, L_max);
	memory_clearing(H_z, L_max);

	memory_clearing(r, L_max);
	memory_clearing(r_z, L_max);
	delete [] R;

	memory_clearing(v_y, L_max);
	memory_clearing(H_y, L_max);

	delete [] z_lst;
	delete [] y_lst;

	memory_clearing(rho_lst, L_max);
	memory_clearing(v_r_lst, L_max);
	memory_clearing(v_phi_lst, L_max);
	memory_clearing(v_z_lst, L_max);
	memory_clearing(e_lst, L_max);
	memory_clearing(p_lst, L_max);
	memory_clearing(P_lst, L_max);
	memory_clearing(H_r_lst, L_max);
	memory_clearing(H_phi_lst, L_max);
	memory_clearing(H_z_lst, L_max);
*/
	return 0;
}
