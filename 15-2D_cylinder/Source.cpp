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

double max_for_dt(double **mu_l, double **mu_m, double dz, double* dr, double L, double M) {
	double value;
	double maxim = 0.0;

	for (int l = 1; l < L; l++) {
		for (int m = 1; m < M; m++) {
			value = mu_l[l][m] / dz + mu_m[l][m] / (dr[l]);

			if (maxim < value) {
				maxim = value;
			}
		}
	}

	return maxim;
}

double max_array(double **array, double L, double M) {
	double maxim = 0.0;

	for (int l = 0; l < L + 1; l++) {
		for (int m = 0; m < M + 1; m++) {
			if (maxim < std::fabs(array[l][m])) {
				maxim = std::fabs(array[l][m]);
			}
		}
	}

	return maxim;
}

double min_array(double *array, double L) {
	double minim = array[0];

	for (int l = 1; l < L; l++) {
		if (minim > array[l]) {
			minim = array[l];
		}
	}

	return minim;
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

	int L_max = 80;
	int M_max = 40;
	int N_max = 15000;
	
	// U = {rho*r*R, rho*u*r*R, rho*v*r*R, rho*w*r*R, rho*energy*r*R, H_phi*R, H_z*r*R, H_y*r}

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

	double **mu_l_left;
	double **mu_l_right;
	double **mu_m_left;
	double **mu_m_right;
	memory_allocation(mu_l_left, L_max, M_max);
	memory_allocation(mu_l_right, L_max, M_max);
	memory_allocation(mu_m_left, L_max, M_max);
	memory_allocation(mu_m_right, L_max, M_max);

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
	double *dr = (double *)malloc((L_max + 1) * sizeof(double));

	double **v_y;
	double **H_y;
	memory_allocation(v_y, L_max, M_max);
	memory_allocation(H_y, L_max, M_max);

	// filling zeros

	for (int l = 0; l < L_max + 1; l++) {
		R[l] = 0;
		dr[l] = 0;
	}

	// grid steps

	double dz = 1.0 / L_max;
	double dy = 1.0 / M_max;
	double dt = 1.0 / N_max;

	// initial condition

	for (int l = 0; l < L_max + 1; l++) {
		R[l] = r2(l * dz) - r1(l * dz);
		dr[l] = R[l] / M_max;
		
		for (int m = 0; m < M_max + 1; m++) {
			r[l][m] = (1 - m * dy) * r1(l * dz) + m * dy * r2(l * dz);
			r_z[l][m] = (1 - m * dy) * der_r1(l * dz) + m * dy * der_r2(l * dz);
		}
	}

	//printf("%lf\n", r[50][50]);

	//for (int m = 0; m < M_max + 1; m++) {
	//	H_z[L_max][m] = 0;
	//}
	//for (int l = 0; l < L_max + 1; l++) {
	//	H_z[l][M_max] = H_z[L_max][M_max] * r[L_max][M_max] * dr[L_max] / r[l][M_max] / dr[l];
	//}
	//for (int l = 0; l < L_max; l++) {
	//	for (int m = M_max; m >= 0; m--) {
	//		H_z[l][m] = 2 * H_z[L_max][m] * ((r[L_max][m] + r[L_max][m + 1]) * dr[L_max] / (r[l][m] + r[l][m + 1]) / dr[l]) - H_z[l][m + 1];
	//	}
	//}


	for (int l = 0; l < L_max + 1; l++) {
		for (int m = 0; m < M_max + 1; m++) {
			rho[l][m] = 1.0;
			v_r[l][m] = 0.1;
			v_phi[l][m] = 0;
			v_z[l][m] = 0.1;
			H_r[l][m] = 0;
			H_phi[l][m] = (1 - 0.9 * l * dz) * r_0 / r[l][m];
			//printf("%lf\n", H_phi[l][m]);
			H_z[l][m] = 0;
			
			e[l][m] = beta / (2.0 * (gamma - 1.0));
			p[l][m] = beta / 2.0;
			P[l][m] = p[l][m] + 1.0 / 2.0 * (pow(H_z[l][m], 2) + pow(H_r[l][m], 2) + pow(H_phi[l][m], 2));

			v_y[l][m] = v_r[l][m] - v_z[l][m] * r_z[l][m];
			H_y[l][m] = H_r[l][m] - H_z[l][m] * r_z[l][m];
		}
	}

	// filling meshes for u

	for (int l = 0; l < L_max + 1; l++) {
		for (int m = 0; m < M_max + 1; m++) {
			u_td_1[l][m] = rho[l][m] * r[l][m];
			u_td_2[l][m] = rho[l][m] * v_z[l][m] * r[l][m];
			u_td_3[l][m] = rho[l][m] * v_r[l][m] * r[l][m];
			u_td_4[l][m] = rho[l][m] * v_phi[l][m] * r[l][m];
			u_td_5[l][m] = rho[l][m] * e[l][m] * r[l][m];
			u_td_6[l][m] = H_phi[l][m];
			u_td_7[l][m] = H_z[l][m] * r[l][m];
			u_td_8[l][m] = H_r[l][m] * r[l][m];
		}
	}

	// evolution in time

	for (int n = 0; n < N_max; n++) {
		
		// left boundary condition for current layer

		for (int m = 0; m < M_max + 1; m++) {
			rho[0][m] = 1.0;
			v_phi[0][m] = 0;
			//v_y[0][m] = 0;
			v_r[0][m] = v_z[0][m] * r_z[0][m];
			H_z[0][m] = 0;
			//H_y[0][m] = 0;
			H_r[0][m] = H_z[0][m] * r_z[0][m];
			H_phi[0][m] = r_0 / r[0][m];
			//printf("%lf\n", H_phi[0][m]);
			e[0][m] = beta / (2.0 * (gamma - 1.0)) * pow(rho[0][m], gamma - 1.0);
		}

		// up an down boundary condition for current layer

		for (int l = 0; l < L_max + 1; l++) {
			/*v_y[l][0] = 0;
			v_y[l][M_max] = 0;
			H_y[l][0] = 0;
			H_y[l][M_max] = 0;*/
			v_r[l][0] = v_z[l][0] * r_z[l][0];
			v_r[l][M_max] = v_z[l][M_max] * r_z[l][M_max];
			H_r[l][0] = H_z[l][0] * r_z[l][0];
			H_r[l][M_max] = H_z[l][M_max] * r_z[l][M_max];
		}
		
		// left boundary condition for current layer for transport part

		for (int m = 0; m < M_max + 1; m++) {
			u_td_1_next[0][m] = rho[0][m] * r[0][m];
			u_td_2_next[0][m] = rho[0][m] * v_z[0][m] * r[0][m];
			u_td_3_next[0][m] = rho[0][m] * v_r[0][m] * r[0][m];
			u_td_4_next[0][m] = rho[0][m] * v_phi[0][m] * r[0][m];
			u_td_5_next[0][m] = rho[0][m] * e[0][m] * r[0][m];
			u_td_6_next[0][m] = H_phi[0][m];
			u_td_7_next[0][m] = H_z[0][m] * r[0][m];
			u_td_8_next[0][m] = H_r[0][m] * r[0][m]; // y -> r
		}

		// up and down boundary condition for current layer for transport part

		for (int l = 0; l < L_max + 1; l++) {
			u_td_1_next[l][0] = rho[l][0] * r[l][0];
			u_td_2_next[l][0] = rho[l][0] * v_z[l][0] * r[l][0];
			u_td_3_next[l][0] = rho[l][0] * v_r[l][0] * r[l][0];
			u_td_4_next[l][0] = rho[l][0] * v_phi[l][0] * r[l][0];
			u_td_5_next[l][0] = rho[l][0] * e[l][0] * r[l][0];
			u_td_6_next[l][0] = H_phi[l][0];
			u_td_7_next[l][0] = H_z[l][0] * r[l][0];
			u_td_8_next[l][0] = H_r[l][0] * r[l][0]; // y -> r

			u_td_1_next[l][M_max] = rho[l][M_max] * r[l][M_max];
			u_td_2_next[l][M_max] = rho[l][M_max] * v_z[l][M_max] * r[l][M_max];
			u_td_3_next[l][M_max] = rho[l][M_max] * v_r[l][M_max] * r[l][M_max];
			u_td_4_next[l][M_max] = rho[l][M_max] * v_phi[l][M_max] * r[l][M_max];
			u_td_5_next[l][M_max] = rho[l][M_max] * e[l][M_max] * r[l][M_max];
			u_td_6_next[l][M_max] = H_phi[l][M_max];
			u_td_7_next[l][M_max] = H_z[l][M_max] * r[l][M_max];
			u_td_8_next[l][M_max] = H_r[l][M_max] * r[l][M_max]; // y -> r
		}

		// constants for current layer
		
		for (int l = 1; l < L_max; l++) {
			for (int m = 1; m < M_max; m++) {
				mu_l_left[l][m] = mu_0_l + (fabs(v_z[l - 1][m]) + fabs(v_z[l][m])) / 4.0;
				mu_l_right[l][m] = mu_0_l + (fabs(v_z[l][m]) + fabs(v_z[l + 1][m])) / 4.0;
				mu_m_left[l][m] = mu_0_m + (fabs(v_y[l][m - 1]) + fabs(v_y[l][m])) / 4.0; // y -> r
				mu_m_right[l][m] = mu_0_m + (fabs(v_y[l][m]) + fabs(v_y[l][m + 1])) / 4.0; // y -> r
			}
		}

		// update tau
		
		double* tmp;
		tmp = new double [L_max + 1];

		for (int l = 0; l < L_max + 1; l++) {
			tmp[l] = dr[l] / (dr[l] * max_array(mu_l_left, L_max, M_max) + dz * max_array(mu_m_left, L_max, M_max));
		}
		dt = k * dz * min_array(tmp, L_max + 1) / 2.0;

		delete [] tmp;
		//printf("%lf\n", dt);
		// found zero value in dt

		// constants for current layer

		for (int l = 1; l < L_max; l++) {
			for (int m = 1; m < M_max; m++) {
				C_0[l][m] = 1 - dt / dz * (v_z[l + 1][m] - v_z[l - 1][m]) / 4.0 - dt / dr[l] * (v_y[l][m + 1] - v_y[l][m - 1]) / 4.0 // y -> r
							  - dt / dz * (mu_l_left[l][m] + mu_l_right[l][m]) - dt / dr[l] * (mu_m_left[l][m] + mu_m_right[l][m]);
				C_l_left[l][m] = dt / dz * ((v_z[l - 1][m] + v_z[l][m]) / 4.0 + mu_l_left[l][m]);
				C_l_right[l][m] = dt / dz * (- (v_z[l][m] + v_z[l + 1][m]) / 4.0 + mu_l_right[l][m]);
				C_m_left[l][m] = dt / dr[l] * ((v_r[l][m - 1] + v_r[l][m]) / 4.0 + mu_m_left[l][m]); // y -> r
				C_m_right[l][m] = dt / dr[l] * (- (v_r[l][m] + v_r[l][m + 1]) / 4.0 + mu_m_right[l][m]); // y -> r
			}
		}


		printf("n=%d   rho[40][20]=%lf   v_z[40][20]=%lf   v_phi[40][20]=%lf\n", n, rho[40][20], v_z[40][20], v_phi[40][20]);

		//printf("mu_l_l=%lf   mu_l_r=%lf   mu_m_l=%lf   mu_m_r=%lf   C0=%lf\n", mu_l_left[50][50], mu_l_right[50][50],
																			   //mu_m_left[50][50], mu_m_right[50][50], C_0[50][50]);
		//printf("C_l_l=%lf   C_l_r=%lf   C_m_l=%lf   C_m_r=%lf   C0=%lf\n", C_l_left[50][50], C_l_right[50][50],  C_m_left[50][50], C_m_right[50][50], C_0[50][50]);
		//printf("v_z=%lf   v_y=%lf   C0=%lf\n", v_z[50][50], v_y[50][50], C_0[50][50]);
		//printf("v_z_l_l=%lf   v_z_l_r=%lf   v_z_m_r=%lf   v_z_m_r=%lf\n", v_z[49][50], v_y[51][50], v_z[50][49], v_z[50][51]);
		//printf("R=%lf\n", R[50]);
		//printf("dt=%lf   dz=%lf   dy=%lf\n", dt, dz, dy);

		//printf("\n\n");




		//printf("v_z[0][20]=%lf   v_z[40][0]=%lf   v_z[40][40]=%lf\n", v_z[0][20], v_z[40][0], v_z[40][40]);
		//printf("rho[0][20]=%lf   rho[40][0]=%lf   rho[40][40]=%lf\n", rho[0][20], rho[40][0], rho[40][40]);





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
									dt / (2.0 * dz) * ((P[l + 1][m] - pow(H_z[l + 1][m], 2)) * r[l + 1][m] - 
													   (P[l - 1][m] - pow(H_z[l - 1][m], 2)) * r[l - 1][m]) +
									dt / (2 * dr[l]) * (H_z[l][m + 1] * H_r[l][m + 1] * r[l][m + 1] - 
														H_z[l][m - 1] * H_r[l][m - 1] * r[l][m - 1]);

				//printf("v_z[0][20]=%lf   v_z[40][0]=%lf   v_z[40][40]=%lf\n", v_z[0][20], v_z[40][0], v_z[40][40]);

				//printf("u_td_2_next[0][20]=%lf   u_td_2_next[40][0]=%lf   u_td_2_next[40][40]=%lf\n", u_td_2_next[0][20], u_td_2_next[40][0], u_td_2_next[40][40]);
				/*printf("u_td_2_next[0][20]=%lf\n", dt / (2.0 * dz) * ((P[l + 1][m] - pow(H_z[l + 1][m], 2)) * r[l + 1][m] - 
													   (P[l - 1][m] - pow(H_z[l - 1][m], 2)) * r[l - 1][m]) +
									dt / (2 * dr[l]) * (H_z[l][m + 1] * H_r[l][m + 1] * r[l][m + 1] - 
														H_z[l][m - 1] * H_r[l][m - 1] * r[l][m - 1]));*/


				u_td_3_next[l][m] = C_0[l][m] * u_td_3[l][m] + 
									C_l_left[l][m] * u_td_3[l - 1][m] + 
									C_l_right[l][m] * u_td_3[l + 1][m] +
									C_m_left[l][m] * u_td_3[l][m - 1] +
									C_m_right[l][m] * u_td_3[l][m + 1] + 
									dt / (2.0 * dz) * (H_z[l + 1][m] * H_r[l + 1][m] * r[l + 1][m] - 
													   H_z[l - 1][m] * H_r[l - 1][m] * r[l - 1][m]) -
									dt / (2.0 * dr[l]) * ((P[l][m + 1] - pow(H_r[l][m + 1], 2)) * r[l][m + 1] - // y -> r
													      (P[l][m - 1] - pow(H_r[l][m - 1], 2)) * r[l][m - 1]) + // y -> r
									dt * (rho[l][m] * pow(v_phi[l][m], 2) + P[l][m] - pow(H_phi[l][m], 2))/* * R[l]*/;

				u_td_4_next[l][m] = C_0[l][m] * u_td_4[l][m] + 
									C_l_left[l][m] * u_td_4[l - 1][m] + 
									C_l_right[l][m] * u_td_4[l + 1][m] +
									C_m_left[l][m] * u_td_4[l][m - 1] +
									C_m_right[l][m] * u_td_4[l][m + 1] + 
									dt / (2.0 * dz) * (H_phi[l + 1][m] * H_z[l + 1][m] * r[l + 1][m] - 
													   H_phi[l - 1][m] * H_z[l - 1][m] * r[l - 1][m]) +
									dt / (2.0 * dr[l]) * (H_phi[l][m + 1] * H_r[l][m + 1] * r[l][m + 1] - // y -> r
													      H_phi[l][m - 1] * H_r[l][m - 1] * r[l][m - 1]) + // y -> r
									dt * (- rho[l][m] * r[l][m] * v_phi[l][m] + H_phi[l][m] * H_r[l][m])/* * R[l]*/;

				u_td_5_next[l][m] = C_0[l][m] * u_td_5[l][m] + 
									C_l_left[l][m] * u_td_5[l - 1][m] + 
									C_l_right[l][m] * u_td_5[l + 1][m] +
									C_m_left[l][m] * u_td_5[l][m - 1] +
									C_m_right[l][m] * u_td_5[l][m + 1] -
									dt * p[l][m] * (1.0 / (2.0 * dz) * (v_z[l + 1][m] * r[l + 1][m] -
																		v_z[l - 1][m] * r[l - 1][m]) + 
													1.0 / (2.0 * dr[l]) * (v_r[l][m + 1] * r[l][m + 1] - // y -> r
																		   v_r[l][m - 1] * r[l][m - 1])); // y -> r

				u_td_6_next[l][m] = C_0[l][m] * u_td_6[l][m] + 
									C_l_left[l][m] * u_td_6[l - 1][m] + 
									C_l_right[l][m] * u_td_6[l + 1][m] +
									C_m_left[l][m] * u_td_6[l][m - 1] +
									C_m_right[l][m] * u_td_6[l][m + 1] +
									dt / (2.0 * dz) * (H_z[l + 1][m] * v_phi[l + 1][m] - 
													   H_z[l - 1][m] * v_phi[l - 1][m]) + 
									dt / (2.0 * dr[l]) * (H_r[l][m + 1] * v_phi[l][m + 1] - // y -> r
													      H_r[l][m - 1] * v_phi[l][m - 1]); // y -> r
				
				u_td_7_next[l][m] = C_0[l][m] * u_td_7[l][m] + 
									C_l_left[l][m] * u_td_7[l - 1][m] + 
									C_l_right[l][m] * u_td_7[l + 1][m] +
									C_m_left[l][m] * u_td_7[l][m - 1] +
									C_m_right[l][m] * u_td_7[l][m + 1] +
									dt / (2.0 * dr[l]) * (H_r[l][m + 1] * v_z[l][m + 1] * r[l][m + 1] - // y -> r
													      H_r[l][m - 1] * v_z[l][m - 1] * r[l][m - 1]); // y -> r

				u_td_8_next[l][m] = C_0[l][m] * u_td_8[l][m] + 
									C_l_left[l][m] * u_td_8[l - 1][m] + 
									C_l_right[l][m] * u_td_8[l + 1][m] +
									C_m_left[l][m] * u_td_8[l][m - 1] +
									C_m_right[l][m] * u_td_8[l][m + 1] + 
									dt / (2.0 * dz) * (H_z[l + 1][m] * v_r[l + 1][m] * r[l + 1][m] - // y -> r
													   H_z[l - 1][m] * v_r[l - 1][m] * r[l - 1][m]); // y -> r
			}
		}

		//printf("%lf %lf\n", u_td_1_next[50][50], C_0[50][50]);

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



		for (int l = 0; l < L_max; l++) {
			u_td_1_next[l][0] = u_td_1_next[l][1];
			u_td_2_next[l][0] = u_td_2_next[l][1];
			u_td_3_next[l][0] = u_td_3_next[l][1];
			u_td_4_next[l][0] = u_td_4_next[l][1];
			u_td_5_next[l][0] = u_td_5_next[l][1];
			u_td_6_next[l][0] = u_td_6_next[l][1];
			u_td_7_next[l][0] = u_td_7_next[l][1];
			u_td_8_next[l][0] = u_td_8_next[l][1];

			u_td_1_next[l][M_max] = u_td_1_next[l][M_max - 1];
			u_td_2_next[l][M_max] = u_td_2_next[l][M_max - 1];
			u_td_3_next[l][M_max] = u_td_3_next[l][M_max - 1];
			u_td_4_next[l][M_max] = u_td_4_next[l][M_max - 1];
			u_td_5_next[l][M_max] = u_td_5_next[l][M_max - 1];
			u_td_6_next[l][M_max] = u_td_6_next[l][M_max - 1];
			u_td_7_next[l][M_max] = u_td_7_next[l][M_max - 1];
			u_td_8_next[l][M_max] = u_td_8_next[l][M_max - 1];
		}

		
		// update parameters of problem

		for (int l = 0; l < L_max + 1; l++) {
			for (int m = 0; m < M_max + 1; m++) {
				rho[l][m] = u_td_1_next[l][m] / r[l][m];
				v_r[l][m] = u_td_3_next[l][m] / u_td_1_next[l][m];
				v_phi[l][m] = u_td_4_next[l][m] / u_td_1_next[l][m];
				v_z[l][m] = u_td_2_next[l][m] / u_td_1_next[l][m];
				v_y[l][m] = v_r[l][m] - v_z[l][m] * r_z[l][m];

				H_y[l][m] = u_td_8_next[l][m] / r[l][m];
				H_phi[l][m] = u_td_6_next[l][m];
				H_z[l][m] = u_td_7_next[l][m] / r[l][m];
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
	}

	// checkout

	/*
	std :: ofstream out0("Data/grid.txt");
	for (int m = M_max; m >= 0; m--) {
		for (int l = 0; l < L_max + 1; l++) {
			out0 << r[l][m] << "; ";
		}
		out0 << "\n\n";
	}
	
	std :: ofstream out1("Data/dr.txt");
	for (int m = M_max; m >= 0; m--) {
		for (int l = 0; l < L_max + 1; l++) {
			out1 << dr[l] << "   ";
		}
		out1 << "\n\n";
	}
	
	out0.close();
	out1.close();
	*/

	// output results in file

	std :: ofstream out2("Res.plt");
	int np = (L_max + 1) * (M_max + 1);
	int ne = L_max * M_max;
	double hfr;
	out2 << "VARIABLES=\n";
	out2 << "\"X\"\n\"Y\"\n\"Rho\"\n\"Vz\"\n\"Vr\"\n\"Vl\"\n\"Vphi\"\n\"Energy\"\n\"Hz\"\n\"Hr\"\n\"Hphi*r\"\n\"Hphi\"\n";
	out2 << "ZONE \n F=FEPOINT, ET=Quadrilateral, N=" << np << " E=" << ne << "\n ";
	for (int m = 0; m < M_max + 1; m++) {
		for (int l = 0; l < L_max + 1; l++) {
			hfr = H_phi[l][m] * r[l][m];
			out2 << l * dz << " " << r[l][m] << " " << rho[l][m] << " " << 
				v_z[l][m] << " " << v_r[l][m] << " " << std::sqrt(v_z[l][m] * v_z[l][m] + v_r[l][m] * v_r[l][m]) << " " << 
				v_phi[l][m] << " " << e[l][m] << " " << H_z[l][m] << " " << H_r[l][m] << " " << hfr << " " << H_phi[l][m] << "\n";
		}
	}

	int i1 = 0;
	int i2 = 0;
	int i3 = 0;
	int i4 = 0;
	
	for (int m = 0; m < M_max; m++) {
		for (int l = 0; l < L_max; l++) {
			i1 = l + m * (L_max + 1) + 1;
			i2 = l + 1 + m * (L_max + 1) + 1;
			i3 = l + 1 + (m + 1) * (L_max + 1) + 1;
			i4 = l + (m + 1) * (L_max + 1) + 1;
			out2 << i1 << " " << i2 << " " << i3 << " " << i4 << "\n";
		}
	}

	out2.close();

	// free memory

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
	
	memory_clearing(mu_l_left, L_max);
	memory_clearing(mu_l_right, L_max);
	memory_clearing(mu_m_left, L_max);
	memory_clearing(mu_m_right, L_max);
	
	memory_clearing(rho, L_max);
	memory_clearing(v_r, L_max);
	memory_clearing(v_phi, L_max);
	memory_clearing(e, L_max);
	memory_clearing(p, L_max);
	memory_clearing(P, L_max);
	memory_clearing(H_r, L_max);
	memory_clearing(H_z, L_max);
	
	memory_clearing(r, L_max);
	memory_clearing(r_z, L_max);
	delete [] R;
	delete [] dr;

	memory_clearing(v_y, L_max);
	memory_clearing(H_y, L_max);
	
	return 0;
}
