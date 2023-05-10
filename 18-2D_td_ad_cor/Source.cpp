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

double max(double a, double b, double c, double d, double e) {
	double res;

	res = max(max(a, b), max(max(c, d), e));

	return res;
}

double min(double a, double b) {
	double res;

	if (a < b) {
		res = a;
	} else {
		res = b;
	}

	return res;
}

double min(double a, double b, double c) {
	double res;

	res = min(min(a, b), c);

	return res;
}

double min(double a, double b, double c, double d, double e) {
	double res;

	res = min(min(a, b), min(min(c, d), e));

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
	int N_max = 10000;
	
	// u = {rho*r, rho*r*v_z, rho*r*v_r, rho*r*v_phi, rho*r*energy,  H_phi, H_z*r, H_r*r}

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

	double **u_ad_1;
	double **u_ad_2;
	double **u_ad_3;
	double **u_ad_4;
	double **u_ad_5;
	double **u_ad_6;
	double **u_ad_7;
	double **u_ad_8;
	memory_allocation(u_ad_1, L_max, M_max);
	memory_allocation(u_ad_2, L_max, M_max);
	memory_allocation(u_ad_3, L_max, M_max);
	memory_allocation(u_ad_4, L_max, M_max);
	memory_allocation(u_ad_5, L_max, M_max);
	memory_allocation(u_ad_6, L_max, M_max);
	memory_allocation(u_ad_7, L_max, M_max);
	memory_allocation(u_ad_8, L_max, M_max);

	double **u_a_1;
	double **u_a_2;
	double **u_a_3;
	double **u_a_4;
	double **u_a_5;
	double **u_a_6;
	double **u_a_7;
	double **u_a_8;
	memory_allocation(u_a_1, L_max, M_max);
	memory_allocation(u_a_2, L_max, M_max);
	memory_allocation(u_a_3, L_max, M_max);
	memory_allocation(u_a_4, L_max, M_max);
	memory_allocation(u_a_5, L_max, M_max);
	memory_allocation(u_a_6, L_max, M_max);
	memory_allocation(u_a_7, L_max, M_max);
	memory_allocation(u_a_8, L_max, M_max);

	double **u_b_1;
	double **u_b_2;
	double **u_b_3;
	double **u_b_4;
	double **u_b_5;
	double **u_b_6;
	double **u_b_7;
	double **u_b_8;
	memory_allocation(u_b_1, L_max, M_max);
	memory_allocation(u_b_2, L_max, M_max);
	memory_allocation(u_b_3, L_max, M_max);
	memory_allocation(u_b_4, L_max, M_max);
	memory_allocation(u_b_5, L_max, M_max);
	memory_allocation(u_b_6, L_max, M_max);
	memory_allocation(u_b_7, L_max, M_max);
	memory_allocation(u_b_8, L_max, M_max);

	double **u_max_1;
	double **u_max_2;
	double **u_max_3;
	double **u_max_4;
	double **u_max_5;
	double **u_max_6;
	double **u_max_7;
	double **u_max_8;
	memory_allocation(u_max_1, L_max, M_max);
	memory_allocation(u_max_2, L_max, M_max);
	memory_allocation(u_max_3, L_max, M_max);
	memory_allocation(u_max_4, L_max, M_max);
	memory_allocation(u_max_5, L_max, M_max);
	memory_allocation(u_max_6, L_max, M_max);
	memory_allocation(u_max_7, L_max, M_max);
	memory_allocation(u_max_8, L_max, M_max);

	double **u_min_1;
	double **u_min_2;
	double **u_min_3;
	double **u_min_4;
	double **u_min_5;
	double **u_min_6;
	double **u_min_7;
	double **u_min_8;
	memory_allocation(u_min_1, L_max, M_max);
	memory_allocation(u_min_2, L_max, M_max);
	memory_allocation(u_min_3, L_max, M_max);
	memory_allocation(u_min_4, L_max, M_max);
	memory_allocation(u_min_5, L_max, M_max);
	memory_allocation(u_min_6, L_max, M_max);
	memory_allocation(u_min_7, L_max, M_max);
	memory_allocation(u_min_8, L_max, M_max);

	double **u_cor_1;
	double **u_cor_2;
	double **u_cor_3;
	double **u_cor_4;
	double **u_cor_5;
	double **u_cor_6;
	double **u_cor_7;
	double **u_cor_8;
	memory_allocation(u_cor_1, L_max, M_max);
	memory_allocation(u_cor_2, L_max, M_max);
	memory_allocation(u_cor_3, L_max, M_max);
	memory_allocation(u_cor_4, L_max, M_max);
	memory_allocation(u_cor_5, L_max, M_max);
	memory_allocation(u_cor_6, L_max, M_max);
	memory_allocation(u_cor_7, L_max, M_max);
	memory_allocation(u_cor_8, L_max, M_max);

	// constants
	
	double k = 0.5;
	double gamma = 1.67;
	double beta = 1.0;
	double H_z0 = 0;

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

	double **mu_ad_l_left;
	double **mu_ad_l_right;
	double **mu_ad_m_left;
	double **mu_ad_m_right;
	memory_allocation(mu_ad_l_left, L_max, M_max);
	memory_allocation(mu_ad_l_right, L_max, M_max);
	memory_allocation(mu_ad_m_left, L_max, M_max);
	memory_allocation(mu_ad_m_right, L_max, M_max);

	double **P_plus_1;
	double **P_plus_2;
	double **P_plus_3;
	double **P_plus_4;
	double **P_plus_5;
	double **P_plus_6;
	double **P_plus_7;
	double **P_plus_8;
	memory_allocation(P_plus_1, L_max, M_max);
	memory_allocation(P_plus_2, L_max, M_max);
	memory_allocation(P_plus_3, L_max, M_max);
	memory_allocation(P_plus_4, L_max, M_max);
	memory_allocation(P_plus_5, L_max, M_max);
	memory_allocation(P_plus_6, L_max, M_max);
	memory_allocation(P_plus_7, L_max, M_max);
	memory_allocation(P_plus_8, L_max, M_max);

	double **P_minus_1;
	double **P_minus_2;
	double **P_minus_3;
	double **P_minus_4;
	double **P_minus_5;
	double **P_minus_6;
	double **P_minus_7;
	double **P_minus_8;
	memory_allocation(P_minus_1, L_max, M_max);
	memory_allocation(P_minus_2, L_max, M_max);
	memory_allocation(P_minus_3, L_max, M_max);
	memory_allocation(P_minus_4, L_max, M_max);
	memory_allocation(P_minus_5, L_max, M_max);
	memory_allocation(P_minus_6, L_max, M_max);
	memory_allocation(P_minus_7, L_max, M_max);
	memory_allocation(P_minus_8, L_max, M_max);

	double **Q_plus_1;
	double **Q_plus_2;
	double **Q_plus_3;
	double **Q_plus_4;
	double **Q_plus_5;
	double **Q_plus_6;
	double **Q_plus_7;
	double **Q_plus_8;
	memory_allocation(Q_plus_1, L_max, M_max);
	memory_allocation(Q_plus_2, L_max, M_max);
	memory_allocation(Q_plus_3, L_max, M_max);
	memory_allocation(Q_plus_4, L_max, M_max);
	memory_allocation(Q_plus_5, L_max, M_max);
	memory_allocation(Q_plus_6, L_max, M_max);
	memory_allocation(Q_plus_7, L_max, M_max);
	memory_allocation(Q_plus_8, L_max, M_max);

	double **Q_minus_1;
	double **Q_minus_2;
	double **Q_minus_3;
	double **Q_minus_4;
	double **Q_minus_5;
	double **Q_minus_6;
	double **Q_minus_7;
	double **Q_minus_8;
	memory_allocation(Q_minus_1, L_max, M_max);
	memory_allocation(Q_minus_2, L_max, M_max);
	memory_allocation(Q_minus_3, L_max, M_max);
	memory_allocation(Q_minus_4, L_max, M_max);
	memory_allocation(Q_minus_5, L_max, M_max);
	memory_allocation(Q_minus_6, L_max, M_max);
	memory_allocation(Q_minus_7, L_max, M_max);
	memory_allocation(Q_minus_8, L_max, M_max);

	double **R_plus_1;
	double **R_plus_2;
	double **R_plus_3;
	double **R_plus_4;
	double **R_plus_5;
	double **R_plus_6;
	double **R_plus_7;
	double **R_plus_8;
	memory_allocation(R_plus_1, L_max, M_max);
	memory_allocation(R_plus_2, L_max, M_max);
	memory_allocation(R_plus_3, L_max, M_max);
	memory_allocation(R_plus_4, L_max, M_max);
	memory_allocation(R_plus_5, L_max, M_max);
	memory_allocation(R_plus_6, L_max, M_max);
	memory_allocation(R_plus_7, L_max, M_max);
	memory_allocation(R_plus_8, L_max, M_max);

	double **R_minus_1;
	double **R_minus_2;
	double **R_minus_3;
	double **R_minus_4;
	double **R_minus_5;
	double **R_minus_6;
	double **R_minus_7;
	double **R_minus_8;
	memory_allocation(R_minus_1, L_max, M_max);
	memory_allocation(R_minus_2, L_max, M_max);
	memory_allocation(R_minus_3, L_max, M_max);
	memory_allocation(R_minus_4, L_max, M_max);
	memory_allocation(R_minus_5, L_max, M_max);
	memory_allocation(R_minus_6, L_max, M_max);
	memory_allocation(R_minus_7, L_max, M_max);
	memory_allocation(R_minus_8, L_max, M_max);

	double **COR_l_left_1;
	double **COR_l_left_2;
	double **COR_l_left_3;
	double **COR_l_left_4;
	double **COR_l_left_5;
	double **COR_l_left_6;
	double **COR_l_left_7;
	double **COR_l_left_8;
	memory_allocation(COR_l_left_1, L_max, M_max);
	memory_allocation(COR_l_left_2, L_max, M_max);
	memory_allocation(COR_l_left_3, L_max, M_max);
	memory_allocation(COR_l_left_4, L_max, M_max);
	memory_allocation(COR_l_left_5, L_max, M_max);
	memory_allocation(COR_l_left_6, L_max, M_max);
	memory_allocation(COR_l_left_7, L_max, M_max);
	memory_allocation(COR_l_left_8, L_max, M_max);

	double **COR_l_right_1;
	double **COR_l_right_2;
	double **COR_l_right_3;
	double **COR_l_right_4;
	double **COR_l_right_5;
	double **COR_l_right_6;
	double **COR_l_right_7;
	double **COR_l_right_8;
	memory_allocation(COR_l_right_1, L_max, M_max);
	memory_allocation(COR_l_right_2, L_max, M_max);
	memory_allocation(COR_l_right_3, L_max, M_max);
	memory_allocation(COR_l_right_4, L_max, M_max);
	memory_allocation(COR_l_right_5, L_max, M_max);
	memory_allocation(COR_l_right_6, L_max, M_max);
	memory_allocation(COR_l_right_7, L_max, M_max);
	memory_allocation(COR_l_right_8, L_max, M_max);

	double **COR_m_left_1;
	double **COR_m_left_2;
	double **COR_m_left_3;
	double **COR_m_left_4;
	double **COR_m_left_5;
	double **COR_m_left_6;
	double **COR_m_left_7;
	double **COR_m_left_8;
	memory_allocation(COR_m_left_1, L_max, M_max);
	memory_allocation(COR_m_left_2, L_max, M_max);
	memory_allocation(COR_m_left_3, L_max, M_max);
	memory_allocation(COR_m_left_4, L_max, M_max);
	memory_allocation(COR_m_left_5, L_max, M_max);
	memory_allocation(COR_m_left_6, L_max, M_max);
	memory_allocation(COR_m_left_7, L_max, M_max);
	memory_allocation(COR_m_left_8, L_max, M_max);

	double **COR_m_right_1;
	double **COR_m_right_2;
	double **COR_m_right_3;
	double **COR_m_right_4;
	double **COR_m_right_5;
	double **COR_m_right_6;
	double **COR_m_right_7;
	double **COR_m_right_8;
	memory_allocation(COR_m_right_1, L_max, M_max);
	memory_allocation(COR_m_right_2, L_max, M_max);
	memory_allocation(COR_m_right_3, L_max, M_max);
	memory_allocation(COR_m_right_4, L_max, M_max);
	memory_allocation(COR_m_right_5, L_max, M_max);
	memory_allocation(COR_m_right_6, L_max, M_max);
	memory_allocation(COR_m_right_7, L_max, M_max);
	memory_allocation(COR_m_right_8, L_max, M_max);

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

	for (int l = 0; l < L_max + 1; l++) {
		for (int m = 0; m < M_max + 1; m++) {
			rho[l][m] = 1.0;
			v_z[l][m] = 0.1;
			v_r[l][m] = 0.1;
			v_phi[l][m] = 0;
			H_phi[l][m] = (1 - 0.9 * l * dz) * r_0 / r[l][m];
			H_z[l][m] = H_z0;
			H_r[l][m] = 0;
			
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

		// constants for current layer
		
		for (int l = 1; l < L_max; l++) {
			for (int m = 1; m < M_max; m++) {
				mu_l_left[l][m] = mu_0_l + (fabs(v_z[l - 1][m]) + fabs(v_z[l][m])) / 4.0;
				mu_l_right[l][m] = mu_0_l + (fabs(v_z[l][m]) + fabs(v_z[l + 1][m])) / 4.0;
				mu_m_left[l][m] = mu_0_m + (fabs(v_y[l][m - 1]) + fabs(v_y[l][m])) / 4.0;
				mu_m_right[l][m] = mu_0_m + (fabs(v_y[l][m]) + fabs(v_y[l][m + 1])) / 4.0;

				mu_ad_l_left[l][m] = mu_l_left[l][m] - fabs(v_z[l - 1][m] + v_z[l][m]) / 4.0;
				mu_ad_l_right[l][m] = mu_l_right[l][m] - fabs(v_z[l][m] + v_z[l + 1][m]) / 4.0;
				mu_ad_m_left[l][m] = mu_m_left[l][m] - fabs(v_y[l][m - 1] + v_y[l][m]) / 4.0;
				mu_ad_m_right[l][m] = mu_m_right[l][m] - fabs(v_y[l][m] + v_y[l][m + 1]) / 4.0;
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

		// constants for current layer

		for (int l = 1; l < L_max; l++) {
			for (int m = 1; m < M_max; m++) {
				C_0[l][m] = 1 - dt / dz * (v_z[l + 1][m] - v_z[l - 1][m]) / 4.0 - dt / dr[l] * (v_y[l][m + 1] - v_y[l][m - 1]) / 4.0
							  - dt / dz * (mu_l_left[l][m] + mu_l_right[l][m]) - dt / dr[l] * (mu_m_left[l][m] + mu_m_right[l][m]);
				C_l_left[l][m] = dt / dz * ((v_z[l - 1][m] + v_z[l][m]) / 4.0 + mu_l_left[l][m]);
				C_l_right[l][m] = dt / dz * (- (v_z[l][m] + v_z[l + 1][m]) / 4.0 + mu_l_right[l][m]);
				C_m_left[l][m] = dt / dr[l] * ((v_r[l][m - 1] + v_r[l][m]) / 4.0 + mu_m_left[l][m]);
				C_m_right[l][m] = dt / dr[l] * (- (v_r[l][m] + v_r[l][m + 1]) / 4.0 + mu_m_right[l][m]);
			}
		}

		// layer checkout

		printf("n=%d   rho[40][20]=%lf   v_z[40][20]=%lf   v_phi[40][20]=%lf\n", n, rho[40][20], v_z[40][20], v_phi[40][20]);

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

				u_td_3_next[l][m] = C_0[l][m] * u_td_3[l][m] + 
									C_l_left[l][m] * u_td_3[l - 1][m] + 
									C_l_right[l][m] * u_td_3[l + 1][m] +
									C_m_left[l][m] * u_td_3[l][m - 1] +
									C_m_right[l][m] * u_td_3[l][m + 1] + 
									dt / (2.0 * dz) * (H_z[l + 1][m] * H_r[l + 1][m] * r[l + 1][m] - 
													   H_z[l - 1][m] * H_r[l - 1][m] * r[l - 1][m]) -
									dt / (2.0 * dr[l]) * ((P[l][m + 1] - pow(H_r[l][m + 1], 2)) * r[l][m + 1] - 
													      (P[l][m - 1] - pow(H_r[l][m - 1], 2)) * r[l][m - 1]) + 
									dt * (rho[l][m] * pow(v_phi[l][m], 2) + P[l][m] - pow(H_phi[l][m], 2));

				u_td_4_next[l][m] = C_0[l][m] * u_td_4[l][m] + 
									C_l_left[l][m] * u_td_4[l - 1][m] + 
									C_l_right[l][m] * u_td_4[l + 1][m] +
									C_m_left[l][m] * u_td_4[l][m - 1] +
									C_m_right[l][m] * u_td_4[l][m + 1] + 
									dt / (2.0 * dz) * (H_phi[l + 1][m] * H_z[l + 1][m] * r[l + 1][m] - 
													   H_phi[l - 1][m] * H_z[l - 1][m] * r[l - 1][m]) +
									dt / (2.0 * dr[l]) * (H_phi[l][m + 1] * H_r[l][m + 1] * r[l][m + 1] - 
													      H_phi[l][m - 1] * H_r[l][m - 1] * r[l][m - 1]) + 
									dt * (- rho[l][m] * r[l][m] * v_phi[l][m] + H_phi[l][m] * H_r[l][m]);

				u_td_5_next[l][m] = C_0[l][m] * u_td_5[l][m] + 
									C_l_left[l][m] * u_td_5[l - 1][m] + 
									C_l_right[l][m] * u_td_5[l + 1][m] +
									C_m_left[l][m] * u_td_5[l][m - 1] +
									C_m_right[l][m] * u_td_5[l][m + 1] -
									dt * p[l][m] * (1.0 / (2.0 * dz) * (v_z[l + 1][m] * r[l + 1][m] -
																		v_z[l - 1][m] * r[l - 1][m]) + 
													1.0 / (2.0 * dr[l]) * (v_r[l][m + 1] * r[l][m + 1] - 
																		   v_r[l][m - 1] * r[l][m - 1]));

				u_td_6_next[l][m] = C_0[l][m] * u_td_6[l][m] + 
									C_l_left[l][m] * u_td_6[l - 1][m] + 
									C_l_right[l][m] * u_td_6[l + 1][m] +
									C_m_left[l][m] * u_td_6[l][m - 1] +
									C_m_right[l][m] * u_td_6[l][m + 1] +
									dt / (2.0 * dz) * (H_z[l + 1][m] * v_phi[l + 1][m] - 
													   H_z[l - 1][m] * v_phi[l - 1][m]) + 
									dt / (2.0 * dr[l]) * (H_r[l][m + 1] * v_phi[l][m + 1] - 
													      H_r[l][m - 1] * v_phi[l][m - 1]);
				
				u_td_7_next[l][m] = C_0[l][m] * u_td_7[l][m] + 
									C_l_left[l][m] * u_td_7[l - 1][m] + 
									C_l_right[l][m] * u_td_7[l + 1][m] +
									C_m_left[l][m] * u_td_7[l][m - 1] +
									C_m_right[l][m] * u_td_7[l][m + 1] +
									dt / (2.0 * dr[l]) * (H_r[l][m + 1] * v_z[l][m + 1] * r[l][m + 1] - 
													      H_r[l][m - 1] * v_z[l][m - 1] * r[l][m - 1]);

				u_td_8_next[l][m] = C_0[l][m] * u_td_8[l][m] + 
									C_l_left[l][m] * u_td_8[l - 1][m] + 
									C_l_right[l][m] * u_td_8[l + 1][m] +
									C_m_left[l][m] * u_td_8[l][m - 1] +
									C_m_right[l][m] * u_td_8[l][m + 1] + 
									dt / (2.0 * dz) * (H_z[l + 1][m] * v_r[l + 1][m] * r[l + 1][m] - 
													   H_z[l - 1][m] * v_r[l - 1][m] * r[l - 1][m]);
			}
		}

		// left boundary condition for current layer

		for (int m = 0; m < M_max + 1; m++) {
			rho[0][m] = 1.0;
			v_phi[0][m] = 0;
			// v_r[0][m] = v_z[0][m] * r_z[0][m];
			H_phi[0][m] = r_0 / r[0][m];
			H_z[0][m] = H_z0;
			// H_r[0][m] = H_z[0][m] * r_z[0][m];
			e[0][m] = beta / (2.0 * (gamma - 1.0)) * pow(rho[0][m], gamma - 1.0);
		}

		// left boundary condition for current layer for transport part

		for (int m = 0; m < M_max + 1; m++) {
			u_td_1_next[0][m] = rho[0][m] * r[0][m];
			u_td_2_next[0][m] = u_td_2_next[1][m];
			// u_td_2_next[0][m] = 2 * u_td_2_next[1][m] - u_td_2_next[2][m];
			u_td_3_next[0][m] = u_td_2_next[0][m] * r_z[0][m];
			u_td_4_next[0][m] = rho[0][m] * v_phi[0][m] * r[0][m];
			u_td_5_next[0][m] = rho[0][m] * e[0][m] * r[0][m];
			u_td_6_next[0][m] = H_phi[0][m];
			u_td_7_next[0][m] = H_z[0][m] * r[0][m];
			u_td_8_next[0][m] = u_td_7_next[0][m] * r_z[0][m];
		}

		// up and down boundary condition for current layer for transport part

		for (int l = 1; l < L_max + 1; l++) {
			u_td_1_next[l][0] = u_td_1_next[l][1];
			u_td_2_next[l][0] = u_td_2_next[l][1];
			u_td_3_next[l][0] = u_td_2_next[l][0] * r_z[l][0];
			u_td_4_next[l][0] = u_td_4_next[l][1];
			u_td_5_next[l][0] = u_td_5_next[l][1];
			u_td_6_next[l][0] = u_td_6_next[l][1];
			u_td_7_next[l][0] = u_td_7_next[l][1];
			u_td_8_next[l][0] = u_td_7_next[l][0] * r_z[l][0];

			u_td_1_next[l][M_max] = u_td_1_next[l][M_max - 1];
			u_td_2_next[l][M_max] = u_td_2_next[l][M_max - 1];
			u_td_3_next[l][M_max] = u_td_2_next[l][M_max] * r_z[l][M_max];
			u_td_4_next[l][M_max] = u_td_4_next[l][M_max - 1];
			u_td_5_next[l][M_max] = u_td_5_next[l][M_max - 1];
			u_td_6_next[l][M_max] = u_td_6_next[l][M_max - 1];
			u_td_7_next[l][M_max] = u_td_7_next[l][M_max - 1];
			u_td_8_next[l][M_max] = u_td_7_next[l][M_max] * r_z[l][M_max];
		}

		// right boundary condition for current layer for transport part

		for (int m = 1; m < M_max; m++) {
			u_td_1_next[L_max][m] = u_td_1_next[L_max - 1][m];
			u_td_2_next[L_max][m] = u_td_2_next[L_max - 1][m];
			u_td_3_next[L_max][m] = u_td_3_next[L_max - 1][m];
			u_td_4_next[L_max][m] = u_td_4_next[L_max - 1][m];
			u_td_5_next[L_max][m] = u_td_5_next[L_max - 1][m];
			u_td_6_next[L_max][m] = u_td_6_next[L_max - 1][m];
			u_td_7_next[L_max][m] = u_td_7_next[L_max - 1][m];
			u_td_8_next[L_max][m] = u_td_8_next[L_max - 1][m];
		}

		// filling central points of grid for anti-diffusion part

		for (int l = 1; l < L_max; l++) {
			for (int m = 1; m < M_max; m++) {
				u_ad_1[l][m] = u_td_1_next[l][m] * (1.0 + dt / dz * (mu_ad_l_left[l][m] + mu_ad_l_right[l][m]) + 
														  dt / dr[l] * (mu_ad_m_left[l][m] + mu_ad_m_right[l][m])) - 
							   dt / dz * (mu_ad_l_left[l][m] * u_td_1_next[l - 1][m] + mu_ad_l_right[l][m] * u_td_1_next[l + 1][m]) - 
							   dt / dr[l] * (mu_ad_m_left[l][m] * u_td_1_next[l][m - 1] + mu_ad_m_right[l][m] * u_td_1_next[l][m + 1]);

				u_ad_2[l][m] = u_td_2_next[l][m] * (1.0 + dt / dz * (mu_ad_l_left[l][m] + mu_ad_l_right[l][m]) + 
														  dt / dr[l] * (mu_ad_m_left[l][m] + mu_ad_m_right[l][m])) - 
							   dt / dz * (mu_ad_l_left[l][m] * u_td_2_next[l - 1][m] + mu_ad_l_right[l][m] * u_td_2_next[l + 1][m]) - 
							   dt / dr[l] * (mu_ad_m_left[l][m] * u_td_2_next[l][m - 1] + mu_ad_m_right[l][m] * u_td_2_next[l][m + 1]);

				u_ad_3[l][m] = u_td_3_next[l][m] * (1.0 + dt / dz * (mu_ad_l_left[l][m] + mu_ad_l_right[l][m]) + 
														  dt / dr[l] * (mu_ad_m_left[l][m] + mu_ad_m_right[l][m])) - 
							   dt / dz * (mu_ad_l_left[l][m] * u_td_3_next[l - 1][m] + mu_ad_l_right[l][m] * u_td_3_next[l + 1][m]) - 
							   dt / dr[l] * (mu_ad_m_left[l][m] * u_td_3_next[l][m - 1] + mu_ad_m_right[l][m] * u_td_3_next[l][m + 1]);

				u_ad_4[l][m] = u_td_4_next[l][m] * (1.0 + dt / dz * (mu_ad_l_left[l][m] + mu_ad_l_right[l][m]) + 
														  dt / dr[l] * (mu_ad_m_left[l][m] + mu_ad_m_right[l][m])) - 
							   dt / dz * (mu_ad_l_left[l][m] * u_td_4_next[l - 1][m] + mu_ad_l_right[l][m] * u_td_4_next[l + 1][m]) - 
							   dt / dr[l] * (mu_ad_m_left[l][m] * u_td_4_next[l][m - 1] + mu_ad_m_right[l][m] * u_td_4_next[l][m + 1]);

				u_ad_5[l][m] = u_td_5_next[l][m] * (1.0 + dt / dz * (mu_ad_l_left[l][m] + mu_ad_l_right[l][m]) + 
														  dt / dr[l] * (mu_ad_m_left[l][m] + mu_ad_m_right[l][m])) - 
							   dt / dz * (mu_ad_l_left[l][m] * u_td_5_next[l - 1][m] + mu_ad_l_right[l][m] * u_td_5_next[l + 1][m]) - 
							   dt / dr[l] * (mu_ad_m_left[l][m] * u_td_5_next[l][m - 1] + mu_ad_m_right[l][m] * u_td_5_next[l][m + 1]);

				u_ad_6[l][m] = u_td_6_next[l][m] * (1.0 + dt / dz * (mu_ad_l_left[l][m] + mu_ad_l_right[l][m]) + 
														  dt / dr[l] * (mu_ad_m_left[l][m] + mu_ad_m_right[l][m])) - 
							   dt / dz * (mu_ad_l_left[l][m] * u_td_6_next[l - 1][m] + mu_ad_l_right[l][m] * u_td_6_next[l + 1][m]) - 
							   dt / dr[l] * (mu_ad_m_left[l][m] * u_td_6_next[l][m - 1] + mu_ad_m_right[l][m] * u_td_6_next[l][m + 1]);

				u_ad_7[l][m] = u_td_7_next[l][m] * (1.0 + dt / dz * (mu_ad_l_left[l][m] + mu_ad_l_right[l][m]) + 
														  dt / dr[l] * (mu_ad_m_left[l][m] + mu_ad_m_right[l][m])) - 
							   dt / dz * (mu_ad_l_left[l][m] * u_td_7_next[l - 1][m] + mu_ad_l_right[l][m] * u_td_7_next[l + 1][m]) - 
							   dt / dr[l] * (mu_ad_m_left[l][m] * u_td_7_next[l][m - 1] + mu_ad_m_right[l][m] * u_td_7_next[l][m + 1]);

				u_ad_8[l][m] = u_td_8_next[l][m] * (1.0 + dt / dz * (mu_ad_l_left[l][m] + mu_ad_l_right[l][m]) + 
														  dt / dr[l] * (mu_ad_m_left[l][m] + mu_ad_m_right[l][m])) - 
							   dt / dz * (mu_ad_l_left[l][m] * u_td_8_next[l - 1][m] + mu_ad_l_right[l][m] * u_td_8_next[l + 1][m]) - 
							   dt / dr[l] * (mu_ad_m_left[l][m] * u_td_8_next[l][m - 1] + mu_ad_m_right[l][m] * u_td_8_next[l][m + 1]);
			}
		}

		// left boundary condition for current layer for anti-diffusion part

		for (int m = 0; m < M_max + 1; m++) {
			//u_ad_1[0][m] = rho[0][m] * r[0][m];
			//u_ad_2[0][m] = u_ad_2[1][m];
			//// u_ad_2[0][m] = 2 * u_ad_2[1][m] - u_ad_2[2][m];
			//u_ad_3[0][m] = u_ad_2[0][m] * r_z[0][m];
			//u_ad_4[0][m] = rho[0][m] * v_phi[0][m] * r[0][m];
			//u_ad_5[0][m] = rho[0][m] * e[0][m] * r[0][m];
			//u_ad_6[0][m] = H_phi[0][m];
			//u_ad_7[0][m] = H_z[0][m] * r[0][m];
			//u_ad_8[0][m] = u_ad_7[0][m] * r_z[0][m];

			u_ad_1[0][m] = u_td_1_next[0][m];
			u_ad_2[0][m] = u_td_2_next[0][m];
			u_ad_3[0][m] = u_td_3_next[0][m];
			u_ad_4[0][m] = u_td_4_next[0][m];
			u_ad_5[0][m] = u_td_5_next[0][m];
			u_ad_6[0][m] = u_td_6_next[0][m];
			u_ad_7[0][m] = u_td_7_next[0][m];
			u_ad_8[0][m] = u_td_8_next[0][m];
		}

		// up and down boundary condition for current layer for anti-diffusion part

		for (int l = 1; l < L_max + 1; l++) {
			/*u_ad_1[l][0] = u_ad_1[l][1];
			u_ad_2[l][0] = u_ad_2[l][1];
			u_ad_3[l][0] = u_ad_2[l][0] * r_z[l][0];
			u_ad_4[l][0] = u_ad_4[l][1];
			u_ad_5[l][0] = u_ad_5[l][1];
			u_ad_6[l][0] = u_ad_6[l][1];
			u_ad_7[l][0] = u_ad_7[l][1];
			u_ad_8[l][0] = u_ad_7[l][0] * r_z[l][0];

			u_ad_1[l][M_max] = u_ad_1[l][M_max - 1];
			u_ad_2[l][M_max] = u_ad_2[l][M_max - 1];
			u_ad_3[l][M_max] = u_ad_2[l][M_max] * r_z[l][M_max];
			u_ad_4[l][M_max] = u_ad_4[l][M_max - 1];
			u_ad_5[l][M_max] = u_ad_5[l][M_max - 1];
			u_ad_6[l][M_max] = u_ad_6[l][M_max - 1];
			u_ad_7[l][M_max] = u_ad_7[l][M_max - 1];
			u_ad_8[l][M_max] = u_ad_7[l][M_max] * r_z[l][M_max];*/

			u_ad_1[l][0] = u_td_1_next[l][0];
			u_ad_2[l][0] = u_td_2_next[l][0];
			u_ad_3[l][0] = u_td_3_next[l][0];
			u_ad_4[l][0] = u_td_4_next[l][0];
			u_ad_5[l][0] = u_td_5_next[l][0];
			u_ad_6[l][0] = u_td_6_next[l][0];
			u_ad_7[l][0] = u_td_7_next[l][0];
			u_ad_8[l][0] = u_td_8_next[l][0];

			u_ad_1[l][M_max] = u_td_1_next[l][M_max];
			u_ad_2[l][M_max] = u_td_2_next[l][M_max];
			u_ad_3[l][M_max] = u_td_3_next[l][M_max];
			u_ad_4[l][M_max] = u_td_4_next[l][M_max];
			u_ad_5[l][M_max] = u_td_5_next[l][M_max];
			u_ad_6[l][M_max] = u_td_6_next[l][M_max];
			u_ad_7[l][M_max] = u_td_7_next[l][M_max];
			u_ad_8[l][M_max] = u_td_8_next[l][M_max];
		}

		// right boundary condition for current layer for anti-diffusion part

		for (int m = 1; m < M_max; m++) {
			/*u_ad_1[L_max][m] = u_ad_1[L_max - 1][m];
			u_ad_2[L_max][m] = u_ad_2[L_max - 1][m];
			u_ad_3[L_max][m] = u_ad_3[L_max - 1][m];
			u_ad_4[L_max][m] = u_ad_4[L_max - 1][m];
			u_ad_5[L_max][m] = u_ad_5[L_max - 1][m];
			u_ad_6[L_max][m] = u_ad_6[L_max - 1][m];
			u_ad_7[L_max][m] = u_ad_7[L_max - 1][m];
			u_ad_8[L_max][m] = u_ad_8[L_max - 1][m];*/

			u_ad_1[L_max][m] = u_td_1_next[L_max][m];
			u_ad_2[L_max][m] = u_td_2_next[L_max][m];
			u_ad_3[L_max][m] = u_td_3_next[L_max][m];
			u_ad_4[L_max][m] = u_td_4_next[L_max][m];
			u_ad_5[L_max][m] = u_td_5_next[L_max][m];
			u_ad_6[L_max][m] = u_td_6_next[L_max][m];
			u_ad_7[L_max][m] = u_td_7_next[L_max][m];
			u_ad_8[L_max][m] = u_td_8_next[L_max][m];
		}

		// constants for current layer for correction part
		
		for (int l = 0; l < L_max + 1; l++) {
			for (int m = 0; m < M_max + 1; m++) {
				u_a_1[l][m] = max(u_td_1[l][m], u_td_1_next[l][m]);
				u_a_2[l][m] = max(u_td_2[l][m], u_td_2_next[l][m]);
				u_a_3[l][m] = max(u_td_3[l][m], u_td_3_next[l][m]);
				u_a_4[l][m] = max(u_td_4[l][m], u_td_4_next[l][m]);
				u_a_5[l][m] = max(u_td_5[l][m], u_td_5_next[l][m]);
				u_a_6[l][m] = max(u_td_6[l][m], u_td_6_next[l][m]);
				u_a_7[l][m] = max(u_td_7[l][m], u_td_7_next[l][m]);
				u_a_8[l][m] = max(u_td_8[l][m], u_td_8_next[l][m]);

				u_b_1[l][m] = min(u_td_1[l][m], u_td_1_next[l][m]);
				u_b_2[l][m] = min(u_td_2[l][m], u_td_2_next[l][m]);
				u_b_3[l][m] = min(u_td_3[l][m], u_td_3_next[l][m]);
				u_b_4[l][m] = min(u_td_4[l][m], u_td_4_next[l][m]);
				u_b_5[l][m] = min(u_td_5[l][m], u_td_5_next[l][m]);
				u_b_6[l][m] = min(u_td_6[l][m], u_td_6_next[l][m]);
				u_b_7[l][m] = min(u_td_7[l][m], u_td_7_next[l][m]);
				u_b_8[l][m] = min(u_td_8[l][m], u_td_8_next[l][m]);
			}
		}

		for (int l = 1; l < L_max; l++) {
			for (int m = 1; m < M_max; m++) {
				u_max_1[l][m] = max(u_a_1[l - 1][m], u_a_1[l + 1][m], u_a_1[l][m - 1], u_a_1[l][m + 1], u_a_1[l][m]);
				u_max_2[l][m] = max(u_a_2[l - 1][m], u_a_2[l + 1][m], u_a_2[l][m - 1], u_a_2[l][m + 1], u_a_2[l][m]);
				u_max_3[l][m] = max(u_a_3[l - 1][m], u_a_3[l + 1][m], u_a_3[l][m - 1], u_a_3[l][m + 1], u_a_3[l][m]);
				u_max_4[l][m] = max(u_a_4[l - 1][m], u_a_4[l + 1][m], u_a_4[l][m - 1], u_a_4[l][m + 1], u_a_4[l][m]);
				u_max_5[l][m] = max(u_a_5[l - 1][m], u_a_5[l + 1][m], u_a_5[l][m - 1], u_a_5[l][m + 1], u_a_5[l][m]);
				u_max_6[l][m] = max(u_a_6[l - 1][m], u_a_6[l + 1][m], u_a_6[l][m - 1], u_a_6[l][m + 1], u_a_6[l][m]);
				u_max_7[l][m] = max(u_a_7[l - 1][m], u_a_7[l + 1][m], u_a_7[l][m - 1], u_a_7[l][m + 1], u_a_7[l][m]);
				u_max_8[l][m] = max(u_a_8[l - 1][m], u_a_8[l + 1][m], u_a_8[l][m - 1], u_a_8[l][m + 1], u_a_8[l][m]);

				u_min_1[l][m] = min(u_b_1[l - 1][m], u_b_1[l + 1][m], u_b_1[l][m - 1], u_b_1[l][m + 1], u_b_1[l][m]);
				u_min_2[l][m] = min(u_b_2[l - 1][m], u_b_2[l + 1][m], u_b_2[l][m - 1], u_b_2[l][m + 1], u_b_2[l][m]);
				u_min_3[l][m] = min(u_b_3[l - 1][m], u_b_3[l + 1][m], u_b_3[l][m - 1], u_b_3[l][m + 1], u_b_3[l][m]);
				u_min_4[l][m] = min(u_b_4[l - 1][m], u_b_4[l + 1][m], u_b_4[l][m - 1], u_b_4[l][m + 1], u_b_4[l][m]);
				u_min_5[l][m] = min(u_b_5[l - 1][m], u_b_5[l + 1][m], u_b_5[l][m - 1], u_b_5[l][m + 1], u_b_5[l][m]);
				u_min_6[l][m] = min(u_b_6[l - 1][m], u_b_6[l + 1][m], u_b_6[l][m - 1], u_b_6[l][m + 1], u_b_6[l][m]);
				u_min_7[l][m] = min(u_b_7[l - 1][m], u_b_7[l + 1][m], u_b_7[l][m - 1], u_b_7[l][m + 1], u_b_7[l][m]);
				u_min_8[l][m] = min(u_b_8[l - 1][m], u_b_8[l + 1][m], u_b_8[l][m - 1], u_b_8[l][m + 1], u_b_8[l][m]);
			}
		}

		// filling P, Q and R constants

		for (int l = 1; l < L_max; l++) {
			for (int m = 1; m < M_max; m++) {
				P_plus_1[l][m] = max(0, mu_ad_l_left[l][m] * (u_td_1_next[l][m] - u_td_1_next[l - 1][m])) - 
								 min(0, mu_ad_l_right[l][m] * (u_td_1_next[l + 1][m] - u_td_1_next[l][m])) + 
								 max(0, mu_ad_m_left[l][m] * (u_td_1_next[l][m] - u_td_1_next[l][m - 1])) - 
								 min(0, mu_ad_m_right[l][m] * (u_td_1_next[l][m + 1] - u_td_1_next[l][m]));
				P_plus_2[l][m] = max(0, mu_ad_l_left[l][m] * (u_td_2_next[l][m] - u_td_2_next[l - 1][m])) - 
								 min(0, mu_ad_l_right[l][m] * (u_td_2_next[l + 1][m] - u_td_2_next[l][m])) + 
								 max(0, mu_ad_m_left[l][m] * (u_td_2_next[l][m] - u_td_2_next[l][m - 1])) - 
								 min(0, mu_ad_m_right[l][m] * (u_td_2_next[l][m + 1] - u_td_2_next[l][m]));
				P_plus_3[l][m] = max(0, mu_ad_l_left[l][m] * (u_td_3_next[l][m] - u_td_3_next[l - 1][m])) - 
								 min(0, mu_ad_l_right[l][m] * (u_td_3_next[l + 1][m] - u_td_3_next[l][m])) + 
								 max(0, mu_ad_m_left[l][m] * (u_td_3_next[l][m] - u_td_3_next[l][m - 1])) - 
								 min(0, mu_ad_m_right[l][m] * (u_td_3_next[l][m + 1] - u_td_3_next[l][m]));
				P_plus_4[l][m] = max(0, mu_ad_l_left[l][m] * (u_td_4_next[l][m] - u_td_4_next[l - 1][m])) - 
								 min(0, mu_ad_l_right[l][m] * (u_td_4_next[l + 1][m] - u_td_4_next[l][m])) + 
								 max(0, mu_ad_m_left[l][m] * (u_td_4_next[l][m] - u_td_4_next[l][m - 1])) - 
								 min(0, mu_ad_m_right[l][m] * (u_td_4_next[l][m + 1] - u_td_4_next[l][m]));
				P_plus_5[l][m] = max(0, mu_ad_l_left[l][m] * (u_td_5_next[l][m] - u_td_5_next[l - 1][m])) - 
								 min(0, mu_ad_l_right[l][m] * (u_td_5_next[l + 1][m] - u_td_5_next[l][m])) + 
								 max(0, mu_ad_m_left[l][m] * (u_td_5_next[l][m] - u_td_5_next[l][m - 1])) - 
								 min(0, mu_ad_m_right[l][m] * (u_td_5_next[l][m + 1] - u_td_5_next[l][m]));
				P_plus_6[l][m] = max(0, mu_ad_l_left[l][m] * (u_td_6_next[l][m] - u_td_6_next[l - 1][m])) - 
								 min(0, mu_ad_l_right[l][m] * (u_td_6_next[l + 1][m] - u_td_6_next[l][m])) + 
								 max(0, mu_ad_m_left[l][m] * (u_td_6_next[l][m] - u_td_6_next[l][m - 1])) - 
								 min(0, mu_ad_m_right[l][m] * (u_td_6_next[l][m + 1] - u_td_6_next[l][m]));
				P_plus_7[l][m] = max(0, mu_ad_l_left[l][m] * (u_td_7_next[l][m] - u_td_7_next[l - 1][m])) - 
								 min(0, mu_ad_l_right[l][m] * (u_td_7_next[l + 1][m] - u_td_7_next[l][m])) + 
								 max(0, mu_ad_m_left[l][m] * (u_td_7_next[l][m] - u_td_7_next[l][m - 1])) - 
								 min(0, mu_ad_m_right[l][m] * (u_td_7_next[l][m + 1] - u_td_7_next[l][m]));
				P_plus_8[l][m] = max(0, mu_ad_l_left[l][m] * (u_td_8_next[l][m] - u_td_8_next[l - 1][m])) - 
								 min(0, mu_ad_l_right[l][m] * (u_td_8_next[l + 1][m] - u_td_8_next[l][m])) + 
								 max(0, mu_ad_m_left[l][m] * (u_td_8_next[l][m] - u_td_8_next[l][m - 1])) - 
								 min(0, mu_ad_m_right[l][m] * (u_td_8_next[l][m + 1] - u_td_8_next[l][m]));

				Q_plus_1[l][m] = u_max_1[l][m] - u_td_1_next[l][m];
				Q_plus_2[l][m] = u_max_2[l][m] - u_td_2_next[l][m];
				Q_plus_3[l][m] = u_max_3[l][m] - u_td_3_next[l][m];
				Q_plus_4[l][m] = u_max_4[l][m] - u_td_4_next[l][m];
				Q_plus_5[l][m] = u_max_5[l][m] - u_td_5_next[l][m];
				Q_plus_6[l][m] = u_max_6[l][m] - u_td_6_next[l][m];
				Q_plus_7[l][m] = u_max_7[l][m] - u_td_7_next[l][m];
				Q_plus_8[l][m] = u_max_8[l][m] - u_td_8_next[l][m];

				P_minus_1[l][m] = max(0, mu_ad_l_right[l][m] * (u_td_1_next[l + 1][m] - u_td_1_next[l][m])) - 
								  min(0, mu_ad_l_left[l][m] * (u_td_1_next[l][m] - u_td_1_next[l - 1][m])) + 
								  max(0, mu_ad_m_right[l][m] * (u_td_1_next[l][m + 1] - u_td_1_next[l][m])) - 
								  min(0, mu_ad_m_left[l][m] * (u_td_1_next[l][m] - u_td_1_next[l][m - 1]));
				P_minus_2[l][m] = max(0, mu_ad_l_right[l][m] * (u_td_2_next[l + 1][m] - u_td_2_next[l][m])) - 
								  min(0, mu_ad_l_left[l][m] * (u_td_2_next[l][m] - u_td_2_next[l - 1][m])) + 
								  max(0, mu_ad_m_right[l][m] * (u_td_2_next[l][m + 1] - u_td_2_next[l][m])) - 
								  min(0, mu_ad_m_left[l][m] * (u_td_2_next[l][m] - u_td_2_next[l][m - 1]));
				P_minus_3[l][m] = max(0, mu_ad_l_right[l][m] * (u_td_3_next[l + 1][m] - u_td_3_next[l][m])) - 
								  min(0, mu_ad_l_left[l][m] * (u_td_3_next[l][m] - u_td_3_next[l - 1][m])) + 
								  max(0, mu_ad_m_right[l][m] * (u_td_3_next[l][m + 1] - u_td_3_next[l][m])) - 
								  min(0, mu_ad_m_left[l][m] * (u_td_3_next[l][m] - u_td_3_next[l][m - 1]));
				P_minus_4[l][m] = max(0, mu_ad_l_right[l][m] * (u_td_4_next[l + 1][m] - u_td_4_next[l][m])) - 
								  min(0, mu_ad_l_left[l][m] * (u_td_4_next[l][m] - u_td_4_next[l - 1][m])) + 
								  max(0, mu_ad_m_right[l][m] * (u_td_4_next[l][m + 1] - u_td_4_next[l][m])) - 
								  min(0, mu_ad_m_left[l][m] * (u_td_4_next[l][m] - u_td_4_next[l][m - 1]));
				P_minus_5[l][m] = max(0, mu_ad_l_right[l][m] * (u_td_5_next[l + 1][m] - u_td_5_next[l][m])) - 
								  min(0, mu_ad_l_left[l][m] * (u_td_5_next[l][m] - u_td_5_next[l - 1][m])) + 
								  max(0, mu_ad_m_right[l][m] * (u_td_5_next[l][m + 1] - u_td_5_next[l][m])) - 
								  min(0, mu_ad_m_left[l][m] * (u_td_5_next[l][m] - u_td_5_next[l][m - 1]));
				P_minus_6[l][m] = max(0, mu_ad_l_right[l][m] * (u_td_6_next[l + 1][m] - u_td_6_next[l][m])) - 
								  min(0, mu_ad_l_left[l][m] * (u_td_6_next[l][m] - u_td_6_next[l - 1][m])) + 
								  max(0, mu_ad_m_right[l][m] * (u_td_6_next[l][m + 1] - u_td_6_next[l][m])) - 
								  min(0, mu_ad_m_left[l][m] * (u_td_6_next[l][m] - u_td_6_next[l][m - 1]));
				P_minus_7[l][m] = max(0, mu_ad_l_right[l][m] * (u_td_7_next[l + 1][m] - u_td_7_next[l][m])) - 
								  min(0, mu_ad_l_left[l][m] * (u_td_7_next[l][m] - u_td_7_next[l - 1][m])) + 
								  max(0, mu_ad_m_right[l][m] * (u_td_7_next[l][m + 1] - u_td_7_next[l][m])) - 
								  min(0, mu_ad_m_left[l][m] * (u_td_7_next[l][m] - u_td_7_next[l][m - 1]));
				P_minus_8[l][m] = max(0, mu_ad_l_right[l][m] * (u_td_8_next[l + 1][m] - u_td_8_next[l][m])) - 
								  min(0, mu_ad_l_left[l][m] * (u_td_8_next[l][m] - u_td_8_next[l - 1][m])) + 
								  max(0, mu_ad_m_right[l][m] * (u_td_8_next[l][m + 1] - u_td_8_next[l][m])) - 
								  min(0, mu_ad_m_left[l][m] * (u_td_8_next[l][m] - u_td_8_next[l][m - 1]));

				Q_minus_1[l][m] = u_td_1_next[l][m] - u_min_1[l][m];
				Q_minus_2[l][m] = u_td_2_next[l][m] - u_min_2[l][m];
				Q_minus_3[l][m] = u_td_3_next[l][m] - u_min_3[l][m];
				Q_minus_4[l][m] = u_td_4_next[l][m] - u_min_4[l][m];
				Q_minus_5[l][m] = u_td_5_next[l][m] - u_min_5[l][m];
				Q_minus_6[l][m] = u_td_6_next[l][m] - u_min_6[l][m];
				Q_minus_7[l][m] = u_td_7_next[l][m] - u_min_7[l][m];
				Q_minus_8[l][m] = u_td_8_next[l][m] - u_min_8[l][m];

				if (P_plus_1[l][m] > 0) {
					R_plus_1[l][m] = min(1, Q_plus_1[l][m] / P_plus_1[l][m]);
				} else {
					R_plus_1[l][m] = 0;
				}
				if (P_plus_2[l][m] > 0) {
					R_plus_2[l][m] = min(1, Q_plus_2[l][m] / P_plus_2[l][m]);
				} else {
					R_plus_2[l][m] = 0;
				}
				if (P_plus_3[l][m] > 0) {
					R_plus_3[l][m] = min(1, Q_plus_3[l][m] / P_plus_3[l][m]);
				} else {
					R_plus_3[l][m] = 0;
				}
				if (P_plus_4[l][m] > 0) {
					R_plus_4[l][m] = min(1, Q_plus_4[l][m] / P_plus_4[l][m]);
				} else {
					R_plus_4[l][m] = 0;
				}
				if (P_plus_5[l][m] > 0) {
					R_plus_5[l][m] = min(1, Q_plus_5[l][m] / P_plus_5[l][m]);
				} else {
					R_plus_5[l][m] = 0;
				}
				if (P_plus_6[l][m] > 0) {
					R_plus_6[l][m] = min(1, Q_plus_6[l][m] / P_plus_6[l][m]);
				} else {
					R_plus_6[l][m] = 0;
				}
				if (P_plus_7[l][m] > 0) {
					R_plus_7[l][m] = min(1, Q_plus_7[l][m] / P_plus_7[l][m]);
				} else {
					R_plus_7[l][m] = 0;
				}
				if (P_plus_8[l][m] > 0) {
					R_plus_8[l][m] = min(1, Q_plus_8[l][m] / P_plus_8[l][m]);
				} else {
					R_plus_8[l][m] = 0;
				}

				if (P_minus_1[l][m] > 0) {
					R_minus_1[l][m] = min(1, Q_minus_1[l][m] / P_minus_1[l][m]);
				} else {
					R_minus_1[l][m] = 0;
				}
				if (P_minus_2[l][m] > 0) {
					R_minus_2[l][m] = min(1, Q_minus_2[l][m] / P_minus_2[l][m]);
				} else {
					R_minus_2[l][m] = 0;
				}
				if (P_minus_3[l][m] > 0) {
					R_minus_3[l][m] = min(1, Q_minus_3[l][m] / P_minus_3[l][m]);
				} else {
					R_minus_3[l][m] = 0;
				}
				if (P_minus_4[l][m] > 0) {
					R_minus_4[l][m] = min(1, Q_minus_4[l][m] / P_minus_4[l][m]);
				} else {
					R_minus_4[l][m] = 0;
				}
				if (P_minus_5[l][m] > 0) {
					R_minus_5[l][m] = min(1, Q_minus_5[l][m] / P_minus_5[l][m]);
				} else {
					R_minus_5[l][m] = 0;
				}
				if (P_minus_6[l][m] > 0) {
					R_minus_6[l][m] = min(1, Q_minus_6[l][m] / P_minus_6[l][m]);
				} else {
					R_minus_6[l][m] = 0;
				}
				if (P_minus_7[l][m] > 0) {
					R_minus_7[l][m] = min(1, Q_minus_7[l][m] / P_minus_7[l][m]);
				} else {
					R_minus_7[l][m] = 0;
				}
				if (P_minus_8[l][m] > 0) {
					R_minus_8[l][m] = min(1, Q_minus_8[l][m] / P_minus_8[l][m]);
				} else {
					R_minus_8[l][m] = 0;
				}
			}
		}

		// constants for correction part

		for (int l = 1; l < L_max; l++) {
			for (int m = 1; m < M_max; m++) {
				if (mu_ad_l_left[l][m] * (u_td_1_next[l][m] - u_td_1_next[l - 1][m]) >= 0) {
					COR_l_left_1[l][m] = min(R_plus_1[l][m], R_minus_1[l - 1][m]);
				} else {
					COR_l_left_1[l][m] = min(R_plus_1[l - 1][m], R_minus_1[l][m]);
				}
				if (mu_ad_l_left[l][m] * (u_td_2_next[l][m] - u_td_2_next[l - 1][m]) >= 0) {
					COR_l_left_2[l][m] = min(R_plus_2[l][m], R_minus_2[l - 1][m]);
				} else {
					COR_l_left_2[l][m] = min(R_plus_2[l - 1][m], R_minus_2[l][m]);
				}
				if (mu_ad_l_left[l][m] * (u_td_3_next[l][m] - u_td_3_next[l - 1][m]) >= 0) {
					COR_l_left_3[l][m] = min(R_plus_3[l][m], R_minus_3[l - 1][m]);
				} else {
					COR_l_left_3[l][m] = min(R_plus_3[l - 1][m], R_minus_3[l][m]);
				}
				if (mu_ad_l_left[l][m] * (u_td_4_next[l][m] - u_td_4_next[l - 1][m]) >= 0) {
					COR_l_left_4[l][m] = min(R_plus_4[l][m], R_minus_4[l - 1][m]);
				} else {
					COR_l_left_4[l][m] = min(R_plus_4[l - 1][m], R_minus_4[l][m]);
				}
				if (mu_ad_l_left[l][m] * (u_td_5_next[l][m] - u_td_5_next[l - 1][m]) >= 0) {
					COR_l_left_5[l][m] = min(R_plus_5[l][m], R_minus_5[l - 1][m]);
				} else {
					COR_l_left_5[l][m] = min(R_plus_5[l - 1][m], R_minus_5[l][m]);
				}
				if (mu_ad_l_left[l][m] * (u_td_6_next[l][m] - u_td_6_next[l - 1][m]) >= 0) {
					COR_l_left_6[l][m] = min(R_plus_6[l][m], R_minus_6[l - 1][m]);
				} else {
					COR_l_left_6[l][m] = min(R_plus_6[l - 1][m], R_minus_6[l][m]);
				}
				if (mu_ad_l_left[l][m] * (u_td_7_next[l][m] - u_td_7_next[l - 1][m]) >= 0) {
					COR_l_left_7[l][m] = min(R_plus_7[l][m], R_minus_7[l - 1][m]);
				} else {
					COR_l_left_7[l][m] = min(R_plus_7[l - 1][m], R_minus_7[l][m]);
				}
				if (mu_ad_l_left[l][m] * (u_td_8_next[l][m] - u_td_8_next[l - 1][m]) >= 0) {
					COR_l_left_8[l][m] = min(R_plus_8[l][m], R_minus_8[l - 1][m]);
				} else {
					COR_l_left_8[l][m] = min(R_plus_8[l - 1][m], R_minus_8[l][m]);
				}

				if (mu_ad_l_right[l][m] * (u_td_1_next[l + 1][m] - u_td_1_next[l][m]) >= 0) {
					COR_l_right_1[l][m] = min(R_plus_1[l + 1][m], R_minus_1[l][m]);
				} else {
					COR_l_right_1[l][m] = min(R_plus_1[l][m], R_minus_1[l + 1][m]);
				}
				if (mu_ad_l_right[l][m] * (u_td_2_next[l + 1][m] - u_td_2_next[l][m]) >= 0) {
					COR_l_right_2[l][m] = min(R_plus_2[l + 1][m], R_minus_2[l][m]);
				} else {
					COR_l_right_2[l][m] = min(R_plus_2[l][m], R_minus_2[l + 1][m]);
				}
				if (mu_ad_l_right[l][m] * (u_td_3_next[l + 1][m] - u_td_3_next[l][m]) >= 0) {
					COR_l_right_3[l][m] = min(R_plus_3[l + 1][m], R_minus_3[l][m]);
				} else {
					COR_l_right_3[l][m] = min(R_plus_3[l][m], R_minus_3[l + 1][m]);
				}
				if (mu_ad_l_right[l][m] * (u_td_4_next[l + 1][m] - u_td_4_next[l][m]) >= 0) {
					COR_l_right_4[l][m] = min(R_plus_4[l + 1][m], R_minus_4[l][m]);
				} else {
					COR_l_right_4[l][m] = min(R_plus_4[l][m], R_minus_4[l + 1][m]);
				}
				if (mu_ad_l_right[l][m] * (u_td_5_next[l + 1][m] - u_td_5_next[l][m]) >= 0) {
					COR_l_right_5[l][m] = min(R_plus_5[l + 1][m], R_minus_5[l][m]);
				} else {
					COR_l_right_5[l][m] = min(R_plus_5[l][m], R_minus_5[l + 1][m]);
				}
				if (mu_ad_l_right[l][m] * (u_td_6_next[l + 1][m] - u_td_6_next[l][m]) >= 0) {
					COR_l_right_6[l][m] = min(R_plus_6[l + 1][m], R_minus_6[l][m]);
				} else {
					COR_l_right_6[l][m] = min(R_plus_6[l][m], R_minus_6[l + 1][m]);
				}
				if (mu_ad_l_right[l][m] * (u_td_7_next[l + 1][m] - u_td_7_next[l][m]) >= 0) {
					COR_l_right_7[l][m] = min(R_plus_7[l + 1][m], R_minus_7[l][m]);
				} else {
					COR_l_right_7[l][m] = min(R_plus_7[l][m], R_minus_7[l + 1][m]);
				}
				if (mu_ad_l_right[l][m] * (u_td_8_next[l + 1][m] - u_td_8_next[l][m]) >= 0) {
					COR_l_right_8[l][m] = min(R_plus_8[l + 1][m], R_minus_8[l][m]);
				} else {
					COR_l_right_8[l][m] = min(R_plus_8[l][m], R_minus_8[l + 1][m]);
				}

				if (mu_ad_m_left[l][m] * (u_td_1_next[l][m] - u_td_1_next[l][m - 1]) >= 0) {
					COR_m_left_1[l][m] = min(R_plus_1[l][m], R_minus_1[l][m - 1]);
				} else {
					COR_m_left_1[l][m] = min(R_plus_1[l][m - 1], R_minus_1[l][m]);
				}
				if (mu_ad_m_left[l][m] * (u_td_2_next[l][m] - u_td_2_next[l][m - 1]) >= 0) {
					COR_m_left_2[l][m] = min(R_plus_2[l][m], R_minus_2[l][m - 1]);
				} else {
					COR_m_left_2[l][m] = min(R_plus_2[l][m - 1], R_minus_2[l][m]);
				}
				if (mu_ad_m_left[l][m] * (u_td_3_next[l][m] - u_td_3_next[l][m - 1]) >= 0) {
					COR_m_left_3[l][m] = min(R_plus_3[l][m], R_minus_3[l][m - 1]);
				} else {
					COR_m_left_3[l][m] = min(R_plus_3[l][m - 1], R_minus_3[l][m]);
				}
				if (mu_ad_m_left[l][m] * (u_td_4_next[l][m] - u_td_4_next[l][m - 1]) >= 0) {
					COR_m_left_4[l][m] = min(R_plus_4[l][m], R_minus_4[l][m - 1]);
				} else {
					COR_m_left_4[l][m] = min(R_plus_4[l][m - 1], R_minus_4[l][m]);
				}
				if (mu_ad_m_left[l][m] * (u_td_5_next[l][m] - u_td_5_next[l][m - 1]) >= 0) {
					COR_m_left_5[l][m] = min(R_plus_5[l][m], R_minus_5[l][m - 1]);
				} else {
					COR_m_left_5[l][m] = min(R_plus_5[l][m - 1], R_minus_5[l][m]);
				}
				if (mu_ad_m_left[l][m] * (u_td_6_next[l][m] - u_td_6_next[l][m - 1]) >= 0) {
					COR_m_left_6[l][m] = min(R_plus_6[l][m], R_minus_6[l][m - 1]);
				} else {
					COR_m_left_6[l][m] = min(R_plus_6[l][m - 1], R_minus_6[l][m]);
				}
				if (mu_ad_m_left[l][m] * (u_td_7_next[l][m] - u_td_7_next[l][m - 1]) >= 0) {
					COR_m_left_7[l][m] = min(R_plus_7[l][m], R_minus_7[l][m - 1]);
				} else {
					COR_m_left_7[l][m] = min(R_plus_7[l][m - 1], R_minus_7[l][m]);
				}
				if (mu_ad_m_left[l][m] * (u_td_8_next[l][m] - u_td_8_next[l][m - 1]) >= 0) {
					COR_m_left_8[l][m] = min(R_plus_8[l][m], R_minus_8[l][m - 1]);
				} else {
					COR_m_left_8[l][m] = min(R_plus_8[l][m - 1], R_minus_8[l][m]);
				}

				if (mu_ad_m_right[l][m] * (u_td_1_next[l][m + 1] - u_td_1_next[l][m]) >= 0) {
					COR_m_right_1[l][m] = min(R_plus_1[l][m + 1], R_minus_1[l][m]);
				} else {
					COR_m_right_1[l][m] = min(R_plus_1[l][m], R_minus_1[l][m + 1]);
				}
				if (mu_ad_m_right[l][m] * (u_td_2_next[l][m + 1] - u_td_2_next[l][m]) >= 0) {
					COR_m_right_2[l][m] = min(R_plus_2[l][m + 1], R_minus_2[l][m]);
				} else {
					COR_m_right_2[l][m] = min(R_plus_2[l][m], R_minus_2[l][m + 1]);
				}
				if (mu_ad_m_right[l][m] * (u_td_3_next[l][m + 1] - u_td_3_next[l][m]) >= 0) {
					COR_m_right_3[l][m] = min(R_plus_3[l][m + 1], R_minus_3[l][m]);
				} else {
					COR_m_right_3[l][m] = min(R_plus_3[l][m], R_minus_3[l][m + 1]);
				}
				if (mu_ad_m_right[l][m] * (u_td_4_next[l][m + 1] - u_td_4_next[l][m]) >= 0) {
					COR_m_right_4[l][m] = min(R_plus_4[l][m + 1], R_minus_4[l][m]);
				} else {
					COR_m_right_4[l][m] = min(R_plus_4[l][m], R_minus_4[l][m + 1]);
				}
				if (mu_ad_m_right[l][m] * (u_td_5_next[l][m + 1] - u_td_5_next[l][m]) >= 0) {
					COR_m_right_5[l][m] = min(R_plus_5[l][m + 1], R_minus_5[l][m]);
				} else {
					COR_m_right_5[l][m] = min(R_plus_5[l][m], R_minus_5[l][m + 1]);
				}
				if (mu_ad_m_right[l][m] * (u_td_6_next[l][m + 1] - u_td_6_next[l][m]) >= 0) {
					COR_m_right_6[l][m] = min(R_plus_6[l][m + 1], R_minus_6[l][m]);
				} else {
					COR_m_right_6[l][m] = min(R_plus_6[l][m], R_minus_6[l][m + 1]);
				}
				if (mu_ad_m_right[l][m] * (u_td_7_next[l][m + 1] - u_td_7_next[l][m]) >= 0) {
					COR_m_right_7[l][m] = min(R_plus_7[l][m + 1], R_minus_7[l][m]);
				} else {
					COR_m_right_7[l][m] = min(R_plus_7[l][m], R_minus_7[l][m + 1]);
				}
				if (mu_ad_m_right[l][m] * (u_td_8_next[l][m + 1] - u_td_8_next[l][m]) >= 0) {
					COR_m_right_8[l][m] = min(R_plus_8[l][m + 1], R_minus_8[l][m]);
				} else {
					COR_m_right_8[l][m] = min(R_plus_8[l][m], R_minus_8[l][m + 1]);
				}
			}
		}

		// filling central points of grid for correction part

		for (int l = 1; l < L_max; l++) {
			for (int m = 1; m < M_max; m++) {
				u_cor_1[l][m] = u_td_1_next[l][m] - 
								dt / dz * (COR_l_right_1[l][m] * mu_ad_l_right[l][m] * (u_td_1_next[l + 1][m] - u_td_1_next[l][m]) - 
										   COR_l_left_1[l][m] * mu_ad_l_left[l][m] * (u_td_1_next[l][m] - u_td_1_next[l - 1][m])) - 
								dt / dr[l] * (COR_m_right_1[l][m] * mu_ad_m_right[l][m] * (u_td_1_next[l][m + 1] - u_td_1_next[l][m]) - 
											  COR_m_left_1[l][m] * mu_ad_m_left[l][m] * (u_td_1_next[l][m] - u_td_1_next[l][m - 1]));
				u_cor_2[l][m] = u_td_2_next[l][m] - 
								dt / dz * (COR_l_right_2[l][m] * mu_ad_l_right[l][m] * (u_td_2_next[l + 1][m] - u_td_2_next[l][m]) - 
										   COR_l_left_2[l][m] * mu_ad_l_left[l][m] * (u_td_2_next[l][m] - u_td_2_next[l - 1][m])) - 
								dt / dr[l] * (COR_m_right_2[l][m] * mu_ad_m_right[l][m] * (u_td_2_next[l][m + 1] - u_td_2_next[l][m]) - 
											  COR_m_left_2[l][m] * mu_ad_m_left[l][m] * (u_td_2_next[l][m] - u_td_2_next[l][m - 1]));
				u_cor_3[l][m] = u_td_3_next[l][m] - 
								dt / dz * (COR_l_right_3[l][m] * mu_ad_l_right[l][m] * (u_td_3_next[l + 1][m] - u_td_3_next[l][m]) - 
										   COR_l_left_3[l][m] * mu_ad_l_left[l][m] * (u_td_3_next[l][m] - u_td_3_next[l - 1][m])) - 
								dt / dr[l] * (COR_m_right_3[l][m] * mu_ad_m_right[l][m] * (u_td_3_next[l][m + 1] - u_td_3_next[l][m]) - 
											  COR_m_left_3[l][m] * mu_ad_m_left[l][m] * (u_td_3_next[l][m] - u_td_3_next[l][m - 1]));
				u_cor_4[l][m] = u_td_4_next[l][m] - 
								dt / dz * (COR_l_right_4[l][m] * mu_ad_l_right[l][m] * (u_td_4_next[l + 1][m] - u_td_4_next[l][m]) - 
										   COR_l_left_4[l][m] * mu_ad_l_left[l][m] * (u_td_4_next[l][m] - u_td_4_next[l - 1][m])) - 
								dt / dr[l] * (COR_m_right_4[l][m] * mu_ad_m_right[l][m] * (u_td_4_next[l][m + 1] - u_td_4_next[l][m]) - 
											  COR_m_left_4[l][m] * mu_ad_m_left[l][m] * (u_td_4_next[l][m] - u_td_4_next[l][m - 1]));
				u_cor_5[l][m] = u_td_5_next[l][m] - 
								dt / dz * (COR_l_right_5[l][m] * mu_ad_l_right[l][m] * (u_td_5_next[l + 1][m] - u_td_5_next[l][m]) - 
										   COR_l_left_5[l][m] * mu_ad_l_left[l][m] * (u_td_5_next[l][m] - u_td_5_next[l - 1][m])) - 
								dt / dr[l] * (COR_m_right_5[l][m] * mu_ad_m_right[l][m] * (u_td_5_next[l][m + 1] - u_td_5_next[l][m]) - 
											  COR_m_left_5[l][m] * mu_ad_m_left[l][m] * (u_td_5_next[l][m] - u_td_5_next[l][m - 1]));
				u_cor_6[l][m] = u_td_6_next[l][m] - 
								dt / dz * (COR_l_right_6[l][m] * mu_ad_l_right[l][m] * (u_td_6_next[l + 1][m] - u_td_6_next[l][m]) - 
										   COR_l_left_6[l][m] * mu_ad_l_left[l][m] * (u_td_6_next[l][m] - u_td_6_next[l - 1][m])) - 
								dt / dr[l] * (COR_m_right_6[l][m] * mu_ad_m_right[l][m] * (u_td_6_next[l][m + 1] - u_td_6_next[l][m]) - 
											  COR_m_left_6[l][m] * mu_ad_m_left[l][m] * (u_td_6_next[l][m] - u_td_6_next[l][m - 1]));
				u_cor_7[l][m] = u_td_7_next[l][m] - 
								dt / dz * (COR_l_right_7[l][m] * mu_ad_l_right[l][m] * (u_td_7_next[l + 1][m] - u_td_7_next[l][m]) - 
										   COR_l_left_7[l][m] * mu_ad_l_left[l][m] * (u_td_7_next[l][m] - u_td_7_next[l - 1][m])) - 
								dt / dr[l] * (COR_m_right_7[l][m] * mu_ad_m_right[l][m] * (u_td_7_next[l][m + 1] - u_td_7_next[l][m]) - 
											  COR_m_left_7[l][m] * mu_ad_m_left[l][m] * (u_td_7_next[l][m] - u_td_7_next[l][m - 1]));
				u_cor_8[l][m] = u_td_8_next[l][m] - 
								dt / dz * (COR_l_right_8[l][m] * mu_ad_l_right[l][m] * (u_td_8_next[l + 1][m] - u_td_8_next[l][m]) - 
										   COR_l_left_8[l][m] * mu_ad_l_left[l][m] * (u_td_8_next[l][m] - u_td_8_next[l - 1][m])) - 
								dt / dr[l] * (COR_m_right_8[l][m] * mu_ad_m_right[l][m] * (u_td_8_next[l][m + 1] - u_td_8_next[l][m]) - 
											  COR_m_left_8[l][m] * mu_ad_m_left[l][m] * (u_td_8_next[l][m] - u_td_8_next[l][m - 1]));
			}
		}

		// left boundary condition for current layer for correction part

		for (int m = 0; m < M_max + 1; m++) {
			u_cor_1[0][m] = rho[0][m] * r[0][m];
			u_cor_2[0][m] = u_cor_2[1][m];
			// u_cor_2[0][m] = 2 * u_cor_2[1][m] - u_cor_2[2][m];
			u_cor_3[0][m] = u_cor_2[0][m] * r_z[0][m];
			u_cor_4[0][m] = rho[0][m] * v_phi[0][m] * r[0][m];
			u_cor_5[0][m] = rho[0][m] * e[0][m] * r[0][m];
			u_cor_6[0][m] = H_phi[0][m];
			u_cor_7[0][m] = H_z[0][m] * r[0][m];
			u_cor_8[0][m] = u_cor_7[0][m] * r_z[0][m];

			/*u_cor_1[0][m] = u_td_1_next[0][m];
			u_cor_2[0][m] = u_td_2_next[0][m];
			u_cor_3[0][m] = u_td_3_next[0][m];
			u_cor_4[0][m] = u_td_4_next[0][m];
			u_cor_5[0][m] = u_td_5_next[0][m];
			u_cor_6[0][m] = u_td_6_next[0][m];
			u_cor_7[0][m] = u_td_7_next[0][m];
			u_cor_8[0][m] = u_td_8_next[0][m];*/
		}

		// up and down boundary condition for current layer for correction part

		for (int l = 1; l < L_max + 1; l++) {
			u_cor_1[l][0] = u_cor_1[l][1];
			u_cor_2[l][0] = u_cor_2[l][1];
			u_cor_3[l][0] = u_cor_2[l][0] * r_z[l][0];
			u_cor_4[l][0] = u_cor_4[l][1];
			u_cor_5[l][0] = u_cor_5[l][1];
			u_cor_6[l][0] = u_cor_6[l][1];
			u_cor_7[l][0] = u_cor_7[l][1];
			u_cor_8[l][0] = u_cor_7[l][0] * r_z[l][0];

			u_cor_1[l][M_max] = u_cor_1[l][M_max - 1];
			u_cor_2[l][M_max] = u_cor_2[l][M_max - 1];
			u_cor_3[l][M_max] = u_cor_2[l][M_max] * r_z[l][M_max];
			u_cor_4[l][M_max] = u_cor_4[l][M_max - 1];
			u_cor_5[l][M_max] = u_cor_5[l][M_max - 1];
			u_cor_6[l][M_max] = u_cor_6[l][M_max - 1];
			u_cor_7[l][M_max] = u_cor_7[l][M_max - 1];
			u_cor_8[l][M_max] = u_cor_7[l][M_max] * r_z[l][M_max];

			/*u_cor_1[l][0] = u_td_1_next[l][0];
			u_cor_2[l][0] = u_td_2_next[l][0];
			u_cor_3[l][0] = u_td_3_next[l][0];
			u_cor_4[l][0] = u_td_4_next[l][0];
			u_cor_5[l][0] = u_td_5_next[l][0];
			u_cor_6[l][0] = u_td_6_next[l][0];
			u_cor_7[l][0] = u_td_7_next[l][0];
			u_cor_8[l][0] = u_td_8_next[l][0];

			u_cor_1[l][M_max] = u_td_1_next[l][M_max];
			u_cor_2[l][M_max] = u_td_2_next[l][M_max];
			u_cor_3[l][M_max] = u_td_3_next[l][M_max];
			u_cor_4[l][M_max] = u_td_4_next[l][M_max];
			u_cor_5[l][M_max] = u_td_5_next[l][M_max];
			u_cor_6[l][M_max] = u_td_6_next[l][M_max];
			u_cor_7[l][M_max] = u_td_7_next[l][M_max];
			u_cor_8[l][M_max] = u_td_8_next[l][M_max];*/
		}

		// right boundary condition for current layer for anti-diffusion part

		for (int m = 1; m < M_max; m++) {
			u_cor_1[L_max][m] = u_cor_1[L_max - 1][m];
			u_cor_2[L_max][m] = u_cor_2[L_max - 1][m];
			u_cor_3[L_max][m] = u_cor_3[L_max - 1][m];
			u_cor_4[L_max][m] = u_cor_4[L_max - 1][m];
			u_cor_5[L_max][m] = u_cor_5[L_max - 1][m];
			u_cor_6[L_max][m] = u_cor_6[L_max - 1][m];
			u_cor_7[L_max][m] = u_cor_7[L_max - 1][m];
			u_cor_8[L_max][m] = u_cor_8[L_max - 1][m];

			/*u_cor_1[L_max][m] = u_td_1_next[L_max][m];
			u_cor_2[L_max][m] = u_td_2_next[L_max][m];
			u_cor_3[L_max][m] = u_td_3_next[L_max][m];
			u_cor_4[L_max][m] = u_td_4_next[L_max][m];
			u_cor_5[L_max][m] = u_td_5_next[L_max][m];
			u_cor_6[L_max][m] = u_td_6_next[L_max][m];
			u_cor_7[L_max][m] = u_td_7_next[L_max][m];
			u_cor_8[L_max][m] = u_td_8_next[L_max][m];*/
		}

		// update parameters of problem

		for (int l = 0; l < L_max + 1; l++) {
			for (int m = 0; m < M_max + 1; m++) {
				rho[l][m] = u_cor_1[l][m] / r[l][m];
				v_z[l][m] = u_cor_2[l][m] / u_cor_1[l][m];
				v_r[l][m] = u_cor_3[l][m] / u_cor_1[l][m];
				v_phi[l][m] = u_cor_4[l][m] / u_cor_1[l][m];
				v_y[l][m] = v_r[l][m] - v_z[l][m] * r_z[l][m];

				H_phi[l][m] = u_cor_6[l][m];
				H_z[l][m] = u_cor_7[l][m] / r[l][m];
				H_r[l][m] = u_cor_8[l][m] / r[l][m];
				H_y[l][m] = H_r[l][m] - H_z[l][m] * r_z[l][m];

				e[l][m] = u_cor_5[l][m] / u_cor_1[l][m];
				p[l][m] = (gamma - 1) * rho[l][m] * e[l][m];
				P[l][m] = p[l][m] + 1.0 / 2.0 * (pow(H_z[l][m], 2) + pow(H_r[l][m], 2) + pow(H_phi[l][m], 2));
			}
		}

		// update u_td

		for (int l = 0; l < L_max + 1; l++) {
			for (int m = 0; m < M_max + 1; m++) {
				u_td_1[l][m] = u_cor_1[l][m];
				u_td_2[l][m] = u_cor_2[l][m];
				u_td_3[l][m] = u_cor_3[l][m];
				u_td_4[l][m] = u_cor_4[l][m];
				u_td_5[l][m] = u_cor_5[l][m];
				u_td_6[l][m] = u_cor_6[l][m];
				u_td_7[l][m] = u_cor_7[l][m];
				u_td_8[l][m] = u_cor_8[l][m];
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

	memory_clearing(u_ad_1, L_max);
	memory_clearing(u_ad_2, L_max);
	memory_clearing(u_ad_3, L_max);
	memory_clearing(u_ad_4, L_max);
	memory_clearing(u_ad_5, L_max);
	memory_clearing(u_ad_6, L_max);
	memory_clearing(u_ad_7, L_max);
	memory_clearing(u_ad_8, L_max);

	memory_clearing(u_a_1, L_max);
	memory_clearing(u_a_2, L_max);
	memory_clearing(u_a_3, L_max);
	memory_clearing(u_a_4, L_max);
	memory_clearing(u_a_5, L_max);
	memory_clearing(u_a_6, L_max);
	memory_clearing(u_a_7, L_max);
	memory_clearing(u_a_8, L_max);

	memory_clearing(u_b_1, L_max);
	memory_clearing(u_b_2, L_max);
	memory_clearing(u_b_3, L_max);
	memory_clearing(u_b_4, L_max);
	memory_clearing(u_b_5, L_max);
	memory_clearing(u_b_6, L_max);
	memory_clearing(u_b_7, L_max);
	memory_clearing(u_b_8, L_max);

	memory_clearing(u_max_1, L_max);
	memory_clearing(u_max_2, L_max);
	memory_clearing(u_max_3, L_max);
	memory_clearing(u_max_4, L_max);
	memory_clearing(u_max_5, L_max);
	memory_clearing(u_max_6, L_max);
	memory_clearing(u_max_7, L_max);
	memory_clearing(u_max_8, L_max);

	memory_clearing(u_min_1, L_max);
	memory_clearing(u_min_2, L_max);
	memory_clearing(u_min_3, L_max);
	memory_clearing(u_min_4, L_max);
	memory_clearing(u_min_5, L_max);
	memory_clearing(u_min_6, L_max);
	memory_clearing(u_min_7, L_max);
	memory_clearing(u_min_8, L_max);

	memory_clearing(u_cor_1, L_max);
	memory_clearing(u_cor_2, L_max);
	memory_clearing(u_cor_3, L_max);
	memory_clearing(u_cor_4, L_max);
	memory_clearing(u_cor_5, L_max);
	memory_clearing(u_cor_6, L_max);
	memory_clearing(u_cor_7, L_max);
	memory_clearing(u_cor_8, L_max);
	
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

	memory_clearing(mu_ad_l_left, L_max);
	memory_clearing(mu_ad_l_right, L_max);
	memory_clearing(mu_ad_m_left, L_max);
	memory_clearing(mu_ad_m_right, L_max);

	memory_clearing(P_plus_1, L_max);
	memory_clearing(P_plus_2, L_max);
	memory_clearing(P_plus_3, L_max);
	memory_clearing(P_plus_4, L_max);
	memory_clearing(P_plus_5, L_max);
	memory_clearing(P_plus_6, L_max);
	memory_clearing(P_plus_7, L_max);
	memory_clearing(P_plus_8, L_max);

	memory_clearing(P_minus_1, L_max);
	memory_clearing(P_minus_2, L_max);
	memory_clearing(P_minus_3, L_max);
	memory_clearing(P_minus_4, L_max);
	memory_clearing(P_minus_5, L_max);
	memory_clearing(P_minus_6, L_max);
	memory_clearing(P_minus_7, L_max);
	memory_clearing(P_minus_8, L_max);

	memory_clearing(Q_plus_1, L_max);
	memory_clearing(Q_plus_2, L_max);
	memory_clearing(Q_plus_3, L_max);
	memory_clearing(Q_plus_4, L_max);
	memory_clearing(Q_plus_5, L_max);
	memory_clearing(Q_plus_6, L_max);
	memory_clearing(Q_plus_7, L_max);
	memory_clearing(Q_plus_8, L_max);

	memory_clearing(Q_minus_1, L_max);
	memory_clearing(Q_minus_2, L_max);
	memory_clearing(Q_minus_3, L_max);
	memory_clearing(Q_minus_4, L_max);
	memory_clearing(Q_minus_5, L_max);
	memory_clearing(Q_minus_6, L_max);
	memory_clearing(Q_minus_7, L_max);
	memory_clearing(Q_minus_8, L_max);

	memory_clearing(R_plus_1, L_max);
	memory_clearing(R_plus_2, L_max);
	memory_clearing(R_plus_3, L_max);
	memory_clearing(R_plus_4, L_max);
	memory_clearing(R_plus_5, L_max);
	memory_clearing(R_plus_6, L_max);
	memory_clearing(R_plus_7, L_max);
	memory_clearing(R_plus_8, L_max);

	memory_clearing(R_minus_1, L_max);
	memory_clearing(R_minus_2, L_max);
	memory_clearing(R_minus_3, L_max);
	memory_clearing(R_minus_4, L_max);
	memory_clearing(R_minus_5, L_max);
	memory_clearing(R_minus_6, L_max);
	memory_clearing(R_minus_7, L_max);
	memory_clearing(R_minus_8, L_max);

	memory_clearing(COR_l_left_1, L_max);
	memory_clearing(COR_l_left_2, L_max);
	memory_clearing(COR_l_left_3, L_max);
	memory_clearing(COR_l_left_4, L_max);
	memory_clearing(COR_l_left_5, L_max);
	memory_clearing(COR_l_left_6, L_max);
	memory_clearing(COR_l_left_7, L_max);
	memory_clearing(COR_l_left_8, L_max);

	memory_clearing(COR_l_right_1, L_max);
	memory_clearing(COR_l_right_2, L_max);
	memory_clearing(COR_l_right_3, L_max);
	memory_clearing(COR_l_right_4, L_max);
	memory_clearing(COR_l_right_5, L_max);
	memory_clearing(COR_l_right_6, L_max);
	memory_clearing(COR_l_right_7, L_max);
	memory_clearing(COR_l_right_8, L_max);

	memory_clearing(COR_m_left_1, L_max);
	memory_clearing(COR_m_left_2, L_max);
	memory_clearing(COR_m_left_3, L_max);
	memory_clearing(COR_m_left_4, L_max);
	memory_clearing(COR_m_left_5, L_max);
	memory_clearing(COR_m_left_6, L_max);
	memory_clearing(COR_m_left_7, L_max);
	memory_clearing(COR_m_left_8, L_max);

	memory_clearing(COR_m_right_1, L_max);
	memory_clearing(COR_m_right_2, L_max);
	memory_clearing(COR_m_right_3, L_max);
	memory_clearing(COR_m_right_4, L_max);
	memory_clearing(COR_m_right_5, L_max);
	memory_clearing(COR_m_right_6, L_max);
	memory_clearing(COR_m_right_7, L_max);
	memory_clearing(COR_m_right_8, L_max);
	
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
