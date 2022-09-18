#include <iostream>
#include <fstream>
#include <cmath>

int round(double d);
double max(double *list, double len);

int main() {
	// mesh options

	int I_max = 100;
	int N_max = 1000;

	// u = {rho*S, rho*v*S, rho*energy*S, H_phi*S}
	// g(u) = {0, -S*d/dz(p + H_phi^2 / 2), -p*d/dz(v*S), 0}

	// u_mesh_1 for rho*S, u_mesh_2 for rho*v*S, u_mesh_3 for rho*energy*S, u_mesh_4 for H_phi*S

	double **u_mesh_1 = (double **)malloc((N_max + 1) * sizeof(double *));
	for (int i = 0; i < N_max + 1; i++) {
		u_mesh_1[i] = (double *)malloc((I_max + 1) * sizeof(double));
	}

	double **u_mesh_2 = (double **)malloc((N_max + 1) * sizeof(double *));
	for (int i = 0; i < N_max + 1; i++) {
		u_mesh_2[i] = (double *)malloc((I_max + 1) * sizeof(double));
	}

	double **u_mesh_3 = (double **)malloc((N_max + 1) * sizeof(double *));
	for (int i = 0; i < N_max + 1; i++) {
		u_mesh_3[i] = (double *)malloc((I_max + 1) * sizeof(double));
	}

	double **u_mesh_4 = (double **)malloc((N_max + 1) * sizeof(double *));
	for (int i = 0; i < N_max + 1; i++) {
		u_mesh_4[i] = (double *)malloc((I_max + 1) * sizeof(double));
	}

	// constants
	
	double k = 0.5;
	double gamma = 1.67;
	double beta = 1.0;
	double mu_0 = 0.7;

	// constants on current layer

	double *mu_left = (double *)malloc((I_max + 1) * sizeof(double));
	double *mu_right = (double *)malloc((I_max + 1) * sizeof(double));
	double *C_left = (double *)malloc((I_max + 1) * sizeof(double));
	double *C_0 = (double *)malloc((I_max + 1) * sizeof(double));
	double *C_right = (double *)malloc((I_max + 1) * sizeof(double));

	// create meshes for parameters of problem

	double *rho_mesh = (double *)malloc((I_max + 1) * sizeof(double));
	double *v_mesh = (double *)malloc((I_max + 1) * sizeof(double));
	double *e_mesh = (double *)malloc((I_max + 1) * sizeof(double));
	double *p_mesh = (double *)malloc((I_max + 1) * sizeof(double));
	double *H_phi_mesh = (double *)malloc((I_max + 1) * sizeof(double));
	double *S_mesh = (double *)malloc((I_max + 1) * sizeof(double));

	// initial condition

	double rho_0 = 1.0;
	double v_0 = 0.1;
	double e_0 = beta / (2.0 * (gamma - 1.0));
	double p_0 = beta / 2.0 * pow(rho_0, gamma);

	// grid steps

	double h = 1.0 / I_max;
	double tau = k * h / (2.0 * mu_0 + v_0);

	// filling meshes for parameters

	for (int i = 0; i < I_max + 1; i++) {
		rho_mesh[i] = rho_0;
		v_mesh[i] = v_0;
		e_mesh[i] = e_0;
		p_mesh[i] = p_0;
		H_phi_mesh[i] = 1.0 - 0.9 * i * h;
		S_mesh[i] = 2.0 * pow((i * h - 0.5), 2) + 0.5;
	}

	// filling meshes for u

	for (int i = 0; i < I_max + 1; i++) {
		u_mesh_1[0][i] = rho_mesh[i] * S_mesh[i];
		u_mesh_2[0][i] = rho_mesh[i] * v_mesh[i] * S_mesh[i];
		u_mesh_3[0][i] = rho_mesh[i] * e_mesh[i] * S_mesh[i];
		u_mesh_4[0][i] = H_phi_mesh[i] * S_mesh[i];
	}

	// evolution in time
	
	for (int n = 0; n < N_max; n++) {

		// left boundary condition for current layer

		rho_mesh[0] = 1.0;
		v_mesh[0] = u_mesh_2[n][1] / (rho_mesh[0] * S_mesh[0]);
		e_mesh[0] = beta / (2.0 * (gamma - 1.0));
		H_phi_mesh[0] = 1.0;
		
		u_mesh_1[n + 1][0] = rho_mesh[0] * S_mesh[0];
		u_mesh_2[n + 1][0] = rho_mesh[0] * v_mesh[0] * S_mesh[0];
		u_mesh_3[n + 1][0] = rho_mesh[0] * e_mesh[0] * S_mesh[0];
		u_mesh_4[n + 1][0] = H_phi_mesh[0] * S_mesh[0];

		// constants for current layer

		for (int i = 1; i < I_max; i++) {
			mu_left[i] = mu_0 + (abs(v_mesh[i - 1]) + abs(v_mesh[i])) / 4.0;
			mu_right[i] = mu_0 + (abs(v_mesh[i]) + abs(v_mesh[i + 1])) / 4.0;
			C_left[i] = tau / h * (v_mesh[i - 1] + v_mesh[i]) / 4.0 + tau / h * mu_left[i];
			C_0[i] = 1 - tau / h * (v_mesh[i + 1] - v_mesh[i - 1]) / 4.0 - tau / h * (mu_left[i] + mu_right[i]);
			C_right[i] = - tau / h * (v_mesh[i] + v_mesh[i + 1]) / 4.0 + tau / h * mu_right[i];
		}

		// filling central points of grid

		for (int i = 1; i < I_max; i++) {
			u_mesh_1[n + 1][i] = C_left[i] * u_mesh_1[n][i - 1] + C_0[i] * u_mesh_1[n][i] + C_right[i] * u_mesh_1[n][i + 1];

			u_mesh_2[n + 1][i] = C_left[i] * u_mesh_2[n][i - 1] + C_0[i] * u_mesh_2[n][i] + 
				C_right[i] * u_mesh_2[n][i + 1] - tau / (2.0 * h) * (p_mesh[i + 1] - p_mesh[i - 1] + 
				pow(H_phi_mesh[i + 1], 2) / 2.0 - pow(H_phi_mesh[i - 1], 2) / 2.0) * S_mesh[i];

	        u_mesh_3[n + 1][i] = C_left[i] * u_mesh_3[n][i - 1] + C_0[i] * u_mesh_3[n][i] + 
				C_right[i] * u_mesh_3[n][i + 1] - tau / (2.0 * h) * p_mesh[i] * (v_mesh[i + 1] * S_mesh[i + 1] - 
				v_mesh[i - 1] * S_mesh[i - 1]);

			u_mesh_4[n + 1][i] = C_left[i] * u_mesh_4[n][i - 1] + C_0[i] * u_mesh_4[n][i] + C_right[i] * u_mesh_4[n][i + 1];
		}

		// right boundary condition for current layer

		u_mesh_1[n + 1][I_max] = u_mesh_1[n + 1][I_max - 1];
		u_mesh_2[n + 1][I_max] = u_mesh_2[n + 1][I_max - 1];
		u_mesh_3[n + 1][I_max] = u_mesh_3[n + 1][I_max - 1];
		u_mesh_4[n + 1][I_max] = u_mesh_4[n + 1][I_max - 1];

		// update parameters of problem

		for (int i = 0; i < I_max + 1; i++) {
			rho_mesh[i] = u_mesh_1[n + 1][i] / S_mesh[i];
			v_mesh[i] = u_mesh_2[n + 1][i] / u_mesh_1[n + 1][i];
			e_mesh[i] = u_mesh_3[n + 1][i] / u_mesh_1[n + 1][i];
			p_mesh[i] = beta / 2 * pow(rho_mesh[i], gamma);
			H_phi_mesh[i] = u_mesh_4[n + 1][i] / S_mesh[i];
		}

		// update tau

		tau = k * h / (2 * mu_0 + max(v_mesh, I_max + 1));
	}

	// checkout

	double *x_lst = (double *)malloc((I_max + 1) * sizeof(double));
	double *rho_lst = (double *)malloc((I_max + 1) * sizeof(double));
	double *v_lst = (double *)malloc((I_max + 1) * sizeof(double));
	double *e_lst = (double *)malloc((I_max + 1) * sizeof(double));
	double *H_lst = (double *)malloc((I_max + 1) * sizeof(double));

	x_lst[0] = 0;
	for (int i = 1; i < I_max + 1; i++) {
		x_lst[i] = x_lst[i - 1] + h;
	}

	// enter the time layer for which you want to display the values

	int n = N_max;

	// for rho

	for (int i = 0; i < I_max + 1; i++) {
		rho_lst[i] = u_mesh_1[n][i] / S_mesh[i];
	}

	// for velocity

	for (int i = 0; i < I_max + 1; i++) {
		v_lst[i] = u_mesh_2[n][i] / u_mesh_1[n][i];
	}

	// for energy

	for (int i = 0; i < I_max + 1; i++) {
		e_lst[i] = u_mesh_3[n][i] / u_mesh_1[n][i];
	}

	// for strength

	for (int i = 0; i < I_max + 1; i++) {
		H_lst[i] = u_mesh_4[n][i] / S_mesh[i];
	}

	// output results in file

	std :: ofstream f("data_MGD_CPP.txt");
	for (int i = 0; i < I_max + 1; i++) {
		f << x_lst[i] << ", " << rho_lst[i] << ", " << v_lst[i] << ", " << e_lst[i] << ", " << H_lst[i] << "\n";
	}
	f.close();

	// free memory
	for (int i = 0; i < N_max + 1; i++) {
		delete [] u_mesh_1[i];
	}
	delete [] u_mesh_1;
	for (int i = 0; i < N_max + 1; i++) {
		delete [] u_mesh_2[i];
	}
	delete [] u_mesh_2;
	for (int i = 0; i < N_max + 1; i++) {
		delete [] u_mesh_3[i];
	}
	delete [] u_mesh_3;
	for (int i = 0; i < N_max + 1; i++) {
		delete [] u_mesh_4[i];
	}
	delete [] u_mesh_4;

	delete [] rho_mesh;
	delete [] v_mesh;
	delete [] e_mesh;
	delete [] p_mesh;
	delete [] H_phi_mesh;
	delete [] S_mesh;

	delete [] x_lst;
	delete [] rho_lst;
	delete [] v_lst;
	delete [] e_lst;
	delete [] H_lst;
	
	return 0;
}

int round(double d) {
	return (int)floor(d + 0.5);
}

double max(double *list, double len) {
	double maxim = 0.0;

	for (int i = 0; i < len; i++) {
		if (maxim < list[i]) {
			maxim = list[i];
		}
	}

	return maxim;
}
