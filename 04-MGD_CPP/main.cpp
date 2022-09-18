#include <iostream>
#include <fstream>
#include <cmath>

int round(double d);

int main() {

	// mesh options
	double time = 1.0;
	int I_max = 50;		// there will be (I_max + 1) values in the mesh (enter the node number starting from 0)
	double tau = 0.01;
	double h = 1.0 / (double) I_max;
	int N_max = round(time / tau);
	std :: cout << N_max;

	// u_mesh_1 for rho*S, u_mesh_2 for rho*v*S, u_mesh_3 for rho*energy*S, u_mesh_4 for H_phi*S
	// u = {rho*S, rho*v*S, rho*energy*S, H_phi*S}

	double **u_mesh_1 = (double **)malloc((N_max + 1) * sizeof(double *));
	for (int i = 0; i < N_max + 1; i++) {
		u_mesh_1[i] = (double *)malloc((I_max + 1) * sizeof(double));
		for (int j = 0; j < I_max + 1; j++) {
			u_mesh_1[i][j] = 0;
		}
	}

	double **u_mesh_2 = (double **)malloc((N_max + 1) * sizeof(double *));
	for (int i = 0; i < N_max + 1; i++) {
		u_mesh_2[i] = (double *)malloc((I_max + 1) * sizeof(double));
		for (int j = 0; j < I_max + 1; j++) {
			u_mesh_2[i][j] = 0;
		}
	}

	double **u_mesh_3 = (double **)malloc((N_max + 1) * sizeof(double *));
	for (int i = 0; i < N_max + 1; i++) {
		u_mesh_3[i] = (double *)malloc((I_max + 1) * sizeof(double));
		for (int j = 0; j < I_max + 1; j++) {
			u_mesh_3[i][j] = 0;
		}
	}

	double **u_mesh_4 = (double **)malloc((N_max + 1) * sizeof(double *));
	for (int i = 0; i < N_max + 1; i++) {
		u_mesh_4[i] = (double *)malloc((I_max + 1) * sizeof(double));
		for (int j = 0; j < I_max + 1; j++) {
			u_mesh_4[i][j] = 0;
		}
	}

	// f_mesh_1 for rho*v*S, f_mesh_2 for rho*v^2*S, f_mesh_3 for rho*v*energy*S, f_mesh_4 for H_phi*v*S
	// f(u) = {rho*v*S, rho*v^2*S, rho*v*energy*S, H_phi*v*S}

	double **f_mesh_1 = (double **)malloc((N_max + 1) * sizeof(double *));
	for (int i = 0; i < N_max + 1; i++) {
		f_mesh_1[i] = (double *)malloc((I_max + 1) * sizeof(double));
		for (int j = 0; j < I_max + 1; j++) {
			f_mesh_1[i][j] = 0;
		}
	}

	double **f_mesh_2 = (double **)malloc((N_max + 1) * sizeof(double *));
	for (int i = 0; i < N_max + 1; i++) {
		f_mesh_2[i] = (double *)malloc((I_max + 1) * sizeof(double));
		for (int j = 0; j < I_max + 1; j++) {
			f_mesh_2[i][j] = 0;
		}
	}
	
	double **f_mesh_3 = (double **)malloc((N_max + 1) * sizeof(double *));
	for (int i = 0; i < N_max + 1; i++) {
		f_mesh_3[i] = (double *)malloc((I_max + 1) * sizeof(double));
		for (int j = 0; j < I_max + 1; j++) {
			f_mesh_3[i][j] = 0;
		}
	}
	
	double **f_mesh_4 = (double **)malloc((N_max + 1) * sizeof(double *));
	for (int i = 0; i < N_max + 1; i++) {
		f_mesh_4[i] = (double *)malloc((I_max + 1) * sizeof(double));
		for (int j = 0; j < I_max + 1; j++) {
			f_mesh_4[i][j] = 0;
		}
	}

	// g(u) = {0, -S*dp/dz, -p*d/dz(v*S), 0}

	// initial condition

	double v_0 = 0.1;
	double gamma = 1.67;
	double beta = 1.0;
	double mu_0 = 0.7;

	double *S_mesh = (double *)malloc((I_max + 1) * sizeof(double));
	for (int i = 0; i < I_max + 1; i++) {
		S_mesh[i] = 2.0 * (i * h - 0.5) * (i * h - 0.5) + 0.5;
	}

	double *v_mesh = (double *)malloc((I_max + 1) * sizeof(double));
	for (int i = 0; i < I_max + 1; i++) {
		v_mesh[i] = v_0;
	}

	double *p_mesh = (double *)malloc((I_max + 1) * sizeof(double));
	for (int i = 0; i < I_max + 1; i++) {
		p_mesh[i] = beta / 2.0 * (1.0 - 0.8 * i * h);
	}

	double *rho_mesh = (double *)malloc((I_max + 1) * sizeof(double));
	for (int i = 0; i < I_max + 1; i++) {
		rho_mesh[i] = 1.0 - 0.8 * i * h;
	}

	double *H_phi_mesh = (double *)malloc((I_max + 1) * sizeof(double));
	for (int i = 0; i < I_max + 1; i++) {
		H_phi_mesh[i] = 1.0 - 0.8 * i * h;
	}

	for (int i = 0; i < I_max + 1; i++) {
		u_mesh_1[0][i] = rho_mesh[i] * S_mesh[i];
		u_mesh_2[0][i] = rho_mesh[i] * v_0 * S_mesh[i];
		u_mesh_3[0][i] = p_mesh[i] * S_mesh[i];
		u_mesh_4[0][i] = H_phi_mesh[i] * S_mesh[i];
		f_mesh_1[0][i] = u_mesh_1[0][i] * v_mesh[i];
		f_mesh_2[0][i] = u_mesh_2[0][i] * v_mesh[i];
		f_mesh_3[0][i] = u_mesh_3[0][i] * v_mesh[i];
		f_mesh_4[0][i] = u_mesh_4[0][i] * v_mesh[i];
	}

	// mesh filling

	double mu_left;
    double mu_right;
    double C_left;
    double C_0;
    double C_right;

	for (int n = 0; n < N_max; n++) {
		// filling mesh for u

		for (int i = 1; i < I_max; i++) {
			mu_left = mu_0 + (abs(v_mesh[i - 1]) + abs(v_mesh[i])) / 4.0;
			mu_right = mu_0 + (abs(v_mesh[i]) + abs(v_mesh[i + 1])) / 4.0;
			C_left = tau / h * (v_mesh[i - 1] + v_mesh[i]) / 4.0 + tau / h * mu_left;
			C_0 = 1.0 - tau / h * (v_mesh[i + 1] - v_mesh[i - 1]) / 4.0 - tau / h * (mu_left + mu_right);
			C_right = - tau / h * (v_mesh[i] + v_mesh[i + 1]) / 4.0 + tau / h * mu_right;

			u_mesh_1[n + 1][i] = C_left * u_mesh_1[n][i - 1] + C_0 * u_mesh_1[n][i] + C_right * u_mesh_1[n][i + 1];
		}
		u_mesh_1[n + 1][0] = rho_mesh[0] * S_mesh[0];
		u_mesh_1[n + 1][I_max] = u_mesh_1[n + 1][I_max - 1];

		for (int i = 1; i < I_max; i++) {
			mu_left = mu_0 + (abs(v_mesh[i - 1]) + abs(v_mesh[i])) / 4.0;
			mu_right = mu_0 + (abs(v_mesh[i]) + abs(v_mesh[i + 1])) / 4.0;
			C_left = tau / h * (v_mesh[i - 1] + v_mesh[i]) / 4.0 + tau / h * mu_left;
			C_0 = 1.0 - tau / h * (v_mesh[i + 1] - v_mesh[i - 1]) / 4.0 - tau / h * (mu_left + mu_right);
			C_right = - tau / h * (v_mesh[i] + v_mesh[i + 1]) / 4.0 + tau / h * mu_right;

			u_mesh_2[n + 1][i] = C_left * u_mesh_2[n][i - 1] + C_0 * u_mesh_2[n][i] + C_right * u_mesh_2[n][i + 1] -
				tau / (2.0 * h) * 1.0 / gamma * S_mesh[i] * (p_mesh[i + 1] - p_mesh[i - 1]);
		}
		u_mesh_2[n + 1][0] = rho_mesh[0] * v_mesh[0] * S_mesh[0];
		u_mesh_2[n + 1][I_max] = u_mesh_2[n + 1][I_max - 1];

		for (int i = 1; i < I_max; i++) {
			mu_left = mu_0 + (abs(v_mesh[i - 1]) + abs(v_mesh[i])) / 4.0;
			mu_right = mu_0 + (abs(v_mesh[i]) + abs(v_mesh[i + 1])) / 4.0;
			C_left = tau / h * (v_mesh[i - 1] + v_mesh[i]) / 4.0 + tau / h * mu_left;
			C_0 = 1.0 - tau / h * (v_mesh[i + 1] - v_mesh[i - 1]) / 4.0 - tau / h * (mu_left + mu_right);
			C_right = - tau / h * (v_mesh[i] + v_mesh[i + 1]) / 4.0 + tau / h * mu_right;

			u_mesh_3[n + 1][i] = C_left * u_mesh_3[n][i - 1] + C_0 * u_mesh_3[n][i] + C_right * u_mesh_3[n][i + 1] + 
				tau / (2.0 * h) * (gamma - 1.0) * v_mesh[i] * S_mesh[i] * (p_mesh[i + 1] - p_mesh[i - 1]);
		}
		u_mesh_3[n + 1][0] = p_mesh[0] * S_mesh[0];
		u_mesh_3[n + 1][I_max] = u_mesh_3[n + 1][I_max - 1];

		for (int i = 1; i < I_max; i++) {
			mu_left = mu_0 + (abs(v_mesh[i - 1]) + abs(v_mesh[i])) / 4.0;
			mu_right = mu_0 + (abs(v_mesh[i]) + abs(v_mesh[i + 1])) / 4.0;
			C_left = tau / h * (v_mesh[i - 1] + v_mesh[i]) / 4.0 + tau / h * mu_left;
			C_0 = 1.0 - tau / h * (v_mesh[i + 1] - v_mesh[i - 1]) / 4.0 - tau / h * (mu_left + mu_right);
			C_right = - tau / h * (v_mesh[i] + v_mesh[i + 1]) / 4.0 + tau / h * mu_right;

			u_mesh_4[n + 1][i] = C_left * u_mesh_4[n][i - 1] + C_0 * u_mesh_4[n][i] + C_right * u_mesh_4[n][i + 1];
		}
		u_mesh_4[n + 1][0] = rho_mesh[0] * S_mesh[0];
		u_mesh_4[n + 1][I_max] = u_mesh_4[n + 1][I_max - 1];

		// calculate the velocity as u_2 / u_1

		v_mesh[0] = u_mesh_2[n + 1][1];	// proven analytically
		for (int i = 1; i < I_max + 1; i++) {
			v_mesh[i] = u_mesh_2[n + 1][i] / u_mesh_1[n + 1][i];
		}

		// filling mesh for f(u)

		for (int i = 0; i < I_max + 1; i++) {
			f_mesh_1[n + 1][i] = u_mesh_1[n + 1][i] * v_mesh[i];
		}
		for (int i = 0; i < I_max + 1; i++) {
			f_mesh_2[n + 1][i] = u_mesh_2[n + 1][i] * v_mesh[i];
		}
		for (int i = 0; i < I_max + 1; i++) {
			f_mesh_3[n + 1][i] = u_mesh_3[n + 1][i] * v_mesh[i];
		}
		for (int i = 0; i < I_max + 1; i++) {
			f_mesh_4[n + 1][i] = u_mesh_4[n + 1][i] * v_mesh[i];
		}
	}

	// checkout

	double *x_lst = (double *)malloc((I_max + 1) * sizeof(double));
	for (int i = 0; i < I_max + 1; i++) {
		x_lst[i] = 0;
	}
	double *rho_lst = (double *)malloc((I_max + 1) * sizeof(double));
	for (int i = 0; i < I_max + 1; i++) {
		rho_lst[i] = 0;
	}
	double *v_lst = (double *)malloc((I_max + 1) * sizeof(double));
	for (int i = 0; i < I_max + 1; i++) {
		v_lst[i] = 0;
	}
	double *e_lst = (double *)malloc((I_max + 1) * sizeof(double));
	for (int i = 0; i < I_max + 1; i++) {
		e_lst[i] = 0;
	}
	double *H_lst = (double *)malloc((I_max + 1) * sizeof(double));
	for (int i = 0; i < I_max + 1; i++) {
		H_lst[i] = 0;
	}

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
		v_lst[i] = f_mesh_1[n][i] / u_mesh_1[n][i];
	}

	// for energy
	for (int i = 0; i < I_max + 1; i++) {
		e_lst[i] = u_mesh_3[n][i] / ((gamma - 1) * u_mesh_1[n][i]);
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
	for(int i = 0; i < N_max + 1; i++) {
		delete [] u_mesh_1[i];
	}
	delete [] u_mesh_1;
	for(int i = 0; i < N_max + 1; i++) {
		delete [] u_mesh_2[i];
	}
	delete [] u_mesh_2;
	for(int i = 0; i < N_max + 1; i++) {
		delete [] u_mesh_3[i];
	}
	delete [] u_mesh_3;
	for(int i = 0; i < N_max + 1; i++) {
		delete [] u_mesh_4[i];
	}
	delete [] u_mesh_4;
	for(int i = 0; i < N_max + 1; i++) {
		delete [] f_mesh_1[i];
	}
	delete [] f_mesh_1;
	for(int i = 0; i < N_max + 1; i++) {
		delete [] f_mesh_2[i];
	}
	delete [] f_mesh_2;
	for(int i = 0; i < N_max + 1; i++) {
		delete [] f_mesh_3[i];
	}
	delete [] f_mesh_3;
	for(int i = 0; i < N_max + 1; i++) {
		delete [] f_mesh_4[i];
	}
	delete [] f_mesh_4;
	delete [] S_mesh;
	delete [] v_mesh;
	delete [] p_mesh;
	delete [] rho_mesh;
	delete [] H_phi_mesh;
	delete [] x_lst;
	delete [] rho_lst;
	delete [] v_lst;
	delete [] e_lst;
	delete [] H_lst;

	return 0;
}

int round(double d) {
	return floor(d + 0.5);
}