#include<stdio.h>
#include<iostream>
#include<fstream>
#include<cmath>
#include<omp.h>

void memory_allocation_2D(double**& array, int rows, int columns) {

	array = new double* [rows];
	for (int i = 0; i < rows; i++) {
		array[i] = new double[columns];
		for (int j = 0; j < columns; j++) {
			array[i][j] = 0;
		}
	}
}

void memory_clearing_2D(double**& array, int rows) {
	for (int i = 0; i < rows; i++) {
		delete[] array[i];
	}
	delete[] array;
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



int main(int argc, char* argv[]) {

	// task data

	double gamma = 1.67;
	double beta = 0.05;

	double H_z0 = 0;

	// discrete solution area

	double T = 10.0;
	double t = 0.0;
	double dt = 0.001;

	int L_max = 80;
	int M_max = 40;

	double dz = 1.0 / L_max;
	double dy = 1.0 / M_max;

	// parallel parameters

	int procs = atoi(argv[1]);
	omp_set_num_threads(procs);

	double begin;
	double end;
	double total;

	// creating arrays for u components

	// u = {rho*r, rho*r*v_z, rho*r*v_r, rho*r*v_phi, rho*r*energy,  H_phi, H_z*r, H_r*r}

	double** u_1;
	double** u_2;
	double** u_3;
	double** u_4;
	double** u_5;
	double** u_6;
	double** u_7;
	double** u_8;
	memory_allocation_2D(u_1, L_max + 1, M_max + 1);
	memory_allocation_2D(u_2, L_max + 1, M_max + 1);
	memory_allocation_2D(u_3, L_max + 1, M_max + 1);
	memory_allocation_2D(u_4, L_max + 1, M_max + 1);
	memory_allocation_2D(u_5, L_max + 1, M_max + 1);
	memory_allocation_2D(u_6, L_max + 1, M_max + 1);
	memory_allocation_2D(u_7, L_max + 1, M_max + 1);
	memory_allocation_2D(u_8, L_max + 1, M_max + 1);

	double** u0_1;
	double** u0_2;
	double** u0_3;
	double** u0_4;
	double** u0_5;
	double** u0_6;
	double** u0_7;
	double** u0_8;
	memory_allocation_2D(u0_1, L_max + 1, M_max + 1);
	memory_allocation_2D(u0_2, L_max + 1, M_max + 1);
	memory_allocation_2D(u0_3, L_max + 1, M_max + 1);
	memory_allocation_2D(u0_4, L_max + 1, M_max + 1);
	memory_allocation_2D(u0_5, L_max + 1, M_max + 1);
	memory_allocation_2D(u0_6, L_max + 1, M_max + 1);
	memory_allocation_2D(u0_7, L_max + 1, M_max + 1);
	memory_allocation_2D(u0_8, L_max + 1, M_max + 1);

	// creating arrays for physical parameters

	double** rho;
	double** v_r;
	double** v_phi;
	double** v_z;
	double** e;
	double** p;
	double** P;
	double** H_r;
	double** H_phi;
	double** H_z;
	memory_allocation_2D(rho, L_max + 1, M_max + 1);
	memory_allocation_2D(v_r, L_max + 1, M_max + 1);
	memory_allocation_2D(v_phi, L_max + 1, M_max + 1);
	memory_allocation_2D(v_z, L_max + 1, M_max + 1);
	memory_allocation_2D(e, L_max + 1, M_max + 1);
	memory_allocation_2D(p, L_max + 1, M_max + 1);
	memory_allocation_2D(P, L_max + 1, M_max + 1);
	memory_allocation_2D(H_r, L_max + 1, M_max + 1);
	memory_allocation_2D(H_phi, L_max + 1, M_max + 1);
	memory_allocation_2D(H_z, L_max + 1, M_max + 1);

	double** r;
	double** r_z;
	memory_allocation_2D(r, L_max + 1, M_max + 1);
	memory_allocation_2D(r_z, L_max + 1, M_max + 1);
	double* R = new double[L_max + 1];
	double* dr = new double[L_max + 1];

	double** v_y;
	double** H_y;
	memory_allocation_2D(v_y, L_max + 1, M_max + 1);
	memory_allocation_2D(H_y, L_max + 1, M_max + 1);

	// filling zeros

	for (int l = 0; l < L_max + 1; l++) {
		R[l] = 0;
		dr[l] = 0;
	}

	// creating arrays for axes

	double r_0 = (r1(0) + r2(0)) / 2.0;

	for (int l = 0; l < L_max + 1; l++) {
		R[l] = r2(l * dz) - r1(l * dz);
		dr[l] = R[l] / M_max;

		for (int m = 0; m < M_max + 1; m++) {
			r[l][m] = (1 - m * dy) * r1(l * dz) + m * dy * r2(l * dz);
			r_z[l][m] = (1 - m * dy) * der_r1(l * dz) + m * dy * der_r2(l * dz);
		}
	}

	// initial condition

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
			u0_1[l][m] = rho[l][m] * r[l][m];
			u0_2[l][m] = rho[l][m] * v_z[l][m] * r[l][m];
			u0_3[l][m] = rho[l][m] * v_r[l][m] * r[l][m];
			u0_4[l][m] = rho[l][m] * v_phi[l][m] * r[l][m];
			u0_5[l][m] = rho[l][m] * e[l][m] * r[l][m];
			u0_6[l][m] = H_phi[l][m];
			u0_7[l][m] = H_z[l][m] * r[l][m];
			u0_8[l][m] = H_r[l][m] * r[l][m];
		}
	}

	// start count time

	begin = omp_get_wtime();

	// start time

	while (t < T) {

		// counting central part

#pragma omp parallel for collapse(2)

		for (int l = 1; l < L_max; l++) {
			for (int m = 1; m < M_max; m++) {
				u_1[l][m] = u0_1[l][m] + 0.25 * (u0_1[l + 1][m] - 2 * u0_1[l][m] + u0_1[l - 1][m]) +
					0.25 * (u0_1[l][m + 1] - 2 * u0_1[l][m] + u0_1[l][m - 1]) +
					dt * (0 -
						(u0_1[l + 1][m] * v_z[l + 1][m] - u0_1[l - 1][m] * v_z[l - 1][m]) / (2 * dz) -
						(u0_1[l][m + 1] * v_r[l][m + 1] - u0_1[l][m - 1] * v_r[l][m - 1]) / (2 * dr[l]));

				u_2[l][m] = u0_2[l][m] + 0.25 * (u0_2[l + 1][m] - 2 * u0_2[l][m] + u0_2[l - 1][m]) +
					0.25 * (u0_2[l][m + 1] - 2 * u0_2[l][m] + u0_2[l][m - 1]) +
					dt * (-((P[l + 1][m] - pow(H_z[l + 1][m], 2)) * r[l + 1][m] -
						(P[l - 1][m] - pow(H_z[l - 1][m], 2)) * r[l - 1][m]) / (2 * dz) +
						(H_z[l][m + 1] * H_r[l][m + 1] * r[l][m + 1] -
							H_z[l][m - 1] * H_r[l][m - 1] * r[l][m - 1]) / (2 * dr[l]) -
						(u0_2[l + 1][m] * v_z[l + 1][m] - u0_2[l - 1][m] * v_z[l - 1][m]) / (2 * dz) -
						(u0_2[l][m + 1] * v_r[l][m + 1] - u0_2[l][m - 1] * v_r[l][m - 1]) / (2 * dr[l]));

				u_3[l][m] = u0_3[l][m] + 0.25 * (u0_3[l + 1][m] - 2 * u0_3[l][m] + u0_3[l - 1][m]) +
					0.25 * (u0_3[l][m + 1] - 2 * u0_3[l][m] + u0_3[l][m - 1]) +
					dt * (rho[l][m] * pow(v_phi[l][m], 2) + P[l][m] - pow(H_phi[l][m], 2) +
						(H_z[l + 1][m] * H_r[l + 1][m] * r[l + 1][m] -
							H_z[l - 1][m] * H_r[l - 1][m] * r[l - 1][m]) / (2 * dz) -
						((P[l][m + 1] - pow(H_r[l][m + 1], 2)) * r[l][m + 1] -
							(P[l][m - 1] - pow(H_r[l][m - 1], 2)) * r[l][m - 1]) / (2 * dr[l]) -
						(u0_3[l + 1][m] * v_z[l + 1][m] - u0_3[l - 1][m] * v_z[l - 1][m]) / (2 * dz) -
						(u0_3[l][m + 1] * v_r[l][m + 1] - u0_3[l][m - 1] * v_r[l][m - 1]) / (2 * dr[l]));


				u_4[l][m] = u0_4[l][m] + 0.25 * (u0_4[l + 1][m] - 2 * u0_4[l][m] + u0_4[l - 1][m]) +
					0.25 * (u0_4[l][m + 1] - 2 * u0_4[l][m] + u0_4[l][m - 1]) +
					dt * (-rho[l][m] * r[l][m] * v_phi[l][m] + H_phi[l][m] * H_r[l][m] +
						(H_phi[l + 1][m] * H_z[l + 1][m] * r[l + 1][m] -
							H_phi[l - 1][m] * H_z[l - 1][m] * r[l - 1][m]) / (2 * dz) +
						(H_phi[l][m + 1] * H_r[l][m + 1] * r[l][m + 1] -
							H_phi[l][m - 1] * H_r[l][m - 1] * r[l][m - 1]) / (2 * dr[l]) -
						(u0_4[l + 1][m] * v_z[l + 1][m] - u0_4[l - 1][m] * v_z[l - 1][m]) / (2 * dz) -
						(u0_4[l][m + 1] * v_r[l][m + 1] - u0_4[l][m - 1] * v_r[l][m - 1]) / (2 * dr[l]));

				u_5[l][m] = u0_5[l][m] + 0.25 * (u0_5[l + 1][m] - 2 * u0_5[l][m] + u0_5[l - 1][m]) +
					0.25 * (u0_5[l][m + 1] - 2 * u0_5[l][m] + u0_5[l][m - 1]) +
					dt * (-p[l][m] * ((v_z[l + 1][m] * r[l + 1][m] -
						v_z[l - 1][m] * r[l - 1][m]) / (2 * dz) +
						(v_r[l][m + 1] * r[l][m + 1] -
							v_r[l][m - 1] * r[l][m - 1]) / (2 * dr[l])) -
						(u0_5[l + 1][m] * v_z[l + 1][m] - u0_5[l - 1][m] * v_z[l - 1][m]) / (2 * dz) -
						(u0_5[l][m + 1] * v_r[l][m + 1] - u0_5[l][m - 1] * v_r[l][m - 1]) / (2 * dr[l]));

				u_6[l][m] = u0_6[l][m] + 0.25 * (u0_6[l + 1][m] - 2 * u0_6[l][m] + u0_6[l - 1][m]) +
					0.25 * (u0_6[l][m + 1] - 2 * u0_6[l][m] + u0_6[l][m - 1]) +
					dt * ((H_z[l + 1][m] * v_phi[l + 1][m] -
						H_z[l - 1][m] * v_phi[l - 1][m]) / (2 * dz) +
						(H_r[l][m + 1] * v_phi[l][m + 1] -
							H_r[l][m - 1] * v_phi[l][m - 1]) / (2 * dr[l]) -
						(u0_6[l + 1][m] * v_z[l + 1][m] - u0_6[l - 1][m] * v_z[l - 1][m]) / (2 * dz) -
						(u0_6[l][m + 1] * v_r[l][m + 1] - u0_6[l][m - 1] * v_r[l][m - 1]) / (2 * dr[l]));

				u_7[l][m] = u0_7[l][m] + 0.25 * (u0_7[l + 1][m] - 2 * u0_7[l][m] + u0_7[l - 1][m]) +
					0.25 * (u0_7[l][m + 1] - 2 * u0_7[l][m] + u0_7[l][m - 1]) +
					dt * ((H_r[l][m + 1] * v_z[l][m + 1] * r[l][m + 1] -
						H_r[l][m - 1] * v_z[l][m - 1] * r[l][m - 1]) / (2 * dr[l]) -
						(u0_7[l][m + 1] * v_r[l][m + 1] - u0_7[l][m - 1] * v_r[l][m - 1]) / (2 * dr[l]));


				u_8[l][m] = u0_8[l][m] + 0.25 * (u0_8[l + 1][m] - 2 * u0_8[l][m] + u0_8[l - 1][m]) +
					0.25 * (u0_8[l][m + 1] - 2 * u0_8[l][m] + u0_8[l][m - 1]) +
					dt * ((H_z[l + 1][m] * v_r[l + 1][m] * r[l + 1][m] -
						H_z[l - 1][m] * v_r[l - 1][m] * r[l - 1][m]) / (2 * dz) -
						(u0_8[l + 1][m] * v_r[l + 1][m] - u0_8[l - 1][m] * v_z[l - 1][m]) / (2 * dz));
			}
		}

		// left boundary condition

#pragma omp parallel for

		for (int m = 0; m < M_max + 1; m++) {
			rho[0][m] = 1.0;
			v_phi[0][m] = 0;
			// v_r[0][m] = v_z[0][m] * r_z[0][m];
			H_phi[0][m] = r_0 / r[0][m];
			H_z[0][m] = H_z0;
			// H_r[0][m] = H_z[0][m] * r_z[0][m];
			e[0][m] = beta / (2.0 * (gamma - 1.0)) * pow(rho[0][m], gamma - 1.0);
		}

#pragma omp parallel for

		for (int m = 0; m < M_max + 1; m++) {
			u_1[0][m] = rho[0][m] * r[0][m];
			u_2[0][m] = u_2[1][m];
			u_3[0][m] = u_2[0][m] * r_z[0][m];
			u_4[0][m] = rho[0][m] * v_phi[0][m] * r[0][m];
			u_5[0][m] = rho[0][m] * e[0][m] * r[0][m];
			u_6[0][m] = H_phi[0][m];
			u_7[0][m] = H_z[0][m] * r[0][m];
			u_8[0][m] = u_7[0][m] * r_z[0][m];
		}

		// up and down boundary condition

#pragma omp parallel for

		for (int l = 1; l < L_max + 1; l++) {
			u_1[l][0] = u_1[l][1];
			u_2[l][0] = u_2[l][1];
			u_3[l][0] = u_2[l][0] * r_z[l][0];
			u_4[l][0] = u_4[l][1];
			u_5[l][0] = u_5[l][1];
			u_6[l][0] = u_6[l][1];
			u_7[l][0] = u_7[l][1];
			u_8[l][0] = u_7[l][0] * r_z[l][0];

			u_1[l][M_max] = u_1[l][M_max - 1];
			u_2[l][M_max] = u_2[l][M_max - 1];
			u_3[l][M_max] = u_2[l][M_max] * r_z[l][M_max];
			u_4[l][M_max] = u_4[l][M_max - 1];
			u_5[l][M_max] = u_5[l][M_max - 1];
			u_6[l][M_max] = u_6[l][M_max - 1];
			u_7[l][M_max] = u_7[l][M_max - 1];
			u_8[l][M_max] = u_7[l][M_max] * r_z[l][M_max];
		}

		// right boundary condition

#pragma omp parallel for

		for (int m = 1; m < M_max; m++) {
			u_1[L_max][m] = u_1[L_max - 1][m];
			u_2[L_max][m] = u_2[L_max - 1][m];
			u_3[L_max][m] = u_3[L_max - 1][m];
			u_4[L_max][m] = u_4[L_max - 1][m];
			u_5[L_max][m] = u_5[L_max - 1][m];
			u_6[L_max][m] = u_6[L_max - 1][m];
			u_7[L_max][m] = u_7[L_max - 1][m];
			u_8[L_max][m] = u_8[L_max - 1][m];
		}

#pragma omp parallel for collapse(2)

		// data update

		for (int l = 0; l < L_max + 1; l++) {
			for (int m = 0; m < M_max + 1; m++) {
				rho[l][m] = u_1[l][m] / r[l][m];
				v_z[l][m] = u_2[l][m] / u_1[l][m];
				v_r[l][m] = u_3[l][m] / u_1[l][m];
				v_phi[l][m] = u_4[l][m] / u_1[l][m];
				v_y[l][m] = v_r[l][m] - v_z[l][m] * r_z[l][m];

				H_phi[l][m] = u_6[l][m];
				H_z[l][m] = u_7[l][m] / r[l][m];
				H_r[l][m] = u_8[l][m] / r[l][m];
				H_y[l][m] = H_r[l][m] - H_z[l][m] * r_z[l][m];

				e[l][m] = u_5[l][m] / u_1[l][m];
				p[l][m] = (gamma - 1) * rho[l][m] * e[l][m];
				P[l][m] = p[l][m] + 1.0 / 2.0 * (pow(H_z[l][m], 2) + pow(H_r[l][m], 2) + pow(H_phi[l][m], 2));
			}
		}

#pragma omp parallel for collapse(2)

		for (int l = 0; l < L_max + 1; l++) {
			for (int m = 0; m < M_max + 1; m++) {
				u0_1[l][m] = u_1[l][m];
				u0_2[l][m] = u_2[l][m];
				u0_3[l][m] = u_3[l][m];
				u0_4[l][m] = u_4[l][m];
				u0_5[l][m] = u_5[l][m];
				u0_6[l][m] = u_6[l][m];
				u0_7[l][m] = u_7[l][m];
				u0_8[l][m] = u_8[l][m];
			}
		}

		// time step

		t += dt;

		// checkout

		//printf("%lf		%lf		%lf		%lf		%lf		%lf\n", t, rho[20][40], v_z[20][40], v_phi[20][40], e[20][40], H_phi[20][40]);

	}

	// finish count time

	end = omp_get_wtime();

	// calculation the time

	total = end - begin;
	printf("Calculation time : %lf sec\n", total);

	// output results in file

	std::ofstream out("21-2D_MHD_viscosity.plt");
	int np = (L_max + 1) * (M_max + 1);
	int ne = L_max * M_max;
	double hfr;
	out << "VARIABLES=\n";
	out << "\"r\"\n\"z\"\n\"Rho\"\n\"Vz\"\n\"Vr\"\n\"Vl\"\n\"Vphi\"\n\"Energy\"\n\"Hz\"\n\"Hr\"\n\"Hphi*r\"\n\"Hphi\"\n";
	out << "ZONE \n F=FEPOINT, ET=Quadrilateral, N=" << np << " E=" << ne << "\n ";
	for (int m = 0; m < M_max + 1; m++) {
		for (int l = 0; l < L_max + 1; l++) {
			hfr = H_phi[l][m] * r[l][m];
			out << l * dz << " " << r[l][m] << " " << rho[l][m] << " " <<
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
			out << i1 << " " << i2 << " " << i3 << " " << i4 << "\n";
		}
	}

	out.close();

	// clear memory

	memory_clearing_2D(u_1, L_max + 1);
	memory_clearing_2D(u_2, L_max + 1);
	memory_clearing_2D(u_3, L_max + 1);
	memory_clearing_2D(u_4, L_max + 1);
	memory_clearing_2D(u_5, L_max + 1);
	memory_clearing_2D(u_6, L_max + 1);
	memory_clearing_2D(u_7, L_max + 1);
	memory_clearing_2D(u_8, L_max + 1);

	memory_clearing_2D(u0_1, L_max + 1);
	memory_clearing_2D(u0_2, L_max + 1);
	memory_clearing_2D(u0_3, L_max + 1);
	memory_clearing_2D(u0_4, L_max + 1);
	memory_clearing_2D(u0_5, L_max + 1);
	memory_clearing_2D(u0_6, L_max + 1);
	memory_clearing_2D(u0_7, L_max + 1);
	memory_clearing_2D(u0_8, L_max + 1);

	memory_clearing_2D(rho, L_max + 1);
	memory_clearing_2D(v_r, L_max + 1);
	memory_clearing_2D(v_phi, L_max + 1);
	memory_clearing_2D(v_z, L_max + 1);
	memory_clearing_2D(e, L_max + 1);
	memory_clearing_2D(p, L_max + 1);
	memory_clearing_2D(P, L_max + 1);
	memory_clearing_2D(H_r, L_max + 1);
	memory_clearing_2D(H_phi, L_max + 1);
	memory_clearing_2D(H_z, L_max + 1);

	memory_clearing_2D(r, L_max + 1);
	memory_clearing_2D(r_z, L_max + 1);
	delete[] R;
	delete[] dr;

	memory_clearing_2D(v_y, L_max + 1);
	memory_clearing_2D(H_y, L_max + 1);

	return 0;
}
