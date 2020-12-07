#include "LBM_D2Q9_KBC.h"

void LBM_D2Q9_KBC::drawing_liquid(void) {
	GLint i = 0, j = 0, ileft = 0, iright = 0, jup = 0, jdown = 0, k = 0;

	glPushMatrix();
	glColor3f(1.0, 0.0, 0.0);
	glPointSize(5.0f);
	glBegin(GL_POINTS);
	double Vort_max = 0.0;
	for (i = 0; i <= l; i++) {
		ileft = (i > 0) ? (i - 1) : (l);
		iright = (i < l) ? (i + 1) : (0);
		for (j = 0; j <= m; j++) {
			jdown = (j > 0) ? (j - 1) : (m);
			jup = (j < m) ? (j + 1) : (0);
			Vorticity[i][j] = 0.5 * (Nodes[i][jup].u[0] - Nodes[i][jdown].u[0]) - \
				0.5 * (Nodes[iright][j].u[1] - Nodes[ileft][j].u[1]);
			if (abs(Vorticity[i][j]) > Vort_max) { Vort_max = Vorticity[i][j]; }
		}
	}
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			if (Nodes[i][j].status == "L") {
				if (Vorticity[i][j] > 0.0) {
					glColor3f(abs(Vorticity[i][j]) / Vort_max, 0.0, 0.0);
				}
				else {
					glColor3f(0.0, 0.0, abs(Vorticity[i][j]) / Vort_max);
				}
				glVertex2f(i * Lx / l, j * Ly / m);
			}
		}
	}
	glEnd();
	glPopMatrix();
}

void LBM_D2Q9_KBC::dump_file() {
	int i = 0, j = 0;
	FILE* fp;
	if (!(fp = fopen("data.plt", "w")))
	{
		printf("Ð´ÈëÎÄ¼þÊ§°Ü!\n");
		exit(1);
	}
	fprintf(fp, "TITLE =\"lattice boltzmann\"\n");
	fprintf(fp, "VARIABLES = \"X\" \"Y\" \"RHO\" \"U\" \"V\" \"VORT\" \n");
	fprintf(fp, "ZONE I=%d, J=%d\n", l + 1, m + 1);
	fprintf(fp, "F=POINT\n");
	for (j = 0; j <= m; j++) {
		for (i = 0; i <= l; i++) {
			fprintf(fp, "%d %d %lf %lf %lf %lf\n", i, j, Nodes[i][j].rho, Nodes[i][j].u[0], Nodes[i][j].u[1], Vorticity[i][j]);
		}
	}
	fclose(fp);
}

void LBM_D2Q9_KBC::Calculate_KBC_params(double rho, double* f, double* feq, double& KBC_gamma, double KBC_Delta_s[], double KBC_Delta_h[]) {
	double Moment[3][3], Moment_eq[3][3];
	vector<vector<vector<int> > > C_matrix(3, vector<vector<int> >(3, vector<int>(Q, 0)));
	C_matrix[0][0] = vector<int>{ 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	C_matrix[1][0] = vector<int>{ 0, 1, 0, -1, 0, 1, -1, -1, 1 };
	C_matrix[0][1] = vector<int>{ 0, 0, 1, 0, -1, 1, 1, -1, -1 };
	C_matrix[2][0] = vector<int>{ 0, 1, 0, 1, 0, 1, 1, 1, 1 };
	C_matrix[0][2] = vector<int>{ 0, 0, 1, 0, 1, 1, 1, 1, 1 };
	C_matrix[1][1] = vector<int>{ 0, 0, 0, 0, 0, 1, -1, 1, -1 };
	C_matrix[2][2] = vector<int>{ 0, 0, 0, 0, 0, 1, 1, 1, 1 };
	C_matrix[2][1] = vector<int>{ 0, 0, 0, 0, 0, 1, 1, -1, -1 };
	C_matrix[1][2] = vector<int>{ 0, 0, 0, 0, 0, 1, -1, -1, 1 };
	for (int M_i = 0; M_i < 3; M_i++) {
		for (int M_j = 0; M_j < 3; M_j++) {
			Moment[M_i][M_j] = 0.0;
			Moment_eq[M_i][M_j] = 0.0;
			for (int k = 0; k < Q; k++) {
				Moment[M_i][M_j] += f[k] * C_matrix[M_i][M_j][k];
				Moment_eq[M_i][M_j] += feq[k] * C_matrix[M_i][M_j][k];
			}
		}
	}

	//KBC - A
	double KBC_T = Moment[2][0] + Moment[0][2];
	double KBC_N = Moment[2][0] - Moment[0][2];
	double KBC_PIxy = Moment[1][1];
	double KBC_Qxxy = Moment[2][1], KBC_Qxyy = Moment[1][2];
	double KBC_A = Moment[2][2];

	KBC_Delta_s[0] = (-KBC_T);

	KBC_Delta_s[1] = 0.5 * (0.5 * (KBC_T + KBC_N));
	KBC_Delta_s[3] = 0.5 * (0.5 * (KBC_T + KBC_N));

	KBC_Delta_s[2] = 0.5 * (0.5 * (KBC_T - KBC_N));
	KBC_Delta_s[4] = 0.5 * (0.5 * (KBC_T - KBC_N));

	KBC_Delta_s[5] = 0.25 * (KBC_PIxy);
	KBC_Delta_s[6] = 0.25 * (-KBC_PIxy);
	KBC_Delta_s[7] = 0.25 * (KBC_PIxy);
	KBC_Delta_s[8] = 0.25 * (-KBC_PIxy); 

	KBC_T = Moment_eq[2][0] + Moment_eq[0][2];
	KBC_N = Moment_eq[2][0] - Moment_eq[0][2];
	KBC_PIxy = Moment_eq[1][1];
	KBC_Qxxy = Moment_eq[2][1], KBC_Qxyy = Moment_eq[1][2];
	KBC_A = Moment_eq[2][2];

	KBC_Delta_s[0] -= (-KBC_T);

	KBC_Delta_s[1] -= 0.5 * (0.5 * (KBC_T + KBC_N));
	KBC_Delta_s[3] -= 0.5 * (0.5 * (KBC_T + KBC_N));

	KBC_Delta_s[2] -= 0.5 * (0.5 * (KBC_T - KBC_N));
	KBC_Delta_s[4] -= 0.5 * (0.5 * (KBC_T - KBC_N));

	KBC_Delta_s[5] -= 0.25 * (KBC_PIxy);
	KBC_Delta_s[6] -= 0.25 * (-KBC_PIxy);
	KBC_Delta_s[7] -= 0.25 * (KBC_PIxy);
	KBC_Delta_s[8] -= 0.25 * (-KBC_PIxy);

	double gamma_nominator = 0.0, gamma_determinator = 0.0;
	KBC_Delta_h[0] = f[0] - feq[0] - KBC_Delta_s[0];
	for (int k = 1; k < Q; k++) {
		KBC_Delta_h[k] = f[k] - feq[k] - KBC_Delta_s[k];
		gamma_nominator += KBC_Delta_s[k] * KBC_Delta_h[k] / (feq[k] + 1e-7);
		gamma_determinator += KBC_Delta_h[k] * KBC_Delta_h[k] / (feq[k] + 1e-7);
	}
	KBC_gamma = gamma_nominator / (gamma_determinator + 1e-7);
}

void LBM_D2Q9_KBC::Collision(void) {
	int i = 0, j = 0, k = 0, rr = 0;

	//Collide process
#pragma omp parallel for private(i,j,rr,k)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			if (Nodes[i][j].status == "L") {
				double u_eq[D] = { 0.0, 0.0 };
				for (rr = 0; rr < D; rr++) {
					u_eq[rr] = Nodes[i][j].u[rr] + tau * (Nodes[i][j].S[rr] + Nodes[i][j].F[rr]) / Nodes[i][j].rho;
				}
				//Already got f_eq
				calculate_feq(Nodes[i][j].rho, u_eq, Nodes[i][j].feq);

				//Calculate result_fi = fi + beta * (2*s_i^{eq} - 2*s_i) + beta*gamma*(h_i^eq - h_i), where,
				//beta*gamma = 1 - (2*beta - 1) * \sum{delta_s_i * delta_h_i / feq_i} / \sum{delta_h_i * delta_h_i / feq_i}
				double KBC_Delta_s[Q], KBC_Delta_h[Q];
				double KBC_gamma = 0.0;

				Calculate_KBC_params(Nodes[i][j].rho, Nodes[i][j].f, Nodes[i][j].feq, KBC_gamma, KBC_Delta_s, KBC_Delta_h);

				KBC_gamma = (1.0 - (2.0 * beta - 1.0) * KBC_gamma) / beta;
				//KBC_gamma = 2.0;

				for (k = 0; k < Q; k++) {
					Nodes[i][j].f[k] = Nodes[i][j].f[k] - 2 * beta * KBC_Delta_s[k] - beta * KBC_gamma * KBC_Delta_h[k];
					//Nodes[i][j].f[k] = Nodes[i][j].f[k] - 2 * beta * (Nodes[i][j].f[k] - Nodes[i][j].feq[k]);
				}
			}
		}
	}
	//the end of 2-D circulation.
}


/*
Nodes[i][j].f[0] = Nodes[i][j].rho * (1.0 - KBC_T + KBC_A);

Nodes[i][j].f[1] = 0.5 * Nodes[i][j].rho * ((0.5 * (KBC_T + KBC_N) + Nodes[i][j].u[0] - KBC_Qxyy - KBC_A));
Nodes[i][j].f[3] = 0.5 * Nodes[i][j].rho * ((0.5 * (KBC_T + KBC_N) - Nodes[i][j].u[0] + KBC_Qxyy - KBC_A));

Nodes[i][j].f[2] = 0.5 * Nodes[i][j].rho * ((0.5 * (KBC_T - KBC_N) + Nodes[i][j].u[1] - KBC_Qxxy - KBC_A));
Nodes[i][j].f[4] = 0.5 * Nodes[i][j].rho * ((0.5 * (KBC_T - KBC_N) - Nodes[i][j].u[1] + KBC_Qxxy - KBC_A));

Nodes[i][j].f[5] = 0.25 * Nodes[i][j].rho * (KBC_A + KBC_PIxy + KBC_Qxyy + KBC_Qxxy);
Nodes[i][j].f[6] = 0.25 * Nodes[i][j].rho * (KBC_A - KBC_PIxy - KBC_Qxyy + KBC_Qxxy);
Nodes[i][j].f[7] = 0.25 * Nodes[i][j].rho * (KBC_A + KBC_PIxy - KBC_Qxyy - KBC_Qxxy);
Nodes[i][j].f[8] = 0.25 * Nodes[i][j].rho * (KBC_A - KBC_PIxy + KBC_Qxyy - KBC_Qxxy);
*/