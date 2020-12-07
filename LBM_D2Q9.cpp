#include "LBM_D2Q9.h"
void LBM_D2Q9::Bounceback_Boundary() {
	GLint i = 0, j = 0;
	if (up_boundary == 'B') {
		for (i = 0; i <= l; i++) {
			Nodes[i][m].status = "Bup";
		}
	}
	if (down_boundary == 'B') {
		for (i = 0; i <= l; i++) {
			Nodes[i][0].status = "Bdown";
		}
	}
	if (left_boundary == 'B') {
		for (j = 0; j <= m; j++) {
			Nodes[0][j].status = "Bleft";
		}
	}
	if (right_boundary == 'B') {
		for (j = 0; j <= m; j++) {
			Nodes[l][j].status = "Bright";
		}
	}
}

double LBM_D2Q9::Cal_error() {
	int i = 0, j = 0, rr = 0;
	double temp_error = 0.0, final_error = 0.0;
#pragma omp parallel for private(i,j,rr,temp_error)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			if (Nodes[i][j].status == "L") {
				temp_error = 0.0;
				for (rr = 0; rr < D; rr++) {
					temp_error += (Nodes[i][j].u[rr] - velocity_before[i][j][rr]) * (Nodes[i][j].u[rr] - velocity_before[i][j][rr]);
				}
				final_error += sqrt(temp_error);
			}
		}
	}
	return final_error;
}

void LBM_D2Q9::calculate_feq(GLdouble rh, GLdouble u[], GLdouble feq_temp[])
{
	GLdouble eu = 0.0, uv = 0.0;
	GLint k = 0;

	/*
	for (GLint rr = 0; rr < D; rr++) {
		uv += (u[rr] * u[rr]);
	}
#pragma omp parallel for private(k,rr,eu)
	for (k = 0; k < Q; k++) {
		eu = 0.0;
		for (rr = 0; rr < D; rr++) {
			eu += (e[k][rr] * u[rr]);
		}
		feq_temp[k] = w[k] * rh * (1.0 + eu / cs2 + 0.5 * eu * eu / cs2 / cs2 - 0.5 * uv / cs2);
	}
	*/
#pragma omp parallel for private(k)
	for (k = 0; k < Q; k++) {
		feq_temp[k] = rh * w[k] * (2.0 - sqrt(1.0 + 3.0 * u[0] * u[0])) * \
			(2.0 - sqrt(1.0 + 3.0 * u[1] * u[1])) * \
			pow((2.0 * u[0] + sqrt(1.0 + 3.0 * u[0] * u[0])) / (1 - u[0]), e[k][0]) * \
			pow((2.0 * u[1] + sqrt(1.0 + 3.0 * u[1] * u[1])) / (1 - u[1]), e[k][1]);
	}
}

void LBM_D2Q9::Record_velocity(void) {
	//Record the velocity
	int i = 0, j = 0, rr = 0;
#pragma omp parallel for private(i,j,rr)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			for (rr = 0; rr < D; rr++) {
				velocity_before[i][j][rr] = Nodes[i][j].u[rr];
			}
		}
	}
}

void LBM_D2Q9::rebounce_f() {
	int i = 0, j = 0, index = 0;
	int n_of_Bounce_index = Bounce_index.size();
#pragma omp parallel for private(i,j,index)
	for (index = 0; index < n_of_Bounce_index; index++)
	{
		i = Bounce_index[index][0];
		j = Bounce_index[index][1];
		if (Nodes[i][j].status[0] == 'B')
		{
			std::swap(Nodes[i][j].f[1], Nodes[i][j].f[3]);
			std::swap(Nodes[i][j].f[2], Nodes[i][j].f[4]);
			std::swap(Nodes[i][j].f[5], Nodes[i][j].f[7]);
			std::swap(Nodes[i][j].f[6], Nodes[i][j].f[8]);
		}
		//end if.
	}
	//The end of 3-D circulation.
	return;
}

void LBM_D2Q9::Streaming() {
	int i = 0, j = 0, ileft = 0, iright = 0, jup = 0, jdown = 0, k = 0;

	//Streaming process
#pragma omp parallel for private(i,j,ileft,iright,jdown,jup)
	for (i = 0; i <= l; i++) {
		ileft = (i > 0) ? (i - 1) : (l);
		iright = (i < l) ? (i + 1) : (0);
		for (j = 0; j <= m; j++) {
			jdown = (j > 0) ? (j - 1) : (m);
			jup = (j < m) ? (j + 1) : (0);

			Nodes[iright][j].fnew[1] = Nodes[i][j].f[1];
			Nodes[i][jup].fnew[2] = Nodes[i][j].f[2];
			Nodes[ileft][j].fnew[3] = Nodes[i][j].f[3];
			Nodes[i][jdown].fnew[4] = Nodes[i][j].f[4];

			Nodes[iright][jup].fnew[5] = Nodes[i][j].f[5];
			Nodes[ileft][jup].fnew[6] = Nodes[i][j].f[6];
			Nodes[ileft][jdown].fnew[7] = Nodes[i][j].f[7];
			Nodes[iright][jdown].fnew[8] = Nodes[i][j].f[8];
		}
	}
	//the end of 2-D circulation.
#pragma omp parallel for private(i,j,k)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			for (k = 1; k < Q; k++) {
				Nodes[i][j].f[k] = Nodes[i][j].fnew[k];
			}
		}
	}
	//the end of 2-D circulation.
}

void LBM_D2Q9::macro_process() {
	GLint i = 0, j = 0, k = 0, rr = 0;
	//macroscopic process
#pragma omp parallel for private(i,j,rr,k)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			if (Nodes[i][j].status == "L") {
				Nodes[i][j].rho = 0.0;
				for (k = 0; k < Q; k++) {
					Nodes[i][j].rho += Nodes[i][j].f[k];
				}

				for (rr = 0; rr < D; rr++) {
					Nodes[i][j].u[rr] = 0.0;
					for (k = 0; k < Q; k++) {
						Nodes[i][j].u[rr] += e[k][rr] * Nodes[i][j].f[k];
					}
					Nodes[i][j].u[rr] *= (c / Nodes[i][j].rho);
				}
			}
			//end if.
		}
	}
	//the end of 2-D circulation.
}

void LBM_D2Q9::Collision(void) {
	int i = 0, j = 0, k = 0, rr = 0;
	GLdouble u_tackle[D] = { 0.0, 0.0 };

	//Collide process
#pragma omp parallel for private(i,j,rr,k,u_tackle)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			if (Nodes[i][j].status == "L") {
				for (rr = 0; rr < D; rr++) {
					u_tackle[rr] = Nodes[i][j].u[rr] + tau * (Nodes[i][j].S[rr] + Nodes[i][j].F[rr]) / Nodes[i][j].rho;
				}
				calculate_feq(Nodes[i][j].rho, u_tackle, Nodes[i][j].feq);
				for (k = 0; k < Q; k++) {
					Nodes[i][j].f[k] = Nodes[i][j].f[k] * (1.0 - 1.0 / tau) + Nodes[i][j].feq[k] / tau;
				}
			}
		}
	}
	//the end of 2-D circulation.
}

void LBM_D2Q9::drawing_liquid(void) {
	GLint i = 0, j = 0;
	GLdouble scale_velo = 10;

	glPushMatrix();
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			if (Nodes[i][j].status == "L") {
				glVertex2f(i * Lx / l, j * Ly / m);
				glVertex2f(i * Lx / l + Nodes[i][j].u[0] * scale_velo, j * Ly / m + Nodes[i][j].u[1] * scale_velo);
			}
		}
	}
	glEnd();
	glPopMatrix();
}

void LBM_D2Q9::dump_file(string filename) {
	int i = 0, j = 0;
	FILE* fp;
	if (!(fp = fopen(filename.c_str(), "w")))
	{
		printf("写入文件失败!\n");
		exit(1);
	}
	fprintf(fp, "TITLE =\"lattice boltzmann\"\n");
	fprintf(fp, "VARIABLES = \"X\" \"Y\" \"RHO\" \"U\" \"V\" \n");
	fprintf(fp, "ZONE I=%d, J=%d\n", l + 1, m + 1);
	fprintf(fp, "F=POINT\n");
	for (j = 0; j <= m; j++) {
		for (i = 0; i <= l; i++) {
			fprintf(fp, "%d %d %lf %lf %lf\n", i, j, Nodes[i][j].rho, Nodes[i][j].u[0], Nodes[i][j].u[1]);
		}
	}
	fclose(fp);
}

void LBM_D2Q9::read_file(string filename) {
	int i = 0, j = 0;
	FILE* fp;
	char Temp_str1 = 'c', Temp_str2 = 'c';
	int Temp_i = 0, Temp_j = 0;

	if (!(fp = fopen(filename.c_str(), "r"))) {
		printf("读入文件失败!\n");
		exit(1);
	}
	while (true) {
		Temp_str2 = Temp_str1;
		fscanf(fp, "%c", &Temp_str1);
		if (Temp_str1 == 'T' && Temp_str2 == 'N') {
			break;
		}
	}
	for (j = 0; j <= m; j++) {
		for (i = 0; i <= l; i++) {
			fscanf(fp, "%d %d %lf %lf %lf\n", &Temp_i, &Temp_j, &Nodes[i][j].rho, &Nodes[i][j].u[0], &Nodes[i][j].u[1]);
		}
	}
	fclose(fp);
}