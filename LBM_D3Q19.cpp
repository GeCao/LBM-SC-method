#include "LBM_D3Q19.h"
void LBM_D3Q19::Bounceback_Boundary()
{
	int i = 0, j = 0, r = 0;
	if (up_boundary == 'B') {
		for (i = 0; i <= l; i++) {
			for (j = 0; j <= m; j++) {
				Nodes[i][j][h].status = "Bup";
			}
		}
	}
	if (down_boundary == 'B') {
		for (i = 0; i <= l; i++) {
			for (j = 0; j <= m; j++) {
				Nodes[i][j][0].status = "Bdown";
			}
		}
	}
	if (front_boundary == 'B') {
		for (i = 0; i <= l; i++) {
			for (r = 0; r <= h; r++) {
				Nodes[i][0][r].status = "Bfront";
			}
		}
	}
	if (behind_boundary == 'B') {
		for (i = 0; i <= l; i++) {
			for (r = 0; r <= h; r++) {
				Nodes[i][m][r].status = "Bbehind";
			}
		}
	}
	if (left_boundary == 'B') {
		for (j = 0; j <= m; j++) {
			for (r = 0; r <= h; r++) {
				Nodes[0][j][r].status = "Bleft";
			}
		}
	}
	if (right_boundary == 'B') {
		for (j = 0; j <= m; j++) {
			for (r = 0; r <= h; r++) {
				Nodes[l][j][r].status = "Bright";
			}
		}
	}
}

double LBM_D3Q19::Cal_error()
{
	int i = 0, j = 0, r = 0, rr = 0;
	double temp_error = 0.0, final_error = 0.0;
#pragma omp parallel for private(i,j,r,rr,temp_error)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			for (r = 0; r < h; r++) {
				if (Nodes[i][j][r].status == "L") {
					temp_error = 0.0;
					for (rr = 0; rr < D; rr++) {
						temp_error += (Nodes[i][j][r].u[rr] - velocity_before[i][j][r][rr]) * (Nodes[i][j][r].u[rr] - velocity_before[i][j][r][rr]);
					}
					final_error += sqrt(temp_error);
				}
			}
		}
	}
	return final_error;
}

void LBM_D3Q19::calculate_feq(GLdouble rh, GLdouble u[], GLdouble feq_temp[])
{
	GLdouble eu = 0.0, uv = 0.0;
	GLint k = 0, rr = 0;

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
	return;
}

void LBM_D3Q19::Record_velocity(void)
{
	//Record the velocity
	int i = 0, j = 0, r = 0, rr = 0;
#pragma omp parallel for private(i,j,r,rr)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			for (r = 0; r <= h; r++) {
				for (rr = 0; rr < D; rr++) {
					velocity_before[i][j][r][rr] = Nodes[i][j][r].u[rr];
				}
			}
		}
	}
}

void LBM_D3Q19::rebounce_f()
{
	int i = 0, j = 0, r = 0, index = 0;
	int n_of_Bounce_index = Bounce_index.size();
#pragma omp parallel for private(i,j,r,index)
	for (index = 0; index < n_of_Bounce_index; index++)
	{
		i = Bounce_index[index][0];
		j = Bounce_index[index][1];
		r = Bounce_index[index][2];
		if (Nodes[i][j][r].status[0] == 'B')
		{
			std::swap(Nodes[i][j][r].f[1], Nodes[i][j][r].f[3]);
			std::swap(Nodes[i][j][r].f[2], Nodes[i][j][r].f[4]);
			std::swap(Nodes[i][j][r].f[5], Nodes[i][j][r].f[7]);
			std::swap(Nodes[i][j][r].f[6], Nodes[i][j][r].f[8]);

			std::swap(Nodes[i][j][r].f[9], Nodes[i][j][r].f[14]);

			std::swap(Nodes[i][j][r].f[10], Nodes[i][j][r].f[17]);
			std::swap(Nodes[i][j][r].f[11], Nodes[i][j][r].f[18]);
			std::swap(Nodes[i][j][r].f[12], Nodes[i][j][r].f[15]);
			std::swap(Nodes[i][j][r].f[13], Nodes[i][j][r].f[16]);
		}
		//end if.
	}
	//The end of 3-D circulation.
}

void LBM_D3Q19::Streaming()
{
	int i = 0, j = 0, r = 0, ileft = 0, iright = 0, jfront = 0, jbehind = 0, rup = 0, rdown = 0, k = 0;

	//Streaming process
#pragma omp parallel for private(i,j,r,ileft,iright,jbehind,jfront,rdown,rup)
	for (i = 0; i <= l; i++) {
		ileft = (i > 0) ? (i - 1) : (l);
		iright = (i < l) ? (i + 1) : (0);
		for (j = 0; j <= m; j++) {
			jfront = (j > 0) ? (j - 1) : (m);
			jbehind = (j < m) ? (j + 1) : (0);
			for (r = 0; r <= h; r++) {
				rdown = (r > 0) ? (r - 1) : (h);
				rup = (r < h) ? (r + 1) : (0);

				Nodes[iright][j][r].fnew[1] = Nodes[i][j][r].f[1];
				Nodes[i][jbehind][r].fnew[2] = Nodes[i][j][r].f[2];
				Nodes[ileft][j][r].fnew[3] = Nodes[i][j][r].f[3];
				Nodes[i][jfront][r].fnew[4] = Nodes[i][j][r].f[4];

				Nodes[iright][jbehind][r].fnew[5] = Nodes[i][j][r].f[5];
				Nodes[ileft][jbehind][r].fnew[6] = Nodes[i][j][r].f[6];
				Nodes[ileft][jfront][r].fnew[7] = Nodes[i][j][r].f[7];
				Nodes[iright][jfront][r].fnew[8] = Nodes[i][j][r].f[8];

				Nodes[i][j][rup].fnew[9] = Nodes[i][j][r].f[9];

				Nodes[iright][j][rup].fnew[10] = Nodes[i][j][r].f[10];
				Nodes[i][jbehind][rup].fnew[11] = Nodes[i][j][r].f[11];
				Nodes[ileft][j][rup].fnew[12] = Nodes[i][j][r].f[12];
				Nodes[i][jfront][rup].fnew[13] = Nodes[i][j][r].f[13];

				Nodes[i][j][rdown].fnew[14] = Nodes[i][j][r].f[14];

				Nodes[iright][j][rdown].fnew[15] = Nodes[i][j][r].f[15];
				Nodes[i][jbehind][rdown].fnew[16] = Nodes[i][j][r].f[16];
				Nodes[ileft][j][rdown].fnew[17] = Nodes[i][j][r].f[17];
				Nodes[i][jfront][rdown].fnew[18] = Nodes[i][j][r].f[18];
			}
		}
	}
	//the end of 3-D circulation.
#pragma omp parallel for private(i,j,r,k)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			for (r = 0; r <= h; r++) {
				for (k = 1; k < Q; k++) {
					Nodes[i][j][r].f[k] = Nodes[i][j][r].fnew[k];
				}
			}
		}
	}
	//the end of 3-D circulation.
}

void LBM_D3Q19::macro_process()
{
	int i = 0, j = 0, r = 0, k = 0, rr = 0;
	//macroscopic process
#pragma omp parallel for private(i,j,r,rr,k)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			for (r = 0; r <= h; r++) {
				if (Nodes[i][j][r].status == "L") {
					Nodes[i][j][r].rho = 0.0;
					for (k = 0; k < Q; k++) {
						Nodes[i][j][r].rho += Nodes[i][j][r].f[k];
					}

					if (Nodes[i][j][r].rho != 0.0) {
						for (rr = 0; rr < D; rr++) {
							Nodes[i][j][r].u[rr] = 0.0;
							for (k = 0; k < Q; k++) {
								Nodes[i][j][r].u[rr] += e[k][rr] * Nodes[i][j][r].f[k] * c;
							}
						}
						for (rr = 0; rr < D; rr++) {
							Nodes[i][j][r].u[rr] /= Nodes[i][j][r].rho;
						}
					}
					//end if.
				}
				//end if.
			}
		}
	}
	//the end of 3-D circulation.
}

void LBM_D3Q19::Collision()
{
	int i = 0, j = 0, r = 0, k = 0, rr = 0;
	GLdouble u_tackle[D] = { 0.0 };

	//Collide process
#pragma omp parallel for private(i,j,r,rr,k,u_tackle)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			for (r = 0; r <= h; r++) {
				if (Nodes[i][j][r].status == "L") {
					for (rr = 0; rr < D; rr++) {
						u_tackle[rr] = Nodes[i][j][r].u[rr] + tau * (Nodes[i][j][r].S[rr] + Nodes[i][j][r].F[rr]) / Nodes[i][j][r].rho;
					}

					calculate_feq(Nodes[i][j][r].rho, u_tackle, Nodes[i][j][r].feq);
					for (k = 0; k < Q; k++) {
						Nodes[i][j][r].f[k] = Nodes[i][j][r].f[k] * (1.0 - 1.0 / tau) + Nodes[i][j][r].feq[k] / tau;
					}
				}
				//end if.
			}
		}
	}
	//the end of 3-D circulation.
}

void LBM_D3Q19::drawing_liquid(void)
{
	GLint i = 0, j = 0, r = 0;
	GLdouble scale_velo = 20;

	glColor3f(0.0, 0.0, 1.0);
	glPointSize(0.4);
	glBegin(GL_POINTS);
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			for (r = 0; r <= h; r++) {
				if (Nodes[i][j][r].status == "L") {
					if (Nodes[i][j][r].rho > 0.5 * (rho_h + rho_l)) {
						glColor3f(0.0, 0.0, 1.0);
						glPushMatrix();
						glTranslated(i * Lx / l, j * Ly / m, r * Lz / h);
						glutSolidSphere(radius, 8, 8);
						glPopMatrix();
					}
				}
			}
		}
	}
	glEnd();
}

void LBM_D3Q19::dump_file(string filename)
{
	int i = 0, j = 0, r = 0;
	FILE* fp;
	if (!(fp = fopen((filename).c_str(), "w")))
	{
		printf("写入文件失败!\n");
		exit(1);
	}
	fprintf(fp, "TITLE =\"lattice boltzmann\"\n");
	fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\" \"U\" \"V\" \"W\" \"RHO\" \n");
	fprintf(fp, "ZONE I=%d, J=%d, K=%d\n", l + 1, m + 1, h + 1);
	fprintf(fp, "F=POINT\n");
	for (r = 0; r <= h; r++)
	{
		for (j = 0; j <= m; j++)
		{
			for (i = 0; i <= l; i++)
			{
				fprintf(fp, "%d %d %d %lf %lf %lf %lf\n", i, j, r, Nodes[i][j][r].u[0], Nodes[i][j][r].u[1], Nodes[i][j][r].u[2], Nodes[i][j][r].rho);
			}
		}
	}
	fclose(fp);
}

void LBM_D3Q19::read_file(string filename)
{
	int i = 0, j = 0, r = 0;
	FILE* fp;
	char Temp_str1 = 'c', Temp_str2 = 'c';
	int Temp_i = 0, Temp_j = 0, Temp_r = 0;

	if (!(fp = fopen((filename).c_str(), "r"))) {
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
	for (r = 0; r <= h; r++) {
		for (j = 0; j <= m; j++) {
			for (i = 0; i <= l; i++) {
				fscanf(fp, "%d %d %d %lf %lf %lf %lf\n", &Temp_i, &Temp_j, &Temp_r, 
					&Nodes[i][j][r].u[0], &Nodes[i][j][r].u[1], &Nodes[i][j][r].u[2],
					&Nodes[i][j][r].rho);
			}
		}
	}
	fclose(fp);
}