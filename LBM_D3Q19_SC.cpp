#include "LBM_D3Q19_SC.h"
void LBM_D3Q19_SC::calcu_upr()
{
	int i = 0, j = 0, r = 0, rr = 0;
#pragma omp parallel for private(i,j,r,rr)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			for (r = 0; r <= h; r++) {
				if (Nodes[i][j][r].status == "L") {
					for (rr = 0; rr < D; rr++) {
						Nodes[i][j][r].upu[rr] = Nodes[i][j][r].u[rr] + (Nodes[i][j][r].F[rr] + Nodes[i][j][r].S[rr]) / (2.0 * Nodes[i][j][r].rho);
					}
				}
				else {
					for (rr = 0; rr < D; rr++) {
						Nodes[i][j][r].upu[rr] = Nodes[i][j][r].u[rr];
					}
				}
				//end if.
			}
		}
	}
	//the end of 3-D circulation.
}

void LBM_D3Q19_SC::calcu_Fxy()
{
	int i = 0, j = 0, r = 0, rr = 0, k = 0;
	int xp = 0, yp = 0, zp = 0;
	double  R = 1.0, a = 12 * cs2, b = 4.0, temp_Fxy = 0.0, temp_rho = 0.0;
	double G1 = -1.0 / 3.0, psx_w = 0.0;
	double Tc = 0.3773 * a / (b * R);
	double TT = TT0 * Tc;
	double temp_sum[D] = { 0.0, 0.0 };

	TT = cs2 / R;
#pragma omp parallel for private(i,j,r,temp_rho,temp_Fxy)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			for (r = 0; r <= h; r++) {
				if (Nodes[i][j][r].status == "L") {
					temp_rho = Nodes[i][j][r].rho * b / 4.0;
					temp_Fxy = R * TT * (1.0 + (4.0 * temp_rho - 2.0 * temp_rho * temp_rho) / pow((1.0 - temp_rho), 3.0)) - a * temp_rho - cs2;

					psx[i][j][r] = sqrt(fabs(2.0 * temp_rho * temp_Fxy / G1 / cs2));
					Nodes[i][j][r].p = temp_rho * cs2 + 0.5 * cs2 * G1 * psx[i][j][r] * psx[i][j][r];
				}
				//end if.
			}
		}
	}
	//the end of 3-D circulation.
	temp_rho = rho_w;
	temp_Fxy = R * TT * (1.0 + (4.0 * temp_rho - 2.0 * temp_rho * temp_rho) / pow((1.0 - temp_rho), 3.0)) - a * temp_rho - cs2;
	psx_w = sqrt(fabs(2.0 * rho_w * temp_Fxy / G1 / cs2));

	G1 = -1.0 / 3.0;
#pragma omp parallel for private(i,j,r,rr,k,xp,yp,zp,temp_sum)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			for (r = 0; r <= h; r++) {
				for (rr = 0; rr < D; rr++) {
					Nodes[i][j][r].F[rr] = 0.0;
					Nodes[i][j][r].S[rr] = 0.0;
					temp_sum[rr] = 0.0;
				}
				if (Nodes[i][j][r].status == "L") {
					for (k = 1; k < Q; k++) {
						xp = i + (int)e[k][0];
						yp = j + (int)e[k][1];
						zp = r + (int)e[k][2];
						if (xp > l) xp = 0;
						else if (xp < 0) xp = l;
						if (yp > m) yp = 0;
						else if (yp < 0) yp = m;
						if (zp > h) zp = 0;
						else if (zp < 0) zp = h;

						if (Nodes[xp][yp][zp].status[0] == 'B') {
							//interact with solid nodes.(prepare)
							for (rr = 0; rr < D; rr++) {
								temp_sum[rr] += w[k] * e[k][rr];
							}
						}
						else {
							//interact with fluid nodes.(prepare)
							for (rr = 0; rr < D; rr++) {
								Nodes[i][j][r].F[rr] += w[k] * e[k][rr] * psx[xp][yp][zp];
							}
						}
						//end if.
					}
					//Final wall-fluid interaction.
					for (rr = 0; rr < D; rr++) {
						Nodes[i][j][r].S[rr] = -G1 * temp_sum[rr] * psx[i][j][r] * psx_w * c;
					}
					//Final fluid-fluid interaction.
					for (rr = 0; rr < D; rr++) {
						Nodes[i][j][r].F[rr] *= -G1 * psx[i][j][r] * c;
					}
				}
				//end if.
			}
		}
	}
	//the end of 3-D circulation.
}

void LBM_D3Q19_SC::compute() {
	//Record the velocity
	Record_velocity();

	//Streaming process
	Streaming();

	//macro process
	macro_process();

	calcu_upr();

	calcu_Fxy();

	rebounce_f();
	//Collision process
	Collision();
}

void LBM_D3Q19_SC::drawing_liquid(void)
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