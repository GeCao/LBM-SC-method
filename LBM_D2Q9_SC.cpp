#include "LBM_D2Q9_SC.h"
void LBM_D2Q9_SC::calcu_upr() {
	int i = 0, j = 0, rr = 0;
#pragma omp parallel for private(i,j,rr)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			if (Nodes[i][j].status == "L") {
				for (rr = 0; rr < D; rr++) {
					Nodes[i][j].upu[rr] = Nodes[i][j].u[rr] + (Nodes[i][j].F[rr] + Nodes[i][j].S[rr]) / (2.0 * Nodes[i][j].rho);
				}
			}
			else {
				for (rr = 0; rr < D; rr++) {
					Nodes[i][j].upu[rr] = Nodes[i][j].u[rr];
				}
			}
			//end if.
		}
	}
	//the end of 2-D circulation.
}

void LBM_D2Q9_SC::calcu_Fxy() {
	int i = 0, j = 0, rr = 0, k = 0;
	int xp = 0, yp = 0;
	double  R = 1.0, a = 12 * cs2, b = 4.0, temp_Fxy = 0.0, temp_rho = 0.0;
	double G1 = -1.0 / 3.0, psx_w = 0.0, psx_up = 0.0;
	double Tc = 0.3773 * a / (b * R);
	double TT = TT0 * Tc;

	TT = cs2 / R;
#pragma omp parallel for private(i,j,temp_rho,temp_Fxy)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			if (Nodes[i][j].status == "L") {
				temp_rho = Nodes[i][j].rho * b / 4.0;
				//temp_Fxy = cs2*(1.0 + (4.0*temp_rho - 2.0*temp_rho*temp_rho) / pow((1.0 - temp_rho), 3.0)) - a*Nodes[i][j].rho - cs2;
				temp_Fxy = R * TT * (1.0 + (4.0 * temp_rho - 2.0 * temp_rho * temp_rho) / pow((1.0 - temp_rho), 3.0)) - a * temp_rho - cs2;

				psx[i][j] = sqrt(2.0 * Nodes[i][j].rho * temp_Fxy / G1 / cs2);
				Nodes[i][j].p = Nodes[i][j].rho * cs2 + 0.5 * cs2 * G1 * psx[i][j] * psx[i][j];
			}
			//end if.
		}
	}
	//the end of 2-D circulation.
	temp_rho = rho_w * b / 4.0;
	temp_Fxy = R * TT * (1.0 + (4.0 * temp_rho - 2.0 * temp_rho * temp_rho) / pow((1.0 - temp_rho), 3.0)) - a * temp_rho - cs2;
	psx_w = sqrt(2.0 * rho_w * temp_Fxy / G1 / cs2);

	psx_up = psx_w;
	G1 = -1.0 / 3.0;
#pragma omp parallel for private(i,j,rr,k,xp,yp)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			for (rr = 0; rr < D; rr++) {
				Nodes[i][j].F[rr] = 0.0;
				Nodes[i][j].S[rr] = 0.0;
			}
			if (Nodes[i][j].status == "L") {
				for (k = 1; k < Q; k++) {
					xp = i + (int)e[k][0];
					yp = j + (int)e[k][1];
					if (xp > l) xp = 0;
					else if (xp < 0) xp = l;
					if (yp > m) yp = 0;
					else if (yp < 0) yp = m;

					if (Nodes[xp][yp].status[0] == 'B' && yp == 0) {
						//interact with solid nodes.(prepare)
						for (rr = 0; rr < D; rr++) {
							Nodes[i][j].S[rr] += w[k] * e[k][rr] * psx_w;
						}
					}
					else if (Nodes[xp][yp].status[0] == 'B' && yp == m) {
						//interact with solid nodes.(prepare)
						for (rr = 0; rr < D; rr++) {
							Nodes[i][j].S[rr] += w[k] * e[k][rr] * psx_up;
						}
					}
					else {
						//interact with fluid nodes.(prepare)
						for (rr = 0; rr < D; rr++)
						{
							Nodes[i][j].F[rr] += w[k] * e[k][rr] * psx[xp][yp];
						}
					}
					//end if.
				}
				//Final wall-fluid interaction.
				for (rr = 0; rr < D; rr++) {
					Nodes[i][j].S[rr] = -G1 * Nodes[i][j].S[rr] * psx[i][j] * c;
				}
				//Final fluid-fluid interaction.
				for (rr = 0; rr < D; rr++) {
					Nodes[i][j].F[rr] = -G1 * Nodes[i][j].F[rr] * psx[i][j] * c;
				}
			}
			//end if.
		}
	}
	//the end of 2-D circulation.
}

void LBM_D2Q9_SC::compute() {
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

void LBM_D2Q9_SC::drawing_liquid(void) {
	GLint i = 0, j = 0;
	GLdouble scale_velo = 20;

	glPushMatrix();
	glColor3f(0.0, 0.0, 1.0);
	glPointSize(4.0);
	glBegin(GL_POINTS);
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			if (Nodes[i][j].status == "L") {
				if (Nodes[i][j].rho > 0.5 * (rho_h + rho_l)) {
					glColor3f(0.0, 0.0, 1.0);
					glVertex2f(i * Lx / l, j * Ly / m);
				}
			}
		}
	}
	glEnd();
	glPopMatrix();
}