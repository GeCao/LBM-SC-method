#pragma once
#include "LBM_D2Q9.h"
#include <vector>
class LBM_D2Q9_KBC : public LBM_D2Q9 {
public:
	LBM_D2Q9_KBC()
	{
		Bounceback_Boundary();

		GLint i = 0, j = 0, k = 0, rr = 0;
		std::vector<GLint> temp_bounce(D, 0);
		for (i = 0; i <= l; i++) {
			for (j = 0; j <= m; j++) {
				for (rr = 0; rr < D; rr++) {
					Nodes[i][j].F[rr] = 0.0;
					Nodes[i][j].S[rr] = 0.0;
					Nodes[i][j].upu[rr] = 0.0;
				}
				Nodes[i][j].p = p0;

				if (Nodes[i][j].status == "L") {
					Nodes[i][j].u[1] = KBC_sigma * Vmax * sin(2.0 * pi * (i / (l + 1) + 0.25));
					if (j <= m / 2) {
						Nodes[i][j].u[0] = Vmax * tanh(KBC_kappa * (j / (m + 1) - 0.25));
					}
					else {
						Nodes[i][j].u[0] = Vmax * tanh(KBC_kappa * (-j / (m + 1) + 0.75));
					}
					Nodes[i][j].rho = rho_h;
				}
				else
				{
					Nodes[i][j].rho = rho_w;
					temp_bounce[0] = 0.0;
					temp_bounce[1] = 0.0;
					Bounce_index.push_back(temp_bounce);
				}

				calculate_feq(Nodes[i][j].rho, Nodes[i][j].u, Nodes[i][j].feq);
				for (k = 0; k < Q; k++)
				{
					Nodes[i][j].f[k] = Nodes[i][j].feq[k];
				}
			}
		}
		//the end of 3-D circulation.
	}

	virtual void Calculate_KBC_params(double rho, double* f, double* feq, double& KBC_gamma, double KBC_Delta_s[], double KBC_Delta_h[]);

	virtual void Collision(void);

	virtual void drawing_liquid(void);

	virtual void dump_file();

	virtual GLint marching_cube(void) {
		throw "Marching cube algorithm: Did not implemented!";
	}
private:
	GLdouble Vorticity[l + 1][m + 1];
};