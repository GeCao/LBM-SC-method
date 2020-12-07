#pragma once
#include "LBM_D2Q9.h"
#include <vector>
#include "lattice_node.h"
#include "marching_cube.h"
class LBM_D2Q9_SC : public LBM_D2Q9 {
public:
	LBM_D2Q9_SC()
	{
		Bounceback_Boundary();

		GLint i = 0, j = 0, k = 0, rr = 0;
		std::vector<GLint> temp_bounce(D, 0);
		for (i = 0; i <= l; i++) {
			for (j = 0; j <= m; j++) {
				for (rr = 0; rr < D; rr++) {
					Nodes[i][j].u[rr] = 0.0;
					Nodes[i][j].F[rr] = 0.0;
					Nodes[i][j].S[rr] = 0.0;
					Nodes[i][j].upu[rr] = 0.0;
				}
				Nodes[i][j].p = p0;
				psx[i][j] = 0.0;

				if (Nodes[i][j].status == "L") {
					Nodes[i][j].rho = rho_l;

					if ((i > 0.2 * l && i < 0.8 * l) && (j > 0.2 * m && j < 0.8 * m)) {
						Nodes[i][j].rho = rho_h;
						Nodes[i][j].u[2] = UU;
					}
				}
				else
				{
					Nodes[i][j].rho = rho_w;
					temp_bounce[0] = i;
					temp_bounce[1] = j;
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

	/*
	* Being used only in Shan-Chen model
	*/
	virtual void calcu_upr();

	/*
	* Being used only in Shan-Chen model
	*/
	virtual void calcu_Fxy();

	virtual void compute();

	virtual void drawing_liquid(void);

	virtual GLint marching_cube(void) {
		throw "Marching cube algorithm: Did not implemented!";
	}

private:
	GLdouble psx[l + 1][m + 1];
};