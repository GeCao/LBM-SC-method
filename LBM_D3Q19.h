#pragma once
#include "LBM.h"
#include <vector>
#include "lattice_node.h"
#include "marching_cube.h"
class LBM_D3Q19 : public LBM {
public:
	LBM_D3Q19()
	{
		Bounceback_Boundary();

		GLint i = 0, j = 0, r = 0, k = 0, rr = 0;
		vector<GLint> temp_bounce(D, 0);
		for (i = 0; i <= l; i++){
			for (j = 0; j <= m; j++){
				for (r = 0; r <= h; r++){
					for (rr = 0; rr < D; rr++){
						Nodes[i][j][r].u[rr] = 0.0;
						Nodes[i][j][r].F[rr] = 0.0;
						Nodes[i][j][r].S[rr] = 0.0;
						Nodes[i][j][r].upu[rr] = 0.0;
					}
					Nodes[i][j][r].p = p0;

					if (Nodes[i][j][r].status == "L"){
						Nodes[i][j][r].rho = rho_l;
						/*
						if ((i - 0.5*l)*(i - 0.5*l) + (j - 0.5*m)*(j - 0.5*m) + (r - 0.5*h)*(r - 0.5*h) < RR*RR)
						{
							Nodes[i][j][r].rho = rho_h;
							Nodes[i][j][r].u[2] = UU;
						}*/
						if ((i > 0.2 * l && i < 0.8 * l) && (j > 0.2 * m && j < 0.8 * m) && (r > 0.4 * h && r < 0.6 * h)){
							Nodes[i][j][r].rho = rho_h;
							Nodes[i][j][r].u[2] = UU;
						}
					}
					else
					{
						Nodes[i][j][r].rho = rho_w;
						temp_bounce[0] = i;
						temp_bounce[1] = j;
						temp_bounce[2] = r;
						Bounce_index.push_back(temp_bounce);
					}

					calculate_feq(Nodes[i][j][r].rho, Nodes[i][j][r].u, Nodes[i][j][r].feq);
					for (k = 0; k < Q; k++)
					{
						Nodes[i][j][r].f[k] = Nodes[i][j][r].feq[k];
					}
				}
			}
		}
		//the end of 3-D circulation.
	}

	/*
	* Set the Boundary condition.
	*/
	virtual void Bounceback_Boundary();

	virtual void calculate_feq(GLdouble rh, GLdouble u[], GLdouble feq_temp[]);

	virtual double Cal_error();

	virtual void Record_velocity(void);

	virtual void rebounce_f();

	virtual void Streaming();

	virtual void macro_process();

	virtual void Collision();

	virtual void dump_file(string filenam="data");
	virtual void read_file(string filenam = "data");

	virtual void drawing_liquid(void);

	virtual GLint marching_cube(void) {
		throw "Marching cube algorithm: Did not implemented!";
	}

protected:
	lattice_node Nodes[l + 1][m + 1][h + 1];
	GLdouble velocity_before[l + 1][m + 1][h + 1][D];
	std::vector<std::vector<GLint> > Bounce_index;
};