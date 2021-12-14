#pragma once
#include "LBM_D3Q19.h"
#include <vector>
#include "lattice_node.h"
#include "marching_cube.h"
#ifdef use_HCZ == 1
class LBM_D3Q19_HCZ : public LBM_D3Q19 {
public:
	LBM_D3Q19_HCZ() : LBM_D3Q19("PPBBPP") {
		l = 40; m = 40; h = 40;
		Lx = l; Ly = m; Lz = h;


		thetaAVE = 120.0;
		thetaR = thetaAVE - 0.0; thetaA = thetaAVE + 0.0;

		b = 4.0;
		a = 12.0 * cs2;
		Kappa = 0.1;
		RT = cs2;
		UU = 0.0;
		DD = 30.0;
		rho_h = 0.251;
		rho_l = 0.024;
		psi_max = 0.251;
		psi_min = 0.024;
		tau_f = 0.7;
		tau_g = 0.7;
		gforce = 0.0;
		Re = UU * DD / (cs2 * tau_f);//Reynolds number
		We = rho_h * DD * UU * UU / (Kappa * 0.0011 / 0.1);
		visc = cs2 * (tau_f - 0.5); //初始粘度 | The initial viscosity of particles.

		BUOYANCY = 0;
	}

	virtual void initialization() {
		psx = allocate_space_3D<double>(l + 1, m + 1, h + 1);
		virtue_down = allocate_space_2D<double>(l + 1, m + 1);
		virtue_up = allocate_space_2D<double>(l + 1, m + 1);
		Nodes = allocate_space_3D<lattice_node_D3Q19_HCZ>(l + 1, m + 1, h + 1);

		int i = 0, j = 0, r = 0, k = 0, rr = 0;
		double temp_rho = 0.0;

		Bounceback_Boundary();

#pragma omp parallel for private(i,j,r,rr,k,temp_rho)
		for (i = 0; i <= l; i++){
			for (j = 0; j <= m; j++){
				for (r = 0; r <= h; r++){
					Nodes[i][j][r].rh = psi_min;
					Nodes[i][j][r].rho = rho_l;

					if (Nodes[i][j][r].status == "L")
					{
						if ((i - l / 2) * (i - l / 2) + (j - m / 2) * (j - m / 2) + (r) * (r) < 0.25 * DD * DD) {
							Nodes[i][j][r].rh = psi_max;
							Nodes[i][j][r].rho = rho_h;
						}
					}

					for (rr = 0; rr < D; rr++)
					{
						Nodes[i][j][r].u[rr] = 0.0;
						Nodes[i][j][r].Force[rr] = 0.0;
					}

					temp_rho = b * Nodes[i][j][r].rho / 4.0;
					Nodes[i][j][r].p = Nodes[i][j][r].rho * RT * temp_rho * (4.0 - 2.0 * temp_rho) / pow((1 - temp_rho), 3) - a * Nodes[i][j][r].rho * Nodes[i][j][r].rho + Nodes[i][j][r].rho * RT;

					calculate_feq_geq(Nodes[i][j][r].rh, Nodes[i][j][r].rho, Nodes[i][j][r].p, Nodes[i][j][r].u, Nodes[i][j][r].feq, Nodes[i][j][r].geq);
					for (k = 0; k < Q; k++)
					{
						Nodes[i][j][r].f[k] = Nodes[i][j][r].feq[k];
						Nodes[i][j][r].g[k] = Nodes[i][j][r].geq[k];
					}
				}
			}
		}
		//the end of 3-D circulation.
		return;
	}

	virtual void Bounceback_Boundary();

	virtual void Dirichlet_Boundary();

	virtual void calculate_feq_geq(GLdouble rh, GLdouble rho, GLdouble p, GLdouble u[], GLdouble feq_temp[], GLdouble geq_temp[]);

	virtual void rebounce_g();

	virtual void Streaming();

	virtual void macro_process();
	virtual void geten();

	virtual void Collision();

	virtual void drawing_liquid();
	virtual void dump_file(string filename);

	virtual void compute() {
		//Record the velocity
		Record_velocity();

		//Streaming process
		Streaming();

		//macro process
		macro_process();
		geten();

		//Collision process
		rebounce_f();
		rebounce_g();
		Collision();
	}
private:
	double*** psx; // (l + 1)*(m + 1)*(h + 1)
	double** virtue_down; // (l + 1)*(m + 1)
	double** virtue_up; // (l + 1)*(m + 1)

	double thetaAVE;
	double thetaR, thetaA;

	double b;
	double a;
	double Kappa;
	double RT;
	double UU;
	double DD;
	double rho_h;
	double rho_l;
	double psi_max;
	double psi_min;
	double tau_f;
	double tau_g;
	double gforce;
	double We;
	int BUOYANCY;//该值设置为1时添加重力，置为0时不添加重力。
};
#endif // 