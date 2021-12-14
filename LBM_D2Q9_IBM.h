#pragma once
#include "LBM_D2Q9_KBC.h"
#include "lagrangian_node.h"
class LBM_D2Q9_IBM : public LBM_D2Q9_KBC {
public:
	LBM_D2Q9_IBM() : LBM_D2Q9_KBC("VOBBPP") { //Boundary: VonNeumann, Outlet, Bounceback, Bounceback, X, X
		/*Set cylinder parameters*/
		DD = 20.0;
		cylinder_center[0] = 10 * DD; cylinder_center[1] = 10 * DD - 10;
		cylinder_area = pi * (0.25 * DD * DD);

		/*Set resolution and size of drawing pictures*/
		l = int(40 * DD); m = int(20 * DD);
		Lx = l; Ly = m;

		/*Set parameters for LBM*/
		rho_h = 1.0;
		Vmax = 0.05;
		Re = 200.0;
		visc = Vmax * DD / Re; //初始粘度 | The initial viscosity of particles.
		tau = visc / cs2 + 0.5;
		beta = 1.0 / 2.0 / tau; //弛豫时间

		Ma = Vmax / cs2;
	}

	virtual void initialization() {
		std::cout << "[Re number]: " << Re << std::endl;
		//Allocate Nodes space:
		Nodes = allocate_space_2D<lattice_node_D2Q9>(l + 1, m + 1);
		Vorticity = allocate_space_2D<double>(l + 1, m + 1);
		

		Vonneumann_Boundary();
		Outlet_Boundary();
		Bounceback_Boundary();

		int i = 0, j = 0, k = 0, rr = 0;
		std::vector<int> temp_bounce(D, 0);

		for (i = 0; i <= l; i++) {
			for (j = 0; j <= m; j++) {
				Vorticity[i][j] = 0.0;
				for (rr = 0; rr < D; rr++) {
					Nodes[i][j].F[rr] = 0.0;
					Nodes[i][j].S[rr] = 0.0;
					Nodes[i][j].upu[rr] = 0.0;
				}
				Nodes[i][j].p = p0;
				double dist2_cylinder = pow(cylinder_center[0] - i, 2) + pow(cylinder_center[1] - j, 2);
				/*
				if(dist2_cylinder > DD * DD / 4.0 && dist2_cylinder <= (DD) * (DD)){
					Nodes[i][j].u[0] = Vmax * tanh(80.0 * (-j * 1.0 / (m + 1) + 0.5));
					Nodes[i][j].u[1] = Vmax * 0.01 * sin(2.0 * pi * (i * 1.0 / (l + 1) + 0.25));
				}*/

				if (dist2_cylinder <= DD * DD / 4.0) {
					Nodes[i][j].status = "BCylinder";
				}

				if (Nodes[i][j].status == "L") {
					Nodes[i][j].rho = rho_h;
				}
				else if (Nodes[i][j].status[0] == 'B') {
					Nodes[i][j].rho = rho_w;
					temp_bounce[0] = i;
					temp_bounce[1] = j;
					Bounce_index.push_back(temp_bounce);
				}
				else {
					Nodes[i][j].rho = rho_h;
				}

				calculate_feq(Nodes[i][j].rho, Nodes[i][j].u, Nodes[i][j].feq);
				for (k = 0; k < Q; k++) {
					Nodes[i][j].f[k] = Nodes[i][j].feq[k];
				}
			}
		}
		//the end of 3-D circulation.

		/*Set points around the cylinder*/
		int population_idx_of_direction[3][3] = { {7, 3, 6}, {4, 0, 2}, {8, 1, 5} };
		for (i = 0; i <= l; i++) {
			for (j = 0; j <= m; j++) {
				bool detected_cylinder = false;
				if (Nodes[i][j].status[0] == 'L') {
					int x_axis[3] = { (i > 0) ? i - 1 : l, i, (i < l) ? i + 1 : 0 };
					int y_axis[3] = { (j > 0) ? j - 1 : m, j, (j < m) ? j + 1 : 0 };
					for (int ii = 0; ii < 3; ii++) {
						for (int jj = 0; jj < 3; jj++) {
							if (Nodes[x_axis[ii]][y_axis[jj]].status == "BCylinder") {
								double xb = i, yb = j;
								double xc = cylinder_center[0], yc = cylinder_center[1];
								double xs = x_axis[ii], ys = y_axis[jj];

								double quad_a = (xs - xb) * (xs - xb) + (ys - yb) * (ys - yb);
								double quad_b = 2 * (xs - xb) * (xb - xc) + 2 * (ys - yb) * (yb - yc);
								double quad_c = (xb - xc) * (xb - xc) + (yb - yc) * (yb - yc) - DD * DD * 0.25;
								double t_ij = (-quad_b + sqrt(quad_b * quad_b - 4.0 * quad_a * quad_c)) /
									(2.0 * quad_a);
								if (t_ij > 1 + 1e-7 || t_ij < 0 - 1e-7) {
									t_ij = (-quad_b - sqrt(quad_b * quad_b - 4.0 * quad_a * quad_c)) /
										(2.0 * quad_a);
								}
								if (t_ij > 1 + 1e-7 || t_ij < 0 - 1e-7) {
									std::cout << "The computation of cylinder might wrong: (xs, ys) = ";
									std::cout << xs << ", " << ys << " || (xb, yb) = " << xb << ", " << yb << std::endl;
								}

								if (!detected_cylinder) {
									idx_near_cylinder_boundary.push_back(std::vector<int>{i, j});
									idx_of_missing_props.push_back(
										std::vector<int>{population_idx_of_direction[ii][jj]});
									
									tij_near_cylinder.push_back(std::vector<double>{t_ij});
									detected_cylinder = true;
									continue;
								}
								else {
									idx_of_missing_props[idx_of_missing_props.size() - 1].push_back(
										population_idx_of_direction[ii][jj]
									);
									tij_near_cylinder[tij_near_cylinder.size() - 1].push_back(t_ij);
								}
							}
						}
					}
				}
			}
		}

		FILE* fp_log_aerodynamics_forces;
		if (!(fp_log_aerodynamics_forces = fopen("log_aerodynamics_forces.txt", "w"))){
			printf("写入文件失败!\n");
			exit(-1);
		}
		fprintf(fp_log_aerodynamics_forces, "Aerodynamics Forces: %dD\n", D);
		fprintf(fp_log_aerodynamics_forces, "Step Force_x Force_y Torque_z\n");
		fclose(fp_log_aerodynamics_forces);
	}

	virtual ~LBM_D2Q9_IBM() {
		delete_space_2D(Vorticity, l + 1, m + 1);
	}

	virtual void drawing_liquid(void);

	virtual void dump_file(string filename);

	virtual void compute() {
		//Record the velocity
		Record_velocity();

		//Collision process
		Collision();

		//Streaming process
		Streaming();
		post_streaming();

		//macro process
		macro_process();
		rebounce_f();
		for (int i = 0; i <= l; i++) { //Free-stream boundary for up/down boundary
			int j = 0;
			if (Nodes[i][j].status[0] == 'B'){
				std::swap(Nodes[i][j].f[1], Nodes[i][j].f[3]);
				std::swap(Nodes[i][j].f[2], Nodes[i][j].f[4]);
				std::swap(Nodes[i][j].f[5], Nodes[i][j].f[8]);
				std::swap(Nodes[i][j].f[6], Nodes[i][j].f[7]);
			}
			j = m;
			if (Nodes[i][j].status[0] == 'B') {
				std::swap(Nodes[i][j].f[1], Nodes[i][j].f[3]);
				std::swap(Nodes[i][j].f[2], Nodes[i][j].f[4]);
				std::swap(Nodes[i][j].f[5], Nodes[i][j].f[8]);
				std::swap(Nodes[i][j].f[6], Nodes[i][j].f[7]);
			}
		}
		Vonneumann_f(Vmax); //left, up, right, down (Corresponding to the streaming index)
		outlet_f();
	}

	virtual void post_streaming();

	
	virtual void Collision(void) {
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

					for (k = 0; k < Q; k++) {
						//Nodes[i][j].f[k] = Nodes[i][j].f[k] - 2 * beta * KBC_Delta_s[k] - beta * KBC_gamma * KBC_Delta_h[k];
						Nodes[i][j].f[k] = Nodes[i][j].f[k] - 2 * beta * (Nodes[i][j].f[k] - Nodes[i][j].feq[k]);
					}
				}
			}
		}
		//the end of 2-D circulation.
	}
	
private:
	bool is_in_array(int k, std::vector<int> array) {
		int n = array.size();
		for (int i = 0; i < n; i++) {
			if (array[i] == k) { return true; }
		}
		return false;
	}

	double DD; //The diameter of cylinder.
	double cylinder_center[2]; //The center of cylinder.
	double cylinder_area;//The area of cylinder.
	
	std::vector<std::vector<int> > idx_near_cylinder_boundary; // nb X 2
	std::vector<std::vector<int> > idx_of_missing_props; // nb X *
	std::vector<std::vector<double> > tij_near_cylinder; // nb X *
};

