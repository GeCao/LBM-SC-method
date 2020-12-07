#pragma once
#include "LBM.h"
#include <vector>
#include "lattice_node.h"
#include "marching_cube.h"
class LBM_D2Q9 : public LBM {
public:
	LBM_D2Q9()
	{
		Bounceback_Boundary();

		//Taylor-Green vortex
		GLint i = 0, j = 0, k = 0, rr = 0;
		std::vector<GLint> temp_bounce(D, 0);

		Analytical_solution(0.0);

		for (i = 0; i <= l; i++) {
			for (j = 0; j <= m; j++) {
				for (rr = 0; rr < D; rr++) {
					//Nodes[i][j].u[rr] = 0.0;
					Nodes[i][j].F[rr] = 0.0;
					Nodes[i][j].S[rr] = 0.0;
					Nodes[i][j].upu[rr] = 0.0;
				}
				Nodes[i][j].p = p0;

				if (Nodes[i][j].status[0] == 'B') {
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

	void Analytical_solution(double t) {
		int i = 0, j = 0;
#pragma omp parallel for private(i,j)
		for (i = 0; i <= l; i++) {
			for (j = 0; j <= m; j++) {
				Nodes[i][j].u[0] = -Vmax * K_y / sqrt(K_x * K_x + K_y * K_y) * exp(-visc * K * K * t) * sin(K_y * j) * cos(K_x * i);
				Nodes[i][j].u[1] = Vmax * K_x / sqrt(K_x * K_x + K_y * K_y) * exp(-visc * K * K * t) * sin(K_x * i) * cos(K_y * j);
				Nodes[i][j].rho = 1 - Ma * Ma / 2.0 / K / K * (K_y * K_y * cos(2.0 * K_x * i) + K_x * K_x * cos(2.0 * K_y * j));
			}
		}
	}

	double Analytical_error(double t) {
		int i = 0, j = 0;
		double final_error = 0.0;
#pragma omp parallel for private(i,j)
		for (i = 0; i <= l; i++) {
			for (j = 0; j <= m; j++) {
				double uxx = -Vmax * K_y / sqrt(K_x * K_x + K_y * K_y) * exp(-visc * K * K * t) * sin(K_y * j) * cos(K_x * i);
				double uyy = Vmax * K_x / sqrt(K_x * K_x + K_y * K_y) * exp(-visc * K * K * t) * sin(K_x * i) * cos(K_y * j);
				double rhoxxyy = 1 - Ma * Ma / 2.0 / K / K * (K_y * K_y * cos(2.0 * K_x * i) + K_x * K_x * cos(2.0 * K_y * j));
				final_error += pow(uxx - Nodes[i][j].u[0], 2) + pow(uyy - Nodes[i][j].u[1], 2);
			}
		}
		final_error = sqrt(final_error);
		return final_error;
	}

	virtual void Bounceback_Boundary();

	virtual void calculate_feq(GLdouble rh, GLdouble u[], GLdouble feq_temp[]);

	virtual double Cal_error();

	virtual void Record_velocity(void);

	virtual void rebounce_f();

	virtual void Streaming();

	virtual void macro_process();

	virtual void Collision(void);

	virtual void dump_file(string filenam = "data");
	virtual void read_file(string filenam = "data");

	virtual void drawing_liquid(void);

	void drawing_analytical_liquid(double t) {
		GLint i = 0, j = 0;
		GLdouble scale_velo = 10;

		glPushMatrix();
		glColor3f(0.0, 0.0, 1.0);
		glBegin(GL_LINES);
		for (i = 0; i <= l; i++) {
			for (j = 0; j <= m; j++) {
				if (Nodes[i][j].status == "L") {
					double uxx = -Vmax * K_y / sqrt(K_x * K_x + K_y * K_y) * exp(-visc * K * K * t) * sin(K_y * j) * cos(K_x * i);
					double uyy = Vmax * K_x / sqrt(K_x * K_x + K_y * K_y) * exp(-visc * K * K * t) * sin(K_x * i) * cos(K_y * j);
					double rhoxxyy = 1 - Ma * Ma / 2.0 / K / K * (K_y * K_y * cos(2.0 * K_x * i) + K_x * K_x * cos(2.0 * K_y * j));
					glVertex2f(Lx + i * Lx / l, j * Ly / m);
					glVertex2f(Lx + i * Lx / l + uxx * scale_velo, j * Ly / m + uyy * scale_velo);
				}
			}
		}
		glEnd();
		glPopMatrix();
	}

	virtual GLint marching_cube(void) {
		throw "Marching cube algorithm: Did not implemented!";
	}

protected:
	lattice_node Nodes[l + 1][m + 1];
	GLdouble velocity_before[l + 1][m + 1][D];
	std::vector<std::vector<GLint> > Bounce_index;
};