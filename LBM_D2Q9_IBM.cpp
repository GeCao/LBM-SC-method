#include "LBM_D2Q9_IBM.h"
void LBM_D2Q9_IBM::dump_file(string filename) {
	int i = 0, j = 0, k = 0, rr = 0;
	FILE* fp;

	//Compute Vorticity.
	int ileft = 0, iright = 0, jup = 0, jdown = 0;
#pragma omp parallel for private(i, j, ileft, iright, jup, jdown)
	for (i = 0; i <= l; i++) {
		ileft = (i > 0) ? (i - 1) : (l);
		iright = (i < l) ? (i + 1) : (0);
		for (j = 0; j <= m; j++) {
			jdown = (j > 0) ? (j - 1) : (m);
			jup = (j < m) ? (j + 1) : (0);
			Vorticity[i][j] = 0.5 * (Nodes[i][jup].u[0] - Nodes[i][jdown].u[0]) - \
				0.5 * (Nodes[iright][j].u[1] - Nodes[ileft][j].u[1]);
			if (Nodes[i][j].status[0] == 'B') { Vorticity[i][j] = 0.0; }
		}
	}

	double aerodynamics_force[D] = { 0.0 };
	double torque = 0.0;
#pragma omp parallel for private(i, j, k, rr)
	for (int idx_point = 0; idx_point < idx_near_cylinder_boundary.size(); idx_point++) {
		i = idx_near_cylinder_boundary[idx_point][0];
		j = idx_near_cylinder_boundary[idx_point][1];
		Vorticity[i][j] = 0.0;
		double force_xb[D] = { 0.0 };
		for (k = 0; k < Q; k++) {
			if (is_in_array(k, idx_of_missing_props[idx_point])) { //K: From Fluid to Obstacle: +
				for (rr = 0; rr < D; rr++) {
					force_xb[rr] += e[k][rr] * Nodes[i][j].f[k];
				}
			}
			else { //K: From Obstacle to Fluid: -
				for (rr = 0; rr < D; rr++) {
					force_xb[rr] -= e[k][rr] * Nodes[i][j].f[k];
				}
			}
		}
		for (rr = 0; rr < D; rr++) {
			aerodynamics_force[rr] += force_xb[rr];
		}
		torque += (i - cylinder_center[0]) * force_xb[1] - (j - cylinder_center[1]) * force_xb[0];
	}

	if (!(fp = fopen(filename.c_str(), "w")))
	{
		printf("写入文件失败!\n");
		exit(1);
	}
	fprintf(fp, "TITLE =\"lattice boltzmann\"\n");
	fprintf(fp, "VARIABLES = \"X\" \"Y\" \"RHO\" \"U\" \"V\" \"VORT\" \n");
	fprintf(fp, "ZONE I=%d, J=%d\n", l + 1, m - 1);
	fprintf(fp, "F=POINT\n");
	for (j = 1; j < m; j++) {
		for (i = 0; i <= l; i++) {
			fprintf(fp, "%d %d %lf %lf %lf %lf\n", i, j, Nodes[i][j].rho, Nodes[i][j].u[0], Nodes[i][j].u[1], Vorticity[i][j]);
		}
	}
	fclose(fp);

	/* Log Aerodynamics Forces*/
	FILE* fp_log_aerodynamics_forces;
	if (!(fp_log_aerodynamics_forces = fopen("log_aerodynamics_forces.txt", "a+"))) {
		printf("写入文件失败!\n");
		exit(-1);
	}
	fseek(fp_log_aerodynamics_forces, 0, SEEK_END);
	fprintf(fp_log_aerodynamics_forces, "%d %lf %lf %lf\n", get_time(), aerodynamics_force[0], aerodynamics_force[1], torque);
	fclose(fp_log_aerodynamics_forces);
}

void LBM_D2Q9_IBM::drawing_liquid(void) {
	GLint i = 0, j = 0, ileft = 0, iright = 0, jup = 0, jdown = 0, k = 0;

	glPushMatrix();
	glColor3f(1.0, 0.0, 0.0);
	glPointSize(5.0f);
	glBegin(GL_POINTS);
	double Vort_max = 0.0;
	for (i = 0; i <= l; i++) {
		ileft = (i > 0) ? (i - 1) : (l);
		iright = (i < l) ? (i + 1) : (0);
		for (j = 0; j <= m; j++) {
			jdown = (j > 0) ? (j - 1) : (m);
			jup = (j < m) ? (j + 1) : (0);
			Vorticity[i][j] = 0.5 * (Nodes[i][jup].u[0] - Nodes[i][jdown].u[0]) - \
				0.5 * (Nodes[iright][j].u[1] - Nodes[ileft][j].u[1]);
			if (abs(Vorticity[i][j]) > Vort_max) { Vort_max = abs(Vorticity[i][j]); }
		}
	}
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			if (Nodes[i][j].status == "L") {
				if (Vorticity[i][j] > 0.0) {
					glColor3f(std::min(1.0, abs(Vorticity[i][j]) * 6.0 / Vort_max), 0.0, 0.0);
				}
				else {
					glColor3f(0.0, 0.0, abs(Vorticity[i][j]) / Vort_max);
				}
				glVertex2f(i * Lx / l, j * Ly / m);
			}
		}
	}
	glEnd();
	glPopMatrix();
	
	/*
	GLint i = 0, j = 0;
	GLdouble scale_velo = 10;

	glPushMatrix();
	glColor3f(1.0, 0.0, 0.0);
	glPointSize(5.2);
	glBegin(GL_POINTS);
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			if (Nodes[i][j].status == "L") {
				glColor3f(Nodes[i][j].rho / 1.2, 0.0, 0.0);
				glVertex2f(i * Lx / l, j * Ly / m);
			}
		}
	}
	for (i = 0; i < idx_near_cylinder_boundary.size(); i++) {
		glColor3f(0.0, 1.0, 0.0);
		glVertex2f(idx_near_cylinder_boundary[i][0] * Lx / l, \
			idx_near_cylinder_boundary[i][1] * Ly / m);
	}
	glEnd();
	glPopMatrix();
	*/
}

void LBM_D2Q9_IBM::post_streaming() {
	/*
	* After streaming, before the bounce-back!
	*/
	int i = 0, j = 0, k = 0, opp_k = 0, idx_point = 0, idx_pop = 0;

#pragma omp parallel for private(i, j, k, opp_k, idx_point, idx_pop)
	for (idx_point = 0; idx_point < idx_near_cylinder_boundary.size(); idx_point++) {
		i = idx_near_cylinder_boundary[idx_point][0];
		j = idx_near_cylinder_boundary[idx_point][1];
		double u_tgt[2] = { 0.0, 0.0 };
		double rho_tgt = 0.0;

		//We did not make any change because we suppose that the reflection from cylinder has been corrrected!
		for (k = 0; k < Q; k++) {
			rho_tgt += Nodes[i][j].f[k];
		}

		for (idx_pop = 0; idx_pop < idx_of_missing_props[idx_point].size(); idx_pop++) {
			k = idx_of_missing_props[idx_point][idx_pop]; //The k is not the missing k! (From fluid to obstacle)
			opp_k = opp_index[k]; //The opp_k is actually the missing k! (From obstacle to fluid)

			double t_ij = tij_near_cylinder[idx_point][idx_pop];
			double xw = e[k][0] * t_ij + i, yw = e[k][1] * t_ij + j;

			//Value the target velocity:
			//The "f" is been interpolated with the help of setting velocity of obstacle as zero.
			u_tgt[0] += (t_ij * Nodes[i + e[opp_k][0]][j + e[opp_k][1]].f[opp_k] + 0.0) / (1.0 + t_ij) * e[opp_k][0];
			u_tgt[1] += (t_ij * Nodes[i + e[opp_k][0]][j + e[opp_k][1]].f[opp_k] + 0.0) / (1.0 + t_ij) * e[opp_k][1];

			//Prepare for next step: Reset the boundary value at Cylinder wall
			//Please aware while in this time, the Bounce-back condition has not been applied yet!
			//So only set the direction of fluid->obstacle in obstacle region.
			Nodes[i + e[k][0]][j + e[k][1]].f[k] = Nodes[i][j].f[k] * (1 - t_ij) + \
				Nodes[i + e[k][0]][j + e[k][1]].f[k];
		}
		
		u_tgt[0] /= idx_of_missing_props[idx_point].size() * rho_tgt;
		u_tgt[1] /= idx_of_missing_props[idx_point].size() * rho_tgt;

		double u_grad[2][2] = { {0.0, 0.0}, {0.0, 0.0} };
		if (is_in_array(1, idx_of_missing_props[idx_point])) { //k = 1: From fluid to obstacle; k = 3: From fluid to fluid
			u_grad[0][0] = Nodes[i][j].u[0] - Nodes[i - 1][j].u[0]; //Partial(u)/Partial(x)
			u_grad[1][0] = Nodes[i][j].u[1] - Nodes[i - 1][j].u[1]; //Partial(v)/Partial(x)
		}
		else if (is_in_array(3, idx_of_missing_props[idx_point])) { //k = 3: From fluid to obstacle; k = 1: From fluid to fluid
			u_grad[0][0] = Nodes[i + 1][j].u[0] - Nodes[i][j].u[0]; //Partial(u)/Partial(x)
			u_grad[1][0] = Nodes[i + 1][j].u[1] - Nodes[i][j].u[1]; //Partial(v)/Partial(x)
		}
		else { 
			u_grad[0][0] = 0.5 * (Nodes[i + 1][j].u[0] - Nodes[i - 1][j].u[0]); //Partial(u)/Partial(x)
			u_grad[1][0] = 0.5 * (Nodes[i + 1][j].u[1] - Nodes[i - 1][j].u[1]); //Partial(v)/Partial(x)
		}

		if (is_in_array(2, idx_of_missing_props[idx_point])) { //k = 2: From fluid to obstacle; k = 4: From fluid to fluid
			u_grad[0][1] = Nodes[i][j].u[0] - Nodes[i][j - 1].u[0]; //Partial(u)/Partial(y)
			u_grad[1][1] = Nodes[i][j].u[1] - Nodes[i][j - 1].u[1]; //Partial(v)/Partial(y)
		}
		else if (is_in_array(4, idx_of_missing_props[idx_point])) { //k = 2: From fluid to obstacle; k = 4: From fluid to fluid
			u_grad[0][1] = Nodes[i][j + 1].u[0] - Nodes[i][j].u[0]; //Partial(u)/Partial(y)
			u_grad[1][1] = Nodes[i][j + 1].u[1] - Nodes[i][j].u[1]; //Partial(v)/Partial(y)
		}
		else {
			u_grad[0][1] = 0.5 * (Nodes[i][j + 1].u[0] - Nodes[i][j - 1].u[0]); //Partial(u)/Partial(y)
			u_grad[1][1] = 0.5 * (Nodes[i][j + 1].u[1] - Nodes[i][j - 1].u[1]); //Partial(v)/Partial(y)
		}

		double P_ab[D][D] = { {0.0, 0.0}, {0.0, 0.0} };
		P_ab[0][0] = rho_tgt * cs2 + rho_tgt * u_tgt[0] * u_tgt[0] - \
			rho_tgt * cs2 * u_grad[0][0] / beta;
		P_ab[1][0] = 0 + rho_tgt * u_tgt[1] * u_tgt[0] - \
			rho_tgt * cs2 * (u_grad[0][1] + u_grad[1][0]) / 2.0 / beta;
		P_ab[0][1] = P_ab[1][0];
		P_ab[1][1] = rho_tgt * cs2 + rho_tgt * u_tgt[1] * u_tgt[1] - \
			rho_tgt * cs2 * u_grad[1][1] / beta;

		for (idx_pop = 0; idx_pop < idx_of_missing_props[idx_point].size(); idx_pop++) {
			k = idx_of_missing_props[idx_point][idx_pop]; //The k is not the missing k! (From fluid to obstacle)
			opp_k = opp_index[k]; //The opp_k is actually the missing k! (From obstacle to fluid)

			Nodes[i][j].f[opp_k] = rho_tgt;
			for (int ii = 0; ii < D; ii++) {
				for (int jj = 0; jj < D; jj++) {
					Nodes[i][j].f[opp_k] += (P_ab[ii][jj] - rho_tgt * cs2 * ((ii == jj) ? 1 : 0)) * \
						(e[opp_k][ii] * e[opp_k][jj] * c2 - cs2 * ((ii == jj) ? 1 : 0)) / (2.0 * cs2 * cs2);
				}
				Nodes[i][j].f[opp_k] += rho_tgt * u_tgt[ii] * e[opp_k][ii] * c / cs2;
			}
			Nodes[i][j].f[opp_k] *= w[opp_k];
		}
		
	}
}