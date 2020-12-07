#pragma once
#include <iostream>
#include<GL/glut.h>
#include "utility.h"
class LBM
{
public:
	virtual void Bounceback_Boundary() {
		throw "Set Bounce back boundary in LBM: Not implemented!";
	}

	virtual void calculate_feq(GLdouble rh, GLdouble u[], GLdouble feq_temp[]) {
		throw "Caluate f_eq: Not implemented!";
	}

	virtual double Cal_error() {
		throw "Calculate the iteration error in LBM: Not implemented!";
	}

	virtual void Record_velocity(void) {
		throw "Record the velocity in last iteration step: Not implemented!";
	}

	virtual void calcu_upr() {
		throw "Calculate new velocity for force and source terms in S-C model: Not implemented!";
	}

	virtual void calcu_Fxy() {
		throw "Calculate the force and source terms in S-C model: Not implemented!";
	}

	virtual void rebounce_f() {
		throw "Rebounce f in Bounceback boundary area: Not implemented!";
	}

	virtual void Streaming() {
		throw "Streaming process: Not implemented!";
	}

	virtual void macro_process() {
		throw "macro process for calculate density and velocity: Not implemented!";
	}

	virtual void Collision(void) {
		throw "Collision step: Not implemented!";
	}

	virtual void drawing_liquid(void) {
		throw "Drawing liquid: Not implemented!";
	}

	virtual void dump_file(std::string filename) {
		throw "Dump file: Not implemented!";
	}

	virtual void read_file(std::string filename) {
		throw "Read saved files: Not implemented!";
	}

	virtual void compute() {
		//Record the velocity
		Record_velocity();

		//Streaming process
		Streaming();

		//macro process
		macro_process();

		rebounce_f();
		//Collision process
		Collision();
	}

	virtual GLint marching_cube(void) {
		throw "Marching cube algorithm: Did not implemented!";
	}

	virtual int get_tmax() { return t_max; }
	virtual bool ifdumpfile() { return if_dump; }

	LBM() : dt(1.0), c(1.0), RR(10.0), UU(-0.0), rho_w(0.08), t_max(4000), if_dump(true), radius(0.15) {
		c2 = c * c;
		cs2 = c2 / 3.0;
		p0 = 0.0; //初始压强 | The initial pressure of particles.
		visc = cs2 * (tau - 0.5); //初始粘度 | The initial viscosity of particles.

		if (exercise_id == 2) {
			Vmax = 0.2; //ex2
			Re = Vmax * l / visc; //ex2
			tau = 1.0; //弛豫时间
			beta = 1.0 / 2.0 / tau;
		}
		else if (exercise_id == 3) {
			beta = 0.85;
			tau = 1.0 / 2.0 / beta; //弛豫时间
			Re = 3000.0; //ex3
			Vmax = Re * visc / l; //ex3
		}
		else {
			Vmax = 0.2; //ex2
			Re = Vmax * l / visc; //ex2
			tau = 1.0; //弛豫时间
			beta = 1.0 / 2.0 / tau;
		}
		Ma = Vmax / cs2;

		t_max = 4000; //最大的迭代计算步数 | The max steps this programm run.
	}

protected:
	double dt; //时间间隔 | The time step
	double c; //格子声速
	double c2;
	double cs2;
	double p0; //初始压强 | The initial pressure of particles.
	double visc; //初始粘度 | The initial viscosity of particles.

	//parameters
	double RR; //液滴半径 | The radius of droplet.
	double UU; //液滴初始速度 | The initial velocity of droplet.
	double rho_w; //壁面密度 | The density of wall.

	double Ma;
	double beta;
	double tau;
	double Re;
	double Vmax;

	//Visualization
	int t_max; //最大的迭代计算步数 | The max steps this programm run.
	bool if_dump; //是否进行文件的输出 | Whether do we need to dump files.
	double radius; //在以小球的形式绘制流场时，小球的半径 | while fluid particle was seen as a little sphere, hereby defined its radius.
};
