#pragma once
#include <GL/glut.h>
#include <string>
using namespace std;
class lattice_node {
public:
	GLdouble f[Q];
	GLdouble fnew[Q];
	GLdouble feq[Q];
	GLdouble u[D];
	GLdouble upu[D];
	GLdouble rho;
	GLdouble p;
	GLdouble F[D];
	GLdouble S[D];
	std::string status;

	lattice_node()
	{
		for (int rr = 0; rr < D; rr++)
		{
			u[rr] = 0.0;
			upu[rr] = 0.0;
			F[rr] = 0.0;
			S[rr] = 0.0;
		}
		for (int k = 0; k < Q; k++)
		{
			f[k] = 0.0;
			fnew[k] = 0.0;
			feq[k] = 0.0;
		}
		rho = 0.0;
		p = 0.0;
		status = "L";
	}
};