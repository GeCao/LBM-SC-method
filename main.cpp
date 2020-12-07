//A program for lattice Boltzmann Method
#pragma comment(lib,"glew32.lib")
#include <GL/glew.h>
#include <GL/glut.h>
#include<fstream>
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include <vector>
#include<omp.h>
#include<math.h>
#include<string.h>
#include "utility.h"
#include "LBM_D3Q19_SC.h"
#include "LBM_D2Q9_KBC.h"
using namespace std;

GLint control = 1; //该值为1时，执行计算 | While this value been set as 1, this programm can run.
GLint case_index = 0; //图层层次
GLint Nwri = 10; //每隔Nwri步，进行一次文件的输出 | Every "Nwri" steps, dump a file.
GLint Ndraw = 1; //每隔Ndraw步，进行一次图形的绘制 | Every "Ndraw" steps, draw a picture.

GLuint vboId;//vertex buffer object句柄    
GLuint vaoId;//vertext array object句柄    
GLuint programId;//shader program 句柄  
GLint size_vertices = 0;

//LBM_D3Q19 has not implemented yet, but LBM_D3Q19_SC did
LBM_D3Q19_SC lbm0;

GLint time = 0; //时间步 | time step.

void init(void) {
	if (D == 3)
	{
		glClearColor(0.0, 0.0, 0.0, 0.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);

		gluPerspective(fovy, aspect, zFar, zNear);
		gluLookAt(eyex, eyey, eyez, centerx, centery, centerz, upx, upy, upz);

		//Create VAO object
		glGenVertexArrays(1, &vaoId);
		glBindVertexArray(vaoId);

		//Create VBO object
		glGenBuffers(1, &vboId);
		glBindBuffer(GL_ARRAY_BUFFER, vboId);
	}
	else if (D == 2) {
		glLoadIdentity();
		glClearColor(1.0, 1.0, 1.0, 0.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glMatrixMode(GL_PROJECTION);
		gluOrtho2D(winx1, winx2, winy1, winy2);
	}
}

void LBM_method()
{
	GLdouble final_error = 1.0;

	omp_set_num_threads(NUM_THREADS);
	if (time == 0) {
		lbm0.drawing_liquid();
		glutSwapBuffers();
		system("pause");
	}

	lbm0.compute();

	if (time % Ndraw == 0)
	{
		if (D == 3) {
			size_vertices = lbm0.marching_cube();
		}
		else if (exercise_id == 2) {
			lbm0.drawing_liquid();
			//lbm0.drawing_analytical_liquid((double)time);
		}
		else if (exercise_id == 3) {
			lbm0.drawing_liquid();
		}
	}

	if (time % Nwri == 0)
	{
		std::cout << "=====================" << endl;
		final_error = lbm0.Cal_error();
		std::cout << "Step: " << time << ";  error: " << final_error << endl;
		if (exercise_id == 2) {
			//std::cout << "Analytical error: " << lbm0.Analytical_error(time) << endl;
		}
		if (lbm0.ifdumpfile()) { 
			char str_step[10];
			itoa(time / Nwri, str_step, 10);
			lbm0.dump_file(string("data") + str_step + ".plt");
		}
	}
}

static void display(void)
{
	init();
	if (control == 1)
	{
		LBM_method();
		time++;
	}
	if (time > lbm0.get_tmax())
	{
		system("pause");
		exit(0);
	}

	glutSwapBuffers();
}

void chooseMode(GLint menuIteemNum)
{
	switch (menuIteemNum)
	{
	case 0:
		case_index = 0; break;
	case 1:
		case_index = 1; break;
	case 2:
		case_index = 2; break;
	default:
		case_index = -1; break;
	}
	glutPostRedisplay();
}

void mouseFunc(GLint button, GLint action, GLint xMouse, GLint yMouse)
{
	glutCreateMenu(chooseMode);
	glutAddMenuEntry("figure : LBM", 0);

	glutAttachMenu(GLUT_RIGHT_BUTTON);
}

static void idle(void)
{
	glutPostRedisplay();
}

int main(int argc, char* argv[])
{
	glutInit(&argc, argv);
	if (exercise_id == 2) {
		glutInitWindowSize(800, 400);
	}
	else {
		glutInitWindowSize(400, 400);
	}
	glutInitWindowPosition(10, 10);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("Fluid dynamics with LBM");

	if(D == 3)
	    glewInit();

	init();

	glutDisplayFunc(display);
	glutMouseFunc(mouseFunc);
	glutIdleFunc(idle);

	/*
	glEnable(GL_CULL_FACE);
	glEnable(GL_LIGHTING);
	const GLfloat light_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
	const GLfloat light_diffuse[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	const GLfloat light_position[] = { Lx + 2.0f, Ly + 5.0f, Lz + 5.0f, 0.0f };

	const GLfloat mat_ambient[] = { 1.0f, 0.0f, 0.31f, 1.0f };
	const GLfloat mat_diffuse[] = { 1.0f, 0.5f, 0.31f, 1.0f };
	const GLfloat mat_specular[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	const GLfloat high_shininess[] = { 32.0f };
	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);
	*/

	glutMainLoop();

	return EXIT_SUCCESS;
}

GLint LBM_D3Q19_SC::marching_cube(void)
{
	GLint i = 0, j = 0, r = 0, k = 0;
	GLint rho_map[l + 1][m + 1][h + 1];
	GLint temp_index = 0, multi = 1;
	GLdouble ave_rho = 0.5 * (rho_h + rho_l);
	GLint index1, index2, index3;
	GLfloat Point1[3], Point2[3], Point3[3];
	GLfloat D_Point12[3], D_Point23[3];
	GLint count_size = 0;

#pragma omp parallel for private(i,j,r)
	for (i = 1; i < l - 1; i++) {
		for (j = 1; j < m - 1; j++) {
			for (r = 1; r < h - 1; r++) {
				if (Nodes[i][j][r].rho > ave_rho) {
					rho_map[i][j][r] = 1;
				}
				else {
					rho_map[i][j][r] = 0;
				}
			}
		}
	}

	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_CULL_FACE);
	glFrontFace(GL_CCW);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glDepthMask(true);
	glColor4f(0.8, 0.8, 0.8, 0.4);
	GLfloat vertices[(l - 3) * (m - 3) * (h - 3) * 9];
	GLint idx_vertices = 0;
	for (i = 1; i < l - 2; i++) {
		for (j = 1; j < m - 2; j++) {
			for (r = 1; r < h - 2; r++) {
				temp_index = 0;
				multi = 1;
				for (k = 0; k < 8; k++) {
					temp_index += rho_map[i + a2fVertexOffset[k][0]][j + a2fVertexOffset[k][1]][r + a2fVertexOffset[k][2]] * multi;
					multi *= 2;
				}
				for (k = 0; k < 16; k += 3) {
					if (temp_index >= 256 || a2iTriangleConnectionTable[temp_index][k] == -1) {
						break;
					}
					index1 = a2iTriangleConnectionTable[temp_index][k];//得到边的索引
					index2 = a2iTriangleConnectionTable[temp_index][k + 1];//得到边的索引
					index3 = a2iTriangleConnectionTable[temp_index][k + 2];//得到边的索引

					//进行三角形的面片表示
					vertices[idx_vertices++] = Point1[0] = i + 0.5 * (a2fVertexOffset[a2iEdgeConnection[index1][0]][0] + a2fVertexOffset[a2iEdgeConnection[index1][1]][0]);
					vertices[idx_vertices++] = Point1[1] = j + 0.5 * (a2fVertexOffset[a2iEdgeConnection[index1][0]][1] + a2fVertexOffset[a2iEdgeConnection[index1][1]][1]);
					vertices[idx_vertices++] = Point1[2] = r + 0.5 * (a2fVertexOffset[a2iEdgeConnection[index1][0]][2] + a2fVertexOffset[a2iEdgeConnection[index1][1]][2]);

					vertices[idx_vertices++] = Point2[0] = i + 0.5 * (a2fVertexOffset[a2iEdgeConnection[index2][0]][0] + a2fVertexOffset[a2iEdgeConnection[index2][1]][0]);
					vertices[idx_vertices++] = Point2[1] = j + 0.5 * (a2fVertexOffset[a2iEdgeConnection[index2][0]][1] + a2fVertexOffset[a2iEdgeConnection[index2][1]][1]);
					vertices[idx_vertices++] = Point2[2] = r + 0.5 * (a2fVertexOffset[a2iEdgeConnection[index2][0]][2] + a2fVertexOffset[a2iEdgeConnection[index2][1]][2]);

					vertices[idx_vertices++] = Point3[0] = i + 0.5 * (a2fVertexOffset[a2iEdgeConnection[index3][0]][0] + a2fVertexOffset[a2iEdgeConnection[index3][1]][0]);
					vertices[idx_vertices++] = Point3[1] = j + 0.5 * (a2fVertexOffset[a2iEdgeConnection[index3][0]][1] + a2fVertexOffset[a2iEdgeConnection[index3][1]][1]);
					vertices[idx_vertices++] = Point3[2] = r + 0.5 * (a2fVertexOffset[a2iEdgeConnection[index3][0]][2] + a2fVertexOffset[a2iEdgeConnection[index3][1]][2]);


					//进行法线的计算和描述
					D_Point12[0] = Point1[0] - Point2[0]; D_Point12[1] = Point1[1] - Point2[1]; D_Point12[2] = Point1[2] - Point2[2];
					D_Point23[0] = Point2[0] - Point3[0]; D_Point23[1] = Point2[1] - Point3[1]; D_Point23[2] = Point2[2] - Point3[2];
					glNormal3f(D_Point12[1] * D_Point23[2] - D_Point12[2] * D_Point23[1],
						-(D_Point12[0] * D_Point23[2] - D_Point12[2] * D_Point23[0]),
						D_Point12[0] * D_Point23[1] - D_Point12[1] * D_Point23[0]);

					/*
					//进行水的贴图和三角形的面片表示
					glBegin(GL_TRIANGLES);
					glTexCoord2f(0.0f, 0.0f);
					glVertex3fv(Point1);
					glTexCoord2f(1.0f, 0.0f);
					glVertex3fv(Point2);
					glTexCoord2f(1.0f, 1.0f);
					glVertex3fv(Point3);
					glEnd();
					*/
					count_size += 3;
					//glTexCoord2f(0.0f, 1.0f); glVertex3d(0.5*(Point3[0] + Point1[0]), 0.5*(Point3[1] + Point1[1]), 0.5*(Point3[2] + Point1[2]));
				}
			}
		}
	}
	//传入VBO数据
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STREAM_DRAW);
	//解除VBO绑定
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glClear(GL_COLOR_BUFFER_BIT);

	//绑定VBO
	glBindBuffer(GL_ARRAY_BUFFER, vboId);
	glEnableVertexAttribArray(0);

	//解释顶点数据方式
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	//绘制模型
	glDrawArrays(GL_TRIANGLES, 0, size_vertices);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glDisableVertexAttribArray(0);
	return count_size;
}