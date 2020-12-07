#pragma once
#include <math.h>
#define NUM_THREADS 4
//定义pi的值 | Define the value of pi
#define pi 3.1415926
//定义D*Q*模型的参数
#define Q 19
#define D 3
#define exercise_id 3
//网格边界面：x=0,x=l,y=0,y=m.z=0,z=h,内部点为(l-1)*(m-1)*(h-1)个.
static char up_boundary = 'P';//上边界为Periodic/Bounceback/Von Neumann(Flux)/Dirichlet boundaries? choose P,B,V,D for your prefer.
static char down_boundary = 'P';//下边界为Periodic/Bounceback/Von Neumann(Flux)/Dirichlet boundaries? choose P,B,V,D for your prefer.
static char left_boundary = 'P';//左边界为Periodic/Bounceback/Von Neumann(Flux)/Dirichlet boundaries? choose P,B,V,D for your prefer.
static char right_boundary = 'P';//右边界为Periodic/Bounceback/Von Neumann(Flux)/Dirichlet boundaries? choose P,B,V,D for your prefer.
static char front_boundary = 'P';//前边界为Periodic/Bounceback/Von Neumann(Flux)/Dirichlet boundaries? choose P,B,V,D for your prefer.
static char behind_boundary = 'P';//后边界为Periodic/Bounceback/Von Neumann(Flux)/Dirichlet boundaries? choose P,B,V,D for your prefer.
/*前处理模块结束*/
#if D == 3
#define l 60
#define m 60
#define h 60
static double w[Q] = { 1.0 / 3.0, //Original point,a=0
1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, //direction: nx,ny,px,py,a=1,2,3,4
1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, //direction: nxny,pxny,pxpy,nxpy,a=5,6,7,8
1.0 / 18.0, //direction: nz,a=9
1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, //direction:  nxnz,nynz,,pxnz,pynz,a=10,11,12,13
1.0 / 18.0, //direction: pz,a=14
1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 //direction:  nxpz,nypz,pxpz,pypz,a=15,16,17,18
};
static double e[Q][D] = { { 0.0, 0.0, 0.0 }, //Original point,a=0
{ 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { -1.0, 0.0, 0.0 }, { 0.0, -1.0, 0.0 }, //direction: nx,ny,px,py,a=1,2,3,4
{ 1.0, 1.0, 0.0 }, { -1.0, 1.0, 0.0 }, { -1.0, -1.0, 0.0 }, { 1.0, -1.0, 0.0 }, //direction: nxny,pxny,pxpy,nxpy,a=5,6,7,8
{ 0.0, 0.0, 1.0 }, //direction: nz,a=9
{ 1.0, 0.0, 1.0 }, { 0.0, 1.0, 1.0 }, { -1.0, 0.0, 1.0 }, { 0.0, -1.0, 1.0 }, //direction: nxnz,nynz,,pxnz,pynz,a=10,11,12,13
{ 0.0, 0.0, -1.0 }, //direction: pz,a=14
{ 1.0, 0.0, -1.0 }, { 0.0, 1.0, -1.0 }, { -1.0, 0.0, -1.0 }, { 0.0, -1.0, -1.0 } //direction: nxpz,nypz,pxpz,pypz,a=15,16,17,18
};
static int opp_index[Q] = { 0, 3, 4, 1, 2, 7, 8, 5, 6, 14, 17, 18, 15, 16, 9, 12, 13, 10, 11 };
#elif D == 2
#define l 128
#define m 128
#define h 4
static double w[Q] = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };
static double e[Q][D] = { { 0.0, 0.0 },
{ 1.0, 0.0 }, { 0.0, 1.0 }, { -1.0, 0.0 }, { 0.0, -1.0 },
{ 1.0, 1.0 }, { -1.0, 1.0 }, { -1.0, -1.0 }, { 1.0, -1.0 } };
static int opp_index[Q] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };
#endif
static double Lx = l, Ly = m, Lz = h; //流场的长与宽

//TT0W is the value of T/T0;   RHW and RLW are the coexisting densities in the sepcified T/T0.
static double TT0W[12] = { 0.975, 0.95, 0.925, 0.9, 0.875, 0.85, 0.825, 0.8, 0.775, 0.75, 0.7, 0.65 };
static double RHW[12] = { 0.16, 0.21, 0.23, 0.247, 0.265, 0.279, 0.29, 0.314, 0.30, 0.33, 0.36, 0.38 };
static double RLW[12] = { 0.08, 0.067, 0.05, 0.0405, 0.038, 0.032, 0.025, 0.0245, 0.02, 0.015, 0.009, 0.006 };
//Please specify which temperature and corresponding \rho_h, \rho_l in above 'data' are chosen.
//Initial T/T0, rho_h, and rho_l for the C-S EOS are listed in above 'data' section
static double TT0 = TT0W[4]; //标准温度，用于气体状态方程 | Standard temprature, been used for EOS.
static double rho_h = RHW[4]; //液滴密度 | The density of droplet.
static double rho_l = RLW[4]; //空气密度 | The density of air.

#if exercise_id == 2
static double winx1 = 0.0, winx2 = 2.0 * Lx, winy1 = 0.0, winy2 = Ly;//Window size
#else
static double winx1 = 0.0, winx2 = Lx, winy1 = 0.0, winy2 = Ly;//Window size
#endif

static double KBC_sigma = 0.05, KBC_kappa = 80.0; //ex3

static double lambda_alpha = 1.0;
static double K_x = 2 * pi / lambda_alpha / l, K_y = 2 * pi / lambda_alpha / m;
static double K = sqrt(K_x * K_x + K_y * K_y);

static GLdouble fovy = 35, aspect = 1, zFar = 1.0 + Lz, zNear = -1.0;//投影
static GLdouble eyex = -0.4 * Lx, eyey = -0.7 * Ly, eyez = 1.1 * Lz, centerx = 0.5 * Lx, centery = 0.5 * Ly, centerz = 0.5 * Lz, upx = 0, upy = 0.0, upz = 1.0;//三维观察