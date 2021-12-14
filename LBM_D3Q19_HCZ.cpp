#include "LBM_D3Q19_HCZ.h"
#ifdef use_HCZ == 1
void LBM_D3Q19_HCZ::Bounceback_Boundary()
{
	GLint i = 0, j = 0, r = 0;

	if (up_boundary == 'B')
	{
#pragma omp parallel for private(i,j)
		for (i = 0; i <= l; i++)
		{
			for (j = 1; j < m; j++)
			{
				Nodes[i][j][h - 1].status = "Bup";
			}
		}
	}
	if (down_boundary == 'B')
	{
#pragma omp parallel for private(i,j)
		for (i = 0; i <= l; i++)
		{
			for (j = 0; j < m; j++)
			{
				Nodes[i][j][1].status = "Bdown";
			}
		}
	}

	if (front_boundary == 'B')
	{
#pragma omp parallel for private(i,r)
		for (i = 0; i <= l; i++)
		{
			for (r = 1; r < h; r++)
			{
				Nodes[i][1][r].status = "Bfront";
			}
		}
	}
	if (behind_boundary == 'B')
	{
#pragma omp parallel for private(i,r)
		for (i = 0; i <= l; i++)
		{
			for (r = 1; r < h; r++)
			{
				Nodes[i][m - 1][r].status = "Bbehind";
			}
		}
	}

	if (left_boundary == 'B')
	{
#pragma omp parallel for private(j,r)
		for (j = 1; j < m; j++)
		{
			for (r = 1; r < h; r++)
			{
				Nodes[0][j][r].status = "Bleft";
			}
		}
	}
	if (right_boundary == 'B')
	{
#pragma omp parallel for private(j,r)
		for (j = 1; j < m; j++)
		{
			for (r = 1; r < h; r++)
			{
				Nodes[l][j][r].status = "Bright";
			}
		}
	}

	return;
}

void LBM_D3Q19_HCZ::Dirichlet_Boundary()
{
	GLint i = 0, j = 0, r = 0;

	if (up_boundary == 'D')
	{
#pragma omp parallel for private(i,j)
		for (i = 0; i <= l; i++)
		{
			for (j = 1; j < m; j++)
			{
				Nodes[i][j][h - 1].status = "Dup";
			}
		}
	}
	if (down_boundary == 'D')
	{
#pragma omp parallel for private(i,j)
		for (i = 0; i <= l; i++)
		{
			for (j = 1; j < m; j++)
			{
				Nodes[i][j][1].status = "Ddown";
			}
		}
	}

	if (front_boundary == 'D')
	{
#pragma omp parallel for private(i,r)
		for (i = 0; i <= l; i++)
		{
			for (r = 1; r < h; r++)
			{
				Nodes[i][1][r].status = "Dfront";
			}
		}
	}
	if (behind_boundary == 'D')
	{
#pragma omp parallel for private(i,r)
		for (i = 0; i <= l; i++)
		{
			for (r = 1; r < h; r++)
			{
				Nodes[i][m - 1][r].status = "Dbehind";
			}
		}
	}

	if (left_boundary == 'D')
	{
#pragma omp parallel for private(j,r)
		for (j = 1; j < m; j++)
		{
			for (r = 1; r < h; r++)
			{
				Nodes[0][j][r].status = "Dleft";
			}
		}
	}
	if (right_boundary == 'D')
	{
#pragma omp parallel for private(j,r)
		for (j = 1; j < m; j++)
		{
			for (r = 1; r < h; r++)
			{
				Nodes[l][j][r].status = "Dright";
			}
		}
	}

	return;
}

void LBM_D3Q19_HCZ::calculate_feq_geq(GLdouble rh, GLdouble rho, GLdouble p, GLdouble u[], GLdouble feq_temp[], GLdouble geq_temp[])
{
	GLdouble eu = 0.0, uv = 0.0;
	GLint k = 0, rr = 0;

	for (rr = 0; rr < D; rr++)
	{
		uv += (u[rr] * u[rr]);
	}
#pragma omp parallel for private(k,rr,eu)
	for (k = 0; k < Q; k++)
	{
		eu = 0.0;
		for (rr = 0; rr < D; rr++)
		{
			eu += (e[k][rr] * u[rr]);
		}
		eu *= c;
		feq_temp[k] = eu / cs2 + 0.5 * eu * eu / cs2 / cs2 - 0.5 * uv / cs2;
		geq_temp[k] = w[k] * (p + cs2 * rho * (feq_temp[k]));
		feq_temp[k] = w[k] * rh * (1.0 + feq_temp[k]);
	}
	return;
}

void LBM_D3Q19_HCZ::rebounce_g()
{
	GLint i = 0, j = 0, r = 0;
	GLdouble swap_temp = 0.0;
#pragma omp parallel for private(i,j,r,swap_temp)
	for (i = 0; i <= l; i++) {
		for (j = 0; j <= m; j++) {
			for (r = 0; r < h; r++) {
				if (Nodes[i][j][r].status[0] == 'B')
				{
					swap_temp = Nodes[i][j][r].g[1]; Nodes[i][j][r].g[1] = Nodes[i][j][r].g[3]; Nodes[i][j][r].g[3] = swap_temp;
					swap_temp = Nodes[i][j][r].g[2]; Nodes[i][j][r].g[2] = Nodes[i][j][r].g[4]; Nodes[i][j][r].g[4] = swap_temp;
					swap_temp = Nodes[i][j][r].g[5]; Nodes[i][j][r].g[5] = Nodes[i][j][r].g[7]; Nodes[i][j][r].g[7] = swap_temp;
					swap_temp = Nodes[i][j][r].g[6]; Nodes[i][j][r].g[6] = Nodes[i][j][r].g[8]; Nodes[i][j][r].g[8] = swap_temp;

					swap_temp = Nodes[i][j][r].g[9]; Nodes[i][j][r].g[9] = Nodes[i][j][r].g[14]; Nodes[i][j][r].g[14] = swap_temp;

					swap_temp = Nodes[i][j][r].g[10]; Nodes[i][j][r].g[10] = Nodes[i][j][r].g[17]; Nodes[i][j][r].g[17] = swap_temp;
					swap_temp = Nodes[i][j][r].g[11]; Nodes[i][j][r].g[11] = Nodes[i][j][r].g[18]; Nodes[i][j][r].g[18] = swap_temp;
					swap_temp = Nodes[i][j][r].g[12]; Nodes[i][j][r].g[12] = Nodes[i][j][r].g[15]; Nodes[i][j][r].g[15] = swap_temp;
					swap_temp = Nodes[i][j][r].g[13]; Nodes[i][j][r].g[13] = Nodes[i][j][r].g[16]; Nodes[i][j][r].g[16] = swap_temp;
				}
				//end if.
			}
		}
	}
	//The end of 3-D circulation.
}

void LBM_D3Q19_HCZ::Streaming()
{
	int i = 0, j = 0, r = 0, ileft = 0, iright = 0, jfront = 0, jbehind = 0, rup = 0, rdown = 0, k = 0;

	//Streaming process::L
#pragma omp parallel for private(i,j,r,ileft,iright,jbehind,jfront,rdown,rup)
	for (i = 0; i <= l; i++)
	{
		ileft = (i > 0) ? (i - 1) : (l);
		iright = (i < l) ? (i + 1) : (0);
		for (j = 1; j < m; j++)
		{
			jfront = (j > 1) ? (j - 1) : (m - 1);
			jbehind = (j < m - 1) ? (j + 1) : (1);
			for (r = 1; r < h; r++)
			{
				rdown = (r > 1) ? (r - 1) : (h - 1);
				rup = (r < h - 1) ? (r + 1) : (1);

				Nodes[iright][j][r].fnew[1] = Nodes[i][j][r].f[1];
				Nodes[i][jbehind][r].fnew[2] = Nodes[i][j][r].f[2];
				Nodes[ileft][j][r].fnew[3] = Nodes[i][j][r].f[3];
				Nodes[i][jfront][r].fnew[4] = Nodes[i][j][r].f[4];

				Nodes[iright][jbehind][r].fnew[5] = Nodes[i][j][r].f[5];
				Nodes[ileft][jbehind][r].fnew[6] = Nodes[i][j][r].f[6];
				Nodes[ileft][jfront][r].fnew[7] = Nodes[i][j][r].f[7];
				Nodes[iright][jfront][r].fnew[8] = Nodes[i][j][r].f[8];

				Nodes[i][j][rup].fnew[9] = Nodes[i][j][r].f[9];

				Nodes[iright][j][rup].fnew[10] = Nodes[i][j][r].f[10];
				Nodes[i][jbehind][rup].fnew[11] = Nodes[i][j][r].f[11];
				Nodes[ileft][j][rup].fnew[12] = Nodes[i][j][r].f[12];
				Nodes[i][jfront][rup].fnew[13] = Nodes[i][j][r].f[13];

				Nodes[i][j][rdown].fnew[14] = Nodes[i][j][r].f[14];

				Nodes[iright][j][rdown].fnew[15] = Nodes[i][j][r].f[15];
				Nodes[i][jbehind][rdown].fnew[16] = Nodes[i][j][r].f[16];
				Nodes[ileft][j][rdown].fnew[17] = Nodes[i][j][r].f[17];
				Nodes[i][jfront][rdown].fnew[18] = Nodes[i][j][r].f[18];
			}
		}
	}
	//the end of 3-D circulation.
#pragma omp parallel for private(i,j,r,k)
	for (i = 0; i <= l; i++)
	{
		for (j = 1; j < m; j++)
		{
			for (r = 1; r < h; r++)
			{
				for (k = 1; k < Q; k++)
				{
					Nodes[i][j][r].f[k] = Nodes[i][j][r].fnew[k];
				}
			}
		}
	}
	//the end of 3-D circulation.

	//Streaming process::L
#pragma omp parallel for private(i,j,r,ileft,iright,jbehind,jfront,rdown,rup)
	for (i = 0; i <= l; i++)
	{
		ileft = (i > 0) ? (i - 1) : (l);
		iright = (i < l) ? (i + 1) : (0);
		for (j = 1; j < m; j++)
		{
			jfront = (j > 1) ? (j - 1) : (m - 1);
			jbehind = (j < m - 1) ? (j + 1) : (1);
			for (r = 1; r < h; r++)
			{
				rdown = (r > 1) ? (r - 1) : (h - 1);
				rup = (r < h - 1) ? (r + 1) : (1);

				Nodes[iright][j][r].fnew[1] = Nodes[i][j][r].g[1];
				Nodes[i][jbehind][r].fnew[2] = Nodes[i][j][r].g[2];
				Nodes[ileft][j][r].fnew[3] = Nodes[i][j][r].g[3];
				Nodes[i][jfront][r].fnew[4] = Nodes[i][j][r].g[4];

				Nodes[iright][jbehind][r].fnew[5] = Nodes[i][j][r].g[5];
				Nodes[ileft][jbehind][r].fnew[6] = Nodes[i][j][r].g[6];
				Nodes[ileft][jfront][r].fnew[7] = Nodes[i][j][r].g[7];
				Nodes[iright][jfront][r].fnew[8] = Nodes[i][j][r].g[8];

				Nodes[i][j][rup].fnew[9] = Nodes[i][j][r].g[9];

				Nodes[iright][j][rup].fnew[10] = Nodes[i][j][r].g[10];
				Nodes[i][jbehind][rup].fnew[11] = Nodes[i][j][r].g[11];
				Nodes[ileft][j][rup].fnew[12] = Nodes[i][j][r].g[12];
				Nodes[i][jfront][rup].fnew[13] = Nodes[i][j][r].g[13];

				Nodes[i][j][rdown].fnew[14] = Nodes[i][j][r].g[14];

				Nodes[iright][j][rdown].fnew[15] = Nodes[i][j][r].g[15];
				Nodes[i][jbehind][rdown].fnew[16] = Nodes[i][j][r].g[16];
				Nodes[ileft][j][rdown].fnew[17] = Nodes[i][j][r].g[17];
				Nodes[i][jfront][rdown].fnew[18] = Nodes[i][j][r].g[18];
			}
		}
	}
	//the end of 3-D circulation.
#pragma omp parallel for private(i,j,r,k)
	for (i = 0; i <= l; i++)
	{
		for (j = 1; j < m; j++)
		{
			for (r = 1; r < h; r++)
			{
				for (k = 1; k < Q; k++)
				{
					Nodes[i][j][r].g[k] = Nodes[i][j][r].fnew[k];
				}
			}
		}
	}
	//the end of 3-D circulation.
	return;
}

void LBM_D3Q19_HCZ::macro_process()
{
	int i = 0, j = 0, r = 0, k = 0, rr = 0;
	int hi = 0, hj = 0;
	int ileft = 0, iright = 0, jfront = 0, jbehind = 0, rup = 0, rdown = 0;
	double psi_min0 = psi_min;
	double psi_max0 = psi_max;
	double temp_rh = 0.0;
	int temp_inter1 = 0, temp_inter2 = 0;
	double temp_angle = 0.0, theta = 0.0;

	//macroscopic process
#pragma omp parallel for private(i,j,r,rr,k,temp_rh)
	for (i = 0; i <= l; i++) {
		for (j = 1; j < m; j++) {
			for (r = 1; r < h; r++) {
				//calculation of rh & rho.
				if (Nodes[i][j][r].status == "L")//压力边界下的rh变量是固定的，固体壁面的rh应作为插值的容器已经得到。
				{
					Nodes[i][j][r].rh = 0.0;
					for (k = 0; k < Q; k++)
					{
						Nodes[i][j][r].rh += Nodes[i][j][r].f[k];
					}
				}
				//计算内部流体和固体边界的rho.
				Nodes[i][j][r].rho = rho_l + (Nodes[i][j][r].rh - psi_min0) / (psi_max0 - psi_min0) * (rho_h - rho_l);

				//calculation of prho & fai.
				if (Nodes[i][j][r].status == "L" || Nodes[i][j][r].status[0] == 'D')//压力边界上的场点也算是一种流体
				{
					temp_rh = b * Nodes[i][j][r].rh / 4.0;
					Nodes[i][j][r].prho = Nodes[i][j][r].p - cs2 * Nodes[i][j][r].rho;
					Nodes[i][j][r].fai = Nodes[i][j][r].rh * RT * (4 * temp_rh - 2 * temp_rh * temp_rh) / pow((1 - temp_rh), 3) - a * Nodes[i][j][r].rh * Nodes[i][j][r].rh;
				}
				//end if.
			}
		}
	}
	//The end of 3-D circulation.

	//========================================================================
	//==                                                                    ==
	//==                         接触角的计算                               ==
	//==                                                                    ==
	//========================================================================

	if (down_boundary == 'B' && up_boundary == 'B') {
#pragma omp parallel for private(i,j,r,ileft,iright,jfront,jbehind,rup,rdown,temp_angle,theta)
		for (i = 0; i <= l; i++)//x方向是显然的周期性边界条件
		{
			ileft = (i > 0) ? (i - 1) : (l);
			iright = (i < l) ? (i + 1) : (0);
			for (j = 1; j < m; j++)
			{
				jfront = j - 1;
				jbehind = j + 1;

				//下地面:全部设置为固壁边界条件
				temp_angle = sqrt(pow(Nodes[iright][j][2].rh - Nodes[ileft][j][2].rh, 2) + pow(Nodes[i][jbehind][2].rh - Nodes[i][jfront][2].rh, 2));
				if (temp_angle == 0 && Nodes[i][j][1].rh >= Nodes[i][j][3].rh) {
					theta = 0;
				}
				else if (temp_angle == 0 && Nodes[i][j][1].rh < Nodes[i][j][3].rh) {
					theta = 180;
				}
				else {
					theta = 90 - atan((Nodes[i][j][1].rh - Nodes[i][j][3].rh) / temp_angle) * 180 / pi;
				}

				if (theta > thetaA) {
					theta = thetaA;
					Nodes[i][j][1].rh = Nodes[i][j][3].rh + tan(pi * (90.0 - theta) / 180.0) * temp_angle;
				}
				else if (theta < thetaR) {
					theta = thetaR;
					Nodes[i][j][1].rh = Nodes[i][j][3].rh + tan(pi * (90.0 - theta) / 180.0) * temp_angle;
				}

				//上地面:全部设置为固壁边界条件
				/*
				temp_angle = sqrt(pow(Nodes[iright][j][h - 2].rh - Nodes[ileft][j][h - 2].rh, 2) + pow(Nodes[i][jbehind][h - 2].rh - Nodes[i][jfront][h - 2].rh, 2));
				if (temp_angle == 0 && Nodes[i][j][h - 1].rh >= Nodes[i][j][h - 3].rh) {
					theta = 0;
				}
				else if (temp_angle == 0 && Nodes[i][j][h - 1].rh < Nodes[i][j][h - 3].rh) {
					theta = 180;
				}
				else {
					theta = 90 - atan((Nodes[i][j][h - 1].rh - Nodes[i][j][h - 3].rh) / temp_angle) * 180 / pi;
				}

				if (theta > thetaA) {
					theta = thetaA;
					Nodes[i][j][h - 1].rh = Nodes[i][j][h - 3].rh + tan(pi * (90.0 - theta) / 180.0) * temp_angle;
				}
				else if (theta < thetaR) {
					theta = thetaR;
					Nodes[i][j][h - 1].rh = Nodes[i][j][h - 3].rh + tan(pi * (90.0 - theta) / 180.0) * temp_angle;
				}
				*/

				Nodes[i][j][h - 1].rh = Nodes[i][j][h - 3].rh;

				//上下地面的一层虚网格
				Nodes[i][j][0].rh = Nodes[i][j][1].rh;
				Nodes[i][j][h].rh = Nodes[i][j][h - 1].rh;

				//辅助变量prho和fai的上下固壁边界插值计算
				Nodes[i][j][h - 1].prho = Nodes[i][j][h - 2].prho;
				Nodes[i][j][1].prho = Nodes[i][j][2].prho;
				Nodes[i][j][h - 1].fai = Nodes[i][j][h - 2].fai;
				Nodes[i][j][1].fai = Nodes[i][j][2].fai;
			}
		}
	}

	if (front_boundary == 'B' && behind_boundary == 'B')
	{
#pragma omp parallel for private(i,j,r,ileft,iright,jfront,jbehind,rup,rdown,temp_angle,theta)
		for (i = 0; i <= l; i++)//x方向是显然的周期性边界条件
		{
			ileft = (i > 0) ? (i - 1) : (l);
			iright = (i < l) ? (i + 1) : (0);
			for (r = 1; r < h; r++)
			{
				rdown = r - 1;
				rup = r + 1;

				/*
				//前表面:全部设置为固壁边界条件
				temp_angle = sqrt(pow(Nodes[iright][2][r].rh - Nodes[ileft][2][r].rh, 2) + pow(Nodes[i][2][rup].rh - Nodes[i][2][rdown].rh, 2));
				if (temp_angle == 0 && Nodes[i][1][r].rh >= Nodes[i][3][r].rh)
				{
				theta = 0;
				}
				else if (temp_angle == 0 && Nodes[i][1][r].rh < Nodes[i][3][r].rh)
				{
				theta = 180;
				}
				else
				{
				theta = 90 - atan((Nodes[i][1][r].rh - Nodes[i][3][r].rh) / temp_angle) * 180 / pi;
				}

				if (theta > thetaA)
				{
				theta = thetaA;
				Nodes[i][1][r].rh = Nodes[i][3][r].rh + tan(pi*(90.0 - theta) / 180.0)*temp_angle;
				}
				else if (theta < thetaR)
				{
				theta = thetaR;
				Nodes[i][1][r].rh = Nodes[i][3][r].rh + tan(pi*(90.0 - theta) / 180.0)*temp_angle;
				}

				//后表面:全部设置为固壁边界条件
				temp_angle = sqrt(pow(Nodes[iright][m - 2][r].rh - Nodes[ileft][m - 2][r].rh, 2) + pow(Nodes[i][m - 2][rup].rh - Nodes[i][m - 2][rdown].rh, 2));
				if (temp_angle == 0 && Nodes[i][m - 1][r].rh >= Nodes[i][m - 3][r].rh)
				{
				theta = 0;
				}
				else if (temp_angle == 0 && Nodes[i][m - 1][r].rh < Nodes[i][m - 3][r].rh)
				{
				theta = 180;
				}
				else
				{
				theta = 90 - atan((Nodes[i][m - 1][r].rh - Nodes[i][m - 3][r].rh) / temp_angle) * 180 / pi;
				}

				if (theta > thetaA)
				{
				theta = thetaA;
				Nodes[i][m - 1][r].rh = Nodes[i][m - 3][r].rh + tan(pi*(90.0 - theta) / 180.0)*temp_angle;
				}
				else if (theta < thetaR)
				{
				theta = thetaR;
				Nodes[i][m - 1][r].rh = Nodes[i][m - 3][r].rh + tan(pi*(90.0 - theta) / 180.0)*temp_angle;
				}
				*/

				//前后表面的一层虚网格
				Nodes[i][0][r].rh = Nodes[i][1][r].rh;
				Nodes[i][m][r].rh = Nodes[i][m - 1][r].rh;

				//辅助变量prho和fai的前后固壁边界插值计算
				Nodes[i][m - 1][r].prho = Nodes[i][m - 2][r].prho;
				Nodes[i][1][r].prho = Nodes[i][2][r].prho;
				Nodes[i][m - 1][r].fai = Nodes[i][m - 2][r].fai;
				Nodes[i][1][r].fai = Nodes[i][2][r].fai;
			}
		}
	}


	if (front_boundary == 'D' && behind_boundary == 'D')
	{
#pragma omp parallel for private(i,r)
		for (i = 0; i <= l; i++)//x方向是显然的周期性边界条件
		{
			for (r = 1; r < h; r++)
			{
				//辅助变量prho和fai的前后固壁边界插值计算
				Nodes[i][m - 1][r].prho = Nodes[i][m - 2][r].prho;
				Nodes[i][1][r].prho = Nodes[i][2][r].prho;
				Nodes[i][m - 1][r].fai = Nodes[i][m - 2][r].fai;
				Nodes[i][1][r].fai = Nodes[i][2][r].fai;
			}
		}
	}


	//虚网格上的线处理：rh
#pragma omp parallel for private(i)
	for (i = 0; i <= l; i++)//x方向是显然的周期性边界条件
	{
		Nodes[i][0][0].rh = 0.5 * (Nodes[i][0][1].rh + Nodes[i][1][0].rh);
		Nodes[i][0][h].rh = 0.5 * (Nodes[i][1][h].rh + Nodes[i][0][h - 1].rh);
		Nodes[i][m][0].rh = 0.5 * (Nodes[i][m][1].rh + Nodes[i][m - 1][0].rh);
		Nodes[i][m][h].rh = 0.5 * (Nodes[i][m][h - 1].rh + Nodes[i][m - 1][h].rh);
	}

	//计算laplace_rh
#pragma omp parallel for private(i,j,r,ileft,iright,jbehind,jfront,rdown,rup)
	for (i = 0; i <= l; i++) {
		ileft = (i > 0) ? (i - 1) : (l);
		iright = (i < l) ? (i + 1) : (0);
		for (j = 1; j < m; j++) {
			jfront = j - 1;
			jbehind = j + 1;
			for (r = 1; r < h; r++) {
				rdown = r - 1;
				rup = r + 1;

				Nodes[i][j][r].laplace_rh = (Nodes[ileft][j][r].rh + Nodes[iright][j][r].rh + Nodes[i][jfront][r].rh + Nodes[i][jbehind][r].rh + Nodes[i][j][rup].rh + Nodes[i][j][rdown].rh) * 2.0 / 6.0;
				Nodes[i][j][r].laplace_rh += (Nodes[ileft][jfront][r].rh + Nodes[iright][jfront][r].rh + Nodes[ileft][jbehind][r].rh + Nodes[iright][jbehind][r].rh) / 6.0;
				Nodes[i][j][r].laplace_rh += (Nodes[ileft][j][rdown].rh + Nodes[iright][j][rdown].rh + Nodes[ileft][j][rup].rh + Nodes[iright][j][rup].rh) / 6.0;
				Nodes[i][j][r].laplace_rh += (Nodes[i][jfront][rdown].rh + Nodes[i][jbehind][rdown].rh + Nodes[i][jfront][rup].rh + Nodes[i][jbehind][rup].rh) / 6.0;
				Nodes[i][j][r].laplace_rh -= 24.0 * Nodes[i][j][r].rh / 6.0;
			}
		}
	}
	//the end of 3-D circulation.

	if (front_boundary == 'D' && behind_boundary == 'D') {
#pragma omp parallel for private(i,r)
		for (i = 0; i <= l; i++) {
			for (r = 1; r < h; r++) {
				Nodes[i][0][r].laplace_rh = Nodes[i][1][r].laplace_rh;
				Nodes[i][m][r].laplace_rh = Nodes[i][m - 1][r].laplace_rh;
			}
		}
	}

	//在前后虚网格上铺开fai和prho
	if (front_boundary == 'D' && behind_boundary == 'D') {
#pragma omp parallel for private(i,r)
		for (i = 0; i <= l; i++) {
			for (r = 1; r < h; r++) {
				Nodes[i][0][r].fai = Nodes[i][1][r].fai;
				Nodes[i][m][r].fai = Nodes[i][m - 1][r].fai;
				Nodes[i][0][r].prho = Nodes[i][1][r].prho;
				Nodes[i][m][r].prho = Nodes[i][m - 1][r].prho;
			}
		}
		//end of 2-Dcirculation.
	}
	//end if.

#pragma omp parallel for private(i,j,r,ileft,iright,jbehind,jfront,rdown,rup)
	for (i = 0; i <= l; i++) {
		ileft = (i > 0) ? (i - 1) : (l);
		iright = (i < l) ? (i + 1) : (0);
		for (j = 1; j < m; j++) {
			jfront = j - 1;
			jbehind = j + 1;
			for (r = 1; r < h; r++) {
				rdown = r - 1;
				rup = r + 1;

				if (Nodes[i][j][r].status[0] == 'B') { continue; }

				Nodes[i][j][r].Force[0] = 2.0 * (Nodes[iright][j][r].laplace_rh - Nodes[ileft][j][r].laplace_rh);
				Nodes[i][j][r].Force[0] += (Nodes[iright][jbehind][r].laplace_rh - Nodes[ileft][jfront][r].laplace_rh);
				Nodes[i][j][r].Force[0] += (Nodes[iright][jfront][r].laplace_rh - Nodes[ileft][jbehind][r].laplace_rh);
				Nodes[i][j][r].Force[0] += (Nodes[iright][j][rup].laplace_rh - Nodes[ileft][j][rdown].laplace_rh);
				Nodes[i][j][r].Force[0] += (Nodes[iright][j][rdown].laplace_rh - Nodes[ileft][j][rup].laplace_rh);
				Nodes[i][j][r].Force[0] *= Kappa * Nodes[i][j][r].rh / 12.0;

				Nodes[i][j][r].Force[1] = 2.0 * (Nodes[i][jbehind][r].laplace_rh - Nodes[i][jfront][r].laplace_rh);
				Nodes[i][j][r].Force[1] += (Nodes[iright][jbehind][r].laplace_rh - Nodes[ileft][jfront][r].laplace_rh);
				Nodes[i][j][r].Force[1] += (Nodes[ileft][jbehind][r].laplace_rh - Nodes[iright][jfront][r].laplace_rh);
				Nodes[i][j][r].Force[1] += (Nodes[i][jbehind][rup].laplace_rh - Nodes[i][jfront][rdown].laplace_rh);
				Nodes[i][j][r].Force[1] += (Nodes[i][jbehind][rdown].laplace_rh - Nodes[i][jfront][rup].laplace_rh);
				Nodes[i][j][r].Force[1] *= Kappa * Nodes[i][j][r].rh / 12.0;

				Nodes[i][j][r].Force[2] = 2.0 * (Nodes[i][j][rup].laplace_rh - Nodes[i][j][rdown].laplace_rh);
				Nodes[i][j][r].Force[2] += (Nodes[iright][j][rup].laplace_rh - Nodes[ileft][j][rdown].laplace_rh);
				Nodes[i][j][r].Force[2] += (Nodes[ileft][j][rup].laplace_rh - Nodes[iright][j][rdown].laplace_rh);
				Nodes[i][j][r].Force[2] += (Nodes[i][jbehind][rup].laplace_rh - Nodes[i][jfront][rdown].laplace_rh);
				Nodes[i][j][r].Force[2] += (Nodes[i][jfront][rup].laplace_rh - Nodes[i][jbehind][rdown].laplace_rh);
				Nodes[i][j][r].Force[2] *= Kappa * Nodes[i][j][r].rh / 12.0;

				if (BUOYANCY == 1) {
					Nodes[i][j][r].Force[2] -= gforce * Nodes[i][j][r].rho;
				}

				Nodes[i][j][r].dfai[0] = 2.0 * (Nodes[iright][j][r].fai - Nodes[ileft][j][r].fai);
				Nodes[i][j][r].dfai[0] += (Nodes[iright][jbehind][r].fai - Nodes[ileft][jfront][r].fai);
				Nodes[i][j][r].dfai[0] += (Nodes[iright][jfront][r].fai - Nodes[ileft][jbehind][r].fai);
				Nodes[i][j][r].dfai[0] += (Nodes[iright][j][rup].fai - Nodes[ileft][j][rdown].fai);
				Nodes[i][j][r].dfai[0] += (Nodes[iright][j][rdown].fai - Nodes[ileft][j][rup].fai);
				Nodes[i][j][r].dfai[0] /= 12.0;

				Nodes[i][j][r].dfai[1] = 2.0 * (Nodes[i][jbehind][r].fai - Nodes[i][jfront][r].fai);
				Nodes[i][j][r].dfai[1] += (Nodes[iright][jbehind][r].fai - Nodes[ileft][jfront][r].fai);
				Nodes[i][j][r].dfai[1] += (Nodes[ileft][jbehind][r].fai - Nodes[iright][jfront][r].fai);
				Nodes[i][j][r].dfai[1] += (Nodes[i][jbehind][rup].fai - Nodes[i][jfront][rdown].fai);
				Nodes[i][j][r].dfai[1] += (Nodes[i][jbehind][rdown].fai - Nodes[i][jfront][rup].fai);
				Nodes[i][j][r].dfai[1] /= 12.0;

				Nodes[i][j][r].dfai[2] = 2.0 * (Nodes[i][j][rup].fai - Nodes[i][j][rdown].fai);
				Nodes[i][j][r].dfai[2] += (Nodes[iright][j][rup].fai - Nodes[ileft][j][rdown].fai);
				Nodes[i][j][r].dfai[2] += (Nodes[ileft][j][rup].fai - Nodes[iright][j][rdown].fai);
				Nodes[i][j][r].dfai[2] += (Nodes[i][jbehind][rup].fai - Nodes[i][jfront][rdown].fai);
				Nodes[i][j][r].dfai[2] += (Nodes[i][jfront][rup].fai - Nodes[i][jbehind][rdown].fai);
				Nodes[i][j][r].dfai[2] /= 12.0;

				Nodes[i][j][r].dprho[0] = 2.0 * (Nodes[iright][j][r].prho - Nodes[ileft][j][r].prho);
				Nodes[i][j][r].dprho[0] += (Nodes[iright][jbehind][r].prho - Nodes[ileft][jfront][r].prho);
				Nodes[i][j][r].dprho[0] += (Nodes[iright][jfront][r].prho - Nodes[ileft][jbehind][r].prho);
				Nodes[i][j][r].dprho[0] += (Nodes[iright][j][rup].prho - Nodes[ileft][j][rdown].prho);
				Nodes[i][j][r].dprho[0] += (Nodes[iright][j][rdown].prho - Nodes[ileft][j][rup].prho);
				Nodes[i][j][r].dprho[0] /= 12.0;

				Nodes[i][j][r].dprho[1] = 2.0 * (Nodes[i][jbehind][r].prho - Nodes[i][jfront][r].prho);
				Nodes[i][j][r].dprho[1] += (Nodes[iright][jbehind][r].prho - Nodes[ileft][jfront][r].prho);
				Nodes[i][j][r].dprho[1] += (Nodes[ileft][jbehind][r].prho - Nodes[iright][jfront][r].prho);
				Nodes[i][j][r].dprho[1] += (Nodes[i][jbehind][rup].prho - Nodes[i][jfront][rdown].prho);
				Nodes[i][j][r].dprho[1] += (Nodes[i][jbehind][rdown].prho - Nodes[i][jfront][rup].prho);
				Nodes[i][j][r].dprho[1] /= 12.0;

				Nodes[i][j][r].dprho[2] = 2.0 * (Nodes[i][j][rup].prho - Nodes[i][j][rdown].prho);
				Nodes[i][j][r].dprho[2] += (Nodes[iright][j][rup].prho - Nodes[ileft][j][rdown].prho);
				Nodes[i][j][r].dprho[2] += (Nodes[ileft][j][rup].prho - Nodes[iright][j][rdown].prho);
				Nodes[i][j][r].dprho[2] += (Nodes[i][jbehind][rup].prho - Nodes[i][jfront][rdown].prho);
				Nodes[i][j][r].dprho[2] += (Nodes[i][jfront][rup].prho - Nodes[i][jbehind][rdown].prho);
				Nodes[i][j][r].dprho[2] /= 12.0;
			}
		}
	}
	//the end of 3-D circulation.
	return;
}

void LBM_D3Q19_HCZ::geten()
{
	GLint i = 0, j = 0, r = 0, rr = 0, k = 0;

#pragma omp parallel for private(i,j,r,rr,k)
	for (i = 0; i <= l; i++)
	{
		for (j = 1; j < m; j++)
		{
			for (r = 1; r < h; r++)
			{
				if (Nodes[i][j][r].status == "L")//压力边界点处的速度由插值得到，故不作计算
				{
					for (rr = 0; rr < D; rr++)
					{
						Nodes[i][j][r].u[rr] = 0.0;
						for (k = 0; k < Q; k++)
						{
							Nodes[i][j][r].u[rr] += e[k][rr] * Nodes[i][j][r].g[k];
						}
						Nodes[i][j][r].u[rr] *= c;
						Nodes[i][j][r].u[rr] += 0.5 * RT * dt * Nodes[i][j][r].Force[rr];
						Nodes[i][j][r].u[rr] /= (Nodes[i][j][r].rho * RT);
					}
				}
				//end if.

				if (Nodes[i][j][r].status == "L" || Nodes[i][j][r].status[0] == 'D')
				{
					Nodes[i][j][r].p = 0.0;
					for (k = 0; k < Q; k++)
					{
						Nodes[i][j][r].p += Nodes[i][j][r].g[k];
					}
					for (rr = 0; rr < D; rr++)
					{
						Nodes[i][j][r].p += (-0.5) * Nodes[i][j][r].u[rr] * Nodes[i][j][r].dprho[rr] * dt;
					}
				}
				//end if.
			}
		}
	}
	return;
}

void LBM_D3Q19_HCZ::Collision()
{
	GLint i = 0, j = 0, r = 0, k = 0, rr = 0;
	GLdouble temp_f = 0.0, temp_g1 = 0.0, temp_g2 = 0.0;

	//Collide process
#pragma omp parallel for private(i,j,r,rr,k,temp_f,temp_g1,temp_g2)
	for (i = 0; i <= l; i++)
	{
		for (j = 1; j < m; j++)
		{
			for (r = 1; r < h; r++)
			{
				if (Nodes[i][j][r].status == "L")// || Nodes[i][j][r].status[0] == 'D')
				{
					calculate_feq_geq(Nodes[i][j][r].rh, Nodes[i][j][r].rho, Nodes[i][j][r].p, Nodes[i][j][r].u, Nodes[i][j][r].feq, Nodes[i][j][r].geq);
					for (k = 0; k < Q; k++)
					{
						Nodes[i][j][r].f[k] = Nodes[i][j][r].f[k] * (1.0 - 1.0 / tau_f) + Nodes[i][j][r].feq[k] / tau_f;
						Nodes[i][j][r].g[k] = Nodes[i][j][r].g[k] * (1.0 - 1.0 / tau_g) + Nodes[i][j][r].geq[k] / tau_g;

						temp_f = 0.0;
						temp_g1 = 0.0;
						temp_g2 = 0.0;

						for (rr = 0; rr < D; rr++)
						{
							temp_f += (e[k][rr] * c - Nodes[i][j][r].u[rr]) * (-Nodes[i][j][r].dfai[rr]) * Nodes[i][j][r].feq[k];
							temp_g1 += (e[k][rr] * c - Nodes[i][j][r].u[rr]) * Nodes[i][j][r].Force[rr] * Nodes[i][j][r].feq[k] / Nodes[i][j][r].rh;
							temp_g2 += (e[k][rr] * c - Nodes[i][j][r].u[rr]) * (-Nodes[i][j][r].dprho[rr]) * ((Nodes[i][j][r].feq[k] / Nodes[i][j][r].rh) - w[k]);
						}
						Nodes[i][j][r].f[k] += dt * (tau_f - 0.5) / tau_f * temp_f / (RT * Nodes[i][j][r].rh);
						Nodes[i][j][r].g[k] += dt * (tau_g - 0.5) / tau_g * (temp_g1 + temp_g2);
					}
				}
				//end if.
			}
		}
	}
	//the end of 3-D circulation.
	return;
}

void LBM_D3Q19_HCZ::drawing_liquid(void)
{
	int i = 0, j = 0, r = 0;
	int hi = 0, hj = 0;
	double rho_ave = 0.5 * (rho_l + rho_h);

	glColor3f(0.0, 0.0, 1.0);
	glPointSize(0.00001);
	glBegin(GL_POINTS);
	for (i = 0; i <= l; i++)
	{
		for (j = 0; j <= m; j++)
		{
			for (r = 0; r <= h; r++)
			{
				if (Nodes[i][j][r].status == "L")
				{
					if (Nodes[i][j][r].rho >= rho_ave)
					{
						glVertex3f(i, j, r);
					}
				}
			}
		}
	}
	glEnd();

	glColor3f(0.0, 1.0, 0.0);
	//绘制流场外围
	glBegin(GL_LINE_STRIP);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(l, 0.0, 0.0);
	glVertex3f(l, m, 0.0);
	glVertex3f(0.0, m, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glEnd();
	glBegin(GL_LINE_STRIP);
	glVertex3f(0.0, 0.0, h);
	glVertex3f(l, 0.0, h);
	glVertex3f(l, m, h);
	glVertex3f(0.0, m, h);
	glVertex3f(0.0, 0.0, h);
	glEnd();
	glBegin(GL_LINES);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, h);
	glVertex3f(l, 0.0, 0.0);
	glVertex3f(l, 0.0, h);
	glVertex3f(0.0, m, 0.0);
	glVertex3f(0.0, m, h);
	glVertex3f(l, m, 0.0);
	glVertex3f(l, m, h);
	glEnd();
	//绘制流场外围结束

	return;
}

void LBM_D3Q19_HCZ::dump_file(string filename)
{
	int i = 0, j = 0, r = 0;
	FILE* fp;
	if (!(fp = fopen((filename).c_str(), "w")))
	{
		printf("写入文件失败!\n");
		exit(1);
	}
	fprintf(fp, "TITLE =\"lattice boltzmann\"\n");
	fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\" \"U\" \"V\" \"W\" \"RHO\" \n");
	fprintf(fp, "ZONE I=%d, J=%d, K=%d\n", l + 1, m - 1, h - 3);
	fprintf(fp, "F=POINT\n");
	for (r = 2; r < h - 1; r++) //B
	{
		for (j = 1; j < m; j++) //P
		{
			for (i = 0; i <= l; i++) //P
			{
				fprintf(fp, "%d %d %d %lf %lf %lf %lf\n", i, j, r, Nodes[i][j][r].u[0], Nodes[i][j][r].u[1], Nodes[i][j][r].u[2], Nodes[i][j][r].rho);
			}
		}
	}
	fclose(fp);
}

#endif // use_HCZ == 1