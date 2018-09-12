//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is a serial C version of the NPB BT code. This C        //
//  version is developed by the Center for Manycore Programming at Seoul   //
//  National University and derived from the serial Fortran versions in    //
//  "NPB3.3-SER" developed by NAS.                                         //
//                                                                         //
//  Permission to use, copy, distribute and modify this software for any   //
//  purpose with or without fee is hereby granted. This software is        //
//  provided "as is" without express or implied warranty.                  //
//                                                                         //
//  Information on NPB 3.3, including the technical report, the original   //
//  specifications, source code, results and information on how to submit  //
//  new results, is available at:                                          //
//                                                                         //
//           http://www.nas.nasa.gov/Software/NPB/                         //
//                                                                         //
//  Send comments or suggestions for this C version to cmp@aces.snu.ac.kr  //
//                                                                         //
//          Center for Manycore Programming                                //
//          School of Computer Science and Engineering                     //
//          Seoul National University                                      //
//          Seoul 151-744, Korea                                           //
//                                                                         //
//          E-mail:  cmp@aces.snu.ac.kr                                    //
//                                                                         //
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
// Authors: Sangmin Seo, Jungwon Kim, Jun Lee, Jeongho Nah, Gangwon Jo,    //
//          and Jaejin Lee                                                 //
//-------------------------------------------------------------------------//

#include "header.h"
#include "work_lhs.h"
#include "timers.h"


double lhs_y[2][3][5][5];
double rhs_y[2][5];
int len_y = 4;
double fjac_y[4][5][5];
double njac_y[4][5][5];

void fn_init_y(int a, int i, int j, int k);
void buf_compute_y(int j, int a);
void binvcrhs_y(int a1, int a2, int b1, int b2, int c2);
void matvec_sub_y(int a1, int a2, int b2, int c2);
void matmul_sub_y(int a1, int a2, int b1, int b2, int c1, int c2);
void binvrhs_y(int a1, int a2, int b2);

//---------------------------------------------------------------------
// Performs line solves in Y direction by first factoring
// the block-tridiagonal matrix into an upper triangular matrix, 
// and then performing back substitution to solve for the unknow
// vectors of each line.  
// 
// Make sure we treat elements zero to cell_size in the direction
// of the sweep.
//---------------------------------------------------------------------
void y_solve()
{
	int i, j, k, m, n, jsize, l;
	double pivot, coeff;
	//int z, a = -1, b = 2;
	//---------------------------------------------------------------------
	//---------------------------------------------------------------------

	if (timeron) timer_start(t_ysolve);

	//---------------------------------------------------------------------
	//---------------------------------------------------------------------

	//---------------------------------------------------------------------
	// This function computes the left hand side for the three y-factors   
	//---------------------------------------------------------------------

	jsize = grid_points[1] - 1;

	//---------------------------------------------------------------------
	// Compute the indices for storing the tri-diagonal matrix;
	// determine a (labeled f) and n jacobians for cell c
	//---------------------------------------------------------------------

	


	for (k = 1; k <= grid_points[2] - 2; k++) {
		for (i = 1; i <= grid_points[0] - 2; i++) {
			for (j = 0; j <= jsize; j++) {
				tmp1 = rho_i[k][j][i];
				tmp2 = tmp1 * tmp1;
				tmp3 = tmp1 * tmp2;

				fjac[j][0][0] = 0.0;
				fjac[j][1][0] = 0.0;
				fjac[j][2][0] = 1.0;
				fjac[j][3][0] = 0.0;
				fjac[j][4][0] = 0.0;

				fjac[j][0][1] = -(u[k][j][i][1] * u[k][j][i][2]) * tmp2;
				fjac[j][1][1] = u[k][j][i][2] * tmp1;
				fjac[j][2][1] = u[k][j][i][1] * tmp1;
				fjac[j][3][1] = 0.0;
				fjac[j][4][1] = 0.0;

				fjac[j][0][2] = -(u[k][j][i][2] * u[k][j][i][2] * tmp2)
					+ c2 * qs[k][j][i];
				fjac[j][1][2] = -c2 *  u[k][j][i][1] * tmp1;
				fjac[j][2][2] = (2.0 - c2) *  u[k][j][i][2] * tmp1;
				fjac[j][3][2] = -c2 * u[k][j][i][3] * tmp1;
				fjac[j][4][2] = c2;

				fjac[j][0][3] = -(u[k][j][i][2] * u[k][j][i][3]) * tmp2;
				fjac[j][1][3] = 0.0;
				fjac[j][2][3] = u[k][j][i][3] * tmp1;
				fjac[j][3][3] = u[k][j][i][2] * tmp1;
				fjac[j][4][3] = 0.0;

				fjac[j][0][4] = (c2 * 2.0 * square[k][j][i] - c1 * u[k][j][i][4])
					* u[k][j][i][2] * tmp2;
				fjac[j][1][4] = -c2 * u[k][j][i][1] * u[k][j][i][2] * tmp2;
				fjac[j][2][4] = c1 * u[k][j][i][4] * tmp1
					- c2 * (qs[k][j][i] + u[k][j][i][2] * u[k][j][i][2] * tmp2);
				fjac[j][3][4] = -c2 * (u[k][j][i][2] * u[k][j][i][3]) * tmp2;
				fjac[j][4][4] = c1 * u[k][j][i][2] * tmp1;

				njac[j][0][0] = 0.0;
				njac[j][1][0] = 0.0;
				njac[j][2][0] = 0.0;
				njac[j][3][0] = 0.0;
				njac[j][4][0] = 0.0;

				njac[j][0][1] = -c3c4 * tmp2 * u[k][j][i][1];
				njac[j][1][1] = c3c4 * tmp1;
				njac[j][2][1] = 0.0;
				njac[j][3][1] = 0.0;
				njac[j][4][1] = 0.0;

				njac[j][0][2] = -con43 * c3c4 * tmp2 * u[k][j][i][2];
				njac[j][1][2] = 0.0;
				njac[j][2][2] = con43 * c3c4 * tmp1;
				njac[j][3][2] = 0.0;
				njac[j][4][2] = 0.0;

				njac[j][0][3] = -c3c4 * tmp2 * u[k][j][i][3];
				njac[j][1][3] = 0.0;
				njac[j][2][3] = 0.0;
				njac[j][3][3] = c3c4 * tmp1;
				njac[j][4][3] = 0.0;

				njac[j][0][4] = -(c3c4
					- c1345) * tmp3 * (u[k][j][i][1] * u[k][j][i][1])
					- (con43 * c3c4
						- c1345) * tmp3 * (u[k][j][i][2] * u[k][j][i][2])
					- (c3c4 - c1345) * tmp3 * (u[k][j][i][3] * u[k][j][i][3])
					- c1345 * tmp2 * u[k][j][i][4];

				njac[j][1][4] = (c3c4 - c1345) * tmp2 * u[k][j][i][1];
				njac[j][2][4] = (con43 * c3c4 - c1345) * tmp2 * u[k][j][i][2];
				njac[j][3][4] = (c3c4 - c1345) * tmp2 * u[k][j][i][3];
				njac[j][4][4] = (c1345)* tmp1;

			}


			for (j = 0; j <= jsize; j++) {
				
				if(j == 0)
				{
					for(n = 0; n < 5; n++)
					{
						for(m = 0; m < 5; m++)
						{
							lhs_buf[k][j][i][AA][n][m] = 0.0;
							lhs_buf[k][j][i][BB][n][m] = 0.0;
							lhs_buf[k][j][i][CC][n][m] = 0.0;
						}
					}
					for(n = 0; n < 5; n++)
					{
						lhs_buf[k][j][i][BB][n][n] = 1.0;
					}
				}
				else if(j == jsize)
				{
					for(n = 0; n < 5; n++)
					{
						for(m = 0; m < 5; m++)
						{
							lhs_buf[k][j][i][AA][n][m] = 0.0;
							lhs_buf[k][j][i][BB][n][m] = 0.0;
							lhs_buf[k][j][i][CC][n][m] = 0.0;
						}
					}
					for(n = 0; n < 5; n++)
					{
						lhs_buf[k][j][i][BB][n][n] = 1.0;
					}
				}
				else
				{
					tmp1 = dt * ty1;
					tmp2 = dt * ty2;

					lhs_buf[k][j][i][AA][0][0] = -tmp2 * fjac[j - 1][0][0]
						- tmp1 * njac[j - 1][0][0]
						- tmp1 * dy1;
					lhs_buf[k][j][i][AA][1][0] = -tmp2 * fjac[j - 1][1][0]
						- tmp1 * njac[j - 1][1][0];
					lhs_buf[k][j][i][AA][2][0] = -tmp2 * fjac[j - 1][2][0]
						- tmp1 * njac[j - 1][2][0];
					lhs_buf[k][j][i][AA][3][0] = -tmp2 * fjac[j - 1][3][0]
						- tmp1 * njac[j - 1][3][0];
					lhs_buf[k][j][i][AA][4][0] = -tmp2 * fjac[j - 1][4][0]
						- tmp1 * njac[j - 1][4][0];

					lhs_buf[k][j][i][AA][0][1] = -tmp2 * fjac[j - 1][0][1]
						- tmp1 * njac[j - 1][0][1];
					lhs_buf[k][j][i][AA][1][1] = -tmp2 * fjac[j - 1][1][1]
						- tmp1 * njac[j - 1][1][1]
						- tmp1 * dy2;
					lhs_buf[k][j][i][AA][2][1] = -tmp2 * fjac[j - 1][2][1]
						- tmp1 * njac[j - 1][2][1];
					lhs_buf[k][j][i][AA][3][1] = -tmp2 * fjac[j - 1][3][1]
						- tmp1 * njac[j - 1][3][1];
					lhs_buf[k][j][i][AA][4][1] = -tmp2 * fjac[j - 1][4][1]
						- tmp1 * njac[j - 1][4][1];

					lhs_buf[k][j][i][AA][0][2] = -tmp2 * fjac[j - 1][0][2]
						- tmp1 * njac[j - 1][0][2];
					lhs_buf[k][j][i][AA][1][2] = -tmp2 * fjac[j - 1][1][2]
						- tmp1 * njac[j - 1][1][2];
					lhs_buf[k][j][i][AA][2][2] = -tmp2 * fjac[j - 1][2][2]
						- tmp1 * njac[j - 1][2][2]
						- tmp1 * dy3;
					lhs_buf[k][j][i][AA][3][2] = -tmp2 * fjac[j - 1][3][2]
						- tmp1 * njac[j - 1][3][2];
					lhs_buf[k][j][i][AA][4][2] = -tmp2 * fjac[j - 1][4][2]
						- tmp1 * njac[j - 1][4][2];

					lhs_buf[k][j][i][AA][0][3] = -tmp2 * fjac[j - 1][0][3]
						- tmp1 * njac[j - 1][0][3];
					lhs_buf[k][j][i][AA][1][3] = -tmp2 * fjac[j - 1][1][3]
						- tmp1 * njac[j - 1][1][3];
					lhs_buf[k][j][i][AA][2][3] = -tmp2 * fjac[j - 1][2][3]
						- tmp1 * njac[j - 1][2][3];
					lhs_buf[k][j][i][AA][3][3] = -tmp2 * fjac[j - 1][3][3]
						- tmp1 * njac[j - 1][3][3]
						- tmp1 * dy4;
					lhs_buf[k][j][i][AA][4][3] = -tmp2 * fjac[j - 1][4][3]
						- tmp1 * njac[j - 1][4][3];

					lhs_buf[k][j][i][AA][0][4] = -tmp2 * fjac[j - 1][0][4]
						- tmp1 * njac[j - 1][0][4];
					lhs_buf[k][j][i][AA][1][4] = -tmp2 * fjac[j - 1][1][4]
						- tmp1 * njac[j - 1][1][4];
					lhs_buf[k][j][i][AA][2][4] = -tmp2 * fjac[j - 1][2][4]
						- tmp1 * njac[j - 1][2][4];
					lhs_buf[k][j][i][AA][3][4] = -tmp2 * fjac[j - 1][3][4]
						- tmp1 * njac[j - 1][3][4];
					lhs_buf[k][j][i][AA][4][4] = -tmp2 * fjac[j - 1][4][4]
						- tmp1 * njac[j - 1][4][4]
						- tmp1 * dy5;

					lhs_buf[k][j][i][BB][0][0] = 1.0
						+ tmp1 * 2.0 * njac[j][0][0]
						+ tmp1 * 2.0 * dy1;
					lhs_buf[k][j][i][BB][1][0] = tmp1 * 2.0 * njac[j][1][0];
					lhs_buf[k][j][i][BB][2][0] = tmp1 * 2.0 * njac[j][2][0];
					lhs_buf[k][j][i][BB][3][0] = tmp1 * 2.0 * njac[j][3][0];
					lhs_buf[k][j][i][BB][4][0] = tmp1 * 2.0 * njac[j][4][0];

					lhs_buf[k][j][i][BB][0][1] = tmp1 * 2.0 * njac[j][0][1];
					lhs_buf[k][j][i][BB][1][1] = 1.0
						+ tmp1 * 2.0 * njac[j][1][1]
						+ tmp1 * 2.0 * dy2;
					lhs_buf[k][j][i][BB][2][1] = tmp1 * 2.0 * njac[j][2][1];
					lhs_buf[k][j][i][BB][3][1] = tmp1 * 2.0 * njac[j][3][1];
					lhs_buf[k][j][i][BB][4][1] = tmp1 * 2.0 * njac[j][4][1];

					lhs_buf[k][j][i][BB][0][2] = tmp1 * 2.0 * njac[j][0][2];
					lhs_buf[k][j][i][BB][1][2] = tmp1 * 2.0 * njac[j][1][2];
					lhs_buf[k][j][i][BB][2][2] = 1.0
						+ tmp1 * 2.0 * njac[j][2][2]
						+ tmp1 * 2.0 * dy3;
					lhs_buf[k][j][i][BB][3][2] = tmp1 * 2.0 * njac[j][3][2];
					lhs_buf[k][j][i][BB][4][2] = tmp1 * 2.0 * njac[j][4][2];

					lhs_buf[k][j][i][BB][0][3] = tmp1 * 2.0 * njac[j][0][3];
					lhs_buf[k][j][i][BB][1][3] = tmp1 * 2.0 * njac[j][1][3];
					lhs_buf[k][j][i][BB][2][3] = tmp1 * 2.0 * njac[j][2][3];
					lhs_buf[k][j][i][BB][3][3] = 1.0
						+ tmp1 * 2.0 * njac[j][3][3]
						+ tmp1 * 2.0 * dy4;
					lhs_buf[k][j][i][BB][4][3] = tmp1 * 2.0 * njac[j][4][3];

					lhs_buf[k][j][i][BB][0][4] = tmp1 * 2.0 * njac[j][0][4];
					lhs_buf[k][j][i][BB][1][4] = tmp1 * 2.0 * njac[j][1][4];
					lhs_buf[k][j][i][BB][2][4] = tmp1 * 2.0 * njac[j][2][4];
					lhs_buf[k][j][i][BB][3][4] = tmp1 * 2.0 * njac[j][3][4];
					lhs_buf[k][j][i][BB][4][4] = 1.0
						+ tmp1 * 2.0 * njac[j][4][4]
						+ tmp1 * 2.0 * dy5;

					lhs_buf[k][j][i][CC][0][0] = tmp2 * fjac[j + 1][0][0]
						- tmp1 * njac[j + 1][0][0]
						- tmp1 * dy1;
					lhs_buf[k][j][i][CC][1][0] = tmp2 * fjac[j + 1][1][0]
						- tmp1 * njac[j + 1][1][0];
					lhs_buf[k][j][i][CC][2][0] = tmp2 * fjac[j + 1][2][0]
						- tmp1 * njac[j + 1][2][0];
					lhs_buf[k][j][i][CC][3][0] = tmp2 * fjac[j + 1][3][0]
						- tmp1 * njac[j + 1][3][0];
					lhs_buf[k][j][i][CC][4][0] = tmp2 * fjac[j + 1][4][0]
						- tmp1 * njac[j + 1][4][0];

					lhs_buf[k][j][i][CC][0][1] = tmp2 * fjac[j + 1][0][1]
						- tmp1 * njac[j + 1][0][1];
					lhs_buf[k][j][i][CC][1][1] = tmp2 * fjac[j + 1][1][1]
						- tmp1 * njac[j + 1][1][1]
						- tmp1 * dy2;
					lhs_buf[k][j][i][CC][2][1] = tmp2 * fjac[j + 1][2][1]
						- tmp1 * njac[j + 1][2][1];
					lhs_buf[k][j][i][CC][3][1] = tmp2 * fjac[j + 1][3][1]
						- tmp1 * njac[j + 1][3][1];
					lhs_buf[k][j][i][CC][4][1] = tmp2 * fjac[j + 1][4][1]
						- tmp1 * njac[j + 1][4][1];

					lhs_buf[k][j][i][CC][0][2] = tmp2 * fjac[j + 1][0][2]
						- tmp1 * njac[j + 1][0][2];
					lhs_buf[k][j][i][CC][1][2] = tmp2 * fjac[j + 1][1][2]
						- tmp1 * njac[j + 1][1][2];
					lhs_buf[k][j][i][CC][2][2] = tmp2 * fjac[j + 1][2][2]
						- tmp1 * njac[j + 1][2][2]
						- tmp1 * dy3;
					lhs_buf[k][j][i][CC][3][2] = tmp2 * fjac[j + 1][3][2]
						- tmp1 * njac[j + 1][3][2];
					lhs_buf[k][j][i][CC][4][2] = tmp2 * fjac[j + 1][4][2]
						- tmp1 * njac[j + 1][4][2];

					lhs_buf[k][j][i][CC][0][3] = tmp2 * fjac[j + 1][0][3]
						- tmp1 * njac[j + 1][0][3];
					lhs_buf[k][j][i][CC][1][3] = tmp2 * fjac[j + 1][1][3]
						- tmp1 * njac[j + 1][1][3];
					lhs_buf[k][j][i][CC][2][3] = tmp2 * fjac[j + 1][2][3]
						- tmp1 * njac[j + 1][2][3];
					lhs_buf[k][j][i][CC][3][3] = tmp2 * fjac[j + 1][3][3]
						- tmp1 * njac[j + 1][3][3]
						- tmp1 * dy4;
					lhs_buf[k][j][i][CC][4][3] = tmp2 * fjac[j + 1][4][3]
						- tmp1 * njac[j + 1][4][3];

					lhs_buf[k][j][i][CC][0][4] = tmp2 * fjac[j + 1][0][4]
						- tmp1 * njac[j + 1][0][4];
					lhs_buf[k][j][i][CC][1][4] = tmp2 * fjac[j + 1][1][4]
						- tmp1 * njac[j + 1][1][4];
					lhs_buf[k][j][i][CC][2][4] = tmp2 * fjac[j + 1][2][4]
						- tmp1 * njac[j + 1][2][4];
					lhs_buf[k][j][i][CC][3][4] = tmp2 * fjac[j + 1][3][4]
						- tmp1 * njac[j + 1][3][4];
					lhs_buf[k][j][i][CC][4][4] = tmp2 * fjac[j + 1][4][4]
						- tmp1 * njac[j + 1][4][4]
						- tmp1 * dy5;
				}
			}

			
		}
	}

	
	for (k = 1; k <= grid_points[2] - 2; k++) {
		for (i = 1; i <= grid_points[0] - 2; i++) {	
			for (j = 0; j <= jsize; j++)
			{
				if(j == 0)
				{
					pivot = 1.00/lhs_buf[k][j][i][BB][0][0];
					lhs_buf[k][j][i][BB][1][0] = lhs_buf[k][j][i][BB][1][0]*pivot;
					lhs_buf[k][j][i][BB][2][0] = lhs_buf[k][j][i][BB][2][0]*pivot;
					lhs_buf[k][j][i][BB][3][0] = lhs_buf[k][j][i][BB][3][0]*pivot;
					lhs_buf[k][j][i][BB][4][0] = lhs_buf[k][j][i][BB][4][0]*pivot;
					lhs_buf[k][j][i][CC][0][0] = lhs_buf[k][j][i][CC][0][0]*pivot;
					lhs_buf[k][j][i][CC][1][0] = lhs_buf[k][j][i][CC][1][0]*pivot;
					lhs_buf[k][j][i][CC][2][0] = lhs_buf[k][j][i][CC][2][0]*pivot;
					lhs_buf[k][j][i][CC][3][0] = lhs_buf[k][j][i][CC][3][0]*pivot;
					lhs_buf[k][j][i][CC][4][0] = lhs_buf[k][j][i][CC][4][0]*pivot;
					rhs[k][j][i][0]   = rhs[k][j][i][0]  *pivot;

					coeff = lhs_buf[k][j][i][BB][0][1];
					lhs_buf[k][j][i][BB][1][1]= lhs_buf[k][j][i][BB][1][1] - coeff*lhs_buf[k][j][i][BB][1][0];
					lhs_buf[k][j][i][BB][2][1]= lhs_buf[k][j][i][BB][2][1] - coeff*lhs_buf[k][j][i][BB][2][0];
					lhs_buf[k][j][i][BB][3][1]= lhs_buf[k][j][i][BB][3][1] - coeff*lhs_buf[k][j][i][BB][3][0];
					lhs_buf[k][j][i][BB][4][1]= lhs_buf[k][j][i][BB][4][1] - coeff*lhs_buf[k][j][i][BB][4][0];
					lhs_buf[k][j][i][CC][0][1] = lhs_buf[k][j][i][CC][0][1] - coeff*lhs_buf[k][j][i][CC][0][0];
					lhs_buf[k][j][i][CC][1][1] = lhs_buf[k][j][i][CC][1][1] - coeff*lhs_buf[k][j][i][CC][1][0];
					lhs_buf[k][j][i][CC][2][1] = lhs_buf[k][j][i][CC][2][1] - coeff*lhs_buf[k][j][i][CC][2][0];
					lhs_buf[k][j][i][CC][3][1] = lhs_buf[k][j][i][CC][3][1] - coeff*lhs_buf[k][j][i][CC][3][0];
					lhs_buf[k][j][i][CC][4][1] = lhs_buf[k][j][i][CC][4][1] - coeff*lhs_buf[k][j][i][CC][4][0];
					rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff*rhs[k][j][i][0];

					coeff = lhs_buf[k][j][i][BB][0][2];
					lhs_buf[k][j][i][BB][1][2]= lhs_buf[k][j][i][BB][1][2] - coeff*lhs_buf[k][j][i][BB][1][0];
					lhs_buf[k][j][i][BB][2][2]= lhs_buf[k][j][i][BB][2][2] - coeff*lhs_buf[k][j][i][BB][2][0];
					lhs_buf[k][j][i][BB][3][2]= lhs_buf[k][j][i][BB][3][2] - coeff*lhs_buf[k][j][i][BB][3][0];
					lhs_buf[k][j][i][BB][4][2]= lhs_buf[k][j][i][BB][4][2] - coeff*lhs_buf[k][j][i][BB][4][0];
					lhs_buf[k][j][i][CC][0][2] = lhs_buf[k][j][i][CC][0][2] - coeff*lhs_buf[k][j][i][CC][0][0];
					lhs_buf[k][j][i][CC][1][2] = lhs_buf[k][j][i][CC][1][2] - coeff*lhs_buf[k][j][i][CC][1][0];
					lhs_buf[k][j][i][CC][2][2] = lhs_buf[k][j][i][CC][2][2] - coeff*lhs_buf[k][j][i][CC][2][0];
					lhs_buf[k][j][i][CC][3][2] = lhs_buf[k][j][i][CC][3][2] - coeff*lhs_buf[k][j][i][CC][3][0];
					lhs_buf[k][j][i][CC][4][2] = lhs_buf[k][j][i][CC][4][2] - coeff*lhs_buf[k][j][i][CC][4][0];
					rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff*rhs[k][j][i][0];

					coeff = lhs_buf[k][j][i][BB][0][3];
					lhs_buf[k][j][i][BB][1][3]= lhs_buf[k][j][i][BB][1][3] - coeff*lhs_buf[k][j][i][BB][1][0];
					lhs_buf[k][j][i][BB][2][3]= lhs_buf[k][j][i][BB][2][3] - coeff*lhs_buf[k][j][i][BB][2][0];
					lhs_buf[k][j][i][BB][3][3]= lhs_buf[k][j][i][BB][3][3] - coeff*lhs_buf[k][j][i][BB][3][0];
					lhs_buf[k][j][i][BB][4][3]= lhs_buf[k][j][i][BB][4][3] - coeff*lhs_buf[k][j][i][BB][4][0];
					lhs_buf[k][j][i][CC][0][3] = lhs_buf[k][j][i][CC][0][3] - coeff*lhs_buf[k][j][i][CC][0][0];
					lhs_buf[k][j][i][CC][1][3] = lhs_buf[k][j][i][CC][1][3] - coeff*lhs_buf[k][j][i][CC][1][0];
					lhs_buf[k][j][i][CC][2][3] = lhs_buf[k][j][i][CC][2][3] - coeff*lhs_buf[k][j][i][CC][2][0];
					lhs_buf[k][j][i][CC][3][3] = lhs_buf[k][j][i][CC][3][3] - coeff*lhs_buf[k][j][i][CC][3][0];
					lhs_buf[k][j][i][CC][4][3] = lhs_buf[k][j][i][CC][4][3] - coeff*lhs_buf[k][j][i][CC][4][0];
					rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff*rhs[k][j][i][0];

					coeff = lhs_buf[k][j][i][BB][0][4];
					lhs_buf[k][j][i][BB][1][4]= lhs_buf[k][j][i][BB][1][4] - coeff*lhs_buf[k][j][i][BB][1][0];
					lhs_buf[k][j][i][BB][2][4]= lhs_buf[k][j][i][BB][2][4] - coeff*lhs_buf[k][j][i][BB][2][0];
					lhs_buf[k][j][i][BB][3][4]= lhs_buf[k][j][i][BB][3][4] - coeff*lhs_buf[k][j][i][BB][3][0];
					lhs_buf[k][j][i][BB][4][4]= lhs_buf[k][j][i][BB][4][4] - coeff*lhs_buf[k][j][i][BB][4][0];
					lhs_buf[k][j][i][CC][0][4] = lhs_buf[k][j][i][CC][0][4] - coeff*lhs_buf[k][j][i][CC][0][0];
					lhs_buf[k][j][i][CC][1][4] = lhs_buf[k][j][i][CC][1][4] - coeff*lhs_buf[k][j][i][CC][1][0];
					lhs_buf[k][j][i][CC][2][4] = lhs_buf[k][j][i][CC][2][4] - coeff*lhs_buf[k][j][i][CC][2][0];
					lhs_buf[k][j][i][CC][3][4] = lhs_buf[k][j][i][CC][3][4] - coeff*lhs_buf[k][j][i][CC][3][0];
					lhs_buf[k][j][i][CC][4][4] = lhs_buf[k][j][i][CC][4][4] - coeff*lhs_buf[k][j][i][CC][4][0];
					rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff*rhs[k][j][i][0];


					pivot = 1.00/lhs_buf[k][j][i][BB][1][1];
					lhs_buf[k][j][i][BB][2][1] = lhs_buf[k][j][i][BB][2][1]*pivot;
					lhs_buf[k][j][i][BB][3][1] = lhs_buf[k][j][i][BB][3][1]*pivot;
					lhs_buf[k][j][i][BB][4][1] = lhs_buf[k][j][i][BB][4][1]*pivot;
					lhs_buf[k][j][i][CC][0][1] = lhs_buf[k][j][i][CC][0][1]*pivot;
					lhs_buf[k][j][i][CC][1][1] = lhs_buf[k][j][i][CC][1][1]*pivot;
					lhs_buf[k][j][i][CC][2][1] = lhs_buf[k][j][i][CC][2][1]*pivot;
					lhs_buf[k][j][i][CC][3][1] = lhs_buf[k][j][i][CC][3][1]*pivot;
					lhs_buf[k][j][i][CC][4][1] = lhs_buf[k][j][i][CC][4][1]*pivot;
					rhs[k][j][i][1]   = rhs[k][j][i][1]  *pivot;

					coeff = lhs_buf[k][j][i][BB][1][0];
					lhs_buf[k][j][i][BB][2][0]= lhs_buf[k][j][i][BB][2][0] - coeff*lhs_buf[k][j][i][BB][2][1];
					lhs_buf[k][j][i][BB][3][0]= lhs_buf[k][j][i][BB][3][0] - coeff*lhs_buf[k][j][i][BB][3][1];
					lhs_buf[k][j][i][BB][4][0]= lhs_buf[k][j][i][BB][4][0] - coeff*lhs_buf[k][j][i][BB][4][1];
					lhs_buf[k][j][i][CC][0][0] = lhs_buf[k][j][i][CC][0][0] - coeff*lhs_buf[k][j][i][CC][0][1];
					lhs_buf[k][j][i][CC][1][0] = lhs_buf[k][j][i][CC][1][0] - coeff*lhs_buf[k][j][i][CC][1][1];
					lhs_buf[k][j][i][CC][2][0] = lhs_buf[k][j][i][CC][2][0] - coeff*lhs_buf[k][j][i][CC][2][1];
					lhs_buf[k][j][i][CC][3][0] = lhs_buf[k][j][i][CC][3][0] - coeff*lhs_buf[k][j][i][CC][3][1];
					lhs_buf[k][j][i][CC][4][0] = lhs_buf[k][j][i][CC][4][0] - coeff*lhs_buf[k][j][i][CC][4][1];
					rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff*rhs[k][j][i][1];

					coeff = lhs_buf[k][j][i][BB][1][2];
					lhs_buf[k][j][i][BB][2][2]= lhs_buf[k][j][i][BB][2][2] - coeff*lhs_buf[k][j][i][BB][2][1];
					lhs_buf[k][j][i][BB][3][2]= lhs_buf[k][j][i][BB][3][2] - coeff*lhs_buf[k][j][i][BB][3][1];
					lhs_buf[k][j][i][BB][4][2]= lhs_buf[k][j][i][BB][4][2] - coeff*lhs_buf[k][j][i][BB][4][1];
					lhs_buf[k][j][i][CC][0][2] = lhs_buf[k][j][i][CC][0][2] - coeff*lhs_buf[k][j][i][CC][0][1];
					lhs_buf[k][j][i][CC][1][2] = lhs_buf[k][j][i][CC][1][2] - coeff*lhs_buf[k][j][i][CC][1][1];
					lhs_buf[k][j][i][CC][2][2] = lhs_buf[k][j][i][CC][2][2] - coeff*lhs_buf[k][j][i][CC][2][1];
					lhs_buf[k][j][i][CC][3][2] = lhs_buf[k][j][i][CC][3][2] - coeff*lhs_buf[k][j][i][CC][3][1];
					lhs_buf[k][j][i][CC][4][2] = lhs_buf[k][j][i][CC][4][2] - coeff*lhs_buf[k][j][i][CC][4][1];
					rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff*rhs[k][j][i][1];

					coeff = lhs_buf[k][j][i][BB][1][3];
					lhs_buf[k][j][i][BB][2][3]= lhs_buf[k][j][i][BB][2][3] - coeff*lhs_buf[k][j][i][BB][2][1];
					lhs_buf[k][j][i][BB][3][3]= lhs_buf[k][j][i][BB][3][3] - coeff*lhs_buf[k][j][i][BB][3][1];
					lhs_buf[k][j][i][BB][4][3]= lhs_buf[k][j][i][BB][4][3] - coeff*lhs_buf[k][j][i][BB][4][1];
					lhs_buf[k][j][i][CC][0][3] = lhs_buf[k][j][i][CC][0][3] - coeff*lhs_buf[k][j][i][CC][0][1];
					lhs_buf[k][j][i][CC][1][3] = lhs_buf[k][j][i][CC][1][3] - coeff*lhs_buf[k][j][i][CC][1][1];
					lhs_buf[k][j][i][CC][2][3] = lhs_buf[k][j][i][CC][2][3] - coeff*lhs_buf[k][j][i][CC][2][1];
					lhs_buf[k][j][i][CC][3][3] = lhs_buf[k][j][i][CC][3][3] - coeff*lhs_buf[k][j][i][CC][3][1];
					lhs_buf[k][j][i][CC][4][3] = lhs_buf[k][j][i][CC][4][3] - coeff*lhs_buf[k][j][i][CC][4][1];
					rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff*rhs[k][j][i][1];

					coeff = lhs_buf[k][j][i][BB][1][4];
					lhs_buf[k][j][i][BB][2][4]= lhs_buf[k][j][i][BB][2][4] - coeff*lhs_buf[k][j][i][BB][2][1];
					lhs_buf[k][j][i][BB][3][4]= lhs_buf[k][j][i][BB][3][4] - coeff*lhs_buf[k][j][i][BB][3][1];
					lhs_buf[k][j][i][BB][4][4]= lhs_buf[k][j][i][BB][4][4] - coeff*lhs_buf[k][j][i][BB][4][1];
					lhs_buf[k][j][i][CC][0][4] = lhs_buf[k][j][i][CC][0][4] - coeff*lhs_buf[k][j][i][CC][0][1];
					lhs_buf[k][j][i][CC][1][4] = lhs_buf[k][j][i][CC][1][4] - coeff*lhs_buf[k][j][i][CC][1][1];
					lhs_buf[k][j][i][CC][2][4] = lhs_buf[k][j][i][CC][2][4] - coeff*lhs_buf[k][j][i][CC][2][1];
					lhs_buf[k][j][i][CC][3][4] = lhs_buf[k][j][i][CC][3][4] - coeff*lhs_buf[k][j][i][CC][3][1];
					lhs_buf[k][j][i][CC][4][4] = lhs_buf[k][j][i][CC][4][4] - coeff*lhs_buf[k][j][i][CC][4][1];
					rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff*rhs[k][j][i][1];


					pivot = 1.00/lhs_buf[k][j][i][BB][2][2];
					lhs_buf[k][j][i][BB][3][2] = lhs_buf[k][j][i][BB][3][2]*pivot;
					lhs_buf[k][j][i][BB][4][2] = lhs_buf[k][j][i][BB][4][2]*pivot;
					lhs_buf[k][j][i][CC][0][2] = lhs_buf[k][j][i][CC][0][2]*pivot;
					lhs_buf[k][j][i][CC][1][2] = lhs_buf[k][j][i][CC][1][2]*pivot;
					lhs_buf[k][j][i][CC][2][2] = lhs_buf[k][j][i][CC][2][2]*pivot;
					lhs_buf[k][j][i][CC][3][2] = lhs_buf[k][j][i][CC][3][2]*pivot;
					lhs_buf[k][j][i][CC][4][2] = lhs_buf[k][j][i][CC][4][2]*pivot;
					rhs[k][j][i][2]   = rhs[k][j][i][2]  *pivot;

					coeff = lhs_buf[k][j][i][BB][2][0];
					lhs_buf[k][j][i][BB][3][0]= lhs_buf[k][j][i][BB][3][0] - coeff*lhs_buf[k][j][i][BB][3][2];
					lhs_buf[k][j][i][BB][4][0]= lhs_buf[k][j][i][BB][4][0] - coeff*lhs_buf[k][j][i][BB][4][2];
					lhs_buf[k][j][i][CC][0][0] = lhs_buf[k][j][i][CC][0][0] - coeff*lhs_buf[k][j][i][CC][0][2];
					lhs_buf[k][j][i][CC][1][0] = lhs_buf[k][j][i][CC][1][0] - coeff*lhs_buf[k][j][i][CC][1][2];
					lhs_buf[k][j][i][CC][2][0] = lhs_buf[k][j][i][CC][2][0] - coeff*lhs_buf[k][j][i][CC][2][2];
					lhs_buf[k][j][i][CC][3][0] = lhs_buf[k][j][i][CC][3][0] - coeff*lhs_buf[k][j][i][CC][3][2];
					lhs_buf[k][j][i][CC][4][0] = lhs_buf[k][j][i][CC][4][0] - coeff*lhs_buf[k][j][i][CC][4][2];
					rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff*rhs[k][j][i][2];

					coeff = lhs_buf[k][j][i][BB][2][1];
					lhs_buf[k][j][i][BB][3][1]= lhs_buf[k][j][i][BB][3][1] - coeff*lhs_buf[k][j][i][BB][3][2];
					lhs_buf[k][j][i][BB][4][1]= lhs_buf[k][j][i][BB][4][1] - coeff*lhs_buf[k][j][i][BB][4][2];
					lhs_buf[k][j][i][CC][0][1] = lhs_buf[k][j][i][CC][0][1] - coeff*lhs_buf[k][j][i][CC][0][2];
					lhs_buf[k][j][i][CC][1][1] = lhs_buf[k][j][i][CC][1][1] - coeff*lhs_buf[k][j][i][CC][1][2];
					lhs_buf[k][j][i][CC][2][1] = lhs_buf[k][j][i][CC][2][1] - coeff*lhs_buf[k][j][i][CC][2][2];
					lhs_buf[k][j][i][CC][3][1] = lhs_buf[k][j][i][CC][3][1] - coeff*lhs_buf[k][j][i][CC][3][2];
					lhs_buf[k][j][i][CC][4][1] = lhs_buf[k][j][i][CC][4][1] - coeff*lhs_buf[k][j][i][CC][4][2];
					rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff*rhs[k][j][i][2];

					coeff = lhs_buf[k][j][i][BB][2][3];
					lhs_buf[k][j][i][BB][3][3]= lhs_buf[k][j][i][BB][3][3] - coeff*lhs_buf[k][j][i][BB][3][2];
					lhs_buf[k][j][i][BB][4][3]= lhs_buf[k][j][i][BB][4][3] - coeff*lhs_buf[k][j][i][BB][4][2];
					lhs_buf[k][j][i][CC][0][3] = lhs_buf[k][j][i][CC][0][3] - coeff*lhs_buf[k][j][i][CC][0][2];
					lhs_buf[k][j][i][CC][1][3] = lhs_buf[k][j][i][CC][1][3] - coeff*lhs_buf[k][j][i][CC][1][2];
					lhs_buf[k][j][i][CC][2][3] = lhs_buf[k][j][i][CC][2][3] - coeff*lhs_buf[k][j][i][CC][2][2];
					lhs_buf[k][j][i][CC][3][3] = lhs_buf[k][j][i][CC][3][3] - coeff*lhs_buf[k][j][i][CC][3][2];
					lhs_buf[k][j][i][CC][4][3] = lhs_buf[k][j][i][CC][4][3] - coeff*lhs_buf[k][j][i][CC][4][2];
					rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff*rhs[k][j][i][2];

					coeff = lhs_buf[k][j][i][BB][2][4];
					lhs_buf[k][j][i][BB][3][4]= lhs_buf[k][j][i][BB][3][4] - coeff*lhs_buf[k][j][i][BB][3][2];
					lhs_buf[k][j][i][BB][4][4]= lhs_buf[k][j][i][BB][4][4] - coeff*lhs_buf[k][j][i][BB][4][2];
					lhs_buf[k][j][i][CC][0][4] = lhs_buf[k][j][i][CC][0][4] - coeff*lhs_buf[k][j][i][CC][0][2];
					lhs_buf[k][j][i][CC][1][4] = lhs_buf[k][j][i][CC][1][4] - coeff*lhs_buf[k][j][i][CC][1][2];
					lhs_buf[k][j][i][CC][2][4] = lhs_buf[k][j][i][CC][2][4] - coeff*lhs_buf[k][j][i][CC][2][2];
					lhs_buf[k][j][i][CC][3][4] = lhs_buf[k][j][i][CC][3][4] - coeff*lhs_buf[k][j][i][CC][3][2];
					lhs_buf[k][j][i][CC][4][4] = lhs_buf[k][j][i][CC][4][4] - coeff*lhs_buf[k][j][i][CC][4][2];
					rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff*rhs[k][j][i][2];


					pivot = 1.00/lhs_buf[k][j][i][BB][3][3];
					lhs_buf[k][j][i][BB][4][3] = lhs_buf[k][j][i][BB][4][3]*pivot;
					lhs_buf[k][j][i][CC][0][3] = lhs_buf[k][j][i][CC][0][3]*pivot;
					lhs_buf[k][j][i][CC][1][3] = lhs_buf[k][j][i][CC][1][3]*pivot;
					lhs_buf[k][j][i][CC][2][3] = lhs_buf[k][j][i][CC][2][3]*pivot;
					lhs_buf[k][j][i][CC][3][3] = lhs_buf[k][j][i][CC][3][3]*pivot;
					lhs_buf[k][j][i][CC][4][3] = lhs_buf[k][j][i][CC][4][3]*pivot;
					rhs[k][j][i][3]   = rhs[k][j][i][3]  *pivot;

					coeff = lhs_buf[k][j][i][BB][3][0];
					lhs_buf[k][j][i][BB][4][0]= lhs_buf[k][j][i][BB][4][0] - coeff*lhs_buf[k][j][i][BB][4][3];
					lhs_buf[k][j][i][CC][0][0] = lhs_buf[k][j][i][CC][0][0] - coeff*lhs_buf[k][j][i][CC][0][3];
					lhs_buf[k][j][i][CC][1][0] = lhs_buf[k][j][i][CC][1][0] - coeff*lhs_buf[k][j][i][CC][1][3];
					lhs_buf[k][j][i][CC][2][0] = lhs_buf[k][j][i][CC][2][0] - coeff*lhs_buf[k][j][i][CC][2][3];
					lhs_buf[k][j][i][CC][3][0] = lhs_buf[k][j][i][CC][3][0] - coeff*lhs_buf[k][j][i][CC][3][3];
					lhs_buf[k][j][i][CC][4][0] = lhs_buf[k][j][i][CC][4][0] - coeff*lhs_buf[k][j][i][CC][4][3];
					rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff*rhs[k][j][i][3];

					coeff = lhs_buf[k][j][i][BB][3][1];
					lhs_buf[k][j][i][BB][4][1]= lhs_buf[k][j][i][BB][4][1] - coeff*lhs_buf[k][j][i][BB][4][3];
					lhs_buf[k][j][i][CC][0][1] = lhs_buf[k][j][i][CC][0][1] - coeff*lhs_buf[k][j][i][CC][0][3];
					lhs_buf[k][j][i][CC][1][1] = lhs_buf[k][j][i][CC][1][1] - coeff*lhs_buf[k][j][i][CC][1][3];
					lhs_buf[k][j][i][CC][2][1] = lhs_buf[k][j][i][CC][2][1] - coeff*lhs_buf[k][j][i][CC][2][3];
					lhs_buf[k][j][i][CC][3][1] = lhs_buf[k][j][i][CC][3][1] - coeff*lhs_buf[k][j][i][CC][3][3];
					lhs_buf[k][j][i][CC][4][1] = lhs_buf[k][j][i][CC][4][1] - coeff*lhs_buf[k][j][i][CC][4][3];
					rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff*rhs[k][j][i][3];

					coeff = lhs_buf[k][j][i][BB][3][2];
					lhs_buf[k][j][i][BB][4][2]= lhs_buf[k][j][i][BB][4][2] - coeff*lhs_buf[k][j][i][BB][4][3];
					lhs_buf[k][j][i][CC][0][2] = lhs_buf[k][j][i][CC][0][2] - coeff*lhs_buf[k][j][i][CC][0][3];
					lhs_buf[k][j][i][CC][1][2] = lhs_buf[k][j][i][CC][1][2] - coeff*lhs_buf[k][j][i][CC][1][3];
					lhs_buf[k][j][i][CC][2][2] = lhs_buf[k][j][i][CC][2][2] - coeff*lhs_buf[k][j][i][CC][2][3];
					lhs_buf[k][j][i][CC][3][2] = lhs_buf[k][j][i][CC][3][2] - coeff*lhs_buf[k][j][i][CC][3][3];
					lhs_buf[k][j][i][CC][4][2] = lhs_buf[k][j][i][CC][4][2] - coeff*lhs_buf[k][j][i][CC][4][3];
					rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff*rhs[k][j][i][3];

					coeff = lhs_buf[k][j][i][BB][3][4];
					lhs_buf[k][j][i][BB][4][4]= lhs_buf[k][j][i][BB][4][4] - coeff*lhs_buf[k][j][i][BB][4][3];
					lhs_buf[k][j][i][CC][0][4] = lhs_buf[k][j][i][CC][0][4] - coeff*lhs_buf[k][j][i][CC][0][3];
					lhs_buf[k][j][i][CC][1][4] = lhs_buf[k][j][i][CC][1][4] - coeff*lhs_buf[k][j][i][CC][1][3];
					lhs_buf[k][j][i][CC][2][4] = lhs_buf[k][j][i][CC][2][4] - coeff*lhs_buf[k][j][i][CC][2][3];
					lhs_buf[k][j][i][CC][3][4] = lhs_buf[k][j][i][CC][3][4] - coeff*lhs_buf[k][j][i][CC][3][3];
					lhs_buf[k][j][i][CC][4][4] = lhs_buf[k][j][i][CC][4][4] - coeff*lhs_buf[k][j][i][CC][4][3];
					rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff*rhs[k][j][i][3];


					pivot = 1.00/lhs_buf[k][j][i][BB][4][4];
					lhs_buf[k][j][i][CC][0][4] = lhs_buf[k][j][i][CC][0][4]*pivot;
					lhs_buf[k][j][i][CC][1][4] = lhs_buf[k][j][i][CC][1][4]*pivot;
					lhs_buf[k][j][i][CC][2][4] = lhs_buf[k][j][i][CC][2][4]*pivot;
					lhs_buf[k][j][i][CC][3][4] = lhs_buf[k][j][i][CC][3][4]*pivot;
					lhs_buf[k][j][i][CC][4][4] = lhs_buf[k][j][i][CC][4][4]*pivot;
					rhs[k][j][i][4]   = rhs[k][j][i][4]  *pivot;

					coeff = lhs_buf[k][j][i][BB][4][0];
					lhs_buf[k][j][i][CC][0][0] = lhs_buf[k][j][i][CC][0][0] - coeff*lhs_buf[k][j][i][CC][0][4];
					lhs_buf[k][j][i][CC][1][0] = lhs_buf[k][j][i][CC][1][0] - coeff*lhs_buf[k][j][i][CC][1][4];
					lhs_buf[k][j][i][CC][2][0] = lhs_buf[k][j][i][CC][2][0] - coeff*lhs_buf[k][j][i][CC][2][4];
					lhs_buf[k][j][i][CC][3][0] = lhs_buf[k][j][i][CC][3][0] - coeff*lhs_buf[k][j][i][CC][3][4];
					lhs_buf[k][j][i][CC][4][0] = lhs_buf[k][j][i][CC][4][0] - coeff*lhs_buf[k][j][i][CC][4][4];
					rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff*rhs[k][j][i][4];

					coeff = lhs_buf[k][j][i][BB][4][1];
					lhs_buf[k][j][i][CC][0][1] = lhs_buf[k][j][i][CC][0][1] - coeff*lhs_buf[k][j][i][CC][0][4];
					lhs_buf[k][j][i][CC][1][1] = lhs_buf[k][j][i][CC][1][1] - coeff*lhs_buf[k][j][i][CC][1][4];
					lhs_buf[k][j][i][CC][2][1] = lhs_buf[k][j][i][CC][2][1] - coeff*lhs_buf[k][j][i][CC][2][4];
					lhs_buf[k][j][i][CC][3][1] = lhs_buf[k][j][i][CC][3][1] - coeff*lhs_buf[k][j][i][CC][3][4];
					lhs_buf[k][j][i][CC][4][1] = lhs_buf[k][j][i][CC][4][1] - coeff*lhs_buf[k][j][i][CC][4][4];
					rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff*rhs[k][j][i][4];

					coeff = lhs_buf[k][j][i][BB][4][2];
					lhs_buf[k][j][i][CC][0][2] = lhs_buf[k][j][i][CC][0][2] - coeff*lhs_buf[k][j][i][CC][0][4];
					lhs_buf[k][j][i][CC][1][2] = lhs_buf[k][j][i][CC][1][2] - coeff*lhs_buf[k][j][i][CC][1][4];
					lhs_buf[k][j][i][CC][2][2] = lhs_buf[k][j][i][CC][2][2] - coeff*lhs_buf[k][j][i][CC][2][4];
					lhs_buf[k][j][i][CC][3][2] = lhs_buf[k][j][i][CC][3][2] - coeff*lhs_buf[k][j][i][CC][3][4];
					lhs_buf[k][j][i][CC][4][2] = lhs_buf[k][j][i][CC][4][2] - coeff*lhs_buf[k][j][i][CC][4][4];
					rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff*rhs[k][j][i][4];

					coeff = lhs_buf[k][j][i][BB][4][3];
					lhs_buf[k][j][i][CC][0][3] = lhs_buf[k][j][i][CC][0][3] - coeff*lhs_buf[k][j][i][CC][0][4];
					lhs_buf[k][j][i][CC][1][3] = lhs_buf[k][j][i][CC][1][3] - coeff*lhs_buf[k][j][i][CC][1][4];
					lhs_buf[k][j][i][CC][2][3] = lhs_buf[k][j][i][CC][2][3] - coeff*lhs_buf[k][j][i][CC][2][4];
					lhs_buf[k][j][i][CC][3][3] = lhs_buf[k][j][i][CC][3][3] - coeff*lhs_buf[k][j][i][CC][3][4];
					lhs_buf[k][j][i][CC][4][3] = lhs_buf[k][j][i][CC][4][3] - coeff*lhs_buf[k][j][i][CC][4][4];
					rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff*rhs[k][j][i][4];
				}
				else if(j == jsize)
				{					
					rhs[k][j][i][0] = rhs[k][j][i][0] - lhs_buf[k][j][i][AA][0][0]*rhs[k][j - 1][i][0]
					- lhs_buf[k][j][i][AA][1][0]*rhs[k][j - 1][i][1]
					- lhs_buf[k][j][i][AA][2][0]*rhs[k][j - 1][i][2]
					- lhs_buf[k][j][i][AA][3][0]*rhs[k][j - 1][i][3]
					- lhs_buf[k][j][i][AA][4][0]*rhs[k][j - 1][i][4];
					rhs[k][j][i][1] = rhs[k][j][i][1] - lhs_buf[k][j][i][AA][0][1]*rhs[k][j - 1][i][0]
					- lhs_buf[k][j][i][AA][1][1]*rhs[k][j - 1][i][1]
					- lhs_buf[k][j][i][AA][2][1]*rhs[k][j - 1][i][2]
					- lhs_buf[k][j][i][AA][3][1]*rhs[k][j - 1][i][3]
					- lhs_buf[k][j][i][AA][4][1]*rhs[k][j - 1][i][4];
					rhs[k][j][i][2] = rhs[k][j][i][2] - lhs_buf[k][j][i][AA][0][2]*rhs[k][j - 1][i][0]
					- lhs_buf[k][j][i][AA][1][2]*rhs[k][j - 1][i][1]
					- lhs_buf[k][j][i][AA][2][2]*rhs[k][j - 1][i][2]
					- lhs_buf[k][j][i][AA][3][2]*rhs[k][j - 1][i][3]
					- lhs_buf[k][j][i][AA][4][2]*rhs[k][j - 1][i][4];
					rhs[k][j][i][3] = rhs[k][j][i][3] - lhs_buf[k][j][i][AA][0][3]*rhs[k][j - 1][i][0]
					- lhs_buf[k][j][i][AA][1][3]*rhs[k][j - 1][i][1]
					- lhs_buf[k][j][i][AA][2][3]*rhs[k][j - 1][i][2]
					- lhs_buf[k][j][i][AA][3][3]*rhs[k][j - 1][i][3]
					- lhs_buf[k][j][i][AA][4][3]*rhs[k][j - 1][i][4];
					rhs[k][j][i][4] = rhs[k][j][i][4] - lhs_buf[k][j][i][AA][0][4]*rhs[k][j - 1][i][0]
					- lhs_buf[k][j][i][AA][1][4]*rhs[k][j - 1][i][1]
					- lhs_buf[k][j][i][AA][2][4]*rhs[k][j - 1][i][2]
					- lhs_buf[k][j][i][AA][3][4]*rhs[k][j - 1][i][3]
					- lhs_buf[k][j][i][AA][4][4]*rhs[k][j - 1][i][4];



					
					lhs_buf[k][j][i][BB][0][0] = lhs_buf[k][j][i][BB][0][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j - 1][i][CC][0][0]
						  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j - 1][i][CC][0][1]
						  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j - 1][i][CC][0][2]
						  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j - 1][i][CC][0][3]
						  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j - 1][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][1] = lhs_buf[k][j][i][BB][0][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j - 1][i][CC][0][0]
						  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j - 1][i][CC][0][1]
						  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j - 1][i][CC][0][2]
						  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j - 1][i][CC][0][3]
						  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j - 1][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][2] = lhs_buf[k][j][i][BB][0][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j - 1][i][CC][0][0]
						  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j - 1][i][CC][0][1]
						  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j - 1][i][CC][0][2]
						  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j - 1][i][CC][0][3]
						  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j - 1][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][3] = lhs_buf[k][j][i][BB][0][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j - 1][i][CC][0][0]
						  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j - 1][i][CC][0][1]
						  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j - 1][i][CC][0][2]
						  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j - 1][i][CC][0][3]
						  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j - 1][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][4] = lhs_buf[k][j][i][BB][0][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j - 1][i][CC][0][0]
						  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j - 1][i][CC][0][1]
						  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j - 1][i][CC][0][2]
						  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j - 1][i][CC][0][3]
						  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j - 1][i][CC][0][4];
					lhs_buf[k][j][i][BB][1][0] = lhs_buf[k][j][i][BB][1][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j - 1][i][CC][1][0]
						  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j - 1][i][CC][1][1]
						  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j - 1][i][CC][1][2]
						  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j - 1][i][CC][1][3]
						  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j - 1][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][1] = lhs_buf[k][j][i][BB][1][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j - 1][i][CC][1][0]
						  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j - 1][i][CC][1][1]
						  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j - 1][i][CC][1][2]
						  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j - 1][i][CC][1][3]
						  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j - 1][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][2] = lhs_buf[k][j][i][BB][1][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j - 1][i][CC][1][0]
						  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j - 1][i][CC][1][1]
						  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j - 1][i][CC][1][2]
						  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j - 1][i][CC][1][3]
						  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j - 1][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][3] = lhs_buf[k][j][i][BB][1][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j - 1][i][CC][1][0]
						  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j - 1][i][CC][1][1]
						  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j - 1][i][CC][1][2]
						  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j - 1][i][CC][1][3]
						  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j - 1][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][4] = lhs_buf[k][j][i][BB][1][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j - 1][i][CC][1][0]
						  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j - 1][i][CC][1][1]
						  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j - 1][i][CC][1][2]
						  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j - 1][i][CC][1][3]
						  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j - 1][i][CC][1][4];
					lhs_buf[k][j][i][BB][2][0] = lhs_buf[k][j][i][BB][2][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j - 1][i][CC][2][0]
						  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j - 1][i][CC][2][1]
						  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j - 1][i][CC][2][2]
						  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j - 1][i][CC][2][3]
						  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j - 1][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][1] = lhs_buf[k][j][i][BB][2][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j - 1][i][CC][2][0]
						  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j - 1][i][CC][2][1]
						  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j - 1][i][CC][2][2]
						  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j - 1][i][CC][2][3]
						  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j - 1][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][2] = lhs_buf[k][j][i][BB][2][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j - 1][i][CC][2][0]
						  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j - 1][i][CC][2][1]
						  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j - 1][i][CC][2][2]
						  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j - 1][i][CC][2][3]
						  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j - 1][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][3] = lhs_buf[k][j][i][BB][2][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j - 1][i][CC][2][0]
						  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j - 1][i][CC][2][1]
						  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j - 1][i][CC][2][2]
						  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j - 1][i][CC][2][3]
						  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j - 1][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][4] = lhs_buf[k][j][i][BB][2][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j - 1][i][CC][2][0]
						  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j - 1][i][CC][2][1]
						  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j - 1][i][CC][2][2]
						  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j - 1][i][CC][2][3]
						  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j - 1][i][CC][2][4];
					lhs_buf[k][j][i][BB][3][0] = lhs_buf[k][j][i][BB][3][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j - 1][i][CC][3][0]
						  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j - 1][i][CC][3][1]
						  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j - 1][i][CC][3][2]
						  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j - 1][i][CC][3][3]
						  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j - 1][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][1] = lhs_buf[k][j][i][BB][3][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j - 1][i][CC][3][0]
						  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j - 1][i][CC][3][1]
						  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j - 1][i][CC][3][2]
						  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j - 1][i][CC][3][3]
						  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j - 1][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][2] = lhs_buf[k][j][i][BB][3][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j - 1][i][CC][3][0]
						  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j - 1][i][CC][3][1]
						  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j - 1][i][CC][3][2]
						  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j - 1][i][CC][3][3]
						  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j - 1][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][3] = lhs_buf[k][j][i][BB][3][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j - 1][i][CC][3][0]
						  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j - 1][i][CC][3][1]
						  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j - 1][i][CC][3][2]
						  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j - 1][i][CC][3][3]
						  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j - 1][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][4] = lhs_buf[k][j][i][BB][3][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j - 1][i][CC][3][0]
						  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j - 1][i][CC][3][1]
						  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j - 1][i][CC][3][2]
						  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j - 1][i][CC][3][3]
						  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j - 1][i][CC][3][4];
					lhs_buf[k][j][i][BB][4][0] = lhs_buf[k][j][i][BB][4][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j - 1][i][CC][4][0]
						  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j - 1][i][CC][4][1]
						  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j - 1][i][CC][4][2]
						  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j - 1][i][CC][4][3]
						  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j - 1][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][1] = lhs_buf[k][j][i][BB][4][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j - 1][i][CC][4][0]
						  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j - 1][i][CC][4][1]
						  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j - 1][i][CC][4][2]
						  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j - 1][i][CC][4][3]
						  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j - 1][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][2] = lhs_buf[k][j][i][BB][4][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j - 1][i][CC][4][0]
						  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j - 1][i][CC][4][1]
						  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j - 1][i][CC][4][2]
						  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j - 1][i][CC][4][3]
						  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j - 1][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][3] = lhs_buf[k][j][i][BB][4][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j - 1][i][CC][4][0]
						  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j - 1][i][CC][4][1]
						  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j - 1][i][CC][4][2]
						  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j - 1][i][CC][4][3]
						  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j - 1][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][4] = lhs_buf[k][j][i][BB][4][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j - 1][i][CC][4][0]
						  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j - 1][i][CC][4][1]
						  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j - 1][i][CC][4][2]
						  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j - 1][i][CC][4][3]
						  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j - 1][i][CC][4][4];
					
					

					//binvrhs( lhs[jsize][BB], rhs[k][jsize][i] );
					//binvrhs(jsize, BB, k, jsize, i);
					pivot = 1.00/lhs_buf[k][j][i][BB][0][0];
					lhs_buf[k][j][i][BB][1][0] = lhs_buf[k][j][i][BB][1][0]*pivot;
					lhs_buf[k][j][i][BB][2][0] = lhs_buf[k][j][i][BB][2][0]*pivot;
					lhs_buf[k][j][i][BB][3][0] = lhs_buf[k][j][i][BB][3][0]*pivot;
					lhs_buf[k][j][i][BB][4][0] = lhs_buf[k][j][i][BB][4][0]*pivot;
					rhs[k][j][i][0]   = rhs[k][j][i][0]  *pivot;

					coeff = lhs_buf[k][j][i][BB][0][1];
					lhs_buf[k][j][i][BB][1][1]= lhs_buf[k][j][i][BB][1][1] - coeff*lhs_buf[k][j][i][BB][1][0];
					lhs_buf[k][j][i][BB][2][1]= lhs_buf[k][j][i][BB][2][1] - coeff*lhs_buf[k][j][i][BB][2][0];
					lhs_buf[k][j][i][BB][3][1]= lhs_buf[k][j][i][BB][3][1] - coeff*lhs_buf[k][j][i][BB][3][0];
					lhs_buf[k][j][i][BB][4][1]= lhs_buf[k][j][i][BB][4][1] - coeff*lhs_buf[k][j][i][BB][4][0];
					rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff*rhs[k][j][i][0];

					coeff = lhs_buf[k][j][i][BB][0][2];
					lhs_buf[k][j][i][BB][1][2]= lhs_buf[k][j][i][BB][1][2] - coeff*lhs_buf[k][j][i][BB][1][0];
					lhs_buf[k][j][i][BB][2][2]= lhs_buf[k][j][i][BB][2][2] - coeff*lhs_buf[k][j][i][BB][2][0];
					lhs_buf[k][j][i][BB][3][2]= lhs_buf[k][j][i][BB][3][2] - coeff*lhs_buf[k][j][i][BB][3][0];
					lhs_buf[k][j][i][BB][4][2]= lhs_buf[k][j][i][BB][4][2] - coeff*lhs_buf[k][j][i][BB][4][0];
					rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff*rhs[k][j][i][0];

					coeff = lhs_buf[k][j][i][BB][0][3];
					lhs_buf[k][j][i][BB][1][3]= lhs_buf[k][j][i][BB][1][3] - coeff*lhs_buf[k][j][i][BB][1][0];
					lhs_buf[k][j][i][BB][2][3]= lhs_buf[k][j][i][BB][2][3] - coeff*lhs_buf[k][j][i][BB][2][0];
					lhs_buf[k][j][i][BB][3][3]= lhs_buf[k][j][i][BB][3][3] - coeff*lhs_buf[k][j][i][BB][3][0];
					lhs_buf[k][j][i][BB][4][3]= lhs_buf[k][j][i][BB][4][3] - coeff*lhs_buf[k][j][i][BB][4][0];
					rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff*rhs[k][j][i][0];

					coeff = lhs_buf[k][j][i][BB][0][4];
					lhs_buf[k][j][i][BB][1][4]= lhs_buf[k][j][i][BB][1][4] - coeff*lhs_buf[k][j][i][BB][1][0];
					lhs_buf[k][j][i][BB][2][4]= lhs_buf[k][j][i][BB][2][4] - coeff*lhs_buf[k][j][i][BB][2][0];
					lhs_buf[k][j][i][BB][3][4]= lhs_buf[k][j][i][BB][3][4] - coeff*lhs_buf[k][j][i][BB][3][0];
					lhs_buf[k][j][i][BB][4][4]= lhs_buf[k][j][i][BB][4][4] - coeff*lhs_buf[k][j][i][BB][4][0];
					rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff*rhs[k][j][i][0];


					pivot = 1.00/lhs_buf[k][j][i][BB][1][1];
					lhs_buf[k][j][i][BB][2][1] = lhs_buf[k][j][i][BB][2][1]*pivot;
					lhs_buf[k][j][i][BB][3][1] = lhs_buf[k][j][i][BB][3][1]*pivot;
					lhs_buf[k][j][i][BB][4][1] = lhs_buf[k][j][i][BB][4][1]*pivot;
					rhs[k][j][i][1]   = rhs[k][j][i][1]  *pivot;

					coeff = lhs_buf[k][j][i][BB][1][0];
					lhs_buf[k][j][i][BB][2][0]= lhs_buf[k][j][i][BB][2][0] - coeff*lhs_buf[k][j][i][BB][2][1];
					lhs_buf[k][j][i][BB][3][0]= lhs_buf[k][j][i][BB][3][0] - coeff*lhs_buf[k][j][i][BB][3][1];
					lhs_buf[k][j][i][BB][4][0]= lhs_buf[k][j][i][BB][4][0] - coeff*lhs_buf[k][j][i][BB][4][1];
					rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff*rhs[k][j][i][1];

					coeff = lhs_buf[k][j][i][BB][1][2];
					lhs_buf[k][j][i][BB][2][2]= lhs_buf[k][j][i][BB][2][2] - coeff*lhs_buf[k][j][i][BB][2][1];
					lhs_buf[k][j][i][BB][3][2]= lhs_buf[k][j][i][BB][3][2] - coeff*lhs_buf[k][j][i][BB][3][1];
					lhs_buf[k][j][i][BB][4][2]= lhs_buf[k][j][i][BB][4][2] - coeff*lhs_buf[k][j][i][BB][4][1];
					rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff*rhs[k][j][i][1];

					coeff = lhs_buf[k][j][i][BB][1][3];
					lhs_buf[k][j][i][BB][2][3]= lhs_buf[k][j][i][BB][2][3] - coeff*lhs_buf[k][j][i][BB][2][1];
					lhs_buf[k][j][i][BB][3][3]= lhs_buf[k][j][i][BB][3][3] - coeff*lhs_buf[k][j][i][BB][3][1];
					lhs_buf[k][j][i][BB][4][3]= lhs_buf[k][j][i][BB][4][3] - coeff*lhs_buf[k][j][i][BB][4][1];
					rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff*rhs[k][j][i][1];

					coeff = lhs_buf[k][j][i][BB][1][4];
					lhs_buf[k][j][i][BB][2][4]= lhs_buf[k][j][i][BB][2][4] - coeff*lhs_buf[k][j][i][BB][2][1];
					lhs_buf[k][j][i][BB][3][4]= lhs_buf[k][j][i][BB][3][4] - coeff*lhs_buf[k][j][i][BB][3][1];
					lhs_buf[k][j][i][BB][4][4]= lhs_buf[k][j][i][BB][4][4] - coeff*lhs_buf[k][j][i][BB][4][1];
					rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff*rhs[k][j][i][1];


					pivot = 1.00/lhs_buf[k][j][i][BB][2][2];
					lhs_buf[k][j][i][BB][3][2] = lhs_buf[k][j][i][BB][3][2]*pivot;
					lhs_buf[k][j][i][BB][4][2] = lhs_buf[k][j][i][BB][4][2]*pivot;
					rhs[k][j][i][2]   = rhs[k][j][i][2]  *pivot;

					coeff = lhs_buf[k][j][i][BB][2][0];
					lhs_buf[k][j][i][BB][3][0]= lhs_buf[k][j][i][BB][3][0] - coeff*lhs_buf[k][j][i][BB][3][2];
					lhs_buf[k][j][i][BB][4][0]= lhs_buf[k][j][i][BB][4][0] - coeff*lhs_buf[k][j][i][BB][4][2];
					rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff*rhs[k][j][i][2];

					coeff = lhs_buf[k][j][i][BB][2][1];
					lhs_buf[k][j][i][BB][3][1]= lhs_buf[k][j][i][BB][3][1] - coeff*lhs_buf[k][j][i][BB][3][2];
					lhs_buf[k][j][i][BB][4][1]= lhs_buf[k][j][i][BB][4][1] - coeff*lhs_buf[k][j][i][BB][4][2];
					rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff*rhs[k][j][i][2];

					coeff = lhs_buf[k][j][i][BB][2][3];
					lhs_buf[k][j][i][BB][3][3]= lhs_buf[k][j][i][BB][3][3] - coeff*lhs_buf[k][j][i][BB][3][2];
					lhs_buf[k][j][i][BB][4][3]= lhs_buf[k][j][i][BB][4][3] - coeff*lhs_buf[k][j][i][BB][4][2];
					rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff*rhs[k][j][i][2];

					coeff = lhs_buf[k][j][i][BB][2][4];
					lhs_buf[k][j][i][BB][3][4]= lhs_buf[k][j][i][BB][3][4] - coeff*lhs_buf[k][j][i][BB][3][2];
					lhs_buf[k][j][i][BB][4][4]= lhs_buf[k][j][i][BB][4][4] - coeff*lhs_buf[k][j][i][BB][4][2];
					rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff*rhs[k][j][i][2];


					pivot = 1.00/lhs_buf[k][j][i][BB][3][3];
					lhs_buf[k][j][i][BB][4][3] = lhs_buf[k][j][i][BB][4][3]*pivot;
					rhs[k][j][i][3]   = rhs[k][j][i][3]  *pivot;

					coeff = lhs_buf[k][j][i][BB][3][0];
					lhs_buf[k][j][i][BB][4][0]= lhs_buf[k][j][i][BB][4][0] - coeff*lhs_buf[k][j][i][BB][4][3];
					rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff*rhs[k][j][i][3];

					coeff = lhs_buf[k][j][i][BB][3][1];
					lhs_buf[k][j][i][BB][4][1]= lhs_buf[k][j][i][BB][4][1] - coeff*lhs_buf[k][j][i][BB][4][3];
					rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff*rhs[k][j][i][3];

					coeff = lhs_buf[k][j][i][BB][3][2];
					lhs_buf[k][j][i][BB][4][2]= lhs_buf[k][j][i][BB][4][2] - coeff*lhs_buf[k][j][i][BB][4][3];
					rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff*rhs[k][j][i][3];

					coeff = lhs_buf[k][j][i][BB][3][4];
					lhs_buf[k][j][i][BB][4][4]= lhs_buf[k][j][i][BB][4][4] - coeff*lhs_buf[k][j][i][BB][4][3];
					rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff*rhs[k][j][i][3];


					pivot = 1.00/lhs_buf[k][j][i][BB][4][4];
					rhs[k][j][i][4]   = rhs[k][j][i][4]  *pivot;

					coeff = lhs_buf[k][j][i][BB][4][0];
					rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff*rhs[k][j][i][4];

					coeff = lhs_buf[k][j][i][BB][4][1];
					rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff*rhs[k][j][i][4];

					coeff = lhs_buf[k][j][i][BB][4][2];
					rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff*rhs[k][j][i][4];

					coeff = lhs_buf[k][j][i][BB][4][3];
					rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff*rhs[k][j][i][4];
				}
				else
				{					
					rhs[k][j][i][0] = rhs[k][j][i][0] - lhs_buf[k][j][i][AA][0][0]*rhs[k][j - 1][i][0]
					- lhs_buf[k][j][i][AA][1][0]*rhs[k][j - 1][i][1]
					- lhs_buf[k][j][i][AA][2][0]*rhs[k][j - 1][i][2]
					- lhs_buf[k][j][i][AA][3][0]*rhs[k][j - 1][i][3]
					- lhs_buf[k][j][i][AA][4][0]*rhs[k][j - 1][i][4];
					rhs[k][j][i][1] = rhs[k][j][i][1] - lhs_buf[k][j][i][AA][0][1]*rhs[k][j - 1][i][0]
					- lhs_buf[k][j][i][AA][1][1]*rhs[k][j - 1][i][1]
					- lhs_buf[k][j][i][AA][2][1]*rhs[k][j - 1][i][2]
					- lhs_buf[k][j][i][AA][3][1]*rhs[k][j - 1][i][3]
					- lhs_buf[k][j][i][AA][4][1]*rhs[k][j - 1][i][4];
					rhs[k][j][i][2] = rhs[k][j][i][2] - lhs_buf[k][j][i][AA][0][2]*rhs[k][j - 1][i][0]
					- lhs_buf[k][j][i][AA][1][2]*rhs[k][j - 1][i][1]
					- lhs_buf[k][j][i][AA][2][2]*rhs[k][j - 1][i][2]
					- lhs_buf[k][j][i][AA][3][2]*rhs[k][j - 1][i][3]
					- lhs_buf[k][j][i][AA][4][2]*rhs[k][j - 1][i][4];
					rhs[k][j][i][3] = rhs[k][j][i][3] - lhs_buf[k][j][i][AA][0][3]*rhs[k][j - 1][i][0]
					- lhs_buf[k][j][i][AA][1][3]*rhs[k][j - 1][i][1]
					- lhs_buf[k][j][i][AA][2][3]*rhs[k][j - 1][i][2]
					- lhs_buf[k][j][i][AA][3][3]*rhs[k][j - 1][i][3]
					- lhs_buf[k][j][i][AA][4][3]*rhs[k][j - 1][i][4];
					rhs[k][j][i][4] = rhs[k][j][i][4] - lhs_buf[k][j][i][AA][0][4]*rhs[k][j - 1][i][0]
					- lhs_buf[k][j][i][AA][1][4]*rhs[k][j - 1][i][1]
					- lhs_buf[k][j][i][AA][2][4]*rhs[k][j - 1][i][2]
					- lhs_buf[k][j][i][AA][3][4]*rhs[k][j - 1][i][3]
					- lhs_buf[k][j][i][AA][4][4]*rhs[k][j - 1][i][4];
					
					


				
					lhs_buf[k][j][i][BB][0][0] = lhs_buf[k][j][i][BB][0][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j - 1][i][CC][0][0]
						  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j - 1][i][CC][0][1]
						  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j - 1][i][CC][0][2]
						  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j - 1][i][CC][0][3]
						  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j - 1][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][1] = lhs_buf[k][j][i][BB][0][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j - 1][i][CC][0][0]
						  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j - 1][i][CC][0][1]
						  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j - 1][i][CC][0][2]
						  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j - 1][i][CC][0][3]
						  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j - 1][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][2] = lhs_buf[k][j][i][BB][0][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j - 1][i][CC][0][0]
						  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j - 1][i][CC][0][1]
						  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j - 1][i][CC][0][2]
						  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j - 1][i][CC][0][3]
						  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j - 1][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][3] = lhs_buf[k][j][i][BB][0][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j - 1][i][CC][0][0]
						  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j - 1][i][CC][0][1]
						  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j - 1][i][CC][0][2]
						  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j - 1][i][CC][0][3]
						  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j - 1][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][4] = lhs_buf[k][j][i][BB][0][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j - 1][i][CC][0][0]
						  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j - 1][i][CC][0][1]
						  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j - 1][i][CC][0][2]
						  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j - 1][i][CC][0][3]
						  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j - 1][i][CC][0][4];
					lhs_buf[k][j][i][BB][1][0] = lhs_buf[k][j][i][BB][1][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j - 1][i][CC][1][0]
						  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j - 1][i][CC][1][1]
						  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j - 1][i][CC][1][2]
						  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j - 1][i][CC][1][3]
						  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j - 1][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][1] = lhs_buf[k][j][i][BB][1][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j - 1][i][CC][1][0]
						  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j - 1][i][CC][1][1]
						  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j - 1][i][CC][1][2]
						  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j - 1][i][CC][1][3]
						  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j - 1][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][2] = lhs_buf[k][j][i][BB][1][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j - 1][i][CC][1][0]
						  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j - 1][i][CC][1][1]
						  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j - 1][i][CC][1][2]
						  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j - 1][i][CC][1][3]
						  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j - 1][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][3] = lhs_buf[k][j][i][BB][1][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j - 1][i][CC][1][0]
						  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j - 1][i][CC][1][1]
						  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j - 1][i][CC][1][2]
						  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j - 1][i][CC][1][3]
						  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j - 1][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][4] = lhs_buf[k][j][i][BB][1][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j - 1][i][CC][1][0]
						  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j - 1][i][CC][1][1]
						  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j - 1][i][CC][1][2]
						  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j - 1][i][CC][1][3]
						  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j - 1][i][CC][1][4];
					lhs_buf[k][j][i][BB][2][0] = lhs_buf[k][j][i][BB][2][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j - 1][i][CC][2][0]
						  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j - 1][i][CC][2][1]
						  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j - 1][i][CC][2][2]
						  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j - 1][i][CC][2][3]
						  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j - 1][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][1] = lhs_buf[k][j][i][BB][2][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j - 1][i][CC][2][0]
						  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j - 1][i][CC][2][1]
						  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j - 1][i][CC][2][2]
						  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j - 1][i][CC][2][3]
						  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j - 1][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][2] = lhs_buf[k][j][i][BB][2][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j - 1][i][CC][2][0]
						  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j - 1][i][CC][2][1]
						  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j - 1][i][CC][2][2]
						  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j - 1][i][CC][2][3]
						  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j - 1][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][3] = lhs_buf[k][j][i][BB][2][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j - 1][i][CC][2][0]
						  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j - 1][i][CC][2][1]
						  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j - 1][i][CC][2][2]
						  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j - 1][i][CC][2][3]
						  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j - 1][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][4] = lhs_buf[k][j][i][BB][2][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j - 1][i][CC][2][0]
						  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j - 1][i][CC][2][1]
						  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j - 1][i][CC][2][2]
						  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j - 1][i][CC][2][3]
						  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j - 1][i][CC][2][4];
					lhs_buf[k][j][i][BB][3][0] = lhs_buf[k][j][i][BB][3][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j - 1][i][CC][3][0]
						  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j - 1][i][CC][3][1]
						  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j - 1][i][CC][3][2]
						  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j - 1][i][CC][3][3]
						  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j - 1][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][1] = lhs_buf[k][j][i][BB][3][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j - 1][i][CC][3][0]
						  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j - 1][i][CC][3][1]
						  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j - 1][i][CC][3][2]
						  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j - 1][i][CC][3][3]
						  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j - 1][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][2] = lhs_buf[k][j][i][BB][3][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j - 1][i][CC][3][0]
						  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j - 1][i][CC][3][1]
						  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j - 1][i][CC][3][2]
						  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j - 1][i][CC][3][3]
						  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j - 1][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][3] = lhs_buf[k][j][i][BB][3][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j - 1][i][CC][3][0]
						  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j - 1][i][CC][3][1]
						  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j - 1][i][CC][3][2]
						  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j - 1][i][CC][3][3]
						  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j - 1][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][4] = lhs_buf[k][j][i][BB][3][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j - 1][i][CC][3][0]
						  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j - 1][i][CC][3][1]
						  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j - 1][i][CC][3][2]
						  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j - 1][i][CC][3][3]
						  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j - 1][i][CC][3][4];
					lhs_buf[k][j][i][BB][4][0] = lhs_buf[k][j][i][BB][4][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j - 1][i][CC][4][0]
						  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j - 1][i][CC][4][1]
						  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j - 1][i][CC][4][2]
						  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j - 1][i][CC][4][3]
						  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j - 1][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][1] = lhs_buf[k][j][i][BB][4][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j - 1][i][CC][4][0]
						  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j - 1][i][CC][4][1]
						  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j - 1][i][CC][4][2]
						  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j - 1][i][CC][4][3]
						  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j - 1][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][2] = lhs_buf[k][j][i][BB][4][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j - 1][i][CC][4][0]
						  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j - 1][i][CC][4][1]
						  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j - 1][i][CC][4][2]
						  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j - 1][i][CC][4][3]
						  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j - 1][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][3] = lhs_buf[k][j][i][BB][4][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j - 1][i][CC][4][0]
						  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j - 1][i][CC][4][1]
						  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j - 1][i][CC][4][2]
						  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j - 1][i][CC][4][3]
						  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j - 1][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][4] = lhs_buf[k][j][i][BB][4][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j - 1][i][CC][4][0]
						  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j - 1][i][CC][4][1]
						  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j - 1][i][CC][4][2]
						  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j - 1][i][CC][4][3]
						  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j - 1][i][CC][4][4];

					
					pivot = 1.00/lhs_buf[k][j][i][BB][0][0];
					lhs_buf[k][j][i][BB][1][0] = lhs_buf[k][j][i][BB][1][0]*pivot;
					lhs_buf[k][j][i][BB][2][0] = lhs_buf[k][j][i][BB][2][0]*pivot;
					lhs_buf[k][j][i][BB][3][0] = lhs_buf[k][j][i][BB][3][0]*pivot;
					lhs_buf[k][j][i][BB][4][0] = lhs_buf[k][j][i][BB][4][0]*pivot;
					lhs_buf[k][j][i][CC][0][0] = lhs_buf[k][j][i][CC][0][0]*pivot;
					lhs_buf[k][j][i][CC][1][0] = lhs_buf[k][j][i][CC][1][0]*pivot;
					lhs_buf[k][j][i][CC][2][0] = lhs_buf[k][j][i][CC][2][0]*pivot;
					lhs_buf[k][j][i][CC][3][0] = lhs_buf[k][j][i][CC][3][0]*pivot;
					lhs_buf[k][j][i][CC][4][0] = lhs_buf[k][j][i][CC][4][0]*pivot;
					rhs[k][j][i][0]   = rhs[k][j][i][0]  *pivot;

					coeff = lhs_buf[k][j][i][BB][0][1];
					lhs_buf[k][j][i][BB][1][1]= lhs_buf[k][j][i][BB][1][1] - coeff*lhs_buf[k][j][i][BB][1][0];
					lhs_buf[k][j][i][BB][2][1]= lhs_buf[k][j][i][BB][2][1] - coeff*lhs_buf[k][j][i][BB][2][0];
					lhs_buf[k][j][i][BB][3][1]= lhs_buf[k][j][i][BB][3][1] - coeff*lhs_buf[k][j][i][BB][3][0];
					lhs_buf[k][j][i][BB][4][1]= lhs_buf[k][j][i][BB][4][1] - coeff*lhs_buf[k][j][i][BB][4][0];
					lhs_buf[k][j][i][CC][0][1] = lhs_buf[k][j][i][CC][0][1] - coeff*lhs_buf[k][j][i][CC][0][0];
					lhs_buf[k][j][i][CC][1][1] = lhs_buf[k][j][i][CC][1][1] - coeff*lhs_buf[k][j][i][CC][1][0];
					lhs_buf[k][j][i][CC][2][1] = lhs_buf[k][j][i][CC][2][1] - coeff*lhs_buf[k][j][i][CC][2][0];
					lhs_buf[k][j][i][CC][3][1] = lhs_buf[k][j][i][CC][3][1] - coeff*lhs_buf[k][j][i][CC][3][0];
					lhs_buf[k][j][i][CC][4][1] = lhs_buf[k][j][i][CC][4][1] - coeff*lhs_buf[k][j][i][CC][4][0];
					rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff*rhs[k][j][i][0];

					coeff = lhs_buf[k][j][i][BB][0][2];
					lhs_buf[k][j][i][BB][1][2]= lhs_buf[k][j][i][BB][1][2] - coeff*lhs_buf[k][j][i][BB][1][0];
					lhs_buf[k][j][i][BB][2][2]= lhs_buf[k][j][i][BB][2][2] - coeff*lhs_buf[k][j][i][BB][2][0];
					lhs_buf[k][j][i][BB][3][2]= lhs_buf[k][j][i][BB][3][2] - coeff*lhs_buf[k][j][i][BB][3][0];
					lhs_buf[k][j][i][BB][4][2]= lhs_buf[k][j][i][BB][4][2] - coeff*lhs_buf[k][j][i][BB][4][0];
					lhs_buf[k][j][i][CC][0][2] = lhs_buf[k][j][i][CC][0][2] - coeff*lhs_buf[k][j][i][CC][0][0];
					lhs_buf[k][j][i][CC][1][2] = lhs_buf[k][j][i][CC][1][2] - coeff*lhs_buf[k][j][i][CC][1][0];
					lhs_buf[k][j][i][CC][2][2] = lhs_buf[k][j][i][CC][2][2] - coeff*lhs_buf[k][j][i][CC][2][0];
					lhs_buf[k][j][i][CC][3][2] = lhs_buf[k][j][i][CC][3][2] - coeff*lhs_buf[k][j][i][CC][3][0];
					lhs_buf[k][j][i][CC][4][2] = lhs_buf[k][j][i][CC][4][2] - coeff*lhs_buf[k][j][i][CC][4][0];
					rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff*rhs[k][j][i][0];

					coeff = lhs_buf[k][j][i][BB][0][3];
					lhs_buf[k][j][i][BB][1][3]= lhs_buf[k][j][i][BB][1][3] - coeff*lhs_buf[k][j][i][BB][1][0];
					lhs_buf[k][j][i][BB][2][3]= lhs_buf[k][j][i][BB][2][3] - coeff*lhs_buf[k][j][i][BB][2][0];
					lhs_buf[k][j][i][BB][3][3]= lhs_buf[k][j][i][BB][3][3] - coeff*lhs_buf[k][j][i][BB][3][0];
					lhs_buf[k][j][i][BB][4][3]= lhs_buf[k][j][i][BB][4][3] - coeff*lhs_buf[k][j][i][BB][4][0];
					lhs_buf[k][j][i][CC][0][3] = lhs_buf[k][j][i][CC][0][3] - coeff*lhs_buf[k][j][i][CC][0][0];
					lhs_buf[k][j][i][CC][1][3] = lhs_buf[k][j][i][CC][1][3] - coeff*lhs_buf[k][j][i][CC][1][0];
					lhs_buf[k][j][i][CC][2][3] = lhs_buf[k][j][i][CC][2][3] - coeff*lhs_buf[k][j][i][CC][2][0];
					lhs_buf[k][j][i][CC][3][3] = lhs_buf[k][j][i][CC][3][3] - coeff*lhs_buf[k][j][i][CC][3][0];
					lhs_buf[k][j][i][CC][4][3] = lhs_buf[k][j][i][CC][4][3] - coeff*lhs_buf[k][j][i][CC][4][0];
					rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff*rhs[k][j][i][0];

					coeff = lhs_buf[k][j][i][BB][0][4];
					lhs_buf[k][j][i][BB][1][4]= lhs_buf[k][j][i][BB][1][4] - coeff*lhs_buf[k][j][i][BB][1][0];
					lhs_buf[k][j][i][BB][2][4]= lhs_buf[k][j][i][BB][2][4] - coeff*lhs_buf[k][j][i][BB][2][0];
					lhs_buf[k][j][i][BB][3][4]= lhs_buf[k][j][i][BB][3][4] - coeff*lhs_buf[k][j][i][BB][3][0];
					lhs_buf[k][j][i][BB][4][4]= lhs_buf[k][j][i][BB][4][4] - coeff*lhs_buf[k][j][i][BB][4][0];
					lhs_buf[k][j][i][CC][0][4] = lhs_buf[k][j][i][CC][0][4] - coeff*lhs_buf[k][j][i][CC][0][0];
					lhs_buf[k][j][i][CC][1][4] = lhs_buf[k][j][i][CC][1][4] - coeff*lhs_buf[k][j][i][CC][1][0];
					lhs_buf[k][j][i][CC][2][4] = lhs_buf[k][j][i][CC][2][4] - coeff*lhs_buf[k][j][i][CC][2][0];
					lhs_buf[k][j][i][CC][3][4] = lhs_buf[k][j][i][CC][3][4] - coeff*lhs_buf[k][j][i][CC][3][0];
					lhs_buf[k][j][i][CC][4][4] = lhs_buf[k][j][i][CC][4][4] - coeff*lhs_buf[k][j][i][CC][4][0];
					rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff*rhs[k][j][i][0];


					pivot = 1.00/lhs_buf[k][j][i][BB][1][1];
					lhs_buf[k][j][i][BB][2][1] = lhs_buf[k][j][i][BB][2][1]*pivot;
					lhs_buf[k][j][i][BB][3][1] = lhs_buf[k][j][i][BB][3][1]*pivot;
					lhs_buf[k][j][i][BB][4][1] = lhs_buf[k][j][i][BB][4][1]*pivot;
					lhs_buf[k][j][i][CC][0][1] = lhs_buf[k][j][i][CC][0][1]*pivot;
					lhs_buf[k][j][i][CC][1][1] = lhs_buf[k][j][i][CC][1][1]*pivot;
					lhs_buf[k][j][i][CC][2][1] = lhs_buf[k][j][i][CC][2][1]*pivot;
					lhs_buf[k][j][i][CC][3][1] = lhs_buf[k][j][i][CC][3][1]*pivot;
					lhs_buf[k][j][i][CC][4][1] = lhs_buf[k][j][i][CC][4][1]*pivot;
					rhs[k][j][i][1]   = rhs[k][j][i][1]  *pivot;

					coeff = lhs_buf[k][j][i][BB][1][0];
					lhs_buf[k][j][i][BB][2][0]= lhs_buf[k][j][i][BB][2][0] - coeff*lhs_buf[k][j][i][BB][2][1];
					lhs_buf[k][j][i][BB][3][0]= lhs_buf[k][j][i][BB][3][0] - coeff*lhs_buf[k][j][i][BB][3][1];
					lhs_buf[k][j][i][BB][4][0]= lhs_buf[k][j][i][BB][4][0] - coeff*lhs_buf[k][j][i][BB][4][1];
					lhs_buf[k][j][i][CC][0][0] = lhs_buf[k][j][i][CC][0][0] - coeff*lhs_buf[k][j][i][CC][0][1];
					lhs_buf[k][j][i][CC][1][0] = lhs_buf[k][j][i][CC][1][0] - coeff*lhs_buf[k][j][i][CC][1][1];
					lhs_buf[k][j][i][CC][2][0] = lhs_buf[k][j][i][CC][2][0] - coeff*lhs_buf[k][j][i][CC][2][1];
					lhs_buf[k][j][i][CC][3][0] = lhs_buf[k][j][i][CC][3][0] - coeff*lhs_buf[k][j][i][CC][3][1];
					lhs_buf[k][j][i][CC][4][0] = lhs_buf[k][j][i][CC][4][0] - coeff*lhs_buf[k][j][i][CC][4][1];
					rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff*rhs[k][j][i][1];

					coeff = lhs_buf[k][j][i][BB][1][2];
					lhs_buf[k][j][i][BB][2][2]= lhs_buf[k][j][i][BB][2][2] - coeff*lhs_buf[k][j][i][BB][2][1];
					lhs_buf[k][j][i][BB][3][2]= lhs_buf[k][j][i][BB][3][2] - coeff*lhs_buf[k][j][i][BB][3][1];
					lhs_buf[k][j][i][BB][4][2]= lhs_buf[k][j][i][BB][4][2] - coeff*lhs_buf[k][j][i][BB][4][1];
					lhs_buf[k][j][i][CC][0][2] = lhs_buf[k][j][i][CC][0][2] - coeff*lhs_buf[k][j][i][CC][0][1];
					lhs_buf[k][j][i][CC][1][2] = lhs_buf[k][j][i][CC][1][2] - coeff*lhs_buf[k][j][i][CC][1][1];
					lhs_buf[k][j][i][CC][2][2] = lhs_buf[k][j][i][CC][2][2] - coeff*lhs_buf[k][j][i][CC][2][1];
					lhs_buf[k][j][i][CC][3][2] = lhs_buf[k][j][i][CC][3][2] - coeff*lhs_buf[k][j][i][CC][3][1];
					lhs_buf[k][j][i][CC][4][2] = lhs_buf[k][j][i][CC][4][2] - coeff*lhs_buf[k][j][i][CC][4][1];
					rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff*rhs[k][j][i][1];

					coeff = lhs_buf[k][j][i][BB][1][3];
					lhs_buf[k][j][i][BB][2][3]= lhs_buf[k][j][i][BB][2][3] - coeff*lhs_buf[k][j][i][BB][2][1];
					lhs_buf[k][j][i][BB][3][3]= lhs_buf[k][j][i][BB][3][3] - coeff*lhs_buf[k][j][i][BB][3][1];
					lhs_buf[k][j][i][BB][4][3]= lhs_buf[k][j][i][BB][4][3] - coeff*lhs_buf[k][j][i][BB][4][1];
					lhs_buf[k][j][i][CC][0][3] = lhs_buf[k][j][i][CC][0][3] - coeff*lhs_buf[k][j][i][CC][0][1];
					lhs_buf[k][j][i][CC][1][3] = lhs_buf[k][j][i][CC][1][3] - coeff*lhs_buf[k][j][i][CC][1][1];
					lhs_buf[k][j][i][CC][2][3] = lhs_buf[k][j][i][CC][2][3] - coeff*lhs_buf[k][j][i][CC][2][1];
					lhs_buf[k][j][i][CC][3][3] = lhs_buf[k][j][i][CC][3][3] - coeff*lhs_buf[k][j][i][CC][3][1];
					lhs_buf[k][j][i][CC][4][3] = lhs_buf[k][j][i][CC][4][3] - coeff*lhs_buf[k][j][i][CC][4][1];
					rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff*rhs[k][j][i][1];

					coeff = lhs_buf[k][j][i][BB][1][4];
					lhs_buf[k][j][i][BB][2][4]= lhs_buf[k][j][i][BB][2][4] - coeff*lhs_buf[k][j][i][BB][2][1];
					lhs_buf[k][j][i][BB][3][4]= lhs_buf[k][j][i][BB][3][4] - coeff*lhs_buf[k][j][i][BB][3][1];
					lhs_buf[k][j][i][BB][4][4]= lhs_buf[k][j][i][BB][4][4] - coeff*lhs_buf[k][j][i][BB][4][1];
					lhs_buf[k][j][i][CC][0][4] = lhs_buf[k][j][i][CC][0][4] - coeff*lhs_buf[k][j][i][CC][0][1];
					lhs_buf[k][j][i][CC][1][4] = lhs_buf[k][j][i][CC][1][4] - coeff*lhs_buf[k][j][i][CC][1][1];
					lhs_buf[k][j][i][CC][2][4] = lhs_buf[k][j][i][CC][2][4] - coeff*lhs_buf[k][j][i][CC][2][1];
					lhs_buf[k][j][i][CC][3][4] = lhs_buf[k][j][i][CC][3][4] - coeff*lhs_buf[k][j][i][CC][3][1];
					lhs_buf[k][j][i][CC][4][4] = lhs_buf[k][j][i][CC][4][4] - coeff*lhs_buf[k][j][i][CC][4][1];
					rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff*rhs[k][j][i][1];


					pivot = 1.00/lhs_buf[k][j][i][BB][2][2];
					lhs_buf[k][j][i][BB][3][2] = lhs_buf[k][j][i][BB][3][2]*pivot;
					lhs_buf[k][j][i][BB][4][2] = lhs_buf[k][j][i][BB][4][2]*pivot;
					lhs_buf[k][j][i][CC][0][2] = lhs_buf[k][j][i][CC][0][2]*pivot;
					lhs_buf[k][j][i][CC][1][2] = lhs_buf[k][j][i][CC][1][2]*pivot;
					lhs_buf[k][j][i][CC][2][2] = lhs_buf[k][j][i][CC][2][2]*pivot;
					lhs_buf[k][j][i][CC][3][2] = lhs_buf[k][j][i][CC][3][2]*pivot;
					lhs_buf[k][j][i][CC][4][2] = lhs_buf[k][j][i][CC][4][2]*pivot;
					rhs[k][j][i][2]   = rhs[k][j][i][2]  *pivot;

					coeff = lhs_buf[k][j][i][BB][2][0];
					lhs_buf[k][j][i][BB][3][0]= lhs_buf[k][j][i][BB][3][0] - coeff*lhs_buf[k][j][i][BB][3][2];
					lhs_buf[k][j][i][BB][4][0]= lhs_buf[k][j][i][BB][4][0] - coeff*lhs_buf[k][j][i][BB][4][2];
					lhs_buf[k][j][i][CC][0][0] = lhs_buf[k][j][i][CC][0][0] - coeff*lhs_buf[k][j][i][CC][0][2];
					lhs_buf[k][j][i][CC][1][0] = lhs_buf[k][j][i][CC][1][0] - coeff*lhs_buf[k][j][i][CC][1][2];
					lhs_buf[k][j][i][CC][2][0] = lhs_buf[k][j][i][CC][2][0] - coeff*lhs_buf[k][j][i][CC][2][2];
					lhs_buf[k][j][i][CC][3][0] = lhs_buf[k][j][i][CC][3][0] - coeff*lhs_buf[k][j][i][CC][3][2];
					lhs_buf[k][j][i][CC][4][0] = lhs_buf[k][j][i][CC][4][0] - coeff*lhs_buf[k][j][i][CC][4][2];
					rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff*rhs[k][j][i][2];

					coeff = lhs_buf[k][j][i][BB][2][1];
					lhs_buf[k][j][i][BB][3][1]= lhs_buf[k][j][i][BB][3][1] - coeff*lhs_buf[k][j][i][BB][3][2];
					lhs_buf[k][j][i][BB][4][1]= lhs_buf[k][j][i][BB][4][1] - coeff*lhs_buf[k][j][i][BB][4][2];
					lhs_buf[k][j][i][CC][0][1] = lhs_buf[k][j][i][CC][0][1] - coeff*lhs_buf[k][j][i][CC][0][2];
					lhs_buf[k][j][i][CC][1][1] = lhs_buf[k][j][i][CC][1][1] - coeff*lhs_buf[k][j][i][CC][1][2];
					lhs_buf[k][j][i][CC][2][1] = lhs_buf[k][j][i][CC][2][1] - coeff*lhs_buf[k][j][i][CC][2][2];
					lhs_buf[k][j][i][CC][3][1] = lhs_buf[k][j][i][CC][3][1] - coeff*lhs_buf[k][j][i][CC][3][2];
					lhs_buf[k][j][i][CC][4][1] = lhs_buf[k][j][i][CC][4][1] - coeff*lhs_buf[k][j][i][CC][4][2];
					rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff*rhs[k][j][i][2];

					coeff = lhs_buf[k][j][i][BB][2][3];
					lhs_buf[k][j][i][BB][3][3]= lhs_buf[k][j][i][BB][3][3] - coeff*lhs_buf[k][j][i][BB][3][2];
					lhs_buf[k][j][i][BB][4][3]= lhs_buf[k][j][i][BB][4][3] - coeff*lhs_buf[k][j][i][BB][4][2];
					lhs_buf[k][j][i][CC][0][3] = lhs_buf[k][j][i][CC][0][3] - coeff*lhs_buf[k][j][i][CC][0][2];
					lhs_buf[k][j][i][CC][1][3] = lhs_buf[k][j][i][CC][1][3] - coeff*lhs_buf[k][j][i][CC][1][2];
					lhs_buf[k][j][i][CC][2][3] = lhs_buf[k][j][i][CC][2][3] - coeff*lhs_buf[k][j][i][CC][2][2];
					lhs_buf[k][j][i][CC][3][3] = lhs_buf[k][j][i][CC][3][3] - coeff*lhs_buf[k][j][i][CC][3][2];
					lhs_buf[k][j][i][CC][4][3] = lhs_buf[k][j][i][CC][4][3] - coeff*lhs_buf[k][j][i][CC][4][2];
					rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff*rhs[k][j][i][2];

					coeff = lhs_buf[k][j][i][BB][2][4];
					lhs_buf[k][j][i][BB][3][4]= lhs_buf[k][j][i][BB][3][4] - coeff*lhs_buf[k][j][i][BB][3][2];
					lhs_buf[k][j][i][BB][4][4]= lhs_buf[k][j][i][BB][4][4] - coeff*lhs_buf[k][j][i][BB][4][2];
					lhs_buf[k][j][i][CC][0][4] = lhs_buf[k][j][i][CC][0][4] - coeff*lhs_buf[k][j][i][CC][0][2];
					lhs_buf[k][j][i][CC][1][4] = lhs_buf[k][j][i][CC][1][4] - coeff*lhs_buf[k][j][i][CC][1][2];
					lhs_buf[k][j][i][CC][2][4] = lhs_buf[k][j][i][CC][2][4] - coeff*lhs_buf[k][j][i][CC][2][2];
					lhs_buf[k][j][i][CC][3][4] = lhs_buf[k][j][i][CC][3][4] - coeff*lhs_buf[k][j][i][CC][3][2];
					lhs_buf[k][j][i][CC][4][4] = lhs_buf[k][j][i][CC][4][4] - coeff*lhs_buf[k][j][i][CC][4][2];
					rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff*rhs[k][j][i][2];


					pivot = 1.00/lhs_buf[k][j][i][BB][3][3];
					lhs_buf[k][j][i][BB][4][3] = lhs_buf[k][j][i][BB][4][3]*pivot;
					lhs_buf[k][j][i][CC][0][3] = lhs_buf[k][j][i][CC][0][3]*pivot;
					lhs_buf[k][j][i][CC][1][3] = lhs_buf[k][j][i][CC][1][3]*pivot;
					lhs_buf[k][j][i][CC][2][3] = lhs_buf[k][j][i][CC][2][3]*pivot;
					lhs_buf[k][j][i][CC][3][3] = lhs_buf[k][j][i][CC][3][3]*pivot;
					lhs_buf[k][j][i][CC][4][3] = lhs_buf[k][j][i][CC][4][3]*pivot;
					rhs[k][j][i][3]   = rhs[k][j][i][3]  *pivot;

					coeff = lhs_buf[k][j][i][BB][3][0];
					lhs_buf[k][j][i][BB][4][0]= lhs_buf[k][j][i][BB][4][0] - coeff*lhs_buf[k][j][i][BB][4][3];
					lhs_buf[k][j][i][CC][0][0] = lhs_buf[k][j][i][CC][0][0] - coeff*lhs_buf[k][j][i][CC][0][3];
					lhs_buf[k][j][i][CC][1][0] = lhs_buf[k][j][i][CC][1][0] - coeff*lhs_buf[k][j][i][CC][1][3];
					lhs_buf[k][j][i][CC][2][0] = lhs_buf[k][j][i][CC][2][0] - coeff*lhs_buf[k][j][i][CC][2][3];
					lhs_buf[k][j][i][CC][3][0] = lhs_buf[k][j][i][CC][3][0] - coeff*lhs_buf[k][j][i][CC][3][3];
					lhs_buf[k][j][i][CC][4][0] = lhs_buf[k][j][i][CC][4][0] - coeff*lhs_buf[k][j][i][CC][4][3];
					rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff*rhs[k][j][i][3];

					coeff = lhs_buf[k][j][i][BB][3][1];
					lhs_buf[k][j][i][BB][4][1]= lhs_buf[k][j][i][BB][4][1] - coeff*lhs_buf[k][j][i][BB][4][3];
					lhs_buf[k][j][i][CC][0][1] = lhs_buf[k][j][i][CC][0][1] - coeff*lhs_buf[k][j][i][CC][0][3];
					lhs_buf[k][j][i][CC][1][1] = lhs_buf[k][j][i][CC][1][1] - coeff*lhs_buf[k][j][i][CC][1][3];
					lhs_buf[k][j][i][CC][2][1] = lhs_buf[k][j][i][CC][2][1] - coeff*lhs_buf[k][j][i][CC][2][3];
					lhs_buf[k][j][i][CC][3][1] = lhs_buf[k][j][i][CC][3][1] - coeff*lhs_buf[k][j][i][CC][3][3];
					lhs_buf[k][j][i][CC][4][1] = lhs_buf[k][j][i][CC][4][1] - coeff*lhs_buf[k][j][i][CC][4][3];
					rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff*rhs[k][j][i][3];

					coeff = lhs_buf[k][j][i][BB][3][2];
					lhs_buf[k][j][i][BB][4][2]= lhs_buf[k][j][i][BB][4][2] - coeff*lhs_buf[k][j][i][BB][4][3];
					lhs_buf[k][j][i][CC][0][2] = lhs_buf[k][j][i][CC][0][2] - coeff*lhs_buf[k][j][i][CC][0][3];
					lhs_buf[k][j][i][CC][1][2] = lhs_buf[k][j][i][CC][1][2] - coeff*lhs_buf[k][j][i][CC][1][3];
					lhs_buf[k][j][i][CC][2][2] = lhs_buf[k][j][i][CC][2][2] - coeff*lhs_buf[k][j][i][CC][2][3];
					lhs_buf[k][j][i][CC][3][2] = lhs_buf[k][j][i][CC][3][2] - coeff*lhs_buf[k][j][i][CC][3][3];
					lhs_buf[k][j][i][CC][4][2] = lhs_buf[k][j][i][CC][4][2] - coeff*lhs_buf[k][j][i][CC][4][3];
					rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff*rhs[k][j][i][3];

					coeff = lhs_buf[k][j][i][BB][3][4];
					lhs_buf[k][j][i][BB][4][4]= lhs_buf[k][j][i][BB][4][4] - coeff*lhs_buf[k][j][i][BB][4][3];
					lhs_buf[k][j][i][CC][0][4] = lhs_buf[k][j][i][CC][0][4] - coeff*lhs_buf[k][j][i][CC][0][3];
					lhs_buf[k][j][i][CC][1][4] = lhs_buf[k][j][i][CC][1][4] - coeff*lhs_buf[k][j][i][CC][1][3];
					lhs_buf[k][j][i][CC][2][4] = lhs_buf[k][j][i][CC][2][4] - coeff*lhs_buf[k][j][i][CC][2][3];
					lhs_buf[k][j][i][CC][3][4] = lhs_buf[k][j][i][CC][3][4] - coeff*lhs_buf[k][j][i][CC][3][3];
					lhs_buf[k][j][i][CC][4][4] = lhs_buf[k][j][i][CC][4][4] - coeff*lhs_buf[k][j][i][CC][4][3];
					rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff*rhs[k][j][i][3];


					pivot = 1.00/lhs_buf[k][j][i][BB][4][4];
					lhs_buf[k][j][i][CC][0][4] = lhs_buf[k][j][i][CC][0][4]*pivot;
					lhs_buf[k][j][i][CC][1][4] = lhs_buf[k][j][i][CC][1][4]*pivot;
					lhs_buf[k][j][i][CC][2][4] = lhs_buf[k][j][i][CC][2][4]*pivot;
					lhs_buf[k][j][i][CC][3][4] = lhs_buf[k][j][i][CC][3][4]*pivot;
					lhs_buf[k][j][i][CC][4][4] = lhs_buf[k][j][i][CC][4][4]*pivot;
					rhs[k][j][i][4]   = rhs[k][j][i][4]  *pivot;

					coeff = lhs_buf[k][j][i][BB][4][0];
					lhs_buf[k][j][i][CC][0][0] = lhs_buf[k][j][i][CC][0][0] - coeff*lhs_buf[k][j][i][CC][0][4];
					lhs_buf[k][j][i][CC][1][0] = lhs_buf[k][j][i][CC][1][0] - coeff*lhs_buf[k][j][i][CC][1][4];
					lhs_buf[k][j][i][CC][2][0] = lhs_buf[k][j][i][CC][2][0] - coeff*lhs_buf[k][j][i][CC][2][4];
					lhs_buf[k][j][i][CC][3][0] = lhs_buf[k][j][i][CC][3][0] - coeff*lhs_buf[k][j][i][CC][3][4];
					lhs_buf[k][j][i][CC][4][0] = lhs_buf[k][j][i][CC][4][0] - coeff*lhs_buf[k][j][i][CC][4][4];
					rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff*rhs[k][j][i][4];

					coeff = lhs_buf[k][j][i][BB][4][1];
					lhs_buf[k][j][i][CC][0][1] = lhs_buf[k][j][i][CC][0][1] - coeff*lhs_buf[k][j][i][CC][0][4];
					lhs_buf[k][j][i][CC][1][1] = lhs_buf[k][j][i][CC][1][1] - coeff*lhs_buf[k][j][i][CC][1][4];
					lhs_buf[k][j][i][CC][2][1] = lhs_buf[k][j][i][CC][2][1] - coeff*lhs_buf[k][j][i][CC][2][4];
					lhs_buf[k][j][i][CC][3][1] = lhs_buf[k][j][i][CC][3][1] - coeff*lhs_buf[k][j][i][CC][3][4];
					lhs_buf[k][j][i][CC][4][1] = lhs_buf[k][j][i][CC][4][1] - coeff*lhs_buf[k][j][i][CC][4][4];
					rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff*rhs[k][j][i][4];

					coeff = lhs_buf[k][j][i][BB][4][2];
					lhs_buf[k][j][i][CC][0][2] = lhs_buf[k][j][i][CC][0][2] - coeff*lhs_buf[k][j][i][CC][0][4];
					lhs_buf[k][j][i][CC][1][2] = lhs_buf[k][j][i][CC][1][2] - coeff*lhs_buf[k][j][i][CC][1][4];
					lhs_buf[k][j][i][CC][2][2] = lhs_buf[k][j][i][CC][2][2] - coeff*lhs_buf[k][j][i][CC][2][4];
					lhs_buf[k][j][i][CC][3][2] = lhs_buf[k][j][i][CC][3][2] - coeff*lhs_buf[k][j][i][CC][3][4];
					lhs_buf[k][j][i][CC][4][2] = lhs_buf[k][j][i][CC][4][2] - coeff*lhs_buf[k][j][i][CC][4][4];
					rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff*rhs[k][j][i][4];

					coeff = lhs_buf[k][j][i][BB][4][3];
					lhs_buf[k][j][i][CC][0][3] = lhs_buf[k][j][i][CC][0][3] - coeff*lhs_buf[k][j][i][CC][0][4];
					lhs_buf[k][j][i][CC][1][3] = lhs_buf[k][j][i][CC][1][3] - coeff*lhs_buf[k][j][i][CC][1][4];
					lhs_buf[k][j][i][CC][2][3] = lhs_buf[k][j][i][CC][2][3] - coeff*lhs_buf[k][j][i][CC][2][4];
					lhs_buf[k][j][i][CC][3][3] = lhs_buf[k][j][i][CC][3][3] - coeff*lhs_buf[k][j][i][CC][3][4];
					lhs_buf[k][j][i][CC][4][3] = lhs_buf[k][j][i][CC][4][3] - coeff*lhs_buf[k][j][i][CC][4][4];
					rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff*rhs[k][j][i][4];
				}
			
			}
		}
	}			


	for (k = 1; k <= grid_points[2] - 2; k++) {
		for (i = 1; i <= grid_points[0] - 2; i++) {				
			for (j = jsize - 1; j >= 0; j--) {
				for (m = 0; m < BLOCK_SIZE; m++) {
					for (n = 0; n < BLOCK_SIZE; n++) {
						rhs[k][j][i][m] = rhs[k][j][i][m]
							- lhs_buf[k][j][i][CC][n][m] * rhs[k][j + 1][i][n];
					}
				}
			}
		}
	}



	if (timeron) timer_stop(t_ysolve);

}

