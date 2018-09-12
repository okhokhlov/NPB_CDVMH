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

double lhs_z[2][3][5][5];
double rhs_z[2][5];
double rhsp_z[5];
int len_z = 4;
double fjac_z[4][5][5];
double njac_z[4][5][5];

void fn_init_z(int a, int i, int j, int k);
void buf_compute_z(int k, int a);
void binvcrhs_z(int a1, int a2, int b1, int b2, int c1);
void matvec_sub_z(int a1, int a2, int b1, int c1 );
void matmul_sub_z(int a1, int a2, int b1, int b2, int c1, int c2 );
void binvrhs_z(int a1, int a2, int b1);

//---------------------------------------------------------------------
// Performs line solves in Z direction by first factoring
// the block-tridiagonal matrix into an upper triangular matrix, 
// and then performing back substitution to solve for the unknow
// vectors of each line.  
// 
// Make sure we treat elements zero to cell_size in the direction
// of the sweep.
//---------------------------------------------------------------------
void z_solve()
{
	int i, j, k, m, n, ksize, l;
	double pivot, coeff;
	//int z, a = -1, b = 2;
	//---------------------------------------------------------------------
	//---------------------------------------------------------------------

	if (timeron) timer_start(t_zsolve);

	//---------------------------------------------------------------------
	//---------------------------------------------------------------------

	//---------------------------------------------------------------------
	// This function computes the left hand side for the three z-factors   
	//---------------------------------------------------------------------

	ksize = grid_points[2] - 1;

	//---------------------------------------------------------------------
	// Compute the indices for storing the block-diagonal matrix;
	// determine c (labeled f) and s jacobians
	//---------------------------------------------------------------------
	

	
	for (j = 1; j <= grid_points[1] - 2; j++) {
		for (i = 1; i <= grid_points[0] - 2; i++) {
			for (k = 0; k <= ksize; k++) {
				tmp1 = 1.0 / u[k][j][i][0];
				tmp2 = tmp1 * tmp1;
				tmp3 = tmp1 * tmp2;

				fjac[k][0][0] = 0.0;
				fjac[k][1][0] = 0.0;
				fjac[k][2][0] = 0.0;
				fjac[k][3][0] = 1.0;
				fjac[k][4][0] = 0.0;

				fjac[k][0][1] = -(u[k][j][i][1] * u[k][j][i][3]) * tmp2;
				fjac[k][1][1] = u[k][j][i][3] * tmp1;
				fjac[k][2][1] = 0.0;
				fjac[k][3][1] = u[k][j][i][1] * tmp1;
				fjac[k][4][1] = 0.0;

				fjac[k][0][2] = -(u[k][j][i][2] * u[k][j][i][3]) * tmp2;
				fjac[k][1][2] = 0.0;
				fjac[k][2][2] = u[k][j][i][3] * tmp1;
				fjac[k][3][2] = u[k][j][i][2] * tmp1;
				fjac[k][4][2] = 0.0;

				fjac[k][0][3] = -(u[k][j][i][3] * u[k][j][i][3] * tmp2)
					+ c2 * qs[k][j][i];
				fjac[k][1][3] = -c2 *  u[k][j][i][1] * tmp1;
				fjac[k][2][3] = -c2 *  u[k][j][i][2] * tmp1;
				fjac[k][3][3] = (2.0 - c2) *  u[k][j][i][3] * tmp1;
				fjac[k][4][3] = c2;

				fjac[k][0][4] = (c2 * 2.0 * square[k][j][i] - c1 * u[k][j][i][4])
					* u[k][j][i][3] * tmp2;
				fjac[k][1][4] = -c2 * (u[k][j][i][1] * u[k][j][i][3]) * tmp2;
				fjac[k][2][4] = -c2 * (u[k][j][i][2] * u[k][j][i][3]) * tmp2;
				fjac[k][3][4] = c1 * (u[k][j][i][4] * tmp1)
					- c2 * (qs[k][j][i] + u[k][j][i][3] * u[k][j][i][3] * tmp2);
				fjac[k][4][4] = c1 * u[k][j][i][3] * tmp1;

				njac[k][0][0] = 0.0;
				njac[k][1][0] = 0.0;
				njac[k][2][0] = 0.0;
				njac[k][3][0] = 0.0;
				njac[k][4][0] = 0.0;

				njac[k][0][1] = -c3c4 * tmp2 * u[k][j][i][1];
				njac[k][1][1] = c3c4 * tmp1;
				njac[k][2][1] = 0.0;
				njac[k][3][1] = 0.0;
				njac[k][4][1] = 0.0;

				njac[k][0][2] = -c3c4 * tmp2 * u[k][j][i][2];
				njac[k][1][2] = 0.0;
				njac[k][2][2] = c3c4 * tmp1;
				njac[k][3][2] = 0.0;
				njac[k][4][2] = 0.0;

				njac[k][0][3] = -con43 * c3c4 * tmp2 * u[k][j][i][3];
				njac[k][1][3] = 0.0;
				njac[k][2][3] = 0.0;
				njac[k][3][3] = con43 * c3 * c4 * tmp1;
				njac[k][4][3] = 0.0;

				njac[k][0][4] = -(c3c4
					- c1345) * tmp3 * (u[k][j][i][1] * u[k][j][i][1])
					- (c3c4 - c1345) * tmp3 * (u[k][j][i][2] * u[k][j][i][2])
					- (con43 * c3c4
						- c1345) * tmp3 * (u[k][j][i][3] * u[k][j][i][3])
					- c1345 * tmp2 * u[k][j][i][4];

				njac[k][1][4] = (c3c4 - c1345) * tmp2 * u[k][j][i][1];
				njac[k][2][4] = (c3c4 - c1345) * tmp2 * u[k][j][i][2];
				njac[k][3][4] = (con43 * c3c4
					- c1345) * tmp2 * u[k][j][i][3];
				njac[k][4][4] = (c1345)* tmp1;
			}

			
			for (k = 0; k <= ksize; k++)
			{
				if(k == 0)
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
				else if(k == ksize)
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
					tmp1 = dt * tz1;
					tmp2 = dt * tz2;

					lhs_buf[k][j ][i][AA][0][0] = -tmp2 * fjac[k - 1][0][0]
						- tmp1 * njac[k - 1][0][0]
						- tmp1 * dz1;
					lhs_buf[k][j ][i][AA][1][0] = -tmp2 * fjac[k - 1][1][0]
						- tmp1 * njac[k - 1][1][0];
					lhs_buf[k][j ][i][AA][2][0] = -tmp2 * fjac[k - 1][2][0]
						- tmp1 * njac[k - 1][2][0];
					lhs_buf[k][j ][i][AA][3][0] = -tmp2 * fjac[k - 1][3][0]
						- tmp1 * njac[k - 1][3][0];
					lhs_buf[k][j ][i][AA][4][0] = -tmp2 * fjac[k - 1][4][0]
						- tmp1 * njac[k - 1][4][0];

					lhs_buf[k][j ][i][AA][0][1] = -tmp2 * fjac[k - 1][0][1]
						- tmp1 * njac[k - 1][0][1];
					lhs_buf[k][j ][i][AA][1][1] = -tmp2 * fjac[k - 1][1][1]
						- tmp1 * njac[k - 1][1][1]
						- tmp1 * dz2;
					lhs_buf[k][j ][i][AA][2][1] = -tmp2 * fjac[k - 1][2][1]
						- tmp1 * njac[k - 1][2][1];
					lhs_buf[k][j ][i][AA][3][1] = -tmp2 * fjac[k - 1][3][1]
						- tmp1 * njac[k - 1][3][1];
					lhs_buf[k][j ][i][AA][4][1] = -tmp2 * fjac[k - 1][4][1]
						- tmp1 * njac[k - 1][4][1];

					lhs_buf[k][j ][i][AA][0][2] = -tmp2 * fjac[k - 1][0][2]
						- tmp1 * njac[k - 1][0][2];
					lhs_buf[k][j ][i][AA][1][2] = -tmp2 * fjac[k - 1][1][2]
						- tmp1 * njac[k - 1][1][2];
					lhs_buf[k][j ][i][AA][2][2] = -tmp2 * fjac[k - 1][2][2]
						- tmp1 * njac[k - 1][2][2]
						- tmp1 * dz3;
					lhs_buf[k][j ][i][AA][3][2] = -tmp2 * fjac[k - 1][3][2]
						- tmp1 * njac[k - 1][3][2];
					lhs_buf[k][j ][i][AA][4][2] = -tmp2 * fjac[k - 1][4][2]
						- tmp1 * njac[k - 1][4][2];

					lhs_buf[k][j ][i][AA][0][3] = -tmp2 * fjac[k - 1][0][3]
						- tmp1 * njac[k - 1][0][3];
					lhs_buf[k][j ][i][AA][1][3] = -tmp2 * fjac[k - 1][1][3]
						- tmp1 * njac[k - 1][1][3];
					lhs_buf[k][j ][i][AA][2][3] = -tmp2 * fjac[k - 1][2][3]
						- tmp1 * njac[k - 1][2][3];
					lhs_buf[k][j ][i][AA][3][3] = -tmp2 * fjac[k - 1][3][3]
						- tmp1 * njac[k - 1][3][3]
						- tmp1 * dz4;
					lhs_buf[k][j ][i][AA][4][3] = -tmp2 * fjac[k - 1][4][3]
						- tmp1 * njac[k - 1][4][3];

					lhs_buf[k][j ][i][AA][0][4] = -tmp2 * fjac[k - 1][0][4]
						- tmp1 * njac[k - 1][0][4];
					lhs_buf[k][j ][i][AA][1][4] = -tmp2 * fjac[k - 1][1][4]
						- tmp1 * njac[k - 1][1][4];
					lhs_buf[k][j ][i][AA][2][4] = -tmp2 * fjac[k - 1][2][4]
						- tmp1 * njac[k - 1][2][4];
					lhs_buf[k][j ][i][AA][3][4] = -tmp2 * fjac[k - 1][3][4]
						- tmp1 * njac[k - 1][3][4];
					lhs_buf[k][j ][i][AA][4][4] = -tmp2 * fjac[k - 1][4][4]
						- tmp1 * njac[k - 1][4][4]
						- tmp1 * dz5;

					lhs_buf[k][j ][i][BB][0][0] = 1.0
						+ tmp1 * 2.0 * njac[k][0][0]
						+ tmp1 * 2.0 * dz1;
					lhs_buf[k][j ][i][BB][1][0] = tmp1 * 2.0 * njac[k][1][0];
					lhs_buf[k][j ][i][BB][2][0] = tmp1 * 2.0 * njac[k][2][0];
					lhs_buf[k][j ][i][BB][3][0] = tmp1 * 2.0 * njac[k][3][0];
					lhs_buf[k][j ][i][BB][4][0] = tmp1 * 2.0 * njac[k][4][0];

					lhs_buf[k][j ][i][BB][0][1] = tmp1 * 2.0 * njac[k][0][1];
					lhs_buf[k][j ][i][BB][1][1] = 1.0
						+ tmp1 * 2.0 * njac[k][1][1]
						+ tmp1 * 2.0 * dz2;
					lhs_buf[k][j ][i][BB][2][1] = tmp1 * 2.0 * njac[k][2][1];
					lhs_buf[k][j ][i][BB][3][1] = tmp1 * 2.0 * njac[k][3][1];
					lhs_buf[k][j ][i][BB][4][1] = tmp1 * 2.0 * njac[k][4][1];

					lhs_buf[k][j ][i][BB][0][2] = tmp1 * 2.0 * njac[k][0][2];
					lhs_buf[k][j ][i][BB][1][2] = tmp1 * 2.0 * njac[k][1][2];
					lhs_buf[k][j ][i][BB][2][2] = 1.0
						+ tmp1 * 2.0 * njac[k][2][2]
						+ tmp1 * 2.0 * dz3;
					lhs_buf[k][j ][i][BB][3][2] = tmp1 * 2.0 * njac[k][3][2];
					lhs_buf[k][j ][i][BB][4][2] = tmp1 * 2.0 * njac[k][4][2];

					lhs_buf[k][j ][i][BB][0][3] = tmp1 * 2.0 * njac[k][0][3];
					lhs_buf[k][j ][i][BB][1][3] = tmp1 * 2.0 * njac[k][1][3];
					lhs_buf[k][j ][i][BB][2][3] = tmp1 * 2.0 * njac[k][2][3];
					lhs_buf[k][j ][i][BB][3][3] = 1.0
						+ tmp1 * 2.0 * njac[k][3][3]
						+ tmp1 * 2.0 * dz4;
					lhs_buf[k][j ][i][BB][4][3] = tmp1 * 2.0 * njac[k][4][3];

					lhs_buf[k][j ][i][BB][0][4] = tmp1 * 2.0 * njac[k][0][4];
					lhs_buf[k][j ][i][BB][1][4] = tmp1 * 2.0 * njac[k][1][4];
					lhs_buf[k][j ][i][BB][2][4] = tmp1 * 2.0 * njac[k][2][4];
					lhs_buf[k][j ][i][BB][3][4] = tmp1 * 2.0 * njac[k][3][4];
					lhs_buf[k][j ][i][BB][4][4] = 1.0
						+ tmp1 * 2.0 * njac[k][4][4]
						+ tmp1 * 2.0 * dz5;

					lhs_buf[k][j ][i][CC][0][0] = tmp2 * fjac[k + 1][0][0]
						- tmp1 * njac[k + 1][0][0]
						- tmp1 * dz1;
					lhs_buf[k][j ][i][CC][1][0] = tmp2 * fjac[k + 1][1][0]
						- tmp1 * njac[k + 1][1][0];
					lhs_buf[k][j ][i][CC][2][0] = tmp2 * fjac[k + 1][2][0]
						- tmp1 * njac[k + 1][2][0];
					lhs_buf[k][j ][i][CC][3][0] = tmp2 * fjac[k + 1][3][0]
						- tmp1 * njac[k + 1][3][0];
					lhs_buf[k][j ][i][CC][4][0] = tmp2 * fjac[k + 1][4][0]
						- tmp1 * njac[k + 1][4][0];

					lhs_buf[k][j ][i][CC][0][1] = tmp2 * fjac[k + 1][0][1]
						- tmp1 * njac[k + 1][0][1];
					lhs_buf[k][j ][i][CC][1][1] = tmp2 * fjac[k + 1][1][1]
						- tmp1 * njac[k + 1][1][1]
						- tmp1 * dz2;
					lhs_buf[k][j ][i][CC][2][1] = tmp2 * fjac[k + 1][2][1]
						- tmp1 * njac[k + 1][2][1];
					lhs_buf[k][j ][i][CC][3][1] = tmp2 * fjac[k + 1][3][1]
						- tmp1 * njac[k + 1][3][1];
					lhs_buf[k][j ][i][CC][4][1] = tmp2 * fjac[k + 1][4][1]
						- tmp1 * njac[k + 1][4][1];

					lhs_buf[k][j ][i][CC][0][2] = tmp2 * fjac[k + 1][0][2]
						- tmp1 * njac[k + 1][0][2];
					lhs_buf[k][j ][i][CC][1][2] = tmp2 * fjac[k + 1][1][2]
						- tmp1 * njac[k + 1][1][2];
					lhs_buf[k][j ][i][CC][2][2] = tmp2 * fjac[k + 1][2][2]
						- tmp1 * njac[k + 1][2][2]
						- tmp1 * dz3;
					lhs_buf[k][j ][i][CC][3][2] = tmp2 * fjac[k + 1][3][2]
						- tmp1 * njac[k + 1][3][2];
					lhs_buf[k][j ][i][CC][4][2] = tmp2 * fjac[k + 1][4][2]
						- tmp1 * njac[k + 1][4][2];

					lhs_buf[k][j ][i][CC][0][3] = tmp2 * fjac[k + 1][0][3]
						- tmp1 * njac[k + 1][0][3];
					lhs_buf[k][j ][i][CC][1][3] = tmp2 * fjac[k + 1][1][3]
						- tmp1 * njac[k + 1][1][3];
					lhs_buf[k][j ][i][CC][2][3] = tmp2 * fjac[k + 1][2][3]
						- tmp1 * njac[k + 1][2][3];
					lhs_buf[k][j ][i][CC][3][3] = tmp2 * fjac[k + 1][3][3]
						- tmp1 * njac[k + 1][3][3]
						- tmp1 * dz4;
					lhs_buf[k][j ][i][CC][4][3] = tmp2 * fjac[k + 1][4][3]
						- tmp1 * njac[k + 1][4][3];

					lhs_buf[k][j ][i][CC][0][4] = tmp2 * fjac[k + 1][0][4]
						- tmp1 * njac[k + 1][0][4];
					lhs_buf[k][j ][i][CC][1][4] = tmp2 * fjac[k + 1][1][4]
						- tmp1 * njac[k + 1][1][4];
					lhs_buf[k][j ][i][CC][2][4] = tmp2 * fjac[k + 1][2][4]
						- tmp1 * njac[k + 1][2][4];
					lhs_buf[k][j ][i][CC][3][4] = tmp2 * fjac[k + 1][3][4]
						- tmp1 * njac[k + 1][3][4];
					lhs_buf[k][j ][i][CC][4][4] = tmp2 * fjac[k + 1][4][4]
						- tmp1 * njac[k + 1][4][4]
						- tmp1 * dz5;
				}
			
			}
		}
	}
			
	for (j = 1; j <= grid_points[1] - 2; j++) {
		for (i = 1; i <= grid_points[0] - 2; i++) {			
			for (k = 0; k <= ksize; k++)
			{
				if(k == 0)
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
				else if(k == ksize)
				{
					rhs[k][j][i][0] = rhs[k][j][i][0] - lhs_buf[k][j][i][AA][0][0]*rhs[k - 1][j][i][0]
									- lhs_buf[k][j][i][AA][1][0]*rhs[k - 1][j][i][1]
									- lhs_buf[k][j][i][AA][2][0]*rhs[k - 1][j][i][2]
									- lhs_buf[k][j][i][AA][3][0]*rhs[k - 1][j][i][3]
									- lhs_buf[k][j][i][AA][4][0]*rhs[k - 1][j][i][4];
					rhs[k][j][i][1] = rhs[k][j][i][1] - lhs_buf[k][j][i][AA][0][1]*rhs[k - 1][j][i][0]
									- lhs_buf[k][j][i][AA][1][1]*rhs[k - 1][j][i][1]
									- lhs_buf[k][j][i][AA][2][1]*rhs[k - 1][j][i][2]
									- lhs_buf[k][j][i][AA][3][1]*rhs[k - 1][j][i][3]
									- lhs_buf[k][j][i][AA][4][1]*rhs[k - 1][j][i][4];
					rhs[k][j][i][2] = rhs[k][j][i][2] - lhs_buf[k][j][i][AA][0][2]*rhs[k - 1][j][i][0]
									- lhs_buf[k][j][i][AA][1][2]*rhs[k - 1][j][i][1]
									- lhs_buf[k][j][i][AA][2][2]*rhs[k - 1][j][i][2]
									- lhs_buf[k][j][i][AA][3][2]*rhs[k - 1][j][i][3]
									- lhs_buf[k][j][i][AA][4][2]*rhs[k - 1][j][i][4];
					rhs[k][j][i][3] = rhs[k][j][i][3] - lhs_buf[k][j][i][AA][0][3]*rhs[k - 1][j][i][0]
									- lhs_buf[k][j][i][AA][1][3]*rhs[k - 1][j][i][1]
									- lhs_buf[k][j][i][AA][2][3]*rhs[k - 1][j][i][2]
									- lhs_buf[k][j][i][AA][3][3]*rhs[k - 1][j][i][3]
									- lhs_buf[k][j][i][AA][4][3]*rhs[k - 1][j][i][4];
					rhs[k][j][i][4] = rhs[k][j][i][4] - lhs_buf[k][j][i][AA][0][4]*rhs[k - 1][j][i][0]
									- lhs_buf[k][j][i][AA][1][4]*rhs[k - 1][j][i][1]
									- lhs_buf[k][j][i][AA][2][4]*rhs[k - 1][j][i][2]
									- lhs_buf[k][j][i][AA][3][4]*rhs[k - 1][j][i][3]
									- lhs_buf[k][j][i][AA][4][4]*rhs[k - 1][j][i][4];

									
									
									
					lhs_buf[k][j][i][BB][0][0] = lhs_buf[k][j][i][BB][0][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k - 1][j][i][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k - 1][j][i][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k - 1][j][i][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k - 1][j][i][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k - 1][j][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][1] = lhs_buf[k][j][i][BB][0][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k - 1][j][i][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k - 1][j][i][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k - 1][j][i][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k - 1][j][i][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k - 1][j][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][2] = lhs_buf[k][j][i][BB][0][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k - 1][j][i][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k - 1][j][i][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k - 1][j][i][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k - 1][j][i][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k - 1][j][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][3] = lhs_buf[k][j][i][BB][0][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k - 1][j][i][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k - 1][j][i][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k - 1][j][i][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k - 1][j][i][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k - 1][j][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][4] = lhs_buf[k][j][i][BB][0][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k - 1][j][i][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k - 1][j][i][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k - 1][j][i][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k - 1][j][i][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k - 1][j][i][CC][0][4];
					lhs_buf[k][j][i][BB][1][0] = lhs_buf[k][j][i][BB][1][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k - 1][j][i][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k - 1][j][i][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k - 1][j][i][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k - 1][j][i][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k - 1][j][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][1] = lhs_buf[k][j][i][BB][1][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k - 1][j][i][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k - 1][j][i][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k - 1][j][i][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k - 1][j][i][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k - 1][j][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][2] = lhs_buf[k][j][i][BB][1][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k - 1][j][i][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k - 1][j][i][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k - 1][j][i][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k - 1][j][i][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k - 1][j][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][3] = lhs_buf[k][j][i][BB][1][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k - 1][j][i][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k - 1][j][i][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k - 1][j][i][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k - 1][j][i][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k - 1][j][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][4] = lhs_buf[k][j][i][BB][1][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k - 1][j][i][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k - 1][j][i][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k - 1][j][i][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k - 1][j][i][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k - 1][j][i][CC][1][4];
					lhs_buf[k][j][i][BB][2][0] = lhs_buf[k][j][i][BB][2][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k - 1][j][i][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k - 1][j][i][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k - 1][j][i][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k - 1][j][i][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k - 1][j][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][1] = lhs_buf[k][j][i][BB][2][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k - 1][j][i][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k - 1][j][i][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k - 1][j][i][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k - 1][j][i][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k - 1][j][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][2] = lhs_buf[k][j][i][BB][2][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k - 1][j][i][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k - 1][j][i][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k - 1][j][i][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k - 1][j][i][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k - 1][j][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][3] = lhs_buf[k][j][i][BB][2][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k - 1][j][i][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k - 1][j][i][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k - 1][j][i][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k - 1][j][i][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k - 1][j][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][4] = lhs_buf[k][j][i][BB][2][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k - 1][j][i][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k - 1][j][i][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k - 1][j][i][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k - 1][j][i][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k - 1][j][i][CC][2][4];
					lhs_buf[k][j][i][BB][3][0] = lhs_buf[k][j][i][BB][3][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k - 1][j][i][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k - 1][j][i][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k - 1][j][i][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k - 1][j][i][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k - 1][j][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][1] = lhs_buf[k][j][i][BB][3][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k - 1][j][i][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k - 1][j][i][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k - 1][j][i][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k - 1][j][i][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k - 1][j][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][2] = lhs_buf[k][j][i][BB][3][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k - 1][j][i][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k - 1][j][i][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k - 1][j][i][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k - 1][j][i][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k - 1][j][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][3] = lhs_buf[k][j][i][BB][3][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k - 1][j][i][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k - 1][j][i][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k - 1][j][i][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k - 1][j][i][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k - 1][j][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][4] = lhs_buf[k][j][i][BB][3][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k - 1][j][i][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k - 1][j][i][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k - 1][j][i][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k - 1][j][i][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k - 1][j][i][CC][3][4];
					lhs_buf[k][j][i][BB][4][0] = lhs_buf[k][j][i][BB][4][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k - 1][j][i][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k - 1][j][i][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k - 1][j][i][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k - 1][j][i][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k - 1][j][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][1] = lhs_buf[k][j][i][BB][4][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k - 1][j][i][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k - 1][j][i][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k - 1][j][i][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k - 1][j][i][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k - 1][j][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][2] = lhs_buf[k][j][i][BB][4][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k - 1][j][i][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k - 1][j][i][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k - 1][j][i][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k - 1][j][i][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k - 1][j][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][3] = lhs_buf[k][j][i][BB][4][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k - 1][j][i][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k - 1][j][i][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k - 1][j][i][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k - 1][j][i][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k - 1][j][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][4] = lhs_buf[k][j][i][BB][4][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k - 1][j][i][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k - 1][j][i][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k - 1][j][i][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k - 1][j][i][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k - 1][j][i][CC][4][4];

					
					
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
					rhs[k][j][i][0] = rhs[k][j][i][0] - lhs_buf[k][j][i][AA][0][0]*rhs[k - 1][j][i][0]
									- lhs_buf[k][j][i][AA][1][0]*rhs[k - 1][j][i][1]
									- lhs_buf[k][j][i][AA][2][0]*rhs[k - 1][j][i][2]
									- lhs_buf[k][j][i][AA][3][0]*rhs[k - 1][j][i][3]
									- lhs_buf[k][j][i][AA][4][0]*rhs[k - 1][j][i][4];
					rhs[k][j][i][1] = rhs[k][j][i][1] - lhs_buf[k][j][i][AA][0][1]*rhs[k - 1][j][i][0]
									- lhs_buf[k][j][i][AA][1][1]*rhs[k - 1][j][i][1]
									- lhs_buf[k][j][i][AA][2][1]*rhs[k - 1][j][i][2]
									- lhs_buf[k][j][i][AA][3][1]*rhs[k - 1][j][i][3]
									- lhs_buf[k][j][i][AA][4][1]*rhs[k - 1][j][i][4];
					rhs[k][j][i][2] = rhs[k][j][i][2] - lhs_buf[k][j][i][AA][0][2]*rhs[k - 1][j][i][0]
									- lhs_buf[k][j][i][AA][1][2]*rhs[k - 1][j][i][1]
									- lhs_buf[k][j][i][AA][2][2]*rhs[k - 1][j][i][2]
									- lhs_buf[k][j][i][AA][3][2]*rhs[k - 1][j][i][3]
									- lhs_buf[k][j][i][AA][4][2]*rhs[k - 1][j][i][4];
					rhs[k][j][i][3] = rhs[k][j][i][3] - lhs_buf[k][j][i][AA][0][3]*rhs[k - 1][j][i][0]
									- lhs_buf[k][j][i][AA][1][3]*rhs[k - 1][j][i][1]
									- lhs_buf[k][j][i][AA][2][3]*rhs[k - 1][j][i][2]
									- lhs_buf[k][j][i][AA][3][3]*rhs[k - 1][j][i][3]
									- lhs_buf[k][j][i][AA][4][3]*rhs[k - 1][j][i][4];
					rhs[k][j][i][4] = rhs[k][j][i][4] - lhs_buf[k][j][i][AA][0][4]*rhs[k - 1][j][i][0]
									- lhs_buf[k][j][i][AA][1][4]*rhs[k - 1][j][i][1]
									- lhs_buf[k][j][i][AA][2][4]*rhs[k - 1][j][i][2]
									- lhs_buf[k][j][i][AA][3][4]*rhs[k - 1][j][i][3]
									- lhs_buf[k][j][i][AA][4][4]*rhs[k - 1][j][i][4];

									
									
									
									
					lhs_buf[k][j][i][BB][0][0] = lhs_buf[k][j][i][BB][0][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k - 1][j][i][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k - 1][j][i][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k - 1][j][i][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k - 1][j][i][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k - 1][j][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][1] = lhs_buf[k][j][i][BB][0][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k - 1][j][i][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k - 1][j][i][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k - 1][j][i][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k - 1][j][i][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k - 1][j][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][2] = lhs_buf[k][j][i][BB][0][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k - 1][j][i][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k - 1][j][i][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k - 1][j][i][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k - 1][j][i][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k - 1][j][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][3] = lhs_buf[k][j][i][BB][0][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k - 1][j][i][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k - 1][j][i][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k - 1][j][i][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k - 1][j][i][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k - 1][j][i][CC][0][4];
					lhs_buf[k][j][i][BB][0][4] = lhs_buf[k][j][i][BB][0][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k - 1][j][i][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k - 1][j][i][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k - 1][j][i][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k - 1][j][i][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k - 1][j][i][CC][0][4];
					lhs_buf[k][j][i][BB][1][0] = lhs_buf[k][j][i][BB][1][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k - 1][j][i][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k - 1][j][i][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k - 1][j][i][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k - 1][j][i][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k - 1][j][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][1] = lhs_buf[k][j][i][BB][1][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k - 1][j][i][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k - 1][j][i][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k - 1][j][i][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k - 1][j][i][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k - 1][j][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][2] = lhs_buf[k][j][i][BB][1][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k - 1][j][i][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k - 1][j][i][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k - 1][j][i][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k - 1][j][i][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k - 1][j][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][3] = lhs_buf[k][j][i][BB][1][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k - 1][j][i][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k - 1][j][i][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k - 1][j][i][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k - 1][j][i][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k - 1][j][i][CC][1][4];
					lhs_buf[k][j][i][BB][1][4] = lhs_buf[k][j][i][BB][1][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k - 1][j][i][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k - 1][j][i][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k - 1][j][i][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k - 1][j][i][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k - 1][j][i][CC][1][4];
					lhs_buf[k][j][i][BB][2][0] = lhs_buf[k][j][i][BB][2][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k - 1][j][i][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k - 1][j][i][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k - 1][j][i][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k - 1][j][i][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k - 1][j][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][1] = lhs_buf[k][j][i][BB][2][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k - 1][j][i][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k - 1][j][i][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k - 1][j][i][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k - 1][j][i][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k - 1][j][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][2] = lhs_buf[k][j][i][BB][2][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k - 1][j][i][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k - 1][j][i][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k - 1][j][i][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k - 1][j][i][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k - 1][j][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][3] = lhs_buf[k][j][i][BB][2][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k - 1][j][i][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k - 1][j][i][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k - 1][j][i][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k - 1][j][i][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k - 1][j][i][CC][2][4];
					lhs_buf[k][j][i][BB][2][4] = lhs_buf[k][j][i][BB][2][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k - 1][j][i][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k - 1][j][i][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k - 1][j][i][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k - 1][j][i][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k - 1][j][i][CC][2][4];
					lhs_buf[k][j][i][BB][3][0] = lhs_buf[k][j][i][BB][3][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k - 1][j][i][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k - 1][j][i][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k - 1][j][i][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k - 1][j][i][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k - 1][j][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][1] = lhs_buf[k][j][i][BB][3][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k - 1][j][i][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k - 1][j][i][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k - 1][j][i][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k - 1][j][i][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k - 1][j][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][2] = lhs_buf[k][j][i][BB][3][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k - 1][j][i][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k - 1][j][i][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k - 1][j][i][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k - 1][j][i][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k - 1][j][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][3] = lhs_buf[k][j][i][BB][3][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k - 1][j][i][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k - 1][j][i][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k - 1][j][i][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k - 1][j][i][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k - 1][j][i][CC][3][4];
					lhs_buf[k][j][i][BB][3][4] = lhs_buf[k][j][i][BB][3][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k - 1][j][i][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k - 1][j][i][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k - 1][j][i][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k - 1][j][i][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k - 1][j][i][CC][3][4];
					lhs_buf[k][j][i][BB][4][0] = lhs_buf[k][j][i][BB][4][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k - 1][j][i][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k - 1][j][i][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k - 1][j][i][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k - 1][j][i][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k - 1][j][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][1] = lhs_buf[k][j][i][BB][4][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k - 1][j][i][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k - 1][j][i][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k - 1][j][i][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k - 1][j][i][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k - 1][j][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][2] = lhs_buf[k][j][i][BB][4][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k - 1][j][i][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k - 1][j][i][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k - 1][j][i][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k - 1][j][i][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k - 1][j][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][3] = lhs_buf[k][j][i][BB][4][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k - 1][j][i][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k - 1][j][i][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k - 1][j][i][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k - 1][j][i][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k - 1][j][i][CC][4][4];
					lhs_buf[k][j][i][BB][4][4] = lhs_buf[k][j][i][BB][4][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k - 1][j][i][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k - 1][j][i][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k - 1][j][i][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k - 1][j][i][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k - 1][j][i][CC][4][4];

					
					
					
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

	
	for (j = 1; j <= grid_points[1] - 2; j++) {
		for (i = 1; i <= grid_points[0] - 2; i++) {
			for (k = ksize - 1; k >= 0; k--) {
				for (m = 0; m < BLOCK_SIZE; m++) {
					for (n = 0; n < BLOCK_SIZE; n++) {
						rhs[k][j][i][m] = rhs[k][j][i][m]
							- lhs_buf[k][j][i][CC][n][m] * rhs[k + 1][j][i][n];
					}
				}
			}
		}
	}
	
	
	if (timeron) timer_stop(t_zsolve);
	

	}







