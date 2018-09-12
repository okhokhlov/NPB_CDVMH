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

double lhs_[2][3][5][5];
double rhs_[2][5];
int len = 4;
//double fjac_[4][5][5];
//double njac_[4][5][5];

void matvec_sub_(int a1, int a2, int b3, int c3);
void matmul_sub_(int a1, int a2, int b1, int b2, int c1, int c2);
void binvcrhs_(int a1, int a2, int b1, int b2, int c3);
void binvrhs_(int a1, int a2, int b3);
void buf_compute_(int i, int a);
void fn_init_(int a, int i, int j, int k);
//---------------------------------------------------------------------
// 
// Performs line solves in X direction by first factoring
// the block-tridiagonal matrix into an upper triangular matrix, 
// and then performing back substitution to solve for the unknow
// vectors of each line.  
// 
// Make sure we treat elements zero to cell_size in the direction
// of the sweep.
// 
//---------------------------------------------------------------------
void x_solve()
{
	int i, j, k, m, n, isize, l = 0;
	
	double t1, t2, t3, tm1, tm2, tm3;
	double pivot, coeff;
	//int z, a = -1, b = 2;
	//---------------------------------------------------------------------
	//---------------------------------------------------------------------

	if (timeron) timer_start(t_xsolve);

	//---------------------------------------------------------------------
	//---------------------------------------------------------------------

	//---------------------------------------------------------------------
	// This function computes the left hand side in the xi-direction
	//---------------------------------------------------------------------

	isize = grid_points[0] - 1;

	//---------------------------------------------------------------------
	// determine a (labeled f) and n jacobians
	//---------------------------------------------------------------------


//#pragma dvm redistribute (lhs_buf[block][block][block][][][])	
	
//#pragma dvm region out(lhs_buf), in(u, qs, square)
{	
////#pragma dvm parallel ([k][j][i] on lhs_buf[k][j][i][][][]) private(n,m,tm1, tm2, tm3, t1, t2, t3 , tmp1, tmp2, tmp3), shadow_renew (u[0:0][0:0][2:2][0:0], qs[0:0][0:0][1:1], square[0:0][0:0][1:1])
	for (k = 1; k <= grid_points[2] - 2; k++) {
		for (j = 1; j <= grid_points[1] - 2; j++) {
			for (i = 0; i <= isize; i++) {
				if(i == 0)
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
				else if(i == isize)
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
					tmp1 = 1.0 / u[k][j][i][0]; //rho_i[k][j][i];
					tmp2 = tmp1 * tmp1;
					tmp3 = tmp1 * tmp2;

					tm1 = 1.0 / u[k][j][i + 1][0];  //rho_i[k][j][i + 1];
					tm2 = tm1 * tm1;
					tm3 = tm1 * tm2;
					
					t1 = 1.0 / u[k][j][i - 1][0];   //rho_i[k][j][i - 1];
					t2 = t1 * t1;
					t3 = t1 * t2;
					
					
					lhs_buf[k][j][i][AA][0][0] = - dt * tx1 * dx1;
					lhs_buf[k][j][i][AA][1][0] = -dt * tx2;
					lhs_buf[k][j][i][AA][2][0] = 0;
					lhs_buf[k][j][i][AA][3][0] = 0;
					lhs_buf[k][j][i][AA][4][0] = 0;

						
					lhs_buf[k][j][i][AA][0][1] = -dt * tx2 * (-(u[k][j][i - 1][1] * t2 * u[k][j][i - 1][1])
						+ c2 * qs[k][j][i - 1]) - dt * tx1 * (-con43 * c3c4 * t2 * u[k][j][i - 1][1]);
					lhs_buf[k][j][i][AA][1][1] = -dt * tx2 * (2.0 - c2) * (u[k][j][i - 1][1] / u[k][j][i - 1][0])
						- dt * tx1 * con43 * c3c4 * t1
						- dt * tx1 * dx2;
					lhs_buf[k][j][i][AA][2][1] = -dt * tx2 * (-c2 * (u[k][j][i - 1][2] * t1));
					lhs_buf[k][j][i][AA][3][1] = -dt * tx2 * (-c2 * (u[k][j][i - 1][3] * t1));
					lhs_buf[k][j][i][AA][4][1] = -dt * tx2 * c2;

						
					lhs_buf[k][j][i][AA][0][2] = -dt * tx2 * (-(u[k][j][i - 1][1] * u[k][j][i - 1][2]) * t2)
						- dt * tx1 * (-c3c4 * t2 * u[k][j][i - 1][2]);
					lhs_buf[k][j][i][AA][1][2] = -dt * tx2 * u[k][j][i - 1][2] * t1;
					lhs_buf[k][j][i][AA][2][2] = -dt * tx2 * u[k][j][i - 1][1] * t1
						- dt * tx1 * c3c4 * t1
						- dt * tx1 * dx3;
					lhs_buf[k][j][i][AA][3][2] = 0;
					lhs_buf[k][j][i][AA][4][2] = 0;

						
					lhs_buf[k][j][i][AA][0][3] = -dt * tx2 * (-(u[k][j][i - 1][1] * u[k][j][i - 1][3]) * t2)
						- dt * tx1 * (-c3c4 * t2 * u[k][j][i - 1][3]);
					lhs_buf[k][j][i][AA][1][3] = -dt * tx2 * u[k][j][i - 1][3] * t1;
					lhs_buf[k][j][i][AA][2][3] = 0;
					lhs_buf[k][j][i][AA][3][3] = -dt * tx2 * u[k][j][i - 1][1] * t1
						- dt * tx1 * c3c4 * t1
						- dt * tx1 * dx4;
					lhs_buf[k][j][i][AA][4][3] = 0;

						
					lhs_buf[k][j][i][AA][0][4] = -dt * tx2 * (c2 * 2.0 * square[k][j][i - 1] - c1 * u[k][j][i - 1][4])
						* (u[k][j][i - 1][1] * t2)
						- dt * tx1 * (-(con43 * c3c4
						- c1345) * t3 * (u[k][j][i - 1][1] * u[k][j][i - 1][1])
						- (c3c4 - c1345) * t3 * (u[k][j][i - 1][2] * u[k][j][i - 1][2])
						- (c3c4 - c1345) * t3 * (u[k][j][i - 1][3] * u[k][j][i - 1][3])
						- c1345 * t2 * u[k][j][i - 1][4]);
					lhs_buf[k][j][i][AA][1][4] = -dt * tx2 * (c1 *  u[k][j][i - 1][4] * t1
						- c2 * (u[k][j][i - 1][1] * u[k][j][i - 1][1] * t2 + qs[k][j][i - 1]))
						- dt * tx1 * (con43 * c3c4 - c1345) * t2 * u[k][j][i - 1][1];
					lhs_buf[k][j][i][AA][2][4] = -dt * tx2 * (-c2 * (u[k][j][i - 1][2] * u[k][j][i - 1][1]) * t2)
						- dt * tx1 * (c3c4 - c1345) * t2 * u[k][j][i - 1][2];
					lhs_buf[k][j][i][AA][3][4] = -dt * tx2 * (-c2 * (u[k][j][i - 1][3] * u[k][j][i - 1][1]) * t2)
						- dt * tx1 * (c3c4 - c1345) * t2 * u[k][j][i - 1][3];
					lhs_buf[k][j][i][AA][4][4] = -dt * tx2 * c1 * (u[k][j][i - 1][1] * t1)
						- dt * tx1 * (c1345)* t1
						- dt * tx1 * dx5;

						
						
						
						
						
					lhs_buf[k][j][i][BB][0][0] = 1.0
						+ dt * tx1 * 2.0 * dx1;
					lhs_buf[k][j][i][BB][1][0] = 0;
					lhs_buf[k][j][i][BB][2][0] = 0;
					lhs_buf[k][j][i][BB][3][0] = 0;
					lhs_buf[k][j][i][BB][4][0] = 0;

					
					lhs_buf[k][j][i][BB][0][1] = dt * tx1 * 2.0 * (-con43 * c3c4 * tmp2 * u[k][j][i][1]);
					lhs_buf[k][j][i][BB][1][1] = 1.0 + dt * tx1 * 2.0 * con43 * c3c4 * tmp1
						+ dt * tx1 * 2.0 * dx2;
					lhs_buf[k][j][i][BB][2][1] = 0;
					lhs_buf[k][j][i][BB][3][1] = 0;
					lhs_buf[k][j][i][BB][4][1] = 0;

					
					lhs_buf[k][j][i][BB][0][2] = dt * tx1 * 2.0 * (-c3c4 * tmp2 * u[k][j][i][2]);
					lhs_buf[k][j][i][BB][1][2] = 0;
					lhs_buf[k][j][i][BB][2][2] = 1.0 + dt * tx1 * 2.0 * c3c4 * tmp1
						+ dt * tx1 * 2.0 * dx3;
					lhs_buf[k][j][i][BB][3][2] = 0;
					lhs_buf[k][j][i][BB][4][2] = 0;

					
					lhs_buf[k][j][i][BB][0][3] = dt * tx1 * 2.0 * (-c3c4 * tmp2 * u[k][j][i][3]);
					lhs_buf[k][j][i][BB][1][3] = 0;
					lhs_buf[k][j][i][BB][2][3] = 0;
					lhs_buf[k][j][i][BB][3][3] = 1.0 + dt * tx1 * 2.0 * c3c4 * tmp1
						+ dt * tx1 * 2.0 * dx4;
					lhs_buf[k][j][i][BB][4][3] = 0;

					
					lhs_buf[k][j][i][BB][0][4] = dt * tx1 * 2.0 * (-(con43 * c3c4
						- c1345) * tmp3 * (u[k][j][i][1] * u[k][j][i][1])
						- (c3c4 - c1345) * tmp3 * (u[k][j][i][2] * u[k][j][i][2])
						- (c3c4 - c1345) * tmp3 * (u[k][j][i][3] * u[k][j][i][3])
						- c1345 * tmp2 * u[k][j][i][4]);
					lhs_buf[k][j][i][BB][1][4] = dt * tx1 * 2.0 * (con43 * c3c4	- c1345) * tmp2 * u[k][j][i][1];
					lhs_buf[k][j][i][BB][2][4] = dt * tx1 * 2.0 * (c3c4 - c1345) * tmp2 * u[k][j][i][2];
					lhs_buf[k][j][i][BB][3][4] = dt * tx1 * 2.0 * (c3c4 - c1345) * tmp2 * u[k][j][i][3];
					lhs_buf[k][j][i][BB][4][4] = 1.0
						+ dt * tx1 * 2.0 * (c1345)* tmp1
						+ dt * tx1 * 2.0 * dx5;

						
						
						
						
					lhs_buf[k][j][i][CC][0][0] = - dt * tx1 * dx1;
					lhs_buf[k][j][i][CC][1][0] = dt * tx2;
					lhs_buf[k][j][i][CC][2][0] = 0;
					lhs_buf[k][j][i][CC][3][0] = 0;
					lhs_buf[k][j][i][CC][4][0] = 0;

					
					lhs_buf[k][j][i][CC][0][1] = dt * tx2 * (-(u[k][j][i + 1][1] * tm2 * u[k][j][i + 1][1])
						+ c2 * qs[k][j][i + 1])
						- dt * tx1 * (-con43 * c3c4 * tm2 * u[k][j][i + 1][1]);
					lhs_buf[k][j][i][CC][1][1] = dt * tx2 * ((2.0 - c2) * (u[k][j][i + 1][1] / u[k][j][i + 1][0]))
						- dt * tx1 * con43 * c3c4 * tm1
						- dt * tx1 * dx2;
					lhs_buf[k][j][i][CC][2][1] = dt * tx2 * (-c2 * (u[k][j][i + 1][2] * tm1));
					lhs_buf[k][j][i][CC][3][1] = dt * tx2 * (-c2 * (u[k][j][i + 1][3] * tm1));
					lhs_buf[k][j][i][CC][4][1] = dt * tx2 * c2;

					
					lhs_buf[k][j][i][CC][0][2] = dt * tx2 * (-(u[k][j][i + 1][1] * u[k][j][i + 1][2]) * tm2)
						- dt * tx1 * (-c3c4 * tm2 * u[k][j][i + 1][2]);
					lhs_buf[k][j][i][CC][1][2] = dt * tx2 * u[k][j][i + 1][2] * tm1;
					lhs_buf[k][j][i][CC][2][2] = dt * tx2 * u[k][j][i + 1][1] * tm1
						- dt * tx1 * c3c4 * tm1
						- dt * tx1 * dx3;
					lhs_buf[k][j][i][CC][3][2] = 0;
					lhs_buf[k][j][i][CC][4][2] = 0;

					
					lhs_buf[k][j][i][CC][0][3] = dt * tx2 * (-(u[k][j][i + 1][1] * u[k][j][i + 1][3]) * tm2)
						- dt * tx1 * (-c3c4 * tm2 * u[k][j][i + 1][3]);
					lhs_buf[k][j][i][CC][1][3] = dt * tx2 * u[k][j][i + 1][3] * tm1;
					lhs_buf[k][j][i][CC][2][3] = 0;
					lhs_buf[k][j][i][CC][3][3] = dt * tx2 * u[k][j][i + 1][1] * tm1
						- dt * tx1 * c3c4 * tm1
						- dt * tx1 * dx4;
					lhs_buf[k][j][i][CC][4][3] = 0;

					
					lhs_buf[k][j][i][CC][0][4] = dt * tx2 * (c2 * 2.0 * square[k][j][i + 1] - c1 * u[k][j][i + 1][4])
						* (u[k][j][i + 1][1] * tm2)
						- dt * tx1 * (-(con43 * c3c4
						- c1345) * tm3 * (u[k][j][i + 1][1] * u[k][j][i + 1][1])
						- (c3c4 - c1345) * tm3 * (u[k][j][i + 1][2] * u[k][j][i + 1][2])
						- (c3c4 - c1345) * tm3 * (u[k][j][i + 1][3] * u[k][j][i + 1][3])
						- c1345 * tm2 * u[k][j][i + 1][4]);
					lhs_buf[k][j][i][CC][1][4] = dt * tx2 * (c1 *  u[k][j][i + 1][4] * tm1
						- c2 * (u[k][j][i + 1][1] * u[k][j][i + 1][1] * tm2 + qs[k][j][i + 1]))
						- dt * tx1 * (con43 * c3c4 - c1345) * tm2 * u[k][j][i + 1][1];
					lhs_buf[k][j][i][CC][2][4] = dt * tx2 * (-c2 * (u[k][j][i + 1][2] * u[k][j][i + 1][1]) * tm2)
						- dt * tx1 * (c3c4 - c1345) * tm2 * u[k][j][i + 1][2];
					lhs_buf[k][j][i][CC][3][4] = dt * tx2 * (-c2 * (u[k][j][i + 1][3] * u[k][j][i + 1][1]) * tm2)
						- dt * tx1 * (c3c4 - c1345) * tm2 * u[k][j][i + 1][3];
					lhs_buf[k][j][i][CC][4][4] = dt * tx2 * c1 * (u[k][j][i + 1][1] * tm1)
						- dt * tx1 * (c1345)* tm1
						- dt * tx1 * dx5;
				}
			}

		}
	}	

}//end region
	
//#pragma dvm redistribute (lhs_buf[][][][][][])	
	
//#pragma dvm redistribute (lhs_buf[block][block][block][][][])

//#pragma dvm parallel([k][j][i] on lhs_buf[k][j][i][][][]) //private(n,m, pivot , coeff)//, across(rhs[0:0][0:0][1:0][0:0], lhs_buf[0:0][0:0][1:0][0:0][0:0][0:0])
	for (k = 1; k <= grid_points[2] - 2; k++) {
		for (j = 1; j <= grid_points[1] - 2; j++) {
			for (i = 0; i <= isize - 1; i++) {
				if(i == 0)
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
				else if (i == isize)
				{
					rhs[k][j][i][0] = rhs[k][j][i][0] - lhs_buf[k][j][i][AA][0][0]*rhs[k][j][i - 1][0]
									- lhs_buf[k][j][i][AA][1][0]*rhs[k][j][i - 1][1]
									- lhs_buf[k][j][i][AA][2][0]*rhs[k][j][i - 1][2]
									- lhs_buf[k][j][i][AA][3][0]*rhs[k][j][i - 1][3]
									- lhs_buf[k][j][i][AA][4][0]*rhs[k][j][i - 1][4];
					rhs[k][j][i][1] = rhs[k][j][i][1] - lhs_buf[k][j][i][AA][0][1]*rhs[k][j][i - 1][0]
									- lhs_buf[k][j][i][AA][1][1]*rhs[k][j][i - 1][1]
									- lhs_buf[k][j][i][AA][2][1]*rhs[k][j][i - 1][2]
									- lhs_buf[k][j][i][AA][3][1]*rhs[k][j][i - 1][3]
									- lhs_buf[k][j][i][AA][4][1]*rhs[k][j][i - 1][4];
					rhs[k][j][i][2] = rhs[k][j][i][2] - lhs_buf[k][j][i][AA][0][2]*rhs[k][j][i - 1][0]
									- lhs_buf[k][j][i][AA][1][2]*rhs[k][j][i - 1][1]
									- lhs_buf[k][j][i][AA][2][2]*rhs[k][j][i - 1][2]
									- lhs_buf[k][j][i][AA][3][2]*rhs[k][j][i - 1][3]
									- lhs_buf[k][j][i][AA][4][2]*rhs[k][j][i - 1][4];
					rhs[k][j][i][3] = rhs[k][j][i][3] - lhs_buf[k][j][i][AA][0][3]*rhs[k][j][i - 1][0]
									- lhs_buf[k][j][i][AA][1][3]*rhs[k][j][i - 1][1]
									- lhs_buf[k][j][i][AA][2][3]*rhs[k][j][i - 1][2]
									- lhs_buf[k][j][i][AA][3][3]*rhs[k][j][i - 1][3]
									- lhs_buf[k][j][i][AA][4][3]*rhs[k][j][i - 1][4];
					rhs[k][j][i][4] = rhs[k][j][i][4] - lhs_buf[k][j][i][AA][0][4]*rhs[k][j][i - 1][0]
									- lhs_buf[k][j][i][AA][1][4]*rhs[k][j][i - 1][1]
									- lhs_buf[k][j][i][AA][2][4]*rhs[k][j][i - 1][2]
									- lhs_buf[k][j][i][AA][3][4]*rhs[k][j][i - 1][3]
									- lhs_buf[k][j][i][AA][4][4]*rhs[k][j][i - 1][4];
					


					
					lhs_buf[k][j][i][BB][0][0] = lhs_buf[k][j][i][BB][0][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j][i - 1][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j][i - 1][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j][i - 1][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j][i - 1][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j][i - 1][CC][0][4];
					lhs_buf[k][j][i][BB][0][1] = lhs_buf[k][j][i][BB][0][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j][i - 1][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j][i - 1][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j][i - 1][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j][i - 1][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j][i - 1][CC][0][4];
					lhs_buf[k][j][i][BB][0][2] = lhs_buf[k][j][i][BB][0][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j][i - 1][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j][i - 1][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j][i - 1][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j][i - 1][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j][i - 1][CC][0][4];
					lhs_buf[k][j][i][BB][0][3] = lhs_buf[k][j][i][BB][0][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j][i - 1][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j][i - 1][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j][i - 1][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j][i - 1][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j][i - 1][CC][0][4];
					lhs_buf[k][j][i][BB][0][4] = lhs_buf[k][j][i][BB][0][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j][i - 1][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j][i - 1][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j][i - 1][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j][i - 1][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j][i - 1][CC][0][4];
					lhs_buf[k][j][i][BB][1][0] = lhs_buf[k][j][i][BB][1][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j][i - 1][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j][i - 1][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j][i - 1][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j][i - 1][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j][i - 1][CC][1][4];
					lhs_buf[k][j][i][BB][1][1] = lhs_buf[k][j][i][BB][1][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j][i - 1][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j][i - 1][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j][i - 1][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j][i - 1][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j][i - 1][CC][1][4];
					lhs_buf[k][j][i][BB][1][2] = lhs_buf[k][j][i][BB][1][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j][i - 1][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j][i - 1][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j][i - 1][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j][i - 1][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j][i - 1][CC][1][4];
					lhs_buf[k][j][i][BB][1][3] = lhs_buf[k][j][i][BB][1][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j][i - 1][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j][i - 1][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j][i - 1][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j][i - 1][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j][i - 1][CC][1][4];
					lhs_buf[k][j][i][BB][1][4] = lhs_buf[k][j][i][BB][1][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j][i - 1][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j][i - 1][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j][i - 1][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j][i - 1][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j][i - 1][CC][1][4];
					lhs_buf[k][j][i][BB][2][0] = lhs_buf[k][j][i][BB][2][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j][i - 1][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j][i - 1][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j][i - 1][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j][i - 1][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j][i - 1][CC][2][4];
					lhs_buf[k][j][i][BB][2][1] = lhs_buf[k][j][i][BB][2][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j][i - 1][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j][i - 1][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j][i - 1][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j][i - 1][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j][i - 1][CC][2][4];
					lhs_buf[k][j][i][BB][2][2] = lhs_buf[k][j][i][BB][2][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j][i - 1][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j][i - 1][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j][i - 1][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j][i - 1][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j][i - 1][CC][2][4];
					lhs_buf[k][j][i][BB][2][3] = lhs_buf[k][j][i][BB][2][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j][i - 1][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j][i - 1][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j][i - 1][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j][i - 1][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j][i - 1][CC][2][4];
					lhs_buf[k][j][i][BB][2][4] = lhs_buf[k][j][i][BB][2][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j][i - 1][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j][i - 1][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j][i - 1][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j][i - 1][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j][i - 1][CC][2][4];
					lhs_buf[k][j][i][BB][3][0] = lhs_buf[k][j][i][BB][3][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j][i - 1][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j][i - 1][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j][i - 1][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j][i - 1][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j][i - 1][CC][3][4];
					lhs_buf[k][j][i][BB][3][1] = lhs_buf[k][j][i][BB][3][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j][i - 1][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j][i - 1][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j][i - 1][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j][i - 1][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j][i - 1][CC][3][4];
					lhs_buf[k][j][i][BB][3][2] = lhs_buf[k][j][i][BB][3][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j][i - 1][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j][i - 1][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j][i - 1][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j][i - 1][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j][i - 1][CC][3][4];
					lhs_buf[k][j][i][BB][3][3] = lhs_buf[k][j][i][BB][3][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j][i - 1][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j][i - 1][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j][i - 1][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j][i - 1][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j][i - 1][CC][3][4];
					lhs_buf[k][j][i][BB][3][4] = lhs_buf[k][j][i][BB][3][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j][i - 1][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j][i - 1][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j][i - 1][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j][i - 1][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j][i - 1][CC][3][4];
					lhs_buf[k][j][i][BB][4][0] = lhs_buf[k][j][i][BB][4][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j][i - 1][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j][i - 1][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j][i - 1][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j][i - 1][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j][i - 1][CC][4][4];
					lhs_buf[k][j][i][BB][4][1] = lhs_buf[k][j][i][BB][4][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j][i - 1][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j][i - 1][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j][i - 1][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j][i - 1][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j][i - 1][CC][4][4];
					lhs_buf[k][j][i][BB][4][2] = lhs_buf[k][j][i][BB][4][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j][i - 1][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j][i - 1][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j][i - 1][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j][i - 1][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j][i - 1][CC][4][4];
					lhs_buf[k][j][i][BB][4][3] = lhs_buf[k][j][i][BB][4][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j][i - 1][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j][i - 1][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j][i - 1][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j][i - 1][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j][i - 1][CC][4][4];
					lhs_buf[k][j][i][BB][4][4] = lhs_buf[k][j][i][BB][4][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j][i - 1][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j][i - 1][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j][i - 1][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j][i - 1][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j][i - 1][CC][4][4];					
					
					
					


					
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
					rhs[k][j][i][0] = rhs[k][j][i][0] - lhs_buf[k][j][i][AA][0][0]*rhs[k][j][i - 1][0]
									- lhs_buf[k][j][i][AA][1][0]*rhs[k][j][i - 1][1]
									- lhs_buf[k][j][i][AA][2][0]*rhs[k][j][i - 1][2]
									- lhs_buf[k][j][i][AA][3][0]*rhs[k][j][i - 1][3]
									- lhs_buf[k][j][i][AA][4][0]*rhs[k][j][i - 1][4];
					rhs[k][j][i][1] = rhs[k][j][i][1] - lhs_buf[k][j][i][AA][0][1]*rhs[k][j][i - 1][0]
									- lhs_buf[k][j][i][AA][1][1]*rhs[k][j][i - 1][1]
									- lhs_buf[k][j][i][AA][2][1]*rhs[k][j][i - 1][2]
									- lhs_buf[k][j][i][AA][3][1]*rhs[k][j][i - 1][3]
									- lhs_buf[k][j][i][AA][4][1]*rhs[k][j][i - 1][4];
					rhs[k][j][i][2] = rhs[k][j][i][2] - lhs_buf[k][j][i][AA][0][2]*rhs[k][j][i - 1][0]
									- lhs_buf[k][j][i][AA][1][2]*rhs[k][j][i - 1][1]
									- lhs_buf[k][j][i][AA][2][2]*rhs[k][j][i - 1][2]
									- lhs_buf[k][j][i][AA][3][2]*rhs[k][j][i - 1][3]
									- lhs_buf[k][j][i][AA][4][2]*rhs[k][j][i - 1][4];
					rhs[k][j][i][3] = rhs[k][j][i][3] - lhs_buf[k][j][i][AA][0][3]*rhs[k][j][i - 1][0]
									- lhs_buf[k][j][i][AA][1][3]*rhs[k][j][i - 1][1]
									- lhs_buf[k][j][i][AA][2][3]*rhs[k][j][i - 1][2]
									- lhs_buf[k][j][i][AA][3][3]*rhs[k][j][i - 1][3]
									- lhs_buf[k][j][i][AA][4][3]*rhs[k][j][i - 1][4];
					rhs[k][j][i][4] = rhs[k][j][i][4] - lhs_buf[k][j][i][AA][0][4]*rhs[k][j][i - 1][0]
									- lhs_buf[k][j][i][AA][1][4]*rhs[k][j][i - 1][1]
									- lhs_buf[k][j][i][AA][2][4]*rhs[k][j][i - 1][2]
									- lhs_buf[k][j][i][AA][3][4]*rhs[k][j][i - 1][3]
									- lhs_buf[k][j][i][AA][4][4]*rhs[k][j][i - 1][4];
					

					
					lhs_buf[k][j][i][BB][0][0] = lhs_buf[k][j][i][BB][0][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j][i - 1][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j][i - 1][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j][i - 1][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j][i - 1][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j][i - 1][CC][0][4];
					lhs_buf[k][j][i][BB][0][1] = lhs_buf[k][j][i][BB][0][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j][i - 1][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j][i - 1][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j][i - 1][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j][i - 1][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j][i - 1][CC][0][4];
					lhs_buf[k][j][i][BB][0][2] = lhs_buf[k][j][i][BB][0][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j][i - 1][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j][i - 1][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j][i - 1][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j][i - 1][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j][i - 1][CC][0][4];
					lhs_buf[k][j][i][BB][0][3] = lhs_buf[k][j][i][BB][0][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j][i - 1][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j][i - 1][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j][i - 1][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j][i - 1][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j][i - 1][CC][0][4];
					lhs_buf[k][j][i][BB][0][4] = lhs_buf[k][j][i][BB][0][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j][i - 1][CC][0][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j][i - 1][CC][0][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j][i - 1][CC][0][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j][i - 1][CC][0][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j][i - 1][CC][0][4];
					lhs_buf[k][j][i][BB][1][0] = lhs_buf[k][j][i][BB][1][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j][i - 1][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j][i - 1][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j][i - 1][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j][i - 1][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j][i - 1][CC][1][4];
					lhs_buf[k][j][i][BB][1][1] = lhs_buf[k][j][i][BB][1][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j][i - 1][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j][i - 1][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j][i - 1][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j][i - 1][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j][i - 1][CC][1][4];
					lhs_buf[k][j][i][BB][1][2] = lhs_buf[k][j][i][BB][1][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j][i - 1][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j][i - 1][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j][i - 1][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j][i - 1][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j][i - 1][CC][1][4];
					lhs_buf[k][j][i][BB][1][3] = lhs_buf[k][j][i][BB][1][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j][i - 1][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j][i - 1][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j][i - 1][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j][i - 1][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j][i - 1][CC][1][4];
					lhs_buf[k][j][i][BB][1][4] = lhs_buf[k][j][i][BB][1][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j][i - 1][CC][1][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j][i - 1][CC][1][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j][i - 1][CC][1][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j][i - 1][CC][1][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j][i - 1][CC][1][4];
					lhs_buf[k][j][i][BB][2][0] = lhs_buf[k][j][i][BB][2][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j][i - 1][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j][i - 1][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j][i - 1][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j][i - 1][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j][i - 1][CC][2][4];
					lhs_buf[k][j][i][BB][2][1] = lhs_buf[k][j][i][BB][2][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j][i - 1][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j][i - 1][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j][i - 1][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j][i - 1][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j][i - 1][CC][2][4];
					lhs_buf[k][j][i][BB][2][2] = lhs_buf[k][j][i][BB][2][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j][i - 1][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j][i - 1][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j][i - 1][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j][i - 1][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j][i - 1][CC][2][4];
					lhs_buf[k][j][i][BB][2][3] = lhs_buf[k][j][i][BB][2][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j][i - 1][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j][i - 1][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j][i - 1][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j][i - 1][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j][i - 1][CC][2][4];
					lhs_buf[k][j][i][BB][2][4] = lhs_buf[k][j][i][BB][2][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j][i - 1][CC][2][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j][i - 1][CC][2][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j][i - 1][CC][2][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j][i - 1][CC][2][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j][i - 1][CC][2][4];
					lhs_buf[k][j][i][BB][3][0] = lhs_buf[k][j][i][BB][3][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j][i - 1][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j][i - 1][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j][i - 1][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j][i - 1][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j][i - 1][CC][3][4];
					lhs_buf[k][j][i][BB][3][1] = lhs_buf[k][j][i][BB][3][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j][i - 1][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j][i - 1][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j][i - 1][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j][i - 1][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j][i - 1][CC][3][4];
					lhs_buf[k][j][i][BB][3][2] = lhs_buf[k][j][i][BB][3][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j][i - 1][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j][i - 1][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j][i - 1][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j][i - 1][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j][i - 1][CC][3][4];
					lhs_buf[k][j][i][BB][3][3] = lhs_buf[k][j][i][BB][3][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j][i - 1][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j][i - 1][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j][i - 1][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j][i - 1][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j][i - 1][CC][3][4];
					lhs_buf[k][j][i][BB][3][4] = lhs_buf[k][j][i][BB][3][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j][i - 1][CC][3][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j][i - 1][CC][3][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j][i - 1][CC][3][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j][i - 1][CC][3][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j][i - 1][CC][3][4];
					lhs_buf[k][j][i][BB][4][0] = lhs_buf[k][j][i][BB][4][0] - lhs_buf[k][j][i][AA][0][0]*lhs_buf[k][j][i - 1][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][0]*lhs_buf[k][j][i - 1][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][0]*lhs_buf[k][j][i - 1][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][0]*lhs_buf[k][j][i - 1][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][0]*lhs_buf[k][j][i - 1][CC][4][4];
					lhs_buf[k][j][i][BB][4][1] = lhs_buf[k][j][i][BB][4][1] - lhs_buf[k][j][i][AA][0][1]*lhs_buf[k][j][i - 1][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][1]*lhs_buf[k][j][i - 1][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][1]*lhs_buf[k][j][i - 1][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][1]*lhs_buf[k][j][i - 1][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][1]*lhs_buf[k][j][i - 1][CC][4][4];
					lhs_buf[k][j][i][BB][4][2] = lhs_buf[k][j][i][BB][4][2] - lhs_buf[k][j][i][AA][0][2]*lhs_buf[k][j][i - 1][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][2]*lhs_buf[k][j][i - 1][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][2]*lhs_buf[k][j][i - 1][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][2]*lhs_buf[k][j][i - 1][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][2]*lhs_buf[k][j][i - 1][CC][4][4];
					lhs_buf[k][j][i][BB][4][3] = lhs_buf[k][j][i][BB][4][3] - lhs_buf[k][j][i][AA][0][3]*lhs_buf[k][j][i - 1][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][3]*lhs_buf[k][j][i - 1][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][3]*lhs_buf[k][j][i - 1][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][3]*lhs_buf[k][j][i - 1][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][3]*lhs_buf[k][j][i - 1][CC][4][4];
					lhs_buf[k][j][i][BB][4][4] = lhs_buf[k][j][i][BB][4][4] - lhs_buf[k][j][i][AA][0][4]*lhs_buf[k][j][i - 1][CC][4][0]
											  - lhs_buf[k][j][i][AA][1][4]*lhs_buf[k][j][i - 1][CC][4][1]
											  - lhs_buf[k][j][i][AA][2][4]*lhs_buf[k][j][i - 1][CC][4][2]
											  - lhs_buf[k][j][i][AA][3][4]*lhs_buf[k][j][i - 1][CC][4][3]
											  - lhs_buf[k][j][i][AA][4][4]*lhs_buf[k][j][i - 1][CC][4][4];

					

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

//#pragma dvm redistribute (lhs_buf[][][][][][])	
	
	
	double rhs_[BLOCK_SIZE], rhs__[BLOCK_SIZE];

//#pragma dvm redistribute (lhs_buf[block][block][block][][][])

	//#pragma dvm region local(lhs_buf)
	//{
//#pragma dvm parallel([k][j][i] on rhs[k][j][i][]) //private(n,m, rhs_, rhs__), across(rhs[0:0][0:0][0:1][0:0])//, shadow_renew (rhs[0:1][0:1][0:1][0:0])//, stage(isize)
		for (k = 1; k <= grid_points[2] - 2; k++) {
			for (j = 1; j <= grid_points[1] - 2; j++) {
				for (i = isize - 1; i >= 0; i--) {
					for (m = 0; m < BLOCK_SIZE; m++) {
						for (n = 0; n < BLOCK_SIZE; n++) {
							rhs[k][j][i][m] = rhs[k][j][i][m]
							- lhs_buf[k][j][i][CC][n][m]*rhs[k][j][i+1][n];
						}
					}
					//for (m = 0; m < BLOCK_SIZE; m++) {
					//	rhs_[m] = rhs[k][j][i][m];
					//	rhs__[m] = rhs[k][j][i + 1][m];
					//}

					//for (m = 0; m < BLOCK_SIZE; m++) {
					//	rhs_[m] = rhs_[m] - lhs_buf[k][j][i][CC][0][m] * rhs__[0];
					//	rhs_[m] = rhs_[m] - lhs_buf[k][j][i][CC][1][m] * rhs__[1];
					//	rhs_[m] = rhs_[m] - lhs_buf[k][j][i][CC][2][m] * rhs__[2];
					//	rhs_[m] = rhs_[m] - lhs_buf[k][j][i][CC][3][m] * rhs__[3];
					//	rhs_[m] = rhs_[m] - lhs_buf[k][j][i][CC][4][m] * rhs__[4];
					//}

					//for (m = 0; m < BLOCK_SIZE; m++) {
					//	rhs[k][j][i][m] = rhs_[m];
					//}
				}
			}
		}
	//}
//#pragma dvm redistribute (lhs_buf[][][][][][])


	if (timeron) timer_stop(t_xsolve);

}



