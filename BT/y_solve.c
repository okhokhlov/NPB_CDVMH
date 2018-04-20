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
	/*
	for (k = 1; k <= grid_points[2]-2; k++) {
	  for (i = 1; i <= grid_points[0]-2; i++) {
		for (j = 0; j <= jsize; j++) {
		  tmp1 = rho_i[k][j][i];
		  tmp2 = tmp1 * tmp1;
		  tmp3 = tmp1 * tmp2;

		  fjac[j][0][0] = 0.0;
		  fjac[j][1][0] = 0.0;
		  fjac[j][2][0] = 1.0;
		  fjac[j][3][0] = 0.0;
		  fjac[j][4][0] = 0.0;

		  fjac[j][0][1] = - ( u[k][j][i][1]*u[k][j][i][2] ) * tmp2;
		  fjac[j][1][1] = u[k][j][i][2] * tmp1;
		  fjac[j][2][1] = u[k][j][i][1] * tmp1;
		  fjac[j][3][1] = 0.0;
		  fjac[j][4][1] = 0.0;

		  fjac[j][0][2] = - ( u[k][j][i][2]*u[k][j][i][2]*tmp2)
			+ c2 * qs[k][j][i];
		  fjac[j][1][2] = - c2 *  u[k][j][i][1] * tmp1;
		  fjac[j][2][2] = ( 2.0 - c2 ) *  u[k][j][i][2] * tmp1;
		  fjac[j][3][2] = - c2 * u[k][j][i][3] * tmp1;
		  fjac[j][4][2] = c2;

		  fjac[j][0][3] = - ( u[k][j][i][2]*u[k][j][i][3] ) * tmp2;
		  fjac[j][1][3] = 0.0;
		  fjac[j][2][3] = u[k][j][i][3] * tmp1;
		  fjac[j][3][3] = u[k][j][i][2] * tmp1;
		  fjac[j][4][3] = 0.0;

		  fjac[j][0][4] = ( c2 * 2.0 * square[k][j][i] - c1 * u[k][j][i][4] )
			* u[k][j][i][2] * tmp2;
		  fjac[j][1][4] = - c2 * u[k][j][i][1]*u[k][j][i][2] * tmp2;
		  fjac[j][2][4] = c1 * u[k][j][i][4] * tmp1
			- c2 * ( qs[k][j][i] + u[k][j][i][2]*u[k][j][i][2] * tmp2 );
		  fjac[j][3][4] = - c2 * ( u[k][j][i][2]*u[k][j][i][3] ) * tmp2;
		  fjac[j][4][4] = c1 * u[k][j][i][2] * tmp1;

		  njac[j][0][0] = 0.0;
		  njac[j][1][0] = 0.0;
		  njac[j][2][0] = 0.0;
		  njac[j][3][0] = 0.0;
		  njac[j][4][0] = 0.0;

		  njac[j][0][1] = - c3c4 * tmp2 * u[k][j][i][1];
		  njac[j][1][1] =   c3c4 * tmp1;
		  njac[j][2][1] =   0.0;
		  njac[j][3][1] =   0.0;
		  njac[j][4][1] =   0.0;

		  njac[j][0][2] = - con43 * c3c4 * tmp2 * u[k][j][i][2];
		  njac[j][1][2] =   0.0;
		  njac[j][2][2] =   con43 * c3c4 * tmp1;
		  njac[j][3][2] =   0.0;
		  njac[j][4][2] =   0.0;

		  njac[j][0][3] = - c3c4 * tmp2 * u[k][j][i][3];
		  njac[j][1][3] =   0.0;
		  njac[j][2][3] =   0.0;
		  njac[j][3][3] =   c3c4 * tmp1;
		  njac[j][4][3] =   0.0;

		  njac[j][0][4] = - (  c3c4
			  - c1345 ) * tmp3 * (u[k][j][i][1]*u[k][j][i][1])
			- ( con43 * c3c4
				- c1345 ) * tmp3 * (u[k][j][i][2]*u[k][j][i][2])
			- ( c3c4 - c1345 ) * tmp3 * (u[k][j][i][3]*u[k][j][i][3])
			- c1345 * tmp2 * u[k][j][i][4];

		  njac[j][1][4] = (  c3c4 - c1345 ) * tmp2 * u[k][j][i][1];
		  njac[j][2][4] = ( con43 * c3c4 - c1345 ) * tmp2 * u[k][j][i][2];
		  njac[j][3][4] = ( c3c4 - c1345 ) * tmp2 * u[k][j][i][3];
		  njac[j][4][4] = ( c1345 ) * tmp1;

		}

		//---------------------------------------------------------------------
		// now joacobians set, so form left hand side in y direction
		//---------------------------------------------------------------------
		lhsinit(lhs, jsize);
		for (j = 1; j <= jsize-1; j++) {
		  tmp1 = dt * ty1;
		  tmp2 = dt * ty2;

		  lhs[j][AA][0][0] = - tmp2 * fjac[j-1][0][0]
			- tmp1 * njac[j-1][0][0]
			- tmp1 * dy1;
		  lhs[j][AA][1][0] = - tmp2 * fjac[j-1][1][0]
			- tmp1 * njac[j-1][1][0];
		  lhs[j][AA][2][0] = - tmp2 * fjac[j-1][2][0]
			- tmp1 * njac[j-1][2][0];
		  lhs[j][AA][3][0] = - tmp2 * fjac[j-1][3][0]
			- tmp1 * njac[j-1][3][0];
		  lhs[j][AA][4][0] = - tmp2 * fjac[j-1][4][0]
			- tmp1 * njac[j-1][4][0];

		  lhs[j][AA][0][1] = - tmp2 * fjac[j-1][0][1]
			- tmp1 * njac[j-1][0][1];
		  lhs[j][AA][1][1] = - tmp2 * fjac[j-1][1][1]
			- tmp1 * njac[j-1][1][1]
			- tmp1 * dy2;
		  lhs[j][AA][2][1] = - tmp2 * fjac[j-1][2][1]
			- tmp1 * njac[j-1][2][1];
		  lhs[j][AA][3][1] = - tmp2 * fjac[j-1][3][1]
			- tmp1 * njac[j-1][3][1];
		  lhs[j][AA][4][1] = - tmp2 * fjac[j-1][4][1]
			- tmp1 * njac[j-1][4][1];

		  lhs[j][AA][0][2] = - tmp2 * fjac[j-1][0][2]
			- tmp1 * njac[j-1][0][2];
		  lhs[j][AA][1][2] = - tmp2 * fjac[j-1][1][2]
			- tmp1 * njac[j-1][1][2];
		  lhs[j][AA][2][2] = - tmp2 * fjac[j-1][2][2]
			- tmp1 * njac[j-1][2][2]
			- tmp1 * dy3;
		  lhs[j][AA][3][2] = - tmp2 * fjac[j-1][3][2]
			- tmp1 * njac[j-1][3][2];
		  lhs[j][AA][4][2] = - tmp2 * fjac[j-1][4][2]
			- tmp1 * njac[j-1][4][2];

		  lhs[j][AA][0][3] = - tmp2 * fjac[j-1][0][3]
			- tmp1 * njac[j-1][0][3];
		  lhs[j][AA][1][3] = - tmp2 * fjac[j-1][1][3]
			- tmp1 * njac[j-1][1][3];
		  lhs[j][AA][2][3] = - tmp2 * fjac[j-1][2][3]
			- tmp1 * njac[j-1][2][3];
		  lhs[j][AA][3][3] = - tmp2 * fjac[j-1][3][3]
			- tmp1 * njac[j-1][3][3]
			- tmp1 * dy4;
		  lhs[j][AA][4][3] = - tmp2 * fjac[j-1][4][3]
			- tmp1 * njac[j-1][4][3];

		  lhs[j][AA][0][4] = - tmp2 * fjac[j-1][0][4]
			- tmp1 * njac[j-1][0][4];
		  lhs[j][AA][1][4] = - tmp2 * fjac[j-1][1][4]
			- tmp1 * njac[j-1][1][4];
		  lhs[j][AA][2][4] = - tmp2 * fjac[j-1][2][4]
			- tmp1 * njac[j-1][2][4];
		  lhs[j][AA][3][4] = - tmp2 * fjac[j-1][3][4]
			- tmp1 * njac[j-1][3][4];
		  lhs[j][AA][4][4] = - tmp2 * fjac[j-1][4][4]
			- tmp1 * njac[j-1][4][4]
			- tmp1 * dy5;

		  lhs[j][BB][0][0] = 1.0
			+ tmp1 * 2.0 * njac[j][0][0]
			+ tmp1 * 2.0 * dy1;
		  lhs[j][BB][1][0] = tmp1 * 2.0 * njac[j][1][0];
		  lhs[j][BB][2][0] = tmp1 * 2.0 * njac[j][2][0];
		  lhs[j][BB][3][0] = tmp1 * 2.0 * njac[j][3][0];
		  lhs[j][BB][4][0] = tmp1 * 2.0 * njac[j][4][0];

		  lhs[j][BB][0][1] = tmp1 * 2.0 * njac[j][0][1];
		  lhs[j][BB][1][1] = 1.0
			+ tmp1 * 2.0 * njac[j][1][1]
			+ tmp1 * 2.0 * dy2;
		  lhs[j][BB][2][1] = tmp1 * 2.0 * njac[j][2][1];
		  lhs[j][BB][3][1] = tmp1 * 2.0 * njac[j][3][1];
		  lhs[j][BB][4][1] = tmp1 * 2.0 * njac[j][4][1];

		  lhs[j][BB][0][2] = tmp1 * 2.0 * njac[j][0][2];
		  lhs[j][BB][1][2] = tmp1 * 2.0 * njac[j][1][2];
		  lhs[j][BB][2][2] = 1.0
			+ tmp1 * 2.0 * njac[j][2][2]
			+ tmp1 * 2.0 * dy3;
		  lhs[j][BB][3][2] = tmp1 * 2.0 * njac[j][3][2];
		  lhs[j][BB][4][2] = tmp1 * 2.0 * njac[j][4][2];

		  lhs[j][BB][0][3] = tmp1 * 2.0 * njac[j][0][3];
		  lhs[j][BB][1][3] = tmp1 * 2.0 * njac[j][1][3];
		  lhs[j][BB][2][3] = tmp1 * 2.0 * njac[j][2][3];
		  lhs[j][BB][3][3] = 1.0
			+ tmp1 * 2.0 * njac[j][3][3]
			+ tmp1 * 2.0 * dy4;
		  lhs[j][BB][4][3] = tmp1 * 2.0 * njac[j][4][3];

		  lhs[j][BB][0][4] = tmp1 * 2.0 * njac[j][0][4];
		  lhs[j][BB][1][4] = tmp1 * 2.0 * njac[j][1][4];
		  lhs[j][BB][2][4] = tmp1 * 2.0 * njac[j][2][4];
		  lhs[j][BB][3][4] = tmp1 * 2.0 * njac[j][3][4];
		  lhs[j][BB][4][4] = 1.0
			+ tmp1 * 2.0 * njac[j][4][4]
			+ tmp1 * 2.0 * dy5;

		  lhs[j][CC][0][0] =  tmp2 * fjac[j+1][0][0]
			- tmp1 * njac[j+1][0][0]
			- tmp1 * dy1;
		  lhs[j][CC][1][0] =  tmp2 * fjac[j+1][1][0]
			- tmp1 * njac[j+1][1][0];
		  lhs[j][CC][2][0] =  tmp2 * fjac[j+1][2][0]
			- tmp1 * njac[j+1][2][0];
		  lhs[j][CC][3][0] =  tmp2 * fjac[j+1][3][0]
			- tmp1 * njac[j+1][3][0];
		  lhs[j][CC][4][0] =  tmp2 * fjac[j+1][4][0]
			- tmp1 * njac[j+1][4][0];

		  lhs[j][CC][0][1] =  tmp2 * fjac[j+1][0][1]
			- tmp1 * njac[j+1][0][1];
		  lhs[j][CC][1][1] =  tmp2 * fjac[j+1][1][1]
			- tmp1 * njac[j+1][1][1]
			- tmp1 * dy2;
		  lhs[j][CC][2][1] =  tmp2 * fjac[j+1][2][1]
			- tmp1 * njac[j+1][2][1];
		  lhs[j][CC][3][1] =  tmp2 * fjac[j+1][3][1]
			- tmp1 * njac[j+1][3][1];
		  lhs[j][CC][4][1] =  tmp2 * fjac[j+1][4][1]
			- tmp1 * njac[j+1][4][1];

		  lhs[j][CC][0][2] =  tmp2 * fjac[j+1][0][2]
			- tmp1 * njac[j+1][0][2];
		  lhs[j][CC][1][2] =  tmp2 * fjac[j+1][1][2]
			- tmp1 * njac[j+1][1][2];
		  lhs[j][CC][2][2] =  tmp2 * fjac[j+1][2][2]
			- tmp1 * njac[j+1][2][2]
			- tmp1 * dy3;
		  lhs[j][CC][3][2] =  tmp2 * fjac[j+1][3][2]
			- tmp1 * njac[j+1][3][2];
		  lhs[j][CC][4][2] =  tmp2 * fjac[j+1][4][2]
			- tmp1 * njac[j+1][4][2];

		  lhs[j][CC][0][3] =  tmp2 * fjac[j+1][0][3]
			- tmp1 * njac[j+1][0][3];
		  lhs[j][CC][1][3] =  tmp2 * fjac[j+1][1][3]
			- tmp1 * njac[j+1][1][3];
		  lhs[j][CC][2][3] =  tmp2 * fjac[j+1][2][3]
			- tmp1 * njac[j+1][2][3];
		  lhs[j][CC][3][3] =  tmp2 * fjac[j+1][3][3]
			- tmp1 * njac[j+1][3][3]
			- tmp1 * dy4;
		  lhs[j][CC][4][3] =  tmp2 * fjac[j+1][4][3]
			- tmp1 * njac[j+1][4][3];

		  lhs[j][CC][0][4] =  tmp2 * fjac[j+1][0][4]
			- tmp1 * njac[j+1][0][4];
		  lhs[j][CC][1][4] =  tmp2 * fjac[j+1][1][4]
			- tmp1 * njac[j+1][1][4];
		  lhs[j][CC][2][4] =  tmp2 * fjac[j+1][2][4]
			- tmp1 * njac[j+1][2][4];
		  lhs[j][CC][3][4] =  tmp2 * fjac[j+1][3][4]
			- tmp1 * njac[j+1][3][4];
		  lhs[j][CC][4][4] =  tmp2 * fjac[j+1][4][4]
			- tmp1 * njac[j+1][4][4]
			- tmp1 * dy5;
		}

		//---------------------------------------------------------------------
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		// performs guaussian elimination on this cell.
		//
		// assumes that unpacking routines for non-first cells
		// preload C' and rhs' from previous cell.
		//
		// assumed send happens outside this routine, but that
		// c'(JMAX) and rhs'(JMAX) will be sent to next cell
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		// multiply c[k][0][i] by b_inverse and copy back to c
		// multiply rhs(0) by b_inverse(0) and copy to rhs
		//---------------------------------------------------------------------

		//binvcrhs( lhs[0][BB], lhs[0][CC], rhs[k][0][i] );
		binvcrhs( 0, BB, 0, CC, k, 0, i );

		//---------------------------------------------------------------------
		// begin inner most do loop
		// do all the elements of the cell unless last
		//---------------------------------------------------------------------
		for (j = 1; j <= jsize-1; j++) {
		  //-------------------------------------------------------------------
		  // subtract A*lhs_vector(j-1) from lhs_vector(j)
		  //
		  // rhs(j) = rhs(j) - A*rhs(j-1)
		  //-------------------------------------------------------------------

		  //matvec_sub(lhs[j][AA], rhs[k][j-1][i], rhs[k][j][i]);
		  matvec_sub(j, AA, k, j-1, i, k, j, i);

		  //-------------------------------------------------------------------
		  // B(j) = B(j) - C(j-1)*A(j)
		  //-------------------------------------------------------------------

		  //matmul_sub(lhs[j][AA], lhs[j-1][CC], lhs[j][BB]);
		  matmul_sub(j, AA, j-1, CC, j, BB);

		  //-------------------------------------------------------------------
		  // multiply c[k][j][i] by b_inverse and copy back to c
		  // multiply rhs[k][0][i] by b_inverse[k][0][i] and copy to rhs
		  //-------------------------------------------------------------------

		  //binvcrhs( lhs[j][BB], lhs[j][CC], rhs[k][j][i] );
		  binvcrhs( j, BB, j, CC, k, j, i );
		}

		//---------------------------------------------------------------------
		// rhs(jsize) = rhs(jsize) - A*rhs(jsize-1)
		//---------------------------------------------------------------------

		//matvec_sub(lhs[jsize][AA], rhs[k][jsize-1][i], rhs[k][jsize][i]);
		matvec_sub(jsize, AA, k, jsize-1, i, k, jsize, i);

		//---------------------------------------------------------------------
		// B(jsize) = B(jsize) - C(jsize-1)*A(jsize)
		// matmul_sub(AA,i,jsize,k,c,
		// $              CC,i,jsize-1,k,c,BB,i,jsize,k)
		//---------------------------------------------------------------------

		//matmul_sub(lhs[jsize][AA], lhs[jsize-1][CC], lhs[jsize][BB]);
		matmul_sub(jsize, AA, jsize-1, CC, jsize, BB);

		//---------------------------------------------------------------------
		// multiply rhs(jsize) by b_inverse(jsize) and copy to rhs
		//---------------------------------------------------------------------

		//binvrhs( lhs[jsize][BB], rhs[k][jsize][i] );
		binvrhs( jsize, BB, k, jsize, i );

		//---------------------------------------------------------------------
		// back solve: if last cell, then generate U(jsize)=rhs(jsize)
		// else assume U(jsize) is loaded in un pack backsub_info
		// so just use it
		// after u(jstart) will be sent to next cell
		//---------------------------------------------------------------------
		for (j = jsize-1; j >= 0; j--) {
		  for (m = 0; m < BLOCK_SIZE; m++) {
			for (n = 0; n < BLOCK_SIZE; n++) {
			  rhs[k][j][i][m] = rhs[k][j][i][m]
				- lhs[j][CC][n][m]*rhs[k][j+1][i][n];
			}
		  }
		}
	  }
	}
	if (timeron) timer_stop(t_ysolve);
	 */

	for (k = 1; k <= grid_points[2] - 2; k++) {
		for (i = 1; i <= grid_points[0] - 2; i++) {




			//lhsinit(lhs, jsize);

			/*for (n = 0; n < 5; n++) {
				for (m = 0; m < 5; m++) {
					lhs[jsize][0][n][m] = 0.0;
					lhs[jsize][1][n][m] = 0.0;
					lhs[jsize][2][n][m] = 0.0;
				}
			}

			for (m = 0; m < 5; m++) {
				lhs[jsize][1][m][m] = 1.0;
			}
			*/
			
			for (j = 1; j <= jsize - 1; j++) {
				for (l = 0; l < 4; l++)
				{
					if (j - 1 + l <= jsize)
						fn_init_y(l, i, j - 1 + l, k);
				}
				if (j == 1)
				{
					rhs_y[0][0] = rhs[k][0][i][0];
					rhs_y[0][1] = rhs[k][0][i][1];
					rhs_y[0][2] = rhs[k][0][i][2];
					rhs_y[0][3] = rhs[k][0][i][3];
					rhs_y[0][4] = rhs[k][0][i][4];

					for (n = 0; n < 5; n++) {
						for (m = 0; m < 5; m++) {
							lhs_y[0][0][n][m] = 0.0;
							lhs_y[0][1][n][m] = 0.0;
							lhs_y[0][2][n][m] = 0.0;
						}
					}
					for (m = 0; m < 5; m++) {
						lhs_y[0][1][m][m] = 1.0;
					}

					binvcrhs_y(0, BB, 0, CC, 0);

					for (n = 0; n < 5; n++) {
						for (m = 0; m < 5; m++) {
							lhs[0][0][n][m] = lhs_y[0][0][n][m];
							lhs[0][1][n][m] = lhs_y[0][1][n][m];
							lhs[0][2][n][m] = lhs_y[0][2][n][m];
						}
					}

					rhs[k][0][i][0] = rhs_y[0][0];
					rhs[k][0][i][1] = rhs_y[0][1];
					rhs[k][0][i][2] = rhs_y[0][2];
					rhs[k][0][i][3] = rhs_y[0][3];
					rhs[k][0][i][4] = rhs_y[0][4];
				}
				else
				{
					for (n = 0; n < 5; n++) {
						for (m = 0; m < 5; m++) {
							lhs_y[0][0][n][m] = lhs_y[1][0][n][m];
							lhs_y[0][1][n][m] = lhs_y[1][1][n][m];
							lhs_y[0][2][n][m] = lhs_y[1][2][n][m];
						}
					}
					//buf_compute_(1,0);
					rhs_y[0][0] = rhs_y[1][0];
					rhs_y[0][1] = rhs_y[1][1];
					rhs_y[0][2] = rhs_y[1][2];
					rhs_y[0][3] = rhs_y[1][3];
					rhs_y[0][4] = rhs_y[1][4];
				}

				

				rhs_y[1][0] = rhs[k][j][i][0];
				rhs_y[1][1] = rhs[k][j][i][1];
				rhs_y[1][2] = rhs[k][j][i][2];
				rhs_y[1][3] = rhs[k][j][i][3];
				rhs_y[1][4] = rhs[k][j][i][4];


				buf_compute_y(1, 1);
				matvec_sub_y(1, AA, 0, 1);
				matmul_sub_y(1, AA, 0, CC, 1, BB);
				binvcrhs_y(1, BB, 1, CC, 1);


				for (n = 0; n < 5; n++) {
					for (m = 0; m < 5; m++) {
						lhs[j][0][n][m] = lhs_y[1][0][n][m];
						lhs[j][1][n][m] = lhs_y[1][1][n][m];
						lhs[j][2][n][m] = lhs_y[1][2][n][m];
					}
				}

				rhs[k][j][i][0] = rhs_y[1][0];
				rhs[k][j][i][1] = rhs_y[1][1];
				rhs[k][j][i][2] = rhs_y[1][2];
				rhs[k][j][i][3] = rhs_y[1][3];
				rhs[k][j][i][4] = rhs_y[1][4];
				
				if(j == jsize-1)
				{
					for (n = 0; n < 5; n++) {
						for (m = 0; m < 5; m++) {
						lhs[jsize][0][n][m] = 0.0;
						lhs[jsize][1][n][m] = 0.0;
						lhs[jsize][2][n][m] = 0.0;
						}
					}

					for (m = 0; m < 5; m++) {
					lhs[jsize][1][m][m] = 1.0;
					}
					
					rhs_y[1][0] = rhs[k][jsize][i][0];
					rhs_y[1][1] = rhs[k][jsize][i][1];
					rhs_y[1][2] = rhs[k][jsize][i][2];
					rhs_y[1][3] = rhs[k][jsize][i][3];
					rhs_y[1][4] = rhs[k][jsize][i][4];
					
					for (n = 0; n < 5; n++) {
						for (m = 0; m < 5; m++) {
							lhs_y[1][0][n][m] = lhs[jsize][0][n][m];
							lhs_y[1][1][n][m] = lhs[jsize][1][n][m];
							lhs_y[1][2][n][m] = lhs[jsize][2][n][m];
						}
					}
					
					for (n = 0; n < 5; n++) {
						for (m = 0; m < 5; m++) {
							lhs_y[0][0][n][m] = lhs[jsize-1][0][n][m];
							lhs_y[0][1][n][m] = lhs[jsize-1][1][n][m];
							lhs_y[0][2][n][m] = lhs[jsize-1][2][n][m];
						}
					}

					matvec_sub_y(1, AA, 0, 1);
					matmul_sub_y(1, AA,0, CC, 1, BB);
					binvrhs_y(1, BB, 1);
					
					for (n = 0; n < 5; n++) {
						for (m = 0; m < 5; m++) {
							lhs[jsize][0][n][m] = lhs_y[1][0][n][m];
							lhs[jsize][1][n][m] = lhs_y[1][1][n][m];
							lhs[jsize][2][n][m] = lhs_y[1][2][n][m];
						}
					}
					rhs[k][jsize][i][0] = rhs_y[1][0];
					rhs[k][jsize][i][1] = rhs_y[1][1];
					rhs[k][jsize][i][2] = rhs_y[1][2];
					rhs[k][jsize][i][3] = rhs_y[1][3];
					rhs[k][jsize][i][4] = rhs_y[1][4];
					
					
				}
				
				for (n = 0; n < 5; n++) {
					for (m = 0; m < 5; m++) {
						lhs_buf[k][j][i][n][m] = lhs[j][CC][n][m];
					}
				}
				if(i == jsize-1)
				{
					for (n = 0; n < 5; n++) {
						for (m = 0; m < 5; m++) {
							lhs_buf[k][0][i][n][m] = lhs[0][CC][n][m];
							lhs_buf[k][jsize][i][n][m] = lhs[jsize][CC][n][m];
						}
					}
				
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
							- lhs_buf[k][j][i][n][m] * rhs[k][j + 1][i][n];
					}
				}
			}
		}
	}

}


void fn_init_y(int a, int i, int j, int k)
{
	tmp1 = rho_i[k][j][i];
	tmp2 = tmp1 * tmp1;
	tmp3 = tmp1 * tmp2;

	fjac_y[a][0][0] = 0.0;
	fjac_y[a][1][0] = 0.0;
	fjac_y[a][2][0] = 1.0;
	fjac_y[a][3][0] = 0.0;
	fjac_y[a][4][0] = 0.0;

	fjac_y[a][0][1] = -(u[k][j][i][1] * u[k][j][i][2]) * tmp2;
	fjac_y[a][1][1] = u[k][j][i][2] * tmp1;
	fjac_y[a][2][1] = u[k][j][i][1] * tmp1;
	fjac_y[a][3][1] = 0.0;
	fjac_y[a][4][1] = 0.0;

	fjac_y[a][0][2] = -(u[k][j][i][2] * u[k][j][i][2] * tmp2)
		+ c2 * qs[k][j][i];
	fjac_y[a][1][2] = -c2 *  u[k][j][i][1] * tmp1;
	fjac_y[a][2][2] = (2.0 - c2) *  u[k][j][i][2] * tmp1;
	fjac_y[a][3][2] = -c2 * u[k][j][i][3] * tmp1;
	fjac_y[a][4][2] = c2;

	fjac_y[a][0][3] = -(u[k][j][i][2] * u[k][j][i][3]) * tmp2;
	fjac_y[a][1][3] = 0.0;
	fjac_y[a][2][3] = u[k][j][i][3] * tmp1;
	fjac_y[a][3][3] = u[k][j][i][2] * tmp1;
	fjac_y[a][4][3] = 0.0;

	fjac_y[a][0][4] = (c2 * 2.0 * square[k][j][i] - c1 * u[k][j][i][4])
		* u[k][j][i][2] * tmp2;
	fjac_y[a][1][4] = -c2 * u[k][j][i][1] * u[k][j][i][2] * tmp2;
	fjac_y[a][2][4] = c1 * u[k][j][i][4] * tmp1
		- c2 * (qs[k][j][i] + u[k][j][i][2] * u[k][j][i][2] * tmp2);
	fjac_y[a][3][4] = -c2 * (u[k][j][i][2] * u[k][j][i][3]) * tmp2;
	fjac_y[a][4][4] = c1 * u[k][j][i][2] * tmp1;

	njac_y[a][0][0] = 0.0;
	njac_y[a][1][0] = 0.0;
	njac_y[a][2][0] = 0.0;
	njac_y[a][3][0] = 0.0;
	njac_y[a][4][0] = 0.0;

	njac_y[a][0][1] = -c3c4 * tmp2 * u[k][j][i][1];
	njac_y[a][1][1] = c3c4 * tmp1;
	njac_y[a][2][1] = 0.0;
	njac_y[a][3][1] = 0.0;
	njac_y[a][4][1] = 0.0;

	njac_y[a][0][2] = -con43 * c3c4 * tmp2 * u[k][j][i][2];
	njac_y[a][1][2] = 0.0;
	njac_y[a][2][2] = con43 * c3c4 * tmp1;
	njac_y[a][3][2] = 0.0;
	njac_y[a][4][2] = 0.0;

	njac_y[a][0][3] = -c3c4 * tmp2 * u[k][j][i][3];
	njac_y[a][1][3] = 0.0;
	njac_y[a][2][3] = 0.0;
	njac_y[a][3][3] = c3c4 * tmp1;
	njac_y[a][4][3] = 0.0;

	njac_y[a][0][4] = -(c3c4
		- c1345) * tmp3 * (u[k][j][i][1] * u[k][j][i][1])
		- (con43 * c3c4
			- c1345) * tmp3 * (u[k][j][i][2] * u[k][j][i][2])
		- (c3c4 - c1345) * tmp3 * (u[k][j][i][3] * u[k][j][i][3])
		- c1345 * tmp2 * u[k][j][i][4];

	njac_y[a][1][4] = (c3c4 - c1345) * tmp2 * u[k][j][i][1];
	njac_y[a][2][4] = (con43 * c3c4 - c1345) * tmp2 * u[k][j][i][2];
	njac_y[a][3][4] = (c3c4 - c1345) * tmp2 * u[k][j][i][3];
	njac_y[a][4][4] = (c1345)* tmp1;

}

void matvec_sub_y(int a1, int a2, int b2, int c2)
{
	//---------------------------------------------------------------------
	// rhs[kc][jc][ic][i] = rhs[kc][jc][ic][i] 
	// $                  - lhs[ia][lhs[a1][a2]][0][i]*
	//---------------------------------------------------------------------
	rhs_y[c2][0] = rhs_y[c2][0] - lhs_y[a1][a2][0][0] * rhs_y[b2][0]
		- lhs_y[a1][a2][1][0] * rhs_y[b2][1]
		- lhs_y[a1][a2][2][0] * rhs_y[b2][2]
		- lhs_y[a1][a2][3][0] * rhs_y[b2][3]
		- lhs_y[a1][a2][4][0] * rhs_y[b2][4];
	rhs_y[c2][1] = rhs_y[c2][1] - lhs_y[a1][a2][0][1] * rhs_y[b2][0]
		- lhs_y[a1][a2][1][1] * rhs_y[b2][1]
		- lhs_y[a1][a2][2][1] * rhs_y[b2][2]
		- lhs_y[a1][a2][3][1] * rhs_y[b2][3]
		- lhs_y[a1][a2][4][1] * rhs_y[b2][4];
	rhs_y[c2][2] = rhs_y[c2][2] - lhs_y[a1][a2][0][2] * rhs_y[b2][0]
		- lhs_y[a1][a2][1][2] * rhs_y[b2][1]
		- lhs_y[a1][a2][2][2] * rhs_y[b2][2]
		- lhs_y[a1][a2][3][2] * rhs_y[b2][3]
		- lhs_y[a1][a2][4][2] * rhs_y[b2][4];
	rhs_y[c2][3] = rhs_y[c2][3] - lhs_y[a1][a2][0][3] * rhs_y[b2][0]
		- lhs_y[a1][a2][1][3] * rhs_y[b2][1]
		- lhs_y[a1][a2][2][3] * rhs_y[b2][2]
		- lhs_y[a1][a2][3][3] * rhs_y[b2][3]
		- lhs_y[a1][a2][4][3] * rhs_y[b2][4];
	rhs_y[c2][4] = rhs_y[c2][4] - lhs_y[a1][a2][0][4] * rhs_y[b2][0]
		- lhs_y[a1][a2][1][4] * rhs_y[b2][1]
		- lhs_y[a1][a2][2][4] * rhs_y[b2][2]
		- lhs_y[a1][a2][3][4] * rhs_y[b2][3]
		- lhs_y[a1][a2][4][4] * rhs_y[b2][4];
}

void buf_compute_y(int j, int a)
{
	tmp1 = dt * ty1;
	tmp2 = dt * ty2;

	lhs_y[a][AA][0][0] = -tmp2 * fjac_y[j - 1][0][0]
		- tmp1 * njac_y[j - 1][0][0]
		- tmp1 * dy1;
	lhs_y[a][AA][1][0] = -tmp2 * fjac_y[j - 1][1][0]
		- tmp1 * njac_y[j - 1][1][0];
	lhs_y[a][AA][2][0] = -tmp2 * fjac_y[j - 1][2][0]
		- tmp1 * njac_y[j - 1][2][0];
	lhs_y[a][AA][3][0] = -tmp2 * fjac_y[j - 1][3][0]
		- tmp1 * njac_y[j - 1][3][0];
	lhs_y[a][AA][4][0] = -tmp2 * fjac_y[j - 1][4][0]
		- tmp1 * njac_y[j - 1][4][0];

	lhs_y[a][AA][0][1] = -tmp2 * fjac_y[j - 1][0][1]
		- tmp1 * njac_y[j - 1][0][1];
	lhs_y[a][AA][1][1] = -tmp2 * fjac_y[j - 1][1][1]
		- tmp1 * njac_y[j - 1][1][1]
		- tmp1 * dy2;
	lhs_y[a][AA][2][1] = -tmp2 * fjac_y[j - 1][2][1]
		- tmp1 * njac_y[j - 1][2][1];
	lhs_y[a][AA][3][1] = -tmp2 * fjac_y[j - 1][3][1]
		- tmp1 * njac_y[j - 1][3][1];
	lhs_y[a][AA][4][1] = -tmp2 * fjac_y[j - 1][4][1]
		- tmp1 * njac_y[j - 1][4][1];

	lhs_y[a][AA][0][2] = -tmp2 * fjac_y[j - 1][0][2]
		- tmp1 * njac_y[j - 1][0][2];
	lhs_y[a][AA][1][2] = -tmp2 * fjac_y[j - 1][1][2]
		- tmp1 * njac_y[j - 1][1][2];
	lhs_y[a][AA][2][2] = -tmp2 * fjac_y[j - 1][2][2]
		- tmp1 * njac_y[j - 1][2][2]
		- tmp1 * dy3;
	lhs_y[a][AA][3][2] = -tmp2 * fjac_y[j - 1][3][2]
		- tmp1 * njac_y[j - 1][3][2];
	lhs_y[a][AA][4][2] = -tmp2 * fjac_y[j - 1][4][2]
		- tmp1 * njac_y[j - 1][4][2];

	lhs_y[a][AA][0][3] = -tmp2 * fjac_y[j - 1][0][3]
		- tmp1 * njac_y[j - 1][0][3];
	lhs_y[a][AA][1][3] = -tmp2 * fjac_y[j - 1][1][3]
		- tmp1 * njac_y[j - 1][1][3];
	lhs_y[a][AA][2][3] = -tmp2 * fjac_y[j - 1][2][3]
		- tmp1 * njac_y[j - 1][2][3];
	lhs_y[a][AA][3][3] = -tmp2 * fjac_y[j - 1][3][3]
		- tmp1 * njac_y[j - 1][3][3]
		- tmp1 * dy4;
	lhs_y[a][AA][4][3] = -tmp2 * fjac_y[j - 1][4][3]
		- tmp1 * njac_y[j - 1][4][3];

	lhs_y[a][AA][0][4] = -tmp2 * fjac_y[j - 1][0][4]
		- tmp1 * njac_y[j - 1][0][4];
	lhs_y[a][AA][1][4] = -tmp2 * fjac_y[j - 1][1][4]
		- tmp1 * njac_y[j - 1][1][4];
	lhs_y[a][AA][2][4] = -tmp2 * fjac_y[j - 1][2][4]
		- tmp1 * njac_y[j - 1][2][4];
	lhs_y[a][AA][3][4] = -tmp2 * fjac_y[j - 1][3][4]
		- tmp1 * njac_y[j - 1][3][4];
	lhs_y[a][AA][4][4] = -tmp2 * fjac_y[j - 1][4][4]
		- tmp1 * njac_y[j - 1][4][4]
		- tmp1 * dy5;

	lhs_y[a][BB][0][0] = 1.0
		+ tmp1 * 2.0 * njac_y[j][0][0]
		+ tmp1 * 2.0 * dy1;
	lhs_y[a][BB][1][0] = tmp1 * 2.0 * njac_y[j][1][0];
	lhs_y[a][BB][2][0] = tmp1 * 2.0 * njac_y[j][2][0];
	lhs_y[a][BB][3][0] = tmp1 * 2.0 * njac_y[j][3][0];
	lhs_y[a][BB][4][0] = tmp1 * 2.0 * njac_y[j][4][0];

	lhs_y[a][BB][0][1] = tmp1 * 2.0 * njac_y[j][0][1];
	lhs_y[a][BB][1][1] = 1.0
		+ tmp1 * 2.0 * njac_y[j][1][1]
		+ tmp1 * 2.0 * dy2;
	lhs_y[a][BB][2][1] = tmp1 * 2.0 * njac_y[j][2][1];
	lhs_y[a][BB][3][1] = tmp1 * 2.0 * njac_y[j][3][1];
	lhs_y[a][BB][4][1] = tmp1 * 2.0 * njac_y[j][4][1];

	lhs_y[a][BB][0][2] = tmp1 * 2.0 * njac_y[j][0][2];
	lhs_y[a][BB][1][2] = tmp1 * 2.0 * njac_y[j][1][2];
	lhs_y[a][BB][2][2] = 1.0
		+ tmp1 * 2.0 * njac_y[j][2][2]
		+ tmp1 * 2.0 * dy3;
	lhs_y[a][BB][3][2] = tmp1 * 2.0 * njac_y[j][3][2];
	lhs_y[a][BB][4][2] = tmp1 * 2.0 * njac_y[j][4][2];

	lhs_y[a][BB][0][3] = tmp1 * 2.0 * njac_y[j][0][3];
	lhs_y[a][BB][1][3] = tmp1 * 2.0 * njac_y[j][1][3];
	lhs_y[a][BB][2][3] = tmp1 * 2.0 * njac_y[j][2][3];
	lhs_y[a][BB][3][3] = 1.0
		+ tmp1 * 2.0 * njac_y[j][3][3]
		+ tmp1 * 2.0 * dy4;
	lhs_y[a][BB][4][3] = tmp1 * 2.0 * njac_y[j][4][3];

	lhs_y[a][BB][0][4] = tmp1 * 2.0 * njac_y[j][0][4];
	lhs_y[a][BB][1][4] = tmp1 * 2.0 * njac_y[j][1][4];
	lhs_y[a][BB][2][4] = tmp1 * 2.0 * njac_y[j][2][4];
	lhs_y[a][BB][3][4] = tmp1 * 2.0 * njac_y[j][3][4];
	lhs_y[a][BB][4][4] = 1.0
		+ tmp1 * 2.0 * njac_y[j][4][4]
		+ tmp1 * 2.0 * dy5;

	lhs_y[a][CC][0][0] = tmp2 * fjac_y[j + 1][0][0]
		- tmp1 * njac_y[j + 1][0][0]
		- tmp1 * dy1;
	lhs_y[a][CC][1][0] = tmp2 * fjac_y[j + 1][1][0]
		- tmp1 * njac_y[j + 1][1][0];
	lhs_y[a][CC][2][0] = tmp2 * fjac_y[j + 1][2][0]
		- tmp1 * njac_y[j + 1][2][0];
	lhs_y[a][CC][3][0] = tmp2 * fjac_y[j + 1][3][0]
		- tmp1 * njac_y[j + 1][3][0];
	lhs_y[a][CC][4][0] = tmp2 * fjac_y[j + 1][4][0]
		- tmp1 * njac_y[j + 1][4][0];

	lhs_y[a][CC][0][1] = tmp2 * fjac_y[j + 1][0][1]
		- tmp1 * njac_y[j + 1][0][1];
	lhs_y[a][CC][1][1] = tmp2 * fjac_y[j + 1][1][1]
		- tmp1 * njac_y[j + 1][1][1]
		- tmp1 * dy2;
	lhs_y[a][CC][2][1] = tmp2 * fjac_y[j + 1][2][1]
		- tmp1 * njac_y[j + 1][2][1];
	lhs_y[a][CC][3][1] = tmp2 * fjac_y[j + 1][3][1]
		- tmp1 * njac_y[j + 1][3][1];
	lhs_y[a][CC][4][1] = tmp2 * fjac_y[j + 1][4][1]
		- tmp1 * njac_y[j + 1][4][1];

	lhs_y[a][CC][0][2] = tmp2 * fjac_y[j + 1][0][2]
		- tmp1 * njac_y[j + 1][0][2];
	lhs_y[a][CC][1][2] = tmp2 * fjac_y[j + 1][1][2]
		- tmp1 * njac_y[j + 1][1][2];
	lhs_y[a][CC][2][2] = tmp2 * fjac_y[j + 1][2][2]
		- tmp1 * njac_y[j + 1][2][2]
		- tmp1 * dy3;
	lhs_y[a][CC][3][2] = tmp2 * fjac_y[j + 1][3][2]
		- tmp1 * njac_y[j + 1][3][2];
	lhs_y[a][CC][4][2] = tmp2 * fjac_y[j + 1][4][2]
		- tmp1 * njac_y[j + 1][4][2];

	lhs_y[a][CC][0][3] = tmp2 * fjac_y[j + 1][0][3]
		- tmp1 * njac_y[j + 1][0][3];
	lhs_y[a][CC][1][3] = tmp2 * fjac_y[j + 1][1][3]
		- tmp1 * njac_y[j + 1][1][3];
	lhs_y[a][CC][2][3] = tmp2 * fjac_y[j + 1][2][3]
		- tmp1 * njac_y[j + 1][2][3];
	lhs_y[a][CC][3][3] = tmp2 * fjac_y[j + 1][3][3]
		- tmp1 * njac_y[j + 1][3][3]
		- tmp1 * dy4;
	lhs_y[a][CC][4][3] = tmp2 * fjac_y[j + 1][4][3]
		- tmp1 * njac_y[j + 1][4][3];

	lhs_y[a][CC][0][4] = tmp2 * fjac_y[j + 1][0][4]
		- tmp1 * njac_y[j + 1][0][4];
	lhs_y[a][CC][1][4] = tmp2 * fjac_y[j + 1][1][4]
		- tmp1 * njac_y[j + 1][1][4];
	lhs_y[a][CC][2][4] = tmp2 * fjac_y[j + 1][2][4]
		- tmp1 * njac_y[j + 1][2][4];
	lhs_y[a][CC][3][4] = tmp2 * fjac_y[j + 1][3][4]
		- tmp1 * njac_y[j + 1][3][4];
	lhs_y[a][CC][4][4] = tmp2 * fjac_y[j + 1][4][4]
		- tmp1 * njac_y[j + 1][4][4]
		- tmp1 * dy5;
}

void matmul_sub_y(int a1, int a2, int b1, int b2, int c1, int c2)
{
	lhs_y[c1][c2][0][0] = lhs_y[c1][c2][0][0] - lhs_y[a1][a2][0][0] * lhs_y[b1][b2][0][0]
		- lhs_y[a1][a2][1][0] * lhs_y[b1][b2][0][1]
		- lhs_y[a1][a2][2][0] * lhs_y[b1][b2][0][2]
		- lhs_y[a1][a2][3][0] * lhs_y[b1][b2][0][3]
		- lhs_y[a1][a2][4][0] * lhs_y[b1][b2][0][4];
	lhs_y[c1][c2][0][1] = lhs_y[c1][c2][0][1] - lhs_y[a1][a2][0][1] * lhs_y[b1][b2][0][0]
		- lhs_y[a1][a2][1][1] * lhs_y[b1][b2][0][1]
		- lhs_y[a1][a2][2][1] * lhs_y[b1][b2][0][2]
		- lhs_y[a1][a2][3][1] * lhs_y[b1][b2][0][3]
		- lhs_y[a1][a2][4][1] * lhs_y[b1][b2][0][4];
	lhs_y[c1][c2][0][2] = lhs_y[c1][c2][0][2] - lhs_y[a1][a2][0][2] * lhs_y[b1][b2][0][0]
		- lhs_y[a1][a2][1][2] * lhs_y[b1][b2][0][1]
		- lhs_y[a1][a2][2][2] * lhs_y[b1][b2][0][2]
		- lhs_y[a1][a2][3][2] * lhs_y[b1][b2][0][3]
		- lhs_y[a1][a2][4][2] * lhs_y[b1][b2][0][4];
	lhs_y[c1][c2][0][3] = lhs_y[c1][c2][0][3] - lhs_y[a1][a2][0][3] * lhs_y[b1][b2][0][0]
		- lhs_y[a1][a2][1][3] * lhs_y[b1][b2][0][1]
		- lhs_y[a1][a2][2][3] * lhs_y[b1][b2][0][2]
		- lhs_y[a1][a2][3][3] * lhs_y[b1][b2][0][3]
		- lhs_y[a1][a2][4][3] * lhs_y[b1][b2][0][4];
	lhs_y[c1][c2][0][4] = lhs_y[c1][c2][0][4] - lhs_y[a1][a2][0][4] * lhs_y[b1][b2][0][0]
		- lhs_y[a1][a2][1][4] * lhs_y[b1][b2][0][1]
		- lhs_y[a1][a2][2][4] * lhs_y[b1][b2][0][2]
		- lhs_y[a1][a2][3][4] * lhs_y[b1][b2][0][3]
		- lhs_y[a1][a2][4][4] * lhs_y[b1][b2][0][4];
	lhs_y[c1][c2][1][0] = lhs_y[c1][c2][1][0] - lhs_y[a1][a2][0][0] * lhs_y[b1][b2][1][0]
		- lhs_y[a1][a2][1][0] * lhs_y[b1][b2][1][1]
		- lhs_y[a1][a2][2][0] * lhs_y[b1][b2][1][2]
		- lhs_y[a1][a2][3][0] * lhs_y[b1][b2][1][3]
		- lhs_y[a1][a2][4][0] * lhs_y[b1][b2][1][4];
	lhs_y[c1][c2][1][1] = lhs_y[c1][c2][1][1] - lhs_y[a1][a2][0][1] * lhs_y[b1][b2][1][0]
		- lhs_y[a1][a2][1][1] * lhs_y[b1][b2][1][1]
		- lhs_y[a1][a2][2][1] * lhs_y[b1][b2][1][2]
		- lhs_y[a1][a2][3][1] * lhs_y[b1][b2][1][3]
		- lhs_y[a1][a2][4][1] * lhs_y[b1][b2][1][4];
	lhs_y[c1][c2][1][2] = lhs_y[c1][c2][1][2] - lhs_y[a1][a2][0][2] * lhs_y[b1][b2][1][0]
		- lhs_y[a1][a2][1][2] * lhs_y[b1][b2][1][1]
		- lhs_y[a1][a2][2][2] * lhs_y[b1][b2][1][2]
		- lhs_y[a1][a2][3][2] * lhs_y[b1][b2][1][3]
		- lhs_y[a1][a2][4][2] * lhs_y[b1][b2][1][4];
	lhs_y[c1][c2][1][3] = lhs_y[c1][c2][1][3] - lhs_y[a1][a2][0][3] * lhs_y[b1][b2][1][0]
		- lhs_y[a1][a2][1][3] * lhs_y[b1][b2][1][1]
		- lhs_y[a1][a2][2][3] * lhs_y[b1][b2][1][2]
		- lhs_y[a1][a2][3][3] * lhs_y[b1][b2][1][3]
		- lhs_y[a1][a2][4][3] * lhs_y[b1][b2][1][4];
	lhs_y[c1][c2][1][4] = lhs_y[c1][c2][1][4] - lhs_y[a1][a2][0][4] * lhs_y[b1][b2][1][0]
		- lhs_y[a1][a2][1][4] * lhs_y[b1][b2][1][1]
		- lhs_y[a1][a2][2][4] * lhs_y[b1][b2][1][2]
		- lhs_y[a1][a2][3][4] * lhs_y[b1][b2][1][3]
		- lhs_y[a1][a2][4][4] * lhs_y[b1][b2][1][4];
	lhs_y[c1][c2][2][0] = lhs_y[c1][c2][2][0] - lhs_y[a1][a2][0][0] * lhs_y[b1][b2][2][0]
		- lhs_y[a1][a2][1][0] * lhs_y[b1][b2][2][1]
		- lhs_y[a1][a2][2][0] * lhs_y[b1][b2][2][2]
		- lhs_y[a1][a2][3][0] * lhs_y[b1][b2][2][3]
		- lhs_y[a1][a2][4][0] * lhs_y[b1][b2][2][4];
	lhs_y[c1][c2][2][1] = lhs_y[c1][c2][2][1] - lhs_y[a1][a2][0][1] * lhs_y[b1][b2][2][0]
		- lhs_y[a1][a2][1][1] * lhs_y[b1][b2][2][1]
		- lhs_y[a1][a2][2][1] * lhs_y[b1][b2][2][2]
		- lhs_y[a1][a2][3][1] * lhs_y[b1][b2][2][3]
		- lhs_y[a1][a2][4][1] * lhs_y[b1][b2][2][4];
	lhs_y[c1][c2][2][2] = lhs_y[c1][c2][2][2] - lhs_y[a1][a2][0][2] * lhs_y[b1][b2][2][0]
		- lhs_y[a1][a2][1][2] * lhs_y[b1][b2][2][1]
		- lhs_y[a1][a2][2][2] * lhs_y[b1][b2][2][2]
		- lhs_y[a1][a2][3][2] * lhs_y[b1][b2][2][3]
		- lhs_y[a1][a2][4][2] * lhs_y[b1][b2][2][4];
	lhs_y[c1][c2][2][3] = lhs_y[c1][c2][2][3] - lhs_y[a1][a2][0][3] * lhs_y[b1][b2][2][0]
		- lhs_y[a1][a2][1][3] * lhs_y[b1][b2][2][1]
		- lhs_y[a1][a2][2][3] * lhs_y[b1][b2][2][2]
		- lhs_y[a1][a2][3][3] * lhs_y[b1][b2][2][3]
		- lhs_y[a1][a2][4][3] * lhs_y[b1][b2][2][4];
	lhs_y[c1][c2][2][4] = lhs_y[c1][c2][2][4] - lhs_y[a1][a2][0][4] * lhs_y[b1][b2][2][0]
		- lhs_y[a1][a2][1][4] * lhs_y[b1][b2][2][1]
		- lhs_y[a1][a2][2][4] * lhs_y[b1][b2][2][2]
		- lhs_y[a1][a2][3][4] * lhs_y[b1][b2][2][3]
		- lhs_y[a1][a2][4][4] * lhs_y[b1][b2][2][4];
	lhs_y[c1][c2][3][0] = lhs_y[c1][c2][3][0] - lhs_y[a1][a2][0][0] * lhs_y[b1][b2][3][0]
		- lhs_y[a1][a2][1][0] * lhs_y[b1][b2][3][1]
		- lhs_y[a1][a2][2][0] * lhs_y[b1][b2][3][2]
		- lhs_y[a1][a2][3][0] * lhs_y[b1][b2][3][3]
		- lhs_y[a1][a2][4][0] * lhs_y[b1][b2][3][4];
	lhs_y[c1][c2][3][1] = lhs_y[c1][c2][3][1] - lhs_y[a1][a2][0][1] * lhs_y[b1][b2][3][0]
		- lhs_y[a1][a2][1][1] * lhs_y[b1][b2][3][1]
		- lhs_y[a1][a2][2][1] * lhs_y[b1][b2][3][2]
		- lhs_y[a1][a2][3][1] * lhs_y[b1][b2][3][3]
		- lhs_y[a1][a2][4][1] * lhs_y[b1][b2][3][4];
	lhs_y[c1][c2][3][2] = lhs_y[c1][c2][3][2] - lhs_y[a1][a2][0][2] * lhs_y[b1][b2][3][0]
		- lhs_y[a1][a2][1][2] * lhs_y[b1][b2][3][1]
		- lhs_y[a1][a2][2][2] * lhs_y[b1][b2][3][2]
		- lhs_y[a1][a2][3][2] * lhs_y[b1][b2][3][3]
		- lhs_y[a1][a2][4][2] * lhs_y[b1][b2][3][4];
	lhs_y[c1][c2][3][3] = lhs_y[c1][c2][3][3] - lhs_y[a1][a2][0][3] * lhs_y[b1][b2][3][0]
		- lhs_y[a1][a2][1][3] * lhs_y[b1][b2][3][1]
		- lhs_y[a1][a2][2][3] * lhs_y[b1][b2][3][2]
		- lhs_y[a1][a2][3][3] * lhs_y[b1][b2][3][3]
		- lhs_y[a1][a2][4][3] * lhs_y[b1][b2][3][4];
	lhs_y[c1][c2][3][4] = lhs_y[c1][c2][3][4] - lhs_y[a1][a2][0][4] * lhs_y[b1][b2][3][0]
		- lhs_y[a1][a2][1][4] * lhs_y[b1][b2][3][1]
		- lhs_y[a1][a2][2][4] * lhs_y[b1][b2][3][2]
		- lhs_y[a1][a2][3][4] * lhs_y[b1][b2][3][3]
		- lhs_y[a1][a2][4][4] * lhs_y[b1][b2][3][4];
	lhs_y[c1][c2][4][0] = lhs_y[c1][c2][4][0] - lhs_y[a1][a2][0][0] * lhs_y[b1][b2][4][0]
		- lhs_y[a1][a2][1][0] * lhs_y[b1][b2][4][1]
		- lhs_y[a1][a2][2][0] * lhs_y[b1][b2][4][2]
		- lhs_y[a1][a2][3][0] * lhs_y[b1][b2][4][3]
		- lhs_y[a1][a2][4][0] * lhs_y[b1][b2][4][4];
	lhs_y[c1][c2][4][1] = lhs_y[c1][c2][4][1] - lhs_y[a1][a2][0][1] * lhs_y[b1][b2][4][0]
		- lhs_y[a1][a2][1][1] * lhs_y[b1][b2][4][1]
		- lhs_y[a1][a2][2][1] * lhs_y[b1][b2][4][2]
		- lhs_y[a1][a2][3][1] * lhs_y[b1][b2][4][3]
		- lhs_y[a1][a2][4][1] * lhs_y[b1][b2][4][4];
	lhs_y[c1][c2][4][2] = lhs_y[c1][c2][4][2] - lhs_y[a1][a2][0][2] * lhs_y[b1][b2][4][0]
		- lhs_y[a1][a2][1][2] * lhs_y[b1][b2][4][1]
		- lhs_y[a1][a2][2][2] * lhs_y[b1][b2][4][2]
		- lhs_y[a1][a2][3][2] * lhs_y[b1][b2][4][3]
		- lhs_y[a1][a2][4][2] * lhs_y[b1][b2][4][4];
	lhs_y[c1][c2][4][3] = lhs_y[c1][c2][4][3] - lhs_y[a1][a2][0][3] * lhs_y[b1][b2][4][0]
		- lhs_y[a1][a2][1][3] * lhs_y[b1][b2][4][1]
		- lhs_y[a1][a2][2][3] * lhs_y[b1][b2][4][2]
		- lhs_y[a1][a2][3][3] * lhs_y[b1][b2][4][3]
		- lhs_y[a1][a2][4][3] * lhs_y[b1][b2][4][4];
	lhs_y[c1][c2][4][4] = lhs_y[c1][c2][4][4] - lhs_y[a1][a2][0][4] * lhs_y[b1][b2][4][0]
		- lhs_y[a1][a2][1][4] * lhs_y[b1][b2][4][1]
		- lhs_y[a1][a2][2][4] * lhs_y[b1][b2][4][2]
		- lhs_y[a1][a2][3][4] * lhs_y[b1][b2][4][3]
		- lhs_y[a1][a2][4][4] * lhs_y[b1][b2][4][4];
}

void binvcrhs_y(int a1, int a2, int b1, int b2, int c2)
{
	double pivot, coeff;

	pivot = 1.00 / lhs_y[a1][a2][0][0];
	lhs_y[a1][a2][1][0] = lhs_y[a1][a2][1][0] * pivot;
	lhs_y[a1][a2][2][0] = lhs_y[a1][a2][2][0] * pivot;
	lhs_y[a1][a2][3][0] = lhs_y[a1][a2][3][0] * pivot;
	lhs_y[a1][a2][4][0] = lhs_y[a1][a2][4][0] * pivot;
	lhs_y[b1][b2][0][0] = lhs_y[b1][b2][0][0] * pivot;
	lhs_y[b1][b2][1][0] = lhs_y[b1][b2][1][0] * pivot;
	lhs_y[b1][b2][2][0] = lhs_y[b1][b2][2][0] * pivot;
	lhs_y[b1][b2][3][0] = lhs_y[b1][b2][3][0] * pivot;
	lhs_y[b1][b2][4][0] = lhs_y[b1][b2][4][0] * pivot;
	rhs_y[c2][0] = rhs_y[c2][0] * pivot;

	coeff = lhs_y[a1][a2][0][1];
	lhs_y[a1][a2][1][1] = lhs_y[a1][a2][1][1] - coeff*lhs_y[a1][a2][1][0];
	lhs_y[a1][a2][2][1] = lhs_y[a1][a2][2][1] - coeff*lhs_y[a1][a2][2][0];
	lhs_y[a1][a2][3][1] = lhs_y[a1][a2][3][1] - coeff*lhs_y[a1][a2][3][0];
	lhs_y[a1][a2][4][1] = lhs_y[a1][a2][4][1] - coeff*lhs_y[a1][a2][4][0];
	lhs_y[b1][b2][0][1] = lhs_y[b1][b2][0][1] - coeff*lhs_y[b1][b2][0][0];
	lhs_y[b1][b2][1][1] = lhs_y[b1][b2][1][1] - coeff*lhs_y[b1][b2][1][0];
	lhs_y[b1][b2][2][1] = lhs_y[b1][b2][2][1] - coeff*lhs_y[b1][b2][2][0];
	lhs_y[b1][b2][3][1] = lhs_y[b1][b2][3][1] - coeff*lhs_y[b1][b2][3][0];
	lhs_y[b1][b2][4][1] = lhs_y[b1][b2][4][1] - coeff*lhs_y[b1][b2][4][0];
	rhs_y[c2][1] = rhs_y[c2][1] - coeff*rhs_y[c2][0];

	coeff = lhs_y[a1][a2][0][2];
	lhs_y[a1][a2][1][2] = lhs_y[a1][a2][1][2] - coeff*lhs_y[a1][a2][1][0];
	lhs_y[a1][a2][2][2] = lhs_y[a1][a2][2][2] - coeff*lhs_y[a1][a2][2][0];
	lhs_y[a1][a2][3][2] = lhs_y[a1][a2][3][2] - coeff*lhs_y[a1][a2][3][0];
	lhs_y[a1][a2][4][2] = lhs_y[a1][a2][4][2] - coeff*lhs_y[a1][a2][4][0];
	lhs_y[b1][b2][0][2] = lhs_y[b1][b2][0][2] - coeff*lhs_y[b1][b2][0][0];
	lhs_y[b1][b2][1][2] = lhs_y[b1][b2][1][2] - coeff*lhs_y[b1][b2][1][0];
	lhs_y[b1][b2][2][2] = lhs_y[b1][b2][2][2] - coeff*lhs_y[b1][b2][2][0];
	lhs_y[b1][b2][3][2] = lhs_y[b1][b2][3][2] - coeff*lhs_y[b1][b2][3][0];
	lhs_y[b1][b2][4][2] = lhs_y[b1][b2][4][2] - coeff*lhs_y[b1][b2][4][0];
	rhs_y[c2][2] = rhs_y[c2][2] - coeff*rhs_y[c2][0];

	coeff = lhs_y[a1][a2][0][3];
	lhs_y[a1][a2][1][3] = lhs_y[a1][a2][1][3] - coeff*lhs_y[a1][a2][1][0];
	lhs_y[a1][a2][2][3] = lhs_y[a1][a2][2][3] - coeff*lhs_y[a1][a2][2][0];
	lhs_y[a1][a2][3][3] = lhs_y[a1][a2][3][3] - coeff*lhs_y[a1][a2][3][0];
	lhs_y[a1][a2][4][3] = lhs_y[a1][a2][4][3] - coeff*lhs_y[a1][a2][4][0];
	lhs_y[b1][b2][0][3] = lhs_y[b1][b2][0][3] - coeff*lhs_y[b1][b2][0][0];
	lhs_y[b1][b2][1][3] = lhs_y[b1][b2][1][3] - coeff*lhs_y[b1][b2][1][0];
	lhs_y[b1][b2][2][3] = lhs_y[b1][b2][2][3] - coeff*lhs_y[b1][b2][2][0];
	lhs_y[b1][b2][3][3] = lhs_y[b1][b2][3][3] - coeff*lhs_y[b1][b2][3][0];
	lhs_y[b1][b2][4][3] = lhs_y[b1][b2][4][3] - coeff*lhs_y[b1][b2][4][0];
	rhs_y[c2][3] = rhs_y[c2][3] - coeff*rhs_y[c2][0];

	coeff = lhs_y[a1][a2][0][4];
	lhs_y[a1][a2][1][4] = lhs_y[a1][a2][1][4] - coeff*lhs_y[a1][a2][1][0];
	lhs_y[a1][a2][2][4] = lhs_y[a1][a2][2][4] - coeff*lhs_y[a1][a2][2][0];
	lhs_y[a1][a2][3][4] = lhs_y[a1][a2][3][4] - coeff*lhs_y[a1][a2][3][0];
	lhs_y[a1][a2][4][4] = lhs_y[a1][a2][4][4] - coeff*lhs_y[a1][a2][4][0];
	lhs_y[b1][b2][0][4] = lhs_y[b1][b2][0][4] - coeff*lhs_y[b1][b2][0][0];
	lhs_y[b1][b2][1][4] = lhs_y[b1][b2][1][4] - coeff*lhs_y[b1][b2][1][0];
	lhs_y[b1][b2][2][4] = lhs_y[b1][b2][2][4] - coeff*lhs_y[b1][b2][2][0];
	lhs_y[b1][b2][3][4] = lhs_y[b1][b2][3][4] - coeff*lhs_y[b1][b2][3][0];
	lhs_y[b1][b2][4][4] = lhs_y[b1][b2][4][4] - coeff*lhs_y[b1][b2][4][0];
	rhs_y[c2][4] = rhs_y[c2][4] - coeff*rhs_y[c2][0];


	pivot = 1.00 / lhs_y[a1][a2][1][1];
	lhs_y[a1][a2][2][1] = lhs_y[a1][a2][2][1] * pivot;
	lhs_y[a1][a2][3][1] = lhs_y[a1][a2][3][1] * pivot;
	lhs_y[a1][a2][4][1] = lhs_y[a1][a2][4][1] * pivot;
	lhs_y[b1][b2][0][1] = lhs_y[b1][b2][0][1] * pivot;
	lhs_y[b1][b2][1][1] = lhs_y[b1][b2][1][1] * pivot;
	lhs_y[b1][b2][2][1] = lhs_y[b1][b2][2][1] * pivot;
	lhs_y[b1][b2][3][1] = lhs_y[b1][b2][3][1] * pivot;
	lhs_y[b1][b2][4][1] = lhs_y[b1][b2][4][1] * pivot;
	rhs_y[c2][1] = rhs_y[c2][1] * pivot;

	coeff = lhs_y[a1][a2][1][0];
	lhs_y[a1][a2][2][0] = lhs_y[a1][a2][2][0] - coeff*lhs_y[a1][a2][2][1];
	lhs_y[a1][a2][3][0] = lhs_y[a1][a2][3][0] - coeff*lhs_y[a1][a2][3][1];
	lhs_y[a1][a2][4][0] = lhs_y[a1][a2][4][0] - coeff*lhs_y[a1][a2][4][1];
	lhs_y[b1][b2][0][0] = lhs_y[b1][b2][0][0] - coeff*lhs_y[b1][b2][0][1];
	lhs_y[b1][b2][1][0] = lhs_y[b1][b2][1][0] - coeff*lhs_y[b1][b2][1][1];
	lhs_y[b1][b2][2][0] = lhs_y[b1][b2][2][0] - coeff*lhs_y[b1][b2][2][1];
	lhs_y[b1][b2][3][0] = lhs_y[b1][b2][3][0] - coeff*lhs_y[b1][b2][3][1];
	lhs_y[b1][b2][4][0] = lhs_y[b1][b2][4][0] - coeff*lhs_y[b1][b2][4][1];
	rhs_y[c2][0] = rhs_y[c2][0] - coeff*rhs_y[c2][1];

	coeff = lhs_y[a1][a2][1][2];
	lhs_y[a1][a2][2][2] = lhs_y[a1][a2][2][2] - coeff*lhs_y[a1][a2][2][1];
	lhs_y[a1][a2][3][2] = lhs_y[a1][a2][3][2] - coeff*lhs_y[a1][a2][3][1];
	lhs_y[a1][a2][4][2] = lhs_y[a1][a2][4][2] - coeff*lhs_y[a1][a2][4][1];
	lhs_y[b1][b2][0][2] = lhs_y[b1][b2][0][2] - coeff*lhs_y[b1][b2][0][1];
	lhs_y[b1][b2][1][2] = lhs_y[b1][b2][1][2] - coeff*lhs_y[b1][b2][1][1];
	lhs_y[b1][b2][2][2] = lhs_y[b1][b2][2][2] - coeff*lhs_y[b1][b2][2][1];
	lhs_y[b1][b2][3][2] = lhs_y[b1][b2][3][2] - coeff*lhs_y[b1][b2][3][1];
	lhs_y[b1][b2][4][2] = lhs_y[b1][b2][4][2] - coeff*lhs_y[b1][b2][4][1];
	rhs_y[c2][2] = rhs_y[c2][2] - coeff*rhs_y[c2][1];

	coeff = lhs_y[a1][a2][1][3];
	lhs_y[a1][a2][2][3] = lhs_y[a1][a2][2][3] - coeff*lhs_y[a1][a2][2][1];
	lhs_y[a1][a2][3][3] = lhs_y[a1][a2][3][3] - coeff*lhs_y[a1][a2][3][1];
	lhs_y[a1][a2][4][3] = lhs_y[a1][a2][4][3] - coeff*lhs_y[a1][a2][4][1];
	lhs_y[b1][b2][0][3] = lhs_y[b1][b2][0][3] - coeff*lhs_y[b1][b2][0][1];
	lhs_y[b1][b2][1][3] = lhs_y[b1][b2][1][3] - coeff*lhs_y[b1][b2][1][1];
	lhs_y[b1][b2][2][3] = lhs_y[b1][b2][2][3] - coeff*lhs_y[b1][b2][2][1];
	lhs_y[b1][b2][3][3] = lhs_y[b1][b2][3][3] - coeff*lhs_y[b1][b2][3][1];
	lhs_y[b1][b2][4][3] = lhs_y[b1][b2][4][3] - coeff*lhs_y[b1][b2][4][1];
	rhs_y[c2][3] = rhs_y[c2][3] - coeff*rhs_y[c2][1];

	coeff = lhs_y[a1][a2][1][4];
	lhs_y[a1][a2][2][4] = lhs_y[a1][a2][2][4] - coeff*lhs_y[a1][a2][2][1];
	lhs_y[a1][a2][3][4] = lhs_y[a1][a2][3][4] - coeff*lhs_y[a1][a2][3][1];
	lhs_y[a1][a2][4][4] = lhs_y[a1][a2][4][4] - coeff*lhs_y[a1][a2][4][1];
	lhs_y[b1][b2][0][4] = lhs_y[b1][b2][0][4] - coeff*lhs_y[b1][b2][0][1];
	lhs_y[b1][b2][1][4] = lhs_y[b1][b2][1][4] - coeff*lhs_y[b1][b2][1][1];
	lhs_y[b1][b2][2][4] = lhs_y[b1][b2][2][4] - coeff*lhs_y[b1][b2][2][1];
	lhs_y[b1][b2][3][4] = lhs_y[b1][b2][3][4] - coeff*lhs_y[b1][b2][3][1];
	lhs_y[b1][b2][4][4] = lhs_y[b1][b2][4][4] - coeff*lhs_y[b1][b2][4][1];
	rhs_y[c2][4] = rhs_y[c2][4] - coeff*rhs_y[c2][1];


	pivot = 1.00 / lhs_y[a1][a2][2][2];
	lhs_y[a1][a2][3][2] = lhs_y[a1][a2][3][2] * pivot;
	lhs_y[a1][a2][4][2] = lhs_y[a1][a2][4][2] * pivot;
	lhs_y[b1][b2][0][2] = lhs_y[b1][b2][0][2] * pivot;
	lhs_y[b1][b2][1][2] = lhs_y[b1][b2][1][2] * pivot;
	lhs_y[b1][b2][2][2] = lhs_y[b1][b2][2][2] * pivot;
	lhs_y[b1][b2][3][2] = lhs_y[b1][b2][3][2] * pivot;
	lhs_y[b1][b2][4][2] = lhs_y[b1][b2][4][2] * pivot;
	rhs_y[c2][2] = rhs_y[c2][2] * pivot;

	coeff = lhs_y[a1][a2][2][0];
	lhs_y[a1][a2][3][0] = lhs_y[a1][a2][3][0] - coeff*lhs_y[a1][a2][3][2];
	lhs_y[a1][a2][4][0] = lhs_y[a1][a2][4][0] - coeff*lhs_y[a1][a2][4][2];
	lhs_y[b1][b2][0][0] = lhs_y[b1][b2][0][0] - coeff*lhs_y[b1][b2][0][2];
	lhs_y[b1][b2][1][0] = lhs_y[b1][b2][1][0] - coeff*lhs_y[b1][b2][1][2];
	lhs_y[b1][b2][2][0] = lhs_y[b1][b2][2][0] - coeff*lhs_y[b1][b2][2][2];
	lhs_y[b1][b2][3][0] = lhs_y[b1][b2][3][0] - coeff*lhs_y[b1][b2][3][2];
	lhs_y[b1][b2][4][0] = lhs_y[b1][b2][4][0] - coeff*lhs_y[b1][b2][4][2];
	rhs_y[c2][0] = rhs_y[c2][0] - coeff*rhs_y[c2][2];

	coeff = lhs_y[a1][a2][2][1];
	lhs_y[a1][a2][3][1] = lhs_y[a1][a2][3][1] - coeff*lhs_y[a1][a2][3][2];
	lhs_y[a1][a2][4][1] = lhs_y[a1][a2][4][1] - coeff*lhs_y[a1][a2][4][2];
	lhs_y[b1][b2][0][1] = lhs_y[b1][b2][0][1] - coeff*lhs_y[b1][b2][0][2];
	lhs_y[b1][b2][1][1] = lhs_y[b1][b2][1][1] - coeff*lhs_y[b1][b2][1][2];
	lhs_y[b1][b2][2][1] = lhs_y[b1][b2][2][1] - coeff*lhs_y[b1][b2][2][2];
	lhs_y[b1][b2][3][1] = lhs_y[b1][b2][3][1] - coeff*lhs_y[b1][b2][3][2];
	lhs_y[b1][b2][4][1] = lhs_y[b1][b2][4][1] - coeff*lhs_y[b1][b2][4][2];
	rhs_y[c2][1] = rhs_y[c2][1] - coeff*rhs_y[c2][2];

	coeff = lhs_y[a1][a2][2][3];
	lhs_y[a1][a2][3][3] = lhs_y[a1][a2][3][3] - coeff*lhs_y[a1][a2][3][2];
	lhs_y[a1][a2][4][3] = lhs_y[a1][a2][4][3] - coeff*lhs_y[a1][a2][4][2];
	lhs_y[b1][b2][0][3] = lhs_y[b1][b2][0][3] - coeff*lhs_y[b1][b2][0][2];
	lhs_y[b1][b2][1][3] = lhs_y[b1][b2][1][3] - coeff*lhs_y[b1][b2][1][2];
	lhs_y[b1][b2][2][3] = lhs_y[b1][b2][2][3] - coeff*lhs_y[b1][b2][2][2];
	lhs_y[b1][b2][3][3] = lhs_y[b1][b2][3][3] - coeff*lhs_y[b1][b2][3][2];
	lhs_y[b1][b2][4][3] = lhs_y[b1][b2][4][3] - coeff*lhs_y[b1][b2][4][2];
	rhs_y[c2][3] = rhs_y[c2][3] - coeff*rhs_y[c2][2];

	coeff = lhs_y[a1][a2][2][4];
	lhs_y[a1][a2][3][4] = lhs_y[a1][a2][3][4] - coeff*lhs_y[a1][a2][3][2];
	lhs_y[a1][a2][4][4] = lhs_y[a1][a2][4][4] - coeff*lhs_y[a1][a2][4][2];
	lhs_y[b1][b2][0][4] = lhs_y[b1][b2][0][4] - coeff*lhs_y[b1][b2][0][2];
	lhs_y[b1][b2][1][4] = lhs_y[b1][b2][1][4] - coeff*lhs_y[b1][b2][1][2];
	lhs_y[b1][b2][2][4] = lhs_y[b1][b2][2][4] - coeff*lhs_y[b1][b2][2][2];
	lhs_y[b1][b2][3][4] = lhs_y[b1][b2][3][4] - coeff*lhs_y[b1][b2][3][2];
	lhs_y[b1][b2][4][4] = lhs_y[b1][b2][4][4] - coeff*lhs_y[b1][b2][4][2];
	rhs_y[c2][4] = rhs_y[c2][4] - coeff*rhs_y[c2][2];


	pivot = 1.00 / lhs_y[a1][a2][3][3];
	lhs_y[a1][a2][4][3] = lhs_y[a1][a2][4][3] * pivot;
	lhs_y[b1][b2][0][3] = lhs_y[b1][b2][0][3] * pivot;
	lhs_y[b1][b2][1][3] = lhs_y[b1][b2][1][3] * pivot;
	lhs_y[b1][b2][2][3] = lhs_y[b1][b2][2][3] * pivot;
	lhs_y[b1][b2][3][3] = lhs_y[b1][b2][3][3] * pivot;
	lhs_y[b1][b2][4][3] = lhs_y[b1][b2][4][3] * pivot;
	rhs_y[c2][3] = rhs_y[c2][3] * pivot;

	coeff = lhs_y[a1][a2][3][0];
	lhs_y[a1][a2][4][0] = lhs_y[a1][a2][4][0] - coeff*lhs_y[a1][a2][4][3];
	lhs_y[b1][b2][0][0] = lhs_y[b1][b2][0][0] - coeff*lhs_y[b1][b2][0][3];
	lhs_y[b1][b2][1][0] = lhs_y[b1][b2][1][0] - coeff*lhs_y[b1][b2][1][3];
	lhs_y[b1][b2][2][0] = lhs_y[b1][b2][2][0] - coeff*lhs_y[b1][b2][2][3];
	lhs_y[b1][b2][3][0] = lhs_y[b1][b2][3][0] - coeff*lhs_y[b1][b2][3][3];
	lhs_y[b1][b2][4][0] = lhs_y[b1][b2][4][0] - coeff*lhs_y[b1][b2][4][3];
	rhs_y[c2][0] = rhs_y[c2][0] - coeff*rhs_y[c2][3];

	coeff = lhs_y[a1][a2][3][1];
	lhs_y[a1][a2][4][1] = lhs_y[a1][a2][4][1] - coeff*lhs_y[a1][a2][4][3];
	lhs_y[b1][b2][0][1] = lhs_y[b1][b2][0][1] - coeff*lhs_y[b1][b2][0][3];
	lhs_y[b1][b2][1][1] = lhs_y[b1][b2][1][1] - coeff*lhs_y[b1][b2][1][3];
	lhs_y[b1][b2][2][1] = lhs_y[b1][b2][2][1] - coeff*lhs_y[b1][b2][2][3];
	lhs_y[b1][b2][3][1] = lhs_y[b1][b2][3][1] - coeff*lhs_y[b1][b2][3][3];
	lhs_y[b1][b2][4][1] = lhs_y[b1][b2][4][1] - coeff*lhs_y[b1][b2][4][3];
	rhs_y[c2][1] = rhs_y[c2][1] - coeff*rhs_y[c2][3];

	coeff = lhs_y[a1][a2][3][2];
	lhs_y[a1][a2][4][2] = lhs_y[a1][a2][4][2] - coeff*lhs_y[a1][a2][4][3];
	lhs_y[b1][b2][0][2] = lhs_y[b1][b2][0][2] - coeff*lhs_y[b1][b2][0][3];
	lhs_y[b1][b2][1][2] = lhs_y[b1][b2][1][2] - coeff*lhs_y[b1][b2][1][3];
	lhs_y[b1][b2][2][2] = lhs_y[b1][b2][2][2] - coeff*lhs_y[b1][b2][2][3];
	lhs_y[b1][b2][3][2] = lhs_y[b1][b2][3][2] - coeff*lhs_y[b1][b2][3][3];
	lhs_y[b1][b2][4][2] = lhs_y[b1][b2][4][2] - coeff*lhs_y[b1][b2][4][3];
	rhs_y[c2][2] = rhs_y[c2][2] - coeff*rhs_y[c2][3];

	coeff = lhs_y[a1][a2][3][4];
	lhs_y[a1][a2][4][4] = lhs_y[a1][a2][4][4] - coeff*lhs_y[a1][a2][4][3];
	lhs_y[b1][b2][0][4] = lhs_y[b1][b2][0][4] - coeff*lhs_y[b1][b2][0][3];
	lhs_y[b1][b2][1][4] = lhs_y[b1][b2][1][4] - coeff*lhs_y[b1][b2][1][3];
	lhs_y[b1][b2][2][4] = lhs_y[b1][b2][2][4] - coeff*lhs_y[b1][b2][2][3];
	lhs_y[b1][b2][3][4] = lhs_y[b1][b2][3][4] - coeff*lhs_y[b1][b2][3][3];
	lhs_y[b1][b2][4][4] = lhs_y[b1][b2][4][4] - coeff*lhs_y[b1][b2][4][3];
	rhs_y[c2][4] = rhs_y[c2][4] - coeff*rhs_y[c2][3];


	pivot = 1.00 / lhs_y[a1][a2][4][4];
	lhs_y[b1][b2][0][4] = lhs_y[b1][b2][0][4] * pivot;
	lhs_y[b1][b2][1][4] = lhs_y[b1][b2][1][4] * pivot;
	lhs_y[b1][b2][2][4] = lhs_y[b1][b2][2][4] * pivot;
	lhs_y[b1][b2][3][4] = lhs_y[b1][b2][3][4] * pivot;
	lhs_y[b1][b2][4][4] = lhs_y[b1][b2][4][4] * pivot;
	rhs_y[c2][4] = rhs_y[c2][4] * pivot;

	coeff = lhs_y[a1][a2][4][0];
	lhs_y[b1][b2][0][0] = lhs_y[b1][b2][0][0] - coeff*lhs_y[b1][b2][0][4];
	lhs_y[b1][b2][1][0] = lhs_y[b1][b2][1][0] - coeff*lhs_y[b1][b2][1][4];
	lhs_y[b1][b2][2][0] = lhs_y[b1][b2][2][0] - coeff*lhs_y[b1][b2][2][4];
	lhs_y[b1][b2][3][0] = lhs_y[b1][b2][3][0] - coeff*lhs_y[b1][b2][3][4];
	lhs_y[b1][b2][4][0] = lhs_y[b1][b2][4][0] - coeff*lhs_y[b1][b2][4][4];
	rhs_y[c2][0] = rhs_y[c2][0] - coeff*rhs_y[c2][4];

	coeff = lhs_y[a1][a2][4][1];
	lhs_y[b1][b2][0][1] = lhs_y[b1][b2][0][1] - coeff*lhs_y[b1][b2][0][4];
	lhs_y[b1][b2][1][1] = lhs_y[b1][b2][1][1] - coeff*lhs_y[b1][b2][1][4];
	lhs_y[b1][b2][2][1] = lhs_y[b1][b2][2][1] - coeff*lhs_y[b1][b2][2][4];
	lhs_y[b1][b2][3][1] = lhs_y[b1][b2][3][1] - coeff*lhs_y[b1][b2][3][4];
	lhs_y[b1][b2][4][1] = lhs_y[b1][b2][4][1] - coeff*lhs_y[b1][b2][4][4];
	rhs_y[c2][1] = rhs_y[c2][1] - coeff*rhs_y[c2][4];

	coeff = lhs_y[a1][a2][4][2];
	lhs_y[b1][b2][0][2] = lhs_y[b1][b2][0][2] - coeff*lhs_y[b1][b2][0][4];
	lhs_y[b1][b2][1][2] = lhs_y[b1][b2][1][2] - coeff*lhs_y[b1][b2][1][4];
	lhs_y[b1][b2][2][2] = lhs_y[b1][b2][2][2] - coeff*lhs_y[b1][b2][2][4];
	lhs_y[b1][b2][3][2] = lhs_y[b1][b2][3][2] - coeff*lhs_y[b1][b2][3][4];
	lhs_y[b1][b2][4][2] = lhs_y[b1][b2][4][2] - coeff*lhs_y[b1][b2][4][4];
	rhs_y[c2][2] = rhs_y[c2][2] - coeff*rhs_y[c2][4];

	coeff = lhs_y[a1][a2][4][3];
	lhs_y[b1][b2][0][3] = lhs_y[b1][b2][0][3] - coeff*lhs_y[b1][b2][0][4];
	lhs_y[b1][b2][1][3] = lhs_y[b1][b2][1][3] - coeff*lhs_y[b1][b2][1][4];
	lhs_y[b1][b2][2][3] = lhs_y[b1][b2][2][3] - coeff*lhs_y[b1][b2][2][4];
	lhs_y[b1][b2][3][3] = lhs_y[b1][b2][3][3] - coeff*lhs_y[b1][b2][3][4];
	lhs_y[b1][b2][4][3] = lhs_y[b1][b2][4][3] - coeff*lhs_y[b1][b2][4][4];
	rhs_y[c2][3] = rhs_y[c2][3] - coeff*rhs_y[c2][4];
}

void binvrhs_y(int a1, int a2, int b2)
{
	double pivot, coeff;

	pivot = 1.00 / lhs_y[a1][a2][0][0];
	lhs_y[a1][a2][1][0] = lhs_y[a1][a2][1][0] * pivot;
	lhs_y[a1][a2][2][0] = lhs_y[a1][a2][2][0] * pivot;
	lhs_y[a1][a2][3][0] = lhs_y[a1][a2][3][0] * pivot;
	lhs_y[a1][a2][4][0] = lhs_y[a1][a2][4][0] * pivot;
	rhs_y[b2][0] = rhs_y[b2][0] * pivot;

	coeff = lhs_y[a1][a2][0][1];
	lhs_y[a1][a2][1][1] = lhs_y[a1][a2][1][1] - coeff*lhs_y[a1][a2][1][0];
	lhs_y[a1][a2][2][1] = lhs_y[a1][a2][2][1] - coeff*lhs_y[a1][a2][2][0];
	lhs_y[a1][a2][3][1] = lhs_y[a1][a2][3][1] - coeff*lhs_y[a1][a2][3][0];
	lhs_y[a1][a2][4][1] = lhs_y[a1][a2][4][1] - coeff*lhs_y[a1][a2][4][0];
	rhs_y[b2][1] = rhs_y[b2][1] - coeff*rhs_y[b2][0];

	coeff = lhs_y[a1][a2][0][2];
	lhs_y[a1][a2][1][2] = lhs_y[a1][a2][1][2] - coeff*lhs_y[a1][a2][1][0];
	lhs_y[a1][a2][2][2] = lhs_y[a1][a2][2][2] - coeff*lhs_y[a1][a2][2][0];
	lhs_y[a1][a2][3][2] = lhs_y[a1][a2][3][2] - coeff*lhs_y[a1][a2][3][0];
	lhs_y[a1][a2][4][2] = lhs_y[a1][a2][4][2] - coeff*lhs_y[a1][a2][4][0];
	rhs_y[b2][2] = rhs_y[b2][2] - coeff*rhs_y[b2][0];

	coeff = lhs_y[a1][a2][0][3];
	lhs_y[a1][a2][1][3] = lhs_y[a1][a2][1][3] - coeff*lhs_y[a1][a2][1][0];
	lhs_y[a1][a2][2][3] = lhs_y[a1][a2][2][3] - coeff*lhs_y[a1][a2][2][0];
	lhs_y[a1][a2][3][3] = lhs_y[a1][a2][3][3] - coeff*lhs_y[a1][a2][3][0];
	lhs_y[a1][a2][4][3] = lhs_y[a1][a2][4][3] - coeff*lhs_y[a1][a2][4][0];
	rhs_y[b2][3] = rhs_y[b2][3] - coeff*rhs_y[b2][0];

	coeff = lhs_y[a1][a2][0][4];
	lhs_y[a1][a2][1][4] = lhs_y[a1][a2][1][4] - coeff*lhs_y[a1][a2][1][0];
	lhs_y[a1][a2][2][4] = lhs_y[a1][a2][2][4] - coeff*lhs_y[a1][a2][2][0];
	lhs_y[a1][a2][3][4] = lhs_y[a1][a2][3][4] - coeff*lhs_y[a1][a2][3][0];
	lhs_y[a1][a2][4][4] = lhs_y[a1][a2][4][4] - coeff*lhs_y[a1][a2][4][0];
	rhs_y[b2][4] = rhs_y[b2][4] - coeff*rhs_y[b2][0];


	pivot = 1.00 / lhs_y[a1][a2][1][1];
	lhs_y[a1][a2][2][1] = lhs_y[a1][a2][2][1] * pivot;
	lhs_y[a1][a2][3][1] = lhs_y[a1][a2][3][1] * pivot;
	lhs_y[a1][a2][4][1] = lhs_y[a1][a2][4][1] * pivot;
	rhs_y[b2][1] = rhs_y[b2][1] * pivot;

	coeff = lhs_y[a1][a2][1][0];
	lhs_y[a1][a2][2][0] = lhs_y[a1][a2][2][0] - coeff*lhs_y[a1][a2][2][1];
	lhs_y[a1][a2][3][0] = lhs_y[a1][a2][3][0] - coeff*lhs_y[a1][a2][3][1];
	lhs_y[a1][a2][4][0] = lhs_y[a1][a2][4][0] - coeff*lhs_y[a1][a2][4][1];
	rhs_y[b2][0] = rhs_y[b2][0] - coeff*rhs_y[b2][1];

	coeff = lhs_y[a1][a2][1][2];
	lhs_y[a1][a2][2][2] = lhs_y[a1][a2][2][2] - coeff*lhs_y[a1][a2][2][1];
	lhs_y[a1][a2][3][2] = lhs_y[a1][a2][3][2] - coeff*lhs_y[a1][a2][3][1];
	lhs_y[a1][a2][4][2] = lhs_y[a1][a2][4][2] - coeff*lhs_y[a1][a2][4][1];
	rhs_y[b2][2] = rhs_y[b2][2] - coeff*rhs_y[b2][1];

	coeff = lhs_y[a1][a2][1][3];
	lhs_y[a1][a2][2][3] = lhs_y[a1][a2][2][3] - coeff*lhs_y[a1][a2][2][1];
	lhs_y[a1][a2][3][3] = lhs_y[a1][a2][3][3] - coeff*lhs_y[a1][a2][3][1];
	lhs_y[a1][a2][4][3] = lhs_y[a1][a2][4][3] - coeff*lhs_y[a1][a2][4][1];
	rhs_y[b2][3] = rhs_y[b2][3] - coeff*rhs_y[b2][1];

	coeff = lhs_y[a1][a2][1][4];
	lhs_y[a1][a2][2][4] = lhs_y[a1][a2][2][4] - coeff*lhs_y[a1][a2][2][1];
	lhs_y[a1][a2][3][4] = lhs_y[a1][a2][3][4] - coeff*lhs_y[a1][a2][3][1];
	lhs_y[a1][a2][4][4] = lhs_y[a1][a2][4][4] - coeff*lhs_y[a1][a2][4][1];
	rhs_y[b2][4] = rhs_y[b2][4] - coeff*rhs_y[b2][1];


	pivot = 1.00 / lhs_y[a1][a2][2][2];
	lhs_y[a1][a2][3][2] = lhs_y[a1][a2][3][2] * pivot;
	lhs_y[a1][a2][4][2] = lhs_y[a1][a2][4][2] * pivot;
	rhs_y[b2][2] = rhs_y[b2][2] * pivot;

	coeff = lhs_y[a1][a2][2][0];
	lhs_y[a1][a2][3][0] = lhs_y[a1][a2][3][0] - coeff*lhs_y[a1][a2][3][2];
	lhs_y[a1][a2][4][0] = lhs_y[a1][a2][4][0] - coeff*lhs_y[a1][a2][4][2];
	rhs_y[b2][0] = rhs_y[b2][0] - coeff*rhs_y[b2][2];

	coeff = lhs_y[a1][a2][2][1];
	lhs_y[a1][a2][3][1] = lhs_y[a1][a2][3][1] - coeff*lhs_y[a1][a2][3][2];
	lhs_y[a1][a2][4][1] = lhs_y[a1][a2][4][1] - coeff*lhs_y[a1][a2][4][2];
	rhs_y[b2][1] = rhs_y[b2][1] - coeff*rhs_y[b2][2];

	coeff = lhs_y[a1][a2][2][3];
	lhs_y[a1][a2][3][3] = lhs_y[a1][a2][3][3] - coeff*lhs_y[a1][a2][3][2];
	lhs_y[a1][a2][4][3] = lhs_y[a1][a2][4][3] - coeff*lhs_y[a1][a2][4][2];
	rhs_y[b2][3] = rhs_y[b2][3] - coeff*rhs_y[b2][2];

	coeff = lhs_y[a1][a2][2][4];
	lhs_y[a1][a2][3][4] = lhs_y[a1][a2][3][4] - coeff*lhs_y[a1][a2][3][2];
	lhs_y[a1][a2][4][4] = lhs_y[a1][a2][4][4] - coeff*lhs_y[a1][a2][4][2];
	rhs_y[b2][4] = rhs_y[b2][4] - coeff*rhs_y[b2][2];


	pivot = 1.00 / lhs_y[a1][a2][3][3];
	lhs_y[a1][a2][4][3] = lhs_y[a1][a2][4][3] * pivot;
	rhs_y[b2][3] = rhs_y[b2][3] * pivot;

	coeff = lhs_y[a1][a2][3][0];
	lhs_y[a1][a2][4][0] = lhs_y[a1][a2][4][0] - coeff*lhs_y[a1][a2][4][3];
	rhs_y[b2][0] = rhs_y[b2][0] - coeff*rhs_y[b2][3];

	coeff = lhs_y[a1][a2][3][1];
	lhs_y[a1][a2][4][1] = lhs_y[a1][a2][4][1] - coeff*lhs_y[a1][a2][4][3];
	rhs_y[b2][1] = rhs_y[b2][1] - coeff*rhs_y[b2][3];

	coeff = lhs_y[a1][a2][3][2];
	lhs_y[a1][a2][4][2] = lhs_y[a1][a2][4][2] - coeff*lhs_y[a1][a2][4][3];
	rhs_y[b2][2] = rhs_y[b2][2] - coeff*rhs_y[b2][3];

	coeff = lhs_y[a1][a2][3][4];
	lhs_y[a1][a2][4][4] = lhs_y[a1][a2][4][4] - coeff*lhs_y[a1][a2][4][3];
	rhs_y[b2][4] = rhs_y[b2][4] - coeff*rhs_y[b2][3];


	pivot = 1.00 / lhs_y[a1][a2][4][4];
	rhs_y[b2][4] = rhs_y[b2][4] * pivot;

	coeff = lhs_y[a1][a2][4][0];
	rhs_y[b2][0] = rhs_y[b2][0] - coeff*rhs_y[b2][4];

	coeff = lhs_y[a1][a2][4][1];
	rhs_y[b2][1] = rhs_y[b2][1] - coeff*rhs_y[b2][4];

	coeff = lhs_y[a1][a2][4][2];
	rhs_y[b2][2] = rhs_y[b2][2] - coeff*rhs_y[b2][4];

	coeff = lhs_y[a1][a2][4][3];
	rhs_y[b2][3] = rhs_y[b2][3] - coeff*rhs_y[b2][4];
}
