#include <cdvmh_helpers.h>

#define DVM0C(n) ((DvmType)(n))

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

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//
//  header.h
//
//---------------------------------------------------------------------
//---------------------------------------------------------------------
 
//---------------------------------------------------------------------
// The following include file is generated automatically by the
// "setparams" utility. It defines 
//      maxcells:      the square root of the maximum number of processors
//      problem_size:  12, 64, 102, 162 (for class T, A, B, C)
//      dt_default:    default time step for this problem size if no
//                     config file
//      niter_default: default number of iterations for this problem size
//---------------------------------------------------------------------

#include "npbparams.h"
#include "../common/type.h"

#define AA            0
#define BB            1
#define CC            2
#define BLOCK_SIZE    5

/* common /global/ */
extern double elapsed_time;
extern int grid_points[3];
extern logical timeron;

/* common /constants/ */
extern double tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3, 
              dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, 
              dy5, dz1, dz2, dz3, dz4, dz5, dssp, dt, 
              ce[5][13], dxmax, dymax, dzmax, xxcon1, xxcon2, 
              xxcon3, xxcon4, xxcon5, dx1tx1, dx2tx1, dx3tx1,
              dx4tx1, dx5tx1, yycon1, yycon2, yycon3, yycon4,
              yycon5, dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1,
              zzcon1, zzcon2, zzcon3, zzcon4, zzcon5, dz1tz1, 
              dz2tz1, dz3tz1, dz4tz1, dz5tz1, dnxm1, dnym1, 
              dnzm1, c1c2, c1c5, c3c4, c1345, conz1, c1, c2, 
              c3, c4, c5, c4dssp, c5dssp, dtdssp, dttx1,
              dttx2, dtty1, dtty2, dttz1, dttz2, c2dttx1, 
              c2dtty1, c2dttz1, comz1, comz4, comz5, comz6, 
              c3c4tx3, c3c4ty3, c3c4tz3, c2iv, con43, con16;

#define IMAX      PROBLEM_SIZE
#define JMAX      PROBLEM_SIZE
#define KMAX      PROBLEM_SIZE
#define IMAXP     IMAX/2*2
#define JMAXP     JMAX/2*2


// to improve cache performance, grid dimensions padded by 1 
// for even number sizes only.
/* common /fields/ */

extern DvmType u[64];

extern DvmType us[64];
extern DvmType vs[64];
extern DvmType ws[64];
extern DvmType qs[64];
extern DvmType rho_i[64];
extern DvmType square[64];
extern DvmType forcing[64];



//#pragma dvm array
extern double rhs    [KMAX][JMAXP+1][IMAXP+1][5];

/* common /work_1d/ */
extern double cuf[PROBLEM_SIZE+1];
extern double q  [PROBLEM_SIZE+1];
extern double ue [PROBLEM_SIZE+1][5];
extern double buf[PROBLEM_SIZE+1][5];
      

//-----------------------------------------------------------------------
// Timer constants
//-----------------------------------------------------------------------
#define t_total     1
#define t_rhsx      2
#define t_rhsy      3
#define t_rhsz      4
#define t_rhs       5
#define t_xsolve    6
#define t_ysolve    7
#define t_zsolve    8
#define t_rdis1     9
#define t_rdis2     10
#define t_add       11
#define t_last      11


void initialize();
void lhsinit(double lhs[][3][5][5], int size);
void exact_solution(double xi, double eta, double zeta, double dtemp[5]);
void exact_rhs();
void set_constants();
void adi();
void compute_rhs();
void x_solve();
void y_solve();

//void matvec_sub(double ablock[5][5], double avec[5], double bvec[5]);
void matvec_sub(int a1, int a2, int b1, int b2, int b3, int c1, int c2, int c3 );

//void matmul_sub(double ablock[5][5], double bblock[5][5], double cblock[5][5]);
void matmul_sub(int a1, int a2, int b1, int b2, int c1, int c2);

//void binvcrhs(double lhs[5][5], double c[5][5], double r[5]);
void binvcrhs(int a1, int a2, int b1, int b2, int c1, int c2, int c3);

//void binvrhs(double lhs[5][5], double r[5]);
void binvrhs(int a1, int a2, int b1, int b2, int b3);

void z_solve();
void add();
void error_norm(double rms[5]);
void rhs_norm(double rms[5]);
void verify(int no_time_steps, char *class, logical *verified);


#include "work_lhs.h"
#include "../common/timers.h"

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
  int i, j, k, m, n, ksize;
  int z, a = -1, b = 2;
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  if (timeron) timer_start(t_zsolve);

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  // This function computes the left hand side for the three z-factors   
  //---------------------------------------------------------------------

  ksize = grid_points[2]-1;

  //---------------------------------------------------------------------
  // Compute the indices for storing the block-diagonal matrix;
  // determine c (labeled f) and s jacobians
  //---------------------------------------------------------------------
  /*for (j = 1; j <= grid_points[1]-2; j++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (k = 0; k <= ksize; k++) {
        tmp1 = 1.0 / u[k][j][i][0];
        tmp2 = tmp1 * tmp1;
        tmp3 = tmp1 * tmp2;

        fjac[k][0][0] = 0.0;
        fjac[k][1][0] = 0.0;
        fjac[k][2][0] = 0.0;
        fjac[k][3][0] = 1.0;
        fjac[k][4][0] = 0.0;

        fjac[k][0][1] = - ( u[k][j][i][1]*u[k][j][i][3] ) * tmp2;
        fjac[k][1][1] = u[k][j][i][3] * tmp1;
        fjac[k][2][1] = 0.0;
        fjac[k][3][1] = u[k][j][i][1] * tmp1;
        fjac[k][4][1] = 0.0;

        fjac[k][0][2] = - ( u[k][j][i][2]*u[k][j][i][3] ) * tmp2;
        fjac[k][1][2] = 0.0;
        fjac[k][2][2] = u[k][j][i][3] * tmp1;
        fjac[k][3][2] = u[k][j][i][2] * tmp1;
        fjac[k][4][2] = 0.0;

        fjac[k][0][3] = - (u[k][j][i][3]*u[k][j][i][3] * tmp2 ) 
          + c2 * qs[k][j][i];
        fjac[k][1][3] = - c2 *  u[k][j][i][1] * tmp1;
        fjac[k][2][3] = - c2 *  u[k][j][i][2] * tmp1;
        fjac[k][3][3] = ( 2.0 - c2 ) *  u[k][j][i][3] * tmp1;
        fjac[k][4][3] = c2;

        fjac[k][0][4] = ( c2 * 2.0 * square[k][j][i] - c1 * u[k][j][i][4] )
          * u[k][j][i][3] * tmp2;
        fjac[k][1][4] = - c2 * ( u[k][j][i][1]*u[k][j][i][3] ) * tmp2;
        fjac[k][2][4] = - c2 * ( u[k][j][i][2]*u[k][j][i][3] ) * tmp2;
        fjac[k][3][4] = c1 * ( u[k][j][i][4] * tmp1 )
          - c2 * ( qs[k][j][i] + u[k][j][i][3]*u[k][j][i][3] * tmp2 );
        fjac[k][4][4] = c1 * u[k][j][i][3] * tmp1;

        njac[k][0][0] = 0.0;
        njac[k][1][0] = 0.0;
        njac[k][2][0] = 0.0;
        njac[k][3][0] = 0.0;
        njac[k][4][0] = 0.0;

        njac[k][0][1] = - c3c4 * tmp2 * u[k][j][i][1];
        njac[k][1][1] =   c3c4 * tmp1;
        njac[k][2][1] =   0.0;
        njac[k][3][1] =   0.0;
        njac[k][4][1] =   0.0;

        njac[k][0][2] = - c3c4 * tmp2 * u[k][j][i][2];
        njac[k][1][2] =   0.0;
        njac[k][2][2] =   c3c4 * tmp1;
        njac[k][3][2] =   0.0;
        njac[k][4][2] =   0.0;

        njac[k][0][3] = - con43 * c3c4 * tmp2 * u[k][j][i][3];
        njac[k][1][3] =   0.0;
        njac[k][2][3] =   0.0;
        njac[k][3][3] =   con43 * c3 * c4 * tmp1;
        njac[k][4][3] =   0.0;

        njac[k][0][4] = - (  c3c4
            - c1345 ) * tmp3 * (u[k][j][i][1]*u[k][j][i][1])
          - ( c3c4 - c1345 ) * tmp3 * (u[k][j][i][2]*u[k][j][i][2])
          - ( con43 * c3c4
              - c1345 ) * tmp3 * (u[k][j][i][3]*u[k][j][i][3])
          - c1345 * tmp2 * u[k][j][i][4];

        njac[k][1][4] = (  c3c4 - c1345 ) * tmp2 * u[k][j][i][1];
        njac[k][2][4] = (  c3c4 - c1345 ) * tmp2 * u[k][j][i][2];
        njac[k][3][4] = ( con43 * c3c4
            - c1345 ) * tmp2 * u[k][j][i][3];
        njac[k][4][4] = ( c1345 )* tmp1;
      }

      //---------------------------------------------------------------------
      // now jacobians set, so form left hand side in z direction
      //---------------------------------------------------------------------
      lhsinit(lhs, ksize);
      for (k = 1; k <= ksize-1; k++) {
        tmp1 = dt * tz1;
        tmp2 = dt * tz2;

        lhs[k][AA][0][0] = - tmp2 * fjac[k-1][0][0]
          - tmp1 * njac[k-1][0][0]
          - tmp1 * dz1; 
        lhs[k][AA][1][0] = - tmp2 * fjac[k-1][1][0]
          - tmp1 * njac[k-1][1][0];
        lhs[k][AA][2][0] = - tmp2 * fjac[k-1][2][0]
          - tmp1 * njac[k-1][2][0];
        lhs[k][AA][3][0] = - tmp2 * fjac[k-1][3][0]
          - tmp1 * njac[k-1][3][0];
        lhs[k][AA][4][0] = - tmp2 * fjac[k-1][4][0]
          - tmp1 * njac[k-1][4][0];

        lhs[k][AA][0][1] = - tmp2 * fjac[k-1][0][1]
          - tmp1 * njac[k-1][0][1];
        lhs[k][AA][1][1] = - tmp2 * fjac[k-1][1][1]
          - tmp1 * njac[k-1][1][1]
          - tmp1 * dz2;
        lhs[k][AA][2][1] = - tmp2 * fjac[k-1][2][1]
          - tmp1 * njac[k-1][2][1];
        lhs[k][AA][3][1] = - tmp2 * fjac[k-1][3][1]
          - tmp1 * njac[k-1][3][1];
        lhs[k][AA][4][1] = - tmp2 * fjac[k-1][4][1]
          - tmp1 * njac[k-1][4][1];

        lhs[k][AA][0][2] = - tmp2 * fjac[k-1][0][2]
          - tmp1 * njac[k-1][0][2];
        lhs[k][AA][1][2] = - tmp2 * fjac[k-1][1][2]
          - tmp1 * njac[k-1][1][2];
        lhs[k][AA][2][2] = - tmp2 * fjac[k-1][2][2]
          - tmp1 * njac[k-1][2][2]
          - tmp1 * dz3;
        lhs[k][AA][3][2] = - tmp2 * fjac[k-1][3][2]
          - tmp1 * njac[k-1][3][2];
        lhs[k][AA][4][2] = - tmp2 * fjac[k-1][4][2]
          - tmp1 * njac[k-1][4][2];

        lhs[k][AA][0][3] = - tmp2 * fjac[k-1][0][3]
          - tmp1 * njac[k-1][0][3];
        lhs[k][AA][1][3] = - tmp2 * fjac[k-1][1][3]
          - tmp1 * njac[k-1][1][3];
        lhs[k][AA][2][3] = - tmp2 * fjac[k-1][2][3]
          - tmp1 * njac[k-1][2][3];
        lhs[k][AA][3][3] = - tmp2 * fjac[k-1][3][3]
          - tmp1 * njac[k-1][3][3]
          - tmp1 * dz4;
        lhs[k][AA][4][3] = - tmp2 * fjac[k-1][4][3]
          - tmp1 * njac[k-1][4][3];

        lhs[k][AA][0][4] = - tmp2 * fjac[k-1][0][4]
          - tmp1 * njac[k-1][0][4];
        lhs[k][AA][1][4] = - tmp2 * fjac[k-1][1][4]
          - tmp1 * njac[k-1][1][4];
        lhs[k][AA][2][4] = - tmp2 * fjac[k-1][2][4]
          - tmp1 * njac[k-1][2][4];
        lhs[k][AA][3][4] = - tmp2 * fjac[k-1][3][4]
          - tmp1 * njac[k-1][3][4];
        lhs[k][AA][4][4] = - tmp2 * fjac[k-1][4][4]
          - tmp1 * njac[k-1][4][4]
          - tmp1 * dz5;

        lhs[k][BB][0][0] = 1.0
          + tmp1 * 2.0 * njac[k][0][0]
          + tmp1 * 2.0 * dz1;
        lhs[k][BB][1][0] = tmp1 * 2.0 * njac[k][1][0];
        lhs[k][BB][2][0] = tmp1 * 2.0 * njac[k][2][0];
        lhs[k][BB][3][0] = tmp1 * 2.0 * njac[k][3][0];
        lhs[k][BB][4][0] = tmp1 * 2.0 * njac[k][4][0];

        lhs[k][BB][0][1] = tmp1 * 2.0 * njac[k][0][1];
        lhs[k][BB][1][1] = 1.0
          + tmp1 * 2.0 * njac[k][1][1]
          + tmp1 * 2.0 * dz2;
        lhs[k][BB][2][1] = tmp1 * 2.0 * njac[k][2][1];
        lhs[k][BB][3][1] = tmp1 * 2.0 * njac[k][3][1];
        lhs[k][BB][4][1] = tmp1 * 2.0 * njac[k][4][1];

        lhs[k][BB][0][2] = tmp1 * 2.0 * njac[k][0][2];
        lhs[k][BB][1][2] = tmp1 * 2.0 * njac[k][1][2];
        lhs[k][BB][2][2] = 1.0
          + tmp1 * 2.0 * njac[k][2][2]
          + tmp1 * 2.0 * dz3;
        lhs[k][BB][3][2] = tmp1 * 2.0 * njac[k][3][2];
        lhs[k][BB][4][2] = tmp1 * 2.0 * njac[k][4][2];

        lhs[k][BB][0][3] = tmp1 * 2.0 * njac[k][0][3];
        lhs[k][BB][1][3] = tmp1 * 2.0 * njac[k][1][3];
        lhs[k][BB][2][3] = tmp1 * 2.0 * njac[k][2][3];
        lhs[k][BB][3][3] = 1.0
          + tmp1 * 2.0 * njac[k][3][3]
          + tmp1 * 2.0 * dz4;
        lhs[k][BB][4][3] = tmp1 * 2.0 * njac[k][4][3];

        lhs[k][BB][0][4] = tmp1 * 2.0 * njac[k][0][4];
        lhs[k][BB][1][4] = tmp1 * 2.0 * njac[k][1][4];
        lhs[k][BB][2][4] = tmp1 * 2.0 * njac[k][2][4];
        lhs[k][BB][3][4] = tmp1 * 2.0 * njac[k][3][4];
        lhs[k][BB][4][4] = 1.0
          + tmp1 * 2.0 * njac[k][4][4] 
          + tmp1 * 2.0 * dz5;

        lhs[k][CC][0][0] =  tmp2 * fjac[k+1][0][0]
          - tmp1 * njac[k+1][0][0]
          - tmp1 * dz1;
        lhs[k][CC][1][0] =  tmp2 * fjac[k+1][1][0]
          - tmp1 * njac[k+1][1][0];
        lhs[k][CC][2][0] =  tmp2 * fjac[k+1][2][0]
          - tmp1 * njac[k+1][2][0];
        lhs[k][CC][3][0] =  tmp2 * fjac[k+1][3][0]
          - tmp1 * njac[k+1][3][0];
        lhs[k][CC][4][0] =  tmp2 * fjac[k+1][4][0]
          - tmp1 * njac[k+1][4][0];

        lhs[k][CC][0][1] =  tmp2 * fjac[k+1][0][1]
          - tmp1 * njac[k+1][0][1];
        lhs[k][CC][1][1] =  tmp2 * fjac[k+1][1][1]
          - tmp1 * njac[k+1][1][1]
          - tmp1 * dz2;
        lhs[k][CC][2][1] =  tmp2 * fjac[k+1][2][1]
          - tmp1 * njac[k+1][2][1];
        lhs[k][CC][3][1] =  tmp2 * fjac[k+1][3][1]
          - tmp1 * njac[k+1][3][1];
        lhs[k][CC][4][1] =  tmp2 * fjac[k+1][4][1]
          - tmp1 * njac[k+1][4][1];

        lhs[k][CC][0][2] =  tmp2 * fjac[k+1][0][2]
          - tmp1 * njac[k+1][0][2];
        lhs[k][CC][1][2] =  tmp2 * fjac[k+1][1][2]
          - tmp1 * njac[k+1][1][2];
        lhs[k][CC][2][2] =  tmp2 * fjac[k+1][2][2]
          - tmp1 * njac[k+1][2][2]
          - tmp1 * dz3;
        lhs[k][CC][3][2] =  tmp2 * fjac[k+1][3][2]
          - tmp1 * njac[k+1][3][2];
        lhs[k][CC][4][2] =  tmp2 * fjac[k+1][4][2]
          - tmp1 * njac[k+1][4][2];

        lhs[k][CC][0][3] =  tmp2 * fjac[k+1][0][3]
          - tmp1 * njac[k+1][0][3];
        lhs[k][CC][1][3] =  tmp2 * fjac[k+1][1][3]
          - tmp1 * njac[k+1][1][3];
        lhs[k][CC][2][3] =  tmp2 * fjac[k+1][2][3]
          - tmp1 * njac[k+1][2][3];
        lhs[k][CC][3][3] =  tmp2 * fjac[k+1][3][3]
          - tmp1 * njac[k+1][3][3]
          - tmp1 * dz4;
        lhs[k][CC][4][3] =  tmp2 * fjac[k+1][4][3]
          - tmp1 * njac[k+1][4][3];

        lhs[k][CC][0][4] =  tmp2 * fjac[k+1][0][4]
          - tmp1 * njac[k+1][0][4];
        lhs[k][CC][1][4] =  tmp2 * fjac[k+1][1][4]
          - tmp1 * njac[k+1][1][4];
        lhs[k][CC][2][4] =  tmp2 * fjac[k+1][2][4]
          - tmp1 * njac[k+1][2][4];
        lhs[k][CC][3][4] =  tmp2 * fjac[k+1][3][4]
          - tmp1 * njac[k+1][3][4];
        lhs[k][CC][4][4] =  tmp2 * fjac[k+1][4][4]
          - tmp1 * njac[k+1][4][4]
          - tmp1 * dz5;
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
      // c'(KMAX) and rhs'(KMAX) will be sent to next cell.
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // outer most do loops - sweeping in i direction
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // multiply c[0][j][i] by b_inverse and copy back to c
      // multiply rhs(0) by b_inverse(0) and copy to rhs
      //---------------------------------------------------------------------
      
	  //binvcrhs( lhs[0][BB], lhs[0][CC], rhs[0][j][i] );
	  binvcrhs( 0, BB, 0, CC, 0, j, i );

      //---------------------------------------------------------------------
      // begin inner most do loop
      // do all the elements of the cell unless last 
      //---------------------------------------------------------------------
      for (k = 1; k <= ksize-1; k++) {
        //-------------------------------------------------------------------
        // subtract A*lhs_vector(k-1) from lhs_vector(k)
        // 
        // rhs(k) = rhs(k) - A*rhs(k-1)
        //-------------------------------------------------------------------
        
		//matvec_sub(lhs[k][AA], rhs[k-1][j][i], rhs[k][j][i]);
		matvec_sub(k, AA, k-1, j, i, k, j, i);

        //-------------------------------------------------------------------
        // B(k) = B(k) - C(k-1)*A(k)
        // matmul_sub(AA,i,j,k,c,CC,i,j,k-1,c,BB,i,j,k)
        //-------------------------------------------------------------------
        
		//matmul_sub(lhs[k][AA], lhs[k-1][CC], lhs[k][BB]);
		matmul_sub(k, AA, k-1, CC, k, BB);

        //-------------------------------------------------------------------
        // multiply c[k][j][i] by b_inverse and copy back to c
        // multiply rhs[0][j][i] by b_inverse[0][j][i] and copy to rhs
        //-------------------------------------------------------------------
        
		//binvcrhs( lhs[k][BB], lhs[k][CC], rhs[k][j][i] );
		binvcrhs( k, BB, k, CC, k, j, i );
      }

      //---------------------------------------------------------------------
      // Now finish up special cases for last cell
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // rhs(ksize) = rhs(ksize) - A*rhs(ksize-1)
      //---------------------------------------------------------------------
      
	  //matvec_sub(lhs[ksize][AA], rhs[ksize-1][j][i], rhs[ksize][j][i]);
	  matvec_sub(ksize, AA, ksize-1, j, i, ksize, j, i);

      //---------------------------------------------------------------------
      // B(ksize) = B(ksize) - C(ksize-1)*A(ksize)
      // matmul_sub(AA,i,j,ksize,c,
      // $              CC,i,j,ksize-1,c,BB,i,j,ksize)
      //---------------------------------------------------------------------
      
	  //matmul_sub(lhs[ksize][AA], lhs[ksize-1][CC], lhs[ksize][BB]);
	  matmul_sub(ksize, AA, ksize-1, CC, ksize, BB);
	  

      //---------------------------------------------------------------------
      // multiply rhs(ksize) by b_inverse(ksize) and copy to rhs
      //---------------------------------------------------------------------
      
	  //binvrhs( lhs[ksize][BB], rhs[ksize][j][i] );
	  binvrhs( ksize, BB, ksize, j, i );

      //---------------------------------------------------------------------
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // back solve: if last cell, then generate U(ksize)=rhs(ksize)
      // else assume U(ksize) is loaded in un pack backsub_info
      // so just use it
      // after u(kstart) will be sent to next cell
      //---------------------------------------------------------------------

      for (k = ksize-1; k >= 0; k--) {
        for (m = 0; m < BLOCK_SIZE; m++) {
          for (n = 0; n < BLOCK_SIZE; n++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] 
              - lhs[k][CC][n][m]*rhs[k+1][j][i][n];
          }
        }
      }
    }
  }
  if (timeron) timer_stop(t_zsolve);
  */
  
  for (j = 1; j <= grid_points[1]-2; j++) {
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (k = 1; k <= ksize-1; k++) {
		
		if(k == 1){
			a = -1;
			b = 2;
		}
		if((k >= 2) && (k <= ksize - 2)){
			a = -2;
			b = 2;
		}
		if(k == ksize - 1){
			a = -2;
			b = 1;
		} 
		  
		for(z = a; z <= b; z++){
			
			tmp1 = 1.0 / (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(0)));
			tmp2 = tmp1 * tmp1;
			tmp3 = tmp1 * tmp2;

			fjac[k + z][0][0] = 0.0;
			fjac[k + z][1][0] = 0.0;
			fjac[k + z][2][0] = 0.0;
			fjac[k + z][3][0] = 1.0;
			fjac[k + z][4][0] = 0.0;

			fjac[k + z][0][1] = - ( (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(1)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3))) ) * tmp2;
			fjac[k + z][1][1] = (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3))) * tmp1;
			fjac[k + z][2][1] = 0.0;
			fjac[k + z][3][1] = (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(1))) * tmp1;
			fjac[k + z][4][1] = 0.0;

			fjac[k + z][0][2] = - ( (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(2)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3))) ) * tmp2;
			fjac[k + z][1][2] = 0.0;
			fjac[k + z][2][2] = (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3))) * tmp1;
			fjac[k + z][3][2] = (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(2))) * tmp1;
			fjac[k + z][4][2] = 0.0;

			fjac[k + z][0][3] = - ((*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3))) * tmp2 ) 
			  + c2 * (*(double *)dvmh_get_element_addr_C(qs, 3, DVM0C(k + z), DVM0C(j), DVM0C(i)));
			fjac[k + z][1][3] = - c2 *  (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(1))) * tmp1;
			fjac[k + z][2][3] = - c2 *  (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(2))) * tmp1;
			fjac[k + z][3][3] = ( 2.0 - c2 ) *  (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3))) * tmp1;
			fjac[k + z][4][3] = c2;

			fjac[k + z][0][4] = ( c2 * 2.0 * (*(double *)dvmh_get_element_addr_C(square, 3, DVM0C(k + z), DVM0C(j), DVM0C(i))) - c1 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(4))) )
			  * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3))) * tmp2;
			fjac[k + z][1][4] = - c2 * ( (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(1)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3))) ) * tmp2;
			fjac[k + z][2][4] = - c2 * ( (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(2)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3))) ) * tmp2;
			fjac[k + z][3][4] = c1 * ( (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(4))) * tmp1 )
			  - c2 * ( (*(double *)dvmh_get_element_addr_C(qs, 3, DVM0C(k + z), DVM0C(j), DVM0C(i))) + (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3))) * tmp2 );
			fjac[k + z][4][4] = c1 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3))) * tmp1;

			njac[k + z][0][0] = 0.0;
			njac[k + z][1][0] = 0.0;
			njac[k + z][2][0] = 0.0;
			njac[k + z][3][0] = 0.0;
			njac[k + z][4][0] = 0.0;

			njac[k + z][0][1] = - c3c4 * tmp2 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(1)));
			njac[k + z][1][1] =   c3c4 * tmp1;
			njac[k + z][2][1] =   0.0;
			njac[k + z][3][1] =   0.0;
			njac[k + z][4][1] =   0.0;

			njac[k + z][0][2] = - c3c4 * tmp2 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(2)));
			njac[k + z][1][2] =   0.0;
			njac[k + z][2][2] =   c3c4 * tmp1;
			njac[k + z][3][2] =   0.0;
			njac[k + z][4][2] =   0.0;

			njac[k + z][0][3] = - con43 * c3c4 * tmp2 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3)));
			njac[k + z][1][3] =   0.0;
			njac[k + z][2][3] =   0.0;
			njac[k + z][3][3] =   con43 * c3 * c4 * tmp1;
			njac[k + z][4][3] =   0.0;

			njac[k + z][0][4] = - (  c3c4
				- c1345 ) * tmp3 * ((*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(1)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(1))))
			  - ( c3c4 - c1345 ) * tmp3 * ((*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(2)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(2))))
			  - ( con43 * c3c4
				  - c1345 ) * tmp3 * ((*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3))))
			  - c1345 * tmp2 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(4)));

			njac[k + z][1][4] = (  c3c4 - c1345 ) * tmp2 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(1)));
			njac[k + z][2][4] = (  c3c4 - c1345 ) * tmp2 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(2)));
			njac[k + z][3][4] = ( con43 * c3c4
				- c1345 ) * tmp2 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k + z), DVM0C(j), DVM0C(i), DVM0C(3)));
			njac[k + z][4][4] = ( c1345 )* tmp1;
		}
      //}

      //---------------------------------------------------------------------
      // now jacobians set, so form left hand side in z direction
      //---------------------------------------------------------------------
        if(k == 1){
		  lhsinit(lhs, ksize);
		}
      //for (k = 1; k <= ksize-1; k++) {
        tmp1 = dt * tz1;
        tmp2 = dt * tz2;

        lhs[k][AA][0][0] = - tmp2 * fjac[k-1][0][0]
          - tmp1 * njac[k-1][0][0]
          - tmp1 * dz1; 
        lhs[k][AA][1][0] = - tmp2 * fjac[k-1][1][0]
          - tmp1 * njac[k-1][1][0];
        lhs[k][AA][2][0] = - tmp2 * fjac[k-1][2][0]
          - tmp1 * njac[k-1][2][0];
        lhs[k][AA][3][0] = - tmp2 * fjac[k-1][3][0]
          - tmp1 * njac[k-1][3][0];
        lhs[k][AA][4][0] = - tmp2 * fjac[k-1][4][0]
          - tmp1 * njac[k-1][4][0];

        lhs[k][AA][0][1] = - tmp2 * fjac[k-1][0][1]
          - tmp1 * njac[k-1][0][1];
        lhs[k][AA][1][1] = - tmp2 * fjac[k-1][1][1]
          - tmp1 * njac[k-1][1][1]
          - tmp1 * dz2;
        lhs[k][AA][2][1] = - tmp2 * fjac[k-1][2][1]
          - tmp1 * njac[k-1][2][1];
        lhs[k][AA][3][1] = - tmp2 * fjac[k-1][3][1]
          - tmp1 * njac[k-1][3][1];
        lhs[k][AA][4][1] = - tmp2 * fjac[k-1][4][1]
          - tmp1 * njac[k-1][4][1];

        lhs[k][AA][0][2] = - tmp2 * fjac[k-1][0][2]
          - tmp1 * njac[k-1][0][2];
        lhs[k][AA][1][2] = - tmp2 * fjac[k-1][1][2]
          - tmp1 * njac[k-1][1][2];
        lhs[k][AA][2][2] = - tmp2 * fjac[k-1][2][2]
          - tmp1 * njac[k-1][2][2]
          - tmp1 * dz3;
        lhs[k][AA][3][2] = - tmp2 * fjac[k-1][3][2]
          - tmp1 * njac[k-1][3][2];
        lhs[k][AA][4][2] = - tmp2 * fjac[k-1][4][2]
          - tmp1 * njac[k-1][4][2];

        lhs[k][AA][0][3] = - tmp2 * fjac[k-1][0][3]
          - tmp1 * njac[k-1][0][3];
        lhs[k][AA][1][3] = - tmp2 * fjac[k-1][1][3]
          - tmp1 * njac[k-1][1][3];
        lhs[k][AA][2][3] = - tmp2 * fjac[k-1][2][3]
          - tmp1 * njac[k-1][2][3];
        lhs[k][AA][3][3] = - tmp2 * fjac[k-1][3][3]
          - tmp1 * njac[k-1][3][3]
          - tmp1 * dz4;
        lhs[k][AA][4][3] = - tmp2 * fjac[k-1][4][3]
          - tmp1 * njac[k-1][4][3];

        lhs[k][AA][0][4] = - tmp2 * fjac[k-1][0][4]
          - tmp1 * njac[k-1][0][4];
        lhs[k][AA][1][4] = - tmp2 * fjac[k-1][1][4]
          - tmp1 * njac[k-1][1][4];
        lhs[k][AA][2][4] = - tmp2 * fjac[k-1][2][4]
          - tmp1 * njac[k-1][2][4];
        lhs[k][AA][3][4] = - tmp2 * fjac[k-1][3][4]
          - tmp1 * njac[k-1][3][4];
        lhs[k][AA][4][4] = - tmp2 * fjac[k-1][4][4]
          - tmp1 * njac[k-1][4][4]
          - tmp1 * dz5;

        lhs[k][BB][0][0] = 1.0
          + tmp1 * 2.0 * njac[k][0][0]
          + tmp1 * 2.0 * dz1;
        lhs[k][BB][1][0] = tmp1 * 2.0 * njac[k][1][0];
        lhs[k][BB][2][0] = tmp1 * 2.0 * njac[k][2][0];
        lhs[k][BB][3][0] = tmp1 * 2.0 * njac[k][3][0];
        lhs[k][BB][4][0] = tmp1 * 2.0 * njac[k][4][0];

        lhs[k][BB][0][1] = tmp1 * 2.0 * njac[k][0][1];
        lhs[k][BB][1][1] = 1.0
          + tmp1 * 2.0 * njac[k][1][1]
          + tmp1 * 2.0 * dz2;
        lhs[k][BB][2][1] = tmp1 * 2.0 * njac[k][2][1];
        lhs[k][BB][3][1] = tmp1 * 2.0 * njac[k][3][1];
        lhs[k][BB][4][1] = tmp1 * 2.0 * njac[k][4][1];

        lhs[k][BB][0][2] = tmp1 * 2.0 * njac[k][0][2];
        lhs[k][BB][1][2] = tmp1 * 2.0 * njac[k][1][2];
        lhs[k][BB][2][2] = 1.0
          + tmp1 * 2.0 * njac[k][2][2]
          + tmp1 * 2.0 * dz3;
        lhs[k][BB][3][2] = tmp1 * 2.0 * njac[k][3][2];
        lhs[k][BB][4][2] = tmp1 * 2.0 * njac[k][4][2];

        lhs[k][BB][0][3] = tmp1 * 2.0 * njac[k][0][3];
        lhs[k][BB][1][3] = tmp1 * 2.0 * njac[k][1][3];
        lhs[k][BB][2][3] = tmp1 * 2.0 * njac[k][2][3];
        lhs[k][BB][3][3] = 1.0
          + tmp1 * 2.0 * njac[k][3][3]
          + tmp1 * 2.0 * dz4;
        lhs[k][BB][4][3] = tmp1 * 2.0 * njac[k][4][3];

        lhs[k][BB][0][4] = tmp1 * 2.0 * njac[k][0][4];
        lhs[k][BB][1][4] = tmp1 * 2.0 * njac[k][1][4];
        lhs[k][BB][2][4] = tmp1 * 2.0 * njac[k][2][4];
        lhs[k][BB][3][4] = tmp1 * 2.0 * njac[k][3][4];
        lhs[k][BB][4][4] = 1.0
          + tmp1 * 2.0 * njac[k][4][4] 
          + tmp1 * 2.0 * dz5;

        lhs[k][CC][0][0] =  tmp2 * fjac[k+1][0][0]
          - tmp1 * njac[k+1][0][0]
          - tmp1 * dz1;
        lhs[k][CC][1][0] =  tmp2 * fjac[k+1][1][0]
          - tmp1 * njac[k+1][1][0];
        lhs[k][CC][2][0] =  tmp2 * fjac[k+1][2][0]
          - tmp1 * njac[k+1][2][0];
        lhs[k][CC][3][0] =  tmp2 * fjac[k+1][3][0]
          - tmp1 * njac[k+1][3][0];
        lhs[k][CC][4][0] =  tmp2 * fjac[k+1][4][0]
          - tmp1 * njac[k+1][4][0];

        lhs[k][CC][0][1] =  tmp2 * fjac[k+1][0][1]
          - tmp1 * njac[k+1][0][1];
        lhs[k][CC][1][1] =  tmp2 * fjac[k+1][1][1]
          - tmp1 * njac[k+1][1][1]
          - tmp1 * dz2;
        lhs[k][CC][2][1] =  tmp2 * fjac[k+1][2][1]
          - tmp1 * njac[k+1][2][1];
        lhs[k][CC][3][1] =  tmp2 * fjac[k+1][3][1]
          - tmp1 * njac[k+1][3][1];
        lhs[k][CC][4][1] =  tmp2 * fjac[k+1][4][1]
          - tmp1 * njac[k+1][4][1];

        lhs[k][CC][0][2] =  tmp2 * fjac[k+1][0][2]
          - tmp1 * njac[k+1][0][2];
        lhs[k][CC][1][2] =  tmp2 * fjac[k+1][1][2]
          - tmp1 * njac[k+1][1][2];
        lhs[k][CC][2][2] =  tmp2 * fjac[k+1][2][2]
          - tmp1 * njac[k+1][2][2]
          - tmp1 * dz3;
        lhs[k][CC][3][2] =  tmp2 * fjac[k+1][3][2]
          - tmp1 * njac[k+1][3][2];
        lhs[k][CC][4][2] =  tmp2 * fjac[k+1][4][2]
          - tmp1 * njac[k+1][4][2];

        lhs[k][CC][0][3] =  tmp2 * fjac[k+1][0][3]
          - tmp1 * njac[k+1][0][3];
        lhs[k][CC][1][3] =  tmp2 * fjac[k+1][1][3]
          - tmp1 * njac[k+1][1][3];
        lhs[k][CC][2][3] =  tmp2 * fjac[k+1][2][3]
          - tmp1 * njac[k+1][2][3];
        lhs[k][CC][3][3] =  tmp2 * fjac[k+1][3][3]
          - tmp1 * njac[k+1][3][3]
          - tmp1 * dz4;
        lhs[k][CC][4][3] =  tmp2 * fjac[k+1][4][3]
          - tmp1 * njac[k+1][4][3];

        lhs[k][CC][0][4] =  tmp2 * fjac[k+1][0][4]
          - tmp1 * njac[k+1][0][4];
        lhs[k][CC][1][4] =  tmp2 * fjac[k+1][1][4]
          - tmp1 * njac[k+1][1][4];
        lhs[k][CC][2][4] =  tmp2 * fjac[k+1][2][4]
          - tmp1 * njac[k+1][2][4];
        lhs[k][CC][3][4] =  tmp2 * fjac[k+1][3][4]
          - tmp1 * njac[k+1][3][4];
        lhs[k][CC][4][4] =  tmp2 * fjac[k+1][4][4]
          - tmp1 * njac[k+1][4][4]
          - tmp1 * dz5;
      //}

        
        if(k == 1){	
	  	  binvcrhs( 0, BB, 0, CC, 0, j, i );
		}

              
		//matvec_sub(lhs[k][AA], rhs[k-1][j][i], rhs[k][j][i]);
		matvec_sub(k, AA, k-1, j, i, k, j, i);

        
		//matmul_sub(lhs[k][AA], lhs[k-1][CC], lhs[k][BB]);/
		matmul_sub(k, AA, k-1, CC, k, BB);

        //-------------------------------------------------------------------
        // multiply c[k][j][i] by b_inverse and copy back to c
        // multiply rhs[0][j][i] by b_inverse[0][j][i] and copy to rhs
        //-------------------------------------------------------------------
        
		//binvcrhs( lhs[k][BB], lhs[k][CC], rhs[k][j][i] );
		binvcrhs( k, BB, k, CC, k, j, i );
      //}

	    if(k == ksize-1){
     	  matvec_sub(ksize, AA, ksize-1, j, i, ksize, j, i);

      	  matmul_sub(ksize, AA, ksize-1, CC, ksize, BB);
	  
	      binvrhs( ksize, BB, ksize, j, i );
		}

		for(m = 0; m < BLOCK_SIZE; m++){
			for(n = 0; n < BLOCK_SIZE; n++){				
				lhs_buf[k][j][i][n][m] = lhs[k][CC][n][m];
			}
		}
      }
    }
  }
  
  for (j = 1; j <= grid_points[1]-2; j++) {
    for (i = 1; i <= grid_points[0]-2; i++) {  
      for (k = ksize-1; k >= 0; k--) {
        for (m = 0; m < BLOCK_SIZE; m++) {
          for (n = 0; n < BLOCK_SIZE; n++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] 
              - lhs_buf[k][j][i][n][m]*rhs[k+1][j][i][n];
          }
        }
      }
    }
  }
  
  
  
}

