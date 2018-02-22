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
  int i, j, k, m, n, isize;
  int z, a = -1, b = 2;
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  if (timeron) timer_start(t_xsolve);

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  // This function computes the left hand side in the xi-direction
  //---------------------------------------------------------------------

  isize = grid_points[0]-1;

  //---------------------------------------------------------------------
  // determine a (labeled f) and n jacobians
  //---------------------------------------------------------------------
  /*for (k = 1; k <= grid_points[2]-2; k++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (i = 0; i <= isize; i++) {
        tmp1 = rho_i[k][j][i];
        tmp2 = tmp1 * tmp1;
        tmp3 = tmp1 * tmp2;
        //-------------------------------------------------------------------
        // 
        //-------------------------------------------------------------------
        fjac[i][0][0] = 0.0;
        fjac[i][1][0] = 1.0;
        fjac[i][2][0] = 0.0;
        fjac[i][3][0] = 0.0;
        fjac[i][4][0] = 0.0;

        fjac[i][0][1] = -(u[k][j][i][1] * tmp2 * u[k][j][i][1])
          + c2 * qs[k][j][i];
        fjac[i][1][1] = ( 2.0 - c2 ) * ( u[k][j][i][1] / u[k][j][i][0] );
        fjac[i][2][1] = - c2 * ( u[k][j][i][2] * tmp1 );
        fjac[i][3][1] = - c2 * ( u[k][j][i][3] * tmp1 );
        fjac[i][4][1] = c2;

        fjac[i][0][2] = - ( u[k][j][i][1]*u[k][j][i][2] ) * tmp2;
        fjac[i][1][2] = u[k][j][i][2] * tmp1;
        fjac[i][2][2] = u[k][j][i][1] * tmp1;
        fjac[i][3][2] = 0.0;
        fjac[i][4][2] = 0.0;

        fjac[i][0][3] = - ( u[k][j][i][1]*u[k][j][i][3] ) * tmp2;
        fjac[i][1][3] = u[k][j][i][3] * tmp1;
        fjac[i][2][3] = 0.0;
        fjac[i][3][3] = u[k][j][i][1] * tmp1;
        fjac[i][4][3] = 0.0;

        fjac[i][0][4] = ( c2 * 2.0 * square[k][j][i] - c1 * u[k][j][i][4] )
          * ( u[k][j][i][1] * tmp2 );
        fjac[i][1][4] = c1 *  u[k][j][i][4] * tmp1 
          - c2 * ( u[k][j][i][1]*u[k][j][i][1] * tmp2 + qs[k][j][i] );
        fjac[i][2][4] = - c2 * ( u[k][j][i][2]*u[k][j][i][1] ) * tmp2;
        fjac[i][3][4] = - c2 * ( u[k][j][i][3]*u[k][j][i][1] ) * tmp2;
        fjac[i][4][4] = c1 * ( u[k][j][i][1] * tmp1 );

        njac[i][0][0] = 0.0;
        njac[i][1][0] = 0.0;
        njac[i][2][0] = 0.0;
        njac[i][3][0] = 0.0;
        njac[i][4][0] = 0.0;

        njac[i][0][1] = - con43 * c3c4 * tmp2 * u[k][j][i][1];
        njac[i][1][1] =   con43 * c3c4 * tmp1;
        njac[i][2][1] =   0.0;
        njac[i][3][1] =   0.0;
        njac[i][4][1] =   0.0;

        njac[i][0][2] = - c3c4 * tmp2 * u[k][j][i][2];
        njac[i][1][2] =   0.0;
        njac[i][2][2] =   c3c4 * tmp1;
        njac[i][3][2] =   0.0;
        njac[i][4][2] =   0.0;

        njac[i][0][3] = - c3c4 * tmp2 * u[k][j][i][3];
        njac[i][1][3] =   0.0;
        njac[i][2][3] =   0.0;
        njac[i][3][3] =   c3c4 * tmp1;
        njac[i][4][3] =   0.0;

        njac[i][0][4] = - ( con43 * c3c4
            - c1345 ) * tmp3 * (u[k][j][i][1]*u[k][j][i][1])
          - ( c3c4 - c1345 ) * tmp3 * (u[k][j][i][2]*u[k][j][i][2])
          - ( c3c4 - c1345 ) * tmp3 * (u[k][j][i][3]*u[k][j][i][3])
          - c1345 * tmp2 * u[k][j][i][4];

        njac[i][1][4] = ( con43 * c3c4
            - c1345 ) * tmp2 * u[k][j][i][1];
        njac[i][2][4] = ( c3c4 - c1345 ) * tmp2 * u[k][j][i][2];
        njac[i][3][4] = ( c3c4 - c1345 ) * tmp2 * u[k][j][i][3];
        njac[i][4][4] = ( c1345 ) * tmp1;
      }
      //---------------------------------------------------------------------
      // now jacobians set, so form left hand side in x direction
      //---------------------------------------------------------------------
      lhsinit(lhs, isize);
      for (i = 1; i <= isize-1; i++) {
        tmp1 = dt * tx1;
        tmp2 = dt * tx2;

        lhs[i][AA][0][0] = - tmp2 * fjac[i-1][0][0]
          - tmp1 * njac[i-1][0][0]
          - tmp1 * dx1; 
        lhs[i][AA][1][0] = - tmp2 * fjac[i-1][1][0]
          - tmp1 * njac[i-1][1][0];
        lhs[i][AA][2][0] = - tmp2 * fjac[i-1][2][0]
          - tmp1 * njac[i-1][2][0];
        lhs[i][AA][3][0] = - tmp2 * fjac[i-1][3][0]
          - tmp1 * njac[i-1][3][0];
        lhs[i][AA][4][0] = - tmp2 * fjac[i-1][4][0]
          - tmp1 * njac[i-1][4][0];

        lhs[i][AA][0][1] = - tmp2 * fjac[i-1][0][1]
          - tmp1 * njac[i-1][0][1];
        lhs[i][AA][1][1] = - tmp2 * fjac[i-1][1][1]
          - tmp1 * njac[i-1][1][1]
          - tmp1 * dx2;
        lhs[i][AA][2][1] = - tmp2 * fjac[i-1][2][1]
          - tmp1 * njac[i-1][2][1];
        lhs[i][AA][3][1] = - tmp2 * fjac[i-1][3][1]
          - tmp1 * njac[i-1][3][1];
        lhs[i][AA][4][1] = - tmp2 * fjac[i-1][4][1]
          - tmp1 * njac[i-1][4][1];

        lhs[i][AA][0][2] = - tmp2 * fjac[i-1][0][2]
          - tmp1 * njac[i-1][0][2];
        lhs[i][AA][1][2] = - tmp2 * fjac[i-1][1][2]
          - tmp1 * njac[i-1][1][2];
        lhs[i][AA][2][2] = - tmp2 * fjac[i-1][2][2]
          - tmp1 * njac[i-1][2][2]
          - tmp1 * dx3;
        lhs[i][AA][3][2] = - tmp2 * fjac[i-1][3][2]
          - tmp1 * njac[i-1][3][2];
        lhs[i][AA][4][2] = - tmp2 * fjac[i-1][4][2]
          - tmp1 * njac[i-1][4][2];

        lhs[i][AA][0][3] = - tmp2 * fjac[i-1][0][3]
          - tmp1 * njac[i-1][0][3];
        lhs[i][AA][1][3] = - tmp2 * fjac[i-1][1][3]
          - tmp1 * njac[i-1][1][3];
        lhs[i][AA][2][3] = - tmp2 * fjac[i-1][2][3]
          - tmp1 * njac[i-1][2][3];
        lhs[i][AA][3][3] = - tmp2 * fjac[i-1][3][3]
          - tmp1 * njac[i-1][3][3]
          - tmp1 * dx4;
        lhs[i][AA][4][3] = - tmp2 * fjac[i-1][4][3]
          - tmp1 * njac[i-1][4][3];

        lhs[i][AA][0][4] = - tmp2 * fjac[i-1][0][4]
          - tmp1 * njac[i-1][0][4];
        lhs[i][AA][1][4] = - tmp2 * fjac[i-1][1][4]
          - tmp1 * njac[i-1][1][4];
        lhs[i][AA][2][4] = - tmp2 * fjac[i-1][2][4]
          - tmp1 * njac[i-1][2][4];
        lhs[i][AA][3][4] = - tmp2 * fjac[i-1][3][4]
          - tmp1 * njac[i-1][3][4];
        lhs[i][AA][4][4] = - tmp2 * fjac[i-1][4][4]
          - tmp1 * njac[i-1][4][4]
          - tmp1 * dx5;

        lhs[i][BB][0][0] = 1.0
          + tmp1 * 2.0 * njac[i][0][0]
          + tmp1 * 2.0 * dx1;
        lhs[i][BB][1][0] = tmp1 * 2.0 * njac[i][1][0];
        lhs[i][BB][2][0] = tmp1 * 2.0 * njac[i][2][0];
        lhs[i][BB][3][0] = tmp1 * 2.0 * njac[i][3][0];
        lhs[i][BB][4][0] = tmp1 * 2.0 * njac[i][4][0];

        lhs[i][BB][0][1] = tmp1 * 2.0 * njac[i][0][1];
        lhs[i][BB][1][1] = 1.0
          + tmp1 * 2.0 * njac[i][1][1]
          + tmp1 * 2.0 * dx2;
        lhs[i][BB][2][1] = tmp1 * 2.0 * njac[i][2][1];
        lhs[i][BB][3][1] = tmp1 * 2.0 * njac[i][3][1];
        lhs[i][BB][4][1] = tmp1 * 2.0 * njac[i][4][1];

        lhs[i][BB][0][2] = tmp1 * 2.0 * njac[i][0][2];
        lhs[i][BB][1][2] = tmp1 * 2.0 * njac[i][1][2];
        lhs[i][BB][2][2] = 1.0
          + tmp1 * 2.0 * njac[i][2][2]
          + tmp1 * 2.0 * dx3;
        lhs[i][BB][3][2] = tmp1 * 2.0 * njac[i][3][2];
        lhs[i][BB][4][2] = tmp1 * 2.0 * njac[i][4][2];

        lhs[i][BB][0][3] = tmp1 * 2.0 * njac[i][0][3];
        lhs[i][BB][1][3] = tmp1 * 2.0 * njac[i][1][3];
        lhs[i][BB][2][3] = tmp1 * 2.0 * njac[i][2][3];
        lhs[i][BB][3][3] = 1.0
          + tmp1 * 2.0 * njac[i][3][3]
          + tmp1 * 2.0 * dx4;
        lhs[i][BB][4][3] = tmp1 * 2.0 * njac[i][4][3];

        lhs[i][BB][0][4] = tmp1 * 2.0 * njac[i][0][4];
        lhs[i][BB][1][4] = tmp1 * 2.0 * njac[i][1][4];
        lhs[i][BB][2][4] = tmp1 * 2.0 * njac[i][2][4];
        lhs[i][BB][3][4] = tmp1 * 2.0 * njac[i][3][4];
        lhs[i][BB][4][4] = 1.0
          + tmp1 * 2.0 * njac[i][4][4]
          + tmp1 * 2.0 * dx5;

        lhs[i][CC][0][0] =  tmp2 * fjac[i+1][0][0]
          - tmp1 * njac[i+1][0][0]
          - tmp1 * dx1;
        lhs[i][CC][1][0] =  tmp2 * fjac[i+1][1][0]
          - tmp1 * njac[i+1][1][0];
        lhs[i][CC][2][0] =  tmp2 * fjac[i+1][2][0]
          - tmp1 * njac[i+1][2][0];
        lhs[i][CC][3][0] =  tmp2 * fjac[i+1][3][0]
          - tmp1 * njac[i+1][3][0];
        lhs[i][CC][4][0] =  tmp2 * fjac[i+1][4][0]
          - tmp1 * njac[i+1][4][0];

        lhs[i][CC][0][1] =  tmp2 * fjac[i+1][0][1]
          - tmp1 * njac[i+1][0][1];
        lhs[i][CC][1][1] =  tmp2 * fjac[i+1][1][1]
          - tmp1 * njac[i+1][1][1]
          - tmp1 * dx2;
        lhs[i][CC][2][1] =  tmp2 * fjac[i+1][2][1]
          - tmp1 * njac[i+1][2][1];
        lhs[i][CC][3][1] =  tmp2 * fjac[i+1][3][1]
          - tmp1 * njac[i+1][3][1];
        lhs[i][CC][4][1] =  tmp2 * fjac[i+1][4][1]
          - tmp1 * njac[i+1][4][1];

        lhs[i][CC][0][2] =  tmp2 * fjac[i+1][0][2]
          - tmp1 * njac[i+1][0][2];
        lhs[i][CC][1][2] =  tmp2 * fjac[i+1][1][2]
          - tmp1 * njac[i+1][1][2];
        lhs[i][CC][2][2] =  tmp2 * fjac[i+1][2][2]
          - tmp1 * njac[i+1][2][2]
          - tmp1 * dx3;
        lhs[i][CC][3][2] =  tmp2 * fjac[i+1][3][2]
          - tmp1 * njac[i+1][3][2];
        lhs[i][CC][4][2] =  tmp2 * fjac[i+1][4][2]
          - tmp1 * njac[i+1][4][2];

        lhs[i][CC][0][3] =  tmp2 * fjac[i+1][0][3]
          - tmp1 * njac[i+1][0][3];
        lhs[i][CC][1][3] =  tmp2 * fjac[i+1][1][3]
          - tmp1 * njac[i+1][1][3];
        lhs[i][CC][2][3] =  tmp2 * fjac[i+1][2][3]
          - tmp1 * njac[i+1][2][3];
        lhs[i][CC][3][3] =  tmp2 * fjac[i+1][3][3]
          - tmp1 * njac[i+1][3][3]
          - tmp1 * dx4;
        lhs[i][CC][4][3] =  tmp2 * fjac[i+1][4][3]
          - tmp1 * njac[i+1][4][3];

        lhs[i][CC][0][4] =  tmp2 * fjac[i+1][0][4]
          - tmp1 * njac[i+1][0][4];
        lhs[i][CC][1][4] =  tmp2 * fjac[i+1][1][4]
          - tmp1 * njac[i+1][1][4];
        lhs[i][CC][2][4] =  tmp2 * fjac[i+1][2][4]
          - tmp1 * njac[i+1][2][4];
        lhs[i][CC][3][4] =  tmp2 * fjac[i+1][3][4]
          - tmp1 * njac[i+1][3][4];
        lhs[i][CC][4][4] =  tmp2 * fjac[i+1][4][4]
          - tmp1 * njac[i+1][4][4]
          - tmp1 * dx5;
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
      // c'(IMAX) and rhs'(IMAX) will be sent to next cell
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // outer most do loops - sweeping in i direction
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // multiply c[k][j][0] by b_inverse and copy back to c
      // multiply rhs(0) by b_inverse(0) and copy to rhs
      //---------------------------------------------------------------------
      
	  //binvcrhs( lhs[0][BB], lhs[0][CC], rhs[k][j][0] );
	  binvcrhs( 0, BB, 0, CC, k, j, 0 );

      //---------------------------------------------------------------------
      // begin inner most do loop
      // do all the elements of the cell unless last 
      //---------------------------------------------------------------------
      for (i = 1; i <= isize-1; i++) {
        //-------------------------------------------------------------------
        // rhs(i) = rhs(i) - A*rhs(i-1)
        //-------------------------------------------------------------------
		
        //matvec_sub(lhs[i][AA], rhs[k][j][i-1], rhs[k][j][i]);
		matvec_sub(i, AA, k, j, i-1, k, j, i);
		
        //-------------------------------------------------------------------
        // B(i) = B(i) - C(i-1)*A(i)
        //-------------------------------------------------------------------
        
		//matmul_sub(lhs[i][AA], lhs[i-1][CC], lhs[i][BB]);
		matmul_sub(i, AA, i-1, CC, i, BB);


        //-------------------------------------------------------------------
        // multiply c[k][j][i] by b_inverse and copy back to c
        // multiply rhs[k][j][0] by b_inverse[k][j][0] and copy to rhs
        //-------------------------------------------------------------------
        
		//binvcrhs( lhs[i][BB], lhs[i][CC], rhs[k][j][i] );
		binvcrhs( i, BB, i, CC, k, j, i );
      }

      //---------------------------------------------------------------------
      // rhs(isize) = rhs(isize) - A*rhs(isize-1)
      //---------------------------------------------------------------------
	  
      //matvec_sub(lhs[isize][AA], rhs[k][j][isize-1], rhs[k][j][isize]);
	  matvec_sub(isize, AA, k, j, isize-1, k, j, isize);
	  
      //---------------------------------------------------------------------
      // B(isize) = B(isize) - C(isize-1)*A(isize)
      //---------------------------------------------------------------------
      
	  //matmul_sub(lhs[isize][AA], lhs[isize-1][CC], lhs[isize][BB]);
	  matmul_sub(isize, AA, isize-1, CC, isize, BB);

      //---------------------------------------------------------------------
      // multiply rhs() by b_inverse() and copy to rhs
      //---------------------------------------------------------------------
      
	  //binvrhs( lhs[isize][BB], rhs[k][j][isize] );
	  binvrhs( isize, BB, k, j, isize );

      //---------------------------------------------------------------------
      // back solve: if last cell, then generate U(isize)=rhs(isize)
      // else assume U(isize) is loaded in un pack backsub_info
      // so just use it
      // after u(istart) will be sent to next cell
      //---------------------------------------------------------------------
      for (i = isize-1; i >=0; i--) {
        for (m = 0; m < BLOCK_SIZE; m++) {
          for (n = 0; n < BLOCK_SIZE; n++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] 
              - lhs[i][CC][n][m]*rhs[k][j][i+1][n];
          }
        }
      }
    }
  }
  if (timeron) timer_stop(t_xsolve);
  */

  for (k = 1; k <= grid_points[2]-2; k++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (i = 1; i <= isize-1; i++) {
		  
		if(i == 1){
			a = -1;
			b = 2;
		}
		if((i >= 2) && (i <= isize - 2)){
			a = -2;
			b = 2;
		}
		if(i == isize - 1){
			a = -2;
			b = 1;
		} 
		  
		for(z = a; z <= b; z++){  
			tmp1 = (*(double *)dvmh_get_element_addr_C(rho_i, 3, DVM0C(k), DVM0C(j), DVM0C(i + z)));
			tmp2 = tmp1 * tmp1;
			tmp3 = tmp1 * tmp2;
			//-------------------------------------------------------------------
			// 
			//-------------------------------------------------------------------
			fjac[i + z][0][0] = 0.0;
			fjac[i + z][1][0] = 1.0;
			fjac[i + z][2][0] = 0.0;
			fjac[i + z][3][0] = 0.0;
			fjac[i + z][4][0] = 0.0;

			fjac[i + z][0][1] = -((*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1))) * tmp2 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1))))
			  + c2 * (*(double *)dvmh_get_element_addr_C(qs, 3, DVM0C(k), DVM0C(j), DVM0C(i + z)));
			fjac[i + z][1][1] = ( 2.0 - c2 ) * ( (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1))) / (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(0))) );
			fjac[i + z][2][1] = - c2 * ( (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(2))) * tmp1 );
			fjac[i + z][3][1] = - c2 * ( (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(3))) * tmp1 );
			fjac[i + z][4][1] = c2;

			fjac[i + z][0][2] = - ( (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(2))) ) * tmp2;
			fjac[i + z][1][2] = (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(2))) * tmp1;
			fjac[i + z][2][2] = (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1))) * tmp1;
			fjac[i + z][3][2] = 0.0;
			fjac[i + z][4][2] = 0.0;

			fjac[i + z][0][3] = - ( (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(3))) ) * tmp2;
			fjac[i + z][1][3] = (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(3))) * tmp1;
			fjac[i + z][2][3] = 0.0;
			fjac[i + z][3][3] = (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1))) * tmp1;
			fjac[i + z][4][3] = 0.0;

			fjac[i + z][0][4] = ( c2 * 2.0 * (*(double *)dvmh_get_element_addr_C(square, 3, DVM0C(k), DVM0C(j), DVM0C(i + z))) - c1 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(4))) )
			  * ( (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1))) * tmp2 );
			fjac[i + z][1][4] = c1 *  (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(4))) * tmp1 
			  - c2 * ( (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1))) * tmp2 + (*(double *)dvmh_get_element_addr_C(qs, 3, DVM0C(k), DVM0C(j), DVM0C(i + z))) );
			fjac[i + z][2][4] = - c2 * ( (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(2)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1))) ) * tmp2;
			fjac[i + z][3][4] = - c2 * ( (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(3)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1))) ) * tmp2;
			fjac[i + z][4][4] = c1 * ( (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1))) * tmp1 );

			njac[i + z][0][0] = 0.0;
			njac[i + z][1][0] = 0.0;
			njac[i + z][2][0] = 0.0;
			njac[i + z][3][0] = 0.0;
			njac[i + z][4][0] = 0.0;

			njac[i + z][0][1] = - con43 * c3c4 * tmp2 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1)));
			njac[i + z][1][1] =   con43 * c3c4 * tmp1;
			njac[i + z][2][1] =   0.0;
			njac[i + z][3][1] =   0.0;
			njac[i + z][4][1] =   0.0;

			njac[i + z][0][2] = - c3c4 * tmp2 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(2)));
			njac[i + z][1][2] =   0.0;
			njac[i + z][2][2] =   c3c4 * tmp1;
			njac[i + z][3][2] =   0.0;
			njac[i + z][4][2] =   0.0;

			njac[i + z][0][3] = - c3c4 * tmp2 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(3)));
			njac[i + z][1][3] =   0.0;
			njac[i + z][2][3] =   0.0;
			njac[i + z][3][3] =   c3c4 * tmp1;
			njac[i + z][4][3] =   0.0;

			njac[i + z][0][4] = - ( con43 * c3c4
				- c1345 ) * tmp3 * ((*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1))))
			  - ( c3c4 - c1345 ) * tmp3 * ((*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(2)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(2))))
			  - ( c3c4 - c1345 ) * tmp3 * ((*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(3)))*(*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(3))))
			  - c1345 * tmp2 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(4)));

			njac[i + z][1][4] = ( con43 * c3c4
				- c1345 ) * tmp2 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(1)));
			njac[i + z][2][4] = ( c3c4 - c1345 ) * tmp2 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(2)));
			njac[i + z][3][4] = ( c3c4 - c1345 ) * tmp2 * (*(double *)dvmh_get_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i + z), DVM0C(3)));
			njac[i + z][4][4] = ( c1345 ) * tmp1;
	    }
      //}
      //---------------------------------------------------------------------
      // now jacobians set, so form left hand side in x direction
      //---------------------------------------------------------------------
          
		if(i == 1){
		  lhsinit(lhs, isize);
		} 
      //for (i = 1; i <= isize-1; i++) {
        tmp1 = dt * tx1;
        tmp2 = dt * tx2;

        lhs[i][AA][0][0] = - tmp2 * fjac[i-1][0][0]
          - tmp1 * njac[i-1][0][0]
          - tmp1 * dx1; 
        lhs[i][AA][1][0] = - tmp2 * fjac[i-1][1][0]
          - tmp1 * njac[i-1][1][0];
        lhs[i][AA][2][0] = - tmp2 * fjac[i-1][2][0]
          - tmp1 * njac[i-1][2][0];
        lhs[i][AA][3][0] = - tmp2 * fjac[i-1][3][0]
          - tmp1 * njac[i-1][3][0];
        lhs[i][AA][4][0] = - tmp2 * fjac[i-1][4][0]
          - tmp1 * njac[i-1][4][0];

        lhs[i][AA][0][1] = - tmp2 * fjac[i-1][0][1]
          - tmp1 * njac[i-1][0][1];
        lhs[i][AA][1][1] = - tmp2 * fjac[i-1][1][1]
          - tmp1 * njac[i-1][1][1]
          - tmp1 * dx2;
        lhs[i][AA][2][1] = - tmp2 * fjac[i-1][2][1]
          - tmp1 * njac[i-1][2][1];
        lhs[i][AA][3][1] = - tmp2 * fjac[i-1][3][1]
          - tmp1 * njac[i-1][3][1];
        lhs[i][AA][4][1] = - tmp2 * fjac[i-1][4][1]
          - tmp1 * njac[i-1][4][1];

        lhs[i][AA][0][2] = - tmp2 * fjac[i-1][0][2]
          - tmp1 * njac[i-1][0][2];
        lhs[i][AA][1][2] = - tmp2 * fjac[i-1][1][2]
          - tmp1 * njac[i-1][1][2];
        lhs[i][AA][2][2] = - tmp2 * fjac[i-1][2][2]
          - tmp1 * njac[i-1][2][2]
          - tmp1 * dx3;
        lhs[i][AA][3][2] = - tmp2 * fjac[i-1][3][2]
          - tmp1 * njac[i-1][3][2];
        lhs[i][AA][4][2] = - tmp2 * fjac[i-1][4][2]
          - tmp1 * njac[i-1][4][2];

        lhs[i][AA][0][3] = - tmp2 * fjac[i-1][0][3]
          - tmp1 * njac[i-1][0][3];
        lhs[i][AA][1][3] = - tmp2 * fjac[i-1][1][3]
          - tmp1 * njac[i-1][1][3];
        lhs[i][AA][2][3] = - tmp2 * fjac[i-1][2][3]
          - tmp1 * njac[i-1][2][3];
        lhs[i][AA][3][3] = - tmp2 * fjac[i-1][3][3]
          - tmp1 * njac[i-1][3][3]
          - tmp1 * dx4;
        lhs[i][AA][4][3] = - tmp2 * fjac[i-1][4][3]
          - tmp1 * njac[i-1][4][3];

        lhs[i][AA][0][4] = - tmp2 * fjac[i-1][0][4]
          - tmp1 * njac[i-1][0][4];
        lhs[i][AA][1][4] = - tmp2 * fjac[i-1][1][4]
          - tmp1 * njac[i-1][1][4];
        lhs[i][AA][2][4] = - tmp2 * fjac[i-1][2][4]
          - tmp1 * njac[i-1][2][4];
        lhs[i][AA][3][4] = - tmp2 * fjac[i-1][3][4]
          - tmp1 * njac[i-1][3][4];
        lhs[i][AA][4][4] = - tmp2 * fjac[i-1][4][4]
          - tmp1 * njac[i-1][4][4]
          - tmp1 * dx5;

        lhs[i][BB][0][0] = 1.0
          + tmp1 * 2.0 * njac[i][0][0]
          + tmp1 * 2.0 * dx1;
        lhs[i][BB][1][0] = tmp1 * 2.0 * njac[i][1][0];
        lhs[i][BB][2][0] = tmp1 * 2.0 * njac[i][2][0];
        lhs[i][BB][3][0] = tmp1 * 2.0 * njac[i][3][0];
        lhs[i][BB][4][0] = tmp1 * 2.0 * njac[i][4][0];

        lhs[i][BB][0][1] = tmp1 * 2.0 * njac[i][0][1];
        lhs[i][BB][1][1] = 1.0
          + tmp1 * 2.0 * njac[i][1][1]
          + tmp1 * 2.0 * dx2;
        lhs[i][BB][2][1] = tmp1 * 2.0 * njac[i][2][1];
        lhs[i][BB][3][1] = tmp1 * 2.0 * njac[i][3][1];
        lhs[i][BB][4][1] = tmp1 * 2.0 * njac[i][4][1];

        lhs[i][BB][0][2] = tmp1 * 2.0 * njac[i][0][2];
        lhs[i][BB][1][2] = tmp1 * 2.0 * njac[i][1][2];
        lhs[i][BB][2][2] = 1.0
          + tmp1 * 2.0 * njac[i][2][2]
          + tmp1 * 2.0 * dx3;
        lhs[i][BB][3][2] = tmp1 * 2.0 * njac[i][3][2];
        lhs[i][BB][4][2] = tmp1 * 2.0 * njac[i][4][2];

        lhs[i][BB][0][3] = tmp1 * 2.0 * njac[i][0][3];
        lhs[i][BB][1][3] = tmp1 * 2.0 * njac[i][1][3];
        lhs[i][BB][2][3] = tmp1 * 2.0 * njac[i][2][3];
        lhs[i][BB][3][3] = 1.0
          + tmp1 * 2.0 * njac[i][3][3]
          + tmp1 * 2.0 * dx4;
        lhs[i][BB][4][3] = tmp1 * 2.0 * njac[i][4][3];

        lhs[i][BB][0][4] = tmp1 * 2.0 * njac[i][0][4];
        lhs[i][BB][1][4] = tmp1 * 2.0 * njac[i][1][4];
        lhs[i][BB][2][4] = tmp1 * 2.0 * njac[i][2][4];
        lhs[i][BB][3][4] = tmp1 * 2.0 * njac[i][3][4];
        lhs[i][BB][4][4] = 1.0
          + tmp1 * 2.0 * njac[i][4][4]
          + tmp1 * 2.0 * dx5;

        lhs[i][CC][0][0] =  tmp2 * fjac[i+1][0][0]
          - tmp1 * njac[i+1][0][0]
          - tmp1 * dx1;
        lhs[i][CC][1][0] =  tmp2 * fjac[i+1][1][0]
          - tmp1 * njac[i+1][1][0];
        lhs[i][CC][2][0] =  tmp2 * fjac[i+1][2][0]
          - tmp1 * njac[i+1][2][0];
        lhs[i][CC][3][0] =  tmp2 * fjac[i+1][3][0]
          - tmp1 * njac[i+1][3][0];
        lhs[i][CC][4][0] =  tmp2 * fjac[i+1][4][0]
          - tmp1 * njac[i+1][4][0];

        lhs[i][CC][0][1] =  tmp2 * fjac[i+1][0][1]
          - tmp1 * njac[i+1][0][1];
        lhs[i][CC][1][1] =  tmp2 * fjac[i+1][1][1]
          - tmp1 * njac[i+1][1][1]
          - tmp1 * dx2;
        lhs[i][CC][2][1] =  tmp2 * fjac[i+1][2][1]
          - tmp1 * njac[i+1][2][1];
        lhs[i][CC][3][1] =  tmp2 * fjac[i+1][3][1]
          - tmp1 * njac[i+1][3][1];
        lhs[i][CC][4][1] =  tmp2 * fjac[i+1][4][1]
          - tmp1 * njac[i+1][4][1];

        lhs[i][CC][0][2] =  tmp2 * fjac[i+1][0][2]
          - tmp1 * njac[i+1][0][2];
        lhs[i][CC][1][2] =  tmp2 * fjac[i+1][1][2]
          - tmp1 * njac[i+1][1][2];
        lhs[i][CC][2][2] =  tmp2 * fjac[i+1][2][2]
          - tmp1 * njac[i+1][2][2]
          - tmp1 * dx3;
        lhs[i][CC][3][2] =  tmp2 * fjac[i+1][3][2]
          - tmp1 * njac[i+1][3][2];
        lhs[i][CC][4][2] =  tmp2 * fjac[i+1][4][2]
          - tmp1 * njac[i+1][4][2];

        lhs[i][CC][0][3] =  tmp2 * fjac[i+1][0][3]
          - tmp1 * njac[i+1][0][3];
        lhs[i][CC][1][3] =  tmp2 * fjac[i+1][1][3]
          - tmp1 * njac[i+1][1][3];
        lhs[i][CC][2][3] =  tmp2 * fjac[i+1][2][3]
          - tmp1 * njac[i+1][2][3];
        lhs[i][CC][3][3] =  tmp2 * fjac[i+1][3][3]
          - tmp1 * njac[i+1][3][3]
          - tmp1 * dx4;
        lhs[i][CC][4][3] =  tmp2 * fjac[i+1][4][3]
          - tmp1 * njac[i+1][4][3];

        lhs[i][CC][0][4] =  tmp2 * fjac[i+1][0][4]
          - tmp1 * njac[i+1][0][4];
        lhs[i][CC][1][4] =  tmp2 * fjac[i+1][1][4]
          - tmp1 * njac[i+1][1][4];
        lhs[i][CC][2][4] =  tmp2 * fjac[i+1][2][4]
          - tmp1 * njac[i+1][2][4];
        lhs[i][CC][3][4] =  tmp2 * fjac[i+1][3][4]
          - tmp1 * njac[i+1][3][4];
        lhs[i][CC][4][4] =  tmp2 * fjac[i+1][4][4]
          - tmp1 * njac[i+1][4][4]
          - tmp1 * dx5;
      //}

      
	    if(i == 1){
		  //binvcrhs( lhs[0][BB], lhs[0][CC], rhs[k][j][0] );
		  binvcrhs( 0, BB, 0, CC, k, j, 0 );
	    }
      //---------------------------------------------------------------------
      // begin inner most do loop
      // do all the elements of the cell unless last 
      //---------------------------------------------------------------------
      //for (i = 1; i <= isize-1; i++) {
        //-------------------------------------------------------------------
        // rhs(i) = rhs(i) - A*rhs(i-1)
        //-------------------------------------------------------------------
		
        //matvec_sub(lhs[i][AA], rhs[k][j][i-1], rhs[k][j][i]);
		matvec_sub(i, AA, k, j, i-1, k, j, i);
		
        //-------------------------------------------------------------------
        // B(i) = B(i) - C(i-1)*A(i)
        //-------------------------------------------------------------------
        
		//matmul_sub(lhs[i][AA], lhs[i-1][CC], lhs[i][BB]);
		matmul_sub(i, AA, i-1, CC, i, BB);


        
        
		//binvcrhs( lhs[i][BB], lhs[i][CC], rhs[k][j][i] );
		binvcrhs( i, BB, i, CC, k, j, i );
		
		
		if(i == isize-1){
			//matvec_sub(lhs[isize][AA], rhs[k][j][isize-1], rhs[k][j][isize]);
			matvec_sub(isize, AA, k, j, isize-1, k, j, isize);
	  
			//matmul_sub(lhs[isize][AA], lhs[isize-1][CC], lhs[isize][BB]);
			matmul_sub(isize, AA, isize-1, CC, isize, BB);

			//binvrhs( lhs[isize][BB], rhs[k][j][isize] );
			binvrhs( isize, BB, k, j, isize );
		}
		
		
		
		for(m = 0; m < BLOCK_SIZE; m++){
			for(n = 0; n < BLOCK_SIZE; n++){				
				lhs_buf[k][j][i][n][m] = lhs[i][CC][n][m];
			}
		}
      }
    }
  }
      

      
  for (k = 1; k <= grid_points[2]-2; k++) {
    for (j = 1; j <= grid_points[1]-2; j++) {     
	  for (i = isize-1; i >=0; i--) {
        for (m = 0; m < BLOCK_SIZE; m++) {
          for (n = 0; n < BLOCK_SIZE; n++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] 
              - lhs_buf[k][j][i][n][m]*rhs[k][j][i+1][n];
          }
        }
      }
    }
  }
  
  
  
  
  
}

