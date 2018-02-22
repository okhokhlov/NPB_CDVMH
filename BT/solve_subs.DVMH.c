#include <cdvmh_helpers.h>

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
// subtracts bvec=bvec - ablock*avec
//---------------------------------------------------------------------
/*void matvec_sub(double ablock[5][5], double avec[5], double bvec[5])
{
  //---------------------------------------------------------------------
  // rhs[kc][jc][ic][i] = rhs[kc][jc][ic][i] 
  // $                  - lhs[ia][ablock][0][i]*
  //---------------------------------------------------------------------
  bvec[0] = bvec[0] - ablock[0][0]*avec[0]
                    - ablock[1][0]*avec[1]
                    - ablock[2][0]*avec[2]
                    - ablock[3][0]*avec[3]
                    - ablock[4][0]*avec[4];
  bvec[1] = bvec[1] - ablock[0][1]*avec[0]
                    - ablock[1][1]*avec[1]
                    - ablock[2][1]*avec[2]
                    - ablock[3][1]*avec[3]
                    - ablock[4][1]*avec[4];
  bvec[2] = bvec[2] - ablock[0][2]*avec[0]
                    - ablock[1][2]*avec[1]
                    - ablock[2][2]*avec[2]
                    - ablock[3][2]*avec[3]
                    - ablock[4][2]*avec[4];
  bvec[3] = bvec[3] - ablock[0][3]*avec[0]
                    - ablock[1][3]*avec[1]
                    - ablock[2][3]*avec[2]
                    - ablock[3][3]*avec[3]
                    - ablock[4][3]*avec[4];
  bvec[4] = bvec[4] - ablock[0][4]*avec[0]
                    - ablock[1][4]*avec[1]
                    - ablock[2][4]*avec[2]
                    - ablock[3][4]*avec[3]
                    - ablock[4][4]*avec[4];
}
*/

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

////////////////////////////////////////////
void matvec_sub(int a1, int a2, int b1, int b2, int b3, int c1, int c2, int c3 )
{
  //---------------------------------------------------------------------
  // rhs[kc][jc][ic][i] = rhs[kc][jc][ic][i] 
  // $                  - lhs[ia][lhs[a1][a2]][0][i]*
  //---------------------------------------------------------------------
  rhs[c1][c2][c3][0] = rhs[c1][c2][c3][0] - lhs[a1][a2][0][0]*rhs[b1][b2][b3][0]
                    - lhs[a1][a2][1][0]*rhs[b1][b2][b3][1]
                    - lhs[a1][a2][2][0]*rhs[b1][b2][b3][2]
                    - lhs[a1][a2][3][0]*rhs[b1][b2][b3][3]
                    - lhs[a1][a2][4][0]*rhs[b1][b2][b3][4];
  rhs[c1][c2][c3][1] = rhs[c1][c2][c3][1] - lhs[a1][a2][0][1]*rhs[b1][b2][b3][0]
                    - lhs[a1][a2][1][1]*rhs[b1][b2][b3][1]
                    - lhs[a1][a2][2][1]*rhs[b1][b2][b3][2]
                    - lhs[a1][a2][3][1]*rhs[b1][b2][b3][3]
                    - lhs[a1][a2][4][1]*rhs[b1][b2][b3][4];
  rhs[c1][c2][c3][2] = rhs[c1][c2][c3][2] - lhs[a1][a2][0][2]*rhs[b1][b2][b3][0]
                    - lhs[a1][a2][1][2]*rhs[b1][b2][b3][1]
                    - lhs[a1][a2][2][2]*rhs[b1][b2][b3][2]
                    - lhs[a1][a2][3][2]*rhs[b1][b2][b3][3]
                    - lhs[a1][a2][4][2]*rhs[b1][b2][b3][4];
  rhs[c1][c2][c3][3] = rhs[c1][c2][c3][3] - lhs[a1][a2][0][3]*rhs[b1][b2][b3][0]
                    - lhs[a1][a2][1][3]*rhs[b1][b2][b3][1]
                    - lhs[a1][a2][2][3]*rhs[b1][b2][b3][2]
                    - lhs[a1][a2][3][3]*rhs[b1][b2][b3][3]
                    - lhs[a1][a2][4][3]*rhs[b1][b2][b3][4];
  rhs[c1][c2][c3][4] = rhs[c1][c2][c3][4] - lhs[a1][a2][0][4]*rhs[b1][b2][b3][0]
                    - lhs[a1][a2][1][4]*rhs[b1][b2][b3][1]
                    - lhs[a1][a2][2][4]*rhs[b1][b2][b3][2]
                    - lhs[a1][a2][3][4]*rhs[b1][b2][b3][3]
                    - lhs[a1][a2][4][4]*rhs[b1][b2][b3][4];
}

////////////////////////////////////////////

//---------------------------------------------------------------------
// subtracts a(i,j,k) X b(i,j,k) from c(i,j,k)
//---------------------------------------------------------------------
/*void matmul_sub(double ablock[5][5], double bblock[5][5], double cblock[5][5])
{
  cblock[0][0] = cblock[0][0] - ablock[0][0]*bblock[0][0]
                              - ablock[1][0]*bblock[0][1]
                              - ablock[2][0]*bblock[0][2]
                              - ablock[3][0]*bblock[0][3]
                              - ablock[4][0]*bblock[0][4];
  cblock[0][1] = cblock[0][1] - ablock[0][1]*bblock[0][0]
                              - ablock[1][1]*bblock[0][1]
                              - ablock[2][1]*bblock[0][2]
                              - ablock[3][1]*bblock[0][3]
                              - ablock[4][1]*bblock[0][4];
  cblock[0][2] = cblock[0][2] - ablock[0][2]*bblock[0][0]
                              - ablock[1][2]*bblock[0][1]
                              - ablock[2][2]*bblock[0][2]
                              - ablock[3][2]*bblock[0][3]
                              - ablock[4][2]*bblock[0][4];
  cblock[0][3] = cblock[0][3] - ablock[0][3]*bblock[0][0]
                              - ablock[1][3]*bblock[0][1]
                              - ablock[2][3]*bblock[0][2]
                              - ablock[3][3]*bblock[0][3]
                              - ablock[4][3]*bblock[0][4];
  cblock[0][4] = cblock[0][4] - ablock[0][4]*bblock[0][0]
                              - ablock[1][4]*bblock[0][1]
                              - ablock[2][4]*bblock[0][2]
                              - ablock[3][4]*bblock[0][3]
                              - ablock[4][4]*bblock[0][4];
  cblock[1][0] = cblock[1][0] - ablock[0][0]*bblock[1][0]
                              - ablock[1][0]*bblock[1][1]
                              - ablock[2][0]*bblock[1][2]
                              - ablock[3][0]*bblock[1][3]
                              - ablock[4][0]*bblock[1][4];
  cblock[1][1] = cblock[1][1] - ablock[0][1]*bblock[1][0]
                              - ablock[1][1]*bblock[1][1]
                              - ablock[2][1]*bblock[1][2]
                              - ablock[3][1]*bblock[1][3]
                              - ablock[4][1]*bblock[1][4];
  cblock[1][2] = cblock[1][2] - ablock[0][2]*bblock[1][0]
                              - ablock[1][2]*bblock[1][1]
                              - ablock[2][2]*bblock[1][2]
                              - ablock[3][2]*bblock[1][3]
                              - ablock[4][2]*bblock[1][4];
  cblock[1][3] = cblock[1][3] - ablock[0][3]*bblock[1][0]
                              - ablock[1][3]*bblock[1][1]
                              - ablock[2][3]*bblock[1][2]
                              - ablock[3][3]*bblock[1][3]
                              - ablock[4][3]*bblock[1][4];
  cblock[1][4] = cblock[1][4] - ablock[0][4]*bblock[1][0]
                              - ablock[1][4]*bblock[1][1]
                              - ablock[2][4]*bblock[1][2]
                              - ablock[3][4]*bblock[1][3]
                              - ablock[4][4]*bblock[1][4];
  cblock[2][0] = cblock[2][0] - ablock[0][0]*bblock[2][0]
                              - ablock[1][0]*bblock[2][1]
                              - ablock[2][0]*bblock[2][2]
                              - ablock[3][0]*bblock[2][3]
                              - ablock[4][0]*bblock[2][4];
  cblock[2][1] = cblock[2][1] - ablock[0][1]*bblock[2][0]
                              - ablock[1][1]*bblock[2][1]
                              - ablock[2][1]*bblock[2][2]
                              - ablock[3][1]*bblock[2][3]
                              - ablock[4][1]*bblock[2][4];
  cblock[2][2] = cblock[2][2] - ablock[0][2]*bblock[2][0]
                              - ablock[1][2]*bblock[2][1]
                              - ablock[2][2]*bblock[2][2]
                              - ablock[3][2]*bblock[2][3]
                              - ablock[4][2]*bblock[2][4];
  cblock[2][3] = cblock[2][3] - ablock[0][3]*bblock[2][0]
                              - ablock[1][3]*bblock[2][1]
                              - ablock[2][3]*bblock[2][2]
                              - ablock[3][3]*bblock[2][3]
                              - ablock[4][3]*bblock[2][4];
  cblock[2][4] = cblock[2][4] - ablock[0][4]*bblock[2][0]
                              - ablock[1][4]*bblock[2][1]
                              - ablock[2][4]*bblock[2][2]
                              - ablock[3][4]*bblock[2][3]
                              - ablock[4][4]*bblock[2][4];
  cblock[3][0] = cblock[3][0] - ablock[0][0]*bblock[3][0]
                              - ablock[1][0]*bblock[3][1]
                              - ablock[2][0]*bblock[3][2]
                              - ablock[3][0]*bblock[3][3]
                              - ablock[4][0]*bblock[3][4];
  cblock[3][1] = cblock[3][1] - ablock[0][1]*bblock[3][0]
                              - ablock[1][1]*bblock[3][1]
                              - ablock[2][1]*bblock[3][2]
                              - ablock[3][1]*bblock[3][3]
                              - ablock[4][1]*bblock[3][4];
  cblock[3][2] = cblock[3][2] - ablock[0][2]*bblock[3][0]
                              - ablock[1][2]*bblock[3][1]
                              - ablock[2][2]*bblock[3][2]
                              - ablock[3][2]*bblock[3][3]
                              - ablock[4][2]*bblock[3][4];
  cblock[3][3] = cblock[3][3] - ablock[0][3]*bblock[3][0]
                              - ablock[1][3]*bblock[3][1]
                              - ablock[2][3]*bblock[3][2]
                              - ablock[3][3]*bblock[3][3]
                              - ablock[4][3]*bblock[3][4];
  cblock[3][4] = cblock[3][4] - ablock[0][4]*bblock[3][0]
                              - ablock[1][4]*bblock[3][1]
                              - ablock[2][4]*bblock[3][2]
                              - ablock[3][4]*bblock[3][3]
                              - ablock[4][4]*bblock[3][4];
  cblock[4][0] = cblock[4][0] - ablock[0][0]*bblock[4][0]
                              - ablock[1][0]*bblock[4][1]
                              - ablock[2][0]*bblock[4][2]
                              - ablock[3][0]*bblock[4][3]
                              - ablock[4][0]*bblock[4][4];
  cblock[4][1] = cblock[4][1] - ablock[0][1]*bblock[4][0]
                              - ablock[1][1]*bblock[4][1]
                              - ablock[2][1]*bblock[4][2]
                              - ablock[3][1]*bblock[4][3]
                              - ablock[4][1]*bblock[4][4];
  cblock[4][2] = cblock[4][2] - ablock[0][2]*bblock[4][0]
                              - ablock[1][2]*bblock[4][1]
                              - ablock[2][2]*bblock[4][2]
                              - ablock[3][2]*bblock[4][3]
                              - ablock[4][2]*bblock[4][4];
  cblock[4][3] = cblock[4][3] - ablock[0][3]*bblock[4][0]
                              - ablock[1][3]*bblock[4][1]
                              - ablock[2][3]*bblock[4][2]
                              - ablock[3][3]*bblock[4][3]
                              - ablock[4][3]*bblock[4][4];
  cblock[4][4] = cblock[4][4] - ablock[0][4]*bblock[4][0]
                              - ablock[1][4]*bblock[4][1]
                              - ablock[2][4]*bblock[4][2]
                              - ablock[3][4]*bblock[4][3]
                              - ablock[4][4]*bblock[4][4];
}
*/

///////////////////////////////////////

void matmul_sub(int a1, int a2, int b1, int b2, int c1, int c2 )
{
  lhs[c1][c2][0][0] = lhs[c1][c2][0][0] - lhs[a1][a2][0][0]*lhs[b1][b2][0][0]
                              - lhs[a1][a2][1][0]*lhs[b1][b2][0][1]
                              - lhs[a1][a2][2][0]*lhs[b1][b2][0][2]
                              - lhs[a1][a2][3][0]*lhs[b1][b2][0][3]
                              - lhs[a1][a2][4][0]*lhs[b1][b2][0][4];
  lhs[c1][c2][0][1] = lhs[c1][c2][0][1] - lhs[a1][a2][0][1]*lhs[b1][b2][0][0]
                              - lhs[a1][a2][1][1]*lhs[b1][b2][0][1]
                              - lhs[a1][a2][2][1]*lhs[b1][b2][0][2]
                              - lhs[a1][a2][3][1]*lhs[b1][b2][0][3]
                              - lhs[a1][a2][4][1]*lhs[b1][b2][0][4];
  lhs[c1][c2][0][2] = lhs[c1][c2][0][2] - lhs[a1][a2][0][2]*lhs[b1][b2][0][0]
                              - lhs[a1][a2][1][2]*lhs[b1][b2][0][1]
                              - lhs[a1][a2][2][2]*lhs[b1][b2][0][2]
                              - lhs[a1][a2][3][2]*lhs[b1][b2][0][3]
                              - lhs[a1][a2][4][2]*lhs[b1][b2][0][4];
  lhs[c1][c2][0][3] = lhs[c1][c2][0][3] - lhs[a1][a2][0][3]*lhs[b1][b2][0][0]
                              - lhs[a1][a2][1][3]*lhs[b1][b2][0][1]
                              - lhs[a1][a2][2][3]*lhs[b1][b2][0][2]
                              - lhs[a1][a2][3][3]*lhs[b1][b2][0][3]
                              - lhs[a1][a2][4][3]*lhs[b1][b2][0][4];
  lhs[c1][c2][0][4] = lhs[c1][c2][0][4] - lhs[a1][a2][0][4]*lhs[b1][b2][0][0]
                              - lhs[a1][a2][1][4]*lhs[b1][b2][0][1]
                              - lhs[a1][a2][2][4]*lhs[b1][b2][0][2]
                              - lhs[a1][a2][3][4]*lhs[b1][b2][0][3]
                              - lhs[a1][a2][4][4]*lhs[b1][b2][0][4];
  lhs[c1][c2][1][0] = lhs[c1][c2][1][0] - lhs[a1][a2][0][0]*lhs[b1][b2][1][0]
                              - lhs[a1][a2][1][0]*lhs[b1][b2][1][1]
                              - lhs[a1][a2][2][0]*lhs[b1][b2][1][2]
                              - lhs[a1][a2][3][0]*lhs[b1][b2][1][3]
                              - lhs[a1][a2][4][0]*lhs[b1][b2][1][4];
  lhs[c1][c2][1][1] = lhs[c1][c2][1][1] - lhs[a1][a2][0][1]*lhs[b1][b2][1][0]
                              - lhs[a1][a2][1][1]*lhs[b1][b2][1][1]
                              - lhs[a1][a2][2][1]*lhs[b1][b2][1][2]
                              - lhs[a1][a2][3][1]*lhs[b1][b2][1][3]
                              - lhs[a1][a2][4][1]*lhs[b1][b2][1][4];
  lhs[c1][c2][1][2] = lhs[c1][c2][1][2] - lhs[a1][a2][0][2]*lhs[b1][b2][1][0]
                              - lhs[a1][a2][1][2]*lhs[b1][b2][1][1]
                              - lhs[a1][a2][2][2]*lhs[b1][b2][1][2]
                              - lhs[a1][a2][3][2]*lhs[b1][b2][1][3]
                              - lhs[a1][a2][4][2]*lhs[b1][b2][1][4];
  lhs[c1][c2][1][3] = lhs[c1][c2][1][3] - lhs[a1][a2][0][3]*lhs[b1][b2][1][0]
                              - lhs[a1][a2][1][3]*lhs[b1][b2][1][1]
                              - lhs[a1][a2][2][3]*lhs[b1][b2][1][2]
                              - lhs[a1][a2][3][3]*lhs[b1][b2][1][3]
                              - lhs[a1][a2][4][3]*lhs[b1][b2][1][4];
  lhs[c1][c2][1][4] = lhs[c1][c2][1][4] - lhs[a1][a2][0][4]*lhs[b1][b2][1][0]
                              - lhs[a1][a2][1][4]*lhs[b1][b2][1][1]
                              - lhs[a1][a2][2][4]*lhs[b1][b2][1][2]
                              - lhs[a1][a2][3][4]*lhs[b1][b2][1][3]
                              - lhs[a1][a2][4][4]*lhs[b1][b2][1][4];
  lhs[c1][c2][2][0] = lhs[c1][c2][2][0] - lhs[a1][a2][0][0]*lhs[b1][b2][2][0]
                              - lhs[a1][a2][1][0]*lhs[b1][b2][2][1]
                              - lhs[a1][a2][2][0]*lhs[b1][b2][2][2]
                              - lhs[a1][a2][3][0]*lhs[b1][b2][2][3]
                              - lhs[a1][a2][4][0]*lhs[b1][b2][2][4];
  lhs[c1][c2][2][1] = lhs[c1][c2][2][1] - lhs[a1][a2][0][1]*lhs[b1][b2][2][0]
                              - lhs[a1][a2][1][1]*lhs[b1][b2][2][1]
                              - lhs[a1][a2][2][1]*lhs[b1][b2][2][2]
                              - lhs[a1][a2][3][1]*lhs[b1][b2][2][3]
                              - lhs[a1][a2][4][1]*lhs[b1][b2][2][4];
  lhs[c1][c2][2][2] = lhs[c1][c2][2][2] - lhs[a1][a2][0][2]*lhs[b1][b2][2][0]
                              - lhs[a1][a2][1][2]*lhs[b1][b2][2][1]
                              - lhs[a1][a2][2][2]*lhs[b1][b2][2][2]
                              - lhs[a1][a2][3][2]*lhs[b1][b2][2][3]
                              - lhs[a1][a2][4][2]*lhs[b1][b2][2][4];
  lhs[c1][c2][2][3] = lhs[c1][c2][2][3] - lhs[a1][a2][0][3]*lhs[b1][b2][2][0]
                              - lhs[a1][a2][1][3]*lhs[b1][b2][2][1]
                              - lhs[a1][a2][2][3]*lhs[b1][b2][2][2]
                              - lhs[a1][a2][3][3]*lhs[b1][b2][2][3]
                              - lhs[a1][a2][4][3]*lhs[b1][b2][2][4];
  lhs[c1][c2][2][4] = lhs[c1][c2][2][4] - lhs[a1][a2][0][4]*lhs[b1][b2][2][0]
                              - lhs[a1][a2][1][4]*lhs[b1][b2][2][1]
                              - lhs[a1][a2][2][4]*lhs[b1][b2][2][2]
                              - lhs[a1][a2][3][4]*lhs[b1][b2][2][3]
                              - lhs[a1][a2][4][4]*lhs[b1][b2][2][4];
  lhs[c1][c2][3][0] = lhs[c1][c2][3][0] - lhs[a1][a2][0][0]*lhs[b1][b2][3][0]
                              - lhs[a1][a2][1][0]*lhs[b1][b2][3][1]
                              - lhs[a1][a2][2][0]*lhs[b1][b2][3][2]
                              - lhs[a1][a2][3][0]*lhs[b1][b2][3][3]
                              - lhs[a1][a2][4][0]*lhs[b1][b2][3][4];
  lhs[c1][c2][3][1] = lhs[c1][c2][3][1] - lhs[a1][a2][0][1]*lhs[b1][b2][3][0]
                              - lhs[a1][a2][1][1]*lhs[b1][b2][3][1]
                              - lhs[a1][a2][2][1]*lhs[b1][b2][3][2]
                              - lhs[a1][a2][3][1]*lhs[b1][b2][3][3]
                              - lhs[a1][a2][4][1]*lhs[b1][b2][3][4];
  lhs[c1][c2][3][2] = lhs[c1][c2][3][2] - lhs[a1][a2][0][2]*lhs[b1][b2][3][0]
                              - lhs[a1][a2][1][2]*lhs[b1][b2][3][1]
                              - lhs[a1][a2][2][2]*lhs[b1][b2][3][2]
                              - lhs[a1][a2][3][2]*lhs[b1][b2][3][3]
                              - lhs[a1][a2][4][2]*lhs[b1][b2][3][4];
  lhs[c1][c2][3][3] = lhs[c1][c2][3][3] - lhs[a1][a2][0][3]*lhs[b1][b2][3][0]
                              - lhs[a1][a2][1][3]*lhs[b1][b2][3][1]
                              - lhs[a1][a2][2][3]*lhs[b1][b2][3][2]
                              - lhs[a1][a2][3][3]*lhs[b1][b2][3][3]
                              - lhs[a1][a2][4][3]*lhs[b1][b2][3][4];
  lhs[c1][c2][3][4] = lhs[c1][c2][3][4] - lhs[a1][a2][0][4]*lhs[b1][b2][3][0]
                              - lhs[a1][a2][1][4]*lhs[b1][b2][3][1]
                              - lhs[a1][a2][2][4]*lhs[b1][b2][3][2]
                              - lhs[a1][a2][3][4]*lhs[b1][b2][3][3]
                              - lhs[a1][a2][4][4]*lhs[b1][b2][3][4];
  lhs[c1][c2][4][0] = lhs[c1][c2][4][0] - lhs[a1][a2][0][0]*lhs[b1][b2][4][0]
                              - lhs[a1][a2][1][0]*lhs[b1][b2][4][1]
                              - lhs[a1][a2][2][0]*lhs[b1][b2][4][2]
                              - lhs[a1][a2][3][0]*lhs[b1][b2][4][3]
                              - lhs[a1][a2][4][0]*lhs[b1][b2][4][4];
  lhs[c1][c2][4][1] = lhs[c1][c2][4][1] - lhs[a1][a2][0][1]*lhs[b1][b2][4][0]
                              - lhs[a1][a2][1][1]*lhs[b1][b2][4][1]
                              - lhs[a1][a2][2][1]*lhs[b1][b2][4][2]
                              - lhs[a1][a2][3][1]*lhs[b1][b2][4][3]
                              - lhs[a1][a2][4][1]*lhs[b1][b2][4][4];
  lhs[c1][c2][4][2] = lhs[c1][c2][4][2] - lhs[a1][a2][0][2]*lhs[b1][b2][4][0]
                              - lhs[a1][a2][1][2]*lhs[b1][b2][4][1]
                              - lhs[a1][a2][2][2]*lhs[b1][b2][4][2]
                              - lhs[a1][a2][3][2]*lhs[b1][b2][4][3]
                              - lhs[a1][a2][4][2]*lhs[b1][b2][4][4];
  lhs[c1][c2][4][3] = lhs[c1][c2][4][3] - lhs[a1][a2][0][3]*lhs[b1][b2][4][0]
                              - lhs[a1][a2][1][3]*lhs[b1][b2][4][1]
                              - lhs[a1][a2][2][3]*lhs[b1][b2][4][2]
                              - lhs[a1][a2][3][3]*lhs[b1][b2][4][3]
                              - lhs[a1][a2][4][3]*lhs[b1][b2][4][4];
  lhs[c1][c2][4][4] = lhs[c1][c2][4][4] - lhs[a1][a2][0][4]*lhs[b1][b2][4][0]
                              - lhs[a1][a2][1][4]*lhs[b1][b2][4][1]
                              - lhs[a1][a2][2][4]*lhs[b1][b2][4][2]
                              - lhs[a1][a2][3][4]*lhs[b1][b2][4][3]
                              - lhs[a1][a2][4][4]*lhs[b1][b2][4][4];
}

///////////////////////////////////////

/*void binvcrhs(double lhs[5][5], double c[5][5], double r[5])
{
  double pivot, coeff;

  pivot = 1.00/lhs[0][0];
  lhs[1][0] = lhs[1][0]*pivot;
  lhs[2][0] = lhs[2][0]*pivot;
  lhs[3][0] = lhs[3][0]*pivot;
  lhs[4][0] = lhs[4][0]*pivot;
  c[0][0] = c[0][0]*pivot;
  c[1][0] = c[1][0]*pivot;
  c[2][0] = c[2][0]*pivot;
  c[3][0] = c[3][0]*pivot;
  c[4][0] = c[4][0]*pivot;
  r[0]   = r[0]  *pivot;

  coeff = lhs[0][1];
  lhs[1][1]= lhs[1][1] - coeff*lhs[1][0];
  lhs[2][1]= lhs[2][1] - coeff*lhs[2][0];
  lhs[3][1]= lhs[3][1] - coeff*lhs[3][0];
  lhs[4][1]= lhs[4][1] - coeff*lhs[4][0];
  c[0][1] = c[0][1] - coeff*c[0][0];
  c[1][1] = c[1][1] - coeff*c[1][0];
  c[2][1] = c[2][1] - coeff*c[2][0];
  c[3][1] = c[3][1] - coeff*c[3][0];
  c[4][1] = c[4][1] - coeff*c[4][0];
  r[1]   = r[1]   - coeff*r[0];

  coeff = lhs[0][2];
  lhs[1][2]= lhs[1][2] - coeff*lhs[1][0];
  lhs[2][2]= lhs[2][2] - coeff*lhs[2][0];
  lhs[3][2]= lhs[3][2] - coeff*lhs[3][0];
  lhs[4][2]= lhs[4][2] - coeff*lhs[4][0];
  c[0][2] = c[0][2] - coeff*c[0][0];
  c[1][2] = c[1][2] - coeff*c[1][0];
  c[2][2] = c[2][2] - coeff*c[2][0];
  c[3][2] = c[3][2] - coeff*c[3][0];
  c[4][2] = c[4][2] - coeff*c[4][0];
  r[2]   = r[2]   - coeff*r[0];

  coeff = lhs[0][3];
  lhs[1][3]= lhs[1][3] - coeff*lhs[1][0];
  lhs[2][3]= lhs[2][3] - coeff*lhs[2][0];
  lhs[3][3]= lhs[3][3] - coeff*lhs[3][0];
  lhs[4][3]= lhs[4][3] - coeff*lhs[4][0];
  c[0][3] = c[0][3] - coeff*c[0][0];
  c[1][3] = c[1][3] - coeff*c[1][0];
  c[2][3] = c[2][3] - coeff*c[2][0];
  c[3][3] = c[3][3] - coeff*c[3][0];
  c[4][3] = c[4][3] - coeff*c[4][0];
  r[3]   = r[3]   - coeff*r[0];

  coeff = lhs[0][4];
  lhs[1][4]= lhs[1][4] - coeff*lhs[1][0];
  lhs[2][4]= lhs[2][4] - coeff*lhs[2][0];
  lhs[3][4]= lhs[3][4] - coeff*lhs[3][0];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][0];
  c[0][4] = c[0][4] - coeff*c[0][0];
  c[1][4] = c[1][4] - coeff*c[1][0];
  c[2][4] = c[2][4] - coeff*c[2][0];
  c[3][4] = c[3][4] - coeff*c[3][0];
  c[4][4] = c[4][4] - coeff*c[4][0];
  r[4]   = r[4]   - coeff*r[0];


  pivot = 1.00/lhs[1][1];
  lhs[2][1] = lhs[2][1]*pivot;
  lhs[3][1] = lhs[3][1]*pivot;
  lhs[4][1] = lhs[4][1]*pivot;
  c[0][1] = c[0][1]*pivot;
  c[1][1] = c[1][1]*pivot;
  c[2][1] = c[2][1]*pivot;
  c[3][1] = c[3][1]*pivot;
  c[4][1] = c[4][1]*pivot;
  r[1]   = r[1]  *pivot;

  coeff = lhs[1][0];
  lhs[2][0]= lhs[2][0] - coeff*lhs[2][1];
  lhs[3][0]= lhs[3][0] - coeff*lhs[3][1];
  lhs[4][0]= lhs[4][0] - coeff*lhs[4][1];
  c[0][0] = c[0][0] - coeff*c[0][1];
  c[1][0] = c[1][0] - coeff*c[1][1];
  c[2][0] = c[2][0] - coeff*c[2][1];
  c[3][0] = c[3][0] - coeff*c[3][1];
  c[4][0] = c[4][0] - coeff*c[4][1];
  r[0]   = r[0]   - coeff*r[1];

  coeff = lhs[1][2];
  lhs[2][2]= lhs[2][2] - coeff*lhs[2][1];
  lhs[3][2]= lhs[3][2] - coeff*lhs[3][1];
  lhs[4][2]= lhs[4][2] - coeff*lhs[4][1];
  c[0][2] = c[0][2] - coeff*c[0][1];
  c[1][2] = c[1][2] - coeff*c[1][1];
  c[2][2] = c[2][2] - coeff*c[2][1];
  c[3][2] = c[3][2] - coeff*c[3][1];
  c[4][2] = c[4][2] - coeff*c[4][1];
  r[2]   = r[2]   - coeff*r[1];

  coeff = lhs[1][3];
  lhs[2][3]= lhs[2][3] - coeff*lhs[2][1];
  lhs[3][3]= lhs[3][3] - coeff*lhs[3][1];
  lhs[4][3]= lhs[4][3] - coeff*lhs[4][1];
  c[0][3] = c[0][3] - coeff*c[0][1];
  c[1][3] = c[1][3] - coeff*c[1][1];
  c[2][3] = c[2][3] - coeff*c[2][1];
  c[3][3] = c[3][3] - coeff*c[3][1];
  c[4][3] = c[4][3] - coeff*c[4][1];
  r[3]   = r[3]   - coeff*r[1];

  coeff = lhs[1][4];
  lhs[2][4]= lhs[2][4] - coeff*lhs[2][1];
  lhs[3][4]= lhs[3][4] - coeff*lhs[3][1];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][1];
  c[0][4] = c[0][4] - coeff*c[0][1];
  c[1][4] = c[1][4] - coeff*c[1][1];
  c[2][4] = c[2][4] - coeff*c[2][1];
  c[3][4] = c[3][4] - coeff*c[3][1];
  c[4][4] = c[4][4] - coeff*c[4][1];
  r[4]   = r[4]   - coeff*r[1];


  pivot = 1.00/lhs[2][2];
  lhs[3][2] = lhs[3][2]*pivot;
  lhs[4][2] = lhs[4][2]*pivot;
  c[0][2] = c[0][2]*pivot;
  c[1][2] = c[1][2]*pivot;
  c[2][2] = c[2][2]*pivot;
  c[3][2] = c[3][2]*pivot;
  c[4][2] = c[4][2]*pivot;
  r[2]   = r[2]  *pivot;

  coeff = lhs[2][0];
  lhs[3][0]= lhs[3][0] - coeff*lhs[3][2];
  lhs[4][0]= lhs[4][0] - coeff*lhs[4][2];
  c[0][0] = c[0][0] - coeff*c[0][2];
  c[1][0] = c[1][0] - coeff*c[1][2];
  c[2][0] = c[2][0] - coeff*c[2][2];
  c[3][0] = c[3][0] - coeff*c[3][2];
  c[4][0] = c[4][0] - coeff*c[4][2];
  r[0]   = r[0]   - coeff*r[2];

  coeff = lhs[2][1];
  lhs[3][1]= lhs[3][1] - coeff*lhs[3][2];
  lhs[4][1]= lhs[4][1] - coeff*lhs[4][2];
  c[0][1] = c[0][1] - coeff*c[0][2];
  c[1][1] = c[1][1] - coeff*c[1][2];
  c[2][1] = c[2][1] - coeff*c[2][2];
  c[3][1] = c[3][1] - coeff*c[3][2];
  c[4][1] = c[4][1] - coeff*c[4][2];
  r[1]   = r[1]   - coeff*r[2];

  coeff = lhs[2][3];
  lhs[3][3]= lhs[3][3] - coeff*lhs[3][2];
  lhs[4][3]= lhs[4][3] - coeff*lhs[4][2];
  c[0][3] = c[0][3] - coeff*c[0][2];
  c[1][3] = c[1][3] - coeff*c[1][2];
  c[2][3] = c[2][3] - coeff*c[2][2];
  c[3][3] = c[3][3] - coeff*c[3][2];
  c[4][3] = c[4][3] - coeff*c[4][2];
  r[3]   = r[3]   - coeff*r[2];

  coeff = lhs[2][4];
  lhs[3][4]= lhs[3][4] - coeff*lhs[3][2];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][2];
  c[0][4] = c[0][4] - coeff*c[0][2];
  c[1][4] = c[1][4] - coeff*c[1][2];
  c[2][4] = c[2][4] - coeff*c[2][2];
  c[3][4] = c[3][4] - coeff*c[3][2];
  c[4][4] = c[4][4] - coeff*c[4][2];
  r[4]   = r[4]   - coeff*r[2];


  pivot = 1.00/lhs[3][3];
  lhs[4][3] = lhs[4][3]*pivot;
  c[0][3] = c[0][3]*pivot;
  c[1][3] = c[1][3]*pivot;
  c[2][3] = c[2][3]*pivot;
  c[3][3] = c[3][3]*pivot;
  c[4][3] = c[4][3]*pivot;
  r[3]   = r[3]  *pivot;

  coeff = lhs[3][0];
  lhs[4][0]= lhs[4][0] - coeff*lhs[4][3];
  c[0][0] = c[0][0] - coeff*c[0][3];
  c[1][0] = c[1][0] - coeff*c[1][3];
  c[2][0] = c[2][0] - coeff*c[2][3];
  c[3][0] = c[3][0] - coeff*c[3][3];
  c[4][0] = c[4][0] - coeff*c[4][3];
  r[0]   = r[0]   - coeff*r[3];

  coeff = lhs[3][1];
  lhs[4][1]= lhs[4][1] - coeff*lhs[4][3];
  c[0][1] = c[0][1] - coeff*c[0][3];
  c[1][1] = c[1][1] - coeff*c[1][3];
  c[2][1] = c[2][1] - coeff*c[2][3];
  c[3][1] = c[3][1] - coeff*c[3][3];
  c[4][1] = c[4][1] - coeff*c[4][3];
  r[1]   = r[1]   - coeff*r[3];

  coeff = lhs[3][2];
  lhs[4][2]= lhs[4][2] - coeff*lhs[4][3];
  c[0][2] = c[0][2] - coeff*c[0][3];
  c[1][2] = c[1][2] - coeff*c[1][3];
  c[2][2] = c[2][2] - coeff*c[2][3];
  c[3][2] = c[3][2] - coeff*c[3][3];
  c[4][2] = c[4][2] - coeff*c[4][3];
  r[2]   = r[2]   - coeff*r[3];

  coeff = lhs[3][4];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][3];
  c[0][4] = c[0][4] - coeff*c[0][3];
  c[1][4] = c[1][4] - coeff*c[1][3];
  c[2][4] = c[2][4] - coeff*c[2][3];
  c[3][4] = c[3][4] - coeff*c[3][3];
  c[4][4] = c[4][4] - coeff*c[4][3];
  r[4]   = r[4]   - coeff*r[3];


  pivot = 1.00/lhs[4][4];
  c[0][4] = c[0][4]*pivot;
  c[1][4] = c[1][4]*pivot;
  c[2][4] = c[2][4]*pivot;
  c[3][4] = c[3][4]*pivot;
  c[4][4] = c[4][4]*pivot;
  r[4]   = r[4]  *pivot;

  coeff = lhs[4][0];
  c[0][0] = c[0][0] - coeff*c[0][4];
  c[1][0] = c[1][0] - coeff*c[1][4];
  c[2][0] = c[2][0] - coeff*c[2][4];
  c[3][0] = c[3][0] - coeff*c[3][4];
  c[4][0] = c[4][0] - coeff*c[4][4];
  r[0]   = r[0]   - coeff*r[4];

  coeff = lhs[4][1];
  c[0][1] = c[0][1] - coeff*c[0][4];
  c[1][1] = c[1][1] - coeff*c[1][4];
  c[2][1] = c[2][1] - coeff*c[2][4];
  c[3][1] = c[3][1] - coeff*c[3][4];
  c[4][1] = c[4][1] - coeff*c[4][4];
  r[1]   = r[1]   - coeff*r[4];

  coeff = lhs[4][2];
  c[0][2] = c[0][2] - coeff*c[0][4];
  c[1][2] = c[1][2] - coeff*c[1][4];
  c[2][2] = c[2][2] - coeff*c[2][4];
  c[3][2] = c[3][2] - coeff*c[3][4];
  c[4][2] = c[4][2] - coeff*c[4][4];
  r[2]   = r[2]   - coeff*r[4];

  coeff = lhs[4][3];
  c[0][3] = c[0][3] - coeff*c[0][4];
  c[1][3] = c[1][3] - coeff*c[1][4];
  c[2][3] = c[2][3] - coeff*c[2][4];
  c[3][3] = c[3][3] - coeff*c[3][4];
  c[4][3] = c[4][3] - coeff*c[4][4];
  r[3]   = r[3]   - coeff*r[4];
}
*/

////////////////////////////////////

void binvcrhs(int a1, int a2, int b1, int b2, int c1, int c2, int c3)
{
  double pivot, coeff;

  pivot = 1.00/lhs[a1][a2][0][0];
  lhs[a1][a2][1][0] = lhs[a1][a2][1][0]*pivot;
  lhs[a1][a2][2][0] = lhs[a1][a2][2][0]*pivot;
  lhs[a1][a2][3][0] = lhs[a1][a2][3][0]*pivot;
  lhs[a1][a2][4][0] = lhs[a1][a2][4][0]*pivot;
  lhs[b1][b2][0][0] = lhs[b1][b2][0][0]*pivot;
  lhs[b1][b2][1][0] = lhs[b1][b2][1][0]*pivot;
  lhs[b1][b2][2][0] = lhs[b1][b2][2][0]*pivot;
  lhs[b1][b2][3][0] = lhs[b1][b2][3][0]*pivot;
  lhs[b1][b2][4][0] = lhs[b1][b2][4][0]*pivot;
  rhs[c1][c2][c3][0]   = rhs[c1][c2][c3][0]  *pivot;

  coeff = lhs[a1][a2][0][1];
  lhs[a1][a2][1][1]= lhs[a1][a2][1][1] - coeff*lhs[a1][a2][1][0];
  lhs[a1][a2][2][1]= lhs[a1][a2][2][1] - coeff*lhs[a1][a2][2][0];
  lhs[a1][a2][3][1]= lhs[a1][a2][3][1] - coeff*lhs[a1][a2][3][0];
  lhs[a1][a2][4][1]= lhs[a1][a2][4][1] - coeff*lhs[a1][a2][4][0];
  lhs[b1][b2][0][1] = lhs[b1][b2][0][1] - coeff*lhs[b1][b2][0][0];
  lhs[b1][b2][1][1] = lhs[b1][b2][1][1] - coeff*lhs[b1][b2][1][0];
  lhs[b1][b2][2][1] = lhs[b1][b2][2][1] - coeff*lhs[b1][b2][2][0];
  lhs[b1][b2][3][1] = lhs[b1][b2][3][1] - coeff*lhs[b1][b2][3][0];
  lhs[b1][b2][4][1] = lhs[b1][b2][4][1] - coeff*lhs[b1][b2][4][0];
  rhs[c1][c2][c3][1]   = rhs[c1][c2][c3][1]   - coeff*rhs[c1][c2][c3][0];

  coeff = lhs[a1][a2][0][2];
  lhs[a1][a2][1][2]= lhs[a1][a2][1][2] - coeff*lhs[a1][a2][1][0];
  lhs[a1][a2][2][2]= lhs[a1][a2][2][2] - coeff*lhs[a1][a2][2][0];
  lhs[a1][a2][3][2]= lhs[a1][a2][3][2] - coeff*lhs[a1][a2][3][0];
  lhs[a1][a2][4][2]= lhs[a1][a2][4][2] - coeff*lhs[a1][a2][4][0];
  lhs[b1][b2][0][2] = lhs[b1][b2][0][2] - coeff*lhs[b1][b2][0][0];
  lhs[b1][b2][1][2] = lhs[b1][b2][1][2] - coeff*lhs[b1][b2][1][0];
  lhs[b1][b2][2][2] = lhs[b1][b2][2][2] - coeff*lhs[b1][b2][2][0];
  lhs[b1][b2][3][2] = lhs[b1][b2][3][2] - coeff*lhs[b1][b2][3][0];
  lhs[b1][b2][4][2] = lhs[b1][b2][4][2] - coeff*lhs[b1][b2][4][0];
  rhs[c1][c2][c3][2]   = rhs[c1][c2][c3][2]   - coeff*rhs[c1][c2][c3][0];

  coeff = lhs[a1][a2][0][3];
  lhs[a1][a2][1][3]= lhs[a1][a2][1][3] - coeff*lhs[a1][a2][1][0];
  lhs[a1][a2][2][3]= lhs[a1][a2][2][3] - coeff*lhs[a1][a2][2][0];
  lhs[a1][a2][3][3]= lhs[a1][a2][3][3] - coeff*lhs[a1][a2][3][0];
  lhs[a1][a2][4][3]= lhs[a1][a2][4][3] - coeff*lhs[a1][a2][4][0];
  lhs[b1][b2][0][3] = lhs[b1][b2][0][3] - coeff*lhs[b1][b2][0][0];
  lhs[b1][b2][1][3] = lhs[b1][b2][1][3] - coeff*lhs[b1][b2][1][0];
  lhs[b1][b2][2][3] = lhs[b1][b2][2][3] - coeff*lhs[b1][b2][2][0];
  lhs[b1][b2][3][3] = lhs[b1][b2][3][3] - coeff*lhs[b1][b2][3][0];
  lhs[b1][b2][4][3] = lhs[b1][b2][4][3] - coeff*lhs[b1][b2][4][0];
  rhs[c1][c2][c3][3]   = rhs[c1][c2][c3][3]   - coeff*rhs[c1][c2][c3][0];

  coeff = lhs[a1][a2][0][4];
  lhs[a1][a2][1][4]= lhs[a1][a2][1][4] - coeff*lhs[a1][a2][1][0];
  lhs[a1][a2][2][4]= lhs[a1][a2][2][4] - coeff*lhs[a1][a2][2][0];
  lhs[a1][a2][3][4]= lhs[a1][a2][3][4] - coeff*lhs[a1][a2][3][0];
  lhs[a1][a2][4][4]= lhs[a1][a2][4][4] - coeff*lhs[a1][a2][4][0];
  lhs[b1][b2][0][4] = lhs[b1][b2][0][4] - coeff*lhs[b1][b2][0][0];
  lhs[b1][b2][1][4] = lhs[b1][b2][1][4] - coeff*lhs[b1][b2][1][0];
  lhs[b1][b2][2][4] = lhs[b1][b2][2][4] - coeff*lhs[b1][b2][2][0];
  lhs[b1][b2][3][4] = lhs[b1][b2][3][4] - coeff*lhs[b1][b2][3][0];
  lhs[b1][b2][4][4] = lhs[b1][b2][4][4] - coeff*lhs[b1][b2][4][0];
  rhs[c1][c2][c3][4]   = rhs[c1][c2][c3][4]   - coeff*rhs[c1][c2][c3][0];


  pivot = 1.00/lhs[a1][a2][1][1];
  lhs[a1][a2][2][1] = lhs[a1][a2][2][1]*pivot;
  lhs[a1][a2][3][1] = lhs[a1][a2][3][1]*pivot;
  lhs[a1][a2][4][1] = lhs[a1][a2][4][1]*pivot;
  lhs[b1][b2][0][1] = lhs[b1][b2][0][1]*pivot;
  lhs[b1][b2][1][1] = lhs[b1][b2][1][1]*pivot;
  lhs[b1][b2][2][1] = lhs[b1][b2][2][1]*pivot;
  lhs[b1][b2][3][1] = lhs[b1][b2][3][1]*pivot;
  lhs[b1][b2][4][1] = lhs[b1][b2][4][1]*pivot;
  rhs[c1][c2][c3][1]   = rhs[c1][c2][c3][1]  *pivot;

  coeff = lhs[a1][a2][1][0];
  lhs[a1][a2][2][0]= lhs[a1][a2][2][0] - coeff*lhs[a1][a2][2][1];
  lhs[a1][a2][3][0]= lhs[a1][a2][3][0] - coeff*lhs[a1][a2][3][1];
  lhs[a1][a2][4][0]= lhs[a1][a2][4][0] - coeff*lhs[a1][a2][4][1];
  lhs[b1][b2][0][0] = lhs[b1][b2][0][0] - coeff*lhs[b1][b2][0][1];
  lhs[b1][b2][1][0] = lhs[b1][b2][1][0] - coeff*lhs[b1][b2][1][1];
  lhs[b1][b2][2][0] = lhs[b1][b2][2][0] - coeff*lhs[b1][b2][2][1];
  lhs[b1][b2][3][0] = lhs[b1][b2][3][0] - coeff*lhs[b1][b2][3][1];
  lhs[b1][b2][4][0] = lhs[b1][b2][4][0] - coeff*lhs[b1][b2][4][1];
  rhs[c1][c2][c3][0]   = rhs[c1][c2][c3][0]   - coeff*rhs[c1][c2][c3][1];

  coeff = lhs[a1][a2][1][2];
  lhs[a1][a2][2][2]= lhs[a1][a2][2][2] - coeff*lhs[a1][a2][2][1];
  lhs[a1][a2][3][2]= lhs[a1][a2][3][2] - coeff*lhs[a1][a2][3][1];
  lhs[a1][a2][4][2]= lhs[a1][a2][4][2] - coeff*lhs[a1][a2][4][1];
  lhs[b1][b2][0][2] = lhs[b1][b2][0][2] - coeff*lhs[b1][b2][0][1];
  lhs[b1][b2][1][2] = lhs[b1][b2][1][2] - coeff*lhs[b1][b2][1][1];
  lhs[b1][b2][2][2] = lhs[b1][b2][2][2] - coeff*lhs[b1][b2][2][1];
  lhs[b1][b2][3][2] = lhs[b1][b2][3][2] - coeff*lhs[b1][b2][3][1];
  lhs[b1][b2][4][2] = lhs[b1][b2][4][2] - coeff*lhs[b1][b2][4][1];
  rhs[c1][c2][c3][2]   = rhs[c1][c2][c3][2]   - coeff*rhs[c1][c2][c3][1];

  coeff = lhs[a1][a2][1][3];
  lhs[a1][a2][2][3]= lhs[a1][a2][2][3] - coeff*lhs[a1][a2][2][1];
  lhs[a1][a2][3][3]= lhs[a1][a2][3][3] - coeff*lhs[a1][a2][3][1];
  lhs[a1][a2][4][3]= lhs[a1][a2][4][3] - coeff*lhs[a1][a2][4][1];
  lhs[b1][b2][0][3] = lhs[b1][b2][0][3] - coeff*lhs[b1][b2][0][1];
  lhs[b1][b2][1][3] = lhs[b1][b2][1][3] - coeff*lhs[b1][b2][1][1];
  lhs[b1][b2][2][3] = lhs[b1][b2][2][3] - coeff*lhs[b1][b2][2][1];
  lhs[b1][b2][3][3] = lhs[b1][b2][3][3] - coeff*lhs[b1][b2][3][1];
  lhs[b1][b2][4][3] = lhs[b1][b2][4][3] - coeff*lhs[b1][b2][4][1];
  rhs[c1][c2][c3][3]   = rhs[c1][c2][c3][3]   - coeff*rhs[c1][c2][c3][1];

  coeff = lhs[a1][a2][1][4];
  lhs[a1][a2][2][4]= lhs[a1][a2][2][4] - coeff*lhs[a1][a2][2][1];
  lhs[a1][a2][3][4]= lhs[a1][a2][3][4] - coeff*lhs[a1][a2][3][1];
  lhs[a1][a2][4][4]= lhs[a1][a2][4][4] - coeff*lhs[a1][a2][4][1];
  lhs[b1][b2][0][4] = lhs[b1][b2][0][4] - coeff*lhs[b1][b2][0][1];
  lhs[b1][b2][1][4] = lhs[b1][b2][1][4] - coeff*lhs[b1][b2][1][1];
  lhs[b1][b2][2][4] = lhs[b1][b2][2][4] - coeff*lhs[b1][b2][2][1];
  lhs[b1][b2][3][4] = lhs[b1][b2][3][4] - coeff*lhs[b1][b2][3][1];
  lhs[b1][b2][4][4] = lhs[b1][b2][4][4] - coeff*lhs[b1][b2][4][1];
  rhs[c1][c2][c3][4]   = rhs[c1][c2][c3][4]   - coeff*rhs[c1][c2][c3][1];


  pivot = 1.00/lhs[a1][a2][2][2];
  lhs[a1][a2][3][2] = lhs[a1][a2][3][2]*pivot;
  lhs[a1][a2][4][2] = lhs[a1][a2][4][2]*pivot;
  lhs[b1][b2][0][2] = lhs[b1][b2][0][2]*pivot;
  lhs[b1][b2][1][2] = lhs[b1][b2][1][2]*pivot;
  lhs[b1][b2][2][2] = lhs[b1][b2][2][2]*pivot;
  lhs[b1][b2][3][2] = lhs[b1][b2][3][2]*pivot;
  lhs[b1][b2][4][2] = lhs[b1][b2][4][2]*pivot;
  rhs[c1][c2][c3][2]   = rhs[c1][c2][c3][2]  *pivot;

  coeff = lhs[a1][a2][2][0];
  lhs[a1][a2][3][0]= lhs[a1][a2][3][0] - coeff*lhs[a1][a2][3][2];
  lhs[a1][a2][4][0]= lhs[a1][a2][4][0] - coeff*lhs[a1][a2][4][2];
  lhs[b1][b2][0][0] = lhs[b1][b2][0][0] - coeff*lhs[b1][b2][0][2];
  lhs[b1][b2][1][0] = lhs[b1][b2][1][0] - coeff*lhs[b1][b2][1][2];
  lhs[b1][b2][2][0] = lhs[b1][b2][2][0] - coeff*lhs[b1][b2][2][2];
  lhs[b1][b2][3][0] = lhs[b1][b2][3][0] - coeff*lhs[b1][b2][3][2];
  lhs[b1][b2][4][0] = lhs[b1][b2][4][0] - coeff*lhs[b1][b2][4][2];
  rhs[c1][c2][c3][0]   = rhs[c1][c2][c3][0]   - coeff*rhs[c1][c2][c3][2];

  coeff = lhs[a1][a2][2][1];
  lhs[a1][a2][3][1]= lhs[a1][a2][3][1] - coeff*lhs[a1][a2][3][2];
  lhs[a1][a2][4][1]= lhs[a1][a2][4][1] - coeff*lhs[a1][a2][4][2];
  lhs[b1][b2][0][1] = lhs[b1][b2][0][1] - coeff*lhs[b1][b2][0][2];
  lhs[b1][b2][1][1] = lhs[b1][b2][1][1] - coeff*lhs[b1][b2][1][2];
  lhs[b1][b2][2][1] = lhs[b1][b2][2][1] - coeff*lhs[b1][b2][2][2];
  lhs[b1][b2][3][1] = lhs[b1][b2][3][1] - coeff*lhs[b1][b2][3][2];
  lhs[b1][b2][4][1] = lhs[b1][b2][4][1] - coeff*lhs[b1][b2][4][2];
  rhs[c1][c2][c3][1]   = rhs[c1][c2][c3][1]   - coeff*rhs[c1][c2][c3][2];

  coeff = lhs[a1][a2][2][3];
  lhs[a1][a2][3][3]= lhs[a1][a2][3][3] - coeff*lhs[a1][a2][3][2];
  lhs[a1][a2][4][3]= lhs[a1][a2][4][3] - coeff*lhs[a1][a2][4][2];
  lhs[b1][b2][0][3] = lhs[b1][b2][0][3] - coeff*lhs[b1][b2][0][2];
  lhs[b1][b2][1][3] = lhs[b1][b2][1][3] - coeff*lhs[b1][b2][1][2];
  lhs[b1][b2][2][3] = lhs[b1][b2][2][3] - coeff*lhs[b1][b2][2][2];
  lhs[b1][b2][3][3] = lhs[b1][b2][3][3] - coeff*lhs[b1][b2][3][2];
  lhs[b1][b2][4][3] = lhs[b1][b2][4][3] - coeff*lhs[b1][b2][4][2];
  rhs[c1][c2][c3][3]   = rhs[c1][c2][c3][3]   - coeff*rhs[c1][c2][c3][2];

  coeff = lhs[a1][a2][2][4];
  lhs[a1][a2][3][4]= lhs[a1][a2][3][4] - coeff*lhs[a1][a2][3][2];
  lhs[a1][a2][4][4]= lhs[a1][a2][4][4] - coeff*lhs[a1][a2][4][2];
  lhs[b1][b2][0][4] = lhs[b1][b2][0][4] - coeff*lhs[b1][b2][0][2];
  lhs[b1][b2][1][4] = lhs[b1][b2][1][4] - coeff*lhs[b1][b2][1][2];
  lhs[b1][b2][2][4] = lhs[b1][b2][2][4] - coeff*lhs[b1][b2][2][2];
  lhs[b1][b2][3][4] = lhs[b1][b2][3][4] - coeff*lhs[b1][b2][3][2];
  lhs[b1][b2][4][4] = lhs[b1][b2][4][4] - coeff*lhs[b1][b2][4][2];
  rhs[c1][c2][c3][4]   = rhs[c1][c2][c3][4]   - coeff*rhs[c1][c2][c3][2];


  pivot = 1.00/lhs[a1][a2][3][3];
  lhs[a1][a2][4][3] = lhs[a1][a2][4][3]*pivot;
  lhs[b1][b2][0][3] = lhs[b1][b2][0][3]*pivot;
  lhs[b1][b2][1][3] = lhs[b1][b2][1][3]*pivot;
  lhs[b1][b2][2][3] = lhs[b1][b2][2][3]*pivot;
  lhs[b1][b2][3][3] = lhs[b1][b2][3][3]*pivot;
  lhs[b1][b2][4][3] = lhs[b1][b2][4][3]*pivot;
  rhs[c1][c2][c3][3]   = rhs[c1][c2][c3][3]  *pivot;

  coeff = lhs[a1][a2][3][0];
  lhs[a1][a2][4][0]= lhs[a1][a2][4][0] - coeff*lhs[a1][a2][4][3];
  lhs[b1][b2][0][0] = lhs[b1][b2][0][0] - coeff*lhs[b1][b2][0][3];
  lhs[b1][b2][1][0] = lhs[b1][b2][1][0] - coeff*lhs[b1][b2][1][3];
  lhs[b1][b2][2][0] = lhs[b1][b2][2][0] - coeff*lhs[b1][b2][2][3];
  lhs[b1][b2][3][0] = lhs[b1][b2][3][0] - coeff*lhs[b1][b2][3][3];
  lhs[b1][b2][4][0] = lhs[b1][b2][4][0] - coeff*lhs[b1][b2][4][3];
  rhs[c1][c2][c3][0]   = rhs[c1][c2][c3][0]   - coeff*rhs[c1][c2][c3][3];

  coeff = lhs[a1][a2][3][1];
  lhs[a1][a2][4][1]= lhs[a1][a2][4][1] - coeff*lhs[a1][a2][4][3];
  lhs[b1][b2][0][1] = lhs[b1][b2][0][1] - coeff*lhs[b1][b2][0][3];
  lhs[b1][b2][1][1] = lhs[b1][b2][1][1] - coeff*lhs[b1][b2][1][3];
  lhs[b1][b2][2][1] = lhs[b1][b2][2][1] - coeff*lhs[b1][b2][2][3];
  lhs[b1][b2][3][1] = lhs[b1][b2][3][1] - coeff*lhs[b1][b2][3][3];
  lhs[b1][b2][4][1] = lhs[b1][b2][4][1] - coeff*lhs[b1][b2][4][3];
  rhs[c1][c2][c3][1]   = rhs[c1][c2][c3][1]   - coeff*rhs[c1][c2][c3][3];

  coeff = lhs[a1][a2][3][2];
  lhs[a1][a2][4][2]= lhs[a1][a2][4][2] - coeff*lhs[a1][a2][4][3];
  lhs[b1][b2][0][2] = lhs[b1][b2][0][2] - coeff*lhs[b1][b2][0][3];
  lhs[b1][b2][1][2] = lhs[b1][b2][1][2] - coeff*lhs[b1][b2][1][3];
  lhs[b1][b2][2][2] = lhs[b1][b2][2][2] - coeff*lhs[b1][b2][2][3];
  lhs[b1][b2][3][2] = lhs[b1][b2][3][2] - coeff*lhs[b1][b2][3][3];
  lhs[b1][b2][4][2] = lhs[b1][b2][4][2] - coeff*lhs[b1][b2][4][3];
  rhs[c1][c2][c3][2]   = rhs[c1][c2][c3][2]   - coeff*rhs[c1][c2][c3][3];

  coeff = lhs[a1][a2][3][4];
  lhs[a1][a2][4][4]= lhs[a1][a2][4][4] - coeff*lhs[a1][a2][4][3];
  lhs[b1][b2][0][4] = lhs[b1][b2][0][4] - coeff*lhs[b1][b2][0][3];
  lhs[b1][b2][1][4] = lhs[b1][b2][1][4] - coeff*lhs[b1][b2][1][3];
  lhs[b1][b2][2][4] = lhs[b1][b2][2][4] - coeff*lhs[b1][b2][2][3];
  lhs[b1][b2][3][4] = lhs[b1][b2][3][4] - coeff*lhs[b1][b2][3][3];
  lhs[b1][b2][4][4] = lhs[b1][b2][4][4] - coeff*lhs[b1][b2][4][3];
  rhs[c1][c2][c3][4]   = rhs[c1][c2][c3][4]   - coeff*rhs[c1][c2][c3][3];


  pivot = 1.00/lhs[a1][a2][4][4];
  lhs[b1][b2][0][4] = lhs[b1][b2][0][4]*pivot;
  lhs[b1][b2][1][4] = lhs[b1][b2][1][4]*pivot;
  lhs[b1][b2][2][4] = lhs[b1][b2][2][4]*pivot;
  lhs[b1][b2][3][4] = lhs[b1][b2][3][4]*pivot;
  lhs[b1][b2][4][4] = lhs[b1][b2][4][4]*pivot;
  rhs[c1][c2][c3][4]   = rhs[c1][c2][c3][4]  *pivot;

  coeff = lhs[a1][a2][4][0];
  lhs[b1][b2][0][0] = lhs[b1][b2][0][0] - coeff*lhs[b1][b2][0][4];
  lhs[b1][b2][1][0] = lhs[b1][b2][1][0] - coeff*lhs[b1][b2][1][4];
  lhs[b1][b2][2][0] = lhs[b1][b2][2][0] - coeff*lhs[b1][b2][2][4];
  lhs[b1][b2][3][0] = lhs[b1][b2][3][0] - coeff*lhs[b1][b2][3][4];
  lhs[b1][b2][4][0] = lhs[b1][b2][4][0] - coeff*lhs[b1][b2][4][4];
  rhs[c1][c2][c3][0]   = rhs[c1][c2][c3][0]   - coeff*rhs[c1][c2][c3][4];

  coeff = lhs[a1][a2][4][1];
  lhs[b1][b2][0][1] = lhs[b1][b2][0][1] - coeff*lhs[b1][b2][0][4];
  lhs[b1][b2][1][1] = lhs[b1][b2][1][1] - coeff*lhs[b1][b2][1][4];
  lhs[b1][b2][2][1] = lhs[b1][b2][2][1] - coeff*lhs[b1][b2][2][4];
  lhs[b1][b2][3][1] = lhs[b1][b2][3][1] - coeff*lhs[b1][b2][3][4];
  lhs[b1][b2][4][1] = lhs[b1][b2][4][1] - coeff*lhs[b1][b2][4][4];
  rhs[c1][c2][c3][1]   = rhs[c1][c2][c3][1]   - coeff*rhs[c1][c2][c3][4];

  coeff = lhs[a1][a2][4][2];
  lhs[b1][b2][0][2] = lhs[b1][b2][0][2] - coeff*lhs[b1][b2][0][4];
  lhs[b1][b2][1][2] = lhs[b1][b2][1][2] - coeff*lhs[b1][b2][1][4];
  lhs[b1][b2][2][2] = lhs[b1][b2][2][2] - coeff*lhs[b1][b2][2][4];
  lhs[b1][b2][3][2] = lhs[b1][b2][3][2] - coeff*lhs[b1][b2][3][4];
  lhs[b1][b2][4][2] = lhs[b1][b2][4][2] - coeff*lhs[b1][b2][4][4];
  rhs[c1][c2][c3][2]   = rhs[c1][c2][c3][2]   - coeff*rhs[c1][c2][c3][4];

  coeff = lhs[a1][a2][4][3];
  lhs[b1][b2][0][3] = lhs[b1][b2][0][3] - coeff*lhs[b1][b2][0][4];
  lhs[b1][b2][1][3] = lhs[b1][b2][1][3] - coeff*lhs[b1][b2][1][4];
  lhs[b1][b2][2][3] = lhs[b1][b2][2][3] - coeff*lhs[b1][b2][2][4];
  lhs[b1][b2][3][3] = lhs[b1][b2][3][3] - coeff*lhs[b1][b2][3][4];
  lhs[b1][b2][4][3] = lhs[b1][b2][4][3] - coeff*lhs[b1][b2][4][4];
  rhs[c1][c2][c3][3]   = rhs[c1][c2][c3][3]   - coeff*rhs[c1][c2][c3][4];
}

////////////////////////////////////


/*void binvrhs(double lhs[5][5], double r[5])
{
  double pivot, coeff;

  pivot = 1.00/lhs[0][0];
  lhs[1][0] = lhs[1][0]*pivot;
  lhs[2][0] = lhs[2][0]*pivot;
  lhs[3][0] = lhs[3][0]*pivot;
  lhs[4][0] = lhs[4][0]*pivot;
  r[0]   = r[0]  *pivot;

  coeff = lhs[0][1];
  lhs[1][1]= lhs[1][1] - coeff*lhs[1][0];
  lhs[2][1]= lhs[2][1] - coeff*lhs[2][0];
  lhs[3][1]= lhs[3][1] - coeff*lhs[3][0];
  lhs[4][1]= lhs[4][1] - coeff*lhs[4][0];
  r[1]   = r[1]   - coeff*r[0];

  coeff = lhs[0][2];
  lhs[1][2]= lhs[1][2] - coeff*lhs[1][0];
  lhs[2][2]= lhs[2][2] - coeff*lhs[2][0];
  lhs[3][2]= lhs[3][2] - coeff*lhs[3][0];
  lhs[4][2]= lhs[4][2] - coeff*lhs[4][0];
  r[2]   = r[2]   - coeff*r[0];

  coeff = lhs[0][3];
  lhs[1][3]= lhs[1][3] - coeff*lhs[1][0];
  lhs[2][3]= lhs[2][3] - coeff*lhs[2][0];
  lhs[3][3]= lhs[3][3] - coeff*lhs[3][0];
  lhs[4][3]= lhs[4][3] - coeff*lhs[4][0];
  r[3]   = r[3]   - coeff*r[0];

  coeff = lhs[0][4];
  lhs[1][4]= lhs[1][4] - coeff*lhs[1][0];
  lhs[2][4]= lhs[2][4] - coeff*lhs[2][0];
  lhs[3][4]= lhs[3][4] - coeff*lhs[3][0];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][0];
  r[4]   = r[4]   - coeff*r[0];


  pivot = 1.00/lhs[1][1];
  lhs[2][1] = lhs[2][1]*pivot;
  lhs[3][1] = lhs[3][1]*pivot;
  lhs[4][1] = lhs[4][1]*pivot;
  r[1]   = r[1]  *pivot;

  coeff = lhs[1][0];
  lhs[2][0]= lhs[2][0] - coeff*lhs[2][1];
  lhs[3][0]= lhs[3][0] - coeff*lhs[3][1];
  lhs[4][0]= lhs[4][0] - coeff*lhs[4][1];
  r[0]   = r[0]   - coeff*r[1];

  coeff = lhs[1][2];
  lhs[2][2]= lhs[2][2] - coeff*lhs[2][1];
  lhs[3][2]= lhs[3][2] - coeff*lhs[3][1];
  lhs[4][2]= lhs[4][2] - coeff*lhs[4][1];
  r[2]   = r[2]   - coeff*r[1];

  coeff = lhs[1][3];
  lhs[2][3]= lhs[2][3] - coeff*lhs[2][1];
  lhs[3][3]= lhs[3][3] - coeff*lhs[3][1];
  lhs[4][3]= lhs[4][3] - coeff*lhs[4][1];
  r[3]   = r[3]   - coeff*r[1];

  coeff = lhs[1][4];
  lhs[2][4]= lhs[2][4] - coeff*lhs[2][1];
  lhs[3][4]= lhs[3][4] - coeff*lhs[3][1];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][1];
  r[4]   = r[4]   - coeff*r[1];


  pivot = 1.00/lhs[2][2];
  lhs[3][2] = lhs[3][2]*pivot;
  lhs[4][2] = lhs[4][2]*pivot;
  r[2]   = r[2]  *pivot;

  coeff = lhs[2][0];
  lhs[3][0]= lhs[3][0] - coeff*lhs[3][2];
  lhs[4][0]= lhs[4][0] - coeff*lhs[4][2];
  r[0]   = r[0]   - coeff*r[2];

  coeff = lhs[2][1];
  lhs[3][1]= lhs[3][1] - coeff*lhs[3][2];
  lhs[4][1]= lhs[4][1] - coeff*lhs[4][2];
  r[1]   = r[1]   - coeff*r[2];

  coeff = lhs[2][3];
  lhs[3][3]= lhs[3][3] - coeff*lhs[3][2];
  lhs[4][3]= lhs[4][3] - coeff*lhs[4][2];
  r[3]   = r[3]   - coeff*r[2];

  coeff = lhs[2][4];
  lhs[3][4]= lhs[3][4] - coeff*lhs[3][2];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][2];
  r[4]   = r[4]   - coeff*r[2];


  pivot = 1.00/lhs[3][3];
  lhs[4][3] = lhs[4][3]*pivot;
  r[3]   = r[3]  *pivot;

  coeff = lhs[3][0];
  lhs[4][0]= lhs[4][0] - coeff*lhs[4][3];
  r[0]   = r[0]   - coeff*r[3];

  coeff = lhs[3][1];
  lhs[4][1]= lhs[4][1] - coeff*lhs[4][3];
  r[1]   = r[1]   - coeff*r[3];

  coeff = lhs[3][2];
  lhs[4][2]= lhs[4][2] - coeff*lhs[4][3];
  r[2]   = r[2]   - coeff*r[3];

  coeff = lhs[3][4];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][3];
  r[4]   = r[4]   - coeff*r[3];


  pivot = 1.00/lhs[4][4];
  r[4]   = r[4]  *pivot;

  coeff = lhs[4][0];
  r[0]   = r[0]   - coeff*r[4];

  coeff = lhs[4][1];
  r[1]   = r[1]   - coeff*r[4];

  coeff = lhs[4][2];
  r[2]   = r[2]   - coeff*r[4];

  coeff = lhs[4][3];
  r[3]   = r[3]   - coeff*r[4];
}
*/

/////////////////////////////////////////

void binvrhs(int a1, int a2, int b1, int b2, int b3)
{
  double pivot, coeff;

  pivot = 1.00/lhs[a1][a2][0][0];
  lhs[a1][a2][1][0] = lhs[a1][a2][1][0]*pivot;
  lhs[a1][a2][2][0] = lhs[a1][a2][2][0]*pivot;
  lhs[a1][a2][3][0] = lhs[a1][a2][3][0]*pivot;
  lhs[a1][a2][4][0] = lhs[a1][a2][4][0]*pivot;
  rhs[b1][b2][b3][0]   = rhs[b1][b2][b3][0]  *pivot;

  coeff = lhs[a1][a2][0][1];
  lhs[a1][a2][1][1]= lhs[a1][a2][1][1] - coeff*lhs[a1][a2][1][0];
  lhs[a1][a2][2][1]= lhs[a1][a2][2][1] - coeff*lhs[a1][a2][2][0];
  lhs[a1][a2][3][1]= lhs[a1][a2][3][1] - coeff*lhs[a1][a2][3][0];
  lhs[a1][a2][4][1]= lhs[a1][a2][4][1] - coeff*lhs[a1][a2][4][0];
  rhs[b1][b2][b3][1]   = rhs[b1][b2][b3][1]   - coeff*rhs[b1][b2][b3][0];

  coeff = lhs[a1][a2][0][2];
  lhs[a1][a2][1][2]= lhs[a1][a2][1][2] - coeff*lhs[a1][a2][1][0];
  lhs[a1][a2][2][2]= lhs[a1][a2][2][2] - coeff*lhs[a1][a2][2][0];
  lhs[a1][a2][3][2]= lhs[a1][a2][3][2] - coeff*lhs[a1][a2][3][0];
  lhs[a1][a2][4][2]= lhs[a1][a2][4][2] - coeff*lhs[a1][a2][4][0];
  rhs[b1][b2][b3][2]   = rhs[b1][b2][b3][2]   - coeff*rhs[b1][b2][b3][0];

  coeff = lhs[a1][a2][0][3];
  lhs[a1][a2][1][3]= lhs[a1][a2][1][3] - coeff*lhs[a1][a2][1][0];
  lhs[a1][a2][2][3]= lhs[a1][a2][2][3] - coeff*lhs[a1][a2][2][0];
  lhs[a1][a2][3][3]= lhs[a1][a2][3][3] - coeff*lhs[a1][a2][3][0];
  lhs[a1][a2][4][3]= lhs[a1][a2][4][3] - coeff*lhs[a1][a2][4][0];
  rhs[b1][b2][b3][3]   = rhs[b1][b2][b3][3]   - coeff*rhs[b1][b2][b3][0];

  coeff = lhs[a1][a2][0][4];
  lhs[a1][a2][1][4]= lhs[a1][a2][1][4] - coeff*lhs[a1][a2][1][0];
  lhs[a1][a2][2][4]= lhs[a1][a2][2][4] - coeff*lhs[a1][a2][2][0];
  lhs[a1][a2][3][4]= lhs[a1][a2][3][4] - coeff*lhs[a1][a2][3][0];
  lhs[a1][a2][4][4]= lhs[a1][a2][4][4] - coeff*lhs[a1][a2][4][0];
  rhs[b1][b2][b3][4]   = rhs[b1][b2][b3][4]   - coeff*rhs[b1][b2][b3][0];


  pivot = 1.00/lhs[a1][a2][1][1];
  lhs[a1][a2][2][1] = lhs[a1][a2][2][1]*pivot;
  lhs[a1][a2][3][1] = lhs[a1][a2][3][1]*pivot;
  lhs[a1][a2][4][1] = lhs[a1][a2][4][1]*pivot;
  rhs[b1][b2][b3][1]   = rhs[b1][b2][b3][1]  *pivot;

  coeff = lhs[a1][a2][1][0];
  lhs[a1][a2][2][0]= lhs[a1][a2][2][0] - coeff*lhs[a1][a2][2][1];
  lhs[a1][a2][3][0]= lhs[a1][a2][3][0] - coeff*lhs[a1][a2][3][1];
  lhs[a1][a2][4][0]= lhs[a1][a2][4][0] - coeff*lhs[a1][a2][4][1];
  rhs[b1][b2][b3][0]   = rhs[b1][b2][b3][0]   - coeff*rhs[b1][b2][b3][1];

  coeff = lhs[a1][a2][1][2];
  lhs[a1][a2][2][2]= lhs[a1][a2][2][2] - coeff*lhs[a1][a2][2][1];
  lhs[a1][a2][3][2]= lhs[a1][a2][3][2] - coeff*lhs[a1][a2][3][1];
  lhs[a1][a2][4][2]= lhs[a1][a2][4][2] - coeff*lhs[a1][a2][4][1];
  rhs[b1][b2][b3][2]   = rhs[b1][b2][b3][2]   - coeff*rhs[b1][b2][b3][1];

  coeff = lhs[a1][a2][1][3];
  lhs[a1][a2][2][3]= lhs[a1][a2][2][3] - coeff*lhs[a1][a2][2][1];
  lhs[a1][a2][3][3]= lhs[a1][a2][3][3] - coeff*lhs[a1][a2][3][1];
  lhs[a1][a2][4][3]= lhs[a1][a2][4][3] - coeff*lhs[a1][a2][4][1];
  rhs[b1][b2][b3][3]   = rhs[b1][b2][b3][3]   - coeff*rhs[b1][b2][b3][1];

  coeff = lhs[a1][a2][1][4];
  lhs[a1][a2][2][4]= lhs[a1][a2][2][4] - coeff*lhs[a1][a2][2][1];
  lhs[a1][a2][3][4]= lhs[a1][a2][3][4] - coeff*lhs[a1][a2][3][1];
  lhs[a1][a2][4][4]= lhs[a1][a2][4][4] - coeff*lhs[a1][a2][4][1];
  rhs[b1][b2][b3][4]   = rhs[b1][b2][b3][4]   - coeff*rhs[b1][b2][b3][1];


  pivot = 1.00/lhs[a1][a2][2][2];
  lhs[a1][a2][3][2] = lhs[a1][a2][3][2]*pivot;
  lhs[a1][a2][4][2] = lhs[a1][a2][4][2]*pivot;
  rhs[b1][b2][b3][2]   = rhs[b1][b2][b3][2]  *pivot;

  coeff = lhs[a1][a2][2][0];
  lhs[a1][a2][3][0]= lhs[a1][a2][3][0] - coeff*lhs[a1][a2][3][2];
  lhs[a1][a2][4][0]= lhs[a1][a2][4][0] - coeff*lhs[a1][a2][4][2];
  rhs[b1][b2][b3][0]   = rhs[b1][b2][b3][0]   - coeff*rhs[b1][b2][b3][2];

  coeff = lhs[a1][a2][2][1];
  lhs[a1][a2][3][1]= lhs[a1][a2][3][1] - coeff*lhs[a1][a2][3][2];
  lhs[a1][a2][4][1]= lhs[a1][a2][4][1] - coeff*lhs[a1][a2][4][2];
  rhs[b1][b2][b3][1]   = rhs[b1][b2][b3][1]   - coeff*rhs[b1][b2][b3][2];

  coeff = lhs[a1][a2][2][3];
  lhs[a1][a2][3][3]= lhs[a1][a2][3][3] - coeff*lhs[a1][a2][3][2];
  lhs[a1][a2][4][3]= lhs[a1][a2][4][3] - coeff*lhs[a1][a2][4][2];
  rhs[b1][b2][b3][3]   = rhs[b1][b2][b3][3]   - coeff*rhs[b1][b2][b3][2];

  coeff = lhs[a1][a2][2][4];
  lhs[a1][a2][3][4]= lhs[a1][a2][3][4] - coeff*lhs[a1][a2][3][2];
  lhs[a1][a2][4][4]= lhs[a1][a2][4][4] - coeff*lhs[a1][a2][4][2];
  rhs[b1][b2][b3][4]   = rhs[b1][b2][b3][4]   - coeff*rhs[b1][b2][b3][2];


  pivot = 1.00/lhs[a1][a2][3][3];
  lhs[a1][a2][4][3] = lhs[a1][a2][4][3]*pivot;
  rhs[b1][b2][b3][3]   = rhs[b1][b2][b3][3]  *pivot;

  coeff = lhs[a1][a2][3][0];
  lhs[a1][a2][4][0]= lhs[a1][a2][4][0] - coeff*lhs[a1][a2][4][3];
  rhs[b1][b2][b3][0]   = rhs[b1][b2][b3][0]   - coeff*rhs[b1][b2][b3][3];

  coeff = lhs[a1][a2][3][1];
  lhs[a1][a2][4][1]= lhs[a1][a2][4][1] - coeff*lhs[a1][a2][4][3];
  rhs[b1][b2][b3][1]   = rhs[b1][b2][b3][1]   - coeff*rhs[b1][b2][b3][3];

  coeff = lhs[a1][a2][3][2];
  lhs[a1][a2][4][2]= lhs[a1][a2][4][2] - coeff*lhs[a1][a2][4][3];
  rhs[b1][b2][b3][2]   = rhs[b1][b2][b3][2]   - coeff*rhs[b1][b2][b3][3];

  coeff = lhs[a1][a2][3][4];
  lhs[a1][a2][4][4]= lhs[a1][a2][4][4] - coeff*lhs[a1][a2][4][3];
  rhs[b1][b2][b3][4]   = rhs[b1][b2][b3][4]   - coeff*rhs[b1][b2][b3][3];


  pivot = 1.00/lhs[a1][a2][4][4];
  rhs[b1][b2][b3][4]   = rhs[b1][b2][b3][4]  *pivot;

  coeff = lhs[a1][a2][4][0];
  rhs[b1][b2][b3][0]   = rhs[b1][b2][b3][0]   - coeff*rhs[b1][b2][b3][4];

  coeff = lhs[a1][a2][4][1];
  rhs[b1][b2][b3][1]   = rhs[b1][b2][b3][1]   - coeff*rhs[b1][b2][b3][4];

  coeff = lhs[a1][a2][4][2];
  rhs[b1][b2][b3][2]   = rhs[b1][b2][b3][2]   - coeff*rhs[b1][b2][b3][4];

  coeff = lhs[a1][a2][4][3];
  rhs[b1][b2][b3][3]   = rhs[b1][b2][b3][3]   - coeff*rhs[b1][b2][b3][4];
}

/////////////////////////////////////////


