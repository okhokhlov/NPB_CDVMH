#include <cdvmh_helpers.h>

//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is a serial C version of the NPB LU code. This C        //
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
//---  applu.incl   
//---------------------------------------------------------------------
//---------------------------------------------------------------------

//---------------------------------------------------------------------
// npbparams.h defines parameters that depend on the class and 
// number of nodes
//---------------------------------------------------------------------

#include "npbparams.h"
#include "../common/type.h"

//---------------------------------------------------------------------
// parameters which can be overridden in runtime config file
// isiz1,isiz2,isiz3 give the maximum size
// ipr = 1 to print out verbose information
// omega = 2.0 is correct for all classes
// tolrsd is tolerance levels for steady state residuals
//---------------------------------------------------------------------
#define IPR_DEFAULT     1
#define OMEGA_DEFAULT   1.2
#define TOLRSD1_DEF     1.0e-08
#define TOLRSD2_DEF     1.0e-08
#define TOLRSD3_DEF     1.0e-08
#define TOLRSD4_DEF     1.0e-08
#define TOLRSD5_DEF     1.0e-08

#define C1              1.40e+00
#define C2              0.40e+00
#define C3              1.00e-01
#define C4              1.00e+00
#define C5              1.40e+00

//---------------------------------------------------------------------
// grid
//---------------------------------------------------------------------
/* common/cgcon/ */
extern double dxi, deta, dzeta;
extern double tx1, tx2, tx3;
extern double ty1, ty2, ty3;
extern double tz1, tz2, tz3;
extern int nx, ny, nz;
extern int nx0, ny0, nz0;
extern int ist, iend;
extern int jst, jend;
extern int ii1, ii2;
extern int ji1, ji2;
extern int ki1, ki2;

//---------------------------------------------------------------------
// dissipation
//---------------------------------------------------------------------
/* common/disp/ */
extern double dx1, dx2, dx3, dx4, dx5;
extern double dy1, dy2, dy3, dy4, dy5;
extern double dz1, dz2, dz3, dz4, dz5;
extern double dssp;

//---------------------------------------------------------------------
// field variables and residuals
// to improve cache performance, second two dimensions padded by 1 
// for even number sizes only.
// Note: corresponding array (called "v") in routines blts, buts, 
// and l2norm are similarly padded
//---------------------------------------------------------------------
/* common/cvar/ */
extern double u    [ISIZ3][ISIZ2/2*2+1][ISIZ1/2*2+1][5];
extern double rsd  [ISIZ3][ISIZ2/2*2+1][ISIZ1/2*2+1][5];
extern double frct [ISIZ3][ISIZ2/2*2+1][ISIZ1/2*2+1][5];
extern double flux [ISIZ1][5];
extern double qs   [ISIZ3][ISIZ2/2*2+1][ISIZ1/2*2+1];
extern double rho_i[ISIZ3][ISIZ2/2*2+1][ISIZ1/2*2+1];

//---------------------------------------------------------------------
// output control parameters
//---------------------------------------------------------------------
/* common/cprcon/ */
extern int ipr, inorm;

//---------------------------------------------------------------------
// newton-raphson iteration control parameters
//---------------------------------------------------------------------
/* common/ctscon/ */
extern double dt, omega, tolrsd[5], rsdnm[5], errnm[5], frc, ttotal;
extern int itmax, invert;

/* common/cjac/ */
extern double a[ISIZ2][ISIZ1/2*2+1][5][5];
extern double b[ISIZ2][ISIZ1/2*2+1][5][5];
extern double c[ISIZ2][ISIZ1/2*2+1][5][5];
extern double d[ISIZ2][ISIZ1/2*2+1][5][5];


//---------------------------------------------------------------------
// coefficients of the exact solution
//---------------------------------------------------------------------
/* common/cexact/ */
extern double ce[5][13];


//---------------------------------------------------------------------
// timers
//---------------------------------------------------------------------
/* common/timer/ */
extern double maxtime;
extern logical timeron;
#define t_total   1
#define t_rhsx    2
#define t_rhsy    3
#define t_rhsz    4
#define t_rhs     5
#define t_jacld   6
#define t_blts    7
#define t_jacu    8
#define t_buts    9
#define t_add     10
#define t_l2norm  11
#define t_last    11


void read_input();
void domain();
void setcoeff();
void setbv();
void exact(int i, int j, int k, double u000ijk[]);
void setiv();
void erhs();
void ssor(int niter);
void rhs();
void l2norm (int ldx, int ldy, int ldz, int nx0, int ny0, int nz0,
     int ist, int iend, int jst, int jend,
     double v[][ldy/2*2+1][ldx/2*2+1][5], double sum[5]);
void jacld(int k);
void blts (int ldmx, int ldmy, int ldmz, int nx, int ny, int nz, int k,
    double omega,
    double v[][ldmy/2*2+1][ldmx/2*2+1][5], 
    double ldz[ldmy][ldmx/2*2+1][5][5],
    double ldy[ldmy][ldmx/2*2+1][5][5],
    double ldx[ldmy][ldmx/2*2+1][5][5],
    double d[ldmy][ldmx/2*2+1][5][5],
    int ist, int iend, int jst, int jend, int nx0, int ny0);
void jacu(int k);
void buts(int ldmx, int ldmy, int ldmz, int nx, int ny, int nz, int k,
    double omega,
    double v[][ldmy/2*2+1][ldmx/2*2+1][5],
    double tv[ldmy][ldmx/2*2+1][5],
    double d[ldmy][ldmx/2*2+1][5][5],
    double udx[ldmy][ldmx/2*2+1][5][5],
    double udy[ldmy][ldmx/2*2+1][5][5],
    double udz[ldmy][ldmx/2*2+1][5][5],
    int ist, int iend, int jst, int jend, int nx0, int ny0);
void error();
void pintgr();
void verify(double xcr[5], double xce[5], double xci, 
            char *Class, logical *verified);


//---------------------------------------------------------------------
//   end of include file
//---------------------------------------------------------------------


//---------------------------------------------------------------------
// set the boundary values of dependent variables
//---------------------------------------------------------------------
void setbv()
{
  //---------------------------------------------------------------------
  // local variables
  //---------------------------------------------------------------------
  int i, j, k, m;
  double temp1[5], temp2[5];

  //---------------------------------------------------------------------
  // set the dependent variable values along the top and bottom faces
  //---------------------------------------------------------------------
  for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++) {
      exact( i, j, 0, temp1 );
      exact( i, j, nz-1, temp2 );
      for (m = 0; m < 5; m++) {
        u[0][j][i][m] = temp1[m];
        u[nz-1][j][i][m] = temp2[m];
      }
    }
  }

  //---------------------------------------------------------------------
  // set the dependent variable values along north and south faces
  //---------------------------------------------------------------------
  for (k = 0; k < nz; k++) {
    for (i = 0; i < nx; i++) {
      exact( i, 0, k, temp1 );
      exact( i, ny-1, k, temp2 );
      for (m = 0; m < 5; m++) {
        u[k][0][i][m] = temp1[m];
        u[k][ny-1][i][m] = temp2[m];
      }
    }
  }

  //---------------------------------------------------------------------
  // set the dependent variable values along east and west faces
  //---------------------------------------------------------------------
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      exact( 0, j, k, temp1 );
      exact( nx-1, j, k, temp2 );
      for (m = 0; m < 5; m++) {
        u[k][j][0][m] = temp1[m];
        u[k][j][nx-1][m] = temp2[m];
      }
    }
  }
}


