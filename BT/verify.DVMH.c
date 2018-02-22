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

#include <stdio.h>
#include <math.h>
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



//---------------------------------------------------------------------
// verification routine                         
//---------------------------------------------------------------------
void verify(int no_time_steps, char *Class, logical *verified)
{
  double xcrref[5], xceref[5], xcrdif[5], xcedif[5]; 
  double epsilon, xce[5], xcr[5], dtref = 0.0;
  int m;

  //---------------------------------------------------------------------
  // tolerance level
  //---------------------------------------------------------------------
  epsilon = 1.0e-08;

  //---------------------------------------------------------------------
  // compute the error norm and the residual norm, and exit if not printing
  //---------------------------------------------------------------------
  error_norm(xce);
  compute_rhs();

  rhs_norm(xcr);

  for (m = 0; m < 5; m++) {
    xcr[m] = xcr[m] / dt;
  }

  *Class = 'U';
  *verified = true;

  for (m = 0; m < 5; m++) {
    xcrref[m] = 1.0;
    xceref[m] = 1.0;
  }

  //---------------------------------------------------------------------
  // reference data for 12X12X12 grids after 60 time steps, with DT = 1.0e-02
  //---------------------------------------------------------------------
  if ( (grid_points[0] == 12) && (grid_points[1] == 12) &&
       (grid_points[2] == 12) && (no_time_steps == 60))  {

    *Class = 'S';
    dtref = 1.0e-2;

    //---------------------------------------------------------------------
    // Reference values of RMS-norms of residual.
    //---------------------------------------------------------------------
    xcrref[0] = 1.7034283709541311e-01;
    xcrref[1] = 1.2975252070034097e-02;
    xcrref[2] = 3.2527926989486055e-02;
    xcrref[3] = 2.6436421275166801e-02;
    xcrref[4] = 1.9211784131744430e-01;

    //---------------------------------------------------------------------
    // Reference values of RMS-norms of solution error.
    //---------------------------------------------------------------------
    xceref[0] = 4.9976913345811579e-04;
    xceref[1] = 4.5195666782961927e-05;
    xceref[2] = 7.3973765172921357e-05;
    xceref[3] = 7.3821238632439731e-05;
    xceref[4] = 8.9269630987491446e-04;

    //---------------------------------------------------------------------
    // reference data for 24X24X24 grids after 200 time steps, 
    // with DT = 0.8e-3
    //---------------------------------------------------------------------
  } else if ( (grid_points[0] == 24) && (grid_points[1] == 24) &&
              (grid_points[2] == 24) && (no_time_steps == 200) ) {

    *Class = 'W';
    dtref = 0.8e-3;
    //---------------------------------------------------------------------
    // Reference values of RMS-norms of residual.
    //---------------------------------------------------------------------
    xcrref[0] = 0.1125590409344e+03;
    xcrref[1] = 0.1180007595731e+02;
    xcrref[2] = 0.2710329767846e+02;
    xcrref[3] = 0.2469174937669e+02;
    xcrref[4] = 0.2638427874317e+03;

    //---------------------------------------------------------------------
    // Reference values of RMS-norms of solution error.
    //---------------------------------------------------------------------
    xceref[0] = 0.4419655736008e+01;
    xceref[1] = 0.4638531260002e+00;
    xceref[2] = 0.1011551749967e+01;
    xceref[3] = 0.9235878729944e+00;
    xceref[4] = 0.1018045837718e+02;

    //---------------------------------------------------------------------
    // reference data for 64X64X64 grids after 200 time steps, 
    // with DT = 0.8e-3
    //---------------------------------------------------------------------
  } else if ( (grid_points[0] == 64) && (grid_points[1] == 64) &&
              (grid_points[2] == 64) && (no_time_steps == 200) ) {

    *Class = 'A';
    dtref = 0.8e-3;
    //---------------------------------------------------------------------
    // Reference values of RMS-norms of residual.
    //---------------------------------------------------------------------
    xcrref[0] = 1.0806346714637264e+02;
    xcrref[1] = 1.1319730901220813e+01;
    xcrref[2] = 2.5974354511582465e+01;
    xcrref[3] = 2.3665622544678910e+01;
    xcrref[4] = 2.5278963211748344e+02;

    //---------------------------------------------------------------------
    // Reference values of RMS-norms of solution error.
    //---------------------------------------------------------------------
    xceref[0] = 4.2348416040525025e+00;
    xceref[1] = 4.4390282496995698e-01;
    xceref[2] = 9.6692480136345650e-01;
    xceref[3] = 8.8302063039765474e-01;
    xceref[4] = 9.7379901770829278e+00;

    //---------------------------------------------------------------------
    // reference data for 102X102X102 grids after 200 time steps,
    // with DT = 3.0e-04
    //---------------------------------------------------------------------
  } else if ( (grid_points[0] == 102) && (grid_points[1] == 102) &&
              (grid_points[2] == 102) && (no_time_steps == 200) ) {

    *Class = 'B';
    dtref = 3.0e-4;

    //---------------------------------------------------------------------
    // Reference values of RMS-norms of residual.
    //---------------------------------------------------------------------
    xcrref[0] = 1.4233597229287254e+03;
    xcrref[1] = 9.9330522590150238e+01;
    xcrref[2] = 3.5646025644535285e+02;
    xcrref[3] = 3.2485447959084092e+02;
    xcrref[4] = 3.2707541254659363e+03;

    //---------------------------------------------------------------------
    // Reference values of RMS-norms of solution error.
    //---------------------------------------------------------------------
    xceref[0] = 5.2969847140936856e+01;
    xceref[1] = 4.4632896115670668e+00;
    xceref[2] = 1.3122573342210174e+01;
    xceref[3] = 1.2006925323559144e+01;
    xceref[4] = 1.2459576151035986e+02;

    //---------------------------------------------------------------------
    // reference data for 162X162X162 grids after 200 time steps,
    // with DT = 1.0e-04
    //---------------------------------------------------------------------
  } else if ( (grid_points[0] == 162) && (grid_points[1] == 162) &&
              (grid_points[2] == 162) && (no_time_steps == 200) ) {

    *Class = 'C';
    dtref = 1.0e-4;

    //---------------------------------------------------------------------
    // Reference values of RMS-norms of residual.
    //---------------------------------------------------------------------
    xcrref[0] = 0.62398116551764615e+04;
    xcrref[1] = 0.50793239190423964e+03;
    xcrref[2] = 0.15423530093013596e+04;
    xcrref[3] = 0.13302387929291190e+04;
    xcrref[4] = 0.11604087428436455e+05;

    //---------------------------------------------------------------------
    // Reference values of RMS-norms of solution error.
    //---------------------------------------------------------------------
    xceref[0] = 0.16462008369091265e+03;
    xceref[1] = 0.11497107903824313e+02;
    xceref[2] = 0.41207446207461508e+02;
    xceref[3] = 0.37087651059694167e+02;
    xceref[4] = 0.36211053051841265e+03;

    //---------------------------------------------------------------------
    // reference data for 408x408x408 grids after 250 time steps,
    // with DT = 0.2e-04
    //---------------------------------------------------------------------
  } else if ( (grid_points[0] == 408) && (grid_points[1] == 408) &&
              (grid_points[2] == 408) && (no_time_steps == 250) ) {

    *Class = 'D';
    dtref = 0.2e-4;

    //---------------------------------------------------------------------
    // Reference values of RMS-norms of residual.
    //---------------------------------------------------------------------
    xcrref[0] = 0.2533188551738e+05;
    xcrref[1] = 0.2346393716980e+04;
    xcrref[2] = 0.6294554366904e+04;
    xcrref[3] = 0.5352565376030e+04;
    xcrref[4] = 0.3905864038618e+05;

    //---------------------------------------------------------------------
    // Reference values of RMS-norms of solution error.
    //---------------------------------------------------------------------
    xceref[0] = 0.3100009377557e+03;
    xceref[1] = 0.2424086324913e+02;
    xceref[2] = 0.7782212022645e+02;
    xceref[3] = 0.6835623860116e+02;
    xceref[4] = 0.6065737200368e+03;

    //---------------------------------------------------------------------
    // reference data for 1020x1020x1020 grids after 250 time steps,
    // with DT = 0.4e-05
    //---------------------------------------------------------------------
  } else if ( (grid_points[0] == 1020) && (grid_points[1] == 1020) &&
              (grid_points[2] == 1020) && (no_time_steps == 250) ) {

    *Class = 'E';
    dtref = 0.4e-5;

    //---------------------------------------------------------------------
    // Reference values of RMS-norms of residual.
    //---------------------------------------------------------------------
    xcrref[0] = 0.9795372484517e+05;
    xcrref[1] = 0.9739814511521e+04;
    xcrref[2] = 0.2467606342965e+05;
    xcrref[3] = 0.2092419572860e+05;
    xcrref[4] = 0.1392138856939e+06;

    //---------------------------------------------------------------------
    // Reference values of RMS-norms of solution error.
    //---------------------------------------------------------------------
    xceref[0] = 0.4327562208414e+03;
    xceref[1] = 0.3699051964887e+02;
    xceref[2] = 0.1089845040954e+03;
    xceref[3] = 0.9462517622043e+02;
    xceref[4] = 0.7765512765309e+03;

  } else {
    *verified = false;
  }

  //---------------------------------------------------------------------
  // verification test for residuals if gridsize is one of 
  // the defined grid sizes above (*Class != 'U')
  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  // Compute the difference of solution values and the known reference values.
  //---------------------------------------------------------------------
  for (m = 0; m < 5; m++) {
    xcrdif[m] = fabs((xcr[m]-xcrref[m])/xcrref[m]);
    xcedif[m] = fabs((xce[m]-xceref[m])/xceref[m]);
  }

  //---------------------------------------------------------------------
  // Output the comparison of computed results to known cases.
  //---------------------------------------------------------------------
  if (*Class != 'U') {
    dvmh_void_printf(" Verification being performed for class %c\n", *Class);
    dvmh_void_printf(" accuracy setting for epsilon = %20.13E\n", epsilon);
    *verified = (fabs(dt-dtref) <= epsilon);
    if (!(*verified)) {  
      *Class = 'U';
      dvmh_void_printf(" DT does not match the reference value of %15.8E\n", dtref);
    }
  } else { 
    dvmh_void_printf(" Unknown class\n");
  }

  if (*Class != 'U') {
    dvmh_void_printf(" Comparison of RMS-norms of residual\n");
  } else {
    dvmh_void_printf(" RMS-norms of residual\n");
  }

  for (m = 0; m < 5; m++) {
    if (*Class == 'U') {
      dvmh_void_printf("          %2d%20.13E\n", m+1, xcr[m]);
    } else if (xcrdif[m] <= epsilon) {
      dvmh_void_printf("          %2d%20.13E%20.13E%20.13E\n", 
          m+1, xcr[m], xcrref[m], xcrdif[m]);
    } else { 
      *verified = false;
      dvmh_void_printf(" FAILURE: %2d%20.13E%20.13E%20.13E\n",
          m+1, xcr[m], xcrref[m], xcrdif[m]);
    }
  }

  if (*Class != 'U') {
    dvmh_void_printf(" Comparison of RMS-norms of solution error\n");
  } else {
    dvmh_void_printf(" RMS-norms of solution error\n");
  }

  for (m = 0; m < 5; m++) {
    if (*Class == 'U') {
      dvmh_void_printf("          %2d%20.13E\n", m+1, xce[m]);
    } else if (xcedif[m] <= epsilon) {
      dvmh_void_printf("          %2d%20.13E%20.13E%20.13E\n", 
          m+1, xce[m], xceref[m], xcedif[m]);
    } else {
      *verified = false;
      dvmh_void_printf(" FAILURE: %2d%20.13E%20.13E%20.13E\n",
          m+1, xce[m], xceref[m], xcedif[m]);
    }
  }

  if (*Class == 'U') {
    dvmh_void_printf(" No reference values provided\n");
    dvmh_void_printf(" No verification performed\n");
  } else if (*verified) {
    dvmh_void_printf(" Verification Successful\n");
  } else {
    dvmh_void_printf(" Verification failed\n");
  }
}

