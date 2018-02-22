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
//   program applu
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//
//   driver for the performance evaluation of the solver for
//   five coupled parabolic/elliptic partial differential equations.
//
//---------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

#include "../common/timers.h"
#include "../common/print_results.h"


//---------------------------------------------------------------------
// grid
//---------------------------------------------------------------------
/* common/cgcon/ */
double dxi, deta, dzeta;
double tx1, tx2, tx3;
double ty1, ty2, ty3;
double tz1, tz2, tz3;
int nx, ny, nz;
int nx0, ny0, nz0;
int ist, iend;
int jst, jend;
int ii1, ii2;
int ji1, ji2;
int ki1, ki2;

//---------------------------------------------------------------------
// dissipation
//---------------------------------------------------------------------
/* common/disp/ */
double dx1, dx2, dx3, dx4, dx5;
double dy1, dy2, dy3, dy4, dy5;
double dz1, dz2, dz3, dz4, dz5;
double dssp;

//---------------------------------------------------------------------
// field variables and residuals
// to improve cache performance, second two dimensions padded by 1 
// for even number sizes only.
// Note: corresponding array (called "v") in routines blts, buts, 
// and l2norm are similarly padded
//---------------------------------------------------------------------
/* common/cvar/ */
double u    [ISIZ3][ISIZ2/2*2+1][ISIZ1/2*2+1][5];
double rsd  [ISIZ3][ISIZ2/2*2+1][ISIZ1/2*2+1][5];
double frct [ISIZ3][ISIZ2/2*2+1][ISIZ1/2*2+1][5];
double flux [ISIZ1][5];
double qs   [ISIZ3][ISIZ2/2*2+1][ISIZ1/2*2+1];
double rho_i[ISIZ3][ISIZ2/2*2+1][ISIZ1/2*2+1];

//---------------------------------------------------------------------
// output control parameters
//---------------------------------------------------------------------
/* common/cprcon/ */
int ipr, inorm;

//---------------------------------------------------------------------
// newton-raphson iteration control parameters
//---------------------------------------------------------------------
/* common/ctscon/ */
double dt, omega, tolrsd[5], rsdnm[5], errnm[5], frc, ttotal;
int itmax, invert;

/* common/cjac/ */
double a[ISIZ2][ISIZ1/2*2+1][5][5];
double b[ISIZ2][ISIZ1/2*2+1][5][5];
double c[ISIZ2][ISIZ1/2*2+1][5][5];
double d[ISIZ2][ISIZ1/2*2+1][5][5];


//---------------------------------------------------------------------
// coefficients of the exact solution
//---------------------------------------------------------------------
/* common/cexact/ */
double ce[5][13];


//---------------------------------------------------------------------
// timers
//---------------------------------------------------------------------
/* common/timer/ */
double maxtime;
logical timeron;


int main(int argc, char *argv[])
{
    dvmh_line_C(129, "lu.c");
#ifdef _OPENMP
    dvmh_init_C(INITFLAG_OPENMP, &argc, &argv);
#else
    dvmh_init_C(0, &argc, &argv);
#endif

  char Class;
  logical verified;
  double mflops;

  double t, tmax, trecs[t_last+1];
  int i;
  char *t_names[t_last+1];

  //---------------------------------------------------------------------
  // Setup info for timers
  //---------------------------------------------------------------------
  FILE *fp;
  if ((fp = dvmh_fopen("timer.flag", "r")) != NULL) {
    timeron = true;
    t_names[t_total] = "total";
    t_names[t_rhsx] = "rhsx";
    t_names[t_rhsy] = "rhsy";
    t_names[t_rhsz] = "rhsz";
    t_names[t_rhs] = "rhs";
    t_names[t_jacld] = "jacld";
    t_names[t_blts] = "blts";
    t_names[t_jacu] = "jacu";
    t_names[t_buts] = "buts";
    t_names[t_add] = "add";
    t_names[t_l2norm] = "l2norm";
    dvmh_fclose(fp);
  } else {
    timeron = false;
  }

  //---------------------------------------------------------------------
  // read input data
  //---------------------------------------------------------------------
  read_input();

  //---------------------------------------------------------------------
  // set up domain sizes
  //---------------------------------------------------------------------
  domain();

  //---------------------------------------------------------------------
  // set up coefficients
  //---------------------------------------------------------------------
  setcoeff();

  //---------------------------------------------------------------------
  // set the boundary values for dependent variables
  //---------------------------------------------------------------------
  setbv();

  //---------------------------------------------------------------------
  // set the initial values for dependent variables
  //---------------------------------------------------------------------
  setiv();

  //---------------------------------------------------------------------
  // compute the forcing term based on prescribed exact solution
  //---------------------------------------------------------------------
  erhs();

  //---------------------------------------------------------------------
  // perform one SSOR iteration to touch all pages
  //---------------------------------------------------------------------
  ssor(1);

  //---------------------------------------------------------------------
  // reset the boundary and initial values
  //---------------------------------------------------------------------
  setbv();
  setiv();

  //---------------------------------------------------------------------
  // perform the SSOR iterations
  //---------------------------------------------------------------------
  ssor(itmax);

  //---------------------------------------------------------------------
  // compute the solution error
  //---------------------------------------------------------------------
  error();

  //---------------------------------------------------------------------
  // compute the surface integral
  //---------------------------------------------------------------------
  pintgr();

  //---------------------------------------------------------------------
  // verification test
  //---------------------------------------------------------------------
  verify ( rsdnm, errnm, frc, &Class, &verified );
  mflops = (double)itmax * (1984.77 * (double)nx0
      * (double)ny0
      * (double)nz0
      - 10923.3 * pow(((double)(nx0+ny0+nz0)/3.0), 2.0) 
      + 27770.9 * (double)(nx0+ny0+nz0)/3.0
      - 144010.0)
    / (maxtime*1000000.0);

  print_results("LU", Class, nx0,
                ny0, nz0, itmax,
                maxtime, mflops, "          floating point", verified, 
                NPBVERSION, COMPILETIME, CS1, CS2, CS3, CS4, CS5, CS6, 
                "(none)");

  //---------------------------------------------------------------------
  // More timers
  //---------------------------------------------------------------------
  if (timeron) {
    for (i = 1; i <= t_last; i++) {
      trecs[i] = timer_read(i);
    }
    tmax = maxtime;
    if (tmax == 0.0) tmax = 1.0;

    dvmh_void_printf("  SECTION     Time (secs)\n");
    for (i = 1; i <= t_last; i++) {
      dvmh_void_printf("  %-8s:%9.3f  (%6.2f%%)\n",
          t_names[i], trecs[i], trecs[i]*100./tmax);
      if (i == t_rhs) {
        t = trecs[t_rhsx] + trecs[t_rhsy] + trecs[t_rhsz];
        dvmh_void_printf("     --> %8s:%9.3f  (%6.2f%%)\n", "sub-rhs", t, t*100./tmax);
        t = trecs[i] - t;
        dvmh_void_printf("     --> %8s:%9.3f  (%6.2f%%)\n", "rest-rhs", t, t*100./tmax);
      }
    }
  }

  dvmh_line_C(258, "lu.c");
  dvmh_exit_C( 0);
    dvmh_exit_C(0);
}



void initCdvmhGlobals_lu_393000618() {
    dvmh_data_enter_C((const void *)&dxi, sizeof(dxi));
    dvmh_data_enter_C((const void *)&deta, sizeof(deta));
    dvmh_data_enter_C((const void *)&dzeta, sizeof(dzeta));
    dvmh_data_enter_C((const void *)&tx1, sizeof(tx1));
    dvmh_data_enter_C((const void *)&tx2, sizeof(tx2));
    dvmh_data_enter_C((const void *)&tx3, sizeof(tx3));
    dvmh_data_enter_C((const void *)&ty1, sizeof(ty1));
    dvmh_data_enter_C((const void *)&ty2, sizeof(ty2));
    dvmh_data_enter_C((const void *)&ty3, sizeof(ty3));
    dvmh_data_enter_C((const void *)&tz1, sizeof(tz1));
    dvmh_data_enter_C((const void *)&tz2, sizeof(tz2));
    dvmh_data_enter_C((const void *)&tz3, sizeof(tz3));
    dvmh_data_enter_C((const void *)&nx, sizeof(nx));
    dvmh_data_enter_C((const void *)&ny, sizeof(ny));
    dvmh_data_enter_C((const void *)&nz, sizeof(nz));
    dvmh_data_enter_C((const void *)&nx0, sizeof(nx0));
    dvmh_data_enter_C((const void *)&ny0, sizeof(ny0));
    dvmh_data_enter_C((const void *)&nz0, sizeof(nz0));
    dvmh_data_enter_C((const void *)&ist, sizeof(ist));
    dvmh_data_enter_C((const void *)&iend, sizeof(iend));
    dvmh_data_enter_C((const void *)&jst, sizeof(jst));
    dvmh_data_enter_C((const void *)&jend, sizeof(jend));
    dvmh_data_enter_C((const void *)&ii1, sizeof(ii1));
    dvmh_data_enter_C((const void *)&ii2, sizeof(ii2));
    dvmh_data_enter_C((const void *)&ji1, sizeof(ji1));
    dvmh_data_enter_C((const void *)&ji2, sizeof(ji2));
    dvmh_data_enter_C((const void *)&ki1, sizeof(ki1));
    dvmh_data_enter_C((const void *)&ki2, sizeof(ki2));
    dvmh_data_enter_C((const void *)&dx1, sizeof(dx1));
    dvmh_data_enter_C((const void *)&dx2, sizeof(dx2));
    dvmh_data_enter_C((const void *)&dx3, sizeof(dx3));
    dvmh_data_enter_C((const void *)&dx4, sizeof(dx4));
    dvmh_data_enter_C((const void *)&dx5, sizeof(dx5));
    dvmh_data_enter_C((const void *)&dy1, sizeof(dy1));
    dvmh_data_enter_C((const void *)&dy2, sizeof(dy2));
    dvmh_data_enter_C((const void *)&dy3, sizeof(dy3));
    dvmh_data_enter_C((const void *)&dy4, sizeof(dy4));
    dvmh_data_enter_C((const void *)&dy5, sizeof(dy5));
    dvmh_data_enter_C((const void *)&dz1, sizeof(dz1));
    dvmh_data_enter_C((const void *)&dz2, sizeof(dz2));
    dvmh_data_enter_C((const void *)&dz3, sizeof(dz3));
    dvmh_data_enter_C((const void *)&dz4, sizeof(dz4));
    dvmh_data_enter_C((const void *)&dz5, sizeof(dz5));
    dvmh_data_enter_C((const void *)&dssp, sizeof(dssp));
    dvmh_data_enter_C((const void *)u, sizeof(u));
    dvmh_data_enter_C((const void *)rsd, sizeof(rsd));
    dvmh_data_enter_C((const void *)frct, sizeof(frct));
    dvmh_data_enter_C((const void *)flux, sizeof(flux));
    dvmh_data_enter_C((const void *)qs, sizeof(qs));
    dvmh_data_enter_C((const void *)rho_i, sizeof(rho_i));
    dvmh_data_enter_C((const void *)&ipr, sizeof(ipr));
    dvmh_data_enter_C((const void *)&inorm, sizeof(inorm));
    dvmh_data_enter_C((const void *)&dt, sizeof(dt));
    dvmh_data_enter_C((const void *)&omega, sizeof(omega));
    dvmh_data_enter_C((const void *)tolrsd, sizeof(tolrsd));
    dvmh_data_enter_C((const void *)rsdnm, sizeof(rsdnm));
    dvmh_data_enter_C((const void *)errnm, sizeof(errnm));
    dvmh_data_enter_C((const void *)&frc, sizeof(frc));
    dvmh_data_enter_C((const void *)&ttotal, sizeof(ttotal));
    dvmh_data_enter_C((const void *)&itmax, sizeof(itmax));
    dvmh_data_enter_C((const void *)&invert, sizeof(invert));
    dvmh_data_enter_C((const void *)a, sizeof(a));
    dvmh_data_enter_C((const void *)b, sizeof(b));
    dvmh_data_enter_C((const void *)c, sizeof(c));
    dvmh_data_enter_C((const void *)d, sizeof(d));
    dvmh_data_enter_C((const void *)ce, sizeof(ce));
    dvmh_data_enter_C((const void *)&maxtime, sizeof(maxtime));
    dvmh_data_enter_C((const void *)&timeron, sizeof(timeron));
}

void finishCdvmhGlobals_lu_326655107() {
    dvmh_data_exit_C((const void *)&timeron, 0);
    dvmh_data_exit_C((const void *)&maxtime, 0);
    dvmh_data_exit_C((const void *)ce, 0);
    dvmh_data_exit_C((const void *)d, 0);
    dvmh_data_exit_C((const void *)c, 0);
    dvmh_data_exit_C((const void *)b, 0);
    dvmh_data_exit_C((const void *)a, 0);
    dvmh_data_exit_C((const void *)&invert, 0);
    dvmh_data_exit_C((const void *)&itmax, 0);
    dvmh_data_exit_C((const void *)&ttotal, 0);
    dvmh_data_exit_C((const void *)&frc, 0);
    dvmh_data_exit_C((const void *)errnm, 0);
    dvmh_data_exit_C((const void *)rsdnm, 0);
    dvmh_data_exit_C((const void *)tolrsd, 0);
    dvmh_data_exit_C((const void *)&omega, 0);
    dvmh_data_exit_C((const void *)&dt, 0);
    dvmh_data_exit_C((const void *)&inorm, 0);
    dvmh_data_exit_C((const void *)&ipr, 0);
    dvmh_data_exit_C((const void *)rho_i, 0);
    dvmh_data_exit_C((const void *)qs, 0);
    dvmh_data_exit_C((const void *)flux, 0);
    dvmh_data_exit_C((const void *)frct, 0);
    dvmh_data_exit_C((const void *)rsd, 0);
    dvmh_data_exit_C((const void *)u, 0);
    dvmh_data_exit_C((const void *)&dssp, 0);
    dvmh_data_exit_C((const void *)&dz5, 0);
    dvmh_data_exit_C((const void *)&dz4, 0);
    dvmh_data_exit_C((const void *)&dz3, 0);
    dvmh_data_exit_C((const void *)&dz2, 0);
    dvmh_data_exit_C((const void *)&dz1, 0);
    dvmh_data_exit_C((const void *)&dy5, 0);
    dvmh_data_exit_C((const void *)&dy4, 0);
    dvmh_data_exit_C((const void *)&dy3, 0);
    dvmh_data_exit_C((const void *)&dy2, 0);
    dvmh_data_exit_C((const void *)&dy1, 0);
    dvmh_data_exit_C((const void *)&dx5, 0);
    dvmh_data_exit_C((const void *)&dx4, 0);
    dvmh_data_exit_C((const void *)&dx3, 0);
    dvmh_data_exit_C((const void *)&dx2, 0);
    dvmh_data_exit_C((const void *)&dx1, 0);
    dvmh_data_exit_C((const void *)&ki2, 0);
    dvmh_data_exit_C((const void *)&ki1, 0);
    dvmh_data_exit_C((const void *)&ji2, 0);
    dvmh_data_exit_C((const void *)&ji1, 0);
    dvmh_data_exit_C((const void *)&ii2, 0);
    dvmh_data_exit_C((const void *)&ii1, 0);
    dvmh_data_exit_C((const void *)&jend, 0);
    dvmh_data_exit_C((const void *)&jst, 0);
    dvmh_data_exit_C((const void *)&iend, 0);
    dvmh_data_exit_C((const void *)&ist, 0);
    dvmh_data_exit_C((const void *)&nz0, 0);
    dvmh_data_exit_C((const void *)&ny0, 0);
    dvmh_data_exit_C((const void *)&nx0, 0);
    dvmh_data_exit_C((const void *)&nz, 0);
    dvmh_data_exit_C((const void *)&ny, 0);
    dvmh_data_exit_C((const void *)&nx, 0);
    dvmh_data_exit_C((const void *)&tz3, 0);
    dvmh_data_exit_C((const void *)&tz2, 0);
    dvmh_data_exit_C((const void *)&tz1, 0);
    dvmh_data_exit_C((const void *)&ty3, 0);
    dvmh_data_exit_C((const void *)&ty2, 0);
    dvmh_data_exit_C((const void *)&ty1, 0);
    dvmh_data_exit_C((const void *)&tx3, 0);
    dvmh_data_exit_C((const void *)&tx2, 0);
    dvmh_data_exit_C((const void *)&tx1, 0);
    dvmh_data_exit_C((const void *)&dzeta, 0);
    dvmh_data_exit_C((const void *)&deta, 0);
    dvmh_data_exit_C((const void *)&dxi, 0);
}

