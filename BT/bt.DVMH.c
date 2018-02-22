#include <cdvmh_helpers.h>

#define DVM0C(n) ((DvmType)(n))
#define DVM0C0 DVM0C(0)
#define DVM0C1 DVM0C(1)

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
// program BT
//---------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

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


#include "../common/timers.h"
#include "../common/print_results.h"

/* common /global/ */
double elapsed_time;
int grid_points[3];
logical timeron;

/* common /constants/ */
double tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3, 
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

// to improve cache performance, grid dimensions padded by 1 
// for even number sizes only.
/* common /fields/ */

DvmType u[64];


DvmType us[64];
DvmType vs[64];
DvmType ws[64];
DvmType qs[64];
DvmType rho_i[64];
DvmType square[64];
DvmType forcing[64];



//#pragma dvm array align ([i][j][k][l] with u[i][j][k][l])
double rhs    [KMAX][JMAXP+1][IMAXP+1][5];

/* common /work_1d/ */
double cuf[PROBLEM_SIZE+1];// как такие выравнивать? по [i] with u[i][0][0][0]
double q  [PROBLEM_SIZE+1];// спросить про разницу запусков. просто файла запуск через dvm run.
double ue [PROBLEM_SIZE+1][5];
double buf[PROBLEM_SIZE+1][5];

/* common /work_lhs/ */
double fjac[PROBLEM_SIZE+1][5][5];
double njac[PROBLEM_SIZE+1][5][5];
double lhs [PROBLEM_SIZE+1][3][5][5];
double tmp1, tmp2, tmp3;

double lhs_buf[PROBLEM_SIZE+1][PROBLEM_SIZE+1][PROBLEM_SIZE+1][5][5];


int main(int argc, char *argv[])
{
    dvmh_line_C(109, "bt.c");
#ifdef _OPENMP
    dvmh_init_C(INITFLAG_OPENMP, &argc, &argv);
#else
    dvmh_init_C(0, &argc, &argv);
#endif

  int i, niter, step;
  double navg, mflops, n3;

  double tmax, t, trecs[t_last+1];
  logical verified;
  char Class;
  char *t_names[t_last+1];

  //---------------------------------------------------------------------
  // Root node reads input file (if it exists) else takes
  // defaults from parameters
  //---------------------------------------------------------------------
  FILE *fp;
  if ((fp = dvmh_fopen("timer.flag", "r")) != NULL) {
    timeron = true;
    t_names[t_total] = "total";
    t_names[t_rhsx] = "rhsx";
    t_names[t_rhsy] = "rhsy";
    t_names[t_rhsz] = "rhsz";
    t_names[t_rhs] = "rhs";
    t_names[t_xsolve] = "xsolve";
    t_names[t_ysolve] = "ysolve";
    t_names[t_zsolve] = "zsolve";
    t_names[t_rdis1] = "redist1";
    t_names[t_rdis2] = "redist2";
    t_names[t_add] = "add";
    dvmh_fclose(fp);
  } else {
    timeron = false;
  }

  dvmh_void_printf("\n\n NAS Parallel Benchmarks (NPB3.3-SER-C) - BT Benchmark\n\n");

  if ((fp = dvmh_fopen("inputbt.data", "r")) != NULL) {
    int result;
    dvmh_void_printf(" Reading from input file inputbt.data\n");
    result = dvmh_fscanf(fp, "%d", &niter);
    while (dvmh_fgetc(fp) != '\n');
    result = dvmh_fscanf(fp, "%lf", &dt);
    while (dvmh_fgetc(fp) != '\n');
    result = dvmh_fscanf(fp, "%d%d%d\n", 
        &grid_points[0], &grid_points[1], &grid_points[2]);
    dvmh_fclose(fp);
  } else {
    dvmh_void_printf(" No input file inputbt.data. Using compiled defaults\n");
    niter = NITER_DEFAULT;
    dt    = DT_DEFAULT;
    grid_points[0] = PROBLEM_SIZE;
    grid_points[1] = PROBLEM_SIZE;
    grid_points[2] = PROBLEM_SIZE;
  }

  dvmh_void_printf(" Size: %4dx%4dx%4d\n",
      grid_points[0], grid_points[1], grid_points[2]);
  dvmh_void_printf(" Iterations: %4d    dt: %10.6f\n", niter, dt);
  dvmh_void_printf("\n");

  if ( (grid_points[0] > IMAX) ||
       (grid_points[1] > JMAX) ||
       (grid_points[2] > KMAX) ) {
    dvmh_void_printf(" %d, %d, %d\n", grid_points[0], grid_points[1], grid_points[2]);
    dvmh_void_printf(" Problem size too big for compiled array sizes\n");
    dvmh_line_C(173, "bt.c");
    dvmh_exit_C( 0);
  }

  set_constants();

  for (i = 1; i <= t_last; i++) {
    timer_clear(i);
  }

  initialize();

  exact_rhs();

  //---------------------------------------------------------------------
  // do one time step to touch all code, and reinitialize
  //---------------------------------------------------------------------
  adi();
  initialize();

  for (i = 1; i <= t_last; i++) {
    timer_clear(i);
  }
  timer_start(1);

  for (step = 1; step <= niter; step++) {
    if ((step % 20) == 0 || step == 1) {
      dvmh_void_printf(" Time step %4d\n", step);
    }

    adi();
  }

  timer_stop(1);
  tmax = timer_read(1);

  verify(niter, &Class, &verified);

  n3 = 1.0*grid_points[0]*grid_points[1]*grid_points[2];
  navg = (grid_points[0]+grid_points[1]+grid_points[2])/3.0;
  if(tmax != 0.0) {
    mflops = 1.0e-6 * (double)niter *
      (3478.8 * n3 - 17655.7 * (navg*navg) + 28023.7 * navg)
      / tmax;
  } else {
    mflops = 0.0;
  }
  print_results("BT", Class, grid_points[0], 
                grid_points[1], grid_points[2], niter,
                tmax, mflops, "          floating point", 
                verified, NPBVERSION,COMPILETIME, CS1, CS2, CS3, CS4, CS5, 
                CS6, "(none)");

  //---------------------------------------------------------------------
  // More timers
  //---------------------------------------------------------------------
  if (timeron) {
    for (i = 1; i <= t_last; i++) {
      trecs[i] = timer_read(i);
    }

    if (tmax == 0.0) tmax = 1.0;
    dvmh_void_printf("  SECTION   Time (secs)\n");
    for (i = 1; i <= t_last; i++) {
      dvmh_void_printf("  %-8s:%9.3f  (%6.2f%%)\n", 
          t_names[i], trecs[i], trecs[i]*100./tmax);
      if (i == t_rhs) {
        t = trecs[t_rhsx] + trecs[t_rhsy] + trecs[t_rhsz];
        dvmh_void_printf("    --> %8s:%9.3f  (%6.2f%%)\n", "sub-rhs", t, t*100./tmax);
        t = trecs[t_rhs] - t;
        dvmh_void_printf("    --> %8s:%9.3f  (%6.2f%%)\n", "rest-rhs", t, t*100./tmax);
      } else if (i==t_zsolve) {
        t = trecs[t_zsolve] - trecs[t_rdis1] - trecs[t_rdis2];
        dvmh_void_printf("    --> %8s:%9.3f  (%6.2f%%)\n", "sub-zsol", t, t*100./tmax);
      } else if (i==t_rdis2) {
        t = trecs[t_rdis1] + trecs[t_rdis2];
        dvmh_void_printf("    --> %8s:%9.3f  (%6.2f%%)\n", "redist", t, t*100./tmax);
      }
    }
  }

  dvmh_line_C(253, "bt.c");
  dvmh_exit_C( 0);
    dvmh_exit_C(0);
}



void initCdvmhGlobals_bt_2144513821() {
    dvmh_line_C(70, "bt.c");
    dvmh_array_create_C(u, 4, -rt_DOUBLE, DVM0C(24), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1, DVM0C(5), DVM0C1, DVM0C1);
    dvmh_distribute_C(u, 4, dvmh_distribution_block_(), dvmh_distribution_replicated_(), dvmh_distribution_replicated_(), dvmh_distribution_replicated_());

    dvmh_line_C(74, "bt.c");
    dvmh_array_create_C(us, 3, -rt_DOUBLE, DVM0C(24), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1);
    dvmh_align_C(us, u, 4, dvmh_alignment_linear_C(1, 1, 0), dvmh_alignment_linear_C(2, 1, 0), dvmh_alignment_linear_C(3, 1, 0), dvmh_alignment_linear_C(-1, 0, 0));

    dvmh_line_C(76, "bt.c");
    dvmh_array_create_C(vs, 3, -rt_DOUBLE, DVM0C(24), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1);
    dvmh_align_C(vs, u, 4, dvmh_alignment_linear_C(1, 1, 0), dvmh_alignment_linear_C(2, 1, 0), dvmh_alignment_linear_C(3, 1, 0), dvmh_alignment_linear_C(-1, 0, 0));

    dvmh_line_C(78, "bt.c");
    dvmh_array_create_C(ws, 3, -rt_DOUBLE, DVM0C(24), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1);
    dvmh_align_C(ws, u, 4, dvmh_alignment_linear_C(1, 1, 0), dvmh_alignment_linear_C(2, 1, 0), dvmh_alignment_linear_C(3, 1, 0), dvmh_alignment_linear_C(-1, 0, 0));

    dvmh_line_C(80, "bt.c");
    dvmh_array_create_C(qs, 3, -rt_DOUBLE, DVM0C(24), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1);
    dvmh_align_C(qs, u, 4, dvmh_alignment_linear_C(1, 1, 0), dvmh_alignment_linear_C(2, 1, 0), dvmh_alignment_linear_C(3, 1, 0), dvmh_alignment_linear_C(-1, 0, 0));

    dvmh_line_C(82, "bt.c");
    dvmh_array_create_C(rho_i, 3, -rt_DOUBLE, DVM0C(24), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1);
    dvmh_align_C(rho_i, u, 4, dvmh_alignment_linear_C(1, 1, 0), dvmh_alignment_linear_C(2, 1, 0), dvmh_alignment_linear_C(3, 1, 0), dvmh_alignment_linear_C(-1, 0, 0));

    dvmh_line_C(84, "bt.c");
    dvmh_array_create_C(square, 3, -rt_DOUBLE, DVM0C(24), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1);
    dvmh_align_C(square, u, 4, dvmh_alignment_linear_C(1, 1, 0), dvmh_alignment_linear_C(2, 1, 0), dvmh_alignment_linear_C(3, 1, 0), dvmh_alignment_linear_C(-1, 0, 0));

    dvmh_line_C(86, "bt.c");
    dvmh_array_create_C(forcing, 4, -rt_DOUBLE, DVM0C(24), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1, DVM0C(25), DVM0C1, DVM0C1, DVM0C(5), DVM0C1, DVM0C1);
    dvmh_align_C(forcing, u, 4, dvmh_alignment_linear_C(1, 1, 0), dvmh_alignment_linear_C(2, 1, 0), dvmh_alignment_linear_C(3, 1, 0), dvmh_alignment_linear_C(-1, 0, 0));

    dvmh_data_enter_C((const void *)&elapsed_time, sizeof(elapsed_time));
    dvmh_data_enter_C((const void *)grid_points, sizeof(grid_points));
    dvmh_data_enter_C((const void *)&timeron, sizeof(timeron));
    dvmh_data_enter_C((const void *)&tx1, sizeof(tx1));
    dvmh_data_enter_C((const void *)&tx2, sizeof(tx2));
    dvmh_data_enter_C((const void *)&tx3, sizeof(tx3));
    dvmh_data_enter_C((const void *)&ty1, sizeof(ty1));
    dvmh_data_enter_C((const void *)&ty2, sizeof(ty2));
    dvmh_data_enter_C((const void *)&ty3, sizeof(ty3));
    dvmh_data_enter_C((const void *)&tz1, sizeof(tz1));
    dvmh_data_enter_C((const void *)&tz2, sizeof(tz2));
    dvmh_data_enter_C((const void *)&tz3, sizeof(tz3));
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
    dvmh_data_enter_C((const void *)&dt, sizeof(dt));
    dvmh_data_enter_C((const void *)ce, sizeof(ce));
    dvmh_data_enter_C((const void *)&dxmax, sizeof(dxmax));
    dvmh_data_enter_C((const void *)&dymax, sizeof(dymax));
    dvmh_data_enter_C((const void *)&dzmax, sizeof(dzmax));
    dvmh_data_enter_C((const void *)&xxcon1, sizeof(xxcon1));
    dvmh_data_enter_C((const void *)&xxcon2, sizeof(xxcon2));
    dvmh_data_enter_C((const void *)&xxcon3, sizeof(xxcon3));
    dvmh_data_enter_C((const void *)&xxcon4, sizeof(xxcon4));
    dvmh_data_enter_C((const void *)&xxcon5, sizeof(xxcon5));
    dvmh_data_enter_C((const void *)&dx1tx1, sizeof(dx1tx1));
    dvmh_data_enter_C((const void *)&dx2tx1, sizeof(dx2tx1));
    dvmh_data_enter_C((const void *)&dx3tx1, sizeof(dx3tx1));
    dvmh_data_enter_C((const void *)&dx4tx1, sizeof(dx4tx1));
    dvmh_data_enter_C((const void *)&dx5tx1, sizeof(dx5tx1));
    dvmh_data_enter_C((const void *)&yycon1, sizeof(yycon1));
    dvmh_data_enter_C((const void *)&yycon2, sizeof(yycon2));
    dvmh_data_enter_C((const void *)&yycon3, sizeof(yycon3));
    dvmh_data_enter_C((const void *)&yycon4, sizeof(yycon4));
    dvmh_data_enter_C((const void *)&yycon5, sizeof(yycon5));
    dvmh_data_enter_C((const void *)&dy1ty1, sizeof(dy1ty1));
    dvmh_data_enter_C((const void *)&dy2ty1, sizeof(dy2ty1));
    dvmh_data_enter_C((const void *)&dy3ty1, sizeof(dy3ty1));
    dvmh_data_enter_C((const void *)&dy4ty1, sizeof(dy4ty1));
    dvmh_data_enter_C((const void *)&dy5ty1, sizeof(dy5ty1));
    dvmh_data_enter_C((const void *)&zzcon1, sizeof(zzcon1));
    dvmh_data_enter_C((const void *)&zzcon2, sizeof(zzcon2));
    dvmh_data_enter_C((const void *)&zzcon3, sizeof(zzcon3));
    dvmh_data_enter_C((const void *)&zzcon4, sizeof(zzcon4));
    dvmh_data_enter_C((const void *)&zzcon5, sizeof(zzcon5));
    dvmh_data_enter_C((const void *)&dz1tz1, sizeof(dz1tz1));
    dvmh_data_enter_C((const void *)&dz2tz1, sizeof(dz2tz1));
    dvmh_data_enter_C((const void *)&dz3tz1, sizeof(dz3tz1));
    dvmh_data_enter_C((const void *)&dz4tz1, sizeof(dz4tz1));
    dvmh_data_enter_C((const void *)&dz5tz1, sizeof(dz5tz1));
    dvmh_data_enter_C((const void *)&dnxm1, sizeof(dnxm1));
    dvmh_data_enter_C((const void *)&dnym1, sizeof(dnym1));
    dvmh_data_enter_C((const void *)&dnzm1, sizeof(dnzm1));
    dvmh_data_enter_C((const void *)&c1c2, sizeof(c1c2));
    dvmh_data_enter_C((const void *)&c1c5, sizeof(c1c5));
    dvmh_data_enter_C((const void *)&c3c4, sizeof(c3c4));
    dvmh_data_enter_C((const void *)&c1345, sizeof(c1345));
    dvmh_data_enter_C((const void *)&conz1, sizeof(conz1));
    dvmh_data_enter_C((const void *)&c1, sizeof(c1));
    dvmh_data_enter_C((const void *)&c2, sizeof(c2));
    dvmh_data_enter_C((const void *)&c3, sizeof(c3));
    dvmh_data_enter_C((const void *)&c4, sizeof(c4));
    dvmh_data_enter_C((const void *)&c5, sizeof(c5));
    dvmh_data_enter_C((const void *)&c4dssp, sizeof(c4dssp));
    dvmh_data_enter_C((const void *)&c5dssp, sizeof(c5dssp));
    dvmh_data_enter_C((const void *)&dtdssp, sizeof(dtdssp));
    dvmh_data_enter_C((const void *)&dttx1, sizeof(dttx1));
    dvmh_data_enter_C((const void *)&dttx2, sizeof(dttx2));
    dvmh_data_enter_C((const void *)&dtty1, sizeof(dtty1));
    dvmh_data_enter_C((const void *)&dtty2, sizeof(dtty2));
    dvmh_data_enter_C((const void *)&dttz1, sizeof(dttz1));
    dvmh_data_enter_C((const void *)&dttz2, sizeof(dttz2));
    dvmh_data_enter_C((const void *)&c2dttx1, sizeof(c2dttx1));
    dvmh_data_enter_C((const void *)&c2dtty1, sizeof(c2dtty1));
    dvmh_data_enter_C((const void *)&c2dttz1, sizeof(c2dttz1));
    dvmh_data_enter_C((const void *)&comz1, sizeof(comz1));
    dvmh_data_enter_C((const void *)&comz4, sizeof(comz4));
    dvmh_data_enter_C((const void *)&comz5, sizeof(comz5));
    dvmh_data_enter_C((const void *)&comz6, sizeof(comz6));
    dvmh_data_enter_C((const void *)&c3c4tx3, sizeof(c3c4tx3));
    dvmh_data_enter_C((const void *)&c3c4ty3, sizeof(c3c4ty3));
    dvmh_data_enter_C((const void *)&c3c4tz3, sizeof(c3c4tz3));
    dvmh_data_enter_C((const void *)&c2iv, sizeof(c2iv));
    dvmh_data_enter_C((const void *)&con43, sizeof(con43));
    dvmh_data_enter_C((const void *)&con16, sizeof(con16));
    dvmh_data_enter_C((const void *)rhs, sizeof(rhs));
    dvmh_data_enter_C((const void *)cuf, sizeof(cuf));
    dvmh_data_enter_C((const void *)q, sizeof(q));
    dvmh_data_enter_C((const void *)ue, sizeof(ue));
    dvmh_data_enter_C((const void *)buf, sizeof(buf));
    dvmh_data_enter_C((const void *)fjac, sizeof(fjac));
    dvmh_data_enter_C((const void *)njac, sizeof(njac));
    dvmh_data_enter_C((const void *)lhs, sizeof(lhs));
    dvmh_data_enter_C((const void *)&tmp1, sizeof(tmp1));
    dvmh_data_enter_C((const void *)&tmp2, sizeof(tmp2));
    dvmh_data_enter_C((const void *)&tmp3, sizeof(tmp3));
    dvmh_data_enter_C((const void *)lhs_buf, sizeof(lhs_buf));
}

void finishCdvmhGlobals_bt_1111161825() {
    dvmh_delete_object_(forcing);
    dvmh_delete_object_(square);
    dvmh_delete_object_(rho_i);
    dvmh_delete_object_(qs);
    dvmh_delete_object_(ws);
    dvmh_delete_object_(vs);
    dvmh_delete_object_(us);
    dvmh_delete_object_(u);
    dvmh_data_exit_C((const void *)lhs_buf, 0);
    dvmh_data_exit_C((const void *)&tmp3, 0);
    dvmh_data_exit_C((const void *)&tmp2, 0);
    dvmh_data_exit_C((const void *)&tmp1, 0);
    dvmh_data_exit_C((const void *)lhs, 0);
    dvmh_data_exit_C((const void *)njac, 0);
    dvmh_data_exit_C((const void *)fjac, 0);
    dvmh_data_exit_C((const void *)buf, 0);
    dvmh_data_exit_C((const void *)ue, 0);
    dvmh_data_exit_C((const void *)q, 0);
    dvmh_data_exit_C((const void *)cuf, 0);
    dvmh_data_exit_C((const void *)rhs, 0);
    dvmh_data_exit_C((const void *)&con16, 0);
    dvmh_data_exit_C((const void *)&con43, 0);
    dvmh_data_exit_C((const void *)&c2iv, 0);
    dvmh_data_exit_C((const void *)&c3c4tz3, 0);
    dvmh_data_exit_C((const void *)&c3c4ty3, 0);
    dvmh_data_exit_C((const void *)&c3c4tx3, 0);
    dvmh_data_exit_C((const void *)&comz6, 0);
    dvmh_data_exit_C((const void *)&comz5, 0);
    dvmh_data_exit_C((const void *)&comz4, 0);
    dvmh_data_exit_C((const void *)&comz1, 0);
    dvmh_data_exit_C((const void *)&c2dttz1, 0);
    dvmh_data_exit_C((const void *)&c2dtty1, 0);
    dvmh_data_exit_C((const void *)&c2dttx1, 0);
    dvmh_data_exit_C((const void *)&dttz2, 0);
    dvmh_data_exit_C((const void *)&dttz1, 0);
    dvmh_data_exit_C((const void *)&dtty2, 0);
    dvmh_data_exit_C((const void *)&dtty1, 0);
    dvmh_data_exit_C((const void *)&dttx2, 0);
    dvmh_data_exit_C((const void *)&dttx1, 0);
    dvmh_data_exit_C((const void *)&dtdssp, 0);
    dvmh_data_exit_C((const void *)&c5dssp, 0);
    dvmh_data_exit_C((const void *)&c4dssp, 0);
    dvmh_data_exit_C((const void *)&c5, 0);
    dvmh_data_exit_C((const void *)&c4, 0);
    dvmh_data_exit_C((const void *)&c3, 0);
    dvmh_data_exit_C((const void *)&c2, 0);
    dvmh_data_exit_C((const void *)&c1, 0);
    dvmh_data_exit_C((const void *)&conz1, 0);
    dvmh_data_exit_C((const void *)&c1345, 0);
    dvmh_data_exit_C((const void *)&c3c4, 0);
    dvmh_data_exit_C((const void *)&c1c5, 0);
    dvmh_data_exit_C((const void *)&c1c2, 0);
    dvmh_data_exit_C((const void *)&dnzm1, 0);
    dvmh_data_exit_C((const void *)&dnym1, 0);
    dvmh_data_exit_C((const void *)&dnxm1, 0);
    dvmh_data_exit_C((const void *)&dz5tz1, 0);
    dvmh_data_exit_C((const void *)&dz4tz1, 0);
    dvmh_data_exit_C((const void *)&dz3tz1, 0);
    dvmh_data_exit_C((const void *)&dz2tz1, 0);
    dvmh_data_exit_C((const void *)&dz1tz1, 0);
    dvmh_data_exit_C((const void *)&zzcon5, 0);
    dvmh_data_exit_C((const void *)&zzcon4, 0);
    dvmh_data_exit_C((const void *)&zzcon3, 0);
    dvmh_data_exit_C((const void *)&zzcon2, 0);
    dvmh_data_exit_C((const void *)&zzcon1, 0);
    dvmh_data_exit_C((const void *)&dy5ty1, 0);
    dvmh_data_exit_C((const void *)&dy4ty1, 0);
    dvmh_data_exit_C((const void *)&dy3ty1, 0);
    dvmh_data_exit_C((const void *)&dy2ty1, 0);
    dvmh_data_exit_C((const void *)&dy1ty1, 0);
    dvmh_data_exit_C((const void *)&yycon5, 0);
    dvmh_data_exit_C((const void *)&yycon4, 0);
    dvmh_data_exit_C((const void *)&yycon3, 0);
    dvmh_data_exit_C((const void *)&yycon2, 0);
    dvmh_data_exit_C((const void *)&yycon1, 0);
    dvmh_data_exit_C((const void *)&dx5tx1, 0);
    dvmh_data_exit_C((const void *)&dx4tx1, 0);
    dvmh_data_exit_C((const void *)&dx3tx1, 0);
    dvmh_data_exit_C((const void *)&dx2tx1, 0);
    dvmh_data_exit_C((const void *)&dx1tx1, 0);
    dvmh_data_exit_C((const void *)&xxcon5, 0);
    dvmh_data_exit_C((const void *)&xxcon4, 0);
    dvmh_data_exit_C((const void *)&xxcon3, 0);
    dvmh_data_exit_C((const void *)&xxcon2, 0);
    dvmh_data_exit_C((const void *)&xxcon1, 0);
    dvmh_data_exit_C((const void *)&dzmax, 0);
    dvmh_data_exit_C((const void *)&dymax, 0);
    dvmh_data_exit_C((const void *)&dxmax, 0);
    dvmh_data_exit_C((const void *)ce, 0);
    dvmh_data_exit_C((const void *)&dt, 0);
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
    dvmh_data_exit_C((const void *)&tz3, 0);
    dvmh_data_exit_C((const void *)&tz2, 0);
    dvmh_data_exit_C((const void *)&tz1, 0);
    dvmh_data_exit_C((const void *)&ty3, 0);
    dvmh_data_exit_C((const void *)&ty2, 0);
    dvmh_data_exit_C((const void *)&ty1, 0);
    dvmh_data_exit_C((const void *)&tx3, 0);
    dvmh_data_exit_C((const void *)&tx2, 0);
    dvmh_data_exit_C((const void *)&tx1, 0);
    dvmh_data_exit_C((const void *)&timeron, 0);
    dvmh_data_exit_C((const void *)grid_points, 0);
    dvmh_data_exit_C((const void *)&elapsed_time, 0);
}

