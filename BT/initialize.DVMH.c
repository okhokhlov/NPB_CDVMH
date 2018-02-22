#include <cdvmh_helpers.h>

#define DVM0C(n) ((DvmType)(n))
#define DVM0C0 DVM0C(0)
#define DVM0C1 DVM0C(1)

void loop_initialize_52_host(DvmType *, DvmType []);
void loop_initialize_203_host(DvmType *, double * DVMH_RESTRICT, double * DVMH_RESTRICT, double * DVMH_RESTRICT, DvmType []);

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



//---------------------------------------------------------------------
// This subroutine initializes the field variable u using 
// tri-linear transfinite interpolation of the boundary values     
//---------------------------------------------------------------------
void initialize()
{
  int i, j, k, m, ix, iy, iz;
  double xi, eta, zeta, Pface[2][3][5], Pxi, Peta, Pzeta, temp[5];

  //---------------------------------------------------------------------
  // Later (in compute_rhs) we compute 1/u for every element. A few of 
  // the corner elements are not used, but it convenient (and faster) 
  // to compute the whole thing with a simple loop. Make sure those 
  // values are nonzero by initializing the whole thing here. 
  //---------------------------------------------------------------------
  do {
  	dvmh_line_C(51, "initialize.c");
  	DvmType cur_loop = dvmh_loop_create_C(0, 3, DVM0C0, DVM0C(grid_points[2] - 1), DVM0C1, DVM0C0, DVM0C(grid_points[1] - 1), DVM0C1, DVM0C0, DVM0C(grid_points[0] - 1), DVM0C1);
  	dvmh_loop_map_C(cur_loop, u, 4, dvmh_alignment_linear_C(1, 1, 0), dvmh_alignment_linear_C(2, 1, 0), dvmh_alignment_linear_C(3, 1, 0), dvmh_alignment_linear_C(-1, 0, 0));
  	dvmh_loop_register_handler_C(cur_loop, DEVICE_TYPE_HOST, HANDLER_TYPE_MASTER, dvmh_handler_func_C((DvmHandlerFunc)loop_initialize_52_host, 1, u));

  	dvmh_line_C(52, "initialize.c");
  	dvmh_loop_perform_C(cur_loop);
  	cur_loop = 0;
} while (0);


  //---------------------------------------------------------------------
  // first store the "interpolated" values everywhere on the grid    
  //---------------------------------------------------------------------
  /*
  for (k = 0; k <= grid_points[2]-1; k++) {
    zeta = (double)(k) * dnzm1;
    for (j = 0; j <= grid_points[1]-1; j++) {
      eta = (double)(j) * dnym1;
      for (i = 0; i <= grid_points[0]-1; i++) {
        xi = (double)(i) * dnxm1;

        for (ix = 0; ix < 2; ix++) {
          exact_solution((double)ix, eta, zeta, &Pface[ix][0][0]);
        }

        for (iy = 0; iy < 2; iy++) {
          exact_solution(xi, (double)iy , zeta, &Pface[iy][1][0]);
        }

        for (iz = 0; iz < 2; iz++) {
          exact_solution(xi, eta, (double)iz, &Pface[iz][2][0]);
        }

        for (m = 0; m < 5; m++) {
          Pxi   = xi   * Pface[1][0][m] + (1.0-xi)   * Pface[0][0][m];
          Peta  = eta  * Pface[1][1][m] + (1.0-eta)  * Pface[0][1][m];
          Pzeta = zeta * Pface[1][2][m] + (1.0-zeta) * Pface[0][2][m];

          u[k][j][i][m] = Pxi + Peta + Pzeta - 
                          Pxi*Peta - Pxi*Pzeta - Peta*Pzeta + 
                          Pxi*Peta*Pzeta;
        }
      }
    }
  }

  //---------------------------------------------------------------------
  // now store the exact values on the boundaries        
  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  // west face                                                  
  //---------------------------------------------------------------------
  i = 0;
  xi = 0.0;
  for (k = 0; k <= grid_points[2]-1; k++) {
    zeta = (double)(k) * dnzm1;
    for (j = 0; j <= grid_points[1]-1; j++) {
      eta = (double)(j) * dnym1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }

  //---------------------------------------------------------------------
  // east face                                                      
  //---------------------------------------------------------------------
  i = grid_points[0]-1;
  xi = 1.0;
  for (k = 0; k <= grid_points[2]-1; k++) {
    zeta = (double)(k) * dnzm1;
    for (j = 0; j <= grid_points[1]-1; j++) {
      eta = (double)(j) * dnym1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }

  //---------------------------------------------------------------------
  // south face                                                 
  //---------------------------------------------------------------------
  j = 0;
  eta = 0.0;
  for (k = 0; k <= grid_points[2]-1; k++) {
    zeta = (double)(k) * dnzm1;
    for (i = 0; i <= grid_points[0]-1; i++) {
      xi = (double)(i) * dnxm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }

  //---------------------------------------------------------------------
  // north face                                    
  //---------------------------------------------------------------------
  j = grid_points[1]-1;
  eta = 1.0;
  for (k = 0; k <= grid_points[2]-1; k++) {
    zeta = (double)(k) * dnzm1;
    for (i = 0; i <= grid_points[0]-1; i++) {
      xi = (double)(i) * dnxm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }

  //---------------------------------------------------------------------
  // bottom face                                       
  //---------------------------------------------------------------------
  k = 0;
  zeta = 0.0;
  for (j = 0; j <= grid_points[1]-1; j++) {
    eta = (double)(j) * dnym1;
    for (i =0; i <= grid_points[0]-1; i++) {
      xi = (double)(i) *dnxm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }

  //---------------------------------------------------------------------
  // top face     
  //---------------------------------------------------------------------
  k = grid_points[2]-1;
  zeta = 1.0;
  for (j = 0; j <= grid_points[1]-1; j++) {
    eta = (double)(j) * dnym1;
    for (i = 0; i <= grid_points[0]-1; i++) {
      xi = (double)(i) * dnxm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }
  
  */
  
  
  
  do {
  	dvmh_line_C(202, "initialize.c");
  	DvmType cur_loop = dvmh_loop_create_C(0, 3, DVM0C0, DVM0C(grid_points[2] - 1), DVM0C1, DVM0C0, DVM0C(grid_points[1] - 1), DVM0C1, DVM0C0, DVM0C(grid_points[0] - 1), DVM0C1);
  	dvmh_loop_map_C(cur_loop, u, 4, dvmh_alignment_linear_C(1, 1, 0), dvmh_alignment_linear_C(2, 1, 0), dvmh_alignment_linear_C(3, 1, 0), dvmh_alignment_linear_C(-1, 0, 0));
  	dvmh_loop_register_handler_C(cur_loop, DEVICE_TYPE_HOST, HANDLER_TYPE_MASTER, dvmh_handler_func_C((DvmHandlerFunc)loop_initialize_203_host, 4, (const void *)&dnxm1, (const void *)&dnym1, (const void *)&dnzm1, u));

  	dvmh_line_C(203, "initialize.c");
  	dvmh_loop_perform_C(cur_loop);
  	cur_loop = 0;
} while (0);



  i = 0;
  xi = 0.0;
  for (k = 0; k <= grid_points[2]-1; k++) {
    for (j = 0; j <= grid_points[1]-1; j++) {
	  zeta = (double)(k) * dnzm1;
      eta = (double)(j) * dnym1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        (dvmh_get_own_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = temp[m]) : ((void)0));
      }
    }
  }

  i = grid_points[0]-1;
  xi = 1.0;
  for (k = 0; k <= grid_points[2]-1; k++) {  
    for (j = 0; j <= grid_points[1]-1; j++) {
      zeta = (double)(k) * dnzm1;
      eta = (double)(j) * dnym1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        (dvmh_get_own_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = temp[m]) : ((void)0));
      }
    }
  }

  j = 0;
  eta = 0.0;
  for (k = 0; k <= grid_points[2]-1; k++) { 
    for (i = 0; i <= grid_points[0]-1; i++) {
	  zeta = (double)(k) * dnzm1;
      xi = (double)(i) * dnxm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        (dvmh_get_own_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = temp[m]) : ((void)0));
      }
    }
  }

  j = grid_points[1]-1;
  eta = 1.0;
  for (k = 0; k <= grid_points[2]-1; k++) { 
    for (i = 0; i <= grid_points[0]-1; i++) {
	  zeta = (double)(k) * dnzm1;
      xi = (double)(i) * dnxm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        (dvmh_get_own_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = temp[m]) : ((void)0));
      }
    }
  }

  k = 0;
  zeta = 0.0;
  for (j = 0; j <= grid_points[1]-1; j++) {
    eta = (double)(j) * dnym1;
    for (i =0; i <= grid_points[0]-1; i++) {
      xi = (double)(i) *dnxm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        (dvmh_get_own_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = temp[m]) : ((void)0));
      }
    }
  }

  k = grid_points[2]-1;
  zeta = 1.0;
  for (j = 0; j <= grid_points[1]-1; j++) {
    for (i = 0; i <= grid_points[0]-1; i++) {
	  eta = (double)(j) * dnym1;
      xi = (double)(i) * dnxm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        (dvmh_get_own_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(u, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = temp[m]) : ((void)0));
      }
    }
  }
}


void lhsinit(double lhs[][3][5][5], int size)
{
  int i, m, n;

  i = size;
  //---------------------------------------------------------------------
  // zero the whole left hand side for starters
  //---------------------------------------------------------------------
  for (n = 0; n < 5; n++) {
    for (m = 0; m < 5; m++) {
      lhs[0][0][n][m] = 0.0;
      lhs[0][1][n][m] = 0.0;
      lhs[0][2][n][m] = 0.0;
      lhs[i][0][n][m] = 0.0;
      lhs[i][1][n][m] = 0.0;
      lhs[i][2][n][m] = 0.0;
    }
  }

  //---------------------------------------------------------------------
  // next, set all diagonal values to 1. This is overkill, but convenient
  //---------------------------------------------------------------------
  for (m = 0; m < 5; m++) {
    lhs[0][1][m][m] = 1.0;
    lhs[i][1][m][m] = 1.0;
  }
}


void loop_initialize_52_host(DvmType *pLoopRef, DvmType u_hdr[]) {
	/* Loop reference and device number */
	DvmType loop_ref = *pLoopRef;
	DvmType device_num = dvmh_loop_get_device_num_C(loop_ref);
	/* Parameters */
	double (* DVMH_RESTRICT u)[u_hdr[1]/u_hdr[2]][u_hdr[2]/u_hdr[3]][u_hdr[3]] = (double (*)[u_hdr[1]/u_hdr[2]][u_hdr[2]/u_hdr[3]][u_hdr[3]])dvmh_get_natural_base_C(device_num, u_hdr);
	/* Supplementary variables for loop handling */
	DvmType boundsLow[3], boundsHigh[3], loopSteps[3];
	/* User variables - loop index variables and other private variables */
	int i;
	int j;
	int k;
	int m;

	dvmh_loop_fill_bounds_C(loop_ref, boundsLow, boundsHigh, loopSteps);

	for (k = boundsLow[0]; k <= boundsHigh[0]; k++)
		for (j = boundsLow[1]; j <= boundsHigh[1]; j++)
			for (i = boundsLow[2]; i <= boundsHigh[2]; i++)
			{
			    for (m = 0; m < 5; m++) {
			        u[k][j][i][m] = 1.;
			    }
			}
}

void loop_initialize_203_host(DvmType *pLoopRef, double * DVMH_RESTRICT dnxm1_ptr, double * DVMH_RESTRICT dnym1_ptr, double * DVMH_RESTRICT dnzm1_ptr, DvmType u_hdr[]) {
	/* Loop reference and device number */
	DvmType loop_ref = *pLoopRef;
	DvmType device_num = dvmh_loop_get_device_num_C(loop_ref);
	/* Parameters */
	double dnxm1 = *dnxm1_ptr;
	double dnym1 = *dnym1_ptr;
	double dnzm1 = *dnzm1_ptr;
	double (* DVMH_RESTRICT u)[u_hdr[1]/u_hdr[2]][u_hdr[2]/u_hdr[3]][u_hdr[3]] = (double (*)[u_hdr[1]/u_hdr[2]][u_hdr[2]/u_hdr[3]][u_hdr[3]])dvmh_get_natural_base_C(device_num, u_hdr);
	/* Supplementary variables for loop handling */
	DvmType boundsLow[3], boundsHigh[3], loopSteps[3];
	/* User variables - loop index variables and other private variables */
	int i;
	int j;
	int k;
	int m;
	int ix;
	int iy;
	int iz;
	double xi;
	double eta;
	double zeta;
	double Pface[2][3][5];
	double Pxi;
	double Peta;
	double Pzeta;
	double temp[5];

	dvmh_loop_fill_bounds_C(loop_ref, boundsLow, boundsHigh, loopSteps);

	for (k = boundsLow[0]; k <= boundsHigh[0]; k++)
		for (j = boundsLow[1]; j <= boundsHigh[1]; j++)
			for (i = boundsLow[2]; i <= boundsHigh[2]; i++)
			{
			    zeta = (double)(k) * dnzm1;
			    eta = (double)(j) * dnym1;
			    xi = (double)(i) * dnxm1;
			    for (ix = 0; ix < 2; ix++) {
			        exact_solution((double)ix, eta, zeta, &Pface[ix][0][0]);
			    }
			    for (iy = 0; iy < 2; iy++) {
			        exact_solution(xi, (double)iy, zeta, &Pface[iy][1][0]);
			    }
			    for (iz = 0; iz < 2; iz++) {
			        exact_solution(xi, eta, (double)iz, &Pface[iz][2][0]);
			    }
			    for (m = 0; m < 5; m++) {
			        Pxi = xi * Pface[1][0][m] + (1. - xi) * Pface[0][0][m];
			        Peta = eta * Pface[1][1][m] + (1. - eta) * Pface[0][1][m];
			        Pzeta = zeta * Pface[1][2][m] + (1. - zeta) * Pface[0][2][m];
			        u[k][j][i][m] = Pxi + Peta + Pzeta - Pxi * Peta - Pxi * Pzeta - Peta * Pzeta + Pxi * Peta * Pzeta;
			    }
			}
}

