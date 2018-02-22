#include <cdvmh_helpers.h>

#define DVM0C(n) ((DvmType)(n))
#define DVM0C0 DVM0C(0)
#define DVM0C1 DVM0C(1)
#define DVM0C2 DVM0C(2)
#define DVM0C3 DVM0C(3)

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
// compute the right hand side based on exact solution
//---------------------------------------------------------------------
void exact_rhs()
{
  double dtemp[5], xi, eta, zeta, dtpp;
  int m, i, j, k, ip1, im1, jp1, jm1, km1, kp1;
  int z, a = 0,b = 0;
  //---------------------------------------------------------------------
  // initialize                                  
  //---------------------------------------------------------------------
  
  //#pragma dvm parallel ([k][j][i] on forcing[k][j][i][])
  for (k = 0; k <= grid_points[2]-1; k++) {
    for (j = 0; j <= grid_points[1]-1; j++) {
      for (i = 0; i <= grid_points[0]-1; i++) {
        for (m = 0; m < 5; m++) {
          (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = 0.0) : ((void)0));
        }
      }
    }
  }
  

  
  

  //---------------------------------------------------------------------
  // xi-direction flux differences                      
  //---------------------------------------------------------------------
  /*for (k = 1; k <= grid_points[2]-2; k++) {
    zeta = (double)(k) * dnzm1;
    for (j = 1; j <= grid_points[1]-2; j++) {
      eta = (double)(j) * dnym1;

      for (i = 0; i <= grid_points[0]-1; i++) {
        xi = (double)(i) * dnxm1;

        exact_solution(xi, eta, zeta, dtemp);
        for (m = 0; m < 5; m++) {
          ue[i][m] = dtemp[m];
        }

        dtpp = 1.0 / dtemp[0];

        for (m = 1; m < 5; m++) {
          buf[i][m] = dtpp * dtemp[m];
        }

        cuf[i]    = buf[i][1] * buf[i][1];
        buf[i][0] = cuf[i] + buf[i][2] * buf[i][2] + buf[i][3] * buf[i][3];
        q[i] = 0.5*(buf[i][1]*ue[i][1] + buf[i][2]*ue[i][2] +
                    buf[i][3]*ue[i][3]);
      }

      for (i = 1; i <= grid_points[0]-2; i++) {
        im1 = i-1;
        ip1 = i+1;

        forcing[k][j][i][0] = forcing[k][j][i][0] -
          tx2*( ue[ip1][1]-ue[im1][1] )+
          dx1tx1*(ue[ip1][0]-2.0*ue[i][0]+ue[im1][0]);

        forcing[k][j][i][1] = forcing[k][j][i][1] - tx2 * (
            (ue[ip1][1]*buf[ip1][1]+c2*(ue[ip1][4]-q[ip1]))-
            (ue[im1][1]*buf[im1][1]+c2*(ue[im1][4]-q[im1])))+
          xxcon1*(buf[ip1][1]-2.0*buf[i][1]+buf[im1][1])+
          dx2tx1*( ue[ip1][1]-2.0* ue[i][1]+ue[im1][1]);

        forcing[k][j][i][2] = forcing[k][j][i][2] - tx2 * (
            ue[ip1][2]*buf[ip1][1]-ue[im1][2]*buf[im1][1])+
          xxcon2*(buf[ip1][2]-2.0*buf[i][2]+buf[im1][2])+
          dx3tx1*( ue[ip1][2]-2.0*ue[i][2] +ue[im1][2]);

        forcing[k][j][i][3] = forcing[k][j][i][3] - tx2*(
            ue[ip1][3]*buf[ip1][1]-ue[im1][3]*buf[im1][1])+
          xxcon2*(buf[ip1][3]-2.0*buf[i][3]+buf[im1][3])+
          dx4tx1*( ue[ip1][3]-2.0* ue[i][3]+ ue[im1][3]);

        forcing[k][j][i][4] = forcing[k][j][i][4] - tx2*(
            buf[ip1][1]*(c1*ue[ip1][4]-c2*q[ip1])-
            buf[im1][1]*(c1*ue[im1][4]-c2*q[im1]))+
          0.5*xxcon3*(buf[ip1][0]-2.0*buf[i][0]+
              buf[im1][0])+
          xxcon4*(cuf[ip1]-2.0*cuf[i]+cuf[im1])+
          xxcon5*(buf[ip1][4]-2.0*buf[i][4]+buf[im1][4])+
          dx5tx1*( ue[ip1][4]-2.0* ue[i][4]+ ue[im1][4]);
      }

      //---------------------------------------------------------------------
      // Fourth-order dissipation                         
      //---------------------------------------------------------------------
      for (m = 0; m < 5; m++) {
        i = 1;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (5.0*ue[i][m] - 4.0*ue[i+1][m] +ue[i+2][m]);
        i = 2;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (-4.0*ue[i-1][m] + 6.0*ue[i][m] -
            4.0*ue[i+1][m] +     ue[i+2][m]);
      }

      for (i = 3; i <= grid_points[0]-4; i++) {
        for (m = 0; m < 5; m++) {
          forcing[k][j][i][m] = forcing[k][j][i][m] - dssp*
            (ue[i-2][m] - 4.0*ue[i-1][m] +
             6.0*ue[i][m] - 4.0*ue[i+1][m] + ue[i+2][m]);
        }
      }

      for (m = 0; m < 5; m++) {
        i = grid_points[0]-3;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (ue[i-2][m] - 4.0*ue[i-1][m] +
           6.0*ue[i][m] - 4.0*ue[i+1][m]);
        i = grid_points[0]-2;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (ue[i-2][m] - 4.0*ue[i-1][m] + 5.0*ue[i][m]);
      }
    }
  }
  */

  for (k = 1; k <= grid_points[2]-2; k++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (i = 1; i <= grid_points[0]-2; i++) {
		zeta = (double)(k) * dnzm1;
		eta = (double)(j) * dnym1;
		
		
		if(i == 1){
			a = -1;
			b = 2;
		}
		if(i == 2){
			a = -2;
			b = 2;
		}
		if((i >= 3)&&(i <= grid_points[0]-4)){
			a = -2;
			b = 2;
		}
		if(i == grid_points[0]-3){
			a = -2;
			b = 2;
		}
		if(i == grid_points[0]-2){
			a = -2;
			b = 1;
		}
		
		for(z = a; z <= b; z++){
			xi = (double)(i+z) * dnxm1;
			
			
			exact_solution(xi, eta, zeta, dtemp);
			for (m = 0; m < 5; m++) {
			  ue[i + z][m] = dtemp[m];
			}

			dtpp = 1.0 / dtemp[0];

			for (m = 1; m < 5; m++) {
			  buf[i + z][m] = dtpp * dtemp[m];
			}

			cuf[i + z]    = buf[i + z][1] * buf[i + z][1];
			buf[i + z][0] = cuf[i + z] + buf[i + z][2] * buf[i + z][2] + buf[i + z][3] * buf[i + z][3];
			q[i + z] = 0.5*(buf[i + z][1]*ue[i + z][1] + buf[i + z][2]*ue[i + z][2] +
						buf[i + z][3]*ue[i + z][3]);
		}			
		///////////////////
		/*
		xi = (double)(0) * dnxm1;
		
		
        exact_solution(xi, eta, zeta, dtemp);
        for (m = 0; m < 5; m++) {
          ue[0][m] = dtemp[m];
        }

        dtpp = 1.0 / dtemp[0];

        for (m = 1; m < 5; m++) {
          buf[0][m] = dtpp * dtemp[m];
        }

        cuf[0]    = buf[0][1] * buf[0][1];
        buf[0][0] = cuf[0] + buf[0][2] * buf[0][2] + buf[0][3] * buf[0][3];
        q[0] = 0.5*(buf[0][1]*ue[0][1] + buf[0][2]*ue[0][2] +
                    buf[0][3]*ue[0][3]);
					
		xi = (double)(grid_points[0]-1) * dnxm1;
		
		
        exact_solution(xi, eta, zeta, dtemp);
        for (m = 0; m < 5; m++) {
          ue[grid_points[0]-1][m] = dtemp[m];
        }

        dtpp = 1.0 / dtemp[0];

        for (m = 1; m < 5; m++) {
          buf[grid_points[0]-1][m] = dtpp * dtemp[m];
        }

        cuf[grid_points[0]-1]    = buf[grid_points[0]-1][1] * buf[grid_points[0]-1][1];
        buf[grid_points[0]-1][0] = cuf[grid_points[0]-1] + buf[grid_points[0]-1][2] * buf[grid_points[0]-1][2] + buf[grid_points[0]-1][3] * buf[grid_points[0]-1][3];
        q[grid_points[0]-1] = 0.5*(buf[grid_points[0]-1][1]*ue[grid_points[0]-1][1] + buf[grid_points[0]-1][2]*ue[grid_points[0]-1][2] +
                    buf[grid_points[0]-1][3]*ue[grid_points[0]-1][3]);
		*/
		////////////////////
					
					
					
					
					
					
					
					
					
        

      //for (i = 1; i <= grid_points[0]-2; i++) {
        im1 = i-1;
        ip1 = i+1;

        (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C0) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C0)) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(0))) -
          tx2*( ue[ip1][1]-ue[im1][1] )+
          dx1tx1*(ue[ip1][0]-2.0*ue[i][0]+ue[im1][0])) : ((void)0));

        (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C1) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C1)) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(1))) - tx2 * (
            (ue[ip1][1]*buf[ip1][1]+c2*(ue[ip1][4]-q[ip1]))-
            (ue[im1][1]*buf[im1][1]+c2*(ue[im1][4]-q[im1])))+
          xxcon1*(buf[ip1][1]-2.0*buf[i][1]+buf[im1][1])+
          dx2tx1*( ue[ip1][1]-2.0* ue[i][1]+ue[im1][1])) : ((void)0));

        (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C2) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C2)) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(2))) - tx2 * (
            ue[ip1][2]*buf[ip1][1]-ue[im1][2]*buf[im1][1])+
          xxcon2*(buf[ip1][2]-2.0*buf[i][2]+buf[im1][2])+
          dx3tx1*( ue[ip1][2]-2.0*ue[i][2] +ue[im1][2])) : ((void)0));

        (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C3) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C3)) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(3))) - tx2*(
            ue[ip1][3]*buf[ip1][1]-ue[im1][3]*buf[im1][1])+
          xxcon2*(buf[ip1][3]-2.0*buf[i][3]+buf[im1][3])+
          dx4tx1*( ue[ip1][3]-2.0* ue[i][3]+ ue[im1][3])) : ((void)0));

        (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(4)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(4))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(4))) - tx2*(
            buf[ip1][1]*(c1*ue[ip1][4]-c2*q[ip1])-
            buf[im1][1]*(c1*ue[im1][4]-c2*q[im1]))+
          0.5*xxcon3*(buf[ip1][0]-2.0*buf[i][0]+
              buf[im1][0])+
          xxcon4*(cuf[ip1]-2.0*cuf[i]+cuf[im1])+
          xxcon5*(buf[ip1][4]-2.0*buf[i][4]+buf[im1][4])+
          dx5tx1*( ue[ip1][4]-2.0* ue[i][4]+ ue[im1][4])) : ((void)0));
      //}

		for(m = 0; m < 5; m++){
		  if(i == 1){
			(dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) - dssp *
			(5.0*ue[i][m] - 4.0*ue[i+1][m] +ue[i+2][m])) : ((void)0));
		  }
		  if(i == 2){
			(dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) - dssp *
			(-4.0*ue[i-1][m] + 6.0*ue[i][m] -
            4.0*ue[i+1][m] +     ue[i+2][m])) : ((void)0)); 
		  }
		  if(i == grid_points[0]-3){
			(dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) - dssp *
			(ue[i-2][m] - 4.0*ue[i-1][m] +
			6.0*ue[i][m] - 4.0*ue[i+1][m])) : ((void)0)); 
		  }
		  if(i == grid_points[0]-2){
			 (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) - dssp *
			 (ue[i-2][m] - 4.0*ue[i-1][m] + 5.0*ue[i][m])) : ((void)0)); 
		  }
		  if((i >= 3)&&(i <= grid_points[0]-4)){
			(dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) - dssp*
            (ue[i-2][m] - 4.0*ue[i-1][m] +
             6.0*ue[i][m] - 4.0*ue[i+1][m] + ue[i+2][m])) : ((void)0));
		  }
		}
	  }
    }
  }
  
  
  
  
  
  
  
  
  //---------------------------------------------------------------------
  // eta-direction flux differences             
  //---------------------------------------------------------------------
  /*for (k = 1; k <= grid_points[2]-2; k++) {
    zeta = (double)(k) * dnzm1;
    for (i = 1; i <= grid_points[0]-2; i++) {
      xi = (double)(i) * dnxm1;

      for (j = 0; j <= grid_points[1]-1; j++) {
        eta = (double)(j) * dnym1;

        exact_solution(xi, eta, zeta, dtemp);
        for (m = 0; m < 5; m++) {
          ue[j][m] = dtemp[m];
        }

        dtpp = 1.0/dtemp[0];

        for (m = 1; m < 5; m++) {
          buf[j][m] = dtpp * dtemp[m];
        }

        cuf[j]    = buf[j][2] * buf[j][2];
        buf[j][0] = cuf[j] + buf[j][1] * buf[j][1] + buf[j][3] * buf[j][3];
        q[j] = 0.5*(buf[j][1]*ue[j][1] + buf[j][2]*ue[j][2] +
                    buf[j][3]*ue[j][3]);
      }

      for (j = 1; j <= grid_points[1]-2; j++) {
        jm1 = j-1;
        jp1 = j+1;

        forcing[k][j][i][0] = forcing[k][j][i][0] -
          ty2*( ue[jp1][2]-ue[jm1][2] )+
          dy1ty1*(ue[jp1][0]-2.0*ue[j][0]+ue[jm1][0]);

        forcing[k][j][i][1] = forcing[k][j][i][1] - ty2*(
            ue[jp1][1]*buf[jp1][2]-ue[jm1][1]*buf[jm1][2])+
          yycon2*(buf[jp1][1]-2.0*buf[j][1]+buf[jm1][1])+
          dy2ty1*( ue[jp1][1]-2.0* ue[j][1]+ ue[jm1][1]);

        forcing[k][j][i][2] = forcing[k][j][i][2] - ty2*(
            (ue[jp1][2]*buf[jp1][2]+c2*(ue[jp1][4]-q[jp1]))-
            (ue[jm1][2]*buf[jm1][2]+c2*(ue[jm1][4]-q[jm1])))+
          yycon1*(buf[jp1][2]-2.0*buf[j][2]+buf[jm1][2])+
          dy3ty1*( ue[jp1][2]-2.0*ue[j][2] +ue[jm1][2]);

        forcing[k][j][i][3] = forcing[k][j][i][3] - ty2*(
            ue[jp1][3]*buf[jp1][2]-ue[jm1][3]*buf[jm1][2])+
          yycon2*(buf[jp1][3]-2.0*buf[j][3]+buf[jm1][3])+
          dy4ty1*( ue[jp1][3]-2.0*ue[j][3]+ ue[jm1][3]);

        forcing[k][j][i][4] = forcing[k][j][i][4] - ty2*(
            buf[jp1][2]*(c1*ue[jp1][4]-c2*q[jp1])-
            buf[jm1][2]*(c1*ue[jm1][4]-c2*q[jm1]))+
          0.5*yycon3*(buf[jp1][0]-2.0*buf[j][0]+
              buf[jm1][0])+
          yycon4*(cuf[jp1]-2.0*cuf[j]+cuf[jm1])+
          yycon5*(buf[jp1][4]-2.0*buf[j][4]+buf[jm1][4])+
          dy5ty1*(ue[jp1][4]-2.0*ue[j][4]+ue[jm1][4]);
      }

      //---------------------------------------------------------------------
      // Fourth-order dissipation                      
      //---------------------------------------------------------------------
      for (m = 0; m < 5; m++) {
        j = 1;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (5.0*ue[j][m] - 4.0*ue[j+1][m] +ue[j+2][m]);
        j = 2;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (-4.0*ue[j-1][m] + 6.0*ue[j][m] -
           4.0*ue[j+1][m] +       ue[j+2][m]);
      }

      for (j = 3; j <= grid_points[1]-4; j++) {
        for (m = 0; m < 5; m++) {
          forcing[k][j][i][m] = forcing[k][j][i][m] - dssp*
            (ue[j-2][m] - 4.0*ue[j-1][m] +
             6.0*ue[j][m] - 4.0*ue[j+1][m] + ue[j+2][m]);
        }
      }

      for (m = 0; m < 5; m++) {
        j = grid_points[1]-3;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (ue[j-2][m] - 4.0*ue[j-1][m] +
           6.0*ue[j][m] - 4.0*ue[j+1][m]);
        j = grid_points[1]-2;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (ue[j-2][m] - 4.0*ue[j-1][m] + 5.0*ue[j][m]);
      }
    }
  }
  */
  
  for (k = 1; k <= grid_points[2]-2; k++) {   
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 1; j <= grid_points[1]-2; j++) {
		xi = (double)(i) * dnxm1;
	    zeta = (double)(k) * dnzm1; 
		
		if(j == 1){
			a = -1;
			b = 2;
		}
		if(j == 2){
			a = -2;
			b = 2;
		}
		if((j >= 3)&&(j <= grid_points[1]-4)){
			a = -2;
			b = 2;
		}
		if(j == grid_points[1]-3){
			a = -2;
			b = 2;
		}
		if(j == grid_points[1]-2){
			a = -2;
			b = 1;
		}
		for(z = a; z <= b; z++){
			eta = (double)(j+z) * dnym1;

			exact_solution(xi, eta, zeta, dtemp);
			for (m = 0; m < 5; m++) {
			  ue[j + z][m] = dtemp[m];
			}

			dtpp = 1.0/dtemp[0];

			for (m = 1; m < 5; m++) {
			  buf[j + z][m] = dtpp * dtemp[m];
			}

			cuf[j + z]    = buf[j + z][2] * buf[j + z][2];
			buf[j + z][0] = cuf[j + z] + buf[j + z][1] * buf[j + z][1] + buf[j + z][3] * buf[j + z][3];
			q[j + z] = 0.5*(buf[j + z][1]*ue[j + z][1] + buf[j + z][2]*ue[j + z][2] +
						buf[j + z][3]*ue[j + z][3]);
        }

      //for (j = 1; j <= grid_points[1]-2; j++) {
        jm1 = j-1;
        jp1 = j+1;

        (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C0) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C0)) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(0))) -
          ty2*( ue[jp1][2]-ue[jm1][2] )+
          dy1ty1*(ue[jp1][0]-2.0*ue[j][0]+ue[jm1][0])) : ((void)0));

        (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C1) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C1)) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(1))) - ty2*(
            ue[jp1][1]*buf[jp1][2]-ue[jm1][1]*buf[jm1][2])+
          yycon2*(buf[jp1][1]-2.0*buf[j][1]+buf[jm1][1])+
          dy2ty1*( ue[jp1][1]-2.0* ue[j][1]+ ue[jm1][1])) : ((void)0));

        (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C2) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C2)) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(2))) - ty2*(
            (ue[jp1][2]*buf[jp1][2]+c2*(ue[jp1][4]-q[jp1]))-
            (ue[jm1][2]*buf[jm1][2]+c2*(ue[jm1][4]-q[jm1])))+
          yycon1*(buf[jp1][2]-2.0*buf[j][2]+buf[jm1][2])+
          dy3ty1*( ue[jp1][2]-2.0*ue[j][2] +ue[jm1][2])) : ((void)0));

        (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C3) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C3)) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(3))) - ty2*(
            ue[jp1][3]*buf[jp1][2]-ue[jm1][3]*buf[jm1][2])+
          yycon2*(buf[jp1][3]-2.0*buf[j][3]+buf[jm1][3])+
          dy4ty1*( ue[jp1][3]-2.0*ue[j][3]+ ue[jm1][3])) : ((void)0));

        (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(4)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(4))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(4))) - ty2*(
            buf[jp1][2]*(c1*ue[jp1][4]-c2*q[jp1])-
            buf[jm1][2]*(c1*ue[jm1][4]-c2*q[jm1]))+
          0.5*yycon3*(buf[jp1][0]-2.0*buf[j][0]+
              buf[jm1][0])+
          yycon4*(cuf[jp1]-2.0*cuf[j]+cuf[jm1])+
          yycon5*(buf[jp1][4]-2.0*buf[j][4]+buf[jm1][4])+
          dy5ty1*(ue[jp1][4]-2.0*ue[j][4]+ue[jm1][4])) : ((void)0));
      //}

            
	    for(m = 0; m < 5; m++){
			if(j == 1){
				(dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) - dssp *
				(5.0*ue[j][m] - 4.0*ue[j+1][m] +ue[j+2][m])) : ((void)0));
			}
			if(j == 2){
				(dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) - dssp *
				(-4.0*ue[j-1][m] + 6.0*ue[j][m] -
				4.0*ue[j+1][m] +       ue[j+2][m])) : ((void)0));
			}
			if(j == grid_points[1]-3){
				(dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) - dssp *
				(ue[j-2][m] - 4.0*ue[j-1][m] +
				6.0*ue[j][m] - 4.0*ue[j+1][m])) : ((void)0));
			}
			if(j == grid_points[1]-2){
				(dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) - dssp *
				(ue[j-2][m] - 4.0*ue[j-1][m] + 5.0*ue[j][m])) : ((void)0));
			}
			if((j >= 3)&&(j <= grid_points[1]-4)){
				(dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) - dssp*
				(ue[j-2][m] - 4.0*ue[j-1][m] +
				6.0*ue[j][m] - 4.0*ue[j+1][m] + ue[j+2][m])) : ((void)0));
			}	
		}
	  
	  
	  
	  
	  
	  }
    }
  }
  
  
  
  //---------------------------------------------------------------------
  // zeta-direction flux differences                      
  //---------------------------------------------------------------------
  /*for (j = 1; j <= grid_points[1]-2; j++) {
    eta = (double)(j) * dnym1;
    for (i = 1; i <= grid_points[0]-2; i++) {
      xi = (double)(i) * dnxm1;

      for (k = 0; k <= grid_points[2]-1; k++) {
        zeta = (double)(k) * dnzm1;

        exact_solution(xi, eta, zeta, dtemp);
        for (m = 0; m < 5; m++) {
          ue[k][m] = dtemp[m];
        }

        dtpp = 1.0/dtemp[0];

        for (m = 1; m < 5; m++) {
          buf[k][m] = dtpp * dtemp[m];
        }

        cuf[k]    = buf[k][3] * buf[k][3];
        buf[k][0] = cuf[k] + buf[k][1] * buf[k][1] + buf[k][2] * buf[k][2];
        q[k] = 0.5*(buf[k][1]*ue[k][1] + buf[k][2]*ue[k][2] +
                    buf[k][3]*ue[k][3]);
      }

      for (k = 1; k <= grid_points[2]-2; k++) {
        km1 = k-1;
        kp1 = k+1;

        forcing[k][j][i][0] = forcing[k][j][i][0] -
          tz2*( ue[kp1][3]-ue[km1][3] )+
          dz1tz1*(ue[kp1][0]-2.0*ue[k][0]+ue[km1][0]);

        forcing[k][j][i][1] = forcing[k][j][i][1] - tz2 * (
            ue[kp1][1]*buf[kp1][3]-ue[km1][1]*buf[km1][3])+
          zzcon2*(buf[kp1][1]-2.0*buf[k][1]+buf[km1][1])+
          dz2tz1*( ue[kp1][1]-2.0* ue[k][1]+ ue[km1][1]);

        forcing[k][j][i][2] = forcing[k][j][i][2] - tz2 * (
            ue[kp1][2]*buf[kp1][3]-ue[km1][2]*buf[km1][3])+
          zzcon2*(buf[kp1][2]-2.0*buf[k][2]+buf[km1][2])+
          dz3tz1*(ue[kp1][2]-2.0*ue[k][2]+ue[km1][2]);

        forcing[k][j][i][3] = forcing[k][j][i][3] - tz2 * (
            (ue[kp1][3]*buf[kp1][3]+c2*(ue[kp1][4]-q[kp1]))-
            (ue[km1][3]*buf[km1][3]+c2*(ue[km1][4]-q[km1])))+
          zzcon1*(buf[kp1][3]-2.0*buf[k][3]+buf[km1][3])+
          dz4tz1*( ue[kp1][3]-2.0*ue[k][3] +ue[km1][3]);

        forcing[k][j][i][4] = forcing[k][j][i][4] - tz2 * (
            buf[kp1][3]*(c1*ue[kp1][4]-c2*q[kp1])-
            buf[km1][3]*(c1*ue[km1][4]-c2*q[km1]))+
          0.5*zzcon3*(buf[kp1][0]-2.0*buf[k][0]
              +buf[km1][0])+
          zzcon4*(cuf[kp1]-2.0*cuf[k]+cuf[km1])+
          zzcon5*(buf[kp1][4]-2.0*buf[k][4]+buf[km1][4])+
          dz5tz1*( ue[kp1][4]-2.0*ue[k][4]+ ue[km1][4]);
      }

      //---------------------------------------------------------------------
      // Fourth-order dissipation                        
      //---------------------------------------------------------------------
      for (m = 0; m < 5; m++) {
        k = 1;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (5.0*ue[k][m] - 4.0*ue[k+1][m] +ue[k+2][m]);
        k = 2;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (-4.0*ue[k-1][m] + 6.0*ue[k][m] -
           4.0*ue[k+1][m] +       ue[k+2][m]);
      }

      for (k = 3; k <= grid_points[2]-4; k++) {
        for (m = 0; m < 5; m++) {
          forcing[k][j][i][m] = forcing[k][j][i][m] - dssp*
            (ue[k-2][m] - 4.0*ue[k-1][m] +
             6.0*ue[k][m] - 4.0*ue[k+1][m] + ue[k+2][m]);
        }
      }

      for (m = 0; m < 5; m++) {
        k = grid_points[2]-3;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (ue[k-2][m] - 4.0*ue[k-1][m] +
           6.0*ue[k][m] - 4.0*ue[k+1][m]);
        k = grid_points[2]-2;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (ue[k-2][m] - 4.0*ue[k-1][m] + 5.0*ue[k][m]);
      }

    }
  }
  */
  
  for (j = 1; j <= grid_points[1]-2; j++) {   
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (k = 1; k <= grid_points[2]-2; k++) {
        xi = (double)(i) * dnxm1;
		eta = (double)(j) * dnym1;
		
		if(k == 1){
			a = -1;
			b = 2;
		}
		if(k == 2){
			a = -2;
			b = 2;
		}
		if((k >= 3)&&(k <= grid_points[1]-4)){
			a = -2;
			b = 2;
		}
		if(k == grid_points[1]-3){
			a = -2;
			b = 2;
		}
		if(k == grid_points[1]-2){
			a = -2;
			b = 1;
		}
		for(z = a; z <= b; z++){
			zeta = (double)(k + z) * dnzm1;

			exact_solution(xi, eta, zeta, dtemp);
			for (m = 0; m < 5; m++) {
			  ue[k + z][m] = dtemp[m];
			}

			dtpp = 1.0/dtemp[0];

			for (m = 1; m < 5; m++) {
			  buf[k + z][m] = dtpp * dtemp[m];
			}

			cuf[k + z]    = buf[k + z][3] * buf[k + z][3];
			buf[k + z][0] = cuf[k + z] + buf[k + z][1] * buf[k + z][1] + buf[k + z][2] * buf[k + z][2];
			q[k + z] = 0.5*(buf[k + z][1]*ue[k + z][1] + buf[k + z][2]*ue[k + z][2] +
						buf[k + z][3]*ue[k + z][3]);
		}

      //for (k = 1; k <= grid_points[2]-2; k++) {
        km1 = k-1;
        kp1 = k+1;

        (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C0) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C0)) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(0))) -
          tz2*( ue[kp1][3]-ue[km1][3] )+
          dz1tz1*(ue[kp1][0]-2.0*ue[k][0]+ue[km1][0])) : ((void)0));

        (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C1) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C1)) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(1))) - tz2 * (
            ue[kp1][1]*buf[kp1][3]-ue[km1][1]*buf[km1][3])+
          zzcon2*(buf[kp1][1]-2.0*buf[k][1]+buf[km1][1])+
          dz2tz1*( ue[kp1][1]-2.0* ue[k][1]+ ue[km1][1])) : ((void)0));

        (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C2) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C2)) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(2))) - tz2 * (
            ue[kp1][2]*buf[kp1][3]-ue[km1][2]*buf[km1][3])+
          zzcon2*(buf[kp1][2]-2.0*buf[k][2]+buf[km1][2])+
          dz3tz1*(ue[kp1][2]-2.0*ue[k][2]+ue[km1][2])) : ((void)0));

        (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C3) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C3)) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(3))) - tz2 * (
            (ue[kp1][3]*buf[kp1][3]+c2*(ue[kp1][4]-q[kp1]))-
            (ue[km1][3]*buf[km1][3]+c2*(ue[km1][4]-q[km1])))+
          zzcon1*(buf[kp1][3]-2.0*buf[k][3]+buf[km1][3])+
          dz4tz1*( ue[kp1][3]-2.0*ue[k][3] +ue[km1][3])) : ((void)0));

        (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(4)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(4))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(4))) - tz2 * (
            buf[kp1][3]*(c1*ue[kp1][4]-c2*q[kp1])-
            buf[km1][3]*(c1*ue[km1][4]-c2*q[km1]))+
          0.5*zzcon3*(buf[kp1][0]-2.0*buf[k][0]
              +buf[km1][0])+
          zzcon4*(cuf[kp1]-2.0*cuf[k]+cuf[km1])+
          zzcon5*(buf[kp1][4]-2.0*buf[k][4]+buf[km1][4])+
          dz5tz1*( ue[kp1][4]-2.0*ue[k][4]+ ue[km1][4])) : ((void)0));
      //}

      
        
		  
		  
		for(m = 0; m < 5; m++){
			if(k == 1){
				(dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) - dssp *
				(5.0*ue[k][m] - 4.0*ue[k+1][m] +ue[k+2][m])) : ((void)0));
			}
			if(k == 2){
				(dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) - dssp *
				(-4.0*ue[k-1][m] + 6.0*ue[k][m] -
				4.0*ue[k+1][m] +       ue[k+2][m])) : ((void)0));
			}
			if(k == grid_points[2]-3){
				(dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) - dssp *
				(ue[k-2][m] - 4.0*ue[k-1][m] +
				6.0*ue[k][m] - 4.0*ue[k+1][m])) : ((void)0));
			}
			if(k == grid_points[2]-2){
				(dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) - dssp *
				(ue[k-2][m] - 4.0*ue[k-1][m] + 5.0*ue[k][m])) : ((void)0));
			}
			if((k >= 3)&&(k <= grid_points[2]-4)){
				(dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) - dssp*
				(ue[k-2][m] - 4.0*ue[k-1][m] +
				6.0*ue[k][m] - 4.0*ue[k+1][m] + ue[k+2][m])) : ((void)0));
			}
		}
      }
    }
  }
  
  
  
  
  //---------------------------------------------------------------------
  // now change the sign of the forcing function, 
  //---------------------------------------------------------------------
  for (k = 1; k <= grid_points[2]-2; k++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (i = 1; i <= grid_points[0]-2; i++) {
        for (m = 0; m < 5; m++) {
          (dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)) ? ((*(double *)dvmh_get_own_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m))) = -1.0 * (*(double *)dvmh_get_element_addr_C(forcing, 4, DVM0C(k), DVM0C(j), DVM0C(i), DVM0C(m)))) : ((void)0));
        }
      }
    }
  }
}







