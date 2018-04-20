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

#include <math.h>
#include "header.h"

//---------------------------------------------------------------------
// this function computes the norm of the difference between the
// computed solution and the exact solution
//---------------------------------------------------------------------
void error_norm(double rms[5])
{
  int i, j, k, m, d;
  double xi, eta, zeta, u_exact[5], add;
  double r0, r1, r2, r3, r4;
  
  r0 = 0;
  r1 = 0;
  r2 = 0;
  r3 = 0;
  r4 = 0;

  for (m = 0; m < 5; m++) {
    rms[m] = 0.0;
  }

  #pragma dvm redistribute (u[block][block][block][])
  #pragma dvm region
  {
	  
  #pragma dvm parallel ([k][j][i] on u[k][j][i][]) reduction(sum(r1), sum(r2), sum(r3), sum(r4), sum(r0)), private(u_exact, xi, eta, zeta, m, add)
  for (k = 0; k <= grid_points[2]-1; k++) {    
    for (j = 0; j <= grid_points[1]-1; j++) {      
      for (i = 0; i <= grid_points[0]-1; i++) {
		zeta = (double)(k) * dnzm1;
		eta = (double)(j) * dnym1;
        xi = (double)(i) * dnxm1;
        exact_solution(xi, eta, zeta, u_exact);

        //for (m = 0; m < 5; m++) {
          //add = u[k][j][i][m]-u_exact[m];
          //rms[m] = rms[m] + add*add;
        //}
		add = u[k][j][i][0] - u_exact[0];
		r0 = r0 + add * add;
		add = u[k][j][i][1] - u_exact[1];
		r1 = r1 + add * add;
		add = u[k][j][i][2] - u_exact[2];
		r2 = r2 + add * add;
		add = u[k][j][i][3] - u_exact[3];
		r3 = r3 + add * add;
		add = u[k][j][i][4] - u_exact[4];
		r4 = r4 + add * add;
		
      }
    }
  }

  }//end region
  #pragma dvm redistribute(u[][][][])
  
  rms[0] = r0;
  rms[1] = r1;
  rms[2] = r2;
  rms[3] = r3;
  rms[4] = r4;
  
  for (m = 0; m < 5; m++) {
    for (d = 0; d < 3; d++) {
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    }
    rms[m] = sqrt(rms[m]);
  }
}


void rhs_norm(double rms[5])
{
  int i, j, k, d, m;
  double add;
  double r0, r1, r2, r3, r4;
  
  r0 = 0;
  r1 = 0;
  r2 = 0;
  r3 = 0;
  r4 = 0;

  for (m = 0; m < 5; m++) {
    rms[m] = 0.0;
  } 

  #pragma dvm redistribute (u[block][block][block][])
  #pragma dvm region
  {
	  
  #pragma dvm parallel ([k][j][i] on u[k][j][i][]) reduction(sum(r1), sum(r2), sum(r3), sum(r4), sum(r0)), private(add)

  for (k = 1; k <= grid_points[2]-2; k++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (i = 1; i <= grid_points[0]-2; i++) {
        //for (m = 0; m < 5; m++) {
          //add = rhs[k][j][i][m];
          //rms[m] = rms[m] + add*add;
        //}
		add = rhs[k][j][i][0];
		r0 = r0 + add * add;
		add = rhs[k][j][i][1];
		r1 = r1 + add * add;
		add = rhs[k][j][i][2];
		r2 = r2 + add * add;
		add = rhs[k][j][i][3];
		r3 = r3 + add * add;
		add = rhs[k][j][i][4];
		r4 = r4 + add * add;
		
      } 
    } 
  }

  }//end region
  #pragma dvm redistribute(u[][][][])
  
  
  rms[0] = r0;
  rms[1] = r1;
  rms[2] = r2;
  rms[3] = r3;
  rms[4] = r4;
  
  for (m = 0; m < 5; m++) {
    for (d = 0; d < 3; d++) {
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    } 
    rms[m] = sqrt(rms[m]);
  } 
}
