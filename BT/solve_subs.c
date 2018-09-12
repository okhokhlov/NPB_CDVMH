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

#include "header.h"
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

