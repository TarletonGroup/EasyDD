/*
 Calculate nodal forces induced by a non-singular straight segment of dislocation on a linear rectangular surface element.

 Notations and details can be found in S. Queyreau, J. Marian, B.D. Wirth, A. Arsenlis, MSMSE, 22(3):035004, (2014).

 Translated from matlab code into C by Daniel Celis Garza.
 Edits: June 12, 2017.
*/
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mex.h>

#include "vector_utils.h"
#include "serial_forces_lin_rect.h"

// Calculate atanh because MS Visual Studio is inferior to GCC.
#ifdef _WIN32
  double atanh( double r ){
    return 0.5 * (log(1+r) - log(1-r));
  }
#endif

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
  // We bundle fx into an array op pointers because it makes it easier to make loops of many operations. We don't do it with the points x1--x6 because we don't gain anything from it aside from shorter argmument lists for functions, makes things harder to read, and would make it easy to mistakenly modify x1 and x2 as a side effect of having to rotate them if the line vector t is parallel to the surface element.
  double *x1_arr, *x2_arr, *x3_arr, *x4_arr, *x5_arr, *x6_arr;
  double *b_arr;
  double mu,nu,a;
  double a_sq, one_m_nu, factor;
  double *fx_arr[4];
  double *ftot_arr;
  int num_surface_elements, num_dislocation_segs;
  double x1[3], x2[3], x3[3], x4[3], x5[3], x6[3];
  double b[3];
  double *fx[4];
  double ftot[3];
  int i, j, k, idx1, idx2;
  //int debug = 1;
  //do {} while( debug == 1 );

  x1_arr = (double *) mxGetPr(prhs[0]);
  x2_arr = (double *) mxGetPr(prhs[1]);
  x3_arr = (double *) mxGetPr(prhs[2]);
  x4_arr = (double *) mxGetPr(prhs[3]);
  x5_arr = (double *) mxGetPr(prhs[4]);
  x6_arr = (double *) mxGetPr(prhs[5]);
  b_arr  = (double *) mxGetPr(prhs[6]);
  mu = mxGetScalar(prhs[7]);
  nu = mxGetScalar(prhs[8]);
  a  = mxGetScalar(prhs[9]);
  num_surface_elements = mxGetScalar(prhs[10]);
  num_dislocation_segs = mxGetScalar(prhs[11]);

  plhs[0] = mxCreateDoubleMatrix(3*num_surface_elements,1,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(3*num_surface_elements,1,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(3*num_surface_elements,1,mxREAL);
  plhs[3] = mxCreateDoubleMatrix(3*num_surface_elements,1,mxREAL);
  plhs[4] = mxCreateDoubleMatrix(3*num_surface_elements,1,mxREAL);

  fx_arr[0] = (double *) mxGetPr(plhs[0]);
  fx_arr[1] = (double *) mxGetPr(plhs[1]);
  fx_arr[2] = (double *) mxGetPr(plhs[2]);
  fx_arr[3] = (double *) mxGetPr(plhs[3]);
  ftot_arr  = (double *) mxGetPr(plhs[4]);
  for (i = 0; i < 4; i++){
      fx[i] = malloc(3*sizeof(double));
  }
  a_sq     = a*a;
  one_m_nu = 1.-nu;
  factor   = 0.25*mu/pi/one_m_nu;
  // Naive case, we repeat many operations.
  // TO DO: Change this from the naive case to something more efficient.
  // Loop through the number of surface elements.
  idx1 = 0;
  for (i = 0; i < num_surface_elements; i++){
    idx2 = 0;
    // Transfer rectangular element i's coordinates into x3--x6.
    for (k = 0; k < 3; k++){
      x3[k] = x3_arr[idx1+k];
      x4[k] = x4_arr[idx1+k];
      x5[k] = x5_arr[idx1+k];
      x6[k] = x6_arr[idx1+k];
    }
    // Loop through the dislocation segments.
    for (j = 0; j < num_dislocation_segs; j++){
      // Transfer dislocation segment j's coordinates and burger's vector into x1--x2 and b
      for (k = 0; k < 3; k++){
        x1[k] = x1_arr[idx2+k];
        x2[k] = x2_arr[idx2+k];
        b [k] = b_arr [idx2+k];
      }
      main_nodal_surface_force_linear_rectangle(x1,x2,x3,x4,x5,x6,b,mu,nu,a,a_sq,one_m_nu,factor,fx,ftot);
      // Add the force contributions for segment j to the surface element i.
      for (k = 0; k < 3; k++){
        fx_arr[0][idx1+k] += fx[0][k];
        fx_arr[1][idx1+k] += fx[1][k];
        fx_arr[2][idx1+k] += fx[2][k];
        fx_arr[3][idx1+k] += fx[3][k];
        ftot_arr [idx1+k] += ftot[k];
      }
      idx2 += 3;
    }
    idx1 += 3;
  }
  for (i = 0; i < 4; i++){
      free(fx[i]);
  }
}

/*
// Testing purposes
int main(void){
  double x1[3], x2[3], x3[3], x4[3], x5[3], x6[3], b[3], mu, nu, a, *nodal_force[4], total_force[3];
  double a_sq, one_m_nu, factor;
  for(int i = 0; i < 4; i++){
    nodal_force[i] = malloc(3*sizeof(double));
  }
  x1[0] = 0.5; x1[1] = 0.; x1[2] = 0.5;
  x2[0] = 0.5; x2[1] = 1.; x2[2] = 1.5;
  x3[0] = 0.; x3[1] = 0.; x3[2] = 0.;
  x4[0] = 1.; x4[1] = 0.; x4[2] = 0.;
  x5[0] = 0.; x5[1] = 1.; x5[2] = 1.;
  x6[0] = 1.; x6[1] = 1.; x6[2] = 1.;
  b[0] = -3./10.; b[1] = 5./10.; b[2] = 4./10.;
  mu = 0.6;
  nu = 0.3;
  a = 0.01;
  a_sq = a*a;
  one_m_nu = 1.-nu;
  factor   = 0.25*mu/pi/one_m_nu;
  main_nodal_surface_force_linear_rectangle(x1, x2, x3, x4, x5, x6, b, mu, nu, a, a_sq, one_m_nu, factor, nodal_force, total_force);
}
*/
