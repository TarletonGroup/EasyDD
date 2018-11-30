/*
 Calculate nodal forces induced by a non-singular straight segment of dislocation on a linear rectangular surface element.

 Notations and details can be found in S. Queyreau, J. Marian, B.D. Wirth, A. Arsenlis, MSMSE, 22(3):035004, (2014).

 Translated from matlab code into C by Daniel Celis Garza.
 Edits: June 1, 2017.
*/
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

// Calculate atanh because MS Visual Studio is inferior to GCC.
#ifdef _WIN32
  double atanh( double r ){
    return 0.5 * (log(1+r) - log(1-r));
  }
#endif
#include "vector_utils.h"
#include "serial_forces_lin_rect.h"

// Testing purposes
int main(int argc, char **argv){
  double x1[3], x2[3], x3[3], x4[3], x5[3], x6[3], b[3], mu, nu, a, *fx[n_nodes], ftot[3];
  double *dln_node_arr[2], *se_node_arr[n_nodes];
  double *b_arr;
  double *fx_arr[n_nodes];
  double *ftot_arr;
  double a_sq, factor, one_m_nu;
  int idx1, idx2, n_se, n_dln;
  for(int i = 0; i < n_nodes; i++){
    fx[i] = malloc(3*sizeof(double));
  }

  // Read input.
  FILE * ptr_file;
  ptr_file = fopen(argv[1], "r");
  if (ptr_file == NULL){
    printf("File does not exist.\n");
    exit(EXIT_FAILURE);
  }

  fscanf(ptr_file, "%i", &n_dln);
  fscanf(ptr_file, "%i", &n_se);

  // Memory allocation
  b_arr = (double *) malloc(3 * n_dln * sizeof(double));
  dln_node_arr[0] = (double *) malloc(3 * n_dln * sizeof(double));
  dln_node_arr[1] = (double *) malloc(3 * n_dln * sizeof(double));
  for (int i = 0; i < n_nodes; i++){
    se_node_arr[i] = (double *) malloc(3 * n_se  * sizeof(double));
    fx[i] = (double *) malloc(3 * sizeof(double));
    fx_arr[i] = (double *) malloc(3 * n_se * sizeof(double));
  }
  ftot_arr = (double *) malloc(3 * n_se * sizeof(double));

  for (int i = 0; i < n_dln*3; i+=3){
    fscanf(ptr_file, "%lf %lf %lf", &dln_node_arr[0][i], &dln_node_arr[0][i+1], &dln_node_arr[0][i+2] );
    fscanf(ptr_file, "%lf %lf %lf", &dln_node_arr[1][i], &dln_node_arr[1][i+1], &dln_node_arr[1][i+2] );
  }
  for (int i = 0; i < n_se*3; i+=3){
    for (int j = 0; j < n_nodes; j++){
      fscanf(ptr_file, "%lf %lf %lf", &se_node_arr[j][i], &se_node_arr[j][i+1], &se_node_arr[j][i+2] );
    }
  }
  for (int i = 0; i < n_dln*3; i+=3){
    fscanf(ptr_file, "%lf %lf %lf", &b_arr[i], &b_arr[i+1], &b_arr[i+2] );
  }
  fscanf(ptr_file, "%lf", &mu );
  fscanf(ptr_file, "%lf", &nu );
  fscanf(ptr_file, "%lf", &a );
  fclose(ptr_file);

  a_sq     = a*a;
  one_m_nu = 1.-nu;
  factor   = 0.25*mu/pi/one_m_nu;

  idx1 = 0;
  for (int i = 0; i < n_se; i++){
    idx2 = 0;
    // Transfer rectangular element i's coordinates into x3--x6.
    for (int k = 0; k < 3; k++){
      x3[k] = se_node_arr[0][idx1+k];
      x4[k] = se_node_arr[1][idx1+k];
      x5[k] = se_node_arr[2][idx1+k];
      x6[k] = se_node_arr[3][idx1+k];
    }
    // Loop through the dislocation segments.
    for (int j = 0; j < n_dln; j++){
      // Transfer dislocation segment j's coordinates and burger's vector into x1--x2 and b
      for (int k = 0; k < 3; k++){
        x1[k] = dln_node_arr[0][idx2+k];
        x2[k] = dln_node_arr[1][idx2+k];
        b [k] = b_arr[idx2+k];
      }
      main_nodal_surface_force_linear_rectangle(x1,x2,x3,x4,x5,x6,b,mu,nu,a,a_sq,one_m_nu,factor,fx,ftot);
      for (int k = 0; k < 3; k++){
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
  #ifdef debug
    for (int i = 0; i < n_se; i++){
      printf("ftot_arr[%d] = [%f, %f, %f]\n", i, ftot_arr[3*i], ftot_arr[3*i+1], ftot_arr[3*i+2]);
    }
  #endif

  //main_nodal_surface_force_linear_rectangle(x1, x2, x3, x4, x5, x6, b, mu, nu, a, a_sq, one_m_nu, factor, fx, ftot);
  for(int i = 0; i < n_nodes; i++){free(fx[i]); free(se_node_arr[i]);}
  free(b_arr); free(dln_node_arr[0]); free(dln_node_arr[1]); free(ftot_arr);
  return 0;
}
