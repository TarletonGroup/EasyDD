/*
 Calculate nodal forces induced by a non-singular straight segment of dislocation on a linear rectangular surface element.

 Notations and details can be found in S. Queyreau, J. Marian, B.D. Wirth, A. Arsenlis, MSMSE, 22(3):035004, (2014).

 Translated from matlab code into C by Daniel Celis Garza.
 Parallelised by Daniel Celis Garza
 Edits: June 1, 2017.
 Parallelisaion: August 31, 2017
*/
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include "helper_cuda.h"
#ifdef debug3
  #include <time.h>
#endif
#include "vector_utils.h"
#include "serial_forces_lin_rect.h"
#include "cuda_vector_map_utils.h"
#include "cuda_forces_lin_rect.h"
/*
#ifdef _WIN32
  double atanh( double r ){
    return 0.5 * (log(1+r) - log(1-r));
  }
#endif
*/

/* Work out best parallelisation.
  // Maximum number of threads per block.
  //cudaDeviceProp deviceProp;
  //cudaGetDeviceProperties(&deviceProp, 1);
  // Number of blocks launched.
  // deviceProp.maxGridSize[0], deviceProp.maxGridSize[1], deviceProp.maxGridSize[2]
  // Size of each block launched
  // deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1], deviceProp.maxThreadsDim[2]
*/
// Testing purposes
int main(int argc, char **argv){
  int n_se, n_dln, threads_per_block;
  FILE * ptr_file;
  printf("Over SEs\n");
  ptr_file = fopen(argv[1], "r");
  if (ptr_file == NULL){
    printf("File does not exist.\n");
    exit(EXIT_FAILURE);
  }
  fscanf(ptr_file, "%i", &n_dln);
  fscanf(ptr_file, "%i", &n_se);
  fscanf(ptr_file, "%i", &threads_per_block);
  fclose(ptr_file);
  main_se_cuda_nodal_surface_force_linear_rectangle(n_se, n_dln, threads_per_block, argv);
  printf("\n");

  printf("Over DLNs\n");
  ptr_file = fopen(argv[1], "r");
  if (ptr_file == NULL){
    printf("File does not exist.\n");
    exit(EXIT_FAILURE);
  }
  fscanf(ptr_file, "%i", &n_dln);
  fscanf(ptr_file, "%i", &n_se);
  fscanf(ptr_file, "%i", &threads_per_block);
  fclose(ptr_file);
  main_dln_cuda_nodal_surface_force_linear_rectangle(n_se, n_dln, threads_per_block, argv);
  printf("\n");

  return 0;
}
