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
#include <mex.h>
#include "vector_utils.h"
#include "serial_forces_lin_rect.h"
#include "cuda_vector_map_utils.h"
#include "cuda_forces_lin_rect.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
  // Node arrays from MATLAB. To be mapped into x_se_arr and then passed to d_x_se_arr.
  double *dln_node_arr[2], *se_node_arr[n_nodes];
  // Variables for the special case where the line segment and surface element are parallel.
  double x1[3], x2[3], x3[3], x4[3], x5[3], x6[3], b[3], p[3], q[3], t[3], n[3], p_norm, q_norm;
  double *fx[n_nodes];
  double ftot[3];
  // Burger's vectors array from MATLAB. To be passed straight into d_b_arr.
  double *b_arr[1];
  // Material properties from MATLAB to be placed into shared memory in device.
  double mu, nu, a, a_sq, one_m_nu, factor;
  // Nodal force array (3 coordinates per node per SE, 3*n_nodes*n_se). To be inversely mapped to *fx_arr[n_nodes].
  double *x_fx_arr;
  // Nodal force array to be sent back to MATLAB.
  double *fx_arr[n_nodes];
  // Total force array on SE (3 coordinates per SE, 3*n_se) to be sent back to MATLAB.
  double *ftot_arr, *x_ftot_arr;
  // Maps of SE and DLN node arrays.
  double *x_se_arr, *x_dln_arr, *x_b_arr;
  // Device arrays.
  double *d_x_b_arr, *d_x_se_arr, *d_x_dln_arr, *d_fx_arr, *d_ftot_arr;
  int threads_per_block, blocks_per_grid, n_se, n_dln, para_scheme;
  //int debug = 1;
  //while(debug == 1){}
  cudaSetDevice(0);
  // Stagger cuda function calls to take advantage of asynchronous calls.
  // If memory becomes an issue, make copying x_dln_arr, x_se_arr and x_b_arr to the device a synchronous operation and free the pointers straight after.
  n_se = (int) mxGetScalar(prhs[10]);
  // Allocate and set forces to 0 in device.
  checkCudaErrors( cudaMalloc( (void **) &d_fx_arr  , 3 * n_se  * n_nodes * sizeof(double) ) );
  checkCudaErrors( cudaMalloc( (void **) &d_ftot_arr, 3 * n_se            * sizeof(double) ) );
  checkCudaErrors( cudaMemsetAsync(d_fx_arr  , 0.0, 3 * n_se * n_nodes * sizeof(double)) );
  checkCudaErrors( cudaMemsetAsync(d_ftot_arr, 0.0, 3 * n_se           * sizeof(double)) );
  // Execute host code while device sets force arrays to zero.
  n_dln = (int) mxGetScalar(prhs[11]);
  dln_node_arr[0] = (double *) mxGetPr(prhs[0]);
  dln_node_arr[1] = (double *) mxGetPr(prhs[1]);
  // Execute host code while copying values from host to device.
  se_node_arr[0] = (double *) mxGetPr(prhs[2]);
  se_node_arr[1] = (double *) mxGetPr(prhs[3]);
  se_node_arr[2] = (double *) mxGetPr(prhs[4]);
  se_node_arr[3] = (double *) mxGetPr(prhs[5]);
  b_arr[0] = (double *) mxGetPr(prhs[6]);
  para_scheme = (int) mxGetScalar(prhs[13]);
  // Map dislocation node arrays to 1D array for parallelisation.
  if(para_scheme == 1){
    x_dln_arr = element_host_device_map(dln_node_arr, n_dln, 2);
    // Allocate and copy values of dislocation nodes to device.
    checkCudaErrors( cudaMalloc     ( (void **) &d_x_dln_arr, 3 * n_dln * 2 * sizeof(double) ) );
    checkCudaErrors( cudaMemcpyAsync(d_x_dln_arr,  x_dln_arr, 3 * n_dln * 2 * sizeof(double), cudaMemcpyHostToDevice) );
    x_se_arr  = se_host_device_map(se_node_arr[0], se_node_arr[1], se_node_arr[2], se_node_arr[3], n_se);
    // Allocate and copy values of surface element nodes to device.
    checkCudaErrors( cudaMalloc     ( (void **) &d_x_se_arr, 3 * n_se * n_nodes * sizeof(double) ) );
    checkCudaErrors( cudaMemcpyAsync(d_x_se_arr,   x_se_arr, 3 * n_se * n_nodes * sizeof(double), cudaMemcpyHostToDevice) );
    // Map Burger's vector array to 1D array for parallelisation.
    x_b_arr  = element_host_device_map(b_arr, n_dln, 1);
    // Allocate and copy values of Burger's vectors to device.
    checkCudaErrors( cudaMalloc     ( (void **) &d_x_b_arr, 3 * n_dln * sizeof(double) ) );
    checkCudaErrors( cudaMemcpyAsync(d_x_b_arr,    x_b_arr, 3 * n_dln * sizeof(double), cudaMemcpyHostToDevice) );
  }
  else{
    x_dln_arr = dln_host_device_map(dln_node_arr[0], dln_node_arr[1], n_dln);
    // Allocate and copy values of dislocation nodes to device.
    checkCudaErrors( cudaMalloc     ( (void **) &d_x_dln_arr, 3 * n_dln * 2 * sizeof(double) ) );
    checkCudaErrors( cudaMemcpyAsync(d_x_dln_arr,  x_dln_arr, 3 * n_dln * 2 * sizeof(double), cudaMemcpyHostToDevice) );
    x_se_arr  = element_host_device_map(se_node_arr, n_se, n_nodes);
    // Allocate and copy values of surface element nodes to device.
    checkCudaErrors( cudaMalloc     ( (void **) &d_x_se_arr, 3 * n_se * n_nodes * sizeof(double) ) );
    checkCudaErrors( cudaMemcpyAsync(d_x_se_arr,   x_se_arr, 3 * n_se * n_nodes * sizeof(double), cudaMemcpyHostToDevice) );
    x_b_arr   = b_host_device_map(b_arr[0], n_dln);
    // Allocate and copy values of Burger's vectors to device.
    checkCudaErrors( cudaMalloc     ( (void **) &d_x_b_arr, 3 * n_dln * sizeof(double) ) );
    checkCudaErrors( cudaMemcpyAsync(d_x_b_arr,    x_b_arr, 3 * n_dln * sizeof(double), cudaMemcpyHostToDevice) );
  }

  // Execute host code while copying values from host to device.
  // Copy constant values to device.
  mu = mxGetScalar(prhs[7]);
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_mu, &mu, sizeof(mu)) );
  nu = mxGetScalar(prhs[8]);
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_nu, &nu, sizeof(nu)) );
  one_m_nu = 1.-nu;
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_one_m_nu, &one_m_nu, sizeof(one_m_nu)) );
  a  = mxGetScalar(prhs[9]);
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_a, &a, sizeof(a)) );
  a_sq     = a*a;
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_a_sq, &a_sq, sizeof(a_sq)) );
  factor   = 0.25*mu/pi/one_m_nu;
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_factor, &factor, sizeof(factor)) );
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_eps   , &eps     , sizeof(eps)) );
  // Link force arrays to MATLAB.
  plhs[0] = mxCreateDoubleMatrix(3 * n_se, 1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(3 * n_se, 1, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(3 * n_se, 1, mxREAL);
  plhs[3] = mxCreateDoubleMatrix(3 * n_se, 1, mxREAL);
  plhs[4] = mxCreateDoubleMatrix(3 * n_se, 1, mxREAL);
  fx_arr[0] = (double *) mxGetPr(plhs[0]);
  fx_arr[1] = (double *) mxGetPr(plhs[1]);
  fx_arr[2] = (double *) mxGetPr(plhs[2]);
  fx_arr[3] = (double *) mxGetPr(plhs[3]);
  ftot_arr  = (double *) mxGetPr(plhs[4]);
  threads_per_block = (int) mxGetScalar(prhs[12]);
  // CUDA
  if(para_scheme == 1){
    blocks_per_grid   = (n_dln + threads_per_block - 1) / threads_per_block;
    dln_cuda_nodal_surface_force_linear_rectangle<<<blocks_per_grid, threads_per_block>>>(d_x_dln_arr, d_x_se_arr, d_x_b_arr, d_fx_arr, d_ftot_arr, n_se, n_dln);
  }
  else{
    blocks_per_grid   = (n_se + threads_per_block - 1) / threads_per_block;
    se_cuda_nodal_surface_force_linear_rectangle<<<blocks_per_grid, threads_per_block>>>(d_x_dln_arr, d_x_se_arr, d_x_b_arr, d_fx_arr, d_ftot_arr, n_se, n_dln);
  }
  // Host code is executed asynchronously from the kernel execution.
  // Free all 1D arrays used to copy data to device.
  free(x_se_arr); free(x_dln_arr); free(x_b_arr);
  x_fx_arr   = (double *) malloc(3 * n_se * n_nodes * sizeof(double));
  // Special case, where dislocation line is parallel with surface element.
  // Default behaviour is we take special cases into account.
  // Add the custom compiler flag -Dsc at the end of the compilation
  // command when wanting to ignore with the special case.
  #ifndef sc
    int idx1, idx2;
    // Initialise forces.
    for (int i = 0; i < n_nodes; i++){
	  fx[i] = (double * ) malloc(3*sizeof(double));
      for (int j = 0; j < 3*n_se; j++){
        fx_arr[i][j] = 0.0;
      }
    }
    for (int i = 0; i < 3*n_se; i++){
      ftot_arr[i] = 0.0;
    }
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
          b [k] = b_arr[0][idx2+k];
        }
        init_vector(x1, x2, 3, t);
        init_vector2(x3, x4, 3, p, &p_norm);
        init_vector2(x3, x5, 3, q, &q_norm);
        cross_product(p, q, n);
        normalise_vector(n, 3, n);
        if (dot_product(t, n, 3) == 0.0){
          nodal_surface_force_linear_rectangle_special(x1, x2, x3, x4, x5, x6, b, t, p, q, n, p_norm, q_norm, mu, nu, a, a_sq, one_m_nu, factor/p_norm/q_norm, fx, ftot);
          // Add the force contributions for segment j to the surface element i.
          for (int k = 0; k < 3; k++){
            fx_arr[0][idx1+k] += fx[0][k];
            fx_arr[1][idx1+k] += fx[1][k];
            fx_arr[2][idx1+k] += fx[2][k];
            fx_arr[3][idx1+k] += fx[3][k];
            ftot_arr [idx1+k] += ftot[k];
          }
        }
        idx2 += 3;
      }
      idx1 += 3;
    }
    x_ftot_arr = (double *) malloc(3 * n_se * sizeof(double));
    // Synchronously copy forces from device to host.
    checkCudaErrors( cudaMemcpy(x_fx_arr, d_fx_arr, 3 * n_se * n_nodes * sizeof(double), cudaMemcpyDeviceToHost) );
    if(para_scheme == 1){
      // Map 1D device array to 2D array for MATLAB.
      dln_add_fx_device_host_map(x_fx_arr, fx_arr, n_se, n_nodes);
      // Synchronously copy forces from device to host.
      checkCudaErrors( cudaMemcpy(x_ftot_arr, d_ftot_arr, 3 * n_se * sizeof(double), cudaMemcpyDeviceToHost) );
      for (int i = 0; i < 3*n_se; i++){
        ftot_arr[i] += x_ftot_arr[i];
      }
    }
    else{
      // Map 1D device array to 2D array for MATLAB.
      add_fx_device_host_map(x_fx_arr, fx_arr, n_se, n_nodes);
      // Synchronously copy forces from device to host.
      checkCudaErrors( cudaMemcpy(x_ftot_arr, d_ftot_arr, 3 * n_se * sizeof(double), cudaMemcpyDeviceToHost) );
      add_ftot_device_host_map(x_ftot_arr, ftot_arr, n_se);
    }
    free(x_ftot_arr);
  #endif

  // This code snippet ignores the special case where dislocation lines are parallel to the surface element.
  // It is only compiled if the flag -Dsc for the special case is present.
  #ifdef sc
  // Synchronously copy forces from device to host.
  checkCudaErrors( cudaMemcpy(x_fx_arr, d_fx_arr, 3 * n_se * n_nodes * sizeof(double), cudaMemcpyDeviceToHost) );
    if(para_scheme == 1){
      // Map 1D device array to 2D array for MATLAB.
      dln_fx_device_host_map(x_fx_arr, fx_arr, n_se, n_nodes);
      checkCudaErrors( cudaMemcpy(ftot_arr, d_ftot_arr, 3 * n_se * sizeof(double), cudaMemcpyDeviceToHost) );
    }
    else{
      x_ftot_arr = (double *) malloc(3 * n_se * sizeof(double));
      fx_device_host_map(x_fx_arr, fx_arr, n_se, n_nodes);
      checkCudaErrors( cudaMemcpy(x_ftot_arr, d_ftot_arr, 3 * n_se * sizeof(double), cudaMemcpyDeviceToHost) );
      ftot_device_host_map(x_ftot_arr, ftot_arr, n_se);
      free(x_ftot_arr);
    }
  #endif
	for (int i = 0; i < n_nodes; i++){
		free(fx[i]);
	}
   free(x_fx_arr);
  // CUDA exit.
  cudaDeviceReset();
}
