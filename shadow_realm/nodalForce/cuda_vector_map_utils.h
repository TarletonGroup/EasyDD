// double precision atomicAdd_dbl
#if __CUDA_ARCH__ < 600
  __device__ double atomicAdd_dbl(double* address, double val)
  {
      unsigned long long int* address_as_ull =
                                (unsigned long long int*)address;
      unsigned long long int old = *address_as_ull, assumed;

      do {
          assumed = old;
          old = atomicCAS(address_as_ull, assumed,
                          __double_as_longlong(val +
                                 __longlong_as_double(assumed)));

      // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
      } while (assumed != old);

      return __longlong_as_double(old);
  }
#endif

double *element_host_device_map(double *i_node_arr[], int i_n_elem_scope, int i_n_nodes){
  /*
    Maps E elements with N nodes from the host to the global device array.
    This ensures coalesced memory accesses when parallelising over elements (be it dislocations or surface elements or both via dynamic parallelism).
    i_node_arr[n][3e + 0:2] = [x_en, y_en, z_en]

      to

    o_g_elem_arr =
      [x_00, x_10, ..., x_(E-1)0, x_01, x_11, ..., x_(E-1)1, ..., x_0N, x_1(N-1), ..., x_(E-1)(N-1),
       y_00, y_10, ..., y_(E-1)0, y_01, y_11, ..., y_(E-1)1, ..., y_0N, y_1(N-1), ..., y_(E-1)(N-1),
       z_00, z_10, ..., z_(E-1)0, z_01, z_11, ..., z_(E-1)1, ..., z_0N, z_1(N-1), ..., z_(E-1)(N-1)]
    i_n_elem_scope  = E
    i_n_nodes       = N
    n_elem_n_nodes  = E*N
    idxi = input (initial) index
    idxt = temp index
    idxf = output (final) index
  */
  int const n_elem_n_nodes = i_n_elem_scope*i_n_nodes;
  double *o_g_elem_arr;
  int idxt = 0, idxf = 0, idxi = 0;
  // Allocate a 1D output array of length 3*E*N.
  o_g_elem_arr = (double *) malloc(3 * n_elem_n_nodes * sizeof(double));
  // Loop over the nodes of the element.
  for (int i = 0; i < i_n_nodes; i++){
    // Reset the output index to point at the i'th node of the first coordinate of the first element.
    idxf = idxt;
    // Loop over coordinates.
    for (int j = 0; j < 3; j++){
      // Set the input index to point at the j'th coordinate of the first element of the i'th node.
      idxi = j;
      // Loop over elements.
      for (int k = 0; k < i_n_elem_scope; k++){
        o_g_elem_arr[idxf + k] = i_node_arr[i][idxi];
        // Advance the input index to point at the j'th coordinate of the (k+1)'th element of the i'th node.
        idxi += 3;
      }
      // Advance the output index to point at the i'th node of the j'th coordinate of the first element.
      idxf += n_elem_n_nodes;
    }
    // Advance the temporary index to point at the (i+1)'th node of the first coordinate of the first element.
    idxt += i_n_elem_scope;
  }
  return o_g_elem_arr;
}

double *se_host_device_map(double *i_x3_arr, double *i_x4_arr, double *i_x5_arr, double *i_x6_arr, int i_n_se){
  /*
    Maps E surface elements with 4 nodes.
    This is for looping through in a parallelisation over dislocation line segments.
    i_xn_arr[3e + 0:2] = [x_en, y_en, z_en]

      to

    o_g_dln_arr =
      [x_00, y_00, z_00, x_01, y_01, z_01, ..., x_0(N-1)    , y_0(N-1)    , z_0(N-1)
       x_10, y_10, z_10, x_11, y_11, z_11, ..., x_1(N-1)    , y_1(N-1)    , z_1(N-1), ...,
       x_(E-1)0, y_(E-1)0, z_(E-1)0      , ..., x_(E-1)(N-1), y_(E-1)(N-1), z_(E-1)(N-1)]
    i_se = E
    n_se_nodes = 4*E
    idxi = input (initial) index
    idxf = output (final) index
  */
  int const n_se_nodes = i_n_se*4;
  double *o_g_se_arr;
  int idxi = 0, idxf = 0;
  // Allocate a 1D output array of length 3*2*E.
  o_g_se_arr = (double *) malloc(3 * n_se_nodes * sizeof(double));
  // Loop over dislocation line segments.
  for (int i = 0; i < i_n_se; i++){
    // Loop over coordinates.
    for (int j = 0; j < 3; j++){
      // Displace the input and output indices to point at the j'th coordinate of the first node of the i'th dislocation line segment.
      o_g_se_arr[idxf + j] = i_x3_arr[idxi + j];
      // Displace the input and output indices to point at the j'th coordinate of the second node of the i'th dislocation line segment.
      o_g_se_arr[idxf + j + 3] = i_x4_arr[idxi + j];
      o_g_se_arr[idxf + j + 6] = i_x5_arr[idxi + j];
      o_g_se_arr[idxf + j + 9] = i_x6_arr[idxi + j];
    }
    // Advance the input index to point at the first coordinate of the (i+1)'th dislocation line segment.
    idxi += 3;
    // Advance the output index to point at the location where the first coordinate of the first node of the (i+1)'th dislocation line segment should go.
    idxf += 12;
  }
  return o_g_se_arr;
}

void fx_device_host_map(double *i_g_fx_arr, double *o_fx_arr[], int i_n_se, int i_n_nodes){
  /*
    Maps the 1D nodal force array for E surface elements of N nodes each to a 2D array usable by MATLAB.
    i_g_fx_arr =
      [x_00, x_10, ..., x_(E-1)0, x_01, x_11, ..., x_(E-1)1, ..., x_0N, x_1(N-1), ..., x_(E-1)3,
       y_00, y_10, ..., y_(E-1)0, y_01, y_11, ..., y_(E-1)1, ..., y_0N, y_1(N-1), ..., y_(E-1)3,
       z_00, z_10, ..., z_(E-1)0, z_01, z_11, ..., z_(E-1)1, ..., z_0N, z_1(N-1), ..., z_(E-1)3]

       to

    o_fx_arr[n][3e + 0:2] = [x_en, y_en, z_en]
    i_n_se = E
    i_n_nodes    = N
    n_se_n_nodes = N*E
    idxi = input (initial) index
    idxf = output (final) index
  */
  int n_se_n_nodes = i_n_se*i_n_nodes;
  int idxi = 0, idxf = 0;
  // Loop over nodes in element.
  for (int i = 0; i < i_n_nodes; i++){
    // Reset the output index to point at the first element.
    idxf = 0;
    // Loop over elements.
    for (int j = 0; j < i_n_se; j++){
      // Loop over coordinates.
      for (int k = 0; k < 3; k++){
        // Displace the input index to point at the k'th coordinate of the i'th node of the j'th element.
        // Displace the output index to point at the k'th coordinate of the i'th node of the j'th element.
        o_fx_arr[i][idxf + k] = i_g_fx_arr[idxi + j + k*n_se_n_nodes];
        //if(i == 0){
        //  o_ftot_arr[idxf + k] = i_g_ftot_arr[j + k*i_n_se];
        //}
        //
      }
      // Advance the output index to point at the first coordinate of the (j+1)'th element.
      idxf += 3;
    }
    // Advance the input index to point at the first coordinate of the (i+1)'th node of the first element.
    idxi += i_n_se;
  }
}

void add_fx_device_host_map(double *i_g_fx_arr, double *o_fx_arr[], int i_n_se, int i_n_nodes){
  /*
    Maps the 1D nodal force array for E surface elements of N nodes each to a 2D array usable by MATLAB.
    i_g_fx_arr =
      [x_00, x_10, ..., x_(E-1)0, x_01, x_11, ..., x_(E-1)1, ..., x_0N, x_1(N-1), ..., x_(E-1)3,
       y_00, y_10, ..., y_(E-1)0, y_01, y_11, ..., y_(E-1)1, ..., y_0N, y_1(N-1), ..., y_(E-1)3,
       z_00, z_10, ..., z_(E-1)0, z_01, z_11, ..., z_(E-1)1, ..., z_0N, z_1(N-1), ..., z_(E-1)3]

       to

    o_fx_arr[n][3e + 0:2] = [x_en, y_en, z_en]
    i_n_se = E
    i_n_nodes    = N
    n_se_n_nodes = N*E
    idxi = input (initial) index
    idxf = output (final) index
  */
  int n_se_n_nodes = i_n_se*i_n_nodes;
  int idxi = 0, idxf = 0;
  // Loop over nodes in element.
  for (int i = 0; i < i_n_nodes; i++){
    // Reset the output index to point at the first element.
    idxf = 0;
    // Loop over elements.
    for (int j = 0; j < i_n_se; j++){
      // Loop over coordinates.
      for (int k = 0; k < 3; k++){
        // Displace the input index to point at the k'th coordinate of the i'th node of the j'th element.
        // Displace the output index to point at the k'th coordinate of the i'th node of the j'th element.
        o_fx_arr[i][idxf + k] += i_g_fx_arr[idxi + j + k*n_se_n_nodes];
        //if(i == 0){
        //  o_ftot_arr[idxf + k] += i_g_ftot_arr[j + k*i_n_se];
        //}
        //
      }
      // Advance the output index to point at the first coordinate of the (j+1)'th element.
      idxf += 3;
    }
    // Advance the input index to point at the first coordinate of the (i+1)'th node of the first element.
    idxi += i_n_se;
  }
}

void ftot_device_host_map(double *i_g_ftot_arr, double *o_ftot_arr, int i_n_se){
  /*
    Maps the 1D total force array for E surface elements to a 1D array usable by MATLAB.
    i_g_ftot_arr =
      [x_0, x_1, ..., x_(E-1)
       y_0, y_1, ..., y_(E-1)
       z_0, z_1, ..., z_(E-1)]

       to

     i_g_ftot_arr =
       [x_0, y_0, z_0, x_1, y_1, z_1, ..., x_(E-1), y_(E-1), z_(E-1)]
    i_n_se = E
    idxf = output (final) index
  */
  int idxf = 0;
  // Loop over surface elements.
  for (int i = 0; i < i_n_se; i++){
    // Loop over coordinates.
    for (int j = 0; j < 3; j++){
      // Displace the input and output indices to point at the j'th coordinate of the i'th surface element.
      o_ftot_arr[idxf + j] = i_g_ftot_arr[i + j*i_n_se];
    }
    // Advance the output index to point at the first coordinate of the (i+1)'th surface element.
    idxf += 3;
  }
}

void add_ftot_device_host_map(double *i_g_ftot_arr, double *o_ftot_arr, int i_n_se){
  /*
    Maps the 1D total force array for E surface elements to a 1D array usable by MATLAB.
    i_g_ftot_arr =
      [x_0, x_1, ..., x_(E-1)
       y_0, y_1, ..., y_(E-1)
       z_0, z_1, ..., z_(E-1)]

       to

     i_g_ftot_arr =
       [x_0, y_0, z_0, x_1, y_1, z_1, ..., x_(E-1), y_(E-1), z_(E-1)]
    i_n_se = E
    idxf = output (final) index
  */
  int idxf = 0;
  // Loop over surface elements.
  for (int i = 0; i < i_n_se; i++){
    // Loop over coordinates.
    for (int j = 0; j < 3; j++){
      // Displace the input and output indices to point at the j'th coordinate of the i'th surface element.
      o_ftot_arr[idxf + j] += i_g_ftot_arr[i + j*i_n_se];
    }
    // Advance the output index to point at the first coordinate of the (i+1)'th surface element.
    idxf += 3;
  }
}

double *dln_host_device_map(double *i_x1_arr, double *i_x2_arr, int i_n_dln){
  /*
    Maps E dislocation line segments with 2 nodes.
    This is for looping through in a parallelisation over surface elements.
    i_xn_arr[3e + 0:2] = [x_en, y_en, z_en]

      to

    o_g_dln_arr =
      [x_00, y_00, z_00, x_01, y_01, z_01, ..., x_10, y_10, z_10, x_11, y_11, z_11, ...,
       x_(E-1)0, y_(E-1)0, z_(E-1)0, x_(E-1)1, y_(E-1)1, z_(E-1)1]
    i_n_dln = E
    n_dln_nodes = 2*E
    idxi = input (initial) index
    idxf = output (final) index
  */
  int const n_dln_nodes = i_n_dln*2;
  double *o_g_dln_arr;
  int idxi = 0, idxf = 0;
  // Allocate a 1D output array of length 3*2*E.
  o_g_dln_arr = (double *) malloc(3 * n_dln_nodes * sizeof(double));
  // Loop over dislocation line segments.
  for (int i = 0; i < i_n_dln; i++){
    // Loop over coordinates.
    for (int j = 0; j < 3; j++){
      // Displace the input and output indices to point at the j'th coordinate of the first node of the i'th dislocation line segment.
      o_g_dln_arr[idxf + j] = i_x1_arr[idxi + j];
      // Displace the input and output indices to point at the j'th coordinate of the second node of the i'th dislocation line segment.
      o_g_dln_arr[idxf + j + 3] = i_x2_arr[idxi + j];
    }
    // Advance the input index to point at the first coordinate of the (i+1)'th dislocation line segment.
    idxi += 3;
    // Advance the output index to point at the location where the first coordinate of the first node of the (i+1)'th dislocation line segment should go.
    idxf += 6;
  }
  return o_g_dln_arr;
}

double *b_host_device_map(double *i_b_arr, int i_n_b){
  /*
    Maps E Burgers vectors.
    This is for looping through in a parallelisation over surface elements.
    i_b_arr[3e + 0:2] = [x_e, y_e, z_e]

      to

    o_g_dln_arr =
      [x_0, y_0, z_0, x_1, y_1, z_1, ..., x_(E-1), y_(E-1), z_(E-1)]
    i_n_b = E
    idx = index
  */
  double *o_g_b_arr;
  int idx = 0;
  // Allocate return array.
  o_g_b_arr = (double *) malloc(i_n_b * 3 * sizeof(double));
  // Loop over dislocation line segments.
  for (int i = 0; i < i_n_b; i++){
    // Loop over coordinates.
    for (int j = 0; j < 3; j++){
      o_g_b_arr[idx + j] = i_b_arr[idx + j];
    }
    // Advance index to find the value of the first coordinate at the (i+1)'th Burgers vector.
    idx += 3;
  }
  return o_g_b_arr;
}

__device__ void cuda_init_force(double nodal_force[][3], double *total_force){
  // Sets forces to zero.
  for (int i = 0; i < n_nodes; i++){
    for (int j = 0; j < 3; j++){
      nodal_force[i][j] = 0.0;
    }
  }
  for (int i = 0; i < 3; i++){
    total_force[i] = 0.0;
  }
}

__device__ void element_device_thread_map(double *i_g_elem_arr, double o_x[][3], int i_n_elem_scope, int i_n_nodes){
  /*
    Maps E elements with N nodes from the global device array to the local node arrays.
    This ensures coalesced memory accesses when parallelising over surface elements.
    o_g_elem_arr =
      [x_00, x_10, ..., x_(E-1)0, x_01, x_11, ..., x_(E-1)1, ..., x_0N, x_1(N-1), ..., x_(E-1)3,
       y_00, y_10, ..., y_(E-1)0, y_01, y_11, ..., y_(E-1)1, ..., y_0N, y_1(N-1), ..., y_(E-1)3,
       z_00, z_10, ..., z_(E-1)0, z_01, z_11, ..., z_(E-1)1, ..., z_0N, z_1(N-1), ..., z_(E-1)3]

      to

    i_x_arr[n][3e + 0:2] = [x_en, y_en, z_en]

    i_g_elem_arr    = E
    i_n_nodes       = N
    n_elem_n_nodes  = E*N
    idx  = unique index for the specific thread.
    idxi = global array index
  */
  int n_elem_n_nodes = i_n_nodes*i_n_elem_scope;
  int idx  = threadIdx.x + blockIdx.x * blockDim.x;
  int idxi = idx;
  // The lhs of this loop is not ideal due to the row-order of C but the rhs allows fully coalesced memory access of global memory.
  // Loop over coordinates.
  for (int i = 0; i < 3; i++){
    // Loop over nodes.
    for (int j = 0; j < i_n_nodes; j++){
      o_x[j][i] = i_g_elem_arr[idxi + j*i_n_elem_scope];
    }
    // Advance the global array index to point at the i'th coordinate of the first node of the idx'th element.
    idxi += n_elem_n_nodes;
  }
  /*
  // Loop over nodes.
  for (int i = 0; i < i_n_nodes; i++){
    // Loop over coordinates.
    for (int j = 0; j < 3; j++){
      // Displace the input index to point at the j'th coordinate of the i'th node of the idx'th element.
      o_x[i][j] = i_g_elem_arr[idxi + j*n_elem_n_nodes];
    }
    // Advance the input index to point at the first coordinate of the (i+1)'th node of the idx'th element.
    idxi += i_n_elem_scope;
  }
  */
}

__device__ void se_device_thread_map(double *i_g_se_arr,
                                     double *o_x3, double *o_x4, double *o_x5, double *o_x6,
                                     int i_n_se, int i_idx){
  /*
    Maps E elements with 4 nodes from the global device array to the local node arrays.
    This ensures coalesced memory accesses when parallelising over surface elements.
    o_g_se_arr =
      [x_00, x_10, ..., x_(E-1)0, x_01, x_11, ..., x_(E-1)1, ..., x_0N, x_1(N-1), ..., x_(E-1)3,
       y_00, y_10, ..., y_(E-1)0, y_01, y_11, ..., y_(E-1)1, ..., y_0N, y_1(N-1), ..., y_(E-1)3,
       z_00, z_10, ..., z_(E-1)0, z_01, z_11, ..., z_(E-1)1, ..., z_0N, z_1(N-1), ..., z_(E-1)3]

      to

    i_xn_arr[3e + 0:2] = [x_en, y_en, z_en]

    i_n_se    = E
    n_elem_n_nodes  = 4*E
    idx  = unique index for the specific thread.
    idxi = global array index
  */
  int n_elem_n_nodes = 4*i_n_se;
  int idxi = i_idx;
  // Loop over coordinates.
  for (int i = 0; i < 3; i++){
    // Displace the global array index to point at the i'th coordinate of the appropriate node of the idx'th surface element.
    o_x3[i] = i_g_se_arr[idxi           ];
    o_x4[i] = i_g_se_arr[idxi +   i_n_se];
    o_x5[i] = i_g_se_arr[idxi + 2*i_n_se];
    o_x6[i] = i_g_se_arr[idxi + 3*i_n_se];
    // Advance the global array index to point at the i'th coordinate of the first node of the idx'th element.
    idxi += n_elem_n_nodes;
  }
}

__device__ void add_force_thread_device(double i_nodal_force[][3], double *i_total_force, double *o_g_fx_arr, double *o_g_ftot_arr, int i_n_se, int i_n_nodes, int idx){
  /*
    Performs atomic addition of nodal and total forces for E surface elements with N nodes from local thread memory to global device memory.
    i_fx_arr[n][3e + 0:2] = [x_en, y_en, z_en]

      added to

    o_g_fx_arr =
      [x_00, x_10, ..., x_(E-1)0, x_01, x_11, ..., x_(E-1)1, ..., x_0N, x_1(N-1), ..., x_(E-1)3,
       y_00, y_10, ..., y_(E-1)0, y_01, y_11, ..., y_(E-1)1, ..., y_0N, y_1(N-1), ..., y_(E-1)3,
       z_00, z_10, ..., z_(E-1)0, z_01, z_11, ..., z_(E-1)1, ..., z_0N, z_1(N-1), ..., z_(E-1)3]

    and

    i_total_force[3e + 0:2] = [x_e, y_e, z_e]

      added to

    o_g_ftot_arr =
      [x_0, x_1, ..., x_(E-1)
       y_0, y_1, ..., y_(E-1)
       z_0, z_1, ..., z_(E-1)]
    i_n_se = E
    i_n_nodes    = N
    n_se_n_nodes = E*N
    idxf = output (final) index
  */
  int const n_se_n_nodes = i_n_se*i_n_nodes;
  int idxf = idx;
  // Nodal force.
  // Loop over nodes.
  for (int i = 0; i < i_n_nodes; i++){
    // Loop over coordinates.
    for (int j = 0; j < 3; j++){
      // Displace the output index to point at the j'th coordinate of the i'th node of the idx'th surface element.
      #if __CUDA_ARCH__ < 600
        atomicAdd_dbl(&o_g_fx_arr[idxf + j*n_se_n_nodes], i_nodal_force[i][j]);
      #else
        atomicAdd(&o_g_fx_arr[idxf + j*n_se_n_nodes], i_nodal_force[i][j]);
      #endif
    }
    // Advance the output index to point at the (i+1)'th node of the first coordinate of the idx'th surface element.
    idxf += i_n_se;
  }
  // Total force per surface element.
  // Loop over coordinates.
  idxf = idx;
  for (int i = 0; i < 3; i++){
    #if __CUDA_ARCH__ < 600
        atomicAdd_dbl(&o_g_ftot_arr[idxf], i_total_force[i]);
    #else
        atomicAdd(&o_g_ftot_arr[idxf], i_total_force[i]);
    #endif
    // Advance the output index to point at the (i+1)'th coordinate of the idx'th surface element.
    idxf += i_n_se;
  }
}

__device__ void dln_add_force_thread_device(double i_nodal_force[][3], double *i_total_force, double *o_g_fx_arr, double *o_g_ftot_arr, int i_n_se, int i_n_nodes, int idx){
  /*
    Performs atomic addition of nodal and total forces for E surface elements with N nodes from local thread memory to global device memory.
    i_fx_arr[n][3e + 0:2] = [x_en, y_en, z_en]

      added to

    o_g_fx_arr =
      [x_00    , y_00    , z_00    , x_01    , y_01    , z_01    , ..., x_0(N-1)    , y_0(N-1)    , z_0(N-1),
       x_10    , y_10    , z_10    , x_11    , y_11    , z_11    , ..., x_1(N-1)    , y_1(N-1)    , z_1(N-1), ...
       x_(E-1)0, y_(E-1)0, z_(E-1)0, x_(E-1)1, y_(E-1)1, z_(E-1)1, ..., x_(E-1)(N-1), y_(E-1)(N-1), z_(E-1)(N-1)]

    and

    i_total_force[3e + 0:2] = [x_e, y_e, z_e]

      added to

    o_g_ftot_arr =
      [x_0, y_0, z_0, x_1 , y_1 , z_1, ..., x_(E-1) , y_(E-1) , z_(E-1)]
    i_n_se = E
    i_n_nodes    = N
    n_se_n_nodes = E*N
    idxf = output (final) index
  */
  int const idxf = 3*i_n_nodes*idx;
  int idxa = 0;
  // Nodal force.
  // Loop over nodes.
  for (int i = 0; i < i_n_nodes; i++){
    // Loop over coordinates.
    for (int j = 0; j < 3; j++){
      // Displace the output index to point at the j'th coordinate of the i'th node of the idx'th surface element.
      #if __CUDA_ARCH__ < 600
        atomicAdd_dbl(&o_g_fx_arr[idxf + idxa + j], i_nodal_force[i][j]);
      #else
        atomicAdd(&o_g_fx_arr[idxf + idxa + j], i_nodal_force[i][j]);
      #endif
    }
    idxa += 3;
  }
  // Total force per surface element.
  // Loop over coordinates.
  idxa = 3*idx;
  for (int i= 0; i < 3; i++){
    #if __CUDA_ARCH__ < 600
        atomicAdd_dbl(&o_g_ftot_arr[idxa + i], i_total_force[i]);
    #else
        atomicAdd(&o_g_ftot_arr[idxa + i], i_total_force[i]);
    #endif
  }
}

void dln_fx_device_host_map(double *i_g_fx_arr, double *o_fx_arr[], int i_n_se, int i_n_nodes){
  /*
    Maps the 1D nodal force array for E surface elements of N nodes each to a 2D array usable by MATLAB.
    i_g_fx_arr =
      [x_00    , y_00    , z_00    , x_01    , y_01    , z_01    , ..., x_0(N-1)    , y_0(N-1)    , z_0(N-1),
       x_10    , y_10    , z_10    , x_11    , y_11    , z_11    , ..., x_1(N-1)    , y_1(N-1)    , z_1(N-1), ...
       x_(E-1)0, y_(E-1)0, z_(E-1)0, x_(E-1)1, y_(E-1)1, z_(E-1)1, ..., x_(E-1)(N-1), y_(E-1)(N-1), z_(E-1)(N-1)]

       to

    o_fx_arr[n][3e + 0:2] = [x_en, y_en, z_en]
    i_n_se = E
    i_n_nodes    = N
    n_se_n_nodes = N*E
    idxi = input (initial) index
    idxf = output (final) index
  */
  int idxi = 0, idxf = 0;
  // Loop over nodes in element.
  for (int i = 0; i < i_n_nodes; i++){
    // Loop over elements.
    for (int j = 0; j < i_n_se; j++){
      // Loop over coordinates.
      for (int k = 0; k < 3; k++){
        // Displace the input index to point at the k'th coordinate of the i'th node of the j'th element.
        // Displace the output index to point at the k'th coordinate of the i'th node of the j'th element.
        o_fx_arr[i][idxf + k] = i_g_fx_arr[idxi + idxf*i_n_nodes + k];
      }
      // Advance the output index to point at the first coordinate of the (j+1)'th element.
      idxf += 3;
    }
    // Advance the input index to point at the first coordinate of the (i+1)'th node of the first element.
    idxi += 3;
    idxf  = 0;
  }
}

void dln_add_fx_device_host_map(double *i_g_fx_arr, double *o_fx_arr[], int i_n_se, int i_n_nodes){
  /*
    Maps the 1D nodal force array for E surface elements of N nodes each to a 2D array usable by MATLAB.
    i_g_fx_arr =
      [x_00    , y_00    , z_00    , x_01    , y_01    , z_01    , ..., x_0(N-1)    , y_0(N-1)    , z_0(N-1),
       x_10    , y_10    , z_10    , x_11    , y_11    , z_11    , ..., x_1(N-1)    , y_1(N-1)    , z_1(N-1), ...
       x_(E-1)0, y_(E-1)0, z_(E-1)0, x_(E-1)1, y_(E-1)1, z_(E-1)1, ..., x_(E-1)(N-1), y_(E-1)(N-1), z_(E-1)(N-1)]

       to

    o_fx_arr[n][3e + 0:2] = [x_en, y_en, z_en]
    i_n_se = E
    i_n_nodes    = N
    n_se_n_nodes = N*E
    idxi = input (initial) index
    idxf = output (final) index
  */
  int idxi = 0, idxf = 0;
  // Loop over nodes in element.
  for (int i = 0; i < i_n_nodes; i++){
    // Loop over elements.
    for (int j = 0; j < i_n_se; j++){
      // Loop over coordinates.
      for (int k = 0; k < 3; k++){
        // Displace the input index to point at the k'th coordinate of the i'th node of the j'th element.
        // Displace the output index to point at the k'th coordinate of the i'th node of the j'th element.
        o_fx_arr[i][idxf + k] += i_g_fx_arr[idxi + idxf*i_n_nodes + k];
      }
      // Advance the output index to point at the first coordinate of the (j+1)'th element.
      idxf += 3;
    }
    // Advance the input index to point at the first coordinate of the (i+1)'th node of the first element.
    idxi += 3;
    idxf  = 0;
  }
}

__device__ void dln_device_thread_map(double *i_g_dln_arr,
                                      double *o_x1, double *o_x2,
                                      int i_dln, int idx){
  /*
    Maps E dislocation line segments with 2 nodes.
    This is for looping through in a parallelisation over surface elements.
    o_g_dln_arr =
      [x_00, y_00, z_00, x_01, y_01, z_01, ..., x_10, y_10, z_10, x_11, y_11, z_11, ...,
       x_(E-1)0, y_(E-1)0, z_(E-1)0, x_(E-1)1, y_(E-1)1, z_(E-1)1]

      to

    i_xn_arr[3e + 0:2] = [x_en, y_en, z_en]

    i_dln = E
    n_dln_nodes   = 2*E
    idxi = input (initial) index
    idxf = output (final) index
  */
  int idxi = idx;
  int n_dln_nodes = 2*i_dln;
  // Loop over coordinates.
  for (int i = 0; i < 3; i++){
    o_x1[i] = i_g_dln_arr[idxi];
    // Displace the input index to point at the i'th coordinat of the second node of the idx'th dislocaiton line segment.
    o_x2[i] = i_g_dln_arr[idxi + i_dln];
    // Advance the index to point at the (i+1)'th coordinate of the first node of the idx'th dislocation line segment.
    idxi += n_dln_nodes;
  }
}

__device__ void b_device_thread_map(double *i_g_b_arr,
                                    double *o_b,
                                    int i_n_b, int idx){
  /*
    Maps E Burgers vectors.
    This is for looping through in a parallelisation over surface elements.
    o_g_dln_arr =
      [x_0, y_0, z_0, x_1, y_1, z_1, ..., x_(E-1), y_(E-1), z_(E-1)]

      to

    i_b_arr[3e + 0:2] = [x_e, y_e, z_e]

    i_n_b = E
    idx = index
  */
  int idxi = idx;
  // Loop over coordinates.
  for (int i = 0; i < 3; i++){
    o_b[i] = i_g_b_arr[idxi];
    // Advance the input index to point at the (i+1)'th coordinate of the idx'th Burgers vector.
    idxi += i_n_b;
  }
}

__device__ double cuda_dot_product(double *i_vec1, double *i_vec2, int i_vec_size){
  // Returns the dot product of i_vec1, i_vec2.
  double result = 0.0;
  for (int i = 0; i < i_vec_size; i++){
    result += i_vec1[i]*i_vec2[i];
  }
  return result;
}

__device__ double *cuda_cross_product(double *i_vec1, double *i_vec2,
                      double *o_vec){
  // Returns the cross product of i_vec1 x i_vec2.
  o_vec[0] = i_vec1[1]*i_vec2[2] - i_vec1[2]*i_vec2[1];
  o_vec[1] = i_vec1[2]*i_vec2[0] - i_vec1[0]*i_vec2[2];
  o_vec[2] = i_vec1[0]*i_vec2[1] - i_vec1[1]*i_vec2[0];
  return o_vec;
}

__device__ double *cuda_cross_product2(double *i_vec1, double *i_vec2){
  double *o_vec;
  o_vec = (double *) malloc(3*sizeof(double));
  // Returns the cross product of i_vec1 x i_vec2.
  o_vec[0] = i_vec1[1]*i_vec2[2] - i_vec1[2]*i_vec2[1];
  o_vec[1] = i_vec1[2]*i_vec2[0] - i_vec1[0]*i_vec2[2];
  o_vec[2] = i_vec1[0]*i_vec2[1] - i_vec1[1]*i_vec2[0];
  return o_vec;
}

__device__ double *cuda_normalise_vector(double *i_vec, int i_vec_size,
                         double *o_vec){
  // Returns a normalised i_vec to o_vec.
  double mag_vec = 0.0;
  mag_vec = cuda_dot_product(i_vec, i_vec, i_vec_size);
  // Check magnitude is not zero.
  if (mag_vec == 0.0){
    //printf("ERROR: nodal_surface_force_linear_rectangle: normalise_vector: mag_vec = 0: A vector cannot have magnitude 0\n");
    //asm("trap;");
  }
  mag_vec = sqrt(mag_vec);
  for(int i = 0; i < i_vec_size; i++){
    o_vec[i] = i_vec[i]/mag_vec;
  }
  return o_vec;
}

__device__ void cuda_normalise_vector2(double *i_vec, int i_vec_size,
                       double *o_vec, double *o_mag_vec){
  // Returns a normalised i_vec to o_vec, and the magnitude of the vector in o_mag_vec.
  // Has to be a void function in order to 'return' two values, the magnitude should be passed as a reference eg:
  /*
    int    size = 3;
    double i_vec[size], o_vec[size], magnitude;
    normalise_vector2(i_vec, size, o_vec, &magnitude);
  */
  *o_mag_vec = cuda_dot_product(i_vec, i_vec, i_vec_size);
  // Check magnitude is not zero.
  if (*o_mag_vec == 0.0){
    //printf("ERROR: nodal_surface_force_linear_rectangle: normalise_vector2: o_mag_vec = 0: A vector cannot have magnitude 0\n");
    //asm("trap;");
  }
  *o_mag_vec = sqrt(*o_mag_vec);
  for(int i = 0; i < i_vec_size; i++){
    o_vec[i] = i_vec[i]/ *o_mag_vec;
  }
}

__device__ double *cuda_build_vector(double *i_x1, double *i_x2, int i_vec_size,
                     double *o_vec){
  // Returns a vector o_vec which translates the point i_x1 to i_x2.
  for (int i = 0; i < i_vec_size; i++){
    o_vec[i] = i_x2[i] - i_x1[i];
  }
  return o_vec;
}

__device__ double *cuda_init_vector(double *i_x1, double *i_x2, int i_vec_size,
                    double *io_vec){
  // Builds and returns a normalised vector io_vect.
  cuda_build_vector(i_x1, i_x2, i_vec_size, io_vec);
  cuda_normalise_vector(io_vec, i_vec_size, io_vec);
  return io_vec;
}

__device__ void cuda_init_vector2(double *i_x1  , double *i_x2, int i_vec_size,
                  double *io_vec, double *o_mag_vec){
  // Builds and returns a normalised vector io_vect, and the magnitude of the vector in o_mag_vec.
  // Has to be a void function in order to 'return' two values, the magnitude should be passed as a reference eg:
  /*
    int    size = 3;
    double i_x1[size], i_x2[size], o_vec[size], magnitude;
    init_vector2(i_x1, i_x2, size, o_vec, &magnitude);
  */
  cuda_build_vector(i_x1, i_x2, i_vec_size, io_vec);
  cuda_normalise_vector2(io_vec, i_vec_size, io_vec, o_mag_vec);
}
