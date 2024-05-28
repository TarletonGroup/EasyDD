double dot_product(double *i_vec1, double *i_vec2, int i_vec_size){
  // Returns the dot product of i_vec1, i_vec2.
  double result = 0.0;
  for (int i = 0; i < i_vec_size; i++){
    result += i_vec1[i]*i_vec2[i];
  }
  return result;
}

double *cross_product(double *i_vec1, double *i_vec2,
                      double *o_vec){
  // Returns the cross product of i_vec1 x i_vec2.
  o_vec[0] = i_vec1[1]*i_vec2[2] - i_vec1[2]*i_vec2[1];
  o_vec[1] = i_vec1[2]*i_vec2[0] - i_vec1[0]*i_vec2[2];
  o_vec[2] = i_vec1[0]*i_vec2[1] - i_vec1[1]*i_vec2[0];
  return o_vec;
}

double *cross_product2(double *i_vec1, double *i_vec2){
  double *o_vec = (double *) malloc(3*sizeof(double));
  // Returns the cross product of i_vec1 x i_vec2.
  o_vec[0] = i_vec1[1]*i_vec2[2] - i_vec1[2]*i_vec2[1];
  o_vec[1] = i_vec1[2]*i_vec2[0] - i_vec1[0]*i_vec2[2];
  o_vec[2] = i_vec1[0]*i_vec2[1] - i_vec1[1]*i_vec2[0];
  return o_vec;
}

double *normalise_vector(double *i_vec, int i_vec_size,
                         double *o_vec){
  // Returns a normalised i_vec to o_vec.
  double mag_vec = 0.0;
  mag_vec = dot_product(i_vec, i_vec, i_vec_size);
  // Check magnitude is not zero.
  if (mag_vec == 0.0){
    fprintf(stderr, "ERROR: nodal_surface_force_linear_rectangle: normalise_vector: mag_vec = 0: A vector cannot have magnitude 0\n");
    //exit(EXIT_FAILURE);
  }
  mag_vec = sqrt(mag_vec);
  for(int i = 0; i < i_vec_size; i++){
    o_vec[i] = i_vec[i]/mag_vec;
  }
  return o_vec;
}

void normalise_vector2(double *i_vec, int i_vec_size,
                       double *o_vec, double *o_mag_vec){
  // Returns a normalised i_vec to o_vec, and the magnitude of the vector in o_mag_vec.
  // Has to be a void function in order to 'return' two values, the magnitude should be passed as a reference eg:
  /*
    int    size = 3;
    double i_vec[size], o_vec[size], magnitude;
    normalise_vector2(i_vec, size, o_vec, &magnitude);
  */
  *o_mag_vec = dot_product(i_vec, i_vec, i_vec_size);
  // Check magnitude is not zero.
  if (*o_mag_vec == 0.0){
    fprintf(stderr, "ERROR: nodal_surface_force_linear_rectangle: normalise_vector2: o_mag_vec = 0: A vector cannot have magnitude 0\n");
    //exit(EXIT_FAILURE);
  }
  *o_mag_vec = sqrt(*o_mag_vec);
  for(int i = 0; i < i_vec_size; i++){
    o_vec[i] = i_vec[i]/ *o_mag_vec;
  }
}

double *arbitrary_rotation_matrix_3d(double i_theta, double *i_rot_centre, double *i_rot_axis, double *i_point, double *o_result){
  // Rotates i_point an angle of i_theta about the unit vector i_rot_axis passing through the point i_rot_centre..
  double u_sq, v_sq, w_sq, au, bv, cw, m_ux_m_vy_m_wz, costheta, one_m_costheta, sintheta;
  double mag_rot_axis;
  // Always assume the user is stupid and check whether i_rot_axis is normalised, if it's not normalise it.
  mag_rot_axis = dot_product(i_rot_axis, i_rot_axis, 3);
  if(mag_rot_axis != 1.0){
    mag_rot_axis = sqrt(mag_rot_axis);
    for (int i = 0; i < 3; i++){
      i_rot_axis[i] /= mag_rot_axis;
    }
  }
  // cos(i_theta), 1 - cos(i_theta), sin(i_theta)
  sintheta = sin(i_theta);
  costheta = cos(i_theta);
  one_m_costheta = 1. - costheta;
  // u^2, v^2, w^2
  u_sq = i_rot_axis[0]*i_rot_axis[0];
  v_sq = i_rot_axis[1]*i_rot_axis[1];
  w_sq = i_rot_axis[2]*i_rot_axis[2];
  // a*u, b*v, c*w, -u*x-v*y-w*z
  au = i_rot_centre[0]*i_rot_axis[0];
  bv = i_rot_centre[1]*i_rot_axis[1];
  cw = i_rot_centre[2]*i_rot_axis[2];
  m_ux_m_vy_m_wz = -(i_point[0]*i_rot_axis[0] + i_point[1]*i_rot_axis[1] + i_point[2]*i_rot_axis[2]);

  o_result[0] = one_m_costheta*(i_rot_centre[0]*(v_sq + w_sq) - i_rot_axis[0]*(bv + cw + m_ux_m_vy_m_wz))
            + costheta*i_point[0]
            + sintheta*(-i_rot_centre[2]*i_rot_axis[1] + i_rot_centre[1]*i_rot_axis[2] - i_rot_axis[2]*i_point[1] + i_rot_axis[1]*i_point[2]);
  o_result[1] = one_m_costheta*(i_rot_centre[1]*(u_sq + w_sq) - i_rot_axis[1]*(au + cw + m_ux_m_vy_m_wz))
            + costheta*i_point[1]
            + sintheta*( i_rot_centre[2]*i_rot_axis[0] - i_rot_centre[0]*i_rot_axis[2] + i_rot_axis[2]*i_point[0] - i_rot_axis[0]*i_point[2]);
  o_result[2] = one_m_costheta*(i_rot_centre[2]*(u_sq + v_sq) - i_rot_axis[2]*(au + bv + m_ux_m_vy_m_wz))
            + costheta*i_point[2]
            + sintheta*(-i_rot_centre[1]*i_rot_axis[0] + i_rot_centre[0]*i_rot_axis[1] - i_rot_axis[1]*i_point[0] + i_rot_axis[0]*i_point[1]);
  return o_result;
}

double *build_vector(double *i_x1, double *i_x2, int i_vec_size,
                     double *o_vec){
  // Returns a vector o_vec which translates the point i_x1 to i_x2.
  for (int i = 0; i < i_vec_size; i++){
    o_vec[i] = i_x2[i] - i_x1[i];
  }
  return o_vec;
}

double *init_vector(double *i_x1, double *i_x2, int i_vec_size,
                    double *io_vec){
  // Builds and returns a normalised vector io_vect.
  build_vector(i_x1, i_x2, i_vec_size, io_vec);
  normalise_vector(io_vec, i_vec_size, io_vec);
  return io_vec;
}

void init_vector2(double *i_x1  , double *i_x2, int i_vec_size,
                  double *io_vec, double *o_mag_vec){
  // Builds and returns a normalised vector io_vect, and the magnitude of the vector in o_mag_vec.
  // Has to be a void function in order to 'return' two values, the magnitude should be passed as a reference eg:
  /*
    int    size = 3;
    double i_x1[size], i_x2[size], o_vec[size], magnitude;
    init_vector2(i_x1, i_x2, size, o_vec, &magnitude);
  */
  build_vector(i_x1, i_x2, i_vec_size, io_vec);
  normalise_vector2(io_vec, i_vec_size, io_vec, o_mag_vec);
}
