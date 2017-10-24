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
// double precision atomicAdd
#if __CUDA_ARCH__ < 600
  __device__ double atomicAdd(double* address, double val)
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
/*
#ifdef _WIN32
  double atanh( double r ){
    return 0.5 * (log(1+r) - log(1-r));
  }
#endif
*/
const int n_nodes = 4;
const int n_limits = 8;
const float pi = 4.0 * atan(1.0);

float dot_product(float *i_vec1, float *i_vec2, int i_vec_size){
  // Returns the dot product of i_vec1, i_vec2.
  float result = 0.0;
  for (int i = 0; i < i_vec_size; i++){
    result += i_vec1[i]*i_vec2[i];
  }
  return result;
}

float *cross_product(float *i_vec1, float *i_vec2,
                      float *o_vec){
  // Returns the cross product of i_vec1 x i_vec2.
  o_vec[0] = i_vec1[1]*i_vec2[2] - i_vec1[2]*i_vec2[1];
  o_vec[1] = i_vec1[2]*i_vec2[0] - i_vec1[0]*i_vec2[2];
  o_vec[2] = i_vec1[0]*i_vec2[1] - i_vec1[1]*i_vec2[0];
  return o_vec;
}

float *cross_product2(float *i_vec1, float *i_vec2){
  float *o_vec = (float *) malloc(3*sizeof(float));
  // Returns the cross product of i_vec1 x i_vec2.
  o_vec[0] = i_vec1[1]*i_vec2[2] - i_vec1[2]*i_vec2[1];
  o_vec[1] = i_vec1[2]*i_vec2[0] - i_vec1[0]*i_vec2[2];
  o_vec[2] = i_vec1[0]*i_vec2[1] - i_vec1[1]*i_vec2[0];
  return o_vec;
}

float *normalise_vector(float *i_vec, int i_vec_size,
                         float *o_vec){
  // Returns a normalised i_vec to o_vec.
  float mag_vec = 0.0;
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

void normalise_vector2(float *i_vec, int i_vec_size,
                       float *o_vec, float *o_mag_vec){
  // Returns a normalised i_vec to o_vec, and the magnitude of the vector in o_mag_vec.
  // Has to be a void function in order to 'return' two values, the magnitude should be passed as a reference eg:
  /*
    int    size = 3;
    float i_vec[size], o_vec[size], magnitude;
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

float *arbitrary_rotation_matrix_3d(float i_theta, float *i_rot_centre, float *i_rot_axis, float *i_point, float *o_result){
  // Rotates i_point an angle of i_theta about the unit vector i_rot_axis passing through the point i_rot_centre..
  float u_sq, v_sq, w_sq, au, bv, cw, m_ux_m_vy_m_wz, costheta, one_m_costheta, sintheta;
  float mag_rot_axis;
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

float *build_vector(float *i_x1, float *i_x2, int i_vec_size,
                     float *o_vec){
  // Returns a vector o_vec which translates the point i_x1 to i_x2.
  for (int i = 0; i < i_vec_size; i++){
    o_vec[i] = i_x2[i] - i_x1[i];
  }
  return o_vec;
}

float *init_vector(float *i_x1, float *i_x2, int i_vec_size,
                    float *io_vec){
  // Builds and returns a normalised vector io_vect.
  build_vector(i_x1, i_x2, i_vec_size, io_vec);
  normalise_vector(io_vec, i_vec_size, io_vec);
  return io_vec;
}

void init_vector2(float *i_x1  , float *i_x2, int i_vec_size,
                  float *io_vec, float *o_mag_vec){
  // Builds and returns a normalised vector io_vect, and the magnitude of the vector in o_mag_vec.
  // Has to be a void function in order to 'return' two values, the magnitude should be passed as a reference eg:
  /*
    int    size = 3;
    float i_x1[size], i_x2[size], o_vec[size], magnitude;
    init_vector2(i_x1, i_x2, size, o_vec, &magnitude);
  */
  build_vector(i_x1, i_x2, i_vec_size, io_vec);
  normalise_vector2(io_vec, i_vec_size, io_vec, o_mag_vec);
}

float init_point(float *i_vec1, float *i_vec2,
                  float *i_vec3, float *i_vec4,
                  int     i_vec_size){
  // Initialises points on the surface element given four vectors.
  float result = 0.0;
  float denom = 0.0;
  denom = dot_product(i_vec3, i_vec4, i_vec_size);
  if (denom == 0.0){
    fprintf(stderr, "nodal_surface_force_linear_rectangle: init_point: division by zero, dot_product(i_vec3, i_vec4) = %2.14f\n", denom);
    //exit(EXIT_FAILURE);
  }
  result = dot_product(i_vec1, i_vec2, i_vec_size)/denom;
  return result;
}

float seed_single_integral(float a, float b){
  // Single integral seed.
  float result = 0.0;
  float arg;
  arg = a+b;
  if (arg <= 0.0){
    fprintf(stderr, "nodal_surface_force_linear_rectangle: seed_single_integral: log(%2.14f) = %2.14f\n", arg, log(arg));
    //exit(EXIT_FAILURE);
  }
  result = log(arg);
  return result;
}

float integral_type_1(float a, float b,
                       float c, float d,
                       float e){
  float result = 0.0;
  result = 0.5*(a*b + (c-d)*e);
  return result;
}

float numer_seed_double_integral(float a,
                                  float b, float c, float d,
                                  float e){
  // Returns the argument of the numerator of the seed integral for float integrals.
  float result = 0.0;
  result = a*(b - c + d) + e;
  return result;
}

float denom_seed_double_integral(float a, float b,
                                  float c, float d){
  // Returns the argument of the denominator of the seed integral for float integrals.
  float result = 0.0;
  result = a*b + c*d;
  if(result == 0.0){
    fprintf(stderr, "nodal_surface_force_linear_rectangle: denom_seed_integral_00m3: division by 0, result = %2.14f\n", result);
    //exit(EXIT_FAILURE);
  }
  return result;
}

float seed_double_integral(float a, float b, float c, float d, float e,
                            float f, float g, float h, float i){
  // Returns the seed integral for float integrals.
  float result      = 0.0;
  float numerator   = 0.0;
  float denominator = 0.0;
  numerator   = numer_seed_double_integral(a, b, c, d, e);
  denominator = denom_seed_double_integral(f, g, h, i);
  if(denominator > 0.){
    denominator = sqrt(denominator);
    result = 2.0/denominator * atan(numerator/denominator);
  }
  else{
    denominator = sqrt(-denominator);
    result = -2.0/denominator * atanh(numerator/denominator);
  }
  return result;
}

float integral_type_2(float a,
                       float b, float c){
  float result = 0.0;
  result = a + b*c;
  return result;
}

float integral_type_3(float a,
                       float b, float c,
                       float d, float e){
  float result = 0.0;
  result = a + b*c + d*e;
  return result;
}

float integral_type_4(float a,
                       float b, float c,
                       float d, float e,
                       float f, float g){
  float result = 0.0;
  result = a + b*c + d*e + f*g;
  return result;
}

float integral_type_5(float a,
                       float b, float c,
                       float d, float e,
                       float f, float g,
                       float h, float i){
  float result = 0.0;
  result = a + b*c + d*e + f*g + h*i;
  return result;
}

float integral_type_6(float a, float b,
                       float c, float d,
                       float e, float f, float g,
                       float h, float i,
                       float j, float k){
  float result = 0.0;
  result = a*b + c*d + (e+f)*g + h*i + j*k;
  return result;
}

float integral_type_7(float a, float b,
                       float c, float d,
                       float e, float f, float g,
                       float h, float i){
  float result = 0.0;
  result = a*b + c*d + (e+f)*g + h*i;
  return result;
}

float integral_type_8(float a, float b,
                       float c, float d,
                       float e, float f,
                       float g, float h){
  float result = 0.0;
  result = a*b + c*d + e*f + g*h;
  return result;
}

float integral_type_9(float a,
                       float b, float c,
                       float d, float e,
                       float f, float g,
                       float h, float i,
                       float j, float k){
  float result = 0.0;
  result = a + b*c + d*e + f*g + h*i + j*k;
  return result;
}

float integral_type_10(float a,
                        float b, float c,
                        float d, float e,
                        float f, float g,
                        float h, float i,
                        float j, float k){
  float result = 0.0;
  result = a + b*c + d*e + f*g + h*i + j*k;
  return result;
}

float integral_type_11(float a,
                        float b, float c,
                        float d, float e,
                        float f, float g,
                        float h, float i){
  float result = 0.0;
  result = a + b*c + d*e + f*g + h*i;
  return result;
}

float integral_type_12(float a, float b,
                        float c, float d,
                        float e, float f,
                        float g, float h,
                        float i, float j,
                        float k, float l, float m){
  float result = 0.0;
  result = a*b + c*d + e*f + g*h + i*j + k*l*m;
  return result;
}

void integral_vector(float *i_p, float *i_q, float *i_b, float *i_t, float *i_n,
                     float i_one_m_nu, float i_a_sq,
                     float o_vec_int[][3]){
  // Calculates the vectors that are multiplied by the integrals.
  // Has to be void in order to return a 2D array, otherwise we'd have to use dynamic memory allocation and therefore pointers. We don't need such flexibility as this is quite specific.
  /*
    (t x b), (p x b)
    (q x b), (b x t)
    t*((p x b) dot n), t*((q x b) dot n)
    t*((t x b) dot n), t*((b x t) dot n)
    (t dot n) * (p x b), (t dot n) * (q x b)
    (t dot n) * (t x b)
    n * ((p x b) dot t), n * ((q x b) dot t)
  */
  float t_x_b[3], p_x_b[3],
         q_x_b[3], b_x_t[3],
         t_p_x_b_dot_n[3], t_q_x_b_dot_n[3],
         t_t_x_b_dot_n[3], t_b_x_t_dot_n[3],
         t_dot_n_p_x_b[3], t_dot_n_q_x_b[3],
         t_dot_n_t_x_b[3],
         n_p_x_b_dot_t[3], n_q_x_b_dot_t[3];
  /*
    t dot n
    ((p x b) dot n), ((q x b) dot n)
    ((t x b) dot n), ((b x t) dot n)
    ((p x b) dot t), ((q x b) dot t)
    (t dot n) * ((p x b) dot t), (t dot n) * ((q x b) dot t)
    3*(t dot n) * ((p x b) dot t), 3*(t dot n) * ((q x b) dot t)
    1.5 * (1-nu) * a^2, a^3
  */
  float t_dot_n,
         p_x_b_dot_n, q_x_b_dot_n,
         t_x_b_dot_n, b_x_t_dot_n,
         p_x_b_dot_t, q_x_b_dot_t,
         t_dot_n_p_x_b_dot_t, t_dot_n_q_x_b_dot_t,
         t_dot_n_p_x_b_dot_t_3, t_dot_n_q_x_b_dot_t_3,
         one_m_nu_1p5_a_sq, a_sq_3;

  // Calculate auxiliary local variables.
  one_m_nu_1p5_a_sq = 1.5 * i_one_m_nu * i_a_sq;
  a_sq_3 = 3.0 * i_a_sq;
  cross_product(i_p, i_b, p_x_b);
  cross_product(i_q, i_b, q_x_b);
  cross_product(i_t, i_b, t_x_b);
  cross_product(i_b, i_t, b_x_t);
  t_dot_n     = dot_product(i_t  , i_n, 3);
  p_x_b_dot_n = dot_product(p_x_b, i_n, 3);
  q_x_b_dot_n = dot_product(q_x_b, i_n, 3);
  t_x_b_dot_n = dot_product(t_x_b, i_n, 3);
  b_x_t_dot_n = dot_product(b_x_t, i_n, 3);
  p_x_b_dot_t = dot_product(p_x_b, i_t, 3);
  q_x_b_dot_t = dot_product(q_x_b, i_t, 3);
  t_dot_n_p_x_b_dot_t   = t_dot_n * p_x_b_dot_t;
  t_dot_n_q_x_b_dot_t   = t_dot_n * q_x_b_dot_t;
  t_dot_n_p_x_b_dot_t_3 = 3.0 * t_dot_n_p_x_b_dot_t;
  t_dot_n_q_x_b_dot_t_3 = 3.0 * t_dot_n_q_x_b_dot_t;

  for (int i=0; i<3; i++)
  {
      // Preparing common array factors.
      t_t_x_b_dot_n[i] = i_t  [i]*t_x_b_dot_n;
      t_b_x_t_dot_n[i] = i_t  [i]*b_x_t_dot_n;
      t_p_x_b_dot_n[i] = i_t  [i]*p_x_b_dot_n;
      t_q_x_b_dot_n[i] = i_t  [i]*q_x_b_dot_n;
      t_dot_n_t_x_b[i] = t_x_b[i]*t_dot_n;
      t_dot_n_q_x_b[i] = q_x_b[i]*t_dot_n;
      t_dot_n_p_x_b[i] = p_x_b[i]*t_dot_n;
      n_q_x_b_dot_t[i] = q_x_b_dot_t*i_n[i];
      n_p_x_b_dot_t[i] = p_x_b_dot_t*i_n[i];
      // Calculating vectors.
      o_vec_int [0][i] = (t_dot_n_t_x_b[i] + t_t_x_b_dot_n[i])*i_one_m_nu        + b_x_t[i]*t_dot_n + t_b_x_t_dot_n[i]; // checked
      o_vec_int [1][i] = (t_dot_n_q_x_b[i] + t_q_x_b_dot_n[i])*i_one_m_nu        - n_q_x_b_dot_t[i] + i_q[i]*b_x_t_dot_n; // checked
      o_vec_int [2][i] = (t_dot_n_p_x_b[i] + t_p_x_b_dot_n[i])*i_one_m_nu        - n_p_x_b_dot_t[i] + i_p[i]*b_x_t_dot_n; // checked
      o_vec_int [3][i] = (t_dot_n_t_x_b[i] + t_t_x_b_dot_n[i])*one_m_nu_1p5_a_sq; // checked
      o_vec_int [4][i] = (t_dot_n_q_x_b[i] + t_q_x_b_dot_n[i])*one_m_nu_1p5_a_sq - n_q_x_b_dot_t[i]*a_sq_3; // checked
      o_vec_int [5][i] = (t_dot_n_p_x_b[i] + t_p_x_b_dot_n[i])*one_m_nu_1p5_a_sq - n_p_x_b_dot_t[i]*a_sq_3; // checked
      o_vec_int [6][i] = - t_dot_n_q_x_b_dot_t_3*i_q[i]; // checked
      o_vec_int [7][i] = - t_dot_n_p_x_b_dot_t_3*i_p[i]; // checked
      o_vec_int [8][i] = - t_dot_n_q_x_b_dot_t_3*i_t[i]; // checked
      o_vec_int [9][i] = - t_dot_n_p_x_b_dot_t_3*i_t[i]; // checked
      o_vec_int[10][i] = - t_dot_n_p_x_b_dot_t_3*i_q[i] - t_dot_n_q_x_b_dot_t_3*i_p[i]; // checked
  }
}

float *vertex_force_linear_rectangle(float *i_sch, float i_vec_int[][3],
                                      float i_r, float i_s, float i_rs,
                                      float i_factor,
                                      float *o_force){
  // Calculates the force on a vertex.
  float f[11];
  f [0] = i_sch [3] - i_s*i_sch [6] - i_r*i_sch [8] + i_rs*i_sch [0];
  f [1] = i_sch [5] - i_s*i_sch [7] - i_r*i_sch[10] + i_rs*i_sch [2];
  f [2] = i_sch [4] - i_s*i_sch [9] - i_r*i_sch [7] + i_rs*i_sch [1];
  f [3] = i_sch[19] - i_s*i_sch[14] - i_r*i_sch[15] + i_rs*i_sch[11];
  f [4] = i_sch[21] - i_s*i_sch[16] - i_r*i_sch[18] + i_rs*i_sch[13];
  f [5] = i_sch[20] - i_s*i_sch[17] - i_r*i_sch[16] + i_rs*i_sch[12];
  f [6] = i_sch[35] - i_s*i_sch[28] - i_r*i_sch[32] + i_rs*i_sch[22];
  f [7] = i_sch[36] - i_s*i_sch[31] - i_r*i_sch[27] + i_rs*i_sch[23];
  f [8] = i_sch[33] - i_s*i_sch[26] - i_r*i_sch[30] + i_rs*i_sch[25];
  f [9] = i_sch[37] - i_s*i_sch[29] - i_r*i_sch[26] + i_rs*i_sch[24];
  f[10] = i_sch[34] - i_s*i_sch[27] - i_r*i_sch[28] + i_rs*i_sch[19];
  for (int i = 0; i < 11; i++){
    for (int j = 0; j < 3; j++){
      o_force[j] += i_vec_int[i][j] * f[i];
    }
  }
  o_force[0] = o_force[0]*i_factor;
  o_force[1] = o_force[1]*i_factor;
  o_force[2] = o_force[2]*i_factor;
  return o_force;
}

float *integrals_linear_rectangle(float *i_r, float *i_s, float *i_y,
                                   float *i_p, float *i_q, float *i_t,
                                   float i_p_dot_t, float i_q_dot_t, float i_a_sq,
                                   float *o_sch, int i_num_integrals){
  // Sign vector for quick evaluation of integrals via the dot product.
  static float signv[8] = {1.0,-1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0};
  static const float third = 1.0/3.0;
  static const float two_third = 2.0/3.0;
  // Integrals.
  float // Single integrals.
         a0m1[n_limits], b0m1[n_limits], c0m1[n_limits], a1m1[n_limits], b1m1[n_limits], c1m1[n_limits],
         a01 [n_limits], b01 [n_limits], c01 [n_limits], a11 [n_limits], b11 [n_limits],
         // float integrals.
         d00m3[n_limits], e00m3[n_limits], f00m3[n_limits], d01m3[n_limits], d10m3[n_limits], d11m3[n_limits], e01m3[n_limits], e10m3[n_limits],
         e11m3[n_limits], f01m3[n_limits], f10m3[n_limits], f11m3[n_limits], d00m1[n_limits], e00m1[n_limits], f00m1[n_limits], d01m1[n_limits],
         e01m1[n_limits], f01m1[n_limits], d10m1[n_limits], e10m1[n_limits], f10m1[n_limits], d11m1[n_limits], e11m1[n_limits], f11m1[n_limits],
         d02m3[n_limits], d20m3[n_limits], e20m3[n_limits], f20m3[n_limits], d001 [n_limits], e001 [n_limits], f001 [n_limits], d02m1[n_limits],
         d20m1[n_limits], e20m1[n_limits], f20m1[n_limits], d12m3[n_limits], d21m3[n_limits], e21m3[n_limits], f21m3[n_limits], d22m3[n_limits],
         // Triple integrals
         h [38][n_limits], h001m1[n_limits], h010m1[n_limits], h100m1[n_limits], h021m3[n_limits], h201m3[n_limits], h000m1[n_limits];
 /*
   y^2, r^2, s^2 coordinates
   y*(p dot t), y*(q dot t), r*(p dot t), r*(q dot t)
    r dot p   ,  r dot q   ,  r dot t
   (r dot p)^2, (r dot q)^2, (r dot t)^2
 */
 float y_sq[n_limits], r_sq[n_limits], s_sq[n_limits],
        y_p_dot_t[n_limits], y_q_dot_t[n_limits], r_p_dot_t[n_limits], s_q_dot_t[n_limits],
        r_dot_p[n_limits], r_dot_q[n_limits], r_dot_t[n_limits],
        r_dot_p_sq[n_limits], r_dot_q_sq[n_limits], r_dot_t_sq[n_limits];
  /*
   ra = sqrt((r_vec dot r_vec) + a**2) from non-singular dislocation theory, there are 8 r_vec, thus 8 ra's.
   ra_sq = ra^2 (element-wise squaring)
   ra_c_o_3 = 1/3 * ra^3 (element-wise cubing and division)
  */
  float ra[n_limits], ra_sq[n_limits], ra_c_o_3[n_limits];
  // Vector R from x1 to x3, x4, x5, x6 and x2 to x3, x4, x5, x6. 8 vectors with 3 components each.
  float r_vec[8][3];
  /*
   Auxiliary constants to reduce computational cost.
   2*(p dot t)
   2*(q dot t)
   1-(p dot t)
   1-(q dot t)
   1-(p dot t)^2
   1-(q dot t)^2
   1/(1-(p dot t)^2)
   1/(1-(q dot t)^2)
   1-(p dot t)^2-(q dot t)^2
   1/(1-(p dot t)^2-(q dot t)^2)
   1/3 * 1/(1-(p dot t)^2-(q dot t)^2)
   (p dot t)*(q dot t)
  */
  float                         two_p_dot_t,
                                 two_q_dot_t,
                               one_m_p_dot_t,
                               one_m_q_dot_t,
                               one_m_p_dot_t_sq,
                               one_m_q_dot_t_sq,
                         one_o_one_m_p_dot_t_sq,
                         one_o_one_m_q_dot_t_sq,
                  one_m_p_dot_t_sq_m_q_dot_t_sq,
            one_o_one_m_p_dot_t_sq_m_q_dot_t_sq,
        trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq,
        p_dot_t_q_dot_t;
  // Auxiliary scalars.
  two_p_dot_t = 2.0*i_p_dot_t;
  two_q_dot_t = 2.0*i_q_dot_t;
  one_m_p_dot_t = 1.0 - i_p_dot_t;
  one_m_q_dot_t = 1.0 - i_q_dot_t;
  one_m_p_dot_t_sq = 1.0 - i_p_dot_t * i_p_dot_t;
  one_m_q_dot_t_sq = 1.0 - i_q_dot_t * i_q_dot_t;
  p_dot_t_q_dot_t  =       i_q_dot_t * i_p_dot_t;
  one_o_one_m_p_dot_t_sq = 1.0 / one_m_p_dot_t_sq;
  one_o_one_m_q_dot_t_sq = 1.0 / one_m_q_dot_t_sq;
  one_m_p_dot_t_sq_m_q_dot_t_sq = one_m_p_dot_t_sq + one_m_q_dot_t_sq - 1.0;
  one_o_one_m_p_dot_t_sq_m_q_dot_t_sq = 1.0/one_m_p_dot_t_sq_m_q_dot_t_sq;
  trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq = third * one_o_one_m_p_dot_t_sq_m_q_dot_t_sq;
  // Auxiliary arrays.
  // Calculate r for all 8 combinations (from x1 to x3, x4, x5, x6 and x2 to x3, x4, x5, x6).
  for (int i = 0; i < n_limits; i++){
    for (int j = 0; j < 3; j++){
      r_vec[i][j] = i_r[i] * i_p[j] + i_s[i] * i_q[j] + i_y[i] * i_t[j];
    }
  }
  for (int i = 0; i < n_limits; i++){
    r_sq      [i] = i_r[i] * i_r[i];
    s_sq      [i] = i_s[i] * i_s[i];
    y_sq      [i] = i_y[i] * i_y[i];
    ra        [i] = sqrt(r_sq[i]+s_sq[i]+y_sq[i] + two_p_dot_t*i_r[i]*i_y[i] + two_q_dot_t*i_s[i]*i_y[i] + i_a_sq);
    ra_sq     [i] = ra   [i] * ra[i];
    ra_c_o_3  [i] = ra_sq[i] * ra[i]/3.0;
    r_dot_p   [i] = dot_product(r_vec[i], i_p, 3);
    r_dot_q   [i] = dot_product(r_vec[i], i_q, 3);
    r_dot_t   [i] = dot_product(r_vec[i], i_t, 3);
    r_dot_p_sq[i] = r_dot_p[i] * r_dot_p[i];
    r_dot_q_sq[i] = r_dot_q[i] * r_dot_q[i];
    r_dot_t_sq[i] = r_dot_t[i] * r_dot_t[i];
    y_p_dot_t [i] = i_y[i] * i_p_dot_t;
    y_q_dot_t [i] = i_y[i] * i_q_dot_t;
    r_p_dot_t [i] = i_r[i] * i_p_dot_t;
    s_q_dot_t [i] = i_s[i] * i_q_dot_t;
  }
  // Calculate integrals.
  for (int i = 0; i < n_limits; i++){
    // Linear seed integrals.
    a0m1[i] = seed_single_integral(ra[i], r_dot_p[i]); // checked
    b0m1[i] = seed_single_integral(ra[i], r_dot_q[i]); // checked
    c0m1[i] = seed_single_integral(ra[i], r_dot_t[i]); // checked
    // Linear integrals.
    // result = 0.5*(a*b + (c-d)*e)
    a01 [i] = integral_type_1(ra[i], r_dot_p[i], ra_sq[i], r_dot_p_sq[i], a0m1[i]); // checked
    b01 [i] = integral_type_1(ra[i], r_dot_q[i], ra_sq[i], r_dot_q_sq[i], b0m1[i]); // checked
    c01 [i] = integral_type_1(ra[i], r_dot_t[i], ra_sq[i], r_dot_t_sq[i], c0m1[i]); // checked
    // type_2 = a + b*c
    a1m1[i] = integral_type_2(ra[i]      , y_p_dot_t[i]             , -a0m1[i]); // checked
    b1m1[i] = integral_type_2(ra[i]      , y_q_dot_t[i]             , -b0m1[i]); // checked
    c1m1[i] = integral_type_2(ra[i]      , r_p_dot_t[i]+s_q_dot_t[i], -c0m1[i]); // checked
    a11 [i] = integral_type_2(ra_c_o_3[i], y_p_dot_t[i]             , -a01 [i]); // checked
    b11 [i] = integral_type_2(ra_c_o_3[i], y_q_dot_t[i]             , -b01 [i]); // checked
    // float seed integrals.
    d00m3[i] = seed_double_integral(1.0, ra[i], r_dot_p[i], r_dot_q[i], 0.0,
                                    1.0             , i_a_sq, one_m_p_dot_t_sq_m_q_dot_t_sq, y_sq[i]); // checked
    e00m3[i] = seed_double_integral(one_m_p_dot_t, ra[i], i_r[i], i_y[i], s_q_dot_t[i],
                                    one_m_p_dot_t_sq, i_a_sq, one_m_p_dot_t_sq_m_q_dot_t_sq, s_sq[i]); // checked
    f00m3[i] = seed_double_integral(one_m_q_dot_t, ra[i], i_s[i], i_y[i], r_p_dot_t[i],
                                    one_m_q_dot_t_sq, i_a_sq, one_m_p_dot_t_sq_m_q_dot_t_sq, r_sq[i]); // checked
    // float integrals.
    // type_2 = a + b*c
    d01m3[i] = integral_type_2(-a0m1[i], -y_q_dot_t[i], d00m3[i]); // checked
    d10m3[i] = integral_type_2(-b0m1[i], -y_p_dot_t[i], d00m3[i]); // checked
    // type_3 = a + b*c + d*e
    e01m3[i] = -one_o_one_m_p_dot_t_sq * integral_type_3(a0m1[i], s_q_dot_t[i], e00m3[i], -i_p_dot_t, c0m1[i]); // checked
    f01m3[i] = -one_o_one_m_q_dot_t_sq * integral_type_3(b0m1[i], r_p_dot_t[i], f00m3[i], -i_q_dot_t, c0m1[i]); // checked
    // type_2 = a + b*c
    e10m3[i] = integral_type_2(-c0m1[i], -i_p_dot_t, e01m3[i]); // checked
    f10m3[i] = integral_type_2(-c0m1[i], -i_q_dot_t, f01m3[i]); // checked
    // type_6 = a*b + c*d + (e+f)*g + h*i + j*k
    d00m1[i] = integral_type_6(i_r[i], b0m1[i], i_s[i], a0m1[i], i_a_sq, y_sq[i], -d00m3[i], -y_p_dot_t[i], d10m3[i], -y_q_dot_t[i], d01m3[i]); // checked
    // type_7 = a*b + c*d + (e+f)*g + h*i
    e00m1[i] = integral_type_7(i_r[i], c0m1[i], i_y[i], a0m1[i], i_a_sq, s_sq[i], -e00m3[i], -s_q_dot_t[i], e01m3[i]); // checked
    f00m1[i] = integral_type_7(i_s[i], c0m1[i], i_y[i], b0m1[i], i_a_sq, r_sq[i], -f00m3[i], -r_p_dot_t[i], f01m3[i]); // checked
    // type_2 = a + b*c
    d11m3[i] = integral_type_2(-a1m1[i], -y_q_dot_t[i], d10m3[i]); // checked
    // type_4 = a + b*c + d*e + f*g
    e11m3[i] = one_o_one_m_p_dot_t_sq * integral_type_4(-a1m1[i], r_p_dot_t[i], c0m1[i], -i_p_dot_t, e00m1[i], -s_q_dot_t[i], e10m3[i]); // checked
    f11m3[i] = one_o_one_m_q_dot_t_sq * integral_type_4(-b1m1[i], s_q_dot_t[i], c0m1[i], -i_q_dot_t, f00m1[i], -r_p_dot_t[i], f10m3[i]); // checked
    // type_3 = a + b*c + d*e
    d20m3[i] = integral_type_3(d00m1[i], -y_p_dot_t[i], d10m3[i], -i_r[i], b0m1[i]); // checked
    d02m3[i] = integral_type_3(d00m1[i], -y_q_dot_t[i], d01m3[i], -i_s[i], a0m1[i]); // checked
    e20m3[i] = integral_type_3(e00m1[i], -i_p_dot_t     , e11m3[i], -i_r[i], c0m1[i]); // checked
    f20m3[i] = integral_type_3(f00m1[i], -i_q_dot_t     , f11m3[i], -i_s[i], c0m1[i]); // checked
    // type_2 = a + b*c
    d01m1[i] = integral_type_2(a01[i], -y_q_dot_t[i], d00m1[i]); // checked
    d10m1[i] = integral_type_2(b01[i], -y_p_dot_t[i], d00m1[i]); // checked
    e01m1[i] = one_o_one_m_p_dot_t_sq * integral_type_3(a01[i], -s_q_dot_t[i], e00m1[i], -i_p_dot_t, c01[i]); // checked
    f01m1[i] = one_o_one_m_q_dot_t_sq * integral_type_3(b01[i], -r_p_dot_t[i], f00m1[i], -i_q_dot_t, c01[i]); // checked
    // type_2 = a + b*c
    e10m1[i] = integral_type_2(c01[i], -i_p_dot_t, e01m1[i]); // checked
    f10m1[i] = integral_type_2(c01[i], -i_q_dot_t, f01m1[i]); // checked
    // type_6 = a*b + c*d + (e+f)*g + h*i + j*k
    d001[i]  = third * integral_type_6(i_r[i], b01[i], i_s[i], a01[i], y_sq[i], i_a_sq, d00m1[i], y_p_dot_t[i], d10m1[i], y_q_dot_t[i], d01m1[i]); // checked
    // type_7 = a*b + c*d + (e+f)*g + h*i
    e001[i]  = third * integral_type_7(i_r[i], c01[i], i_y[i], a01[i], s_sq[i], i_a_sq, e00m1[i], s_q_dot_t[i], e01m1[i]); // checked
    f001[i]  = third * integral_type_7(i_s[i], c01[i], i_y[i], b01[i], r_sq[i], i_a_sq, f00m1[i], r_p_dot_t[i], f01m1[i]); // checked
    // type_2 = a + b*c
    d11m1[i] = integral_type_2(a11[i], -y_q_dot_t[i], d10m1[i]); // checked
    // type_4 = a + b*c + d*e + f*g
    e11m1[i] = one_o_one_m_p_dot_t_sq * integral_type_4(a11[i], -r_p_dot_t[i], c01[i], i_p_dot_t , e001[i], -s_q_dot_t[i], e10m1[i]); // checked
    f11m1[i] = one_o_one_m_q_dot_t_sq * integral_type_4(b11[i], -s_q_dot_t[i], c01[i], i_q_dot_t , f001[i], -r_p_dot_t[i], f10m1[i]); // checked
    // type_3 = a + b*c + d*e
    d02m1[i] = integral_type_3(-d001[i], -y_q_dot_t[i], d01m1[i],  i_s[i], a01 [i]); // checked
    d20m1[i] = integral_type_3(-d001[i], -y_p_dot_t[i], d10m1[i],  i_r[i], b01 [i]); // checked
    e20m1[i] = integral_type_3(-e001[i], -i_p_dot_t     , e11m1[i],  i_r[i], c01 [i]); // checked
    f20m1[i] = integral_type_3(-f001[i], -i_q_dot_t     , f11m1[i],  i_s[i], c01 [i]); // checked
    d12m3[i] = integral_type_3(d10m1[i], -y_q_dot_t[i], d11m3[i], -i_s[i], a1m1[i]); // checked
    d21m3[i] = integral_type_3(d01m1[i], -y_p_dot_t[i], d11m3[i], -i_r[i], b1m1[i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    e21m3[i] = one_o_one_m_p_dot_t_sq * integral_type_5(e01m1[i], y_p_dot_t[i], a1m1[i], -i_p_dot_t, e10m1[i], p_dot_t_q_dot_t*i_s[i], e11m3[i], -i_r[i], c1m1[i]); // checked
    f21m3[i] = one_o_one_m_q_dot_t_sq * integral_type_5(f01m1[i], y_q_dot_t[i], b1m1[i], -i_q_dot_t, f10m1[i], p_dot_t_q_dot_t*i_r[i], f11m3[i], -i_s[i], c1m1[i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    // type_8 = a*b + c*d + e*f + g*h
    d22m3[i] = integral_type_5(-d001[i], p_dot_t_q_dot_t*y_sq[i], d11m3[i], i_r[i], b01[i], -i_r[i]*i_s[i], ra[i], i_s[i], a01[i])
             + i_y[i] * integral_type_8(-i_p_dot_t, d10m1[i], -i_q_dot_t, d01m1[i], i_p_dot_t*i_s[i], a1m1[i], i_q_dot_t*i_r[i], b1m1[i]); // checked
    // Triple integrals
    // type_3 = a + b*c + d*e
    h[0][i] = -one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_3(d00m1[i], -i_q_dot_t, e00m1[i], -i_p_dot_t, f00m1[i]); // checked
    // type_2 = a + b*c
    h[1][i] = integral_type_2(-f00m1[i], -i_p_dot_t, h[0][i]); // checked
    h[2][i] = integral_type_2(-e00m1[i], -i_q_dot_t, h[0][i]); // checked
    // type_3 = a + b*c + d*e
    h001m1[i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_3(d001[i], -i_q_dot_t, e001[i], -i_p_dot_t, f001[i]); // checked
    // type_2 = a + b*c
    h100m1[i] = integral_type_2(f001[i], -i_p_dot_t, h001m1[i]); // checked
    h010m1[i] = integral_type_2(e001[i], -i_q_dot_t, h001m1[i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    h[3][i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_5(-d11m1[i], r_p_dot_t[i], f10m1[i], -i_p_dot_t, h010m1[i], s_q_dot_t[i], e10m1[i], -i_q_dot_t, h100m1[i]); // checked
    // type_3 = a + b*c + d*e
    h[4][i] = integral_type_3(h010m1[i], -i_p_dot_t, h[3][i], -i_r[i], f10m1[i]); // checked
    h[5][i] = integral_type_3(h100m1[i], -i_q_dot_t, h[3][i], -i_s[i], e10m1[i]); // checked
    // type_4 = a + b*c + d*e + f*g
    h021m3[i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d02m1[i], -two_q_dot_t, h010m1[i], i_p_dot_t, f20m1[i], i_q_dot_t*s_sq[i], e00m1[i]); // checked
    h201m3[i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d20m1[i], -two_p_dot_t, h100m1[i], i_q_dot_t, e20m1[i], i_p_dot_t*r_sq[i], f00m1[i]);
    // type_8 = a*b + c*d + e*f + g*h
    h000m1[i] = 0.5 * integral_type_8(-i_a_sq, 0.0, i_y[i], d00m1[i], i_s[i], e00m1[i], i_r[i], f00m1[i]); // checked
    // type_4 = a + b*c + d*e + f*g
    h[6][i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d10m1[i], r_p_dot_t[i], f00m1[i], i_q_dot_t, e10m1[i], -i_p_dot_t, h000m1[i]); // checked
    h[8][i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d01m1[i], s_q_dot_t[i], e00m1[i], i_p_dot_t, f10m1[i], -i_q_dot_t, h000m1[i]); // checked
    // type_2 = a + b*c
    h[7][i] = integral_type_2(-e10m1[i], -i_q_dot_t, h[6][i]); // checked
    // type_3 = a + b*c + d*e
    h[9] [i] = integral_type_3(h000m1[i], -i_p_dot_t, h[6][i], -i_r[i], f00m1[i]); // checked
    h[10][i] = integral_type_3(h000m1[i], -i_q_dot_t, h[8][i], -i_s[i], e00m1[i]); // checked
    h[11][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_3(-d00m3[i], i_p_dot_t, f00m3[i], i_q_dot_t, e00m3[i]); // checked
    // type_2 = a + b*c
    h[12][i] = integral_type_2(-third * f00m3[i], -i_p_dot_t, h[11][i]); // checked
    h[13][i] = integral_type_2(-third * e00m3[i], -i_q_dot_t, h[11][i]); // checked
    // type_4 = a + b*c + d*e + f*g
    h[14][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d10m3[i], r_p_dot_t[i], f00m3[i], i_q_dot_t, e10m3[i], -i_p_dot_t, 0.0); // checked
    h[15][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d01m3[i], s_q_dot_t[i], e00m3[i], i_p_dot_t, f10m3[i], -i_q_dot_t, 0.0); // checked
    // type_2 = a + b*c
    h[16][i] = integral_type_2(-third * f10m3[i], -i_p_dot_t, h[15][i]); // checked
    // type_3 = a + b*c + d*e
    h[17][i] = integral_type_3(third * 0.0, -third*i_r[i], f00m3[i], -i_p_dot_t, h[14][i]); // checked
    h[18][i] = integral_type_3(third * 0.0, -third*i_s[i], e00m3[i], -i_q_dot_t, h[15][i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    h[19][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_5(-d11m3[i], r_p_dot_t[i], f10m3[i], -i_p_dot_t, h[2][i], s_q_dot_t[i], e10m3[i], -i_q_dot_t, h[1][i]); // checked
    // type_3 = a + b*c + d*e
    h[20][i] = integral_type_3(third * h[2][i], -third*i_r[i], f10m3[i], -i_p_dot_t, h[19][i]); // checked
    h[21][i] = integral_type_3(third * h[1][i], -third*i_s[i], e10m3[i], -i_q_dot_t, h[19][i]); // checked
    // type_4 = a + b*c + d*e + f*g
    h[23][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d20m3[i], -two_p_dot_t, h[1][i], i_q_dot_t, e20m3[i], i_p_dot_t*r_sq[i], f00m3[i]); // checked
    h[22][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_4(-d02m3[i], -two_q_dot_t, h[2][i], i_p_dot_t, f20m3[i], i_q_dot_t*s_sq[i], e00m3[i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    h[25][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_5(h[2][i], s_q_dot_t[i], e01m3[i], -i_y[i], d01m3[i], i_p_dot_t, f11m3[i], -i_q_dot_t, h[0][i]); // checked
    h[24][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_5(h[1][i], r_p_dot_t[i], f01m3[i], -i_y[i], d10m3[i], i_q_dot_t, e11m3[i], -i_p_dot_t, h[0][i]); // checked
    // type_9 = a + b*c + d*e + f*g + h*i + j*k
    h[26][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_9(h[7][i], r_p_dot_t[i], f11m3[i], -i_y[i], d11m3[i], s_q_dot_t[i], e11m3[i], -i_p_dot_t, h[8][i], -i_q_dot_t, h[6][i]); // checked
    // type_3 = a + b*c + d*e
    h[27][i] = integral_type_3(third * h[8][i], -i_p_dot_t, h[26][i], -third*i_r[i], f11m3[i]); // checked
    h[28][i] = integral_type_3(third * h[6][i], -i_q_dot_t, h[26][i], -third*i_s[i], e11m3[i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    h[29][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_5(h[9] [i], -two_p_dot_t, h[6][i], i_q_dot_t, e21m3[i], i_p_dot_t*r_sq[i], f01m3[i], -i_y[i], d20m3[i]); // checked
    h[30][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_5(h[10][i], -two_q_dot_t, h[8][i], i_p_dot_t, f21m3[i], i_q_dot_t*s_sq[i], e01m3[i], -i_y[i], d02m3[i]); // checked
    // type_3 = a + b*c + d*e
    h[31][i] = integral_type_3(two_third * h[6][i], -i_p_dot_t, h[29][i], -third*r_sq[i], f01m3[i]); // checked
    h[32][i] = integral_type_3(two_third * h[8][i], -i_q_dot_t, h[30][i], -third*s_sq[i], e01m3[i]); // checked
    // type_10 = a + b*c + d*e + f*g + h*i + j*k
    h[33][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_10(h[5][i], -two_q_dot_t, h[3][i], -i_y[i], d12m3[i], r_p_dot_t[i], f21m3[i], -i_p_dot_t, h021m3[i], i_q_dot_t*s_sq[i], e11m3[i]); // checked
    h[37][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_10(h[4][i], -two_p_dot_t, h[3][i], -i_y[i], d21m3[i], s_q_dot_t[i], e21m3[i], -i_q_dot_t, h201m3[i], i_p_dot_t*r_sq[i], f11m3[i]); // checked
    // type_11 = a + b*c + d*e + f*g + h*i
    h[34][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_11(-d22m3[i], i_p_dot_t*r_sq[i], f20m3[i], i_q_dot_t*s_sq[i], e20m3[i], -i_p_dot_t*2.0, h[5][i], -i_q_dot_t*2.0, h[4][i]);
    // type_12 = a*b + c*d + e*f + g*h + i*j + k*l*m;
    h[36][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_12(2.0*one_m_q_dot_t_sq, h[3][i], -i_p_dot_t, h[4][i], p_dot_t_q_dot_t, h201m3[i], y_p_dot_t[i], d21m3[i], -p_dot_t_q_dot_t*i_s[i], e21m3[i], -one_m_q_dot_t_sq, r_sq[i], f11m3[i]); // checked
    h[35][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * integral_type_12(2.0*one_m_p_dot_t_sq, h[3][i], -i_q_dot_t, h[5][i], p_dot_t_q_dot_t, h021m3[i], y_q_dot_t[i], d12m3[i], -p_dot_t_q_dot_t*i_r[i], f21m3[i], -one_m_p_dot_t_sq, s_sq[i], e11m3[i]); // checked
  }
  // Evaluating the integrals.
  for (int i = 0; i<i_num_integrals; i++){
    o_sch[i] = dot_product(h[i], signv, n_limits);
  }
  return o_sch;
}

void compute_forces_linear_rectangle(float *i_sch, float i_vec_int[][3],
                                     float *i_rp, float *i_sp, float i_factor,
                                     float *o_nodal_force[n_nodes], float *o_total_force){
  // Calculating nodal forces
  // x3.
  vertex_force_linear_rectangle(i_sch, i_vec_int, i_rp[1], i_sp[1], i_rp[1]*i_sp[1],  i_factor, o_nodal_force[0]);
  // x4.
  vertex_force_linear_rectangle(i_sch, i_vec_int, i_rp[0], i_sp[1], i_rp[0]*i_sp[1], -i_factor, o_nodal_force[1]);
  // x5.
  vertex_force_linear_rectangle(i_sch, i_vec_int, i_rp[1], i_sp[0], i_rp[1]*i_sp[0], -i_factor, o_nodal_force[2]);
  // x5
  vertex_force_linear_rectangle(i_sch, i_vec_int, i_rp[0], i_sp[0], i_rp[0]*i_sp[0],  i_factor, o_nodal_force[3]);
  for (int i = 0; i < n_nodes; i++){
    o_total_force[0] += o_nodal_force[i][0];
    o_total_force[1] += o_nodal_force[i][1];
    o_total_force[2] += o_nodal_force[i][2];
  }
}

void init_force(float *nodal_force[n_nodes], float *total_force){
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

void add_force(float *p_nodal_force[n_nodes], float *p_total_force, float *nodal_force[n_nodes], float *total_force){
  // Adds forces for averaging purposes later on.
  for (int i = 0; i < n_nodes; i++){
    for (int j = 0; j < 3; j++){
      nodal_force[i][j] += p_nodal_force[i][j];
    }
  }
  for (int i = 0; i < 3; i++){
    total_force[i] += p_total_force[i];
  }
}

void mean_force(float *nodal_force[n_nodes], float *total_force, int n_samples){
  // Sets forces to zero.
  for (int i = 0; i < n_nodes; i++){
    for (int j = 0; j < 3; j++){
      nodal_force[i][j] = nodal_force[i][j]/n_samples;
    }
  }
  for (int i = 0; i < 3; i++){
    total_force[i] = total_force[i]/n_samples;
  }
}

void nodal_surface_force_linear_rectangle(float *x1, float *x2, float *x3, float *x4, float *x5, float *x6, float *b, float *p, float *q, float *n, float p_norm, float q_norm, float mu, float nu, float a, float a_sq, float one_m_nu, float factor, float *nodal_force[n_nodes], float *total_force){
  // Vertices of the dislocation segment and a linear rectangular surface.
  /*
          x5 -------------------- x6
             |                  |                 t
             |                  |               ---->           b
             |                  |           x1 -------- x2    ^>
             |                  |             /        |     /
      ^      |                  |            /         \
    q |      |                  |          ...        ...
      |   X3 -------------------- X4

          ---->
            p
  */
  // Characteristic vectors.
  float t[3];
  // Basis vectors (unitary).
  float p_x_t[3], q_x_t[3];
  // Limits of the distance vector from the plane to dislocation line segment.
  // r_lim[0][] = vector from x1 to x3, r_lim[1][] = vector from x2 to x6.
  float r_lim[2][3];
  // r, s limits
  float rp[2], sp[2];
  //  y, r, s coordinates
  float y[n_limits], r[n_limits], s[n_limits];
  // Vectors for the integrals.
  float vec_int[11][3];
  // Scalar value of the integrals.
  float sch[38];
  //  Auxiliary constants to reduce computational cost.
  //  p dot t, q dot t
  float p_dot_t, q_dot_t;
  // Set forces to zero.
  init_force(nodal_force, total_force);
  // Build unit vectors t.
  init_vector (x1, x2, 3, t);
  // Dot products.
  p_dot_t = dot_product(p, t, 3);
  q_dot_t = dot_product(q, t, 3);
  //*******************WARNING*******************//
  // This formulation assumes x3-x6 is diagonal! //
  //*******************WARNING*******************//
  cross_product(p, t, p_x_t);
  cross_product(q, t, q_x_t);
  // Vectors between x3 and x1, and x6 and x2.
  for (int i = 0; i < 3; i++){
    r_lim[0][i] = x3[i] - x1[i];
    r_lim[1][i] = x6[i] - x2[i];
  }
  // Integral bounds for y, r, s.
  for (int i = 0; i < 2; i++){
    rp[i] = init_point(r_lim[i], q_x_t, p, q_x_t, 3);
    sp[i] = init_point(r_lim[i], p_x_t, q, p_x_t, 3);
  }
  // Assign coordinates for the evaluation of the integrals.
  y[0] = y[2] = y[4] = y[6] = init_point(r_lim[1], n    , t, n    , 3);
  y[1] = y[3] = y[5] = y[7] = init_point(r_lim[0], n    , t, n    , 3);
  r[0] = r[1] = r[2] = r[3] = rp[1];
  r[4] = r[5] = r[6] = r[7] = rp[0];
  s[0] = s[1] = s[4] = s[5] = sp[1];
  s[2] = s[3] = s[6] = s[7] = sp[0];
  // Calculate vectors for integrals.
  integral_vector(p, q, b, t, n, one_m_nu, a_sq, vec_int);
  // Calculate integrals.
  integrals_linear_rectangle(r, s, y, p, q, t, p_dot_t, q_dot_t, a_sq, sch, 38);
  // Calculate nodal forces.
  compute_forces_linear_rectangle(sch, vec_int, rp, sp, factor, nodal_force, total_force);
  //printf("total_force[x, y, z] = [%2.14f, %2.14f, %2.14f]\n", total_force[0], total_force[1], total_force[2]);
}

void nodal_surface_force_linear_rectangle_special(float *x1, float *x2, float *x3, float *x4, float *x5, float *x6, float *b, float *t, float *p, float *q, float *n, float p_norm, float q_norm, float mu, float nu, float a, float a_sq, float one_m_nu, float factor, float *nodal_force[n_nodes], float *total_force){
  /*
    Forces
    nodal_force[0][] = F_x3[x, y, z], nodal_force[1][] = F_x4[x, y, z],
    nodal_force[2][] = F_x5[x, y, z], nodal_force[3][] = F_x6[x, y, z]
    total_force[x, y, z] = F_x3[x, y, z] + F_x4[x, y, z] + F_x5[x, y, z] + F_x6[x, y, z]
  */
  // Modulus of p and q.
  float rot_centre[3], rot_x1[3], rot_x2[3];
  float t_x_n[3], mag_t_x_n, p_total_force[3], *p_nodal_force[n_nodes];
  int rotation;
  rotation = 3;
  // Initialise force to zero.
  init_force(nodal_force, total_force);

  for (int i = 0; i < n_nodes; i++){
    p_nodal_force[i] = (float *) malloc(3 * sizeof(float));
  }
  float angle = 0.01*pi/180.;
  int j;
  cross_product(t, n, t_x_n);
  mag_t_x_n = sqrt(dot_product(t_x_n, t_x_n, 3));
  for (int i = 0; i < 3; i++){
    // Halfway between x1 and x2. x1 + (x2-x1)/2
    rot_centre[i] = 0.5*(x1[i] + x2[i]);
    t_x_n[i] = t_x_n[i]/mag_t_x_n;
  }
  //FILE *fp;
  //fp = fopen("./tests/test2.txt", "w");
  for (int i = 0; i < rotation; i++){
    j = i+1;
    arbitrary_rotation_matrix_3d(j*angle, rot_centre, t_x_n, x1, rot_x1);
    arbitrary_rotation_matrix_3d(j*angle, rot_centre, t_x_n, x2, rot_x2);
    nodal_surface_force_linear_rectangle(rot_x1, rot_x2, x3, x4, x5, x6, b, p, q, n, p_norm, q_norm, mu, nu, a, a_sq, one_m_nu, factor, p_nodal_force, p_total_force);
    //printf("theta  = %f\nrot_x1 = [%f, %f, %f]\nrot_x2 = [%f, %f, %f]\nx1 = [%f, %f, %f]\nx2 = [%f, %f, %f]\n", j*angle, rot_x1[0], rot_x1[1], rot_x1[2], rot_x2[0], rot_x2[1], rot_x2[2], x1[0], x1[1], x1[2], x2[0], x2[1], x2[2]);
    //fprintf(fp, "%3f %2.14f %2.14f %2.14f\n", j*angle*180./pi, p_total_force[0], p_total_force[1], p_total_force[2]);
    add_force(p_nodal_force, p_total_force, nodal_force, total_force);
    arbitrary_rotation_matrix_3d(-j*angle, rot_centre, t_x_n, x1, rot_x1);
    arbitrary_rotation_matrix_3d(-j*angle, rot_centre, t_x_n, x2, rot_x2);
    nodal_surface_force_linear_rectangle(rot_x1, rot_x2, x3, x4, x5, x6, b, p, q, n, p_norm, q_norm, mu, nu, a, a_sq, one_m_nu, factor, p_nodal_force, p_total_force);
    //printf("theta  = %f\nrot_x1 = [%f, %f, %f]\nrot_x2 = [%f, %f, %f]\nx1 = [%f, %f, %f]\nx2 = [%f, %f, %f]\n", j*angle, rot_x1[0], rot_x1[1], rot_x1[2], rot_x2[0], rot_x2[1], rot_x2[2], x1[0], x1[1], x1[2], x2[0], x2[1], x2[2]);
    //fprintf(fp, "%3.6f %2.14f %2.14f %2.14f\n", -j*angle*180./pi, p_total_force[0], p_total_force[1], p_total_force[2]);
    add_force(p_nodal_force, p_total_force, nodal_force, total_force);
  }
  //fclose(fp);
  mean_force(nodal_force, total_force, rotation*2);
  for (int i = 0; i < n_nodes; i++){
    free(p_nodal_force[i]);
  }
  //printf("total_force[x, y, z] = [%2.14f, %2.14f, %2.14f]\n", total_force[0], total_force[1], total_force[2]);
}

// Parallel functions.
float *element_host_device_map(float *i_node_arr[], int i_n_elem_scope, int i_n_nodes){
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
  float *o_g_elem_arr;
  int idxt = 0, idxf = 0, idxi = 0;
  // Allocate a 1D output array of length 3*E*N.
  o_g_elem_arr = (float *) malloc(3 * n_elem_n_nodes * sizeof(float));
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

float *se_host_device_map(float *i_x3_arr, float *i_x4_arr, float *i_x5_arr, float *i_x6_arr, int i_n_se){
  /*
    Maps E surface elements with 4 nodes.
    This ensures coalesced memory accesses when parallelising over surface elements.
    i_xn_arr[3e + 0:2] = [x_en, y_en, z_en]

      to

    o_g_se_arr =
      [x_00, x_10, ..., x_(E-1)0, x_01, x_11, ..., x_(E-1)1, ..., x_0N, x_1(N-1), ..., x_(E-1)3,
       y_00, y_10, ..., y_(E-1)0, y_01, y_11, ..., y_(E-1)1, ..., y_0N, y_1(N-1), ..., y_(E-1)3,
       z_00, z_10, ..., z_(E-1)0, z_01, z_11, ..., z_(E-1)1, ..., z_0N, z_1(N-1), ..., z_(E-1)3]
    i_n_se    = E
    n_elem_n_nodes  = 4*E
    idxi = input (initial) index
    idxf = output (final) index
  */
  /*
  int const n_se_n_nodes = i_n_se*4;
  float *o_g_se_arr;
  int idxi = 0, idxf = 0;
  // Allocate a 1D output array of length 3*4*E.
  o_g_se_arr = (float *) malloc(3 * n_se_n_nodes * sizeof(float));
  // Loop over coordinates.
  for (int i = 0; i < 3; i++){
    // Reset the input index to point at the i'th coordinate of the first node of the j'th element.
    idxi = i;
    // Loop over surface elements.
    for (int j = 0; j < i_n_se; j++){
      // Displace the ouptut index to point at the next node of the j'th element.
      o_g_se_arr[idxf + j           ] = i_x3_arr[idxi];
      o_g_se_arr[idxf + j +   i_n_se] = i_x4_arr[idxi];
      o_g_se_arr[idxf + j + 2*i_n_se] = i_x5_arr[idxi];
      o_g_se_arr[idxf + j + 3*i_n_se] = i_x6_arr[idxi];
      // Advance the input index to point at the i'th coordinate of the first node of the (j+1)'th element.
      idxi += 3;
    }
    // Advance the output index to point at the (i+1)'th coordinate of the first node of the first element.
    idxf += n_se_n_nodes;
  }
  return o_g_se_arr;
  */
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
  float *o_g_se_arr;
  int idxi = 0, idxf = 0;
  // Allocate a 1D output array of length 3*2*E.
  o_g_se_arr = (float *) malloc(3 * n_se_nodes * sizeof(float));
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

void fx_device_host_map(float *i_g_fx_arr, float *o_fx_arr[], int i_n_se, int i_n_nodes){
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
      }
      // Advance the output index to point at the first coordinate of the (j+1)'th element.
      idxf += 3;
    }
    // Advance the input index to point at the first coordinate of the (i+1)'th node of the first element.
    idxi += i_n_se;
  }
}

void ftot_device_host_map(float *i_g_ftot_arr, float *o_ftot_arr, int i_n_se){
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

float *dln_host_device_map(float *i_x1_arr, float *i_x2_arr, int i_n_dln){
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
  float *o_g_dln_arr;
  int idxi = 0, idxf = 0;
  // Allocate a 1D output array of length 3*2*E.
  o_g_dln_arr = (float *) malloc(3 * n_dln_nodes * sizeof(float));
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

float *b_host_device_map(float *i_b_arr, int i_n_b){
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
  float *o_g_b_arr;
  int idx = 0;
  // Allocate return array.
  o_g_b_arr = (float *) malloc(i_n_b * 3 * sizeof(float));
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

#if debug
  void print_fmap(float *i_x_arr, int i_n_elements,
                  int i_n_nodes, int i_node_label){
    // Print nodes for forward map.
    int const n_elements_n_nodes = i_n_elements * i_n_nodes;
    // Save initial node label.
    int init_node_label = i_node_label;
    // Loop over all coordinates of all nodes in the scope.
    for (int i = 0; i < n_elements_n_nodes * 3; i++){
      // If we're onto a new block of nodes.
      if (i%i_n_elements == 0){
        // If we're onto a new coordinate block and we're not at the first entry.
        if (i%n_elements_n_nodes == 0 && i!=0){
          // Skip a line.
          printf("\n");
          // Reset the node label back to the first node of the element.
          i_node_label = init_node_label;
        }
        // Print the node label and skip a line.
        printf("X%d\n", i_node_label);
        // Increse the node label.
        i_node_label++;
      }
      // Print the array entry.
      printf("x_arr[%d] = %f\n", i, i_x_arr[i]);
    }
  }
#endif

__device__ void cuda_init_force(float nodal_force[][3], float *total_force){
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

__device__ void se_device_thread_map(float *i_g_se_arr,
                                     float *o_x3, float *o_x4, float *o_x5, float *o_x6,
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

__device__ void add_force_thread_device(float i_nodal_force[][3], float *i_total_force, float *o_g_fx_arr, float *o_g_ftot_arr, int i_n_se, int i_n_nodes, int idx){
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
      atomicAdd(&o_g_fx_arr[idxf + j*n_se_n_nodes], i_nodal_force[i][j]);
    }
    // Advance the output index to point at the (i+1)'th node of the first coordinate of the idx'th surface element.
    idxf += i_n_se;
  }
  // Total force per surface element.
  // Loop over coordinates.
  idxf = idx;
  for (int i= 0; i < 3; i++){
    atomicAdd(&o_g_ftot_arr[idxf], i_total_force[i]);
    // Advance the output index to point at the (i+1)'th coordinate of the idx'th surface element.
    idxf += i_n_se;
  }
}

__device__ void dln_device_thread_map(float *i_g_dln_arr,
                                      float *o_x1, float *o_x2,
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

__device__ void b_device_thread_map(float *i_g_b_arr,
                                    float *o_b,
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

__device__ float cuda_dot_product(float *i_vec1, float *i_vec2, int i_vec_size){
  // Returns the dot product of i_vec1, i_vec2.
  float result = 0.0;
  for (int i = 0; i < i_vec_size; i++){
    result += i_vec1[i]*i_vec2[i];
  }
  return result;
}

__device__ float *cuda_cross_product(float *i_vec1, float *i_vec2,
                      float *o_vec){
  // Returns the cross product of i_vec1 x i_vec2.
  o_vec[0] = i_vec1[1]*i_vec2[2] - i_vec1[2]*i_vec2[1];
  o_vec[1] = i_vec1[2]*i_vec2[0] - i_vec1[0]*i_vec2[2];
  o_vec[2] = i_vec1[0]*i_vec2[1] - i_vec1[1]*i_vec2[0];
  return o_vec;
}

__device__ float *cuda_cross_product2(float *i_vec1, float *i_vec2){
  float *o_vec;
  o_vec = (float *) malloc(3*sizeof(float));
  // Returns the cross product of i_vec1 x i_vec2.
  o_vec[0] = i_vec1[1]*i_vec2[2] - i_vec1[2]*i_vec2[1];
  o_vec[1] = i_vec1[2]*i_vec2[0] - i_vec1[0]*i_vec2[2];
  o_vec[2] = i_vec1[0]*i_vec2[1] - i_vec1[1]*i_vec2[0];
  return o_vec;
}

__device__ float *cuda_normalise_vector(float *i_vec, int i_vec_size,
                         float *o_vec){
  // Returns a normalised i_vec to o_vec.
  float mag_vec = 0.0;
  mag_vec = cuda_dot_product(i_vec, i_vec, i_vec_size);
  // Check magnitude is not zero.
  if (mag_vec == 0.0){
    //printf("ERROR: nodal_surface_force_linear_rectangle: normalise_vector: mag_vec = 0: A vector cannot have magnitude 0\n");
    asm("trap;");
  }
  mag_vec = sqrt(mag_vec);
  for(int i = 0; i < i_vec_size; i++){
    o_vec[i] = i_vec[i]/mag_vec;
  }
  return o_vec;
}

__device__ void cuda_normalise_vector2(float *i_vec, int i_vec_size,
                       float *o_vec, float *o_mag_vec){
  // Returns a normalised i_vec to o_vec, and the magnitude of the vector in o_mag_vec.
  // Has to be a void function in order to 'return' two values, the magnitude should be passed as a reference eg:
  /*
    int    size = 3;
    float i_vec[size], o_vec[size], magnitude;
    normalise_vector2(i_vec, size, o_vec, &magnitude);
  */
  *o_mag_vec = cuda_dot_product(i_vec, i_vec, i_vec_size);
  // Check magnitude is not zero.
  if (*o_mag_vec == 0.0){
    //printf("ERROR: nodal_surface_force_linear_rectangle: normalise_vector2: o_mag_vec = 0: A vector cannot have magnitude 0\n");
    asm("trap;");
  }
  *o_mag_vec = sqrt(*o_mag_vec);
  for(int i = 0; i < i_vec_size; i++){
    o_vec[i] = i_vec[i]/ *o_mag_vec;
  }
}

__device__ float *cuda_build_vector(float *i_x1, float *i_x2, int i_vec_size,
                     float *o_vec){
  // Returns a vector o_vec which translates the point i_x1 to i_x2.
  for (int i = 0; i < i_vec_size; i++){
    o_vec[i] = i_x2[i] - i_x1[i];
  }
  return o_vec;
}

__device__ float *cuda_init_vector(float *i_x1, float *i_x2, int i_vec_size,
                    float *io_vec){
  // Builds and returns a normalised vector io_vect.
  cuda_build_vector(i_x1, i_x2, i_vec_size, io_vec);
  cuda_normalise_vector(io_vec, i_vec_size, io_vec);
  return io_vec;
}

__device__ void cuda_init_vector2(float *i_x1  , float *i_x2, int i_vec_size,
                  float *io_vec, float *o_mag_vec){
  // Builds and returns a normalised vector io_vect, and the magnitude of the vector in o_mag_vec.
  // Has to be a void function in order to 'return' two values, the magnitude should be passed as a reference eg:
  /*
    int    size = 3;
    float i_x1[size], i_x2[size], o_vec[size], magnitude;
    init_vector2(i_x1, i_x2, size, o_vec, &magnitude);
  */
  cuda_build_vector(i_x1, i_x2, i_vec_size, io_vec);
  cuda_normalise_vector2(io_vec, i_vec_size, io_vec, o_mag_vec);
}

__device__ float cuda_init_point(float *i_vec1, float *i_vec2,
                  float *i_vec3, float *i_vec4,
                  int     i_vec_size){
  // Initialises points on the surface element given four vectors.
  float result = 0.0;
  float denom = 0.0;
  denom = cuda_dot_product(i_vec3, i_vec4, i_vec_size);
  if (denom == 0.0){
    //printf("nodal_surface_force_linear_rectangle: init_point: division by zero, dot_product(i_vec3, i_vec4) = %2.14f\n", denom);
    asm("trap;");
  }
  result = cuda_dot_product(i_vec1, i_vec2, i_vec_size)/denom;
  return result;
}

__device__ void cuda_integral_vector(float *i_p, float *i_q, float *i_b, float *i_t, float *i_n,
                     float i_one_m_nu, float i_a_sq,
                     float o_vec_int[][3]){
  // Calculates the vectors that are multiplied by the integrals.
  // Has to be void in order to return a 2D array, otherwise we'd have to use dynamic memory allocation and therefore pointers. We don't need such flexibility as this is quite specific.
  /*
    (t x b), (p x b)
    (q x b), (b x t)
    t*((p x b) dot n), t*((q x b) dot n)
    t*((t x b) dot n), t*((b x t) dot n)
    (t dot n) * (p x b), (t dot n) * (q x b)
    (t dot n) * (t x b)
    n * ((p x b) dot t), n * ((q x b) dot t)
  */
  float t_x_b[3], p_x_b[3],
         q_x_b[3], b_x_t[3],
         t_p_x_b_dot_n[3], t_q_x_b_dot_n[3],
         t_t_x_b_dot_n[3], t_b_x_t_dot_n[3],
         t_dot_n_p_x_b[3], t_dot_n_q_x_b[3],
         t_dot_n_t_x_b[3],
         n_p_x_b_dot_t[3], n_q_x_b_dot_t[3];
  /*
    t dot n
    ((p x b) dot n), ((q x b) dot n)
    ((t x b) dot n), ((b x t) dot n)
    ((p x b) dot t), ((q x b) dot t)
    (t dot n) * ((p x b) dot t), (t dot n) * ((q x b) dot t)
    3*(t dot n) * ((p x b) dot t), 3*(t dot n) * ((q x b) dot t)
    1.5 * (1-nu) * a^2, a^3
  */
  float t_dot_n,
         p_x_b_dot_n, q_x_b_dot_n,
         t_x_b_dot_n, b_x_t_dot_n,
         p_x_b_dot_t, q_x_b_dot_t,
         t_dot_n_p_x_b_dot_t, t_dot_n_q_x_b_dot_t,
         t_dot_n_p_x_b_dot_t_3, t_dot_n_q_x_b_dot_t_3,
         one_m_nu_1p5_a_sq, a_sq_3;

  // Calculate auxiliary local variables.
  one_m_nu_1p5_a_sq = 1.5 * i_one_m_nu * i_a_sq;
  a_sq_3 = 3.0 * i_a_sq;
  cuda_cross_product(i_p, i_b, p_x_b);
  cuda_cross_product(i_q, i_b, q_x_b);
  cuda_cross_product(i_t, i_b, t_x_b);
  cuda_cross_product(i_b, i_t, b_x_t);
  t_dot_n     = cuda_dot_product(i_t  , i_n, 3);
  p_x_b_dot_n = cuda_dot_product(p_x_b, i_n, 3);
  q_x_b_dot_n = cuda_dot_product(q_x_b, i_n, 3);
  t_x_b_dot_n = cuda_dot_product(t_x_b, i_n, 3);
  b_x_t_dot_n = cuda_dot_product(b_x_t, i_n, 3);
  p_x_b_dot_t = cuda_dot_product(p_x_b, i_t, 3);
  q_x_b_dot_t = cuda_dot_product(q_x_b, i_t, 3);
  t_dot_n_p_x_b_dot_t   = t_dot_n * p_x_b_dot_t;
  t_dot_n_q_x_b_dot_t   = t_dot_n * q_x_b_dot_t;
  t_dot_n_p_x_b_dot_t_3 = 3.0 * t_dot_n_p_x_b_dot_t;
  t_dot_n_q_x_b_dot_t_3 = 3.0 * t_dot_n_q_x_b_dot_t;

  for (int i=0; i<3; i++){
      // Preparing common array factors.
      t_t_x_b_dot_n[i] = i_t  [i]*t_x_b_dot_n;
      t_b_x_t_dot_n[i] = i_t  [i]*b_x_t_dot_n;
      t_p_x_b_dot_n[i] = i_t  [i]*p_x_b_dot_n;
      t_q_x_b_dot_n[i] = i_t  [i]*q_x_b_dot_n;
      t_dot_n_t_x_b[i] = t_x_b[i]*t_dot_n;
      t_dot_n_q_x_b[i] = q_x_b[i]*t_dot_n;
      t_dot_n_p_x_b[i] = p_x_b[i]*t_dot_n;
      n_q_x_b_dot_t[i] = q_x_b_dot_t*i_n[i];
      n_p_x_b_dot_t[i] = p_x_b_dot_t*i_n[i];
      // Calculating vectors.
      o_vec_int [0][i] = (t_dot_n_t_x_b[i] + t_t_x_b_dot_n[i])*i_one_m_nu        + b_x_t[i]*t_dot_n + t_b_x_t_dot_n[i]; // checked
      o_vec_int [1][i] = (t_dot_n_q_x_b[i] + t_q_x_b_dot_n[i])*i_one_m_nu        - n_q_x_b_dot_t[i] + i_q[i]*b_x_t_dot_n; // checked
      o_vec_int [2][i] = (t_dot_n_p_x_b[i] + t_p_x_b_dot_n[i])*i_one_m_nu        - n_p_x_b_dot_t[i] + i_p[i]*b_x_t_dot_n; // checked
      o_vec_int [3][i] = (t_dot_n_t_x_b[i] + t_t_x_b_dot_n[i])*one_m_nu_1p5_a_sq; // checked
      o_vec_int [4][i] = (t_dot_n_q_x_b[i] + t_q_x_b_dot_n[i])*one_m_nu_1p5_a_sq - n_q_x_b_dot_t[i]*a_sq_3; // checked
      o_vec_int [5][i] = (t_dot_n_p_x_b[i] + t_p_x_b_dot_n[i])*one_m_nu_1p5_a_sq - n_p_x_b_dot_t[i]*a_sq_3; // checked
      o_vec_int [6][i] = - t_dot_n_q_x_b_dot_t_3*i_q[i]; // checked
      o_vec_int [7][i] = - t_dot_n_p_x_b_dot_t_3*i_p[i]; // checked
      o_vec_int [8][i] = - t_dot_n_q_x_b_dot_t_3*i_t[i]; // checked
      o_vec_int [9][i] = - t_dot_n_p_x_b_dot_t_3*i_t[i]; // checked
      o_vec_int[10][i] = - t_dot_n_p_x_b_dot_t_3*i_q[i] - t_dot_n_q_x_b_dot_t_3*i_p[i]; // checked
  }
}

__device__ float cuda_seed_single_integral(float a, float b){
  // Single integral seed.
  float result = 0.0;
  float arg;
  arg = a+b;
  if (arg <= 0.0){
    //printf("nodal_surface_force_linear_rectangle: seed_single_integral: log(%2.14f) = %2.14f\n", arg, log(arg));
    asm("trap;");
  }
  result = log(arg);
  return result;
}

__device__ float cuda_integral_type_1(float a, float b,
                       float c, float d,
                       float e){
  float result = 0.0;
  result = 0.5*(a*b + (c-d)*e);
  return result;
}

__device__ float cuda_numer_seed_double_integral(float a,
                                  float b, float c, float d,
                                  float e){
  // Returns the argument of the numerator of the seed integral for float integrals.
  float result = 0.0;
  result = a*(b - c + d) + e;
  return result;
}

__device__ float cuda_denom_seed_double_integral(float a, float b,
                                  float c, float d){
  // Returns the argument of the denominator of the seed integral for float integrals.
  float result = 0.0;
  result = a*b + c*d;
  if(result == 0.0){
    //printf("nodal_surface_force_linear_rectangle: denom_seed_integral_00m3: division by 0, result = %2.14f\n", result);
    asm("trap;");
  }
  return result;
}

__device__ float cuda_seed_double_integral(float a, float b, float c, float d, float e,
                            float f, float g, float h, float i){
  // Returns the seed integral for float integrals.
  float result      = 0.0;
  float numerator   = 0.0;
  float denominator = 0.0;
  numerator   = cuda_numer_seed_double_integral(a, b, c, d, e);
  denominator = cuda_denom_seed_double_integral(f, g, h, i);
  if(denominator > 0.){
    denominator = sqrt(denominator);
    result = 2.0/denominator * atan(numerator/denominator);
  }
  else{
    denominator = sqrt(-denominator);
    result = -2.0/denominator * atanh(numerator/denominator);
  }
  return result;
}

__device__ float cuda_integral_type_2(float a,
                       float b, float c){
  float result = 0.0;
  result = a + b*c;
  return result;
}

__device__ float cuda_integral_type_3(float a,
                       float b, float c,
                       float d, float e){
  float result = 0.0;
  result = a + b*c + d*e;
  return result;
}

__device__ float cuda_integral_type_4(float a,
                       float b, float c,
                       float d, float e,
                       float f, float g){
  float result = 0.0;
  result = a + b*c + d*e + f*g;
  return result;
}

__device__ float cuda_integral_type_5(float a,
                       float b, float c,
                       float d, float e,
                       float f, float g,
                       float h, float i){
  float result = 0.0;
  result = a + b*c + d*e + f*g + h*i;
  return result;
}

__device__ float cuda_integral_type_6(float a, float b,
                       float c, float d,
                       float e, float f, float g,
                       float h, float i,
                       float j, float k){
  float result = 0.0;
  result = a*b + c*d + (e+f)*g + h*i + j*k;
  return result;
}

__device__ float cuda_integral_type_7(float a, float b,
                       float c, float d,
                       float e, float f, float g,
                       float h, float i){
  float result = 0.0;
  result = a*b + c*d + (e+f)*g + h*i;
  return result;
}

__device__ float cuda_integral_type_8(float a, float b,
                       float c, float d,
                       float e, float f,
                       float g, float h){
  float result = 0.0;
  result = a*b + c*d + e*f + g*h;
  return result;
}

__device__ float cuda_integral_type_9(float a,
                       float b, float c,
                       float d, float e,
                       float f, float g,
                       float h, float i,
                       float j, float k){
  float result = 0.0;
  result = a + b*c + d*e + f*g + h*i + j*k;
  return result;
}

__device__ float cuda_integral_type_10(float a,
                        float b, float c,
                        float d, float e,
                        float f, float g,
                        float h, float i,
                        float j, float k){
  float result = 0.0;
  result = a + b*c + d*e + f*g + h*i + j*k;
  return result;
}

__device__ float cuda_integral_type_11(float a,
                        float b, float c,
                        float d, float e,
                        float f, float g,
                        float h, float i){
  float result = 0.0;
  result = a + b*c + d*e + f*g + h*i;
  return result;
}

__device__ float cuda_integral_type_12(float a, float b,
                        float c, float d,
                        float e, float f,
                        float g, float h,
                        float i, float j,
                        float k, float l, float m){
  float result = 0.0;
  result = a*b + c*d + e*f + g*h + i*j + k*l*m;
  return result;
}

__device__ float *cuda_integrals_linear_rectangle(float *i_r, float *i_s, float *i_y,
                                   float *i_p, float *i_q, float *i_t,
                                   float i_p_dot_t, float i_q_dot_t, float i_a_sq,
                                   float *o_sch){
  // Sign vector for quick evaluation of integrals via the dot product.
  float signv[8] = {1.0,-1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0};
  const float third = 1.0/3.0;
  const float two_third = 2.0/3.0;
  // Integrals.
  float // Single integrals.
         a0m1[n_limits], b0m1[n_limits], c0m1[n_limits], a1m1[n_limits], b1m1[n_limits], c1m1[n_limits],
         a01 [n_limits], b01 [n_limits], c01 [n_limits], a11 [n_limits], b11 [n_limits],
         // float integrals.
         d00m3[n_limits], e00m3[n_limits], f00m3[n_limits], d01m3[n_limits], d10m3[n_limits], d11m3[n_limits], e01m3[n_limits], e10m3[n_limits],
         e11m3[n_limits], f01m3[n_limits], f10m3[n_limits], f11m3[n_limits], d00m1[n_limits], e00m1[n_limits], f00m1[n_limits], d01m1[n_limits],
         e01m1[n_limits], f01m1[n_limits], d10m1[n_limits], e10m1[n_limits], f10m1[n_limits], d11m1[n_limits], e11m1[n_limits], f11m1[n_limits],
         d02m3[n_limits], d20m3[n_limits], e20m3[n_limits], f20m3[n_limits], d001 [n_limits], e001 [n_limits], f001 [n_limits], d02m1[n_limits],
         d20m1[n_limits], e20m1[n_limits], f20m1[n_limits], d12m3[n_limits], d21m3[n_limits], e21m3[n_limits], f21m3[n_limits], d22m3[n_limits],
         // Triple integrals
         h [38][n_limits], h001m1[n_limits], h010m1[n_limits], h100m1[n_limits], h021m3[n_limits], h201m3[n_limits], h000m1[n_limits];
 /*
   y^2, r^2, s^2 coordinates
   y*(p dot t), y*(q dot t), r*(p dot t), r*(q dot t)
    r dot p   ,  r dot q   ,  r dot t
   (r dot p)^2, (r dot q)^2, (r dot t)^2
 */
 float y_sq[n_limits], r_sq[n_limits], s_sq[n_limits],
        y_p_dot_t[n_limits], y_q_dot_t[n_limits], r_p_dot_t[n_limits], s_q_dot_t[n_limits],
        r_dot_p[n_limits], r_dot_q[n_limits], r_dot_t[n_limits],
        r_dot_p_sq[n_limits], r_dot_q_sq[n_limits], r_dot_t_sq[n_limits];
  /*
   ra = sqrt((r_vec dot r_vec) + a**2) from non-singular dislocation theory, there are 8 r_vec, thus 8 ra's.
   ra_sq = ra^2 (element-wise squaring)
   ra_c_o_3 = 1/3 * ra^3 (element-wise cubing and division)
  */
  float ra[n_limits], ra_sq[n_limits], ra_c_o_3[n_limits];
  // Vector R from x1 to x3, x4, x5, x6 and x2 to x3, x4, x5, x6. 8 vectors with 3 components each.
  float r_vec[n_limits][3];
  /*
   Auxiliary constants to reduce computational cost.
   2*(p dot t)
   2*(q dot t)
   1-(p dot t)
   1-(q dot t)
   1-(p dot t)^2
   1-(q dot t)^2
   1/(1-(p dot t)^2)
   1/(1-(q dot t)^2)
   1-(p dot t)^2-(q dot t)^2
   1/(1-(p dot t)^2-(q dot t)^2)
   1/3 * 1/(1-(p dot t)^2-(q dot t)^2)
   (p dot t)*(q dot t)
  */
  float                         two_p_dot_t,
                                 two_q_dot_t,
                               one_m_p_dot_t,
                               one_m_q_dot_t,
                               one_m_p_dot_t_sq,
                               one_m_q_dot_t_sq,
                         one_o_one_m_p_dot_t_sq,
                         one_o_one_m_q_dot_t_sq,
                  one_m_p_dot_t_sq_m_q_dot_t_sq,
            one_o_one_m_p_dot_t_sq_m_q_dot_t_sq,
        trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq,
        p_dot_t_q_dot_t;
  // Auxiliary scalars.
  two_p_dot_t = 2.0*i_p_dot_t;
  two_q_dot_t = 2.0*i_q_dot_t;
  one_m_p_dot_t = 1.0 - i_p_dot_t;
  one_m_q_dot_t = 1.0 - i_q_dot_t;
  one_m_p_dot_t_sq = 1.0 - i_p_dot_t * i_p_dot_t;
  one_m_q_dot_t_sq = 1.0 - i_q_dot_t * i_q_dot_t;
  p_dot_t_q_dot_t  =       i_q_dot_t * i_p_dot_t;
  one_o_one_m_p_dot_t_sq = 1.0 / one_m_p_dot_t_sq;
  one_o_one_m_q_dot_t_sq = 1.0 / one_m_q_dot_t_sq;
  one_m_p_dot_t_sq_m_q_dot_t_sq = one_m_p_dot_t_sq + one_m_q_dot_t_sq - 1.0;
  one_o_one_m_p_dot_t_sq_m_q_dot_t_sq = 1.0/one_m_p_dot_t_sq_m_q_dot_t_sq;
  trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq = third * one_o_one_m_p_dot_t_sq_m_q_dot_t_sq;
  // Auxiliary arrays.
  // Calculate r for all 8 combinations (from x1 to x3, x4, x5, x6 and x2 to x3, x4, x5, x6).
  for (int i = 0; i < n_limits; i++){
    for (int j = 0; j < 3; j++){
      r_vec[i][j] = i_r[i] * i_p[j] + i_s[i] * i_q[j] + i_y[i] * i_t[j];
    }
  }
  for (int i = 0; i < n_limits; i++){
    r_sq      [i] = i_r[i] * i_r[i];
    s_sq      [i] = i_s[i] * i_s[i];
    y_sq      [i] = i_y[i] * i_y[i];
    ra        [i] = sqrt(r_sq[i]+s_sq[i]+y_sq[i] + two_p_dot_t*i_r[i]*i_y[i] + two_q_dot_t*i_s[i]*i_y[i] + i_a_sq);
    ra_sq     [i] = ra   [i] * ra[i];
    ra_c_o_3  [i] = ra_sq[i] * ra[i]/3.0;
    r_dot_p   [i] = cuda_dot_product(r_vec[i], i_p, 3);
    r_dot_q   [i] = cuda_dot_product(r_vec[i], i_q, 3);
    r_dot_t   [i] = cuda_dot_product(r_vec[i], i_t, 3);
    r_dot_p_sq[i] = r_dot_p[i] * r_dot_p[i];
    r_dot_q_sq[i] = r_dot_q[i] * r_dot_q[i];
    r_dot_t_sq[i] = r_dot_t[i] * r_dot_t[i];
    y_p_dot_t [i] = i_y[i] * i_p_dot_t;
    y_q_dot_t [i] = i_y[i] * i_q_dot_t;
    r_p_dot_t [i] = i_r[i] * i_p_dot_t;
    s_q_dot_t [i] = i_s[i] * i_q_dot_t;
  }
  // Calculate integrals.
  for (int i = 0; i < n_limits; i++){
    // Linear seed integrals.
    a0m1[i] = cuda_seed_single_integral(ra[i], r_dot_p[i]); // checked
    b0m1[i] = cuda_seed_single_integral(ra[i], r_dot_q[i]); // checked
    c0m1[i] = cuda_seed_single_integral(ra[i], r_dot_t[i]); // checked
    // Linear integrals.
    // result = 0.5*(a*b + (c-d)*e)
    a01 [i] = cuda_integral_type_1(ra[i], r_dot_p[i], ra_sq[i], r_dot_p_sq[i], a0m1[i]); // checked
    b01 [i] = cuda_integral_type_1(ra[i], r_dot_q[i], ra_sq[i], r_dot_q_sq[i], b0m1[i]); // checked
    c01 [i] = cuda_integral_type_1(ra[i], r_dot_t[i], ra_sq[i], r_dot_t_sq[i], c0m1[i]); // checked
    // type_2 = a + b*c
    a1m1[i] = cuda_integral_type_2(ra[i]      , y_p_dot_t[i]             , -a0m1[i]); // checked
    b1m1[i] = cuda_integral_type_2(ra[i]      , y_q_dot_t[i]             , -b0m1[i]); // checked
    c1m1[i] = cuda_integral_type_2(ra[i]      , r_p_dot_t[i]+s_q_dot_t[i], -c0m1[i]); // checked
    a11 [i] = cuda_integral_type_2(ra_c_o_3[i], y_p_dot_t[i]             , -a01 [i]); // checked
    b11 [i] = cuda_integral_type_2(ra_c_o_3[i], y_q_dot_t[i]             , -b01 [i]); // checked
    // float seed integrals.
    d00m3[i] = cuda_seed_double_integral(1.0, ra[i], r_dot_p[i], r_dot_q[i], 0.0,
                                    1.0             , i_a_sq, one_m_p_dot_t_sq_m_q_dot_t_sq, y_sq[i]); // checked
    e00m3[i] = cuda_seed_double_integral(one_m_p_dot_t, ra[i], i_r[i], i_y[i], s_q_dot_t[i],
                                    one_m_p_dot_t_sq, i_a_sq, one_m_p_dot_t_sq_m_q_dot_t_sq, s_sq[i]); // checked
    f00m3[i] = cuda_seed_double_integral(one_m_q_dot_t, ra[i], i_s[i], i_y[i], r_p_dot_t[i],
                                    one_m_q_dot_t_sq, i_a_sq, one_m_p_dot_t_sq_m_q_dot_t_sq, r_sq[i]); // checked
    // float integrals.
    // type_2 = a + b*c
    d01m3[i] = cuda_integral_type_2(-a0m1[i], -y_q_dot_t[i], d00m3[i]); // checked
    d10m3[i] = cuda_integral_type_2(-b0m1[i], -y_p_dot_t[i], d00m3[i]); // checked
    // type_3 = a + b*c + d*e
    e01m3[i] = -one_o_one_m_p_dot_t_sq * cuda_integral_type_3(a0m1[i], s_q_dot_t[i], e00m3[i], -i_p_dot_t, c0m1[i]); // checked
    f01m3[i] = -one_o_one_m_q_dot_t_sq * cuda_integral_type_3(b0m1[i], r_p_dot_t[i], f00m3[i], -i_q_dot_t, c0m1[i]); // checked
    // type_2 = a + b*c
    e10m3[i] = cuda_integral_type_2(-c0m1[i], -i_p_dot_t, e01m3[i]); // checked
    f10m3[i] = cuda_integral_type_2(-c0m1[i], -i_q_dot_t, f01m3[i]); // checked
    // type_6 = a*b + c*d + (e+f)*g + h*i + j*k
    d00m1[i] = cuda_integral_type_6(i_r[i], b0m1[i], i_s[i], a0m1[i], i_a_sq, y_sq[i], -d00m3[i], -y_p_dot_t[i], d10m3[i], -y_q_dot_t[i], d01m3[i]); // checked
    // type_7 = a*b + c*d + (e+f)*g + h*i
    e00m1[i] = cuda_integral_type_7(i_r[i], c0m1[i], i_y[i], a0m1[i], i_a_sq, s_sq[i], -e00m3[i], -s_q_dot_t[i], e01m3[i]); // checked
    f00m1[i] = cuda_integral_type_7(i_s[i], c0m1[i], i_y[i], b0m1[i], i_a_sq, r_sq[i], -f00m3[i], -r_p_dot_t[i], f01m3[i]); // checked
    // type_2 = a + b*c
    d11m3[i] = cuda_integral_type_2(-a1m1[i], -y_q_dot_t[i], d10m3[i]); // checked
    // type_4 = a + b*c + d*e + f*g
    e11m3[i] = one_o_one_m_p_dot_t_sq * cuda_integral_type_4(-a1m1[i], r_p_dot_t[i], c0m1[i], -i_p_dot_t, e00m1[i], -s_q_dot_t[i], e10m3[i]); // checked
    f11m3[i] = one_o_one_m_q_dot_t_sq * cuda_integral_type_4(-b1m1[i], s_q_dot_t[i], c0m1[i], -i_q_dot_t, f00m1[i], -r_p_dot_t[i], f10m3[i]); // checked
    // type_3 = a + b*c + d*e
    d20m3[i] = cuda_integral_type_3(d00m1[i], -y_p_dot_t[i], d10m3[i], -i_r[i], b0m1[i]); // checked
    d02m3[i] = cuda_integral_type_3(d00m1[i], -y_q_dot_t[i], d01m3[i], -i_s[i], a0m1[i]); // checked
    e20m3[i] = cuda_integral_type_3(e00m1[i], -i_p_dot_t     , e11m3[i], -i_r[i], c0m1[i]); // checked
    f20m3[i] = cuda_integral_type_3(f00m1[i], -i_q_dot_t     , f11m3[i], -i_s[i], c0m1[i]); // checked
    // type_2 = a + b*c
    d01m1[i] = cuda_integral_type_2(a01[i], -y_q_dot_t[i], d00m1[i]); // checked
    d10m1[i] = cuda_integral_type_2(b01[i], -y_p_dot_t[i], d00m1[i]); // checked
    e01m1[i] = one_o_one_m_p_dot_t_sq * cuda_integral_type_3(a01[i], -s_q_dot_t[i], e00m1[i], -i_p_dot_t, c01[i]); // checked
    f01m1[i] = one_o_one_m_q_dot_t_sq * cuda_integral_type_3(b01[i], -r_p_dot_t[i], f00m1[i], -i_q_dot_t, c01[i]); // checked
    // type_2 = a + b*c
    e10m1[i] = cuda_integral_type_2(c01[i], -i_p_dot_t, e01m1[i]); // checked
    f10m1[i] = cuda_integral_type_2(c01[i], -i_q_dot_t, f01m1[i]); // checked
    // type_6 = a*b + c*d + (e+f)*g + h*i + j*k
    d001[i]  = third * cuda_integral_type_6(i_r[i], b01[i], i_s[i], a01[i], y_sq[i], i_a_sq, d00m1[i], y_p_dot_t[i], d10m1[i], y_q_dot_t[i], d01m1[i]); // checked
    // type_7 = a*b + c*d + (e+f)*g + h*i
    e001[i]  = third * cuda_integral_type_7(i_r[i], c01[i], i_y[i], a01[i], s_sq[i], i_a_sq, e00m1[i], s_q_dot_t[i], e01m1[i]); // checked
    f001[i]  = third * cuda_integral_type_7(i_s[i], c01[i], i_y[i], b01[i], r_sq[i], i_a_sq, f00m1[i], r_p_dot_t[i], f01m1[i]); // checked
    // type_2 = a + b*c
    d11m1[i] = cuda_integral_type_2(a11[i], -y_q_dot_t[i], d10m1[i]); // checked
    // type_4 = a + b*c + d*e + f*g
    e11m1[i] = one_o_one_m_p_dot_t_sq * cuda_integral_type_4(a11[i], -r_p_dot_t[i], c01[i], i_p_dot_t , e001[i], -s_q_dot_t[i], e10m1[i]); // checked
    f11m1[i] = one_o_one_m_q_dot_t_sq * cuda_integral_type_4(b11[i], -s_q_dot_t[i], c01[i], i_q_dot_t , f001[i], -r_p_dot_t[i], f10m1[i]); // checked
    // type_3 = a + b*c + d*e
    d02m1[i] = cuda_integral_type_3(-d001[i], -y_q_dot_t[i], d01m1[i],  i_s[i], a01 [i]); // checked
    d20m1[i] = cuda_integral_type_3(-d001[i], -y_p_dot_t[i], d10m1[i],  i_r[i], b01 [i]); // checked
    e20m1[i] = cuda_integral_type_3(-e001[i], -i_p_dot_t     , e11m1[i],  i_r[i], c01 [i]); // checked
    f20m1[i] = cuda_integral_type_3(-f001[i], -i_q_dot_t     , f11m1[i],  i_s[i], c01 [i]); // checked
    d12m3[i] = cuda_integral_type_3(d10m1[i], -y_q_dot_t[i], d11m3[i], -i_s[i], a1m1[i]); // checked
    d21m3[i] = cuda_integral_type_3(d01m1[i], -y_p_dot_t[i], d11m3[i], -i_r[i], b1m1[i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    e21m3[i] = one_o_one_m_p_dot_t_sq * cuda_integral_type_5(e01m1[i], y_p_dot_t[i], a1m1[i], -i_p_dot_t, e10m1[i], p_dot_t_q_dot_t*i_s[i], e11m3[i], -i_r[i], c1m1[i]); // checked
    f21m3[i] = one_o_one_m_q_dot_t_sq * cuda_integral_type_5(f01m1[i], y_q_dot_t[i], b1m1[i], -i_q_dot_t, f10m1[i], p_dot_t_q_dot_t*i_r[i], f11m3[i], -i_s[i], c1m1[i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    // type_8 = a*b + c*d + e*f + g*h
    d22m3[i] = cuda_integral_type_5(-d001[i], p_dot_t_q_dot_t*y_sq[i], d11m3[i], i_r[i], b01[i], -i_r[i]*i_s[i], ra[i], i_s[i], a01[i])
             + i_y[i] * cuda_integral_type_8(-i_p_dot_t, d10m1[i], -i_q_dot_t, d01m1[i], i_p_dot_t*i_s[i], a1m1[i], i_q_dot_t*i_r[i], b1m1[i]); // checked
    // Triple integrals
    // type_3 = a + b*c + d*e
    h[0][i] = -one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_3(d00m1[i], -i_q_dot_t, e00m1[i], -i_p_dot_t, f00m1[i]); // checked
    // type_2 = a + b*c
    h[1][i] = cuda_integral_type_2(-f00m1[i], -i_p_dot_t, h[0][i]); // checked
    h[2][i] = cuda_integral_type_2(-e00m1[i], -i_q_dot_t, h[0][i]); // checked
    // type_3 = a + b*c + d*e
    h001m1[i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_3(d001[i], -i_q_dot_t, e001[i], -i_p_dot_t, f001[i]); // checked
    // type_2 = a + b*c
    h100m1[i] = cuda_integral_type_2(f001[i], -i_p_dot_t, h001m1[i]); // checked
    h010m1[i] = cuda_integral_type_2(e001[i], -i_q_dot_t, h001m1[i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    h[3][i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_5(-d11m1[i], r_p_dot_t[i], f10m1[i], -i_p_dot_t, h010m1[i], s_q_dot_t[i], e10m1[i], -i_q_dot_t, h100m1[i]); // checked
    // type_3 = a + b*c + d*e
    h[4][i] = cuda_integral_type_3(h010m1[i], -i_p_dot_t, h[3][i], -i_r[i], f10m1[i]); // checked
    h[5][i] = cuda_integral_type_3(h100m1[i], -i_q_dot_t, h[3][i], -i_s[i], e10m1[i]); // checked
    // type_4 = a + b*c + d*e + f*g
    h021m3[i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_4(-d02m1[i], -two_q_dot_t, h010m1[i], i_p_dot_t, f20m1[i], i_q_dot_t*s_sq[i], e00m1[i]); // checked
    h201m3[i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_4(-d20m1[i], -two_p_dot_t, h100m1[i], i_q_dot_t, e20m1[i], i_p_dot_t*r_sq[i], f00m1[i]);
    // type_8 = a*b + c*d + e*f + g*h
    h000m1[i] = 0.5 * cuda_integral_type_8(-i_a_sq, 0.0, i_y[i], d00m1[i], i_s[i], e00m1[i], i_r[i], f00m1[i]); // checked
    // type_4 = a + b*c + d*e + f*g
    h[6][i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_4(-d10m1[i], r_p_dot_t[i], f00m1[i], i_q_dot_t, e10m1[i], -i_p_dot_t, h000m1[i]); // checked
    h[8][i] = one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_4(-d01m1[i], s_q_dot_t[i], e00m1[i], i_p_dot_t, f10m1[i], -i_q_dot_t, h000m1[i]); // checked
    // type_2 = a + b*c
    h[7][i] = cuda_integral_type_2(-e10m1[i], -i_q_dot_t, h[6][i]); // checked
    // type_3 = a + b*c + d*e
    h[9] [i] = cuda_integral_type_3(h000m1[i], -i_p_dot_t, h[6][i], -i_r[i], f00m1[i]); // checked
    h[10][i] = cuda_integral_type_3(h000m1[i], -i_q_dot_t, h[8][i], -i_s[i], e00m1[i]); // checked
    h[11][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_3(-d00m3[i], i_p_dot_t, f00m3[i], i_q_dot_t, e00m3[i]); // checked
    // type_2 = a + b*c
    h[12][i] = cuda_integral_type_2(-third * f00m3[i], -i_p_dot_t, h[11][i]); // checked
    h[13][i] = cuda_integral_type_2(-third * e00m3[i], -i_q_dot_t, h[11][i]); // checked
    // type_4 = a + b*c + d*e + f*g
    h[14][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_4(-d10m3[i], r_p_dot_t[i], f00m3[i], i_q_dot_t, e10m3[i], -i_p_dot_t, 0.0); // checked
    h[15][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_4(-d01m3[i], s_q_dot_t[i], e00m3[i], i_p_dot_t, f10m3[i], -i_q_dot_t, 0.0); // checked
    // type_2 = a + b*c
    h[16][i] = cuda_integral_type_2(-third * f10m3[i], -i_p_dot_t, h[15][i]); // checked
    // type_3 = a + b*c + d*e
    h[17][i] = cuda_integral_type_3(third * 0.0, -third*i_r[i], f00m3[i], -i_p_dot_t, h[14][i]); // checked
    h[18][i] = cuda_integral_type_3(third * 0.0, -third*i_s[i], e00m3[i], -i_q_dot_t, h[15][i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    h[19][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_5(-d11m3[i], r_p_dot_t[i], f10m3[i], -i_p_dot_t, h[2][i], s_q_dot_t[i], e10m3[i], -i_q_dot_t, h[1][i]); // checked
    // type_3 = a + b*c + d*e
    h[20][i] = cuda_integral_type_3(third * h[2][i], -third*i_r[i], f10m3[i], -i_p_dot_t, h[19][i]); // checked
    h[21][i] = cuda_integral_type_3(third * h[1][i], -third*i_s[i], e10m3[i], -i_q_dot_t, h[19][i]); // checked
    // type_4 = a + b*c + d*e + f*g
    h[23][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_4(-d20m3[i], -two_p_dot_t, h[1][i], i_q_dot_t, e20m3[i], i_p_dot_t*r_sq[i], f00m3[i]); // checked
    h[22][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_4(-d02m3[i], -two_q_dot_t, h[2][i], i_p_dot_t, f20m3[i], i_q_dot_t*s_sq[i], e00m3[i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    h[25][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_5(h[2][i], s_q_dot_t[i], e01m3[i], -i_y[i], d01m3[i], i_p_dot_t, f11m3[i], -i_q_dot_t, h[0][i]); // checked
    h[24][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_5(h[1][i], r_p_dot_t[i], f01m3[i], -i_y[i], d10m3[i], i_q_dot_t, e11m3[i], -i_p_dot_t, h[0][i]); // checked
    // type_9 = a + b*c + d*e + f*g + h*i + j*k
    h[26][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_9(h[7][i], r_p_dot_t[i], f11m3[i], -i_y[i], d11m3[i], s_q_dot_t[i], e11m3[i], -i_p_dot_t, h[8][i], -i_q_dot_t, h[6][i]); // checked
    // type_3 = a + b*c + d*e
    h[27][i] = cuda_integral_type_3(third * h[8][i], -i_p_dot_t, h[26][i], -third*i_r[i], f11m3[i]); // checked
    h[28][i] = cuda_integral_type_3(third * h[6][i], -i_q_dot_t, h[26][i], -third*i_s[i], e11m3[i]); // checked
    // type_5 = a + b*c + d*e + f*g + h*i
    h[29][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_5(h[9] [i], -two_p_dot_t, h[6][i], i_q_dot_t, e21m3[i], i_p_dot_t*r_sq[i], f01m3[i], -i_y[i], d20m3[i]); // checked
    h[30][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_5(h[10][i], -two_q_dot_t, h[8][i], i_p_dot_t, f21m3[i], i_q_dot_t*s_sq[i], e01m3[i], -i_y[i], d02m3[i]); // checked
    // type_3 = a + b*c + d*e
    h[31][i] = cuda_integral_type_3(two_third * h[6][i], -i_p_dot_t, h[29][i], -third*r_sq[i], f01m3[i]); // checked
    h[32][i] = cuda_integral_type_3(two_third * h[8][i], -i_q_dot_t, h[30][i], -third*s_sq[i], e01m3[i]); // checked
    // type_10 = a + b*c + d*e + f*g + h*i + j*k
    h[33][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_10(h[5][i], -two_q_dot_t, h[3][i], -i_y[i], d12m3[i], r_p_dot_t[i], f21m3[i], -i_p_dot_t, h021m3[i], i_q_dot_t*s_sq[i], e11m3[i]); // checked
    h[37][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_10(h[4][i], -two_p_dot_t, h[3][i], -i_y[i], d21m3[i], s_q_dot_t[i], e21m3[i], -i_q_dot_t, h201m3[i], i_p_dot_t*r_sq[i], f11m3[i]); // checked
    // type_11 = a + b*c + d*e + f*g + h*i
    h[34][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_11(-d22m3[i], i_p_dot_t*r_sq[i], f20m3[i], i_q_dot_t*s_sq[i], e20m3[i], -i_p_dot_t*2.0, h[5][i], -i_q_dot_t*2.0, h[4][i]);
    // type_12 = a*b + c*d + e*f + g*h + i*j + k*l*m;
    h[36][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_12(2.0*one_m_q_dot_t_sq, h[3][i], -i_p_dot_t, h[4][i], p_dot_t_q_dot_t, h201m3[i], y_p_dot_t[i], d21m3[i], -p_dot_t_q_dot_t*i_s[i], e21m3[i], -one_m_q_dot_t_sq, r_sq[i], f11m3[i]); // checked
    h[35][i] = trd_one_o_one_m_p_dot_t_sq_m_q_dot_t_sq * cuda_integral_type_12(2.0*one_m_p_dot_t_sq, h[3][i], -i_q_dot_t, h[5][i], p_dot_t_q_dot_t, h021m3[i], y_q_dot_t[i], d12m3[i], -p_dot_t_q_dot_t*i_r[i], f21m3[i], -one_m_p_dot_t_sq, s_sq[i], e11m3[i]); // checked
  }
  // Evaluating the integrals.
  for (int i = 0; i < 38; i++){
    o_sch[i] = cuda_dot_product(h[i], signv, n_limits);
  }
  return o_sch;
}

__device__ float *cuda_vertex_force_linear_rectangle(float *i_sch, float i_vec_int[][3],
                                      float i_r, float i_s, float i_rs,
                                      float i_factor,
                                      float *o_force){
  // Calculate the force on a vertex.
  float f[11];
  f [0] = i_sch [3] - i_s*i_sch [6] - i_r*i_sch [8] + i_rs*i_sch [0];
  f [1] = i_sch [5] - i_s*i_sch [7] - i_r*i_sch[10] + i_rs*i_sch [2];
  f [2] = i_sch [4] - i_s*i_sch [9] - i_r*i_sch [7] + i_rs*i_sch [1];
  f [3] = i_sch[19] - i_s*i_sch[14] - i_r*i_sch[15] + i_rs*i_sch[11];
  f [4] = i_sch[21] - i_s*i_sch[16] - i_r*i_sch[18] + i_rs*i_sch[13];
  f [5] = i_sch[20] - i_s*i_sch[17] - i_r*i_sch[16] + i_rs*i_sch[12];
  f [6] = i_sch[35] - i_s*i_sch[28] - i_r*i_sch[32] + i_rs*i_sch[22];
  f [7] = i_sch[36] - i_s*i_sch[31] - i_r*i_sch[27] + i_rs*i_sch[23];
  f [8] = i_sch[33] - i_s*i_sch[26] - i_r*i_sch[30] + i_rs*i_sch[25];
  f [9] = i_sch[37] - i_s*i_sch[29] - i_r*i_sch[26] + i_rs*i_sch[24];
  f[10] = i_sch[34] - i_s*i_sch[27] - i_r*i_sch[28] + i_rs*i_sch[19];
  for (int i = 0; i < 11; i++){
    for (int j = 0; j < 3; j++){
      o_force[j] += i_vec_int[i][j] * f[i];
    }
  }
  o_force[0] = o_force[0]*i_factor;
  o_force[1] = o_force[1]*i_factor;
  o_force[2] = o_force[2]*i_factor;
  return o_force;
}

__device__ void cuda_compute_forces_linear_rectangle(float *i_sch, float i_vec_int[][3],
                                     float *i_rp, float *i_sp, float i_factor,
                                     float o_nodal_force[][3], float *o_total_force){
  // Calculating nodal forces
  // x3.
  cuda_vertex_force_linear_rectangle(i_sch, i_vec_int, i_rp[1], i_sp[1], i_rp[1]*i_sp[1],  i_factor, o_nodal_force[0]);
  // x4.
  cuda_vertex_force_linear_rectangle(i_sch, i_vec_int, i_rp[0], i_sp[1], i_rp[0]*i_sp[1], -i_factor, o_nodal_force[1]);
  // x5.
  cuda_vertex_force_linear_rectangle(i_sch, i_vec_int, i_rp[1], i_sp[0], i_rp[1]*i_sp[0], -i_factor, o_nodal_force[2]);
  // x5
  cuda_vertex_force_linear_rectangle(i_sch, i_vec_int, i_rp[0], i_sp[0], i_rp[0]*i_sp[0],  i_factor, o_nodal_force[3]);
  for (int i = 0; i < n_nodes; i++){
    o_total_force[0] += o_nodal_force[i][0];
    o_total_force[1] += o_nodal_force[i][1];
    o_total_force[2] += o_nodal_force[i][2];
  }
}

// Device constants.
__constant__ float d_mu, d_nu, d_a, d_a_sq, d_one_m_nu, d_factor;

__global__ void se_cuda_nodal_surface_force_linear_rectangle(float *g_dln_arr, float *g_se_arr, float *g_b_arr, float *g_fx_arr, float *g_ftot_arr, int n_se, int n_dln){
  float x1[3], x2[3], x3[3], x4[3], x5[3], x6[3], b[3];
  // Characteristic vectors.
  float t[3], p[3], q[3], n[3];
  // Basis vectors (unitary).
  float p_x_t[3], q_x_t[3];
  // Limits of the distance vector from the plane to dislocation line segment.
  // r_lim[0][] = vector from x1 to x3, r_lim[1][] = vector from x2 to x6.
  float r_lim[2][3];
  //r, s limits
  float rp[2], sp[2];
  //  y, r, s coordinates
  float y[n_limits], r[n_limits], s[n_limits];
  // Vectors for the integrals.
  float vec_int[11][3];
  // Scalar value of the integrals.
  float sch[38];
  //  Auxiliary constants to reduce computational cost.
  //  p dot t, q dot t
  float p_dot_t, q_dot_t;
  // Auxiliary variables for nodal force calculation.
  float l_factor;
  float p_norm, q_norm;
  float nodal_force[n_nodes][3], total_force[3];
  int idx  = threadIdx.x + blockIdx.x * blockDim.x;
  int idxd = 0;
  int idxb = 0;
  // Because each index corresponds to a surface element, we have to make sure that even if we have more threads than surface elements we won't go out of bounds in our arrays.
  if (idx < n_se){
    se_device_thread_map(g_se_arr, x3, x4, x5, x6, n_se, idx);
    for (int j = 0; j < n_dln; j++){
      cuda_init_force(nodal_force, total_force);
      x1[0] = g_dln_arr[idxd  ];
      x1[1] = g_dln_arr[idxd+1];
      x1[2] = g_dln_arr[idxd+2];
      x2[0] = g_dln_arr[idxd+3];
      x2[1] = g_dln_arr[idxd+4];
      x2[2] = g_dln_arr[idxd+5];
      b[0] = g_b_arr[idxb];
      b[1] = g_b_arr[idxb+1];
      b[2] = g_b_arr[idxb+2];
      idxd += 6;
      idxb += 3;
      // Vectors
      cuda_init_vector(x1, x2, 3, t);
      cuda_init_vector2(x3, x4, 3, p, &p_norm);
      cuda_init_vector2(x3, x5, 3, q, &q_norm);
      cuda_cross_product(p, q, n);
      // Only register results where t dot n is different to zero. The special case will be treated outside the parallel code outside in the serial portion of the code. This is due to warp branching, which means all code branches are evaluated but changes in memory are modified only for those instructions whose branch conditions have been met. Having the special case here means having to execute every instruction of every branch, severly hindering performance. The branching is up here in case this changes in the future.
      if (cuda_dot_product(t, n, 3) != 0.0){
        cuda_normalise_vector(n, 3, n);
        // Local factor.
        l_factor = d_factor/p_norm/q_norm;
        // Dot products.
        p_dot_t = cuda_dot_product(p, t, 3);
        q_dot_t = cuda_dot_product(q, t, 3);
        //*******************WARNING*******************//
        // This formulation assumes x3-x6 is diagonal! //
        //*******************WARNING*******************//
        cuda_cross_product(p, t, p_x_t);
        cuda_cross_product(q, t, q_x_t);
        // Vectors between x3 and x1, and x6 and x2.
        for (int i = 0; i < 3; i++){
          r_lim[0][i] = x3[i] - x1[i];
          r_lim[1][i] = x6[i] - x2[i];
        }
        // Integral bounds for y, r, s.
        for (int i = 0; i < 2; i++){
          rp[i] = cuda_init_point(r_lim[i], q_x_t, p, q_x_t, 3);
          sp[i] = cuda_init_point(r_lim[i], p_x_t, q, p_x_t, 3);
        }
        // Assign coordinates for the evaluation of the integrals.
        y[0] = y[2] = y[4] = y[6] = cuda_init_point(r_lim[1], n, t, n, 3);
        y[1] = y[3] = y[5] = y[7] = cuda_init_point(r_lim[0], n, t, n, 3);
        r[0] = r[1] = r[2] = r[3] = rp[1];
        r[4] = r[5] = r[6] = r[7] = rp[0];
        s[0] = s[1] = s[4] = s[5] = sp[1];
        s[2] = s[3] = s[6] = s[7] = sp[0];
        // Calculate vectors for integrals.
        cuda_integral_vector(p, q, b, t, n, d_one_m_nu, d_a_sq, vec_int);
        // Calculate integrals.
        cuda_integrals_linear_rectangle(r, s, y, p, q, t, p_dot_t, q_dot_t, d_a_sq, sch);
        // Calculate nodal forces.
        cuda_compute_forces_linear_rectangle(sch, vec_int, rp, sp, l_factor, nodal_force,
           total_force);
        add_force_thread_device(nodal_force, total_force, g_fx_arr, g_ftot_arr, n_se, n_nodes, idx);
        #ifdef debug
        /*
          printf("SE = %d\n", (threadIdx.x + blockIdx.x * blockDim.x)%n_se+1);
          printf("t = %f %f %f\n", t[0], t[1], t[2]);
          printf("x1 = %f %f %f\n", x1[0], x1[1], x1[2]);
          printf("p = %f %f %f\n", p[0], p[1], p[2]);
          printf("q = %f %f %f\n", q[0], q[1], q[2]);
          printf("n = %f %f %f\n", n[0], n[1], n[2]);
          printf("b = %f %f %f\n", b[0], b[1], b[2]);
          for (int i = 0; i < 11; i++){printf("vec_int[%d] = %1.8f %1.8f %1.8f\n", i, vec_int[i][0],vec_int[i][1],vec_int[i][2]);}
          for (int i = 0; i < n_nodes; i++){printf("nodal_force[%d] = %1.8f %1.8f %1.8f\n", i, nodal_force[i][0], nodal_force[i][1], nodal_force[i][2]);}
        */
        /*
        for (int i = 0; i < 4; i++){
          printf("SE = %d\nnodal_force[%d] = %2.14f %2.14f %2.14f\n",(threadIdx.x + blockIdx.x * blockDim.x)%n_se+1, i, nodal_force[i][0], nodal_force[i][1], nodal_force[i][2]);
        }
        */
        printf("dln = %d, SE = %d\ntotal_force = %2.14f %2.14f %2.14f\n", j+1, (threadIdx.x + blockIdx.x * blockDim.x)%n_se+1, total_force[0], total_force[1], total_force[2]);
        #endif
      }
    }
  }
}

__global__ void dln_cuda_nodal_surface_force_linear_rectangle(float *g_dln_arr, float *g_se_arr, float *g_b_arr, float *g_fx_arr, float *g_ftot_arr, int n_se, int n_dln){
  float x1[3], x2[3], x3[3], x4[3], x5[3], x6[3], b[3];
  // Characteristic vectors.
  float t[3], p[3], q[3], n[3];
  // Basis vectors (unitary).
  float p_x_t[3], q_x_t[3];
  // Limits of the distance vector from the plane to dislocation line segment.
  // r_lim[0][] = vector from x1 to x3, r_lim[1][] = vector from x2 to x6.
  float r_lim[2][3];
  //r, s limits
  float rp[2], sp[2];
  //  y, r, s coordinates
  float y[n_limits], r[n_limits], s[n_limits];
  // Vectors for the integrals.
  float vec_int[11][3];
  // Scalar value of the integrals.
  float sch[38];
  //  Auxiliary constants to reduce computational cost.
  //  p dot t, q dot t
  float p_dot_t, q_dot_t;
  // Auxiliary variables for nodal force calculation.
  float l_factor;
  float p_norm, q_norm;
  float nodal_force[n_nodes][3], total_force[3];
  int idx   = threadIdx.x + blockIdx.x * blockDim.x;
  int idxse = 0;
  // Because each index corresponds to a dislocation line segment, we have to make sure that even if we have more threads than surface elements we won't go out of bounds in our arrays.
  if (idx < n_dln){
    dln_device_thread_map(g_dln_arr, x1, x2, n_dln, idx);
    b_device_thread_map(g_b_arr, b, n_dln, idx);
    for (int j = 0; j < n_se; j++){
      cuda_init_force(nodal_force, total_force);
      x3[0] = g_se_arr[idxse  ];
      x3[1] = g_se_arr[idxse+1];
      x3[2] = g_se_arr[idxse+2];
      x4[0] = g_se_arr[idxse+3];
      x4[1] = g_se_arr[idxse+4];
      x4[2] = g_se_arr[idxse+5];
      x5[0] = g_se_arr[idxse+6];
      x5[1] = g_se_arr[idxse+7];
      x5[2] = g_se_arr[idxse+8];
      x6[0] = g_se_arr[idxse+9];
      x6[1] = g_se_arr[idxse+10];
      x6[2] = g_se_arr[idxse+11];
      idxse += 12;
      // Vectors
      cuda_init_vector(x1, x2, 3, t);
      cuda_init_vector2(x3, x4, 3, p, &p_norm);
      cuda_init_vector2(x3, x5, 3, q, &q_norm);
      cuda_cross_product(p, q, n);
      // Only register results where t dot n is different to zero. The special case will be treated outside the parallel code outside in the serial portion of the code. This is due to warp branching, which means all code branches are evaluated but changes in memory are modified only for those instructions whose branch conditions have been met. Having the special case here means having to execute every instruction of every branch, severly hindering performance. The branching is up here in case this changes in the future.
      if (cuda_dot_product(t, n, 3) != 0.0){
        cuda_normalise_vector(n, 3, n);
        // Local factor.
        l_factor = d_factor/p_norm/q_norm;
        // Dot products.
        p_dot_t = cuda_dot_product(p, t, 3);
        q_dot_t = cuda_dot_product(q, t, 3);
        //*******************WARNING*******************//
        // This formulation assumes x3-x6 is diagonal! //
        //*******************WARNING*******************//
        cuda_cross_product(p, t, p_x_t);
        cuda_cross_product(q, t, q_x_t);
        // Vectors between x3 and x1, and x6 and x2.
        for (int i = 0; i < 3; i++){
          r_lim[0][i] = x3[i] - x1[i];
          r_lim[1][i] = x6[i] - x2[i];
        }
        // Integral bounds for y, r, s.
        for (int i = 0; i < 2; i++){
          rp[i] = cuda_init_point(r_lim[i], q_x_t, p, q_x_t, 3);
          sp[i] = cuda_init_point(r_lim[i], p_x_t, q, p_x_t, 3);
        }
        // Assign coordinates for the evaluation of the integrals.
        y[0] = y[2] = y[4] = y[6] = cuda_init_point(r_lim[1], n, t, n, 3);
        y[1] = y[3] = y[5] = y[7] = cuda_init_point(r_lim[0], n, t, n, 3);
        r[0] = r[1] = r[2] = r[3] = rp[1];
        r[4] = r[5] = r[6] = r[7] = rp[0];
        s[0] = s[1] = s[4] = s[5] = sp[1];
        s[2] = s[3] = s[6] = s[7] = sp[0];
        // Calculate vectors for integrals.
        cuda_integral_vector(p, q, b, t, n, d_one_m_nu, d_a_sq, vec_int);
        // Calculate integrals.
        cuda_integrals_linear_rectangle(r, s, y, p, q, t, p_dot_t, q_dot_t, d_a_sq, sch);
        // Calculate nodal forces.
        cuda_compute_forces_linear_rectangle(sch, vec_int, rp, sp, l_factor, nodal_force,
           total_force);
        add_force_thread_device(nodal_force, total_force, g_fx_arr, g_ftot_arr, n_se, n_nodes, j);
        #ifdef debug
        /*
          printf("SE = %d\n", (threadIdx.x + blockIdx.x * blockDim.x)%n_se+1);
          printf("t = %f %f %f\n", t[0], t[1], t[2]);
          printf("x1 = %f %f %f\n", x1[0], x1[1], x1[2]);
          printf("p = %f %f %f\n", p[0], p[1], p[2]);
          printf("q = %f %f %f\n", q[0], q[1], q[2]);
          printf("n = %f %f %f\n", n[0], n[1], n[2]);
          printf("b = %f %f %f\n", b[0], b[1], b[2]);
          for (int i = 0; i < 11; i++){printf("vec_int[%d] = %1.8f %1.8f %1.8f\n", i, vec_int[i][0],vec_int[i][1],vec_int[i][2]);}
          for (int i = 0; i < n_nodes; i++){printf("nodal_force[%d] = %1.8f %1.8f %1.8f\n", i, nodal_force[i][0], nodal_force[i][1], nodal_force[i][2]);}
        */
        /*
        for (int i = 0; i < 4; i++){
          printf("SE = %d\nnodal_force[%d] = %2.14f %2.14f %2.14f\n",(threadIdx.x + blockIdx.x * blockDim.x)%n_se+1, i, nodal_force[i][0], nodal_force[i][1], nodal_force[i][2]);
        }
        */
        printf("dln = %d, SE = %d\ntotal_force = %2.14f %2.14f %2.14f\n", (threadIdx.x + blockIdx.x * blockDim.x)%n_dln+1, j+1, total_force[0], total_force[1], total_force[2]);
        #endif
      }
    }
  }
}


/* Work out best parallelisation.
// Maximum number of threads per block.
//cudaDeviceProp deviceProp;
//cudaGetDeviceProperties(&deviceProp, 1);
// Number of blocks launched.
// deviceProp.maxGridSize[0], deviceProp.maxGridSize[1], deviceProp.maxGridSize[2]
// Size of each block launched
// deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1], deviceProp.maxThreadsDim[2]
*/
/*
void main_se_cuda_nodal_surface_force_linear_rectangle(int n_se, int n_dln){
  // Node arrays from MATLAB. To be mapped into x_se_arr and then passed to d_x_se_arr.
  float *dln_node_arr[2], *se_node_arr[n_nodes];
  // Variables for the special case where the line segment and surface element are parallel.
  float x1[3], x2[3], x3[3], x4[3], x5[3], x6[3], b[3], p[3], q[3], t[3], n[3], p_norm, q_norm;
  float *fx[n_nodes];
  float ftot[3];
  // Burger's vectors array from MATLAB. To be passed straight into d_b_arr.
  float *b_arr[1];
  // Material properties from MATLAB to be placed into shared memory in device.
  float mu, nu, a, a_sq, one_m_nu, factor;
  // Nodal force array (3 coordinates per node per SE, 3*n_nodes*n_se). To be inversely mapped to *fx_arr[n_nodes].
  float *x_fx_arr;
  // Nodal force array to be sent back to MATLAB.
  float *fx_arr[n_nodes];
  // Total force array on SE (3 coordinates per SE, 3*n_se) to be sent back to MATLAB.
  float *ftot_arr, *x_ftot_arr;
  // Maps of SE and DLN node arrays.
  float *x_se_arr, *x_dln_arr, *x_b_arr;
  // Device arrays.
  float *d_x_b_arr, *d_x_se_arr, *d_x_dln_arr, *d_fx_arr, *d_ftot_arr;
  int idx1, idx2;

  // Memory allocation
  b_arr[0] = (float *) malloc(3 * n_dln * sizeof(float));
  dln_node_arr[0] = (float *) malloc(3 * n_dln * sizeof(float));
  dln_node_arr[1] = (float *) malloc(3 * n_dln * sizeof(float));
  for (int i = 0; i < n_nodes; i++){
    se_node_arr[i] = (float *) malloc(3 * n_se  * sizeof(float));
    fx[i] = (float *) malloc(3 * sizeof(float));
  }
  // Read input.
  FILE * ptr_file;
  ptr_file = fopen("cuda_input.txt", "r");
  if (ptr_file == NULL){
    printf("File does not exist.\n");
  }
  // Skip two lines.
  fscanf(ptr_file, "%*[^\n]\n");
  fscanf(ptr_file, "%*[^\n]\n");
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
    fscanf(ptr_file, "%lf %lf %lf", &b_arr[0][i], &b_arr[0][i+1], &b_arr[0][i+2] );
  }
  fscanf(ptr_file, "%lf", &mu );
  fscanf(ptr_file, "%lf", &nu );
  fscanf(ptr_file, "%lf", &a );
  fclose(ptr_file);
  // Auxiliary constants.
  a_sq     = a*a;
  one_m_nu = 1.-nu;
  factor   = 0.25*mu/pi/one_m_nu;
  // Forward data map.
  x_b_arr   = b_host_device_map(b_arr[0], n_dln);
  x_dln_arr = dln_host_device_map(dln_node_arr[0], dln_node_arr[1], n_dln);
  x_se_arr  = element_host_device_map(se_node_arr, n_se, n_nodes);
  #ifdef debug
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
  #endif
  // Allocate device memory.
  checkCudaErrors( cudaMalloc( (void **) &d_x_dln_arr, 3 * n_dln * 2       * sizeof(float) ) );
  checkCudaErrors( cudaMalloc( (void **) &d_x_b_arr  , 3 * n_dln           * sizeof(float) ) );
  checkCudaErrors( cudaMalloc( (void **) &d_x_se_arr , 3 * n_se  * n_nodes * sizeof(float) ) );
  checkCudaErrors( cudaMalloc( (void **) &d_fx_arr   , 3 * n_se  * n_nodes * sizeof(float) ) );
  checkCudaErrors( cudaMalloc( (void **) &d_ftot_arr , 3 * n_se            * sizeof(float) ) );
  // Copy host to device.
  // Only pass node coordinates to device.
  checkCudaErrors( cudaMemcpy(d_x_se_arr , x_se_arr , 3*n_se *n_nodes*sizeof(float), cudaMemcpyHostToDevice) );
  free(x_se_arr);
  checkCudaErrors( cudaMemcpy(d_x_dln_arr, x_dln_arr, 3*n_dln*2      *sizeof(float), cudaMemcpyHostToDevice) );
  free(x_dln_arr);
  checkCudaErrors( cudaMemcpy(d_x_b_arr  , x_b_arr, 3*n_dln        *sizeof(float), cudaMemcpyHostToDevice) );
  free(x_b_arr);
  // Initialising force arrays to zero.
  checkCudaErrors( cudaMemset(d_fx_arr  , 0.0, 3*n_se *n_nodes*sizeof(float)) );
  checkCudaErrors( cudaMemset(d_ftot_arr, 0.0, 3*n_se         *sizeof(float)) );
  checkCudaErrors( cudaMemcpyToSymbol(d_mu      , &mu      , sizeof(mu)) );
  checkCudaErrors( cudaMemcpyToSymbol(d_nu      , &nu      , sizeof(nu)) );
  checkCudaErrors( cudaMemcpyToSymbol(d_a       , &a       , sizeof(a)) );
  checkCudaErrors( cudaMemcpyToSymbol(d_a_sq    , &a_sq    , sizeof(a_sq)) );
  checkCudaErrors( cudaMemcpyToSymbol(d_one_m_nu, &one_m_nu, sizeof(one_m_nu)) );
  checkCudaErrors( cudaMemcpyToSymbol(d_factor  , &factor  , sizeof(factor)) );
  // CUDA
  //nodal_surface_force_linear_rectangle<<<1,n_se>>>(d_x_se_arr);
  //nodal_surface_force_linear_rectangle<<<n_se,1>>>(d_x_se_arr);
  #ifdef debug
    cudaEventRecord(start);
  #endif
  se_cuda_nodal_surface_force_linear_rectangle<<<6,6>>>(d_x_dln_arr, d_x_se_arr, d_x_b_arr, d_fx_arr, d_ftot_arr, n_se, n_dln);
  #ifdef debug
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    printf("Time spent in SE parallelisation: %f (ms) \n", milliseconds);
  #endif
  // Copy device to host
  // Only pass forces to host.
  // x_fx_arr treated the same way as x_se_arr.
  x_fx_arr   = (float *) malloc(3 * n_se * n_nodes * sizeof(float));
  x_ftot_arr = (float *) malloc(3 * n_se           * sizeof(float));
  checkCudaErrors( cudaMemcpy(x_fx_arr  , d_fx_arr  , 3 * n_se * n_nodes * sizeof(float), cudaMemcpyDeviceToHost) );
  checkCudaErrors( cudaMemcpy(x_ftot_arr, d_ftot_arr, 3 * n_se           * sizeof(float), cudaMemcpyDeviceToHost) );
  ftot_arr = (float *) malloc(3 * n_se * sizeof(float));
  ftot_device_host_map(x_ftot_arr, ftot_arr, n_se);
  free(x_ftot_arr);
  for (int i = 0; i < n_nodes; i++){
    fx_arr[i] = (float *) malloc(3 * n_se * sizeof(float));
  }
  fx_device_host_map(x_fx_arr, fx_arr, n_se, n_nodes);
  free(x_fx_arr);
  // Hunt for the special case of parallel line segments to surface elements.
  cudaDeviceSynchronize();
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
        nodal_surface_force_linear_rectangle_special(x1, x2, x3, x4, x5, x6, b, t, p, q, n, p_norm, q_norm, mu, nu, a, fx, ftot);
        // Add the force contributions for segment j to the surface element i.
        for (int k = 0; k < 3; k++){
          //printf("dln = %d, se = %d, ftot[%d] = %f\n", j, i, k, ftot[k]);
          //printf("x1 = %f %f %f\n", x1[0], x1[1], x1[2]);
          //printf("x2 = %f %f %f\n", x2[0], x2[1], x2[2]);
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
  #ifdef debug
    //for (int i = 0; i < 3*n_se*n_nodes; i++){
    //  printf("x_fx_arr[%d] = %2.14f\n", i, x_fx_arr[i]);
    //}
    for (int i = 0; i < 3*n_se; i++){
      printf("x_ftot_arr[%d] = %2.14f\n", i, ftot_arr[i]);
    }
    //for(int i = 0; i < n_nodes; i++){
    //  for(int j = 0; j < 3*n_se; j+=3){
    //    printf("fx_arr[%d] = %2.18f %2.18f %2.18f\n", i, fx_arr[i][j], fx_arr[i][j+1], fx_arr[i][j+2]);
    //  }
    //}
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
  #endif
  // CUDA exit.
  cudaDeviceReset();
  for(int i = 0; i < n_nodes; i++){free(fx[i]);}
  // Don't free these when using MATLAB. It silently crashes the program. Free them when using straight C.
  //free(b_arr[0]); free(dln_node_arr[0]); free(dln_node_arr[1]); //free(x3_arr); free(x4_arr); free(x5_arr); free(x6_arr);
  //for(int i=0; i < n_nodes; i++){free(se_node_arr[i]);}
  //free(ftot_arr);
  //for(int i=0; i < n_nodes; i++){free(fx_arr[i]);}
}
*/

// main_dln_cuda_nodal_surface_force_linear_rectangle(int n_se, int n_dln)
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
  // Node arrays from MATLAB. To be mapped into x_se_arr and then passed to d_x_se_arr.
  float *dln_node_arr[2], *se_node_arr[n_nodes];
  // Variables for the special case where the line segment and surface element are parallel.
  float x1[3], x2[3], x3[3], x4[3], x5[3], x6[3], b[3], p[3], q[3], t[3], n[3], p_norm, q_norm;
  float *fx[n_nodes];
  float ftot[3];
  // Burger's vectors array from MATLAB. To be passed straight into d_b_arr.
  float *b_arr[1];
  // Material properties from MATLAB to be placed into shared memory in device.
  float mu, nu, a, a_sq, one_m_nu, factor;
  // Nodal force array (3 coordinates per node per SE, 3*n_nodes*n_se). To be inversely mapped to *fx_arr[n_nodes].
  float *x_fx_arr;
  // Nodal force array to be sent back to MATLAB.
  float *fx_arr[n_nodes];
  // Total force array on SE (3 coordinates per SE, 3*n_se) to be sent back to MATLAB.
  float *ftot_arr, *x_ftot_arr;
  // Maps of SE and DLN node arrays.
  float *x_se_arr, *x_dln_arr, *x_b_arr;
  // Device arrays.
  float *d_x_b_arr, *d_x_se_arr, *d_x_dln_arr, *d_fx_arr, *d_ftot_arr;
  int threads_per_block, blocks_per_grid, n_se, n_dln;
  int idx1, idx2;
  //int debug = 1;
  //while(debug == 1){}
  // Stagger cuda function calls to take advantage of asynchronous calls.
  // If memory becomes an issue, make copying x_dln_arr, x_se_arr and x_b_arr to the device a synchronous operation and free the pointers straight after.
  n_se = (int) mxGetScalar(prhs[10]);
  // Allocate and set forces to 0 in device.
  checkCudaErrors( cudaMalloc( (void **) &d_fx_arr  , 3 * n_se  * n_nodes * sizeof(float) ) );
  checkCudaErrors( cudaMalloc( (void **) &d_ftot_arr, 3 * n_se            * sizeof(float) ) );
  checkCudaErrors( cudaMemsetAsync(d_fx_arr  , 0.0, 3 * n_se * n_nodes * sizeof(float)) );
  checkCudaErrors( cudaMemsetAsync(d_ftot_arr, 0.0, 3 * n_se           * sizeof(float)) );
  // Execute host code while device sets force arrays to zero.
  n_dln = (int) mxGetScalar(prhs[11]);
  dln_node_arr[0] = (float *) mxGetPr(prhs[0]);
  dln_node_arr[1] = (float *) mxGetPr(prhs[1]);
  // Map dislocation node arrays to 1D array for parallelisation.
  x_dln_arr       = element_host_device_map(dln_node_arr, n_dln, 2);
  // Allocate and copy values of dislocation nodes to device.
  checkCudaErrors( cudaMalloc     ( (void **) &d_x_dln_arr, 3 * n_dln * 2 * sizeof(float) ) );
  checkCudaErrors( cudaMemcpyAsync(d_x_dln_arr,  x_dln_arr, 3 * n_dln * 2 * sizeof(float), cudaMemcpyHostToDevice) );
  // Execute host code while copying values from host to device.
  se_node_arr[0] = (float *) mxGetPr(prhs[2]);
  se_node_arr[1] = (float *) mxGetPr(prhs[3]);
  se_node_arr[2] = (float *) mxGetPr(prhs[4]);
  se_node_arr[3] = (float *) mxGetPr(prhs[5]);
  // Map surface element node arrays to 1D array for parallelisation.
  x_se_arr       = se_host_device_map(se_node_arr[0], se_node_arr[1], se_node_arr[2], se_node_arr[3], n_se);
  // Allocate and copy values of surface element nodes to device.
  checkCudaErrors( cudaMalloc     ( (void **) &d_x_se_arr, 3 * n_se * n_nodes * sizeof(float) ) );
  checkCudaErrors( cudaMemcpyAsync(d_x_se_arr,   x_se_arr, 3 * n_se * n_nodes * sizeof(float), cudaMemcpyHostToDevice) );
  // Execute host code while copying values from host to device.
  b_arr[0] = (float *) mxGetPr(prhs[6]);
  // Map Burger's vector array to 1D array for parallelisation.
  x_b_arr  = element_host_device_map(b_arr, n_dln, 1);
  // Allocate and copy values of Burger's vectors to device.
  checkCudaErrors( cudaMalloc     ( (void **) &d_x_b_arr, 3 * n_dln * sizeof(float) ) );
  checkCudaErrors( cudaMemcpyAsync(d_x_b_arr,    x_b_arr, 3 * n_dln * sizeof(float), cudaMemcpyHostToDevice) );
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
  // Link force arrays to MATLAB.
  plhs[0] = mxCreateDoubleMatrix(3 * n_se, 1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(3 * n_se, 1, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(3 * n_se, 1, mxREAL);
  plhs[3] = mxCreateDoubleMatrix(3 * n_se, 1, mxREAL);
  plhs[4] = mxCreateDoubleMatrix(3 * n_se, 1, mxREAL);
  fx_arr[0] = (float *) mxGetPr(plhs[0]);
  fx_arr[1] = (float *) mxGetPr(plhs[1]);
  fx_arr[2] = (float *) mxGetPr(plhs[2]);
  fx_arr[3] = (float *) mxGetPr(plhs[3]);
  ftot_arr  = (float *) mxGetPr(plhs[4]);
  threads_per_block = (int) mxGetScalar(prhs[12]);
  blocks_per_grid   = (n_dln + threads_per_block - 1) / threads_per_block;
  // CUDA
  dln_cuda_nodal_surface_force_linear_rectangle<<<blocks_per_grid, threads_per_block>>>(d_x_dln_arr, d_x_se_arr, d_x_b_arr, d_fx_arr, d_ftot_arr, n_se, n_dln);
  // Host code is executed asynchronously from the kernel execution.
  // Free all 1D arrays used to copy data to device.
  free(x_se_arr); free(x_dln_arr); free(x_b_arr);
  x_fx_arr   = (float *) malloc(3 * n_se * n_nodes * sizeof(float));
  x_ftot_arr = (float *) malloc(3 * n_se           * sizeof(float));
  // Synchronously copy forces from device to host.
  checkCudaErrors( cudaMemcpy(x_fx_arr, d_fx_arr, 3 * n_se * n_nodes * sizeof(float), cudaMemcpyDeviceToHost) );
  // Map 1D device array to 2D array for MATLAB.
  fx_device_host_map(x_fx_arr, fx_arr, n_se, n_nodes);
  free(x_fx_arr);
  checkCudaErrors( cudaMemcpy(x_ftot_arr, d_ftot_arr, 3 * n_se * sizeof(float), cudaMemcpyDeviceToHost) );
  // Map 1D device array to 1D array for MATLAB.
  ftot_device_host_map(x_ftot_arr, ftot_arr, n_se);
  free(x_ftot_arr);
  // CUDA exit.
  cudaDeviceReset();
}
