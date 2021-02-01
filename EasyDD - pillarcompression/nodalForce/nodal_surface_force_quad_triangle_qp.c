//
//  SegQuadTrianForces.c
//
//
//  Created by Sylvain Queyreau on 12/12/14.
//  Refactored and modified by Daniel Celis Garza 20/06/2017
//
//
//  contains: SegQuadTrianForces (arbitrary dislocation and surface orientations)
//            SpecialSegQuadTrianForces (parallel case)
//            Main (read input data from a file)
//
//
//  These routines evaluate the nodal forces associated with the traction field induced
//  by non-singular straight dislocation onto a linear triangular element.
//  The non-singular stress field expression for isotropic elasticity is provided
//  in [Cai et al. Journal of the Mechanics and Physics of Solids 54 (2006) 561–587].
//  The surface element is a three-node linear triangular element of arbitrary
//  shape.
//
//  The nodal force evaluation is associated to a triple integration along the dislocation
//  and along the surface element. This can be achieved by a sequence of three integrations
//  by part that present recurrence relationships.
//
//  See Latex file for more details.
//
#include <quadmath.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <tgmath.h>
//#include <mex.h>
/*-------------------------------------------------------------------------
 *
 *      Function:     SegQuadTrianForces
 *      Description:  Function for calculating nodal forces induced by a
 *                    non-singular segment on a quadratic triangular surface
 *                    element. Parallel and non-parallel case.
 *
 *      Arguments:
 *          p1*,p2*      endpoints for dislocation segment beginning
 *                         at point p1 and ending at point p2. the line
 *                         vector t = p1p2
 *          p3*,p4*,p5*  endpoints for an arbitrary shaped triangular element,
 *                         3 additional mid points p6, p7, p8 are added at the
 *                         centre of each side of the triangle.
 *                         the surface normal n is defined by the cross product
 *                         n = p3p4 x p3p5.
 *          bx,by,bz     Burgers vector for segment p1p2
 *          a            core parameter
 *          mu           shear modulus
 *          nu           poisson ratio
 *          fp3* - fp8*  pointers to locations in which to return forces
 *                         on nodes p3 thru p8 respectively. Units are
 *                         defined by the units used for xi, b and mu.
 *
 *-----------------------------------------------------------------------*/
const int n_nodes = 6;
const int n_limits = 8;
const int n_limits_o_2 = 4;
 /*
  Auxiliary functions.
 */
 __float128 dot_product(__float128 *i_vec1, __float128 *i_vec2, int i_vec_size){
   // Returns the dot product of i_vec1, i_vec2.
   __float128 result = 0.0q;
   for (int i = 0; i < i_vec_size; i++){
     result += i_vec1[i]*i_vec2[i];
   }
   return result;
 }

 __float128 *cross_product(__float128 *i_vec1, __float128 *i_vec2,
                       __float128 *o_vec){
   // Returns the cross product of i_vec1 x i_vec2.
   o_vec[0] = i_vec1[1]*i_vec2[2] - i_vec1[2]*i_vec2[1];
   o_vec[1] = i_vec1[2]*i_vec2[0] - i_vec1[0]*i_vec2[2];
   o_vec[2] = i_vec1[0]*i_vec2[1] - i_vec1[1]*i_vec2[0];
   return o_vec;
 }

 __float128 *cross_product2(__float128 *i_vec1, __float128 *i_vec2){
   __float128 *o_vec = malloc(3*sizeof(__float128));
   // Returns the cross product of i_vec1 x i_vec2.
   o_vec[0] = i_vec1[1]*i_vec2[2] - i_vec1[2]*i_vec2[1];
   o_vec[1] = i_vec1[2]*i_vec2[0] - i_vec1[0]*i_vec2[2];
   o_vec[2] = i_vec1[0]*i_vec2[1] - i_vec1[1]*i_vec2[0];
   return o_vec;
 }

 __float128 *normalise_vector(__float128 *i_vec, int i_vec_size,
                          __float128 *o_vec){
   // Returns a normalised i_vec to o_vec.
   __float128 mag_vec = 0.0q;
   mag_vec = dot_product(i_vec, i_vec, i_vec_size);
   // Check magnitude is not zero.
   if (mag_vec == 0.0q){
     fprintf(stderr, "ERROR: nodal_surface_force_quad_triangle: normalise_vector: mag_vec = 0: A vector cannot have magnitude 0\n");
     exit(EXIT_FAILURE);
   }
   mag_vec = sqrtq(mag_vec);
   for(int i = 0; i < i_vec_size; i++){
     o_vec[i] = i_vec[i]/mag_vec;
   }
   return o_vec;
 }

 void normalise_vector2(__float128 *i_vec, int i_vec_size,
                        __float128 *o_vec, __float128 *o_mag_vec){
   // Returns a normalised i_vec to o_vec, and the magnitude of the vector in o_mag_vec.
   // Has to be a void function in order to 'return' two values, the magnitude should be passed as a reference eg:
   /*
     int    size = 3;
     __float128 i_vec[size], o_vec[size], magnitude;
     normalise_vector2(i_vec, size, o_vec, &magnitude);
   */
   *o_mag_vec = dot_product(i_vec, i_vec, i_vec_size);
   // Check magnitude is not zero.
   if (*o_mag_vec == 0.0q){
     fprintf(stderr, "ERROR: nodal_surface_force_quad_triangle: normalise_vector2: o_mag_vec = 0: A vector cannot have magnitude 0\n");
     exit(EXIT_FAILURE);
   }
   *o_mag_vec = sqrtq(*o_mag_vec);
   for(int i = 0; i < i_vec_size; i++){
     o_vec[i] = i_vec[i]/ *o_mag_vec;
   }
 }

 __float128 *arbitrary_rotation_matrix_3d(__float128 i_theta, __float128 *i_rot_centre, __float128 *i_rot_axis, __float128 *i_point, __float128 *o_result){
   // Rotates i_point an angle of i_theta about the unit vector i_rot_axis passing through the point i_rot_centre..
   __float128 u_sq, v_sq, w_sq, au, bv, cw, m_ux_m_vy_m_wz, costheta, one_m_costheta, sintheta;
   __float128 mag_rot_axis;
   // Always assume the user is stupid and check whether i_rot_axis is normalised, if it's not normalise it.
   mag_rot_axis = dot_product(i_rot_axis, i_rot_axis, 3);
   if(mag_rot_axis != 1.0q){
     mag_rot_axis = sqrtq(mag_rot_axis);
     for (int i = 0; i < 3; i++){
       i_rot_axis[i] /= mag_rot_axis;
     }
   }
   // cosq(i_theta), 1 - cosq(i_theta), sinq(i_theta)
   sintheta = sinq(i_theta);
   costheta = cosq(i_theta);
   one_m_costheta = 1.0q - costheta;
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
/*
 __float128 *arbitrary_rotation_matrix_3d(__float128 theta, __float128 *abc, __float128 *uvw, __float128 *xyz, __float128 *o_xpypzp){
   // Rotates point xyz an angle of theta about the unit vector uvw passing through the point abc.
   __float128 costheta       = cosq(theta);
   __float128 one_m_costheta = 1. - costheta;
   __float128 mag_uvw = dot_product(uvw, uvw, 3);
   if(mag_uvw != 1.0q){
     normalise_vector(uvw, 3, uvw);
   }
   o_xpypzp[0] = (abc[0]*(uvw[1]*uvw[1] + uvw[2]*uvw[2]) - uvw[0]*(abc[1]*uvw[1] + abc[2]*uvw[2] - uvw[0]*xyz[0] - uvw[1]*xyz[1] - uvw[2]*xyz[2]))*one_m_costheta + xyz[0]*costheta + (-abc[2]*uvw[1] + abc[1]*uvw[2] - uvw[2]*xyz[1] + uvw[1]*xyz[2])*sinq(theta);
   o_xpypzp[1] = (abc[1]*(uvw[0]*uvw[0] + uvw[2]*uvw[2]) - uvw[1]*(abc[0]*uvw[0] + abc[2]*uvw[2] - uvw[0]*xyz[0] - uvw[1]*xyz[1] - uvw[2]*xyz[2]))*one_m_costheta + xyz[1]*costheta + ( abc[2]*uvw[0] - abc[0]*uvw[2] + uvw[2]*xyz[0] - uvw[0]*xyz[2])*sinq(theta);
   o_xpypzp[2] = (abc[2]*(uvw[0]*uvw[0] + uvw[1]*uvw[1]) - uvw[2]*(abc[0]*uvw[0] + abc[1]*uvw[1] - uvw[0]*xyz[0] - uvw[1]*xyz[1] - uvw[2]*xyz[2]))*one_m_costheta + xyz[2]*costheta + (-abc[1]*uvw[0] + abc[0]*uvw[1] - uvw[1]*xyz[0] + uvw[0]*xyz[1])*sinq(theta);
   //printf("%Qf, %Qf, %Qf\n", xpypzp[0], xpypzp[1], xpypzp[2]);
   return o_xpypzp;
 }*/

 __float128 *build_vector(__float128 *i_x1, __float128 *i_x2, int i_vec_size,
                      __float128 *o_vec){
   // Returns a vector o_vec which translates the point i_x1 to i_x2.
   for (int i = 0; i < i_vec_size; i++){
     o_vec[i] = i_x2[i] - i_x1[i];
   }
   return o_vec;
 }

 __float128 *init_vector(__float128 *i_x1, __float128 *i_x2, int i_vec_size,
                     __float128 *io_vec){
   // Builds and returns a normalised vector io_vect.
   build_vector(i_x1, i_x2, i_vec_size, io_vec);
   normalise_vector(io_vec, i_vec_size, io_vec);
   return io_vec;
 }

 void init_vector2(__float128 *i_x1  , __float128 *i_x2, int i_vec_size,
                   __float128 *io_vec, __float128 *o_mag_vec){
   // Builds and returns a normalised vector io_vect, and the magnitude of the vector in o_mag_vec.
   // Has to be a void function in order to 'return' two values, the magnitude should be passed as a reference eg:
   /*
     int    size = 3;
     __float128 i_x1[size], i_x2[size], o_vec[size], magnitude;
     init_vector2(i_x1, i_x2, size, o_vec, &magnitude);
   */
   build_vector(i_x1, i_x2, i_vec_size, io_vec);
   normalise_vector2(io_vec, i_vec_size, io_vec, o_mag_vec);
 }

 __float128 init_point(__float128 *i_vec1, __float128 *i_vec2,
                   __float128 *i_vec3, __float128 *i_vec4,
                   int     i_vec_size){
   // Initialises points on the surface element given four vectors.
   __float128 result = 0.0q;
   __float128 denom = 0.0q;
   denom = dot_product(i_vec3, i_vec4, i_vec_size);
   if (denom == 0.0q){
     fprintf(stderr, "nodal_surface_force_quad_triangle: init_point: division by zero, dot_product(i_vec3, i_vec4) = %Qf\n", denom);
     exit(EXIT_FAILURE);
   }
   result = dot_product(i_vec1, i_vec2, i_vec_size)/denom;
   return result;
 }

 __float128 seed_single_integral(__float128 a, __float128 b){
   // Single integral seed.
   __float128 result = 0.0q;
   __float128 arg;
   arg = a+b;
   if (arg <= 0.0q){
     fprintf(stderr, "nodal_surface_force_quad_triangle: seed_single_integral: logq(%Qf) = %Qf\n", arg, logq(arg));
     exit(EXIT_FAILURE);
   }
   result = logq(arg);
   return result;
 }

 void init_force(__float128 *nodal_force[n_nodes], __float128 *total_force){
   // Sets forces to zero.
   for (int i = 0; i < n_nodes; i++){
     for (int j = 0; j < 3; j++){
       nodal_force[i][j] = 0.0q;
     }
   }
   for (int i = 0; i < 3; i++){
     total_force[i] = 0.0q;
   }
 }

 void add_force(__float128 *p_nodal_force[n_nodes], __float128 *p_total_force, __float128 *nodal_force[n_nodes], __float128 *total_force){
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

 void mean_force(__float128 *nodal_force[n_nodes], __float128 *total_force, int n_samples){
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

 void integral_vector(__float128 *i_p, __float128 *i_q, __float128 *i_b, __float128 *i_t, __float128 *i_n,
                      __float128 i_one_m_nu, __float128 i_a_sq,
                      __float128 o_vec_int[][3]){
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
  __float128 t_x_b[3], p_x_b[3],
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
  __float128 t_dot_n,
         p_x_b_dot_n, q_x_b_dot_n,
         t_x_b_dot_n, b_x_t_dot_n,
         p_x_b_dot_t, q_x_b_dot_t,
         t_dot_n_p_x_b_dot_t, t_dot_n_q_x_b_dot_t,
         t_dot_n_p_x_b_dot_t_3, t_dot_n_q_x_b_dot_t_3,
         one_m_nu_1p5_a_sq, a_sq_3;

   // Calculate auxiliary local variables.
   one_m_nu_1p5_a_sq = 1.5q * i_one_m_nu * i_a_sq;
   a_sq_3 = 3.0q * i_a_sq;
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
   t_dot_n_p_x_b_dot_t_3 = 3.0q * t_dot_n_p_x_b_dot_t;
   t_dot_n_q_x_b_dot_t_3 = 3.0q * t_dot_n_q_x_b_dot_t;

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

 __float128 *integrals_quad_triangle(__float128 *rv , __float128 *sv , __float128 *yv,
                                 __float128 *rfs, __float128 *sfr,
                                 __float128 r1, __float128 s1,
                                 __float128 g, __float128 m,
                                 __float128 f, __float128 h,
                                 __float128 c, __float128 d, __float128 e,
                                 __float128 a2,
                                 __float128 *o_sch, int i_num_integrals){
   // Sign vector for quick evaluation of integrals via the dot product.
   static __float128 signv[8] = {1.0,-1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0};
   __float128 rv2[8], sv2[8], yv2[8];
   __float128 ra[8], rra[8], sra[8], rar1[8], rar2[8], ras1[8], ras2[8];
   __float128 rdotp[8], rdotq[8], rrdott[8], srdott[8];
   __float128 rardotp[8], rardotq[8], rrardott[8], srardott[8];
   __float128 sfr_sq;
   // integrals
   __float128 a01_s1[8], a01_s2[8], b01_r1[8], b01_r2[8];
   __float128 c01_s1[8], c01_s2[8], c01_r1[8], c01_r2[8];
   __float128 a11_s1[8], a11_s2[8], b11_r1[8], b11_r2[8];
   __float128 c11_s1[8], c11_s2[8], c11_r1[8], c11_r2[8];
   __float128 a0m1_s1[8], a0m1_s2[8], b0m1_r1[8], b0m1_r2[8];
   __float128 c0m1_s1[8], c0m1_s2[8], c0m1_r1[8], c0m1_r2[8];
   __float128 a1m1_s1[8], a1m1_s2[8], b1m1_r1[8], b1m1_r2[8];
   __float128 aa, bb, cc;
   __float128 a21_s1[8], a21_s2[8], b21_r1[8], b21_r2[8];
   __float128 root, tp1, tp2, tp3, tp4, tp5;
   __float128 e003_s1[8], e003_s2[8], f003_r1[8], f003_r2[8];
   __float128 d003[8], d003_r1[8], d003_r2[8], d003_s1[8], d003_s2[8];
   __float128 e103_s1[8], e103_s2[8], f103_r1[8], f103_r2[8];
   __float128 d103[8],d013[8];
   __float128 aas001[8], bbr001[8];
   __float128 e001_s1[8], e001_s2[8], e101_s1[8], e101_s2[8], e011_s1[8], e011_s2[8];
   __float128 f001_r1[8], f001_r2[8], f101_r1[8], f101_r2[8], f011_r1[8], f011_r2[8];
   __float128 d001[8], d011[8], d101[8];
   __float128 aas00m1[8], bbr00m1[8];
   __float128 e00m1_s1[8], e00m1_s2[8], f00m1_r1[8], f00m1_r2[8];
   __float128 d00m1[8];
   __float128 f201_r1[8], f201_r2[8], e201_s1[8], e201_s2[8];
   __float128 f111_r1[8], f111_r2[8], e111_s1[8], e111_s2[8];
   __float128 aas10m1[8], aas01m1[8];
   __float128 bbr01m1[8], bbr10m1[8];
   __float128 d201[8], d021[8], d111[8], d111_1[8], d111_2[8];
   __float128 e203_s1[8], e203_s2[8], f203_r1[8], f203_r2[8];
   __float128 aas101[8], aas011[8];
   __float128 bbr101[8], bbr011[8];
   __float128 d203[8], d023[8], d113[8], d113_1[8], d113_2[8];
   __float128 e013_s1[8], e013_s2[8], f013_r1[8], f013_r2[8];
   __float128 e113_s1[8], e113_s2[8], f113_r1[8], f113_r2[8];
   __float128 e213_s1[8], e213_s2[8], f213_r1[8], f213_r2[8];
   __float128 aas111[8], bbr111[8], d213[8], d123[8];
   __float128 f313_r1[8], f313_r2[8];
   __float128 e313_s1[8], e313_s2[8];
   __float128 aas201[8], aas021[8] ,bbr201[8], bbr021[8] ;
   __float128 d303[8], d033[8];
   __float128 h0001[8], ffr1001[8], ees0101[8];
   __float128 ees000m1[8], ffr000m1[8], h0011[8], h1001[8], h0101[8];
   __float128 h1103_1[8], h1103_2[8];
   __float128 ffr0001[8], ees0001[8], ffr0101[8], ees1001[8];
   __float128 ees2001[8], ffr0201[8], ffr2001[8], ees0201[8], ees1101[8], ffr1101[8];
   __float128 h1105_1[8], h1105_2[8];
   __float128 ees0003[8], ffr0003[8], ees1003[8], ffr1003[8], ees0103[8], ffr0103[8];
   __float128 ees1103[8], ffr1103[8], ees1013[8], ffr1013[8], ees0113[8];
   __float128 ffr0113[8], ees2013[8], ffr2013[8], ees0213[8], ffr0213[8];
   __float128 ees1113[8], ffr1113[8];
   __float128 h2215_1[8], h2215_2[8];
   __float128 ees2003[8], ffr0203[8], ees0203[8], ffr2003[8], ees3013[8], ffr0313[8], ffr3013[8], ees0313[8];
   __float128 ees2113[8], ffr1213[8], ffr2113[8], ees1213[8];
   // Triple integrals.
   __float128 hijkl[i_num_integrals][n_limits];
   // Indices for the evaluation of integrals.
   int ir1[4],ir2[4],is1[4],is2[4];
   int i;
   // Auxiliary scalars.
   __float128 c2, d2, f2, g2, e2;
   __float128 ce, de, ef, fg, eg, cf, cg, onepf2p2ef, oneph2p2eh;
   __float128 h2, m2;
   __float128 cd, hm, dh, eh, em, dm;
   __float128 f2g, fg2, h2m, hm2, f3, h3, g3, m3, cde;
   __float128 r12, r13, s12, s13;
   // Auxiliary scalars
   c2 = c*c; d2 = d*d; e2 = e*e; h2=h*h; f2=f*f; g2=g*g; m2=m*m; eh=e*h; ef=e*f;
   dh=h*d; cf=f*c; cd=c*d; cg=c*g; fg=f*g; eg=e*g; em=e*m; hm=h*m; dm=d*m;
   de=d*e; ce=c*e; f3=f2*f; h3=h2*h; m3=m2*m; g3=g2*g; h2m=h2*m; hm2=h*m2;
   f2g=f2*g; fg2=f*g2; cde= cd*e;
   r12 = r1*r1; r13 = r12*r1; s12 = s1*s1; s13 = s12*s1;
   onepf2p2ef = 1.0q+f2+2.0q*ef; oneph2p2eh = 1.0q+h2+2.0q*eh;
   for(int i = 0; i < n_limits; i++){
     rv2[i] = rv[i]*rv[i];
     sv2[i] = sv[i]*sv[i];
     yv2[i] = yv[i]*yv[i];
     // ra = sqrtq(r.r)
     ra      [i] = sqrtq(yv2[i]+rv2[i]+sv2[i]+2.0q*c*rv[i]*yv[i]+2.0q*d*sv[i]*yv[i]+2.0q*e*rv[i]*sv[i]+a2);
     rra     [i] = sqrtq(yv2[i]+rv2[i]+sfr[i]*sfr[i]+2.0q*c*rv[i]*yv[i]+2.0q*d*sfr[i]*yv[i]+2.0q*e*rv[i]*sfr[i]+a2);
     sra     [i] = sqrtq(yv2[i]+rfs[i]*rfs[i]+sv2[i]+2.0q*c*rfs[i]*yv[i]+2.0q*d*sv[i]*yv[i]+2.0q*e*rfs[i]*sv[i]+a2);
     rar1    [i] = sqrtq(yv2[i]+r1*r1+sv2[i]+2.0q*c*r1*yv[i]+2.0q*d*sv[i]*yv[i]+2.0q*e*r1*sv[i]+a2);
     rar2    [i] = sqrtq(yv2[i]+(onepf2p2ef)*sv2[i]+2.0q*sv[i]*yv[i]*(cf+d)+2.0q*sv[i]*(fg+eg)+2.0q*cg*yv[i]+g2+a2);
     ras1    [i] = sqrtq(yv2[i]+rv2[i]+s1*s1+2.0q*c*rv[i]*yv[i]+2.0q*d*s1*yv[i]+2.0q*e*rv[i]*s1+a2);
     ras2    [i] = sqrtq(yv2[i]+(oneph2p2eh)*rv2[i]+2.0q*rv[i]*yv[i]*(dh+c)+2.0q*rv[i]*(hm+em)+2.0q*dm*yv[i]+m2+a2);
     rdotp   [i] = rv [i] + c*yv [i] + e*sv[i];
     rdotq   [i] = sv [i] + d*yv [i] + e*rv[i];
     rrdott  [i] = yv [i] + c*rv [i] + d*sfr[i];
     srdott  [i] = yv [i] + c*rfs[i] + d*sv[i];
     rardotp [i] = ra [i] + rdotp[i];
     rardotq [i] = ra [i] + rdotq[i];
     rrardott[i] = rra[i] + rrdott[i];
     srardott[i] = sra[i] + srdott[i];
}
   //
   // LINEAR INTEGRALS
   //
   /* The following shorthand is introduced to simplify notation:
      Ail = /int r^i Ra^-l dr
      Bjl = /int s^j Ra^-l ds
      Ckl = /int y^k Ra^-l dy
      Linear integrals are the last object created by the sequence of
      integration by part. since r2 depends upon s and s2 depends up on
      r, 2 different sollutions must be provided for Ail_s1, Ail_s2  and
      Bjl_r1 and Bjl_r2 and 4 different solutions for Ckl.
   */

   /* Since integration paths are different for r1 and r2,
      and for s1 and s2, we need to evaluate integrals Ail_r1
      at vector components corresponding to r1 */
   ir1[0]=4;ir1[1]=5;ir1[2]=6;ir1[3]=7;
   is1[0]=2;is1[1]=3;is1[2]=6;is1[3]=7;
   ir2[0]=0;ir2[1]=1;ir2[2]=2;ir2[3]=3;
   is2[0]=0;is2[1]=1;is2[2]=4;is2[3]=5;

   //
   // LINEAR SEEDS ABC01
   //
   for(int j = 0; j < 4; j++) {
     // r1
     i=ir1[j];
     b01_r1[i]= logq( fabsq(rardotq[i]));
     c01_r1[i]= logq( fabsq(rrardott[i]));
     // abc11
     b11_r1[i]=sra[i]-(d*yv[i]+e*rv[i])*b01_r1[i];
     c11_r1[i]=rar1[i]-(c*rv[i]+d*sfr[i])*c01_r1[i];
     // abc1m1
     b0m1_r1[i]=0.5q*(rdotq[i]*ra[i] +(ra[i]*ra[i]-(rdotq[i]*rdotq[i]))*b01_r1[i]);
     b1m1_r1[i]=1.0q/3.0q*(ra[i]*ra[i]*ra[i])-(d*yv[i]+e*rv[i])*b0m1_r1[i];
     c0m1_r1[i]=0.5q*(rrdott[i]*rra[i] +(rra[i]*rra[i]-(rrdott[i]*rrdott[i]))*c01_r1[i]);
     // abc0m3
     sfr_sq=sfr[i]*sfr[i];
     // abc21
     b21_r1[i]= -1.0q*(-sv[i]*rar1[i]+b0m1_r1[i]-(d*yv[i]+e*r1)*(-1.0q)*b11_r1[i]);
     // abc31
     root= 1.0q/sqrtq((1.0q-d2)*a2+((1.0q-d2)*(1.0q-e2)-(c-de)*(c-de))*rv2[i]);
     f003_r1[i]= 2.0q*root*atanq(((1.0q-d)*(rra[i]-sv[i]-d*yv[i]-e*rv[i])+(1.0q-d2)*yv[i]+(c-de)*rv[i])*root);
     root= 1.0q/sqrtq((1.0q-e2)*a2+(1.0q-c2-d2-e2+2.0q*c*de)*yv2[i]);
     d003_r1[i]=2.0q*root*atanq(((1.0q-e)*(rar1[i]+rv[i]-sv[i])+(c-d)*yv[i])*root);
     f103_r1[i]=1.0q/(1.0q-d2)*((+d*b01_r1[i] -c01_r1[i]) +(cd-e)*rv[i]*f003_r1[i]);
     f013_r1[i]=1.0q/(1.0q-d2)*((+d*c01_r1[i] -b01_r1[i]) +(de-c)*rv[i]*f003_r1[i]);
     f001_r1[i]=1.0q/(1.0q-d2)/(1.0q-2.0q)*(((a2+rv2[i])*(1.0q-d2)+2.0q*cd*e*rv2[i]-1.0q*rv2[i]*(e2+c2))*f003_r1[i]+(-(sfr[i]*(1.0q-d2)-c*rv[i]*d+e*rv[i])*c01_r1[i]-(yv[i]*(1.0q-d2)+c*rv[i]-e*d*rv[i])*b01_r1[i]));
     f101_r1[i]=1.0q/(1.0q-d2)*(1.0q/(1.0q-2.0q)*(+d*b0m1_r1[i]-c0m1_r1[i])+(cd-e)*rv[i]*f001_r1[i]);
     f011_r1[i]=1.0q/(1.0q-d2)*(1.0q/(1.0q-2.0q)*(+d*c0m1_r1[i]-b0m1_r1[i])+(de-c)*rv[i]*f001_r1[i]);
     f00m1_r1[i]=1.0q/(1.0q-d2)/(-3.0q)*((-(a2+rv2[i])*(1.0q-d2)+2.0q*cd*e*-rv2[i]+rv2[i]*(e2+c2))*f001_r1[i]+(-(sfr[i]*(1.0q-d2)-c*rv[i]*d+e*rv[i])*c0m1_r1[i]-(yv[i]*(1.0q-d2)+c*rv[i]-e*d*rv[i])*b0m1_r1[i]));
     f201_r1[i]=1.0q/(1.0q-d2)*(1.0q/(-1.0q)*(f00m1_r1[i]+d*b1m1_r1[i]-sfr[i]*c0m1_r1[i])+(cd-e)*rv[i]*f101_r1[i]);
     f111_r1[i]=1.0q/(1.0q-d2)*(1.0q/(-1.0q)*(-d*f00m1_r1[i]+d*sfr[i]*c0m1_r1[i]-b1m1_r1[i])+(de-c)*rv[i]*f101_r1[i]);
     f203_r1[i]=1.0q/(1.0q-d2)*((f001_r1[i]+d*b11_r1[i]-sfr[i]*c01_r1[i])+(cd-e)*rv[i]*f103_r1[i]);
     f113_r1[i]=1.0q/(1.0q-d2)*((-d*f001_r1[i]+d*yv[i]*b01_r1[i]-c11_r1[i])+(cd-e)*rv[i]*f013_r1[i]);
     f213_r1[i]=1.0q/(1.0q-d2)*((f011_r1[i]-d*f101_r1[i]+d*yv[i]*b11_r1[i]-sfr[i]*c11_r1[i])+(cd-e)*rv[i]*f113_r1[i]);
     f313_r1[i]=1.0q/(1.0q-d2)*((2.0q*f111_r1[i]-d*f201_r1[i]+d*yv[i]*b21_r1[i]-sfr_sq*c11_r1[i])+(cd-e)*rv[i]*f213_r1[i]);
     // r2
     aa=(onepf2p2ef);
     i=ir2[j];
     bb=(yv[i]*(cf+d)+g*(f+e));
     cc=(yv2[i]+2.0q*cg*yv[i]+g2+a2);
     b01_r2[i]=logq(2.0q*sqrtq(aa)*sqrtq(aa*sv2[i]+2.0q*bb*sv[i]+cc)+2.0q*(aa*sv[i]+bb))/sqrtq(aa);
     bb=cg+(cf+d)*sv[i];
     cc=sqrtq(-bb*bb+sv2[i]*(onepf2p2ef)+2.0q*sv[i]*(fg+eg)+g2+a2);
     c01_r2[i]=logq(2.0q*sqrtq((yv[i]+bb)*(yv[i]+bb)+cc*cc)+2.0q*(yv[i]+bb));
     // linear integrals.
     b11_r2[i]=1.0q/(onepf2p2ef)*(sra[i]-(yv[i]*(cf+d)+fg+eg)*b01_r2[i]);
     c11_r2[i]=rar2[i]-(sv[i]*(cf+d)+cg)*c01_r2[i];
     b0m1_r2[i]= 0.5q*(sv[i]*sra[i]+(yv[i]*(cf+d)+g*(f+e))*b11_r2[i]+(yv2[i]+2.0q*g*c*yv[i]+g2+a2)*b01_r2[i]);
     b1m1_r2[i]=1.0q/(onepf2p2ef)*(1.0q/3.0q*sra[i]*sra[i]*sra[i]-(yv[i]*(cf+d)+g*(f+e))*b0m1_r2[i]);
     c0m1_r2[i]=0.5q*(+rar2[i]*(yv[i]+sv[i]*(cf+d)+cg)+((onepf2p2ef-(cf+d)*(cf+d))*sv2[i]+2.0q*sv[i]*(fg+eg-cg*(cf+d))+(1.0q-c2)*g2+a2)*c01_r2[i]);
     b21_r2[i]= -1.0q/(onepf2p2ef)*(-sv[i]*rar2[i]+b0m1_r2[i]-(yv[i]*(cf+d)+fg+eg)*(-1.0q)*b11_r2[i]);
     // Double integrals
     tp1=sqrtq(onepf2p2ef);
     tp2=(yv[i]*(cf+d)+fg+eg)/tp1;
     tp3=sqrtq(tp1-cf-d);
     tp4=(tp1*(yv[i]+cg)-tp2*(cf+d))/tp3;
     tp5=sqrtq(-tp4*tp4-(tp1+cf+d)*tp2*tp2+(tp1+cf+d)*(a2+g2+yv2[i]+2.0q*cg*yv[i]));
     f003_r2[i]=2.0q/tp3/tp5*atanq((tp3*(rar2[i]-tp1*sv[i]-tp2)+tp4)/tp5);
     tp1=sqrtq(onepf2p2ef);
     tp2=sqrtq(tp1-f-e);
     root=1.0q/sqrtq(yv2[i]*(1.0q-c2-d2-e2+2.0q*c*de)+a2*(1.0q-e2));
     d003_r2[i]=2.0q*root*atanq((tp2*tp2*rar2[i]+(tp1*c-cf-d)*yv[i]+tp1*g-fg-eg-tp1*tp2*tp2*sv[i])*root);
     f103_r2[i]=1.0q/(onepf2p2ef -(cf+d)*(cf+d))*((cf+d)*b01_r2[i] -c01_r2[i] -(fg+eg-cg*(cf+d))*f003_r2[i]);
     f013_r2[i]=1.0q/(onepf2p2ef -(cf+d)*(cf+d))*((cf+d)*c01_r2[i] -(onepf2p2ef)*b01_r2[i] +((cf+d)*(fg+eg)-(onepf2p2ef)*cg)*f003_r2[i]);
     tp1=onepf2p2ef-(cf+d)*(cf+d);
     f001_r2[i] = 1.0q/(1.0q-2.0q)*(((fg+eg)*(cf+d)-cg*(onepf2p2ef))/tp1*b01_r2[i]-yv[i]*b01_r2[i] +(cg*(cf+d)-(fg+eg))/tp1*c01_r2[i]-sv[i]*c01_r2[i] +(a2+g2)*f003_r2[i] - ((fg+eg)*(fg+eg-cg*(cf+d))+cg*(cg*(onepf2p2ef)-(cf+d)*(fg+eg)))/tp1*f003_r2[i]);
     f101_r2[i]=1.0q/(1.0q-2.0q)/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*b0m1_r2[i]-c0m1_r2[i]-(1.0q-2.0q)*(fg+eg-cg*(cf+d))*f001_r2[i]);
     f011_r2[i]=1.0q/(1.0q-2.0q)/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*c0m1_r2[i]-(onepf2p2ef)*b0m1_r2[i]+(1.0q-2.0q)*((cf+d)*(fg+eg)-(onepf2p2ef)*cg)*f001_r2[i]);
     tp1=onepf2p2ef-(cf+d)*(cf+d);
     f00m1_r2[i] = 1.0q/(-3.0q)*(((fg+eg)*(cf+d)-cg*(onepf2p2ef))/tp1*b0m1_r2[i]-yv[i]*b0m1_r2[i] +(cg*(cf+d)-(fg+eg))/tp1*c0m1_r2[i]-sv[i]*c0m1_r2[i] -1.0q*(a2+g2)*f001_r2[i] +1.0q*((fg+eg)*(fg+eg-cg*(cf+d))+cg*(cg*(onepf2p2ef)-(cf+d)*(fg+eg)))/tp1*f001_r2[i]);
     tp1=onepf2p2ef-(cf+d)*(cf+d);
     f201_r2[i]=1.0q/(-1.0q)/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*b1m1_r2[i]-sv[i]*c0m1_r2[i]+1.0q*f00m1_r2[i]-(-1.0q)*(fg+eg-cg*(cf+d))*f101_r2[i]);
     f111_r2[i]=1.0q/(-1.0q)/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*sv[i]*c0m1_r2[i]-(onepf2p2ef)*b1m1_r2[i]-(cf+d)*f00m1_r2[i]+(-1.0q)*((cf+d)*(fg+eg)-(onepf2p2ef)*cg)*f101_r2[i]);
     f203_r2[i]=1.0q/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*b11_r2[i]-sv[i]*c01_r2[i]+f001_r2[i]-(fg+eg-cg*(cf+d))*f103_r2[i]);
     f113_r2[i]=1.0q/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*yv[i]*b01_r2[i]-c11_r2[i]-(cf+d)*f001_r2[i]-(fg+eg-cg*(cf+d))*f013_r2[i]);
     f213_r2[i]=1.0q/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*yv[i]*b11_r2[i]-sv[i]*c11_r2[i]+f011_r2[i]-(cf+d)*f101_r2[i]-(fg+eg-cg*(cf+d))*f113_r2[i]);
     f313_r2[i]=1.0q/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*yv[i]*b21_r2[i]-powq(sv[i],2.0q)*c11_r2[i]+2.0q*f111_r2[i]-(cf+d)*f201_r2[i]-(fg+eg-cg*(cf+d))*f213_r2[i]);
     // s1
     i=is1[j];
     a01_s1[i]= logq(fabsq(rardotp[i]));
     c01_s1[i]= logq(fabsq(srardott[i]));
     a11_s1[i]=rra[i]-(c*yv[i]+e*sv[i])*a01_s1[i];
     c11_s1[i]=ras1[i]-(d*sv[i]+c*rfs[i])*c01_s1[i];
     a0m1_s1[i]=0.5q*(rdotp[i]*ra[i] +(ra[i]*ra[i]-(rdotp[i]*rdotp[i]))*a01_s1[i]);
     a1m1_s1[i]=-1.0q/3.0q*(-ras1[i]*ras1[i]*ras1[i])-(c*yv[i]+e*sv[i])*a0m1_s1[i];
     c0m1_s1[i]=0.5q*(srdott[i]*sra[i] +(sra[i]*ra[i]-(srdott[i]*srdott[i]))*c01_s1[i]);
     a21_s1[i]= -1.0q*(-rv[i]*ras1[i]+a0m1_s1[i]-(c*yv[i]+e*s1)*(-1.0q)*a11_s1[i]);
     root= 1.0q/sqrtq((1.0q-c2)*a2+((1.0q-c2)*(1.0q-e2)-(d-ce)*(d-ce))*sv2[i]);
     e003_s1[i]=2.0q*root*atanq(((1.0q-c)*(sra[i]-rv[i]-c*yv[i]-e*sv[i])+(1.0q-c2)*yv[i]+(d-c*e)*sv[i])*root);
     root= 1.0q/sqrtq((1.0q-e2)*a2+(1.0q-c2-d2-e2+2.0q*c*de)*yv2[i]);
     d003_s1[i]=2.0q*root*atanq(((1.0q-e)*(ras1[i]+sv[i]-rv[i])+(d-c)*yv[i])*root);
     e103_s1[i]=1.0q/(1.0q-c2)*((+c*a01_s1[i] -c01_s1[i]) +(cd-e)*sv[i]*e003_s1[i]);
     e013_s1[i]=1.0q/(1.0q-c2)*((+c*c01_s1[i] -a01_s1[i]) +(ce-d)*sv[i]*e003_s1[i]);
     e001_s1[i] = 1.0q/(1.0q-c2)/(1.0q-2.0q)*(((a2+sv2[i])*(1.0q-c2)+2.0q*c*d*e*sv2[i]-sv2[i]*(e2+d2))*e003_s1[i]+(-(rfs[i]*(1.0q-c2)-d*sv[i]*c+e*sv[i])*c01_s1[i]-(yv[i]*(1.0q-c2)+d*sv[i]-e*c*sv[i])*a01_s1[i]));
     e101_s1[i]=1.0q/(1.0q-c2)*(1.0q/(1.0q-2.0q)*(+c*a0m1_s1[i]-c0m1_s1[i])+(cd-e)*sv[i]*e001_s1[i]);
     e011_s1[i]=1.0q/(1.0q-c2)*(1.0q/(1.0q-2.0q)*(+c*c0m1_s1[i]-a0m1_s1[i])+(c*e-d)*sv[i]*e001_s1[i]);
     e00m1_s1[i] = 1.0q/(1.0q-c2)/(-3.0q)*((-(a2+sv2[i])*(1.0q-c2)+2.0q*cd*e*-sv2[i]+sv2[i]*(e2+d2))*e001_s1[i]+(-(rfs[i]*(1.0q-c2)-d*sv[i]*c+e*sv[i])*c0m1_s1[i]-(yv[i]*(1.0q-c2)+d*sv[i]-e*c*sv[i])*a0m1_s1[i]));
     e201_s1[i]=1.0q/(1.0q-c2)*(1.0q/(-1.0q)*(e00m1_s1[i]+c*a1m1_s1[i]-rfs[i]*c0m1_s1[i])+(cd-e)*sv[i]*e101_s1[i]);
     e111_s1[i]=1.0q/(1.0q-c2)*(1.0q/(1.0q-2.0q)*(-c*e00m1_s1[i]+c*rfs[i]*c0m1_s1[i]-a1m1_s1[i])+(ce-d)*sv[i]*e101_s1[i]);
     e203_s1[i]=1.0q/(1.0q-c2)*((e001_s1[i]+c*a11_s1[i]-rfs[i]*c01_s1[i])+(cd-e)*sv[i]*e103_s1[i]);
     e113_s1[i]=1.0q/(1.0q-c2)*((-c*e001_s1[i]+c*yv[i]*a01_s1[i]-c11_s1[i])+(cd-e)*sv[i]*e013_s1[i]);
     e213_s1[i]=1.0q/(1.0q-c2)*((e011_s1[i]-c*e101_s1[i]+c*yv[i]*a11_s1[i]-rfs[i]*c11_s1[i])+(cd-e)*sv[i]*e113_s1[i]);
     e313_s1[i]=1.0q/(1.0q-c2)*((2.0q*e111_s1[i]-c*e201_s1[i]+c*yv[i]*a21_s1[i]-powq(rfs[i],2.0q)*c11_s1[i])+(cd-e)*sv[i]*e213_s1[i]);
     // s2
     aa=(oneph2p2eh);
     i=is2[j];
     bb=(yv[i]*(h*d+c)+m*(h+e));
     cc=(yv2[i]+2.0q*m*d*yv[i]+m2+a2);
     a01_s2[i]=logq(2.0q*sqrtq(aa)*sqrtq(aa*rv2[i]+2.0q*bb*rv[i]+cc)+2.0q*(aa*rv[i]+bb))/sqrtq(aa);
     bb=dm+(dh+c)*rv[i];
     cc=sqrtq(-bb*bb+rv2[i]*(oneph2p2eh)+2.0q*rv[i]*(hm+em)+m2+a2);
     c01_s2[i]=logq(2.0q*sqrtq((yv[i]+bb)*(yv[i]+bb)+cc*cc)+2.0q*(yv[i]+bb));
     a11_s2[i]=1.0q/(oneph2p2eh)*(rra[i]-(yv[i]*(dh+c)+hm+em)*a01_s2[i]);
     c11_s2[i]=ras2[i]-(rv[i]*(dh+c)+dm)*c01_s2[i];
     a0m1_s2[i]= 0.5q*(rv[i]*rra[i]+(yv[i]*(dh+c)+m*(h+e))*a11_s2[i]+(yv2[i]+2.0q*m*d*yv[i]+m2+a2)*a01_s2[i]);
     a1m1_s2[i]=1.0q/(oneph2p2eh)*(1.0q/3.0q*rra[i]*rra[i]*rra[i]-(yv[i]*(dh+c)+m*(h+e))*a0m1_s2[i]);
     c0m1_s2[i]=0.5q*(+ras2[i]*(yv[i]+rv[i]*(dh+c)+dm)+((oneph2p2eh-(dh+c)*(dh+c))*rv2[i]+2.0q*rv[i]*(hm+em-dm*(dh+c))+(1.0q-d2)*m2+a2)*c01_s2[i]);
     a21_s2[i]= -1.0q/(oneph2p2eh)*(-rv[i]*ras2[i]+a0m1_s2[i]-(yv[i]*(dh+c)+hm+em)*(-1.0q)*a11_s2[i]);
     //
     // DOUBLE INTEGRALS
     //
     /* DOUBLE integrals are built from previous linear integrals and a set
        of 3 seed functions.
        Eikl = /int/int r^i y^k Ra^-l dr dy
        Fjkl = /int/int s^j y^k Ra^-l ds dy
        Dijl = /int/int r^i s^j Ra^-l dr ds
        For Eikl and Fjkl integrals, integration paths are different for r1,r2 or
        s1,s2 as consequence of s2 = h*r+m and r2= f*s+g.
        We introduce also the following linear integrals:
        AASijl = /int r^i s(r)^j Ra^-l dr = /int r^i (h*r+m)^j Ra^-l dr - /int r^i s1^j Ra^-l dr
        BBrijl = /int r(s)^i s^j Ra^-l ds = /int (f*s+g)^i s^j Ra^-l ds - /int r1^i s^j Ra^-l ds
     */
     //
     // SEED FUNCTIONS: DEF003
     tp1=sqrtq(oneph2p2eh);
     tp2=(yv[i]*(dh+c)+hm+em)/tp1;
     tp3=sqrtq(tp1-dh-c);
     tp4=(tp1*(yv[i]+dm)-tp2*(dh+c))/tp3;
     tp5=sqrtq(-tp4*tp4-(tp1+dh+c)*tp2*tp2+(tp1+dh+c)*(a2+m2+yv2[i]+2.0q*dm*yv[i]));
     e003_s2[i]=2.0q/tp3/tp5*atanq((tp3*(ras2[i]-tp1*rv[i]-tp2)+tp4)/tp5);
     tp1=sqrtq(oneph2p2eh);
     tp2=sqrtq(tp1-h-e);
     root=1.0q/sqrtq(yv2[i]*(1.0q-c2-d2-e2+2.0q*c*de)+a2*(1.0q-e2));
     d003_s2[i]=2.0q*root*atanq((tp2*tp2*ras2[i]+(tp1*d-dh-c)*yv[i]+tp1*m-hm-em-tp1*tp2*tp2*rv[i])*root);
     e103_s2[i]=1.0q/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*a01_s2[i] -c01_s2[i] -(hm+em-dm*(dh+c))*e003_s2[i]);
     e013_s2[i]=1.0q/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*c01_s2[i] -(oneph2p2eh)*a01_s2[i] +((dh+c)*(hm+em) -(oneph2p2eh)*dm)*e003_s2[i]);
     tp1=oneph2p2eh-(dh+c)*(dh+c);
     e001_s2[i] = 1.0q/(1.0q-2.0q)*(((hm+em)*(dh+c)-dm*(oneph2p2eh))/tp1*a01_s2[i]-yv[i]*a01_s2[i] +(dm*(dh+c)-(hm+em))/tp1*c01_s2[i]-rv[i]*c01_s2[i] +(a2+m2)*e003_s2[i] - ((hm+em)*(hm+em-dm*(dh+c))+dm*(dm*(oneph2p2eh)-(dh+c)*(hm+em)))/tp1*e003_s2[i]);
     e101_s2[i]=1.0q/(1.0q-2.0q)/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*a0m1_s2[i]-c0m1_s2[i]-(1.0q-2.0q)*(hm+em-dm*(dh+c))*e001_s2[i]);
     e011_s2[i]=1.0q/(1.0q-2.0q)/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*c0m1_s2[i]-(oneph2p2eh)*a0m1_s2[i]+(1.0q-2.0q)*((dh+c)*(hm+em)-(oneph2p2eh)*dm)*e001_s2[i]);
     tp1=oneph2p2eh-(dh+c)*(dh+c);
     e00m1_s2[i] = 1.0q/(-3.0q)*(((hm+em)*(dh+c)-dm*(oneph2p2eh))/tp1*a0m1_s2[i]-yv[i]*a0m1_s2[i] +(dm*(dh+c)-(hm+em))/tp1*c0m1_s2[i]-rv[i]*c0m1_s2[i] -1.0q*(a2+m2)*e001_s2[i] +1.0q*((hm+em)*(hm+em-dm*(dh+c))+dm*(dm*(oneph2p2eh)-(dh+c)*(hm+em)))/tp1*e001_s2[i]);
     tp1=oneph2p2eh-(dh+c)*(dh+c);
     e201_s2[i]=1.0q/(-1.0q)/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*a1m1_s2[i]-rv[i]*c0m1_s2[i]+e00m1_s2[i]-(-1.0q)*(hm+em-dm*(dh+c))*e101_s2[i]);
     e111_s2[i]=1.0q/(1.0q-2.0q)/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*rv[i]*c0m1_s2[i]-(oneph2p2eh)*a1m1_s2[i]-(dh+c)*e00m1_s2[i]+(-1.0q)*((dh+c)*(hm+em)-(oneph2p2eh)*dm)*e101_s2[i]);
     e203_s2[i]=1.0q/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*a11_s2[i]-rv[i]*c01_s2[i]+e001_s2[i]-(hm+em-dm*(dh+c))*e103_s2[i]);
     e113_s2[i]=1.0q/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*yv[i]*a01_s2[i]-c11_s2[i]-(dh+c)*e001_s2[i]-(hm+em-dm*(dh+c))*e013_s2[i]);
     e213_s2[i]=1.0q/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*yv[i]*a11_s2[i]-rv[i]*c11_s2[i]+e011_s2[i]-(dh+c)*e101_s2[i]-(hm+em-dm*(dh+c))*e113_s2[i]);
     e313_s2[i]=1.0q/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*yv[i]*a21_s2[i]-powq(rv[i],2.0q)*c11_s2[i]+2.0q*e111_s2[i]-(dh+c)*e201_s2[i]-(hm+em-dm*(dh+c))*e213_s2[i]);
   }

   for(int i = 0; i < 8; i++) {
     d003[i]= 0.5q*d003_r1[i] +0.5q*d003_r2[i] +0.5q*d003_s1[i] +0.5q*d003_s2[i];
     aas001[i]= a01_s1[i]+a01_s2[i];
     bbr001[i]= b01_r1[i]+b01_r2[i];
     d103[i]=1.0q/(1.0q-e2)*((+e*aas001[i] -bbr001[i]) -(c-de)*yv[i]*d003[i]);
     d013[i]=1.0q/(1.0q-e2)*((+e*bbr001[i] -aas001[i]) -(d-ce)*yv[i]*d003[i]);
     d001[i] =1.0q/((1.0q-e2)*(-1.0q))*( ((a2+yv2[i]) * (1.0q-e2) +2.0q*cd*e*1. * yv2[i] - 1.0q*yv2[i]*(c2+d2)) * d003[i]-(yv[i]*(d-ce) +m*(1.0q-e2)) * a01_s2[i]-(h*(1.0q-e2))*a11_s2[i]-(yv[i]*(d-ce)+sv[i]*(1.0q-e2)) * a01_s1[i]-(yv[i]*(c-de)+(1.0q-e2)*g)*b01_r2[i]-(f*(1.0q-e2))*b11_r2[i]-(yv[i]*(c-e*d)+rv[i]*(1.0q-e2)) * b01_r1[i]);
     aas00m1[i]= a0m1_s1[i]+a0m1_s2[i];
     bbr00m1[i]= b0m1_r1[i]+b0m1_r2[i];
     d101[i]=1.0q/(1.0q-e2)*(1.0q/(1.0q-2.0q)*(+e*aas00m1[i]-bbr00m1[i])-(c-de)*yv[i]*d001[i]);
     d011[i]=1.0q/(1.0q-e2)*(1.0q/(1.0q-2.0q)*(+e*bbr00m1[i]-aas00m1[i])-(d-ce)*yv[i]*d001[i]);
     d00m1[i] =1.0q/((1.0q-e2)*(-3.0q))*( (-(a2+yv2[i]) * (1.0q-e2) +2.0q*cd*e*-1. * yv2[i] +1.0q*yv2[i]*(c2+d2)) * d001[i]-(yv[i]*(d-e*c) +m*(1.0q-e2)) * a0m1_s2[i]-(h*(1.0q-e2))*a1m1_s2[i]-(yv[i]*(d-e*c)+sv[i]*(1.0q-e2)) * a0m1_s1[i]-(yv[i]*(c-e*d)+g*(1.0q-e2)) * b0m1_r2[i]-(f*(1.0q-e2))*b1m1_r2[i]-(yv[i]*(c-e*d)+rv[i]*(1.0q-e2)) * b0m1_r1[i]);
     aas10m1[i]= a1m1_s1[i]+a1m1_s2[i];
     aas01m1[i]= h*a1m1_s2[i] +m* a0m1_s2[i] +s1*a0m1_s1[i];
     bbr01m1[i]= +b1m1_r2[i] +b1m1_r1[i];
     bbr10m1[i]= f*b1m1_r2[i] +g* b0m1_r2[i] +r1*b0m1_r1[i];
     d201[i]= 1.0q/(1.0q-e2)*(1.0q/(-1.0q)*(d00m1[i] +e*aas10m1[i] -bbr10m1[i]) -(c-de)*yv[i]*d101[i]);
     d111_1[i]= 1.0q/(1.0q-e2)*(1.0q/(-1.0q)*(-e*d00m1[i] +e*aas01m1[i] -bbr01m1[i]) -(c-de)*yv[i]*d011[i]);
     d111_2[i]= 1.0q/(1.0q-e2)*(1.0q/(-1.0q)*(-e*d00m1[i] +e*bbr10m1[i] -aas10m1[i]) -(d-ce)*yv[i]*d101[i]);
     d021[i]=1.0q/(1.0q-e2)*(1.0q/(-1.0q)*(d00m1[i]+e*bbr01m1[i]-aas01m1[i])-(d-ce)*yv[i]*d011[i]);
     d111[i]= 0.5q*d111_1[i] +0.5q*d111_2[i];
     aas101[i]= a11_s1[i]+a11_s2[i];
     aas011[i]= h*a11_s2[i] +m* a01_s2[i] +s1*a01_s1[i];
     bbr101[i]= f*b11_r2[i] +g* b01_r2[i] +r1*b01_r1[i];
     bbr011[i]= +b11_r2[i] +b11_r1[i];
     d203[i]=1.0q/(1.0q-e2)*((d001[i]+e*aas101[i]-bbr101[i])-(c-de)*yv[i]*d103[i]);
     d023[i]=1.0q/(1.0q-e2)*((d001[i]+e*bbr011[i]-aas011[i])-(d-ce)*yv[i]*d013[i]);
     d113_1[i]=1.0q/(1.0q-e2)*((-e*d001[i]+e*aas011[i]-bbr011[i])-(c-de)*yv[i]*d013[i]);
     d113_2[i]=1.0q/(1.0q-e2)*((-e*d001[i]+e*bbr101[i]-aas101[i])-(d-ce)*yv[i]*d103[i]);
     d113[i]= 0.5q*d113_1[i] +0.5q*d113_2[i];
     aas111[i] = h*a21_s2[i] +m* a11_s2[i] +s1*a11_s1[i];
     bbr111[i] = f*b21_r2[i] +g* b11_r2[i] +r1*b11_r1[i];
     d213[i]=1.0q/(1.0q-e2)*((d011[i]-e*d101[i]+e*aas111[i]-bbr111[i])-(c-de)*yv[i]*d113[i]);
     d123[i]=1.0q/(1.0q-e2)*((d101[i]-e*d011[i]+e*bbr111[i]-aas111[i])-(d-ce)*yv[i]*d113[i]);
     aas201[i] = a21_s1[i]+a21_s2[i];
     aas021[i] = h2*a21_s2[i] +2.0q*hm*a11_s2[i] +m2* a01_s2[i] +s1*s1*a01_s1[i];
     bbr201[i] = f2*b21_r2[i] +2.0q*fg*b11_r2[i] +g2* b01_r2[i] +r1*r1*b01_r1[i];
     bbr021[i] = +b21_r2[i] +b21_r1[i];
     d303[i]=1.0q/(1.0q-e2)*((2.0q*d101[i]+e*aas201[i]-bbr201[i])-(c-de)*yv[i]*d203[i]);
     d033[i]=1.0q/(1.0q-e2)*((2.0q*d011[i]+e*bbr021[i]-aas021[i])-(d-ce)*yv[i]*d023[i]);
     //
     // TRIPLE INTEGRALS
     //
     /* triple integrals are built from previous __float128 integrals.
        Hijkl = /int/int/int r^i s^j y^k Ra^-l dr ds dy
        Since r2 is a fct(s) and s2 is a fct(r) in the case of a triangle, we introduce
        additional __float128 integrals:
        EESijkl = /int/int r^i s(r)^j y^k Ra^-l dr dy = /int/int r^i (hr+m)^j y^k Ra^-l drdy - /int/int r^i sa^j y^k Ra^-l drdy
        FFRijkl = /int/int r(s)^i s^j y^k Ra^-l dr dy = /int/int (fs+g)^i s^j y^k Ra^-l drdy - /int/int r1^i s^j y^k Ra^-l drdy
     */

     /* SEED FUNCTION */
     /* since no analytical was found for the only tripple seed integral,
        we perform instead an adaptative Simpson quadrature */
     //  functionPtr = &FDr1s1;
     //  Dr1s1 = adaptiveSimpsons(functionPtr, y1, y2, 1e-9, 20q);
     //  functionPtr = &FDr1s2;
     //  Dr1s2 = adaptiveSimpsons(functionPtr, y1, y2, 1e-9, 20q);
     //  functionPtr = &FDr2s1;
     //  Dr2s1 = adaptiveSimpsons(functionPtr, y1, y2, 1e-9, 20q);
     //  functionPtr = &FDr2s2;
     //  Dr2s2 = adaptiveSimpsons(functionPtr, y1, y2, 1e-9, 20q);
     //
     //  H0003[0] = Dr2s2; H0003[1] = 0;
     //  H0003[2] = Dr2s1; H0003[3] = 0;
     //  H0003[4] = Dr1s2; H0003[5] = 0;
     //  H0003[6] = Dr1s1; H0003[7] = 0;
     ffr1001[i] = f*f101_r2[i] +g*f001_r2[i] +r1*f001_r1[i];
     ees0101[i] = h*e101_s2[i] +m*e001_s2[i] +s1*e001_s1[i];
     h0001[i] = 1.0q/(-2.0q)*(a2*0.0q-ffr1001[i]-ees0101[i]-yv[i]*d001[i]);
     ees000m1[i] = e00m1_s1[i] +e00m1_s2[i];
     ffr000m1[i] = f00m1_r2[i] +f00m1_r1[i];
     h0011[i] = (1.0q/-1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -d00m1[i]) +(d-ce)*ees000m1[i] +(c-de)*ffr000m1[i]);
     h1001[i] = 1.0q/-1.0q/(1.0q-e2)*(+e*ees000m1[i] -ffr000m1[i] -1.0q*(de-c)*h0011[i]);
     h0101[i] = 1.0q/-1.0q/(1.0q-e2)*(-ees000m1[i] +e*ffr000m1[i] -1.0q*(ce-d)*h0011[i]);
     ees0001[i] = e001_s1[i] +e001_s2[i];
     ffr0001[i] = f001_r2[i] +f001_r1[i];
     ees1001[i] = e101_s1[i] +e101_s2[i];
     ffr0101[i] = f101_r2[i] +f101_r1[i];
     hijkl[0][i] = (1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -d001[i]) +(d-ce)*ees0001[i] +(c-de)*ffr0001[i]);
     hijkl[1][i] = 1.0q/(1.0q-e2)*( +e*ees0001[i] -ffr0001[i] +(de-c)*hijkl[0][i]);
     hijkl[2][i] = 1.0q/(1.0q-e2)*( -ees0001[i] +e*ffr0001[i] +(ce-d)*hijkl[0][i]);
     hijkl[5][i] = (1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -d101[i]) -(c-de)*h0001[i] +(d-ce)*ees1001[i] +(c-de)*ffr1001[i]);
     hijkl[6][i] = (1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -d011[i])  -(d-ce)*h0001[i] +(d-ce)*ees0101[i] +(c-de)*ffr0101[i]);
     h1103_1[i] = 1.0q/(1.0q-e2)*( -e*h0001[i] +e*ees0101[i] -ffr0101[i] +(de-c)*hijkl[6][i]);
     h1103_2[i] = 1.0q/(1.0q-e2)*( -e*h0001[i] +e*ffr1001[i] -ees1001[i] +(ce-d)*hijkl[5][i]);
     hijkl[7][i]= 0.5q*h1103_1[i] +0.5q*h1103_2[i];
     hijkl[3][i] = 1.0q/(1.0q-e2)*(h0001[i] +e*ees1001[i] -ffr1001[i] +(de-c)*hijkl[5][i]);
     hijkl[4][i] = 1.0q/(1.0q-e2)*(h0001[i] -ees0101[i] +e*ffr0101[i] +(ce-d)*hijkl[6][i]);
     ees2001[i]= e201_s1[i] +e201_s2[i];
     ffr0201[i]= f201_r2[i] +f201_r1[i];
     ffr2001[i]= f2*f201_r2[i] +2.0q*fg*f101_r2[i] +g2*f001_r2[i] +r12*f001_r1[i];
     ees0201[i]= h2*e201_s2[i] +2.0q*hm*e101_s2[i] +m2*e001_s2[i] +s12*e001_s1[i];
     hijkl[8][i] = (1.0q/1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -d201[i]) -2.0q*(c-de)*h1001[i] +(d-ce)*ees2001[i] +(c-de)*ffr2001[i]);
     hijkl[9][i] = (1.0q/1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -d021[i]) -2.0q*(d-ce)*h0101[i] +(d-ce)*ees0201[i] +(c-de)*ffr0201[i]);
     ees1101[i]= h*e201_s2[i] +m*e101_s2[i] +s1*e101_s1[i];
     ffr1101[i]= f*f201_r2[i] +g*f101_r2[i] +r1*f101_r1[i];
     hijkl[12][i] = (1.0q/1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -d111[i]) -(c-de)*h0101[i] -(d-ce)*h1001[i] +(d-ce)*ees1101[i] +(c-de)*ffr1101[i]);
     hijkl[10][i] = 1.0q/(1.0q-e2)*(h0101[i] -e*h1001[i] +e*ees1101[i] -ffr1101[i] +(de-c)*hijkl[12][i]);
     hijkl[11][i] = 1.0q/(1.0q-e2)*(h1001[i] -e*h0101[i] +e*ffr1101[i] -ees1101[i]  +(ce-d)*hijkl[12][i]);
     hijkl[13][i] = 1.0q/(1.0q-e2)*(2.0q*h1001[i] +e*ees2001[i] -ffr2001[i] +(de-c)*hijkl[8][i]);
     hijkl[14][i] = 1.0q/(1.0q-e2)*(2.0q*h0101[i] -ees0201[i] +e*ffr0201[i] +(ce-d)*hijkl[9][i]);
     ees0003[i] = e003_s1[i] +e003_s2[i];
     ffr0003[i] = f003_r2[i] +f003_r1[i];
     hijkl[15][i] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -d003[i]) +(d-ce)*ees0003[i] +(c-de)*ffr0003[i]);
     hijkl[16][i] = 1.0q/3.0q/(1.0q-e2)*( +e*ees0003[i] -ffr0003[i] +3.0q*(de-c)*hijkl[15][i]);
     hijkl[17][i] = 1.0q/3.0q/(1.0q-e2)*( -ees0003[i] +e*ffr0003[i] +3.0q*(ce-d)*hijkl[15][i]);
     ees1003[i] = e103_s1[i] +e103_s2[i];
     ffr1003[i] = f*f103_r2[i] +g*f003_r2[i] +r1*f003_r1[i];
     ees0103[i] = h*e103_s2[i] +m*e003_s2[i] +s1*e003_s1[i];
     ffr0103[i] = f103_r2[i] +f103_r1[i];
     hijkl[20][i] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -d103[i]) -(c-de)*0. +(d-ce)*ees1003[i] +(c-de)*ffr1003[i]);
     hijkl[21][i] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -d013[i]) -(d-ce)*0. +(d-ce)*ees0103[i] +(c-de)*ffr0103[i]);
     h1105_1[i] = 1.0q/3.0q/(1.0q-e2)*( -e*0. +e*ees0103[i] -ffr0103[i] +3.0q*(de-c)*hijkl[21][i]);
     h1105_2[i] = 1.0q/3.0q/(1.0q-e2)*( -e*0. +e*ffr1003[i] -ees1003[i] +3.0q*(ce-d)*hijkl[20][i]);
     hijkl[22][i]= 0.5q*h1105_1[i] +0.5q*h1105_2[i];
     hijkl[18][i] = 1.0q/3.0q/(1.0q-e2)*(0. +e*ees1003[i] -ffr1003[i] +3.0q*(de-c)*hijkl[20][i]);
     hijkl[19][i] = 1.0q/3.0q/(1.0q-e2)*(0. -ees0103[i] +e*ffr0103[i] +3.0q*(ce-d)*hijkl[21][i]);
     ees1103[i] = h*e203_s2[i] +m*e103_s2[i] +s1*e103_s1[i];
     ffr1103[i] = f*f203_r2[i] +g*f103_r2[i] +r1*f103_r1[i];
     hijkl[23][i] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -d113[i]) -(c-de)*hijkl[2][i] -(d-ce)*hijkl[1][i] +(d-ce)*ees1103[i] +(c-de)*ffr1103[i]);
     hijkl[35][i] = 1.0q/3.0q/(1.0q-e2)*(hijkl[2][i] -e*hijkl[1][i] +e*ees1103[i] -ffr1103[i] +3.0q*(de-c)*hijkl[23][i]);
     hijkl[36][i] = 1.0q/3.0q/(1.0q-e2)*(hijkl[1][i] -e*hijkl[2][i] -ees1103[i] +e*ffr1103[i] +3.0q*(ce-d)*hijkl[23][i]);
     ees1013[i] = e113_s1[i] +e113_s2[i];
     ffr1013[i] = f*f113_r2[i] +g*f013_r2[i] +r1*f013_r1[i];
     hijkl[27][i] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(hijkl[1][i] -yv[i]*d103[i]) -(c-de)*hijkl[0][i] +(d-ce)*ees1013[i] +(c-de)*ffr1013[i]);
     hijkl[29][i] = 1.0q/3.0q/(1.0q-e2)*(hijkl[0][i] +e*ees1013[i] -ffr1013[i] +3.0q*(de-c)*hijkl[27][i]);
     ees0113[i] = h*e113_s2[i] +m*e013_s2[i] +s1*e013_s1[i];
     ffr0113[i] = f113_r2[i] +f113_r1[i];
     hijkl[28][i] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(hijkl[2][i] -yv[i]*d013[i]) -(d-ce)*hijkl[0][i] +(d-ce)*ees0113[i] +(c-de)*ffr0113[i]);
     hijkl[30][i] = 1.0q/3.0q/(1.0q-e2)*(hijkl[0][i] -ees0113[i] +e*ffr0113[i] +3.0q*(ce-d)*hijkl[28][i]);
     ees2013[i] = e213_s1[i] +e213_s2[i];
     ffr2013[i] = f2*f213_r2[i] +2.0q*fg*f113_r2[i] +g2*f013_r2[i] +r12*f013_r1[i];
     ees0213[i] = h2*e213_s2[i] +2.0q*hm*e113_s2[i] +m2*e013_s2[i] +s12*e013_s1[i];
     ffr0213[i] = f213_r2[i] +f213_r1[i];
     hijkl[31][i] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(hijkl[3][i] -yv[i]*d203[i]) -2.0q*(c-de)*hijkl[5][i] +(d-ce)*ees2013[i] +(c-de)*ffr2013[i]);
     hijkl[32][i] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(hijkl[4][i] -yv[i]*d023[i]) -2.0q*(d-ce)*hijkl[6][i] +(d-ce)*ees0213[i] +(c-de)*ffr0213[i]);
     ees1113[i] = h*e213_s2[i] +m*e113_s2[i] +s1*e113_s1[i];
     ffr1113[i] = f*f213_r2[i] +g*f113_r2[i] +r1*f113_r1[i];
     hijkl[24][i] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(hijkl[7][i] -yv[i]*d113[i]) -(c-de)*hijkl[6][i] -(d-ce)*hijkl[5][i] +(d-ce)*ees1113[i] +(c-de)*ffr1113[i]);
     hijkl[25][i] = 1.0q/3.0q/(1.0q-e2)*(hijkl[6][i] -e*hijkl[5][i] +e*ees1113[i] -ffr1113[i] +3.0q*(de-c)*hijkl[24][i]);
     hijkl[26][i] = 1.0q/3.0q/(1.0q-e2)*(hijkl[5][i] -e*hijkl[6][i] -ees1113[i] +e*ffr1113[i] +3.0q*(ce-d)*hijkl[24][i]);
     hijkl[33][i] = 1.0q/3.0q/(1.0q-e2)*(2.0q*hijkl[5][i] +e*ees2013[i] -ffr2013[i] +3.0q*(de-c)*hijkl[31][i]);
     hijkl[34][i] = 1.0q/3.0q/(1.0q-e2)*(2.0q*hijkl[6][i] -ees0213[i] +e*ffr0213[i] +3.0q*(ce-d)*hijkl[32][i]);
     ees2003[i] = e203_s1[i] +e203_s2[i];
     ffr0203[i] = f203_r2[i] +f203_r1[i];
     ees0203[i] = h2*e203_s2[i] +2.0q*hm*e103_s2[i] +m2*e003_s2[i] +s12*e003_s1[i];
     ffr2003[i] = f2*f203_r2[i] +2.0q*fg*f103_r2[i] +g2*f003_r2[i] +r12*f003_r1[i];
     hijkl[37][i] = 1.0q/3.0q/(1.0q-e2)*(2.0q*hijkl[1][i] +e*ees2003[i] -ffr2003[i] +3.0q*(de-c)*hijkl[29][i]);
     hijkl[38][i] = 1.0q/3.0q/(1.0q-e2)*(2.0q*hijkl[2][i] -ees0203[i] +e*ffr0203[i] +3.0q*(ce-d)*hijkl[30][i]);
     ees3013[i] = e313_s1[i] +e313_s2[i];
     ffr0313[i] = f313_r2[i] +f313_r1[i];
     ffr3013[i] = f3*f313_r2[i] +3.0q*f2g*f213_r2[i] +3.0q*fg2*f113_r2[i] +g3*f013_r2[i] +r13*f013_r1[i];
     ees0313[i] = h3*e313_s2[i] +3.0q*h2m*e213_s2[i] +3.0q*hm2*e113_s2[i] +m3*e013_s2[i] +s13*e013_s1[i];
     hijkl[44][i] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(hijkl[13][i] -yv[i]*d303[i]) -3.0q*(c-de)*hijkl[8][i] +(d-ce)*ees3013[i] +(c-de)*ffr3013[i]);
     hijkl[45][i] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(hijkl[14][i] -yv[i]*d033[i]) -3.0q*(d-ce)*hijkl[9][i] +(d-ce)*ees0313[i] +(c-de)*ffr0313[i]);
     hijkl[46][i] = 1.0q/3.0q/(1.0q-e2)*(-e*3.0q*hijkl[8][i] -ees3013[i] +e*ffr3013[i] +3.0q*(ce-d)*hijkl[44][i]);
     hijkl[47][i] = 1.0q/3.0q/(1.0q-e2)*(-e*3.0q*hijkl[9][i] +e*ees0313[i] -ffr0313[i] +3.0q*(de-c)*hijkl[45][i]);
     ees2113[i] = h*e313_s2[i] +m*e213_s2[i] +s1*e213_s1[i];
     ffr1213[i] = f*f313_r2[i] +g*f213_r2[i] +r1*f213_r1[i];
     ffr2113[i] = f2*f313_r2[i] +2.0q*fg*f213_r2[i] +g2*f113_r2[i] +r12*f113_r1[i];
     ees1213[i] = h2*e313_s2[i] +2.0q*hm*e213_s2[i] +m2*e113_s2[i] +s12*e113_s1[i];
     hijkl[40][i] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(hijkl[10][i] -yv[i]*d213[i]) -2.0q*(c-de)*hijkl[12][i] -(d-ce)*hijkl[8][i] +(d-ce)*ees2113[i] +(c-de)*ffr2113[i]);
     hijkl[41][i] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(hijkl[11][i] -yv[i]*d123[i]) -(c-de)*hijkl[9][i] -2.0q*(d-ce)*hijkl[12][i] +(d-ce)*ees1213[i] +(c-de)*ffr1213[i]);
     h2215_1[i] = 1.0q/3.0q/(1.0q-e2)*(hijkl[9][i] -e*2.0q*hijkl[12][i] +e*ees1213[i] -ffr1213[i] +3.0q*(de-c)*hijkl[41][i]);
     h2215_2[i] = 1.0q/3.0q/(1.0q-e2)*(hijkl[8][i] -e*2.0q*hijkl[12][i] +e*ffr2113[i] -ees2113[i] +3.0q*(ce-d)*hijkl[40][i]);
     hijkl[39][i]= 0.5q*h2215_1[i] +0.5q*h2215_2[i];
     hijkl[42][i] = 1.0q/3.0q/(1.0q-e2)*(3.0q*hijkl[8][i] +e*ees3013[i] -ffr3013[i] +3.0q*(de-c)*hijkl[44][i]);
     hijkl[43][i] = 1.0q/3.0q/(1.0q-e2)*(3.0q*hijkl[9][i] -ees0313[i] +e*ffr0313[i] +3.0q*(ce-d)*hijkl[45][i]);
   }

   // Evaluating the integrals.
   for (int i = 0; i < i_num_integrals; i++){
     o_sch[i] = dot_product(hijkl[i], signv, n_limits);
   }
  return o_sch;
}

__float128 *node_force_quad_triangle(__float128 *i_sch, __float128 i_vec_int[][3],
                                 __float128 i_a, __float128 i_b, __float128 i_c,
                                 __float128 i_d, __float128 i_e, __float128 i_f,
                                 __float128 i_factor,
                                 __float128 *o_force){
  // Calculates the force on a node.
  __float128 f[11];
  f [0] = i_a*i_sch [8] + i_b*i_sch [5] + i_c*i_sch [9] + i_d*i_sch [6] + i_e*i_sch[12] + i_f*i_sch [0];
  f [1] = i_a*i_sch[10] + i_b*i_sch [7] + i_c*i_sch[14] + i_d*i_sch [4] + i_e*i_sch[11] + i_f*i_sch [2];
  f [2] = i_a*i_sch[13] + i_b*i_sch [3] + i_c*i_sch[11] + i_d*i_sch [7] + i_e*i_sch[10] + i_f*i_sch [1];
  f [3] = i_a*i_sch[29] + i_b*i_sch[20] + i_c*i_sch[30] + i_d*i_sch[21] + i_e*i_sch[23] + i_f*i_sch[15];
  f [4] = i_a*i_sch[35] + i_b*i_sch[22] + i_c*i_sch[38] + i_d*i_sch[19] + i_e*i_sch[36] + i_f*i_sch[17];
  f [5] = i_a*i_sch[37] + i_b*i_sch[18] + i_c*i_sch[36] + i_d*i_sch[22] + i_e*i_sch[35] + i_f*i_sch[16];
  f [6] = i_a*i_sch[39] + i_b*i_sch[26] + i_c*i_sch[43] + i_d*i_sch[34] + i_e*i_sch[47] + i_f*i_sch[30];
  f [7] = i_a*i_sch[42] + i_b*i_sch[33] + i_c*i_sch[39] + i_d*i_sch[25] + i_e*i_sch[46] + i_f*i_sch[29];
  f [8] = i_a*i_sch[40] + i_b*i_sch[24] + i_c*i_sch[45] + i_d*i_sch[32] + i_e*i_sch[41] + i_f*i_sch[28];
  f [9] = i_a*i_sch[44] + i_b*i_sch[31] + i_c*i_sch[41] + i_d*i_sch[24] + i_e*i_sch[40] + i_f*i_sch[27];
  f[10] = i_a*i_sch[46] + i_b*i_sch[25] + i_c*i_sch[47] + i_d*i_sch[26] + i_e*i_sch[39] + i_f*i_sch[23];
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

void compute_forces_quad_triangle(__float128 *i_sch, __float128 i_vec_int[][3],
                                     __float128 i_r1, __float128 i_r2,
                                     __float128 i_s1, __float128 i_s2,
                                     __float128 i_one_o_dr_sq, __float128 i_one_o_ds_sq,
                                     __float128 i_one_o_dr, __float128 i_one_o_ds,
                                     __float128 i_one_o_dr_one_o_ds,
                                     __float128 i_factor,
                                     __float128 *o_nodal_force[n_nodes], __float128 *o_total_force){
  // Calculating nodal forces
  __float128 a, b, c, d, e, f;
  // x3
  a = 2.0q*i_one_o_dr_sq;
  b = -4.0q*i_r1*i_one_o_dr_sq-3.0q*i_one_o_dr-4.0q*i_s1*i_one_o_dr_one_o_ds;
  c = 2.0q*i_one_o_ds_sq;
  d = -4.0q*i_s1*i_one_o_ds_sq-3.0q*i_one_o_ds-4.0q*i_r1*i_one_o_dr_one_o_ds;
  e = 4.0q*i_one_o_dr_one_o_ds;
  f = 1.0q+4.0q*i_r1*i_s1*i_one_o_dr_one_o_ds+3.0q*i_r1*i_one_o_dr+3.0q*i_s1*i_one_o_ds+2.0q*i_r1*i_r1*i_one_o_dr_sq+2.0q*i_s1*i_s1*i_one_o_ds_sq;
  node_force_quad_triangle(i_sch, i_vec_int, a, b, c, d, e, f, i_factor, o_nodal_force[0]);
  // x4
  a = 2.0q;
  b = -3.0q*i_r1 - i_r2;
  //c = d = e = 0.0q;
  f = i_r1*i_r1 + i_r1*i_r2;
  node_force_quad_triangle(i_sch, i_vec_int, a, b, 0., 0., 0., f, i_factor*i_one_o_dr_sq, o_nodal_force[1]);
  // x5
  //a = b = e = 0.0q;
  b = 0.0q;
  c = 2.0q;
  d = -3.0q*i_s1 - i_s2;
  f = i_s1*i_s1 + i_s1*i_s2;
  node_force_quad_triangle(i_sch, i_vec_int, 0., 0., c, d, 0., f, i_factor*i_one_o_ds_sq, o_nodal_force[2]);
  // x6
  a = -4.0q*i_one_o_dr_sq;
  b = 4.0q*i_one_o_dr+8.0q*i_r1*i_one_o_dr_sq+4.0q*i_s1*i_one_o_dr_one_o_ds;
  //c = 0.0q;
  d = 4.0q*i_r1*i_one_o_dr_one_o_ds;
  e = -4.0q*i_one_o_dr_one_o_ds;
  f = -4.0q*i_r1*i_one_o_dr-4.0q*i_r1*i_r1*i_one_o_dr_sq-4.0q*i_r1*i_s1*i_one_o_dr_one_o_ds;
  node_force_quad_triangle(i_sch, i_vec_int, a, b, 0., d, e, f, i_factor, o_nodal_force[3]);
  // x7
  //a = 0.0q;
  b = 4.0q*i_s1*i_one_o_dr_one_o_ds;
  c = -4.0q*i_one_o_ds_sq;
  d = 4.0q*i_one_o_ds+8.0q*i_s1*i_one_o_ds_sq+4.0q*i_r1*i_one_o_dr_one_o_ds;
  e = -4.0q*i_one_o_dr_one_o_ds;
  f = -4.0q*i_s1*i_one_o_ds-4.0q*i_s1*i_s1*i_one_o_ds_sq-4.0q*i_r1*i_s1*i_one_o_dr_one_o_ds;
  node_force_quad_triangle(i_sch, i_vec_int, 0., b, c, d, e, f, i_factor, o_nodal_force[4]);
  // x8
  // a = c = 0.0q;
  b = -4.0q*i_s1;
  d = -4.0q*i_r1;
  e = 4.0q;
  f = 4.0q*i_r1*i_s1;
  node_force_quad_triangle(i_sch, i_vec_int, 0., b, 0., d, e, f, i_factor*i_one_o_dr_one_o_ds, o_nodal_force[5]);
  //[8] [5] [9] [6] [12] [0]

  for (int i = 0; i < n_nodes; i++){
    o_total_force[0] += o_nodal_force[i][0];
    o_total_force[1] += o_nodal_force[i][1];
    o_total_force[2] += o_nodal_force[i][2];
  }
}

 void nodal_surface_force_quad_triangle(__float128 *x1, __float128 *x2, __float128 *x3, __float128 *x4, __float128 *x5, __float128 *b, __float128 *p, __float128 *q, __float128 *n, __float128 p_dot_q, __float128 p_dot_q_sq, __float128 mu, __float128 nu, __float128 a, __float128 a_sq, __float128 one_m_nu, __float128 factor, __float128 *nodal_force[n_nodes], __float128 *total_force){
   // Characteristic vectors.
   __float128 t[3];
   // Basis vectors (unitary).
   __float128 p_x_t[3], q_x_t[3];
   //  Auxiliary constants to reduce computational cost.
   //  p dot t, q dot t
   __float128 p_dot_t, q_dot_t;
   // Limits of the distance vector from the plane to dislocation line segment.
   // r_lim[0][] = vector from x1 to x3, r_lim[1][] = vector from x1 to x4, r_lim[2][] = vector from x1 to x5, r_lim[3][] = vector from x2 to x3.
   __float128 r_lim[4][3];
   // ds = (sp[2] - sp[1]), rs = (rp[2] - rp[1]), one_o_ds = 1/(sp[2] - sp[1]), one_o_rs = 1/(rp[2] - rp[1])
   __float128 ds, dr, one_o_dr, one_o_ds, one_o_dr_sq, one_o_ds_sq, one_o_dr_one_o_ds, dr_o_ds, ds_o_dr;
   //r, s limits, the points where r and s become functions of each other are r[2] and s[2].
   __float128 rp[3], sp[3];
   //  y, r, s coordinates.
   __float128 y[n_limits], r[n_limits], s[n_limits];
   // r(s), s(r) coordinates.
   __float128 rfs[n_limits], sfr[n_limits];
   // Vectors for the integrals.
   __float128 vec_int[11][3];
   // Scalar value of the integrals.
   __float128 sch[48];
   // Set forces to zero.
   init_force(nodal_force, total_force);
   // Build unit vector t.
   init_vector (x1, x2, 3, t);
   // Dot products.
   p_dot_t = dot_product(p, t, 3);
   q_dot_t = dot_product(q, t, 3);
   //*******************WARNING*******************//
   // This formulation assumes x3-x6 is diagonal! //
   //*******************WARNING*******************//
   cross_product(p, t, p_x_t);
   cross_product(q, t, q_x_t);
   // Vectors between x3 and x1, x4 and x1, x5 and x1, x3 and x2.
   for (int i = 0; i < 3; i++){
     r_lim[0][i] = x3[i] - x1[i];
     r_lim[1][i] = x4[i] - x1[i];
     r_lim[2][i] = x5[i] - x1[i];
     r_lim[3][i] = x3[i] - x2[i];
   }
   // Integral bounds for y, r, s.
   // Assign coordinates for the evaluation of the integrals.
   /*
    *  Integrals are functions of y,r and s. But ultimately, the nodal force evaluation
    *  will be the sum of evaluations of the antiderivative for the various bounds.
    *  It is therefore more convenient to organize integrals and variables in the shape of
    *  vectors of 8 components. The different components correspond to permutation
    * of the 2 bounds per variables: rp[0], r2, sp[0], s2, y1 and y2.
    */
  rp[0] = init_point(r_lim[0], q_x_t, p, q_x_t, 3);
  rp[1] = init_point(r_lim[1], q_x_t, p, q_x_t, 3);
  sp[0] = init_point(r_lim[0], p_x_t, q, p_x_t, 3);
  sp[1] = init_point(r_lim[2], p_x_t, q, p_x_t, 3);
  y[0] = y[2] = y[4] = y[6] = init_point(r_lim[3], n, t, n    , 3);
  y[1] = y[3] = y[5] = y[7] = init_point(r_lim[0], n, t, n    , 3);
  r[0] = r[1] = r[2] = r[3] = rp[1];
  r[4] = r[5] = r[6] = r[7] = rp[0];
  s[0] = s[1] = s[4] = s[5] = sp[1];
  s[2] = s[3] = s[6] = s[7] = sp[0];
  // in the case of triangular element and depending on the integration order,
  // the upper bound of r can be a function of s. These vectors exprese this dependence.
         ds = sp[1] - sp[0];
         dr = rp[1] - rp[0];
  one_o_dr = 1.0q/dr;
  one_o_ds = 1.0q/ds;
  one_o_dr_sq = one_o_dr * one_o_dr;
  one_o_ds_sq = one_o_ds * one_o_ds;
  one_o_dr_one_o_ds = one_o_dr * one_o_ds;
  dr_o_ds = dr * one_o_ds;
  ds_o_dr = ds * one_o_dr;
  sp [2] = sp[1] + rp[0] * ds_o_dr;
  rp [2] = rp[1] + sp[0] * dr_o_ds;
  rfs[0] = rfs[1] = rp[2] - dr_o_ds * sp[1];
  rfs[2] = rfs[3] = rp[2] - dr_o_ds * sp[0];
  rfs[4] = rfs[5] = rfs[6] = rfs[7] = rp[0];
  sfr[0] = sfr[1] = sp[2] - ds_o_dr * rp[1];
  sfr[2] = sfr[3] = sfr[6] = sfr[7] = sp[0];
  sfr[4] = sfr[5] = sp[2] - ds_o_dr * rp[0];

  // Calculate vectors for integrals.
  integral_vector(p, q, b, t, n, one_m_nu, a_sq, vec_int);
  // Calculate integrals.
  integrals_quad_triangle(r, s, y, rfs, sfr, rp[0], sp[0], rp[2], sp[2], -dr_o_ds, -ds_o_dr, p_dot_t, q_dot_t, p_dot_q, a_sq, sch, 48);
  // Calculate nodal forces.
  compute_forces_quad_triangle(sch, vec_int, rp[0], rp[1], sp[0], sp[1], one_o_dr_sq, one_o_ds_sq, one_o_dr, one_o_ds, one_o_dr_one_o_ds, factor, nodal_force, total_force);
 }

 void main_nodal_surface_force_quad_triangle(__float128 *x1, __float128 *x2, __float128 *x3, __float128 *x4, __float128 *x5, __float128 *b, __float128 mu, __float128 nu, __float128 a, __float128 *nodal_force[n_nodes], __float128 *total_force){
   /*
     Forces
     nodal_force[0][] = F_x3[x, y, z], nodal_force[1][] = F_x4[x, y, z],
     nodal_force[2][] = F_x5[x, y, z], nodal_force[3][] = F_x6[x, y, z]
     total_force[x, y, z] = F_x3[x, y, z] + F_x4[x, y, z] + F_x5[x, y, z] + F_x6[x, y, z]
   */
   // Characteristic vectors.
   __float128 p[3], q[3], n[3], t[3];
   // t dot q, p dot q, (p dot q)^2, alpha = sine(angle between p and q).
   __float128 t_dot_n, p_dot_q, p_dot_q_sq, alpha;
   __float128 a_sq, factor, one_m_nu;
   __float128 rot_centre[3], rot_x1[3], rot_x2[3];
   __float128 t_x_n[3], mag_t_x_n, p_total_force[3], *p_nodal_force[n_nodes];
   int rotation;
   rotation = 3;
   // Vectors
   init_vector(x1, x2, 3, t);
   init_vector(x3, x4, 3, p);
   init_vector(x3, x5, 3, q);
   cross_product(p, q, n);
   normalise_vector(n, 3, n);
   t_dot_n = dot_product(t, n, 3);
   p_dot_q = dot_product(p, q, 3);
   p_dot_q_sq = p_dot_q * p_dot_q;
   // Auxiliary variables.
   a_sq     = a*a;
   one_m_nu = 1.0q-nu;
   alpha    = sinq(acosq(p_dot_q));
   factor   = 0.25q*alpha*mu/M_PI/one_m_nu;
   if(t_dot_n != 0.0q){
     nodal_surface_force_quad_triangle(x1, x2, x3, x4, x5, b, p, q, n, p_dot_q, p_dot_q_sq, mu, nu, a, a_sq, one_m_nu, factor, nodal_force, total_force);
   }
   else{
     for (int i = 0; i < n_nodes; i++){
       p_nodal_force[i] = malloc(3*sizeof(__float128));
     }
     __float128 angle = 0.01q*M_PI/180.0q;
     int j;
     cross_product(t, n, t_x_n);
     mag_t_x_n = sqrtq(dot_product(t_x_n, t_x_n, 3));
     for (int i = 0; i < 3; i++){
       // Halfway between x1 and x2. x1 + (x2-x1)/2
       rot_centre[i] = 0.5q*(x1[i] + x2[i]);
       t_x_n[i] = t_x_n[i]/mag_t_x_n;
     }
     //FILE *fp;
     //fp = fopen(".0q/tests/test2.txt", "w");
     for (int i = 0; i < rotation; i++){
       j = i+1;
       arbitrary_rotation_matrix_3d(j*angle, rot_centre, t_x_n, x1, rot_x1);
       arbitrary_rotation_matrix_3d(j*angle, rot_centre, t_x_n, x2, rot_x2);
       nodal_surface_force_quad_triangle(rot_x1, rot_x2, x3, x4, x5, b, p, q, n, p_dot_q, p_dot_q_sq, mu, nu, a, a_sq, one_m_nu, factor, p_nodal_force, p_total_force);
       add_force(p_nodal_force, p_total_force, nodal_force, total_force);
       arbitrary_rotation_matrix_3d(-j*angle, rot_centre, t_x_n, x1, rot_x1);
       arbitrary_rotation_matrix_3d(-j*angle, rot_centre, t_x_n, x2, rot_x2);
       nodal_surface_force_quad_triangle(rot_x1, rot_x2, x3, x4, x5, b, p, q, n, p_dot_q, p_dot_q_sq, mu, nu, a, a_sq, one_m_nu, factor, p_nodal_force, p_total_force);
       add_force(p_nodal_force, p_total_force, nodal_force, total_force);
     }
     //fclose(fp);
     mean_force(nodal_force, total_force, rotation*2);
     for (int i = 0; i < n_nodes; i++){
       free(p_nodal_force[i]);
     }
   }
 }

 // Testing purposes
 int main(void){
   __float128 x1[3], x2[3], x3[3], x4[3], x5[3], b[3], mu, nu, a, *nodal_force[n_nodes], total_force[3];
   char x11[256], x12[256], x13[256]; // String size
   char x21[256], x22[256], x23[256];
   char x31[256], x32[256], x33[256];
   char x41[256], x42[256], x43[256];
   char x51[256], x52[256], x53[256];
   char b1[256], b2[256], b3[256];
   char c_mu[256], c_nu[256], c_a[256];
   for(int i = 0; i < n_nodes; i++){
     nodal_force[i] = malloc(3*sizeof(__float128));
   }
   FILE * ptr_file;
   ptr_file =fopen("input.txt", "r");
   if (ptr_file == NULL){
     printf("File does not exist.\n");
   }
   fscanf(ptr_file, "%s %s %s \n", x11, x12, x13 );
   fscanf(ptr_file, "%s %s %s \n", x21, x22, x23 );
   fscanf(ptr_file, "%s %s %s \n", x31, x32, x33 );
   fscanf(ptr_file, "%s %s %s \n", x41, x42, x43 );
   fscanf(ptr_file, "%s %s %s \n", x51, x52, x53 );
   fscanf(ptr_file, "%s %s %s \n", b1, b2, b3 );
   fscanf(ptr_file, "%s \n", c_mu );
   fscanf(ptr_file, "%s \n", c_nu );
   fscanf(ptr_file, "%s \n", c_a );
   fclose(ptr_file);
   // Converstion Char into quad prec floating point numbers
   x1[0] = strtoflt128(x11, NULL); x1[1] = strtoflt128(x12, NULL); x1[2] = strtoflt128(x13, NULL);
   x2[0] = strtoflt128(x21, NULL); x2[1] = strtoflt128(x22, NULL); x2[2] = strtoflt128(x23, NULL);
   x3[0] = strtoflt128(x31, NULL); x3[1] = strtoflt128(x32, NULL); x3[2] = strtoflt128(x33, NULL);
   x4[0] = strtoflt128(x41, NULL); x4[1] = strtoflt128(x42, NULL); x4[2] = strtoflt128(x43, NULL);
   x5[0] = strtoflt128(x51, NULL); x5[1] = strtoflt128(x52, NULL); x5[2] = strtoflt128(x53, NULL);
   b [0] = strtoflt128(b1, NULL) ; b [1] = strtoflt128(b2, NULL) ; b [2] = strtoflt128(b3 , NULL);
   mu = strtoflt128(c_mu, NULL);
   nu = strtoflt128(c_nu, NULL);
   a  = strtoflt128(c_a , NULL);
   // Nodal force calculation.
   main_nodal_surface_force_quad_triangle(x1, x2, x3, x4, x5, b, mu, nu, a, nodal_force, total_force);
   for(int i=0;i<3;i++) { printf("%Qf\t", nodal_force[0][i]*10000.q); }; printf("\n");
   for(int i=0;i<3;i++) { printf("%Qf\t", nodal_force[1][i]*10000.q); }; printf("\n");
   for(int i=0;i<3;i++) { printf("%Qf\t", nodal_force[2][i]*10000.q); }; printf("\n");
   for(int i=0;i<3;i++) { printf("%Qf\t", nodal_force[3][i]*10000.q); }; printf("\n");
   for(int i=0;i<3;i++) { printf("%Qf\t", nodal_force[4][i]*10000.q); }; printf("\n");
   for(int i=0;i<3;i++) { printf("%Qf\t", nodal_force[5][i]*10000.q); }; printf("\n");
   for(int i = 0; i < n_nodes; i++){
     free(nodal_force[i]);
   }
   return 0;
 }
