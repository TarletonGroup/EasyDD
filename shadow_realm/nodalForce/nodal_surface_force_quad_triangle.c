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
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
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
const double pi = 4.0 * atan(1.0);
 /*
  Auxiliary functions.
 */
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
   double *o_vec = malloc(3*sizeof(double));
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
     fprintf(stderr, "ERROR: nodal_surface_force_quad_triangle: normalise_vector: mag_vec = 0: A vector cannot have magnitude 0\n");
     exit(EXIT_FAILURE);
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
     fprintf(stderr, "ERROR: nodal_surface_force_quad_triangle: normalise_vector2: o_mag_vec = 0: A vector cannot have magnitude 0\n");
     exit(EXIT_FAILURE);
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
/*
 double *arbitrary_rotation_matrix_3d(double theta, double *abc, double *uvw, double *xyz, double *o_xpypzp){
   // Rotates point xyz an angle of theta about the unit vector uvw passing through the point abc.
   double costheta       = cos(theta);
   double one_m_costheta = 1. - costheta;
   double mag_uvw = dot_product(uvw, uvw, 3);
   if(mag_uvw != 1.){
     normalise_vector(uvw, 3, uvw);
   }
   o_xpypzp[0] = (abc[0]*(uvw[1]*uvw[1] + uvw[2]*uvw[2]) - uvw[0]*(abc[1]*uvw[1] + abc[2]*uvw[2] - uvw[0]*xyz[0] - uvw[1]*xyz[1] - uvw[2]*xyz[2]))*one_m_costheta + xyz[0]*costheta + (-abc[2]*uvw[1] + abc[1]*uvw[2] - uvw[2]*xyz[1] + uvw[1]*xyz[2])*sin(theta);
   o_xpypzp[1] = (abc[1]*(uvw[0]*uvw[0] + uvw[2]*uvw[2]) - uvw[1]*(abc[0]*uvw[0] + abc[2]*uvw[2] - uvw[0]*xyz[0] - uvw[1]*xyz[1] - uvw[2]*xyz[2]))*one_m_costheta + xyz[1]*costheta + ( abc[2]*uvw[0] - abc[0]*uvw[2] + uvw[2]*xyz[0] - uvw[0]*xyz[2])*sin(theta);
   o_xpypzp[2] = (abc[2]*(uvw[0]*uvw[0] + uvw[1]*uvw[1]) - uvw[2]*(abc[0]*uvw[0] + abc[1]*uvw[1] - uvw[0]*xyz[0] - uvw[1]*xyz[1] - uvw[2]*xyz[2]))*one_m_costheta + xyz[2]*costheta + (-abc[1]*uvw[0] + abc[0]*uvw[1] - uvw[1]*xyz[0] + uvw[0]*xyz[1])*sin(theta);
   //printf("%f, %f, %f\n", xpypzp[0], xpypzp[1], xpypzp[2]);
   return o_xpypzp;
 }*/

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

 double init_point(double *i_vec1, double *i_vec2,
                   double *i_vec3, double *i_vec4,
                   int     i_vec_size){
   // Initialises points on the surface element given four vectors.
   double result = 0.0;
   double denom = 0.0;
   denom = dot_product(i_vec3, i_vec4, i_vec_size);
   if (denom == 0.0){
     fprintf(stderr, "nodal_surface_force_quad_triangle: init_point: division by zero, dot_product(i_vec3, i_vec4) = %2.14f\n", denom);
     exit(EXIT_FAILURE);
   }
   result = dot_product(i_vec1, i_vec2, i_vec_size)/denom;
   return result;
 }

 double seed_single_integral(double a, double b){
   // Single integral seed.
   double result = 0.0;
   double arg;
   arg = a+b;
   if (arg <= 0.0){
     fprintf(stderr, "nodal_surface_force_quad_triangle: seed_single_integral: log(%2.14f) = %2.14f\n", arg, log(arg));
     exit(EXIT_FAILURE);
   }
   result = log(arg);
   return result;
 }

 void init_force(double *nodal_force[n_nodes], double *total_force){
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

 void add_force(double *p_nodal_force[n_nodes], double *p_total_force, double *nodal_force[n_nodes], double *total_force){
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

 void mean_force(double *nodal_force[n_nodes], double *total_force, int n_samples){
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

 void integral_vector(double *i_p, double *i_q, double *i_b, double *i_t, double *i_n,
                      double i_one_m_nu, double i_a_sq,
                      double o_vec_int[][3]){
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
  double t_x_b[3], p_x_b[3],
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
  double t_dot_n,
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

 double *integrals_quad_triangle(double *rv , double *sv , double *yv,
                                 double *rfs, double *sfr,
                                 double r1, double s1,
                                 double g, double m,
                                 double f, double h,
                                 double c, double d, double e,
                                 double a2,
                                 double *o_sch, int i_num_integrals){
   // Sign vector for quick evaluation of integrals via the dot product.
   static double signv[8] = {1.0,-1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0};
   double rv2[8], sv2[8], yv2[8];
   double ra[8], rra[8], sra[8], rar1[8], rar2[8], ras1[8], ras2[8];
   double rdotp[8], rdotq[8], rrdott[8], srdott[8];
   double rardotp[8], rardotq[8], rrardott[8], srardott[8];
   // integrals
   double a01_s1[8], a01_s2[8], b01_r1[8], b01_r2[8];
   double c01_s1[8], c01_s2[8], c01_r1[8], c01_r2[8];
   double a11_s1[8], a11_s2[8], b11_r1[8], b11_r2[8];
   double c11_s1[8], c11_s2[8], c11_r1[8], c11_r2[8];
   double a0m1_s1[8], a0m1_s2[8], b0m1_r1[8], b0m1_r2[8];
   double c0m1_s1[8], c0m1_s2[8], c0m1_r1[8], c0m1_r2[8];
   double a1m1_s1[8], a1m1_s2[8], b1m1_r1[8], b1m1_r2[8];
   double c1m1_s1[8], c1m1_s2[8], c1m1_r1[8], c1m1_r2[8];
   double aa, bb, cc;
   double vtmp1[8], vtmp2[8];
   double a0m3_s1[8], a0m3_s2[8], b0m3_r1[8], b0m3_r2[8];
   double c0m3_s1[8], c0m3_s2[8], c0m3_r1[8], c0m3_r2[8];
   double a1m3_s1[8], a1m3_s2[8], b1m3_r1[8], b1m3_r2[8];
   double a2m1_s1[8], a2m1_s2[8], b2m1_r1[8], b2m1_r2[8];
   double a21_s1[8], a21_s2[8], b21_r1[8], b21_r2[8];
   double a31_s1[8], a31_s2[8], b31_r1[8], b31_r2[8];
   double root, tp1, tp2, tp3, tp4, tp5;
   double e003_s1[8], e003_s2[8], f003_r1[8], f003_r2[8];
   double d003[8], d003_r1[8], d003_r2[8], d003_s1[8], d003_s2[8];
   double e103_s1[8], e103_s2[8], f103_r1[8], f103_r2[8];
   double d103[8],d013[8];
   double aas001[8], bbr001[8];
   double e001_s1[8], e001_s2[8], e101_s1[8], e101_s2[8], e011_s1[8], e011_s2[8];
   double f001_r1[8], f001_r2[8], f101_r1[8], f101_r2[8], f011_r1[8], f011_r2[8];
   double d001[8], d011[8], d101[8];
   double aas00m1[8], bbr00m1[8];
   double e00m1_s1[8], e00m1_s2[8], f00m1_r1[8], f00m1_r2[8];
   double d00m1[8], d00m3[8];
   double e00m3_s1[8], e00m3_s2[8], f00m3_r1[8], f00m3_r2[8];
   double f10m1_r1[8], f10m1_r2[8], e10m1_s1[8], e10m1_s2[8];
   double f201_r1[8], f201_r2[8], e201_s1[8], e201_s2[8];
   double f111_r1[8], f111_r2[8], e111_s1[8], e111_s2[8];
   double aas00m3[8], bbr00m3[8];
   double d10m1[8], d01m1[8];
   double aas10m1[8], aas01m1[8], aas11m1[8];
   double bbr01m1[8], bbr10m1[8], bbr11m1[8];
   double d201[8], d021[8], d111[8], d111_1[8], d111_2[8];
   double e203_s1[8], e203_s2[8], f203_r1[8], f203_r2[8];
   double aas101[8], aas011[8];
   double bbr101[8], bbr011[8];
   double d203[8], d023[8], d113[8], d113_1[8], d113_2[8];
   double e013_s1[8], e013_s2[8], f013_r1[8], f013_r2[8];
   double e023_s1[8], e023_s2[8], f023_r1[8], f023_r2[8];
   double e113_s1[8], e113_s2[8], f113_r1[8], f113_r2[8];
   double e123_s1[8], e123_s2[8], f123_r1[8], f123_r2[8];
   double e213_s1[8], e213_s2[8], f213_r1[8], f213_r2[8];
   double aas111[8], bbr111[8], d213[8], d123[8];
   double f313_r1[8], f223_r1[8], f313_r2[8], f223_r2[8];
   double e313_s1[8],   e223_s1[8], e313_s2[8], e223_s2[8];
   double aas211[8], aas121[8], bbr211[8], bbr121[8];
   double aas201[8], aas021[8] ,bbr201[8], bbr021[8] ;
   double d313[8], d223[8], d223_1[8], d223_2[8], d133[8], d303[8], d033[8];
   double h0001[8], ffr1001[8], ees0101[8];
   double h000m1[8], ffr100m1[8], ees010m1[8];
   double ees000m1[8], ffr000m1[8], h0011[8], h1001[8], h0101[8];
   double ees100m1[8], ffr010m1[8], h1011[8], h0111[8], h1101[8], h2001[8], h0201[8];
   double h1103_1[8], h1103_2[8];
   double ffr0001[8], ees0001[8], ffr0101[8], ees1001[8];
   double h0023[8];
   double h0123[8], ees0111[8], ffr0111[8], h1023[8], ffr1011[8], ees1011[8];
   double ees0011[8], ffr0011[8], ees2001[8], ffr0201[8], ffr2001[8], ees0201[8], ees1101[8], ffr1101[8];
   double h1105_1[8], h1105_2[8], h0025[8];
   double ees0003[8], ffr0003[8], ees1003[8], ffr1003[8], ees0103[8], ffr0103[8], ees0013[8], ffr0013[8];
   double ees1103[8], ffr1103[8], ees1013[8], ffr1013[8], ees0113[8];
   double ffr0113[8], ees2013[8], ffr2013[8], ees0213[8], ffr0213[8];
   double ees1113[8], ffr1113[8];
   double h1135[8], h2215_1[8], h2215_2[8];
   double ees2003[8], ffr0203[8], ees0203[8], ffr2003[8], ees3013[8], ffr0313[8], ffr3013[8], ees0313[8];
   double ees1123[8], ffr1123[8], ees2113[8], ffr1213[8], ffr2113[8], ees1213[8];
   // Triple integrals.
   double hijkl[i_num_integrals][n_limits];
   // Indices for the evaluation of integrals.
   int ir1[4],ir2[4],is1[4],is2[4];
   // Auxiliary scalars.
   double c2, d2, f2, g2, e2;
   double ce, de, ef, fg, eg, cf, cg, onepf2p2ef, oneph2p2eh;
   double h2, m2;
   double cd, hm, dh, eh, em, dm;
   double f2g, fg2, h2m, hm2, f3, h3, g3, m3, cde;
   double r12, r13, s12, s13;
   // Auxiliary scalars
   c2 = c*c; d2 = d*d; e2 = e*e; h2=h*h; f2=f*f; g2=g*g; m2=m*m; eh=e*h; ef=e*f;
   dh=h*d; cf=f*c; cd=c*d; cg=c*g; fg=f*g; eg=e*g; em=e*m; hm=h*m; dm=d*m;
   de=d*e; ce=c*e; f3=f2*f; h3=h2*h; m3=m2*m; g3=g2*g; h2m=h2*m; hm2=h*m2;
   f2g=f2*g; fg2=f*g2; cde= cd*e;
   r12 = r1*r1; r13 = r12*r1; s12 = s1*s1; s13 = s12*s1;
   onepf2p2ef = 1.+f2+2.*ef; oneph2p2eh = 1.+h2+2.*eh;
   for(int i = 0; i < n_limits; i++){
     rv2[i] = rv[i]*rv[i];
     sv2[i] = sv[i]*sv[i];
     yv2[i] = yv[i]*yv[i];
     // ra = sqrt(r.r)
     ra      [i] = sqrt(yv2[i]+rv2[i]+sv2[i]+2.*c*rv[i]*yv[i]+2.*d*sv[i]*yv[i]+2.*e*rv[i]*sv[i]+a2);
     rra     [i] = sqrt(yv2[i]+rv2[i]+sfr[i]*sfr[i]+2.*c*rv[i]*yv[i]+2.*d*sfr[i]*yv[i]+2.*e*rv[i]*sfr[i]+a2);
     sra     [i] = sqrt(yv2[i]+rfs[i]*rfs[i]+sv2[i]+2.*c*rfs[i]*yv[i]+2.*d*sv[i]*yv[i]+2.*e*rfs[i]*sv[i]+a2);
     rar1    [i] = sqrt(yv2[i]+r1*r1+sv2[i]+2.*c*r1*yv[i]+2.*d*sv[i]*yv[i]+2.*e*r1*sv[i]+a2);
     rar2    [i] = sqrt(yv2[i]+(onepf2p2ef)*sv2[i]+2.*sv[i]*yv[i]*(cf+d)+2.*sv[i]*(fg+eg)+2.*cg*yv[i]+g2+a2);
     ras1    [i] = sqrt(yv2[i]+rv2[i]+s1*s1+2.*c*rv[i]*yv[i]+2.*d*s1*yv[i]+2.*e*rv[i]*s1+a2);
     ras2    [i] = sqrt(yv2[i]+(oneph2p2eh)*rv2[i]+2.*rv[i]*yv[i]*(dh+c)+2.*rv[i]*(hm+em)+2.*dm*yv[i]+m2+a2);
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
   int i;
   for(int j = 0; j < 4; j++) {
     // r1
     i=ir1[j];
     b01_r1[i]= log( fabs(rardotq[i]));
     c01_r1[i]= log( fabs(rrardott[i]));
     // r2
     aa=(onepf2p2ef);
     i=ir2[j];
     bb=(yv[i]*(cf+d)+g*(f+e));
     cc=(yv2[i]+2.*cg*yv[i]+g2+a2);
     b01_r2[i]=log(2.*sqrt(aa)*sqrt(aa*sv2[i]+2.*bb*sv[i]+cc)+2.*(aa*sv[i]+bb))/sqrt(aa);
     bb=cg+(cf+d)*sv[i];
     cc=sqrt(-bb*bb+sv2[i]*(onepf2p2ef)+2.*sv[i]*(fg+eg)+g2+a2);
     c01_r2[i]=log(2.*sqrt((yv[i]+bb)*(yv[i]+bb)+cc*cc)+2.*(yv[i]+bb));
     // s1
     i=is1[j];
     a01_s1[i]= log(fabs(rardotp[i]));
     c01_s1[i]= log(fabs(srardott[i])) ;
     // s2
     aa=(oneph2p2eh);
     i=is2[j];
     bb=(yv[i]*(h*d+c)+m*(h+e));
     cc=(yv2[i]+2.*m*d*yv[i]+m2+a2);
     a01_s2[i]=log(2.*sqrt(aa)*sqrt(aa*rv2[i]+2.*bb*rv[i]+cc)+2.*(aa*rv[i]+bb))/sqrt(aa);
     bb=dm+(dh+c)*rv[i];
     cc=sqrt(-bb*bb+rv2[i]*(oneph2p2eh)+2.*rv[i]*(hm+em)+m2+a2);
     c01_s2[i]=log(2.*sqrt((yv[i]+bb)*(yv[i]+bb)+cc*cc)+2.*(yv[i]+bb));
     // linear integrals.
     // abc11
     // r1
     i=ir1[j];
     b11_r1[i]=sra[i]-(d*yv[i]+e*rv[i])*b01_r1[i];
     c11_r1[i]=rar1[i]-(c*rv[i]+d*sfr[i])*c01_r1[i];
     // r2
     i=ir2[j];
     b11_r2[i]=1./(onepf2p2ef)*(sra[i]-(yv[i]*(cf+d)+fg+eg)*b01_r2[i]);
     c11_r2[i]=rar2[i]-(sv[i]*(cf+d)+cg)*c01_r2[i];
     // s1
     i=is1[j];
     a11_s1[i]=rra[i]-(c*yv[i]+e*sv[i])*a01_s1[i];
     c11_s1[i]=ras1[i]-(d*sv[i]+c*rfs[i])*c01_s1[i];
     // s2
     i=is2[j];
     a11_s2[i]=1./(oneph2p2eh)*(rra[i]-(yv[i]*(dh+c)+hm+em)*a01_s2[i]);
     c11_s2[i]=ras2[i]-(rv[i]*(dh+c)+dm)*c01_s2[i];
     // abc1m1
     // r1
     i=ir1[j];
     b0m1_r1[i]=0.5*(rdotq[i]*ra[i] +(ra[i]*ra[i]-(rdotq[i]*rdotq[i]))*b01_r1[i]);
     b1m1_r1[i]=1./3.*(ra[i]*ra[i]*ra[i])-(d*yv[i]+e*rv[i])*b0m1_r1[i];
     c0m1_r1[i]=0.5*(rrdott[i]*rra[i] +(rra[i]*rra[i]-(rrdott[i]*rrdott[i]))*c01_r1[i]);
     c1m1_r1[i]=-1./3.*(-sra[i]*sra[i]*sra[i])-(c*rfs[i]+d*sv[i])*c0m1_r1[i];
     // r2
     i=ir2[j];
     b0m1_r2[i]= 0.5*(sv[i]*sra[i]+(yv[i]*(cf+d)+g*(f+e))*b11_r2[i]+(yv2[i]+2.*g*c*yv[i]+g2+a2)*b01_r2[i]);
     b1m1_r2[i]=1./(onepf2p2ef)*(1./3.*sra[i]*sra[i]*sra[i]-(yv[i]*(cf+d)+g*(f+e))*b0m1_r2[i]);
     c0m1_r2[i]=0.5*(+rar2[i]*(yv[i]+sv[i]*(cf+d)+cg)+((onepf2p2ef-(cf+d)*(cf+d))*sv2[i]+2.*sv[i]*(fg+eg-cg*(cf+d))+(1.-c2)*g2+a2)*c01_r2[i]);
     c1m1_r2[i]=1./3.*rar2[i]*rar2[i]*rar2[i]-(sv[i]*(cf+d)+cg)*c0m1_r2[i];
     // s1
     i=is1[j];
     a0m1_s1[i]=0.5*(rdotp[i]*ra[i] +(ra[i]*ra[i]-(rdotp[i]*rdotp[i]))*a01_s1[i]);
     a1m1_s1[i]=-1./3.*(-ras1[i]*ras1[i]*ras1[i])-(c*yv[i]+e*sv[i])*a0m1_s1[i];
     c0m1_s1[i]=0.5*(srdott[i]*sra[i] +(sra[i]*ra[i]-(srdott[i]*srdott[i]))*c01_s1[i]);
     c1m1_s1[i]=-1./3.*(-rra[i]*rra[i]*rra[i])-(c*rv[i]+d*sfr[i])*c0m1_s1[i];
     // s2
     i=is2[j];
     a0m1_s2[i]= 0.5*(rv[i]*rra[i]+(yv[i]*(dh+c)+m*(h+e))*a11_s2[i]+(yv2[i]+2.*m*d*yv[i]+m2+a2)*a01_s2[i]);
     a1m1_s2[i]=1./(oneph2p2eh)*(1./3.*rra[i]*rra[i]*rra[i]-(yv[i]*(dh+c)+m*(h+e))*a0m1_s2[i]);
     c0m1_s2[i]=0.5*(+ras2[i]*(yv[i]+rv[i]*(dh+c)+dm)+((oneph2p2eh-(dh+c)*(dh+c))*rv2[i]+2.*rv[i]*(hm+em-dm*(dh+c))+(1.-d2)*m2+a2)*c01_s2[i]);
     c1m1_s2[i]=1./3.*ras2[i]*ras2[i]*ras2[i]-(rv[i]*(dh+c)+dm)*c0m1_s2[i];
     // abc0m3
     // r1
     i=ir1[j];
     b0m3_r1[i]=1./4.*(rdotq[i]*pow(rar1[i],3.) +3.*(pow(rar1[i],2.)-(pow(rdotq[i],2.)))*b0m1_r1[i]);
     c0m3_r1[i]=1./4.*(rrdott[i]*pow(rar1[i],3.) +3.*(pow(rar1[i],2.)-pow(rrdott[i],2.))*c0m1_r1[i]);
     // r2
     i=ir2[j];
     vtmp1[i]=((cf+d)*yv[i]+fg+eg)/(onepf2p2ef);
     vtmp2[i]= pow(((cf+d)*yv[i]+fg+eg),2.)/(onepf2p2ef);
     b0m3_r2[i]=-1./4.*(-(sv[i]+vtmp1[i])*pow(rar2[i],3.) -3.*(yv2[i]+g2+2.*cg*yv[i]+a2-vtmp2[i])*b0m1_r2[i]);
     c0m3_r2[i]=1./4.*(pow(rar2[i],3.)*(yv[i]+sv[i]*(cf+d)+cg)+3.*((onepf2p2ef-pow((cf+d),2.))*sv2[i]+2.*sv[i]*(fg+eg-cg*(cf+d))+(1.-c2)*g2+a2)*c0m1_r2[i]);
     // s1
     i=is1[j];
     a0m3_s1[i]=1./4.*(rdotp[i]*pow(ras1[i],3.) +3.*(pow(ras1[i],2.)-(pow(rdotp[i],2.)))*a0m1_s1[i]);
     c0m3_s1[i]=1./4.*(srdott[i]*pow(ras1[i],3.) +3.*(pow(ras1[i],2.)-pow(srdott[i],2.))*c0m1_s1[i]);
     // s2
     i=is2[j];
     vtmp1[i]=((dh+c)*yv[i]+hm+em)/(oneph2p2eh);
     vtmp2[i]= pow((dh+c)*yv[i]+hm+em,2.)/(oneph2p2eh);
     a0m3_s2[i]=-1./4.*(-(rv[i]+vtmp1[i])*pow(ras2[i],3.) -3.*(yv2[i]+m2+2.*dm*yv[i]+a2-vtmp2[i])*a0m1_s2[i]);
     c0m3_s2[i]=1./4.*(pow(ras2[i],3.)*(yv[i]+rv[i]*(dh+c)+dm)+3.*((oneph2p2eh-pow((dh+c),2.))*rv2[i]+2.*rv[i]*(hm+em-dm*(dh+c))+(1.-d2)*m2+a2)*c0m1_s2[i]);
     // abc1m3
     // r1
     i=ir1[j];
     b1m3_r1[i]= 1./5.*pow(rar1[i],5.)-(d*yv[i]+e*r1)*b0m3_r1[i];
     // r2
     i=ir2[j];
     b1m3_r2[i]= 1./5./(onepf2p2ef)*pow(rar2[i],5.)-1./(onepf2p2ef)*(yv[i]*(cf+d)+fg+eg)*b0m3_r2[i];
     // s1
     i=is1[j];
     a1m3_s1[i]= 1./5.*pow(ras1[i],5.)-(c*yv[i]+e*s1)*a0m3_s1[i];
     // s2
     i=is2[j];
     a1m3_s2[i]= -1./5./(oneph2p2eh)*(-pow(ras2[i],5.)+5.*(yv[i]*(dh+c)+hm+em)*a0m3_s2[i]);
     // abc2m1
     // r1
     i=ir1[j];
     b2m1_r1[i]= -1./3.*(-sv[i]*pow(rar1[i],3.)+b0m3_r1[i]-(d*yv[i]+e*r1)*(-3.)*b1m1_r1[i]);
     // r2
     i=ir2[j];
     b2m1_r2[i]= 1./(-3.)/(onepf2p2ef)*(-sv[i]*pow(rar2[i],3.)+b0m3_r2[i]-(yv[i]*(cf+d)+fg+eg)*(-3.)*b1m1_r2[i]);
     // s1
     i=is1[j];
     a2m1_s1[i]= -1./3.*(-rv[i]*pow(ras1[i],3.)+a0m3_s1[i]-(c*yv[i]+e*s1)*(-3.)*a1m1_s1[i]);
     // s2
     i=is2[j];
     a2m1_s2[i]= 1./(-3.)/(oneph2p2eh)*(-rv[i]*pow(ras2[i],3.)+a0m3_s2[i]-(yv[i]*(dh+c)+hm+em)*(-3.)*a1m1_s2[i]);
     // abc21
     // r1
     i=ir1[j];
     b21_r1[i]= -1.*(-sv[i]*pow(rar1[i],1.)+b0m1_r1[i]-(d*yv[i]+e*r1)*(-1.)*b11_r1[i]);
     // r2
     i=ir2[j];
     b21_r2[i]= -1./(onepf2p2ef)*(-sv[i]*pow(rar2[i],1.)+b0m1_r2[i]-(yv[i]*(cf+d)+fg+eg)*(-1.)*b11_r2[i]);
     // s1
     i=is1[j];
     a21_s1[i]= -1.*(-rv[i]*pow(ras1[i],1.)+a0m1_s1[i]-(c*yv[i]+e*s1)*(-1.)*a11_s1[i]);
     // s2
     i=is2[j];
     a21_s2[i]= -1./(oneph2p2eh)*(-rv[i]*pow(ras2[i],1.)+a0m1_s2[i]-(yv[i]*(dh+c)+hm+em)*(-1.)*a11_s2[i]);
     // abc31
     // r1
     i=ir1[j];
     b31_r1[i]= -1.*(-sv[i]*sv[i]*pow(rar1[i],1.)+2.*b1m1_r1[i]-(d*yv[i]+e*r1)*(-1.)*b21_r1[i]);
     // r2
     i=ir2[j];
     b31_r2[i]= -1./(onepf2p2ef)*(-sv[i]*sv[i]*pow(rar2[i],1.)+2.*b1m1_r2[i]-(yv[i]*(cf+d)+fg+eg)*(-1.)*b21_r2[i]);
     // s1
     i=is1[j];
     a31_s1[i]= -1.*(-rv[i]*rv[i]*pow(ras1[i],1.)+2.*a1m1_s1[i]-(c*yv[i]+e*s1)*(-1.)*a21_s1[i]);
     // s2
     i=is2[j];
     a31_s2[i]= -1./(oneph2p2eh)*(-rv[i]*rv[i]*pow(ras2[i],1.)+2.*a1m1_s2[i]-(yv[i]*(dh+c)+hm+em)*(-1.)*a21_s2[i]);
     //
     // DOUBLE INTEGRALS
     //
     /* Double integrals are built from previous linear integrals and a set
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
     // r1
     i=ir1[j];
     root= 1./sqrt((1.-d2)*a2+((1.-d2)*(1.-e2)-(c-de)*(c-de))*rv2[i]);
     f003_r1[i]= 2.*root*atan(((1.-d)*(rra[i]-sv[i]-d*yv[i]-e*rv[i])+(1.-d2)*yv[i]+(c-de)*rv[i])*root);
     root= 1./sqrt((1.-e2)*a2+(1.-c2-d2-e2+2.*c*de)*yv2[i]);
     d003_r1[i]=2.*root*atan(((1.-e)*(rar1[i]+rv[i]-sv[i])+(c-d)*yv[i])*root);
     // r2
     i=ir2[j];
     tp1=sqrt(onepf2p2ef);
     tp2=(yv[i]*(cf+d)+fg+eg)/tp1;
     tp3=sqrt(tp1-cf-d);
     tp4=(tp1*(yv[i]+cg)-tp2*(cf+d))/tp3;
     tp5=sqrt(-tp4*tp4-(tp1+cf+d)*tp2*tp2+(tp1+cf+d)*(a2+g2+yv2[i]+2.*cg*yv[i]));
     f003_r2[i]=2./tp3/tp5*atan((tp3*(rar2[i]-tp1*sv[i]-tp2)+tp4)/tp5);
     tp1=sqrt(onepf2p2ef);
     tp2=sqrt(tp1-f-e);
     root=1./sqrt(yv2[i]*(1.-c2-d2-e2+2.*c*de)+a2*(1.-e2));
     d003_r2[i]=2.*root*atan((tp2*tp2*rar2[i]+(tp1*c-cf-d)*yv[i]+tp1*g-fg-eg-tp1*tp2*tp2*sv[i])*root);
     // s1
     i=is1[j];
     root= 1./sqrt((1.-c2)*a2+((1.-c2)*(1.-e2)-(d-ce)*(d-ce))*sv2[i]);
     e003_s1[i]=2.*root*atan(((1.-c)*(sra[i]-rv[i]-c*yv[i]-e*sv[i])+(1.-c2)*yv[i]+(d-c*e)*sv[i])*root);
     root= 1./sqrt((1.-e2)*a2+(1.-c2-d2-e2+2.*c*de)*yv2[i]);
     d003_s1[i]=2.*root*atan(((1.-e)*(ras1[i]+sv[i]-rv[i])+(d-c)*yv[i])*root);
     // s2
     i=is2[j];
     tp1=sqrt(oneph2p2eh);
     tp2=(yv[i]*(dh+c)+hm+em)/tp1;
     tp3=sqrt(tp1-dh-c);
     tp4=(tp1*(yv[i]+dm)-tp2*(dh+c))/tp3;
     tp5=sqrt(-tp4*tp4-(tp1+dh+c)*tp2*tp2+(tp1+dh+c)*(a2+m2+yv2[i]+2.*dm*yv[i]));
     e003_s2[i]=2./tp3/tp5*atan((tp3*(ras2[i]-tp1*rv[i]-tp2)+tp4)/tp5);
     tp1=sqrt(oneph2p2eh);
     tp2=sqrt(tp1-h-e);
     root=1./sqrt(yv2[i]*(1.-c2-d2-e2+2.*c*de)+a2*(1.-e2));
     d003_s2[i]=2.*root*atan((tp2*tp2*ras2[i]+(tp1*d-dh-c)*yv[i]+tp1*m-hm-em-tp1*tp2*tp2*rv[i])*root);
     // r1
     i=ir1[j];
     f103_r1[i]=1./(1.-d2)*((+d*b01_r1[i] -c01_r1[i]) +(cd-e)*rv[i]*f003_r1[i]);
     f013_r1[i]=1./(1.-d2)*((+d*pow(sfr[i],0.)*c01_r1[i] -pow(yv[i],0.)*b01_r1[i]) +(de-c)*rv[i]*f003_r1[i]);
     // r2
     i=ir2[j];
     f103_r2[i]=1./(onepf2p2ef -(cf+d)*(cf+d))*((cf+d)*pow(yv[i],0.)*b01_r2[i] -c01_r2[i] -(fg+eg-cg*(cf+d))*f003_r2[i]);
     f013_r2[i]=1./(onepf2p2ef -(cf+d)*(cf+d))*((cf+d)*pow(sv[i],0.)*c01_r2[i] -(onepf2p2ef)*pow(yv[i],0.)*b01_r2[i] +((cf+d)*(fg+eg)-(onepf2p2ef)*cg)*f003_r2[i]);
     // s1
     i=is1[j];
     e103_s1[i]=1./(1.-c2)*((+c*a01_s1[i] -c01_s1[i]) +(cd-e)*sv[i]*e003_s1[i]);
     e013_s1[i]=1./(1.-c2)*((+c*pow(rfs[i],0.)*c01_s1[i] -pow(yv[i],0.)*a01_s1[i]) +(ce-d)*sv[i]*e003_s1[i]);
     // s2
     i=is2[j];
     e103_s2[i]=1./(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*pow(yv[i],0.)*a01_s2[i] -c01_s2[i] -(hm+em-dm*(dh+c))*e003_s2[i]);
     e013_s2[i]=1./(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*pow(rv[i],0.)*c01_s2[i] -(oneph2p2eh)*pow(yv[i],0.)*a01_s2[i] +((dh+c)*(hm+em) -(oneph2p2eh)*dm)*e003_s2[i]);
     // r1
     i=ir1[j];
     f001_r1[i] = 1./(1.-d2)/(1.-2.)*(((a2+rv2[i])*(1.-d2)+2.*cd*e*rv2[i]-1.*rv2[i]*(e2+c2))*f003_r1[i]+(-pow(sfr[i],0.)*(sfr[i]*(1.-d2)-c*rv[i]*d+e*rv[i])*c01_r1[i]-pow(yv[i],0.)*(yv[i]*(1.-d2)+c*rv[i]-e*d*rv[i])*b01_r1[i]));
     f101_r1[i]=1./(1.-d2)*(1./(1.-2.)*(+d*pow(yv[i],0.)*b0m1_r1[i]-pow(sfr[i],0.)*c0m1_r1[i])+(cd-e)*rv[i]*f001_r1[i]);
     f011_r1[i]=1./(1.-d2)*(1./(1.-2.)*(+d*pow(sfr[i],0.)*c0m1_r1[i]-pow(yv[i],0.)*b0m1_r1[i])+(de-c)*rv[i]*f001_r1[i]);
     // r2
     i=ir2[j];
     tp1=onepf2p2ef-(cf+d)*(cf+d);
     f001_r2[i] = 1./(1.-2.)*(((fg+eg)*(cf+d)-cg*(onepf2p2ef))/tp1*b01_r2[i]-pow(yv[i],(1.))*b01_r2[i] +(cg*(cf+d)-(fg+eg))/tp1*c01_r2[i]-pow(sv[i],(1.))*c01_r2[i] +(a2+g2)*f003_r2[i] - ((fg+eg)*(fg+eg-cg*(cf+d))+cg*(cg*(onepf2p2ef)-(cf+d)*(fg+eg)))/tp1*f003_r2[i]);
     f101_r2[i]=1./(1.-2.)/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*pow(yv[i],0.)*b0m1_r2[i]-pow(sv[i],0.)*c0m1_r2[i]-(1.-2.)*(fg+eg-cg*(cf+d))*f001_r2[i]);
     f011_r2[i]=1./(1.-2.)/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*pow(sv[i],0.)*c0m1_r2[i]-(onepf2p2ef)*pow(yv[i],0.)*b0m1_r2[i]+(1.-2.)*((cf+d)*(fg+eg)-(onepf2p2ef)*cg)*f001_r2[i]);
     // s1
     i=is1[j];
     e001_s1[i] = 1./(1.-c2)/(1.-2.)*(((a2+sv2[i])*(1.-c2)+2.*c*d*e*sv2[i]-sv2[i]*(e2+d2))*e003_s1[i]+(-pow(rfs[i],0.)*(rfs[i]*(1.-c2)-d*sv[i]*c+e*sv[i])*c01_s1[i]-pow(yv[i],0.)*(yv[i]*(1.-c2)+d*sv[i]-e*c*sv[i])*a01_s1[i]));
     e101_s1[i]=1./(1.-c2)*(1./(1.-2.)*(+c*pow(yv[i],0.)*a0m1_s1[i]-pow(rfs[i],0.)*c0m1_s1[i])+(cd-e)*sv[i]*e001_s1[i]);
     e011_s1[i]=1./(1.-c2)*(1./(1.-2.)*(+c*pow(rfs[i],0.)*c0m1_s1[i]-pow(yv[i],0.)*a0m1_s1[i])+(c*e-d)*sv[i]*e001_s1[i]);
     // s2
     i=is2[j];
     tp1=oneph2p2eh-(dh+c)*(dh+c);
     e001_s2[i] = 1./(1.-2.)*(((hm+em)*(dh+c)-dm*(oneph2p2eh))/tp1*a01_s2[i]-pow(yv[i],(1.))*a01_s2[i] +(dm*(dh+c)-(hm+em))/tp1*c01_s2[i]-pow(rv[i],(1.))*c01_s2[i] +(a2+m2)*e003_s2[i] - ((hm+em)*(hm+em-dm*(dh+c))+dm*(dm*(oneph2p2eh)-(dh+c)*(hm+em)))/tp1*e003_s2[i]);
     e101_s2[i]=1./(1.-2.)/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*pow(yv[i],0.)*a0m1_s2[i]-pow(rv[i],0.)*c0m1_s2[i]-(1.-2.)*(hm+em-dm*(dh+c))*e001_s2[i]);
     e011_s2[i]=1./(1.-2.)/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*pow(rv[i],0.)*c0m1_s2[i]-(oneph2p2eh)*pow(yv[i],0.)*a0m1_s2[i]+(1.-2.)*((dh+c)*(hm+em)-(oneph2p2eh)*dm)*e001_s2[i]);
     // r1
     i=ir1[j];
     f00m1_r1[i] = 1./(1.-d2)/(-3.)*((-(a2+rv2[i])*(1.-d2)+2.*cd*e*-rv2[i]+rv2[i]*(e2+c2))*f001_r1[i]+(-pow(sfr[i],0.)*(sfr[i]*(1.-d2)-c*rv[i]*d+e*rv[i])*c0m1_r1[i]-pow(yv[i],0.)*(yv[i]*(1.-d2)+c*rv[i]-e*d*rv[i])*b0m1_r1[i]));
     // r2
     i=ir2[j];
     tp1=onepf2p2ef-(cf+d)*(cf+d);
     f00m1_r2[i] = 1./(-3.)*(((fg+eg)*(cf+d)-cg*(onepf2p2ef))/tp1*b0m1_r2[i]-pow(yv[i],(1.))*b0m1_r2[i] +(cg*(cf+d)-(fg+eg))/tp1*c0m1_r2[i]-pow(sv[i],(1.))*c0m1_r2[i] -1.*(a2+g2)*f001_r2[i] +1.*((fg+eg)*(fg+eg-cg*(cf+d))+cg*(cg*(onepf2p2ef)-(cf+d)*(fg+eg)))/tp1*f001_r2[i]);
     // s1
     i=is1[j];
     e00m1_s1[i] = 1./(1.-c2)/(-3.)*((-(a2+sv2[i])*(1.-c2)+2.*cd*e*-sv2[i]+sv2[i]*(e2+d2))*e001_s1[i]+(-pow(rfs[i],0.)*(rfs[i]*(1.-c2)-d*sv[i]*c+e*sv[i])*c0m1_s1[i]-pow(yv[i],0.)*(yv[i]*(1.-c2)+d*sv[i]-e*c*sv[i])*a0m1_s1[i]));
     // s2
     i=is2[j];
     tp1=oneph2p2eh-(dh+c)*(dh+c);
     e00m1_s2[i] = 1./(-3.)*(((hm+em)*(dh+c)-dm*(oneph2p2eh))/tp1*a0m1_s2[i]-pow(yv[i],(1.))*a0m1_s2[i] +(dm*(dh+c)-(hm+em))/tp1*c0m1_s2[i]-pow(rv[i],(1.))*c0m1_s2[i] -1.*(a2+m2)*e001_s2[i] +1.*((hm+em)*(hm+em-dm*(dh+c))+dm*(dm*(oneph2p2eh)-(dh+c)*(hm+em)))/tp1*e001_s2[i]);
     // r1
     i=ir1[j];
     f00m3_r1[i] = 1./(1.-d2)/(-5.)*((-3.*(a2+rv2[i])*(1.-d2)+2.*cd*e*-3.*rv2[i]+3.*rv2[i]*(e2+c2))*f00m1_r1[i]+(-pow(sfr[i],0.)*(sfr[i]*(1.-d2)-cd*rv[i]+e*rv[i])*c0m3_r1[i]-pow(yv[i],0.)*(yv[i]*(1.-d2)+c*rv[i]-e*d*rv[i])*b0m3_r1[i]));
     // r2
     i=ir2[j];
     tp1=onepf2p2ef-(cf+d)*(cf+d);
     f00m3_r2[i] = 1./(-5.)*(((fg+eg)*(cf+d)-cg*(onepf2p2ef))/tp1*b0m3_r2[i]-pow(yv[i],(1.))*b0m3_r2[i] +(cg*(cf+d)-(fg+eg))/tp1*c0m3_r2[i]-pow(sv[i],(1.))*c0m3_r2[i] -3.*(a2+g2)*f00m1_r2[i] +3.*((fg+eg)*(fg+eg-cg*(cf+d))+cg*(cg*(onepf2p2ef)-(cf+d)*(fg+eg)))/tp1*f00m1_r2[i]);
     // s1
     i=is1[j];
     e00m3_s1[i] = 1./(1.-c2)/(-5.)*((-3.*(a2+sv2[i])*(1.-c2)+2.*cd*e*-3.*sv2[i]+3.*sv2[i]*(e2+d2))*e00m1_s1[i]+(-pow(rfs[i],0.)*(rfs[i]*(1.-c2)-cd*sv[i]+e*sv[i])*c0m3_s1[i]-pow(yv[i],0.)*(yv[i]*(1.-c2)+d*sv[i]-e*c*sv[i])*a0m3_s1[i]));
     // s2
     i=is2[j];
     tp1=oneph2p2eh-(dh+c)*(dh+c);
     e00m3_s2[i] = 1./(-5.)*(((hm+em)*(dh+c)-dm*(oneph2p2eh))/tp1*a0m3_s2[i]-pow(yv[i],(1.))*a0m3_s2[i] +(dm*(dh+c)-(hm+em))/tp1*c0m3_s2[i]-pow(rv[i],(1.))*c0m3_s2[i] -3.*(a2+m2)*e00m1_s2[i] +3.*((hm+em)*(hm+em-dm*(dh+c))+dm*(dm*(oneph2p2eh)-(dh+c)*(hm+em)))/tp1*e00m1_s2[i]);
     // r1
     i=ir1[j];
     f10m1_r1[i]=1./(1.-d2)*(1./(-3.)*(+d*pow(yv[i],0.)*b0m3_r1[i]-pow(sfr[i],0.)*c0m3_r1[i])+(cd-e)*rv[i]*f00m1_r1[i]);
     // r2
     i=ir2[j];
     f10m1_r2[i]=1./(-3.)/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*pow(yv[i],0.)*b0m3_r2[i]-pow(sv[i],0.)*c0m3_r2[i] -(-3.)*(fg+eg-cg*(cf+d))*f00m1_r2[i]);
     // s1
     i=is1[j];
     e10m1_s1[i]=1./(1.-c2)*(1./(-3.)*(+c*pow(yv[i],0.)*a0m3_s1[i]-pow(rfs[i],0.)*c0m3_s1[i])+(cd-e)*sv[i]*e00m1_s1[i]);
     // s2
     i=is2[j];
     e10m1_s2[i]=1./(-3.)/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*pow(yv[i],0.)*a0m3_s2[i]-pow(rv[i],0.)*c0m3_s2[i] -(-3.)*(hm+em-dm*(dh+c))*e00m1_s2[i]);
     // r1
     i=ir1[j];
     f201_r1[i]=1./(1.-d2)*(1./(-1.)*(f00m1_r1[i]+d*pow(yv[i],0.)*b1m1_r1[i]-pow(sfr[i],1.)*c0m1_r1[i])+(cd-e)*rv[i]*f101_r1[i]);
     f111_r1[i]=1./(1.-d2)*(1./(-1.)*(-d*f00m1_r1[i]+d*pow(sfr[i],1.)*c0m1_r1[i]-pow(yv[i],0.)*b1m1_r1[i])+(de-c)*rv[i]*f101_r1[i]);
     // r2
     i=ir2[j];
     f201_r2[i]=1./(-1.)/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*pow(yv[i],0.)*b1m1_r2[i]-pow(sv[i],1.)*c0m1_r2[i]+1.*f00m1_r2[i]-(-1.)*(fg+eg-cg*(cf+d))*f101_r2[i]);
     f111_r2[i]=1./(-1.)/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*pow(sv[i],1.)*c0m1_r2[i]-(onepf2p2ef)*pow(yv[i],0.)*b1m1_r2[i]-(cf+d)*f00m1_r2[i]+(-1.)*((cf+d)*(fg+eg)-(onepf2p2ef)*cg)*f101_r2[i]);
     // s1
     i=is1[j];
     e201_s1[i]=1./(1.-c2)*(1./(-1.)*(e00m1_s1[i]+c*pow(yv[i],0.)*a1m1_s1[i]-pow(rfs[i],1.)*c0m1_s1[i])+(cd-e)*sv[i]*e101_s1[i]);
     e111_s1[i]=1./(1.-c2)*(1./(1.-2.)*(-c*e00m1_s1[i]+c*pow(rfs[i],1.)*c0m1_s1[i]-pow(yv[i],0.)*a1m1_s1[i])+(ce-d)*sv[i]*e101_s1[i]);
     // s2
     i=is2[j];
     e201_s2[i]=1./(-1.)/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*pow(yv[i],0.)*a1m1_s2[i]-pow(rv[i],1.)*c0m1_s2[i]+e00m1_s2[i]-(-1.)*(hm+em-dm*(dh+c))*e101_s2[i]);
     e111_s2[i]=1./(1.-2.)/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*pow(rv[i],1.)*c0m1_s2[i]-(oneph2p2eh)*pow(yv[i],0.)*a1m1_s2[i]-(dh+c)*e00m1_s2[i]+(-1.)*((dh+c)*(hm+em)-(oneph2p2eh)*dm)*e101_s2[i]);
     // r1
     i=ir1[j];
     f203_r1[i]=1./(1.-d2)*((f001_r1[i]+d*pow(yv[i],0.)*b11_r1[i]-pow(sfr[i],1.)*c01_r1[i])+(cd-e)*rv[i]*f103_r1[i]);
     // r2
     i=ir2[j];
     f203_r2[i]=1./(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*pow(yv[i],0.)*b11_r2[i]-pow(sv[i],1.)*c01_r2[i]+f001_r2[i]-(fg+eg-cg*(cf+d))*f103_r2[i]);
     // s1
     i=is1[j];
     e203_s1[i]=1./(1.-c2)*((e001_s1[i]+c*pow(yv[i],0.)*a11_s1[i]-pow(rfs[i],1.)*c01_s1[i])+(cd-e)*sv[i]*e103_s1[i]);
     // s2
     i=is2[j];
     e203_s2[i]=1./(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*pow(yv[i],0.)*a11_s2[i]-pow(rv[i],1.)*c01_s2[i]+e001_s2[i]-(hm+em-dm*(dh+c))*e103_s2[i]);
     // r1
     i=ir1[j];
     f023_r1[i]=1./(1.-d2)*((f001_r1[i]+d*pow(sfr[i],0.)*c11_r1[i]-pow(yv[i],1.)*b01_r1[i])+(de-c)*rv[i]*f013_r1[i]);
     // r2
     i=ir2[j];
     f023_r2[i]=1./(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*pow(sv[i],0.)*c11_r2[i]-(onepf2p2ef)*pow(yv[i],1.)*b01_r2[i]+(onepf2p2ef)*f001_r2[i]+((cf+d)*(fg+eg)-(onepf2p2ef)*cg)*f013_r2[i]);
     // s1
     i=is1[j];
     e023_s1[i]=1./(1.-c2)*((e001_s1[i]+c*pow(rfs[i],0.)*c11_s1[i]-pow(yv[i],1.)*a01_s1[i])+(ce-d)*sv[i]*e013_s1[i]);
     // s2
     i=is2[j];
     e023_s2[i]=1./(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*pow(rv[i],0.)*c11_s2[i]-(oneph2p2eh)*pow(yv[i],1.)*a01_s2[i]+(oneph2p2eh)*e001_s2[i]+((dh+c)*(hm+em)-(oneph2p2eh)*dm)*e013_s2[i]);
     // r1
     i=ir1[j];
     f113_r1[i]=1./(1.-d2)*((-d*f001_r1[i]+d*pow(yv[i],1.)*b01_r1[i]-pow(sfr[i],0.)*c11_r1[i])+(cd-e)*rv[i]*f013_r1[i]);
     // r2
     i=ir2[j];
     f113_r2[i]=1./(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*pow(yv[i],1.)*b01_r2[i]-pow(sv[i],0.)*c11_r2[i]-(cf+d)*f001_r2[i]-(fg+eg-cg*(cf+d))*f013_r2[i]);
     // s1
     i=is1[j];
     e113_s1[i]=1./(1.-c2)*((-c*e001_s1[i]+c*pow(yv[i],1.)*a01_s1[i]-pow(rfs[i],0.)*c11_s1[i])+(cd-e)*sv[i]*e013_s1[i]);
     // s2
     i=is2[j];
     e113_s2[i]=1./(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*pow(yv[i],1.)*a01_s2[i]-pow(rv[i],0.)*c11_s2[i]-(dh+c)*e001_s2[i]-(hm+em-dm*(dh+c))*e013_s2[i]);
     // r1
     i=ir1[j];
     f213_r1[i]=1./(1.-d2)*((f011_r1[i]-d*f101_r1[i]+d*pow(yv[i],1.)*b11_r1[i]-pow(sfr[i],1.)*c11_r1[i])+(cd-e)*rv[i]*f113_r1[i]);
     f123_r1[i]=1./(1.-d2)*((f101_r1[i]-d*f011_r1[i]+d*pow(sfr[i],1.)*c11_r1[i]-pow(yv[i],1.)*b11_r1[i])+(de-c)*rv[i]*f113_r1[i]);
     // r2
     i=ir2[j];
     f213_r2[i]=1./(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*pow(yv[i],1.)*b11_r2[i]-pow(sv[i],1.)*c11_r2[i]+f011_r2[i]-(cf+d)*f101_r2[i]-(fg+eg-cg*(cf+d))*f113_r2[i]);
     f123_r2[i]=1./(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*pow(sv[i],1.)*c11_r2[i]-(onepf2p2ef)*pow(yv[i],1.)*b11_r2[i]+(onepf2p2ef)*f101_r2[i]-(cf+d)*f011_r2[i]+((cf+d)*(fg+eg)-(onepf2p2ef)*cg)*f113_r2[i]);
     // s1
     i=is1[j];
     e213_s1[i]=1./(1.-c2)*((e011_s1[i]-c*e101_s1[i]+c*pow(yv[i],1.)*a11_s1[i]-pow(rfs[i],1.)*c11_s1[i])+(cd-e)*sv[i]*e113_s1[i]);
     e123_s1[i]=1./(1.-c2)*((e101_s1[i]-c*e011_s1[i]+c*pow(rfs[i],1.)*c11_s1[i]-pow(yv[i],1.)*a11_s1[i])+(ce-d)*sv[i]*e113_s1[i]);
     // s2
     i=is2[j];
     e213_s2[i]=1./(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*pow(yv[i],1.)*a11_s2[i]-pow(rv[i],1.)*c11_s2[i]+e011_s2[i]-(dh+c)*e101_s2[i]-(hm+em-dm*(dh+c))*e113_s2[i]);
     e123_s2[i]=1./(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*pow(rv[i],1.)*c11_s2[i]-(oneph2p2eh)*pow(yv[i],1.)*a11_s2[i]+(oneph2p2eh)*e101_s2[i]-(dh+c)*e011_s2[i]+((dh+c)*(hm+em)-(oneph2p2eh)*dm)*e113_s2[i]);
     // r1
     i=ir1[j];
     f313_r1[i]=1./(1.-d2)*((2.*f111_r1[i]-d*f201_r1[i]+d*pow(yv[i],1.)*b21_r1[i]-pow(sfr[i],2.)*c11_r1[i])+(cd-e)*rv[i]*f213_r1[i]);
     f223_r1[i]=1./(1.-d2)*((f201_r1[i]-d*2.*f111_r1[i]+d*pow(sfr[i],2.)*c11_r1[i]-pow(yv[i],1.)*b21_r1[i])+(de-c)*rv[i]*f213_r1[i]);
     // r2
     i=ir2[j];
     f313_r2[i]=1./(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*pow(yv[i],1.)*b21_r2[i]-pow(sv[i],2.)*c11_r2[i]+2.*f111_r2[i]-(cf+d)*f201_r2[i]-(fg+eg-cg*(cf+d))*f213_r2[i]);
     f223_r2[i]=1./(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*pow(sv[i],2.)*c11_r2[i]-(onepf2p2ef)*pow(yv[i],1.)*b21_r2[i]+(onepf2p2ef)*f201_r2[i]-2.*(cf+d)*f111_r2[i]+((cf+d)*(fg+eg)-(onepf2p2ef)*cg)*f213_r2[i]);
     // s1
     i=is1[j];
     e313_s1[i]=1./(1.-c2)*((2.*e111_s1[i]-c*e201_s1[i]+c*pow(yv[i],1.)*a21_s1[i]-pow(rfs[i],2.)*c11_s1[i])+(cd-e)*sv[i]*e213_s1[i]);
     e223_s1[i]=1./(1.-c2)*((e201_s1[i]-c*2.*e111_s1[i]+c*pow(rfs[i],2.)*c11_s1[i]-pow(yv[i],1.)*a21_s1[i])+(c*e-d)*sv[i]*e213_s1[i]);
     // s2
     i=is2[j];
     e313_s2[i]=1./(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*pow(yv[i],1.)*a21_s2[i]-pow(rv[i],2.)*c11_s2[i]+2.*e111_s2[i]-(dh+c)*e201_s2[i]-(hm+em-dm*(dh+c))*e213_s2[i]);
     e223_s2[i]=1./(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*pow(rv[i],2.)*c11_s2[i]-(oneph2p2eh)*pow(yv[i],1.)*a21_s2[i]+(oneph2p2eh)*e201_s2[i]-2.*(dh+c)*e111_s2[i]+((dh+c)*(hm+em)-(oneph2p2eh)*dm)*e213_s2[i]);
   }

   for(int i = 0; i < 8; i++) {
     d003[i]= 0.5*d003_r1[i] +0.5*d003_r2[i] +0.5*d003_s1[i] +0.5*d003_s2[i];
     aas001[i]= a01_s1[i]+a01_s2[i];
     bbr001[i]= b01_r1[i]+b01_r2[i];
     d103[i]=1./(1.-e2)*((+e*aas001[i] -bbr001[i]) -(c-de)*yv[i]*d003[i]);
     d013[i]=1./(1.-e2)*((+e*bbr001[i] -aas001[i]) -(d-ce)*yv[i]*d003[i]);
     d001[i] =1./((1.-e2)*(-1.))*( ((a2+yv2[i]) * (1.-e2) +2.*cd*e*1. * yv2[i] - 1.*yv2[i]*(c2+d2)) * d003[i]-(yv[i]*(d-ce) +pow(m,(1.))*(1.-e2)) * a01_s2[i]-(h*(1.-e2))*a11_s2[i]-(yv[i]*(d-ce)+pow(sv[i],(1.))*(1.-e2)) * a01_s1[i]-(yv[i]*(c-de)+(1.-e2)*g)*b01_r2[i]-(f*(1.-e2))*b11_r2[i]-(yv[i]*(c-e*d)+pow(rv[i],(1.))*(1.-e2)) * b01_r1[i]);
     aas00m1[i]= a0m1_s1[i]+a0m1_s2[i];
     bbr00m1[i]= b0m1_r1[i]+b0m1_r2[i];
     d101[i]=1./(1.-e2)*(1./(1.-2.)*(+e*aas00m1[i]-bbr00m1[i])-(c-de)*yv[i]*d001[i]);
     d011[i]=1./(1.-e2)*(1./(1.-2.)*(+e*bbr00m1[i]-aas00m1[i])-(d-ce)*yv[i]*d001[i]);
     d00m1[i] =1./((1.-e2)*(-3.))*( (-(a2+yv2[i]) * (1.-e2) +2.*cd*e*-1. * yv2[i] +1.*yv2[i]*(c2+d2)) * d001[i]-(yv[i]*(d-e*c) +pow(m,(1.))*(1.-e2)) * a0m1_s2[i]-(h*(1.-e2))*a1m1_s2[i]-(yv[i]*(d-e*c)+pow(sv[i],(1.))*(1.-e2)) * a0m1_s1[i]-(yv[i]*(c-e*d)+pow(g,(1.))*(1.-e2)) * b0m1_r2[i]-(f*(1.-e2))*b1m1_r2[i]-(yv[i]*(c-e*d)+pow(rv[i],(1.))*(1.-e2)) * b0m1_r1[i]);
     d00m3[i] =1./((1.-e2)*(-5.))*( (-3.*(a2+yv2[i]) * (1.-e2) +2.*cd*e*-3. * yv2[i] +3.*yv2[i]*(c2+d2)) * d00m1[i]-(yv[i]*(d-e*c) +pow(m,(1.))*(1.-e2)) * a0m3_s2[i]-(h*(1.-e2))*a1m3_s2[i]-(yv[i]*(d-e*c)+pow(sv[i],(1.))*(1.-e2)) * a0m3_s1[i]-(yv[i]*(c-e*d)+pow(g,(1.))*(1.-e2)) * b0m3_r2[i]-(f*(1.-e2))*b1m3_r2[i]-(yv[i]*(c-e*d)+pow(rv[i],(1.))*(1.-e2)) * b0m3_r1[i]);
     aas00m3[i]= a0m3_s1[i]+a0m3_s2[i];
     bbr00m3[i]= b0m3_r1[i]+b0m3_r2[i];
     d10m1[i]=1./(1.-e2)*(1./(-3.)*(+e*aas00m3[i]-bbr00m3[i])-(c-de)*yv[i]*d00m1[i]);
     d01m1[i]=1./(1.-e2)*(1./(-3.)*(+e*bbr00m3[i]-aas00m3[i])-(d-ce)*yv[i]*d00m1[i]);
     aas10m1[i]= a1m1_s1[i]+a1m1_s2[i];
     aas01m1[i]= h*a1m1_s2[i] +m* a0m1_s2[i] +s1*a0m1_s1[i];
     aas11m1[i]= h*a2m1_s2[i] +m* a1m1_s2[i] +s1*a1m1_s1[i];
     bbr01m1[i]= +b1m1_r2[i] +b1m1_r1[i];
     bbr10m1[i]= f*b1m1_r2[i] +g* b0m1_r2[i] +r1*b0m1_r1[i];
     bbr11m1[i]= f*b2m1_r2[i] +g* b1m1_r2[i] +r1*b1m1_r1[i];
     d201[i]= 1./(1.-e2)*(1./(-1.)*(d00m1[i] +e*aas10m1[i] -bbr10m1[i]) -(c-de)*yv[i]*d101[i]);
     d111_1[i]= 1./(1.-e2)*(1./(-1.)*(-e*d00m1[i] +e*aas01m1[i] -bbr01m1[i]) -(c-de)*yv[i]*d011[i]);
     d111_2[i]= 1./(1.-e2)*(1./(-1.)*(-e*d00m1[i] +e*bbr10m1[i] -aas10m1[i]) -(d-ce)*yv[i]*d101[i]);
     d021[i]=1./(1.-e2)*(1./(-1.)*(d00m1[i]+e*bbr01m1[i]-aas01m1[i])-(d-ce)*yv[i]*d011[i]);
     d111[i]= 0.5*d111_1[i] +0.5*d111_2[i];
     aas101[i]= a11_s1[i]+a11_s2[i];
     aas011[i]= h*a11_s2[i] +m* a01_s2[i] +s1*a01_s1[i];
     bbr101[i]= f*b11_r2[i] +g* b01_r2[i] +r1*b01_r1[i];
     bbr011[i]= +b11_r2[i] +b11_r1[i];
     d203[i]=1./(1.-e2)*((d001[i]+e*aas101[i]-bbr101[i])-(c-de)*yv[i]*d103[i]);
     d023[i]=1./(1.-e2)*((d001[i]+e*bbr011[i]-aas011[i])-(d-ce)*yv[i]*d013[i]);
     d113_1[i]=1./(1.-e2)*((-e*d001[i]+e*aas011[i]-bbr011[i])-(c-de)*yv[i]*d013[i]);
     d113_2[i]=1./(1.-e2)*((-e*d001[i]+e*bbr101[i]-aas101[i])-(d-ce)*yv[i]*d103[i]);
     d113[i]= 0.5*d113_1[i] +0.5*d113_2[i];
     aas111[i] = h*a21_s2[i] +m* a11_s2[i] +s1*a11_s1[i];
     bbr111[i] = f*b21_r2[i] +g* b11_r2[i] +r1*b11_r1[i];
     d213[i]=1./(1.-e2)*((d011[i]-e*d101[i]+e*aas111[i]-bbr111[i])-(c-de)*yv[i]*d113[i]);
     d123[i]=1./(1.-e2)*((d101[i]-e*d011[i]+e*bbr111[i]-aas111[i])-(d-ce)*yv[i]*d113[i]);
     aas211[i] = h*a31_s2[i] +m* a21_s2[i] +s1*a21_s1[i];
     aas121[i] = h2*a31_s2[i] +2.*hm*a21_s2[i] +m2* a11_s2[i] +s1*s1*a11_s1[i];
     bbr211[i] = f2*b31_r2[i] +2.*fg*b21_r2[i] +g2* b11_r2[i] +r1*r1*b11_r1[i];
     bbr121[i] = f*b31_r2[i] +g* b21_r2[i] +r1*b21_r1[i];
     d313[i]=1./(1.-e2)*((2.*d111[i]-e*d201[i]+e*aas211[i]-bbr211[i])-(c-de)*yv[i]*d213[i]);
     d223_1[i]=1./(1.-e2)*((d201[i]-e*2.*d111[i]+e*bbr211[i]-aas211[i])-(d-ce)*yv[i]*d213[i]);
     d223_2[i]=1./(1.-e2)*((d021[i]-e*2.*d111[i]+e*aas121[i]-bbr121[i])-(c-de)*yv[i]*d123[i]);
     d133[i]=1./(1.-e2)*((2.*d111[i]-e*d021[i]+e*bbr121[i]-aas121[i])-(d-ce)*yv[i]*d123[i]);
     aas201[i] = a21_s1[i]+a21_s2[i];
     aas021[i] = h2*a21_s2[i] +2.*hm*a11_s2[i] +m2* a01_s2[i] +s1*s1*a01_s1[i];
     bbr201[i] = f2*b21_r2[i] +2.*fg*b11_r2[i] +g2* b01_r2[i] +r1*r1*b01_r1[i];
     bbr021[i] = +b21_r2[i] +b21_r1[i];
     d303[i]=1./(1.-e2)*((2.*d101[i]+e*aas201[i]-bbr201[i])-(c-de)*yv[i]*d203[i]);
     d033[i]=1./(1.-e2)*((2.*d011[i]+e*bbr021[i]-aas021[i])-(d-ce)*yv[i]*d023[i]);
     d223[i]= 0.5*d223_1[i] +0.5*d223_2[i];
     //
     // TRIPLE INTEGRALS
     //
     /* triple integrals are built from previous double integrals.
        Hijkl = /int/int/int r^i s^j y^k Ra^-l dr ds dy
        Since r2 is a fct(s) and s2 is a fct(r) in the case of a triangle, we introduce
        additional double integrals:
        EESijkl = /int/int r^i s(r)^j y^k Ra^-l dr dy = /int/int r^i (hr+m)^j y^k Ra^-l drdy - /int/int r^i sa^j y^k Ra^-l drdy
        FFRijkl = /int/int r(s)^i s^j y^k Ra^-l dr dy = /int/int (fs+g)^i s^j y^k Ra^-l drdy - /int/int r1^i s^j y^k Ra^-l drdy
     */

     /* SEED FUNCTION */
     /* since no analytical was found for the only tripple seed integral,
        we perform instead an adaptative Simpson quadrature */
     //  functionPtr = &FDr1s1;
     //  Dr1s1 = adaptiveSimpsons(functionPtr, y1, y2, 1e-9, 20);
     //  functionPtr = &FDr1s2;
     //  Dr1s2 = adaptiveSimpsons(functionPtr, y1, y2, 1e-9, 20);
     //  functionPtr = &FDr2s1;
     //  Dr2s1 = adaptiveSimpsons(functionPtr, y1, y2, 1e-9, 20);
     //  functionPtr = &FDr2s2;
     //  Dr2s2 = adaptiveSimpsons(functionPtr, y1, y2, 1e-9, 20);
     //
     //  H0003[0] = Dr2s2; H0003[1] = 0;
     //  H0003[2] = Dr2s1; H0003[3] = 0;
     //  H0003[4] = Dr1s2; H0003[5] = 0;
     //  H0003[6] = Dr1s1; H0003[7] = 0;
     ffr1001[i] = f*f101_r2[i] +g*f001_r2[i] +r1*f001_r1[i];
     ees0101[i] = h*e101_s2[i] +m*e001_s2[i] +s1*e001_s1[i];
     h0001[i] = 1./(-2.)*(a2*0.-ffr1001[i]-ees0101[i]-pow(yv[i],1.)*d001[i]);
     ffr100m1[i] = f*f10m1_r2[i] +g*f00m1_r2[i] +r1*f00m1_r1[i];
     ees010m1[i] = h*e10m1_s2[i] +m*e00m1_s2[i] +s1*e00m1_s1[i];
     h000m1[i] = 1./(-4.)*(-a2*h0001[i]-ffr100m1[i]-ees010m1[i]-pow(yv[i],1.)*d00m1[i]);
     ees000m1[i] = e00m1_s1[i] +e00m1_s2[i];
     ffr000m1[i] = f00m1_r2[i] +f00m1_r1[i];
     h0011[i] = (1./-1./(1.-e2-c2-d2+2.*cde))*((1.-e2)*( -pow(yv[i],0.)*d00m1[i]) +(d-ce)*ees000m1[i] +(c-de)*ffr000m1[i]);
     h1001[i] = 1./-1./(1.-e2)*(+e*ees000m1[i] -ffr000m1[i] -1.*(de-c)*h0011[i]);
     h0101[i] = 1./-1./(1.-e2)*(-ees000m1[i] +e*ffr000m1[i] -1.*(ce-d)*h0011[i]);
     ees100m1[i] = e10m1_s1[i] +e10m1_s2[i];
     ffr010m1[i] = f10m1_r2[i] +f10m1_r1[i];
     h1011[i] = (1./-1./(1.-e2-c2-d2+2.*cde))*((1.-e2)*( -pow(yv[i],0.)*d10m1[i]) -(c-de)*h000m1[i] +(d-ce)*ees100m1[i] +(c-de)*ffr100m1[i]);
     h0111[i] = (1./-1./(1.-e2-c2-d2+2.*cde))*((1.-e2)*( -pow(yv[i],0.)*d01m1[i]) -(d-ce)*h000m1[i] +(d-ce)*ees010m1[i] +(c-de)*ffr010m1[i]);
     h1101[i] = -1./(1.-e2)*( -e*h000m1[i] -ees100m1[i] +e*ffr100m1[i] -(ce-d)*h1011[i]);
     //  h1101_2[i] = 1./-1./(1.-e2)*( -e*h000m1[i] +e*ees010m1[i] -ffr010m1[i] -(de-c)*h0111[i]);
     h2001[i] = 1./-1./(1.-e2)*(1.*h000m1[i] +e*ees100m1[i] -ffr100m1[i] -(de-c)*h1011[i]);
     h0201[i] = 1./-1./(1.-e2)*(1.*h000m1[i] -ees010m1[i] +e*ffr010m1[i] -(ce-d)*h0111[i]);
     ees0001[i] = e001_s1[i] +e001_s2[i];
     ffr0001[i] = f001_r2[i] +f001_r1[i];
     ees1001[i] = e101_s1[i] +e101_s2[i];
     ffr0101[i] = f101_r2[i] +f101_r1[i];
     hijkl[0][i] = (1./(1.-e2-c2-d2+2.*cde))*((1.-e2)*( -pow(yv[i],0.)*d001[i]) +(d-ce)*ees0001[i] +(c-de)*ffr0001[i]);
     hijkl[1][i] = 1./(1.-e2)*( +e*ees0001[i] -ffr0001[i] +(de-c)*hijkl[0][i]);
     hijkl[2][i] = 1./(1.-e2)*( -ees0001[i] +e*ffr0001[i] +(ce-d)*hijkl[0][i]);
     hijkl[5][i] = (1./(1.-e2-c2-d2+2.*cde))*((1.-e2)*( -pow(yv[i],0.)*d101[i]) -(c-de)*h0001[i] +(d-ce)*ees1001[i] +(c-de)*ffr1001[i]);
     hijkl[6][i] = (1./(1.-e2-c2-d2+2.*cde))*((1.-e2)*( -pow(yv[i],0.)*d011[i])  -(d-ce)*h0001[i] +(d-ce)*ees0101[i] +(c-de)*ffr0101[i]);
     h1103_1[i] = 1./(1.-e2)*( -e*h0001[i] +e*ees0101[i] -ffr0101[i] +(de-c)*hijkl[6][i]);
     h1103_2[i] = 1./(1.-e2)*( -e*h0001[i] +e*ffr1001[i] -ees1001[i] +(ce-d)*hijkl[5][i]);
     hijkl[7][i]= 0.5*h1103_1[i] +0.5*h1103_2[i];
     ees0011[i]= e011_s1[i] +e011_s2[i];
     ffr0011[i]= f011_r2[i] +f011_r1[i];
     h0023[i] = (1./(1.-e2-c2-d2+2.*cde))*((1.-e2)*(h0001[i] -pow(yv[i],1.)*d001[i]) +(d-ce)*ees0011[i] +(c-de)*ffr0011[i]);
     hijkl[3][i] = 1./(1.-e2)*(h0001[i] +e*ees1001[i] -ffr1001[i] +(de-c)*hijkl[5][i]);
     hijkl[4][i] = 1./(1.-e2)*(h0001[i] -ees0101[i] +e*ffr0101[i] +(ce-d)*hijkl[6][i]);
     ees2001[i]= e201_s1[i] +e201_s2[i];
     ffr0201[i]= f201_r2[i] +f201_r1[i];
     ffr2001[i]= f2*f201_r2[i] +2.*fg*f101_r2[i] +g2*f001_r2[i] +r12*f001_r1[i];
     ees0201[i]= h2*e201_s2[i] +2.*hm*e101_s2[i] +m2*e001_s2[i] +s12*e001_s1[i];
     hijkl[8][i] = (1./1./(1.-e2-c2-d2+2.*cde))*((1.-e2)*( -pow(yv[i],0.)*d201[i]) -2.*(c-de)*h1001[i] +(d-ce)*ees2001[i] +(c-de)*ffr2001[i]);
     hijkl[9][i] = (1./1./(1.-e2-c2-d2+2.*cde))*((1.-e2)*( -pow(yv[i],0.)*d021[i]) -2.*(d-ce)*h0101[i] +(d-ce)*ees0201[i] +(c-de)*ffr0201[i]);
     ees1101[i]= h*e201_s2[i] +m*e101_s2[i] +s1*e101_s1[i];
     ffr1101[i]= f*f201_r2[i] +g*f101_r2[i] +r1*f101_r1[i];
     hijkl[12][i] = (1./1./(1.-e2-c2-d2+2.*cde))*((1.-e2)*( -pow(yv[i],0.)*d111[i]) -(c-de)*h0101[i] -(d-ce)*h1001[i] +(d-ce)*ees1101[i] +(c-de)*ffr1101[i]);
     hijkl[10][i] = 1./(1.-e2)*(h0101[i] -e*h1001[i] +e*ees1101[i] -ffr1101[i] +(de-c)*hijkl[12][i]);
     hijkl[11][i] = 1./(1.-e2)*(h1001[i] -e*h0101[i] +e*ffr1101[i] -ees1101[i]  +(ce-d)*hijkl[12][i]);
     ees0111[i] = h*e111_s2[i] +m*e011_s2[i] +s1*e011_s1[i];
     ffr1011[i] = f*f111_r2[i] +g*f011_r2[i] +r1*f011_r1[i];
     ffr0111[i] = f111_r2[i] +f111_r1[i];
     ees1011[i] = e111_s1[i] +e111_s2[i];
     h0123[i] = (1./(1.-e2-c2-d2+2.*cde))*((1.-e2)*(h0101[i] -pow(yv[i],1.)*d011[i]) -(d-ce)*h0011[i] +(d-ce)*ees0111[i] +(c-de)*ffr0111[i]);
     h1023[i] = (1./(1.-e2-c2-d2+2.*cde))*((1.-e2)*(h1001[i] -pow(yv[i],1.)*d101[i]) -(c-de)*h0011[i] +(d-ce)*ees1011[i] +(c-de)*ffr1011[i]);
     hijkl[13][i] = 1./(1.-e2)*(2.*h1001[i] +e*ees2001[i] -ffr2001[i] +(de-c)*hijkl[8][i]);
     hijkl[14][i] = 1./(1.-e2)*(2.*h0101[i] -ees0201[i] +e*ffr0201[i] +(ce-d)*hijkl[9][i]);
     ees0003[i] = e003_s1[i] +e003_s2[i];
     ffr0003[i] = f003_r2[i] +f003_r1[i];
     hijkl[15][i] = (1./3./(1.-e2-c2-d2+2.*cde))*((1.-e2)*( -pow(yv[i],0.)*d003[i]) +(d-ce)*ees0003[i] +(c-de)*ffr0003[i]);
     hijkl[16][i] = 1./3./(1.-e2)*( +e*ees0003[i] -ffr0003[i] +3.*(de-c)*hijkl[15][i]);
     hijkl[17][i] = 1./3./(1.-e2)*( -ees0003[i] +e*ffr0003[i] +3.*(ce-d)*hijkl[15][i]);
     ees1003[i] = e103_s1[i] +e103_s2[i];
     ffr1003[i] = f*f103_r2[i] +g*f003_r2[i] +r1*f003_r1[i];
     ees0103[i] = h*e103_s2[i] +m*e003_s2[i] +s1*e003_s1[i];
     ffr0103[i] = f103_r2[i] +f103_r1[i];
     hijkl[20][i] = (1./3./(1.-e2-c2-d2+2.*cde))*((1.-e2)*( -pow(yv[i],0.)*d103[i]) -(c-de)*0. +(d-ce)*ees1003[i] +(c-de)*ffr1003[i]);
     hijkl[21][i] = (1./3./(1.-e2-c2-d2+2.*cde))*((1.-e2)*( -pow(yv[i],0.)*d013[i]) -(d-ce)*0. +(d-ce)*ees0103[i] +(c-de)*ffr0103[i]);
     h1105_1[i] = 1./3./(1.-e2)*( -e*0. +e*ees0103[i] -ffr0103[i] +3.*(de-c)*hijkl[21][i]);
     h1105_2[i] = 1./3./(1.-e2)*( -e*0. +e*ffr1003[i] -ees1003[i] +3.*(ce-d)*hijkl[20][i]);
     hijkl[22][i]= 0.5*h1105_1[i] +0.5*h1105_2[i];
     ees0013[i] = e013_s1[i] +e013_s2[i];
     ffr0013[i] = f013_r2[i] +f013_r1[i];
     h0025[i] = (1./3./(1.-e2-c2-d2+2.*cde))*((1.-e2)*(0. - pow(yv[i],1.)*d003[i]) +(d-ce)*ees0013[i] +(c-de)*ffr0013[i]);
     hijkl[18][i] = 1./3./(1.-e2)*(0. +e*ees1003[i] -ffr1003[i] +3.*(de-c)*hijkl[20][i]);
     hijkl[19][i] = 1./3./(1.-e2)*(0. -ees0103[i] +e*ffr0103[i] +3.*(ce-d)*hijkl[21][i]);
     ees1103[i] = h*e203_s2[i] +m*e103_s2[i] +s1*e103_s1[i];
     ffr1103[i] = f*f203_r2[i] +g*f103_r2[i] +r1*f103_r1[i];
     hijkl[23][i] = (1./3./(1.-e2-c2-d2+2.*cde))*((1.-e2)*( -pow(yv[i],0.)*d113[i]) -(c-de)*hijkl[2][i] -(d-ce)*hijkl[1][i] +(d-ce)*ees1103[i] +(c-de)*ffr1103[i]);
     hijkl[35][i] = 1./3./(1.-e2)*(hijkl[2][i] -e*hijkl[1][i] +e*ees1103[i] -ffr1103[i] +3.*(de-c)*hijkl[23][i]);
     hijkl[36][i] = 1./3./(1.-e2)*(hijkl[1][i] -e*hijkl[2][i] -ees1103[i] +e*ffr1103[i] +3.*(ce-d)*hijkl[23][i]);
     ees1013[i] = e113_s1[i] +e113_s2[i];
     ffr1013[i] = f*f113_r2[i] +g*f013_r2[i] +r1*f013_r1[i];
     hijkl[27][i] = (1./3./(1.-e2-c2-d2+2.*cde))*((1.-e2)*(hijkl[1][i] -pow(yv[i],1.)*d103[i]) -(c-de)*hijkl[0][i] +(d-ce)*ees1013[i] +(c-de)*ffr1013[i]);
     hijkl[29][i] = 1./3./(1.-e2)*(hijkl[0][i] +e*ees1013[i] -ffr1013[i] +3.*(de-c)*hijkl[27][i]);
     ees0113[i] = h*e113_s2[i] +m*e013_s2[i] +s1*e013_s1[i];
     ffr0113[i] = f113_r2[i] +f113_r1[i];
     hijkl[28][i] = (1./3./(1.-e2-c2-d2+2.*cde))*((1.-e2)*(hijkl[2][i] -pow(yv[i],1.)*d013[i]) -(d-ce)*hijkl[0][i] +(d-ce)*ees0113[i] +(c-de)*ffr0113[i]);
     hijkl[30][i] = 1./3./(1.-e2)*(hijkl[0][i] -ees0113[i] +e*ffr0113[i] +3.*(ce-d)*hijkl[28][i]);
     ees2013[i] = e213_s1[i] +e213_s2[i];
     ffr2013[i] = f2*f213_r2[i] +2.*fg*f113_r2[i] +g2*f013_r2[i] +r12*f013_r1[i];
     ees0213[i] = h2*e213_s2[i] +2.*hm*e113_s2[i] +m2*e013_s2[i] +s12*e013_s1[i];
     ffr0213[i] = f213_r2[i] +f213_r1[i];
     hijkl[31][i] = (1./3./(1.-e2-c2-d2+2.*cde))*((1.-e2)*(hijkl[3][i] -pow(yv[i],1.)*d203[i]) -2.*(c-de)*hijkl[5][i] +(d-ce)*ees2013[i] +(c-de)*ffr2013[i]);
     hijkl[32][i] = (1./3./(1.-e2-c2-d2+2.*cde))*((1.-e2)*(hijkl[4][i] -pow(yv[i],1.)*d023[i]) -2.*(d-ce)*hijkl[6][i] +(d-ce)*ees0213[i] +(c-de)*ffr0213[i]);
     ees1113[i] = h*e213_s2[i] +m*e113_s2[i] +s1*e113_s1[i];
     ffr1113[i] = f*f213_r2[i] +g*f113_r2[i] +r1*f113_r1[i];
     hijkl[24][i] = (1./3./(1.-e2-c2-d2+2.*cde))*((1.-e2)*(hijkl[7][i] -pow(yv[i],1.)*d113[i]) -(c-de)*hijkl[6][i] -(d-ce)*hijkl[5][i] +(d-ce)*ees1113[i] +(c-de)*ffr1113[i]);
     hijkl[25][i] = 1./3./(1.-e2)*(hijkl[6][i] -e*hijkl[5][i] +e*ees1113[i] -ffr1113[i] +3.*(de-c)*hijkl[24][i]);
     hijkl[26][i] = 1./3./(1.-e2)*(hijkl[5][i] -e*hijkl[6][i] -ees1113[i] +e*ffr1113[i] +3.*(ce-d)*hijkl[24][i]);
     hijkl[33][i] = 1./3./(1.-e2)*(2.*hijkl[5][i] +e*ees2013[i] -ffr2013[i] +3.*(de-c)*hijkl[31][i]);
     hijkl[34][i] = 1./3./(1.-e2)*(2.*hijkl[6][i] -ees0213[i] +e*ffr0213[i] +3.*(ce-d)*hijkl[32][i]);
     ees2003[i] = e203_s1[i] +e203_s2[i];
     ffr0203[i] = f203_r2[i] +f203_r1[i];
     ees0203[i] = h2*e203_s2[i] +2.*hm*e103_s2[i] +m2*e003_s2[i] +s12*e003_s1[i];
     ffr2003[i] = f2*f203_r2[i] +2.*fg*f103_r2[i] +g2*f003_r2[i] +r12*f003_r1[i];
     hijkl[37][i] = 1./3./(1.-e2)*(2.*hijkl[1][i] +e*ees2003[i] -ffr2003[i] +3.*(de-c)*hijkl[29][i]);
     hijkl[38][i] = 1./3./(1.-e2)*(2.*hijkl[2][i] -ees0203[i] +e*ffr0203[i] +3.*(ce-d)*hijkl[30][i]);
     ees3013[i] = e313_s1[i] +e313_s2[i];
     ffr0313[i] = f313_r2[i] +f313_r1[i];
     ffr3013[i] = f3*f313_r2[i] +3.*f2g*f213_r2[i] +3.*fg2*f113_r2[i] +g3*f013_r2[i] +r13*f013_r1[i];
     ees0313[i] = h3*e313_s2[i] +3.*h2m*e213_s2[i] +3.*hm2*e113_s2[i] +m3*e013_s2[i] +s13*e013_s1[i];
     hijkl[44][i] = (1./3./(1.-e2-c2-d2+2.*cde))*((1.-e2)*(hijkl[13][i] -pow(yv[i],1.)*d303[i]) -3.*(c-de)*hijkl[8][i] +(d-ce)*ees3013[i] +(c-de)*ffr3013[i]);
     hijkl[45][i] = (1./3./(1.-e2-c2-d2+2.*cde))*((1.-e2)*(hijkl[14][i] -pow(yv[i],1.)*d033[i]) -3.*(d-ce)*hijkl[9][i] +(d-ce)*ees0313[i] +(c-de)*ffr0313[i]);
     ees1123[i] = h*e223_s2[i] +m*e123_s2[i] +s1*e123_s1[i];
     ffr1123[i] = f*f223_r2[i] +g*f123_r2[i] +r1*f123_r1[i];
     hijkl[46][i] = 1./3./(1.-e2)*(-e*3.*hijkl[8][i] -ees3013[i] +e*ffr3013[i] +3.*(ce-d)*hijkl[44][i]);
     hijkl[47][i] = 1./3./(1.-e2)*(-e*3.*hijkl[9][i] +e*ees0313[i] -ffr0313[i] +3.*(de-c)*hijkl[45][i]);
     h1135[i] = (1./3./(1.-e2-c2-d2+2.*cde))*((1.-e2)*(2.*hijkl[12][i] -pow(yv[i],2.)*d113[i]) -(c-de)*h0123[i] -(d-ce)*h1023[i] +(d-ce)*ees1123[i] +(c-de)*ffr1123[i]);
     ees2113[i] = h*e313_s2[i] +m*e213_s2[i] +s1*e213_s1[i];
     ffr1213[i] = f*f313_r2[i] +g*f213_r2[i] +r1*f213_r1[i];
     ffr2113[i] = f2*f313_r2[i] +2.*fg*f213_r2[i] +g2*f113_r2[i] +r12*f113_r1[i];
     ees1213[i] = h2*e313_s2[i] +2.*hm*e213_s2[i] +m2*e113_s2[i] +s12*e113_s1[i];
     hijkl[40][i] = (1./3./(1.-e2-c2-d2+2.*cde))*((1.-e2)*(hijkl[10][i] -pow(yv[i],1.)*d213[i]) -2.*(c-de)*hijkl[12][i] -(d-ce)*hijkl[8][i] +(d-ce)*ees2113[i] +(c-de)*ffr2113[i]);
     hijkl[41][i] = (1./3./(1.-e2-c2-d2+2.*cde))*((1.-e2)*(hijkl[11][i] -pow(yv[i],1.)*d123[i]) -(c-de)*hijkl[9][i] -2.*(d-ce)*hijkl[12][i] +(d-ce)*ees1213[i] +(c-de)*ffr1213[i]);
     h2215_1[i] = 1./3./(1.-e2)*(hijkl[9][i] -e*2.*hijkl[12][i] +e*ees1213[i] -ffr1213[i] +3.*(de-c)*hijkl[41][i]);
     h2215_2[i] = 1./3./(1.-e2)*(hijkl[8][i] -e*2.*hijkl[12][i] +e*ffr2113[i] -ees2113[i] +3.*(ce-d)*hijkl[40][i]);
     hijkl[39][i]= 0.5*h2215_1[i] +0.5*h2215_2[i];
     hijkl[42][i] = 1./3./(1.-e2)*(3.*hijkl[8][i] +e*ees3013[i] -ffr3013[i] +3.*(de-c)*hijkl[44][i]);
     hijkl[43][i] = 1./3./(1.-e2)*(3.*hijkl[9][i] -ees0313[i] +e*ffr0313[i] +3.*(ce-d)*hijkl[45][i]);
   }

   // Evaluating the integrals.
   for (int i = 0; i < i_num_integrals; i++){
     o_sch[i] = dot_product(hijkl[i], signv, n_limits);
   }
  return o_sch;
}

double *node_force_quad_triangle(double *i_sch, double i_vec_int[][3],
                                 double i_a, double i_b, double i_c,
                                 double i_d, double i_e, double i_f,
                                 double i_factor,
                                 double *o_force){
  // Calculates the force on a node.
  double f[11];
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

void compute_forces_quad_triangle(double *i_sch, double i_vec_int[][3],
                                     double i_r1, double i_r2,
                                     double i_s1, double i_s2,
                                     double i_one_o_dr_sq, double i_one_o_ds_sq,
                                     double i_one_o_dr, double i_one_o_ds,
                                     double i_one_o_dr_one_o_ds,
                                     double i_factor,
                                     double *o_nodal_force[n_nodes], double *o_total_force){
  // Calculating nodal forces
  double a, b, c, d, e, f;
  // x3
  a = 2.*i_one_o_dr_sq;
  b = -4.*i_r1*i_one_o_dr_sq-3.*i_one_o_dr-4.*i_s1*i_one_o_dr_one_o_ds;
  c = 2.*i_one_o_ds_sq;
  d = -4.*i_s1*i_one_o_ds_sq-3.*i_one_o_ds-4.*i_r1*i_one_o_dr_one_o_ds;
  e = 4.*i_one_o_dr_one_o_ds;
  f = 1.+4.*i_r1*i_s1*i_one_o_dr_one_o_ds+3.*i_r1*i_one_o_dr+3.*i_s1*i_one_o_ds+2.*i_r1*i_r1*i_one_o_dr_sq+2.*i_s1*i_s1*i_one_o_ds_sq;
  node_force_quad_triangle(i_sch, i_vec_int, a, b, c, d, e, f, i_factor, o_nodal_force[0]);
  // x4
  a = 2.;
  b = -3.*i_r1 - i_r2;
  //c = d = e = 0.;
  f = i_r1*i_r1 + i_r1*i_r2;
  node_force_quad_triangle(i_sch, i_vec_int, a, b, 0., 0., 0., f, i_factor*i_one_o_dr_sq, o_nodal_force[1]);
  // x5
  //a = b = e = 0.;
  b = 0.;
  c = 2.;
  d = -3.*i_s1 - i_s2;
  f = i_s1*i_s1 + i_s1*i_s2;
  node_force_quad_triangle(i_sch, i_vec_int, 0., 0., c, d, 0., f, i_factor*i_one_o_ds_sq, o_nodal_force[2]);
  // x6
  a = -4.*i_one_o_dr_sq;
  b = 4.*i_one_o_dr+8.*i_r1*i_one_o_dr_sq+4.*i_s1*i_one_o_dr_one_o_ds;
  //c = 0.;
  d = 4.*i_r1*i_one_o_dr_one_o_ds;
  e = -4.*i_one_o_dr_one_o_ds;
  f = -4.*i_r1*i_one_o_dr-4.*i_r1*i_r1*i_one_o_dr_sq-4.*i_r1*i_s1*i_one_o_dr_one_o_ds;
  node_force_quad_triangle(i_sch, i_vec_int, a, b, 0., d, e, f, i_factor, o_nodal_force[3]);
  // x7
  //a = 0.;
  b = 4.*i_s1*i_one_o_dr_one_o_ds;
  c = -4.*i_one_o_ds_sq;
  d = 4.*i_one_o_ds+8.*i_s1*i_one_o_ds_sq+4.*i_r1*i_one_o_dr_one_o_ds;
  e = -4.*i_one_o_dr_one_o_ds;
  f = -4.*i_s1*i_one_o_ds-4.*i_s1*i_s1*i_one_o_ds_sq-4.*i_r1*i_s1*i_one_o_dr_one_o_ds;
  node_force_quad_triangle(i_sch, i_vec_int, 0., b, c, d, e, f, i_factor, o_nodal_force[4]);
  // x8
  // a = c = 0.;
  b = -4.*i_s1;
  d = -4.*i_r1;
  e = 4.;
  f = 4.*i_r1*i_s1;
  node_force_quad_triangle(i_sch, i_vec_int, 0., b, 0., d, e, f, i_factor*i_one_o_dr_one_o_ds, o_nodal_force[5]);
  //[8] [5] [9] [6] [12] [0]

  for (int i = 0; i < n_nodes; i++){
    o_total_force[0] += o_nodal_force[i][0];
    o_total_force[1] += o_nodal_force[i][1];
    o_total_force[2] += o_nodal_force[i][2];
  }
}

 void nodal_surface_force_quad_triangle(double *x1, double *x2, double *x3, double *x4, double *x5, double *b, double *p, double *q, double *n, double p_dot_q, double p_dot_q_sq, double mu, double nu, double a, double a_sq, double one_m_nu, double factor, double *nodal_force[n_nodes], double *total_force){
   // Characteristic vectors.
   double t[3];
   // Basis vectors (unitary).
   double p_x_t[3], q_x_t[3];
   //  Auxiliary constants to reduce computational cost.
   //  p dot t, q dot t
   double p_dot_t, q_dot_t;
   // Limits of the distance vector from the plane to dislocation line segment.
   // r_lim[0][] = vector from x1 to x3, r_lim[1][] = vector from x1 to x4, r_lim[2][] = vector from x1 to x5, r_lim[3][] = vector from x2 to x3.
   double r_lim[4][3];
   // ds = (sp[2] - sp[1]), rs = (rp[2] - rp[1]), one_o_ds = 1/(sp[2] - sp[1]), one_o_rs = 1/(rp[2] - rp[1])
   double ds, dr, one_o_dr, one_o_ds, one_o_dr_sq, one_o_ds_sq, one_o_dr_one_o_ds, dr_o_ds, ds_o_dr;
   //r, s limits, the points where r and s become functions of each other are r[2] and s[2].
   double rp[3], sp[3];
   //  y, r, s coordinates.
   double y[n_limits], r[n_limits], s[n_limits];
   // r(s), s(r) coordinates.
   double rfs[n_limits], sfr[n_limits];
   // Vectors for the integrals.
   double vec_int[11][3];
   // Scalar value of the integrals.
   double sch[48];
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
   /*
   yp[0] = init_point(r_lim[0], n    , t, n    , 3);
   yp[1] = init_point(r_lim[3], n    , t, n    , 3);
   rp[0] = init_point(r_lim[0], q_x_t, p, q_x_t, 3);
   rp[1] = init_point(r_lim[1], q_x_t, p, q_x_t, 3);
   sp[0] = init_point(r_lim[0], p_x_t, q, p_x_t, 3);
   sp[1] = init_point(r_lim[2], p_x_t, q, p_x_t, 3);
   */
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
  one_o_dr = 1./dr;
  one_o_ds = 1./ds;
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

 void main_nodal_surface_force_quad_triangle(double *x1, double *x2, double *x3, double *x4, double *x5, double *b, double mu, double nu, double a, double *nodal_force[n_nodes], double *total_force){
   /*
     Forces
     nodal_force[0][] = F_x3[x, y, z], nodal_force[1][] = F_x4[x, y, z],
     nodal_force[2][] = F_x5[x, y, z], nodal_force[3][] = F_x6[x, y, z]
     total_force[x, y, z] = F_x3[x, y, z] + F_x4[x, y, z] + F_x5[x, y, z] + F_x6[x, y, z]
   */
   // Characteristic vectors.
   double p[3], q[3], n[3], t[3];
   // t dot q, p dot q, (p dot q)^2, alpha = sine(angle between p and q).
   double t_dot_n, p_dot_q, p_dot_q_sq, alpha;
   double a_sq, factor, one_m_nu;
   double rot_centre[3], rot_x1[3], rot_x2[3];
   double t_x_n[3], mag_t_x_n, p_total_force[3], *p_nodal_force[n_nodes];
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
   one_m_nu = 1.-nu;
   alpha    = sin(acos(p_dot_q));
   factor   = 0.25*alpha*mu/pi/one_m_nu;
   if(t_dot_n != 0.0){
     nodal_surface_force_quad_triangle(x1, x2, x3, x4, x5, b, p, q, n, p_dot_q, p_dot_q_sq, mu, nu, a, a_sq, one_m_nu, factor, nodal_force, total_force);
   }
   else{
     for (int i = 0; i < n_nodes; i++){
       p_nodal_force[i] = malloc(3*sizeof(double));
     }
     double angle = 0.01*pi/180;
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
   double x1[3], x2[3], x3[3], x4[3], x5[3], b[3], mu, nu, a, *nodal_force[n_nodes], total_force[3];
   for(int i = 0; i < n_nodes; i++){
     nodal_force[i] = malloc(3*sizeof(double));
   }
   FILE * ptr_file;
   ptr_file =fopen("input.txt", "r");
   if (ptr_file == NULL){
     printf("File does not exist.\n");
   }
   fscanf(ptr_file, "%lf %lf %lf", &x1[0], &x1[1], &x1[2] );
   fscanf(ptr_file, "%lf %lf %lf", &x2[0], &x2[1], &x2[2] );
   fscanf(ptr_file, "%lf %lf %lf", &x3[0], &x3[1], &x3[2] );
   fscanf(ptr_file, "%lf %lf %lf", &x4[0], &x4[1], &x4[2] );
   fscanf(ptr_file, "%lf %lf %lf", &x5[0], &x5[1], &x5[2] );
   fscanf(ptr_file, "%lf %lf %lf", &b[0], &b[1], &b[2] );
   fscanf(ptr_file, "%lf", &mu );
   fscanf(ptr_file, "%lf", &nu );
   fscanf(ptr_file, "%lf", &a );
   fclose(ptr_file);
   main_nodal_surface_force_quad_triangle(x1, x2, x3, x4, x5, b, mu, nu, a, nodal_force, total_force);
   for(int i=0;i<3;i++) { printf("%f\t", nodal_force[0][i]); }; printf("\n");
   for(int i=0;i<3;i++) { printf("%f\t", nodal_force[1][i]); }; printf("\n");
   for(int i=0;i<3;i++) { printf("%f\t", nodal_force[2][i]); }; printf("\n");
   for(int i=0;i<3;i++) { printf("%f\t", nodal_force[3][i]); }; printf("\n");
   for(int i=0;i<3;i++) { printf("%f\t", nodal_force[4][i]); }; printf("\n");
   for(int i=0;i<3;i++) { printf("%f\t", nodal_force[5][i]); }; printf("\n");
   for(int i = 0; i < n_nodes; i++){
     free(nodal_force[i]);
   }
   return 0;
 }
