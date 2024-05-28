const double eps = 1e-14;

__device__ double cuda_init_point(double *i_vec1, double *i_vec2,
                  double *i_vec3, double *i_vec4,
                  int     i_vec_size){
  // Initialises points on the surface element given four vectors.
  double result = 0.0;
  double denom = 0.0;
  denom = cuda_dot_product(i_vec3, i_vec4, i_vec_size);
  if (denom == 0.0){
    //printf("nodal_surface_force_linear_rectangle: init_point: division by zero, dot_product(i_vec3, i_vec4) = %2.14f\n", denom);
    asm("trap;");
  }
  result = cuda_dot_product(i_vec1, i_vec2, i_vec_size)/denom;
  return result;
}

__device__ void cuda_integral_vector(double *i_p, double *i_q, double *i_b, double *i_t, double *i_n,
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

__device__ double cuda_seed_single_integral(double a, double b){
  // Single integral seed.
  double result = 0.0;
  double arg;
  arg = a+b;
  if (arg <= 0.0){
    //printf("nodal_surface_force_linear_rectangle: seed_single_integral: log(%2.14f) = %2.14f\n", arg, log(arg));
    asm("trap;");
  }
  result = log(arg);
  return result;
}

__device__ double cuda_integral_type_1(double a, double b,
                       double c, double d,
                       double e){
  double result = 0.0;
  result = 0.5*(a*b + (c-d)*e);
  return result;
}

__device__ double cuda_numer_seed_double_integral(double a,
                                  double b, double c, double d,
                                  double e){
  // Returns the argument of the numerator of the seed integral for double integrals.
  double result = 0.0;
  result = a*(b - c + d) + e;
  return result;
}

__device__ double cuda_denom_seed_double_integral(double a, double b,
                                  double c, double d){
  // Returns the argument of the denominator of the seed integral for double integrals.
  double result = 0.0;
  result = a*b + c*d;
  if(result == 0.0){
    //printf("nodal_surface_force_linear_rectangle: denom_seed_integral_00m3: division by 0, result = %2.14f\n", result);
    asm("trap;");
  }
  return result;
}

__device__ double cuda_seed_double_integral(double a, double b, double c, double d, double e,
                            double f, double g, double h, double i){
  // Returns the seed integral for double integrals.
  double result      = 0.0;
  double numerator   = 0.0;
  double denominator = 0.0;
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

__device__ double cuda_integral_type_2(double a,
                       double b, double c){
  double result = 0.0;
  result = a + b*c;
  return result;
}

__device__ double cuda_integral_type_3(double a,
                       double b, double c,
                       double d, double e){
  double result = 0.0;
  result = a + b*c + d*e;
  return result;
}

__device__ double cuda_integral_type_4(double a,
                       double b, double c,
                       double d, double e,
                       double f, double g){
  double result = 0.0;
  result = a + b*c + d*e + f*g;
  return result;
}

__device__ double cuda_integral_type_5(double a,
                       double b, double c,
                       double d, double e,
                       double f, double g,
                       double h, double i){
  double result = 0.0;
  result = a + b*c + d*e + f*g + h*i;
  return result;
}

__device__ double cuda_integral_type_6(double a, double b,
                       double c, double d,
                       double e, double f, double g,
                       double h, double i,
                       double j, double k){
  double result = 0.0;
  result = a*b + c*d + (e+f)*g + h*i + j*k;
  return result;
}

__device__ double cuda_integral_type_7(double a, double b,
                       double c, double d,
                       double e, double f, double g,
                       double h, double i){
  double result = 0.0;
  result = a*b + c*d + (e+f)*g + h*i;
  return result;
}

__device__ double cuda_integral_type_8(double a, double b,
                       double c, double d,
                       double e, double f,
                       double g, double h){
  double result = 0.0;
  result = a*b + c*d + e*f + g*h;
  return result;
}

__device__ double cuda_integral_type_9(double a,
                       double b, double c,
                       double d, double e,
                       double f, double g,
                       double h, double i,
                       double j, double k){
  double result = 0.0;
  result = a + b*c + d*e + f*g + h*i + j*k;
  return result;
}

__device__ double cuda_integral_type_10(double a,
                        double b, double c,
                        double d, double e,
                        double f, double g,
                        double h, double i,
                        double j, double k){
  double result = 0.0;
  result = a + b*c + d*e + f*g + h*i + j*k;
  return result;
}

__device__ double cuda_integral_type_11(double a,
                        double b, double c,
                        double d, double e,
                        double f, double g,
                        double h, double i){
  double result = 0.0;
  result = a + b*c + d*e + f*g + h*i;
  return result;
}

__device__ double cuda_integral_type_12(double a, double b,
                        double c, double d,
                        double e, double f,
                        double g, double h,
                        double i, double j,
                        double k, double l, double m){
  double result = 0.0;
  result = a*b + c*d + e*f + g*h + i*j + k*l*m;
  return result;
}

__device__ double *cuda_integrals_linear_rectangle(double *i_r, double *i_s, double *i_y,
                                   double *i_p, double *i_q, double *i_t,
                                   double i_p_dot_t, double i_q_dot_t, double i_a_sq,
                                   double *o_sch){
  // Sign vector for quick evaluation of integrals via the dot product.
  double signv[8] = {1.0,-1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0};
  const double third = 1.0/3.0;
  const double two_third = 2.0/3.0;
  // Integrals.
  double // Single integrals.
         a0m1[n_limits], b0m1[n_limits], c0m1[n_limits], a1m1[n_limits], b1m1[n_limits], c1m1[n_limits],
         a01 [n_limits], b01 [n_limits], c01 [n_limits], a11 [n_limits], b11 [n_limits],
         // Double integrals.
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
 double y_sq[n_limits], r_sq[n_limits], s_sq[n_limits],
        y_p_dot_t[n_limits], y_q_dot_t[n_limits], r_p_dot_t[n_limits], s_q_dot_t[n_limits],
        r_dot_p[n_limits], r_dot_q[n_limits], r_dot_t[n_limits],
        r_dot_p_sq[n_limits], r_dot_q_sq[n_limits], r_dot_t_sq[n_limits];
  /*
   ra = sqrt((r_vec dot r_vec) + a**2) from non-singular dislocation theory, there are 8 r_vec, thus 8 ra's.
   ra_sq = ra^2 (element-wise squaring)
   ra_c_o_3 = 1/3 * ra^3 (element-wise cubing and division)
  */
  double ra[n_limits], ra_sq[n_limits], ra_c_o_3[n_limits];
  // Vector R from x1 to x3, x4, x5, x6 and x2 to x3, x4, x5, x6. 8 vectors with 3 components each.
  double r_vec[n_limits][3];
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
  double                         two_p_dot_t,
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
    // Double seed integrals.
    d00m3[i] = cuda_seed_double_integral(1.0, ra[i], r_dot_p[i], r_dot_q[i], 0.0,
                                    1.0             , i_a_sq, one_m_p_dot_t_sq_m_q_dot_t_sq, y_sq[i]); // checked
    e00m3[i] = cuda_seed_double_integral(one_m_p_dot_t, ra[i], i_r[i], i_y[i], s_q_dot_t[i],
                                    one_m_p_dot_t_sq, i_a_sq, one_m_p_dot_t_sq_m_q_dot_t_sq, s_sq[i]); // checked
    f00m3[i] = cuda_seed_double_integral(one_m_q_dot_t, ra[i], i_s[i], i_y[i], r_p_dot_t[i],
                                    one_m_q_dot_t_sq, i_a_sq, one_m_p_dot_t_sq_m_q_dot_t_sq, r_sq[i]); // checked
    // Double integrals.
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

__device__ double *cuda_vertex_force_linear_rectangle(double *i_sch, double i_vec_int[][3],
                                      double i_r, double i_s, double i_rs,
                                      double i_factor,
                                      double *o_force){
  // Calculate the force on a vertex.
  double f[11];
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

__device__ void cuda_compute_forces_linear_rectangle(double *i_sch, double i_vec_int[][3],
                                     double *i_rp, double *i_sp, double i_factor,
                                     double o_nodal_force[][3], double *o_total_force){
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
__constant__ double d_mu, d_nu, d_a, d_a_sq, d_one_m_nu, d_factor, d_eps;

__global__ void se_cuda_nodal_surface_force_linear_rectangle(double *g_dln_arr, double *g_se_arr, double *g_b_arr, double *g_fx_arr, double *g_ftot_arr, int n_se, int n_dln){
  double x1[3], x2[3], x3[3], x4[3], x5[3], x6[3], b[3];
  // Characteristic vectors.
  double t[3], p[3], q[3], n[3];
  // Basis vectors (unitary).
  double p_x_t[3], q_x_t[3];
  // Limits of the distance vector from the plane to dislocation line segment.
  // r_lim[0][] = vector from x1 to x3, r_lim[1][] = vector from x2 to x6.
  double r_lim[2][3];
  //r, s limits
  double rp[2], sp[2];
  //  y, r, s coordinates
  double y[n_limits], r[n_limits], s[n_limits];
  // Vectors for the integrals.
  double vec_int[11][3];
  // Scalar value of the integrals.
  double sch[38];
  //  Auxiliary constants to reduce computational cost.
  //  p dot t, q dot t
  double p_dot_t, q_dot_t;
  // Auxiliary variables for nodal force calculation.
  double l_factor;
  double p_norm, q_norm;
  double nodal_force[n_nodes][3], total_force[3];
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
      if (abs(cuda_dot_product(t, n, 3)) > d_eps){
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
      }
    }
  }
}

__global__ void dln_cuda_nodal_surface_force_linear_rectangle(double *g_dln_arr, double *g_se_arr, double *g_b_arr, double *g_fx_arr, double *g_ftot_arr, int n_se, int n_dln){
  double x1[3], x2[3], x3[3], x4[3], x5[3], x6[3], b[3];
  // Characteristic vectors.
  double t[3], p[3], q[3], n[3];
  // Basis vectors (unitary).
  double p_x_t[3], q_x_t[3];
  // Limits of the distance vector from the plane to dislocation line segment.
  // r_lim[0][] = vector from x1 to x3, r_lim[1][] = vector from x2 to x6.
  double r_lim[2][3];
  //r, s limits
  double rp[2], sp[2];
  //  y, r, s coordinates
  double y[n_limits], r[n_limits], s[n_limits];
  // Vectors for the integrals.
  double vec_int[11][3];
  // Scalar value of the integrals.
  double sch[38];
  //  Auxiliary constants to reduce computational cost.
  //  p dot t, q dot t
  double p_dot_t, q_dot_t;
  // Auxiliary variables for nodal force calculation.
  double l_factor;
  double p_norm, q_norm;
  double nodal_force[n_nodes][3], total_force[3];
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
      if (abs(cuda_dot_product(t, n, 3)) > d_eps){
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
        dln_add_force_thread_device(nodal_force, total_force, g_fx_arr, g_ftot_arr, n_se, n_nodes, j);
      }
    }
  }
}

void main_se_cuda_nodal_surface_force_linear_rectangle(int n_se, int n_dln, int threads_per_block, char **argv){
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
  int blocks_per_grid;
  int idx1, idx2;
  cudaSetDevice(0);
  // Memory allocation
  b_arr[0] = (double *) malloc(3 * n_dln * sizeof(double));
  dln_node_arr[0] = (double *) malloc(3 * n_dln * sizeof(double));
  dln_node_arr[1] = (double *) malloc(3 * n_dln * sizeof(double));
  for (int i = 0; i < n_nodes; i++){
    se_node_arr[i] = (double *) malloc(3 * n_se  * sizeof(double));
    fx[i] = (double *) malloc(3 * sizeof(double));
    fx_arr[i] = (double *) malloc(3 * n_se * sizeof(double));
  }
  // Read input.
  FILE * ptr_file;
  ptr_file = fopen(argv[1], "r");
  if (ptr_file == NULL){
    printf("File does not exist.\n");
  }
  // Skip two lines.
  fscanf(ptr_file, "%*[^\n]\n");
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
  #ifdef debug2
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
  #endif
  #ifdef debug3
    clock_t begin = clock();
  #endif
  #ifdef debug2
    cudaEventRecord(start);
  #endif
  // Auxiliary constants.
  a_sq     = a*a;
  one_m_nu = 1.-nu;
  factor   = 0.25*mu/pi/one_m_nu;
  // Forward data map.
  x_b_arr   = b_host_device_map(b_arr[0], n_dln);
  x_dln_arr = dln_host_device_map(dln_node_arr[0], dln_node_arr[1], n_dln);
  x_se_arr  = element_host_device_map(se_node_arr, n_se, n_nodes);
  // Allocate device memory.
  checkCudaErrors( cudaMalloc( (void **) &d_x_dln_arr, 3 * n_dln * 2       * sizeof(double) ) );
  checkCudaErrors( cudaMalloc( (void **) &d_x_b_arr  , 3 * n_dln           * sizeof(double) ) );
  checkCudaErrors( cudaMalloc( (void **) &d_x_se_arr , 3 * n_se  * n_nodes * sizeof(double) ) );
  checkCudaErrors( cudaMalloc( (void **) &d_fx_arr   , 3 * n_se  * n_nodes * sizeof(double) ) );
  checkCudaErrors( cudaMalloc( (void **) &d_ftot_arr , 3 * n_se            * sizeof(double) ) );
  // Copy host to device.
  // Only pass node coordinates to device.
  checkCudaErrors( cudaMemcpyAsync(d_x_se_arr , x_se_arr , 3*n_se *n_nodes*sizeof(double), cudaMemcpyHostToDevice) );
  checkCudaErrors( cudaMemcpyAsync(d_x_dln_arr, x_dln_arr, 3*n_dln*2      *sizeof(double), cudaMemcpyHostToDevice) );
  checkCudaErrors( cudaMemcpyAsync(d_x_b_arr  , x_b_arr, 3*n_dln        *sizeof(double), cudaMemcpyHostToDevice) );
  // Initialising force arrays to zero.
  checkCudaErrors( cudaMemsetAsync(d_fx_arr  , 0.0, 3*n_se *n_nodes*sizeof(double)) );
  checkCudaErrors( cudaMemsetAsync(d_ftot_arr, 0.0, 3*n_se         *sizeof(double)) );
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_mu      , &mu      , sizeof(mu)) );
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_nu      , &nu      , sizeof(nu)) );
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_a       , &a       , sizeof(a)) );
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_a_sq    , &a_sq    , sizeof(a_sq)) );
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_one_m_nu, &one_m_nu, sizeof(one_m_nu)) );
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_factor  , &factor  , sizeof(factor)) );
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_eps     , &eps     , sizeof(eps)) );
  blocks_per_grid = (n_se + threads_per_block - 1)/threads_per_block;
  se_cuda_nodal_surface_force_linear_rectangle<<<blocks_per_grid, threads_per_block>>>(d_x_dln_arr, d_x_se_arr, d_x_b_arr, d_fx_arr, d_ftot_arr, n_se, n_dln);
  #ifdef debug2
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    printf("Time spent in SE parallelisation: %f (ms) \n", milliseconds);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
  #endif
  // Host code is executed asynchronously from the kernel execution.
  // Free all 1D arrays used to copy data to device.
  free(x_se_arr); free(x_dln_arr); free(x_b_arr);
  // Special case, where dislocation line is parallel with surface element.
  // Initialise forces.
  ftot_arr = (double *) malloc(3 * n_se * sizeof(double));
  for (int i = 0; i < n_nodes; i++){
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
      if (abs(dot_product(t, n, 3)) > eps){
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
  // Copy device to host
  // Only pass forces to host.
  x_fx_arr   = (double *) malloc(3 * n_se * n_nodes * sizeof(double));
  x_ftot_arr = (double *) malloc(3 * n_se           * sizeof(double));
  checkCudaErrors( cudaMemcpy(x_fx_arr  , d_fx_arr  , 3 * n_se * n_nodes * sizeof(double), cudaMemcpyDeviceToHost) );
  add_fx_device_host_map(x_fx_arr, fx_arr, n_se, n_nodes);
  free(x_fx_arr);
  checkCudaErrors( cudaMemcpy(x_ftot_arr, d_ftot_arr, 3 * n_se           * sizeof(double), cudaMemcpyDeviceToHost) );
  add_ftot_device_host_map(x_ftot_arr, ftot_arr, n_se);
  free(x_ftot_arr);
  #ifdef debug
    for (int i = 0; i < n_se; i++){
      printf("ftot_arr[%d] = [%2.14f, %2.14f, %2.14f]\n", i, ftot_arr[3*i], ftot_arr[3*i+1], ftot_arr[3*i+2]);
    }
    //for(int i = 0; i < n_nodes; i++){
    //  for(int j = 0; j < 3*n_se; j+=3){
    //    printf("fx_arr[%d] = %2.18f %2.18f %2.18f\n", i, fx_arr[i][j], fx_arr[i][j+1], fx_arr[i][j+2]);
    //  }
    //}
  #endif
  #ifdef debug3
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time elapsed = %f ms\n", time_spent*1000);
  #endif
  // CUDA exit.
  cudaDeviceReset();
  for(int i = 0; i < n_nodes; i++){free(fx[i]);}
  // Don't free these when using MATLAB. It silently crashes the program. Free them when using straight C.
  free(b_arr[0]); free(dln_node_arr[0]); free(dln_node_arr[1]); //free(x3_arr); free(x4_arr); free(x5_arr); free(x6_arr);
  for(int i=0; i < n_nodes; i++){free(se_node_arr[i]);}
  free(ftot_arr);
}

void main_dln_cuda_nodal_surface_force_linear_rectangle(int n_se, int n_dln, int threads_per_block, char **argv){
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
  int blocks_per_grid;
  //int n_se, n_dln;
  int idx1, idx2;
  cudaSetDevice(0);
  // Memory allocation
  b_arr[0] = (double *) malloc(3 * n_dln * sizeof(double));
  dln_node_arr[0] = (double *) malloc(3 * n_dln * sizeof(double));
  dln_node_arr[1] = (double *) malloc(3 * n_dln * sizeof(double));
  for (int i = 0; i < n_nodes; i++){
    se_node_arr[i] = (double *) malloc(3 * n_se  * sizeof(double));
    fx[i] = (double *) malloc(3 * sizeof(double));
    fx_arr[i] = (double *) malloc(3 * n_se * sizeof(double));
  }
  // Read input.
  FILE * ptr_file;
  ptr_file = fopen(argv[1], "r");
  if (ptr_file == NULL){
    printf("File does not exist.\n");
  }
  // Skip two lines.
  fscanf(ptr_file, "%*[^\n]\n");
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
  #ifdef debug2
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
  #endif
  #ifdef debug3
    clock_t begin = clock();
  #endif
  #ifdef debug2
    cudaEventRecord(start);
  #endif
  // Auxiliary constants.
  a_sq     = a*a;
  one_m_nu = 1.-nu;
  factor   = 0.25*mu/pi/one_m_nu;
  // Forward data map.
  x_b_arr   = element_host_device_map(b_arr, n_dln, 1);
  x_dln_arr = element_host_device_map(dln_node_arr, n_dln, 2);
  x_se_arr  = se_host_device_map(se_node_arr[0], se_node_arr[1], se_node_arr[2], se_node_arr[3], n_se);
  // Allocate device memory.
  checkCudaErrors( cudaMalloc( (void **) &d_x_dln_arr, 3 * n_dln * 2       * sizeof(double) ) );
  checkCudaErrors( cudaMalloc( (void **) &d_x_b_arr  , 3 * n_dln           * sizeof(double) ) );
  checkCudaErrors( cudaMalloc( (void **) &d_x_se_arr , 3 * n_se  * n_nodes * sizeof(double) ) );
  checkCudaErrors( cudaMalloc( (void **) &d_fx_arr   , 3 * n_se  * n_nodes * sizeof(double) ) );
  checkCudaErrors( cudaMalloc( (void **) &d_ftot_arr , 3 * n_se            * sizeof(double) ) );
  // Copy host to device.
  // Only pass node coordinates to device.
  checkCudaErrors( cudaMemcpyAsync(d_x_se_arr , x_se_arr , 3*n_se*n_nodes*sizeof(double), cudaMemcpyHostToDevice) );
  checkCudaErrors( cudaMemcpyAsync(d_x_dln_arr, x_dln_arr, 3*n_dln*2     *sizeof(double), cudaMemcpyHostToDevice) );
  checkCudaErrors( cudaMemcpyAsync(d_x_b_arr  , x_b_arr  , 3*n_dln       *sizeof(double), cudaMemcpyHostToDevice) );
  // Initialising force arrays to zero.
  checkCudaErrors( cudaMemsetAsync(d_fx_arr  , 0.0, 3*n_se*n_nodes*sizeof(double)) );
  checkCudaErrors( cudaMemsetAsync(d_ftot_arr, 0.0, 3*n_se        *sizeof(double)) );
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_mu      , &mu      , sizeof(mu)) );
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_nu      , &nu      , sizeof(nu)) );
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_a       , &a       , sizeof(a)) );
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_a_sq    , &a_sq    , sizeof(a_sq)) );
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_one_m_nu, &one_m_nu, sizeof(one_m_nu)) );
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_factor  , &factor  , sizeof(factor)) );
  checkCudaErrors( cudaMemcpyToSymbolAsync(d_eps     , &eps     , sizeof(eps)) );
  // Change this to be automated with CUDA's automatic estimation.
  blocks_per_grid = (n_dln + threads_per_block - 1)/threads_per_block;
  // CUDA
  dln_cuda_nodal_surface_force_linear_rectangle<<<blocks_per_grid, threads_per_block>>>(d_x_dln_arr, d_x_se_arr, d_x_b_arr, d_fx_arr, d_ftot_arr, n_se, n_dln);
  #ifdef debug2
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    printf("Time spent in DLN parallelisation: %f (ms) \n", milliseconds);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
  #endif

  // Host code is executed asynchronously from the kernel execution.
  // Free all 1D arrays used to copy data to device.
  free(x_se_arr); free(x_dln_arr); free(x_b_arr);
  // Special case, where dislocation line is parallel with surface element.
  // Initialise forces.
  ftot_arr = (double *) malloc(3 * n_se * sizeof(double));
  for (int i = 0; i < n_nodes; i++){
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
      if (abs(dot_product(t, n, 3)) > eps){
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
  x_fx_arr   = (double *) malloc(3 * n_se * n_nodes * sizeof(double));
  x_ftot_arr = (double *) malloc(3 * n_se * sizeof(double));
  // Synchronously copy forces from device to host.
  checkCudaErrors( cudaMemcpy(x_fx_arr, d_fx_arr, 3 * n_se * n_nodes * sizeof(double), cudaMemcpyDeviceToHost) );
  // Map 1D device array to 2D array for MATLAB.
  dln_add_fx_device_host_map(x_fx_arr, fx_arr, n_se, n_nodes);
  free(x_fx_arr);
  checkCudaErrors( cudaMemcpy(x_ftot_arr, d_ftot_arr, 3 * n_se * sizeof(double), cudaMemcpyDeviceToHost) );
  for (int i = 0; i < 3*n_se; i++){
    ftot_arr[i] += x_ftot_arr[i];
  }
  free(x_ftot_arr);

  #ifdef debug
    for (int i = 0; i < n_se; i++){
      printf("ftot_arr[%d] = [%2.14f, %2.14f, %2.14f]\n", i, ftot_arr[3*i], ftot_arr[3*i+1], ftot_arr[3*i+2]);
    }
    //for(int i = 0; i < n_nodes; i++){
    //  for(int j = 0; j < 3*n_se; j+=3){
    //    printf("fx_arr[%d] = %2.18f %2.18f %2.18f\n", i, fx_arr[i][j], fx_arr[i][j+1], fx_arr[i][j+2]);
    //  }
    //}
  #endif
  #ifdef debug3
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time elapsed = %f ms\n", time_spent*1000);
  #endif
  // CUDA exit.
  cudaDeviceReset();
  for(int i = 0; i < n_nodes; i++){free(fx[i]); free(se_node_arr[i]);}
  free(b_arr[0]); free(dln_node_arr[0]); free(dln_node_arr[1]); free(ftot_arr);
}
