/**************************************************************************
 *
 *      Module:  This module contains the functions needed for
 *               calculating interactions between dislocation
 *               segments.  See Tom Arsenlis for details on the
 *               method used to do the calculations.
 *
 *      Included functions:
 *               SegSegForceIntegrals()
 *               SpecialSegSegForce()
 *               SegSegForce()
 *               LocalSegForces()
 *               ComputeForces()
 *
 *************************************************************************/
//#include <mex.h>
//#include <cuda.h>
//#include <cuda_runtime.h>

#include <math.h>
#include <string.h>
#include <stdio.h>

/**************************************************************************
 * The following section defines the variables used later on:
 *
 *          p1*,p2*      endpoints for dislocation segment beginning
 *                       at point p1 and ending at point p2
 *          p3*,p4*      endpoints for dislocation segment beginning
 *                       at point p3 and ending at point p4
 *          bpx,bpy,bpz  burgers vector for segment p1->p2
 *          bx,by,bz     burgers vector for segment p3->p4
 *          a            core parameter
 *          MU           shear modulus
 *          NU           poisson ratio
 *          seg12Local   1 if either node of segment p1->p2 is local to
 *                       the current domain, zero otherwise.
 *          seg34Local   1 if either node of segment p3->p4 is local to
 *                       the current domain, zero otherwise.
 *          fp1*, fp2*,  pointers to locations in which to return forces
 *          fp3*, fp4*   on nodes p1 thru p4 respectively
 *
 *
 *          With regards to other variables, the following notation has
 *          used: t stands for the line vector, p stands for the normalised
 *          of a vector, c stands for a cross product operation, d stands 
 *          for a dot product operation. For example the variable,
 *          
 *          tctpctx is "t" "c"ross product with "t" "p"rime "c"ross product
 *                      with "t" "x" component
 *
 *************************************************************************/

__global__ void SegForceNBodyCUDA(double const *SoA,
                    double const a, double const MU, double const NU,
                    int const S,
                    double *f0x, double *f0y, double *f0z,
                    double *f1x, double *f1y, double *f1z);

__device__ void SegSegForce(double p1x, double p1y, double p1z,
                        double p2x, double p2y, double p2z,
                        double p3x, double p3y, double p3z,
                        double p4x, double p4y, double p4z,
                        double bpx, double bpy, double bpz,
                        double bx, double by, double bz,
                        double a, double MU, double NU,
                        int seg12Local, int seg34Local,
                        double *fp1x, double *fp1y, double *fp1z,
                        double *fp2x, double *fp2y, double *fp2z,
                        double *fp3x, double *fp3y, double *fp3z,
                        double *fp4x, double *fp4y, double *fp4z);

__device__ void SpecialSegSegForce(double p1x, double p1y, double p1z,
                               double p2x, double p2y, double p2z,
                               double p3x, double p3y, double p3z,
                               double p4x, double p4y, double p4z,
                               double bpx, double bpy, double bpz,
                               double bx, double by, double bz,
                               double a, double MU, double NU, double ecrit,
                               int seg12Local, int seg34Local,
                               double *fp1x, double *fp1y, double *fp1z,
                               double *fp2x, double *fp2y, double *fp2z,
                               double *fp3x, double *fp3y, double *fp3z,
                               double *fp4x, double *fp4y, double *fp4z);

__global__ void SegForceNBodyCUDA(double const *SoA,
                    double const a, double const MU, double const NU,
                    int const S,
                    double *f0x, double *f0y, double *f0z,
                    double *f1x, double *f1y, double *f1z)
{
    int j;
    double p1x, p1y, p1z;
    double p2x, p2y, p2z;
    double p3x, p3y, p3z;
    double p4x, p4y, p4z;
    double b1x, b1y, b1z;
    double b2x, b2y, b2z;
    int seg12Local, seg34Local;
    double fp1x, fp1y, fp1z;
    double fp2x, fp2y, fp2z;
    double fp3x, fp3y, fp3z;
    double fp4x, fp4y, fp4z;
        
        seg12Local=1;

        //pure N-body formulation. only calculated fp1 and fp2, and call function N^2 times.
        
        int i = threadIdx.x + blockIdx.x * blockDim.x;

        if (i < S) {
            
            b1x = SoA[9*i];
            b1y = SoA[9*i+1];
            b1z = SoA[9*i+2];
            p1x = SoA[9*i+3];
            p1y = SoA[9*i+4];
            p1z = SoA[9*i+5];
            p2x = SoA[9*i+6];
            p2y = SoA[9*i+7];
            p2z = SoA[9*i+8];

            for (j=0; j<S; j++) { //use blocks of B_vec, P1_vec, P2_vec in shared memory for efficient data reuse?

                    b2x = SoA[9*j];
                    b2y = SoA[9*j+1];
                    b2z = SoA[9*j+2];
                    p3x = SoA[9*j+3];
                    p3y = SoA[9*j+4];
                    p3z = SoA[9*j+5];
                    p4x = SoA[9*j+6];
                    p4y = SoA[9*j+7];
                    p4z = SoA[9*j+8];
                    
                    SegSegForce(p1x, p1y, p1z,
                                p2x, p2y, p2z,
                                p3x, p3y, p3z,
                                p4x, p4y, p4z,
                                b1x, b1y, b1z,
                                b2x, b2y, b2z,
                                a, MU, NU,
                                seg12Local, seg34Local,
                                &fp1x, &fp1y, &fp1z,
                                &fp2x, &fp2y, &fp2z,
                                &fp3x, &fp3y, &fp3z,
                                &fp4x, &fp4y, &fp4z);
                    
                    f0x[i] += fp1x;
                    f0y[i] += fp1y;
                    f0z[i] += fp1z;
                    
                    f1x[i] += fp2x;
                    f1y[i] += fp2y;
                    f1z[i] += fp2z;
            }
        }
    return;
}

/*
 *      SpecialSegSegForce() and SegSegForce() each reference the
 *      other function, so we need prototypes for these functions...
 */

/*-------------------------------------------------------------------------
 *
 *      Function:     SpecialSegSegForceIntegrals
 *      Description:  Calculates the integrals required for the
 *                    force calculation for near-parallel dislocation
 *                    segment pairs.
 *-----------------------------------------------------------------------*/
__device__ void SpecialSegSegForceIntegrals(double a2,double d2,double yin,double zin,
                                  double *f_003, double *f_103, double *f_013,
                                 double *f_113, double *f_213, double *f_123,
                                 double *f_005, double *f_105, double *f_015,
                                 double *f_115, double *f_215, double *f_125)
{
        double a2_d2, a2d2inv, ypz, ymz, Ra, Rainv, Log_Ra_ypz, common1;

        a2_d2 = a2 + d2;   
        a2d2inv = 1 / a2_d2;
        ypz = yin + zin;
        ymz = yin - zin;
        Ra = sqrt(a2_d2 + ypz*ypz);
        Rainv = 1 / Ra;
        Log_Ra_ypz = log(Ra + ypz);
        
        common1 = ymz * Ra * a2d2inv;

        *f_003 = Ra * a2d2inv;
        *f_103 = -0.5 * (Log_Ra_ypz - common1);
        *f_013 = -0.5 * (Log_Ra_ypz + common1);
        *f_113 = -Log_Ra_ypz;
        *f_213 = zin*Log_Ra_ypz - Ra;
        *f_123 = yin*Log_Ra_ypz - Ra;
        
        *f_005 =  a2d2inv * (2*a2d2inv*Ra - Rainv);
        *f_105 =  a2d2inv * (common1 - yin*Rainv);
        *f_015 = -a2d2inv * (common1 + zin*Rainv);
        *f_115 = -a2d2inv * ypz * Rainv;
        *f_215 =  Rainv - zin * *f_115;
        *f_125 =  Rainv - yin * *f_115;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     SegSegForceIntegrals
 *      Description:  Calculates the integrals required for the
 *                    force calculation
 *
 *-----------------------------------------------------------------------*/
__device__ void SegSegForceIntegrals(double a, double d, double c,double yin,double zin,
                        double *f_003, double *f_103, double *f_013,
                        double *f_113, double *f_203, double *f_023,
                        double *f_005, double *f_105, double *f_015,
                        double *f_115, double *f_205, double *f_025,
                        double *f_215, double *f_125, double *f_225,
                        double *f_305, double *f_035, double *f_315,
                        double *f_135)
{
        double c2, onemc2, onemc2inv, denom;
        double a2_d2, y2, z2, Ra;
        double Ra_Rdot_t, log_Ra_Rdot_t, zlog_Ra_Rdot_t;
        double Ra_Rdot_tp, log_Ra_Rdot_tp, ylog_Ra_Rdot_tp;
        double Rainv, Ra2_R_tinv, zRa2_R_tinv, z2Ra2_R_tinv;
        double Ra2_R_tpinv, yRa2_R_tpinv, y2Ra2_R_tpinv;
        double adf_003, tf_113;
        double commonf025, commonf035, commonf205;
        double commonf223, commonf225, commonf305;
        double ycommonf025, zcommonf305, zcommonf205;

        c2 = c*c;
        onemc2 = 1-c2; 
        onemc2inv = 1/onemc2;         
        a2_d2 = a*a+d*d*onemc2;
        y2    = yin*yin;
        z2    = zin*zin;
        Ra    = sqrt(a2_d2 + y2 + z2 + 2*yin*zin*c);       
        Rainv = 1/Ra;
        
        Ra_Rdot_tp = Ra+zin+yin*c;       
        Ra_Rdot_t  = Ra+yin+zin*c;       
         
         log_Ra_Rdot_tp =     log(Ra_Rdot_tp); 
        ylog_Ra_Rdot_tp = yin*log_Ra_Rdot_tp;
        
         log_Ra_Rdot_t =     log(Ra_Rdot_t); 
        zlog_Ra_Rdot_t = zin*log_Ra_Rdot_t; 
        
          Ra2_R_tpinv = Rainv/Ra_Rdot_tp; 
         yRa2_R_tpinv = yin*  Ra2_R_tpinv;
        y2Ra2_R_tpinv = yin* yRa2_R_tpinv;
        
          Ra2_R_tinv = Rainv/Ra_Rdot_t; 
         zRa2_R_tinv = zin* Ra2_R_tinv;
        z2Ra2_R_tinv = zin*zRa2_R_tinv;
        
        denom = 1/sqrt(onemc2*a2_d2);
        *f_003 = -2*denom*atan((1+c)*(Ra+yin+zin)*denom);
        
        adf_003 = a2_d2**f_003;
        commonf223 = (c*Ra - adf_003)*onemc2inv;

        *f_103 = (c*log_Ra_Rdot_t  - log_Ra_Rdot_tp)*onemc2inv; 
        *f_013 = (c*log_Ra_Rdot_tp - log_Ra_Rdot_t )*onemc2inv;
        *f_113 = (c*adf_003 - Ra)*onemc2inv;
        *f_203 =  zlog_Ra_Rdot_t  + commonf223;
        *f_023 =  ylog_Ra_Rdot_tp + commonf223;

         commonf225 = *f_003 - c*Rainv;
         commonf025 = c*yRa2_R_tpinv - Rainv;
        ycommonf025 = yin*commonf025;
         commonf205 = c*zRa2_R_tinv  - Rainv;
        zcommonf205 = zin*commonf205;
         commonf305 = log_Ra_Rdot_t  -(yin-c*zin)*Rainv - c2*z2Ra2_R_tinv;
        zcommonf305 = zin*commonf305;
         commonf035 = log_Ra_Rdot_tp -(zin-c*yin)*Rainv - c2*y2Ra2_R_tpinv;

        tf_113 = 2**f_113;
        
        *f_005 = (*f_003 - yRa2_R_tpinv - zRa2_R_tinv)/(a2_d2);
        *f_105 = (Ra2_R_tpinv - c*Ra2_R_tinv)*onemc2inv;
        *f_015 = (Ra2_R_tinv  - c*Ra2_R_tpinv)*onemc2inv;
        *f_115 = (Rainv - c*(yRa2_R_tpinv + zRa2_R_tinv + *f_003))*onemc2inv;
        *f_205 = (yRa2_R_tpinv + c2*zRa2_R_tinv  + commonf225)*onemc2inv;
        *f_025 = (zRa2_R_tinv  + c2*yRa2_R_tpinv + commonf225)*onemc2inv;
        *f_215 = (*f_013 - ycommonf025 + c*(zcommonf205-*f_103))*onemc2inv; 
        *f_125 = (*f_103 - zcommonf205 + c*(ycommonf025 - *f_013))*onemc2inv; 
        *f_225 = (*f_203 - zcommonf305 + c*(y2*commonf025 - tf_113))*onemc2inv;
        *f_305 = (y2Ra2_R_tpinv + c*commonf305 + 2**f_103)*onemc2inv;
        *f_035 = (z2Ra2_R_tinv  + c*commonf035 + 2**f_013)*onemc2inv;
        *f_315 = (tf_113 - y2*commonf025 + c*(zcommonf305 - *f_203))*onemc2inv;
        *f_135 = (tf_113 - z2*commonf205 + c*(yin*commonf035-*f_023))*onemc2inv;
       
        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     SpecialSegSegForce
 *      Description:  Special function for calculating forces between
 *                    dislocation segments too close to parallel to be
 *                    calculated via the function used for regular
 *                    segment/segment forces.
 *      Arguments:
 *          p1*,p2*      endpoints for dislocation segment beginning
 *                       at point p1 and ending at point p2
 *          p3*,p4*      endpoints for dislocation segment beginning
 *                       at point p3 and ending at point p4
 *          bpx,bpy,bpz  burgers vector for segment p1->p2
 *          bx,by,bz     burgers vector for segment p3->p4
 *          a            core parameter
 *          MU           shear modulus
 *          NU           poisson ratio
 *          seg12Local   1 if either node of segment p1->p2 is local to
 *                       the current domain, zero otherwise.
 *          seg34Local   1 if either node of segment p3->p4 is local to
 *                       the current domain, zero otherwise.
 *          fp1*, fp2*,  pointers to locations in which to return forces
 *          fp3*, fp4*   on nodes p1 thru p4 respectively
 *
 *
 *Refer to Appendix A.2. in "Enabling Strain Hardening Simulations with
 *Dislocation Dynamics" by A. Arsenlis et al. for mathematical treatment
 *-----------------------------------------------------------------------*/
__device__ void SpecialSegSegForce(double p1x, double p1y, double p1z,
                                double p2x, double p2y, double p2z,
                                double p3x, double p3y, double p3z,
                                double p4x, double p4y, double p4z,
                                double bpx, double bpy, double bpz,
                                double bx, double by, double bz,
                                double a, double MU, double NU, double ecrit,
                                int seg12Local, int seg34Local,
                                double *fp1x, double *fp1y, double *fp1z,
                                double *fp2x, double *fp2y, double *fp2z,
                                double *fp3x, double *fp3y, double *fp3z,
                                double *fp4x, double *fp4y, double *fp4z)
{
        double eps, c, a2, d2, flip;
        double Rx, Ry, Rz, Rdtp;
        //double Rdt;
        double oneoverL, oneoverLp;
        double temp, tempx, tempy, tempz, tempx2, tempy2, tempz2;
        double common1, common4;
        double common2x, common2y, common2z;
        double common3x, common3y, common3z;
        //double p1modx, p1mody, p1modz;
        //double p2modx, p2mody, p2modz;
        double p3modx, p3mody, p3modz;
        double p4modx, p4mody, p4modz;
        double vec1x, vec1y, vec1z;
        double tx, ty, tz;
        double tpx, tpy, tpz;
        double diffx, diffy, diffz, magdiff;
        double ndx, ndy, ndz;
        double wx, wy, wz;
        double qx, qy, qz;
        double ya, yb, za, zb;
        double fp1xcor, fp1ycor, fp1zcor;
        double fp2xcor, fp2ycor, fp2zcor;
        //double fp3xcor, fp3ycor, fp3zcor;
        //double fp4xcor, fp4ycor, fp4zcor;
        double f_003,  f_103,  f_013,  f_113,  f_213,  f_123,  f_005,  f_105;
        double f_003a, f_103a, f_013a, f_113a, f_213a, f_123a, f_005a, f_105a;
        double f_015,  f_115,  f_215,  f_125;
        double f_015a, f_115a, f_215a, f_125a;
        double Fint_003, Fint_113, Fint_005, Fint_115;
        double I_003x, I_003y, I_003z;
        double I_113x, I_113y, I_113z;
        double I_005x, I_005y, I_005z;
        double I_115x, I_115y, I_115z;
        double m4p, m8p, m4pn, a2m4pn, a2m8p;
        //double tdb, tdbp, nddb;
        //double bctx, bcty, bctz;
        //double bpctx, bpcty, bpctz;
        //double ndctx, ndcty, ndctz;
        //double bpctdb, double bpctdnd;
        //double bpctctx, bpctcty, bpctctz;
        double tpdb, tpdbp, nddbp;
        double bctpx, bctpy, bctpz;
        double bpctpx, bpctpy, bpctpz;
        double ndctpx, ndctpy, ndctpz;
        double bctpdbp, bctpdnd;
        double bctpctpx, bctpctpy, bctpctpz;
        double diffMag2;
        //double p1modMag2, p2modMag2;
        double p3modMag2, p4modMag2;
        double cotanthetac;
        double pivalue=3.141592653589793;


        cotanthetac = sqrt((1 - ecrit*1.01) / (ecrit*1.01));
        
        /* bunch of pre-factors that appear in equations to avoid
         * rewriting them all the time */
        eps    = 1e-16;
        a2     = a*a;
        m4p    = 0.25 * MU / pivalue;
        m8p    = 0.5 * m4p;
        m4pn   = m4p / ( 1 - NU ); /* factor mu/[8*pi*(1-nu)] */
        a2m4pn = a2 * m4pn;
        a2m8p  = a2 * m8p;

        /* setting pointers to zero */
        *fp1x = 0.0;
        *fp1y = 0.0;
        *fp1z = 0.0;
            
        *fp2x = 0.0;
        *fp2y = 0.0;
        *fp2z = 0.0;
            
        *fp3x = 0.0;
        *fp3y = 0.0;
        *fp3z = 0.0;
            
        *fp4x = 0.0;
        *fp4y = 0.0;
        *fp4z = 0.0;
        
        /* calculating line vector x,y,z-components for segment 4-3*/
        vec1x = p4x - p3x;
        vec1y = p4y - p3y;
        vec1z = p4z - p3z;
        
        /* line vector magnitude */
        oneoverL = 1/sqrt(vec1x*vec1x + vec1y*vec1y + vec1z*vec1z);
        
        /* unit line vector x,y,z-components for segment 4-3*/
        tx = vec1x*oneoverL;
        ty = vec1y*oneoverL;
        tz = vec1z*oneoverL;
        
        /* doing same thing, except for segment 2-1 */
        vec1x = p2x - p1x;
        vec1y = p2y - p1y;
        vec1z = p2z - p1z;
            
        oneoverLp = 1/sqrt(vec1x*vec1x + vec1y*vec1y + vec1z*vec1z);
        
        /* unit line vector x,y,z-components for segment 2-1*/
        tpx = vec1x*oneoverLp;
        tpy = vec1y*oneoverLp;
        tpz = vec1z*oneoverLp;
        
        /*dot product of unit line vectors of segments, equivalent to cosine
         *of theta angle between the two vectors*/
        c = tx*tpx + ty*tpy + tz*tpz; 
             
        /*just making sure notation is consistent and you don't get negative
         *angles*/
        flip = 0;
        if (c < 0) {
            flip = 1;
            tempx = p2x;
            tempy = p2y;
            tempz = p2z;
            p2x = p1x;
            p2y = p1y;
            p2z = p1z;
            p1x = tempx;
            p1y = tempy;
            p1z = tempz;
            tpx = -tpx;
            tpy = -tpy;
            tpz = -tpz;
            bpx = -bpx;
            bpy = -bpy;
            bpz = -bpz;
        } 
             
/*
 *      Find f1 and f2, but only if at least one of the endpoints
 *      is local to the domain.
 */
        if (seg12Local) {
            temp = (p4x-p3x)*tpx + (p4y-p3y)*tpy + (p4z-p3z)*tpz;
             
            p4modx = p3x + temp*tpx;
            p4mody = p3y + temp*tpy;
            p4modz = p3z + temp*tpz;
             
            diffx = p4x - p4modx;
            diffy = p4y - p4mody;
            diffz = p4z - p4modz;
             
            magdiff = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

            tempx2 = (0.5 * cotanthetac) * magdiff * tpx;
            tempy2 = (0.5 * cotanthetac) * magdiff * tpy;
            tempz2 = (0.5 * cotanthetac) * magdiff * tpz;

            p3modx = p3x + 0.5*diffx + tempx2;
            p3mody = p3y + 0.5*diffy + tempy2;
            p3modz = p3z + 0.5*diffz + tempz2;
             
            p4modx = p4modx + 0.5*diffx - tempx2;
            p4mody = p4mody + 0.5*diffy - tempy2;
            p4modz = p4modz + 0.5*diffz - tempz2;
             
            Rx = (p3modx - p1x);
            Ry = (p3mody - p1y);
            Rz = (p3modz - p1z);
             
            Rdtp = Rx*tpx + Ry*tpy + Rz*tpz;
             
            ndx = Rx - Rdtp*tpx;
            ndy = Ry - Rdtp*tpy;
            ndz = Rz - Rdtp*tpz;
             
            d2 = ndx*ndx + ndy*ndy + ndz*ndz;
             
            yb = p4modx*tpx + p4mody*tpy + p4modz*tpz;
            ya = p3modx*tpx + p3mody*tpy + p3modz*tpz;
            za = -(p1x*tpx + p1y*tpy + p1z*tpz);
            zb = -(p2x*tpx + p2y*tpy + p2z*tpz);
             
            SpecialSegSegForceIntegrals(a2, d2, ya, za,
                                        &f_003, &f_103, &f_013, &f_113,
                                        &f_213, &f_123, &f_005, &f_105,
                                        &f_015, &f_115, &f_215, &f_125);
        
            SpecialSegSegForceIntegrals(a2, d2, ya, zb,
                                        &f_003a, &f_103a, &f_013a, &f_113a,
                                        &f_213a, &f_123a, &f_005a, &f_105a,
                                        &f_015a, &f_115a, &f_215a, &f_125a);
        
            f_003 = f_003 - f_003a;
            f_103 = f_103 - f_103a;
            f_013 = f_013 - f_013a;
            f_113 = f_113 - f_113a;
            f_213 = f_213 - f_213a;
            f_123 = f_123 - f_123a;
            f_005 = f_005 - f_005a;
            f_105 = f_105 - f_105a;
            f_015 = f_015 - f_015a;
            f_115 = f_115 - f_115a;
            f_215 = f_215 - f_215a;
            f_125 = f_125 - f_125a;
        
            SpecialSegSegForceIntegrals(a2, d2, yb, za,
                                        &f_003a, &f_103a, &f_013a, &f_113a,
                                        &f_213a, &f_123a, &f_005a, &f_105a,
                                        &f_015a, &f_115a, &f_215a, &f_125a);
        
            f_003 = f_003 - f_003a;
            f_103 = f_103 - f_103a;
            f_013 = f_013 - f_013a;
            f_113 = f_113 - f_113a;
            f_213 = f_213 - f_213a;
            f_123 = f_123 - f_123a;
            f_005 = f_005 - f_005a;
            f_105 = f_105 - f_105a;
            f_015 = f_015 - f_015a;
            f_115 = f_115 - f_115a;
            f_215 = f_215 - f_215a;
            f_125 = f_125 - f_125a;
        
            SpecialSegSegForceIntegrals(a2, d2, yb, zb,
                                        &f_003a, &f_103a, &f_013a, &f_113a,
                                        &f_213a, &f_123a, &f_005a, &f_105a,
                                        &f_015a, &f_115a, &f_215a, &f_125a);
        
            f_003 = f_003 + f_003a;
            f_103 = f_103 + f_103a;
            f_013 = f_013 + f_013a;
            f_113 = f_113 + f_113a;
            f_213 = f_213 + f_213a;
            f_123 = f_123 + f_123a;
            f_005 = f_005 + f_005a;
            f_105 = f_105 + f_105a;
            f_015 = f_015 + f_015a;
            f_115 = f_115 + f_115a;
            f_215 = f_215 + f_215a;
            f_125 = f_125 + f_125a;
             
            tpdb = tpx*bx + tpy*by + tpz*bz;
            tpdbp = tpx*bpx + tpy*bpy + tpz*bpz;
            nddbp = ndx*bpx + ndy*bpy + ndz*bpz;
        
            bctpx = by*tpz - bz*tpy; 
            bctpy = bz*tpx - bx*tpz; 
            bctpz = bx*tpy - by*tpx;
            
            bpctpx = bpy*tpz - bpz*tpy; 
            bpctpy = bpz*tpx - bpx*tpz; 
            bpctpz = bpx*tpy - bpy*tpx;
            

            ndctpx = ndy*tpz - ndz*tpy; 
            ndctpy = ndz*tpx - ndx*tpz; 
            ndctpz = ndx*tpy - ndy*tpx;
            
            bctpdbp = bctpx*bpx + bctpy*bpy + bctpz*bpz;
            bctpdnd = bctpx*ndx + bctpy*ndy + bctpz*ndz;
            
            bctpctpx = tpdb*tpx - bx;
            bctpctpy = tpdb*tpy - by;
            bctpctpz = tpdb*tpz - bz;
            
            common1 = tpdbp*tpdb;
            
            common2x = common1*ndx;
            common2y = common1*ndy;
            common2z = common1*ndz;
            
            common3x = bctpdnd*bpctpx;
            common3y = bctpdnd*bpctpy;
            common3z = bctpdnd*bpctpz;
            
            I_003x = m4pn*(nddbp*bctpctpx+bctpdbp*ndctpx-common3x) -
                     m4p*common2x; 
            I_003y = m4pn*(nddbp*bctpctpy+bctpdbp*ndctpy-common3y) -
                     m4p*common2y; 
            I_003z = m4pn*(nddbp*bctpctpz+bctpdbp*ndctpz-common3z) -
                     m4p*common2z; 
            
            common1 = (m4pn-m4p)*tpdbp;
            
            I_113x =  common1*bctpctpx;
            I_113y =  common1*bctpctpy;
            I_113z =  common1*bctpctpz;
            
            common1 = m4pn*bctpdnd*nddbp;
            
            I_005x = -a2m8p*common2x - a2m4pn*common3x - common1*ndctpx;
            I_005y = -a2m8p*common2y - a2m4pn*common3y - common1*ndctpy;
            I_005z = -a2m8p*common2z - a2m4pn*common3z - common1*ndctpz;
            
            common1 = a2m8p*tpdbp;
            common4 = m4pn*bctpdnd*tpdbp;
            
            I_115x = -common1*bctpctpx - common4*ndctpx;
            I_115y = -common1*bctpctpy - common4*ndctpy;
            I_115z = -common1*bctpctpz - common4*ndctpz;
            
            Fint_003 = f_013 - za*f_003;
            Fint_113 = f_123 - za*f_113;
            Fint_005 = f_015 - za*f_005;
            Fint_115 = f_125 - za*f_115;
             
            *fp2x = (I_003x*Fint_003 + I_113x*Fint_113 + I_005x*Fint_005 +
                     I_115x*Fint_115) * oneoverLp;
            *fp2y = (I_003y*Fint_003 + I_113y*Fint_113 + I_005y*Fint_005 +
                     I_115y*Fint_115) * oneoverLp;
            *fp2z = (I_003z*Fint_003 + I_113z*Fint_113 + I_005z*Fint_005 +
                     I_115z*Fint_115) * oneoverLp;
             
            Fint_003 = zb*f_003 - f_013;
            Fint_113 = zb*f_113 - f_123;
            Fint_005 = zb*f_005 - f_015;
            Fint_115 = zb*f_115 - f_125;
             
            *fp1x = (I_003x*Fint_003 + I_113x*Fint_113 + I_005x*Fint_005 +
                     I_115x*Fint_115) * oneoverLp;
            *fp1y = (I_003y*Fint_003 + I_113y*Fint_113 + I_005y*Fint_005 +
                     I_115y*Fint_115) * oneoverLp;
            *fp1z = (I_003z*Fint_003 + I_113z*Fint_113 + I_005z*Fint_005 +
                     I_115z*Fint_115) * oneoverLp;
             

            diffMag2 = (diffx*diffx + diffy*diffy + diffz*diffz);
            p3modMag2 = (p3modx*p3modx + p3mody*p3mody + p3modz*p3modz);
            p4modMag2 = (p4modx*p4modx + p4mody*p4mody + p4modz*p4modz);

            if (diffMag2 > (eps * (p3modMag2+p4modMag2))) {
        
                SegSegForce(p3x, p3y, p3z, p3modx, p3mody, p3modz,
                            p1x, p1y, p1z, p2x, p2y, p2z,
                            bx, by, bz, bpx, bpy, bpz, a, MU, NU,
                            seg12Local, seg34Local,
                            &wx, &wy, &wz, &qx, &qy, &qz,
                            &fp1xcor, &fp1ycor, &fp1zcor,
                            &fp2xcor, &fp2ycor, &fp2zcor);
        
                *fp1x = *fp1x + fp1xcor;
                *fp1y = *fp1y + fp1ycor;
                *fp1z = *fp1z + fp1zcor;
                *fp2x = *fp2x + fp2xcor;
                *fp2y = *fp2y + fp2ycor;
                *fp2z = *fp2z + fp2zcor;
        
                SegSegForce(p4modx, p4mody, p4modz, p4x, p4y, p4z,
                            p1x, p1y, p1z, p2x, p2y, p2z,
                            bx, by, bz, bpx, bpy, bpz, a, MU, NU,
                            seg12Local, seg34Local,
                            &wx, &wy, &wz, &qx, &qy, &qz,
                            &fp1xcor, &fp1ycor, &fp1zcor,
                            &fp2xcor, &fp2ycor, &fp2zcor);
        
                *fp1x = *fp1x + fp1xcor;
                *fp1y = *fp1y + fp1ycor;
                *fp1z = *fp1z + fp1zcor;
                *fp2x = *fp2x + fp2xcor;
                *fp2y = *fp2y + fp2ycor;
                *fp2z = *fp2z + fp2zcor;
            }
             
/*
 *          If we flipped points 1 and 2 earlier, we have to compensate
 *          again here, but all that really needs to be switched are the
 *          forces.
 */
            if (flip == 1) {
                tempx = *fp2x;
                tempy = *fp2y;
                tempz = *fp2z;
                *fp2x = *fp1x;
                *fp2y = *fp1y;
                *fp2z = *fp1z;
                *fp1x = tempx;
                *fp1y = tempy;
                *fp1z = tempz;
            }
        } /* if segment p1->p2 is local */

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       SegSegForce
 *      Description:    Used to calculate the interaction forces between
 *                      dislocation segments analytically.
 *
 *      Arguments:
 *              p1*,p2*      endpoints for first dislocation segment starting
 *                           at p1x,p1y,p1z and ending at p2x,p2y,p2z
 *              p3*,p4*      endpoints for seond dislocation segment starting
 *                           at p3x,p3y,p3z and ending at p4x,p4y,p4z
 *              bxp,byp,bzp  burgers vector for segment p1 to p2
 *              bx,by,bz     burgers vector for segment p3 to p4
 *              a            core parameter
 *              MU           shear modulus
 *              NU           poisson ratio
 *              seg12Local   1 if either node of segment p1->p2 is local to
 *                           the current domain, zero otherwise.
 *              seg34Local   1 if either node of segment p3->p4 is local to
 *                           the current domain, zero otherwise.
 *              fp1*,fp2*,   pointers to locations in which to return
 *              fp3*,fp4*    forces on nodes located at p1, p2, p3 and
 *                           p4 respectively
 *            
 *Appendix A.1. in "Enabling Strain Hardening Simulations with
 *Dislocation Dynamics" by A. Arsenlis et al. for mathematical treatment          
 *-----------------------------------------------------------------------*/
__device__ void SegSegForce(double p1x, double p1y, double p1z,
                        double p2x, double p2y, double p2z,
                        double p3x, double p3y, double p3z,
                        double p4x, double p4y, double p4z,
                        double bpx, double bpy, double bpz,
                        double bx, double by, double bz,
                        double a, double MU, double NU,
                        int seg12Local, int seg34Local,
                        double *fp1x, double *fp1y, double *fp1z,
                        double *fp2x, double *fp2y, double *fp2z,
                        double *fp3x, double *fp3y, double *fp3z,
                        double *fp4x, double *fp4y, double *fp4z)
{
        double eps, d, c, c2, onemc2, onemc2inv, oneoverL, oneoverLp;
        double temp1, temp2, temp1a, temp1b, temp2a, temp2b;
        double R1x, R1y, R1z, R2x, R2y, R2z;
        double a2, m4p, m4pd, m8p, m8pd, m4pn, m4pnd, m4pnd2, m4pnd3;
        double a2m4pnd, a2m8pd, a2m4pn, a2m8p;
        double vec1x, vec1y, vec1z;
        double tx, ty, tz;
        double tpx, tpy, tpz;
        double tctpx, tctpy, tctpz;
        double ya, yb, za, zb;
        double f_003a, f_103a, f_013a, f_113a, f_203a, f_023a, f_005a, f_105a;
        double f_003,  f_103,  f_013,  f_113,  f_203,  f_023,  f_005,  f_105;
        double f_015a, f_115a, f_205a, f_025a, f_215a, f_125a, f_225a, f_305a;
        double f_015,  f_115,  f_205,  f_025,  f_215,  f_125,  f_225,  f_305;
        double f_035a, f_315a, f_135a;
        double f_035,  f_315,  f_135;
        double Fint_003, Fint_005, Fint_013, Fint_015, Fint_025, Fint_103;
        double Fint_105, Fint_115, Fint_125, Fint_205, Fint_215;
        double I_003x, I_003y, I_003z, I_005x, I_005y, I_005z;
        double I_013x, I_013y, I_013z, I_015x, I_015y, I_015z;
        double I_025x, I_025y, I_025z, I_103x, I_103y, I_103z;
        double I_105x, I_105y, I_105z, I_115x, I_115y, I_115z;
        double I_125x, I_125y, I_125z, I_205x, I_205y, I_205z;
        double I_215x, I_215y, I_215z;
        double I00ax, I00ay, I00az, I01ax, I01ay, I01az;
        double I10ax, I10ay, I10az, I00bx, I00by, I00bz;
        double I01bx, I01by, I01bz;
        //double I10bx, I10by, I10bz;
        double bctctpx, bctctpy, bctctpz;
        double bctdbp;
        double bctx, bcty, bctz;
        //double bpctpctx, bpctpcty, bpctpctz;
        double bpctpdb;
        double bpctpx, bpctpy, bpctpz;
        //double tcbpctx, tcbpcty, tcbpctz;
        //double tcbpdb;
        //double tcbpdtp;
        //double tcbpx, tcbpy, tcbpz;
        //double tctpcbpctx, tctpcbpcty, tctpcbpctz;
        double tctpcbpdb;
        //double tctpcbpdtp;
        //double tctpcbpx, tctpcbpy, tctpcbpz;
        //double tctpctx, tctpcty, tctpctz;
        double tctpdb;
        double tdb, tdbp;
        double tpcbctpx, tpcbctpy, tpcbctpz;
        double tpcbdbp;
        double tpcbdt;
        //double tpcbx, tpcby, tpcbz;
        //double tpctcbctpx, tpctcbctpy, tpctcbctpz;
        double tpctcbdbp;
        double tpctcbdt;
        //double tpctcbx, tpctcby, tpctcbz;
        double tpctctpx, tpctctpy, tpctctpz;
        double tpctdbp;
        double tpctx, tpcty, tpctz;
        double tpdb, tpdbp;
        double pivalue=3.141592653589793;

        eps = 1e-6;            

        *fp1x = 0.0;
        *fp1y = 0.0;
        *fp1z = 0.0;

        *fp2x = 0.0;
        *fp2y = 0.0;
        *fp2z = 0.0;

        *fp3x = 0.0;
        *fp3y = 0.0;
        *fp3z = 0.0;

        *fp4x = 0.0;
        *fp4y = 0.0;
        *fp4z = 0.0;

        vec1x = p4x - p3x; /*x-component of line vector of segment 4-3*/
        vec1y = p4y - p3y; /*y-component of line vector of segment 4-3*/
        vec1z = p4z - p3z; /*z-component of line vector of segment 4-3*/

        /*inverse magnitude of vector for segment 4-3*/
        oneoverL = 1/sqrt(vec1x*vec1x+vec1y*vec1y+vec1z*vec1z); 
        
        /*unit line vector x,y,z-components for segment 4-3*/
        tx = vec1x*oneoverL;
        ty = vec1y*oneoverL;
        tz = vec1z*oneoverL;
        
        vec1x = p2x - p1x; /*x-component of line vector of segment 2-1*/
        vec1y = p2y - p1y; /*y-component of line vector of segment 2-1*/
        vec1z = p2z - p1z; /*z-component of line vector of segment 2-1*/

        /*inverse magnitude of vector for segment 2-1*/
        oneoverLp = 1/sqrt(vec1x*vec1x+vec1y*vec1y+vec1z*vec1z);

        /*unit line vector x,y,z-components for segment 2-1*/
        tpx = vec1x*oneoverLp;
        tpy = vec1y*oneoverLp;
        tpz = vec1z*oneoverLp;
        
        /*tctpx, tctpy, tctpz are the determinants for submatrices of 
         *cross-product of line vectors of segments 2-1 and 4-3:
         *det|i j k ; tx ty tz ; tpx tpy tpz|*/
        tctpx = ty*tpz - tz*tpy;
        tctpy = tz*tpx - tx*tpz;
        tctpz = tx*tpy - ty*tpx;
        
        /* c is the dot product of the unit line vectors of segments 
         * 2-1 and 4-3. c2 is the magnitude of dot product of the line*/
        c = tx*tpx + ty*tpy + tz*tpz;
        c2 = c*c;
        onemc2 = 1 - c2;
        
        /*In order to check whether the two segments are parallel or not
         *an "if" statemement is used comparing the eps, an arbitrarily
         *small number close to zero, to the onemc2 variable. If the latter
         *is larger than eps, it will treat it as non-parallel and run this
         *loop, otherwise (else statement found at the bottom of script) 
         *it will use "SpecialSegSegForces" function which is defined above*/
        if (onemc2 > eps) {

            onemc2inv = 1/onemc2;

            R1x = p3x - p1x;
            R1y = p3y - p1y;
            R1z = p3z - p1z;

            R2x = p4x - p2x;
            R2y = p4y - p2y;
            R2z = p4z - p2z;

            d = (R2x*tctpx + R2y*tctpy + R2z*tctpz) * onemc2inv;

            temp1a = R1x*tx + R1y*ty + R1z*tz;
            temp1b = R2x*tx + R2y*ty + R2z*tz;

            temp2a = R1x*tpx + R1y*tpy + R1z*tpz;
            temp2b = R2x*tpx + R2y*tpy + R2z*tpz;

            ya = (temp1a - c*temp2a) * onemc2inv;
            yb = (temp1b - c*temp2b) * onemc2inv;

            za = (temp2a - c*temp1a) * onemc2inv;
            zb = (temp2b - c*temp1b) * onemc2inv;

/*
 *          For this first call to SegSegForceIntegrals() we can
 *          just pass the addresses of f_nnn variables rather than use
 *          the f_nnna variables and then copy the values.
 */
            
            /* ~f(x3-x2)? */
            SegSegForceIntegrals(a, d, c, ya, za, &f_003, &f_103,
                                 &f_013, &f_113, &f_203, &f_023,
                                 &f_005, &f_105, &f_015, &f_115,
                                 &f_205, &f_025, &f_215, &f_125,
                                 &f_225, &f_305, &f_035, &f_315,
                                 &f_135);
            /* ~f(x3-x2)? */
            SegSegForceIntegrals(a, d, c, ya, zb, &f_003a, &f_103a,
                                 &f_013a, &f_113a, &f_203a, &f_023a,
                                 &f_005a, &f_105a, &f_015a, &f_115a,
                                 &f_205a, &f_025a, &f_215a, &f_125a,
                                 &f_225a, &f_305a, &f_035a, &f_315a,
                                 &f_135a);

            f_003 = f_003 - f_003a;
            f_103 = f_103 - f_103a;
            f_013 = f_013 - f_013a;
            f_113 = f_113 - f_113a;
            f_203 = f_203 - f_203a;
            f_023 = f_023 - f_023a;
            f_005 = f_005 - f_005a;
            f_105 = f_105 - f_105a;
            f_015 = f_015 - f_015a;
            f_115 = f_115 - f_115a;
            f_205 = f_205 - f_205a;
            f_025 = f_025 - f_025a;
            f_215 = f_215 - f_215a;
            f_125 = f_125 - f_125a;
            f_225 = f_225 - f_225a;
            f_305 = f_305 - f_305a;
            f_035 = f_035 - f_035a;
            f_315 = f_315 - f_315a;
            f_135 = f_135 - f_135a;        

            /* ~f(x4-x1)? */
            SegSegForceIntegrals(a, d, c, yb, za, &f_003a, &f_103a,
                                 &f_013a, &f_113a, &f_203a, &f_023a,
                                 &f_005a, &f_105a, &f_015a, &f_115a,
                                 &f_205a, &f_025a, &f_215a, &f_125a,
                                 &f_225a, &f_305a, &f_035a, &f_315a,
                                 &f_135a); 
            

            f_003 = f_003 - f_003a;
            f_103 = f_103 - f_103a;
            f_013 = f_013 - f_013a;
            f_113 = f_113 - f_113a;
            f_203 = f_203 - f_203a;
            f_023 = f_023 - f_023a;
            f_005 = f_005 - f_005a;
            f_105 = f_105 - f_105a;
            f_015 = f_015 - f_015a;
            f_115 = f_115 - f_115a;
            f_205 = f_205 - f_205a;
            f_025 = f_025 - f_025a;
            f_215 = f_215 - f_215a;
            f_125 = f_125 - f_125a;
            f_225 = f_225 - f_225a;
            f_305 = f_305 - f_305a;
            f_035 = f_035 - f_035a;
            f_315 = f_315 - f_315a;
            f_135 = f_135 - f_135a;

            /* ~f(x3-x1)? */
            SegSegForceIntegrals(a, d, c, yb, zb, &f_003a, &f_103a,
                                 &f_013a, &f_113a, &f_203a, &f_023a,
                                 &f_005a, &f_105a, &f_015a, &f_115a,
                                 &f_205a, &f_025a, &f_215a, &f_125a,
                                 &f_225a, &f_305a, &f_035a, &f_315a,
                                 &f_135a);

            f_003 = f_003 + f_003a;
            f_103 = f_103 + f_103a;
            f_013 = f_013 + f_013a;
            f_113 = f_113 + f_113a;
            f_203 = f_203 + f_203a;
            f_023 = f_023 + f_023a;
            f_005 = f_005 + f_005a;
            f_105 = f_105 + f_105a;
            f_015 = f_015 + f_015a;
            f_115 = f_115 + f_115a;
            f_205 = f_205 + f_205a;
            f_025 = f_025 + f_025a;
            f_215 = f_215 + f_215a;
            f_125 = f_125 + f_125a;
            f_225 = f_225 + f_225a;
            f_305 = f_305 + f_305a;
            f_035 = f_035 + f_035a;
            f_315 = f_315 + f_315a;
            f_135 = f_135 + f_135a;


            a2 = a*a;
            m4p = 0.25 * MU / pivalue;
            m4pd =  m4p * d;
            m8p = 0.5 * m4p;
            m8pd = m8p * d;
            m4pn = m4p / ( 1 - NU );
            m4pnd = m4pn * d;
            m4pnd2 = m4pnd * d;
            m4pnd3 = m4pnd2 * d;
            a2m4pnd = a2 * m4pnd;
            a2m8pd = a2 * m8pd;
            a2m4pn = a2 * m4pn;
            a2m8p = a2 * m8p;

            tpctx = -tctpx;
            tpcty = -tctpy;
            tpctz = -tctpz;

            //tcbpx = ty*bpz - tz*bpy;
            //tcbpy = tz*bpx - tx*bpz;
            //tcbpz = tx*bpy - ty*bpx;

            //tpcbx = tpy*bz - tpz*by;
            //tpcby = tpz*bx - tpx*bz;
            //tpcbz = tpx*by - tpy*bx;

            bctx = by*tz - bz*ty;
            bcty = bz*tx - bx*tz;
            bctz = bx*ty - by*tx;


            bpctpx = bpy*tpz - bpz*tpy;
            bpctpy = bpz*tpx - bpx*tpz;
            bpctpz = bpx*tpy - bpy*tpx;

            tdb = tx*bx + ty*by + tz*bz;
            tdbp = tx*bpx + ty*bpy + tz*bpz;
            tpdb = tpx*bx + tpy*by + tpz*bz;
            tpdbp = tpx*bpx + tpy*bpy + tpz*bpz;

            tctpdb =  tctpx*bx + tctpy*by + tctpz*bz;
            tpctdbp = tpctx*bpx + tpcty*bpy + tpctz*bpz;
            //tcbpdtp = tpctdbp;
            tpcbdt = tctpdb;

            bpctpdb = bpctpx*bx + bpctpy*by + bpctpz*bz;
            bctdbp = bctx*bpx + bcty*bpy + bctz*bpz;
            //tcbpdb = bctdbp;
            tpcbdbp = bpctpdb;

            //tctpctx = tpx - c*tx;
            //tctpcty = tpy - c*ty;
            //tctpctz = tpz - c*tz;


            tpctctpx = tx - c*tpx;
            tpctctpy = ty - c*tpy;
            tpctctpz = tz - c*tpz;

            //tctpcbpx = tdbp*tpx - tpdbp*tx;
            //tctpcbpy = tdbp*tpy - tpdbp*ty;
            //tctpcbpz = tdbp*tpz - tpdbp*tz;

            //tpctcbx = tpdb*tx - tdb*tpx;
            //tpctcby = tpdb*ty - tdb*tpy;
            //tpctcbz = tpdb*tz - tdb*tpz;

            //tcbpctx = bpx - tdbp*tx;
            //tcbpcty = bpy - tdbp*ty;
            //tcbpctz = bpz - tdbp*tz;

            tpcbctpx = bx - tpdb*tpx;
            tpcbctpy = by - tpdb*tpy;
            tpcbctpz = bz - tpdb*tpz;
   
            //bpctpctx = tdbp*tpx - c*bpx;
            //bpctpcty = tdbp*tpy - c*bpy;
            //bpctpctz = tdbp*tpz - c*bpz;

            bctctpx = tpdb*tx - c*bx;
            bctctpy = tpdb*ty - c*by;
            bctctpz = tpdb*tz - c*bz;

            //tctpcbpctx = tdbp*tpctx;
            //tctpcbpcty = tdbp*tpcty;
            //tctpcbpctz = tdbp*tpctz;

            //tpctcbctpx = tpdb*tctpx;
            //tpctcbctpy = tpdb*tctpy;
            //tpctcbctpz = tpdb*tctpz;

            //tctpcbpdtp = tdbp - tpdbp*c;
            tpctcbdt = tpdb - tdb*c;
            tctpcbpdb =  tdbp*tpdb - tpdbp*tdb;
            tpctcbdbp = tctpcbpdb;
/*
 *          Only calculate the forces for segment p1->p2 if at least one
 *          of the segment's nodes is local to the current domain.
 */
            if (seg12Local) {

                temp1 = tpdb*tdbp + tpctcbdbp;

                I00ax = temp1 * tctpx;
                I00ay = temp1 * tctpy;
                I00az = temp1 * tctpz;

                I00bx = bpctpx * tpctcbdt;
                I00by = bpctpy * tpctcbdt;
                I00bz = bpctpz * tpctcbdt;

                temp1 = m4pnd * tpctdbp;
                temp2 = m4pnd * bctdbp;

                I_003x = m4pd*I00ax - m4pnd*I00bx + temp1*bctctpx +
                         temp2*tpctctpx;
                I_003y = m4pd*I00ay - m4pnd*I00by + temp1*bctctpy +
                         temp2*tpctctpy;
                I_003z = m4pd*I00az - m4pnd*I00bz + temp1*bctctpz +
                         temp2*tpctctpz;

                temp1 = m4pnd3 * tpctcbdt * tpctdbp;

                I_005x = a2m8pd*I00ax - a2m4pnd*I00bx - temp1*tpctctpx; 
                I_005y = a2m8pd*I00ay - a2m4pnd*I00by - temp1*tpctctpy; 
                I_005z = a2m8pd*I00az - a2m4pnd*I00bz - temp1*tpctctpz; 

                I01ax = tpctx*tpcbdbp - tpcbctpx*tdbp;
                I01ay = tpcty*tpcbdbp - tpcbctpy*tdbp;
                I01az = tpctz*tpcbdbp - tpcbctpz*tdbp;

                I01bx = -bpctpx * tpcbdt;
                I01by = -bpctpy * tpcbdt;
                I01bz = -bpctpz * tpcbdt;

                temp1 = m4pn * tpdbp;

                I_013x = -temp1 * bctctpx + m4p*I01ax - m4pn*I01bx;
                I_013y = -temp1 * bctctpy + m4p*I01ay - m4pn*I01by;
                I_013z = -temp1 * bctctpz + m4p*I01az - m4pn*I01bz;

                temp1 = m4pnd2 * (tpcbdt*tpctdbp + tpctcbdt*tpdbp);

                I_015x = a2m8p*I01ax - a2m4pn*I01bx + temp1*tpctctpx;
                I_015y = a2m8p*I01ay - a2m4pn*I01by + temp1*tpctctpy;
                I_015z = a2m8p*I01az - a2m4pn*I01bz + temp1*tpctctpz;

                I10ax = bctctpx*tdbp - tpctx*bctdbp;
                I10ay = bctctpy*tdbp - tpcty*bctdbp;
                I10az = bctctpz*tdbp - tpctz*bctdbp;

                temp1 = m4pn * tdbp; 
                temp2 = m4pn * bctdbp;

                I_103x = m4p*I10ax - temp1*bctctpx + temp2*tpctx;
                I_103y = m4p*I10ay - temp1*bctctpy + temp2*tpcty;
                I_103z = m4p*I10az - temp1*bctctpz + temp2*tpctz;

                temp1 = m4pnd2 * tpctcbdt * tdbp;
                temp2 = m4pnd2 * tpctcbdt * tpctdbp;

                I_105x = a2m8p*I10ax + temp1*tpctctpx - temp2*tpctx;
                I_105y = a2m8p*I10ay + temp1*tpctctpy - temp2*tpcty;
                I_105z = a2m8p*I10az + temp1*tpctctpz - temp2*tpctz;

                temp1 = (m4pnd * tpcbdt * tpdbp);

                I_025x = -temp1 * tpctctpx;
                I_025y = -temp1 * tpctctpy;
                I_025z = -temp1 * tpctctpz;

                temp1 = (m4pnd * tpctcbdt * tdbp);

                I_205x = temp1 * tpctx;
                I_205y = temp1 * tpcty;
                I_205z = temp1 * tpctz;

                temp1 = m4pnd * (tpctcbdt*tpdbp + tpcbdt*tpctdbp);
                temp2 = m4pnd * tpcbdt * tdbp;

                I_115x = temp1*tpctx - temp2*tpctctpx;
                I_115y = temp1*tpcty - temp2*tpctctpy;
                I_115z = temp1*tpctz - temp2*tpctctpz;

                temp1 = (m4pn * tpcbdt * tpdbp);

                I_125x = -temp1 * tpctx;
                I_125y = -temp1 * tpcty;
                I_125z = -temp1 * tpctz;

                temp1 = (m4pn * tpcbdt * tdbp);

                I_215x = -temp1 * tpctx;
                I_215y = -temp1 * tpcty;
                I_215z = -temp1 * tpctz;

                Fint_003 = f_013 - zb*f_003;
                Fint_103 = f_113 - zb*f_103;
                Fint_013 = f_023 - zb*f_013;
                Fint_005 = f_015 - zb*f_005;
                Fint_105 = f_115 - zb*f_105;
                Fint_015 = f_025 - zb*f_015;
                Fint_115 = f_125 - zb*f_115;
                Fint_205 = f_215 - zb*f_205;
                Fint_025 = f_035 - zb*f_025;
                Fint_215 = f_225 - zb*f_215;
                Fint_125 = f_135 - zb*f_125;

                *fp1x = (I_003x*Fint_003 + I_103x*Fint_103 + I_013x*Fint_013 +
                         I_005x*Fint_005 + I_105x*Fint_105 + I_015x*Fint_015 +
                         I_115x*Fint_115 + I_205x*Fint_205 + I_025x*Fint_025 +
                         I_215x*Fint_215 + I_125x*Fint_125) * oneoverLp;

                *fp1y = (I_003y*Fint_003 + I_103y*Fint_103 + I_013y*Fint_013 +
                         I_005y*Fint_005 + I_105y*Fint_105 + I_015y*Fint_015 +
                         I_115y*Fint_115 + I_205y*Fint_205 + I_025y*Fint_025 +
                         I_215y*Fint_215 + I_125y*Fint_125) * oneoverLp;

                *fp1z = (I_003z*Fint_003 + I_103z*Fint_103 + I_013z*Fint_013 +
                         I_005z*Fint_005 + I_105z*Fint_105 + I_015z*Fint_015 +
                         I_115z*Fint_115 + I_205z*Fint_205 + I_025z*Fint_025 +
                         I_215z*Fint_215 + I_125z*Fint_125) * oneoverLp;
   
                Fint_003 = za*f_003 - f_013;
                Fint_103 = za*f_103 - f_113;
                Fint_013 = za*f_013 - f_023;
                Fint_005 = za*f_005 - f_015;
                Fint_105 = za*f_105 - f_115;
                Fint_015 = za*f_015 - f_025;
                Fint_115 = za*f_115 - f_125;
                Fint_205 = za*f_205 - f_215;
                Fint_025 = za*f_025 - f_035;
                Fint_215 = za*f_215 - f_225;
                Fint_125 = za*f_125 - f_135;

                *fp2x = (I_003x*Fint_003 + I_103x*Fint_103 + I_013x*Fint_013 +
                         I_005x*Fint_005 + I_105x*Fint_105 + I_015x*Fint_015 +
                         I_115x*Fint_115 + I_205x*Fint_205 + I_025x*Fint_025 +
                         I_215x*Fint_215 + I_125x*Fint_125) * oneoverLp;

                *fp2y = (I_003y*Fint_003 + I_103y*Fint_103 + I_013y*Fint_013 +
                         I_005y*Fint_005 + I_105y*Fint_105 + I_015y*Fint_015 +
                         I_115y*Fint_115 + I_205y*Fint_205 + I_025y*Fint_025 +
                         I_215y*Fint_215 + I_125y*Fint_125) * oneoverLp;

                *fp2z = (I_003z*Fint_003 + I_103z*Fint_103 + I_013z*Fint_013 +
                         I_005z*Fint_005 + I_105z*Fint_105 + I_015z*Fint_015 +
                         I_115z*Fint_115 + I_205z*Fint_205 + I_025z*Fint_025 +
                         I_215z*Fint_215 + I_125z*Fint_125) * oneoverLp;
   
            } /* if segment p1->p2 is "local" */

        } else {
/*
 *          The two lines are parallel, so we have to use a special
 *          lower dimensional function
 */
            SpecialSegSegForce(p1x, p1y, p1z, p2x, p2y, p2z,
                               p3x, p3y, p3z, p4x, p4y, p4z,
                               bpx, bpy, bpz, bx, by, bz, a, MU, NU,
                               eps, seg12Local, seg34Local,
                               fp1x, fp1y, fp1z, fp2x, fp2y, fp2z,
                               fp3x, fp3y, fp3z, fp4x, fp4y, fp4z);
       }

       return;
}


/* --------------------------------------------------------------*/
/* ------------ the MEX driver runs on the CPU ------------------*/
/* --------------------------------------------------------------*/

/*
void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) 
{
    //Get the total memory needed for the various arrays on the GPU (bytes)
    mwSize m = mxGetM(prhs[0]);
    mwSize n = mxGetN(prhs[0]);
    mwSize arraySize = m * n;
    mem_size = sizeof(double) * arraySize;
    //all the other arrays are the same size (except a,MU,NU,linkid,S)
    
    // Define variables
    double *P1x_vec, *P1y_vec, *P1z_vec;
    double *P2x_vec, *P2y_vec, *P2z_vec;
    double *Bx_vec, *By_vec, *Bz_vec;
    double a,MU,NU;
    int linkid,S;
    double *f0x, *f0y, *f0z, *f1x, *f1y, *f1z;

    // allocate the memory on the GPU
    // memory for input vectors
    HANDLE_ERROR ( cudaMalloc( &P1x_vec, mem_size ) );
    HANDLE_ERROR ( cudaMalloc( &P1y_vec, mem_size ) );
    HANDLE_ERROR ( cudaMalloc( &P1z_vec, mem_size ) );
    HANDLE_ERROR ( cudaMalloc( &P2x_vec, mem_size ) );
    HANDLE_ERROR ( cudaMalloc( &P2y_vec, mem_size ) );
    HANDLE_ERROR ( cudaMalloc( &P2z_vec, mem_size ) );
    HANDLE_ERROR ( cudaMalloc( &Bx_vec, mem_size ) );
    HANDLE_ERROR ( cudaMalloc( &By_vec, mem_size ) );
    HANDLE_ERROR ( cudaMalloc( &Bz_vec, mem_size ) );
    // memory for constants
    HANDLE_ERROR ( cudaMalloc( &a, sizeof(double) ) );
    HANDLE_ERROR ( cudaMalloc( &MU, sizeof(double) ) );
    HANDLE_ERROR ( cudaMalloc( &NU, sizeof(double) ) );
    HANDLE_ERROR ( cudaMalloc( &linkid, sizeof(int) ) );
    HANDLE_ERROR ( cudaMalloc( &S, sizeof(int) ) );
    // memory for output vectors
    HANDLE_ERROR ( cudaMalloc( &f0x, mem_size ) );
    HANDLE_ERROR ( cudaMalloc( &f0y, mem_size ) );
    HANDLE_ERROR ( cudaMalloc( &f0z, mem_size ) );
    HANDLE_ERROR ( cudaMalloc( &f1x, mem_size ) );
    HANDLE_ERROR ( cudaMalloc( &f1y, mem_size ) );
    HANDLE_ERROR ( cudaMalloc( &f1z, mem_size ) );
    
    // Copy input data across to the GPU
    // memory for input vectors
    HANDLE_ERROR ( cudaMemcpy( P1x_vec, (double*) mexGetData(prhs[0]), mem_size, cudaMemcpyHostToDevice ) );
    HANDLE_ERROR ( cudaMemcpy( P1y_vec, (double*) mexGetData(prhs[1]), mem_size, cudaMemcpyHostToDevice ) );
    HANDLE_ERROR ( cudaMemcpy( P1z_vec, (double*) mexGetData(prhs[2]), mem_size, cudaMemcpyHostToDevice ) );
    HANDLE_ERROR ( cudaMemcpy( P2x_vec, (double*) mexGetData(prhs[3]), mem_size, cudaMemcpyHostToDevice ) );
    HANDLE_ERROR ( cudaMemcpy( P2y_vec, (double*) mexGetData(prhs[4]), mem_size, cudaMemcpyHostToDevice ) );
    HANDLE_ERROR ( cudaMemcpy( P2z_vec, (double*) mexGetData(prhs[5]), mem_size, cudaMemcpyHostToDevice ) );
    HANDLE_ERROR ( cudaMemcpy( Bx_vec, (double*) mexGetData(prhs[6]), mem_size, cudaMemcpyHostToDevice ) );
    HANDLE_ERROR ( cudaMemcpy( By_vec, (double*) mexGetData(prhs[7]), mem_size, cudaMemcpyHostToDevice ) );
    HANDLE_ERROR ( cudaMemcpy( Bz_vec, (double*) mexGetData(prhs[8]), mem_size, cudaMemcpyHostToDevice ) );
    // memory for constants
    HANDLE_ERROR ( cudaMemcpy( a, (double) mexGetData(prhs[9]), sizeof(double), cudaMemcpyHostToDevice ) );
    HANDLE_ERROR ( cudaMemcpy( MU, (double) mexGetData(prhs[10]), sizeof(double), cudaMemcpyHostToDevice ) );
    HANDLE_ERROR ( cudaMemcpy( NU, (double) mexGetData(prhs[11]), sizeof(double), cudaMemcpyHostToDevice ) );
    HANDLE_ERROR ( cudaMemcpy( linkid, (int) mexGetData(prhs[12]), sizeof(int), cudaMemcpyHostToDevice ) );
    HANDLE_ERROR ( cudaMemcpy( S, (int) mexGetData(prhs[13]), sizeof(int), cudaMemcpyHostToDevice ) );
    
    // Define the block and grid size - make this dynamic at later stage !!!!
    int blockSize = 1024;
    dim3 block(blockSize);
    dim 3 grid(ceil(arraySize/(double)blockSize));
    
    // Use the CUDA runtime to run the kernel
    SegForceNBodyCUDA <<< grid, block >>> (P1x_vec, P1y_vec, P1z_vec,
                                           P2x_vec, P2y_vec, P2z_vec,
                                           Bx_vec, By_vec, Bz_vec,
                                           a, MU, NU,
                                           linkid, S,
                                           f0x, f0y, f0z,
                                           f1x, f1y, f1z);
                                           
    // Create the output arrays for MATLAB
    plhs[0]=mxCreateNumericMatrix(m,n,mxSINGLE_CLASS,mxREAL);
    plhs[1]=mxCreateNumericMatrix(m,n,mxSINGLE_CLASS,mxREAL);
    plhs[2]=mxCreateNumericMatrix(m,n,mxSINGLE_CLASS,mxREAL);
    plhs[3]=mxCreateNumericMatrix(m,n,mxSINGLE_CLASS,mxREAL);
    plhs[4]=mxCreateNumericMatrix(m,n,mxSINGLE_CLASS,mxREAL);
    plhs[5]=mxCreateNumericMatrix(m,n,mxSINGLE_CLASS,mxREAL);

    //Copy the data from the card into the MATLAB array
    cudaMemcpy( (double*)mxGetData(plhs[0]), f0x, mem_size, cudaMemcpyDeviceToHost );
    cudaMemcpy( (double*)mxGetData(plhs[1]), f0y, mem_size, cudaMemcpyDeviceToHost );
    cudaMemcpy( (double*)mxGetData(plhs[2]), f0z, mem_size, cudaMemcpyDeviceToHost );
    cudaMemcpy( (double*)mxGetData(plhs[3]), f1x, mem_size, cudaMemcpyDeviceToHost );
    cudaMemcpy( (double*)mxGetData(plhs[4]), f1y, mem_size, cudaMemcpyDeviceToHost );
    cudaMemcpy( (double*)mxGetData(plhs[5]), f1z, mem_size, cudaMemcpyDeviceToHost );

    //Free the data on the card
    cudaFree( P1x_vec );
    cudaFree( P1y_vec );
    cudaFree( P1z_vec );
    cudaFree( P2x_vec );
    cudaFree( P2y_vec );
    cudaFree( P2z_vec );
    cudaFree( Bx_vec );
    cudaFree( By_vec );
    cudaFree( Bz_vec );
    cudaFree( a );
    cudaFree( MU );
    cudaFree( NU );
    cudaFree( linkid );
    cudaFree( S );
}
*/

    
    
