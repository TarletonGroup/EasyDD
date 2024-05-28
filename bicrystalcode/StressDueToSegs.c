/**************************************************************************
 *
 *      Function:    StressDueToSeg
 *      Description: Calculate the stress at point p from the segment
 *                   starting at point p1 and ending at point p2.
 *
 *      Arguments:
 *         px, py, pz     coordinates of field point at which stress is to
 *                        be evaluated
 *         p1x, p1y, p1z  starting position of the dislocation segment
 *         p2x, p2y, p2z  ending position of the dislocation segment
 *         bx, by, bz     burgers vector associated with segment going
 *                        from p1 to p2
 *         a              core value
 *         MU             shear modulus
 *         NU             poisson ratio
 *         stress         array of stresses form the indicated segment
 *                        at the field point requested
 *                            [0] = stressxx
 *                            [1] = stressyy
 *                            [2] = stresszz
 *                            [3] = stressxy
 *                            [4] = stressyz
 *                            [5] = stressxz
 *
 *************************************************************************/
#define _USE_MATH_DEFINES
#include <math.h>
#include <mex.h>

static void StressDueToSeg(int N, int S, double const *px, double const *py, double const *pz,
                    double const *p1x, double const *p1y, double const *p1z,
                    double const *p2x, double const *p2y, double const *p2z,
                    double const *bx, double const *by, double const *bz,
                    double a, double MU, double NU,
                    double *sxx, double *syy, double *szz,
                    double *sxy, double *syz, double *sxz);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
    {
    int N,S;
    double const *px,*py,*pz;
    double const *p1x,*p1y,*p1z,*p2x,*p2y,*p2z;
    double const *bx,*by,*bz;
    double a,MU,NU;
    double *sxx,*syy,*szz,*sxy,*syz,*sxz;
    int i;
    
    /*get the scalar inputs N,S,a,MU,NU
    create pointers to the input matrices*/
    N=mxGetScalar(prhs[0]);
    S=mxGetScalar(prhs[1]);
    
    px=mxGetPr(prhs[2]);
    py=mxGetPr(prhs[3]);
    pz=mxGetPr(prhs[4]);
    
    p1x=mxGetPr(prhs[5]);
    p1y=mxGetPr(prhs[6]);
    p1z=mxGetPr(prhs[7]);
    
    p2x=mxGetPr(prhs[8]);
    p2y=mxGetPr(prhs[9]);
    p2z=mxGetPr(prhs[10]); 
    
    bx=mxGetPr(prhs[11]);
    by=mxGetPr(prhs[12]);
    bz=mxGetPr(prhs[13]);
    
    a=mxGetScalar(prhs[14]);
    MU=mxGetScalar(prhs[15]);
    NU=mxGetScalar(prhs[16]);
     
    /*call the output pointer to the output matrix*/
    plhs[0]=mxCreateDoubleMatrix(N,1,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(N,1,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(N,1,mxREAL);
    plhs[3]=mxCreateDoubleMatrix(N,1,mxREAL);
    plhs[4]=mxCreateDoubleMatrix(N,1,mxREAL);
    plhs[5]=mxCreateDoubleMatrix(N,1,mxREAL);

    /*create a C pointer to a copy of the output matrix*/
	sxx= mxGetPr(plhs[0]);
    syy= mxGetPr(plhs[1]);
    szz= mxGetPr(plhs[2]);
    sxy= mxGetPr(plhs[3]);
    syz= mxGetPr(plhs[4]);
    sxz= mxGetPr(plhs[5]);
    
    /*call the C subroutine*/
    StressDueToSeg(N,S, px, py, pz,p1x, p1y, p1z,
                p2x, p2y, p2z, bx, by, bz,
                a, MU, NU,
                sxx, syy, szz, sxy, syz, sxz);
    
    }

void StressDueToSeg(int N, int S, double const *px, double const *py, double const *pz,
                    double const *p1x, double const *p1y, double const *p1z,
                    double const *p2x, double const *p2y, double const *p2z,
                    double const *bx, double const *by, double const *bz,
                    double a, double MU, double NU,
                    double *sxx, double *syy, double *szz,
                    double *sxy, double *syz, double *sxz)
{
    double   oneoverLp, common;
    double   vec1x, vec1y, vec1z;
    double   tpx, tpy, tpz;
    double   Rx, Ry, Rz, Rdt;
    double   ndx, ndy, ndz;
    double   d2, s1, s2, a2, a2_d2, a2d2inv;
    double   Ra, Rainv, Ra3inv, sRa3inv;
    double   s_03a, s_13a, s_05a, s_15a, s_25a;
    double   s_03b, s_13b, s_05b, s_15b, s_25b;
    double   s_03, s_13, s_05, s_15, s_25;
    double   m4p, m8p, m4pn, mn4pn, a2m8p;
    double   txbx, txby, txbz;
    double   dxbx, dxby, dxbz;
    double   dxbdt, dmdxx, dmdyy, dmdzz, dmdxy, dmdyz, dmdxz;
    double   tmtxx, tmtyy, tmtzz, tmtxy, tmtyz, tmtxz;
    double   tmdxx, tmdyy, tmdzz, tmdxy, tmdyz, tmdxz;
    double   tmtxbxx, tmtxbyy, tmtxbzz, tmtxbxy, tmtxbyz, tmtxbxz;
    double   dmtxbxx, dmtxbyy, dmtxbzz, dmtxbxy, dmtxbyz, dmtxbxz;
    double   tmdxbxx, tmdxbyy, tmdxbzz, tmdxbxy, tmdxbyz, tmdxbxz;
    double   I_03xx, I_03yy, I_03zz, I_03xy, I_03yz, I_03xz;
    double   I_13xx, I_13yy, I_13zz, I_13xy, I_13yz, I_13xz;
    double   I_05xx, I_05yy, I_05zz, I_05xy, I_05yz, I_05xz;
    double   I_15xx, I_15yy, I_15zz, I_15xy, I_15yz, I_15xz;
    double   I_25xx, I_25yy, I_25zz, I_25xy, I_25yz, I_25xz;
    int     coord, seg;
    
    /*precompute some constants*/
    m4p = 0.25 * MU / M_PI;
    m8p = 0.5 * m4p;
    m4pn = m4p / (1 - NU);
    mn4pn = m4pn * NU;
    a2 = a * a;
    a2m8p = a2 * m8p;

    for (coord=0; coord<N; coord++) { /*loop over the grid points*/
        
        for (seg=0; seg<S; seg++) { /*loop over the segments*/
        
        vec1x = p2x[seg] - p1x[seg];
        vec1y = p2y[seg] - p1y[seg];
        vec1z = p2z[seg] - p1z[seg];

        oneoverLp = 1.0 / sqrt(vec1x*vec1x + vec1y*vec1y + vec1z*vec1z);

        tpx = vec1x * oneoverLp;
        tpy = vec1y * oneoverLp;
        tpz = vec1z * oneoverLp;

        Rx = px[coord] - p1x[seg];
        Ry = py[coord] - p1y[seg];
        Rz = pz[coord] - p1z[seg];

        Rdt = Rx*tpx + Ry*tpy + Rz*tpz;

        ndx = Rx - Rdt*tpx;
        ndy = Ry - Rdt*tpy;
        ndz = Rz - Rdt*tpz;

        d2 = ndx*ndx + ndy*ndy + ndz*ndz;

        s1 = -Rdt;
        s2 = -((px[coord]-p2x[seg])*tpx + (py[coord]-p2y[seg])*tpy + (pz[coord]-p2z[seg])*tpz);
        a2_d2 = a2 + d2;
        a2d2inv = 1.0 / a2_d2;

        Ra = sqrt(a2_d2 + s1*s1);
        Rainv = 1.0 / Ra;
        Ra3inv = Rainv * Rainv * Rainv;
        sRa3inv = s1 * Ra3inv;

        s_03a = s1 * Rainv * a2d2inv;
        s_13a = -Rainv;
        s_05a = (2.0*s_03a + sRa3inv) * a2d2inv;
        s_15a = -Ra3inv;
        s_25a = s_03a - sRa3inv;

        Ra = sqrt(a2_d2 + s2*s2);
        Rainv = 1.0 / Ra;
        Ra3inv = Rainv * Rainv * Rainv;
        sRa3inv = s2 * Ra3inv;

        s_03b = s2 * Rainv * a2d2inv;
        s_13b = -Rainv;
        s_05b = (2.0*s_03b + sRa3inv) * a2d2inv;
        s_15b = -Ra3inv;
        s_25b = s_03b - sRa3inv;

        s_03 = s_03b - s_03a;
        s_13 = s_13b - s_13a;
        s_05 = s_05b - s_05a;
        s_15 = s_15b - s_15a;
        s_25 = s_25b - s_25a;

        txbx = tpy*bz[seg] - tpz*by[seg];
        txby = tpz*bx[seg] - tpx*bz[seg];
        txbz = tpx*by[seg] - tpy*bx[seg];

        dxbx = ndy*bz[seg] - ndz*by[seg];
        dxby = ndz*bx[seg] - ndx*bz[seg];
        dxbz = ndx*by[seg] - ndy*bx[seg];

        dxbdt = dxbx*tpx + dxby*tpy + dxbz*tpz;

        dmdxx = ndx * ndx;
        dmdyy = ndy * ndy;
        dmdzz = ndz * ndz;
        dmdxy = ndx * ndy;
        dmdyz = ndy * ndz;
        dmdxz = ndx * ndz;

        tmtxx = tpx * tpx;
        tmtyy = tpy * tpy;
        tmtzz = tpz * tpz;
        tmtxy = tpx * tpy;
        tmtyz = tpy * tpz;
        tmtxz = tpx * tpz;

        tmdxx = 2.0 * tpx * ndx;
        tmdyy = 2.0 * tpy * ndy;
        tmdzz = 2.0 * tpz * ndz;
        tmdxy = tpx*ndy + tpy*ndx;
        tmdyz = tpy*ndz + tpz*ndy;
        tmdxz = tpx*ndz + tpz*ndx;

        tmtxbxx = 2.0 * tpx * txbx;
        tmtxbyy = 2.0 * tpy * txby;
        tmtxbzz = 2.0 * tpz * txbz;
        tmtxbxy = tpx*txby + tpy*txbx;
        tmtxbyz = tpy*txbz + tpz*txby;
        tmtxbxz = tpx*txbz + tpz*txbx;

        dmtxbxx = 2.0 * ndx * txbx;
        dmtxbyy = 2.0 * ndy * txby;
        dmtxbzz = 2.0 * ndz * txbz;
        dmtxbxy = ndx*txby + ndy*txbx;
        dmtxbyz = ndy*txbz + ndz*txby;
        dmtxbxz = ndx*txbz + ndz*txbx;

        tmdxbxx = 2.0 * tpx * dxbx;
        tmdxbyy = 2.0 * tpy * dxby;
        tmdxbzz = 2.0 * tpz * dxbz;
        tmdxbxy = tpx*dxby + tpy*dxbx;
        tmdxbyz = tpy*dxbz + tpz*dxby;
        tmdxbxz = tpx*dxbz + tpz*dxbx;

        common = m4pn * dxbdt;

        I_03xx = common + m4pn*dmtxbxx - m4p*tmdxbxx;
        I_03yy = common + m4pn*dmtxbyy - m4p*tmdxbyy;
        I_03zz = common + m4pn*dmtxbzz - m4p*tmdxbzz;
        I_03xy = m4pn*dmtxbxy - m4p*tmdxbxy;
        I_03yz = m4pn*dmtxbyz - m4p*tmdxbyz;
        I_03xz = m4pn*dmtxbxz - m4p*tmdxbxz;

        I_13xx = -mn4pn * tmtxbxx;
        I_13yy = -mn4pn * tmtxbyy;
        I_13zz = -mn4pn * tmtxbzz;
        I_13xy = -mn4pn * tmtxbxy;
        I_13yz = -mn4pn * tmtxbyz;
        I_13xz = -mn4pn * tmtxbxz;

        I_05xx = common*(a2+dmdxx) - a2m8p*tmdxbxx;
        I_05yy = common*(a2+dmdyy) - a2m8p*tmdxbyy;
        I_05zz = common*(a2+dmdzz) - a2m8p*tmdxbzz;
        I_05xy = common*dmdxy - a2m8p*tmdxbxy;
        I_05yz = common*dmdyz - a2m8p*tmdxbyz;
        I_05xz = common*dmdxz - a2m8p*tmdxbxz;

        I_15xx = a2m8p*tmtxbxx - common*tmdxx;
        I_15yy = a2m8p*tmtxbyy - common*tmdyy;
        I_15zz = a2m8p*tmtxbzz - common*tmdzz;
        I_15xy = a2m8p*tmtxbxy - common*tmdxy;
        I_15yz = a2m8p*tmtxbyz - common*tmdyz;
        I_15xz = a2m8p*tmtxbxz - common*tmdxz;

        I_25xx = common * tmtxx;
        I_25yy = common * tmtyy;
        I_25zz = common * tmtzz;
        I_25xy = common * tmtxy;
        I_25yz = common * tmtyz;
        I_25xz = common * tmtxz;

        sxx[coord] += I_03xx*s_03 + I_13xx*s_13 + I_05xx*s_05 +
        I_15xx*s_15 + I_25xx*s_25;

        syy[coord] += I_03yy*s_03 + I_13yy*s_13 + I_05yy*s_05 +
        I_15yy*s_15 + I_25yy*s_25;

        szz[coord] += I_03zz*s_03 + I_13zz*s_13 + I_05zz*s_05 +
        I_15zz*s_15 + I_25zz*s_25;

        sxy[coord] += I_03xy*s_03 + I_13xy*s_13 + I_05xy*s_05 +
        I_15xy*s_15 + I_25xy*s_25;

        syz[coord] += I_03yz*s_03 + I_13yz*s_13 + I_05yz*s_05 +
        I_15yz*s_15 + I_25yz*s_25;

        sxz[coord] += I_03xz*s_03 + I_13xz*s_13 + I_05xz*s_05 +
        I_15xz*s_15 + I_25xz*s_25;
        
        }
    }
}