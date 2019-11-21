#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <mex.h>
#include "UtildaMex.h"

/*#ifdef _WIN32
    #define isnan(x) (_isnan(x))
    #define isfinite(x) (_finite(x))
    // this function does rounding as MS visual studio can't do it!
    int round( double r ) {
    return (r > 0.0) ? (r + 0.5) : (r - 0.5); 
    }
#endif
*/

/************************** MEX gateway function ***********************/

void mexFunction(int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
    double *x0,*y0,*z0;
    double *x1,*y1,*z1;
    double *x2,*y2,*z2;
    double *bx,*by,*bz;
    double *nx,*ny,*nz;
    double NU;
    int point_array_length;
    int segments_array_length;
    double *Ux,*Uy,*Uz;
    
    /* field points */
    x0 = (double *) mxGetPr(prhs[0]);
    y0 = (double *) mxGetPr(prhs[1]);
    z0 = (double *) mxGetPr(prhs[2]);
    
    /* segments burgers vector*/
    bx = (double *) mxGetPr(prhs[3]);
    by = (double *) mxGetPr(prhs[4]);
    bz = (double *) mxGetPr(prhs[5]);
    
    /* segments starting nodes */
    x1 = (double *) mxGetPr(prhs[6]);
    y1 = (double *) mxGetPr(prhs[7]);
    z1 = (double *) mxGetPr(prhs[8]);
    
    /* segments end nodes */
    x2 = (double *) mxGetPr(prhs[9]);
    y2 = (double *) mxGetPr(prhs[10]);
    z2 = (double *) mxGetPr(prhs[11]);
    
    /*slip planes*/
    nx = (double *) mxGetPr(prhs[12]);
    ny = (double *) mxGetPr(prhs[13]);
    nz = (double *) mxGetPr(prhs[14]);

    /* material constants */
    NU = mxGetScalar(prhs[15]);
    
    /* array sizes */
    point_array_length = mxGetScalar(prhs[16]);
    segments_array_length = mxGetScalar(prhs[17]);
    
    /* pre-allocate memory for results */
    plhs[0] = mxCreateDoubleMatrix(point_array_length,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(point_array_length,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(point_array_length,1,mxREAL);
    
    Ux = (double *) mxGetPr(plhs[0]);
    Uy = (double *) mxGetPr(plhs[1]);
    Uz = (double *) mxGetPr(plhs[2]);
    
    UtildaMex(x0,y0,z0,bx,by,bz,x1,y1,z1,x2,y2,z2,nx,ny,nz,
              NU,point_array_length,segments_array_length,Ux,Uy,Uz);
}

/**************************************************************************/

void UtildaMex(double *x0, double *y0, double *z0, double *bx, double *by, double *bz,
               double *x1, double *y1, double *z1, double *x2, double *y2, double *z2,
               double *nx, double *ny, double *nz,
               double NU, int point_array_length, int segments_array_length,
               double *Ux, double *Uy, double *Uz)
{
    int i,j,k;
    double coord[3];
    double lvec[3];
    double bvec[3];
    double nvec[3];
    double nnorm;
    double bdotbinv,TminusAdotb,TminusBdotb;
    double Apfactor,Bpfactor;
    double A[3], B[3], Ap[3], Bp[3], T[3];
    double utilda1[3], utilda2[3], utilda[3], utilda_tot[3];
    int singularity_check;
    const double P0[3] = {0,0,0};
    const double P1[3] = {sqrt(3),0.5,2+sqrt(3)}; /*random - more robust way??*/
    
    /*make sure U array is zeroed!*/
    for (i=0; i<point_array_length; i++)
    {
        Ux[i] = 0;
        Uy[i] = 0;
        Uz[i] = 0;
    }
    
    /*start nested segment-coordinate loops*/
    for (i=0; i<segments_array_length; i++)
    {
        /*assign b-vector*/
        bvec[0] = bx[i];
        bvec[1] = by[i];
        bvec[2] = bz[i];
        
        /*assign A,B*/
        A[0] = x1[i];
        A[1] = y1[i];
        A[2] = z1[i];
        
        B[0] = x2[i];
        B[1] = y2[i];
        B[2] = z2[i];
        
        /*calculate segment vector*/
        lvec[0] = B[0] - A[0];
        lvec[1] = B[1] - A[1];
        lvec[2] = B[2] - A[2];
        
        /*calculate slip plane normal, n=cross(l,b) */
        nvec[0] = lvec[1]*bvec[2] - lvec[2]*bvec[1];
        nvec[1] = -(lvec[0]*bvec[2] - lvec[2]*bvec[0]);
        nvec[2] = lvec[0]*bvec[1] - lvec[1]*bvec[0];
        
        nnorm = sqrt(nvec[0]*nvec[0] + nvec[1]*nvec[1] + nvec[2]*nvec[2]);
        nvec[0] = nvec[0]/nnorm;
        nvec[1] = nvec[1]/nnorm;
        nvec[2] = nvec[2]/nnorm;
        
        /*checks on slip plane to guarantee numerical stability*/
        SlipPlaneCheck(nvec,nnorm,nx[i],ny[i],nz[i]);
        
        /*printf("n[0] = %f, n[1] = %f, n[2] = %f \n",nvec[0],nvec[1],nvec[2]);*/
        
        /*find intersection point between slip plane and arbitrary vector*/
        /*the arbitrary vector is hard-coded! Has warnings if something is wrong though*/
        PlaneLineIntersection(nvec,A,P0,P1,T);
        
        /*printf("Cp[0] = %f, Cp[1] = %f, Cp[2] = %f \n",T[0],T[1],T[2]);*/
        
        /*find orthogonal projection of closure point Cp on the line collinear with
         b which passes through a node*/
        bdotbinv = 1.0/(bvec[0]*bvec[0] + bvec[1]*bvec[1] + bvec[2]*bvec[2]);
        TminusAdotb = (T[0]-A[0])*bvec[0] + (T[1]-A[1])*bvec[1] + (T[2]-A[2])*bvec[2];
        TminusBdotb = (T[0]-B[0])*bvec[0] + (T[1]-B[1])*bvec[1] + (T[2]-B[2])*bvec[2];
        Apfactor = TminusAdotb*bdotbinv;
        Bpfactor = TminusBdotb*bdotbinv;
        
        Ap[0] = A[0] + Apfactor*bvec[0];
        Ap[1] = A[1] + Apfactor*bvec[1];
        Ap[2] = A[2] + Apfactor*bvec[2];
        
        Bp[0] = B[0] + Bpfactor*bvec[0];
        Bp[1] = B[1] + Bpfactor*bvec[1];
        Bp[2] = B[2] + Bpfactor*bvec[2];
        
        /*printf("Ap[0] = %f, Ap[1] = %f, Ap[2] = %f \n",Ap[0],Ap[1],Ap[2]);
        printf("Bp[0] = %f, Bp[1] = %f, Bp[2] = %f \n",Bp[0],Bp[1],Bp[2]);*/
        
        for (j=0; j<point_array_length; j++)
        {
            coord[0] = x0[j];
            coord[1] = y0[j];
            coord[2] = z0[j];
            
            /* Calculate displacement due to the two subtriangles associated with each
            * segment
            *
            *   A --> B
            *   ^ \   |
            *   |  \  |
            *   |   \ |
            *   Ap<--Bp      . coord where displacement is needed
            *
            * In other words, of triangles Bp->A->B and Ap->A->Bp with point of interest*/
        
            BarnettTriangle(coord,Bp,A,B,bvec,NU,utilda1);
            BarnettTriangle(coord,Ap,A,Bp,bvec,NU,utilda2);
            utilda[0] = utilda1[0] + utilda2[0];
            utilda[1] = utilda1[1] + utilda2[1];
            utilda[2] = utilda1[2] + utilda2[2];
            
            if ((!isfinite(utilda[0])) || (!isfinite(utilda[1])) || (!isfinite(utilda[2])))
            {
                utilda[0] = 0;
                utilda[1] = 0;
                utilda[2] = 0;
            }
            
            /*printf("ux = %f, uy = %f, uz = %f \n",utilda[0],utilda[1],utilda[2]);*/
            
            Ux[j] += utilda[0];
            Uy[j] += utilda[1];
            Uz[j] += utilda[2];
        }
    }
}

/**************************************************************************/
                            
void SlipPlaneCheck(double *nvec, double nnorm, double nx, double ny, double nz)
{
    int check1,check2,check3;
    
    check1 = (isnan(nvec[0]) || isnan(nvec[1]) || isnan(nvec[2])); /*NaN*/
    check2 = (nvec[0]==0 && nvec[1]==0 && nvec[2]==0); /*Undefined*/
    check3 = (nnorm < 1E-10); /*Inf*/
    
    if (check1 || check2 || check3)
    {
        /*take slip plane given in input file*/
        nvec[0] = nx;
        nvec[1] = ny;
        nvec[2] = nz;
        nnorm = sqrt(nvec[0]*nvec[0] + nvec[1]*nvec[1] + nvec[2]*nvec[2]);
        nvec[0] = nvec[0]/nnorm;
        nvec[1] = nvec[1]/nnorm;
        nvec[2] = nvec[2]/nnorm;
        /*check input values as well*/
        check1 = (isnan(nvec[0]) || isnan(nvec[1]) || isnan(nvec[2])); /*NaN*/
        check2 = (nvec[0]==0 && nvec[1]==0 && nvec[2]==0); /*Undefined*/
        check3 = (nnorm < 1E-10); /*Inf*/
        if (check1 || check2 || check3)
        {
            /*Calculation cannot yield answer, and input is also incorrect. The world hates you*/
            printf("UtildaMex.c: Incorrect input slip plane vector in links array! Put dummy for now. Please fix! \n");
            nvec[0] = 1;
            nvec[1] = 1;
            nvec[2] = 1;
        }
    }
}
    
/**************************************************************************/

void PlaneLineIntersection(double *n, double *V0, double const *P0, double const *P1, double *I)
{
    /* Inputs:
     * n - normal vector of the plane
     * V0 - point that belongs to the plane
     * P0 - end point 1 of arbitrary segment P0P1
     * P1 - end point 2 of arbitrary segment P0P1
     * 
     * P0P1 is kept constant for all slip planes. Not sure if this is the best
     * way to do it, as it entails some checks, which may be inefficient.
     *
     * Outputs:
     * I - point of intersection
     * Check indicator of:
     *              0 => disjoint (no intersection) or lies in the plane
     *              1 => the plane intersects P0P1 in the unique point I
     */
    
    double u[3],w[3];
    double D,N,sI;
    
    /*preallocate I*/
    I[0]=0;
    I[1]=0;
    I[2]=0;
    
    u[0] = P1[0] - P0[0];
    u[1] = P1[1] - P0[1];
    u[2] = P1[2] - P0[2];
    D = n[0]*u[0] + n[1]*u[1] + n[2]*u[2]; /*dot product*/
    
    /* check if line lies in plane */
    if (fabs(D) < 10E-7)
    {
        printf("UtildaMex.c: Choose different arbitrary vector for displacement calculations \n");
        return;
    }
    
    w[0] = P0[0] - V0[0];
    w[1] = P0[1] - V0[1];
    w[2] = P0[2] - V0[2];
    N = -(n[0]*w[0] + n[1]*w[1] + n[2]*w[2]); /*dot product*/
    
    sI = N/D;
    I[0] = P0[0] + sI*u[0];
    I[1] = P0[1] + sI*u[1];
    I[2] = P0[2] + sI*u[2];
}


/**************************************************************************/

void BarnettTriangle(double *point, double *A, double *B, double *C, double *b, double NU, double *utilda)
{
    const double C1 = (1.0-2.0*NU)/(8.0*M_PI*(1.0-NU));
    const double C2 = 1.0/(8.0*M_PI*(1.0-NU));
    
    /* a -----> b
       ^      /
        \    /
         \  /
          c   */
    
    /* A is p1, B is p2, C is centroid */
    
    double RA[3], RB[3], RC[3];
    double vecAB[3], vecBC[3], vecCA[3];
    int i;
    double normAB, normBC, normCA;
    double modRA, modRB, modRC;
    double modRAinv, modRBinv, modRCinv;
    
    double tAB[3], tBC[3], tCA[3];
    double lamA[3], lamB[3], lamC[3];
    
    double fAB[3], fBC[3], fCA[3];
    double gAB[3], gBC[3], gCA[3];
    
    double omega;
    double eps;
    double n[3];
    
    for (i=0; i<3; i++){
        RA[i] = A[i] - point[i];
        RB[i] = B[i] - point[i];
        RC[i] = C[i] - point[i];
        
        vecAB[i] = B[i] - A[i];
        vecBC[i] = C[i] - B[i];
        vecCA[i] = A[i] - C[i];
    }
    
    safenorm(vecAB, tAB);
    safenorm(vecBC, tBC);
    safenorm(vecCA, tCA);
    
    modRA = sqrt(RA[0]*RA[0] + RA[1]*RA[1] + RA[2]*RA[2]);
    modRB = sqrt(RB[0]*RB[0] + RB[1]*RB[1] + RB[2]*RB[2]);
    modRC = sqrt(RC[0]*RC[0] + RC[1]*RC[1] + RC[2]*RC[2]);
    
    /*unit tangents along the directed segments AB, BC, and CA*/
    safenorm(RA, lamA);
    safenorm(RB, lamB);
    safenorm(RC, lamC);

    /*calculate fAB, fBC, and fCA*/
    n[0] = tAB[1]*tBC[2] - tAB[2]*tBC[1];
    n[1] = tAB[2]*tBC[0] - tAB[0]*tBC[2];
    n[2] = tAB[0]*tBC[1] - tAB[1]*tBC[0];
    
    fab(b,tAB,lamA,lamB,modRA,modRB,fAB);
    fab(b,tBC,lamB,lamC,modRB,modRC,fBC);
    fab(b,tCA,lamC,lamA,modRC,modRA,fCA);
    
    /*calculate gAB, gBC, and gCA*/
    gab(b,lamA,lamB,gAB);
    gab(b,lamB,lamC,gBC);
    gab(b,lamC,lamA,gCA);
    
    /*calculate solid angle*/
    omega = solang(lamA,lamB,lamC,point,n);
    
    /*printf("Omega = %f \n",omega);*/
    for (i=0; i<3; i++){
        utilda[i] = -b[i]*omega/(4.0*M_PI) - C1*(fAB[i] + fBC[i] + fCA[i]) + C2*(gAB[i] + gBC[i] + gCA[i]);
    }
    
}

/*****************************Auxiliary functions******************************/

/*f vector calculation*/
void fab(double *b, double *tAB, double *lamA, double *lamB, double RA, double RB, double *f_vec)
{
    double numerator ;
    double denominator;
    double logarithm;
    double eps;
    
    eps = 1E-10;
    numerator = RB*(1.0 + lamB[0]*tAB[0] + lamB[1]*tAB[1] + lamB[2]*tAB[2]);
    denominator = RA*(1.0 + lamA[0]*tAB[0] + lamA[1]*tAB[1] + lamA[2]*tAB[2]);
    
    /*this can happen if field point is on the dislocation segment*/
    if (fabs(denominator) < eps)
        logarithm =0;
    else if (fabs(numerator) < eps)
        logarithm =0;
    else
        logarithm = log(numerator/denominator);
    
    f_vec[0] = logarithm*(b[1]*tAB[2] - b[2]*tAB[1]);
    f_vec[1] = logarithm*(b[2]*tAB[0] - b[0]*tAB[2]);
    f_vec[2] = logarithm*(b[0]*tAB[1] - b[1]*tAB[0]);
    
}

/*g vector calculation*/
void gab(double *b, double *lamA, double *lamB, double *g_vec)
{
    /*const double invdenom = 1.0/(1.0 + (lamA[0]*lamB[0] + lamA[1]*lamB[1] + lamA[2]*lamB[2]));*/
    double denominator;
    
    double lamAxlamB[3];
    double bdlamAxlamB;
    double numerator[3];
    const double eps = 1E-12;
    
    lamAxlamB[0] = lamA[1]*lamB[2] - lamA[2]*lamB[1];
    lamAxlamB[1] = lamA[2]*lamB[0] - lamA[0]*lamB[2];
    lamAxlamB[2] = lamA[0]*lamB[1] - lamA[1]*lamB[0];
    
    bdlamAxlamB = b[0]*lamAxlamB[0] + b[1]*lamAxlamB[1] + b[2]*lamAxlamB[2];
    
    numerator[0] = bdlamAxlamB*(lamA[0] + lamB[0]);
    numerator[1] = bdlamAxlamB*(lamA[1] + lamB[1]);
    numerator[2] = bdlamAxlamB*(lamA[2] + lamB[2]);
    
    denominator = 1.0 + (lamA[0]*lamB[0] + lamA[1]*lamB[1] + lamA[2]*lamB[2]);
    
    if (fabs(denominator) < eps)
    {
        g_vec[0] = 0;
        g_vec[1] = 0;
        g_vec[2] = 0;
    }
    else
    {
        g_vec[0] = numerator[0]/denominator;
        g_vec[1] = numerator[1]/denominator;
        g_vec[2] = numerator[2]/denominator;
    }
    
}

/*solid angle calculation*/
double solang(double *lamA, double *lamB, double *lamC, double *p, double *n)
{
    double omega;
    const double a = acos( lamB[0]*lamC[0] + lamB[1]*lamC[1] + lamB[2]*lamC[2] ); /*a,b or c can be NaN when Ap ~ Bp*/
    const double b = acos( lamC[0]*lamA[0] + lamC[1]*lamA[1] + lamC[2]*lamA[2] );
    const double c = acos( lamA[0]*lamB[0] + lamA[1]*lamB[1] + lamA[2]*lamB[2] );
    
    const double s = (a+b+c)/2;
    double svec[4];
    int i;
    double temp, lamAdotn, eps;
    double sign;
    
    svec[0] = s*0.5;
    svec[1] = (s-a)*0.5;
    svec[2] = (s-b)*0.5;
    svec[3] = (s-c)*0.5;
    
    /*numerical accuracy may produce very small -ve numbers*/
    for (i=0; i<4; i++){
        if (svec[i] < 0)
            svec[i] = fabs(svec[i]);
    }
    temp = tan(svec[0])*tan(svec[1])*tan(svec[2])*tan(svec[3]);
    omega = 4*atan(sqrt(temp));
    
    lamAdotn = lamA[0]*n[0]+lamA[1]*n[1]+lamA[2]*n[2];
    
    eps = 1E-10;
    if (lamAdotn > eps)
    {
        sign =  1;
        omega = -sign*omega;
    }
    else if (lamAdotn < -eps)
    {
        sign = -1;
        omega = -sign*omega;
    }
    else
    {
        sign = 0;
        omega =0;
    }
    
    return omega;
}

/* This function exchange two rows of a matrix */
void safenorm( double R[3],double unitR[3])
{
    double normR;
    const double eps= 1E-10;
    
    normR = sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]);
    
    if (normR > eps)
    {     
        unitR[0] = R[0]/normR;	
        unitR[1] = R[1]/normR;	
        unitR[2] = R[2]/normR;
    }
    else
    {	
        unitR[0] = 0;	
        unitR[1] = 0;	
        unitR[2] = 0;
    }
}