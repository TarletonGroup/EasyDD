#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <mex.h>

void triangleDisplacement(double *point, double *A, double *B, double *C, double *b, double *n, double NU, double *utilda);
/* windows VS compiler requires all variables be declared before use */
/* Declare auxiliary functions */
void fab(double *b, double *t, double *lamA, double *lamB, double RA, double RB, double *f_vec);
void gab(double *b, double *lamA, double *lamB, double *g_vec);
double solang(double *lamA, double *lamB, double *lamC, double *p, double *plane_n);

/************************** MEX gateway function ***********************/

void mexFunction(int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
    double *point;
    double *A, *B, *C;
    double *b, *n;
    double *utilda;
    double NU;
    
    point = (double *) mxGetPr(prhs[0]);
    A = (double *) mxGetPr(prhs[1]);
    B = (double *) mxGetPr(prhs[2]);
    C = (double *) mxGetPr(prhs[3]);
    b = (double *) mxGetPr(prhs[4]);
    n = (double *) mxGetPr(prhs[5]);
    NU = mxGetScalar(prhs[6]);
    
    plhs[0] = mxCreateDoubleMatrix(3,1,mxREAL);
    utilda = (double *) mxGetPr(plhs[0]);
    
    triangleDisplacement(point,A,B,C,b,n,NU,utilda);
}

/**************************************************************************/

void triangleDisplacement(double *point, double *A, double *B, double *C, double *b, double *plane_n, double NU, double *utilda)
{
    const double C1 = (1.0-2.0*NU)/(8.0*M_PI*(1.0-NU));
    const double C2 = 1.0/(8.0*M_PI*(1.0-NU));
    
    /*     a -----> b
            ^      /
             \    /
              \  /
               c          */
    
    /* A is p1, B is p2, C is centroid */
    
    double RA[3], RB[3], RC[3];
    double vecAB[3], vecBC[3], vecCA[3];
    int i;
    double modABinv, modBCinv, modCAinv;
    double modRA, modRB, modRC;
    double modRAinv, modRBinv, modRCinv; 
         
    double tAB[3], tBC[3], tCA[3];
    double lamA[3], lamB[3], lamC[3];
    
    double fAB[3], fBC[3], fCA[3];    
    double gAB[3], gBC[3], gCA[3];
    
    double omega;
     
    for (i=0; i<3; i++){
        RA[i] = A[i] - point[i];
        RB[i] = B[i] - point[i];
        RC[i] = C[i] - point[i];
        
        /*RA[i] = point[i] - A[i];
        RB[i] = point[i] - B[i];
        RC[i] = point[i] - C[i];*/
        
        vecAB[i] = B[i] - A[i];
        vecBC[i] = C[i] - B[i];
        vecCA[i] = A[i] - C[i];
    }
    
    modABinv = 1.0/sqrt(vecAB[0]*vecAB[0] + vecAB[1]*vecAB[1] + vecAB[2]*vecAB[2]);
    modBCinv = 1.0/sqrt(vecBC[0]*vecBC[0] + vecBC[1]*vecBC[1] + vecBC[2]*vecBC[2]);
    modCAinv = 1.0/sqrt(vecCA[0]*vecCA[0] + vecCA[1]*vecCA[1] + vecCA[2]*vecCA[2]);
    
    modRA = sqrt(RA[0]*RA[0] + RA[1]*RA[1] + RA[2]*RA[2]);
    modRB = sqrt(RB[0]*RB[0] + RB[1]*RB[1] + RB[2]*RB[2]);
    modRC = sqrt(RC[0]*RC[0] + RC[1]*RC[1] + RC[2]*RC[2]);
    modRAinv = 1.0/modRA;
    modRBinv = 1.0/modRB;
    modRCinv = 1.0/modRC;
    
    /*unit tangents along the directed segments AB, BC, and CA*/

    for (i=0; i<3; i++){
        tAB[i] = vecAB[i]*modABinv;
        tBC[i] = vecBC[i]*modBCinv;
        tCA[i] = vecCA[i]*modCAinv;
        
        lamA[i] = RA[i]*modRAinv;
        lamB[i] = RB[i]*modRBinv;
        lamC[i] = RC[i]*modRCinv;
    }
     /*
     printf("tAB \n");
     for (int i=0; i<3; i++){
         printf("%f \n",tAB[i]);
     } 
     
     printf("lamA \n");
     for (int i=0; i<3; i++){
         printf("%f \n",lamA[i]);
     }
     printf("lamB \n");
     for (int i=0; i<3; i++){
         printf("%f \n",lamB[i]);
     }
     printf("lamC \n");
     for (int i=0; i<3; i++){
         printf("%f \n",lamC[i]);
     }
     */
    
    /*calculate fAB, fBC, and fCA*/
   
    fab(b,tAB,lamA,lamB,modRA,modRB,fAB);
    fab(b,tBC,lamB,lamC,modRB,modRC,fBC);
    fab(b,tCA,lamC,lamA,modRC,modRA,fCA);
    
    /*
     printf("fAB \n");
     for (int i=0; i<3; i++){
         printf("%f \n",fAB[i]);
     }
     printf("fBC \n");
     for (int i=0; i<3; i++){
         printf("%f \n",fBC[i]);
     }
     printf("fCA \n");
     for (int i=0; i<3; i++){
         printf("%f \n",fCA[i]);
     }
    */
    
    /*calculate gAB, gBC, and gCA*/
    
    gab(b,lamA,lamB,gAB);
    gab(b,lamB,lamC,gBC);
    gab(b,lamC,lamA,gCA);
    
    /*
    printf("gAB \n");
     for (int i=0; i<3; i++){
        printf("%f \n",gAB[i]);
    }
    printf("gBC \n");
    for (int i=0; i<3; i++){
        printf("%f \n",gBC[i]);
    }
    printf("gCA \n");
    for (int i=0; i<3; i++){
        printf("%f \n",gCA[i]);
    }
     */
    
    /*calculate solid angle*/
    omega = solang(lamA,lamB,lamC,point,plane_n);

    /*
    printf("Omega = %f \n",omega);
    */
    
    for (i=0; i<3; i++){
        utilda[i] = -b[i]*omega/(4.0*M_PI) - C1*(fAB[i] + fBC[i] + fCA[i]) + C2*(gAB[i] + gBC[i] + gCA[i]);
        /*printf("f-sum = %f \n",fAB[i] + fBC[i] + fCA[i]);*/
        /*printf("g-sum = %f \n",gAB[i] + gBC[i] + gCA[i]);*/
    }
    
}

/*****************************Auxiliary functions******************************/

/*f vector calculation*/
void fab(double *b, double *tAB, double *lamA, double *lamB, double RA, double RB, double *f_vec)
{
    const double numerator = RB*(1.0 + lamB[0]*tAB[0] + lamB[1]*tAB[1] + lamB[2]*tAB[2]);
    const double denominator = RA*(1.0 + lamA[0]*tAB[0] + lamA[1]*tAB[1] + lamA[2]*tAB[2]);
    const double logarithm = log(numerator/denominator);
    
    f_vec[0] = logarithm*(b[1]*tAB[2] - b[2]*tAB[1]);
    f_vec[1] = logarithm*(b[2]*tAB[0] - b[0]*tAB[2]);
    f_vec[2] = logarithm*(b[0]*tAB[1] - b[1]*tAB[0]);

/* Fivel paper, might have typo?    
//     double bxt[3];
//     const double lamBdtAB = lamB[0]*tAB[0] + lamB[1]*tAB[1] + lamB[2]*tAB[2]; 
//     const double lamAdtAB = lamA[0]*tAB[0] + lamA[1]*tAB[1] + lamA[2]*tAB[2];
//     const double lnRARB = log(RB/RA);
//     
//     bxt[0] = b[1]*tAB[2] - b[2]*tAB[1];
//     bxt[1] = b[2]*tAB[0] - b[0]*tAB[2];
//     bxt[2] = b[0]*tAB[1] - b[1]*tAB[0];
//     
//     f_vec[0] = bxt[0] * lnRARB * (lamBdtAB / lamAdtAB);
//     f_vec[1] = bxt[1] * lnRARB * (lamBdtAB / lamAdtAB);
//     f_vec[2] = bxt[2] * lnRARB * (lamBdtAB / lamAdtAB);*/

}

/*g vector calculation*/
void gab(double *b, double *lamA, double *lamB, double *g_vec)
{
    const double invdenom = 1.0/(1.0 + (lamA[0]*lamB[0] + lamA[1]*lamB[1] + lamA[2]*lamB[2]));

    double lamAxlamB[3];
    double bdlamAxlamB;
    
    double numerator[3];
    
    lamAxlamB[0] = lamA[1]*lamB[2] - lamA[2]*lamB[1];
    lamAxlamB[1] = lamA[2]*lamB[0] - lamA[0]*lamB[2];
    lamAxlamB[2] = lamA[0]*lamB[1] - lamA[1]*lamB[0];

    bdlamAxlamB = b[0]*lamAxlamB[0] + b[1]*lamAxlamB[1] + b[2]*lamAxlamB[2];

    
    numerator[0] = bdlamAxlamB*(lamA[0] + lamB[0]);
    numerator[1] = bdlamAxlamB*(lamA[1] + lamB[1]);
    numerator[2] = bdlamAxlamB*(lamA[2] + lamB[2]);

    g_vec[0] = numerator[0]*invdenom;
    g_vec[1] = numerator[1]*invdenom;
    g_vec[2] = numerator[2]*invdenom;
    
}

/*solid angle calculation*/
double solang(double *lamA, double *lamB, double *lamC, double *p, double *plane_n)
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
  
    lamAdotn = lamA[0]*plane_n[0]+lamA[1]*plane_n[1]+lamA[2]*plane_n[2];
    eps = 1E-10;
    if (lamAdotn > 0)
       sign =  1;
    if (lamAdotn < 0)
       sign = -1;
    if (fabs(lamAdotn) < eps)
       sign = 0;
    omega = -sign*omega;  
    
     printf("sign = %f \n",sign);
    return omega;
}

/* void solang(double *lamA, double *lamB, double *lamC, double *p, double *plane_n, double *omega)
 {
     const double cross1 = lamB[1]*lamC[2] - lamC[1]*lamB[2];
     const double cross2 = lamB[0]*lamC[2] - lamC[0]*lamB[2];
     const double cross3 = lamB[0]*lamC[1] - lamC[0]*lamB[1];
     const double detLam = lamA[0]*cross1 + lamA[1]*cross2 + lamA[2]*cross3;
    
     const double al = sqrt(lamA[0]*lamA[0] + lamA[1]*lamA[1] + lamA[2]*lamA[2]);
     const double bl = sqrt(lamB[0]*lamB[0] + lamB[1]*lamB[1] + lamB[2]*lamB[2]);
     const double cl = sqrt(lamC[0]*lamC[0] + lamC[1]*lamC[1] + lamC[2]*lamC[2]);
     
     const double AdotB = lamA[0]*lamB[0] + lamA[1]*lamB[1] + lamA[2]*lamB[2];
     const double AdotC = lamA[0]*lamC[0] + lamA[1]*lamC[1] + lamA[2]*lamC[2];
     const double BdotC = lamB[0]*lamC[0] + lamB[1]*lamC[1] + lamB[2]*lamC[2];
     
     const double div = al*bl*cl + AdotB*cl + AdotC*bl + BdotC*al;
     //double at = atan(div/detLam);
     double at = atan2(detLam,div);
     
     if (at<0.0)
         at = at + M_PI;
     
     double sign;
     if (lamA[0]*plane_n[0]+lamA[1]*plane_n[1]+lamA[2]*plane_n[2] > 0)
        sign =  1;
     if (lamA[0]*plane_n[0]+lamA[1]*plane_n[1]+lamA[2]*plane_n[2] < 0)
        sign = -1;
     
      omega[0] = 2.0*at;
      omega[0] = -sign*omega[0];
             
      const double x = p[0]*plane_n[0] + p[1]*plane_n[1] + p[2]*plane_n[2] + plane_n[3];
      if (x<0.0)
          omega[0] = -omega[0];
     
}*/