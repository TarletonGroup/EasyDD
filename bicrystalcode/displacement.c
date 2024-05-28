#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>

void triangleDisplacement(double *point, double *A, double *B, double *C, double *b, double *n, double NU, double *utilda);

/* Declare auxiliary functions */
void fab(double *b, double *t, double *lamA, double *lamB, double RA, double RB, double *f_vec);
void gab(double *b, double *lamA, double *lamB, double *g_vec);
void solang(double *lamA, double *lamB, double *lamC, double *p, double *plane_n, double *omega);

/************************** MEX gateway function ***********************/

int  main() {
    
    double point[] = {0,0,0};
    double A[] = {10,0,-100};
    double B[] = {10,0,100};
    double C[] = {10,10,0};
    double b[] = {0.5,1,1};
    double n[] = {-1,0,0};
    double utilda[3];
    double NU=0.305;
    
    triangleDisplacement(point,A,B,C,b,n,NU,utilda);
    
    //for (int i=0; i<3; i++){
    //    printf("%f \n",utilda[i]);
    //}

return 0;
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
    for (int i=0; i<3; i++){
        RA[i] = A[i] - point[i];
        RB[i] = B[i] - point[i];
        RC[i] = C[i] - point[i];
        
        vecAB[i] = B[i] - A[i];
        vecBC[i] = C[i] - B[i];
        vecCA[i] = A[i] - C[i];
    }
    
    const double modABinv = 1.0/sqrt(vecAB[0]*vecAB[0] + vecAB[1]*vecAB[1] + vecAB[2]*vecAB[2]);
    const double modBCinv = 1.0/sqrt(vecBC[0]*vecBC[0] + vecBC[1]*vecBC[1] + vecBC[2]*vecBC[2]);
    const double modCAinv = 1.0/sqrt(vecCA[0]*vecCA[0] + vecCA[1]*vecCA[1] + vecCA[2]*vecCA[2]);
    
    const double modRA = sqrt(RA[0]*RA[0] + RA[1]*RA[1] + RA[2]*RA[2]);
    const double modRB = sqrt(RB[0]*RB[0] + RB[1]*RB[1] + RB[2]*RB[2]);
    const double modRC = sqrt(RC[0]*RC[0] + RC[1]*RC[1] + RC[2]*RC[2]);
    const double modRAinv = 1.0/modRA;
    const double modRBinv = 1.0/modRB;
    const double modRCinv = 1.0/modRC;
    
    /*unit tangents along the directed segments AB, BC, and CA*/
    double tAB[3], tBC[3], tCA[3];
    double lamA[3], lamB[3], lamC[3];
    for (int i=0; i<3; i++){
        tAB[i] = vecAB[i]*modABinv;
        tBC[i] = vecBC[i]*modBCinv;
        tCA[i] = vecCA[i]*modCAinv;
        
        lamA[i] = RA[i]*modRAinv;
        lamB[i] = RB[i]*modRBinv;
        lamC[i] = RC[i]*modRCinv;
    }
    
//     printf("lamA \n");
//     for (int i=0; i<3; i++){
//         printf("%f \n",lamA[i]);
//     }
//     printf("lamB \n");
//     for (int i=0; i<3; i++){
//         printf("%f \n",lamB[i]);
//     }
//     printf("lamC \n");
//     for (int i=0; i<3; i++){
//         printf("%f \n",lamC[i]);
//     }
    
    /*calculate fAB, fBC, and fCA*/
    double fAB[3], fBC[3], fCA[3];
    fab(b,tAB,lamA,lamB,modRA,modRB,fAB);
    fab(b,tBC,lamB,lamC,modRB,modRC,fBC);
    fab(b,tCA,lamC,lamA,modRC,modRA,fCA);
    
//     printf("fAB \n");
//     for (int i=0; i<3; i++){
//         printf("%f \n",fAB[i]);
//     }
//     printf("fBC \n");
//     for (int i=0; i<3; i++){
//         printf("%f \n",fBC[i]);
//     }
//     printf("fCA \n");
//     for (int i=0; i<3; i++){
//         printf("%f \n",fCA[i]);
//     }
    
    /*calculate gAB, gBC, and gCA*/
    double gAB[3], gBC[3], gCA[3];
    gab(b,lamA,lamB,gAB);
    gab(b,lamB,lamC,gBC);
    gab(b,lamC,lamA,gCA);
    
//     printf("gAB \n");
//     for (int i=0; i<3; i++){
//         printf("%f \n",gAB[i]);
//     }
//     printf("gBC \n");
//     for (int i=0; i<3; i++){
//         printf("%f \n",gBC[i]);
//     }
//     printf("gCA \n");
//     for (int i=0; i<3; i++){
//         printf("%f \n",gCA[i]);
//     }
    
    /*calculate solid angle*/
    double omega[1];
    solang(lamA,lamB,lamC,point,plane_n,omega);
    
    //printf("Solid Angle \n");
    //printf("%f \n",omega[0]);
    
    for (int i=0; i<3; i++){
        utilda[i] = -b[i]*omega[0]/(4.0*M_PI) - C1*(fAB[i] + fBC[i] + fCA[i]) + C2*(gAB[i] + gBC[i] + gCA[i]);
    }
    
}

/*****************************Auxiliary functions******************************/

/*f vector calculation*/
void fab(double *b, double *t, double *lamA, double *lamB, double RA, double RB, double *f_vec)
{
    const double numerator = RB*(1.0 + lamB[0]*t[0] + lamB[1]*t[1] + lamB[2]*t[2]);
    const double denominator = RA*(1.0 + lamA[0]*t[0] + lamA[1]*t[1] + lamA[2]*t[2]);
    const double logarithm = log(numerator/denominator);
    
    f_vec[0] = logarithm*(b[1]*t[2] - b[2]*t[1]);
    f_vec[1] = logarithm*(b[2]*t[0] - b[0]*t[2]);
    f_vec[2] = logarithm*(b[0]*t[1] - b[1]*t[0]);

}

/*g vector calculation*/
void gab(double *b, double *lamA, double *lamB, double *g_vec)
{
    const double invdenom = 1.0/(1.0 + (lamA[0]*lamB[0] + lamA[1]*lamB[1] + lamA[2]*lamB[2]));
    
    double lamAxlamB[3];
    lamAxlamB[0] = lamA[1]*lamB[2] - lamA[2]*lamB[1];
    lamAxlamB[1] = lamA[2]*lamB[0] - lamA[0]*lamB[2];
    lamAxlamB[2] = lamA[0]*lamB[1] - lamA[1]*lamB[0];
    
    const double bdlamAxlamB = b[0]*lamAxlamB[0] + b[1]*lamAxlamB[1] + b[2]*lamAxlamB[2];
    
    double numerator[3];
    numerator[0] = bdlamAxlamB*(lamA[0] + lamB[0]);
    numerator[1] = bdlamAxlamB*(lamA[1] + lamB[1]);
    numerator[2] = bdlamAxlamB*(lamA[2] + lamB[2]);
    
    g_vec[0] = numerator[0]*invdenom;
    g_vec[1] = numerator[1]*invdenom;
    g_vec[2] = numerator[2]*invdenom;
    
}

/*solid angle calculation*/
void solang(double *lamA, double *lamB, double *lamC, double *p, double *plane_n, double *omega)
{
    /*determinant of [lamA;lamB;lamC]
     | lamA[1] lamA[2] lamA[3] |
     | lamB[1] lamB[2] lamB[3] |
     | lamC[1] lamC[2] lamC[3] | */
    const double cross1 = lamB[1]*lamC[2] - lamC[1]*lamB[2];
    const double cross2 = lamB[0]*lamC[2] - lamC[0]*lamB[2];
    const double cross3 = lamB[0]*lamC[1] - lamC[0]*lamB[1];
    const double detLam = lamA[0]*cross1 + lamA[1]*cross2 + lamA[2]*cross3;
    
    /*magnitude of lamA, lamB, lamC*/
    const double al = sqrt(lamA[0]*lamA[0] + lamA[1]*lamA[1] + lamA[2]*lamA[2]);
    const double bl = sqrt(lamB[0]*lamB[0] + lamB[1]*lamB[1] + lamB[2]*lamB[2]);
    const double cl = sqrt(lamC[0]*lamC[0] + lamC[1]*lamC[1] + lamC[2]*lamC[2]);
    
    /*dot(A,B) , dot(A,C) , dot(B,C)*/
    const double AdotB = lamA[0]*lamB[0] + lamA[1]*lamB[1] + lamA[2]*lamB[2];
    const double AdotC = lamA[0]*lamC[0] + lamA[1]*lamC[1] + lamA[2]*lamC[2];
    const double BdotC = lamB[0]*lamC[0] + lamB[1]*lamC[1] + lamB[2]*lamC[2];
    
    /*main formula*/
    const double div = al*bl*cl + AdotB*cl + AdotC*bl + BdotC*al;
    double at = atan(div/detLam);
    
    if (at<0.0)
        at = at + M_PI;
    
    omega[0] = 2.0*at;
    
    /*Test whether point is below or above plane!*/
    const double x = p[0]*plane_n[0] + p[1]*plane_n[1] + p[2]*plane_n[2] + plane_n[3];
    if (x<0.0)
        omega[0] = -omega[0];
}
