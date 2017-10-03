#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <mex.h>

#ifdef _WIN32
/* this function calculates atanh since MS visual studio can't do it!*/
double atanh( double r ) {
    return 0.5*(log(1+r) - log(1-r));
}
#endif

void NodalSurfForceLinearRectangle(double *x1, double *x2, double *x3,
 	double *x4, double *x5, double *x6, double *b, double mu, double nu, double a,
    double *fx3, double *fx4, double *fx5, double *fx6, double *ftot);
 double dot8(double a[8],double b[8]);
 double dot3(double a[3],double b[3]);
 void cross3(double a[3], double b[3], double c[3]);

void mexFunction(int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
    double *x1, *x2, *x3, *x4, *x5, *x6;
    double *b;
    double mu,nu,a;
    double *fx3, *fx4, *fx5, *fx6, *ftot;

    x1 = (double *) mxGetPr(prhs[0]);
    x2 = (double *) mxGetPr(prhs[1]);
    x3 = (double *) mxGetPr(prhs[2]);
    x4 = (double *) mxGetPr(prhs[3]);
    x5 = (double *) mxGetPr(prhs[4]);
    x6 = (double *) mxGetPr(prhs[5]);
    b = (double *) mxGetPr(prhs[6]);
    mu = mxGetScalar(prhs[7]);
    nu = mxGetScalar(prhs[8]);
    a = mxGetScalar(prhs[9]);
    /*t = (double *) mxGetPr(prhs[10]);*/

    plhs[0] = mxCreateDoubleMatrix(3,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(3,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(3,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(3,1,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(3,1,mxREAL);

    fx3 = (double *) mxGetPr(plhs[0]);
    fx4 = (double *) mxGetPr(plhs[1]);
    fx5 = (double *) mxGetPr(plhs[2]);
    fx6 = (double *) mxGetPr(plhs[3]);
    ftot = (double *) mxGetPr(plhs[4]);
    
    NodalSurfForceLinearRectangle(x1,x2,x3,x4,x5,x6,b,mu,nu,a,fx3,fx4,fx5,fx6,ftot);
    
}

 void NodalSurfForceLinearRectangle(double *x1, double *x2, double *x3,
 	double *x4, double *x5, double *x6, double *b, double mu, double nu, double a,
 	double *fx3, double *fx4, double *fx5, double *fx6, double *ftot)
 {
    /* declaration of variables */
    int i,j;
    double p[3], q[3], t[3];
    double L, Lprime,tnorm;
    double n[3],pxt[3],qxt[3];
    double Rlim1[3],Rlim2[3];
    double y1,y2,r1,r2,s1,s2;
    double rv[8],sv[8],yv[8],signv[8];
    double pdott,qdott;
    double R[8][3];
    double Ra[8],Rdotp[8],Rdotq[8],Rdott[8];
    double A0m1[8],B0m1[8],C0m1[8];
    double A1m1[8],B1m1[8],C1m1[8];
    double A01[8],B01[8],C01[8];
    double A11[8],B11[8],C11[8];
    double udotu,vdotv,pdott2,qdott2;
    double D00m3[8],E00m3[8],F00m3[8];
    double temp1,cond1;
    double D01m3[8],D10m3[8];
    double E01m3[8],E10m3[8];
    double F01m3[8],F10m3[8];
    double D00m1[8],E00m1[8],F00m1[8];
    double D11m3[8],E11m3[8],F11m3[8];
    double D20m3[8],D02m3[8];
    double E20m3[8],F20m3[8];
    double D01m1[8],D10m1[8],E01m1[8],F01m1[8],E10m1[8],F10m1[8];
    double D001[8],E001[8],F001[8];
    double D11m1[8],E11m1[8],F11m1[8],D02m1[8],D20m1[8],E20m1[8],F20m1[8];
    double D12m3[8],D21m3[8],E21m3[8],F21m3[8],D22m3[8];
    double H001m3[8],H100m3[8],H010m3[8],H001m1[8],H100m1[8],H010m1[8];
    double H111m3[8],H210m3[8],H120m3[8],H021m3[8],H201m3[8];
    double H000m3[8],H000m1[8];
    double AA;
    double H101m3[8],H011m3[8],H110m3[8],H200m3[8],H020m3[8];
    double H001m5[8],H100m5[8],H010m5[8];
    double H101m5[8],H011m5[8],H110m5[8],H200m5[8],H020m5[8];
    double H111m5[8],H210m5[8],H120m5[8],H201m5[8],H021m5[8],H012m5[8],H102m5[8];
    double H112m5[8],H211m5[8],H121m5[8],H202m5[8],H022m5[8],H301m5[8],H031m5[8];
    double H122m5[8],H212m5[8],H221m5[8],H311m5[8],H131m5[8];
    double scH001m3,scH100m3,scH010m3;
    double scH111m3,scH210m3,scH120m3;
    double scH101m3,scH110m3,scH011m3,scH200m3,scH020m3;
    double scH001m5,scH100m5,scH010m5,scH101m5,scH011m5,scH110m5,scH200m5;
    double scH020m5,scH111m5,scH210m5,scH120m5,scH021m5,scH201m5,scH102m5;
    double scH012m5,scH112m5,scH211m5,scH121m5,scH202m5,scH022m5,scH301m5;
    double scH031m5,scH122m5,scH221m5,scH131m5,scH311m5,scH212m5;
    double tcrossb[3],pcrossb[3],qcrossb[3],bcrosst[3];
    double tdotn,factor,a2;
    double ttbn[3],tbtn[3],tpbn[3],tqbn[3];
    double ri,si;
    double I111m3[3],I120m3[3],I210m3[3],I111m5[3],I120m5[3],I210m5[3];
    double I131m5[3],I311m5[3],I122m5[3],I212m5[3],I221m5[3];
    double F111m3,F120m3,F210m3,F111m5,F120m5;
    double F210m5,F131m5,F311m5,F122m5,F212m5,F221m5;
    double fLLprime[3];
    
    /*calculate unit vectors p,q,t*/
    for (i=0;i<3;i++){
        t[i]=x2[i]-x1[i];
        p[i]=x4[i]-x3[i];
        q[i]=x5[i]-x3[i];
    }
    L=1.0/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    Lprime=1.0/sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
    tnorm=1.0/sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
    for (i=0;i<3;i++){
        p[i]=p[i]*L;
        q[i]=q[i]*Lprime;
        t[i]=t[i]*tnorm;
    }
    
    /*calculate n,pxt,qxt*/
    cross3(p,q,n);
    cross3(p,t,pxt);
    cross3(q,t,qxt);
    
    /*calculate y,r,s*/
    for (i=0;i<3;i++){
        Rlim1[i]=x3[i]-x1[i];
        Rlim2[i]=x6[i]-x2[i];
    }
    y1=dot3(Rlim1,n)/dot3(t,n);
    y2=dot3(Rlim2,n)/dot3(t,n);
    r1=dot3(Rlim1,qxt)/dot3(p,qxt);
    r2=dot3(Rlim2,qxt)/dot3(p,qxt);
    s1=dot3(Rlim1,pxt)/dot3(q,pxt);
    s2=dot3(Rlim2,pxt)/dot3(q,pxt);
    
    /*assign rv,sv,yv,signv*/
    for (i=0;i<4;i++)
        rv[i]=r2;
    for (i=4;i<8;i++)
        rv[i]=r1;
    
    sv[0]=s2;
    sv[1]=s2;
    sv[2]=s1;
    sv[3]=s1;
    sv[4]=s2;
    sv[5]=s2;
    sv[6]=s1;
    sv[7]=s1;
    
    yv[0]=y2;
    yv[1]=y1;
    yv[2]=y2;
    yv[3]=y1;
    yv[4]=y2;
    yv[5]=y1;
    yv[6]=y2;
    yv[7]=y1;
    
    signv[0]=1;
    signv[1]=-1;
    signv[2]=-1;
    signv[3]=1;
    signv[4]=-1;
    signv[5]=1;
    signv[6]=1;
    signv[7]=-1;
    
    /*calculate R*/
    pdott = dot3(p,t);
    qdott = dot3(q,t);
    
    for (i=0;i<8;i++)
    {
        for (j=0;j<3;j++)
        {
            R[i][j] = rv[i]*p[j]+sv[i]*q[j]+yv[i]*t[j];
        }
    }
    
    /*calculate Ra, Rdotp, Rdotq, Rdott*/
    
    for (i=0;i<8;i++)
    {
        Ra[i] = sqrt(rv[i]*rv[i] + sv[i]*sv[i] + yv[i]*yv[i] +
                     2.0*pdott*rv[i]*yv[i] + 2.0*qdott*sv[i]*yv[i]+a*a);
        Rdotp[i] = R[i][0]*p[0] + R[i][1]*p[1] + R[i][2]*p[2];
        Rdotq[i] = R[i][0]*q[0] + R[i][1]*q[1] + R[i][2]*q[2];
        Rdott[i] = R[i][0]*t[0] + R[i][1]*t[1] + R[i][2]*t[2];
    }
    
    /*calculate A0m1,B0m1,C0m1 linear seed integrals*/
    for (i=0;i<8;i++)
    {
        A0m1[i] = log(Ra[i]+Rdotp[i]);
        B0m1[i] = log(Ra[i]+Rdotq[i]);
        C0m1[i] = log(Ra[i]+Rdott[i]);
    }
    
    /*calculate linear integrals*/
    for (i=0;i<8;i++)
    {
        A1m1[i] = Ra[i] - pdott*yv[i]*A0m1[i];
        B1m1[i] = Ra[i] - qdott*yv[i]*B0m1[i];
        C1m1[i] = Ra[i] - (pdott*rv[i]+qdott*sv[i])*C0m1[i];
    }
    
    for (i=0;i<8;i++)
    {
        A01[i]= 0.5*(Rdotp[i]*Ra[i] + (Ra[i]*Ra[i]-(Rdotp[i]*Rdotp[i]))*A0m1[i]);
        B01[i]= 0.5*(Rdotq[i]*Ra[i] + (Ra[i]*Ra[i]-(Rdotq[i]*Rdotq[i]))*B0m1[i]);
        C01[i]= 0.5*(Rdott[i]*Ra[i] + (Ra[i]*Ra[i]-(Rdott[i]*Rdott[i]))*C0m1[i]);
    }
    
    for (i=0;i<8;i++)
    {
        A11[i]= 1.0/3.0*Ra[i]*Ra[i]*Ra[i] - pdott*yv[i]*A01[i];
        B11[i]= 1.0/3.0*Ra[i]*Ra[i]*Ra[i] - qdott*yv[i]*B01[i];
        C11[i]= 1.0/3.0*Ra[i]*Ra[i]*Ra[i] - (pdott*rv[i]+qdott*sv[i])*C01[i];
    }
    
    /*calculate double seed integrals */
    a2=a*a;
    pdott2=pdott*pdott;
    qdott2=qdott*qdott;
    udotu=(1.0-pdott2);
    vdotv=(1.0-qdott2);
    
    for (i=0;i<8;i++)
    {
        D00m3[i]=0;
        E00m3[i]=0;
        F00m3[i]=0;
    }
    
    for (i=0; i<8; i++)
    {
        cond1=(a2+(1-qdott2-pdott2)*yv[i]*yv[i]);
        if (cond1 > 0)
        {
            temp1=sqrt(cond1);
            D00m3[i] = (2.0/temp1 * atan((Ra[i] - Rdotp[i] + Rdotq[i])/temp1));
        }
        else if (cond1 < 0)
        {
            temp1=sqrt(fabs(cond1));
            D00m3[i]=(-2.0/temp1 * atanh((Ra[i] - Rdotp[i] + Rdotq[i])/temp1));
        }
    }
    
    for (i=0; i<8; i++)
    {
        cond1=((1-pdott2)*a2+(1.0-qdott2-pdott2)*sv[i]*sv[i]);
        if (cond1 > 0)
        {
            temp1=sqrt(cond1);
            E00m3[i] = (2.0/temp1 * atan(((1.0-pdott)*(Ra[i]-rv[i]+yv[i])+qdott*sv[i])/temp1));
        }
        else if (cond1 < 0)
        {
            temp1=sqrt(fabs(cond1));
            E00m3[i]=(-2.0/temp1 * atanh(((1.0-pdott)*(Ra[i]-rv[i]+yv[i])+qdott*sv[i])/temp1));
        }
    }
    
    for (i=0; i<8; i++)
    {
        cond1=((1-qdott2)*a2+(1-qdott2-pdott2)*rv[i]*rv[i]);
        if (cond1 > 0)
        {
            temp1=sqrt(cond1);
            F00m3[i] = (2.0/temp1 * atan(((1.0-qdott)*(Ra[i]-sv[i]+yv[i])+pdott*rv[i])/temp1));
        }
        else if (cond1 < 0)
        {
            temp1=sqrt(fabs(cond1));
            F00m3[i] = (-2.0/temp1 * atanh(((1.0-qdott)*(Ra[i]-sv[i]+yv[i])+pdott*rv[i])/temp1));
        }
    }
    
    /*calculate double integrals*/
    /* D_{ijl} = \iint r^i s^j Ra^l dr ds % m stands for '-'
     * E_{ikl} = \iint r^i y^k Ra^l dr dy
     * F_{jkl} = \iint s^j y^k Ra^l ds dy */
    for (i=0;i<8;i++)
    {
        D01m3[i]= -A0m1[i] - qdott*yv[i]*D00m3[i];
        D10m3[i]= -B0m1[i] - pdott*yv[i]*D00m3[i];
        E01m3[i]= -1.0/udotu * (A0m1[i] - pdott*C0m1[i] + qdott*sv[i]*E00m3[i]);
        F01m3[i]= -1.0/vdotv * (B0m1[i] - qdott*C0m1[i] + pdott*rv[i]*F00m3[i]);
        E10m3[i]= -C0m1[i]-pdott*E01m3[i];
        F10m3[i]= -C0m1[i]-qdott*F01m3[i];
        
        D00m1[i]= rv[i]*B0m1[i] + sv[i]*A0m1[i] - pdott*yv[i]*D10m3[i] - qdott*yv[i]*D01m3[i] - (a2+yv[i]*yv[i])*D00m3[i];
        E00m1[i]= rv[i]*C0m1[i] + yv[i]*A0m1[i] - qdott*sv[i]*E01m3[i] - (a2+sv[i]*sv[i])*E00m3[i] ;
        F00m1[i]= sv[i]*C0m1[i] + yv[i]*B0m1[i] - pdott*rv[i]*F01m3[i] - (a2+rv[i]*rv[i])*F00m3[i] ;
        D11m3[i]= -A1m1[i] - qdott*yv[i]*D10m3[i];
        E11m3[i]= 1.0/udotu *(-A1m1[i] + pdott*rv[i]*C0m1[i] - pdott*E00m1[i]-qdott*sv[i]*E10m3[i]);
        F11m3[i]= 1.0/vdotv *(-B1m1[i] + qdott*sv[i]*C0m1[i] - qdott*F00m1[i]-pdott*rv[i]*F10m3[i]);
        
        D20m3[i]= D00m1[i] - pdott*yv[i]*D10m3[i] - rv[i]*B0m1[i];
        D02m3[i]= D00m1[i] - qdott*yv[i]*D01m3[i] - sv[i]*A0m1[i];
        E20m3[i]= E00m1[i] - pdott*E11m3[i] - rv[i]*C0m1[i];
        F20m3[i]= F00m1[i] - qdott*F11m3[i] - sv[i]*C0m1[i];
        
        D01m1[i]= A01[i] - qdott*yv[i]*D00m1[i];
        D10m1[i]= B01[i] - pdott*yv[i]*D00m1[i];
        E01m1[i]= 1.0/udotu * (A01[i] - pdott*C01[i] - qdott*sv[i]*E00m1[i]) ;
        F01m1[i]= 1.0/vdotv * (B01[i] - qdott*C01[i] - pdott*rv[i]*F00m1[i]) ;
        E10m1[i]= C01[i] - pdott*E01m1[i];
        F10m1[i]= C01[i] - qdott*F01m1[i];
        
        D001[i]= 1.0/3.0*( (yv[i]*yv[i]+a2)*D00m1[i] + rv[i]*B01[i] + sv[i]*A01[i] + yv[i]*pdott*D10m1[i] +yv[i]*qdott*D01m1[i] ) ;
        E001[i]= 1.0/3.0*( (sv[i]*sv[i]+a2)*E00m1[i] + rv[i]*C01[i] + yv[i]*A01[i] + qdott*sv[i]*E01m1[i]);
        F001[i]= 1.0/3.0*( (rv[i]*rv[i]+a2)*F00m1[i] + sv[i]*C01[i] + yv[i]*B01[i] + pdott*rv[i]*F01m1[i]);

        D11m1[i]=A11[i]-qdott*yv[i]*D10m1[i];
        E11m1[i]=(1.0/udotu)*(A11[i]-pdott*rv[i]*C01[i] - qdott*sv[i]*E10m1[i]+pdott*E001[i]);
        F11m1[i]=(1.0/vdotv)*(B11[i]-qdott*sv[i]*C01[i] - pdott*rv[i]*F10m1[i]+qdott*F001[i]);
        D02m1[i]=sv[i]*A01[i]-D001[i]-qdott*yv[i]*D01m1[i];
        D20m1[i]=rv[i]*B01[i]-D001[i]-pdott*yv[i]*D10m1[i];
        E20m1[i]=rv[i]*C01[i]-E001[i]-pdott*E11m1[i];
        F20m1[i]=sv[i]*C01[i]-F001[i]-qdott*F11m1[i];
        
        D12m3[i]= D10m1[i] - qdott*yv[i]*D11m3[i] - sv[i]*A1m1[i];
        D21m3[i]= D01m1[i] - pdott*yv[i]*D11m3[i] - rv[i]*B1m1[i];
        E21m3[i]= 1.0/udotu * (E01m1[i] - pdott*E10m1[i] - rv[i]*C1m1[i] + pdott*yv[i]*A1m1[i] + qdott*pdott*sv[i]*E11m3[i]) ;
        F21m3[i]= 1.0/vdotv * (F01m1[i] - qdott*F10m1[i] - sv[i]*C1m1[i] + qdott*yv[i]*B1m1[i] + qdott*pdott*rv[i]*F11m3[i]) ;
        D22m3[i]=-D001[i]-qdott*yv[i]*D01m1[i]-pdott*yv[i]*D10m1[i]+pdott*qdott*yv[i]*yv[i]*D11m3[i]-rv[i]*sv[i]*Ra[i]+rv[i]*B01[i]+sv[i]*A01[i]+qdott*rv[i]*yv[i]*B1m1[i]+pdott*sv[i]*yv[i]*A1m1[i];
    }
    
    /*calculate triple integrals*/
    /* H_{ijkl} = \iiint r^i s^j y^k Ra^l drdsdy
     * m stands for '-' */
    AA = 1.0/(1-pdott2-qdott2);
    for (i=0; i<8; i++)
    {
        H001m3[i]= -1.0 * (D00m1[i] - qdott*E00m1[i] - pdott*F00m1[i]) /(1.0-pdott2-qdott2);
        H100m3[i]= -pdott*H001m3[i] - F00m1[i];
        H010m3[i]= -qdott*H001m3[i] - E00m1[i];
        
        H001m1[i]= 1.0/(1.0-pdott2-qdott2) * (D001[i] - pdott*F001[i] - qdott*E001[i]);
        H100m1[i]= F001[i] - pdott*H001m1[i];
        H010m1[i]= E001[i] - qdott*H001m1[i];
        
        H111m3[i]=(1.0/(1.0-pdott2-qdott2))*(pdott*rv[i]*F10m1[i]+qdott*sv[i]*E10m1[i]-pdott*H010m1[i]-qdott*H100m1[i]-D11m1[i]);
        H210m3[i]= H010m1[i] - rv[i]*F10m1[i] - pdott*H111m3[i];
        H120m3[i]= H100m1[i] - sv[i]*E10m1[i] - qdott*H111m3[i];
        H021m3[i]=1.0/(1.0-pdott2-qdott2)*(-2.0*qdott*H010m1[i]-D02m1[i]+qdott*sv[i]*sv[i]*E00m1[i]+pdott*F20m1[i]);
        H201m3[i]=1.0/(1.0-pdott2-qdott2)*(-2.0*pdott*H100m1[i]-D20m1[i]+pdott*rv[i]*rv[i]*F00m1[i]+qdott*E20m1[i]);

        /*H000m3 is the only integral for which there is no analytical
        * antiderivative. This integral could be very easily estimated by a Gauss
        * quadrature. This is not actually required at all as all the H000m3 terms cancel
        * out exactly in the final nodal force calculation. H000m3 is set to 0 for
        * convenience.*/
        H000m3[i]= 0;
        /*This means that all antiderivatives in relation with H000m3 are
        * modified and their value would disagree with a Gauss quadrature if H000m3
        * is not reset to its regular value. However the final sum and nodal forces
        * remain unchanged.*/
        
        H000m1[i]= 1.0/2.0*(rv[i]*F00m1[i] + sv[i]*E00m1[i] + yv[i]*D00m1[i] - a2*H000m3[i]);

        H101m3[i]= (-D10m1[i]+pdott*rv[i]*F00m1[i]+qdott*E10m1[i]-pdott*H000m1[i])*AA ;
        H011m3[i]= (-D01m1[i]+qdott*sv[i]*E00m1[i]+pdott*F10m1[i]-qdott*H000m1[i])*AA ;
        H110m3[i]= -E10m1[i] - qdott*H101m3[i];
        H200m3[i]= H000m1[i] - rv[i]*F00m1[i] - pdott*H101m3[i];
        H020m3[i]= H000m1[i] - sv[i]*E00m1[i] - qdott*H011m3[i];
        
        H001m5[i]=  (pdott*F00m3[i] + qdott*E00m3[i] - D00m3[i])/3.0 * AA ;
        H100m5[i]= -1.0/3.0*F00m3[i] - pdott*H001m5[i];
        H010m5[i]= -1.0/3.0*E00m3[i] - qdott*H001m5[i];
        
        H101m5[i]= (-pdott*H000m3[i]-D10m3[i]+pdott*rv[i]*F00m3[i]+qdott*E10m3[i])/3.0 * AA;
        H011m5[i]= (-qdott*H000m3[i]-D01m3[i]+qdott*sv[i]*E00m3[i]+pdott*F10m3[i])/3.0 * AA;
        H110m5[i]= -1.0/3.0*F10m3[i] - pdott*H011m5[i];
        H200m5[i]= 1.0/3.0*H000m3[i] - 1.0/3.0*rv[i]*F00m3[i] - pdott*H101m5[i];
        H020m5[i]= 1.0/3.0*H000m3[i] - 1.0/3.0*sv[i]*E00m3[i] - qdott*H011m5[i];
                          
        H111m5[i]=  (-D11m3[i] + pdott*rv[i]*F10m3[i] + qdott*sv[i]*E10m3[i] - pdott*H010m3[i] - qdott*H100m3[i])/3.0 * AA ;
        H210m5[i]= -1.0/3.0*rv[i]*F10m3[i] + 1.0/3.0*H010m3[i] - pdott*H111m5[i];
        H120m5[i]= -1.0/3.0*sv[i]*E10m3[i] + 1.0/3.0*H100m3[i] - qdott*H111m5[i];
        H201m5[i]= (-2.0*pdott*H100m3[i]-D20m3[i]+pdott*rv[i]*rv[i]*F00m3[i]+qdott*E20m3[i])/3.0 * AA;
        H021m5[i]= (-2.0*qdott*H010m3[i]-D02m3[i]+qdott*sv[i]*sv[i]*E00m3[i]+pdott*F20m3[i])/3.0 * AA;
        H012m5[i]= (H010m3[i] -qdott*H001m3[i]-yv[i]*D01m3[i]+qdott*sv[i]*E01m3[i]+pdott*F11m3[i])/3.0 * AA;
        H102m5[i]= (H100m3[i] -pdott*H001m3[i]-yv[i]*D10m3[i]+pdott*rv[i]*F01m3[i]+qdott*E11m3[i])/3.0 * AA;
                                            
        H112m5[i]= (-yv[i]*D11m3[i] + pdott*rv[i]*F11m3[i] + qdott*sv[i]*E11m3[i] + H110m3[i] - pdott*H011m3[i] - qdott*H101m3[i])/3 * AA;
        H211m5[i]= -1.0/3.0*rv[i]*F11m3[i] + 1.0/3.0*H011m3[i] - pdott*H112m5[i];
        H121m5[i]= -1.0/3.0*sv[i]*E11m3[i] + 1.0/3.0*H101m3[i] - qdott*H112m5[i];
        H202m5[i]= (H200m3[i]-2.0*pdott*H101m3[i]-yv[i]*D20m3[i]+pdott*rv[i]*rv[i]*F01m3[i]+qdott*E21m3[i])/3.0 * AA;
        H022m5[i]= (H020m3[i]-2.0*qdott*H011m3[i]-yv[i]*D02m3[i]+qdott*sv[i]*sv[i]*E01m3[i]+pdott*F21m3[i])/3.0 * AA;
        H301m5[i]= -1.0/3.0*rv[i]*rv[i]*F01m3[i] + 2.0/3.0*H101m3[i] - pdott*H202m5[i];
        H031m5[i]= -1.0/3.0*sv[i]*sv[i]*E01m3[i] + 2.0/3.0*H011m3[i] - qdott*H022m5[i];
                 
        H122m5[i]= (H120m3[i]-pdott*H021m3[i]-2.0*qdott*H111m3[i]-yv[i]*D12m3[i]+pdott*rv[i]*F21m3[i]+qdott*sv[i]*sv[i]*E11m3[i])/3.0 *AA;
        H212m5[i]= (H210m3[i]-qdott*H201m3[i]-2.0*pdott*H111m3[i]-yv[i]*D21m3[i]+qdott*sv[i]*E21m3[i]+pdott*rv[i]*rv[i]*F11m3[i])/3.0 * AA;
        H221m5[i]= (pdott*rv[i]*rv[i]*F20m3[i]+qdott*sv[i]*sv[i]*E20m3[i]-2.0*pdott*H120m3[i]-2.0*qdott*H210m3[i]-D22m3[i])/3.0 * AA;
        H311m5[i]= (2.0*(1.0-qdott*qdott)*H111m3[i]+qdott*pdott*H201m3[i]-pdott*H210m3[i]-(1.0-qdott*qdott)*rv[i]*rv[i]*F11m3[i]-pdott*qdott*sv[i]*E21m3[i]+pdott*yv[i]*D21m3[i])/3.0 * AA;
        H131m5[i]= (2.0*(1.0-pdott*pdott)*H111m3[i]+qdott*pdott*H021m3[i]-qdott*H120m3[i]-(1.0-pdott*pdott)*sv[i]*sv[i]*E11m3[i]-pdott*qdott*rv[i]*F21m3[i]+qdott*yv[i]*D12m3[i])/3.0 * AA;

    }
    
    /* scalarization of integrals */
    /*Every triple integrals are now evaluated.
    * if terms are not functions of r,s and y they disappear in the summation
    * only functions of the three variable matter */
    scH001m3= dot8(H001m3,signv);
    scH100m3= dot8(H100m3,signv);
    scH010m3= dot8(H010m3,signv);
    scH111m3= dot8(H111m3,signv);
    scH210m3= dot8(H210m3,signv);
    scH120m3= dot8(H120m3,signv);
    scH101m3= dot8(H101m3,signv);
    scH110m3= dot8(H110m3,signv);
    scH011m3= dot8(H011m3,signv);
    scH200m3= dot8(H200m3,signv);
    scH020m3= dot8(H020m3,signv);
    scH001m5= dot8(H001m5,signv);
    scH100m5= dot8(H100m5,signv);
    scH010m5= dot8(H010m5,signv);
    scH101m5= dot8(H101m5,signv);
    scH011m5= dot8(H011m5,signv);
    scH110m5= dot8(H110m5,signv);
    scH200m5= dot8(H200m5,signv);
    scH020m5= dot8(H020m5,signv);
    scH111m5= dot8(H111m5,signv);
    scH210m5= dot8(H210m5,signv);
    scH120m5= dot8(H120m5,signv);
    scH021m5= dot8(H021m5,signv);
    scH201m5= dot8(H201m5,signv);
    scH102m5= dot8(H102m5,signv);
    scH012m5= dot8(H012m5,signv);
    scH112m5= dot8(H112m5,signv);
    scH211m5= dot8(H211m5,signv);
    scH121m5= dot8(H121m5,signv);
    scH202m5= dot8(H202m5,signv);
    scH022m5= dot8(H022m5,signv);
    scH301m5= dot8(H301m5,signv);
    scH031m5= dot8(H031m5,signv);
    scH122m5= dot8(H122m5,signv);
    scH221m5= dot8(H221m5,signv); 
    scH131m5= dot8(H131m5,signv); 
    scH311m5= dot8(H311m5,signv); 
    scH212m5= dot8(H212m5,signv); 

    /*final sum*/
    /* the force is now calculated for each of the 4 patch points */
    cross3(t,b,tcrossb);
    cross3(p,b,pcrossb);
    cross3(q,b,qcrossb);
    cross3(b,t,bcrosst);
    tdotn = dot3(t,n);
    for (i=0;i<3;i++)
    {
        ttbn[i] = t[i]*dot3(tcrossb,n);
        tbtn[i] = t[i]*dot3(bcrosst,n);
        tpbn[i] = t[i]*dot3(pcrossb,n);
        tqbn[i] = t[i]*dot3(qcrossb,n);
    }
    factor=mu/4.0/3.141592653589793238462643383279502884197169399375105820974944592307816/(1.0-nu)*L*Lprime;

    /* Vectors (identical for all xi nodes of the surface element) */
    for (i=0;i<3;i++)
    {
        I111m3[i]=(tcrossb[i]*tdotn+t[i]*(dot3(tcrossb,n)))*(1.0-nu) + (bcrosst[i]*tdotn+tbtn[i]);
        I120m3[i]=(qcrossb[i]*tdotn+tqbn[i])*(1.0-nu)-(dot3(qcrossb,t)*n[i])+(q[i]*(dot3(bcrosst,n)));
        I210m3[i]=(pcrossb[i]*tdotn+tpbn[i])*(1.0-nu)-(dot3(pcrossb,t)*n[i])+(p[i]*(dot3(bcrosst,n)));
        I111m5[i]=(tcrossb[i]*tdotn+ttbn[i])*1.5*(1.0-nu)*a2;
        I120m5[i]=(qcrossb[i]*tdotn+tqbn[i])*1.5*(1.0-nu)*a2-(dot3(qcrossb,t)*n[i]) * 3.0*a2;
        I210m5[i]=(pcrossb[i]*tdotn+tpbn[i])*1.5*(1.0-nu)*a2-(dot3(pcrossb,t)*n[i]) * 3.0*a2;
        I131m5[i]=-(dot3(qcrossb,t)*tdotn)*q[i] * 3.0;
        I311m5[i]=-(dot3(pcrossb,t)*tdotn)*p[i] * 3.0;
        I122m5[i]=-(dot3(qcrossb,t)*tdotn)*t[i] * 3.0;
        I212m5[i]=-(dot3(pcrossb,t)*tdotn)*t[i] * 3.0;
        I221m5[i]=-((dot3(pcrossb,t)*tdotn)*q[i]+(dot3(qcrossb,t)*tdotn)*p[i])*3.0;
    }
        
    /* Force in x3 */
    ri= r2;
    si= s2;
    F111m3=scH111m3-si*scH101m3-ri*scH011m3+ri*si*scH001m3;
    F120m3=scH120m3-si*scH110m3-ri*scH020m3+ri*si*scH010m3;
    F210m3=scH210m3-si*scH200m3-ri*scH110m3+ri*si*scH100m3;
    F111m5=scH111m5-si*scH101m5-ri*scH011m5+ri*si*scH001m5;
    F120m5=scH120m5-si*scH110m5-ri*scH020m5+ri*si*scH010m5;
    F210m5=scH210m5-si*scH200m5-ri*scH110m5+ri*si*scH100m5;
    F131m5=scH131m5-si*scH121m5-ri*scH031m5+ri*si*scH021m5;
    F311m5=scH311m5-si*scH301m5-ri*scH211m5+ri*si*scH201m5;
    F122m5=scH122m5-si*scH112m5-ri*scH022m5+ri*si*scH012m5;
    F212m5=scH212m5-si*scH202m5-ri*scH112m5+ri*si*scH102m5;
    F221m5=scH221m5-si*scH211m5-ri*scH121m5+ri*si*scH111m5;
    for (i=0;i<3;i++)
    {
        fLLprime[i]= I111m3[i] * F111m3 + I210m3[i] * F210m3 + I120m3[i] * F120m3
        + I111m5[i] * F111m5 + I210m5[i] * F210m5 + I120m5[i] * F120m5
        + I212m5[i] * F212m5 + I122m5[i] * F122m5 + I221m5[i] * F221m5
        + I311m5[i] * F311m5 + I131m5[i] * F131m5;
        
        fx3[i] = fLLprime[i]*factor;
    }
    
    /* Force in x4 */
    ri= r1;
    si= s2;
    F111m3=-(scH111m3-si*scH101m3-ri*scH011m3+ri*si*scH001m3);
    F120m3=-(scH120m3-si*scH110m3-ri*scH020m3+ri*si*scH010m3);
    F210m3=-(scH210m3-si*scH200m3-ri*scH110m3+ri*si*scH100m3);
    F111m5=-(scH111m5-si*scH101m5-ri*scH011m5+ri*si*scH001m5);
    F120m5=-(scH120m5-si*scH110m5-ri*scH020m5+ri*si*scH010m5);
    F210m5=-(scH210m5-si*scH200m5-ri*scH110m5+ri*si*scH100m5);
    F131m5=-(scH131m5-si*scH121m5-ri*scH031m5+ri*si*scH021m5);
    F311m5=-(scH311m5-si*scH301m5-ri*scH211m5+ri*si*scH201m5);
    F122m5=-(scH122m5-si*scH112m5-ri*scH022m5+ri*si*scH012m5);
    F212m5=-(scH212m5-si*scH202m5-ri*scH112m5+ri*si*scH102m5);
    F221m5=-(scH221m5-si*scH211m5-ri*scH121m5+ri*si*scH111m5);
    /* Vectors (identical for all xi nodes of the surface element) */
    for (i=0;i<3;i++)
    {
        fLLprime[i]= I111m3[i] * F111m3 + I210m3[i] * F210m3 + I120m3[i] * F120m3
        + I111m5[i] * F111m5 + I210m5[i] * F210m5 + I120m5[i] * F120m5
        + I212m5[i] * F212m5 + I122m5[i] * F122m5 + I221m5[i] * F221m5
        + I311m5[i] * F311m5 + I131m5[i] * F131m5;
        
        fx4[i] = fLLprime[i]*factor;
    }
    
    /* Force in x5 */
    ri= r2;
    si= s1;
    F111m3=-(scH111m3-si*scH101m3-ri*scH011m3+ri*si*scH001m3);
    F120m3=-(scH120m3-si*scH110m3-ri*scH020m3+ri*si*scH010m3);
    F210m3=-(scH210m3-si*scH200m3-ri*scH110m3+ri*si*scH100m3);
    F111m5=-(scH111m5-si*scH101m5-ri*scH011m5+ri*si*scH001m5);
    F120m5=-(scH120m5-si*scH110m5-ri*scH020m5+ri*si*scH010m5);
    F210m5=-(scH210m5-si*scH200m5-ri*scH110m5+ri*si*scH100m5);
    F131m5=-(scH131m5-si*scH121m5-ri*scH031m5+ri*si*scH021m5);
    F311m5=-(scH311m5-si*scH301m5-ri*scH211m5+ri*si*scH201m5);
    F122m5=-(scH122m5-si*scH112m5-ri*scH022m5+ri*si*scH012m5);
    F212m5=-(scH212m5-si*scH202m5-ri*scH112m5+ri*si*scH102m5);
    F221m5=-(scH221m5-si*scH211m5-ri*scH121m5+ri*si*scH111m5);
    /* Vectors (identical for all xi nodes of the surface element) */
    for (i=0;i<3;i++)
    {
        fLLprime[i]= I111m3[i] * F111m3 + I210m3[i] * F210m3 + I120m3[i] * F120m3
        + I111m5[i] * F111m5 + I210m5[i] * F210m5 + I120m5[i] * F120m5
        + I212m5[i] * F212m5 + I122m5[i] * F122m5 + I221m5[i] * F221m5
        + I311m5[i] * F311m5 + I131m5[i] * F131m5;
        
        fx5[i] = fLLprime[i]*factor;
    }
    
    /* Force in x6 */
    ri= r1;
    si= s1;
    F111m3=scH111m3-si*scH101m3-ri*scH011m3+ri*si*scH001m3;
    F120m3=scH120m3-si*scH110m3-ri*scH020m3+ri*si*scH010m3;
    F210m3=scH210m3-si*scH200m3-ri*scH110m3+ri*si*scH100m3;
    F111m5=scH111m5-si*scH101m5-ri*scH011m5+ri*si*scH001m5;
    F120m5=scH120m5-si*scH110m5-ri*scH020m5+ri*si*scH010m5;
    F210m5=scH210m5-si*scH200m5-ri*scH110m5+ri*si*scH100m5;
    F131m5=scH131m5-si*scH121m5-ri*scH031m5+ri*si*scH021m5;
    F311m5=scH311m5-si*scH301m5-ri*scH211m5+ri*si*scH201m5;
    F122m5=scH122m5-si*scH112m5-ri*scH022m5+ri*si*scH012m5;
    F212m5=scH212m5-si*scH202m5-ri*scH112m5+ri*si*scH102m5;
    F221m5=scH221m5-si*scH211m5-ri*scH121m5+ri*si*scH111m5;
    /* Vectors (identical for all xi nodes of the surface element) */
    for (i=0;i<3;i++)
    {
        fLLprime[i]= I111m3[i] * F111m3 + I210m3[i] * F210m3 + I120m3[i] * F120m3
        + I111m5[i] * F111m5 + I210m5[i] * F210m5 + I120m5[i] * F120m5
        + I212m5[i] * F212m5 + I122m5[i] * F122m5 + I221m5[i] * F221m5
        + I311m5[i] * F311m5 + I131m5[i] * F131m5;
        
        fx6[i] = fLLprime[i]*factor;
    }

    /*total force*/
    for (i=0;i<3;i++)
    {
    ftot[i]=fx3[i]+fx4[i]+fx5[i]+fx6[i];
    }
 
}

/* Auxiliary functions */
double dot3(double a[3],double b[3])
{
    double c;
    c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    return c;
}

double dot8(double a[8],double b[8])
{
    int i;
    double c;
    c=0;
    for (i=0;i<8;i++)
    {
        c += a[i]*b[i];
    }
    return c;
}

void cross3(double a[3], double b[3], double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}