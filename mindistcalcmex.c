#define _USE_MATH_DEFINES
#include <math.h>
#include <mex.h>
#include <matrix.h>

void MinDistCalc(double *x0vx0, double *x1vx1, double *y0vy0, double *y1vy1,
    int x0vx0_length, double *dist2, double *ddist2dt, double *L1, double *L2);
void minimum(double *array, int size, double *minimum, int *location);

/************************** MEX gateway function ***********************/

void mexFunction(int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
    double *x0vx0, *x1vx1, *y0vy0, *y1vy1;
    double *dist2, *ddist2dt, *L1, *L2;
    int x0vx0_length;
    
    x0vx0 = (double *) mxGetPr(prhs[0]);
    x1vx1 = (double *) mxGetPr(prhs[1]);
    y0vy0 = (double *) mxGetPr(prhs[2]);
    y1vy1 = (double *) mxGetPr(prhs[3]);
    x0vx0_length = mxGetNumberOfElements(prhs[0]); 
    
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
    
    dist2= (double *) mxGetPr(plhs[0]);
    ddist2dt = (double *) mxGetPr(plhs[1]);
    L1 = (double *) mxGetPr(plhs[2]);
    L2 = (double *) mxGetPr(plhs[3]);
    
    MinDistCalc(x0vx0,x1vx1,y0vy0,y1vy1,x0vx0_length,dist2,ddist2dt,L1,L2);
}

/**************************************************************************/

void MinDistCalc(double *x0vx0, double *x1vx1, double *y0vy0, double *y1vy1,
                 int x0vx0_length, double *dist2, double *ddist2dt, double *L1, double *L2)
{
    
    /* [dist2,ddist2dt,L1,L2]=mindistcalc(x0vx0,x1vx1,y0vy0,y1vy1)
     * this function finds the minimum distance bewtween two line segments
     * seg1=x0->x1 seg2=y0->y1
     * dist2 = square of the minimum distance between the two points
     * L1 = normalize position on seg1 that is closest to seg2
     * L2 = normalized position on seg2 that is closest to seg1
     * ddist2dt = time rate of change of the distance between L1 and L2*/
    
    double x0[3], x1[3], y0[3], y1[3];
    double vx0[3], vx1[3], vy0[3], vy1[3];
    double seg1[3], seg2[3];
    double vseg1[3], vseg2[3];
    double dist[4];
	double mindist2[1]={0};
    double A, B, C, D, E, F, G, eps;
	int i,L1_check,L2_check,pos[1]={0};
    
    for (i=0;i<4;i++)
        dist[i]=0.0;
    
    A=0.0;
    B=0.0;
    C=0.0;
    D=0.0;
    E=0.0;
    F=0.0;
    G=0.0;
    
    for (i=0; i<3; i++)
    {
        x0[i]=x0vx0[i];
        x1[i]=x1vx1[i];
        y0[i]=y0vy0[i];
        y1[i]=y1vy1[i];
    }
    
    if (x0vx0_length==6){
        for (i=0; i<3; i++)
        {
            vx0[i]=x0vx0[i+3];
            vx1[i]=x1vx1[i+3];
            vy0[i]=y0vy0[i+3];
            vy1[i]=y1vy1[i+3];
        }
    }
    else {
        for (i=0; i<3; i++)
        {
            vx0[i]=0.0;
            vx1[i]=0.0;
            vy0[i]=0.0;
            vy1[i]=0.0;
        }
    }
    
    for (i=0; i<3; i++)
    {
        seg1[i]=x1[i]-x0[i];
        seg2[i]=y1[i]-y0[i];
        vseg1[i]=vx1[i]-vx0[i];
        vseg2[i]=vy1[i]-vy0[i];
        
        A += seg1[i]*seg1[i];
        B += 2.0*seg1[i]*(x0[i]-y0[i]);
        C += 2.0*seg1[i]*seg2[i];
        D += 2.0*seg2[i]*(y0[i]-x0[i]);
        E += seg2[i]*seg2[i];
        F += x0[i]*x0[i] + y0[i]*y0[i];
    }
    
    G=C*C-4.0*A*E;
    eps=1E-12;
    
    if (A<eps) {
        *L1=0.0;
        if (E<eps)
            *L2=0.0;
        else
            *L2=-0.5*D/E;
    }
    else if (E<eps) {
        *L2=0.0;
        if (A<eps)
            *L1=0.0;
        else
            *L1=-0.5*B/A;
    }
    else if (fabs(G)<eps) {
        for (i=0; i<3; i++){
            dist[0]+=(y0[i]-x0[i])*(y0[i]-x0[i]);
            dist[1]+=(y1[i]-x0[i])*(y1[i]-x0[i]);
            dist[2]+=(y0[i]-x1[i])*(y0[i]-x1[i]);
            dist[3]+=(y1[i]-x1[i])*(y1[i]-x1[i]);
        }
        
        minimum(dist,4,mindist2,pos);
        *L1=floor(*pos*0.5);
        *L2= (*pos-1)%2;
    }
    else {
        *L2=(2.0*A*D+B*C)/G;
        *L1=0.5*(C*L2[0]-B)/A;
    }
    
    /*now check to make sure that L2 and L1 are betwen 0 and 1*/
    if (*L1<0)
        *L1=0;
    else if (*L1>1)
        *L1=1;
    
    if (*L2<0)
        *L2=0;
    else if (*L2>1)
        *L2=1;
    
    *dist2=0;
    *ddist2dt=0;
    
    /* now calculate the distance^2 and the time rate of change of the distance between the points at L1 and L2*/
    for (i=0; i<3; i++){
        *dist2+=(x0[i]+seg1[i]*L1[0]-y0[i]-seg2[i]*L2[0])*(x0[i]+seg1[i]*L1[0]-y0[i]-seg2[i]*L2[0]);
        *ddist2dt+=2.0*((vx0[i]+vseg1[i]*L1[0]-vy0[i]-vseg2[i]*L2[0])*(x0[i]+seg1[i]*L1[0]-y0[i]-seg2[i]*L2[0]));
    }
}

void minimum(double *array, int size, double *min, int *location)
{
    int c;
    *min = array[0];
    for ( c = 1 ; c < size ; c++ )
    {
        if ( array[c] < *min )
        {
            *min = array[c];
            *location = c+1;
        }
    }
}