#define _USE_MATH_DEFINES
#include <math.h>
#include <mex.h>
#include <matrix.h>

void MinDistCalc(double *x0, double *x1, double *y0, double *y1, double *vx0, double *vx1, double *vy0, double *vy1, double *dist2, double *ddist2dt, double *L1, double *L2);

/************************** MEX gateway function ***********************/

void mexFunction(int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
    double *seg_x0, *seg_y0, *seg_z0;
    double *seg_vx0, *seg_vy0, *seg_vz0;
    double *seg_x1, *seg_y1, *seg_z1;
    double *seg_vx1, *seg_vy1, *seg_vz1;
    
    double x0[3], x1[3], y0[3], y1[3];
    double vx0[3], vx1[3], vy0[3], vy1[3];
    
    double dist2[1]={0}, ddist2dt[1]={0}, L1[1]={0}, L2[1]={0}, logic[1]={0};
    double *dist2_array, *ddist2dt_array, *L1_array, *L2_array;
    int xlen,xlen2,p,j,k;
    
    /* node 1 */
    seg_x0 = (double *) mxGetPr(prhs[0]);
    seg_y0 = (double *) mxGetPr(prhs[1]);
    seg_z0 = (double *) mxGetPr(prhs[2]);
    seg_vx0 = (double *) mxGetPr(prhs[3]);
    seg_vy0 = (double *) mxGetPr(prhs[4]);
    seg_vz0 = (double *) mxGetPr(prhs[5]);
    
    /* node 2 */
    seg_x1 = (double *) mxGetPr(prhs[6]);
    seg_y1 = (double *) mxGetPr(prhs[7]);
    seg_z1 = (double *) mxGetPr(prhs[8]);
    seg_vx1 = (double *) mxGetPr(prhs[9]);
    seg_vy1 = (double *) mxGetPr(prhs[10]);
    seg_vz1 = (double *) mxGetPr(prhs[11]);
    
    xlen = mxGetNumberOfElements(prhs[0]);
    xlen2 = xlen*xlen;
    
    plhs[0] = mxCreateDoubleMatrix(xlen2,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(xlen2,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(xlen2,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(xlen2,1,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);
    
    dist2_array = (double *) mxGetPr(plhs[0]);
    ddist2dt_array = (double *) mxGetPr(plhs[1]);
    L1_array = (double *) mxGetPr(plhs[2]);
    L2_array = (double *) mxGetPr(plhs[3]);
    
    for (p=0; p<xlen; p++){
        
        x0[0]= seg_x0[p];
        x0[1]= seg_y0[p];
        x0[2]= seg_z0[p];
        
        x1[0]= seg_x1[p];
        x1[1]= seg_y1[p];
        x1[2]= seg_z1[p];
        
        vx0[0]= seg_vx0[p];
        vx0[1]= seg_vy0[p];
        vx0[2]= seg_vz0[p];
        
        vx1[0]= seg_vx1[p];
        vx1[1]= seg_vy1[p];
        vx1[2]= seg_vz1[p];
        
        for (j=0; j<xlen; j++){
            
            y0[0]= seg_x0[j];
            y0[1]= seg_y0[j];
            y0[2]= seg_z0[j];
            
            y1[0]= seg_x1[j];
            y1[1]= seg_y1[j];
            y1[2]= seg_z1[j];
            
            vy0[0]= seg_vx0[j];
            vy0[1]= seg_vy0[j];
            vy0[2]= seg_vz0[j];
            
            vy1[0]= seg_vx1[j];
            vy1[1]= seg_vy1[j];
            vy1[2]= seg_vz1[j];
            
            MinDistCalc(x0,x1,y0,y1,vx0,vx1,vy0,vy1,dist2,ddist2dt,L1,L2);
            
            k=p+j*xlen;
            dist2_array[k] = *dist2;
            ddist2dt_array[k] = *ddist2dt;
            L1_array[k] = *L1;
            L2_array[k] = *L2;

            }
        }
}

/**************************************************************************/

void MinDistCalc(double *x0, double *x1, double *y0, double *y1, double *vx0, double *vx1, double *vy0, double *vy1, double *dist2,double *ddist2dt, double *L1, double *L2)
{
    /* [dist2,ddist2dt,L1,L2]=mindistcalc(x0vx0,x1vx1,y0vy0,y1vy1)
     * this function finds the minimum distance bewtween two line segments
     * seg1=x0->x1 seg2=y0->y1
     * dist2 = square of the minimum distance between the two points
     * L1 = normalize position on seg1 that is closest to seg2
     * L2 = normalized position on seg2 that is closest to seg1
     * ddist2dt = time rate of change of the distance between L1 and L2*/
    
    double seg1[3], seg2[3];
    double vseg1[3], vseg2[3];
    double dist[4];
    double mindist2;
    double A, B, C, D, E, F, G, eps;
    int i,c,pos=0;
    
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
        
        mindist2 = dist[0];
        for ( c = 1 ; c < 4 ; c++ )
        {
            if ( dist[c] < mindist2 )
            {
                mindist2 = dist[c];
                pos = c+1;
            }
        }
        
        *L1=floor(pos*0.5);
        *L2= (pos-1)%2;
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