#define _USE_MATH_DEFINES
#include <math.h>
#include <mex.h>
#include <matrix.h>

/*#ifdef _WIN32
    // this function does rounding as MS visual studio can't do it!
    int round( double r ) {
    return (r > 0.0) ? (r + 0.5) : (r - 0.5); 
    }
#endif
*/
void CreateInputMex (double *rn_x,double *rn_y,double *rn_z, double *v_x,double *v_y,double *v_z, double *flag,
       double *links_c1, double *links_c2,int links_length,double *count,
       double *r11x,double *r11y,double *r11z, double *r21x,double *r21y,double *r21z, double *r12x,double *r12y,double *r12z, double *r22x,double *r22y,double *r22z,
       double *v11x,double *v11y,double *v11z, double *v21x,double *v21y,double *v21z, double *v12x,double *v12y,double *v12z, double *v22x,double *v22y,double *v22z,
       double *in1s1, double *in2s1, double *in1s2, double *n2s2,double *ii, double *ij);

/************************** MEX gateway function ***********************/

void mexFunction(int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
    
    double *rn_x, *rn_y, *rn_z;
    double *v_x, *v_y, *v_z;
    double *flag;
    double *links_c1, *links_c2;
    int links_length;
    double *r11x, *r11y, *r11z, *r21x, *r21y, *r21z, *r12x, *r12y, *r12z, *r22x, *r22y, *r22z;
    double *v11x, *v11y, *v11z, *v21x, *v21y, *v21z, *v12x, *v12y, *v12z, *v22x, *v22y, *v22z;
    double *in1s1, *in2s1, *in1s2, *in2s2, *ii, *ij; 
    int g;
    double *count=0;
    /*printf("Number of possibilities=%i \n",g);*/
    
    rn_x = (double *) mxGetPr(prhs[0]);
    rn_y = (double *) mxGetPr(prhs[1]);
    rn_z = (double *) mxGetPr(prhs[2]);
    v_x = (double *) mxGetPr(prhs[3]);
    v_y = (double *) mxGetPr(prhs[4]);
    v_z = (double *) mxGetPr(prhs[5]);
    flag = (double *) mxGetPr(prhs[6]);
    links_c1 = (double *) mxGetPr(prhs[7]);
    links_c2 = (double *) mxGetPr(prhs[8]);
    
    links_length = mxGetNumberOfElements(prhs[7]); 
    
    g=0.5*(links_length)*(links_length-1);
    
     plhs[0] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[1] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[2] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[3] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[4] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[5] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[6] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[7] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[8] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[9] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[10] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[11] = mxCreateDoubleMatrix(g,1,mxREAL);
     
     plhs[12] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[13] = mxCreateDoubleMatrix(g,1,mxREAL);     
     plhs[14] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[15] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[16] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[17] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[18] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[19] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[20] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[21] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[22] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[23] = mxCreateDoubleMatrix(g,1,mxREAL);
     
     plhs[24] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[25] = mxCreateDoubleMatrix(g,1,mxREAL);    
     plhs[26] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[27] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[28] = mxCreateDoubleMatrix(g,1,mxREAL);
     plhs[29] = mxCreateDoubleMatrix(g,1,mxREAL);
     
     plhs[30] = mxCreateDoubleMatrix(1,1,mxREAL);
     
     r11x = (double *) mxGetPr(plhs[0]);
     r11y = (double *) mxGetPr(plhs[1]);
     r11z = (double *) mxGetPr(plhs[2]);
     r21x = (double *) mxGetPr(plhs[3]);
     r21y = (double *) mxGetPr(plhs[4]);
     r21z = (double *) mxGetPr(plhs[5]);
     r12x = (double *) mxGetPr(plhs[6]);
     r12y = (double *) mxGetPr(plhs[7]);
     r12z = (double *) mxGetPr(plhs[8]);
     r22x = (double *) mxGetPr(plhs[9]);
     r22y = (double *) mxGetPr(plhs[10]);
     r22z = (double *) mxGetPr(plhs[11]);

     v11x = (double *) mxGetPr(plhs[12]);
     v11y = (double *) mxGetPr(plhs[13]);
     v11z = (double *) mxGetPr(plhs[14]);
     v21x = (double *) mxGetPr(plhs[15]);
     v21y = (double *) mxGetPr(plhs[16]);
     v21z = (double *) mxGetPr(plhs[17]);
     v12x = (double *) mxGetPr(plhs[18]);
     v12y = (double *) mxGetPr(plhs[19]);
     v12z = (double *) mxGetPr(plhs[20]);
     v22x = (double *) mxGetPr(plhs[21]);
     v22y = (double *) mxGetPr(plhs[22]);
     v22z = (double *) mxGetPr(plhs[23]);
     
     in1s1 = (double *) mxGetPr(plhs[24]);
     in2s1 = (double *) mxGetPr(plhs[25]);
     in1s2 = (double *) mxGetPr(plhs[26]);
     in2s2 = (double *) mxGetPr(plhs[27]);
     ii = (double *) mxGetPr(plhs[28]);
     ij = (double *) mxGetPr(plhs[29]);
     
     count = (double *) mxGetPr(plhs[30]);
     
    CreateInputMex (rn_x,rn_y,rn_z,v_x,v_y,v_z,flag,links_c1,links_c2,links_length,count,
            r11x,r11y,r11z,r21x,r21y,r21z,r12x,r12y,r12z,r22x,r22y,r22z,
            v11x,v11y,v11z,v21x,v21y,v21z,v12x,v12y,v12z,v22x,v22y,v22z,
            in1s1, in2s1, in1s2, in2s2, ii, ij);
}

/**************************************************************************/

void CreateInputMex (double *rn_x,double *rn_y,double *rn_z, double *v_x,double *v_y,double *v_z, double *flag,
       double *links_c1, double *links_c2,int links_length,double *count,
       double *r11x,double *r11y,double *r11z, double *r21x,double *r21y,double *r21z, double *r12x,double *r12y,double *r12z, double *r22x,double *r22y,double *r22z,
       double *v11x,double *v11y,double *v11z, double *v21x,double *v21y,double *v21z, double *v12x,double *v12y,double *v12z, double *v22x,double *v22y,double *v22z,
       double *in1s1, double *in2s1, double *in1s2, double *in2s2, double *ii, double *ij)
{
    int i=0; 
    int j;
    int p=0;
    int n1s1_int,n2s1_int,n1s2_int,n2s2_int;
    int flag1;
    int flag2;
    n1s1_int =0;
    n2s1_int =0;
    n1s2_int =0;
    n2s2_int= 0;
    *count=0;
    
    while (i<(links_length-1))
    {
        flag1=(int)round(flag[(int)round(links_c1[i])-1]);
        flag2=(int)round(flag[(int)round(links_c2[i])-1]);
        /*printf("i=%i, flag1=%i, flag2=%i \n",i,*flag1,*flag2);*/
        if (flag1==67 || flag2==67)
        {
            i=i+1;
            continue;
        }
        j=i+1;
        /*printf("i=%i \n",i);*/
        while (j<(links_length))
        {
            n1s1_int = (int)round(links_c1[i]-1); /*correct for matlab indexing*/
            n2s1_int = (int)round(links_c2[i]-1); /*correct for matlab indexing*/
            n1s2_int = (int)round(links_c1[j]-1); /*correct for matlab indexing*/
            n2s2_int = (int)round(links_c2[j]-1); /*correct for matlab indexing*/
            /*printf("n1s1=%i, n2s1=%i, n1s2=%i, n2s2=%i \n",n1s1_int,n2s1_int,n1s2_int,n2s2_int);*/
            
            if ((n1s1_int!=n1s2_int)&&(n1s1_int!=n2s2_int)&&(n2s1_int!=n1s2_int)&&(n2s1_int!=n2s2_int))
            {
                /*printf("p=%i \n",p);*/
                r11x[p] = rn_x[n1s1_int];
                r11y[p] = rn_y[n1s1_int];
                r11z[p] = rn_z[n1s1_int];
                r21x[p] = rn_x[n2s1_int];
                r21y[p] = rn_y[n2s1_int];
                r21z[p] = rn_z[n2s1_int];
                r12x[p] = rn_x[n1s2_int];
                r12y[p] = rn_y[n1s2_int];
                r12z[p] = rn_z[n1s2_int];
                r22x[p] = rn_x[n2s2_int];
                r22y[p] = rn_y[n2s2_int];
                r22z[p] = rn_z[n2s2_int];
                
                v11x[p] = v_x[n1s1_int];
                v11y[p] = v_y[n1s1_int];
                v11z[p] = v_z[n1s1_int];
                v21x[p] = v_x[n2s1_int];
                v21y[p] = v_y[n2s1_int];
                v21z[p] = v_z[n2s1_int];
                v12x[p] = v_x[n1s2_int];
                v12y[p] = v_y[n1s2_int];
                v12z[p] = v_z[n1s2_int];
                v22x[p] = v_x[n2s2_int];
                v22y[p] = v_y[n2s2_int];
                v22z[p] = v_z[n2s2_int];
                
                /*printf("n1s1=%i, n2s1=%i, n1s2=%i, n2s2=%i \n",n1s1_int,n2s1_int,n1s2_int,n2s2_int);*/
                in1s1[p] = n1s1_int +1;
                in2s1[p] = n2s1_int +1;
                in1s2[p] = n1s2_int +1;
                in2s2[p] = n2s2_int +1;
                ii[p] = i+1;
                ij[p] = j+1; 
                
                *count=*count+1;
                p=p+1;
            }
            
         j=j+1;
        }
    i=i+1;    
    }
    
}