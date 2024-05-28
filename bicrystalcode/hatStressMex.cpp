#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <mex.h>

void hatStressMex(double *uhat, int *nc, double *x, double *y, double *z, double *D, int mx, int mz, double w, double h, double d, double *x0, double *sigma);


/************************** MEX gateway function ***********************/

void mexFunction(int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
    double *uhat;
    int *nc;
    double *x, *y, *z, *D;
    double *x0, *sigma;
    int mx,mz;
    double w,h,d;
    
    uhat = (double *) mxGetPr(prhs[0]);
    nc = (int *) mxGetPr(prhs[1]);
    x = (double *) mxGetPr(prhs[2]);
    y = (double *) mxGetPr(prhs[3]);
    z = (double *) mxGetPr(prhs[4]);
    D = (double *) mxGetPr(prhs[5]);
    
    mx = mxGetScalar(prhs[6]);
    mz = mxGetScalar(prhs[7]);
   	w = mxGetScalar(prhs[8]);
    h = mxGetScalar(prhs[9]);
    d = mxGetScalar(prhs[10]);
	x0 = (double *) mxGetPr(prhs[11]);

    plhs[0] = mxCreateDoubleMatrix(6,1,mxREAL);
    sigma = (double *) mxGetPr(plhs[0]);
    
    hatStressMex(uhat,nc,x,y,z,D,mx,mz,w,h,d,x0,sigma);
}

/**************************************************************************/

void hatStressMex(double *uhat, int *nc, double *x, double *y, double *z, double *D, int mx, int mz, double w, double h, double d, double *x0, double *sigma)
{
	int i,j,k,p;
	double xc[3];
    
	i = ceil(x0[0]/w);
	j = ceil(x0[1]/w);
	k = ceil(x0[2]/w);
    
	p = i + (k-1)*mx + (j-1)*mx*mz - 1;
    
    //std::cout << i << ' ' << j << ' ' << k << ' ' << p << std::endl;
    
    xc[0] = 0.5*( x[nc[8*p]-1] + x[nc[8*p+6]-1] );
    xc[1] = 0.5*( y[nc[8*p]-1] + y[nc[8*p+6]-1] );
    xc[2] = 0.5*( z[nc[8*p]-1] + z[nc[8*p+6]-1] );
    
    //std::cout << xc[0] << ' ' << xc[1] << ' ' << xc[2] << std::endl;
    
	double const ds1dx = 1.0/(w*0.5);
	double const ds2dy = 1.0/(h*0.5);
	double const ds3dz = 1.0/(d*0.5);
    
    //std::cout << ds1dx << ' ' << ds2dy << ' ' << ds3dz << std::endl;
    
	double const s1 = (x0[0] - xc[0])*ds1dx;
	double const s2 = (x0[1] - xc[1])*ds2dy;
	double const s3 = (x0[2] - xc[2])*ds3dz;
    
    //std::cout << s1 << ' ' << s2 << ' ' << s3 << std::endl;
    
	double const pm1[8] = {-1,1,1,-1,-1,1,1,-1};
	double const pm2[8] = {1,1,1,1,-1,-1,-1,-1};
	double const pm3[8] = {-1,-1,1,1,-1,-1,1,1};
    
	double dNds1[8], dNds2[8], dNds3[8];
	double B[6][24] = {{0}};
    
	for (int a=0; a<8; a++){
		dNds1[a] = 0.125*pm1[a]*(1+pm2[a]*s2)*(1+pm3[a]*s3);
		dNds2[a] = 0.125*(1+pm1[a]*s1)*pm2[a]*(1+pm3[a]*s3);
		dNds3[a] = 0.125*(1+pm1[a]*s1)*(1+pm2[a]*s2)*pm3[a];
        
        //std::cout << dNds1[a] << std::endl;
        //std::cout << dNds2[a] << std::endl;
        //std::cout << dNds3[a] << std::endl;
        
		//In Matlab 3*(a-1)+1 is equivalent to 3*a for C/C++
					//3*(a-1)+2 is equivalent to 3*a+1
					//3*(a-1)+3 is equivalent to 3*a+2
					//3*a 	  is equivalent to 3*(a+1)-1
        
		B[0][3*a] = dNds1[a]*ds1dx;
		B[1][3*a+1] = dNds2[a]*ds2dy;
		B[2][3*a+2] = dNds3[a]*ds3dz;
        
		B[3][3*a] = B[1][3*a+1];
		B[3][3*a+1] = B[0][3*a];
        
		B[4][3*a] = B[2][3*(a+1)-1];
		B[4][3*(a+1)-1] = B[0][3*a];
        
		B[5][3*a+1] = B[2][3*(a+1)-1];
		B[5][3*(a+1)-1] = B[1][3*a+1];
	}
    
    //for (int a=0; a<8; a++){
    //    for (int b=0; b<24; b++){
    //        std::cout << B[a][b] << ' ';
    //    }
    //    std::cout << ';' << std::endl;
    //}
    
	double U[24] = {0};
	for (int a=0; a<8; a++){
		int loc = nc[8*p+a]-1;
        //std::cout << loc << std::endl;
		U[3*a] = uhat[3*loc];
		U[3*a+1] = uhat[3*loc+1];
		U[3*a+2] = uhat[3*loc+2];
	}
    
    //for (int a=0; a<1; a++){
    //    std::cout << U[a] << std::endl;
    //}
    
	double BtimesU[6] = {0};
	for (int a=0; a<6; a++)
	{
		for (int b=0; b<24; b++)
		{
			BtimesU[a] += B[a][b]*U[b];
		}
	}

    //for (int a=0; a<6; a++){
    //    std::cout << BtimesU[a] << std::endl;
    //}
    
	for (int a=0; a<6; a++)
	{
		for (int b=0; b<6; b++){
			sigma[a] += D[6*a+b]*BtimesU[b];
		}
	}
    
    //for (int a=0; a<6; a++){
    //    std::cout << sigma[a] << std::endl;
    //}
}
