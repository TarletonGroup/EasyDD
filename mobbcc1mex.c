#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <mex.h>

/*#ifdef _WIN32
    // this function does rounding as MS visual studio can't do it!
    int round( double r ) {
    return (r > 0.0) ? (r + 0.5) : (r - 0.5); 
    }
#endif
*/

void cross3(double a[3], double b[3], double c[3]);
void matrixinverse(double m[3][3],double minv[3][3]);
double determinant(double m[3][3]);

/************************** MEX gateway function ***********************/

/*(fseg,rn,links,connectivity,nodelist,conlist)*/
void mexFunction(int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
    /*********************************************************************/
    /*MEX memory management*/
    /*********************************************************************/
    double *nodelist;
    double **connlist;
    double *connlist_pointer;
    double **connectivity;
    double *connectivity_pointer;
    double **links;
    double *links_pointer;
    double **rn;
    double *rn_pointer;
    double **fseg;
    double *fseg_pointer;
    double Beclimb,Bclimb,Bline,Bglide,Bscrew,Bedge;
    int nodelist_length;
    int i,j,M,N;
    
     /*variables internal to code*/
    int n,n0,numNbrs;
    int ii,linkid,posinlink,n1;
    double *fx,*fy,*fz;
    /*double *vx,*vy,*vz;*/
    double *f;
    double Btotal[3][3]={0};
    double Btotalinv[3][3]={0};
    double eye[3][3];
    double L,cth2,sth2;
    double rt[3],fsegn0[3],burgv[3];
    double linedir[3];
    double mdir[3],ndir[3];
    double ldotl,bdotb,ldotb,factor;
    double ltl[3][3],ntn[3][3],mtm[3][3];
    int x,y;
    double eps=1E-12;
    double *B11,*B22,*B33,*B12,*B13,*B23;
    
    /*arrays*/
    fseg_pointer = (double *) mxGetPr(prhs[0]);
    rn_pointer = (double *) mxGetPr(prhs[1]);
    links_pointer = (double *) mxGetPr(prhs[2]);
    connectivity_pointer = (double *) mxGetPr(prhs[3]);
    nodelist = (double *) mxGetPr(prhs[4]);
    connlist_pointer = (double *) mxGetPr(prhs[5]);
    /*B parameters*/
    Beclimb = mxGetScalar(prhs[6]);
    Bline = mxGetScalar(prhs[7]);
    Bscrew = mxGetScalar(prhs[8]);
    Bedge = mxGetScalar(prhs[9]);
    
    /*allocate memory for each 2D matrix*/
    /*fseg*/
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    fseg = mxCalloc(N,sizeof(double*));
    for (j=0;j<N;j++){
        fseg[j]=fseg_pointer+M*j;
    }
    /*rn*/
    M = mxGetM(prhs[1]);
    N = mxGetN(prhs[1]);
    rn = mxCalloc(N,sizeof(double*));
    for (j=0;j<N;j++){
        rn[j]=rn_pointer+M*j;
    }
    /*links*/
    M = mxGetM(prhs[2]);
    N = mxGetN(prhs[2]);
    links = mxCalloc(N,sizeof(double*));
    for (j=0;j<N;j++){
        links[j]=links_pointer+M*j;
    }
    /*connectivity*/
    M = mxGetM(prhs[3]);
    N = mxGetN(prhs[3]);
    connectivity = mxCalloc(N,sizeof(double*));
    for (j=0;j<N;j++){
        connectivity[j]=connectivity_pointer+M*j;
    }
    /*connlist*/
    M = mxGetM(prhs[5]);
    N = mxGetN(prhs[5]);
    connlist = mxCalloc(N,sizeof(double*));
    for (j=0;j<N;j++){
        connlist[j]=connlist_pointer+M*j;
    }
    
    /*allocate memory of output arrays (and temporary ones)*/
    M = mxGetM(prhs[4]);
    N = mxGetN(prhs[4]);
    if (N<M){
        nodelist_length=M;
    }
    else{
        nodelist_length=N;
    }
    
    plhs[0] = mxCreateDoubleMatrix(nodelist_length,1,mxREAL); /*fx*/
    plhs[1] = mxCreateDoubleMatrix(nodelist_length,1,mxREAL); /*fy*/
    plhs[2] = mxCreateDoubleMatrix(nodelist_length,1,mxREAL); /*fz*/
    fx = (double *) mxGetPr(plhs[0]);
    fy = (double *) mxGetPr(plhs[1]);
    fz = (double *) mxGetPr(plhs[2]);
    
    plhs[3] = mxCreateDoubleMatrix(nodelist_length,1,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(nodelist_length,1,mxREAL);
    plhs[5] = mxCreateDoubleMatrix(nodelist_length,1,mxREAL);
    plhs[6] = mxCreateDoubleMatrix(nodelist_length,1,mxREAL);
    plhs[7] = mxCreateDoubleMatrix(nodelist_length,1,mxREAL);
    plhs[8] = mxCreateDoubleMatrix(nodelist_length,1,mxREAL);

    B11 = (double *) mxGetPr(plhs[3]);
    B22 = (double *) mxGetPr(plhs[4]);
    B33 = (double *) mxGetPr(plhs[5]);
    B12 = (double *) mxGetPr(plhs[6]);
    B13 = (double *) mxGetPr(plhs[7]);
    B23 = (double *) mxGetPr(plhs[8]);
    
    /*initialize eye(3,3)*/
    for (x=0;x<3;x++){
        for(y=0;y<3;y++){
            if (x==y){
                eye[x][y]=1.0;
            }
            else{
                eye[x][y]=0.0;
            }
        }
    }
    
    /*********************************************************************/
    /*start main loop**/
    /*********************************************************************/
    for (n=0; n<nodelist_length; n++){
        n0=(int)round(nodelist[n])-1;
        numNbrs=(int)round(connlist[0][n]);
        /*zero Btotal matrix*/
        for (x=0;x<3;x++){
            for (y=0;y<3;y++){
                Btotal[x][y]=0.0;
            }
        }
        for (i=0;i<numNbrs;i++)
        {
            ii = (int)round(connlist[i+1][n])-1;
            linkid = (int)round(connectivity[2*ii+1][n0])-1;
            posinlink=(int)round(connectivity[2*ii+2][n0])-1;
            n1 = (int)round(links[1-posinlink][linkid])-1;
            L=0.0;
            for (j=0;j<3;j++){
                rt[j]=rn[j][n1]-rn[j][n0];
                L=L+rt[j]*rt[j];
            }
            L=sqrt(L);
            if (L>0.0){
                ldotl=0;
                bdotb=0;
                ldotb=0;
                for (j=0;j<3;j++){
                    fsegn0[j]=fseg[3*posinlink+j][linkid];
                    burgv[j] = links[2+j][linkid];
                    linedir[j] = rt[j]/L;                
                    bdotb = bdotb + burgv[j]*burgv[j];
                    ldotb = ldotb + burgv[j]*linedir[j];
                }
                fx[n] = fx[n] + fsegn0[0];
                fy[n] = fy[n] + fsegn0[1];
                fz[n] = fz[n] + fsegn0[2];
                
                /*if (fabs(burgv[0]*burgv[1]*burgv[2])<eps){
                    for (x=0;x<3;x++){
                        for (y=0;y<3;y++){
                            Btotal[x][y]=Btotal[x][y] + (0.5*L*(Beclimb*eye[x][y]+(Bline-Beclimb)*linedir[x]*linedir[y]));
                        }
                    }
                }
                else{*/
                
                cth2 = (ldotb*ldotb) / bdotb;
                sth2 = 1.0-cth2;
                for (x=0;x<3;x++){
                    for (y=0;y<3;y++){
                        Btotal[x][y] = Btotal[x][y] + (0.5*L*(Bscrew*eye[x][y]+(Bline-Bscrew)*linedir[x]*linedir[y]));
                    }
                }
                if (sth2 > eps) {
                    cross3(burgv,linedir,ndir);
                    factor = 1.0/(sqrt(bdotb*sth2));
                    for (x=0;x<3;x++){
                        ndir[x]=ndir[x]*factor;
                    }
                    cross3(ndir,linedir,mdir);
                    Bglide=1.0 / sqrt( sth2/(Bedge*Bedge) + cth2/(Bscrew*Bscrew) );
                    Bclimb=sqrt( Beclimb*Beclimb * sth2 + Bscrew*Bscrew * cth2 );
                    for (x=0;x<3;x++){
                        for (y=0;y<3;y++){
                            Btotal[x][y]=Btotal[x][y]+((0.5*L)*((Bglide-Bscrew)*(mdir[x]*mdir[y])+(Bclimb-Bscrew)*(ndir[x]*ndir[y])));
                        }
                    }
                }
                /*}*/
            }
        }
        
        /*printf("determinant=%lf \n",determinant(Btotal));
        printf("Btotal \n");
        for (x=0;x<3;x++){
            printf("%lf %lf %lf \n",Btotal[x][0],Btotal[x][1],Btotal[x][2]);
        }
        printf("\n");*/

        /*problem with matrix inverse*/
        /*why is determinant different for same matrix???*/
        /*matrixinverse(Btotal,Btotalinv);*/

        /*vx[n]=0;
        vy[n]=0;
        vz[n]=0;
        fsegn0[0]=fx[n];
        fsegn0[1]=fy[n];
        fsegn0[2]=fz[n];
        for (x=0;x<3;x++){
            vx[n]=vz[n]+Btotalinv[x][0]*fsegn0[x];
            vy[n]=vz[n]+Btotalinv[x][1]*fsegn0[x];
            vz[n]=vz[n]+Btotalinv[x][2]*fsegn0[x];
        }*/

        B11[n]=Btotal[0][0];
        B22[n]=Btotal[1][1];
        B33[n]=Btotal[2][2];
        B12[n]=Btotal[0][1];
        B13[n]=Btotal[0][2];
        B23[n]=Btotal[1][2];
    }
}

void cross3(double a[3], double b[3], double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

void matrixinverse(double m[3][3],double minv[3][3])
{
    double invdet = 1.0 / determinant(m);
    
    minv[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
    minv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
    minv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
    minv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
    minv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
    minv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
    minv[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
    minv[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
    minv[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;
}

double determinant(double m[3][3])
{
    double det;    
    
    return det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
    m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
    m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}