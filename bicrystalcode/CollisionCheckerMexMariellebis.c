#define _USE_MATH_DEFINES
#include <math.h>
#include <mex.h>
#include <matrix.h>

/* Old windows compilers do not have a round function. If this does not compile on your machine as
a result of undefined round function, uncomment this section C-preprosessor function.
*/
/*#ifdef _WIN32
    // this function does rounding as MS visual studio can't do it!
    int round( double r ) {
    return (r > 0.0) ? (r + 0.5) : (r - 0.5);
    }
#endif
*/


void MinDistCalc(double x0[3], double x1[3], double y0[3], double y1[3], double vx0[3], double vx1[3], double vy0[3], double vy1[3], double dist2[1], double ddist2dt[1], double L1[1], double L2[1]);

/************************** MEX gateway function ***********************/

void mexFunction(int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
    /********* variable definitions *********/
    double *rn_x, *rn_y, *rn_z;
    double *v_x, *v_y, *v_z;
    double *flag;
    double *links_c1, *links_c2;
    double *connectivity_pointer, **connectivity;
    int rn_length;
    int links_length;
    double mindist;
    double mindist2;
    double *n1s1,*n2s1,*n1s2,*n2s2,*s1, *s2;
    int n1s1_int,n2s1_int,n1s2_int,n2s2_int;
    double norm2t0, norm2t1, tmp1d, tmp2d, tmp3d;
    double *colliding_segments;
    double x0[3], x1[3], y0[3], y1[3];
    double vx0[3], vx1[3], vy0[3], vy1[3];
    double dist2[1]={0}, ddist2dt[1]={0}, L1[1]={0}, L2[1]={0};
    int logic,flag1,flag2;
    int i,j,k;
    int link_col,link_row,nodenoti,linkid,tmp;
    const double eps=1E-12;
    int connectivity_M, connectivity_N;
    double *floop;
    double * segpair;
    /********* MEX memory management *********/
    rn_x = (double *) mxGetPr(prhs[0]);
    rn_y = (double *) mxGetPr(prhs[1]);
    rn_z = (double *) mxGetPr(prhs[2]);
    flag = (double *) mxGetPr(prhs[3]);
    v_x = (double *) mxGetPr(prhs[4]);
    v_y = (double *) mxGetPr(prhs[5]);
    v_z = (double *) mxGetPr(prhs[6]);
    links_c1 = (double *) mxGetPr(prhs[7]);
    links_c2 = (double *) mxGetPr(prhs[8]);
    connectivity_pointer = (double *) mxGetPr(prhs[9]);
    mindist = mxGetScalar(prhs[10]);
    rn_length = mxGetNumberOfElements(prhs[0]);
    links_length = mxGetNumberOfElements(prhs[7]);
    /*printf("links_length=%i \n",links_length);*/
    /*create 2D array for connectivity */
    connectivity_M = mxGetM(prhs[9]);
    connectivity_N = mxGetN(prhs[9]);
    connectivity = mxCalloc(connectivity_N,sizeof(double*));
    for (j=0;j<connectivity_N;j++){
        connectivity[j]=connectivity_pointer+connectivity_M*j;
    }
    /*printf("%i x %i \n",connectivity_M,connectivity_N);
    for (i=0;i<connectivity_M;i++){
        printf("\n");
        for (j=0;j<connectivity_N;j++){
            printf("%f ",connectivity[j][i]);
        }
    }*/
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[5] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[6] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[7] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[8] = mxCreateDoubleMatrix(1,1,mxREAL);

    colliding_segments = (double *) mxGetPr(plhs[0]);
    n1s1 = (double *) mxGetPr(plhs[1]);
    n2s1 = (double *) mxGetPr(plhs[2]);
    n1s2 = (double *) mxGetPr(plhs[3]);
    n2s2 = (double *) mxGetPr(plhs[4]);
    floop = (double *) mxGetPr(plhs[5]);
    s1 = (double *) mxGetPr(plhs[6]);
    s2 = (double *) mxGetPr(plhs[7]);
    segpair = (double *) mxGetPr(plhs[8]);

    /********* collision checker routine *********/
    mindist2 = mindist*mindist;
    i=0;
    segpair[0]=0;
    colliding_segments[0] = 0;

    /*look for unconnected links*/
    while (i<(links_length-1)){
        flag1=(int)round(flag[(int)round(links_c1[i])-1]);
        flag2=(int)round(flag[(int)round(links_c2[i])-1]);
        if (flag1==67 || flag2==67)
        {
            i=i+1;
            continue;
        }
        j=i+1;
        while (j<(links_length)){
            n1s1_int = (int)round(links_c1[i])-1; /*correct for matlab indexing*/
            n2s1_int = (int)round(links_c2[i])-1; /*correct for matlab indexing*/
            n1s2_int = (int)round(links_c1[j])-1; /*correct for matlab indexing*/
            n2s2_int = (int)round(links_c2[j])-1; /*correct for matlab indexing*/
            /*printf("n1s1=%i, n2s1=%i, n1s2=%i, n2s2=%i \n",n1s1_int,n2s1_int,n1s2_int,n2s2_int);*/
            if ((n1s1_int!=n1s2_int)&&(n1s1_int!=n2s2_int)&&(n2s1_int!=n1s2_int)&&(n2s1_int!=n2s2_int)){
                /*uncomment to compare with matlab script to check n1s1 and n2s1 - checked*/
                /*printf("n1s1=%i, n2s1=%i, n1s2=%i, n2s2=%i \n",n1s1_int,n2s1_int,n1s2_int,n2s2_int);*/
                segpair[0]=segpair[0]+1;
                x0[0] = rn_x[n1s1_int];
                x0[1] = rn_y[n1s1_int];
                x0[2] = rn_z[n1s1_int];
                vx0[0] = v_x[n1s1_int];
                vx0[1] = v_y[n1s1_int];
                vx0[2] = v_z[n1s1_int];

                x1[0] = rn_x[n2s1_int];
                x1[1] = rn_y[n2s1_int];
                x1[2] = rn_z[n2s1_int];
                vx1[0] = v_x[n2s1_int];
                vx1[1] = v_y[n2s1_int];
                vx1[2] = v_z[n2s1_int];

                y0[0] = rn_x[n1s2_int];
                y0[1] = rn_y[n1s2_int];
                y0[2] = rn_z[n1s2_int];
                vy0[0] = v_x[n1s2_int];
                vy0[1] = v_y[n1s2_int];
                vy0[2] = v_z[n1s2_int];

                y1[0] = rn_x[n2s2_int];
                y1[1] = rn_y[n2s2_int];
                y1[2] = rn_z[n2s2_int];
                vy1[0] = v_x[n2s2_int];
                vy1[1] = v_y[n2s2_int];
                vy1[2] = v_z[n2s2_int];

                MinDistCalc(x0,x1,y0,y1,vx0,vx1,vy0,vy1,dist2,ddist2dt,L1,L2);

                logic = ((dist2[0]<mindist2)&&(ddist2dt[0]<-eps))||(dist2[0]<eps);
                if (logic == 1){
                    colliding_segments[0]=1;
                    n1s1[0] = (double)(n1s1_int+1); /*correct for matlab indexing*/
                    n1s2[0] = (double)(n1s2_int+1); /*correct for matlab indexing*/
                    n2s1[0] = (double)(n2s1_int+1); /*correct for matlab indexing*/
                    n2s2[0] = (double)(n2s2_int+1); /*correct for matlab indexing*/
                    s1[0] = (double) (i+1);
                    s2[0] = (double) (j+1);
                    floop[0] = 1;
                    printf("Unconnected links found... Running collision correction.");
                    /*remember to de-allocate 2D connectivity array*/
                    mxFree(connectivity);
                    return;
                }
            }
            j=j+1;
        }
        i=i+1;
    }

    i=0;
    while (i<rn_length){
        if ((int)round(flag[i])==67){
            i=i+1;
            continue;
        }

        // Bruce Bromage and Daniel Celis 22/06/2020.
        // Node i is the hinge node. We only want to collide it if it is a junction. Else it gets
        // remeshed by remesh all. This was causing problems by colliding and remeshing the same
        // segment, thus undoing any changes and preventing the evolution of the network.
        if ((int)round(connectivity[0][i] < 3)){
            i = i + 1;
            continue;
        }

        j=0;
        while (j<=(int)round(connectivity[0][i])-1 /*correct for matlab indexing*/){
            link_row = (int)round(connectivity[2*j+1][i])-1/*correct for matlab indexing*/;
            link_col = 3-(int)round(connectivity[2*j+2][i]);
            if (link_col == 1){
                nodenoti = (int)round(links_c1[link_row])-1/*correct for matlab indexing*/;
            }
            else if (link_col == 2){
                nodenoti = (int)round(links_c2[link_row])-1/*correct for matlab indexing*/;
            }
            else {
                printf("Error in accessing connectivity. See collision_checker or links array");
                nodenoti=0; /*dummy value*/
            }
            /*printf("n1s1=%i, n2s1=%i, nodenoti=%i \n",n1s1_int,n2s1_int,nodenoti);*/

            // Find the squared length of link 0.
            tmp1d = rn_x[i] - rn_x[nodenoti];
            tmp2d = rn_y[i] - rn_y[nodenoti];
            tmp3d = rn_z[i] - rn_z[nodenoti];
            norm2t0 = tmp1d*tmp1d + tmp2d*tmp2d + tmp3d*tmp3d;

            k=0;
            while (k<=(int)round(connectivity[0][i])-1 /*correct for matlab indexing*/){
                linkid = (int)round(connectivity[2*k+1][i])-1;/*correct for matlab indexing*/
                if (j!=k){
                    n1s1_int = (int)round(links_c1[linkid])-1;/*correct for matlab indexing*/
                    n2s1_int = (int)round(links_c2[linkid])-1;/*correct for matlab indexing*/

                    // Find squared length of link 1.
                    tmp1d = rn_x[n1s1_int] - rn_x[n2s1_int];
                    tmp2d = rn_y[n1s1_int] - rn_y[n2s1_int];
                    tmp3d = rn_z[n1s1_int] - rn_z[n2s1_int];
                    norm2t1 = tmp1d*tmp1d + tmp2d*tmp2d + tmp3d*tmp3d;

                    /* MinDistCalc calculates the minimum distance between two lines. In this case,
                       one of the lines is collapsed into a point. This makes it so there are two
                       possible distance calculations. By swapping the node labels, we collapse the
                       shortest segment into a point and measure the distance from this point to the
                       other segment. This is a shorter distance than if we were to collapse the
                       longer segment.
                    */
                    if (norm2t0 > norm2t1){
                        tmp = linkid;
                        linkid = link_row;
                        link_row = tmp;

                        if (n1s1_int == i){
                            // tmp = n2s1_int;
                            // n2s1_int = nodenoti;
                            // nodenoti = tmp;
                            nodenoti = n2s1_int;
                        }
                        else {
                            // tmp = n1s1_int;
                            // n1s1_int = nodenoti;
                            // nodenoti = tmp;
                            nodenoti = n1s1_int;
                        }

                        n1s1_int = (int)round(links_c1[linkid])-1;/*correct for matlab indexing*/
                        n2s1_int = (int)round(links_c2[linkid])-1;/*correct for matlab indexing*/
                    }

                    /*uncomment to compare with matlab script to check n1s1 and n2s1 - checked*/
                    /*printf("n1s1=%i, n2s1=%i, nodenoti=%i \n",n1s1_int,n2s1_int,nodenoti);*/

                    x0[0] = rn_x[n1s1_int];
                    x0[1] = rn_y[n1s1_int];
                    x0[2] = rn_z[n1s1_int];
                    vx0[0] = v_x[n1s1_int];
                    vx0[1] = v_y[n1s1_int];
                    vx0[2] = v_z[n1s1_int];

                    x1[0] = rn_x[n2s1_int];
                    x1[1] = rn_y[n2s1_int];
                    x1[2] = rn_z[n2s1_int];
                    vx1[0] = v_x[n2s1_int];
                    vx1[1] = v_y[n2s1_int];
                    vx1[2] = v_z[n2s1_int];

                    y0[0] = rn_x[nodenoti];
                    y0[1] = rn_y[nodenoti];
                    y0[2] = rn_z[nodenoti];
                    vy0[0] = v_x[nodenoti];
                    vy0[1] = v_y[nodenoti];
                    vy0[2] = v_z[nodenoti];

                    y1[0] = rn_x[nodenoti];
                    y1[1] = rn_y[nodenoti];
                    y1[2] = rn_z[nodenoti];
                    vy1[0] = v_x[nodenoti];
                    vy1[1] = v_y[nodenoti];
                    vy1[2] = v_z[nodenoti];

                    MinDistCalc(x0,x1,y0,y1,vx0,vx1,vy0,vy1,dist2,ddist2dt,L1,L2);
                    logic=(dist2[0]<mindist2)&(ddist2dt[0]<-eps);

                    if (logic == 1){
                        colliding_segments[0]=1;
                        n1s1[0] = (double)(n1s1_int+1); /*correct for matlab indexing*/
                        n1s2[0] = (double)(nodenoti+1); /*correct for matlab indexing*/
                        n2s1[0] = (double)(n2s1_int+1); /*correct for matlab indexing*/
                        n2s2[0] = (double)(nodenoti+1); /*correct for matlab indexing*/
                        s1[0] = (double) linkid + 1;
                        s2[0] = (double) link_row + 1;

                        floop[0] = 2;
                        printf("Hinge condition found... Running collision correction");
                        /*remember to de-allocate 2D connectivity array*/
                        mxFree(connectivity);
                        return;
                    }
                }
                k=k+1;
            }
            j=j+1;
        }
        i=i+1;
    }

    /*printf("n1s2_int = %i, n2s2_int = %i",n1s2_int,n2s2_int);*/
    n1s1[0] = (double)(n1s1_int+1); /*correct for matlab indexing*/
    n2s1[0] = (double)(n2s1_int+1); /*correct for matlab indexing*/
    n1s2[0] = (double)(n1s2_int+1); /*correct for matlab indexing*/
    n2s2[0] = (double)(n2s2_int+1); /*correct for matlab indexing*/
    /*remember to de-allocate 2D connectivity array*/
    mxFree(connectivity);
}

/**************************************************************************/

void MinDistCalc(double x0[3], double x1[3], double y0[3], double y1[3], double vx0[3], double vx1[3], double vy0[3], double vy1[3], double dist2[1],double ddist2dt[1], double L1[1], double L2[1])
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
