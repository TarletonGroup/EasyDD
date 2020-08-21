//
//  SegQuadTrianForces.c
//
//
//  Created by Sylvain Queyreau on 12/12/14.
//
//
//  contains: SegQuadTrianForces (arbitrary dislocation and surface orientations)
//            SpecialSegQuadTrianForces (parallel case)
//            Main (read input data from a file)
//
//
//  These routines evaluate the nodal forces associated with the traction field induced
//  by non-singular straight dislocation onto a quadratic triangular element.
//  The non-singular stress field expression for isotropic elasticity is provided
//  in [Cai et al. Journal of the Mechanics and Physics of Solids 54 (2006) 561–587].
//  The surface element is a six-node quadratic triangular element of arbitrary
//  shape. The surface element has to remain planar.
//
//  The nodal force evaluation is associated to a triple integration along the dislocation
//  and along the surface element. This can be achieved by a sequence of three integrations
//  by part that present recurrence relationships.
//
//  See Latex file for more details.
//

#include <quadmath.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <tgmath.h>

/*-------------------------------------------------------------------------
 *
 *      Function:     SegQuadTrianForces
 *      Description:  Function for calculating nodal forces induced by a
 *                    non-singular segment on a quadratic triangular surface
 *                    element. Non-parallel case.
 *
 *      Arguments:
 *          p1*,p2*      endpoints for dislocation segment beginning
 *                         at point p1 and ending at point p2. the line
 *                         vector t = p1p2
 *          p3*,p4*,p5*  endpoints for an arbitrary shaped triangular element,
 *                         3 additional mid points p6, p7, p8 are added at the
 *                         centre of each side of the triangle.
 *                         the surface normal n is defined by the cross product
 *                         n = p3p4 x p3p5.
 *          bx,by,bz     Burgers vector for segment p1p2
 *          a            core parameter
 *          MU           shear modulus
 *          NU           poisson ratio
 *          fp3* - fp8*  pointers to locations in which to return forces
 *                         on nodes p3 thru p8 respectively. Units are
 *                         defined by the units used for xi, b and Mu.
 *
 *-----------------------------------------------------------------------*/
void SegQuadTrianForces(__float128 p1x,__float128 p1y,__float128 p1z,
                        __float128 p2x,__float128 p2y,__float128 p2z,
                        __float128 p3x,__float128 p3y,__float128 p3z,
                        __float128 p4x,__float128 p4y,__float128 p4z,
                        __float128 p5x,__float128 p5y,__float128 p5z,
                        __float128 bx, __float128 by, __float128 bz,
                        __float128 MU, __float128 NU, __float128 acore,
                        __float128 *fp3x, __float128 *fp3y, __float128 *fp3z,
                        __float128 *fp4x, __float128 *fp4y, __float128 *fp4z,
                        __float128 *fp5x, __float128 *fp5y, __float128 *fp5z,
                        __float128 *fp6x, __float128 *fp6y, __float128 *fp6z,
                        __float128 *fp7x, __float128 *fp7y, __float128 *fp7z,
                        __float128 *fp8x, __float128 *fp8y, __float128 *fp8z)
{

  //
  // VARIABLE DECLARATION
  //
  int i,j,w;
  int ir1[4],ir2[4],is1[4],is2[4];
  __float128 r1, r2, s1, s2;
  __float128 r12, r22, s12, s22;
  __float128 c, d, e, f, g, c2, d2, f2, g2, e2, a2;
  __float128 ce, de, ef, fg, eg, cf, cg, onepf2p2ef, oneph2p2eh;
  __float128 x1[3], x2[3], x3[3], x4[3], x5[3];
  __float128 x1x2[3], x3x4[3], x3x5[3];
  __float128 oneoverL, alpha, tdotn, dummy;
  __float128 p[3], q[3], t[3], n[3], pxt[3], qxt[3];
  __float128 dummyvec[3], x1x3[3], x2x3[3], x1x4[3], x1x5[3];
  __float128 pdotqxt, qdotpxt, ndotx1x3, ndotx2x3;
  __float128 qxtdotx1x3, qxtdotx1x4, pxtdotx1x3, pxtdotx1x5;
  __float128 y1, y2, r13, s13;
  __float128 h, m,   h2, m2;
  __float128 cd, hm, dh, eh, em, dm;
  __float128 f2g, fg2, h2m, hm2, f3, h3, g3, m3, cde;
  __float128 dr, ds, dr_o_ds, ds_o_dr, r0, s0;
  __float128 rv[8], sv[8], yv[8], rfs[8], sfr[8];
  __float128 rv2[8], sv2[8], yv2[8];
  __float128 Ra[8], rRa[8], sRa[8], Rar1[8], Rar2[8], Ras1[8], Ras2[8];
  __float128 Rdotp[8], Rdotq[8], rRdott[8], sRdott[8];
  __float128 RaRdotp[8], RaRdotq[8], rRaRdott[8], sRaRdott[8];
  __float128 A01_s1[8], A01_s2[8], B01_r1[8], B01_r2[8];
  __float128 C01_s1[8], C01_s2[8], C01_r1[8], C01_r2[8];
  __float128 A11_s1[8], A11_s2[8], B11_r1[8], B11_r2[8];
  __float128 C11_s1[8], C11_s2[8], C11_r1[8], C11_r2[8];
  __float128 A0m1_s1[8], A0m1_s2[8], B0m1_r1[8], B0m1_r2[8];
  __float128 C0m1_s1[8], C0m1_s2[8], C0m1_r1[8], C0m1_r2[8];
  __float128 A1m1_s1[8], A1m1_s2[8], B1m1_r1[8], B1m1_r2[8];
  __float128 C1m1_s1[8], C1m1_s2[8], C1m1_r1[8], C1m1_r2[8];
  __float128 AA, BB, CC, DD, EE, FF;
  __float128 vtmp1[8], vtmp2[8];
  __float128 A0m3_s1[8], A0m3_s2[8], B0m3_r1[8], B0m3_r2[8];
  __float128 C0m3_s1[8], C0m3_s2[8], C0m3_r1[8], C0m3_r2[8];
  __float128 A1m3_s1[8], A1m3_s2[8], B1m3_r1[8], B1m3_r2[8];
  __float128 A2m1_s1[8], A2m1_s2[8], B2m1_r1[8], B2m1_r2[8];
  __float128 A21_s1[8], A21_s2[8], B21_r1[8], B21_r2[8];
  __float128 A31_s1[8], A31_s2[8], B31_r1[8], B31_r2[8];
  __float128 root, tp1, tp2, tp3, tp4, tp5;
  __float128 E003_s1[8], E003_s2[8], F003_r1[8], F003_r2[8];
  __float128 D003[8], D003_r1[8], D003_r2[8], D003_s1[8], D003_s2[8];
  __float128 E103_s1[8], E103_s2[8], F103_r1[8], F103_r2[8];
  __float128 D103[8],D013[8];
  __float128 AAS001[8], BBR001[8];
  __float128 E001_s1[8], E001_s2[8], E101_s1[8], E101_s2[8], E011_s1[8], E011_s2[8];
  __float128 F001_r1[8], F001_r2[8], F101_r1[8], F101_r2[8], F011_r1[8], F011_r2[8];
  __float128 D001[8], D011[8], D101[8];
  __float128 AAS00m1[8], BBR00m1[8];
  __float128 E00m1_s1[8], E00m1_s2[8], F00m1_r1[8], F00m1_r2[8];
  __float128 D00m1[8], D00m3[8];
  __float128 E00m3_s1[8], E00m3_s2[8], F00m3_r1[8], F00m3_r2[8];
  __float128 F10m1_r1[8], F10m1_r2[8], E10m1_s1[8], E10m1_s2[8];
  __float128 F201_r1[8], F201_r2[8], E201_s1[8], E201_s2[8];
  __float128 F111_r1[8], F111_r2[8], E111_s1[8], E111_s2[8];
  __float128 AAS00m3[8], BBR00m3[8];
  __float128 D10m1[8], D01m1[8];
  __float128 AAS10m1[8], AAS01m1[8], AAS11m1[8];
  __float128 BBR01m1[8], BBR10m1[8], BBR11m1[8];
  __float128 D201[8], D021[8], D111[8], D111_1[8], D111_2[8];
  __float128 E203_s1[8], E203_s2[8], F203_r1[8], F203_r2[8];
  __float128 AAS101[8], AAS011[8];
  __float128 BBR101[8], BBR011[8];
  __float128 D203[8], D023[8], D113[8], D113_1[8], D113_2[8];
  __float128 E013_s1[8], E013_s2[8], F013_r1[8], F013_r2[8];
  __float128 E023_s1[8], E023_s2[8], F023_r1[8], F023_r2[8];
  __float128 E113_s1[8], E113_s2[8], F113_r1[8], F113_r2[8];
  __float128 E123_s1[8], E123_s2[8], F123_r1[8], F123_r2[8];
  __float128 E213_s1[8], E213_s2[8], F213_r1[8], F213_r2[8];
  __float128 AAS111[8], BBR111[8], D213[8], D123[8];
  __float128 F313_r1[8], F223_r1[8], F313_r2[8], F223_r2[8];
  __float128 E313_s1[8],   E223_s1[8], E313_s2[8], E223_s2[8];
  __float128 AAS211[8], AAS121[8], BBR211[8], BBR121[8];
  __float128 AAS201[8], AAS021[8] ,BBR201[8], BBR021[8] ;
  __float128 D313[8], D223[8], D223_1[8], D223_2[8], D133[8], D303[8], D033[8];
  __float128 (*functionPtr)(__float128);
  __float128 Dr1s1, Dr2s1, Dr1s2, Dr2s2;
  __float128 H0003[8];
  __float128 H0001[8], FFR1001[8], EES0101[8];
  __float128 H000m1[8], FFR100m1[8], EES010m1[8];
  __float128 EES000m1[8], FFR000m1[8], H0011[8], H1001[8], H0101[8];
  __float128 EES100m1[8], FFR010m1[8], H1011[8], H0111[8], H1101[8], H2001[8], H0201[8];
  __float128 H0013[8], H1003[8], H0103[8], H1013[8], H0113[8], H1103[8], H1103_1[8], H1103_2[8];
  __float128 FFR0001[8], EES0001[8], FFR0101[8], EES1001[8];
  __float128 H0023[8], H2003[8], H0203[8], H2013[8], H0213[8], H1113[8], H2103[8], H1203[8], H3003[8], H0303[8];
  __float128 H0123[8], EES0111[8], FFR0111[8], H1023[8], FFR1011[8], EES1011[8];
  __float128 EES0011[8], FFR0011[8], EES2001[8], FFR0201[8], FFR2001[8], EES0201[8], EES1101[8], FFR1101[8];
  __float128 H0015[8], H1005[8], H0105[8], H1015[8], H0115[8], H1105[8], H1105_1[8], H1105_2[8], H0025[8], H2005[8], H0205[8];
  __float128 EES0003[8], FFR0003[8], EES1003[8], FFR1003[8], EES0103[8], FFR0103[8], EES0013[8], FFR0013[8];
  __float128 H1115[8], H2105[8], H1205[8], H0125[8], H0215[8];
  __float128 H1025[8], H2015[8];
  __float128 H2025[8], H0225[8], H1125[8], H2115[8], H1215[8];
  __float128 EES1103[8], FFR1103[8], EES1013[8], FFR1013[8], EES0113[8];
  __float128 FFR0113[8], EES2013[8], FFR2013[8], EES0213[8], FFR0213[8];
  __float128 EES1113[8], FFR1113[8];
  __float128 H3015[8], H0315[8], H3005[8], H0305[8], H3025[8], H0325[8], H3115[8], H1315[8];
  __float128 H1135[8], H2125[8], H1225[8], H2215[8], H2215_1[8], H2215_2[8], H4015[8], H0415[8];
  __float128 EES2003[8], FFR0203[8], EES0203[8], FFR2003[8], EES3013[8], FFR0313[8], FFR3013[8], EES0313[8];
  __float128 EES1123[8], FFR1123[8], EES2113[8], FFR1213[8], FFR2113[8], EES1213[8];
  __float128 signv[8];
  __float128 sch[0]=0.0q, sch[1]=0.0q, sch[2]=0.0q, sch[3]=0.0q;
  __float128 sch[4]=0.0q, sch[5]=0.0q, sch[6]=0.0q, sch[7]=0.0q;
  __float128 sch[8]=0.0q, sch[9]=0.0q, sch[10]=0.0q, sch[11]=0.0q;
  __float128 sch[13]=0.0q, sch[14]=0.0q, sch[15]=0.0q, sch[16]=0.0q;
  __float128 sch[17]=0.0q, sch[18]=0.0q, sch[19]=0.0q, sch[20]=0.0q;
  __float128 sch[21]=0.0q, sch[22]=0.0q, sch[23]=0.0q, sch[24]=0.0q;
  __float128 sch[25]=0.0q, sch[26]=0.0q, sch[27]=0.0q, sch[28]=0.0q;
  __float128 sch[29]=0.0q, sch[30]=0.0q, sch[31]=0.0q, sch[32]=0.0q;
  __float128 sch[33]=0.0q, sch[34]=0.0q, sch[12]=0.0q;
  __float128 sch[35]=0.0q,   sch[36]=0.0q;
  __float128 sch[37]=0.0q,   sch[38]=0.0q,   sch[39]=0.0q,   sch[40]=0.0q;
  __float128 sch[41]=0.0q,   sch[42]=0.0q,   sch[43]=0.0q,   sch[44]=0.0q;
  __float128 sch[45]=0.0q,   sch[46]=0.0q,   sch[47]=0.0q;
  __float128 txb[3], pxb[3], qxb[3], bxt[3];
  __float128 ttbn[3], tbtn[3], tpbn[3], tqbn[3];
  __float128 txbdotn, bxtdotn, pxbdotn, qxbdotn, factor;
  __float128 b[3];
  __float128 I0013[3], I0103[3], I1003[3], I0015[3], I0105[3], I1005[3], I0215[3], I2015[3], I0125[3], I1025[3], I1115[3];
  __float128 pxbdott, qxbdott;
  __float128 F0013, F0103, F1003, F0015;
  __float128 F0105, F1005, F0215, F2015;
  __float128 F0125, F1025, F1115;
  __float128 fLLprime[3], fx4[3], fx5[3], fx6[3], fx7[3], fx8[3], fx3[3], ftot[3];
  char csca[256];
  char cv3[3][256];
  char cv8[8][256];

  //
  // VARIABLE INITIALIZATION
  //
  for(i=0;i<8;i++) {
    A01_s1[i] = 0.0q; A01_s2[i] = 0.0q;
    B01_r1[i] = 0.0q; B01_r2[i] = 0.0q;
    C01_s1[i] = 0.0q; C01_s2[i] = 0.0q;
    C01_r1[i] = 0.0q; C01_r2[i] = 0.0q;
    A11_s1[i] = 0.0q; A11_s2[i] = 0.0q;
    B11_r1[i] = 0.0q; B11_r2[i] = 0.0q;
    C11_s1[i] = 0.0q; C11_s2[i] = 0.0q;
    C11_r1[i] = 0.0q; C11_r2[i] = 0.0q;
    A0m1_s1[i] = 0.0q; A0m1_s2[i] = 0.0q;
    B0m1_r1[i] = 0.0q; B0m1_r2[i] = 0.0q;
    C0m1_s1[i] = 0.0q; C0m1_s2[i] = 0.0q;
    C0m1_r1[i] = 0.0q; C0m1_r2[i] = 0.0q;
    A1m1_s1[i] = 0.0q; A1m1_s2[i] = 0.0q;
    B1m1_r1[i] = 0.0q; B1m1_r2[i] = 0.0q;
    C1m1_s1[i] = 0.0q; C1m1_s2[i] = 0.0q;
    C1m1_r1[i] = 0.0q; C1m1_r2[i] = 0.0q;
    A0m3_s1[i] = 0.0q; A0m3_s2[i] = 0.0q;
    B0m3_r1[i] = 0.0q; B0m3_r2[i] = 0.0q;
    C0m3_s1[i] = 0.0q; C0m3_s2[i] = 0.0q;
    C0m3_r1[i] = 0.0q; C0m3_r2[i] = 0.0q;
    A1m3_s1[i] = 0.0q; A1m3_s2[i] = 0.0q;
    B1m3_r1[i] = 0.0q; B1m3_r2[i] = 0.0q;
    A2m1_s1[i] = 0.0q; A2m1_s2[i] = 0.0q;
    B2m1_r1[i] = 0.0q; B2m1_r2[i] = 0.0q;
    A21_s1[i] = 0.0q; A21_s2[i] = 0.0q;
    B21_r1[i] = 0.0q; B21_r2[i] = 0.0q;
    A31_s1[i] = 0.0q; A31_s2[i] = 0.0q;
    B31_r1[i] = 0.0q; B31_r2[i] = 0.0q;
    E003_s1[i] = 0.0q; E003_s2[i] = 0.0q;
    F003_r1[i] = 0.0q; F003_r2[i] = 0.0q;
    D003[i] = 0.0q;
    D003_r1[i] = 0.0q; D003_r2[i] = 0.0q;
    D003_s1[i] = 0.0q; D003_s2[i] = 0.0q;
    E103_s1[i] = 0.0q; E103_s2[i] = 0.0q;
    F103_r1[i] = 0.0q; F103_r2[i] = 0.0q;
    E001_s1[i] = 0.0q; E001_s2[i] = 0.0q;
    E101_s1[i] = 0.0q; E101_s2[i] = 0.0q;
    E011_s1[i] = 0.0q; E011_s2[i] = 0.0q;
    F001_r1[i] = 0.0q; F001_r2[i] = 0.0q;
    F101_r1[i] = 0.0q; F101_r2[i] = 0.0q;
    F011_r1[i] = 0.0q; F011_r2[i] = 0.0q;
    D001[i] = 0.0q; D101[i] = 0.0q; D011[i] = 0.0q;
    AAS001[i] = 0.0q; BBR001[i] = 0.0q;
    AAS00m1[i] = 0.0q; BBR00m1[i] = 0.0q;
    E00m1_s1[i] = 0.0q; E00m1_s2[i] = 0.0q;
    F00m1_r1[i] = 0.0q; F00m1_r2[i] = 0.0q;
    D00m1[i] = 0.0q; D00m3[i] = 0.0q;
    E00m3_s1[i] = 0.0q; E00m3_s2[i] = 0.0q;
    F00m3_r1[i] = 0.0q; F00m3_r2[i] = 0.0q;
    E10m1_s1[i] = 0.0q; E10m1_s2[i] = 0.0q;
    F10m1_r1[i] = 0.0q; F10m1_r2[i] = 0.0q;
    E201_s1[i] = 0.0q; E201_s2[i] = 0.0q;
    F201_r1[i] = 0.0q; F201_r2[i] = 0.0q;
    E111_s1[i] = 0.0q; E111_s2[i] = 0.0q;
    F111_r1[i] = 0.0q; F111_r2[i] = 0.0q;
    AAS00m3[i] = 0.0q; BBR00m3[i] = 0.0q;
    D10m1[i] = 0.0q; D01m1[i] = 0.0q;
    AAS10m1[i] = 0.0q; AAS01m1[i] = 0.0q; AAS11m1[i] = 0.0q;
    BBR01m1[i] = 0.0q; BBR10m1[i] = 0.0q; BBR11m1[i] = 0.0q;
    D201[i] = 0.0q; D021[i] = 0.0q;
    D111[i] = 0.0q; D111_1[i] = 0.0q; D111_2[i] = 0.0q;
    E203_s1[i] = 0.0q; E203_s2[i] = 0.0q;
    F203_r1[i] = 0.0q; F203_r2[i] = 0.0q;
    AAS101[i] = 0.0q; AAS011[i] = 0.0q;
    BBR101[i] = 0.0q; BBR011[i] = 0.0q;
    D203[i] = 0.0q; D023[i] = 0.0q;
    D113[i] = 0.0q; D113_1[i] = 0.0q; D113_2[i] = 0.0q;
    E013_s1[i] = 0.0q; E013_s2[i] = 0.0q;
    F013_r1[i] = 0.0q; F013_r2[i] = 0.0q;
    E023_s1[i] = 0.0q; E023_s2[i] = 0.0q;
    F023_r1[i] = 0.0q; F023_r2[i] = 0.0q;
    E113_s1[i] = 0.0q; E113_s2[i] = 0.0q;
    F113_r1[i] = 0.0q; F113_r2[i] = 0.0q;
    E123_s1[i] = 0.0q; E123_s2[i] = 0.0q; F123_r1[i] = 0.0q; F123_r2[i] = 0.0q;
    E213_s1[i] = 0.0q; E213_s2[i] = 0.0q; F213_r1[i] = 0.0q; F213_r2[i] = 0.0q;
    AAS111[i] = 0.0q; BBR111[i] = 0.0q; D213[i] = 0.0q; D123[i] = 0.0q;
    F313_r1[i] = 0.0q; F223_r1[i] = 0.0q; F313_r2[i] = 0.0q; F223_r2[i] = 0.0q;
    E313_s1[i] = 0.0q; E223_s1[i] = 0.0q; E313_s2[i] = 0.0q; E223_s2[i] = 0.0q;

    AAS211[i] = 0.0q; AAS121[i] = 0.0q; BBR211[i] = 0.0q; BBR121[i] = 0.0q;
    AAS201[i] = 0.0q; AAS021[i] = 0.0q; BBR201[i] = 0.0q; BBR021[i] = 0.0q;
    D313[i] = 0.0q; D133[i] = 0.0q; D303[i] = 0.0q; D033[i] = 0.0q;
    D223[i] = 0.0q; D223_1[i] = 0.0q; D223_2[i] = 0.0q;

    H0003[i] = 0.0q; H0001[i] = 0.0q; EES0101[i] = 0.0q; FFR1001[i] = 0.0q;
    H000m1[i] = 0.0q; EES010m1[i] = 0.0q; FFR100m1[i] = 0.0q;

    EES000m1[i] = 0.0q; FFR000m1[i] = 0.0q; H0011[i] = 0.0q; H1001[i] = 0.0q;
    H0101[i] = 0.0q; EES100m1[i] = 0.0q; FFR010m1[i] = 0.0q; H1011[i] = 0.0q;
    H0111[i] = 0.0q; H1101[i] = 0.0q; H2001[i] = 0.0q; H0201[i] = 0.0q;

    H0013[i] = 0.0q; H1003[i] = 0.0q; H0103[i] = 0.0q; H1013[i] = 0.0q;
    H0113[i] = 0.0q; H1103[i] = 0.0q; H1103_1[i] = 0.0q; H1103_2[i] = 0.0q;
    FFR0001[i] = 0.0q; EES0001[i] = 0.0q; FFR0101[i] = 0.0q; EES1001[i] = 0.0q;

    H0023[i] = 0.0q; H2003[i] = 0.0q; H0203[i] = 0.0q; H2013[i] = 0.0q;
    H0213[i] = 0.0q; H1113[i] = 0.0q; H2103[i] = 0.0q; H1203[i] = 0.0q;
    H3003[i] = 0.0q; H0303[i] = 0.0q;
    H0123[i] = 0.0q; EES0111[i] = 0.0q; FFR0111[i] = 0.0q;
    H1023[i] = 0.0q; FFR1011[i] = 0.0q; EES1011[i] = 0.0q;
    EES0011[i] = 0.0q; FFR0011[i] = 0.0q; EES2001[i] = 0.0q; FFR0201[i] = 0.0q;
    FFR2001[i] = 0.0q; EES0201[i] = 0.0q; EES1101[i] = 0.0q; FFR1101[i] = 0.0q;

    H0015[i] = 0.0q; H1005[i] = 0.0q; H0105[i] = 0.0q;
    H1015[i] = 0.0q; H0115[i] = 0.0q; H1105[i] = 0.0q; H1105_1[i] = 0.0q; H1105_2[i] = 0.0q;
    H0025[i] = 0.0q; H2005[i] = 0.0q; H0205[i] = 0.0q;
    EES0003[i] = 0.0q; FFR0003[i] = 0.0q; EES1003[i] = 0.0q;
    FFR1003[i] = 0.0q; EES0103[i] = 0.0q; FFR0103[i] = 0.0q;
    EES0013[i] = 0.0q; FFR0013[i] = 0.0q;

    H1115[i] = 0.0q; H2105[i] = 0.0q; H1205[i] = 0.0q; H0125[i] = 0.0q; H0215[i] = 0.0q;
    H2025[i] = 0.0q; H0225[i] = 0.0q; H1125[i] = 0.0q; H2115[i] = 0.0q; H1215[i] = 0.0q;
    H1025[i] = 0.0q; H2015[i] = 0.0q;
    EES1103[i] = 0.0q; FFR1103[i] = 0.0q; EES1013[i] = 0.0q; FFR1013[i] = 0.0q;
    EES0113[i] = 0.0q; FFR0113[i] = 0.0q; EES2013[i] = 0.0q; FFR2013[i] = 0.0q;
    EES0213[i] = 0.0q; FFR0213[i] = 0.0q; EES1113[i] = 0.0q; FFR1113[i] = 0.0q;

    H3015[i] = 0.0q; H0315[i] = 0.0q; H3005[i] = 0.0q; H0305[i] = 0.0q;
    H3025[i] = 0.0q; H0325[i] = 0.0q; H3115[i] = 0.0q; H1315[i] = 0.0q;
    H1135[i] = 0.0q; H2125[i] = 0.0q; H1225[i] = 0.0q; H2215[i] = 0.0q; H2215_1[i] = 0.0q; H2215_2[i] = 0.0q;
    H4015[i] = 0.0q; H0415[i] = 0.0q;
    EES2003[i] = 0.0q; FFR0203[i] = 0.0q; EES0203[i] = 0.0q; FFR2003[i] = 0.0q;
    EES3013[i] = 0.0q; FFR0313[i] = 0.0q; FFR3013[i] = 0.0q; EES0313[i] = 0.0q;
    EES1123[i] = 0.0q; FFR1123[i] = 0.0q; EES2113[i] = 0.0q; FFR1213[i] = 0.0q;
    FFR2113[i] = 0.0q; EES1213[i] = 0.0q;

  }
/*  TAKEN CARE OF THIS
  *fp3x = 0.0q; *fp3y = 0.0q; *fp3z = 0.0q;
  *fp4x = 0.0q; *fp4y = 0.0q; *fp4z = 0.0q;
  *fp5x = 0.0q; *fp5y = 0.0q; *fp5z = 0.0q;
  *fp6x = 0.0q; *fp6y = 0.0q; *fp6z = 0.0q;
  *fp7x = 0.0q; *fp7y = 0.0q; *fp7z = 0.0q;
  *fp8x = 0.0q; *fp8y = 0.0q; *fp8z = 0.0q;

  x1[0] = p1x; x1[1] = p1y; x1[2] = p1z;
  x2[0] = p2x; x2[1] = p2y; x2[2] = p2z;
  x3[0] = p3x; x3[1] = p3y; x3[2] = p3z;
  x4[0] = p4x; x4[1] = p4y; x4[2] = p4z;
  x5[0] = p5x; x5[1] = p5y; x5[2] = p5z;
*/
  //
  // DEFINITION OF GEOMETRIC VARIABLES
  //
  printf("Beginning of SegQuadTrianForces Function \n");
  /*  TAKEN CARE OF THIS
  // Vectors connecting the differents nodes
  for(i=0;i<3;i++) {
    x1x2[i] = x2[i] - x1[i];
    x3x4[i] = x4[i] - x3[i];
    x3x5[i] = x5[i] - x3[i];
    x1x3[i] = x3[i] - x1[i];
    x1x4[i] = x4[i] - x1[i];
    x1x5[i] = x5[i] - x1[i];
    x2x3[i] = x3[i] - x2[i];
  }
  */

  /*
   *  The vector Ra spanning the dislocation and the surface is decomposed as follows:
   *  R = y*t +r*p +s*q
   *  y is scalar spanning the dislocation, and r and s span the surface element
   *  The tripple integration required to express the nodal forces are actually
   *  done along y, r and s
   */
   /* TAKEN CARE OF THIS
  quadmath_snprintf(cv3[0], 256,  "%*.30Qf", bx); quadmath_snprintf(cv3[1], 256,  "%*.30Qf", by); quadmath_snprintf(cv3[2], 256,  "%*.30Qf", bz);
  printf("b = %s %s %s \n", cv3[0], cv3[1], cv3[2] );

  oneoverL = 1.0q / sqrtq(x3x4[0]*x3x4[0] +x3x4[1]*x3x4[1] +x3x4[2]*x3x4[2]);
  quadmath_snprintf(cv3[0], 256,  "%*.30Qf", oneoverL);
  p[0]=x3x4[0]*oneoverL;
  p[1]=x3x4[1]*oneoverL;
  p[2]=x3x4[2]*oneoverL;

  oneoverL = 1.0q / sqrtq(x3x5[0]*x3x5[0] +x3x5[1]*x3x5[1] +x3x5[2]*x3x5[2]);

  q[0]=x3x5[0]*oneoverL;
  q[1]=x3x5[1]*oneoverL;
  q[2]=x3x5[2]*oneoverL;

  oneoverL = 1.0q / sqrtq(x1x2[0]*x1x2[0] +x1x2[1]*x1x2[1] +x1x2[2]*x1x2[2]);

  t[0]=x1x2[0]*oneoverL;
  t[1]=x1x2[1]*oneoverL;
  t[2]=x1x2[2]*oneoverL;
  // Important Dot products
  c = t[0]*p[0] +t[1]*p[1] +t[2]*p[2];
  d = t[0]*q[0] +t[1]*q[1] +t[2]*q[2];
  e = q[0]*p[0] +q[1]*p[1] +q[2]*p[2];

  alpha = sinq(acosq(e));

  //
  // DEFINITION OF THE SURFACE NORMAL
  //
  n[0] = p[1]*q[2] - p[2]*q[1];
  n[1] = p[2]*q[0] - p[0]*q[2];
  n[2] = p[0]*q[1] - p[1]*q[0];

  oneoverL = 1.0q / sqrtq(n[0]*n[0] +n[1]*n[1] +n[2]*n[2]);

  n[0] = n[0]*oneoverL;
  n[1] = n[1]*oneoverL;
  n[2] = n[2]*oneoverL;

  // Important cross products
  pxt[0] = p[1]*t[2] - p[2]*t[1];
  pxt[1] = p[2]*t[0] - p[0]*t[2];
  pxt[2] = p[0]*t[1] - p[1]*t[0];

  qxt[0] = q[1]*t[2] - q[2]*t[1];
  qxt[1] = q[2]*t[0] - q[0]*t[2];
  qxt[2] = q[0]*t[1] - q[1]*t[0];
  */
  /* TAKEN CARE OF THIS
  tdotn = t[0]*n[0] +t[1]*n[1] +t[2]*n[2];
  pdotqxt = p[0]*qxt[0] +p[1]*qxt[1] +p[2]*qxt[2];
  qdotpxt = q[0]*pxt[0] +q[1]*pxt[1] +q[2]*pxt[2];

  ndotx1x3 = n[0]*x1x3[0] +n[1]*x1x3[1] +n[2]*x1x3[2];
  ndotx2x3 = n[0]*x2x3[0] +n[1]*x2x3[1] +n[2]*x2x3[2];

  qxtdotx1x3 = qxt[0]*x1x3[0] +qxt[1]*x1x3[1] +qxt[2]*x1x3[2];
  qxtdotx1x4 = qxt[0]*x1x4[0] +qxt[1]*x1x4[1] +qxt[2]*x1x4[2];

  pxtdotx1x3 = pxt[0]*x1x3[0] +pxt[1]*x1x3[1] +pxt[2]*x1x3[2];
  pxtdotx1x5 = pxt[0]*x1x5[0] +pxt[1]*x1x5[1] +pxt[2]*x1x5[2];
  */
  //
  // INTEGRAL BOUNDS
  //
  // correspond to min and max values of variables r,s,y
  /* TAKEN CARE OF THIS
  y1 =   ndotx1x3 / tdotn;
  y2 =   ndotx2x3 / tdotn;
  r1 = qxtdotx1x3 / pdotqxt;
  r2 = qxtdotx1x4 / pdotqxt;
  s1 = pxtdotx1x3 / qdotpxt;
  s2 = pxtdotx1x5 / qdotpxt;
  */
  /*
  r12 = r1*r1;
  s12 = s1*s1;
  */
  r22 = r2*r2;
  s22 = s2*s2;
  r13 = r1*r12;
  s13 = s1*s12;

  /* Constants built from previous dot and cross products */
  /* TAKEN CARE OF THIS
  a2 = acore*acore;
  e2 = e*e;
  */
  c2 = c*c;
  d2 = d*d;
  /* TAKEN CARE OF THIS
  ds=s2-s1;
  dr=r2-r1;
  dr_o_ds=dr*one_o_ds;
  ds_o_dr=ds*one_o_dr;
  */
  a2 = acore*acore;
  c2 = c*c;
  d2 = d*d;
  e2 = e*e;
  ds=s2-s1;
  dr=r2-r1;
  dr_o_ds=dr/ds;
  ds_o_dr=ds/dr;
  s0=s2+ds_o_dr*r1;
  r0=r2+dr_o_ds*s1;
  f=-dr_o_ds;
  g=r0;
  h=-ds_o_dr;
  m=s0;
  h2=h*h; f2=f*f;
  g2=g*g; m2=m*m;
  eh=e*h; ef=e*f;
  dh=h*d;
  cf=f*c;
  cd=c*d;
  cg=c*g;
  fg=f*g;
  eg=e*g;
  em=e*m;
  hm=h*m;
  dm=d*m;
  de=d*e;
  ce=c*e;
  f3=f2*f;
  h3=h2*h;
  m3=m2*m;
  g3=g2*g;
  h2m=h2*m;
  hm2=h*m2;
  f2g=f2*g;
  fg2=f*g2;
  cde= cd*e;
  onepf2p2ef = 1.0q+f2+2.0q*ef;
  oneph2p2eh = 1.0q+h2+2.0q*eh;
  /*
  cde= cd*e;
  onepf2p2ef = 1.0q+f2+2.0q*ef;
  oneph2p2eh = 1.0q+h2+2.0q*eh;
  */

  /*
   *  Integrals are functions of y,r and s. But ultimately, the nodal force evaluation
   *  will be the sum of evaluations of the antiderivative for the various bounds.
   *  It is therefore more convenient to organize integrals and variables in the shape of
   *  vectors of 8 components. The different components correspond to permutation
   * of the 2 bounds per variables: r1, r2, s1, s2, y1 and y2.
   */
   /* TAKEN CARE OF THIS
  rv[0]=r2; rv[1]=r2; rv[2]=r2; rv[3]=r2; rv[4]=r1; rv[5]=r1; rv[6]=r1; rv[7]=r1;
  sv[0]=s2; sv[1]=s2; sv[2]=s1; sv[3]=s1; sv[4]=s2; sv[5]=s2; sv[6]=s1; sv[7]=s1;
  yv[0]=y2; yv[1]=y1; yv[2]=y2; yv[3]=y1; yv[4]=y2; yv[5]=y1; yv[6]=y2; yv[7]=y1;
  // in the case of triangular element and depending on the integration order,
  // the upper bound of r can be a function of s. These vectors expresse this dependence.
  rfs[0]=f*s2+r0; rfs[1]=f*s2+r0; rfs[2]=f*s1+r0; rfs[3]=f*s1+r0; rfs[4]=r1; rfs[5]=r1; rfs[6]=r1; rfs[7]=r1;
  sfr[0]=h*r2+s0; sfr[1]=h*r2+s0; sfr[2]=s1; sfr[3]=s1; sfr[4]=h*r1+s0; sfr[5]=h*r1+s0; sfr[6]=s1; sfr[7]=s1;
  */
  /*
  for(i=0;i<8;i++) {
    rv2[i] = rv[i]*rv[i];
    sv2[i] = sv[i]*sv[i];
    yv2[i] = yv[i]*yv[i];
  }
  */
/*
  // Ra = sqrt(R.R)
  for(i=0;i<8;i++) {
    Ra[i]=sqrtq(yv2[i]+rv2[i]+sv2[i]+2.0q*c*rv[i]*yv[i]+2.0q*d*sv[i]*yv[i]+2.0q*e*rv[i]*sv[i]+a2);
    rRa[i]=sqrtq(yv2[i]+rv2[i]+sfr[i]*sfr[i]+2.0q*c*rv[i]*yv[i]+2.0q*d*sfr[i]*yv[i]+2.0q*e*rv[i]*sfr[i]+a2);
    sRa[i]=sqrtq(yv2[i]+rfs[i]*rfs[i]+sv2[i]+2.0q*c*rfs[i]*yv[i]+2.0q*d*sv[i]*yv[i]+2.0q*e*rfs[i]*sv[i]+a2);
    Rar1[i]=sqrtq(yv2[i]+r1*r1+sv2[i]+2.0q*c*r1*yv[i]+2.0q*d*sv[i]*yv[i]+2.0q*e*r1*sv[i]+a2);
    Rar2[i]=sqrtq(yv2[i]+(1.0q+f2+2.0q*e*f)*sv2[i]+2.0q*sv[i]*yv[i]*(cf+d)+2.0q*sv[i]*(fg+eg)+2.0q*cg*yv[i]+g2+a2);
    Ras1[i]=sqrtq(yv2[i]+rv2[i]+s1*s1+2.0q*c*rv[i]*yv[i]+2.0q*d*s1*yv[i]+2.0q*e*rv[i]*s1+a2);
    Ras2[i]=sqrtq(yv2[i]+(1.0q+h2+2.0q*e*h)*rv2[i]+2.0q*rv[i]*yv[i]*(dh+c)+2.0q*rv[i]*(hm+em)+2.0q*dm*yv[i]+m2+a2);
  }
*/
  for(i=0;i<8;i++) {
    Rdotp[i] = rv[i]+c*yv[i]+e*sv[i];
    Rdotq[i] = sv[i]+d*yv[i]+e*rv[i];
    rRdott[i] = yv[i]+c*rv[i]+d*sfr[i];
    sRdott[i] = yv[i]+c*rfs[i]+d*sv[i];
    RaRdotp[i] = Ra[i] +Rdotp[i];
    RaRdotq[i] = Ra[i] +Rdotq[i];
    rRaRdott[i] = rRa[i] +rRdott[i];
    sRaRdott[i] = sRa[i] +sRdott[i];
  }

  //
  // LINEAR INTEGRALS
  //
  /* The following shorthand is introduced to simplify notation:
     Ail = /int r^i Ra^-l dr
     Bjl = /int s^j Ra^-l ds
     Ckl = /int y^k Ra^-l dy
     Linear integrals are the last object created by the sequence of
     integration by part. since r2 depends upon s and s2 depends up on
     r, 2 different sollutions must be provided for Ail_s1, Ail_s2  and
     Bjl_r1 and Bjl_r2 and 4 different solutions for Ckl.
  */

  /* Since integration paths are different for r1 and r2,
     and for s1 and s2, we need to evaluate integrals Ail_r1
     at vector components corresponding to r1 */
  ir1[0]=4;ir1[1]=5;ir1[2]=6;ir1[3]=7;
  is1[0]=2;is1[1]=3;is1[2]=6;is1[3]=7;
  ir2[0]=0;ir2[1]=1;ir2[2]=2;ir2[3]=3;
  is2[0]=0;is2[1]=1;is2[2]=4;is2[3]=5;

  //
  // LINEAR SEEDS ABC01
  //
  // printf("Linear Integrals \n");
  for(j=0;j<4;j++) {
    i=ir1[j];
    B01_r1[i]= logq( fabsq(RaRdotq[i]));
    C01_r1[i]= logq( fabsq(rRaRdott[i]));
  }

  AA=(onepf2p2ef);
  for(j=0;j<4;j++) {
    i=ir2[j];
    BB=(yv[i]*(cf+d)+g*(f+e));
    CC=(yv2[i]+2.0q*cg*yv[i]+g2+a2);
    B01_r2[i]=logq(2.0q*sqrtq(AA)*sqrtq(AA*sv2[i]+2.0q*BB*sv[i]+CC)+2.0q*(AA*sv[i]+BB))/sqrtq(AA);
    BB=cg+(cf+d)*sv[i];
    CC=sqrtq(-BB*BB+sv2[i]*(onepf2p2ef)+2.0q*sv[i]*(fg+eg)+g2+a2);
    C01_r2[i]=logq(2.0q*sqrtq((yv[i]+BB)*(yv[i]+BB)+CC*CC)+2.0q*(yv[i]+BB));
  }

  for(j=0;j<4;j++) {
    i=is1[j];
    A01_s1[i]= logq(fabsq(RaRdotp[i]));
    C01_s1[i]= logq(fabsq(sRaRdott[i])) ;
  }

  AA=(oneph2p2eh);
  for(j=0;j<4;j++) {
    i=is2[j];
    BB=(yv[i]*(h*d+c)+m*(h+e));
    CC=(yv2[i]+2.0q*m*d*yv[i]+m2+a2);
    A01_s2[i]=logq(2.0q*sqrtq(AA)*sqrtq(AA*rv2[i]+2.0q*BB*rv[i]+CC)+2.0q*(AA*rv[i]+BB))/sqrtq(AA);
    BB=dm+(dh+c)*rv[i];
    CC=sqrtq(-BB*BB+rv2[i]*(oneph2p2eh)+2.0q*rv[i]*(hm+em)+m2+a2);
    C01_s2[i]=logq(2.0q*sqrtq((yv[i]+BB)*(yv[i]+BB)+CC*CC)+2.0q*(yv[i]+BB));
  }

  // ABC11
  for(j=0;j<4;j++) {
    i=ir1[j];
    B11_r1[i]=sRa[i]-(d*yv[i]+e*rv[i])*B01_r1[i];
    C11_r1[i]=Rar1[i]-(c*rv[i]+d*sfr[i])*C01_r1[i];

    i=ir2[j];
    B11_r2[i]=1.0q/(onepf2p2ef)*(sRa[i]-(yv[i]*(cf+d)+fg+eg)*B01_r2[i]);
    C11_r2[i]=Rar2[i]-(sv[i]*(cf+d)+cg)*C01_r2[i];

    i=is1[j];
    A11_s1[i]=rRa[i]-(c*yv[i]+e*sv[i])*A01_s1[i];
    C11_s1[i]=Ras1[i]-(d*sv[i]+c*rfs[i])*C01_s1[i];

    i=is2[j];
    A11_s2[i]=1.0q/(oneph2p2eh)*(rRa[i]-(yv[i]*(dh+c)+hm+em)*A01_s2[i]);
    C11_s2[i]=Ras2[i]-(rv[i]*(dh+c)+dm)*C01_s2[i];
  }

  // ABC1m1
  for(j=0;j<4;j++) {
    i=ir1[j];
    B0m1_r1[i]=0.5q*(Rdotq[i]*Ra[i] +(Ra[i]*Ra[i]-(Rdotq[i]*Rdotq[i]))*B01_r1[i]);
    B1m1_r1[i]=1.0q/3.0q*(Ra[i]*Ra[i]*Ra[i])-(d*yv[i]+e*rv[i])*B0m1_r1[i];
    C0m1_r1[i]=0.5q*(rRdott[i]*rRa[i] +(rRa[i]*rRa[i]-(rRdott[i]*rRdott[i]))*C01_r1[i]);
    C1m1_r1[i]=-1.0q/3.0q*(-sRa[i]*sRa[i]*sRa[i])-(c*rfs[i]+d*sv[i])*C0m1_r1[i];

    i=ir2[j];
    B0m1_r2[i]= 0.5q*(sv[i]*sRa[i]+(yv[i]*(cf+d)+g*(f+e))*B11_r2[i]+(yv2[i]+2.0q*g*c*yv[i]+g2+a2)*B01_r2[i]);
    B1m1_r2[i]=1.0q/(onepf2p2ef)*(1.0q/3.0q*sRa[i]*sRa[i]*sRa[i]-(yv[i]*(cf+d)+g*(f+e))*B0m1_r2[i]);
    C0m1_r2[i]=0.5q*(+Rar2[i]*(yv[i]+sv[i]*(cf+d)+cg)+((onepf2p2ef-(cf+d)*(cf+d))*sv2[i]+2.0q*sv[i]*(fg+eg-cg*(cf+d))+(1.0q-c2)*g2+a2)*C01_r2[i]);
    C1m1_r2[i]=1.0q/3.0q*Rar2[i]*Rar2[i]*Rar2[i]-(sv[i]*(cf+d)+cg)*C0m1_r2[i];

    i=is1[j];
    A0m1_s1[i]=0.5q*(Rdotp[i]*Ra[i] +(Ra[i]*Ra[i]-(Rdotp[i]*Rdotp[i]))*A01_s1[i]);
    A1m1_s1[i]=-1.0q/3.0q*(-Ras1[i]*Ras1[i]*Ras1[i])-(c*yv[i]+e*sv[i])*A0m1_s1[i];
    C0m1_s1[i]=0.5q*(sRdott[i]*sRa[i] +(sRa[i]*Ra[i]-(sRdott[i]*sRdott[i]))*C01_s1[i]);
    C1m1_s1[i]=-1.0q/3.0q*(-rRa[i]*rRa[i]*rRa[i])-(c*rv[i]+d*sfr[i])*C0m1_s1[i];

    i=is2[j];
    A0m1_s2[i]= 0.5q*(rv[i]*rRa[i]+(yv[i]*(dh+c)+m*(h+e))*A11_s2[i]+(yv2[i]+2.0q*m*d*yv[i]+m2+a2)*A01_s2[i]);
    A1m1_s2[i]=1.0q/(oneph2p2eh)*(1.0q/3.0q*rRa[i]*rRa[i]*rRa[i]-(yv[i]*(dh+c)+m*(h+e))*A0m1_s2[i]);
    C0m1_s2[i]=0.5q*(+Ras2[i]*(yv[i]+rv[i]*(dh+c)+dm)+((oneph2p2eh-(dh+c)*(dh+c))*rv2[i]+2.0q*rv[i]*(hm+em-dm*(dh+c))+(1.0q-d2)*m2+a2)*C01_s2[i]);
    C1m1_s2[i]=1.0q/3.0q*Ras2[i]*Ras2[i]*Ras2[i]-(rv[i]*(dh+c)+dm)*C0m1_s2[i];
  }

  // ABC0m3
  for(j=0;j<4;j++) {
    i=ir1[j];
    B0m3_r1[i]=1.0q/4.0q*(Rdotq[i]*powq(Rar1[i],3.0q) +3.0q*(powq(Rar1[i],2.0q)-(powq(Rdotq[i],2.0q)))*B0m1_r1[i]);
    C0m3_r1[i]=1.0q/4.0q*(rRdott[i]*powq(Rar1[i],3.0q) +3.0q*(powq(Rar1[i],2.0q)-powq(rRdott[i],2.0q))*C0m1_r1[i]);

    i=ir2[j];
    vtmp1[i]=((cf+d)*yv[i]+fg+eg)/(onepf2p2ef);
    vtmp2[i]= powq(((cf+d)*yv[i]+fg+eg),2.0q)/(onepf2p2ef);
    B0m3_r2[i]=-1.0q/4.0q*(-(sv[i]+vtmp1[i])*powq(Rar2[i],3.0q) -3.0q*(yv2[i]+g2+2.0q*cg*yv[i]+a2-vtmp2[i])*B0m1_r2[i]);
    C0m3_r2[i]=1.0q/4.0q*(powq(Rar2[i],3.0q)*(yv[i]+sv[i]*(cf+d)+cg)+3.0q*((onepf2p2ef-powq((cf+d),2.0q))*sv2[i]+2.0q*sv[i]*(fg+eg-cg*(cf+d))+(1.0q-c2)*g2+a2)*C0m1_r2[i]);

    i=is1[j];
    A0m3_s1[i]=1.0q/4.0q*(Rdotp[i]*powq(Ras1[i],3.0q) +3.0q*(powq(Ras1[i],2.0q)-(powq(Rdotp[i],2.0q)))*A0m1_s1[i]);
    C0m3_s1[i]=1.0q/4.0q*(sRdott[i]*powq(Ras1[i],3.0q) +3.0q*(powq(Ras1[i],2.0q)-powq(sRdott[i],2.0q))*C0m1_s1[i]);

    i=is2[j];
    vtmp1[i]=((dh+c)*yv[i]+hm+em)/(oneph2p2eh);
    vtmp2[i]= powq((dh+c)*yv[i]+hm+em,2.0q)/(oneph2p2eh);
    A0m3_s2[i]=-1.0q/4.0q*(-(rv[i]+vtmp1[i])*powq(Ras2[i],3.0q) -3.0q*(yv2[i]+m2+2.0q*dm*yv[i]+a2-vtmp2[i])*A0m1_s2[i]);
    C0m3_s2[i]=1.0q/4.0q*(powq(Ras2[i],3.0q)*(yv[i]+rv[i]*(dh+c)+dm)+3.0q*((oneph2p2eh-powq((dh+c),2.0q))*rv2[i]+2.0q*rv[i]*(hm+em-dm*(dh+c))+(1.0q-d2)*m2+a2)*C0m1_s2[i]);
  }

  // ABC1m3
  for(j=0;j<4;j++) {
    i=ir1[j];
    B1m3_r1[i]= 1.0q/5.0q*powq(Rar1[i],5.0q)-(d*yv[i]+e*r1)*B0m3_r1[i];

    i=ir2[j];
    B1m3_r2[i]= 1.0q/5.0q/(onepf2p2ef)*powq(Rar2[i],5.0q)-1.0q/(onepf2p2ef)*(yv[i]*(cf+d)+fg+eg)*B0m3_r2[i];

    i=is1[j];
    A1m3_s1[i]= 1.0q/5.0q*powq(Ras1[i],5.0q)-(c*yv[i]+e*s1)*A0m3_s1[i];

    i=is2[j];
    A1m3_s2[i]= -1.0q/5.0q/(oneph2p2eh)*(-powq(Ras2[i],5.0q)+5.0q*(yv[i]*(dh+c)+hm+em)*A0m3_s2[i]);
  }

  // ABC2m1
  for(j=0;j<4;j++) {
    i=ir1[j];
    B2m1_r1[i]= -1.0q/3.0q*(-sv[i]*powq(Rar1[i],3.0q)+B0m3_r1[i]-(d*yv[i]+e*r1)*(-3.0q)*B1m1_r1[i]);

    i=ir2[j];
    B2m1_r2[i]= 1.0q/(-3.0q)/(onepf2p2ef)*(-sv[i]*powq(Rar2[i],3.0q)+B0m3_r2[i]-(yv[i]*(cf+d)+fg+eg)*(-3.0q)*B1m1_r2[i]);

    i=is1[j];
    A2m1_s1[i]= -1.0q/3.0q*(-rv[i]*powq(Ras1[i],3.0q)+A0m3_s1[i]-(c*yv[i]+e*s1)*(-3.0q)*A1m1_s1[i]);

    i=is2[j];
    A2m1_s2[i]= 1.0q/(-3.0q)/(oneph2p2eh)*(-rv[i]*powq(Ras2[i],3.0q)+A0m3_s2[i]-(yv[i]*(dh+c)+hm+em)*(-3.0q)*A1m1_s2[i]);
  }

  // ABC21
  for(j=0;j<4;j++) {
    i=ir1[j];
    B21_r1[i]= -1.0q*(-sv[i]*powq(Rar1[i],1.0q)+B0m1_r1[i]-(d*yv[i]+e*r1)*(-1.0q)*B11_r1[i]);

    i=ir2[j];
    B21_r2[i]= -1.0q/(onepf2p2ef)*(-sv[i]*powq(Rar2[i],1.0q)+B0m1_r2[i]-(yv[i]*(cf+d)+fg+eg)*(-1.0q)*B11_r2[i]);

    i=is1[j];
    A21_s1[i]= -1.0q*(-rv[i]*powq(Ras1[i],1.0q)+A0m1_s1[i]-(c*yv[i]+e*s1)*(-1.0q)*A11_s1[i]);

    i=is2[j];
    A21_s2[i]= -1.0q/(oneph2p2eh)*(-rv[i]*powq(Ras2[i],1.0q)+A0m1_s2[i]-(yv[i]*(dh+c)+hm+em)*(-1.0q)*A11_s2[i]);
  }

  // ABC31
  for(j=0;j<4;j++) {
    i=ir1[j];
    B31_r1[i]= -1.0q*(-sv[i]*sv[i]*powq(Rar1[i],1.0q)+2.0q*B1m1_r1[i]-(d*yv[i]+e*r1)*(-1.0q)*B21_r1[i]);

    i=ir2[j];
    B31_r2[i]= -1.0q/(onepf2p2ef)*(-sv[i]*sv[i]*powq(Rar2[i],1.0q)+2.0q*B1m1_r2[i]-(yv[i]*(cf+d)+fg+eg)*(-1.0q)*B21_r2[i]);

    i=is1[j];
    A31_s1[i]= -1.0q*(-rv[i]*rv[i]*powq(Ras1[i],1.0q)+2.0q*A1m1_s1[i]-(c*yv[i]+e*s1)*(-1.0q)*A21_s1[i]);

    i=is2[j];
    A31_s2[i]= -1.0q/(oneph2p2eh)*(-rv[i]*rv[i]*powq(Ras2[i],1.0q)+2.0q*A1m1_s2[i]-(yv[i]*(dh+c)+hm+em)*(-1.0q)*A21_s2[i]);
  }

  //
  // DOUBLE INTEGRALS
  //
  /* Double integrals are built from previous linear integrals and a set
     of 3 seed functions.
     Eikl = /int/int r^i y^k Ra^-l dr dy
     Fjkl = /int/int s^j y^k Ra^-l ds dy
     Dijl = /int/int r^i s^j Ra^-l dr ds
     For Eikl and Fjkl integrals, integration paths are different for r1,r2 or
     s1,s2 as consequence of s2 = h*r+m and r2= f*s+g.
     We introduce also the following linear integrals:
     AASijl = /int r^i s(r)^j Ra^-l dr = /int r^i (h*r+m)^j Ra^-l dr - /int r^i s1^j Ra^-l dr
     BBrijl = /int r(s)^i s^j Ra^-l ds = /int (f*s+g)^i s^j Ra^-l ds - /int r1^i s^j Ra^-l ds
  */

  //
  // SEED FUNCTIONS: DEF003
  //
  for(j=0;j<4;j++) {
    i=ir1[j];
    root= 1.0q/sqrtq((1.0q-d2)*a2+((1.0q-d2)*(1.0q-e2)-(c-de)*(c-de))*rv2[i]);
    F003_r1[i]= 2.0q*root*atanq(((1.0q-d)*(rRa[i]-sv[i]-d*yv[i]-e*rv[i])+(1.0q-d2)*yv[i]+(c-de)*rv[i])*root);
    root= 1.0q/sqrtq((1.0q-e2)*a2+(1.0q-c2-d2-e2+2.0q*c*de)*yv2[i]);
    D003_r1[i]=2.0q*root*atanq(((1.0q-e)*(Rar1[i]+rv[i]-sv[i])+(c-d)*yv[i])*root);

    i=ir2[j];
    tp1=sqrtq(onepf2p2ef);
    tp2=(yv[i]*(cf+d)+fg+eg)/tp1;
    tp3=sqrtq(tp1-cf-d);
    tp4=(tp1*(yv[i]+cg)-tp2*(cf+d))/tp3;
    tp5=sqrtq(-tp4*tp4-(tp1+cf+d)*tp2*tp2+(tp1+cf+d)*(a2+g2+yv2[i]+2.0q*cg*yv[i]));
    F003_r2[i]=2.0q/tp3/tp5*atanq((tp3*(Rar2[i]-tp1*sv[i]-tp2)+tp4)/tp5);

    tp1=sqrtq(onepf2p2ef);
    tp2=sqrtq(tp1-f-e);
    root=1.0q/sqrtq(yv2[i]*(1.0q-c2-d2-e2+2.0q*c*de)+a2*(1.0q-e2));
    D003_r2[i]=2.0q*root*atanq((tp2*tp2*Rar2[i]+(tp1*c-cf-d)*yv[i]+tp1*g-fg-eg-tp1*tp2*tp2*sv[i])*root);

    i=is1[j];
    root= 1.0q/sqrtq((1.0q-c2)*a2+((1.0q-c2)*(1.0q-e2)-(d-ce)*(d-ce))*sv2[i]);
    E003_s1[i]=2.0q*root*atanq(((1.0q-c)*(sRa[i]-rv[i]-c*yv[i]-e*sv[i])+(1.0q-c2)*yv[i]+(d-c*e)*sv[i])*root);
    root= 1.0q/sqrtq((1.0q-e2)*a2+(1.0q-c2-d2-e2+2.0q*c*de)*yv2[i]);
    D003_s1[i]=2.0q*root*atanq(((1.0q-e)*(Ras1[i]+sv[i]-rv[i])+(d-c)*yv[i])*root);

    i=is2[j];
    tp1=sqrtq(oneph2p2eh);
    tp2=(yv[i]*(dh+c)+hm+em)/tp1;
    tp3=sqrtq(tp1-dh-c);
    tp4=(tp1*(yv[i]+dm)-tp2*(dh+c))/tp3;
    tp5=sqrtq(-tp4*tp4-(tp1+dh+c)*tp2*tp2+(tp1+dh+c)*(a2+m2+yv2[i]+2.0q*dm*yv[i]));
    E003_s2[i]=2.0q/tp3/tp5*atanq((tp3*(Ras2[i]-tp1*rv[i]-tp2)+tp4)/tp5);

    tp1=sqrtq(oneph2p2eh);
    tp2=sqrtq(tp1-h-e);
    root=1.0q/sqrtq(yv2[i]*(1.0q-c2-d2-e2+2.0q*c*de)+a2*(1.0q-e2));
    D003_s2[i]=2.0q*root*atanq((tp2*tp2*Ras2[i]+(tp1*d-dh-c)*yv[i]+tp1*m-hm-em-tp1*tp2*tp2*rv[i])*root);
  }

  for(i=0;i<8;i++) {
    D003[i]= 0.5q*D003_r1[i] +0.5q*D003_r2[i] +0.5q*D003_s1[i] +0.5q*D003_s2[i];
    AAS001[i]= A01_s1[i]+A01_s2[i];
    BBR001[i]= B01_r1[i]+B01_r2[i];
    D103[i]=1.0q/(1.0q-e2)*((+e*AAS001[i] -BBR001[i]) -(c-de)*yv[i]*D003[i]);
    D013[i]=1.0q/(1.0q-e2)*((+e*BBR001[i] -AAS001[i]) -(d-ce)*yv[i]*D003[i]);
  }

  for(j=0;j<4;j++) {
    i=ir1[j];
    F103_r1[i]=1.0q/(1.0q-d2)*((+d*B01_r1[i] -C01_r1[i]) +(cd-e)*rv[i]*F003_r1[i]);
    F013_r1[i]=1.0q/(1.0q-d2)*((+d*powq(sfr[i],0.0q)*C01_r1[i] -powq(yv[i],0.0q)*B01_r1[i]) +(de-c)*rv[i]*F003_r1[i]);
    i=ir2[j];
    F103_r2[i]=1.0q/(onepf2p2ef -(cf+d)*(cf+d))*((cf+d)*powq(yv[i],0.0q)*B01_r2[i] -C01_r2[i] -(fg+eg-cg*(cf+d))*F003_r2[i]);
    F013_r2[i]=1.0q/(onepf2p2ef -(cf+d)*(cf+d))*((cf+d)*powq(sv[i],0.0q)*C01_r2[i] -(onepf2p2ef)*powq(yv[i],0.0q)*B01_r2[i] +((cf+d)*(fg+eg)-(onepf2p2ef)*cg)*F003_r2[i]);
    i=is1[j];
    E103_s1[i]=1.0q/(1.0q-c2)*((+c*A01_s1[i] -C01_s1[i]) +(cd-e)*sv[i]*E003_s1[i]);
    E013_s1[i]=1.0q/(1.0q-c2)*((+c*powq(rfs[i],0.0q)*C01_s1[i] -powq(yv[i],0.0q)*A01_s1[i]) +(ce-d)*sv[i]*E003_s1[i]);
    i=is2[j];
    E103_s2[i]=1.0q/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*powq(yv[i],0.0q)*A01_s2[i] -C01_s2[i] -(hm+em-dm*(dh+c))*E003_s2[i]);
    E013_s2[i]=1.0q/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*powq(rv[i],0.0q)*C01_s2[i] -(oneph2p2eh)*powq(yv[i],0.0q)*A01_s2[i] +((dh+c)*(hm+em) -(oneph2p2eh)*dm)*E003_s2[i]);
  }

  for(j=0;j<4;j++) {
    i=ir1[j];
    F001_r1[i] = 1.0q/(1.0q-d2)/(1.0q-2.0q)*(((a2+rv2[i])*(1.0q-d2)+2.0q*cd*e*rv2[i]-1.0q*rv2[i]*(e2+c2))*F003_r1[i]+(-powq(sfr[i],0.0q)*(sfr[i]*(1.0q-d2)-c*rv[i]*d+e*rv[i])*C01_r1[i]-powq(yv[i],0.0q)*(yv[i]*(1.0q-d2)+c*rv[i]-e*d*rv[i])*B01_r1[i]));
    F101_r1[i]=1.0q/(1.0q-d2)*(1.0q/(1.0q-2.0q)*(+d*powq(yv[i],0.0q)*B0m1_r1[i]-powq(sfr[i],0.0q)*C0m1_r1[i])+(cd-e)*rv[i]*F001_r1[i]);
    F011_r1[i]=1.0q/(1.0q-d2)*(1.0q/(1.0q-2.0q)*(+d*powq(sfr[i],0.0q)*C0m1_r1[i]-powq(yv[i],0.0q)*B0m1_r1[i])+(de-c)*rv[i]*F001_r1[i]);
    i=ir2[j];
    tp1=onepf2p2ef-(cf+d)*(cf+d);
    F001_r2[i] = 1.0q/(1.0q-2.0q)*(((fg+eg)*(cf+d)-cg*(onepf2p2ef))/tp1*B01_r2[i]-powq(yv[i],(1.0q))*B01_r2[i] +(cg*(cf+d)-(fg+eg))/tp1*C01_r2[i]-powq(sv[i],(1.0q))*C01_r2[i] +(a2+g2)*F003_r2[i] - ((fg+eg)*(fg+eg-cg*(cf+d))+cg*(cg*(onepf2p2ef)-(cf+d)*(fg+eg)))/tp1*F003_r2[i]);
    F101_r2[i]=1.0q/(1.0q-2.0q)/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*powq(yv[i],0.0q)*B0m1_r2[i]-powq(sv[i],0.0q)*C0m1_r2[i]-(1.0q-2.0q)*(fg+eg-cg*(cf+d))*F001_r2[i]);
    F011_r2[i]=1.0q/(1.0q-2.0q)/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*powq(sv[i],0.0q)*C0m1_r2[i]-(onepf2p2ef)*powq(yv[i],0.0q)*B0m1_r2[i]+(1.0q-2.0q)*((cf+d)*(fg+eg)-(onepf2p2ef)*cg)*F001_r2[i]);
    i=is1[j];
    E001_s1[i] = 1.0q/(1.0q-c2)/(1.0q-2.0q)*(((a2+sv2[i])*(1.0q-c2)+2.0q*c*d*e*sv2[i]-sv2[i]*(e2+d2))*E003_s1[i]+(-powq(rfs[i],0.0q)*(rfs[i]*(1.0q-c2)-d*sv[i]*c+e*sv[i])*C01_s1[i]-powq(yv[i],0.0q)*(yv[i]*(1.0q-c2)+d*sv[i]-e*c*sv[i])*A01_s1[i]));

    E101_s1[i]=1.0q/(1.0q-c2)*(1.0q/(1.0q-2.0q)*(+c*powq(yv[i],0.0q)*A0m1_s1[i]-powq(rfs[i],0.0q)*C0m1_s1[i])+(cd-e)*sv[i]*E001_s1[i]);
    E011_s1[i]=1.0q/(1.0q-c2)*(1.0q/(1.0q-2.0q)*(+c*powq(rfs[i],0.0q)*C0m1_s1[i]-powq(yv[i],0.0q)*A0m1_s1[i])+(c*e-d)*sv[i]*E001_s1[i]);
    i=is2[j];
    tp1=oneph2p2eh-(dh+c)*(dh+c);
    E001_s2[i] = 1.0q/(1.0q-2.0q)*(((hm+em)*(dh+c)-dm*(oneph2p2eh))/tp1*A01_s2[i]-powq(yv[i],(1.0q))*A01_s2[i] +(dm*(dh+c)-(hm+em))/tp1*C01_s2[i]-powq(rv[i],(1.0q))*C01_s2[i] +(a2+m2)*E003_s2[i] - ((hm+em)*(hm+em-dm*(dh+c))+dm*(dm*(oneph2p2eh)-(dh+c)*(hm+em)))/tp1*E003_s2[i]);
    E101_s2[i]=1.0q/(1.0q-2.0q)/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*powq(yv[i],0.0q)*A0m1_s2[i]-powq(rv[i],0.0q)*C0m1_s2[i]-(1.0q-2.0q)*(hm+em-dm*(dh+c))*E001_s2[i]);
    E011_s2[i]=1.0q/(1.0q-2.0q)/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*powq(rv[i],0.0q)*C0m1_s2[i]-(oneph2p2eh)*powq(yv[i],0.0q)*A0m1_s2[i]+(1.0q-2.0q)*((dh+c)*(hm+em)-(oneph2p2eh)*dm)*E001_s2[i]);
  }

  for(w=0;w<8;w++) {
    D001[w] =1.0q/((1.0q-e2)*(-1.0q))*( ((a2+yv2[w]) * (1.0q-e2) +2.0q*cd*e*1.0q * yv2[w] - 1.0q*yv2[w]*(c2+d2)) * D003[w]-(yv[w]*(d-ce) +powq(m,(1.0q))*(1.0q-e2)) * A01_s2[w]-(h*(1.0q-e2))*A11_s2[w]-(yv[w]*(d-ce)+powq(sv[w],(1.0q))*(1.0q-e2)) * A01_s1[w]-(yv[w]*(c-de)+(1.0q-e2)*g)*B01_r2[w]-(f*(1.0q-e2))*B11_r2[w]-(yv[w]*(c-e*d)+powq(rv[w],(1.0q))*(1.0q-e2)) * B01_r1[w]);
    AAS00m1[w]= A0m1_s1[w]+A0m1_s2[w];
    BBR00m1[w]= B0m1_r1[w]+B0m1_r2[w];
    D101[w]=1.0q/(1.0q-e2)*(1.0q/(1.0q-2.0q)*(+e*AAS00m1[w]-BBR00m1[w])-(c-de)*yv[w]*D001[w]);
    D011[w]=1.0q/(1.0q-e2)*(1.0q/(1.0q-2.0q)*(+e*BBR00m1[w]-AAS00m1[w])-(d-ce)*yv[w]*D001[w]);
  }

  for(j=0;j<4;j++) {
    w=ir1[j];
    F00m1_r1[w] = 1.0q/(1.0q-d2)/(-3.0q)*((-(a2+rv2[w])*(1.0q-d2)+2.0q*cd*e*-rv2[w]+rv2[w]*(e2+c2))*F001_r1[w]+(-powq(sfr[w],0.0q)*(sfr[w]*(1.0q-d2)-c*rv[w]*d+e*rv[w])*C0m1_r1[w]-powq(yv[w],0.0q)*(yv[w]*(1.0q-d2)+c*rv[w]-e*d*rv[w])*B0m1_r1[w]));
    w=ir2[j];
    tp1=onepf2p2ef-(cf+d)*(cf+d);
    F00m1_r2[w] = 1.0q/(-3.0q)*(((fg+eg)*(cf+d)-cg*(onepf2p2ef))/tp1*B0m1_r2[w]-powq(yv[w],(1.0q))*B0m1_r2[w] +(cg*(cf+d)-(fg+eg))/tp1*C0m1_r2[w]-powq(sv[w],(1.0q))*C0m1_r2[w] -1.0q*(a2+g2)*F001_r2[w] +1.0q*((fg+eg)*(fg+eg-cg*(cf+d))+cg*(cg*(onepf2p2ef)-(cf+d)*(fg+eg)))/tp1*F001_r2[w]);
    w=is1[j];
    E00m1_s1[w] = 1.0q/(1.0q-c2)/(-3.0q)*((-(a2+sv2[w])*(1.0q-c2)+2.0q*cd*e*-sv2[w]+sv2[w]*(e2+d2))*E001_s1[w]+(-powq(rfs[w],0.0q)*(rfs[w]*(1.0q-c2)-d*sv[w]*c+e*sv[w])*C0m1_s1[w]-powq(yv[w],0.0q)*(yv[w]*(1.0q-c2)+d*sv[w]-e*c*sv[w])*A0m1_s1[w]));
    w=is2[j];
    tp1=oneph2p2eh-(dh+c)*(dh+c);
    E00m1_s2[w] = 1.0q/(-3.0q)*(((hm+em)*(dh+c)-dm*(oneph2p2eh))/tp1*A0m1_s2[w]-powq(yv[w],(1.0q))*A0m1_s2[w] +(dm*(dh+c)-(hm+em))/tp1*C0m1_s2[w]-powq(rv[w],(1.0q))*C0m1_s2[w] -1.0q*(a2+m2)*E001_s2[w] +1.0q*((hm+em)*(hm+em-dm*(dh+c))+dm*(dm*(oneph2p2eh)-(dh+c)*(hm+em)))/tp1*E001_s2[w]);
  }

  for(w=0;w<8;w++) {
    D00m1[w] =1.0q/((1.0q-e2)*(-3.0q))*( (-(a2+yv2[w]) * (1.0q-e2) +2.0q*cd*e*-1.0q * yv2[w] +1.0q*yv2[w]*(c2+d2)) * D001[w]-(yv[w]*(d-e*c) +powq(m,(1.0q))*(1.0q-e2)) * A0m1_s2[w]-(h*(1.0q-e2))*A1m1_s2[w]-(yv[w]*(d-e*c)+powq(sv[w],(1.0q))*(1.0q-e2)) * A0m1_s1[w]-(yv[w]*(c-e*d)+powq(g,(1.0q))*(1.0q-e2)) * B0m1_r2[w]-(f*(1.0q-e2))*B1m1_r2[w]-(yv[w]*(c-e*d)+powq(rv[w],(1.0q))*(1.0q-e2)) * B0m1_r1[w]);
    D00m3[w] =1.0q/((1.0q-e2)*(-5.0q))*( (-3.0q*(a2+yv2[w]) * (1.0q-e2) +2.0q*cd*e*-3.0q * yv2[w] +3.0q*yv2[w]*(c2+d2)) * D00m1[w]-(yv[w]*(d-e*c) +powq(m,(1.0q))*(1.0q-e2)) * A0m3_s2[w]-(h*(1.0q-e2))*A1m3_s2[w]-(yv[w]*(d-e*c)+powq(sv[w],(1.0q))*(1.0q-e2)) * A0m3_s1[w]-(yv[w]*(c-e*d)+powq(g,(1.0q))*(1.0q-e2)) * B0m3_r2[w]-(f*(1.0q-e2))*B1m3_r2[w]-(yv[w]*(c-e*d)+powq(rv[w],(1.0q))*(1.0q-e2)) * B0m3_r1[w]);
  }

  for(j=0;j<4;j++) {
    w=ir1[j];
    F00m3_r1[w] = 1.0q/(1.0q-d2)/(-5.0q)*((-3.0q*(a2+rv2[w])*(1.0q-d2)+2.0q*cd*e*-3.0q*rv2[w]+3.0q*rv2[w]*(e2+c2))*F00m1_r1[w]+(-powq(sfr[w],0.0q)*(sfr[w]*(1.0q-d2)-cd*rv[w]+e*rv[w])*C0m3_r1[w]-powq(yv[w],0.0q)*(yv[w]*(1.0q-d2)+c*rv[w]-e*d*rv[w])*B0m3_r1[w]));
    w=ir2[j];
    tp1=onepf2p2ef-(cf+d)*(cf+d);
    F00m3_r2[w] = 1.0q/(-5.0q)*(((fg+eg)*(cf+d)-cg*(onepf2p2ef))/tp1*B0m3_r2[w]-powq(yv[w],(1.0q))*B0m3_r2[w] +(cg*(cf+d)-(fg+eg))/tp1*C0m3_r2[w]-powq(sv[w],(1.0q))*C0m3_r2[w] -3.0q*(a2+g2)*F00m1_r2[w] +3.0q*((fg+eg)*(fg+eg-cg*(cf+d))+cg*(cg*(onepf2p2ef)-(cf+d)*(fg+eg)))/tp1*F00m1_r2[w]);
    w=is1[j];
    E00m3_s1[w] = 1.0q/(1.0q-c2)/(-5.0q)*((-3.0q*(a2+sv2[w])*(1.0q-c2)+2.0q*cd*e*-3.0q*sv2[w]+3.0q*sv2[w]*(e2+d2))*E00m1_s1[w]+(-powq(rfs[w],0.0q)*(rfs[w]*(1.0q-c2)-cd*sv[w]+e*sv[w])*C0m3_s1[w]-powq(yv[w],0.0q)*(yv[w]*(1.0q-c2)+d*sv[w]-e*c*sv[w])*A0m3_s1[w]));
    w=is2[j];
    tp1=oneph2p2eh-(dh+c)*(dh+c);
    E00m3_s2[w] = 1.0q/(-5.0q)*(((hm+em)*(dh+c)-dm*(oneph2p2eh))/tp1*A0m3_s2[w]-powq(yv[w],(1.0q))*A0m3_s2[w] +(dm*(dh+c)-(hm+em))/tp1*C0m3_s2[w]-powq(rv[w],(1.0q))*C0m3_s2[w] -3.0q*(a2+m2)*E00m1_s2[w] +3.0q*((hm+em)*(hm+em-dm*(dh+c))+dm*(dm*(oneph2p2eh)-(dh+c)*(hm+em)))/tp1*E00m1_s2[w]);
  }

  for(j=0;j<4;j++) {
    w=ir1[j];
    F10m1_r1[w]=1.0q/(1.0q-d2)*(1.0q/(-3.0q)*(+d*powq(yv[w],0.0q)*B0m3_r1[w]-powq(sfr[w],0.0q)*C0m3_r1[w])+(cd-e)*rv[w]*F00m1_r1[w]);
    w=ir2[j];
    F10m1_r2[w]=1.0q/(-3.0q)/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*powq(yv[w],0.0q)*B0m3_r2[w]-powq(sv[w],0.0q)*C0m3_r2[w] -(-3.0q)*(fg+eg-cg*(cf+d))*F00m1_r2[w]);
    w=is1[j];
    E10m1_s1[w]=1.0q/(1.0q-c2)*(1.0q/(-3.0q)*(+c*powq(yv[w],0.0q)*A0m3_s1[w]-powq(rfs[w],0.0q)*C0m3_s1[w])+(cd-e)*sv[w]*E00m1_s1[w]);
    w=is2[j];
    E10m1_s2[w]=1.0q/(-3.0q)/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*powq(yv[w],0.0q)*A0m3_s2[w]-powq(rv[w],0.0q)*C0m3_s2[w] -(-3.0q)*(hm+em-dm*(dh+c))*E00m1_s2[w]);
  }

  for(w=0;w<8;w++) {
    AAS00m3[w]= A0m3_s1[w]+A0m3_s2[w];
    BBR00m3[w]= B0m3_r1[w]+B0m3_r2[w];
    D10m1[w]=1.0q/(1.0q-e2)*(1.0q/(-3.0q)*(+e*AAS00m3[w]-BBR00m3[w])-(c-de)*yv[w]*D00m1[w]);
    D01m1[w]=1.0q/(1.0q-e2)*(1.0q/(-3.0q)*(+e*BBR00m3[w]-AAS00m3[w])-(d-ce)*yv[w]*D00m1[w]);
  }

  for(j=0;j<4;j++) {
    w=ir1[j];
    F201_r1[w]=1.0q/(1.0q-d2)*(1.0q/(-1.0q)*(F00m1_r1[w]+d*powq(yv[w],0.0q)*B1m1_r1[w]-powq(sfr[w],1.0q)*C0m1_r1[w])+(cd-e)*rv[w]*F101_r1[w]);
    F111_r1[w]=1.0q/(1.0q-d2)*(1.0q/(-1.0q)*(-d*F00m1_r1[w]+d*powq(sfr[w],1.0q)*C0m1_r1[w]-powq(yv[w],0.0q)*B1m1_r1[w])+(de-c)*rv[w]*F101_r1[w]);
    w=ir2[j];
    F201_r2[w]=1.0q/(-1.0q)/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*powq(yv[w],0.0q)*B1m1_r2[w]-powq(sv[w],1.0q)*C0m1_r2[w]+1.0q*F00m1_r2[w]-(-1.0q)*(fg+eg-cg*(cf+d))*F101_r2[w]);
    F111_r2[w]=1.0q/(-1.0q)/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*powq(sv[w],1.0q)*C0m1_r2[w]-(onepf2p2ef)*powq(yv[w],0.0q)*B1m1_r2[w]-(cf+d)*F00m1_r2[w]+(-1.0q)*((cf+d)*(fg+eg)-(onepf2p2ef)*cg)*F101_r2[w]);
    w=is1[j];
    E201_s1[w]=1.0q/(1.0q-c2)*(1.0q/(-1.0q)*(E00m1_s1[w]+c*powq(yv[w],0.0q)*A1m1_s1[w]-powq(rfs[w],1.0q)*C0m1_s1[w])+(cd-e)*sv[w]*E101_s1[w]);
    E111_s1[w]=1.0q/(1.0q-c2)*(1.0q/(1.0q-2.0q)*(-c*E00m1_s1[w]+c*powq(rfs[w],1.0q)*C0m1_s1[w]-powq(yv[w],0.0q)*A1m1_s1[w])+(ce-d)*sv[w]*E101_s1[w]);
    w=is2[j];
    E201_s2[w]=1.0q/(-1.0q)/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*powq(yv[w],0.0q)*A1m1_s2[w]-powq(rv[w],1.0q)*C0m1_s2[w]+E00m1_s2[w]-(-1.0q)*(hm+em-dm*(dh+c))*E101_s2[w]);
    E111_s2[w]=1.0q/(1.0q-2.0q)/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*powq(rv[w],1.0q)*C0m1_s2[w]-(oneph2p2eh)*powq(yv[w],0.0q)*A1m1_s2[w]-(dh+c)*E00m1_s2[w]+(-1.0q)*((dh+c)*(hm+em)-(oneph2p2eh)*dm)*E101_s2[w]);
  }

  for(w=0;w<8;w++) {
    AAS10m1[w]= A1m1_s1[w]+A1m1_s2[w];
    AAS01m1[w]= h*A1m1_s2[w] +m* A0m1_s2[w] +s1*A0m1_s1[w];
    AAS11m1[w]= h*A2m1_s2[w] +m* A1m1_s2[w] +s1*A1m1_s1[w];
    BBR01m1[w]= +B1m1_r2[w] +B1m1_r1[w];
    BBR10m1[w]= f*B1m1_r2[w] +g* B0m1_r2[w] +r1*B0m1_r1[w];
    BBR11m1[w]= f*B2m1_r2[w] +g* B1m1_r2[w] +r1*B1m1_r1[w];
    D201[w]= 1.0q/(1.0q-e2)*(1.0q/(-1.0q)*(D00m1[w] +e*AAS10m1[w] -BBR10m1[w]) -(c-de)*yv[w]*D101[w]);

    //
    D111_1[w]= 1.0q/(1.0q-e2)*(1.0q/(-1.0q)*(-e*D00m1[w] +e*AAS01m1[w] -BBR01m1[w]) -(c-de)*yv[w]*D011[w]);
    D111_2[w]= 1.0q/(1.0q-e2)*(1.0q/(-1.0q)*(-e*D00m1[w] +e*BBR10m1[w] -AAS10m1[w]) -(d-ce)*yv[w]*D101[w]);
    //

    D021[w]=1.0q/(1.0q-e2)*(1.0q/(-1.0q)*(D00m1[w]+e*BBR01m1[w]-AAS01m1[w])-(d-ce)*yv[w]*D011[w]);
  }

 for(w=0;w<8;w++) {
       D111[w]= 0.5q*D111_1[w] +0.5q*D111_2[w];
 }


  for(j=0;j<4;j++) {
    w=ir1[j];
    F203_r1[w]=1.0q/(1.0q-d2)*((F001_r1[w]+d*powq(yv[w],0.0q)*B11_r1[w]-powq(sfr[w],1.0q)*C01_r1[w])+(cd-e)*rv[w]*F103_r1[w]);
    w=ir2[j];
    F203_r2[w]=1.0q/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*powq(yv[w],0.0q)*B11_r2[w]-powq(sv[w],1.0q)*C01_r2[w]+F001_r2[w]-(fg+eg-cg*(cf+d))*F103_r2[w]);
    w=is1[j];
    E203_s1[w]=1.0q/(1.0q-c2)*((E001_s1[w]+c*powq(yv[w],0.0q)*A11_s1[w]-powq(rfs[w],1.0q)*C01_s1[w])+(cd-e)*sv[w]*E103_s1[w]);
    w=is2[j];
    E203_s2[w]=1.0q/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*powq(yv[w],0.0q)*A11_s2[w]-powq(rv[w],1.0q)*C01_s2[w]+E001_s2[w]-(hm+em-dm*(dh+c))*E103_s2[w]);
  }

  for(w=0;w<8;w++) {
    AAS101[w]= A11_s1[w]+A11_s2[w];
    AAS011[w]= h*A11_s2[w] +m* A01_s2[w] +s1*A01_s1[w];
    BBR101[w]= f*B11_r2[w] +g* B01_r2[w] +r1*B01_r1[w];
    BBR011[w]= +B11_r2[w] +B11_r1[w];
    D203[w]=1.0q/(1.0q-e2)*((D001[w]+e*AAS101[w]-BBR101[w])-(c-de)*yv[w]*D103[w]);
    D023[w]=1.0q/(1.0q-e2)*((D001[w]+e*BBR011[w]-AAS011[w])-(d-ce)*yv[w]*D013[w]);
  }

  for(j=0;j<4;j++) {
    w=ir1[j];
    F023_r1[w]=1.0q/(1.0q-d2)*((F001_r1[w]+d*powq(sfr[w],0.0q)*C11_r1[w]-powq(yv[w],1.0q)*B01_r1[w])+(de-c)*rv[w]*F013_r1[w]);
    w=ir2[j];
    F023_r2[w]=1.0q/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*powq(sv[w],0.0q)*C11_r2[w]-(onepf2p2ef)*powq(yv[w],1.0q)*B01_r2[w]+(onepf2p2ef)*F001_r2[w]+((cf+d)*(fg+eg)-(onepf2p2ef)*cg)*F013_r2[w]);
    w=is1[j];
    E023_s1[w]=1.0q/(1.0q-c2)*((E001_s1[w]+c*powq(rfs[w],0.0q)*C11_s1[w]-powq(yv[w],1.0q)*A01_s1[w])+(ce-d)*sv[w]*E013_s1[w]);
    w=is2[j];
    E023_s2[w]=1.0q/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*powq(rv[w],0.0q)*C11_s2[w]-(oneph2p2eh)*powq(yv[w],1.0q)*A01_s2[w]+(oneph2p2eh)*E001_s2[w]+((dh+c)*(hm+em)-(oneph2p2eh)*dm)*E013_s2[w]);
  }

  for(j=0;j<4;j++) {
    w=ir1[j];
    F113_r1[w]=1.0q/(1.0q-d2)*((-d*F001_r1[w]+d*powq(yv[w],1.0q)*B01_r1[w]-powq(sfr[w],0.0q)*C11_r1[w])+(cd-e)*rv[w]*F013_r1[w]);
    w=ir2[j];
    F113_r2[w]=1.0q/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*powq(yv[w],1.0q)*B01_r2[w]-powq(sv[w],0.0q)*C11_r2[w]-(cf+d)*F001_r2[w]-(fg+eg-cg*(cf+d))*F013_r2[w]);
    w=is1[j];
    E113_s1[w]=1.0q/(1.0q-c2)*((-c*E001_s1[w]+c*powq(yv[w],1.0q)*A01_s1[w]-powq(rfs[w],0.0q)*C11_s1[w])+(cd-e)*sv[w]*E013_s1[w]);
    w=is2[j];
    E113_s2[w]=1.0q/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*powq(yv[w],1.0q)*A01_s2[w]-powq(rv[w],0.0q)*C11_s2[w]-(dh+c)*E001_s2[w]-(hm+em-dm*(dh+c))*E013_s2[w]);
  }

  for(w=0;w<8;w++) {
    D113_1[w]=1.0q/(1.0q-e2)*((-e*D001[w]+e*AAS011[w]-BBR011[w])-(c-de)*yv[w]*D013[w]);
    D113_2[w]=1.0q/(1.0q-e2)*((-e*D001[w]+e*BBR101[w]-AAS101[w])-(d-ce)*yv[w]*D103[w]);
  }
  for(w=0;w<8;w++) {
    D113[w]= 0.5q*D113_1[w] +0.5q*D113_2[w];
  }

  for(j=0;j<4;j++) {
    w=ir1[j];
    F213_r1[w]=1.0q/(1.0q-d2)*((F011_r1[w]-d*F101_r1[w]+d*powq(yv[w],1.0q)*B11_r1[w]-powq(sfr[w],1.0q)*C11_r1[w])+(cd-e)*rv[w]*F113_r1[w]);
    F123_r1[w]=1.0q/(1.0q-d2)*((F101_r1[w]-d*F011_r1[w]+d*powq(sfr[w],1.0q)*C11_r1[w]-powq(yv[w],1.0q)*B11_r1[w])+(de-c)*rv[w]*F113_r1[w]);
    w=ir2[j];
    F213_r2[w]=1.0q/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*powq(yv[w],1.0q)*B11_r2[w]-powq(sv[w],1.0q)*C11_r2[w]+F011_r2[w]-(cf+d)*F101_r2[w]-(fg+eg-cg*(cf+d))*F113_r2[w]);
    F123_r2[w]=1.0q/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*powq(sv[w],1.0q)*C11_r2[w]-(onepf2p2ef)*powq(yv[w],1.0q)*B11_r2[w]+(onepf2p2ef)*F101_r2[w]-(cf+d)*F011_r2[w]+((cf+d)*(fg+eg)-(onepf2p2ef)*cg)*F113_r2[w]);
    w=is1[j];
    E213_s1[w]=1.0q/(1.0q-c2)*((E011_s1[w]-c*E101_s1[w]+c*powq(yv[w],1.0q)*A11_s1[w]-powq(rfs[w],1.0q)*C11_s1[w])+(cd-e)*sv[w]*E113_s1[w]);
    E123_s1[w]=1.0q/(1.0q-c2)*((E101_s1[w]-c*E011_s1[w]+c*powq(rfs[w],1.0q)*C11_s1[w]-powq(yv[w],1.0q)*A11_s1[w])+(ce-d)*sv[w]*E113_s1[w]);
    w=is2[j];
    E213_s2[w]=1.0q/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*powq(yv[w],1.0q)*A11_s2[w]-powq(rv[w],1.0q)*C11_s2[w]+E011_s2[w]-(dh+c)*E101_s2[w]-(hm+em-dm*(dh+c))*E113_s2[w]);
    E123_s2[w]=1.0q/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*powq(rv[w],1.0q)*C11_s2[w]-(oneph2p2eh)*powq(yv[w],1.0q)*A11_s2[w]+(oneph2p2eh)*E101_s2[w]-(dh+c)*E011_s2[w]+((dh+c)*(hm+em)-(oneph2p2eh)*dm)*E113_s2[w]);
  }

  for(w=0;w<8;w++) {
    AAS111[w] = h*A21_s2[w] +m* A11_s2[w] +s1*A11_s1[w];
    BBR111[w] = f*B21_r2[w] +g* B11_r2[w] +r1*B11_r1[w];
    D213[w]=1.0q/(1.0q-e2)*((D011[w]-e*D101[w]+e*AAS111[w]-BBR111[w])-(c-de)*yv[w]*D113[w]);
    D123[w]=1.0q/(1.0q-e2)*((D101[w]-e*D011[w]+e*BBR111[w]-AAS111[w])-(d-ce)*yv[w]*D113[w]);
  }

  for(j=0;j<4;j++) {
    w=ir1[j];
    F313_r1[w]=1.0q/(1.0q-d2)*((2.0q*F111_r1[w]-d*F201_r1[w]+d*powq(yv[w],1.0q)*B21_r1[w]-powq(sfr[w],2.0q)*C11_r1[w])+(cd-e)*rv[w]*F213_r1[w]);
    F223_r1[w]=1.0q/(1.0q-d2)*((F201_r1[w]-d*2.0q*F111_r1[w]+d*powq(sfr[w],2.0q)*C11_r1[w]-powq(yv[w],1.0q)*B21_r1[w])+(de-c)*rv[w]*F213_r1[w]);
    w=ir2[j];
    F313_r2[w]=1.0q/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*powq(yv[w],1.0q)*B21_r2[w]-powq(sv[w],2.0q)*C11_r2[w]+2.0q*F111_r2[w]-(cf+d)*F201_r2[w]-(fg+eg-cg*(cf+d))*F213_r2[w]);
    F223_r2[w]=1.0q/(onepf2p2ef-(cf+d)*(cf+d))*((cf+d)*powq(sv[w],2.0q)*C11_r2[w]-(onepf2p2ef)*powq(yv[w],1.0q)*B21_r2[w]+(onepf2p2ef)*F201_r2[w]-2.0q*(cf+d)*F111_r2[w]+((cf+d)*(fg+eg)-(onepf2p2ef)*cg)*F213_r2[w]);
    w=is1[j];
    E313_s1[w]=1.0q/(1.0q-c2)*((2.0q*E111_s1[w]-c*E201_s1[w]+c*powq(yv[w],1.0q)*A21_s1[w]-powq(rfs[w],2.0q)*C11_s1[w])+(cd-e)*sv[w]*E213_s1[w]);
    E223_s1[w]=1.0q/(1.0q-c2)*((E201_s1[w]-c*2.0q*E111_s1[w]+c*powq(rfs[w],2.0q)*C11_s1[w]-powq(yv[w],1.0q)*A21_s1[w])+(c*e-d)*sv[w]*E213_s1[w]);
    w=is2[j];
    E313_s2[w]=1.0q/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*powq(yv[w],1.0q)*A21_s2[w]-powq(rv[w],2.0q)*C11_s2[w]+2.0q*E111_s2[w]-(dh+c)*E201_s2[w]-(hm+em-dm*(dh+c))*E213_s2[w]);
    E223_s2[w]=1.0q/(oneph2p2eh-(dh+c)*(dh+c))*((dh+c)*powq(rv[w],2.0q)*C11_s2[w]-(oneph2p2eh)*powq(yv[w],1.0q)*A21_s2[w]+(oneph2p2eh)*E201_s2[w]-2.0q*(dh+c)*E111_s2[w]+((dh+c)*(hm+em)-(oneph2p2eh)*dm)*E213_s2[w]);
  }

  for(w=0;w<8;w++) {
    AAS211[w] = h*A31_s2[w] +m* A21_s2[w] +s1*A21_s1[w];
    AAS121[w] = h2*A31_s2[w] +2.0q*hm*A21_s2[w] +m2* A11_s2[w] +s1*s1*A11_s1[w];
    BBR211[w] = f2*B31_r2[w] +2.0q*fg*B21_r2[w] +g2* B11_r2[w] +r1*r1*B11_r1[w];
    BBR121[w] = f*B31_r2[w] +g* B21_r2[w] +r1*B21_r1[w];
    D313[w]=1.0q/(1.0q-e2)*((2.0q*D111[w]-e*D201[w]+e*AAS211[w]-BBR211[w])-(c-de)*yv[w]*D213[w]);

    D223_1[w]=1.0q/(1.0q-e2)*((D201[w]-e*2.0q*D111[w]+e*BBR211[w]-AAS211[w])-(d-ce)*yv[w]*D213[w]);
    D223_2[w]=1.0q/(1.0q-e2)*((D021[w]-e*2.0q*D111[w]+e*AAS121[w]-BBR121[w])-(c-de)*yv[w]*D123[w]);

    D133[w]=1.0q/(1.0q-e2)*((2.0q*D111[w]-e*D021[w]+e*BBR121[w]-AAS121[w])-(d-ce)*yv[w]*D123[w]);
    AAS201[w] = A21_s1[w]+A21_s2[w];
    AAS021[w] = h2*A21_s2[w] +2.0q*hm*A11_s2[w] +m2* A01_s2[w] +s1*s1*A01_s1[w];
    BBR201[w] = f2*B21_r2[w] +2.0q*fg*B11_r2[w] +g2* B01_r2[w] +r1*r1*B01_r1[w];
    BBR021[w] = +B21_r2[w] +B21_r1[w];
    D303[w]=1.0q/(1.0q-e2)*((2.0q*D101[w]+e*AAS201[w]-BBR201[w])-(c-de)*yv[w]*D203[w]);
    D033[w]=1.0q/(1.0q-e2)*((2.0q*D011[w]+e*BBR021[w]-AAS021[w])-(d-ce)*yv[w]*D023[w]);
  }

  for(w=0;w<8;w++) {
    D223[w]= 0.5q*D223_1[w] +0.5q*D223_2[w];
  }

  //
  // TRIPLE INTEGRALS
  //
  /* triple integrals are built from previous double integrals.
     Hijkl = /int/int/int r^i s^j y^k Ra^-l dr ds dy
     Since r2 is a fct(s) and s2 is a fct(r) in the case of a triangle, we introduce
     additional double integrals:
     EESijkl = /int/int r^i s(r)^j y^k Ra^-l dr dy = /int/int r^i (hr+m)^j y^k Ra^-l drdy - /int/int r^i sa^j y^k Ra^-l drdy
     FFRijkl = /int/int r(s)^i s^j y^k Ra^-l dr dy = /int/int (fs+g)^i s^j y^k Ra^-l drdy - /int/int r1^i s^j y^k Ra^-l drdy
  */

  /* SEED FUNCTION */
  /* since no analytical was found for the only tripple seed integral,
     we perform instead an adaptative Simpson quadrature */
  //  functionPtr = &FDr1s1;
  //  Dr1s1 = adaptiveSimpsons(functionPtr, y1, y2, 1e-9, 20);
  //  functionPtr = &FDr1s2;
  //  Dr1s2 = adaptiveSimpsons(functionPtr, y1, y2, 1e-9, 20);
  //  functionPtr = &FDr2s1;
  //  Dr2s1 = adaptiveSimpsons(functionPtr, y1, y2, 1e-9, 20);
  //  functionPtr = &FDr2s2;
  //  Dr2s2 = adaptiveSimpsons(functionPtr, y1, y2, 1e-9, 20);
  //
  //  H0003[0] = Dr2s2; H0003[1] = 0;
  //  H0003[2] = Dr2s1; H0003[3] = 0;
  //  H0003[4] = Dr1s2; H0003[5] = 0;
  //  H0003[6] = Dr1s1; H0003[7] = 0;

  H0003[0] = 0.0q; H0003[1] = 0.0q;
  H0003[2] = 0.0q; H0003[3] = 0.0q;
  H0003[4] = 0.0q; H0003[5] = 0.0q;
  H0003[6] = 0.0q; H0003[7] = 0.0q;

  for(w=0;w<8;w++) {
    FFR1001[w] = f*F101_r2[w] +g*F001_r2[w] +r1*F001_r1[w];
    EES0101[w] = h*E101_s2[w] +m*E001_s2[w] +s1*E001_s1[w];
    H0001[w] = 1.0q/(-2.0q)*(a2*H0003[w]-FFR1001[w]-EES0101[w]-powq(yv[w],1.0q)*D001[w]);
    FFR100m1[w] = f*F10m1_r2[w] +g*F00m1_r2[w] +r1*F00m1_r1[w];
    EES010m1[w] = h*E10m1_s2[w] +m*E00m1_s2[w] +s1*E00m1_s1[w];
    H000m1[w] = 1.0q/(-4.0q)*(-a2*H0001[w]-FFR100m1[w]-EES010m1[w]-powq(yv[w],1.0q)*D00m1[w]);
  }

  for(w=0;w<8;w++) {
    EES000m1[w] = E00m1_s1[w] +E00m1_s2[w];
    FFR000m1[w] = F00m1_r2[w] +F00m1_r1[w];
    H0011[w] = (1.0q/-1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -powq(yv[w],0.0q)*D00m1[w]) +(d-ce)*EES000m1[w] +(c-de)*FFR000m1[w]);
    H1001[w] = 1.0q/-1.0q/(1.0q-e2)*(+e*EES000m1[w] -FFR000m1[w] -1.0q*(de-c)*H0011[w]);
    H0101[w] = 1.0q/-1.0q/(1.0q-e2)*(-EES000m1[w] +e*FFR000m1[w] -1.0q*(ce-d)*H0011[w]);
  }

  for(w=0;w<8;w++) {
    EES100m1[w] = E10m1_s1[w] +E10m1_s2[w];
    FFR010m1[w] = F10m1_r2[w] +F10m1_r1[w];
    H1011[w] = (1.0q/-1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -powq(yv[w],0.0q)*D10m1[w]) -(c-de)*H000m1[w] +(d-ce)*EES100m1[w] +(c-de)*FFR100m1[w]);
    H0111[w] = (1.0q/-1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -powq(yv[w],0.0q)*D01m1[w]) -(d-ce)*H000m1[w] +(d-ce)*EES010m1[w] +(c-de)*FFR010m1[w]);
    H1101[w] = -1.0q/(1.0q-e2)*( -e*H000m1[w] -EES100m1[w] +e*FFR100m1[w] -(ce-d)*H1011[w]);
    //  H1101_2[w] = 1.0q/-1.0q/(1.0q-e2)*( -e*H000m1[w] +e*EES010m1[w] -FFR010m1[w] -(de-c)*H0111[w]);
    H2001[w] = 1.0q/-1.0q/(1.0q-e2)*(1.0q*H000m1[w] +e*EES100m1[w] -FFR100m1[w] -(de-c)*H1011[w]);
    H0201[w] = 1.0q/-1.0q/(1.0q-e2)*(1.0q*H000m1[w] -EES010m1[w] +e*FFR010m1[w] -(ce-d)*H0111[w]);
  }

  /*  for(w=0;w<8;w++) {
    H1101[w]= 0.5q*H1101_1[w] +0.5q*H1101_2[w];
  } */

  for(w=0;w<8;w++) {
    EES0001[w] = E001_s1[w] +E001_s2[w];
    FFR0001[w] = F001_r2[w] +F001_r1[w];
    EES1001[w] = E101_s1[w] +E101_s2[w];
    FFR0101[w] = F101_r2[w] +F101_r1[w];
    H0013[w] = (1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -powq(yv[w],0.0q)*D001[w]) +(d-ce)*EES0001[w] +(c-de)*FFR0001[w]);
    H1003[w] = 1.0q/(1.0q-e2)*( +e*EES0001[w] -FFR0001[w] +(de-c)*H0013[w]);
    H0103[w] = 1.0q/(1.0q-e2)*( -EES0001[w] +e*FFR0001[w] +(ce-d)*H0013[w]);
    H1013[w] = (1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -powq(yv[w],0.0q)*D101[w]) -(c-de)*H0001[w] +(d-ce)*EES1001[w] +(c-de)*FFR1001[w]);
    H0113[w] = (1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -powq(yv[w],0.0q)*D011[w])  -(d-ce)*H0001[w] +(d-ce)*EES0101[w] +(c-de)*FFR0101[w]);
    H1103_1[w] = 1.0q/(1.0q-e2)*( -e*H0001[w] +e*EES0101[w] -FFR0101[w] +(de-c)*H0113[w]);
    H1103_2[w] = 1.0q/(1.0q-e2)*( -e*H0001[w] +e*FFR1001[w] -EES1001[w] +(ce-d)*H1013[w]);
  }

  for(w=0;w<8;w++) {
    H1103[w]= 0.5q*H1103_1[w] +0.5q*H1103_2[w];
  }

  for(w=0;w<8;w++) {
    EES0011[w]= E011_s1[w] +E011_s2[w];
    FFR0011[w]= F011_r2[w] +F011_r1[w];
    H0023[w] = (1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(H0001[w] -powq(yv[w],1.0q)*D001[w]) +(d-ce)*EES0011[w] +(c-de)*FFR0011[w]);
    H2003[w] = 1.0q/(1.0q-e2)*(H0001[w] +e*EES1001[w] -FFR1001[w] +(de-c)*H1013[w]);
    H0203[w] = 1.0q/(1.0q-e2)*(H0001[w] -EES0101[w] +e*FFR0101[w] +(ce-d)*H0113[w]);
  }

  for(w=0;w<8;w++) {
    EES2001[w]= E201_s1[w] +E201_s2[w];
    FFR0201[w]= F201_r2[w] +F201_r1[w];
    FFR2001[w]= f2*F201_r2[w] +2.0q*fg*F101_r2[w] +g2*F001_r2[w] +r12*F001_r1[w];
    EES0201[w]= h2*E201_s2[w] +2.0q*hm*E101_s2[w] +m2*E001_s2[w] +s12*E001_s1[w];
    H2013[w] = (1.0q/1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -powq(yv[w],0.0q)*D201[w]) -2.0q*(c-de)*H1001[w] +(d-ce)*EES2001[w] +(c-de)*FFR2001[w]);
    H0213[w] = (1.0q/1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -powq(yv[w],0.0q)*D021[w]) -2.0q*(d-ce)*H0101[w] +(d-ce)*EES0201[w] +(c-de)*FFR0201[w]);
//    H0213[w] = (1.0q/1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -powq(yv[w],0.0q)*D021[w]) -2.0q*(d-ce)*H0101[w] +(c-de)*EES0201[w] +(d-ce)*FFR0201[w]);
  }

  for(w=0;w<8;w++) {
    EES1101[w]= h*E201_s2[w] +m*E101_s2[w] +s1*E101_s1[w];
    FFR1101[w]= f*F201_r2[w] +g*F101_r2[w] +r1*F101_r1[w];
    H1113[w] = (1.0q/1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -powq(yv[w],0.0q)*D111[w]) -(c-de)*H0101[w] -(d-ce)*H1001[w] +(d-ce)*EES1101[w] +(c-de)*FFR1101[w]);
    H2103[w] = 1.0q/(1.0q-e2)*(H0101[w] -e*H1001[w] +e*EES1101[w] -FFR1101[w] +(de-c)*H1113[w]);
    H1203[w] = 1.0q/(1.0q-e2)*(H1001[w] -e*H0101[w] +e*FFR1101[w] -EES1101[w]  +(ce-d)*H1113[w]);
    // H1203[w] = 1.0q/(1.0q-e2)*( -e*2.0q*H0101[w] +e*EES0201[w] -FFR0201[w] +(de-c)*H0213[w]);
  }

  for(w=0;w<8;w++) {
    EES0111[w] = h*E111_s2[w] +m*E011_s2[w] +s1*E011_s1[w];
    FFR1011[w] = f*F111_r2[w] +g*F011_r2[w] +r1*F011_r1[w];
    FFR0111[w] = F111_r2[w] +F111_r1[w];
    EES1011[w] = E111_s1[w] +E111_s2[w];
    H0123[w] = (1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(H0101[w] -powq(yv[w],1.0q)*D011[w]) -(d-ce)*H0011[w] +(d-ce)*EES0111[w] +(c-de)*FFR0111[w]);
    H1023[w] = (1.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(H1001[w] -powq(yv[w],1.0q)*D101[w]) -(c-de)*H0011[w] +(d-ce)*EES1011[w] +(c-de)*FFR1011[w]);
  }

  for(w=0;w<8;w++) {
    H3003[w] = 1.0q/(1.0q-e2)*(2.0q*H1001[w] +e*EES2001[w] -FFR2001[w] +(de-c)*H2013[w]);
    H0303[w] = 1.0q/(1.0q-e2)*(2.0q*H0101[w] -EES0201[w] +e*FFR0201[w] +(ce-d)*H0213[w]);
  }

  for(w=0;w<8;w++) {
    EES0003[w] = E003_s1[w] +E003_s2[w];
    FFR0003[w] = F003_r2[w] +F003_r1[w];
    H0015[w] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -powq(yv[w],0.0q)*D003[w]) +(d-ce)*EES0003[w] +(c-de)*FFR0003[w]);
    H1005[w] = 1.0q/3.0q/(1.0q-e2)*( +e*EES0003[w] -FFR0003[w] +3.0q*(de-c)*H0015[w]);
    H0105[w] = 1.0q/3.0q/(1.0q-e2)*( -EES0003[w] +e*FFR0003[w] +3.0q*(ce-d)*H0015[w]);
  }

  for(w=0;w<8;w++) {
    EES1003[w] = E103_s1[w] +E103_s2[w];
    FFR1003[w] = f*F103_r2[w] +g*F003_r2[w] +r1*F003_r1[w];
    EES0103[w] = h*E103_s2[w] +m*E003_s2[w] +s1*E003_s1[w];
    FFR0103[w] = F103_r2[w] +F103_r1[w];
    H1015[w] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -powq(yv[w],0.0q)*D103[w]) -(c-de)*H0003[w] +(d-ce)*EES1003[w] +(c-de)*FFR1003[w]);
    H0115[w] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -powq(yv[w],0.0q)*D013[w]) -(d-ce)*H0003[w] +(d-ce)*EES0103[w] +(c-de)*FFR0103[w]);
    H1105_1[w] = 1.0q/3.0q/(1.0q-e2)*( -e*H0003[w] +e*EES0103[w] -FFR0103[w] +3.0q*(de-c)*H0115[w]);
    H1105_2[w] = 1.0q/3.0q/(1.0q-e2)*( -e*H0003[w] +e*FFR1003[w] -EES1003[w] +3.0q*(ce-d)*H1015[w]);
  }

  for(w=0;w<8;w++) {
    H1105[w]= 0.5q*H1105_1[w] +0.5q*H1105_2[w];
  }

  for(w=0;w<8;w++) {
    EES0013[w] = E013_s1[w] +E013_s2[w];
    FFR0013[w] = F013_r2[w] +F013_r1[w];
    H0025[w] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(H0003[w] - powq(yv[w],1.0q)*D003[w]) +(d-ce)*EES0013[w] +(c-de)*FFR0013[w]);
    H2005[w] = 1.0q/3.0q/(1.0q-e2)*(H0003[w] +e*EES1003[w] -FFR1003[w] +3.0q*(de-c)*H1015[w]);
    H0205[w] = 1.0q/3.0q/(1.0q-e2)*(H0003[w] -EES0103[w] +e*FFR0103[w] +3.0q*(ce-d)*H0115[w]);
  }

  for(w=0;w<8;w++) {
    EES1103[w] = h*E203_s2[w] +m*E103_s2[w] +s1*E103_s1[w];
    FFR1103[w] = f*F203_r2[w] +g*F103_r2[w] +r1*F103_r1[w];
    H1115[w] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*( -powq(yv[w],0.0q)*D113[w]) -(c-de)*H0103[w] -(d-ce)*H1003[w] +(d-ce)*EES1103[w] +(c-de)*FFR1103[w]);
  }

  for(w=0;w<8;w++) {
    H2105[w] = 1.0q/3.0q/(1.0q-e2)*(H0103[w] -e*H1003[w] +e*EES1103[w] -FFR1103[w] +3.0q*(de-c)*H1115[w]);
    H1205[w] = 1.0q/3.0q/(1.0q-e2)*(H1003[w] -e*H0103[w] -EES1103[w] +e*FFR1103[w] +3.0q*(ce-d)*H1115[w]);
  }

  for(w=0;w<8;w++) {
    EES1013[w] = E113_s1[w] +E113_s2[w];
    FFR1013[w] = f*F113_r2[w] +g*F013_r2[w] +r1*F013_r1[w];
    H1025[w] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(H1003[w] -powq(yv[w],1.0q)*D103[w]) -(c-de)*H0013[w] +(d-ce)*EES1013[w] +(c-de)*FFR1013[w]);
    H2015[w] = 1.0q/3.0q/(1.0q-e2)*(H0013[w] +e*EES1013[w] -FFR1013[w] +3.0q*(de-c)*H1025[w]);

    EES0113[w] = h*E113_s2[w] +m*E013_s2[w] +s1*E013_s1[w];
    FFR0113[w] = F113_r2[w] +F113_r1[w];
    H0125[w] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(H0103[w] -powq(yv[w],1.0q)*D013[w]) -(d-ce)*H0013[w] +(d-ce)*EES0113[w] +(c-de)*FFR0113[w]);
    H0215[w] = 1.0q/3.0q/(1.0q-e2)*(H0013[w] -EES0113[w] +e*FFR0113[w] +3.0q*(ce-d)*H0125[w]);
  }

  for(w=0;w<8;w++) {
    EES2013[w] = E213_s1[w] +E213_s2[w];
    FFR2013[w] = f2*F213_r2[w] +2.0q*fg*F113_r2[w] +g2*F013_r2[w] +r12*F013_r1[w];
    EES0213[w] = h2*E213_s2[w] +2.0q*hm*E113_s2[w] +m2*E013_s2[w] +s12*E013_s1[w];
    FFR0213[w] = F213_r2[w] +F213_r1[w];
    H2025[w] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(H2003[w] -powq(yv[w],1.0q)*D203[w]) -2.0q*(c-de)*H1013[w] +(d-ce)*EES2013[w] +(c-de)*FFR2013[w]);
    H0225[w] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(H0203[w] -powq(yv[w],1.0q)*D023[w]) -2.0q*(d-ce)*H0113[w] +(d-ce)*EES0213[w] +(c-de)*FFR0213[w]);
  }

  for(w=0;w<8;w++) {
    EES1113[w] = h*E213_s2[w] +m*E113_s2[w] +s1*E113_s1[w];
    FFR1113[w] = f*F213_r2[w] +g*F113_r2[w] +r1*F113_r1[w];
    H1125[w] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(H1103[w] -powq(yv[w],1.0q)*D113[w]) -(c-de)*H0113[w] -(d-ce)*H1013[w] +(d-ce)*EES1113[w] +(c-de)*FFR1113[w]);
    H2115[w] = 1.0q/3.0q/(1.0q-e2)*(H0113[w] -e*H1013[w] +e*EES1113[w] -FFR1113[w] +3.0q*(de-c)*H1125[w]);
    H1215[w] = 1.0q/3.0q/(1.0q-e2)*(H1013[w] -e*H0113[w] -EES1113[w] +e*FFR1113[w] +3.0q*(ce-d)*H1125[w]);
  }

  for(w=0;w<8;w++) {
    H3015[w] = 1.0q/3.0q/(1.0q-e2)*(2.0q*H1013[w] +e*EES2013[w] -FFR2013[w] +3.0q*(de-c)*H2025[w]);
    H0315[w] = 1.0q/3.0q/(1.0q-e2)*(2.0q*H0113[w] -EES0213[w] +e*FFR0213[w] +3.0q*(ce-d)*H0225[w]);
  }

  for(w=0;w<8;w++) {
    EES2003[w] = E203_s1[w] +E203_s2[w];
    FFR0203[w] = F203_r2[w] +F203_r1[w];
    EES0203[w] = h2*E203_s2[w] +2.0q*hm*E103_s2[w] +m2*E003_s2[w] +s12*E003_s1[w];
    FFR2003[w] = f2*F203_r2[w] +2.0q*fg*F103_r2[w] +g2*F003_r2[w] +r12*F003_r1[w];
    H3005[w] = 1.0q/3.0q/(1.0q-e2)*(2.0q*H1003[w] +e*EES2003[w] -FFR2003[w] +3.0q*(de-c)*H2015[w]);
    H0305[w] = 1.0q/3.0q/(1.0q-e2)*(2.0q*H0103[w] -EES0203[w] +e*FFR0203[w] +3.0q*(ce-d)*H0215[w]);
  }

  for(w=0;w<8;w++) {
    EES3013[w] = E313_s1[w] +E313_s2[w];
    FFR0313[w] = F313_r2[w] +F313_r1[w];
    FFR3013[w] = f3*F313_r2[w] +3.0q*f2g*F213_r2[w] +3.0q*fg2*F113_r2[w] +g3*F013_r2[w] +r13*F013_r1[w];
    EES0313[w] = h3*E313_s2[w] +3.0q*h2m*E213_s2[w] +3.0q*hm2*E113_s2[w] +m3*E013_s2[w] +s13*E013_s1[w];
    H3025[w] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(H3003[w] -powq(yv[w],1.0q)*D303[w]) -3.0q*(c-de)*H2013[w] +(d-ce)*EES3013[w] +(c-de)*FFR3013[w]);
    H0325[w] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(H0303[w] -powq(yv[w],1.0q)*D033[w]) -3.0q*(d-ce)*H0213[w] +(d-ce)*EES0313[w] +(c-de)*FFR0313[w]);
  }

  for(w=0;w<8;w++) {
    EES1123[w] = h*E223_s2[w] +m*E123_s2[w] +s1*E123_s1[w];
    FFR1123[w] = f*F223_r2[w] +g*F123_r2[w] +r1*F123_r1[w];
    H3115[w] = 1.0q/3.0q/(1.0q-e2)*(-e*3.0q*H2013[w] -EES3013[w] +e*FFR3013[w] +3.0q*(ce-d)*H3025[w]);
    H1315[w] = 1.0q/3.0q/(1.0q-e2)*(-e*3.0q*H0213[w] +e*EES0313[w] -FFR0313[w] +3.0q*(de-c)*H0325[w]);
    H1135[w] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(2.0q*H1113[w] -powq(yv[w],2.0q)*D113[w]) -(c-de)*H0123[w] -(d-ce)*H1023[w] +(d-ce)*EES1123[w] +(c-de)*FFR1123[w]);
  }

  for(w=0;w<8;w++) {
    EES2113[w] = h*E313_s2[w] +m*E213_s2[w] +s1*E213_s1[w];
    FFR1213[w] = f*F313_r2[w] +g*F213_r2[w] +r1*F213_r1[w];
    FFR2113[w] = f2*F313_r2[w] +2.0q*fg*F213_r2[w] +g2*F113_r2[w] +r12*F113_r1[w];
    EES1213[w] = h2*E313_s2[w] +2.0q*hm*E213_s2[w] +m2*E113_s2[w] +s12*E113_s1[w];
    H2125[w] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(H2103[w] -powq(yv[w],1.0q)*D213[w]) -2.0q*(c-de)*H1113[w] -(d-ce)*H2013[w] +(d-ce)*EES2113[w] +(c-de)*FFR2113[w]);
    H1225[w] = (1.0q/3.0q/(1.0q-e2-c2-d2+2.0q*cde))*((1.0q-e2)*(H1203[w] -powq(yv[w],1.0q)*D123[w]) -(c-de)*H0213[w] -2.0q*(d-ce)*H1113[w] +(d-ce)*EES1213[w] +(c-de)*FFR1213[w]);
    H2215_1[w] = 1.0q/3.0q/(1.0q-e2)*(H0213[w] -e*2.0q*H1113[w] +e*EES1213[w] -FFR1213[w] +3.0q*(de-c)*H1225[w]);
    H2215_2[w] = 1.0q/3.0q/(1.0q-e2)*(H2013[w] -e*2.0q*H1113[w] +e*FFR2113[w] -EES2113[w] +3.0q*(ce-d)*H2125[w]);
  }

  for(w=0;w<8;w++) {
    H2215[w]= 0.5q*H2215_1[w] +0.5q*H2215_2[w];
  }

  for(w=0;w<8;w++) {
    H4015[w] = 1.0q/3.0q/(1.0q-e2)*(3.0q*H2013[w] +e*EES3013[w] -FFR3013[w] +3.0q*(de-c)*H3025[w]);
    H0415[w] = 1.0q/3.0q/(1.0q-e2)*(3.0q*H0213[w] -EES0313[w] +e*FFR0313[w] +3.0q*(ce-d)*H0325[w]);
  }

  for(i=0;i<3;i++) { quadmath_snprintf(cv3[i], 256,  "%*.30Qf", p[i]); } printf("p =  %s %s %s  \n", cv3[0], cv3[1], cv3[2]);
  for(i=0;i<3;i++) { quadmath_snprintf(cv3[i], 256,  "%*.30Qf", q[i]); } printf("q =  %s %s %s  \n", cv3[0], cv3[1], cv3[2]);
  for(i=0;i<3;i++) { quadmath_snprintf(cv3[i], 256,  "%*.30Qf", t[i]); } printf("t =  %s %s %s  \n", cv3[0], cv3[1], cv3[2]);

  quadmath_snprintf(cv3[0], 256,  "%*.30Qf", alpha); printf("alpha = %s \n\n ", cv3[0]);

  for(i=0;i<3;i++) { quadmath_snprintf(cv3[i], 256,  "%*.30Qf", n[i]); } printf("n =  %s %s %s  \n", cv3[0], cv3[1], cv3[2]);
  for(i=0;i<3;i++) { quadmath_snprintf(cv3[i], 256,  "%*.30Qf", pxt[i]); } printf("pxt =  %s %s %s  \n", cv3[0], cv3[1], cv3[2]);
  for(i=0;i<3;i++) { quadmath_snprintf(cv3[i], 256,  "%*.30Qf", qxt[i]); } printf("qxt =  %s %s %s  \n\n", cv3[0], cv3[1], cv3[2]);

  quadmath_snprintf(cv3[0], 256,  "%*.30Qf", tdotn); printf("tdotn = %s \n", cv3[0]);
  quadmath_snprintf(cv3[0], 256,  "%*.30Qf", a2); printf("a2 = %s \n", cv3[0]);
  quadmath_snprintf(cv3[0], 256,  "%*.30Qf", c2);
  quadmath_snprintf(cv3[1], 256,  "%*.30Qf", d2);
  quadmath_snprintf(cv3[2], 256,  "%*.30Qf", e2); printf("c2 d2 e2 =  %s %s %s  \n\n", cv3[0], cv3[1], cv3[2]);

  quadmath_snprintf(cv3[0], 256,  "%*.30Qf", y1); quadmath_snprintf(cv3[1], 256,  "%*.30Qf", y2); printf("y1 y2 =  %s %s  \n", cv3[0], cv3[1]);
  quadmath_snprintf(cv3[0], 256,  "%*.30Qf", r1); quadmath_snprintf(cv3[1], 256,  "%*.30Qf", r2); printf("r1 r2 =  %s %s  \n", cv3[0], cv3[1]);
  quadmath_snprintf(cv3[0], 256,  "%*.30Qf", s1); quadmath_snprintf(cv3[1], 256,  "%*.30Qf", s2); printf("s1 s2 =  %s %s  \n", cv3[0], cv3[1]);

  printf("Ra =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", Ra[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("rRa =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", rRa[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("sRa =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", sRa[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("Rar1 =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", Rar1[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("Rar2 =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", Rar2[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("Ras1 =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", Ras1[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("Ras2 =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", Ras2[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");


  printf("Rdotp =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", Rdotp[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("Rdotq =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", Rdotq[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("rRdott =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", rRdott[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("sRdott =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", sRdott[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("RaRdotp =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", RaRdotp[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("RaRdotq =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", RaRdotq[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("rRaRdott =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", rRaRdott[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("sRaRdott =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", sRaRdott[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");

  printf("A01 =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", A01_s1[i]+A01_s2[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("B01 =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", B01_r1[i]+B01_r2[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("C01_s =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", C01_s1[i]+C01_s2[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("C01_r =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", C01_r1[i]+C01_r2[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");

  printf("A11 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",A11_s1[i]+A11_s2[i]);} printf("\n");
  printf("B11 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",B11_r1[i]+B11_r2[i]);} printf("\n");
  printf("C11_s =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",C11_s1[i]+C11_s2[i]);} printf("\n");
  printf("C11_r =\n");for(i=0;i<8;i++) {printf("  %19.15e\n",C11_r1[i]+C11_r2[i]);} printf("\n\n");

  printf("A0m1 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",A0m1_s1[i]+A0m1_s2[i]);} printf("\n");
  printf("B0m1 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",B0m1_r1[i]+B0m1_r2[i]);} printf("\n");
  printf("C0m1_s =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",C0m1_s1[i]+C0m1_s2[i]);} printf("\n");
  printf("C0m1_r =\n");for(i=0;i<8;i++) {printf("  %19.15e\n",C0m1_r1[i]+C0m1_r2[i]);} printf("\n\n");

  printf("A1m1 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",A1m1_s1[i]+A1m1_s2[i]);} printf("\n");
  printf("B1m1 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",B1m1_r1[i]+B1m1_r2[i]);} printf("\n");
  printf("C1m1_s =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",C1m1_s1[i]+C1m1_s2[i]);} printf("\n");
  printf("C1m1_r =\n");for(i=0;i<8;i++) {printf("  %19.15e\n",C1m1_r1[i]+C1m1_r2[i]);} printf("\n\n");

  printf("A0m3 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",A0m3_s1[i]+A0m3_s2[i]);} printf("\n");
  printf("B0m3 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",B0m3_r1[i]+B0m3_r2[i]);} printf("\n");
  printf("C0m3_s =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",C0m3_s1[i]+C0m3_s2[i]);} printf("\n");
  printf("C0m3_r =\n");for(i=0;i<8;i++) {printf("  %19.15e\n",C0m3_r1[i]+C0m3_r2[i]);} printf("\n\n");

  printf("A1m3 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",A1m3_s1[i]+A1m3_s2[i]);} printf("\n");
  printf("B1m3 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",B1m3_r1[i]+B1m3_r2[i]);} printf("\n\n");

  printf("A2m1 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",A2m1_s1[i]+A2m1_s2[i]);} printf("\n");
  printf("B2m1 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",B2m1_r1[i]+B2m1_r2[i]);} printf("\n\n");

  printf("A21 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",A21_s1[i]+A21_s2[i]);} printf("\n");
  printf("B21 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",B21_r1[i]+B21_r2[i]);} printf("\n\n");

  printf("A31 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",A31_s1[i]+A31_s2[i]);} printf("\n");
  printf("B31 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",B31_r1[i]+B31_r2[i]);} printf("\n\n");

  printf("D003 =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", D003[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("E003 =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", E003_s1[i]+E003_s2[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("F003 =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", F003_r1[i]+F003_r2[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");

  printf("D103 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D103[i]);} printf("\n");
  printf("D013 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D013[i]);} printf("\n");
  printf("E103 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E103_s1[i]+E103_s2[i]);} printf("\n");
  printf("F103 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F103_r1[i]);} printf("\n");
  printf("F103 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F103_r2[i]);} printf("\n");

  printf("E001 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E001_s1[i]+E001_s2[i]);} printf("\n");
  printf("E101 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E101_s1[i]+E101_s2[i]);} printf("\n");
  printf("E011 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E011_s1[i]+E011_s2[i]);} printf("\n");

  printf("F001 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F001_r1[i]+F001_r2[i]);} printf("\n");
  printf("F101 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F101_r1[i]+F101_r2[i]);} printf("\n");
  printf("F011 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F011_r1[i]+F011_r2[i]);} printf("\n\n");

  printf("D001 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D001[i]);} printf("\n");
  printf("D101 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D101[i]);} printf("\n");
  printf("D011 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D011[i]);} printf("\n");

  printf("D00m1 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D00m1[i]);} printf("\n");
  printf("E00m1 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E00m1_s1[i]+E00m1_s2[i]);} printf("\n");
  printf("F00m1 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F00m1_r1[i]+F00m1_r2[i]);} printf("\n\n");

  printf("D00m3 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D00m3[i]);} printf("\n");
  printf("E00m3 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E00m3_s1[i]+E00m3_s2[i]);} printf("\n");
  printf("F00m3 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F00m3_r1[i]+F00m3_r2[i]);} printf("\n\n");

  printf("E10m1 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E10m1_s1[i]+E10m1_s2[i]);} printf("\n");
  printf("F10m1 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F10m1_r1[i]+F10m1_r2[i]);} printf("\n\n");

  printf("E201 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E201_s1[i]+E201_s2[i]);} printf("\n");
  printf("F201 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F201_r1[i]+F201_r2[i]);} printf("\n\n");

  printf("E111 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E111_s1[i]+E111_s2[i]);} printf("\n");
  printf("F111 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F111_r1[i]+F111_r2[i]);} printf("\n\n");

  printf("D10m1 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D10m1[i]);} printf("\n");
  printf("D01m1 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D01m1[i]);} printf("\n");

  printf("D201 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D201[i]);} printf("\n");
  printf("D021 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D021[i]);} printf("\n");
  printf("D111 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D111[i]);} printf("\n\n");

  printf("E203 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E203_s1[i]+E203_s2[i]);} printf("\n");
  printf("F203 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F203_r1[i]+F203_r2[i]);} printf("\n\n");

  printf("D203 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D203[i]);} printf("\n");
  printf("D023 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D023[i]);} printf("\n");

  printf("E013 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E013_s1[i]+E013_s2[i]);} printf("\n");
  printf("F013 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F013_r1[i]+F013_r2[i]);} printf("\n\n");
  printf("E023 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E023_s1[i]+E023_s2[i]);} printf("\n");
  printf("F023 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F023_r1[i]+F023_r2[i]);} printf("\n\n");

  printf("E113 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E113_s1[i]+E113_s2[i]);} printf("\n");
  printf("F113 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F113_r1[i]+F113_r2[i]);} printf("\n\n");
  printf("D113 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D113[i]);} printf("\n");

  printf("E123 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E123_s1[i]+E123_s2[i]);} printf("\n");
  printf("F123 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F123_r1[i]+F123_r2[i]);} printf("\n\n");
  printf("E213 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E213_s1[i]+E213_s2[i]);} printf("\n");
  printf("F213 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F213_r1[i]+F213_r2[i]);} printf("\n\n");

  printf("D123 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D123[i]);} printf("\n");
  printf("D213 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D213[i]);} printf("\n");

  printf("E313 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E313_s1[i]+E313_s2[i]);} printf("\n");
  printf("F313 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F313_r1[i]+F313_r2[i]);} printf("\n\n");
  printf("E223 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",E223_s1[i]+E223_s2[i]);} printf("\n");
  printf("F223 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",F223_r1[i]+F223_r2[i]);} printf("\n\n");

  printf("D313 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D313[i]);} printf("\n");
  printf("D223 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D223[i]);} printf("\n");
  printf("D133 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D133[i]);} printf("\n");
  printf("D303 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D303[i]);} printf("\n");
  printf("D033 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",D033[i]);} printf("\n");

  printf("H0003 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0003[i]);} printf("\n");

  printf("H0001 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0001[i]);} printf("\n");
  printf("H000m1 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H000m1[i]);} printf("\n");

  printf("H0011 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0011[i]);} printf("\n");
  printf("H1001 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1001[i]);} printf("\n");
  printf("H0101 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0101[i]);} printf("\n");

  printf("H1011 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1011[i]);} printf("\n");
  printf("H0111 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0111[i]);} printf("\n");
  /* printf("H1101 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1101[i]);} printf("\n"); */
  printf("H2001 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H2001[i]);} printf("\n");
  printf("H0201 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0201[i]);} printf("\n");

  printf("H0013 =\n");for(i=0;i<8;i++) {quadmath_snprintf(cv8[i], 256,  "%*.30Qf", H0013[i]); }; for(i=0;i<8;i++) {printf(" %s \n",cv8[i]);} printf("\n");
  printf("H0013 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0013[i]);} printf("\n");
  printf("H1003 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1003[i]);} printf("\n");
  printf("H0103 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0103[i]);} printf("\n");
  printf("H1013 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1013[i]);} printf("\n");
  printf("H0113 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0113[i]);} printf("\n");
  printf("H1103 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1103[i]);} printf("\n");

  printf("H0023 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0023[i]);} printf("\n");
  printf("H2003 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H2003[i]);} printf("\n");
  printf("H0203 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0203[i]);} printf("\n");

  printf("H2013 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H2013[i]);} printf("\n");
  printf("H0213 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0213[i]);} printf("\n");
  printf("H1113 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1113[i]);} printf("\n");

  printf("H2103 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H2103[i]);} printf("\n");
  printf("H1203 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1203[i]);} printf("\n");

  printf("H3003 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H3003[i]);} printf("\n");
  printf("H0303 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0303[i]);} printf("\n");

  printf("H0015 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0015[i]);} printf("\n");
  printf("H1005 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1005[i]);} printf("\n");
  printf("H0105 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0105[i]);} printf("\n");

  printf("H1015 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1015[i]);} printf("\n");
  printf("H0115 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0115[i]);} printf("\n");
  printf("H1105 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1105[i]);} printf("\n");

  printf("H0025 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0025[i]);} printf("\n");
  printf("H2005 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H2005[i]);} printf("\n");
  printf("H0205 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0205[i]);} printf("\n");

  printf("H0123 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0123[i]);} printf("\n");
  printf("H1023 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1023[i]);} printf("\n");

  printf("H1115 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1115[i]);} printf("\n");
  printf("H2105 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H2105[i]);} printf("\n");
  printf("H1205 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1205[i]);} printf("\n");

  printf("H1025 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1025[i]);} printf("\n");
  printf("H2015 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H2015[i]);} printf("\n");
  printf("H0125 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0125[i]);} printf("\n");

  printf("H0215 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0215[i]);} printf("\n");
  printf("H2025 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H2025[i]);} printf("\n");
  printf("H0225 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H0225[i]);} printf("\n");

  printf("H1125 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1125[i]);} printf("\n");
  printf("H2115 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H2115[i]);} printf("\n");
  printf("H1215 =\n");for(i=0;i<8;i++) {printf("  %19.15e \n",H1215[i]);} printf("\n");


  signv[0] =  1.0q;
  signv[1] = -1.0q;
  signv[2] = -1.0q;
  signv[3] =  1.0q;
  signv[4] = -1.0q;
  signv[5] =  1.0q;
  signv[6] =  1.0q;
  signv[7] = -1.0q;

  //
  // FINAL EVALUATION OF TRIPLE INTEGRALS
  //
  for(w=0;w<8;w++) {
    sch[0]= sch[0] +H0013[w]*signv[w];
    sch[1]= sch[1] +H1003[w]*signv[w];
    sch[2]= sch[2] +H0103[w]*signv[w];
    sch[3]= sch[3] +H2003[w]*signv[w];
    sch[4]= sch[4] +H0203[w]*signv[w];

    sch[5]= sch[5] +H1013[w]*signv[w];
    sch[6]= sch[6] +H0113[w]*signv[w];
    sch[7]= sch[7] +H1103[w]*signv[w];
    sch[8]= sch[8] +H2013[w]*signv[w];
    sch[9]= sch[9] +H0213[w]*signv[w];

    sch[10]= sch[10] +H2103[w]*signv[w];
    sch[11]= sch[11] +H1203[w]*signv[w];
    sch[12]= sch[12] +H1113[w]*signv[w];
    sch[13]= sch[13] +H3003[w]*signv[w];
    sch[14]= sch[14] +H0303[w]*signv[w];

    sch[15]= sch[15] +H0015[w]*signv[w];
    sch[16]= sch[16] +H1005[w]*signv[w];
    sch[17]= sch[17] +H0105[w]*signv[w];
    sch[18]= sch[18] +H2005[w]*signv[w];
    sch[19]= sch[19] +H0205[w]*signv[w];

    sch[20]= sch[20] +H1015[w]*signv[w];
    sch[21]= sch[21] +H0115[w]*signv[w];
    sch[22]= sch[22] +H1105[w]*signv[w];
    sch[23]= sch[23] +H1115[w]*signv[w];
    sch[24]= sch[24] +H1125[w]*signv[w];

    sch[25]= sch[25] +H2115[w]*signv[w];
    sch[26]= sch[26] +H1215[w]*signv[w];
    sch[27]= sch[27] +H1025[w]*signv[w];
    sch[28]= sch[28] +H0125[w]*signv[w];
    sch[29]= sch[29] +H2015[w]*signv[w];

    sch[30]= sch[30] +H0215[w]*signv[w];
    sch[31]= sch[31] +H2025[w]*signv[w];
    sch[32]= sch[32] +H0225[w]*signv[w];
    sch[33]= sch[33] +H3015[w]*signv[w];
    sch[34]= sch[34] +H0315[w]*signv[w];

    sch[35]= sch[35] +H2105[w]*signv[w];
    sch[36]= sch[36] +H1205[w]*signv[w];
    sch[37]= sch[37] +H3005[w]*signv[w];
    sch[38]= sch[38] +H0305[w]*signv[w];
    sch[39]= sch[39] +H2215[w]*signv[w];

    sch[40]= sch[40] +H2125[w]*signv[w];
    sch[41]= sch[41] +H1225[w]*signv[w];
    sch[42]= sch[42] +H4015[w]*signv[w];
    sch[43]= sch[43] +H0415[w]*signv[w];
    sch[44]= sch[44] +H3025[w]*signv[w];

    sch[45]= sch[45] +H0325[w]*signv[w];
    sch[46]= sch[46] +H3115[w]*signv[w];
    sch[47]= sch[47] +H1315[w]*signv[w];
  }

  quadmath_snprintf(csca, 256,  "%*.30Qf", sch[0]); printf("sch[0] = %s \n\n ", csca);
  printf("sch[0] =\n"); printf("  %19.15e \n",sch[0]);
  printf("sch[1] =\n"); printf("  %19.15e \n",sch[1]);
  printf("sch[2] =\n"); printf("  %19.15e \n",sch[2]);

  printf("sch[3] =\n"); printf("  %19.15e \n",sch[3]);
  printf("sch[4] =\n"); printf("  %19.15e \n",sch[4]);

  printf("sch[5] =\n"); printf("  %19.15e \n",sch[5]);
  printf("sch[6] =\n"); printf("  %19.15e \n",sch[6]);
  printf("sch[7] =\n"); printf("  %19.15e \n",sch[7]);

  printf("sch[8] =\n"); printf("  %19.15e \n",sch[8]);
  printf("sch[9] =\n"); printf("  %19.15e \n",sch[9]);
  printf("sch[10] =\n"); printf("  %19.15e \n",sch[10]);
  printf("sch[11] =\n"); printf("  %19.15e \n",sch[11]);

  printf("sch[12] =\n"); printf("  %19.15e \n",sch[12]);

  printf("sch[13] =\n"); printf("  %19.15e \n",sch[13]);
  printf("sch[14] =\n"); printf("  %19.15e \n",sch[14]);

  printf("sch[15] =\n"); printf("  %19.15e \n",sch[15]);
  printf("sch[16] =\n"); printf("  %19.15e \n",sch[16]);
  printf("sch[17] =\n"); printf("  %19.15e \n",sch[17]);

  printf("sch[18] =\n"); printf("  %19.15e \n",sch[18]);
  printf("sch[19] =\n"); printf("  %19.15e \n",sch[19]);

  printf("sch[20] =\n"); printf("  %19.15e \n",sch[20]);
  printf("sch[21] =\n"); printf("  %19.15e \n",sch[21]);
  printf("sch[22] =\n"); printf("  %19.15e \n",sch[22]);

  printf("sch[23] =\n"); printf("  %19.15e \n",sch[23]);

  printf("sch[24] =\n"); printf("  %19.15e \n",sch[24]);
  printf("sch[25] =\n"); printf("  %19.15e \n",sch[25]);
  printf("sch[26] =\n"); printf("  %19.15e \n",sch[26]);

  printf("sch[27] =\n"); printf("  %19.15e \n",sch[27]);
  printf("sch[28] =\n"); printf("  %19.15e \n",sch[28]);
  printf("sch[29] =\n"); printf("  %19.15e \n",sch[29]);
  printf("sch[30] =\n"); printf("  %19.15e \n",sch[30]);

  printf("sch[31] =\n"); printf("  %19.15e \n",sch[31]);
  printf("sch[32] =\n"); printf("  %19.15e \n",sch[32]);
  printf("sch[33] =\n"); printf("  %19.15e \n",sch[33]);
  printf("sch[34] =\n"); printf("  %19.15e \n",sch[34]);

  printf("sch[29] =\n"); printf("  %19.15e \n",sch[29]);
  printf("sch[30] =\n"); printf("  %19.15e \n",sch[30]);
  printf("sch[35] =\n"); printf("  %19.15e \n",sch[35]);
  printf("sch[36] =\n"); printf("  %19.15e \n",sch[36]);

  printf("sch[37] =\n"); printf("  %19.15e \n",sch[37]);
  printf("sch[38] =\n"); printf("  %19.15e \n",sch[38]);

  printf("sch[39] =\n"); printf("  %19.15e \n",sch[39]);
  printf("sch[40] =\n"); printf("  %19.15e \n",sch[40]);
  printf("sch[41] =\n"); printf("  %19.15e \n",sch[41]);

  printf("sch[42] =\n"); printf("  %19.15e \n",sch[42]);
  printf("sch[43] =\n"); printf("  %19.15e \n",sch[43]);

  printf("sch[44] =\n"); printf("  %19.15e \n",sch[44]);
  printf("sch[45] =\n"); printf("  %19.15e \n",sch[45]);

  printf("sch[46] =\n"); printf("  %19.15e \n",sch[46]);
  printf("sch[47] =\n"); printf("  %19.15e \n",sch[47]);

  //
  // NODAL FORCES FINAL EVALUATION
  //
  /* TAKEN CARE OF THIS
  // Addition vector definitions
  b[0] = bx ;
  b[1] = by ;
  b[2] = bz ;

  txb[0] = t[1]*b[2] - t[2]*b[1];
  txb[1] = t[2]*b[0] - t[0]*b[2];
  txb[2] = t[0]*b[1] - t[1]*b[0];

  pxb[0] = p[1]*b[2] - p[2]*b[1];
  pxb[1] = p[2]*b[0] - p[0]*b[2];
  pxb[2] = p[0]*b[1] - p[1]*b[0];

  qxb[0] = q[1]*b[2] - q[2]*b[1];
  qxb[1] = q[2]*b[0] - q[0]*b[2];
  qxb[2] = q[0]*b[1] - q[1]*b[0];

  bxt[0] = b[1]*t[2] - b[2]*t[1];
  bxt[1] = b[2]*t[0] - b[0]*t[2];
  bxt[2] = b[0]*t[1] - b[1]*t[0];

  tdotn= t[0]*n[0] +t[1]*n[1] +t[2]*n[2];
  txbdotn= txb[0]*n[0] +txb[1]*n[1] +txb[2]*n[2];
  bxtdotn= bxt[0]*n[0] +bxt[1]*n[1] +bxt[2]*n[2];
  pxbdotn= pxb[0]*n[0] +pxb[1]*n[1] +pxb[2]*n[2];
  qxbdotn= qxb[0]*n[0] +qxb[1]*n[1] +qxb[2]*n[2];

  pxbdott= pxb[0]*t[0] +pxb[1]*t[1] +pxb[2]*t[2];
  qxbdott= qxb[0]*t[0] +qxb[1]*t[1] +qxb[2]*t[2];

  for(i=0;i<3;i++) {
    ttbn[i]= t[i]*txbdotn;
    tbtn[i]= t[i]*bxtdotn;
    tpbn[i]= t[i]*pxbdotn;
    tqbn[i]= t[i]*qxbdotn;
  }
  // Important constant
  factor=MU*alpha*0.25q/M_PI/(1.0q-NU);

  // Vectors that are common to all nodal forces
  for(i=0;i<3;i++) {
    I0013[i]=(txb[i]*tdotn+ttbn[i])*(1.0q-NU) + (bxt[i]*tdotn+tbtn[i]);
    I0103[i]=(qxb[i]*tdotn+tqbn[i])*(1.0q-NU)-(qxbdott*n[i])+(q[i]*bxtdotn);
    I1003[i]=(pxb[i]*tdotn+tpbn[i])*(1.0q-NU)-(pxbdott*n[i])+(p[i]*bxtdotn);

    I0015[i]=(txb[i]*tdotn+ttbn[i])*1.5q*(1.0q-NU)*a2;
    I0105[i]=(qxb[i]*tdotn+tqbn[i])*1.5q*(1.0q-NU)*a2-(qxbdott*n[i]) * 3.0q*a2;
    I1005[i]=(pxb[i]*tdotn+tpbn[i])*1.5q*(1.0q-NU)*a2-(pxbdott*n[i]) * 3.0q*a2;

    I0215[i]=-(qxbdott*(tdotn))*q[i] * 3.0q;
    I2015[i]=-(pxbdott*(tdotn))*p[i] * 3.0q;
    I0125[i]=-(qxbdott*tdotn)*t[i] * 3.0q;
    I1025[i]=-(pxbdott*tdotn)*t[i] * 3.0q;

    I1115[i]=-((pxbdott*tdotn)*q[i] +(qxbdott*tdotn)*p[i])*3.0q;
  }
  */
  //
  // Fx4
  //
  // constants associated to the node shape function
  AA=2.0q;
  BB=-(3.0q*r1+r2);
  CC=r1*r1+r1*r2;

  // triple integrals required at that node
  F0013=AA*sch[8]+BB*sch[5]+CC*sch[0];
  F0103=AA*sch[10]+BB*sch[7]+CC*sch[2];
  F1003=AA*sch[13]+BB*sch[3]+CC*sch[1];

  F0015=AA*sch[29]+BB*sch[20]+CC*sch[15];
  F0105=AA*sch[35]+BB*sch[22]+CC*sch[17];
  F1005=AA*sch[37]+BB*sch[18]+CC*sch[16];

  F0215=AA*sch[39]+BB*sch[26]+CC*sch[30];
  F2015=AA*sch[42]+BB*sch[33]+CC*sch[29];
  F0125=AA*sch[40]+BB*sch[24]+CC*sch[28];
  F1025=AA*sch[44]+BB*sch[31]+CC*sch[27];

  F1115=AA*sch[46]+BB*sch[25]+CC*sch[23];

  for(i=0;i<3;i++) {
    fLLprime[i]= I0013[i] * F0013 + I1003[i] * F1003 + I0103[i] * F0103
      + I0015[i] * F0015 + I1005[i] * F1005 + I0105[i] * F0105
      + I1025[i] * F1025 + I0125[i] * F0125 + I1115[i] * F1115
      + I2015[i] * F2015 + I0215[i] * F0215;

    fx4[i]= fLLprime[i]*factor*one_o_dr*one_o_dr;
  }
  //
  // Fx5
  //
  // constants associated to the node shape function
  AA=2.0q;
  BB=-(3.0q*s1+s2);
  CC=s1*s1+s1*s2;

  // triple integrals required at that node
  F0013=AA*sch[9]+BB*sch[6]+CC*sch[0];
  F0103=AA*sch[14]+BB*sch[4]+CC*sch[2];
  F1003=AA*sch[11]+BB*sch[7]+CC*sch[1];

  F0015=AA*sch[30]+BB*sch[21]+CC*sch[15];
  F0105=AA*sch[38]+BB*sch[19]+CC*sch[17];
  F1005=AA*sch[36]+BB*sch[22]+CC*sch[16];
  F0215=AA*sch[43]+BB*sch[34]+CC*sch[30];
  F2015=AA*sch[39]+BB*sch[25]+CC*sch[29];
  F0125=AA*sch[45]+BB*sch[32]+CC*sch[28];
  F1025=AA*sch[41]+BB*sch[24]+CC*sch[27];
  F1115=AA*sch[47]+BB*sch[26]+CC*sch[23];

  for(i=0;i<3;i++) {
    fLLprime[i]= I0013[i] * F0013 + I1003[i] * F1003 + I0103[i] * F0103
      + I0015[i] * F0015 + I1005[i] * F1005 + I0105[i] * F0105
      + I1025[i] * F1025 + I0125[i] * F0125 + I1115[i] * F1115
      + I2015[i] * F2015 + I0215[i] * F0215;

    fx5[i]= fLLprime[i]*factor*one_o_ds*one_o_ds;
  }
  //
  // Fx3
  //
  // constants associated to the node shape function
  AA=2.0q*one_o_dr*one_o_dr;
  BB=2.0q*one_o_ds*one_o_ds;
  CC=4.0q*one_o_dr*one_o_ds;
  DD=-4.0q*r1*one_o_dr*one_o_dr-3.0q*one_o_dr-4.0q*s1*one_o_dr*one_o_ds;
  EE=-4.0q*s1*one_o_ds*one_o_ds-3.0q*one_o_ds-4.0q*r1*one_o_dr*one_o_ds;
  FF=1.0q+4.0q*r1*s1*one_o_dr*one_o_ds+3.0q*r1*one_o_dr+3.0q*s1*one_o_ds+2.0q*r1*r1*one_o_dr*one_o_dr+2.0q*s1*s1*one_o_ds*one_o_ds;

  // triple integrals required at that node
  F0013=AA*sch[8]+DD*sch[5]+BB*sch[9]+EE*sch[6]+CC*sch[12]+FF*sch[0];
  F0103=AA*sch[10]+DD*sch[7]+BB*sch[14]+EE*sch[4]+CC*sch[11]+FF*sch[2];
  F1003=AA*sch[13]+DD*sch[3]+BB*sch[11]+EE*sch[7]+CC*sch[10]+FF*sch[1];

  F0015=AA*sch[29]+DD*sch[20]+BB*sch[30]+EE*sch[21]+CC*sch[23]+FF*sch[15];
  F0105=AA*sch[35]+DD*sch[22]+BB*sch[38]+EE*sch[19]+CC*sch[36]+FF*sch[17];
  F1005=AA*sch[37]+DD*sch[18]+BB*sch[36]+EE*sch[22]+CC*sch[35]+FF*sch[16];
  F0215=AA*sch[39]+DD*sch[26]+BB*sch[43]+EE*sch[34]+CC*sch[47]+FF*sch[30];
  F2015=AA*sch[42]+DD*sch[33]+BB*sch[39]+EE*sch[25]+CC*sch[46]+FF*sch[29];
  F0125=AA*sch[40]+DD*sch[24]+BB*sch[45]+EE*sch[32]+CC*sch[41]+FF*sch[28];
  F1025=AA*sch[44]+DD*sch[31]+BB*sch[41]+EE*sch[24]+CC*sch[40]+FF*sch[27];
  F1115=AA*sch[46]+DD*sch[25]+BB*sch[47]+EE*sch[26]+CC*sch[39]+FF*sch[23];

  for(i=0;i<3;i++) {
    fLLprime[i]= I0013[i] * F0013 + I1003[i] * F1003 + I0103[i] * F0103
      + I0015[i] * F0015 + I1005[i] * F1005 + I0105[i] * F0105
      + I1025[i] * F1025 + I0125[i] * F0125 + I1115[i] * F1115
      + I2015[i] * F2015 + I0215[i] * F0215;

    fx3[i]= fLLprime[i]*factor;
  }

  for(i=0;i<3;i++) { quadmath_snprintf(cv3[i], 256,  "%*.30Qf", fx3[i]); } printf("fx3 = \n  %s %s %s  \n", cv3[0], cv3[1], cv3[2]);

  //
  // Fx6
  //
  // constants associated to the node shape function
  AA=-4.0q*one_o_dr*one_o_dr;
  BB=-4.0q*one_o_dr*one_o_ds;
  CC=4.0q*one_o_dr+8.0*r1*one_o_dr*one_o_dr+4.0q*s1*one_o_dr*one_o_ds;
  DD=4.0q*r1*one_o_dr*one_o_ds;
  EE=-4.0q*r1*one_o_dr-4.0q*r1*r1*one_o_dr*one_o_dr-4.0q*r1*s1*one_o_dr*one_o_ds;

  // triple integrals required at that node
  F0013=AA*sch[8]+BB*sch[12]+CC*sch[5]+DD*sch[6]+EE*sch[0];
  F0103=AA*sch[10]+BB*sch[11]+CC*sch[7]+DD*sch[4]+EE*sch[2];
  F1003=AA*sch[13]+BB*sch[10]+CC*sch[3]+DD*sch[7]+EE*sch[1];

  F0015=AA*sch[29]+BB*sch[23]+CC*sch[20]+DD*sch[21]+EE*sch[15];
  F0105=AA*sch[35]+BB*sch[36]+CC*sch[22]+DD*sch[19]+EE*sch[17];
  F1005=AA*sch[37]+BB*sch[35]+CC*sch[18]+DD*sch[22]+EE*sch[16];
  F0215=AA*sch[39]+BB*sch[47]+CC*sch[26]+DD*sch[34]+EE*sch[30];
  F2015=AA*sch[42]+BB*sch[46]+CC*sch[33]+DD*sch[25]+EE*sch[29];
  F0125=AA*sch[40]+BB*sch[41]+CC*sch[24]+DD*sch[32]+EE*sch[28];
  F1025=AA*sch[44]+BB*sch[40]+CC*sch[31]+DD*sch[24]+EE*sch[27];
  F1115=AA*sch[46]+BB*sch[39]+CC*sch[25]+DD*sch[26]+EE*sch[23];

  for(i=0;i<3;i++) {
    fLLprime[i]= I0013[i] * F0013 + I1003[i] * F1003 + I0103[i] * F0103
      + I0015[i] * F0015 + I1005[i] * F1005 + I0105[i] * F0105
      + I1025[i] * F1025 + I0125[i] * F0125 + I1115[i] * F1115
      + I2015[i] * F2015 + I0215[i] * F0215;

    fx6[i]= fLLprime[i]*factor;
  }
  //
  // Fx7
  //
  // constants associated to the node shape function
  AA=-4.0q*one_o_ds*one_o_ds;
  BB=-4.0q*one_o_dr*one_o_ds;
  CC=4.0q*s1*one_o_dr*one_o_ds;
  DD=4.0q*one_o_ds+8.0*s1*one_o_ds*one_o_ds+4.0q*r1*one_o_dr*one_o_ds;
  EE=-4.0q*s1*one_o_ds-4.0q*s1*s1*one_o_ds*one_o_ds-4.0q*r1*s1*one_o_dr*one_o_ds;

  // triple integrals required at that node
  F0013=AA*sch[9]+BB*sch[12]+CC*sch[5]+DD*sch[6]+EE*sch[0];
  F0103=AA*sch[14]+BB*sch[11]+CC*sch[7]+DD*sch[4]+EE*sch[2];
  F1003=AA*sch[11]+BB*sch[10]+CC*sch[3]+DD*sch[7]+EE*sch[1];

  F0015=AA*sch[30]+BB*sch[23]+CC*sch[20]+DD*sch[21]+EE*sch[15];
  F0105=AA*sch[38]+BB*sch[36]+CC*sch[22]+DD*sch[19]+EE*sch[17];
  F1005=AA*sch[36]+BB*sch[35]+CC*sch[18]+DD*sch[22]+EE*sch[16];
  F0215=AA*sch[43]+BB*sch[47]+CC*sch[26]+DD*sch[34]+EE*sch[30];
  F2015=AA*sch[39]+BB*sch[46]+CC*sch[33]+DD*sch[25]+EE*sch[29];
  F0125=AA*sch[45]+BB*sch[41]+CC*sch[24]+DD*sch[32]+EE*sch[28];
  F1025=AA*sch[41]+BB*sch[40]+CC*sch[31]+DD*sch[24]+EE*sch[27];
  F1115=AA*sch[47]+BB*sch[39]+CC*sch[25]+DD*sch[26]+EE*sch[23];

  for(i=0;i<3;i++) {
    fLLprime[i]= I0013[i] * F0013 + I1003[i] * F1003 + I0103[i] * F0103
      + I0015[i] * F0015 + I1005[i] * F1005 + I0105[i] * F0105
      + I1025[i] * F1025 + I0125[i] * F0125 + I1115[i] * F1115
      + I2015[i] * F2015 + I0215[i] * F0215;

    fx7[i]= fLLprime[i]*factor;
  }
  //
  // Fx8
  //
  // constants associated to the node shape function
  BB=4.0q;
  CC=-4.0q*s1;
  DD=-4.0q*r1;
  EE=4.0q*r1*s1;

  // triple integrals required at that node
  F0013=BB*sch[12]+CC*sch[5]+DD*sch[6]+EE*sch[0];
  F0103=BB*sch[11]+CC*sch[7]+DD*sch[4]+EE*sch[2];
  F1003=BB*sch[10]+CC*sch[3]+DD*sch[7]+EE*sch[1];

  F0015=BB*sch[23]+CC*sch[20]+DD*sch[21]+EE*sch[15];
  F0105=BB*sch[36]+CC*sch[22]+DD*sch[19]+EE*sch[17];
  F1005=BB*sch[35]+CC*sch[18]+DD*sch[22]+EE*sch[16];
  F0215=BB*sch[47]+CC*sch[26]+DD*sch[34]+EE*sch[30];
  F2015=BB*sch[46]+CC*sch[33]+DD*sch[25]+EE*sch[29];
  F0125=BB*sch[41]+CC*sch[24]+DD*sch[32]+EE*sch[28];
  F1025=BB*sch[40]+CC*sch[31]+DD*sch[24]+EE*sch[27];
  F1115=BB*sch[39]+CC*sch[25]+DD*sch[26]+EE*sch[23];

  for(i=0;i<3;i++) {
    fLLprime[i]= I0013[i] * F0013 + I1003[i] * F1003 + I0103[i] * F0103
      + I0015[i] * F0015 + I1005[i] * F1005 + I0105[i] * F0105
      + I1025[i] * F1025 + I0125[i] * F0125 + I1115[i] * F1115
      + I2015[i] * F2015 + I0215[i] * F0215;

    fx8[i]= fLLprime[i]*factor*one_o_dr*one_o_ds;
  }

  for(i=0;i<3;i++) {
    ftot[i]= fx6[i] +fx7[i] +fx3[i] +fx4[i] +fx5[i] +fx8[i];
  }

  *fp3x = fx3[0]; *fp3y = fx3[1]; *fp3z = fx3[2];
  *fp4x = fx4[0]; *fp4y = fx4[1]; *fp4z = fx4[2];
  *fp5x = fx5[0]; *fp5y = fx5[1]; *fp5z = fx5[2];
  *fp6x = fx6[0]; *fp6y = fx6[1]; *fp6z = fx6[2];
  *fp7x = fx7[0]; *fp7y = fx7[1]; *fp7z = fx7[2];
  *fp8x = fx8[0]; *fp8y = fx8[1]; *fp8z = fx8[2];

printf("End of SegQuadTrianForces Function \n");
}


/*
 *
 *  MAIN,
 *  Main function is simply used to read segment and surface data
 *  from file 'input.txt'
 *
 */
int main()
{
  __float128 x1[3], x2[3], x3[3], x4[3], x5[3];
  __float128 b[3], MU, NU, acore;
  __float128 fx3[3], fx4[3], fx5[3];
  __float128 fx6[3], fx7[3], fx8[3];
  char x11[256], x12[256], x13[256]; // String size
  char x21[256], x22[256], x23[256];
  char x31[256], x32[256], x33[256];
  char x41[256], x42[256], x43[256];
  char x51[256], x52[256], x53[256];
  char b1[256], b2[256], b3[256];
  char c_mu[256], c_nu[256], c_a[256];
  char cv3[3][256];
  int i;

  FILE * ptr_file;
  printf("Reading input.txt... \n");

  ptr_file =fopen("input.txt", "r");
  fscanf(ptr_file, "%s %s %s \n", x11, x12, x13 );
  fscanf(ptr_file, "%s %s %s \n", x21, x22, x23 );
  fscanf(ptr_file, "%s %s %s \n", x31, x32, x33 );
  fscanf(ptr_file, "%s %s %s \n", x41, x42, x43 );
  fscanf(ptr_file, "%s %s %s \n", x51, x52, x53 );
  fscanf(ptr_file, "%s %s %s \n", b1, b2, b3 );
  fscanf(ptr_file, "%s \n", c_mu );
  fscanf(ptr_file, "%s \n", c_nu );
  fscanf(ptr_file, "%s \n", c_a );

  // Converstion Char into quad prec floating point numbers
  x1[0] = strtoflt128(x11, NULL); x1[1] = strtoflt128(x12, NULL); x1[2] = strtoflt128(x13, NULL);
  x2[0] = strtoflt128(x21, NULL); x2[1] = strtoflt128(x22, NULL); x2[2] = strtoflt128(x23, NULL);
  x3[0] = strtoflt128(x31, NULL); x3[1] = strtoflt128(x32, NULL); x3[2] = strtoflt128(x33, NULL);
  x4[0] = strtoflt128(x41, NULL); x4[1] = strtoflt128(x42, NULL); x4[2] = strtoflt128(x43, NULL);
  x5[0] = strtoflt128(x51, NULL); x5[1] = strtoflt128(x52, NULL); x5[2] = strtoflt128(x53, NULL);
  b[0] = strtoflt128(b1, NULL); b[1] = strtoflt128(b2, NULL); b[2] = strtoflt128(b3, NULL);
  MU = strtoflt128(c_mu, NULL);
  NU = strtoflt128(c_nu, NULL);
  acore = strtoflt128(c_a, NULL);

  // Converstion Char into quad prec floating point numbers
  quadmath_snprintf(x11, 256, "%*.30Qf", x1[0]); quadmath_snprintf(x12, 256,  "%*.30Qf", x1[1]); quadmath_snprintf(x13, 256,  "%*.30Qf", x1[2]); // String size must correspond to the one used for the reading
  quadmath_snprintf(x21, 256,  "%*.30Qf", x2[0]); quadmath_snprintf(x22, 256,  "%*.30Qf", x2[1]); quadmath_snprintf(x23, 256,  "%*.30Qf", x2[2]);
  quadmath_snprintf(x31, 256,  "%*.30Qf", x3[0]); quadmath_snprintf(x32, 256,  "%*.30Qf", x3[1]); quadmath_snprintf(x33, 256,  "%*.30Qf", x3[2]);

  quadmath_snprintf(x41, 256,  "%*.30Qf", x4[0]); quadmath_snprintf(x42, 256,  "%*.30Qf", x4[1]); quadmath_snprintf(x43, 256,  "%*.30Qf", x4[2]);
  quadmath_snprintf(x51, 256,  "%*.30Qf", x5[0]); quadmath_snprintf(x52, 256,  "%*.30Qf", x5[1]); quadmath_snprintf(x53, 256,  "%*.30Qf", x5[2]);
  quadmath_snprintf(b1 , 256,  "%*.30Qf", b[0] ); quadmath_snprintf(b2, 256,  "%*.30Qf", b[1]); quadmath_snprintf(b3, 256,  "%*.30Qf", b[2]);

  quadmath_snprintf(c_mu, 256,  "%*.30Qf",MU); quadmath_snprintf(c_nu, 256,  "%*.30Qf", NU); quadmath_snprintf(c_a, 256,  "%*.30Qf", acore);

  printf("Done reading input.txt... \n");
  printf("x1 = %s %s %s \n", x11, x12, x13 );
  printf("x2 = %s %s %s \n", x21, x22, x23 );
  printf("x3 = %s %s %s \n", x31, x32, x33 );
  printf("x4 = %s %s %s \n", x41, x42, x43 );
  printf("x5 = %s %s %s \n", x51, x52, x53 );
  printf("b = %s %s %s \n", b1, b2, b3 );
  printf("MU = %s NU = %s acore = %s \n", c_mu, c_nu, c_a );

  /* Call Nodal Force Function */
  SegQuadTrianForces(x1[0], x1[1], x1[2],
                     x2[0], x2[1], x2[2],
                     x3[0], x3[1], x3[2],
                     x4[0], x4[1], x4[2],
                     x5[0], x5[1], x5[2],
                     b[0],  b[1],  b[2],
                     MU, NU, acore,
                     &fx3[0], &fx3[1], &fx3[2],
                     &fx4[0], &fx4[1], &fx4[2],
                     &fx5[0], &fx5[1], &fx5[2],
                     &fx6[0], &fx6[1], &fx6[2],
                     &fx7[0], &fx7[1], &fx7[2],
                     &fx8[0], &fx8[1], &fx8[2]);

  for(i=0;i<3;i++) { quadmath_snprintf(cv3[i], 256,  "%*.30Qf", fx3[i]); } printf("fx3 = \n  %s %s %s  \n", cv3[0], cv3[1], cv3[2]);
  for(i=0;i<3;i++) { quadmath_snprintf(cv3[i], 256,  "%*.30Qf", fx4[i]); } printf("fx4 = \n  %s %s %s  \n", cv3[0], cv3[1], cv3[2]);
  for(i=0;i<3;i++) { quadmath_snprintf(cv3[i], 256,  "%*.30Qf", fx5[i]); } printf("fx5 = \n  %s %s %s  \n", cv3[0], cv3[1], cv3[2]);
  for(i=0;i<3;i++) { quadmath_snprintf(cv3[i], 256,  "%*.30Qf", fx6[i]); } printf("fx6 = \n  %s %s %s  \n", cv3[0], cv3[1], cv3[2]);
  for(i=0;i<3;i++) { quadmath_snprintf(cv3[i], 256,  "%*.30Qf", fx7[i]); } printf("fx7 = \n  %s %s %s  \n", cv3[0], cv3[1], cv3[2]);
  for(i=0;i<3;i++) { quadmath_snprintf(cv3[i], 256,  "%*.30Qf", fx8[i]); } printf("fx8 = \n  %s %s %s  \n", cv3[0], cv3[1], cv3[2]);

  return 0;
}
