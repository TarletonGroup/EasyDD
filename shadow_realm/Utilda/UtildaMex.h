/* global function */
void UtildaMex(double *x0, double *y0, double *z0, double *bx, double *by, double *bz,
               double *x1, double *y1, double *z1, double *x2, double *y2, double *z2,
               double *nx, double *ny, double *nz,
               double NU, int point_array_length, int segments_array_length,
               double *Ux, double *Uy, double *Uz);

/* auxiliary functions to global */
void SlipPlaneCheck(double *nvec, double normn, double nx, double ny, double nz);
void PlaneLineIntersection(double *n, double *V0, double const *P0, double const *P1, double *I);

/* barnett triangle function*/
void BarnettTriangle(double *point, double *A, double *B, double *C, double *b, double NU, double *utilda);

/* auxiliary functions to barnett triangle function */
void fab(double *b, double *t, double *lamA, double *lamB, double RA, double RB, double *f_vec);
void gab(double *b, double *lamA, double *lamB, double *g_vec);
double solang(double *lamA, double *lamB, double *lamC, double *p, double *plane_n);

/*random functions*/
void safenorm( double R[3],double unitR[3]);