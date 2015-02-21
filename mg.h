#ifndef __MG__H_
#define __MG__H_

#include "types.h"
#include "includes.h"


// generate mm number of random global indeces
void zran3(REAL ***z,int n1,int n2,int n3,int nx,int ny,int* j1,int* j2,int* j3,int *m1, int *m0, int mm, grid_t* grid);
// generate the v vector of A*u = v
void gen_v(REAL ***z,int n1,int n2,int n3,int nx,int ny, grid_t* grid);
// The original zran3 function.
void gen_v_orig(REAL ***z,int n1,int n2,int n3,int nx,int ny, grid_t* grid);

int getInputPars();
void resid(REAL ***u, REAL*** v, REAL*** r, int n1,int n2,int n3, double a[4]);
void mg3P(REAL ****u, REAL*** v, REAL**** r, double a[4], double c[4], int n1,int n2,int n3, int restriction);
double norm2u3(REAL*** r,int n1,int n2,int n3, int nx,int ny,int nz);

void exchange(REAL ***r, int n1,int n2,int n3 );
void interp_mpi(REAL ***z, int mm1, int mm2, int mm3, REAL ***u, int n1,int n2,int n3 );
void interp(REAL*** z,int mm1,int mm2,int mm3, REAL*** u, int n1,int n2,int n3 );
void psinv(REAL*** r, REAL ***u,int n1,int n2,int n3, double c[4]);
void rprj3_mpi(REAL*** r, int m1k,int m2k,int m3k, REAL*** s,int m1j,int m2j,int m3j);
void rprj3(REAL*** r, int m1k,int m2k,int m3k, REAL*** s,int m1j,int m2j,int m3j);

void comm3(REAL ***u,int n1,int n2,int n3);

#endif
