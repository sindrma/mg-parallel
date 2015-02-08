#ifndef __MG__H_
#define __MG__H_

#include "types.h"

typedef struct grid_s
{
    int is1, is2, is3, ie1, ie2, ie3;
} grid_t;


void zran3(REAL ***z,int n1,int n2,int n3,int nx,int ny,int* j1,int* j2,int* j3,int *m1, int *m0, int mm, grid_t* grid);
void gen_v(REAL ***z,int n1,int n2,int n3,int nx,int ny, grid_t* grid);
void gen_v_orig(REAL ***z,int n1,int n2,int n3,int nx,int ny, grid_t* grid);
void setup(int *n1, int *n2, int *n3, grid_t* grid);
int getInputPars();
void resid(REAL ***u, REAL*** v, REAL*** r, int n1,int n2,int n3, double a[4]);
void mg3P(REAL ****u, REAL*** v, REAL**** r, double a[4], double c[4], int n1,int n2,int n3);
double norm2u3(REAL*** r,int n1,int n2,int n3, int nx,int ny,int nz);

void interp(REAL*** z,int mm1,int mm2,int mm3, REAL*** u, int n1,int n2,int n3 );
void psinv(REAL*** r, REAL ***u,int n1,int n2,int n3, double c[4]);
void rprj3(REAL*** r, int m1k,int m2k,int m3k, REAL*** s,int m1j,int m2j,int m3j);
void bubble(double ten[],int j1[],int j2[],int j3[],int m,int ind );

void comm3(REAL ***u,int n1,int n2,int n3);

#endif
