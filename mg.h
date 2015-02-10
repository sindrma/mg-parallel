#ifndef __MG__H_
#define __MG__H_

#include "types.h"
#include "utility.h"

void setup(int *n1, int *n2, int *n3, grid_t* grid);
int getInputPars();
void resid(REAL ***u, REAL*** v, REAL*** r, int n1,int n2,int n3, double a[4]);
void mg3P(REAL ****u, REAL*** v, REAL**** r, double a[4], double c[4], int n1,int n2,int n3);
double norm2u3(REAL*** r,int n1,int n2,int n3, int nx,int ny,int nz);

void interp(REAL*** z,int mm1,int mm2,int mm3, REAL*** u, int n1,int n2,int n3 );
void psinv(REAL*** r, REAL ***u,int n1,int n2,int n3, double c[4]);
void rprj3(REAL*** r, int m1k,int m2k,int m3k, REAL*** s,int m1j,int m2j,int m3j);

void comm3(REAL ***u,int n1,int n2,int n3);

#endif
