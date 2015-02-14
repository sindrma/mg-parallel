#ifndef _UTILITY_H_
#define _UTILITY_H_

#include "types.h"
#include "includes.h"

#define ERROR(s) printf s; printf("\n"); fflush(stdout);
#define DEBUG(s) printf s; printf("\n"); fflush(stdout);

//bubble does a bubble sort in direction dir
void bubble(double ten[],int j1[],int j2[],int j3[],int m,int ind );

REAL ** exchange_data(REAL** data,int size);
REAL ** getGhostCells(REAL*** mat,int x,int y,int z);
REAL *** splitMatrix(REAL*** mat,int x,int y,int z, int processorID);
REAL * flattenMatrix(REAL*** mat,int x,int y,int z);

double TestNorm(double r[],int n1,int n2,int n3);

//Fill the 3-dimensional array z with zeroes
void zero3old(REAL z[],int off,int n1,int n2,int n3);
void zero3(REAL ***z,int n1,int n2,int n3);

REAL ***alloc3D(int n, int m,int k);
void free3D(REAL*** arr);

//Allocate course and fine grids
//Finest grid has index depth-1 and has dimensions (n+2)*(m+2)*(k+2)
REAL **** allocGrids(size_t depth, size_t n1, size_t n2, size_t n3, size_t pad);
void      freeGrids(REAL**** grids);

#endif
