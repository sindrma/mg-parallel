#ifndef _UTILITY_H_
#define _UTILITY_H_

#include "types.h"
#include "includes.h"

#define ERROR(s) printf s; printf("\n"); fflush(stdout);
#define DEBUG(s) printf s; printf("\n"); fflush(stdout);

//bubble does a bubble sort in direction dir
void bubble(double ten[],int j1[],int j2[],int j3[],int m,int ind );

void printMatrix(REAL *** mat,int x,int y,int z);
REAL *** merge_matrices(REAL*** mat1, REAL *** mat2,int x1, int y1, int z1, int x2, int y2, int z2, bool buffer);
REAL ** exchange_data(REAL** data,int size);
REAL ** getGhostCells(REAL*** mat,int x,int y,int z);
REAL *** splitMatrix(REAL*** mat,int x,int y,int z, int processorID, int size,bool addBuffer);
REAL *** unflattenMatrix(REAL* mat,int x,int y,int z);
REAL * flattenMatrix(REAL*** mat,int x,int y,int z);

double TestNorm(double r[],int n1,int n2,int n3);

REAL **alloc2D(int n, int m);
void free2D(REAL** buffer, int n);

//Fill the 3-dimensional array z with zeroes
void zero3old(REAL z[],int off,int n1,int n2,int n3);
void zero3(REAL ***z,int n1,int n2,int n3);

REAL ***alloc3D(int n, int m,int k);
// does not free up all pointers 
void free3D(REAL*** arr);
void free3D_old(REAL*** arr, int n, int m);

//Allocate course and fine grids
//Finest grid has index depth-1 and has dimensions (n+2)*(m+2)*(k+2)
REAL **** allocGrids(size_t depth, size_t n1, size_t n2, size_t n3, size_t pad);
void      freeGrids(REAL**** grids);

#endif
