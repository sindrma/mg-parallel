#ifndef _UTILITY_H_
#define _UTILITY_H_

#include "types.h"
#include "includes.h"

#define ERROR(s) printf s; printf("\n"); fflush(stdout);
#define DEBUG(s) printf s; printf("\n"); fflush(stdout);

typedef struct grid_s
{
    int is1, is2, is3, ie1, ie2, ie3;
} grid_t;

// generate mm number of random global indeces
void zran3(REAL ***z,int n1,int n2,int n3,int nx,int ny,int* j1,int* j2,int* j3,int *m1, int *m0, int mm, grid_t* grid);
// generate the v vector of A*u = v
void gen_v(REAL ***z,int n1,int n2,int n3,int nx,int ny, grid_t* grid);
// The original zran3 function.
void gen_v_orig(REAL ***z,int n1,int n2,int n3,int nx,int ny, grid_t* grid);

//bubble does a bubble sort in direction dir
void bubble(double ten[],int j1[],int j2[],int j3[],int m,int ind );

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
