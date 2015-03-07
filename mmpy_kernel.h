#ifndef __kernel_H_
#define __kernel_H_
//
///opt/apps/cuda/5.5/bin/nvcc  -ftz=true -m64 -O3 -gencode arch=compute_35,code=\"sm_35,compute_35\" -Xcompiler -fopenmp -Xptxas -dlcm=ca -c -D_DOUBLE --ptxas-options=-v -DBLOCKDIM_X=1 -DBLOCKDIM_Y=1 -DBLOCKDIM_Z=1 mmpy_kernel.cu 

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

void interp_kernel_wrapper(REAL *z, int mm1, int mm2, int mm3, REAL *u, int n1,int n2,int n3 );
void comm3_kernel_wrapper(REAL* u,int n1,int n2,int n3, int I, int J, int K);
void resid_kernel_wrapper(REAL ***u, REAL*** v, REAL*** r, int n1,int n2,int n3);
void psinv_kernel_wrapper(REAL* r, REAL* u, int n1,int n2,int n3, double c[4]);
void rprj3_kernel_wrapper(REAL*** r, int m1k,int m2k,int m3k, REAL*** s,int m1j,int m2j,int m3j);
REAL* copy_3d(REAL *** mat, int n_1, int n_2, int n_3);

#ifdef __cplusplus
}
#endif

#endif