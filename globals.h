/*
---------------------------------------------------------------------
  Parameter lm (declared and set in "npbparams.h") is the log-base2 of 
  the edge size max for the partition on a given node, so must be changed 
  either to save space (if running a small case) or made bigger for larger 
  cases, for example, 512^3. Thus lm=7 means that the largest dimension 
  of a partition that can be solved on a node is 2^7 = 128. lm is set 
  automatically in npbparams.h
  Parameters ndim1, ndim2, ndim3 are the local problem dimensions. 
---------------------------------------------------------------------
*/
#include <stdbool.h>
#include "npbparams.h"
#include "includes.h"

/*
integer nm      ! actual dimension including ghost cells for communications
     >      , nv      ! size of rhs array
     >      , nr      ! size of residual array
     >      , nm2     ! size of communication buffer
     >      , maxlevel! maximum number of levels
*/
int nm;  //=n+2
int nv;  //=(n+2)^3
int nm2; //2D-array 2*(n+2)^2
int nr;  //8/7 * nv + 8/7*nm^2 + 8*5/7*nm + 8*lm 
int m;
inline void init_globals()
{
    nm = 2 + (1 << lm);
    nv = (2 + (1 << ndim1)) * (2+ (1<<ndim2)) * (2+(1<<ndim3));
    nm2 = 2*nm*nm;
    nr = 8 * (nv+nm*nm+5*nm+7*lm)/7; // = upper bound for sum_i (n_i+2)^3
    m = nm+1;
}
#define maxlevel 11

//---------------------------------------------------------------------
int nx[maxlevel];
int ny[maxlevel];
int nz[maxlevel];

char Class;// = CLASS;
int debug_vec[8];
int ir[maxlevel];
int m1[maxlevel];
int m2[maxlevel];
int m3[maxlevel];
int lt, lb;

//---------------------------------------------------------------------
//  Set at m=1024, can handle cases up to 1024^3 case
//---------------------------------------------------------------------
bool timeron = false;

//Timer IDs
const int T_total=0,  T_init=1,  T_bench=2,  T_mg3P=3,
      T_psinv=4,  T_resid=5, T_resid2=6, T_rprj3=7,
      T_interp=8, T_norm2=9, T_last=9;

