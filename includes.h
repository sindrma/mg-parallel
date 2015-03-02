//C standard library include files
#include "ptools_ppf.h"
//Need this defined in order to get the function definitions of some non-standard functions (e.g. strdup)
#define _XOPEN_SOURCE 1000

#ifndef __INCLUDES_H__
#define __INCLUDES_H__
	#include <stdio.h>
	#include <stdlib.h>
	#include <math.h>
	#include <string.h>
	#include <strings.h>
	#include <stdbool.h>
	#include <ctype.h>
	#include <sys/time.h>
#endif

//NOTE: actual values are in setup.c- setting the values in a header file will cause duplicate definitions!!!
#ifndef __STRUCTS_AND_VARIABLES__
#define __STRUCTS_AND_VARIABLES__
	const int nx_default; 
	const int ny_default; 
	const int nz_default; 
	const int nit_default; 
	const int lm; 
	const int lt_default; 
	const int debug_default; 
	const int ndim1; 
	const int ndim2; 
	const int ndim3; 
	const bool convertdouble;
	const char* compiletime;
	const char* npbversion;
	const char* cs1;
	const char* cs2;
	const char* cs3;
	const char* cs4;
	const char* cs5;
	const char* cs6;
	const char* cs7;
	int nm;  //=n+2
	int nv;  //=(n+2)^3
	int nm2; //2D-array 2*(n+2)^2
	int nr;  //8/7 * nv + 8/7*nm^2 + 8*5/7*nm + 8*lm 
	int m;
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
    int *split_levels;

	//---------------------------------------------------------------------
	//  Set at m=1024, can handle cases up to 1024^3 case
	//---------------------------------------------------------------------

	//Timer IDs
	const int T_total,  T_init,  T_bench,  T_mg3P,
	      T_psinv,  T_resid, T_resid2, T_rprj3,
	      T_interp, T_norm2, T_last;

	//Used to store variables that are the same accross processors
	struct params{
		int n_size;	//size of the cube
		int n_it;	//number of iterations
		int lt;		//threads requested? not sure if this is correct
		int mpi_size;	//# of processors available
		int mpi_rank;	//processor id
        int mpi_orig_size;
        bool active;
		double seed;	//seed used to generate random indices
		char class;	//class of test- legacy code
		int geometry[3];//geometry of how the data is split - *UNIMPLEMENTED*
	};
	
	struct params *global_params;
	
	//Represents the benchmark results
	typedef struct results_s
	{
	    char* name;
	    char* MachineName;
	    char* PrLang;
	    char  clss;
	    int   n1, n2, n3, niter;
	    double time,acctime,wctime,mops;
	    double tmSent, tmReceived;
	    int   RecArrSize;
	    char* optype;
	    int   numthreads;
	    int   verified;
	} results_t;
	
	typedef struct grid_s
	{
	    int is1, is2, is3, ie1, ie2, ie3;
	} grid_t;
#endif