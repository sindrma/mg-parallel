/*
	Setup.c: contains any functions used in setting up the multigrid program
	-Working on:
*/
#include "setup.h"
#include "../npbparams.h"

//Timer IDs-
const int T_total=0,  T_init=1,  T_bench=2,  T_mg3P=3,
      T_psinv=4,  T_resid=5, T_resid2=6, T_rprj3=7,
      T_interp=8, T_norm2=9, T_last=9;

//original setup function- allocates the grid variables
void setup(int *n1, int *n2, int *n3, grid_t* grid)
{
	int j, k;
	
	int ax;
	int size1=3, size2=10;
	int *mi = malloc(sizeof(int)*size1*size2);
	int *ng = malloc(sizeof(int)*size1*size2);
	
	//old init globals function
	nm = 2 + (1 << lm);
	nv = (2 + (1 << ndim1)) * (2+ (1<<ndim2)) * (2+(1<<ndim3));
	nm2 = 2*nm*nm;
	nr = 8 * (nv+nm*nm+5*nm+7*lm)/7; // = upper bound for sum_i (n_i+2)^3
	m = nm+1;
	
	ng[  (lt-1)*size1]=nx[lt-1];
	ng[1+(lt-1)*size1]=ny[lt-1];
	ng[2+(lt-1)*size1]=nz[lt-1]  / global_params->mpi_size;

	for(ax=0;ax<size1;ax++)
		for(k=lt-2;k>=0;k--)
			ng[ax+k*size1]=ng[ax+(k+1)*size1]/2;
	
	for(k=lt-2;k>=0;k--) {
		nx[k]=ng[  k*size1];
		ny[k]=ng[1+k*size1];
		nz[k]=ng[2+k*size1];
	}
	
	for(k=lt-1;k>=0;k--) {
		for(ax=0;ax<size1;ax++) {
			mi[ax+k*size1] = 2 + ng[ax+k*size1];
		}
		m1[k]=mi[  k*size1];
		m2[k]=mi[1+k*size1];
		m3[k]=(mi[2+k*size1] -2 ) / global_params->mpi_size + 2;
	}
	
	k = lt-1;
	grid->is1 = 2 + ng[k*size1] - ng[k*size1];
	grid->ie1 = 1 + ng[k*size1];
	*n1= 3 + grid->ie1 - grid->is1;
	grid->is2 = 2 + ng[1+k*size1] - ng[1+k*size1];
	grid->ie2 = 1 + ng[1+k*size1]; 
	*n2= 3 + grid->ie2 - grid->is2;
	grid->is3 = 2 + ng[2+k*size1] - ng[2+k*size1];
	grid->ie3 = 1 + ng[2+k*size1];
	//NOTE: we are splitting the data into strips
	//global_params->mpi_size
	// printf ("IE3: %d, IS3: %d \n",grid->ie3, grid->is3);
	*n3= (3 + grid->ie3 - grid->is3);
	// printf ("N3: %d \n",*n3);
	
	ir[lt-1]=0;
	for(j = lt-2;j>=0;j--) {
		ir[j]=ir[j+1]+m1[j+1]*m2[j+1]*m3[j+1];
	}

	free(mi);
	free(ng);
}

struct params* setup_local(int argc, const char **argv)
{
	struct params parameters;
	struct params *p = (struct params *) malloc(sizeof(struct params));
	int nArg;
	
	//setup default struct
	p->n_it		= 4;
	p->n_size	= 256;
	p->lt 		= 8;
	p->seed		= 314159265.0;
	p->class	= 'U';
	MPI_Comm_rank( MPI_COMM_WORLD, &p->mpi_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &p->mpi_size );
	
	//check command line input for manual entry
	for (nArg=1; nArg < argc; nArg+=2){
		// printf("%d : %s\n",nArg,argv[nArg]);
		if (strcmp(argv[nArg],"-nit") == 0) {
			if (sscanf (argv[nArg + 1], "%i", &p->n_it)!=1) {
				printf ("error - not an integer");
			}
		}
		if (strcmp(argv[nArg],"-n") == 0) {
			if (sscanf (argv[nArg + 1], "%i", &p->n_size)!=1) {
				printf ("error - not an integer");
			}
		}
		
		if (strcmp(argv[nArg],"-lt") == 0) {
			if (sscanf (argv[nArg + 1], "%i", &p->lt)!=1) {
				printf ("error - not an integer");
			} else {
				//TODO : add lt>maxlevel check (it's in the original code)
				
			}
		}
		if (strcmp(argv[nArg],"-s") == 0) {
			if (sscanf (argv[nArg + 1], "%f", &p->seed)!=1) {
				printf ("error - not a float");
			}
		}
		
	}
	
	//determine class of test
	if( p->n_size==32 && p->n_it==4 )
		p->class = 'S';
	else if( p->n_size==64 && p->n_it==40 )
		p->class = 'W';
	else if( p->n_size==256 && p->n_it==20 )
		p->class= 'B';
	else if( p->n_size==512 && p->n_it==20 )
		p->class = 'C';
	else if( p->n_size==256 && p->n_it==4 )
		p->class = 'A';
	
	//print function-
	if(p->mpi_rank == 0){
		printf("\n NAS Parallel Benchmarks C version\n");
		printf(" Multithreaded Version %s.%c np=%d\n", "MG", p->class, omp_get_max_threads());
		printf(" Size:  %dx%dx%d\n Iterations:   %d\n", p->n_size, p->n_size, p->n_size, p->n_it );
	}
	
	return p;
}

//sets a and c matrices
//i decided to make these functions in the case that we decide to change it- using set coefficients based on class seems weird to me
void set_a(double *a ,struct params* global){
	a[0] = -8.0/3.0; 
	a[1] =  0.0;
	a[2] =  1.0/6.0; 
	a[3] =  1.0/12.0;
}
void set_c(double *c ,struct params* global){
	if(global->class=='A' || global->class=='S' || global->class=='W') {
		//---------------------------------------------------------------------
		//     Coefficients for the S(a) smoother
		//---------------------------------------------------------------------
		c[0] =  -3.0/8.0;
		c[1] =  +1.0/32.0;
		c[2] =  -1.0/64.0;
		c[3] =   0.0;
	} else {
		//---------------------------------------------------------------------
		//     Coefficients for the S(b) smoother
		//---------------------------------------------------------------------
		c[0] =  -3.0/17.0;
		c[1] =  +1.0/33.0;
		c[2] =  -1.0/61.0;
		c[3] =   0.0;
	}
}
