/*
	-=Current Goals=-
	-start parallelizing the "v" cycle

	-=Command Line=-
	mpirun -np 2 ./mg -n 64
	parameter options:
	-s	: seed
	-n	: problem size
	-nit	: number of iterations
	-lt	: sets the lt... i think this defines how large the "W" cycle is?
*/

#include <omp.h>
#include "functions/setup.h"
#include "functions/results.h"
#include "random.h"
#include "utility.h"
#include "timer.h"
#include "mg.h"

//Some global constants
const char * BMName="MG"; //Benchmark name

int main(int argc, const char **argv)
{	
	int argc_test,mpi_rank,mpi_size;
	char ** argv_test;
	MPI_Init( &argc_test, &argv_test );
	/*local processor setup- 
		-initializes variables that will be the same for every processor
		-parses command line options
	*/
	global_params = setup_local(argc,argv);
	
	//k is the current level. It is passed down through subroutine args and is NOT global. 
	//it is the current iteration.
	int k, it;
	//pointers of pointers in 3-D
	REAL**** u, //approximation matrix
		**** r, //residual error matrix
		***  v, //values matrix
		a[4], c[4]; //what are these used for?
	
	double rnm2, tinit;
	int n1, n2, n3, nit;
	int i;
	
	//use global parameters to set these
	lt = global_params->lt;
	nit = global_params->n_it;
	nx[lt-1] = global_params->n_size;
	ny[lt-1] = global_params->n_size;
	nz[lt-1] = global_params->n_size;
	Class = global_params->class;
	
	//set a and b matrices for smoother using the class value
	set_a(a,global_params);
	set_c(c,global_params);
	
	k = lt;
	
	//setup initializes some of the variables used...
	grid_t grid;
	setup(&n1, &n2, &n3, &grid);
	
	//processor 0 tracks time and sets up solution matrix
	init_timers();
	timer_start(T_init);
	
	v = alloc3D(n1, n2, n3);
	gen_v(v,n1,n2,n3,nx[lt-1],ny[lt-1], &grid);
	
	//Initialize arrays- currently does this in strips
	u = allocGrids(lt, n1-2, n2-2, (n3-2) / global_params->mpi_size, 2);
	r = allocGrids(lt, n1-2, n2-2, (n3-2) / global_params->mpi_size, 2);
	zero3(u[0],n1,n2,n3 / global_params->mpi_size); //zero-out all of u
	
	//create local v (for residual)
	REAL ***  local_v = splitMatrix(
		v, //matrix to split
		n1,n2,n3, //matrix size
		global_params->mpi_rank,
		global_params->mpi_size,
		true //result matrix should have a 1 cell boundary buffer
	);
	
	timer_stop(T_init);
	timer_start(T_bench);
	resid(u[0],local_v,r[0],n1,n2,(n3-2)/global_params->mpi_size + 2,a);
	
	//each processor runs the multigrid algorithm
	for(it=1;it<=nit;it++) {
		//actual call to multigrid
		mg3P(u,local_v,r,a,c,n1,n2,(n3-2)/global_params->mpi_size + 2,0);
		//compute the residual error here...
		//only pass in the spliced portion of v...
		resid(u[0],local_v,r[0],n1,n2,(n3-2)/global_params->mpi_size + 2,a);
	} 
	
	//processor 1 interprets results
	if(global_params->mpi_rank == 0){
		timer_stop(T_bench);
		tinit = timer_elapsed(T_init);
		// printf(" Initialization time: %f seconds\n", tinit);
		
		// //combine r results from each processor
		// REAL * message = (REAL*) malloc(sizeof(REAL*)*(n1)*(n2)*(n3));
		// MPI_Recv(message,n1*n2*n3,MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// 
		// //unflatten message into 3d matrix
		// REAL *** strip = unflattenMatrix(message,n1,n2,(n3));
		// 
		// //combine results
		// REAL *** result = merge_matrices(
		// 	r[0],strip,
		// 	n1,n2,(n3),
		// 	n1,n2,(n3),
		// 	true	//1 cell of boundary data in each matrix
		// );
		
		// rnm2=norm2u3(result,n1,n2,n3,nx[lt-1],ny[lt-1],nz[lt-1]);
		rnm2=norm2u3(r[0],n1,n2,n3,nx[lt-1],ny[lt-1],nz[lt-1]);
		double tm = timer_elapsed(T_bench);

		//validates the results and prints to console
		interpret_results(rnm2, global_params, tm);
		
	} else {
	// 	//send residual result
	// 	REAL * message = flattenMatrix(r[0],n1,n2,(n3));
	// 	MPI_Send(message, (n1)*(n2)*(n3), MPI_DOUBLE, 0, global_params->mpi_rank, MPI_COMM_WORLD);
	}
	free(v);
	freeGrids(u);
	freeGrids(r);
	
	MPI_Barrier(MPI_COMM_WORLD);
	PPF_Print( MPI_COMM_WORLD, "Message from %N - finished\n" );
	
	MPI_Finalize();
	return 0;
}

void gen_v(REAL*** z,int n1,int n2,int n3,int nx,int ny, grid_t* grid)
{
	//---------------------------------------------------------------------
	//      Generate righthandside of the equation A*u = v
	//      for the serial version of the code.
	//---------------------------------------------------------------------
	int m0, m1, mm=10, i1, i2, i3, i;
	int *j1 = malloc(sizeof(int)*mm*2), 
		*j2 = malloc(sizeof(int)*mm*2),
		*j3 = malloc(sizeof(int)*mm*2);
	
	zran3(z,n1,n2,n3,nx,ny,j1,j2,j3, &m1, &m0, mm, grid);
	#pragma omp parallel for private(i1,i2,i3)
	for(i3=0;i3<n3;i3++)
		for(i2=0;i2<n2;i2++)
			for(i1=0;i1<n1;i1++)
				z[i3][i2][i1] = 0.0;
	for(i=mm;i>=m0;i--)
		z[j3[i-1]][j2[i-1]][j1[i-1]] = -1.0;
	for(i=mm;i>=m1;i--)
		z[j3[i-1+mm]][j2[i-1+mm]][j1[i-1+mm]] = 1.0;

	free(j1);
	free(j2);
	free(j3);

	comm3(z,n1,n2,n3);
}

void zran3(REAL ***z,int n1,int n2,int n3,int nx,int ny,int* j1,int* j2,int* j3,int *m1, int *m0, int mm, grid_t* grid){
	int is1 = grid->is1, is2 = grid->is2, is3 = grid->is3, ie1 = grid->ie1, ie2 = grid->ie2, ie3 = grid->ie3;
	int i, i0, i1, i2, i3, d1, e1, e2, e3;
	int *jg = malloc(sizeof(int)*4*mm*2);
	double xx, x0, x1, a1, a2, ai;
	double best;
	double *ten= malloc(sizeof(double)*mm*2);
	
	zero3(z,n1,n2,n3);
	i = is1-2+nx*(is2-2+ny*(is3-2));
	
	d1 = ie1 - is1 + 1;
	e1 = ie1 - is1 + 2;
	e2 = ie2 - is2 + 2;
	e3 = ie3 - is3 + 2;
	
	double seed=314159265.0, a=pow(5.0,13);
	//double rng = drand48();
	a1 = rnd_power( a, nx );
	a2 = rnd_power( a, nx*ny );
	ai = rnd_power( a, i );
	x0 = rnd_randlc( seed, ai );
	
	for(i3=2;i3<=e3;i3++) {
		x1 = x0;
		for(i2 = 2;i2<=e2;i2++) {
			xx = x1;
			rnd_vranlc( d1, xx, a,z[0][0],(1+n1*(i2-1+n2*(i3-1))));
			x1 = rnd_randlc( x1, a1 );
		}
		x0 = rnd_randlc( x0, a2 );
	}
	
	for(i=0;i<mm;i++) {
		ten[i+mm] = 0.0;
		j1[i+mm] = 0;
		j2[i+mm] = 0;
		j3[i+mm] = 0;
		ten[i] = 1.0;
		j1[i] = 0;
		j2[i] = 0;
		j3[i] = 0;
	}
	
	for(i3=1;i3<n3-1;i3++) {
		for(i2=1;i2<n2-1;i2++) {
			for(i1=1;i1<n1-1;i1++) {
				if( z[i3][i2][i1] > ten[mm] ) {
					ten[mm] = z[i3][i2][i1];
					j1[mm] = i1;
					j2[mm] = i2;
					j3[mm] = i3;
					bubble( ten, j1, j2, j3, mm, 1 );
				}
				if( z[i3][i2][i1] < ten[0] ) {
					ten[0] = z[i3][i2][i1];
					j1[0] = i1;
					j2[0] = i2;
					j3[0] = i3;
					bubble( ten, j1, j2, j3, mm, 0 );
				}
			}
		}
	}
	
	//---------------------------------------------------------------------
	//     Now which of these are globally best?
	//---------------------------------------------------------------------
	i1 = mm;
	i0 = mm;
	for(i=mm-1;i>=0;i--)
	{
		//best = z[0][0][j1[i1-1+mm]+n1*(j2[i1-1+mm]+n2*(j3[i1-1+mm]))];
		best = z[j1[i1-1+mm]][j2[i1-1+mm]][j3[i1-1+mm]];
		if(best==z[j1[i1-1+mm]][j2[i1-1+mm]][j3[i1-1+mm]])
		{
			jg[4*(i+mm)] = 0;
			jg[1+4*(i+mm)] = is1 - 2 + j1[i1-1+mm];
			jg[2+4*(i+mm)] = is2 - 2 + j2[i1-1+mm];
			jg[3+4*(i+mm)] = is3 - 2 + j3[i1-1+mm];
			i1 = i1-1;
		} else {
			jg[4*(i+mm)] = 0;
			jg[1+4*(i+mm)] = 0; 
			jg[2+4*(i+mm)] = 0; 
			jg[3+4*(i+mm)] = 0; 
		}         
		ten[i+mm] = best;
		
		best = z[j3[i0-1]][j2[i0-1]][j1[i0-1]];
		if(best==z[j3[i0-1]][j2[i0-1]][j1[i0-1]]) {
			jg[4*i] = 0;
			jg[1+4*i] = is1 - 2 + j1[i0-1];
			jg[2+4*i] = is2 - 2 + j2[i0-1];
			jg[3+4*i] = is3 - 2 + j3[i0-1];
			i0 = i0-1;
		} else {
			jg[4*i] = 0;
			jg[1+4*i] = 0; 
			jg[2+4*i] = 0; 
			jg[3+4*i] = 0; 
		}
		ten[i] = best;
	}
	
	free(jg);
	free(ten);
	
	*m1 = i1+1;
	*m0 = i0+1;
}

void gen_v_orig(REAL ***z,int n1,int n2,int n3,int nx,int ny, grid_t* grid){
	int is1 = grid->is1, is2 = grid->is2, is3 = grid->is3, ie1 = grid->ie1, ie2 = grid->ie2, ie3 = grid->ie3;
	int i0, m0, m1;
	
	int mm=10, i1, i2, i3, d1, e1, e2, e3;
	double xx, x0, x1, a1, a2, ai;
	double best;
	double *ten= malloc(sizeof(double)*mm*2);
	int i;
	int *j1 = malloc(sizeof(int)*mm*2),
		*j2 = malloc(sizeof(int)*mm*2),
		*j3 = malloc(sizeof(int)*mm*2);
	int *jg = malloc(sizeof(int)*4*mm*2);
	
	zero3(z,n1,n2,n3);
	i = is1-2+nx*(is2-2+ny*(is3-2));
	
	d1 = ie1 - is1 + 1;
	e1 = ie1 - is1 + 2;
	e2 = ie2 - is2 + 2;
	e3 = ie3 - is3 + 2;
	
	double seed=314159265.0, a=pow(5.0,13);
	//double rng = drand48();
	a1 = rnd_power( a, nx );
	a2 = rnd_power( a, nx*ny );
	ai = rnd_power( a, i );
	x0 = rnd_randlc( seed, ai );
	
	for(i3=2;i3<=e3;i3++)
	{
		x1 = x0;
		for(i2 = 2;i2<=e2;i2++)
		{
			xx = x1;
			rnd_vranlc( d1, xx, a,z[0][0],(1+n1*(i2-1+n2*(i3-1))));
			x1 = rnd_randlc( x1, a1 );
		}
		x0 = rnd_randlc( x0, a2 );
	}
	
	for(i=0;i<mm;i++) {
		ten[i+mm] = 0.0;
		j1[i+mm] = 0;
		j2[i+mm] = 0;
		j3[i+mm] = 0;
		ten[i] = 1.0;
		j1[i] = 0;
		j2[i] = 0;
		j3[i] = 0;
	}
	
	for(i3=1;i3<n3-1;i3++) {
		for(i2=1;i2<n2-1;i2++) {
			for(i1=1;i1<n1-1;i1++) {
				if( z[i3][i2][i1] > ten[mm] ) {
					ten[mm] = z[i3][i2][i1]; 
					j1[mm] = i1;
					j2[mm] = i2;
					j3[mm] = i3;
					bubble( ten, j1, j2, j3, mm, 1 );
				}
				if( z[i3][i2][i1] < ten[0] ) {
					ten[0] = z[i3][i2][i1]; 
					j1[0] = i1;
					j2[0] = i2;
					j3[0] = i3;
					bubble( ten, j1, j2, j3, mm, 0 );
				}
			}
		}
	}
	
	//---------------------------------------------------------------------
	//     Now which of these are globally best?
	//---------------------------------------------------------------------
	i1 = mm;
	i0 = mm;
	for(i=mm-1;i>=0;i--)
	{
		//best = z[0][0][j1[i1-1+mm]+n1*(j2[i1-1+mm]+n2*(j3[i1-1+mm]))];
		best = z[j1[i1-1+mm]][j2[i1-1+mm]][j3[i1-1+mm]];
		if(best==z[j1[i1-1+mm]][j2[i1-1+mm]][j3[i1-1+mm]]) {
			jg[4*(i+mm)] = 0;
			jg[1+4*(i+mm)] = is1 - 2 + j1[i1-1+mm]; 
			jg[2+4*(i+mm)] = is2 - 2 + j2[i1-1+mm]; 
			jg[3+4*(i+mm)] = is3 - 2 + j3[i1-1+mm]; 
			i1 = i1-1;
		} else {
			jg[4*(i+mm)] = 0;
			jg[1+4*(i+mm)] = 0; 
			jg[2+4*(i+mm)] = 0; 
			jg[3+4*(i+mm)] = 0; 
		}         
		ten[i+mm] = best;
		
		best = z[j3[i0-1]][j2[i0-1]][j1[i0-1]];
		if(best==z[j3[i0-1]][j2[i0-1]][j1[i0-1]]) {
			jg[4*i] = 0;
			jg[1+4*i] = is1 - 2 + j1[i0-1]; 
			jg[2+4*i] = is2 - 2 + j2[i0-1]; 
			jg[3+4*i] = is3 - 2 + j3[i0-1]; 
			i0 = i0-1;
		} else {
			jg[4*i] = 0;
			jg[1+4*i] = 0; 
			jg[2+4*i] = 0; 
			jg[3+4*i] = 0; 
		}
		ten[i] = best;
	}
	
	
	free(jg);
	free(ten);
	
	m1 = i1+1;
	m0 = i0+1;
	#pragma omp parallel for private(i1,i2,i3)
	for(i3=0;i3<n3;i3++)
		for(i2=0;i2<n2;i2++)
			for(i1=0;i1<n1;i1++)
				z[i3][i2][i1] = 0.0;
	for(i=mm;i>=m0;i--)
		z[j3[i-1]][j2[i-1]][j1[i-1]] = -1.0;
	for(i=mm;i>=m1;i--)
		z[j3[i-1+mm]][j2[i-1+mm]][j1[i-1+mm]] = 1.0;
	free(j1);
	free(j2);
	free(j3);
	
	comm3(z,n1,n2,n3);
}

double norm2u3(REAL*** r,int n1,int n2,int n3, int nx,int ny,int nz)
{
	//---------------------------------------------------------------------
	//     norm2u3 evaluates approximations to the L2 norm and the
	//     uniform (or L-infinity or Chebyshev) norm, under the
	//     assumption that the boundaries are periodic or zero.  Add the
	//     boundaries in with half weight (quarter weight on the edges
	//     and eighth weight at the corners) for inhomogeneous boundaries.
	//---------------------------------------------------------------------
	//      double precision r(n1,n2,n3)
	double rnmu = 0.0;
	double rnm2 = 0.0;
	int i1,i2,i3;
	double localmax;
	
	localmax = 0;
	for(i3=1;i3<n3-1;i3++)
	{
		double localmax = 0;
		for(i2=1;i2<n2-1;i2++)
		{
			for(i1=1;i1<n1-1;i1++)
			{
				rnm2 += r[i3][i2][i1]*r[i3][i2][i1];
				//absolute value of floating point
				double a = fabs(r[i3][i2][i1]);
				//max of floating point
				localmax = fmax(localmax,a);
			}
		}
	}
		
	rnmu=fmax(rnmu,localmax);
	
	//mpi reduce
	rnm2 = sqrt( rnm2 / ((double) nx*ny*nz ));
	
	return rnm2;
}

void resid(REAL ***u, REAL*** v, REAL*** r,
           int n1,int n2,int n3, double a[4])
{
    //NOTE: Au = v - r
    //---------------------------------------------------------------------
    //     resid computes the residual:  r = v - Au
    //   
    //     This  implementation costs  15A + 4M per result, where
    //     A and M denote the costs of Addition (or Subtraction) and 
    //     Multiplication, respectively. 
    //     Presuming coefficient a(1) is zero (the NPB assumes this,
    //     but it is thus not a general case), 3A + 1M may be eliminated,
    //     resulting in 12A + 3M.
    //     Note that this vectorizes, and is also fine for cache 
    //     based machines.  
    //---------------------------------------------------------------------
    int i3, i2, i1;

    static bool init = false;

    //Private arrays for each thread
    static double **_u1;
    static double **_u2;
    
    double *u1, *u2;

    if (!init) {
        _u1 = (REAL**)malloc(sizeof(REAL*)*omp_get_max_threads()); 
        _u2 = (double**)malloc(sizeof(REAL*)*omp_get_max_threads());
        for (i1 = 0; i1 < omp_get_max_threads(); i1++) {
            _u1[i1] = malloc(sizeof(REAL)*(nm+1));
            _u2[i1] = malloc(sizeof(REAL)*(nm+1));
        }
        init = true;
    }

    #pragma omp parallel private(i1,i2,i3,u1,u2)
    {
        u1 = _u1[omp_get_thread_num()];
        u2 = _u2[omp_get_thread_num()];

	//cycles through each dimension
        #pragma omp for
        for(i3=1;i3<n3-1;i3++) {
            for(i2=1;i2<n2-1;i2++) {
                for(i1=0;i1<n1;i1++) {
                    u1[i1] = u[i3][i2-1][i1]   + u[i3][i2+1][i1] 
                           + u[i3-1][i2][i1]   + u[i3+1][i2][i1];
                    u2[i1] = u[i3-1][i2-1][i1] + u[i3-1][i2+1][i1]
                           + u[i3+1][i2-1][i1] + u[i3+1][i2+1][i1];
                }
                for(i1=1;i1<n1-1;i1++) {
                    r[i3][i2][i1] = v[i3][i2][i1] 
                        - a[0] * u[i3][i2][i1]
                        //---------------------------------------------------------------------
                        //  Assume a(1) = 0      (Enable 2 lines below if a(1) not= 0)
                        //---------------------------------------------------------------------
                        //    >                     - a[1] * ( u(i1-1,i2,i3) + u(i1+1,i2,i3)
                        //    >                              + u1(i1) )
                        //---------------------------------------------------------------------
                        - a[2] * ( u2[i1] + u1[i1-1] + u1[i1+1] )
                        - a[3] * ( u2[i1-1] + u2[i1+1] );
                }
            }
        }
    }
    //---------------------------------------------------------------------
    //     exchange boundary data
    //---------------------------------------------------------------------
    comm3(r,n1,n2,n3);
}

//Multigrid 3-Dimensions function
//Call: mg3P(u,v,r,a,c,n1,n2,n3);
void mg3P(REAL**** u, REAL*** v, REAL**** r, double a[4], double c[4], int n1,int n2,int n3,int restriction)
{
	//---------------------------------------------------------------------
	//     multigrid V-cycle routine
	//---------------------------------------------------------------------
	//      double precision u(nr),v(nv),r(nr)
	int j,k;
	//---------------------------------------------------------------------
	//     down cycle.
	//     restrict the residual from the fine grid to the coarse
	//---------------------------------------------------------------------
	for(k=lt-1-restriction;k>=1;k--) {
		j = k-1;
		// rprj3(
		// 	r[lt-1-k],		//current level residual
		// 	m1[k],m2[k],m3[k],	//current level m-values
		// 	r[lt-1-j],		//next level residual
		// 	m1[j],m2[j],m3[j]	//next level m-values
		// );
		printf ("M1k: %d, M2k: %d, M3k: %d \n",m1[k],m2[k],m3[k]);
		rprj3_mpi(r[lt-1-k],m1[k],m2[k],m3[k],r[lt-1-j],m1[j],m2[j],m3[j]);
	}
	
	k = 0;
	//---------------------------------------------------------------------
	//     compute an approximate solution on the coarsest grid
	//---------------------------------------------------------------------
	zero3(u[lt-1-k],m1[k],m2[k],m3[k]);
	psinv(r[lt-1-k],u[lt-1-k],m1[k],m2[k],m3[k], c);
	
	//	m1,m2,m3 are ints of maxlevel size- they are defined in globals
	//	n1,n2,n3 are the size of the matrix
	for(k=1;k<lt-1-restriction;k++) {
		j = k-1;
		//---------------------------------------------------------------------
		//        prolongate from level k-1  to k
		//---------------------------------------------------------------------
		zero3(u[lt-1-k],m1[k],m2[k],m3[k]);
		
		interp_mpi(u[lt-1-j],m1[j],m2[j],m3[j],u[lt-1-k], m1[k],m2[k],m3[k]);
		//---------------------------------------------------------------------
		//        compute residual for level k
		//---------------------------------------------------------------------
		resid(u[lt-1-k],r[lt-1-k],r[lt-1-k],m1[k],m2[k],m3[k], a);
		//---------------------------------------------------------------------
		//        apply smoother
		//---------------------------------------------------------------------
		psinv(r[lt-1-k],u[lt-1-k],m1[k],m2[k],m3[k],c);
	}
	j = lt - 2;
	k = lt - 1;
	
	interp_mpi(u[lt-1-j],m1[j],m2[j],m3[j],u[0], n1,n2,n3);
	resid(u[0],v,r[0],n1,n2,n3, a);
	psinv(r[0],u[0],n1,n2,n3,c);
}
/*NOTE: rprj3 projects onto the next coarser grid*/
void rprj3_mpi(REAL*** r, int m1k,int m2k,int m3k, REAL*** s,int m1j,int m2j,int m3j) {
	
	int j3, j2, j1, i3, i2, i1, d1, d2, d3;
	
	double x2,y2;
	double *x1,*y1;
	
	bool init = false;
	//Private arrays for each thread
	static double **_x1;
	static double **_y1;
	//keeps track of the number of ghost cell regions the processor should wait for- first and last  strips require less
	int receiveCount = 2;
	
	
	// REAL ** ghost_data = getGhostCells(r,m1k-2,m2k-2,(m3k-2));
	// int messageSize = (m1k-2)*(m2k-2);
	// 
	// //send the ghost cell data to neighbor processors
	// REAL ** results = exchange_data(ghost_data,messageSize);
	// //TODO : does not seem correct?
	// for (i2 = 1; i2 < m1k-2; i2++) {
	// 	for (i1 = 1; i1 < m2k-2; i1++) {
	// 		r[0][i2][i1] = results[0][(m1k-2)*i2 + i1];
	// 		r[m3k][i2][i1] = results[1][(m1k-2)*i2 + i1];
	// 	}
	// }
	
	//this seems to be incorrect - (it's larger than necessary)
	//it initializes for the largest possible array needed (finest level), instead of adjusting to the current level size
	x1 = malloc(sizeof(double)*(nm+1));
	y1 = malloc(sizeof(double)*(nm+1));
	d1 = (m1k==3) ? 2 : 1;
	d2 = (m2k==3) ? 2 : 1;
	d3 = (m3k==3) ? 2 : 1;
	
	for(j3=2;j3<=m3j-1;j3++) {
		i3 = 2*j3-d3-1;
		for(j2=2;j2<=m2j-1;j2++) {
			i2 = 2*j2-d2-1;
			for(j1=2;j1<=m1j;j1++) {
				i1 = 2*j1-d1-1;
				x1[i1-1] = r[i3][i2-1][i1-1] + r[i3][i2+1][i1-1]
					+ r[i3-1][i2][i1-1] + r[i3+1][i2][i1-1];
				y1[i1-1] = r[i3-1][i2-1][i1-1] + r[i3+1][i2-1][i1-1] 
					+ r[i3-1][i2+1][i1-1] + r[i3+1][i2+1][i1-1];
			}
			
			for(j1=2;j1<=m1j-1;j1++) {
				i1 = 2*j1-d1-1;
				y2 = r[i3-1][i2-1][i1]  + r[i3+1][i2-1][i1] 
					+ r[i3-1][i2+1][i1] + r[i3+1][i2+1][i1];
				x2 = r[i3][i2-1][i1]    + r[i3][i2+1][i1]
					+ r[i3-1][i2][i1]   + r[i3+1][i2][i1];
				s[j3-1][j2-1][j1-1] =
					0.5      *   r[i3][i2][i1]
					+ 0.25   * ( r[i3][i2][i1-1]+r[i3][i2][i1+1]+x2)
					+ 0.125  * ( x1[i1-1] + x1[i1+1] + y2)
					+ 0.0625 * ( y1[i1-1] + y1[i1+1] );
			}
		}
	}
	comm3(s,m1j,m2j,m3j);
}

//r is the finer level of residual error, and s is the coarser level
void rprj3(REAL*** r, int m1k,int m2k,int m3k,
			REAL*** s,int m1j,int m2j,int m3j)
{
	//---------------------------------------------------------------------
	//     rprj3 projects onto the next coarser grid,
	//     using a trilinear Finite Element projection:  s = r' = P r
	//     
	//     This  implementation costs  20A + 4M per result, where
	//     A and M denote the costs of Addition and Multiplication.
	//     Note that this vectorizes, and is also fine for cache 
	//     based machines.  
	//---------------------------------------------------------------------
	//     double precision r(m1k,m2k,m3k), s(m1j,m2j,m3j)
	int j3, j2, j1, i3, i2, i1, d1, d2, d3;
	
	double x2,y2;
	double *x1,*y1;
	
	bool init = false;
	//Private arrays for each thread
	static double **_x1;
	static double **_y1;
	
	//used for openmp- should be re-written and simplified...
	if (!init) {
		_x1 = (double**)malloc(sizeof(double*)*omp_get_max_threads());
		_y1 = (double**)malloc(sizeof(double*)*omp_get_max_threads());
		
		for (i1 = 0; i1 < omp_get_max_threads(); i1++) {
			_x1[i1] = malloc(sizeof(double)*(nm+1));
			_y1[i1] = malloc(sizeof(double)*(nm+1));
		}
		init = true;
	}
	
	d1 = (m1k==3) ? 2 : 1;
	d2 = (m2k==3) ? 2 : 1;
	d3 = (m3k==3) ? 2 : 1;
	
	#pragma omp parallel private(j1,j2,j3,i1,i2,i3,x1,y1,x2,y2)
	{
		x1 = _x1[omp_get_thread_num()];
		y1 = _y1[omp_get_thread_num()];
		
		#pragma omp for
		for(j3=2;j3<=m3j-1;j3++) {
			i3 = 2*j3-d3-1;
			for(j2=2;j2<=m2j-1;j2++) {
				i2 = 2*j2-d2-1;
				for(j1=2;j1<=m1j;j1++) {
					i1 = 2*j1-d1-1;
					x1[i1-1] = r[i3][i2-1][i1-1] + r[i3][i2+1][i1-1]
						+ r[i3-1][i2][i1-1] + r[i3+1][i2][i1-1];
					y1[i1-1] = r[i3-1][i2-1][i1-1] + r[i3+1][i2-1][i1-1] 
						+ r[i3-1][i2+1][i1-1] + r[i3+1][i2+1][i1-1];
				}
				
				for(j1=2;j1<=m1j-1;j1++) {
					i1 = 2*j1-d1-1;
					y2 = r[i3-1][i2-1][i1]  + r[i3+1][i2-1][i1] 
						+ r[i3-1][i2+1][i1] + r[i3+1][i2+1][i1];
					x2 = r[i3][i2-1][i1]    + r[i3][i2+1][i1]
						+ r[i3-1][i2][i1]   + r[i3+1][i2][i1];
					s[j3-1][j2-1][j1-1] =
						0.5      *   r[i3][i2][i1]
						+ 0.25   * ( r[i3][i2][i1-1]+r[i3][i2][i1+1]+x2)
						+ 0.125  * ( x1[i1-1] + x1[i1+1] + y2)
						+ 0.0625 * ( y1[i1-1] + y1[i1+1] );
				}
			}
		}
	}
	comm3(s,m1j,m2j,m3j);
}

//TODO : 
//	-remove omp
//	-exchange ghost cell data
//	-processor 1 compiles results
//z is the coarser level, and u is the finer level (currently zeroed out)
void interp_mpi(REAL ***z, int mm1, int mm2, int mm3, REAL ***u,
        int n1,int n2,int n3 ){
	int i3, i2, i1, d1, d2, d3, t1, t2, t3;
	
	// note that m = 1037 in globals.h but for this only need to be
	// 535 to handle up to 1024^3
	//      integer m
	//      parameter( m=535 )
	int m=535;
	double *z1,*z2,*z3;
	
	static bool init = false;
	
	//Private arrays for each thread
	static double **_z1;
	static double **_z2;
	static double **_z3;
	
	if (!init) {
		_z1 = (double**)malloc(sizeof(double*)*omp_get_max_threads()); 
		_z2 = (double**)malloc(sizeof(double*)*omp_get_max_threads());
		_z3 = (double**)malloc(sizeof(double*)*omp_get_max_threads());
		
		for (i1 = 0; i1 < omp_get_max_threads(); i1++) {
			_z1[i1] = malloc(sizeof(double)*m);
			_z2[i1] = malloc(sizeof(double)*m);
			_z3[i1] = malloc(sizeof(double)*m);
		}
		init = true;
	}
	
	if( n1 != 3 && n2 != 3 && n3 != 3 ) {
		z1 = _z1[omp_get_thread_num()];
		z2 = _z2[omp_get_thread_num()];
		z3 = _z3[omp_get_thread_num()];
		
		for(i3=1;i3<=mm3-1;i3++) {
			for(i2=1;i2<=mm2-1;i2++) {
				//grab values from z
				for(i1=1;i1<=mm1;i1++) {
					z1[i1-1] = z[i3-1][i2][i1-1] + z[i3-1][i2-1][i1-1];
					z2[i1-1] = z[i3][i2-1][i1-1] + z[i3-1][i2-1][i1-1];
					z3[i1-1] = z[i3][i2][i1-1]   + z[i3][i2-1][i1-1] + z1[i1-1];
				}
				//
				for(i1=1;i1<=mm1-1;i1++) {
					u[2*i3-2][2*i2-2][2*i1-2]  += z[i3-1][i2-1][i1-1];
					u[2*i3-2][2*i2-2][2*i1-1]  +=
						0.5*(z[i3-1][i2-1][i1] +  z[i3-1][i2-1][i1-1]);
				}
				for(i1=1;i1<=mm1-1;i1++) {
					u[2*i3-2][2*i2-1][2*i1-2] += 0.5  *  z1[i1-1];
					u[2*i3-2][2*i2-1][2*i1-1] += 0.25 * (z1[i1-1] + z1[i1] );
				}
				for(i1=1;i1<=mm1-1;i1++) {
					u[2*i3-1][2*i2-2][2*i1-2] += 0.5  * z2[i1-1];
					u[2*i3-1][2*i2-2][2*i1-1] += 0.25 *(z2[i1-1] + z2[i1] );
				}
				for(i1=1;i1<=mm1-1;i1++) {
					u[2*i3-1][2*i2-1][2*i1-2] += 0.25*z3[i1-1];
					u[2*i3-1][2*i2-1][2*i1-1] += 0.125*( z3[i1-1] + z3[i1] );
				}
			}
		}
	} else {
		if(n1==3) {
			d1 = 2;
			t1 = 1;
		} else {
			d1 = 1;
			t1 = 0;
		}
		
		if(n2==3) {
			d2 = 2;
			t2 = 1;
		} else {
			d2 = 1;
			t2 = 0;
		}
		
		if(n3==3) {
			d3 = 2;
			t3 = 1;
		} else {
			d3 = 1;
			t3 = 0;
		}
		
		for(i3=1;i3<=mm3-1;i3++) {
			for(i2=1;i2<=mm2-1;i2++) {
				for(i1=1;i1<=mm1-1;i1++) {
					u[2*i3-1-d3][2*i2-1-d2][2*i1-1-d1] +=
						z[i3-1][i2-1][i1-1];
				}
				for(i1=1;i1<=mm1-1;i1++) {
					u[2*i3-1-d3][2*i2-1-d2][2*i1-1-t1] +=
						0.5*(z[i3-1][i2-1][i1] + z[i3-1][i2-1][i1-1]);
				}
			}
			for(i2=1;i2<=mm2-1;i2++) {
				for(i1=1;i1<=mm1-1;i1++) {
					u[2*i3-1-d3][2*i2-1-t2][2*i1-1-d1] +=
						0.5*(z[i3-1][i2][i1-1] + z[i3-1][i2-1][i1-1]);
				}
				for(i1=1;i1<=mm1-1;i1++) {
					u[2*i3-1-d3][2*i2-1-t2][2*i1-1-t1] +=
					0.25*(z[i3-1][i2][i1] + z[i3-1][i2-1][i1]
						+z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
				}
			}
			for(i3=1;i3<=mm3-1;i3++) {
				for(i2=1;i2<=mm2-1;i2++) {
					for(i1=1;i1<=mm1-1;i1++) {
						u[2*i3-1-t3][2*i2-1-d2][2*i1-1-d1] =
							0.5*(z[i3][i2-1][i1-1]+z[i3-1][i2-1][i1-1]);
					}
					for(i1=1;i1<=mm1-1;i1++) {
						u[2*i3-1-t3][2*i2-1-d2][2*i1-1-t1] +=
							0.25*(z[i3][i2-1][i1] + z[i3][i2-1][i1-1]
							+z[i3-1][i2-1][i1]+z[i3-1][i2-1][i1-1]);
					}
				}
				for(i2=1;i2<=mm2-1;i2++) {
					for(i1=1;i1<=mm1-1;i1++) {
						u[2*i3-1-t3][2*i2-1-t2][2*i1-1-d1] +=
							0.25*(z[i3][i2][i1-1]+z[i3][i2-1][i1-1]
							+z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
					}
					for(i1=1;i1<=mm1-1;i1++) {
						u[2*i3-1-t3][2*i2-1-t2][2*i1-1-t1] +=
							0.125*(z[i3][i2][i1]+z[i3][i2-1][i1]
							+z[i3][i2][i1-1]+z[i3][i2-1][i1-1]
							+z[i3-1][i2][i1]+z[i3-1][i2-1][i1]
							+z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
					}
				}
			}
		}
	}
	
	//TODO : 
	// -exchange ghost cell data
	// -remove openmp code
	// REAL ** ghost_data = getGhostCells(u,n1-2,n2-2,(n3-2));
	// int messageSize = (n1-2)*(n2-2);
	// REAL ** results = exchange_data(ghost_data,messageSize);
	// printf(" mm1:  %d\n mm2:   %d\n", mm1-2, mm2-2 );
	// printf(" n1:  %d\n n2:   %d\n", n1, n2 );
	// printf(" -------\n" );
	//TODO : need to look over the mechanics of this- the boundaries are not 0!
	// for (i2 = 1; i2 < n1-2; i2++) {
	// 	for (i1 = 1; i1 < n2-2; i1++) {
	// 		// z[0][i2][i1] = 0.0; - returns bad- shouldn't be 0!
	// 		z[0][i2][i1] += results[0][(mm1-2)*i2 + i1];
	// 		z[mm3-2][i2][i1] += results[1][(mm2-2)*i2 + i1];
	// 	}
	// }
}

//Upward cycle / prolongation
void interp(REAL ***z, int mm1, int mm2, int mm3, REAL ***u,
        int n1,int n2,int n3 )
{
    //---------------------------------------------------------------------
    //     interp adds the trilinear interpolation of the correction
    //     from the coarser grid to the current approximation:  u = u + Qu'
    //     
    //     Observe that this  implementation costs  16A + 4M, where
    //     A and M denote the costs of Addition and Multiplication.  
    //     Note that this vectorizes, and is also fine for cache 
    //     based machines.  Vector machines may get slightly better 
    //     performance however, with 8 separate "do i1" loops, rather than 4.
    //---------------------------------------------------------------------
    //      double precision z(mm1,mm2,mm3),u(n1,n2,n3)
    int i3, i2, i1, d1, d2, d3, t1, t2, t3;

    // note that m = 1037 in globals.h but for this only need to be
    // 535 to handle up to 1024^3
    //      integer m
    //      parameter( m=535 )
    int m=535;
    double *z1,*z2,*z3;

    static bool init = false;

    //Private arrays for each thread
    static double **_z1;
    static double **_z2;
    static double **_z3;

    if (!init) {
        _z1 = (double**)malloc(sizeof(double*)*omp_get_max_threads()); 
        _z2 = (double**)malloc(sizeof(double*)*omp_get_max_threads());
        _z3 = (double**)malloc(sizeof(double*)*omp_get_max_threads());

        for (i1 = 0; i1 < omp_get_max_threads(); i1++) {
            _z1[i1] = malloc(sizeof(double)*m);
            _z2[i1] = malloc(sizeof(double)*m);
            _z3[i1] = malloc(sizeof(double)*m);
        }
        init = true;
    }

    if( n1 != 3 && n2 != 3 && n3 != 3 ) {
        #pragma omp parallel private(i1,i2,i3,z1,z2,z3)
        {
            z1 = _z1[omp_get_thread_num()];
            z2 = _z2[omp_get_thread_num()];
            z3 = _z3[omp_get_thread_num()];

            #pragma omp for
            for(i3=1;i3<=mm3-1;i3++) {
                for(i2=1;i2<=mm2-1;i2++) {
                    for(i1=1;i1<=mm1;i1++) {
                        z1[i1-1] = z[i3-1][i2][i1-1] + z[i3-1][i2-1][i1-1];
                        z2[i1-1] = z[i3][i2-1][i1-1] + z[i3-1][i2-1][i1-1];
                        z3[i1-1] = z[i3][i2][i1-1]   + z[i3][i2-1][i1-1] + z1[i1-1];
                    }

                    for(i1=1;i1<=mm1-1;i1++) {
                        u[2*i3-2][2*i2-2][2*i1-2]  += z[i3-1][i2-1][i1-1];
                        u[2*i3-2][2*i2-2][2*i1-1]  +=
                            0.5*(z[i3-1][i2-1][i1] +  z[i3-1][i2-1][i1-1]);
                    }
                    for(i1=1;i1<=mm1-1;i1++) {
                        u[2*i3-2][2*i2-1][2*i1-2] += 0.5  *  z1[i1-1];
                        u[2*i3-2][2*i2-1][2*i1-1] += 0.25 * (z1[i1-1] + z1[i1] );
                    }
                    for(i1=1;i1<=mm1-1;i1++) {
                        u[2*i3-1][2*i2-2][2*i1-2] += 0.5  * z2[i1-1];
                        u[2*i3-1][2*i2-2][2*i1-1] += 0.25 *(z2[i1-1] + z2[i1] );
                    }
                    for(i1=1;i1<=mm1-1;i1++) {
                        u[2*i3-1][2*i2-1][2*i1-2] += 0.25*z3[i1-1];
                        u[2*i3-1][2*i2-1][2*i1-1] += 0.125*( z3[i1-1] + z3[i1] );
                    }
                }
            }
        }
    } else {
        if(n1==3) {
            d1 = 2;
            t1 = 1;
        } else {
            d1 = 1;
            t1 = 0;
        }

        if(n2==3) {
            d2 = 2;
            t2 = 1;
        } else {
            d2 = 1;
            t2 = 0;
        }

        if(n3==3) {
            d3 = 2;
            t3 = 1;
        } else {
            d3 = 1;
            t3 = 0;
        }

        #pragma omp parallel private(i1,i2,i3)
        {
            #pragma omp for
            for(i3=1;i3<=mm3-1;i3++) {
                for(i2=1;i2<=mm2-1;i2++) {
                    for(i1=1;i1<=mm1-1;i1++) {
                        u[2*i3-1-d3][2*i2-1-d2][2*i1-1-d1] +=
                            z[i3-1][i2-1][i1-1];
                    }
                    for(i1=1;i1<=mm1-1;i1++) {
                        u[2*i3-1-d3][2*i2-1-d2][2*i1-1-t1] +=
                            0.5*(z[i3-1][i2-1][i1] + z[i3-1][i2-1][i1-1]);
                    }
                }
                for(i2=1;i2<=mm2-1;i2++) {
                    for(i1=1;i1<=mm1-1;i1++) {
                        u[2*i3-1-d3][2*i2-1-t2][2*i1-1-d1] +=
                            0.5*(z[i3-1][i2][i1-1] + z[i3-1][i2-1][i1-1]);
                    }
                    for(i1=1;i1<=mm1-1;i1++) {
                        u[2*i3-1-d3][2*i2-1-t2][2*i1-1-t1] +=
                            0.25*(z[i3-1][i2][i1] + z[i3-1][i2-1][i1]
                                    +z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
                    }
                }
            }
            #pragma omp for nowait
            for(i3=1;i3<=mm3-1;i3++) {
                for(i2=1;i2<=mm2-1;i2++) {
                    for(i1=1;i1<=mm1-1;i1++) {
                        u[2*i3-1-t3][2*i2-1-d2][2*i1-1-d1] =
                            0.5*(z[i3][i2-1][i1-1]+z[i3-1][i2-1][i1-1]);
                    }
                    for(i1=1;i1<=mm1-1;i1++) {
                        u[2*i3-1-t3][2*i2-1-d2][2*i1-1-t1] +=
                            0.25*(z[i3][i2-1][i1] + z[i3][i2-1][i1-1]
                                    +z[i3-1][i2-1][i1]+z[i3-1][i2-1][i1-1]);
                    }
                }
                for(i2=1;i2<=mm2-1;i2++) {
                    for(i1=1;i1<=mm1-1;i1++) {
                        u[2*i3-1-t3][2*i2-1-t2][2*i1-1-d1] +=
                            0.25*(z[i3][i2][i1-1]+z[i3][i2-1][i1-1]
                                    +z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
                    }
                    for(i1=1;i1<=mm1-1;i1++) {
                        u[2*i3-1-t3][2*i2-1-t2][2*i1-1-t1] +=
                            0.125*(z[i3][i2][i1]+z[i3][i2-1][i1]
                                    +z[i3][i2][i1-1]+z[i3][i2-1][i1-1]
                                    +z[i3-1][i2][i1]+z[i3-1][i2-1][i1]
                                    +z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
                    }
                }
            }
        }
    }
}

//smoother
void psinv(REAL*** r, REAL*** u, int n1,int n2,int n3, double c[4])
{
    //---------------------------------------------------------------------
    //     psinv applies an approximate inverse as smoother:  u = u + Cr
    //
    //     This  implementation costs  15A + 4M per result, where
    //     A and M denote the costs of Addition and Multiplication.  
    //     Presuming coefficient c(3) is zero (the NPB assumes this,
    //     but it is thus not a general case), 2A + 1M may be eliminated,
    //     resulting in 13A + 3M.
    //     Note that this vectorizes, and is also fine for cache 
    //     based machines.  
    //---------------------------------------------------------------------
    //      double precision u(n1,n2,n3),r(n1,n2,n3),c(0:3)
    int i3, i2, i1;

    double *r1, *r2;

    static bool init = false;

    //Private arrays for each thread
    static double **_r1;
    static double **_r2;
    
    if (!init) {
        _r1 = (double**)malloc(sizeof(double*)*omp_get_max_threads()); 
        _r2 = (double**)malloc(sizeof(double*)*omp_get_max_threads());
        for (i1 = 0; i1 < omp_get_max_threads(); i1++) {
            _r1[i1] = malloc(sizeof(double)*(nm+1));
            _r2[i1] = malloc(sizeof(double)*(nm+1));
        }
        init = true;
    }

    #pragma omp parallel private(i1,i2,i3,r1,r2)
    {
        r1 = _r1[omp_get_thread_num()];
        r2 = _r2[omp_get_thread_num()];

        #pragma omp for
        for(i3=1;i3<n3-1;i3++) {
            for(i2=1;i2<n2-1;i2++) {
                for(i1=0;i1<n1;i1++) {
                    r1[i1] = r[i3][i2-1][i1]+ r[i3][i2+1][i1]
                        + r[i3-1][i2][i1] + r[i3+1][i2][i1];
                    r2[i1] = r[i3-1][i2-1][i1] + r[i3-1][i2+1][i1]
                        + r[i3+1][i2-1][i1] + r[i3+1][i2+1][i1];
                }
                for(i1=1;i1<n1-1;i1++) {
                    u[i3][i2][i1] += 
                        c[0] * r[i3][i2][i1]
                        + c[1] * ( r[i3][i2][i1-1] + r[i3][i2][i1+1]
                                + r1[i1] )
                        + c[2] * ( r2[i1] + r1[i1-1] + r1[i1+1] );
                    //---------------------------------------------------------------------
                    //  Assume c(3) = 0    (Enable line below if c(3) not= 0)
                    //---------------------------------------------------------------------
                    //    >                     + c(3) * ( r2(i1-1) + r2(i1+1) )
                    //---------------------------------------------------------------------
                }
            }
        }
    }
    //---------------------------------------------------------------------
    //     exchange boundary points
    //---------------------------------------------------------------------
    comm3(u,n1,n2,n3);
}

void comm3(REAL*** u,int n1,int n2,int n3)
{
//---------------------------------------------------------------------
//     comm3 organizes the communication on all borders 
//---------------------------------------------------------------------
    int i1, i2, i3;

    #pragma omp parallel private(i1,i2,i3)
    {
        #pragma omp for
        for(i3=1;i3<n3-1;i3++) {
            for(i2=1;i2<n2-1;i2++) {
                u[i3][i2][0] = u[i3][i2][n1-2];
                u[i3][i2][n1-1] = u[i3][i2][1];
            }
        }

        for(i3=1;i3<n3-1;i3++)
        {
            for(i1=0;i1<n1;i1++)
            {
                u[i3][0][i1] = u[i3][n2-2][i1];
                u[i3][n2-1][i1]  = u[i3][1][i1];
            }
        }

		//TODO : this doesn't look like it'll work in multiple processors
        #pragma omp for nowait
        for(i2=0;i2<n2;i2++)
        {
            for(i1=0;i1<n1;i1++)
            {
                u[0][i2][i1] = u[n3-2][i2][i1];
                u[n3-1][i2][i1] = u[1][i2][i1];
            }
        }

    }
}
