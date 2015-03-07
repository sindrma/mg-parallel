/*
   !-------------------------------------------------------------------------!
   !                                                                         !
   !         N  A  S     P A R A L L E L         B E N C H M A R K S  3.0    !
   !                                                                         !
   !                        J A V A         V E R S I O N                    !
   !                                                                         !
   !                                  MG                                     !
   !                                                                         !
   !-------------------------------------------------------------------------!
   !                                                                         !
   !    This benchmark is a serial/multithreaded version of the              !
   !    NPB3_0_JAV MG code.                                                  !
   !                                                                         !
   !    Permission to use, copy, distribute and modify this software         !
   !    for any purpose with or without fee is hereby granted.  We           !
   !    request, however, that all derived work reference the NAS            !
   !    Parallel Benchmarks 3.0. This software is provided "as is"           !
   !    without express or implied warranty.                                 !
   !                                                                         !
   !    Information on NPB 3.0, including the Technical Report NAS-02-008    !
   !    "Implementation of the NAS Parallel Benchmarks in Java",             !
   !    original specifications, source code, results and information        !
   !    on how to submit new results, is available at:                       !
   !                                                                         !
   !            http://www.nas.nasa.gov/Software/NPB/                        !
   !                                                                         !
   !    Send comments or suggestions to  npb@nas.nasa.gov                    !
   !                                                                         !
   !           NAS Parallel Benchmarks Group                                 !
   !           NASA Ames Research Center                                     !
   !           Mail Stop: T27A-1                                             !
   !           Moffett Field, CA   94035-1000                                !
   !                                                                         !
   !           E-mail:  npb@nas.nasa.gov                                     !
   !           Fax:     (650) 604-3957                                       !
   !                                                                         !
   !-------------------------------------------------------------------------!
   ! Authors: E. Barszcz                                                     !
   !          P. Frederickson                                                !
   !          A. Woo                                                         !
   !          M. Yarrow                                                      !
   ! Translation to Java and MultiThreaded Code                              !
   !           M. Frumkin                                                    !
   !           M. Schultz                                                    !
   !-------------------------------------------------------------------------!
 */

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
	
	check functions:
		-resid
		-rprj
		-psinv
		-interp
		-comm3
		-exchange
		
	hypothesis: whatever data is in u needs to be checked!
		-copy serial code
			-print r and u at each level
				//n=8 nit=1
				//2.643068e-01 - 4 processors
				//1.042504e-01 - 1 procesor
				
	2 theories on the mpi stampede error-
		1) an omp thread uses a matrix before it has been zeroed out
		2) an mpi buffer is used before it has been allocated
*/

#include <omp.h>
#include "functions/setup.h"
#include "functions/results.h"
#include "random.h"
#include "utility.h"
#include "timer.h"
#include "mg.h"
#include "mmpy_kernel.h"

#define MIN_ELEMS_PER_PROC 2

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
<<<<<<< Updated upstream
	//PPF_Print( MPI_COMM_WORLD, "n1=%d, n2=%d, n3=%d, nxk=%d, nyk=%d\n", n1, n2, n3, nx[lt-1], ny[lt-1]);
	//processor 0 tracks time and sets up solution matrix
    
    if(global_params->mpi_rank==0){
        init_timers();
        timer_start(T_init);
    }
=======
	// PPF_Print( MPI_COMM_WORLD, "n1=%d, n2=%d, n3=%d, nxk=%d, nyk=%d\n", n1, n2, n3, nx[lt-1], ny[lt-1]);
	//processor 0 tracks time and sets up solution matrix
    // if(global_params->mpi_rank==0){
    //     init_timers();
    //     timer_start(T_init);
    // }
>>>>>>> Stashed changes
	
	v = alloc3D(n1, n2, n3);
    
	gen_v(v,n1,n2,n3,nx[lt-1],ny[lt-1], &grid);
	
	//Initialize arrays- currently does this in strips
<<<<<<< Updated upstream
	// printf("N1= %d, N3=%d\n",n1-2,((n3-2) / global_params->mpi_size + 2 + 2));
    //PPF_Print( MPI_COMM_WORLD, "n3=%d\n", n3 );
    // how the matrixes should be split up from the start.
    // TODO: write a new alloc function that is smarter than this one....
    
	//u = allocGrids(lt, (n1-2)/global_params->mpi_size, n2-2, n3-2, 2);
	//r = allocGrids(lt, (n1-2)/global_params->mpi_size, n2-2, n3-2, 2);
    
    // Just for testing, should be removed in the finished design
    u = allocGrids(lt, n1-2, n2-2, n3-2, 2);
    //PPF_Print( MPI_COMM_WORLD, "n1=%d, n2=%d, n3=%d\n", n1, n2, n3 );
	r = allocGrids(lt, n1-2, n2-2, n3-2, 2);
	
	zero3(u[0],n1,n2,(n3-2) / global_params->mpi_size + 2); //zero-out all of u
=======
	printf("N1= %d, N3=%d\n",n1-2,((n3-2) / global_params->mpi_size + 2 + 2));
	u = allocGrids(lt, n1-2, n2-2, (n3-2) / global_params->mpi_size, 2);
	r = allocGrids(lt, n1-2, n2-2, (n3-2) / global_params->mpi_size, 2);
>>>>>>> Stashed changes
	
	//create local v (different strip per processor)
	REAL ***  local_v = splitMatrix(
		v, //matrix to split
		n1,n2,n3, //matrix size
		global_params->mpi_rank,
    	global_params->mpi_size,
		false //don't put a buffer on this!
	);
	//exchange local_v data
<<<<<<< Updated upstream
    if(global_params->mpi_size > 1){
        exchange(local_v,n1,n2,(n3-2)/global_params->mpi_size+2);
    }
	comm3(local_v,n1,n2,(n3-2)/global_params->mpi_size+2);
    /*
    int i1,i2,i3, x, y, z;
	for(i3=0;i3<n3;i3++){
		for(i2=0;i2<n2;i2++){
			for(i1=0;i1<n1;i1++){
                r[0][i3][i2][i1] = i3;
			}
		}
	}
    
    for(i3=0;i3<(n3-2)/2 + 2;i3++){
		for(i2=0;i2<(n2-2)/2 + 2;i2++){
			for(i1=0;i1<(n1-2)/2 + 2;i1++){
                r[1][i3][i2][i1] = i3;
			}
		}
	}
    */
    //if(global_params->mpi_rank==0){
    //    printMatrix(r[0], n1, n2, (n3-2)/global_params->mpi_size + 2);
        //printMatrix(r[0], n1, n2, n3);
    //}
    /*
    int split = 1;
    for(k=lt-1;k>=1;k--) {
        MPI_Barrier(MPI_COMM_WORLD);
        if((m3[k] - 2) < MIN_ELEMS_PER_PROC && global_params->active){
                MPI_Request req = MPI_REQUEST_NULL;
                MPI_Status status;
                REAL* message;
                global_params->mpi_size = global_params->mpi_size - split;
                // if not divisible by split, send of your matrix
                if(global_params->mpi_rank % (split*2) != 0){
                    global_params->active = false;
                    //printf("m1=%d, m2=%d, m3=%d\n", m1[lt-1], m2[lt-1], m3[lt-1]);
                    //printMatrix(r[0], m1[lt-1], m2[lt-1], m3[lt-1]);
                    message = flattenMatrix(r[lt-1-k], m1[k], m2[k], m3[k]);
                    MPI_Isend(message, m1[k]*m2[k]*m3[k], MPI_DOUBLE, 
                              global_params->mpi_rank-split, 4, MPI_COMM_WORLD, &req);
                    MPI_Wait(&req, &status);
                }
                // else recv elements to work on
                else{
                    message = (REAL*) malloc(sizeof(REAL)*m1[k]*m2[k]*m3[k]);
                    MPI_Irecv(message, m1[k]*m2[k]*m3[k]*split, MPI_DOUBLE, 
                              global_params->mpi_rank + split, 4,MPI_COMM_WORLD, &req);
                    MPI_Wait(&req, &status);
                    for(z=1; z<m3[k]; z++){
                        for(y=0; y<m2[k]; y++){
                            for(x=0; x<m1[k]; x++){
                                r[lt-1-k][z + m3[k] - 2][y][x] = message[m2[k]*m1[k]*z + m1[k]*y + x];  
                            }
                        }
                    }
                   
                }
                split *= 2;
                free(message);
            }
         if(global_params->mpi_rank==0) {
             printf("m1=%d, m2=%d, m3=%d\n", m1[k], m2[k], m3[k]);
             printMatrix(r[lt-1-k], m1[k], m2[k], m3[k]*2);
         }
    }
    */
    /*
    if(m3[0]*split > MIN_ELEMS_PER_PROC){
            MPI_Request req = MPI_REQUEST_NULL;
            MPI_Status status;
            REAL* message;
            global_params->mpi_size = global_params->mpi_size + split;
            if(global_params->mpi_rank % (split/2) == 0 && global_params->active){
                // send data
                REAL* data
                data = flattenMatrix(u[lt-1], m1[0], m2[0], m3[0]);
                message = (REAL*)malloc(sizeof(REAL)*m1[0]*m2[0]*(m3[0]/2));
                for(z=1; z<m3[0]+1; z++){
                    for(y=0; y<m2[0]; y++){
                        for(x=0; x<m1[0]; x++){
                            message[m2[0]*m1[0]*z m1[0]*y + x] = data[m2[0]*m1[0]*(z+m3[0]-1) + m1[0]*y + x];  
                        }
                    }
                }
                free(data);
                MPI_Isend(message, m1[0]*m2[0]*(m3[0]/2), MPI_DOUBLE, 
                          global_params->mpi_rank+split, 4, MPI_COMM_WORLD, &req);
                MPI_Wait(&req, &status);
                free(message);
            } else if(global_params->mpi_rank % split == 0 && !global_params->active){
                // receive data
                global_params->active = true;
                message = (REAL*)malloc(sizeof(REAL)*m1[0]*m2[0]*(m3[0]/2));
                MPI_Irecv(message, m1[0]*m2[0]*(m3[0]/2), MPI_DOUBLE, 
                          global_params->mpi_rank - split, 4,MPI_COMM_WORLD, &req);
                MPI_Wait(&req, &status);
                free(message);
            }
            split /= 2;
    }
	*/
	// if(global_params->mpi_rank == 0){
	// 	printMatrix(local_v,n1,n2,(n3-2)/global_params->mpi_size+2);
	// }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(global_params->mpi_rank==0){
        timer_stop(T_init);
        timer_start(T_bench);
    }
    
    
    //PPF_Print( MPI_COMM_WORLD, "n1=%d, n2=%d, n3=%d\n", n1, n2, n3 );
	resid(u[0],local_v,r[0],n1,n2,(n3-2)/global_params->mpi_size + 2,a);
	//     if(global_params->mpi_rank == 1){
	// 	printMatrix(r[0],n1,n2,(n3-2)/global_params->mpi_size + 2);
=======
	if(global_params->mpi_size > 1){
		exchange(local_v,n1,n2,(n3-2)/global_params->mpi_size+2);
	}
	comm3(local_v,n1,n2,(n3-2)/global_params->mpi_size+2);
	    
	//MPI_Barrier(MPI_COMM_WORLD);
	// if(global_params->mpi_rank==0){
	// 	timer_stop(T_init);
	// 	timer_start(T_bench);
>>>>>>> Stashed changes
	// }
	
	resid(u[0],local_v,r[0],n1,n2,(n3-2)/global_params->mpi_size + 2,a);
	// resid_kernel_wrapper(u[0],local_v,r[0],n1,n2,(n3-2)/global_params->mpi_size + 2);
	if(global_params->mpi_rank == 0) {
		printf("\nRESIDUAL:\n");
		printMatrix(r[0],n1,n2,(n3-2)/global_params->mpi_size + 2);
	}
	// error is sometime after here-
	//each processor runs the multigrid algorithm
	for(it=1;it<=nit;it++) {
		//actual call to multigrid
		mg3P(u,local_v,r,a,c,n1,n2,(n3-2)/global_params->mpi_size + 2,0);
		//compute the residual error here...
		//only pass in the spliced portion of v...
		resid(u[0],local_v,r[0],n1,n2,(n3-2)/global_params->mpi_size + 2,a);
	}
<<<<<<< Updated upstream
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(global_params->mpi_rank==0){
        timer_stop(T_bench);
        tinit = timer_elapsed(T_init);
        printf(" Initialization time: %f seconds\n", tinit);
    }
    
    
    //PPF_Print( MPI_COMM_WORLD, "n1=%d, n2=%d, n3=%d, nx=%d, ny=%d, nz=%d\n", n1, n2, n3, nx[lt-1],ny[lt-1],nz[lt-1] );
    rnm2 = norm2u3(r[0], n1, n2, (n3-2)/global_params->mpi_size+2, nx[lt-1],ny[lt-1],nz[lt-1]);
    
    MPI_Barrier(MPI_COMM_WORLD);
    //validates the results and prints to console
    if(global_params->mpi_rank==0){
        double tm = timer_elapsed(T_bench);
        interpret_results(rnm2, global_params, tm);
    }
    
=======
	//     
	    MPI_Barrier(MPI_COMM_WORLD);
	//     if(global_params->mpi_rank==0){
	// //         timer_stop(T_bench);
	// 		tinit = 1.0;
	// //         tinit = timer_elapsed(T_init);
	// //         printf(" Initialization time: %f seconds\n", tinit);
	//     }
	//     
	//     // PPF_Print( MPI_COMM_WORLD, "n1=%d, n2=%d, n3=%d, nx=%d, ny=%d, nz=%d\n", n1, n2, n3, nx[lt-1],ny[lt-1],nz[lt-1] );
	rnm2 = norm2u3(r[0], n1, n2, (n3-2)/global_params->mpi_size+2, nx[lt-1],ny[lt-1],nz[lt-1]);
	// //     
	// //     MPI_Barrier(MPI_COMM_WORLD);
	//     //validates the results and prints to console
	if(global_params->mpi_rank == 0){
	// //         double tm = timer_elapsed(T_bench);
		double tm = 0.0;
		interpret_results(rnm2, global_params, tm, omp_get_max_threads());
	}
>>>>>>> Stashed changes
	free3D(v);
	freeGrids(u);
	freeGrids(r);
	
	MPI_Barrier(MPI_COMM_WORLD);
	//PPF_Print( MPI_COMM_WORLD, "Message from %N - finished\n" );
	
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
	
	for(i3=0;i3<mm*2;i3++){
		j1[i3] = 0;
		j2[i3] = 0;
		j3[i3] = 0;
	}
	
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
	
	REAL local_rnm2 = 0.0, global_rnm2=0.0;
	int i1,i2,i3;
    
    #pragma omp parallel private(i1,i2,i3)
    {
        #pragma omp for reduction (+:local_rnm2)
        for(i3=1;i3<n3-1;i3++)
        {
            for(i2=1;i2<n2-1;i2++)
            {
                for(i1=1;i1<n1-1;i1++)
                {
                    local_rnm2+=r[i3][i2][i1]*r[i3][i2][i1];
                }
            }

        }

    }
	
    MPI_Reduce(&local_rnm2, &global_rnm2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(global_params->mpi_rank==0){
		// printf("\nNZ=%d\nN3=%d",nz,n3);
	    global_rnm2 = sqrt( global_rnm2 / ((double) nx*ny*nz ));
    }
	
	return global_rnm2;
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
	
	static bool resid_init = false;
	
	//Private arrays for each thread
	static double **_u1;
	static double **_u2;
	
	double *u1, *u2;
	
	if (!resid_init) {
		_u1 = (REAL**)malloc(sizeof(REAL*)*omp_get_max_threads()); 
		_u2 = (double**)malloc(sizeof(REAL*)*omp_get_max_threads());
		for (i1 = 0; i1 < omp_get_max_threads(); i1++) {
			_u1[i1] = malloc(sizeof(REAL)*(nm+1));
			_u2[i1] = malloc(sizeof(REAL)*(nm+1));
		}
		resid_init = true;
	}
	//cycles through each dimension
	#pragma omp parallel private(i1,i2,i3,u1,u2)
    {
        u1 = _u1[omp_get_thread_num()];
        u2 = _u2[omp_get_thread_num()];

        #pragma omp for
        for(i3=1;i3<n3-1;i3++)
        {
            for(i2=1;i2<n2-1;i2++)
            {
                for(i1=0;i1<n1;i1++)
                {
                    u1[i1] = u[i3][i2-1][i1]   + u[i3][i2+1][i1] 
                           + u[i3-1][i2][i1]   + u[i3+1][i2][i1];
                    u2[i1] = u[i3-1][i2-1][i1] + u[i3-1][i2+1][i1]
                           + u[i3+1][i2-1][i1] + u[i3+1][i2+1][i1];
                }
                for(i1=1;i1<n1-1;i1++)
                {
<<<<<<< Updated upstream
=======
					//issue: r[1][0][1] = v[1][0][1] 
					//		- a[0] * u[1][0][1]
					//		- a[2] * edges
					//		- a[3] * corners;
>>>>>>> Stashed changes
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
	//		1- accross processors
	//		2- accross strips
	//---------------------------------------------------------------------
	//exchange(r,n1,n2,n3);
	comm3(r,n1,n2,n3);
}

//Multigrid 3-Dimensions function
void mg3P(REAL**** u, REAL*** v, REAL**** r, double a[4], double c[4], int n1,int n2,int n3,int restriction)
{
	//---------------------------------------------------------------------
	//     multigrid V-cycle routine
	//---------------------------------------------------------------------
	//      double precision u(nr),v(nv),r(nr)
	int j,k,x,y,z;
    int split;
    
    split = 1;
	
	//---------------------------------------------------------------------
	//     down cycle.
	//     restrict the residual from the fine grid to the coarse
	//---------------------------------------------------------------------
	//MEMORY TODO : check m1,m2,m3
	for(k=lt-1-restriction;k>=1;k--) {
		j = k-1;
<<<<<<< Updated upstream
        
        // if to few elements per procs in Z, reduce proc count by two.
        /*
        if((m3[k] - 2)*split < MIN_ELEMS_PER_PROC && global_params->active){
            MPI_Request req = MPI_REQUEST_NULL;
            MPI_Status status;
            REAL* message;
            global_params->mpi_size = global_params->mpi_size - split;
            // if not divisible by split, send of your matrix
            if(global_params->mpi_rank % (split*2) != 0){
                global_params->active = false;
                message = flattenMatrix(r[lt-1-k], m1[k], m2[k], m3[k]);
                MPI_Isend(message, m1[k]*m2[k]*m3[k], MPI_DOUBLE, 
                          global_params->mpi_rank-split, 4, MPI_COMM_WORLD, &req);
                MPI_Wait(&req, &status);
            }
            // else recv elements to work on
            else{
                message = (REAL*) malloc(sizeof(REAL)*m1[k]*m2[k]*m3[k]);
                MPI_Irecv(message, m1[k]*m2[k]*m3[k], MPI_DOUBLE, 
                          global_params->mpi_rank + split, 4,MPI_COMM_WORLD, &req);
                MPI_Wait(&req, &status);
                for(z=1; z<m3[k]; z++){
                    for(y=0; y<m2[k]; y++){
                        for(x=0; x<m1[k]; x++){
                            r[lt-1-k][z + m3[k] - 2][y][x] = message[m2[k]*m1[k]*z + m1[k]*y + x];  
                        }
                    }
                }
            }
            split *= 2;
            free(message);
        }
        */
        if(global_params->active){
            rprj3_mpi(r[lt-1-k],m1[k],m2[k],m3[k]*split,r[lt-1-j],m1[j],m2[j],m3[j]);  
        }
		//ERROR HAS HAPPENED BY THIS TIME- it looks like the last matrix is not filled out
		// if(global_params->mpi_rank == 0 && k==lt-1){
		// 	printMatrix(r[lt-1-k],m1[k],m2[k],m3[k]);
		// }
=======
		rprj3_mpi(r[lt-1-k],m1[k],m2[k],m3[k],r[lt-1-j],m1[j],m2[j],m3[j]);
>>>>>>> Stashed changes
	}
	
	k = 0;
	//---------------------------------------------------------------------
	//     compute an approximate solution on the coarsest grid
	//---------------------------------------------------------------------
	// printf("\nY SIZE=%d\nZ SIZE=%d\n",temp,temp);
	zero3(u[lt-1-k],m1[k],m2[k],m3[k]);
	psinv(r[lt-1-k],u[lt-1-k],m1[k],m2[k],m3[k], c);
	
	for(k=1;k<lt-1-restriction;k++) {
		j = k-1;
        // split to more processors when level is reached again
        /*
        if(m3[j]*split > MIN_ELEMS_PER_PROC){
            MPI_Request req = MPI_REQUEST_NULL;
            MPI_Status status;
            REAL* message;
            global_params->mpi_size = global_params->mpi_size + split;
            if(global_params->mpi_rank % (split/2) == 0 && global_params->active){
                // send data
                REAL* data
                data = flattenMatrix(u[lt-1-j], m1[j], m2[j], m3[j]);
                message = (REAL*)malloc(sizeof(REAL)*m1[j]*m2[j]*(m3[j]/2));
                for(z=1; z<m3[k]+1; z++){
                    for(y=0; y<m2[k]; y++){
                        for(x=0; x<m1[k]; x++){
                            message[m2[j]*m1[j]*z m1[j]*y + x] = data[m2[k]*m1[k]*(z+m3[j]-1) + m1[k]*y + x];  
                        }
                    }
                }
                free(data);
                MPI_Isend(message, m1[k]*m2[k]*(m3[k]/2), MPI_DOUBLE, 
                          global_params->mpi_rank+split, 4, MPI_COMM_WORLD, &req);
                MPI_Wait(&req, &status);
                free(message);
            } else if(global_params->mpi_rank % split == 0 && !global_params->active){
                // receive data
                global_params->active = true;
                message = (REAL*)malloc(sizeof(REAL)*m1[j]*m2[j]*(m3[j]/2));
                MPI_Irecv(message, m1[k]*m2[k]*(m3[k]/2), MPI_DOUBLE, 
                          global_params->mpi_rank - split, 4,MPI_COMM_WORLD, &req);
                MPI_Wait(&req, &status);
                free(message);
            }
            split /= 2;
        }
        */
        
		//---------------------------------------------------------------------
		//        prolongate from level k-1  to k
		//---------------------------------------------------------------------
        if(global_params->active){
            zero3(u[lt-1-k],m1[k],m2[k],m3[k]);
            
            interp_mpi(u[lt-1-j],m1[j],m2[j],m3[j]*split,u[lt-1-k], m1[k],m2[k],m3[k]);
            //---------------------------------------------------------------------
            //        compute residual for level k
            //---------------------------------------------------------------------
            resid(u[lt-1-k],r[lt-1-k],r[lt-1-k],m1[k],m2[k],m3[k], a);
            //---------------------------------------------------------------------
            //        apply smoother
            //---------------------------------------------------------------------
            psinv(r[lt-1-k],u[lt-1-k],m1[k],m2[k],m3[k],c);
        }
	}
	j = lt - 2;
	k = lt - 1;
	// 
	interp_mpi(u[lt-1-j],m1[j],m2[j],m3[j],u[0], n1,n2,n3);
	resid(u[0],v,r[0],n1,n2,n3, a);
	psinv(r[0],u[0],n1,n2,n3,c);
}
/*NOTE: rprj3 projects onto the next coarser grid*/
void rprj3_mpi(REAL*** r, int m1k,int m2k,int m3k, REAL*** s,int m1j,int m2j,int m3j) {
	
	int j3, j2, j1, i3, i2, i1, d1, d2, d3;
	
	double x2,y2;
	double *x1,*y1;
	
	static bool rprj3_init = false;
	//Private arrays for each thread
	static double **_x1;
	static double **_y1;
	//keeps track of the number of ghost cell regions the processor should wait for- first and last  strips require less
	int receiveCount = 2;
	
<<<<<<< Updated upstream
	//Exchange boundary data accross processors
	//exchange(r,m1k,m2k,m3k);
	
=======
>>>>>>> Stashed changes
	//this seems to be incorrect - (it's larger than necessary)
	//it initializes for the largest possible array needed (finest level), instead of adjusting to the current level size
	if (!rprj3_init)
    {
        _x1 = (double**)malloc(sizeof(double*)*omp_get_max_threads()); 
        _y1 = (double**)malloc(sizeof(double*)*omp_get_max_threads());

        for (i1 = 0; i1 < omp_get_max_threads(); i1++)
        {
            _x1[i1] = malloc(sizeof(double)*(nm+1));
            _y1[i1] = malloc(sizeof(double)*(nm+1));
<<<<<<< Updated upstream
=======
			for(j1=0;j3<=nm+1;j1++){
				_x1[i1][j1] = 0.0;
			}
>>>>>>> Stashed changes
        }
        rprj3_init = true;
    }
	d1 = (m1k==3) ? 2 : 1;
	d2 = (m2k==3) ? 2 : 1;
	d3 = (m3k==3) ? 2 : 1;
	
	#pragma omp parallel private(j1,j2,j3,i1,i2,i3,x1,y1,x2,y2)
    {
        x1 = _x1[omp_get_thread_num()];
        y1 = _y1[omp_get_thread_num()];

        #pragma omp for
        for(j3=2;j3<=m3j-1;j3++)
        {
            i3 = 2*j3-d3-1;
            for(j2=2;j2<=m2j-1;j2++)
            {
                i2 = 2*j2-d2-1;
                for(j1=2;j1<=m1j;j1++)
                {
                    i1 = 2*j1-d1-1;
                    x1[i1-1] = r[i3][i2-1][i1-1] + r[i3][i2+1][i1-1]
                             + r[i3-1][i2][i1-1] + r[i3+1][i2][i1-1];
                    y1[i1-1] = r[i3-1][i2-1][i1-1] + r[i3+1][i2-1][i1-1] 
                             + r[i3-1][i2+1][i1-1] + r[i3+1][i2+1][i1-1];
                }

                for(j1=2;j1<=m1j-1;j1++)
                {
                    i1 = 2*j1-d1-1;
                    y2 = r[i3-1][i2-1][i1] + r[i3+1][i2-1][i1] 
                       + r[i3-1][i2+1][i1] + r[i3+1][i2+1][i1];
                    x2 = r[i3][i2-1][i1]   + r[i3][i2+1][i1]
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

void exchange(REAL ***r, int n1,int n2,int n3 ){
	REAL ** ghost_data = getGhostCells(r,n1,n2,n3);
	//just send the x*y planes- not including the buffer
	int messageSize = n1*n2;
	int i1,i2;
	REAL ** results = exchange_data(ghost_data,messageSize);
	
	// put results into the r matrix
	for (i2 = 0; i2 < n2; i2++) {
		for (i1 = 0; i1 < n1; i1++) {
			if(global_params->mpi_rank != 0 ){
				r[0][i2][i1] = results[0][(n1)*i2 + i1];
			}
			if(global_params->mpi_rank != global_params->mpi_size - 1 ){
				r[n3-1][i2][i1] = results[1][(n1)*i2 + i1];
			}
		}
	}
    free2D(ghost_data, 2);
    free2D(results, 2);
}

void exchange_add(REAL ***r, int n1,int n2,int n3 ){
	REAL ** ghost_data = getGhostCells(r,n1,n2,n3);
	//just send the x*y planes- not including the buffer
	int messageSize = n1*n2;
	int i1,i2;
	REAL ** results = exchange_data(ghost_data,messageSize);
	
	// put results into the r matrix
	for (i2 = 0; i2 < n2; i2++) {
		for (i1 = 0; i1 < n1; i1++) {
			if(global_params->mpi_rank != 0 ){
				r[0][i2][i1] += results[0][(n1)*i2 + i1];
			}
			if(global_params->mpi_rank != global_params->mpi_size - 1 ){
				r[n3-1][i2][i1] += results[1][(n1)*i2 + i1];
			}
		}
	}
    free2D(ghost_data, 2);
    free2D(results, 2);
}

//z is the coarser level, and u is the finer level
void interp_mpi(REAL ***z, int mm1, int mm2, int mm3, REAL ***u,
        int n1,int n2,int n3 ){
	int i3, i2, i1, d1, d2, d3, t1, t2, t3;
	
	// note that m = 1037 in globals.h but for this only need to be
	// 535 to handle up to 1024^3
	//      integer m
	//      parameter( m=535 )
	int m=535;
	double *z1,*z2,*z3;
	
	static bool interp_init = false;
<<<<<<< Updated upstream
	
	//Exchange boundary data accross processors
	//exchange(z,mm1,mm2,mm3);
=======
>>>>>>> Stashed changes
	
	//Private arrays for each thread
	static double **_z1;
	static double **_z2;
	static double **_z3;
	
	if (!interp_init) {
		_z1 = (double**)malloc(sizeof(double*)*omp_get_max_threads()); 
		_z2 = (double**)malloc(sizeof(double*)*omp_get_max_threads());
		_z3 = (double**)malloc(sizeof(double*)*omp_get_max_threads());
		
		for (i1 = 0; i1 < omp_get_max_threads(); i1++) {
			_z1[i1] = malloc(sizeof(double)*m);
			_z2[i1] = malloc(sizeof(double)*m);
			_z3[i1] = malloc(sizeof(double)*m);
			for (i2 = 0; i2 < m; i2++) {
				_z1[i1][i2] = 0.0;
				_z2[i1][i2] = 0.0;
				_z3[i1][i2] = 0.0;
			}
		}
		interp_init = true;
	}
	
	if( n1 != 3 && n2 != 3 && n3 != 3 ) {
        #pragma omp parallel private(i1,i2,i3,z1,z2,z3)
        {
            z1 = _z1[omp_get_thread_num()];
            z2 = _z2[omp_get_thread_num()];
            z3 = _z3[omp_get_thread_num()];

            #pragma omp for
            for(i3=1;i3<=mm3-1;i3++)
            {
                for(i2=1;i2<=mm2-1;i2++)
                {

                    for(i1=1;i1<=mm1;i1++)
                    {
                        z1[i1-1] = z[i3-1][i2][i1-1] + z[i3-1][i2-1][i1-1];
                        z2[i1-1] = z[i3][i2-1][i1-1] + z[i3-1][i2-1][i1-1];
                        z3[i1-1] = z[i3][i2][i1-1]   + z[i3][i2-1][i1-1] + z1[i1-1];
                    }

                    for(i1=1;i1<=mm1-1;i1++)
                    {
                        u[2*i3-2][2*i2-2][2*i1-2]  += z[i3-1][i2-1][i1-1];
                        u[2*i3-2][2*i2-2][2*i1-1]  +=
                            0.5*(z[i3-1][i2-1][i1] +  z[i3-1][i2-1][i1-1]);
                    }
                    for(i1=1;i1<=mm1-1;i1++)
                    {
                        u[2*i3-2][2*i2-1][2*i1-2] += 0.5  *  z1[i1-1];
                        u[2*i3-2][2*i2-1][2*i1-1] += 0.25 * (z1[i1-1] + z1[i1] );
                    }
                    for(i1=1;i1<=mm1-1;i1++)
                    {
                        u[2*i3-1][2*i2-2][2*i1-2] += 0.5  * z2[i1-1];
                        u[2*i3-1][2*i2-2][2*i1-1] += 0.25 *(z2[i1-1] + z2[i1] );
                    }
                    for(i1=1;i1<=mm1-1;i1++)
                    {
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
            for(i3=1;i3<=mm3-1;i3++)
            {
                for(i2=1;i2<=mm2-1;i2++)
                {
                    for(i1=1;i1<=mm1-1;i1++)
                    {
                        u[2*i3-1-d3][2*i2-1-d2][2*i1-1-d1] +=
                            z[i3-1][i2-1][i1-1];
                    }
                    for(i1=1;i1<=mm1-1;i1++)
                    {
                        u[2*i3-1-d3][2*i2-1-d2][2*i1-1-t1] +=
                            0.5*(z[i3-1][i2-1][i1] + z[i3-1][i2-1][i1-1]);
                    }
                }
                for(i2=1;i2<=mm2-1;i2++)
                {
                    for(i1=1;i1<=mm1-1;i1++)
                    {
                        u[2*i3-1-d3][2*i2-1-t2][2*i1-1-d1] +=
                            0.5*(z[i3-1][i2][i1-1] + z[i3-1][i2-1][i1-1]);
                    }
                    for(i1=1;i1<=mm1-1;i1++)
                    {
                        u[2*i3-1-d3][2*i2-1-t2][2*i1-1-t1] +=
                            0.25*(z[i3-1][i2][i1] + z[i3-1][i2-1][i1]
                                    +z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
                    }
                }
            }
            #pragma omp for nowait
            for(i3=1;i3<=mm3-1;i3++)
            {
                for(i2=1;i2<=mm2-1;i2++)
                {
                    for(i1=1;i1<=mm1-1;i1++)
                    {
                        u[2*i3-1-t3][2*i2-1-d2][2*i1-1-d1] =
                            0.5*(z[i3][i2-1][i1-1]+z[i3-1][i2-1][i1-1]);
                    }
                    for(i1=1;i1<=mm1-1;i1++)
                    {
                        u[2*i3-1-t3][2*i2-1-d2][2*i1-1-t1] +=
                            0.25*(z[i3][i2-1][i1] + z[i3][i2-1][i1-1]
                                    +z[i3-1][i2-1][i1]+z[i3-1][i2-1][i1-1]);
                    }
                }
                for(i2=1;i2<=mm2-1;i2++)
                {
                    for(i1=1;i1<=mm1-1;i1++)
                    {
                        u[2*i3-1-t3][2*i2-1-t2][2*i1-1-d1] +=
                            0.25*(z[i3][i2][i1-1]+z[i3][i2-1][i1-1]
                                    +z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
                    }
                    for(i1=1;i1<=mm1-1;i1++)
                    {
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
	// exchange(u,n1,n2,n3);
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

    static bool psinv_init = false;

    //Private arrays for each thread
    static double **_r1;
    static double **_r2;
    
    if (!psinv_init) {
        _r1 = (double**)malloc(sizeof(double*)*omp_get_max_threads()); 
        _r2 = (double**)malloc(sizeof(double*)*omp_get_max_threads());
        for (i1 = 0; i1 < omp_get_max_threads(); i1++) {
            _r1[i1] = malloc(sizeof(double)*(nm+1));
            _r2[i1] = malloc(sizeof(double)*(nm+1));
        }
        psinv_init = true;
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
	//Exchange boundary data accross processors
<<<<<<< Updated upstream
	//exchange(u,n1,n2,n3);
=======
>>>>>>> Stashed changes
    comm3(u,n1,n2,n3);
}

void comm3(REAL*** u,int n1,int n2,int n3)
{
//---------------------------------------------------------------------
//     comm3 organizes the communication on all borders 
//---------------------------------------------------------------------
	//exchange around the x-axis
	int i1, i2, i3;
	int isLast = (global_params->mpi_rank == global_params->mpi_size - 1 ? 1:0);
	int isFirst = (global_params->mpi_rank == 0 ? 1:0);
    
    if(global_params->mpi_size > 1){
        exchange(u,n1,n2,n3);
    }
    #pragma omp parallel private(i1,i2,i3, isLast, isFirst)
    {
        #pragma omp for
        for(i3=isFirst;i3<n3-isLast;i3++) {
            for(i2=isFirst;i2<n2-isLast;i2++) {
                u[i3][i2][0] = u[i3][i2][n1-2];
                u[i3][i2][n1-1] = u[i3][i2][1];
            }
        }
        for(i3=isFirst;i3<n3-isLast;i3++)
        {
            for(i1=0;i1<n1;i1++)
            {
                u[i3][0][i1] = u[i3][n2-2][i1];
                u[i3][n2-1][i1]  = u[i3][1][i1];
            }
        }
    }
	//Exchange first and last xy-planes between processors
	if(global_params->mpi_size == 1){
		for(i2=0;i2<n2;i2++) {
			for(i1=0;i1<n1;i1++) {
				u[0][i2][i1] = u[n3-2][i2][i1];
				u[n3-1][i2][i1] = u[1][i2][i1];
			}
		}
	} else {
		//use mpi to exchange this data
		if(global_params->mpi_rank == 0 || global_params->mpi_rank == global_params->mpi_size - 1){
			//Similar call to getGhostCells- might be able to reuse it (offset is not the same though, needs to be tested)
            MPI_Request r_req = MPI_REQUEST_NULL, s_req = MPI_REQUEST_NULL; 
            MPI_Status status;
            
			REAL * ghost_cells = (REAL*) malloc(sizeof(REAL)*n1*n2);
			REAL * message = (REAL*) malloc(sizeof(REAL*)*n1*n2);
			
			if(global_params->mpi_rank == 0){
				//get first plane
				for(i2=0;i2<n2;i2++) {
					for(i1=0;i1<n1;i1++) {
						ghost_cells[(i2)*(n1) + i1] = u[1][i2][i1];
					}
				}
				//send/receive front plane
				//MPI_Send(ghost_cells, n1*n2, MPI_DOUBLE, global_params->mpi_size - 1, 3, MPI_COMM_WORLD);
				//MPI_Recv(message,n1*n2,MPI_DOUBLE, global_params->mpi_size - 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                MPI_Irecv(message, n1*n2, MPI_DOUBLE, global_params->mpi_size - 1, 3,MPI_COMM_WORLD, &r_req);
                MPI_Isend(ghost_cells, n1*n2, MPI_DOUBLE, global_params->mpi_size - 1, 3, MPI_COMM_WORLD, &s_req);
				
                MPI_Wait(&r_req, &status);
                MPI_Wait(&s_req, &status);
                
				for(i2=0;i2<n2;i2++) {
					for(i1=0;i1<n1;i1++) {
						u[0][i2][i1] = message[(i2)*(n1) + i1];
					}
				}
			}
            if(global_params->mpi_rank == global_params->mpi_size-1){	
				//get last plane
				for(i2=0;i2<n2;i2++) {
					for(i1=0;i1<n1;i1++) {
						ghost_cells[(i2)*(n1) + i1] = u[n3-2][i2][i1];
					}
				}
				//receive/send back plane
				//MPI_Recv(message,n1*n2,MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				//MPI_Send(ghost_cells, n1*n2, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
                
                MPI_Irecv(message,n1*n2, MPI_DOUBLE, 0, 3,MPI_COMM_WORLD, &r_req);
                MPI_Isend(ghost_cells, n1*n2, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &s_req);
				
                MPI_Wait(&r_req, &status);
                MPI_Wait(&s_req, &status);
                
				for(i2=0;i2<n2;i2++) {
					for(i1=0;i1<n1;i1++) {
						u[n3-1][i2][i1] = message[(i2)*(n1) + i1];
					}
				}
			}
            free(ghost_cells);
            free(message);
		}
        
	}
}
