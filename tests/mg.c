#include <omp.h>
#include "../functions/setup.h"
#include "../functions/results.h"
#include "../random.h"
#include "../utility.h"
#include "../timer.h"
#include "mg.h"
#include "../mg.h"

//mpirun -np 2 ./mg
//TODO : testExchange isn't correct yet
int main(int argc, const char **argv)
{	
	// test_Flatten_AND_Unflatten(4,4,4);
	
	// testSplit(
	// 	4+2,4+2,4+2,	//size of matrix (x,y,z) plus boundary buffer
	// 	2,		//number of processors
	// 	true 	//use some of the matrix as boundary points
	// );
	// testSplit(
	// 	4,4,4,	//size of matrix (x,y,z) plus boundary buffer
	// 	2,		//number of processors
	// 	false 	//use some of the matrix as boundary points
	// );
	
	// testMerge(
	// 	4+2,4+2,4+2,	//size of matrix (x,y,z) plus boundary buffer
	// 	2,		//number of processors
	// 	true 	//use some of the matrix as boundary points
	// );
	// testMerge(
	// 	4,4,4,	//size of matrix (x,y,z) plus boundary buffer
	// 	2,		//number of processors
	// 	false 	//use some of the matrix as boundary points
	// );
	
	//testGhostCell(4,4,4);
	// testExchange(4,4,4);
	
    
    
	// PPF_Print( MPI_COMM_WORLD, "\n  Test Complete\n" );
	// printf("\n  Test Complete\n");
	
	return 0;
}

void testResid(int argc, const char **argv){
	MPI_Init( &argc, &argv);
	/*local processor setup- 
		-initializes variables that will be the same for every processor
		-parses command line options
	*/
	global_params = setup_local(argc,argv);
	
	//k is the current level. It is passed down through subroutine args and is NOT global. 
	//it is the current iteration.
	int k, it;
	//pointers of pointers in 3-D
	REAL*** u, //approximation matrix
		*** r, //residual error matrix
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
	
	set_a(a,global_params);
	set_c(c,global_params);
	
	k = lt;
	
	grid_t grid;
	setup(&n1, &n2, &n3, &grid);
	
	v = alloc3D(n1, n2, n3);
	gen_v(v,n1,n2,n3,nx[lt-1],ny[lt-1], &grid);
	
	u = alloc3D(n1-2, n2-2, ((n3-2) / global_params->mpi_size + 2 + 2));
	r = alloc3D(n1-2, n2-2, ((n3-2) / global_params->mpi_size + 2 + 2));
	
	zero3(u[0],n1,n2,n3 / global_params->mpi_size); 
	
	//create local v (for residual)
	REAL ***  local_v = splitMatrix(
		v, //matrix to split
		n1,n2,n3, //matrix size
		global_params->mpi_rank,
		global_params->mpi_size,
		false //don't put a buffer on this!
	);
	
	resid(u,local_v,r,n1,n2,(n3-2)/global_params->mpi_size + 2,a);
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void testComm3(int argc, const char **argv){
    int mpi_rank,mpi_size;
	MPI_Init( &argc, &argv);
	/*local processor setup- 
		-initializes variables that will be the same for every processor
		-parses command line options
	*/
	global_params = setup_local(argc,argv);
    
    REAL*** data_mpi = alloc3D(4,4,2);
    REAL*** data_orig = alloc3D(4,4,2);
	//populate matrix (with a buffer)
	for(i3=0;i3<z;i3++){
		for(i2=0;i2<y;i2++){
			for(i1=0;i1<x;i1++){
				data_mpi[i3][i2][i1] = i1+i2+i3;
                data_orig[i3][i2][i1] = i1+i2+i3;
			}
		}
	}
    
    comm3(data_mpi,4,4,4);
    
    bool differ = false;
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void test_Flatten_AND_Unflatten(int x,int y,int z){
	int i1,i2,i3;
	REAL *** data = alloc3D(x,y,z);
	//populate matrix (with a buffer)
	for(i3=0;i3<z;i3++){
		for(i2=0;i2<y;i2++){
			for(i1=0;i1<x;i1++){
				data[i3][i2][i1] = i1+i2+i3;
			}
		}
	}
	printf("  \nOriginal Matrix:\n");
	printMatrix(data,x,y,z);
	
	REAL * message = flattenMatrix(data,x,y,z);
	//print matrix as 1D
	printf("  \nFlattened Matrix:\n");
	for(i1=0;i1<x*y*z;i1++){
		printf(" %f:",message[i1]);
	}
	
	REAL *** strip = unflattenMatrix(message,x,y,z);
	//print 3d matrix
	printf("  \nFlattened/Unflattened Matrix:\n");
	printMatrix(strip,x,y,z);
}

//grabs the first and last z-level of the matrix- this is the ghost cell data
void testGhostCell(int x,int y,int z){
	int i1,i2,i3;
	REAL *** data = generateMatrix(x,y,z,false);
	
	printf("  \nOriginal Matrix:\n");
	printMatrix(data,x,y,z);
	
	//get ghost cell data
	REAL ** ghost_data = getGhostCells(data,x,y,z);
	//print data
	printf("  \nFirst Plane:\n");
	for(i1=0;i1<(x-2)*(y-2);i1++){
		printf(" %f:",ghost_data[0][i1]);
	}
	printf("  \nSecond Plane:\n");
	for(i1=0;i1<(x-2)*(y-2);i1++){
		printf(" %f:",ghost_data[1][i1]);
	}
}

void testExchange(int x,int y,int z){
	int i1,i2,i3;
	REAL *** data = generateMatrix(x,y,z,false);
	
	//setup global params (used for exchanging data)
	int argc_test;
	char ** argv_test;
	const char ** argv;
	MPI_Init( &argc_test, &argv_test );
	global_params = setup_local(argc_test,argv);
	
	if(global_params->mpi_rank == 0){
		printf("  \nOriginal Matrix:\n");
		printMatrix(data,x,y,z);
	}
	
	//get the slice for the matrix
	REAL *** strip = splitMatrix(
		data,x,y,z,
		global_params->mpi_rank,
		global_params->mpi_size,
		false
	);
	
	//get ghost cell data
	REAL ** ghost_data = getGhostCells(strip,x,y,z/global_params->mpi_size);
	
	//exchange data
	int messageSize = (x)*(y);
	REAL ** results = exchange_data(ghost_data,messageSize);
	
	//copy over results
	for (i2 = 0; i2 < x; i2++) {
		for (i1 = 0; i1 < y; i1++) {
			if(global_params->mpi_rank != 0){
				strip[0][i2][i1] = results[0][(x)*i2 + i1];
			}
			if(global_params->mpi_rank != global_params->mpi_size - 1){
				strip[z/global_params->mpi_size - 1][i2][i1] = results[1][(x)*i2 + i1];
			}
		}
	}
	
	//print matrix
	printf("  \nResult Matrix: %d\n",global_params->mpi_rank);
	printMatrix(strip,x,y,z/global_params->mpi_size);
	printf("\n\n");
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}

/*
	creates a matrix serially
	splits it along the z-axis, (dependent on splitting being correct)
	merges results back together again
*/
void testMerge(int x,int y,int z,int num_processors,bool buffered){
	REAL *** data = alloc3D(x,y,z);
	
	int i1,i2,i3;
	int start = buffered ? 1 : 0;
	//populate matrix (with a buffer)
	for(i3=start;i3<z-start;i3++){
		for(i2=start;i2<y-start;i2++){
			for(i1=start;i1<x-start;i1++){
				data[i3][i2][i1] = i1+i2+i3;
			}
		}
	}
	
	printf("  \nOriginal Matrix:\n");
	printMatrix(data,x,y,z);
	
	REAL ****  mat = (REAL****) malloc(sizeof(REAL***)*num_processors);
	//split matrices
	for(i1=0;i1<num_processors;i1++){
		mat[i1] = splitMatrix(data,x,y,z,i1,num_processors,buffered);
	}
	
	//merge matrices
	REAL *** result = mat[0]; //set to first matrix
	for(i1=1;i1<num_processors;i1++){
		result = merge_matrices(
			result,mat[i1],
			x,y,(z/num_processors)*i1,
			x,y,z/num_processors,
			buffered	//1 cell of boundary data in each matrix
		);
	}
	
	//print result
	printf("  \nSplit/Merged Matrix:\n");
	printMatrix(result,x,y,z);
}
/*
testSplit-
	creates a matrix serially, 
	splits it along the z-axis, 
	then examines the pieces to make sure they are correct
*/
void testSplit(int x,int y,int z,int num_processors,bool buffered){
	//create an x by y by z matrix
	REAL *** data = alloc3D(x,y,z);
	
	int i1,i2,i3;
	int start = buffered ? 1 : 0;
	//populate matrix (with a buffer)
	for(i3=start;i3<z-start;i3++){
		for(i2=start;i2<y-start;i2++){
			for(i1=start;i1<x-start;i1++){
				data[i3][i2][i1] = i1+i2+i3;
			}
		}
	}
	
	printf("  Original Matrix:\n");
	printMatrix(data,x,y,z);
	REAL ****  mat = (REAL****) malloc(sizeof(REAL***)*num_processors);
	
	//splt the matrices
	for(i1=0;i1<num_processors;i1++){
		mat[i1] = splitMatrix(data,x,y,z,i1,num_processors,buffered);
	}
	//print split results
	for(i1=0;i1<num_processors;i1++){
		printf("  \n--------Matrix: %d--------\n",i1);
		printMatrix(mat[i1],x,y,z/num_processors + start);
	}
	
	//TODO : compare values instead of printing them out
	
}

//creates a generic matrix
REAL *** generateMatrix(int x,int y,int z, bool buffered){
	int i1,i2,i3;
	REAL *** data = alloc3D(x,y,z);
	int start = buffered ? 1 : 0;
	//populate matrix (with a buffer)
	for(i3=start;i3<z-start;i3++){
		for(i2=start;i2<y-start;i2++){
			for(i1=start;i1<x-start;i1++){
				data[i3][i2][i1] = i1+i2+i3;
			}
		}
	}
	return data;
}