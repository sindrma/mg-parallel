#include <assert.h>
#include "utility.h"
#include "ptools_ppf.h"
#include <string.h>

//prints out a 3d matrix
void printMatrix(REAL *** mat,int x,int y,int z){
	int i1,i2,i3;
	for(i3=0;i3<z;i3++){
		printf("\n  z=%d", i3);
		for(i2=0;i2<y;i2++){
			printf("\n");
			for(i1=0;i1<x;i1++){
				printf("%f ", mat[i3][i2][i1]);
			}
		}
	}
}

//merges striped matrices together
REAL *** merge_matrices(REAL*** mat1, REAL *** mat2,int x1, int y1, int z1, int x2, int y2, int z2, bool buffer){
	int start = buffer ? 1 : 0; //offfset by 1 if there is a buffer
	//only expands in the z dimension
	REAL *** result = alloc3D((x1),(y1),(z1+z2));
	int i1,i2,i3;
	//copy matrix 1
	for(i3=buffer;i3<(z1);i3++){
		for(i2=buffer;i2<(y1-buffer);i2++){
			for(i1=buffer;i1<(x1-buffer);i1++){
				result[i3][i2][i1] = mat1[i3][i2][i1];
			}
		}
	}
	//copy matrix 2
	for(i3=z1;i3<(z1+z2-buffer);i3++){
		for(i2=buffer;i2<(y1-buffer);i2++){
			for(i1=buffer;i1<(x1-buffer);i1++){
				result[i3][i2][i1] = mat2[i3 - z1 + buffer][i2][i1];
			}
		}
	}
	return result;
}

//even processors send data first and then receive while odd ones are in the other order.
REAL ** exchange_data(REAL** data,int size){
<<<<<<< Updated upstream
	//REAL ** messages = (REAL**) malloc(sizeof(REAL**)*2);
	//messages[0] = (REAL*) malloc(sizeof(REAL*)*size);
	//messages[1] = (REAL*) malloc(sizeof(REAL*)*size);
=======
>>>>>>> Stashed changes
    
    REAL **messages = alloc2D(2, size);
    
    MPI_Request lr_req = MPI_REQUEST_NULL, rr_req = MPI_REQUEST_NULL; 
    MPI_Request ls_req = MPI_REQUEST_NULL, rs_req = MPI_REQUEST_NULL;
    MPI_Status status;
    
    // initiate receive from the left and send to the left
    if(global_params->mpi_rank != 0){
        //MPI_Request lr_req, rr_req;
        MPI_Irecv(messages[0],size, MPI_DOUBLE, global_params->mpi_rank - 1, 1,MPI_COMM_WORLD, &lr_req);
        MPI_Isend(data[0], size, MPI_DOUBLE, global_params->mpi_rank - 1, 2, MPI_COMM_WORLD, &ls_req);
    }
    
    // initiate receive from the right and send to the right
    if(global_params->mpi_rank != global_params->mpi_size-1){
        //MPI_Request ls_req, rs_req;
        MPI_Irecv(messages[1], size, MPI_DOUBLE, global_params->mpi_rank + 1, 2,MPI_COMM_WORLD,&rr_req);
        MPI_Isend(data[1], size, MPI_DOUBLE, global_params->mpi_rank + 1, 1, MPI_COMM_WORLD,&rs_req);
    }
    
    MPI_Wait(&lr_req, &status);
    MPI_Wait(&ls_req, &status);
    MPI_Wait(&rr_req, &status);
    MPI_Wait(&rs_req, &status);
    
    
	return messages;
}

//returns the ghost cells for a matrix (so that they can be sent to neigbor processors)
//currently returns front and back planes-
REAL ** getGhostCells(REAL*** mat,int x,int y, int z){
	int i2,i1;
	int offset = (x)*(y); //plane size
    
    REAL **ghost_cells = alloc2D(2, offset);
	
	//note: the matrix is already padded with boundary data so it must be offset by 1
	for(i2=0;i2<y;i2++){
		for(i1=0;i1<x;i1++){
			//get first plane
			ghost_cells[0][(i2)*(x) + i1] = mat[0+1][i2][i1];
			//get last plane
			ghost_cells[1][(i2)*(x) + i1] = mat[z-1-1][i2][i1];
		}
	}
	
	return ghost_cells;
}

// REAL ** splitMatrix
//given a 3d matrix, returns the strip that a processor should have
//has the option to maintain boundary points on the result matrix
//TODO : current split values: 258x258x129- should this be 258x258x130?=
//n1,n2,n3 are input
REAL *** splitMatrix(REAL*** mat,int x,int y,int z,int processorID,int size,bool addBuffer){
	int i3,i2,i1;
	int z_size = (z - 2) / size + 2; //130
	int offset = (z_size - 2) * processorID; //128 offset per processor
	//allocate 3D strip matrix
	
	REAL*** strip = alloc3D(x,y,z_size);
	
	for(i3=1;i3<z_size-1;i3++){
		for(i2=1;i2<y-1;i2++){
			for(i1=1;i1<x-1;i1++){
				strip[i3][i2][i1] = mat[i3 + offset][i2][i1];
			}
		}
	}
	return strip;
}
//opposite of flatten matrix
REAL *** unflattenMatrix(REAL* mat,int x,int y,int z){
	REAL *** data = alloc3D(x,y,z);
	int i3,i2,i1;
	// printf(" unflat Values:  %dx%dx%d\n ", x, y, z );
	//map 3D matrix to 1D
	for(i3=0;i3<z;i3++){
		for(i2=0;i2<y;i2++){
			for(i1=0;i1<x;i1++){
				data[i3][i2][i1] = mat[i3*x*y + i2*x + i1];
			}
		}
	}
	return data;
}

/*flatten function: 
	Given matrix pointer, x,y,z sizes, return 1D matrix
*/
REAL * flattenMatrix(REAL*** mat,int x,int y,int z){
	REAL * data = (REAL*) malloc(sizeof(REAL)*x*y*z);
	// printf(" flat Values:  %dx%dx%d\n ", x, y, z );
	int i3,i2,i1;
	//map 3D matrix to 1D
	for(i3=0;i3<z;i3++){
		for(i2=0;i2<y;i2++){
			for(i1=0;i1<x;i1++){
				data[i3*x*y + i2*x + i1] = mat[i3][i2][i1];
			}
		}
	}
	return data;
}

//Fill the 3-dimensional array z with zeroes
void zero3(double ***z,int n1,int n2,int n3)
{
    int i1, i2, i3;

    #pragma omp parallel for private(i1,i2,i3)
    for(i3=0;i3<n3;i3++)
        for(i2=0;i2<n2;i2++)
            for(i1=0;i1<n1;i1++)
                z[i3][i2][i1] = 0.0; //[off+i1+n1*i2+n1*n2*i3]=0.0;

}


//Fill the 3-dimensional array z with zeroes
void zero3old(double z[],int off,int n1,int n2,int n3)
{
    int i1, i2, i3;

    #pragma omp parallel for private(i1,i2,i3)
    for(i3=0;i3<n3;i3++)
        for(i2=0;i2<n2;i2++)
            for(i1=0;i1<n1;i1++)
                z[off+i1+n1*i2+n1*n2*i3]=0.0;
}

//bubble does a bubble sort in direction dir
void bubble(double ten[],int j1[],int j2[],int j3[],int m,int ind )
{
    double temp;
    int i, j_temp=0;

    if( ind == 1 )
    {
        for(i=0;i<m-1;i++)
        {
            if( ten[i+m*ind] > ten[i+1+m*ind] )
            {
                temp = ten[i+1+m*ind];
                ten[i+1+m*ind] = ten[i+m*ind];
                ten[i+m*ind] = temp;

                j_temp           = j1[i+1+m*ind];
                j1[i+1+m*ind] = j1[i+m*ind];
                j1[i+m*ind] = j_temp;

                j_temp           = j2[i+1+m*ind];
                j2[i+1+m*ind] = j2[i+m*ind];
                j2[i+m*ind] = j_temp;

                j_temp           = j3[ i+1+m*ind ];
                j3[i+1+m*ind] = j3[ i+m*ind ];
                j3[i+m*ind] = j_temp;
            }
            else 
            {
                return;
            }
        }
    }
    else
    {
        for(i=0;i<m-1;i++)
        {
            if( ten[i+m*ind] < ten[i+1+m*ind] )
            {
                temp = ten[i+1+m*ind];
                ten[i+1+m*ind] = ten[i+m*ind];
                ten[i+m*ind] = temp;

                j_temp           = j1[i+1+m*ind];
                j1[i+1+m*ind] = j1[i+m*ind];
                j1[i+m*ind] = j_temp;

                j_temp           = j2[i+1+m*ind];
                j2[i+1+m*ind] = j2[i+m*ind];
                j2[i+m*ind] = j_temp;

                j_temp           = j3[ i+1+m*ind ];
                j3[i+1+m*ind] = j3[ i+m*ind ];
                j3[i+m*ind] = j_temp;
            }
            else 
            {
                return;
            }
        }
    }
}

REAL **alloc2D(int n, int m){
    int i,j;
    REAL ** buffer = (REAL**) malloc(sizeof(REAL*)*n);
    for(i=0; i<n; i++){
        buffer[i] = (REAL*) malloc(sizeof(REAL)*m);
		for(j=0; j<m; j++){
			buffer[i][m] = 0.0;
		}
    }
	return buffer;
}

void free2D(REAL** buffer, int n){
    int i;
    for(i=0; i<n; i++){
        free(buffer[i]);
    }
    free(buffer);
}


REAL ***alloc3D(int n, int m,int k)
{
    REAL ***m_buffer=NULL;

    int nx=n, ny=m, nk = k;

    m_buffer = (REAL***)malloc(sizeof(REAL**)* nk);
    assert(m_buffer);  

    REAL** m_tempzy = (REAL**)malloc(sizeof(REAL*)* nk * ny);
    REAL *m_tempzyx = (REAL*)malloc(sizeof(REAL)* nx * ny * nk );
    int z,y; 
    for ( z = 0 ; z < nk ; z++, m_tempzy += ny ) {    
        m_buffer[z] = m_tempzy;
        for ( y = 0 ; y < ny ; y++, m_tempzyx += nx ) {
            m_buffer[z][y] = m_tempzyx;
        }
    }
	zero3(m_buffer,n,m,k);
    return m_buffer;
}

void free3D(REAL*** arr)
{
    free(arr[0][0]);
    free(arr[0]);
    free(arr);
}

void free3D_old(REAL*** arr, int n, int m){
    int i,j;
    for(i=0; i<n; i++){
        for(j=0; j<m; j++){
            free(arr[i][j]);   
        }
        free(arr[i]);
    }
    free(arr);
}


REAL **** allocGrids(size_t depth, size_t n3, size_t n2, size_t n1, size_t pad)
{
    size_t _n1 = n1, _n2 = n2, _n3 = n3;
    long total = 0;
    size_t indexes = 0;
    size_t zsize = 0, zysize = 0;
    long i,d,z,y,x;


    for (i = 0; i < depth; i++)
    {
        total  += (_n1+pad)*(_n2+pad)*(_n3+pad);
        zsize  += (_n1+pad);
        zysize += (_n1+pad)*(_n2+pad);
        _n1 /= 2; _n2 /= 2; _n3 /= 2;
    }


    REAL**** buffer   = (REAL****) malloc(sizeof(REAL***)* depth);
    REAL***  tempdz   = (REAL***)  malloc(sizeof(REAL**) * zsize);
    REAL**   tempdzy  = (REAL**)   malloc(sizeof(REAL*)  * zysize);
    REAL*    tempdzyx = (REAL*)    malloc(sizeof(REAL)   * total);

    _n1 = n1; _n2 = n2; _n3 = n3;

    for (d = 0; d < depth; d++)
    {
        buffer[d] = tempdz;
        for (z = 0; z < _n1+pad; z++, tempdzy += _n2+pad)
        {
            buffer[d][z] = tempdzy;
            for (y = 0; y < _n2+pad; y++, tempdzyx += _n3+pad)
            {
                buffer[d][z][y] = tempdzyx;
            }
        }
        tempdz += _n1+pad;
        _n1 /= 2; _n2 /= 2; _n3 /= 2;
    }
	//zero out arrays
	_n1 = n1; _n2 = n2; _n3 = n3;
	for (d = 0; d < depth; d++){
		for (z = 0; z < _n1+pad; z++){
			for (y = 0; y < _n2+pad; y++){
				for (x = 0; x < _n3+pad; x++){
					buffer[d][z][y][x] = 0.0;
				}
			}
		}
		// zero3(buffer[d],_n3+pad,_n2+pad,_n1+pad);
		_n1 /= 2; _n2 /= 2; _n3 /= 2;
	}

    return buffer;
}

void freeGrids(REAL**** grids)
{
    free(grids[0][0][0]);
    free(grids[0][0]);
    free(grids[0]);
    free(grids);
}

