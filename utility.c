#include <assert.h>
#include "utility.h"

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
	REAL ** messages = (REAL**) malloc(sizeof(REAL**)*2);
	messages[0] = (REAL*) malloc(sizeof(REAL*)*size);
	messages[1] = (REAL*) malloc(sizeof(REAL*)*size);
	if(global_params->mpi_rank%2 == 0){
		if(global_params->mpi_rank != 0){
			//send left data
			MPI_Send(data[0], size, MPI_DOUBLE, global_params->mpi_rank - 1, 1, MPI_COMM_WORLD);
			//receive left data
			MPI_Recv(messages[0],size,MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		if(global_params->mpi_rank != (global_params->mpi_size - 1)){
			//send right data
			MPI_Send(data[1], size, MPI_DOUBLE, global_params->mpi_rank + 1, 1, MPI_COMM_WORLD);
			//receive right data
			MPI_Recv(messages[1],size,MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	} else {
		//receive left data
		MPI_Recv(messages[0],size,MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//send left data
		MPI_Send(data[0], size, MPI_DOUBLE, global_params->mpi_rank - 1, 1, MPI_COMM_WORLD);
		if(global_params->mpi_rank != (global_params->mpi_size - 1)){
			//receive right data
			MPI_Recv(messages[1],size,MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//send right data
			MPI_Send(data[1], size, MPI_DOUBLE, global_params->mpi_rank + 1, 1, MPI_COMM_WORLD);
		}
	}
	return messages;
}

//returns the ghost cells for a matrix (so that they can be sent to neigbor processors)
//currently returns front and back planes-
REAL ** getGhostCells(REAL*** mat,int x,int y, int z){
	int i2,i1;
	int offset = x*y; //plane size
	//in strips the ghost cells consists of two planes
	REAL ** ghost_cells = (REAL**) malloc(sizeof(REAL*)*2);
	ghost_cells[0] = (REAL*) malloc(sizeof(REAL)*offset);
	ghost_cells[1] = (REAL*) malloc(sizeof(REAL)*offset);
	
	//note: the matrix is already padded with boundary data so it must be offset by 1
	for(i2=0;i2<y;i2++){
		for(i1=0;i1<x;i1++){
			//get first plane
			ghost_cells[0][i2*x + i1] = mat[0][i2][i1];
			//get last plane
			ghost_cells[1][i2*x + i1] = mat[z-1][i2][i1];
		}
	}
	
	return ghost_cells;
}

// REAL ** splitMatrix
//given a 3d matrix, returns the strip that a processor should have
//has the option to maintain boundary points on the result matrix
//TODO : current split values: 258x258x129- should this be 258x258x130?=
REAL *** splitMatrix(REAL*** mat,int x,int y,int z,int processorID,int size,bool addBuffer){
	int i3,i2,i1;
	int offset = (z/size) * processorID;
	//allocate 3D strip matrix
	// printf(" split Values:  %dx%dx%d\n ", x, y, z );
	//add padding on both indexes if required
	//don't need to add padding for first and last indices (already there)
	int endOffset = ((processorID != (size-1) || !addBuffer) ? 0 : 1);
	int startOffset = (processorID != 0 || !addBuffer ? 0 : 1);
	int padding = addBuffer ?
		(
			((processorID != (size-1)) ? 1 : 0)
			+ (processorID != 0 ? 1 : 0)
		) : 0;
	
	REAL*** strip = alloc3D(x,y,(z/size + padding));
	
	for(i3=startOffset;i3<(z/size)-endOffset;i3++){
		for(i2=0;i2<y;i2++){
			for(i1=0;i1<x;i1++){
				strip[i3 + (addBuffer ? 1-startOffset : 0)][i2][i1] = mat[i3 + offset][i2][i1];
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

double TestNorm(double r[],int n1,int n2,int n3)
{
    double rnm2=0.0;
    int i1,i2,i3;
    for(i3=1;i3<n3-1;i3++)
        for(i2=1;i2<n2-1;i2++)
            for(i1=1;i1<n1-1;i1++)
                rnm2+=r[i1+n1*(i2+n2*i3)]*r[i1+n1*(i2+n2*i3)];

    rnm2 = sqrt( rnm2 / ((double)n1*n2*n3));
    printf("*****TestNorm  %f\n", rnm2);
    return rnm2;
}

//Fill the 3-dimensional array z with zeroes
void zero3(double ***z,int n1,int n2,int n3)
{
    int i1, i2, i3;

    #pragma omp parallel for private(i1,i2,i3)
    for(i3=0;i3<n3;i3++)
        for(i2=0;i2<n2;i2++)
            for(i1=0;i1<n1;i1++)
                z[i3][i2][i1] = 0.; //[off+i1+n1*i2+n1*n2*i3]=0.0;

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

    return m_buffer;
}

void free3D(REAL*** arr)
{
    free(arr[0][0]);
    free(arr[0]);
    free(arr);
}


REAL **** allocGrids(size_t depth, size_t n1, size_t n2, size_t n3, size_t pad)
{
    size_t _n1 = n1, _n2 = n2, _n3 = n3;
    long total = 0;
    size_t indexes = 0;
    size_t zsize = 0, zysize = 0;
    long i,d,z,y;


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

    return buffer;
}

void freeGrids(REAL**** grids)
{
    free(grids[0][0][0]);
    free(grids[0][0]);
    free(grids[0]);
    free(grids);
}

