// Matrix multiply device code
#include "mmpy_kernel.h"
#include <assert.h>
#include <math.h>
// #include "utils.h"
#include "types.h"
using namespace std;

/*
	make clean
	make bx=1 by=1 bz=1
	./mmpy
	
	CURRENTLY WORKING ON:
		-comm3_kernel
			-first and last z-indices (corners, rows, & columns)
	
	const unsigned int bx = BLOCKDIM_X,by = BLOCKDIM_Y,bz = BLOCKDIM_Y;
	nvcc   -arch=sm_21 -c --compiler-options -fno-strict-aliasing  -O3 -I/opt/nvidia/latest/cuda/include   -D_DOUBLE --ptxas-options=-v -DBLOCKDIM_X=16 -DBLOCKDIM_Y=16 setGrid.cu
	
	
	
	TODO : Naive versions of
		-interp
			
		-test different block sizes
*/
/*
	Interpolation uses mm1-1 x mm2-1 x mm3-1 threads
	-each thread writes to 8 different pixels (as part of interpolation)
*/
__global__ void interp_kernel(REAL *z, int mm1, int mm2, int mm3, REAL *u, int n1,int n2,int n3 ){
	const unsigned int bx = BLOCKDIM_X,by = BLOCKDIM_Y,bz = BLOCKDIM_Z;
	const unsigned int tx = threadIdx.x,ty = threadIdx.y,tz = threadIdx.z;
	const unsigned int I =  blockIdx.x*bx + tx,
		J =  blockIdx.y*by + ty,
		K =  blockIdx.z*bz + tz;
	
	if((I < mm1) && (J < mm2) && (K < mm3)){
		int j3 = (K+1);
		int j2 = (J+1);
		int j1 = (I+1);
			
		int i3 = 2*j3;
		int i2 = 2*j2;
		int i1 = 2*j1;
	
		//set offset index for i and j
		j3 *= (mm1*mm2);
		j2 *= (mm1);
		j1 *= 1;
		i3 *= (n1*n2);
		i2 *= (n1);
		i1 *= 1;
		
		//offsets for each direction (by 1)
		int i_off_z = (n1*n2);
		int i_off_y = (n1);
		int i_off_x = 1;
		int off_z = (mm1*mm2);
		int off_y = (mm1);
		int off_x = 1;
		
		/*
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
		*/
		
		
		
		u[(i3-2*i_off_z) + (i2-2*i_off_y) + (i1-2*i_off_x)] += z[(j3-off_z) + (j2-off_y) + (j1-off_x)];
		u[(i3-2*i_off_z) + (i2-2*i_off_y) + (i1-i_off_x)] += 0.5 * (z[(j3-off_z) + (j2-off_y) + (j1-off_x)] + z[(j3-off_z) + (j2-off_y) + j1]);
		u[(i3-2*i_off_z) + (i2-i_off_y) + (i1-2*i_off_x)] += 0.5  * (z[(j3-off_z) + j2 + (j1-off_x)] + z[(j3-off_z) + (j2-off_y) + (j1-off_x)]);
		//OLD (ABOVE ^): u[(i3-2*i_off_z) + (i2-i_off_y) + (i1-2*i_off_x)] += 0.5  * (z[j3 + j2 + (j1-off_x)]   + z[j3 + (j2-off_y) + (j1-off_x)] + z[(j3-off_z) + j2 + (j1-off_x)] + z[(j3-off_z) + (j2-off_y) + (j1-off_x)]);
		u[(i3-2*i_off_z) + (i2-i_off_y) + (i1-i_off_x)] += 0.25 * (z[(j3-off_z) + j2 + (j1-off_x)] + z[(j3-off_z) + (j2-off_y) + (j1-off_x)] + z[(j3-off_z)+ j2 + j1] + z[(j3-off_z) + (j2-off_y) + j1]);
		u[(i3-i_off_z) + (i2-2*i_off_y) + (i1-2*i_off_x)] += 0.5  * (z[j3 + (j2-off_y) + (j1-off_x)] + z[(j3-off_z) + (j2-off_y) + (j1-off_x)]);
		// OLD (ABOVE ^):u[(i3-i_off_z) + (i2-2*i_off_y) + (i1-2*i_off_x)] += 0.5  * (z[j3 + (j2-off_y) + (j1-off_x)] + z[(j3-off_z) + (j2-off_y) + (j1-off_x)]);
		u[(i3-i_off_z) + (i2-2*i_off_y) + (i1-i_off_x)] += 0.25 * (z[j3 + (j2-off_y) + (j1-off_x)] + z[(j3-off_z) + (j2-off_y) + (j1-off_x)] 
										+ z[j3 + (j2-off_y) + j1] + z[(j3-off_z) + (j2-off_y) + j1]);
		
		u[(i3-i_off_z) + (i2-i_off_y) + (i1-2*i_off_x)] += 0.25* (z[j3 + j2 + (j1-off_x)] + z[j3 + (j2-off_y) + (j1-off_x)] + z[(j3-off_z) + j2 + (j1-off_x)] + z[(j3-off_z) + (j2-off_y) + (j1-off_x)]);
		u[(i3-i_off_z) + (i2-i_off_y) + (i1-i_off_x)] += 0.125* (z[j3 + j2 + (j1-off_x)] + z[j3 + (j2-off_y) + (j1-off_x)] + z[(j3-off_z) + j2 + (j1-off_x)] 
												+ z[(j3-off_z) + (j2-off_y) + (j1-off_x)] + z[j3 + j2 + j1] + z[i3 + (i2-off_y) + i1]
												+ z[(j3-off_z) + j2 + j1] + z[(j3-off_z) + (j2-off_y) + j1] );
	}
}

__device__ void comm3_kernel(REAL* u,int n1,int n2,int n3, int I, int J, int K) {
	//---------------------------------------------------------------------
	//     comm3_kernel organizes the communication on all borders 
	//---------------------------------------------------------------------
	I++;J++;K++;

	if(I == 1){
		u[K*(n1*n2) + J*(n1) + n1-1] = u[K*(n1*n2) + J*(n1) + 1];
	}
	if(I == n1-2){
		u[K*(n1*n2) + J*(n1)] = u[K*(n1*n2) + J*(n1) + n1-2];
	}

	if(J==1){
		u[K*(n1*n2) + (n1-1)*n2 + I] = u[K*(n1*n2) + (n1) + I];
		if(I == 1){
			u[K*(n1*n2) + (n2-1)*n1 + n1 - 1] = u[K*(n1*n2) + (n2-1)*n1 + 1];
		}
		if(I == n1-2){
			u[K*(n1*n2) + (n2-1)*n1] = u[K*(n1*n2) + (n2-1)*n1 + n1 - 2];
		}
	}
	if(J==n2-2){
		u[K*(n1*n2) + I] = u[K*(n1*n2) + (n1-2)*n2 + I];
		if(I == 1){
			u[K*(n1*n2) + n1 - 1] = u[K*(n1*n2) + 1];
		}
		if(I == n1-2){
			u[K*(n1*n2)] = u[K*(n1*n2) + n1 - 2];
		}
	}

	if(K==1){
		u[(n3-1)*(n1*n2) + J*(n1) + I] = u[(n1*n2) + J*(n1) + I];
		if(I == 1){
			u[(n3-1)*(n1*n2) + J*(n1) + n1 - 1] = u[(n3-1)*(n1*n2) + J*(n1) + 1];
		}
		if(I == n1-2){
			u[(n3-1)*(n1*n2) + J*(n1)] = u[(n3-1)*(n1*n2) + J*(n1) + n1 - 2];
		}
		if(J==1){
			//todo: should switch n2 and n1 in the n1-1 part?
			u[(n3-1)*(n1*n2) + (n1-1)*n2 + I] = u[(n3-1)*(n1*n2) + (n1) + I];
			if(I == 1){
				u[(n3-1)*(n1*n2) + (n2-1)*n1 + n1 - 1] = u[(n3-1)*(n1*n2) + (n2-1)*n1 + 1];
			}
			if(I == n1-2){
				u[(n3-1)*(n1*n2) + (n2-1)*n1] = u[(n3-1)*(n1*n2) + (n2-1)*n1 + n1 - 2];
			}
		}
		if(J==n2-2){
			u[(n3-1)*(n1*n2) + I] = u[(n3-1)*(n1*n2) + (n1-2)*n2 + I];
			if(I == 1){
				u[(n3-1)*(n1*n2) + n1 - 1] = u[(n3-1)*(n1*n2) + 1];
			}
			if(I == n1-2){
				u[(n3-1)*(n1*n2)] = u[(n3-1)*(n1*n2) +  + n1 - 2];
			}
		}
	}
	if(K==n2-2){
		u[J*(n1) + I] = u[(n3-2)*(n1*n2) + J*(n1) + I];
		if(I == 1){
			u[J*(n1) + n1 - 1] = u[J*(n1) + 1];
		}
		if(I == n1-2){
			u[J*(n1)] = u[J*(n1) + n1 - 2];
		}
		if(J==1){
			//todo: should switch n2 and n1 in the n1-1 part?
			u[(n1-1)*n2 + I] = u[(n1) + I];
			if(I == 1){
				u[(n2-1)*n1 + n1 - 1] = u[(n2-1)*n1 + 1];
			}
			if(I == n1-2){
				u[(n2-1)*n1] = u[ (n2-1)*n1 + n1 - 2];
			}
		}
		if(J==n2-2){
			u[I] = u[(n1-2)*n2 + I];
			if(I == 1){
				u[n1 - 1] = u[1];
			}
			if(I == n1-2){
				u[0] = u[n1 - 2];
			}
		}
	}
}

__global__ void resid_kernel(REAL *u, REAL* v, REAL* r, int n1,int n2,int n3, double a[4]){
	const unsigned int bx = BLOCKDIM_X,by = BLOCKDIM_Y,bz = BLOCKDIM_Z;
	const unsigned int tx = threadIdx.x,ty = threadIdx.y,tz = threadIdx.z;
	const unsigned int I =  blockIdx.x*bx + tx,
		J =  blockIdx.y*by + ty,
		K =  blockIdx.z*bz + tz;
	
	if((I < n1) && (J < n2) && (K < n3)){
		int i3 = (K+1)*(n1*n2);
		int i2 = (J+1)*(n1);
		int i1 = (I+1);
		
		int off_z = (n1*n2);
		int off_y = (n1);
		int off_x = 1;
		
		int center = u[i3 + i2 + i1];
		
		int edges =   u[(i3-off_z)+(i2-off_y)+(i1)]	+ u[(i3+off_z)+(i2-off_y)+(i1)]
					+ u[(i3-off_z)+(i2+off_y)+(i1)]	+ u[(i3+off_z)+(i2+off_y)+(i1)]
					+ u[(i3)+(i2-off_y)+(i1-off_x)]	+ u[(i3)+(i2+off_y)+(i1-off_x)]
					+ u[(i3-off_z)+(i2)+(i1-off_x)]	+ u[(i3+off_z)+(i2)+(i1-off_x)]
					+ u[(i3)+(i2-off_y)+(i1+off_x)]	+ u[(i3)+(i2+off_y)+(i1+off_x)]
					+ u[(i3-off_z)+(i2)+(i1+off_x)]	+ u[(i3+off_z)+(i2)+(i1+off_x)];
		
		int corners =     u[(i3-off_z)+(i2-off_y)+(i1-off_x)]	+ u[(i3+off_z)+(i2-off_y)+(i1-off_x)] 
						+ u[(i3-off_z)+(i2+off_y)+(i1-off_x)]	+ u[(i3+off_z)+(i2+off_y)+(i1-off_x)]
						+ u[(i3-off_z)+(i2-off_y)+(i1+off_x)]	+ u[(i3+off_z)+(i2-off_y)+(i1+off_x)] 
						+ u[(i3-off_z)+(i2+off_y)+(i1+off_x)]	+ u[(i3+off_z)+(i2+off_y)+(i1+off_x)];
		//TODO: copy this array over
		// a[0] = -8.0/3.0; 
		// a[1] =  0.0;
		// a[2] =  1.0/6.0; 
		// a[3] =  1.0/12.0;
		
		r[i3 + i2 + i1] = v[i3 + i2 + i1]
			- (-8.0/3.0) * center
			- (1.0/6.0) * edges
			- (1.0/12.0) * corners;
	}
	
	comm3_kernel(r,n1,n2,n3,I,J,K);
}

//TODO : need to copy over c4 instead of using fixed values
__global__ void psinv_kernel(REAL* r, REAL* u, int n1,int n2,int n3, double c[4]){
	const unsigned int bx = BLOCKDIM_X,by = BLOCKDIM_Y,bz = BLOCKDIM_Z;
	const unsigned int tx = threadIdx.x,ty = threadIdx.y,tz = threadIdx.z;
	const unsigned int I =  blockIdx.x*bx + tx,
		J =  blockIdx.y*by + ty,
		K =  blockIdx.z*bz + tz;
	
	if((I < n1) && (J < n2) && (K < n3)){
		int i3 = (K+1)*(n1*n2);
		int i2 = (J+1)*(n1);
		int i1 = (I+1);
		
		int off_z = (n1*n2);
		int off_y = (n1);
		int off_x = 1;
		
		int center = r[i3 + i2 + i1];
		
		int faces =   r[(i3)+(i2-off_y)+(i1)]	+ r[(i3)+(i2+off_y)+(i1)]
					+ r[(i3-off_z)+(i2)+(i1)]	+ r[(i3+off_z)+(i2)+(i1)]
					+ r[(i3)+(i2)+(i1-off_x)]	+ r[(i3)+(i2)+(i1+off_x)];
		
		int edges =   r[(i3-off_z)+(i2-off_y)+(i1)]	+ r[(i3+off_z)+(i2-off_y)+(i1)]
					+ r[(i3-off_z)+(i2+off_y)+(i1)]	+ r[(i3+off_z)+(i2+off_y)+(i1)]
					+ r[(i3)+(i2-off_y)+(i1-off_x)]	+ r[(i3)+(i2+off_y)+(i1-off_x)]
					+ r[(i3-off_z)+(i2)+(i1-off_x)]	+ r[(i3+off_z)+(i2)+(i1-off_x)]
					+ r[(i3)+(i2-off_y)+(i1+off_x)]	+ r[(i3)+(i2+off_y)+(i1+off_x)]
					+ r[(i3-off_z)+(i2)+(i1+off_x)]	+ r[(i3+off_z)+(i2)+(i1+off_x)];
		
		u[i3 + i2 + i1] += (-3.0/17.0) * center
						+  (+1.0/33.0) * faces
						+  (-1.0/61.0) * edges;
		
		// u[i3 + i2 + i1] += c[0] * center
		// 				+c[1] * faces
		// 				+c[2] * edges;
		comm3_kernel(u,n1,n2,n3,I,J,K);
	}
}

__global__ void rprj3_kernel(REAL* r, int m1k,int m2k,int m3k, REAL* s,int m1j,int m2j,int m3j) {
	const unsigned int bx = BLOCKDIM_X,by = BLOCKDIM_Y,bz = BLOCKDIM_Z;
	const unsigned int tx = threadIdx.x,ty = threadIdx.y,tz = threadIdx.z;
	const unsigned int I =  blockIdx.x*bx + tx,
		J =  blockIdx.y*by + ty,
		K =  blockIdx.z*bz + tz;
	
	if((I < m1k) && (J < m2k) && (K < m3k)){
		int j3 = (K+1);
		int j2 = (J+1);
		int j1 = (I+1);
	
		int i3 = 2*j3;
		int i2 = 2*j2;
		int i1 = 2*j1;
	
		//set offset index for i and j
		j3 *= (m1j*m2j);
		j2 *= (m1j);
		j1 *= 1;
		i3 *= (m1k*m2k);
		i2 *= (m1k);
		i1 *= 1;
	
		//offsets for each direction (by 1)
		int off_z = (m1k*m2k);
		int off_y = (m1k);
		int off_x = 1;
	
		int center = r[i3 + i2 + i1];
	
		int faces =   r[(i3)+(i2-off_y)+(i1)]	+ r[(i3)+(i2+off_y)+(i1)]
					+ r[(i3-off_z)+(i2)+(i1)]	+ r[(i3+off_z)+(i2)+(i1)]
					+ r[(i3)+(i2)+(i1-off_x)]	+ r[(i3)+(i2)+(i1+off_x)];
	
		int corners =     r[(i3-off_z)+(i2-off_y)+(i1-off_x)]	+ r[(i3+off_z)+(i2-off_y)+(i1-off_x)] 
						+ r[(i3-off_z)+(i2+off_y)+(i1-off_x)]	+ r[(i3+off_z)+(i2+off_y)+(i1-off_x)]
						+ r[(i3-off_z)+(i2-off_y)+(i1+off_x)]	+ r[(i3+off_z)+(i2-off_y)+(i1+off_x)] 
						+ r[(i3-off_z)+(i2+off_y)+(i1+off_x)]	+ r[(i3+off_z)+(i2+off_y)+(i1+off_x)];
	
		int edges =   r[(i3-off_z)+(i2-off_y)+(i1)]	+ r[(i3+off_z)+(i2-off_y)+(i1)]
					+ r[(i3-off_z)+(i2+off_y)+(i1)]	+ r[(i3+off_z)+(i2+off_y)+(i1)]
					+ r[(i3)+(i2-off_y)+(i1-off_x)]	+ r[(i3)+(i2+off_y)+(i1-off_x)]
					+ r[(i3-off_z)+(i2)+(i1-off_x)]	+ r[(i3+off_z)+(i2)+(i1-off_x)]
					+ r[(i3)+(i2-off_y)+(i1+off_x)]	+ r[(i3)+(i2+off_y)+(i1+off_x)]
					+ r[(i3-off_z)+(i2)+(i1+off_x)]	+ r[(i3+off_z)+(i2)+(i1+off_x)];
	
		s[(j3) + (j2) + (j1)] = 
				0.5	*	center
			+ 0.25	*	faces
			+ 0.125	*	edges
			+ 0.0625*	corners;
	
		comm3_kernel(s,m1j,m2j,m3j,I,J,K);
	}
}

//WRAPPER FUNCTIONS TO BE USED:
void interp_kernel_wrapper(REAL *z, int mm1, int mm2, int mm3, REAL *u, int n1,int n2,int n3 ) {
	// getMean<<<nBlocks, blockSize>>>(devData, pitch, width,height);
	
}

void comm3_kernel_wrapper(REAL* u,int n1,int n2,int n3, int I, int J, int K) {
	
}

void resid_kernel_wrapper(REAL ***u, REAL*** v, REAL*** r, int n1,int n2,int n3) {
	REAL a[4];
	int _ntx = BLOCKDIM_X;
	int _nty = BLOCKDIM_Y;
	int _ntz = BLOCKDIM_Z;
	
	dim3 threads(_ntx,_nty,_ntz);
	//TODO:  -2 should probably go before division
	int numblocksX = n1 / threads.x - 2;
	int numblocksY = n2 / threads.y - 2;
	int numblocksZ = n3 / threads.z - 2;
	
	if( n1 % _ntx != 0  )
		numblocksX++;
	if( n2 % _nty != 0  )
		numblocksY++;
	if( n3 % _ntz != 0  )
		numblocksZ++;
	// printf("NUMBLOCKSX/Y/Z: %d/%d/%d\n",numblocksX,numblocksX,numblocksZ);
	dim3 grid(numblocksX, numblocksY, numblocksZ);
	
	// print configurations
	// printf("n: %d, tx: %d, ty: %d, gridX: %d, gridY: %d\n\n", n, threads.x, threads.y, grid.x, grid.y);
	
	//r will get overwritten
	REAL * d_r_0 = copy_3d(r,n1,n2,n3);
	REAL * d_v_0 = copy_3d(v,n1,n2,n3);
	REAL * d_u_0 = copy_3d(u,n1,n2,n3);
	
	a[0] = -8.0/3.0; 
	a[1] =  0.0;
	a[2] =  1.0/6.0; 
	a[3] =  1.0/12.0;
	
	// (REAL *u, REAL* v, REAL* r, int n1,int n2,int n3, double a[4])
	resid_kernel<<< grid, threads >>>(d_u_0, d_v_0, d_r_0, n1, n2, n3, a);
	// checkCUDAError("Error in rprj3 kernel");
	
	int return_size = (n1)*(n2)*(n3)*sizeof(REAL);
	REAL * result = (REAL*) malloc(return_size);
	
	// copy result from device to host
	cudaMemcpy(result, d_r_0, return_size, cudaMemcpyDeviceToHost);
	// checkCUDAError("Unable to retrieve result from device");
	//unflatten result:
	int i3,i2,i1;
	//map 3D matrix to 1D
	for(i3=0;i3<n3;i3++){
		for(i2=0;i2<n2;i2++){
			for(i1=0;i1<n1;i1++){
				r[i3][i2][i1] = result[i3*n1*n2 + i2*n1 + i1];
			}
		}
	}
	free(result);
	
	assert(cudaSuccess == cudaFree(d_r_0));
	assert(cudaSuccess == cudaFree(d_u_0));
	
	cudaThreadExit();
}

void psinv_kernel_wrapper(REAL* r, REAL* u, int n1,int n2,int n3, double c[4]) {
	
}

void rprj3_kernel_wrapper(REAL*** r, int m1k,int m2k,int m3k, REAL*** s,int m1j,int m2j,int m3j) {
	// setup execution configurations
	int _ntx = BLOCKDIM_X;
	int _nty = BLOCKDIM_Y;
	int _ntz = BLOCKDIM_Z;
	
	dim3 threads(_ntx,_nty,_ntz);
	int numblocksX = (m1k / threads.x - 2) / 2;
	int numblocksY = (m2k / threads.y - 2) / 2;
	int numblocksZ = (m3k / threads.z - 2) / 2;
	
	if( m1k % _ntx != 0  )
		numblocksX++;
	if( m2k % _nty != 0  )
		numblocksY++;
	if( m3k % _ntz != 0  )
		numblocksZ++;
	// printf("NUMBLOCKSX/Y/Z: %d/%d/%d\n",numblocksX,numblocksX,numblocksZ);
	
	dim3 grid(numblocksX, numblocksY, numblocksZ);
	
	// print configurations
	// printf("n: %d, tx: %d, ty: %d, gridX: %d, gridY: %d\n\n", m1k, threads.x, threads.y, grid.x, grid.y);
	
	int n_temp = (m1k - 2) / 2 + 2;
	int next_size = (n_temp)*(n_temp)*(n_temp)*sizeof(REAL);
	
	REAL * d_r_0 = copy_3d(r,m1k,m1k,m1k);
	REAL * d_r_1 = copy_3d(s,n_temp,n_temp,n_temp);
	
	// printf("\ng.x=%d, g.y=%d, g.z=%d\n",grid.x,grid.y,grid.z);
	// printf("t.x=%d, t.y=%d, t.z=%d\n",threads.x,threads.y,threads.z);
	
	rprj3_kernel<<< grid, threads >>>(d_r_0,m1k,m2k,m3k, d_r_1,n_temp,n_temp,n_temp);
	// checkCUDAError("Error in rprj3 kernel");
	
	REAL * result = (REAL*) malloc(next_size);
	
	// copy result from device to host
	cudaMemcpy(result, d_r_1, next_size, cudaMemcpyDeviceToHost);
	// checkCUDAError("Unable to retrieve result from device");
	// REAL *** r_0 = unflattenMatrix(result,n_temp,n_temp,n_temp);
	
	// printMat(r_0,n_temp,n_temp,n_temp);
	// printf("\n");
	// 
	assert(cudaSuccess == cudaFree(d_r_0));
	assert(cudaSuccess == cudaFree(d_r_1));
	
	cudaThreadExit();
}

REAL* copy_3d(REAL *** mat, int n1, int n2, int n3){
	REAL *handle;
	int total_size = (n1)*(n2)*(n3)*sizeof(REAL);
	
	// REAL *temp = flattenMatrix(mat,n_1,n_2,n_3);
	REAL * data = (REAL*) malloc(sizeof(REAL)*n1*n2*n3);
	int i3,i2,i1;
	//map 3D matrix to 1D
	for(i3=0;i3<n3;i3++){
		for(i2=0;i2<n2;i2++){
			for(i1=0;i1<n1;i1++){
				data[i3*n1*n2 + i2*n1 + i1] = mat[i3][i2][i1];
			}
		}
	}
	
	cudaMalloc((void**) &handle, total_size);
	// checkCUDAError("Error allocating device memory in copy_3d");
	cudaMemset((void **) handle,0,total_size);
	// checkCUDAError("Error initializing device memory matrix r[0]");
	cudaMemcpy(handle, data, total_size, cudaMemcpyHostToDevice);
	// checkCUDAError("Error copying matrix to device");
	return handle;
}
