#include <assert.h>
#include "utility.h"

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

