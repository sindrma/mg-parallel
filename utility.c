#include <assert.h>
#include "utility.h"


//setup function - FLAGGED
//NOTE: this is what creates the error in the first place

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

    for(i=0;i<mm;i++)
    {
        ten[i+mm] = 0.0;
        j1[i+mm] = 0;
        j2[i+mm] = 0;
        j3[i+mm] = 0;
        ten[i] = 1.0;
        j1[i] = 0;
        j2[i] = 0;
        j3[i] = 0;
    }

    for(i3=1;i3<n3-1;i3++)
    {
        for(i2=1;i2<n2-1;i2++)
        {
            for(i1=1;i1<n1-1;i1++)
            {
                if( z[i3][i2][i1] > ten[mm] )
                {
                    ten[mm] = z[i3][i2][i1]; 
                    j1[mm] = i1;
                    j2[mm] = i2;
                    j3[mm] = i3;
                    bubble( ten, j1, j2, j3, mm, 1 );
                }
                if( z[i3][i2][i1] < ten[0] )
                {
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
        }
        else
        {
            jg[4*(i+mm)] = 0;
            jg[1+4*(i+mm)] = 0; 
            jg[2+4*(i+mm)] = 0; 
            jg[3+4*(i+mm)] = 0; 
        }         
        ten[i+mm] = best;

        best = z[j3[i0-1]][j2[i0-1]][j1[i0-1]];
        if(best==z[j3[i0-1]][j2[i0-1]][j1[i0-1]])
        {
            jg[4*i] = 0;
            jg[1+4*i] = is1 - 2 + j1[i0-1]; 
            jg[2+4*i] = is2 - 2 + j2[i0-1]; 
            jg[3+4*i] = is3 - 2 + j3[i0-1]; 
            i0 = i0-1;
        }
        else
        {
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

    for(i=0;i<mm;i++)
    {
        ten[i+mm] = 0.0;
        j1[i+mm] = 0;
        j2[i+mm] = 0;
        j3[i+mm] = 0;
        ten[i] = 1.0;
        j1[i] = 0;
        j2[i] = 0;
        j3[i] = 0;
    }

    for(i3=1;i3<n3-1;i3++)
    {
        for(i2=1;i2<n2-1;i2++)
        {
            for(i1=1;i1<n1-1;i1++)
            {
                if( z[i3][i2][i1] > ten[mm] )
                {
                    ten[mm] = z[i3][i2][i1]; 
                    j1[mm] = i1;
                    j2[mm] = i2;
                    j3[mm] = i3;
                    bubble( ten, j1, j2, j3, mm, 1 );
                }
                if( z[i3][i2][i1] < ten[0] )
                {
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
        }
        else
        {
            jg[4*(i+mm)] = 0;
            jg[1+4*(i+mm)] = 0; 
            jg[2+4*(i+mm)] = 0; 
            jg[3+4*(i+mm)] = 0; 
        }         
        ten[i+mm] = best;

        best = z[j3[i0-1]][j2[i0-1]][j1[i0-1]];
        if(best==z[j3[i0-1]][j2[i0-1]][j1[i0-1]])
        {
            jg[4*i] = 0;
            jg[1+4*i] = is1 - 2 + j1[i0-1]; 
            jg[2+4*i] = is2 - 2 + j2[i0-1]; 
            jg[3+4*i] = is3 - 2 + j3[i0-1]; 
            i0 = i0-1;
        }
        else
        {
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

