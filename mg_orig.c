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

#include <omp.h>
#include "globals.h"
#include "results.h"
#include "includes.h"
#include "utility.h"
#include "timer.h"
#include "random.h"
#include "mg.h"


//Some global constants
const char * BMName="MG"; //Benchmark name

results_t results;

int zoff,zsize3,zsize2,zsize1;
int uoff,usize1,usize2,usize3;
int roff,rsize1,rsize2,rsize3;

void print_timers(char **t_names)
{ 
    printf("  SECTION   Time (secs)\n");

    double tmax = timer_elapsed(T_bench);
    if (tmax == 0.0) tmax = 1.0;

    int i;
    for (i=T_bench;i<=T_last;i++)
    {
        double t = timer_elapsed(i);
        if (i==T_resid2)
        {
            t = timer_elapsed(T_resid);
            printf("          --> total mg-resid %6.4f (%6.4f %%)\n", t, t*100./tmax);
        }
        else
        {
            printf("        %s  %6.4f (%6.4f%%)\n", t_names[i], t, t*100./tmax);
        }
    }
}



int main(int argc, const char **argv)
{

    init_globals();
    //k is the current level. It is passed down through subroutine args and is NOT global. 
    //it is the current iteration.
    int k, it;
    
    double tinit, mflops = 0;
    
    REAL**** u,
        **** r,
        ***  v,
            ///*v = malloc(sizeof(double)*nv),
            //*r = malloc(sizeof(double)*nr), 
            a[4], c[4];
    
    double rnm2, epsilon;
    int n1, n2, n3, nit;
    double verify_value;
    bool verified;
    int i;
    
    char *t_names[16]; //Timer names
    
    init_timers();
    timer_start(T_init);

    //Initialize the timer names
    FILE* f1;
    if( (f1 = fopen("timer.flag", "r")) )
    {
        timeron = true;
        t_names[T_init] = strdup("init");
        t_names[T_bench] = strdup("benchmark");
        t_names[T_mg3P] = strdup("mg3P");
        t_names[T_psinv] = strdup("psinv");
        t_names[T_resid] = strdup("resid");
        t_names[T_rprj3] = strdup("rprj3");
        t_names[T_interp] = strdup("interp");
        t_names[T_norm2] = strdup("norm2");
        fclose(f1);
    }
    else
        timeron = false;


    FILE* f2;
    if ( (f2 = fopen("mg.input", "r")) )
    {
        printf("Reading from input file mg.input\n");

        int ret;

        ret = fscanf(f2, "%d", &lt);
        if(lt>maxlevel) {
            printf("lt=%d Maximum allowable=%d\n", lt, maxlevel);
            exit(0);
        }
        ret = fscanf(f2, "%d", &nx[lt-1]);
        ret = fscanf(f2, "%d", &ny[lt-1]);
        ret = fscanf(f2, "%d", &nz[lt-1]);
        ret = fscanf(f2, "%d", &nit);
        printf("lt %d x %d y %d z %d nit %d\n", lt, nx[lt-1], ny[lt-1], nz[lt-1], nit);

        fflush(stdout);

        for (i = 0; i < 8; i++)
            debug_vec[i] = debug_default;
        fclose(f2);
    }
    else
    {
        printf("No input file. Using compiled defaults.\n");

        lt = lt_default;
        nit = nit_default;
        nx[lt-1] = nx_default;
        ny[lt-1] = ny_default;
        nz[lt-1] = nz_default;
        for (i = 0; i < 8; i++)
            debug_vec[i] = debug_default;
    }  

    if (nx[lt-1] != ny[lt-1] || nx[lt-1] != nz[lt-1])
        Class = 'U';
    else if(nx[lt-1]==32&&nit==4 )
        Class = 'S';
    else if( nx[lt-1]==64&&nit==40 )
        Class = 'W';
    else if( nx[lt-1]==256&&nit==20 )
        Class= 'B';
    else if( nx[lt-1]==512&&nit==20 )
        Class = 'C';
    else if( nx[lt-1]==256&&nit==4 )
        Class = 'A';
    else{
        Class = 'U';
    }

    a[0] = -8.0/3.0; 
    a[1] =  0.0;
    a[2] =  1.0/6.0; 
    a[3] =  1.0/12.0;

    if(Class=='A'||Class=='S'||Class=='W') 
    {
        //---------------------------------------------------------------------
        //     Coefficients for the S(a) smoother
        //---------------------------------------------------------------------
        c[0] =  -3.0/8.0;
        c[1] =  +1.0/32.0;
        c[2] =  -1.0/64.0;
        c[3] =   0.0;
    }
    else
    {
        //---------------------------------------------------------------------
        //     Coefficients for the S(b) smoother
        //---------------------------------------------------------------------
        c[0] =  -3.0/17.0;
        c[1] =  +1.0/33.0;
        c[2] =  -1.0/61.0;
        c[3] =   0.0;
    }

    k = lt;

    printf(" NAS Parallel Benchmarks C version\n");
    printf(" Multithreaded Version %s.%c np=%d\n", BMName, CLASS, omp_get_max_threads());

    printf(" Size:  %dx%dx%d\n Iterations:   %d\n", nx[lt-1], ny[lt-1], nz[lt-1], nit );

    //Initialize arrays
    grid_t grid;
    setup(&n1, &n2, &n3, &grid);
    u = allocGrids(lt, n1-2, n2-2, n3-2, 2);
    r = allocGrids(lt, n1-2, n2-2, n3-2, 2);
    v = alloc3D(n1, n2, n3);

    zero3(u[0],n1,n2,n3);
    zran3(v,n1,n2,n3,nx[lt-1],ny[lt-1], &grid);

    resid(u[0],v,r[0],n1,n2,n3,a);

    //--------------------------------------------------------------------
    //    One iteration for startup
    //--------------------------------------------------------------------
    mg3P(u,v,r,a,c,n1,n2,n3);
    resid(u[0],v,r[0],n1,n2,n3,a);
    zero3(u[0],n1,n2,n3);
    zran3(v,n1,n2,n3,nx[lt-1],ny[lt-1], &grid);

    timer_stop(T_init);
    timer_start(T_bench);     

    if (timeron) timer_start(T_resid2);
    resid(u[0],v,r[0],n1,n2,n3,a);
    if (timeron) timer_stop(T_resid2);
    for(it=1;it<=nit;it++)
    {
        if (timeron) timer_start(T_mg3P);
        mg3P(u,v,r,a,c,n1,n2,n3);
        if (timeron) timer_stop(T_mg3P);

        if (timeron) timer_start(T_resid2);
        resid(u[0],v,r[0],n1,n2,n3,a);
        if (timeron) timer_stop(T_resid2);
    }
    timer_stop(T_bench);

    tinit = timer_elapsed(T_init);
    printf(" Initialization time: %f seconds\n", tinit);

    rnm2=norm2u3(r[0],n1,n2,n3,nx[lt-1],ny[lt-1],nz[lt-1]);
    double tm = timer_elapsed(T_bench);

    verify_value=0.0;
    epsilon = 1.0E-8;
    if (CLASS != 'U') 
    {
        if(CLASS=='S') 
            verify_value = 0.530770700573E-4;
        else if(CLASS=='W') 
            verify_value = 0.250391406439E-17; 
        else if(CLASS=='A') 
            verify_value = 0.2433365309E-5;
        else if(CLASS=='B') 
            verify_value = 0.180056440132E-5;
        else if(CLASS=='C') 
            verify_value = 0.570674826298E-6;
        printf("class %c\n", CLASS);
        printf(" L2 Norm is %e\n", rnm2);
        if(fabs( rnm2 - verify_value ) < epsilon ) 
        {
            verified = 1;
            printf(" Deviation is   %e\n", (rnm2 - verify_value));
        }
        else
        {
            verified = 0;
            printf(" The correct L2 Norm is %e\n", verify_value);
        }
    }
    else
    {
        verified = -1;
    }

    print_verification(CLASS,verified,BMName); 

    if( tm > 0.0 ) 
    {
        mflops = 58.0*nx[lt-1]*ny[lt-1]*nz[lt-1];
        mflops *= nit / (tm*1000000.0);
    }

    set_results(&results, "MG",
            CLASS,
            nx[lt-1],
            ny[lt-1],
            nz[lt-1],
            nit,
            tm,
            mflops,
            (const char*)"floating point",
            verified,
            omp_get_max_threads());

    print_results(&results, stdout);
    if (timeron) print_timers(t_names);

    freeGrids(u);
    freeGrids(r);
    free(v);
//    free(r);

    return 0;
}

void setup(int *n1, int *n2, int *n3, grid_t* grid)
{
    int j, k;

    int ax;
    int size1=3, size2=10;
    int *mi = malloc(sizeof(int)*size1*size2);
    int *ng = malloc(sizeof(int)*size1*size2);

    ng[  (lt-1)*size1]=nx[lt-1];
    ng[1+(lt-1)*size1]=ny[lt-1];
    ng[2+(lt-1)*size1]=nz[lt-1];

    for(ax=0;ax<size1;ax++)
        for(k=lt-2;k>=0;k--)
            ng[ax+k*size1]=ng[ax+(k+1)*size1]/2;

    for(k=lt-2;k>=0;k--)
    {
        nx[k]=ng[  k*size1];
        ny[k]=ng[1+k*size1];
        nz[k]=ng[2+k*size1];
    }

    for(k=lt-1;k>=0;k--)
    {
        for(ax=0;ax<size1;ax++)
        {
            mi[ax+k*size1] = 2 + ng[ax+k*size1];
        }
        m1[k]=mi[k*size1];
        m2[k]=mi[1+k*size1];
        m3[k]=mi[2+k*size1];
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
    *n3= 3 + grid->ie3 - grid->is3;

    ir[lt-1]=0;
    for(j = lt-2;j>=0;j--)
    {
        ir[j]=ir[j+1]+m1[j+1]*m2[j+1]*m3[j+1];
    }

    free(mi);
    free(ng);
}

//---------------------------------------------------------------------
//     zran3  loads +1 at ten randomly chosen points,
//     loads -1 at a different ten random points,
//     and zero elsewhere.
//---------------------------------------------------------------------
void zran3(REAL*** z,int n1,int n2,int n3,int nx,int ny, grid_t* grid)
{
    int mm=10;
    double *ten= malloc(sizeof(double)*mm*2);
    
    int *j1 = malloc(sizeof(int)*mm*2), 
        *j2 = malloc(sizeof(int)*mm*2),
        *j3 = malloc(sizeof(int)*mm*2);
    int *jg = malloc(sizeof(int)*4*mm*2);

    int is1 = grid->is1, is2 = grid->is2, is3 = grid->is3, ie1 = grid->ie1, ie2 = grid->ie2, ie3 = grid->ie3;
    int i0, m0, m1, i1, i2, i3, d1, e1, e2, e3;
    double xx, x0, x1, a1, a2, ai, best;
    int i;
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
    free(jg);
    free(ten);

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
    if (timeron) timer_start(T_norm2);      
    double rnmu = 0.0;
    double rnm2=0.0;
    int i1,i2,i3;

    double localmax;

    #pragma omp parallel private(i1,i2,i3,localmax)
    {
        localmax = 0;

        #pragma omp for reduction (+:rnm2)
        for(i3=1;i3<n3-1;i3++)
        {
            double localmax = 0;
            for(i2=1;i2<n2-1;i2++)
            {
                for(i1=1;i1<n1-1;i1++)
                {
                    rnm2+=r[i3][i2][i1]*r[i3][i2][i1];
                    double a=fabs(r[i3][i2][i1]);
                    localmax=fmax(localmax,a);
                }
            }

        }

        #pragma omp critical
        {
            rnmu=fmax(rnmu,localmax);
        }
    }

    rnm2=sqrt( rnm2 / ((double) nx*ny*nz ));
    if (timeron) timer_stop(T_norm2);

    return rnm2;
}

void resid(REAL ***u, REAL*** v, REAL*** r,
           int n1,int n2,int n3, double a[4])
{
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

    if (!init)
    {
        _u1 = (REAL**)malloc(sizeof(REAL*)*omp_get_max_threads()); 
        _u2 = (double**)malloc(sizeof(REAL*)*omp_get_max_threads());
        for (i1 = 0; i1 < omp_get_max_threads(); i1++)
        {
            _u1[i1] = malloc(sizeof(REAL)*(nm+1));
            _u2[i1] = malloc(sizeof(REAL)*(nm+1));
        }
        init = true;
    }

    if (timeron) timer_start(T_resid);
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
    if (timeron) timer_stop(T_resid);
}

void mg3P(REAL**** u, REAL*** v, REAL**** r, double a[4], double c[4], int n1,int n2,int n3)
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
    for(k=lt-1;k>=1;k--)
    {
        j = k-1;
        rprj3(r[lt-1-k],m1[k],m2[k],m3[k],r[lt-1-j],m1[j],m2[j],m3[j]);
    }
    k = 0;
    //---------------------------------------------------------------------
    //     compute an approximate solution on the coarsest grid
    //---------------------------------------------------------------------
    zero3(u[lt-1-k],m1[k],m2[k],m3[k]);
    psinv(r[lt-1-k],u[lt-1-k],m1[k],m2[k],m3[k], c);
    for(k=1;k<lt-1;k++)
    {     
        j = k-1;
        //---------------------------------------------------------------------
        //        prolongate from level k-1  to k
        //---------------------------------------------------------------------
        zero3(u[lt-1-k],m1[k],m2[k],m3[k]);
        interp(u[lt-1-j],m1[j],m2[j],m3[j],u[lt-1-k], m1[k],m2[k],m3[k]);
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
    k = lt-1;
    interp(u[lt-1-j],m1[j],m2[j],m3[j],u[0], n1,n2,n3);
    resid(u[0],v,r[0],n1,n2,n3, a);
    psinv(r[0],u[0],n1,n2,n3,c);
}

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
    //      double precision r(m1k,m2k,m3k), s(m1j,m2j,m3j)
    int j3, j2, j1, i3, i2, i1, d1, d2, d3;

    double x2,y2;
    double *x1,*y1;
    
    bool init = false;
    //Private arrays for each thread
    static double **_x1;
    static double **_y1;

    if (!init)
    {
        _x1 = (double**)malloc(sizeof(double*)*omp_get_max_threads()); 
        _y1 = (double**)malloc(sizeof(double*)*omp_get_max_threads());

        for (i1 = 0; i1 < omp_get_max_threads(); i1++)
        {
            _x1[i1] = malloc(sizeof(double)*(nm+1));
            _y1[i1] = malloc(sizeof(double)*(nm+1));
        }
        init = true;
    }

    if (timeron) timer_start(T_rprj3);
    if(m1k==3)
        d1 = 2;
    else
        d1 = 1;

    if(m2k==3)
        d2 = 2;
    else
        d2 = 1;

    if(m3k==3)
        d3 = 2;
    else
        d3 = 1;

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
    if (timeron) timer_stop(T_rprj3);
}

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

    if (!init)
    {
        _z1 = (double**)malloc(sizeof(double*)*omp_get_max_threads()); 
        _z2 = (double**)malloc(sizeof(double*)*omp_get_max_threads());
        _z3 = (double**)malloc(sizeof(double*)*omp_get_max_threads());

        for (i1 = 0; i1 < omp_get_max_threads(); i1++)
        {
            _z1[i1] = malloc(sizeof(double)*m);
            _z2[i1] = malloc(sizeof(double)*m);
            _z3[i1] = malloc(sizeof(double)*m);
        }
        init = true;
    }

    if (timeron) timer_start(T_interp);
    if( n1 != 3 && n2 != 3 && n3 != 3 )
    {
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
    }
    else
    {

        if(n1==3)
        {
            d1 = 2;
            t1 = 1;
        }else{
            d1 = 1;
            t1 = 0;
        }

        if(n2==3)
        {
            d2 = 2;
            t2 = 1;
        }else{
            d2 = 1;
            t2 = 0;
        }

        if(n3==3)
        {
            d3 = 2;
            t3 = 1;
        }else{
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
    if (timeron) timer_stop(T_interp);
}

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
    
    if (!init)
    {
        _r1 = (double**)malloc(sizeof(double*)*omp_get_max_threads()); 
        _r2 = (double**)malloc(sizeof(double*)*omp_get_max_threads());
        for (i1 = 0; i1 < omp_get_max_threads(); i1++)
        {
            _r1[i1] = malloc(sizeof(double)*(nm+1));
            _r2[i1] = malloc(sizeof(double)*(nm+1));
        }
        init = true;
    }
    if (timeron) timer_start(T_psinv);

    #pragma omp parallel private(i1,i2,i3,r1,r2)
    {
        r1 = _r1[omp_get_thread_num()];
        r2 = _r2[omp_get_thread_num()];

        #pragma omp for
        for(i3=1;i3<n3-1;i3++)
        {
            for(i2=1;i2<n2-1;i2++)
            {
                for(i1=0;i1<n1;i1++)
                {
                    r1[i1] = r[i3][i2-1][i1]+ r[i3][i2+1][i1]
                        + r[i3-1][i2][i1] + r[i3+1][i2][i1];
                    r2[i1] = r[i3-1][i2-1][i1] + r[i3-1][i2+1][i1]
                        + r[i3+1][i2-1][i1] + r[i3+1][i2+1][i1];
                }
                for(i1=1;i1<n1-1;i1++)
                {
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
    if (timeron) timer_stop(T_psinv);
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
        for(i3=1;i3<n3-1;i3++)
        {
            for(i2=1;i2<n2-1;i2++)
            {
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
