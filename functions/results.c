#include "results.h"
#include "../utility.h"

void init_results(results_t* res)
{
    memset(res, 0, sizeof(results_t));
    res->clss = 'S';
    res->optype = strdup("floating point");
}

void print_results(results_t* res)
{
    printf("***** NAS Parallel Benchmarks C version %s ****\n", res->name);
    printf("* Class             = %c\n", res->clss);

    if(res->n2 == 0 && res->n3 == 0 ) {
        printf("* Size              = %d\n", res->n1);
    } else {
        printf("* Size              = %d X %d X %d\n", res->n1, res->n2, res->n3);
    }
    printf("* Iterations        = %d\n", res->niter);
    printf("* Time in seconds   = %6.4f\n", res->time);
    printf("* ACCTime           = %e\n", res->acctime);
    printf("* Mops total        = %6.4f\n", res->mops);
    printf("* Operation type    = %s\n", res->optype);
    
    printf("* Verification      = ");
    if (res->verified==1) printf("Successful\n");
    else if(res->verified==0) printf("Failed\n");
    else printf("Not Performed\n");

    printf("* Threads requested = %d\n", res->numthreads);
    
    printf("* Please send all errors/feedbacks to:\n");
    printf("* NPB Working Team\n");
    printf("* npb@nas.nasa.gov\n");
    printf("***************************************************************\n");
}

void print_verification(char clss, int verified, const char* BMName)
{
    if (clss == 'U' || verified == -1) {
        verified = -1;
        printf(" Problem size unknown\n");
        printf("%s.%c: Verification Not Performed\n", BMName, clss);
    } else if (verified==1) 
        printf("%s.%c: Verification Successful\n", BMName, clss);
    else
        printf("%s.%c: Verification Failed\n", BMName, clss);
}

void set_results(results_t* res, const char* bname, char CLSS, int bn1, int bn2, int bn3, int bniter, double btime, double bmops, const char *boptype, int passed_verification, int num_threads)
{
	res->name=strdup(bname);
	res->clss=CLSS;
	res->n1=bn1;
	res->n2=bn2;
	res->n3=bn3;
	res->niter=bniter;
	res->time=btime;
	res->mops=bmops;
	res->optype=strdup(boptype);
	res->verified=passed_verification;
	res->numthreads=num_threads;
}

void interpret_results(double rnm2,struct params* global, double tm, int numThreads){
	results_t res;
	
	bool verified = 0;
	double verify_value = 0.0;
	double epsilon = 1.0E-8;
	
	printf(" L2 Norm is %e\n", rnm2);
	if (global->class != 'U') {
		if(global->class == 'S') 
			verify_value = 0.530770700573E-4;
		else if(global->class == 'W') 
			verify_value = 0.250391406439E-17; 
		else if(global->class == 'A') 
			verify_value = 0.2433365309E-5;
		else if(global->class == 'B') 
			verify_value = 0.180056440132E-5;
		else if(global->class == 'C') 
			verify_value = 0.570674826298E-6;
		if(fabs( rnm2 - verify_value ) < epsilon ) {
			verified = 1;
			printf(" Deviation is   %e\n", (rnm2 - verify_value));
		} else {
			verified = 0;
			printf(" The correct L2 Norm is %e\n", verify_value);
		}
	} else {
		verified = -1;
	}
	
	print_verification(global->class,verified,"BM");
	
	double mflops = 0.0;
	if( tm > 0.0 ) {
		mflops = 58.0*global->n_size*global->n_size*global->n_size;
		mflops *= global->n_it / (tm*1000000.0);
	}
	
	set_results(
		&res, 
		"MG",
		global->class,
		global->n_size,
		global->n_size,
		global->n_size,
		global->n_it,
		tm,
		mflops,
		(const char*)"floating point",
		verified,
		numThreads
	);

<<<<<<< Updated upstream
	print_results(&res, stdout);
}

//OLD MPI CODE

/*
	// REAL number;
	// REAL * tempa = (REAL*) malloc(sizeof(REAL)*10);
	REAL * tempdzyx = (REAL*) malloc(sizeof(REAL)*n1*n2*n3);
	//assign values to temp
	int i3, i2, i1;
	double ***z = u[0];
	//maps 3-D to 1-D array for transmission
	for(i3=0;i3<n3;i3++){
	        for(i2=0;i2<n2;i2++){
	            for(i1=0;i1<n1;i1++){
				tempdzyx[i3*n1*n2 + i2*n1 + i1] = z[i3][i2][i1];
			}
		}
	}
	//making the first element 0.1 because it will probably be 0.0- which makes it hard to tell if it actually transmitted correctly
	tempdzyx[0] = 0.1;
	REAL * tempa = (REAL*) malloc(sizeof(REAL)*10);
	for(i1=0;i1<10;i1++){
		tempa[i1] = tempdzyx[i1] + 1.0;
	}

	number = -1.5;
	// MPI_Send(&number, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
	// MPI_Send(tempa, 10, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
	MPI_Send(tempdzyx, n1*n2*n3, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
*/

/*
} else {
	// PPF_Print( MPI_COMM_WORLD, "Its the other processor!" );
	// REAL number;
	// REAL * tempa = (REAL*) malloc(sizeof(REAL)*10);
	REAL * tempdzyx = (REAL*) malloc(sizeof(REAL)*258*258*258);
	// MPI_Recv(&number, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// MPI_Recv(tempa, 10, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(tempdzyx, 258*258*258, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	PPF_Print( MPI_COMM_WORLD,"Process 1 received number %f from process 0\n", tempdzyx[0]);
	// PPF_Print( MPI_COMM_WORLD,"Processor receive succes\n");
}
/*
PPF_Print( MPI_COMM_WORLD, "Message from %N\n" );
PPF_Print( MPI_COMM_WORLD, (mpi_rank < (mpi_size/2) )
			? "Message from first half of nodes (%N)\n"
			: "Message from second half of nodes\n" );
PPF_Print( MPI_COMM_WORLD, (mpi_rank % 2)
			? "[%N] Message from odd numbered nodes\n"
			: "[%N] Message from even numbered nodes\n" );

*/
=======
	print_results(&res);
}
>>>>>>> Stashed changes
