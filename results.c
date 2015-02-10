#include "results.h"

void init_results(results_t* res)
{
    memset(res, 0, sizeof(results_t));
    res->clss = 'S';
    res->optype = strdup("floating point");
}

void print_results(results_t* res, FILE *out)
{
    fprintf(out, "***** NAS Parallel Benchmarks C version %s ****\n", res->name);
    fprintf(out, "* Class             = %c\n", res->clss);

    if(res->n2 == 0 && res->n3 == 0 ){
        fprintf(out, "* Size              = %d\n", res->n1);
    }else{
        fprintf(out, "* Size              = %d X %d X %d\n", res->n1, res->n2, res->n3);
    }
    fprintf(out, "* Iterations        = %d\n", res->niter);
    fprintf(out, "* Time in seconds   = %6.4f\n", res->time);
    fprintf(out, "* ACCTime           = %e\n", res->acctime);
    fprintf(out, "* Mops total        = %6.4f\n", res->mops);
    fprintf(out, "* Operation type    = %s\n", res->optype);
    
    fprintf(out, "* Verification      = ");
    if (res->verified==1) fprintf(out, "Successful\n");
    else if(res->verified==0) fprintf(out, "Failed\n");
    else fprintf(out, "Not Performed\n");

    fprintf(out, "* Threads requested = %d\n", res->numthreads);
    
    fprintf(out, "* Please send all errors/feedbacks to:\n");
    fprintf(out, "* NPB Working Team\n");
    fprintf(out, "* npb@nas.nasa.gov\n");
    fprintf(out, "***************************************************************\n");
}

void print_verification(char clss, int verified, const char* BMName)
{
    if (clss == 'U'||verified == -1) {
        verified = -1;
        printf(" Problem size unknown\n");
        printf("%s.%c: Verification Not Performed\n", BMName, clss);
    }
    else if (verified==1) 
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
