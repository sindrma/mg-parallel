#ifndef __RESULTS_H__
#define __RESULTS_H__

#include "includes.h"

//Represents the benchmark results
typedef struct results_s
{
    char* name;
    char* MachineName;
    char* PrLang;
    char  clss;
    int   n1, n2, n3, niter;
    double time,acctime,wctime,mops;
    double tmSent, tmReceived;
    int   RecArrSize;
    char* optype;
    int   numthreads;
    int   verified;
} results_t;

void print_results(results_t* res, FILE *out);
void init_results(results_t* res);

void set_results(results_t* res, const char* bname, char CLSS, int bn1, int bn2, int bn3, int bniter, double btime, double bmops, const char*boptype, int passed_verification, int num_threads);

void print_verification(char clss, int verified, const char* BMName);
#endif
