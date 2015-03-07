#ifndef __RESULTS_H__
#define __RESULTS_H__

#include "../includes.h"

void print_results(results_t* res);
void init_results(results_t* res);
void set_results(results_t* res, const char* bname, char CLSS, int bn1, int bn2, int bn3, int bniter, double btime, double bmops, const char*boptype, int passed_verification, int num_threads);
void print_verification(char clss, int verified, const char* BMName);
void interpret_results(double rnm2, struct params* global, double tm, int numThreads);

#endif
