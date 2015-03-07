#include "../types.h"
#include "../includes.h"

#define CLASS 'A'
#include <stdbool.h>

/*
This file is generated automatically by the setparams utility.
  It sets the number of processors and the class of the NPB
in this directory. Do not modify it by hand.
*/
void setup(int *n1, int *n2, int *n3, grid_t* grid);
struct params* setup_local(int argc, char **argv);
void set_a(double *a ,struct params* global);
void set_c(double *c ,struct params* global);
