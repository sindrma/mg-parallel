#ifndef __RANDOM_H__
#define __RANDOM_H__

double rnd_randlc(double x, double a);
//double rnd_randlc(double a);
double rnd_vranlc(double n, double x, double a, double *y, int offset);
double rnd_ipow46(double a, int exponent );
double rnd_power( double a, int n );

#endif
