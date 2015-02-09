/*
---------------------------------------------------------------------
  Parameter lm (declared and set in "npbparams.h") is the log-base2 of 
  the edge size max for the partition on a given node, so must be changed 
  either to save space (if running a small case) or made bigger for larger 
  cases, for example, 512^3. Thus lm=7 means that the largest dimension 
  of a partition that can be solved on a node is 2^7 = 128. lm is set 
  automatically in npbparams.h
  Parameters ndim1, ndim2, ndim3 are the local problem dimensions. 
---------------------------------------------------------------------
*/
#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#include <stdbool.h>
#include "includes.h"



#endif
