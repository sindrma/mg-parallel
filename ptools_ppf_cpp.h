/* ptools_ppf.h */
/* John May -- 23 Oct 96 */
/* Definititions and declarations for the parallel print function */
/* 
 *	 Copyright (c) 1996-1998 The Regents of the Univerity of California
 *			and Lawrence Livermore National Laboratory
 *       All rights reserved.
 */

/* Modification history */
/*
 *	13 Jan 98	johnmay	Change names to follow Ptools conventions
 */

#ifndef _PTOOLS_PPF_H
#define _PTOOLS_PPF_H

#include <mpi.h>

#define PPF_MAXLINE 256

extern "C" {
int PPF_Print( MPI_Comm comm, const char * string, ... );
}

#endif /* _PTOOLS_PPF_H */

