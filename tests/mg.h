#ifndef __MG__H_
#define __MG__H_

#include "../types.h"
#include "../includes.h"

void test_Flatten_AND_Unflatten(int x,int y,int z);

void testGhostCell(int x,int y,int z);
void testExchange(int x,int y,int z);


void testMerge(int x,int y,int z,int num_processors,bool buffered);
void testSplit(int x,int y,int z,int num_processors,bool buffered);

REAL *** generateMatrix(int x,int y,int z, bool buffered);

#endif
