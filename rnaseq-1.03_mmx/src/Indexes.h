#ifndef __INDEXES_GUARD__
#define __INDEXES_GUARD__
#include "stdio.h"
#include "stdlib.h"
#include "file.h"
#include "batlib.h"
extern "C" 
{
        #include "iniparser.h"
        #include <time.h>
        #include "MemManager.h"
        #include "MiscUtilities.h"
        #include "TextConverter.h"
        #include "BWT.h"
}

struct SA
{

	unsigned Start;
	unsigned End;
	unsigned Start_Location;// exact location of the first occurance...
	unsigned End_Location;

};

struct RANGEINDEX
{
	FILE* Index;
	FILE* Blocks;
	SA *SA_Index;
	int *SA_Blocks;
	unsigned Hits;
	char COMPRESS;
};

BWT* Load_Indexes(char *BWTINDEX,char *OCCFILE, char *SAFILE, MMPool* & mmPool);
void Load_Range_Index(char* INDFILE, char* BLKFILE,RANGEINDEX & Range_Index);
void UnLoad_Indexes(BWT* revfmi,BWT *fwfmi,MMPool* mmPool,RANGEINDEX Range_Index);
#endif
