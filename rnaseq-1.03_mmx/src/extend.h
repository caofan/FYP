#ifndef __EXTEND_GUARD__
#define __EXTEND_GUARD__

#include "file.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct Junction {
	unsigned p;
	unsigned q;
	unsigned r;
	char *Chrom;
	int Mismatches;
	int misL;
	char signal[5];
	bool isCanonical();
};
extern unsigned char* Original_Text;
Junction* extend(char* R, unsigned x, unsigned y, unsigned p, unsigned q, bool partion=false);
void loadPac(char* filename);
void Get_Bases_ASCII (unsigned Location,int StringLength,char* Org_String);
#endif
