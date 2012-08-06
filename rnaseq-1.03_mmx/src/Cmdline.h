#ifndef __CMDLINE_GUARD__
#define __CMDLINE_GUARD__

#include <getopt.h>
#include <stdlib.h>

class Parameters
{
	public:
	int MAX_TAGS_TO_PROCESS;	
	int Patternfile_Count;
	char* PATTERNFILE;	
	char* PATTERNFILE1;	
	char* JUNCTIONFILE;
	int EXONGAP;
	Parameters():MAX_TAGS_TO_PROCESS(0),Patternfile_Count(0),EXONGAP(20000),JUNCTIONFILE(NULL),PATTERNFILE(NULL),PATTERNFILE1(NULL)
	{
	}
};

struct Index_Info
{
	char* BWTFILE ; 
	char* OCCFILE ;
	char* SAFILE;
	char* REVBWTINDEX;
	char* REVOCCFILE;
	char* REVSAFILE;
	char* LOCATIONFILE;
	char* INDFILE;
	char* BLKFILE;
	char* BINFILE;
	char* SORTEDRANGEFILE;
	char* RANGEFILE;
	char* INFOFILE;
	unsigned char* Original_Text;//encoded genome file...
};

void Parse_Command_line(int argc, char* argv[],Index_Info & Ind,Parameters & CL);

#endif
