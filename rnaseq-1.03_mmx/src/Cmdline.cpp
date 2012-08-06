#include <Cmdline.h>
#include "stdlib.h"
#include <stdio.h> 
#include "string.h"
#include "const.h"
#define FALSE 0
option Long_Options[]=
{
{"help",0,NULL,'h'},
{"query",1,NULL,'q'},
{0,0,0,0}
};


void Parse_Command_line(int argc, char* argv[],Index_Info & Ind,Parameters & CL)
{
	int Current_Option=0;
	char Short_Options[] ="jM:bhq:t:g:G:n:N:o:w:mprO::";//allowed options....
	char* This_Program = argv[0];//Current program name....
	char Help_String[]=
"Parameters:\n"
" --help | -h\t\t\t\t Print help\n"
" --query | -m \t\t map only\n"
" --query | -M \t\t Store <number> junctions..\n"
" --query | -p \t\t process only\n"
" --query | -q <filename>\t\t Query file(File of Tags)\n"
" --query | -r \t\t mark exons in refGene..\n"
" --query | -o <filename>\t\t junction file...\n"
" --query | -w <filename>\t\t Wiggle file...\n"
" --query | -O <filename>\t\t log read mapping info...\n"
" --query | -G <number>\t\t maximum gap between two introns...\n"
" --query | -n <number>\t\t number of mismatches in non spliced junctions...\n"
" --query | -N <number>\t\t number of mismatches in splices...\n"
" --query | -b \t\t build files from refgene\n"
" --query | -j \t\t output all junctions...\n"
;

	if(argc == 1) {printf("%s \n",Help_String);exit(0);}
	char *Source=(char*)malloc(sizeof(char)*6500);//create space for file names...
	char *options, *value; 
	char* Name;int Last_Dash;char* Genome_Name;
	char MARKEX=FALSE;

	for(;;)	
	{
		Current_Option=getopt_long(argc, argv, Short_Options, Long_Options, NULL);
		if (Current_Option == -1 ) break;
		switch(Current_Option)
		{
			case 'h':
				printf("%s \n",Help_String);exit(0);
			case 't':
				CL.MAX_TAGS_TO_PROCESS=atoi(optarg);
				break;
			case 'q':
				if(!CL.Patternfile_Count){CL.PATTERNFILE=optarg;}
				else CL.PATTERNFILE1=optarg;
				CL.Patternfile_Count++;
				break;
			case 'G':
				CL.EXONGAP=atoi(optarg);
				break;
			case 'o':
				CL.JUNCTIONFILE=optarg;
				break;
			/*case 'j':
				DUMP_ALL_JUNC=TRUE;
				break;
			case'b':
				MARKEX=TRUE;
				break;
			case 'm':
				MAPMODE=TRUE;PROCESSMODE=FALSE;
				break;
			case 'M':
				MAX_HITS_TO_STORE=atoi(optarg);
				break;
			case 'p':
				MAPMODE=FALSE;PROCESSMODE=TRUE;
				break;
			case 'n':
				MIS_IN_INITMAP=atoi(optarg);
				break;
			case 'N':
				COUNT=atoi(optarg);
				break;
			case 'r':
				USEREFGENE=TRUE;
				break;
			case 'w':
				WIGGLEFILE=optarg;
				break;
			case 'O':
				WRITE_SPLITREAD=TRUE;
				if (optarg) MAPFILE=optarg;
				break;*/
			case 'g':
				Name=optarg;Last_Dash=0;Genome_Name=optarg;
				for(;Name[0]!=0;Name++)
				{
					if (Name[0]=='/') 
					{
						Last_Dash++;Genome_Name=Name;
					}
				}

				Ind.REVBWTINDEX = (char*)Source;
				if(Last_Dash) Last_Dash=Genome_Name-optarg+1; else Genome_Name--;
				strncpy(Ind.REVBWTINDEX,optarg,Last_Dash);
				Ind.REVBWTINDEX[Last_Dash+0]='r';Ind.REVBWTINDEX[Last_Dash+1]='e';Ind.REVBWTINDEX[Last_Dash+2]='v';
				strcpy(Ind.REVBWTINDEX+Last_Dash+3,Genome_Name+1);
				strcat(Ind.REVBWTINDEX+Last_Dash+3,".bwt"); 

				Ind.BWTFILE=Ind.REVBWTINDEX+500;
				strncpy(Ind.BWTFILE,optarg,Last_Dash);
				strcpy(Ind.BWTFILE+Last_Dash,Genome_Name+1);
				strcat(Ind.BWTFILE+Last_Dash,".bwt"); 


				Ind.REVOCCFILE = Ind.BWTFILE+500;
				strncpy(Ind.REVOCCFILE,optarg,Last_Dash);
				Ind.REVOCCFILE[Last_Dash+0]='r';Ind.REVOCCFILE[Last_Dash+1]='e';Ind.REVOCCFILE[Last_Dash+2]='v';
				strcpy(Ind.REVOCCFILE+Last_Dash+3,Genome_Name+1);
				strcat(Ind.REVOCCFILE+Last_Dash+3,".fmv"); 


				Ind.OCCFILE=Ind.REVOCCFILE+500;			
				strncpy(Ind.OCCFILE,optarg,Last_Dash);
				strcpy(Ind.OCCFILE+Last_Dash,Genome_Name+1);
				strcat(Ind.OCCFILE+Last_Dash,".fmv"); 

				Ind.SAFILE=Ind.OCCFILE+500;			
				strncpy(Ind.SAFILE,optarg,Last_Dash);
				strcpy(Ind.SAFILE+Last_Dash,Genome_Name+1);
				strcat(Ind.SAFILE+Last_Dash,".sa");

				Ind.REVSAFILE = Ind.SAFILE+500;
				strncpy(Ind.REVSAFILE,optarg,Last_Dash);
				Ind.REVSAFILE[Last_Dash+0]='r';Ind.REVSAFILE[Last_Dash+1]='e';Ind.REVSAFILE[Last_Dash+2]='v';
				strcpy(Ind.REVSAFILE+Last_Dash+3,Genome_Name+1);
				strcat(Ind.REVSAFILE+Last_Dash+3,".sa"); 

				Ind.BINFILE=Ind.REVSAFILE+500;			
				strncpy(Ind.BINFILE,optarg,Last_Dash);
				strcpy(Ind.BINFILE+Last_Dash,Genome_Name+1);
				strcat(Ind.BINFILE+Last_Dash,".pac");

				Ind.LOCATIONFILE=Ind.BINFILE+500;			
				strncpy(Ind.LOCATIONFILE,optarg,Last_Dash);
				strcpy(Ind.LOCATIONFILE+Last_Dash,Genome_Name+1);
				strcat(Ind.LOCATIONFILE+Last_Dash,".ann.location");

                                Ind.BLKFILE = Ind.LOCATIONFILE+500;
                                strncpy(Ind.BLKFILE,optarg,Last_Dash);
                                strcpy(Ind.BLKFILE+Last_Dash,Genome_Name+1);
                                strcat(Ind.BLKFILE+Last_Dash,".blk.");
				sprintf(Ind.BLKFILE+strlen(Ind.BLKFILE),"%d",RQFACTOR);

                                Ind.INDFILE = Ind.BLKFILE+500;
                                strncpy(Ind.INDFILE,optarg,Last_Dash);
                                strcpy(Ind.INDFILE+Last_Dash,Genome_Name+1);
                                strcat(Ind.INDFILE+Last_Dash,".ind.");
				sprintf(Ind.INDFILE+strlen(Ind.INDFILE),"%d",RQFACTOR);

                                Ind.RANGEFILE = Ind.INDFILE+500;
                                strncpy(Ind.RANGEFILE,optarg,Last_Dash);
                                strcpy(Ind.RANGEFILE+Last_Dash,Genome_Name+1);
                                strcat(Ind.RANGEFILE+Last_Dash,".range");

                                Ind.SORTEDRANGEFILE = Ind.RANGEFILE+500;
                                strncpy(Ind.SORTEDRANGEFILE,optarg,Last_Dash);
                                strcpy(Ind.SORTEDRANGEFILE+Last_Dash,Genome_Name+1);
                                strcat(Ind.SORTEDRANGEFILE+Last_Dash,".sort");

                                Ind.INFOFILE = Ind.SORTEDRANGEFILE+500;
                                strncpy(Ind.INFOFILE,optarg,Last_Dash);
                                strcpy(Ind.INFOFILE+Last_Dash,Genome_Name+1);
                                strcat(Ind.INFOFILE+Last_Dash,".info");
				break;
			default:
				printf("%s \n",Help_String);
				exit(0);
		}
	}	
	CL.Patternfile_Count--;
	/*if (MARKEX) {Mark_Exons();exit(0);}*/
}
