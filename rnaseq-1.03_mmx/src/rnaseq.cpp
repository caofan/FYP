//BATMAN 1.1 - added Substring search..
//BATMAN 1.10 - enhanced default output processing...
//		handle N's.
//		print blanks...
//		--maxhits bug fixes
//{-----------------------------  INCLUDE FILES  -------------------------------------------------/
#include <stdio.h>
#include <string>
#include <limits.h>
#include <inttypes.h>
#include <string.h>
#include <stdlib.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include "bfix.h"
#include <getopt.h>
#include "zlib.h"
#include "const.h"
#include "assert.h"
#include "Hash.h"
#include "rqindex.h"
extern "C"
{
	#include "iniparser.h"
	#include <time.h>
	#include "MemManager.h"
	#include "MiscUtilities.h"
	#include "TextConverter.h"
	#include "BWT.h"
}
//}-----------------------------  INCLUDE FILES  -------------------------------------------------/

using namespace std;
#include "batlib.h"
#include "Cmdline.h"
#include "Indexes.h"
#include "file.h"
#include "extend.h"
//#include "batlibX.cpp"
//#include "rqindex.cpp"
//{-----------------------------  DEFINES  -------------------------------------------------/

//
//{---------------------------- GLOBAL VARIABLES -------------------------------------------------
int INDEX_RESOLUTION=30000;
int EXONGAP;
unsigned CONVERSION_FACTOR;
char* LOG_SUCCESS_FILE=NULL;
FILE *Log_SFile;
int MAX_MISMATCHES=1;
BWT *fwfmi,*revfmi;
unsigned Conversion_Factor;
map <unsigned, Ann_Info> Annotations;
unsigned Location_Array[80];
unsigned char* Original_Text;
char Char_To_CodeC[256];
char Char_To_Code[256];

unsigned RQ_Hits;
SA* SA_Index;
char COMPRESS;
int *SA_Blocks;
//}---------------------------- GLOBAL VARIABLES -------------------------------------------------


//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*
int Find_Single_Junc(char* Read,char* Converted_Read,MEMX & MF_Pre,MEMX & MF_Suf,int STRINGLENGTH,LEN & L,PAIR* & Pairs,Junction *Final_Juncs,int & Err);
void Load_All_Indexes(Index_Info Genome_Files,BWT* & fwfmi,BWT* & revfmi,MMPool* & mmPool,RANGEINDEX & Range_Index);
void Init(BWT *revfmi,unsigned & SOURCELENGTH,PAIR* & Pairs,gzFile & Input_File,gzFile & Mate_File,FILETYPE & File_Info,Parameters & CL,FILE* & OUT,Index_Info & Genome_Files);
bool  Progress_Bar(Parameters & CL,unsigned & Number_of_Tags,unsigned & Progress,unsigned & Tag_Count,FILETYPE & File_Info);
int Find_Pairings(int & Pairs_Index,SARange* MF_Pre_Hits,SARange* MF_Suf_Hits,PAIR* &  Pairs, char* Read,int STRINGLENGTH,int & Final_Juncs_Ptr,Junction *Final_Juncs,int & Min_Mismatch,int Mis_In_AnchorL,int Mis_In_AnchorR);
int Find_Two_Junc(char* Read,char* Converted_Read,MEMX & MF_Pre,MEMX & MF_Suf,MEMX & MF_Mid,int STRINGLENGTH,LEN & L,PAIR* & Pairs,Junction *Final_Juncs,int & Err);
//}-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------/*

int main(int argc, char* argv[])
{
	unsigned Total_Hits=0,Tags_Processed=0,Tag_Count=0;
	time_t Start_Time,End_Time;
	unsigned Number_of_Tags=1000;
	unsigned Progress=0;
	FILE* OUT;

//------------------- INIT -------------------------------------------------------------
	Index_Info Genome_Files;
	Parameters CL;
	MMPool *mmPool;
	RANGEINDEX Range_Index;
	unsigned SOURCELENGTH;
	PAIR *Pairs;
	Junction Final_Juncs[MAX_JUNCS_TO_STORE];
	gzFile Input_File,Mate_File;
	FILETYPE File_Info;

	Parse_Command_line(argc,argv,Genome_Files,CL);
	Load_All_Indexes(Genome_Files,fwfmi,revfmi,mmPool,Range_Index);
	Init(revfmi,SOURCELENGTH,Pairs,Input_File,Mate_File,File_Info,CL,OUT,Genome_Files);
//------------------- INIT -------------------------------------------------------------

	time(&Start_Time);
	READ Head,Tail;
	int LOOKUPSIZE=3;
	MEMLOOK MLook;MLook.Lookupsize=3;
	Build_Tables(fwfmi,revfmi,MLook);
	LEN L;L.IGNOREHEAD=0;
	Split_Read(RQFACTOR,L);//we are scanning 18-mers...
	Conversion_Factor=revfmi->textLength-RQFACTOR;
//--------------------- Setup Data Structure for Batman Prefix ----------------------------------------
	MEMX MF_Pre,MF_Suf,MF_Mid;//MemX is the data structure for doing Batman alignment. MF_Pre is for the prefix, MC for suffix..
	Init_Batman(MF_Pre,L,MLook,MAX_MISMATCHES);
	Init_Batman(MF_Suf,L,MLook,MAX_MISMATCHES);
	Init_Batman(MF_Mid,L,MLook,MAX_MISMATCHES);
//--------------------- Setup Data Structure for Batman End----------------------------------------

	time_t Maptime;
	{
		fprintf(stderr,"======================]\r[");//progress bar....
		while (Read_Tag(Head,Tail,Input_File,Mate_File,File_Info))
		{
			if(!Progress_Bar(CL,Number_of_Tags,Progress,Tag_Count,File_Info)) break;

			int Err;
			if(Find_Single_Junc(Head.Tag_Copy,Head.Tag,MF_Pre,MF_Suf,File_Info.STRINGLENGTH,L,Pairs,Final_Juncs,Err))
			{
				fprintf(OUT,"%s",Head.Description);
				for(int i=0;Final_Juncs[i].p!=INT_MAX;i++)
				{
					fprintf(OUT,"%s\t%d\t%d\t%d\t%d\n",Final_Juncs[i].Chrom,Final_Juncs[i].p,Final_Juncs[i].q,Final_Juncs[i].Mismatches,Tag_Count);
				}
			}
			else//Read is not mapped
			{
				Find_Two_Junc(Head.Tag_Copy,Head.Tag,MF_Pre,MF_Suf,MF_Mid,File_Info.STRINGLENGTH,L,Pairs,Final_Juncs,Err);

				if(!Err) printf(">%d-%s%s",Tag_Count,Head.Description,Head.Tag_Copy);
			}

		}
		UnLoad_Indexes(fwfmi,revfmi,mmPool,Range_Index);
		fprintf(stderr,"\r[++++++++100%%+++++++++]\n");//progress bar....
		time(&End_Time);
		Maptime=difftime(End_Time,Start_Time);
	}
	/*printf("Average %d\n",Sum);*/
	fprintf(stderr,"%u / %u Tags/Hits",Tag_Count,Total_Hits);
	time(&End_Time);fprintf(stderr,"\n Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
	fprintf(stderr,"\n Time Taken  - %.0lf Seconds ..\n ",Maptime);
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  Find_Two_Junc
 *  Description:  Find Two junctions in the read Read
 *  		  Returns least mismatch junctions in Final_Juncs, terminated by sentinel Final_Junc[Last].p=INT_MAX
 *  		  Err is set if there were overflows in hits,
 *  		  returns zero if no junctions found..
 * =====================================================================================
 */
int Find_Two_Junc(char* Read,char* Converted_Read,MEMX & MF_Pre,MEMX & MF_Suf,MEMX & MF_Mid,int STRINGLENGTH,LEN & L,PAIR* & Pairs,Junction *Final_Juncs,int & Err)
{
	int Pairs_Index=0;
	int Final_Juncs_Ptr=0;
	int Min_Mismatch=INT_MAX;
	Err=0;

	int preMis[2],sufMis[2];
	int preSplit[2], sufSplit[2];

	MF_Pre.Hits=0;MF_Pre.Hit_Array_Ptr=0;MF_Pre.Current_Tag=Converted_Read;MF_Pre.Hit_Array[0].Start=0;//setup read details to alignmentstructure..
	int Last_Mis_Suf=-1;
	int Last_Mis_Pre=Scan(MF_Pre,MAX_MISMATCHES,L,fwfmi,revfmi,0,UINT_MAX);
	if(Last_Mis_Pre < 0)
		return 0;
	preMis[0] = Last_Mis_Pre;
	preSplit[0] = 0; preSplit[1] = MF_Pre.Hit_Array_Ptr;


	MF_Suf.Hits=0;MF_Suf.Hit_Array_Ptr=0;MF_Suf.Current_Tag=Converted_Read+(STRINGLENGTH-RQFACTOR);MF_Suf.Hit_Array[0].Start=0;//setup read details to alignmentstructure..
	Last_Mis_Suf=Scan(MF_Suf,MAX_MISMATCHES,L,fwfmi,revfmi,0,UINT_MAX);
	if(Last_Mis_Suf < 0)
		return 0;
	sufMis[0] = Last_Mis_Suf;
	sufSplit[0] = 0; sufSplit[1] = MF_Suf.Hit_Array_Ptr;


	Last_Mis_Pre=Scan(MF_Pre,MAX_MISMATCHES,L,fwfmi,revfmi,Last_Mis_Pre+1,UINT_MAX);
	preMis[1] = Last_Mis_Pre;
	Last_Mis_Suf=Scan(MF_Suf,MAX_MISMATCHES,L,fwfmi,revfmi,Last_Mis_Suf+1,UINT_MAX);
	sufMis[1] = Last_Mis_Suf;

	for(int i = 0; i < 2 && preMis[i] >=0; i++)
	{
		for(int j = 0; j < 2 && sufMis[j] >=0; j++)
		{
			Err+=Pair_Reads(MF_Pre.Hit_Array+preSplit[i], MF_Suf.Hit_Array+sufSplit[j],Pairs,Pairs_Index);
			if(!Pairs_Index)
				return Err;
			while(--Pairs_Index >= 0)
			{
				Ann_Info A,A1;
				int L1=Pairs[Pairs_Index].Head;
				int L2=Pairs[Pairs_Index].Tail;

				Location_To_Genome(Pairs[Pairs_Index].Head,A);
				Location_To_Genome(Pairs[Pairs_Index].Tail,A1);
				if (Pairs[Pairs_Index].Head+STRINGLENGTH > A.Size||Pairs[Pairs_Index].Tail+STRINGLENGTH > A1.Size)//check for a Boundary Hit..
				{
					continue;
				}
				Junction* midPart = extend(Read,1,STRINGLENGTH-RQFACTOR,L1,L2,true);//partially extend the two anchors.
				for(int i=0;midPart[i].p!=INT_MAX;i++)
				{
					if(true)//check for mismatches
					{
						LEN midL; midL.IGNOREHEAD=0;
						int midLength = midPart[i].q - midPart[i].p + 1;
						char Converted_Mid_Tag[midLength+1];
						strncpy(Converted_Mid_Tag, &Converted_Read[midPart[i].p],midLength);
						Split_Read(midLength,midL);
						MF_Mid.Hits=0; MF_Mid.Hit_Array_Ptr=0;MF_Mid.Current_Tag=Converted_Mid_Tag; MF_Mid.Hit_Array[0].Start=0;
						int last_mid_mis = Scan(MF_Mid,MAX_MISMATCHES,midL,fwfmi,revfmi,0,UINT_MAX);//Scan for positions for the mid part;
						//Junction* leftJuncs = extend();
						//Junction* rightJuncs = extend();
					}
				}
			}
		}
	}


	if(Err)
		return 0;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  Find_Single_Junc
 *  Description:  Find Single junctions in the read Read
 *  		  Returns least mismatch junctions in Final_Juncs, terminated by sentinel Final_Junc[Last].p=INT_MAX
 *  		  Err is set if there were overflows in hits,
 *  		  returns zero if no junctions found..
 * =====================================================================================
 */
int Find_Single_Junc(char* Read,char* Converted_Read,MEMX & MF_Pre,MEMX & MF_Suf,int STRINGLENGTH,LEN & L,PAIR* & Pairs,Junction *Final_Juncs,int & Err)
{
	int Pairs_Index=0;
	int Final_Juncs_Ptr=0;
	int Ret_Value=0;//if Any err in pairing, gets set..
	int Min_Mismatch=INT_MAX;
	Err=0;
	//Get Hits for head..
	MF_Pre.Hits=0;MF_Pre.Hit_Array_Ptr=0;MF_Pre.Current_Tag=Converted_Read;MF_Pre.Hit_Array[0].Start=0;//setup read details to alignmentstructure..
	int Last_Mis_Suf= -1;
	int Last_Mis_Pre=Scan(MF_Pre,MAX_MISMATCHES,L,fwfmi,revfmi,0,UINT_MAX);
	int Top_Hit_End=MF_Pre.Hit_Array_Ptr;//Save the end point top hits..
	//------------------------------------ starts Lowest Mismatch Pairing of Suf/Pref -------------------------------------------------
	if(Last_Mis_Pre>=0)//Did we get a hit?
	{
		//Get Hits for suffix..
		MF_Suf.Hits=0;MF_Suf.Hit_Array_Ptr=0;MF_Suf.Current_Tag=Converted_Read+(STRINGLENGTH-RQFACTOR);MF_Suf.Hit_Array[0].Start=0;//setup read details to alignmentstructure..
		Last_Mis_Suf=Scan(MF_Suf,MAX_MISMATCHES,L,fwfmi,revfmi,0,UINT_MAX);
		int Top_Hit_End_Suf=MF_Suf.Hit_Array_Ptr;//Save the end point top hits..
	}
	if(Last_Mis_Pre >=0 && Last_Mis_Suf>=0)
	{
		Err+=Find_Pairings(Pairs_Index,MF_Pre.Hit_Array,MF_Suf.Hit_Array,Pairs,Read,STRINGLENGTH,Final_Juncs_Ptr,Final_Juncs,Min_Mismatch,Last_Mis_Pre,Last_Mis_Suf);
	}
	else return 0;//no anchors
	/*if (Final_Juncs_Ptr || Err)
	{
		Final_Juncs[Final_Juncs_Ptr].p=INT_MAX;
		return Final_Juncs_Ptr;
	}*/
	//------------------------------------ End Lowest Mismatch Pairing of Suf/Pref -------------------------------------------------
	//------------------------------------ Start Higest/Lowest Mismatch Pairing of Suf/Pref -------------------------------------------------
	SARange* Sub_Opt_Pre_Start=MF_Pre.Hit_Array+MF_Pre.Hit_Array_Ptr;
	SARange* Sub_Opt_Suf_Start=MF_Suf.Hit_Array+MF_Suf.Hit_Array_Ptr;
	int M1=Last_Mis_Pre,M2=Last_Mis_Suf;

	Last_Mis_Suf=Scan(MF_Suf,MAX_MISMATCHES,L,fwfmi,revfmi,Last_Mis_Suf+1,UINT_MAX);
	Last_Mis_Pre=Scan(MF_Pre,MAX_MISMATCHES,L,fwfmi,revfmi,Last_Mis_Pre+1,UINT_MAX);
	if(Last_Mis_Pre >=0)
	{
		Err+=Find_Pairings(Pairs_Index,Sub_Opt_Pre_Start,MF_Suf.Hit_Array,Pairs,Read,STRINGLENGTH,Final_Juncs_Ptr,Final_Juncs,Min_Mismatch,Last_Mis_Pre,M2);
	}
	if(Last_Mis_Suf >=0)
	{
		Err+=Find_Pairings(Pairs_Index,MF_Pre.Hit_Array,Sub_Opt_Suf_Start,Pairs,Read,STRINGLENGTH,Final_Juncs_Ptr,Final_Juncs,Min_Mismatch,M1,Last_Mis_Suf);
	}
	//------------------------------------ End Higest/Lowest Mismatch Pairing of Suf/Pref -------------------------------------------------
	//------------------------------------ End Higest Mismatch Pairing of Suf/Pref -------------------------------------------------
	if(Last_Mis_Pre >=0 && Last_Mis_Suf>=0)
	{
		Err+=Find_Pairings(Pairs_Index,Sub_Opt_Pre_Start,Sub_Opt_Suf_Start,Pairs,Read,STRINGLENGTH,Final_Juncs_Ptr,Final_Juncs,Min_Mismatch,Last_Mis_Pre,Last_Mis_Suf);
	}
	//------------------------------------ End Higest Mismatch Pairing of Suf/Pref -------------------------------------------------

	Final_Juncs[Final_Juncs_Ptr].p=INT_MAX;
	if (Err) return 0;
	return Final_Juncs_Ptr;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  Find_Pairings
 *  Description:  Find Junctions by pairing SA-Ranges in MF_Pre and MF_Suf..
 *  		  Returns least mismatch junctions, with mismatches not exceeding Min_Mismatch in Final_Juncs. Last junc at Final_Junc[Final_Juncs_Ptr-1]
 *  		  Err is set if there were overflows in hits, and returns zero..
 * =====================================================================================
 */
int Find_Pairings(int & Pairs_Index,SARange* MF_Pre_Hits,SARange* MF_Suf_Hits,PAIR* &  Pairs, char* Read,int STRINGLENGTH,int & Final_Juncs_Ptr,Junction *Final_Juncs,int & Min_Mismatch,int Mis_In_AnchorL, int Mis_In_AnchorR)
{

	int Ret_Value=0;//if Any err in pairing, gets set..
	Pairs_Index=0;
	int Err=Pair_Reads(MF_Pre_Hits,MF_Suf_Hits,Pairs,Pairs_Index);
	if(!Pairs_Index) return Err;
	while(--Pairs_Index>=0)
	{
		Ann_Info A,A1;
		int L1=Pairs[Pairs_Index].Head;
		int L2=Pairs[Pairs_Index].Tail;

		Location_To_Genome(Pairs[Pairs_Index].Head,A);
		Location_To_Genome(Pairs[Pairs_Index].Tail,A1);
		if (Pairs[Pairs_Index].Head+STRINGLENGTH > A.Size||Pairs[Pairs_Index].Tail+STRINGLENGTH > A1.Size)//check for a Boundary Hit..
		{
			continue;
		}
		if (A1.ID == A.ID)//in the same chromosome?
		{
			{
				Junction* junctions = extend(Read,1,STRINGLENGTH-RQFACTOR,L1,L2);
				if(junctions[0].Mismatches<=Min_Mismatch)
				{
					if(Min_Mismatch>junctions[0].Mismatches)
					{
						Min_Mismatch=junctions[0].Mismatches;Final_Juncs_Ptr=0;//Store min mismatch junctions...
					}
					for(int i=0;junctions[i].p!=INT_MAX;i++)
					{
						if(Final_Juncs_Ptr>=MAX_JUNCS_TO_STORE-1) {Err=1;break;}
						Location_To_Genome(junctions[i].p,A);
						Location_To_Genome(junctions[i].q,A); //A1?
						junctions[i].Chrom=A.Name;
						Final_Juncs[Final_Juncs_Ptr++]=junctions[i];
					}
				}
				delete [] junctions;
			}
		}
	}
	return Err;
}
/*int Find_Single_Junc(char* Read,char* Converted_Read,MEMX & MF,MEMX & MF_Suf,int STRINGLENGTH,LEN & L,PAIR* & Pairs,Junction *Final_Juncs,int & Err)
{
	int Pairs_Index=0;
	int Final_Juncs_Ptr=0;
	int Ret_Value=0;
	Err=0;
	//Get Hits for head..
	MF.Hits=0;MF.Hit_Array_Ptr=0;MF.Current_Tag=Converted_Read;MF.Hit_Array[0].Start=0;//setup read details to alignmentstructure..
	int Last_Mis_Suf= -1;
	int Last_Mis=Scan(MF,MAX_MISMATCHES,L,fwfmi,revfmi,0,UINT_MAX);
	int Top_Hit_End=MF.Hit_Array_Ptr;//Save the end point top hits..
	if(Last_Mis>=0)//Did we get a hit?
	{
		//Get Hits for suffix..
		MF_Suf.Hits=0;MF_Suf.Hit_Array_Ptr=0;MF_Suf.Current_Tag=Converted_Read+(STRINGLENGTH-RQFACTOR);MF_Suf.Hit_Array[0].Start=0;//setup read details to alignmentstructure..
		Last_Mis_Suf=Scan(MF_Suf,MAX_MISMATCHES,L,fwfmi,revfmi,0,UINT_MAX);
		int Top_Hit_End_Suf=MF_Suf.Hit_Array_Ptr;//Save the end point top hits..
	}
	int Suf_Last_Hits=0;//MF_Suf.Hit_Array_Ptr;
	int Pre_Last_Hits=0;//MF.Hit_Array_Ptr;
	int Pass=0;
	if(Last_Mis >=0 && Last_Mis_Suf>=0)
	{
		int Min_Mismatch=INT_MAX;
		do
		{
			Pairs_Index=0;
			//Err=Pair_Reads(MF.Hit_Array,MF_Suf.Hit_Array,Pairs,Pairs_Index);
			if(Last_Mis != -1) Err+=Pair_Reads(MF.Hit_Array+Pre_Last_Hits,MF_Suf.Hit_Array,Pairs,Pairs_Index);
			if (Pass++)//not the initial poass..
				Err+=Pair_Reads(MF.Hit_Array,MF_Suf.Hit_Array+Suf_Last_Hits,Pairs,Pairs_Index);
			//if(!Pairs_Index) return 0;
			while(--Pairs_Index>=0)
			{
				Ann_Info A,A1;
				int L1=Pairs[Pairs_Index].Head;
				int L2=Pairs[Pairs_Index].Tail;
				Junction* junctions = extend(Read,1,STRINGLENGTH-RQFACTOR,L1,L2);

				Location_To_Genome(Pairs[Pairs_Index].Head,A);
				Location_To_Genome(Pairs[Pairs_Index].Tail,A1);
				if (Pairs[Pairs_Index].Head+STRINGLENGTH > A.Size||Pairs[Pairs_Index].Tail+STRINGLENGTH > A.Size)
				{
					continue;//Boundary Hit..
				}
				if (A1.ID == A.ID)//in the same chromosome?
				{
					if(junctions[0].Mismatches<=Min_Mismatch)
					{
						if(Min_Mismatch>junctions[0].Mismatches)
						{
							Min_Mismatch=junctions[0].Mismatches;Final_Juncs_Ptr=0;//Store min mismatch junctions...
						}
						for(int i=0;junctions[i].p!=INT_MAX;i++)
						{
							if(Final_Juncs_Ptr>=MAX_JUNCS_TO_STORE-1) break;
							Location_To_Genome(junctions[i].p,A);
							Location_To_Genome(junctions[i].q,A);
							junctions[i].Chrom=A.Name;
							Final_Juncs[Final_Juncs_Ptr++]=junctions[i];
						}
					}
				}
				delete [] junctions;
			}
			Suf_Last_Hits=MF_Suf.Hit_Array_Ptr;
			Pre_Last_Hits=MF_Suf.Hit_Array_Ptr;
			if(Last_Mis_Suf != -1)
			{
				Last_Mis_Suf=Scan(MF_Suf,MAX_MISMATCHES,L,fwfmi,revfmi,Last_Mis_Suf+1,UINT_MAX);
			}
			if(Last_Mis != -1)
			{
				Last_Mis=Scan(MF,MAX_MISMATCHES,L,fwfmi,revfmi,Last_Mis+1,UINT_MAX);
			}
		}
		while (Last_Mis_Suf != -1 || Last_Mis != -1);
		Ret_Value+=Final_Juncs_Ptr;
	}
	Final_Juncs[Final_Juncs_Ptr].p=INT_MAX;
	return Ret_Value;
}*/

void Load_All_Indexes(Index_Info Genome_Files,BWT* & fwfmi,BWT* & revfmi,MMPool* & mmPool,RANGEINDEX & Range_Index)
{
	if(mkdir("Raw_Out",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) && errno != EEXIST) {fprintf(stderr,"Cannot create Temp directory..\n");exit(-1);};
	fprintf (stderr,"Loading index %s\n",Genome_Files.BWTFILE);
	fwfmi= Load_Indexes(Genome_Files.BWTFILE,Genome_Files.OCCFILE,Genome_Files.SAFILE,mmPool);
	fprintf (stderr,"Loading index %s\n",Genome_Files.REVBWTINDEX);
	revfmi= Load_Indexes(Genome_Files.REVBWTINDEX,Genome_Files.REVOCCFILE,Genome_Files.REVSAFILE,mmPool);
	fwfmi->saInterval=revfmi->saInterval;
	fprintf (stderr,"Loading Location %s\n",Genome_Files.LOCATIONFILE);
	int Genome_Count= Load_Location(Genome_Files.LOCATIONFILE,Annotations,Location_Array);
	fprintf (stderr,"Loading packed genome %s\n",Genome_Files.BINFILE);
	loadPac(Genome_Files.BINFILE);
	fprintf (stderr,"Loading index %s\n",Genome_Files.INDFILE);
	Load_Range_Index(Genome_Files.INDFILE,Genome_Files.BLKFILE,Range_Index);
	SA_Index=Range_Index.SA_Index;
	SA_Blocks=Range_Index.SA_Blocks;
	COMPRESS=Range_Index.COMPRESS;
	RQ_Hits=Range_Index.Hits;
	fprintf(stderr,"Done...\n");
}

void Init(BWT *revfmi,unsigned & SOURCELENGTH,PAIR* & Pairs,gzFile & Input_File,gzFile & Mate_File,FILETYPE & File_Info,Parameters & CL,FILE* & OUT,Index_Info & Genome_Files)
{
	EXONGAP=CL.EXONGAP;
	SOURCELENGTH = revfmi->textLength;
	CONVERSION_FACTOR=revfmi->textLength-RQFACTOR;//+1;
	Char_To_Code['N']=0;Char_To_Code['n']=0;
	Char_To_Code['A']=0;Char_To_Code['C']=1;Char_To_Code['G']=2;Char_To_Code['T']=3;
	Char_To_Code['a']=0;Char_To_Code['c']=1;Char_To_Code['g']=2;Char_To_Code['t']=3;
	Char_To_Code['+']='+';Char_To_Code['-']='-';//we are using character count to store the fmicode for acgt
	Char_To_CodeC['N']=3;Char_To_CodeC['n']=3;
	Char_To_CodeC[0]=3;Char_To_CodeC[1]=2;Char_To_CodeC[2]=1;Char_To_CodeC[3]=0;
	Char_To_CodeC['a']=3;Char_To_CodeC['c']=2;Char_To_CodeC['g']=1;Char_To_CodeC['t']=0;
	Char_To_CodeC['-']='-';Char_To_CodeC['+']='+';//we are using character count to store the fmicode for acgt

	if (!(Pairs=(PAIR*)malloc(sizeof(PAIR)*(MAX_HITS_TO_STORE+10)))) {fprintf(stderr,"Allocate_Memory():malloc error...\n");exit(100);}
	Open_Files(Input_File,Mate_File,CL);
	Detect_Input(File_Info,Input_File,Mate_File);
	if (CL.JUNCTIONFILE) OUT=File_Open(CL.JUNCTIONFILE,"w"); else OUT=stdout;
	FILE* Inf_File;
	if((Inf_File=File_Exist_Open(Genome_Files.INFOFILE)))
	{
		char String[40];
		int Num;
		fscanf(Inf_File,"%s%d",String,&Num);
		fscanf(Inf_File,"%s%d",String,&INDEX_RESOLUTION);
	}
	else
	{
		INDEX_RESOLUTION=30000;
	}
	//fclose(Inf_File);
}

bool  Progress_Bar(Parameters & CL,unsigned & Number_of_Tags,unsigned & Progress,unsigned & Tag_Count,FILETYPE & File_Info)
{
	Tag_Count++;
	Progress++;
	if (CL.MAX_TAGS_TO_PROCESS && Tag_Count >= CL.MAX_TAGS_TO_PROCESS) return false;

	if (Progress==Number_of_Tags)
	{
		if (CL.MAX_TAGS_TO_PROCESS)
		{
			Number_of_Tags=(CL.MAX_TAGS_TO_PROCESS)/20;
			Progress=0;
			Show_Progress(Tag_Count*100/CL.MAX_TAGS_TO_PROCESS);
		}
		else
		{
			off64_t Current_Pos=ftello64(File_Info.Org_File);
			unsigned Average_Length=Current_Pos/Tag_Count+1;//+1 avoids divide by zero..
			Number_of_Tags=(File_Info.File_Size/Average_Length)/20;
			Progress=0;
			Show_Progress(Current_Pos*100/File_Info.File_Size);
		}
	}
	return true;
}
