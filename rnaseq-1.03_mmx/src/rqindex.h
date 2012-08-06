#ifndef __RQINDEX_GUARD__
#define __RQINDEX_GUARD__
#include "common.h"
#define dist(x,y) (((x)>(y)) ? (x)-(y): (y)-(x))

extern int INDEX_RESOLUTION;
const int  SAGAP_CUTOFF=5;
const int  SAGAP_CUTOFF_INDEX=2;
const int SAGAP_CUTOFF_T=SAGAP_CUTOFF_INDEX, SAGAP_CUTOFF_H= SAGAP_CUTOFF;
const int MAX_JUNCS_TO_STORE=500;
const int MAX_HITS_TO_STORE=MAX_JUNCS_TO_STORE;//200;//later resolve this and MAXCOUNT conflict..
const int MAXCOUNT=MAX_HITS_TO_STORE; //Maximum number of pairs to enum...
struct PAIR
{
	unsigned Head;
	unsigned Tail;
	int HLevel;
	int TLevel;
};

struct TAG_INFO
{
	unsigned SA_Start;//start of the Sa range
	unsigned Gap;//length of SA range
	unsigned Block_Start;//Start of block info..
	unsigned Index;//Index to SA_index
	unsigned First;//first location of hit
	unsigned Last;//last location of hit
	unsigned Field_Length;
};

int Scan(MEMX & MF,int MAX_MISMATCHES, LEN & L,BWT* fwfmi,BWT* revfmi,int Next_Mis,int Max_Hits);
unsigned log2(unsigned v); // 32-bit word to find the log of
unsigned Get_Block_Start(unsigned SAValue,unsigned & M);
void Load_Info( TAG_INFO & Tag, SARange & Head);
inline unsigned Get_Location(TAG_INFO & Tag, unsigned Offset);
void Search_Small_Gap(SARange & Head, SARange & Tail, unsigned d,PAIR* Pairs,int & Pairs_Index);
void Get_Head_Tail(SARange & Head, SARange & Tail,unsigned d,PAIR* Pairs,int & Pairs_Index);
int Pair_Reads(SARange *Head_Hits,SARange *Tail_Hits,PAIR *Pairs,int & Pairs_Index);
#endif
