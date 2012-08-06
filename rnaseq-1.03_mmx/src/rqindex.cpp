#include "rqindex.h"
#include "Indexes.h"
#include "zlib.h"
#include "common.h"
#include "assert.h"
#include "batlib.h"
#include "bfix.h"
using namespace std;

extern unsigned CONVERSION_FACTOR;
extern int EXONGAP;

extern unsigned RQ_Hits;
extern SA* SA_Index;
extern char COMPRESS;
extern int *SA_Blocks;
extern BWT *revfmi;
extern unsigned Conversion_Factor;
//int Scan(MEMX & MF,MEMX & MC,int MAX_MISMATCHES, LEN & L,BWT* fwfmi,BWT* revfmi,int Next_Mis,int & Top,int Max_Hits)
int Scan(MEMX & MF,int MAX_MISMATCHES, LEN & L,BWT* fwfmi,BWT* revfmi,int Next_Mis,int Max_Hits)
{
	if(MAX_MISMATCHES < Next_Mis) return -1;
	assert(Next_Mis >=0);assert(MAX_MISMATCHES >= Next_Mis);assert (Next_Mis <= 5);
	int In_Mis=0,Hits=0;MF.Hits=0;//MC.Hits=0;
	if (Next_Mis == 0) goto Zero; else if (Next_Mis ==1) goto One;else if (Next_Mis ==2) goto Two;else if (Next_Mis ==3) goto Three;else if (Next_Mis ==4) goto Four; else goto Five;
Zero:
	Hits+=Zero_Mismatch(MF.Current_Tag,L,revfmi,MF);
	//Hits+=Zero_Mismatch(MC.Current_Tag,L,revfmi,MC);
One:
	if (!Hits && MAX_MISMATCHES >0)
	{
		In_Mis=1;
		Hits+=One_Mismatch(MF.Current_Tag,L,Max_Hits,fwfmi,revfmi,MF);
		//Hits+=One_Mismatch(MC.Current_Tag,L,Max_Hits,fwfmi,revfmi,MC);
	}
Two:
	if (!Hits && MAX_MISMATCHES >1)
	{
		In_Mis=2;
		Hits+=Two_Mismatch(MF.Current_Tag,L,Max_Hits,fwfmi,revfmi,MF);
		//Hits+=Two_Mismatch(MC.Current_Tag,L,Max_Hits,fwfmi,revfmi,MC);
	}
Three:
	if (!Hits && MAX_MISMATCHES >2)
	{
		In_Mis=3;
		Hits+=Three_Mismatch(MF.Current_Tag,L,Max_Hits,fwfmi,revfmi,MF);
		//Hits+=Three_Mismatch(MC.Current_Tag,L,Max_Hits,fwfmi,revfmi,MC);
	}
Four:
	if (!Hits && MAX_MISMATCHES >3)
	{
		In_Mis=4;
		Hits+=Four_Mismatch(MF.Current_Tag,L,Max_Hits,fwfmi,revfmi,MF);
		//Hits+=Four_Mismatch(MC.Current_Tag,L,Max_Hits,fwfmi,revfmi,MC);
	}
Five:
	if (!Hits && MAX_MISMATCHES >4)
	{
		In_Mis=5;
		Hits+=Five_Mismatch(MF.Current_Tag,L,Max_Hits,fwfmi,revfmi,MF);
		//Hits+=Five_Mismatch(MC.Current_Tag,L,Max_Hits,fwfmi,revfmi,MC);
	}

	MF.Hit_Array[MF.Hit_Array_Ptr].Start=0;//MC.Hit_Array[MC.Hit_Array_Ptr].Start=0;//tag sentinels to sa lists..
	MF.Hit_Array_Ptr++;//MC.Hit_Array_Ptr++;//Setup for suboptimal hits..
	assert(In_Mis <= MAX_MISMATCHES);
	//Top=Hits;
	return (Hits ? In_Mis : -1) ;
}

//{-----------------------------  Misc  -------------------------------------------------

static const char LogTable256[] = 
{

  1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
  
};

unsigned log2(unsigned v) // 32-bit word to find the log of
{
#ifdef BUILTIN_LOG
	return (32 -__builtin_clz(v));
#else
	unsigned r;     // r will be lg(v)
	register unsigned int t, tt; // temporaries

	if (tt = v >> 16)
	{
		r= (t = tt >> 8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	}
	else 
	{
		r= (t = v >> 8) ? 8 + LogTable256[t] : LogTable256[v];
	}
	return r;
#endif
}

//}-----------------------------  Misc  -------------------------------------------------


//{--------------------------------  Pair reads -------------------------------------------------------

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Get_Block_Start
 *  Description:  Finds the Block Value of the SAValue by bin search...
 *  		  Assumes SAValue exists...
 * =====================================================================================
 */
unsigned Get_Block_Start(unsigned SAValue,unsigned & M)
{
	unsigned L=0;
	unsigned H=RQ_Hits;

	while (L < H)
	{
		M=(L+H)/2;
		if (SA_Index[M].Start < SAValue)
		{
			L=M+1;
		}
		else
		{
			H=M;
		}
	}
	if (L==H) M=H; 

#ifdef DEBUG
	if (SA_Index[M].Start!=SAValue) 
	{
		printf("Get_Block_Start(): Bad SAValue !\n");
		exit(0);
	}
#endif

	return SA_Index[M].End;
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Load_Info
 *  Description:  Loads information necessary for reading coded info.... 
 * =====================================================================================
 */

void Load_Info( TAG_INFO & Tag, SARange & Head)
{
	Tag.Gap=Head.End-Head.Start+1;
	Tag.SA_Start=Head.Start;
	Tag.Block_Start=Get_Block_Start(Head.Start, Tag.Index);//get info block..
	Tag.First=SA_Index[Tag.Index].Start_Location;
	Tag.Last=SA_Index[Tag.Index].End_Location;
	Tag.Field_Length=log2(Tag.Gap);
}


inline unsigned Get_Location(TAG_INFO & Tag, unsigned Offset)
{
	unsigned SAPos;
	if(0==Offset) return Tag.First;
	if (Offset==Tag.Gap-1) return Tag.Last;
	Offset--;
	if(COMPRESS)
	{
		SAPos=Tag.SA_Start + bfx((unsigned char*)SA_Blocks,Tag.Block_Start+(Offset*Tag.Field_Length),Tag.Field_Length);
		return CONVERSION_FACTOR-BWTSaValue(revfmi,SAPos);
	}
	else
	{
		return (unsigned)SA_Blocks[Tag.Block_Start+Offset];
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Search_Small_Gap
 *  Description:  Do the searching when at least one of the pairs have small SA range... 
 * 		  Store the result in Pairs, starting from Pairs[Pairs_Index]
 * 		  Terminates hit list when encountering Pairs[x].Head=0... 
 * =====================================================================================
 */
void Search_Small_Gap(SARange & Head, SARange & Tail, unsigned d,PAIR* Pairs,int & Pairs_Index)
{
	unsigned H1,H2,T1,T2,L,H,M;
	TAG_INFO Head_Info,Tail_Info;
	assert(Pairs_Index<MAX_HITS_TO_STORE);


	if(Head.Start==Head.End)//Head is unique...
	{

		H1=Head.Start;//Conversion_Factor-BWTSaValue(revfmi,Head.Start);
		if (Tail.Start==Tail.End)// both hits are unique...
		{
			T2=Tail.Start;//Conversion_Factor-BWTSaValue(revfmi,Tail.Start);
			if (T2>H1 && T2 < H1+d)//modify for multiple hits...
			{
				Pairs[Pairs_Index].Head=H1;
				Pairs[Pairs_Index].HLevel=Head.Level;
				Pairs[Pairs_Index].Tail=T2;
				Pairs[Pairs_Index].TLevel=Tail.Level;
				Pairs_Index++;
				if (Pairs_Index == MAX_HITS_TO_STORE){Pairs[Pairs_Index].Head=0;return;}
#ifdef DEBUG
				if(H1>T2) 
				{
					printf("Search_Small_Gap(6):Enum error...\n");
					exit(0);
				}
#endif
			}
			Pairs[Pairs_Index].Head=0;
			return;
		}
		else //tail has multiples...
		{
			if (Tail.End-Tail.Start > SAGAP_CUTOFF)
			{
				Load_Info(Tail_Info,Tail);
				T1=Tail_Info.First;T2=Tail_Info.Last;


				if(H1<T1 )//Possible case for T1>H1 
				{
					if(H1+d>T1)
					{
						M=0;
						while (H1<T1 && H1+d>T1)//enumerate hits...
						{
							Pairs[Pairs_Index].Head=H1;
							Pairs[Pairs_Index].HLevel=Head.Level;
							Pairs[Pairs_Index].Tail=T1;
							Pairs[Pairs_Index].TLevel=Tail.Level;
							Pairs_Index++;
							if (Pairs_Index == MAX_HITS_TO_STORE){Pairs[Pairs_Index].Head=0;return;}
#ifdef DEBUG
						if(H1>T1) 
						{
							printf("Search_Small_Gap(5):Enum error...\n");
							exit(0);
						}
#endif
							if (Pairs_Index>MAXCOUNT) 
							{
								Pairs[Pairs_Index].Head=0;
								return;
							}
							M++;
							T1=Get_Location(Tail_Info,M);
							if(M>=Tail_Info.Gap) break;
						}
					}
				}
				else if(T2>H1) //H1 inside tail gaps..
				{

					L=0;
					H=Tail_Info.Gap;
					while (L < H)
					{
						M=(L+H)/2;
						if (Get_Location(Tail_Info,M) > H1)
						{
							H=M;
						}
						else
						{
							L=M+1;
						}
					}
					if (L==H) M=H;//find tail position closest to unique head...

					T1=Get_Location(Tail_Info,M);
					while (H1<T1 && H1+d>T1)//enumerate hits...
					{
						Pairs[Pairs_Index].Head=H1;
						Pairs[Pairs_Index].HLevel=Head.Level;
						Pairs[Pairs_Index].Tail=T1;
						Pairs[Pairs_Index].TLevel=Tail.Level;
						Pairs_Index++;
						if (Pairs_Index == MAX_HITS_TO_STORE){Pairs[Pairs_Index].Head=0;return;}
#ifdef DEBUG
						if(H1>T1) 
						{
							printf("Search_Small_Gap(4):Enum error...\n");
							exit(0);
						}
#endif
						if (Pairs_Index>MAXCOUNT) 
						{
							Pairs[Pairs_Index].Head=0;
							return;
						}
						M++;
						T1=Get_Location(Tail_Info,M);
						if(M>=Tail_Info.Gap) break;
					}
				}

				Pairs[Pairs_Index].Head=0;
				return;
			}
			else//Unique head and tail with gap below cutoff...
			{
				for (unsigned i=Tail.Start;i<=Tail.End;i++)
				{
					//unsigned Hit=(Tail.Conversion_Factor-BWTSaValue(revfmi,i));
					unsigned Hit=(Conversion_Factor-BWTSaValue(revfmi,i));
					if(Hit > H1 && Hit < H1+d)
					{
						Pairs[Pairs_Index].Head=H1;
						Pairs[Pairs_Index].HLevel=Head.Level;
						Pairs[Pairs_Index].Tail=Hit;
						Pairs[Pairs_Index].TLevel=Tail.Level;
						Pairs_Index++;
						if (Pairs_Index == MAX_HITS_TO_STORE){Pairs[Pairs_Index].Head=0;return;}
#ifdef DEBUG
						if(H1>Hit) 
						{
							printf("Search_Small_Gap(3):Enum error...\n");
							exit(0);
						}
#endif

						if(Pairs_Index >MAXCOUNT) break;
					}
				}
				Pairs[Pairs_Index].Head=0;
				return;
			}
		}
	}
	else //if(Tail.End==Tail.Start)//Tail is unique...
	{
		T1=Tail.Start;//Conversion_Factor-BWTSaValue(revfmi,Tail.Start);
		if(Head.End-Head.Start>SAGAP_CUTOFF)//Unique tail, but with multiple possible heads...
		{
			Load_Info(Head_Info,Head);
			H1=Head_Info.First;H2=Head_Info.Last;
			if(T1 > H1) //Head should not be after T1...
			{
				if(H2+d >T1)//Tail not too far away from heads...
				{
					//T1-d is between H1 and H2, search for the closest hit...
					L=0;
					H=Head_Info.Gap-1;
					T1-=d;
					while (L < H)
					{
						M=(L+H)/2;
						if (Get_Location(Head_Info,M) < T1)
						{
							L=M+1;
						}
						else
						{
							H=M;
						}
					}
					if (L==H) M=L;
					H1=Get_Location(Head_Info,M);

					T1+=d;
					while (H1<T1 && H1+d>T1)//enumerate hits...
					{
						Pairs[Pairs_Index].Head=H1;
						Pairs[Pairs_Index].HLevel=Head.Level;
						Pairs[Pairs_Index].Tail=T1;
						Pairs[Pairs_Index].TLevel=Tail.Level;
						Pairs_Index++;
						if (Pairs_Index == MAX_HITS_TO_STORE){Pairs[Pairs_Index].Head=0;return;}
#ifdef DEBUG
						if(H1>T1) 
						{
							printf("Search_Small_Gap(2):Enum error...\n");
							exit(0);
						}
#endif
						if (Pairs_Index>MAXCOUNT) 
						{
							break;	
						}
						M++;
						H1=Get_Location(Head_Info,M);
						if(M>=Head_Info.Gap) break;
					}

				}
			}
			Pairs[Pairs_Index].Head=0;
			return;
		}
		else//Unique tail and head with gap below cutoff...
		{
			for (unsigned i=Head.Start;i<=Head.End;i++)//try all hits...
			{
				//unsigned Hit=(Head.Conversion_Factor-BWTSaValue(revfmi,i));
				unsigned Hit=(Conversion_Factor-BWTSaValue(revfmi,i));
				if(Hit < T1 && T1 < Hit+d)
				{
					Pairs[Pairs_Index].Head=Hit;
					Pairs[Pairs_Index].HLevel=Head.Level;
					Pairs[Pairs_Index].Tail=T1;
					Pairs[Pairs_Index].TLevel=Tail.Level;
					Pairs_Index++;
					if (Pairs_Index == MAX_HITS_TO_STORE){Pairs[Pairs_Index].Head=0;return;}
#ifdef DEBUG
				if(Hit>T1) 
				{
					printf("Search_Small_Gap(1):Enum error...\n");
					exit(0);
				}
#endif
					if(Pairs_Index >MAXCOUNT) break;
				}
			}
			Pairs[Pairs_Index].Head=0;
			return;
		}
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Get_Head_Tail...
 *  Description:  find the possible pairs for head and tail that are a distance d apart.. 
 *  		  Pairs are stored in the structure Pair..
 *  		  Pairs_Index is the offset of the sentinel indicating last hit...
 * =====================================================================================
 */

void Get_Head_Tail(SARange & Head, SARange & Tail,unsigned d,PAIR* Pairs,int & Pairs_Index)
{
	assert(Pairs_Index<MAX_HITS_TO_STORE);
	unsigned H1,H2,T1,T2;
	unsigned H,M,L;
	unsigned Location=0,Hd,Tl;

	TAG_INFO Head_Info,Tail_Info;


	Load_Info(Head_Info,Head);
	Load_Info(Tail_Info,Tail);
	H1=Head_Info.First;H2=Head_Info.Last + d;
	T1=Tail_Info.First;T2=Tail_Info.Last;


	if (T1>H2 || T2<H1) //disjoint head and pair...
	{
		Pairs[Pairs_Index].Head=0;//cannot be paired...
		return;
	}

	if (T1>H1) M=0;//Get in M the index of the closest value of Tail larger than H1
	else
	{
		L=0;
		H=Tail_Info.Gap-1;
		while (L < H)
		{
			M=(L+H)/2;
			if (Get_Location(Tail_Info,M) < H1)
			{
				L=M+1;
			}
			else
			{
				H=M;
			}
		}
		if (L==H) M=H;
	}


	Hd=H1;H=1;
	Tl=Get_Location(Tail_Info,M);M++;

#ifdef DEBUG
	if(Tl<Hd) {printf("Get_Head_Tail(): bin search error..\n");exit(0);}
#endif

	int MHead=1;//Array pos. of head...
	int TempM;
	unsigned TempTl;

	while (Hd<T2)//Tail must always be after Head...
	{
		if (Hd<Tl)
		{
			TempM=M;TempTl=Tl;
			while (Hd+d>Tl)//enum. Hits...
			{
				Pairs[Pairs_Index].Head=Hd;
				Pairs[Pairs_Index].HLevel=Head.Level;
				Pairs[Pairs_Index].Tail=Tl;
				Pairs[Pairs_Index].TLevel=Tail.Level;
				Pairs_Index++;
				if (Pairs_Index == MAX_HITS_TO_STORE){Pairs[Pairs_Index].Head=0;return;}
#ifdef DEBUG
				if(Hd>Tl) 
				{
					printf("Get_Head_Tail():Enum error...\n");
					exit(0);
				}
#endif
				if(M >= Tail_Info.Gap) break;else Tl=Get_Location(Tail_Info,M);
				M++;
			}
			M=TempM;Tl=TempTl;

			if (MHead>=Head_Info.Gap) break; else Hd=Get_Location(Head_Info,MHead);
			MHead++;
		}
		else
		{
			if(M >= Tail_Info.Gap) break;else Tl=Get_Location(Tail_Info,M);
			M++;
		}
	}
	Pairs[Pairs_Index].Head=0;
}


//Will not execute if Pairs found is >500
int Pair_Reads(SARange *Head_Hits,SARange *Tail_Hits,PAIR *Pairs,int & Pairs_Index)
{
	if (Pairs_Index >=MAX_HITS_TO_STORE) {return 1;}
	SARange Head,Tail;
	int TGap,HGap;
	struct Valid_Hed
	{
		unsigned Location;
		bool Valid;
	}
	Multi_Head[SAGAP_CUTOFF+1];
	static int Sum=0;
	int Err=0;

	for(int i=0;Head_Hits[i].Start;i++)//Iterate Head Hits
	{
		for(int j=0;Tail_Hits[j].Start;j++)//With Tail Hits
		{
			Head=Head_Hits[i];
			Tail=Tail_Hits[j];
			TGap=Tail.End-Tail.Start;
			HGap=Head.End-Head.Start;

			if(TGap>INDEX_RESOLUTION || HGap > INDEX_RESOLUTION)
			{
				Err=1;
				continue;
			}
			if (TGap<= SAGAP_CUTOFF_T || HGap <= SAGAP_CUTOFF_H)//Small sa ranges
			{
				if(TGap && HGap)//Head and tail multi hits...
				{
					unsigned Start;
					if(HGap<TGap)// || HGap > SAGAP_CUTOFF_INDEX)///debug=1
					{
						Start=Head.Start;
						for (int k=0;k<=HGap;k++)
						{
							Head.Start=Start;
							//unsigned Loc=Head.End=Head.Start=Head.Conversion_Factor-BWTSaValue(revfmi,Head.Start);
							unsigned Loc=Head.End=Head.Start=Conversion_Factor-BWTSaValue(revfmi,Head.Start);
							//Search_Small_Gap(Head,Tail,EXONGAP,Pairs,Pairs_Index);
							Multi_Head[k].Location=Loc;Multi_Head[k].Valid=true;
							/*for (int j=0;j<k;j++)
							{
								if (dist(Multi_Head[j].Location,Head.Start) <400000)
								{
									Multi_Head[j].Valid=false;Multi_Head[k].Valid=false;
								}
							}*/
							Start++;
						}
						for (int k=0;k<=HGap;k++)
						{
							if (Multi_Head[k].Valid)
							{
								Head.Start=Head.End=Multi_Head[k].Location;
								Search_Small_Gap(Head,Tail,EXONGAP,Pairs,Pairs_Index);
							}
						}
					}
					else
					{
						Start=Tail.Start;
						for (int k=0;k<=TGap;k++)
						{
							Tail.Start=Start;
							//unsigned Loc=Tail.End=Tail.Start=Tail.Conversion_Factor-BWTSaValue(revfmi,Tail.Start);
							unsigned Loc=Tail.End=Tail.Start=Conversion_Factor-BWTSaValue(revfmi,Tail.Start);
							//Search_Small_Gap(Head,Tail,EXONGAP,Pairs,Pairs_Index);
							Multi_Head[k].Location=Loc;Multi_Head[k].Valid=true;
							/*for (int j=0;j<k;j++)
							{
								if (dist(Multi_Head[j].Location,Tail.Start) <400000)
								{
									Multi_Head[j].Valid=false;Multi_Head[k].Valid=false;
								}
							}*/
							Start++;
						}
						for (int k=0;k<=TGap;k++)
						{
							if (Multi_Head[k].Valid)
							{
								Tail.Start=Tail.End=Multi_Head[k].Location;
								Search_Small_Gap(Head,Tail,EXONGAP,Pairs,Pairs_Index);
							}
						}
					}
				}
				else//Head or Tail unique..
				{
					Search_Small_Gap(Head,Tail,EXONGAP,Pairs,Pairs_Index);
				}
			}
			else//Two big SA gaps...
			{
				Sum++;int T=1;
				//printf("%d:",Sum);
				//if (HGap > 500) {printf("HGAP:%d\n",HGap);continue;}
//Sum++;T=1;
				Get_Head_Tail(Head, Tail,EXONGAP,Pairs,Pairs_Index);
			}
			if(Pairs_Index>=MAX_HITS_TO_STORE-1) {Err=1;break;}
			assert(Pairs_Index<=MAX_HITS_TO_STORE);
		}
		if(Pairs_Index>=MAX_HITS_TO_STORE-1) {Err=1;break;}
		assert(Pairs_Index<=MAX_HITS_TO_STORE);
	}
	return Err;

}

//}--------------------------------  Pair reads -------------------------------------------------------
