#include "Indexes.h"

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Load_Indexes 
 *  Description:  Opens FM index 
 * =====================================================================================
 */

BWT* Load_Indexes(char *BWTINDEX,char *OCCFILE, char *SAFILE, MMPool* & mmPool)
{
        int PoolSize = 524288;
	MMMasterInitialize(3, 0, FALSE, NULL);
	mmPool = MMPoolCreate(PoolSize);
	return BWTLoad(mmPool, BWTINDEX, OCCFILE, SAFILE, NULL, NULL, NULL);//Load FM index
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  UnLoad_Indexes 
 *  Description:  Opens FM index 
 * =====================================================================================
 */

void UnLoad_Indexes(BWT* revfmi,BWT *fwfmi,MMPool* mmPool,RANGEINDEX Range_Index)
{
	BWTFree(mmPool,fwfmi);
	BWTFree(mmPool,revfmi);
	MMPoolFree(mmPool);
	free(Range_Index.SA_Index);
	free(Range_Index.SA_Blocks);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Load_Range_Index
 *  Description:  Loads the range index sturcture..
 * =====================================================================================
 */
void Load_Range_Index(char* INDFILE, char* BLKFILE,RANGEINDEX & Range_Index)
{
	unsigned Index_Size,Block_Size;

	Range_Index.Index=File_Open(INDFILE,"rb");
	Range_Index.Blocks=File_Open(BLKFILE,"rb");
	Range_Index.SA_Index=(SA*)malloc((Index_Size=Get_File_Size(Range_Index.Index)));//contains the index to blocks...
	Range_Index.SA_Blocks=(int*)malloc((Block_Size=Get_File_Size(Range_Index.Blocks)));
	if (!Range_Index.SA_Index || !Range_Index.SA_Blocks)
	{
		printf ("Load_Range_Index(): Memory allocation failed!..\n");
		exit(0);
	}
	fread(&Range_Index.COMPRESS,1,1,Range_Index.Blocks);
	fread(&Range_Index.Hits,1,sizeof(Range_Index.Hits),Range_Index.Blocks);
	fread(Range_Index.SA_Index,Index_Size,1,Range_Index.Index);
	fread(Range_Index.SA_Blocks,Block_Size,1,Range_Index.Blocks);
}

