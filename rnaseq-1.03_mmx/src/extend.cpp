#include "extend.h"
#include "limits.h"
#include "file.h"
using namespace std;


bool Junction::isCanonical(){
	return !(strcmp(signal, "GTAG") && strcmp(signal, "GCAG") && strcmp(signal, "ATAC") && strcmp(signal, "CTAC") && strcmp(signal, "CTGC") && strcmp(signal, "GTAT"));
}

/*int main(int argc, char* argv[]){
	loadPac("~/Documents/test.fasta.pac");
	//printf("%s\n",Original_Text);
	char* R = "ACGTTTTGCGTGAGTGCTGCTAGTGGTACGTGTGTGTACGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
	Junction* junctions = extend(R, 20, 120, 20,980,false);

	delete [] junctions;
	return 0;
}*/

void loadPac(char* filename){
	FILE* Original_File=File_Open(filename,"rb");
	Original_Text=(unsigned char*) malloc(Get_File_Size(Original_File));
	fread(Original_Text,Get_File_Size(Original_File),1,Original_File);
}

void Get_Bases_ASCII (unsigned Location,int StringLength,char* Org_String)
{
	for (int i=0;i<StringLength;i++)
	{
		unsigned char L= (unsigned char)(Original_Text[(Location+i)/4]<< (((Location+i) % 4) * 2)) >>6;
		Org_String[i]="ACGT"[L];
	}
}

/*
class Extension{
	private:
		char* R;
		unsigned x, y ,p ,q;
	public:
		Extension(char* new_R, unsigned new_x, unsigned new_y, unsigned new_p, unsigned new_q);
		Junction* extend();
		Junction* partialExtend();
}
*/

int findFirst(int array[], int v, int len, int c) {
	int l = 0, r = len-1;
	int mid = (l+r)/2;

	while(l+1!=r){
		if(c*array[mid] < c*v)
			l = mid;
		else
			r = mid;
		mid = (l+r)/2;
	}
	if(r>len || array[r]!=v)
		return -1;
	else
		return r;
}

int findLast(int array[], int v, int len, int c){
	int l = 0, r = len-1;
	int mid = (l+r)/2;

	while(l+1!=r){
		if(c*array[mid] > c*v)
			r = mid;
		else
			l = mid;
		mid = (l+r)/2;
	}
	if(r>len || array[l]!=v)
		return -1;
	else
		return l;
}

//Find full partitions
Junction* findFullPar(int misL[], int misR[], int size, unsigned x, unsigned p, unsigned q){
	//Find the partitions that give the minimum mismatches.
	//Partitions store the size of the portion of the left extension.
	int min= misL[0] + misR[0];
	int parCount = 1;
	int partitions[size+1];
	partitions[0] = 0;
	for(int i=1; i<size+1; i++)
	{
		//printf("%d\n",min);
		int temp = misL[i]+misR[i];
		if(temp<min)
		{
			parCount = 1;
			min = temp;
			partitions[0] = i;
		}
		else if(temp == min)
		{
			partitions[parCount] = i;
			++parCount;
		}
	}

	//Find the junctions based on the partitions obtained
	Junction *junctions = new Junction[parCount+1];
	printf("Min mismatches: %d\n",min);
	for(int i = 0; i<parCount; i++)
	{
		junctions[i].p = (unsigned)(p + x + partitions[i] + 1);
		junctions[i].q = (unsigned)(q - size + partitions[i] - 1);
		junctions[i].r = x + partitions[i] + 1;
		junctions[i].Mismatches = min;
		char donor[5], acceptor[3];
		donor[2] = 0;
		acceptor[2] = 0;
		Get_Bases_ASCII(junctions[i].p,2,donor);
		Get_Bases_ASCII(junctions[i].q-1,2,acceptor);
		strcat(donor,acceptor);
		//strcpy(donor,"GTAG");
		strcpy(junctions[i].signal, donor);
		//printf("%d\t%d\t%d\t%s\n",junctions[i].p,junctions[i].q,junctions[i].r,junctions[i].signal);
		//if(junctions[i].isCanonical())
		//	printf("canonical\n");
		//else
		//	printf("non-canonical\n");
	}
	junctions[parCount].p = INT_MAX;
	return junctions;
}

Junction* findPartialPar(int misL[], int misR[], int size, unsigned x, int maxMisL, int maxMisR){
	if(misL[size-1] < maxMisL)
		maxMisL = misL[size-1];
	if(misR[0] < maxMisR)
		maxMisR = misR[0];
	Junction *junctions = new Junction[(maxMisL+1)*(maxMisR+1)+1];
	int juncPtr = 0;
	for(int i=maxMisL; i>=0; i--)
	{
		int posL = findLast(misL,i,size+1,1);
		for(int j=maxMisR; j>=0; j--)
		{
			int posR = findFirst(misR,j,size+1,-1);
			if(posR-posL-1 > 18)
			{
				junctions[juncPtr].p = (unsigned)(x+posL+1);
				junctions[juncPtr].q = (unsigned)(x+posR);
				junctions[juncPtr].Mismatches = i+j;
				junctions[juncPtr].misL = i;
				juncPtr++;
			}
		}
	}
	junctions[juncPtr].p = INT_MAX;
	return junctions;
}

//Assume the coordinates are 0-based.
//All intervals used are closed interval.
Junction* extend(char* R, unsigned x, unsigned y, unsigned p, unsigned q, bool partial /*=false*/){
	unsigned size = y-x-1;
	int misL[size+1], misR[size+1];
	char basesL[size+1], basesR[size+1];
	basesL[size] = 0;
	basesR[size] = 0;

	Get_Bases_ASCII(p+x+1, size, basesL);
	Get_Bases_ASCII(q-size, size, basesR);

	//char str1[602];
	//str1[600]=0;
	//Get_Bases_ASCII(0, 600, str1);
	//printf("size: %d, %d, %s\n",size,q-size, str1);
	//printf("%c\t%d\t%s\n",basesL[0],strlen(basesL),basesL);
	//printf("%c\t%d\t%s\n",basesR[0],strlen(basesR),basesR);
	int countL=0, countR=0;
	misL[0] = misR[size] = 0;
	for(int i=1; i<size+1; i++)
	{
		if(R[x+i] != basesL[i-1])
			++countL;
		if(R[x+size+1-i] != basesR[size-i])
			++countR;
		//printf("%c\t%c\n",R[x+i], basesL[i-1]);
		misL[i] = countL;
		misR[size-i] = countR;
		//printf("misL: %d, %d\n",i,misL[i]);
		//printf("misR: %d, %d\n",i,misR[size-i]);
	}

	if(partial)
		return findPartialPar(misL, misR, size, x, 2, 2);//!Need to change these two to global variables.
	else
		return findFullPar(misL, misR, size, x, p, q);


}
