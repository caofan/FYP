/**
 * @author 
 * @version 2012-05-16
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
using namespace std;

struct Junc {
	unsigned p;
	unsigned q;
	unsigned r;
	char signal[5];
	bool isCanonical();
};
unsigned char* Original_Text;

void loadPac(char* filename);
Junc* extend(char* R, unsigned x, unsigned y, unsigned p, unsigned q);
void Get_Bases (unsigned Location,int StringLength,char* Org_String);
unsigned Get_File_Size(FILE* File);
FILE* File_Open(const char* File_Name,const char* Mode);

bool Junc::isCanonical(){
	return !(strcmp(signal, "GTAG") && strcmp(signal, "GCAG") && strcmp(signal, "ATAC") && strcmp(signal, "CTAC") && strcmp(signal, "CTGC") && strcmp(signal, "GTAT"));
}

int main(int argc, char* argv[]){
	loadPac("test.fasta.pac");
	//printf("%s\n",Original_Text);
	char* R = "GCATCGATCAGCATGCATCGATGGCATGCGTGGGAGAATGTTTTTTTTTTTTTTTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
	Junc* junctions = extend(R, 20, 100, 560,680);
	
	delete [] junctions;
	return 0;
}

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

unsigned Get_File_Size(FILE* File)
{
	fseek (File , 0 , SEEK_END);
	unsigned Size = ftell (File);
	rewind (File);
	return Size;
}

FILE* File_Open(const char* File_Name,const char* Mode)
{
	FILE* Handle;
	Handle=fopen64(File_Name,Mode);
	if (Handle==NULL)
	{
		printf("File %s Cannot be opened ....\n",File_Name);
		exit(1);
	}
	else return Handle;
}

//Assume the coordinates are 0-based.
//All intervals used are closed interval.
Junc* extend(char* R, unsigned x, unsigned y, unsigned p, unsigned q){
	unsigned size = y-x-1;
	int misL[size+1], misR[size+1];
	char basesL[size+1], basesR[size+1];
	basesL[size] = 0;
	basesR[size] = 0;
	char str1[602];
	str1[600]=0;

	Get_Bases_ASCII(p+x+1, size, basesL);
	Get_Bases_ASCII(q-size, size, basesR);
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

	//Find the partitions that give the minimum mismatches.
	//Partitions store the size of the portion of the left extension.
	int min= misL[0] + misR[0];
	int parCount = 1;
	int partitions[size+1];
	for(int i=1; i<size+1; i++)
	{
		printf("%d\n",min);
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
	Junc *junctions = new Junc[parCount];
	for(int i = 0; i<parCount; i++)
	{
		junctions[i].p = p + x + partitions[i] + 1;
		junctions[i].q = q - size + partitions[i] - 1;
		junctions[i].r = x + partitions[i] + 1;
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

	return junctions;
}
