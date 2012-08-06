#ifndef __HASH_GUARD__
#define __HASH_GUARD__
#include <map>
using namespace std;

struct OP//Ordered pairs..
{
	unsigned x;
	unsigned y;
	//int Motif;
};

struct OPX//Ordered pairs..
{
	unsigned x;
	unsigned y;
	int Motif;
};

struct JStat//Junction stats..
{
	int Count;
	int Junc_Type;
};

struct OP_Cmp
{
	bool operator()( OP OP1, OP OP2)
	{
		if (OP1.x == OP2.x)
		{
			return OP1.y< OP2.y;
		}
		else return OP1.x < OP2.x;
	} 
};

typedef map <OP,JStat>::iterator map_it; 
class Hash
{

	public:
	map <OP,JStat,OP_Cmp> Junctions;
	map_it Junc_I,JJ;
	map_it Last;

	void Insert (OP Location,int Paring);
	void Delete (OP Location);
	bool Init_Iterate(OP & Location,JStat & Data);
	bool Iterate(OP & Location,JStat & Data);
	map_it Begin();


};
#endif
