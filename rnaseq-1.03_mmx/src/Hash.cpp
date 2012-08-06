#include "Hash.h"
#include "assert.h"

map_it Hash::Begin ()
{
	return Junctions.begin();
}

void Hash::Insert (OP Location,int Paring)
{
	JStat T;
	assert(Paring >=255 && Paring <=262); assert(Location.x <= Location.y);
	Junc_I=Junctions.find(Location);
	if (Junc_I == Junctions.end()) //New entry..
	{
		T.Count=1;T.Junc_Type=Paring;
		Junctions[Location]=T;
	}
	else 
	{
		assert((Junc_I->second).Count>0 && (Junc_I->second).Junc_Type==Paring);
		(Junc_I->second).Count++;
	}
}


void Hash::Delete (OP Location)
{
	Junc_I=Junctions.find(Location);
	if (Junc_I != Junctions.end()) Junctions.erase(Junc_I);
}

bool Hash::Init_Iterate(OP & Location,JStat & Data)
{
	Junc_I=Junctions.begin();
	if (Junc_I == Junctions.end()) return false; 
	else 
	{
		Last=Junc_I;//save last pos...
		Location = Junc_I->first;
		Data = Junc_I->second;
		Junc_I++;
		return true;
	}
}

bool Hash::Iterate(OP & Location,JStat & Data)
{
	if(Junc_I == Junctions.end()) return false;
	Last=Junc_I;//save last pos...
	Location = Junc_I->first;
	Data = Junc_I->second;
	Junc_I++;
	return true;
}

//}-----------------------------  Classes  -------------------------------------------------/
