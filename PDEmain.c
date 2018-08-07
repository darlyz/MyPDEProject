#include "PDE.h"
void ReadJson(char*);
void ReadNodeMesh(char*);
void WriteMesh   (char*);
//int WriteResult (Node_Unknow*);
//int InitialPDE  (Coor, Node_Mesh*, Node_Unknow*, Materail, Time_List, Equation_Set);
//int AssembleMX  (Coor, Node_Mesh*, Node_Unknow*, Materail, Time_List, Equation_Set);
//int SolveUnknow (Equation_Set);
//int PackResault (Coor);

int main (int argc, char* argv[])
{
	ReadJson("testpro.json");
	FieldNum=2;
	ReadNodeMesh("test.dat");
	//InitialPDE  (Coor0, NodMsh, NodF, Mate0, TList0, EqSet0);
	//AssembleMX  (Coor0, NodMsh, NodF, Mate0, TList0, EqSet0);
	//SolveUnknow (EqSet0);
	//PackResault (NodF);
	WriteMesh ("test.post.msh");
	//WriteResult (Coor0, NodMsh);
	
	return 1;
}