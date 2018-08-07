#include "PDE.h"

void ReadJson(char*);
void ReadNodeMesh(char*);

int main()
{
	ReadJson("testpro.json");
	
	if(NodMsh==NULL) printf("NodMsh is not creat\n");
	if(EdgMsh==NULL) printf("EdgMsh is not creat\n");
	if(FacMsh==NULL) printf("FacMsh is not creat\n");
	if(VolMsh==NULL) printf("VolMsh is not creat\n");
	
	FieldNum=2;
	ReadNodeMesh("test.dat");
	printf("%d\n",NodMsh[3].Mesh_Topo[(NodMsh[3].Mesh_Scale*NodMsh[3].Node_Count)-2]);
	
	//if(NodMsh!=NULL) free(NodMsh);
	//if(EdgMsh!=NULL) free(EdgMsh);
	//if(FacMsh!=NULL) free(FacMsh);
	//if(VolMsh!=NULL) free(VolMsh);
	return 1;
}
