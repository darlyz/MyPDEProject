#include "PDE.h"

void jumprow(int n,FILE *fp);

void ReadNodeMesh(char* pre_file, char* dat_file)
//void main(int argv, char* argc[])
{
	char   temp_char;
	int    temp_int;
	double temp_double;
	
	//---------- need to be determined by pre file
	Nodmsh0.Mesh_Type = 1;
	FieldNum    = 1;
	ElemType[1] = Q4;
	//ElemType[2] = W4;
	//ElemType[3] = Q4;
	//ElemType[4] = T3;
	//---------

	FILE *ReadData;
	if((ReadData=fopen(dat_file,"r"))==NULL)
	//if((ReadData=fopen(argc[2],"r"))==NULL)
		return;
	
	jumprow(2,ReadData);

	// read coor
	int ReadCoor=1;
	if(ReadCoor){
	fscanf(ReadData,"%d",&Coor0.Total_Nodes);
	fscanf(ReadData,"%d",&Coor0.Dim);
	Coor0.Coor =(double*)malloc(sizeof(double)*Coor0.Total_Nodes*Coor0.Dim);
	memset      (Coor0.Coor, 0, sizeof(double)*Coor0.Total_Nodes*Coor0.Dim);
	for (int i=0; i<Coor0.Total_Nodes; i++)
	{
		fscanf(ReadData,"%d",&temp_int);
		for (int j=0; j<Coor0.Dim; j++)
			fscanf(ReadData,"%lf",&(Coor0.Coor[i + j*Coor0.Total_Nodes]));
	}
	}

	// read ubf, usually we only care for Field A, but we need to jump the data
	int ReadBoundaryForce=1;
	if(ReadBoundaryForce){
	dof = (int*)malloc(sizeof(int)*FieldNum);          // DOF of Fields
	memset    (dof, 0, sizeof(int)*FieldNum);

	IDNodeNum = (int*)malloc(sizeof(int)*FieldNum);    // the total number of Nodes which assigned ID
	memset    (IDNodeNum, 0, sizeof(int)*FieldNum);

	IDNodeSN = (int**)malloc(sizeof(int*)*FieldNum);   // the Serial Number of Nodes which assigned ID
	memset     (IDNodeSN, 0, sizeof(int*)*FieldNum);

	IDs = (int**)malloc(sizeof(int*)*FieldNum);        // the ID of DOFs
	memset     (IDs, 0, sizeof(int*)*FieldNum);

	Ubf = (double**)malloc(sizeof(double*)*FieldNum);        // the value of IDs
	memset        (Ubf, 0, sizeof(double*)*FieldNum);

	for (int k=0; k<FieldNum; k++)
	{
		jumprow(2,ReadData);
		fscanf(ReadData,"%d",&(IDNodeNum[k]));
		fscanf(ReadData,"%d",&(dof[k]));

		if(IDNodeNum[k]!=0)
		{
			IDNodeSN[k] = (int*)malloc(sizeof(int)*IDNodeNum[k]);
			memset    (IDNodeSN[k], 0, sizeof(int)*IDNodeNum[k]);

			IDs[k] = (int*)malloc(sizeof(int)*IDNodeNum[k]*dof[k]);
			memset   (IDs[k], 0,  sizeof(int)*IDNodeNum[k]*dof[k]);

			Ubf[k] = (double*)malloc(sizeof(double)*IDNodeNum[k]*dof[k]);
			memset   (Ubf[k], 0,     sizeof(double)*IDNodeNum[k]*dof[k]);
		

			for (int i=0; i<IDNodeNum[k]; i++)
			{
				fscanf(ReadData,"%d",&(IDNodeSN[k][i]));
				for (int j=0; j<dof[k]; j++)
					fscanf(ReadData,"%d",&(IDs[k][i*dof[k] + j]));
			}
		
			jumprow(3,ReadData);
			
			for (int i=0; i<IDNodeNum[k]; i++)
			{
				fscanf(ReadData,"%d",&temp_int);
				for (int j=0; j<dof[k]; j++)
					fscanf(ReadData,"%lf",&(Ubf[k][i*dof[k] + j]));
			}
		}
		else
			jumprow(3,ReadData);
		
	}
	}

	// read elem, usually we only care for Field A, so we ignore others
	int ReadMesh=1;
	// method 1
	if(ReadMesh){
	long long position;
	position=ftell(ReadData);
	int NodeSum=0;
	for (int k=1; k<=Nodmsh0.Mesh_Type; k++)
	{
		jumprow(2,ReadData);
		fscanf(ReadData,"%d",&(Nodmsh0.Mesh_Scale[k]));
		fscanf(ReadData,"%d",&(Nodmsh0.Node_Count[k]));
		jumprow(Nodmsh0.Mesh_Scale[k]+1,ReadData);
		NodeSum += Nodmsh0.Mesh_Scale[k]*Nodmsh0.Node_Count[k];
	}

	int aa = fseek(ReadData, position, SEEK_SET);

	Nodmsh0.Mesh_Topo = (int*)malloc(sizeof(int)*(NodeSum+1));
	memset    (Nodmsh0.Mesh_Topo, 0, sizeof(int)*(NodeSum+1));
	
	NodeSum=0;
	for (int k=1; k<=Nodmsh0.Mesh_Type; k++)
	{
		jumprow(3,ReadData);
		for (int i=1; i<=Nodmsh0.Mesh_Scale[k]; i++)
		{
			fscanf(ReadData,"%d",&temp_int);
			for (int j=1; j<=Nodmsh0.Node_Count[k]; j++)
				fscanf(ReadData,"%d",&(Nodmsh0.Mesh_Topo[NodeSum + (i-1)*Nodmsh0.Node_Count[k] + j]));
		}
		jumprow(1,ReadData);
		NodeSum += Nodmsh0.Mesh_Scale[k]*Nodmsh0.Node_Count[k];
	}
	
	// method 2
	/*
	int NodeSum=0;
	int **Mesh_Topo;
	Mesh_Topo = (int**)malloc(sizeof(int*)*Nodmsh0.Mesh_Type);
	memset     (Mesh_Topo, 0, sizeof(int*)*Nodmsh0.Mesh_Type);
	for (int k=1; k<=Nodmsh0.Mesh_Type; k++)
	{
		jumprow(2,ReadData);
		fscanf(ReadData,"%d",&(Nodmsh0.Mesh_Scale[k]));
		fscanf(ReadData,"%d",&(Nodmsh0.Node_Count[k]));
		Mesh_Topo[k-1] = (int*)malloc(sizeof(int)*Nodmsh0.Mesh_Scale[k]*Nodmsh0.Node_Count[k]);
		memset    (Mesh_Topo[k-1], 0, sizeof(int)*Nodmsh0.Mesh_Scale[k]*Nodmsh0.Node_Count[k]);
		for (int i=1; i<=Nodmsh0.Mesh_Scale[k]; i++)
		{
			fscanf(ReadData,"%d",&temp_int);
			for (int j=1; j<=Nodmsh0.Node_Count[k]; j++)
				fscanf(ReadData,"%d",&(Mesh_Topo[k-1][(i-1)*Nodmsh0.Node_Count[k] + j-1]));
		}
		jumprow(1,ReadData);
		NodeSum+=Nodmsh0.Mesh_Scale[k]*Nodmsh0.Node_Count[k];
	}
	
	Nodmsh0.Mesh_Topo = (int*)malloc(sizeof(int)*(NodeSum+1));
	
	NodeSum=0;
	for (int k=1; k<=Nodmsh0.Mesh_Type; k++)
	{
		memcpy(&(Nodmsh0.Mesh_Topo[NodeSum+1]),Mesh_Topo[k-1],sizeof(int)*Nodmsh0.Mesh_Scale[k]*Nodmsh0.Node_Count[k]);
		NodeSum+=Nodmsh0.Mesh_Scale[k]*Nodmsh0.Node_Count[k];
	}

	for (int k=1; k<=Nodmsh0.Mesh_Type; k++)
	{
		free(Mesh_Topo[k-1]);
		Mesh_Topo[k-1]=NULL;
	}
	free(Mesh_Topo);
	Mesh_Topo=NULL;
	*/
	}
	
	return;
}


void jumprow(int n,FILE *fp) // jump n rows
{
	char voidchar[255];
	int i;
	for (i=1;i<=n;i++)
		fgets(voidchar,255,fp);
	return;
}