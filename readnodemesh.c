#include "PDE.h"

void jumprow(int n,FILE *fp);

void ReadNodeMesh(char* dat_file)
//void main(int argv, char* argc[])
{
	char   temp_char;
	int    temp_int;
	double temp_double;
	
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
	int NodeSum=0;
	for (int k=1; k<=TypesNum; k++)
	{
		jumprow(2,ReadData);
		fscanf(ReadData,"%d",&(NodMsh[k-1].Mesh_Scale));
		fscanf(ReadData,"%d",&(NodMsh[k-1].Elem_CompN));
		NodMsh[k-1].Mesh_Topo = (int*)calloc(sizeof(int),NodMsh[k-1].Mesh_Scale*NodMsh[k-1].Elem_CompN);
		for (int i=1; i<=NodMsh[k-1].Mesh_Scale; i++)
		{
			fscanf(ReadData,"%d",&temp_int);
			for (int j=1; j<=NodMsh[k-1].Elem_CompN; j++)
				fscanf(ReadData,"%d",&(NodMsh[k-1].Mesh_Topo[(i-1)*NodMsh[k-1].Elem_CompN + j-1]));
		}
		jumprow(1,ReadData);
	}
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