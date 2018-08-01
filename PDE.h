#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

typedef struct Coor
{
	int Dim, Total_Nodes;
	double *Coor;                                      // *Coor : [x1,x2,x3,...,y1,y2,y3,...,z1,z2,z3,...]
}Coor;
Coor Coor0;

typedef struct Node_Unknow
{
	int DoF;
	int *Unknowns_Status;
	double *Status_Value;
	char **Unknow_Name;
}Node_Unknow;
Node_Unknow *NodF;

typedef struct Edge_Unknow
{
	int DoF, Total_Edge;
	int *Unknowns_Status, *Edge_Dirt, *Edge_List;
	double *Status_Value;
	char **Unknow_Name;
}Edge_Unknow;
Edge_Unknow *EdgF;

typedef struct Face_Unknow
{
	int DoF, Total_Face;
	int *Unknowns_Status, *Face_Norm, *Face_List;
	double *Status_Value;
	char **Unknow_Name;
}Face_Unknow;
Face_Unknow *FacF;

typedef struct Volm_Unknow
{
	int DoF, Total_Volm;
	int *Unknown_Status, *Volm_List;
	double *Status_Value;
	char **Unknow_Name;
}Volm_Unknow;
Volm_Unknow *VolF;

typedef struct Node_Mesh
{
	int Mesh_Type;
	int Node_Count,Mesh_Scale,*Mesh_Topo;
}Node_Mesh;
Node_Mesh *NodMsh;

typedef struct Edge_Mesh
{
	int Mesh_Type;
	int Edge_Count,Mesh_Scale,*Mesh_Topo;
}Edge_Mesh;
Edge_Mesh *EdgMsh;

typedef struct Face_Mesh
{
	int Mesh_Type;
	int Face_Count,Mesh_Scale,*Mesh_Topo;
}Face_Mesh;
Face_Mesh *FacMsh;

typedef struct Volm_Mesh
{
	int Mesh_Type;
	int Mesh_Scale,*Mesh_Topo;
}Volm_Mesh;
Volm_Mesh VolMsh;

typedef struct Materail
{
	int Mate_Type;
	double* Mate;
}Materail;
Materail Mate0;

typedef struct Time_List
{
	int Time_Intervals;
	double Inti_Time,End_Time,Current_Time;
	double* Time_Series;
}Time_List;
Time_List TList0;

typedef struct Equation_Set
{
	int Total_Equations,Total_Nontriaval;
	int *Triaval_Per_Row;
	int **Column_Index;
	double **Matrix;
}Equation_Set;
Equation_Set EqSet0;

int FieldNum;
int ElemType[10];
int *dof;         // DOF of Fields
int *IDNodeNum;   // the total number of Nodes which assigned ID
int **IDNodeSN;   // the Serial Number of Nodes which assigned ID
int **IDs;        // the ID of DOFs
double **Ubf;     // the value of IDs	
enum Mesh_Type{P1=1, L2, L3, T3, T6, Q4, Q8, Q9, W4, W10, C8, C20, C27, H6, H15, H18};