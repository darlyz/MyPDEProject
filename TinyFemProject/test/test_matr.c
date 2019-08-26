/*
 Copyright: Copyright (c) 2019
 Created: 2019-4-19
 Author: Zhang_Licheng
 All rights reserved
*/
#include "fem.h"

void matrcalc(Coor_Info, Node_Mesh, Dof_Tag, Equat_Set*, Materail, int);

int main (int argc, char* argv[])
{
	Coor_Info Coor;
	Node_Mesh Mesh;
	Dof_Tag   ID;
	Materail  Mate;
	Equat_Set Equa;
	int node_dof = 1;
	int *dof_tag;
	double *bnd_constr;
	
	Coor.dim = 2;
	Coor.total_nodes = 8;
	Coor.coordinate = (double*) malloc (Coor.dim*Coor.total_nodes*sizeof(double));
	Coor.coordinate[ 0] = 0.000000e+00;    Coor.coordinate[ 1] = 1.000000e+00; 
    Coor.coordinate[ 2] = 1.000000e+00;    Coor.coordinate[ 3] = 1.000000e+00; 
    Coor.coordinate[ 4] = 0.000000e+00;    Coor.coordinate[ 5] = 0.000000e+00; 
    Coor.coordinate[ 6] = 1.000000e+00;    Coor.coordinate[ 7] = 0.000000e+00; 
    Coor.coordinate[ 8] = 2.000000e+00;    Coor.coordinate[ 9] = 1.000000e+00; 
    Coor.coordinate[10] = 2.000000e+00;    Coor.coordinate[11] = 0.000000e+00; 
    Coor.coordinate[12] = 3.000000e+00;    Coor.coordinate[13] = 1.000000e+00; 
    Coor.coordinate[14] = 3.000000e+00;    Coor.coordinate[15] = 0.000000e+00;
	
	//dof_tag    = (int*)    malloc (Coor.total_nodes*node_dof*sizeof(int));
	//bnd_constr = (double*) calloc (Coor.total_nodes*node_dof,sizeof(double));
	//for (int i=0; i<Coor.total_nodes*node_dof; i++) dof_tag[i] = 1;
	//dof_tag[0] = -1;    bnd_constr[0] = 0.0;
    //dof_tag[2] = -1;    bnd_constr[2] = 0.0;
    //dof_tag[6] = -1;    bnd_constr[6] = 1;
    //dof_tag[7] = -1;    bnd_constr[7] = 1;
	
	Mesh.typeN = 1;

	Mesh.type = (Mesh_Type* )malloc(sizeof(Mesh_Type)*16);
    Mesh.elem_nodeN = (int* )malloc(sizeof(int )*16);
    Mesh.mesh_scale = (int* )malloc(sizeof(int )*16);
    Mesh.mesh_topo  = (int**)malloc(sizeof(int*)*16);

	Mesh.elem_nodeN[0] = 5;
	Mesh.mesh_scale[0] = 3;
	int total = 0;
	for (int i=0; i<Mesh.typeN; i++)
		total += Mesh.elem_nodeN[i] * Mesh.mesh_scale[i];
	Mesh.mesh_topo[0] = (int*) malloc (total*sizeof(int));
	Mesh.mesh_topo[0][ 0] = 4;    Mesh.mesh_topo[0][ 1] = 2;    Mesh.mesh_topo[0][ 2] = 1;    Mesh.mesh_topo[0][ 3] = 3;    Mesh.mesh_topo[0][ 4] = 1; 
    Mesh.mesh_topo[0][ 5] = 6;    Mesh.mesh_topo[0][ 6] = 5;    Mesh.mesh_topo[0][ 7] = 2;    Mesh.mesh_topo[0][ 8] = 4;    Mesh.mesh_topo[0][ 9] = 1; 
    Mesh.mesh_topo[0][10] = 8;    Mesh.mesh_topo[0][11] = 7;    Mesh.mesh_topo[0][12] = 5;    Mesh.mesh_topo[0][13] = 6;    Mesh.mesh_topo[0][14] = 1; 
	
	Mate.mate_varN = 2;
	Mate.mate_type = 1;
	Mate.mate = (double*) malloc (Mate.mate_varN*Mate.mate_type*sizeof(double));
	Mate.mate[0] = 8.8541878e-12;
	Mate.mate[1] = 0.0;
	
	Equa.total_equations = Coor.total_nodes * node_dof;
	Equa.matrix = (double**) malloc (Equa.total_equations * sizeof(double*));
	Equa.vector = (double* ) calloc (Equa.total_equations , sizeof(double));
	for (int i=0; i<Equa.total_equations; i++)
		Equa.matrix[i] = (double*) calloc (Equa.total_equations, sizeof(double));
	
	printf("begin!\n");

	matrcalc(Coor, Mesh, ID, &Equa, Mate, node_dof);
	
	//for (int i=0; i<Equa.total_equations; i++)
	//{
	//	for (int j=0; j<Equa.total_equations; j++)
	//		printf("%lf  ",Equa.matrix[i][j]);
	//	printf("%lf\n",Equa.vector[i]);
	//}
	
	
	printf("done!\n");
	
	return 1;
}