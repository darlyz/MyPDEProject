/*
 Copyright: Copyright (c) 2019
 Created: 2019-4-19
 Author: Zhang_Licheng
 All rights reserved
*/
#include "fem.h"

void matrcalc(Coor_Info, Node_Mesh, Materail, Equat_Set, int, int*, double*);

int main (int argc, char* argv[])
{
	Coor_Info Coor;
	Node_Mesh Mesh;
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
	
	dof_tag    = (int*)    malloc (Coor.total_nodes*node_dof*sizeof(int));
	bnd_constr = (double*) calloc (Coor.total_nodes*node_dof,sizeof(double));
	for (int i=0; i<Coor.total_nodes*node_dof; i++) dof_tag[i] = 1;
	dof_tag[0] = -1;    bnd_constr[0] = 0.0;
    dof_tag[2] = -1;    bnd_constr[2] = 0.0;
    dof_tag[6] = -1;    bnd_constr[6] = 1;
    dof_tag[7] = -1;    bnd_constr[7] = 1;
	
	Mesh.typeN = 1;
	Mesh.elem_nodeN = (int*) malloc (Mesh.typeN*sizeof(int));
	Mesh.mesh_scale = (int*) malloc (Mesh.typeN*sizeof(int));
	Mesh.elem_nodeN[0] = 5;
	Mesh.mesh_scale[0] = 3;
	int total = 0;
	for (int i=0; i<Mesh.typeN; i++)
		total += Mesh.elem_nodeN[i] * Mesh.mesh_scale[i];
	Mesh.mesh_topo = (int*) malloc (total*sizeof(int));
	Mesh.mesh_topo[ 0] = 4;    Mesh.mesh_topo[ 1] = 2;    Mesh.mesh_topo[ 2] = 1;    Mesh.mesh_topo[ 3] = 3;    Mesh.mesh_topo[ 4] = 1; 
    Mesh.mesh_topo[ 5] = 6;    Mesh.mesh_topo[ 6] = 5;    Mesh.mesh_topo[ 7] = 2;    Mesh.mesh_topo[ 8] = 4;    Mesh.mesh_topo[ 9] = 1; 
    Mesh.mesh_topo[10] = 8;    Mesh.mesh_topo[11] = 7;    Mesh.mesh_topo[12] = 5;    Mesh.mesh_topo[13] = 6;    Mesh.mesh_topo[14] = 1; 
	
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

	matrcalc(Coor, Mesh, Mate, Equa, node_dof, dof_tag, bnd_constr);
	
	//for (int i=0; i<Equa.total_equations; i++)
	//{
	//	for (int j=0; j<Equa.total_equations; j++)
	//		printf("%lf  ",Equa.matrix[i][j]);
	//	printf("%lf\n",Equa.vector[i]);
	//}
	
	
	printf("done!\n");
	
	return 1;
}