/*
 Copyright: Copyright (c) 2019
 Created: 2019-4-19
 Author: Zhang_Licheng
 All rights reserved
*/
#include "fem.h"

void matrcalc(Coor_Info, Node_Mesh, Dof_Tag, Equat_Set*, Materail, double*, int);

int main (int argc, char* argv[])
{
	Coor_Info Coor;
	Node_Mesh Mesh;
	Dof_Tag   ID;
	Materail  Mate;
	Equat_Set Equa;
	int node_dof = 1;
	
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

	ID.tag_nodn = 4;
	ID.dof_num  = node_dof;
	ID.tag_nods = (int*)    malloc (ID.tag_nodn*ID.dof_num*sizeof(int));
	ID.dof_tag  = (int*)    malloc (ID.tag_nodn*ID.dof_num*sizeof(int));
	ID.dof_val  = (double*) calloc (ID.tag_nodn*ID.dof_num,sizeof(double));
	ID.tag_nods[0] = 1 ;    ID.dof_tag[0] = -1;    ID.dof_val[0] = 0.0;
    ID.tag_nods[1] = 3 ;    ID.dof_tag[1] = -1;    ID.dof_val[1] = 0.0;
    ID.tag_nods[2] = 7 ;    ID.dof_tag[2] = -1;    ID.dof_val[2] = 1;
    ID.tag_nods[3] = 8 ;    ID.dof_tag[3] = -1;    ID.dof_val[3] = 1;

	double *result = (double*)calloc(Coor.total_nodes*node_dof,sizeof(double));
	
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

	Equa.node_equa_index = (int**)malloc(Coor.total_nodes * sizeof(int*));
	for (int i=0; i<Coor.total_nodes; i++)
        Equa.node_equa_index[i] = (int*)malloc(node_dof*sizeof(int));
	
	Equa.node_equa_index[0][0] = -1;
	Equa.node_equa_index[1][0] =  1;
	Equa.node_equa_index[2][0] = -1;
	Equa.node_equa_index[3][0] =  2;
	Equa.node_equa_index[4][0] =  3;
	Equa.node_equa_index[5][0] =  4;
	Equa.node_equa_index[6][0] = -1;
	Equa.node_equa_index[7][0] = -1;

	Equa.total_equations = 4;
	Equa.matrix = (double**) malloc (Equa.total_equations * sizeof(double*));
	Equa.vector = (double* ) calloc (Equa.total_equations , sizeof(double));
	Equa.row_nontriaval = (int* )malloc(Equa.total_equations * sizeof(int));
    Equa.column_index   = (int**)malloc(Equa.total_equations * sizeof(int*));

	Equa.row_nontriaval[0] = 4;  Equa.column_index[0] = (int*)malloc(4 * sizeof(int));  Equa.matrix[0] = (double*)calloc(4 , sizeof(double));
	Equa.row_nontriaval[1] = 4;  Equa.column_index[1] = (int*)malloc(4 * sizeof(int));  Equa.matrix[1] = (double*)calloc(4 , sizeof(double));
	Equa.row_nontriaval[2] = 4;  Equa.column_index[2] = (int*)malloc(4 * sizeof(int));  Equa.matrix[2] = (double*)calloc(4 , sizeof(double));
	Equa.row_nontriaval[3] = 4;  Equa.column_index[3] = (int*)malloc(4 * sizeof(int));  Equa.matrix[3] = (double*)calloc(4 , sizeof(double));
	
	Equa.column_index[0][0] = 1;	Equa.column_index[0][1] = 2; Equa.column_index[0][2] = 3;	Equa.column_index[0][3] = 4;
	Equa.column_index[1][0] = 1;	Equa.column_index[1][1] = 2; Equa.column_index[1][2] = 3;	Equa.column_index[1][3] = 4;
	Equa.column_index[2][0] = 1;	Equa.column_index[2][1] = 2; Equa.column_index[2][2] = 3;	Equa.column_index[2][3] = 4;
	Equa.column_index[3][0] = 1;	Equa.column_index[3][1] = 2; Equa.column_index[3][2] = 3;	Equa.column_index[3][3] = 4;
	
	printf("begin!\n");

	matrcalc(Coor, Mesh, ID, &Equa, Mate, result, node_dof);
	
	//for (int i=0; i<Equa.total_equations; i++)
	//{
	//	for (int j=0; j<Equa.total_equations; j++)
	//		printf("%lf  ",Equa.matrix[i][j]);
	//	printf("%lf\n",Equa.vector[i]);
	//}
	
	
	printf("done!\n");
	
	return 1;
}

int Binary_Search(int* dest, int dest_len, int key) {
    int idx_start = 0;
    int idx_end   = dest_len-1;
    int idx_mid   = -1;
    int insert    = 0;

    if (dest_len == 0)
    	return 0;

    while(1)
    {
        idx_mid = (idx_start + idx_end)/2;

        if (key < dest[idx_mid])
            idx_end   = idx_mid - 1;

        else if (key > dest[idx_mid])
            idx_start = idx_mid + 1;

        else
        	return idx_mid;
    }
}