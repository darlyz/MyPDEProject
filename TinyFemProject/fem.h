/*
 Copyright: Copyright (c) 2019
 Created: 2019-4-19
 Author: Zhang_Licheng
 All rights reserved
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

typedef enum Mesh_Type {P1=1, L2, L3, T3, T6, Q4, Q8, Q9, W4, W10, C8, C20, C27, H6, H15, H18} Mesh_Type;
typedef enum Matr_Type {lump=1, dist} Matr_Type;

typedef struct Elem_Matr
{
    double *matr_0, *matr_1, *matr_2, *matr_3;
    double *left_matr, *righ_vect;
}Elem_Matr;

typedef struct Gaus_Info
{
    int gaus_num, refc_num;
    double *gaus_coor, *gaus_weig;
}Gaus_Info;

typedef struct Elem_Info
{
    int dim, node_cont, *elem_node;
	double *node_coor, *coup_valu, *elem_mate, **refr_shap;
}Elem_Info;

typedef struct Coor_Info
{
    int dim, total_nodes;
    double *coordinate;
}Coor_Info;

typedef struct Node_Mesh
{
	Mesh_Type *type;
	int typeN, *elem_nodeN, *mesh_scale, *mesh_topo;
}Node_Mesh;

typedef struct Materail
{
	int mate_type, mate_varN;
	double* mate;
}Materail;

typedef struct Equation_Set
{
	int total_equations, total_nontriaval, *triaval_per_row, **column_index;
	double **matrix,*vector;
}Equat_Set;

typedef struct Dirichlet
{
    int  dirichlet_nodn;
    int *dirichlet_nods;
    int    *dof_tag;
    double *dof_val;
}Dirichlet;