/*
 Copyright: Copyright (c) 2018
 Created: 2018-4-19
 Author: Zhang_Licheng
 Title: complex tensor expression in dummy index summation
 All rights reserved
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include"cJSON.h"

typedef enum Mesh_Type {P1=1, L2, L3, T3, T6, Q4, Q8, Q9, W4, W10, C8, C20, C27, H6, H15, H18} Mesh_Type;
typedef enum Matr_Type {lump=1, dist} Matr_Type;

typedef struct Coor
{
	int dim, total_nodes;
	double *coor;                                      // *coor : [x1,x2,x3,...,y1,y2,y3,...,z1,z2,z3,...]
}Coor;

typedef struct Gaus_Info
{
    int gaus_num, refc_num;
    double *gaus_coor, *gaus_weig;
}Gaus_Info;

typedef struct Elem_Info
{
    int dim, node_cont, *elem_node;
	double *node_coor, *coup_valu, *elem_mate **refr_shap;
}Elem_Info;

typedef struct Node_Unknow
{
	int dof, *unknowns_status;
	double *status_value;
	char **unknow_name;
}Node_Unknow;

typedef struct Edge_Unknow
{
	int dof, total_edge, *unknowns_status, *edge_dirt, *edge_list;
	char **unknow_name;
	double *status_value;
}Edge_Unknow;

typedef struct Face_Unknow
{
	int dof, total_face, *unknowns_status, *face_norm, *face_list;
	char **unknow_name;
	double *status_value;
}Face_Unknow;

typedef struct Volm_Unknow
{
	int dof, total_volm, *unknown_status, *volm_list;
	char **unknow_name;
	double *status_value;
}Volm_Unknow;

typedef struct Node_Mesh
{
	Mesh_Type *type;
	int typeN, *elem_nodeN, *mesh_scale, *mesh_topo;
}Node_Mesh;

typedef struct Edge_Mesh
{
	Mesh_Type *type;
	int typeN, *elem_edgeN, *mesh_scale, *mesh_topo;
}Edge_Mesh;

typedef struct Face_Mesh
{
	Mesh_Type *type;
	int typeN, *elem_faceN, *mesh_scale, *mesh_topo;
}Face_Mesh;

typedef struct Volm_Mesh
{
	Mesh_Type *type;
	int typeN, *mesh_scale, *mesh_topo;
}Volm_Mesh;

typedef struct Materail
{
	int mate_type, mate_varN;
	double* mate;
}Materail;

typedef struct Time_List
{
	int time_intervals;
	double inti_time, end_time, current_time, *time_series;
}Time_List;

typedef struct Equation_Set
{
	int total_equations, total_nontriaval, *triaval_per_row, **column_index;
	double **matrix;
}Equation_Set;

typedef struct Elem_Matrix
{
    double *matr_0, *matr_1, *matr_2, *matr_3;
    double *left_matr, *righ_vec;
}Elem_Matr;

coor coor0;
node_unknow *nodf;
edge_unknow *edgf;
face_unknow *facf;
volm_unknow *volf;
node_mesh *nodmsh;
edge_mesh *edgmsh;
face_mesh *facmsh;
volm_mesh *volmsh;
materail mate0;
time_list tlist0;
equation_set eqset0;

int fieldnum;
int typesnum;
int elemtype[10];
int *dof;         // dof of fields
int *idnodenum;   // the total number of nodes which assigned id
int **idnodesn;   // the serial number of nodes which assigned id
int **ids;        // the id of dofs
double **ubf;     // the value of ids	
