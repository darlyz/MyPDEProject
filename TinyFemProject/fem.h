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
#include <ctype.h>
#include <stdbool.h>

typedef enum Mesh_Type {P1=1, L2, L3, T3, T6, Q4, Q8, Q9, W4, W10, C8, C20, C27, H6, H15, H18, P5, P13, P14} Mesh_Type;
typedef enum Matr_Type {lump=1, dist} Matr_Type;

typedef struct Elem_Matr
{
    double *matr_0, *matr_1, *matr_2, *matr_3;   // stif, damp, mass, high order mass
    double *left_matr, *righ_vect;
}Elem_Matr;

typedef struct Gaus_Info
{
    int gaus_num, refc_num;         // gaussian node number, reference dim
    double *gaus_coor, *gaus_weig;  // gaussian node coordinates, weight
}Gaus_Info;

typedef struct Elem_Info
{
    int dim, node_cont, *elem_node; // element dim, node count, nodes S.N. and materail S.N.
	double *node_coor, *coup_valu;  // element nodes coordinates, (XXcoupled value)
    double *elem_mate, **refr_shap; // element materail parameter values, reference shap values
}Elem_Info;

typedef struct Coor_Info
{
    int dim, total_nodes;           // coordinates dim, total nodes count
    double *coordinate;             // all nodes coordinates
}Coor_Info;

typedef struct Node_Mesh
{
	Mesh_Type *type;                // element types
	int typeN, *elem_nodeN;         // type count, element nodes count per type
    int *mesh_scale, **mesh_topo;   // element count per type, mesh topology of all types
}Node_Mesh;

typedef struct Materail
{
	int mate_type, mate_varN;   // materail count, parameter count
	double* mate;               // materail parameter values
}Materail;

typedef struct Equation_Set
{
	int      total_equations;   // total number of equations
    int      total_nontriaval;  // total number of nontrivials
    int     *row_nontriaval;    // number of nontrivials per row        [eq   S.N.]           = <count>
    int    **column_index;      // nontrivials column index per row     [eq   S.N.][index]    = <clm S.N.>
    int    **node_equa_index;   // equations S.N. of dofs on every node [node S.N.][dof S.N.] = <eq  S.N.>
	double **matrix,*vector;
}Equat_Set;

typedef struct Dof_Tag
{
    int     tag_nodn;  // total number of nodes has constraint dof
    int    *tag_nods;  // S.N. of nodes has constraint dof          [index]    = <node S.N.>
    int     dof_num;   // dof number a node has
    int    *dof_tag;   // ID of dof                                 [i*dofn+j] = <ID S.N.  >
    double *dof_val;   // value of dof                              [i*dofn+j] = <ID value >
}Dof_Tag;

typedef struct init_data
{
    double *init_0, *init_1, *init_2;
    int dof_num;
    int total_nodes;
    int init_oder;
}Init_Data;