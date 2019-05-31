/*
 Copyright: Copyright (c) 2019
 Created: 2019-4-19
 Author: Zhang_Licheng
 All rights reserved
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

struct elem_matr
{
    double *matr_0;
    double *matr_1;
    double *matr_2;
    double *matr_3;
    double *left_matr;
    double *righ_vec;
};
typedef struct elem_matr elem_matr;

struct gaus_info
{
    int gaus_num;
    int refc_num;
    double *gaus_coor;
    double *gaus_weig;
};
typedef struct gaus_info gaus_info;

struct elem_info
{
    int node_num;
    int dim;
	double *node_coor;
	double *refc_coor;
};
typedef struct elem_info elem_info;

struct coor_info
{
    int dim;
    int node_num;
    double *coordinate;
};
typedef struct coor_info coor_info;

struct Node_Mesh
{
	int Mesh_TypeN;
	int *Elem_NodeN;
    int *Mesh_Scale;
    int *Mesh_Topo;
};
typedef struct Node_Mesh Node_Mesh;

typedef enum matr_type {lump=1, dist} matr_type;