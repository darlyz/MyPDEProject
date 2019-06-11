/*
 Copyright: Copyright (c) 2018
 Created: 2018-4-19
 Author: Zhang_Licheng
 Title: complex tensor expression in dummy index summation
 All rights reserved
*/
#include "PDE.h"
//#include <omp.h>

void set_gaus(Gaus_Info*);
void set_refr_shap(Elem_Info, Gaus_Info);
void SetEmatr(Elem_Matr*, int, Matr_Type*);
void SetEinfo(Elem_Info*, int, int, int, int);

void ComposeEquation(
	Coor         Coor0,
	Node_Mesh    Nmesh,
	Node_Unknow  Nuknw,
	Materail     Mate0,
	Time_List    Tlist,
	Equation_Set Eqset
){
	
	
	int pre_node = 0
	for (int type_i=1; type_i<=Nmesh.typeN; type_i++)
	{
		int elem_nodeN = Nmesh.elem_nodeN[type_i-1];
        int mesh_scale = Nmesh.mesh_scale[type_i-1];
        int ematr_size = (elem_nodeN-1) * Nuknw.dof;
		
		Matr_Type M_type[4]={dist,lump,lump,lump};
		
		Gaus_Info Gauss;
		Elem_Info Einfo;
		Elem_Matr Ematr;
		
		set_gaus(&Gauss); // !!!!
		SetEmatr(&Ematr, ematr_size, M_type);
		SetEinfo(&Einfo, Coor0.dim, elem_nodeN, Mate0.mate_varN, Gauss.gaus_num);
		
		for (int elem_i; elem_i<=mesh_scale; elem_i++)
		{
			memcpy(Einfo.elem_node,
				   &Nmesh.mesh_topo[pre_node+(elem_i-1)*elem_nodeN],
				   elem_nodeN*sizeof(int));
			memcpy(Einfo.elem_mate,
				   &Mate0.mate[Einfo.elem_node[Einfo.node_cont]-1],
				   Mate0.mate_varN*sizeof(double))
			for (int node_i=1; node_i<elem_nodeN; node_i++)
				memcpy(Einfo.node_coor,
					   &Coor0.coor[(Einfo.elem_node[node_i-1]-1)*Einfo.dim],
					   Einfo.dim*sizeof(double));
			set_refr_shap(Einfo,Gauss);
			
			ElemCalculate(Gauss, Einfo, &Ematr, elem_i);
			
			memset(Ematr.left_matr,0.0,ematr_size*ematr_size*sizeof(double));
			
			Compose
			for (int i=1; i<=ematr_size; ++i)
            {
                for (int j=1; j<=ematr_size; ++j)
                {
                    Ematr.left_matr[(i-1)*ematr_size+j-1] += 
                                          +es[(i-1)*ematr_size+j-1];
                }
                Ematr.left_matr[(i-1)*ematr_size+i-1] += 0
                                      ;
            }
			
			// ..........
		
		}
		
		ClearEmatr(&Ematr);
		ClearEinfo(&Einfo, Gauss.gaus_num)
		pre_node += elem_nodeN * mesh_scale;
	}
}

void SetEmatr(Elem_Matr &Ematr, int elem_dof, Matr_Type *M_type)
{
    Ematr->matr_0    = (double*)malloc(pow(elem_dof,M_type[0])*sizeof(double));
    Ematr->matr_1    = (double*)malloc(pow(elem_dof,M_type[1])*sizeof(double));
    Ematr->matr_2    = (double*)malloc(pow(elem_dof,M_type[2])*sizeof(double));
    Ematr->matr_3    = (double*)malloc(pow(elem_dof,M_type[3])*sizeof(double));
    Ematr->left_matr = (double*)malloc(pow(elem_dof,M_type[0])*sizeof(double));
    Ematr->righ_vec  = (double*)malloc(elem_dof*sizeof(double));
}
    
void ClearEmatr(Elem_Matr *Ematr)
{
    free(Ematr->matr_0);
    free(Ematr->matr_1);
    free(Ematr->matr_2);
    free(Ematr->matr_3);
	free(Ematr->left_matr);
	free(Ematr->righ_vec );
}

void SetEinfo(Elem_Info *Einfo, int dim, int elem_nodeN, int mate_varN, int gaus_num)
{
	Einfo->dim = dim;
	Einfo->node_cont  = elem_nodeN - 1;
	
	Einfo->elem_node = (int*    )malloc(elem_nodeN * sizeof(int))
	Einfo->node_coor = (double* )malloc(Einfo->node_cont * Einfo->dim * sizeof(double));
	//Einfo->coup_valu = (double*)malloc(Einfo->node_cont * coupN * sizeof(double));
	Einfo->elem_mate = (double* )malloc(mate_varN * sizeof(double));
	Einfo->refr_shap = (double**)malloc(gaus_num * sizeof(double*));
	for (int gaus_i=1; gaus_i<=gaus_num; gaus_i++)
		Einfo->refr_shap[gaus_i-1] = (double*)malloc(Einfo->node_cont * (Einfo->dim+1) * sizeof(double));
}

void ClearEinfo(Elem_Info *Einfo, int gaus_num)
{
	free(Einfo->elem_node);
    free(Einfo->node_coor);
    //free(Einfo->coup_valu);
    free(Einfo->elem_mate);
	for (int gaus_i=1; gaus_i<=gaus_num; gaus_i++)
		free(Einfo->refr_shap[gaus_i-1]);
	free(Einfo->refr_shap);
}