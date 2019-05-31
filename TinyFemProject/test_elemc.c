/*
 Copyright: Copyright (c) 2019
 Created: 2019-4-19
 Author: Zhang_Licheng
 All rights reserved
*/
#include "fem.h"

void set_gaus(gaus_calc*);
void set_elem(elem_info*);
void set_refr_shap(elem_info, gaus_calc,double**);
void elemcalc ();

int main (int argc, char* argv[])
{
	elem_info e_info;
	gaus_calc g_calc;
	elem_matr ematr;
	set_gaus(&g_calc);
	set_elem(&e_info);

	int node_dof = 1;
	int elem_dof = e_info.node_num*node_dof;
	ematr.matr_0 = (double*)calloc(elem_dof*elem_dof,sizeof(double));
	ematr.matr_1 = (double*)calloc(elem_dof*elem_dof,sizeof(double));
	ematr.matr_2 = (double*)calloc(elem_dof*elem_dof,sizeof(double));
	ematr.righ_vec = (double*)calloc(elem_dof,sizeof(double));
	
	double** refr_shap;
	refr_shap = (double**)malloc(g_calc.gaus_num*sizeof(double*));
	for (int gaus_i = 1; gaus_i <= g_calc.gaus_num; gaus_i ++)
		refr_shap[gaus_i-1] = (double*)malloc(e_info.node_num*(e_info.dim+1)*sizeof(double));
	set_refr_shap(e_info,g_calc,refr_shap);
	
	double node_coor[]={0.,0.};
	double coup_valu[]={0.,0.};
	double elem_mate[]={8.8541878e-12,0.};
	int elem_n = 2;

	elemcalc(node_coor,coup_valu,elem_mate,elem_n,refr_shap,g_calc,e_info,&ematr);

	for (int i=0; i<elem_dof*elem_dof; i++)
		printf("%e\n",ematr.matr_0[i]);

	for (int i=0; i<elem_dof; i++)
		printf("%e\n",ematr.righ_vec[i]);
	
	
	
	
	printf("done!\n");
	
	return 1;
}

