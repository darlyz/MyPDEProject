/*
 Copyright: Copyright (c) 2019
 Created: 2019-4-19
 Author: Zhang_Licheng
 All rights reserved
*/
#include "fem.h"

void set_gaus(Gaus_Info*);
void set_refr_shap(double**, double*, int, int, int);
void show_elem(Elem_Info, int, int, int, int);
void elemcalc ();

int main (int argc, char* argv[])
{
	Elem_Info E_info;
	Gaus_Info G_info;
	Elem_Matr E_matr;
	set_gaus(&G_info);
	
	E_info.node_cont = 4;
	E_info.dim = 2;
	E_info.node_coor = (double*) malloc (E_info.node_cont * E_info.dim * sizeof(double));
	
	E_info.node_coor[(1-1)*(4)+1-1]= 2.; E_info.node_coor[(2-1)*(4)+1-1]= 0.;
	E_info.node_coor[(1-1)*(4)+2-1]= 2.; E_info.node_coor[(2-1)*(4)+2-1]= 1.;
	E_info.node_coor[(1-1)*(4)+3-1]= 1.; E_info.node_coor[(2-1)*(4)+3-1]= 1.;
	E_info.node_coor[(1-1)*(4)+4-1]= 1.; E_info.node_coor[(2-1)*(4)+4-1]= 0.;
	
	E_info.elem_mate = (double*) malloc (2 * sizeof(double));
	E_info.elem_mate[0] = 8.8541878e-12;
	E_info.elem_mate[1] = 0.0;


	int node_dof = 1;
	int elem_dof = E_info.node_cont * node_dof;
	E_matr.matr_0    = (double*) calloc (elem_dof * elem_dof, sizeof(double));
	E_matr.matr_1    = (double*) calloc (elem_dof * elem_dof, sizeof(double));
	E_matr.matr_2    = (double*) calloc (elem_dof * elem_dof, sizeof(double));
	E_matr.left_matr = (double*) calloc (elem_dof * elem_dof, sizeof(double));
	E_matr.righ_vect = (double*) calloc (elem_dof, sizeof(double));
	
	E_info.refr_shap = (double**) malloc (G_info.gaus_num * sizeof(double*));
	for (int gaus_i = 1; gaus_i <= G_info.gaus_num; gaus_i ++)
		E_info.refr_shap[gaus_i-1] = (double*) malloc (E_info.node_cont * (E_info.dim+1) * sizeof(double));
	
	set_refr_shap(E_info.refr_shap, G_info.gaus_coor, G_info.gaus_num, E_info.node_cont, E_info.dim);

	//for (int i=0; i<G_info.gaus_num; i++){
    //    for (int j=0; j<E_info.node_cont * (E_info.dim+1); j++)
    //        printf("%e ",E_info.refr_shap[i][j]);
    //    printf("\n");
    //}
	
	double coup_valu[]={0.,0.};
	int elem_n = 2;
	show_elem(E_info, 2, 0, 2, 4);

	elemcalc(elem_n, G_info, E_info, &E_matr);

	//for (int i=0; i<elem_dof*elem_dof; i++)
	//	printf("%e\n",E_matr.matr_0[i]);
//
	//for (int i=0; i<elem_dof; i++)
	//	printf("%e\n",E_matr.righ_vect[i]);
	
	
	
	
	printf("done!\n");
	
	return 1;
}
void show_elem(Elem_Info E_info, int dim, int elem_nodeN, int mate_varN, int gaus_num)
{
	printf("dim: %d\n",E_info.dim);

	printf("node_cont: %d\n",E_info.node_cont);

	printf("elem_node:\n");
	for (int i=0; i<elem_nodeN; i++)
		printf("%d ",E_info.elem_node[i]);
	printf("\n");

	printf("node_coor:\n");
	for (int i=0; i<E_info.node_cont; i++){
		for (int j=0; j<dim; j++)
			printf("%e ",E_info.node_coor[j*E_info.node_cont+i]);
		printf("\n");
	}

	printf("elem_mate:\n");
	for (int i=0; i<mate_varN; i++)
		printf("%e ",E_info.elem_mate[i]);
	printf("\n");

	printf("refr_shap:\n");
	for (int i=0; i<gaus_num; i++){
		for (int j=0; j<E_info.node_cont*(dim+1); j++)
			printf("%e ",E_info.refr_shap[i][j]);
		printf("\n");
	}

}