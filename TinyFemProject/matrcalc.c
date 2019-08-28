/*
 Copyright: Copyright (c) 2019
 Created: 2019-4-19
 Author: Zhang_Licheng
 All rights reserved
*/
#include "fem.h"

void set_gaus  (Gaus_Info*);
void set_matr  (Elem_Matr*, int, Matr_Type*);
void set_elem  (Elem_Info*, int, int, int, int);
void reset_matr(Elem_Matr*, int, Matr_Type*);
void show_elem (Elem_Info,  int, int, int, int);
void clear_matr(Elem_Matr*);
void clear_elem(Elem_Info*, int);
void set_refr_shap(double**, double*, int, int, int);
void elemcalc(int, Gaus_Info, Elem_Info, Elem_Matr*);

void matrcalc(
    Coor_Info  Coor,
    Node_Mesh  Mesh,
	Dof_Tag    ID,
	Equat_Set* Equa,
    Materail   Mate,
    int node_dof
){
    int pre_node = 0;
    for (int type_i=1; type_i<=Mesh.typeN; type_i++)
    {
        int elem_nodeN = Mesh.elem_nodeN[type_i-1];
        int mesh_scale = Mesh.mesh_scale[type_i-1];
        int ematr_size = node_dof*(elem_nodeN-1);

        Matr_Type M_type[4]={dist,lump,lump,lump};
		
		Gaus_Info G_info;
		Elem_Info E_info;
		Elem_Matr E_matr;
		
		set_gaus(&G_info);
		set_matr(&E_matr, ematr_size, M_type);
		set_elem(&E_info, Coor.dim, elem_nodeN, Mate.mate_varN, G_info.gaus_num);
        set_refr_shap(E_info.refr_shap, G_info.gaus_coor, G_info.gaus_num, E_info.node_cont, E_info.dim);

        for (int elem_i=1; elem_i<=mesh_scale; elem_i++)
        {

            memcpy(E_info.elem_node, &Mesh.mesh_topo[type_i-1][(elem_i-1)*elem_nodeN], elem_nodeN*sizeof(int));
			memcpy(E_info.elem_mate, &Mate.mate[E_info.elem_node[E_info.node_cont]-1], Mate.mate_varN*sizeof(double));
            
			for (int node_i=1; node_i<elem_nodeN; node_i++) {
                for (int dim_i=1; dim_i<=E_info.dim; dim_i++)

                    E_info.node_coor[(dim_i-1)*(elem_nodeN-1)+node_i-1] = 
                    Coor.coordinate[(E_info.elem_node[node_i-1]-1)*E_info.dim + dim_i-1];
            }
            
            reset_matr(&E_matr, ematr_size, M_type);

			elemcalc(elem_i, G_info, E_info, &E_matr);
        }

        clear_matr(&E_matr);
        clear_elem(&E_info, G_info.gaus_num);
    }
}

void set_matr(Elem_Matr* E_matr, int elem_dof, Matr_Type *M_type)
{
    int size[4];
    for (int i=0; i<4; i++) {
        if     (M_type[i] == lump) size[i] = elem_dof;
        else if(M_type[i] == dist) size[i] = elem_dof*elem_dof;
    }
    E_matr->matr_0    = (double*)malloc(size[0]          *sizeof(double));
    E_matr->matr_1    = (double*)malloc(size[1]          *sizeof(double));
    E_matr->matr_2    = (double*)malloc(size[2]          *sizeof(double));
    E_matr->matr_3    = (double*)malloc(size[3]          *sizeof(double));
    E_matr->left_matr = (double*)malloc(elem_dof*elem_dof*sizeof(double));
    E_matr->righ_vect = (double*)malloc(elem_dof         *sizeof(double));
}

void reset_matr(Elem_Matr* E_matr, int elem_dof, Matr_Type *M_type)
{
    int size[4];
    for (int i=0; i<4; i++){
        if     (M_type[i] == lump) size[i] = elem_dof;
        else if(M_type[i] == dist) size[i] = elem_dof*elem_dof;
    }
    memset(E_matr->matr_0    , 0, size[0]          *sizeof(double));
    memset(E_matr->matr_1    , 0, size[1]          *sizeof(double));
    memset(E_matr->matr_2    , 0, size[2]          *sizeof(double));
    memset(E_matr->matr_3    , 0, size[3]          *sizeof(double));
    memset(E_matr->left_matr , 0, elem_dof*elem_dof*sizeof(double));
    memset(E_matr->righ_vect , 0, elem_dof         *sizeof(double));
}

    
void clear_matr(Elem_Matr* E_matr)
{
    free(E_matr->matr_0);
    free(E_matr->matr_1);
    free(E_matr->matr_2);
    free(E_matr->matr_3);
    free(E_matr->left_matr);
    free(E_matr->righ_vect);
}

void set_elem(Elem_Info *E_info, int dim, int elem_nodeN, int mate_varN, int gaus_num)
{
	E_info->dim = dim;
	E_info->node_cont  = elem_nodeN - 1;
	
	E_info->elem_node = (int*    )malloc(elem_nodeN * sizeof(int));
	E_info->node_coor = (double* )malloc(E_info->node_cont * E_info->dim * sizeof(double));
	//E_info->coup_valu = (double*)malloc(E_info->node_cont * coupN * sizeof(double));
	E_info->elem_mate = (double* )malloc(mate_varN * sizeof(double));
	E_info->refr_shap = (double**)malloc(gaus_num * sizeof(double*));
	for (int gaus_i=1; gaus_i<=gaus_num; gaus_i++)
		E_info->refr_shap[gaus_i-1] = (double*)malloc(E_info->node_cont * (E_info->dim+1) * sizeof(double));
}

void clear_elem(Elem_Info *E_info, int gaus_num)
{
	free(E_info->elem_node);
    free(E_info->node_coor);
    //free(E_info->coup_valu);
    free(E_info->elem_mate);
	for (int gaus_i=1; gaus_i<=gaus_num; gaus_i++)
		free(E_info->refr_shap[gaus_i-1]);
	free(E_info->refr_shap);
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