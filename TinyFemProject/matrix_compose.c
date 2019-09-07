/*
 Copyright: Copyright (c) 2019
 Created: 2019-4-19
 Author: Zhang_Licheng
 All rights reserved
*/
#include "fem.h"

void set_gaus  ();
void set_matr  ();
void set_elem  ();
void reset_matr();
void show_elem ();
void show_matr ();
void show_elem_matr();
void clear_matr();
void clear_elem();
void set_refr_shap();
void elemcalc();
void show_elem_stif();
int  Binary_Search_();
int  traversal_search_();

void matrcalc(
    Coor_Info  Coor,
    Node_Mesh  Mesh,
	Dof_Tag    ID,
	Equat_Set* Equa,
    Materail   Mate,
    double *result,
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

            if (E_info.elem_node[E_info.node_cont] <= 0)
                continue;

			memcpy(E_info.elem_mate, &Mate.mate[E_info.elem_node[E_info.node_cont]-1], Mate.mate_varN*sizeof(double));
            
			for (int node_i=1; node_i<elem_nodeN; node_i++) {
                for (int dim_i=1; dim_i<=E_info.dim; dim_i++)

                    E_info.node_coor[(dim_i-1)*(elem_nodeN-1)+node_i-1] = 
                    Coor.coordinate[(E_info.elem_node[node_i-1]-1)*E_info.dim + dim_i-1];
            }
            
            reset_matr(&E_matr, ematr_size, M_type);

			elemcalc(elem_i, G_info, E_info, &E_matr);

            //show_elem(E_info, elem_nodeN, Mate.mate_varN, G_info.gaus_num);
            //show_elem_stif(E_info.node_cont, E_matr);
            //show_elem_matr(E_matr, ematr_size, M_type);

            // construct element left matrix
            for (int dof_i=0; dof_i<ematr_size; dof_i++) {
                
                for (int dof_j=0; dof_j<ematr_size; dof_j++) {
                    // distributed matrix
                    E_matr.left_matr[dof_i*ematr_size + dof_j] += E_matr.matr_0[dof_i*ematr_size + dof_j];
                }
                // lumped matrix
                //E_matr.left_matr[dof_i*ematr_size + dof_i] += E_matr.matr_1[dof_i]*dt/2
                //                                           +  E_matr.matr_2[dof_i];
            }

            // construct total left matrix and right vector
            for (int nod_i=0; nod_i<elem_nodeN-1; nod_i++) {
                
                int nodi_SN = E_info.elem_node[nod_i];

                for (int dof_i=0; dof_i<node_dof; dof_i++) {

                    int ID_i = Equa->node_equa_index[nodi_SN-1][dof_i];

                    int MD_i = nod_i*node_dof + dof_i;

                    if (ID_i == 0)
                        continue;

                    if (ID_i  > 0) {

                        // lumped matrix for right vector
                        Equa->vector[ID_i-1] += E_matr.righ_vect[MD_i];
                        //                   +  E_matr.righ_vect[MD_i] * dt*dt/2
                        //                   +  E_matr.matr_2   [MD_i] * u_n[(nodi_SN-1)*node_dof + dof_i]
                        //                   +  E_matr.matr_2   [MD_i] * v_n[(nodi_SN-1)*node_dof + dof_i] * dt
                        //                   +  E_matr.matr_1   [MD_i] * u_n[(nodi_SN-1)*node_dof + dof_i] * dt/2;
                    }

                    int idx = 0;

                    for (int nod_j=0; nod_j<elem_nodeN-1; nod_j++) {

                        int nodj_SN = E_info.elem_node[nod_j];

                        for (int dof_j=0; dof_j<node_dof; dof_j++) {

                            int ID_j = Equa->node_equa_index[nodj_SN-1][dof_j];

                            int MD_j = nod_j*node_dof + dof_j;

                            if (ID_j == 0)
                                continue;

                            if (ID_j  > 0) {

                                // distributed matrix for right vector
                                //Equa->vector[ID_j-1] += -E_matr.matr_0[MD_i*ematr_size + MD_j] * u_n[(nodi_SN-1)*node_dof + dof_i] *dt*dt/4
                                //                     +   E_matr.matr_0[MD_i*ematr_size + MD_j] * v_n[(nodi_SN-1)*node_dof + dof_i] *dt*dt/2;
                            
                                // deal with direclet boundary conditions
                                if (ID_i < 0)
                                    Equa->vector[ID_j-1] -= E_matr.left_matr[MD_i*ematr_size + MD_j] * result[(nodi_SN-1)*node_dof+dof_i];

                                // compose whole matrix
                                else if (ID_i > 0) {

                                    idx = Binary_Search_(Equa->column_index[ID_i-1], Equa->row_nontriaval[ID_i-1], ID_j);

                                    if (idx <= 0)
                                        continue;

                                    Equa->matrix[ID_i-1][idx] += E_matr.left_matr[MD_i*ematr_size + MD_j];
                                    //idx ++;
                                }
                            }
                        }
                    }
                }
            }
        }
        clear_matr(&E_matr);
        clear_elem(&E_info, G_info.gaus_num);
    }

    //show_matr(*Equa);
    printf("Compose whole matrix done!\n");
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