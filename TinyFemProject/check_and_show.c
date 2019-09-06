#include "fem.h"

// used it in elemcalc.c
void show_elem_stif(int node_cont, Elem_Matr E_matr) {
	for (int i=0; i<node_cont*node_cont; i++)
	    printf("%e\n",E_matr.matr_0[i]);
    
	for (int i=0; i<node_cont; i++)
	    printf("%e\n",E_matr.righ_vect[i]);

	printf("\n");
}

// used it in initial.c
void show_adj(int* adj_nodn, int** adj_topo, int total_nodes) {
    for (int i=0; i<total_nodes; i++) {
        printf("%d, %d : ",i+1,adj_nodn[i]);
        for (int j=0; j<adj_nodn[i]; j++)
            printf("%d ",adj_topo[i][j]);
        printf("\n");
    }
}

// used it in initial.c
void show_node_eq_index(Equat_Set Equa, int total_nodes, int node_dof) {
    printf("total_equations = %d\n",Equa.total_equations);
    for (int i=0; i<total_nodes; i++) {
        for (int j=0; j<node_dof; j++)
            printf("\t%d",Equa.node_equa_index[i][j]);
        printf("\n");
    }
}

// used it in initial.c
void show_non_trivial(Equat_Set Equa) {
    printf("non triaval: %d\n",Equa.total_nontriaval);
    for (int i=0; i<Equa.total_equations; i++) {
        printf("%d, ",i+1);
        for (int j=0; j<Equa.row_nontriaval[i]; j++)
            printf("%d ",Equa.column_index[i][j]);
        printf("\n");
    }
}

// used it in matrcalc.c
void show_elem(Elem_Info E_info, int elem_nodeN, int mate_varN, int gaus_num) {

    int dim = E_info.dim;

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

// used it in matrcalc.c
void show_matr(Equat_Set Equa) {
    for (int i=0; i<Equa.total_equations; i++){
        printf("%d: ",i+1);
        for (int j=0; j<Equa.row_nontriaval[i]; j++)
            printf("%d,%e ",Equa.column_index[i][j],Equa.matrix[i][j]);
        printf("-->%e\n",Equa.vector[i]);
    }
}

// used it in matrcalc.c
void show_elem_matr(Elem_Matr E_matr, int ematr_size, Matr_Type *M_type) {
    int size[4];
    for (int i=0; i<4; i++){
        if     (M_type[i] == lump) size[i] = ematr_size;
        else if(M_type[i] == dist) size[i] = ematr_size*ematr_size;
    }

    printf("\nmatr0:\n");
    for (int i=0; i<size[0]; i++)
        printf("%e\n",E_matr.matr_0[i]);

    printf("\nforc:\n");
    for (int i=0; i<ematr_size; i++)
        printf("%e\n",E_matr.righ_vect[i]);
}

void show_mesh(Node_Mesh Mesh) {
    for (int i=0; i<Mesh.typeN; i++) {
        printf("type_i=%d\n",i+1);
        for(int j=0; j<Mesh.mesh_scale[i];j++) {
            printf("elem_i=%d: ",j+1);
            for (int k=0; k<Mesh.elem_nodeN[i];k++)
                printf("%d ",Mesh.mesh_topo[i][j*Mesh.elem_nodeN[i]+k]);
            printf("\n");
        }
        printf("\n");
    }
}