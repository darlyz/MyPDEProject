#include "fem.h"

#define init_adj_num 20
#define minsize_int32 -2147483648

void free_adj();
void int_qsort();
void insert_node();
void insert_node_();
void show_adj();
void show_node_eq_index();
void show_non_trivial();

void initial(
    Coor_Info Coor,
    Node_Mesh Mesh,
    Dof_Tag   ID,
    Equat_Set *Equa,
    double **result,
    int node_dof
){
    // ---------------------------------------- initial solution space ----------------------------------------
    (*result) = (double *)calloc(Coor.total_nodes*node_dof, sizeof(int));
    for (int i=0; i<ID.tag_nodn; i++) {
        int node_SN = ID.tag_nods[i] - 1;
        for (int j=0; j<ID.dof_num; j++)
            if (ID.dof_tag[i*ID.dof_num + j] == -1 )
                (*result)[node_SN*ID.dof_num + j] = ID.dof_val[i*ID.dof_num + j];
    }

    // ---------------------------------------- convert mesh to graph ----------------------------------------
    int   nodeN;
    int  *node;
    int  *adj_nodn;
    int **adj_topo;

    adj_nodn = (int* )calloc(Coor.total_nodes,sizeof(int));
    adj_topo = (int**)calloc(Coor.total_nodes,sizeof(int*));
    for (int i=0; i<Coor.total_nodes; i++)
        adj_topo[i] = (int*)malloc(init_adj_num*sizeof(int));

    for (int type_i=0; type_i<Mesh.typeN; type_i++) {

        nodeN = Mesh.elem_nodeN[type_i] -1;

        for (int elem_i=0; elem_i<Mesh.mesh_scale[type_i]; elem_i++) {

            node = &Mesh.mesh_topo[type_i][elem_i*Mesh.elem_nodeN[type_i]];

            for (int node_i=0; node_i<nodeN; node_i++) {

                for (int node_j=0; node_j<nodeN; node_j++) {

                    if (node_i == node_j)
                        continue;

                    insert_node(adj_topo + node[node_i] - 1, adj_nodn + node[node_i] - 1, init_adj_num, node[node_j]);
                    //insert_node_(node[node_i], node[node_j], adj_nodn, adj_topo);
                }
            }
        }
    }

    //show_adj(adj_nodn, adj_topo, Coor.total_nodes);

    // --------------------------------------- initial Equation Set ---------------------------------------
    //int constraint_count = 0;
    //for (int i=0; i<ID.tag_nodn; i++)
    //    for (int j=0; j<ID.dof_num; j++)
    //        if (ID.dof_tag[i*ID.dof_num + j] <= 0)
    //            constraint_count ++ ;

    //Equa->total_equations = Coor.total_nodes * node_dof - constraint_count;
    //printf("constraint_count = %d\n",constraint_count);

    int min_int32 = -2147483648;

    Equa->node_equa_index = (int**)malloc(Coor.total_nodes * sizeof(int*));
    for (int i=0; i<Coor.total_nodes; i++) {
        Equa->node_equa_index[i] = (int*)malloc(node_dof*sizeof(int));
        memcpy(Equa->node_equa_index[i], &min_int32, node_dof*sizeof(int));
    }

    for (int i=0; i<ID.tag_nodn; i++) 
        for (int j=0; j<ID.dof_num; j++)
            Equa->node_equa_index[ID.tag_nods[i]-1][j] = ID.dof_tag[i*ID.dof_num + j];

    int equation_index = 0;
    for (int i=0; i<Coor.total_nodes; i++)
        for (int j=0; j<node_dof; j++)
            if (Equa->node_equa_index[i][j] == min_int32)
                Equa->node_equa_index[i][j] = (++equation_index);

    Equa->total_equations = equation_index;
    //show_node_eq_index(*Equa, Coor.total_nodes, node_dof);
    //printf("equation_index = %d\n",equation_index);
    
    Equa->row_nontriaval = (int* )malloc(Equa->total_equations * sizeof(int));
    Equa->column_index   = (int**)malloc(Equa->total_equations * sizeof(int*));
    Equa->matrix      = (double**)malloc(Equa->total_equations * sizeof(double*));
    Equa->vector      = (double* )malloc(Equa->total_equations * sizeof(double));

    for (int i=0; i<Coor.total_nodes; i++)
        for (int j=0; j<node_dof; j++) {

            int row_index   = Equa->node_equa_index[i][j];

            if (row_index > 0) {
                Equa->column_index[row_index-1] = (int*   )malloc((adj_nodn[i]+1)*node_dof * sizeof(int));
                Equa->matrix      [row_index-1] = (double*)malloc((adj_nodn[i]+1)*node_dof * sizeof(double));
            }
        }


    // --------------------------------------- convert graph to matrix ---------------------------------------
    Equa->total_nontriaval = 0;
    for (int i=0; i<Coor.total_nodes; i++)

        for (int j=0; j<node_dof; j++) {

            int non_trivial = 0;
            int row_index   = Equa->node_equa_index[i][j];

            if (row_index > 0) {

                for (int jj=0; jj<node_dof; jj++) {

                    int clm_index = Equa->node_equa_index[i][jj];

                    if (clm_index >0)
                        Equa->column_index[row_index-1][non_trivial++] = clm_index;
                    
                    else
                        continue;
                }

                for (int ii=0; ii<adj_nodn[i]; ii++)

                    for (int jj=0; jj<node_dof; jj++) {

                        int clm_index = Equa->node_equa_index[ adj_topo[i][ii] -1 ][ jj ];

                        if (clm_index > 0)
                            Equa->column_index[row_index-1][non_trivial++] = clm_index;

                        else
                            continue;
                    }
            }

            else
                continue;
            
            Equa->row_nontriaval[row_index-1] = non_trivial;
            Equa->total_nontriaval += non_trivial;
        }


    for (int i=0; i<Equa->total_equations; i++)
        int_qsort(Equa->column_index[i], Equa->row_nontriaval[i]);

    //show_non_trivial(*Equa);

    free_adj(Coor.total_nodes, adj_nodn, adj_topo);

    printf("Initial done!\n");
}

void free_adj(int total_nodes, int* adj_nodn, int** adj_topo) {
    for (int i=0; i<total_nodes; i++)
        free(adj_topo[i]);
    free(adj_topo);
    free(adj_nodn);
}