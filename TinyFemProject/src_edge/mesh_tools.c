#include "fem.h"

#define init_adj_num 20

void insert_node();
void insert_node_();
void remove_node();
int Binary_Search_();

void mesh2graph(
    int  *adj_nodn,
    int **adj_topo,
    Node_Mesh Mesh
){
    int   nodeN;
    int  *node;

    for (int type_i=0; type_i<Mesh.typeN; type_i++) {

        nodeN = Mesh.nodeN[type_i];

        for (int elem_i=0; elem_i<Mesh.scale[type_i]; elem_i++) {

            node = &(Mesh.topo[type_i][elem_i*nodeN]);

            for (int node_i=0; node_i<nodeN; node_i++) {

                for (int node_j=0; node_j<nodeN; node_j++) {

                    if (node_i == node_j)
                        continue;

                    insert_node(adj_topo + node[node_i] - 1,
                                adj_nodn + node[node_i] - 1,
                                init_adj_num,
                                node[node_j] );
                    //insert_node_(node[node_i], node[node_j], adj_nodn, adj_topo);
                }
            }
        }
    }
}

void graph_contraction(
    int   *adj_nodn,
    int  **adj_topo,
    int  **adj_nodn_ctr,
    int ***adj_topo_ctr,
    int total_nodes
){
    int  *temp_adj_nodn = (int* )malloc(total_nodes*sizeof(int));
    int **temp_adj_topo = (int**)malloc(total_nodes*sizeof(int*));

    memcpy(temp_adj_nodn, adj_nodn, total_nodes*sizeof(int));
    
    for (int i=0; i<total_nodes; i++) {
        temp_adj_topo[i] = (int*)malloc(adj_nodn[i]*sizeof(int));
        memcpy(temp_adj_topo[i], adj_topo[i], adj_nodn[i]*sizeof(int));
    }

    int nodeSN;

    for (int i=1; i<total_nodes; i++) {
        
        int j = 0;

        while (j < temp_adj_nodn[i]) {

            nodeSN = temp_adj_topo[i][j];

            j++;

            if (nodeSN < i+1) {

                if ( Binary_Search_(temp_adj_topo[nodeSN-1], temp_adj_nodn[nodeSN-1], i+1) != -1) {
                    
                    remove_node(&(temp_adj_topo[i]), &(temp_adj_nodn[i]), nodeSN);

                    j--;
                }
            }
        }
    }

    (*adj_nodn_ctr) = (int* )malloc(total_nodes*sizeof(int));
    (*adj_topo_ctr) = (int**)malloc(total_nodes*sizeof(int*));

    for (int i = 0; i < total_nodes; i++) {

        if (temp_adj_nodn[i] != 0) {

            (*adj_topo_ctr)[i] = (int*)malloc(temp_adj_nodn[i]*sizeof(int));

            memcpy((*adj_topo_ctr)[i], temp_adj_topo[i], temp_adj_nodn[i]*sizeof(int));
        }

        if (i > 0)
            temp_adj_nodn[i] += temp_adj_nodn[i-1];

        (*adj_nodn_ctr)[i] = temp_adj_nodn[i];
    }

    for (int i = 0; i <total_nodes; i++)
        free(temp_adj_topo[i]);
    free(temp_adj_topo);
    free(temp_adj_nodn);
}

void edge_census(
    int  *adj_nodn,
    int **adj_topo,
    int **edge_nodes,
    int *edgeN,
    int total_nodes
){
    (*edgeN) = adj_nodn[total_nodes-1];
    (*edge_nodes) = (int*)malloc(2*(*edgeN)*sizeof(int));

    int edge_i = 0;
    int nodnum = 0;
    for (int i = 0; i < total_nodes; i++) {

        nodnum = (i>0) ? (adj_nodn[i]-adj_nodn[i-1]) : (adj_nodn[i]);

        for (int j = 0; j <nodnum; j++) {

            (*edge_nodes)[edge_i*2 + 0] = i+1;
            (*edge_nodes)[edge_i*2 + 1] = adj_topo[i][j];
            edge_i ++ ;
        }
    }
}