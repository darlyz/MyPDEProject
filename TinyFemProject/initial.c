#include "fem.h"

#define init_adj_num 20

void initial(
    Coor_Info Coor,
    Node_Mesh Mesh,
    Equat_Set Equa,
    Dof_Tag   ID,
    int node_dof
    
){

    int  *adj_nodn;
    int **adj_topo;

    adj_nodn = (int* )calloc(Coor.total_nodes,sizeof(int));
    adj_topo = (int**)calloc(Coor.total_nodes,sizeof(int*));
    for (int i=0; i<Coor.total_nodes; i++)
        adj_topo[i] = (int*)malloc(init_adj_num*sizeof(int));

    for (int type_i=0; type_i<Mesh.typeN; type_i++) {

        int *node;
        int nodeN = Mesh.elem_nodeN[type_i] -1;
        node = (int*)malloc(nodeN*sizeof(int));

        for (int elem_i=0; elem_i<Mesh.mesh_scale[type_i]; elem_i++) {

            memcpy(node, &Mesh.mesh_topo[type_i][elem_i*Mesh.elem_nodeN[type_i]], nodeN*sizeof(int));

            for (int node_i; node_i<nodeN; node_i++) {

                int nodeSN = node[node_i];

                for (int node_j=0; node_j<nodeN; node_j++) {

                    if (node_i == node_j)
                        continue;

                    int insert_SN = node[node_j];

                    insert_node(nodeSN, insert_SN, adj_nodn, adj_topo);
                }
            }
        }

        free(node);
    }


    printf("initial done!\n");
}

int Binary_Search(int* dest, int dest_len, int key)
{
    int idx_start = 0;
    int idx_end   = dest_len-1;
    int idx_mid   = -1;
    int insert    = 0;

    if (dest_len == 0)
    	return 0;

    while(1)
    {
        idx_mid = (idx_start + idx_end)/2;

        if      (key < dest[idx_mid])
            idx_end   = idx_mid - 1;

        else if (key > dest[idx_mid])
            idx_start = idx_mid + 1;

        else
        	return -1;

        if (idx_start >= idx_end) {

        	if (idx_mid == idx_start)
        		return idx_mid;

        	else if (idx_mid == idx_end)
        		return idx_mid+1;
        }
    }
}

void insert_node(int nodeSN, int insert_SN, int* adj_nodn, int** adj_topo)
{
    if (adj_nodn[nodeSN] == 0) {

        adj_topo[nodeSN][0] = insert_SN;
        adj_nodn[nodeSN] ++;
    }

    else {

        int insert_idx = Binary_Search(adj_topo[nodeSN], adj_nodn[nodeSN], insert_SN);

        if (insert_SN != -1) {

            if (adj_nodn[nodeSN] % init_adj_num == 0) {

                int *temp_topo;
                temp_topo = (int*)calloc( (adj_nodn[nodeSN]/init_adj_num + 1)*init_adj_num*sizeof(int) );
                
                memcpy(temp_topo, adj_topo[nodeSN], insert_idx*sizeof(int));
                temp_topo[insert_idx] = insert_SN;
                memcpy(temp_topo + insert_SN + 1,
                       adj_topo[nodeSN] + insert_SN,
                      (adj_nodn[nodeSN] - insert_idx)*sizeof(int) );
                
                free(adj_topo[nodeSN]);
                adj_topo[nodeSN] = temp_topo;
                temp_topo = NULL;

            }

            else {

                for (int i=adj_nodn[nodeSN]; i>insert_idx; i--)
                    adj_topo[nodeSN][i] = adj_topo[nodeSN][i-1];

                adj_topo[nodeSN][insert_idx] = insert_SN;
            }

            adj_nodn[nodeSN] ++;
        }
    }
}

int compare_int(const void *value1, const void *value2) 
{
    return *(int*)value1 - *(int*)value2;
}

void int_qsort(int* array, int array_len)
{
    qsort(array, array_len, sizeof(int), compare_int);
}