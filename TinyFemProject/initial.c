#include "fem.h"

#define init_adj_num 20

void insert_node(int, int, int*, int**);
void show_adj(int*, int**, int);

void initial(
    Coor_Info Coor,
    Node_Mesh Mesh,
    Equat_Set Equa,
    Dof_Tag   ID,
    int node_dof
    
){

    int   nodeN;
    int  *node;
    int  *node_SN;
    int  *insert_SN;
    int  *adj_nodn;
    int **adj_topo;

    adj_nodn = (int* )calloc(Coor.total_nodes,sizeof(int));
    adj_topo = (int**)calloc(Coor.total_nodes,sizeof(int*));
    for (int i=0; i<Coor.total_nodes; i++)
        adj_topo[i] = (int*)malloc(init_adj_num*sizeof(int));

    for (int type_i=0; type_i<Mesh.typeN; type_i++) {

        nodeN = Mesh.elem_nodeN[type_i] -1;
        //node = (int*)malloc(nodeN*sizeof(int));

        for (int elem_i=0; elem_i<Mesh.mesh_scale[type_i]; elem_i++) {

            //memcpy(node, &Mesh.mesh_topo[type_i][elem_i*Mesh.elem_nodeN[type_i]], nodeN*sizeof(int));
            node = &Mesh.mesh_topo[type_i][elem_i*Mesh.elem_nodeN[type_i]];

            for (int node_i=0; node_i<nodeN; node_i++) {

                node_SN = &node[node_i];

                for (int node_j=0; node_j<nodeN; node_j++) {

                    if (node_i == node_j)
                        continue;

                    insert_SN = &node[node_j];

                    insert_node(*node_SN, *insert_SN, adj_nodn, adj_topo);
                }
            }
        }
    }

    show_adj(adj_nodn, adj_topo, Coor.total_nodes);

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

int search(int* dest, int dest_len, int key)
{
    for (int i=0; i<dest_len; i++) {
        if ( dest[0] > key )
            return 0;

        else if ( dest[dest_len-1] < key )
            return dest_len;

        else if ( dest[i] < key && dest[i+1] > key)
            return i+1;
    }

    return -1;
}

void insert_node(int node_SN, int insert_SN, int* adj_nodn, int** adj_topo)
{
    node_SN --;

    if (adj_nodn[node_SN] == 0) {

        adj_topo[node_SN][0] = insert_SN;
        adj_nodn[node_SN] ++;
    }

    else {

        int insert_idx = Binary_Search(adj_topo[node_SN], adj_nodn[node_SN], insert_SN);

        if (insert_idx != -1) {

            if (adj_nodn[node_SN] % init_adj_num == 0) {

                int *temp_topo;
                temp_topo = (int*)calloc( (adj_nodn[node_SN]/init_adj_num + 1)*init_adj_num, sizeof(int) );
                
                memcpy(temp_topo, adj_topo[node_SN], insert_idx*sizeof(int));
                temp_topo[insert_idx] = insert_SN;
                memcpy(temp_topo + insert_SN + 1,
                       adj_topo[node_SN] + insert_SN,
                      (adj_nodn[node_SN] - insert_idx)*sizeof(int) );
                
                free(adj_topo[node_SN]);
                adj_topo[node_SN] = temp_topo;
                temp_topo = NULL;

            }

            else {

                for (int i=adj_nodn[node_SN]; i>insert_idx; i--)
                    adj_topo[node_SN][i] = adj_topo[node_SN][i-1];

                adj_topo[node_SN][insert_idx] = insert_SN;
            }

            adj_nodn[node_SN] ++;
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

void show_adj(int* adj_nodn, int** adj_topo, int total_nodes)
{
    for (int i=0; i<total_nodes; i++) {
        printf("%d, ",i+1);
        for (int j=0; j<adj_nodn[i]; j++)
            printf("%d ",adj_topo[i][j]);
        printf("\n");
    }
}