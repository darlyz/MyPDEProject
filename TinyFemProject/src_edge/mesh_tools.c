#include "fem.h"

#define init_adj_num 20

void insert_node();
void insert_node_();
void remove_node();
int Binary_Search_();

void Mesh2Graph_Neighbour(
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
                }
            }
        }
    }
}

int CalcEMeshEdgeN(Mesh_Type type) {
    switch (type) {
        case P1 : return  0; break;
        case L2 : return  1; break;
        case L3 : return  2; break;
        case T3 : return  3; break;
        case T6 : return  6; break;
        case Q4 : return  4; break;
        case Q8 : return  8; break;
        case Q9 : return  8; break;
        case W4 : return  6; break;
        case W10: return 12; break;
        case C8 : return 12; break;
        case C20: return 24; break;
        case C27: return 24; break;
        case H6 : return  9; break;
        case H15: return 18; break;
        case H18: return 18; break;
        case P5 : return  8; break;
        case P13: return 16; break;
        case P14: return 16; break;
        default : return  0;
    }
}

int LocalEMeshTopo(Mesh_Type type, int **local_topo) {

    int edgeN = CalcEMeshEdgeN(type);

    (*local_topo) = (int*)malloc(2*edgeN*sizeof(int));
    switch (type) {

        case L2:
            (*local_topo)[0] = 1; (*local_topo)[1] = 2;
            break;

        case L3:
            (*local_topo)[0] = 2; (*local_topo)[1] = 1;
            (*local_topo)[2] = 2; (*local_topo)[3] = 3;
            break;

        case T3:
            (*local_topo)[0] = 1; (*local_topo)[1] = 2;
            (*local_topo)[2] = 2; (*local_topo)[3] = 3;
            (*local_topo)[4] = 3; (*local_topo)[5] = 1;
            break;


        case Q4:
            (*local_topo)[0] = 1; (*local_topo)[1] = 2;
            (*local_topo)[2] = 4; (*local_topo)[3] = 3;
            (*local_topo)[4] = 1; (*local_topo)[5] = 4;
            (*local_topo)[6] = 2; (*local_topo)[7] = 3;
            break;

        
        case W4:
            (*local_topo)[ 0] = 1; (*local_topo)[ 1] = 2;
            (*local_topo)[ 2] = 1; (*local_topo)[ 3] = 3;
            (*local_topo)[ 4] = 1; (*local_topo)[ 5] = 4;
            (*local_topo)[ 6] = 2; (*local_topo)[ 7] = 3;
            (*local_topo)[ 8] = 4; (*local_topo)[ 9] = 2;
            (*local_topo)[10] = 3; (*local_topo)[11] = 4;
            break;

    }

    return edgeN;
}

void Mesh2Graph_Connected(
    int  *adj_nodn,
    int **adj_topo,
    Node_Mesh Mesh
){
    int   nodeN;
    int  *node;

    for (int type_i=0; type_i<Mesh.typeN; type_i++) {

        nodeN = Mesh.nodeN[type_i];

        int *local_topo;
        int edgeN = LocalEMeshTopo(Mesh.type[type_i], &local_topo);

        for (int elem_i=0; elem_i<Mesh.scale[type_i]; elem_i++) {

            node = &(Mesh.topo[type_i][elem_i*nodeN]);

            for (int edge_i=0; edge_i<edgeN; edge_i++) {

                insert_node(adj_topo + node[local_topo[edge_i*2+0] - 1] - 1,
                            adj_nodn + node[local_topo[edge_i*2+0] - 1] - 1,
                            init_adj_num,
                            node[local_topo[edge_i*2+1] - 1]);

                insert_node(adj_topo + node[local_topo[edge_i*2+1] - 1] - 1,
                            adj_nodn + node[local_topo[edge_i*2+1] - 1] - 1,
                            init_adj_num,
                            node[local_topo[edge_i*2+0] - 1]);
            }
        }
    }
}

void GraphContraction(
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

void EdgeCensus(
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

void NodeMesh2Edge(
    Node_Mesh  NMesh,
    Edge_Mesh *EMesh,
    Edge_Info  Edges
){
    int typeN = NMesh.typeN;

    EMesh->typeN = typeN;
    EMesh->l_dim = (int* )malloc(typeN * sizeof(int));
    EMesh->edgeN = (int* )malloc(typeN * sizeof(int));
    EMesh->scale = (int* )malloc(typeN * sizeof(int));
    EMesh->topo  = (int**)malloc(typeN * sizeof(int));
    EMesh->drict = (int**)malloc(typeN * sizeof(int));
    EMesh->type  = (int* )malloc(typeN * sizeof(Mesh_Type));

    memcpy(EMesh->l_dim, NMesh.l_dim, typeN*sizeof(int));
    memcpy(EMesh->scale, NMesh.scale, typeN*sizeof(int));
    memcpy(EMesh->type , NMesh.type , typeN*sizeof(Mesh_Type));

    // treat EMesh->edgeN

    for (int type_i = 0; type_i < typeN; type_i++) {

        int *local_topo;

        int mesh_scale = EMesh->scale[type_i];
        int elem_edgeN = LocalEMeshTopo(NMesh.type[type_i], &local_topo);
        EMesh->topo [type_i] = (int*)malloc(elem_edgeN*mesh_scale*sizeof(int));
        EMesh->drict[type_i] = (int*)malloc(elem_edgeN*mesh_scale*sizeof(int));
        
        EMesh->edgeN[type_i] = elem_edgeN;

        for (int elem_i=0; elem_i<mesh_scale; elem_i++) {

            int *node = &NMesh.topo[type_i][elem_i*NMesh.nodeN[type_i]];

            for (int edge_i=0; edge_i<elem_edgeN; edge_i++) {

                int Estart = node[local_topo[edge_i*2+0]-1];
                int Eend   = node[local_topo[edge_i*2+1]-1];

                EMesh->drict[type_i][elem_i*elem_edgeN + edge_i] = Estart > Eend ? -1 : 1;

                if ( Estart > Eend ) { 
                    int temp = Eend;
                    Eend = Estart;
                    Estart = temp;
                }

                int pre_edgeN = ( Estart == 1 ? 0 : Edges.adj_nodn[Estart-2] );
                EMesh->topo[type_i][elem_i*elem_edgeN + edge_i] 
                = pre_edgeN + Binary_Search_(Edges.adj_topo[Estart-1], Edges.adj_nodn[Estart-1]-pre_edgeN, Eend) + 1;
            }
        }

        free(local_topo);
        local_topo = NULL;
    }

    // ...
}
