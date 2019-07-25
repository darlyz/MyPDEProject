#include "fem.h"
#include "metis.h"
#define list_size 1000000

void nontrival_graph();

void initial(
    Coor_Info Coor,
    Node_Mesh Mesh,
    Equat_Set Equa,
    Dirichlet D1st,
    int node_dof
    
){
    int initial_tag; memset(&initial_tag, 127, sizeof(int));
    int dirichlet_tag = -1;
    int empty_dof_tag =  0;
    int empty_elm_tag =  0;

    int *eq_id;
    eq_id = (int*)malloc(Coor.total_nodes * node_dof * sizeof(int));
    memset(eq_id,  127,  Coor.total_nodes * node_dof * sizeof(int));

    for (int i=0; i<D1st.dirichlet_nodn; i++) {
        int node_i = D1st.dirichlet_nods[i] - 1;
        memcpy(&eq_id[node_i*node_dof], &D1st.dof_tag[i*node_dof], node_dof*sizeof(int));
    }

    int pre_nods = 0;
    for (int type_i=1; type_i<=Mesh.typeN; type_i++) {
        int elem_nodeN = Mesh.elem_nodeN[type_i-1];
        int mesh_scale = Mesh.mesh_scale[type_i-1];

        for (int elem_i=1; elem_i<=mesh_scale; elem_i++) {
            int elem_node[elem_nodeN];
            memcpy(elem_node, &Mesh.mesh_topo[pre_nods + (elem_i-1)*elem_nodeN], elem_nodeN*sizeof(int));
            
            if (elem_node[elem_nodeN-1] <= 0)
                for (int node_i=1; node_i<elem_nodeN; node_i++) 
                    for (int dof_i=1; dof_i<=node_dof; dof_i++)
                        eq_id[(elem_node[node_i-1]-1)*node_dof + dof_i-1] = elem_node[elem_nodeN-1];
        }
        pre_nods += elem_nodeN*mesh_scale;
    }
    
    int eq_serial_num = 0;
    for (int i=0; i<Coor.total_nodes*node_dof; i++) {
        if (eq_id[i] == initial_tag)
            eq_id[i] = ++eq_serial_num;
    }

    pre_nods = 0;
    for (int type_i=1; type_i<=Mesh.typeN; type_i++) {
        int elem_nodeN = Mesh.elem_nodeN[type_i-1];
        int mesh_scale = Mesh.mesh_scale[type_i-1];

        for (int elem_i=1; elem_i<=mesh_scale; elem_i++) {
            int elem_node[elem_nodeN];
            memcpy(elem_node, &Mesh.mesh_topo[pre_nods + (elem_i-1)*elem_nodeN], elem_nodeN*sizeof(int));

            if (elem_node[elem_nodeN-1] <= 0) continue;

            int eq_count = 0;
            int id_set[500];
            for (int node_i=1; node_i<elem_nodeN; node_i++) 
                for (int dof_i=1; dof_i<=node_dof; dof_i++) {
                    if (eq_id[(elem_node[node_i-1]-1)*node_dof + dof_i-1] > 0)
                        id_set[eq_count++] = eq_id[(elem_node[node_i-1]-1)*node_dof + dof_i-1];
                }

            //if (eq_count > 0)
               // nontrival_graph(); // need to design
        }
        pre_nods += elem_nodeN*mesh_scale;

    }


    printf("initial done!\n");
}