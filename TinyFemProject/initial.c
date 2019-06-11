#include "fem.h"

void initial(
    Coor_Info Coor,
    Node_Mesh Mesh,
    Equat_Set Equa,
    int node_dof,
    int *dof_tag,
    
){
    int *eq_id = (int*)malloc(Coor.total_nodes * node_dof * sizeof(int));
    for (int i=1; i<=Coor.total_nodes * node_dof; i++)
        eq_id[i] = i;

}