#include "fem.h"

void result_compose( Equat_Set Equa, Dof_Tag ID, double* result, int total_nodes, int dof_num ) {

    for (int i=0; i<total_nodes; i++) {

        for (int j=0; j<dof_num; j++) {

            int id_SN = Equa.node_equa_index[i][j];

            if (id_SN > 0)
                result[i*dof_num + j] = Equa.vector[id_SN - 1];
        }
    }

    printf("result composed done!\n");
}