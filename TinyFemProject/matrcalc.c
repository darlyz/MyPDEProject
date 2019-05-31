#include "fem.h"

void matrcalc(
    coor_info Coor,
    Node_Mesh Mesh,
    int node_dof,
    int *dof_tag,
    double *bnd_constr
){
    int pre_node = 0;
    for (int type_i=1; type_i<=Mesh.Mesh_TypeN; type_i++)
    {
        int NodeN = Mesh.Elem_NodeN[type_i-1];
        int MeshN = Mesh.Mesh_Scale[type_i-1];
        int elem_dof = node_dof*(NodeN-1);

        matr_type matr_type[4]={dist,lump,lump,lump};

        elem_matr E_matr;
        set_matr(&E_matr, elem_dof, matr_type);

        for (int elem_i=1; elem_i<MeshN; elem_i++)
        {
            int node[NodeN-1];
            memcpy(node, &Mesh.Mesh_Topo[pre_node+(elem_i-1)*NodeN], (NodeN-1)*sizeof(int));

            double E_coor[(NodeN-1)*Coor.dim];
            for (int node_i=1; node_i<NodeN; node_i++)
            {
                int node_num = node[node_i-1];
                memcpy(E_coor,&Coor.coordinate[(node_num-1)*Coor.dim],Coor.dim*sizeof(double));
                elemclac(E_coor, coup_valu, elem_mate, elem_i, refr_shap, gauss, einfo, &E_matr)
            }
        }
        pre_node += NodeN * MeshN;
    }
}

void set_matr(elem_matr &E_matr, int elem_dof, matr_type *M_type)
{
    E_matr->matr_0    = (double*)malloc(pow(elem_dof,M_type[0])*sizeof(double));
    E_matr->matr_1    = (double*)malloc(pow(elem_dof,M_type[1])*sizeof(double));
    E_matr->matr_2    = (double*)malloc(pow(elem_dof,M_type[2])*sizeof(double));
    E_matr->matr_3    = (double*)malloc(pow(elem_dof,M_type[3])*sizeof(double));
    E_matr->left_matr = (double*)malloc(pow(elem_dof,M_type[0])*sizeof(double));
    E_matr->righ_vec  = (double*)malloc(elem_dof*sizeof(double));
}
    
void clear_matr(elem_matr &E_matr)
{
    free(E_matr->matr_0);
    free(E_matr->matr_1);
    free(E_matr->matr_2);
    free(E_matr->matr_3);
}