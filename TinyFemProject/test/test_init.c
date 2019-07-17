#include "fem.h"

void ReadNodeMesh();
void print_mesh();
void print_coor();
void initial();

int main (int argc, char* argv[])
{
    Coor_Info Coor;
    Node_Mesh Mesh;
    Equat_Set Equa;
    Dirichlet D1st;
    int *dof;
    int *IDNodeNum;
    int **IDNodeSN;
    int **IDs;
    double **Ubf;

    Mesh.typeN = 2;
    Mesh.mesh_scale = (int*)malloc(Mesh.typeN*sizeof(int));
    Mesh.elem_nodeN = (int*)malloc(Mesh.typeN*sizeof(int));
    ReadNodeMesh(argv[1], &Coor, &Mesh, 2, dof, IDNodeNum, IDNodeSN, IDs, Ubf);
    
    //print_mesh(Mesh);
    //print_coor(Coor);
    
    initial(Coor, Mesh, &Equa, &D1st, dof[0]);

}


void print_mesh(Node_Mesh Mesh)
{
    int NodeSum = 0;
    for (int k=0; k<Mesh.typeN; k++)
    {
        printf("type %d:\n scale:%d, %d per elem\n",k+1, Mesh.mesh_scale[k], Mesh.elem_nodeN[k]);
        for (int i=0; i<Mesh.mesh_scale[k]; i++)
        {
            printf("\t%d",i+1);
            for (int j=0; j<Mesh.elem_nodeN[k]; j++)
                printf("\t%d",Mesh.mesh_topo[NodeSum + i*Mesh.elem_nodeN[k] + j]);
            printf("\n");
        }
        NodeSum += Mesh.mesh_scale[k]*Mesh.elem_nodeN[k];
    }
}
void print_coor(Coor_Info Coor)
{
    for (int i=0; i<Coor.total_nodes; i++)
    {
        printf("\t%d",i+1);
        for (int j=0; j<Coor.dim; j++)
            printf("\t%lf",Coor.coordinate[i + j*Coor.total_nodes]);
        printf("\n");
    }
}
