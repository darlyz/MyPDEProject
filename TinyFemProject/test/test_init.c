#include "fem.h"

void readmesh();
void print_mesh();
void print_coor();
void initial();

int main (int argc, char* argv[])
{
    char data_file[]="../mesh/exam.gid/exam.dat";
    Coor_Info Coor;
    Node_Mesh Mesh;
    Init_Data Init;
    Equat_Set Equa;
    Dof_Tag   ID;
    int node_dof = 1;

    readmesh(data_file, &Coor, &Mesh, &Init, &ID);

    //print_mesh(Mesh);
    //print_coor(Coor);

    clock_t start, finish;
    double  duration;

    start  = clock();
    initial(Coor, Mesh, Equa, ID, node_dof);
    finish = clock();
    
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "%f seconds\n", duration );
}


void print_mesh(Node_Mesh Mesh)
{
    for (int k=0; k<Mesh.typeN; k++) {
        printf("Mesh Type %d: %d\n", k, Mesh.type[k]);
        for (int i=0; i<Mesh.mesh_scale[k]; i++) {
            for (int j=0; j<Mesh.elem_nodeN[k]; j++) 
                printf("%d ",Mesh.mesh_topo[k][i*Mesh.elem_nodeN[k]+j]);
            printf("\n");
        }
    }
}
void print_coor(Coor_Info Coor)
{
    printf("%d %d\n", Coor.total_nodes, Coor.dim);
    for (int i=0; i<Coor.total_nodes; i++) {
        for (int j=0; j<Coor.dim; j++)
            printf("%lg ",Coor.coordinate[i*Coor.dim + j]);
        printf("\n");
    }
}
