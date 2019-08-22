/*
 Copyright: Copyright (c) 2019
 Created: 2019-4-19
 Author: Zhang_Licheng
 All rights reserved
*/
#include "fem.h"

void readmesh(char*, Coor_Info*, Node_Mesh*, Init_Data*, Dof_Tag*);
void show_coor (Coor_Info);
void show_ID   (Dof_Tag  );
void show_value(Dof_Tag  );
void show_init (Init_Data);
void show_mesh (Node_Mesh);

int main (int argc, char* argv[])
{
    Coor_Info Coor;
	Node_Mesh Mesh;
	Materail  Mate;
	Equat_Set Equa;
    Init_Data Init;
    Dof_Tag   ID;
	int node_dof = 1;

    if (argc > 1) {
        readmesh(argv[1], &Coor, &Mesh, &Init, &ID);
    }

    if (argc == 3) {

        if (strncmp(argv[2],"-coor" ,5) == 0) {
            show_coor(Coor);
        }
        else if (strncmp(argv[2],"-ID"   ,3) == 0) {
            show_ID(ID);
        }
        else if (strncmp(argv[2],"-init" ,5) == 0) {
            show_init(Init);
        }
        else if (strncmp(argv[2],"-mesh" ,5) == 0) {
            show_mesh(Mesh);
        }
    }

}

void show_coor(Coor_Info Coor)
{
    printf("%d %d\n", Coor.total_nodes, Coor.dim);
    for (int i=0; i<Coor.total_nodes; i++) {
        for (int j=0; j<Coor.dim; j++)
            printf("%lg ",Coor.coordinate[i*Coor.dim + j]);
        printf("\n");
    }
}

void show_ID(Dof_Tag ID)
{
    printf("%d %d\n", ID.tag_nodn, ID.dof_num);
    for (int i=0; i<ID.tag_nodn; i++) {
        printf("%d ",ID.tag_nods[i]);
        for (int j=0; j<ID.dof_num; j++)
            printf("%d ",ID.dof_tag[i*ID.dof_num + j]);
        printf("\n");
    }
    printf("-------------------------\n");
    for (int i=0; i<ID.tag_nodn; i++) {
        printf("%d ",ID.tag_nods[i]);
        for (int j=0; j<ID.dof_num; j++)
            printf("%lg ",ID.dof_val[i*ID.dof_num + j]);
        printf("\n");
    }
}

void show_init(Init_Data Init)
{
    printf("%d %d\n", Init.total_nodes, Init.dof_num);
    for (int k=0; k<=Init.init_oder; k++) {
        printf("Init_%d\n",k);
        for (int i=0; i<Init.total_nodes; i++) {
            for (int j=0; j<Init.dof_num; j++) {
                if (k==0)
                    printf("%lg ",Init.init_0[i*Init.dof_num + j]);
                if (k==1)
                    printf("%lg ",Init.init_1[i*Init.dof_num + j]);
                if (k==2)
                    printf("%lg ",Init.init_2[i*Init.dof_num + j]);
            }
            printf("\n");
        }
    }
}

void show_mesh(Node_Mesh Mesh)
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