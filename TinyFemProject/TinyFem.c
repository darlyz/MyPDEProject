#include "fem.h"

void readmesh(Coor_Info*, Node_Mesh*, Dof_Tag*, Init_Data*, char*);
void initial (Coor_Info,  Node_Mesh,  Dof_Tag,  Equat_Set*, int);
void matrcalc(Coor_Info,  Node_Mesh,  Dof_Tag,  Equat_Set*, Materail, int);

int main(int argc, char* argv[])
{

    char data_file[]="./mesh/exam.gid/exam.dat";
    Coor_Info Coor;
    Node_Mesh Mesh;
    Init_Data Init;
    Equat_Set Equa;
    Materail  Mate;
    Dof_Tag   ID;
    int node_dof = 1;

    readmesh(&Coor, &Mesh, &ID, &Init, data_file);
    initial ( Coor,  Mesh,  ID, &Equa, node_dof);
    matrcalc( Coor,  Mesh,  ID, &Equa, Mate, node_dof);
}