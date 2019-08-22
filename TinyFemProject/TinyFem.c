#include "fem.h"

int main(int argc, char* argv[])
{

    char data_file[]="./mesh/exam.gid/exam.dat";
    Coor_Info Coor;
    Node_Mesh Mesh;
    Init_Data Init;
    Equat_Set Equa;
    Dof_Tag   ID;
    int node_dof = 1;

    readmesh(data_file, &Coor, &Mesh, &Init, &ID);
    initial(Coor, Mesh, Equa, ID, node_dof);

}