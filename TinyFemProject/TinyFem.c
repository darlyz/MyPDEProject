#include "fem.h"

void readmesh();
void initial ();
void matrcalc();

int main(int argc, char* argv[])
{

    char data_file[]="./mesh/exam2.gid/exam2.dat";
    Coor_Info Coor;
    Node_Mesh Mesh;
    Init_Data Init;
    Equat_Set Equa;
    Materail  Mate;
    Dof_Tag   ID;
    int node_dof = 1;
    double* result;

    readmesh(&Coor, &Mesh, &ID, &Init, data_file);
    initial ( Coor,  Mesh,  ID, &Equa, &result, node_dof);
    matrcalc( Coor,  Mesh,  ID, &Equa, Mate, result, node_dof);
}