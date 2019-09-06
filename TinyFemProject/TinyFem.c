#include "fem.h"

void readmesh();
void initial ();
void matrcalc();
void show_mesh();

int main(int argc, char* argv[])
{

    char data_file[]="./mesh/exam3.gid/exam3.dat";
    Coor_Info Coor;
    Node_Mesh Mesh;
    Init_Data Init;
    Equat_Set Equa;
    Materail  Mate;
    Dof_Tag   ID;
    int node_dof = 1;
    double* result;

    readmesh(&Coor, &Mesh, &ID, &Init, data_file);
    //show_mesh(Mesh);
    //Mesh.typeN = 1;
    initial ( Coor,  Mesh,  ID, &Equa, &result, node_dof);
    matrcalc( Coor,  Mesh,  ID, &Equa, Mate, result, node_dof);
    matrsolv( Equa );
}