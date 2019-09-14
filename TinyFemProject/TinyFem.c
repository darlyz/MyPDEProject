#include "fem.h"

void readmesh();
void readmate();
void initial ();
void matrcalc();
void matrsolv();
void show_mesh();

int main(int argc, char* argv[])
{

    char data_file[]="./mesh/exam3.gid/exam3.dat";
    char mate_file[]="./mate/exam3.mat";
    Coor_Info Coor;
    Node_Mesh Mesh;
    Init_Data Init;
    Equat_Set Equa;
    Materail  Mate;
    Dof_Tag   ID;
    int node_dof = 1;
    double* result;

    readmesh(&Coor, &Mesh, &ID, &Init, data_file);
    readmate(&Mate, mate_file);
    //show_mesh(Mesh);
    //Mesh.typeN = 1;
    initial ( Coor,  Mesh,  ID, &Equa, &result, node_dof);
    matrcalc( Coor,  Mesh,  ID, &Equa, Mate, result, node_dof);
    matrsolv( Equa );
}