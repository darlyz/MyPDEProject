#include "fem.h"

int main(int argc, char* argv[])
{

    char data_file[]="./mesh/exam.gid/exam.dat"
    Coor_Info Coor;
    Node_Mesh Mesh;
    Init_Data Init;
    Dof_Tag   ID;

    readmesh(data_file, &Coor, &Mesh, &Init, &ID);

}