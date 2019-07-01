/*
 Copyright: Copyright (c) 2018
 Created: 2018-4-19
 Author: Zhang_Licheng
 Title: complex tensor expression in dummy index summation
 All rights reserved
*/
#include "fem.h"

void jumprow(int n,FILE *fp);

void ReadNodeMesh(
    char* dat_file,
    Coor_Info  *Coor,
    Node_Mesh  *Mesh,
    int FieldNum,
    int *dof,
    int *IDNodeNum,
    int **IDNodeSN,
    int **IDs,
    double **Ubf
){
    char   temp_char;
    int    temp_int;
    double temp_double;
    
    FILE *ReadData;
    if((ReadData=fopen(dat_file,"r"))==NULL)
        return;
    
    jumprow(2,ReadData);

    // read coor
    fscanf(ReadData,"%d",&Coor->total_nodes);
    fscanf(ReadData,"%d",&Coor->dim);
    Coor->coordinate =(double*)calloc(Coor->total_nodes*Coor->dim, sizeof(double));
    for (int i=0; i<Coor->total_nodes; i++)
    {
        fscanf(ReadData,"%d",&temp_int);
        for (int j=0; j<Coor->dim; j++)
            fscanf(ReadData,"%lf",&(Coor->coordinate[i + j*Coor->total_nodes]));
    }

    // read ubf, usually we only care for Field A, but we need to jump the data
    dof       = (int* )calloc(FieldNum,sizeof(int));          // DOF of Fields
    IDNodeNum = (int* )calloc(FieldNum,sizeof(int));          // the total number of Nodes which assigned ID
    IDNodeSN  = (int**)calloc(FieldNum,sizeof(int*));         // the Serial Number of Nodes which assigned ID
    IDs       = (int**)calloc(FieldNum,sizeof(int*));         // the ID of DOFs
    Ubf       = (double**)calloc(FieldNum,sizeof(double*));   // the value of IDs

    for (int k=0; k<FieldNum; k++)
    {
        jumprow(2,ReadData);
        fscanf(ReadData,"%d",&(IDNodeNum[k]));
        fscanf(ReadData,"%d",&(dof[k]));

        if(IDNodeNum[k]!=0)
        {
            IDNodeSN[k] = (int*)calloc(IDNodeNum[k]       ,sizeof(int));
            IDs[k]   = (int   *)calloc(IDNodeNum[k]*dof[k],sizeof(int));
            Ubf[k]   = (double*)calloc(IDNodeNum[k]*dof[k],sizeof(double));
        
            for (int i=0; i<IDNodeNum[k]; i++)
            {
                fscanf(ReadData,"%d",&(IDNodeSN[k][i]));
                for (int j=0; j<dof[k]; j++)
                    fscanf(ReadData,"%d",&(IDs[k][i*dof[k] + j]));
            }
        
            jumprow(3,ReadData);
            
            for (int i=0; i<IDNodeNum[k]; i++)
            {
                fscanf(ReadData,"%d",&temp_int);
                for (int j=0; j<dof[k]; j++)
                    fscanf(ReadData,"%lf",&(Ubf[k][i*dof[k] + j]));
            }
        }
        else
            jumprow(3,ReadData);
        
    }
    
    // read elem, usually we only care for Field A, so we ignore others
    int NodeSum = 0;
    long int pos = ftell(ReadData);
    for (int k=0; k<Mesh->typeN; k++)
    {
        jumprow(2,ReadData);
        fscanf(ReadData,"%d",&(Mesh->mesh_scale[k]));
        fscanf(ReadData,"%d",&(Mesh->elem_nodeN[k]));
        printf("%d\t%d\n",Mesh->mesh_scale[k], Mesh->elem_nodeN[k]);
        jumprow(Mesh->mesh_scale[k]+1,ReadData);
        NodeSum += Mesh->mesh_scale[k]*Mesh->elem_nodeN[k];
    }
    Mesh->mesh_topo = (int*)calloc(NodeSum, sizeof(int));

    fseek(ReadData, pos, SEEK_SET);
    NodeSum = 0;
    for (int k=0; k<Mesh->typeN; k++)
    {
        jumprow(3,ReadData);
        for (int i=0; i<Mesh->mesh_scale[k]; i++)
        {
            fscanf(ReadData,"%d",&temp_int);
            for (int j=0; j<Mesh->elem_nodeN[k]; j++)
                fscanf(ReadData,"%d",&(Mesh->mesh_topo[NodeSum + i*Mesh->elem_nodeN[k] + j]));
        }
        jumprow(1,ReadData);
        NodeSum += Mesh->mesh_scale[k]*Mesh->elem_nodeN[k];
    }
    
    return;
}


void jumprow(int n,FILE *fp) // jump n rows
{
    char voidchar[255];
    int i;
    for (i=1;i<=n;i++)
        fgets(voidchar,255,fp);
    return;
}