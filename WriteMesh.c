#include "PDE.h"

void WriteMesh (char *PostMeshFile)
{
    FILE *fp;
    if ((fp = fopen(PostMeshFile,"w"))==NULL)
        return;

    for (int k=1; k<=TypesNum; k++)
    {
        char EType[5];
        char TypeName[20];
        switch(ElemType[k-1])
        {
            case P1 : strcpy(EType,"P1" ); strcpy(TypeName,"Point"); break;
            case L2 : strcpy(EType,"L2" ); strcpy(TypeName,"Line" ); break;
            case L3 : strcpy(EType,"L3" ); strcpy(TypeName,"Line" ); break;
            case T3 : strcpy(EType,"T3" ); strcpy(TypeName,"Triangle"); break;
            case T6 : strcpy(EType,"T6" ); strcpy(TypeName,"Triangle"); break;
            case Q4 : strcpy(EType,"Q4" ); strcpy(TypeName,"Quadrilateral"); break;
            case Q8 : strcpy(EType,"Q8" ); strcpy(TypeName,"Quadrilateral"); break;
            case Q9 : strcpy(EType,"Q9" ); strcpy(TypeName,"Quadrilateral"); break;
            case W4 : strcpy(EType,"W4" ); strcpy(TypeName,"Tetrahedra"); break;
            case W10: strcpy(EType,"W10"); strcpy(TypeName,"Tetrahedra"); break;
            case C8 : strcpy(EType,"C8" ); strcpy(TypeName,"Hexahedra"); break;
            case C20: strcpy(EType,"C20"); strcpy(TypeName,"Hexahedra"); break;
            case C27: strcpy(EType,"C27"); strcpy(TypeName,"Hexahedra"); break;
            case H6 : strcpy(EType,"H6" ); strcpy(TypeName,"Prism"); break;
            case H15: strcpy(EType,"H15"); strcpy(TypeName,"Prism"); break;
            case H18: strcpy(EType,"H18"); strcpy(TypeName,"Prism"); break;
            default: break;
        }        
        
        fprintf(fp,"%s %s %s %d %s %s %s %d\n","Mesh",EType,"Dimension",Coor0.Dim,"Elemtype",TypeName,"Nnode",NodMsh[k-1].Elem_CompN-1);
        
        if(k==1)
        {
            fprintf(fp,"%s\n","Coordinates");
            for (int j=1; j<=Coor0.Total_Nodes; ++j)
            {
                fprintf(fp,"%d",j);
                for (int i=1; i<=Coor0.Dim; ++i)
                {
                    fprintf(fp," %e",Coor0.Coor[(i-1)*Coor0.Total_Nodes+j-1]);
                }
                fprintf(fp,"\n");
            }
            fprintf(fp,"%s\n","End coordinates");
        }
        
        
        fprintf(fp,"%s\n","Elements");
        for (int i=1; i<=NodMsh[k-1].Mesh_Scale; ++i)
        {
            fprintf(fp,"%d",i);
            for (int j=1; j<=NodMsh[k-1].Elem_CompN; ++j)
            {
                fprintf(fp," %d",NodMsh[k-1].Mesh_Topo[(i-1)*NodMsh[k-1].Elem_CompN+j-1]);
            }
            fprintf(fp,"\n");
        }
        fprintf(fp,"%s\n","End elements");
    }
    
    fclose(fp);
    return;
}