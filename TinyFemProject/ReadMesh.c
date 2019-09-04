/*
 Copyright: Copyright (c) 2018
 Created: 2018-4-19
 Author: Zhang_Licheng
 Title: complex tensor expression in dummy index summation
 All rights reserved
*/
#include "fem.h"

void readmesh(
    Coor_Info *Coor,
    Node_Mesh *Mesh,
    Dof_Tag   *ID,
    Init_Data *Init,
    char      *dat_file
){
    char   temp_char;
    char   temp_str[255];
    int    temp_int;
    double temp_double;

    FILE *ReadData;
    if((ReadData=fopen(dat_file,"r"))==NULL)
        return;

    typedef enum{coor, id, value, init, mesh} paragh;
    paragh dataparagh;
    int line_count;

    Init->init_oder = -1;
    Mesh->typeN     =  0;
    Mesh->type = (Mesh_Type* )malloc(sizeof(Mesh_Type)*16);
    Mesh->elem_nodeN = (int* )malloc(sizeof(int )*16);
    Mesh->mesh_scale = (int* )malloc(sizeof(int )*16);
    Mesh->mesh_topo  = (int**)malloc(sizeof(int*)*16);

    while(fgets(temp_str, 255, ReadData) != NULL) {

        line_count ++;

        if (strncmp(temp_str,"-coordinate" ,3) == 0) {
            dataparagh = coor;
            line_count = 0;
            continue;
        }

        else if (strncmp(temp_str,"-ID-"   ,3) == 0) {
            dataparagh = id;
            line_count = 0;
            continue;
        }

        else if (strncmp(temp_str,"-value-",3) == 0) {
            dataparagh = value;
            line_count = 0;
            continue;
        }

        else if (strncmp(temp_str,"-init-" ,3) == 0) {
            dataparagh = init;
            line_count = 0;
            Init->init_oder ++ ;
            continue;
        }

        else if (strncmp(temp_str,"-mesh-" ,3) == 0) {
            dataparagh = mesh;
            line_count = 0;
            Mesh->typeN ++ ;
            continue;
        }

        if (dataparagh == coor) {

            if (line_count == 1) {

                sscanf(temp_str, "%d %d", &Coor->total_nodes, &Coor->dim);

                Coor->coordinate = (double*)malloc(Coor->total_nodes
                                                  *Coor->dim
                                                  *sizeof(double));
            }
            else if (line_count < Coor->total_nodes + 2) {

                int  word_count = 0;
                int  word_begin = 0;
                bool word_find  = true;
                for (int i=0; temp_str[i] != '\0'; i++) {

                    if ( isspace(temp_str[i]) && word_find ) {
                        word_find = false;
                        temp_str[i] = '\0';
                        sscanf(temp_str + word_begin, "%lg", &Coor->coordinate[(line_count-2)*Coor->dim + word_count]);
                        word_count ++;
                    }

                    else if ( !isspace(temp_str[i]) && !word_find ) {
                        word_find  = true;
                        word_begin = i;
                    }
                }
            }
        }

        else if (dataparagh == id) {

            if (line_count == 1) {

                sscanf(temp_str, "%d %d", &ID->tag_nodn, &ID->dof_num);

                ID->tag_nods = (int*   )malloc(ID->tag_nodn*sizeof(int));
                ID->dof_tag  = (int*   )malloc(ID->tag_nodn*ID->dof_num*sizeof(int));
                ID->dof_val  = (double*)malloc(ID->tag_nodn*ID->dof_num*sizeof(double));
            }
            else if (line_count < ID->tag_nodn + 2) {

                int  word_count = 0;
                int  word_begin = 0;
                bool word_find  = true;

                for (int i=0; temp_str[i] != '\0'; i++) {

                    if ( isspace(temp_str[i]) && word_find ) {

                        word_find = false;
                        temp_str[i] = '\0';

                        if (word_count == 0)
                            sscanf(temp_str + word_begin, "%d", &ID->tag_nods[line_count-2]);
                        
                        else
                            sscanf(temp_str + word_begin, "%d", &ID->dof_tag[(line_count-2)*ID->dof_num + word_count-1]);
                        
                        word_count ++;
                    }

                    else if ( !isspace(temp_str[i]) && !word_find ) {
                        word_find  = true;
                        word_begin = i;
                    }
                }
            }
        }

        else if (dataparagh == value) {

            if (line_count == 1) {
                continue;
            }

            else if (line_count < ID->tag_nodn + 2) {

                int  word_count = 0;
                int  word_begin = 0;
                bool word_find  = true;

                for (int i=0; temp_str[i] != '\0'; i++) {

                    if ( isspace(temp_str[i]) && word_find ) {

                        word_find = false;
                        temp_str[i] = '\0';

                        if (word_count != 0)
                            sscanf(temp_str + word_begin, "%lg", &ID->dof_val[(line_count-2)*ID->dof_num + word_count-1]);
                        
                        word_count ++;
                    }

                    else if ( !isspace(temp_str[i]) && !word_find ) {
                        word_find  = true;
                        word_begin = i;
                    }
                }
            }
        }

        else if (dataparagh == init) {

            if (line_count == 1) {

                if (Init->init_oder == 0) {
                    sscanf(temp_str, "%d %d", &Init->total_nodes, &Init->dof_num);
                    Init->init_0 = (double*)malloc(Init->total_nodes*Init->dof_num*sizeof(double));
                }

                else if (Init->init_oder == 1)
                    Init->init_1 = (double*)malloc(Init->total_nodes*Init->dof_num*sizeof(double));

                else if (Init->init_oder == 2)
                    Init->init_2 = (double*)malloc(Init->total_nodes*Init->dof_num*sizeof(double));

            }

            else if (line_count < Init->total_nodes + 2) {

                int  word_count = 0;
                int  word_begin = 0;
                bool word_find  = true;

                for (int i=0; temp_str[i] != '\0'; i++) {

                    if ( isspace(temp_str[i]) && word_find ) {

                        word_find = false;
                        temp_str[i] = '\0';

                        if (word_count != 0) {

                            if (Init->init_oder == 0)
                                sscanf(temp_str + word_begin, "%lg", &Init->init_0[(line_count-2)*Init->dof_num + word_count-1]);

                            else if (Init->init_oder == 1)
                                sscanf(temp_str + word_begin, "%lg", &Init->init_1[(line_count-2)*Init->dof_num + word_count-1]);

                            else if (Init->init_oder == 2)
                                sscanf(temp_str + word_begin, "%lg", &Init->init_2[(line_count-2)*Init->dof_num + word_count-1]);

                        }

                        word_count ++;
                    }

                    else if ( !isspace(temp_str[i]) && !word_find ) {
                        word_find  = true;
                        word_begin = i;
                    }
                }
            }
        }

        else if (dataparagh == mesh) {

            if (line_count == 1) {

                /*
                if (strncmp(temp_str, "l2", 2) == 0)
                    Mesh->type[Mesh->typeN-1] = L2;

                else if (strncmp(temp_str, "l3", 2) == 0)
                    Mesh->type[Mesh->typeN-1] = L3;

                else if (strncmp(temp_str, "t3", 2) == 0)
                    Mesh->type[Mesh->typeN-1] = T3;

                else if (strncmp(temp_str, "t6", 2) == 0)
                    Mesh->type[Mesh->typeN-1] = T6;

                else if (strncmp(temp_str, "q4", 2) == 0)
                    Mesh->type[Mesh->typeN-1] = Q4;

                else if (strncmp(temp_str, "q8", 2) == 0)
                    Mesh->type[Mesh->typeN-1] = Q8;

                else if (strncmp(temp_str, "q9", 2) == 0)
                    Mesh->type[Mesh->typeN-1] = Q9;

                else if (strncmp(temp_str, "w4", 2) == 0)
                    Mesh->type[Mesh->typeN-1] = W4;

                else if (strncmp(temp_str, "w10", 3) == 0)
                    Mesh->type[Mesh->typeN-1] = W10;

                else if (strncmp(temp_str, "c8", 2) == 0)
                    Mesh->type[Mesh->typeN-1] = C8;

                else if (strncmp(temp_str, "c20", 3) == 0)
                    Mesh->type[Mesh->typeN-1] = C20;

                else if (strncmp(temp_str, "c27", 3) == 0)
                    Mesh->type[Mesh->typeN-1] = C27;

                else if (strncmp(temp_str, "h6", 3) == 0)
                    Mesh->type[Mesh->typeN-1] = H6;

                else if (strncmp(temp_str, "h15", 3) == 0)
                    Mesh->type[Mesh->typeN-1] = H15;

                else if (strncmp(temp_str, "h18", 3) == 0)
                    Mesh->type[Mesh->typeN-1] = H18;
                */
                
                char shap_type = temp_str[0];
                int  shap_nodn;
                sscanf(temp_str + 1, "%d", &shap_nodn);

                switch (shap_type) {
                    case 'l':
                        switch (shap_nodn) {
                            case 2:
                                Mesh->type[Mesh->typeN-1] = L2;
                                break;
                            case 3:
                                Mesh->type[Mesh->typeN-1] = L3;
                                break;
                            default:
                                printf("Unknown shap type: %s!\n",temp_str);
                                return;
                        }
                        break;
                    case 't':
                        switch (shap_nodn) {
                            case 3:
                                Mesh->type[Mesh->typeN-1] = T3;
                                break;
                            case 6:
                                Mesh->type[Mesh->typeN-1] = T6;
                                break;
                            default:
                                printf("Unknown shap type: %s!\n",temp_str);
                                return;
                        }
                        break;
                    case 'q':
                        switch (shap_nodn) {
                            case 4:
                                Mesh->type[Mesh->typeN-1] = Q4;
                                break;
                            case 8:
                                Mesh->type[Mesh->typeN-1] = Q8;
                                break;
                            case 9:
                                Mesh->type[Mesh->typeN-1] = Q9;
                                break;
                            default:
                                printf("Unknown shap type: %s!\n",temp_str);
                                return;
                        }
                        break;
                    case 'w':
                        switch (shap_nodn) {
                            case 4:
                                Mesh->type[Mesh->typeN-1] = W4;
                                break;
                            case 10:
                                Mesh->type[Mesh->typeN-1] = W10;
                                break;
                            default:
                                printf("Unknown shap type: %s!\n",temp_str);
                                return;
                        }
                        break;
                    case 'c':
                        switch (shap_nodn) {
                            case 8:
                                Mesh->type[Mesh->typeN-1] = C8;
                                break;
                            case 20:
                                Mesh->type[Mesh->typeN-1] = C20;
                                break;
                            case 27:
                                Mesh->type[Mesh->typeN-1] = C27;
                                break;
                            default:
                                printf("Unknown shap type: %s!\n",temp_str);
                                return;
                        }
                        break;
                    case 'h':
                        switch (shap_nodn) {
                            case 6:
                                Mesh->type[Mesh->typeN-1] = H6;
                                break;
                            case 15:
                                Mesh->type[Mesh->typeN-1] = H15;
                                break;
                            case 18:
                                Mesh->type[Mesh->typeN-1] = H18;
                                break;
                            default:
                                printf("Unknown shap type: %s!\n",temp_str);
                                return;
                        }
                        break;
                    case 'p':
                        switch (shap_nodn) {
                            case 1:
                                Mesh->type[Mesh->typeN-1] = P1;
                                break;
                            case 5:
                                Mesh->type[Mesh->typeN-1] = P5;
                                break;
                            case 13:
                                Mesh->type[Mesh->typeN-1] = P13;
                                break;
                            case 14:
                                Mesh->type[Mesh->typeN-1] = P14;
                                break;
                            default:
                                printf("Unknown shap type: %s!\n",temp_str);
                                return;
                        }
                        break;
                    default:
                        printf("Unknown shap type: %s!\n",temp_str);
                        return;
                }


            }

            else if (line_count == 2) {

                sscanf(temp_str, "%d %d", &Mesh->mesh_scale[Mesh->typeN-1],
                                          &Mesh->elem_nodeN[Mesh->typeN-1]);

                Mesh->mesh_topo[Mesh->typeN-1]  = (int*)malloc(sizeof(int)*
                                           Mesh->mesh_scale[Mesh->typeN-1]*
                                           Mesh->elem_nodeN[Mesh->typeN-1]);
            }

            else if (line_count < Mesh->mesh_scale[Mesh->typeN-1] + 3) {

                int  word_count = 0;
                int  word_begin = 0;
                bool word_find  = true;

                for (int i=0; temp_str[i] != '\0'; i++) {

                    if ( isspace(temp_str[i]) && word_find ) {

                        word_find = false;
                        temp_str[i] = '\0';

                        sscanf(temp_str + word_begin, "%d", &Mesh->mesh_topo[Mesh->typeN-1][(line_count-3)*Mesh->elem_nodeN[Mesh->typeN-1] + word_count]);

                        word_count ++;
                    }

                    else if ( !isspace(temp_str[i]) && !word_find ) {
                        word_find  = true;
                        word_begin = i;
                    }
                }
            }
        }
    }
    printf("Read mesh done!\n");
}