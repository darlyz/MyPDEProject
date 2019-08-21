/*
 Copyright: Copyright (c) 2018
 Created: 2018-4-19
 Author: Zhang_Licheng
 Title: complex tensor expression in dummy index summation
 All rights reserved
*/
#include "fem.h"

void readmesh(
    char* dat_file,
    Coor_Info *Coor,
    Node_Mesh *Mesh,
    Init_Data *Init,
    Dof_Tag   *ID
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
    }
}