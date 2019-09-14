#include "fem.h"

void readmate(Materail *Mate, char* mat_file) {
    
    FILE *ReadMate;
    if((ReadMate = fopen(mat_file,"r"))==NULL) {
        printf("Read materials failed!\n");
        return;
    }

    fscanf(ReadMate, "%d", &(Mate->field_num));
    Mate->mate_cont = (int*    )malloc(Mate->field_num * sizeof(int));
    Mate->mate_varN = (int*    )malloc(Mate->field_num * sizeof(int));
    Mate->mate      = (double**)malloc(Mate->field_num * sizeof(double*));
    
    for (int field_i=0; field_i<Mate->field_num; field_i++) {

        fscanf(ReadMate, "%d", &(Mate->mate_cont[field_i]));
        fscanf(ReadMate, "%d", &(Mate->mate_varN[field_i]));

        Mate->mate[field_i] = (double*)malloc(Mate->mate_cont[field_i] * Mate->mate_varN[field_i] * sizeof(double));

        for (int mat_i=0; mat_i< Mate->mate_cont[field_i]; mat_i++)
            for (int var_i; var_i<Mate->mate_varN[field_i]; var_i++)
                fscanf(ReadMate, "%le", &(Mate->mate[field_i][mat_i * (Mate->mate_varN[field_i]) + var_i]));
    }

    printf("Read materials done!\n");
}