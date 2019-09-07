#include "fem.h"
#include "az_aztec.h"

void matrsolv(Equat_Set Equa) {

    int total_nontriaval = Equa.total_nontriaval;
    int total_equations  = Equa.total_equations;
    int *row_nontriaval  = Equa.row_nontriaval;
    int **column_index   = Equa.column_index;
    double **matrix      = Equa.matrix;

    int    *bindx = (int   *)calloc(total_nontriaval+1, sizeof(int));
    double *val   = (double*)calloc(total_nontriaval+1, sizeof(double));

    // convert matrix to MSR format
    int nontrivial_count = 0, temp_count = 0;
    for (int eq_i=0; eq_i<total_equations; eq_i++) {

        if (eq_i==0)
            bindx[0] = total_equations + 1;
        
        else {
            nontrivial_count += row_nontriaval[eq_i-1];
            bindx[eq_i] = nontrivial_count;
        }

        for (int clm_i=0; clm_i<row_nontriaval[eq_i]; clm_i++) {

            if (column_index[eq_i][clm_i] == eq_i) {

                val[eq_i] = matrix[eq_i][clm_i];
                continue;
            }

            temp_count ++;

            bindx[total_equations + temp_count] = column_index[eq_i][clm_i];
            val  [total_equations + temp_count] = matrix      [eq_i][clm_i];
        }
    }
    bindx[total_equations] = nontrivial_count + row_nontriaval[total_equations-1];

    // 
    int N_update = total_equations;
    int *update  = (int*)malloc(N_update * sizeof(int));
    for (int i = 0; i < N_update; i++)
        update [i] = i;

    double *b=(double *) calloc(N_update,sizeof(double));
    double *x=(double *) calloc(N_update,sizeof(double));

    int  noptions[AZ_OPTIONS_SIZE];
    double params[AZ_PARAMS_SIZE];

    AZ_defaults(noptions, params);
    choose_solver_options(noptions, params);

    params  [AZ_tol      ] = 1e-8;
    //params  [Az_dropout  ] = 0.  ;
    params  [AZ_omega    ] = 1.  ;
    noptions[AZ_max_iter ] = 1000;
    noptions[AZ_solver   ] = AZ_bicgstab ;
    noptions[AZ_precond  ] = AZ_ls ;
    noptions[AZ_poly_ord ] = 7 ;
    noptions[AZ_scaling  ] = AZ_none ;
    noptions[AZ_output   ] = 9 ;


    //----------------------------------------------------------------

    int    nproc_config[AZ_PROC_SIZE];
    double status[AZ_STATUS_SIZE];

    int *ndata_org;
    int *nexternal;
    int *nupdate_index;
    int *nextern_index;

    AZ_set_proc_config(nproc_config, AZ_NOT_MPI);

    int MSRorVBR = AZ_MSR_MATRIX;

    AZ_solve(x,b, noptions, params,NULL,bindx,NULL,NULL,
             NULL,val, ndata_org, status, nproc_config);

}