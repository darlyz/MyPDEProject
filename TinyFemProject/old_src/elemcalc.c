/*
 Copyright: Copyright (c) 2019
 Created: 2019-4-19
 Author: Zhang_Licheng
 All rights reserved
*/
#include "fem.h"

double transe_coor(double*,double*,double*,double*,int,int);
double inv(double*,double*,int);
void calc_real_shap(int, int, double*, double*, double*);
void show_elem_stif(int, Elem_Matr);

void elemcalc (
    int elem_n,
    Gaus_Info  G_info,
	Elem_Info  E_info,
    Elem_Matr *E_matr
){
    double x,y,z;
    double rx,ry,rz;
	double tempvar;

    double  u[E_info.node_cont];
    double ux[E_info.node_cont];
	double uy[E_info.node_cont];

    double eps,rho;
    eps = E_info.elem_mate[0];
    rho = E_info.elem_mate[1];

	double real_coor[3];
	double refr_coor[3];
	double jacb_matr[9];
	double invt_jacb[9];
	double real_shap[E_info.node_cont*(E_info.dim+1)];

    for (int gaus_i = 1; gaus_i <= G_info.gaus_num; gaus_i ++)
    {
        for (int i=1;i<=E_info.dim;i++)
		{
            refr_coor[i-1] = G_info.gaus_coor[(i-1)*G_info.gaus_num + gaus_i-1];
			real_coor[i-1] = 0.;
			for (int j=1;j<=E_info.dim;j++)
				jacb_matr[(i-1)*E_info.dim+j-1]=0.;
		}

        double det = transe_coor(real_coor, refr_coor, E_info.node_coor, jacb_matr, E_info.dim, E_info.node_cont);

        double invt_det = inv(invt_jacb, jacb_matr, E_info.dim);

        x  = real_coor[0];
        y  = real_coor[1];

        rx = refr_coor[0];
        ry = refr_coor[1];

		double weight = det*G_info.gaus_weig[gaus_i-1];
		calc_real_shap(E_info.dim, E_info.node_cont, E_info.refr_shap[gaus_i-1], real_shap, invt_jacb);

		for (int i=1; i<=E_info.node_cont; ++i)
            u[i-1]  = +real_shap[(i-1)*(E_info.dim+1)+1-1];

		for (int i=1; i<=E_info.node_cont; ++i)
            ux[i-1] = +real_shap[(i-1)*(E_info.dim+1)+2-1];

        for (int i=1; i<=E_info.node_cont; ++i)
            uy[i-1] = +real_shap[(i-1)*(E_info.dim+1)+3-1];

		// stiffness calculation
		for (int i=1; i<=E_info.node_cont; ++i)
			for (int j=1; j<=E_info.node_cont; ++j)
			{
				tempvar = ux[i-1]*ux[j-1]*eps
						+ uy[i-1]*uy[j-1]*eps;
				E_matr->matr_0[(i-1)*E_info.node_cont + j-1] += tempvar*weight;
			}

		for (int i=1; i<=E_info.node_cont; ++i)
			for (int j=1; j<=E_info.node_cont; ++j)
			{
				tempvar = u[i-1]*u[j-1];
				E_matr->matr_0[(i-1)*E_info.node_cont + j-1] += tempvar*weight;
			}
		
		// loadvector calculation
		for (int i=1; i<=E_info.node_cont; ++i)
		{
			tempvar = u[i-1]*rho;
			E_matr->righ_vect[i-1] += tempvar*weight;
		}
    }

	//show_elem_stif(E_info.node_cont, *E_matr);

}

void show_elem_stif(int node_cont, Elem_Matr E_matr)
{
	for (int i=0; i<node_cont*node_cont; i++)
	    printf("%e\n",E_matr.matr_0[i]);
    
	for (int i=0; i<node_cont; i++)
	    printf("%e\n",E_matr.righ_vect[i]);

	printf("\n");
}