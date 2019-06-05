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

        double det = transe_coor(real_coor,refr_coor,E_info.node_coor,jacb_matr,E_info.dim,E_info.node_cont);

        double invt_det = inv(invt_jacb,jacb_matr,E_info.dim);

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

	for (int i=0; i<E_info.node_cont*E_info.node_cont; i++)
	    printf("%e\n",E_matr->matr_0[i]);
    
	for (int i=0; i<E_info.node_cont; i++)
	    printf("%e\n",E_matr->righ_vect[i]);

	printf("\n");

}


void set_gaus(Gaus_Info* G_info)
{
	G_info->gaus_num = 4;
	G_info->refc_num = 2;
	G_info->gaus_coor = (double*)malloc(G_info->gaus_num*G_info->refc_num*sizeof(double));
	G_info->gaus_weig = (double*)malloc(G_info->gaus_num*sizeof(double));
	
	G_info->gaus_coor[(1-1)*(4)+1-1] =  5.773502692e-001;
    G_info->gaus_coor[(2-1)*(4)+1-1] =  5.773502692e-001;
	
    G_info->gaus_weig[1-1]=1.000000000e+000;
	
    G_info->gaus_coor[(1-1)*(4)+2-1] =  5.773502692e-001;
    G_info->gaus_coor[(2-1)*(4)+2-1] = -5.773502692e-001;
	
    G_info->gaus_weig[2-1]=1.000000000e+000;
	
    G_info->gaus_coor[(1-1)*(4)+3-1] = -5.773502692e-001;
    G_info->gaus_coor[(2-1)*(4)+3-1] =  5.773502692e-001;
	
    G_info->gaus_weig[3-1]=1.000000000e+000;
	
    G_info->gaus_coor[(1-1)*(4)+4-1] = -5.773502692e-001;
    G_info->gaus_coor[(2-1)*(4)+4-1] = -5.773502692e-001;
	
    G_info->gaus_weig[4-1]=1.000000000e+000;
}

double lagrange_shapfunc(double *shap_coor, int node_i) //quadrilateral
{
	double u,v,shap_val;
	u = shap_coor[0];
	v = shap_coor[1];
	switch (node_i)
    {
    case 1: shap_val = +(+1.-u)*(+1.-v)/4.;
        break;
    case 2: shap_val = +(+1.+u)*(+1.-v)/4.;
        break;
    case 3: shap_val = +(+1.+u)*(+1.+v)/4.;
        break;
    case 4: shap_val = +(+1.-u)*(+1.+v)/4.;
        break;
	}
	return shap_val;
}
double lagrange_deriva_shapfunc(double *shap_coor, int node_i, int coor_i) // quadrilateral
{
	double u,v,shap_val;
	u = shap_coor[0];
	v = shap_coor[1];
	switch (node_i)
    {
    case 1: 
		switch(coor_i){
			case 1: shap_val = -(+1.-v)/4.; break;
			case 2: shap_val = -(+1.-u)/4.; break;
		}
        break;
    case 2: 
		switch(coor_i){
			case 1: shap_val = +(+1.-v)/4.; break;
			case 2: shap_val = -(+1.+u)/4.; break;
		}
        break;
    case 3: 
		switch(coor_i){
			case 1: shap_val = +(+1.+v)/4.; break;
			case 2: shap_val = +(+1.+u)/4.; break;
		}
        break;
    case 4:
		switch(coor_i){
			case 1: shap_val = -(+1.+v)/4.; break;
			case 2: shap_val = +(+1.-u)/4.; break;
		}
        break;
	}
	return shap_val;
}

double det(double* matr, int dim)
{
	if      (dim == 1) return matr[0];
	else if (dim == 2) return matr[(1-1)*dim + 1-1]*matr[(2-1)*dim + 2-1] - matr[(1-1)*dim + 2-1]*matr[(2-1)*dim + 1-1];
	else
	{
		double  sdet  = 0.0;
		double *sub_matr = (double*)malloc((dim-1)*(dim-1)*sizeof(double)); 
		for (int i=1; i<=dim; i++)
		{
			int sign = (i%2==0)?-1:1;
			for (int r=2,rr=1; r<=dim; r++,rr++)
				for (int c=1,cc=0; c<=dim; c++)
				{
					if (c == i) continue;
					sub_matr[(rr-1)*(dim-1) + (cc++)] = matr[(r-1)*dim + c-1];
				}
			sdet += sign*matr[i-1]*det(sub_matr,dim-1);
		}
		free(sub_matr);
		return sdet;
	}
}

double inv(double* inv_matr, double* matr, int dim)
{
	double sdet = det(matr,dim);
	if (fabs(sdet) < 1e-30) return 1e30;
	double *sub_matr = (double*)malloc((dim-1)*(dim-1)*sizeof(double));
	for (int i=1; i<=dim; i++)
		for (int j=1; j<=dim; j++)
		{
			int sign = ((i+j)%2==0)?1:-1;
			for (int r=1,rr=0; r<=dim; r++)
			{
				if(r == i) continue;
				rr++;
				for (int c=1,cc=0; c<=dim; c++)
				{
					if(c == j) continue;
					cc++;
					sub_matr[(rr-1)*(dim-1) + cc-1] = matr[(r-1)*dim + c-1];
				}
			}
			inv_matr[(j-1)*dim + i-1] = sign*det(sub_matr,dim-1)/sdet;
		}
	free(sub_matr);
	return det(inv_matr,dim);
}

double jacobi(
	double* node_coor,
	double* refr_coor,
	double* jacb_matr,
	int dim,
	int shap_nodn
){
	for (int i=1; i<=dim; ++i)
		for (int j=1; j<=dim; ++j)
			for (int node_i=1; node_i<=shap_nodn; ++node_i)
				jacb_matr[(j-1)*dim+i-1] += node_coor[(j-1)*shap_nodn + node_i-1]*lagrange_deriva_shapfunc(refr_coor,node_i,i);
	return det(jacb_matr,dim);
}

double transe_coor(
	double* real_coor,
	double* refr_coor,
	double* node_coor,
	double* jacb_matr,
	int dim,
	int shap_nodn
){
    for (int i=1; i<=dim; i++)
		for (int node_i=1; node_i<=shap_nodn; node_i++)
			real_coor[i-1] += node_coor[(i-1)*shap_nodn + node_i-1]*lagrange_shapfunc(refr_coor,node_i);
	return jacobi(node_coor,refr_coor,jacb_matr,dim,shap_nodn);
}

void calc_refr_shap( double* refr_shap, double* refr_coor,int node_cont, int dim)
{
	for (int node_i=1; node_i<=node_cont; node_i++)
	{
		refr_shap[(node_i-1)*(dim+1)] = lagrange_shapfunc(refr_coor, node_i);
		for (int dim_i=1; dim_i<=(dim+1); dim_i++)
			refr_shap[(node_i-1)*(dim+1) + dim_i] = lagrange_deriva_shapfunc(refr_coor,node_i,dim_i);
	}
}

void set_refr_shap(double** refr_shap, double* gaus_coor, int gaus_num, int node_cont, int dim)
{
	double refr_coor[]={0.0, 0.0, 0.0};
	for (int gaus_i = 1; gaus_i <= gaus_num; gaus_i ++)
    {
		for (int i = 1; i <= dim; i ++)
            refr_coor[i-1] = gaus_coor[(i-1)*gaus_num + gaus_i-1];
		calc_refr_shap( refr_shap[gaus_i-1], refr_coor, node_cont, dim);
	}
}

void calc_real_shap(int dim, int node_cont, double* refr_shap, double* real_shap, double* invt_jacb)
{
	for (int node_i=1; node_i<=node_cont; node_i++)
	{
		real_shap[(node_i-1)*(dim+1)] = refr_shap[(node_i-1)*(dim+1)];
		for (int dim_i=1; dim_i<=dim; dim_i++)
		{
			double temp = 0.0;
			for (int k=1; k<=dim; ++k)
				temp += invt_jacb[(k-1)*dim + dim_i-1] * refr_shap[(node_i-1)*(dim+1) + k];
			real_shap[(node_i-1)*(dim+1) + dim_i] = temp;
		}
	}
}
