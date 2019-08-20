#include "fem.h"
double lagrange_shapfunc(double*, int);
double lagrange_deriva_shapfunc(double*, int, int);

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
