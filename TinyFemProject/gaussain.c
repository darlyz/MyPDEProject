#include "fem.h"

void set_gaus(Gaus_Info* G_info)
{
	G_info->gaus_num = 4;
	G_info->dim = 2;
	G_info->gaus_coor = (double*)malloc(G_info->gaus_num*G_info->dim*sizeof(double));
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