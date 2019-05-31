#include "PDE.h"

int compare_int(const void *a, const void *b) // for qsort(int*, int size, sizeof(int), compare_int)
{
	return (*(int*)a-*(int*)b);//a>b 返回正值
}

int compare_char(const void *a, const void *b) // for qsort(char*, int size, sizeof(char), compare_char)
{
    return(*(char*)a-*(char*)b);
}

int compare_double(const void *a, const void *b) // for qsort(double*, int size, sizeof(double), compare_double)
{
    if(fabs(*(double*)a-*(double *)b)<1*exp(-20))
        return 0;
    else
        return(((*(double*)a-*(double*)b)>0)?1:-1);
}

