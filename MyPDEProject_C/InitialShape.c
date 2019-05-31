/*
 Copyright: Copyright (c) 2018
 Created: 2018-4-19
 Author: Zhang_Licheng
 Title: complex tensor expression in dummy index summation
 All rights reserved
*/
#include "PDE.h"

void SetTetraShap4 (int NumNode, double *val, double *GausCoor)
{
	double u = GausCoor[0];
	double v = GausCoor[1];
	double w = GausCoor[2];
	switch (NumNode)
	{
		case 1: *val =    +u     ;	break;
		case 2: *val =      +v   ;	break;
		case 3: *val =        +w ;	break;
		case 4: *val = +1.-u-v-w ;	break;
		default: break;
	}
	return;
}
void SetTetraShap10(int NumNode, double *val, double *GausCoor)
{
	double u = GausCoor[0];
	double v = GausCoor[1];
	double w = GausCoor[2];
	switch (NumNode)
	{
		case  1: *val = +(-1.+2.*u)          *     u         ; break;
		case  2: *val =         +u           *       v   *4. ; break;
		case  3: *val = +(-1.     +2.*v)     *       v       ; break;
		case  4: *val =              +v      *         w *4. ; break;
		case  5: *val = +(-1.          +2.*w)*         w     ; break;
		case  6: *val =                   +w *     u     *4. ; break;
		case  7: *val =         +u           *(+1.-u-v-w)*4. ; break;
		case  8: *val =              +v      *(+1.-u-v-w)*4. ; break;
		case  9: *val =                   +w *(+1.-u-v-w)*4. ; break;
		case 10: *val = +(+1.-2.*u-2.*v-2.*w)*(+1.-u-v-w)    ; break;
		default: break;
	}
	return;
}
void SetPrismShap6 (int NumNode, double *val, double *GausCoor)
{
	double u = GausCoor[0];
	double v = GausCoor[1];
	double w = GausCoor[2];
	switch (NumNode)
	{
		case 1: *val =      +u   *(+1.-w)/2. ; break;
		case 2: *val =        +v *(+1.-w)/2. ; break;
		case 3: *val = +(+1.-u-v)*(+1.-w)/2. ; break;
		case 4: *val =      +u   *(+1.+w)/2. ; break;
		case 5: *val =        +v *(+1.+w)/2. ; break;
		case 6: *val = +(+1.-u-v)*(+1.+w)/2. ; break;
		default: break;
	}
	return;
}
void SetPrismShap18(int NumNode, double *val, double *GausCoor)
{
	double u = GausCoor[0];
	double v = GausCoor[1];
	double w = GausCoor[2];
	switch (NumNode)
	{
		case  1: *val = +(-1.+2.*u     )*     u   *   w*(-1.+w)/2. ; break;
		case  2: *val =         +u      *       v *4.*w*(-1.+w)/2. ; break;
		case  3: *val = +(-1.     +2.*v)*       v *   w*(-1.+w)/2. ; break;
		case  4: *val =              +v *(+1.-u-v)*4.*w*(-1.+w)/2. ; break;
		case  5: *val = +(+1.-2.*u-2.*v)*(+1.-u-v)   *w*(-1.+w)/2. ; break;
		case  6: *val =         +u      *(+1.-u-v)*4.*w*(-1.+w)/2. ; break;
		case  7: *val = +(-1.+2.*u     )*     u        *(+1.-w*w)  ; break;
		case  8: *val =         +u      *       v *4.  *(+1.-w*w)  ; break;
		case  9: *val = +(-1.     +2.*v)*       v      *(+1.-w*w)  ; break;
		case 10: *val =              +v *(+1.-u-v)*4.  *(+1.-w*w)  ; break;
		case 11: *val = +(+1.-2.*u-2.*v)*(+1.-u-v)     *(+1.-w*w)  ; break;
		case 12: *val =         +u      *(+1.-u-v)*4.  *(+1.-w*w)  ; break;
		case 13: *val = +(-1.+2.*u     )*     u      *w*(+1.+w)/2. ; break;
		case 14: *val =         +u      *       v *4.*w*(+1.+w)/2. ; break;
		case 15: *val = +(-1.     +2.*v)*       v    *w*(+1.+w)/2. ; break;
		case 16: *val =              +v *(+1.-u-v)*4.*w*(+1.+w)/2. ; break;
		case 17: *val = +(+1.-2.*u-2.*v)*(+1.-u-v)   *w*(+1.+w)/2. ; break;
		case 18: *val =         +u      *(+1.-u-v)*4.*w*(+1.+w)/2. ; break;
		default: break;
	}
	return;
}