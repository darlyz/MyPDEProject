#include "PDE.h"
#include <omp.h>

void ComposeMatrix(double**, double**, double**, double**,
		              int*, int*, double*, double*, double*);

//int InitialPDE (Coor0, NodMsh, NodF, Mate0, TList0, EqSet0)
Coor Coor0;
Node_Mesh *NodMsh;
Node_Unknow *NodF;
Materail Mate0;
Time_List TList0;
Equation_Set EqSet0;
{
	int **RowMatrix;
	int *tempRow;
	int tempCount = 0;
	int RowCount = 0;
	
	for (int k=1; k<=TypesNum; k++)
		RowCount += NodMsh[k-1].Mesh_Scale;
	
	RowMatrix = (int**)malloc(sizeof(int*)*RowCount);
	for (int k=1; k<=RowCount; k++)
		RowMatrix[k-1] = (int*)malloc(sizeof(int)*101);
	
	for (int k=1; k<=TypesNum; k++)
	{
		for (int ElemIndex=1; ElemIndex<=NodMsh[k-1].Mesh_Scale; ElemIndex++)
		{
			
	
		}
	}
	
	for (int k=1; k<=TypesNum; k++)
	{

		int Node_Count   = NodMsh[k-1].Elem_CompN-1;
		int Elem_CompN   = NodMsh[k-1].Elem_CompN;
		int ElemMtrxSize = NodF[k-1].DoF*Node_Count;

		double **TDerivative0ro = (double**)malloc(sizeof(double*)*NodMsh[k-1].Mesh_Scale);
		double **TDerivative1st = (double**)malloc(sizeof(double*)*NodMsh[k-1].Mesh_Scale);
		double **TDerivative2nd = (double**)malloc(sizeof(double*)*NodMsh[k-1].Mesh_Scale);
		double **RHSData        = (double**)malloc(sizeof(double*)*NodMsh[k-1].Mesh_Scale);

		#pragma omp parallel for
		for(int i=0; i<NodMsh[k-1].Mesh_Scale; i++)
		{
			double *TDerivative0ro[i] = (double*)malloc(sizeof(double)*ElemMtrxSize*ElemMtrxSize);
			double *TDerivative1st[i] = (double*)malloc(sizeof(double)*ElemMtrxSize*ElemMtrxSize);
			double *TDerivative2nd[i] = (double*)malloc(sizeof(double)*ElemMtrxSize*ElemMtrxSize);
			double *RHSData[i]        = (double*)malloc(sizeof(double)*ElemMtrxSize);
		}

		//int ElemMatix2Globle = (int*)malloc(sizeof(int)*size*size);
		//int ElemRHS2Globle   = (int*)malloc(sizeof(int)*size);

		#pragma omp parallel for
		for (int ElemIndex=1; ElemIndex<=NodMsh[k-1].Mesh_Scale; ElemIndex++)
		{
			double *ElemNodeCoor = (double*)malloc(sizeof(double)*Coor0.Dim*Node_Count);
			int MateIndex = NodMsh[k-1].Mesh_Topo[ElemIndex*Elem_CompN-1];
			int NodeIndex;

			for(int i=1; i<=Node_Count; i++)
			{
				NodeIndex = NodMsh[k-1].Mesh_Topo[(ElemIndex-1)*Elem_CompN+i-1];
				for(int j=1; j<=Coor0.Dim; j++)
					ElemNodeCoor[(j-1)*Node_Count+i-1] = Coor0.Coor[(j-1)*NodMsh[k-1].Mesh_Scale+NodeIndex-1];
			}

			// there need to be calculate couple Data as decouple mode
			double *CoupData;
			
			ElemCalculate(TDerivative0ro[ElemIndex-1], TDerivative1st[ElemIndex-1], 
			              TDerivative2nd[ElemIndex-1], RHSData[ElemIndex-1], 
			              ElemIndex, MateIndex, CoupData, ElemNodeCoor);

			free(ElemNodeCoor);
		}
		
		ComposeMatrix(TDerivative0ro, TDerivative1st, TDerivative2nd, RHSData,
		              IDNodeSN[0], IDs[0], Ubf[0], Matrix, RHS); // free ElemMatrix as Composing
		
		
		#pragma omp parallel for
		for(int i=0; i<NodMsh[k-1].Mesh_Scale; i++)
		{
			if(TDerivative0ro[i]!=NULL) free(TDerivative0ro[i]);
			if(TDerivative1st[i]!=NULL) free(TDerivative1st[i]);
			if(TDerivative2nd[i]!=NULL) free(TDerivative2nd[i]);
			if(RHSData[i]!=NULL) free(RHSData[i]);
		}
		if(TDerivative0ro!=NULL) free(TDerivative0ro);
		if(TDerivative1st!=NULL) free(TDerivative1st);
		if(TDerivative2nd!=NULL) free(TDerivative2nd);
		if(RHSData!=NULL) free(RHSData);
	}
}

void ComposeMatrix(double **TDerivative0ro, double **TDerivative1st, double **TDerivative2nd, double **RHSData,
		           int *IDNodeSN, int *IDs, double *Ubf, double *Matrix, double *RHS);
{
	for(int ElemIndex
}

void ElemCalculate(double *TDerivative0ro, double *TDerivative1st, double *TDerivative2nd, double *RHSData, 
			       int ElemIndex, int MateIndex, double *CoupData, double *ElemNodeCoor);
{
	
}
