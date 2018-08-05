#include"PDE.h"

char* ReadJsonFile  (char*, char*);
int   SetElemType   (char*, int);
void  CreatMeshSpace(char*, int);

void ReadJson(char *JsonFile)
{
	char *JsonData;
	JsonData = ReadJsonFile(JsonFile,JsonData);
	
	cJSON *json = cJSON_Parse(JsonData);
	
	cJSON *RegionMesh = cJSON_GetObjectItem(json,"RegionMesh");
	cJSON *BndaryMesh = cJSON_GetObjectItem(json,"BndaryMesh");
	char  *RegMesh    = cJSON_Print(RegionMesh);
	char  *BndMesh    = cJSON_Print(BndaryMesh);
	
	int NType = SetElemType(RegMesh, 0);
	    NType = SetElemType(BndMesh, NType);
	
	TypesNum = NType;
	
	cJSON *ElemtType  = cJSON_GetObjectItem(json,"ElemtType");
	char  *ElemtT     = cJSON_Print(ElemtType);
	
	CreatMeshSpace(ElemtT,NType);
	
	
	free(JsonData);
	cJSON_Delete(json);
	return;
}


char* ReadJsonFile(char *FileName, char *JsonData)
{
	FILE *fp;
	if((fp=fopen(FileName,"r"))==NULL)
		return NULL;
	char c,str[255];
	int i=0;
	while(fgets(str,255,fp) != NULL)
	{
		i+=sizeof(str);  //printf("%s\n",str);
	}//printf("%d\n",i);
	
	rewind(fp);

	JsonData = (char*)calloc(sizeof(char),i+1);
	while(fgets(str,255,fp) != NULL)
		strcat(JsonData,str);
	fclose(fp);
	
	return JsonData;
}

int SetElemType(char *MeshName, int NType)
{
	int kk=0;
	while(1)  // clear '\"' at the head and tail of the json key value
	{
		MeshName[kk]=MeshName[kk+1];
		if(MeshName[kk]=='\"')
		{
			MeshName[kk]='\0';
			break;
		}
		kk++;
	}
	
	char Temp_Type[4];
	int i=0,j=0,k=0;
	while(1)
	{
		Temp_Type[i]=MeshName[j];
		if(MeshName[j]==',' || MeshName[j]=='\0')
		{
			if(Temp_Type[0]=='P'&&Temp_Type[1]=='1')
				ElemType[k]=P1;
			else if(Temp_Type[0]=='L')
			{
				if     (Temp_Type[1]=='2')  ElemType[k+NType]=L2;
				else if(Temp_Type[1]=='3')  ElemType[k+NType]=L3;
				else return 0;
			}
			else if(Temp_Type[0]=='T')
			{
				if     (Temp_Type[1]=='3')  ElemType[k+NType]=T3;
				else if(Temp_Type[1]=='6')  ElemType[k+NType]=T6;
				else return 0;
			}
			else if(Temp_Type[0]=='Q')
			{
				if     (Temp_Type[1]=='4')  ElemType[k+NType]=Q4;
				else if(Temp_Type[1]=='8')  ElemType[k+NType]=Q8;
				else if(Temp_Type[1]=='9')  ElemType[k+NType]=Q9;
				else return 0;
			}
			else if(Temp_Type[0]=='W')
			{
				if     (Temp_Type[1]=='4')  ElemType[k+NType]=W4;
				else if(Temp_Type[1]=='1'
				      &&Temp_Type[2]=='0')  ElemType[k+NType]=W10;
				else return 0;
			}
			else if(Temp_Type[0]=='C')
			{
				if     (Temp_Type[1]=='8')  ElemType[k+NType]=C8;
				else if(Temp_Type[1]=='2'
				      &&Temp_Type[2]=='0')  ElemType[k+NType]=C20;
				else if(Temp_Type[1]=='2'
				      &&Temp_Type[2]=='7')  ElemType[k+NType]=C27;
				else return 0;
			}
			else if(Temp_Type[0]=='H')
			{
				if     (Temp_Type[1]=='6')  ElemType[k+NType]=H6;
				else if(Temp_Type[1]=='1'
				      &&Temp_Type[2]=='5')  ElemType[k+NType]=H15;
				else if(Temp_Type[1]=='1'
				      &&Temp_Type[2]=='8')  ElemType[k+NType]=H18;
				else return 0;
			}
			else return 0;
			i=0;
			k++;
		}
		if(MeshName[j]=='\0')
			break;
		if(MeshName[j]!=',')
			i++;
		j++;
	}
	
	return k+NType;
}

void CreatMeshSpace(char *ElemtType, int NType)
{
	int kk=0;
	while(1)  // clear '\"' at the head and tail of the json key value
	{
		ElemtType[kk]=ElemtType[kk+1];
		if(ElemtType[kk]=='\"')
		{
			ElemtType[kk]='\0';
			break;
		}
		kk++;
	}
	
	char Temp_Type[5];
	int i=0,j=0,k=0;
	while(1)
	{
		Temp_Type[i]=ElemtType[j];
		if(ElemtType[j]==',' || ElemtType[j]=='\0')
		{
			if     (Temp_Type[0]=='N'
				  &&Temp_Type[1]=='o'
				  &&Temp_Type[2]=='d'
				  &&Temp_Type[3]=='e'&&NodMsh!=NULL)
						NodMsh = (Node_Mesh*)malloc(sizeof(Node_Mesh)*NType);
			else if(Temp_Type[0]=='E'
				  &&Temp_Type[1]=='d'
				  &&Temp_Type[2]=='g'
				  &&Temp_Type[3]=='e'&&EdgMsh!=NULL)
						EdgMsh = (Edge_Mesh*)malloc(sizeof(Edge_Mesh)*NType);
			else if(Temp_Type[0]=='F'
				  &&Temp_Type[1]=='a'
				  &&Temp_Type[2]=='c'
				  &&Temp_Type[3]=='e'&&FacMsh!=NULL)
						FacMsh = (Face_Mesh*)malloc(sizeof(Face_Mesh)*NType);
			else if(Temp_Type[0]=='V'
				  &&Temp_Type[1]=='o'
				  &&Temp_Type[2]=='l'
				  &&Temp_Type[3]=='m'&&VolMsh!=NULL)
						VolMsh = (Volm_Mesh*)malloc(sizeof(Volm_Mesh)*NType);
			else return;
			
			i=0;
			k++;
		}
		if(ElemtType[j]=='\0')
			break;
		if(ElemtType[j]!=',')
			i++;
		j++;
	}
	return;
	
}

