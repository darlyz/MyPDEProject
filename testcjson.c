#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<regex.h>
#include"cJSON.h"

int main(int argc, char* argv[])
{
	FILE *fp;
	if((fp=fopen(argv[1],"r"))==NULL)
		return -1;
	
	char str[255];
	int i=0;
	while(fgets(str,255,fp) != NULL)
	{
		i+=sizeof(str);  //printf("%s\n",str);
		
	}//printf("%d\n",i);
	
	rewind(fp);

	char *data;
	data = (char*)calloc(sizeof(char),i+1);
	while(fgets(str,255,fp) != NULL)
		strcat(data,str);
	fclose(fp);
	//printf("%d\t%s\n",i,data);
	
    //char *data = "{\"love\":[\"LOL\",\"Go shopping\"]}";
    //从缓冲区中解析出JSON结构
    cJSON * json= cJSON_Parse(data);
    
    //将传入的JSON结构转化为字符串 并打印
    char *json_data = NULL;
    printf("data:%s\n",json_data = cJSON_Print(json));
	
	cJSON *RegionMesh = cJSON_GetObjectItem(json,"RegionMesh");
	cJSON *BndaryMesh = cJSON_GetObjectItem(json,"BndaryMesh");
	char  *RegMesh    = cJSON_Print(RegionMesh);
	char  *BndMesh    = cJSON_Print(BndaryMesh);
	printf("RegionMesh:%s\n",RegMesh);
	printf("BndaryMesh:%s\n",BndMesh);
	
	char *pattern="[A-Z][0-9]+";
	
    
    free(json_data);
    //将JSON结构所占用的数据空间释放
    cJSON_Delete(json);
    return 0;
}