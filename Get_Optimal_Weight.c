#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<getopt.h>
#include<time.h>

void GetNormalBiopsySetInfo(FILE*, FILE*);
void GetOddBiopsySetInfo(FILE*, FILE*);
void GetOptimalWeightFileInfo(char*);

int main(int argc,char *argv[])
{
        clock_t begin,end;
        int i,option;
        char *filename=NULL;
        extern char *optarg;
        
        while((option=getopt(argc,argv,"f:h"))!=EOF)
	{
		switch(option)
		{
			case 'f': filename=optarg;
				  printf("\n\tFile Name <%s>\n",filename);
				  break;
			case '?': printf("Unrecognised option\n");
				  exit(1);
                }
        }
	
	begin=clock();
	GetOptimalWeightFileInfo(filename);
	end=clock();
        return 0;
}

void GetNormalBiopsySetInfo(FILE *fr, FILE *fw)
{
	int i;
	char *imagename=(char*)calloc(100,sizeof(char));
	//fprintf(fw,"6\n");
	for(i=0;i<6;i++)
	{
		fscanf(fr,"%s",imagename);
		fprintf(fw,"%s\n",imagename);
	}
	fclose(fw);
	free(imagename);
}

void GetOddBiopsySetInfo(FILE *fr, FILE *fw)
{
	int i;
	char *imagename=(char*)calloc(100,sizeof(char));
	//fprintf(fw,"4\n");
	for(i=0;i<4;i++)
	{
		fscanf(fr,"%s",imagename);
		fprintf(fw,"%s\n",imagename);
	}
	fclose(fw);
	free(imagename);
}

void GetOptimalWeightFileInfo(char *filename)
{
        int i,count;
        FILE *f_read,*f_write,*f_biopsy;
          
	char *path=(char*)calloc(100, sizeof(char));
	sprintf(path,"Optimal_Silhouette.txt");
	f_read=fopen(path,"r");
            
	f_biopsy=fopen(filename,"r");
	
	//fscanf(f_read,"%d",&count);
                
	char *imagepath=(char*)malloc(200 * sizeof(char));
	char *sourcefile=(char*)malloc(200 * sizeof(char));
	char *link=(char*)malloc(200 * sizeof(char));
                
	fscanf(f_biopsy,"%s",sourcefile);
	strcpy(imagepath,sourcefile);
	sprintf(link,"/weights.txt");
	strcat(imagepath,link);
	f_write=fopen(imagepath,"w");
	GetNormalBiopsySetInfo(f_read, f_write);
	
	fscanf(f_biopsy,"%s",sourcefile);
	strcpy(imagepath,sourcefile);
	sprintf(link,"/weights.txt");
	strcat(imagepath,link);
	f_write=fopen(imagepath,"w");
	GetNormalBiopsySetInfo(f_read, f_write);
	
	fscanf(f_biopsy,"%s",sourcefile);
	strcpy(imagepath,sourcefile);
	sprintf(link,"/weights.txt");
	strcat(imagepath,link);
	f_write=fopen(imagepath,"w");
	GetNormalBiopsySetInfo(f_read, f_write);
	
	fscanf(f_biopsy,"%s",sourcefile);
	strcpy(imagepath,sourcefile);
	sprintf(link,"/weights.txt");
	strcat(imagepath,link);
	f_write=fopen(imagepath,"w");
	GetNormalBiopsySetInfo(f_read, f_write);
	
	fscanf(f_biopsy,"%s",sourcefile);
	strcpy(imagepath,sourcefile);
	sprintf(link,"/weights.txt");
	strcat(imagepath,link);
	f_write=fopen(imagepath,"w");
	GetNormalBiopsySetInfo(f_read, f_write);
	
	fscanf(f_biopsy,"%s",sourcefile);
	strcpy(imagepath,sourcefile);
	sprintf(link,"/weights.txt");
	strcat(imagepath,link);
	f_write=fopen(imagepath,"w");
	GetNormalBiopsySetInfo(f_read, f_write);
	
	fscanf(f_biopsy,"%s",sourcefile);
	strcpy(imagepath,sourcefile);
	sprintf(link,"/weights.txt");
	strcat(imagepath,link);
	f_write=fopen(imagepath,"w");
	GetOddBiopsySetInfo(f_read, f_write);
	
	fscanf(f_biopsy,"%s",sourcefile);
	strcpy(imagepath,sourcefile);
	sprintf(link,"/weights.txt");
	strcat(imagepath,link);
	f_write=fopen(imagepath,"w");
	GetNormalBiopsySetInfo(f_read, f_write);
	
	fscanf(f_biopsy,"%s",sourcefile);
	strcpy(imagepath,sourcefile);
	sprintf(link,"/weights.txt");
	strcat(imagepath,link);
	f_write=fopen(imagepath,"w");
	GetNormalBiopsySetInfo(f_read, f_write);
	
	fscanf(f_biopsy,"%s",sourcefile);
	strcpy(imagepath,sourcefile);
	sprintf(link,"/weights.txt");
	strcat(imagepath,link);
	f_write=fopen(imagepath,"w");
	GetNormalBiopsySetInfo(f_read, f_write);
		
                
	fclose(f_biopsy);
	fclose(f_read);
	free(path);
	free(link);
	free(imagepath);
	free(sourcefile);
}