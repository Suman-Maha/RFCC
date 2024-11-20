#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<assert.h>
#include<time.h>
#include<sys/timeb.h>
#include<getopt.h>

#define PI 3.14159265

void write_instruction(void);
void cut_string(char*,int,char*);
int** int_matrix_allocation(int,int);
double** double_matrix_allocation(int,int);
void GetAngularValues(FILE*, FILE*, int*,int, int*);
double GetScalingFactorValue(double, double, double, double);
void Normalizefuzzy_membership(double**, int, int);
double getThresholdValue(double**, int, int, int*, double*, int, double);
void ModifyMembership(double**, int, int, double);
void Updatefuzzy_membership(double*,int*,double*,int,double**,int,double,int*,double*,int, double);
void assign_class(double**,int*,int,int);
void int_matrix_deallocation(int**,int);
void double_matrix_deallocation(double**,int);
double compute_Silhoutte_Index(int*,int*,double*,int*,int,int);
void test_silhouette_index_computation(FILE*,int,int,double,double*,int);

int main(int argc,char *argv[])
{
        clock_t begin,end;
        int i,option,num_of_clusters,maxcount,parameter_length;
        double fuzzifier,elapsed_time;
        extern char *optarg;
	char *filename=NULL;
	char *parameterfile=NULL;
        FILE *f_imagenames,*f_parameter,*f_time;
        
        num_of_clusters=3;
        fuzzifier=2.0;
        
        while((option=getopt(argc,argv,"f:p:c:m:h"))!=EOF)
	{
		switch(option)
		{
			case 'f': filename=optarg;
				  printf("\n\tImage File Name <%s>\n",filename);
				  break;
			case 'p': parameterfile=optarg;
				  printf("\n\tParameter File Name <%s>\n",parameterfile);
				  break;
                        case 'c': num_of_clusters=atoi(optarg);
				  printf("\tNumber of Clusters: %d\n",num_of_clusters);
				  break;
			case 'm': fuzzifier=atof(optarg);
				  printf("\tProbabilistic Weighing Exponent: %lf\n",fuzzifier);
				  break;
			case 'h': write_instruction();
				  break;
			case '?': printf("Unrecognised option\n");
				  exit(1);
		}
	}
	if(filename==NULL)
        {
		write_instruction();
        }
	f_imagenames=fopen(filename,"r");
	if(f_imagenames==NULL)
	{
		printf("\nError: Error in Input Image Details File.\n");
		exit(0);
	}
	fscanf(f_imagenames,"%d",&maxcount);
	
	f_parameter=fopen(parameterfile,"r");
        fscanf(f_parameter,"%d",&parameter_length);
	double *weights=(double*)calloc(parameter_length,sizeof(double));
	
	for(i=0;i<parameter_length;i++)
	{
		 fscanf(f_parameter,"%lf",&weights[i]);
	}

	fclose(f_parameter);
	
        begin=clock();
	test_silhouette_index_computation(f_imagenames,maxcount,num_of_clusters,fuzzifier,weights,parameter_length);
	end=clock();
	free(weights);
	elapsed_time=(double)(end-begin)/CLOCKS_PER_SEC;
	f_time=fopen("Execution_time_Silhouette.txt","w");
	printf("\n******************** TOTAL TIME ELAPSED : %lf SECONDS ********************\n\n",elapsed_time);
	fprintf(f_time,"\n\n********** TIME ELAPSED FOR SILHOUETTE INDEX COMPUTATION: %lf SECONDS **********\n\n",elapsed_time);
	fclose(f_time);
	return 0;
}

void write_instruction(void)
{
	system("clear");
	printf("f:\tInput Image Details File\n");
	printf("p:\tParameter Set File\n");
	printf("c:\tNumber of Clusters\n");
	printf("m:\tProbabilistic Weighting Exponent\n");
	printf("h:\tHelp\n");
	exit(1);
}

void cut_string(char *str1,int index,char *str2)
{
	int i;
	for(i=0;i<index;i++)
	{
		str2[i]=str1[i];
	}
	str2[i]='\0';
}

int** int_matrix_allocation(int row,int column)
{
        int i,**data;
        data=(int**)calloc(row,sizeof(int*));
        assert(data!=NULL);
        for(i=0;i<row;i++)
        {
                data[i]=(int*)calloc(column,sizeof(int));
                assert(data[i]!=NULL);
        }
        return(data);
}

double** double_matrix_allocation(int row,int column)
{
        int i;
        double **data;
        data=(double**)calloc(row,sizeof(double*));
        assert(data!=NULL);
        for(i=0;i<row;i++)
        {
                data[i]=(double*)calloc(column,sizeof(double));
                assert(data[i]!=NULL);
        }
        return(data);
}

void GetAngularValues(FILE *f_SW, FILE *f_Modified, int *AngularData,int AngleLength, int *AngleFreq)
{
        int i,freq1;
	double freq2;
        for(i=0;i<AngleLength;i++)
        {
		AngularData[i]=i;
		fscanf(f_Modified,"%d",&freq1);
		fscanf(f_SW,"%lf",&freq2);
                AngleFreq[i]=(int)((double)freq1+freq2+0.5);
        }
        fclose(f_SW);
	fclose(f_Modified);
}

double GetScalingFactorValue(double lower, double upper, double interval, double concentration)
{
        int i,n;
        double so,se,ans;
        double *x=(double*)calloc(5000, sizeof(double));
        double *y=(double*)calloc(5000, sizeof(double));
        n=(upper-lower)/interval;
        if(n%2==1)
        {
                n=n+1;
        }
        interval=(upper-lower)/n;
        for(i=0; i<=n; i++)
        {
                x[i]=lower+(i*interval);
                y[i]=exp(concentration*cos(x[i]));
        }
        so=0;
        se=0;
        for(i=1; i<n; i++)
        {
                if(i%2==1)
                {
                        so=so+y[i];
                }
                else
                {
                        se=se+y[i];
                }
 
        }
        ans=interval/3*(y[0]+y[n]+4*so+2*se);
        free(x);
        free(y);
        return(ans);
}

void Normalizefuzzy_membership(double **fuzzy_membership, int num_of_clusters, int AngleLength)
{
        int i,j;
        double sum;
        for(j=0;j<AngleLength;j++)
        {
                sum=0.00;
                for(i=0;i<num_of_clusters;i++)
                {
                        sum+=fuzzy_membership[i][j];
                }
                for(i=0;i<num_of_clusters;i++)
                {
                        fuzzy_membership[i][j]/=sum;
                }
        }
}

double getThresholdValue(double **fuzzy_membership, int num_of_clusters, int AngleLength, int *Modified_hue_histogram, double *SW_hue_histogram, int AngleFreqSum, double SWFreqSum)
{
	int i,j,index;
	double max1,max2,threshold,sum;
	sum=0.00;
	for(j=0;j<AngleLength;j++)
	{
		max1=0.00;
		for(i=0;i<num_of_clusters;i++)
                {
                        if(fuzzy_membership[i][j]>=max1)
                        {
                                max1=fuzzy_membership[i][j];
                                index=i;
                        }
                }
                max2=0.00;
                for(i=0;i<num_of_clusters;i++)
                {
                        if(fuzzy_membership[i][j]>=max2 && i!=index)
                        {
                                max2=fuzzy_membership[i][j];
                        }
                }
                sum+=(fabs(max1-max2)*((double)Modified_hue_histogram[j]+SW_hue_histogram[j]));
	}
	threshold=((double)sum/((double)AngleFreqSum+SWFreqSum));
	return threshold;
}

void ModifyMembership(double **fuzzy_membership, int num_of_clusters, int AngleLength, double threshold)
{
        int i,j,index;
        double max1, max2;
        for(j=0;j<AngleLength;j++)
        {
                max1=0.00;
                for(i=0;i<num_of_clusters;i++)
                {
                        if(fuzzy_membership[i][j]>=max1)
                        {
                                max1=fuzzy_membership[i][j];
                                index=i;
                        }
                }
                max2=0.00;
                for(i=0;i<num_of_clusters;i++)
                {
                        if(fuzzy_membership[i][j]>=max2 && i!=index)
                        {
                                max2=fuzzy_membership[i][j];
                        }
                }
                if((max1-max2)>threshold)
                {
                        fuzzy_membership[index][j]=1.00;
                        for(i=0;i<num_of_clusters;i++)
                        {
                                if(i!=index)
                                {
                                        fuzzy_membership[i][j]=0.00;
                                }
                        }
                }
        }
}

void Updatefuzzy_membership(double *Concentration,int *AngularData,double *ClusterCentroid,int num_of_clusters,double **fuzzy_membership,int AngleLength,double fuzzifier,int *Modified_hue_histogram,double *SW_hue_histogram,int AngleFreqSum, double SWFreqSum)
{
        int i,j,count;
        double MembershipDenominator, MembershipNumerator, factor, threshold, lower, upper, interval, sum;
        factor=((double)1/fuzzifier);
        lower=0.0, upper=(2*PI), interval=0.1;
        for(i=0;i<num_of_clusters;i++)
        {
                MembershipDenominator=GetScalingFactorValue(lower, upper, interval, Concentration[i]);
                for(j=0;j<AngleLength;j++)
                {
                        MembershipNumerator=exp(Concentration[i]*cos((AngularData[j]-ClusterCentroid[i])*(PI/180.0)));
                        fuzzy_membership[i][j]=pow(((double)MembershipNumerator/MembershipDenominator),factor);
                }
        }
        Normalizefuzzy_membership(fuzzy_membership, num_of_clusters, AngleLength);
	
	for(j=0;j<AngleLength;j++)
	{
		count=0;
		for(i=0;i<num_of_clusters;i++)
		{
			if(fuzzy_membership[i][j]==1.0)
			{
				count++;
			}
		}
		
		if(count>1)
		{
			printf("\nError: Error in Program\n");
			exit(0);
		}
		
		else if(count==1)
		{
			for(i=0;i<num_of_clusters;i++)
			{
				if(fuzzy_membership[i][j]==1.0)
				{
					fuzzy_membership[i][j]=1.0;
				}
				else
				{
					fuzzy_membership[i][j]=0.0;
				}
			}
		}
		
		else if(count==0)
		{
			sum=0.0;
			for(i=0;i<num_of_clusters;i++)
			{
				sum+=fuzzy_membership[i][j];
			}
			for(i=0;i<num_of_clusters;i++)
			{
				fuzzy_membership[i][j]=(double)fuzzy_membership[i][j]/sum;
			}
		}
	}
}

void assign_class(double **fuzzy_membership,int *ClusterLabel,int num_of_clusters,int AngleLength)
{
        int i,j,index;
        double max;
        
        for(j=0;j<AngleLength;j++)
	{
		max=0.0;
                index=0;
		for(i=0;i<num_of_clusters;i++)
                {
                        if(fuzzy_membership[i][j]>=max)
                        {
                                max=fuzzy_membership[i][j];
                                index=i;
                        }
                }
                ClusterLabel[j]=index;
        }
}

void int_matrix_deallocation(int **matrix,int row)
{
	int i;
	for(i=0;i<row;i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}

void double_matrix_deallocation(double **matrix,int row)
{
	int i;
	for(i=0;i<row;i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}

double compute_Silhoutte_Index(int *AngularData,int *ClusterLabel,double *Concentration,int *AngleFreq,int AngleLength,int num_of_clusters)
{
        int i,j,k,*count,**freq;
        double max,*avg,*a,*b,*silhouette_width,*silhouette,*silhouette_cluster,total_silhouette,global_silhouette,lower, upper, interval, ScaleFactor;
        
        freq=int_matrix_allocation(num_of_clusters,AngleLength);
        a=(double*)calloc(AngleLength,sizeof(double));
        b=(double*)calloc(AngleLength,sizeof(double));
        silhouette_width=(double*)calloc(AngleLength,sizeof(double));
        
	lower=0, upper=(2*PI), interval=0.1;
	
        for(i=0;i<AngleLength;i++)
                freq[ClusterLabel[i]][AngularData[i]]=AngleFreq[i];
        
        for(i=0;i<AngleLength;i++)
        {
                avg=(double*)calloc(num_of_clusters,sizeof(double));
                count=(int*)calloc(num_of_clusters,sizeof(int));
                
                for(j=0;j<num_of_clusters;j++)
                {
			ScaleFactor=GetScalingFactorValue(lower, upper, interval, Concentration[j]);
                        if(ClusterLabel[i]==j)
                        {
                                for(k=0;k<AngleLength;k++)
                                {
                                        if(AngularData[i]==k)
                                        {
						avg[j]+=fabs(log(ScaleFactor)-(Concentration[j]*cos(((AngularData[i]-k)*PI)/180.0)))*(freq[j][k]-1);
                                                count[j]+=(freq[j][k]-1);
                                        }
                                        else
                                        {
                                                avg[j]+=fabs(log(ScaleFactor)-(Concentration[j]*cos(((AngularData[i]-k)*PI)/180.0)))*freq[j][k];
                                                count[j]+=freq[j][k];
                                        }
                                }
                        }
                        else
                        {
                                for(k=0;k<=AngleLength;k++)
                                {
                                        avg[j]+=fabs(log(ScaleFactor)-(Concentration[j]*cos(((AngularData[i]-k)*PI)/180.0)))*freq[j][k];
                                        count[j]+=freq[j][k];
                                }
                        }
                        
                        if(count[j]>0)
                                avg[j]=(avg[j]/(double)count[j]);
                        else
                                avg[j]=1;
                }
                
                a[i]=avg[ClusterLabel[i]];
                
                if(ClusterLabel[i]==0)
                        b[i]=avg[1];
                else
                        b[i]=avg[0];
                
                for(j=0;j<num_of_clusters;j++)
                {
                        if(j!=ClusterLabel[i])
                        {
                                if(b[i]>avg[j])
                                {
                                        b[i]=avg[j];
                                }
                        }
                }
                                
                if(a[i]>=b[i])
                        max=a[i];
                else
                        max=b[i];
                
                silhouette_width[i]=(b[i]-a[i])/(double)max;
                
                free(avg);
                free(count);
        }
        silhouette=(double*)calloc(num_of_clusters,sizeof(double));
	silhouette_cluster=(double*)calloc(num_of_clusters,sizeof(double));
	count=(int*)calloc(num_of_clusters,sizeof(int));
        
        for(i=0;i<AngleLength;i++)
        {
                silhouette[ClusterLabel[i]]+=(silhouette_width[i]*AngleFreq[i]);
                count[ClusterLabel[i]]+=AngleFreq[i];
        }
	
	total_silhouette=0.0;
        for(j=0;j<num_of_clusters;j++)
	{
		if(count[j]>0)
			silhouette_cluster[j]=(silhouette[j]/(double)count[j]);
		total_silhouette+=silhouette_cluster[j];
	}
	global_silhouette=((double)total_silhouette/num_of_clusters);
        
	free(silhouette_cluster);
        free(silhouette);
	free(count);
        free(a);
        free(b);
        free(silhouette_width);
        int_matrix_deallocation(freq,num_of_clusters);
        
        return global_silhouette;
}

void test_silhouette_index_computation(FILE *f_imagenames,int maxcount,int num_of_clusters,double fuzzifier,double *weights,int parameter_length)
{
        int i,j,num_of_samples,count, AngleLength, AngleFreqSum ;
        int m,n;
        double global_silhouette, center_achrom, conc_achrom, maximum_silhouette, weight, SWFreqSum;
	int *Modified_hue_histogram;
	double *SW_hue_histogram;
        double **fuzzy_membership;
        count=0,AngleLength=360;
        FILE *f_Silhoette,*f_upto, *f_Silhouette_optimize;
        f_Silhoette=fopen("Global_Silhouette.txt","w");
	f_Silhouette_optimize=fopen("Optimal_Silhouette.txt","w");
	
        while(count<maxcount)
        {
                count++;
                printf("\n\n********************Opeartion on Image %d is started********************\n", count);
                char *imagename=(char*)calloc(100,sizeof(char));
                char *imagepath=(char*)calloc(200,sizeof(char));
                char *sourcefile=(char*)calloc(100,sizeof(char));
		Modified_hue_histogram=(int*)calloc(AngleLength,sizeof(int));
                SW_hue_histogram=(double*)calloc(AngleLength,sizeof(double));
		
		FILE *f_hueHist,*f_center,*f_s_index, *f_conc, *f_Mod_hist, *f_SWHue, *f_ModHue;
		
                fscanf(f_imagenames,"%s",sourcefile);
                
                cut_string(sourcefile, (strlen(sourcefile)-4), imagename);
        
                strcpy(imagepath,imagename);
                strcat(imagepath, "/SW_hue_hist.txt");
                f_hueHist=fopen(imagepath, "r");
		
		strcpy(imagepath,imagename);
                strcat(imagepath, "/Modified_hue_hist.txt");
                f_Mod_hist=fopen(imagepath, "r");
		
                int *AngularData=(int*)calloc(AngleLength,sizeof(int));
                int *AngleFreq=(int*)calloc(AngleLength,sizeof(int));
		GetAngularValues(f_hueHist, f_Mod_hist, AngularData, AngleLength, AngleFreq);
		
		strcpy(imagepath,imagename);
                strcat(imagepath, "/SW_hue_hist.txt");
                f_SWHue=fopen(imagepath, "r");
		for(i=0;i<AngleLength;i++)
		{
			fscanf(f_SWHue,"%lf",&SW_hue_histogram[i]);
		}
		fclose(f_SWHue);
		SWFreqSum=0.0;
		for(i=0;i<AngleLength;i++)
		{
			SWFreqSum+=SW_hue_histogram[i];
		}

		strcpy(imagepath,imagename);
                strcat(imagepath, "/Modified_hue_hist.txt");
                f_ModHue=fopen(imagepath, "r");
		for(i=0;i<AngleLength;i++)
		{
			fscanf(f_ModHue,"%d",&Modified_hue_histogram[i]);
		}
		fclose(f_ModHue);
		AngleFreqSum=0;
		for(i=0;i<AngleLength;i++)
		{
			AngleFreqSum+=Modified_hue_histogram[i];
		}

		strcpy(imagepath,imagename);
                strcat(imagepath,"/SILHOUETTE.txt");
                f_s_index=fopen(imagepath, "w");

		maximum_silhouette=0.0;
		
		for(j=0;j<parameter_length;j++)
		{
                        fuzzy_membership=double_matrix_allocation(num_of_clusters,AngleLength);
                                
                        char *path=(char*)calloc(200,sizeof(char));
                        double *ClusterCentroid=(double*)calloc(num_of_clusters,sizeof(double));
			double *Concentration=(double*)calloc(num_of_clusters, sizeof(double));
                        int *ClusterLabel=(int*)calloc(AngleLength,sizeof(int));
                                
                        sprintf(path,"/Centroid_w_%0.2lf.txt",weights[j]);
                        strcpy(imagepath,imagename);
                        strcat(imagepath,path);
                        f_center=fopen(imagepath,"r");
			
			sprintf(path,"/kappa_w_%0.2lf.txt",weights[j]);
                        strcpy(imagepath,imagename);
                        strcat(imagepath,path);
                        f_conc=fopen(imagepath,"r");
			
			fscanf(f_center,"%lf  %lf  %lf",&ClusterCentroid[0],&ClusterCentroid[1],&ClusterCentroid[2]);
			fclose(f_center);

			fscanf(f_conc,"%lf  %lf  %lf",&Concentration[0],&Concentration[1],&Concentration[2]);
			fclose(f_conc);
				
			Updatefuzzy_membership(Concentration,AngularData,ClusterCentroid,num_of_clusters,fuzzy_membership,AngleLength,fuzzifier,Modified_hue_histogram,SW_hue_histogram,AngleFreqSum,SWFreqSum);
                        assign_class(fuzzy_membership,ClusterLabel,num_of_clusters,AngleLength);
			global_silhouette=compute_Silhoutte_Index(AngularData,ClusterLabel,Concentration,AngleFreq,AngleLength,num_of_clusters);
			
			if(global_silhouette>=maximum_silhouette)
			{
				maximum_silhouette=global_silhouette;
				weight=weights[j];
			}
			
                        fprintf(f_s_index,"%lf ",global_silhouette);
                        fprintf(f_Silhoette,"%lf ",global_silhouette);
                                
                        free(path);
                        free(ClusterLabel);
                        free(ClusterCentroid);
			free(Concentration);
                        double_matrix_deallocation(fuzzy_membership,num_of_clusters);
		}
		
		fprintf(f_Silhoette,"\n");
		fprintf(f_Silhouette_optimize,"%0.2lf\n",weight);
		fclose(f_s_index);
		free(AngularData);
		free(AngleFreq);
		free(imagename);
		free(imagepath);
		free(sourcefile);
		free(SW_hue_histogram);
		free(Modified_hue_histogram);
                printf("\n\n********************Opeartion on Image %d is completed********************\n", count);
		
		f_upto=fopen("Upto_Image_Silhoette.txt","w");
		fprintf(f_upto,"Silhoette Computation done Upto Image: %d", count);
		fclose(f_upto);
        }
        fclose(f_Silhouette_optimize);
        fclose(f_Silhoette);
}