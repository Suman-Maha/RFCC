#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<assert.h>
#include<time.h>
#include<getopt.h>
#include<sys/stat.h>

#define PI 3.14159265

struct RGB
{
        double red;
        double green;
        double blue;
};

void write_instruction(void);
double FindMaximum(double, double);
double FindMinimum(double, double);
double findMinimum(double, double, double);
int** integer_matrix_allocation(int, int);
double** double_matrix_allocation(int, int);
void integer_matrix_deallocation(int**, int);
void double_matrix_deallocation(double**, int);
void cut_string(char*, int, char*);
void skip_comments(FILE*);
void read_ppm_write_pgmMatrix(FILE*, int**, int**, int**, double**, double**, double**, int, int);
void nonlinearRGB_to_linearRGB(double**, double**, double**, int, int, double);
void RGB_to_HSI_conversion(double**, double**, double**, double**, double**, double**, int, int);
void compute_SW_hue_histograms(double**, double**, int**, double**, int, int, double*);
void GetSmoothHistogram(FILE*, double*, double*);
void ComputeNormalizedHistogram(double*, double*, int, int);
double FindMedianValue(double*, int);
int ThresholdingBackward(double*, int);
int ThresholdingForward(double*, int);
double ComputeInitialCentroid(double*, int, int);
double ComputeMeanResultantDirection(double*, int, int);
double ComputeInitialConcentration(double);
void GetAngularValues(double*, int);
void CopyOriginalArray(double*, double*, int);
void CopyOriginalMatrix(double**, double**, int, int);
double GetScalingFactorValue(double, double, double, double);
void NormalizeMembershipValue(double**, int, int);
double getThresholdValue(double**, int, int, int*, double*, int, double);
void UpdateMembershipValue(double*, double*, double*, int, double**,int, double, int*, double*, int, double);
void find_roughsets(double**,int,int,int*,int**,double);
double GetActualAngle(double,double);
void estimate_centroid(double**,double*,double*,int,int,double,int*,int**,int*,double*,double);
double extract_kappa_value(double);
void estimate_concentration_kappa(double**,double*,double*,double*,int,int,double,int*,int**,int*,double*,int);
void getHSIInformation(FILE*, int, int, double**);
void findClusterCenterSaturation(double**, double*, int, int, double*, double*);
struct RGB RG_Sector_Conversion(double, double, double);
struct RGB GB_Sector_Conversion(double, double, double);
struct RGB BR_Sector_Conversion(double, double, double);
void HSI_to_RGB_conversion(double*, double*, double*, double*, double*, double*, int);
void convert_to_log_domain(double*, double*, double*, int);
double zero_check_normalization(int);
void read_ppm_write_pgmMatrix_Mapping(FILE*, int**, int**, int**, double**, double**, double**, int, int);
void GetResizedImageMatrix(double**, double**, double**, int, int, double**);
double** GetMinor(double**,int,int,int);
double Determinant_value(double**,int);
void transpose(double**,double**, int, int);
void GetInverse(double**,double**,int);
void MatrixMultiplication(double**, double**, int, int, int, double**);
void UpdateFirstMatrix(double**, double**, double**, int, int, int);
void UpdateSecondMatrix(double**, double**, double**, int, int, int);
int checkifLargeDifferenceArray(double*, double*, int, double);
int checkifLargeDifferenceMatrix(double**, double**, int, int, double);
void MapConcentrationDensity(double**,double**,int,int,int,int,int,double**);
void ReformImageChannels(double**,int**,int**,int**,int,int);
void read_pgmMatrix_write_ppm(int**, int**, int**, FILE*, int, int, int);
void Test_SW_hue_histogram(FILE*, char*, int, int, double);
void TestInitialization(FILE*,int);
void TestCircularClustering(FILE*, int, int, double*, int, double);
void TestfindClusterHSI(FILE*,int,int,double*,int);
void TestMappingProcedure(FILE*,char*,int,double*,int,int,char*);

int main(int argc,char *argv[])
{
        clock_t begin,end;
        int i,option, num_of_clusters, max_iteration, maxcount,parameter_length,rectsize;
        double fuzzifier, elapsed_time,concentration,alpha;
	double *SW_hue_histogram;
        extern char *optarg;
	char *filename=NULL,*parameterfile=NULL,*TargetImage=NULL, *filepath;
        FILE *f_imagenames,*f_parameter, *f_time;
        num_of_clusters=3;
	max_iteration=100;
        fuzzifier=2.0;
	rectsize=3;
        while((option=getopt(argc,argv,"f:c:p:m:s:a:r:h"))!=EOF)
	{
		switch(option)
		{
			case 'f': filename=optarg;
				  printf("\n\tFile Name <%s>\n",filename);
				  break;
			case 'c': num_of_clusters=atoi(optarg);
				  printf("\n\tNumber of Clusters: %d\n",num_of_clusters);
				  break;
                        case 'p': parameterfile=optarg;
				  printf("\n\tParameter File Name <%s>\n",parameterfile);
				  break;
			case 'm': fuzzifier=atof(optarg);
				  printf("\n\tFuzzifier Value: %lf\n",fuzzifier);
				  break;
			case 's': rectsize=atoi(optarg);
			          printf("\n\tRectangle Size: %dx%d\n", rectsize,rectsize);
				  break;
			case 'a': alpha=atof(optarg);
			          printf("\n\tAlpha Parameter Value: %lf\n", alpha);
				  break;
                        case 'r': TargetImage=optarg;
                                  printf("\n\tReference Imagename <%s>\n",TargetImage);
                                  break;
			case 'h': write_instruction();
				  break;
			case '?': printf("Unrecognised option\n");
				  exit(1);
		}
	}
        f_parameter=fopen(parameterfile,"r");
        fscanf(f_parameter,"%d",&parameter_length);
	double *weights=(double*)calloc(parameter_length,sizeof(double));
	for(i=0;i<parameter_length;i++)
	{
		 fscanf(f_parameter,"%lf",&weights[i]);
	}
	fclose(f_parameter);
	
	filepath=(char*)calloc(200, sizeof(char));
	filepath="/home/suman/Suman_Backup/IMPT/Different_Formulations/UCSB_Data/UCSB_Images/";
	filepath="./ImageFolder/";
	
        begin=clock();
	f_imagenames=fopen(filename,"r");
        fscanf(f_imagenames,"%d",&maxcount);
	Test_SW_hue_histogram(f_imagenames, filepath, maxcount, rectsize, alpha);
	fclose(f_imagenames);
	f_imagenames=fopen(filename,"r");
        fscanf(f_imagenames,"%d",&maxcount);
	TestInitialization(f_imagenames, maxcount);
	f_imagenames=fopen(filename,"r");
        fscanf(f_imagenames,"%d",&maxcount);
	TestCircularClustering(f_imagenames, maxcount, num_of_clusters, weights, parameter_length, fuzzifier);
        f_imagenames=fopen(filename,"r");
        fscanf(f_imagenames,"%d",&maxcount);
	TestfindClusterHSI(f_imagenames,maxcount,num_of_clusters,weights,parameter_length);	
        fclose(f_imagenames);
        f_imagenames=fopen(filename,"r");
        fscanf(f_imagenames,"%d",&maxcount);
	TestMappingProcedure(f_imagenames,filepath,maxcount,weights,parameter_length,num_of_clusters,TargetImage);
        fclose(f_imagenames);
	end=clock();
        free(weights);
	elapsed_time=(double)(end-begin)/CLOCKS_PER_SEC;
	f_time=fopen("Execution_time.txt","w");
	printf("\n******************** TOTAL TIME ELAPSED : %lf SECONDS ********************\n\n",elapsed_time);
	fprintf(f_time,"\n\n********** TOTAL TIME ELAPSED FOR COLOR NORMALIZATION: %lf SECONDS **********\n\n",elapsed_time);
	fclose(f_time);
	return 0;
}

void write_instruction(void)
{
	printf("\n\tf: Image File Name\n");
	printf("\tc: Number of Clusters\n");
	printf("\tp: Parameter File Name\n");
	printf("\tm: Fuzzifier Under Consideration\n");
	printf("\tr: Reference Image Filename\n");
	printf("\th: Help for the instructions\n");
}

double FindMaximum(double a, double b)
{
	if(a>=b)
		return a;
	else 
		return b;
}

double FindMinimum(double a, double b)
{
	if(a<=b)
		return a;
	else 
		return b;
}

double findMinimum(double a, double b, double c)
{
	double min;
	if(a<=b)
		min=a;
	else 
		min=b;
	if(c<=min)
		min=c;
	return min;
}

int** integer_matrix_allocation(int row,int column)
{
        int i;
        int **matrix;
        matrix=(int**)calloc(row, sizeof(int*));
        assert(matrix!=NULL);
        for(i=0;i<row;i++)
        {
                matrix[i]=(int*)calloc(column, sizeof(int));
                assert(matrix[i]!=NULL);
        }
        return matrix;
}

double** double_matrix_allocation(int row, int col)
{
	int i;
	double **matrix;
        matrix=(double**)calloc(row, sizeof(double*));
        assert(matrix!=NULL);
	for(i=0;i<row;i++)
	{
		matrix[i]=(double*)calloc(col, sizeof(double));
                assert(matrix[i]!=NULL);
	}
	return matrix;
}

void integer_matrix_deallocation(int **matrix, int row)
{
	int i;
	for(i=0;i<row;i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}

void double_matrix_deallocation(double **matrix, int row)
{
	int i;
	for(i=0;i<row;i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}

void cut_string(char *str1, int index, char *str2)
{
	int i;
	for(i=0;i<index;i++)
	{
		str2[i]=str1[i];
	}
	str2[i]='\0';
}

void skip_comments(FILE *fp)
{
	char c;
	c=getc(fp);
	if(c!='#')
	{
		ungetc(c,fp);
	}
	else
	{
		do
		{
			c=getc(fp);
		}while(c!='\n');
		skip_comments(fp);		
	}
}

void read_ppm_write_pgmMatrix(FILE *fp, int **image_red, int **image_green, int **image_blue, double **normalized_image_red, double **normalized_image_green, double **normalized_image_blue, int col, int row)
{
	int i,j,red,green,blue;
	char ch;
	fscanf(fp,"%c",&ch);
	for(i = 0; i < row; i++)
	{
		for(j=0;j< col;j++)
		{
			red=fgetc(fp);
			image_red[i][j]=red;
			green=fgetc(fp);
			image_green[i][j]=green;
			blue=fgetc(fp);
			image_blue[i][j]=blue;
		}
	}
        for(i=0;i<row;i++)
        {
                for(j=0;j<col;j++)
                {
                        normalized_image_red[i][j]=((double)image_red[i][j]/255.0);
                        normalized_image_green[i][j]=((double)image_green[i][j]/255.0);
                        normalized_image_blue[i][j]=((double)image_blue[i][j]/255.0);
                }
        }
	fclose(fp);
}

void nonlinearRGB_to_linearRGB(double **normalized_image_red, double **normalized_image_green, double **normalized_image_blue, int col, int row, double gamma)
{
	int i,j;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			if(normalized_image_red[i][j]<=0.018)
			{
				normalized_image_red[i][j]=((double)normalized_image_red[i][j]/4.5);
			}
			else
			{
				normalized_image_red[i][j]=pow(((normalized_image_red[i][j]+0.099)/1.099),gamma);
			}
			if(normalized_image_green[i][j]<=0.018)
			{
				normalized_image_green[i][j]=((double)normalized_image_green[i][j]/4.5);
			}
			else
			{
				normalized_image_green[i][j]=pow(((normalized_image_green[i][j]+0.099)/1.099),gamma);
			}
			if(normalized_image_blue[i][j]<=0.018)
			{
				normalized_image_blue[i][j]=((double)normalized_image_blue[i][j]/4.5);
			}
			else
			{
				normalized_image_blue[i][j]=pow(((normalized_image_blue[i][j]+0.099)/1.099),gamma);
			}
		}
	}
}

void RGB_to_HSI_conversion(double **normalized_image_red, double **normalized_image_green, double **normalized_image_blue, double **theta, double **saturation, double **intensity, int col, int row)
{
	int i,j;
	double hue_numer,hue_deno,minvalue,theta_value;
	double pi=3.1415926535;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			if((normalized_image_red[i][j]==normalized_image_green[i][j]) && (normalized_image_red[i][j]==normalized_image_blue[i][j]))
			{
				theta[i][j]=((double)90/360);
				continue;
			}
			hue_numer=0.5*((normalized_image_red[i][j]-normalized_image_green[i][j])+(normalized_image_red[i][j]-normalized_image_blue[i][j]));
			hue_deno=sqrt((pow((normalized_image_red[i][j]-normalized_image_green[i][j]),2)+((normalized_image_red[i][j]-normalized_image_blue[i][j])*(normalized_image_green[i][j]-normalized_image_blue[i][j]))));
			theta_value=((double)(acos(hue_numer/hue_deno)*180)/pi);
			if(normalized_image_blue[i][j]<=normalized_image_green[i][j])
			{
				theta[i][j]=((double)theta_value/360);
			}
			else
			{
				theta[i][j]=((double)(360-theta_value)/360);
			}
			minvalue=findMinimum(normalized_image_red[i][j],normalized_image_green[i][j],normalized_image_blue[i][j]);
			saturation[i][j]=(1.00-((3.00*minvalue)/(normalized_image_red[i][j]+normalized_image_green[i][j]+normalized_image_blue[i][j])));
			intensity[i][j]=((normalized_image_red[i][j]+normalized_image_green[i][j]+normalized_image_blue[i][j])/3.00);
		}
	}
}

void Compute_Neighborhood_Average(double **hue_actual, int row, int col, int rectsize, double alpha, int **Modified_hue, int *Modified_hue_histogram)
{
	int i,j,rr,cc,count;
	double neighbor_gross_hue,neighbor_average_hue,hue_average;
	count=0;
	neighbor_gross_hue=0.0;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			for(rr=(i-(rectsize/2));rr<((i-(rectsize/2))+rectsize);rr++)
			{
				for(cc=(j-(rectsize/2));cc<((j-(rectsize/2))+rectsize);cc++)
				{
					if(rr>=0 && rr<row && cc>=0 && cc<col)
					{
						if(rr==i && cc==j);
						else
						{
							count++;
							neighbor_gross_hue+=hue_actual[rr][cc];
						}
					}
				}
			}
			neighbor_average_hue=((((double)alpha/(double)count)*neighbor_gross_hue)+hue_actual[i][j]);
			hue_average=((double)neighbor_average_hue/(double)(1.0+alpha));
			Modified_hue[i][j]=((int)(hue_average+0.5));
		}
	}
	
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			Modified_hue_histogram[Modified_hue[i][j]]++;
		}
	}
}

void compute_SW_hue_histograms(double **theta, double **saturation, int **hue, double **hue_actual, int row, int col, double *SW_hue_histogram)
{
	int i,j,bin=360;
        for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			hue_actual[i][j]=(bin*theta[i][j]);
		}
	}
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			hue[i][j]=(int)(hue_actual[i][j]+0.5);
		}
	}
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			SW_hue_histogram[hue[i][j]]=SW_hue_histogram[hue[i][j]]+saturation[i][j];
		}
	}
}

void GetSmoothHistogram(FILE *fp, double *SW_hue_histogram, double *histogram_smooth)
{
	int i,size;
	double *array;
	size=3;
	array=(double*)calloc(size,sizeof(double));
	histogram_smooth[0]=FindMinimum(SW_hue_histogram[0],SW_hue_histogram[1]);
	histogram_smooth[359]=FindMaximum(SW_hue_histogram[358],SW_hue_histogram[359]);
	for(i=1;i<359;i++)
	{
		array[0]=SW_hue_histogram[i-1];
		array[1]=SW_hue_histogram[i];
		array[2]=SW_hue_histogram[i+1];
		histogram_smooth[i]=FindMedianValue(array, size);
	}
	free(array);
	for(i=0;i<360;i++)
	{
		fprintf(fp,"%d\t%lf\n",i,histogram_smooth[i]);
	}
	fclose(fp);
}

void ComputeNormalizedHistogram(double *SW_hue_histogram, double *histogram_norm, int start_angle, int end_angle)
{
        int i;
	double sum;
        sum=0.0;
        for(i=start_angle;i<end_angle;i++)
        {
                sum+=SW_hue_histogram[i];
        }
        for(i=start_angle;i<end_angle;i++)
        {
                histogram_norm[i-start_angle]=((double)SW_hue_histogram[i]/sum);
        }
}

double FindMedianValue(double *array, int size)
{
	int i,j;
	double temp;
	for(i=0;i<size;i++)
	{
		for(j=i+1;j<size;j++)
		{
			if(array[i]>array[j])
			{
				temp=array[i];
				array[i]=array[j];
				array[j]=temp;
			}
		}
	}
	return array[1];
}

int ThresholdingBackward(double *hist, int end_angle)
{
        int i,t,tt,thr;
        double min_sigma_W, sum_sqr_1, sum_sqr_2, sum_1, sum_2, sigma_1, omega_1, mean_1, sigma_2, omega_2, mean_2, sigma_W;
        thr = -1;
        min_sigma_W = 9e20;
        sum_sqr_1 = sum_1 = 0;
        sum_sqr_2 = sum_2 = 0;
        for (t = (end_angle-1); t >= (end_angle/2); t--) 
        {
                if (t == (end_angle-1)) 
                {
                        sigma_1 = omega_1 = mean_1 = 0;
                        for (i = (end_angle-1); i >= (end_angle/2); i--) 
                        {
                                omega_1 += hist[i];
                                mean_1 += i * hist[i];
                                sum_sqr_1 += pow(i,2) * hist[i];
                                sum_1 += i * hist[i];
                        }
                        mean_1 /= omega_1;
                        for (i = (end_angle-1); i >= (end_angle/2); i--)
                        sigma_1 += hist[i] * pow((i - mean_1),2);

                        sigma_2 = omega_2 = mean_2 = 0;
                        for (i = ((end_angle/2)-1); i >= 0; i--) 
                        {
                                omega_2 += hist[i];
                                mean_2 += i * hist[i];
                                sum_sqr_2 += pow(i,2) * hist[i];
                                sum_2 += i * hist[i];
                        }
                        mean_2 /= omega_2;
                        for (i = ((end_angle/2)-1); i >= 0; i--)
                                sigma_2 += hist[i] * pow((i - mean_2),2);

                }
                else 
                {
                        tt = t + 1;
                        omega_1 -= hist[tt];
                        sum_sqr_1 -= pow(tt,2) * hist[tt];
                        sum_1 -= tt * hist[tt];

                        tt = ((t - (end_angle/2))+1);
                        omega_1 += hist[tt];
                        sum_sqr_1 += pow(tt,2) * hist[tt];
                        sum_1 += tt * hist[tt];

                        mean_1 = sum_1 / omega_1;

                        sigma_1 = sum_sqr_1 - pow(sum_1,2)/omega_1;

                        tt = ((t - (end_angle/2))+1);
                        omega_2 -= hist[tt];
                        sum_sqr_2 -= pow(tt,2) * hist[tt];
                        sum_2 -= tt * hist[tt];

                        tt = ((t - end_angle)+1);
                        omega_2 += hist[t-1];
                        sum_sqr_2 += pow(tt,2) * hist[t-1];
                        sum_2 += tt * hist[t-1];

                        mean_2 = sum_2 / omega_2;

                        sigma_2 = sum_sqr_2 - pow(sum_2,2)/omega_2;
                }
                if (omega_1 == 0) continue;
                if (omega_2 == 0) continue;
                sigma_W = sigma_1 + sigma_2;
                if (sigma_W < min_sigma_W) 
                {
                        min_sigma_W = sigma_W;
                        thr = t;
                }
        }
        return thr;
}

int ThresholdingForward(double *hist, int end_angle)
{
        int i,t,tt,thr;
        double min_sigma_W, sum_sqr_1, sum_sqr_2, sum_1, sum_2, sigma_1, omega_1, mean_1, sigma_2, omega_2, mean_2, sigma_W;
        thr = -1;
        min_sigma_W = 9e20;
        sum_sqr_1 = sum_1 = 0;
        sum_sqr_2 = sum_2 = 0;
        for (t = 0; t < (end_angle/2); t++) 
        {
                if (t == 0) 
                {
                        sigma_1 = omega_1 = mean_1 = 0;
                        for (i = 0; i <= ((end_angle/2)-1); i++) 
                        {
                                omega_1 += hist[i];
                                mean_1 += i * hist[i];
                                sum_sqr_1 += pow(i,2) * hist[i];
                                sum_1 += i * hist[i];
                        }
                        mean_1 /= omega_1;
                        for (i = 0; i <= ((end_angle/2)-1); i++)
                        sigma_1 += hist[i] * pow((i - mean_1),2);

                        sigma_2 = omega_2 = mean_2 = 0;
                        for (i = (end_angle/2); i <= (end_angle-1); i++) 
                        {
                                omega_2 += hist[i];
                                mean_2 += i * hist[i];
                                sum_sqr_2 += pow(i,2) * hist[i];
                                sum_2 += i * hist[i];
                        }
                        mean_2 /= omega_2;
                        for (i = (end_angle/2); i <= (end_angle-1); i++)
                                sigma_2 += hist[i] * pow((i - mean_2),2);

                }
                else 
                {
          
                        tt = t - 1;
                        omega_1 -= hist[tt];
                        sum_sqr_1 -= pow(tt,2) * hist[tt];
                        sum_1 -= tt * hist[tt];

                        tt = t + ((end_angle/2)-1);
                        omega_1 += hist[tt];
                        sum_sqr_1 += pow(tt,2) * hist[tt];
                        sum_1 += tt * hist[tt];

                        mean_1 = sum_1 / omega_1;

                        sigma_1 = sum_sqr_1 - pow(sum_1,2)/omega_1;

                        tt = t + ((end_angle/2)-1);
                        omega_2 -= hist[tt];
                        sum_sqr_2 -= pow(tt,2) * hist[tt];
                        sum_2 -= tt * hist[tt];

                        tt = t + (end_angle-1);
                        omega_2 += hist[t-1];
                        sum_sqr_2 += pow(tt,2) * hist[t-1];
                        sum_2 += tt * hist[t-1];

                        mean_2 = sum_2 / omega_2;

                        sigma_2 = sum_sqr_2 - pow(sum_2,2)/omega_2;
                }
                if (omega_1 == 0) continue;
                if (omega_2 == 0) continue;
                sigma_W = sigma_1 + sigma_2;
                if (sigma_W < min_sigma_W) 
                {
                        min_sigma_W = sigma_W;
                        thr = t;
                }
        }
        return thr;
}

double ComputeInitialCentroid(double *histogram_norm, int start_angle, int end_angle)
{
        int i;
        double sum1, sum2, centroid;
        sum1=0.0, sum2=0.0;
        for(i=start_angle;i<end_angle;i++)
        {
                sum1+=(histogram_norm[i]*sin(((double)i*PI)/180.0));
                sum2+=(histogram_norm[i]*cos(((double)i*PI)/180.0));
        }
        centroid=((atan2(sum1,sum2))*(180.0/PI));
        if(centroid<0)
        {
                centroid+=360.0;
        }
        return centroid;
}

double ComputeMeanResultantDirection(double *histogram_norm, int start_angle, int end_angle)
{
        int i;
        double sum, cos_sum, sin_sum, resultant_direction, mean_resultant_direction;
        sum=0.0, cos_sum=0.0, sin_sum=0.0;
        for(i=start_angle;i<end_angle;i++)
        {
                sum+=histogram_norm[i];
                cos_sum+=(histogram_norm[i]*cos(((double)i*PI)/180.0));
                sin_sum+=(histogram_norm[i]*sin(((double)i*PI)/180.0));
        }
        resultant_direction=sqrt(pow(cos_sum,2.0)+pow(sin_sum,2.0));
        mean_resultant_direction=((double)resultant_direction/sum);
        return mean_resultant_direction;
}

double ComputeInitialConcentration(double mean_resultant_direction)
{
        double kappa;
        
        if(mean_resultant_direction>=0 && mean_resultant_direction<0.53)
        {
                kappa=(double)((2 * mean_resultant_direction) + pow(mean_resultant_direction,3) + ((double)(5 * pow(mean_resultant_direction,5))/6));
        }
        else if(mean_resultant_direction>=0.53 && mean_resultant_direction<0.85)
        {
                kappa=(double)(((1.39 * mean_resultant_direction) + ((double)0.43/(1-mean_resultant_direction)))-0.4);
        }
        else
        {
                kappa=(double)((double)1/(((3 * mean_resultant_direction) + pow(mean_resultant_direction,3)) - (4 * pow(mean_resultant_direction,2))));
        }
        if(kappa<0)
        {
                kappa=fabs(kappa);
        }
        return kappa;
}

void GetAngularValues(double *AngleValue, int AngleLength)
{
        int i;
        for(i=0;i<AngleLength;i++)
        {
                AngleValue[i]=(((double)i*PI)/180.0);
        }
}

void CopyOriginalArray(double *array1, double *array2, int Numofelements)
{
	int i;
	for(i=0;i<Numofelements;i++)
	{
                array1[i]=array2[i];
	}
}

void CopyOriginalMatrix(double **matrix1, double **matrix2, int row, int col)
{
        int i,j;
        for(i=0;i<row;i++)
        {
                for(j=0;j<col;j++)
                {
                        matrix1[i][j]=matrix2[i][j];
                }
        }
}

double GetScalingFactorValue(double lower, double upper, double interval, double concentration)
{
        int i,n;
        double sum_odd,sum_even,ans;
        double *x, *y;
        x=(double*)calloc(5000, sizeof(double));
        y=(double*)calloc(5000, sizeof(double));
        n=(int)((upper-lower)/interval);
        if(n%2==1)
        {
                n=n+1;
        }
        interval=((double)(upper-lower)/n);
        for(i=0; i<=n; i++)
        {
                x[i]=lower+(i*interval);
                y[i]=exp(concentration * cos(x[i]));
        }
        sum_odd=0.0;
        sum_even=0.0;
        for(i=1; i<n; i++)
        {
                if(i%2==1)
                {
                        sum_odd=sum_odd+y[i];
                }
                else
                {
                        sum_even=sum_even+y[i];
                }
 
        }
        ans=((double)interval/3)*((y[0]+y[n])+(4*sum_odd)+(2*sum_even));
        free(x);
        free(y);
        return(ans);
}

double getThresholdValue(double **MembershipValue, int num_of_clusters, int AngleLength, int *Modified_hue_histogram, double *SW_hue_histogram, int AngleFreqSum, double SWFreqSum)
{
	int i,j,index;
	double max1,max2,threshold,sum;
	sum=0.00;
	for(j=0;j<AngleLength;j++)
	{
		max1=0.00;
		for(i=0;i<num_of_clusters;i++)
                {
                        if(MembershipValue[i][j]>=max1)
                        {
                                max1=MembershipValue[i][j];
                                index=i;
                        }
                }
                max2=0.00;
                for(i=0;i<num_of_clusters;i++)
                {
                        if(MembershipValue[i][j]>=max2 && i!=index)
                        {
                                max2=MembershipValue[i][j];
                        }
                }
                sum+=(fabs(max1-max2)*((double)Modified_hue_histogram[j]+SW_hue_histogram[j]));
	}
	threshold=((double)sum/((double)AngleFreqSum+SWFreqSum));
	return threshold;
}

void NormalizeMembershipValue(double **MembershipValue, int num_of_clusters, int AngleLength)
{
        int i,j;
        double sum;
        for(j=0;j<AngleLength;j++)
        {
                sum=0.00;
                for(i=0;i<num_of_clusters;i++)
                {
                        sum+=MembershipValue[i][j];
                }
                for(i=0;i<num_of_clusters;i++)
                {
                        MembershipValue[i][j]=((double)MembershipValue[i][j]/sum);
                }
        }
}

void UpdateMembershipValue(double *Concentration,double *AngleValue,double *ClusterCenter,int num_of_clusters,double **MembershipValue,int AngleLength,double fuzzifier, int *Modified_hue_histogram, double *SW_hue_histogram, int AngleFreqSum, double SWFreqSum)
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
                        MembershipNumerator=exp(Concentration[i]*cos(AngleValue[j]-ClusterCenter[i]));
                        MembershipValue[i][j]=pow(((double)MembershipNumerator/MembershipDenominator),factor);
                }
        }
        NormalizeMembershipValue(MembershipValue, num_of_clusters, AngleLength);
	
	for(j=0;j<AngleLength;j++)
	{
		count=0;
		for(i=0;i<num_of_clusters;i++)
		{
			if(MembershipValue[i][j]==1.0)
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
				if(MembershipValue[i][j]==1.0)
				{
					MembershipValue[i][j]=1.0;
				}
				else
				{
					MembershipValue[i][j]=0.0;
				}
			}
		}
		
		else if(count==0)
		{
			sum=0.0;
			for(i=0;i<num_of_clusters;i++)
			{
				sum+=MembershipValue[i][j];
			}
			for(i=0;i<num_of_clusters;i++)
			{
				MembershipValue[i][j]=(double)MembershipValue[i][j]/sum;
			}
		}
	}
}

void find_roughsets(double **MembershipValue,int AngleLength,int num_of_clusters,int *lower,int **boundary,double threshold)
{
        int i,j,k;
        int tempclass;
        double temp;
        int *class;
        double *array;
        class=(int*)calloc(num_of_clusters, sizeof(int));
        array=(double*)calloc(num_of_clusters, sizeof(double));
        for(i=0;i<AngleLength;i++)
        {
                lower[i]=-1;
                for(j=0;j<num_of_clusters;j++)
                {
                        boundary[i][j]=-1;
                }
        }
        for(i=0;i<AngleLength;i++)
        {
                for(j=0;j<num_of_clusters;j++)
                {
                        array[j]=MembershipValue[j][i];
                        class[j]=j;
                }
                for(j=0;j<num_of_clusters;j++)
                {
                        for(k=j+1;k<num_of_clusters;k++)
                        {
                                if(array[j]<array[k])
                                {
                                        temp=array[j];
                                        array[j]=array[k];
                                        array[k]=temp;
                                        tempclass=class[j];
                                        class[j]=class[k];
                                        class[k]=tempclass;
                                }
                        }
                }
                if((array[0]-array[1])>threshold)
                {
                        lower[i]=class[0];
                }
                else
                {
                        boundary[i][class[0]]=1;
                        for(j=1;j<num_of_clusters;j++)
                        {
                                if((array[0]-array[j])<=threshold)
                                {
                                        boundary[i][class[j]]=1;
                                }
                        }
                }

        }
        free(class);
        free(array);
}

double GetActualAngle(double numerator,double denominator)
{
        double theta,angle;
        theta=atan2(numerator,denominator);
        if(theta<0)
        {
                angle=theta+(2*PI);
        }
        else
        {
                angle=theta;
        }
        return angle;
}

void estimate_centroid(double **MembershipValue,double *AngleValue,double *ClusterCenter,int num_of_clusters,int AngleLength,double fuzzifier,int *lower,int **boundary,int *Modified_hue_histogram,double *SW_hue_histogram,double weightage)
{
        int i,j;
        double sin_lower,cos_lower,sin_boundary,cos_boundary,centroid_numerator,centroid_denominator;
        for(i=0;i<num_of_clusters;i++)
        {
                sin_lower=0.0, cos_lower=0.0, sin_boundary=0.0, cos_boundary=0.0;
                for(j=0;j<AngleLength;j++)
                {
                        if(lower[j]==i)
                        {
                                sin_lower+=(sin(AngleValue[j])*((double)Modified_hue_histogram[j]+SW_hue_histogram[j]));
                                cos_lower+=(cos(AngleValue[j])*((double)Modified_hue_histogram[j]+SW_hue_histogram[j]));
                        }
                        if(boundary[j][i]==1)
                        {
                                sin_boundary+=((pow(MembershipValue[i][j],fuzzifier)*sin(AngleValue[j]))*((double)Modified_hue_histogram[j]+SW_hue_histogram[j]));
                                cos_boundary+=((pow(MembershipValue[i][j],fuzzifier)*cos(AngleValue[j]))*((double)Modified_hue_histogram[j]+SW_hue_histogram[j]));
                        }
                }
                centroid_numerator=((weightage*sin_lower)+((1.0-weightage)*sin_boundary));
		centroid_denominator=((weightage*cos_lower)+((1.0-weightage)*cos_boundary));
		ClusterCenter[i]=GetActualAngle(centroid_numerator, centroid_denominator);
        }
}

double extract_kappa_value(double bessel_ratio)
{
        double kappa;
        if(bessel_ratio>=0 && bessel_ratio<0.53)
        {
                kappa=(double)((2 * bessel_ratio) + pow(bessel_ratio,3) + ((double)(5 * pow(bessel_ratio,5))/6));
        }
        else if(bessel_ratio>=0.53 && bessel_ratio<0.85)
        {
                kappa=(double)(((1.39 * bessel_ratio) + ((double)0.43/(1-bessel_ratio)))-0.4);
        }
        else
        {
                kappa=(double)((double)1/(((3 * bessel_ratio) + pow(bessel_ratio,3)) - (4 * pow(bessel_ratio,2))));
        }
        
        if(kappa<0)
        {
                kappa=fabs(kappa);
        }
        return kappa;
}

void estimate_concentration_kappa(double **MembershipValue,double *AngleValue,double *ClusterCenter,double *Concentration,int num_of_clusters,int AngleLength,double fuzzifier,int *lower,int **boundary,int *Modified_hue_histogram,double *SW_hue_histogram,int weightage)
{
	int i,j;
        double cos_lower,member_boundary,cos_boundary,kappa_numerator,kappa_denominator,bessel_ratio,cardinality_lower;
        for(i=0;i<num_of_clusters;i++)
        {
		cardinality_lower=0.0;
		cos_lower=0.0, member_boundary=0.0, cos_boundary=0.0;
                for(j=0;j<AngleLength;j++)
                {
                        if(lower[j]==i)
                        {
                                cardinality_lower+=((double)Modified_hue_histogram[j]+SW_hue_histogram[j]);
				cos_lower+=((cos(AngleValue[j]-ClusterCenter[i]))*((double)Modified_hue_histogram[j]+SW_hue_histogram[j]));
                                
                        }
                        if(boundary[j][i]==1)
                        {
				member_boundary+=(pow(MembershipValue[i][j],fuzzifier)*((double)Modified_hue_histogram[j]+SW_hue_histogram[j]));
				cos_boundary+=((pow(MembershipValue[i][j],fuzzifier))*(cos(AngleValue[j]-ClusterCenter[i]))*((double)Modified_hue_histogram[j]+SW_hue_histogram[j]));
                        }
                }
                kappa_numerator=((weightage*cos_lower)+((1.0-weightage)*cos_boundary));
		kappa_denominator=((weightage*cardinality_lower)+((1.0-weightage)*member_boundary));
		bessel_ratio=((double)kappa_numerator/kappa_denominator);
		Concentration[i]=extract_kappa_value(bessel_ratio);
        }
}

void getHSIInformation(FILE *fp, int num_of_clusters, int AngleLength, double **HSI)
{
        int i,j;
        for(i=0;i<AngleLength;i++)
        {
                for(j=0;j<num_of_clusters;j++)
                {
                        fscanf(fp,"%lf",&HSI[i][j]);
                }
        }
        fclose(fp);
}

void findClusterCenterSaturation(double **HSI, double *ClusterCenter, int num_of_clusters, int TotalPixels, double *saturation_center, double *intensity_center)
{
        int i,j;
        double saturation_sum, intensity_weighted_sum, saturation_weighted_sum;
        for(i=0;i<num_of_clusters;i++)
        {
                saturation_sum=0.00, intensity_weighted_sum=0.00, saturation_weighted_sum=0.00;
                for(j=0;j<TotalPixels;j++)
                {
                        if(fabs(ClusterCenter[i]-HSI[j][0])<1.00)
                        {
                                saturation_sum+=HSI[j][1];
                                saturation_weighted_sum+=pow(HSI[j][1],2);
                                intensity_weighted_sum+=(HSI[j][1] * HSI[j][2]);
                        }
                }
                if(saturation_sum==0)
                {
                        saturation_center[i]=0.00;
                        intensity_center[i]=0.00;
                }
                else
                {
                        saturation_center[i]=(double)(saturation_weighted_sum/saturation_sum);
                        intensity_center[i]=(double)(intensity_weighted_sum/saturation_sum);
                }
        }
}

struct RGB RG_Sector_Conversion(double hue, double saturation, double intensity)
{
        struct RGB color;
        double red,green,blue;
        blue=(intensity*(1.00-saturation));
        red=(intensity*(1.00+(saturation*cos(hue*PI/180))/cos((60-hue)*PI/180)));
        green=((3.00*intensity)-(red+blue));
        color.red=red;
        color.green=green;
        color.blue=blue;
        return color;
}

struct RGB GB_Sector_Conversion(double hue, double saturation, double intensity)
{
        struct RGB color;
        double red,green,blue;
        hue=(hue-120);
        red=(intensity*(1.00-saturation));
        green=(intensity*(1.00+(saturation*cos(hue*PI/180))/cos((60-hue)*PI/180)));
        blue=((3.00*intensity)-(red+green));
        color.red=red;
        color.green=green;
        color.blue=blue;
        return color;
}

struct RGB BR_Sector_Conversion(double hue, double saturation, double intensity)
{
        struct RGB color;
        double red,green,blue;
        hue=(hue-240);
        green=(intensity*(1.00-saturation));
        blue=(intensity*(1.00+(saturation*cos(hue*PI/180))/cos((60-hue)*PI/180)));
        red=((3.00*intensity)-(green+blue));
        color.red=red;
        color.green=green;
        color.blue=blue;
        return color;
}

void HSI_to_RGB_conversion(double *hue_center, double *saturation_center, double *intensity_center, double *normalized_red, double *normalized_green, double *normalized_blue, int num_of_clusters)
{
        int i;
        struct RGB color_component;
        for(i=0;i<num_of_clusters;i++)
        {
                if(hue_center[i]>0 && hue_center[i]<=120)
                {
                        color_component=RG_Sector_Conversion(hue_center[i], saturation_center[i], intensity_center[i]);
                        normalized_red[i]=color_component.red;
                        normalized_green[i]=color_component.green;
                        normalized_blue[i]=color_component.blue;
                }
                else if(hue_center[i]>120 && hue_center[i]<=240)
                {
                        color_component=GB_Sector_Conversion(hue_center[i], saturation_center[i], intensity_center[i]);
                        normalized_red[i]=color_component.red;
                        normalized_green[i]=color_component.green;
                        normalized_blue[i]=color_component.blue;
                }
                else
                {
                        color_component=BR_Sector_Conversion(hue_center[i], saturation_center[i], intensity_center[i]);
                        normalized_red[i]=color_component.red;
                        normalized_green[i]=color_component.green;
                        normalized_blue[i]=color_component.blue;
                }
        }
        
}

void convert_to_log_domain(double *normalized_red, double *normalized_green, double *normalized_blue, int num_of_clusters)
{
        int i;
        double sum_stain,sqrt_sum_stain;
        for(i=0;i<num_of_clusters;i++)
        {
                if(normalized_red[i]!=0)
                {
			normalized_red[i]=log((double)1/normalized_red[i]);
                }
                
                if(normalized_green[i]!=0)
                {
			normalized_green[i]=log((double)1/normalized_green[i]);
                }
                
                if(normalized_blue[i]!=0)
                {
			normalized_blue[i]=log((double)1/normalized_blue[i]);
                }
        }
        for(i=0;i<num_of_clusters;i++)
        {   
                sum_stain=pow(normalized_red[i],2.0)+pow(normalized_green[i],2.0)+pow(normalized_blue[i],2.0);
                sqrt_sum_stain=fabs(sqrt(sum_stain));
                normalized_red[i]=(double)normalized_red[i]/sqrt_sum_stain;
                normalized_green[i]=(double)normalized_green[i]/sqrt_sum_stain;
                normalized_blue[i]=(double)normalized_blue[i]/sqrt_sum_stain;
        }
}

double zero_check_normalization(int val)
{
        double value;
        if(val==0)
            value=((double)255/(val+1));
        else
            value=(double)255/val;
        
        return (log(value));
}

void read_ppm_write_pgmMatrix_Mapping(FILE *fp, int **image_red, int **image_green, int **image_blue, double **normalized_image_red, double **normalized_image_green, double **normalized_image_blue, int col, int row)
{
	int i,j,red,green,blue;
	char ch;
	fscanf(fp,"%c",&ch);
	for(i = 0; i < row; i++)
	{
		for(j=0;j< col;j++)
		{
			red=fgetc(fp);
			image_red[i][j]=red;
			green=fgetc(fp);
			image_green[i][j]=green;
			blue=fgetc(fp);
			image_blue[i][j]=blue;
		}
	}
        for(i=0;i<row;i++)
        {
                for(j=0;j<col;j++)
                {
                    
                        normalized_image_red[i][j]=zero_check_normalization(image_red[i][j]);
                        normalized_image_green[i][j]=zero_check_normalization(image_green[i][j]);
                        normalized_image_blue[i][j]=zero_check_normalization(image_blue[i][j]);
                }
        }	
	fclose(fp);
}

void GetResizedImageMatrix(double **normalized_image_red, double **normalized_image_green, double **normalized_image_blue, int row, int col, double **ResizedImage)
{
        int i,j,k;
        int size=(row*col);
        k=0;
        while(k<size)
        {
                for(i=0;i<row;i++)
                {
                        for(j=0;j<col;j++)
                        {
                                ResizedImage[k][0]=normalized_image_red[i][j];
                                ResizedImage[k][1]=normalized_image_green[i][j];
                                ResizedImage[k][2]=normalized_image_blue[i][j];
                                k++;
                        }
                }
        }
}

double** GetMinor(double **matrix,int order,int row_index,int col_index)
{
        int p,q,flag,i,j;
        p=0;
        q=0;
        double **minor_of_element;
        minor_of_element=double_matrix_allocation(order-1,order-1);
	for(i=0;i<order;i++)
	{
		for(j=0;j<order;j++)
		{
			if((i==row_index && col_index!=(order-1)) || (j==col_index && col_index!=(order-1))  || (i==row_index && j==col_index && col_index!=(order-1)))
			{flag=0;}
			else if((i==row_index && col_index==(order-1)) || (j==col_index && col_index==(order-1))  || (i==row_index && j==col_index && col_index==(order-1)))
			{}
			else
			{
				minor_of_element[p][q]=matrix[i][j];
				if(q<(order-2))
				{
					q++;
				}
				flag=1;	
			}
		}
		q=0;
		if((p<(order-2)) && flag==1)
		{
			p++;
		}
	}
	return minor_of_element;
}

double Determinant_value(double **matrix,int order)
{
        int i;
        double sum=0.0;

        if(order==1)
        {
                return(matrix[0][0]);
        }
        else
        {
                for(i=0;i<order;i++)
                {
                        sum+=(matrix[0][i]*pow((-1),(0+i))*Determinant_value(GetMinor(matrix,order,0,i),order-1));
                }
        }
        return sum;
}

void transpose(double **matrix,double **transposed_matrix, int row, int col)
{
        int i,j;
        for(i=0;i<row;i++)
        {
                for(j=0;j<col;j++)
                {
                        transposed_matrix[j][i]=matrix[i][j];
                }
        }
}

void GetInverse(double **matrix,double **inverse_matrix,int order)
{
        int i,j;
        double value=Determinant_value(matrix,order);
        if(value==0)
        {
                printf("\nThe given matrix is a Singular Matrix  - Inverse of the matrix does not exist\n\n");
        }
        else
        {
                double **cofactor_matrix,**cofactor_matrix_transpose;
                cofactor_matrix=double_matrix_allocation(order,order);
                cofactor_matrix_transpose=double_matrix_allocation(order,order);

	        for(i=0;i<order;i++)
	        {
	                for(j=0;j<order;j++)
	                {
	                        cofactor_matrix[i][j]=(pow((-1),(i+j))*Determinant_value(GetMinor(matrix,order,i,j),order-1));
	                }
	        }  
	        transpose(cofactor_matrix,cofactor_matrix_transpose,order,order);
	        for(i=0;i<order;i++)
	        {
	                for(j=0;j<order;j++)
	                {
	                        inverse_matrix[i][j]=((double)cofactor_matrix_transpose[i][j]/value);
	                }
	        }  
	        double_matrix_deallocation(cofactor_matrix,order);
                double_matrix_deallocation(cofactor_matrix_transpose,order);
	 } 
}

void MatrixMultiplication(double **firstMatrix, double **secondMatrix, int row, int order, int col, double **multipliedmatrix)
{
	int i,j,k;
	double sum=0.00;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			for(k=0;k<order;k++)
			{
				sum+=(firstMatrix[i][k] * secondMatrix[k][j]);
			}
			multipliedmatrix[i][j]=FindMaximum(sum,1e-10);
			sum=0.00;
		}
	}
}

void UpdateFirstMatrix(double **firstMatrix, double **secondMatrix, double **ActualMatrix, int row, int col, int order)
{
	int i,j;
	double **intermediateMatrix1, **intermediateMatrix2, **intermediateMatrix3, **secondMatrixTranspose;
	intermediateMatrix1=double_matrix_allocation(row, order);
	intermediateMatrix2=double_matrix_allocation(row, col);
	intermediateMatrix3=double_matrix_allocation(row, order);
	secondMatrixTranspose=double_matrix_allocation(col, order);
	transpose(secondMatrix, secondMatrixTranspose, order, col);
	MatrixMultiplication(ActualMatrix, secondMatrixTranspose, row, col, order, intermediateMatrix1);
	MatrixMultiplication(firstMatrix, secondMatrix, row, order, col, intermediateMatrix2);
	MatrixMultiplication(intermediateMatrix2, secondMatrixTranspose, row, col, order, intermediateMatrix3);
	for(i=0;i<row;i++)
	{
		for(j=0;j<order;j++)
		{
			firstMatrix[i][j]=(firstMatrix[i][j] * (double)(intermediateMatrix1[i][j]/intermediateMatrix3[i][j]));
		}
	}
	
	double_matrix_deallocation(intermediateMatrix1, row);
	double_matrix_deallocation(intermediateMatrix2, row);
	double_matrix_deallocation(intermediateMatrix3, row);
	double_matrix_deallocation(secondMatrixTranspose, col);
}

void UpdateSecondMatrix(double **firstMatrix, double **secondMatrix, double **ActualMatrix, int row, int col, int order)
{
	int i,j;
	double **intermediateMatrix1, **intermediateMatrix2, **intermediateMatrix3, **firstMatrixTranspose;
	intermediateMatrix1=double_matrix_allocation(order, col);
	intermediateMatrix2=double_matrix_allocation(order, order);
	intermediateMatrix3=double_matrix_allocation(order, col);
	firstMatrixTranspose=double_matrix_allocation(order, row);
	transpose(firstMatrix, firstMatrixTranspose, row, order);
	MatrixMultiplication(firstMatrixTranspose, ActualMatrix, order, row, col, intermediateMatrix1);
	MatrixMultiplication(firstMatrixTranspose, firstMatrix, order, row, order, intermediateMatrix2);
	MatrixMultiplication(intermediateMatrix2, secondMatrix, order, order, col, intermediateMatrix3);
	for(i=0;i<order;i++)
	{
		for(j=0;j<col;j++)
		{
			secondMatrix[i][j]=(secondMatrix[i][j] * (double)(intermediateMatrix1[i][j]/intermediateMatrix3[i][j]));
		}
	}
	double_matrix_deallocation(intermediateMatrix1, order);
	double_matrix_deallocation(intermediateMatrix2, order);
	double_matrix_deallocation(intermediateMatrix3, order);
	double_matrix_deallocation(firstMatrixTranspose, order);
}

int checkifLargeDifferenceArray(double *array1, double *array2, int Numofelements, double epsilon)
{
	int i,flag;
        flag=0;
        for(i=0;i<Numofelements;i++)
        {
                if(fabs(array1[i]-array2[i])>=epsilon)
                {
                        flag=1;
                }
        }
	return flag;
}

int checkifLargeDifferenceMatrix(double **matrix1, double **matrix2, int row, int col, double eps)
{
	int i,j, flag=0;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			if(fabs(matrix1[i][j]-matrix2[i][j])>eps)
			{
				flag=1;
			}
		}
	}
	return flag;
}

void MapConcentrationDensity(double **centermatrix,double **stain_density,int row,int col,int num_of_clusters,int num_of_channels,int maxval,double **MappedMatrix)
{
        int i,j;
	MatrixMultiplication(stain_density, centermatrix, (row*col), num_of_clusters, num_of_channels, MappedMatrix);
        for(i=0;i<(row*col);i++)
        {
                for(j=0;j<num_of_channels;j++)
                {
                        MappedMatrix[i][j]=maxval*exp(-MappedMatrix[i][j]);
                }
        }
}

void ReformImageChannels(double **MappedMatrix,int **image_red,int **image_green,int **image_blue,int row,int col)
{
        int i,j,k;
        int size=(row*col);
        k=0;
        while(k<size)
        {
                for(i=0;i<row;i++)
                {
                        for(j=0;j<col;j++)
                        {
                                image_red[i][j]=(int)(MappedMatrix[k][0]+0.5);
                                image_green[i][j]=(int)(MappedMatrix[k][1]+0.5);
                                image_blue[i][j]=(int)(MappedMatrix[k][2]+0.5);
                                k++;
                        }
                }
        }
}

void read_pgmMatrix_write_ppm(int **image_red, int **image_green, int **image_blue, FILE *ft, int col, int row, int maxval)
{
	int i,j,red,green,blue;
	fprintf(ft,"P6\n%d %d\n%d\n",col,row,maxval);
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			red=image_red[i][j];
			fputc(red,ft);
			green=image_green[i][j];
			fputc(green,ft);
			blue=image_blue[i][j];
			fputc(blue,ft);
		}
	}
}

void Test_SW_hue_histogram(FILE *f_imagenames, char *filepath, int maxcount, int rectsize, double alpha)
{
        int i, j, row, col, maxval, count=0;
        double gamma=2.22;
	FILE *f_source, *f_hueHist, *f_HSI, *f_Modified_hist;
	double **normalized_image_red, **normalized_image_green, **normalized_image_blue, **theta, **saturation, **intensity, **hue_actual;
	int **image_red, **image_green, **image_blue, **hue, **Modified_hue;
    
        while(count<maxcount)
        {
                int i, row, col, maxval,AngleLength;
                char ch1, ch2, ch3;
                char *imagename=(char*)malloc(100 * sizeof(char));
                char *imagepath=(char*)malloc(200 * sizeof(char));
                char *sourcefile=(char*)malloc(100 * sizeof(char));
		char *data_path=(char*)calloc(200, sizeof(char));
		double *SW_hue_histogram=(double*)calloc(360, sizeof(double));
		int *Modified_hue_histogram=(int*)calloc(360, sizeof(double));
		
                fscanf(f_imagenames,"%s",sourcefile);
		strcpy(data_path, filepath);
		strcat(data_path, sourcefile);
                f_source=fopen(data_path,"r");
                fscanf(f_source,"%c%c%c",&ch1,&ch2,&ch3);
                skip_comments(f_source);
                fscanf(f_source,"%d %d",&col,&row);
                fscanf(f_source,"%d",&maxval);
                cut_string(sourcefile, (strlen(sourcefile)-4), imagename);
                image_red=integer_matrix_allocation(row, col);
                image_green=integer_matrix_allocation(row, col);
                image_blue=integer_matrix_allocation(row, col);	
                normalized_image_red=double_matrix_allocation(row, col);
                normalized_image_green=double_matrix_allocation(row, col);
                normalized_image_blue=double_matrix_allocation(row, col);
                read_ppm_write_pgmMatrix(f_source, image_red, image_green, image_blue, normalized_image_red, normalized_image_green, normalized_image_blue, col, row);
                nonlinearRGB_to_linearRGB(normalized_image_red, normalized_image_green, normalized_image_blue, col, row, gamma);
                hue_actual=double_matrix_allocation(row, col);
		Modified_hue=integer_matrix_allocation(row, col);
                theta=double_matrix_allocation(row, col);
                saturation=double_matrix_allocation(row, col);
                intensity=double_matrix_allocation(row, col);
                RGB_to_HSI_conversion(normalized_image_red, normalized_image_green, normalized_image_blue, theta, saturation, intensity, col, row);
                hue=integer_matrix_allocation(row, col);
                compute_SW_hue_histograms(theta, saturation, hue, hue_actual, row, col, SW_hue_histogram);
		Compute_Neighborhood_Average(hue_actual, row, col, rectsize, alpha, Modified_hue, Modified_hue_histogram);
                mkdir(imagename,0777);
                strcpy(imagepath,imagename);
                strcat(imagepath, "/HSI_values.txt");
                f_HSI=fopen(imagepath, "w");
                fprintf(f_HSI,"%d\n",(row*col));
                for(i=0;i<row;i++)
                {
                        for(j=0;j<col;j++)
                        {
                                fprintf(f_HSI,"%lf %lf %lf\n",hue_actual[i][j],saturation[i][j],intensity[i][j]);
                        }
                }
                fclose(f_HSI);
		strcpy(imagepath,imagename);
                strcat(imagepath, "/SW_hue_hist.txt");
                f_hueHist=fopen(imagepath, "w");
		for(i=0;i<360;i++)
		{
			fprintf(f_hueHist,"%lf\n",SW_hue_histogram[i]);
		}
		fclose(f_hueHist);
		strcpy(imagepath,imagename);
                strcat(imagepath, "/Modified_hue_hist.txt");
                f_Modified_hist=fopen(imagepath, "w");
		for(i=0;i<360;i++)
		{
			fprintf(f_Modified_hist,"%d\n",Modified_hue_histogram[i]);
		}
		fclose(f_Modified_hist);
                free(imagename);
                free(imagepath);
                free(sourcefile);
		free(data_path);
                double_matrix_deallocation(normalized_image_red,row);
                double_matrix_deallocation(normalized_image_green,row);
                double_matrix_deallocation(normalized_image_blue,row);
                double_matrix_deallocation(hue_actual,row);
                double_matrix_deallocation(theta,row);
                double_matrix_deallocation(saturation,row);
                double_matrix_deallocation(intensity,row);
                integer_matrix_deallocation(image_red,row);
                integer_matrix_deallocation(image_green,row);
                integer_matrix_deallocation(image_blue,row);
                integer_matrix_deallocation(hue,row);
		integer_matrix_deallocation(Modified_hue, row);
                free(SW_hue_histogram);
		free(Modified_hue_histogram);
                count++;
                printf("\nSW Histogram computation on Image %d is completed\n", count);
        }
}

void TestInitialization(FILE *f_imagenames,int maxcount)
{
        int num_of_samples, thres_E_H, thres_H_Ach, thres_Start_Ach, threshold1, threshold2, count;
	double center_Ach, center_H, center_E, kappa_Ach, kappa_H, kappa_E, MeanLength_Ach, MeanLength_H, MeanLength_E;
	char ch1, ch2, ch3;
	double *histogram_smooth, *histogram_norm, *histogram_intermediate1, *histogram_intermediate2, *histogram_intermediate3, *SW_hue_histogram;
	FILE *f_source, *f_smooth_hist, *f_conc, *f_center, *f_hueHist;
	count=0;
        while(count<maxcount)
        {
                int i, row, col, maxval,AngleLength;
                char *hue_file;
		FILE *f_center, *f_conc;
                char *imagename=(char*)malloc(50 * sizeof(char));
                char *imagepath=(char*)malloc(200 * sizeof(char));
                char *sourcefile=(char*)malloc(75 * sizeof(char));
		SW_hue_histogram=(double*)calloc(360, sizeof(double));
		histogram_smooth=(double*)calloc(360, sizeof(double));
		histogram_norm=(double*)calloc(360, sizeof(double));
                fscanf(f_imagenames,"%s",sourcefile);
                f_source=fopen(sourcefile,"r");
                cut_string(sourcefile, (strlen(sourcefile)-4), imagename);
		strcpy(imagepath,imagename);
                strcat(imagepath, "/SW_hue_hist.txt");
                f_hueHist=fopen(imagepath, "r");
		for(i=0;i<360;i++)
		{
			fscanf(f_hueHist,"%lf",&SW_hue_histogram[i]);
		}
		fclose(f_hueHist);
		
		strcpy(imagepath,imagename);
                strcat(imagepath, "/Smooth_Histogram.txt");
                f_smooth_hist=fopen(imagepath, "w");
		GetSmoothHistogram(f_smooth_hist, SW_hue_histogram, histogram_smooth);
		ComputeNormalizedHistogram(histogram_smooth, histogram_norm, 0, 360);
		thres_Start_Ach=ThresholdingForward(histogram_norm, 360);
		histogram_intermediate1=(double*)calloc((360-(thres_Start_Ach+1)),sizeof(double));
		ComputeNormalizedHistogram(histogram_smooth, histogram_intermediate1, (thres_Start_Ach+1), 360);
		thres_E_H=thres_Start_Ach+ThresholdingBackward(histogram_intermediate1, (360-(thres_Start_Ach+1)));
		histogram_intermediate2=(double*)calloc((thres_E_H-(thres_Start_Ach+1)),sizeof(double));
		ComputeNormalizedHistogram(histogram_smooth, histogram_intermediate2, (thres_Start_Ach+1), thres_E_H);
		thres_H_Ach=thres_Start_Ach+ThresholdingBackward(histogram_intermediate2, (thres_E_H-(thres_Start_Ach+1)));
		histogram_intermediate3=(double*)calloc((thres_H_Ach-(thres_Start_Ach+1)),sizeof(double));
		ComputeNormalizedHistogram(histogram_smooth, histogram_intermediate3, (thres_Start_Ach+1), thres_H_Ach);
		threshold1=thres_Start_Ach+ThresholdingBackward(histogram_intermediate3, (thres_H_Ach-(thres_Start_Ach+1)));
		center_Ach=ComputeInitialCentroid(histogram_norm, 0, thres_Start_Ach);
		center_H=ComputeInitialCentroid(histogram_norm, (thres_Start_Ach+1), thres_E_H);
		center_E=ComputeInitialCentroid(histogram_norm, (thres_E_H+1), 360);
		MeanLength_Ach=ComputeMeanResultantDirection(histogram_norm, 0, thres_Start_Ach);
		MeanLength_H=ComputeMeanResultantDirection(histogram_norm, (thres_Start_Ach+1), thres_E_H);
		MeanLength_E=ComputeMeanResultantDirection(histogram_norm, (thres_E_H+1), 360);
		kappa_Ach=ComputeInitialConcentration(MeanLength_Ach);
		kappa_H=ComputeInitialConcentration(MeanLength_H);
		kappa_E=ComputeInitialConcentration(MeanLength_E);
		strcpy(imagepath,imagename);
                strcat(imagepath, "/Initial_Centroid.txt");
		f_center=fopen(imagepath,"w");
		fprintf(f_center,"%lf\t%lf\t%lf\n",center_Ach,center_H,center_E);
		fclose(f_center);
		strcpy(imagepath,imagename);
                strcat(imagepath, "/Initial_Concentration.txt");
		f_conc=fopen(imagepath,"w");
		fprintf(f_conc,"%lf\t%lf\t%lf\n",kappa_Ach,kappa_H,kappa_E);
		fclose(f_conc);
		free(imagename);
                free(imagepath);
                free(sourcefile);
		free(SW_hue_histogram);
		free(histogram_smooth);
		free(histogram_norm);
		free(histogram_intermediate1);
		free(histogram_intermediate2);
		free(histogram_intermediate3);
                count++;
                printf("\nOpeartion on Image %d has been accomplished\n", count);
        }
        fclose(f_imagenames);
}

void TestCircularClustering(FILE *f_imagenames, int maxcount, int num_of_clusters, double *weightage, int parameter_length, double fuzzifier)
{
        int i,j,flag,AngleLength, objective_value, objective_new, count, itercount, maxiter, AngleFreqSum;
        double threshold, epsilon, SWFreqSum;
	int *Modified_hue_histogram;
	double *SW_hue_histogram;
        double **MembershipValue;
	FILE *f_threshold;
	count=0, itercount=0, maxiter=100, AngleLength=360, epsilon=0.99;
	
	f_threshold=fopen("Threshold_Value.txt","w");
        while(count<maxcount)
        {
                count++;
                printf("\n\n********************Opeartion on Image %d has been started********************\n", count);
                char *imagename=(char*)calloc(200, sizeof(char));
                char *imagepath=(char*)calloc(200, sizeof(char));
                char *sourcefile=(char*)calloc(200, sizeof(char));
                char *path=(char*)calloc(200, sizeof(char));
		double *AngleValue=(double*)calloc(360, sizeof(double));
                FILE *f_hueHist,*f_iter,*f_centers, *f_kappa;
                fscanf(f_imagenames,"%s",sourcefile);
                cut_string(sourcefile, (strlen(sourcefile)-4), imagename);
		GetAngularValues(AngleValue, AngleLength);
                for(j=0;j<parameter_length;j++)
                {
                        double *ClusterCenter=(double*)calloc(num_of_clusters, sizeof(double));
                        double *Concentration=(double*)calloc(num_of_clusters, sizeof(double));
                        double *ClusterCenterOld=(double*)calloc(num_of_clusters, sizeof(double));
			Modified_hue_histogram=(int*)calloc(AngleLength, sizeof(int));
			SW_hue_histogram=(double*)calloc(AngleLength, sizeof(double));
                        int *lower=(int*)calloc(AngleLength, sizeof(int));
                        int **boundary=integer_matrix_allocation(AngleLength,num_of_clusters);
                        MembershipValue=double_matrix_allocation(num_of_clusters,AngleLength);
			FILE *f_init_center, *f_threshold;
			strcpy(imagepath,imagename);
			strcat(imagepath, "/SW_hue_hist.txt");
			f_hueHist=fopen(imagepath, "r");
			for(i=0;i<360;i++)
			{
				fscanf(f_hueHist,"%lf",&SW_hue_histogram[i]);
			}
			fclose(f_hueHist);
			strcpy(imagepath,imagename);
			strcat(imagepath, "/Modified_hue_hist.txt");
			f_hueHist=fopen(imagepath, "r");
			for(i=0;i<360;i++)
			{
				fscanf(f_hueHist,"%d",&Modified_hue_histogram[i]);
			}
			fclose(f_hueHist);
			strcpy(imagepath,imagename);
			strcat(imagepath, "/Initial_Centroid.txt");
			f_init_center=fopen(imagepath,"r");
			fscanf(f_init_center,"%lf  %lf  %lf",&ClusterCenter[0],&ClusterCenter[1],&ClusterCenter[2]);
			fclose(f_init_center);
			FILE *f_init_conc;
			strcpy(imagepath,imagename);
			strcat(imagepath, "/Initial_Concentration.txt");
			f_init_conc=fopen(imagepath,"r");
			fscanf(f_init_conc,"%lf  %lf  %lf",&Concentration[0],&Concentration[1],&Concentration[2]);
			fclose(f_init_conc);
			for(i=0;i<num_of_clusters;i++)
			{
				ClusterCenter[i]=((double)(ClusterCenter[i]*PI)/180.0);
			}
                        strcpy(imagepath,imagename);
                        sprintf(path,"/Iteration_w_%0.2lf.txt",weightage[j]);
                        strcat(imagepath,path);
                        f_iter=fopen(imagepath,"w");
                        strcpy(imagepath,imagename);
                        sprintf(path,"/Centroid_w_%0.2lf.txt",weightage[j]);
                        strcat(imagepath,path);
                        f_centers=fopen(imagepath,"w");
                        strcpy(imagepath,imagename);
                        sprintf(path,"/kappa_w_%0.2lf.txt",weightage[j]);
                        strcat(imagepath,path);
                        f_kappa=fopen(imagepath,"w");
                        fprintf(f_iter,"Initial Cluster Parameters: ");
                        for(i=0;i<num_of_clusters;i++)
                        {
                                fprintf(f_iter,"%0.8lf  ",((ClusterCenter[i]*180.00)/PI));
                        }
                        fprintf(f_iter,"\n");
                        for(i=0;i<num_of_clusters;i++)
                        {
                                fprintf(f_iter,"%0.8lf  ",Concentration[i]);
                        }
                        fprintf(f_iter,"\n");
                        printf("\n*************** AT WEIGHTAGE : %lf ***************\n",weightage[j]);
			AngleFreqSum=0;
			for(i=0;i<AngleLength;i++)
			{
				AngleFreqSum+=Modified_hue_histogram[i];
			}
			SWFreqSum=0.0;
			for(i=0;i<AngleLength;i++)
			{
				SWFreqSum+=SW_hue_histogram[i];
			}
			
			UpdateMembershipValue(Concentration, AngleValue, ClusterCenter, num_of_clusters, MembershipValue, AngleLength, fuzzifier, Modified_hue_histogram, SW_hue_histogram, AngleFreqSum, SWFreqSum);
			threshold=getThresholdValue(MembershipValue, num_of_clusters, AngleLength, Modified_hue_histogram, SW_hue_histogram, AngleFreqSum, SWFreqSum);
			
			printf("\nKI PROBLEM?\n");
			
                        do
                        {
                                flag=0;
                                itercount++;
                                CopyOriginalArray(ClusterCenterOld, ClusterCenter, num_of_clusters);
				UpdateMembershipValue(Concentration,AngleValue,ClusterCenter,num_of_clusters,MembershipValue,AngleLength,fuzzifier, Modified_hue_histogram, SW_hue_histogram, AngleFreqSum, SWFreqSum);
				find_roughsets(MembershipValue, AngleLength, num_of_clusters, lower, boundary, threshold);
				estimate_centroid(MembershipValue,AngleValue,ClusterCenter,num_of_clusters,AngleLength,fuzzifier,lower,boundary,Modified_hue_histogram,SW_hue_histogram,weightage[j]);
				estimate_concentration_kappa(MembershipValue,AngleValue,ClusterCenter,Concentration,num_of_clusters,AngleLength,fuzzifier,lower,boundary,Modified_hue_histogram,SW_hue_histogram,weightage[j]);
                                flag=checkifLargeDifferenceArray(ClusterCenterOld, ClusterCenter, num_of_clusters, epsilon);
                                fprintf(f_iter,"\n\n*************** Iteration Number : %d ***************\n\n",itercount);
                                for(i=0;i<num_of_clusters;i++)
                                {
                                        fprintf(f_iter,"%0.8lf  ",((ClusterCenter[i]*180.00)/PI));
                                }
                                fprintf(f_iter,"\n");
                                for(i=0;i<num_of_clusters;i++)
                                {
                                        fprintf(f_iter,"%0.8lf  ",Concentration[i]);
                                }
                                fprintf(f_iter,"\n");

                        }while(flag==1 && itercount<maxiter);
                        fclose(f_iter);
                        for(i=0;i<num_of_clusters;i++)
                        {
                                fprintf(f_centers,"%0.8lf  ",((ClusterCenter[i]*180.00)/PI));
                        }
                        fclose(f_centers);
                        for(i=0;i<num_of_clusters;i++)
                        {
                                fprintf(f_kappa,"%0.4lf  ",Concentration[i]);
                        }
                        fclose(f_kappa);
                        double_matrix_deallocation(MembershipValue, num_of_clusters);
                        free(ClusterCenter);
                        free(ClusterCenterOld);
                        free(Concentration);
			free(Modified_hue_histogram);
                        free(lower);
                        integer_matrix_deallocation(boundary, AngleLength);
                        itercount=0;
                }
                
                fprintf(f_threshold,"%0.6lf\n",threshold);
                
		free(AngleValue);
                free(imagename);
                free(imagepath);
                free(sourcefile);
                free(path);
        }
        fclose(f_imagenames);
	fclose(f_threshold);
}

void TestfindClusterHSI(FILE *f_imagenames,int maxcount,int num_of_clusters,double *weights,int parameter_length)
{
        int i,j,TotalPixels,count=0;
        double **HSI;
        while(count<maxcount)
        {
                count++;
                printf("\n\n********************RGB domain Mean computation on Image %d is started********************\n", count);
                char *imagename=(char*)malloc(100 * sizeof(char));
                char *imagepath=(char*)malloc(200 * sizeof(char));
                char *sourcefile=(char*)malloc(100 * sizeof(char));
                fscanf(f_imagenames,"%s",sourcefile);
                FILE *f_hueHist,*f_centers,*f_RGB;
                cut_string(sourcefile, (strlen(sourcefile)-4), imagename);
                strcpy(imagepath,imagename);
                strcat(imagepath, "/HSI_values.txt");
                f_hueHist=fopen(imagepath, "r");
                fscanf(f_hueHist,"%d",&TotalPixels);
                HSI=double_matrix_allocation(TotalPixels,num_of_clusters);
                getHSIInformation(f_hueHist, num_of_clusters, TotalPixels, HSI);
                for(i=0;i<parameter_length;i++)
                {
			char *path=(char*)calloc(200,sizeof(char));
			double *hue_center=(double*)calloc(num_of_clusters, sizeof(double));
			double *saturation_center=(double*)calloc(num_of_clusters, sizeof(double));
			double *intensity_center=(double*)calloc(num_of_clusters, sizeof(double));
			double *normalized_red=(double*)calloc(num_of_clusters, sizeof(double));
			double *normalized_green=(double*)calloc(num_of_clusters, sizeof(double));
			double *normalized_blue=(double*)calloc(num_of_clusters, sizeof(double));
			strcpy(imagepath,imagename);
			sprintf(path,"/Centroid_w_%0.2lf.txt",weights[i]);
			strcat(imagepath,path);
			f_centers=fopen(imagepath, "r");
			fscanf(f_centers,"%lf  %lf  %lf",&hue_center[0],&hue_center[1],&hue_center[2]);
			fclose(f_centers);
			printf("\n\tWeight=%0.2lf\n",weights[i]);
			findClusterCenterSaturation(HSI, hue_center, num_of_clusters, TotalPixels, saturation_center, intensity_center);
			HSI_to_RGB_conversion(hue_center, saturation_center, intensity_center, normalized_red, normalized_green, normalized_blue, num_of_clusters);
			convert_to_log_domain(normalized_red, normalized_green, normalized_blue, num_of_clusters);
			free(hue_center);
			strcpy(imagepath,imagename);
			sprintf(path,"/RGB_Mean_w_%0.2lf.txt",weights[i]);
			strcat(imagepath,path);
			f_RGB=fopen(imagepath, "w");
			fprintf(f_RGB,"%lf %lf %lf\n",normalized_red[0],normalized_red[1],normalized_red[2]);
			fprintf(f_RGB,"%lf %lf %lf\n",normalized_green[0],normalized_green[1],normalized_green[2]);
			fprintf(f_RGB,"%lf %lf %lf",normalized_blue[0],normalized_blue[1],normalized_blue[2]);
			fclose(f_RGB);
			free(saturation_center);
			free(intensity_center);
			free(normalized_red);
			free(normalized_green);
			free(normalized_blue);
			free(path);
                }
                double_matrix_deallocation(HSI,TotalPixels);
                free(imagename);
                free(imagepath);
                free(sourcefile);
        }
}

void TestMappingProcedure(FILE *f_imagenames,char *filepath,int maxcount,double *weights,int parameter_length,int num_of_clusters,char *TargetImage)
{
        int i,j,k,l,row,col,maxval,flag,count,num_of_channels, itercount;
	double eps,weight;
        double **targetcenter,**Inverse_centermatrix,**sourcecenter, **sourcecenter_transpose, **sourcecenter_transpose_Old,**temp;
        char ch1,ch2,ch3;
	FILE *f_upto;
	count=0, num_of_channels=3;
	eps=1e-4;
	
        while(count<maxcount)
        {
                count++;
                printf("\n\n********************Mapping Procedure of Image %d is started********************\n", count);
                char *imagename=(char*)malloc(100 * sizeof(char));
                char *imagepath=(char*)malloc(200 * sizeof(char));
                char *sourcefile=(char*)malloc(100 * sizeof(char));
		char *data_path=(char*)calloc(200, sizeof(char));
                FILE *fimage,*fsource,*ftarget,*fnormalized,*fstain,*f_nmf_center;
                fscanf(f_imagenames,"%s",sourcefile);
		strcpy(data_path, filepath);
		strcat(data_path, sourcefile);
                fimage=fopen(data_path,"r");
		fscanf(fimage,"%c%c%c",&ch1,&ch2,&ch3);
		skip_comments(fimage);
		fscanf(fimage,"%d %d",&col,&row);
		fscanf(fimage,"%d",&maxval);
		int **image_red,**image_green,**image_blue;
		double **normalized_image_red,**normalized_image_green,**normalized_image_blue,**ResizedImage,**stain_density,**MappedMatrix,**sourcecenterOld,**stain_densityOld;
		image_red=integer_matrix_allocation(row,col);
		image_green=integer_matrix_allocation(row,col);
		image_blue=integer_matrix_allocation(row,col);
		normalized_image_red=double_matrix_allocation(row,col);
		normalized_image_green=double_matrix_allocation(row,col);
		normalized_image_blue=double_matrix_allocation(row,col);
		ResizedImage=double_matrix_allocation((row*col),num_of_channels);
		read_ppm_write_pgmMatrix_Mapping(fimage,image_red,image_green,image_blue,normalized_image_red,normalized_image_green,normalized_image_blue,col,row);
		GetResizedImageMatrix(normalized_image_red,normalized_image_green,normalized_image_blue,row,col,ResizedImage);
		for(i=0;i<parameter_length;i++)
		{
			char *path=(char*)calloc(200,sizeof(char));
			cut_string(sourcefile, (strlen(sourcefile)-4), imagename);
			strcpy(imagepath,imagename);
			sprintf(path,"/RGB_Mean_w_%0.2lf.txt",weights[i]);
			strcat(imagepath,path);
			fsource=fopen(imagepath, "r");
			sourcecenter=double_matrix_allocation(num_of_channels,num_of_clusters);
			for(k=0;k<num_of_channels;k++)
			{
				for(l=0;l<num_of_clusters;l++)
				{
					fscanf(fsource,"%lf",&sourcecenter[k][l]);
				}
			}
			fclose(fsource);
			sourcecenter_transpose=double_matrix_allocation(num_of_clusters, num_of_channels);
			transpose(sourcecenter, sourcecenter_transpose, num_of_channels, num_of_clusters);
			Inverse_centermatrix=double_matrix_allocation(num_of_channels,num_of_clusters);
			GetInverse(sourcecenter_transpose,Inverse_centermatrix,num_of_clusters);
			stain_density=double_matrix_allocation((row*col),num_of_clusters);
			MatrixMultiplication(ResizedImage, Inverse_centermatrix, (row*col), num_of_channels, num_of_clusters, stain_density);
			stain_densityOld=double_matrix_allocation((row*col),num_of_clusters);
			sourcecenter_transpose_Old=double_matrix_allocation(num_of_clusters,num_of_channels);
			itercount=0;
			do
			{
				flag=0;
				itercount++;
				CopyOriginalMatrix(sourcecenter_transpose_Old, sourcecenter_transpose, num_of_clusters, num_of_channels);
				UpdateFirstMatrix(stain_density, sourcecenter_transpose, ResizedImage, (row*col), num_of_channels, num_of_clusters);
				UpdateSecondMatrix(stain_density, sourcecenter_transpose, ResizedImage, (row*col), num_of_channels, num_of_clusters);
				flag=checkifLargeDifferenceMatrix(sourcecenter_transpose, sourcecenter_transpose_Old, num_of_clusters, num_of_channels, eps);
			}while(flag==1 && itercount<100);
			printf("\nIteration Required : %d\n",itercount);
			strcpy(imagepath,imagename);
			sprintf(path,"/RGB_NMF_Mean_w_%0.2lf.txt",weights[i]);
			strcat(imagepath,path);
			f_nmf_center=fopen(imagepath, "w");
			for(k=0;k<num_of_clusters;k++)
			{
				for(l=0;l<num_of_channels;l++)
				{
					fprintf(f_nmf_center,"%lf ",sourcecenter_transpose[k][l]);
					if(l==num_of_channels-1)
					{
						fprintf(f_nmf_center,"\n");
					}
				}
			}
			fclose(f_nmf_center);
			printf("\n............... IMAGE = %d : NMF based Computation is done for weight=%lf ...............\n", count, weights[i]);
			cut_string(TargetImage,(strlen(TargetImage)-4),imagename);
			strcpy(imagepath,imagename);
			sprintf(path,"/RGB_NMF_Mean_w_%0.2lf.txt",weights[i]);
			strcat(imagepath,path);
			ftarget=fopen(imagepath, "r");
			targetcenter=double_matrix_allocation(num_of_clusters,num_of_channels);
			for(k=0;k<num_of_clusters;k++)
			{
				for(l=0;l<num_of_channels;l++)
				{
					fscanf(ftarget,"%lf",&targetcenter[k][l]);
				}
			}
			fclose(ftarget);
			MappedMatrix=double_matrix_allocation((row*col),num_of_channels);
			MapConcentrationDensity(targetcenter,stain_density,row,col,num_of_clusters,num_of_channels,maxval,MappedMatrix);
			ReformImageChannels(MappedMatrix,image_red,image_green,image_blue,row,col);
			cut_string(sourcefile, (strlen(sourcefile)-4), imagename);
			sprintf(path,"_VMM_w_%0.2lf.ppm",weights[i]);
			strcat(imagename,path);
			sprintf(path,"Normalized_UCSB");
			mkdir(path,0777);
			strcat(path,"/");
			strcpy(imagepath,path);
			strcat(imagepath, imagename);
			fnormalized=fopen(imagepath, "w");
			read_pgmMatrix_write_ppm(image_red,image_green,image_blue,fnormalized,col,row,maxval);
			fclose(fnormalized);
			double_matrix_deallocation(sourcecenter,num_of_channels);
			double_matrix_deallocation(sourcecenter_transpose,num_of_clusters);
			double_matrix_deallocation(sourcecenter_transpose_Old,num_of_clusters);
			double_matrix_deallocation(targetcenter,num_of_clusters);
			double_matrix_deallocation(Inverse_centermatrix,num_of_channels);
			double_matrix_deallocation(stain_density,(row*col));
			double_matrix_deallocation(stain_densityOld,(row*col));
			double_matrix_deallocation(MappedMatrix,(row*col));
			free(path);
			printf("\n***************Mapping based on weight=%lf is done***************\n",weights[i]);
		}
		free(imagename);
		free(imagepath);
		free(sourcefile);
		free(data_path);
		integer_matrix_deallocation(image_red,row);
		integer_matrix_deallocation(image_green,row);
		integer_matrix_deallocation(image_blue,row);
		double_matrix_deallocation(normalized_image_red,row);
		double_matrix_deallocation(normalized_image_green,row);
		double_matrix_deallocation(normalized_image_blue,row);
		double_matrix_deallocation(ResizedImage,(row*col));
		
		f_upto=fopen("Upto_Image.txt","w");
		fprintf(f_upto,"***** Mapping is Done Upto Image : %d *****", count);
		fclose(f_upto);
	}
}
