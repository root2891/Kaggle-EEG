
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>


#define VERY_BIG  (1e30)
void dtw(int* size__of_file, int *no_of_params, double* result)

{
  Rprintf("file loaded\n");
int size_file = *size__of_file;
int params = *no_of_params;
double **globdist;
double **Dist;

double  total;

unsigned int i, j, k;

FILE *file1;
double **x;
if ((file1=fopen("D:/Work Related/Kaggle EEG/table2.txt","rb"))==NULL)
{fprintf(stderr,"File table1 cannot be opened\n");
exit(1);
}
if ((x = malloc(size_file * sizeof(double *))) == NULL)
     fprintf(stderr,"Memory allocation error (x)\n");

for (i=0; i < size_file; i++)
     if ((x[i] = malloc(params * sizeof(double))) == NULL)
     fprintf(stderr,"Memory allocation error (x)\n");
double myvariable;
for(i = 0; i < size_file; i++)
  {
    for (j = 0 ; j < params; j++)
    {
      fscanf(file1,"%lf",&myvariable);
      x[i][j] = myvariable;
      Rprintf("val: %.0f\n",x[i][j]);
    }
  }
     /* allocate memory for Dist */
if ((Dist = malloc(params * sizeof(double *))) == NULL)
     fprintf(stderr,"Memory allocation error (Dist)\n");

for (i=0; i < params; i++)
if ((Dist[i] = malloc(params * sizeof(double))) == NULL)
     fprintf(stderr,"Memory allocation error (Dist)\n");

     /* allocate memory for globdist */
if ((globdist = malloc(params * sizeof(double *))) == NULL)
     fprintf(stderr,"Memory allocation error (globdist)\n");

for (i=0; i < params; i++)
if ((globdist[i] = malloc(params * sizeof(double))) == NULL)
     fprintf(stderr,"Memory allocation error (globdist)\n");
int a = 0, b= 0;
int count = 0;
for(a = 0; a<size_file; a++)
{
  for(b = 0; b<size_file; b++)
  {
    if(a<b)
    {

      /*Compute distance matrix*/

      for(i=0;i<params;i++) 
      {
        for(j=0;j<params;j++) 
        {
          total = 0;
            total = total + abs((x[a][i] - x[b][j]));
          Dist[i][j] = total;
        }
      }


      globdist[0][0] = Dist[0][0];

      for(i=1; i<params; i++)
      { 
          globdist[i][0] = Dist[i][0] + globdist[i-1][0];
      }

      for(k=1; k<params; k++) 
      {
          globdist[0][k] = Dist[0][k] + globdist[0][k-1];
      }
        
          
      for(i=1; i<params; i++)
      {
          for(k=1; k<params; k++)
          {
             if(globdist[i-1][k] < globdist[i][k-1])
             {
                    globdist[i][k] = Dist[i][k] +    globdist[i-1][k];         
             }
             else
             {
                  globdist[i][k] = Dist[i][k] +   globdist[i][k-1];
             }
          }
      }
      *(result+count) = globdist[params-1][params-1];
      count++;
    }
  }
}
Rprintf("result: %.0f\n",result[0]);
fclose(file1);
for(i =0; i<size_file; i++)
      free(x[i]);
   free(x);   

for(i =0; i<params; i++)
      free(Dist[i]);
      free(Dist);     
for(i =0; i<params; i++)
      free(globdist[i]);
      free(globdist);
}
