#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrutil.h"

void bDeviation(char **geno, double *baf, int *numprobes, double *out)
{
  int j;

  for (j=0; j<*numprobes; j++) 
  { 

    if(geno[j][0]==geno[j][1])
    {	
      out[j]=DMIN(baf[j],1-baf[j]);
    }
    else
    {
      if((geno[j][0]!=geno[j][1]) & (geno[j][0]!='N'))
      {
        out[j]=fabs(baf[j]-0.5);
      }
      else
      {
         out[j]=-9;
      }
    }
  }
  return;
}
