#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include"myConst.h"
#include"fold.h"
#include"energy_const.h"
#include"fold_vars.h"
#include"pair_mat.h"
#include"convert_Vienna.h"
#include"params.h"
#include"ParCal.h"
#include<limits.h>
#include"misc.h"
#include"RNAbor.h"

/*Partition function calculation for the RNAhairpin program*/
int GetZB(double ***ZB, double ***ZM, double ***ZM1, int i, int j, int h, char sequence[MAXSIZE]){
  int l,r,k;
  if(h==1)
    {ZB[i][j][h]=exp(-HP_Energy(i,j,S0,sequence+1)/kT);
      ZB[j][i][h]=1;
      for (l=i+1;l<min(i+30,j-5)+1;++l)
        {for(r=max(l+4,j-(30-(l-i)));r<j;++r)
            {if(BP(l,r,sequence))
               {ZB[i][j][h]+=ZB[l][r][h]*exp(-IL_Energy(i,j,l,r,S0)/kT);
                 ZB[j][i][h]+=ZB[r][l][h];
               }
            }
        }
    }
  else
    {for (l=i+1;l<min(i+30,j-5)+1;++l)
        {for(r=max(l+4,j-(30-(l-i)));r<j;++r)
            {if(BP(l,r,sequence))
               {ZB[i][j][h]+=exp(-IL_Energy(i,j,l,r,S0)/kT)*ZB[l][r][h];
                 ZB[j][i][h]+=ZB[r][l][h];
               }
            }
        } 
      for(r=i+6;r<j-4;++r)
        {for(k=1;k<h;++k)
            {ZB[i][j][h]+=exp(-(ML_close+MLbasepairAndAUpenalty(j,i,S0))/kT)*ZM[i+1][r-1][k]*ZM1[r][j-1][h-k];
              ZB[j][i][h]+=ZM[r-1][i+1][k]*ZM1[j-1][r][h-k];
            }
        }
    }
}

int GetZM1(double ***ZB, double ***ZM, double ***ZM1, int i, int j, int h, char sequence[MAXSIZE])
{ int r;
  for (r=i+4;r<j+1;++r)
    { ZM1[i][j][h]+=ZB[i][r][h]*exp(-(ML_base*(j-r)+MLbasepairAndAUpenalty(i,r,S0))/kT);
      ZM1[j][i][h]+=ZB[r][i][h];
    }
}

int GetZM(double ***ZM, double ***ZM1, int i, int j, int h, char sequence[MAXSIZE])
{int r,k;
  for(r=i;r<j-3;++r)
    {ZM[i][j][h]+=ZM1[r][j][h]*exp(-ML_base*(r-i)/kT);
      ZM[j][i][h]+=ZM1[j][r][h];}
  for(r=i+5;r<j-3;++r)
    { for(k=1;k<h;++k)
        {ZM[i][j][h]+=ZM[i][r-1][k]*ZM1[r][j][h-k];
          ZM[j][i][h]+=ZM[r-1][i][k]*ZM1[j][r][h-k];
        }
    }
}

int GetZ(int j, int h, double ***ZB,double **Z, double **N, char sequence[MAXSIZE])
{int r,k;
 if(j<5)
    {Z[j][h]=0;    // the empty structure is no longer counted in when h!=0 since it doesn't contain hairpins.
      N[j][h]=0;
    }
  else
    { Z[j][h]+=Z[j-1][h];
      N[j][h]+=N[j-1][h];
      for(r=1;r<j-3;++r)
        {if(BP(r,j,sequence))
           {if(r==1)
               {Z[j][h]+=ZB[r][j][h]*exp(-AU_Penalty(1,j,S0)/kT);
                 N[j][h]+=ZB[j][r][h];}
             else
               {for(k=0;k<h;++k)
                   {Z[j][h]+=Z[r-1][k]*ZB[r][j][h-k]*exp(-AU_Penalty(r,j,S0)/kT);
                       N[j][h]+=N[r-1][k]*ZB[j][r][h-k];
                   }
               }
           }
        }
    }
}



