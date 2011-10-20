#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h> 
#include"fold.h"
#include"energy_const.h"
#include"fold_vars.h"
#include"pair_mat.h"
#include"convert_Vienna.h"
#include"params.h"
#include"myConst.h"
#include"misc.h"
#include"RNAhairpin.h"
#include"McCaskill.h"
#include"limits.h"

PRIVATE short *encode_seq(const char *seq);
int MaxNumHairpin(char sequence[MAXSIZE], double **HPMAX);
double HairpinPartition(double *HP,double *HN, int H_MAX, char sequence[MAXSIZE]);
/*This is an algorithm to calculate the partition function with respect to the number of hairpins under Turner Energy model.*/
/*Author: Yang Ding; Date: 07/04/2009*/

int main(int argc, char *argv[]){
  char sequence[MAXSIZE],filename[FILENAMELENGTH];
  FILE *inputfile;
  int i,j,h;
  double **HPMAX;
  double HP[MAXSIZE],HN[MAXSIZE];
  int H_MAX=INT_MAX;
  double totalpar; //total partition function
  if (argc==3||argc==5)
    {
      i=1;
      while(i<argc)
	{
	  if(!strcmp(argv[i], "-s"))
	    {
	      sequence[0]='@';
	      strncpy(sequence+1,argv[i+1],strlen(argv[i+1])); // the sequence starts from index 1
	      sequence[strlen(argv[i+1])+1]='\0'; //remember to end the sequence
	      i+=2;
	    }
	  else if(!strcmp(argv[i], "-n"))
	    {
	      H_MAX=atoi(argv[i+1]);
	      i+=2;
	    }
	  else if(argc==3)
	    {
	     i+=2;
	    }
	  else
	    {
	      printf("Error: Unrecognizable Flag!\n");
	      printf("Usage:");
	      printf("%s -s SEQUENCE [-n NUMBER]\n",argv[0]);
	      printf("SEQUENCE is the RNA sequence.\n");
	      printf("NUMBER is the maximum number of hairpin that one considers.\n");
	      printf("If not specified, the program will calculate the maximum possible number of hairpins given the current RNA sequence, and use that number as an upper bound.\n");
	      printf("The length of RNA sequence should be less than %d and larger than or equal to 5.\n",MAXSIZE);
	      exit(1);
	    }
	}
    CheckSequence(sequence);
    S0=encode_seq(sequence+1);
    seqlen=strlen(sequence)-1;
    Initialize_Params();
    make_pair_matrix();//needed for pair matching
    kT = (temperature+K0)*GASCONST/1000.0;
    //printf("%.15f\n",kT);
    ML_base=(double)P->MLbase/100;
    ML_close=(double)P->MLclosing/100;
    HPMAX=Allocate2DMatrix( seqlen+1, seqlen+1);
    for(i=1;i<seqlen+1;++i)
      {for (j=1;j<seqlen+1;++j)
          HPMAX[i][j]=0;
      }
    MaxNumHairpin(sequence,HPMAX);
    if(H_MAX>HPMAX[1][seqlen])
	H_MAX=HPMAX[1][seqlen];
    for (i=0;i<H_MAX+1;++i)
      {HP[i]=0;
        HN[i]=0;}
    HairpinPartition(HP,HN,H_MAX,sequence);
    
    char testString[5];
    testString[0] = '.';
    testString[1] = '.';
    testString[2] = '.';
    testString[3] = '.';
    testString[4] = '.';
    
    double** McCaskillZ;
    McCaskillZ=runMcCaskill(sequence);
    totalpar=McCaskillZ[1][seqlen];
    printf("The RNA sequence is %s:\n",sequence+1);
    printf("The total number of structures is %.0f.\n", McCaskillZ[seqlen][1]);
    printf("The first column is the number of hairpins in a secondary structure.\n");
    printf("The partition function and number of secondary structure with k-hairpin is:\n");
    printf("%-20s%-30s%-30s%-30s\n","K-hairpin","Number","Partition Function","Probability");
    for (h=0;h<H_MAX+1;++h)
      printf("%-8d\t%-30.15lf\t%-30.15lf\t%-30.15lf\n",h,HN[h],HP[h],HP[h]/totalpar);
    }
  else{
    printf("Usage:");
    printf("%s -s SEQUENCE [-n NUMBER]\n",argv[0]);
    printf("SEQUENCE is the RNA sequence\n");
    printf("NUMBER is the maximum number of hairpin that one considers.\n");
    printf("If not specified, the program will calculate the maximum possible number of hairpins given the current RNA sequence, and use that number as an upper bound\n");
    printf("The length of RNA sequence should be less than %d and larger than or equal to 5\n",MAXSIZE);
    exit(1);
  }
  return 0;
}

/*Below modified from fold.c*/
PRIVATE short *encode_seq(const char *seq) {
  unsigned int k,l;
  short *S0_out;
  int i;
  l = strlen(seq);
  S0_out = (short *) space(sizeof(short)*(l+2));
  S0_out[0]=l;
  for (k=1; k<=l; k++) { /* make numerical encoding of seq */
    S0_out[k]= (short) encode_char(toupper(seq[k-1]));
  }
  return S0_out;
}

/*This function calculates the maximum number of hairpins given a sequence.*/
int MaxNumHairpin(char sequence[MAXSIZE], double **HPMAX){
  int d,i,j,r;
  for (d=4;d<seqlen;++d)
    {for(i=1;i<seqlen-d+1;++i)
        {j=i+d;
          HPMAX[i][j]=HPMAX[i][j-1];
          for(r=i;r<j-3;++r)
            {if (BP(r,j,sequence))
                { if(r==i)
                    {double S1[3]={1.,HPMAX[i][j],HPMAX[i+1][j-1]};
                      HPMAX[i][j]=MaxInArray(S1,3);}
                  else
                    {double S2[2]={1.,HPMAX[r+1][j-1]};
		      double S3[3]={1.,HPMAX[i][j],HPMAX[i][r-1]+MaxInArray(S2,2)};
                      HPMAX[i][j]=MaxInArray(S3,3);
                    }
                }
            }
        }
    }
  return 0;
}

/*The following function is the main recursion*/
double HairpinPartition(double *HP,double *HN, int H_MAX, char sequence[MAXSIZE])
{ double ***ZB;
  double ***ZM;
  double ***ZM1;
  double **Z;
  double **N;
  Z=Allocate2DMatrix( seqlen+1, H_MAX+1);
  N=Allocate2DMatrix( seqlen+1, H_MAX+1);
  ZB=Allocate3DMatrix( seqlen+1,seqlen+1, H_MAX+1);
  ZM=Allocate3DMatrix( seqlen+1,seqlen+1, H_MAX+1);
  ZM1=Allocate3DMatrix( seqlen+1,seqlen+1, H_MAX+1);
  int i,j,r,d,h;
  for(i=1;i<seqlen+1;++i)
    {for (j=1;j<seqlen+1;++j)
	{for(r=1;r<H_MAX+1;++r)
	    { ZB[i][j][r]=0;
	      ZM[i][j][r]=0;
	      ZM1[i][j][r]=0;
	    }
	}
    }
  for (i=1;i<seqlen+1;++i)
    {for(j=1;j<H_MAX+1;++j)
	{Z[i][j]=0;
	  N[i][j]=0; 
	}
    }
  for(h=0;h<H_MAX+1;++h)
    {if(h==0)
	{
	  HP[h]=1.0;
	  HN[h]=1.0;
	  for(i=0;i<seqlen+1;++i)
	    {Z[i][0]=1;
	      N[i][0]=1;
	    }
	}
      else
	{ 
	  for(d=4;d<seqlen;++d)
	    {for(i=1;i<seqlen-d+1;++i)
		{j=i+d;
		  if(BP(i,j,sequence))
		    GetZB(ZB,ZM,ZM1,i,j,h,sequence);
		  GetZM1(ZB,ZM,ZM1,i,j,h,sequence);
		  GetZM(ZM,ZM1,i,j,h,sequence);
		}
	    }
	  for(j=1;j<seqlen+1;++j)
	    {GetZ(j,h,ZB,Z,N,sequence);
	    }
	  HP[h]=Z[seqlen][h];
	  HN[h]=N[seqlen][h];
	}
    }
}

