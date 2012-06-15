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
#include"misc.h"
#include<limits.h>

/*The following function checks if the ith and jth positions in the sequence could base pair*/
int BP(int i,int j,char sequence[MAXSIZE])  //safe
{if(j-i<4)
    return 0;
  else if ((sequence[i]=='A'&& sequence[j]=='U')||(sequence[i]=='U'&& sequence[j]=='A'))
    return 1;
  else if ((sequence[i]=='C'&& sequence[j]=='G')||(sequence[i]=='G'&& sequence[j]=='C'))
    return 1;
  else if ((sequence[i]=='U'&& sequence[j]=='G')||(sequence[i]=='G'&& sequence[j]=='U'))
    return 1;
  else
    return 0;
}

int ElementInList(char subshape[SHAPELENGTH],char shapeList[SHAPELENGTH][SHAPELENGTH],int listlen) //check if a subshape is in list, if yes,return the index of the subshape,if not return -1
{int i;
 for(i=0;i<listlen;++i)
    {if(strcmp(shapeList[i],subshape)==0)
        return i;
    }
 return -1;
}

int NumInList(int a, int L[SHAPELENGTH],int listlen) //check if a number is in this list,if not return -1
{int i;
 for(i=0;i<listlen;++i)
    {if(a==L[i])
        return i;
    }
 return -1;
}

/*This function returns the largest value in an array.*/
double MaxInArray(double array[],int arraylen)   //safe
{ double largest=array[0];
  int i;
  for (i=0;i<arraylen;++i)
    {if (array[i]>largest)
        largest=array[i];
    }
  return largest;
}

/* The following function check if the input strings are valid RNA sequence, structure */
int CheckSequenceAndStructure(char sequence[MAXSIZE], char structure[MAXSIZE]) {
  int i;
  
  // Sequence is 1 indexed, structure is 0 indexed
  if (strlen(sequence) != strlen(structure) + 1) {
    printf("The RNA sequence and structure must be of identical length.\n");
    exit(1);
  }
  
  for (i = 1; i < strlen(sequence); ++i) {
    // Sequence validation
    if (
      toupper(sequence[i]) != 'A' && 
      toupper(sequence[i]) != 'U' && 
      toupper(sequence[i]) != 'C' && 
      toupper(sequence[i]) != 'G' && 
      toupper(sequence[i]) != 'T'
    ) {
      printf("The input sequence should only contain A, U, T, C, G!\n");
      exit(1);
    } else if (toupper(sequence[i]) == 'T') {
      sequence[i] = 'U';
    } else {
      sequence[i] = toupper(sequence[i]);
    }
    
    // Structure validation
    if (
      structure[i - 1] != '(' && 
      structure[i - 1] != ')' && 
      structure[i - 1] != '.'
    ) {
      printf("The input structure should only contain '.', '(' and ')'!\n");
      exit(1);
    }
  }
  
  return 0;
}

double ***Allocate3DMatrix(int a, int b, int c)
// allocate a matrix of size a*b*c
{ int i;
  int j;
  double ***Matrix;
  Matrix=(double ***) malloc((a) * sizeof (double **));
  if(Matrix == NULL)
        {printf("out of memory\n");
         exit(1);}
  for (i=0; i<a; i++){
    Matrix[i]=(double **) malloc((b)*sizeof(double *));
    if(Matrix[i] == NULL)
        {printf("out of memory\n");
         exit(1);}
  }
  for(i=0;i<a;i++)
    {for(j=0;j<b;j++)
        {
        Matrix[i][j]=(double *) malloc((c)*sizeof(double));
        if(Matrix[i][j] == NULL)
        {printf("out of memory\n");
         exit(1);}
        }
    }
  return Matrix;
}

int Free3DMatrix(double ***Matrix, int a, int b, int c)
// Free a matrix of size a*b*c
{ int i;
  int j;
  for(i=0;i<a;i++)
	{
	for(j=0;j<b;j++)
	free(Matrix[i][j]);
	}
  for(i=0;i<a;i++)
	free(Matrix[i]);
  free(Matrix);
  return 0;
}

double **Allocate2DMatrix(int a, int b)
// allocate a matrix of size a*b
{ int i;
  double **Matrix;
  Matrix=(double **) calloc(a, sizeof (double *));
  if(Matrix == NULL)
        {printf("out of memory\n");
         exit(1);}
  for (i=0; i<a; i++){
    Matrix[i]=(double *) calloc(b, sizeof(double));
    if(Matrix[i] == NULL)
        {printf("out of memory\n");
         exit(1);}
  }
  return Matrix;
}

int Free2DMatrix(double **Matrix,int a, int b)
// free a matrix of size a*b
{ int i;
  for (i=0; i<a; i++){
    free(Matrix[i]);
  }
  free(Matrix);
  return 0;
}





