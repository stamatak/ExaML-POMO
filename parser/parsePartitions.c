/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 * 
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */


#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h> 
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <inttypes.h>



#include "axml.h"

/*****************************FUNCTIONS FOR READING MULTIPLE MODEL SPECIFICATIONS************************************************/


extern char modelFileName[1024];
extern char excludeFileName[1024];
extern char proteinModelFileName[1024];
extern char secondaryStructureFileName[1024];


extern char seq_file[1024];

extern char *protModels[NUM_PROT_MODELS];

boolean lineContainsOnlyWhiteChars(char *line)
{
  size_t 
    i, 
    n = strlen(line);

  if(n == 0)
    return TRUE;

  for(i = 0; i < n; i++)
    {
      if(!whitechar(line[i]))
	return FALSE;
    }
  return TRUE;
}


static int isNum(int c)
{
  
  return (c == '0' || c == '1' || c == '2' || c == '3' || c == '4' ||
	  c == '5' || c == '6' || c == '7' || c == '8' || c == '9');
}


static void skipWhites(char **ch)
{
  while(**ch == ' ' || **ch == '\t')
    *ch = *ch + 1;
}

static void analyzeIdentifier(char **ch, int modelNumber, tree *tr)
{
  char
    *start = *ch,
    ident[2048] = "";
  char model[128] = "";  
  char thisModel[1024];
  size_t 
    i = 0, n, j;
  int containsComma = 0;

  while(**ch != '=')
    {
      if(**ch == '\n' || **ch == '\r')
	{
	  printf("\nPartition file parsing error!\n");
	  printf("Each line must contain a \"=\" character\n");
	  printf("Offending line: %s\n", start);
	  printf("ExaML will exit now.\n\n");
	  exit(-1);
	}

      if(**ch != ' ' && **ch != '\t')
	{
	  ident[i] = **ch;      
	  i++;
	}
      *ch = *ch + 1;
    }
  
  n = i;
  i = 0;
  
  for(i = 0; i < n; i++)
    if(ident[i] == ',') 
      containsComma = 1;

  if(!containsComma)
    {
      printf("Error, model file must have format: BIN, DNA, AA, or POMO model, then a comma, and then the partition name\n");
      exit(-1);
    }
  else
    {
      boolean found = FALSE;
      i = 0;
      while(ident[i] != ',')
	{
	  model[i] = ident[i];
	  i++;
	}      
      
      /* AA */

      for(i = 0; i < NUM_PROT_MODELS && !found; i++)
	{	
	  strcpy(thisModel, protModels[i]);
	  
	  if(strcasecmp(model, thisModel) == 0)
	    {	      	      
	      tr->initialPartitionData[modelNumber].protModels = (int)i;		  
	      tr->initialPartitionData[modelNumber].protFreqs  = 0;
	      tr->initialPartitionData[modelNumber].dataType   = AA_DATA;
	      found = TRUE;
	    }
	  	  
	  strcpy(thisModel, protModels[i]);
	  strcat(thisModel, "F");
	  
	  if(strcasecmp(model, thisModel) == 0)
	    {	      
	      tr->initialPartitionData[modelNumber].protModels = (int)i;		  
	      tr->initialPartitionData[modelNumber].protFreqs  = 1;
	      tr->initialPartitionData[modelNumber].dataType   = AA_DATA;
	      found = TRUE;

	      if(tr->initialPartitionData[modelNumber].protModels == AUTO)
		{
		  printf("\nError: Option AUTOF has been deprecated, exiting\n\n");
		  errorExit(-1);
		}
	      
	      if(tr->initialPartitionData[modelNumber].protModels == LG4M || tr->initialPartitionData[modelNumber].protModels == LG4X)
		{
		  printf("\nError: Options LG4MF and LG4XF have been deprecated.\n");
		  printf("They shall only be used with the given base frequencies of the model, exiting\n\n");
		  errorExit(-1);
		}
	    }	

	  strcpy(thisModel, protModels[i]);
	  strcat(thisModel, "X");

	  if(strcasecmp(model, thisModel) == 0)
	    {	      
	      tr->initialPartitionData[modelNumber].protModels = (int)i;		  
	      tr->initialPartitionData[modelNumber].protFreqs  = 0;
	      tr->initialPartitionData[modelNumber].optimizeBaseFrequencies = TRUE;
	      tr->initialPartitionData[modelNumber].dataType   = AA_DATA;
	      found = TRUE;

	      if(tr->initialPartitionData[modelNumber].protModels == AUTO)
		{
		  printf("\nError: Option AUTOX has been deprecated, exiting\n\n");
		  errorExit(-1);
		}
	      
	      if(tr->initialPartitionData[modelNumber].protModels == LG4M || tr->initialPartitionData[modelNumber].protModels == LG4X)
		{
		  printf("\nError: Options LG4MX and LG4XX have been deprecated.\n");
		  printf("They shall only be used with the given base frequencies of the model, exiting\n\n");
		  errorExit(-1);
		}

	    }	

	  /*if(found)
	    printf("%s %d\n", model, i);*/
	}
      
      if(!found)
	{		  	  
	  if(strcasecmp(model, "DNA") == 0)
	    {	     	      
	      tr->initialPartitionData[modelNumber].protModels = -1;		  
	      tr->initialPartitionData[modelNumber].protFreqs  = -1;
	      tr->initialPartitionData[modelNumber].dataType   = DNA_DATA;
	      tr->initialPartitionData[modelNumber].optimizeBaseFrequencies = FALSE;
	      found = TRUE;
	    }
	  else
	    {
	      if(strcasecmp(model, "DNAX") == 0)
		{	     	      
		  tr->initialPartitionData[modelNumber].protModels = -1;		  
		  tr->initialPartitionData[modelNumber].protFreqs  = -1;
		  tr->initialPartitionData[modelNumber].dataType   = DNA_DATA;
		  tr->initialPartitionData[modelNumber].optimizeBaseFrequencies = TRUE;
		  found = TRUE;
		}	      	    
	      else
		{	    	  
		  if(strcasecmp(model, "BIN") == 0)
		    {	     	      
		      tr->initialPartitionData[modelNumber].protModels = -1;		  
		      tr->initialPartitionData[modelNumber].protFreqs  = -1;
		      tr->initialPartitionData[modelNumber].dataType   = BINARY_DATA;
		      tr->initialPartitionData[modelNumber].optimizeBaseFrequencies = FALSE;
		      found = TRUE;
		    }
		  else
		    {
		      if(strcasecmp(model, "BINX") == 0)
			{
			  tr->initialPartitionData[modelNumber].protModels = -1;		  
			  tr->initialPartitionData[modelNumber].protFreqs  = -1;
			  tr->initialPartitionData[modelNumber].dataType   = BINARY_DATA;
			  tr->initialPartitionData[modelNumber].optimizeBaseFrequencies = TRUE;
			  found = TRUE;
			}
		      else
			{
			  if(strcasecmp(model, "POMO16") == 0)
			    {
			      tr->initialPartitionData[modelNumber].protModels = -1;		  
			      tr->initialPartitionData[modelNumber].protFreqs  = -1;
			      tr->initialPartitionData[modelNumber].dataType   = POMO_16;
			      tr->initialPartitionData[modelNumber].optimizeBaseFrequencies = FALSE;
			      
			      found = TRUE;
			    }
			  else
			    {
			      if(strcasecmp(model, "POMO16X") == 0)
				{
				  tr->initialPartitionData[modelNumber].protModels = -1;		  
				  tr->initialPartitionData[modelNumber].protFreqs  = -1;
				  tr->initialPartitionData[modelNumber].dataType   = POMO_16;
				  tr->initialPartitionData[modelNumber].optimizeBaseFrequencies = TRUE;
				  
				  found = TRUE;
				}
			      else
				{
				  if(strcasecmp(model, "POMO64") == 0)
				    {
				      tr->initialPartitionData[modelNumber].protModels = -1;		  
				      tr->initialPartitionData[modelNumber].protFreqs  = -1;
				      tr->initialPartitionData[modelNumber].dataType   = POMO_16;
				      tr->initialPartitionData[modelNumber].optimizeBaseFrequencies = FALSE;				      
				      found = TRUE;
				    }
				  else
				    {
				      if(strcasecmp(model, "POMO64X") == 0)
					{
					   tr->initialPartitionData[modelNumber].protModels = -1;		  
					   tr->initialPartitionData[modelNumber].protFreqs  = -1;
					   tr->initialPartitionData[modelNumber].dataType   = POMO_64;
					   tr->initialPartitionData[modelNumber].optimizeBaseFrequencies = TRUE;
					   
					  found = TRUE;
					}
				      else
					{			      
					  if(strcasecmp(model, "MULTI") == 0)
					    {	     	      
					      tr->initialPartitionData[modelNumber].protModels = -1;		  
					      tr->initialPartitionData[modelNumber].protFreqs  = -1;
					      tr->initialPartitionData[modelNumber].dataType   = GENERIC_32;
					      
					      found = TRUE;
					    }
					  else
					    {
					      if(strcasecmp(model, "CODON") == 0)
						{	     	      
						  tr->initialPartitionData[modelNumber].protModels = -1;		  
						  tr->initialPartitionData[modelNumber].protFreqs  = -1;
						  tr->initialPartitionData[modelNumber].dataType   = GENERIC_64;
						  
						  found = TRUE;
						}
					    }
					}
				    }				    				
				}
			    }
			}
		    }
		}
	    }
	}

      if(!found)
	{
	  printf("ERROR: you specified the unknown model %s for partition %d\n", model, modelNumber);
	  exit(-1);
	}
           

      i = 0;
      while(ident[i++] != ',');      

      tr->initialPartitionData[modelNumber].partitionName = (char*)malloc((n - i + 1) * sizeof(char));          

      j = 0;
      while(i < n)	
	tr->initialPartitionData[modelNumber].partitionName[j++] =  ident[i++];

      tr->initialPartitionData[modelNumber].partitionName[j] = '\0';                      
    }
}



static void setModel(int model, int64_t position, int *a)
{
  if(a[position] == -1)
    a[position] = model;
  else
    {
      printf("ERROR trying to assign model %d to position %" PRId64 "\n", model, position);
      printf("while already model %d has been assigned to this position\n", a[position]);
      exit(-1);
    }      
}




void parsePartitions(analdef *adef, rawdata *rdta, tree *tr)
{
  FILE *f; 
  int numberOfModels = 0; 
  size_t nbytes = 0;
  char *ch;
  char *cc = (char *)NULL;
  char **p_names;
  size_t n, l;
  int64_t lower, upper, modulo;
  char buf[256];
  int64_t **partitions;
  int pairsCount;
  int64_t i, as;
  int64_t j;
  int64_t k; 
  char* endptr;

  f = myfopen(modelFileName, "rb");   

 
  while(getline(&cc, &nbytes, f) > -1)
    {     
      if(!lineContainsOnlyWhiteChars(cc))	
	numberOfModels++;
       
      if(cc)
	free(cc);
      
      cc = (char *)NULL;
    }     
  
  rewind(f);
      
  p_names = (char **)malloc(sizeof(char *) * (size_t)numberOfModels);
  partitions = (int64_t **)malloc(sizeof(int64_t *) * (size_t)numberOfModels);     
  
  tr->initialPartitionData = (pInfo*)malloc(sizeof(pInfo) * (size_t)numberOfModels);
      
  for(i = 0; i < numberOfModels; i++) 
    {     
      tr->initialPartitionData[i].protModels = adef->proteinMatrix;
      tr->initialPartitionData[i].protFreqs  = adef->protEmpiricalFreqs;
      tr->initialPartitionData[i].dataType   = -1;
    }

  for(i = 0; i < (int64_t)numberOfModels; i++)    
    partitions[i] = (int64_t *)NULL;
    
  i = 0;
  while(getline(&cc, &nbytes, f) > -1)
    {          
      if(!lineContainsOnlyWhiteChars(cc))
	{
	  n = strlen(cc);	 
	  p_names[i] = (char *)malloc(sizeof(char) * ((size_t)n + 1));
	  strcpy(&(p_names[i][0]), cc);
	  i++;
	}
      if(cc)
	free(cc);
      cc = (char *)NULL;
    }         

  for(i = 0; i < numberOfModels; i++)
    {           
      ch = p_names[i];     
      pairsCount = 0;
      skipWhites(&ch);
      
      if(*ch == '=')
	{
	  printf("Identifier missing prior to '=' in %s\n", p_names[i]);
	  exit(-1);
	}
      
      analyzeIdentifier(&ch, (int)i, tr);
      ch++;
            
    numberPairs:
      pairsCount++;
      partitions[i] = (int64_t *)realloc((void *)partitions[i], (1 + 3 * (size_t)pairsCount) * sizeof(int64_t));
      partitions[i][0] = pairsCount;
      partitions[i][3 + 3 * (pairsCount - 1)] = -1; 	
      
      skipWhites(&ch);
      
      if(!isNum(*ch))
	{
	  printf("%c Number expected in %s\n", *ch, p_names[i]);
	  exit(-1);
	}   
      
      l = 0;
      while(isNum(*ch))		 
	{       
	  buf[l] = *ch;
	  ch++;	
	  l++;
	}
      buf[l] = '\0';
      lower = strtoimax(buf, &endptr, 10);      
      partitions[i][1 + 3 * (pairsCount - 1)] = lower;   
      
      skipWhites(&ch);
      
      /* NEW */
      
      if((*ch != '-') && (*ch != ','))
	{
	  if(*ch == '\0' || *ch == '\n' || *ch == '\r')
	    {
	      upper = lower;
	      goto SINGLE_NUMBER;
	    }
	  else
	    {
	      printf("'-' or ',' expected in %s\n", p_names[i]);
	      exit(-1);
	    }
	}	 
      
      if(*ch == ',')
	{	     
	  upper = lower;
	  goto SINGLE_NUMBER;
	}
      
      /* END NEW */
      
      ch++;   
      
      skipWhites(&ch);
      
      if(!isNum(*ch))
	{
	  printf("%c Number expected in %s\n", *ch, p_names[i]);
	  exit(-1);
	}    
      
      l = 0;
      while(isNum(*ch))
	{    
	  buf[l] = *ch;
	  ch++;	
	  l++;
	}
      buf[l] = '\0';
      upper = strtoimax(buf, &endptr, 10); 
    SINGLE_NUMBER:
      partitions[i][2 + 3 * (pairsCount - 1)] = upper;        	  
      
      if(upper < lower)
	{
	  printf("Upper bound %" PRId64 " smaller than lower bound %" PRId64 " for this partition: %s\n", upper, lower,  p_names[i]);
	  exit(-1);
	}
      
      skipWhites(&ch);
      
      if(*ch == '\0' || *ch == '\n' || *ch == '\r') /* PC-LINEBREAK*/
	{    
	  goto parsed;
	}
      
      if(*ch == ',')
	{	 
	  ch++;
	  goto numberPairs;
	}
      
      if(*ch == '\\')
	{
	  ch++;
	  skipWhites(&ch);	 	

	  if(!isNum(*ch))
	    {
	      printf("%c Number expected in %s\n", *ch, p_names[i]);
	      exit(-1);
	    }   

	  if(adef->compressPatterns == FALSE)
	    {
	      printf("\nError: You are not allowed to use interleaved partitions, that is, assign non-contiguous sites\n");
	      printf("to the same partition model, when pattern compression is disabled via the -c flag!\n\n");
	      exit(-1);
	    }
	  
	  l = 0;
	  while(isNum(*ch))
	    {
	      buf[l] = *ch;
	      ch++;	
	      l++;
	    }
	  buf[l] = '\0';
	  modulo = strtoimax(buf, &endptr, 10);
	  partitions[i][3 + 3 * (pairsCount - 1)] = modulo; 	
	  
	  skipWhites(&ch);
	  if(*ch == '\0' || *ch == '\n' || *ch == '\r')
	    {	     
	      goto parsed;
	    }
	  if(*ch == ',')
	    {	       
	      ch++;
	      goto numberPairs;
	    }
	}  
      
      
      printf("\nError: You may be using \"/\" for specifying interleaved partitions in the model file, while it should be \"\\\" !\n\n");
      assert(0);
       
    parsed:
      ;
    }
  
  fclose(f);
 
  /*********************************************************************************************************************/ 

  for(i = 0; i <= rdta->sites; i++)
    tr->model[i] = -1;
  
  for(i = 0; i < (int64_t)numberOfModels; i++)
    {   
      as = partitions[i][0];     
      
      for(j = 0; j < as; j++)
	{
	  lower = partitions[i][1 + j * 3];
	  upper = partitions[i][2 + j * 3]; 
	  modulo = partitions[i][3 + j * 3];	
	 
	  if(modulo == -1)
	    {
	      for(k = lower; k <= upper; k++)
		setModel((int)i, k, tr->model);
	    }
	  else
	    {
	      for(k = lower; k <= upper; k += modulo)
		{
		  if(k <= rdta->sites)
		    setModel((int)i, k, tr->model);	      
		}
	    }
	}        
    }


  for(i = 1; i < rdta->sites + 1; i++)
    {      
      if(tr->model[i] == -1)
	{
	  printf("ERROR: Alignment Position %zu has not been assigned any model\n", i);
	  exit(-1);
	}      
    }  

  {
    int 
      pomo16 = 0,
      pomo64 = 0;

    for(i = 0; i < numberOfModels; i++)
      {
	if(tr->initialPartitionData[i].dataType == POMO_16)
	  pomo16++;
	if(tr->initialPartitionData[i].dataType == POMO_64)
	  pomo64++;
      }    

    if(pomo16 > 0)
      {
	if(pomo16 < numberOfModels)
	  {
	    printf("\nError: When using POMO all partitions either need to use POMO_16 or POMO_64\n\n");
	    errorExit(-1);	    
	  }

	if(pomo16 == numberOfModels && adef->model !=  M_POMOGAMMA_16)
	  {
	    printf("\nError, for using a partitioned POMO16 model you also need to specify POMO in the command line via -m POMO16\n\n");
	    errorExit(-1);
	  }
      }

    if(pomo64 > 0)
      {
	if(pomo64 < numberOfModels)
	  {
	    printf("\nError: When using POMO all partitions either need to use POMO_16 or POMO_64\n\n");
	    errorExit(-1);	    
	  }

	if(pomo64 == numberOfModels && adef->model !=  M_POMOGAMMA_64)
	  {
	    printf("\nError, for using a partitioned POMO64 model you also need to specify POMO in the command line via -m POMO64\n\n");
	    errorExit(-1);
	  }
      }

  }

  for(i = 0; i < numberOfModels; i++)
    {
      free(partitions[i]);
      free(p_names[i]);
    }
  
  free(partitions);
  free(p_names);    
    
  tr->NumberOfModels = numberOfModels;     
  
  
}

/*******************************************************************************************************************************/


void parseProteinModel(analdef *adef)
{
  FILE *f; 
  int doublesRead = 0;
  int result = 0;
  int i, j;
  double acc = 0.0;

  assert(adef->userProteinModel);
  printf("User-defined prot mod %s\n", proteinModelFileName);

  adef->externalAAMatrix = (double*)malloc(420 * sizeof(double));

  f = myfopen(proteinModelFileName, "rb");
  
 

  while(doublesRead < 420)
    {     
      result = fscanf(f, "%lf", &(adef->externalAAMatrix[doublesRead++]));           

      if(result == EOF)
	{
	  printf("Error protein model file must consist of exactly 420 entries \n");
	  printf("The first 400 entries are for the rates of the AA matrix, while the\n");
	  printf("last 20 should contain the empirical base frequencies\n");
	  printf("Reached End of File after %d entries\n", (doublesRead - 1));
	  exit(-1);
	}    
    }
       
  fclose(f);

  /* CHECKS */
  for(i = 0; i < 20; i++)
    for(j = 0; j < 20; j++)
      {
	if(i != j)
	  {
	    if(adef->externalAAMatrix[i * 20 + j] != adef->externalAAMatrix[j * 20 + i])
	      {
		printf("Error user-defined Protein model matrix must be symmetric\n");
		printf("Entry P[%d][%d]=%f at position %d is not equal to P[%d][%d]=%f at position %d\n", 
		       i, j,  adef->externalAAMatrix[i * 20 + j], (i * 20 + j),
		       j, i,  adef->externalAAMatrix[j * 20 + i], (j * 20 + i));
		exit(-1);
	      }
	  }
      }

  acc = 0.0;

  for(i = 400; i < 420; i++)    
    acc += adef->externalAAMatrix[i];         

  if((acc > 1.0 + 1.0E-6) || (acc <  1.0 - 1.0E-6))
    {
      printf("Base frequencies in user-defined AA substitution matrix do not sum to 1.0\n");
      printf("the sum is %1.80f\n", acc);
      exit(-1);
    }

}




