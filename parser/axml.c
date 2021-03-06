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

#ifdef WIN32
#include <direct.h>
#endif

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
#include <stdarg.h>
#include <limits.h>
#include <inttypes.h>

#ifdef  _FINE_GRAIN_MPI
#include <mpi.h>
#endif



#ifdef _USE_PTHREADS
#include <pthread.h>

#endif

#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
#include <xmmintrin.h>
/*
  special bug fix, enforces denormalized numbers to be flushed to zero,
  without this program is a tiny bit faster though.
  #include <emmintrin.h> 
  #define MM_DAZ_MASK    0x0040
  #define MM_DAZ_ON    0x0040
  #define MM_DAZ_OFF    0x0000
*/
#endif

#include "axml.h"
#include "globalVariables.h"







/***************** UTILITY FUNCTIONS **************************/


void myBinFwrite(const void *ptr, size_t size, size_t nmemb)
{ 
  size_t  
    bytes_written = fwrite(ptr, size, nmemb, byteFile);
  
  assert(bytes_written == nmemb);
}





void *malloc_aligned(size_t size) 
{
  void 
    *ptr = (void *)NULL;
 
  int 
    res;
  

#if defined (__APPLE__)
  /* 
     presumably malloc on MACs always returns 
     a 16-byte aligned pointer
  */

  ptr = malloc(size);
  
  if(ptr == (void*)NULL) 
   assert(0);
  
#ifdef __AVX
  assert(0);
#endif


#else
  res = posix_memalign( &ptr, BYTE_ALIGNMENT, size );

  if(res != 0) 
    assert(0);
#endif 
   
  return ptr;
}








void printBothOpen(const char* format, ... )
{
  FILE *f = myfopen(infoFileName, "ab");

  va_list args;
  va_start(args, format);
  vfprintf(f, format, args );
  va_end(args);

  va_start(args, format);
  vprintf(format, args );
  va_end(args);

  fclose(f);
}

void printBothOpenMPI(const char* format, ... )
{
#ifdef _WAYNE_MPI
  if(processID == 0)
#endif
    {
      FILE *f = myfopen(infoFileName, "ab");

      va_list args;
      va_start(args, format);
      vfprintf(f, format, args );
      va_end(args);
      
      va_start(args, format);
      vprintf(format, args );
      va_end(args);
      
      fclose(f);
    }
}


boolean getSmoothFreqs(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].smoothFrequencies;
}

const unsigned int *getBitVector(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].bitVector;
}


int getStates(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].states;
}

unsigned char getUndetermined(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].undetermined;
}





partitionLengths *getPartitionLengths(pInfo *p)
{
  int 
    dataType  = p->dataType,
    states    = p->states,
    tipLength = p->maxTipStates;

  assert(states != -1 && tipLength != -1);

  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  pLength.leftLength = pLength.rightLength = states * states;
  pLength.eignLength = states;
  pLength.evLength   = states * states;
  pLength.eiLength   = states * states;
  pLength.substRatesLength = (states * states - states) / 2;
  pLength.frequenciesLength = states;
  pLength.tipVectorLength   = tipLength * states;
  pLength.symmetryVectorLength = (states * states - states) / 2;
  pLength.frequencyGroupingLength = states;  
  pLength.nonGTR = FALSE;  
  pLength.optimizeBaseFrequencies = FALSE;

  return (&pLengths[dataType]); 
}







double gettime(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 60.0*localtm.tm_min + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
#endif
}



double randum (long  *seed)
{
  long  sum, mult0, mult1, seed0, seed1, seed2, newseed0, newseed1, newseed2;
  double res;

  mult0 = 1549;
  seed0 = *seed & 4095;
  sum  = mult0 * seed0;
  newseed0 = sum & 4095;
  sum >>= 12;
  seed1 = (*seed >> 12) & 4095;
  mult1 =  406;
  sum += mult0 * seed1 + mult1 * seed0;
  newseed1 = sum & 4095;
  sum >>= 12;
  seed2 = (*seed >> 24) & 255;
  sum += mult0 * seed2 + mult1 * seed1;
  newseed2 = sum & 255;

  *seed = newseed2 << 24 | newseed1 << 12 | newseed0;
  res = 0.00390625 * (newseed2 + 0.000244140625 * (newseed1 + 0.000244140625 * newseed0));

  return res;
}

static int filexists(char *filename)
{
  FILE *fp;
  int res;
  fp = fopen(filename,"rb");

  if(fp)
    {
      res = 1;
      fclose(fp);
    }
  else
    res = 0;

  return res;
}


FILE *myfopen(const char *path, const char *mode)
{
  FILE *fp = fopen(path, mode);

  if(strcmp(mode,"r") == 0 || strcmp(mode,"rb") == 0)
    {
      if(fp)
	return fp;
      else
	{
	  if(processID == 0)
	    printf("\n Error: the file %s you want to open for reading does not exist, exiting ...\n\n", path);
	  errorExit(-1);
	  return (FILE *)NULL;
	}
    }
  else
    {
      if(fp)
	return fp;
      else
	{
	  if(processID == 0)
	    printf("\n Error: the file %s you want to open for writing or appending can not be opened [mode: %s], exiting ...\n\n",
		   path, mode);
	  errorExit(-1);
	  return (FILE *)NULL;
	}
    }


}





/********************* END UTILITY FUNCTIONS ********************/


/******************************some functions for the likelihood computation ****************************/










/***********************reading and initializing input ******************/

static void getnums (rawdata *rdta)
{
  if (fscanf(INFILE, "%d %" PRId64, & rdta->numsp, & rdta->sites) != 2)
    {
      if(processID == 0)
	printf("\n Error: problem reading number of species and sites\n\n");
      errorExit(-1);
    }

  if (rdta->numsp < 4)
    {
      if(processID == 0)
	printf("\n Error: too few species\n\n");
      errorExit(-1);
    }

  if (rdta->sites < 1)
    {
      if(processID == 0)
	printf("\n Error: too few sites\n\n");
      errorExit(-1);
    }

  return;
}





boolean whitechar (int ch)
{
  return (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r');
}


static void uppercase (int *chptr)
{
  int  ch;

  ch = *chptr;
  if ((ch >= 'a' && ch <= 'i') || (ch >= 'j' && ch <= 'r')
      || (ch >= 's' && ch <= 'z'))
    *chptr = ch + 'A' - 'a';
}




static void getyspace (rawdata *rdta)
{
  size_t size = 4 * ((size_t)(rdta->sites / 4 + 1));
  
 

  int    i;
  unsigned char *y0;

  rdta->y = (unsigned char **) malloc(((size_t)rdta->numsp + 1) * sizeof(unsigned char *));
  assert(rdta->y);   

  y0 = (unsigned char *)calloc(((size_t)(rdta->numsp + 1)) * size, sizeof(unsigned char));

  /*
    printf("Raw alignment data Assigning %Zu bytes\n", ((size_t)(rdta->numsp + 1)) * size * sizeof(unsigned char));

  */

  assert(y0);   

  rdta->y0 = y0;

  for (i = 0; i <= rdta->numsp; i++)
    {
      rdta->y[i] = y0;
      y0 += size;
    }

  return;
}




static boolean setupTree (tree *tr, analdef *adef)
{
  nodeptr  
    p0;
  
  int
    tips,
    inter; 

  if(!adef->readTaxaOnly)
    {
      /*tr->bigCutoff = FALSE;*/

      tr->patternPosition = (int64_t*)NULL;
      tr->columnPosition = (int64_t*)NULL;

      /*tr->maxCategories = MAX(4, adef->categories);*/

      /*tr->partitionContributions = (double *)malloc(sizeof(double) * tr->NumberOfModels);

      for(i = 0; i < tr->NumberOfModels; i++)
	tr->partitionContributions[i] = -1.0;

      tr->perPartitionLH = (double *)malloc(sizeof(double) * tr->NumberOfModels);
      

      for(i = 0; i < tr->NumberOfModels; i++)
	{
	  tr->perPartitionLH[i] = 0.0;	 
	}

      if(adef->grouping)
	tr->grouped = TRUE;
      else
	tr->grouped = FALSE;

      if(adef->constraint)
	tr->constrained = TRUE;
      else
	tr->constrained = FALSE;

	tr->treeID = 0;*/
    }

  tips  = tr->mxtips;
  inter = tr->mxtips - 1;

  if(!adef->readTaxaOnly)
    {
      tr->yVector      = (unsigned char **)  malloc(((size_t)tr->mxtips + 1) * sizeof(unsigned char *));

      /*      tr->fracchanges  = (double *)malloc(tr->NumberOfModels * sizeof(double));
	      tr->likelihoods  = (double *)malloc(adef->multipleRuns * sizeof(double));*/
    }

  /*tr->numberOfTrees = -1;

 

  tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;

  tr->tree_string  = (char*)calloc(tr->treeStringLength, sizeof(char)); 
  tr->tree0 = (char*)calloc(tr->treeStringLength, sizeof(char));
  tr->tree1 = (char*)calloc(tr->treeStringLength, sizeof(char));*/


  /*TODO, must that be so long ?*/

  if(!adef->readTaxaOnly)
    {
            
      /*tr->td[0].count = 0;
      tr->td[0].ti    = (traversalInfo *)malloc(sizeof(traversalInfo) * tr->mxtips);
      tr->td[0].executeModel = (boolean *)malloc(sizeof(boolean) * tr->NumberOfModels);
      tr->td[0].parameterValues = (double *)malloc(sizeof(double) * tr->NumberOfModels);
       
      for(i = 0; i < tr->NumberOfModels; i++)
	tr->fracchanges[i] = -1.0;
      tr->fracchange = -1.0;

      tr->constraintVector = (int *)malloc((2 * tr->mxtips) * sizeof(int));*/

      tr->nameList = (char **)malloc(sizeof(char *) * ((size_t)tips + 1));
    }

  if (!(p0 = (nodeptr) malloc(((size_t)tips + 3 * (size_t)inter) * sizeof(node))))
    {
      printf("\n Error: unable to obtain sufficient tree memory\n\n");
      return  FALSE;
    }
  
 

  

  tr->vLength = 0;

  tr->h = (hashtable*)NULL;


  return TRUE;
}


static void checkTaxonName(char *buffer, int len)
{
  int i;

  for(i = 0; i < len - 1; i++)
    {
      boolean valid;

      switch(buffer[i])
	{
	case '\0':
	case '\t':
	case '\n':
	case '\r':
	case ' ':
	case ':':
	case ',':
	case '(':
	case ')':
	case ';':
	case '[':
	case ']':
	  valid = FALSE;
	  break;
	default:
	  valid = TRUE;
	}

      if(!valid)
	{
	  printf("\n Error: Taxon Name \"%s\" is invalid at position %d, it contains illegal character %c\n\n", buffer, i, buffer[i]);
	  printf(" Illegal characters in taxon-names are: tabulators, carriage returns, spaces, \":\", \",\", \")\", \"(\", \";\", \"]\", \"[\"\n");
	  printf(" Exiting\n");
	  exit(-1);
	}

    }
  assert(buffer[len - 1] == '\0');
}

static boolean getdata(analdef *adef, rawdata *rdta, tree *tr)
{
  int64_t
    i, 
    j, 
    basesread, 
    basesnew;
   
  int
    ch, my_i, meaning,
    len,
    meaningAA[256], 
    meaningDNA[256], 
    meaningBINARY[256],
    meaningGeneric32[256],
    meaningGeneric64[256];
  
  boolean  
    allread, 
    firstpass;
  
  char 
    buffer[nmlngth + 2];
  
  unsigned char
    genericChars32[32] = {'0', '1', '2', '3', '4', '5', '6', '7', 
			  '8', '9', 'A', 'B', 'C', 'D', 'E', 'F',
			  'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
			  'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V'};  
  unsigned long 
    total = 0,
    gaps  = 0;

  for (i = 0; i < 256; i++)
    {      
      meaningAA[i]          = -1;
      meaningDNA[i]         = -1;
      meaningBINARY[i]      = -1;
      meaningGeneric32[i]   = -1;
      meaningGeneric64[i]   = -1;      
    }

  /* generic 32 data */

  for(i = 0; i < 32; i++)
    meaningGeneric32[genericChars32[i]] = (int)i;
  meaningGeneric32['-'] = getUndetermined(GENERIC_32);
  meaningGeneric32['?'] = getUndetermined(GENERIC_32);

  /* AA data */

  meaningAA['A'] =  0;  /* alanine */
  meaningAA['R'] =  1;  /* arginine */
  meaningAA['N'] =  2;  /*  asparagine*/
  meaningAA['D'] =  3;  /* aspartic */
  meaningAA['C'] =  4;  /* cysteine */
  meaningAA['Q'] =  5;  /* glutamine */
  meaningAA['E'] =  6;  /* glutamic */
  meaningAA['G'] =  7;  /* glycine */
  meaningAA['H'] =  8;  /* histidine */
  meaningAA['I'] =  9;  /* isoleucine */
  meaningAA['L'] =  10; /* leucine */
  meaningAA['K'] =  11; /* lysine */
  meaningAA['M'] =  12; /* methionine */
  meaningAA['F'] =  13; /* phenylalanine */
  meaningAA['P'] =  14; /* proline */
  meaningAA['S'] =  15; /* serine */
  meaningAA['T'] =  16; /* threonine */
  meaningAA['W'] =  17; /* tryptophan */
  meaningAA['Y'] =  18; /* tyrosine */
  meaningAA['V'] =  19; /* valine */
  meaningAA['B'] =  20; /* asparagine, aspartic 2 and 3*/
  meaningAA['Z'] =  21; /*21 glutamine glutamic 5 and 6*/

  meaningAA['X'] = 
    meaningAA['?'] = 
    meaningAA['*'] = 
    meaningAA['-'] = 
    getUndetermined(AA_DATA);

  /* DNA data */

  meaningDNA['A'] =  1;
  meaningDNA['B'] = 14;
  meaningDNA['C'] =  2;
  meaningDNA['D'] = 13;
  meaningDNA['G'] =  4;
  meaningDNA['H'] = 11;
  meaningDNA['K'] = 12;
  meaningDNA['M'] =  3;  
  meaningDNA['R'] =  5;
  meaningDNA['S'] =  6;
  meaningDNA['T'] =  8;
  meaningDNA['U'] =  8;
  meaningDNA['V'] =  7;
  meaningDNA['W'] =  9; 
  meaningDNA['Y'] = 10;

  meaningDNA['N'] = 
    meaningDNA['O'] = 
    meaningDNA['X'] = 
    meaningDNA['-'] = 
    meaningDNA['?'] = 
    getUndetermined(DNA_DATA);

  /* BINARY DATA */

  meaningBINARY['0'] = 1;
  meaningBINARY['1'] = 2;
  
  meaningBINARY['-'] = 
    meaningBINARY['?'] = 
    getUndetermined(BINARY_DATA);

 

  /*******************************************************************/

  basesread = basesnew = 0;

  allread = FALSE;
  firstpass = TRUE;
  ch = ' ';

  while (! allread)
    {
      for (i = 1; i <= tr->mxtips; i++)
	{
	  if (firstpass)
	    {
	      ch = getc(INFILE);
	      while(ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r')
		ch = getc(INFILE);

	      my_i = 0;

	      do
		{
		  buffer[my_i] = (char)ch;
		  ch = getc(INFILE);
		  my_i++;
		  if(my_i >= nmlngth)
		    {
		      if(processID == 0)
			{
			  printf("Taxon Name to long at taxon %" PRId64 ", adapt constant nmlngth in\n", i);
			  printf("axml.h, current setting %d\n", nmlngth);
			}
		      errorExit(-1);
		    }
		}
	      while(ch !=  ' ' && ch != '\n' && ch != '\t' && ch != '\r');

	      while(ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r')
		ch = getc(INFILE);
	      
	      ungetc(ch, INFILE);

	      buffer[my_i] = '\0';
	      len = (int)strlen(buffer) + 1;
	      checkTaxonName(buffer, len);
	      tr->nameList[i] = (char *)malloc(sizeof(char) * (size_t)len);
	      strcpy(tr->nameList[i], buffer);
	    }

	  j = basesread;

	  while ((j < rdta->sites) && ((ch = getc(INFILE)) != EOF) && (ch != '\n') && (ch != '\r'))
	    {
	      uppercase(& ch);

	      assert(tr->dataVector[j + 1] != -1);

	      switch(tr->dataVector[j + 1])
		{
		case BINARY_DATA:
		  meaning = meaningBINARY[ch];
		  break;
		case DNA_DATA:
		case SECONDARY_DATA:
		case SECONDARY_DATA_6:
		case SECONDARY_DATA_7:
		  /*
		     still dealing with DNA/RNA here, hence just act if as they where DNA characters
		     corresponding column merging for sec struct models will take place later
		  */
		  meaning = meaningDNA[ch];
		  break;
		case AA_DATA:
		  meaning = meaningAA[ch];
		  break;
		case GENERIC_32:
		  meaning = meaningGeneric32[ch];
		  break;
		case GENERIC_64:
		  meaning = meaningGeneric64[ch];
		  break;
		  //mth now map the character to the state
		case POMO_16:
		case POMO_64:
		  meaning = meaningDNA[ch];		
		  break;
		default:
		  assert(0);
		}

	      if (meaning != -1)
		{
		  j++;
		  rdta->y[i][j] = (unsigned char)ch;		 
		}
	      else
		{
		  if(!whitechar(ch))
		    {
		      printf("\n Error: bad base (%c) at site %" PRId64 " of sequence %" PRId64 "\n\n",
			     ch, j + 1, i);
		      return FALSE;
		    }
		}
	    }

	  if (ch == EOF)
	    {
	      printf("\n Error: end-of-file at site %"  PRId64 " of sequence %" PRId64 "\n\n", j + 1, i);
	      return  FALSE;
	    }

	  if (! firstpass && (j == basesread))
	    i--;
	  else
	    {
	      if (i == 1)
		basesnew = j;
	      else
		if (j != basesnew)
		  {
		    printf("\n Error: sequences out of alignment\n");
		    printf("%"  PRId64 " (instead of %"  PRId64 ") residues read in sequence %"  PRId64  " %s\n",
			   j - basesread, basesnew - basesread, i, tr->nameList[i]);
		    return  FALSE;
		  }
	    }
	  while (ch != '\n' && ch != EOF && ch != '\r') ch = getc(INFILE);  /* flush line *//* PC-LINEBREAK*/
	}

      firstpass = FALSE;
      basesread = basesnew;
      allread = (basesread >= rdta->sites);
    }

  for(j = 1; j <= tr->mxtips; j++)
    for(i = 1; i <= rdta->sites; i++)
      {
	assert(tr->dataVector[i] != -1);

	switch(tr->dataVector[i])
	  {
	  case BINARY_DATA:
	    meaning = meaningBINARY[rdta->y[j][i]];
	    if(meaning == getUndetermined(BINARY_DATA))
	      gaps++;
	    break;

	  case SECONDARY_DATA:
	  case SECONDARY_DATA_6:
	  case SECONDARY_DATA_7:
	    assert(tr->secondaryStructurePairs[i - 1] != -1);
	    assert(i - 1 == tr->secondaryStructurePairs[tr->secondaryStructurePairs[i - 1]]);
	    /*
	       don't worry too much about undetermined column count here for sec-struct, just count
	       DNA/RNA gaps here and worry about the rest later-on, falling through to DNA again :-)
	    */
	  case DNA_DATA:
	    meaning = meaningDNA[rdta->y[j][i]];
	    if(meaning == getUndetermined(DNA_DATA))
	      gaps++;
	    break;

	  case AA_DATA:
	    meaning = meaningAA[rdta->y[j][i]];
	    if(meaning == getUndetermined(AA_DATA))
	      gaps++;
	    break;

	  case GENERIC_32:
	    meaning = meaningGeneric32[rdta->y[j][i]];
	    if(meaning == getUndetermined(GENERIC_32))
	      gaps++;
	    break;

	  case GENERIC_64:
	    meaning = meaningGeneric64[rdta->y[j][i]];
	    if(meaning == getUndetermined(GENERIC_64))
	      gaps++;
	    break;
	    //mth count gaps in alignment
	  case POMO_16:
	  case POMO_64:
	    meaning = meaningDNA[rdta->y[j][i]];
	    if(meaning == getUndetermined(DNA_DATA))
	      gaps++;	   
	    break;
	  default:
	    assert(0);
	  }

	total++;
	rdta->y[j][i] = (unsigned char)meaning;
      }

  adef->gapyness = (double)gaps / (double)total;
    
  /*myBinFwrite(&(adef->gapyness), sizeof(double), 1);*/

  printf("\n\ngappyness: %f\n", adef->gapyness);
  
  /*for(i = 1; i <= tr->mxtips; i++)
    {
      int 
	len = strlen(tr->nameList[i]) + 1;
      
      myBinFwrite(&len, sizeof(int), 1);
      myBinFwrite(tr->nameList[i], sizeof(char), len);
      
      printf("%d %s\n", len, tr->nameList[i]);
      }     */
  
  return  TRUE;
}




static hashNumberType hashString(char *p, hashNumberType tableSize)
{
  hashNumberType 
    h = 0;
  
  for(; *p; p++)
    h = 31 * h + (hashNumberType)*p;
  
  return (h % tableSize);
}

static void addword(char *s, stringHashtable *h, int nodeNumber)
{
  hashNumberType position = hashString(s, h->tableSize);
  stringEntry *p = h->table[position];
  
  for(; p!= NULL; p = p->next)
    {
      if(strcmp(s, p->word) == 0)		 
	return;	  	
    }

  p = (stringEntry *)malloc(sizeof(stringEntry));

  assert(p);
  
  p->nodeNumber = nodeNumber;
  p->word = (char *)malloc((strlen(s) + 1) * sizeof(char));

  strcpy(p->word, s);
  
  p->next =  h->table[position];
  
  h->table[position] = p;
}


static stringHashtable *initStringHashTable(hashNumberType n)
{
  /* 
     init with primes 
  */
    
  static const hashNumberType initTable[] = {53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
					     196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
					     50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};
 

  /* init with powers of two

  static const  hashNumberType initTable[] = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,
					      32768, 65536, 131072, 262144, 524288, 1048576, 2097152,
					      4194304, 8388608, 16777216, 33554432, 67108864, 134217728,
					      268435456, 536870912, 1073741824, 2147483648U};
  */
  
  stringHashtable *h = (stringHashtable*)malloc(sizeof(stringHashtable));
  
  hashNumberType
    tableSize,
    i,
    primeTableLength = sizeof(initTable)/sizeof(initTable[0]),
    maxSize = (hashNumberType)-1;    

  assert(n <= maxSize);

  i = 0;

  while(initTable[i] < n && i < primeTableLength)
    i++;

  assert(i < primeTableLength);

  tableSize = initTable[i];  

  h->table = (stringEntry**)calloc(tableSize, sizeof(stringEntry*));
  h->tableSize = tableSize;    

  return h;
}



static void getinput(analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr)
{
  int64_t 
    i;

  INFILE = myfopen(seq_file, "rb");
  
  getnums(rdta);
  
     
  /*myBinFwrite(&(rdta->sites), sizeof(int), 1);
  myBinFwrite(&(rdta->numsp), sizeof(int), 1);  

  printf("%d %d\n", rdta->sites, rdta->numsp);*/
    

  tr->mxtips            = rdta->numsp;
  
  
  rdta->wgt             = (int *)    malloc(((size_t)rdta->sites + 1) * sizeof(int));
  cdta->alias           = (int64_t *)    malloc(((size_t)rdta->sites + 1) * sizeof(int64_t));
  cdta->aliaswgt        = (int *)    malloc(((size_t)rdta->sites + 1) * sizeof(int)); 
  tr->model             = (int *)    calloc(((size_t)rdta->sites + 1), sizeof(int));
  tr->initialDataVector  = (int *)    malloc(((size_t)rdta->sites + 1) * sizeof(int));
  tr->extendedDataVector = (int *)    malloc(((size_t)rdta->sites + 1) * sizeof(int));         
  
  assert(!adef->useWeightFile);
    
  for (i = 1; i <= rdta->sites; i++)
    rdta->wgt[i] = 1;
  
  if(adef->useMultipleModel)
    {
      int ref;
      
      parsePartitions(adef, rdta, tr);
      
      for(i = 1; i <= rdta->sites; i++)
	{
	  ref = tr->model[i];
	  tr->initialDataVector[i] = tr->initialPartitionData[ref].dataType;
	}
    }
  else
    {
      int dataType = -1;
	              
      tr->initialPartitionData  = (pInfo*)malloc(sizeof(pInfo));
      tr->initialPartitionData->optimizeBaseFrequencies = FALSE;
      
      
      tr->initialPartitionData[0].partitionName = (char*)malloc(128 * sizeof(char));
      strcpy(tr->initialPartitionData[0].partitionName, "No Name Provided");
      
      tr->initialPartitionData[0].protModels = adef->proteinMatrix;
      tr->initialPartitionData[0].protFreqs  = adef->protEmpiricalFreqs;
      
      
      tr->NumberOfModels = 1;
           
      
      if(adef->model == M_PROTCAT || adef->model == M_PROTGAMMA)
	dataType = AA_DATA;
      if(adef->model == M_GTRCAT || adef->model == M_GTRGAMMA)
	dataType = DNA_DATA;
      if(adef->model == M_BINCAT || adef->model == M_BINGAMMA)
	dataType = BINARY_DATA;
      if(adef->model == M_32CAT || adef->model == M_32GAMMA)
	dataType = GENERIC_32;
      if(adef->model == M_64CAT || adef->model == M_64GAMMA)
	dataType = GENERIC_64;
      //mth set data type to POMO
      if(adef->model == M_POMOGAMMA_16)
	dataType = POMO_16;
       if(adef->model == M_POMOGAMMA_64)
	 dataType = POMO_64;
      
      
      assert(dataType == BINARY_DATA || dataType == DNA_DATA || dataType == AA_DATA || 
	     dataType == GENERIC_32  || dataType == GENERIC_64 || dataType == POMO_16 || 
	     dataType == POMO_64);
      
      tr->initialPartitionData[0].dataType = dataType;
      
      for(i = 0; i <= rdta->sites; i++)
	{
	  tr->initialDataVector[i] = dataType;
	  tr->model[i]      = 0;
	}
    }
  
 
  tr->dataVector    = tr->initialDataVector;
  tr->partitionData = tr->initialPartitionData;
 
  
 
  
  getyspace(rdta);

  setupTree(tr, adef);

      
	
  if(!getdata(adef, rdta, tr))
    {
      printf("Problem reading alignment file \n");
      errorExit(1);
    }
      
  tr->nameHash = initStringHashTable(10 * (unsigned int)tr->mxtips);
  tr->speciesHash = initStringHashTable(10 * (unsigned int)tr->mxtips);

  for(i = 1; i <= tr->mxtips; i++)
    addword(tr->nameList[i], tr->nameHash, (int)i);
      
  fclose(INFILE);
}







static void sitesort(rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)
{
  int64_t  
    gap, i, j, jj, jg, k, n, nsp;
  
  int64_t  
    *index; 
    
  int
    *category = (int*)NULL;

  boolean  
    flip, 
    tied;
  
  unsigned char  
    **data;

  if(adef->useMultipleModel)    
    category      = tr->model;
  

  index    = cdta->alias;
  data     = rdta->y;
  n        = rdta->sites;
  nsp      = rdta->numsp;
  index[0] = -1;


  if(adef->compressPatterns)
    {
      for (gap = n / 2; gap > 0; gap /= 2)
	{
	  for (i = gap + 1; i <= n; i++)
	    {
	      j = i - gap;

	      do
		{
		  jj = index[j];
		  jg = index[j+gap];
		  if(adef->useMultipleModel)
		    {		     		      
		      assert(category[jj] != -1 &&
			     category[jg] != -1);
		     
		      flip = (category[jj] > category[jg]);
		      tied = (category[jj] == category[jg]);		     

		    }
		  else
		    {
		      flip = 0;
		      tied = 1;
		    }

		  for (k = 1; (k <= nsp) && tied; k++)
		    {
		      flip = (data[k][jj] >  data[k][jg]);
		      tied = (data[k][jj] == data[k][jg]);
		    }

		  if (flip)
		    {
		      index[j]     = jg;
		      index[j+gap] = jj;
		      j -= gap;
		    }
		}
	      while (flip && (j > 0));
	    }
	}
    }
}


static void sitecombcrunch (rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)
{
  
  boolean  
    tied;
  
  int64_t
    i,
    sitei, 
    j, 
    sitej, 
    k;

  int
    *aliasModel = (int*)NULL,
    *aliasSuperModel = (int*)NULL;

  int64_t
    undeterminedSites = 0;

  if(adef->useMultipleModel)
    {
      aliasSuperModel = (int*)malloc(sizeof(int) * ((size_t)rdta->sites + 1));
      aliasModel      = (int*)malloc(sizeof(int) * ((size_t)rdta->sites + 1));
    } 

  i = 0;
  cdta->alias[0]    = cdta->alias[1];
  cdta->aliaswgt[0] = 0;

  if(adef->mode == PER_SITE_LL)
    {      
      assert(0);

      /*
      tr->patternPosition = (int*)malloc(sizeof(int) * rdta->sites);
      tr->columnPosition  = (int*)malloc(sizeof(int) * rdta->sites);

      for(i = 0; i < rdta->sites; i++)
	{
	  tr->patternPosition[i] = -1;
	  tr->columnPosition[i]  = -1;
	}
      */
    }

  

  i = 0;
  for (j = 1; j <= rdta->sites; j++)
    {
      int 
	allGap = TRUE;

      unsigned char 
	undetermined;

      sitei = cdta->alias[i];
      sitej = cdta->alias[j];      

      undetermined = getUndetermined(tr->dataVector[sitej]);
      
      for(k = 1; k <= rdta->numsp; k++)
	{	 
	  if(rdta->y[k][sitej] != undetermined)
	    {
	      allGap = FALSE;
	      break;
	    }
	}
      
      if(allGap)      
      	undeterminedSites++;	 

#ifdef _DEBUG_UNDET_REMOVAL
      if(allGap)
	printf("Skipping gap site %d\n", sitej);
#endif
          
      if(!adef->compressPatterns)
	tied = 0;
      else
	{
	  if(adef->useMultipleModel)
	    {	     
	      tied = (tr->model[sitei] == tr->model[sitej]);
	      if(tied)
		assert(tr->dataVector[sitei] == tr->dataVector[sitej]);
	    }
	  else
	    tied = 1;
	}
      
      for (k = 1; tied && (k <= rdta->numsp); k++)
	tied = (rdta->y[k][sitei] == rdta->y[k][sitej]);
	      
      assert(!(tied && allGap));

      if(tied && !allGap)
	{
	  if(adef->mode == PER_SITE_LL)
	    {
	      tr->patternPosition[j - 1] = i;
	      tr->columnPosition[j - 1] = sitej;
	      /*printf("Pattern %d from column %d also at site %d\n", i, sitei, sitej);*/
	    }


	  cdta->aliaswgt[i] += rdta->wgt[sitej];
	  if(adef->useMultipleModel)
	    {
	      aliasModel[i]      = tr->model[sitej];
	      aliasSuperModel[i] = tr->dataVector[sitej];
	    }
	}
      else
	{
	  if(!allGap)
	    {
	      if(cdta->aliaswgt[i] > 0) 
		i++;
	      
	      if(adef->mode == PER_SITE_LL)
		{
		  tr->patternPosition[j - 1] = i;
		  tr->columnPosition[j - 1] = sitej;
		  /*printf("Pattern %d is from cloumn %d\n", i, sitej);*/
		}
	      
	      cdta->aliaswgt[i] = rdta->wgt[sitej];
	      cdta->alias[i] = sitej;
	      if(adef->useMultipleModel)
		{
		  aliasModel[i]      = tr->model[sitej];
		  aliasSuperModel[i] = tr->dataVector[sitej];
		}
	    }	  
	}
    }

  cdta->endsite = (size_t)i;
  if (cdta->aliaswgt[i] > 0) 
    cdta->endsite++;

#ifdef _DEBUG_UNDET_REMOVAL
  printf("included sites: %d\n", cdta->endsite);
#endif

  if(adef->mode == PER_SITE_LL)
    {
      assert(0);

      for(i = 0; i < rdta->sites; i++)
	{
	  int64_t p  = tr->patternPosition[i];
	  int64_t c  = tr->columnPosition[i];

	  assert(p >= 0 && (size_t) p < cdta->endsite);
	  assert(c >= 1 && c <= rdta->sites);
	}
    }


  if(adef->useMultipleModel)
    {
      for(i = 0; i <= rdta->sites; i++)
	{
	  tr->model[i]      = aliasModel[i];
	  tr->dataVector[i] = aliasSuperModel[i];	  
	}
    }

  if(adef->useMultipleModel)
    {
      free(aliasModel);
      free(aliasSuperModel);
    }     

  if(undeterminedSites > 0)    
    printBothOpen("\nAlignment has %d completely undetermined sites that will be automatically removed from the binary alignment file\n\n", undeterminedSites);
}


static boolean makeweights (analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr)
{
  int  i;

 
    
  for (i = 1; i <= rdta->sites; i++)
    cdta->alias[i] = i;

  sitesort(rdta, cdta, tr, adef);
  sitecombcrunch(rdta, cdta, tr, adef);  

  return TRUE;
}



static boolean makevalues(rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)
{
  int  
    i, 
    model, 
    modelCounter;
  
  size_t 
    j;

  unsigned char
    *y    = (unsigned char *)malloc(((size_t)rdta->numsp) * ((size_t)cdta->endsite) * sizeof(unsigned char));
  

  /*

  printf("compressed data Assigning %Zu bytes\n", ((size_t)rdta->numsp) * ((size_t)cdta->endsite) * sizeof(unsigned char));

  */
  
  
    {
      for (i = 1; i <= rdta->numsp; i++)
	for (j = 0; j < cdta->endsite; j++)   
	  y[(((size_t)(i - 1)) * ((size_t)cdta->endsite)) + j] = rdta->y[i][cdta->alias[j]];
      
      /*
	printf("Free on raw data\n");
      */

      free(rdta->y0);
      free(rdta->y);
      
    }

  rdta->y0 = y;
 
  if(!adef->useMultipleModel)
    tr->NumberOfModels = 1;

#ifdef _DEBUG_UNDET_REMOVAL
  for(i = 0; i < cdta->endsite; i++)
    printf("%d ", tr->model[i]);

  printf("\n");
#endif

  if(adef->useMultipleModel)
    {
      tr->partitionData[0].lower = 0;

      model        = tr->model[0];
      modelCounter = 0;
     
      i            = 1;

      while((size_t) i <  cdta->endsite)
	{
	  if(tr->model[i] != model)
	    {
	      tr->partitionData[modelCounter].upper     = (size_t)i;
	      tr->partitionData[modelCounter + 1].lower = (size_t)i;

	      model = tr->model[i];	     
	      modelCounter++;
	    }
	  i++;
	}

      if(modelCounter <  tr->NumberOfModels - 1)
	{
	  printf("\nYou specified %d partitions, but after parsing and pre-processing ExaML only found %d partitions\n", tr->NumberOfModels, modelCounter + 1);
	  printf("Presumably one or more partitions vanished because they consisted entirely of undetermined characters.\n");
	  printf("Please fix your data!\n\n");
	  exit(-1);
	}


      tr->partitionData[tr->NumberOfModels - 1].upper = (size_t)cdta->endsite;      
    
      for(i = 0; i < tr->NumberOfModels; i++)		  
	tr->partitionData[i].width      = tr->partitionData[i].upper -  tr->partitionData[i].lower;
	 
      model        = tr->model[0];
      modelCounter = 0;
      tr->model[0] = modelCounter;
      i            = 1;
	
      while((size_t) i < cdta->endsite)
	{	 
	  if(tr->model[i] != model)
	    {
	      model = tr->model[i];
	      modelCounter++;
	      tr->model[i] = modelCounter;
	    }
	  else
	    tr->model[i] = modelCounter;
	  i++;
	}      
    }
  else
    {
      tr->partitionData[0].lower = 0;
      tr->partitionData[0].upper = (size_t)cdta->endsite;
      tr->partitionData[0].width =  tr->partitionData[0].upper -  tr->partitionData[0].lower;
    }

  tr->rdta       = rdta;
  tr->cdta       = cdta; 

  tr->originalCrunchedLength = tr->cdta->endsite;
    
  for(i = 0; i < rdta->numsp; i++)
    tr->yVector[i + 1] = &(rdta->y0[(tr->originalCrunchedLength) * ((size_t)i)]);

  return TRUE;
}



static void initAdef(analdef *adef)
{  
  adef->useSecondaryStructure  = FALSE;
  adef->bootstrapBranchLengths = FALSE;
  adef->model                  = M_GTRCAT;
  adef->max_rearrange          = 21;
  adef->stepwidth              = 5;
  adef->initial                = adef->bestTrav = 10;
  adef->initialSet             = FALSE;
  adef->restart                = FALSE;
  adef->mode                   = BIG_RAPID_MODE;
  adef->categories             = 25;
  adef->boot                   = 0;
  adef->rapidBoot              = 0;
  adef->useWeightFile          = FALSE;
  adef->checkpoints            = 0;
  adef->startingTreeOnly       = 0;
  adef->multipleRuns           = 1;
  adef->useMultipleModel       = FALSE;
  adef->likelihoodEpsilon      = 0.1;
  adef->constraint             = FALSE;
  adef->grouping               = FALSE;
  adef->randomStartingTree     = FALSE;
  adef->parsimonySeed          = 0;
  adef->proteinMatrix          = JTT;
  adef->protEmpiricalFreqs     = 0;  
  adef->useInvariant           = FALSE;
  adef->permuteTreeoptimize    = FALSE;
  adef->useInvariant           = FALSE;
  adef->allInOne               = FALSE;
  adef->likelihoodTest         = FALSE;
  adef->perGeneBranchLengths   = FALSE;
  adef->generateBS             = FALSE;
  adef->bootStopping           = FALSE;
  adef->gapyness               = 0.0;
  adef->similarityFilterMode   = 0;
  adef->useExcludeFile         = FALSE;
  adef->userProteinModel       = FALSE;
  adef->externalAAMatrix       = (double*)NULL;
  adef->computeELW             = FALSE;
  adef->computeDistance        = FALSE;
  adef->thoroughInsertion      = FALSE;
  adef->compressPatterns       = TRUE; 
  adef->readTaxaOnly           = FALSE;
  adef->meshSearch             = 0;
  adef->useCheckpoint          = FALSE;
  adef->leaveDropMode          = FALSE;
  adef->slidingWindowSize      = 100;
#ifdef _BAYESIAN 
  adef->bayesian               = FALSE;
#endif

}




static int dataExists(char *model, analdef *adef)
{
  /********** BINARY ********************/

   if(strcmp(model, "BIN\0") == 0)
    {
      adef->model = M_BINGAMMA;      
      return 1;
    }  

  /*********** DNA **********************/

  if(strcmp(model, "DNA\0") == 0)
    {
      adef->model = M_GTRGAMMA;     
      return 1;
    }

  /*************** AA GTR ********************/

  if(strcmp(model, "PROT\0") == 0)
    {
      adef->model = M_PROTGAMMA;     
      return 1;
    } 

  /************* POMO16 **********************/

  //mth allow users to specifi -m POMO16 as data model
  if(strcmp(model, "POMO16\0") == 0)
    {
      adef->model = M_POMOGAMMA_16;     
      return 1;
    } 
  
  /************** POMO64 ****************/

  //mth allow users to specifi -m POMO16 as data model
  if(strcmp(model, "POMO64\0") == 0)
    {
      adef->model = M_POMOGAMMA_64;     
      return 1;
    } 


  return 0;
}

/*********************************************************************************************/

static void printVersionInfo(void)
{
  printf("\n\nThis is the parse-examl version %s released by Alexandros Stamatakis, Andre J. Aberer, and Alexey Kozlov in %s.\n\n",  programVersion, programDate); 
}

static void printREADME(void)
{
  printVersionInfo();
  printf("\n");  
  printf("\nTo report bugs use the RAxML google group\n");
  printf("Please send us all input files, the exact invocation, details of the HW and operating system,\n");
  printf("as well as all error messages printed to screen.\n\n\n");

  printf("parse-examl\n");
  printf("      -s sequenceFileName\n");
  printf("      -n outputFileName\n");
  printf("      -m substitutionModel\n");
  printf("      -p pomoMapFile\n");
  printf("      [-c]\n");
  printf("      [-q]\n");
  printf("      [-h]\n");
  printf("\n"); 
  printf("      -m type of data to be parsed:\n");
  printf("\n"); 
  printf("              For Binary data use: BIN\n");
  printf("              For DNA data use:    DNA\n");	
  printf("              For AA data use:     PROT\n");			   
  printf("              For POMO data use:   POMO16 or POMO64\n");
  printf("\n"); 
  printf("      -p      Specify the name of the POMO species name to taxon names mapping of corresponding individuals.\n");
  printf("              The mapping file needs to be a plain text file containing one line per species.\n");
  printf("              Each species line needs to contain the species name followed by the taxon names of the corresponding\n");
  printf("              individuals from the DNA input alignment separated by whitespaces.\n");
  printf("\n");
  printf("      -c      disable site pattern compression\n");
  printf("\n");
  printf("      -q      Specify the file name which contains the assignment of models to alignment\n");
  printf("              partitions for multiple models of substitution. For the syntax of this file\n");
  printf("              please consult the manual.\n");  
  printf("\n");
  printf("      -h      Display this help message.\n");
  printf("\n");
  printf("\n\n\n\n");

}



/*********************************************************************************************/









static void analyzeRunId(char id[128])
{
  int i = 0;

  while(id[i] != '\0')
    {    
      if(i >= 128)
	{
	  printf("\n Error: run id after \"-n\" is too long, it has %d characters please use a shorter one\n\n", i);
	  assert(0);
	}
      
      if(id[i] == '/')
	{
	  printf("\n Error character %c not allowed in run ID\n\n", id[i]);
	  assert(0);
	}


      i++;
    }

  if(i == 0)
    {
      printf("\n Error: please provide a string for the run id after \"-n\" \n\n");
      assert(0);
    }

}


static void get_args(int argc, char *argv[], analdef *adef, tree *tr)
{
  boolean
    pomoMapSet = FALSE,
    bad_opt    = FALSE;

  char   
    model[2048] = "";
  
  int  
    c,
    nameSet = 0,
    alignmentSet = 0,     
    modelSet = 0;


  run_id[0] = 0; 
  seq_file[0] = 0;
  model[0] = 0;
  weightFileName[0] = 0;
  modelFileName[0] = 0;

  /*********** tr inits **************/
 
  tr->bootStopCriterion = -1;
  tr->wcThreshold = 0.03;
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = FALSE;
  tr->catOnly = FALSE;
 
  tr->multiStateModel  = GTR_MULTI_STATE;
  tr->useGappedImplementation = FALSE;
  tr->saveMemory = FALSE;
    
  /********* tr inits end*************/


  while( !bad_opt && ( ( c = getopt(argc,argv,"q:s:n:m:p:hc") ) != -1 ) )
    {
    switch(c)
      {                
      case 'c':
	adef->compressPatterns = FALSE;
	break;
      case 'h':
        printREADME();
	errorExit(0);
	break;
      case 'q':
	strcpy(modelFileName,optarg);
	adef->useMultipleModel = TRUE;
        break;                 
      case 'n':
        strcpy(run_id,optarg);
	analyzeRunId(run_id);
	nameSet = 1;
        break;     
      case 's':
	strcpy(seq_file, optarg);
	alignmentSet = 1;
	break;
      case 'm':
	strcpy(model,optarg);
	if(dataExists(model, adef) == 0)
	  {
	    printf("\n Error: model %s does not exist\n\n", model);               
	    errorExit(-1);
	  }
	else
	  modelSet = 1;
	break; 
      case 'p':
	strcpy(pomoMapFileName, optarg);
	pomoMapSet = TRUE;
	break;
	/*case 'r':
	adef->randomSeed;
	break;
      case 'b':
	adef->bootstrapReplicates;
	break;*/
      default:
	errorExit(-1);
    }
  }  

  if((adef->model == M_POMOGAMMA_16 || adef->model == M_POMOGAMMA_64) && !pomoMapSet)
    {
      if(processID == 0)
        {
          printREADME();	    
	  printf("\n Error, for the POMO models you need to specify a species mapping with the \"-p\" option\n\n");
        }
      errorExit(-1);
    }
    
  if(!adef->useMultipleModel && !modelSet)
    {
      if(processID == 0)
        {
          printREADME();	    
	  printf("\n Error, you must specify a data type for unpartitioned alignment with the \"-m\" option\n\n");
        }
      errorExit(-1);
    }

  if(!nameSet)
    {
      if(processID == 0)
        {
          printREADME();	    
	  printf("\n Error: please specify a name for this run with -n\n\n");
        }
      errorExit(-1);
    }


  if(!alignmentSet)
    {
      if(processID == 0)
        {
          printREADME();	    
	  printf("\n Error: please specify an alignment for this run with -s\n\n");
        } 
      errorExit(-1);
    }
  
  
   strcat(infoFileName,         "RAxML_info."); 
   strcat(infoFileName,         run_id);
  
   if(processID == 0)
     {
       int infoFileExists = 0;
       
       infoFileExists = filexists(infoFileName);
       
       if(infoFileExists)
	 {
	   printf("\n Error: output files with the run ID <%s> already exist... exiting\n\n", run_id);
	   exit(-1);
	 }
     }

  strcat(byteFileName, run_id);
  strcat(byteFileName, ".binary");
  
  if(filexists(byteFileName))
    {
      printf("\n Error: binary compressed file %s you want to generate already exists... exiting\n\n", byteFileName);
      exit(0);
    }

  byteFile = fopen(byteFileName, "wb");

  if ( !byteFile )  
    printf("%s\n", byteFileName);

  return;
}




void errorExit(int e)
{

#ifdef _WAYNE_MPI
  MPI_Finalize();
#endif

  exit(e);

}





 




/***********************reading and initializing input ******************/


/********************PRINTING various INFO **************************************/





void getDataTypeString(tree *tr, int model, char typeOfData[1024])
{
  switch(tr->partitionData[model].dataType)
    {
    case AA_DATA:
      strcpy(typeOfData,"AA");
      break;
    case DNA_DATA:
      strcpy(typeOfData,"DNA");
      break;
    case BINARY_DATA:
      strcpy(typeOfData,"BINARY/MORPHOLOGICAL");
      break;   
    case GENERIC_32:
      strcpy(typeOfData,"Multi-State");
      break;
    case GENERIC_64:
      strcpy(typeOfData,"Codon"); 
      break;
      //mth copy POMO to a string
    case POMO_16:
      strcpy(typeOfData, "POMO_16");
      break;
    case POMO_64:
      strcpy(typeOfData, "POMO_64");
      break;
    default:
      assert(0);
    }
}





/************************************************************************************/





  









static void smoothFreqs(const int n, double *pfreqs, double *dst, pInfo *partitionData)
{
  int 
    countScale = 0, 
    l,
    loopCounter = 0;  
  

 
  for(l = 0; l < n; l++)
    if(pfreqs[l] < FREQ_MIN)
      countScale++;
 

  /* for(l = 0; l < n; l++)
    if(pfreqs[l] == 0.0)
    countScale++;*/

  if(countScale > 0)
    {	     
      while(countScale > 0)
	{
	  double correction = 0.0;
	  double factor = 1.0;
	  
	  for(l = 0; l < n; l++)
	    {
	      if(pfreqs[l] == 0.0)		  
		correction += FREQ_MIN;		   		  
	      else
		if(pfreqs[l] < FREQ_MIN)		    
		  {
		    correction += (FREQ_MIN - pfreqs[l]);
		    factor -= (FREQ_MIN - pfreqs[l]);
		  }
	    }		      	    	    
	  
	  countScale = 0;
	  
	  for(l = 0; l < n; l++)
	    {		    
	      if(pfreqs[l] >= FREQ_MIN)		      
		pfreqs[l] = pfreqs[l] - (pfreqs[l] * correction * factor);	
	      else
		pfreqs[l] = FREQ_MIN;
	      
	      if(pfreqs[l] < FREQ_MIN)
		countScale++;
	    }
	  assert(loopCounter < 100);
	  loopCounter++;
	}		    
    }

  for(l = 0; l < n; l++)
    dst[l] = pfreqs[l];

  
  if(partitionData->nonGTR)
    {
      int k;

      assert(partitionData->dataType == SECONDARY_DATA_7 || partitionData->dataType == SECONDARY_DATA_6 || partitionData->dataType == SECONDARY_DATA);
       
      for(l = 0; l < n; l++)
	{
	  int count = 1;	
	  
	  for(k = 0; k < n; k++)
	    {
	      if(k != l && partitionData->frequencyGrouping[l] == partitionData->frequencyGrouping[k])
		{
		  count++;
		  dst[l] += pfreqs[k];
		}
	    }
	  dst[l] /= ((double)count);
	}            
     }  
}
	    


static void genericBaseFrequencies(tree *tr, const int numFreqs, rawdata *rdta, cruncheddata *cdta, size_t lower, size_t upper, int model, boolean smoothFrequencies,
				   const unsigned int *bitMask)
{
  double 
    wj, 
    acc,
    pfreqs[64], 
    sumf[64],   
    temp[64];
 
  size_t
    j,
    i;

  int     
    countStatesPresent = 0,
    statesPresent[64], 
    k, 
    l;

  unsigned char  
    *yptr;  
	  
  for(l = 0; l < numFreqs; l++)	    
    {
      pfreqs[l] = 1.0 / ((double)numFreqs);     
      statesPresent[l] = 0;
    }

#ifdef _DEBUG_UNDET_REMOVAL	  
  printf("bounds %d %d\n", lower, upper);

  for(j = lower; j < upper; j++) 
    {
      for(i = 0; i < rdta->numsp; i++)
	{
	  unsigned int 
	    code;

	  yptr = &(rdta->y0[((size_t)i) * (tr->originalCrunchedLength)]);
	  
	  code = yptr[j];

	  printf("%c",  inverseMeaningPROT[code]);
	}
      printf("\n");
    }
  
  printf("\n\n");
#endif

  for(i = 0; i < (size_t)rdta->numsp; i++) 
    {
      yptr = &(rdta->y0[((size_t)i) * (tr->originalCrunchedLength)]);
      
      for(j = lower; j < upper; j++) 
	{
	  unsigned int 	      
	    code = bitMask[yptr[j]];
	  
	  switch(numFreqs)
	    {
	    case 2:
	      switch(code)
		{
		case 1:
		  statesPresent[0] = 1;
		  break;
		case 2:
		  statesPresent[1] = 1;
		  break;
		default:
		  ;
		}
	      break;
	    case 4:	  
	      switch(code)
		{
		case 1:
		  statesPresent[0] = 1;
		  break;
		case 2:
		  statesPresent[1] = 1;
		  break;
		case 4:
		  statesPresent[2] = 1;
		  break;
		case 8:
		  statesPresent[3] = 1;
		  break;
		default:
		  ;
		}
	      break;	       
	    case 20:	      
	      if(yptr[j] >= 0 && yptr[j] < 20)
		statesPresent[yptr[j]] = 1;
	      break;	    	    
	    default:
	      assert(0);
	    }
	}
    }
	      
  for(i = 0, countStatesPresent = 0; i < (size_t)numFreqs; i++)
    {      
      if(statesPresent[i] == 1)
	countStatesPresent++;
    }  

  for (k = 1; k <= 8; k++) 
    {	     	   	    	      			    
      for(l = 0; l < numFreqs; l++)
	sumf[l] = 0.0;
	      
      for(i = 0; i < (size_t)rdta->numsp; i++) 
	{		 
	  yptr = &(rdta->y0[((size_t)i) * (tr->originalCrunchedLength)]);
	  
	  for(j = lower; j < upper; j++) 
	    {
	      unsigned int 
		code = bitMask[yptr[j]];
	      
	      assert(code >= 1);
	      
	      for(l = 0; l < numFreqs; l++)
		{
		  if((code >> l) & 1)
		    temp[l] = pfreqs[l];
		  else
		    temp[l] = 0.0;
		}		      	      
	      
	      for(l = 0, acc = 0.0; l < numFreqs; l++)
		{
		  if(temp[l] != 0.0)
		    acc += temp[l];
		}
	      
	      wj = ((double)cdta->aliaswgt[j]) / acc;
	      
	      for(l = 0; l < numFreqs; l++)
		{
		  if(temp[l] != 0.0)		    
		    sumf[l] += wj * temp[l];			     				   			     		   
		}
	    }
	}	    	      
      
      for(l = 0, acc = 0.0; l < numFreqs; l++)
	{
	  if(sumf[l] != 0.0)
	    acc += sumf[l];
	}
	      
      for(l = 0; l < numFreqs; l++)
	pfreqs[l] = sumf[l] / acc;	     
    }

  
  
  if(countStatesPresent < numFreqs)
    {           
      printf("Partition %s number %d has a problem, the number of expected states is %d the number of states that are present is %d.\n", 
	     tr->partitionData[model].partitionName, model, numFreqs, countStatesPresent);
      printf("Please go and fix your data!\n\n");   
    }

  if(smoothFrequencies)         
    {          
      smoothFreqs(numFreqs, pfreqs,  tr->partitionData[model].frequencies, &(tr->partitionData[model]));	   
    }
  else    
    {
      boolean 
	zeroFreq = FALSE;

      char 
	typeOfData[1024];

      getDataTypeString(tr, model, typeOfData);  

      for(l = 0; l < numFreqs; l++)
	{
	  if(pfreqs[l] == 0.0)
	    {	      
	      printBothOpen("Empirical base frequency for state number %d is equal to zero in %s data partition %s\n", l, typeOfData, tr->partitionData[model].partitionName);
	      printBothOpen("Since this is probably not what you want to do, RAxML will soon exit.\n\n");
	      zeroFreq = TRUE;	   
	    }      
	}
      
      if(zeroFreq)
	exit(-1);
	

      for(l = 0; l < numFreqs; l++)
	{ 	 
	  assert(pfreqs[l] > 0.0);
	  tr->partitionData[model].frequencies[l] = pfreqs[l];
	}          
    }  
}







static void baseFrequenciesGTR(rawdata *rdta, cruncheddata *cdta, tree *tr)
{  
  int
    model;

  size_t
    lower,
    upper;
  
  int
    states;

  for(model = 0; model < tr->NumberOfModels; model++)
    {      
      lower = tr->partitionData[model].lower;
      upper = tr->partitionData[model].upper;	  	 
      states = tr->partitionData[model].states;
	
      switch(tr->partitionData[model].dataType)
	{
	case GENERIC_32:
	  switch(tr->multiStateModel)
	    {
	    case ORDERED_MULTI_STATE:
	    case MK_MULTI_STATE:	   
	      {	       
		int 
		  i;
		double 
		  freq = 1.0 / (double)states,
		  acc = 0.0;

		for(i = 0; i < states; i++)
		  {
		    acc += freq;
		    tr->partitionData[model].frequencies[i] = freq;
		    /*printf("%f \n", freq);*/
		  }
		/*printf("Frequency Deviation: %1.60f\n", acc);*/
	      }
	      break;
	     case GTR_MULTI_STATE:
	      genericBaseFrequencies(tr, states, rdta, cdta, lower, upper, model, TRUE,
				     bitVector32);
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case GENERIC_64:	 
	  assert(0);
	  break;
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case SECONDARY_DATA:
	case AA_DATA:
	case DNA_DATA:
	case BINARY_DATA:	  
	  //mth calc empirical base freqs to string

	  genericBaseFrequencies(tr, states, rdta, cdta, lower, upper, model, 
				 getSmoothFreqs(tr->partitionData[model].dataType),
				 getBitVector(tr->partitionData[model].dataType));
	  break;	
	case POMO_16:  
	case POMO_64:	 
	  memset(tr->partitionData[model].frequencies, 0, sizeof(double) * (size_t)tr->partitionData[model].states);
	  genericBaseFrequencies(tr, 4, rdta, cdta, lower, upper, model, 
				 getSmoothFreqs(DNA_DATA),
				 getBitVector(DNA_DATA));	  	 		  
	  break;
	default:
	  assert(0);     
	}      
    }
  
  return;
}

int lookupWord(char *s, stringHashtable *h)
{
  hashNumberType position = hashString(s, h->tableSize);
  stringEntry *p = h->table[position];

  for(; p!= NULL; p = p->next)
    {
      if(strcmp(s, p->word) == 0)		 
	return p->nodeNumber;	  	
    }

  return -1;
}

//mth 

static void calculatePomoMap(tree *tr, analdef *adef)
{
  if(adef->model == M_POMOGAMMA_16 || adef->model == M_POMOGAMMA_64)
    {
      //mth open the mapping file 

      FILE 
	*mapFile = myfopen(pomoMapFileName, "r");      

      int 
	taxaCount = 0,
	i, 
	j,
	speciesCount = 0,
	individualsCount = 0;
      
      char 
	*cc = (char *)NULL;
      
      size_t 
	nbytes = 0;

      //mth count lines in mapping file to obtain the number of species 
     
      while(getline(&cc, &nbytes, mapFile) > -1)
	if(!lineContainsOnlyWhiteChars(cc))
	  speciesCount++;

      //mth set number of POMO species 

      tr->numberOfPomoSpecies = speciesCount;

      //mth less than 4 species -> return an error 

      if(tr->numberOfPomoSpecies < 4)
	{
	  printf("\nError, we need at least 4 species in the map file %s for building a POMO-based tree, exiting\n\n", pomoMapFileName);
	  errorExit(-1);
	}

      //mth set up some arrays for the individual to species mapping 
      
      tr->pomoIndex = (pomoInd *)malloc(sizeof(pomoInd) * ((size_t)speciesCount));
      tr->pomoMap = (int *)malloc(sizeof(int) * ((size_t)tr->mxtips + 1));
      tr->pomoSpeciesNameList= (char **)malloc(sizeof(char *) * ((size_t)speciesCount));

      printBothOpen("\nNumber of POMO species: %d\n\n", speciesCount);

      speciesCount = 0;

      //mth go back to the beginning of the file to start parsing taxon names 

      rewind(mapFile);

      while(getline(&cc, &nbytes, mapFile) > -1)
	{ 
	  if(!lineContainsOnlyWhiteChars(cc))
	    {
	      const char 
		delimiters[] = " \n\r";
	      
	      char 
		*token = strtok(cc, delimiters);   		 

	      boolean 
		speciesNameSet = FALSE;
	      
	      individualsCount = 0;

	      //mth loop over string tokens in a line 

	      while(token != (char*)NULL) 
		{
		  //mth first name in list is the species name 		  
		  
		  if(!speciesNameSet)
		    {		      
		      tr->pomoSpeciesNameList[speciesCount] = (char *)malloc((strlen(token) + 1) * sizeof(char));
		      strcpy(tr->pomoSpeciesNameList[speciesCount], token);

		      if(lookupWord(token, tr->speciesHash) != -1)
			{
			  printf("\nError: duplicate POMO Species name %s in file %s!\n\n", token, pomoMapFileName);
			  errorExit(-1);
			}
		      
		      addword(token, tr->speciesHash, (int)(speciesCount + 1));
		      
		      printBothOpen("POMO Species name: %s\n", tr->pomoSpeciesNameList[speciesCount]);

		      speciesNameSet = TRUE; 
		    }
		  else
		    {
		      //mth check if individual taxon name is contained in the MSA 
		      
		      int 
			lookup = lookupWord(token, tr->nameHash);
		      
		      //mth if not, exit 
		      if(lookup <= 0)
			{
			  printf("\n Taxon %s from pomo map file %s can not be found in alignment, exiting ... \n\n", token, pomoMapFileName);
			  errorExit(-1);
			}
		      else
			individualsCount++;
		      
		      //mth set up mapping array that tells us that individual with the taxon identifier lookup belongs to 
		      //mth species speciesCount 

		      tr->pomoMap[lookup] = speciesCount;
		    }

		  token = strtok((char*)NULL, delimiters);		    
		}
	      
	      //mth set up a struct data type that contains the individuals to species mapping

	      //mth number of individuals per species 
	      tr->pomoIndex[speciesCount].indCount = individualsCount;

	      //mth allocate array for storing individual taxon indices 
	      tr->pomoIndex[speciesCount].indMap = (int *)malloc(sizeof(int) * (size_t)individualsCount);

	      taxaCount += individualsCount;

	      //mth we need at least two individuals per species 
	      if(individualsCount < 2)
		{
		  printf("\nFor a POMO model the number of individuals per species needs to be at least 2!\n");
		  printf("The species map file line: \n\n%s \n\ncontains only one, exiting ....\n\n", cc);
		  errorExit(-1);
		}

	      speciesCount++;
	    }
	}               
      
      printBothOpen("\n\n");
      
      fclose(mapFile);

      assert(taxaCount == tr->mxtips);      

      //mth figure in the individual taxon indices to the struct that defines the POMO species 

      for(i = 0; i < tr->numberOfPomoSpecies; i++)
	{
	  int
	    counter;
	  
	  for(j = 1, counter = 0; j <= tr->mxtips; j++)	
	    if(i == tr->pomoMap[j])	    
	      tr->pomoIndex[i].indMap[counter++] = j;	           
	}             
    }
}


//mth functions for building POMO CLV

/* Redargless of the number of bins of frequency differences, there are only
    10 classes of Pomo states: the 4 monomorphic, and the 6 diallelic.
*/
typedef enum PomoStateClass {
  MONO_A_STATE_CLASS = 0,
  MONO_C_STATE_CLASS = 1,
  MONO_G_STATE_CLASS = 2,
  MONO_T_STATE_CLASS = 3,
  LAST_MONO_STATE_CLASS = 3,
  DIALLELIC_AC_STATE_CLASS = 4,
  DIALLELIC_AG_STATE_CLASS = 5,
  DIALLELIC_AT_STATE_CLASS = 6,
  DIALLELIC_CG_STATE_CLASS = 7,
  DIALLELIC_CT_STATE_CLASS = 8,
  DIALLELIC_GT_STATE_CLASS = 9,
  LAST_POMO_STATE_CLASS = 9
} PomoStateClass_t;

//MTH obsToPomoCounts maps a binary representation of a data code to:
//  an array of it's effects on the set of model state classes
                                     //    A,  C,  G,  T, AC, AG, AT, CG, CT, GT
static const int obsToPomoCounts[16][10] = {{15, 15, 15, 15, 15, 15, 15, 15, 15, 15},
					    { 4,  1,  1,  1,  4,  4,  4,  1,  1,  1}, // A
					    { 1,  4,  1,  1,  2,  1,  1,  4,  4,  1}, // C
					    { 4,  4,  1,  1,  0,  4,  4,  4,  4,  1}, // M = {AC}
					    { 1,  1,  4,  1,  1,  2,  1,  2,  1,  4}, // G
					    { 4,  1,  4,  1,  4,  0,  4,  2,  1,  4}, // R = {AG}
					    { 1,  4,  4,  1,  2,  2,  1,  0,  4,  4}, // S = {CG}
					    { 4,  4,  4,  1,  0,  0,  4,  0,  4,  4}, // V = {ACG}
					    { 1,  1,  1,  4,  1,  1,  2,  1,  2,  2}, // T
					    { 4,  1,  1,  4,  4,  4,  0,  1,  2,  2}, // W = {AT}
					    { 1,  4,  1,  4,  2,  1,  2,  4,  0,  2}, // Y = {CT}
					    { 4,  4,  1,  4,  0,  4,  0,  4,  0,  2}, // H = {ACT}
					    { 1,  1,  4,  4,  1,  2,  2,  2,  2,  0}, // K = {GT}
					    { 4,  1,  4,  4,  4,  0,  0,  2,  2,  0}, // D = {AGT}
					    { 1,  4,  4,  4,  2,  2,  2,  0,  0,  0}, // B = {CGT}
					    { 4,  4,  4,  4,  0,  0,  0,  0,  0,  0}}; // N = {ACGT}

static double calcBinomProb(unsigned numFirst, unsigned numSecond, double probFirst);

static double logBinomCoefficient(unsigned numFirst, unsigned numSecond);

static double logBinomCoefficient(unsigned numFirst, unsigned numSecond)
{
  const unsigned 
    n = numFirst + numSecond,
    larger = (numFirst > numSecond ? numFirst : numSecond);
  
  if(larger == n)    
    return 0.0;
  else
    {
      unsigned 
	i;
      
      double 
	logp = 0.0;
      
      for(i = 0; i + larger < n; ++i)
	{
	  double 
	    numeratorFactor = (double)(n - i),
	    denominatorFactor = (double)(1 + i);

	  logp += log(numeratorFactor) - log(denominatorFactor);
	}
      
      return logp;
    }
}

static double calcBinomProb(unsigned int numFirst, unsigned int numSecond, double probFirst)
{
  if(numSecond == 0 && numFirst == 0)    
    return 1.0;
  else
    {
      double 
	logp = 0.0;
      
      if(numFirst > 0)
	{
	  assert(probFirst > 0.0);
	  logp += numFirst*log(probFirst);
	} 
      
      if(numSecond > 0)
	{
	  const double 
	    probSecond = 1.0 - probFirst;

	  assert(probSecond > 0.0);
	  logp += numSecond*log(probSecond);
	}
      
      logp += logBinomCoefficient(numFirst, numSecond);
      
      return exp(logp);
    }
}
//mth function for building POMO CLV, to be extended 

static void buildPomoCLV(size_t i, double *pomoBuffer, tree *tr, pInfo *p, unsigned char *y0)
{
  size_t 
    j,
    site;
  
  const size_t 
    numDialleleFreqBins = ((size_t)p->states - 4) / 6 ; // p->states-4 is the number of states for diallelic conditions. There are 6 diallelic combinations.
  
  const double 
    binWidth = 1.0 / ((double)(1 + numDialleleFreqBins));
  
  //mth we allow only for max two states per inds in a species; how do we handle ambiguous characters?
  //mth right now they if we have a site with A and a site with C and a third site with an ambiguous character representing A or G
  //mth this is still accepted by the code

  //mth loop over partition sites 

  for(site = p->lower; site < p->upper; site++)    
    {   
      size_t
	sc;
      
      double 
	cl;
      //mth loop over DNA sequences of individuals for POMO species i
      
      int 
	stillValid[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
	numStillValid = 10;
	// Same dim of full PomoStateClass_t, but only the diallelic cases will be filled
      unsigned int
	diallelicCounts[10][2] = {{0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}};
      
      for(j = 0; j < (size_t)tr->pomoIndex[i].indCount; j++)	
	{
	  const int 
	    taxonIndex = tr->pomoIndex[i].indMap[j];
	  
	  const unsigned char 
	    tipValue = (y0 + sizeof(unsigned char)
			* ((((size_t)taxonIndex - 1) *  tr->originalCrunchedLength)  + site))[0];

	  //printf("i=%lu j=%lu site=%lu tipValue=%u\n", i, j, site, tipValue); //DEBUG

	  if((tipValue < 1) || (tipValue > 15))
	    {
	      printf("\n Invalid code for a DNA state!\n");
	      errorExit(-1);
	    }
	  
	  const int 
	    *effectRow = obsToPomoCounts[tipValue];
	  
	  for (sc = 0; sc <= LAST_MONO_STATE_CLASS; ++sc)
	    {
	      if (effectRow[sc] & 1)
		{
		  if (stillValid[sc])
		    {
		      if(numStillValid == 1)
			{
			  printf("\n Column %d of species %d cannot be explained by PoMo - more than 2 alleles/species required!\n", (int)j, (int)i);
			  errorExit(-1);
			}
		      numStillValid -= 1;
		      stillValid[sc] = 0;
		    }
		}
	    }
	  
	  for(; sc <= LAST_POMO_STATE_CLASS; ++sc)
	    {
	      if(effectRow[sc] & 1)
		{
		  if (stillValid[sc])
		    {
		      if (numStillValid == 1)
			{
			  printf("\n Column %d of species %d cannot be explained by PoMo - more than 2 alleles/species required!\n", (int)j, (int)i);
			  errorExit(-1);
			}
		      numStillValid -= 1;
		      stillValid[sc] = 0;
		    }
		}
	      else
		{
		  if (effectRow[sc] & 4) // signal that first allele is the only one compat w/ this state		    
		    diallelicCounts[sc][0] += 1;		   
		  else 
		    {
		      if (effectRow[sc] & 2) // signal that second allele is the only one compat w/ this state		    
			diallelicCounts[sc][1] += 1;
		    }
		}
	    }
	}
      
      size_t 
	pbIndex;
      //mth set pomo CLV buffer entries for a specific site (indexing of the buffer starts at 0) to some value 
      
      for (sc = 0; sc <= LAST_MONO_STATE_CLASS; ++sc)
	{
	  pbIndex = (site - p->lower) * (size_t)p->states + sc;
	  cl = (stillValid[sc] ? 1.0 : 0.0);
	  pomoBuffer[pbIndex] = cl;
	  ///printf("    pomoBuffer[%lu] = %lf\n", pbIndex, cl); //DEBUGGING
	}
      
      size_t 
	diallelicOffset = 1 + LAST_MONO_STATE_CLASS;

      for (; sc <= LAST_POMO_STATE_CLASS; ++sc, diallelicOffset += numDialleleFreqBins)
	{
	  if (stillValid[sc])
	    {
	      for (j = 0 ; j < numDialleleFreqBins; ++j)
		{
		  double
		    secondAlleleFreq = binWidth * ((double)(1 + j)),
		    firstAlleleFreq = 1.0 - secondAlleleFreq;
		  
		  pbIndex = (site - p->lower) * (size_t)p->states + diallelicOffset + j;
		  cl = calcBinomProb(diallelicCounts[sc][0], diallelicCounts[sc][1], firstAlleleFreq);
		  pomoBuffer[pbIndex] = cl;
		  //printf("    pomoBuffer[%lu] = %lf\n", pbIndex, cl); //DEBUGGING
		}
	    }
	  else
	    {
	      pbIndex = (site - p->lower) * (size_t)p->states + diallelicOffset;
	      
	      for (j = 0 ; j < numDialleleFreqBins; ++j)		
		pomoBuffer[pbIndex + j] = 0.0;	       
	      //printf("    pomoBuffer[%lu] -> pomoBuffer[%lu] = 0.0\n", pbIndex, pbIndex + numDialleleFreqBins - 1); //DEBUGGING
	    }
	}
    }     
}



// #define OLD_LAYOUT 


int main (int argc, char *argv[])
{
  size_t 
    model;

  rawdata      *rdta;
  cruncheddata *cdta;
  tree         *tr;
  analdef      *adef;
  
  /* get the start time */

  masterTime = gettime();

  /* get some memory for the basic data structures */

  adef = (analdef *)malloc(sizeof(analdef));
  rdta = (rawdata *)malloc(sizeof(rawdata));
  cdta = (cruncheddata *)malloc(sizeof(cruncheddata));
  tr   = (tree *)malloc(sizeof(tree));

  /* initialize the analysis parameters in struct adef to default values */

  initAdef(adef);

  /* parse command line arguments: this has a side effect on tr struct and adef struct variables */

  get_args(argc,argv, adef, tr); 
            
  /* parse the phylip file: this should probably be re-done, perhaps using the relatively flexible parser 
     written in C++ by Marc Holder */
  
  getinput(adef, rdta, cdta, tr);  

  printBothOpen("Pattern compression: %s\n", (adef->compressPatterns)?"ON":"OFF");

  makeweights(adef, rdta, cdta, tr);         
      
  makevalues(rdta, cdta, tr, adef);                 
        
  //mth parses the POMO individuals to species mapping file and stores the mapping 

  calculatePomoMap(tr, adef);
                  
  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {	      
      tr->partitionData[model].states = getStates(tr->partitionData[model].dataType);
      tr->partitionData[model].maxTipStates = getUndetermined(tr->partitionData[model].dataType) + 1;  	      
      tr->partitionData[model].nonGTR = FALSE;
      
      partitionLengths 
	*pl = getPartitionLengths(&(tr->partitionData[model]));
      
      tr->partitionData[model].frequencies       = (double*)malloc((size_t)pl->frequenciesLength * sizeof(double));      
    }   
  
  baseFrequenciesGTR(tr->rdta, tr->cdta, tr); 

  {
    int 
      sizeOfSizeT = sizeof(size_t),
      version = (int)programVersionInt,
      magicNumber = 6517718;
    
    size_t 
      i;   
    
    /* NEW, we firstly write, how many bytes size_t comprises */
    
    myBinFwrite(&(sizeOfSizeT),                sizeof(sizeOfSizeT), 1); 

    //error checking for parser!
    myBinFwrite(&version,     sizeof(int), 1);
    myBinFwrite(&magicNumber, sizeof(int), 1);
    //error checking for correct parser end

    
    //mth write the number of POMO species to file instead of the number of taxa/individuals 
    if(adef->model == M_POMOGAMMA_16 || adef->model == M_POMOGAMMA_64)
      myBinFwrite(&(tr->numberOfPomoSpecies),                 sizeof(int), 1);
    else
      myBinFwrite(&(tr->mxtips),                 sizeof(int), 1);

    myBinFwrite(&(tr->originalCrunchedLength), sizeof(size_t), 1);
    myBinFwrite(&(tr->NumberOfModels),         sizeof(int), 1);
    myBinFwrite(&(adef->gapyness),             sizeof(double), 1);
    
    myBinFwrite(tr->cdta->aliaswgt,               sizeof(int), tr->originalCrunchedLength);	  	  	       	
	
    if(adef->model == M_POMOGAMMA_16 || adef->model == M_POMOGAMMA_64)
      {
	for(i = 0; i < (size_t)tr->numberOfPomoSpecies; i++)
	  {
	    int
	      len = (int)strlen(tr->pomoSpeciesNameList[i]) + 1;
	    
	    myBinFwrite(&len, sizeof(int), 1);
	    myBinFwrite(tr->pomoSpeciesNameList[i], sizeof(char), (size_t)len);	
	  }  
      }
    else
      {
	for(i = 1; i <= (size_t)tr->mxtips; i++)
	  {
	    size_t 
	      len = strlen(tr->nameList[i]) + 1;
	    
	    myBinFwrite(&len, sizeof(int), 1);
	    myBinFwrite(tr->nameList[i], sizeof(char), (size_t)len);	
	  }  
      }
    
    // mth: here we loop over the partitions and print information on partitions (not the actual seq data) to the binary alignment file
    // mth that is read by ExaML 

    for(model = 0; model < (size_t)tr->NumberOfModels; model++)
      {
	size_t 
	  len;
	
	pInfo 
	  *p = &(tr->partitionData[model]);
	
	//mth number of states
	myBinFwrite(&(p->states),             sizeof(int), 1);

	//mth number of max tips states including ambiguous characters 
	myBinFwrite(&(p->maxTipStates),       sizeof(int), 1);

	//mth lower and upper index of this partition
	myBinFwrite(&(p->lower),              sizeof(size_t), 1);
	myBinFwrite(&(p->upper),              sizeof(size_t), 1);

	//mth width of the partition in terms of # patterns
	myBinFwrite(&(p->width),              sizeof(size_t), 1);

	//mth datatype of this partition: BIN, DNA, PROT 
	myBinFwrite(&(p->dataType),           sizeof(int), 1);

	//mth prot-specific stuff 
	myBinFwrite(&(p->protModels),         sizeof(int), 1);	
	myBinFwrite(&(p->protFreqs),          sizeof(int), 1);	

	//mth are we using a standard GTR model for this partition?
	myBinFwrite(&(p->nonGTR),                      sizeof(boolean), 1); 	

	//mth do we want to do a ML estimate of the base freqs? 
	myBinFwrite(&(p->optimizeBaseFrequencies),     sizeof(boolean), 1);
				
	//mth name of the partition 
	len = strlen(p->partitionName) + 1;
	myBinFwrite(&len, sizeof(int), 1);
	myBinFwrite(p->partitionName, sizeof(char), (size_t)len);	    

	//mth empirical base frequencies for this partition 
	//mth right now for POMO the base freqs are partially undefined, we need to figure out how to best solve this 
	//mth the first 4 entries of the freq vector contain the empirical DNA base frequencies 
	myBinFwrite(tr->partitionData[model].frequencies, sizeof(double), (size_t)tr->partitionData[model].states);		



      }	            

#ifdef OLD_LAYOUT
    myBinFwrite(rdta->y0, sizeof(unsigned char), (tr->originalCrunchedLength) * ((size_t)tr->mxtips)); 
#else 
    /* 
       Write each partition, taxon by taxon. Thus, if unpartitioned,
       nothing changes.
    */   

    //mth write sequence data to the binary file ....

    size_t
      pomo_multiplier = 0,
      mem_reqs_cat = 0,
      mem_reqs_gamma = 0,
      unique_patterns = 0;

    for(model = 0; model < (size_t) tr->NumberOfModels; ++model )
      {
        pInfo
          *p  = &(tr->partitionData[model]); 
	
        size_t 
          width = p->upper - p->lower; 

	unique_patterns += width;

	//multiply partition width with number of states we need to store in each CLV entry

	mem_reqs_cat += (size_t)tr->partitionData[model].states * width;	

	//mth if we use a POMO model write a CLV to file intsead of the raw alignment sequence 

	if(p->dataType == POMO_16 || p->dataType == POMO_64)
	  {
	    //mth pomoBuffer assumes no rate heterogeneity, in fact even under GAMMA this is not needed 
	    //mth it suffices to store one value per site and per state at the tips 
	    //mth buffer for storing CLV we want write to the binary file 
	    double 
	      *pomoBuffer = (double*)malloc(sizeof(double) * (size_t)tr->partitionData[model].states * width);

	    //mth loop overd POMO species 
	    for(i = 0; i < (size_t)tr->numberOfPomoSpecies; i++)
	      {
		int 
		  j;			

		//mth some verbatim output 
		printBothOpen("\nBuilding CLV for POMO species %zu comprising the following individuals:\n", i);

		for(j = 0; j < tr->pomoIndex[i].indCount; j++)	
		  printBothOpen("%s ", tr->nameList[tr->pomoIndex[i].indMap[j]]);
		printBothOpen("\n");

		//mth the function below builds a CLV for the specific partition and and POMO species 
		//mth that shall be stored in pomoBuffer

		buildPomoCLV(i, pomoBuffer, tr, p, rdta->y0);
		

		//mth write POMO tip CLV to binary alignment file 		
		myBinFwrite(pomoBuffer,  sizeof(double), (size_t)tr->partitionData[model].states * width);		
	      }

	    free(pomoBuffer);
	  }
	else
	  {
	    for(i = 0; i < (size_t)tr->mxtips; ++i)	      
	      myBinFwrite(rdta->y0
			  + sizeof(unsigned char) * (  (i *  tr->originalCrunchedLength)  + p->lower   ) 
			  , sizeof(unsigned char), width); 	     
	  }
      }

    printBothOpen("\n\nYour alignment has %zu %s\n", unique_patterns, (adef->compressPatterns == TRUE)?"unique patterns":"sites");

    //multiply CLV vector length with number of tips and 8, since b bytes are needed to store an inner conditional probability vector    
    mem_reqs_cat *= (size_t)tr->mxtips * sizeof(double);

    //mem reqs for gamma are 4 times higher than for CAT
    mem_reqs_gamma = mem_reqs_cat * 4;

    //now add the space for storing the tips:
    //take care of pomo CLV tips first 
    
    switch(adef->model)
      {
      case M_POMOGAMMA_16:
	pomo_multiplier = sizeof(double) * 16 * 2; //need to store two CLC copies at tips 
	break;
      case M_POMOGAMMA_64:
	pomo_multiplier = sizeof(double) * 64 * 2; //need to store two CLC copies at tips 
	break;
      default:
	pomo_multiplier = sizeof(unsigned char);
      }    

    mem_reqs_cat   += (size_t)tr->mxtips * unique_patterns * pomo_multiplier;
    mem_reqs_gamma += (size_t)tr->mxtips * unique_patterns * pomo_multiplier;
     
    printBothOpen("\n\nUnder CAT the memory required by ExaML for storing CLVs and tip vectors will be\n%zu bytes\n%zu kiloBytes\n%zu MegaBytes\n%zu GigaBytes\n", 
		  mem_reqs_cat, 
		  mem_reqs_cat / 1024 , 
		  mem_reqs_cat / (1024 * 1024),
		  mem_reqs_cat / (1024 * 1024 * 1024));
    
    printBothOpen("\n\nUnder GAMMA the memory required by ExaML for storing CLVs and tip vectors will be\n%zu bytes\n%zu kiloBytes\n%zu MegaBytes\n%zu GigaBytes\n", 
		  mem_reqs_gamma, 
		  mem_reqs_gamma / 1024 , 
		  mem_reqs_gamma / (1024 * 1024),
		  mem_reqs_gamma / (1024 * 1024 * 1024));

    printBothOpen("\nPlease note that, these are just the memory requirements for doing likelihood calculations!\n");
    printBothOpen("To be on the safe side, we recommend that you execute ExaML on a system with twice that memory.\n");

#endif
  }

  fclose(byteFile);  

  printBothOpen("\n\nBinary and compressed alignment file written to file %s\n\n", byteFileName);
  printBothOpen("Parsing completed, exiting now ... \n\n");

  free(adef);

  free(tr->model);
  free(tr->initialDataVector);
  free(tr->extendedDataVector);
  free(tr->initialPartitionData);
  
  free(rdta->wgt);
  free(cdta->alias);          
  free(cdta->aliaswgt);

  free(rdta);
  free(cdta); 
  free(tr);
  
  return 0;
}
