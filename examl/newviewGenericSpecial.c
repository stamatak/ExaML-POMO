/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *
 *  and
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
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
#include <unistd.h>
#endif

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include "axml.h"

#ifdef __SIM_SSE3

#include <stdint.h>
#include <xmmintrin.h>
#include <pmmintrin.h>

/* required to compute the absoliute values of double precision numbers with SSE3 */

typedef const union __attribute__ ((aligned (BYTE_ALIGNMENT)))
{
       uint64_t i[2];
       __m128d m;
} aMask;

extern const aMask absMask;

const aMask absMask = {{0x7fffffffffffffffULL , 0x7fffffffffffffffULL }};

/*const union __attribute__ ((aligned (BYTE_ALIGNMENT)))
{
       uint64_t i[2];
       __m128d m;
       } absMask = {{0x7fffffffffffffffULL , 0x7fffffffffffffffULL }};*/



#endif

/*** generic n-state model ***/
#if defined(__SIM_SSE3) && !defined(__AVX)
#include <xmmintrin.h>
#include <pmmintrin.h>


#define VECTOR_STORE_LEFT _mm_storel_pd


static const union __attribute__ ((aligned (BYTE_ALIGNMENT)))
{
       uint64_t i[2];
       __m128d m;
} absMaskGeneric = {{0x7fffffffffffffffULL , 0x7fffffffffffffffULL }};

#endif



#ifdef __AVX
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <immintrin.h>

static const union __attribute__ ((aligned (BYTE_ALIGNMENT)))
{
  int32_t i[8];
  __m256i m;  
} bitmask = {{0, 0, 0, 0, 0, 0, -1, -1}};


#define VECTOR_STORE_LEFT(x,y) _mm256_maskstore_pd(x, bitmask.m, y)


static const union __attribute__ ((aligned (BYTE_ALIGNMENT)))
{
  uint64_t i[4];
  __m256d m;
  
} absMaskGeneric = {{0x7fffffffffffffffULL, 0x7fffffffffffffffULL, 0x7fffffffffffffffULL, 0x7fffffffffffffffULL}};

#endif

static double haddScalar(VECTOR_REGISTER v)
{
  double 
    result;

#if defined(__SIM_SSE3) && !defined(__AVX)
  
  v = _mm_hadd_pd(v, v);
  
  _mm_storel_pd(&result, v);
#endif
  
#ifdef __AVX 
  double 
    ra[4] __attribute__ ((aligned (BYTE_ALIGNMENT)));

  v = _mm256_hadd_pd(v, v);

  _mm256_store_pd(ra, v);

  result = ra[0] + ra[2];
#endif

  return result;
}

static VECTOR_REGISTER haddBroadCast(VECTOR_REGISTER v)
{
#if defined(__SIM_SSE3) && !defined(__AVX)
  
  return _mm_hadd_pd(v, v);  
   
#endif

#ifdef __AVX
  __m256d
    a;
  
  v = _mm256_hadd_pd(v, v);
  a = _mm256_permute2f128_pd(v, v, 1);
  v = _mm256_add_pd(a, v);
  
  return v;
#endif
}

static boolean scaleEntry(const size_t stride, const size_t i, double *x3, const size_t scalingLoopLength)
{
  double 
    *v = &(x3[stride * i]);
  
  VECTOR_REGISTER 
    minlikelihood_vector = VECTOR_SET_ONE( minlikelihood );

  size_t
    l;

  int
    scale = 1;
  
  for(l = 0; scale && (l < scalingLoopLength); l += VECTOR_WIDTH)
    {
      VECTOR_REGISTER 
	vv = VECTOR_LOAD(&v[l]);
      
      VECTOR_REGISTER 
	v1 = VECTOR_AND(vv, absMaskGeneric.m);

#if defined(__SIM_SSE3) && !defined(__AVX)
      v1 = _mm_cmplt_pd(v1,  minlikelihood_vector);
      if(_mm_movemask_pd( v1 ) != 3)
	scale = 0;
#endif
      
#ifdef __AVX
      v1 = _mm256_cmp_pd(v1,  minlikelihood_vector, _CMP_LT_OS);
      if(_mm256_movemask_pd( v1 ) != 15)
	scale = 0;
#endif
      
    }	    	  
	      
  for(;scale && (l < stride); l++)
    scale = (ABS(v[l]) < minlikelihood);

  if(scale)
    {
      VECTOR_REGISTER 
	twoto = VECTOR_SET_ONE(twotothe256);	       

      for(l = 0; l < scalingLoopLength; l += VECTOR_WIDTH)
	{
	  VECTOR_REGISTER 
	    ex3v = VECTOR_LOAD(&v[l]);		  
	  
	  VECTOR_STORE(&v[l], VECTOR_MUL(ex3v,twoto));	
	}		   		  

      for(;l < stride; l++)
	v[l] *= twotothe256;
      
      //      printf("scale\n");
      return TRUE;
    }
  else
    {
      //printf("no scale\n");
      return FALSE;
    }
}



/* includes MIC-optimized functions */

#ifdef __MIC_NATIVE
#include "mic_native.h"
#endif

extern int processID;

/* bit mask */

extern const unsigned int mask32[32];


/* generic function for computing the P matrices, for computing the conditional likelihood at a node p, given child nodes q and r 
   we compute P(z1) and P(z2) here */

static void makeP(double z1, double z2, double *rptr, double *EI,  double *EIGN, size_t numberOfCategories, double *left, double *right, boolean saveMem, int maxCat, const size_t states)
{
  size_t
    i, 
    j, 
    k,
    /* square of the number of states = P-matrix size */
    statesSquare = states * states;
  
  /* assign some space for pre-computing and later re-using functions */

  double 
    *lz1 = (double*)malloc(sizeof(double) * states),
    *lz2 = (double*)malloc(sizeof(double) * states),
    *d1 = (double*)malloc(sizeof(double) * states),
    *d2 = (double*)malloc(sizeof(double) * states);

  /* multiply branch lengths with eigenvalues */

  for(i = 1; i < states; i++)
    {
      lz1[i] = EIGN[i] * z1;
      lz2[i] = EIGN[i] * z2;
    }


  /* loop over the number of rate categories, this will be 4 for the GAMMA model and 
     variable for the CAT model */

  for(i = 0; i < numberOfCategories; i++)
    {
      /* exponentiate the rate multiplied by the branch */

      for(j = 1; j < states; j++)
	{
	  d1[j] = EXP(rptr[i] * lz1[j]);
	  d2[j] = EXP(rptr[i] * lz2[j]);
	}

      /* now fill the P matrices for the two branch length values */

      for(j = 0; j < states; j++)
	{
	  /* left and right are pre-allocated arrays */

	  left[statesSquare * i  + states * j] = 1.0;
	  right[statesSquare * i + states * j] = 1.0;	  

	  for(k = 1; k < states; k++)
	    {
	      left[statesSquare * i + states * j + k]  = d1[k] * EI[states * j + k];
	      right[statesSquare * i + states * j + k] = d2[k] * EI[states * j + k];
	    }
	}
    }


  /* if memory saving is enabled and we are using CAT we need to do one additional P matrix 
     calculation for a rate of 1.0 to compute the entries of a column/tree site comprising only gaps */


  if(saveMem)
    {
      i = (size_t)maxCat;
      
      for(j = 1; j < states; j++)
	{
	  d1[j] = EXP (lz1[j]);
	  d2[j] = EXP (lz2[j]);
	}

      for(j = 0; j < states; j++)
	{
	  left[statesSquare * i  + states * j] = 1.0;
	  right[statesSquare * i + states * j] = 1.0;

	  for(k = 1; k < states; k++)
	    {
	      left[statesSquare * i + states * j + k]  = d1[k] * EI[states * j + k];
	      right[statesSquare * i + states * j + k] = d2[k] * EI[states * j + k];
	    }
	}
    }
  
  /* free the temporary buffers */

  free(lz1);
  free(lz2);
  free(d1);
  free(d2);
}

static void makeP_FlexLG4(double z1, double z2, double *rptr, double *EI[4],  double *EIGN[4], size_t numberOfCategories, double *left, double *right, const size_t numStates)
{
  size_t
    i,
    j,
    k;
  
  const size_t
    statesSquare = numStates * numStates;

  double    
    d1[64],  
    d2[64];

  assert(numStates <= 64);
       
  for(i = 0; i < numberOfCategories; i++)
    {
      for(j = 1; j < numStates; j++)
	{
	  d1[j] = EXP (rptr[i] * EIGN[i][j] * z1);
	  d2[j] = EXP (rptr[i] * EIGN[i][j] * z2);
	}

      for(j = 0; j < numStates; j++)
	{
	  left[statesSquare * i  + numStates * j] = 1.0;
	  right[statesSquare * i + numStates * j] = 1.0;

	  for(k = 1; k < numStates; k++)
	    {
	      left[statesSquare * i + numStates * j + k]  = d1[k] * EI[i][numStates * j + k];
	      right[statesSquare * i + numStates * j + k] = d2[k] * EI[i][numStates * j + k];
	    }
	}
    }  
}



/* The functions here are organized in a similar way as in evaluateGenericSpecial.c 
   I provide generic, slow but readable function implementations for computing the 
   conditional likelihood arrays at p, given child nodes q and r. Once again we need 
   two generic function implementations, one for CAT and one for GAMMA */

#ifndef _OPTIMIZED_FUNCTIONS

static void newviewCAT_FLEX(int tipCase, double *extEV,
			    int *cptr,
			    double *x1, double *x2, double *x3, double *tipVector,
			    unsigned char *tipX1, unsigned char *tipX2,
			    size_t n, double *left, double *right, int *wgt, int *scalerIncrement, const int states)
{
  double
    *le, 
    *ri, 
    *v, 
    *vl, 
    *vr,
    ump_x1, 
    ump_x2, 
    x1px2;

  size_t
    i, 
    l, 
    j;

  int
    scale, 
    addScale = 0;

  const size_t
    statesSquare = states * states;


  /* here we switch over the different cases for efficiency, but also because 
     each case accesses different data types.

     We consider three cases: either q and r are both tips, q or r are tips, and q and r are inner 
     nodes.
  */
     

  switch(tipCase)
    {
      
      /* both child nodes of p weher we want to update the conditional likelihood are tips */
    case TIP_TIP:     
      /* loop over sites */
      for (i = 0; i < n; i++)
	{
	  /* set a pointer to the P-Matrices for the rate category of this site */
	  le = &left[cptr[i] * statesSquare];
	  ri = &right[cptr[i] * statesSquare];
	  
	  /* pointers to the likelihood entries of the tips q (vl) and r (vr) 
	     We will do reading accesses to these values only.
	   */
	  vl = &(tipVector[states * tipX1[i]]);
	  vr = &(tipVector[states * tipX2[i]]);
	  
	  /* address of the conditional likelihood array entres at site i. This is 
	     a writing access to v */
	  v  = &x3[states * i];
	  
	  /* initialize v */
	  for(l = 0; l < states; l++)
	    v[l] = 0.0;
	  	  
	  /* loop over states to compute the cond likelihoods at p (v) */

	  for(l = 0; l < states; l++)
	    {	      
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;
	      
	      /* le and ri are the P-matrices */

	      for(j = 0; j < states; j++)
		{
		  ump_x1 += vl[j] * le[l * states + j];
		  ump_x2 += vr[j] * ri[l * states + j];
		}
	      
	      x1px2 = ump_x1 * ump_x2;
	      
	      /* multiply with matrix of eigenvectors extEV */

	      for(j = 0; j < states; j++)
		v[j] += x1px2 * extEV[l * states + j];
	    }	   
	}    
      break;
    case TIP_INNER:      

      /* same as above, only that now vl is a tip and vr is the conditional probability vector 
	 at an inner node. Note that, if we have the case that either q or r is a tip, the 
	 nodes will be flipped to ensure that tipX1 always points to the sequence at the tip.
      */

      for (i = 0; i < n; i++)
	{
	  le = &left[cptr[i] * statesSquare];
	  ri = &right[cptr[i] * statesSquare];
	  
	  /* access tip vector lookup table */
	  vl = &(tipVector[states * tipX1[i]]);

	  /* access conditional likelihoo arrays */
	  /* again, vl and vr are reading accesses, while v is a writing access */
	  vr = &x2[states * i];
	  v  = &x3[states * i];
	  
	  /* same as in the loop above */

	  for(l = 0; l < states; l++)
	    v[l] = 0.0;
	  
	  for(l = 0; l < states; l++)
	    {
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;
	      
	      for(j = 0; j < states; j++)
		{
		  ump_x1 += vl[j] * le[l * states + j];
		  ump_x2 += vr[j] * ri[l * states + j];
		}
	      
	      x1px2 = ump_x1 * ump_x2;
	      
	      for(j = 0; j < states; j++)
		v[j] += x1px2 * extEV[l * states + j];
	    }
	  
	  /* now let's check for numerical scaling. 
	     The maths in RAxML are a bit non-standard to avoid/economize on arithmetic operations 
	     at the virtual root and for branch length optimization and hence values stored 
	     in the conditional likelihood vectors can become negative.
	     Below we check if all absolute values stored at position i of v are smaller 
	     than a pre-defined value in axml.h. If they are all smaller we can then safely 
	     multiply them by a large, constant number twotothe256 (without numerical overflow) 
	     that is also speced in axml.h */

	  scale = 1;
	  for(l = 0; scale && (l < states); l++)
	    scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));	   
	  
	  if(scale)
	    {
	      for(l = 0; l < states; l++)
		v[l] *= twotothe256;
	      
	      /* if we have scaled the entries to prevent underflow, we need to keep track of how many scaling 
		 multiplications we did per node such as to undo them at the virtual root, e.g., in 
		 evaluateGeneric() 
		 Note here, that, if we scaled the site we need to increment the scaling counter by the wieght, i.e., 
		 the number of sites this potentially compressed pattern represents ! */ 

	      addScale += wgt[i];	  
	    }
	}   
      break;
    case INNER_INNER:
      
      /* same as above, only that the two child nodes q and r are now inner nodes */

      for(i = 0; i < n; i++)
	{
	  le = &left[cptr[i] * statesSquare];
	  ri = &right[cptr[i] * statesSquare];

	  /* index conditional likelihood vectors of inner nodes */

	  vl = &x1[states * i];
	  vr = &x2[states * i];
	  v = &x3[states * i];

	  for(l = 0; l < states; l++)
	    v[l] = 0.0;
	 
	  for(l = 0; l < states; l++)
	    {
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;

	      for(j = 0; j < states; j++)
		{
		  ump_x1 += vl[j] * le[l * states + j];
		  ump_x2 += vr[j] * ri[l * states + j];
		}

	      x1px2 =  ump_x1 * ump_x2;

	      for(j = 0; j < states; j++)
		v[j] += x1px2 * extEV[l * states + j];	      
	    }

	   scale = 1;
	   for(l = 0; scale && (l < states); l++)
	     scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));
  
	   if(scale)
	     {
	       for(l = 0; l < states; l++)
		 v[l] *= twotothe256;

	       addScale += wgt[i];	   
	     }
	}
      break;
    default:
      assert(0);
    }
   
  /* increment the scaling counter by the additional scalings done at node p */

  *scalerIncrement = addScale;
}


static void newviewGAMMA_FLEX(int tipCase,
			      double *x1, double *x2, double *x3, double *extEV, double *tipVector,
			      unsigned char *tipX1, unsigned char *tipX2,
			      size_t n, double *left, double *right, int *wgt, int *scalerIncrement, const int states, const int maxStateValue)
{
  double  
    *uX1, 
    *uX2, 
    *v, 
    x1px2, 
    *vl, 
    *vr, 
    al, 
    ar;
  
  int  
    i, 
    j, 
    l, 
    k, 
    scale, 
    addScale = 0;

  const int     
    statesSquare = states * states,
    span = states * 4,
    /* this is required for doing some pre-computations that help to save 
       numerical operations. What we are actually computing here are additional lookup tables 
       for each possible state a certain data-type can assume.
       for DNA with ambuguity coding this is 15, for proteins this is 22 or 23, since there 
       also exist one or two amibguity codes for protein data.
       Essentially this is very similar to the tip vectors which we also use as lookup tables */
    precomputeLength = maxStateValue * span;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	/* allocate pre-compute memory space */

	double 
	  *umpX1 = (double*)malloc(sizeof(double) * precomputeLength),
	  *umpX2 = (double*)malloc(sizeof(double) * precomputeLength);

	/* multiply all possible tip state vectors with the respective P-matrices 
	 */

	for(i = 0; i < maxStateValue; i++)
	  {
	    v = &(tipVector[states * i]);

	    for(k = 0; k < span; k++)
	      {

		umpX1[span * i + k] = 0.0;
		umpX2[span * i + k] = 0.0;

		for(l = 0; l < states; l++)
		  {
		    umpX1[span * i + k] +=  v[l] *  left[k * states + l];
		    umpX2[span * i + k] +=  v[l] * right[k * states + l];
		  }

	      }
	  }

	for(i = 0; i < n; i++)
	  {
	    /* access the precomputed arrays (pre-computed multiplication of conditional with the tip state) 
	     */

	    uX1 = &umpX1[span * tipX1[i]];
	    uX2 = &umpX2[span * tipX2[i]];

	    /* loop over discrete GAMMA rates */

	    for(j = 0; j < 4; j++)
	      {
		/* the rest is the same as for CAT */
		v = &x3[i * span + j * states];

		for(k = 0; k < states; k++)
		  v[k] = 0.0;

		for(k = 0; k < states; k++)
		  {		   
		    x1px2 = uX1[j * states + k] * uX2[j * states + k];
		   
		    for(l = 0; l < states; l++)		      					
		      v[l] += x1px2 * extEV[states * k + l];		     
		  }

	      }	   
	  }
	
	/* free precomputed vectors */

	free(umpX1);
	free(umpX2);
      }
      break;
    case TIP_INNER:
      {
	/* we do analogous pre-computations as above, with the only difference that we now do them 
	   only for one tip vector */

	double 
	  *umpX1 = (double*)malloc(sizeof(double) * precomputeLength),
	  *ump_x2 = (double*)malloc(sizeof(double) * states);

	/* precompute P and left tip vector product */

	for(i = 0; i < maxStateValue; i++)
	  {
	    v = &(tipVector[states * i]);

	    for(k = 0; k < span; k++)
	      {
  
		umpX1[span * i + k] = 0.0;

		for(l = 0; l < states; l++)
		  umpX1[span * i + k] +=  v[l] * left[k * states + l];


	      }
	  }

	for (i = 0; i < n; i++)
	  {
	    /* access pre-computed value based on the raw sequence data tipX1 that is used as an index */

	    uX1 = &umpX1[span * tipX1[i]];

	    /* loop over discrete GAMMA rates */

	    for(k = 0; k < 4; k++)
	      {
		v = &(x2[span * i + k * states]);

		for(l = 0; l < states; l++)
		  {
		    ump_x2[l] = 0.0;

		    for(j = 0; j < states; j++)
		      ump_x2[l] += v[j] * right[k * statesSquare + l * states + j];
		  }

		v = &(x3[span * i + states * k]);

		for(l = 0; l < states; l++)
		  v[l] = 0;

		for(l = 0; l < states; l++)
		  {
		    x1px2 = uX1[k * states + l]  * ump_x2[l];
		    for(j = 0; j < states; j++)
		      v[j] += x1px2 * extEV[l * states  + j];
		  }
	      }
	   
	    /* also do numerical scaling as above. Note that here we need to scale 
	       4 * 4 values for DNA or 4 * 20 values for protein data.
	       If they are ALL smaller than our threshold, we scale. Note that,
	       this can cause numerical problems with GAMMA, if the values generated 
	       by the four discrete GAMMA rates are too different.

	       For details, see: 
	       
	       F. Izquierdo-Carrasco, S.A. Smith, A. Stamatakis: "Algorithms, Data Structures, and Numerics for Likelihood-based Phylogenetic Inference of Huge Trees"

	    */
	    

	    v = &x3[span * i];
	    scale = 1;
	    for(l = 0; scale && (l < span); l++)
	      scale = (ABS(v[l]) <  minlikelihood);


	    if (scale)
	      {
		for(l = 0; l < span; l++)
		  v[l] *= twotothe256;
	
		addScale += wgt[i];		    
	      }
	  }

	free(umpX1);
	free(ump_x2);
      }
      break;
    case INNER_INNER:

      /* same as above, without pre-computations */

      for (i = 0; i < n; i++)
       {
	 for(k = 0; k < 4; k++)
	   {
	     vl = &(x1[span * i + states * k]);
	     vr = &(x2[span * i + states * k]);
	     v =  &(x3[span * i + states * k]);


	     for(l = 0; l < states; l++)
	       v[l] = 0;


	     for(l = 0; l < states; l++)
	       {		 

		 al = 0.0;
		 ar = 0.0;

		 for(j = 0; j < states; j++)
		   {
		     al += vl[j] * left[k * statesSquare + l * states + j];
		     ar += vr[j] * right[k * statesSquare + l * states + j];
		   }

		 x1px2 = al * ar;

		 for(j = 0; j < states; j++)
		   v[j] += x1px2 * extEV[states * l + j];

	       }
	   }
	 
	 v = &(x3[span * i]);
	 scale = 1;
	 for(l = 0; scale && (l < span); l++)
	   scale = ((ABS(v[l]) <  minlikelihood));

	 if(scale)
	   {  
	     for(l = 0; l < span; l++)
	       v[l] *= twotothe256;
	     
	     addScale += wgt[i];	    	  
	   }
       }
      break;
    default:
      assert(0);
    }

  /* as above, increment the global counter that counts scaling multiplications by the scaling multiplications 
     carried out for computing the likelihood array at node p */

  *scalerIncrement = addScale;
}

#endif


    
/* The function below computes partial traversals only down to the point/node in the tree where the 
   conditional likelihhod vector summarizing a subtree is already oriented in the correct direction */

void computeTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches, boolean partialTraversal)
{
  /* if it's a tip we don't do anything */

  if(isTip(p->number, maxTips))
    return;

  {
    int 
      i;
    
    /* get the left and right descendants */

    nodeptr 
      q = p->next->back,
      r = p->next->next->back;   

    /* if the left and right children are tips there is not that much to do */

    if(isTip(r->number, maxTips) && isTip(q->number, maxTips))
      {
	/* fix the orientation of p->x */
	
	if (! p->x)
	  getxnode(p);	
	assert(p->x);
	  
	/* add the current node triplet p,q,r to the traversal descriptor */

	ti[*counter].tipCase = TIP_TIP;
	ti[*counter].pNumber = p->number;
	ti[*counter].qNumber = q->number;
	ti[*counter].rNumber = r->number;

	/* copy branches to traversal descriptor */

	for(i = 0; i < numBranches; i++)
	  {	    
	    ti[*counter].qz[i] = q->z[i];
	    ti[*counter].rz[i] = r->z[i];
	  }

	/* increment length counter */

	*counter = *counter + 1;
      }
    else
      {
	/* if either r or q are tips, flip them to make sure that the tip data is stored 
	   for q */

	if(isTip(r->number, maxTips) || isTip(q->number, maxTips))
	  {	    
	    if(isTip(r->number, maxTips))
	      {
		nodeptr 
		  tmp = r;
		r = q;
		q = tmp;
	      }
	   
	    /* if the orientation of the liklihood vector at r is not correct we need to re-compute it 
	       and descend into its subtree to figure out if there are more vrctors in there to re-compute and 
	       re-orient */

	    if(!r->x || !partialTraversal)
	      computeTraversalInfo(r, ti, counter, maxTips, numBranches, partialTraversal);
	    if(! p->x)
	      getxnode(p);	 
	    
	    /* make sure that everything is consistent now */

	    assert(p->x && r->x);

	    /* store data for p, q, r in the traversal descriptor */

	    ti[*counter].tipCase = TIP_INNER;
	    ti[*counter].pNumber = p->number;
	    ti[*counter].qNumber = q->number;
	    ti[*counter].rNumber = r->number;

	    for(i = 0; i < numBranches; i++)
	      {	
		ti[*counter].qz[i] = q->z[i];
		ti[*counter].rz[i] = r->z[i];
	      }

	    *counter = *counter + 1;
	  }
	else
	  {
	    /* same as above, only now q and r are inner nodes. Hence if they are not 
	       oriented correctly they will need to be recomputed and we need to descend into the 
	       respective subtrees to check if everything is consistent in there, potentially expanding 
	       the traversal descriptor */
	   
	    if(! q->x || !partialTraversal)
	      computeTraversalInfo(q, ti, counter, maxTips, numBranches, partialTraversal);
	    if(! r->x || !partialTraversal)
	      computeTraversalInfo(r, ti, counter, maxTips, numBranches, partialTraversal);
	    if(! p->x)
	      getxnode(p);
	     
	    /* check that the vector orientations are consistent now */

	    assert(p->x && r->x && q->x);

	    ti[*counter].tipCase = INNER_INNER;
	    ti[*counter].pNumber = p->number;
	    ti[*counter].qNumber = q->number;
	    ti[*counter].rNumber = r->number;

	    for(i = 0; i < numBranches; i++)
	      {	
		ti[*counter].qz[i] = q->z[i];
		ti[*counter].rz[i] = r->z[i];
	      }

	    *counter = *counter + 1;
	  }
      }
  }
}

#ifdef _OPTIMIZED_FUNCTIONS

//mth generic vectorized N-state function for CLV updates

static void newviewGTRGAMMA_NSTATES(int tipCase,
				    double *x1, double *x2, double *x3, double *extEV, double *tipVector,
				    unsigned char *tipX1, unsigned char *tipX2,
				    size_t n, double *left, double *right, int *wgt, int *scalerIncrement, const size_t numberOfAllCharacters, const size_t numberOfStates, 
				    const size_t gammaRates);

#endif

/* below are the optimized unrolled, and vectorized versions of the above generi cfunctions 
   for computing the conditional likelihood at p given child nodes q and r. The actual implementation is located at the end/bottom of this 
   file.
*/

#if (defined(_OPTIMIZED_FUNCTIONS) && !defined(__AVX))



static void newviewGTRGAMMAPROT_LG4(int tipCase,
				    double *x1, double *x2, double *x3, double *extEV[4], double *tipVector[4],
				    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				    size_t n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling);

static void newviewGTRGAMMA_GAPPED_SAVE(int tipCase,
					double *x1_start, double *x2_start, double *x3_start,
					double *EV, double *tipVector,
					unsigned char *tipX1, unsigned char *tipX2,
					const size_t n, double *left, double *right, int *wgt, int *scalerIncrement, 
					unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap, 
					double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn);

static void newviewGTRGAMMA(int tipCase,
			    double *x1_start, double *x2_start, double *x3_start,
			    double *EV, double *tipVector,
			    unsigned char *tipX1, unsigned char *tipX2,
			    const size_t n, double *left, double *right, int *wgt, int *scalerIncrement
			    );

static void newviewGTRCAT( int tipCase,  double *EV,  int *cptr,
			   double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
			   unsigned char *tipX1, unsigned char *tipX2,
			   size_t n,  double *left, double *right, int *wgt, int *scalerIncrement);


static void newviewGTRCAT_SAVE( int tipCase,  double *EV,  int *cptr,
				double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
				unsigned char *tipX1, unsigned char *tipX2,
				size_t n,  double *left, double *right, int *wgt, int *scalerIncrement,
				unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
				double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats);

static void newviewGTRGAMMAPROT_GAPPED_SAVE(int tipCase,
					    double *x1, double *x2, double *x3, double *extEV, double *tipVector,
					    unsigned char *tipX1, unsigned char *tipX2,
					    size_t n, double *left, double *right, int *wgt, int *scalerIncrement, 
					    unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,  
					    double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn
					    );



static void newviewGTRGAMMAPROT(int tipCase,
				double *x1, double *x2, double *x3, double *extEV, double *tipVector,
				unsigned char *tipX1, unsigned char *tipX2,
				size_t n, double *left, double *right, int *wgt, int *scalerIncrement);
static void newviewGTRCATPROT(int tipCase, double *extEV,
			      int *cptr,
			      double *x1, double *x2, double *x3, double *tipVector,
			      unsigned char *tipX1, unsigned char *tipX2,
			      size_t n, double *left, double *right, int *wgt, int *scalerIncrement );

static void newviewGTRCATPROT_SAVE(int tipCase, double *extEV,
				   int *cptr,
				   double *x1, double *x2, double *x3, double *tipVector,
				   unsigned char *tipX1, unsigned char *tipX2,
				   size_t n, double *left, double *right, int *wgt, int *scalerIncrement,
				   unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
				   double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats);

#endif

#ifdef _OPTIMIZED_FUNCTIONS

static void newviewGTRCAT_BINARY( int tipCase,  double *EV,  int *cptr,
                                  double *x1_start,  double *x2_start,  double *x3_start,  double *tipVector,
                                  int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                  size_t n,  double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling);

static void newviewGTRGAMMA_BINARY(int tipCase,
				   double *x1_start, double *x2_start, double *x3_start,
				   double *EV, double *tipVector,
				   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   const size_t n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling
				   );

#endif

boolean isGap(unsigned int *x, size_t pos)
{
  return (boolean)((x[pos / 32] & mask32[pos % 32]));
}

boolean noGap(unsigned int *x, size_t pos)
{
  return (boolean)(!(x[pos / 32] & mask32[pos % 32]));
}

/* now this is the function that just iterates over the length of the traversal descriptor and 
   just computes the conditional likelihhod arrays in the order given by the descriptor.
   So in a sense, this function has no clue that there is any tree-like structure 
   in the traversal descriptor, it just operates on an array of structs of given length */ 


extern const char inverseMeaningDNA[16]; 

void newviewIterative (tree *tr, int startIndex)
{
  traversalInfo 
    *ti   = tr->td[0].ti;

  int 
    i;

    /* loop over traversal descriptor length. Note that on average we only re-compute the conditionals on 3 -4
       nodes in RAxML */

  for(i = startIndex; i < tr->td[0].count; i++)
    {
      traversalInfo 
	*tInfo = &ti[i];
      
      {
	int 
	  model;
      
#ifdef _USE_OMP
#pragma omp parallel for
#endif
	for(model = 0; model < tr->NumberOfModels; model++)
	  {
	    /* check if this partition has to be processed now - otherwise no need to compute P matrix */
	    if(!tr->td[0].executeModel[model] || tr->partitionData[model].width == 0)
	      continue;
	    
	    size_t
	      categories,
	      states = (size_t)tr->partitionData[model].states;
	    
	    double
	      qz,
	      rz,
	      *rateCategories,
	      *left = tr->partitionData[model].left,
	      *right = tr->partitionData[model].right,
	      plain[1] = {1.0};
	    
	    /* figure out what kind of rate heterogeneity approach we are using */
	    switch(tr->rateHetModel)
	      {
	      case CAT:
		rateCategories = tr->partitionData[model].perSiteRates;
		categories = (size_t)tr->partitionData[model].numberOfCategories;
		break;
	      case GAMMA:
		rateCategories = tr->partitionData[model].gammaRates;
		categories = 4;
		break;
	      case PLAIN:
		rateCategories = plain;
		categories = 1;
		break;
	      default:
		assert(0);
	      }
	    
	    /* if we use per-partition branch length optimization
	       get the branch length of partition model and take the log otherwise
	       use the joint branch length among all partitions that is always stored
	       at index [0] */
	    if(tr->numBranches > 1)
	      {
		qz = tInfo->qz[model];
		rz = tInfo->rz[model];
	      }
	    else
	      {
		qz = tInfo->qz[0];
		rz = tInfo->rz[0];
	      }
	    
	    qz = (qz > zmin) ? log(qz) : log(zmin);
	    rz = (rz > zmin) ? log(rz) : log(zmin);
	    
	    /* compute the left and right P matrices */
#ifdef __MIC_NATIVE
	    switch (tr->partitionData[model].states)
	      {
	      case 2: /* BINARY data */
		assert(0 && "Binary data model is not implemented on Intel MIC");
		break;
	      case 4: /* DNA data */
		{
		  makeP_DNA_MIC(qz, rz, rateCategories,   tr->partitionData[model].EI,
				tr->partitionData[model].EIGN, categories,
				left, right, tr->saveMemory, tr->maxCategories);
		  
		  precomputeTips_DNA_MIC(tInfo->tipCase, tr->partitionData[model].tipVector,
					 left, right,
					 tr->partitionData[model].mic_umpLeft, tr->partitionData[model].mic_umpRight,
					 categories);
		} 
		break;
	      case 20: /* AA data */
		{
		  if(tr->partitionData[model].protModels == LG4M || tr->partitionData[model].protModels == LG4X)
		    {
		      makeP_PROT_LG4_MIC(qz, rz, tr->partitionData[model].gammaRates,
					 tr->partitionData[model].EI_LG4, tr->partitionData[model].EIGN_LG4,
					 4, left, right);
		      
		      precomputeTips_PROT_LG4_MIC(tInfo->tipCase, tr->partitionData[model].tipVector_LG4,
						  left, right,
						  tr->partitionData[model].mic_umpLeft, tr->partitionData[model].mic_umpRight,
						  categories);
		    }
		  else
		    {
		      makeP_PROT_MIC(qz, rz, rateCategories, tr->partitionData[model].EI,
				     tr->partitionData[model].EIGN, categories,
				     left, right, tr->saveMemory, tr->maxCategories);
		      
		      precomputeTips_PROT_MIC(tInfo->tipCase, tr->partitionData[model].tipVector,
					      left, right,
					      tr->partitionData[model].mic_umpLeft, tr->partitionData[model].mic_umpRight,
					      categories);
		    }
		} 
		break;
	      default:
		assert(0);
	      }
#else
	    if(tr->partitionData[model].protModels == LG4M || tr->partitionData[model].protModels == LG4X)
	      makeP_FlexLG4(qz, rz, tr->partitionData[model].gammaRates,
			    tr->partitionData[model].EI_LG4,
			    tr->partitionData[model].EIGN_LG4,
			    4, left, right, 20);
	    else
	      makeP(qz, rz, rateCategories,   tr->partitionData[model].EI,
		    tr->partitionData[model].EIGN, categories,
		    left, right, tr->saveMemory, tr->maxCategories, states);
#endif
	  } // for model
      }

      /* now loop over all partitions for nodes p, q, and r of the current traversal vector entry */
#ifdef _USE_OMP
#pragma omp parallel
#endif
      {
	int
	  m,
	  model,
	  maxModel;
	
#ifdef _USE_OMP
	maxModel = tr->maxModelsPerThread;
#else
	maxModel = tr->NumberOfModels;
#endif

	for(m = 0; m < maxModel; m++)
	  {
	    size_t
	      width  = 0,
	      offset = 0;
	    
	    double
	      *left     = (double*)NULL,
	      *right    = (double*)NULL;
	    
	    unsigned int
	      *globalScaler = (unsigned int*)NULL;

#ifdef _USE_OMP
	    int
	      tid = omp_get_thread_num();

	    /* check if this thread should process this partition */
	    Assign* 
	      pAss = tr->threadPartAssigns[tid * tr->maxModelsPerThread + m];

	    if(pAss)
	      {
		assert(tid == pAss->procId);
		
		model  = pAss->partitionId;
		width  = pAss->width;
		offset = pAss->offset;
		
		left  = tr->partitionData[model].left;
		right = tr->partitionData[model].right;
		globalScaler = tr->partitionData[model].threadGlobalScaler[tid];
	      }
	    else
	      break;
#else
	    model = m;	    

	    /* number of sites in this partition */
	    width  = (size_t)tr->partitionData[model].width;
	    offset = 0;

	    /* set the pointers to the left and right P matrices to the pre-allocated memory space for storing them */
	    
	    left  = tr->partitionData[model].left;
	    right = tr->partitionData[model].right;
	    globalScaler = tr->partitionData[model].globalScaler;
#endif

	    /* this conditional statement is exactly identical to what we do in evaluateIterative */
	    if(tr->td[0].executeModel[model] && width > 0)
	      {	      
		double
		  *x1_start = (double*)NULL,
		  *x2_start = (double*)NULL,		 
		  *x3_start = (double*)NULL, //tr->partitionData[model].xVector[tInfo->pNumber - tr->mxtips - 1],
		  *x1_gapColumn = (double*)NULL,
		  *x2_gapColumn = (double*)NULL,
		  *x3_gapColumn = (double*)NULL;

		int
		  genericTipCase  = -1,
		  scalerIncrement = 0,
		
		  /* integer wieght vector with pattern compression weights */
		  *wgt = tr->partitionData[model].wgt + offset,

		  /* integer rate category vector (for each pattern, _number_ of PSR category assigned to it, NOT actual rate!) */
		  *rateCategory = tr->partitionData[model].rateCategory + offset;

		unsigned int
		  *x1_gap = (unsigned int*)NULL,
		  *x2_gap = (unsigned int*)NULL,
		  *x3_gap = (unsigned int*)NULL;

		unsigned char
		  *tipX1 = (unsigned char *)NULL,
		  *tipX2 = (unsigned char *)NULL;	
	      
		size_t
		  gapOffset = 0,
		  rateHet = discreteRateCategories(tr->rateHetModel),
		  
		  /* get the number of states in the data stored in partition model */
		  
		  states = (size_t)tr->partitionData[model].states,	
		  
		  /* span for single alignment site (in doubles!) */
		  span = rateHet * states,
		  x_offset = offset * (size_t)span,
		  
		  
		  /* get the length of the current likelihood array stored at node p. This is 
		     important mainly for the SEV-based memory saving option described in here:
		     
		     F. Izquierdo-Carrasco, S.A. Smith, A. Stamatakis: "Algorithms, Data Structures, and Numerics for Likelihood-based Phylogenetic Inference of Huge Trees".
		     
		     So tr->partitionData[model].xSpaceVector[i] provides the length of the allocated conditional array of partition model 
		     and node i 
		  */
		
		  availableLength = tr->partitionData[model].xSpaceVector[(tInfo->pNumber - tr->mxtips - 1)],
		  requiredLength = 0;	     

		x3_start = tr->partitionData[model].xVector[tInfo->pNumber - tr->mxtips - 1] + x_offset;

		/* memory saving stuff, not important right now, but if you are interested ask Fernando */
		if(tr->saveMemory)
		  {
		    size_t
		      j,
		      setBits = 0;		  
		    
		    gapOffset = states * (size_t)getUndetermined(tr->partitionData[model].dataType);
		    
		    x1_gap = &(tr->partitionData[model].gapVector[tInfo->qNumber * tr->partitionData[model].gapVectorLength]);
		    x2_gap = &(tr->partitionData[model].gapVector[tInfo->rNumber * tr->partitionData[model].gapVectorLength]);
		    x3_gap = &(tr->partitionData[model].gapVector[tInfo->pNumber * tr->partitionData[model].gapVectorLength]);		      		  
		    
		    for(j = 0; j < (size_t)tr->partitionData[model].gapVectorLength; j++)
		      {		     
			x3_gap[j] = x1_gap[j] & x2_gap[j];
			setBits += (size_t)__builtin_popcount(x3_gap[j]);		      
		      }
		    
		    requiredLength = (width - setBits)  * rateHet * states * sizeof(double);		
		  }
		else
		  /* if we are not trying to save memory the space required to store an inner likelihood array 
		     is the number of sites in the partition times the number of states of the data type in the partition 
		     times the number of discrete GAMMA rates (1 for CAT essentially) times 8 bytes */
		  requiredLength  =  width * rateHet * states * sizeof(double);
		
		/* Initially, even when not using memory saving no space is allocated for inner likelihood arrats hence 
		   availableLength will be zero at the very first time we traverse the tree.
		   Hence we need to allocate something here */
#ifndef _USE_OMP
		if(requiredLength != availableLength)
		  {
		    /* if there is a vector of incorrect length assigned here i.e., x3 != NULL we must free
		       it first */
		    if(x3_start)
		      free(x3_start);
		    
		    /* allocate memory: note that here we use a byte-boundary aligned malloc, because we need the vectors
		       to be aligned at 16 BYTE (SSE3) or 32 BYTE (AVX) boundaries! */
		    
		    x3_start = (double*)malloc_aligned(requiredLength);
		    
		    /* update the data structures for consistent bookkeeping */
		    tr->partitionData[model].xVector[tInfo->pNumber - tr->mxtips - 1] = x3_start;
		    tr->partitionData[model].xSpaceVector[(tInfo->pNumber - tr->mxtips - 1)] = requiredLength;
		  }
#endif

		/* now just set the pointers for data accesses in the newview() implementations above to the corresponding values 
		   according to the tip case */
		
		switch(tInfo->tipCase)
		  {
		  case TIP_TIP:		  
		    if(isPomo(tr->partitionData[model].dataType))
		      {
			//mth add appropriate offset for MIC version, note that, we just count the number of double entries irrespctive of the number of rate cats!
			assert(offset == 0 && x_offset == 0);
			x1_start = tr->partitionData[model].xTipVector[tInfo->qNumber];						
			x2_start = tr->partitionData[model].xTipVector[tInfo->rNumber];

			genericTipCase = TIP_TIP_CLV;
		      }
		    else
		      {
			tipX1    = tr->partitionData[model].yVector[tInfo->qNumber] + offset;
			tipX2    = tr->partitionData[model].yVector[tInfo->rNumber] + offset;

			genericTipCase = TIP_TIP;
		      }
		    
		    if(tr->saveMemory)
		      {
			assert(tInfo->pNumber - tr->mxtips - 1 >= 0);

			x1_gapColumn   = &(tr->partitionData[model].tipVector[gapOffset]);
			x2_gapColumn   = &(tr->partitionData[model].tipVector[gapOffset]);		    
			x3_gapColumn   = &tr->partitionData[model].gapColumn[((size_t)tInfo->pNumber - (size_t)tr->mxtips - 1) * states * rateHet];		    
		      }
		    
		    break;
		  case TIP_INNER:
		    if(isPomo(tr->partitionData[model].dataType))
		      {	
			//mth add appropriate offset for MIC version, note that, we just count the number of double entries irrespctive of the number of rate cats!
			assert(offset == 0 && x_offset == 0);
			x1_start = tr->partitionData[model].xTipVector[tInfo->qNumber];
			
			genericTipCase = TIP_INNER_CLV;
		      }
		    else
		      {
			tipX1    =  tr->partitionData[model].yVector[tInfo->qNumber] + offset;
			
			genericTipCase = TIP_INNER;
		      }
		    
		    x2_start = tr->partitionData[model].xVector[tInfo->rNumber - tr->mxtips - 1] + x_offset;
		    
		    if(tr->saveMemory)
		      {	
			assert(tInfo->rNumber - tr->mxtips - 1 >= 0 &&
			       tInfo->pNumber - tr->mxtips - 1 >= 0);
			       

			x1_gapColumn   = &(tr->partitionData[model].tipVector[gapOffset]);	     
			x2_gapColumn   = &tr->partitionData[model].gapColumn[((size_t)tInfo->rNumber - (size_t)tr->mxtips - 1) * states * rateHet];
			x3_gapColumn   = &tr->partitionData[model].gapColumn[((size_t)tInfo->pNumber - (size_t)tr->mxtips - 1) * states * rateHet];
		      }
		    
		    break;
		  case INNER_INNER:	 
		    if(isPomo(tr->partitionData[model].dataType))	 		 
		      genericTipCase = INNER_INNER;
		    
		    x1_start       = tr->partitionData[model].xVector[tInfo->qNumber - tr->mxtips - 1] + x_offset;
		    x2_start       = tr->partitionData[model].xVector[tInfo->rNumber - tr->mxtips - 1] + x_offset;
		    
		    assert(tInfo->rNumber - tr->mxtips - 1 >= 0 &&
			   tInfo->pNumber - tr->mxtips - 1 >= 0 && 
			   tInfo->qNumber - tr->mxtips - 1 >= 0);

		    if(tr->saveMemory)
		      {
			x1_gapColumn   = &tr->partitionData[model].gapColumn[((size_t)tInfo->qNumber - (size_t)tr->mxtips - 1) * states * rateHet];
			x2_gapColumn   = &tr->partitionData[model].gapColumn[((size_t)tInfo->rNumber - (size_t)tr->mxtips - 1) * states * rateHet];
			x3_gapColumn   = &tr->partitionData[model].gapColumn[((size_t)tInfo->pNumber - (size_t)tr->mxtips - 1) * states * rateHet];
		      }
		    
		    break;
		  default:
		    assert(0);
		  }
		
#ifndef _OPTIMIZED_FUNCTIONS

	      /* memory saving not implemented */

	      assert(!tr->saveMemory);

	      assert(tr->rateHetModel != PLAIN);

	      /* figure out if we need to compute the CAT or GAMMA model of rate heterogeneity */

	      if(tr->rateHetModel == CAT)
		newviewCAT_FLEX(tInfo->tipCase,  tr->partitionData[model].EV, rateCategory,
				x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
				tipX1, tipX2,
				width, left, right, wgt, &scalerIncrement, states);
	      else
		newviewGAMMA_FLEX(tInfo->tipCase,
				  x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
				  tipX1, tipX2,
				  width, left, right, wgt, &scalerIncrement, states, getUndetermined(tr->partitionData[model].dataType) + 1);

#else
	      /* dedicated highly optimized functions. Analogously to the functions in evaluateGeneric() 
		 we also siwtch over the state number */

	      switch(states)
		{		
		case 2:
#ifdef __MIC_NATIVE
 	      assert(0 && "Binary data model is not implemented on Intel MIC");
#else
		  assert(!tr->saveMemory);

		  assert(tr->rateHetModel != PLAIN);
		  
		  if(tr->rateHetModel == CAT)
		    newviewGTRCAT_BINARY(tInfo->tipCase,  tr->partitionData[model].EV,  rateCategory,
					 x1_start,  x2_start,  x3_start, tr->partitionData[model].tipVector,
					 (int*)NULL, tipX1, tipX2,
					 width, left, right, wgt, &scalerIncrement, TRUE);
		  else
		    newviewGTRGAMMA_BINARY(tInfo->tipCase,
					   x1_start, x2_start, x3_start,
					   tr->partitionData[model].EV, tr->partitionData[model].tipVector,
					   (int *)NULL, tipX1, tipX2,
					   width, left, right, wgt, &scalerIncrement, TRUE);		 
#endif
		  break;
		case 4:	/* DNA */
		  assert(tr->rateHetModel != PLAIN);

		  if(tr->rateHetModel == CAT)
		    {		    		     
		      if(tr->saveMemory)
#ifdef __MIC_NATIVE
		     assert(0 && "Neither CAT model of rate heterogeneity nor memory saving are implemented on Intel MIC");
#elif __AVX
			newviewGTRCAT_AVX_GAPPED_SAVE(tInfo->tipCase,  tr->partitionData[model].EV, rateCategory,
						      x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
						      (int*)NULL, tipX1, tipX2,
						      width, left, right, wgt, &scalerIncrement, TRUE, x1_gap, x2_gap, x3_gap,
						      x1_gapColumn, x2_gapColumn, x3_gapColumn, tr->maxCategories);
#else
			newviewGTRCAT_SAVE(tInfo->tipCase,  tr->partitionData[model].EV, rateCategory,
					   x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					   tipX1, tipX2,
					   width, left, right, wgt, &scalerIncrement, x1_gap, x2_gap, x3_gap,
					   x1_gapColumn, x2_gapColumn, x3_gapColumn, tr->maxCategories);
#endif
		      else
#ifdef __MIC_NATIVE
		     assert(0 && "CAT model of rate heterogeneity is not implemented on Intel MIC");
#elif __AVX
			newviewGTRCAT_AVX(tInfo->tipCase,  tr->partitionData[model].EV, rateCategory,
					  x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					  tipX1, tipX2,
					  width, left, right, wgt, &scalerIncrement);
#else
			newviewGTRCAT(tInfo->tipCase,  tr->partitionData[model].EV, rateCategory,
				      x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
				      tipX1, tipX2,
				      width, left, right, wgt, &scalerIncrement);
#endif
		    }
		  else
		    {
		      
		       
		       if(tr->saveMemory)
#ifdef __MIC_NATIVE
		     assert(0 && "Memory saving is not implemented on Intel MIC");
#elif __AVX
			 newviewGTRGAMMA_AVX_GAPPED_SAVE(tInfo->tipCase,
							 x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector, (int*)NULL,
							 tipX1, tipX2,
							 width, left, right, wgt, &scalerIncrement, TRUE,
							 x1_gap, x2_gap, x3_gap, 
							 x1_gapColumn, x2_gapColumn, x3_gapColumn);
#else
		       newviewGTRGAMMA_GAPPED_SAVE(tInfo->tipCase,
						   x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
						   tipX1, tipX2,
						   width, left, right, wgt, &scalerIncrement, 
						   x1_gap, x2_gap, x3_gap, 
						   x1_gapColumn, x2_gapColumn, x3_gapColumn);
#endif
		       else
#ifdef __MIC_NATIVE
			 newviewGTRGAMMA_MIC(tInfo->tipCase,
				  x1_start, x2_start, x3_start, tr->partitionData[model].mic_EV, tr->partitionData[model].tipVector,
				  tipX1, tipX2,
				  width, left, right, wgt, &scalerIncrement,
				  tr->partitionData[model].mic_umpLeft, tr->partitionData[model].mic_umpRight);
#elif __AVX
			 newviewGTRGAMMA_AVX(tInfo->tipCase,
					     x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
					     tipX1, tipX2,
					     width, left, right, wgt, &scalerIncrement);
#else
		       newviewGTRGAMMA(tInfo->tipCase,
					 x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
					 tipX1, tipX2,
					 width, left, right, wgt, &scalerIncrement);
#endif
		    }
		
		  break;		    
		case 20: /* proteins */
		  assert(tr->rateHetModel != PLAIN);
		  if(tr->rateHetModel == CAT)
		    {		     
		      if(tr->saveMemory)
			{
#ifdef __MIC_NATIVE
		     assert(0 && "Neither CAT model of rate heterogeneity nor memory saving are implemented on Intel MIC");
#elif __AVX
			  newviewGTRCATPROT_AVX_GAPPED_SAVE(tInfo->tipCase,  tr->partitionData[model].EV, rateCategory,
							    x1_start, x2_start, x3_start, tr->partitionData[model].tipVector, (int*)NULL,
							    tipX1, tipX2, width, left, right, wgt, &scalerIncrement, TRUE, x1_gap, x2_gap, x3_gap,
							    x1_gapColumn, x2_gapColumn, x3_gapColumn, tr->maxCategories);
#else
			  newviewGTRCATPROT_SAVE(tInfo->tipCase,  tr->partitionData[model].EV, rateCategory,
						 x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
						 tipX1, tipX2, width, left, right, wgt, &scalerIncrement, x1_gap, x2_gap, x3_gap,
						 x1_gapColumn, x2_gapColumn, x3_gapColumn, tr->maxCategories);
#endif
			}
		      else
			{			 			
#ifdef __MIC_NATIVE
		     assert(0 && "CAT model of rate heterogeneity is not implemented on Intel MIC");
#elif __AVX
			  newviewGTRCATPROT_AVX(tInfo->tipCase,  tr->partitionData[model].EV, rateCategory,
						x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
						tipX1, tipX2, width, left, right, wgt, &scalerIncrement);
#else
			  newviewGTRCATPROT(tInfo->tipCase,  tr->partitionData[model].EV, rateCategory,
					    x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					    tipX1, tipX2, width, left, right, wgt, &scalerIncrement);			
#endif
			}
		    }
		  else
		    {		    			 			  
		      if(tr->saveMemory)
			{
#ifdef __MIC_NATIVE
		     assert(0 && "Memory saving is not implemented on Intel MIC");
#elif __AVX
			  newviewGTRGAMMAPROT_AVX_GAPPED_SAVE(tInfo->tipCase,
							      x1_start, x2_start, x3_start,
							      tr->partitionData[model].EV,
							      tr->partitionData[model].tipVector, (int*)NULL,
							      tipX1, tipX2,
							      width, left, right, wgt, &scalerIncrement, TRUE,
							      x1_gap, x2_gap, x3_gap,
							      x1_gapColumn, x2_gapColumn, x3_gapColumn);
#else
			  newviewGTRGAMMAPROT_GAPPED_SAVE(tInfo->tipCase,
							  x1_start, x2_start, x3_start,
							  tr->partitionData[model].EV,
							  tr->partitionData[model].tipVector,
							  tipX1, tipX2,
							  width, left, right, wgt, &scalerIncrement,
							  x1_gap, x2_gap, x3_gap,
							  x1_gapColumn, x2_gapColumn, x3_gapColumn);
#endif
			}
		      else
			{
			  if(tr->partitionData[model].protModels == LG4M || tr->partitionData[model].protModels == LG4X)
			    {
#ifdef __MIC_NATIVE
			      newviewGTRGAMMAPROT_LG4_MIC(tInfo->tipCase,
							x1_start, x2_start, x3_start, tr->partitionData[model].mic_EV, tr->partitionData[model].mic_tipVector,
							tipX1, tipX2,
							width, left, right, wgt, &scalerIncrement,
							tr->partitionData[model].mic_umpLeft, tr->partitionData[model].mic_umpRight);
#elif __AVX
			      newviewGTRGAMMAPROT_AVX_LG4(tInfo->tipCase,
							  x1_start, x2_start, x3_start,
							  tr->partitionData[model].EV_LG4,
							  tr->partitionData[model].tipVector_LG4,
							  (int*)NULL, tipX1, tipX2,
							  width, left, right, wgt, &scalerIncrement, TRUE);
#else
			      newviewGTRGAMMAPROT_LG4(tInfo->tipCase,
						      x1_start, x2_start, x3_start,
						      tr->partitionData[model].EV_LG4,
						      tr->partitionData[model].tipVector_LG4,
						      (int*)NULL, tipX1, tipX2,
						      width, left, right, 
						      wgt, &scalerIncrement, TRUE);
#endif			    
			    }
			  else
			    {
#ifdef __MIC_NATIVE
			      newviewGTRGAMMAPROT_MIC(tInfo->tipCase,
							x1_start, x2_start, x3_start, tr->partitionData[model].mic_EV, tr->partitionData[model].mic_tipVector,
							tipX1, tipX2,
							width, left, right, wgt, &scalerIncrement,
							tr->partitionData[model].mic_umpLeft, tr->partitionData[model].mic_umpRight);
#elif __AVX
			     
			      
			      /*newviewGTRGAMMA_NSTATES(tInfo->tipCase,
						       x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
						       tipX1, tipX2,
						       width, left, right, wgt, &scalerIncrement, 23, 20, 4);*/

			       newviewGTRGAMMAPROT_AVX(tInfo->tipCase,
						      x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
						      tipX1, tipX2,
						      width, left, right, wgt, &scalerIncrement);
#else
			       

			      newviewGTRGAMMAPROT(tInfo->tipCase,
						  x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
						  tipX1, tipX2,
						  width, left, right, wgt, &scalerIncrement);			     					      
#endif
			    }
			}
		    }	
		  break;	
		case 16:
		  //mth POMO16
		  assert(!tr->saveMemory);
		  
		  switch(tr->rateHetModel)
		    {
		    case GAMMA:
		      newviewGTRGAMMA_NSTATES(genericTipCase,
					      x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
					      tipX1, tipX2,
					      width, left, right, wgt, &scalerIncrement, 0, 16, 4);
		      break;
		    case PLAIN:
		      newviewGTRGAMMA_NSTATES(genericTipCase,
					      x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
					      tipX1, tipX2,
					      width, left, right, wgt, &scalerIncrement, 0, 16, 1);
		      break;
		    default:
		      assert(0);
		    }
		  break;
		case 64:
		  //mth POMO64
		  assert(!tr->saveMemory); 
		  switch(tr->rateHetModel)
		    {
		    case GAMMA:
		      newviewGTRGAMMA_NSTATES(genericTipCase,
					      x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
					      tipX1, tipX2,
					      width, left, right, wgt, &scalerIncrement, 0, 64, 4);
		      break;
		    case PLAIN:
		      newviewGTRGAMMA_NSTATES(genericTipCase,
					      x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
					      tipX1, tipX2,
					      width, left, right, wgt, &scalerIncrement, 0, 64, 1);
		      break;
		    default:
		      assert(0);
		    }
		  break;
		default:		  
		  assert(0);
		}
#endif

	      /* important step, here we essentiallt recursively compute the number of scaling multiplications 
		 at node p: it's the sum of the number of scaling multiplications already conducted 
		 for computing nodes q and r plus the scaling multiplications done at node p */

	      globalScaler[tInfo->pNumber] =
		globalScaler[tInfo->qNumber] +
		globalScaler[tInfo->rNumber] +
		(unsigned int)scalerIncrement;

	      /* check that we are not getting an integer overflow ! */

	      assert(globalScaler[tInfo->pNumber] < INT_MAX);
	    }	
	} // for model
    }  // omp parallel block
  }  // for traversal
}


/* here is the generic function that could be called from the user program 
   it re-computes the vector at node p (regardless of whether it's orientation is 
   correct and then it also re-computes reciursively the likelihood arrays 
   in the subtrees of p as needed and if needed */

void newviewGeneric (tree *tr, nodeptr p, boolean masked)
{  
  /* if it's a tip there is nothing to do */

  if(isTip(p->number, tr->mxtips))
    return;
  
  /* the first entry of the traversal descriptor is always reserved for evaluate or branch length optimization calls,
     hence we start filling the array at the second entry with index one. This is not very nice and should be fixed 
     at some point */

  tr->td[0].count = 0;

  /* compute the traversal descriptor */
  computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches, TRUE);

  /* the traversal descriptor has been recomputed -> not sure if it really always changes, something to 
     optimize in the future */
  tr->td[0].traversalHasChanged = TRUE;
  
  /* We do a masked newview, i.e., do not execute newvies for each partition, when for example 
     doing a branch length optimization on the entire tree when branches are estimated on a per partition basis.

     you may imagine that for partition 5 the branch length optimization has already converged whereas 
     for partition 6 we still need to go over the tree again.

     This is explained in more detail in:

     A. Stamatakis, M. Ott: "Load Balance in the Phylogenetic Likelihood Kernel". Proceedings of ICPP 2009

     The external boolean array tr->partitionConverged[] contains exactly that information and is copied 
     to executeModel and subsequently to the executeMask of the traversal descriptor 

  */


  if(masked)
    {
      int model;
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  if(tr->partitionConverged[model])
	    tr->executeModel[model] = FALSE;
	  else
	    tr->executeModel[model] = TRUE;
	}
    }

  /* if there is something to re-compute */

  if(tr->td[0].count > 0)
    {
      /* store execute mask in traversal descriptor */

      storeExecuteMaskInTraversalDescriptor(tr);           
      newviewIterative(tr, 0);
    }

  /* clean up */

  if(masked)
    {
      int model;
      
      for(model = 0; model < tr->NumberOfModels; model++)
	tr->executeModel[model] = TRUE;
    }

  tr->td[0].traversalHasChanged = FALSE;
}


/* optimized function implementations */

#if (defined(_OPTIMIZED_FUNCTIONS) && !defined(__AVX))

static void newviewGTRGAMMA_GAPPED_SAVE(int tipCase,
					double *x1_start, double *x2_start, double *x3_start,
					double *EV, double *tipVector,
					unsigned char *tipX1, unsigned char *tipX2,
					const size_t n, double *left, double *right, int *wgt, int *scalerIncrement, 
					unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap, 
					double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn)
{
  size_t
    i, 
    j, 
    k, 
    l;

  int
    addScale = 0, 
    scaleGap = 0;
  
  double
    *x1,
    *x2,
    *x3,
    *x1_ptr = x1_start,
    *x2_ptr = x2_start,       
    max,
    maxima[2] __attribute__ ((aligned (BYTE_ALIGNMENT))),        
    EV_t[16] __attribute__ ((aligned (BYTE_ALIGNMENT)));      
    
  __m128d 
    values[8],
    EVV[8];  

  for(k = 0; k < 4; k++)
    for (l=0; l < 4; l++)
      EV_t[4 * l + k] = EV[4 * k + l];

  for(k = 0; k < 8; k++)
    EVV[k] = _mm_load_pd(&EV_t[k * 2]);      
 
  

  switch(tipCase)
    {
    case TIP_TIP:
      {
	double *uX1, umpX1[256] __attribute__ ((aligned (BYTE_ALIGNMENT))), *uX2, umpX2[256] __attribute__ ((aligned (BYTE_ALIGNMENT)));


	for (i = 1; i < 16; i++)
	  {	    
	    __m128d x1_1 = _mm_load_pd(&(tipVector[i*4]));
	    __m128d x1_2 = _mm_load_pd(&(tipVector[i*4 + 2]));	   

	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{			 	 
		  __m128d left1 = _mm_load_pd(&left[j*16 + k*4]);
		  __m128d left2 = _mm_load_pd(&left[j*16 + k*4 + 2]);
		  
		  __m128d acc = _mm_setzero_pd();

		  acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
		  acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));
		  		  
		  acc = _mm_hadd_pd(acc, acc);
		  _mm_storel_pd(&umpX1[i*16 + j*4 + k], acc);
		}
	  
	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{
		  __m128d left1 = _mm_load_pd(&right[j*16 + k*4]);
		  __m128d left2 = _mm_load_pd(&right[j*16 + k*4 + 2]);
		  
		  __m128d acc = _mm_setzero_pd();

		  acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
		  acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));
		  		  
		  acc = _mm_hadd_pd(acc, acc);
		  _mm_storel_pd(&umpX2[i*16 + j*4 + k], acc);
		 
		}
	  }   		  
	
	uX1 = &umpX1[240];
	uX2 = &umpX2[240];	   	    	    
	
	for (j = 0; j < 4; j++)
	  {				 		  		  		   
	    __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
	    __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
	    	    
	    __m128d uX2_k0_sse = _mm_load_pd( &uX2[j * 4] );
	    __m128d uX2_k2_sse = _mm_load_pd( &uX2[j * 4 + 2] );
	    
	    __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, uX2_k0_sse );
	    __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, uX2_k2_sse );		    		    		   
	    
	    __m128d EV_t_l0_k0 = EVV[0];
	    __m128d EV_t_l0_k2 = EVV[1];
	    __m128d EV_t_l1_k0 = EVV[2];
	    __m128d EV_t_l1_k2 = EVV[3];
	    __m128d EV_t_l2_k0 = EVV[4];
	    __m128d EV_t_l2_k2 = EVV[5];
	    __m128d EV_t_l3_k0 = EVV[6]; 
	    __m128d EV_t_l3_k2 = EVV[7];
	    
	    EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	    EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	    
	    EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	    EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	    
	    EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	    
	    EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	    EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	    
	    EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	    EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	    EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	    
	    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
	    	  
	    _mm_store_pd( &x3_gapColumn[j * 4 + 0], EV_t_l0_k0 );
	    _mm_store_pd( &x3_gapColumn[j * 4 + 2], EV_t_l2_k0 );	   
	  }  
	
       
	x3 = x3_start;
	
	for (i = 0; i < n; i++)
	  {	    
	    if(!(x3_gap[i / 32] & mask32[i % 32]))	     
	      {
		uX1 = &umpX1[16 * tipX1[i]];
		uX2 = &umpX2[16 * tipX2[i]];	   	    	    		
		
		for (j = 0; j < 4; j++)
		  {				 		  		  		   
		    __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
		    __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
		    
		    
		    __m128d uX2_k0_sse = _mm_load_pd( &uX2[j * 4] );
		    __m128d uX2_k2_sse = _mm_load_pd( &uX2[j * 4 + 2] );
		    
		    
		    //
		    // multiply left * right
		    //
		    
		    __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, uX2_k0_sse );
		    __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, uX2_k2_sse );
		    
		    
		    //
		    // multiply with EV matrix (!?)
		    //
		    
		    __m128d EV_t_l0_k0 = EVV[0];
		    __m128d EV_t_l0_k2 = EVV[1];
		    __m128d EV_t_l1_k0 = EVV[2];
		    __m128d EV_t_l1_k2 = EVV[3];
		    __m128d EV_t_l2_k0 = EVV[4];
		    __m128d EV_t_l2_k2 = EVV[5];
		    __m128d EV_t_l3_k0 = EVV[6]; 
		    __m128d EV_t_l3_k2 = EVV[7];
		    
		    EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
		    EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
		    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
		    
		    EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
		    EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
		    
		    EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
		    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
		    
		    EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
		    EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
		    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
		    
		    EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
		    EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
		    EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
		    
		    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
		    
		    _mm_store_pd( &x3[j * 4 + 0], EV_t_l0_k0 );
		    _mm_store_pd( &x3[j * 4 + 2], EV_t_l2_k0 );
		  }
		
		x3 += 16;
	      }
	  }
      }
      break;
    case TIP_INNER:
      {	
	double 
	  *uX1, 
	  umpX1[256] __attribute__ ((aligned (BYTE_ALIGNMENT)));		 

	for (i = 1; i < 16; i++)
	  {
	    __m128d x1_1 = _mm_load_pd(&(tipVector[i*4]));
	    __m128d x1_2 = _mm_load_pd(&(tipVector[i*4 + 2]));	   

	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m128d left1 = _mm_load_pd(&left[j*16 + k*4]);
		  __m128d left2 = _mm_load_pd(&left[j*16 + k*4 + 2]);
		  
		  __m128d acc = _mm_setzero_pd();

		  acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
		  acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));
		  		  
		  acc = _mm_hadd_pd(acc, acc);
		  _mm_storel_pd(&umpX1[i*16 + j*4 + k], acc);		 
		}
	  }

	{
	  __m128d maxv =_mm_setzero_pd();
	  
	  scaleGap = 0;
	  
	  x2 = x2_gapColumn;			 
	  x3 = x3_gapColumn;
	  
	  uX1 = &umpX1[240];	     
	  
	  for (j = 0; j < 4; j++)
	    {		     		   
	      double *x2_p = &x2[j*4];
	      double *right_k0_p = &right[j*16];
	      double *right_k1_p = &right[j*16 + 1*4];
	      double *right_k2_p = &right[j*16 + 2*4];
	      double *right_k3_p = &right[j*16 + 3*4];
	      __m128d x2_0 = _mm_load_pd( &x2_p[0] );
	      __m128d x2_2 = _mm_load_pd( &x2_p[2] );
	      
	      __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
	      __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
	      __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
	      __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
	      __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
	      __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
	      __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
	      __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );
	      	      
	      right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	      right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	      
	      right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	      right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	      
	      right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	      right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	      right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	      	       
	      right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	      right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	      	       
	      right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	      right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	      
	      right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	      right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	      right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);
	      
	      __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
	      __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
	      
	      __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, right_k0_0 );
	      __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, right_k2_0 );
	      
	      __m128d EV_t_l0_k0 = EVV[0];
	      __m128d EV_t_l0_k2 = EVV[1];
	      __m128d EV_t_l1_k0 = EVV[2];
	      __m128d EV_t_l1_k2 = EVV[3];
	      __m128d EV_t_l2_k0 = EVV[4];
	      __m128d EV_t_l2_k2 = EVV[5];
	      __m128d EV_t_l3_k0 = EVV[6]; 
	      __m128d EV_t_l3_k2 = EVV[7];
	      
	      EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	      EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	      EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	      
	      EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	      EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	      
	      EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	      EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	      
	      EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	      EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	      EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	      
	      EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	      EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	      EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	      
	      EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
	      
	      values[j * 2]     = EV_t_l0_k0;
	      values[j * 2 + 1] = EV_t_l2_k0;		   		   
	      
	      maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
	      maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));		   	     		   
	    }

	  
	  _mm_store_pd(maxima, maxv);
		 
	  max = MAX(maxima[0], maxima[1]);
	  
	  if(max < minlikelihood)
	    {
	      scaleGap = 1;
	      
	      __m128d sv = _mm_set1_pd(twotothe256);
	      
	      _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
	      _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
	      _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
	      _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
	      _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
	      _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
	      _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
	      _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     	      	     
	    }
	  else
	    {
	      _mm_store_pd(&x3[0], values[0]);	   
	      _mm_store_pd(&x3[2], values[1]);
	      _mm_store_pd(&x3[4], values[2]);
	      _mm_store_pd(&x3[6], values[3]);
	      _mm_store_pd(&x3[8], values[4]);	   
	      _mm_store_pd(&x3[10], values[5]);
	      _mm_store_pd(&x3[12], values[6]);
	      _mm_store_pd(&x3[14], values[7]);
	    }
	}		       	
      	
	x3 = x3_start;

	for (i = 0; i < n; i++)
	   {
	     if((x3_gap[i / 32] & mask32[i % 32]))
	       {	       
		 if(scaleGap)
		   {		    
		       addScale += wgt[i];		     
		   }
	       }
	     else
	       {				 
		 __m128d maxv =_mm_setzero_pd();		 
		 
		 if(x2_gap[i / 32] & mask32[i % 32])
		   x2 = x2_gapColumn;
		 else
		   {
		     x2 = x2_ptr;
		     x2_ptr += 16;
		   }
		 		 		 
		 uX1 = &umpX1[16 * tipX1[i]];	     
		 
		 
		 for (j = 0; j < 4; j++)
		   {		     		   
		     double *x2_p = &x2[j*4];
		     double *right_k0_p = &right[j*16];
		     double *right_k1_p = &right[j*16 + 1*4];
		     double *right_k2_p = &right[j*16 + 2*4];
		     double *right_k3_p = &right[j*16 + 3*4];
		     __m128d x2_0 = _mm_load_pd( &x2_p[0] );
		     __m128d x2_2 = _mm_load_pd( &x2_p[2] );
		     
		     __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
		     __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
		     __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
		     __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
		     __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
		     __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
		     __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
		     __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );
		     
		     		     
		     right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
		     right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
		     
		     right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
		     right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
		     
		     right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
		     right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
		     right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
		     
		     
		     right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
		     right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
		     
		     
		     right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
		     right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
		     
		     right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
		     right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
		     right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);
		     
		     {
		       //
		       // load left side from tip vector
		       //
		       
		       __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
		       __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
		       
		       
		       //
		       // multiply left * right
			   //
		       
		       __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, right_k0_0 );
		       __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, right_k2_0 );
		       
		       
		       //
		       // multiply with EV matrix (!?)
		       //		   		  
		       
		       __m128d EV_t_l0_k0 = EVV[0];
		       __m128d EV_t_l0_k2 = EVV[1];
		       __m128d EV_t_l1_k0 = EVV[2];
		       __m128d EV_t_l1_k2 = EVV[3];
		       __m128d EV_t_l2_k0 = EVV[4];
		       __m128d EV_t_l2_k2 = EVV[5];
		       __m128d EV_t_l3_k0 = EVV[6]; 
		       __m128d EV_t_l3_k2 = EVV[7];
		       
		       
		       EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
		       EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
		       EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
		       
		       EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
		       EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
		       
		       EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
		       EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
		       
		       EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
		       EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
		       EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
		       
		       EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
		       EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
		       EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
		       
		       EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
		       
		       values[j * 2]     = EV_t_l0_k0;
		       values[j * 2 + 1] = EV_t_l2_k0;		   		   
			   
		       maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
		       maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));		   
		     }		   
		   }

	     
		 _mm_store_pd(maxima, maxv);
		 
		 max = MAX(maxima[0], maxima[1]);
		 
		 if(max < minlikelihood)
		   {
		     __m128d sv = _mm_set1_pd(twotothe256);
		     
		     _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
		     _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
		     _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
		     _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
		     _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
		     _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
		     _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
		     _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     
		     
		     
		     addScale += wgt[i];
		    
		   }
		 else
		   {
		     _mm_store_pd(&x3[0], values[0]);	   
		     _mm_store_pd(&x3[2], values[1]);
		     _mm_store_pd(&x3[4], values[2]);
		     _mm_store_pd(&x3[6], values[3]);
		     _mm_store_pd(&x3[8], values[4]);	   
		     _mm_store_pd(&x3[10], values[5]);
		     _mm_store_pd(&x3[12], values[6]);
		     _mm_store_pd(&x3[14], values[7]);
		   }		 
		 
		 x3 += 16;
	       }
	   }
      }
      break;
    case INNER_INNER:         
      {
	__m128d maxv =_mm_setzero_pd();
	
	scaleGap = 0;
	
	x1 = x1_gapColumn;	     	    
	x2 = x2_gapColumn;	    
	x3 = x3_gapColumn;
	
	for (j = 0; j < 4; j++)
	  {
	    
	    double *x1_p = &x1[j*4];
	    double *left_k0_p = &left[j*16];
	    double *left_k1_p = &left[j*16 + 1*4];
	    double *left_k2_p = &left[j*16 + 2*4];
	    double *left_k3_p = &left[j*16 + 3*4];
	    
	    __m128d x1_0 = _mm_load_pd( &x1_p[0] );
	    __m128d x1_2 = _mm_load_pd( &x1_p[2] );
	    
	    __m128d left_k0_0 = _mm_load_pd( &left_k0_p[0] );
	    __m128d left_k0_2 = _mm_load_pd( &left_k0_p[2] );
	    __m128d left_k1_0 = _mm_load_pd( &left_k1_p[0] );
	    __m128d left_k1_2 = _mm_load_pd( &left_k1_p[2] );
	    __m128d left_k2_0 = _mm_load_pd( &left_k2_p[0] );
	    __m128d left_k2_2 = _mm_load_pd( &left_k2_p[2] );
	    __m128d left_k3_0 = _mm_load_pd( &left_k3_p[0] );
	    __m128d left_k3_2 = _mm_load_pd( &left_k3_p[2] );
	    
	    left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
	    left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
	    
	    left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
	    left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
	    
	    left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
	    left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
	    left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
	    
	    left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
	    left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
	    
	    left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
	    left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
	    
	    left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
	    left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
	    left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
	    
	    
	    double *x2_p = &x2[j*4];
	    double *right_k0_p = &right[j*16];
	    double *right_k1_p = &right[j*16 + 1*4];
	    double *right_k2_p = &right[j*16 + 2*4];
	    double *right_k3_p = &right[j*16 + 3*4];
	    __m128d x2_0 = _mm_load_pd( &x2_p[0] );
	    __m128d x2_2 = _mm_load_pd( &x2_p[2] );
	    
	    __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
	    __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
	    __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
	    __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
	    __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
	    __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
	    __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
	    __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );
	    
	    right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	    right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	    
	    right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	    right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	    
	    right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	    right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	    right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	    
	    right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	    right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	    	    
	    right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	    right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	    
	    right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	    right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	    right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   		 		
	    
	    __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
	    __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );		 		 	   
	    
	    __m128d EV_t_l0_k0 = EVV[0];
	    __m128d EV_t_l0_k2 = EVV[1];
	    __m128d EV_t_l1_k0 = EVV[2];
	    __m128d EV_t_l1_k2 = EVV[3];
	    __m128d EV_t_l2_k0 = EVV[4];
	    __m128d EV_t_l2_k2 = EVV[5];
	    __m128d EV_t_l3_k0 = EVV[6]; 
	    __m128d EV_t_l3_k2 = EVV[7];
	    
	    EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	    EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	    
	    EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	    EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	    
	    EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	    
	    EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	    EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	    
	    EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	    EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	    EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	    
	    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
	    
	    
	    values[j * 2] = EV_t_l0_k0;
	    values[j * 2 + 1] = EV_t_l2_k0;            	   	    
	    
	    maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
	    maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));
	  }
		     
	_mm_store_pd(maxima, maxv);
	
	max = MAX(maxima[0], maxima[1]);
	
	if(max < minlikelihood)
	  {
	    __m128d sv = _mm_set1_pd(twotothe256);
	    
	    scaleGap = 1;
	    
	    _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
	    _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
	    _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
	    _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
	    _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
	    _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
	    _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
	    _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     	    	 
	  }
	else
	  {
	    _mm_store_pd(&x3[0], values[0]);	   
	    _mm_store_pd(&x3[2], values[1]);
	    _mm_store_pd(&x3[4], values[2]);
	    _mm_store_pd(&x3[6], values[3]);
	    _mm_store_pd(&x3[8], values[4]);	   
	    _mm_store_pd(&x3[10], values[5]);
	    _mm_store_pd(&x3[12], values[6]);
	    _mm_store_pd(&x3[14], values[7]);
	  }
      }

     
      x3 = x3_start;

     for (i = 0; i < n; i++)
       { 
	 if(x3_gap[i / 32] & mask32[i % 32])
	   {	     
	     if(scaleGap)
	       {		 
		 addScale += wgt[i];		 	       
	       }
	   }
	 else
	   {
	     __m128d maxv =_mm_setzero_pd();	     	    
	     
	     if(x1_gap[i / 32] & mask32[i % 32])
	       x1 = x1_gapColumn;
	     else
	       {
		 x1 = x1_ptr;
		 x1_ptr += 16;
	       }
	     
	     if(x2_gap[i / 32] & mask32[i % 32])
	       x2 = x2_gapColumn;
	     else
	       {
		 x2 = x2_ptr;
		 x2_ptr += 16;
	       }
	     
	     
	     for (j = 0; j < 4; j++)
	       {
		 
		 double *x1_p = &x1[j*4];
		 double *left_k0_p = &left[j*16];
		 double *left_k1_p = &left[j*16 + 1*4];
		 double *left_k2_p = &left[j*16 + 2*4];
		 double *left_k3_p = &left[j*16 + 3*4];
		 
		 __m128d x1_0 = _mm_load_pd( &x1_p[0] );
		 __m128d x1_2 = _mm_load_pd( &x1_p[2] );
		 
		 __m128d left_k0_0 = _mm_load_pd( &left_k0_p[0] );
		 __m128d left_k0_2 = _mm_load_pd( &left_k0_p[2] );
		 __m128d left_k1_0 = _mm_load_pd( &left_k1_p[0] );
		 __m128d left_k1_2 = _mm_load_pd( &left_k1_p[2] );
		 __m128d left_k2_0 = _mm_load_pd( &left_k2_p[0] );
		 __m128d left_k2_2 = _mm_load_pd( &left_k2_p[2] );
		 __m128d left_k3_0 = _mm_load_pd( &left_k3_p[0] );
		 __m128d left_k3_2 = _mm_load_pd( &left_k3_p[2] );
		 
		 left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
		 left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
		 
		 left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
		 left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
		 
		 left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
		 left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
		 left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
		 
		 left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
		 left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
		 
		 left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
		 left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
		 
		 left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
		 left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
		 left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
		 
		 
		 //
		 // multiply/add right side
		 //
		 double *x2_p = &x2[j*4];
		 double *right_k0_p = &right[j*16];
		 double *right_k1_p = &right[j*16 + 1*4];
		 double *right_k2_p = &right[j*16 + 2*4];
		 double *right_k3_p = &right[j*16 + 3*4];
		 __m128d x2_0 = _mm_load_pd( &x2_p[0] );
		 __m128d x2_2 = _mm_load_pd( &x2_p[2] );
		 
		 __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
		 __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
		 __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
		 __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
		 __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
		 __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
		 __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
		 __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );
		 
		 right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
		 right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
		 
		 right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
		 right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
		 
		 right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
		 right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
		 right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
		 
		 right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
		 right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
		 
		 
		 right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
		 right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
		 
		 right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
		 right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
		 right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   
		 
		 //
		 // multiply left * right
		 //
		 
		 __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
		 __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );
		 
		 
		 //
		 // multiply with EV matrix (!?)
		 //	     
		 
		 __m128d EV_t_l0_k0 = EVV[0];
		 __m128d EV_t_l0_k2 = EVV[1];
		 __m128d EV_t_l1_k0 = EVV[2];
		 __m128d EV_t_l1_k2 = EVV[3];
		 __m128d EV_t_l2_k0 = EVV[4];
		 __m128d EV_t_l2_k2 = EVV[5];
		 __m128d EV_t_l3_k0 = EVV[6]; 
		 __m128d EV_t_l3_k2 = EVV[7];
		 
		 
		 EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
		 EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
		 EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
		 
		 EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
		 EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
		 
		 EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
		 EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
		 
		 EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
		 EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
		 EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
		 
		 
		 EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
		 EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
		 EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
		 
		 EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
		 
		 
		 values[j * 2] = EV_t_l0_k0;
		 values[j * 2 + 1] = EV_t_l2_k0;            	   	    
		 
		 maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
		 maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));
	       }
	     
	     
	     _mm_store_pd(maxima, maxv);
	     
	     max = MAX(maxima[0], maxima[1]);
	     
	     if(max < minlikelihood)
	       {
		 __m128d sv = _mm_set1_pd(twotothe256);
		 
		 _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
		 _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
		 _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
		 _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
		 _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
		 _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
		 _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
		 _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     
		 
		 
		 addScale += wgt[i];
		
	       }
	     else
	       {
		 _mm_store_pd(&x3[0], values[0]);	   
		 _mm_store_pd(&x3[2], values[1]);
		 _mm_store_pd(&x3[4], values[2]);
		 _mm_store_pd(&x3[6], values[3]);
		 _mm_store_pd(&x3[8], values[4]);	   
		 _mm_store_pd(&x3[10], values[5]);
		 _mm_store_pd(&x3[12], values[6]);
		 _mm_store_pd(&x3[14], values[7]);
	       }	 

	    
		 
	     x3 += 16;

	   }
       }
     break;
    default:
      assert(0);
    }
  
 
  *scalerIncrement = addScale;
}


static void newviewGTRGAMMA(int tipCase,
			    double *x1_start, double *x2_start, double *x3_start,
			    double *EV, double *tipVector,
			    unsigned char *tipX1, unsigned char *tipX2,
			    const size_t n, double *left, double *right, int *wgt, int *scalerIncrement
			    )
{
  size_t
    i, 
    j, 
    k, 
    l;

  int
    addScale = 0;
  
  double
    *x1,
    *x2,
    *x3,
    max,
    maxima[2] __attribute__ ((aligned (BYTE_ALIGNMENT))),       
    EV_t[16] __attribute__ ((aligned (BYTE_ALIGNMENT)));      
    
  __m128d 
    values[8],
    EVV[8];  

  for(k = 0; k < 4; k++)
    for (l=0; l < 4; l++)
      EV_t[4 * l + k] = EV[4 * k + l];

  for(k = 0; k < 8; k++)
    EVV[k] = _mm_load_pd(&EV_t[k * 2]);
   
  switch(tipCase)
    {
    case TIP_TIP:
      {
	double *uX1, umpX1[256] __attribute__ ((aligned (BYTE_ALIGNMENT))), *uX2, umpX2[256] __attribute__ ((aligned (BYTE_ALIGNMENT)));


	for (i = 1; i < 16; i++)
	  {
	    __m128d x1_1 = _mm_load_pd(&(tipVector[i*4]));
	    __m128d x1_2 = _mm_load_pd(&(tipVector[i*4 + 2]));	   

	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m128d left1 = _mm_load_pd(&left[j*16 + k*4]);
		  __m128d left2 = _mm_load_pd(&left[j*16 + k*4 + 2]);
		  
		  __m128d acc = _mm_setzero_pd();

		  acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
		  acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));
		  		  
		  acc = _mm_hadd_pd(acc, acc);
		  _mm_storel_pd(&umpX1[i*16 + j*4 + k], acc);
		}
	  
	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{
		  __m128d left1 = _mm_load_pd(&right[j*16 + k*4]);
		  __m128d left2 = _mm_load_pd(&right[j*16 + k*4 + 2]);
		  
		  __m128d acc = _mm_setzero_pd();

		  acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
		  acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));
		  		  
		  acc = _mm_hadd_pd(acc, acc);
		  _mm_storel_pd(&umpX2[i*16 + j*4 + k], acc);
		 
		}
	  }   	
	  
	for (i = 0; i < n; i++)
	  {
	    x3 = &x3_start[i * 16];

	    
	    uX1 = &umpX1[16 * tipX1[i]];
	    uX2 = &umpX2[16 * tipX2[i]];	   	    	    
	    
	    for (j = 0; j < 4; j++)
	       {				 		  		  		   
		 __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
		 __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
		 				  
		   
		 __m128d uX2_k0_sse = _mm_load_pd( &uX2[j * 4] );
		 __m128d uX2_k2_sse = _mm_load_pd( &uX2[j * 4 + 2] );
 		 

		 //
		 // multiply left * right
		 //
		 
		 __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, uX2_k0_sse );
		 __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, uX2_k2_sse );
		 
		 
		 //
		 // multiply with EV matrix (!?)
		 //
		 
		 __m128d EV_t_l0_k0 = EVV[0];
		 __m128d EV_t_l0_k2 = EVV[1];
		 __m128d EV_t_l1_k0 = EVV[2];
		 __m128d EV_t_l1_k2 = EVV[3];
		 __m128d EV_t_l2_k0 = EVV[4];
		 __m128d EV_t_l2_k2 = EVV[5];
		 __m128d EV_t_l3_k0 = EVV[6]; 
		 __m128d EV_t_l3_k2 = EVV[7];
		 
		 EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
		 EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
		 EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
		 
		 EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
		 EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
		 
		 EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
		 EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
		 
		 EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
		 EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
		 EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
		 
		 EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
		 EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
		 EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
		 
		 EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
		 
		 _mm_store_pd( &x3[j * 4 + 0], EV_t_l0_k0 );
		 _mm_store_pd( &x3[j * 4 + 2], EV_t_l2_k0 );
	       }
	  }
      }
      break;
    case TIP_INNER:
      {	
	double *uX1, umpX1[256] __attribute__ ((aligned (BYTE_ALIGNMENT)));


	for (i = 1; i < 16; i++)
	  {
	    __m128d x1_1 = _mm_load_pd(&(tipVector[i*4]));
	    __m128d x1_2 = _mm_load_pd(&(tipVector[i*4 + 2]));	   

	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m128d left1 = _mm_load_pd(&left[j*16 + k*4]);
		  __m128d left2 = _mm_load_pd(&left[j*16 + k*4 + 2]);
		  
		  __m128d acc = _mm_setzero_pd();

		  acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
		  acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));
		  		  
		  acc = _mm_hadd_pd(acc, acc);
		  _mm_storel_pd(&umpX1[i*16 + j*4 + k], acc);		 
		}
	  }

	 for (i = 0; i < n; i++)
	   {
	     __m128d maxv =_mm_setzero_pd();
	     
	     x2 = &x2_start[i * 16];
	     x3 = &x3_start[i * 16];

	     uX1 = &umpX1[16 * tipX1[i]];	     

	     for (j = 0; j < 4; j++)
	       {

		 //
		 // multiply/add right side
		 //
		 double *x2_p = &x2[j*4];
		 double *right_k0_p = &right[j*16];
		 double *right_k1_p = &right[j*16 + 1*4];
		 double *right_k2_p = &right[j*16 + 2*4];
		 double *right_k3_p = &right[j*16 + 3*4];
		 __m128d x2_0 = _mm_load_pd( &x2_p[0] );
		 __m128d x2_2 = _mm_load_pd( &x2_p[2] );

		 __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
		 __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
		 __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
		 __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
		 __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
		 __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
		 __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
		 __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );



		 right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
		 right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

		 right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
		 right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

		 right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
		 right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
		 right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);


		 right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
		 right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);


		 right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
		 right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

		 right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
		 right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
		 right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);

		 {
		   //
		   // load left side from tip vector
		   //
		   
		   __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
		   __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
		 
		 
		   //
		   // multiply left * right
		   //
		   
		   __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, right_k0_0 );
		   __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, right_k2_0 );
		   
		   
		   //
		   // multiply with EV matrix (!?)
		   //		   		  

		   __m128d EV_t_l0_k0 = EVV[0];
		   __m128d EV_t_l0_k2 = EVV[1];
		   __m128d EV_t_l1_k0 = EVV[2];
		   __m128d EV_t_l1_k2 = EVV[3];
		   __m128d EV_t_l2_k0 = EVV[4];
		   __m128d EV_t_l2_k2 = EVV[5];
		   __m128d EV_t_l3_k0 = EVV[6]; 
		   __m128d EV_t_l3_k2 = EVV[7];

		   
		   EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
		   EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
		   EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
		   
		   EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
		   EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
		   
		   EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
		   EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
		   
		   EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
		   EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
		   EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
		   		   
		   EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
		   EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
		   EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
		   
		   EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
		   
		   values[j * 2]     = EV_t_l0_k0;
		   values[j * 2 + 1] = EV_t_l2_k0;		   		   
		   
		   maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
		   maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));		   
		 }
	       }

	     
	     _mm_store_pd(maxima, maxv);

	     max = MAX(maxima[0], maxima[1]);

	     if(max < minlikelihood)
	       {
		 __m128d sv = _mm_set1_pd(twotothe256);
	       		       	   	 	     
		 _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
		 _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
		 _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
		 _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
		 _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
		 _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
		 _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
		 _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     
		 
		 
		 addScale += wgt[i];
		 
	       }
	     else
	       {
		 _mm_store_pd(&x3[0], values[0]);	   
		 _mm_store_pd(&x3[2], values[1]);
		 _mm_store_pd(&x3[4], values[2]);
		 _mm_store_pd(&x3[6], values[3]);
		 _mm_store_pd(&x3[8], values[4]);	   
		 _mm_store_pd(&x3[10], values[5]);
		 _mm_store_pd(&x3[12], values[6]);
		 _mm_store_pd(&x3[14], values[7]);
	       }
	   }
      }
      break;
    case INNER_INNER:     
     for (i = 0; i < n; i++)
       {
	 __m128d maxv =_mm_setzero_pd();
	 

	 x1 = &x1_start[i * 16];
	 x2 = &x2_start[i * 16];
	 x3 = &x3_start[i * 16];
	 
	 for (j = 0; j < 4; j++)
	   {
	     
	     double *x1_p = &x1[j*4];
	     double *left_k0_p = &left[j*16];
	     double *left_k1_p = &left[j*16 + 1*4];
	     double *left_k2_p = &left[j*16 + 2*4];
	     double *left_k3_p = &left[j*16 + 3*4];
	     
	     __m128d x1_0 = _mm_load_pd( &x1_p[0] );
	     __m128d x1_2 = _mm_load_pd( &x1_p[2] );
	     
	     __m128d left_k0_0 = _mm_load_pd( &left_k0_p[0] );
	     __m128d left_k0_2 = _mm_load_pd( &left_k0_p[2] );
	     __m128d left_k1_0 = _mm_load_pd( &left_k1_p[0] );
	     __m128d left_k1_2 = _mm_load_pd( &left_k1_p[2] );
	     __m128d left_k2_0 = _mm_load_pd( &left_k2_p[0] );
	     __m128d left_k2_2 = _mm_load_pd( &left_k2_p[2] );
	     __m128d left_k3_0 = _mm_load_pd( &left_k3_p[0] );
	     __m128d left_k3_2 = _mm_load_pd( &left_k3_p[2] );
	     
	     left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
	     left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
	     
	     left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
	     left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
	     
	     left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
	     left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
	     left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
	     
	     left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
	     left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
	     
	     left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
	     left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
	     
	     left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
	     left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
	     left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
	     
	     
	     //
	     // multiply/add right side
	     //
	     double *x2_p = &x2[j*4];
	     double *right_k0_p = &right[j*16];
	     double *right_k1_p = &right[j*16 + 1*4];
	     double *right_k2_p = &right[j*16 + 2*4];
	     double *right_k3_p = &right[j*16 + 3*4];
	     __m128d x2_0 = _mm_load_pd( &x2_p[0] );
	     __m128d x2_2 = _mm_load_pd( &x2_p[2] );
	     
	     __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
	     __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
	     __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
	     __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
	     __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
	     __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
	     __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
	     __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );
	     
	     right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	     right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	     
	     right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	     right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	     
	     right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	     right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	     right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	     
	     right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	     right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	     
	     
	     right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	     right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	     
	     right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	     right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	     right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   

             //
             // multiply left * right
             //

	     __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
	     __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );


             //
             // multiply with EV matrix (!?)
             //	     

	     __m128d EV_t_l0_k0 = EVV[0];
	     __m128d EV_t_l0_k2 = EVV[1];
	     __m128d EV_t_l1_k0 = EVV[2];
	     __m128d EV_t_l1_k2 = EVV[3];
	     __m128d EV_t_l2_k0 = EVV[4];
	     __m128d EV_t_l2_k2 = EVV[5];
	     __m128d EV_t_l3_k0 = EVV[6]; 
	     __m128d EV_t_l3_k2 = EVV[7];


	    EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	    EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

	    EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	    EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

	    EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

	    EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	    EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );


	    EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
            EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
            EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

            EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );

	    
	    values[j * 2] = EV_t_l0_k0;
	    values[j * 2 + 1] = EV_t_l2_k0;            	   	    

	    maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
	    maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));
           }
	 	 
	 
	 _mm_store_pd(maxima, maxv);
	 
	 max = MAX(maxima[0], maxima[1]);
	 
	 if(max < minlikelihood)
	   {
	     __m128d sv = _mm_set1_pd(twotothe256);
	       		       	   	 	     
	     _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
	     _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
	     _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
	     _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
	     _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
	     _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
	     _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
	     _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     
	     
	    
	     addScale += wgt[i];
	    
	   }
	 else
	   {
	     _mm_store_pd(&x3[0], values[0]);	   
	     _mm_store_pd(&x3[2], values[1]);
	     _mm_store_pd(&x3[4], values[2]);
	     _mm_store_pd(&x3[6], values[3]);
	     _mm_store_pd(&x3[8], values[4]);	   
	     _mm_store_pd(&x3[10], values[5]);
	     _mm_store_pd(&x3[12], values[6]);
	     _mm_store_pd(&x3[14], values[7]);
	   }	 
       }
   
     break;
    default:
      assert(0);
    }
  
 
  *scalerIncrement = addScale;

}
static void newviewGTRCAT( int tipCase,  double *EV,  int *cptr,
			   double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
			   unsigned char *tipX1, unsigned char *tipX2,
			   size_t n,  double *left, double *right, int *wgt, int *scalerIncrement)
{
  double
    *le,
    *ri,
    *x1,
    *x2, 
    *x3, 
    EV_t[16] __attribute__ ((aligned (BYTE_ALIGNMENT)));
    
  size_t
    i, 
    j;
  
  int
    scale,
    addScale = 0;
   
  __m128d
    minlikelihood_sse = _mm_set1_pd( minlikelihood ),
    sc = _mm_set1_pd(twotothe256),
    EVV[8];  
  
  for(i = 0; i < 4; i++)
    for (j=0; j < 4; j++)
      EV_t[4 * j + i] = EV[4 * i + j];
  
  for(i = 0; i < 8; i++)
    EVV[i] = _mm_load_pd(&EV_t[i * 2]);
  
  switch(tipCase)
    {
    case TIP_TIP:      
      for (i = 0; i < n; i++)
	{	 
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &(tipVector[4 * tipX2[i]]);
	  
	  x3 = &x3_start[i * 4];
	  
	  le =  &left[cptr[i] * 16];
	  ri =  &right[cptr[i] * 16];
	  
	  __m128d x1_0 = _mm_load_pd( &x1[0] );
	  __m128d x1_2 = _mm_load_pd( &x1[2] );
	  
	  __m128d left_k0_0 = _mm_load_pd( &le[0] );
	  __m128d left_k0_2 = _mm_load_pd( &le[2] );
	  __m128d left_k1_0 = _mm_load_pd( &le[4] );
	  __m128d left_k1_2 = _mm_load_pd( &le[6] );
	  __m128d left_k2_0 = _mm_load_pd( &le[8] );
	  __m128d left_k2_2 = _mm_load_pd( &le[10] );
	  __m128d left_k3_0 = _mm_load_pd( &le[12] );
	  __m128d left_k3_2 = _mm_load_pd( &le[14] );
	  
	  left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
	  left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
	  
	  left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
	  left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
	  
	  left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
	  left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
	  left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
	  
	  left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
	  left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
	  
	  left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
	  left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
	  
	  left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
	  left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
	  left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
	  
	  __m128d x2_0 = _mm_load_pd( &x2[0] );
	  __m128d x2_2 = _mm_load_pd( &x2[2] );
	  
	  __m128d right_k0_0 = _mm_load_pd( &ri[0] );
	  __m128d right_k0_2 = _mm_load_pd( &ri[2] );
	  __m128d right_k1_0 = _mm_load_pd( &ri[4] );
	  __m128d right_k1_2 = _mm_load_pd( &ri[6] );
	  __m128d right_k2_0 = _mm_load_pd( &ri[8] );
	  __m128d right_k2_2 = _mm_load_pd( &ri[10] );
	  __m128d right_k3_0 = _mm_load_pd( &ri[12] );
	  __m128d right_k3_2 = _mm_load_pd( &ri[14] );
	  
	  right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	  right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	  
	  right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	  right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	  
	  right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	  right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	  right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	  
	  right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	  right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	  
	  right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	  right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	  
	  right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	  right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	  right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   
	  
	  __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
	  __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );	  	  

	  __m128d EV_t_l0_k0 = EVV[0];
	  __m128d EV_t_l0_k2 = EVV[1];
	  __m128d EV_t_l1_k0 = EVV[2];
	  __m128d EV_t_l1_k2 = EVV[3];
	  __m128d EV_t_l2_k0 = EVV[4];
	  __m128d EV_t_l2_k2 = EVV[5];
	  __m128d EV_t_l3_k0 = EVV[6];
	  __m128d EV_t_l3_k2 = EVV[7];
	  
	  EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	  EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	  EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	  
	  EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	  EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	  
	  EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	  EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	  
	  EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	  EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	  EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	  	  
	  EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	  EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	  EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	  
	  EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );	 
	  	  
	  _mm_store_pd(x3, EV_t_l0_k0);
	  _mm_store_pd(&x3[2], EV_t_l2_k0);	  	 	   	    
	}
      break;
    case TIP_INNER:      
      for (i = 0; i < n; i++)
	{
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &x2_start[4 * i];
	  x3 = &x3_start[4 * i];
	  
	  le =  &left[cptr[i] * 16];
	  ri =  &right[cptr[i] * 16];

	  __m128d x1_0 = _mm_load_pd( &x1[0] );
	  __m128d x1_2 = _mm_load_pd( &x1[2] );
	  
	  __m128d left_k0_0 = _mm_load_pd( &le[0] );
	  __m128d left_k0_2 = _mm_load_pd( &le[2] );
	  __m128d left_k1_0 = _mm_load_pd( &le[4] );
	  __m128d left_k1_2 = _mm_load_pd( &le[6] );
	  __m128d left_k2_0 = _mm_load_pd( &le[8] );
	  __m128d left_k2_2 = _mm_load_pd( &le[10] );
	  __m128d left_k3_0 = _mm_load_pd( &le[12] );
	  __m128d left_k3_2 = _mm_load_pd( &le[14] );
	  
	  left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
	  left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
	  
	  left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
	  left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
	  
	  left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
	  left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
	  left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
	  
	  left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
	  left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
	  
	  left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
	  left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
	  
	  left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
	  left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
	  left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
	  
	  __m128d x2_0 = _mm_load_pd( &x2[0] );
	  __m128d x2_2 = _mm_load_pd( &x2[2] );
	  
	  __m128d right_k0_0 = _mm_load_pd( &ri[0] );
	  __m128d right_k0_2 = _mm_load_pd( &ri[2] );
	  __m128d right_k1_0 = _mm_load_pd( &ri[4] );
	  __m128d right_k1_2 = _mm_load_pd( &ri[6] );
	  __m128d right_k2_0 = _mm_load_pd( &ri[8] );
	  __m128d right_k2_2 = _mm_load_pd( &ri[10] );
	  __m128d right_k3_0 = _mm_load_pd( &ri[12] );
	  __m128d right_k3_2 = _mm_load_pd( &ri[14] );
	  
	  right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	  right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	  
	  right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	  right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	  
	  right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	  right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	  right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	  
	  right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	  right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	  
	  right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	  right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	  
	  right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	  right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	  right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   
	  
	  __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
	  __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );
	  
	  __m128d EV_t_l0_k0 = EVV[0];
	  __m128d EV_t_l0_k2 = EVV[1];
	  __m128d EV_t_l1_k0 = EVV[2];
	  __m128d EV_t_l1_k2 = EVV[3];
	  __m128d EV_t_l2_k0 = EVV[4];
	  __m128d EV_t_l2_k2 = EVV[5];
	  __m128d EV_t_l3_k0 = EVV[6];
	  __m128d EV_t_l3_k2 = EVV[7];
	 
	  
	  EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	  EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	  EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	  
	  EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	  EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	  
	  EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	  EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	  
	  EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	  EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	  EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	  	  
	  EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	  EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	  EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	  
	  EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );	  	 	    		  
	 
	  scale = 1;
	  	  	  	    
	  __m128d v1 = _mm_and_pd(EV_t_l0_k0, absMask.m);
	  v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	  if(_mm_movemask_pd( v1 ) != 3)
	    scale = 0;
	  else
	    {
	      v1 = _mm_and_pd(EV_t_l2_k0, absMask.m);
	      v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	      if(_mm_movemask_pd( v1 ) != 3)
		scale = 0;
	    }
	  	  
	  if(scale)
	    {		      
	      _mm_store_pd(&x3[0], _mm_mul_pd(EV_t_l0_k0, sc));
	      _mm_store_pd(&x3[2], _mm_mul_pd(EV_t_l2_k0, sc));	      	      
	      
	      
	      addScale += wgt[i];	  
	    }	
	  else
	    {
	      _mm_store_pd(x3, EV_t_l0_k0);
	      _mm_store_pd(&x3[2], EV_t_l2_k0);
	    }
	 
	  	  
	}
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];
	  x3 = &x3_start[4 * i];
	  
	  le =  &left[cptr[i] * 16];
	  ri =  &right[cptr[i] * 16];

	  __m128d x1_0 = _mm_load_pd( &x1[0] );
	  __m128d x1_2 = _mm_load_pd( &x1[2] );
	  
	  __m128d left_k0_0 = _mm_load_pd( &le[0] );
	  __m128d left_k0_2 = _mm_load_pd( &le[2] );
	  __m128d left_k1_0 = _mm_load_pd( &le[4] );
	  __m128d left_k1_2 = _mm_load_pd( &le[6] );
	  __m128d left_k2_0 = _mm_load_pd( &le[8] );
	  __m128d left_k2_2 = _mm_load_pd( &le[10] );
	  __m128d left_k3_0 = _mm_load_pd( &le[12] );
	  __m128d left_k3_2 = _mm_load_pd( &le[14] );
	  
	  left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
	  left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
	  
	  left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
	  left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
	  
	  left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
	  left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
	  left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
	  
	  left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
	  left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
	  
	  left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
	  left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
	  
	  left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
	  left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
	  left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
	  
	  __m128d x2_0 = _mm_load_pd( &x2[0] );
	  __m128d x2_2 = _mm_load_pd( &x2[2] );
	  
	  __m128d right_k0_0 = _mm_load_pd( &ri[0] );
	  __m128d right_k0_2 = _mm_load_pd( &ri[2] );
	  __m128d right_k1_0 = _mm_load_pd( &ri[4] );
	  __m128d right_k1_2 = _mm_load_pd( &ri[6] );
	  __m128d right_k2_0 = _mm_load_pd( &ri[8] );
	  __m128d right_k2_2 = _mm_load_pd( &ri[10] );
	  __m128d right_k3_0 = _mm_load_pd( &ri[12] );
	  __m128d right_k3_2 = _mm_load_pd( &ri[14] );
	  
	  right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	  right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	  
	  right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	  right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	  
	  right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	  right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	  right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	  
	  right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	  right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	  
	  right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	  right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	  
	  right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	  right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	  right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   
	  
	  __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
	  __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );
	  
	  __m128d EV_t_l0_k0 = EVV[0];
	  __m128d EV_t_l0_k2 = EVV[1];
	  __m128d EV_t_l1_k0 = EVV[2];
	  __m128d EV_t_l1_k2 = EVV[3];
	  __m128d EV_t_l2_k0 = EVV[4];
	  __m128d EV_t_l2_k2 = EVV[5];
	  __m128d EV_t_l3_k0 = EVV[6];
	  __m128d EV_t_l3_k2 = EVV[7];
	 
	  
	  EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	  EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	  EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	  
	  EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	  EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	  
	  EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	  EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	  
	  EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	  EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	  EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	  	  
	  EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	  EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	  EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	  
	  EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );	  	 	    		  	 

	  scale = 1;
	  	  
	  __m128d v1 = _mm_and_pd(EV_t_l0_k0, absMask.m);
	  v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	  if(_mm_movemask_pd( v1 ) != 3)
	    scale = 0;
	  else
	    {
	      v1 = _mm_and_pd(EV_t_l2_k0, absMask.m);
	      v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	      if(_mm_movemask_pd( v1 ) != 3)
		scale = 0;
	    }
	  	  
	  if(scale)
	    {		      
	      _mm_store_pd(&x3[0], _mm_mul_pd(EV_t_l0_k0, sc));
	      _mm_store_pd(&x3[2], _mm_mul_pd(EV_t_l2_k0, sc));	      	      
	      
	      
	      addScale += wgt[i];	  
	    }	
	  else
	    {
	      _mm_store_pd(x3, EV_t_l0_k0);
	      _mm_store_pd(&x3[2], EV_t_l2_k0);
	    }
	  	  
	}
      break;
    default:
      assert(0);
    }

  
  *scalerIncrement = addScale;
}



static void newviewGTRCAT_SAVE( int tipCase,  double *EV,  int *cptr,
				double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
				unsigned char *tipX1, unsigned char *tipX2,
				size_t n,  double *left, double *right, int *wgt, int *scalerIncrement,
				unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
				double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats)
{
  double
    *le,
    *ri,
    *x1,
    *x2,
    *x3,
    *x1_ptr = x1_start,
    *x2_ptr = x2_start, 
    *x3_ptr = x3_start, 
    EV_t[16] __attribute__ ((aligned (BYTE_ALIGNMENT)));
    
  size_t
    i, 
    j;

  int
    scale, 
    scaleGap = 0,
    addScale = 0;
   
  __m128d
    minlikelihood_sse = _mm_set1_pd( minlikelihood ),
    sc = _mm_set1_pd(twotothe256),
    EVV[8];  
  
  for(i = 0; i < 4; i++)
    for (j=0; j < 4; j++)
      EV_t[4 * j + i] = EV[4 * i + j];
  
  for(i = 0; i < 8; i++)
    EVV[i] = _mm_load_pd(&EV_t[i * 2]);
  
  {
    x1 = x1_gapColumn;	      
    x2 = x2_gapColumn;
    x3 = x3_gapColumn;
    
    le =  &left[maxCats * 16];	     	 
    ri =  &right[maxCats * 16];		   	  	  	  	         

    __m128d x1_0 = _mm_load_pd( &x1[0] );
    __m128d x1_2 = _mm_load_pd( &x1[2] );
    
    __m128d left_k0_0 = _mm_load_pd( &le[0] );
    __m128d left_k0_2 = _mm_load_pd( &le[2] );
    __m128d left_k1_0 = _mm_load_pd( &le[4] );
    __m128d left_k1_2 = _mm_load_pd( &le[6] );
    __m128d left_k2_0 = _mm_load_pd( &le[8] );
    __m128d left_k2_2 = _mm_load_pd( &le[10] );
    __m128d left_k3_0 = _mm_load_pd( &le[12] );
    __m128d left_k3_2 = _mm_load_pd( &le[14] );
    
    left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
    left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
    
    left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
    left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
    
    left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
    left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
    left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
    
    left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
    left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
    
    left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
    left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
    
    left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
    left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
    left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
    
    __m128d x2_0 = _mm_load_pd( &x2[0] );
    __m128d x2_2 = _mm_load_pd( &x2[2] );
    
    __m128d right_k0_0 = _mm_load_pd( &ri[0] );
    __m128d right_k0_2 = _mm_load_pd( &ri[2] );
    __m128d right_k1_0 = _mm_load_pd( &ri[4] );
    __m128d right_k1_2 = _mm_load_pd( &ri[6] );
    __m128d right_k2_0 = _mm_load_pd( &ri[8] );
    __m128d right_k2_2 = _mm_load_pd( &ri[10] );
    __m128d right_k3_0 = _mm_load_pd( &ri[12] );
    __m128d right_k3_2 = _mm_load_pd( &ri[14] );
    
    right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
    right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
    
    right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
    right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
    
    right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
    right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
    right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
    
    right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
    right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
    
    right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
    right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
    
    right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
    right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
    right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   
    
    __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
    __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );
    
    __m128d EV_t_l0_k0 = EVV[0];
    __m128d EV_t_l0_k2 = EVV[1];
    __m128d EV_t_l1_k0 = EVV[2];
    __m128d EV_t_l1_k2 = EVV[3];
    __m128d EV_t_l2_k0 = EVV[4];
    __m128d EV_t_l2_k2 = EVV[5];
    __m128d EV_t_l3_k0 = EVV[6];
    __m128d EV_t_l3_k2 = EVV[7];
        
    EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
    EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
    
    EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
    EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
    
    EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
    
    EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
    EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
    
    EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
    EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
    EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
    
    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );	  	 	    		  
	
    if(tipCase != TIP_TIP)
      {    
	scale = 1;
	      
	__m128d v1 = _mm_and_pd(EV_t_l0_k0, absMask.m);
	v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	if(_mm_movemask_pd( v1 ) != 3)
	  scale = 0;
	else
	  {
	    v1 = _mm_and_pd(EV_t_l2_k0, absMask.m);
	    v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	    if(_mm_movemask_pd( v1 ) != 3)
	      scale = 0;
	  }
	
	if(scale)
	  {		      
	    _mm_store_pd(&x3[0], _mm_mul_pd(EV_t_l0_k0, sc));
	    _mm_store_pd(&x3[2], _mm_mul_pd(EV_t_l2_k0, sc));	      	      
	    
	    scaleGap = TRUE;	   
	  }	
	else
	  {
	    _mm_store_pd(x3, EV_t_l0_k0);
	    _mm_store_pd(&x3[2], EV_t_l2_k0);
	  }
      }
    else
      {
	_mm_store_pd(x3, EV_t_l0_k0);
	_mm_store_pd(&x3[2], EV_t_l2_k0);
      }
  }
  

  switch(tipCase)
    {
    case TIP_TIP:      
      for (i = 0; i < n; i++)
	{
	  if(noGap(x3_gap, i))
	    {
	      x1 = &(tipVector[4 * tipX1[i]]);
	      x2 = &(tipVector[4 * tipX2[i]]);
	  
	      x3 = x3_ptr;
	  
	      if(isGap(x1_gap, i))
		le =  &left[maxCats * 16];
	      else	  	  
		le =  &left[cptr[i] * 16];	  
	  
	      if(isGap(x2_gap, i))
		ri =  &right[maxCats * 16];
	      else	 	  
		ri =  &right[cptr[i] * 16];
	  
	      __m128d x1_0 = _mm_load_pd( &x1[0] );
	      __m128d x1_2 = _mm_load_pd( &x1[2] );
	      
	      __m128d left_k0_0 = _mm_load_pd( &le[0] );
	      __m128d left_k0_2 = _mm_load_pd( &le[2] );
	      __m128d left_k1_0 = _mm_load_pd( &le[4] );
	      __m128d left_k1_2 = _mm_load_pd( &le[6] );
	      __m128d left_k2_0 = _mm_load_pd( &le[8] );
	      __m128d left_k2_2 = _mm_load_pd( &le[10] );
	      __m128d left_k3_0 = _mm_load_pd( &le[12] );
	      __m128d left_k3_2 = _mm_load_pd( &le[14] );
	  
	      left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
	      left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
	      
	      left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
	      left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
	      
	      left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
	      left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
	      left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
	      
	      left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
	      left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
	      
	      left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
	      left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
	      
	      left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
	      left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
	      left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
	      
	      __m128d x2_0 = _mm_load_pd( &x2[0] );
	      __m128d x2_2 = _mm_load_pd( &x2[2] );
	      
	      __m128d right_k0_0 = _mm_load_pd( &ri[0] );
	      __m128d right_k0_2 = _mm_load_pd( &ri[2] );
	      __m128d right_k1_0 = _mm_load_pd( &ri[4] );
	      __m128d right_k1_2 = _mm_load_pd( &ri[6] );
	      __m128d right_k2_0 = _mm_load_pd( &ri[8] );
	      __m128d right_k2_2 = _mm_load_pd( &ri[10] );
	      __m128d right_k3_0 = _mm_load_pd( &ri[12] );
	      __m128d right_k3_2 = _mm_load_pd( &ri[14] );
	      
	      right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	      right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	      
	      right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	      right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	      
	      right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	      right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	      right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	      
	      right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	      right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	      
	      right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	      right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	      
	      right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	      right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	      right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   
	      
	      __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
	      __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );	  	  
	      
	      __m128d EV_t_l0_k0 = EVV[0];
	      __m128d EV_t_l0_k2 = EVV[1];
	      __m128d EV_t_l1_k0 = EVV[2];
	      __m128d EV_t_l1_k2 = EVV[3];
	      __m128d EV_t_l2_k0 = EVV[4];
	      __m128d EV_t_l2_k2 = EVV[5];
	      __m128d EV_t_l3_k0 = EVV[6];
	      __m128d EV_t_l3_k2 = EVV[7];
	      
	      EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	      EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	      EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	      
	      EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	      EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	      
	      EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	      EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	      
	      EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	      EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	      EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	      
	      EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	      EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	      EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	      
	      EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );	 
	  	  
	      _mm_store_pd(x3, EV_t_l0_k0);
	      _mm_store_pd(&x3[2], EV_t_l2_k0);	  	 	   	    

	      x3_ptr += 4;
	    }
	}
      break;
    case TIP_INNER:      
      for (i = 0; i < n; i++)
	{ 
	  if(isGap(x3_gap, i))
	    {
	      if(scaleGap)		   		    
		addScale += wgt[i];
	    }
	  else
	    {	      
	      x1 = &(tipVector[4 * tipX1[i]]);
	      	   
	      x3 = x3_ptr;

	      if(isGap(x1_gap, i))
		le =  &left[maxCats * 16];
	      else
		le =  &left[cptr[i] * 16];

	      if(isGap(x2_gap, i))
		{		 
		  ri =  &right[maxCats * 16];
		  x2 = x2_gapColumn;
		}
	      else
		{
		  ri =  &right[cptr[i] * 16];
		  x2 = x2_ptr;
		  x2_ptr += 4;
		}	  	  	  	  

	      __m128d x1_0 = _mm_load_pd( &x1[0] );
	      __m128d x1_2 = _mm_load_pd( &x1[2] );
	      
	      __m128d left_k0_0 = _mm_load_pd( &le[0] );
	      __m128d left_k0_2 = _mm_load_pd( &le[2] );
	      __m128d left_k1_0 = _mm_load_pd( &le[4] );
	      __m128d left_k1_2 = _mm_load_pd( &le[6] );
	      __m128d left_k2_0 = _mm_load_pd( &le[8] );
	      __m128d left_k2_2 = _mm_load_pd( &le[10] );
	      __m128d left_k3_0 = _mm_load_pd( &le[12] );
	      __m128d left_k3_2 = _mm_load_pd( &le[14] );
	      
	      left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
	      left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
	      
	      left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
	      left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
	      
	      left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
	      left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
	      left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
	      
	      left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
	      left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
	      
	      left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
	      left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
	      
	      left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
	      left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
	      left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
	      
	      __m128d x2_0 = _mm_load_pd( &x2[0] );
	      __m128d x2_2 = _mm_load_pd( &x2[2] );
	      
	      __m128d right_k0_0 = _mm_load_pd( &ri[0] );
	      __m128d right_k0_2 = _mm_load_pd( &ri[2] );
	      __m128d right_k1_0 = _mm_load_pd( &ri[4] );
	      __m128d right_k1_2 = _mm_load_pd( &ri[6] );
	      __m128d right_k2_0 = _mm_load_pd( &ri[8] );
	      __m128d right_k2_2 = _mm_load_pd( &ri[10] );
	      __m128d right_k3_0 = _mm_load_pd( &ri[12] );
	      __m128d right_k3_2 = _mm_load_pd( &ri[14] );
	      
	      right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	      right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	  
	      right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	      right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	      
	      right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	      right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	      right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	      
	      right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	      right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	      
	      right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	      right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	      
	      right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	      right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	      right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   
	      
	      __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
	      __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );
	      
	      __m128d EV_t_l0_k0 = EVV[0];
	      __m128d EV_t_l0_k2 = EVV[1];
	      __m128d EV_t_l1_k0 = EVV[2];
	      __m128d EV_t_l1_k2 = EVV[3];
	      __m128d EV_t_l2_k0 = EVV[4];
	      __m128d EV_t_l2_k2 = EVV[5];
	      __m128d EV_t_l3_k0 = EVV[6];
	      __m128d EV_t_l3_k2 = EVV[7];
	      
	      
	      EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	      EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	      EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	      
	      EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	      EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	      
	      EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	      EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	      
	      EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	      EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	      EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	      
	      EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	      EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	      EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	      
	      EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );	  	 	    		  
	      
	      scale = 1;
	      
	      __m128d v1 = _mm_and_pd(EV_t_l0_k0, absMask.m);
	      v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	      if(_mm_movemask_pd( v1 ) != 3)
		scale = 0;
	      else
		{
		  v1 = _mm_and_pd(EV_t_l2_k0, absMask.m);
		  v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		  if(_mm_movemask_pd( v1 ) != 3)
		    scale = 0;
		}
	  	  
	      if(scale)
		{		      
		  _mm_store_pd(&x3[0], _mm_mul_pd(EV_t_l0_k0, sc));
		  _mm_store_pd(&x3[2], _mm_mul_pd(EV_t_l2_k0, sc));	      	      
		  		  
		  addScale += wgt[i];	  
		}	
	      else
		{
		  _mm_store_pd(x3, EV_t_l0_k0);
		  _mm_store_pd(&x3[2], EV_t_l2_k0);
		}

	      x3_ptr += 4;
	    }
	  	  
	}
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{ 
	  if(isGap(x3_gap, i))
	    {
	      if(scaleGap)		   		    
		addScale += wgt[i];
	    }
	  else
	    {	     
	      x3 = x3_ptr;
	  	  
	      if(isGap(x1_gap, i))
		{
		  x1 = x1_gapColumn;
		  le =  &left[maxCats * 16];
		}
	      else
		{
		  le =  &left[cptr[i] * 16];
		  x1 = x1_ptr;
		  x1_ptr += 4;
		}

	      if(isGap(x2_gap, i))	
		{
		  x2 = x2_gapColumn;
		  ri =  &right[maxCats * 16];	    
		}
	      else
		{
		  ri =  &right[cptr[i] * 16];
		  x2 = x2_ptr;
		  x2_ptr += 4;
		}	 	  	  	  

	      __m128d x1_0 = _mm_load_pd( &x1[0] );
	      __m128d x1_2 = _mm_load_pd( &x1[2] );
	      
	      __m128d left_k0_0 = _mm_load_pd( &le[0] );
	      __m128d left_k0_2 = _mm_load_pd( &le[2] );
	      __m128d left_k1_0 = _mm_load_pd( &le[4] );
	      __m128d left_k1_2 = _mm_load_pd( &le[6] );
	      __m128d left_k2_0 = _mm_load_pd( &le[8] );
	      __m128d left_k2_2 = _mm_load_pd( &le[10] );
	      __m128d left_k3_0 = _mm_load_pd( &le[12] );
	      __m128d left_k3_2 = _mm_load_pd( &le[14] );
	      
	      left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
	      left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
	      
	      left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
	      left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
	      
	      left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
	      left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
	      left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
	      
	      left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
	      left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
	      
	      left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
	      left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
	      
	      left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
	      left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
	      left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
	      
	      __m128d x2_0 = _mm_load_pd( &x2[0] );
	      __m128d x2_2 = _mm_load_pd( &x2[2] );
	      
	      __m128d right_k0_0 = _mm_load_pd( &ri[0] );
	      __m128d right_k0_2 = _mm_load_pd( &ri[2] );
	      __m128d right_k1_0 = _mm_load_pd( &ri[4] );
	      __m128d right_k1_2 = _mm_load_pd( &ri[6] );
	      __m128d right_k2_0 = _mm_load_pd( &ri[8] );
	      __m128d right_k2_2 = _mm_load_pd( &ri[10] );
	      __m128d right_k3_0 = _mm_load_pd( &ri[12] );
	      __m128d right_k3_2 = _mm_load_pd( &ri[14] );
	      
	      right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	      right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	      
	      right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	      right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	      
	      right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	      right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	      right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	      
	      right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	      right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	      
	      right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	      right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	      
	      right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	      right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	      right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   
	      
	      __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
	      __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );
	      
	      __m128d EV_t_l0_k0 = EVV[0];
	      __m128d EV_t_l0_k2 = EVV[1];
	      __m128d EV_t_l1_k0 = EVV[2];
	      __m128d EV_t_l1_k2 = EVV[3];
	      __m128d EV_t_l2_k0 = EVV[4];
	      __m128d EV_t_l2_k2 = EVV[5];
	      __m128d EV_t_l3_k0 = EVV[6];
	      __m128d EV_t_l3_k2 = EVV[7];
	      
	      
	      EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	      EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	      EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	      
	      EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	      EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	      
	      EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	      EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	      
	      EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	      EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	      EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	      
	      EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	      EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	      EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	      
	      EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );	  	 	    		  	 
	      
	      scale = 1;
	      
	      __m128d v1 = _mm_and_pd(EV_t_l0_k0, absMask.m);
	      v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	      if(_mm_movemask_pd( v1 ) != 3)
		scale = 0;
	      else
		{
		  v1 = _mm_and_pd(EV_t_l2_k0, absMask.m);
		  v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		  if(_mm_movemask_pd( v1 ) != 3)
		    scale = 0;
		}
	  	  
	      if(scale)
		{		      
		  _mm_store_pd(&x3[0], _mm_mul_pd(EV_t_l0_k0, sc));
		  _mm_store_pd(&x3[2], _mm_mul_pd(EV_t_l2_k0, sc));	      	      
		  	      
		  addScale += wgt[i];	  
		}	
	      else
		{
		  _mm_store_pd(x3, EV_t_l0_k0);
		  _mm_store_pd(&x3[2], EV_t_l2_k0);
		}
	     
	      x3_ptr += 4;
	    }
	}
      break;
    default:
      assert(0);
    }

  
  *scalerIncrement = addScale;
}

static void newviewGTRGAMMAPROT_GAPPED_SAVE(int tipCase,
					    double *x1, double *x2, double *x3, double *extEV, double *tipVector,
					    unsigned char *tipX1, unsigned char *tipX2,
					    size_t n, double *left, double *right, int *wgt, int *scalerIncrement, 
					    unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,  
					    double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn
					    )
{
  double  *uX1, *uX2, *v;
  double x1px2;
  size_t  i, j, l, k;
  
  int
    scale, 
    addScale = 0,   
    gapScaling = 0;
  
  double 
    *vl, *vr, *x1v, *x2v,
    *x1_ptr = x1,
    *x2_ptr = x2,
    *x3_ptr = x3;

  

  switch(tipCase)
    {
    case TIP_TIP:
      {
	double umpX1[1840], umpX2[1840];

	for(i = 0; i < 23; i++)
	  {
	    v = &(tipVector[20 * i]);

	    for(k = 0; k < 80; k++)
	      {
		double *ll =  &left[k * 20];
		double *rr =  &right[k * 20];
		
		__m128d umpX1v = _mm_setzero_pd();
		__m128d umpX2v = _mm_setzero_pd();

		for(l = 0; l < 20; l+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[l]);
		    umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));
		    umpX2v = _mm_add_pd(umpX2v, _mm_mul_pd(vv, _mm_load_pd(&rr[l])));					
		  }
		
		umpX1v = _mm_hadd_pd(umpX1v, umpX1v);
		umpX2v = _mm_hadd_pd(umpX2v, umpX2v);
		
		_mm_storel_pd(&umpX1[80 * i + k], umpX1v);
		_mm_storel_pd(&umpX2[80 * i + k], umpX2v);
	      }
	  }

	{
	  uX1 = &umpX1[1760];
	  uX2 = &umpX2[1760];

	  for(j = 0; j < 4; j++)
	    {
	      v = &x3_gapColumn[j * 20];

	      __m128d zero =  _mm_setzero_pd();
	      for(k = 0; k < 20; k+=2)		  		    
		_mm_store_pd(&v[k], zero);

	      for(k = 0; k < 20; k++)
		{ 
		  double *eev = &extEV[k * 20];
		  x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
		  __m128d x1px2v = _mm_set1_pd(x1px2);
		  
		  for(l = 0; l < 20; l+=2)
		    {
		      __m128d vv = _mm_load_pd(&v[l]);
		      __m128d ee = _mm_load_pd(&eev[l]);
		      
		      vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
		      
		      _mm_store_pd(&v[l], vv);
		    }
		}
	    }	   
	}	

	for(i = 0; i < n; i++)
	  {
	    if(!(x3_gap[i / 32] & mask32[i % 32]))
	      {
		uX1 = &umpX1[80 * tipX1[i]];
		uX2 = &umpX2[80 * tipX2[i]];
		
		for(j = 0; j < 4; j++)
		  {
		    v = &x3_ptr[j * 20];
		    
		    
		    __m128d zero =  _mm_setzero_pd();
		    for(k = 0; k < 20; k+=2)		  		    
		      _mm_store_pd(&v[k], zero);
		    
		    for(k = 0; k < 20; k++)
		      { 
			double *eev = &extEV[k * 20];
			x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
			__m128d x1px2v = _mm_set1_pd(x1px2);
			
			for(l = 0; l < 20; l+=2)
			  {
			    __m128d vv = _mm_load_pd(&v[l]);
			    __m128d ee = _mm_load_pd(&eev[l]);
			    
			    vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
			    
			    _mm_store_pd(&v[l], vv);
			  }
		      }
		  }	   
		x3_ptr += 80;
	      }
	  }
      }
      break;
    case TIP_INNER:
      {
	double umpX1[1840], ump_x2[20];


	for(i = 0; i < 23; i++)
	  {
	    v = &(tipVector[20 * i]);

	    for(k = 0; k < 80; k++)
	      {
		double *ll =  &left[k * 20];
				
		__m128d umpX1v = _mm_setzero_pd();
		
		for(l = 0; l < 20; l+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[l]);
		    umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));		    					
		  }
		
		umpX1v = _mm_hadd_pd(umpX1v, umpX1v);				
		_mm_storel_pd(&umpX1[80 * i + k], umpX1v);		

	      }
	  }

	{
	  uX1 = &umpX1[1760];

	  for(k = 0; k < 4; k++)
	    {
	      v = &(x2_gapColumn[k * 20]);
	       
	      for(l = 0; l < 20; l++)
		{		   
		  double *r =  &right[k * 400 + l * 20];
		  __m128d ump_x2v = _mm_setzero_pd();	    
		  
		  for(j = 0; j < 20; j+= 2)
		    {
		      __m128d vv = _mm_load_pd(&v[j]);
		      __m128d rr = _mm_load_pd(&r[j]);
		      ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(vv, rr));
		    }
		  
		  ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);
		  
		  _mm_storel_pd(&ump_x2[l], ump_x2v);		   		     
		}

	      v = &(x3_gapColumn[20 * k]);

	      __m128d zero =  _mm_setzero_pd();
	      for(l = 0; l < 20; l+=2)		  		    
		_mm_store_pd(&v[l], zero);
		  
	      for(l = 0; l < 20; l++)
		{
		  double *eev = &extEV[l * 20];
		  x1px2 = uX1[k * 20 + l]  * ump_x2[l];
		  __m128d x1px2v = _mm_set1_pd(x1px2);
		  
		  for(j = 0; j < 20; j+=2)
		    {
		      __m128d vv = _mm_load_pd(&v[j]);
		      __m128d ee = _mm_load_pd(&eev[j]);
		      
		      vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
		      
		      _mm_store_pd(&v[j], vv);
		    }		     		    
		}			
	      
	    }
	  
	  { 
	    v = x3_gapColumn;
	    __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	    
	    scale = 1;
	    for(l = 0; scale && (l < 80); l += 2)
	      {
		__m128d vv = _mm_load_pd(&v[l]);
		__m128d v1 = _mm_and_pd(vv, absMask.m);
		v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		if(_mm_movemask_pd( v1 ) != 3)
		  scale = 0;
	      }	    	  
	  }


	  if (scale)
	    {
	      gapScaling = 1;
	      __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	      
	      for(l = 0; l < 80; l+=2)
		{
		  __m128d ex3v = _mm_load_pd(&v[l]);		  
		  _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
		}		   		  	      	    	       
	    }
	}

	for (i = 0; i < n; i++)
	  {	    
	    if((x3_gap[i / 32] & mask32[i % 32]))
	       {	       
		 if(gapScaling)
		   {		     
		     addScale += wgt[i];		     
		   }
	       }
	     else
	       {
		 uX1 = &umpX1[80 * tipX1[i]];

		  if(x2_gap[i / 32] & mask32[i % 32])
		   x2v = x2_gapColumn;
		  else
		    {
		      x2v = x2_ptr;
		      x2_ptr += 80;
		    }
		 
		 for(k = 0; k < 4; k++)
		   {
		     v = &(x2v[k * 20]);
		     
		     for(l = 0; l < 20; l++)
		       {		   
			 double *r =  &right[k * 400 + l * 20];
			 __m128d ump_x2v = _mm_setzero_pd();	    
			 
			 for(j = 0; j < 20; j+= 2)
			   {
			     __m128d vv = _mm_load_pd(&v[j]);
			     __m128d rr = _mm_load_pd(&r[j]);
			     ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(vv, rr));
			   }
			 
			 ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);
			 
			 _mm_storel_pd(&ump_x2[l], ump_x2v);		   		     
		       }
		     
		     v = &x3_ptr[20 * k];
		     
		     __m128d zero =  _mm_setzero_pd();
		     for(l = 0; l < 20; l+=2)		  		    
		       _mm_store_pd(&v[l], zero);
		     
		     for(l = 0; l < 20; l++)
		       {
			 double *eev = &extEV[l * 20];
			 x1px2 = uX1[k * 20 + l]  * ump_x2[l];
			 __m128d x1px2v = _mm_set1_pd(x1px2);
			 
			 for(j = 0; j < 20; j+=2)
			   {
			     __m128d vv = _mm_load_pd(&v[j]);
			     __m128d ee = _mm_load_pd(&eev[j]);
			     
			     vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
			     
			     _mm_store_pd(&v[j], vv);
			   }		     		    
		       }			
		     
		   }
		 
		 
		 { 
		   v = x3_ptr;
		   __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
		   
		   scale = 1;
		   for(l = 0; scale && (l < 80); l += 2)
		     {
		       __m128d vv = _mm_load_pd(&v[l]);
		       __m128d v1 = _mm_and_pd(vv, absMask.m);
		       v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		       if(_mm_movemask_pd( v1 ) != 3)
			 scale = 0;
		     }	    	  
		 }
		 
		 
		 if (scale)
		   {
		     __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
		     
		     for(l = 0; l < 80; l+=2)
		       {
			 __m128d ex3v = _mm_load_pd(&v[l]);		  
			 _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
		       }		   		  
		     		    
		     addScale += wgt[i];		      
		   }
		 
		 x3_ptr += 80;
	       }
	  }
      }
      break;
    case INNER_INNER:
      {
	for(k = 0; k < 4; k++)
	   {
	     vl = &(x1_gapColumn[20 * k]);
	     vr = &(x2_gapColumn[20 * k]);
	     v =  &(x3_gapColumn[20 * k]);

	     __m128d zero =  _mm_setzero_pd();
	     for(l = 0; l < 20; l+=2)		  		    
	       _mm_store_pd(&v[l], zero);
	     
	     for(l = 0; l < 20; l++)
	       {		 
		 {
		   __m128d al = _mm_setzero_pd();
		   __m128d ar = _mm_setzero_pd();

		   double *ll   = &left[k * 400 + l * 20];
		   double *rr   = &right[k * 400 + l * 20];
		   double *EVEV = &extEV[20 * l];
		   
		   for(j = 0; j < 20; j+=2)
		     {
		       __m128d lv  = _mm_load_pd(&ll[j]);
		       __m128d rv  = _mm_load_pd(&rr[j]);
		       __m128d vll = _mm_load_pd(&vl[j]);
		       __m128d vrr = _mm_load_pd(&vr[j]);
		       
		       al = _mm_add_pd(al, _mm_mul_pd(vll, lv));
		       ar = _mm_add_pd(ar, _mm_mul_pd(vrr, rv));
		     }  		 
		       
		   al = _mm_hadd_pd(al, al);
		   ar = _mm_hadd_pd(ar, ar);
		   
		   al = _mm_mul_pd(al, ar);

		   for(j = 0; j < 20; j+=2)
		     {
		       __m128d vv  = _mm_load_pd(&v[j]);
		       __m128d EVV = _mm_load_pd(&EVEV[j]);

		       vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));

		       _mm_store_pd(&v[j], vv);
		     }		  		   		  
		 }		 

	       }
	   }
	 

	{ 
	   v = x3_gapColumn;
	   __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	   
	   scale = 1;
	   for(l = 0; scale && (l < 80); l += 2)
	     {
	       __m128d vv = _mm_load_pd(&v[l]);
	       __m128d v1 = _mm_and_pd(vv, absMask.m);
	       v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	       if(_mm_movemask_pd( v1 ) != 3)
		 scale = 0;
	     }	    	  
	 }

	 if (scale)
	   {
	     gapScaling = 1;
	     __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	     
	     for(l = 0; l < 80; l+=2)
	       {
		 __m128d ex3v = _mm_load_pd(&v[l]);		  
		 _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
	       }		   		  
	     
	    	  
	   }
      }

      for (i = 0; i < n; i++)
       {
	  if(x3_gap[i / 32] & mask32[i % 32])
	   {	     
	     if(gapScaling)
	       {		
		 addScale += wgt[i];			       
	       }
	   }
	 else
	   {
	     if(x1_gap[i / 32] & mask32[i % 32])
	       x1v = x1_gapColumn;
	     else
	       {
		 x1v = x1_ptr;
		 x1_ptr += 80;
	       }

	     if(x2_gap[i / 32] & mask32[i % 32])
	       x2v = x2_gapColumn;
	     else
	       {
		 x2v = x2_ptr;
		 x2_ptr += 80;
	       }

	     for(k = 0; k < 4; k++)
	       {
		 vl = &(x1v[20 * k]);
		 vr = &(x2v[20 * k]);
		 v =  &x3_ptr[20 * k];
		 		 
		 __m128d zero =  _mm_setzero_pd();
		 for(l = 0; l < 20; l+=2)		  		    
		   _mm_store_pd(&v[l], zero);
		 		 
		 for(l = 0; l < 20; l++)
		   {		 
		     {
		       __m128d al = _mm_setzero_pd();
		       __m128d ar = _mm_setzero_pd();
		       
		       double *ll   = &left[k * 400 + l * 20];
		       double *rr   = &right[k * 400 + l * 20];
		       double *EVEV = &extEV[20 * l];
		       
		       for(j = 0; j < 20; j+=2)
			 {
			   __m128d lv  = _mm_load_pd(&ll[j]);
			   __m128d rv  = _mm_load_pd(&rr[j]);
			   __m128d vll = _mm_load_pd(&vl[j]);
			   __m128d vrr = _mm_load_pd(&vr[j]);
			   
			   al = _mm_add_pd(al, _mm_mul_pd(vll, lv));
			   ar = _mm_add_pd(ar, _mm_mul_pd(vrr, rv));
			 }  		 
		       
		       al = _mm_hadd_pd(al, al);
		       ar = _mm_hadd_pd(ar, ar);
		       
		       al = _mm_mul_pd(al, ar);
		       
		       for(j = 0; j < 20; j+=2)
			 {
			   __m128d vv  = _mm_load_pd(&v[j]);
			   __m128d EVV = _mm_load_pd(&EVEV[j]);
			   
			   vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
			   
			   _mm_store_pd(&v[j], vv);
			 }		  		   		  
		     }		 
		     
		   }
	       }
	     

	     
	     { 
	       v = x3_ptr;
	       __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	       
	       scale = 1;
	       for(l = 0; scale && (l < 80); l += 2)
		 {
		   __m128d vv = _mm_load_pd(&v[l]);
		   __m128d v1 = _mm_and_pd(vv, absMask.m);
		   v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		   if(_mm_movemask_pd( v1 ) != 3)
		     scale = 0;
		 }	    	  
	     }
	     
	     
	     if (scale)
	       {
		 __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
		 
		 for(l = 0; l < 80; l+=2)
		   {
		     __m128d ex3v = _mm_load_pd(&v[l]);		  
		     _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
		   }		   		  
		 		
		 addScale += wgt[i];		 	  
	       }
	     x3_ptr += 80;
	   }
       }
      break;
    default:
      assert(0);
    }

 
  *scalerIncrement = addScale;  
}




static void newviewGTRGAMMAPROT(int tipCase,
				double *x1, double *x2, double *x3, double *extEV, double *tipVector,
				unsigned char *tipX1, unsigned char *tipX2,
				size_t n, double *left, double *right, int *wgt, int *scalerIncrement)
{
  double  *uX1, *uX2, *v;
  double x1px2;
  size_t  i, j, l, k;
  int 
    scale, addScale = 0;
  double *vl, *vr;



  switch(tipCase)
    {
    case TIP_TIP:
      {
	double umpX1[1840], umpX2[1840];

	for(i = 0; i < 23; i++)
	  {
	    v = &(tipVector[20 * i]);

	    for(k = 0; k < 80; k++)
	      {
		double *ll =  &left[k * 20];
		double *rr =  &right[k * 20];
		
		__m128d umpX1v = _mm_setzero_pd();
		__m128d umpX2v = _mm_setzero_pd();

		for(l = 0; l < 20; l+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[l]);
		    umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));
		    umpX2v = _mm_add_pd(umpX2v, _mm_mul_pd(vv, _mm_load_pd(&rr[l])));					
		  }
		
		umpX1v = _mm_hadd_pd(umpX1v, umpX1v);
		umpX2v = _mm_hadd_pd(umpX2v, umpX2v);
		
		_mm_storel_pd(&umpX1[80 * i + k], umpX1v);
		_mm_storel_pd(&umpX2[80 * i + k], umpX2v);

	      }
	  }

	for(i = 0; i < n; i++)
	  {
	    uX1 = &umpX1[80 * tipX1[i]];
	    uX2 = &umpX2[80 * tipX2[i]];

	    for(j = 0; j < 4; j++)
	      {
		v = &x3[i * 80 + j * 20];


		__m128d zero =  _mm_setzero_pd();
		for(k = 0; k < 20; k+=2)		  		    
		  _mm_store_pd(&v[k], zero);

		for(k = 0; k < 20; k++)
		  { 
		    double *eev = &extEV[k * 20];
		    x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
		    __m128d x1px2v = _mm_set1_pd(x1px2);

		    for(l = 0; l < 20; l+=2)
		      {
		      	__m128d vv = _mm_load_pd(&v[l]);
			__m128d ee = _mm_load_pd(&eev[l]);

			vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
			
			_mm_store_pd(&v[l], vv);
		      }
		  }


	      }	   
	  }
      }
      break;
    case TIP_INNER:
      {
	double umpX1[1840], ump_x2[20];


	for(i = 0; i < 23; i++)
	  {
	    v = &(tipVector[20 * i]);

	    for(k = 0; k < 80; k++)
	      {
		double *ll =  &left[k * 20];
				
		__m128d umpX1v = _mm_setzero_pd();
		
		for(l = 0; l < 20; l+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[l]);
		    umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));		    					
		  }
		
		umpX1v = _mm_hadd_pd(umpX1v, umpX1v);				
		_mm_storel_pd(&umpX1[80 * i + k], umpX1v);		


	      }
	  }

	for (i = 0; i < n; i++)
	  {
	    uX1 = &umpX1[80 * tipX1[i]];

	    for(k = 0; k < 4; k++)
	      {
		v = &(x2[80 * i + k * 20]);
	       
		for(l = 0; l < 20; l++)
		  {		   
		    double *r =  &right[k * 400 + l * 20];
		    __m128d ump_x2v = _mm_setzero_pd();	    
		    
		    for(j = 0; j < 20; j+= 2)
		      {
			__m128d vv = _mm_load_pd(&v[j]);
			__m128d rr = _mm_load_pd(&r[j]);
			ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(vv, rr));
		      }
		     
		    ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);
		    
		    _mm_storel_pd(&ump_x2[l], ump_x2v);		   		     
		  }

		v = &(x3[80 * i + 20 * k]);

		__m128d zero =  _mm_setzero_pd();
		for(l = 0; l < 20; l+=2)		  		    
		  _mm_store_pd(&v[l], zero);
		  
		for(l = 0; l < 20; l++)
		  {
		    double *eev = &extEV[l * 20];
		    x1px2 = uX1[k * 20 + l]  * ump_x2[l];
		    __m128d x1px2v = _mm_set1_pd(x1px2);
		  
		    for(j = 0; j < 20; j+=2)
		      {
			__m128d vv = _mm_load_pd(&v[j]);
			__m128d ee = _mm_load_pd(&eev[j]);
			
			vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
			
			_mm_store_pd(&v[j], vv);
		      }		     		    
		  }			

	      }
	   

	    { 
	      v = &(x3[80 * i]);
	      __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	      
	      scale = 1;
	      for(l = 0; scale && (l < 80); l += 2)
		{
		  __m128d vv = _mm_load_pd(&v[l]);
		  __m128d v1 = _mm_and_pd(vv, absMask.m);
		  v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		  if(_mm_movemask_pd( v1 ) != 3)
		    scale = 0;
		}	    	  
	    }


	    if (scale)
	      {

	       __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	       
	       for(l = 0; l < 80; l+=2)
		 {
		   __m128d ex3v = _mm_load_pd(&v[l]);		  
		   _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
		 }		   		  


	
		addScale += wgt[i];
		       
	      }
	  }
      }
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
       {
	 for(k = 0; k < 4; k++)
	   {
	     vl = &(x1[80 * i + 20 * k]);
	     vr = &(x2[80 * i + 20 * k]);
	     v =  &(x3[80 * i + 20 * k]);


	     __m128d zero =  _mm_setzero_pd();
	     for(l = 0; l < 20; l+=2)		  		    
	       _mm_store_pd(&v[l], zero);


	     for(l = 0; l < 20; l++)
	       {		 

		 {
		   __m128d al = _mm_setzero_pd();
		   __m128d ar = _mm_setzero_pd();

		   double *ll   = &left[k * 400 + l * 20];
		   double *rr   = &right[k * 400 + l * 20];
		   double *EVEV = &extEV[20 * l];
		   
		   for(j = 0; j < 20; j+=2)
		     {
		       __m128d lv  = _mm_load_pd(&ll[j]);
		       __m128d rv  = _mm_load_pd(&rr[j]);
		       __m128d vll = _mm_load_pd(&vl[j]);
		       __m128d vrr = _mm_load_pd(&vr[j]);
		       
		       al = _mm_add_pd(al, _mm_mul_pd(vll, lv));
		       ar = _mm_add_pd(ar, _mm_mul_pd(vrr, rv));
		     }  		 
		       
		   al = _mm_hadd_pd(al, al);
		   ar = _mm_hadd_pd(ar, ar);
		   
		   al = _mm_mul_pd(al, ar);

		   for(j = 0; j < 20; j+=2)
		     {
		       __m128d vv  = _mm_load_pd(&v[j]);
		       __m128d EVV = _mm_load_pd(&EVEV[j]);

		       vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));

		       _mm_store_pd(&v[j], vv);
		     }		  		   		  
		 }		 

	       }
	   }
	 


	 { 
	   v = &(x3[80 * i]);
	   __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	   
	   scale = 1;
	   for(l = 0; scale && (l < 80); l += 2)
	     {
	       __m128d vv = _mm_load_pd(&v[l]);
	       __m128d v1 = _mm_and_pd(vv, absMask.m);
	       v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	       if(_mm_movemask_pd( v1 ) != 3)
		 scale = 0;
	     }	    	  
	 }


	 if (scale)
	   {

	       __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	       
	       for(l = 0; l < 80; l+=2)
		 {
		   __m128d ex3v = _mm_load_pd(&v[l]);		  
		   _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
		 }		   		  


	    
	     addScale += wgt[i];
	      
	   }
       }
      break;
    default:
      assert(0);
    }

  
  *scalerIncrement = addScale;

}


     
static void newviewGTRCATPROT(int tipCase, double *extEV,
			      int *cptr,
			      double *x1, double *x2, double *x3, double *tipVector,
			      unsigned char *tipX1, unsigned char *tipX2,
			      size_t n, double *left, double *right, int *wgt, int *scalerIncrement )
{
  double
    *le, *ri, *v, *vl, *vr;

  size_t i, l, j;
  int 
    scale, addScale = 0;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for (i = 0; i < n; i++)
	  {
	    le = &left[cptr[i] * 400];
	    ri = &right[cptr[i] * 400];

	    vl = &(tipVector[20 * tipX1[i]]);
	    vr = &(tipVector[20 * tipX2[i]]);
	    v  = &x3[20 * i];

	    for(l = 0; l < 20; l+=2)
	      _mm_store_pd(&v[l], _mm_setzero_pd());	      		


	    for(l = 0; l < 20; l++)
	      {
		__m128d x1v = _mm_setzero_pd();
		__m128d x2v = _mm_setzero_pd();	 
		double 
		  *ev = &extEV[l * 20],
		  *lv = &le[l * 20],
		  *rv = &ri[l * 20];

		for(j = 0; j < 20; j+=2)
		  {
		    x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));		    
		    x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
		  }

		x1v = _mm_hadd_pd(x1v, x1v);
		x2v = _mm_hadd_pd(x2v, x2v);

		x1v = _mm_mul_pd(x1v, x2v);
		
		for(j = 0; j < 20; j+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[j]);
		    vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
		    _mm_store_pd(&v[j], vv);
		  }		    

	      }	   
	  }
      }
      break;
    case TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    le = &left[cptr[i] * 400];
	    ri = &right[cptr[i] * 400];

	    vl = &(tipVector[20 * tipX1[i]]);
	    vr = &x2[20 * i];
	    v  = &x3[20 * i];

	    for(l = 0; l < 20; l+=2)
	      _mm_store_pd(&v[l], _mm_setzero_pd());	      		

	   

	    for(l = 0; l < 20; l++)
	      {

		__m128d x1v = _mm_setzero_pd();
		__m128d x2v = _mm_setzero_pd();	
		double 
		  *ev = &extEV[l * 20],
		  *lv = &le[l * 20],
		  *rv = &ri[l * 20];

		for(j = 0; j < 20; j+=2)
		  {
		    x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));		    
		    x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
		  }

		x1v = _mm_hadd_pd(x1v, x1v);
		x2v = _mm_hadd_pd(x2v, x2v);

		x1v = _mm_mul_pd(x1v, x2v);
		
		for(j = 0; j < 20; j+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[j]);
		    vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
		    _mm_store_pd(&v[j], vv);
		  }		    

	      }

	    { 	    
	      __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	      
	      scale = 1;
	      for(l = 0; scale && (l < 20); l += 2)
		{
		  __m128d vv = _mm_load_pd(&v[l]);
		  __m128d v1 = _mm_and_pd(vv, absMask.m);
		  v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		  if(_mm_movemask_pd( v1 ) != 3)
		    scale = 0;
		}	    	  
	    }


	    if(scale)
	      {

		__m128d twoto = _mm_set_pd(twotothe256, twotothe256);

		for(l = 0; l < 20; l+=2)
		  {
		    __m128d ex3v = _mm_load_pd(&v[l]);
		    _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));		    
		  }
	
		addScale += wgt[i];	  
	      }
	  }
      }
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  le = &left[cptr[i] * 400];
	  ri = &right[cptr[i] * 400];

	  vl = &x1[20 * i];
	  vr = &x2[20 * i];
	  v = &x3[20 * i];


	    for(l = 0; l < 20; l+=2)
	      _mm_store_pd(&v[l], _mm_setzero_pd());	      		

	 
	  for(l = 0; l < 20; l++)
	    {

		__m128d x1v = _mm_setzero_pd();
		__m128d x2v = _mm_setzero_pd();
		double 
		  *ev = &extEV[l * 20],
		  *lv = &le[l * 20],
		  *rv = &ri[l * 20];


		for(j = 0; j < 20; j+=2)
		  {
		    x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));		    
		    x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
		  }

		x1v = _mm_hadd_pd(x1v, x1v);
		x2v = _mm_hadd_pd(x2v, x2v);

		x1v = _mm_mul_pd(x1v, x2v);
		
		for(j = 0; j < 20; j+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[j]);
		    vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
		    _mm_store_pd(&v[j], vv);
		  }		    

	    }

	    { 	    
	      __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	      
	      scale = 1;
	      for(l = 0; scale && (l < 20); l += 2)
		{
		  __m128d vv = _mm_load_pd(&v[l]);
		  __m128d v1 = _mm_and_pd(vv, absMask.m);
		  v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		  if(_mm_movemask_pd( v1 ) != 3)
		    scale = 0;
		}	    	  
	    }
   

	   if(scale)
	     {

	       __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	       
	       for(l = 0; l < 20; l+=2)
		 {
		   __m128d ex3v = _mm_load_pd(&v[l]);		  
		   _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
		 }		   		  


	       
	       addScale += wgt[i];	   
	     }
	}
      break;
    default:
      assert(0);
    }
  
 
  *scalerIncrement = addScale;

}

static void newviewGTRCATPROT_SAVE(int tipCase, double *extEV,
				   int *cptr,
				   double *x1, double *x2, double *x3, double *tipVector,
				   unsigned char *tipX1, unsigned char *tipX2,
				   size_t n, double *left, double *right, int *wgt, int *scalerIncrement,
				   unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
				   double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats)
{
  double
    *le, 
    *ri, 
    *v, 
    *vl, 
    *vr,
    *x1_ptr = x1,
    *x2_ptr = x2, 
    *x3_ptr = x3;

  size_t
    i, 
    l, 
    j;

  int
    scale, 
    scaleGap = 0,
    addScale = 0;

  {
    vl = x1_gapColumn;	      
    vr = x2_gapColumn;
    v = x3_gapColumn;

    le = &left[maxCats * 400];
    ri = &right[maxCats * 400];	  

    for(l = 0; l < 20; l+=2)
      _mm_store_pd(&v[l], _mm_setzero_pd());	      		
	 
    for(l = 0; l < 20; l++)
      {
	__m128d x1v = _mm_setzero_pd();
	__m128d x2v = _mm_setzero_pd();
	double 
	  *ev = &extEV[l * 20],
	  *lv = &le[l * 20],
	  *rv = &ri[l * 20];


	for(j = 0; j < 20; j+=2)
	  {
	    x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));		    
	    x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
	  }
	
	x1v = _mm_hadd_pd(x1v, x1v);
	x2v = _mm_hadd_pd(x2v, x2v);
	
	x1v = _mm_mul_pd(x1v, x2v);
	
	for(j = 0; j < 20; j+=2)
	  {
	    __m128d vv = _mm_load_pd(&v[j]);
	    vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
	    _mm_store_pd(&v[j], vv);
	  }		    	
      }
    
    if(tipCase != TIP_TIP)
      { 	    
	__m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	      
	scale = 1;
	for(l = 0; scale && (l < 20); l += 2)
	  {
	    __m128d vv = _mm_load_pd(&v[l]);
	    __m128d v1 = _mm_and_pd(vv, absMask.m);
	    v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	    if(_mm_movemask_pd( v1 ) != 3)
	      scale = 0;
	  }	    	        
  
	if(scale)
	  {
	    __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	    
	    for(l = 0; l < 20; l+=2)
	      {
		__m128d ex3v = _mm_load_pd(&v[l]);		  
		_mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
	      }		   		  
	       
	    scaleGap = TRUE;	   
	  }
      }
  }
  
  switch(tipCase)
    {
    case TIP_TIP:
      {
	for (i = 0; i < n; i++)
	  {
	    if(noGap(x3_gap, i))
	      {		
		vl = &(tipVector[20 * tipX1[i]]);
		vr = &(tipVector[20 * tipX2[i]]);
		v  = x3_ptr;

		if(isGap(x1_gap, i))
		  le =  &left[maxCats * 400];
		else	  	  
		  le =  &left[cptr[i] * 400];	  
	  
		if(isGap(x2_gap, i))
		  ri =  &right[maxCats * 400];
		else	 	  
		  ri =  &right[cptr[i] * 400];

		for(l = 0; l < 20; l+=2)
		  _mm_store_pd(&v[l], _mm_setzero_pd());	      		
		
		for(l = 0; l < 20; l++)
		  {
		    __m128d x1v = _mm_setzero_pd();
		    __m128d x2v = _mm_setzero_pd();	 
		    double 
		      *ev = &extEV[l * 20],
		      *lv = &le[l * 20],
		      *rv = &ri[l * 20];
		    
		    for(j = 0; j < 20; j+=2)
		      {
			x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));		    
			x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
		      }
		    
		    x1v = _mm_hadd_pd(x1v, x1v);
		    x2v = _mm_hadd_pd(x2v, x2v);
		    
		    x1v = _mm_mul_pd(x1v, x2v);
		    
		    for(j = 0; j < 20; j+=2)
		      {
			__m128d vv = _mm_load_pd(&v[j]);
			vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
			_mm_store_pd(&v[j], vv);
		      }		   
		  }

		x3_ptr += 20;

	      }   
	  }
      }
      break;
    case TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    if(isGap(x3_gap, i))
	      {
		if(scaleGap)		   		    
		  addScale += wgt[i];
	      }
	    else
	      {	 
		vl = &(tipVector[20 * tipX1[i]]);
	      	       
		v = x3_ptr;

		if(isGap(x1_gap, i))
		  le =  &left[maxCats * 400];
		else
		  le =  &left[cptr[i] * 400];

		if(isGap(x2_gap, i))
		  {		 
		    ri =  &right[maxCats * 400];
		    vr = x2_gapColumn;
		  }
		else
		  {
		    ri =  &right[cptr[i] * 400];
		    vr = x2_ptr;
		    x2_ptr += 20;
		  }	  	  	  	  		  

		for(l = 0; l < 20; l+=2)
		  _mm_store_pd(&v[l], _mm_setzero_pd());	      			   

		for(l = 0; l < 20; l++)
		  {
		    __m128d x1v = _mm_setzero_pd();
		    __m128d x2v = _mm_setzero_pd();	
		    double 
		      *ev = &extEV[l * 20],
		      *lv = &le[l * 20],
		      *rv = &ri[l * 20];
		    
		    for(j = 0; j < 20; j+=2)
		      {
			x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));		    
			x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
		      }
		    
		    x1v = _mm_hadd_pd(x1v, x1v);
		    x2v = _mm_hadd_pd(x2v, x2v);
		    
		    x1v = _mm_mul_pd(x1v, x2v);
		    
		    for(j = 0; j < 20; j+=2)
		      {
			__m128d vv = _mm_load_pd(&v[j]);
			vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
			_mm_store_pd(&v[j], vv);
		      }		    
		  }
		
		{ 	    
		  __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
		  
		  scale = 1;
		  for(l = 0; scale && (l < 20); l += 2)
		    {
		      __m128d vv = _mm_load_pd(&v[l]);
		      __m128d v1 = _mm_and_pd(vv, absMask.m);
		      v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		      if(_mm_movemask_pd( v1 ) != 3)
			scale = 0;
		    }	    	  
		}
		
		
		if(scale)
		  {
		    __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
		    
		    for(l = 0; l < 20; l+=2)
		      {
			__m128d ex3v = _mm_load_pd(&v[l]);
			_mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));		    
		      }
		    
		    addScale += wgt[i];	  
		  }
		x3_ptr += 20;
	      }
	  }
      }
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
	{ 
	  if(isGap(x3_gap, i))
	    {
	      if(scaleGap)		   		    
		addScale += wgt[i];
	    }
	  else
	    {	  	     
	      v = x3_ptr;
	  	  
	      if(isGap(x1_gap, i))
		{
		  vl = x1_gapColumn;
		  le =  &left[maxCats * 400];
		}
	      else
		{
		  le =  &left[cptr[i] * 400];
		  vl = x1_ptr;
		  x1_ptr += 20;
		}

	      if(isGap(x2_gap, i))	
		{
		  vr = x2_gapColumn;
		  ri =  &right[maxCats * 400];	    
		}
	      else
		{
		  ri =  &right[cptr[i] * 400];
		  vr = x2_ptr;
		  x2_ptr += 20;
		}	 	  	  	  

	      for(l = 0; l < 20; l+=2)
		_mm_store_pd(&v[l], _mm_setzero_pd());	      		
	 
	      for(l = 0; l < 20; l++)
		{
		  __m128d x1v = _mm_setzero_pd();
		  __m128d x2v = _mm_setzero_pd();
		  double 
		    *ev = &extEV[l * 20],
		    *lv = &le[l * 20],
		    *rv = &ri[l * 20];
		  		  
		  for(j = 0; j < 20; j+=2)
		    {
		      x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));		    
		      x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
		    }
		  
		  x1v = _mm_hadd_pd(x1v, x1v);
		  x2v = _mm_hadd_pd(x2v, x2v);
		  
		  x1v = _mm_mul_pd(x1v, x2v);
		  
		  for(j = 0; j < 20; j+=2)
		    {
		      __m128d vv = _mm_load_pd(&v[j]);
		      vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
		      _mm_store_pd(&v[j], vv);
		    }		    
		  
		}
	      
	      { 	    
		__m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
		
		scale = 1;
		for(l = 0; scale && (l < 20); l += 2)
		  {
		    __m128d vv = _mm_load_pd(&v[l]);
		    __m128d v1 = _mm_and_pd(vv, absMask.m);
		    v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		    if(_mm_movemask_pd( v1 ) != 3)
		      scale = 0;
		  }	    	  
	      }
  
	      if(scale)
		{
		  __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
		  
		  for(l = 0; l < 20; l+=2)
		    {
		      __m128d ex3v = _mm_load_pd(&v[l]);		  
		      _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
		    }		   		  
		  
		  addScale += wgt[i];	   
		}
	      x3_ptr += 20;
	    }
	}
      break;
    default:
      assert(0);
    }
  
 
  *scalerIncrement = addScale;

}

static void newviewGTRGAMMAPROT_LG4(int tipCase,
				    double *x1, double *x2, double *x3, double *extEV[4], double *tipVector[4],
				    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				    size_t n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double  *uX1, *uX2, *v;
  double x1px2;
  size_t  i, j, l, k;
  
  int 
    scale, addScale = 0;
  double *vl, *vr;
#ifndef __SIM_SSE3
  double al, ar;
#endif



  switch(tipCase)
    {
    case TIP_TIP:
      {
	double umpX1[1840], umpX2[1840];

	for(i = 0; i < 23; i++)
	  {
	   

	    for(k = 0; k < 80; k++)
	      {
		
		v = &(tipVector[k / 20][20 * i]);
#ifdef __SIM_SSE3
		double *ll =  &left[k * 20];
		double *rr =  &right[k * 20];
		
		__m128d umpX1v = _mm_setzero_pd();
		__m128d umpX2v = _mm_setzero_pd();

		for(l = 0; l < 20; l+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[l]);
		    umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));
		    umpX2v = _mm_add_pd(umpX2v, _mm_mul_pd(vv, _mm_load_pd(&rr[l])));					
		  }
		
		umpX1v = _mm_hadd_pd(umpX1v, umpX1v);
		umpX2v = _mm_hadd_pd(umpX2v, umpX2v);
		
		_mm_storel_pd(&umpX1[80 * i + k], umpX1v);
		_mm_storel_pd(&umpX2[80 * i + k], umpX2v);
#else
		umpX1[80 * i + k] = 0.0;
		umpX2[80 * i + k] = 0.0;

		for(l = 0; l < 20; l++)
		  {
		    umpX1[80 * i + k] +=  v[l] *  left[k * 20 + l];
		    umpX2[80 * i + k] +=  v[l] * right[k * 20 + l];
		  }
#endif
	      }
	  }

	for(i = 0; i < n; i++)
	  {
	    uX1 = &umpX1[80 * tipX1[i]];
	    uX2 = &umpX2[80 * tipX2[i]];

	    for(j = 0; j < 4; j++)
	      {
		v = &x3[i * 80 + j * 20];

#ifdef __SIM_SSE3
		__m128d zero =  _mm_setzero_pd();
		for(k = 0; k < 20; k+=2)		  		    
		  _mm_store_pd(&v[k], zero);

		for(k = 0; k < 20; k++)
		  { 
		    double *eev = &extEV[j][k * 20];
		    x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
		    __m128d x1px2v = _mm_set1_pd(x1px2);

		    for(l = 0; l < 20; l+=2)
		      {
		      	__m128d vv = _mm_load_pd(&v[l]);
			__m128d ee = _mm_load_pd(&eev[l]);

			vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
			
			_mm_store_pd(&v[l], vv);
		      }
		  }

#else

		for(k = 0; k < 20; k++)
		  v[k] = 0.0;

		for(k = 0; k < 20; k++)
		  {		   
		    x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
		   
		    for(l = 0; l < 20; l++)		      					
		      v[l] += x1px2 * extEV[j][20 * k + l];		     
		  }
#endif
	      }	   
	  }
      }
      break;
    case TIP_INNER:
      {
	double umpX1[1840], ump_x2[20];


	for(i = 0; i < 23; i++)
	  {
	   

	    for(k = 0; k < 80; k++)
	      { 
		v = &(tipVector[k / 20][20 * i]);
#ifdef __SIM_SSE3
		double *ll =  &left[k * 20];
				
		__m128d umpX1v = _mm_setzero_pd();
		
		for(l = 0; l < 20; l+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[l]);
		    umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));		    					
		  }
		
		umpX1v = _mm_hadd_pd(umpX1v, umpX1v);				
		_mm_storel_pd(&umpX1[80 * i + k], umpX1v);		
#else	    
		umpX1[80 * i + k] = 0.0;

		for(l = 0; l < 20; l++)
		  umpX1[80 * i + k] +=  v[l] * left[k * 20 + l];
#endif

	      }
	  }

	for (i = 0; i < n; i++)
	  {
	    uX1 = &umpX1[80 * tipX1[i]];

	    for(k = 0; k < 4; k++)
	      {
		v = &(x2[80 * i + k * 20]);
#ifdef __SIM_SSE3	       
		for(l = 0; l < 20; l++)
		  {		   
		    double *r =  &right[k * 400 + l * 20];
		    __m128d ump_x2v = _mm_setzero_pd();	    
		    
		    for(j = 0; j < 20; j+= 2)
		      {
			__m128d vv = _mm_load_pd(&v[j]);
			__m128d rr = _mm_load_pd(&r[j]);
			ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(vv, rr));
		      }
		     
		    ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);
		    
		    _mm_storel_pd(&ump_x2[l], ump_x2v);		   		     
		  }

		v = &(x3[80 * i + 20 * k]);

		__m128d zero =  _mm_setzero_pd();
		for(l = 0; l < 20; l+=2)		  		    
		  _mm_store_pd(&v[l], zero);
		  
		for(l = 0; l < 20; l++)
		  {
		    double *eev = &extEV[k][l * 20];
		    x1px2 = uX1[k * 20 + l]  * ump_x2[l];
		    __m128d x1px2v = _mm_set1_pd(x1px2);
		  
		    for(j = 0; j < 20; j+=2)
		      {
			__m128d vv = _mm_load_pd(&v[j]);
			__m128d ee = _mm_load_pd(&eev[j]);
			
			vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
			
			_mm_store_pd(&v[j], vv);
		      }		     		    
		  }			
#else
		for(l = 0; l < 20; l++)
		  {
		    ump_x2[l] = 0.0;

		    for(j = 0; j < 20; j++)
		      ump_x2[l] += v[j] * right[k * 400 + l * 20 + j];
		  }

		v = &(x3[80 * i + 20 * k]);

		for(l = 0; l < 20; l++)
		  v[l] = 0;

		for(l = 0; l < 20; l++)
		  {
		    x1px2 = uX1[k * 20 + l]  * ump_x2[l];
		    for(j = 0; j < 20; j++)
		      v[j] += x1px2 * extEV[k][l * 20  + j];
		  }
#endif
	      }
	   
#ifdef __SIM_SSE3
	    { 
	      v = &(x3[80 * i]);
	      __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	      
	      scale = 1;
	      for(l = 0; scale && (l < 80); l += 2)
		{
		  __m128d vv = _mm_load_pd(&v[l]);
		  __m128d v1 = _mm_and_pd(vv, absMask.m);
		  v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		  if(_mm_movemask_pd( v1 ) != 3)
		    scale = 0;
		}	    	  
	    }
#else
	    v = &x3[80 * i];
	    scale = 1;
	    for(l = 0; scale && (l < 80); l++)
	      scale = (ABS(v[l]) <  minlikelihood);
#endif

	    if (scale)
	      {
#ifdef __SIM_SSE3
	       __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	       
	       for(l = 0; l < 80; l+=2)
		 {
		   __m128d ex3v = _mm_load_pd(&v[l]);		  
		   _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
		 }		   		  
#else
		for(l = 0; l < 80; l++)
		  v[l] *= twotothe256;
#endif

		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i]  += 1;	       
	      }
	  }
      }
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
       {
	 for(k = 0; k < 4; k++)
	   {
	     vl = &(x1[80 * i + 20 * k]);
	     vr = &(x2[80 * i + 20 * k]);
	     v =  &(x3[80 * i + 20 * k]);

#ifdef __SIM_SSE3
	     __m128d zero =  _mm_setzero_pd();
	     for(l = 0; l < 20; l+=2)		  		    
	       _mm_store_pd(&v[l], zero);
#else
	     for(l = 0; l < 20; l++)
	       v[l] = 0;
#endif

	     for(l = 0; l < 20; l++)
	       {		 
#ifdef __SIM_SSE3
		 {
		   __m128d al = _mm_setzero_pd();
		   __m128d ar = _mm_setzero_pd();

		   double *ll   = &left[k * 400 + l * 20];
		   double *rr   = &right[k * 400 + l * 20];
		   double *EVEV = &extEV[k][20 * l];
		   
		   for(j = 0; j < 20; j+=2)
		     {
		       __m128d lv  = _mm_load_pd(&ll[j]);
		       __m128d rv  = _mm_load_pd(&rr[j]);
		       __m128d vll = _mm_load_pd(&vl[j]);
		       __m128d vrr = _mm_load_pd(&vr[j]);
		       
		       al = _mm_add_pd(al, _mm_mul_pd(vll, lv));
		       ar = _mm_add_pd(ar, _mm_mul_pd(vrr, rv));
		     }  		 
		       
		   al = _mm_hadd_pd(al, al);
		   ar = _mm_hadd_pd(ar, ar);
		   
		   al = _mm_mul_pd(al, ar);

		   for(j = 0; j < 20; j+=2)
		     {
		       __m128d vv  = _mm_load_pd(&v[j]);
		       __m128d EVV = _mm_load_pd(&EVEV[j]);

		       vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));

		       _mm_store_pd(&v[j], vv);
		     }		  		   		  
		 }		 
#else
		 al = 0.0;
		 ar = 0.0;

		 for(j = 0; j < 20; j++)
		   {
		     al += vl[j] * left[k * 400 + l * 20 + j];
		     ar += vr[j] * right[k * 400 + l * 20 + j];
		   }

		 x1px2 = al * ar;

		 for(j = 0; j < 20; j++)
		   v[j] += x1px2 * extEV[k][20 * l + j];
#endif
	       }
	   }
	 

#ifdef __SIM_SSE3
	 { 
	   v = &(x3[80 * i]);
	   __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	   
	   scale = 1;
	   for(l = 0; scale && (l < 80); l += 2)
	     {
	       __m128d vv = _mm_load_pd(&v[l]);
	       __m128d v1 = _mm_and_pd(vv, absMask.m);
	       v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	       if(_mm_movemask_pd( v1 ) != 3)
		 scale = 0;
	     }	    	  
	 }
#else
	 v = &(x3[80 * i]);
	 scale = 1;
	 for(l = 0; scale && (l < 80); l++)
	   scale = ((ABS(v[l]) <  minlikelihood));
#endif

	 if (scale)
	   {
#ifdef __SIM_SSE3
	       __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	       
	       for(l = 0; l < 80; l+=2)
		 {
		   __m128d ex3v = _mm_load_pd(&v[l]);		  
		   _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
		 }		   		  
#else	     
	     for(l = 0; l < 80; l++)
	       v[l] *= twotothe256;
#endif

	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	  
	   }
       }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}

#endif

#ifdef _OPTIMIZED_FUNCTIONS

/*** BINARY DATA functions *****/

static void newviewGTRCAT_BINARY( int tipCase,  double *EV,  int *cptr,
                                  double *x1_start,  double *x2_start,  double *x3_start,  double *tipVector,
                                  int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                  size_t n,  double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double
    *le,
    *ri,
    *x1, *x2, *x3;
  size_t i, l;
  int
    scale, addScale = 0;

  switch(tipCase)
    {
    case TIP_TIP:
      {
        for(i = 0; i < n; i++)
          {
            x1 = &(tipVector[2 * tipX1[i]]);
            x2 = &(tipVector[2 * tipX2[i]]);
            x3 = &x3_start[2 * i];         

            le =  &left[cptr[i] * 4];
            ri =  &right[cptr[i] * 4];

            _mm_store_pd(x3, _mm_setzero_pd());     
                     
            for(l = 0; l < 2; l++)
              {                                                                                                                          
                __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&le[l * 2]));
                __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&ri[l * 2]));
                
                al = _mm_hadd_pd(al, al);
                ar = _mm_hadd_pd(ar, ar);
                
                al = _mm_mul_pd(al, ar);
                
                __m128d vv  = _mm_load_pd(x3);
                __m128d EVV = _mm_load_pd(&EV[2 * l]);
                
                vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
                
                _mm_store_pd(x3, vv);                                                     
              }            
          }
      }
      break;
    case TIP_INNER:
      {
        for (i = 0; i < n; i++)
          {
            x1 = &(tipVector[2 * tipX1[i]]);
            x2 = &x2_start[2 * i];
            x3 = &x3_start[2 * i];
            
            le =  &left[cptr[i] * 4];
            ri =  &right[cptr[i] * 4];

            _mm_store_pd(x3, _mm_setzero_pd());     
                     
            for(l = 0; l < 2; l++)
              {                                                                                                                          
                __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&le[l * 2]));
                __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&ri[l * 2]));
                
                al = _mm_hadd_pd(al, al);
                ar = _mm_hadd_pd(ar, ar);
                
                al = _mm_mul_pd(al, ar);
                
                __m128d vv  = _mm_load_pd(x3);
                __m128d EVV = _mm_load_pd(&EV[2 * l]);
                
                vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
                
                _mm_store_pd(x3, vv);                                                     
              }  
            
            __m128d minlikelihood_sse = _mm_set1_pd(minlikelihood);
         
            scale = 1;
            
            __m128d v1 = _mm_and_pd(_mm_load_pd(x3), absMask.m);
            v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
            if(_mm_movemask_pd( v1 ) != 3)
              scale = 0;                         
            
            if(scale)
              {
                __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
                
                __m128d ex3v = _mm_load_pd(x3);           
                _mm_store_pd(x3, _mm_mul_pd(ex3v,twoto));                                                 
                
                if(useFastScaling)
                  addScale += wgt[i];
                else
                  ex3[i]  += 1;   
              }                    
          }
      }
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
        {
          x1 = &x1_start[2 * i];
          x2 = &x2_start[2 * i];
          x3 = &x3_start[2 * i];

          le = &left[cptr[i] * 4];
          ri = &right[cptr[i] * 4];

          _mm_store_pd(x3, _mm_setzero_pd());       
          
          for(l = 0; l < 2; l++)
            {                                                                                                                            
              __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&le[l * 2]));
              __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&ri[l * 2]));
              
              al = _mm_hadd_pd(al, al);
              ar = _mm_hadd_pd(ar, ar);
              
              al = _mm_mul_pd(al, ar);
              
              __m128d vv  = _mm_load_pd(x3);
              __m128d EVV = _mm_load_pd(&EV[2 * l]);
              
              vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
              
              _mm_store_pd(x3, vv);                                                       
            }                             

          __m128d minlikelihood_sse = _mm_set1_pd(minlikelihood);
         
          scale = 1;
                  
          __m128d v1 = _mm_and_pd(_mm_load_pd(x3), absMask.m);
          v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
          if(_mm_movemask_pd( v1 ) != 3)
            scale = 0;                   
         
          if(scale)
            {
              __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
                    
              __m128d ex3v = _mm_load_pd(x3);             
              _mm_store_pd(x3, _mm_mul_pd(ex3v,twoto));                                           
             
              if(useFastScaling)
                addScale += wgt[i];
              else
                ex3[i]  += 1;     
           }             
        }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}

static void newviewGTRGAMMA_BINARY(int tipCase,
				   double *x1_start, double *x2_start, double *x3_start,
				   double *EV, double *tipVector,
				   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   const size_t n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling
				   )
{
  double
    *x1, *x2, *x3;
 
  size_t i, k, l;
  
  int
    scale, addScale = 0; 

  switch(tipCase)
    {
    case TIP_TIP:
      for (i = 0; i < n; i++)
       {
	 x1  = &(tipVector[2 * tipX1[i]]);
	 x2  = &(tipVector[2 * tipX2[i]]);
	 
	 for(k = 0; k < 4; k++)
	   {	     	     	    
	     x3 = &(x3_start[8 * i + 2 * k]);	     
	    	         
	     _mm_store_pd(x3, _mm_setzero_pd());	    
	    	     
	     for(l = 0; l < 2; l++)
	       {		 		 						   		  		 		 
		 __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&left[k * 4 + l * 2]));
		 __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&right[k * 4 + l * 2]));
		 		       
		 al = _mm_hadd_pd(al, al);
		 ar = _mm_hadd_pd(ar, ar);
		   
		 al = _mm_mul_pd(al, ar);
		   
		 __m128d vv  = _mm_load_pd(x3);
		 __m128d EVV = _mm_load_pd(&EV[2 * l]);
		 
		 vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
		 
		 _mm_store_pd(x3, vv);		     	  		   		  
	       }	     	    
	   }
       }
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
       {
	 x1  = &(tipVector[2 * tipX1[i]]);
	 
	 for(k = 0; k < 4; k++)
	   {	     	     
	     x2 = &(x2_start[8 * i + 2 * k]);
	     x3 = &(x3_start[8 * i + 2 * k]);	     
	    	         
	     _mm_store_pd(x3, _mm_setzero_pd());	    
	    	     
	     for(l = 0; l < 2; l++)
	       {		 		 						   		  		 		 
		 __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&left[k * 4 + l * 2]));
		 __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&right[k * 4 + l * 2]));
		 		       
		 al = _mm_hadd_pd(al, al);
		 ar = _mm_hadd_pd(ar, ar);
		   
		 al = _mm_mul_pd(al, ar);
		   
		 __m128d vv  = _mm_load_pd(x3);
		 __m128d EVV = _mm_load_pd(&EV[2 * l]);
		 
		 vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
		 
		 _mm_store_pd(x3, vv);		     	  		   		  
	       }	     	    
	   }
	
	 x3 = &(x3_start[8 * i]);
	 __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	 
	 scale = 1;
	 for(l = 0; scale && (l < 8); l += 2)
	   {
	     __m128d vv = _mm_load_pd(&x3[l]);
	     __m128d v1 = _mm_and_pd(vv, absMask.m);
	     v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	     if(_mm_movemask_pd( v1 ) != 3)
	       scale = 0;
	   }	    	         
	 
	 if(scale)
	   {
	     __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	     
	     for(l = 0; l < 8; l+=2)
	       {
		 __m128d ex3v = _mm_load_pd(&x3[l]);		  
		 _mm_store_pd(&x3[l], _mm_mul_pd(ex3v,twoto));	
	       }		   		  
	     
	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	  
	   }	 
       }      
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
       {	 
	 for(k = 0; k < 4; k++)
	   {	     
	     x1 = &(x1_start[8 * i + 2 * k]);
	     x2 = &(x2_start[8 * i + 2 * k]);
	     x3 = &(x3_start[8 * i + 2 * k]);	     
	    	         
	     _mm_store_pd(x3, _mm_setzero_pd());	    
	    	     
	     for(l = 0; l < 2; l++)
	       {		 		 						   		  		 		 
		 __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&left[k * 4 + l * 2]));
		 __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&right[k * 4 + l * 2]));
		 		       
		 al = _mm_hadd_pd(al, al);
		 ar = _mm_hadd_pd(ar, ar);
		   
		 al = _mm_mul_pd(al, ar);
		   
		 __m128d vv  = _mm_load_pd(x3);
		 __m128d EVV = _mm_load_pd(&EV[2 * l]);
		 
		 vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
		 
		 _mm_store_pd(x3, vv);		     	  		   		  
	       }	     	    
	   }
	
	 x3 = &(x3_start[8 * i]);
	 __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	 
	 scale = 1;
	 for(l = 0; scale && (l < 8); l += 2)
	   {
	     __m128d vv = _mm_load_pd(&x3[l]);
	     __m128d v1 = _mm_and_pd(vv, absMask.m);
	     v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	     if(_mm_movemask_pd( v1 ) != 3)
	       scale = 0;
	   }	    	         
	 
	 if(scale)
	   {
	     __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	     
	     for(l = 0; l < 8; l+=2)
	       {
		 __m128d ex3v = _mm_load_pd(&x3[l]);		  
		 _mm_store_pd(&x3[l], _mm_mul_pd(ex3v,twoto));	
	       }		   		  
	     
	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	  
	   }	 
       }
      break;

    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}


/**** BINARY DATA functions end ****/



#endif

#ifdef _OPTIMIZED_FUNCTIONS

/**** CLV computation at inner node of the tree *********************/


static void newviewGTRGAMMA_NSTATES(int tipCase,
				    double *x1, double *x2, double *x3, double *extEV, double *tipVector,
				    unsigned char *tipX1, unsigned char *tipX2,
				    size_t n, double *left, double *right, int *wgt, int *scalerIncrement, const size_t numberOfAllCharacters, const size_t numberOfStates, 
				    const size_t gammaRates)
{
  double  
    *uX1, 
    *uX2, 
    *v, 
    x1px2,
    *vl, 
    *vr;
  
  size_t
    i, 
    j, 
    l, 
    k;

  int
    addScale = 0;

  const size_t   
    loopLength = numberOfStates - (numberOfStates % VECTOR_WIDTH), //or 18 for testing!
    scalingLoopLength = loopLength * gammaRates,
    statesSquare = numberOfStates * numberOfStates,
    stride = numberOfStates * gammaRates,
    umpLength = numberOfAllCharacters * numberOfStates * gammaRates;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	double 
	  umpX1[umpLength], 
	  umpX2[umpLength];     

	//printf("TT\n");

	for(i = 0; i < numberOfAllCharacters; i++)
	  {
	    v = &(tipVector[numberOfStates * i]);	    

	    for(k = 0; k < stride; k++)
	      {
		double *ll =  &left[k * numberOfStates];
		double *rr =  &right[k * numberOfStates];
		
		VECTOR_REGISTER umpX1v = VECTOR_SET_ZERO();
		VECTOR_REGISTER umpX2v = VECTOR_SET_ZERO();	       			

		for(l = 0; l < loopLength; l += VECTOR_WIDTH)
		  {		   		    
		    VECTOR_REGISTER vv = VECTOR_LOAD(&v[l]);
		    		   
		    umpX1v = VECTOR_ADD(umpX1v, VECTOR_MUL(vv, VECTOR_LOAD(&ll[l])));
		    umpX2v = VECTOR_ADD(umpX2v, VECTOR_MUL(vv, VECTOR_LOAD(&rr[l])));					
		  }							
		
		umpX1[stride * i + k] = haddScalar(umpX1v);
		umpX2[stride * i + k] = haddScalar(umpX2v);
	
		for(;l < numberOfStates; l++)
		  {
		    umpX1[stride * i + k] += v[l] * ll[l]; 
		    umpX2[stride * i + k] += v[l] * rr[l];
		  }
		
	      }
	  }       

	for(i = 0; i < n; i++)
	  {
	    uX1 = &umpX1[stride * tipX1[i]];
	    uX2 = &umpX2[stride * tipX2[i]];

	    for(j = 0; j < gammaRates; j++)
	      {
		v = &x3[i * stride + j * numberOfStates];

		VECTOR_REGISTER zero =  VECTOR_SET_ZERO();
	       
		for(k = 0; k < loopLength; k += VECTOR_WIDTH)		  		    	      
		  VECTOR_STORE(&v[k], zero);
	       
		for(;k < numberOfStates; k++)
		  v[k] = 0.0;
	       
		for(k = 0; k < numberOfStates; k++)
		  { 
		    double 
		      *eev = &extEV[k * numberOfStates];
		    
		    x1px2 = uX1[j * numberOfStates + k] * uX2[j * numberOfStates + k];
		    
		    VECTOR_REGISTER 
		      x1px2v = VECTOR_SET_ONE(x1px2);

		    for(l = 0; l < loopLength; l += VECTOR_WIDTH)
		      {
		      	VECTOR_REGISTER vv = VECTOR_LOAD(&v[l]);
			VECTOR_REGISTER ee = VECTOR_LOAD(&eev[l]);

			vv = VECTOR_ADD(vv, VECTOR_MUL(x1px2v,ee));
			
			VECTOR_STORE(&v[l], vv);
		      }

		    for(;l < numberOfStates; l++)
		      v[l] += x1px2 * eev[l];
		  }
	      }	   
	  }
      }
      break;
    case TIP_INNER:
      {
	double 
	  umpX1[umpLength], 
	  ump_x2[numberOfStates];

	//printf("TI\n");

	for(i = 0; i < numberOfAllCharacters; i++)
	  {
	    v = &(tipVector[numberOfStates * i]);

	    for(k = 0; k < stride; k++)
	      {
		double *ll =  &left[k * numberOfStates];
				
		VECTOR_REGISTER umpX1v = VECTOR_SET_ZERO();
		
		for(l = 0; l < loopLength; l += VECTOR_WIDTH)
		  {
		    VECTOR_REGISTER vv = VECTOR_LOAD(&v[l]);
		    umpX1v = VECTOR_ADD(umpX1v, VECTOR_MUL(vv, VECTOR_LOAD(&ll[l])));		    					
		  }					

		umpX1[stride * i + k] = haddScalar(umpX1v);
		
		for(;l < numberOfStates; l++)		  
		  umpX1[stride * i + k] += v[l] * ll[l]; 	       
	      }
	  }

	for (i = 0; i < n; i++)
	  {
	    uX1 = &umpX1[stride * tipX1[i]];

	    for(k = 0; k < gammaRates; k++)
	      {
		v = &(x2[stride * i + k * numberOfStates]);
	       
		for(l = 0; l < numberOfStates; l++)
		  {		   
		    double *r =  &right[k * statesSquare + l * numberOfStates];
		    VECTOR_REGISTER ump_x2v = VECTOR_SET_ZERO();	    
		    
		    for(j = 0; j < loopLength; j+= VECTOR_WIDTH)
		      {
			VECTOR_REGISTER vv = VECTOR_LOAD(&v[j]);
			VECTOR_REGISTER rr = VECTOR_LOAD(&r[j]);
			ump_x2v = VECTOR_ADD(ump_x2v, VECTOR_MUL(vv, rr));
		      }
		     
		    ump_x2[l] = haddScalar(ump_x2v);

		    for(;j < numberOfStates; j++)
		      ump_x2[l] += v[j] * r[j];
		  }

		v = &(x3[stride * i + numberOfStates * k]);

		VECTOR_REGISTER zero =  VECTOR_SET_ZERO();
		
		for(l = 0; l < loopLength; l += VECTOR_WIDTH)		  		    
		  VECTOR_STORE(&v[l], zero);

		for(;l < numberOfStates; l++)
		  v[l] = 0.0;
		  
		for(l = 0; l < numberOfStates; l++)
		  {
		    double *eev = &extEV[l * numberOfStates];
		    x1px2 = uX1[k * numberOfStates + l]  * ump_x2[l];
		    VECTOR_REGISTER x1px2v = VECTOR_SET_ONE(x1px2);
		  
		    for(j = 0; j < loopLength; j += VECTOR_WIDTH)
		      {
			VECTOR_REGISTER vv = VECTOR_LOAD(&v[j]);
			VECTOR_REGISTER ee = VECTOR_LOAD(&eev[j]);
			
			vv = VECTOR_ADD(vv, VECTOR_MUL(x1px2v,ee));
			
			VECTOR_STORE(&v[j], vv);
		      }		     		    

		    for(;j < numberOfStates; j++)
		      v[j] += x1px2 * eev[j];
		  }			
	      }
	   	   
	    if(scaleEntry(stride, i, x3, scalingLoopLength))
	      addScale += wgt[i];		       	      	
	  }
      }
      break;
    case INNER_INNER:      	      
      for (i = 0; i < n; i++)
       {
	 for(k = 0; k < gammaRates; k++)
	   {
	     vl = &(x1[stride * i + numberOfStates * k]);
	     vr = &(x2[stride * i + numberOfStates * k]);
	     v =  &(x3[stride * i + numberOfStates * k]);

	     VECTOR_REGISTER zero =  VECTOR_SET_ZERO();
	     
	     for(l = 0; l < loopLength; l += VECTOR_WIDTH)		  		    
	       VECTOR_STORE(&v[l], zero);

	     for(;l < numberOfStates; l++)
	       v[l] = 0.0;

	     for(l = 0; l < numberOfStates; l++)
	       {		 		 
		 VECTOR_REGISTER al = VECTOR_SET_ZERO();
		 VECTOR_REGISTER ar = VECTOR_SET_ZERO();		  
		 
		 double *ll   = &left[k * statesSquare + l * numberOfStates];
		 double *rr   = &right[k * statesSquare + l * numberOfStates];
		 double *EVEV = &extEV[numberOfStates * l];
		 
		 double 
		   sal = 0.0,
		   sar = 0.0;
		 
		 for(j = 0; j < loopLength; j += VECTOR_WIDTH)
		   {
		     VECTOR_REGISTER lv  = VECTOR_LOAD(&ll[j]);
		     VECTOR_REGISTER rv  = VECTOR_LOAD(&rr[j]);
		     VECTOR_REGISTER vll = VECTOR_LOAD(&vl[j]);
		     VECTOR_REGISTER vrr = VECTOR_LOAD(&vr[j]);
		     
		     al = VECTOR_ADD(al, VECTOR_MUL(vll, lv));
		     ar = VECTOR_ADD(ar, VECTOR_MUL(vrr, rv));
		   }  		 
		       
		 //Hadd with broadcast!

		 //al = _mm_hadd_pd(al, al);
		 //ar = _mm_hadd_pd(ar, ar);

		 al = haddBroadCast(al);
		 ar = haddBroadCast(ar);
		 
		 
		 if(j < numberOfStates)
		   {
		     for(;j < numberOfStates; j++)
		       {
			 sal += (ll[j] * vl[j]);
			 sar += (rr[j] * vr[j]);
		       }
		     
		     al = VECTOR_ADD(al, VECTOR_SET_ONE(sal));
		     ar = VECTOR_ADD(ar, VECTOR_SET_ONE(sar));
		   }
		 
		 al = VECTOR_MUL(al, ar);
		 
		 for(j = 0; j < loopLength; j += VECTOR_WIDTH)
		   {
		     VECTOR_REGISTER vv  = VECTOR_LOAD(&v[j]);
		     VECTOR_REGISTER EVV = VECTOR_LOAD(&EVEV[j]);
		     
		     vv = VECTOR_ADD(vv, VECTOR_MUL(al, EVV));
		     
		     VECTOR_STORE(&v[j], vv);
		   }	
		 
		 if(j < numberOfStates)
		   {		       
		     VECTOR_STORE_LEFT(&sal, al);
		     for(;j < numberOfStates; j++)
		       v[j] += (sal * EVEV[j]);
		   }		   
	       }
		 
	   }	   

	 if(scaleEntry(stride, i, x3, scalingLoopLength))	   	    
	   addScale += wgt[i];		  
       }
      break;
    case TIP_TIP_CLV:

      //printf("TTC\n");      

      for (i = 0; i < n; i++)
	{
	  for(k = 0; k < gammaRates; k++)
	    {
	      vl = &(x1[numberOfStates * i]);
	      vr = &(x2[numberOfStates * i]);
	      v =  &(x3[stride * i + numberOfStates * k]);
	      
	      VECTOR_REGISTER zero =  VECTOR_SET_ZERO();
	      
	      for(l = 0; l < loopLength; l += VECTOR_WIDTH)		  		    
		VECTOR_STORE(&v[l], zero);
	      
	      for(;l < numberOfStates; l++)
		v[l] = 0.0;
	      
	      for(l = 0; l < numberOfStates; l++)
		{		 		 
		  VECTOR_REGISTER al = VECTOR_SET_ZERO();
		  VECTOR_REGISTER ar = VECTOR_SET_ZERO();		  
		  
		  double *ll   = &left[k * statesSquare + l * numberOfStates];
		  double *rr   = &right[k * statesSquare + l * numberOfStates];
		  double *EVEV = &extEV[numberOfStates * l];
		  
		  double 
		    sal = 0.0,
		    sar = 0.0;
		  
		  for(j = 0; j < loopLength; j += VECTOR_WIDTH)
		    {
		      VECTOR_REGISTER lv  = VECTOR_LOAD(&ll[j]);
		      VECTOR_REGISTER rv  = VECTOR_LOAD(&rr[j]);
		      VECTOR_REGISTER vll = VECTOR_LOAD(&vl[j]);
		      VECTOR_REGISTER vrr = VECTOR_LOAD(&vr[j]);
		      
		      al = VECTOR_ADD(al, VECTOR_MUL(vll, lv));
		      ar = VECTOR_ADD(ar, VECTOR_MUL(vrr, rv));		     
		    }  		 
		  
		  //Hadd with broadcast!
		  
		  //al = _mm_hadd_pd(al, al);
		  //ar = _mm_hadd_pd(ar, ar);
		  
		  al = haddBroadCast(al);
		  ar = haddBroadCast(ar);
		  
		  
		  if(j < numberOfStates)
		    {
		      for(;j < numberOfStates; j++)
			{
			  sal += (ll[j] * vl[j]);
			  sar += (rr[j] * vr[j]);
			}
		      
		      al = VECTOR_ADD(al, VECTOR_SET_ONE(sal));
		      ar = VECTOR_ADD(ar, VECTOR_SET_ONE(sar));
		    }
		  
		  al = VECTOR_MUL(al, ar);
		  
		  for(j = 0; j < loopLength; j += VECTOR_WIDTH)
		    {
		      VECTOR_REGISTER vv  = VECTOR_LOAD(&v[j]);
		      VECTOR_REGISTER EVV = VECTOR_LOAD(&EVEV[j]);
		      
		      vv = VECTOR_ADD(vv, VECTOR_MUL(al, EVV));
		      
		      VECTOR_STORE(&v[j], vv);
		    }	
		  
		  if(j < numberOfStates)
		    {		       
		      VECTOR_STORE_LEFT(&sal, al);
		      for(;j < numberOfStates; j++)
			v[j] += (sal * EVEV[j]);
		    }		   
		}
	    }	  	  
	  
	  //mth todo need scaling here?
	  if(scaleEntry(stride, i, x3, scalingLoopLength))
	    addScale += wgt[i];		  	  
	}         

       
  
      break;
    case TIP_INNER_CLV:
      
      //printf("TIC\n");
      
      for (i = 0; i < n; i++)
	{
	  for(k = 0; k < gammaRates; k++)
	    {
	      vl = &(x1[numberOfStates * i]);
	      vr = &(x2[stride * i + numberOfStates * k]);
	      v =  &(x3[stride * i + numberOfStates * k]);

	     VECTOR_REGISTER zero =  VECTOR_SET_ZERO();
	     
	     for(l = 0; l < loopLength; l += VECTOR_WIDTH)		  		    
	       VECTOR_STORE(&v[l], zero);

	     for(;l < numberOfStates; l++)
	       v[l] = 0.0;

	     for(l = 0; l < numberOfStates; l++)
	       {		 		 
		 VECTOR_REGISTER al = VECTOR_SET_ZERO();
		 VECTOR_REGISTER ar = VECTOR_SET_ZERO();		  
		 
		 double *ll   = &left[k * statesSquare + l * numberOfStates];
		 double *rr   = &right[k * statesSquare + l * numberOfStates];
		 double *EVEV = &extEV[numberOfStates * l];
		 
		 double 
		   sal = 0.0,
		   sar = 0.0;
		 
		 for(j = 0; j < loopLength; j += VECTOR_WIDTH)
		   {
		     VECTOR_REGISTER lv  = VECTOR_LOAD(&ll[j]);
		     VECTOR_REGISTER rv  = VECTOR_LOAD(&rr[j]);
		     VECTOR_REGISTER vll = VECTOR_LOAD(&vl[j]);
		     VECTOR_REGISTER vrr = VECTOR_LOAD(&vr[j]);
		     
		     al = VECTOR_ADD(al, VECTOR_MUL(vll, lv));
		     ar = VECTOR_ADD(ar, VECTOR_MUL(vrr, rv));
		   }  		 
		       
		 //Hadd with broadcast!

		 //al = _mm_hadd_pd(al, al);
		 //ar = _mm_hadd_pd(ar, ar);

		 al = haddBroadCast(al);
		 ar = haddBroadCast(ar);
		 		 
		 if(j < numberOfStates)
		   {
		     for(;j < numberOfStates; j++)
		       {
			 sal += (ll[j] * vl[j]);
			 sar += (rr[j] * vr[j]);
		       }
		     
		     al = VECTOR_ADD(al, VECTOR_SET_ONE(sal));
		     ar = VECTOR_ADD(ar, VECTOR_SET_ONE(sar));
		   }
		 
		 al = VECTOR_MUL(al, ar);
		 
		 for(j = 0; j < loopLength; j += VECTOR_WIDTH)
		   {
		     VECTOR_REGISTER vv  = VECTOR_LOAD(&v[j]);
		     VECTOR_REGISTER EVV = VECTOR_LOAD(&EVEV[j]);
		     
		     vv = VECTOR_ADD(vv, VECTOR_MUL(al, EVV));
		     
		     VECTOR_STORE(&v[j], vv);
		   }	
		 
		 if(j < numberOfStates)
		   {		       
		     VECTOR_STORE_LEFT(&sal, al);
		     for(;j < numberOfStates; j++)
		       v[j] += (sal * EVEV[j]);
		   }		   
	       }
		 
	    }	   

	  if(scaleEntry(stride, i, x3, scalingLoopLength))
	    addScale += wgt[i];	
	}
      break;
    default:
      assert(0);
    }

  
  *scalerIncrement = addScale;
}

#endif
