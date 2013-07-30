/*
 *  cteqpdf-f77.c
 *  cteq-pdf
 *
 *  Created by Zoltan Nagy on 3/24/08.
 *  Copyright 2008 Zoltan Nagy. All rights reserved.
 *
 */

/* Standard C includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* lhpdf includes */
#include "cteqpdf.h"



/*   kind of virtual memory to create pdf object from f77 code  */
static size_t __f77_cteq_pdf_array_size = 0;
static cteq_pdf_t **__f77_cteq_pdf_array = 0;


static
void __f77_cteq_pdf_array_alloc(size_t n)
{
  if(!(__f77_cteq_pdf_array = realloc(__f77_cteq_pdf_array, n*sizeof(cteq_pdf_t *)))) {
   fprintf(stderr, "F77 CTEQ PDF: Unable to allocate memory for the array of pdf objects.");
   exit(1);
  }
  
  /*   set all the new elements to 0  */
  size_t i;
  for(i = __f77_cteq_pdf_array_size; i < n; ++i) 
	__f77_cteq_pdf_array[i] = 0;
  
  __f77_cteq_pdf_array_size = n;
}


int cteq_pdf_alloc_id_(int *id) 
{
  /*   Find an available pointer. */
  size_t i = 0;
  
  if(__f77_cteq_pdf_array)
	while(__f77_cteq_pdf_array[i] && (i < __f77_cteq_pdf_array_size)) ++i;

  /*   Increases the size of the array by 50  */
  if(i == __f77_cteq_pdf_array_size)
	__f77_cteq_pdf_array_alloc(__f77_cteq_pdf_array_size + 50);
  
  /*   Create the pdf   */
  __f77_cteq_pdf_array[i] = cteq_pdf_alloc_id(*id);
  
  return (__f77_cteq_pdf_array[i] ? (int) i : -1);
}

int cteq_pdf_alloc_name_(char *name, long int len) 
{
  /*   Find an available pointer. */
  size_t i = 0;
  
  if(__f77_cteq_pdf_array)
	while(__f77_cteq_pdf_array[i] && (i < __f77_cteq_pdf_array_size)) ++i;
  
  /*   Increases the size of the array by 50  */
  if(i == __f77_cteq_pdf_array_size)
	__f77_cteq_pdf_array_alloc(__f77_cteq_pdf_array_size + 50);
  
  /*   Create the pdf   */
  char *__name = malloc((len+1)*sizeof(char));
  
  if(!__name) return -1;
  strncpy(__name, name, len); __name[len]='\0';
  
  __f77_cteq_pdf_array[i] = cteq_pdf_alloc_name(__name);
  free(__name);
  
  return (__f77_cteq_pdf_array[i] ? (int) i : -1);
}

void cteq_pdf_free_(int *id) 
{
  if(0 < *id && *id < __f77_cteq_pdf_array_size)
    if(__f77_cteq_pdf_array[*id]) {
      cteq_pdf_free(__f77_cteq_pdf_array[*id]);
      __f77_cteq_pdf_array[*id] = 0;
    }
  
  /*  Set the fortran id to -1.   */
  *id = -1;
}

double cteq_pdf_evolvepdf_(int *pdf, int *iprtn, double *x, double *q)
{
  return cteq_pdf_evolvepdf(__f77_cteq_pdf_array[*pdf], *iprtn, *x, *q);
}

double cteq_pdf_evolveas_(int *pdf, double *q)
{
  return cteq_pdf_evolveas(__f77_cteq_pdf_array[*pdf], *q);
}


/**   Returns the evolution order of the pdf set  */
unsigned int cteq_pdf_orderpdf_(int *pdf) {
  return __f77_cteq_pdf_array[*pdf]->order;
}

/**   Returns the evolution order of the \f$\alpha_s\f$ */
unsigned int cteq_pdf_orderas_(int *pdf) {
  return __f77_cteq_pdf_array[*pdf]->order;
}

/**   Returns the number of the active flavours  */
unsigned int cteq_pdf_nfmax_(int *pdf) {
  return __f77_cteq_pdf_array[*pdf]->nfmx;
}

/**   Returns the fitting scale  */
double cteq_pdf_scale_(int *pdf) {
  return __f77_cteq_pdf_array[*pdf]->qini;
}

/**   Returns the alphas at the fitting scale  */
double cteq_pdf_alfas_(int *pdf) {
  return cteq_pdf_alfas(__f77_cteq_pdf_array[*pdf]);
}

/**   Returns the masses  */
double cteq_pdf_mass_(int *pdf, int *i) {
  return __f77_cteq_pdf_array[*pdf]->mass[abs(*i)];
}

/**   Returns the flavour threshold in the evolution */
double cteq_pdf_threshold_(int *pdf, int *nf) {
  return (*nf < 4 ? 0.0 : __f77_cteq_pdf_array[*pdf]->mass[*nf]);
}





