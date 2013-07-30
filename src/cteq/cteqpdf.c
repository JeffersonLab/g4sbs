/* Standard C includes */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* lhpdf includes */
#include "cteqpdf.h"

/*    some technical constatnts  */
static const double xpow = 0.3;
static const double onep = 1.00001;
static const unsigned int maxval = 4;


static
double __cteq_pdf_as(unsigned int ord, unsigned int nf, double q, double lmd)
{
  double b0 = 11.0 - 2.0/3.0*nf;
  double t = log(q/lmd);
  double as = 1.0/(b0*t);
  
  /*  it is a leading order evolution  */ 
  if(ord <= 1) return as;
  
  /*  at NLO or higer order level returns with the NLO alpha_s  */
  double b1 = 51.0 - 19.0/3.0*nf;
  return as*(1.0-b1/b0*as*log(2.0*t));
}

static 
double __cteq_pdf_astolmd(unsigned int ord, unsigned int nf, double as, double q)
{
  double b0 = 11.0 - 2.0/3.0*nf;
  double t = 1.0/(b0*as);
  
  /*  it is a leading order evolution  */ 
  if(ord <= 1) return q*exp(-t);
  
  /*  at NLO or higer order level returns with the NLO alpha_s  */
  double as0, as1, ot, lt, br = (51.0 - 19.0/3.0*nf)/(b0*b0);
  
  do {
	lt = log(2.0*t)/t;
    ot = t;
	
    as0 = (1.0 - br*lt)/(b0*t);
    as1 = (-1.0 - br*(1.0/t-2.0*lt))/(b0*t*t);
    t += (as - as0)/as1;
  } while(fabs(ot-t)/ot > 1e-5);
  
  return q*exp(-t);
}

static
void __cteq_pdf_setlmd(cteq_pdf_t *pdf)
{
  double as;
  unsigned int nf;
  
  for(nf = pdf->nf+1; nf <= 6; ++nf) {
	as = __cteq_pdf_as(pdf->order, nf-1, pdf->mass[nf], pdf->lambda[nf-1]);
	pdf->lambda[nf] = __cteq_pdf_astolmd(pdf->order, nf, as, pdf->mass[nf]);
  }
  
  /*   Under the charm mass every quark is considered as massless.  */
  for(nf = pdf->nf-1; nf > 2; --nf) {
	as = __cteq_pdf_as(pdf->order, nf+1, pdf->mass[nf+1], pdf->lambda[nf+1]);
	pdf->lambda[nf] = __cteq_pdf_astolmd(pdf->order, nf, as, pdf->mass[nf+1]);
  }
}



static
double __cteq_pdf_polint4f(double *xa, double *ya, double x)
{
  double c1, d1, d2, c2, d3, h1, h2, h3, h4, c3, cc1, cd1, 
  cd2, cc2, dd1, dc1, den;
  
  /* Function Body */
  h1 = xa[0] - x; h2 = xa[1] - x;
  h3 = xa[2] - x; h4 = xa[3] - x;
  
  den = (ya[1] - ya[0])/(h1 - h2);
  d1 = h2*den; c1 = h1*den;
  
  den = (ya[2] - ya[1])/(h2 - h3);
  d2 = h3*den; c2 = h2*den;
  
  den = (ya[3] - ya[2])/(h3 - h4);
  d3 = h4*den; c3 = h3*den;
  
  den = (c2 - d1)/(h1 - h3);
  cd1 = h3*den; cc1 = h1*den;
  
  den = (c3 - d2)/(h2 - h4);
  cd2 = h4*den; cc2 = h2*den;
  
  den = (cc2 - cd1)/(h1 - h4);
  dd1 = h4*den; dc1 = h1*den;
  
  if(h3 + h4 < 0.0) return ya[3] + d3 + cd2 + dd1;
  if(h2 + h3 < 0.0) return ya[2] + d2 + cd1 + dc1;
  if(h1 + h2 < 0.0) return ya[1] + c2 + cd1 + dc1;
  
  return ya[0] + c1 + cc1 + dc1;
} 


static 
double __cteq_pdf_pardis(const cteq_pdf_t *pdf, int iprtn, double x, double q)
{
  int jm, jx, jlx = -1, ju, jq, jlq = -1, j1, ip, jtmp, it;
  unsigned int nx = pdf->nx, nq = pdf->nt;
  double tt, ss, fvec[4];
  
  if(iprtn != 0)
    if(q <= pdf->mass[abs(iprtn)]) 
      return 0.0;
  
  tt = log(log(q/pdf->lambda[pdf->nf]));
  ss = pow(x, xpow);
  
  ju = nx + 1;
  while(ju - jlx > 1) {
    jm = (ju + jlx)/2;
    if(x >= pdf->xv[jm]) jlx = jm;
    else ju = jm;
  }
  
  if(jlx <= -1) {
    fprintf(stderr, "cteq-pdf: Severe error x <= 0! x = %g \n", x);
    exit(-1);
  } else if(jlx == 0) jx = 0;
  else if(jlx <= (int) nx-2) jx = jlx-1;
  else if(jlx == (int) nx-1 || x < onep) jx = jlx-2;
  else {
    fprintf(stderr, "cteq-pdf: Severe error: x > 1!  x = %g", x);
    exit(-1);
  }
  
  ju = nq + 1;
  while(ju - jlq > 1) {
    jm = (ju + jlq)/2;
    if (tt >= pdf->tv[jm]) jlq = jm;
    else ju = jm;
  }
  
  if(jlq <= 0) jq = 0;
  else if(jlq <= (int) nq-2) jq = jlq-1;
  else jq = nq-3;
  
  /* get the pdf function values at the lattice points... */
  ip = (iprtn > (int) pdf->mxval ? -iprtn : iprtn);
  jtmp = ((ip + pdf->nfmx)*(nq+1) + (jq-1))*(nx+1) + jx + 1;
  
  if(jx == 0) 
  {
	double fij[4];	
	for(it = 0; it < 4; ++it) {
	  j1 = jtmp + (it+1)*(nx+1);
	  
	  fij[0] = 0.0;
	  fij[1] = (pdf->upd[j1  ])*(pdf->xv[1])*(pdf->xv[1]);
	  fij[2] = (pdf->upd[j1+1])*(pdf->xv[2])*(pdf->xv[2]);
	  fij[3] = (pdf->upd[j1+2])*(pdf->xv[3])*(pdf->xv[3]);
	  
	  /*  Use Polint which allows x to be anywhere w.r.t. the grid */
	  fvec[it] = __cteq_pdf_polint4f(pdf->xvpow, fij, ss)/(x*x);
	}	
  } 
  else if(jlx == nx-1) 
  {
	for(it = 0; it < 4; ++it)
	  fvec[it] = __cteq_pdf_polint4f(pdf->xvpow+nx-3, pdf->upd+jtmp+(it+1)*(nx+1)-1, ss);
  } 
  else 
  {
	double *svec = pdf->xvpow + jx-1;
	double s12 = svec[1]-svec[2], s13 = svec[1]-svec[3], s23 = svec[2]-svec[3],
	s24 = svec[2]-svec[4], s34 = svec[3]-svec[4], sy2 = ss-svec[2], sy3 = ss-svec[3];
	
	double const1 = s13/s23, const2 = s12/s23, const3 = s34/s23, const4 = s24/s23;
	double s1213 = s12 + s13, s2434 = s24 + s34;
	double sdet = s12*s34 - s1213*s2434;
	double tmp = sy2*sy3/sdet; 
	double const5 = (s34*sy2 - s2434*sy3)*tmp/s12, const6 = (s1213*sy2 - s12*sy3)*tmp/s34;
	
	for(it = 0; it < 4; ++it) {
	  j1 = jtmp + (it+1)*(nx+1);
	  
	  double sf2 = pdf->upd[j1], sf3 = pdf->upd[j1+1];
	  double g1 = sf2*const1 - sf3*const2, g4 = sf3*const4 - sf2*const3;
	  fvec[it] = (const5*(pdf->upd[j1-1]-g1) + const6*(pdf->upd[j1+2]-g4) + sf2*sy3-sf3*sy2)/s23;
	}
  }
  
  /*   interpolate in t... */
  if(jlq <= 0) return __cteq_pdf_polint4f(pdf->tv, fvec, tt);
  if(jlq >= nq-1) return __cteq_pdf_polint4f(pdf->tv+nq-3, fvec, tt);
  
  double *tvec = pdf->tv + jq-1;
  double t12 = tvec[1]-tvec[2], t13 = tvec[1]-tvec[3], t23 = tvec[2]-tvec[3],
  t24 = tvec[2]-tvec[4], t34 = tvec[3]-tvec[4], ty2 = tt-tvec[2], ty3 = tt-tvec[3];
  
  double tmp1 = t12 + t13, tmp2 = t24 + t34;
  double tdet = t12*t34 - tmp1*tmp2;

  double g1 = (fvec[1]*t13 - fvec[2]*t12)/t23, g4 = (fvec[2]*t24 - fvec[1]*t34)/t23;
  double h00 = (t34*ty2 - tmp2*ty3)*(fvec[0]-g1)/t12 + (tmp1*ty2 - t12*ty3)*(fvec[3]-g4)/t34;
  
  return (h00*ty2*ty3/tdet + fvec[1]*ty3 - fvec[2]*ty2)/t23;
}

  
static 
cteq_pdf_t * __cteq_pdf_alloc_read_tbl(FILE *file)
{
  unsigned int i;
  
  /*   Allocate the pdf set   */
  cteq_pdf_t *pdf = malloc(sizeof(cteq_pdf_t));
  if(!pdf) return 0;
  
  /*   Reading the first two line    */
  int del;
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');
  
  /*   Reading ord, nfl  */
  float __ord, __nf; 
  fscanf(file, "%f%f", &__ord, &__nf);
  
  pdf->order = (unsigned int) __ord;
  pdf->nf = (unsigned int) __nf;
  pdf->mxval = 2;

  /*   Reading ord, nfl, al mass[1-6]   */
  pdf->mass[0] = 0.0;     /*   gluon mass  */
  fscanf(file, "%lf%lf%lf%lf%lf%lf%lf", pdf->lambda + pdf->nf, 
		 pdf->mass+1,pdf->mass+2,pdf->mass+3,pdf->mass+4,pdf->mass+5,pdf->mass+6);

  /*   Calculating the lambda values at the thresholds  */
  __cteq_pdf_setlmd(pdf);
  
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');

  /*   Reading nx, nt, nfmx   */
  fscanf(file, "%u%u%u", &(pdf->nx), &(pdf->nt), &(pdf->nfmx));
  
  /*   Allocating memory for the grid   */
  unsigned int nupd = (pdf->nx+1)*(pdf->nt+1)*(pdf->nfmx+1+pdf->mxval);

  if(!(pdf->xv = (double *) malloc((pdf->nx + 1)*sizeof(double)))) goto label_free_pdf;
  if(!(pdf->xvpow = (double *) malloc((pdf->nx +1)*sizeof(double)))) goto label_free_xv;
  if(!(pdf->tv = (double *) malloc((pdf->nt + 1)*sizeof(double)))) goto label_free_xvpow;
  if(!(pdf->upd = (double *) malloc(nupd*sizeof(double)))) goto label_free_tv;
  
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');

  /*   Reading qini nad qmax   */
  fscanf(file, "%le%le", &(pdf->qini), &(pdf->qmax));
  
  /*   Reading q values and converting to t  */
  double q;
  for(i = 0; i <= pdf->nt; i++) {
	fscanf(file, "%le", &q);
	pdf->tv[i] = log(log(q/(pdf->lambda[pdf->nf])));
  }

  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');
  
  /*   Reading xmin   */
  fscanf(file, "%le", &(pdf->xmin));

  /*   Reading x values and calculating xvpow values  */
  double x;
  for(i = 0; i <= pdf->nx; i++) {
	fscanf(file, "%le", &x);
	pdf->xv[i] = x;
	pdf->xvpow[i] = pow(x, xpow);
  }
  
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');

  /*   Reading the grid   */
  for(i = 0; i < nupd; i++)
	fscanf(file, "%le", pdf->upd+i);
  
  return pdf;

  /*  garbage collection   */
label_free_tv: free(pdf->tv);
label_free_xvpow: free(pdf->xvpow);
label_free_xv: free(pdf->xv);
label_free_pdf: free(pdf);
  
  return 0;
}

static 
cteq_pdf_t * __cteq_pdf_alloc_read_pds(FILE *file)
{
  unsigned int i;
  
  /*   Allocate the pdf set   */
  cteq_pdf_t *pdf = malloc(sizeof(cteq_pdf_t));
  if(!pdf) return 0;
  
  /*   Reading the first two line    */
  int del;
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');
  
  /*   Reading ord, nfl   */
  float __ord, __nf;
  fscanf(file, "%f%f", &__ord, &__nf);
  
  pdf->order = (unsigned int) __ord;
  pdf->nf = (unsigned int) __nf;
  
  /*   Reading ord, nfl, al mass[1-6]   */
  pdf->mass[0] = 0.0;     /*   gluon mass  */
  fscanf(file, "%lf%lf%lf%lf%lf%lf%lf", pdf->lambda + pdf->nf, 
		 pdf->mass+1,pdf->mass+2,pdf->mass+3,pdf->mass+4,pdf->mass+5,pdf->mass+6);
  
  /*   Calculating the lambda values at the thresholds  */
  __cteq_pdf_setlmd(pdf);
  
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');
  
  /*   Reading nfmx and mxval   */
  fscanf(file, "%d%d%d%u%u", &del, &del, &del, &(pdf->nfmx), &(pdf->mxval));
  if(pdf->mxval > maxval) pdf->mxval = 3;
  
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');

  /*   Reading nx, n and ng   */
  unsigned int ng;
  fscanf(file, "%u%u%d%u", &(pdf->nx), &(pdf->nt), &del, &ng);
 
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');
  
  /*   Skipping the next ng lines and the comment starts QINI, QMAX,...,   */
  for(i = 0; i <= ng; i++) 
	do del = fgetc(file); while(del != '\n');
  
  /*   Allocating memory for the grid   */
  unsigned int nupd = (pdf->nx+1)*(pdf->nt+1)*(pdf->nfmx+1+pdf->mxval);

  if(!(pdf->xv = (double *) malloc((pdf->nx + 1)*sizeof(double)))) goto label_free_pdf;
  if(!(pdf->xvpow = (double *) malloc((pdf->nx +1)*sizeof(double)))) goto label_free_xv;
  if(!(pdf->tv = (double *) malloc((pdf->nt + 1)*sizeof(double)))) goto label_free_xvpow;
  if(!(pdf->upd = (double *) malloc(nupd*sizeof(double)))) goto label_free_tv;
    
  /*   Reading qini nad qmax   */
  fscanf(file, "%le%le", &(pdf->qini), &(pdf->qmax));
  
  /*   Reading t values */
  double tmp;
  for(i = 0; i <= pdf->nt; i++) 
	fscanf(file, "%le%le", &tmp, &(pdf->tv[i])); 
  
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');
  
  /*   Reading xmin   */
  fscanf(file, "%le%le", &(pdf->xmin), &tmp);
  
  /*   Reading x values and calculating xvpow values  */
  pdf->xv[i] = pdf->xvpow[i] = 0.0;
  
  for(i = 1; i <= pdf->nx; i++) {
	fscanf(file, "%le", &(pdf->xv[i]));
	pdf->xvpow[i] = pow(pdf->xv[i], xpow);
  }
  
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');
  
  /*   Reading the grid   */
  for(i = 0; i < nupd; i++)
	fscanf(file, "%le", pdf->upd+i);
  
  return pdf;
  
  /*  garbage collection   */
label_free_tv: free(pdf->tv);
label_free_xvpow: free(pdf->xvpow);
label_free_xv: free(pdf->xv);
label_free_pdf: free(pdf);
  
  return 0;
}


/*    Exported functions    */
cteq_pdf_t * cteq_pdf_alloc(const cteq_pdfset_t *pdfset)
{
  /*   return value  */
  cteq_pdf_t *pdf = 0;
  
  /*   we need the full filename with the absolute path */
  size_t len = strlen(pdfset->path) + strlen(pdfset->filename);
  char *filename = malloc((len+1)*sizeof(char));
  
  if(!filename) {
	fprintf(stderr, "cteq_pdf: Unable to determine the full filename of the table file!\n");
	fprintf(stderr, "    path     : %s\n", pdfset->path);
	fprintf(stderr, "    filename : %s\n", pdfset->filename);
	return pdf;
  }

  strcpy(filename, pdfset->path);
  strcat(filename, pdfset->filename);
  
  /*   Openning the table file  */
  FILE *file = fopen(filename, "r");
  
  if(!file) {
	fprintf(stderr, "cteq_pdf: Unable to open file %s\n", filename);
	free(filename);
	return pdf;
  }
  
  /*    Creating the pdf set   */
  switch(pdfset->itbl) {
	case 1: pdf = __cteq_pdf_alloc_read_tbl(file); break;
	case 2: pdf = __cteq_pdf_alloc_read_pds(file); break;
	default: pdf = 0;
  }
  
  /*    Closing the table file  */
  fclose(file);
  free(filename);
  
  return pdf;
}


cteq_pdf_t * cteq_pdf_alloc_name(const char *name)
{
  /*   First we try to find the pdfset in the database. */
  const cteq_pdfset_t *pdfset = cteq_pdfset_find(cteq_pdfset_database, name);
  if(pdfset) return cteq_pdf_alloc(pdfset);
  
  /*   Otherwise we try to open it as a table file. */
  /*   return value  */
  cteq_pdf_t *pdf = 0;
  
  /*   The table files must have .tbl or .pds extension. 
   Unfortunately there is no infor about the type of the table in the file,
   so we have to use the extension.   */
  unsigned int itbl = 0;
  const char *ext = name;
  while(*ext != '\0') ++ext; 
  
  if(strcmp(ext-4, ".tbl") == 0) itbl = 1;
  else if(strcmp(ext-4, ".pds") == 0) itbl = 2;
  else {
	fprintf(stderr, "cteq_pdf: Unable to identify the type of the table. Please use extension .tbl or .pds!\n");
	return pdf;
  }
  
  /*   Openning the table file  */
  FILE *file = fopen(name, "r");
  
  if(!file) {
	fprintf(stderr, "cteq_pdf: Unable to open file %s\n", name);
	return pdf;
  }
  
  /*    Creating the pdf set   */
  switch(itbl) {
	case 1: pdf = __cteq_pdf_alloc_read_tbl(file); break;
	case 2: pdf = __cteq_pdf_alloc_read_pds(file); break;
	default: pdf = 0;
  }
  
  /*    Closing the table file  */
  fclose(file);
  
  return pdf;
}


cteq_pdf_t * cteq_pdf_alloc_id(int id)
{
  /*   First we try to find the pdfset in the database. */
  const cteq_pdfset_t *pdfset = cteq_pdfset_find_id(cteq_pdfset_database, id);
  if(pdfset) return cteq_pdf_alloc(pdfset);
  
  return 0;
}


void cteq_pdf_free(cteq_pdf_t *pdf) 
{
  if(pdf == 0) return;
  
  if(pdf->xv) free(pdf->xv);
  if(pdf->tv) free(pdf->tv);
  if(pdf->xvpow) free(pdf->xvpow);
  if(pdf->upd) free(pdf->upd);
  
  free(pdf);
}

double cteq_pdf_evolvepdf(const cteq_pdf_t *pdf, int iprtn, double x, double q)
{
  if(abs(iprtn) > pdf->nfmx) return 0.0;  
  double ff = __cteq_pdf_pardis(pdf, iprtn, x, q);
  return ff <= 0.0 ? 0.0 : ff;
}


double cteq_pdf_evolveas(const cteq_pdf_t *pdf, double q)
{
  unsigned int nf = 6;
  while(q < pdf->mass[nf] && nf > 3) --nf;
  return __cteq_pdf_as(pdf->order, nf, q, pdf->lambda[nf]);
}

/**   Returns the evolution order of the pdf set  */
unsigned int cteq_pdf_orderpdf(const cteq_pdf_t *pdf) {
  return pdf->order;
}

/**   Returns the evolution order of the \f$\alpha_s\f$ */
unsigned int cteq_pdf_orderas(const cteq_pdf_t *pdf) {
  return pdf->order;
}

/**   Returns the number of the active flavours  */
unsigned int cteq_pdf_nfmax(const cteq_pdf_t *pdf) {
  return pdf->nfmx;
}

/**   Returns the fitting scale  */
double cteq_pdf_scale(const cteq_pdf_t *pdf) {
  return pdf->qini;
}

/**   Returns the alphas at the fitting scale  */
double cteq_pdf_alfas(const cteq_pdf_t *pdf) {
  return cteq_pdf_evolveas(pdf, pdf->qini);
}

/**   Returns the masses  */
double cteq_pdf_mass(const cteq_pdf_t *pdf, int i) {
  return pdf->mass[abs(i)];
}

/**   Returns the flavour threshold in the evolution */
double cteq_pdf_threshold(const cteq_pdf_t *pdf, unsigned int nf) {
  return (nf < 4 ? 0.0 : pdf->mass[nf]);
}


