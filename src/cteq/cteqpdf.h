/**
 * @file lhpdf.h
 */
#ifndef __cteqpdf_h__
#define __cteqpdf_h__ 1


#ifdef __cplusplus
extern "C" {
#endif

    
  /** There are several available pdf sets. We have to build a database and 
   this class helps to do this. */ 
  typedef struct {
	/*  The ID of the pdf set. This is the same what is used by CTEQ. */
	const int id;
	
	/*  Name of the pdf set. This is the same what is used by CTEQ. E.g.: CTEQ6M */
	const char *name;
	
	/*  Here we can provide some breif description of the pdf set.  */
	const char *description;
	
	/*  This is the path of the table file. It MUST BE absolute path 
	 otherwise it is considered as relative to the working directory. */
	const char *path;

	/*  Name of the pdf table file */
	const char *filename;
	
	/*   Type of the pdf table: 1 ==> .tbl and 2 ==> .pds  */
	const unsigned int itbl;
  } cteq_pdfset_t;
  
  
  /**   The list of the available pdf sets */
  extern const cteq_pdfset_t *cteq_pdfset_database;
 
  /**  Get the database entry using the name of the pdfset. 
   It returns the ponter to the coorespond entry, other it returs null pointer. */ 
  const cteq_pdfset_t * cteq_pdfset_find(const cteq_pdfset_t *, const char *);
   
  /**  Get the database entry using the CTEQ ID of the pdfset. 
   It returns the ponter to the coorespond entry, other it returs null pointer. */ 
  const cteq_pdfset_t * cteq_pdfset_find_id(const cteq_pdfset_t *, int);
  
  
  /**  The structure to store the tabele and related parameters. */ 
  typedef struct 
  {
	/*  alphas parameters  */
	unsigned int order, nf;
	
	/* Lambda(nf=5,MSbar), quark masses */
	double lambda[7], mass[7];

	/* max Nf MxVal  */
	unsigned int nfmx, mxval;
    
	/*  Qini, Qmax, xmin  */
	double qini, qmax, xmin;
	
	/* xv, tv arrays and the pdf grid */
	unsigned int nx, nt;
	double *xv, *xvpow, *tv, *upd;
  } cteq_pdf_t;
  
  
  /**  Creates a pdf set using the database entry */
  cteq_pdf_t * cteq_pdf_alloc(const cteq_pdfset_t *);
  
  /**  Creates a pdf set using its name. It will look up it 
   in the database if it doesn't find the wanted pdf then it try to 
   open it as a table file. */
  cteq_pdf_t * cteq_pdf_alloc_name(const char *);
 
  /**  Creates a pdf using the CTEQ ID of the pdfset. */
  cteq_pdf_t * cteq_pdf_alloc_id(int);
  
  /**  Destructor */
  void cteq_pdf_free(cteq_pdf_t *);

  /**  Evolve the parton distribution function  */
  double cteq_pdf_evolvepdf(const cteq_pdf_t *, int, double, double);
  
  /**  Evolve the \f$\alpha_s\f$  */
  double cteq_pdf_evolveas(const cteq_pdf_t *pdf, double);
  
  /**   Returns the evolution order of the pdf set  */
  unsigned int cteq_pdf_orderpdf(const cteq_pdf_t *);
  
  /**   Returns the evolution order of the \f$\alpha_s\f$ */
  unsigned int cteq_pdf_orderas(const cteq_pdf_t *);
  
  /**   Returns the number of the active flavours  */
  unsigned int cteq_pdf_nfmax(const cteq_pdf_t *);

  /**   Returns the fitting scale  */
  double cteq_pdf_scale(const cteq_pdf_t *);
   
  /**   Returns the alphas at the fitting scale  */
  double cteq_pdf_alfas(const cteq_pdf_t *);
  
  /**   Returns the quark masses  */
  double cteq_pdf_mass(const cteq_pdf_t *, int);
  
  /**   Returns the flavour threshold in the evolution */
  double cteq_pdf_threshold(const cteq_pdf_t *, unsigned int);
  
  
  
  /*****
   *****  FIX THIS PART
   *****/
  
  /**   Returns the description  */
  const char ** cteq_pdf_desc(const cteq_pdf_t *);

   
  

#ifdef __cplusplus
}
#endif


#endif
