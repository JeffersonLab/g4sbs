/*
 *  cteq-database.c
 *  
 *
 *  Created by Zoltan Nagy on 3/22/08.
 *  Copyright 2008 Zoltan Nagy. All rights reserved.
 *
 *  Modified by Seamus Riordan
 *
 */

/* Standard C includes */
#include <string.h>


/* cteq-pdf includes */
#include "cteqpdf.h"

//Add by Jixie: in case the makefile does not define CTEQ_TBL_PATH
//use this definition
#ifndef CTEQ_TBL_PATH
#define  CTEQ_TBL_PATH "cteq-tbls"
#endif

#define CTEQ6STD_TBL_PATH    CTEQ_TBL_PATH"/cteq6std/"
#define CTEQ6_TBL_PATH       CTEQ_TBL_PATH"/cteq6/"
#define CTEQ65S_TBL_PATH     CTEQ_TBL_PATH"/ctq65s/"
#define CTEQ65C_TBL_PATH     CTEQ_TBL_PATH"/ctq65c/"
#define CTEQ6M_TBL_PATH      CTEQ_TBL_PATH"/cteq6m/"
#define CTEQ61_TBL_PATH      CTEQ_TBL_PATH"/cteq61/"
#define CTEQ65_PDS_TBL_PATH  CTEQ_TBL_PATH"/ctq65-pds/"
#define CTEQ66A_TBL_PATH     CTEQ_TBL_PATH"/ctq66a/"
#define CTEQ66C_TBL_PATH     CTEQ_TBL_PATH"/ctq66c/"
#define CTEQ66M_TBL_PATH     CTEQ_TBL_PATH"/ctq66m/"



static const cteq_pdfset_t __cteq_pdfset_database[] = { 
  {400, "CTEQ66.00",  "description", CTEQ66M_TBL_PATH, "ctq66.00.pds", 2},
  
  /*  End of the list */
  {0,0,0,0,0,2}
};



const cteq_pdfset_t *cteq_pdfset_database = __cteq_pdfset_database;


const cteq_pdfset_t * 
cteq_pdfset_find(const cteq_pdfset_t *pdflist, const char *name)
{
  while(pdflist->name) {
	if(strcmp(name, pdflist->name) == 0) return pdflist;
	++pdflist;
  }

  return 0;
}

const cteq_pdfset_t * 
cteq_pdfset_find_id(const cteq_pdfset_t *pdflist, int id)
{
  while(pdflist->name) {
	if(id == pdflist->id) return pdflist;
	++pdflist;
  }
  
  return 0;
}







