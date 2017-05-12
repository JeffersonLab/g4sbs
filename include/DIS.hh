#ifndef G4SBSDIS_HH
#define G4SBSDIS_HH

#include "cteqpdf.h"
//#include "sbstypes.hh"
#include <assert.h>
// Use CTEQ6 parameterization

cteq_pdf_t *__dis_pdf;

void initcteqpdf(){
     __dis_pdf = cteq_pdf_alloc_id(400); // mode 400 = cteq6.6?
     assert(__dis_pdf);
}

double F2N(double x, double Q2,  int nucl){
  
    double qu = cteq_pdf_evolvepdf(__dis_pdf, 1, x, sqrt(Q2) );
    double qd = cteq_pdf_evolvepdf(__dis_pdf, 2, x, sqrt(Q2) );
    double qubar = cteq_pdf_evolvepdf(__dis_pdf, -1, x, sqrt(Q2) );
    double qdbar = cteq_pdf_evolvepdf(__dis_pdf, -2, x, sqrt(Q2) );

    double quv = qu-qubar;
    double qdv = qd-qdbar;

    double qs = cteq_pdf_evolvepdf(__dis_pdf, 3, x, sqrt(Q2) );

    double F2 = 0.0; 
    double e_u =  2.0/3.0;
    double e_d = -1.0/3.0;

    if( !nucl ){
    	F2 += x*( e_u*e_u*quv + e_d*e_d*qdv ); 
    }
    else {
	F2 += x*( e_u*e_u*qdv + e_d*e_d*quv ); 
    }
    // Sea quarks
    F2  += x*(2.0*e_u*e_u*qubar + 2.0*e_d*e_d*(qdbar + qs));

    return F2;

}

double dissigma( double ebeam, double th, double eprime, int nucl ){
    // Return in nb

    double Q2 = 2.0*eprime*ebeam*(1.0-cos(th));
    double nu = ebeam-eprime;
    double Mp = 0.938;

    double x = Q2/(2.0*Mp*nu);
    double y = nu/ebeam;

    if( ! (0.0 < x && x < 1.0 && 0.0 < y && y < 1.0) ){
      //printf("WARNING %s line %d  x = %f, y = %f -> eprime = %f GeV   th = %f deg  ebeam = %f GeV\n", __FILE__,
      //__LINE__, x, y, eprime, th*180/3.14159, ebeam );
      //exit(1);
	return 0.0;;

    }

    double qu = cteq_pdf_evolvepdf(__dis_pdf, 1, x, sqrt(Q2) );
    double qd = cteq_pdf_evolvepdf(__dis_pdf, 2, x, sqrt(Q2) );
    double qubar = cteq_pdf_evolvepdf(__dis_pdf, -1, x, sqrt(Q2) );
    double qdbar = cteq_pdf_evolvepdf(__dis_pdf, -2, x, sqrt(Q2) );

    double quv = qu-qubar;
    double qdv = qd-qdbar;

    double qs = cteq_pdf_evolvepdf(__dis_pdf, 3, x, sqrt(Q2) );

    double F2 = 0.0; 
    double e_u =  2.0/3.0;
    double e_d = -1.0/3.0;

    if( !nucl ){
	F2 += x*( e_u*e_u*quv + e_d*e_d*qdv ); 
    }
    else {
	F2 += x*( e_u*e_u*qdv + e_d*e_d*quv ); 
    }
    // Sea quarks
    F2  += x*(2.0*e_u*e_u*qubar + 2.0*e_d*e_d*(qdbar + qs));
    double F1 = F2/(2.0*x);

    // From PDG
    double ds_dxdy
	= 4.0*3.14159*((1.0-y-pow(x*y*Mp,2.0)/Q2)*F2+y*y*x*F1)
	          /(x*y*Q2*137.0*137.0);

    // In GeV^-2
    double ds_dOmega_dE = ds_dxdy*eprime/(2.0*3.14159*Mp*nu);

    return ds_dOmega_dE*0.197*0.197*1e7; // GeV2 -> nb
}

double dissigma_p(double eb, double th, double ep){
    return dissigma( eb, th, ep, 0);
}
double dissigma_n(double eb, double th, double ep){
    return dissigma( eb, th, ep, 1 );
}

#endif//G4SBSDIS_HH
