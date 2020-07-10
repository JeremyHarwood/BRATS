#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#define GSL_EPS 0.5e-4
#define GSL_WSIZE 10000
#define GSL_KEY GSL_INTEG_GAUSS21


extern double gsl_epsilon;
extern int gsl_key;
extern int gsl_quiet;
	
//#pragma omp threadprivate(gsl_epsilon, gsl_key, gsl_quiet)


double gsl_integ(double (),double,double,gsl_integration_workspace *);
double gsl_integ_i(double (),gsl_integration_workspace *);
double gsl_integ_iu(double (),double,gsl_integration_workspace *);
//double gsl_integ_qng(double (),double,double);

