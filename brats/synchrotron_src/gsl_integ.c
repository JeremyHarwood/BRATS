/* wrapper code for GSL integration */

/* this makes the GSL function look a little more like our original qromo calls */

#include "gsl_integ.h"
#include <stdio.h>

double gsl_epsilon=GSL_EPS;
int gsl_key=GSL_KEY;
int gsl_quiet=0;

double gsl_integ(double f(),double x1,double x2,gsl_integration_workspace *w) {

  double res,err;
  gsl_function F;
  int r;

  F.function=f;
  F.params=NULL;
  //printf("gsl_integ called: %g %g\n",x1,x2);
  r=gsl_integration_qag(&F,x1,x2,0,gsl_epsilon,GSL_WSIZE,gsl_key,w,&res,&err);
  if (!gsl_quiet && r) fprintf(stderr,"gsl_integ error %i: returning %g, err is %g\n",r,res,err);
  return res;
}

double gsl_integ_i(double f(),gsl_integration_workspace *w) {

  double res,err;
  gsl_function F;
  int r;

  F.function=f;
  F.params=NULL;
  //printf("gsl_integ called: %g %g\n",x1,x2);
  r=gsl_integration_qagi(&F,0,GSL_EPS,GSL_WSIZE,w,&res,&err);
  if (!gsl_quiet && r) fprintf(stderr,"gsl_integ_i error %i: returning %g, err is %g\n",r,res,err);
  return res;
}

double gsl_integ_iu(double f(),double x1,gsl_integration_workspace *w) {

  double res,err;
  gsl_function F;
  int r;

  F.function=f;
  F.params=NULL;
  //printf("gsl_integ called: %g %g\n",x1,x2);
  r=gsl_integration_qagiu(&F,x1,0,GSL_EPS,GSL_WSIZE,w,&res,&err);
  if (!gsl_quiet && r) fprintf(stderr,"gsl_integ_iu error %i: returning %g, err is %g\n",r,res,err);
  return res;
}
/*
double gsl_integ_qng(double f(),double x1,double x2) {

  double res,err,n;
  gsl_function F;
  int r;

  F.function=f;
  F.params=NULL;
  //printf("gsl_integ called: %g %g\n",x1,x2);
  r=gsl_integration_qng(&F,x1,x2,0,gsl_epsilon,&res,&err,&n);
  if (!gsl_quiet && r) fprintf(stderr,"gsl_integ error %i: returning %g, err is %g\n",r,res,err);
  return res;
}
*/
