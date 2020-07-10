#include <stdio.h>
#include <math.h>
#include <gsl_integ.h>
//#include <gsl_spline.h>
#include <stdlib.h>
#include <constants.h>
#include <bessel.h>
#include <omp.h>
#include <gsl/gsl_cdf.h>

/* synchrotron code comes from SoS */

#define MAXDATA 100
#define XF1 1.0e-4
#define XF2 0.22e+2
#define YF1 1.0e-3
#define YF2 0.20e+2
#define INFIN  1.0e2

#ifndef INTCO
#define INTCO (3*ECF)/(4*PI*M_EL*M_EL*M_EL*V_C*V_C*V_C*V_C)
#define EMCO (sqrt(3)*ECF*ECF*ECF)/(4.0*PI*EPS_0*V_C*M_EL)
#endif

static double axf1, axf2, ayf1, ayf2;
static double *fx, *fy;
//static int xls=100;

double fr[MAXDATA],em[MAXDATA],eme[MAXDATA];
double (*ne)(double);
double (*ng)(double);
double (*age_gmax)(double);

double gmin,gmax,power,gbreak,gbvalue;
double GMIN,GMAX,EBREAK;
double ageb;
double ebrd,gbrd;

/* the above variables are OK as globals because they are not changed by threads. On the other hand, the following ones are! */

double n0_ext;
double BFIELD, nu;
double sinalpha;
double mb_a;
double yy;
double age;
double t_off;
//t_off = 2*(1e6*365*86400);

#pragma omp threadprivate(n0_ext,BFIELD,nu,sinalpha,mb_a, age, yy,t_off)

gsl_integration_workspace *w1,*w2,*w3,*w4;

#pragma omp threadprivate(w1,w2,w3,w4)


double gbreak_kg(double maxg) {
    
  double gammabreak;

  gammabreak = 1 / (((4*THOMSON/(6*M_EL*V_C*MU_0))*(ageb*ageb)) * (age+t_off));
  if (gammabreak>maxg) gammabreak=maxg;
  //printf("Break at: %.4e (nu of: %.4e)\n", gammabreak, (ECF*BFIELD*gammabreak*gammabreak)/(2*PI*M_EL));

  return(gammabreak);
}


double ne_age_kg(double g) {
    
  /* ne for a finite age continuous injection model */
  /* loss is the term that corrects for loss. */
  /* KG version */

  double gdepterm = 0.0;
  double loss = 0.0;
  
  if ((g<GMIN) || (g>GMAX)) {
    return(0.0);
  }
  else loss = (4*THOMSON/(6*M_EL*V_C*MU_0))*(ageb*ageb); // C1, b and E [sic: gamma] of LG. Note we dont include the time component here compared to standard JP as there are 2 t terms needed later (age and t_off)

    
  if (g < 1/(loss*(age+t_off))) {
    gdepterm = pow(1.0-(loss*g*t_off), (power-1.0)) - pow(1.0-(loss*g*(age+t_off)), (power-1.0));
  }
  else if (g < 1/(loss*t_off) ) {
    gdepterm = pow(1.0-(loss*g*t_off), (power-1.0));
  }
    
    
  //printf("edepterm: %.4e n0_ext: %.4e ageb %.4e e: %.4e exponent: %.4e returnval: %.4e\n", edepterm, n0_ext, ageb, e, (-power-1)/(loss*(power-1)),((n0_ext * edepterm * pow(e, (-power-1)))/(loss*(power-1))) );
    
  return( n0_ext * gdepterm * pow(g, (-power-1))/(loss*(power-1)) );
    
}

double age_gmax_kg(double limit) {
  
  double gm;
  double loss = (4*THOMSON/(6*M_EL*V_C*MU_0))*(ageb*ageb);
    
  gm = 1/(loss*t_off);

  if (gm>limit) gm=limit;
  return(gm);
}


double ffy_kg(double y) {
  /* use the tabulated values of F(y) to return F(y) at any y. For y
     outside the table, use the asymptote of Pacholczyk (1970) for small
     y and 0 for large y. Otherwise return a log-linear two-point
     interpolation. */
    
  int i;
  double av,yv1,yv2, c1, c2, c3, c4;

  c1 = 1.808418021;
  c2 = -1.813799364;
  c3 = 0.8476959474;
  c4 = -0.510131;
    
  if (y>=YF2)
    return 0.0;
  else if (y<=YF1)
    return (c1 * pow(y,(1/3))) + (c2 * y) + (c3 * pow(y,(7/3))) + (c4 * pow(y,3))  ; // Asymptote for small values of y
  else {
    av=log(y);
    yv1=(av-ayf1)*((double)xls-1.0)/(ayf2-ayf1);
    yv2=floor(yv1);
    i=yv2;
    yv1-=yv2;
    return exp(((1.0-yv1)*log(fy[i]))+(yv1*log(fy[i+1])));
  }
}					       


double ff_y_kg(double t) {
    
  float bessk, dummy1, dummy2, dummy3, ff_y_res;
    
  /* given x, returns the Bessel function of order 5/3 */
  bessik(t,5.0/3.0,&dummy1,&bessk,&dummy2,&dummy3);
  /*  printf("x was %g, Bessel result was %g\n",x,bessk); */

  ff_y_res = sqrt(1 - (pow(yy,2)/pow(t,2))) * bessk;
  return ff_y_res;
}


void makefy_kg(void) {
    
  int i;
  double yv, av;
    
  //  printf("Making F(y) lookup table:\n");

  ayf1=log(YF1);
  ayf2=log(YF2);
#pragma omp parallel for private(i,av,yv)
  for (i=0; i<xls; i++) {
    /* compute the y-value */
    av=ayf1+((ayf2-ayf1)*((double)i)/(double)(xls-1));
    yy=exp(av);
    /* integrate ff from xv to infinity */
    yv=gsl_integ_iu(ff_y_kg,yy,w1);
    fy[i]=yy*yv;
  }
}

double integr_kg(double g) {
    
  /* the synchrotron emission integrand */
  double omega_0 = (ECF*BFIELD) / M_EL;
  return ne(g)*ffy_kg((4*PI*nu)/(3*omega_0*g*g));
}


double synchin_kg() {
    
  double gmin,gbreak,gmax,resI,resII;
  double omega_0 = (ECF*BFIELD) / M_EL;
  double classical_e_radius = 2.8179403267e-15;
       
  gmin=GMIN;
  gbreak=gbreak_kg(GMAX);
  gmax=GMAX;
  //printf("gbreak: %.4e \n", gbreak);

  if (gmin<GMIN) gmin=GMIN;
  if (gmax>GMAX) gmax=GMAX;
    
  if ((age+t_off)>0.0) gmax=age_gmax_kg(GMAX);

  if (gmin>gbreak) {
    resI=0.0;
  }
  else {
    resI=gsl_integ(integr_kg,gmin,gbreak,w1);
  }

  if (resI<0) {
    printf("synchrotron integral convergence problem (synchin_kg_I)\n");
    resI=0.0;
  }

  if (gbreak>=gmax) {
    resII=0;
  }
  else {
    resII=gsl_integ(integr_kg,gbreak,gmax,w2);
  }
  //printf("resII: %.4e\n", resI);

  if (resII<0) {
    printf("synchrotron integral convergence problem (synchin_kg_II)\n");
    resII=0.0;
  }
    
  return ( (sqrt(3.0)*M_EL*V_C*classical_e_radius) / (4*PI) )*omega_0*(resI+resII);
}


double emiss_nkg(double n0, double b, double nu_arg) {
    
  double emiss;
    
  n0_ext=n0;
  BFIELD=b;
  nu=nu_arg;

  emiss = 1e26*synchin_kg();

  //printf("age: %.4e frequency: %.4e totalemission: %.4e\n", age, nu, emiss);
    
  return emiss;
}




int ciageingmodels(float **regflux, int imgnum, int *fluxorder, float *frequency, int passingreg, float *chisquared, float *bestage, float *bestoff, float **fluxerror, float *modelalpha, float **modelflux, double usr_gmin, double usr_gmax, int myears, int minmyears, int minoff, int maxoff, int ageresolution, int levels, float *bestnorm, int modelres, int printresults, double inject, double fieldstrength, int model, int printreject, float redshift, float *offerrorsplus, float *offerrorsminus, float *onerrorsplus, float *onerrorsminus, int ul, int *upperlimits, float *bcmb, float *passageingb, float *breakon, float *breakoff, float *conflvl, int *suppress, char *index, int export, int suppresscdf, int extrapolatemodel, float extrapolationfrequency, float fixedage) {


  //#define STEPS 50

  // Declare the function prototypes
  double normgoldensearch_func(double normflux, int imgnum, float **fluxerror, float **regflux, int currentregion, int *fluxorder, float *frequency, double magfield, double *emibymapfreq, int ul, int *upperlimits); 

  int goldensearch( double (*f)(double, int, float **, float **, int, int *, float *, double, double *, int, int *), double* a, double *fa, double* b, double* fb, double tolerance, int imgnum, float **fluxerror, float **regflux, int currentregion, int *fluxorder, float *frequency, double magfield, double *emibymapfreq, int ul, int *upperlimits);


  double b, emi, interval, siglevel68, siglevel90, siglevel95, siglevel99, siglevel995, siglevel999, siglevel9999, dof, c, d, fc, fd, ***emibymapfreq;
  float nu, normfactor, worstfit, normmax, normmin, tmp_maxage, tmp_minage, tmp_maxage_off, tmp_minage_off, confpercentage, maxfitage, maxfitage_off;
  int i, a, j, m, l, s, u, g, currentregion, *medianbin, median, rejectbin, nonrejectbin, oldestregion, freeparams;


  // Setup some variables
  int numregions = 1; // This is for the dummy array and purely for legacy purposes
  int resolution = ageresolution;
  int offresolution = 0; // Varies later dependent on the model used
  int regcounter = 0;
  int hundreds = 0;
  int parallelabort = 0;
  int onegoodfit = 0;
  float steps = 0.0;
  float tmp_steps = 0.0;
  float offsteps = 0.0;
  float tmp_offsteps = 0.0;
  float tmp_bestchi = 1e22;
  float tmp_bestage = 0.0;
  float tmp_bestoffchi = 1e22;
  float tmp_bestoffage = 0.0;
  float avgchisquared = 0.0;
  float avgnormalisation = 0.0;
  float sumchisquared = 0.0;
  float freqsteps = 0.0;
  //float rangecorrection = 0.0;

  // These are changed again later in the loop
  float maxage = (float)myears;
  float minage = (float)minmyears;
  float maxage_off = (float)maxoff;
  float minage_off = (float)minoff;

  float fixbybreak = -1.0;


  double tolerance = -1;

  double *tmp_chibyon, *tmp_chion, *tmp_chibyoff, *tmp_chioff, tmp_deltachiof1, tmp_sort, slope, intercept, tmp_chisquared;
  int cub_cnt = 0;

  int beststepindex = 0;
  int tmp_beststepindex = 0;
  int tmp_bestoffindex = 0;


  t_off = 1e-10*(1e6*365*86400); // Stops overflow when too close to 0. This equates to less than an hour switched off in real time.

  medianbin = (int *)calloc(5, sizeof(int));

  // If we are using a fixed age set the search range to make just 1 attempt. The variables should remain local in scope.
  if (fixedage > 0.0) {
    levels = 1;
    resolution = 0;
  }

  // Get the confidence level ready for later
  if (model == 5) { // CI Off Model
    freeparams = 3;
  }
  else if (model == 6) { // CI Model
    freeparams = 2;
  }
  else {
    fprintf(stderr,"\nError: Model number is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");
    return(6);
  }
  
  dof = (imgnum - ul - freeparams);

  if (suppresscdf !=1) {
    // Calculate the confidence levels
    interval = 0.68;
    siglevel68 = gsl_cdf_chisq_Pinv(interval, dof);
    interval = 0.90;
    siglevel90 = gsl_cdf_chisq_Pinv(interval, dof);
    interval = 0.95;
    siglevel95 = gsl_cdf_chisq_Pinv(interval, dof);
    interval = 0.99;
    siglevel99 = gsl_cdf_chisq_Pinv(interval, dof);
    interval = 0.995;
    siglevel995 = gsl_cdf_chisq_Pinv(interval, dof);
    interval = 0.999;
    siglevel999 = gsl_cdf_chisq_Pinv(interval, dof);
    interval = 0.9999;
    siglevel9999 = gsl_cdf_chisq_Pinv(interval, dof);
  }

  // Start the fitting process
  if (export == 0) {
    for(j=0; j<imgnum; j++) {
      printf("Rest frequency: %.3e\n", frequency[fluxorder[j]]);
    }
  }

#pragma omp parallel 
  {
    w1=gsl_integration_workspace_alloc(GSL_WSIZE);
    w2=gsl_integration_workspace_alloc(GSL_WSIZE);
    w3=gsl_integration_workspace_alloc(GSL_WSIZE);
  }

  gsl_set_error_handler_off(); // to avoid abort
  gsl_epsilon=1e-3;
  gsl_quiet=1;

  fy=calloc(xls,sizeof(double));
  makefy_kg();

  // printf("Ready to go!\n");

  // Now user set, default 10 and 100000
  
  gmin=usr_gmin;
  gmax=usr_gmax;
  
  // Set which model energy distribution to use (always KG at the moment, but allows for additional similar models later).
  if (model == 5) { // CI Off Model
    ne=ne_age_kg;
    age_gmax=age_gmax_kg;

  }
  else if (model == 6) { // CI Model
    ne=ne_age_kg;
    age_gmax=age_gmax_kg;
  }
  else {
    fprintf(stderr,"\nError: Model number is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");
    return 6;
  }
  
  // Default of 2.2 gives a spetcral index of 0.6
  power=(2*inject)+1; // Now a user input

  GMIN=gmin;
  GMAX=gmax;

  b=fieldstrength;
  bcmb[passingreg]=(0.318 * (pow(1+redshift,2))) * 1e-9;
  ageb=sqrt(pow(b,2) + pow(bcmb[passingreg],2));
  passageingb[passingreg]=ageb; // Pass back the ageing field to the calling module (ageb needs to be local)
  //ageb=b;

  if (export == 0) {
    printf("Magnetic field strength: %.3e T\n", b);
    printf("B_CMB: %.3e T\n", bcmb[passingreg]);
    printf("Ageing field strength: %.3e T\n", ageb);
  }

  // Get the model values for only those points which have real data available and carry out a chi-square test

  // printf("Calculating model fits (this may take some time)\n Currently processing region: 1 of %d ...\n", numregions);
  
  emibymapfreq = (double ***)calloc(resolution+1, sizeof(double **));

  for (i=0; i<=resolution; i++){
    emibymapfreq[i] = (double **)calloc(resolution+1, sizeof(double *));
    for (j=0; j<=resolution; j++){
      emibymapfreq[i][j] = (double *)calloc(imgnum, sizeof(double));
    }
  }

  b=fieldstrength;
  bcmb[passingreg]=(0.318 * (pow(1+redshift,2))) * 1e-9;
  ageb=sqrt(pow(b,2) + pow(bcmb[passingreg],2));
  //ageb=b;
 
  for(a=1; a<=numregions; a++) {

    tmp_chibyon = (double *)calloc((resolution+1)*levels, sizeof(double));
    tmp_chion = (double *)calloc((resolution+1)*levels, sizeof(double));
    tmp_chibyoff = (double *)calloc((resolution+1)*levels, sizeof(double));
    tmp_chioff = (double *)calloc((resolution+1)*levels, sizeof(double));

    onegoodfit = 0;
  
    if (parallelabort != 1) { // If we have bad values, stop trying

      // Set the default chisquared to an unlikely values
      chisquared[a] = 5e22;
      worstfit = -1e22;

      if (fixedage > 0.0) {
	minage = fixedage;
	maxage = fixedage;

	minage_off = (float)minoff; // These can be changed to allow for a fixed off time later if there's a demand for it.
	maxage_off = (float)maxoff;
      }
      else {
	if (minmyears == 0) {
	  minage=0.01; // For axis and normalization factor later on
	}
	else {
	  minage = (float)minmyears;
	}
	minage_off = (float)minoff;
	maxage = (float)myears;
	maxage_off = (float)maxoff;

	// Allow varying by break. Currently not implemented, debugging only.
	if (fixbybreak > 0) {

	  //double fixb, rootdenom;
	  //fixb = (4*THOMSON/(6*M_EL*V_C*MU_0))*(ageb*ageb);
	  //rootdenom = sqrt((M_EL*PI*fixbybreak)/(b*ECF));
	  //minage = sqrt(2)/(2*fixb*rootdenom);
	  minage = sqrt(2)/(2*((4*THOMSON/(6*M_EL*V_C*MU_0))*(ageb*ageb))*sqrt((M_EL*PI*fixbybreak)/(b*ECF)));
	  minage /= (1e6*365*86400);
	  maxage = minage;
	  //printf("fixbybreak: %.4e ageb: %.4e fixb: %.4e rootdenom: %.4e fixedage: %.4e\n", fixbybreak, ageb, fixb, rootdenom, minage);
	}
      }

      // For each level of increasing accuracy (changes done at closing bracket)
      for(l=0; l<levels; l++) {

	if (parallelabort != 1) { // If we have bad values, stop trying
	  
#pragma omp parallel for private(j, m, g, steps, offsteps, normfactor, normmax, normmin, tmp_minage, tmp_maxage, tmp_minage_off, tmp_maxage_off,c,d,fc,fd) firstprivate(minage, maxage, minage_off, maxage_off) schedule(dynamic)

	  for(j=0; j<=resolution; j++ ) {

	    if (parallelabort != 1) { // If we have bad values, stop trying

	      if (j == 0) {
		steps = minage;
	      }
	      else {
		steps = minage + ( ((maxage-minage)/resolution) * (float)j );
	      }

	      // Age in seconds
	      age=1e6*steps*365*86400;

	      // Go through each map and find the range of possible normalisation values
	      normmax = -1e36;
	      normmin = 1e36;

	      if (model == 6) {
		offresolution = 0;
	      }
	      else if (model == 5) {
		offresolution = resolution;
	      }

	      //#pragma omp parallel for private(g, t_off, offsteps) schedule(dynamic)
	      for(g=0; g<=offresolution; g++) {

		if (parallelabort != 1) {

		  if (model == 5) {

		    if (g == 0) {
		      offsteps = minage_off;
		    }
		    else {
		      offsteps = minage_off + ( ((maxage_off-minage_off)/offresolution) * (float)g );
		    }
		    // Age in seconds
		    t_off=1e6*offsteps*365*86400;
		  }
		  else if (model == 6) {
		    t_off=0.0;
		    offsteps = 0.0; // Needs initializing to keep valgrind and the compiler happy.
		  }
		
#pragma omp parallel for private(m) schedule(dynamic)
		  for(m=0; m<imgnum; m++) {

		    // Stop trying if we have bad values
		    if (parallelabort != 1) {
		    
		      if ( (model == 5) || (model == 6) ) {
			emibymapfreq[j][g][m]=emiss_nkg(1,b,frequency[fluxorder[m]]);
		      }
		      else {
			emibymapfreq[j][g][m]=emiss_nkg(1,b,frequency[fluxorder[m]]);
		      }
		    		    
		      // Check is the flux value is bad and if so make sure the message is only displayed once
		      if ( ( (emibymapfreq[j][g][m] <= 1e-33) || (emibymapfreq[j][g][m] >= 1e33) ) && (parallelabort!=1) ) {   
			if (*suppress == 0) {
			  printf("\nError: Model has reached its precision limit (model flux effectively 0 or unphysically high) at an age of %.2f (on) %.2f (off) Megayears. Please run plot(ci/cioff) to determine a suitable age range and / or model paramaters (e.g. bfield)\n", steps, offsteps);
			  parallelabort = 1;
			}
			else {
			  // Check this isnt the very first try!
			  if ( (g == 0) && (j == 0) ) {
			    printf("\nError: Model has reached its precision limit (model flux effectively 0 or unphysically high) on the very first attempt! This cannot be suppressed. No output will made for source with index %s.\n", index);
			    *suppress = 0;
			    parallelabort = 1;
			  }
			  else {     
			  // Increment to indicate we had to suppress something (doesnt matter if we have mutiple increments as we are only interested if it is > 1)
			    (*suppress)++;
			    parallelabort = 1;
			  }
			}
		      }

		      normfactor = regflux[fluxorder[m]][a]/(emibymapfreq[j][g][m]);
		      normfactor = log10(normfactor);
		      //rangecorrection = normfactor * (fluxerror[m][a]/regflux[fluxorder[m]][a]);
		      
		      if ( normfactor > normmax ) {
			normmax = normfactor;
		      }
		      
		      if ( normfactor < normmin ) {
			normmin = normfactor;
		      }
		    }
		  }
		  
		  normmax = pow(10,normmax) + (pow(10,normmax)*0.5);
		  normmin = pow(10,normmin) - (pow(10,normmin)*0.5);
		  normmax = log10(normmax);
		  normmin = log10(normmin);
		  
		  if (parallelabort != 1) {
		    
		    currentregion = a;
		    c = normmin;                                                 
		    d = normmax;
		    
		    fc = normgoldensearch_func(c, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, b, emibymapfreq[j][g], ul, upperlimits);
		    fd = normgoldensearch_func(d, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, b, emibymapfreq[j][g], ul, upperlimits);
		    
		    goldensearch(normgoldensearch_func, &c, &fc, &d, &fd, tolerance, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, b, emibymapfreq[j][g], ul, upperlimits);
		    
		    // For testing vs Leith's fit. Also change the loop in the normgoldenrationsearch from main synch module to imgnum-1

		    //printf("emibymapfreq[j][g][4]: %.4e *c: %.4e *d: %.4e\n",emibymapfreq[j][g][4],emibymapfreq[j][g][4]*pow(10, c),emibymapfreq[j][g][4]*pow(10, d));
		    
		    /*if ( (fc < fd) && (emibymapfreq[j][g][4]*pow(10, c) >= 0.03) ) {
		      fc = 1e33;
		      fd = 1e33;
		      //printf("Over UL (c)\n");
		      }
		      else if ( (fd < fc) && (emibymapfreq[j][g][4]*pow(10, d) >= 0.03) ) {
		      fc = 1e33;
		      fd = 1e33;
		      //printf("Over UL (d)\n");
		      }
		      else if (onegoodfit == 0) {
		      onegoodfit = 1;
		      }*/
		  
		    if ((ul > 0) && (onegoodfit == 0)) {
		      if ( (fc < fd) && (emibymapfreq[j][g][4]*pow(10, c) < 0.03) ) {
			onegoodfit = 1;
		      }
		      else if ( (fd < fc) && (emibymapfreq[j][g][4]*pow(10, d) < 0.03) ) {
			onegoodfit = 1;
		      }
		    }
		    else if (onegoodfit == 0) {
		      onegoodfit = 1;
		    }
		  
		    //If the chi squared is less than the previous best, set its new value and age
		    if ( (fc < chisquared[a]) || (fd < chisquared[a]) ) {
		    
		      if (fc < fd) {
			tmp_chisquared = (float)fc;
			bestnorm[a] = pow(10, c);
		      }
		      else {
			tmp_chisquared = (float)fd;
			bestnorm[a] = pow(10, d);
		      }
		    
		      chisquared[a] = tmp_chisquared;
		      bestage[a] = steps;
		      bestoff[a] = offsteps;
		      beststepindex = j+(l*(resolution+1));
		    }
		    
		    if ( (fc > worstfit) || (fd > worstfit) ) {
		      
		      if (fc > fd) {
			tmp_chisquared = (float)fc;
		      }
		      else {
			tmp_chisquared = (float)fd;
		      }
		    
		      worstfit = tmp_chisquared;
		      //worstage = steps;
		    
		    }
		  }
		}
	      }
	    }
	  }

	  if (parallelabort != 1) {
	    
	    // Set the new min and max ages
	    tmp_maxage = maxage;
	    tmp_minage = minage;

	    tmp_maxage_off = maxage_off;
	    tmp_minage_off = minage_off;
	  
	    // Set the new min and max levels
	    if ( ( bestage[a] - ( (  ( (tmp_maxage-tmp_minage) / resolution ) /2)*pow(10, -l) ) ) <= 0.01) { // Prevent an impossible age
	      minage = 0.01;
	    }
	    else {
	      minage = bestage[a] - ( (  ( (tmp_maxage-tmp_minage) / resolution ) /2)*pow(10, -l) );
	    }
	    maxage = bestage[a] + ( (  ( (tmp_maxage-tmp_minage) / resolution ) /2)*pow(10, -l) );	  

	    // Set the new min and max levels for t_off
	    if ( ( bestoff[a] - ( (  ( (tmp_maxage_off-tmp_minage_off) / resolution ) /2)*pow(10, -l) ) ) <= 0) { // Prevent an impossible age
	      minage_off = 0.0;
	    }
	    else {
	      minage_off = bestoff[a] - ( (  ( (tmp_maxage_off-tmp_minage_off) / resolution ) /2)*pow(10, -l) );
	    }
	    maxage_off = bestoff[a] + ( (  ( (tmp_maxage_off-tmp_minage_off) / resolution ) /2)*pow(10, -l) );
	  }
	}
      
	// For use with strict upper limts (Leith's method)
	if ( (onegoodfit != 1) && (ul != 0) ) { // if we are always excluded by the upper limit, escape

	  printf("\n*** Error: Unable to find a model for the given parameters which fits below the all the upper limits. ***\n--> Skipping data set with index %s...\n\n", index);
	  parallelabort = 1;
	  return 5;
	}
      }
     
      // If we have bad values and we are not suppressing, kick us out
      if ( (*suppress == 0) && (parallelabort == 1) ) {
	return 5;
      }
  
      if (printresults == 1) {
	printf("Best fit for %d: %.2f Myr (On), %.2f Myr (Off) Reduced X^2 %.2f and Best Norm: %.2e \n", a, bestage[a], bestoff[a], ( chisquared[a] / (imgnum - ul - freeparams) ), bestnorm[a]);
      }
  	
      //printf("local minimum is between (%12.6e,%12.6e) and (%12.6e,%12.6e)\n", c,fc, d,fd);
  
      regcounter++;
      avgchisquared += chisquared[a];
      avgnormalisation += bestnorm[a];

      if (suppresscdf !=1) {
	// Bin the chisquared region values by confidence level
	if (chisquared[a] < siglevel68) {
	  conflvl[passingreg] = 1.0000;
	  medianbin[0]++;
	}
	else if ( (chisquared[a] >= siglevel68) && (chisquared[a] < siglevel90) ) {
	  conflvl[passingreg] = 68.0000;
	  medianbin[1]++;
	}
	else if ( (chisquared[a] >= siglevel90) && (chisquared[a] < siglevel95) ) {
	  conflvl[passingreg] = 90.0000;
	  medianbin[2]++;
	}
	else if ( (chisquared[a] >= siglevel95) && (chisquared[a] < siglevel99) ) {
	  conflvl[passingreg] = 95.0000;
	  medianbin[3]++;
	}
	else if (chisquared[a] >= siglevel99) {
	  conflvl[passingreg] = 99.0000;
	  medianbin[4]++;
	}
      }

      if (regcounter >= 100) {

	hundreds++;
	printf("%d00...\n", hundreds);

	regcounter = 0;
      }
    }
    
    
    // If we have bad values, kick us out
    if ( (*suppress == 0) && (parallelabort == 1) ) {
      return 5;
    }
  
    // Find the on and off errors. For this we need to start looping through again to get the individual X^2 curves.
    g=0;
    t_off=1e6*bestoff[a]*365*86400; // Fix the t_off to the best age
    tmp_bestchi = 1e26;
    tmp_bestage = 0.0;

    if (minmyears == 0) {
      minage=0.01; // For axis and normalization factor later on
    }
    else {
      minage = (float)minmyears;
    }
    maxage = (float)myears;

    // Allow varying by break. Currently not implemented, debugging only.
    if (fixbybreak > 0) {
      minage = sqrt(2)/(2*((4*THOMSON/(6*M_EL*V_C*MU_0))*(ageb*ageb))*sqrt((M_EL*PI*fixbybreak)/(b*ECF)));
      minage /= (1e6*365*86400);
      maxage = minage;
      //printf("fixbybreak: %.4e ageb: %.4e fixb: %.4e rootdenom: %.4e fixedage: %.4e\n", fixbybreak, ageb, fixb, rootdenom, minage);
    }

    // Find the X^2 curves for on ages
    for(l=0; l<levels; l++) {

#pragma omp parallel for private(j, m, steps, tmp_steps, normfactor, normmax, normmin, tmp_minage, tmp_maxage, c,d,fc,fd) firstprivate(minage, maxage) schedule(dynamic)
      for(j=0; j<=resolution; j++) {

	normmax = -1e36;
	normmin = 1e36;

	if (parallelabort != 1) {

	  if (j == 0) {
	    tmp_steps = minage;
	  }
	  else {
	    tmp_steps = minage + ( ((maxage-minage)/resolution) * (float)j );
	  }
		    
	  // Age in seconds
	  age=1e6*tmp_steps*365*86400;

	  for(m=0; m<imgnum; m++) {
		    
	    // Stop trying if we have bad values
	    if (parallelabort != 1) {
	      
	      if ( (model == 5) || (model == 6) ) {
		emibymapfreq[j][g][m]=emiss_nkg(1,b,frequency[fluxorder[m]]);
	      }
	      else {
		emibymapfreq[j][g][m]=emiss_nkg(1,b,frequency[fluxorder[m]]);
	      }

	      // Check is the flux value is bad and if so make sure the message is only displayed once
	      if ( ( (emibymapfreq[j][g][m] <= 1e-33) || (emibymapfreq[j][g][m] >= 1e33) ) && (parallelabort!=1) ) {

		if (*suppress == 0) {
		  printf("\nError: Model has reached its precision limit (model flux effectively 0 or unphysically high) at an age of %.2f (on) %.2f (off) Megayears. Please run plot(ci/cioff) to determine a suitable age range and / or model paramaters (e.g. bfield)\n", tmp_steps, tmp_offsteps);
		  parallelabort = 1;
		}
		else {
		  // Check this isnt the very first try!
		  if (j == 0) {
		    printf("\nError: Model has reached its precision limit (model flux effectively 0 or unphysically high) on the very first attempt! This cannot be suppressed. No output will made for source with index %s.\n", index);
		    *suppress = 0;
		    parallelabort = 1;
		  }
		  else {     
		    // Increment to indicate we had to suppress something (doesnt matter if we have mutiple increments as we are only interested if it is > 1)
		    (*suppress)++;
		    parallelabort = 1;
		  }
		}		
	      }

	      normfactor = regflux[fluxorder[m]][a]/(emibymapfreq[j][g][m]);
	      normfactor = log10(normfactor);
		      		      
	      if ( normfactor > normmax ) {
		normmax = normfactor;
	      }
		      
	      if ( normfactor < normmin ) {
		normmin = normfactor;
	      }
	    }
	  }

	  normmax = pow(10,normmax) + (pow(10,normmax)*0.5);
	  normmin = pow(10,normmin) - (pow(10,normmin)*0.5);
	  normmax = log10(normmax);
	  normmin = log10(normmin);
		  
	  if (parallelabort != 1) {
		    
	    currentregion = a;
	    c = normmin;                                                 
	    d = normmax;
		    
	    fc = normgoldensearch_func(c, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, b, emibymapfreq[j][g], ul, upperlimits);
	    fd = normgoldensearch_func(d, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, b, emibymapfreq[j][g], ul, upperlimits);
		    
	    goldensearch(normgoldensearch_func, &c, &fc, &d, &fd, tolerance, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, b, emibymapfreq[j][g], ul, upperlimits);
	
	    if (fc < fd) {
	      tmp_chibyon[j+(l*(resolution+1))] = (float)fc;
	      tmp_chion[j+(l*(resolution+1))] = tmp_steps;
	    }
	    else {
	      tmp_chibyon[j+(l*(resolution+1))] = (float)fd;
	      tmp_chion[j+(l*(resolution+1))] = tmp_steps;
	    }

	    if (tmp_chibyon[j+(l*(resolution+1))] < tmp_bestchi) {
		      	    
	      tmp_bestchi = tmp_chibyon[j+(l*(resolution+1))];
	      tmp_bestage = tmp_steps;
	      tmp_beststepindex = (j+(l*(resolution+1)));
	    }
	  }
	}
      }

      // If we have bad values, kick us out
      if ( (*suppress == 0) && (parallelabort == 1) ) {
	return 5;
      }

      // Set the new min and max ages
      tmp_maxage = maxage;
      tmp_minage = minage;
	  
      // Set the new min and max levels
      if ( ( tmp_bestage - ( (  ( (tmp_maxage-tmp_minage) / resolution ) /2)*pow(10, -l) ) ) <= 0.01) { // Prevent an impossible age
	minage = 0.01;
      }
      else {
	minage = tmp_bestage - ( (  ( (tmp_maxage-tmp_minage) / resolution ) /2)*pow(10, -l) );
      }
      maxage = tmp_bestage + ( (  ( (tmp_maxage-tmp_minage) / resolution ) /2)*pow(10, -l) );
    }

    // Now do the same for the off time...
    if (model == 5) {
      j=1;
      age=1e6*bestage[a]*365*86400; // Fix the t_off to the best age
      tmp_bestoffchi = 1e26;
      tmp_bestoffage = 0.0;

      minage_off = (float)minoff;
      maxage_off = (float)maxoff;

      // Find the X^2 curve for off ages
      for(l=0; l<levels; l++) {

	for(g=0; g<=offresolution; g++) {

	  // Go through each map and find the range of possible normalisation values
	  normmax = -1e36;
	  normmin = 1e36;

	  if (parallelabort != 1) {

	    for(m=0; m<imgnum; m++) {
		    
	      // Stop trying if we have bad values
	      if (parallelabort != 1) {
  
		if (g == 0) {
		  tmp_offsteps = minage_off;
		  //printf("Off Steps: %.4f\n", offsteps);
		}
		else {
		  tmp_offsteps = minage_off + ( ((maxage_off-minage_off)/offresolution) * (float)g );
		  //printf("Off Steps: %.4f\n", offsteps);
		}
		    
		// Age in seconds
		t_off=1e6*tmp_offsteps*365*86400;
	    
  
		if ( (model == 5) || (model == 6) ) {
		  emibymapfreq[j][g][m]=emiss_nkg(1,b,frequency[fluxorder[m]]);
		}
		else {
		  emibymapfreq[j][g][m]=emiss_nkg(1,b,frequency[fluxorder[m]]);
		}

		// Check is the flux value is bad and if so make sure the message is only displayed once
		if ( ( (emibymapfreq[j][g][m] <= 1e-33) || (emibymapfreq[j][g][m] >= 1e33) ) && (parallelabort!=1) ) { 
		  if (*suppress == 0) {
		    printf("\nError: Model has reached its precision limit (model flux effectively 0 or unphysically high) at an age of %.2f (on) %.2f (off) Megayears. Please run plot(ci/cioff) to determine a suitable age range and / or model paramaters (e.g. bfield)\n", steps, offsteps);
		    parallelabort = 1;
		  }
		  else {
		    // Check this isnt the very first try!
		    if ( (g == 0) && (j == 0) ) {
		      printf("\nError: Model has reached its precision limit (model flux effectively 0 or unphysically high) on the very first attempt! This cannot be suppressed. No output will made for source with index %s.\n", index);
		      *suppress = 0;
		      parallelabort = 1;
		    }
		    else {     
		      // Increment to indicate we had to suppress something (doesnt matter if we have mutiple increments as we are only interested if it is > 1)
		      (*suppress)++;
		      parallelabort = 1;
		    }
		  }
		}

		normfactor = regflux[fluxorder[m]][a]/(emibymapfreq[j][g][m]);
		normfactor = log10(normfactor);
		      		      
		if ( normfactor > normmax ) {
		  normmax = normfactor;
		}
		      
		if ( normfactor < normmin ) {
		  normmin = normfactor;
		}
	      }
	    }

	    normmax = pow(10,normmax) + (pow(10,normmax)*0.5);
	    normmin = pow(10,normmin) - (pow(10,normmin)*0.5);
	    normmax = log10(normmax);
	    normmin = log10(normmin);
		  
	    if (parallelabort != 1) {
		    
	      currentregion = a;
	      c = normmin;                                                 
	      d = normmax;
		    
	      fc = normgoldensearch_func(c, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, b, emibymapfreq[j][g], ul, upperlimits);
	      fd = normgoldensearch_func(d, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, b, emibymapfreq[j][g], ul, upperlimits);
		    
	      goldensearch(normgoldensearch_func, &c, &fc, &d, &fd, tolerance, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, b, emibymapfreq[j][g], ul, upperlimits);
	
	      if (fc < fd) {
		tmp_chibyoff[g+(l*(resolution+1))] = (float)fc;
		tmp_chioff[g+(l*(resolution+1))] = tmp_offsteps;
	      }
	      else {
		tmp_chibyoff[g+(l*(resolution+1))] = (float)fd;
		tmp_chioff[g+(l*(resolution+1))] = tmp_offsteps;
	      }

	      if (tmp_chibyoff[g+(l*(resolution+1))] < tmp_bestoffchi) {
		      	    
		tmp_bestoffchi = tmp_chibyoff[g+(l*(offresolution+1))];
		tmp_bestoffage = tmp_offsteps;
		tmp_bestoffindex = (g+(l*(resolution+1)));
	      }
	    }
	  }
	}

	// Set the new min and max ages
	tmp_maxage_off = maxage_off;
	tmp_minage_off = minage_off;
	
	// Set the new min and max levels for t_off
	if ( ( tmp_bestoffage - ( (  ( (tmp_maxage_off-tmp_minage_off) / resolution ) /2)*pow(10, -l) ) ) <= 0) { // Prevent an impossible age
	  minage_off = 0.0;
	}
	else {
	  minage_off = tmp_bestoffage - ( (  ( (tmp_maxage_off-tmp_minage_off) / resolution ) /2)*pow(10, -l) );
	}
	maxage_off = tmp_bestoffage + ( (  ( (tmp_maxage_off-tmp_minage_off) / resolution ) /2)*pow(10, -l) );
      }
    }

    // Order the arrays by age for curve fitting later
    for(s=0; s<(resolution+1)*levels; s++) { // /get the on error for CI and CI off
      for(u=s; u<(resolution+1)*levels; u++) {
	if(tmp_chion[s] > tmp_chion[u]) {

	  tmp_sort=tmp_chion[s];
	  tmp_chion[s]=tmp_chion[u];
	  tmp_chion[u]=tmp_sort;

	  tmp_sort=tmp_chibyon[s];
	  tmp_chibyon[s]=tmp_chibyon[u];
	  tmp_chibyon[u]=tmp_sort;

	  if (u == tmp_beststepindex) {
	    tmp_beststepindex=s;
	  }
	}
      }
    }
    
    if (model == 5) { // Get the off error for CI off
      for(s=0; s<(resolution+1)*levels; s++) {
	for(u=s; u<(resolution+1)*levels; u++) {
	  if(tmp_chioff[s] > tmp_chioff[u]) {

	    tmp_sort=tmp_chioff[s];
	    tmp_chioff[s]=tmp_chioff[u];
	    tmp_chioff[u]=tmp_sort;

	    tmp_sort=tmp_chibyoff[s];
	    tmp_chibyoff[s]=tmp_chibyoff[u];
	    tmp_chibyoff[u]=tmp_sort;

	    if (u == tmp_bestoffindex) {
	      tmp_bestoffindex=s;
	    }
	  }
	}
      }
    }
    

    // Currently only linear interpolation works, but this allows for the option of adding a cspline later.
    int interptype = 0;

    if (interptype == 0) { // Linearlly interpolate between the clostest points to a delta chi^2 of 1
      for (cub_cnt=tmp_beststepindex; cub_cnt>=0; cub_cnt--) {

	if ( (tmp_chibyon[cub_cnt]-tmp_chibyon[tmp_beststepindex] >= 1.000) && (tmp_chion[cub_cnt] < tmp_chion[tmp_beststepindex]) ) {

	  tmp_deltachiof1 = tmp_chibyon[tmp_beststepindex]+1;
	  slope = (tmp_chibyon[cub_cnt+1]-tmp_chibyon[cub_cnt]) / (tmp_chion[cub_cnt+1]-tmp_chion[cub_cnt]);
	  intercept = tmp_chibyon[cub_cnt+1] - (slope * tmp_chion[cub_cnt+1]);
	  onerrorsminus[a] = fabs( tmp_chion[tmp_beststepindex] - ((tmp_deltachiof1-intercept) / slope) );
	  break;
	}

	// Else if we have hit 0 and still not got a delta X^2 of 1...
	else if ( (cub_cnt == 0) && (tmp_chibyon[cub_cnt]-tmp_chibyon[beststepindex] < 1.000) ) {
	  onerrorsminus[a] = tmp_chion[tmp_beststepindex];
	}


      }

      for (cub_cnt=tmp_beststepindex; cub_cnt<(resolution+1)*levels; cub_cnt++) {

	if ( (tmp_chibyon[cub_cnt]-tmp_chibyon[tmp_beststepindex] >= 1.000)  && (tmp_chion[cub_cnt] > tmp_chion[tmp_beststepindex]) ) {

	  tmp_deltachiof1 = tmp_chibyon[tmp_beststepindex]+1;
	  slope = (tmp_chibyon[cub_cnt-1]-tmp_chibyon[cub_cnt]) / (tmp_chion[cub_cnt-1]-tmp_chion[cub_cnt]);
	  intercept = tmp_chibyon[cub_cnt-1] - (slope * tmp_chion[cub_cnt-1]);
	  onerrorsplus[a] = fabs( tmp_chion[tmp_beststepindex] - ((tmp_deltachiof1-intercept) / slope) );

	  break;
	}

	// Else if we have hit the end of our data points and still not got a delta X^2 of 1...
	else if ( (cub_cnt == ((resolution+1)*levels)) && (tmp_chibyon[cub_cnt]-tmp_chibyon[tmp_beststepindex] < 1.000) ) {
	  onerrorsplus[a] = 0.000000;
	  printf("*** Warning: Could not find an upper error for the on component of region %s! Error has been set to 0. Consider expanding the age search range or changing the model parameters. ***\n\n", index);
	}

      }
    

      if (model == 5) {
	for (cub_cnt=tmp_bestoffindex; cub_cnt>=0; cub_cnt--) {

	  if ( (tmp_chibyoff[cub_cnt]-tmp_chibyoff[tmp_bestoffindex] >= 1.000) && (tmp_chioff[cub_cnt] < tmp_chioff[tmp_bestoffindex]) ) {

	    tmp_deltachiof1 = tmp_chibyoff[tmp_bestoffindex]+1;
	    slope = (tmp_chibyoff[cub_cnt+1]-tmp_chibyoff[cub_cnt]) / (tmp_chioff[cub_cnt+1]-tmp_chioff[cub_cnt]);
	    intercept = tmp_chibyoff[cub_cnt+1] - (slope * tmp_chioff[cub_cnt+1]);
	    offerrorsminus[a] = fabs( tmp_chioff[tmp_bestoffindex] - ((tmp_deltachiof1-intercept) / slope) );

	    break;
	  }
	  // Else if we have hit 0 and still not got a delta X^2 of 1...
	  else if ( (cub_cnt == 0) && (tmp_chibyoff[cub_cnt]-tmp_chibyoff[tmp_bestoffindex] < 1.000) ) {
	    offerrorsminus[a] = tmp_chioff[tmp_bestoffindex];
	  }
	}

	for (cub_cnt=tmp_bestoffindex; cub_cnt<(resolution+1)*levels; cub_cnt++) {
	  if ( (tmp_chibyoff[cub_cnt]-tmp_chibyoff[tmp_bestoffindex] >= 1.000)  && (tmp_chioff[cub_cnt] > tmp_chioff[tmp_bestoffindex]) ) {

	    if (tmp_chibyoff[cub_cnt] > 1e32) {
	      printf("\n*** Warning: The off component of the chi-squared curve appears to contain a discontinuity! This is normally caused by upper limits. Setting the positive error to 0. ***\n\n");

	      offerrorsplus[a] = 0.0;
	    }

	    tmp_deltachiof1 = tmp_chibyoff[tmp_bestoffindex]+1;
	    slope = (tmp_chibyoff[cub_cnt-1]-tmp_chibyoff[cub_cnt]) / (tmp_chioff[cub_cnt-1]-tmp_chioff[cub_cnt]);
	    intercept = tmp_chibyoff[cub_cnt-1] - (slope * tmp_chioff[cub_cnt-1]);
	    offerrorsplus[a] = fabs( tmp_chioff[tmp_bestoffindex] - ((tmp_deltachiof1-intercept) / slope) );

	    break;
	  }
	  // Else if we have hit the end of our data points and still not got a delta X^2 of 1...
	  else if ( (cub_cnt == ((resolution+1)*levels)) && (tmp_chibyoff[cub_cnt]-tmp_chibyoff[tmp_bestoffindex] < 1.000) ) {
	    offerrorsplus[a] = 0.000000;
	    printf("*** Warning: Could not find an upper error for the off component of region %s! Error has been set to 0. Consider expanding the age search range or changing the model parameters. ***\n", index);
	  }

	}
      }

    }

    /*
      free(tmp_chibyon);
      free(tmp_chion);
      free(tmp_chibyoff);
      free(tmp_chioff);*/
  }

  sumchisquared = avgchisquared;
  avgchisquared /= numregions;
  avgnormalisation /= numregions;

  // Get the model fluxes and spectral indices
  regcounter = 0;
  hundreds = 0;

  // Check min and max frequencies aren't the same
  if (fluxorder[0] != fluxorder[imgnum-1]) {

    // printf("Calculating model fluxes and spectral indices\n Currently processing region: 1 of %d ...\n", numregions);

    steps = 0;

    // Get the predicted flux for each age at myears*resolution intervals
#pragma omp parallel for private(i, j, nu, emi, freqsteps) schedule(dynamic)
    for(i=1; i<=numregions; i++) {
 
      //Best model age fit
      age=1e6*bestage[i]*365*86400;
      t_off=1e6*bestoff[i]*365*86400;


      // No need for normalisation here as all alpha and curvature are relative
      //#pragma omp parallel for private(j, nu, emi, freqsteps) firstprivate(i) schedule(dynamic)
      for (j=0; j<=modelres; j++) {

	if ( (extrapolatemodel == 1) && (extrapolationfrequency > frequency[fluxorder[imgnum-1]]) ) {
	freqsteps = ( (extrapolationfrequency - frequency[fluxorder[0]]) / modelres ) * j;
	}
	else if ( (extrapolatemodel == 1) && (extrapolationfrequency < frequency[fluxorder[0]]) ){
	freqsteps = ( (frequency[fluxorder[imgnum-1]] - extrapolationfrequency) / modelres ) * j;
	}
	else {
	freqsteps = ( (frequency[fluxorder[imgnum-1]] - frequency[fluxorder[0]]) / modelres ) * j;
	}

	// Region Flux
	if ( (extrapolatemodel == 1) && (extrapolationfrequency < frequency[fluxorder[0]]) ){
	nu=extrapolationfrequency + freqsteps;
	}
	else { // This works for both normal fitting and extrapolation to higher frequencies
	nu=frequency[fluxorder[0]] + freqsteps;
	}

	if ( (model == 5) || (model == 6) ) {
	  emi=emiss_nkg(1,b,nu);
	}
	else {
	  emi=emiss_nkg(1,b,nu);
	}
      
	modelflux[passingreg][j] = emi;
	modelflux[passingreg][j] *= bestnorm[i];
      }
    
      // Find alpha for the best model fit of each region  
      modelalpha[i] = -log(modelflux[passingreg][0]/modelflux[passingreg][modelres]) / log(frequency[fluxorder[0]]/frequency[fluxorder[imgnum-1]]);

      regcounter++;

      if (regcounter >= 100) {

	hundreds++;
	printf("%d00...\n", hundreds);
	regcounter = 0;
      }
    }    
  }

  else {
    fprintf(stderr,"\nError: Cannot calculate spectral index. Are both maps the same frequency???\nExiting...\n\n");
    return(500);
  }

  // Print out the confidence levels
  if ( (export == 0) && (suppresscdf != 1) ) {
    printf("\nModel can be rejected with 68 per cent confidence at X^2 > %.2f (Reduced: %.2f)\n", siglevel68, siglevel68 / dof);
    printf("Model can be rejected with 90 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel90, siglevel90 / dof);
    printf("Model can be rejected with 95 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel95, siglevel95 / dof);
    printf("Model can be rejected with 99 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel99, siglevel99 / dof);
    printf("Model can be rejected with 99.5 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel995, siglevel995 / dof);
    printf("Model can be rejected with 99.9 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel999, siglevel999 / dof);
    printf("Model can be rejected with 99.99 per cent confidence X^2 > %.2f (Reduced: %.2f)\n\n", siglevel9999, siglevel9999 / dof);
  }

  // Find the minimum and maximum ages

  maxfitage = -1e22;
  maxfitage_off = -1e22;
  oldestregion = -1;

  for(a=1; a<=numregions; a++) {
    if ((bestage[a]+bestoff[a]) > (maxfitage+maxfitage_off)) {
      maxfitage = bestage[a];
      maxfitage_off = bestoff[a];
      oldestregion = a;
    }
  }


  // (ECF*b*pow(1 / (((4*THOMSON/(6*M_EL*V_C*MU_0))*(ageb*ageb)) * (bestage[a]+bestoffage[a])),2);)/(2*PI*M_EL) 

  // Print out the results and statistics

  // Print out the confidence levels
  if (export == 0) {
    printf("Sum of X^2 %.2f \n", sumchisquared);
    printf("Average X^2 %.2f \n", avgchisquared);
    printf("Average normalisation %.2e \n", avgnormalisation);
    printf("Average Reduced X^2 %.2f for %d degrees of freedom\n", ( avgchisquared / (imgnum - ul - freeparams)), (imgnum - ul - freeparams) );
  }

  a = 1; // This needs changing for when we have more than 1 source / region
  if ((numregions == 1) && (model ==5)) {

    // Set the break to variables so we can pass them back
    breakon[passingreg] = (ECF*b*pow(1 / (((4*THOMSON/(6*M_EL*V_C*MU_0))*(ageb*ageb)) * (1e6*365*86400*(bestage[a]+bestoff[a]))),2))/(2*PI*M_EL);
    breakoff[passingreg] = (ECF*b*pow(1 / (((4*THOMSON/(6*M_EL*V_C*MU_0))*(ageb*ageb)) * (bestoff[a]*1e6*365*86400)),2))/(2*PI*M_EL);

    if (export == 0) {
      printf("Break frequency is %.2e Hz with an off component break at %.2e Hz\n", breakon[passingreg], breakoff[passingreg]);
    }
  }
  else if ((numregions == 1) && (model ==6)) {
    breakon[passingreg] = (ECF*b*pow(1 / (((4*THOMSON/(6*M_EL*V_C*MU_0))*(ageb*ageb)) * (1e6*365*86400*(bestage[a]+bestoff[a]))),2))/(2*PI*M_EL);
    breakoff[passingreg] = (-1); // Indicate not used in model to the caller

    if (export == 0) {
      printf("Break frequency is %.2e Hz\n", breakon[passingreg]);
    }
  }

  if (export == 0) {
    printf("Total Age: %.2f +%.2f -%.2f | On: %.2f +%.2f -%.2f | Off: %.2f +%.2f -%.2f Megayears (Region %d)\n\n", maxfitage+maxfitage_off, sqrt(pow(onerrorsplus[oldestregion],2)+pow(offerrorsplus[oldestregion],2)), sqrt(pow(onerrorsminus[oldestregion],2)+pow(offerrorsminus[oldestregion],2)), maxfitage, onerrorsplus[oldestregion], onerrorsminus[oldestregion], maxfitage_off, offerrorsplus[oldestregion], offerrorsminus[oldestregion], oldestregion);
  }

  if (suppresscdf != 1) {

    // Find the number of rejected regions
    nonrejectbin = 0;
    rejectbin = 0;

    for (i=0; i<=2; i++) {
      nonrejectbin += medianbin[i];
    }
    for (i=3; i<=4; i++) {
      rejectbin += medianbin[i];
    }

    if (export == 0) {

      if (numregions > 1) {
	printf("Bin Statistics:\n");
	printf("===========================================================\n\n");
	printf(" Index  |   Confidence Level   |   Regions   |    Fraction   \n\n");
	printf("   0    |    < 68 per cent     |     %d      |      %.2f     \n", medianbin[0], (float)medianbin[0]/(float)numregions);
	printf("   1    |   68 - 90 per cent   |     %d      |      %.2f     \n", medianbin[1], (float)medianbin[1]/(float)numregions);
	printf("   2    |   90 - 95 per cent   |     %d      |      %.2f     \n", medianbin[2], (float)medianbin[2]/(float)numregions);
	printf("   3    |   95 - 99 per cent   |     %d      |      %.2f     \n", medianbin[3], (float)medianbin[3]/(float)numregions);
	printf("   4    |    > 99 per cent     |     %d      |      %.2f     \n", medianbin[4], (float)medianbin[4]/(float)numregions);
	printf("\n===========================================================\n\n");
      }
      else {

	char confbuff[16];
	
	if (medianbin[0] > 0) {
	      sprintf(confbuff,"< 68");
	}
	else if (medianbin[1] > 0) {
	      sprintf(confbuff,"68 - 90");
	}
	else if (medianbin[2] > 0) {
	      sprintf(confbuff,"90 - 95");
	}
	else if (medianbin[3] > 0) {
	      sprintf(confbuff,"95 - 99");
	}
	else if (medianbin[4] > 0) {
	      sprintf(confbuff,"> 99");
	}
	else  {
	      sprintf(confbuff,"*** Unknown ***");
	}

	printf("This source falls within the %s per cent confidence bin\n\n", confbuff);
	
      }
    }

    // Find the median of the binned chi-squared values and output the result 
    if (nonrejectbin > rejectbin) {

      median = 0;
      for (i=1; i<=2; i++) {
	if (medianbin[i] > medianbin[median]) {
	  median = i;
	}
	if ( (medianbin[i] == medianbin[median]) && (i != median ) && (export == 0) ) {
	  printf("\n *** Warning: Median bin %d and %d contain the same number of regions! *** \n\n", i, median);
	}
      }

      if (export == 0) {
	if (median == 0) {
	  confpercentage = 68.0;
	}
	else if (median == 1) {
	  confpercentage = 90.0;
	}
	else if (median == 2) {
	  confpercentage = 95.0;
	}
  
	printf("Model cannot be rejected at the %.0f per cent significance level based on median binning\n\n", confpercentage);
      }
    }
    else {

      median = 3;
      for (i=4; i<=4; i++) {
	if (medianbin[i] > medianbin[median]) {
	  median = i;
	}
	if ( (medianbin[i] == medianbin[median]) && (i != median ) && (export == 0) ) {
	  printf("\n *** Warning: Median bin %d and %d contain the same number of regions! *** \n\n", i, median);
	}
      }

      if (median == 3) {
	confpercentage = 95.0;
      }
      else if (median == 4) {
	confpercentage = 99.0;
      }

      if (export == 0) {
	printf("Model is rejected at the %.0f per cent significance level based on median binning\n\n", confpercentage);
      }
    }
  }

  // Find if the model can be rejected or not using the average chi-squared

  if ( (printreject == 1) && (export == 0) && (suppresscdf !=1) ) {
    if (avgchisquared > siglevel9999) {
      printf("Model can be rejected at the 99.99 per cent confidence level!\n\n");
    }
    else if (avgchisquared > siglevel999) {
      printf("Model can be rejected at the 99.9 per cent confidence level!\n\n");
    }
    else if (avgchisquared > siglevel995) {
      printf("Model can be rejected at the 99.5 per cent confidence level!\n\n");
    }
    else if (avgchisquared > siglevel99) {
      printf("Model can be rejected at the 99.0 per cent confidence level!)\n\n");
    }
    else if (avgchisquared > siglevel95) {
      printf("Model can be rejected at the 95.0 per cent confidence level!\n\n");
    }
    else if (avgchisquared > siglevel90) {
      printf("Model can be rejected at the 90.0 per cent confidence level!\n\n");
    }
    /*else if (avgchisquared > siglevel68) {
      printf("Model can be rejected at the 68.0 per cent confidence level!\n\n");
      }*/
    else {
      printf("Model is a good fit and cannot be rejected above a 90 per cent confidence!\n\n");
    }
  }

/*
  for (i=0; i<=resolution; i++){
    free(emibymapfreq[i]);
    }
    free(emibymapfreq);*/
  free(medianbin);

  return 0;
}





int plotcimodel(float minf, float maxf, double inject, double fieldstrength, int model, double usr_gmin, double usr_gmax, int modelmyears, int minmodelmyears, char *output, int titles, float redshift, int minoff, int maxoff, int skip, int varyoffage) {
    
  // Set up and plot a set of aged spectra
    
#define STEPSKG 100
    
  // Declaring PGPLOT prototypes
  int cpgopen(const char *device);
  void cpgclos(void);
  void cpgenv(float xmin, float xmax, float ymin, float ymax, int just, int axis);
  void cpglab(const char *xlbl, const char *ylbl, const char *toplbl);
  void cpgdraw(float x, float y);
  void cpgmove(float x, float y);
  void cpgsci(int ci);
  void cpgscr(int ci, float cr, float cg, float cb);
    
  int i,j;
  double b, emi, modnormfact, lowagenorm;
  float delt, yaxismax, yaxismin, nu, bcmb, **plotarray;
  char modeltitle[128];
  int cycleoff = 0; // Sets whether we cycle the on time or off time


  if (model == 6) {
    t_off = 0.0;
    age=1e6*minmodelmyears*365*86400; // For axis and normalization factor later on

    if (minmodelmyears == 0) {
      age=1e6*0.01*365*86400; // For axis and normalization factor later on
    }

    plotarray = (float **)calloc(((modelmyears-minmodelmyears)/skip)+1, sizeof(float *));
    for (i=0; i<=((modelmyears-minmodelmyears)/skip); i++) {
      plotarray[i] = (float *)calloc(STEPSKG+1, sizeof(float));
    }
  }
  //else if ( ((maxoff-minoff) == 0) && (model == 5) ) {
  else if ( (varyoffage == 0) && (model == 5) ) {
    t_off = maxoff*(1e6*365*86400);
    age=1e6*minmodelmyears*365*86400;
    if (minmodelmyears == 0) {
      age=1e6*0.01*365*86400; // For axis and normalization factor later on
    }
    cycleoff = 0;
    plotarray = (float **)calloc(((modelmyears-minmodelmyears)/skip)+1, sizeof(float *));
    for (i=0; i<=((modelmyears-minmodelmyears)/skip); i++) {
      plotarray[i] = (float *)calloc(STEPSKG+1, sizeof(float));
    }
  }
  else if ( (varyoffage == 1) && (model == 5) ) {
    cycleoff = 1;
    t_off = minoff*(1e6*365*86400); // Set this now so we can find a sensible y axis
    age=1e6*modelmyears*365*86400; // As we fix the on age to modelmyears 
    if (modelmyears == 0) {
      age=1e6*0.01*365*86400; // For axis and normalization factor later on
    }
    plotarray = (float **)calloc(((maxoff-minoff)/skip)+1, sizeof(float *));
    for (i=0; i<=((maxoff-minoff)/skip); i++) {
      plotarray[i] = (float *)calloc(STEPSKG+1, sizeof(float));
    }

  }
  else {
    printf("\nError: On / off combination is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");
    return(404);
  }

#pragma omp parallel
  {
    w1=gsl_integration_workspace_alloc(GSL_WSIZE);
    w2=gsl_integration_workspace_alloc(GSL_WSIZE);
    w3=gsl_integration_workspace_alloc(GSL_WSIZE);
    w4=gsl_integration_workspace_alloc(GSL_WSIZE);
  }
    
  gsl_set_error_handler_off(); // to avoid abort
  gsl_epsilon=1e-3;
  gsl_quiet=1;
    
  //fx=calloc(xls,sizeof(double));
  //makefx_kg();

  fy=calloc(xls,sizeof(double));
  makefy_kg();
    
  printf("Ready to go!\n");
    
  gmin=usr_gmin;
  gmax=usr_gmax;


  // Default of 2.2 gives a spetcral index of 0.6
  power=(2*inject)+1; //Now a user input
    
  // Set which model energy distribution to use (always KG at the moment, but allows for additional similar models later).
  if (model == 5) { // CI off model
    ne=ne_age_kg;
    age_gmax=age_gmax_kg;
  }
  else if (model == 6) { // CI model
    ne=ne_age_kg;
    age_gmax=age_gmax_kg;
  }
  else {
    printf(" *** Warning: Unknown Model. Defaulting to CI off ***\n");
    age_gmax=age_gmax_kg;
  }
    
  GMIN=gmin;
  GMAX=gmax;
    
  b=fieldstrength;
  bcmb=(0.318 * (pow(1+redshift,2))) * 1e-9;
  ageb=sqrt(pow(b,2) + pow(bcmb,2));
  //ageb=b;
    
  printf("Magnetic field strength: %.3e T\n", b);
  printf("B_CMB: %.3e T\n", bcmb);
  printf("Ageing field strength: %.3e T\n", ageb);
  
    
  cpgopen(output);  
    
  // Swap to background and forground colours to something sensible
  cpgscr(0, 1, 1, 1);
  cpgscr(1, 0, 0, 0);

  // Get the flux of the lowest age to normalize everything to later for an easier to understand plot. Ages and off times are set further up the start as they are dependent on the plot type
  nu=minf;
  lowagenorm = emiss_nkg(1,b,nu);
    
  // Calculate a sensible Y axis scale
    
  if ( (model == 5) || (model == 6) ) {
    emi = emiss_nkg(1,b,nu);
  }
  else {
    //sprintf(" *** Warning: Unknown Model. Defaulting to CI off ***\n");
    emi=emiss_nkg(1,b,nu);
  }
    
  if(log10(emi) >= 1) {
    yaxismax = log10(emi) * 1.05;
  }
  else {
    yaxismax = log10(emi) * 0.95;
  }

  //printf("maxemi: %.4e yaxismax: %.4e ", emi, yaxismax);
    
  nu=maxf;
  age=1e6*modelmyears*365*86400;

  if (modelmyears == 0) {
    age=1e6*0.01*365*86400; // For axis and normalization factor later on
  }
    
  if ( (model == 5) || (model == 6) ) {
    emi=emiss_nkg(1,b,nu);
  }
  else {
    //sprintf(" *** Warning: Unknown Model. Defaulting to CI off ***\n");
    emi=emiss_nkg(1,b,nu);
  }

  modnormfact=emiss_nkg(1,b,minf);
  emi*=(lowagenorm/modnormfact);

  if ( emi <= 1e-33 ) {
    yaxismin = yaxismax - 12;
  }
  else{
    if(log10(emi) >= 1) {
      yaxismin = log10(emi) * 0.95;
    }
    else {
      yaxismin = log10(emi) * 1.05;
    }
  }

  cpgenv(log10(minf),log10(maxf), yaxismin, yaxismax, 0, 30);
    
  // Set the title if required
    
  if (titles == 0) { sprintf(modeltitle," ");  }
  else {
    if (model == 5) { // CI off model
      if (varyoffage == 0) {
	if (minmodelmyears == 0) {
	  sprintf(modeltitle,"CI_Off model between 0.01 and %d Myrs (ON) and %d Myrs (Off) (arbitrary normalisation)", modelmyears, maxoff);
	}
	else {
	  sprintf(modeltitle,"CI_Off model between %d and %d Myrs (ON) and %d Myrs (Off) (arbitrary normalisation)", minmodelmyears, modelmyears, maxoff);
	}
      }
      else{
	if (modelmyears == 0) {
	  sprintf(modeltitle,"CI_Off model between %d and %d Myrs (Off) and 0.01 Myrs (On) (arbitrary normalisation)", minoff, maxoff);
	}
	else {
	  sprintf(modeltitle,"CI_Off model between %d and %d Myrs (Off) and %d Myrs (On) (arbitrary normalisation)", minoff, maxoff, modelmyears);
	}
      }
    }
    else if (model == 6) { // CI model
      if (minmodelmyears == 0) {
	sprintf(modeltitle,"CI model between 0.01 and %d Myrs (arbitrary normalisation)", modelmyears);
      }
      else {
	sprintf(modeltitle,"CI model between %d and %d Myrs (arbitrary normalisation)", minmodelmyears, modelmyears);
      }
    }
    else {
      sprintf(modeltitle,"Unknown Model Between %d and %d Myrs (arbitrary normalisation)", minmodelmyears, modelmyears);
    }
  }
    

  cpglab("Frequency / Hz", "Flux / Jy", modeltitle);
    
  delt=log10(maxf/minf)/(float)STEPSKG;
    
  //firstprivate(plotarray)
    
  //#pragma omp parallel for private(j, i, emi, nu) schedule(dynamic)
  if (cycleoff != 1) {
    for(j=0; j<=((modelmyears-minmodelmyears)/skip); j++) {
        
      age=1e6*((j*skip)+minmodelmyears)*365*86400;

      if ( (minmodelmyears == 0) && (j == 0)) {
	age=1e6*0.01*365*86400; // For axis and normalization factor later on
      }

      if (j != 0) {
	modnormfact=emiss_nkg(1,b,minf);
	modnormfact=lowagenorm/modnormfact;
      }
      else {
	modnormfact = 1;
      }

      //#pragma omp parallel for private(i, nu, emi) schedule(dynamic)
      for(i=0; i<=STEPSKG; i++) {
	// x-axis frequency for each step
	nu=minf*pow(10,delt*i);
            
	// y-axis emission for each step
	if ( (model == 5) || (model == 6) ) {
	  plotarray[j][i] = emiss_nkg(1,b,nu)*modnormfact;
	  //if ((j==1) || (j==10) ) {
	  //  printf("plotarray[%d][%d] = %.4e\n", j, i, plotarray[j][i]);
	  //}
	}
	else {
	  printf(" *** Warning: Unknown Model. Defaulting to CI off ***\n");
	  plotarray[j][i] = emiss_nkg(1,b,nu)*modnormfact;
	}
      }
    }

    // Easy way to avoid the race condition on colours
    //#pragma omp barrier
    
      for(j=0; j<=((modelmyears-minmodelmyears)/skip); j++) {
	cpgsci(j+2); //The +2 avoid the background and gives the standard 0 to 10 start and end distinctive colours.
        
	for(i=0; i<=STEPSKG; i++) {
	  nu=minf*pow(10,delt*i);
	  (i?cpgdraw:cpgmove)(log10(nu),log10(plotarray[j][i]));
	}
      }
  }
  else {

    age=1e6*modelmyears*365*86400;

    if (modelmyears == 0) {
      age=1e6*0.01*365*86400; // For axis and normalization factor later on
    }

    for(j=0; j<=((maxoff-minoff)/skip); j++) {
        
      t_off=1e6*((j*skip)+minoff)*365*86400;

      if (j != 0) {
	modnormfact=emiss_nkg(1,b,minf);
	modnormfact=lowagenorm/modnormfact;
      }
      else {
	modnormfact = 1;
      }

      //#pragma omp parallel for private(i, nu, emi) schedule(dynamic)
      for(i=0; i<=STEPSKG; i++) {
	// x-axis frequency for each step
	nu=minf*pow(10,delt*i);
            
	// y-axis emission for each step
	if ( (model == 5) || (model == 6) ) {
	  plotarray[j][i] = emiss_nkg(1,b,nu)*modnormfact;
	}
	else {
	  printf(" *** Warning: Unknown Model. Defaulting to CI off Model: %d ***\n", model);
	  plotarray[j][i] = emiss_nkg(1,b,nu)*modnormfact;
	}
      }
    }

    // Easy way to avoid the race condition on colours
    //#pragma omp barrier
    
    for(j=0; j<=((maxoff-minoff)/skip); j++) {
      cpgsci(j+2); //The +2 avoid the background and gives the standard 0 to 10 start and end distinctive colours.
        
      for(i=0; i<=STEPSKG; i++) {
	nu=minf*pow(10,delt*i);
	(i?cpgdraw:cpgmove)(log10(nu),log10(plotarray[j][i]));
      }
    }

  }
    
  cpgclos();
    
  return 0;
}





int cimodeldata(float minf, float maxf, double inject, double fieldstrength, int model, double usr_gmin, double usr_gmax, int minmodelmyears, int modelmyears, char *filename, float redshift, int minoff, int maxoff, int skip, int varyoffage, int exactage, float exactmyears, float exactoff, int dataintervals) {


  // If the min in greater than the max, return an error
  if ( (minmodelmyears > modelmyears) && (exactage != 1) ) {
    printf("\nError: The minimum on time (%d Myr) must be less than the maximum on time (%d Myr). Please use the minmodelmyears and/or the modelmyears commands to set a new range.\n\n", minmodelmyears, modelmyears);
    return 6;
  }

  if (model == 5) {
    // If the min off in greater than the maxoff and we are using the CI off model, return an error
    if ( (minoff > maxoff) && (exactage != 1) ) {
      printf("\nError: The minimum off time (%d Myr) must be less than the maximum off time (%d Myr). Please use the modelminoff and/or the modelmaxoff commands to set a new range.\n\n", minoff, maxoff);
      return 6;
    }
  }

  // Declaring PGPLOT prototypes
  int cpgopen(const char *device);
  void cpgclos(void);
  void cpgenv(float xmin, float xmax, float ymin, float ymax, int just, int axis);
  void cpglab(const char *xlbl, const char *ylbl, const char *toplbl);
  void cpgdraw(float x, float y);
  void cpgmove(float x, float y);
  void cpgsci(int ci);
  void cpgscr(int ci, float cr, float cg, float cb);
  
  int i,j;
  double b, emi, modnormfact, lowagenorm;
  float delt, nu, bcmb, **plotarray;
  int cycleoff = 0; // Sets whether we cycle the on time or off time


  if (model == 6) {
    if (exactage != 1 ) {

      plotarray = (float **)calloc(((modelmyears-minmodelmyears)/skip)+1, sizeof(float *));
      
      for (i=0; i<=((modelmyears-minmodelmyears)/skip); i++) {
	plotarray[i] = (float *)calloc(dataintervals+1, sizeof(float));
      }
    }
    else {
      plotarray = (float **)calloc(1, sizeof(float *)); 
      plotarray[0] = (float *)calloc(dataintervals+1, sizeof(float)); 
    }
  }
  else if ( (varyoffage == 0) && (model == 5) ) {

    cycleoff = 0;

    if (exactage != 1 ) {
      plotarray = (float **)calloc(((modelmyears-minmodelmyears)/skip)+1, sizeof(float *));
      for (i=0; i<=((modelmyears-minmodelmyears)/skip); i++) {
	plotarray[i] = (float *)calloc(dataintervals+1, sizeof(float));
      }
    }
    else {
      plotarray = (float **)calloc(1, sizeof(float *)); 
      plotarray[0] = (float *)calloc(dataintervals+1, sizeof(float)); 
    }
  }
  else if ( (varyoffage == 1) && (model == 5) ) {

    cycleoff = 1;

    if (exactage != 1 ) {
      plotarray = (float **)calloc(((maxoff-minoff)/skip)+1, sizeof(float *));
      for (i=0; i<=((maxoff-minoff)/skip); i++) {
	plotarray[i] = (float *)calloc(dataintervals+1, sizeof(float));
      }
    }
    else {
      plotarray = (float **)calloc(1, sizeof(float *)); 
      plotarray[0] = (float *)calloc(dataintervals+1, sizeof(float)); 
    }
  }
  else {
    printf("\nError: On / off combination is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");
    return(404);
  }

#pragma omp parallel
  {
    w1=gsl_integration_workspace_alloc(GSL_WSIZE);
    w2=gsl_integration_workspace_alloc(GSL_WSIZE);
    w3=gsl_integration_workspace_alloc(GSL_WSIZE);
    w4=gsl_integration_workspace_alloc(GSL_WSIZE);
  }
    
  gsl_set_error_handler_off(); // to avoid abort
  gsl_epsilon=1e-3;
  gsl_quiet=1;

  fy=calloc(xls,sizeof(double));
  makefy_kg();
    
  printf("Ready to go!\n");
    
  gmin=usr_gmin;
  gmax=usr_gmax;


  // Default of 2.2 gives a spetcral index of 0.6
  power=(2*inject)+1; //Now a user input
    
  // Set which model energy distribution to use (always KG at the moment, but allows for additional similar models later).
  if (model == 5) { // CI off model
    ne=ne_age_kg;
    age_gmax=age_gmax_kg;
  }
  else if (model == 6) { // CI model
    ne=ne_age_kg;
    age_gmax=age_gmax_kg;
  }
  else {
    printf(" *** Warning: Unknown model. Defaulting to CI off ***\n");
    age_gmax=age_gmax_kg;
  }
    
  GMIN=gmin;
  GMAX=gmax;
    
  b=fieldstrength;
  bcmb=(0.318 * (pow(1+redshift,2))) * 1e-9;
  ageb=sqrt(pow(b,2) + pow(bcmb,2));
    
  printf("Magnetic field strength: %.3e T\n", b);
  printf("B_CMB: %.3e T\n", bcmb);
  printf("Ageing field strength: %.3e T\n", ageb);

  // Get the flux of the lowest age to normalize everything to later for an easier to understand plot. Ages and off times are set further up the start as they are dependent on the plt type
  nu=minf;

  if (exactage != 1 ) {
    if (model == 6) {
      t_off = 0.0;
      age=1e6*minmodelmyears*365*86400; // For axis and normalization factor later on

      if (minmodelmyears == 0) {
	age=1e6*0.01*365*86400; // For axis and normalization factor later on
      }
    }
    else if ( (varyoffage == 0) && (model == 5) ) {
      t_off = maxoff*(1e6*365*86400);
      age=1e6*minmodelmyears*365*86400;

      if (minmodelmyears == 0) {
	age=1e6*0.01*365*86400;
      }
    }

    else if ( (varyoffage == 1) && (model == 5) ) {
      t_off = minoff*(1e6*365*86400);
      age=1e6*modelmyears*365*86400;

      if (modelmyears == 0) {
	age=1e6*0.01*365*86400;
      }
    }
    else {
      t_off = 0.0;
      age=1e6*minmodelmyears*365*86400;
    }
  }
  else { // If we are using exact times
    if (model == 6) {
      t_off = 0.0;
      age=1e6*exactmyears*365*86400; // For axis and normalization factor later on

      if (exactmyears == 0) {
	age=1e6*0.01*365*86400; // For axis and normalization factor later on
      }
    }
    else if (model == 5) {
      t_off = exactoff*(1e6*365*86400);
      age=1e6*exactmyears*365*86400;

      if (exactmyears == 0) {
	age=1e6*0.01*365*86400;
      }
    }
    else {
      t_off = 0.0;
      age=1e6*minmodelmyears*365*86400;
    }
  }

  lowagenorm = emiss_nkg(1,b,nu);


  // Now do the same but at the highest frequency and ages
  if (exactage != 1 ) {
    if (model == 6) {
      t_off = 0.0;
      age=1e6*modelmyears*365*86400; // For axis and normalization factor later on

      if (minmodelmyears == 0) {
	age=1e6*0.01*365*86400; // For axis and normalization factor later on
      }
    }
    else if ( (varyoffage == 0) && (model == 5) ) {
      t_off = maxoff*(1e6*365*86400);
      age=1e6*modelmyears*365*86400;

      if (minmodelmyears == 0) {
	age=1e6*0.01*365*86400;
      }
    }

    else if ( (varyoffage == 1) && (model == 5) ) {
      t_off = maxoff*(1e6*365*86400);
      age=1e6*modelmyears*365*86400;

      if (modelmyears == 0) {
	age=1e6*0.01*365*86400;
      }
    }
    else {
      t_off = 0.0;
      age=1e6*modelmyears*365*86400;
    }
  }
  else { // If we are using exact times
    if (model == 6) {
      t_off = 0.0;
      age=1e6*exactmyears*365*86400; // For axis and normalization factor later on

      if (exactmyears == 0) {
	age=1e6*0.01*365*86400; // For axis and normalization factor later on
      }
    }
    else if (model == 5) {
      t_off = exactoff*(1e6*365*86400);
      age=1e6*exactmyears*365*86400;

      if (exactmyears == 0) {
	age=1e6*0.01*365*86400;
      }
    }
    else {
      t_off = 0.0;
      age=1e6*minmodelmyears*365*86400;
    }
  }

  nu=maxf;
    
  if ( (model == 5) || (model == 6) ) {
    emi=emiss_nkg(1,b,nu);
  }
  else {
    //sprintf(" *** Warning: Unknown model. Defaulting to CI off ***\n");
    emi=emiss_nkg(1,b,nu);
  }

  modnormfact=emiss_nkg(1,b,minf);
  emi*=(lowagenorm/modnormfact);

  delt=log10(maxf/minf)/(float)dataintervals;
    
  if (exactage != 1) {
 
    if (cycleoff != 1) {
      for(j=0; j<=((modelmyears-minmodelmyears)/skip); j++) {
        
	age=1e6*((j*skip)+minmodelmyears)*365*86400;

	if (model == 5) {
	  t_off=1e6*maxoff*365*86400;
	}
	else {
	  t_off=0.0;
	}

	if ( (minmodelmyears == 0) && (j == 0)) {
	  age=1e6*0.01*365*86400; // For axis and normalization factor later on
	}

	if (j != 0) {
	  modnormfact=emiss_nkg(1,b,minf);
	  modnormfact=lowagenorm/modnormfact;
	}
	else {
	  modnormfact = 1;
	}

	for(i=0; i<=dataintervals; i++) {
	  // x-axis frequency for each step
	  nu=minf*pow(10,delt*i);
            
	  // y-axis emission for each step
	  if ( (model == 5) || (model == 6) ) {
	    plotarray[j][i] = emiss_nkg(1,b,nu)*modnormfact;
	  }
	  else {
	    printf(" *** Warning: Unknown model. Defaulting to CI off ***\n");
	    plotarray[j][i] = emiss_nkg(1,b,nu)*modnormfact;
	  }
	}
      }
    }
    else {

      age=1e6*modelmyears*365*86400;

      if (modelmyears == 0) {
	age=1e6*0.01*365*86400; // For axis and normalization factor later on
      }

      for(j=0; j<=((maxoff-minoff)/skip); j++) {
        
	t_off=1e6*((j*skip)+minoff)*365*86400;

	if (j != 0) {
	  modnormfact=emiss_nkg(1,b,minf);
	  modnormfact=lowagenorm/modnormfact;
	}
	else {
	  modnormfact = 1;
	}
      
	//#pragma omp parallel for private(i, nu, emi) schedule(dynamic)
	for(i=0; i<=dataintervals; i++) {
	  // x-axis frequency for each step
	  nu=minf*pow(10,delt*i);
            
	  // y-axis emission for each step
	  if ( (model == 5) || (model == 6) ) {
	    plotarray[j][i] = emiss_nkg(1,b,nu)*modnormfact;
	  }
	  else {
	    printf(" *** Warning: Unknown model. Defaulting to model: %d ***\n", model);
	    plotarray[j][i] = emiss_nkg(1,b,nu)*modnormfact;
	  }
	}
      }
    }
  }
  else {
    if (model == 5) {
      age=1e6*exactmyears*365*86400;
      t_off=1e6*exactoff*365*86400;

      modnormfact=emiss_nkg(1,b,minf);
      modnormfact=lowagenorm/modnormfact;
    }
    else if (model == 6) {
      age=1e6*exactmyears*365*86400;
      t_off=0.0;

      modnormfact=emiss_nkg(1,b,minf);
      modnormfact=lowagenorm/modnormfact;
    }
    else {
      printf(" *** Warning: Unknown model. Defaulting to CI off model ***\n");
      age=1e6*exactmyears*365*86400;
      t_off=1e6*exactoff*365*86400;

      modnormfact=emiss_nkg(1,b,minf);
      modnormfact=lowagenorm/modnormfact;
    }
    for(i=0; i<=dataintervals; i++) {
      // x-axis frequency for each step
      nu=minf*pow(10,delt*i);
      plotarray[0][i] = emiss_nkg(1,b,nu)*modnormfact;
    }
  }

#pragma omp barrier

  FILE *filestream;

  if ( (filestream = fopen(filename,"w") ) == NULL ) {
    printf("\n*** Error: Cannot open specified directory to output the data file (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in dataloc ***\n\n", filename);
    return 404;
  }

  if (exactage != 1 ) {

    if ( (model == 5) && (varyoffage == 1) ) {

      for(j=0; j<=((maxoff-minoff)/skip); j++) {
	for(i=0; i<=dataintervals; i++) {
	  nu=minf*pow(10,delt*i);
	  fprintf(filestream, "%d %d %.4e %.4e\n", modelmyears, (j*skip)+minoff, nu, plotarray[j][i]);
	}
      }
    }
    else {
      for(j=0; j<=((modelmyears-minmodelmyears)/skip); j++) {
	for(i=0; i<=dataintervals; i++) {
	  nu=minf*pow(10,delt*i);
	  if (model == 5) {
	    fprintf(filestream, "%d %d %.4e %.4e\n", (j*skip)+minmodelmyears, maxoff, nu, plotarray[j][i]);
	  }
	  else if (model == 6) {
	    fprintf(filestream, "%d %.4e %.4e\n", (j*skip)+minmodelmyears, nu, plotarray[j][i]);
	    //plotarray = (float **)calloc(((modelmyears-minmodelmyears)/skip)+1, sizeof(float *));
	  }
	  else {
	    printf(" *** Warning: Unknown model. Defaulting to CI off model ***\n");
	    fprintf(filestream, "%d %d %.4e %.4e\n", (j*skip)+minmodelmyears, (j*skip)+minoff, nu, plotarray[j][i]);
	  }
	}
      }
    }
  }
 else {
   for(i=0; i<=dataintervals; i++) {
     nu=minf*pow(10,delt*i);
     if (model == 5) {
       fprintf(filestream, "%.2f %.2f %.4e %.4e\n", exactmyears, exactoff, nu, plotarray[0][i]);
     }
     else if (model == 6) {
       fprintf(filestream, "%.2f %.4e %.4e\n", exactmyears, nu, plotarray[0][i]);
     }
     else {
	printf(" *** Warning: Unknown model. Defaulting to CI off model ***\n");
	fprintf(filestream, "%.2f %.2f %.4e %.4e\n", exactmyears, exactoff, nu, plotarray[0][i]);
      }
    }
  }

  fclose(filestream);

    
  return 0;
}

