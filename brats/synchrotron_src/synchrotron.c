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
#define INFIN  1.0e2

#ifndef INTCO
#define INTCO (3*ECF)/(4*PI*M_EL*M_EL*M_EL*V_C*V_C*V_C*V_C)
#define EMCO (sqrt(3)*ECF*ECF*ECF)/(4.0*PI*EPS_0*V_C*M_EL)
#endif

static double axf1, axf2;
static double *fx;
static int xls=100;

double fr[MAXDATA],em[MAXDATA],eme[MAXDATA];
double (*ne)(double);
double (*ng)(double);
double (*age_emax)(double);

double gmin,gmax,power,gbreak,gbvalue;
double EMIN,EMAX,EBREAK;
double ageb,age;
double ebrd,gbrd;

/* the above variables are OK as globals because they are not changed by threads. On the other hand, the following ones are! */

double n0_ext;
double BFIELD, nu;
double sinalpha;
double mb_a;

#pragma omp threadprivate(n0_ext,BFIELD,nu,sinalpha,mb_a, age)

gsl_integration_workspace *w1,*w2,*w3;

#pragma omp threadprivate(w1,w2,w3)


double ng_pow(double g) {

  if (g<gmin || g>gmax) return 0.0;
  else return pow(g,-power);
}

double ne_pow(double e) {

  /* ne(E) dE is the number of electrons with energy between E and E+dE */

  /* delta is the power-law index */
  /* normalise to 1 at some arbitrary energy */

  if ((e<EMIN) || (e>EMAX)) return 0.0;
  else return n0_ext*pow(e,-power);
}

double ne_age_jp(double e) {

  /* as above, but with a finite age single-burst model */
  /* loss is the term that corrects for loss. */
  /* JP version */

  double loss;

  if ((e<EMIN) || (e>EMAX)) return(0.0);
  else {
    loss=(4*THOMSON/(6*M_EL*M_EL*V_C*V_C*V_C*MU_0))*ageb*ageb*e*age;
    if (loss>=1.0) return(0.0);
    else return n0_ext*pow(e,-power)*pow((1.0-loss),power-2);
  }
}

double ne_age_kp(double e) {

  /* as above, but with a finite age single-burst model */
  /* loss is the term that corrects for loss. */
  /* KP version */

  double loss;

  if ((e<EMIN) || (e>EMAX)) return(0.0);
  else {
    loss=(THOMSON*sinalpha*sinalpha/(M_EL*M_EL*V_C*V_C*V_C*MU_0))*ageb*ageb*e*age;
    if (loss>=1.0) return(0.0);
    else return n0_ext*pow(e,-power)*pow((1.0-loss),power-2);
  }
}

double ng_age(double g) {

  /* as above, but with a finite age single-burst model */
  /* loss is the term that corrects for loss. */

  double loss;

  if (g<gmin || g>gmax) return 0.0;
  else {
    loss=(4*THOMSON/(6*M_EL*V_C*MU_0))*ageb*ageb*g*age;
    if (loss>=1.0) return(0.0);
    else return pow(g,-power)*pow((1.0-loss),power-2);
  }
}

double age_emax_jp(double limit) {

  double em;

  em=6.0*M_EL*M_EL*V_C*V_C*V_C*MU_0/(4*THOMSON*ageb*ageb*age);
  if (em>limit) em=limit;
  return(em);
}

double age_emax_kp(double limit) {

  double em;

  em=M_EL*M_EL*V_C*V_C*V_C*MU_0/(THOMSON*ageb*ageb*age*sinalpha*sinalpha);
  if (em>limit) em=limit;
  return(em);
}

double ng_break(double g) {
  if (g<gmin || g>gmax) return 0.0;
  else if (g>=gbreak) return gbrd*pow(g,-(power+gbvalue));
  else return pow(g,-power);
}

double ne_break(double e) {

  /* as above, but with a break */
  
  if ((e<EMIN) || (e>EMAX)) return 0.0;
  else if (e>=EBREAK) return n0_ext*ebrd*pow(e,-(power+gbvalue));
  else return n0_ext*pow(e,-power);
}

double ene(double e) {
  return e*ne(e);
}

double ffx(double x)

{
  /* use the tabulated values of F(x) to return F(x) at any x. For x
     outside the table, use the asymptote of Pacholczyk (1970) for small
     x and 0 for large x. Otherwise return a log-linear two-point
     interpolation. */

  int i;
  double av,xv1,xv2;

  if (x>=XF2)
    return 0.0;
  else if (x<=XF1)
    return 4.0*PI*pow((0.5*x),(1.0/3.0))/(sqrt(3.0)*2.68357);
  /* this is a cheat because I'm having problems with the gamma fn */
  else {
    av=log(x);
    xv1=(av-axf1)*((double)xls-1.0)/(axf2-axf1);
    xv2=floor(xv1);
    i=xv2;
    xv1-=xv2;
    return exp(((1.0-xv1)*log(fx[i]))+(xv1*log(fx[i+1])));
  }
}

double ff(double x) {

  float bessk, dummy1, dummy2, dummy3;
  
  /* given x, returns the Bessel function of order 5/3 */
  bessik(x,5.0/3.0,&dummy1,&bessk,&dummy2,&dummy3);
  /*  printf("x was %g, Bessel result was %g\n",x,bessk); */
  return bessk;
}

void makefx(void) {

  int i;
  double xv, yv, av;

  printf("Making F(x) lookup table:\n"); 
  axf1=log(XF1);
  axf2=log(XF2);
#pragma omp parallel for private(i,av,xv,yv)
  for (i=0; i<xls; i++) {
    /* compute the x-value */
    av=axf1+((axf2-axf1)*((double)i)/(double)(xls-1));
    xv=exp(av);
    /* integrate ff from xv to infinity */
    yv=gsl_integ_iu(ff,xv,w1);
    fx[i]=xv*yv;
  }
}

double integr(double e) {

  /* the synchrotron emission integrand */
  return ne(e)*ffx(nu/(INTCO*BFIELD*sinalpha*e*e));
}

double synchin(double alpha) {

  double emin,emax,res;

  sinalpha=sin(alpha);

  emin=sqrt(nu/(BFIELD*sinalpha*INTCO*XF2));
  emax=3*sqrt(nu/(BFIELD*sinalpha*INTCO*XF1));
  if (emin<EMIN) emin=EMIN;
  if (emax>EMAX) emax=EMAX;
  if (age>0.0) emax=age_emax(emax);
  //  else emax=EMAX;
  /* printf("Using min energy %lg\n",emin); */
  /*  printf("in synchpl, emin = %g, emax = %g\n",emin,EMAX); */
  if (emin>emax) {
    /* printf("Synch. integration would be void; no suitable electrons.\n"); */
    res=0;
  } else
    res=gsl_integ(integr,emin,emax,w1);
  if (res<0) {
    printf("synchrotron integral convergence problem\n");
    res=0.0;
  }
  //  printf("returning %g\n",res);
  return EMCO*BFIELD*sinalpha*sinalpha*res;
}

double emiss_n(double n0, double b, double nu_arg) {

  n0_ext=n0;
  BFIELD=b;
  nu=nu_arg;
  return 0.5*gsl_integ(synchin,0.0,PI,w2);;
}

double synchmb(double b) {
  BFIELD=b;

  return b*b*exp(-b*b/(2*mb_a*mb_a))*0.5*gsl_integ(synchin,0.0,PI,w2);
}

double emiss_nmb(double n0, double b, double nu_arg) {

  /* emiss_n integrated over the Maxwell-Boltzmann distribution */

  /* the factor 6 used here is experimentally enough to get the flux to 4 sf */

  n0_ext=n0;
  nu=nu_arg;
  mb_a=b/sqrt(3.0);
  return sqrt(2.0/PI)*gsl_integ(synchmb,0,mb_a*6.0,w3)/pow(mb_a,3.0);
}


int spectralageingmodels(float **regflux, int imgnum, int *fluxorder, float *frequency, int numregions, float *chisquared, float *bestage, float **fluxerror, float *modelalpha, float **modelflux, double usr_gmin, double usr_gmax, int minmyears, int myears, int ageresolution, int levels, float *bestnorm, int modelres, int printresults, double inject, double fieldstrength, int model, int printreject, float redshift, float *ageerrorsplus, float *ageerrorsminus, int suppresscdf) {


  //#define STEPS 50

  // Declare the function prototypes

  double normgoldensearch_func(double normflux, int imgnum, float **fluxerror, float **regflux, int currentregion, int *fluxorder, float *frequency, double magfield, double *emibymapfreq, int ul, int *upperlimits); 

  int goldensearch( double (*f)(double, int, float **, float **, int, int *, float *, double, double *, int, int *), double* a, double *fa, double* b, double* fb, double tolerance, int imgnum, float **fluxerror, float **regflux, int currentregion, int *fluxorder, float *frequency, double magfield, double *emibymapfreq, int ul, int *upperlimits);


  double b, emi, interval, siglevel68, siglevel90, siglevel95, siglevel99, siglevel995, siglevel999, siglevel9999, dof, c, d, fc, fd, **emibymapfreq;
  float nu, normfactor, worstfit, normmax, normmin, tmp_maxage, tmp_minage, tmp_chisquared, confpercentage, maxfitage, minfitage, bcmb;
  int i, a, j, m, l, s, u, currentregion, *medianbin, median, rejectbin, nonrejectbin, oldestregion, youngestregion;

  // Setup some variables
  int ul = 0; // This is the number of datapoints which are strict upper limits (not currently implimented in single injection fits).
  int resolution = ageresolution;
  int regcounter = 0;
  int hundreds = 0;
  int parallelabort = 0;
  float steps = 0.0;
  float avgchisquared = 0.0;
  float avgnormalisation = 0.0;
  float sumchisquared = 0.0;
  float freqsteps = 0.0;
  //float rangecorrection = 0.0;

  float maxage = (float)myears;
  float minage = (float)minmyears;

  double tolerance = -1;

  double *tmp_chibyage, *tmp_chiage, tmp_deltachiof1, tmp_sort, slope, intercept;
  int cub_cnt, beststepindex;

  // This allows for upper limits to be applied to single injection models later on
  int ulindex = 0;
  int *upperlimits = (int *)calloc(ul+1, sizeof(int));

  for (i=0;i<ul;i++) {
    upperlimits[0]=ulindex;
  }


  medianbin = (int *)calloc(5, sizeof(int));


  // Get the confidence level ready for later

  dof = (imgnum - 2);

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
  for(j=0; j<imgnum; j++) {

    printf("Rest frequency: %.3e\n", frequency[fluxorder[j]]);

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

  fx=calloc(xls,sizeof(double));
  makefx();

  printf("Ready to go!\n");

  // Now user set, default 10 and 100000
  
  gmin=usr_gmin;
  gmax=usr_gmax;
  
  // Default of 2.2 gives a spetcral index of 0.6

  power=(2*inject)+1; // Now a user input

  // Set which model to use

  if (model == 1) { // JP Model
    ne=ne_age_jp;
    age_emax=age_emax_jp;
  }

  else if (model == 2) { // KP Model
    ne=ne_age_kp;
    age_emax=age_emax_kp;
  }

  else if (model == 3) { // Tribble model with JP electron energy distribution
    ne=ne_age_jp;
    age_emax=age_emax_jp;
  }

  else if (model == 4) { // Tribble model with KP electron energy distribution
    ne=ne_age_kp;
    age_emax=age_emax_kp;
  }


  else {
    fprintf(stderr,"\nError: Model number is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");

    return 6;
  }

  EMIN=gmin*M_EL*V_C*V_C;
  EMAX=gmax*M_EL*V_C*V_C;

  b=fieldstrength;
  bcmb=(0.318 * (pow(1+redshift,2))) * 1e-9;
  ageb=sqrt(pow(b,2) + pow(bcmb,2));
  //ageb=b;

  printf("Magnetic field strength: %.3e T\n", b);
  printf("B_CMB: %.3e T\n", bcmb);
  printf("Ageing field strength: %.3e T\n", ageb);

  // Get the model values for only those points which have real data available and carry out a chi-square test

  printf("Calculating model fits (this may take some time)\n Currently processing region: 1 of %d ...\n", numregions);
  
  emibymapfreq = (double **)calloc(resolution+1, sizeof(double *));

  for (i=0; i<=resolution; i++){
    emibymapfreq[i] = (double *)calloc(imgnum, sizeof(double));
  }

  b=fieldstrength;
  bcmb=(0.318 * (pow(1+redshift,2))) * 1e-9;
  ageb=sqrt(pow(b,2) + pow(bcmb,2));
  //ageb=b;
 
  for(a=1; a<=numregions; a++) {

    tmp_chibyage = (double *)calloc((resolution+1)*levels, sizeof(double));
    tmp_chiage = (double *)calloc((resolution+1)*levels, sizeof(double));

  
    if (parallelabort != 1) { // If we have bad values, stop trying

      // Set the default chi squared to an unlikely values
      chisquared[a] = 1e22;
      worstfit = -1e22;
      minage = (float)minmyears;
      maxage = (float)myears;

      // For each level of increasing accuracy (changes done at closing bracket)
      for(l=0; l<levels; l++) {

	if (parallelabort != 1) { // If we have bad values, stop trying

#pragma omp parallel for private(j, m, steps, normfactor, normmax, normmin, tmp_minage, tmp_maxage,c,d,fc,fd,) firstprivate(minage, maxage) schedule(dynamic)
	  for(j=0; j<=resolution; j++ ) {

	    if (parallelabort != 1) { // If we have bad values, stop trying

	      if (j == 0) {
		steps = minage;
		//printf("Steps: %.4f\n", steps);
	      }
	      else {
		steps = minage + ( ((maxage-minage)/resolution) * (float)j );
		//printf("Steps: %.4f\n", steps);
	      }

	      // Age in seconds
	      age=1e6*steps*365*86400;

	      //printf("j = %d Age = %.2e steps %.1f minage: %.2e maxage: %.2e resolution: %d myears: %d\n", j, age, steps, minage, maxage, resolution, myears);

	      // Go through each map and find the range of possible normalisation values

	      normmax = -1e36;
	      normmin = 1e36;

#pragma omp parallel for private(m) schedule(dynamic)
	      for(m=0; m<imgnum; m++) {

		// Stop trying if we have bad values
		if (parallelabort != 1) {

		  if ( (model == 3) || (model == 4) ) {
		    emibymapfreq[j][m]=emiss_nmb(1,b,frequency[fluxorder[m]]);
		  }
		  else {
		    emibymapfreq[j][m]=emiss_n(1,b,frequency[fluxorder[m]]);
		  }

		  //printf("emibymapfreq[%d] = %.2e\n", m, emibymapfreq[j][m]);

		  // Check is the flux value is bad and if so make sure the message is only displayed once
		  if ( ( (emibymapfreq[j][m] <= 1e-33) || (emibymapfreq[j][m] >= 1e33) ) && (parallelabort!=1) ) { 

		    printf("\nError: Model has reached its precision limit (model flux effectively 0 or unphysically high) at an age of %.2f Megayears. Please run plot(jp/kp/tribble)model to determine a suitable age range and / or model paramaters (e.g. bfield)\n\n", steps);

		    parallelabort = 1;

		  }
		
		  //printf("level: %d emibymapfreq[%d]: %.4e regflux[%d][%d]: %.4e\n", l, m, emibymapfreq[m], fluxorder[m], a, regflux[fluxorder[m]][a]);
	
		  normfactor = regflux[fluxorder[m]][a]/(emibymapfreq[j][m]);
		
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
                              
		fc = normgoldensearch_func(c, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, b, emibymapfreq[j], ul, upperlimits);
		fd = normgoldensearch_func(d, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, b, emibymapfreq[j], ul, upperlimits);                                          
	      
	    
		goldensearch(normgoldensearch_func, &c, &fc, &d, &fd, tolerance, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, b, emibymapfreq[j], ul, upperlimits);
	    
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
		  beststepindex = j+(l*(resolution+1));
		  //printf("beststepindex: %d\n",beststepindex);

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
	  
	    // Set the best errors for this age (regardless if best overall)
	    if (fc < fd) {
	      tmp_chibyage[j+(l*(resolution+1))] = (float)fc;
	      tmp_chiage[j+(l*(resolution+1))] = steps;
	    }
	    else {
	      tmp_chibyage[j+(l*(resolution+1))] = (float)fd;
	      tmp_chiage[j+(l*(resolution+1))] = steps;
	    }

	  }
	  //printf("Bestage[%d] = %.2f l=%d resolution = %d\n", a, bestage[a], l, resolution);

	  if (parallelabort != 1) {

	    // Set the new min and max ages
	    tmp_maxage = maxage;
	    tmp_minage = minage;

	    // Set the new min and max levels
	    if ( ( bestage[a] - ( (  ( (tmp_maxage-tmp_minage) / resolution ) /2)*pow(10, -l) ) ) <= 0) { // Prevent an impossible age
	      minage = 0.0;
	      //printf("Minage: %.4e\n", minage);
	    }
	    else {
	      minage = bestage[a] - ( (  ( (tmp_maxage-tmp_minage) / resolution ) /2)*pow(10, -l) );
	      //printf("Minage: %.4e\n", minage);
	    }
	    maxage = bestage[a] + ( (  ( (tmp_maxage-tmp_minage) / resolution ) /2)*pow(10, -l) );
	    //printf("Maxage: %.4e\n", maxage);
	  }
	}
      }
      // If we have bad values, kick us out
      if (parallelabort == 1) {
	return 5;
      }
  
      if (printresults == 1) {
	printf("Best fit for %d: %.2f Myr, Reduced X^2 %.2f and Best Norm: %.2e \n", a, bestage[a], ( chisquared[a] / (imgnum -2) ), bestnorm[a]);
      }
  	
      //printf("local minimum is between (%12.6e,%12.6e) and (%12.6e,%12.6e)\n", c,fc, d,fd);
  
      regcounter++;
      avgchisquared += chisquared[a];
      avgnormalisation += bestnorm[a];

      if (suppresscdf !=1) {
	// Bin the chisquared region values by confidence level
	if (chisquared[a] < siglevel68) {
	  medianbin[0]++;
	}
	else if ( (chisquared[a] >= siglevel68) && (chisquared[a] < siglevel90) ) {
	  medianbin[1]++;
	}
	else if ( (chisquared[a] >= siglevel90) && (chisquared[a] < siglevel95) ) {
	  medianbin[2]++;
	}
	else if ( (chisquared[a] >= siglevel95) && (chisquared[a] < siglevel99) ) {
	  medianbin[3]++;
	}
	else if (chisquared[a] >= siglevel99) {
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
    if (parallelabort == 1) {
      return 5;
    }


    // Order the arrays by age for curve fitting later

    for(s=0; s<(resolution+1)*levels; s++) {
      for(u=s; u<(resolution+1)*levels; u++) {
	if(tmp_chiage[s] > tmp_chiage[u]) {

	  tmp_sort=tmp_chiage[s];
	  tmp_chiage[s]=tmp_chiage[u];
	  tmp_chiage[u]=tmp_sort;

	  tmp_sort=tmp_chibyage[s];
	  tmp_chibyage[s]=tmp_chibyage[u];
	  tmp_chibyage[u]=tmp_sort;

	  if (u == beststepindex) {
	    beststepindex=s;
	  }
	}
      }
    }
  
    // Currently only linear interpolation works, but this allows for the option of adding a cspline later.

    int interptype = 0;

    if (interptype == 0) { // Linearlly interpolate between the clostest points to a delta chi^2 of 1
      for (cub_cnt=beststepindex; cub_cnt>=0; cub_cnt--) {

	if ( (tmp_chibyage[cub_cnt]-tmp_chibyage[beststepindex] >= 1.000) && (tmp_chiage[cub_cnt] < tmp_chiage[beststepindex]) ) {

	  // printf("Age: %.4f Delta Chi Lower: %.4e - %.4e = %.4e\n", tmp_chiage[cub_cnt], tmp_chibyage[cub_cnt], tmp_chibyage[beststepindex], tmp_chibyage[cub_cnt]-tmp_chibyage[beststepindex]);
	  tmp_deltachiof1 = tmp_chibyage[beststepindex]+1;

	  slope = (tmp_chibyage[cub_cnt+1]-tmp_chibyage[cub_cnt]) / (tmp_chiage[cub_cnt+1]-tmp_chiage[cub_cnt]);
	  intercept = tmp_chibyage[cub_cnt+1] - (slope * tmp_chiage[cub_cnt+1]);

	  //printf("Minus - Slope: %.4f Intercept: %.4f\n", slope, intercept);

	  ageerrorsminus[a] = fabs( tmp_chiage[beststepindex] - ((tmp_deltachiof1-intercept) / slope) );

	  break;
	}

	// Else if we have hit 0 and still not got a delta X^2 of 1...
	else if ( (cub_cnt == 0) && (tmp_chibyage[cub_cnt]-tmp_chibyage[beststepindex] < 1.000) ) {
	  ageerrorsminus[a] = tmp_chiage[beststepindex];
	}


      }

      for (cub_cnt=beststepindex; cub_cnt<=(resolution+1)*levels; cub_cnt++) {

	if ( (tmp_chibyage[cub_cnt]-tmp_chibyage[beststepindex] >= 1.000)  && (tmp_chiage[cub_cnt] > tmp_chiage[beststepindex]) ) {

	  //printf("Age: %.4f Delta Chi Upper: %.4e - %.4e = %.4e\n", tmp_chiage[cub_cnt], tmp_chibyage[cub_cnt], tmp_chibyage[beststepindex], tmp_chibyage[cub_cnt]-tmp_chibyage[beststepindex]);

	  tmp_deltachiof1 = tmp_chibyage[beststepindex]+1;

	  slope = (tmp_chibyage[cub_cnt-1]-tmp_chibyage[cub_cnt]) / (tmp_chiage[cub_cnt-1]-tmp_chiage[cub_cnt]);
	  intercept = tmp_chibyage[cub_cnt-1] - (slope * tmp_chiage[cub_cnt-1]);

	  //printf("Plus - Slope: %.4f Intercept: %.4f tmp_chiage[beststepindex]: %.4f tmp_deltachiof1: %.4f\n", slope, intercept, tmp_chiage[beststepindex], tmp_deltachiof1);

	  ageerrorsplus[a] = fabs( tmp_chiage[beststepindex] - ((tmp_deltachiof1-intercept) / slope) );

	  break;
	}

	// Else if we have hit the end of our data points and still not got a delta X^2 of 1...
	else if ( (cub_cnt == ((resolution+1)*levels)) && (tmp_chibyage[cub_cnt]-tmp_chibyage[beststepindex] < 1.000) ) {
	  ageerrorsplus[a] = 0.000000;
	  printf("*** Warning: Could not find an upper error for region %d! Error has been set to 0. Consider expanding the age search range (myears). ***\n", a);
	}

      }
    }
      /*
    else { // For fitting a c spline. Currently the data points can be so close you run in to a division by 0. Add this when have more time to think about it.
      
      gsl_interp_accel *acc = gsl_interp_accel_alloc ();
      gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, ((resolution+1)*levels) );

      gsl_spline_init (spline, tmp_chibyage, tmp_chiage,((resolution+1)*levels) );

      tmp_deltachiof1 = tmp_chibyage[beststepindex]+1;

      ageerrorsplus[a] = gsl_spline_eval (spline, tmp_deltachiof1, acc);
      ageerrorsplus[a] = tmp_chiage[beststepindex] - ageerrorsplus[a];

      gsl_spline_free (spline);
      gsl_interp_accel_free (acc);
      
      }*/
  
    // printf("Region: %d Age: %.4f Plus: %.4e Minus: %.4e\n", a, bestage[a], ageerrorsplus[a], ageerrorsminus[a]);

    free(tmp_chibyage);
    free(tmp_chiage);

  }

  sumchisquared = avgchisquared;
  avgchisquared /= numregions;
  avgnormalisation /= numregions;

  // Get the model fluxes and spectral indices
  regcounter = 0;
  hundreds = 0;

  // Check min and max frequencies aren't the same
  if (fluxorder[0] != fluxorder[imgnum-1]) {

    printf("Calculating model fluxes and spectral indices\n Currently processing region: 1 of %d ...\n", numregions);

    steps = 0;

    // Get the predicted flux for each age at myears*resolution intervals
#pragma omp parallel for private(i, j, nu, emi, freqsteps) schedule(dynamic)
    for(i=1; i<=numregions; i++) {
 
      //Best model age fit
      age=1e6*bestage[i]*365*86400;

      // No need for normalisation here as all alpha and curvature are relative
      //#pragma omp parallel for private(j, nu, emi, freqsteps) firstprivate(i) schedule(dynamic)
      for (j=0; j<=modelres; j++) {

	freqsteps = ( (frequency[fluxorder[imgnum-1]] - frequency[fluxorder[0]]) / modelres ) * j;

	// Region Flux
	nu=frequency[fluxorder[0]] + freqsteps;

	// printf("fluxorder[0]= %d   frequency[fluxorder[0]] = %.4e   freqsteps=%.4e \n", fluxorder[0], frequency[fluxorder[0]], freqsteps);

	if ( (model == 3) || (model == 4) ) {
	  emi=emiss_nmb(1,b,nu);
	}
	else {
	  emi=emiss_n(1,b,nu);
	}

	modelflux[i][j] = emi;
	modelflux[i][j] *= bestnorm[i];
      }
      
      // Find alpha for the best model fit of each region
    
      modelalpha[i] = -log(modelflux[i][0]/modelflux[i][modelres]) / log(frequency[fluxorder[0]]/frequency[fluxorder[imgnum-1]]);

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
    return 0;
    //exit(0);
  }


  if (suppresscdf !=1) {
    // Print out the confidence levels
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
  minfitage = 1e22;
  oldestregion = -1;
  youngestregion = -1;

  for(a=1; a<=numregions; a++) {
    if (bestage[a] > maxfitage) {
      maxfitage = bestage[a];
      oldestregion = a;
    }

    if (bestage[a] < minfitage) {
      minfitage = bestage[a];
      youngestregion = a;
    }
  }

  // Print out the average chi-squared value
  printf("Sum of X^2 %.2f \n", sumchisquared);
  printf("Average X^2 %.2f \n", avgchisquared);
  printf("Average normalisation %.2e \n", avgnormalisation);
  printf("Average Reduced X^2 %.2f for %d degrees of freedom\n", ( avgchisquared / (imgnum -2)), (imgnum -2) );
  printf("Minimum Age: %.2f +%.2f -%.2f Megayears (Region %d)\n", minfitage, ageerrorsplus[youngestregion], ageerrorsminus[youngestregion], youngestregion);
  printf("Maximum Age: %.2f +%.2f -%.2f Megayears (Region %d)\n\n", maxfitage, ageerrorsplus[oldestregion], ageerrorsminus[oldestregion], oldestregion);


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


    printf("Bin Statistics:\n");
    printf("===========================================================\n\n");
    printf(" Index  |   Confidence Level   |   Regions   |    Fraction   \n\n");
    printf("   0    |    < 68 per cent     |     %d      |      %.2f     \n", medianbin[0], (float)medianbin[0]/(float)numregions);
    printf("   1    |   68 - 90 per cent   |     %d      |      %.2f     \n", medianbin[1], (float)medianbin[1]/(float)numregions);
    printf("   2    |   90 - 95 per cent   |     %d      |      %.2f     \n", medianbin[2], (float)medianbin[2]/(float)numregions);
    printf("   3    |   95 - 99 per cent   |     %d      |      %.2f     \n", medianbin[3], (float)medianbin[3]/(float)numregions);
    printf("   4    |    > 99 per cent     |     %d      |      %.2f     \n", medianbin[4], (float)medianbin[4]/(float)numregions);
    printf("\n===========================================================\n\n");

    // Find the median of the binned chi-squared values and output the result 
    if (nonrejectbin > rejectbin) {

      median = 0;
      for (i=1; i<=2; i++) {
	if (medianbin[i] > medianbin[median]) {
	  median = i;
	}
	if ( (medianbin[i] == medianbin[median]) && (i != median ) ) {
	  printf("\n *** Warning: Median bin %d and %d contain the same number of regions! *** \n\n", i, median);
	}
      }

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

    else {

      median = 3;
      for (i=4; i<=4; i++) {
	if (medianbin[i] > medianbin[median]) {
	  median = i;
	}
	if ( (medianbin[i] == medianbin[median]) && (i != median ) ) {
	  printf("\n *** Warning: Median bin %d and %d contain the same number of regions! *** \n\n", i, median);
	}
      }

      if (median == 3) {
	confpercentage = 95.0;
      }
      else if (median == 4) {
	confpercentage = 99.0;
      }

      printf("Model is rejected at the %.0f per cent significance level based on median binning\n\n", confpercentage);

    }
  }

  // Find if the model can be rejected or not using the average chi-squared

  if ( (printreject == 1) && (suppresscdf !=1) ) {

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

  for (i=0; i<=resolution; i++){
    free(emibymapfreq[i]);
  }
  free(emibymapfreq);
  free(medianbin);

  return 0;
}




int plotmodel(float minf, float maxf, double inject, double fieldstrength, int model, double usr_gmin, double usr_gmax, int minmodelmyears, int modelmyears, char *output, int titles, float redshift, int skip) {

  // Set up and plot a set of aged spectra

#define STEPS 100

  // If the min in greater than the max, return an error
  if (minmodelmyears > modelmyears) {
    printf("\nError: The minimum on time (%d Myr) must be less than the maximum on time (%d Myr). Please use the minmodelmyears and/or the modelmyears commands to set a new range.\n\n", minmodelmyears, modelmyears);
    return 6;
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
  double b, emi;
  float delt, yaxismax, yaxismin, nu,bcmb;
  float plotarray[((modelmyears-minmodelmyears)/skip)+1][STEPS+1];
  char modeltitle[128];

#pragma omp parallel 
  {
    w1=gsl_integration_workspace_alloc(GSL_WSIZE);
    w2=gsl_integration_workspace_alloc(GSL_WSIZE);
    w3=gsl_integration_workspace_alloc(GSL_WSIZE);
  }

  gsl_set_error_handler_off(); // to avoid abort
  gsl_epsilon=1e-3;
  gsl_quiet=1;

  fx=calloc(xls,sizeof(double));
  makefx();

  printf("Ready to go!\n");

  gmin=usr_gmin;
  gmax=usr_gmax;

  // Default of 2.2 gives a spetcral index of 0.6
  power=(2*inject)+1; //Now a user input

  // Set which model to use
  if (model == 1) { // JP Model
    ne=ne_age_jp;
    age_emax=age_emax_jp;
  }

  else if (model == 2) { // KP Model
    ne=ne_age_kp;
    age_emax=age_emax_kp;
  }

  else if (model == 3) { // Tribble model with JP electron energy distribution
    ne=ne_age_jp;
    age_emax=age_emax_jp;
  }

  else if (model == 4) { // Tribble model with KP electron energy distribution
    ne=ne_age_kp;
    age_emax=age_emax_kp;
  }

  EMIN=gmin*M_EL*V_C*V_C;
  EMAX=gmax*M_EL*V_C*V_C;

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

  // Calculate a sensible Y axis scale
  nu=minf;
  age=0.0;

  if ( (model == 3) || (model == 4) ) {
    emi=emiss_nmb(1,b,nu);
  }
  else {
    emi=emiss_n(1,b,nu);
  }

  if(log10(emi) >= 1) {
    yaxismax = log10(emi) * 1.05;
  }
  else {
    yaxismax = log10(emi) * 0.95;
  }

  nu=maxf;
  age=1e6*modelmyears*365*86400;

  if ( (model == 3) || (model == 4) ) {
    emi=emiss_nmb(1,b,nu);
  }
  else {
    emi=emiss_n(1,b,nu);
  }

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

  // yaxismin = log10(emi) - 12;

  cpgenv(log10(minf),log10(maxf), yaxismin, yaxismax, 0, 30);

  // Set the title if required

  if (titles == 0) { sprintf(modeltitle," ");  }
  else {

    // Set which model to use
    if (model == 1) { // JP Model
      sprintf(modeltitle,"JP Model Between %d and %d MYears (Arbitrary Normalisation)", minmodelmyears, modelmyears);
    }
    else if (model == 2) { // KP Model
      sprintf(modeltitle,"KP Model Between %d and %d MYears (Arbitrary Normalisation)", minmodelmyears, modelmyears);
    }
    else if (model == 3) { // Tribble model with JP electron energy distribution
      sprintf(modeltitle,"Tribble (JP) Model Between %d and %d MYears (Arbitrary Normalisation)", minmodelmyears, modelmyears);
    }
    else if (model == 4) { // Tribble model with KP electron energy distribution
      sprintf(modeltitle,"Tribble (KP) Model Between %d and %d MYears (Arbitrary Normalisation)", minmodelmyears, modelmyears);
    }
    else {
      sprintf(modeltitle,"Unknown Model Between %d and %d MYears (Arbitrary Normalisation)", minmodelmyears, modelmyears);
    }
  }


  cpglab("Frequency / Hz", "Flux / Jy", modeltitle);

  delt=log10(maxf/minf)/(float)STEPS;

  //firstprivate(plotarray)

#pragma omp parallel for private(j, i, emi, nu) schedule(dynamic)
  for(j=0; j<=((modelmyears-minmodelmyears)/skip); j++) {
    
    age=1e6*(minmodelmyears+(j*skip))*365*86400;

    //printf("j = %d Age = %.2e\n",j , age);

    //#pragma omp parallel for private(i, nu, emi) schedule(dynamic)
    for(i=0; i<=STEPS; i++) {
      // x-axis frequency for each step
      nu=minf*pow(10,delt*i);

      // y-axis emission for each step
      if ( (model == 1) || (model == 2) ) {
	plotarray[j][i] = emiss_n(1,b,nu);
      }
      else if ( (model == 3) || (model == 4) ) {
	plotarray[j][i] = emiss_nmb(1,b,nu);
      }
    }
  }

  // Easy way to avoid the race condition on colours
#pragma omp barrier

  for(j=0; j<=((modelmyears-minmodelmyears)/skip); j++) {
    cpgsci(j+2); //The +2 avoid the background and gives the standard 0 to 10 start and end distinctive colours.

    for(i=0; i<=STEPS; i++) {
      nu=minf*pow(10,delt*i);
      (i?cpgdraw:cpgmove)(log10(nu),log10(plotarray[j][i]));
    }
  }
 
  cpgclos();

  return 0;
}





int modeldata(float minf, float maxf, double inject, double fieldstrength, int model, double usr_gmin, double usr_gmax, int minmodelmyears, int modelmyears, char *filename, float redshift, int exactage, float exactmyears, int dataintervals) {


  // If the min in greater than the max, return an error
  if ( (minmodelmyears > modelmyears) && (exactage != 1) ) {
    printf("\nError: The minimum on time (%d Myr) must be less than the maximum on time (%d Myr). Please use the minmodelmyears and/or the modelmyears commands to set a new range.\n\n", minmodelmyears, modelmyears);
    return 6;
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
  double b;
  float delt, nu, bcmb, **plotarray;

  if (exactage != 1 ) {
    plotarray = (float **)calloc((modelmyears-minmodelmyears)+1, sizeof(float *)); 
    for (i=0; i<=(modelmyears-minmodelmyears); i++) {
      plotarray[i] = (float *)calloc(dataintervals+1, sizeof(float)); 
    }
  }
  else {
    plotarray = (float **)calloc(1, sizeof(float *)); 
    plotarray[0] = (float *)calloc(dataintervals+1, sizeof(float)); 
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

  fx=calloc(xls,sizeof(double));
  makefx();

  printf("Ready to go!\n");

  gmin=usr_gmin;
  gmax=usr_gmax;

  // Default of 2.2 gives a spetcral index of 0.6
  power=(2*inject)+1; //Now a user input

  // Set which model to use
  if (model == 1) { // JP Model
    ne=ne_age_jp;
    age_emax=age_emax_jp;
  }

  else if (model == 2) { // KP Model
    ne=ne_age_kp;
    age_emax=age_emax_kp;
  }

  else if (model == 3) { // Tribble model with JP electron energy distribution
    ne=ne_age_jp;
    age_emax=age_emax_jp;
  }

  else if (model == 4) { // Tribble model with KP electron energy distribution
    ne=ne_age_kp;
    age_emax=age_emax_kp;
  }

  EMIN=gmin*M_EL*V_C*V_C;
  EMAX=gmax*M_EL*V_C*V_C;

  b=fieldstrength;
  bcmb=(0.318 * (pow(1+redshift,2))) * 1e-9;
  ageb=sqrt(pow(b,2) + pow(bcmb,2));
  //ageb=b;

  printf("Magnetic field strength: %.3e T\n", b);
  printf("B_CMB: %.3e T\n", bcmb);
  printf("Ageing field strength: %.3e T\n", ageb);


  delt=log10(maxf/minf)/(float)dataintervals;

  if (exactage != 1) {

#pragma omp parallel for private(j, i, nu) schedule(dynamic)
    for(j=0; j<=(modelmyears-minmodelmyears); j++) {
    
      age=1e6*(minmodelmyears+j)*365*86400;

      for(i=0; i<=dataintervals; i++) {
	// x-axis frequency for each step
	nu=minf*pow(10,delt*i);

	// y-axis emission for each step
	if ( (model == 1) || (model == 2) ) {
	  plotarray[j][i] = emiss_n(1,b,nu);
	}
	else if ( (model == 3) || (model == 4) ) {
	  plotarray[j][i] = emiss_nmb(1,b,nu);
	}
      }
    }
  }
  else {
  
    age=1e6*exactmyears*365*86400;

#pragma omp parallel for private(i, nu) schedule(dynamic)
    for(i=0; i<=dataintervals; i++) {
      // x-axis frequency for each step
      nu=minf*pow(10,delt*i);
      age=1e6*exactmyears*365*86400;

      // y-axis emission for each step
      if ( (model == 1) || (model == 2) ) {
	plotarray[0][i] = emiss_n(1,b,nu);
      }
      else if ( (model == 3) || (model == 4) ) {
	plotarray[0][i] = emiss_nmb(1,b,nu);
      }
    }
  }
  // Easy way to avoid the race condition on colours
#pragma omp barrier

  FILE *filestream;

  if ( (filestream = fopen(filename,"w") ) == NULL ) {
    printf("\n*** Error: Cannot open specified directory to output the data file (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in dataloc ***\n\n", filename);
    return 404;
  }

  if (exactage != 1 ) {
    for(j=0; j<=(modelmyears-minmodelmyears); j++) {
      for(i=0; i<=dataintervals; i++) {
	nu=minf*pow(10,delt*i);
	fprintf(filestream, "%d %.4e %.4e\n", minmodelmyears+j, nu, plotarray[j][i]);
      }
    }
  }
  else {
    for(i=0; i<=dataintervals; i++) {
      nu=minf*pow(10,delt*i);
      fprintf(filestream, "%.2f %.4e %.4e\n", exactmyears, nu, plotarray[0][i]);
    }
  }

  fclose(filestream);

  return 0;
}





// f(x) for golden search of minumum chi-squared

double normgoldensearch_func(double normflux, int imgnum, float **fluxerror, float **regflux, int currentregion, int *fluxorder, float *frequency, double magfield, double *emibymapfreq, int ul, int *upperlimits) {

  int i;
  double chisquared = 0.0;
  float testflux = 0.0;

  if (ul > 0) { // If we have strict upper limits, make sure we dont go over them
    for(i=0; i<ul; i++) {
      if ((pow(10, normflux) * emibymapfreq[upperlimits[i]]) >= regflux[upperlimits[i]][currentregion])  { // If the normalised flux is over the upper limit, send back a high chisquared value
	return 1e33;
      }
    }
  }

  for(i=0; i<imgnum-ul; i++) {

    // Get the normalised model flux value
    testflux = pow(10, normflux) * emibymapfreq[i];

    //Calculate the chi squared value
    chisquared += (pow(fabs((regflux[fluxorder[i]][currentregion])-testflux)/(fluxerror[fluxorder[i]][currentregion]), 2) );

  }
 
  
  return chisquared;
}

//Internally Defined Routines 
static int Stopping_Rule(double x0, double x1, double tolerance);

#define sqrt5 2.236067977499789696

int goldensearch( double (*f)(double, int, float **, float **, int, int *, float *, double, double *, int, int*), double *a, double *fa, double* b, double* fb, double tolerance, int imgnum, float **fluxerror, float **regflux, int currentregion, int *fluxorder, float *frequency, double magfield, double *emibymapfreq, int ul, int *upperlimits)
{
  static const double lambda = 0.5 * (sqrt5 - 1.0);
  static const double mu = 0.5 * (3.0 - sqrt5);         // = 1 - lambda
  double x1;
  double x2;
  double fx1;
  double fx2;

  // Find first two internal points and evaluate 
  // the function at the two internal points.

  x1 = *b - lambda * (*b - *a);
  x2 = *a + lambda * (*b - *a);

  fx1 = normgoldensearch_func(x1, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, magfield, emibymapfreq, ul, upperlimits);
  fx2 = normgoldensearch_func(x2, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, magfield, emibymapfreq, ul, upperlimits);
 
  // Verify that the tolerance is an acceptable number

  if (tolerance <= 0.0) tolerance = sqrt(DBL_EPSILON) * (*b - *a);

  // Loop by excluding segments from current endpoints a, b
  // to current internal points x1, x2 and then calculating
  // a new internal point until the length of the interval
  // is less than or equal to the tolerance.

  while ( ! Stopping_Rule( *a, *b, tolerance) ) {
    if (fx1 > fx2) {
      *a = x1;
      *fa = fx1;
      if ( Stopping_Rule( *a, *b, tolerance) ) break;
      x1 = x2;
      fx1 = fx2;
      x2 = *b - mu * (*b - *a);
      fx2 = normgoldensearch_func(x2, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, magfield, emibymapfreq, ul, upperlimits);
    } else {
      *b = x2;
      *fb = fx2;
      if ( Stopping_Rule( *a, *b, tolerance) ) break;
      x2 = x1;
      fx2 = fx1;
      x1 = *a + mu * (*b - *a);
      fx1 = normgoldensearch_func(x1, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, magfield, emibymapfreq, ul, upperlimits);
    }
  }

  return 0;
}



static int Stopping_Rule(double x0, double x1, double tolerance) 
{
  double xm = 0.5 * fabs( x1 + x0 );

  if ( xm <= 1.0 ) return ( fabs( x1 - x0 ) < tolerance ) ? 1 : 0;
  return ( fabs( x1 - x0 ) < tolerance * xm ) ? 1 : 0;
}


int spectralageing_min(float **regflux, int imgnum, int *fluxorder, float *frequency, int numregions, float **fluxerror, double usr_gmin, double usr_gmax, int minmyears, int myears, int ageresolution, int levels, double inject_loop, double fieldstrength, int model, float *sumchisquared, float *regchisquared, float *regbestinject, float redshift) {

  // Declare the function prototypes

  double normgoldensearch_func(double normflux, int imgnum, float **fluxerror, float **regflux, int currentregion, int *fluxorder, float *frequency, double magfield, double *emibymapfreq, int ul, int *upperlimits); 

  int goldensearch( double (*f)(double, int, float **, float **, int, int *, float *, double, double *, int, int *), double* a, double *fa, double* b, double* fb, double tolerance, int imgnum, float **fluxerror, float **regflux, int currentregion, int *fluxorder, float *frequency, double magfield, double *emibymapfreq, int ul, int *upperlimits);


  double b, c, d, fc, fd, **emibymapfreq;
  float normfactor, normmax, normmin, bestage, tmp_maxage, tmp_minage, tmp_chisquared;
  int a, i, j, m, l, currentregion;

  // Setup some variables
  int ul = 0; // This is for strict upper limits. Not currently implimented in singhle injection fits.
  int resolution = ageresolution;
  int regcounter = 0;
  int hundreds = 0;
  int parallelabort = 0;
  float steps = 0.0;
  //float rangecorrection = 0.0;
  float maxage = (float)myears;
  float minage = (float)minmyears;
  float bcmb;

  double tolerance = -1;

  float *mininj_chisquared;

  mininj_chisquared = (float *)calloc(numregions+1, sizeof(float)); 


  // This allows for upper limits to be applied to single injection models later on
  int ulindex = 0;
  int *upperlimits = (int *)calloc(ul+1, sizeof(int));

  for (i=0;i<ul;i++) {
    upperlimits[0]=ulindex;
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

  fx=calloc(xls,sizeof(double));
  makefx();

  // Now user set, default 10 and 100000
  gmin=usr_gmin;
  gmax=usr_gmax;

  // Default of 2.2 gives a spetcral index of 0.6
  power=(2*inject_loop)+1; // Now a user input

  // Set which model to use

  if (model == 1) { // JP Model
    ne=ne_age_jp;
    age_emax=age_emax_jp;
  }

  else if (model == 2) { // KP Model
    ne=ne_age_kp;
    age_emax=age_emax_kp;
  }

  else if (model == 3) { // Tribble model with JP electron energy distribution
    ne=ne_age_jp;
    age_emax=age_emax_jp;
  }

  else if (model == 4) { // Tribble model with KP electron energy distribution
    ne=ne_age_kp;
    age_emax=age_emax_kp;
  }


  else {
    fprintf(stderr,"\nError: Model number is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");

    return 6;
  }

  EMIN=gmin*M_EL*V_C*V_C;
  EMAX=gmax*M_EL*V_C*V_C;

  b=fieldstrength;
  bcmb=(0.318 * (pow(1+redshift,2))) * 1e-9;
  ageb=sqrt(pow(b,2) + pow(bcmb,2));
  //ageb=b;

  // Get the model values for only those points which have real data available and carry out a chi-square test

  printf("Currently processing region: 1 of %d ...\n", numregions);

  emibymapfreq = (double **)calloc(resolution+1, sizeof(double *));

  for (i=0; i<=resolution; i++){
    emibymapfreq[i] = (double *)calloc(imgnum, sizeof(double));
  }

  b=fieldstrength;
  bcmb=(0.318 * (pow(1+redshift,2))) * 1e-9;
  ageb=sqrt(pow(b,2) + pow(bcmb,2));
  //ageb=b;

  for(a=1; a<=numregions; a++) {

    if (parallelabort != 1) { // If we have bad values, stop trying

      // Set the default chi squared to an unlikely values
      mininj_chisquared[a] = 1e+22;
      minage = (float)minmyears;
      maxage = (float)myears;

      // For each level of increasing accuracy (changes done at closing bracket)
      for(l=0; l<levels; l++) {

	if (parallelabort != 1) {

#pragma omp parallel for private(j, m, steps, normfactor, normmax, normmin, tmp_minage, tmp_maxage,c,d,fc,fd) firstprivate(minage, maxage) schedule(dynamic)
	  for(j=0; j<=resolution; j++ ) {

	    if (parallelabort != 1) {

	      if (j == 0) {
		steps = minage;
		//printf("Steps: %.4f\n", steps);
	      }
	      else {
		steps = minage + ( ((maxage-minage)/resolution) * (float)j );
		//printf("Steps: %.4f\n", steps);
	      }

	      // Age in seconds
	      age=1e6*steps*365*86400;

	      // Go through each map and find the range of possible normalisation values

	      // (Re)set the normalisations range to unlikely values
	      normmax = -1e36;
	      normmin = 1e36;

#pragma omp parallel for private(m) schedule(dynamic)
	      for(m=0; m<imgnum; m++) {

		if (parallelabort != 1) {
  
		  if ( (model == 3) || (model == 4) ) {
		    emibymapfreq[j][m]=emiss_nmb(1,b,frequency[fluxorder[m]]);
		  }
		  else {
		    emibymapfreq[j][m]=emiss_n(1,b,frequency[fluxorder[m]]);
		  }

		  if ( ( (emibymapfreq[j][m] <= 1e-33) || (emibymapfreq[j][m] >= 1e33) ) && (parallelabort != 1) ) {
		    printf("\nError: Model has reached its precision limit (model flux effectively 0 or unphysically high) at an age of %.2f Megayears. Please run plot(JP/KP)model to determine a suitable age range\n\n", steps); 

		    parallelabort = 1;

		  }

		  //printf("level: %d emibymapfreq[%d]: %.4e regflux[%d][%d]: %.4e\n", l, m, emibymapfreq[m], fluxorder[m], a, regflux[fluxorder[m]][a]);
	
		  normfactor = regflux[fluxorder[m]][a]/(emibymapfreq[j][m]);
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

		fc = normgoldensearch_func(c, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, b, emibymapfreq[j], ul, upperlimits);
		fd = normgoldensearch_func(d, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, b, emibymapfreq[j], ul, upperlimits);                                          
		goldensearch(normgoldensearch_func, &c, &fc, &d, &fd, tolerance, imgnum, fluxerror, regflux, currentregion, fluxorder, frequency, b, emibymapfreq[j], ul, upperlimits);
      
		//If the chi squared is less than the previous best set its new value and age
		if ( (fc < mininj_chisquared[a]) || (fd < mininj_chisquared[a]) ) {

		  if (fc < fd) {
		    tmp_chisquared = (float)fc;
		  }
		  else {
		    tmp_chisquared = (float)fd;
		  }

		  mininj_chisquared[a] = tmp_chisquared;
		  bestage = steps;
		}
	      }
	    }
	  }

	  if (parallelabort != 1) {
	    // Set the next min and max ages
	    tmp_maxage = maxage;
	    tmp_minage = minage;

	    // Set the new min and max levels
	    if (bestage - ( (  ( (tmp_maxage-tmp_minage) / resolution ) /2)*pow(10, -l) ) <= 0.0) { // Prevent an impossible age
	      minage = 0.0;
	      //printf("Minage: %.4e\n", minage);
	    }
	    else {
	      minage = bestage - ( (  ( (tmp_maxage-tmp_minage) / resolution ) /2)*pow(10, -l) );
	      //printf("Minage: %.4e\n", minage);
	    }
	    maxage = bestage + ( (  ( (tmp_maxage-tmp_minage) / resolution ) /2)*pow(10, -l) );
	    //printf("Maxage: %.4e\n", maxage);
	  }
	}
      }

      // If we have bad values, kick us out
      if (parallelabort == 1) {
	return 5;
      }
    
      if (mininj_chisquared[a] < regchisquared[a]) {
	regchisquared[a] = mininj_chisquared[a];
	regbestinject[a] = inject_loop;
      }

      regcounter++;
      *sumchisquared+= mininj_chisquared[a];
  
      if (regcounter >= 100) {

	hundreds++;
	printf("%d00...\n", hundreds);

	regcounter = 0;
      }
    }
  }

  // If we have bad values, kick us out
  if (parallelabort == 1) {
    return 5;
  }


  // Print out the average chi-squared value
  printf("Sum of X^2 %.2f for model %d with injection index %.2f\n", *sumchisquared, model, inject_loop);
  printf("Average X^2 %.2f\n", *sumchisquared/numregions);
  // Free up the memory
  free(mininj_chisquared);

  return 0;

}




int findmodelvalues(double usr_gmin, double usr_gmax, double fieldstrength, int model, float *reconmapflux, int xdim, int ydim, float usr_frequency, float *bestnorm, float *bestage, int **regionarray, float inject, float beamarea, float redshift) {


  double b;
  int i, j;
  float bcmb;

  // Setup some variables
  int regcounter = 0;
  int thousands = 0;
  int parallelabort = 0;

#pragma omp parallel 
  {
    w1=gsl_integration_workspace_alloc(GSL_WSIZE);
    w2=gsl_integration_workspace_alloc(GSL_WSIZE);
    w3=gsl_integration_workspace_alloc(GSL_WSIZE);
  }

  gsl_set_error_handler_off(); // to avoid abort
  gsl_epsilon=1e-3;
  gsl_quiet=1;

  fx=calloc(xls,sizeof(double));
  makefx();

  // Now user set, default 10 and 100000
  gmin=usr_gmin;
  gmax=usr_gmax;

  // Default of 2.2 gives a spetcral index of 0.6
  power=(2*inject)+1; // Now a user input

  // Set which model to use

  if (model == 1) { // JP Model
    ne=ne_age_jp;
    age_emax=age_emax_jp;
  }

  else if (model == 2) { // KP Model
    ne=ne_age_kp;
    age_emax=age_emax_kp;
  }

  else if (model == 3) { // Tribble model with JP electron energy distribution
    ne=ne_age_jp;
    age_emax=age_emax_jp;
  }

  else if (model == 4) { // Tribble model with KP electron energy distribution
    ne=ne_age_kp;
    age_emax=age_emax_kp;
  }


  else {
    fprintf(stderr,"\nError: Model number is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");

    return 6;
  }

  EMIN=gmin*M_EL*V_C*V_C;
  EMAX=gmax*M_EL*V_C*V_C;

  b=fieldstrength;
  bcmb=(0.318 * (pow(1+redshift,2))) * 1e-9;
  ageb=sqrt(pow(b,2) + pow(bcmb,2));
  //ageb=b;

  // Get the model values for only those points which have real data available and carry out a chi-square test

  printf("Currently determining flux for pixel: 1 of %d ...\n", xdim*ydim);
  
  b=fieldstrength;
  bcmb=(0.318 * (pow(1+redshift,2))) * 1e-9;
  ageb=sqrt(pow(b,2) + pow(bcmb,2));
  //ageb=b;


  for (i=0; i<xdim; i++){  

    for(j=0; j<ydim; j++) {

      if (regionarray[i][j] <= 0) {

	reconmapflux[i*xdim+j] = 0.0;
      }
      else {

	// Age in seconds
	age=1e6*bestage[regionarray[i][j]]*365*86400;

	if ( (model == 1) || (model == 2) ) {
	  reconmapflux[i*xdim+j] = emiss_n(1, b, usr_frequency);

	  reconmapflux[i*xdim+j] *= bestnorm[regionarray[i][j]]; // Scale by the normalisation factor
	  reconmapflux[i*xdim+j] *= beamarea; // Rescale back to Jy/beam
	}

	else if ( (model == 3) || (model == 4) ) {

	  reconmapflux[i*xdim+j] = emiss_nmb(1, b, usr_frequency);

	  reconmapflux[i*xdim+j] *= bestnorm[regionarray[i][j]]; // Scale by the normalisation factor
	  reconmapflux[i*xdim+j] *= beamarea; // Rescale back to Jy/beam

	}
      }

      regcounter++;
  
      if (regcounter >= 10000) {

	thousands++;
	printf("%d0000...\n", thousands);

	regcounter = 0;
      }
    }
  }

  // If we have bad values, kick us out
  if (parallelabort == 1) {
    return 5;
  }


  return 0;

}

