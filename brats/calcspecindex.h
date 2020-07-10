/* Find the polynomial regression of a data set.

   What needs to be passed:

   float **regflux - An array containing fluxes as flux[map][region]. A single region is allowable but the second dimension of the array must still be included.
   float *specindex - An empty array to return the spectral indices. Must be of format specindex[region]
   int numregions - The number of regions in the dataset.
   float freq1 - Frequency of maps 1
   float freq2 - Frequency of maps 2
   float *maxsi - Holder for the maximum spectral index
   float *minsi - Holder for the maximum spectral index
   int printindex - Print the spectral index to the terminal?


Basic usage. We assume imgnum is already set to be an integer number of data points (maps), flux is an array of fluxes by region, and map1 and maps 2 are set to the desired map values for flux[map] (normally min and max frequencies):

  float **alpha;
  alpha = calloc(imgnum, sizeof(float *));

  ...

  calcspecindex(flux, alpha, numregions, map1, map2, freq1, freq2, 1);

Plain values (not logs) should be passed. Multiple spectal indices can be obtained using a simple loop.
This script is designed to be combined with the spectra index mapping routine.

   If you have any question please email jeremy.harwood@physics.org and I will try and answer them as soon as I can
*/


#include <gsl/gsl_fit.h>


int calcspecindex(float **regflux,  int **regionarray, float *specindex, int numregions, int imgnum, float *frequency, int printindex, float *minsi, float *maxsi) {

  // Declaring variables

  int i, j;
  double sumX, sumY, sumXX, sumXY;

  *maxsi = -1e32;
  *minsi = 1e32;

 
  // Determine the spectral indices between the maximum and minimum data points (maps)

  for (i = 1; i <=numregions; i++) {

      sumX = 0;
      sumY = 0;
      sumXX = 0;
      sumXY = 0;

    for (j = 0; j < imgnum; j++) {
      sumX += log(frequency[j]);
      sumY += log(regflux[j][i]);
      sumXX += log(frequency[j]) * log(frequency[j]);
      sumXY += log(frequency[j]) * log(regflux[j][i]);
    }



    // printf("sumX: %.2e sumY: %.2e sumXX: %.2e sumXY: %.2e\n",sumX,sumY,sumXX,sumXY);


    specindex[i] = -( ((imgnum * sumXY) - (sumX * sumY)) / ((imgnum * sumXX) - (sumX * sumX)) );

	// Find the max and min spectral indices as we go
	if (specindex[i] > *maxsi) {
	  *maxsi = specindex[i];
	}
	if (specindex[i] < *minsi) {
	  *minsi = specindex[i];
	}


    if (printindex == 1) {

      printf("Spectal Index of region %d: %.2f\n", i, specindex[i]);
    }
  }

  /*
  if (warningcount > 0) {
    fprintf(stderr,"\n*** WARNING: Regions had a flux less than the RMS %d times! This normally occurs when fixed regions are used and they are not a good match ***\n\n", warningcount);
  }
  */

    printf("Maximum Spectal Index: %.2f\n", *maxsi);
    printf("Minimum Spectal Index: %.2f\n", *minsi);

  return 0;

}



// Calculate the spec index using only 2 points e.g. for colour-colour plots 

int calcspecindexp2p(float **regflux, float *specindex, int numregions, int map1, int map2, float freq1, float freq2, int printindex, float *minsi, float *maxsi, float *rms) {

  // Declaring variables

  int i;
  int warningcount = 0;
  *maxsi = -1e32;
  *minsi = 1e32;

 
  // Determine the spectral indices between the maximum and minimum data points (maps)


  for (i = 1; i <=numregions; i++) {

    // Count how many times a regions has a flux below the rms
    if ( (regflux[map1][i] <= rms[map1])  || (regflux[map2][i] <= rms[map2]) ) {
      warningcount++;
      specindex[i] = 1e34; // These values will be ignored later 
      }

    else {
      if (freq1 > freq2) {
	specindex[i] = (-log(regflux[map2][i]/regflux[map1][i])) / (log(freq2/freq1));

	// Find the max and min spectral indices as we go
	if (specindex[i] > *maxsi) {
	  *maxsi = specindex[i];
	}
	if (specindex[i] < *minsi) {
	  *minsi = specindex[i];
	}
      }

      else if (freq2 > freq1) {
	specindex[i] = (-log(regflux[map1][i]/regflux[map2][i])) / (log(freq1/freq2));

	// Find the max and min spectral indices as we go
	if (specindex[i] > *maxsi) {
	  *maxsi = specindex[i];
	}
	if (specindex[i] < *minsi) {
	  *minsi = specindex[i];
	}
      }

      else {
	fprintf(stderr,"\nError: Cannot calculate spectral index. Are both maps the same frequency???\nExiting...\n\n");
	return 404;
      }
    }


    if (printindex == 1) {

      printf("Spectal Index of region %d: %.2f\n", i, specindex[i]);
    }
  }

  /*
  if (warningcount > 0) {
    fprintf(stderr,"\n*** WARNING: Regions had a flux less than the RMS %d times! This normally occurs when fixed regions are used and they are not a good match ***\n\n", warningcount);
  }
  */

    printf("Maximum Spectal Index: %.2f\n", *maxsi);
    printf("Minimum Spectal Index: %.2f\n", *minsi);

  return 0;

}



int calcspecindexweightedgsl(float **regflux, int **regionarray, float *specindex, double *specindexintercept, float *specindexchisquared, double **specindex_modelflux, double **specindex_modelerror, float **fluxerror, int numregions, int imgnum, float *frequency, int printindex, float *minsi, float *maxsi, int modelres, int *fluxorder, float *specindexerror, int *specindex_errortype, int forceerrortype) {

  // Declaring variables
  int i, j, minsiref, maxsiref;
  double c0, c1, cov00, cov01, cov11, chisq, tmp_modelfrequency;
  double **switchedfluxarray, *logfreq, *weights, fracerror, sow, xw, x2w;

  double lntol10 = log10(10)/log(10);

  // Declare the GSL functions
  int gsl_fit_wlinear (const double * x, const size_t xstride, const double * w, const size_t wstride, const double * y, const size_t ystride, size_t n, double * c0, double * c1, double * cov00, double * cov01, double * cov11, double * chisq);

  int gsl_fit_linear_est(double x, double c0, double c1, double cov00, double cov01, double cov11, double * y, double * y_err);

  *maxsi = -1e32;
  *minsi = 1e32;
  minsiref = 0;
  maxsiref = 0;

  // GSL expects the array to be on a map by map basis, so we have to swap the around here (and take the log while we are at it)

  switchedfluxarray = (double **)calloc(numregions+1, sizeof(double *));

  for (i=1; i<=numregions; i++) {
    switchedfluxarray[i] = (double *)calloc(imgnum, sizeof(double));

    for (j=0; j<imgnum; j++) {
      switchedfluxarray[i][j] = log(regflux[j][i]);
    }
  }

  // Take the log of the frequency for each map
  logfreq = (double *)calloc(imgnum, sizeof(double));

  for (j = 0; j<imgnum; j++) {
    logfreq[j] = log(frequency[j]);
  }


  // Determine the spectral indices between the maximum and minimum data points (maps)

  for (i = 1; i <=numregions; i++) {

    // Calculate the weighting of each map for the region
    weights = (double *)calloc(imgnum, sizeof(double));
  
    // Set the weights
    for (j = 0; j<imgnum; j++) {
      fracerror = (fluxerror[j][i] / regflux[j][i]);
      weights[j] = 1/(pow( (switchedfluxarray[i][j]*fracerror),2));
    }

    gsl_fit_wlinear(logfreq, 1, weights, 1, switchedfluxarray[i], 1, imgnum, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
    
    specindex[i] = -c1;
    specindexintercept[i] = c0;
    specindexchisquared[i] = (float)chisq; // This has to be a float to be accepted by the generic chi squared mapping function


    if ( (imgnum == 2) && (forceerrortype == 0) ) { // Unless we are forcing StdDev, if we only have 2 images just propagate the errors...

      float sumfracerror = 0;
      float logfreqratio = 0;

      for (j = 0; j<imgnum; j++) { // Calculate the sum of the squares.
	sumfracerror += pow((fluxerror[j][i] / regflux[j][i]), 2);
      }
	
      logfreqratio = log(frequency[fluxorder[0]]/frequency[fluxorder[1]]); // Find the log if the ratio of frequnecies
      
      specindexerror[i] = fabs((1/logfreqratio) * sqrt(sumfracerror)); // Determine the error
      
      *specindex_errortype = 0;

    } else { // ... but if we have more, use the standard deviation for WLS
      // Calculate the error for each alpha
      sow = 0.0;
      xw = 0.0;
      x2w = 0.0;
    
      for (j = 0; j<imgnum; j++) { // Calculate the components
	sow += weights[j];
        xw += logfreq[j] * weights[j];
	x2w += pow(logfreq[j],2) * weights[j];
      }

      specindexerror[i] = sqrt(sow/((x2w*sow)-pow(xw, 2))); // Determine the error
      *specindex_errortype = 1;
    } //YYYYYYYY


    //printf("SI: %.3f +- %.3f\n", specindex[i], specindexerror[i]);

    //if ( (imgnum != 2) || ( (imgnum == 2) && (forceerrortype != 0) ) ) { // Check we have the right conditions
    
      for (j = 0; j<=modelres; j++) {

	tmp_modelfrequency = logfreq[fluxorder[0]] + ( ( (logfreq[fluxorder[imgnum-1]] - logfreq[fluxorder[0]]) / modelres ) * j);

	gsl_fit_linear_est(tmp_modelfrequency, c0, c1, cov00, cov01, cov11, &specindex_modelflux[i][j], &specindex_modelerror[i][j]);

	// Change from the ln space needed for determining the weights (for symmetric error bar approximation, to log10 space needed for plotting.
	specindex_modelflux[i][j] = lntol10 * specindex_modelflux[i][j];
	specindex_modelerror[i][j] = lntol10 * specindex_modelerror[i][j];

	// printf("i: %d j: %d fluxorder[0]: %d logfreq[fluxorder[0]]: %.4e tmp_modelfrequency: %.4e c0: %.4e c1 : %.4f cov00: %.4e cov01 : %.4e cov11: %.4e specindex_modelflux[i][j]: %.4e specindex_modelerror[i][j]: %.4e\n", i, j, fluxorder[0], logfreq[fluxorder[0]], tmp_modelfrequency, c0, c1, cov00, cov01, cov11, specindex_modelflux[i][j], specindex_modelerror[i][j]);

      }

      //}
    
    // Find the max and min spectral indices as we go
    if (specindex[i] > *maxsi) { // This is just for setting (legacy?) scales so we dont need to store the error on min/max
      *maxsi = specindex[i];
      minsiref = i;
    }
    if (specindex[i] < *minsi) {
      *minsi = specindex[i];
      maxsiref = i;
    }
      
    
    if (printindex == 1) {
	
      printf("Spectal Index of region %d: %.2f +- %2f\n", i, specindex[i], specindexerror[i]);
    }
    
  }
  
  /*
    if (warningcount > 0) {
    fprintf(stderr,"\n*** WARNING: Regions had a flux less than the RMS %d times! This normally occurs when fixed regions are used and they are not a good match ***\n\n", warningcount);
    }
  */
  
  printf("Maximum Spectal Index: %.2f +- %.2f\n", *maxsi, specindexerror[maxsiref]);
  printf("Minimum Spectal Index: %.2f +- %.2f\n", *minsi, specindexerror[minsiref]);
  
  return 0;
  
}
