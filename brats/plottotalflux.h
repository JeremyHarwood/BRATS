/* 
   Plots the total flux of the source as a function of frequency
*/


int plottotalflux(float **regflux, float *frequency, int *fluxorder, float *rms, int imgnum, int numregions, int symbol, char *target, int extwincontrol, float *fluxcalerror, float bestguess, float specrange, float fluxrange, int resolution, int specindexres, int fluxresolution, int *regionsize, float onsourcenoise, char *settarget, int titles, int labels, int axismarks, char *output) {


  // Declaring variables
  int a, i, m, n, totalpixels;

  float *totalflux, maxflux, minflux, maxfreq, minfreq, averageflux, averagefreq, *modelflux, *modelfreq, powerlaw, freqsteps, rangemax, rangemin, specsteps, specindex, tempchi, pointerror, tmp_modelflux, flux_rangemin, flux_rangemax, fluxsteps, normflux, best_normflux, *tmp_frequency, *pointerrorplus, *pointerrorminus, tmp_error, nonlogflux;

  char maptitle[64];
  char xlabel[64];
  char ylabel[64];

  float bestchi = 1e34;
  float bestspec = 1e34;

  tmp_frequency = (float *)calloc(imgnum, sizeof(float));

  // PGPLOT vars
  int cpgopen(const char *device);
  void cpgclos(void);
  void cpgenv(float xmin, float xmax, float ymin, float ymax, int just, int axis);
  void cpglab(const char *xlbl, const char *ylbl, const char *toplbl);
  void cpgsci(int ci);
  void cpgline(int n, const float *xpts, const float *ypts);
  void cpgpt(int n, const float *xpts, const float *ypts, int symbol);
  void cpgerry(int n, const float *x, const float *y1, const float *y2, float t);
  void cpgscr(int ci, float cr, float cg, float cb);

  totalflux = (float *)calloc(imgnum, sizeof(float));
  
  // Determine the spectral indices between the maximum and minimum data points (maps)

  maxflux = -1e36;
  minflux = 1e36;


  printf("\nTotal flux and statistics for %s:\n\n", target);


  for (i = 0; i < imgnum; i++) {

    for(a = 1; a <= numregions; a++) {
      totalflux[i] += regflux[i][a];
    }

    printf("Totalflux[%.2e Hz]: %.2e Jy\n", frequency[i], totalflux[i]);

    if (totalflux[i] > maxflux) {
      maxflux = totalflux[i];
    }

    if (totalflux[i] < minflux) {
      minflux = totalflux[i];
    }

  }
 

  minflux *= 0.9;
  maxflux *= 1.1;

  minfreq = frequency[fluxorder[0]]*0.9;
  maxfreq = frequency[fluxorder[imgnum-1]]*1.1;
 
  minflux = log10(minflux);
  maxflux = log10(maxflux);


  minfreq = log10(minfreq);
  maxfreq = log10(maxfreq);


  // Convert values to logs
  for (i = 0; i < imgnum; i++) {
    tmp_frequency[i] = log10(frequency[i]);
    totalflux[i] = log10(totalflux[i]);
  }
  
  // Plot the fluxes as a function of frequency

if (extwincontrol != 1) {
    cpgopen(output);
  }

  // Swap to background and forground colours to something sensible
  cpgscr(0, 1, 1, 1);
  cpgscr(1, 0, 0, 0);
  
  cpgsci(1);

  if (axismarks == 1) {
    cpgenv(minfreq, maxfreq, minflux, maxflux, 0, 0);
  }
  else {
    cpgenv(minfreq, maxfreq, minflux, maxflux, 0, -1);
  }


  if (titles == 0) {
    sprintf(maptitle, " ");
  }
  else {
    sprintf(maptitle, "Total Flux vs Frequency for %s", settarget);
  }

  if (labels == 1) {
    sprintf(xlabel,"log(Frequency / Hz)");
    sprintf(ylabel,"log(Total Flux / Jy)");
  }
  else {
    sprintf(xlabel," ");
    sprintf(ylabel," ");
  }

  cpglab(xlabel, ylabel, maptitle);


  cpgpt(imgnum, tmp_frequency, totalflux, symbol);

  totalpixels = 0;

  for(a = 1; a <= numregions; a++) {
      totalpixels +=  regionsize[a];
    }


pointerrorplus = (float *)calloc(imgnum, sizeof(float));
pointerrorminus = (float *)calloc(imgnum, sizeof(float));

// Calculate the range of the errorbars in world coordinates

  for (i = 0; i < imgnum; i++) {
  
    nonlogflux = pow(10, totalflux[i]);

    tmp_error = sqrt( pow(nonlogflux * fluxcalerror[i], 2) + (totalpixels * (pow((rms[i] * onsourcenoise), 2) ) ) );

    pointerrorplus[i] = log10(nonlogflux + tmp_error);

    pointerrorminus[i] = log10(nonlogflux - tmp_error);

    //printf("i: %d nonlogflux: %.3e (log)totalflux: %.3e tmp_error: %.3e fluxcalerror: %.3f totalpixels: %d rms: %.3e onsourcenoise: %.3f pointerrorplus: %.3e pointerrorminus: %.3e\n", i, nonlogflux, totalflux[i], tmp_error, fluxcalerror, totalpixels, rms[i], onsourcenoise, pointerrorplus[i], pointerrorminus[i]);

  }

  cpgerry(imgnum, tmp_frequency, pointerrorplus, pointerrorminus, 0.5);

  free(pointerrorplus);
  free(pointerrorminus);

  averageflux = 0;
  averagefreq = 0;

  // Get the averages

  for (i=0; i < imgnum; i++) {

    averageflux += totalflux[i];
    averagefreq += tmp_frequency[i];

  }

  averageflux /= imgnum;
  averagefreq /= imgnum;

  averageflux = pow(10, averageflux);
  averagefreq = pow(10, averagefreq);


  rangemax = bestguess + specrange;
  rangemin  = bestguess - specrange;
  flux_rangemax = averageflux - (averageflux * fluxrange);
  flux_rangemin = averageflux + (averageflux * fluxrange);

  //Find the steps size across the ranges
  fluxsteps = (flux_rangemax - flux_rangemin) / fluxresolution;
  specsteps = (rangemax - rangemin) / specindexres;

  for(n=0; n<=fluxresolution; n++) {

    normflux = flux_rangemin + (fluxsteps * n);

    for(m=0; m<=specindexres; m++) {

      // Calculate the spectral index for this loop
      specindex = rangemin + (specsteps * m);

      // Reset the temporary value for each normalisation
      tempchi = 0;


      for(i=0; i<imgnum; i++) {

	tmp_modelflux = pow( averagefreq / frequency[i], specindex) * normflux;

	//pointerror = sqrt( pow(pow(10, totalflux[i]) * fluxcalerror, 2) + (totalpixels * (pow(rms[i] * onsourcenoise, 2) ) ) );

	pointerror = pow(10, totalflux[i]) * sqrt( pow(fluxcalerror[i], 2) + (totalpixels * (pow((rms[i] * onsourcenoise)/pow(10, totalflux[i]), 2) ) ) );

	//Calculate the chi squared value
	tempchi += (pow(fabs(pow(10, totalflux[i]) - tmp_modelflux)/pointerror,2));
      }

      //If the chi squared is less than the previous best, set its new value and age

      if (tempchi < bestchi) {

	bestchi = tempchi;
	bestspec = specindex;
	best_normflux = normflux;

      }
    }

  }
  
  //This -1 needs to be set to allow changing of free parameters

  printf("\nBest Reduced X^2 is %.2f with a spectral index of %.2f\n\n", ( bestchi / (imgnum -2) ), bestspec);


  powerlaw = bestspec;

  freqsteps = (frequency[fluxorder[imgnum-1]] - frequency[fluxorder[0]] ) / resolution;

  // Calulate the model fluxes

  modelflux = (float*)calloc(resolution+1, sizeof(float));
  modelfreq = (float*)calloc(resolution+1, sizeof(float));

  for (i=0; i<=resolution; i++) {

    modelfreq[i] =  frequency[fluxorder[0]]  + (freqsteps * i);
    modelflux[i] = pow( averagefreq / modelfreq[i], powerlaw) * best_normflux;


    modelfreq[i] = log10(modelfreq[i]);
    modelflux[i] = log10(modelflux[i]);

  }

  cpgsci(2);
  cpgline(resolution+1, modelfreq,  modelflux);
  cpgsci(1);


if (extwincontrol != 1) {
    cpgclos();
  }


return 0;

}
