/* Plots flux against frequency */


int plotfluxbyregion(float ***regflux, int *imgnum, int **fluxorder, float **frequency, float minfreq, float maxfreq, float **minflux, float **maxflux, int *regnumber, float ***fluxerror, char maptitle[], int titles, int labels, int axismarks, int extwincontrol, int currentset, int uselogx, int uselogy, int autoscalex, int autoscaley, float userminx, float usermaxx, float userminy, float usermaxy, int printincreasing, int skip, int symbol, char *output) {


  int cpgopen(const char *device);
  void cpgclos(void);
  void cpgenv(float xmin, float xmax, float ymin, float ymax, int just, int axis);
  void cpglab(const char *xlbl, const char *ylbl, const char *toplbl);
  void cpgpt1(float xpt, float ypt, int symbol);
  void cpgsci(int ci);
  void cpgask(int flag);
  void cpgscr(int ci, float cr, float cg, float cb);
  void cpgline(int n, const float *xpts, const float *ypts);
  void cpgerry(int n, const float *x, const float *y1, const float *y2, float t);


  int a, b, i, intenscale;

  float *fluxerrorplus, *fluxerrorminus, *plotarray, *frequencyarray;

  double xaxismin, xaxismax;

  double yaxismin = 1e+32;
  double yaxismax = -1e+32;

  char xlabel[64];
  char ylabel[64];


  // Find the overall min and max values

  if (autoscalex == 1) {
    xaxismin = 0.9 * minfreq;
    xaxismax = 1.1 * maxfreq;
  }
  else {
    xaxismin = userminx;
    xaxismax = usermaxx;
  }
  
  // Find the new min max for the autoscale if using skip else save time and use the total min max values
  if (autoscaley == 1) {
      for (a=0; a<imgnum[currentset]; a++){
	for (b=skip; b<=regnumber[currentset]; b+=skip){
	  if (regflux[currentset][fluxorder[currentset][a]][b] > yaxismax) {
	    yaxismax = regflux[currentset][fluxorder[currentset][a]][b];
	  }
	  if (regflux[currentset][fluxorder[currentset][a]][b] < yaxismin) {
	    yaxismin = regflux[currentset][fluxorder[currentset][a]][b];
	  }
	}
      }


    yaxismin *= 0.7;
    yaxismax *= 1.3;
  }

  // Else use the user defined min max values

  else {
    yaxismin = userminy;
    yaxismax = usermaxy;
  }

  if (extwincontrol == 1) {
    cpgask(1);
  }
  else {
    cpgopen(output);
  }

 // Swap to background and forground colours to something sensible
  cpgscr(0, 1, 1, 1);
  cpgscr(1, 0, 0, 0);


  // Setup the plotting area

  cpgsci(1);

  if ( (uselogx == 1) && (uselogy == 1) && (axismarks == 1) ) {
    cpgenv(log10(xaxismin), log10(xaxismax), log10(yaxismin), log10(yaxismax), 0, 0);
  }

  else if ( (uselogx == 1) && (uselogy == 1) && (axismarks != 1) ) {
    cpgenv(log10(xaxismin), log10(xaxismax), log10(yaxismin), log10(yaxismax), 0, -1);
  }

  else if ( (uselogx == 1) && (uselogy == 0) && (axismarks == 1) ) {
    cpgenv(log10(xaxismin), log10(xaxismax), yaxismin, yaxismax, 0, 0);
  }

  else if ( (uselogx == 1) && (uselogy == 0) && (axismarks != 1) ) {
    cpgenv(log10(xaxismin), log10(xaxismax), yaxismin, yaxismax, 0, -1);
  }

  else if ( (uselogx == 0) && (uselogy == 1) && (axismarks == 1) ) {
    cpgenv(xaxismin, xaxismax, log10(yaxismin), log10(yaxismax), 0, 0);
  }

  else if ( (uselogx == 0) && (uselogy == 1) && (axismarks != 1) ) {
    cpgenv(xaxismin, xaxismax, log10(yaxismin), log10(yaxismax), 0, -1);
  }

  else if ( (uselogx == 0) && (uselogy == 0) && (axismarks == 1) ) {
    cpgenv(xaxismin,xaxismax,yaxismin,yaxismax, 0, 0);
  }

  else if ((uselogx == 0) && (uselogy == 0) && (axismarks != 1) ) {
    cpgenv(xaxismin,xaxismax,yaxismin,yaxismax, 0, -1);
  }

  else {
    fprintf(stderr,"Error: Unhandled value of uselogx and uselogy combination or axismarks value! Exiting command...\n");
    return 4;
  }

  if (titles == 0) { sprintf(maptitle," ");  }


  if (labels != 1) {
    sprintf(xlabel," ");
    sprintf(ylabel," ");
  }
  else if ( (uselogx == 1) && (uselogy == 1) && (labels == 1) ) {
    sprintf(xlabel,"log(Frequency / Hz)");
    sprintf(ylabel,"log(Flux / Jy)");
  }
  else if ( (uselogx == 1) && (uselogy == 0)  && (labels == 1) ) {
    sprintf(xlabel,"log(Frequency / Hz)");
    sprintf(ylabel,"Flux / Jy");
  }

  else if ( (uselogx == 0) && (uselogy == 1)  && (labels == 1) ) {
    sprintf(xlabel,"Frequency / Hz");
    sprintf(ylabel,"log(Flux / Jy)");
  }

  else if ( (uselogx == 0) && (uselogy == 0) && (labels == 1) ) {
    sprintf(xlabel,"Frequency (Hz)");
    sprintf(ylabel,"Flux (Jy)");
  }

  else {
    fprintf(stderr,"Error: Unhandled value of uselogx and uselogy combination or labels values! Exiting command...\n");
    return 4;
  }

  cpglab(xlabel, ylabel, maptitle);

  // Start the plotting

  intenscale = 0.0;

  plotarray = (float *)calloc(imgnum[currentset],sizeof(float));
  frequencyarray = (float *)calloc(imgnum[currentset],sizeof(float));
  fluxerrorplus = (float *)calloc(imgnum[currentset],sizeof(float));
  fluxerrorminus  = (float *)calloc(imgnum[currentset],sizeof(float));

  // Get the frequency in order to plot


  for (b=skip; b<=regnumber[currentset]; b+=skip){
    intenscale++;

    // Setup the line colours

    if (intenscale > 15) {
      intenscale = 1;
    }

    cpgsci(intenscale);


    for (a=0; a<imgnum[currentset]; a++) {

      if ( (uselogx == 1) && (uselogy == 1) ) {

	if ( (a == 0) && (b == skip) ) {

	  for (i=0; i<imgnum[currentset]; i++) {

	    frequencyarray[i] = log10(frequency[currentset][fluxorder[currentset][i]]);
	  }
	}

	plotarray[a] = log10(regflux[currentset][fluxorder[currentset][a]][b]);

	fluxerrorplus[a] = log10(regflux[currentset][fluxorder[currentset][a]][b] + fluxerror[currentset][fluxorder[currentset][a]][b]);

	fluxerrorminus[a] = log10(regflux[currentset][fluxorder[currentset][a]][b] - fluxerror[currentset][fluxorder[currentset][a]][b]);


	if ( (printincreasing == 1) && (a > 0 ) ) {

	  if (plotarray[a - 1]  < plotarray[a] ) {
	    printf("Warning: Region %d has an increasing flux! (%.3e vs %.3e) in maps %d and %d\n", b, regflux[currentset][fluxorder[currentset][a-1]][b], regflux[currentset][fluxorder[currentset][a]][b], a-1, a);
	  }
	}

      }
   

      else if ( (uselogx == 0) && (uselogy == 1) ) {
	   

	if ( (a == 0) && (b == skip) ) {

	  for (i=0; i<imgnum[currentset]; i++) {

	    frequencyarray[i] = frequency[currentset][fluxorder[currentset][i]];
	  }
	}

	plotarray[a] = log10(regflux[currentset][fluxorder[currentset][a]][b]);

	fluxerrorplus[a] = log10(regflux[currentset][fluxorder[currentset][a]][b] + fluxerror[currentset][fluxorder[currentset][a]][b]);

	fluxerrorminus[a] = log10(regflux[currentset][fluxorder[currentset][a]][b] - fluxerror[currentset][fluxorder[currentset][a]][b]);


	if ( (printincreasing == 1) && (a > 0 ) ) {

	  if (plotarray[a - 1]  < plotarray[a] ) {
	    printf("Warning: Region %d has an increasing flux! (%.3e vs %.3e) in maps %d and %d\n", b, regflux[currentset][fluxorder[currentset][a-1]][b], regflux[currentset][fluxorder[currentset][a]][b], a-1, a);
	  }
	}

      }


      else if ( (uselogx == 1) && (uselogy == 0) ){
	
	if ( (a == 0) && (b == skip) ) {

	  for (i=0; i<imgnum[currentset]; i++) {

	    frequencyarray[i] = log10(frequency[currentset][fluxorder[currentset][i]]);
	  }
	}

	plotarray[a] = regflux[currentset][fluxorder[currentset][a]][b];

	fluxerrorplus[a] = regflux[currentset][fluxorder[currentset][a]][b] + fluxerror[currentset][fluxorder[currentset][a]][b];

	fluxerrorminus[a] = regflux[currentset][fluxorder[currentset][a]][b] - fluxerror[currentset][fluxorder[currentset][a]][b];


	if ( (printincreasing == 1) && (a > 0 ) ) {

	  if (plotarray[a - 1]  < plotarray[a] ) {
	    printf("Warning: Region %d has an increasing flux! (%.3e vs %.3e) in maps %d and %d\n", b, regflux[currentset][fluxorder[currentset][a-1]][b], regflux[currentset][fluxorder[currentset][a]][b], a-1, a);
	  }
	}

      }

      else {

	
	if ( (a == 0) && (b == skip) ) {

	  for (i=0; i<imgnum[currentset]; i++) {

	    frequencyarray[i] = frequency[currentset][fluxorder[currentset][i]];
	  }
	}

	plotarray[a] = regflux[currentset][fluxorder[currentset][a]][b];

	fluxerrorplus[a] = regflux[currentset][fluxorder[currentset][a]][b] + fluxerror[currentset][fluxorder[currentset][a]][b];

	fluxerrorminus[a] = regflux[currentset][fluxorder[currentset][a]][b] - fluxerror[currentset][fluxorder[currentset][a]][b];


	if ( (printincreasing == 1) && (a > 0 ) ) {

	  if (plotarray[a - 1]  < plotarray[a] ) {
	    printf("Warning: Region %d has an increasing flux! (%.3e vs %.3e) in maps %d and %d\n", b, regflux[currentset][fluxorder[currentset][a-1]][b], regflux[currentset][fluxorder[currentset][a]][b], a-1, a);
	  }
	}

      }

    }


    cpgline(imgnum[currentset], frequencyarray, plotarray);  
    cpgerry(imgnum[currentset], frequencyarray, fluxerrorplus, fluxerrorminus, 0.5);

  }

  if (extwincontrol != 1) {
    cpgclos();
  }


  return 0;

}
