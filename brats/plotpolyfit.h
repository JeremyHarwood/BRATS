/* Plot the poly fits */

int plotpolyfit(int *regnumber, int currentset, int *imgnum, float **frequency, int **fluxorder, int plotres, float ***curvearray, float *mincurve, float *maxcurve, int poly, float userminx, float usermaxx, float userminy, float usermaxy, int autoscalex, int autoscaley, int skip, int titles, int labels, int axismarks, char *target, int extwincontrol, char *output) {

  int i, j, g, intenscale;
  float **xcurve, **ycurve, yaxismin, yaxismax, xaxismin, xaxismax, divsize;

  char xlabel[64];
  char ylabel[64];
  char maptitle[64];

  int cpgopen(const char *device);
  void cpgclos(void);
  void cpgenv(float xmin, float xmax, float ymin, float ymax, int just, int axis);
  void cpglab(const char *xlbl, const char *ylbl, const char *toplbl);
  void cpgsci(int ci);
  void cpgline(int n, const float *xpts, const float *ypts);
  void cpgscr(int ci, float cr, float cg, float cb);


  // Allocate memory to the curve plotting arrays
  xcurve = (float **)calloc(regnumber[currentset]+1,sizeof(float *));
  ycurve = (float **)calloc(regnumber[currentset]+1,sizeof(float *));

  if (autoscaley == 1) {
    yaxismin = 1e32;
    yaxismax = -1e32;
  }
  else {
    yaxismin = userminy;
    yaxismax = usermaxy;
  }

  if (autoscalex == 1) {
    xaxismin = 0.9 * frequency[currentset][fluxorder[currentset][0]];
    xaxismax = 1.1 * frequency[currentset][fluxorder[currentset][imgnum[currentset]-1]];
  }
  else {
    xaxismin = userminx;
    xaxismax = usermaxx;
  }


  divsize = ((log10(frequency[currentset][fluxorder[currentset][imgnum[currentset]-1]]) - log10(frequency[currentset][fluxorder[currentset][0]])) / plotres);
 

  // Organise the values in to approriate arrays
	 
  for(i=(skip); i<=regnumber[currentset]; i+=skip) {

    xcurve[i] = (float *)calloc(plotres,sizeof(float));
    ycurve[i] = (float *)calloc(plotres,sizeof(float));

    for (g = 0; g < plotres; g++){

      if (g == 0) {
	xcurve[i][g] = log10(frequency[currentset][fluxorder[currentset][0]]);

	for(j=0; j<=poly; j++) {
	  ycurve[i][g] += curvearray[currentset][i][j] * pow(xcurve[i][g], j);
	}
	if ( (ycurve[i][g] < yaxismin)  && (autoscaley == 1) ) {
	  yaxismin = ycurve[i][g];
	}

	if ( (ycurve[i][g] > yaxismax) && (autoscaley == 1) ){
	  yaxismax = ycurve[i][g];
	}
      }

      else {
	xcurve[i][g] = xcurve[i][g-1] + divsize;

	for(j=0; j<=poly; j++) {
	  ycurve[i][g] += curvearray[currentset][i][j] * pow(xcurve[i][g], j);
	}

	if ( (ycurve[i][g] < yaxismin) && (autoscaley ==1) ) {
	  yaxismin = ycurve[i][g];
	}

	if ( (ycurve[i][g] > yaxismax) && (autoscaley ==1) ) {
	  yaxismax = ycurve[i][g];
	}
      }
    }
  }


  // Now plot the values
  if (autoscaley == 1) {
    if (yaxismax >= 0) {
      yaxismax *= 1.1;
    }
    else {
      yaxismax *= 0.9;
    }

    if (yaxismax >= 0) {
      yaxismin  *= 0.9;
    }
    else {
      yaxismin  *= 1.1;
    }
  }

  //printf("\nMaximum: %.2f\nMinimum: %.2f \n",yaxismax, yaxismin);

  if (extwincontrol != 1) {
    cpgopen(output);
  }

  // Swap to background and forground colours to something sensible
  cpgscr(0, 1, 1, 1);
  cpgscr(1, 0, 0, 0);

  //Ensure it is reset to the default colours
  intenscale = 1;
  cpgsci(intenscale);

  //Setup the plotting enviroment

  if (titles != 1) {
    sprintf(maptitle," "); 
  }
  else {
    if (poly == 2) {
      sprintf(maptitle,"Spectral Curvature of %s (Quadradtic Best Fit)", target);
    }
    if (poly == 3) {
      sprintf(maptitle,"Spectral Curvature of %s (Cubic Best Fit)", target);
    }
    if (poly > 3) {
      sprintf(maptitle,"Spectral Curvature of %s (Best Fit %d th Order Polynomial)", target, poly);
    }
  }

  if (labels != 1) {
    sprintf(xlabel," ");
    sprintf(ylabel," ");
  }
  else {
    sprintf(xlabel,"log(Frequency / Hz)");
    sprintf(ylabel,"log(Flux / Jy)");
  }

  if (axismarks != 1) {
    cpgenv(log10(xaxismin),log10(xaxismax), yaxismin, yaxismax, 0, -1);

  }
  else{
    cpgenv(log10(xaxismin),log10(xaxismax), yaxismin, yaxismax, 0, 0);
  }

  cpglab(xlabel, ylabel, maptitle);


  for(i=(skip); i<=regnumber[currentset]; i+=skip) {

    if (intenscale > 15) {
      intenscale = 1;
    }

    cpgsci(intenscale);
  

    cpgline(plotres, xcurve[i], ycurve[i]);

    intenscale ++;

  }

  if (extwincontrol != 1) {
    cpgclos();
  }

  // Free up the memory

  for (i=(skip); i<=regnumber[currentset]; i+=skip) {
    free(xcurve[i]);
    free(ycurve[i]);
  }

  free(xcurve);
  free(ycurve);

  return 0;

}
