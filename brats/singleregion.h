/* Selects regions based on a minimum average flux level */

//#include "mapscale.h"
//#include "drawbeam.h"

int singleregion(int *imgnum, int *xdim, int *ydim, float ****flux, float minflux, float **rms, int maxarea, int ***regionarray, float ***regflux, int *regnumber, float ***fluxerror, int averaged, float **regionlocx, float **regionlocy, int **fluxorder, int currentset, float m2m, float hotpixels, float sigmalevel, int titles, int labels, int axismarks, int **regionsize, float **regmaxflux, float **regminflux, float *fluxcalerror, float onsourcenoise, int regionsset, int border, float zoom, float ra, float dec, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, int xshift, int yshift, int largetxt) {

  int mapscale(int scaletype, int istart, int iend, int jstart, int jend, int border, float cellsize, int axismarks, int labels, float ra, float dec, char *maptitle, int xshift, int yshift, int largetxt, int export, int titles);

  int drawbeam(float xpos, float ypos, float bmaj, float bmin, float bpa);

  int a, i, j, m, *reginten, tmp_regsize, istart, jstart, iend, jend, *colourarray, isize, jsize;

  //int intencount = 0;
  int xaxismax = xdim[currentset];
  int yaxismax = ydim[currentset];
  int maxregsize = -1e9;
  int minregsize = 1e9;

  char specmaptitle[128];
  float intenscale = 0.0;
  float zoomedx, zoomedy, **fracerror;


  // Declaring PGPLOT prototypes

  int cpgopen(const char *device);
  void cpgclos(void);
  void cpgpt1(float xpt, float ypt, int symbol);
  void cpgsci(int ci);
  void cpgscr(int ci, float cr, float cg, float cb);
  void cpgpixl(const int *ia, int idim, int jdim, int i1, int i2, int j1, int j2, float x1, float x2, float y1, float y2);

  //Do the memory allocation
  regnumber[currentset] = 1;
  //intencount = 1;

  reginten = (int *)calloc(xdim[currentset]*ydim[currentset], sizeof(int));
  fracerror = (float **)calloc(xdim[currentset]*ydim[currentset], sizeof(float *));

  if (regionsset == 1) {

    free(regflux[currentset]);
    free(regionarray[currentset]);
    free(fluxerror[currentset]);
    free(regionlocx[currentset]);
    free(regionlocy[currentset]);
    free(regionsize[currentset]);
    free(regmaxflux[currentset]);
    free(regminflux[currentset]);

    regionsset = 0;

  }

  regflux[currentset] = (float **)calloc(imgnum[currentset], sizeof(float *));
  regionarray[currentset] = (int **)calloc(xdim[currentset], sizeof(int *));
  fluxerror[currentset] = (float **)calloc(imgnum[currentset], sizeof(float *));
  regionlocx[currentset] = (float *)calloc(xdim[currentset]*ydim[currentset], sizeof(float));
  regionlocy[currentset] = (float *)calloc(ydim[currentset]*xdim[currentset], sizeof(float));
  regionsize[currentset] = (int *)calloc(xdim[currentset]*ydim[currentset], sizeof(int));
  regmaxflux[currentset] = (float *)calloc(imgnum[currentset], sizeof(float));
  regminflux[currentset] = (float *)calloc(imgnum[currentset], sizeof(float));

  for (a=0; a<imgnum[currentset]; a++) {
    regflux[currentset][a] = (float *)calloc(xdim[currentset]*ydim[currentset], sizeof(float));
    fluxerror[currentset][a] = (float *)calloc(xdim[currentset]*ydim[currentset], sizeof(float));
    fracerror[a] = (float *)calloc(xdim[currentset]*ydim[currentset], sizeof(float));
    regmaxflux[currentset][a] = -1e+32;
    regminflux[currentset][a] = 1e+32;
  }
  

  printf("Memory assigned...\n");

  printf("Starting to make the region...\n");

  // Allocate the array memory, zero and take out anything lower than the RMS

  for (i=0; i<xdim[currentset]; i++){
    regionarray[currentset][i] = (int *)calloc(ydim[currentset], sizeof(int));

    for (j=0; j<ydim[currentset]; j++){

      for (a=0; a<imgnum[currentset]; a++){

	if (flux[currentset][a][i][j] <= sigmalevel * rms[currentset][a] ) {

	  regionarray[currentset][i][j] = -1;
	  break;
	}

	else if (a == (imgnum[currentset] -1) ) {
	  regionarray[currentset][i][j] = 0;
	}
      }
    }
  }

  // Get the flux for the whole regions so long as it has not been flagged below RMS

  for (i=0; i<xdim[currentset]; i++){

    for (j=0; j<ydim[currentset]; j++){

      // This skips pixels flagged as below the RMS (-2)
      if (regionarray[currentset][i][j] < 0) {
	;
      }

      else if (regionarray[currentset][i][j] == 0) {
	 
	for (a=0; a<imgnum[currentset]; a++) {

	    regflux[currentset][a][regnumber[currentset]] += flux[currentset][a][i][j];

	    if ((a == (imgnum[currentset]-1)) && (regionarray[currentset][i][j] == 0)) {

	      regionlocx[currentset][regnumber[currentset]] = i;
	      regionlocy[currentset][regnumber[currentset]] = j;
	      regionarray[currentset][i][j] = regnumber[currentset];
	      reginten[regnumber[currentset]] = 2;
	      regionsize[currentset][regnumber[currentset]]++;


	    }
	    
	}
      }
    }
  }

  // Get the min / max flux values for the region
  for (a=0; a<imgnum[currentset]; a++) {
    if (regflux[currentset][a][regnumber[currentset]] > regmaxflux[currentset][a]) {
      regmaxflux[currentset][a] = regflux[currentset][a][regnumber[currentset]];
    }
    if (regflux[currentset][a][regnumber[currentset]] < regminflux[currentset][a]){
      regminflux[currentset][a] = regflux[currentset][a][regnumber[currentset]];
    }
  }


  // Calculate the error of each region

  for (m=1; m<=regnumber[currentset]; m++){

    for (a=0; a<imgnum[currentset]; a++){

      fluxerror[currentset][a][m] = sqrt( pow(regflux[currentset][a][m] * fluxcalerror[a], 2) + (regionsize[currentset][m] * (pow(rms[currentset][a] * onsourcenoise, 2) ) ) );
      
      fracerror[a][m] = fluxerror[currentset][a][m] / regflux[currentset][a][m];
    }

    // Find the max and min regions sizes while we are at it 
    if (regionsize[currentset][m] > maxregsize) {
      maxregsize = regionsize[currentset][m];
    }

    if (regionsize[currentset][m] < minregsize){ 
      minregsize = regionsize[currentset][m];
    }

  }
	 
  //Average the flux if set to do so

  if (averaged == 1) {

    // Reset the max and min values
    for (a=0; a<imgnum[currentset]; a++) {
      regmaxflux[currentset][a] = -1e32;
      regminflux[currentset][a] = 1e32;
    }

    // Average the fluxes
    for (i=1; i<=regnumber[currentset]; i++){

      for (j=0; j<imgnum[currentset]; j++){

	if (regionsize[currentset][i] <= 0) {
	  tmp_regsize = 1; 
	}
	else{
	  tmp_regsize = regionsize[currentset][i];
	}

	regflux[currentset][j][i] /= tmp_regsize;
	fluxerror[currentset][j][i] = regflux[currentset][j][i] * fracerror[j][i];

	// Find the new max and mix values as we go
	if (regflux[currentset][j][i] > regmaxflux[currentset][j]) {
	  regmaxflux[currentset][j] = regflux[currentset][j][i];
	}

	if (regflux[currentset][j][i] < regminflux[currentset][j]) {
	  regminflux[currentset][j] = regflux[currentset][j][i];
	}
      }
    }
  }



  // Get the fluxes for each region (Testing only)
  /*
    for (i=0; i<regnumber[currentset]; i++){

    for (j=0; j<imgnum[currentset]; j++){

    printf("regflux[imgnum: %d][regnum: %d] %.2e\n", j, i, regflux[j][i]);

    }
    }
  */

  // Setup the plot area

  printf("Mapping the regions...");

  cpgopen("/xs");

  // Swap to background and forground colours to something sensible
  cpgscr(0, 1, 1, 1);
  cpgscr(1, 0, 0, 0);

  cpgsci(1);

  // Calculate the zoomed start and end points
  zoomedx = xaxismax / zoom;
  zoomedy = yaxismax / zoom;

  istart = xshift + border + ( (xaxismax - (int)zoomedx) / 2);
  jstart = yshift + border + ( (yaxismax - (int)zoomedy) / 2);

  iend = xshift + xaxismax - ( ( (xaxismax - (int)zoomedx) / 2) + border);
  jend = yshift + yaxismax - ( ( (yaxismax - (int)zoomedy) / 2) + border);

  isize = (iend-istart) + (2*border);
  jsize = (jend-jstart) + (2*border);


  if (titles == 1) { sprintf(specmaptitle,"Adaptive Regions"); }
  else { sprintf(specmaptitle," ");  }

  mapscale(scaletype, istart, iend, jstart, jend, border, cellsize, axismarks, labels, ra, dec, specmaptitle, xshift, yshift, largetxt, 0, titles);

  /*
  //For testing flux values
  for (i=0; i<imgnum; i++){
  for (j=1; j<regnumber[0]; j++){

  printf("Map: %d Region: %d Flux: %.3e \n", i, j, regflux[i][j]);
  
  }
  }
  */

  // Plot the regions to make a map

  colourarray = (int *)calloc(2, sizeof(int));

  for (i=istart; i<iend; i++){
    for (j=jstart; j<jend; j++){
      if (regionarray[currentset][i][j] == -1) {
	intenscale = 0;
      }
      else if (regionarray[currentset][i][j] == -2) {
	intenscale = 0;
      }
      else if (regionarray[currentset][i][j] == 0) {
	intenscale = 0;
      }
      else {
	intenscale = reginten[regionarray[currentset][i][j]];
      }

      colourarray[0] = intenscale;
      colourarray[1] = intenscale;
      colourarray[2] = intenscale;
      colourarray[3] = intenscale;

      cpgpixl(colourarray,2,2,1,2,1,2,j,j+1,i,i+1);
    }
  }

  // Plot the beam if the option is set to on
  if (plotbeam == 1) {
    double tmp_beamposx = ((double)beamposx)/100;
    double tmp_beamposy = (double)beamposy/100;
    int beamposi = (istart+0.5-border) + (tmp_beamposx*(isize-(bmaj/cellsize))) + ((bmaj/cellsize)/2);
    int beamposj = (jstart+0.5-border) + (tmp_beamposy*(jsize-(bmaj/cellsize))) + ((bmaj/cellsize)/2);

    // This must be after the main plotting to create an overlay
    drawbeam(beamposi, beamposj, bmaj/cellsize, bmin/cellsize, bpa);
  }

  // Close off the plotting device so other functions can call it
  cpgclos();

  // Output some info

  printf("\n\nTotal number of regions: %d\n", regnumber[currentset]);
  printf("Largest Region: %d Pixels\n", maxregsize);
  printf("Smallest Region: %d Pixels\n\n", minregsize);


  // Free up redundant memory usage

  for (a=0; a<imgnum[currentset]; a++){
    free(fracerror[a]);
  }

  free(fracerror);
  free(reginten);

  return 0;

}
