/* Selects regions based on a minimum average flux level */

#include "mapscale.h"
#include "drawbeam.h"

int adaptiveregions(int *imgnum, int *xdim, int *ydim, float ****flux, float minflux, float **rms, int maxarea, int ***regionarray, float ***regflux, int *regnumber, float ***fluxerror, int averaged, float **regionlocx, float **regionlocy, int **fluxorder, int currentset, float m2m, float hotpixels, float sigmalevel, int titles, int labels, int axismarks, int **regionsize, float **regmaxflux, float **regminflux, float *fluxcalerror, float onsourcenoise, int regionsset, int border, float zoom, float ra, float dec, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, int xshift, int yshift, int smoothmaps, int largetxt) {

#include "smoothmaps.h"

  int mapscale(int scaletype, int istart, int iend, int jstart, int jend, int border, float cellsize, int axismarks, int labels, float ra, float dec, char *maptitle, int xshift, int yshift, int largetxt, int export, int titles);

  int drawbeam(float xpos, float ypos, float bmaj, float bmin, float bpa);

  int a, i, b, j, x, y, e, f, m, n, z, count, loop, plus, countsinregion, *reginten, countflux, fluxloop, fluxplus, surroundingcount, surroundloop, tmp_regcount, tmp_regsize, searchtotal, istart, jstart, iend, jend, *colourarray, tmp_xshift, tmp_yshift, isize, jsize;

  int intencount = 0;
  int printreg = -1;
  int xaxismax = xdim[currentset];
  int yaxismax = ydim[currentset];
  int maxregsize = -1e9;
  int minregsize = 1e9;

  char specmaptitle[128];
  float intenscale = 0.0;
  float zoomedx, zoomedy, *tmp_smootharray, minrange, maxrange, adj_maxrange;

  int mapnumber, maxcount, maptomap, exceptions, searchcounter;
  float testflux, uppersurround, lowersurround, **fracerror, m2mupper, m2mlower, surroundingflux;

  float wedgemultiplier = 1.05; // Fixed as we always want to be just over the max number of regions. This value seems to give a nice range in most cases.


  // Declaring PGPLOT prototypes

  int cpgopen(const char *device);
  void cpgclos(void);
  void cpgpt1(float xpt, float ypt, int symbol);
  void cpgsci(int ci);
  void cpgscr(int ci, float cr, float cg, float cb);
  void cpgpixl(const int *ia, int idim, int jdim, int i1, int i2, int j1, int j2, float x1, float x2, float y1, float y2);
  void cpgimag(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float a1, float a2, const float *tr);
  void cpgctab(const float *l, const float *r, const float *g, const float *b, int nc, float contra, float bright);
  void cpgwedg(const char *side, float disp, float width, float fg, float bg, const char *label);
  void cpgcont(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, const float *c, int nc, const float *tr);
  void cpgsls(int ls);
  void cpgmtxt(const char *side, float disp, float coord, float fjust, const char *text);
  void cpgsch(float size);

  // If we are not using smooth maps then disable shifting due to a bug
  if (smoothmaps == 1) {
    tmp_xshift = xshift;
    tmp_yshift = yshift;
  }
  else {
    tmp_xshift = 0;
    tmp_yshift = 0;
  }

  //Do the memory allocation

  regnumber[currentset] = 0;
  intencount = 1;

  reginten = (int *)calloc(xdim[currentset]*ydim[currentset], sizeof(int));
  fracerror = (float **)calloc(xdim[currentset]*ydim[currentset], sizeof(float *));

  regflux[currentset] = (float **)calloc(imgnum[currentset], sizeof(float *));
  regionarray[currentset] = (int **)calloc(xdim[currentset]+1, sizeof(int *));
  fluxerror[currentset] = (float **)calloc(imgnum[currentset], sizeof(float *));
  regionlocx[currentset] = (float *)calloc((xdim[currentset]*ydim[currentset])+1, sizeof(float));
  regionlocy[currentset] = (float *)calloc((ydim[currentset]*xdim[currentset])+1, sizeof(float));
  regionsize[currentset] = (int *)calloc((xdim[currentset]*ydim[currentset])+1, sizeof(int));
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

  // Set if we should flag large map to map variations

  if (m2m <= 0) {

    maptomap = 0;
    printf("Map to map variation flagging is set to off\n");
  }
  else if ( (m2m > 0.8) || (m2m < 0.2) ) {

    printf("*** Warning: Map to map variation flagging is set to a very high or low value (suggested range 0.2 - 0.8)!!! ***\nContinuing...\n");
    maptomap = 1;
    m2mupper = 1 + m2m;
    m2mlower = 1 - m2m;
  }

  else {
    printf("Map to map variation flagging is set to on\n");
    maptomap = 1;
    m2mupper = 1 + m2m;
    m2mlower = 1 - m2m;
  }

  // Set if we should use surrounding pixel technique to flag hot pixels

  if (hotpixels <= 0.01) {

    exceptions = 0;
    printf("Surrounding pixel technique will not be used to detect hot and cold pixels\n");
  }
  else if ((hotpixels > 0.8) || (hotpixels < 0.2)  ) {

    printf("*** Warning: Hot pixel detection multiplier is set to a very high or low value (suggested range 0.2 - 0.8)!!! ***\nContinuing...\n");
    exceptions = 1;
    uppersurround = 1 + hotpixels;
    lowersurround = 1 - hotpixels;

  }
  else {
    printf("Surrounding pixel technique for hot pixel detection is set to on\n");
    exceptions = 1;
    uppersurround = 1 + hotpixels;
    lowersurround = 1 - hotpixels;
  }


  printf("Starting to make adaptive regions...\n");

  // Allocate the array memory, zero and take out anything lower than the RMS

  for (i=0; i<xdim[currentset]; i++){
    regionarray[currentset][i] = (int *)calloc(ydim[currentset]+1, sizeof(int));

    for (j=0; j<ydim[currentset]; j++){

      for (a=0; a<imgnum[currentset]; a++){

	if (flux[currentset][a][i][j] <= sigmalevel * rms[currentset][a]) {
	  regionarray[currentset][i][j] = -1;
	  break;
	}
	else if (a == (imgnum[currentset] -1) ) {
	  regionarray[currentset][i][j] = 0;
	}
      }
    }
  }


  // Go through the maps, find the lowest flux map for each region, and remove any unusable pixels

  for (i=0; i<xdim[currentset]; i++){

    for (j=0; j<ydim[currentset]; j++){

      // This skips pixels flagged as below the RMS (-1)
      if (regionarray[currentset][i][j] < 0) {
	;
      }

      else if (regionarray[currentset][i][j] == 0) {

	maxcount = 0;
	mapnumber = 0;
	 
	for (a=0; a<imgnum[currentset]; a++) {

	  testflux = 0.0;
	  loop = 1;
	  count = 1;
	  plus = 1;
	  //stopwhileloop = 0;
	  countsinregion = 1;
	  searchcounter = 1;
	  x = i;
	  y = j;

	  // This one breaks out if a certain map is bad
	  if (regionarray[currentset][i][j] != 0) {
	    break;
	  }
	  else if (regionarray[currentset][i][j] == 0) {

	    if (maptomap == 1) {
	      for (z=1;z<imgnum[currentset];z++) {
	      
		if (flux[currentset][fluxorder[currentset][z]][i][j] > (flux[currentset][fluxorder[currentset][z-1]][i][j]*m2mupper) || ( flux[currentset][fluxorder[currentset][z]][i][j] < (flux[currentset][fluxorder[currentset][z-1]][i][j]*m2mlower) ) ) {
		  regionarray[currentset][i][j] = -2;
		  break;
		}
	      }
	    }
	    
	    while ( testflux < minflux * (sqrt(countsinregion) * rms[currentset][a] * onsourcenoise) ) {

	      if (searchcounter > maxarea) {
		regionarray[currentset][i][j] = -2;
		break;
	      }

	      // See if this pixels value is an exception
	      
	      if (exceptions == 1) {

		//surroundplus = 1;
		surroundingcount = 0;
		surroundingflux = 0;
		m=(x-1);
		n=(y+1);

		for (surroundloop = 0; surroundloop < 8; surroundloop++) {

		  if (surroundloop < 2) {
		    if ((m < xdim[currentset]) && (n < ydim[currentset]) && (m >= 0) && (n >= 0)) {
		      if ((flux[currentset][a][m][n] < (flux[currentset][a][m][n]*uppersurround)) && (flux[currentset][a][m][n] > (flux[currentset][a][m][n] * lowersurround))) {
			surroundingflux +=  flux[currentset][a][m][n];
			surroundingcount++;
			m++;
		      }
		    }
		  }

		  if ((surroundloop >= 2) && (surroundloop < 4))  {
		    if ((m < xdim[currentset]) && (n < ydim[currentset]) && (m >= 0) && (n >= 0)) {
		      if ((flux[currentset][a][m][n] < (flux[currentset][a][m][n]*uppersurround)) && (flux[currentset][a][m][n] > (flux[currentset][a][m][n] * lowersurround))) {
			surroundingflux +=  flux[currentset][a][m][n];
			surroundingcount++;
			n++;
		      }
		    }
		  }
		  
		  if ((surroundloop >= 4) && (surroundloop < 6))  {
		    if ((m < xdim[a]) && (n < ydim[a]) && (m >= 0) && (n >= 0)) {
		      if ((flux[currentset][a][m][n] < (flux[currentset][a][m][n]*uppersurround)) && (flux[currentset][a][m][n] > (flux[currentset][a][m][n] * lowersurround))) {
			surroundingflux += flux[currentset][a][m][n];
			surroundingcount++;
			m--;
		      }
		    }
		  }

		  if ((surroundloop >= 6) && (surroundloop < 8))  {
		    if ((m < xdim[currentset]) && (n < ydim[currentset]) && (m >= 0) && (n >= 0)) {
		      if ((flux[currentset][a][m][n] < (flux[currentset][a][m][n]*uppersurround)) && (flux[currentset][a][m][n] > (flux[currentset][a][m][n] * lowersurround ))) {
			surroundingflux +=  flux[currentset][a][m][n];
			surroundingcount++;
			m--;
		      }
		    }
		  }
		}
	      
		// If there arent enough surrounding values to check the validity, assume it is bad
		if ((x < xdim[currentset]) && (y < ydim[currentset]) && (x >= 0) && (y >= 0)) {

		  if ((surroundingcount <= 3) && (regionarray[currentset][x][y] == 0)) {
		    regionarray[currentset][x][y] = -2;
		  }
		  else if(regionarray[currentset][x][y] == 0)  {
		    surroundingflux /= surroundingcount;
		  }

		  // If the pixel a flux of a set amount greater those surrounding it, mark it as bad

		  /*if (flux[a][x][y] > (surroundmultiplier*surroundingflux)) {
		    regionarray[x][y] = -2;
		    }
		  */
		}
	      }

	      if ( testflux < (minflux *  (sqrt(countsinregion) * rms[currentset][a] * onsourcenoise) ) && (count <= loop) ) {
		if ((x < xdim[currentset]) && (y < ydim[currentset]) && (x >= 0) && (y >= 0)) {
		  if (regionarray[currentset][x][y] == 0) {
		    testflux += flux[currentset][a][x][y];
		    countsinregion++;
		  }
		}
		if (plus == 1) {
		  x++;
		}
		else {
		  x--;
		}
	      }

	      else if ( testflux < (minflux * (sqrt(countsinregion) * rms[currentset][a] * onsourcenoise) ) && (count > loop) && (count <= (2*loop) ) ) {
		if ((x < xdim[currentset]) && (y < ydim[currentset]) && (x >= 0) && (y >= 0)) {
		  if (regionarray[currentset][x][y] == 0) {
		    testflux += flux[currentset][a][x][y];
		    countsinregion++;
		  }
		}
		if (plus == 1) {
		  y++;
		}
		else {
		  y--;
		}
	      }
	    
	      if (count == (2*loop)) {
		loop++;
		count = 0;
		if (plus == 1) {
		  plus = 0;
		}
		else {
		  plus = 1;
		}
	      }

	      if ( testflux < minflux * (sqrt(countsinregion) * rms[currentset][a] * onsourcenoise) ) {
		if (countsinregion > maxcount) {
		  maxcount = countsinregion;
		  mapnumber = a;
		}
	      }
	      // These need to be at the end of the while loop to prevent a division by 0
	      count++;
	      searchcounter++;

	    }
	  }

	  // If all maps reach the minimum flux, set the region and get the flux
	  if (regionarray[currentset][i][j] != 0) {
	    ;
	  }

	  else if (regionarray[currentset][i][j] == 0) {

	    countflux = 0;
	    fluxloop = 1;
	    fluxplus = 1;
	    tmp_regcount = 1;

	  
	    if ((a == (imgnum[currentset]-1)) && (regionarray[currentset][i][j] == 0)) {
	      regnumber[currentset]++;
	      regionsize[currentset][regnumber[currentset]] = 0;
	      e = i;
	      f = j;

	      //fluxaddedcount = 0;

	      regionlocx[currentset][regnumber[currentset]] = i;
	      regionlocy[currentset][regnumber[currentset]] = j;

	      reginten[regnumber[currentset]] = intencount;

	      if (intencount == 15) {
		intencount = 1;
	      }

	      intencount++;

	      searchtotal = 0;
	    
	      while ( regflux[currentset][mapnumber][regnumber[currentset]] < minflux * (sqrt(tmp_regcount) * rms[currentset][a] * onsourcenoise) ) {

		// Break if we go over the maximum search area
		if (searchtotal >= maxarea) {
		  break;
		}

		searchtotal++;
		countflux++;
	     
		//printf("Region: %d tmp_regcount = %d\n", regnumber[currentset], tmp_regcount);
	      
		if (countflux <= fluxloop) {

		  if ((e < xdim[currentset]) && (f < ydim[currentset]) && (e > 0) && (f > 0)) {

		    if (regionarray[currentset][e][f] == 0) {

		      for (b = 0; b<imgnum[currentset]; b++) {

			if ( (regionsize[currentset][regnumber[currentset]] > 0) && (b == 0) ){
			  tmp_regcount++;			}

			if (b == 0) {
			  regionsize[currentset][regnumber[currentset]]++;
			}

			regflux[currentset][b][regnumber[currentset]] += flux[currentset][b][e][f];
			regionarray[currentset][e][f] = regnumber[currentset];

			// Find the min / max fluxes for each map in the dataset when ragions have been considered

			if (regflux[currentset][b][regnumber[currentset]] > regmaxflux[currentset][b]) {
			  regmaxflux[currentset][b] = regflux[currentset][b][regnumber[currentset]];
			}
			else if (regflux[currentset][b][regnumber[currentset]] < regminflux[currentset][b]){
			  regminflux[currentset][b] = regflux[currentset][b][regnumber[currentset]];
			}

			//Outputs .reg data for testing bad regions
			if (printreg >= 0) {
			  if (b == 1) {
			    if (regnumber[currentset] == printreg) {
			      printf("box(%d,%d,1,1,0)\n", f, e);
			    }
			  }
			}
		      }
		    }
		  }
		  if (fluxplus == 1) {
		    e++;
		  }
		  else {
		    e--;
		  }
		}

		else if ((countflux > fluxloop) && (countflux <= (2*fluxloop))) {
		  if ((e < xdim[currentset]) && (f < ydim[currentset]) && (e > 0) && (f > 0)) {
		    if (regionarray[currentset][e][f] == 0) {
		      for (b=0; b<imgnum[currentset]; b++) {

			if ( (regionsize[currentset][regnumber[currentset]] > 0) && (b == 0) ) {
			  tmp_regcount++;
			}

			if (b == 0) {
			  regionsize[currentset][regnumber[currentset]]++;
			}

			regflux[currentset][b][regnumber[currentset]] += flux[currentset][b][e][f];
			regionarray[currentset][e][f] = regnumber[currentset];

			// Find the min / max fluxes for each map in the dataset when ragions have been considered
			if (regflux[currentset][b][regnumber[currentset]] > regmaxflux[currentset][b]) {
			  regmaxflux[currentset][b] = regflux[currentset][b][regnumber[currentset]];
			}
			else if (regflux[currentset][b][regnumber[currentset]] < regminflux[currentset][b]){
			  regminflux[currentset][b] = regflux[currentset][b][regnumber[currentset]];
			}


			//Outputs .reg data for testing bad regions
			if (printreg >= 0) {
			  if (b == 1) {
			    if (regnumber[currentset] == printreg) {
			      printf("box(%d,%d,1,1,0)\n", f, e);
			    }
			  }
			}
		      }
		    }
		  }

		  if (fluxplus == 1) {
		    f++;
		  }
		  else {
		    f--;
		  }
		}

		//Snake round the initial point
		if (countflux >= (2*fluxloop)) {
		  fluxloop++;
		  countflux = 0;
		  if (fluxplus == 1) {
		    fluxplus = 0;
		  }
		  else {
		    fluxplus = 1;
		  }
		}

	      }
	    }
	  }
	}
      }
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

    else if (regionsize[currentset][m] < minregsize){ 
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

	if (regflux[currentset][j][i] < regminflux[currentset][j]){
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

  // If we havent been able to set any regions, kick back to main
  if (regnumber[currentset] < 1) {
    printf("\n*** Error: unable to define any regions for the given parameters ***\n\n");
    return 404;
  }

  // Setup the plot area

  printf("Mapping the regions...");

  cpgopen("/xs");

  // Swap to background and forground colours to something sensible
  cpgscr(0, 1, 1, 1);
  cpgscr(1, 0, 0, 0);

  cpgsci(1);

  // Calculate the zoomed start and end points
  zoomedx = (xaxismax / zoom);
  zoomedy = (yaxismax / zoom);

  istart = tmp_xshift + border + ( (xaxismax - (int)zoomedx) / 2);
  jstart = tmp_yshift + border + ( (yaxismax - (int)zoomedy) / 2);

  iend = tmp_xshift + xaxismax - ( ( (xaxismax - (int)zoomedx) / 2) + border);
  jend = tmp_yshift + yaxismax - ( ( (yaxismax - (int)zoomedy) / 2) + border);

  isize = (iend-istart) + (2*border);
  jsize = (jend-jstart) + (2*border);

  //printf("\nistart: %d jstart: %d iend: %d jend: %d\n\n", istart, jstart, iend, jend);


  if (titles == 1) { sprintf(specmaptitle,"Adaptive Regions"); }
  else { sprintf(specmaptitle," ");  }

  mapscale(scaletype, istart, iend, jstart, jend, border, cellsize, axismarks, labels, ra, dec, specmaptitle, tmp_xshift, tmp_yshift, largetxt, 0, titles);

  /*
  //For testing flux values
  for (i=0; i<imgnum; i++){
  for (j=1; j<regnumber[0]; j++){

  printf("Map: %d Region: %d Flux: %.3e \n", i, j, regflux[i][j]);
  
  }
  }
  */

  // Plot the regions to make a map

  if (smoothmaps == 1)  {

    maxrange = -1e23;
      
    minrange = 0; // Fixed as we will always have a region 0
   
    tmp_smootharray = (float *)calloc(isize*jsize, sizeof(float));

    for (i=0; i<isize; i++) {
      for (j=0; j<jsize; j++) {
	if ( ((i+jstart-border) < 0) || ((j+istart-border) < 0) || ((i+jstart-border) >= xaxismax) || ((j+istart-border) >= yaxismax) ) {
	  tmp_smootharray[i*jsize+j] = 1e24;
	}
	else if (regionarray[currentset][i+(jstart-border)][j+(istart-border)] < 0) { // Set void regions to a high value to create a white background

	  tmp_smootharray[i*jsize+j] = 1e24;
	}
	else {

	  tmp_smootharray[i*jsize+j] = regionarray[currentset][i+(jstart-border)][j+(istart-border)];

	  if (tmp_smootharray[i*jsize+j] < minrange) {
	    minrange = tmp_smootharray[i*jsize+j];
	  }

	  if (tmp_smootharray[i*jsize+j] > maxrange) {
	    maxrange = tmp_smootharray[i*jsize+j];
	  }
	}
      }
    }

    adj_maxrange = maxrange * wedgemultiplier;

    if (adj_maxrange < 0.1) { // this prevents a 0 range, especially when using one region
      printf("\n*** Warning: Range is too small, setting maxrange to a value of 0.1\n\n"); 
      adj_maxrange = 0.1;
    }

    cpgctab(rl,rr,rg,rb,11,contra,bright);
    if (largetxt == 1) { cpgsch(1.30); }
    cpgwedg("RI", 1.5, 5.0, minrange, adj_maxrange, "Regions");
    if (largetxt == 1) { cpgsch(1.00); } // Revert to default

    float transshift[6]={(istart-border), 1.0, 0.0, (jstart-border), 0.0, 1.0};
    cpgimag(tmp_smootharray,isize,jsize,border,isize-border,border,jsize-border,minrange,adj_maxrange,transshift);

    free(tmp_smootharray);
  }

  else {

    colourarray = (int *)calloc(4, sizeof(int));

    for (i=istart; i<iend; i++){
      for (j=jstart; j<jend; j++){
	if ( (i >= xaxismax) || (j >= yaxismax) || (i < 0) || (j < 0) ) {
	  intenscale = 0;
	  continue;
	}
	else if (regionarray[currentset][i][j] == -1) {
	  intenscale = 0;
	  continue;
	}
	else if (regionarray[currentset][i][j] == -2) {
	  intenscale = 0;
	  continue;
	}
	else if (regionarray[currentset][i][j] == 0) {
	  intenscale = 0;
	  continue;
	}
	else {
	  intenscale = reginten[regionarray[currentset][i][j]];
	}

	colourarray[0] = intenscale;
	colourarray[1] = intenscale;
	colourarray[2] = intenscale;
	colourarray[3] = intenscale;

	//cpgpixl(colourarray,2,2,1,2,1,2,j,j+1,i,i+1);
	cpgpixl(colourarray,2,2,1,2,1,2,j+1,j+2,i+1.5,i+2.5); // Corrected for the offset between the two plotting methods
	//cpgpixl(colourarray,2,2,1,2,1,2,xcount+1,ycount+2,xcount+1.5,ycount+2.5); // Corrected for the offset between the two plotting methods
      }
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
