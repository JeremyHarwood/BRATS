/* Selects regions based on a minimum average flux level */


int applyfixedregions(int *imgnum, int *xdim, int *ydim, float ****flux, float **rms, int ***regionarray, float ***regflux, int *regnumber, float ***fluxerror, int averaged, float **regionlocx, float **regionlocy, int titles, int labels, int axismarks, int fixedregnum, int fixedregref, int *regionsset, float *fluxcalerror, float onsourcenoise, float ra, float dec, float cellsize, int border, float zoom, int xshift, int yshift, int smoothmaps, int scaletype, int largetxt, int export) {

#include "smoothmaps.h"

  int mapscale(int scaletype, int istart, int iend, int jstart, int jend, int border, float cellsize, int axismarks, int labels, float ra, float dec, char *maptitle, int xshift, int yshift, int largetxt, int export, int titles);


  int a, i, j, m, tmp_regnum, *reginten, *regionsize, istart, jstart, iend, jend, *colourarray, tmp_xshift, tmp_yshift, isize, jsize;
  int intenscale = 0;
  int intencount = 0;
  int regcount = 0;
  float **fracerror, zoomedx, zoomedy, *tmp_smootharray, minrange, maxrange, adj_maxrange;

  int xaxismax = xdim[fixedregnum];
  int yaxismax = ydim[fixedregnum];

  float wedgemultiplier = 1.05; // Fixed as we always want to be just over the max number of regions. This value seems to give a nice range in most cases.

  char specmaptitle[128];

 // Declaring PGPLOT prototypes

  int cpgopen(const char *device);
  void cpgclos(void);
  void cpgenv(float xmin, float xmax, float ymin, float ymax, int just, int axis);
  void cpglab(const char *xlbl, const char *ylbl, const char *toplbl);
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

  // Check both maps are the same size

  for (i=0; i<regnumber[fixedregref]; i++){

    if ( (regionlocx[fixedregref][i] >= xdim[fixedregnum]) || (regionlocy[fixedregref][i] >= ydim[fixedregnum]) ) {
      continue;
    }
    else{
      regcount++;
    }
  }

  if ( (xdim[fixedregnum] != xdim[fixedregref]) || (ydim[fixedregnum] != ydim[fixedregref]) ) {
    fprintf(stderr,"\nError: Map dimensions do not match. If the images are centrally aligned (as is usually the case), the resize images command can be used to create a set of images that meet this requirement. Any existing regions will be kept. Exiting...\n\n");
    return 404;
  }


  // Set the number of regions to be the same
  regnumber[fixedregnum] = regnumber[fixedregref];
  intencount = 1;

  // Allocating memory

  reginten = (int *)calloc(xdim[fixedregnum]*ydim[fixedregnum], sizeof(int));
  regionsize = (int *)calloc(xdim[fixedregnum]*ydim[fixedregnum], sizeof(int));
  fracerror = (float **)calloc(xdim[fixedregnum]*ydim[fixedregnum], sizeof(float *));

  regflux[fixedregnum] = (float **)calloc(imgnum[fixedregnum], sizeof(float *));
  regionarray[fixedregnum] = (int **)calloc(xdim[fixedregnum], sizeof(int *));
  fluxerror[fixedregnum] = (float **)calloc(imgnum[fixedregnum], sizeof(float *));

  // If this set already has regions, free them up
  if (regionsset[fixedregnum] == 0) {
    free(regionlocx[fixedregnum]);
    free(regionlocx[fixedregnum]);
  }

  regionlocx[fixedregnum] = (float *)calloc(xdim[fixedregref]*ydim[fixedregref], sizeof(float));
  regionlocy[fixedregnum] = (float *)calloc(ydim[fixedregref]*xdim[fixedregref], sizeof(float));

  for (a=0; a<imgnum[fixedregnum]; a++) {
    regflux[fixedregnum][a] = (float *)calloc(xdim[fixedregnum]*ydim[fixedregnum], sizeof(float));
    fluxerror[fixedregnum][a] = (float *)calloc(xdim[fixedregnum]*ydim[fixedregnum], sizeof(float));
    fracerror[a] = (float *)calloc(xdim[fixedregnum]*ydim[fixedregnum], sizeof(float));
  }


  // Copy the array from reference to target
  for (i=0; i<xdim[fixedregref]; i++){
    regionarray[fixedregnum][i] = (int *)calloc(ydim[fixedregnum], sizeof(int));

    for (j=0; j<ydim[fixedregref]; j++){
      regionarray[fixedregnum][i][j] = regionarray[fixedregref][i][j];
    }
  }


  // Copy over the location of the regions center pixel and the colours for mapping
  for (i=0; i<regnumber[fixedregref]; i++){
    regionlocx[fixedregnum][i] = regionlocx[fixedregref][i];
    regionlocy[fixedregnum][i] = regionlocy[fixedregref][i];
  }

  for (a=0; a<imgnum[fixedregnum]; a++){
    for (i=0; i<xdim[fixedregnum]; i++){
      for (j=0; j<ydim[fixedregnum]; j++){

	tmp_regnum = regionarray[fixedregnum][i][j];

	if ( (i >= xdim[fixedregref]) || (j >= ydim[fixedregref]) ) {
	  continue; // If we are outside the bounds of the map, skip the regions
	}

	if (tmp_regnum >= 0) {

	    regflux[fixedregnum][a][tmp_regnum] += flux[fixedregnum][a][i][j];
	    regionsize[tmp_regnum]++;

	    reginten[regionarray[fixedregnum][i][j]] = intencount;

	    if (intencount == 15) {
	      intencount = 1;
	    }

	    intencount++;
	}
      }
    }
  }

  // Calculate the error of each region

  for (m=1; m<=regnumber[fixedregnum]; m++){
    for (a=0; a<imgnum[fixedregnum]; a++){

      fluxerror[fixedregnum][a][m] = sqrt( pow(regflux[fixedregnum][a][m] * fluxcalerror[a], 2) + (regionsize[m] * (pow(rms[fixedregnum][a] * onsourcenoise, 2) ) ) );

      fracerror[a][m] = fluxerror[fixedregnum][a][m] / regflux[fixedregnum][a][m];
    }
  }

  // Average the flux if set to do so

  if (averaged == 1) {
    for (i=0; i<regnumber[fixedregnum]; i++){
      for (j=0; j<imgnum[fixedregnum]; j++){
	regflux[fixedregnum][j][i] /= regionsize[i];
	fluxerror[fixedregnum][j][i] = regflux[fixedregnum][j][i] * fracerror[j][i];
      }
    }  
  }

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

  if (titles == 1) { sprintf(specmaptitle,"Map of Adaptive Regions"); }
  else { sprintf(specmaptitle," ");  }

  mapscale(scaletype, istart, iend, jstart, jend, border, cellsize, axismarks, labels, ra, dec, specmaptitle, tmp_xshift, tmp_yshift, largetxt, export, titles);

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
	else if (regionarray[fixedregnum][i+(jstart-border)][j+(istart-border)] < 0) { // Set void regions to a high value to create a white background
	  tmp_smootharray[i*jsize+j] = 1e24;
	}
	else {
	  tmp_smootharray[i*jsize+j] = regionarray[fixedregnum][i+(jstart-border)][j+(istart-border)];

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

    //cpgimag(tmp_smootharray,xaxismax,yaxismax,istart,iend,jstart,jend,minrange,adj_maxrange,trans);

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
       if (regionarray[fixedregnum][i][j] == -1) {
	 intenscale = 0;
	 continue;
       }
       else if (regionarray[fixedregnum][i][j] == -2) {
	 intenscale = 0;
	 continue;
       }
       else if (regionarray[fixedregnum][i][j] == 0) {
	 intenscale = 0;
	 continue;
       }
       else {
	 intenscale = reginten[regionarray[fixedregnum][i][j]];
       }

	colourarray[0] = intenscale;
	colourarray[1] = intenscale;
	colourarray[2] = intenscale;
	colourarray[3] = intenscale;

	cpgpixl(colourarray,2,2,1,2,1,2,j+1,j+2,i+1.5,i+2.5); // Corrected for the offset between the two plotting methods
     }
   }

 }

   // Close off the plotting device so other functions can call it
   cpgclos();

  printf("\n\nTotal number of regions: %d\n\n", regnumber[fixedregnum]);


// Free up redundant memory usage

  for (a=0; a<imgnum[fixedregnum]; a++){
    free(fracerror[a]);
  }

  free(fracerror);
  free(reginten);
  free(regionsize);

  return 0;

}
