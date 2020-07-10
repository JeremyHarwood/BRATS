// Plot the spectral index error map of a given data set

int specindexerrormap(float **specindexerror, int numregions, int *xdim, int *ydim, char *settarget, float *minsi, float *maxsi, int ***regionarray, int currentset, char *maptitle, float zoom, int border, int titles, int labels, int axismarks, int symbol, int smoothmaps, float **flux, int extwincontrol, int contours, int contourlevels, int contourcolour, int contourlinestyle, int firstcontour, int contourlogs, char *output, float ra, float dec, float delt1, float delt2, float crpix1, float crpix2, float crota1, float crota2, float eq, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, float wedgemultiplier, int xshift, int yshift, float minfreq, float maxfreq, int export, int exportfits, int largetxt) {

#include "smoothmaps.h"

  int mapscale(int scaletype, int istart, int iend, int jstart, int jend, int border, float cellsize, int axismarks, int labels, float ra, float dec, char *maptitle, int xshift, int yshift, int largetxt, int export, int titles);

  int drawbeam(float xpos, float ypos, float bmaj, float bmin, float bpa);

  printf("Creating a spectral index error map...\n");

  // Declaring variables

  int i, j, istart, jstart, iend, jend, isize, jsize;
  float zoomedx, zoomedy, adj_maxrange, minrange, maxrange, *tmp_smootharray;
  char datebuff[32], historytext[MAXCMDLENGTH], flttochar[MAXCMDLENGTH];
  int xaxismax = xdim[currentset];
  int yaxismax = ydim[currentset];

  int intenscale = 0;

  // Declaring PGPLOT prototypes
  int cpgopen(const char *device);
  void cpgclos(void);
  void cpgenv(float xmin, float xmax, float ymin, float ymax, int just, int axis);
  void cpgpt1(float xpt, float ypt, int symbol);
  void cpgsci(int ci);
  void cpgimag(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float a1, float a2, const float *tr);
  void cpgctab(const float *l, const float *r, const float *g, const float *b, int nc, float contra, float bright);
  void cpgwedg(const char *side, float disp, float width, float fg, float bg, const char *label);
  void cpgscr(int ci, float cr, float cg, float cb);
  void cpgcont(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, const float *c, int nc, const float *tr);
  void cpgsls(int ls);
  void cpgsch(float size);

  // Declare funtools variables for FITS export
  Fun *funconn;

  // Set up variables and structs for times and dates
  time_t currenttime;
  struct tm * time_struct;


  if ( (exportfits !=1) || (export != 1) ) { // Do the standard plotting

    if (extwincontrol != 1) { 
      cpgopen(output);
    }

    // Swap to background and forground colours to something sensible
    cpgscr(0, 1, 1, 1);
    cpgscr(1, 0, 0, 0);

    cpgsci(1);

    // Calculate the zoomed start and end points
    zoomedx =  xdim[currentset] / zoom;
    zoomedy = ydim[currentset] / zoom;

    istart = xshift + border + ( (xdim[currentset] - (int)zoomedx) / 2);
    jstart = yshift + border + ( (ydim[currentset]  - (int)zoomedy) / 2);

    iend = xshift + xdim[currentset] - ( ( (xdim[currentset] - (int)zoomedx) / 2) + border);
    jend = yshift + ydim[currentset] - ( ( (ydim[currentset] - (int)zoomedy) / 2) + border);


    if (titles == 0) { sprintf(maptitle," ");  }

    mapscale(scaletype, istart, iend, jstart, jend, border, cellsize, axismarks, labels, ra, dec, maptitle, xshift, yshift, largetxt, export, titles);

 
    if (smoothmaps == 1) {

      minrange = 1e23;
      maxrange = -1e23;

      isize = (iend-istart) + (2*border);
      jsize = (jend-jstart) + (2*border);

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

	    tmp_smootharray[i*jsize+j] = specindexerror[currentset][regionarray[currentset][i+(jstart-border)][j+(istart-border)]];

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

      float transshift[6]={(istart-border), 1.0, 0.0, (jstart-border), 0.0, 1.0};

      cpgctab(rl,rr,rg,rb,11,contra,bright);
      if (largetxt == 1) { cpgsch(1.30); }
      cpgwedg("RI", 1.5, 5.0, minrange, adj_maxrange, "Spectral Index Error");
      if (largetxt == 1) { cpgsch(1.00); } // Revert to default

      cpgimag(tmp_smootharray,isize,jsize,border,isize-border,border,jsize-border,minrange,adj_maxrange,transshift);

      // printf("isize: %d jsize: %d border: %d minrange: %.3f adj_maxrange: %.3f\n", isize,jsize,border, minrange,adj_maxrange);
      
      free(tmp_smootharray);
    }

    else {

      // Use the original mapping style

      // Loop through setting and plotting each region

      for (i=border; i<(xdim[currentset]-border); i++){
	for (j=border; j<(ydim[currentset]-border); j++){

	  if (regionarray[currentset][i][j] == -1) {
	    intenscale = 0;   
	  }

	  else if (regionarray[currentset][i][j] == -2) {
	    intenscale = 0;
	  }

	  else if (regionarray[currentset][i][j] == 0) {
	    intenscale = 0;
	  }

	  else if (specindexerror[currentset][regionarray[currentset][i][j]] > 1e33) {
	    intenscale = 0;
	  }

	  else {

	    if (specindexerror[currentset][regionarray[currentset][i][j]] >= (maxsi[currentset] - ((maxsi[currentset]-minsi[currentset])/15)) ) {
	      intenscale = 1;
	    }
	    else if ( (specindexerror[currentset][regionarray[currentset][i][j]] < (maxsi[currentset] - ((maxsi[currentset]-minsi[currentset])/15))) && (specindexerror[currentset][regionarray[currentset][i][j]] >= (maxsi[currentset] - (2*((maxsi[currentset]-minsi[currentset])/15)))) ) {
	      intenscale = 2;
	    }
	    else if ( (specindexerror[currentset][regionarray[currentset][i][j]] < (maxsi[currentset] - (2*((maxsi[currentset]-minsi[currentset])/15)))) && (specindexerror[currentset][regionarray[currentset][i][j]] >= (maxsi[currentset] - (3*((maxsi[currentset]-minsi[currentset])/15)))) ) {
	      intenscale = 3;
	    }
	    else if ( (specindexerror[currentset][regionarray[currentset][i][j]] < (maxsi[currentset] - (3*((maxsi[currentset]-minsi[currentset])/15)))) && (specindexerror[currentset][regionarray[currentset][i][j]] >= (maxsi[currentset] - (4*((maxsi[currentset]-minsi[currentset])/15)))) ) {
	      intenscale = 4;
	    }
	    else if ( (specindexerror[currentset][regionarray[currentset][i][j]] < (maxsi[currentset] - (4*((maxsi[currentset]-minsi[currentset])/15)))) && (specindexerror[currentset][regionarray[currentset][i][j]] >= (maxsi[currentset] - (5*((maxsi[currentset]-minsi[currentset])/15)))) ) {
	      intenscale = 5;
	    }
	    else if ( (specindexerror[currentset][regionarray[currentset][i][j]] < (maxsi[currentset] - (5*((maxsi[currentset]-minsi[currentset])/15)))) && (specindexerror[currentset][regionarray[currentset][i][j]] >= (maxsi[currentset] - (6*((maxsi[currentset]-minsi[currentset])/15)))) ) {
	      intenscale = 6;
	    }
	    else if ( (specindexerror[currentset][regionarray[currentset][i][j]] < (maxsi[currentset] - (6*((maxsi[currentset]-minsi[currentset])/15)))) && (specindexerror[currentset][regionarray[currentset][i][j]] >= (maxsi[currentset] - (7*((maxsi[currentset]-minsi[currentset])/15)))) ) {
	      intenscale = 7;
	    }
	    else if ( (specindexerror[currentset][regionarray[currentset][i][j]] < (maxsi[currentset] - (7*((maxsi[currentset]-minsi[currentset])/15)))) && (specindexerror[currentset][regionarray[currentset][i][j]] >= (maxsi[currentset] - (8*((maxsi[currentset]-minsi[currentset])/15)))) ) {
	      intenscale = 8;
	    }
	  
	    else if ( (specindexerror[currentset][regionarray[currentset][i][j]] < (maxsi[currentset] - (8*((maxsi[currentset]-minsi[currentset])/15)))) && (specindexerror[currentset][regionarray[currentset][i][j]] >= (maxsi[currentset] - (9*((maxsi[currentset]-minsi[currentset])/15)))) ) {
	      intenscale = 9;
	    }
	    else if ( (specindexerror[currentset][regionarray[currentset][i][j]] < (maxsi[currentset] - (9*((maxsi[currentset]-minsi[currentset])/15)))) && (specindexerror[currentset][regionarray[currentset][i][j]] >= (maxsi[currentset] - (10*((maxsi[currentset]-minsi[currentset])/15)))) ) {
	      intenscale = 10;
	    }
	    else if ( (specindexerror[currentset][regionarray[currentset][i][j]] < (maxsi[currentset] - (10*((maxsi[currentset]-minsi[currentset])/15)))) && (specindexerror[currentset][regionarray[currentset][i][j]] >= (maxsi[currentset] - (11*((maxsi[currentset]-minsi[currentset])/15)))) ) {
	      intenscale = 11;
	    }

	    else if ( (specindexerror[currentset][regionarray[currentset][i][j]] < (maxsi[currentset] - (11*((maxsi[currentset]-minsi[currentset])/15)))) && (specindexerror[currentset][regionarray[currentset][i][j]] >= (maxsi[currentset] - (12*((maxsi[currentset]-minsi[currentset])/15)))) ) {
	      intenscale = 12;
	    }
	  
	    else if ( (specindexerror[currentset][regionarray[currentset][i][j]] < (maxsi[currentset] - (12*((maxsi[currentset]-minsi[currentset])/15)))) && (specindexerror[currentset][regionarray[currentset][i][j]] >= (maxsi[currentset] - (13*((maxsi[currentset]-minsi[currentset])/15)))) ) {
	      intenscale = 13;
	    }
	    else if ( (specindexerror[currentset][regionarray[currentset][i][j]] < (maxsi[currentset] - (13*((maxsi[currentset]-minsi[currentset])/15)))) && (specindexerror[currentset][regionarray[currentset][i][j]] >= (maxsi[currentset] - (14*((maxsi[currentset]-minsi[currentset])/15)))) ) {
	      intenscale = 14;
	    }	 


	    else if ( (specindexerror[currentset][regionarray[currentset][i][j]] < (maxsi[currentset] - (14*((maxsi[currentset]-minsi[currentset])/15)))) && (specindexerror[currentset][regionarray[currentset][i][j]] >= minsi[currentset]) ) {
	      intenscale = 15;
	    }
	  
	    else {
	      intenscale = 0;
	    }
	  }
	  cpgsci(intenscale);
	  cpgpt1(j,i,symbol);
	}
      }

    }


    if (contours == 1) {

      minrange = 1e32;
      maxrange = -1e32;
      float currentcont, *tmp_contourarray;

      cpgsci(contourcolour);
      cpgsls(contourlinestyle);

      tmp_contourarray = (float *)calloc(isize*jsize, sizeof(float));

      for (i=0; i<isize; i++) {
	for (j=0; j<jsize; j++) {
	  if ( ((i+jstart-border) < 0) || ((j+istart-border) < 0) || ((i+jstart-border) >= xaxismax) || ((j+istart-border) >= yaxismax) ) {
	    tmp_contourarray[i*jsize+j] = 0.0;
	  }
	  else if (regionarray[currentset][i+(jstart-border)][j+(istart-border)] < 0) { // Set void regions to a high value to create a white background
	    tmp_contourarray[i*jsize+j] = 0.0;
	  }
	  else {

	    tmp_contourarray[i*jsize+j] = flux[i+(jstart-border)][j+(istart-border)];

	  
	    if (tmp_contourarray[i*jsize+j] < minrange) {
	      minrange = tmp_contourarray[i*jsize+j];
	    }
	  
	    if (tmp_contourarray[i*jsize+j] > maxrange) {
	      maxrange = tmp_contourarray[i*jsize+j];
	    }
	  }
	}
      }

      for (i=firstcontour; i<=contourlevels; i++) {

	if (contourlogs == 1) {
	  currentcont = pow(10, log10(minrange) + (i*(log10(maxrange)-log10(minrange))/contourlevels) );
	}
	else {
	  currentcont = minrange + i*(maxrange-minrange)/contourlevels;
	}

	float transshift[6]={(istart-border), 1.0, 0.0, (jstart-border), 0.0, 1.0};
	cpgcont(tmp_contourarray, isize, jsize, border, isize-border, border, jsize-border, &currentcont, -1, transshift);
      }

      free(tmp_contourarray);
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


    if (extwincontrol != 1) {
      cpgclos();
    }
  }

  else if ( (exportfits == 1) && (export == 1) ) { // If we are exporting as fits
      maxrange = -1e23;
      minrange = 1e23;

      tmp_smootharray = (float *)calloc(xaxismax*yaxismax, sizeof(float));

      for (i=0; i<xaxismax; i++) {
	for (j=0; j<yaxismax; j++) {
	
	  if (regionarray[currentset][i][j] < 0) { // Set void regions to a high value to create a white background
	    tmp_smootharray[i*xaxismax+j] = NAN;
	  }
	  else {
	    tmp_smootharray[i*xaxismax+j] = specindexerror[currentset][regionarray[currentset][i][j]];

	    if (tmp_smootharray[i*xaxismax+j] < minrange) {
	      minrange = tmp_smootharray[i*xaxismax+j];
	    }

	    if (tmp_smootharray[i*xaxismax+j] > maxrange) {
	      maxrange = tmp_smootharray[i*xaxismax+j];
	    }
	  }
	}
      }

      // Create the file and write the header
      funconn = FunOpen(output, "w", NULL);

      // Check we can access the file
      if (funconn == NULL){
	fprintf(stderr,"\n*** Error: Unable to open %s, please check the file exists and access permissions are correct. ***\n\n", output);
	return 404;
      }
 
      // Get the current date and time and format it for the header
      time( &currenttime );
      time_struct = localtime ( &currenttime );
      strftime(datebuff,32,"%Y-%m-%d", time_struct);


      FunParamPutb(funconn, "SIMPLE", 0, 1, NULL, 1);
      FunParamPuti(funconn, "BITPIX", 0, -32, NULL, 1);
      FunParamPuti(funconn, "NAXIS", 0, 2, "", 1);
      FunParamPuti(funconn, "NAXIS1", 0, xaxismax, "X dimension (px)", 1);
      FunParamPuti(funconn, "NAXIS2", 0, yaxismax, "Y dimension (px)", 1);
      FunParamPutb(funconn, "EXTEND", 0, 1, "Tables following main image", 1);
      FunParamPutb(funconn, "BLOCKED", 0, 1, "Tape may be blocked", 1);
      FunParamPuts(funconn, "OBJECT", 0,  settarget, "Source name", 1);
      FunParamPuts(funconn, "TELESCOPE", 0, "BRATS", "Broadband Radio Astronomy Tools", 1);
      FunParamPuts(funconn, "INSTRUME", 0, "", NULL, 1);
      FunParamPuts(funconn, "OBSERVER", 0, "BRATS", "Broadband Radio Astronomy Tools", 1);
      //FunParamPuts(funconn,"DATE-OBS", 0, "", "Obs start date YYYY-MM-DD", 1);
      FunParamPuts(funconn, "DATE-MAP", 0, datebuff, "Last processing date YYYY-MM-DD", 1);
      FunParamPutd(funconn, "BSCALE", 0, 1.0, 11, "REAL = TAPE * BSCALE + BZERO", 1);
      FunParamPutd(funconn, "BZERO", 0, 0.0, 9, NULL, 1);
      FunParamPuts(funconn, "BUNIT", 0,  "", "Units of spectral index", 1);
      FunParamPutd(funconn, "EQUINOX", 0, eq, 11, "Epoch of RA DEC", 1);
      FunParamPutd(funconn, "OBSRA", 0, ra, 11, "Antenna pointing RA", 1);
      FunParamPutd(funconn, "OBSDEC", 0, dec, 11, "Antenna pointing DEC", 1);
      FunParamPutd(funconn, "DATAMAX", 0, maxrange, 11, "Maximum pixel value", 1);
      FunParamPutd(funconn, "DATAMIN", 0, minrange, 11, "Minimum pixel value", 1);
      FunParamPuts(funconn, "CTYPE1", 0, "RA---SIN", NULL, 1);
      FunParamPutd(funconn, "CRVAL1", 0, ra, 11, "RA", 1);
      FunParamPutd(funconn, "CDELT1", 0, delt1, 11, NULL, 1);
      FunParamPutd(funconn, "CRPIX1", 0, crpix1, 11, NULL, 1);
      FunParamPutd(funconn, "CROTA1", 0, crota1, 11, NULL, 1);
      FunParamPuts(funconn, "CTYPE2", 0, "DEC--SIN", NULL, 1);
      FunParamPutd(funconn, "CRVAL2", 0, dec, 11, "Dec", 1);
      FunParamPutd(funconn, "CDELT2", 0, delt2, 11, NULL, 1);
      FunParamPutd(funconn, "CRPIX2", 0, crpix2, 11, NULL, 1);
      FunParamPutd(funconn, "CROTA2", 0, crota2, 11, NULL, 1);
      FunParamPuts(funconn, "ORIGIN", 0, "BRATS", "Broadband Radio Astronomy Tools", 1);


      strcpy(historytext, "--------------------------------------------------------------------");
      FunParamPuts(funconn, "HISTORY", 0, historytext, NULL, 1);
      strcpy(historytext, "BRATS  Created by the Broadband Radio Astronomy Tools");
      FunParamPuts(funconn, "HISTORY", 0, historytext, NULL, 1);
      strcpy(historytext, "BRATS  http://www.askanastronomer.co.uk/brats");
      FunParamPuts(funconn, "HISTORY", 0, historytext, NULL, 1);
      strcpy(historytext, "BRATS  Spectral index error map between ");
      sprintf(flttochar, "%.4e", minfreq);
      strcat(historytext, flttochar);
      strcat(historytext, " and ");
      sprintf(flttochar, "%.4e", maxfreq);
      strcat(historytext, flttochar);
      strcat(historytext, " Hz.");

      FunParamPuts(funconn, "HISTORY", 0, historytext, NULL, 1);


      // Output the image
      FunImagePut(funconn, tmp_smootharray, xaxismax, yaxismax, -32, NULL);
      FunClose(funconn);

      free(tmp_smootharray);

  }
  else {
    fprintf(stderr,"\nError: Unknown display or export type, please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");
    return 404;
  }

  return 0;

}
