 
// Outputs a simple flux map of anything above the RMS noise using the raw data

//#include "mapscale.h"
//#include "drawbeam.h"

int fluxmap(float ****flux, float **rms, int *xdim, int *ydim, int mapnumber, float **maxflux, float **minflux, int currentset, int border, float zoom, float sigmalevel, char fluxmaptitle[], int titles, int labels, int axismarks, int extwincontrol, int smoothmaps, int fluxlogs, int fluxcut, float usercut, float wedgemultiplier, char *output, float ra, float dec, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, int xshift, int yshift, int largetxt, int export) { 

#include "smoothmaps.h"

  int mapscale(int scaletype, int istart, int iend, int jstart, int jend, int border, float cellsize, int axismarks, int labels, float ra, float dec, char *maptitle, int xshift, int yshift, int largetxt, int export, int titles);

  int drawbeam(float xpos, float ypos, float bmaj, float bmin, float bpa);

  int cpgopen(const char *device);
  void cpgclos(void);
  void cpgpt1(float xpt, float ypt, int symbol);
  void cpgsci(int ci);
  void cpgask(int flag);
  void cpgimag(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float a1, float a2, const float *tr);
  void cpgctab(const float *l, const float *r, const float *g, const float *b, int nc, float contra, float bright);
  void cpgwedg(const char *side, float disp, float width, float fg, float bg, const char *label);
  void cpgscr(int ci, float cr, float cg, float cb);
  void cpgsch(float size);

  int i, j, isize, jsize;
  int intenscale = 0;
  float scaleref = (maxflux[currentset][mapnumber] - minflux[currentset][mapnumber]);
  float zoomedx, zoomedy;
  int istart, jstart, iend, jend;

 
  if (extwincontrol == 1) {
    cpgask(1);
  }
  else {
    cpgopen(output);
  }

  // Swap to background and forground colours to something sensible
  cpgscr(0, 1, 1, 1);
  cpgscr(1, 0, 0, 0);

  cpgsci(1);

  // Calculate the zoomed start and end points
  zoomedx = xdim[currentset] / zoom;
  zoomedy = ydim[currentset] / zoom;

  istart = xshift + border + ( (xdim[currentset] - (int)zoomedx) / 2);
  jstart = yshift + border + ( (ydim[currentset]  - (int)zoomedy) / 2);

  iend = xshift + xdim[currentset] - ( ( (xdim[currentset] - (int)zoomedx) / 2) + border);
  jend = yshift + ydim[currentset] - ( ( (ydim[currentset] - (int)zoomedy) / 2) + border);

  isize = (iend-istart) + (2*border);
  jsize = (jend-jstart) + (2*border);

  if (titles == 0) { sprintf(fluxmaptitle," ");  }

  mapscale(scaletype, istart, iend, jstart, jend, border, cellsize, axismarks, labels, ra, dec, fluxmaptitle, xshift, yshift, largetxt, export, titles);

  
  if (smoothmaps == 1) {

    float *tmp_smootharray;
    float minrange = 1e23;
    float maxrange = -1e23;
    int xaxismax = xdim[currentset];
    int yaxismax = ydim[currentset];
  
    tmp_smootharray = (float *)calloc(isize*jsize, sizeof(float));

    for (i=0; i<isize; i++) {
      for (j=0; j<jsize; j++) {
	if ( ((i+jstart-border) < 0) || ((j+istart-border) < 0) || ((i+jstart-border) >= xaxismax) || ((j+istart-border) >= yaxismax) ) {
	  tmp_smootharray[i*jsize+j] = 1e24;
	}
	else if ( (flux[currentset][mapnumber][i+(jstart-border)][j+(istart-border)] < (sigmalevel * rms[currentset][mapnumber])) && (fluxcut != 1) ) { // Set void regions to a high value to create a white background
	  tmp_smootharray[i*jsize+j] = 1e24;
	}
	else if ( (fluxcut == 1) && (flux[currentset][mapnumber][i+(jstart-border)][j+(istart-border)] < usercut) ) {
	  tmp_smootharray[i*jsize+j] = 1e24;
	}
	else {

	  tmp_smootharray[i*jsize+j] = flux[currentset][mapnumber][i+(jstart-border)][j+(istart-border)];

	  if (tmp_smootharray[i*jsize+j] < minrange) {
	    minrange = tmp_smootharray[i*jsize+j];
	  }

	  if (tmp_smootharray[i*jsize+j] > maxrange) {
	    maxrange = tmp_smootharray[i*jsize+j];
	  }

	  if (fluxlogs == 1) {
	    tmp_smootharray[i*jsize+j] = log10(tmp_smootharray[i*jsize+j]);
	  }
	}
      }
    }

    if (fluxlogs == 1) {
      maxrange *= wedgemultiplier;
      minrange = log10(minrange);
      maxrange = log10(maxrange);
      //maxrange = maxrange+((maxrange-minrange)*wedgemultiplier);
    }
    else {
      maxrange *= wedgemultiplier;
    }

    if (fluxcut == 1) {
      if (fluxlogs == 1) {
	minrange = log10(usercut);
      }
      else {
	minrange = usercut;
	printf("\n\nminrange: %.2e\n\n", minrange);
      }
    }

    cpgctab(rl,rr,rg,rb,11,contra,bright);

    if (largetxt == 1) { cpgsch(1.30); }

    if (fluxlogs == 1) {
      cpgwedg("RI", 1.5, 5.0, minrange, maxrange, "log(Flux (Jy))");
    }
    else {
      cpgwedg("RI", 1.5, 5.0, minrange, maxrange, "Flux (Jy)");
    }

    if (largetxt == 1) { cpgsch(1.00); } // Revert to default

    float transshift[6]={(istart-border), 1.0, 0.0, (jstart-border), 0.0, 1.0};
    cpgimag(tmp_smootharray,isize,jsize,border,isize-border,border,jsize-border,minrange,maxrange,transshift);

    free(tmp_smootharray);
  }

  else {

    for (i=istart; i<iend; i++){

      for (j=jstart; j<jend; j++){

	if (flux[currentset][mapnumber][i][j] > (sigmalevel * rms[currentset][mapnumber]) ) {

	  intenscale = ((((flux[currentset][mapnumber][i][j] - minflux[currentset][mapnumber]) / scaleref ) * 15) + 1);

	  if (intenscale > 15) {
	    intenscale = 15;
	  }

	  cpgsci(intenscale);
	  cpgpt1(j, i, 1);
	}
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

  if (extwincontrol != 1) {
    cpgclos();
  }
  
  return 0;

}



// produce a flux map based on the adaptive regions

int regionfluxmap(float ***regflux, int *xdim, int *ydim, int mapnumber, float **maxflux, float **minflux, int currentset, int border, float zoom, char fluxmaptitle[], int titles, int labels, int axismarks, int ***regionarray, int extwincontrol, int **regionsize, int *ext_regnumber, int doregaveraging, int smoothmaps, char *output, float ra, float dec, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, int xshift, int yshift, int largetxt, int export) { 

#include "smoothmaps.h"

  int mapscale(int scaletype, int istart, int iend, int jstart, int jend, int border, float cellsize, int axismarks, int labels, float ra, float dec, char *maptitle, int xshift, int yshift, int largetxt, int export, int titles);

  int drawbeam(float xpos, float ypos, float bmaj, float bmin, float bpa);

  int cpgopen(const char *device);
  void cpgclos(void);
  void cpgpt1(float xpt, float ypt, int symbol);
  void cpgsci(int ci);
  void cpgask(int flag);
  void cpgimag(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float a1, float a2, const float *tr);
  void cpgctab(const float *l, const float *r, const float *g, const float *b, int nc, float contra, float bright);
  void cpgwedg(const char *side, float disp, float width, float fg, float bg, const char *label);
  void cpgscr(int ci, float cr, float cg, float cb);
  void cpgsch(float size);

  int i, j, isize, jsize;
  int regnumber = 0;
  int intenscale = 0;
  float zoomedx, zoomedy, tmp_maxflux, tmp_minflux, testflux, scaleref;
  int istart, jstart, iend, jend;

  // Swap to background and forground colours to something sensible
  cpgscr(0, 1, 1, 1);
  cpgscr(1, 0, 0, 0);
  cpgsci(1);

  // If averaging is set, find the new max and min values

  if (doregaveraging == 1) {

    tmp_maxflux = -1e32;
    tmp_minflux = 1e32;
 
    for (i=0; i<ext_regnumber[currentset]; i++) {

      testflux = regflux[currentset][mapnumber][i] / regionsize[currentset][i];
	
      if (testflux > tmp_maxflux) {
	tmp_maxflux = testflux;
      }

      if (testflux < tmp_minflux) {
	tmp_minflux = testflux;
      }
    }

    scaleref = (tmp_maxflux - tmp_minflux);

  }

  else { // Set a normal scale ref
    scaleref = (maxflux[currentset][mapnumber] - minflux[currentset][mapnumber]);
  }

  
 // Calculate the zoomed start and end points
  zoomedx = xdim[currentset] / zoom;
  zoomedy = ydim[currentset] / zoom;
  
  istart = xshift + border + ( (xdim[currentset] - (int)zoomedx) / 2);
  jstart = yshift + border + ( (ydim[currentset]  - (int)zoomedy) / 2);

  iend = xshift + xdim[currentset] - ( ( (xdim[currentset] - (int)zoomedx) / 2) + border);
  jend = yshift + ydim[currentset] - ( ( (ydim[currentset] - (int)zoomedy) / 2) + border);

  isize = (iend-istart) + (2*border);
  jsize = (jend-jstart) + (2*border);

  if (titles == 0) { sprintf(fluxmaptitle," ");  }

  mapscale(scaletype, istart, iend, jstart, jend, border, cellsize, axismarks, labels, ra, dec, fluxmaptitle, xshift, yshift, largetxt, export, titles);


 if (smoothmaps == 1) {

    float *tmp_smootharray;
    float minrange = 1e23;
    float maxrange = -1e23;
    int xaxismax = xdim[currentset];
    int yaxismax = ydim[currentset];

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

	  regnumber = regionarray[currentset][i+(jstart-border)][j+(istart-border)];

	  tmp_smootharray[i*jsize+j] = regflux[currentset][mapnumber][regnumber];

	  if (tmp_smootharray[i*jsize+j] < minrange) {
	    minrange = tmp_smootharray[i*jsize+j];
	  }

	  if (tmp_smootharray[i*jsize+j] > maxrange) {
	    maxrange = tmp_smootharray[i*jsize+j];
	  }
	}
      }
    }

    cpgctab(rl,rr,rg,rb,11,contra,bright);
    if (largetxt == 1) { cpgsch(1.30); }
    cpgwedg("RI", 1.5, 5.0, minrange, maxrange, "Flux (Jy)");
    if (largetxt == 1) { cpgsch(1.00); } // Revert to default

    float transshift[6]={(istart-border), 1.0, 0.0, (jstart-border), 0.0, 1.0};
    cpgimag(tmp_smootharray,isize,jsize,border,isize-border,border,jsize-border,minrange,maxrange,transshift);

    free(tmp_smootharray);
  }

 else {

   for (i=istart; i<iend; i++){

     for (j=jstart; j<jend; j++){

       regnumber = regionarray[currentset][i][j];

       if (regnumber > 0) {

	 if (doregaveraging == 1) {
	   intenscale = (((( (regflux[currentset][mapnumber][regnumber] / regionsize[currentset][regnumber]) - tmp_minflux) / scaleref ) * 15) + 1);
	 }

	 else {
	   intenscale = ((((regflux[currentset][mapnumber][regnumber] - minflux[currentset][mapnumber]) / scaleref ) * 15) + 1);
	 }
       }
       else {
	 intenscale = 0;
       }

       if (intenscale > 15) {
	 intenscale = 15;
       }

       cpgsci(intenscale);
       cpgpt1(j, i, 1);
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


 if (extwincontrol != 1) {
   cpgclos();
 }

 return 0;

}
