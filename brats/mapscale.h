/* 
  Set the map scale and its labels.
*/

#include <stdio.h>
#include <stdlib.h>

int mapscale(int scaletype, int istart, int iend, int jstart, int jend, int border, float cellsize, int axismarks, int labels, float ra, float dec, char *maptitle, int xshift, int yshift, int largetxt, int export, int titles) {

  void cpgenv(float xmin, float xmax, float ymin, float ymax, int just, int axis);
  void cpgaxis(const char *opt, float x1, float y1, float x2, float y2, float v1, float v2, float step, int nsub, float dmajl, float dmajr, float fmin, float disp, float orient);
  void cpgswin(float x1, float x2, float y1, float y2);
  void cpgtbox(const char *xopt, float xtick, int nxsub, const char *yopt, float ytick, int nysub);
  void cpglab(const char *xlbl, const char *ylbl, const char *toplbl);
  void cpgsch(float size);
  void cpgsvp(float xleft, float xright, float ybot, float ytop);
  void cpgpage(void);
  void cpgpap(float width, float aspect);
  void cpgbox(const char *xopt, float xtick, int nxsub, const char *yopt, float ytick, int nysub);
  void cpgdraw(float x, float y);
  void cpgmove(float x, float y);

  char xlabel[64];
  char ylabel[64];

  float conv_ra, conv_dec, xstartval, xendval, ystartval, yendval;

  float degtorad = 3.14159/180;
  float wholecirclesecs= 86400.0;

  //int largetxt = 1;

  // Enlarge the environment if large text mode is set
  if (largetxt == 1) {
    cpgsch(1.30);
  }
  else {
    cpgsch(1.00); // Check we are reset
  }


  if (scaletype == 1) { // Arcsec

    // Determing the axis values for the given cellsize

    xstartval = 0 - ( ( ( (iend+border) - (istart-border) ) / 2 ) * cellsize);
    xendval = 0 + ( ( ( (iend+border) - (istart-border) ) / 2 ) * cellsize);

    ystartval = 0 - ( ( ( (jend+border) - (jstart-border) ) / 2 ) * cellsize);
    yendval = 0 + ( ( ( (jend+border) - (jstart-border) ) / 2 ) * cellsize);

    if (export == 1) {
      if ( (largetxt == 1) && (labels == 1) ) {
    	if (titles == 1) {
    	  cpgpap(10,0.9);
    	} else {
    	  cpgpap(10,0.85);
    	}
      } else {
    	if (titles == 1) {
    	  cpgpap(10,0.95);
    	} else {
    	  cpgpap(10,0.88);
    	}
      }
    }

    cpgenv(istart-border, iend+border, jstart-border, jend+border, 1, -2);

    if (export == 1) {
      if ( (largetxt == 1) && (labels == 1) ) {
    	if (titles == 1) {
    	  cpgsvp(0.09, 0.8055, 0.115, 0.91);
    	} else {
    	  cpgsvp(0.09, 0.83375, 0.115, 0.99);
    	}
      } else {
    	if (titles == 1) {
    	  cpgsvp(0.08, 0.859, 0.09, 0.91);
    	} else {
    	  cpgsvp(0.08, 0.872, 0.09, 0.99);
    	}
      }
    }

    if (largetxt == 1) {
      cpgsch(1.30);
    }
    else {
      cpgsch(1.00); // Check we are reset
    }

    // Draw the X axis
    cpgaxis("N", istart-border, jstart-border, iend+border, jstart-border, xstartval, xendval, 0.0, 10, 0.7, 0.0, 0.4, 0.2, 0.0);

    // Now the Y axis
    cpgaxis("N", istart-border, jstart-border, istart-border, jend+border, ystartval, yendval, 0.0, 10, 0.0, 0.7, 0.4, -1.2, 0.0);
  
    // There is some strange depenence between the pgplot task box compared to the tbox used later. Just drawing the other axis to save my sanity...
    cpgmove(istart-border, jend+border);
    cpgdraw(iend+border, jend+border);

    cpgmove(iend+border, jend+border);
    cpgdraw(iend+border, jstart-border);
  }


  else if (scaletype == 2) { // DMS

    // Determing the axis values for the given cellsize
    xstartval = 0 - ( ( ( (iend+border) - (istart-border) ) / 2 ) * cellsize);
    xendval = 0 + ( ( ( (iend+border) - (istart-border) ) / 2 ) * cellsize);

    ystartval = 0 - ( ( ( (jend+border) - (jstart-border) ) / 2 ) * cellsize);
    yendval = 0 + ( ( ( (jend+border) - (jstart-border) ) / 2 ) * cellsize);


    // Setup the enviroment depending on if we are exporting
    if (export == 0) {
      cpgenv(istart-border, iend+border, jstart-border, jend+border, 1, -1);
    }
    else {
      cpgpage();
    }

    // Set to range in arcseconds
    cpgswin(xstartval, xendval, ystartval, yendval);

    if (export == 1) {
      if ( (largetxt == 1) && (labels == 1) ) {
	if (titles == 1) {
	  cpgpap(10,0.9);
	  cpgsvp(0.09, 0.8055, 0.115, 0.91);
	} else {
	  cpgpap(10,0.85);
	  cpgsvp(0.09, 0.83375, 0.115, 0.99);
	}
      } else {
	if (titles == 1) {
	  cpgpap(10,0.95);
	  cpgsvp(0.08, 0.859, 0.09, 0.91);
	} else {
	  cpgpap(10,0.88);
	  cpgsvp(0.08, 0.872, 0.09, 0.99);
	}
      }
    }

    // Draw the axis box
    if ( (axismarks == 1)){
      if (largetxt == 1) {
	cpgsch(1.30);
      }
      else {
	cpgsch(1.00); // Check we are reset
      }
      cpgtbox("BCZDNTSXYO", 0.0, 0.0, "BCZDNTSXYO", 0.0, 0.0);
    }
    else {
      cpgtbox("BCZDNXYO", 0.0, 0.0, "BCZDNXYO", 0.0, 0.0);
    }

  }

  else if (scaletype == 3) { // WCS

    float minmaxshift;

   // Convert the RA from decimal degrees to seconds. This will catch both positive values RA e.g. 240 deg and negative ones as well e.g. -120 deg as we are limiting to the 24 hours later.

   conv_dec = dec * 3600; // arcsec
   conv_dec += (yshift*cellsize);

   conv_ra = 240*(360 + ra); // arcsec

   conv_ra +=  ((wholecirclesecs*cos(degtorad * dec) ) / 360) * ((0-xshift) * (cellsize/3600)); // We change the sign of xshift here to account for the fact RA goes right to left
   

    // Determing the axis values for the given cellsize
    minmaxshift =  ((wholecirclesecs*cos(degtorad * dec) ) / 360) * ((((iend+border) - (istart-border) ) / 2 )  * (cellsize/3600) );

    xstartval = conv_ra + minmaxshift;
    xendval = conv_ra - minmaxshift;

    ystartval = conv_dec - ( ( ( (jend+border) - (jstart-border) ) / 2 ) * cellsize);
    yendval = conv_dec + ( ( ( (jend+border) - (jstart-border) ) / 2 ) * cellsize);

    // Setup the enviroment depending on if we are exporting
    if (export == 0) {
      cpgenv(istart-border, iend+border, jstart-border, jend+border, 1, -2);
    }
    else {
      cpgpage();
    }

    cpgswin(xstartval, xendval, ystartval, yendval);

    // Publication format. Maybe at a later date if I have time
    /* if (pubformat == 1) { */
    /*   if (wedgepos == 1) { */
    /* 	cpgpap(10,0.9); */
    /* 	cpgsvp(0.09, 0.8055, 0.115, 0.91); */
    /*   } */
    /*   else if (wedgepos == 2) { */
    /* 	cpgpap(10,0.85); */
    /* 	cpgsvp(0.09, 0.83375, 0.115, 0.99); */
    /*   } else { // If we get lost just stick it in the default right position */
    /* 	cpgpap(10,0.85); */
    /* 	cpgsvp(0.09, 0.83375, 0.115, 0.99); */
    /*   } */

    /* } else { */
      // Adapt the sizes depending on what options are set
      if (export == 1) {
	if ( (largetxt == 1) && (labels == 1) ) {
	  if (titles == 1) {
	    cpgpap(10,0.9);
	    cpgsvp(0.09, 0.8055, 0.115, 0.91);
	  } else {
	    cpgpap(10,0.85);
	    cpgsvp(0.09, 0.83375, 0.115, 0.99);
	  }
	} else {
	  if (titles == 1) {
	    cpgpap(10,0.95);
	    cpgsvp(0.08, 0.859, 0.09, 0.91);
	  } else {
	    cpgpap(10,0.88);
	    cpgsvp(0.08, 0.872, 0.09, 0.99);
	  }
	}
      }
      ///}
   /* else { */
   /*   if (largetxt == 1) { */
   /*     cpgsvp(0.08, 0.85, 0.09, 0.91); */
   /*   } else { */
   /*     printf("Here 4"); */
   /*     cpgsvp(0.08, 0.85, 0.09, 0.91); */
   /*   } */
   /* } */

    // Draw the axis box and set the scale in arcseconds

    if ( (axismarks == 1)){
      if (largetxt == 1) {
	cpgsch(1.30);
      }
      else {
	cpgsch(1.00); // Check we are reset
      }
      cpgtbox("XBCZHNTSYO", 0.0, 0.0, "BCZDNTSYO", 0.0, 0.0);
    }
    else {
      cpgtbox("BCZHNYO", 0.0, 0.0, "BCZDNYO", 0.0, 0.0);
    }

  }
  else { // Pixels

    if ( (axismarks == 1)){
      cpgenv(istart-border, iend+border, jstart-border, jend+border, 1, 0);
    }

    else {
      cpgenv(istart-border, iend+border, jstart-border, jend+border, 1, -1);
    }
  }


  if (labels == 1) {

    if (scaletype == 1) {
      sprintf(xlabel,"Angular Size (Arcseconds)");
      sprintf(ylabel,"Angular Size (Arcseconds)");
    }

    else if (scaletype == 2) {
      sprintf(xlabel,"Angular Size");
      sprintf(ylabel,"Angular Size");
    }

    else if (scaletype == 3) {
      sprintf(xlabel,"Right Ascension");
      sprintf(ylabel,"Declination");
    }

    else {
      sprintf(xlabel,"Pixels");
      sprintf(ylabel,"Pixels");
    }

  }
  else {
    sprintf(xlabel," ");
    sprintf(ylabel," ");
  }

  cpglab(xlabel, ylabel, maptitle);
 
 // If the enviroment is set to use seconds, change it back to pixels before returning
  if ( (scaletype == 2) || (scaletype == 3) ) {
    cpgswin(istart-border, iend+border, jstart-border, jend+border);
  }

  // Reset size to the the default
  cpgsch(1.00);


  //printf("cpgswin(istart-border: %d, iend+border: %d, jstart-border: %d, jend+border: %d)\n",istart-border, iend+border, jstart-border, jend+border);

  return 0;

}
