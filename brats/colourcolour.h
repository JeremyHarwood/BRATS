/* 
   Plots the total flux of the source as a function of frequency
*/


int colourcolour(float *colour1, float *colour2, float mincolour1, float maxcolour1, float mincolour2, float maxcolour2, float map1freq, float map2freq, float map3freq, char *settarget, int regnumber, int axismarks, int titles, int labels, int symbol, int extwincontrol, char *output, int swapsign) {


  // Declaring variables

  int i;

  float minx, miny, maxx, maxy;

  // const float c = 2.99792458e8; //Only needed when axis are in terms of wavelength rather than frequency

  char maptitle[64];
  char xlabel[64];
  char ylabel[64];

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
  void cpgmove(float x, float y);
  void cpgdraw(float x, float y);
  

  // Set the axis ranges

  if (swapsign == 1) {

    maxx = 0;
    maxy = 0;

    if (maxcolour1 > maxcolour2) {
      minx = miny = -(maxcolour1 * 1.2);
    }
    else {
      minx = miny = -(maxcolour2 * 1.2);
    }

  }

  else {

    minx = 0;
    miny = 0;

    if (maxcolour1 > maxcolour2) {
      maxx = maxy = maxcolour1 * 1.2;
    }
    else {
      maxx = maxy = maxcolour2 * 1.2;
    }

  }


  if (extwincontrol != 1) {
    cpgopen(output);
  }

  // Swap to background and forground colours to something sensible
  cpgscr(0, 1, 1, 1);
  cpgscr(1, 0, 0, 0);
  
  cpgsci(1);

  if (axismarks == 1) {
    cpgenv(minx, maxx, miny, maxy, 0, 0);
    //cpgenv(0, 2.5, 0, 2.5, 0, 0);
  }
  else {
    cpgenv(minx, maxx, miny, maxy, 0, -1);
  }


  if (titles == 0) {
    sprintf(maptitle, " ");
  }
  else {
    sprintf(maptitle, "Colour-Colour Plot for %s", settarget);
  }

  if (labels == 1) {
    //sprintf(xlabel,"\\ga\\d%.0f\\u\\u%.0f", (c/map1freq)*100, (c/map2freq)*100); // In cm
    //sprintf(ylabel,"\\ga\\d%.0f\\u\\u%.0f", (c/map2freq)*100, (c/map3freq)*100);
    sprintf(xlabel,"\\ga\\d%.0f\\u\\u%.0f", (map1freq/1e6), (map2freq/1e6)); // in MHz
    sprintf(ylabel,"\\ga\\d%.0f\\u\\u%.0f", (map2freq/1e6), (map3freq/1e6));
  }
  else {
    sprintf(xlabel," ");
    sprintf(ylabel," ");
  }

  cpglab(xlabel, ylabel, maptitle);

  if (swapsign == 1) {

    for (i=0; i < regnumber; i++) {

      colour1[i] = -colour1[i];
      colour2[i] = -colour2[i];
    }
  }

  cpgpt(regnumber, colour1, colour2, symbol);

  // Draw index a = index b

  cpgsci(2);
  cpgmove(minx, miny);
  cpgdraw(maxx, maxy);
  cpgsci(1);

if (extwincontrol != 1) {
    cpgclos();
  }


return 0;

}
