/* Simple module to plot spectral index against curvature

   What needs to be passed:

   float *curvature - Array of curvature values indexed in the same way as the spectral index array
   float *specindex - Array of spectral index values
   int regnumber - The total number of regions
   int symbol - Symbol to use for points (1 for dots). Se pgplot symbols for more options.
   float xmin, xmax, ymin, ymax - Minimum and maximum axis values

Basic usage:

  #include "plotcurvagainstspec.c"

  ...

  plotcurvagainstspec(cuvature, specindex, numberofregions, 1, 0, 10, 0 24);


   If you have any question please email jeremy.harwood@physics.org and I will try and answer them as soon as I can
*/

#include <gsl/gsl_fit.h>

int plotcurvagainstspec(float **curvature, float *specindex, int regnumber, int symbol, float xmin, float xmax, float ymin, float ymax, int poly, int autoscalex, int autoscaley, float userminx, float usermaxx, float userminy, float usermaxy, char *settarget, int titles, int labels, int axismarks, int extwincontrol) {

  // Declaring variables
  int i;
  char maptitle[128];

  // Delcare pgpplot prototypes
  int cpgopen(const char *device);
  void cpgclos(void);
  void cpgenv(float xmin, float xmax, float ymin, float ymax, int just, int axis);
  void cpglab(const char *xlbl, const char *ylbl, const char *toplbl);
  void cpgpt1(float xpt, float ypt, int symbol);
  void cpgsci(int ci);
  void cpgscr(int ci, float cr, float cg, float cb);

  // Plot the spectral index map
  
  if (autoscalex == 1) {
    if (xmin < 0) {
      xmin *= 1.2;
    }
    else {
      xmin *= 0.8;
    }

    if (xmax < 0) {
      xmax *= 0.8;
    }
    else {
      xmax *= 1.2;
    }
  }
  else {
    xmin = userminx;
    xmax = usermaxx;
  }


  if (autoscaley == 1) {
    if (ymin < 0) {
      ymin *= 1.2;
    }
    else {
      ymin *= 0.8;
    }

    if (ymax < 0) {
      ymax *= 0.8;
    }
    else {
      ymax *= 1.2;
    }
  }
  else {

    ymin = userminy;
    ymax = usermaxy;
  }


  if (extwincontrol != 1) {
    cpgopen("/xs");
  }


  // Swap to background and forground colours to something sensible
  cpgscr(0, 1, 1, 1);
  cpgscr(1, 0, 0, 0);
 
  // Layout information
  cpgsci(1);

  if (axismarks == 1) {
    cpgenv(xmin,xmax,ymin,ymax, 0, 0);
  }
  else {
    cpgenv(xmin,xmax,ymin,ymax, 0, -1);
  }


 // Set plot title
  if (titles == 0) {
    sprintf(maptitle, " ");
  }
  else {
    sprintf(maptitle,"Spectral Index vs Spectral Curvature for %s", settarget);
  }
  
  if (labels == 0) {
    cpglab("", "", maptitle);
  }
  else {
    cpglab("Spectral Index", "Spectral Curvature", maptitle);
  }


  // Plotting each region

  for(i=0; i<regnumber; i++) {

    cpgpt1(specindex[i], curvature[i][poly], symbol);

  }

if (extwincontrol != 1) {
  cpgclos();
}

  return 0;

}




int plotmodelvsflux(float *modelflux, float *regflux, float *fluxerrorplus, float *fluxerrorminus, int symbol, int nummodelpts, int imgnum, float *fluxfrequency, float *modelfrequency, int autoscalex, int autoscaley, float userminx, float usermaxx, float userminy, float usermaxy, float minfreq, float maxfreq, float minflux, float maxflux, char *settarget, int titles, int labels, int axismarks, double inject, double fieldstrength, float bestage, float bestoff,float chisquared, int paramlabels, int model, int usereduced, int ul, int *upperlimits, int largetxt) {


  // Declaring variables

  int i;
  double xaxismin, xaxismax, yaxismin, yaxismax;
  char str_inject[64], str_field[64], str_age[64], str_on[64], str_off[64], str_chi[64], ylabel[64], xlabel[64], maptitle[128];

  // Delcare pgpplot prototypes
  int cpgopen(const char *device);
  void cpgclos(void);
  void cpgenv(float xmin, float xmax, float ymin, float ymax, int just, int axis);
  void cpglab(const char *xlbl, const char *ylbl, const char *toplbl);
  void cpgpt(int n, const float *xpts, const float *ypts, int symbol);
  void cpgsci(int ci);
  void cpgerrb(int dir, int n, const float *x, const float *y, const float *e, float t);
  void cpgline(int n, const float *xpts, const float *ypts);
  void cpgerry(int n, const float *x, const float *y1, const float *y2, float t);
  void cpgscr(int ci, float cr, float cg, float cb);
  void cpgmtxt(const char *side, float disp, float coord, float fjust, const char *text);
  void cpgarro(float x1, float y1, float x2, float y2);
  void cpgsah(int fs, float angle, float barb);
  void cpgsch(float size);

  // Swap to background and forground colours to something sensible
  cpgscr(0, 1, 1, 1);
  cpgscr(1, 0, 0, 0);

  // Find the overall min and max values

  if (autoscalex == 1) {
    xaxismin = 0.9 * minfreq;
    xaxismax = 1.1 * maxfreq;
  }
  else {
    xaxismin = userminx;
    xaxismax = usermaxx;

  }

  if (autoscaley == 1) {
    yaxismin = 0.9 * minflux;
    yaxismax = 1.1 * maxflux;
  }
  else {
    yaxismin = userminy;
    yaxismax = usermaxy;
  }


  xaxismin = log10(xaxismin);
  xaxismax = log10(xaxismax);
  yaxismin = log10(yaxismin);
  yaxismax = log10(yaxismax);

  // Layout information

  cpgsci(1);

  if (titles == 0) { sprintf(maptitle," ");  }

  if (labels == 1) {
    sprintf(xlabel,"log(Frequency / Hz)");
    sprintf(ylabel,"log(Flux / Jy)");
  }
  else {
    sprintf(xlabel," ");
    sprintf(ylabel," ");
  }


  if (titles == 0) {
    sprintf(maptitle," "); 
  }
  else {
    sprintf(maptitle,"Model vs Observed Flux as a Function of Frequency");
  }

  if (axismarks == 1) {
    // Enlarge the environment if large text mode is set
    if (largetxt == 1) {
      cpgsch(1.30);
      cpgenv(xaxismin, xaxismax, yaxismin, yaxismax, 0, 0);
      cpgsch(1.00);
    }
    else {
      cpgenv(xaxismin, xaxismax, yaxismin, yaxismax, 0, 0);
    }
  }
  else {
    if (largetxt == 1) {
      cpgsch(1.30);
      cpgenv(xaxismin, xaxismax, yaxismin, yaxismax, 0, -1);
      cpgsch(1.00);
    }
    else {
      cpgenv(xaxismin, xaxismax, yaxismin, yaxismax, 0, -1);
    }
  }

  if (largetxt == 1) {
    cpgsch(1.40);
    cpglab(xlabel, ylabel, maptitle);
    cpgsch(1.00);
  }
  else {
    cpglab(xlabel, ylabel, maptitle);
  }

  // Plot each region
  cpgpt(imgnum, fluxfrequency, regflux, symbol);

  if (ul > 0) {
    cpgsah(2, 45.0, 1.0); // Make the arrow more suitable to the plot style

    float *tmp_fluxerrorplus = (float *)calloc(imgnum-ul,sizeof(float));
    float *tmp_fluxerrorminus  = (float *)calloc(imgnum-ul,sizeof(float));
    float *tmp_fluxerrorfrequency = (float *)calloc(imgnum-ul,sizeof(float));

    int j;
    int tmp_imgnum = 0;

    for (i=0; i<imgnum; i++) {
      for (j=0; j<ul; j++) {
	if ((i != upperlimits[j]) && (j == ul-1)){
	  tmp_fluxerrorfrequency[tmp_imgnum] = fluxfrequency[i];
	  tmp_fluxerrorplus[tmp_imgnum] = fluxerrorplus[i];
	  tmp_fluxerrorminus[tmp_imgnum] = fluxerrorminus[i];

	  tmp_imgnum++;

	  //printf("Freq: %.4e log(Freq): %.3f Flux: %.3f log(Flux): %.3f\n", pow(10,fluxfrequency[i]), fluxfrequency[i], pow(10, regflux[i]), regflux[i]);
	}
      }
    }

    cpgerry(imgnum-ul, tmp_fluxerrorfrequency, tmp_fluxerrorplus, tmp_fluxerrorminus, 0.5);
  
    cpgsch(0.65);

    for (i=0; i<ul; i++) {
      cpgarro(fluxfrequency[upperlimits[i]],regflux[upperlimits[i]],fluxfrequency[upperlimits[i]],regflux[upperlimits[i]]+(regflux[4]*0.08));
    }
      cpgsch(1.0);
  }
  else {
    cpgerry(imgnum, fluxfrequency, fluxerrorplus, fluxerrorminus, 0.5); // This array needs adapting to only include the non-upperlimits
  }

  cpgsci(2);
  cpgline(nummodelpts+1, modelfrequency, modelflux);

  cpgsci(1);
  if (paramlabels == 1) {

  if (largetxt == 1) {
    cpgsch(1.40);
  }

    if (model == 1) {
      cpgmtxt("B", -5.5, 0.05, 0.0, "JP Model");
    }
    if (model == 2) {
      cpgmtxt("B", -5.5, 0.05, 0.0, "KP Model");
    }
    if (model == 3) {
      cpgmtxt("B", -5.5, 0.05, 0.0, "Tribble Model (JP)");
    }
    else if (model == 4) {
      cpgmtxt("B", -5.5, 0.05, 0.0, "Tribble Model (KP)");
    }

    sprintf(str_inject, "Injection index: %.2f", inject);
    sprintf(str_field, "B-Field: %.2e T", fieldstrength);


    if (model == 5) {
      sprintf(str_age, "Spectral Age: %.2f Myr", bestage+bestoff);
      sprintf(str_on, "On time: %.2f", bestage);
      sprintf(str_off, "Off time: %.2f", bestoff);
    }
    else {
      sprintf(str_age, "Spectral Age: %.2f Myr", bestage);
    }


    if (usereduced == 1) {
      sprintf(str_chi, "Reduced Chi-Squared: %.2f", chisquared);
    }
    else {
      sprintf(str_chi, "Chi-Squared: %.2f", chisquared);
    }

    if (model == 5) {
      cpgmtxt("B", -6.5 , 0.05, 0.0, str_inject);
      cpgmtxt("B", -5.5 , 0.05, 0.0, str_field);
      cpgmtxt("B", -4.5 , 0.05, 0.0, str_age);
      cpgmtxt("B", -3.5 , 0.05, 0.0, str_on);
      cpgmtxt("B", -2.5 , 0.05, 0.0, str_off);
      cpgmtxt("B", -1.5 , 0.05, 0.0, str_chi);
    }
    else {
      cpgmtxt("B", -4.5 , 0.05, 0.0, str_inject);
      cpgmtxt("B", -3.5 , 0.05, 0.0, str_field);
      cpgmtxt("B", -2.5 , 0.05, 0.0, str_age);
      cpgmtxt("B", -1.5 , 0.05, 0.0, str_chi);
    }

    if (largetxt == 1) {
      cpgsch(1.00);
    }
  }

  return 0;

}



int plotfluxvsspecindex(float *regflux, float *fluxerrorplus, float *fluxerrorminus, int symbol, int imgnum, float *fluxfrequency, int autoscalex, int autoscaley, float userminx, float usermaxx, float userminy, float usermaxy, float minfreq, float maxfreq, float minflux, float maxflux, char *settarget, int titles, int labels, int axismarks, int paramlabels, double *specindex_modelflux, double *specindex_modelerror, int modelres, int *fluxorder, float specindex, int largetxt, float specindexerror, int specindex_errortype) {


  // Declaring variables
  int j;
  double xaxismin, xaxismax, yaxismin, yaxismax;
  float *tmp_modelfrequency, *tmp_modelflux, *tmp_posmodelerror, *tmp_negmodelerror;
  char str_alpha[64], ylabel[64], xlabel[64], maptitle[128];

  
  // Declare pgpplot prototypes
  int cpgopen(const char *device);
  void cpgclos(void);
  void cpgenv(float xmin, float xmax, float ymin, float ymax, int just, int axis);
  void cpglab(const char *xlbl, const char *ylbl, const char *toplbl);
  void cpgpt(int n, const float *xpts, const float *ypts, int symbol);
  void cpgsci(int ci);
  void cpgerrb(int dir, int n, const float *x, const float *y, const float *e, float t);
  void cpgline(int n, const float *xpts, const float *ypts);
  void cpgerry(int n, const float *x, const float *y1, const float *y2, float t);
  void cpgscr(int ci, float cr, float cg, float cb);
  void cpgmtxt(const char *side, float disp, float coord, float fjust, const char *text);
  void cpgsch(float size);
  void cpgpap(float width, float aspect);

  
  // Allocate memory for the spectral index data
  tmp_modelfrequency = (float *)calloc(modelres+1, sizeof(float));
  tmp_modelflux = (float *)calloc(modelres+1, sizeof(float));
  tmp_posmodelerror = (float *)calloc(modelres+1, sizeof(float));
  tmp_negmodelerror = (float *)calloc(modelres+1, sizeof(float));

  // Calculate the spectral index line to be plotted
  for (j=0; j<=modelres; j++) {

    tmp_modelfrequency[j] = fluxfrequency[0] + ( ( (fluxfrequency[imgnum-1] - fluxfrequency[0]) / modelres ) * j);

    // This is due to one function needing a double and the other a float (GSL vs PGPLOT)!
    tmp_modelflux[j] = specindex_modelflux[j];
    tmp_posmodelerror[j] = specindex_modelflux[j] + specindex_modelerror[j];
    tmp_negmodelerror[j] = specindex_modelflux[j] - specindex_modelerror[j];


    //printf("j: %d index: %.4f freq: %.4e flux: %.4e error %.4e\n", j, specindex, tmp_modelfrequency[j], tmp_modelflux[j], tmp_modelerror[j]);

  }

  // Swap to background and forground colours to something sensible
  cpgscr(0, 1, 1, 1);
  cpgscr(1, 0, 0, 0);


  // Find the overall min and max values

  if (autoscalex == 1) {
    xaxismin = 0.9 * minfreq;
    xaxismax = 1.1 * maxfreq;
  }
  else {
    xaxismin = userminx;
    xaxismax = usermaxx;
  }

  if (autoscaley == 1) {
    /*if ( (ul > 1) && (upperlimits == ) ) {
      yaxismin = 0.8 * minflux;
      }*/
    yaxismin = 0.9 * minflux;
    yaxismax = 1.1 * maxflux;
  }
  else {
    yaxismin = userminy;
    yaxismax = usermaxy;
  }


  xaxismin = log10(xaxismin);
  xaxismax = log10(xaxismax);
  yaxismin = log10(yaxismin);
  yaxismax = log10(yaxismax);

  // Layout information
  cpgsci(1);

  if (titles == 0) { sprintf(maptitle," ");  }

  if (largetxt == 1) { cpgsch(1.30); }

  if (labels == 1) {
    sprintf(xlabel,"log(Frequency / Hz)");
    sprintf(ylabel,"log(Flux / Jy)");
  }
  else {
    sprintf(xlabel," ");
    sprintf(ylabel," ");
  }


  if (titles == 0) {
    sprintf(maptitle," "); 
  }
  else {
    sprintf(maptitle,"Region Flux With Best Fitting Spectral Index Overlaid");
  }


  if (axismarks == 1) {
    cpgenv(xaxismin, xaxismax, yaxismin, yaxismax, 0, 0);
  }
  else {
    cpgenv(xaxismin, xaxismax, yaxismin, yaxismax, 0, -1);
  }


  cpglab(xlabel, ylabel, maptitle);


  // Plot each region

  cpgpt(imgnum, fluxfrequency, regflux, symbol);
  cpgerry(imgnum, fluxfrequency, fluxerrorplus, fluxerrorminus, 0.5);

  cpgsci(2);
  cpgline(modelres+1, tmp_modelfrequency, tmp_modelflux);

  if (specindex_errortype == 1) {
    cpgsci(3);
    cpgline(modelres+1, tmp_modelfrequency, tmp_posmodelerror);

    cpgsci(3);
    cpgline(modelres+1, tmp_modelfrequency, tmp_negmodelerror);
  }

  cpgsci(1);
  if (paramlabels == 1) {
    sprintf(str_alpha, "\\ga : %.2f +- %.2f", specindex, specindexerror);
    cpgmtxt("B", -2.2 , 0.05, 0.0, str_alpha);
    // Can add the 1 sigma error here later
  }

  if (largetxt == 1) { cpgsch(1.00); } // Revert to default

  return 0;

}



