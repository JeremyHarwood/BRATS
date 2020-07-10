/*
Known Issues

1. If the fits header is too long, the funtools library will cause the program to crash out with the error 'no WCS information in file while parsing filter at: XX:XX:XX.XXX'. This normally occurs due to a large history accumulated during reduction and can be fixed by running the AIPS Stalin command or history=False when exporting in CASA.

2. There is a known bug in some versions of GSL that causes the chi-squared CDF function to fail, particularly for a large number of DoF or very low chi-squared values. The suppressconf command has been included as a work around but this can also sometimes throw up oddities when mapping later on. In such cases exporting the data via exportdata and mapping in an external program is recommended or, ideally, update GSL to a version where this issue has been resolved.

3. The default GCC compiler that comes preinstalled with MacOS does not support multicore processing. This will cause an error during installation similar to: "FATAL:/somelocalpath/as/x86_64/as: I don't understand 'm' flag!". This can be resolved by installing the latest version of GCC.

4. The new fixed width when exporting png images sometimes cuts off long title names and/or the beginning of the xaxis scale label. This will be addressed at a later date but is currently worth sacraficing for the sake of easier to produce publication images. The xaxis scale issue usually only occurs on high zooms, and can be fixed by slightly adjusting the zoom value.

5. When outputting in FITS format, certain header parameters (e.g. CTYPE) are grouped together rather than set out in the traditional format. This appears to be a FUNTOOLS issue.

6. When using linear plots, the colour wedge scale is limited to 4 orders of magnitude. This appears to be a PGPLOT issue. The default is currently set to using a log scale which works around this problem for now.
*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <dirent.h>
#include <funtools.h>
#include <time.h>
#include <ctype.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <readline/readline.h>
#include <readline/history.h>

//Global command over the colour scheme
int pltcolour = 0;
int invertcolours = 0;

#include "chomp.h"
#include "help.h"
#include "load.h"
#include "syscom.h"
//#include "allocate.h"
#include "list.h"
#include "adaptiveregions.h"
#include "fluxmap.h"
#include "calcspecindex.h"
#include "plotspecindex.h"
#include "applyfixedregions.h"
#include "plotflux.h"
#include "curvature.h"
#include "plotpolyfit.h"
#include "plottotalflux.h"
#include "./synchrotron_src/synchrotron.c"
#include "comparisonplots.h"
#include "spectralagemap.h"
#include "chisquaredmap.h"
#include "colourcolour.h"
#include "singleregion.h"
#include "fluxerrors.h"
#include "mergeimages.h"
#include "mapfrommodel.h"
#include "diffmap.h"
#include "scalefluxes.h"
#include "errormap.h"
#include "fullexport.h"
#include "fullimport.h"
#include "resizeimages.h"
#include "./synchrotron_src/ci_model.h"
#include "injectbyregion.h"
#include "injectchibyregion.h"
#include "specindexerrormap.h"


// Definitions
#define VERSION "2.6.3.3"
#define MAXCMDLENGTH 1024
#define MAXDATASETS 20
#define MAXNUMMAPS 100
#define NUMBEROFMODELS 4
#define DEFAULTSIGMA 5.00
#define DEFAULTZOOM 1.0
#define DEFAULTBORDER 4
#define DEFAULTFRACBORDER 3.0 // This must be a decimal e.g. 3.0, not 3
#define DEFAULTITLES 1
#define DEFAULTLABELS 1
#define DEFAULTAXISMARKS 1
#define DEFAULTMAPTOMAP -1.0
#define DEFAULTHOTPIXELS 0.2
#define DEFAULTAVEREGIONS 0
#define DEFAULTSIGNOISE 1
#define DEFAULTSEARCHAREA 1
#define DEFAULTPRINTINDEX 0
#define MAXTITLELENGTH 128
#define DEFAULTUSELOGX 1
#define DEFAULTUSELOGY 1
#define DEFAULTAUTOSCALEX 1
#define DEFAULTAUTOSCALEY 1
#define DEFAULTPRINTINC 0
#define DEFAULTPOLY 2
#define DEFAULTPLOTRES 1000
#define DEFAULTREVCURVE 0
#define DEFAULYPRINTPOLY 0
#define DEFAULTERRFLAT 0.03
#define DEFAULTINSTRUMENT 1
#define DEFAULTBESTGUESS 1.0
#define DEFAULTSPECRANGE 0.5
#define DEFAULTNORMRANGE 0.2
#define DEFAULTSYMBOL 1
#define DEFAULTSPECRES 100
#define DEFAULTFLUXRES 100
#define DEFAULTMINMODELFREQ 1e+7
#define DEFAULTMAXMODELFREQ 1e+13
#define DEFAULTAGERES 10
#define DEFAULTMYEARS 100
#define DEFAULTNORMRES 10
#define DEFAULTONSOURCENOISE 3.0
#define DEFAULTMODELRES 100
#define DEFAULTCASADATA 1
#define DEFAULTLEVELS 3
#define DEFAULTPRINTRESULTS 0
#define DEFAULTSMOOTHMAPS 1
#define DEFAULTCONLEVELS 8
#define DEFAULTCONCOLOUR 1
#define DEFAULTCONLINESTYLE 1
#define DEFAULTCONTOURS 0
#define DEFAULTFIRSTCONTOUR 1
#define DEFAULTCONLOGS 1
#define DEFAULTIMAGELOC "./images"
#define DEFAULTDATALOC "./data"
#define DEFAULTIMAGETYPE "png"
#define DEFAULTDATAEXPORTEXT "txt"
#define DEFAULTEXPORT 0
#define DEFAULTFLUXLOGS 1
#define DEFAULTFLUXCUT 0
#define DEFAULTUSERCUT 1e-6
#define DEFAULTWEDGE 1.15
#define DEFAULTFIELDSTRENGTH 1e-9
#define DEFAULTINJECT 0.6
#define DEFAULTCOLOURSIGN 0
#define DEFAULTALWAYSZEROAGE 1
#define DEFAULTPARAMLABELS 1
#define DEFAULTPRINTREJECT 0
#define DEFAULTGMIN 10.0
#define DEFAULTGMAX 1000000.0
#define DEFAULTPOSMAP 0
#define DEFAULTSCALETYPE 3
#define DEFAULTBEAMPOSX 95
#define DEFAULTBEAMPOSY 5
#define DEFAULTPLOTBEAM 1
#define DEFAULTUSEREDUCED 1
#define DEFAULTMODELMYEARS 10
#define DEFAULTMININJECT 0.5
#define DEFAULTMAXINJECT 1.0
#define DEFAULTINJECTINTERVAL 10
#define DEFAULTXSHIFT 0 // These values should be in the code format i.e. opposite signs to the standard user input
#define DEFAULTYSHIFT 0
#define DEFAULTFLUXCALERR 3.0
#define DEFAULTEXAMPLEREDSHIFT 0.2
#define DEFAULTEXPORTCOMPRESSION 1
#define DEFAULTEXPECTEDHEADERS 7
#define DEFAULTEXPECTEDMODELHEADERS 3
#define DEFAULTINDEXCALCTYPE 2
#define DEFAULTMINMODELMYEARS 0
#define DEFAULTMODELMINOFF 0
#define DEFAULTMODELMAXOFF 20
#define DEFAULTMINOFF 0
#define DEFAULTMAXOFF 20
#define DEFAULTMINMYEARS 0
#define DEFAULTVARYOFFAGE 1
#define DEFAULTMINNUMMAPS 3
#define DEFAULTEXPORTFITS 0
#define DEFAULTDATAINTERVALS 100
#define DEFAULTSUPRESSCDF 0
#define DEFAULTEXTRAPOLATEMODEL 0
#define DEFAULTEXTRAPOLATIONFREQUENCY 1e10
#define DEFAULTLARGETXT 0
#define DEFAULTFORCEERRORTYPE 0

char **commandcompletion(const char *, int, int);
char *commandgenerator(const char *, int );

// Array of commands for autocomplete
#include "commandlist.h"
#include "commandlist_functions.h"


int main(int argc, char *argv[]) {

  // Define variables

  int c, i, a, j, o, p, q, r, *imgnum, *xdim, *ydim, **fluxorder, ***regionarray, *regnumber, tmp_currentset, tmp_signaltonoise, tmp_searcharea, *regionsset, dataselect, mapselect, *setaveraged, **regionsize, yesno, setvalues, selectlog, tmp_skip, *curvatureset, looplimit, tmp_plotres, tmp_poly, tmp_symbol, *jpfitmemset, tmp_ageresolution, tmp_myears, tmp_normresolution, *specindexset, numbermodelfluxes, *jpmodelset, *kpmodelset, tmp_modelres, tmp_levels, currentmap, tmp_contourlevels , tmp_contourcolour, tmp_contourlinestyle, tmp_firstcontour, *allocarray, map1, map2, map3, *kpfitmemset, tmp_posmap, tmp_scaletype, tmp_beamposx, tmp_beamposy, dataexporttype, confirm, tmp_modelmyears, modelloop, tmp_model, **mininjectset, tmp_injectinterval, gtlt, changearange, looplimit2, doallmaps, tmp_instrument, *jptribfitmemset, *jptribmodelset, *kptribfitmemset, *kptribmodelset, errorpm, **setmodelres, checkthejpmemory, checkthekpmemory, checkthejptribmemory, headerloop, strfound, line_num, foundheader, importcompression, tmp_indexcalctype, tmp_numdatapoints, *tmp_fluxorder, *specindex_modelres, *specindextype, tmp_minmodelmyears, tmp_minmyears, tmp_modelminoff, tmp_modelmaxoff, tmp_minoff, tmp_maxoff, *imgnum_CI, **fluxorder_CI, regnumber_CI, **upperlimits, *numupperlimits, mapopt, scaleopt, singlemap, t, freeparams, exactage, tmp_dataintervals, count, tmp_extrapolatemodel, *border, istart, iend, tmp_pltcolour, *specindex_errortype;


  int cmdnum = -1;
  int mainloop = 1;
  int numdatasets = 0;
  int arraysallocated = 0;
  int currentset = 0;
  int setnumber = 9999;
  int titles = DEFAULTITLES;
  int labels = DEFAULTLABELS;
  int axismarks = DEFAULTAXISMARKS;
  int regmemallocated = 0;
  int averegions = DEFAULTAVEREGIONS;
  int signaltonoise = DEFAULTSIGNOISE;
  int searcharea = DEFAULTSEARCHAREA;
  int printindex = DEFAULTPRINTINDEX;
  int specallocated = 0;
  int fixedregnum = 0;
  int fixedregref = 0;
  int extwincontrol = 0;
  int validentry = 0;
  int validentry2 = 0;
  int escape = 0;
  int escape2 = 0;
  int doregaveraging = 0;
  int checkuseall = 0;
  int uselogx = DEFAULTUSELOGX;
  int uselogy = DEFAULTUSELOGY;
  int autoscalex = DEFAULTAUTOSCALEX;
  int autoscaley = DEFAULTAUTOSCALEY;
  int skip = 1;
  int changevalues = 0;
  int model = 1;
  int printinc = DEFAULTPRINTINC;
  int curvearrayset = 0;
  int exportversion = 0;
  int poly = DEFAULTPOLY;
  int plotres = DEFAULTPLOTRES;
  int revcurve = DEFAULTREVCURVE;
  int printpoly = DEFAULYPRINTPOLY;
  int symbol = DEFAULTSYMBOL;
  int specres = DEFAULTSPECRES;
  int fluxres = DEFAULTFLUXRES;
  int ageresolution = DEFAULTAGERES;
  int myears = DEFAULTMYEARS;
  int normresolution = DEFAULTNORMRES;
  int jpmodelmemoryset = 0;
  int kpmodelmemoryset = 0;
  int jptribmodelmemoryset = 0;
  int kptribmodelmemoryset = 0;
  int modelres = DEFAULTMODELRES;
  int casadata = DEFAULTCASADATA;
  int levels = DEFAULTLEVELS;
  int printresults = DEFAULTPRINTRESULTS;
  int smoothmaps = DEFAULTSMOOTHMAPS;
  int contourlevels = DEFAULTCONLEVELS;
  int contourcolour = DEFAULTCONCOLOUR;
  int contourlinestyle = DEFAULTCONLINESTYLE;
  int contours = DEFAULTCONTOURS;
  int export = DEFAULTEXPORT;
  int warningcheck = 0;
  int firstcontour = DEFAULTFIRSTCONTOUR;
  int contourlogs = DEFAULTCONLOGS;
  int fluxlogs = DEFAULTFLUXLOGS;
  int fluxcut = DEFAULTFLUXCUT;
  int swapcoloursign = DEFAULTCOLOURSIGN;
  int alwayszeroage = DEFAULTALWAYSZEROAGE;
  int paramlabels = DEFAULTPARAMLABELS;
  int printreject = DEFAULTPRINTREJECT;
  int posmap = DEFAULTPOSMAP;
  int scaletype = DEFAULTSCALETYPE;
  int beamposx = DEFAULTBEAMPOSX;
  int beamposy = DEFAULTBEAMPOSY;
  int plotbeam = DEFAULTPLOTBEAM;
  int usereduced = DEFAULTUSEREDUCED;
  int modelmyears = DEFAULTMODELMYEARS;
  int injectinterval = DEFAULTINJECTINTERVAL;
  int instrument = DEFAULTINSTRUMENT;
  int xshift = DEFAULTXSHIFT;
  int yshift = DEFAULTYSHIFT;
  int exportcompression = DEFAULTEXPORTCOMPRESSION;
  int expectedheaders = DEFAULTEXPECTEDHEADERS;
  int expectedmodelheaders = DEFAULTEXPECTEDMODELHEADERS;
  int indexcalctype = DEFAULTINDEXCALCTYPE;
  int minmodelmyears = DEFAULTMINMODELMYEARS;
  int modelminoff = DEFAULTMODELMINOFF;
  int modelmaxoff = DEFAULTMODELMAXOFF;
  int minoff = DEFAULTMINOFF;
  int maxoff = DEFAULTMAXOFF;
  int minmyears = DEFAULTMINMYEARS;
  int varyoffage = DEFAULTVARYOFFAGE;
  int minnummaps = DEFAULTMINNUMMAPS;
  int exportfits = DEFAULTEXPORTFITS;
  int dataintervals = DEFAULTDATAINTERVALS;
  int firstplot = 1;
  int suppresscdf = DEFAULTSUPRESSCDF;
  int extrapolatemodel = DEFAULTEXTRAPOLATEMODEL;
  int largetxt = DEFAULTLARGETXT;
  int forceerrortype = DEFAULTFORCEERRORTYPE;


  char cmdtxt[MAXCMDLENGTH], setname[MAXDATASETS][MAXCMDLENGTH], setreg[MAXDATASETS][MAXCMDLENGTH], setbg[MAXDATASETS][MAXCMDLENGTH], dirname[MAXCMDLENGTH], regname[MAXCMDLENGTH], bgname[MAXCMDLENGTH], maptitle[MAXTITLELENGTH], **settarget, tmp_maptarget[MAXCMDLENGTH], output[MAXCMDLENGTH], imageloc[MAXCMDLENGTH], imagetype[32], timebuff[32], modelbuff[32], posnegbuff[32], tmp_imageloc[MAXCMDLENGTH], tmp_imagetype[MAXCMDLENGTH], tmp_filename[MAXCMDLENGTH], filename[MAXCMDLENGTH], dataloc[MAXCMDLENGTH], dataexportext[16], datasettype[128], tmp_dataloc[MAXCMDLENGTH], importfilename[MAXCMDLENGTH], tmp_filename2[MAXCMDLENGTH], *tmp_char, citarget[MAXCMDLENGTH], dataoutput[MAXCMDLENGTH], dataoutput_ext[MAXCMDLENGTH], *cmdbuffer, *endptr;


  // Set the default image export location and type

  sprintf(imageloc, DEFAULTIMAGELOC);
  sprintf(imagetype, DEFAULTIMAGETYPE);
  sprintf(dataloc, DEFAULTDATALOC);
  sprintf(dataexportext, DEFAULTDATAEXPORTEXT);

  float **frequency, **frequencyobs, *bmaj, *bmin, *bpa, *beamarea, **rms, ****flux, ***regflux, ***fluxerror, **regionlocx, **regionlocy, tmp_zoom, tmp_maptomap, tmp_hotpixels,  *minsi, *maxsi, **alpha, **maxflux, **minflux, **regmaxflux, **regminflux, ***curvearray, *mincurve, *maxcurve, *x, *y, *coeff, tmp_specrange, tmp_normrange, tmp_specres, tmp_fluxres, **jpchisquared, **jpmodelalpha, ***jpmodelflux, **jpbestnorm, **jpbestage, tmp_onsourcenoise, tmp_minmodelfreq, tmp_maxmodelfreq, *tmp_modelflux, *tmp_regflux, *tmp_modelfreq, *tmp_fluxfreq, *tmp_fluxerror, tmp_minflux, tmp_maxflux, *tmp_fluxerrorplus, *tmp_fluxerrorminus, **realloc_jpchisquared, **realloc_jpbestage, **realloc_jpmodelalpha, **realloc_jpbestnorm, ***realloc_jpmodelflux, tmp_usercut, tmp_wedgemultiplier, *colour1, *colour2, mincolour1, maxcolour1, mincolour2, maxcolour2, **kpchisquared, **kpmodelalpha, ***kpmodelflux, **kpbestnorm, **kpbestage, **realloc_kpchisquared, **realloc_kpbestage, **realloc_kpmodelalpha, **realloc_kpbestnorm, ***realloc_kpmodelflux, **ra, **dec, *cellsize, tmp_flt_inject, tmp_field, tmp_bestage, tmp_chisquared, ***inj_sumchisquared, ***inj_regchisquared,  ***inj_regbestinject, tmp_inject, bestinject1, bestinject2, bestchi1, bestchi2, tmp_mininject, tmp_maxinject, **mininjectstore, **maxinjectstore, **fluxcalerror, errorchange, tmp_flaterror, **jptribchisquared, **jptribmodelalpha, ***jptribmodelflux, **jptribbestnorm, **jptribbestage, **realloc_jptribchisquared, **realloc_jptribbestage, **realloc_jptribmodelalpha, **realloc_jptribbestnorm, ***realloc_jptribmodelflux, **kptribchisquared, **kptribmodelalpha, ***kptribmodelflux, **kptribbestnorm, **kptribbestage, **realloc_kptribchisquared, **realloc_kptribbestage, **realloc_kptribmodelalpha, **realloc_kptribbestnorm, ***realloc_kptribmodelflux, ***bgbuff, **jpageerrorsplus, **kpageerrorsplus, **jptribbleageerrorsplus, **realloc_jpageerrorsplus, **realloc_kpageerrorsplus, **realloc_jptribbleageerrorsplus, *tmp_ageerrorsplus, **jpageerrorsminus, **kpageerrorsminus, **jptribbleageerrorsminus, **realloc_jpageerrorsminus, **realloc_kpageerrorsminus, **realloc_jptribbleageerrorsminus, *tmp_ageerrorsminus, **specindexchisquared, tmp_modfreq, tmp_exampleredshift, *cibcmb, *ciageingb, *cibreakon, *cibreakoff, *ciconflvl, **tmp_integfluxes, *tmp_integfreq, **tmp_fluxerr, tmp_holdintflux, tmp_holdintfreq, tmp_holdfluxerr, tmp_curminfreq, tmp_curmaxfreq, *tmp_chisq, tmp_redshift, *tmp_bage, *tmp_modalpha, **tmp_modflux, *tmp_bnorm, *redshift, *fluxerrorplus_CI, *fluxerrorminus_CI, *tmp_regflux_CI, *modelfreq_CI, *fluxfreq_CI, ***regflux_CI, **frequency_CI, *cichisquared, *cibestage, *cibestoff, ***fluxerror_CI, *cimodelalpha, **cimodelflux, **log_cimodelflux, *cibestnorm, *CI_redshift, *tmp_offerrorsplus_CI, *tmp_offerrorsminus_CI, *tmp_onerrorsplus_CI, *tmp_onerrorsminus_CI, minfluxscalefreq, maxfluxscalefreq, factor1, factor2, singlefactor, factorfreq1, factorfreq2, **delt1, **delt2, **crpix1, **crpix2, **crota1, **crota2, **eq, exactmyears, exactoff, tmp_sigma, tmp_extrapolationfrequency, tmp_fracborder, zoomedx, **specindexerror;
  

  // float tmp_userminx, tmp_usermaxx, tmp_userminy, tmp_usermaxy, tmp_bestguess; // From old input method

  float zoom = DEFAULTZOOM;
  float maptomap = DEFAULTMAPTOMAP;
  float hotpixels = DEFAULTHOTPIXELS;
  float bestguess = DEFAULTBESTGUESS;
  float specrange = DEFAULTSPECRANGE;
  float normrange = DEFAULTNORMRANGE;
  float minmodelfreq = DEFAULTMINMODELFREQ;
  float maxmodelfreq = DEFAULTMAXMODELFREQ;
  float onsourcenoise = DEFAULTONSOURCENOISE;
  float usercut = DEFAULTUSERCUT;
  float wedgemultiplier = DEFAULTWEDGE;
  float userminx = 0.0;
  float usermaxx = 0.0;
  float userminy = 0.0;
  float usermaxy = 0.0;
  float mininject = DEFAULTMININJECT;
  float maxinject = DEFAULTMAXINJECT;
  float flaterror = DEFAULTERRFLAT;
  float exampleredshift = DEFAULTEXAMPLEREDSHIFT;
  float sigmalevel = DEFAULTSIGMA;
  float extrapolationfrequency = DEFAULTEXTRAPOLATIONFREQUENCY;
  float frac_border = DEFAULTFRACBORDER/100; // Default is in percentage format


  double tmp_fieldstrength, tmp_dbl_inject, tmp_power, *jpinject, *kpinject, *jpfield, *kpfield, interval, siglevel68, siglevel90, siglevel95, siglevel99, siglevel995, siglevel999, siglevel9999, dof, tmp_usr_gmin, tmp_usr_gmax, *jptribinject, *kptribinject, *jptribfield, *kptribfield, **specindexintercept, ***specindex_modelflux, ***specindex_modelerror;

  double fieldstrength = DEFAULTFIELDSTRENGTH;
  double inject = DEFAULTINJECT;
  double power = (2*DEFAULTINJECT) + 1;
  double usr_gmin = DEFAULTGMIN;
  double usr_gmax = DEFAULTGMAX;

  FILE *filestream, *importfile;
  DIR *dir;

  // Set up variables and structs for times and dates

  time_t currenttime;
  struct tm * time_struct;

  // Allocate memory needed from the outset
  regionsset = (int *)calloc(MAXDATASETS, sizeof(int));
  settarget = (char **)calloc(MAXDATASETS,sizeof(char *));
  jpfitmemset = (int *)calloc(MAXDATASETS, sizeof(int));
  curvatureset = (int *)calloc(MAXDATASETS, sizeof(int));
  specindexset = (int *)calloc(MAXDATASETS, sizeof(int));
  jpmodelset = (int *)calloc(MAXDATASETS, sizeof(int));
  kpmodelset = (int *)calloc(MAXDATASETS, sizeof(int));
  kpfitmemset = (int *)calloc(MAXDATASETS, sizeof(int));
  jpinject = (double *)calloc(MAXDATASETS, sizeof(double));
  kpinject = (double *)calloc(MAXDATASETS, sizeof(double));
  jpfield = (double *)calloc(MAXDATASETS, sizeof(double));
  kpfield = (double *)calloc(MAXDATASETS, sizeof(double));
  inj_sumchisquared = (float ***)calloc(MAXDATASETS, sizeof(float **));
  inj_regchisquared = (float ***)calloc(MAXDATASETS, sizeof(float ***));
  inj_regbestinject = (float ***)calloc(MAXDATASETS, sizeof(float **));
  mininjectstore = (float **)calloc(MAXDATASETS, sizeof(float *));
  maxinjectstore = (float **)calloc(MAXDATASETS, sizeof(float *));
  mininjectset = (int **)calloc(MAXDATASETS, sizeof(int *));
  fluxcalerror = (float **)calloc(MAXDATASETS, sizeof(float *));
  jptribfitmemset = (int *)calloc(MAXDATASETS, sizeof(int));
  jptribmodelset = (int *)calloc(MAXDATASETS, sizeof(int));
  jptribinject = (double *)calloc(MAXDATASETS, sizeof(double));
  jptribfield = (double *)calloc(MAXDATASETS, sizeof(double));
  kptribfitmemset = (int *)calloc(MAXDATASETS, sizeof(int));
  kptribmodelset = (int *)calloc(MAXDATASETS, sizeof(int));
  kptribinject = (double *)calloc(MAXDATASETS, sizeof(double));
  kptribfield = (double *)calloc(MAXDATASETS, sizeof(double));
  setmodelres = (int **)calloc(MAXDATASETS, sizeof(int *));
  specindex_modelres = (int *)calloc(MAXDATASETS, sizeof(int));
  specindextype = (int *)calloc(MAXDATASETS, sizeof(int));
  border = (int *)calloc(MAXDATASETS, sizeof(int));



  for (i=0; i < MAXDATASETS; i++) {
    inj_sumchisquared[i] = (float **)calloc(NUMBEROFMODELS+1, sizeof(float *));
    inj_regchisquared[i] = (float **)calloc(NUMBEROFMODELS+1, sizeof(float **));
    inj_regbestinject[i] = (float **)calloc(NUMBEROFMODELS+1, sizeof(float *));
    mininjectstore[i] = (float *)calloc(NUMBEROFMODELS+1, sizeof(float));
    maxinjectstore[i] = (float *)calloc(NUMBEROFMODELS+1, sizeof(float));
    mininjectset[i] = (int *)calloc(NUMBEROFMODELS+1, sizeof(int));
    fluxcalerror[i] = (float *)calloc(MAXNUMMAPS, sizeof(float));
    setmodelres[i] = (int *)calloc(NUMBEROFMODELS+1, sizeof(int));
    border[i] = DEFAULTBORDER; // Set border to default value as 0 just in case. 0 causes a plot failure.
  }

  // Declare the function prototypes

  void cpgclos(void);
  int cpgopen(const char *device);
  void cpgask(int);

  int chomp(char extstr[]);
  int help();
  int syscom(const char command[]);

  int list(const int indexnumber, const int imgnum, const char dirname[], const char regname[], const char bgname[]);

  int listprops(const int indexnumber, const int imgnum, const char dirname[], const char regname[], const char bgname[], const int xdim, const int ydim, const float bmaj, const float bmin, const float beamarea);

  int load(char dirname[], int currentset, int numdatasets, int *imgnum, int *xdim, int *ydim, float *bmaj, float *bmin, float *bpa, float *beamarea, float **frequency, float **frequencyobs, char regname[], char bgname[], float ****flux, float **rms, float sigmalevel, int **fluxorder, float zoom, int *border, int titles, int labels, int axismarks, char tmp_maptarget[], float **yaxismax, float **yaxismin, int casadata, int fluxlogs, int fluxcut, float usercut, float wedgemultiplier, float **ra, float **dec, float **delt1, float **delt2, float **crpix1, float **crpix2, float **crota1, float **crota2, float **eq, float *cellsize, float scaletype, int posmap, int beamposx, int beamposy, int plotbeam, float ***bgbuff, int xshift, int yshift, float *redshift, int largetxt, float frac_border);

  int adaptiveregions(int *imgnum, int *xdim, int *ydim, float ****flux, float minflux, float **rms, int maxarea, int ***regionarray, float ***regflux, int *regnumber, float ***fluxerror, int averaged, float **regionlocx, float **regionlocy, int **fluxorder, int currentset, float m2m, float hotpixels, float sigmalevel, int titles, int labels, int axismarks, int **regionsize, float **regmaxflux, float **regminflux, float *fluxcalerror, float onsourcenoise, int regionsset, int border, float zoom, float ra, float dec, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, int xshift, int yshift, int smoothmaps, int largetxt);

  int fluxmap(float ****flux, float **rms, int *xdim, int *ydim, int mapnumber, float **maxflux, float **minflux, int currentset, int border, float zoom, float sigmalevel, char fluxmaptitle[128], int titles, int labels, int axismarks, int extwincontrol, int smoothmaps, int fluxlogs, int fluxcut, float usercut, float wedgemultiplier, char *output, float ra, float dec, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, int xshift, int yshift, int largetxt, int export);

  int regionfluxmap(float ***regflux, int *xdim, int *ydim, int mapnumber, float **maxflux, float **minflux, int currentset, int border, float zoom, char fluxmaptitle[], int titles, int labels, int axismarks, int ***regionarray, int extwincontrol, int **regionsize, int *regnumber, int doregaveraging, int smoothmaps, char *output, float ra, float dec, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, int xshift, int yshift, int largetxt, int export);

  int calcspecindexp2p(float **regflux, float *specindex, int numregions, int map1, int map2, float freq1, float freq2, int printindex, float *minsi, float *maxsi, float *rms);

  int calcspecindex(float **regflux,  int **regionarray, float *specindex, int numregions, int imgnum, float *frequency, int printindex, float *minsi, float *maxsi);

  int calcspecindexweightedgsl(float **regflux, int **regionarray, float *specindex, double *specindexintercept, float *specindexchisquared, double **specindex_modelflux, double **specindex_modelerror, float **fluxerror, int numregions, int imgnum, float *frequency, int printindex, float *minsi, float *maxsi, int modelres, int *fluxorder, float *specindexerror, int *specindex_errortype, int forceerrortype);

  int applyfixedregions(int *imgnum, int *xdim, int *ydim, float ****flux, float **rms, int ***regionarray, float ***regflux, int *regnumber, float ***fluxerror, int averaged, float **regionlocx, float **regionlocy, int titles, int labels, int axismarks, int fixedregnum, int fixedregref, int *regionsset, float *fluxcalerror, float onsourcenoise, float ra, float dec, float cellsize, int border, float zoom, int xshift, int yshift, int smoothmaps, int scaletype, int largetxt, int export);

  int plotpolyfit(int *regnumber, int currentset, int *imgnum, float **frequency, int **fluxorder, int plotres, float ***curvearray, float *mincurve, float *maxcurve, int poly, float userminx, float usermaxx, float userminy, float usermaxy, int autoscalex, int autoscaley, int skip, int titles, int labels, int axismarks, char *target, int extwincontrol, char *output);

  int plotfluxbyregion(float ***regflux, int *imgnum, int **fluxorder, float **frequency, float minfreq, float maxfreq, float **minflux, float **maxflux, int *regnumber, float ***fluxerror, char maptitle[], int titles, int labels, int axismarks, int extwincontrol, int currentset, int uselogx, int uselogy, int autoscalex, int autoscaley, float userminx, float usermaxx, float userminy, float usermaxy, int printincreasing, int skip, int symbol, char *output);

  int plottotalflux(float **regflux, float *frequency, int *fluxorder, float *rms, int imgnum, int numregions, int symbol, char *target, int extwincontrol, float *fluxcalerror, float bestguess, float specrange, float fluxrange, int resolution, int specindexres, int fluxresolution, int *regionsize, float onsourcenoise, char *settarget, int titles, int labels, int axismarks, char *output);

  int plotspecindex(float **specindex, int numregions, int *xdim, int *ydim, char *settarget, float *minsi, float *maxsi, int ***regionarray, int currentset, char *specmaptitle, float zoom, int border, int titles, int labels, int axismarks, int symbol, int smoothmaps, float **flux, int extwincontrol, int contours, int contourlevels, int contourcolour, int contourlinestyle, int firstcontour, int contourlogs, char *output, float ra, float dec, float delt1, float delt2, float crpix1, float crpix2, float crota1, float crota2, float eq, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, float wedgemultiplier, int xshift, int yshift, float minfreq, float maxfreq, int export, int exportfits, int largetxt);

  int curvature(int n, float *x, float *y, int poly, float *coeff, int printpoly);

  int plotmodel(float minf, float maxf, double inject, double fieldstrength, int model, double usr_gmin, double usr_gmax, int minmodelmyears, int modelmyears, char *output, int titles, float redshift, int skip);

  int spectralageingmodels(float **regflux, int imgnum, int *fluxorder, float *frequency, int numregions, float *chisquared, float *bestage, float **fluxerror, float *modelalpha, float **modelflux, double usr_gmin, double usr_gmax, int minmyears, int myears, int ageresolution, int levels, float *bestnorm, int modelres, int printresults, double inject, double fieldstrength, int model, int printreject, float redshift, float *ageerrorsplus, float *ageerrorsminus, int suppresscdf);

  int plotcurvagainstspec(float **curvature, float *specindex, int regnumber, int symbol, float xmin, float xmax, float ymin, float ymax, int poly, int autoscalex, int autoscaley, float userminx, float usermaxx, float userminy, float usermaxy, char *settarget, int titles, int labels, int axismarks, int extwincontrol);

  int plotmodelvsflux(float *modelflux, float *regflux, float *fluxerrorplus, float *fluxerrorminus, int symbol, int nummodelpts, int imgnum, float *fluxfrequency, float *modelfrequency, int autoscalex, int autoscaley, float userminx, float usermaxx, float userminy, float usermaxy, float minfreq, float maxfreq, float minflux, float maxflux, char *settarget, int titles, int labels, int axismarks, double inject, double fieldstrength, float bestage, float bestoff, float chisquared, int paramlabels, int model, int usereduced, int ul, int *upperlimits, int largetxt);

  int spectralagemap(float *bestage, int **regionarray, float **flux, int xaxismax, int yaxismax, char *settarget, int border, float zoom, int titles, int labels, int axismarks, int extwincontrol, int smoothmaps, int contours, int contourlevels, int contourcolour, int contourlinestyle, int firstcontour, int contourlogs, char *output, int alwayszeroage, double inject, double fieldstrength, int paramlabels, int model, float ra, float dec, float delt1, float delt2, float crpix1, float crpix2, float crota1, float crota2, float eq, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, float wedgemultiplier, int xshift, int yshift, int export, int exportfits, int largetxt);

  int chisquaredmap(float *chisquared, int **regionarray, int xaxismax, int yaxismax, char *settarget, int imgnum, int border, float zoom, char maptitle[], int titles, int labels, int axismarks, int extwincontrol, int smoothmaps, int contours, int contourlevels, int contourcolour, int contourlinestyle, float **flux, int firstcontour, int contourlogs, char *output, double inject, double fieldstrength, int paramlabels, int model, float ra, float dec, float delt1, float delt2, float crpix1, float crpix2, float crota1, float crota2, float eq, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, float wedgemultiplier, int usereduced, int xshift, int yshift, int export, int exportfits, int largetxt);

  int colourcolour(float *colour1, float *colour2, float mincolour1, float maxcolour1, float mincolour2, float maxcolour2, float map1freq, float map2freq, float map3freq, char *settarget, int regnumber, int axismarks, int titles, int labels, int symbol, int extwincontrol, char *output, int swapsign);

  int singleregion(int *imgnum, int *xdim, int *ydim, float ****flux, float minflux, float **rms, int maxarea, int ***regionarray, float ***regflux, int *regnumber, float ***fluxerror, int averaged, float **regionlocx, float **regionlocy, int **fluxorder, int currentset, float m2m, float hotpixels, float sigmalevel, int titles, int labels, int axismarks, int **regionsize, float **regmaxflux, float **regminflux, float *fluxcalerror, float onsourcenoise, int regionsset, int border, float zoom, float ra, float dec, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, int xshift, int yshift, int largetxt);

  int spectralageing_min(float **regflux, int imgnum, int *fluxorder, float *frequency, int numregions, float **fluxerror, double usr_gmin, double usr_gmax, int minmyears, int myears, int ageresolution, int levels, double inject_loop, double fieldstrength, int model, float *sumchisquared, float *regchisquared, float *regbestinject, float redshift);

  int fluxerrors(float *fluxcalerror, float *frequency, int imgnum, int instrument, float flaterror);

  int mergeimages(int casadata, char *imageloc);

  int diffmap(int casadata, char *imageloc, int ***regionarray, int *regionsset, int *reg_xdim, int *reg_ydim, int numdatasets, int *imgnum, char setname[][MAXCMDLENGTH], char setreg[][MAXCMDLENGTH], char setbg[][MAXCMDLENGTH]);

  int scalefluxes(int currentset, int mapopt, float minfreq, float maxfreq, int singlemap, int scaleopt, float factorfreq1, float factorfreq2, float factor1, float factor2, float singlefactor, float ****flux, int xdim, int ydim, int *fluxorder, int imgnum, float *frequency, float ***bgbuff, float **rms, float beamarea);

  int mapfrommodel(char *imageloc, char *dirname, int xdim, int ydim, char *maptarget, double usr_gmin, double usr_gmax, double fieldstrength, int model, float *bestnorm, float *bestage, int **regionarray, float inject, float beamarea, float redshift);

  int errormap(float *ageerrors, int **regionarray, float **flux, int xaxismax, int yaxismax, char *settarget, int border, float zoom, int titles, int labels, int axismarks, int extwincontrol, int smoothmaps, int contours, int contourlevels, int contourcolour, int contourlinestyle, int firstcontour, int contourlogs, char *output, double inject, double fieldstrength, int paramlabels, int model, float ra, float dec, float delt1, float delt2, float crpix1, float crpix2, float crota1, float crota2, float eq, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, float wedgemultiplier, int xshift, int yshift, int errorpm, int export, int exportfits, int largetxt);

  int fullexport(char *filename, int exportcompression, char *settarget, int imgnum, float redshift, int xdim, int ydim, float bmaj, float bmin, float bpa, float beamarea, float *ra, float *dec, float *eq, float *delt1, float *delt2, float *crpix1, float *crpix2, float *crota1, float *crota2, float cellsize, char *setname, char *setbg, char *setreg, int posmap, char *dataloc, int regionsset, int jpmodelset, int kpmodelset, int jptribmodelset, int *mininjectset, float *frequency, float *frequencyobs, int *fluxorder, float *rms, float ***flux, float *minflux, float *maxflux, float **bgbuff, float *fluxcalerror, int setaveraged, int regnumber, float *regionlocx, float *regionlocy, float *regminflux, float *regmaxflux, float **fluxerror, int **regionarray, float **regflux, float *mininjectstore, float *maxinjectstore, float **inj_sumchisquared, int *setmodelres, float jpfield, float jpinject, float *jpbestage, float *jpbestnorm, float *jpchisquared, float *jpmodelalpha, float **jpmodelflux, float *jpageerrorsplus, float *jpageerrorsminus, float kpfield, float kpinject, float *kpbestage, float *kpbestnorm, float *kpchisquared, float *kpmodelalpha, float **kpmodelflux, float *kpageerrorsplus, float *kpageerrorsminus, float jptribfield, float jptribinject, float *jptribbestage, float *jptribbestnorm, float *jptribchisquared, float *jptribmodelalpha, float **jptribmodelflux, float *jptribbleageerrorsplus, float *jptribbleageerrorsminus, int *regionsize, float **inj_regbestinject, float **inj_regchisquared);

  int ciageingmodels(float **regflux, int imgnum, int *fluxorder, float *frequency, int passingreg, float *chisquared, float *bestage, float *bestoff, float **fluxerror, float *modelalpha, float **modelflux, double usr_gmin, double usr_gmax, int myears, int minmyears, int minoff, int maxoff, int ageresolution, int levels, float *bestnorm, int modelres, int printresults, double inject, double fieldstrength, int model, int printreject, float redshift, float *offerrorsplus, float *offerrorsminus, float *onerrorsplus, float *onerrorsminus, int ul, int *upperlimits, float *bcmb, float *passageingb, float *breakon, float *breakoff, float *conflvl, int *suppress, char *index, int export, int suppresscdf, int extrapolatemodel, float extrapolationfrequency, float fixedage);

  int injectchibyregion(float *inj_regchisquared, int **regionarray, float **flux, int xaxismax, int yaxismax, char *settarget, int border, float zoom, int titles, int labels, int axismarks, int extwincontrol, int smoothmaps, int contours, int contourlevels, int contourcolour, int contourlinestyle, int firstcontour, int contourlogs, char *output, int model, float ra, float dec, float delt1, float delt2, float crpix1, float crpix2, float crota1, float crota2, float eq, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, float wedgemultiplier, int xshift, int yshift, int export, int exportfits, int usereduced, int dof, int largetxt);

  int injectbyregion(float *inj_regbestinject, int **regionarray, float **flux, int xaxismax, int yaxismax, char *settarget, int border, float zoom, int titles, int labels, int axismarks, int extwincontrol, int smoothmaps, int contours, int contourlevels, int contourcolour, int contourlinestyle, int firstcontour, int contourlogs, char *output, int model, float ra, float dec, float delt1, float delt2, float crpix1, float crpix2, float crota1, float crota2, float eq, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, float wedgemultiplier, int xshift, int yshift, int export, int exportfits, int largetxt);

  int plotfluxvsspecindex(float *regflux, float *fluxerrorplus, float *fluxerrorminus, int symbol, int imgnum, float *fluxfrequency, int autoscalex, int autoscaley, float userminx, float usermaxx, float userminy, float usermaxy, float minfreq, float maxfreq, float minflux, float maxflux, char *settarget, int titles, int labels, int axismarks, int paramlabels, double *specindex_modelflux, double *specindex_modelerror, int modelres, int *fluxorder, float specindex, int largetxt, float specindexerror, int specindex_errortype);

  int specindexerrormap(float **specindexerror, int numregions, int *xdim, int *ydim, char *settarget, float *minsi, float *maxsi, int ***regionarray, int currentset, char *maptitle, float zoom, int border, int titles, int labels, int axismarks, int symbol, int smoothmaps, float **flux, int extwincontrol, int contours, int contourlevels, int contourcolour, int contourlinestyle, int firstcontour, int contourlogs, char *output, float ra, float dec, float delt1, float delt2, float crpix1, float crpix2, float crota1, float crota2, float eq, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, float wedgemultiplier, int xshift, int yshift, float minfreq, float maxfreq, int export, int exportfits, int largetxt);

  int importcoords(char *importfilename, int currentset, float **ra, float **dec, int imgnum);

  void cpgscr(int ci, float cr, float cg, float cb);

  char *usage="Syntax: %s [-h help] [-v version]\n\n";
  
  // Process the command line options

  while ((c=getopt(argc, argv, "hv"))!=EOF) {
    switch (c) {
    case 'h':
      printf("\n");
      fprintf(stderr,usage,argv[0]);
      printf("Details of all functions and their usage can be found in the BRATS cookbook and within BRATS itself via the \"help\" command.\nFor further information, support, or to download the latest version, please visit: http://www.askanastronomer.co.uk/brats\n\n");
      exit(0);
      break;
    case 'v':
      printf("\nBroadband Radio Astronomy ToolS (BRATS) - Version %s\n", VERSION);
      printf("Written by: Jeremy J. Harwood (Jeremy.Harwood@physics.org)\n");
      printf("For further information and support, please visit:\n");
      printf("http://www.askanastronomer.co.uk/brats\n");
      printf("\nIf you have made use of this software, please cite:\n");
      printf(" - Harwood et al., 2013, MNRAS, 435, 3353\n");
      printf(" - Harwood et al., 2015, MNRAS, 454, 3403\n\n");
      exit(0);
      break;
    default:
      fprintf(stderr,usage,argv[0]);
      exit(1);
      break;
    }
  }


  // Output the startup details and version info
  printf("\n");
  printf("##############################################################################\n");
  printf("#        Welcome to the Broadband Radio Analysis ToolS (BRATS) Software      #\n");
  printf("#                              Version %s                               #\n", VERSION);
  printf("#                                                                            #\n");
  printf("#            For further information and support, please visit:              #\n");
  printf("#                 http://www.askanastronomer.co.uk/brats                     #\n");
  printf("#                                                                            #\n");
  printf("#         Developed by Jeremy J. Harwood (Jeremy.Harwood@physics.org)        #\n");
  printf("#                                                                            #\n");
  printf("#             If you have made use of this software please cite              #\n");
  printf("#                  Harwood et al., 2013, MNRAS, 435, 3353                    #\n");
  printf("#                  Harwood et al., 2015, MNRAS, 454, 3403                    #\n");
  printf("##############################################################################\n");
  printf("\n");

  // Auto complete commands, paths, and file names
  rl_attempted_completion_function = commandcompletion;

  // Main program command loop
  while (mainloop == 1) {

    // Set the command number to 0
    cmdnum = 0;

    // Print the command promt and capture any input
    //fputs("BRATS: ", stdout);// Save in case people have problems with readline
    //fflush(stdout);
    //fgets(cmdtxt, sizeof(cmdtxt), stdin); // Save in case people have problems with readline

    cmdbuffer = readline("BRATS: ");

    // If it's not blank, add to the history
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    strcpy(cmdtxt, cmdbuffer);

    // Format out any return trailing white space
    chomp(cmdtxt);

    for (tmp_char = cmdtxt; *tmp_char != '\0'; ++tmp_char) {
      *tmp_char = tolower(*tmp_char);
    }

    if (strlen(cmdtxt) == 0) {
      cmdnum = 0;
    }
    else if ((strstr(cmdtxt, "help") != NULL) && (strlen(cmdtxt) == 4)) {
      cmdnum = 1;
    }
    else if ((strstr(cmdtxt, "quit") != NULL) && (strlen(cmdtxt) == 4)) {
      cmdnum = 2;
    }
    else if ((strstr(cmdtxt, "exit") != NULL) && (strlen(cmdtxt) == 4)) {
      printf("\nPlease use 'quit' to leave BRATS. This is to avoid accidently closing a session and losing data when attempting to exit a seperate terminal window.\n\n");
    }
    else if ((strstr(cmdtxt, "load") != NULL) && (strlen(cmdtxt) == 4)) {
      cmdnum = 3;
    }
    else if ((strstr(cmdtxt, "ls") != NULL) && (strlen(cmdtxt) == 2)) {
      cmdnum = 4;
    }
    else if ((strstr(cmdtxt, "shell") != NULL) && (strlen(cmdtxt) == 5)) {
      cmdnum = 5;
    }
    else if ((strstr(cmdtxt, "list") != NULL) && (strlen(cmdtxt) == 4)) {
      cmdnum = 6;
    }
    else if ((strstr(cmdtxt, "props") != NULL) && (strlen(cmdtxt) == 5)) {
      cmdnum = 7;
    }
    else if ((strstr(cmdtxt, "setregions") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 8;
    }
    else if ((strstr(cmdtxt, "zoom") != NULL) && (strlen(cmdtxt) == 4)) {
      cmdnum = 9;
    }
    else if ((strstr(cmdtxt, "border") != NULL) && (strlen(cmdtxt) == 6)) {
      cmdnum = 10;
    }
    else if ((strstr(cmdtxt, "sigma") != NULL) && (strlen(cmdtxt) == 5)) {
      cmdnum = 11;
    }
    else if ((strstr(cmdtxt, "titles") != NULL) && (strlen(cmdtxt) == 6)) {
      cmdnum = 12;
    }
    else if ((strstr(cmdtxt, "labels") != NULL) && (strlen(cmdtxt) == 6)) {
      cmdnum = 13;
    }
    else if ((strstr(cmdtxt, "axismarks") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 14;
    }
    else if ((strstr(cmdtxt, "maptomap") != NULL) && (strlen(cmdtxt) == 8)) {
      cmdnum = 15;
    }
    else if ((strstr(cmdtxt, "hotpixels") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 16;
    }
    else if ((strstr(cmdtxt, "averegions") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 17;
    }
    else if ((strstr(cmdtxt, "signaltonoise") != NULL) && (strlen(cmdtxt) == 13)) {
      cmdnum = 18;
    }
    else if ((strstr(cmdtxt, "searcharea") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 19;
    }
    else if ((strstr(cmdtxt, "specindex") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 20;
    }
    else if ((strstr(cmdtxt, "printindex") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 21;
    }
    else if ((strstr(cmdtxt, "fixedregions") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 22;
    }
    else if ((strstr(cmdtxt, "fluxmap") != NULL) && (strlen(cmdtxt) == 7)) {
      cmdnum = 23;
    }
    else if ((strstr(cmdtxt, "plotflux") != NULL) && (strlen(cmdtxt) == 8)) {
      cmdnum = 24;
    }
    else if ((strstr(cmdtxt, "autoscale") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 25;
    }
    else if ((strstr(cmdtxt, "printinc") != NULL) && (strlen(cmdtxt) == 8)) {
      cmdnum = 26;
    }
    else if ((strstr(cmdtxt, "uselogs") != NULL) && (strlen(cmdtxt) == 7)) {
      cmdnum = 27;
    }
    else if ((strstr(cmdtxt, "skip") != NULL) && (strlen(cmdtxt) == 4)) {
      cmdnum = 28;
    }
    else if ((strstr(cmdtxt, "fitpoly") != NULL) && (strlen(cmdtxt) == 7)) {
      cmdnum = 29;
    }
    else if ((strstr(cmdtxt, "plotres") != NULL) && (strlen(cmdtxt) == 7)) {
      cmdnum = 30;
    }
    else if ((strstr(cmdtxt, "setpoly") != NULL) && (strlen(cmdtxt) == 7)) {
      cmdnum = 31;
    }
    else if ((strstr(cmdtxt, "curvecon") != NULL) && (strlen(cmdtxt) == 8)) {
      cmdnum = 32;
    }
    else if ((strstr(cmdtxt, "totalflux") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 33;
    }
    else if ((strstr(cmdtxt, "fluxcalerror") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 34;
    }
    else if ((strstr(cmdtxt, "bestguess") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 35;
    }
    else if ((strstr(cmdtxt, "specrange") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 36;
    }
    else if ((strstr(cmdtxt, "normrange") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 37;
    }
    else if ((strstr(cmdtxt, "symbol") != NULL) && (strlen(cmdtxt) == 6)) {
      cmdnum = 38;
    }
    else if ((strstr(cmdtxt, "specres") != NULL) && (strlen(cmdtxt) == 7)) {
      cmdnum = 39;
    }
    else if ((strstr(cmdtxt, "fluxres") != NULL) && (strlen(cmdtxt) == 7)) {
      cmdnum = 40;
    }
    else if ((strstr(cmdtxt, "plotjpmodel") != NULL) && (strlen(cmdtxt) == 11)) {
      cmdnum = 41;
    }
    else if ((strstr(cmdtxt, "fitjpmodel") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 42;
    }
    else if ((strstr(cmdtxt, "onsource") != NULL) && (strlen(cmdtxt) == 8)) {
      cmdnum = 43;
    }
    else if ((strstr(cmdtxt, "minmodelfreq") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 44;
    }
    else if ((strstr(cmdtxt, "maxmodelfreq") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 45;
    }
    else if ((strstr(cmdtxt, "ageres") != NULL) && (strlen(cmdtxt) == 6)) {
      cmdnum = 46;
    }
    else if ((strstr(cmdtxt, "myears") != NULL) && (strlen(cmdtxt) == 6)) {
      cmdnum = 47;
    }
    else if ((strstr(cmdtxt, "normres") != NULL) && (strlen(cmdtxt) == 7)) {
      printf("\n This command is now redundant due to the use of a golden ratio search to full float accurancy\n\n");
      //cmdnum = 48;
    }
    else if ((strstr(cmdtxt, "plotcurvespec") != NULL) && (strlen(cmdtxt) == 13)) {
      cmdnum = 49;
    }
    else if ((strstr(cmdtxt, "plotmodelobs") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 50;
    }
    else if ((strstr(cmdtxt, "modelres") != NULL) && (strlen(cmdtxt) == 8)) {
      cmdnum = 51;
    }
    else if ((strstr(cmdtxt, "specagemap") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 52;
    }
    else if ((strstr(cmdtxt, "chisquaredmap") != NULL) && (strlen(cmdtxt) == 13)) {
      cmdnum = 53;
    }
    else if ((strstr(cmdtxt, "casadata") != NULL) && (strlen(cmdtxt) == 8)) {
      cmdnum = 54;
    }
    else if ((strstr(cmdtxt, "levels") != NULL) && (strlen(cmdtxt) == 6)) {
      cmdnum = 55;
    }
    else if ((strstr(cmdtxt, "printresults") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 56;
    }
    else if ((strstr(cmdtxt, "smoothmaps") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 57;
    }
    else if ((strstr(cmdtxt, "contours") != NULL) && (strlen(cmdtxt) == 8)) {
      cmdnum = 58;
    }
    else if ((strstr(cmdtxt, "contourlevels") != NULL) && (strlen(cmdtxt) == 13)) {
      cmdnum = 59;
    }
    else if ((strstr(cmdtxt, "contourcolour") != NULL) && (strlen(cmdtxt) == 13)) {
      cmdnum = 60;
    }
    else if ((strstr(cmdtxt, "contourcolor") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 60;
    }
    else if ((strstr(cmdtxt, "contourlinestyle") != NULL) && (strlen(cmdtxt) == 16)) {
      cmdnum = 61;
    }
    else if ((strstr(cmdtxt, "export") != NULL) && (strlen(cmdtxt) == 6)) {
      cmdnum = 62;
    }
    else if ((strstr(cmdtxt, "imageloc") != NULL) && (strlen(cmdtxt) == 8)) {
      cmdnum = 64;
    }
    else if ((strstr(cmdtxt, "imagetype") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 65;
    }
    else if ((strstr(cmdtxt, "firstcontour") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 66;
    }
    else if ((strstr(cmdtxt, "contourlogs") != NULL) && (strlen(cmdtxt) == 11)) {
      cmdnum = 67;
    }
    else if ((strstr(cmdtxt, "fluxlogs") != NULL) && (strlen(cmdtxt) == 8)) {
      cmdnum = 68;
    }
    else if ((strstr(cmdtxt, "fluxcut") != NULL) && (strlen(cmdtxt) == 7)) {
      cmdnum = 69;
    }
    else if ((strstr(cmdtxt, "usercut") != NULL) && (strlen(cmdtxt) == 7)) {
      cmdnum = 70;
    }
    else if ((strstr(cmdtxt, "wedgemultiplier") != NULL) && (strlen(cmdtxt) == 15)) {
      cmdnum = 71;
    }
    else if ((strstr(cmdtxt, "bfield") != NULL) && (strlen(cmdtxt) == 6)) {
      cmdnum = 72;
    }
    else if ((strstr(cmdtxt, "injectionindex") != NULL) && (strlen(cmdtxt) == 14)) {
      cmdnum = 73;
    }
    else if ((strstr(cmdtxt, "colourcolour") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 74;
    }
    else if ((strstr(cmdtxt, "colorcolor") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 74;
    }
    else if ((strstr(cmdtxt, "modelpower") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 75;
    }
    else if ((strstr(cmdtxt, "alwayszeroage") != NULL) && (strlen(cmdtxt) == 13)) {
      cmdnum = 76;
    }
    else if ((strstr(cmdtxt, "plotkpmodel") != NULL) && (strlen(cmdtxt) == 11)) {
      cmdnum = 77;
    }
    else if ((strstr(cmdtxt, "fitkpmodel") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 78;
    }
    else if ((strstr(cmdtxt, "paramlabels") != NULL) && (strlen(cmdtxt) == 11)) {
      cmdnum = 79;
    }
    else if ((strstr(cmdtxt, "printreject") != NULL) && (strlen(cmdtxt) == 11)) {
      cmdnum = 80;
    }
    else if ((strstr(cmdtxt, "gmin") != NULL) && (strlen(cmdtxt) == 4)) {
      cmdnum = 81;
    }
    else if ((strstr(cmdtxt, "gmax") != NULL) && (strlen(cmdtxt) == 4)) {
      cmdnum = 82;
    }
    else if ((strstr(cmdtxt, "posmap") != NULL) && (strlen(cmdtxt) == 6)) {
      cmdnum = 83;
    }
    else if ((strstr(cmdtxt, "scaletype") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 84;
    }
    else if ((strstr(cmdtxt, "beampos") != NULL) && (strlen(cmdtxt) == 7)) {
      cmdnum = 85;
    }
    else if ((strstr(cmdtxt, "beam") != NULL) && (strlen(cmdtxt) == 4)) {
      cmdnum = 86;
    }
    else if ((strstr(cmdtxt, "exportdata") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 87;
    }
    else if ((strstr(cmdtxt, "dataloc") != NULL) && (strlen(cmdtxt) == 7)) {
      cmdnum = 88;
    }
    else if ((strstr(cmdtxt, "usereduced") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 89;
    }
    else if ((strstr(cmdtxt, "conflevels") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 90;
    }
    else if ((strstr(cmdtxt, "modelmyears") != NULL) && (strlen(cmdtxt) == 11)) {
      cmdnum = 91;
    }
    else if ((strstr(cmdtxt, "setsingleregion") != NULL) && (strlen(cmdtxt) == 15)) {
      cmdnum = 92;
    }
    else if ((strstr(cmdtxt, "findinject") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 93;
    }
    else if ((strstr(cmdtxt, "mininject") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 94;
    }
    else if ((strstr(cmdtxt, "maxinject") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 95;
    }
    else if ((strstr(cmdtxt, "injectintervals") != NULL) && (strlen(cmdtxt) == 15)) {
      cmdnum = 96;
    }
    else if ((strstr(cmdtxt, "viewerrors") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 97;
    }
    else if ((strstr(cmdtxt, "telescope") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 98;
    }
    else if ((strstr(cmdtxt, "fitjptribble") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 99;
    }
    /*else if ((strstr(cmdtxt, "fitkptribble") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 100;
      }*/
    else if ((strstr(cmdtxt, "plotjptribble") != NULL) && (strlen(cmdtxt) == 13)) {
      cmdnum = 101;
    }
    else if ((strstr(cmdtxt, "plotkptribble") != NULL) && (strlen(cmdtxt) == 13)) {
      cmdnum = 102;
    }
    else if ((strstr(cmdtxt, "combineimages") != NULL) && (strlen(cmdtxt) == 13)) {
      cmdnum = 103;
    }
    else if ((strstr(cmdtxt, "mapfrommodel") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 104;
    }
    else if ((strstr(cmdtxt, "diffmap") != NULL) && (strlen(cmdtxt) == 7)) {
      cmdnum = 105;
    }
    else if ((strstr(cmdtxt, "scaleflux") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 106;
    }
    else if ((strstr(cmdtxt, "shiftimage") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 107;
    }
    else if ((strstr(cmdtxt, "fitintegrated") != NULL) && (strlen(cmdtxt) == 13)) {
      cmdnum = 108;
    }
    else if ((strstr(cmdtxt, "errormap") != NULL) && (strlen(cmdtxt) == 8)) {
      cmdnum = 109;
    }
    else if ((strstr(cmdtxt, "fullexport") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 110;
    }
    else if ((strstr(cmdtxt, "exportcompression") != NULL) && (strlen(cmdtxt) == 17)) {
      cmdnum = 111;
    }
    else if ((strstr(cmdtxt, "fullimport") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 112;
    }
    else if ((strstr(cmdtxt, "specindexcalctype") != NULL) && (strlen(cmdtxt) == 17)) {
      cmdnum = 113;
    }
    else if ((strstr(cmdtxt, "plotspecindex") != NULL) && (strlen(cmdtxt) == 13)) {
      cmdnum = 114;
    }
    else if ((strstr(cmdtxt, "specchisquared") != NULL) && (strlen(cmdtxt) == 14)) {
      cmdnum = 115;
    }
    else if ((strstr(cmdtxt, "resizeimage") != NULL) && (strlen(cmdtxt) == 11)) {
      cmdnum = 116;
    }
    else if ((strstr(cmdtxt, "resizeimages") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 116;
    }
    else if ((strstr(cmdtxt, "plotcioff") != NULL) && (strlen(cmdtxt) == 9)) {
      model = 5;
      cmdnum = 117;
    }
    else if ((strstr(cmdtxt, "plotcimodel") != NULL) && (strlen(cmdtxt) == 11)) {
      model = 6;
      cmdnum = 117;
    }
    /*else if ((strstr(cmdtxt, "plotciana") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 118;
      }*/
    else if ((strstr(cmdtxt, "fitcioff") != NULL) && (strlen(cmdtxt) == 8)) {
      model = 5;
      cmdnum = 119;
    }
    else if ((strstr(cmdtxt, "fitcimodel") != NULL) && (strlen(cmdtxt) == 10)) {
      model = 6;
      cmdnum = 119;
    }
    else if ((strstr(cmdtxt, "minmodelmyears") != NULL) && (strlen(cmdtxt) == 14)) {
      cmdnum = 120;
    }
    else if ((strstr(cmdtxt, "modelminoff") != NULL) && (strlen(cmdtxt) == 11)) {
      cmdnum = 121;
    }
    else if ((strstr(cmdtxt, "minmodeloff") != NULL) && (strlen(cmdtxt) == 11)) {
      cmdnum = 121;
    }
    else if ((strstr(cmdtxt, "modelmaxoff") != NULL) && (strlen(cmdtxt) == 11)) {
      cmdnum = 122;
    }
    else if ((strstr(cmdtxt, "maxmodeloff") != NULL) && (strlen(cmdtxt) == 11)) {
      cmdnum = 122;
    }
    else if ((strstr(cmdtxt, "varyoffage") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 123;
    }
    else if ((strstr(cmdtxt, "modelredshift") != NULL) && (strlen(cmdtxt) == 13)) {
      cmdnum = 124;
    }
    else if ((strstr(cmdtxt, "minoff") != NULL) && (strlen(cmdtxt) == 6)) {
      cmdnum = 125;
    }
    else if ((strstr(cmdtxt, "maxoff") != NULL) && (strlen(cmdtxt) == 6)) {
      cmdnum = 126;
    }
    else if ((strstr(cmdtxt, "minmyears") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 127;
    }
    else if ((strstr(cmdtxt, "jpdata") != NULL) && (strlen(cmdtxt) == 6)) {
      model = 1;
      cmdnum = 128;
    }
    else if ((strstr(cmdtxt, "kpdata") != NULL) && (strlen(cmdtxt) == 6)) {
      model = 2;
      cmdnum = 128;
    }
    else if ((strstr(cmdtxt, "tribbledata") != NULL) && (strlen(cmdtxt) == 11)) {
      model = 3;
      cmdnum = 128;
    }
    else if ((strstr(cmdtxt, "cioffdata") != NULL) && (strlen(cmdtxt) == 9)) {
      model = 5;
      cmdnum = 128;
    }
    else if ((strstr(cmdtxt, "cidata") != NULL) && (strlen(cmdtxt) == 6)) {
      model = 6;
      cmdnum = 128;
    }
    else if ((strstr(cmdtxt, "exportasfits") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 129;
    }
    else if ((strstr(cmdtxt, "dataintervals") != NULL) && (strlen(cmdtxt) == 13)) {
      cmdnum = 130;
    }
    else if ((strstr(cmdtxt, "targetname") != NULL) && (strlen(cmdtxt) == 10)) {
      cmdnum = 131;
    }
    else if ((strstr(cmdtxt, "rmsnoise") != NULL) && (strlen(cmdtxt) == 8)) {
      cmdnum = 132;
    }
    else if ((strstr(cmdtxt, "suppressconf") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 133;
    }
    else if ((strstr(cmdtxt, "injectionmap") != NULL) && (strlen(cmdtxt) == 12)) {
      cmdnum = 134;
    }
    else if ((strstr(cmdtxt, "injectionchimap") != NULL) && (strlen(cmdtxt) == 15)) {
      cmdnum = 135;
    }
    else if ((strstr(cmdtxt, "extendmodel") != NULL) && (strlen(cmdtxt) == 11)) {
      cmdnum = 136;
    }
    else if ((strstr(cmdtxt, "largetext") != NULL) && (strlen(cmdtxt) == 9)) {
      cmdnum = 137;
    }
    else if ( ((strstr(cmdtxt, "specindexerrors") != NULL) && (strlen(cmdtxt) == 15)) || ((strstr(cmdtxt, "specindexerror") != NULL) && (strlen(cmdtxt) == 14)) ){
      cmdnum = 138;
    }
    else if ( ((strstr(cmdtxt, "colourscheme") != NULL) && (strlen(cmdtxt) == 12)) || ((strstr(cmdtxt, "colorscheme") != NULL) && (strlen(cmdtxt) == 11)) ){
      cmdnum = 139;
    }
    else if ( ((strstr(cmdtxt, "invertcolours") != NULL) && (strlen(cmdtxt) == 13)) || ((strstr(cmdtxt, "invertcolors") != NULL) && (strlen(cmdtxt) == 12)) ){
      cmdnum = 140;
    }
    else if ((strstr(cmdtxt, "forceerrortype") != NULL) && (strlen(cmdtxt) == 14)) {
      cmdnum = 141;
    }
    else {
      cmdnum = 999;
    }


    switch (cmdnum) {
      
    case 0:
      ;
      break;

    case 1:
      help();
      break;

    case 2:
      mainloop = 0;
      break;

    case 3: // Load a data set

      if (arraysallocated == 0) {

	// Can these be reduced to MAXDATASETS and MAXNUMMAPS?
	
	xdim = (int *)calloc(MAXDATASETS,sizeof(int));
	ydim = (int *)calloc(MAXDATASETS,sizeof(int));
	bmaj = (float *)calloc(MAXDATASETS,sizeof(float));
	bmin = (float *)calloc(MAXDATASETS,sizeof(float));
	bpa = (float *)calloc(MAXDATASETS,sizeof(float));
	beamarea = (float *)calloc(MAXDATASETS,sizeof(float));
	imgnum = (int *)calloc(MAXDATASETS,sizeof(int));
	maxflux = (float **)calloc(MAXDATASETS,sizeof(float *));
	minflux = (float **)calloc(MAXDATASETS,sizeof(float *));
	
	rms = (float **)calloc(MAXDATASETS,sizeof(float *));
	flux = (float ****)calloc(MAXDATASETS,sizeof(float ***));
	fluxorder =  (int **)calloc(MAXDATASETS,sizeof(int *));
	bgbuff = (float ***)calloc(MAXDATASETS, sizeof(float **));
	
	frequency = (float **)calloc(MAXDATASETS,sizeof(float *));
	frequencyobs = (float **)calloc(MAXDATASETS,sizeof(float *));
	ra = (float **)calloc(MAXDATASETS,sizeof(float *));
	dec = (float **)calloc(MAXDATASETS,sizeof(float *));
	delt1 = (float **)calloc(MAXDATASETS,sizeof(float *));
	delt2 = (float **)calloc(MAXDATASETS,sizeof(float *));
	crpix1 = (float **)calloc(MAXDATASETS,sizeof(float *));
	crpix2 = (float **)calloc(MAXDATASETS,sizeof(float *));
	crota1 = (float **)calloc(MAXDATASETS,sizeof(float *));
        crota2 = (float **)calloc(MAXDATASETS,sizeof(float *));
	eq = (float **)calloc(MAXDATASETS,sizeof(float *));
	cellsize = (float *)calloc(MAXDATASETS,sizeof(float));
	
	redshift = (float *)calloc(MAXDATASETS,sizeof(float));
      
	//imgbuff = (float ***)calloc((MAXDATASETS * MAXNUMMAPS),sizeof(float **));
	//bgbuff = (float ***)calloc((MAXDATASETS * MAXNUMMAPS),sizeof(float **));
      
	arraysallocated = 1;

	currentset = numdatasets;
      }

      else if (numdatasets >= MAXDATASETS) {
	fprintf(stderr,"Error: Maximum number of datasets has been reached. Please delete or replace one!\n");
	break;
      }

      // *** This needs sorting to dynamically allocate the memory!!! ***

      /*else if (arraysallocated > 0) {

      // Else expand the arrays
      // numdatasets = total number of active datasets; currentset = set being modified; arraysallocated = how many arrays have been allocated memory

      currentset = numdatasets;

      int *pnt_arraysallocated;

      pnt_arraysallocated=&arraysallocated;

      allocate(xdim, ydim, bmaj, bmin, beamarea, imgnum, frequency, numdatasets, currentset, pnt_arraysallocated);
      //printf("*** Arrays allocated: %d ***\n", arraysallocated);
      }*/

      currentset = numdatasets; // This allows later addition or modification and deleting arrays
    
      if (load(dirname, currentset, numdatasets, imgnum, xdim, ydim, bmaj, bmin, bpa, beamarea, frequency, frequencyobs, regname, bgname, flux, rms, sigmalevel, fluxorder, zoom, border, titles, labels, axismarks, tmp_maptarget, maxflux, minflux, casadata, fluxlogs, fluxcut, usercut, wedgemultiplier, ra, dec, delt1, delt2, crpix1, crpix2, crota1, crota2, eq, cellsize, scaletype, posmap, beamposx, beamposy, plotbeam, bgbuff, xshift, yshift, redshift, largetxt, frac_border) == 0) {

	// Assign memory dynamically for settarget
	settarget[currentset] = (char *)malloc(sizeof(tmp_maptarget)+1);

	// This avoids the compications of sending over a 2 dimensional array of strings
	strcpy(setname[currentset], dirname);
	strcpy(setreg[currentset], regname);
	strcpy(setbg[currentset], bgname);
	strcpy(settarget[currentset], tmp_maptarget);

	// Set the flux cal errors
	fluxerrors(fluxcalerror[currentset], frequency[currentset], imgnum[currentset], instrument, flaterror);

	// Update the number of datasets (note this is done AFTER the data is loaded!)
	numdatasets++;
      }

      else {
	fprintf(stderr,"Error: Data failed to load!\n");
      }

      break;


    case 4:
      //For ls
      syscom(cmdtxt);
      break;


    case 5:
      //For generic commands
      cmdbuffer = readline("Enter system command (esc to exit): ");
      if (cmdbuffer[0] != 0) {
	add_history(cmdbuffer);
      }

      strcpy(cmdtxt, cmdbuffer);

      chomp(cmdtxt);

      if ((strstr(cmdtxt, "esc") != NULL) && (strlen(cmdtxt) == 3)) {
	printf("Escaping command...\n");
	break;
      }

      syscom(cmdtxt);
      break;


    case 6: // List all loaded datasets

      // This is easier than sending over string arrays, although make list() rather redundant
      printf("========================================================================\n\n");
      printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");
      // Check at least 1 dataset is loaded
      if (numdatasets > 0) {

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);
	}
      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      printf("\n========================================================================\n");

      break;


    case 7: // List properties

      // Check at least 1 dataset is loaded
      if (numdatasets > 0) {

	// This is easier than sending over string arrays, although make list() rather redundant
	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");
	// Check at least 1 dataset is loaded
	if (numdatasets > 0) {

	  // Loop through each of the sets
	  for (i=0;i<numdatasets;i++) {

	    list(i, imgnum[i], setname[i], setreg[i], setbg[i]);
	  }
	}
	else {
	  fprintf(stderr," Error: No datasets have yet been loaded!\n");
	}

	printf("\n========================================================================\n");


	validentry = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("\nEnter dataset number (-1 to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  setnumber = strtol(cmdbuffer, &endptr, 10);

	  if ( (setnumber == 0) && (cmdbuffer == endptr) ) {

	    printf("\nInvalid input, please try again...\n");
	    continue;
	    //printf("Invalid input, escaping command to protect the existing data...\n");
	    // break;
	  }

	  if (setnumber < 0) {
	    printf("Escaping command...\n");
	    break;
	  }

	  else if (setnumber >= numdatasets) {
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again.\n\n");
	    continue;
	  }

	  else {
	    printf("\n");
	    i = setnumber; // Keeps the command a bit cleaner
	    listprops(i, imgnum[i], setname[i], setreg[i], setbg[i], xdim[i], ydim[i], bmaj[i], bmin[i], beamarea[i]);
	    validentry = 1;
	  }
	}
      }

      else {
	fprintf(stderr,"\nError: No datasets have yet been loaded!\n\n");
      }

      break;


    case 8: // Adaptive regions
      
      if (numdatasets <= 0) {
	fprintf(stderr,"\nError: No datasets have yet been loaded!\n\n");
      }

      else {

	if (regmemallocated == 0 ) {

	  regflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  regionarray = (int ***)calloc(MAXDATASETS, sizeof(int *));
	  regnumber = (int *)calloc(MAXDATASETS, sizeof(int));
	  fluxerror = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  regionlocx = (float **)calloc(MAXDATASETS, sizeof(float *));
	  regionlocy = (float **)calloc(MAXDATASETS, sizeof(float *));
	  setaveraged = (int *)calloc(MAXDATASETS, sizeof(int));
	  regionsize = (int **)calloc(MAXDATASETS, sizeof(int *));
	  regmaxflux = (float **)calloc(MAXDATASETS,sizeof(float *));
	  regminflux = (float **)calloc(MAXDATASETS,sizeof(float *));

	  regmemallocated = 1;
	}


	// This is easier than sending over string arrays, although make list() rather redundant
	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");
	// Check at least 1 dataset is loaded

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}
  
	printf("\n========================================================================\n");


	validentry = 0;

	while (validentry == 0) {

	// printf("Enter the index number of the data set to create regions for (-1 to escape, 888 for all): ");

	tmp_currentset = currentset;

	cmdbuffer = readline("Enter the index number of the data set to create regions for (-1 to escape, 888 for all): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	currentset = strtol(cmdbuffer, &endptr, 10);

	//if (scanf("%d", &currentset) == 0) {

	if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (currentset < 0) {
	  printf("Escaping command...\n");
	  currentset = tmp_currentset;
	  break;
	}
	else if ( (currentset >= numdatasets) && (currentset != 888) ) {
	  fprintf(stderr,"\nError: This dataset does not exist! Please try again.\n\n");
	  continue;
	}

	else {
	  validentry = 1;
	  if (currentset == 888) {
	    for (i=0; i<numdatasets; i++){

	      adaptiveregions(imgnum, xdim, ydim, flux, signaltonoise, rms, searcharea, regionarray, regflux, regnumber, fluxerror, averegions, regionlocx, regionlocy, fluxorder, i, maptomap, hotpixels, sigmalevel, titles, labels, axismarks, regionsize, regmaxflux, regminflux, fluxcalerror[i], onsourcenoise, regionsset[i], border[i], zoom, ra[i][posmap], dec[i][posmap], cellsize[i], scaletype, bmaj[i], bmin[i], bpa[i], beamposx, beamposy, plotbeam, xshift, yshift, smoothmaps, largetxt);

	      if(regnumber[i] > 0) { // If we have valid regions, mark as set
		regionsset[i] = 1;
		setaveraged[i] = averegions;
	      }
	      else {
		regionsset[i] = 0;
	      }

	    }
	    currentset = 0; // Put the current set back to 0 to save any mix ups later
	  }

	  else {
	    adaptiveregions(imgnum, xdim, ydim, flux, signaltonoise, rms, searcharea, regionarray, regflux, regnumber, fluxerror, averegions, regionlocx, regionlocy, fluxorder, currentset, maptomap, hotpixels, sigmalevel, titles, labels, axismarks, regionsize, regmaxflux, regminflux, fluxcalerror[currentset], onsourcenoise, regionsset[currentset], border[currentset], zoom, ra[currentset][posmap], dec[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, xshift, yshift, smoothmaps, largetxt);


	      if(regnumber[currentset] > 0) { // If we have valid regions, mark as set
		regionsset[currentset] = 1;
		setaveraged[currentset] = averegions;
	      }
	      else {
		regionsset[currentset] = 0;
	      }
	    

	  }
	}
      }
    }
    
      break;


    case 9: //Set the zoom level for maps

      validentry = 0;

      printf("Zoom currently set to %.2f (DEFAULT: %.2f)\n", zoom, DEFAULTZOOM);

      while (validentry == 0) {

	//Protect the current value incase we need to escape
	tmp_zoom = zoom;

	cmdbuffer = readline("Please enter a value to zoom (-1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	zoom = strtod(cmdbuffer, &endptr);

	if ( (zoom == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else if ( (zoom < 1.0) && (zoom >= -0.001) ) {
	  printf("Value must be greater than or equal to 1. Please try again.\n");
	  continue;
	}

	else if (zoom < 0.001) {
	  zoom = tmp_zoom;
	  printf("Escaping command...\n");
	  break;
	}

	else {
	  // Change the absolute zoom values to new fractional zoom value
	  for (i=0;i<numdatasets;i++) {
	    zoomedx =  xdim[i] / zoom;
	    istart = xshift + ( (xdim[i] - (int)zoomedx) / 2);
	    iend = xshift + xdim[i] - ( (xdim[i] - (int)zoomedx) / 2);
	    border[i] = ceil((iend-istart) * frac_border); // Always round up
	  }

	  printf("Zoom level for maps is now set to %.2f\n", zoom);
	  validentry = 1;
	}
      }

      break;


    case 10: //Set the border level for maps
      //Protect the current value incase we need to escape
      tmp_fracborder = frac_border;

      printf("Border currently set to %.2f per cent (DEFAULT: %d)\n", frac_border*100, DEFAULTBORDER);

      validentry = 0;

      while (validentry == 0) {

	cmdbuffer = readline("Please enter a border as a percentage of the axis (-1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	frac_border = strtod(cmdbuffer, &endptr) / 100;

	if ( ((frac_border < 0.001) && (frac_border > -1e-23)) && (cmdbuffer == endptr) ) {	
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else if (frac_border < -1e-23) {
	  frac_border = tmp_fracborder;
	  printf("Escaping command...\n");
	  break;
	}
	else {

	  if (frac_border < 1e-23) {
	    frac_border = 1e-23; // Ensure we get at least a 1 pixel border to prvent an invalid output. Note <0 already caught
	  }
	
	  for (i=0;i<numdatasets;i++) {
	    zoomedx =  xdim[i] / zoom;
	    istart = xshift + ( (xdim[i] - (int)zoomedx) / 2);
	    iend = xshift + xdim[i] - ( (xdim[i] - (int)zoomedx) / 2);
	    border[i] = ceil((iend-istart) * frac_border); // Always round up
	  }

	  printf("Border for maps is now set to %.2f per cent\n", (frac_border*100));
	  validentry = 1;
	}
      }

      break;


    case 11: //Set the sigma value


      //Protect the current value incase we need to escape
      tmp_sigma = sigmalevel;

      validentry = 0;

      printf("Sigma level currently set to %.2f (DEFAULT: %.2f)\n", sigmalevel, DEFAULTSIGMA);

      while (validentry == 0) {

	cmdbuffer = readline("Please enter a new value for sigma (-1 to escape, integers only): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	sigmalevel = strtod(cmdbuffer, &endptr);

	if ( (sigmalevel == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (sigmalevel < 0) {
	  sigmalevel = tmp_sigma;
	  printf("Escaping command...\n");
	  break;
	}

	else {
	  printf("Sigma level for detection is now set to %.2f\n", sigmalevel);
	  validentry = 1;
	}
      }

      break;


    case 12: //Turn titles on and off for plots and maps

      if (titles == 1) {
	titles = 0;
	printf("Plot and map titles are now turned off!\n");
      }
      else if (titles == 0) {
	titles = 1;
	printf("Plot and map titles are now turned on!\n");
      }
      else {
	printf("\nThe value of titles is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", titles);
      }

      break;


    case 13: //Turn the labels on and off for plots and maps

      if (labels == 1) {
	labels = 0;
	printf("Plot and map labels are now turned off!\n");
      }
      else if (labels == 0) {
	labels = 1;
	printf("Plot and map labels are now turned on!\n");
      }
      else {
	printf("\nThe value of labels is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", labels);
      }

      break;


    case 14: //Turn the axis marking on and off for plots and maps

      if (axismarks == 1) {
	axismarks = 0;
	printf("Plot and map axis markings are now turned off!\n");
      }
      else if (axismarks == 0) {
	axismarks = 1;
	printf("Plot and map axis marking are now turned on!\n");
      }
      else {
	printf("\nThe value of axismarks is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", axismarks);
      }

      break;


    case 15: // Change the maptomap value

      //Protect the current value incase we need to escape
      tmp_maptomap = maptomap;

      printf("maptomap currently set to %.2f (DEFAULT: %.2f)\n", maptomap, DEFAULTMAPTOMAP);

      validentry = 0;

      while (validentry == 0) {

	cmdbuffer = readline("Please enter a value for map to map variations (999 to escape, -1 to turn off): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	maptomap = strtod(cmdbuffer, &endptr);

	if ( (maptomap == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (maptomap == 999) {
	  maptomap = tmp_maptomap;
	  printf("Escaping command...\n");
	  break;
	}
	else if (maptomap <= 0) {
	  printf("Map to map variation flagging is now turned off!\n");
	  break;
	}

	else {
	  printf("Maximum allowable map to map variation is now set to %.2f\n", maptomap);
	  validentry = 1;
	}
      }

      break;


    case 16: // Change the hotpixels value

      //Protect the current value incase we need to escape
      tmp_hotpixels = hotpixels;

      printf("hotpixels currently set to %.2f (DEFAULT: %.2f)\n", hotpixels, DEFAULTHOTPIXELS);

      validentry = 0;

      while (validentry == 0) {

	cmdbuffer = readline("Please enter a value for hotpixel variations (-1 to escape, 0 to turn off): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	hotpixels = strtod(cmdbuffer, &endptr);

	if ( (hotpixels == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (hotpixels < (-0.05)) {
	  hotpixels = tmp_hotpixels;
	  printf("Escaping command...\n");
	  break;
	}
	else if (hotpixels <= 0.01) {
	  printf("Hot and cold pixel error checking is now turned off\n");
	  break;
	}
	else {
	  printf("Hot and cold pixel value now set to %.2f\n", hotpixels);
	  validentry = 1;
	}
      }

      break;


    case 17: //Change the averaging of regions value

      if (averegions == 1) {
	averegions = 0;
	printf("Averaging of regions is now turned off!\n");
      }
      else if (averegions == 0) {
	averegions = 1;
	printf("Averaging of regions is now turned on!\n");
      }
      else {
	printf("\nThe value of averegions is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", averegions);
      }

      break;

    case 18: //Set the signal to noise limit

      //Protect the current value incase we need to escape
      tmp_signaltonoise = signaltonoise;

      printf("Signal to noise currently set to %d (DEFAULT: %d)\n", signaltonoise, DEFAULTSIGNOISE);


      validentry = 0;

      while (validentry == 0) {

	cmdbuffer = readline("Enter a new signal to noise limit (-1 to escape, integers only): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	signaltonoise = strtol(cmdbuffer, &endptr, 10);

	if ( (signaltonoise == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (signaltonoise < 0) {
	  signaltonoise = tmp_signaltonoise;
	  printf("Escaping command...\n");
	  break;
	}

	if (signaltonoise == 0) {
	  //signaltonoise = tmp_signaltonoise;
	  printf("\nThe signal to noise must be greater than 0, please try again...\n\n");
	  continue;
	}

	else {
	  printf("Signal to noise limit is now set to %d\n", signaltonoise);
	  validentry = 1;
	}
      }

      break;


    case 19: // Set the search area limit

      //Protect the current value incase we need to escape
      tmp_searcharea = searcharea;

      printf("Search area currently set to %d (DEFAULT: %d)\n", searcharea, DEFAULTSEARCHAREA);

      validentry = 0;

      while (validentry == 0) {

	cmdbuffer = readline("Enter a new maximum search area (-1 to escape, integers only): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	searcharea = strtol(cmdbuffer, &endptr, 10);

	if ( (searcharea == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (searcharea < 0) {
	  searcharea = tmp_searcharea;
	  printf("Escaping command...\n");
	  break;
	}
	else if (searcharea == 0) {
	  //searcharea = tmp_searcharea;
	  printf("\nThe maximum search area must be greater than 0. Please try again...\n\n");
	  continue;
	}
	else {
	  validentry = 1;
	  printf("Maximum search area is now set to %d\n", searcharea);
	}
      }

      break;


    case 20: // Calculate the spectral index of each region

      if (numdatasets > 0) {

	if (export == 1) {

	  // Open the directory that is currently set and check it exists
	  dir = opendir(imageloc);

	  if (dir != NULL) {
	    closedir (dir); // If we find it move on without comment
	  }
	  else {
	    printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	    break;
	  }
	
	}

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");
	// Check at least 1 dataset is loaded
	if (numdatasets > 0) {

	  // Loop through each of the sets
	  for (i=0;i<numdatasets;i++) {

	    list(i, imgnum[i], setname[i], setreg[i], setbg[i]);
	  }
	}
	else {
	  fprintf(stderr," Error: No datasets have yet been loaded!\n");
	}

	printf("\n========================================================================\n");


	tmp_currentset = currentset;

	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the number of the data set to determine the spectral index map (-1 to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    currentset = tmp_currentset;
	    escape = 1;
	    break;
	  }
	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    printf("\nThis value is greater than the number of data sets that exists! Please try again...\n\n");
	    continue;
	  }
	  else if ( (imgnum[currentset] <= 1) && (currentset != 888) ) { // Check the set has at least two maps
	    fprintf(stderr,"\nError: The data set must contain a minimum of two images for fitting. Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    //escape = 1;
	    continue;
	  }

	  else {
	    validentry = 1;
	  }
	}
      
	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if (contours == 1) {

	  if ( (exportfits != 1)  || (export != 1) ) {

	    validentry = 0;

	    printf("============================================\n\n");
	    printf(" Index  |  Frequency (Obs.)  |  Frequency (Rest)  \n\n");

	    // Loop through each of the sets
	    for (a=0; a<imgnum[currentset]; a++) {

	      printf("   %d       %.2e Hz             %.2e Hz\n", fluxorder[currentset][a], frequencyobs[currentset][fluxorder[currentset][a]], frequency[currentset][fluxorder[currentset][a]]);

	    }

	    printf("\n============================================\n");
    

	    while (validentry == 0) {

	      cmdbuffer = readline("Select a maps frequency for the contours (-1 to escape): ");
	      if (cmdbuffer[0] != 0) {
		add_history(cmdbuffer);
	      }

	      chomp(cmdbuffer);

	      currentmap = strtol(cmdbuffer, &endptr, 10);

	      if ( (currentmap == 0) && (cmdbuffer == endptr) ) {
		printf("\nInvalid input, please try again...\n\n");
		continue;
	      }

	      if (currentmap < 0) {
		printf("Escaping command...\n");
		escape = 1;
		break;
	      }
	
	      else if (currentmap >= imgnum[currentset]) { // Check the set exists
		printf("\nInvalid selection, please choose again...\n\n");
		continue;
	      }
	      else {
		validentry = 1;
	      }

	    }
	  }
	}
	
	else if ( (contours == 1) && (currentset == 888) ) { // Multimaps always uses map 0 as reference. This can easily be changed if needed
	  currentmap = 0;
	}
	else {
	  currentmap = 0;
	}
      
	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if ( (currentset == 888) && (export == 0) ) {
	  looplimit = numdatasets;
	  extwincontrol = 1;
	  sprintf(output,"/xs");
	  cpgopen(output);
	}
	else if ( (currentset == 888) && (export == 1) ) {
	  looplimit = numdatasets;
	  extwincontrol = 0;
	}
	else {
	  looplimit = 1;
	  extwincontrol = 0;
	}


	if (specallocated == 0) {

	  alpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	  minsi = (float *)calloc(MAXDATASETS, sizeof(float));
	  maxsi = (float *)calloc(MAXDATASETS, sizeof(float));
	  specindexintercept = (double **)calloc(MAXDATASETS, sizeof(double *));
	  specindexchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	  specindexerror = (float **)calloc(MAXDATASETS, sizeof(float *)); // This is the error on the actual index
	  specindex_modelflux = (double ***)calloc(MAXDATASETS, sizeof(double **));
	  specindex_modelerror = (double ***)calloc(MAXDATASETS, sizeof(double **)); // These are the error points for plotting
	  specindex_errortype = (int *)calloc(MAXDATASETS, sizeof(int)); // 0 for propagation, 1 for stddev wls
	  
	  specallocated = 1;
	}
      

	for (j=0; j<looplimit; j++){

	  if ( (looplimit > 1) || ( (looplimit == 1) && (currentset == 888) ) ) {
	    currentset = j;
	  }
	  
	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets
	  if (regionsset[currentset] != 1) {
	    printf("\n*** Warning: The regions for data set %d have not been set! Skipping this set... ***\n\n", j);
	    continue;	  
	  }
	  else if (imgnum[currentset] <= 1) { // Check the set exists
	    fprintf(stderr,"\n*** Warning: Data set %d does not contain the minimum two images required for fitting. Skipping this set... ***\n\n", j);
	    continue;
	  }

	  alpha[currentset] = (float *)calloc(xdim[currentset]*ydim[currentset], sizeof(float));
	  specindexintercept[currentset] = (double *)calloc(xdim[currentset]*ydim[currentset], sizeof(double));
	  specindexchisquared[currentset] = (float *)calloc(xdim[currentset]*ydim[currentset], sizeof(float));

	  specindexerror[currentset] = (float *)calloc(regnumber[currentset]+1, sizeof(float));
	  specindex_modelflux[currentset] = (double **)calloc(regnumber[currentset]+1, sizeof(double *));
	  specindex_modelerror[currentset] = (double **)calloc(regnumber[currentset]+1, sizeof(double *));
	  
	  for (j=1; j<=regnumber[currentset]; j++){
	    specindex_modelflux[currentset][j] = (double *)calloc(modelres+1, sizeof(double));
	    specindex_modelerror[currentset][j] = (double *)calloc(modelres+1, sizeof(double));
	  }

	  if (indexcalctype == 2) {
	    calcspecindexweightedgsl(regflux[currentset], regionarray[currentset], alpha[currentset], specindexintercept[currentset], specindexchisquared[currentset], specindex_modelflux[currentset], specindex_modelerror[currentset], fluxerror[currentset], regnumber[currentset], imgnum[currentset], frequency[currentset], printindex, &minsi[currentset], &maxsi[currentset], modelres, fluxorder[currentset], specindexerror[currentset], &specindex_errortype[currentset], forceerrortype);

	    specindex_modelres[currentset] = modelres;
	    specindextype[currentset] = 2;
	  }
	  else if (indexcalctype == 1) {
	    calcspecindex(regflux[currentset], regionarray[currentset], alpha[currentset], regnumber[currentset], imgnum[currentset], frequency[currentset], printindex, &minsi[currentset], &maxsi[currentset]);

	    specindextype[currentset] = 1;
	  }
	  else {
	    printf("\n*** Warning: Unable to determine indexcalctype. Reverting to non-weighted inbuilt method ***\n\n");
	    calcspecindex(regflux[currentset], regionarray[currentset], alpha[currentset], regnumber[currentset], imgnum[currentset], frequency[currentset], printindex, &minsi[currentset], &maxsi[currentset]);

	    specindextype[currentset] = 1;
	  }

	  specindexset[currentset] = 1;

	  if (export == 1) {

	    extwincontrol = 1; // We always want external control if exporting

	    time( &currenttime );
	    time_struct = localtime ( &currenttime );

	    strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	    if ( ( export == 1 ) && (exportfits != 1) ) {

	      if (contours == 1) {
		sprintf(output,"%s/%s_SpectralIndex_Contours_%s.%s/%s", imageloc, settarget[currentset], timebuff, imagetype, imagetype);
	      }
	      else {
		sprintf(output,"%s/%s_SpectralIndex_%s.%s/%s", imageloc, settarget[currentset], timebuff, imagetype, imagetype);
	      }
	      printf("\nExporting as %s\n", output);
	      cpgopen(output);
	    }
	    else if ( ( export == 1 ) && (exportfits == 1) ) {
	      sprintf(output,"%s/%s_SpectralIndex_%s.fits", imageloc, settarget[currentset], timebuff);
	      printf("\nExporting as %s\n", output);
	    }
	  }
	  else {
	    sprintf(output,"/xs");
	  }

	  sprintf(maptitle,"     Spectral Index Map of %s Between %.2e and %.2e Hz", settarget[currentset],frequencyobs[currentset][fluxorder[currentset][0]], frequencyobs[currentset][fluxorder[currentset][imgnum[currentset]-1]]);

	  plotspecindex(alpha, regnumber[currentset], xdim, ydim, settarget[currentset], minsi, maxsi, regionarray, currentset, maptitle, zoom, border[currentset], titles, labels, axismarks, symbol, smoothmaps, flux[currentset][currentmap], extwincontrol, contours, contourlevels, contourcolour, contourlinestyle, firstcontour, contourlogs, output, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, xshift, yshift, frequency[currentset][fluxorder[currentset][0]], frequency[currentset][fluxorder[currentset][imgnum[currentset]-1]], export, exportfits, largetxt);

	  if ( (export == 1) && (exportfits != 1) && (extwincontrol == 1) ){
	    cpgclos();
	    }
	  
	}

	if (extwincontrol == 1) {
	  extwincontrol = 0;
	}
      
	/*if ( (export != 1) && (extwincontrol == 1) ) {
	  printf("Here 2\n");
	  cpgclos();
	}*/
    
      }
      else {
	fprintf(stderr," Error: No data sets have yet been loaded!\n");
      }


      break;


    case 21: //Turn the printing of all spectral index values on and off

      if (printindex == 1) {
	printindex = 0;
	printf("Printing of spectral index values is now turned off!\n");
      }
      else if (printindex == 0) {
	printindex = 1;
	printf("Printing of spectral index values is now turned on!\n");
      }
      else {
	printf("\nThe value of printindex is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", printindex);
      }

      break;


    case 22: //Use a fixed regions from a given data set rather than each using its own

      if (export == 1) {

	// Open the directory that is currently set and check it exists
	dir = opendir(imageloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	  break;
	}
      }
      

      if (numdatasets > 1) { 

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");
	// Check at least 1 dataset is loaded
	if (numdatasets > 0) {

	  // Loop through each of the sets
	  for (i=0;i<numdatasets;i++) {

	    list(i, imgnum[i], setname[i], setreg[i], setbg[i]);
	  }
	}
	else {
	  fprintf(stderr," Error: No datasets have yet been loaded!\n");
	}

	printf("\n========================================================================\n");



	validentry = 0;
	escape =0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set of the regions to be used as reference (-1 to escape, integers only): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  fixedregref = strtol(cmdbuffer, &endptr, 10);

	  // Get the reference data set
	  if ( (fixedregref == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (fixedregref < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if (fixedregref >= numdatasets) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again.\n\n");
	    continue;
	  }

	  else if ( regionsset[fixedregref] != 1 ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for the reference data set! Please run setregions first.\nExiting command...\n\n");
	    escape = 1;
	    break;
	  }
	  else {
	    validentry = 1;
	  }
	}

	if (escape == 1) {
	  break;
	}

	// Get the data set to apply it to
	validentry = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the number the data set to apply the reference regions (-1 to escape, integers only): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  fixedregnum = strtol(cmdbuffer, &endptr, 10);

	  if ( (fixedregnum == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (fixedregnum < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if (fixedregnum >= numdatasets) { 
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	if (escape == 1) {
	  break;
	}

	// Now apply one to the other

	printf("Fixing the regions of data set %d to data set %d...\n", fixedregref, fixedregnum);

	if (applyfixedregions(imgnum, xdim, ydim, flux, rms, regionarray, regflux, regnumber, fluxerror, averegions, regionlocx, regionlocy, titles, labels, axismarks, fixedregnum, fixedregref, regionsset, fluxcalerror[fixedregnum], onsourcenoise, ra[fixedregnum][posmap], dec[fixedregnum][posmap], cellsize[fixedregnum], border[currentset], zoom, xshift, yshift, smoothmaps, scaletype, largetxt, export) == 0) {
	  
	  regionsset[fixedregnum] = 1; // Mark the set as having regions defined

	  printf("*** Warning: Fixed regions may cause unexpected results if the two datasets are not well matched or you are using the setsingleregion command ***\n\n");
	}

      }

      else {
	fprintf(stderr," Error: At least 2 datasets must be loaded to use this command!\n");
      }

      break;


    case 23: // Flux maps of a given data set

      if (numdatasets > 0) {

	if (export == 1) {

	  // Open the directory that is currently set and check it exists
	  dir = opendir(imageloc);

	  if (dir != NULL) {
	    closedir (dir); // If we find it move on without comment
	  }
	  else {
	    printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	    break;
	  }
	}

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to produce the flux maps (-1 to escape, integers only): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if (currentset >= numdatasets) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}
      
	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}
      
	// Raw data or adaptive regions?

	validentry = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Use raw data (0) or adaptive regions (1)? (-1 to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  dataselect = strtol(cmdbuffer, &endptr, 10);

	  if ( (dataselect == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (dataselect < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if (dataselect > 1) { // Check the set exists
	    fprintf(stderr,"\nError: Invalid selection, please try again.\n\n");
	    continue;
	  }

	  else if ( (dataselect == 1) && (regionsset[currentset] != 1) ) { // Check the regions have been set if mapping adaptive
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first or use the raw data option.\n\n");
	    continue;
	  }

	  else {
	    validentry = 1;
	  }
	}

	if (escape == 1) {
	  break;
	}

	// All maps or just one?
	validentry = 0;


	printf("============================================\n\n");
	printf(" Index  |  Frequency (Obs.)  |  Frequency (Rest)  \n\n");

	// Loop through each of the sets
	for (a=0; a<imgnum[currentset]; a++) {

	  printf("   %d       %.2e Hz             %.2e Hz\n", fluxorder[currentset][a], frequencyobs[currentset][fluxorder[currentset][a]], frequency[currentset][fluxorder[currentset][a]]);

	}

	printf("\n============================================\n");
    

   	while (validentry == 0) {

	  cmdbuffer = readline("Select a frequency to produce a flux map for (-1 to escape, 888 all maps): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  mapselect = strtol(cmdbuffer, &endptr, 10);

	  if ( (mapselect == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (mapselect < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }	
	  else if ( (mapselect >= imgnum[currentset]) && (mapselect != 888) ) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again.\n\n");
	    continue;
	  }
	  else if (mapselect == 888) { // Confirm the use of all maps

	    printf("\nWarning: You cannot escape the production of flux maps until the last one has been viewed. If you have a large maps either in number or size this may take some time!\n\n");

	    validentry2 = 0;

	    while (validentry2 == 0) {

	      escape2 = 0;

	      cmdbuffer = readline("Are you sure you want to continue using all maps? (No [0], Yes [1], -1 to escape): ");
	      if (cmdbuffer[0] != 0) {
		add_history(cmdbuffer);
	      }

	      chomp(cmdbuffer);

	     checkuseall = strtol(cmdbuffer, &endptr, 10);

	      if ( (checkuseall == 0) && (cmdbuffer == endptr) ) {
		printf("\nInvalid input, please try again...\n\n");
		continue;
	      }

	      if (checkuseall < 0) {
		printf("Escaping command...\n");
		escape2 = 1;
		escape = 1;
		break;
	      }
	      else if (checkuseall > 1) {
		printf("Invalid selection, please choose again...\n");
		continue;
	      }

	      else if (checkuseall == 0) {
		mapselect = 0;
		validentry2 = 1;
	      }

	      else if (checkuseall == 1){
		validentry = 1;
		validentry2 = 1;
	      }

	      else {
		printf("Invalid entry, please try again...\n");
		continue;
	      }
	    }
	    if (escape2 == 1) {
	      break;
	    }

	  }	
	  else {
	    validentry = 1;
	  }
	}
      
	if (escape == 1) {
	  break;
	}


	if (dataselect == 0) {

	  if (mapselect == 888) {

	    // Set window control to external and open one up. This gives a stop after each map.
	    extwincontrol = 1;

	    // If we are not exporting, open one window and scroll through
	    if ( export == 0 ) {

	      sprintf(output,"/xs");
	      cpgopen(output);
	    }


	    for(a=0; a<imgnum[currentset]; a++) {

	      // If exporting, output open a new session for each with a different name
	      if ( export == 1 ) {

		time( &currenttime );
		time_struct = localtime ( &currenttime );

		strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

		sprintf(output,"%s/%s_%.2eHz_%s.%s/%s", imageloc, settarget[currentset], frequencyobs[currentset][a], timebuff, imagetype, imagetype);

		printf("\nExporting as %s\n", output);

		cpgopen(output);
	      }

	      
	      sprintf(maptitle,"Flux of %s at %.2e Hz (%.2f Sigma)", settarget[currentset], frequencyobs[currentset][a], sigmalevel);

	      fluxmap(flux, rms, xdim, ydim, a, maxflux, minflux, currentset, border[currentset], zoom, sigmalevel, maptitle, titles, labels, axismarks, extwincontrol, smoothmaps, fluxlogs, fluxcut, usercut, wedgemultiplier, output, ra[currentset][posmap], dec[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, xshift, yshift, largetxt, export);

	      // Close for each loop if we are exporting
	      if ( export == 1 ) {
		cpgclos();
	      }

	    }
	    extwincontrol = 0; //Set back to 0 to be safe
	  }

	  else { // Do the single selected map

	    extwincontrol = 1;

	    // This is fine as we are only ever opening or outputting one image
	    if ( export == 1 ) {

	      time( &currenttime );
	      time_struct = localtime ( &currenttime );

	      strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	      sprintf(output,"%s/%s_%.2eHz_%s.%s/%s", imageloc, settarget[currentset], frequencyobs[currentset][mapselect], timebuff, imagetype, imagetype);

	      printf("\nExporting as %s\n", output);

	    }
	    else {
	      sprintf(output,"/xs");
	    }
	    
	    cpgopen(output);

	    sprintf(maptitle,"Flux of %s at %.2e Hz (%.2f Sigma)", settarget[currentset], frequencyobs[currentset][mapselect], sigmalevel);

	    fluxmap(flux, rms, xdim, ydim, mapselect, maxflux, minflux, currentset, border[currentset], zoom, sigmalevel, maptitle, titles, labels, axismarks, extwincontrol, smoothmaps, fluxlogs, fluxcut, usercut, wedgemultiplier, output, ra[currentset][posmap], dec[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, xshift, yshift, largetxt, export);

	    extwincontrol = 0;

	    // Close off the output
	    cpgclos();

	  }
	}

	else if (dataselect == 1) {

	  // Check there is more than one region
	  if (regnumber[currentset] <= 1) {
	    printf("\n Error: Dataset must have more than 1 region to map the adaptive regions\n\n");
	    break;
	  }

	  // Average the flux over the regions?
	  if (setaveraged[currentset] == 1) {
	    printf("Data set is already averaged by region, skipping option.");
	  }
	  else {

	    validentry = 0;

	    while (validentry == 0) {

	      cmdbuffer = readline("Average the flux by region size? (No [0], Yes [1], -1 to escape): ");
	      if (cmdbuffer[0] != 0) {
		add_history(cmdbuffer);
	      }

	      chomp(cmdbuffer);

	      doregaveraging = strtol(cmdbuffer, &endptr, 10);

	      if ( (doregaveraging == 0) && (cmdbuffer == endptr) ) {
		printf("\nInvalid input, please try again...\n\n");
		continue;
	      }

	      if (doregaveraging < 0) {
		printf("Escaping command...\n");
		escape = 1;
		break;
	      }
	      else if (doregaveraging > 1) {
		printf("\nInvalid selection, please try again...\n\n");
		continue;
	      }
	      else {
		validentry = 1;
	      }
	    }

	    if (escape == 1) {
	      break;
	    }
	  }


	  if (mapselect == 888) {

	    extwincontrol = 1;

	    if ( export == 0 ) {
	      sprintf(output,"/xs");
	      cpgopen(output);
	    }

	    for(a=0; a<imgnum[currentset]; a++) {

	      // If exporting, output open a new session for each with differing names
	      if ( export == 1 ) {

		time( &currenttime );
		time_struct = localtime ( &currenttime );

		strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

		sprintf(output,"%s/%s_%.2eHz_%s.%s/%s", imageloc, settarget[currentset], frequencyobs[currentset][a], timebuff, imagetype, imagetype);

		printf("\nExporting as %s\n", output);

		cpgopen(output);
	      }

	      sprintf(maptitle,"Flux of %s at %.2e Hz", settarget[currentset], frequencyobs[currentset][a]);

	      regionfluxmap(regflux, xdim, ydim, a, regmaxflux, regminflux, currentset, border[currentset], zoom, maptitle, titles, labels, axismarks, regionarray, extwincontrol, regionsize, regnumber, doregaveraging, smoothmaps, output, ra[currentset][posmap], dec[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, xshift, yshift, largetxt, export);

	      // Close for each loop if we are exporting
	      if ( export == 1 ) {
		cpgclos();
	      }
	    }
	    extwincontrol = 0;
	  }
	
	  else {

	    extwincontrol = 1;

	    if ( export == 1 ) {

	      time( &currenttime );
	      time_struct = localtime ( &currenttime );

	      strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	      sprintf(output,"%s/%s_%.2eHz_%s.%s/%s", imageloc, settarget[currentset], frequencyobs[currentset][mapselect], timebuff, imagetype, imagetype);

	    }
	    else {
	      sprintf(output,"/xs");
	      cpgopen(output);
	    }

	    cpgopen(output);

	    sprintf(maptitle,"Flux of %s at %.2e Hz", settarget[currentset], frequencyobs[currentset][mapselect]);

	    regionfluxmap(regflux, xdim, ydim, mapselect, regmaxflux, regminflux, currentset, border[currentset], zoom, maptitle, titles, labels, axismarks, regionarray, extwincontrol, regionsize, regnumber, doregaveraging, smoothmaps, output, ra[currentset][posmap], dec[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, xshift, yshift, largetxt, export);

	    extwincontrol = 0;

	    cpgclos();
	  }
	}
      }

      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      /*if (export == 0) {
	cpgclos();
	}*/
      
      
      break;


    case 24: // Plot the flux of a data set

 
      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to plot the flux (-1 to escape, integers only): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if (currentset >= numdatasets) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}
      
	// Check the regions have been set
	if (regionsset[currentset] != 1) {
	  fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	  continue;
	}

	if ( (autoscalex == 0) && (autoscaley == 0) ) {

	  printf("Auto scaling of both axes is currently set to off,  check the values are correct for this plot...\n");
	}

	else if ( (autoscalex == 1) && (autoscaley == 0) ) {

	  printf("Auto scaling of the Y axis is currently set to off (X axis is on), check the values are correct for this plot...\n");
	}

	else if ( (autoscalex == 0) && (autoscaley == 1) ) {

	  printf("Auto scaling of the X axis is currently set to off (Y axis is on), the values are correct for this plot...\n");
	}

	else {
	  printf("Auto scaling of axis is currently set to on for both axes\n");
	}

	// If all seems fine, do the plot

	if ( export == 1 ) {

	  time( &currenttime );
	  time_struct = localtime ( &currenttime );

	  strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	  sprintf(output,"%s/%s_FluxByRegion_%s.%s/%s", imageloc, settarget[currentset], timebuff, imagetype, imagetype);

	  printf("\nExporting as %s\n", output);

	}
	else {
	  sprintf(output,"/xs");
	}

	sprintf(maptitle,"%s between %.2e Hz and %.2e Hz", settarget[currentset], frequencyobs[currentset][fluxorder[currentset][0]], frequencyobs[currentset][fluxorder[currentset][imgnum[currentset]-1]]);

	plotfluxbyregion(regflux, imgnum, fluxorder, frequencyobs, frequencyobs[currentset][fluxorder[currentset][0]], frequencyobs[currentset][fluxorder[currentset][imgnum[currentset]-1]], regminflux, regmaxflux, regnumber, fluxerror, maptitle, titles, labels, axismarks, 0, currentset, uselogx, uselogy, autoscalex, autoscaley, userminx, usermaxx, userminy, usermaxy, printinc, skip, symbol, output);

      }

      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;


    case 25: // Setup auto scaling

      changevalues = 0;
      escape = 0;
      setvalues = 1;

      // If it's already off, turn it back on again

      if (autoscalex == 0) {

	validentry = 0;
	yesno = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Turn on autoscaling of the X axes? (No [0], Yes [1], -1 to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  yesno = strtol(cmdbuffer, &endptr, 10);

	  if ( (yesno == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (yesno < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if (yesno > 1) {
	    printf("\nInvalid selection, please choose again...\n\n");
	    continue;
	  }

	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if (yesno == 1) {
	  autoscalex = 1;
	  setvalues = 0;
	  printf("Auto scaling of the X axis is now on\n");
	}
	else {
	  changevalues = 1;
	  setvalues = 1;
	}

      }

      // Turn off the X axis if needed then get the user defined values

      if ( ( (autoscalex == 1) || (changevalues == 1) ) && (setvalues == 1) ) {

	validentry = 0;
	yesno = 1;

	if (changevalues == 0) {

	  while (validentry == 0) {

	    cmdbuffer = readline("Turn off autoscaling of the X axes? (No [0], Yes [1], -1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    yesno = strtol(cmdbuffer, &endptr, 10);

	    if ( (yesno == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (yesno < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }

	    else if (yesno > 1) {
	      printf("\nInvalid selection, please choose again...\n\n");
	      continue;
	    }

	    else {
	      validentry = 1;
	    }
	  }

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }

	  if (yesno == 0) {
	    printf("Auto scaling of the X axis remains on\n");
	  }
	}

	if (yesno == 1) {

	  // Get the min X axis values

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Please enter a minimum X axis value: ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    userminx = strtod(cmdbuffer, &endptr);

	    if ( (userminx == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }
	    else {
	      validentry = 1;
	    }
	  }

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Please enter a maximum X axis value: ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    usermaxx = strtod(cmdbuffer, &endptr);

	    if ( (usermaxx == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }
	    else {
	      validentry = 1;
	    }
	  }

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }

	  autoscalex = 0;

	  printf("Now using user defined values for the X axis!\n");
      
	}
      }


      // Now for the Y axis

      changevalues = 0;
      escape = 0;
      setvalues = 1;

      if (autoscaley == 0) {

	validentry = 0;
	yesno = 0;
	  
	while (validentry == 0) {

	  cmdbuffer = readline("Turn on autoscaling of the Y axes? (No [0], Yes [1], -1 to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  yesno = strtol(cmdbuffer, &endptr, 10);

	  if ( (yesno == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (yesno < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if (yesno > 1) {
	    printf("\nInvalid selection, please choose again...\n\n");
	    continue;
	  }

	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if (yesno == 1) {
	  printf("Auto scaling of the Y axis is now on\n");
	  autoscaley = 1;
	  setvalues = 0;
	}

	else {
	  changevalues = 1;
	  setvalues = 1;
	}
      }

      // And the other half...

      if ( ( (autoscaley == 1) || (changevalues == 1) ) && (setvalues == 1)  ) {

	validentry = 0;
	yesno = 1;

	if (changevalues == 0) {

	  while (validentry == 0) {

	    cmdbuffer = readline("Turn off autoscaling of the Y axes? (No [0], Yes [1], -1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    yesno = strtol(cmdbuffer, &endptr, 10);

	    if ( (yesno == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (yesno < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else if (yesno > 1) {
	      printf("\nInvalid selection, please choose again...\n\n");
	      continue;
	    }
	    else {
	      validentry = 1;
	    }
	  }

	  if (escape == 1) {
	    break;
	  }

	  if (yesno == 0) {
	    printf("Auto scaling of the Y axis remains on\n");
	  }
	}

	if (yesno == 1) {

	  // Get the min Y axis value
	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Please enter a minimum Y axis value: ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    userminy = strtod(cmdbuffer, &endptr);

	    if ( (userminy == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }
	    else {
	      validentry = 1;
	    }
	  }

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Please enter a maximum Y axis value: ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    usermaxy = strtod(cmdbuffer, &endptr);

	    if ( (usermaxy == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }
	    else {
	      validentry = 1;
	    }
	  }

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }

	  autoscaley = 0;

	  printf("Now using user defined values for the Y axis!\n");
	}
      }

      break;


    case 26: // Turn on and off printing of regions with an increasing flux

      if (printinc == 1) {
	printinc = 0;
	printf("Printing of regions with an increasing flux with incresing frequency is now turned off!\n");
      }
      else if (printinc == 0) {
	printinc = 1;
	printf("Printing of regions with an increasing flux with incresing frequency is now turned on!\n");
      }
      else {
	printf("\nThe value of printinc is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", printinc);
      }

      break;
      

    case 27: // Turn on and off the use of logarithmic values for plots

      validentry = 0;
      escape = 0;

      while (validentry == 0) {

	cmdbuffer = readline("Use logarithmic values for the X axis? (No [0], Yes [1], Both No [2], Both Yes [3], -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	selectlog = strtol(cmdbuffer, &endptr, 10);

	if ( (selectlog == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (selectlog < 0) {
	  printf("Escaping command...\n");
	  escape = 1;
	  break;
	}

	else if (selectlog > 3) { // Check the set exists
	  printf("\nInvalid selection, please try again...\n\n");
	  continue;
	}
	else {
	  validentry = 1;
	}
      }

      // Break the second loop if they want to escape
      if (escape == 1) {
	break;
      }

      // Set the value of uselogx
      if (selectlog == 0) {
	uselogx = 0;
	printf("Use of logarithmic values is now turned off for the X axis!\n");
      }
      else if (selectlog == 1) {
	uselogx = 1;
	printf("Use of logarithmic values is now turned on for the X axis!\n");
      }
      else if (selectlog == 2) {
	uselogx = 0;
	uselogy = 0;
	printf("Use of logarithmic values is now turned off for both the X and Y axes!\n");
	break;
      }
      else if (selectlog == 3) {
	uselogx = 1;
	uselogy = 1;
	printf("Use of logarithmic values is now turned on for both the X and Y axes!\n");
	break;
      }
      else {
	printf("\nThe value of selectlog is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", selectlog);
      }


      // Now for the Y axis if not already set

      validentry = 0;
      escape = 0;

      while (validentry == 0) {

	cmdbuffer = readline("Use logarithmic values for the Y axis? (No [0], Yes [1], -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	selectlog = strtol(cmdbuffer, &endptr, 10);

	if ( (selectlog == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (selectlog < 0) {
	  printf("Escaping command...\n");
	  escape = 1;
	  break;
	}

	else if (selectlog > 1) { // Check the set exists
	  printf("\nInvalid selection, please try again...\n\n");
	  continue;
	}
	else {
	  validentry = 1;
	}
      }

      if (escape == 1) {
	break;
      }

      // Set the value of uselogy
      if (selectlog == 0) {
	uselogy = 0;
	printf("Use of logarithmic values is now turned off for the Y axis!\n");
      }
      else if (selectlog == 1) {
	uselogy = 1;
	printf("Use of logarithmic values is now turned on for the X axis!\n");
      }
      else {
	printf("\nThe value of selectlog is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", selectlog);
      }

      break;


    case 28: // Set the value of how many regions to skip in a plot

      validentry = 0;
      tmp_skip = skip;

      printf("Skip currently set to %d (DEFAULT: 1)\n", skip);


      while (validentry == 0) {

	cmdbuffer = readline("How many regions should be skipped? (1 to plot all regions, -1 to escape, integers only): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	skip = strtol(cmdbuffer, &endptr, 10);

	if ( (skip == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (skip < 0) {
	  printf("Escaping command...\n");
	  skip = tmp_skip;
	  break;
	}

	if (skip == 0) {
	  printf("\nDividing by 0 may cause an end to our Universe. Please try another value...\n\n");
	  continue;
	}
	else {
	  printf("Plots will now only output once every %d regions\n", skip);
	  validentry = 1;
	}
      }

      break;


    case 29: // Calculate curvature values via a polynomial fit

      if (numdatasets > 0) {
	
	if (export == 1) {

	  // Open the directory that is currently set and check it exists
	  dir = opendir(imageloc);

	  if (dir != NULL) {
	    closedir (dir); // If we find it move on without comment
	  }
	  else {
	    printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	    break;
	  }
	}

	
	// List the available datasets
	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to fit the polynomial (-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	 currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
	  else if ( (imgnum[currentset] <= 1) && (currentset != 888) ) { // Check the set has at least two maps
	    fprintf(stderr,"\nError: The data set must contain a minimum of two images for fitting. Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    continue;
	  }

	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	// Set the array memory

	if (curvearrayset == 0) {

	  curvearray = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  mincurve = (float *)calloc(MAXDATASETS, sizeof(float));
	  maxcurve = (float *)calloc(MAXDATASETS, sizeof(float));
	  curvearrayset = 1;
	}

	if (currentset == 888) {

	  looplimit = numdatasets;
	  extwincontrol = 1;

	  if (export != 1) {
	    sprintf(output,"/xs");
	    cpgopen(output);
	  }

	  // If all maps check for memory assignments then free and allocate as appropriate

	  for (j=0; j<looplimit; j++){

	    if (curvatureset[j] == 1) {
	      free(curvearray[currentset]);
	    }

	    if (curvatureset[j] != 1) {
	      curvearray[j] = (float **)calloc(regnumber[j], sizeof(float *));

	      for(i=0;i<=regnumber[j];i++) {
		curvearray[j][i] = (float *)calloc(poly+1, sizeof(float));
	      }
	    }
	  }
	}

	else {
	  // Else if just 1 map check for memoery assignments then free and allocate as appropriate
    
	  looplimit = 1;

	  if (curvatureset[currentset] == 1) {
	    free(curvearray[currentset]);
	  }

	  curvearray[currentset] = (float **)calloc(regnumber[currentset], sizeof(float *));

	  for(i=0;i<=regnumber[currentset];i++) {
	    curvearray[currentset][i] = (float *)calloc(poly+1, sizeof(float));
	  }
	}
    

	for (j=0; j<looplimit; j++){

	  if ( (looplimit > 1) || ( (looplimit == 1) && (currentset == 888) ) ) {
	    currentset = j;
	  }


	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets

	  if (regionsset[currentset] != 1) {
	    printf("\nThe regions for data set %d have not been set! Skipping this set...\n", j);
	    continue;	  
	  }
	  else if (imgnum[currentset] <= 1) { // Check the set exists
	    fprintf(stderr,"\n*** Warning: Data set %d does not contain the minimum two images required for fitting. Skipping this set... ***\n\n", j);
	    continue;
	  }
	  
	  // Create the temporary arrays

	  x = (float *)calloc(imgnum[currentset], sizeof(float));
	  y = (float *)calloc(imgnum[currentset], sizeof(float));
	  coeff = (float *)calloc(poly+1, sizeof(float));

	  mincurve[currentset] = 1e20;
	  maxcurve[currentset] = -1e20;


	  for(i=1;i<=regnumber[currentset];i++) {

	    for (a = 0; a < imgnum[currentset]; a++){
    
	      x[a] = log10(frequency[currentset][fluxorder[currentset][a]]);
	      y[a] = log10(regflux[currentset][fluxorder[currentset][a]][i]);
	    }

	    curvature(imgnum[currentset], x, y, poly, coeff, printpoly);
  
	    //Find the  min and max curvature values
	    //Reversing the axis, may want to add an option for this at a later date


	    if (revcurve == 1) {

	      
	      if (-coeff[poly] < mincurve[currentset]) {

		mincurve[currentset] = (-coeff[poly]);
		//printf("Mincurve is now: %.2f \n", mincurve[currentset]);
	      }

	      if (-coeff[poly] > maxcurve[currentset]) {
		maxcurve[currentset] = (-coeff[poly]);
	
		//printf("Maxcurve is now: %.2f \n", maxcurve[currentset]);
	      }
	    }

	    else {

	      if (coeff[poly] < mincurve[currentset]) {

		mincurve[currentset] = (coeff[poly]);
		//printf("Mincurve is now: %.2f \n", mincurve[currentset]);
	      }

	      if (coeff[poly] > maxcurve[currentset]) {
		maxcurve[currentset] = (coeff[poly]);
	
		//printf("Maxcurve is now: %.2f \n", maxcurve[currentset]);
	      }
	    }
	
	    // Set the curvature array values

	    for (a=0; a<=poly; a++){

	      if (revcurve == 1) {

		curvearray[currentset][i][a] = -coeff[a];
	      }

	      else {

		curvearray[currentset][i][a] = coeff[a];
	      }
	    }

	  }
   
	  curvatureset[currentset] = 1;

	  // Free up the temporary arrays
	  free(x);
	  free(y);
	  free(coeff);

	  // If exporting, output open a new session for each with a different name
	  if (export == 1) {

	    extwincontrol = 1; // We always want external control if exporting

	    time( &currenttime );
	    time_struct = localtime ( &currenttime );

	    strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	    sprintf(output,"%s/%s_PolyFit_%s.%s/%s", imageloc, settarget[currentset], timebuff, imagetype, imagetype);

	    printf("\nExporting as %s\n", output);

	    cpgopen(output);
	  }
	  else {
	    sprintf(output,"/xs");
	  }
      
	  plotpolyfit(regnumber, currentset, imgnum, frequencyobs, fluxorder, plotres, curvearray, mincurve, maxcurve, poly, userminx, usermaxx, userminy, usermaxy, autoscalex, autoscaley, skip, titles, labels, axismarks, settarget[currentset], extwincontrol, output);

	  if ( export == 1 ) {
	    cpgclos();
	  }

	}


	if (extwincontrol == 1) {
	  extwincontrol = 0;
	}

	if ( (export != 1) && (extwincontrol == 1) ) {
	  cpgclos();
	}

      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;


    case 30: // Set the value of the plot resolution

      validentry = 0;
      tmp_plotres = plotres;

      printf("Plotting resolution currently set to %d (DEFAULT: %d)\n", plotres, DEFAULTPLOTRES);

      while (validentry == 0) {

	cmdbuffer = readline("Enter the number of points to use when plotting curves (integers only): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	plotres = strtol(cmdbuffer, &endptr, 10);

	if ( (plotres == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (plotres < 0) {
	  printf("Escaping command...\n");
	  plotres = tmp_plotres;
	  break;
	}

	if ( (plotres >= 0) && (plotres < 10) ) {
	  printf("\nValue must be greater than 10. Please try again...\n\n");
	  continue;
	}
	else {
	  printf("Plotting of curves will now use %d points\n", plotres);
	  validentry = 1;
	}
      }

      break;


    case 31: // Set the order polynomial to use

      validentry = 0;
      tmp_poly = poly;

      printf("Polynomial order currently set to %d (DEFAULT: %d)\n", poly, DEFAULTPOLY);

      while (validentry == 0) {

	cmdbuffer = readline("Enter the order polynomial to use for fitting (integers only, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	poly = strtol(cmdbuffer, &endptr, 10);

	if ( (poly == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else if (poly < 0) {
	  printf("Escaping command...\n");
	  poly = tmp_poly;
	  break;
	}
	else if (poly < 2) {
	  printf("\nPolynomial must by atleast quadratic (poly = 2). Please try again...\n\n");
	  //poly = tmp_poly;
	  continue;
	}
	else if (poly == 2) {
	  printf("Fitpoly command will now fit a quadratic to the data.\n");
	  validentry = 1;
	}
	else if (poly == 3) {
	  printf("Fitpoly command will now fit a cubic to the data.\n");
	  validentry = 1;
	}
	else if (poly > 3) {
	  printf("Fitpoly command will now fit a %d th order polynomial to the data.\n", poly);
	  validentry = 1;
	}
      }

      break;


    case 32: // Change the sign convention use when fitting curves

      if (revcurve == 1) {
	revcurve = 0;
	printf("Switching of the sign convention used when fitting curves has now been changed to OFF!\n");
      }
      else if (revcurve == 0) {
	revcurve = 1;
	printf("Switching of the sign convention used when fitting curves has now been changed to ON!\n");
      }
      else {
	printf("\nThe value of revcurves is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", revcurve);
      }

      break;


    case 33: // Plot the total flux of a source

      if (numdatasets > 0) {

	if (export == 1) {

	  // Open the directory that is currently set and check it exists
	  dir = opendir(imageloc);

	  if (dir != NULL) {
	    closedir (dir); // If we find it move on without comment
	  }
	  else {
	    printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	    break;
	  }
	}

	// List the available datasets
	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to plot the total flux (-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset= strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	// Setup the loops

	if (currentset == 888) {
	  looplimit = numdatasets;
	  extwincontrol = 1;

	  if (export != 1) {
	    sprintf(output,"/xs");
	    cpgopen(output);
	  }

	}
	else {
	  looplimit = 1;
	}

	for (j=0; j<looplimit; j++){

	  if ( (looplimit > 1) || ( (looplimit == 1) && (currentset == 888) ) ) {
	    currentset = j;
	  }

	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets

	  if (regionsset[currentset] != 1) {
	    printf("\nThe regions for data set %d have not been set! Skipping this set...\n", j);
	    continue;	  
	  }

	  // If exporting, output open a new session for each with a different name
	  if (export == 1) {

	    extwincontrol = 1; // We always want external control if exporting

	    time( &currenttime );
	    time_struct = localtime ( &currenttime );

	    strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	    sprintf(output,"%s/%s_TotalFlux_%s.%s/%s", imageloc, settarget[currentset], timebuff, imagetype, imagetype);

	    printf("\nExporting as %s\n", output);

	    cpgopen(output);
	  }
	  else {
	    sprintf(output,"/xs");
	  }

	  plottotalflux(regflux[currentset], frequencyobs[currentset], fluxorder[currentset], rms[currentset], imgnum[currentset], regnumber[currentset], symbol, settarget[currentset], extwincontrol, fluxcalerror[currentset], bestguess, specrange, normrange, plotres, specres, fluxres, regionsize[currentset], onsourcenoise, settarget[currentset], titles, labels, axismarks, output);

	  if ( export == 1 ) {
	    cpgclos();
	  }
	}

      
	if (extwincontrol == 1) {
	  extwincontrol = 0;
	}

	if ( (export != 1) && (extwincontrol == 1) ) {
	  cpgclos();
	}

      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;


    case 34: // Set the values for the flux calibration errors

      changearange = 0;
      doallmaps = 0;

      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to set the flux calibration error (-1 to escape, integers only, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    //printf("This value is greater than the number of data sets that exists! Please try again...\n");
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}
      
	// All maps or just one?
	validentry = 0;
     
	if (currentset == 888) {
	  printf("===========================================================\n\n");
	  printf(" Dataset |  Index  |  Frequency (obs.)  |  Frequency (Rest)  \n\n");

	  // Loop through each of the sets
	  for (i=0; i<numdatasets; i++) {
	    for (a=0; a<imgnum[i]; a++) {

	      printf("   %d          %d         %.2e Hz          %.2e Hz\n", i, fluxorder[i][a], frequencyobs[i][fluxorder[i][a]], frequency[i][fluxorder[i][a]]);

	    }
	  }
	  printf("\n===========================================================\n");

	}

	else {

	  printf("================================================\n\n");
	  printf(" Index  |  Frequency (Obs.)  |  Frequency (Rest)  \n\n");

	  // Loop through each of the sets
	  for (a=0; a<imgnum[currentset]; a++) {

	    printf("   %d         %.2e Hz          %.2e Hz\n", fluxorder[currentset][a], frequencyobs[currentset][fluxorder[currentset][a]], frequency[currentset][fluxorder[currentset][a]]);

	  }

	  printf("\n================================================\n");
    
	}
    

	while (validentry == 0) {

	  cmdbuffer = readline("Which map(s) would you like to change the flux calibration error for? (-1 to escape, 888 all maps, 777 to select a frequency range): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  mapselect = strtol(cmdbuffer, &endptr, 10);

	  if ( (mapselect == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (mapselect < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  /*else if ( (mapselect >= imgnum[currentset]) && (mapselect != 888) && (mapselect != 777)) { // Check the map exists
	    printf("\nInvalid selection, please try again...\n\n");
	    continue;
	  }*/

	  else if ( (mapselect >= imgnum[currentset]) && (currentset != 888) && (mapselect != 888) && (mapselect != 777)) { // Check the map exists
	    printf("\nInvalid selection, please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset == 888) && (mapselect != 888) && (mapselect != 777) ){ // Check the map exists
	    count = 0;
	    for (i=0; i<numdatasets; i++) {
	      if (mapselect >= imgnum[i]) {
		printf("\n*** Warning: Data set %d does not contain an image with index %d. RMS will not be updated for this data set... ***\n", i, mapselect);
	      }
	      else {
		count++;
	      }
	    }
	    if (count == 0) {
	      printf("\nError: No loaded data set contains a map with index %d, please try again...\n\n", mapselect);
	      continue;
	    }
	    else {
	      validentry = 1;
	    }
	  }

	  else if (mapselect == 888) {
	    doallmaps = 1;
	    validentry = 1;
	  }
	  else if (mapselect == 777) {
	    changearange = 1;
	    doallmaps = 1;
	    validentry = 1;
	  }
	  else {
	    validentry = 1;
	  }
	}
      
	if (escape == 1) {
	  break;
	}

	// Get the frequency range to change if required

	if (mapselect == 777) {

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Greater than (0) or less than (1) a frequency? (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    gtlt = strtol(cmdbuffer, &endptr, 10);

	    if ( (gtlt == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (gtlt < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }

	    else if (gtlt > 1) { // Check the entry is valid
	      printf("\nInvalid selection, please try again...\n\n");
	      continue;
	    }

	    else {
	      validentry = 1;
	    }
	  }
      
	  if (escape == 1) {
	    break;
	  }


	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Enter the (observed) frequency above or below which this error occurs in Hz (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    errorchange = strtod(cmdbuffer, &endptr);

	    if ( (errorchange == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (errorchange < 0.0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else {
	      validentry = 1;
	    }
	  }
      
	  if (escape == 1) {
	    break;
	  }

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Enter the percentage error for this selection in decimal form (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    flaterror = strtod(cmdbuffer, &endptr);

	    if ( (flaterror == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (flaterror < 0.0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else if (flaterror > 0.201) { // Check the entry is valid
	      printf("\n ***WARNING: Error is above 20%%. This seems very large!*** \n");
	      validentry = 1;
	    }
	    else {
	      validentry = 1;
	    }
	  }
      
	  if (escape == 1) {
	    break;
	  }

	}

	else {

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Enter the percentage error to apply in decimal form (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    flaterror = strtod(cmdbuffer, &endptr);

	    if ( (flaterror == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (flaterror < 0.0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }

	    else if (flaterror > 0.201) { // Check the entry is valid
	      printf("\n ***WARNING: Error is above 20%%. This seems very large!*** \n");
	      validentry = 1;
	    }

	    else {
	      validentry = 1;
	    }
	  }
      
	  if (escape == 1) {
	    break;
	  }

	}


	// Setup the loops
	if (currentset == 888) {
	  looplimit = numdatasets;
	}
	else {
	  looplimit = 1;
	}


	for (i=0; i<looplimit; i++){

	  if (looplimit > 1) {
	    currentset = i;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = i;
	  }


	  if (doallmaps == 1) {
	    looplimit2 = imgnum[currentset];
	  }
	  else {
	    looplimit2 = 1;
	  }


	  for (j=0; j<looplimit2; j++){

	    if (looplimit2 > 1) {
	      mapselect  = j;
	    }
	    else if ( (looplimit2 == 1) && ( (mapselect == 888) || (mapselect == 777) ) ) {
	      mapselect  = j;
	    }


	    // Do the changes if a range is requested
	    if (changearange == 1) {

	      if ( (gtlt == 0) && (frequencyobs[currentset][fluxorder[currentset][mapselect]] > errorchange) ) {

		fluxcalerror[currentset][fluxorder[currentset][mapselect]] = flaterror;
	      }
	      else if ( (gtlt == 1) && (frequencyobs[currentset][fluxorder[currentset][mapselect]] < errorchange) ) {

		fluxcalerror[currentset][fluxorder[currentset][mapselect]] = flaterror;
	      }
	    }

	    else {
	      fluxcalerror[currentset][mapselect] = flaterror;
	    }

	    //printf("Flux calibration error for dataset %d at %.2e Hz now set to %.2f%%\n", currentset, frequency[currentset][mapselect], flaterror*100);

	  }

	}

	printf("\nFlux calibration errors now set. Please re-run the appropriate region selection command to update the region errors.\n\n");

      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded! If you wish to change the default flux calibration error used on loading of a dataset, please use the instrument command.\n");
      }
    
      break;


    case 35: // Set the value of the plot resolutionbest guess spectral index

      validentry = 0;

      printf("Spectral index best guess currently set to %.2f (DEFAULT: %.2f)\n", bestguess, DEFAULTBESTGUESS);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the spectral index best guess: ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	bestguess = strtod(cmdbuffer, &endptr);

	if ( (bestguess == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  printf("Best guess now set to %.2f \n", bestguess);
	  validentry = 1;
	}
      }

      break;


    case 36: // Set the search range for the spectral index

      validentry = 0;
      tmp_specrange = specrange;

      printf("Spectral idnex search range currently set to %.2f (DEFAULT: %.2f)\n", specrange, DEFAULTSPECRANGE);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the spectral index range (+- absolute value, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	specrange = strtod(cmdbuffer, &endptr);

	if ( (specrange == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (specrange < 0) {
	  printf("Escaping command...\n");
	  specrange  = tmp_specrange;
	  break;
	}

	else {
	  printf("Spectral index range now set to %.2f \n", specrange);
	  validentry = 1;
	}
      }

      break;


    case 37: // Set the value of the plot resolution

      validentry = 0;
      tmp_normrange = normrange;

      printf("Normalisation search range currently set to %.2f (DEFAULT: %.2f)\n", normrange, DEFAULTNORMRANGE);


      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the flux normalisation range (+- percentage, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	normrange = strtod(cmdbuffer, &endptr);

	if ( (normrange == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (normrange < 0) {
	  printf("Escaping command...\n");
	  normrange  = tmp_normrange;
	  break;
	}

	else {
	  printf("Flux normalisation range now set to %.2f \n", normrange);
	  validentry = 1;
	}
      }

      break;


    case 38: // Set the value of the plot resolution

      validentry = 0;
      tmp_symbol = symbol;

      printf("Symbol currently set to %d (DEFAULT: %d)\n", symbol, DEFAULTSYMBOL);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the symbol to be used in plots (-1 to escape, 0 to 126 inclusive): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	symbol = strtol(cmdbuffer, &endptr, 10);

	if ( (symbol == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (symbol < 0) {
	  printf("Escaping command...\n");
	  symbol  = tmp_symbol;
	  break;
	}

	else if (symbol > 126) {
	  printf("\nValue must be between 0 and 126 inclusive. Please try again...\n\n");
	  continue;
	}

	else {
	  printf("Symbol value now set to %d \n", symbol);
	  validentry = 1;
	}
      }

      break;


    case 39: // Set the value of the spectral index fitting resolution

      validentry = 0;
      tmp_specres = specres;

      printf("Spectral index fitting resolution currently set to %d (DEFAULT: %d)\n", specres, DEFAULTSPECRES);


      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the best fit spectral index search resolution (-1 to escape, min = 10): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	specres = strtol(cmdbuffer, &endptr, 10);

	if ( (specres == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (specres < 0) {
	  printf("Escaping command...\n");
	  specres  = tmp_specres;
	  break;
	}

	if ( (specres >= 0) && (specres < 10) ) {
	  printf("\nValue must be greater than or equal to 10\n\n");
	  continue;
	}

	else {
	  printf("Spectral index search resolution now set to %d \n", specres);
	  validentry = 1;
	}
      }

      break;


    case 40: // Set the value of the normalisation search range for the total flux

      validentry = 0;
      tmp_fluxres = fluxres;

      printf("Total flux normalisation search resolution currently set to %d (DEFAULT: %d)\n", fluxres, DEFAULTFLUXRES);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the total flux normalisation search resolution (-1 to escape, min = 10): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	fluxres = strtol(cmdbuffer, &endptr, 10);

	if ( (fluxres == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (fluxres < 0) {
	  printf("Escaping command...\n");
	  fluxres  = tmp_fluxres;
	  break;
	}

	if ( (fluxres >= 0) && (fluxres < 10) ) {
	  printf("\nValue must be greater than or equal to 10\n\n");
	  continue;
	}
	else {
	  printf("Total flux normalisation search resolution now set to %d \n", fluxres);
	  validentry = 1;
	}
      }

      break;


    case 41: // Plot the JP model for an arbitary normalisation

      printf("\nPlotting an example JP model between %.2e Hz and %.2e Hz\n\n", minmodelfreq, maxmodelfreq);

      model = 1;

      if (export == 1) {

	// Open the directory that is currently set and check it exists
	dir = opendir(imageloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	  break;
	}

	time( &currenttime );
	time_struct = localtime ( &currenttime );

	strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	sprintf(output,"%s/ExJPModel_%s.%s/%s", imageloc,  timebuff, imagetype, imagetype);

	printf("\nExporting as %s\n", output);

	cpgopen(output);
      }
      else {
      sprintf(output,"/xs");
      }

      plotmodel(minmodelfreq, maxmodelfreq, inject, fieldstrength, model, usr_gmin, usr_gmax, minmodelmyears, modelmyears, output, titles, exampleredshift, skip);

      break;


    case 42: // Fit the JP model to a dataset
  
      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to fit the JP model (-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }
	  
	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	printf("Escaping command...\n");
	escape = 1;
	break;
	  }
	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
	  else if ( (imgnum[currentset] <= 1) && (currentset != 888) ) { // Check the set exists
	    fprintf(stderr,"\nError: The data set must contain a minimum of two images for fitting. Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    continue;
	  }

	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	// Setup the loops
	if (currentset == 888) {
	  looplimit = numdatasets;
	}
	else {
	  looplimit = 1;
	}

	for (j=0; j<looplimit; j++){

	  if (looplimit > 1) {
	    currentset = j;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = j;
	  }

	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets

	  if (regionsset[currentset] != 1) {
	    printf("\nThe regions for data set %d have not been set! Skipping this set...\n", j);
	    continue;	  
	  }
	  else if (imgnum[currentset] <= 1) { // Check the set has at least two maps
	    fprintf(stderr,"\n*** Warning: Data set %d does not contain the minimum two images required for fitting. Skipping this set... ***\n\n", j);
	    continue;
	  }
	  
	  setmodelres[currentset][1] = modelres; // Set the resolution for this dataset and model

	  // Setup the memory
	  if (jpmodelmemoryset == 0) {

	    jpchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jpbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jpmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jpmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	    jpbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jpageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jpageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));

	    jpmodelmemoryset = 1;
	  }
	
	  if (numdatasets > 1) {

	    allocarray = (int *)calloc(MAXDATASETS, sizeof(int));

	    realloc_jpchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_jpbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_jpmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_jpmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	    realloc_jpbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_jpageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_jpageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));

	    for (o=0; o<numdatasets; o++) {

	      if (jpfitmemset[o] == 1) {

		realloc_jpchisquared[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_jpbestage[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_jpmodelalpha[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_jpmodelflux[o] = (float **)calloc(regnumber[o]+1, sizeof(float *));
		realloc_jpbestnorm[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_jpageerrorsplus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_jpageerrorsminus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));

		allocarray[o] = 1;

	      }
	    }
	  
	    for (o=0; o<numdatasets; o++) {

	      if (jpfitmemset[o] == 1) {

		for (q=0; q<=regnumber[o]; q++) {

		  //printf("q=%d\n", q);
		  realloc_jpmodelflux[o][q] = (float *)calloc(setmodelres[o][1]+1, sizeof(float));
		}
	      }
	    }

	    for (o=0; o<numdatasets-1; o++) { // -1 as we are yet to allocate the new dataset memory (or we will overwrite it with the new fitting anyway)

	      for (p=0; p<=regnumber[o]; p++) {

		if (jpfitmemset[o] == 1) {

		  realloc_jpchisquared[o][p] = jpchisquared[o][p];
		  realloc_jpbestage[o][p] = jpbestage[o][p];
		  realloc_jpmodelalpha[o][p] = jpmodelalpha[o][p];
		  realloc_jpbestnorm[o][p] = jpbestnorm[o][p];
		  realloc_jpageerrorsplus[o][p] = jpageerrorsplus[o][p];
		  realloc_jpageerrorsminus[o][p] = jpageerrorsminus[o][p];

		  for (q=0; q<=setmodelres[o][1]; q++) {
		    realloc_jpmodelflux[o][p][q] = jpmodelflux[o][p][q];
		  }

		}
	      }
	    }
	  
	    for (o=0; o<numdatasets-1; o++) { // -1 as we are yet to allocate the new dataset memory

	      free(jpchisquared[o]);
	      free(jpbestage[o]);
	      free(jpmodelalpha[o]);
	      free(jpbestnorm[o]);
	      free(jpageerrorsplus[o]);
	      free(jpageerrorsminus[o]);
	      free(jpmodelflux[o]);

	    }
	    
	    free(jpchisquared);
	    free(jpbestage);
	    free(jpmodelalpha);
	    free(jpbestnorm);
	    free(jpageerrorsplus);
	    free(jpageerrorsminus);
	    free(jpmodelflux);
	    
	    jpfitmemset[o] = 0;
	    jpmodelmemoryset = 0;


	    // Make the first level of the array again
	    jpchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jpbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jpmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jpmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	    jpbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jpageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jpageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));

	    jpmodelmemoryset = 1;

	  }
	
	  for (o=0; o<numdatasets; o++) {

	    jpchisquared[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    jpbestage[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    jpmodelalpha[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    jpmodelflux[o] = (float **)calloc(regnumber[o]+1, sizeof(float *));
	    jpbestnorm[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    jpageerrorsplus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    jpageerrorsminus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));

	    for (p=0; p<=regnumber[o]; p++) {
	      jpmodelflux[o][p] = (float *)calloc(setmodelres[o][1]+1, sizeof(float));
	    }
	    
	    jpfitmemset[o] = 1;

	  }
	
	
	  if (numdatasets > 1) {
	    for (o=0; o<numdatasets; o++) {

	      if ( (o != currentset) && (allocarray[o] == 1) ) {

		for (p=0; p<=regnumber[o]; p++) {

		  jpchisquared[o][p] = realloc_jpchisquared[o][p];
		  jpbestage[o][p] = realloc_jpbestage[o][p];
		  jpmodelalpha[o][p] = realloc_jpmodelalpha[o][p];
		  jpbestnorm[o][p] = realloc_jpbestnorm[o][p];
		  jpageerrorsplus[o][p] = realloc_jpageerrorsplus[o][p];
		  jpageerrorsminus[o][p] = realloc_jpageerrorsminus[o][p];

		  for (q=0; q<=setmodelres[o][1]; q++) {

		    jpmodelflux[o][p][q] = realloc_jpmodelflux[o][p][q];
		  }
		}
	      }
	    }
	  }
	
	  if (numdatasets > 1) {

	    // Free the reallocation arrays if they are set

	    for (o=0; o<numdatasets-1; o++) {

	      free(realloc_jpchisquared[o]);
	      free(realloc_jpbestage[o]);
	      free(realloc_jpmodelalpha[o]);
	      free(realloc_jpbestnorm[o]);
	      free(realloc_jpageerrorsplus[o]);
	      free(realloc_jpageerrorsminus[o]);
	      //free(realloc_jpmodelflux[o]);

	    }

	    free(realloc_jpchisquared);
	    free(realloc_jpbestage);
	    free(realloc_jpmodelalpha);
	    free(realloc_jpbestnorm);
	    free(realloc_jpageerrorsplus);
	    free(realloc_jpageerrorsminus);
	    free(realloc_jpmodelflux);

	    free(allocarray);

	  }

	  model = 1;

	  // Using a temp array here to allow all models without the need for dummy arrays
	  tmp_ageerrorsplus = (float *)calloc(regnumber[currentset]+1, sizeof(float));
	  tmp_ageerrorsminus = (float *)calloc(regnumber[currentset]+1, sizeof(float));


	  spectralageingmodels(regflux[currentset], imgnum[currentset], fluxorder[currentset], frequency[currentset], regnumber[currentset], jpchisquared[currentset], jpbestage[currentset], fluxerror[currentset], jpmodelalpha[currentset], jpmodelflux[currentset], usr_gmin, usr_gmax, minmyears, myears, ageresolution, levels, jpbestnorm[currentset], setmodelres[currentset][1], printresults, inject, fieldstrength, model, printreject, redshift[currentset], tmp_ageerrorsplus, tmp_ageerrorsminus, suppresscdf);

	  jpfield[currentset] = fieldstrength;
	  jpinject[currentset] = inject;
	  jpmodelset[currentset] = 1;

	  // Transfer the age errors over to the perminant array
	  for (r=1; r<=regnumber[currentset]; r++) {
	    jpageerrorsplus[currentset][r] = tmp_ageerrorsplus[r];
	    jpageerrorsminus[currentset][r] = tmp_ageerrorsminus[r];
	  }

	  free(tmp_ageerrorsplus);
	  free(tmp_ageerrorsminus);  

	}

      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;


    case 43: // Set the value of the multiplier for the onsource noise

      validentry = 0;
      tmp_onsourcenoise = onsourcenoise;

      printf("On-source noise multiplier currently set to %.2f (DEFAULT: %.2f)\n", onsourcenoise, DEFAULTONSOURCENOISE);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a multiplier for the on-source noise (RMS x n) (-1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	onsourcenoise = strtod(cmdbuffer, &endptr);

	if ( (onsourcenoise == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (onsourcenoise < 0) {
	    printf("Escaping command...\n");
	    onsourcenoise = tmp_onsourcenoise;
	    break;
	  }
	  else {
	    printf("\nOn-source noise multiplier now set to %.2f. Regions will need to be remade for this to be applied! (excluding plottotalflux)\n\n", onsourcenoise);
	    validentry = 1;
	  }
	}
      }

      break;


    case 44: // Set the value of the minimum model frequency

      validentry = 0;
      tmp_minmodelfreq = minmodelfreq;

      printf("Minimum example model frequency currently set to %.2e Hz (DEFAULT: %.2e Hz)\n", minmodelfreq, DEFAULTMINMODELFREQ);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a minimum frequency in Hz for the plotting of example models (-1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	minmodelfreq = strtod(cmdbuffer, &endptr);

	if ( (minmodelfreq == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {

	  if (minmodelfreq < 0) {
	    printf("Escaping command...\n");
	    minmodelfreq = tmp_minmodelfreq;
	    break;
	  }
	  else {
	    printf("\nMinimum model plotting frequency now set to %.2e Hz\n\n", minmodelfreq);
	    validentry = 1;
	  }
	}
      }

      break;


    case 45: // Set the value of the maximum model frequency

      validentry = 0;
      tmp_maxmodelfreq = maxmodelfreq;

      printf("Maximum example model frequency currently set to %.2e (DEFAULT: %.2e)\n", maxmodelfreq, DEFAULTMAXMODELFREQ);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a maximum frequency in Hz for the plotting of example models (-1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	maxmodelfreq = strtod(cmdbuffer, &endptr);

	if ( (maxmodelfreq == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}	
	else {

	  if (maxmodelfreq < 0) {
	    printf("Escaping command...\n");
	    maxmodelfreq = tmp_maxmodelfreq;
	    break;
	  }
	  else {
	    printf("\nMaximum model plotting frequency now set to %.2e Hz\n\n", maxmodelfreq);
	    validentry = 1;
	  }
	}
      }

      break;


    case 46: // Set the value of the age resolution

      validentry = 0;
      tmp_ageresolution = ageresolution;

      printf("Age resolution currently set to %d (DEFAULT: %d)\n", ageresolution, DEFAULTAGERES);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the model age resolution (integers only, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	ageresolution = strtol(cmdbuffer, &endptr, 10);

	if ( (ageresolution == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {

	  if (ageresolution < 0) {
	    printf("Escaping command...\n");
	    ageresolution = tmp_ageresolution;
	    break;
	  }
	  else {
	    printf("\nAge resolution is now set to %d\n\n", ageresolution);
	    validentry = 1;
	  }
	}
      }

      break;


    case 47: // Set the value of the model ages to be tested

      validentry = 0;
      tmp_myears = myears;

      printf("MYears currently set to %d (DEFAULT: %d)\n", myears, DEFAULTMYEARS);

      while (validentry == 0) {

	cmdbuffer = readline("Enter the number of MYears over which to test models (integers only, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	myears = strtol(cmdbuffer, &endptr, 10);

	if ( (myears == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {

	  if (myears < 0) {
	    printf("Escaping command...\n");
	    myears = tmp_myears;
	    break;
	  }
	  else {
	    printf("\nMYears is now set to %d\n\n", myears);
	    validentry = 1;
	  }
	}
      }

      break;


    case 48: // Set the value of the normalisation resolution (NOW REDUNDANT! GOLDEN RATIO USED INSTEAD)

      validentry = 0;

      while (validentry == 0) {

	printf("Enter the number of normalisation values over which to test models (-1 to escape, current = %d, DEFAULT = %d): ", normresolution, DEFAULTNORMRES);

	tmp_normresolution = normresolution;

	if (scanf("%d", &normresolution) == 0) {
		 
	  printf("Invalid input, escaping command to protect the existing data...\n");
	  normresolution = tmp_normresolution;
	  break;
	}
	else {

	  if (normresolution < 0) {
	    printf("Escaping command...\n");
	    normresolution = tmp_normresolution;
	    break;
	  }
	  else {
	    printf("\nNormalisation resolution is now set to %d Hz\n\n", normresolution);
	    validentry = 1;
	  }
	}
      }

      break;


    case 49: // Plot curvature against spectral index for a given dataset
      
      if (export == 1) {

	// Open the directory that is currently set and check it exists
	dir = opendir(imageloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	  break;
	}
      }
 
      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to plot curvature against spectral index (-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }
	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (curvatureset[currentset] != 1) ) { // Check the curvature has been set
	    fprintf(stderr,"\nError: Curvature have not yet been calculated for this data set! Please run fitpoly first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (specindexset[currentset] != 1) ) { // Check the spectral index has been set
	    fprintf(stderr,"\nError: Spectral indices have not yet been calculated for this data set! Please run specindex first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (regnumber[currentset] <= 1) ) {
	    printf("\n Error: Dataset must have more than 1 region to produce anything useful!\nExiting command...\n\n");
	    escape = 1;
	    break;
	    fprintf(stderr,"\nError: Dataset must have more than 1 region to produce anything useful!\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if (currentset == 888) {
	  looplimit = numdatasets;
	  extwincontrol = 1;

	  if (export != 1) {
	    sprintf(output,"/xs");
	    cpgopen(output);
	  }

	}
	else {
	  looplimit = 1;
	}
      

	for (j=0; j<looplimit; j++){

	  if (looplimit > 1) {
	    currentset = j;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = j;
	  }

	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets

	  if (regionsset[currentset] != 1) {
	    printf("\nThe regions for data set %d have not been set! Skipping this set...\n", j);
	    continue;	  
	  }
	  if (curvatureset[currentset] != 1) {
	    printf("\nThe curvature for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }
	  if (specindexset[currentset] != 1) {
	    printf("\nThe spectral indices for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;
	  }
	  if ( regnumber[currentset] <= 1 ) {
	    printf("\nDataset %d has 1 or less regions. This is not enough to produce anything useful! Skipping this set...\n", j);
	    continue;
	  }

	  // If exporting, output open a new session for each with a different name
	  if (export == 1) {

	    extwincontrol = 1; // We always want external control if exporting

	    time( &currenttime );
	    time_struct = localtime ( &currenttime );

	    strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	    sprintf(output,"%s/%s_CurveSI_%s.%s/%s", imageloc, settarget[currentset], timebuff, imagetype, imagetype);

	    printf("\nExporting as %s\n", output);

	    cpgopen(output);
	  }
	  else {
	    sprintf(output,"/xs");
	  }

	  plotcurvagainstspec(curvearray[currentset], alpha[currentset], regnumber[currentset], symbol, minsi[currentset], maxsi[currentset], mincurve[currentset], maxcurve[currentset], poly, autoscalex, autoscaley, userminx, usermaxx, userminy, usermaxy, settarget[currentset], titles, labels, axismarks, extwincontrol);

	  if ( export == 1 ) {
	    cpgclos();
	  }

	}

	if (extwincontrol == 1) {
	  extwincontrol = 0;
	}

	if ( (export != 1) && (extwincontrol == 1) ) {
	  cpgclos();
	}

      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;


    case 50: // Plot model flux against observed flux as a function of frequency

      if (export == 1) {

	// Open the directory that is currently set and check it exists
	dir = opendir(imageloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	  break;
	}
      }
      
      
      if (numdatasets > 0) {

	// List the available datasets
	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to compare model and observed values (-1 to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }
	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if (currentset >= numdatasets) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }

	  else if (regionsset[currentset] != 1) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    continue;
	  }

	  else if ( (jpmodelset[currentset] != 1) && (kpmodelset[currentset] != 1) && (jptribmodelset[currentset] != 1) && (kptribmodelset[currentset] != 1) ) { // Check the curvature has been set
	    fprintf(stderr,"\nError: A model has not yet been fitted to this data set! Please run a model fitting command first.\n\n");
	    continue;
	  }

	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}


	// Ask which model to plot and check its validity
	validentry = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Which model would you like to plot? (-1 to escape, 1 JP model, 2 KP model, 3 Tribble (JP), 4 Tribble (KP)): ");

	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  model = strtol(cmdbuffer, &endptr, 10);

	  if ( (model == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (model < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (model > NUMBEROFMODELS) ||  (model == 0) ) { // Check the set exists
	    printf("\nInvalid model number. Please try again...\n\n");
	    continue;
	  }
	  else if ( (model == 1) && (jpmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: JP model vales have not yet been calculated for this data set! Please run fitjpmodel first.\n\n");
	    continue;
	  }
	  else if ( (model == 2) && (kpmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: KP model vales have not yet been calculated for this data set! Please run fitjpmodel first.\n\n");
	    continue;
	  }
	  else if ( (model == 3) && (jptribmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Tribble model (JP) vales have not yet been calculated for this data set! Please run fitjpmodel first.\n\n");
	    continue;
	  }
	  else if ( (model == 4) && (kptribmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Tribble model (JP) vales have not yet been calculated for this data set! Please run fitjpmodel first.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	// Open a window if we are not exporting
	if (export != 1) {
	  sprintf(output,"/xs");
	  cpgopen(output);
	}

	numbermodelfluxes = setmodelres[currentset][model];

	tmp_regflux = (float *)calloc(imgnum[currentset], sizeof(float));
	tmp_modelflux = (float *)calloc(setmodelres[currentset][model]+1, sizeof(float));
	tmp_fluxerror = (float *)calloc(imgnum[currentset], sizeof(float));
	tmp_fluxfreq = (float *)calloc(imgnum[currentset], sizeof(float));
	tmp_modelfreq = (float *)calloc(setmodelres[currentset][model]+1, sizeof(float));
	tmp_fluxerrorplus = (float *)calloc(imgnum[currentset], sizeof(float));
	tmp_fluxerrorminus = (float *)calloc(imgnum[currentset], sizeof(float));

	tmp_minflux = 1e32;
	tmp_maxflux = -1e32;


	// Check that they really want to export this many maps if it is over 10
	if ( (export == 1) && ( ((float)regnumber[currentset]/(float)skip)>10) ) {

	  validentry = 0;

	  printf("Warning: You are about to export %d plots!\n", (int)(regnumber[currentset]/(float)skip) );

	  while (validentry == 0) {


	    cmdbuffer = readline("Are you sure you want to continue? (1 Yes, 0 No): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    warningcheck = strtol(cmdbuffer, &endptr, 10);

	    if ( (warningcheck == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }
		
	    if ( (warningcheck < 0) || (warningcheck > 1) ) { // Check the value is either 1 or 0
	      printf("Value must be either Yes (1) or no (0). Please try again...\n");
	      continue;
	    }

	    else if (warningcheck == 0) { // If 0, exit
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }

	    else if (warningcheck == 1) { // If 1, continue
	      printf("Proceeding...\n");
	      validentry = 1;
	    }

	    else { // Catch anthing odd
	      printf("\nThe value of warningcheck  is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", warningcheck);
	      break;
	    }
	  }

	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}
      
	for (i=skip; i<=regnumber[currentset]; i+=skip) {
	  
	  tmp_minflux = 1e32;
	  tmp_maxflux = -1e32;

	  for (j=0; j<imgnum[currentset]; j++) {

	    tmp_regflux[j] = log10(regflux[currentset][fluxorder[currentset][j]][i]);
	    tmp_fluxfreq[j] = log10(frequencyobs[currentset][fluxorder[currentset][j]]);
	  
	    tmp_fluxerrorplus[j] = log10(regflux[currentset][fluxorder[currentset][j]][i] + fluxerror[currentset][fluxorder[currentset][j]][i]);
	  
	    tmp_fluxerrorminus[j] = log10(regflux[currentset][fluxorder[currentset][j]][i] - fluxerror[currentset][fluxorder[currentset][j]][i]);
	

	    if (pow(10,tmp_fluxerrorminus[j]) < tmp_minflux) {
	      tmp_minflux = pow(10, tmp_fluxerrorminus[j]);
	    }

	    if (pow(10,tmp_fluxerrorplus[j]) > tmp_maxflux) {
	      tmp_maxflux =  pow(10, tmp_fluxerrorplus[j]);
	    }

	  }

	  for (j=0; j<=setmodelres[currentset][model]; j++) {

	    if (model == 1) {
	      tmp_modelflux[j] = log10( jpmodelflux[currentset][i][j] );
	    }

	    else if (model == 2) {
	      tmp_modelflux[j] = log10( kpmodelflux[currentset][i][j] );
	      //printf("tmp_modelflux[j]: %.4e\n", tmp_modelflux[j]);
	    }

	    else if (model == 3) {
	      tmp_modelflux[j] = log10( jptribmodelflux[currentset][i][j] );
	    }

	    else if (model == 4) {
	      tmp_modelflux[j] = log10( kptribmodelflux[currentset][i][j] );
	    }
	    else {
	      fprintf(stderr,"\nError: Model number is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");

	      // Free the temporary variables before leaving
	      free(tmp_fluxfreq);
	      free(tmp_regflux);
	      //free(tmp_modelflux);
	      free(tmp_fluxerror);
	      free(tmp_modelfreq);
	      free(tmp_fluxerrorplus);
	      free(tmp_fluxerrorminus);

	      break;
	    }

	    // Check if the model is the min/max value we need to use for plotting

	    if (tmp_modelflux[j] < tmp_minflux) {
	      tmp_minflux = pow(10, tmp_modelflux[j]);
	    }
	    if (tmp_modelflux[j] > tmp_maxflux) {
	      tmp_maxflux = pow(10, tmp_modelflux[j]);
	    }
	  
	    tmp_modelfreq[j] = log10(frequencyobs[currentset][fluxorder[currentset][0]] + ( ( (frequencyobs[currentset][fluxorder[currentset][imgnum[currentset]-1]] - frequencyobs[currentset][fluxorder[currentset][0]]) / setmodelres[currentset][model] ) * j ) );

	    //printf("tmp_modelfreq[j]: %.4e\n", tmp_modelfreq[j]);

	  }

	  // If exporting, output open a new session for each with a different name
	  if (export == 1) {

	    time( &currenttime );
	    time_struct = localtime ( &currenttime );

	    strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	    if (model == 1) {

	      sprintf(output,"%s/%s_Reg%d_ModObs_JP_%s.%s/%s", imageloc, settarget[currentset], i, timebuff, imagetype, imagetype);
	    }

	    else  if (model == 2) {

	      sprintf(output,"%s/%s_Reg%d_ModObs_KP_%s.%s/%s", imageloc, settarget[currentset], i, timebuff, imagetype, imagetype);
	    }

	    else  if (model == 3) {

	      sprintf(output,"%s/%s_Reg%d_ModObs_TribbleJP_%s.%s/%s", imageloc, settarget[currentset], i, timebuff, imagetype, imagetype);
	    }

	    else  if (model == 4) {

	      sprintf(output,"%s/%s_Reg%d_ModObs_TribbleKP_%s.%s/%s", imageloc, settarget[currentset], i, timebuff, imagetype, imagetype);
	    }

	    else {

	      sprintf(output,"%s/%s_Reg%d_ModObs_UnknownModel_%s.%s/%s", imageloc, settarget[currentset], i, timebuff, imagetype, imagetype);
	    }

	    printf("\nExporting as %s\n", output);

	    cpgopen(output);
	  }
	  else {
	    sprintf(output,"/xs");
	  }

	  if (model == 1) {
	    tmp_flt_inject = jpinject[currentset];
	    tmp_field = jpfield[currentset];
	    tmp_bestage = jpbestage[currentset][i];
	    tmp_chisquared = jpchisquared[currentset][i]/(imgnum[currentset]-2);

	    if (usereduced == 1) {
	      tmp_chisquared = jpchisquared[currentset][i]/(imgnum[currentset]-2);
	    }
	    else {
	      tmp_chisquared = jpchisquared[currentset][i];
	    }
	  }

	  else if (model == 2) {
	    tmp_flt_inject = kpinject[currentset];
	    tmp_field = kpfield[currentset];
	    tmp_bestage = kpbestage[currentset][i];

	    if (usereduced == 1) {
	      tmp_chisquared = kpchisquared[currentset][i]/(imgnum[currentset]-2);
	    }
	    else {
	      tmp_chisquared = kpchisquared[currentset][i];
	    }
	  }

	  else if (model == 3) {
	    tmp_flt_inject = jptribinject[currentset];
	    tmp_field = jptribfield[currentset];
	    tmp_bestage = jptribbestage[currentset][i];

	    if (usereduced == 1) {
	      tmp_chisquared = jptribchisquared[currentset][i]/(imgnum[currentset]-2);
	    }
	    else {
	      tmp_chisquared =jptribchisquared[currentset][i];
	    }
	  }
	  else if (model == 4) {
	    tmp_flt_inject = kptribinject[currentset];
	    tmp_field = kptribfield[currentset];
	    tmp_bestage = kptribbestage[currentset][i];

	    if (usereduced == 1) {
	      tmp_chisquared = kptribchisquared[currentset][i]/(imgnum[currentset]-2);
	    }
	    else {
	      tmp_chisquared = kptribchisquared[currentset][i];
	    }
	  }



	  else {
	    fprintf(stderr,"\nError: Model number is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");
	    break;
	  }


	  // Hack to output model values to a text file
	  /*	  
		  char tmp_filename[256];
		  int t;
		  sprintf(tmp_filename,"%s/modelvalues.txt", dataloc);

		  filestream = fopen(tmp_filename,"w");
		  for (t=0; t<=setmodelres[currentset][1]; t++) {

		  fprintf(filestream, "%.6e %.6e", tmp_modelfreq[t], tmp_modelflux[t]);
		  fprintf(filestream, "\n");
 
		  }
		  fclose(filestream);
	  */

	  // This allows for upper limits to be applied to single injection models later on
	  numupperlimits = (int *)calloc(1, sizeof(int));
	  numupperlimits[0] = 0;
	  upperlimits = (int **)calloc(1, sizeof(int *));
	  upperlimits[0] = (int *)calloc(numupperlimits[0]+1, sizeof(int));
	  upperlimits[0][0] = 0;
	
	  
	  plotmodelvsflux(tmp_modelflux, tmp_regflux, tmp_fluxerrorplus, tmp_fluxerrorminus, symbol, numbermodelfluxes, imgnum[currentset], tmp_fluxfreq,  tmp_modelfreq, autoscalex, autoscaley, userminx, usermaxx, userminy, usermaxy, frequencyobs[currentset][fluxorder[currentset][0]], frequencyobs[currentset][fluxorder[currentset][imgnum[currentset]-1]], tmp_minflux, tmp_maxflux, settarget[currentset], titles, labels, axismarks, tmp_flt_inject, tmp_field, tmp_bestage, 0.0, tmp_chisquared, paramlabels, model, usereduced, numupperlimits[0], upperlimits[0], largetxt);
	  
	  if ( export == 1 ) {
	    cpgclos();
	  }
	
	  if (model == 1) {

	    if (usereduced == 1) {
	      printf("Region: %d Best Age: %.2f MYears Reduced Chi-Squared: %.2f\n", i, jpbestage[currentset][i], tmp_chisquared);
	    }
	    else {
	      printf("Region: %d Best Age: %.2f MYears Chi-Squared: %.2f\n", i, jpbestage[currentset][i], tmp_chisquared);
	    }
	  }
	  else if (model == 2) {

	    if (usereduced == 1) {
	      printf("Region: %d Best Age: %.2f MYears Reduced Chi-Squared: %.2f\n", i, kpbestage[currentset][i], tmp_chisquared);
	    }
	    else {
	      printf("Region: %d Best Age: %.2f MYears Chi-Squared: %.2f\n", i, kpbestage[currentset][i], tmp_chisquared);
	    }
	  }
	  else if (model == 3) {

	    if (usereduced == 1) {
	      printf("Region: %d Best Age: %.2f MYears Reduced Chi-Squared: %.2f\n", i, jptribbestage[currentset][i], tmp_chisquared);
	    }
	    else {
	      printf("Region: %d Best Age: %.2f MYears Chi-Squared: %.2f\n", i, jptribbestage[currentset][i], tmp_chisquared);
	    }
	  }

	  else if (model == 4) {

	    if (usereduced == 1) {
	      printf("Region: %d Best Age: %.2f MYears Reduced Chi-Squared: %.2f\n", i, kptribbestage[currentset][i], tmp_chisquared);
	    }
	    else {
	      printf("Region: %d Best Age: %.2f MYears Chi-Squared: %.2f\n", i, kptribbestage[currentset][i], tmp_chisquared);
	    }
	  }

	}

	if (model == 1) { 
	  printf("JP Model parameters - Injection index: %.2f Magnectic field strength: %.2e T\n\n", jpinject[currentset], jpfield[currentset]);
	}
	else if (model == 2) {
	  printf("KP Model parameters - Injection index: %.2f Magnectic field strength: %.2e T\n\n", kpinject[currentset], kpfield[currentset]);
	}
	else if (model == 3) {
	  printf("Tribble (JP) Model parameters - Injection index: %.2f Magnectic field strength: %.2e T\n\n", jptribinject[currentset], jptribfield[currentset]);
	}
	else if (model == 4) {
	  printf("Tribble (KP) Model parameters - Injection index: %.2f Magnectic field strength: %.2e T\n\n", kptribinject[currentset], kptribfield[currentset]);
	}

	// Free the temporary variables
	free(tmp_fluxfreq);
	free(tmp_regflux);
	//free(tmp_modelflux);
	free(tmp_fluxerror);
	free(tmp_modelfreq);
	free(tmp_fluxerrorplus);
	free(tmp_fluxerrorminus);
	free(upperlimits);
	free(numupperlimits);

	if (export != 1) {
	  cpgclos();
	}
      
      }

      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;


    case 51: // Set the value of the model resolution

      validentry = 0;
      tmp_modelres = modelres;

      printf("Model resolution currently set to %d (DEFAULT: %d)\n", modelres, DEFAULTMODELRES);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the model resolution (integers only, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	modelres = strtol(cmdbuffer, &endptr, 10);

	if ( (modelres == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {

	  if (modelres < 0) {
	    printf("Escaping command...\n");
	    modelres = tmp_modelres;
	    break;
	  }
	  else {
	    printf("\nModel resolution is now set to %d. Model fitting should now be rerun.\n\n", modelres);
	    validentry = 1;
	  }
	}
      }

      break;



    case 52: // Make a spectral age map for a given dataset

      if (export == 1) {

	// Open the directory that is currently set and check it exists
	dir = opendir(imageloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	  break;
	}
      }

      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to map the spectral age (-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (jpmodelset[currentset] != 1) && (kpmodelset[currentset] != 1) && (jptribmodelset[currentset] != 1) && (kptribmodelset[currentset] != 1)  ) { // Check the curvature has been set
	    fprintf(stderr,"\nError: A model has not yet been fitted to this data set! Please run a model fitting command first.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	// Ask which model to map and check its validity

	validentry = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Which model would you like to map? (-1 to escape, 1 for JP model, 2 for KP model, 3 Tribble (JP), 4 Tribble (KP)): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  model = strtol(cmdbuffer, &endptr, 10);

	  if ( (model == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (model < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (model > NUMBEROFMODELS) ||  (model == 0) ) { // Check the set exists
	    printf("\nInvalid model number. Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (model == 1) && (jpmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: JP model values have not yet been calculated for this data set! Please run fitjpmodel first.\n\n");
	    continue;
	  }

	  else if ( (currentset != 888) && (model == 2) && (kpmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: KP model values have not yet been calculated for this data set! Please run fitkpmodel first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (model == 3) && (jptribmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Tribble (JP) model values have not yet been calculated for this data set! Please run fitjptribble first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (model == 4) && (kptribmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Tribble (KP) model values have not yet been calculated for this data set! Please run fitkptribble first.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}
	
	if ( (contours == 1) && (currentset != 888) ) {

	  if ( (exportfits != 1)  || (export != 1) ) {

	    validentry = 0;

	    printf("============================================\n\n");
	    printf(" Index  |  Frequency (Obs.)  |  Frequency (Rest)  \n\n");

	    // Loop through each of the sets
	    for (a=0; a<imgnum[currentset]; a++) {

	      printf("   %d       %.2e Hz             %.2e Hz\n", fluxorder[currentset][a], frequencyobs[currentset][fluxorder[currentset][a]], frequency[currentset][fluxorder[currentset][a]]);

	    }

	    printf("\n============================================\n");
    

    

	    while (validentry == 0) {

	      cmdbuffer = readline("Select a map frequency for the contours (-1 to escape): ");
	      if (cmdbuffer[0] != 0) {
		add_history(cmdbuffer);
	      }

	      chomp(cmdbuffer);

	      currentmap = strtol(cmdbuffer, &endptr, 10);

	      if ( (currentmap == 0) && (cmdbuffer == endptr) ) {
		printf("\nInvalid input, please try again...\n\n");
		continue;
	      }

	      if (currentmap < 0) {
		printf("Escaping command...\n");
		escape = 1;
		break;
	      }   
	      else if (currentmap >= imgnum[currentset]) { // Check the set exists
		printf("Invalid selection, please try again...\n");
		continue;
	      }
	      else {
		validentry = 1;
	      }
	    }


	  }
	}
	
	else if ( (contours == 1) && (currentset == 888) ) { // Multimaps always uses map 0 as reference. This can easily be changed if needed
	  currentmap = 0;
	}
	else { //If we are not using contours, make sure its still a sensible value
	  currentmap = 0;
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if ( (currentset == 888) && (export == 0) ) {
	  looplimit = numdatasets;
	  extwincontrol = 1;
	  sprintf(output,"/xs");
	  cpgopen(output);
	}
	else if ( (currentset == 888) && (export == 1) ) {
	  looplimit = numdatasets;
	  extwincontrol = 0;
	}
	else {
	  looplimit = 1;
	  extwincontrol = 0;
	}
      

	for (j=0; j<looplimit; j++) {

	  if (looplimit > 1) {
	    currentset = j;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = j;
	  }

	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets
	  if (regionsset[currentset] != 1) {
	    printf("\nThe regions for data set %d have not been set! Skipping this set...\n", j);
	    continue;	  
	  }

	  if ( (model == 1) && (jpmodelset[currentset] != 1) ) {
	    printf("\nThe JP model values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  else if ( (model == 2) && (kpmodelset[currentset] != 1) ) {
	    printf("\nThe KP model values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  else if ( (model == 3) && (jptribmodelset[currentset] != 1) ) {
	    printf("\nThe Tribble (JP) model values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  else if ( (model == 4) && (kptribmodelset[currentset] != 1) ) {
	    printf("\nThe Tribble (KP) model values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  time( &currenttime );
	  time_struct = localtime ( &currenttime );
	    
	  strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	  if ( ( export == 1 ) && (exportfits != 1) ) {
	    if (model == 1) {
	      sprintf(modelbuff,"JP");
	    }
	    else if (model == 2) {
	      sprintf(modelbuff,"KP");
	    }
	    else if (model == 3) {
	      sprintf(modelbuff,"TribbleJP");
	    }
	    else if (model == 4) {
	      sprintf(modelbuff,"TribbleKP");
	    }
	    else {
	      sprintf(modelbuff,"Unknown");
	    }


	    if (contours == 1) {
	      sprintf(output,"%s/%s_SpectralAgeing_Contours_%s_%s.%s/%s", imageloc, settarget[currentset], modelbuff, timebuff, imagetype, imagetype);
	      }
	      else {
		sprintf(output,"%s/%s_SpectralAgeing_%s_%s.%s/%s", imageloc, settarget[currentset], modelbuff, timebuff, imagetype, imagetype);
	      }

	    printf("\nExporting as %s\n", output);
	  }

	  else if ( ( export == 1 ) && (exportfits == 1) ) {
	    if (model == 1) {
	      sprintf(output,"%s/%s_SpectralAgeing_JP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else if (model == 2) {
	      sprintf(output,"%s/%s_SpectralAgeing_KP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else if (model == 3) {
	      sprintf(output,"%s/%s_SpectralAgeing_TribbleJP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else if (model == 4) {
	      sprintf(output,"%s/%s_SpectralAgeing_TribbleKP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else {
	      sprintf(output,"%s/%s_SpectralAgeing_UnknownModel_%s.fits", imageloc, settarget[currentset], timebuff);
	    }

	    printf("\nExporting as %s\n", output);
	  }
	  else {
	    sprintf(output,"/xs");
	  }
	
	  if (model == 1) {
      
	    spectralagemap(jpbestage[currentset], regionarray[currentset], flux[currentset][currentmap], xdim[currentset], ydim[currentset], settarget[currentset], border[currentset], zoom, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, firstcontour, contourlogs, output, alwayszeroage, jpinject[currentset], jpfield[currentset], paramlabels, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, xshift, yshift, export, exportfits, largetxt);
	  }
	  else if (model == 2) {
	    spectralagemap(kpbestage[currentset], regionarray[currentset], flux[currentset][currentmap], xdim[currentset], ydim[currentset], settarget[currentset], border[currentset], zoom, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, firstcontour, contourlogs, output, alwayszeroage, kpinject[currentset], kpfield[currentset], paramlabels, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, xshift, yshift, export, exportfits, largetxt);
	  }
	  else if (model == 3) {
      
	    spectralagemap(jptribbestage[currentset], regionarray[currentset], flux[currentset][currentmap], xdim[currentset], ydim[currentset], settarget[currentset], border[currentset], zoom, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, firstcontour, contourlogs, output, alwayszeroage, jptribinject[currentset], jptribfield[currentset], paramlabels, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, xshift, yshift, export, exportfits, largetxt);
	  }
	  else if (model == 4) {
      
	    spectralagemap(kptribbestage[currentset], regionarray[currentset], flux[currentset][currentmap], xdim[currentset], ydim[currentset], settarget[currentset], border[currentset], zoom, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, firstcontour, contourlogs, output, alwayszeroage, kptribinject[currentset], kptribfield[currentset], paramlabels, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, xshift, yshift, export, exportfits, largetxt);
	  }
	  else {
	    fprintf(stderr,"\nError: Model number is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");
	    break;
	  }

	}

	if (extwincontrol == 1) {
	  cpgclos();
	  extwincontrol = 0;
	}

      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;


    case 53: // Make a chi-squared map for a given dataset

      
      if (export == 1) {

	// Open the directory that is currently set and check it exists
	dir = opendir(imageloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	  break;
	}
      }

      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to map the chi-squared values (-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (jpmodelset[currentset] != 1) && (kpmodelset[currentset] != 1) && (jptribmodelset[currentset] != 1) && (kptribmodelset[currentset] != 1) ) { // Check the curvature has been set
	    fprintf(stderr,"\nError: A model has not yet been fitted to this data set! Please run a model fitting command first.\n\n");
	    continue;
	  }

	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}


	// Ask which model to map and check its validity
	validentry = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Which model would you like to map? (-1 to escape, 1 JP model, 2 KP model, 3 Tribble (JP) model, 4 Tribble (KP)): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  model = strtol(cmdbuffer, &endptr, 10);

	  if ( (model == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (model < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (model > NUMBEROFMODELS) ||  (model == 0) ) { // Check the set exists
	    printf("\nInvalid model number. Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (model == 1) && (jpmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: JP model vales have not yet been calculated for this data set! Please run fitjpmodel first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (model == 2) && (kpmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: KP model vales have not yet been calculated for this data set! Please run fitkpmodel first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (model == 3) && (jptribmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Tribble (JP) model vales have not yet been calculated for this data set! Please run fitjptribble first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (model == 4) && (kptribmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Tribble (KP) model vales have not yet been calculated for this data set! Please run fitkptribble first.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}


	if ( (contours == 1) && (currentset != 888) ) {

	  if ( (exportfits != 1)  || (export != 1) ) {

	    validentry = 0;

	    printf("============================================\n\n");
	    printf(" Index  |  Frequency (Obs.)  |  Frequency (Rest)  \n\n");

	    // Loop through each of the sets
	    for (a=0; a<imgnum[currentset]; a++) {

	      printf("   %d       %.2e Hz             %.2e Hz\n", fluxorder[currentset][a], frequencyobs[currentset][fluxorder[currentset][a]], frequency[currentset][fluxorder[currentset][a]]);

	    }

	    printf("\n============================================\n");
    

	    while (validentry == 0) {

	      cmdbuffer = readline("Select a maps frequency for the contours (-1 to escape): ");
	      if (cmdbuffer[0] != 0) {
		add_history(cmdbuffer);
	      }

	      chomp(cmdbuffer);

	      currentmap = strtol(cmdbuffer, &endptr, 10);

	      if ( (currentmap == 0) && (cmdbuffer == endptr) ) {
		printf("\nInvalid input, please try again...\n\n");
		continue;
	      }

	      if (currentmap < 0) {
		printf("Escaping command...\n");
		escape = 1;
		break;
	      }
	      else if (currentmap >= imgnum[currentset]) { // Check the set exists
		printf("\nInvalid selection, please try again...\n\n");
		continue;
	      }
	      else {
		validentry = 1;
	      }
	    }
	  }
	}
	else if ( (contours == 1) && (currentset == 888) ) { // Multimaps always uses map 0 as reference. This can easily be changed if needed
	  currentmap = 0;
	}
	else {
	  currentmap = 0;
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if ( (currentset == 888) && (export == 0) ) {
	  looplimit = numdatasets;
	  extwincontrol = 1;
	  cpgopen("/xs");
	}
	else if ( (currentset == 888) && (export == 1) ) {
	  looplimit = numdatasets;
	  extwincontrol = 0;
	}
	else {
	  looplimit = 1;
	  extwincontrol = 0;
	}

	for (j=0; j<looplimit; j++){

	  if (looplimit > 1) {
	    currentset = j;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = j;
	  }

	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets

	  if (regionsset[currentset] != 1) {
	    printf("\nThe regions for data set %d have not been set! Skipping this set...\n", j);
	    continue;	  
	  }

	  if ( (model == 1) && (jpmodelset[currentset] != 1) ) {
	    printf("\nThe JP model values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  else if ( (model == 2) && (kpmodelset[currentset] != 1) ) {
	    printf("\nThe KP model values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  else if ( (model == 3) && (jptribmodelset[currentset] != 1) ) {
	    printf("\nThe Tribble (JP) model values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  else if ( (model == 4) && (kptribmodelset[currentset] != 1) ) {
	    printf("\nThe Tribble (JP) model values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  dof = (imgnum[currentset] - 2);

	  if (suppresscdf !=1) {
	    // print out the confidence levels
	    interval = 0.68;
	    siglevel68 = gsl_cdf_chisq_Pinv(interval, dof);
	    printf("\nModel can be rejected with 68 per cent confidence at X^2 > %.2f (Reduced: %.2f)\n", siglevel68, siglevel68 / dof);
	    interval = 0.90;
	    siglevel90 = gsl_cdf_chisq_Pinv(interval, dof);
	    printf("Model can be rejected with 90 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel90, siglevel90 / dof);
	    interval = 0.95;
	    siglevel95 = gsl_cdf_chisq_Pinv(interval, dof);
	    printf("Model can be rejected with 95 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel95, siglevel95 / dof);
	    interval = 0.99;
	    siglevel99 = gsl_cdf_chisq_Pinv(interval, dof);
	    printf("Model can be rejected with 99 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel99, siglevel99 / dof);
	    interval = 0.995;
	    siglevel995 = gsl_cdf_chisq_Pinv(interval, dof);
	    printf("Model can be rejected with 99.5 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel995, siglevel995 / dof);
	    interval = 0.999;
	    siglevel999 = gsl_cdf_chisq_Pinv(interval, dof);
	    printf("Model can be rejected with 99.9 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel999, siglevel999 / dof);
	    interval = 0.9999;
	    siglevel9999 = gsl_cdf_chisq_Pinv(interval, dof);
	    printf("Model can be rejected with 99.99 per cent confidence X^2 > %.2f (Reduced: %.2f)\n\n", siglevel9999, siglevel9999 / dof);
	  }

	    time( &currenttime );
	    time_struct = localtime ( &currenttime );

	    strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	  if ( ( export == 1 ) && (exportfits != 1) ) {
	    	    
	    if (model == 1) {
	      sprintf(modelbuff,"JP");
	    }
	    else if (model == 2) {
	      sprintf(modelbuff,"KP");
	    }
	    else if (model == 3) {
	      sprintf(modelbuff,"TribbleJP");
	    }
	    else if (model == 4) {
	      sprintf(modelbuff,"TribbleKP");
	    }
	    else {
	      sprintf(modelbuff,"Unknown");
	    }


	    if (contours == 1) {
	      sprintf(output,"%s/%s_ChiSquared_Contours_%s_%s.%s/%s", imageloc, settarget[currentset], modelbuff, timebuff, imagetype, imagetype);
	    }
	    else {
	      sprintf(output,"%s/%s_ChiSquared_%s_%s.%s/%s", imageloc, settarget[currentset], modelbuff, timebuff, imagetype, imagetype);
	    }

	    printf("\nExporting as %s\n", output);
	  }

	  else if ( ( export == 1 ) && (exportfits == 1) ) {

	    if (model == 1) {
	      sprintf(output,"%s/%s_ChiSquared_JP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else if (model == 2) {
	      sprintf(output,"%s/%s_ChiSquared_KP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else if (model == 3) {
	      sprintf(output,"%s/%s_ChiSquared_TribbleJP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else if (model == 4) {
	      sprintf(output,"%s/%s_ChiSquared_TribbleKP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else {
	      sprintf(output,"%s/%s_ChiSquared_UnknownModel_%s.fits", imageloc, settarget[currentset], timebuff);
	    }

	    printf("\nExporting as %s\n", output);
	  }
	  else {
	    sprintf(output,"/xs");
	  }

	  if (usereduced == 1) {
	    sprintf(maptitle,"Map of Reduced Chi-Squared Values for %s", settarget[currentset]);
	  }
	  else {
	    sprintf(maptitle,"Map of Chi-Squared Values for %s", settarget[currentset]);
	  }
	

	  if (model == 1) {
	    chisquaredmap(jpchisquared[currentset], regionarray[currentset], xdim[currentset], ydim[currentset], settarget[currentset], imgnum[currentset], border[currentset], zoom, maptitle, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, flux[currentset][currentmap], firstcontour, contourlogs, output, jpinject[currentset], jpfield[currentset], paramlabels, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, usereduced, xshift, yshift, export, exportfits, largetxt);
	  }
	  else if (model == 2) {
	    chisquaredmap(kpchisquared[currentset], regionarray[currentset], xdim[currentset], ydim[currentset], settarget[currentset], imgnum[currentset], border[currentset], zoom, maptitle, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, flux[currentset][currentmap], firstcontour, contourlogs, output, kpinject[currentset], kpfield[currentset], paramlabels, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, usereduced, xshift, yshift, export, exportfits, largetxt);
	  }
	  else if (model == 3) {
	    chisquaredmap(jptribchisquared[currentset], regionarray[currentset], xdim[currentset], ydim[currentset], settarget[currentset], imgnum[currentset], border[currentset], zoom, maptitle, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, flux[currentset][currentmap], firstcontour, contourlogs, output, jptribinject[currentset], jptribfield[currentset], paramlabels, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, usereduced, xshift, yshift, export, exportfits, largetxt);
	  }
	  else if (model == 4) {
	    chisquaredmap(kptribchisquared[currentset], regionarray[currentset], xdim[currentset], ydim[currentset], settarget[currentset], imgnum[currentset], border[currentset], zoom, maptitle, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, flux[currentset][currentmap], firstcontour, contourlogs, output, kptribinject[currentset], kptribfield[currentset], paramlabels, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, usereduced, xshift, yshift, export, exportfits, largetxt);
	  }
	  else {
	    fprintf(stderr,"\nError: Model number is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");
	    break;
	  }

	}

	if (extwincontrol == 1) {
	  cpgclos();
	  extwincontrol = 0;
	}

      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;


    case 54: // Change between CASA and AIPS format DATA

      if (casadata == 1) {
	casadata = 0;
	printf("BRATS will now default to AIPS formatted headers\n");
      }
      else if (casadata == 0) {
	casadata = 1;
	printf("BRATS will now default to CASA formatted headers\n");
      }
      else {
	printf("\nThe value of casadata is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", casadata);
      }

      break;


    case 55: // Set number of (age) levels deep to go when model fitting

      validentry = 0;
      tmp_levels = levels;

      printf("Levels currently set to %d (DEFAULT: %d)\n", levels, DEFAULTLEVELS);

      while (validentry == 0) {

	cmdbuffer = readline("Enter the number of levels to apply to model fitting (Integers > 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	levels = strtol(cmdbuffer, &endptr, 10);

	if ( (levels == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else if (levels == 0) {
	  printf("Value must be greater than or equal to 1. Please try again.\n");
	  continue;
	}
	else {
	  if (levels < 0) {
	    printf("Escaping command...\n");
	    levels = tmp_levels;
	    break;
	  }
	  else {
	    printf("\nNumber of levels to apply to model fitting now set to %d\n\n", levels);
	    validentry = 1;
	  }
	}
      }

      break;


    case 56: // Turn on and off the printing of results during model fitting

      if (printresults == 1) {
	printresults = 0;
	printf("Printing of individual results whilst model fitting is now turned off!\n");
      }
      else if (printresults == 0) {
	printresults = 1;
	printf("Printing of individual results whilst model fitting is now turned on!\n");
      }
      else {
	printf("\nThe value of printresults is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", printresults);
      }

      break;

    case 57: // Turn on and off smooth maps

      if (smoothmaps == 1) {
	smoothmaps = 0;
	printf("Plotting of smooth maps is now turned off!\n");
      }
      else if (smoothmaps == 0) {
	smoothmaps = 1;
	printf("Plotting of smooth maps is now turned on!\n");
      }
      else {
	printf("\nThe value of smoothmaps is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", smoothmaps);
      }

      break;


    case 58: // Turn on and off smooth maps

      if (contours == 1) {
	contours = 0;
	printf("Plotting of contours is now turned off!\n");
      }
      else if (contours == 0) {
	contours = 1;
	printf("Plotting of contours is now turned on!\n");
      }
      else {
	printf("\nThe value of contours is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", contours);
      }

      break;



    case 59: // Set number of contour levels to apply to maps

      validentry = 0;
      tmp_contourlevels = contourlevels;

      printf("Number of contours levels currently set to %d (DEFAULT: %d)\n", contourlevels, DEFAULTCONLEVELS);

      while (validentry == 0) {

	cmdbuffer = readline("Enter the number of contour levels to apply (integers > 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	contourlevels = strtol(cmdbuffer, &endptr, 10);

	if ( (contourlevels == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else if (contourlevels == 0) {
	  printf("\nValue must be greater than or equal to 1. Please try again.\n\n");
	  //contourlevels = tmp_contourlevels;
	  continue;
	}
	else {

	  if (contourlevels < 0) {
	    printf("Escaping command...\n");
	    contourlevels = tmp_contourlevels;
	    break;
	  }
	  else {
	    printf("\nNumber of contour levels to apply now set to %d\n\n", contourlevels);
	    validentry = 1;
	  }
	}
      }

      break;


    case 60: // Set colour of the contour lines

      validentry = 0;
      tmp_contourcolour = contourcolour;

      printf("Contour line colour currently set to %d (DEFAULT: %d)\n", contourcolour, DEFAULTCONCOLOUR);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the colour of the contour lines (integers 0 - 15, -1 to escape, 999 for a list of colours): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	contourcolour= strtol(cmdbuffer, &endptr, 10);

	if ( (contourcolour == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else if ( (contourcolour > 15) && (contourcolour != 999) ) {
	  printf("\nValue must be between 0 and 15 (inclusive). Please try again.\n\n");
	  //contourcolour = tmp_contourcolour;
	  continue;
	}
	else if (contourcolour == 999) {
	  printf("\n0: White, 1: Black, 2: Red, 3: Green, 4: Blue, 5: Light blue, 6: Pink, 7: Yellow, 8: Orange, 9: Lime green, 10: Pale green, 11: Pale blue, 12: Purple, 13: Dark pink, 14: Brown grey, 15: Light grey\n\n");
	  //contourcolour = tmp_contourcolour;
	  continue;
	}
	else {
	  if (contourcolour < 0) {
	    printf("Escaping command...\n");
	    contourcolour = tmp_contourcolour;
	    break;
	  }
	  else {
	    printf("\nThe colour of the contour lines is now set to %d\n\n", contourcolour);
	    validentry = 1;
	  }
	}
      }

      break;


    case 61: // Set colour of the contour lines

      validentry = 0;
      tmp_contourlinestyle = contourlinestyle;

      printf("Contour line style currently set to %d (DEFAULT: %d)\n", contourlinestyle, DEFAULTCONLINESTYLE);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the contour line style (integers 1 - 5, -1 to escape, 999 for a list of styles): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	contourlinestyle = strtol(cmdbuffer, &endptr, 10);

	if ( (contourlinestyle == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else if ( ( (contourlinestyle == 0) || (contourlinestyle > 5) ) && (contourlinestyle != 999) ) {
	  printf("\nValue must be between 1 and 5 (inclusive). Please try again.\n\n");
	  //contourlinestyle = tmp_contourlinestyle;
	  continue;
	}
	else if (contourlinestyle == 999) {
	  printf("\n1: Solid, 2: Dashed, 3: Dot-dash-dot-dash, 4: Dotted, 5: Dash-dot-dot-dot\n\n");
	  //contourlinestyle = tmp_contourlinestyle;
	  continue;
	}
	else {
	  if (contourlinestyle < 0) {
	    printf("Escaping command...\n");
	    contourlinestyle = tmp_contourlinestyle;
	    break;
	  }
	  else {
	    printf("\nThe contour line style is now set to %d\n\n", contourlinestyle);
	    validentry = 1;
	  }
	}
      }
    
      break;


    case 62: // Switch between screen display and exporting maps

      if (export == 1) {
	export = 0;
	printf("Exporting of plots and maps is now turned off! Output will be to screen\n");
      }
      else if (export == 0) {
	export = 1;
	printf("Exporting of plots and maps is now turned on! Output type will be that chosen using the device command (CURRENT: %s DEFAULT: %s)\n", imagetype, DEFAULTIMAGETYPE);
      }
      else {
	printf("\nThe value of export is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", export);
      }

      break;


      // Skip case 63 as this relates to '?'


    case 64: // Set where exported data will be saved

      validentry = 0;
      strcpy(tmp_imageloc, imageloc);

      printf("Currently exporting images to %s (DEFAULT: %s)\n", imageloc, DEFAULTIMAGELOC);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a folder location to store exported images. Absolute or relative paths allowed. No training slash (esc to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	strcpy(imageloc, cmdbuffer);

	if ((strstr(imageloc, "esc") != NULL) && (strlen(imageloc) == 3) ) {
	  strcpy(imageloc, tmp_imageloc);
	  printf("Escaping command...\n");
	  break;
	}
	else {

	  // Open the directory that us currently set and check it exists
	  dir = opendir(imageloc);

	  if (dir != NULL) {
	    printf("\nExported data will now be saved to %s\n\n", imageloc);
	    validentry = 1;
	    closedir (dir);
	  }
	  else {
	    strcpy(imageloc, tmp_imageloc);
	    printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists and the correct permission are set. ***\n\n", imageloc);
	    continue;
	  }
	}
      }
    
      break;


    case 65: // Set the export image format type

      validentry = 0;
      strcpy(tmp_imagetype, imagetype);

      printf("Image export format currently set to %s (DEFAULT: %s)\n", imagetype, DEFAULTIMAGETYPE);


      while (validentry == 0) {

	cmdbuffer = readline("Enter an image export format type e.g. png. For a full list option type 'all'. DO NOT INCLUDE THE LEADING SLASH. Warning: Some option may cause the program to crash if they have not been setup correctly! (String, esc to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	for (tmp_char = cmdbuffer; *tmp_char != '\0'; ++tmp_char) {
	  *tmp_char = toupper(*tmp_char);
	}

	strcpy(imagetype, cmdbuffer);

	if ( (strcmp(imagetype,"GIF") != 0) && (strcmp(imagetype,"VGIF") != 0) && (strcmp(imagetype,"NULL") != 0) && (strcmp(imagetype,"PNG") != 0) && (strcmp(imagetype,"TPNG") != 0) && (strcmp(imagetype,"PS") != 0) && (strcmp(imagetype,"VPS") != 0) && (strcmp(imagetype,"VCPS") != 0) && (strcmp(imagetype,"CPS") != 0) && (strcmp(imagetype,"ALL") != 0) && (strcmp(imagetype,"ESC") != 0) ) {
	  printf("\nInvalid format type, please try again...\n\n");
	  continue;
	}

	if ((strstr(imagetype, "ESC") != NULL) && (strlen(imagetype) == 3) ) { // Caps here as we have to shift for the iamge type
	  strcpy(imagetype, tmp_imagetype);
	  printf("Escaping command...\n");
	  break;
	}

	else if ((strstr(imagetype, "ALL") != NULL) && (strlen(imagetype) == 3) ) {
	  strcpy(imagetype, tmp_imagetype);

	  printf("    GIF       (Graphics Interchange Format file, landscape orientation)\n");
	  printf("    VGIF      (Graphics Interchange Format file, portrait orientation)\n");
	  printf("    NULL      (Null device, no output)\n");
	  printf("    PNG       (Portable Network Graphics file)\n");
	  printf("    TPNG      (Portable Network Graphics file - transparent background)\n");
	  printf("    PS        (PostScript file, landscape orientation)\n");
	  printf("    VPS       (PostScript file, portrait orientation)\n");
	  printf("    CPS       (Colour PostScript file, landscape orientation)\n");
	  printf("    VCPS      (Colour PostScript file, portrait orientation)\n");
	  continue;
	}
	else {
	  printf("\nImages will now be exported in %s format\n\n", imagetype);
	  validentry = 1;
	}
      }
    
      break;


    case 66: // Set number of the lowest contour levels to apply to maps

      validentry = 0;
      tmp_firstcontour = firstcontour;

      printf("The first contour levelt o be plotted is currently set to %d (DEFAULT: %d)\n", firstcontour, DEFAULTFIRSTCONTOUR);


      while (validentry == 0) {

	cmdbuffer = readline("Enter the number of the first (lowest flux) contour to plot (Integers >= 0 and < contourlevels-1, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	firstcontour = strtol(cmdbuffer, &endptr, 10);

	if ( (firstcontour == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else if (firstcontour > contourlevels-1) {
	  printf("\nValue must be greater than contourlevels - 1. Please try again or change the value of contourlevels\n\n");
	  continue;
	}
	else {
	  if (firstcontour < 0) {
	    printf("Escaping command...\n");
	    firstcontour = tmp_firstcontour;
	    break;
	  }
	  else {
	    printf("\nFirst contour now set to %d\n\n", firstcontour);
	    validentry = 1;
	  }
	}
      }

      break;


    case 67: // Turn on and off use of log spacing for contours

      if (contourlogs == 1) {
	contourlogs = 0;
	printf("Contour will now be linearly spaced\n");
      }
      else if (contourlogs == 0) {
	contourlogs = 1;
	printf("Contours will now be logarithmically spaced\n");
      }
      else {
	printf("\nThe value of contourlogs is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", contourlogs);
      }

      break;



    case 68: // Turn on and off use of log spacing for flux maps

      if (fluxlogs == 1) {
	fluxlogs = 0;
	printf("Flux will now be mapped linearly (raw data only!)\n");
      }
      else if (fluxlogs == 0) {
	fluxlogs = 1;
	printf("Flux will now be mapped logarithmically (raw data only!)\n");
      }
      else {
	printf("\nThe value of fluxlogs is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", fluxlogs);
      }

      break;


    case 69: // Turn on and off use of log spacing for fluxes

      if (fluxcut == 1) {
	fluxcut = 0;
	printf("Minimum mapped flux will now be the detection value set by 'sigma' (raw data only!)\n");
      }
      else if (fluxcut == 0) {
	fluxcut = 1;
	printf("Minimum mapped flux will now be cut to the user specified value set by 'usercut' (raw data only!)\n");
      }
      else {
	printf("\nThe value of fluxcut is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", fluxcut);
      }

      break;


    case 70: // Set the cut value for flux maps

      validentry = 0;
      tmp_usercut = usercut;

      printf("Cut value for flux maps currently set to %.2f (DEFAULT: %.2f)\n", usercut, DEFAULTUSERCUT);

      while (validentry == 0) {

	cmdbuffer = readline("Enter the desired cut value for flux maps. This must also be turned on using the fluxcut command (raw data only!) (Float >= 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	usercut = strtod(cmdbuffer, &endptr);

	if ( (usercut == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (usercut < 0) {
	    printf("Escaping command...\n");
	    usercut = tmp_usercut;
	    break;
	  }
	  else {
	    printf("\nCut now set to %.2e\n\n", usercut);
	    validentry = 1;
	  }
	}
      }

      break;


    case 71: // Set the wedge multiplier to give the best colour range on maps

      validentry = 0;
      tmp_wedgemultiplier = wedgemultiplier;

      printf("Multiplier to apply to the maximum wedge range currently set to %.2f (DEFAULT: %.2f)\n", wedgemultiplier, DEFAULTWEDGE);


      while (validentry == 0) {

	cmdbuffer = readline("Enter the multiplier to apply to the maximum wedge range (Float > 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	wedgemultiplier = strtod(cmdbuffer, &endptr);

	if ( (wedgemultiplier == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (wedgemultiplier < 0) {
	    printf("Escaping command...\n");
	    wedgemultiplier = tmp_wedgemultiplier;
	    break;
	  }
	  else {
	    printf("\nThe wedge multiplier is now set to %.2f\n\n", wedgemultiplier);
	    validentry = 1;
	  }
	}
      }

      break;


    case 72: // Set the magnetic field strength for models

      validentry = 0;
      tmp_fieldstrength = fieldstrength;
      printf("Magnetic field strength currently set to %.2e T (DEFAULT: %.2e T)\n",fieldstrength , DEFAULTFIELDSTRENGTH);


      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the magnetic field strength in Tesla (Double > 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	fieldstrength = strtod(cmdbuffer, &endptr);

	if ( (fieldstrength == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (fieldstrength < 0) {
	    printf("Escaping command...\n");
	    fieldstrength = tmp_fieldstrength;
	    break;
	  }
	  else {
	    printf("\nThe magnetic field strength is now set to %.2e T\n\n", fieldstrength);
	    validentry = 1;
	  }
	}
      }

      break;


    case 73: // Set the injection index for models

      validentry = 0;
      tmp_dbl_inject = inject;

      printf("Injection index currently set to %.2f (DEFAULT: %.2f)\n", inject, DEFAULTINJECT);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the injection index (Double > 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	inject = strtod(cmdbuffer, &endptr);

	if ( (inject == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (inject < 0) {
	    printf("Escaping command...\n");
	    inject = tmp_dbl_inject;
	    break;
	  }
	  else {
	    power = (2*inject) + 1; // Update the power equivilant
	    printf("\nThe model injection index is now set to %.2f (model power of %.1f)\n\n", inject, power);
	    validentry = 1;
	  }
	}
      }

      break;


    case 74: // Create a colour colour plot of a data set


      if (numdatasets > 0) {

	// Open the directory that is currently set and check it exists
	dir = opendir(imageloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	  break;
	}
      
      
	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to create the colour-colour plot (-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    printf("\nThis value is greater than the number of data sets that exists! Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    //escape = 1;
	    //break;
	    continue;
	  }

	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}


	if ( (currentset == 888) && (export == 0) ) {
	  looplimit = numdatasets;
	  extwincontrol = 1;
	  cpgopen("/xs");
	}
	else if ( (currentset == 888) && (export == 1) ) {
	  looplimit = numdatasets;
	  extwincontrol = 0;
	}
	else {
	  looplimit = 1;
	  extwincontrol = 0;
	}
      

	for (j=0; j<looplimit; j++){

	  if (looplimit > 1) {
	    currentset = j;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = j;
	  }

	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets

	  if (regionsset[currentset] != 1) {
	    printf("\nThe regions for data set %d have not been set! Skipping this set...\n", j);
	    continue;	  
	  }


	  if ( export == 1 ) {

	    time( &currenttime );
	    time_struct = localtime ( &currenttime );

	    strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	    sprintf(output,"%s/%s_colourcolourplot_%s.%s/%s", imageloc, settarget[currentset], timebuff, imagetype, imagetype);

	    printf("\nExporting as %s\n", output);

	  }
	  else {
	    sprintf(output,"/xs");
	  }


	  // Get the maps to use

	  validentry = 0;


	  printf("============================================\n\n");
	  printf(" Index  |  Frequency (Obs.)  |  Frequency (Rest)  \n\n");

	  // Loop through each of the sets
	  for (a=0; a<imgnum[currentset]; a++) {

	    printf("   %d       %.2e Hz             %.2e Hz\n", fluxorder[currentset][a], frequencyobs[currentset][fluxorder[currentset][a]], frequency[currentset][fluxorder[currentset][a]]);

	  }

	  printf("\n============================================\n");
    
	
	  while (validentry == 0) {

	    cmdbuffer = readline("Select the lowest frequency map to be used (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    map1 = strtol(cmdbuffer, &endptr, 10);

	    if ( (map1 == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (map1 < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else if (map1 >= imgnum[currentset]) { // Check the set exists
	      printf("\nInvalid selection, please choose again...\n\n");
	      continue;
	    }
	    else {
	      validentry = 1;
	    }

	  }	

	  if (escape == 1) {
	    break;
	  }

      
      
	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Select the mid frequency map to be used (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    map2 = strtol(cmdbuffer, &endptr, 10);

	    if ( (map2 == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (map2 < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	
	    else if (map2 >= imgnum[currentset]) { // Check the set exists
	      printf("\nInvalid selection, please choose again...\n\n");
	      continue;
	    }
	    else {
	      validentry = 1;
	    }

	  }

	  if (escape == 1) {
	    break;
	  }

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Select the highest frequency map to be used (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    map3 = strtol(cmdbuffer, &endptr, 10);

	    if ( (map3 == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (map3 < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	
	    else if (map3 >= imgnum[currentset]) { // Check the set exists
	      printf("\nInvalid selection, please choose again...\n\n");
	      continue;
	    }
	    else {
	      validentry = 1;
	    }

	  }

	  if (escape == 1) {
	    break;
	  }


	  // Setup the default values and memory
	  mincolour1=1e34;
	  maxcolour1=-1e34;
	  mincolour2=1e34;
	  maxcolour2=-1e34;

	  printf("map1 [%d]: %.2e map2 [%d]: %.2e map3 [%d]: %.2e\n", map1, frequency[currentset][map1], map2, frequency[currentset][map2], map3, frequency[currentset][map3]);

	  colour1 = (float *)calloc(regnumber[currentset], sizeof(float));
	  colour2 = (float *)calloc(regnumber[currentset], sizeof(float));

	  // Calculate the spectral index between the first 2 maps
	  calcspecindexp2p(regflux[currentset], colour1, regnumber[currentset], map1, map2, frequency[currentset][map1], frequency[currentset][map2], printindex, &mincolour1, &maxcolour1, rms[currentset]);

	  calcspecindexp2p(regflux[currentset], colour2, regnumber[currentset], map2, map3, frequency[currentset][map2], frequency[currentset][map3], printindex, &mincolour2, &maxcolour2, rms[currentset]);

	  // Plot the points

	  colourcolour(colour1, colour2, mincolour1, maxcolour1, mincolour2, maxcolour2, frequency[currentset][map1], frequency[currentset][map2], frequency[currentset][map3], settarget[currentset], regnumber[currentset], axismarks, titles, labels, symbol, extwincontrol, output, swapcoloursign);

	  //free(colour1);
	  //free(colour2);


	  if (extwincontrol == 1) {
	    cpgclos();
	    extwincontrol = 0;
	  }
	}

	if (extwincontrol == 1) {
	  cpgclos();
	  extwincontrol = 0;
	}

      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;



    case 75: // Set the injection index for models in terms of power

      validentry = 0;
      tmp_power = power;

      printf("Model power currently set to %.2f (DEFAULT: %.2f)\n", power, (2*DEFAULTINJECT) +1);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the model injection index in terms of power [alpha = (power-1)/2] (Double > 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	power  = strtod(cmdbuffer, &endptr);

	if ( (power == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (power < 0) {
	    printf("Escaping command...\n");
	    power = tmp_power;
	    break;
	  }
	  else {
	    inject = (power-1)/2; // Update the power equivilant
	    printf("\nThe model power is now set to %.2f (injection index of %.1f)\n\n", power, inject);
	    validentry = 1;
	  }
	}
      }

      break;


    case 76: // Turn on and off if the minimum age is always zero for spectral ageing maps

      if (alwayszeroage == 1) {
	alwayszeroage = 0;
	printf("The minimum age range value for spectral ageing maps will now be automatically determined\n");
      }
      else if (alwayszeroage == 0) {
	alwayszeroage = 1;
	printf("The minimum age range value  will now always be zero for spectral ageing maps\n");
      }
      else {
	printf("\nThe value of alwayszeroage is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", alwayszeroage);
      }

      break;


    case 77: // Plot the KP model for an arbitary normalisation

      printf("\nPlotting an example KP model between %.2e Hz and %.2e Hz\n\n", minmodelfreq, maxmodelfreq);


      if (export == 1) {

	// Open the directory that is currently set and check it exists
	dir = opendir(imageloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	  break;
	}

	time( &currenttime );
	time_struct = localtime ( &currenttime );

	strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	sprintf(output,"%s/ExKPModel_%s.%s/%s", imageloc,  timebuff, imagetype, imagetype);

	printf("\nExporting as %s\n", output);

	cpgopen(output);
      }
      else {
	sprintf(output,"/xs");
      }


      model = 2;

      plotmodel(minmodelfreq, maxmodelfreq, inject, fieldstrength, model, usr_gmin, usr_gmax, minmodelmyears, modelmyears, output, titles, exampleredshift, skip);

      break;


    case 78: //Fit the KP model to a dataset
  
      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to fit the KP model(-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    printf("\nThis value is greater than the number of data sets that exists! Please try again...\n\n");
	    continue;
	  }
	  else if ( (imgnum[currentset] <= 1) && (currentset != 888) ) { // Check the set has at least two maps
	    fprintf(stderr,"\nError: The data set must contain a minimum of two images for fitting. Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	// Setup the loops
	if (currentset == 888) {
	  looplimit = numdatasets;
	}
	else {
	  looplimit = 1;
	}

	for (j=0; j<looplimit; j++){

	  if (looplimit > 1) {
	    currentset = j;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = j;
	  }

	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets
	  if (regionsset[currentset] != 1) {
	    printf("\nThe regions for data set %d have not been set! Skipping this set...\n", j);
	    continue;	  
	  }
	  else if (imgnum[currentset] <= 1) { // Check the set has at least two maps
	    fprintf(stderr,"\n*** Warning: Data set %d does not contain the minimum two images required for fitting. Skipping this set... ***\n\n", j);
	    continue;
	  }

	  setmodelres[currentset][2] = modelres; // Set the resolution for this dataset and model

	  // Setup the memory
	  if (kpmodelmemoryset == 0) {

	    kpchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kpbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kpmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kpmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	    kpbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kpageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kpageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));

	    kpmodelmemoryset = 1;
	  }

	  if ( (numdatasets > 1)) {

	    allocarray = (int *)calloc(MAXDATASETS, sizeof(int));

	    realloc_kpchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_kpbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_kpmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_kpmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	    realloc_kpbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_kpageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_kpageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));


	    for (o=0; o<numdatasets; o++) {

	      if (kpfitmemset[o] == 1) {

		realloc_kpchisquared[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_kpbestage[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_kpmodelalpha[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_kpmodelflux[o] = (float **)calloc(regnumber[o]+1, sizeof(float *));
		realloc_kpbestnorm[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_kpageerrorsplus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_kpageerrorsminus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));

		allocarray[o] = 1;

	      }

	    }

	    for (o=0; o<numdatasets; o++) {

	      if (kpfitmemset[o] == 1) {

		for (q=0; q<=regnumber[o]; q++) {

		  //printf("q=%d\n", q);
		  realloc_kpmodelflux[o][q] = (float *)calloc(setmodelres[currentset][2]+1, sizeof(float));
		}
	      }
	    }

	    for (o=0; o<numdatasets; o++) {

	      for (p=0; p<=regnumber[o]; p++) {

		if (kpfitmemset[o] == 1) {

		  realloc_kpchisquared[o][p] = kpchisquared[o][p];
		  realloc_kpbestage[o][p] = kpbestage[o][p];
		  realloc_kpmodelalpha[o][p] = kpmodelalpha[o][p];
		  realloc_kpbestnorm[o][p] = kpbestnorm[o][p];
		  realloc_kpageerrorsplus[o][p] = kpageerrorsplus[o][p];
		  realloc_kpageerrorsminus[o][p] = kpageerrorsminus[o][p];

		  for (q=0; q<=setmodelres[currentset][2]; q++) {
		    realloc_kpmodelflux[o][p][q] = kpmodelflux[o][p][q];
		  }

		}
	      }

	    }

	    for (o=0; o<numdatasets-1; o++) { // -1 as we are yet to allocate the new dataset memory

	      free(kpchisquared[o]);
	      free(kpbestage[o]);
	      free(kpmodelalpha[o]);
	      free(kpbestnorm[o]);
	      free(kpageerrorsplus[o]);
	      free(kpageerrorsminus[o]);
	      free(kpmodelflux[o]);

	    }
	    
	    free(kpchisquared);
	    free(kpbestage);
	    free(kpmodelalpha);
	    free(kpbestnorm);
	    free(kpageerrorsplus);
	    free(kpageerrorsminus);
	    free(kpmodelflux);
	    
	    kpfitmemset[o] = 0;
	    kpmodelmemoryset = 0;


	    // Make the first level of the array again
	    kpchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kpbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kpmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kpmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	    kpbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kpageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kpageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));

	    kpmodelmemoryset = 1;

	  }
	
	  for (o=0; o<numdatasets; o++) {

	    kpchisquared[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    kpbestage[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    kpmodelalpha[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    kpmodelflux[o] = (float **)calloc(regnumber[o]+1, sizeof(float *));
	    kpbestnorm[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    kpageerrorsplus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    kpageerrorsminus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));

	    for (p=0; p<=regnumber[o]; p++) {
	      kpmodelflux[o][p] = (float *)calloc(setmodelres[currentset][2]+1, sizeof(float));
	    }
	    
	    kpfitmemset[o] = 1;

	  }
	

	  if (numdatasets > 1) {
	    for (o=0; o<numdatasets; o++) {

	      if ( (o != currentset) && (allocarray[o] == 1) ) {

		for (p=0; p<=regnumber[o]; p++) {

		  kpchisquared[o][p] = realloc_kpchisquared[o][p];
		  kpbestage[o][p] = realloc_kpbestage[o][p];
		  kpmodelalpha[o][p] = realloc_kpmodelalpha[o][p];
		  kpbestnorm[o][p] = realloc_kpbestnorm[o][p];
		  kpageerrorsplus[o][p] = realloc_kpageerrorsplus[o][p];
		  kpageerrorsminus[o][p] = realloc_kpageerrorsminus[o][p];

		  for (q=0; q<=setmodelres[currentset][2]; q++) {

		    kpmodelflux[o][p][q] = realloc_kpmodelflux[o][p][q];
		  }
		}
	      }
	    }
	  }
	
	  if (numdatasets > 1) {

	    // Free the reallocation arrays if they are set

	    for (o=0; o<numdatasets-1; o++) {

	      free(realloc_kpchisquared[o]);
	      free(realloc_kpbestage[o]);
	      free(realloc_kpmodelalpha[o]);
	      free(realloc_kpbestnorm[o]);
	      free(realloc_kpageerrorsplus[o]);
	      free(realloc_kpageerrorsminus[o]);
	      //free(realloc_kpmodelflux[o]);

	    }

	    free(realloc_kpchisquared);
	    free(realloc_kpbestage);
	    free(realloc_kpmodelalpha);
	    free(realloc_kpbestnorm);
	    free(realloc_kpageerrorsplus);
	    free(realloc_kpageerrorsminus);
	    free(realloc_kpmodelflux);

	    free(allocarray);

	  }

	  model = 2;

	  // Using a temp array here to allow all models without the need for dummy arrays
	  tmp_ageerrorsplus = (float *)calloc(regnumber[currentset]+1, sizeof(float));
	  tmp_ageerrorsminus = (float *)calloc(regnumber[currentset]+1, sizeof(float));


	  spectralageingmodels(regflux[currentset], imgnum[currentset], fluxorder[currentset], frequency[currentset], regnumber[currentset], kpchisquared[currentset], kpbestage[currentset], fluxerror[currentset], kpmodelalpha[currentset], kpmodelflux[currentset], usr_gmin, usr_gmax, minmyears, myears, ageresolution, levels, kpbestnorm[currentset], setmodelres[currentset][2], printresults, inject, fieldstrength, model, printreject, redshift[currentset], tmp_ageerrorsplus, tmp_ageerrorsminus, suppresscdf);


	  kpfield[currentset] = fieldstrength;
	  kpinject[currentset] = inject;
	  kpmodelset[currentset] = 1;

	  // Transfer the age errors over to the perminant array
	  for (r=1; r<=regnumber[currentset]; r++) {
	    kpageerrorsplus[currentset][r] = tmp_ageerrorsplus[r];
	    kpageerrorsminus[currentset][r] = tmp_ageerrorsminus[r];
	  }

	  free(tmp_ageerrorsplus);
 	  free(tmp_ageerrorsminus);  

	}

      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;


    case 79: // Turn on and off if the minimum age is always zero for spectral ageing maps

      if (paramlabels == 1) {
	paramlabels = 0;
	printf("\nParameter labels for model maps and plot are now turned off\n\n");
      }
      else if (paramlabels == 0) {
	paramlabels = 1;
	printf("\nParameter labels for model maps and plot are now turned on\n\n");
      }
      else {
	printf("\nThe value of paramlabels is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", paramlabels);
      }

      break;


    case 80: // Turn on and off if the minimum age is always zero for spectral ageing maps

      if (printreject == 1) {
	printreject = 0;
	printf("\nOutput of whether a model should be rejected based on average chi-squared is now turned off.\n\n");
      }
      else if (printreject == 0) {
	printreject = 1;
	printf("\nOutput of whether a model should be rejected based on average chi-squared is now turned on. Note: If there are a number of regions where you expect a bad fit e.g. from dynamic range issues, the average may not be a suitable value to use.\n\n");
      }
      else {
	printf("\nThe value of printreject is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", printreject);
      }

      break;


    case 81: // Set gmin for model fitting

      validentry = 0;
      tmp_usr_gmin = usr_gmin;

      printf("Minimum gamma value currently set to %.2f (DEFAULT: %.2f)\n", usr_gmin, DEFAULTGMIN);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the minimum gamma value (Double > 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	usr_gmin = strtod(cmdbuffer, &endptr);

	if ( (usr_gmin == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {

	  if (usr_gmin < 0) {
	    printf("Escaping command...\n");
	    usr_gmin = tmp_usr_gmin;
	    break;
	  }
	  else {
	    printf("\nThe minimum gammma is now set to %.2f\n\n", usr_gmin);
	    validentry = 1;
	  }
	}
      }

      break;


    case 82: // Set gmax for model fitting

      validentry = 0;
      tmp_usr_gmax = usr_gmax;

      printf("Maximum gamma value currently set to %.2e (DEFAULT: %.2e)\n", usr_gmax, DEFAULTGMAX);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the maximum gamma value (Double > 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	usr_gmax = strtod(cmdbuffer, &endptr);

	if ( (usr_gmax == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (usr_gmax < 0) {
	    printf("Escaping command...\n");
	    usr_gmax = tmp_usr_gmax;
	    break;
	  }
	  else {
	    printf("\nThe maximum gammma is now set to %.2e\n\n", usr_gmax);
	    validentry = 1;
	  }
	}
      }

      break;


    case 83: // Set which map should be used to define the position when plotting maps in WCS

      validentry = 0;
      tmp_posmap = posmap;

      printf("WCS reference map currently set to %d (DEFAULT: %d)\n", posmap, DEFAULTPOSMAP);

      while (validentry == 0) {

	cmdbuffer = readline("Select the map should to be used as reference for coordinates when plotting maps in WCS. If the number entered is greater than the number of maps in a given dataset, the coordinates will revert to 0 at the centre of the map (-1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	posmap = strtol(cmdbuffer, &endptr, 10);

	if ( (posmap == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (posmap < 0) {
	    printf("Escaping command...\n");
	    posmap = tmp_posmap;
	    break;
	  }
	  else {
	    printf("\nReference map for co-ordinates is now set to %d\n\n", posmap);
	    validentry = 1;
	  }
	}
      }

      break;


    case 84: // Set which map should be used to define the position when plotting maps in WCS

      validentry = 0;
      tmp_scaletype = scaletype;

      printf("Scale type currently set to %d (DEFAULT: %d)\n", scaletype, DEFAULTSCALETYPE);

      while (validentry == 0) {

	cmdbuffer = readline("Select a scale type to use when mapping. 0 = Pixels, 1 = Arcsec, 2 = DMS, 3 = WCS (-1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	scaletype= strtol(cmdbuffer, &endptr, 10);

	if ( (scaletype == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (scaletype < 0) {
	    printf("Escaping command...\n");
	    scaletype = tmp_scaletype;
	    break;
	  }
	  else if (scaletype > 3) {
	    printf("\nInvalid entry, please try again...\n\n");
	    continue;
	  }
	  else {
	    printf("\nScale type is now set to %d\n\n", scaletype);
	    validentry = 1;
	  }
	}
      }

      break;


    case 85: // Set the position to plot the beam when making maps

      validentry = 0;
      escape = 0;
      tmp_beamposx = beamposx;

      printf("Percentage along the X axis currently set to %d (DEFAULT: %d)\n", beamposx, DEFAULTBEAMPOSX);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the beam location along the X axis as a percentage of the total length (integer >= 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	beamposx = strtol(cmdbuffer, &endptr, 10);

	if ( (beamposx == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (beamposx < 0) {
	    printf("Escaping command...\n");
	    beamposx = tmp_beamposx;
	    escape = 1;
	    break;
	  }
	  else {
	    validentry = 1;
	  }
	}
      }

      // Break the second loop if they want to escape
      if (escape == 1) {
	break;
      }

      validentry = 0;
      tmp_beamposy = beamposy;

      printf("Percentage along the Y axis currently set to %d (DEFAULT: %d)\n", beamposy, DEFAULTBEAMPOSY);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the beam location along the Y axis as a percentage (integer >= 0, -1 to escape (both values reset)): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	beamposy = strtol(cmdbuffer, &endptr, 10);

	if ( (beamposy == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (beamposy < 0) {
	    printf("Escaping command...\n");
	    beamposy = tmp_beamposy;
	    beamposx = tmp_beamposx;
	    break;
	  }
	  else {
	    printf("\nThe beam will now be drawn %d per cent  the X axis and %d pixels from the Y axis\n\n", beamposx, beamposy);
	    validentry = 1;
	  }
	}
      }

      break;


    case 86: // Turn on and off the plotting of the beam to maps

      if (plotbeam == 1) {
	plotbeam = 0;
	printf("\nPlotting of the beam to maps has now been turned off.\n\n");
      }
      else if (plotbeam == 0) {
	plotbeam = 1;
	printf("\nPlotting of the beam to maps has now been turned on.\n\n");
      }
      else {
	printf("\nThe value of plotbeam is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", plotbeam);
      }

      break;


    case 87: // Export data to a text file
      
      if (numdatasets > 0) {

	// Open the directory that us currently set and check it exists
	dir = opendir(dataloc);

	if (dir != NULL) {
	  printf("\nOutput directory exists. All data will be output to %s.\n", dataloc);
	  closedir (dir);
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output data (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in dataloc ***\n\n", dataloc);
	  break;
	}

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);
	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which you wish to export data (-1 to escape, integers only): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if (currentset >= numdatasets) { // Check the set exists
	    printf("\nThis value is greater than the number of data sets that exists! Please try again...\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	// Get which data needs exporting
	validentry = 0;
	escape = 0;

	// List the available choices
	printf("======================================================\n");
	printf("  Selection  |  Data type                             \n");
	printf("------------------------------------------------------\n");
	printf("      0      |  JP ages                               \n");
	printf("      1      |  KP ages                               \n");
	printf("      2      |  Tribble (JP) ages                     \n");
	printf("      3      |  JP chi-squared                        \n");
	printf("      4      |  KP chi-squared                        \n");
	printf("      5      |  Tribble (JP) chi-squared              \n");
	printf("      6      |  JP errors                             \n");
	printf("      7      |  KP errors                             \n");
	printf("      8      |  Tribble (JP) errors                   \n");
	printf("      9      |  JP normalization                      \n");
	printf("      10     |  KP normalization                      \n");
	printf("      11     |  Tribble (JP) normalization            \n");
	printf("      12     |  Spectral index                        \n");
	printf("      13     |  Injection index                       \n");
	printf("      14     |  Injection index by region             \n");
	printf("      15     |  Injection index chisquared by region  \n");
	printf("      16     |  Region flux                           \n");
	printf("      17     |  Region array                          \n");
	printf("      -1     |  Escape command                        \n");
	printf("======================================================\n");
    

	while (validentry == 0) {

	  /*cmdbuffer = readline("Which data would you like to export? (-1 to escape, 0: JP Ages, 1: KP Ages, 2: JP X^2, 3: KP X^2, 4: JP Norm, 5:KP Norm, 6: Spectral Index, 7: Injection Index, 8: Tribble(JP) Ages, 9: Tribble(KP) Ages, 10: Tribble(JP) X^2, 11: Tribble(KP) X^2, 12: Tribble(JP) Norm, 13: Tribble(KP) Norm, 14: Region Flux 15: Region Array 16: JP Errors 17: KP Errors 18: Tribble(JP) Errors): ");*/

	  cmdbuffer = readline("Which data would you like to export? (see table above): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  dataexporttype = strtol(cmdbuffer, &endptr, 10);

	  if ( (dataexporttype == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }


	  if (dataexporttype < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if (dataexporttype > 17) {
	    printf("\nInvalid entry, please try again...\n\n");
	    continue;
	  }
	  // Check the set exists
	  else if ( ( (dataexporttype == 0) || (dataexporttype == 3) || (dataexporttype == 6) || (dataexporttype == 9) ) && (jpmodelset[currentset] != 1) ) {
	    printf("\nJP model has not yet been fit for this dataset, please try again...\n\n");
	    continue;
	  }
	  else if ( ( (dataexporttype == 1) || (dataexporttype == 4) || (dataexporttype == 7) || (dataexporttype == 10) ) && (kpmodelset[currentset] != 1) ) {
	    printf("\nKP model has not yet been fit for this dataset, please try again...\n\n");
	    continue;
	  }
	  else if ( ( (dataexporttype == 2) || (dataexporttype == 5) || (dataexporttype == 8) || (dataexporttype == 11) ) && (jptribmodelset[currentset] != 1) ) {
	    printf("\nTribble (JP) model has not yet been fit for this dataset, please try again...\n\n");
	    continue;
	  }
	  // Model no longer in use
	  /*else if ( ( (dataexporttype == 9) || (dataexporttype == 11) || (dataexporttype == 13) ) && (kptribmodelset[currentset] != 1) ) {
	    printf("\nTribble (KP) model has not yet been fit for this dataset, please try again...\n\n");
	    continue;
	    }*/
	  else if ( ( (dataexporttype == 16) || (dataexporttype == 17) ) && (regionsset[currentset] != 1) ) {
	    printf("\nsetregions has not yet been run for this dataset, please try again...\n\n");
	    continue;
	  }
	  else if ( (dataexporttype == 12)  && (specindexset[currentset] != 1) ) {
	    printf("\nSpectral indices has not yet been fit for this dataset, please try again...\n\n");
	    continue;
	  }
	  else if ( ( (dataexporttype == 13) || (dataexporttype == 14) || (dataexporttype == 15) ) && ( (mininjectset[currentset][0] < 1) && (mininjectset[currentset][1] < 1) && (mininjectset[currentset][2] < 1) && (mininjectset[currentset][3] < 1) ) ) {
	    printf("\nA valid injection index minimization has not yet been run for this dataset, please try again...\n\n");
	    continue;
	  }
	  else if ( (dataexporttype == 13) || (dataexporttype == 14) || (dataexporttype == 15) ) {

	    validentry2 = 0;

	    while (validentry2 == 0) {

	      cmdbuffer = readline("Which model would you like to export data for? (-1 to escape, 1: JP, 2: KP, 3: Tribble (JP), 4: Tribble(KP)): ");
	      if (cmdbuffer[0] != 0) {
		add_history(cmdbuffer);
	      }

	      chomp(cmdbuffer);

	      model = strtol(cmdbuffer, &endptr, 10);

	      if ( (model == 0) && (cmdbuffer == endptr) ) {
		printf("\nInvalid input, please try again...\n\n");
		continue;
	      }

	      if (model < 0) {
		printf("Escaping command...\n");
		escape = 1;
		break;
	      }
	      else if ( (model > NUMBEROFMODELS) || (model == 0) ) {
		printf("\nThis model number does not exist. Please try again...\n\n");
		continue;
	      }
	      else if (mininjectset[currentset][model] < 1) {
		printf("\nA valid injection index minimization has not yet been run for this model. Please try again...\n\n");
		continue;
	      }
	      else {
		validentry2 = 1;
		validentry = 1;
	      }
	    }

	    // Break the second loop if they want to escape
	    if (escape == 1) {
	      break;
	    }
	  }

	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	// Set which dataset we are exporting for naming purposes

	if (dataexporttype == 0) {
	  sprintf(datasettype,"JP");
	}
	else if (dataexporttype == 1) {
	  sprintf(datasettype,"KP");
	}
	else if (dataexporttype == 3) {
	  sprintf(datasettype,"JPX2");
	}
	else if (dataexporttype == 4) {
	  sprintf(datasettype,"KPX2");
	}
	else if (dataexporttype == 9) {
	  sprintf(datasettype,"JPNORM");
	}
	else if (dataexporttype == 10) {
	  sprintf(datasettype,"KPNORM");
	}
	else if (dataexporttype == 12) {
	  sprintf(datasettype,"SPECIND");
	}
	else if (dataexporttype == 13) {
	  sprintf(datasettype,"INJIND_MOD%d", model);
	}
	else if (dataexporttype == 2) {
	  sprintf(datasettype,"TRIBJP");
	}
	// Model no longer in use
	/*else if (dataexporttype == 9) {
	  sprintf(datasettype,"TRIBKP");
	  }*/
	else if (dataexporttype == 5) {
	  sprintf(datasettype,"TRIBJPX2");
	}
	/*else if (dataexporttype == 11) {
	  sprintf(datasettype,"TRIBKPX2");
	  }*/
	else if (dataexporttype == 11) {
	  sprintf(datasettype,"TRIBJPNORM");
	}
	/*else if (dataexporttype == 13) {
	  sprintf(datasettype,"TRIBKPNORM");
	  }*/
	else if (dataexporttype == 16) {
	  sprintf(datasettype,"FLUX");
	}
	else if (dataexporttype == 17) {
	  sprintf(datasettype,"REGARRAY");
	}
	else if (dataexporttype == 6) {
	  sprintf(datasettype,"JPERRORS");
	}
	else if (dataexporttype == 7) {
	  sprintf(datasettype,"KPERRORS");
	}
	else if (dataexporttype == 8) {
	  sprintf(datasettype,"TRIBJPERRORS");
	}
	else if (dataexporttype == 14) {
	  sprintf(datasettype,"INJECTBYREGION");
	}
	else if (dataexporttype == 15) {
	  sprintf(datasettype,"INJECTCHISQUAREDBYREGION");
	}
	else {
	  sprintf(datasettype,"Unknown");
	}


	// Get how the filename should be appended
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Please enter an identifier for the file output e.g. entering new will result in the file name 3C436_JP_new.txt for JP model data of 3C436 (esc to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  strcpy(tmp_filename, cmdbuffer);

	  if ( (tmp_filename == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }
	  else if ((strstr(tmp_filename, "esc") != NULL) && (strlen(tmp_filename) == 3)) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}


	// If the selction requires looping through the maps...
	if (dataexporttype == 16) {

	  printf("Exporting data, please wait... ");

	  for (i=0; i<imgnum[currentset]; i++) {

	    // Create the filename
	    sprintf(filename,"%s/%s_%.2e_%s_%s.%s", dataloc, settarget[currentset], frequencyobs[currentset][i], datasettype, tmp_filename, dataexportext);

	    // Check is the file exists or not

	    if ((filestream = fopen(filename, "r")) != NULL) {

	      validentry = 0;
	      escape = 0;

	      while (validentry == 0) {

		cmdbuffer = readline("Warning: File already exists, are you sure you want to overwrite (1: Yes 0: No): ");
		if (cmdbuffer[0] != 0) {
		  add_history(cmdbuffer);
		}

		chomp(cmdbuffer);

		confirm = strtol(cmdbuffer, &endptr, 10);

		if ( (confirm == 0) && (cmdbuffer == endptr) ) {
		  printf("\nInvalid input, please try again...\n\n");
		  continue;
		}
		else if ((confirm != 0) && (confirm != 1)) {
		  printf("\nInvalid entry, please try again...\n\n");
		  continue;
		}
		else if (confirm == 0) {
		  printf("Escaping command...\n");
		  escape = 1;
		  break;
		}
		else {
		  validentry = 1;
		}
	      }

	      fclose(filestream);

	      // Break the second loop if they want to escape
	      if (escape == 1) {
		break;
	      }
	    }
	  
	    if ( (filestream = fopen(filename,"w") ) == NULL ) {
	      printf("\n*** Error: Cannot open specified directory to output data (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in dataloc ***\n\n", dataloc);
	      break;
	    }
	    else {

	      for(j=1; j<=regnumber[currentset]; j++) {
		fprintf(filestream, "%.6e", (regflux[currentset][i][j] / regionsize[currentset][j]) * beamarea[currentset]);
		fprintf(filestream, "\n");
	      }
	      fclose(filestream);
	    }
	  }
	}
	else { // If we just need to go through it once...

	  // Create the filename
	  sprintf(filename,"%s/%s_%s_%s.%s", dataloc, settarget[currentset], datasettype, tmp_filename, dataexportext);

	  // Check is the file exists or not
	  if ((filestream = fopen(filename, "r")) != NULL) {

	    validentry = 0;
	    escape = 0;

	    while (validentry == 0) {

	      cmdbuffer = readline("Warning: File already exists, are you sure you want to overwrite (1: Yes 0: No): ");
	      if (cmdbuffer[0] != 0) {
		add_history(cmdbuffer);
	      }

	      chomp(cmdbuffer);

	      confirm = strtol(cmdbuffer, &endptr, 10);

	      if ( (confirm == 0) && (cmdbuffer == endptr) ) {
		printf("\nInvalid input, please try again...\n\n");
		continue;
	      }
	      else if ((confirm != 0) && (confirm != 1)) {
		printf("\nInvalid entry, please try again...\n\n");
		continue;
	      }
	      else if (confirm == 0) {
		printf("Escaping command...\n");
		escape = 1;
		break;
	      }
	      else {
		validentry = 1;
	      }
	    }

	    fclose(filestream);

	    // Break the second loop if they want to escape
	    if (escape == 1) {
	      break;
	    }

	  }

	  printf("Exporting data, please wait... ");

	  if ( (filestream = fopen(filename,"w") ) == NULL ) {
	    printf("\n*** Error: Cannot open specified directory to output data (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in dataloc ***\n\n", dataloc);
	    break;
	  }
	  else {

	    if (dataexporttype == 0) {
	      for(i=1; i<=regnumber[currentset]; i++) {

		fprintf(filestream, "%f", jpbestage[currentset][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    else if (dataexporttype == 1) {
	      for(i=1; i<=regnumber[currentset]; i++) {

		fprintf(filestream, "%f", kpbestage[currentset][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    else if (dataexporttype == 3) {
	      for(i=1; i<=regnumber[currentset]; i++) {

		fprintf(filestream, "%f", jpchisquared[currentset][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    else if (dataexporttype == 4) {
	      for(i=1; i<=regnumber[currentset]; i++) {

		fprintf(filestream, "%f", kpchisquared[currentset][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    else if (dataexporttype == 9) {
	      for(i=1; i<=regnumber[currentset]; i++) {

		fprintf(filestream, "%.6e", jpbestnorm[currentset][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    else if (dataexporttype == 10) {
	      for(i=1; i<=regnumber[currentset]; i++) {

		fprintf(filestream, "%.6e", kpbestnorm[currentset][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    else if (dataexporttype == 12) {
	      for(i=1; i<=regnumber[currentset]; i++) {

		fprintf(filestream, "%f %f", alpha[currentset][i], specindexerror[currentset][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    else if (dataexporttype == 13) {
	      	      
	      for(i=0; i<=mininjectset[currentset][model]; i++) {
		fprintf(filestream, "%f %f", mininjectstore[currentset][model] + ( ((maxinjectstore[currentset][model]-mininjectstore[currentset][model]) / mininjectset[currentset][model] ) * i ), inj_sumchisquared[currentset][model][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    else if (dataexporttype == 2) {
	      for(i=1; i<=regnumber[currentset]; i++) {

		fprintf(filestream, "%f", jptribbestage[currentset][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    // Model no longer in use
	    /*else if (dataexporttype == 9) {
	      for(i=1; i<=regnumber[currentset]; i++) {

		fprintf(filestream, "%f", kptribbestage[currentset][i]);
		fprintf(filestream, "\n");
	      }
	      }*/
	    else if (dataexporttype == 5) {
	      for(i=1; i<=regnumber[currentset]; i++) {

		fprintf(filestream, "%f", jptribchisquared[currentset][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    /*else if (dataexporttype == 11) {
	      for(i=1; i<=regnumber[currentset]; i++) {

		fprintf(filestream, "%f", kptribchisquared[currentset][i]);
		fprintf(filestream, "\n");
	      }
	      }*/
	    else if (dataexporttype == 11) {
	      for(i=1; i<=regnumber[currentset]; i++) {

		fprintf(filestream, "%.6e", jptribbestnorm[currentset][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    /*else if (dataexporttype == 13) {
	      for(i=1; i<=regnumber[currentset]; i++) {

		fprintf(filestream, "%.6e", kptribbestnorm[currentset][i]);
		fprintf(filestream, "\n");
	      }
	      }*/
	    

	    // The x/y --> i/j mapping has become crossed over somewhere. General mapping works so have just swapped them over for now but need to find the root cause of this problem.
	    else if (dataexporttype == 17) {
	      for(j=0; j<xdim[currentset]; j++) {
		for(i=0; i<ydim[currentset]; i++) {
		  if (regionarray[currentset][i][j] > 0) {
		    fprintf(filestream, "%d %d %d", j, i, regionarray[currentset][i][j]);
		    fprintf(filestream, "\n");
		  }
		}
	      }
	    }
	    else if (dataexporttype == 6) {
	      for(i=1; i<=regnumber[currentset]; i++) {
		fprintf(filestream, "+%f -%f", jpageerrorsplus[currentset][i], jpageerrorsminus[currentset][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    else if (dataexporttype == 7) {
	      for(i=1; i<=regnumber[currentset]; i++) {
		fprintf(filestream, "+%f -%f", kpageerrorsplus[currentset][i], kpageerrorsminus[currentset][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    else if (dataexporttype == 8) {
	      for(i=1; i<=regnumber[currentset]; i++) {
		fprintf(filestream, "+%f -%f", jptribbleageerrorsplus[currentset][i], jptribbleageerrorsminus[currentset][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    else if (dataexporttype == 14) {
	      for(i=1; i<=regnumber[currentset]; i++) {
		fprintf(filestream, "%d %f", i, inj_regbestinject[currentset][model][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    else if (dataexporttype == 15) {
	      for(i=1; i<=regnumber[currentset]; i++) {
		fprintf(filestream, "%d %f", i, inj_regchisquared[currentset][model][i]);
		fprintf(filestream, "\n");
	      }
	    }
	    else {
	      printf("\n*** Error: Unable to output data! Datatype of %d is unknown. Exiting... ***\n\n", dataexporttype);
	      break;
	    }
	  }
	  fclose(filestream);
	}

	printf("Done\n");
      }

      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;



    case 88: // Set where exported data will be saved

      validentry = 0;
      strcpy(tmp_dataloc, dataloc);

      printf("Currently exporting data to %s (DEFAULT: %s)\n", dataloc, DEFAULTDATALOC);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a folder location to store exported data. Absolute or relative paths allowed. No training slash (esc to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	strcpy(dataloc, cmdbuffer);

	if ((strstr(dataloc, "esc") != NULL) && (strlen(dataloc) == 3) ) {
	  strcpy(dataloc, tmp_dataloc);
	  printf("Escaping command...\n");
	  break;
	}
	else {
	  // Open the directory that us currently set and check it exists
	  dir = opendir(dataloc);

	  if (dir != NULL) {
	    printf("\nExported data will now be saved to %s\n\n", dataloc);
	    validentry = 1;
	    closedir (dir);
	  }
	  else {
	    strcpy(dataloc, tmp_dataloc);
	    printf("\n*** Error: Cannot open specified directory to output data (%s)! Please check the directory exists and the correct permission are set. ***\n\n", dataloc);
	    continue;
	  }
	}
      }
    
      break;


    case 89: // Turn on and off the plotting of the beam to maps

      if (usereduced == 1) {
	usereduced = 0;
	printf("\nPlots and maps will now output standard X^2 values.\n\n");
      }
      else if (usereduced == 0) {
	usereduced = 1;
	printf("\nPlots and maps will now output reduced X^2 values.\n\n");
      }
      else {
	printf("\nThe value of usereduced is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", usereduced);
      }

      break;


    case 90: // Print out the confidence levels

 
      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);
	}

	printf("\n========================================================================\n");
    

	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which you wish to find the chi-squared confidence intervals (-1 to escape, 888 for fixed DoF, integers only): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    printf("\nThis value is greater than the number of data sets that exists! Please try again...\n\n");
	    continue;
	  }
	  else {

	    if (currentset != 888) {
	      dof = (imgnum[currentset] - 2);
	    }

	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}
      }
      else {
	printf("\nNo data sets loaded. Entering fixed DoF mode...\n\n");
      }
    
      if ( (numdatasets == 0) || (currentset == 888) ) {

	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the number of degrees of freedom for which you wish to find the chi-squared confidence intervals (Integer > 0, -1 to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  dof = strtod(cmdbuffer, &endptr);

	  if ( (dof == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }
	  else {

	    if (dof < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else if (dof > 4.7e4) {
	    printf("\nUnable to calculate a CDF for that many DoF (max 4.7e4), please try again...\n\n");
	    continue;
	    }
	    else {
	      validentry = 1;
	    }
	  }
	}

	if (escape == 1) {
	  break;
	}
      }

      // print out the confidence levels
      interval = 0.68;
      siglevel68 = gsl_cdf_chisq_Pinv(interval, dof);
      printf("\nModel can be rejected with 68 per cent confidence at X^2 > %.2f (Reduced: %.2f)\n", siglevel68, siglevel68 / dof);
      interval = 0.90;
      siglevel90 = gsl_cdf_chisq_Pinv(interval, dof);
      printf("Model can be rejected with 90 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel90, siglevel90 / dof);
      interval = 0.95;
      siglevel95 = gsl_cdf_chisq_Pinv(interval, dof);
      printf("Model can be rejected with 95 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel95, siglevel95 / dof);
      interval = 0.99;
      siglevel99 = gsl_cdf_chisq_Pinv(interval, dof);
      printf("Model can be rejected with 99 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel99, siglevel99 / dof);
      interval = 0.995;
      siglevel995 = gsl_cdf_chisq_Pinv(interval, dof);
      printf("Model can be rejected with 99.5 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel995, siglevel995 / dof);
      interval = 0.999;
      siglevel999 = gsl_cdf_chisq_Pinv(interval, dof);
      printf("Model can be rejected with 99.9 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel999, siglevel999 / dof);
      interval = 0.9999;
      siglevel9999 = gsl_cdf_chisq_Pinv(interval, dof);
      printf("Model can be rejected with 99.99 per cent confidence X^2 > %.2f (Reduced: %.2f)\n\n", siglevel9999, siglevel9999 / dof);

      break;


    case 91: // Set the number of years to plot for models

      validentry = 0;
      tmp_modelmyears = modelmyears;

      printf("Myr to plot when outputting models currently set to %d (DEFAULT: %d)\n", modelmyears, DEFAULTMODELMYEARS);

      while (validentry == 0) {

	// Input
	cmdbuffer = readline("Enter the maximum age in Myr to plot when outputting models (e.g. plotjpmodel). (Integer > 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	modelmyears = strtol(cmdbuffer, &endptr, 10);

	if ( (modelmyears == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {

	  if (modelmyears < 0) {
	    printf("Escaping command...\n");
	    modelmyears = tmp_modelmyears;
	    break;
	  }
	  else {
	    printf("\nAges up to %d Myr will now be output at intervals of 1 Myr\n\n", modelmyears);
	    validentry = 1;
	  }
	}
      }

      break;


    case 92: // Create a single region for a dataset
      
      if (numdatasets <= 0) {
	fprintf(stderr,"\nError: No datasets have yet been loaded!\n\n");
      }

      else {

	if (regmemallocated == 0 ) {

	  regflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  regionarray = (int ***)calloc(MAXDATASETS, sizeof(int *));
	  regnumber = (int *)calloc(MAXDATASETS, sizeof(int));
	  fluxerror = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  regionlocx = (float **)calloc(MAXDATASETS, sizeof(float *));
	  regionlocy = (float **)calloc(MAXDATASETS, sizeof(float *));
	  setaveraged = (int *)calloc(MAXDATASETS, sizeof(int));
	  regionsize = (int **)calloc(MAXDATASETS, sizeof(int *));
	  regmaxflux = (float **)calloc(MAXDATASETS,sizeof(float *));
	  regminflux = (float **)calloc(MAXDATASETS,sizeof(float *));

	  regmemallocated = 1;
	}


	// This is easier than sending over string arrays, although make list() rather redundant
	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");
	// Check at least 1 dataset is loaded

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}
  
	printf("\n========================================================================\n");

	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the index number of the data set to create regions for (-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }
	
	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    currentset = tmp_currentset;
	    escape = 1;
	    break;
	  }
	  else if ( (currentset >= numdatasets) && (currentset != 888) ) {
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
      
	  else {
	    if (currentset == 888) {
	      for (i=0; i<numdatasets; i++){

		singleregion(imgnum, xdim, ydim, flux, signaltonoise, rms, searcharea, regionarray, regflux, regnumber, fluxerror, averegions, regionlocx, regionlocy, fluxorder, i, maptomap, hotpixels, sigmalevel, titles, labels, axismarks, regionsize, regmaxflux, regminflux, fluxcalerror[i], onsourcenoise, regionsset[i], border[currentset], zoom, ra[i][posmap], dec[i][posmap], cellsize[i], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, xshift, yshift, largetxt);

		regionsset[i] = 1; // Mark the set as having regions defined
		setaveraged[i] = averegions;

	      }
	      currentset = 0; // Put the current set back to 0 to save any mix ups later
	      validentry = 1;
	    }
	
	    else {
	      singleregion(imgnum, xdim, ydim, flux, signaltonoise, rms, searcharea, regionarray, regflux, regnumber, fluxerror, averegions, regionlocx, regionlocy, fluxorder, currentset, maptomap, hotpixels, sigmalevel, titles, labels, axismarks, regionsize, regmaxflux, regminflux, fluxcalerror[currentset], onsourcenoise, regionsset[currentset], border[currentset], zoom, ra[currentset][posmap], dec[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, xshift, yshift, largetxt);

	      regionsset[currentset] = 1;
	      setaveraged[currentset] = averegions;
	      validentry = 1;
	    }
	  }
	}

	if (escape != 1) {
	  printf("*** Warning: It is left to the user to ensure the region used is large enough, sits in a sensible location and reaches the desired signal to noise ratio. The only cut made is that the flux per pixel is greater than sigma * RMS. This should not be an issue when region size >= beamsize as it is comparable to the old manual method of determining spectral ages and indices. *** \n\n");
	}
      }
    
      break;


    case 93: // Minimise the chisquared by varying the injection index.

      // Check if the min and max inject are set correctly
      if (mininject > maxinject) {
	fprintf(stderr,"\n *** Error: Invalid search range. Minimum injection index must smaller than the maximum. Please set a valid range via the mininject and maxinject commands and try again ***\n\n");
	break;
      }
  
      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set to evaluate (-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if (currentset == 888) {
	    
	    validentry2 = 0;

	    while (validentry2 == 0) {

	      cmdbuffer = readline("Computing all dataset can take a VERY long time. Are you sure you wish to continue? Yes (1) or no (0): ");
	      if (cmdbuffer[0] != 0) {
		add_history(cmdbuffer);
	      }

	      chomp(cmdbuffer);

	      warningcheck= strtol(cmdbuffer, &endptr, 10);

	      if ( (warningcheck == 0) && (cmdbuffer == endptr) ) {
		printf("\nInvalid input, please try again...\n\n");
		continue;
	      }

	      if ( (warningcheck < 0) || (warningcheck > 1) ) { // Check the value is either 1 or 0
		printf("\nValue must be either Yes (1) or no (0). Please try again...\n\n");
		continue;
	      }

	      else if (warningcheck == 0) { // If 0, exit
		break;
	      }

	      else if (warningcheck == 1) { // If 1, continue
		printf("Proceeding...\n");
		validentry2 = 1;
		validentry = 1;
	      }

	      else { // Catch anthing odd
		printf("\nThe value of warningcheck  is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", warningcheck);
		break;
	      }
	    }
	  }

	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    printf("\nThis value is greater than the number of data sets that exists! Please try again...\n\n");
	    continue;
	  }
	  else if ( (imgnum[currentset] <= 1) && (currentset != 888) ) { // Check the set has at least two maps
	    fprintf(stderr,"\nError: The data set must contain a minimum of two images for fitting. Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    continue;
	  }

	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Which model would you like to minimise the injection index for? (-1 to escape, 1: JP, 2: KP, 3: Tribble (JP), 4: Tribble (KP), 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  tmp_model = strtol(cmdbuffer, &endptr, 10);

	  if ( (tmp_model == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (tmp_model < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if ( (tmp_model > NUMBEROFMODELS) && (tmp_model != 888) ) { 
	    printf("\nInvalid model number, please try again...\n\n");
	    escape = 1;
	    continue;
	  }

	  else if (tmp_model == 888) {

	    validentry2 = 0;

	    while (validentry2 == 0) {

	      cmdbuffer = readline("Computing all models can take a VERY long time. Are you sure you wish to continue? Yes (1) or no (0): ");
	      if (cmdbuffer[0] != 0) {
		add_history(cmdbuffer);
	      }

	      chomp(cmdbuffer);

	      warningcheck = strtol(cmdbuffer, &endptr, 10);

	      if ( (warningcheck == 0) && (cmdbuffer == endptr) ) {
		printf("\nInvalid input, please try again...\n\n");
		continue;
	      }

	      if ( (warningcheck < 0) || (warningcheck > 1) ) { // Check the value is either 1 or 0
		printf("\nValue must be either Yes (1) or no (0). Please try again...\n\n");
		continue;
	      }
	      else if (warningcheck == 0) { // If 0, exit
		break;
	      }
	      else if (warningcheck == 1) { // If 1, continue
		printf("Proceeding...\n");
		validentry2 = 1;
		validentry = 1;
	      }
	      else { // Catch anthing odd
		printf("\nThe value of warningcheck  is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", warningcheck);
		break;
	      }
	    }
	  }

	  else if ( ( (tmp_model > NUMBEROFMODELS) && (tmp_model != 888) ) || (tmp_model == 0) ) {
	    printf("\nInvalid entry, please try again...\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	// Setup the loops
	if (currentset == 888) {
	  looplimit = numdatasets;
	}
	else {
	  looplimit = 1;
	}

	escape = 0;

	for (j=0; j<looplimit; j++){

	  if (looplimit > 1) {
	    currentset = j;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = j;
	  }

	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets

	  if (regionsset[currentset] != 1) {
	    printf("\nThe regions for data set %d have not been set! Skipping this set...\n", j);
	    continue;	  
	  }
	  else if (imgnum[currentset] <= 1) { // Check the set exists
	    fprintf(stderr,"\n*** Warning: Data set %d does not contain the minimum two images required for fitting. Skipping this set... ***\n\n", j);
	    continue;
	  }

	  
	  if (tmp_model == 888) {
	    modelloop = NUMBEROFMODELS;
	  }
	  else {
	    modelloop = 1;
	  }

	  for (a=1; a<=modelloop; a++) { // Models are indexed from 1

	    if (tmp_model == 888) {
	      model = a;
	    }
	    else {
	      model = tmp_model;
	    }
	  
	    // Setup the memory
	    inj_sumchisquared[currentset][model] = (float *)realloc(inj_sumchisquared[currentset][model], (injectinterval+1) * sizeof(float));
	    inj_regchisquared[currentset][model] = (float *)realloc(inj_regchisquared[currentset][model], (regnumber[currentset]+1) * sizeof(float));
	    inj_regbestinject[currentset][model] = (float *)realloc(inj_regbestinject[currentset][model], (regnumber[currentset]+1) * sizeof(float));

	    // Set the default values
	    for (i=0; i<=injectinterval; i++) {
	      inj_sumchisquared[currentset][model][i] = 0.0;
	    }
	    for (i=0; i<=regnumber[currentset]; i++) {
	      inj_regchisquared[currentset][model][i] = 1e23;
	      inj_regbestinject[currentset][model][i] = 1e23;
	    }


	    for (i=0; i<=injectinterval; i++) {

	      tmp_inject = mininject + ( ((maxinject-mininject) / (float)injectinterval) * i );

	      printf("Calculating the best chisquared for dataset %d using model %d with injection index %.2f\n", currentset, model, tmp_inject);


	      if (spectralageing_min(regflux[currentset], imgnum[currentset], fluxorder[currentset], frequencyobs[currentset], regnumber[currentset], fluxerror[currentset], usr_gmin, usr_gmax, minmyears, myears, ageresolution, levels, tmp_inject, fieldstrength, model, &inj_sumchisquared[currentset][model][i], inj_regchisquared[currentset][model], inj_regbestinject[currentset][model], redshift[currentset]) == 5) {

		printf("\n *** Error: Unable to calculate minimum injection index. Please review any errors and try again ***\n\n");
		escape = 1;
		break;
	      }

	      mininjectset[currentset][model] = injectinterval;
	      maxinjectstore[currentset][model] = maxinject;
	      mininjectstore[currentset][model] = mininject;
	    }
	    if (escape == 1) {
	      break;
	    }
	  }
	  if (escape == 1) {
	    break;
	  }
	}
	if (escape == 1) {
	  break;
	}

	// Output a summary
	printf("\nInjection index fit summary:\n\n");

	for (j=0; j<looplimit; j++) {

	  if (looplimit > 1) { // If it is only one dataset, we can keep the currentset number from before
	    currentset = j;
	  }
	 

	  for (a=1; a<=modelloop; a++) { // Models are indexed from 1

	    if (tmp_model == 888) {
	      model = a;
	    }
	    else {
	      model = tmp_model;
	    }

	    bestchi1 = 1e+24;
	    bestchi2 = 1e+24;

	    for (i=0; i<=injectinterval; i++) {

	      if (inj_sumchisquared[currentset][model][i] < bestchi1) {
		bestchi2 = bestchi1;
		bestinject2 = bestinject1;
		bestchi1 = inj_sumchisquared[currentset][model][i];
		bestinject1 = mininject + ( ((maxinject-mininject) / injectinterval) * i );
	      }
	      else if (inj_sumchisquared[currentset][model][i] < bestchi2) {
		bestchi2 = inj_sumchisquared[currentset][model][i];
		bestinject2 = mininject + ( ((maxinject-mininject) / injectinterval) * i );
	      }
	    }

	    printf("The injection index for dataset %d when fitting model %d has a best fit between %.2f and %.2f\n", currentset, model, bestinject2, bestinject1);

	  }
	}
      
	/*for (i=1; i<=regnumber[currentset]; i++) {
	  printf("Region: %d Inject: %.4f Chi: %.4f\n", i, inj_regbestinject[currentset][model][i], inj_regchisquared[currentset][model][i]);
	  }*/

	printf("\nInjection index fitting complete. These results can be output using the exportdata command for the purposes of analysis.\n\n");

      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;


    case 94: // Set the minimum injection index to be tested

      validentry = 0;
      tmp_mininject = mininject;

      printf("Minimum injection index currently set to %.2f (DEFAULT: %.2f)\n", mininject, DEFAULTMININJECT);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the minimum injection index for findinject fit (Double > 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	mininject = strtod(cmdbuffer, &endptr);

	if ( (mininject == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (mininject < 0) {
	    printf("Escaping command...\n");
	    mininject = tmp_mininject;
	    break;
	  }/* // This is now done checked at run time
	      else if (mininject >= maxinject) {
	      printf("Mininject must be smaller than maxinject. Please either change maxinject first or try again...\n");
	      mininject = tmp_mininject;
	      }*/
	  else {
	    printf("\nThe minimum injection index is now set to %.2f\n\n", mininject);
	    validentry = 1;
	  }
	}
      }

      break;


    case 95: // Set the minimum injection index to be tested

      validentry = 0;
      tmp_maxinject = maxinject;

      printf("Maximum injection index currently set to %.2f (DEFAULT: %.2f)\n", maxinject, DEFAULTMAXINJECT);

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value for the maximum injection index for findinject fit (Double > 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	maxinject = strtod(cmdbuffer, &endptr);

	if ( (maxinject == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (maxinject < 0) {
	    printf("Escaping command...\n");
	    maxinject = tmp_maxinject;
	    break;
	  }/* // This is now done checked at run time
	      else if (mininject >= maxinject) {
	      printf("Maxinject must be greater than mininject. Please either change mininject first or try again...\n");
	      maxinject = tmp_maxinject;
	      }*/
	  else {
	    printf("\nThe maximum injection index is now set to %.2f\n\n", maxinject);
	    validentry = 1;
	  }
	}
      }

      break;


    case 96: // Set the number of injection index intervals to be tested

      validentry = 0;
      tmp_injectinterval = injectinterval;

      printf("Number of intervals currently set to %d (DEFAULT: %d)\n", injectinterval, DEFAULTINJECTINTERVAL);

      while (validentry == 0) {

	cmdbuffer = readline("Enter the number of intervals to fit between mininject and maxinject. n + 1 intevals are tested to find the best fit e.g. with mininject 0.5 and maxinject 1.0, n = 10 will result in steps of 0.05 (integer > 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	injectinterval = strtol(cmdbuffer, &endptr, 10);

	if ( (injectinterval == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (injectinterval < 0) {
	    printf("Escaping command...\n");
	    injectinterval = tmp_injectinterval;
	    break;
	  }
	  else {
	    printf("\nThe injection index will now be fitted over %d intervals\n\n", injectinterval);
	    validentry = 1;
	  }
	}
      }

      break;


    case 97: // View fluxcalerrors

      printf("=====================================================================\n\n");
      printf("Dataset |  Index  |   Frequency (Obs.)   |   Flux Calibration Error \n\n");

      // Loop through each of the sets
      for (i=0; i<numdatasets; i++) {

	for (a=0; a<imgnum[i]; a++) {

	  printf("   %d         %d       %.2e Hz      %.1f%%\n", i, fluxorder[i][a], frequencyobs[i][fluxorder[i][a]], fluxcalerror[i][fluxorder[i][a]]*100);

	}
      }

      printf("\n=====================================================================\n\n");

      break;



    case 98: // Select the telescope to set default values for

      validentry = 0;
      tmp_instrument = instrument;

      printf("Default telescope currently set to %d (DEFAULT: %d)\n", instrument, DEFAULTINSTRUMENT);

      while (validentry == 0) {

	cmdbuffer = readline("Select a telescope to use for the default values on load (0 = Flat Errors, 1 = JVLA, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	instrument = strtol(cmdbuffer, &endptr, 10);

	if ( (instrument == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {

	  if (instrument < 0) {
	    printf("Escaping command...\n");
	    instrument = tmp_instrument;
	    break;
	  }
	  else if (instrument > 1) {
	    printf("\nInvalid selection. Please try again...\n\n");
	    continue;
	  }
	  else {
	    if (instrument > 0) {
	      printf("\nDefault values on load will now be set to telescope %d\n\n", instrument);
	    }
	    validentry = 1;
	  }
	}
      }


      if (instrument == 0) {

	validentry = 0;
	tmp_flaterror = flaterror;

	printf("Percentage error currently set to %.2f (DEFAULT: %.2f)\n", flaterror, DEFAULTERRFLAT);

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the percentage error to apply on load in decimal form (-1 to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	 flaterror = strtod(cmdbuffer, &endptr);

	  if ( (flaterror == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }
	  else {
	    if (flaterror < 0) {
	      printf("Escaping command...\n");
	      flaterror = tmp_flaterror;
	      break;
	    }
	    else {
	      if (flaterror > 0) {
		printf("\nDefault error will now be set to %.2f for all maps on load\n\n", flaterror);
	      }
	      validentry = 1;
	    }
	  }
	}
      }


      break;


    case 99: //Fit the JPTRIB model to a dataset
  
      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to fit the Tribble (JP) model (-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    printf("\nThis value is greater than the number of data sets that exists! Please try again...\n\n");
	    continue;
	  }
	  else if ( (imgnum[currentset] <= 1) && (currentset != 888) ) { // Check the set has at least two maps
	    fprintf(stderr,"\nError: The data set must contain a minimum of two images for fitting. Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    continue;
	  }

	  else {
	    validentry = 1;
	  }
	}
      
	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	// Setup the loops
	if (currentset == 888) {
	  looplimit = numdatasets;
	}
	else {
	  looplimit = 1;
	}

	for (j=0; j<looplimit; j++){

	  if (looplimit > 1) {
	    currentset = j;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = j;
	  }

	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets

	  if (regionsset[currentset] != 1) {
	    printf("\nThe regions for data set %d have not been set! Skipping this set...\n", j);
	    continue;	  
	  }
	  else if (imgnum[currentset] <= 1) { // Check the set exists
	    fprintf(stderr,"\n*** Warning: Data set %d does not contain the minimum two images required for fitting. Skipping this set... ***\n\n", j);
	    continue;
	  }
	  
	  setmodelres[currentset][3] = modelres; // Set the resolution for this dataset and model

	  // Setup the memory
	  if (jptribmodelmemoryset == 0) {

	    jptribchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jptribbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jptribmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jptribmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	    jptribbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jptribbleageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jptribbleageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));

	    jptribmodelmemoryset = 1;
	  }
	

	  if (numdatasets > 1) {

	    allocarray = (int *)calloc(MAXDATASETS, sizeof(int));

	    realloc_jptribchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_jptribbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_jptribmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_jptribmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	    realloc_jptribbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_jptribbleageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_jptribbleageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));


	    for (o=0; o<numdatasets; o++) {

	      if (jptribfitmemset[o] == 1) {

		realloc_jptribchisquared[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_jptribbestage[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_jptribmodelalpha[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_jptribmodelflux[o] = (float **)calloc(regnumber[o]+1, sizeof(float *));
		realloc_jptribbestnorm[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_jptribbleageerrorsplus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_jptribbleageerrorsminus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));

		allocarray[o] = 1;

	      }
	    }
	  
	    for (o=0; o<numdatasets; o++) {

	      if (jptribfitmemset[o] == 1) {

		for (q=0; q<=regnumber[o]; q++) {

		  //printf("q=%d\n", q);
		  realloc_jptribmodelflux[o][q] = (float *)calloc(setmodelres[currentset][3]+1, sizeof(float));
		}
	      }
	    }

	    for (o=0; o<numdatasets; o++) {

	      for (p=0; p<=regnumber[o]; p++) {

		if (jptribfitmemset[o] == 1) {

		  realloc_jptribchisquared[o][p] = jptribchisquared[o][p];
		  realloc_jptribbestage[o][p] = jptribbestage[o][p];
		  realloc_jptribmodelalpha[o][p] = jptribmodelalpha[o][p];
		  realloc_jptribbestnorm[o][p] = jptribbestnorm[o][p];
		  realloc_jptribbleageerrorsplus[o][p] = jptribbleageerrorsplus[o][p];
		  realloc_jptribbleageerrorsminus[o][p] = jptribbleageerrorsminus[o][p];

		  for (q=0; q<=setmodelres[currentset][3]; q++) {
		    realloc_jptribmodelflux[o][p][q] = jptribmodelflux[o][p][q];
		  }

		}
	      }

	    }
	  
	    for (o=0; o<numdatasets-1; o++) { // -1 as we are yet to allocate the new dataset memory

	      free(jptribchisquared[o]);
	      free(jptribbestage[o]);
	      free(jptribmodelalpha[o]);
	      free(jptribbestnorm[o]);
	      free(jptribmodelflux[o]);
	      free(jptribbleageerrorsplus[o]);
	      free(jptribbleageerrorsminus[o]);

	    }
	    
	    free(jptribchisquared);
	    free(jptribbestage);
	    free(jptribmodelalpha);
	    free(jptribbestnorm);
	    free(jptribbleageerrorsplus);
	    free(jptribbleageerrorsminus);
	    free(jptribmodelflux);
	    
	    jptribfitmemset[o] = 0;
	    jptribmodelmemoryset = 0;


	    // Make the first level of the array again
	    jptribchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jptribbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jptribmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jptribmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	    jptribbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jptribbleageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	    jptribbleageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));

	    jptribmodelmemoryset = 1;

	  }
	
	  for (o=0; o<numdatasets; o++) {

	    jptribchisquared[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    jptribbestage[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    jptribmodelalpha[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    jptribmodelflux[o] = (float **)calloc(regnumber[o]+1, sizeof(float *));
	    jptribbestnorm[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    jptribbleageerrorsplus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    jptribbleageerrorsminus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));


	    for (p=0; p<=regnumber[o]; p++) {
	      jptribmodelflux[o][p] = (float *)calloc(setmodelres[currentset][3]+1, sizeof(float));
	    }
	    
	    jptribfitmemset[o] = 1;

	  }
	
	
	  if (numdatasets > 1) {
	    for (o=0; o<numdatasets; o++) {

	      if ( (o != currentset) && (allocarray[o] == 1) ) {

		for (p=0; p<=regnumber[o]; p++) {

		  jptribchisquared[o][p] = realloc_jptribchisquared[o][p];
		  jptribbestage[o][p] = realloc_jptribbestage[o][p];
		  jptribmodelalpha[o][p] = realloc_jptribmodelalpha[o][p];
		  jptribbestnorm[o][p] = realloc_jptribbestnorm[o][p];
		  jptribbleageerrorsplus[o][p] = realloc_jptribbleageerrorsplus[o][p];
		  jptribbleageerrorsminus[o][p] = realloc_jptribbleageerrorsminus[o][p];

		  for (q=0; q<=setmodelres[currentset][3]; q++) {

		    jptribmodelflux[o][p][q] = realloc_jptribmodelflux[o][p][q];
		  }
		}
	      }
	    }
	  }
	
	  if (numdatasets > 1) {

	    // Free the reallocation arrays if they are set

	    for (o=0; o<numdatasets-1; o++) {

	      free(realloc_jptribchisquared[o]);
	      free(realloc_jptribbestage[o]);
	      free(realloc_jptribmodelalpha[o]);
	      free(realloc_jptribbestnorm[o]);
	      free(realloc_jptribbleageerrorsplus[o]);
	      free(realloc_jptribbleageerrorsminus[o]);
	      //free(realloc_jptribmodelflux[o]);

	    }

	    free(realloc_jptribchisquared);
	    free(realloc_jptribbestage);
	    free(realloc_jptribmodelalpha);
	    free(realloc_jptribbestnorm);
	    free(realloc_jptribbleageerrorsplus);
	    free(realloc_jptribbleageerrorsminus);
	    free(realloc_jptribmodelflux);

	    free(allocarray);

	  }

	  model = 3;

	  // Using a temp array here to allow all models without the need for dummy arrays
	  tmp_ageerrorsplus = (float *)calloc(regnumber[currentset]+1, sizeof(float));
	  tmp_ageerrorsminus = (float *)calloc(regnumber[currentset]+1, sizeof(float));

	  spectralageingmodels(regflux[currentset], imgnum[currentset], fluxorder[currentset], frequency[currentset], regnumber[currentset], jptribchisquared[currentset], jptribbestage[currentset], fluxerror[currentset], jptribmodelalpha[currentset], jptribmodelflux[currentset], usr_gmin, usr_gmax, minmyears, myears, ageresolution, levels, jptribbestnorm[currentset], setmodelres[currentset][3], printresults, inject, fieldstrength, model, printreject, redshift[currentset], tmp_ageerrorsplus, tmp_ageerrorsminus, suppresscdf);

	  jptribfield[currentset] = fieldstrength;
	  jptribinject[currentset] = inject;
	  jptribmodelset[currentset] = 1;

	  // Transfer the age errors over to the perminant array
	  for (r=1; r<=regnumber[currentset]; r++) {
	    jptribbleageerrorsplus[currentset][r] = tmp_ageerrorsplus[r];
	    jptribbleageerrorsminus[currentset][r] = tmp_ageerrorsminus[r];
	  }

	  free(tmp_ageerrorsplus);
 	  free(tmp_ageerrorsminus); 
	}

      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;

      
    case 100: //Fit the KPTRIB model to a dataset ***WARNING: THIS IS NO LONGER IN USE!!!! CODE MAY NOT RUN AS EXPECTED AND IS ONLY HERE FOR LEGACY PURPOSES. TO BE REMOVED***
  
      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  printf("Enter the data set for which to fit the Tribble (KP) model (-1 to escape, 888 for all): ");

	  if (scanf("%d", &currentset) == 0) {
		 
	    printf("Invalid input, escaping command to protect the existing data...\n");
	    escape = 1;
	    break;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    printf("This value is greater than the number of data sets that exists! Please try again...\n");
	  }

	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\nExiting command...\n\n");
	    escape = 1;
	    break;
	  }

	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	// Setup the loops
	if (currentset == 888) {
	  looplimit = numdatasets;
	}
	else {
	  looplimit = 1;
	}

	for (j=0; j<looplimit; j++){

	  if (looplimit > 1) {
	    currentset = j;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = j;
	  }

	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets

	  if (regionsset[currentset] != 1) {
	    printf("\nThe regions for data set %d have not been set! Skipping this set...\n", j);
	    continue;	  
	  }

	  setmodelres[currentset][4] = modelres; // Set the resolution for this dataset and model

	  // Setup the memory
	  if (kptribmodelmemoryset == 0) {

	    kptribchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kptribbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kptribmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kptribmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	    kptribbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));

	    kptribmodelmemoryset = 1;
	  }
	

	  if (numdatasets > 1) {

	    allocarray = (int *)calloc(MAXDATASETS, sizeof(int));

	    realloc_kptribchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_kptribbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_kptribmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	    realloc_kptribmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	    realloc_kptribbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));


	    for (o=0; o<numdatasets; o++) {

	      if (kptribfitmemset[o] == 1) {

		realloc_kptribchisquared[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_kptribbestage[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_kptribmodelalpha[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
		realloc_kptribmodelflux[o] = (float **)calloc(regnumber[o]+1, sizeof(float *));
		realloc_kptribbestnorm[o] = (float *)calloc(regnumber[o]+1, sizeof(float));

		allocarray[o] = 1;

	      }
	    }
	  
	    for (o=0; o<numdatasets; o++) {

	      if (kptribfitmemset[o] == 1) {

		for (q=0; q<=regnumber[o]; q++) {

		  //printf("q=%d\n", q);
		  realloc_kptribmodelflux[o][q] = (float *)calloc(setmodelres[currentset][4]+1, sizeof(float));
		}
	      }
	    }

	    for (o=0; o<numdatasets; o++) {

	      for (p=0; p<=regnumber[o]; p++) {

		if (kptribfitmemset[o] == 1) {

		  realloc_kptribchisquared[o][p] = kptribchisquared[o][p];
		  realloc_kptribbestage[o][p] = kptribbestage[o][p];
		  realloc_kptribmodelalpha[o][p] = kptribmodelalpha[o][p];
		  realloc_kptribbestnorm[o][p] = kptribbestnorm[o][p];

		  for (q=0; q<=setmodelres[currentset][4]; q++) {
		    realloc_kptribmodelflux[o][p][q] = kptribmodelflux[o][p][q];
		  }

		}
	      }

	    }
	  
	    for (o=0; o<numdatasets-1; o++) { // -1 as we are yet to allocate the new dataset memory

	      free(kptribchisquared[o]);
	      free(kptribbestage[o]);
	      free(kptribmodelalpha[o]);
	      free(kptribbestnorm[o]);
	      free(kptribmodelflux[o]);

	    }
	    
	    free(kptribchisquared);
	    free(kptribbestage);
	    free(kptribmodelalpha);
	    free(kptribbestnorm);
	    free(kptribmodelflux);
	    
	    kptribfitmemset[o] = 0;
	    kptribmodelmemoryset = 0;


	    // Make the first level of the array again
	    kptribchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kptribbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kptribmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	    kptribmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	    kptribbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));

	    kptribmodelmemoryset = 1;

	  }
	
	  for (o=0; o<numdatasets; o++) {

	    kptribchisquared[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    kptribbestage[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    kptribmodelalpha[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	    kptribmodelflux[o] = (float **)calloc(regnumber[o]+1, sizeof(float *));
	    kptribbestnorm[o] = (float *)calloc(regnumber[o]+1, sizeof(float));


	    for (p=0; p<=regnumber[o]; p++) {
	      kptribmodelflux[o][p] = (float *)calloc(setmodelres[currentset][4]+1, sizeof(float));
	    }
	    
	    kptribfitmemset[o] = 1;

	  }
	
	
	  if (numdatasets > 1) {
	    for (o=0; o<numdatasets; o++) {

	      if ( (o != currentset) && (allocarray[o] == 1) ) {

		for (p=0; p<=regnumber[o]; p++) {

		  kptribchisquared[o][p] = realloc_kptribchisquared[o][p];
		  kptribbestage[o][p] = realloc_kptribbestage[o][p];
		  kptribmodelalpha[o][p] = realloc_kptribmodelalpha[o][p];
		  kptribbestnorm[o][p] = realloc_kptribbestnorm[o][p];

		  for (q=0; q<=setmodelres[currentset][4]; q++) {

		    kptribmodelflux[o][p][q] = realloc_kptribmodelflux[o][p][q];
		  }
		}
	      }
	    }
	  }
	
	  if (numdatasets > 1) {

	    // Free the reallocation arrays if they are set

	    for (o=0; o<numdatasets-1; o++) {

	      free(realloc_kptribchisquared[o]);
	      free(realloc_kptribbestage[o]);
	      free(realloc_kptribmodelalpha[o]);
	      free(realloc_kptribbestnorm[o]);
	      //free(realloc_kptribmodelflux[o]);

	    }

	    free(realloc_kptribchisquared);
	    free(realloc_kptribbestage);
	    free(realloc_kptribmodelalpha);
	    free(realloc_kptribbestnorm);
	    free(realloc_kptribmodelflux);

	    free(allocarray);

	  }

	  model = 4;

	  // Using a temp array here to allow all models without the need for dummy arrays
	  tmp_ageerrorsplus = (float *)calloc(regnumber[currentset]+1, sizeof(float));
	  tmp_ageerrorsminus = (float *)calloc(regnumber[currentset]+1, sizeof(float));

	  spectralageingmodels(regflux[currentset], imgnum[currentset], fluxorder[currentset], frequency[currentset], regnumber[currentset], kptribchisquared[currentset], kptribbestage[currentset], fluxerror[currentset], kptribmodelalpha[currentset], kptribmodelflux[currentset], usr_gmin, usr_gmax, minmyears, myears, ageresolution, levels, kptribbestnorm[currentset], setmodelres[currentset][4], printresults, inject, fieldstrength, model, printreject, redshift[currentset], tmp_ageerrorsplus, tmp_ageerrorsminus, suppresscdf);

	  kptribfield[currentset] = fieldstrength;
	  kptribinject[currentset] = inject;
	  kptribmodelset[currentset] = 1;


	  free(tmp_ageerrorsplus);
	  free(tmp_ageerrorsminus);

	}

      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;
      

    case 101: // Plot the Tribble (JP) model for an arbitary normalisation

      printf("\nPlotting an example Tribble (JP) model between %.2e Hz and %.2e Hz\n\n", minmodelfreq, maxmodelfreq);

      if (export == 1) {

	// Open the directory that is currently set and check it exists
	dir = opendir(imageloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	  break;
	}

	time( &currenttime );
	time_struct = localtime ( &currenttime );

	strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	sprintf(output,"%s/ExTribJPModel_%s.%s/%s", imageloc,  timebuff, imagetype, imagetype);

	printf("\nExporting as %s\n", output);

	cpgopen(output);
      }
      else {
	sprintf(output,"/xs");
      }

      model = 3;

      plotmodel(minmodelfreq, maxmodelfreq, inject, fieldstrength, model, usr_gmin, usr_gmax, minmodelmyears, modelmyears, output, titles, exampleredshift, skip);

      break;


    case 102: // Plot the Tribble (KP) model for an arbitary normalisation

      printf("\nPlotting an example Tribble (KP) model between %.2e Hz and %.2e Hz\n\n", minmodelfreq, maxmodelfreq);

      if (export == 1) {

	// Open the directory that is currently set and check it exists
	dir = opendir(imageloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	  break;
	}

	time( &currenttime );
	time_struct = localtime ( &currenttime );

	strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	sprintf(output,"%s/ExTribKPModel_%s.%s/%s", imageloc,  timebuff, imagetype, imagetype);

	printf("\nExporting as %s\n", output);

	cpgopen(output);
      }
      else {
	sprintf(output,"/xs");
      }

      model = 4;

      plotmodel(minmodelfreq, maxmodelfreq, inject, fieldstrength, model, usr_gmin, usr_gmax, minmodelmyears, modelmyears, output, titles, exampleredshift, skip);

      break;


    case 103: // Merge maps in the image plane

      mergeimages(casadata, imageloc);

      break;


    case 104: // Reconstruct a map from the best values of a model at any frequency

      printf("\n ************************************************************************* \n");
      printf(" *                                                                       *\n");
      printf(" *   Warning: If parameters have changed since the original model fit    * \n");
      printf(" *   (e.g. bfield, injectionindex...), the reconstructed radio map may   *\n");
      printf(" *   not produced the expected result. This allows for more flexibility  *\n");
      printf(" *   but with an increased chance of user error. You have been warned!   *\n");
      printf(" *                                                                       *\n");
      printf(" ************************************************************************* \n\n");

      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to produce the reconstructed map (-1 to escape, integers only): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if (currentset >= numdatasets) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}
      
	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	validentry = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Which model fit should be used? (-1 to escape, 1 JP model, 2 KP model, 3 Tribble (JP), 4 Tribble (KP)): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  model = strtol(cmdbuffer, &endptr, 10);

	  if ( (model == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (model < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if ( (model > NUMBEROFMODELS) ||  (model == 0) ) { // Check the set exists
	    printf("\nInvalid model number. Please try again...\n\n");
	    continue;
	  }

	  else if ( (model == 1) && (jpmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: JP model vales have not yet been calculated for this data set! Please run fitjpmodel first.\n\n");
	    continue;
	  }

	  else if ( (model == 2) && (kpmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: KP model vales have not yet been calculated for this data set! Please run fitjpmodel first.\n\n");
	    continue;
	  }

	  else if ( (model == 3) && (jptribmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Tribble (KP) model vales have not yet been calculated for this data set! Please run fitjpmodel first.\n\n");
	    continue;
	  }

	  else if ( (model == 4) && (kptribmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Tribble (KP) model vales have not yet been calculated for this data set! Please run fitjpmodel first.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if (model == 1) {
	  mapfrommodel(imageloc, setname[currentset], xdim[currentset], ydim[currentset], settarget[currentset], usr_gmin, usr_gmax, fieldstrength, model, jpbestnorm[currentset], jpbestage[currentset], regionarray[currentset], inject, beamarea[currentset], redshift[currentset]);
	}
	if (model == 2) {
	  mapfrommodel(imageloc, setname[currentset], xdim[currentset], ydim[currentset], settarget[currentset], usr_gmin, usr_gmax, fieldstrength, model, kpbestnorm[currentset], kpbestage[currentset], regionarray[currentset], inject, beamarea[currentset], redshift[currentset]);
	}
	if (model == 3) {
	  mapfrommodel(imageloc, setname[currentset], xdim[currentset], ydim[currentset], settarget[currentset], usr_gmin, usr_gmax, fieldstrength, model, jptribbestnorm[currentset], jptribbestage[currentset], regionarray[currentset], inject, beamarea[currentset], redshift[currentset]);
	}
	if (model == 4) {
	  mapfrommodel(imageloc, setname[currentset], xdim[currentset], ydim[currentset], settarget[currentset], usr_gmin, usr_gmax, fieldstrength, model, kptribbestnorm[currentset], kptribbestage[currentset], regionarray[currentset], inject, beamarea[currentset], redshift[currentset]);
	}

      }

      else {
	fprintf(stderr,"\nError: No datasets have yet been loaded!\n\n");
      }

      break;


    case 105: // Subtract one map from another and output the difference to a fits file

      diffmap(casadata, imageloc, regionarray, regionsset, xdim, ydim, numdatasets, imgnum, setname, setreg, setbg);

      break;


    case 106: // Scale the raw flux of a data set

      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to scale the fluxes (-1 to escape, integers only): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if (currentset >= numdatasets) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}
      
	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	validentry = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Which maps should be scaled? (-1 to escape, 1 Single map, 2 Frequency range, 3 All maps): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  mapopt = strtol(cmdbuffer, &endptr, 10);

	  if ( (mapopt == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (mapopt < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (mapopt > 3) || (mapopt == 0) ) { // Check the set exists
	    printf("\nInvalid option. Please try again...\n\n");
	    continue;
	  }
	  else if (mapopt == 1) { // Check the set exists
	    scaleopt = 1;
	    validentry = 1;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if (mapopt == 1) {


	  printf("============================================\n\n");
	  printf(" Index  |  Frequency (Obs.)  |  Frequency (Rest)  \n\n");

	  // Loop through each of the sets
	  for (a=0; a<imgnum[currentset]; a++) {

	    printf("   %d       %.2e Hz             %.2e Hz\n", fluxorder[currentset][a], frequencyobs[currentset][fluxorder[currentset][a]], frequency[currentset][fluxorder[currentset][a]]);

	  }

	  printf("\n============================================\n");
    

	  minfluxscalefreq = -1.0;
	  maxfluxscalefreq = -1.0;

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Please enter the index of the map to be scaled (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    singlemap = strtol(cmdbuffer, &endptr, 10);

	    if ( (singlemap == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (singlemap >= imgnum[currentset] ) { // Check the map exists
	      printf("\nInvalid selection. Please try again...\n\n");
	      continue;
	    }

	    else if (singlemap < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }

	    else {
	      validentry = 1;
	    }
	  }

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }

	}
     
	else if (mapopt == 2) {


	  printf("============================================\n\n");
	  printf(" Index  |  Frequency (Obs.)  |  Frequency (Rest)  \n\n");

	  // Loop through each of the sets
	  for (a=0; a<imgnum[currentset]; a++) {

	    printf("   %d       %.2e Hz             %.2e Hz\n", fluxorder[currentset][a], frequencyobs[currentset][fluxorder[currentset][a]], frequency[currentset][fluxorder[currentset][a]]);

	  }

	  printf("\n============================================\n");
   

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Please enter the minimum frequency to be scaled (in Hz) (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    minfluxscalefreq = strtod(cmdbuffer, &endptr);

	    if ( (minfluxscalefreq == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (minfluxscalefreq < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }

	    else {
	      validentry = 1;
	    }
	  }

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Please enter the maximum frequency to be scaled (in Hz)? (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    maxfluxscalefreq = strtod(cmdbuffer, &endptr);

	    if ( (maxfluxscalefreq == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (maxfluxscalefreq < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }

	    else {
	      validentry = 1;
	    }
	  }

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }
	}
	else if (mapopt == 3) {
	  minfluxscalefreq = -1e27;
	  maxfluxscalefreq = 1e27;
	}
	else {
	  minfluxscalefreq = -1.0;
	  maxfluxscalefreq = -1.0;
	}
     

	if ( (mapopt == 2) || (mapopt == 3) ) {

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("How should the flux be scaled? (-1 to escape, 1 Single value, 2 Linear interpolation): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    scaleopt = strtol(cmdbuffer, &endptr, 10);

	    if ( (scaleopt == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (scaleopt < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else if ( (scaleopt > 2) || (scaleopt == 0) ) { // Check the set exists
	      printf("\nInvalid option. Please try again...\n\n");
	      continue;
	    }
	    else {
	      validentry = 1;
	    }
	  }

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }
	}
	else {
	  scaleopt = 1;
	}


	if (scaleopt == 1 ) {

	  factorfreq1 = -1e27;
	  factorfreq2 = 1e27;
	  factor1 = -1e27;
	  factor2 = 1e27;


	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Please enter the scaling value as a decimal value: ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    singlefactor = strtod(cmdbuffer, &endptr);

	    if ( (singlefactor == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }
	    else {
	      validentry = 1;
	    }
	  }

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }
	}

	else {

	  singlefactor = 1e27;

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Please enter the frequency of the first scaling factor (in Hz) (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    factorfreq1 = strtod(cmdbuffer, &endptr);

	    if ( (factorfreq1 == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (factorfreq1 < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else {
	      validentry = 1;
	    }
	  }

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Please enter the first scaling value as a decimal value: ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    factor1 = strtol(cmdbuffer, &endptr, 10);

	    if ( (factor1 == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }
	    else {
	      validentry = 1;
	    }
	  }

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Please enter the frequency of the second scaling factor (in Hz) (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    factorfreq2 = strtod(cmdbuffer, &endptr);

	    if ( (factorfreq2 == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    /*printf("Please enter the frequency of the second scaling factor (in Hz) (-1 to escape): ");
	    if (scanf("%f", &factorfreq2) == 0) {
	      printf("Invalid input, escaping command to protect the existing data...\n");
	      escape = 1;
	      break;
	    }*/

	    if (factorfreq2 < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }

	    else {
	      validentry = 1;
	    }
	  }

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }


	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Please enter the second scaling value as a decimal value: ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    factor2 = strtol(cmdbuffer, &endptr, 10);

	    if ( (factor2 == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }
	    else {
	      validentry = 1;
	    }
	  }

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }
	}

	validentry = 0;
	yesno = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Scaling of the raw fluxes will now be applied. This will not change the initial FITS maps, only the BRATS data. Are you sure you want to continue? (No [0], Yes [1]): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  yesno = strtol(cmdbuffer, &endptr, 10);

	  if ( (yesno == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if ( (yesno > 1) || (yesno < 0) ) {
	    printf("\nInvalid selection, please choose again...\n\n");
	    continue;
	  }

	  if (yesno == 0) {
	    printf("Exiting command...\n");
	    escape = 1;
	    break;
	  }

	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	scalefluxes(currentset, mapopt, minfluxscalefreq, maxfluxscalefreq, singlemap, scaleopt, factorfreq1, factorfreq2, factor1, factor2, singlefactor, flux, xdim[currentset], ydim[currentset], fluxorder[currentset], imgnum[currentset], frequencyobs[currentset], bgbuff, rms, beamarea[currentset]);


      }

      else {
	fprintf(stderr,"\nError: No datasets have yet been loaded!\n\n");
      }

      break;



    case 107: // Set shifting of images

      validentry = 0;

      printf("X-axis shift currently set to %d (DEFAULT: %d)\n", (0-xshift), (0-DEFAULTXSHIFT));

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value to shift the X-axis by relative to original position in pixels (integers only): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	xshift = strtol(cmdbuffer, &endptr, 10);

	if ( (xshift == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  xshift = 0 - xshift; // Swap the sign so it is in the correct format to be shifted
	  validentry = 1;
	}
      }

      validentry = 0;

      printf("Y-axis shift currently set to %d (DEFAULT: %d)\n", (0-yshift), (0-DEFAULTYSHIFT));

      while (validentry == 0) {

	cmdbuffer = readline("Enter a value to shift the Y-axis by relative to original position in pixels (integers only): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	yshift = strtol(cmdbuffer, &endptr, 10);

	if ( (yshift == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  yshift = 0 - yshift; // Swap the sign so it is in the correct format to be shifted
	  printf("\nImages will now be shifted by x: %d pixels y: %d pixels\n\n", (0-xshift), (0-yshift));
	  validentry = 1;
	}
      }

      break;


    case 108: // Fit a single injection model to the integrated flux

      validentry = 0;
      escape = 0;
      tmp_numdatapoints = 0;
      model = 0;

      
      if (export == 1) {

	// Open the directory that us currently set and check it exists
	dir = opendir(dataloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output data (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in dataloc ***\n\n", dataloc);
	  break;
	}
      }


      while (validentry == 0) {

	cmdbuffer = readline("Which model would you like to fit? (-1 to escape, 1 JP model, 2 KP model, 3 Tribble (JP), 4 Tribble (KP)): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	model = strtol(cmdbuffer, &endptr, 10);

	if ( (model == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}

	if (model < 0) {
	  printf("Escaping command...\n");
	  escape = 1;
	  break;
	}
	else if ( (model > NUMBEROFMODELS) ||  (model == 0) ) { // Check the set exists
	  printf("\nInvalid model number. Please try again...\n\n");
	  continue;
	}
	else {
	  validentry = 1;
	}
      }

      // Break the if they want to escape
      if (escape == 1) {
	break;
      }

      validentry = 0;
      
      while (validentry == 0) {

	cmdbuffer = readline("What is the redshift of your source? (-1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	tmp_redshift = strtod(cmdbuffer, &endptr);

	if ( (tmp_redshift == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (tmp_redshift < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else {
	    validentry = 1;
	  }
	}
      }
      
      // Break the if they want to escape
      if (escape == 1) {
	break;
      }


      validentry = 0;

      while (validentry == 0) {

	cmdbuffer = readline("How many flux values are there to be fitted to? (-1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	tmp_numdatapoints = strtol(cmdbuffer, &endptr, 10);

	if ( (tmp_numdatapoints == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {

	  if (tmp_numdatapoints < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (tmp_numdatapoints < 3) && (tmp_numdatapoints >= 0) ) {
	    printf("\nNumber of fluxes must be at least 3, please try again...\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}
      }
      
      // Break the if they want to escape
      if (escape == 1) {
	break;
      }
    
      validentry = 0;
    
      tmp_integfluxes = (float **)calloc(tmp_numdatapoints, sizeof(float *));
      tmp_integfreq = (float *)calloc(tmp_numdatapoints, sizeof(float));
      tmp_fluxerr = (float **)calloc(tmp_numdatapoints, sizeof(float *));

      for (i=0; i<=tmp_numdatapoints; i++) {
	tmp_integfluxes[i] = (float *)calloc(2, sizeof(float));
	tmp_fluxerr[i] = (float *)calloc(2, sizeof(float));
      }

      for (a=0; a<tmp_numdatapoints; a++) {

	validentry = 0;

	printf("*** Data point %d ***\n", a+1);

	while (validentry == 0) {

	  cmdbuffer = readline("Please enter a frequency (Hz, -1 to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  tmp_holdintfreq = strtod(cmdbuffer, &endptr);

	  if ( (tmp_holdintfreq == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }
	  else {
	    if (tmp_holdintfreq < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else {
	      tmp_integfreq[a] = tmp_holdintfreq;
	      validentry = 1;
	    }
	  }
	}


	if (escape == 1) {
	  break;
	}

	validentry = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Please enter a flux (Jy, -1 to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  tmp_holdintflux = strtod(cmdbuffer, &endptr);

	  if ( (tmp_holdintflux == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }
	  else {
	    if (tmp_holdintflux < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else {
	      tmp_integfluxes[a][1] = tmp_holdintflux;
	      validentry = 1;
	    }
	  }
	}

	if (escape == 1) {
	  break;
	}
      
	validentry = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Please enter the total error as a decimal (-1 to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  tmp_holdfluxerr = strtod(cmdbuffer, &endptr);

	  if ( (tmp_holdfluxerr == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }
	  else {
	    if (tmp_holdfluxerr < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else {
	      tmp_fluxerr[a][1] = tmp_holdfluxerr;
	      validentry = 1;
	    }
	  }
	}

	if (escape == 1) {
	  break;
	}
      }
    
      if (escape == 1) {
	break;
      }
      
      // Get how the filename should be appended
      validentry = 0;

      while (validentry == 0) {

	cmdbuffer = readline("Please enter an identifier for the file output e.g. entering new will result in the file name IntergratedFit_Mod1_new.txt for JP model data (esc to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	strcpy(tmp_filename, cmdbuffer);

	if ( (tmp_filename == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else if ((strstr(tmp_filename, "esc") != NULL) && (strlen(tmp_filename) == 3)) {
	  printf("Escaping command...\n");
	  escape = 1;
	  break;
	}
	else {
	  validentry = 1;
	}
      }
      
      // Break the second loop if they want to escape
      if (escape == 1) {
	break;
      }

      tmp_fluxorder =  (int *)calloc((tmp_numdatapoints+1),sizeof(int));

      for (i=0; i < tmp_numdatapoints; i++) {
	if (i == 0) {
	  tmp_curminfreq = 1e-32;
	}
	else {
	  tmp_curminfreq = tmp_integfreq[tmp_fluxorder[i-1]];
	}

	tmp_curmaxfreq = 1e+32;
      
	for (a=0; a < tmp_numdatapoints; a++) {
	  if (a==0) {
	    tmp_fluxorder[i] = a;
	  }
	  if ( (tmp_integfreq[a] > tmp_curminfreq) && (tmp_integfreq[a] < tmp_curmaxfreq) ) {
	    tmp_curmaxfreq = tmp_integfreq[a];
	    tmp_fluxorder[i] = a;
	  }
	}
      }

      //setmodelres[currentset][model] = modelres;

      tmp_chisq = (float *)calloc(2, sizeof(float));
      tmp_bage = (float *)calloc(2, sizeof(float));
      tmp_modalpha  = (float *)calloc(2, sizeof(float));
      tmp_modflux = (float **)calloc(2, sizeof(float *));
      tmp_modflux[1] = (float *)calloc(modelres+1, sizeof(float));
      tmp_bnorm = (float *)calloc(2, sizeof(float));

      // Using a temp array here to allow all models without the need for dummy arrays
      tmp_ageerrorsplus = (float *)calloc(2, sizeof(float));
      tmp_ageerrorsminus = (float *)calloc(2, sizeof(float));

      /*
	if (redshiftallocated == 0) {
	redshift = (float *)calloc(MAXDATASETS,sizeof(float));
	redshift[0] = 0.2;
	redshiftallocated = 1;
	}*/

      spectralageingmodels(tmp_integfluxes, tmp_numdatapoints, tmp_fluxorder, tmp_integfreq, 1, tmp_chisq, tmp_bage, tmp_fluxerr, tmp_modalpha, tmp_modflux, usr_gmin, usr_gmax, minmyears, myears, ageresolution, levels, tmp_bnorm, modelres, printresults, inject, fieldstrength, model, printreject, tmp_redshift, tmp_ageerrorsplus, tmp_ageerrorsminus, suppresscdf);


      // Integrated errors need weighting by the beam area. To do if anyone needs it...

      // Hack to output model values to a text file. This need reworking for adaptive names.

      //tmp_modfreq = (float *)calloc(modelres+1, sizeof(float));

      /*
	time( &currenttime );
	time_struct = localtime ( &currenttime );

	strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);
      */

      sprintf(tmp_filename2,"%s/IntergratedFit_Mod%d_%s.txt", dataloc, model, tmp_filename);

      printf("\nExporting model values as %s\n", tmp_filename2);


      //sprintf(tmp_filename,"%s/intergratedmodelvalues.txt", dataloc);


      filestream = fopen(tmp_filename2,"w");
      for (t=0; t<=modelres; t++) {

	tmp_modfreq = log10(tmp_integfreq[tmp_fluxorder[0]] + ( ( (tmp_integfreq[tmp_fluxorder[tmp_numdatapoints-1]] - tmp_integfreq[tmp_fluxorder[0]]) / modelres ) * t ) );

	//printf("%.6e \n", tmp_modflux[1][t]);

	fprintf(filestream, "%.6e %.6e", pow(10, tmp_modfreq), tmp_modflux[1][t]);
	fprintf(filestream, "\n");
 
      }
      fclose(filestream);
	  
      // Free up the memory now incase we want to run it again
      free(tmp_fluxorder);
      free(tmp_chisq);
      free(tmp_bage);
      free(tmp_modalpha);
      free(tmp_modflux);
      free(tmp_bnorm);
      free(tmp_ageerrorsplus);
      free(tmp_ageerrorsminus);


      break;


    case 109: // Make an error map for a given dataset

      if (export == 1) {

	// Open the directory that is currently set and check it exists
	dir = opendir(imageloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	  break;
	}
      }

      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to map the spectral age (-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    printf("\nThis value is greater than the number of data sets that exists! Please try again...\n\n");
	    continue;
	  }

	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    continue;
	  }

	  else if ( (currentset != 888) && (jpmodelset[currentset] != 1) && (kpmodelset[currentset] != 1) && (jptribmodelset[currentset] != 1) && (kptribmodelset[currentset] != 1)  ) { // Check the curvature has been set
	    fprintf(stderr,"\nError: A model has not yet been set for this data set! Please run a model fitting command first.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}


	// Ask which model to map and check its validity
	validentry = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Which model would you like to map? (-1 to escape, 1 for JP model, 2 for KP model, 3 Tribble (JP), 4 Tribble (KP)): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  model = strtol(cmdbuffer, &endptr, 10);

	  if ( (model == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (model < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if ( (model > NUMBEROFMODELS) ||  (model == 0) ) { // Check the set exists
	    printf("\nInvalid model number. Please try again...\n\n");
	    continue;
	  }

	  else if ( (currentset != 888) && (model == 1) && (jpmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: JP model values have not yet been calculated for this data set! Please run fitjpmodel first.\n\n");
	    continue;
	  }

	  else if ( (currentset != 888) && (model == 2) && (kpmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: KP model values have not yet been calculated for this data set! Please run fitkpmodel first.\n\n");
	    continue;
	  }

	  else if ( (currentset != 888) && (model == 3) && (jptribmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Tribble (JP) model values have not yet been calculated for this data set! Please run fitjptribble first.\n\n");
	    continue;
	  }

	  else if ( (currentset != 888) && (model == 4) && (kptribmodelset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Tribble (KP) model values have not yet been calculated for this data set! Please run fitkptribble first.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}
      
	if ( (contours == 1) && (currentset != 888) ) {

	  if ( (exportfits != 1)  || (export != 1) ) {

	    validentry = 0;

	    printf("============================================\n\n");
	    printf(" Index  |  Frequency (Obs.)  |  Frequency (Rest)  \n\n");

	    // Loop through each of the sets
	    for (a=0; a<imgnum[currentset]; a++) {

	      printf("   %d       %.2e Hz             %.2e Hz\n", fluxorder[currentset][a], frequencyobs[currentset][fluxorder[currentset][a]], frequency[currentset][fluxorder[currentset][a]]);

	    }

	    printf("\n============================================\n");



	    while (validentry == 0) {

	      cmdbuffer = readline("Select a map frequency for the contours (-1 to escape): ");
	      if (cmdbuffer[0] != 0) {
		add_history(cmdbuffer);
	      }

	      chomp(cmdbuffer);

	      currentmap = strtol(cmdbuffer, &endptr, 10);

	      if ( (currentmap == 0) && (cmdbuffer == endptr) ) {
		printf("\nInvalid input, please try again...\n\n");
		continue;
	      }

	      if (currentmap < 0) {
		printf("Escaping command...\n");
		escape = 1;
		break;
	      }
	      else if (currentmap >= imgnum[currentset]) { // Check the set exists
		printf("\nInvalid selection, please choose again...\n\n");
		continue;
	      }
	      else {
		validentry = 1;
	      }
	    }
	  }
	}
	else if ( (contours == 1) && (currentset == 888) ) { // Multimaps always uses map 0 as reference. This can easily be changed if needed
	  currentmap = 0;
	}
	else { //If we are not using contours, make sure its still a sensible value
	  currentmap = 0;
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}
      
	if ( (currentset == 888) && (export == 0) ) {
	  looplimit = numdatasets;
	  extwincontrol = 1;
	  sprintf(output,"/xs");
	  cpgopen(output);
	}
	else if ( (currentset == 888) && (export == 1) ) {
	  looplimit = numdatasets;
	  extwincontrol = 0;
	}
	else {
	  looplimit = 1;
	  extwincontrol = 0;
	}
      
	validentry = 0;

	while (validentry == 0) {

	  errorpm = 0;

	  cmdbuffer = readline("Use positive or negative errors? (0 for positive, 1 for negative, -1 to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  errorpm = strtol(cmdbuffer, &endptr, 10);

	  if ( (errorpm == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (errorpm < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if (errorpm > 1) { // Make sure it is 0 or 1
	    printf("\nInvalid selection, please choose again...\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}


	for (j=0; j<looplimit; j++) {

	  if (looplimit > 1) {
	    currentset = j;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = j;
	  }

	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets

	  if (regionsset[currentset] != 1) {
	    printf("\nThe regions for data set %d have not been set! Skipping this set...\n", j);
	    continue;	  
	  }

	  if ( (model == 1) && (jpmodelset[currentset] != 1) ) {
	    printf("\nThe JP model values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  else if ( (model == 2) && (kpmodelset[currentset] != 1) ) {
	    printf("\nThe KP model values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  else if ( (model == 3) && (jptribmodelset[currentset] != 1) ) {
	    printf("\nThe Tribble (JP) model values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  else if ( (model == 4) && (kptribmodelset[currentset] != 1) ) {
	    printf("\nThe Tribble (KP) model values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  time( &currenttime );
	  time_struct = localtime ( &currenttime );
	  
	  strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	  if ( ( export == 1 ) && (exportfits != 1) ) {
	    
	    if (model == 1) {
	      sprintf(modelbuff,"JP");
	    }
	    else if (model == 2) {
	      sprintf(modelbuff,"KP");
	    }
	    else if (model == 3) {
	      sprintf(modelbuff,"TribbleJP");
	    }
	    else if (model == 4) {
	      sprintf(modelbuff,"TribbleKP");
	    }
	    else {
	      sprintf(modelbuff,"Unknown");
	    }

	    if (errorpm == 0) {
	      sprintf(posnegbuff,"Pos");
	    }
	    else if (model == 1) {
	      sprintf(posnegbuff,"Neg");
	    }
	    else {
	      sprintf(posnegbuff,"Unknown");
	    }

	    if (contours == 1) {
		sprintf(output,"%s/%s_%sErrorMap_Contours_%s_%s.%s/%s", imageloc, settarget[currentset], posnegbuff, modelbuff, timebuff, imagetype, imagetype);

	      }
	      else {
		sprintf(output,"%s/%s_%sErrorMap_%s_%s.%s/%s", imageloc, settarget[currentset], posnegbuff, modelbuff, timebuff, imagetype, imagetype);
	      }



	    printf("\nExporting as %s\n", output);
	  }

	  else if ( ( export == 1 ) && (exportfits == 1) ) {

	    if (model == 1) {
	      if (errorpm == 0) {
		sprintf(output,"%s/%s_PosErrorMap_JP_%s.fits", imageloc, settarget[currentset], timebuff);
	      }
	      else if (errorpm == 1) {
		sprintf(output,"%s/%s_NegErrorMap_JP_%s.fits", imageloc, settarget[currentset], timebuff);
	      }
	      else {
		sprintf(output,"%s/%s_UnknownErrorMap_JP_%s.fits", imageloc, settarget[currentset], timebuff);
	      }
	    }

	    else if (model == 2) {
	      if (errorpm == 0) {
		sprintf(output,"%s/%s_PosErrorMap_KP_%s.fits", imageloc, settarget[currentset], timebuff);
	      }
	      else if (errorpm == 1) {
		sprintf(output,"%s/%s_NegErrorMap_KP_%s.fits", imageloc, settarget[currentset], timebuff);
	      }
	      else {
		sprintf(output,"%s/%s_UnknownErrorMap_KP_%s.fits", imageloc, settarget[currentset], timebuff);
	      }
	    }

	    else if (model == 3) {
	      if (errorpm == 0) {
		sprintf(output,"%s/%s_PosErrorMap_TribbleJP_%s.fits", imageloc, settarget[currentset], timebuff);
	      }
	      else if (errorpm == 1) {
		sprintf(output,"%s/%s_NegErrorMap_TribbleJP_%s.fits", imageloc, settarget[currentset], timebuff);
	      }
	      else {
		sprintf(output,"%s/%s_UnknownErrorMap_TribbleJP_%s.fits", imageloc, settarget[currentset], timebuff);
	      }
	    }

	    else if (model == 4) {
	      if (errorpm == 0) {
		sprintf(output,"%s/%s_PosErrorMap_TribbleKP_%s.fits", imageloc, settarget[currentset], timebuff);
	      }
	      else if (errorpm == 1) {
		sprintf(output,"%s/%s_NegErrorMap_TribbleKP_%s.fits", imageloc, settarget[currentset], timebuff);
	      }
	      else {
		sprintf(output,"%s/%s_UnknownErrorMap_TribbleKP_%s.fits", imageloc, settarget[currentset], timebuff);
	      }
	    }
	    else {
	      sprintf(output,"%s/%s_ErrorMap_UnknownModel_%s.fits", imageloc, settarget[currentset], timebuff);
	    }

	    printf("\nExporting as %s\n", output);

	  }
	  else {
	    sprintf(output,"/xs");
	  }

	
	  if (model == 1) {
      
	    if (errorpm == 0) {
	      errormap(jpageerrorsplus[currentset], regionarray[currentset], flux[currentset][currentmap], xdim[currentset], ydim[currentset], settarget[currentset], border[currentset], zoom, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, firstcontour, contourlogs, output, jpinject[currentset], jpfield[currentset], paramlabels, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, xshift, yshift, errorpm, export, exportfits, largetxt);

	    }

	    else {
	      errormap(jpageerrorsminus[currentset], regionarray[currentset], flux[currentset][currentmap], xdim[currentset], ydim[currentset], settarget[currentset], border[currentset], zoom, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, firstcontour, contourlogs, output, jpinject[currentset], jpfield[currentset], paramlabels, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, xshift, yshift, errorpm, export, exportfits, largetxt);
	    }

	  }

	  else if (model == 2) {

	    if (errorpm == 0) {
	      errormap(kpageerrorsplus[currentset], regionarray[currentset], flux[currentset][currentmap], xdim[currentset], ydim[currentset], settarget[currentset], border[currentset], zoom, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, firstcontour, contourlogs, output, kpinject[currentset], kpfield[currentset], paramlabels, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, xshift, yshift, errorpm, export, exportfits, largetxt);

	    }

	    else {
	      errormap(kpageerrorsminus[currentset], regionarray[currentset], flux[currentset][currentmap], xdim[currentset], ydim[currentset], settarget[currentset], border[currentset], zoom, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, firstcontour, contourlogs, output, kpinject[currentset], kpfield[currentset], paramlabels, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, xshift, yshift, errorpm, export, exportfits, largetxt);
	    }

	  }

	  else if (model == 3) {
      
	    if (errorpm == 0) {
	      errormap(jptribbleageerrorsplus[currentset], regionarray[currentset], flux[currentset][currentmap], xdim[currentset], ydim[currentset], settarget[currentset], border[currentset], zoom, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, firstcontour, contourlogs, output, jptribinject[currentset], jptribfield[currentset], paramlabels, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, xshift, yshift, errorpm, export, exportfits, largetxt);

	    }

	    else {
	      errormap(jptribbleageerrorsminus[currentset], regionarray[currentset], flux[currentset][currentmap], xdim[currentset], ydim[currentset], settarget[currentset], border[currentset], zoom, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, firstcontour, contourlogs, output, jptribinject[currentset], jptribfield[currentset], paramlabels, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, xshift, yshift, errorpm, export, exportfits, largetxt);
	    }

	  }

	  else if (model == 4) {
	    if (errorpm == 0) {
	      printf("*** Error: KP Tribble no longer in use ***\n");
	    }
	    else {
	      printf("*** Error: KP Tribble no longer in use ***\n");
	    }
	  }

	  else {
	    fprintf(stderr,"\nError: Model number is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");
	    break;
	  }

	}

	if (extwincontrol == 1) {
	  cpgclos();
	  extwincontrol = 0;
	}

      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;


    case 110: // Exports all information associated with a given dataset in .brats format

      if (export  == 1) {

	// Open the directory that us currently set and check it exists
	dir = opendir(dataloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output data (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in dataloc ***\n\n", dataloc);
	  break;
	}
      }

      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set which you would like to export(-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }

	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	// Set the time before the loop to give similar names (sets are distigished by set number)
	time( &currenttime );
	time_struct = localtime ( &currenttime );
	  
	strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

      
	// Assign some place holder memory if fitting hasnt been done

	checkthejpmemory = 0;
	/*	
		for (i=0; i<MAXDATASETS; i++) {
		if (jpfitmemset[i] != 0) {
		checkthejpmemory = 1;
		}
		}*/
	

	if (jpmodelmemoryset == 0) {
	  jpbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jpbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jpchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jpmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jpmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  jpageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jpageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));

	  checkthejpmemory = 1;
	    
	}

	  
	checkthekpmemory = 0;
      
	/*
	  for (i=0; i<MAXDATASETS; i++) {
	  if (kpfitmemset[i] != 0) {
	  checkthekpmemory = 1;
	  }
	  }*/

	if (kpmodelmemoryset == 0) {
	  kpbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	  kpbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	  kpchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	  kpmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	  kpmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  kpageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  kpageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));

	  checkthekpmemory = 1;
	}

	  
	checkthejptribmemory = 0;
      
	/*
	  for (i=0; i<MAXDATASETS; i++) {
	  if (jptribfitmemset[i] != 0) {
	  checkthejptribmemory = 1;
	  }
	  }*/

	if (jptribmodelmemoryset == 0) {
	  jptribbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jptribbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jptribchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jptribmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jptribmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  jptribbleageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jptribbleageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));

	  checkthejptribmemory = 1;

	}
	  

	if (currentset == 888) {
	  looplimit = numdatasets;
	}
	else {
	  looplimit = 1;
	}

	for (j=0; j<looplimit; j++) {

	  if (looplimit > 1) {
	    currentset = j;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = j;
	  }
	  
	  if (j == 0) {
	    // Make the file name
	    validentry = 0;
	    escape = 0;

	    while (validentry == 0) {
    
	      cmdbuffer = readline("Please enter an identifier for the file output [Target_TimeStamp_SetNumber_YourEntry.brats](esc to escape): ");
	      if (cmdbuffer[0] != 0) {
		add_history(cmdbuffer);
	      }

	      chomp(cmdbuffer);

	      strcpy(tmp_filename, cmdbuffer);

	      if ( (tmp_filename == 0) && (cmdbuffer == endptr) ) {
		printf("\nInvalid input, please try again...\n\n");
		continue;
	      }
	      else if ((strstr(tmp_filename, "esc") != NULL) && (strlen(tmp_filename) == 3)) {
		printf("Escaping command...\n");
		escape = 1;
		break;
	      }
	      else {
		validentry = 1;
	      }
	    }
	  
	    // Break the second loop if they want to escape
	    if (escape == 1) {
	      break;
	    }
	  }
	  
	  sprintf(filename,"%s/%s_%s_%d_%s.brats", dataloc, settarget[currentset], timebuff, currentset, tmp_filename);


	  fullexport(filename, exportcompression, settarget[currentset], imgnum[currentset], redshift[currentset], xdim[currentset], ydim[currentset], bmaj[currentset], bmin[currentset], bpa[currentset], beamarea[currentset], ra[currentset], dec[currentset], eq[currentset], delt1[currentset], delt2[currentset], crpix1[currentset], crpix2[currentset], crota1[currentset], crota2[currentset], cellsize[currentset], setname[currentset], setbg[currentset], setreg[currentset], posmap, dataloc, regionsset[currentset], jpmodelset[currentset], kpmodelset[currentset], jptribmodelset[currentset], mininjectset[currentset], frequency[currentset], frequencyobs[currentset], fluxorder[currentset], rms[currentset], flux[currentset], minflux[currentset], maxflux[currentset], bgbuff[currentset], fluxcalerror[currentset], setaveraged[currentset], regnumber[currentset], regionlocx[currentset], regionlocy[currentset], regminflux[currentset], regmaxflux[currentset], fluxerror[currentset], regionarray[currentset], regflux[currentset], mininjectstore[currentset], maxinjectstore[currentset], inj_sumchisquared[currentset], setmodelres[currentset], jpfield[currentset], jpinject[currentset], jpbestage[currentset], jpbestnorm[currentset], jpchisquared[currentset], jpmodelalpha[currentset], jpmodelflux[currentset], jpageerrorsplus[currentset], jpageerrorsminus[currentset], kpfield[currentset], kpinject[currentset], kpbestage[currentset], kpbestnorm[currentset], kpchisquared[currentset], kpmodelalpha[currentset], kpmodelflux[currentset], kpageerrorsplus[currentset], kpageerrorsminus[currentset], jptribfield[currentset], jptribinject[currentset], jptribbestage[currentset], jptribbestnorm[currentset], jptribchisquared[currentset], jptribmodelalpha[currentset], jptribmodelflux[currentset], jptribbleageerrorsplus[currentset], jptribbleageerrorsminus[currentset], regionsize[currentset], inj_regbestinject[currentset], inj_regchisquared[currentset]);

	}

	// Clear any temporary memory assigned
	
	if ( (jpmodelmemoryset == 0) && (checkthejpmemory == 1) ) {
	  free(jpbestage);
	  free(jpbestnorm);
	  free(jpchisquared);
	  free(jpmodelalpha);
	  free(jpmodelflux);
	  free(jpageerrorsplus);
	  free(jpageerrorsminus);
	  checkthejpmemory = 0;
	}
	
	if ( (kpmodelmemoryset == 0) && (checkthekpmemory == 1) ) {
	  
	  free(kpbestage);
	  free(kpbestnorm);
	  free(kpchisquared);
	  free(kpmodelalpha);
	  free(kpmodelflux);
	  free(kpageerrorsplus);
	  free(kpageerrorsminus);
	  checkthekpmemory = 0;
	}
	
	
	if ( (jptribmodelmemoryset == 0) && (checkthejptribmemory == 1) ) {
	  free(jptribbestage);
	  free(jptribbestnorm);
	  free(jptribchisquared);
	  free(jptribmodelalpha);
	  free(jptribmodelflux);
	  free(jptribbleageerrorsplus);
	  free(jptribbleageerrorsminus);
	  checkthejptribmemory = 0;
	}
	
      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;


    case 111: // Turn export compression on and off

      if (exportcompression == 1) {
	exportcompression = 0;
	printf("BRATS will now export datasets in their uncompressed state\n");
      }
      else if (exportcompression == 0) {
	exportcompression = 1;
	printf("BRATS will now export datasets in their compressed state\n");
      }
      else {
	printf("\nThe value of exportcompression is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", exportcompression);
      }

      break;


    case 112: // Import a full dataset

      exportversion = 0; // Reset the export version
      
      if (arraysallocated == 0) {

	xdim = (int *)calloc(MAXDATASETS,sizeof(int));
	ydim = (int *)calloc(MAXDATASETS,sizeof(int));
	bmaj = (float *)calloc(MAXDATASETS,sizeof(float));
	bmin = (float *)calloc(MAXDATASETS,sizeof(float));
	bpa = (float *)calloc(MAXDATASETS,sizeof(float));
	beamarea = (float *)calloc(MAXDATASETS,sizeof(float));
	imgnum = (int *)calloc(MAXDATASETS,sizeof(int));
	maxflux = (float **)calloc(MAXDATASETS,sizeof(float *));
	minflux = (float **)calloc(MAXDATASETS,sizeof(float *));

	rms = (float **)calloc(MAXDATASETS,sizeof(float *));
	flux = (float ****)calloc(MAXDATASETS,sizeof(float ***));
	fluxorder =  (int **)calloc(MAXDATASETS,sizeof(int *));
	bgbuff = (float ***)calloc(MAXDATASETS, sizeof(float **));

	frequency = (float **)calloc(MAXDATASETS,sizeof(float *));
	frequencyobs = (float **)calloc(MAXDATASETS,sizeof(float *));
	ra = (float **)calloc(MAXDATASETS,sizeof(float *));
	dec = (float **)calloc(MAXDATASETS,sizeof(float *));
	eq = (float **)calloc(MAXDATASETS,sizeof(float *));
	delt1 = (float **)calloc(MAXDATASETS,sizeof(float *));
	delt2 = (float **)calloc(MAXDATASETS,sizeof(float *));
	crpix1 = (float **)calloc(MAXDATASETS,sizeof(float *));
	crpix2 = (float **)calloc(MAXDATASETS,sizeof(float *));
	crota1 = (float **)calloc(MAXDATASETS,sizeof(float *));
	crota2 = (float **)calloc(MAXDATASETS,sizeof(float *));
	cellsize = (float *)calloc(MAXDATASETS,sizeof(float));
	redshift = (float *)calloc(MAXDATASETS,sizeof(float));
	//imgbuff = (float ***)calloc((MAXDATASETS * MAXNUMMAPS),sizeof(float **));
	//bgbuff = (float ***)calloc((MAXDATASETS * MAXNUMMAPS),sizeof(float **));

	arraysallocated = 1;
      }
      
      else if (numdatasets >= MAXDATASETS) {
	fprintf(stderr,"Error: Maximum number of datasets has been reached. Please delete or replace one!\n");
	break;
      }

      currentset = numdatasets; // This allows later addition or modification and deleting arrays

      // Get set the .brats file to be loaded and read the header

      validentry = 0;
      escape = 0;

      while (validentry == 0) {

	strcpy(importfilename, "");

	cmdbuffer = readline("Enter the location of a file to import (ls to list current directory, esc to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	strcpy(importfilename, cmdbuffer);

	chomp(importfilename);

	// Check a directory name has been enetered isn't over the maximum character limit
	if (strlen(importfilename) >= MAXCMDLENGTH-1) {
	  fprintf(stderr,"\n*** Error: Folder path must be less than %d characters long! ***\n\n", MAXCMDLENGTH);
	  break;
	}
	else if (strlen(importfilename) == 0) {
	  continue;
	}
	else if ((strstr(importfilename, "esc") != NULL) && (strlen(importfilename) == 3)) {
	  printf("Escaping command...\n");
	  escape = 1;
	  break;
	}
	else if ((strstr(importfilename, "ls") != NULL) && (strlen(importfilename) == 2)) {
	  syscom("ls");
	}
	else if( access(importfilename, R_OK) != -1 ) {
	  validentry = 1;
	}
	else {
	  fprintf(stderr,"\n *** Unable to read %s. Check the file exists and the read permissions are set correctly *** \n\n", importfilename); 
	}
      }

      // Break the second loop if they want to escape
      if (escape == 1) {
	break;
      }
      
      // Check the file exists and that it is in a compatable format
      char line_storage[1024], buffer[1024], compstr[100], readstr1[100], readstr2[100], readstr3[100], readstr4[100], readstr5[100], readstr6[100], readstr7[100], readstr8[100], readstr9[100], readstr10[100], readstr11[100], readstr12[100], readstr13[100];

      strfound = 0;
      line_num = 1;
      foundheader = 0;	
      
      importfile = fopen(importfilename, "r");
      
      if (importfile != NULL) {
	
	while( fgets(line_storage, sizeof(line_storage), importfile) != NULL )  {
	  strfound = 0;
	  sscanf(line_storage,"%s", buffer);
	    
	  if (foundheader == 3) {

	    sscanf(buffer,"%100[^,]", compstr);

	    // printf("buffer: %s\n", buffer);
	    // printf("compstr: %s\n", compstr);

	    if (strcmp(compstr,"1HEAD1") != 0) {
	      printf("\n*** ERROR: Unable to read or locate the header data. The data file may be corrupt! ***\n\n");
	      fclose(importfile);
	      break;
	    }
	    else {
	      printf("Reading the header...\n");
	      //sscanf(buffer,"%*s100[^,],%s100[^,],%f100[^,],%d100[^,],%d100[^,],%d100[^,],%f100[^,],%f100[^,],%f100[^,],%f100[^,],%f100[^,],%s100[^,],%s100[^,],%s", settarget[currentset], &redshift[currentset], &imgnum[currentset], &xdim[currentset], &ydim[currentset], &bmaj[currentset], &bmin[currentset], &bpa[currentset], &beamarea[currentset], &cellsize[currentset], setname[currentset], setreg[currentset], setbg[currentset]);

	      sscanf(buffer,"%*100[^,],%100[^,],%100[^,],%100[^,],%100[^,],%100[^,],%100[^,],%100[^,],%100[^,],%100[^,],%100[^,],%100[^,],%100[^,],%s", readstr1, readstr2, readstr3, readstr4, readstr5, readstr6, readstr7, readstr8, readstr9, readstr10, readstr11, readstr12, readstr13);

	      //printf("READ STRINGS: %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",  readstr1, readstr2, readstr3, readstr4, readstr5, readstr6, readstr7, readstr8, readstr9, readstr10, readstr11, readstr12, readstr13);

	      // Assign memory dynamically for settarget
	      settarget[currentset] = (char *)malloc(sizeof(tmp_maptarget)+1);

	      // Copy over the read strings to their apppropriate variable and format (there is probably a much cleaner way of doing all this!!!)

	      strcpy(settarget[currentset], readstr1);
	      redshift[currentset] = atof(readstr2);
	      imgnum[currentset] = atoi(readstr3);
	      xdim[currentset] = atoi(readstr4);
	      ydim[currentset] = atoi(readstr5);
	      bmaj[currentset] = atof(readstr6);
	      bmin[currentset] = atof(readstr7);
	      bpa[currentset] = atof(readstr8);
	      beamarea[currentset] = atof(readstr9);
	      cellsize[currentset] = atof(readstr10);
	      strcpy(setname[currentset], readstr11);
	      strcpy(setreg[currentset], readstr12);
	      strcpy(setbg[currentset], readstr13);

	      //printf("CONVERTED: %s,%.6f,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%s,%s\n", settarget[currentset], redshift[currentset], imgnum[currentset], xdim[currentset], ydim[currentset], bmaj[currentset], bmin[currentset], bpa[currentset], beamarea[currentset], cellsize[currentset], setname[currentset], setreg[currentset], setbg[currentset]);
	      fclose(importfile);
	      printf("Header imported successfully\n");
	      break;
	    }
	  }
	    
	  if (foundheader == 2) {
	    if ( (strcmp(buffer,"Compact") != 0) && (strcmp(buffer,"Full") != 0) ) {
	      printf("\n*** ERROR: File is not in a known BRATS format (format type missing or incorrect). The data file may be corrupt! ***\n\n");
	      fclose(importfile);
	      break;
	    }
	    else {
	      printf("Data format is of type %s\n", buffer);

	      if (strcmp(buffer,"Compact") == 0) {
		importcompression = 1;
	      }
	      else if (strcmp(buffer,"Full") == 0 ) {
		importcompression = 2;
	      }

	      foundheader++;
	    }
	  }
	    
	  if (foundheader == 1) {
	    if ( (strcmp(buffer,"BRATSV1") != 0) && (strcmp(buffer,"BRATSV2") != 0) ) {
	      printf("\n*** ERROR: File is not in a known BRATS format (version type missing or incorrect). The data file may be corrupt! ***\n\n");
	      fclose(importfile);
	      break;
	    }
	    else {
	      if (strcmp(buffer,"BRATSV1") == 0) {
		exportversion = 1;
	      }
	      else if (strcmp(buffer,"BRATSV2") == 0) {
		exportversion = 2;
	      }
	      else { // Catch all for completeness, but this should never be triggered
		printf("\n *** Warning: Unknown export version. Defaulting to the latest type. ***\n\n");
		exportversion = 2;
	      }
	      printf("Data is of type %s\n", buffer);
	      foundheader++;
	    }
	  }
	    
	  if(strcmp(buffer,"###GENHEAD###") == 0)  {
	    strfound = 1;
	  }
	  if(strfound == 1) {
	    printf("GENHEAD match found on line %d\n", line_num);
	    foundheader = 1;
	  }

	  if (strcmp(buffer,"###ENDOFIMPORTFILE###") == 0)  {
	    printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EM1001) ***\n\n");
	    fclose(importfile);
	    break;
	  }
	    
	  line_num++;
	}
	
      }
      // If we can't find or access the file, kick back to the command prompt
      else {
	fprintf(stderr,"*** Error: Unable to locate or access the .brats file requested ***\n\n");
	break;
      }

      // If we find nothing that matches, the kick out to command promt
      if (foundheader == 0) {
	printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EM1002) ***\n\n");
	fclose(importfile);
	break;
      }

      // Import the RA and DEC (all these functions can be found in fullimport.h)
      if(importcoords(importfilename, currentset, ra, dec, imgnum[currentset]) != 0) {
	break; // Escape if there is a problem found in the function
      }

      // Import the map properties is we are greater than version 2, else calculate them or add dummy entries
      if (exportversion >= 2) {
	if(importmapprops(importfilename, currentset, eq, delt1, delt2, crpix1, crpix2, crota1, crota2, imgnum[currentset]) != 0) {
	  break; // Escape if there is a problem found in the function
	}
      }
      else if (exportversion == 1) {
	printf("\n*** Warning: The file to be imported is in an older format. This will be automatically converted but the results should be checked carefully. ***\n\n");
	if(convertmapprops(importfilename, currentset, exportversion, eq, delt1, delt2, crpix1, crpix2, crota1, crota2, imgnum[currentset], cellsize[currentset], xdim[currentset], ydim[currentset]) != 0) {
	  break; // Escape if there is a problem found in the function
	}
      }
      else { // Asuumer we are using the latest version as before for a catch all
	if(importmapprops(importfilename, currentset, eq, delt1, delt2, crpix1, crpix2, crota1, crota2, imgnum[currentset]) != 0) {
	  break; // Escape if there is a problem found in the function
	}
      }

      // See what has been run in SETHEAD
      importfile = fopen(importfilename, "r");

      if (importfile != NULL) {
	
	headerloop = 0;
	line_num = 1;
	foundheader = 0;

	while( fgets(line_storage, sizeof(line_storage), importfile) != NULL )  {

	  sscanf(line_storage,"%s", buffer);
	  // printf("line_storage: %s buffer: %s\n", line_storage, buffer);

	  
	  if (strcmp(buffer,"###ENDOFIMPORTFILE###") == 0)  {
	    printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EM1003) ***\n\n");
	    fclose(importfile);
	    break;
	  }

	  if (foundheader == 1) {
	    if ( (strcmp(buffer,"###END###") == 0) && (headerloop != expectedheaders) ) {
	      printf("\n*** ERROR: Finished reading SETHEAD before the expected number of iterations had been reached! The data may be corrupt ***\n\n");
	      fclose(importfile);
	      break;
	    }
	    else if ( (strcmp(buffer,"###END###") == 0) && (headerloop == expectedheaders) ) {
	      printf("SETHEAD imported successfully\n");
	      fclose(importfile);
	      break;
	    }
	    else {
	      if (headerloop == 0) {
		sscanf(buffer,"%d", &regionsset[currentset]);
	      }
	      else if (headerloop == 1) {
		sscanf(buffer,"%d", &jpmodelset[currentset]);
	      }
	      else if (headerloop == 2) {
		sscanf(buffer,"%d", &kpmodelset[currentset]);
	      }
	      else if (headerloop == 3) {
		sscanf(buffer,"%d", &jptribmodelset[currentset]);
	      }
	      else if (headerloop == 4) {
		sscanf(buffer,"%d", &mininjectset[currentset][1]);
	      }
	      else if (headerloop == 5) {
		sscanf(buffer,"%d", &mininjectset[currentset][2]);
	      }
	      else if (headerloop == 6) {
		sscanf(buffer,"%d", &mininjectset[currentset][3]);
	      }

	      headerloop++;
	    }	    
	  }

	  if (strcmp(buffer,"###SETHEAD###") == 0)  {
	    printf("SETHEAD found on line %d\n", line_num);
	    foundheader = 1;
	  }

	  line_num++;
	}

	// If we find nothing that matches, the kick out to command promt
	if (foundheader == 0) {
	  printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EM1004) ***\n\n");
	  fclose(importfile);
	  break;
	}

	// printf("SETHEAD: %d %d %d %d %d %d %d\n", regionsset[currentset], jpmodelset[currentset], kpmodelset[currentset], jptribmodelset[currentset], mininjectset[currentset][1], mininjectset[currentset][2], mininjectset[currentset][3]);

      }
      // If we can't find or access the file, kick back to the command prompt
      else {
	fprintf(stderr,"*** Error: Unable to locate or access the .brats file requested ***\n\n");
	break;
      }
      
      // Import the frequency data (all these functions can be found in fullimport.h)
      if(importfrequencies(importfilename, currentset, frequency, frequencyobs, fluxorder, imgnum[currentset]) != 0) {
	break; // Escape if there is a problem found in the function
      }


      if(importrawfluxdata(importfilename, currentset, imgnum[currentset], rms, fluxcalerror, minflux, maxflux, flux, bgbuff, xdim[currentset], ydim[currentset]) != 0) {
	break; // Escape if there is a problem found in the function
      }

      // Up the number of datasets once the raw data has been loaded
      numdatasets++;
      
      if (regionsset[currentset] == 1) {

	if (regmemallocated == 0 ) {
	  
	  regflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  regionarray = (int ***)calloc(MAXDATASETS, sizeof(int *));
	  regnumber = (int *)calloc(MAXDATASETS, sizeof(int));
	  fluxerror = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  regionlocx = (float **)calloc(MAXDATASETS, sizeof(float *));
	  regionlocy = (float **)calloc(MAXDATASETS, sizeof(float *));
	  setaveraged = (int *)calloc(MAXDATASETS, sizeof(int));
	  regionsize = (int **)calloc(MAXDATASETS, sizeof(int *));
	  regmaxflux = (float **)calloc(MAXDATASETS,sizeof(float *));
	  regminflux = (float **)calloc(MAXDATASETS,sizeof(float *));
	  
	  regmemallocated = 1;
	}

	importnumberofregions(importfilename, currentset, regnumber);
	
	if(importregiondata(importfilename, importcompression, currentset, imgnum[currentset], xdim[currentset], ydim[currentset], regionarray, regflux, regnumber, fluxerror, setaveraged, regionlocx, regionlocy, regionsize, regmaxflux, regminflux) != 0) {
	  regionsset[currentset] = 0;
	  break; // Escape if there is a problem found in the function
	}
      }
      

      if (jpmodelset[currentset] == 1) {

	// See what has been run in SETHEAD
	importfile = fopen(importfilename, "r");

	if (importfile != NULL) {
	
	  headerloop = 0;
	  line_num = 1;
	  foundheader = 0;
	  
	  while( fgets(line_storage, sizeof(line_storage), importfile) != NULL )  {

	    sscanf(line_storage,"%s", buffer);
	    // printf("line_storage: %s buffer: %s\n", line_storage, buffer);
	    	    
	    if (strcmp(buffer,"###ENDOFIMPORTFILE###") == 0)  {
	      printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EM1010) ***\n\n");
	      fclose(importfile);
	      jpmodelset[currentset] = 0;
	      break;
	    }

	    if (foundheader == 1) {
	      if ( (strcmp(buffer,"###END###") == 0) && (headerloop != expectedmodelheaders) ) {
		printf("\n*** ERROR: Finished reading JPFITHEADER before the expected number of iterations had been reached! The data may be corrupt ***\n\n");
		fclose(importfile);
		jpmodelset[currentset] = 0;
		break;
	      }
	      else if ( (strcmp(buffer,"###END###") == 0) && (headerloop == expectedmodelheaders) ) {
		printf("JPFITHEAD imported successfully\n");
		fclose(importfile);
		break;
	      }
	      else {
		if (headerloop == 0) {
		  sscanf(buffer,"%lf", &jpfield[currentset]);
		}
		else if (headerloop == 1) {
		  sscanf(buffer,"%lf", &jpinject[currentset]);
		}
		else if (headerloop == 2) {
		  sscanf(buffer,"%d", &setmodelres[currentset][1]);
		}

		headerloop++;
	      }	      
	    }
	    
	    if (strcmp(buffer,"###JPFITHEADER###") == 0)  {
	      printf("JPFITHEADER found on line %d\n", line_num);
	      foundheader = 1;
	    }
	    
	    line_num++;	    
	  }
	  
	  // If we find nothing that matches, the kick out to command promt
	  if (foundheader == 0) {
	    printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EM1011) ***\n\n");
	    fclose(importfile);
	    jpmodelset[currentset] = 0;
	    break;
	  }	  
	}
	// If we can't find or access the file, kick back to the command prompt
	else {
	  fprintf(stderr,"*** Error: Unable to locate or access the .brats file requested ***\n\n");
	  jpmodelset[currentset] = 0;
	  break;
	}
	
	
	// Setup the memory
	if (jpmodelmemoryset == 0) {

	  jpchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jpbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jpmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jpmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  jpbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jpageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jpageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  
	  jpmodelmemoryset = 1;
	}
	
	
	if (numdatasets > 1) {
	  
	  allocarray = (int *)calloc(MAXDATASETS, sizeof(int));
	  
	  realloc_jpchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	  realloc_jpbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	  realloc_jpmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	  realloc_jpmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  realloc_jpbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	  realloc_jpageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  realloc_jpageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  
	
	  for (o=0; o<numdatasets; o++) {
	    
	    if (jpfitmemset[o] == 1) {
	      
	      realloc_jpchisquared[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      realloc_jpbestage[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      realloc_jpmodelalpha[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      realloc_jpmodelflux[o] = (float **)calloc(regnumber[o]+1, sizeof(float *));
	      realloc_jpbestnorm[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      realloc_jpageerrorsplus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      realloc_jpageerrorsminus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      
	      allocarray[o] = 1;
	      
	    }
	  }
	  
	  for (o=0; o<numdatasets; o++) {
	    
	    if (jpfitmemset[o] == 1) {
	      
	      for (q=0; q<=regnumber[o]; q++) {
		
		//printf("q=%d\n", q);
		realloc_jpmodelflux[o][q] = (float *)calloc(setmodelres[o][1]+1, sizeof(float));
	      }
	    }
	  }
	  
	  for (o=0; o<numdatasets; o++) {
	    
	    for (p=0; p<=regnumber[o]; p++) {
	      
	      if (jpfitmemset[o] == 1) {
		
		realloc_jpchisquared[o][p] = jpchisquared[o][p];
		realloc_jpbestage[o][p] = jpbestage[o][p];
		realloc_jpmodelalpha[o][p] = jpmodelalpha[o][p];
		realloc_jpbestnorm[o][p] = jpbestnorm[o][p];
		realloc_jpageerrorsplus[o][p] = jpageerrorsplus[o][p];
		realloc_jpageerrorsminus[o][p] = jpageerrorsminus[o][p];
		
		for (q=0; q<=setmodelres[o][1]; q++) {
		  realloc_jpmodelflux[o][p][q] = jpmodelflux[o][p][q];
		}
		
	      }
	    }
	  }
	  
	  for (o=0; o<numdatasets-1; o++) { // -1 as we are yet to allocate the new dataset memory
	    
	    free(jpchisquared[o]);
	    free(jpbestage[o]);
	    free(jpmodelalpha[o]);
	    free(jpbestnorm[o]);
	    free(jpageerrorsplus[o]);
	    free(jpageerrorsminus[o]);
	    free(jpmodelflux[o]);	    
	  }
	  
	  free(jpchisquared);
	  free(jpbestage);
	  free(jpmodelalpha);
	  free(jpbestnorm);
	  free(jpageerrorsplus);
	  free(jpageerrorsminus);
	  free(jpmodelflux);
	  
	  jpfitmemset[o] = 0;
	  jpmodelmemoryset = 0;
	  
	  
	  // Make the first level of the array again
	  jpchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jpbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jpmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jpmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  jpbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jpageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jpageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  
	  jpmodelmemoryset = 1;
	}
	
	for (o=0; o<numdatasets; o++) {

	  jpchisquared[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	  jpbestage[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	  jpmodelalpha[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	  jpmodelflux[o] = (float **)calloc(regnumber[o]+1, sizeof(float *));
	  jpbestnorm[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	  jpageerrorsplus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	  jpageerrorsminus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));

	  
	  for (p=0; p<=regnumber[o]; p++) {
	    jpmodelflux[o][p] = (float *)calloc(setmodelres[o][1]+1, sizeof(float));
	  }
	  
	  jpfitmemset[o] = 1;
	}
	
	
	if (numdatasets > 1) {
	  for (o=0; o<numdatasets; o++) {
	    
	    if ( (o != currentset) && (allocarray[o] == 1) ) {
	      
	      for (p=0; p<=regnumber[o]; p++) {
		
		jpchisquared[o][p] = realloc_jpchisquared[o][p];
		jpbestage[o][p] = realloc_jpbestage[o][p];
		jpmodelalpha[o][p] = realloc_jpmodelalpha[o][p];
		jpbestnorm[o][p] = realloc_jpbestnorm[o][p];
		jpageerrorsplus[o][p] = realloc_jpageerrorsplus[o][p];
		jpageerrorsminus[o][p] = realloc_jpageerrorsminus[o][p];
		
		for (q=0; q<=setmodelres[o][1]; q++) {
		  
		  jpmodelflux[o][p][q] = realloc_jpmodelflux[o][p][q];
		}
	      }
	    }
	  }
	}
	
	if (numdatasets > 1) {
	  
	  // Free the reallocation arrays if they are set
	  
	  for (o=0; o<numdatasets-1; o++) {
	    
	    free(realloc_jpchisquared[o]);
	    free(realloc_jpbestage[o]);
	    free(realloc_jpmodelalpha[o]);
	    free(realloc_jpbestnorm[o]);
	    free(realloc_jpageerrorsplus[o]);
	    free(realloc_jpageerrorsminus[o]);
	    //free(realloc_jpmodelflux[o]);
	  }
	  
	  free(realloc_jpchisquared);
	  free(realloc_jpbestage);
	  free(realloc_jpmodelalpha);
	  free(realloc_jpbestnorm);
	  free(realloc_jpageerrorsplus);
	  free(realloc_jpageerrorsminus);
	  free(realloc_jpmodelflux);
	  
	  free(allocarray);
	}
	

	if(importjpmodel(importfilename, currentset, imgnum[currentset], xdim[currentset], ydim[currentset], regnumber[currentset], jpchisquared, jpbestage, jpmodelalpha, jpmodelflux, jpbestnorm, jpageerrorsplus, jpageerrorsminus, setmodelres[currentset][1]) != 0) {
	  jpmodelset[currentset] = 0;
	  break; // Escape if there is a problem found in the function
	}
      }
      

      if (kpmodelset[currentset] == 1) {

	// See what has been run in SETHEAD
	importfile = fopen(importfilename, "r");

	if (importfile != NULL) {
	
	  headerloop = 0;
	  line_num = 1;
	  foundheader = 0;
	  
	  while( fgets(line_storage, sizeof(line_storage), importfile) != NULL )  {

	    sscanf(line_storage,"%s", buffer);
	    // printf("line_storage: %s buffer: %s\n", line_storage, buffer);
	    
	    
	    if (strcmp(buffer,"###ENDOFIMPORTFILE###") == 0)  {
	      printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EM1010) ***\n\n");
	      fclose(importfile);
	      kpmodelset[currentset] = 0;
	      break;
	    }

	    if (foundheader == 1) {
	      if ( (strcmp(buffer,"###END###") == 0) && (headerloop != expectedmodelheaders) ) {
		printf("\n*** ERROR: Finished reading KPFITHEADER before the expected number of iterations had been reached! The data may be corrupt ***\n\n");
		fclose(importfile);
		kpmodelset[currentset] = 0;
		break;
	      }
	      else if ( (strcmp(buffer,"###END###") == 0) && (headerloop == expectedmodelheaders) ) {
		printf("KPFITHEAD imported successfully\n");
		fclose(importfile);
		break;
	      }
	      else {
		if (headerloop == 0) {
		  sscanf(buffer,"%lf", &kpfield[currentset]);
		}
		else if (headerloop == 1) {
		  sscanf(buffer,"%lf", &kpinject[currentset]);
		}
		else if (headerloop == 2) {
		  sscanf(buffer,"%d", &setmodelres[currentset][2]);
		}

		headerloop++;
	      }
	    }
	    
	    if (strcmp(buffer,"###KPFITHEADER###") == 0)  {
	      printf("KPFITHEADER found on line %d\n", line_num);
	      foundheader = 1;
	    }
	    
	    line_num++;
	  }
	  
	  // If we find nothing that matches, the kick out to command promt
	  if (foundheader == 0) {
	    printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EM1011) ***\n\n");
	    fclose(importfile);
	    kpmodelset[currentset] = 0;
	    break;
	  }
	  
	}
	// If we can't find or access the file, kick back to the command prompt
	else {
	  fprintf(stderr,"*** Error: Unable to locate or access the .brats file requested ***\n\n");
	  kpmodelset[currentset] = 0;
	  break;
	}
	
	
	// Setup the memory
	if (kpmodelmemoryset == 0) {

	  kpchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	  kpbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	  kpmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	  kpmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  kpbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	  kpageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  kpageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  
	  kpmodelmemoryset = 1;
	}
	
	
	if (numdatasets > 1) {
	  
	  allocarray = (int *)calloc(MAXDATASETS, sizeof(int));
	  
	  realloc_kpchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	  realloc_kpbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	  realloc_kpmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	  realloc_kpmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  realloc_kpbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	  realloc_kpageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  realloc_kpageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  
	  
	  for (o=0; o<numdatasets; o++) {
	    
	    if (kpfitmemset[o] == 1) {
	      
	      realloc_kpchisquared[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      realloc_kpbestage[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      realloc_kpmodelalpha[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      realloc_kpmodelflux[o] = (float **)calloc(regnumber[o]+1, sizeof(float *));
	      realloc_kpbestnorm[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      realloc_kpageerrorsplus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      realloc_kpageerrorsminus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      
	      allocarray[o] = 1;
	    }
	  }
	  
	  for (o=0; o<numdatasets; o++) {
	    
	    if (kpfitmemset[o] == 1) {
	      
	      for (q=0; q<=regnumber[o]; q++) {
		
		//printf("q=%d\n", q);
		realloc_kpmodelflux[o][q] = (float *)calloc(setmodelres[o][2]+1, sizeof(float));
	      }
	    }
	  }
	  
	  for (o=0; o<numdatasets; o++) {
	    
	    for (p=0; p<=regnumber[o]; p++) {
	      
	      if (kpfitmemset[o] == 1) {
		
		realloc_kpchisquared[o][p] = kpchisquared[o][p];
		realloc_kpbestage[o][p] = kpbestage[o][p];
		realloc_kpmodelalpha[o][p] = kpmodelalpha[o][p];
		realloc_kpbestnorm[o][p] = kpbestnorm[o][p];
		realloc_kpageerrorsplus[o][p] = kpageerrorsplus[o][p];
		realloc_kpageerrorsminus[o][p] = kpageerrorsminus[o][p];
		
		for (q=0; q<=setmodelres[o][2]; q++) {
		  realloc_kpmodelflux[o][p][q] = kpmodelflux[o][p][q];
		}
		
	      }
	    }
	    
	  }
	  
	  for (o=0; o<numdatasets-1; o++) { // -1 as we are yet to allocate the new dataset memory
	    
	    free(kpchisquared[o]);
	    free(kpbestage[o]);
	    free(kpmodelalpha[o]);
	    free(kpbestnorm[o]);
	    free(kpageerrorsplus[o]);
	    free(kpageerrorsminus[o]);
	    free(kpmodelflux[o]);
	  }
	  
	  free(kpchisquared);
	  free(kpbestage);
	  free(kpmodelalpha);
	  free(kpbestnorm);
	  free(kpageerrorsplus);
	  free(kpageerrorsminus);
	  free(kpmodelflux);
	  
	  kpfitmemset[o] = 0;
	  kpmodelmemoryset = 0;
	  
	  
	  // Make the first level of the array again
	  kpchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	  kpbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	  kpmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	  kpmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  kpbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	  kpageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  kpageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  
	  kpmodelmemoryset = 1;
	    
	}
	
	for (o=0; o<numdatasets; o++) {

	  kpchisquared[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	  kpbestage[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	  kpmodelalpha[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	  kpmodelflux[o] = (float **)calloc(regnumber[o]+1, sizeof(float *));
	  kpbestnorm[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	  kpageerrorsplus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	  kpageerrorsminus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));

	  
	  for (p=0; p<=regnumber[o]; p++) {
	    kpmodelflux[o][p] = (float *)calloc(setmodelres[o][2]+1, sizeof(float));
	  }
	  
	  kpfitmemset[o] = 1;
	  
	}
	
	
	if (numdatasets > 1) {
	  for (o=0; o<numdatasets; o++) {
	    
	    if ( (o != currentset) && (allocarray[o] == 1) ) {
	      
	      for (p=0; p<=regnumber[o]; p++) {
		
		kpchisquared[o][p] = realloc_kpchisquared[o][p];
		kpbestage[o][p] = realloc_kpbestage[o][p];
		kpmodelalpha[o][p] = realloc_kpmodelalpha[o][p];
		kpbestnorm[o][p] = realloc_kpbestnorm[o][p];
		kpageerrorsplus[o][p] = realloc_kpageerrorsplus[o][p];
		kpageerrorsminus[o][p] = realloc_kpageerrorsminus[o][p];
		
		for (q=0; q<=setmodelres[o][2]; q++) {
		  
		  kpmodelflux[o][p][q] = realloc_kpmodelflux[o][p][q];
		}
	      }
	    }
	  }
	}
	
	if (numdatasets > 1) {
	  
	  // Free the reallocation arrays if they are set
	  
	  for (o=0; o<numdatasets-1; o++) {
	    
	    free(realloc_kpchisquared[o]);
	    free(realloc_kpbestage[o]);
	    free(realloc_kpmodelalpha[o]);
	    free(realloc_kpbestnorm[o]);
	    free(realloc_kpageerrorsplus[o]);
	    free(realloc_kpageerrorsminus[o]);
	    //free(realloc_kpmodelflux[o]);
	    
	  }
	  
	  free(realloc_kpchisquared);
	  free(realloc_kpbestage);
	  free(realloc_kpmodelalpha);
	  free(realloc_kpbestnorm);
	  free(realloc_kpageerrorsplus);
	  free(realloc_kpageerrorsminus);
	  free(realloc_kpmodelflux);
	  
	  free(allocarray);
	  
	}
	

	if(importkpmodel(importfilename, currentset, imgnum[currentset], xdim[currentset], ydim[currentset], regnumber[currentset], kpchisquared, kpbestage, kpmodelalpha, kpmodelflux, kpbestnorm, kpageerrorsplus, kpageerrorsminus, setmodelres[currentset][2]) != 0) {
	  kpmodelset[currentset] = 0;
	  break; // Escape if there is a problem found in the function
	}

      }
      

      if (jptribmodelset[currentset] == 1) {

	// See what has been run in SETHEAD
	importfile = fopen(importfilename, "r");

	if (importfile != NULL) {
	
	  headerloop = 0;
	  line_num = 1;
	  foundheader = 0;
	  
	  while( fgets(line_storage, sizeof(line_storage), importfile) != NULL )  {

	    sscanf(line_storage,"%s", buffer);
	    // printf("line_storage: %s buffer: %s\n", line_storage, buffer);
	    
	    
	    if (strcmp(buffer,"###ENDOFIMPORTFILE###") == 0)  {
	      printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EM1010) ***\n\n");
	      fclose(importfile);
	      jptribmodelset[currentset] = 0;
	      break;
	    }

	    if (foundheader == 1) {
	      if ( (strcmp(buffer,"###END###") == 0) && (headerloop != expectedmodelheaders) ) {
		printf("\n*** ERROR: Finished reading JPTRIBFITHEADER before the expected number of iterations had been reached! The data may be corrupt ***\n\n");
		fclose(importfile);
		jptribmodelset[currentset] = 0;
		break;
	      }
	      else if ( (strcmp(buffer,"###END###") == 0) && (headerloop == expectedmodelheaders) ) {
		printf("JPTRIBFITHEAD imported successfully\n");
		fclose(importfile);
		break;
	      }
	      else {
		if (headerloop == 0) {
		  sscanf(buffer,"%lf", &jptribfield[currentset]);
		}
		else if (headerloop == 1) {
		  sscanf(buffer,"%lf", &jptribinject[currentset]);
		}
		else if (headerloop == 2) {
		  sscanf(buffer,"%d", &setmodelres[currentset][3]);
		}

		headerloop++;
	      }
	      
	    }
	    
	    if (strcmp(buffer,"###JPTRIBFITHEADER###") == 0)  {
	      printf("JPTRIBFITHEADER found on line %d\n", line_num);
	      foundheader = 1;
	    }
	    
	    line_num++;
	    
	  }
	  
	  // If we find nothing that matches, the kick out to command promt
	  if (foundheader == 0) {
	    printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EM1011) ***\n\n");
	    fclose(importfile);
	    jptribmodelset[currentset] = 0;
	    break;
	  }
	  
	}
	// If we can't find or access the file, kick back to the command prompt
	else {
	  fprintf(stderr,"*** Error: Unable to locate or access the .brats file requested ***\n\n");
	  jptribmodelset[currentset] = 0;
	  break;
	}
	
	
	// Setup the memory
	if (jptribmodelmemoryset == 0) {

	  jptribchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jptribbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jptribmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jptribmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  jptribbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jptribbleageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jptribbleageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  
	  jptribmodelmemoryset = 1;
	}
	
	
	if (numdatasets > 1) {
	  
	  allocarray = (int *)calloc(MAXDATASETS, sizeof(int));
	  
	  realloc_jptribchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	  realloc_jptribbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	  realloc_jptribmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	  realloc_jptribmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  realloc_jptribbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	  realloc_jptribbleageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  realloc_jptribbleageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  
	  
	  for (o=0; o<numdatasets; o++) {
	    
	    if (jptribfitmemset[o] == 1) {
	      
	      realloc_jptribchisquared[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      realloc_jptribbestage[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      realloc_jptribmodelalpha[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      realloc_jptribmodelflux[o] = (float **)calloc(regnumber[o]+1, sizeof(float *));
	      realloc_jptribbestnorm[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      realloc_jptribbleageerrorsplus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      realloc_jptribbleageerrorsminus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	      
	      allocarray[o] = 1;
	      
	    }
	  }
	  
	  for (o=0; o<numdatasets; o++) {
	    
	    if (jptribfitmemset[o] == 1) {
	      
	      for (q=0; q<=regnumber[o]; q++) {
		
		//printf("q=%d\n", q);
		realloc_jptribmodelflux[o][q] = (float *)calloc(setmodelres[o][3]+1, sizeof(float));
	      }
	    }
	  }
	  
	  for (o=0; o<numdatasets; o++) {
	    
	    for (p=0; p<=regnumber[o]; p++) {
	      
	      if (jptribfitmemset[o] == 1) {
		
		realloc_jptribchisquared[o][p] = jptribchisquared[o][p];
		realloc_jptribbestage[o][p] = jptribbestage[o][p];
		realloc_jptribmodelalpha[o][p] = jptribmodelalpha[o][p];
		realloc_jptribbestnorm[o][p] = jptribbestnorm[o][p];
		realloc_jptribbleageerrorsplus[o][p] = jptribbleageerrorsplus[o][p];
		realloc_jptribbleageerrorsminus[o][p] = jptribbleageerrorsminus[o][p];
		
		for (q=0; q<=setmodelres[o][3]; q++) {
		  realloc_jptribmodelflux[o][p][q] = jptribmodelflux[o][p][q];
		}
		
	      }
	    }
	    
	  }
	  
	  for (o=0; o<numdatasets-1; o++) { // -1 as we are yet to allocate the new dataset memory
	    
	    free(jptribchisquared[o]);
	    free(jptribbestage[o]);
	    free(jptribmodelalpha[o]);
	    free(jptribbestnorm[o]);
	    free(jptribbleageerrorsplus[o]);
	    free(jptribbleageerrorsminus[o]);
	    free(jptribmodelflux[o]);
	    
	  }
	  
	  free(jptribchisquared);
	  free(jptribbestage);
	  free(jptribmodelalpha);
	  free(jptribbestnorm);
	  free(jptribbleageerrorsplus);
	  free(jptribbleageerrorsminus);
	  free(jptribmodelflux);
	  
	  jptribfitmemset[o] = 0;
	  jptribmodelmemoryset = 0;
	  
	  
	  // Make the first level of the array again
	  jptribchisquared = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jptribbestage = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jptribmodelalpha = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jptribmodelflux = (float ***)calloc(MAXDATASETS, sizeof(float **));
	  jptribbestnorm = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jptribbleageerrorsplus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  jptribbleageerrorsminus = (float **)calloc(MAXDATASETS, sizeof(float *));
	  
	  jptribmodelmemoryset = 1;
	    
	}
	
	for (o=0; o<numdatasets; o++) {

	  jptribchisquared[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	  jptribbestage[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	  jptribmodelalpha[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	  jptribmodelflux[o] = (float **)calloc(regnumber[o]+1, sizeof(float *));
	  jptribbestnorm[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	  jptribbleageerrorsplus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));
	  jptribbleageerrorsminus[o] = (float *)calloc(regnumber[o]+1, sizeof(float));

	  
	  for (p=0; p<=regnumber[o]; p++) {
	    jptribmodelflux[o][p] = (float *)calloc(setmodelres[o][3]+1, sizeof(float));
	  }
	  
	  jptribfitmemset[o] = 1;
	  
	}
	
	
	if (numdatasets > 1) {
	  for (o=0; o<numdatasets; o++) {
	    
	    if ( (o != currentset) && (allocarray[o] == 1) ) {
	      
	      for (p=0; p<=regnumber[o]; p++) {
		
		jptribchisquared[o][p] = realloc_jptribchisquared[o][p];
		jptribbestage[o][p] = realloc_jptribbestage[o][p];
		jptribmodelalpha[o][p] = realloc_jptribmodelalpha[o][p];
		jptribbestnorm[o][p] = realloc_jptribbestnorm[o][p];
		jptribbleageerrorsplus[o][p] = realloc_jptribbleageerrorsplus[o][p];
		jptribbleageerrorsminus[o][p] = realloc_jptribbleageerrorsminus[o][p];
		
		for (q=0; q<=setmodelres[o][3]; q++) {
		  
		  jptribmodelflux[o][p][q] = realloc_jptribmodelflux[o][p][q];
		}
	      }
	    }
	  }
	}
	
	if (numdatasets > 1) {
	  
	  // Free the reallocation arrays if they are set
	  
	  for (o=0; o<numdatasets-1; o++) {
	    
	    free(realloc_jptribchisquared[o]);
	    free(realloc_jptribbestage[o]);
	    free(realloc_jptribmodelalpha[o]);
	    free(realloc_jptribbestnorm[o]);
	    free(realloc_jptribbleageerrorsplus[o]);
	    free(realloc_jptribbleageerrorsminus[o]);
	    //free(realloc_jptribmodelflux[o]);
	    
	  }
	  
	  free(realloc_jptribchisquared);
	  free(realloc_jptribbestage);
	  free(realloc_jptribmodelalpha);
	  free(realloc_jptribbestnorm);
	  free(realloc_jptribbleageerrorsplus);
	  free(realloc_jptribbleageerrorsminus);
	  free(realloc_jptribmodelflux);
	  
	  free(allocarray);
	  
	}
	

	if(importjptribmodel(importfilename, currentset, imgnum[currentset], xdim[currentset], ydim[currentset], regnumber[currentset], jptribchisquared, jptribbestage, jptribmodelalpha, jptribmodelflux, jptribbestnorm, jptribbleageerrorsplus, jptribbleageerrorsminus, setmodelres[currentset][3]) != 0) {
	  jptribmodelset[currentset] = 0;
	  break; // Escape if there is a problem found in the function
	}

      }

      for (q=1; q<NUMBEROFMODELS; q++) {

	model = q;

	if (mininjectset[currentset][model] > 0) {

	  // Setup the memory
	  inj_sumchisquared[currentset][model] = (float *)realloc(inj_sumchisquared[currentset][model], (mininjectset[currentset][model]+1) * sizeof(float));
	  inj_regchisquared[currentset][model] = (float *)realloc(inj_regchisquared[currentset][model], (regnumber[currentset]+1) * sizeof(float));
	  inj_regbestinject[currentset][model] = (float *)realloc(inj_regbestinject[currentset][model], (regnumber[currentset]+1) * sizeof(float));

	  // Set the default values
	  for (i=0; i<=mininjectset[currentset][model]; i++) {
	    inj_sumchisquared[currentset][model][i] = 0.0;
	  }

	  for (i=0; i<=regnumber[currentset]; i++) {
	    inj_regchisquared[currentset][model][i] = 1e23;
	    inj_regbestinject[currentset][model][i] = 1e23;
	  }


	  if(importfindinject(importfilename, inj_sumchisquared[currentset][model], inj_regbestinject[currentset][model], inj_regchisquared[currentset][model], mininjectset[currentset][model], &mininjectstore[currentset][model], &maxinjectstore[currentset][model], model, regnumber[currentset], exportversion) != 0) {
	    mininjectset[currentset][model] = 0;
	    break; // Escape if there is a problem found in the function
	  }
	}
      }

      break;

      
    case 113: // Set what type of fitting should be done for the spectral index

      validentry = 0;
      tmp_indexcalctype = indexcalctype;

      printf("Scale type currently set to %d (DEFAULT: %d)\n", indexcalctype, DEFAULTINDEXCALCTYPE);

      while (validentry == 0) {

	cmdbuffer = readline("Select a scale type to use when mapping. 1 = Least squares (unweighted, no errors, inbuilt), 2 = Weighted least squares (with errors, GSL): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	indexcalctype = strtol(cmdbuffer, &endptr, 10);

	if ( (indexcalctype == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (indexcalctype < 0) {
	    printf("Escaping command...\n");
	    indexcalctype = tmp_indexcalctype;
	    break;
	  }
	  else if (indexcalctype > 2) {
	    printf("\nInvalid entry, please try again...\n\n");
	    continue;
	  }
	  else {
	    printf("\nSpectral index fitting type is now set to %d\n\n", indexcalctype);
	    validentry = 1;
	  }
	}
      }

      break;


    case 114: // Plot the spectral indices for any GSL fit

      if (numdatasets > 0) {

	// Open the directory that is currently set and check it exists
	dir = opendir(imageloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	  break;
	}
	
	// List the available datasets
	
	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");
	
	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {
	  
	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);
	  
	}
	
	printf("\n========================================================================\n");
    
	
	// Get the reference data set
	validentry = 0;
	escape = 0;
	
	while (validentry == 0) {
	  
	  cmdbuffer = readline("Enter the data set for which to plot the spectral index (-1 to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if (currentset >= numdatasets) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
	  else if (regionsset[currentset] != 1) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    continue;
	  }
	  else if (specindexset[currentset] != 1) { // Check the spectral index has been set
	    fprintf(stderr,"\nError: Spectral indices have not yet been fitted to this data set! Please run the specindex command first.\n\n");
	    continue;
	  }
	  else if (specindextype[currentset] != 2) { // Check the spectral index has been set
	    fprintf(stderr,"\nError: Invalid fitting method used. Weighted (GSL) least squares fitting must be used to enable this function.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}
	
	
	tmp_regflux = (float *)calloc(imgnum[currentset], sizeof(float));
	tmp_fluxerror = (float *)calloc(imgnum[currentset], sizeof(float));
	tmp_fluxfreq = (float *)calloc(imgnum[currentset], sizeof(float));
	tmp_fluxerrorplus = (float *)calloc(imgnum[currentset], sizeof(float));
	tmp_fluxerrorminus = (float *)calloc(imgnum[currentset], sizeof(float));

	tmp_minflux = 1e32;
	tmp_maxflux = -1e32;
      
	// Check that they really want to export this many maps if it is over 10
	if ( ((float)regnumber[currentset]/(float)skip)>10 )  {

	  if (export == 1) {
	    printf("\nWarning: You are about to export %d plots! This cannot be escaped early. ", (int)(regnumber[currentset]/(float)skip) );
	  }
	  else {
	    printf("\nWarning: You are about to view %d plots! This cannot be escaped early. ", (int)(regnumber[currentset]/(float)skip) );
	  }

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Are you sure you want to continue? (1 Yes, 0 No): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    warningcheck = strtol(cmdbuffer, &endptr, 10);

	    if ( (warningcheck == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }
		
	    if ( (warningcheck < 0) || (warningcheck > 1) ) { // Check the value is either 1 or 0
	      printf("\nValue must be either Yes (1) or no (0). Please try again...\n\n");
	      continue;
	    }
	    else if (warningcheck == 0) { // If 0, exit
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else if (warningcheck == 1) { // If 1, continue
	      printf("Proceeding...\n");
	      validentry = 1;
	    }
	    else { // Catch anthing odd
	      printf("\nThe value of warningcheck is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", warningcheck);
	      break;
	    }
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}
      
	// Open a window if we are not exporting
	if (export != 1) {
	  sprintf(output,"/xs");
	  cpgopen(output);
	}

	firstplot = 1;

	for (i=skip; i<=regnumber[currentset]; i+=skip) {
	  // Double plot the first cycle to stop it skipping the carriage return
	  if ( (i != skip) && (firstplot == 1) && (export != 1) ) {
	    i = skip;
	    firstplot = 0;
	  }

	  tmp_minflux = 1e32; // Reset the min max after every plot
	  tmp_maxflux = -1e32;

	  for (j=0; j<imgnum[currentset]; j++) {
	    
	    tmp_regflux[j] = log10(regflux[currentset][fluxorder[currentset][j]][i]);
	    tmp_fluxfreq[j] = log10(frequencyobs[currentset][fluxorder[currentset][j]]);
	    
	    tmp_fluxerrorplus[j] = log10(regflux[currentset][fluxorder[currentset][j]][i] + fluxerror[currentset][fluxorder[currentset][j]][i]);
	    
	    tmp_fluxerrorminus[j] = log10(regflux[currentset][fluxorder[currentset][j]][i] - fluxerror[currentset][fluxorder[currentset][j]][i]);
	

	    if (pow(10, tmp_fluxerrorminus[j]) < tmp_minflux) {
	      tmp_minflux = pow(10, tmp_fluxerrorminus[j]);
	    }

	    if (pow(10, tmp_fluxerrorplus[j]) > tmp_maxflux) {
	      tmp_maxflux =  pow(10, tmp_fluxerrorplus[j]);
	    }

	    for (a=0; a<=specindex_modelres[currentset]; a++) {
	      if (specindex_modelflux[currentset][i][a] < tmp_minflux) {
		tmp_minflux = pow(10, specindex_modelflux[currentset][i][a]);
	      }
	      if (specindex_modelflux[currentset][i][a] > tmp_maxflux) {
		tmp_maxflux = pow(10, specindex_modelflux[currentset][i][a]);
	      }

	      if (pow(10,(specindex_modelflux[currentset][i][a] - specindex_modelerror[currentset][i][a])) < tmp_minflux) {
		tmp_minflux = pow(10, specindex_modelflux[currentset][i][a] - specindex_modelerror[currentset][i][a]);
	      }
	      if (pow(10,(specindex_modelflux[currentset][i][a] + specindex_modelerror[currentset][i][a])) > tmp_maxflux) {
		tmp_maxflux = pow(10, specindex_modelflux[currentset][i][a] + specindex_modelerror[currentset][i][a]);
	      }
	    }
	  }
	  
	  // If exporting, output open a new session for each with a different name
	  if (export == 1) {
	    
	    time( &currenttime );
	    time_struct = localtime ( &currenttime );
	    
	    strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);
	    
	    sprintf(output,"%s/%s_Reg%d_ModObs_SpecIndex_%s.%s/%s", imageloc, settarget[currentset], i, timebuff, imagetype, imagetype);
	    
	    printf("\nExporting as %s\n", output);
	    
	    cpgopen(output);
	  }
	  else {
	    // Export output is aloready set, just give the region output
	    if (firstplot != 1) {
	      printf("Current plot is of region %d\n", i);
	    }
	  }
		
	  plotfluxvsspecindex(tmp_regflux, tmp_fluxerrorplus, tmp_fluxerrorminus, symbol, imgnum[currentset], tmp_fluxfreq, autoscalex, autoscaley, userminx, usermaxx, userminy, usermaxy, frequencyobs[currentset][fluxorder[currentset][0]], frequencyobs[currentset][fluxorder[currentset][imgnum[currentset]-1]], tmp_minflux, tmp_maxflux, settarget[currentset], titles, labels, axismarks, paramlabels, specindex_modelflux[currentset][i], specindex_modelerror[currentset][i], specindex_modelres[currentset], fluxorder[currentset], alpha[currentset][i], largetxt, specindexerror[currentset][i], specindex_errortype[currentset]);
	
	
	  if ( export == 1 ) {
	    cpgclos();
	  }

	}
	 
	// Free the temporary variables
	free(tmp_fluxfreq);
	free(tmp_regflux);
	free(tmp_fluxerror);
	free(tmp_fluxerrorplus);
	free(tmp_fluxerrorminus);

	if (export != 1) {
	  cpgclos();
	}
	
      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;

      
    case 115: // Make a chi-squared map for the spectral indices of a given dataset


      if (numdatasets > 0) {

	if (export == 1) {

	  // Open the directory that is currently set and check it exists
	  dir = opendir(imageloc);

	  if (dir != NULL) {
	    closedir (dir); // If we find it move on without comment
	  }
	  else {
	    printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	    break;
	  }
	}

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to map the chi-squared values (-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (specindexset[currentset] != 1) ) { // Check the spectral indices have been set
	    fprintf(stderr,"\nError: Spectral indices have not yet been fitted to this data set! Please run the specindex command first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (imgnum[currentset] <= 2) ) { // Check the fitting type used
	    fprintf(stderr,"\nError: Chi-squared values are not available for spectral index fitting using only two data points.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (specindextype[currentset] == 1) ) { // Check its a GSL fitting type
	    fprintf(stderr,"\nError: Chi-squared map is only valid for spectral index fitting which uses the GSL functions! Spectral index calculation type can be set using the specindexcalctype command.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if ( (contours == 1) && (currentset != 888) ) {

	  if ( (exportfits != 1)  || (export != 1) ) {

	    validentry = 0;

	    printf("============================================\n\n");
	    printf(" Index  |  Frequency (Obs.)  |  Frequency (Rest)  \n\n");

	    // Loop through each of the sets
	    for (a=0; a<imgnum[currentset]; a++) {

	      printf("   %d       %.2e Hz             %.2e Hz\n", fluxorder[currentset][a], frequencyobs[currentset][fluxorder[currentset][a]], frequency[currentset][fluxorder[currentset][a]]);

	    }

	    printf("\n============================================\n");
    

	    while (validentry == 0) {

	      cmdbuffer = readline("Select a maps frequency for the contours (-1 to escape): ");
	      if (cmdbuffer[0] != 0) {
		add_history(cmdbuffer);
	      }

	      chomp(cmdbuffer);

	      currentmap = strtol(cmdbuffer, &endptr, 10);

	      if ( (currentmap == 0) && (cmdbuffer == endptr) ) {
		printf("\nInvalid input, please try again...\n\n");
		continue;
	      }

	      if (currentmap < 0) {
		printf("Escaping command...\n");
		escape = 1;
		break;
	      }
	      else if (currentmap >= imgnum[currentset]) { // Check the set exists
		printf("\nInvalid selection, please choose again...\n\n");
		continue;
	      }
	      else {
		validentry = 1;
	      }
	    }
	  }
	}
	else if ( (contours == 1) && (currentset == 888) ) { // Multimaps always uses map 0 as reference. This can easily be changed if needed
	  currentmap = 0;
	}
	else {
	  currentmap = 0;
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if ( (currentset == 888) && (export == 0) ) {
	  looplimit = numdatasets;
	  extwincontrol = 1;
	  cpgopen("/xs");
	}
	else if ( (currentset == 888) && (export == 1) ) {
	  looplimit = numdatasets;
	  extwincontrol = 0;
	}
	else {
	  looplimit = 1;
	  extwincontrol = 0;
	}
      

	for (j=0; j<looplimit; j++){

	  if (looplimit > 1) {
	    currentset = j;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = j;
	  }

	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets
	  if (regionsset[currentset] != 1) {
	    printf("\nThe regions for data set %d have not been set! Skipping this set...\n", j);
	    continue;	  
	  }

	  if (specindexset[currentset] != 1) {
	    printf("\nThe spectral index values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;
	  }

	  if (specindextype[currentset] == 1) { // Check its a GSL fitting type
	    printf("\nThe spectral index values for data set %d do not use the correct fitting type! Skipping this set...\n\n", j);
	    continue;
	  }

	  if (imgnum[currentset] <= 2) { // Check the fitting type used
	    printf("\nChi-squared values are not available for spectral index fitting using only two data points. Skipping this set...\n\n");
	    continue;
	  }

	  dof = (imgnum[currentset] - 2);

	  if (suppresscdf !=1) {
	    // print out the confidence levels
	    interval = 0.68;
	    siglevel68 = gsl_cdf_chisq_Pinv(interval, dof);
	    printf("\nModel can be rejected with 68 per cent confidence at X^2 > %.2f (Reduced: %.2f)\n", siglevel68, siglevel68 / dof);
	    interval = 0.90;
	    siglevel90 = gsl_cdf_chisq_Pinv(interval, dof);
	    printf("Model can be rejected with 90 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel90, siglevel90 / dof);
	    interval = 0.95;
	    siglevel95 = gsl_cdf_chisq_Pinv(interval, dof);
	    printf("Model can be rejected with 95 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel95, siglevel95 / dof);
	    interval = 0.99;
	    siglevel99 = gsl_cdf_chisq_Pinv(interval, dof);
	    printf("Model can be rejected with 99 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel99, siglevel99 / dof);
	    interval = 0.995;
	    siglevel995 = gsl_cdf_chisq_Pinv(interval, dof);
	    printf("Model can be rejected with 99.5 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel995, siglevel995 / dof);
	    interval = 0.999;
	    siglevel999 = gsl_cdf_chisq_Pinv(interval, dof);
	    printf("Model can be rejected with 99.9 per cent confidence X^2 > %.2f (Reduced: %.2f)\n", siglevel999, siglevel999 / dof);
	    interval = 0.9999;
	    siglevel9999 = gsl_cdf_chisq_Pinv(interval, dof);
	    printf("Model can be rejected with 99.99 per cent confidence X^2 > %.2f (Reduced: %.2f)\n\n", siglevel9999, siglevel9999 / dof);
	  }

	  time( &currenttime );
	  time_struct = localtime ( &currenttime );

	  strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	  if ( ( export == 1 ) && (exportfits != 1) ) {

	    sprintf(output,"%s/%s_ChiSquared_SpecIndex_%s.%s/%s", imageloc, settarget[currentset], timebuff, imagetype, imagetype);
	    printf("\nExporting as %s\n", output);
	  }
	  else if ( ( export == 1 ) && (exportfits == 1) ) {
	    sprintf(output,"%s/%s_ChiSquared_SpecIndex_%s.fits", imageloc, settarget[currentset], timebuff);
	    printf("\nExporting as %s\n", output);
	  }
	  else {
	    sprintf(output,"/xs");
	  }

	  if (usereduced == 1) {
	    sprintf(maptitle,"Map of Reduced Spectral Index Chi-Squared Values for %s", settarget[currentset]);
	  }
	  else {
	    sprintf(maptitle,"Map of Spectral Index Chi-Squared Values for %s", settarget[currentset]);
	  }
	
	  model = 999;

	  chisquaredmap(specindexchisquared[currentset], regionarray[currentset], xdim[currentset], ydim[currentset], settarget[currentset], imgnum[currentset], border[currentset], zoom, maptitle, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, flux[currentset][currentmap], firstcontour, contourlogs, output, jpinject[currentset], jpfield[currentset], paramlabels, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, usereduced, xshift, yshift, export, exportfits, largetxt);

	  /*if (extwincontrol == 1) {
	    cpgclos();
	    extwincontrol = 0;
	  }*/
	}

	if (extwincontrol == 1) {
	  cpgclos();
	  extwincontrol = 0;
	}
      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;


    case 116: // Resize all images in a specific folder

      resizeimages(imageloc);

      break;


    case 117: // Plot the CI (off) model for an arbitary normalisation the model (off vs standard CI) is set at the initial call

      if (export == 1) {

	// Open the directory that is currently set and check it exists
	dir = opendir(imageloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	  break;
	}

	time( &currenttime );
	time_struct = localtime ( &currenttime );

	strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);
	if (model == 5) {
	  sprintf(output,"%s/ExCIOFFModel_%s.%s/%s", imageloc,  timebuff, imagetype, imagetype);
	}
	else if (model == 6) {
	  sprintf(output,"%s/ExCIModel_%s.%s/%s", imageloc,  timebuff, imagetype, imagetype);
	}
	else{
	  sprintf(output,"%s/ExUnknownCIModel_%s.%s/%s", imageloc,  timebuff, imagetype, imagetype);
	}

	printf("\nExporting as %s\n", output);

	cpgopen(output);
      }
      else {
	sprintf(output,"/xs");
      }

      if (model == 5) {

	if (minmodelmyears < 1) {
	  printf("\nThe on time cannot be exactly 0 for a CI off model as the source must be on for at least some period to be visible. Setting the minimum on time to 0.01 Myr.\n");
	  tmp_minmodelmyears = 0; // Changed to 0.01 in module
	}
	else {
	  tmp_minmodelmyears = minmodelmyears;
	}

	if (varyoffage != 1) {

	  if (minmodelmyears > modelmyears) {
	    printf("\nError: The minimum on time (%d Myrs) must be less than the maximum on time (%d Myrs). Please use the minmodelmyears and/or the modelmyears commands to set a new range.\n\n", minmodelmyears, modelmyears);
	    break;
	  }
	  tmp_modelminoff = modelmaxoff; // Temporarily set min and max to be the same. This means we never loop the off age.
	  tmp_modelmyears = modelmyears; // We need to set the tmp_ version to something (check already done above).
	}
	else {
	  tmp_modelminoff = modelminoff;
	
	  if (modelmyears < 1) {
	    printf("\nThe on time cannot be exactly 0 for a CI off model as the source must be on for at least some period to be visible. Setting the maximum on time to 0.1 Myr.\n");
	    tmp_modelmyears = 0;
	  }
	  else {
	    tmp_modelmyears = modelmyears;
	  }

	  if (modelminoff > modelmaxoff) {
	    printf("\nError: The minimum off time (%d Myrs) must be less than the maximum off time (%d Myrs). Please use the modelminoff and/or the modelmaxoff commands to set a new range.\n\n", modelminoff, modelmaxoff);
	    break;
	  }
	}

	if ((modelmaxoff-tmp_modelminoff) == 0) {
	  printf("\nPlotting an example CI off model between %.2e Hz and %.2e Hz with an on time of between %d and %d Myrs and an off time of %d Myrs\n\n", minmodelfreq, maxmodelfreq, minmodelmyears, modelmyears, modelmaxoff);
	}
	else {
	  printf("\nPlotting an example CI off model between %.2e Hz and %.2e Hz with an on time of %d Myrs and an off time between %d and %d Myrs\n\n", minmodelfreq, maxmodelfreq, modelmyears, modelminoff, modelmaxoff);
	}
      }
      else if (model == 6) {

	if (minmodelmyears > modelmyears) {
	  printf("\nError: The minimum on time (%d Myrs) must be less than the maximum on time (%d Myrs). Please use the minmodelmyears and/or the modelmyears commands to set a new range.\n\n", minmodelmyears, modelmyears);
	  break;
	}

	if (minmodelmyears < 1) {
	  printf("\nThe on time cannot be less that 1 for a CI model as the source must be on for at least some period to be visible. Setting the minimum on time to 1 Myr.\n");
	  tmp_minmodelmyears = 0;
	}
	else {
	  tmp_minmodelmyears = minmodelmyears;
	}

	if (modelmyears < 1) {
	  printf("\nThe on time cannot be less that 1 for a CI model as the source must be on for at least some period to be visible. Setting the maximum on time to 1 Myr.\n");
	  tmp_modelmyears = 0;
	}
	else {
	  tmp_modelmyears = modelmyears;
	}

	printf("\nPlotting an example CI model between %.2e Hz and %.2e Hz with an age of between %d and %d Myrs\n\n", minmodelfreq, maxmodelfreq, minmodelmyears, modelmyears);
	tmp_modelminoff = 0;
      }

      plotcimodel(minmodelfreq, maxmodelfreq, inject, fieldstrength, model, usr_gmin, usr_gmax, tmp_modelmyears, tmp_minmodelmyears, output, titles, exampleredshift, tmp_modelminoff, modelmaxoff, skip, varyoffage);

      break;


      /*case 118: // Plot the CI Off analystical model for an arbitary normalisation (t_off = 0 in KGJP)

      // *** No longer in use ***

      printf("\nPlotting an example CI model between %.2e Hz and %.2e Hz\n\n", minmodelfreq, maxmodelfreq);

      model = 5;


      if (export == 1) {

	time( &currenttime );
	time_struct = localtime ( &currenttime );

	strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	sprintf(output,"%s/ExJPModel_%s.%s/%s", imageloc,  timebuff, imagetype, imagetype);

	printf("\nExporting as %s\n", output);

	cpgopen(output);
      }
      else {
	sprintf(output,"/xs");
      }

      plotcimodel_analytical(minmodelfreq, maxmodelfreq, inject, fieldstrength, model, usr_gmin, usr_gmax, modelmyears, output, titles, exampleredshift);

      break;*/
    

    case 119: // Fit CI models to a data set

      
      if (export  == 1) {

	// Open the directory that is currently set and check it exists
	dir = opendir(dataloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output data (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in dataloc ***\n\n", dataloc);
	  break;
	}

	dir = opendir(imageloc);

	if (dir != NULL) {
	  closedir (dir);
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	  break;
	}
	
      }

      // Check that the search range parameters are valid
      if ((model == 5) && (minoff > maxoff)) {
	printf("\nError: The minimum off time (%d Myrs) must be less than the maximum off time (%d Myrs). Please use the minoff and/or the maxoff commands to set a new range.\n\n", minoff, maxoff);
	break;
      }

      if (minmyears > myears) {
	printf("\nError: The minimum on time (%d Myrs) must be less than the maximum on time (%d Myrs). Please use the minmyears and/or the myears commands to set a new range.\n\n", minmyears, myears);
	break;
      }

      char ciheaderfilename[1024];
      char cidatafilename[1024];

      int tmp_upperlimit = 0;

      int CI_linenum = 1;

      char CI_line_storage[1024];

      char tmp_index[33], **CI_index;
      float tmp_redshift, tmp_bfield, tmp_inject, *CI_bfield, *CI_inject, tmp_frequency, tmp_flux, tmp_error;

      FILE *ciheaderfile;
      FILE *cidatafile;


      // Get the header filename
      validentry = 0;
      escape = 0;

      while (validentry == 0) {

	cmdbuffer = readline("Please enter the name and location of the header file (esc to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	strcpy(tmp_filename, cmdbuffer);

	if ( (tmp_filename == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else if ((strstr(tmp_filename, "esc") != NULL) && (strlen(tmp_filename) == 3)) {
	  printf("Escaping command...\n");
	  escape = 1;
	  break;
	}
	else if ((strstr(tmp_filename, "ls") != NULL) && (strlen(tmp_filename) == 2)) {
	  syscom("ls");
	  continue;
	}
	else {
	  if ((ciheaderfile = fopen(tmp_filename, "r")) != NULL) {
	    sprintf(ciheaderfilename,"%s", tmp_filename);
	    fclose(ciheaderfile);
	    validentry = 1;
	  }
	  else {
	    printf("\n *** Unable to locate or access the given header file. Please check the file exists and the correct permission are set. ***\n\n");
	    continue;
	  }
	}
      }

      // Break the second loop if they want to escape
      if (escape == 1) {
	break;
      }

      // Get the data filename
      validentry = 0;
      escape = 0;

      while (validentry == 0) {

	cmdbuffer = readline("Please enter the name and location of the data file (esc to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	strcpy(tmp_filename, cmdbuffer);

	if ( (tmp_filename == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else if ((strstr(tmp_filename, "esc") != NULL) && (strlen(tmp_filename) == 3)) {
	  printf("Escaping command...\n");
	  escape = 1;
	  break;
	}
	else if ((strstr(tmp_filename, "ls") != NULL) && (strlen(tmp_filename) == 2)) {
	  syscom("ls");
	  continue;
	}
	else {
	  if ((cidatafile = fopen(tmp_filename, "r")) != NULL) {
	    sprintf(cidatafilename,"%s", tmp_filename);
	    fclose(cidatafile);
	    validentry = 1;
	  }
	  else {
	    printf("\n *** Unable to locate or access the given data file. Please check the file exists and the correct permission are set. ***\n\n");
	    continue;
	  }
	}
      }

      // Break the second loop if they want to escape
      if (escape == 1) {
	break;
      }

	// Setup the naming structure for the data output file
      if (export == 1) {

	// Get the export filename
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Please enter an name for the file output (no file extension required). Entering an existing file will cause the data to be appeneded (esc to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  strcpy(tmp_filename, cmdbuffer);

	  if ( (tmp_filename == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }
	  else if ((strstr(tmp_filename, "esc") != NULL) && (strlen(tmp_filename) == 3)) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if (model == 5) {
	  strcpy(tmp_filename2, "CIOff_Fitting_Results");
	}

	else if (model == 6) {
	  strcpy(tmp_filename2, "CI_Fitting_Results");
	}
	else {
	  strcpy(tmp_filename2, "UnknownModel_Results");
	}

	if (strcmp(tmp_filename,"") != 0) { // If the user input isnt blank, append an underscore to link the text
	  strcat(tmp_filename2, "_");
	}
	

	sprintf(dataoutput,"%s/%s%s.dat", dataloc, tmp_filename2, tmp_filename);

	if (exportcompression == 0) {
	  sprintf(dataoutput_ext,"%s/%s_MODELDATA.dat", dataloc, tmp_filename);
	}

	// Check if the files already exists
	if ((filestream = fopen(dataoutput, "r")) != NULL) {

	  validentry = 0;
	  escape = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Warning: File already exists, are you sure you want to append to the end of the current data file? (1: Yes 0: No): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    confirm = strtol(cmdbuffer, &endptr, 10);

	    if ( (confirm == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }
	    else if ((confirm != 0) && (confirm != 1)) {
	      printf("\nInvalid entry, please try again...\n\n");
	      continue;
	    }
	    else if (confirm == 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else {
	      validentry = 1;
	    }
	  }

	  fclose(filestream);

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }
	}
	// Now check the model file if we are using uncompressed data export
	else if ( (exportcompression == 0) && ((filestream = fopen(dataoutput_ext, "r")) != NULL) ) {

	  validentry = 0;
	  escape = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Warning: A data file is not present with that name, but a model data file already exists. Are you sure you want to append to the end of this file? (1: Yes 0: No): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    confirm = strtol(cmdbuffer, &endptr, 10);

	    if ( (confirm == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }
	    else if ((confirm != 0) && (confirm != 1)) {
	      printf("\nInvalid entry, please try again...\n\n");
	      continue;
	    }
	    else if (confirm == 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else {
	      validentry = 1;
	      printf("\nExporting model data as %s\n", dataoutput);
	    }
	  }

	  fclose(filestream);

	  // Break the second loop if they want to escape
	  if (escape == 1) {
	    break;
	  }
	}
	printf("Exporting data as %s\n", dataoutput);
	printf("Exporting model data as %s\n", dataoutput);
      }

      // Read in the header file and determine how many entries are expect
      escape = 0;
      regnumber_CI = 0;

      printf("Reading header file...\n");

      int numlines = 0;

      if ((ciheaderfile = fopen(ciheaderfilename, "r")) != NULL) {

	// See what the maximum number of regions / targets is from the header
	while( fgets(CI_line_storage, sizeof(CI_line_storage), ciheaderfile) != NULL )  {
	  numlines++;
	}

	// Allocate memory for the header data
	CI_index = (char **)calloc(numlines+1, sizeof(char *));
	CI_redshift = (float *)calloc(numlines+1, sizeof(float));
	CI_bfield = (float *)calloc(numlines+1, sizeof(float));
	CI_inject = (float *)calloc(numlines+1, sizeof(float));

	for (i=0; i<=numlines; i++) {
	  CI_index[i] = (char *)calloc(33, sizeof(char));
	}

	rewind(ciheaderfile);

	while( fgets(CI_line_storage, sizeof(CI_line_storage), ciheaderfile) != NULL )  {

	  if (sscanf(CI_line_storage,"%100[^,],%f,%f,%f", tmp_index, &tmp_redshift, &tmp_bfield, &tmp_inject) != 4) {
	    printf("*** Warning: Missing entry or incorrect data format found on line %d. Skipping this entry... ***\n", CI_linenum);
	  }
	  else {

	    // Check if the new index is a duplicate key. We are fine to use <= here without risk of checking its own value as we have not yet assigned to the perminant array or incremented regnumber.
	    for (i=1; i<=regnumber_CI; i++) {
	      if (strcasecmp(tmp_index, CI_index[i]) == 0) {
		fprintf(stderr,"\n*** Error: Unable to determine data parameters for each region/target due to a duplicate index of %s on lines %d and %d of the header file. Please ensure each entry has a unique reference index and try again. ***\n\n", tmp_index, i, CI_linenum);
		escape = 1;
		break;
	      }
	    }

	    if (escape == 1) {
	      break;
	    }
	  
	    //Now we are confident the entry is fine, assign a new regnumber	   
	    regnumber_CI++;

	    //Copy the values over to the perminant arrays
	    strcpy(CI_index[regnumber_CI],tmp_index);
	    CI_redshift[regnumber_CI] = tmp_redshift;
	    CI_bfield[regnumber_CI] = tmp_bfield;
	    CI_inject[regnumber_CI] = tmp_inject;
	  }
	  CI_linenum++;
	}
      }
      else {
	fprintf(stderr,"*** Error: Unable to locate or access the given header file. Please check the file exists and the correct permission are set. ***\n\n");
        break;
      }

      if (escape == 1) {
	break;
      }
    
      printf("Found headers for %d regions / target sources\n", regnumber_CI);

      // Reallocate the memory for the header data now we know exactly how much we need (trim)
      CI_redshift = (float *)realloc(CI_redshift, (regnumber_CI+1) * sizeof(float));
      CI_bfield = (float *)realloc(CI_bfield, (regnumber_CI+1) * sizeof(float));
      CI_inject = (float *)realloc(CI_inject, (regnumber_CI+1) * sizeof(float));
    
      cichisquared = (float *)calloc(regnumber_CI+1, sizeof(float));
      cibestage = (float *)calloc(regnumber_CI+1, sizeof(float));
      cibestoff = (float *)calloc(regnumber_CI+1, sizeof(float));
      cimodelalpha = (float *)calloc(regnumber_CI+1, sizeof(float));
      cimodelflux = (float **)calloc(regnumber_CI+1, sizeof(float *));
      cibestnorm = (float *)calloc(regnumber_CI+1, sizeof(float));
      log_cimodelflux = (float **)calloc(regnumber_CI+1, sizeof(float *));

      cibcmb = (float *)calloc(regnumber_CI+1, sizeof(float));
      ciageingb = (float *)calloc(regnumber_CI+1, sizeof(float));
      cibreakon = (float *)calloc(regnumber_CI+1, sizeof(float));
      cibreakoff = (float *)calloc(regnumber_CI+1, sizeof(float));
      ciconflvl = (float *)calloc(regnumber_CI+1, sizeof(float));


      for (i=1; i<=regnumber_CI; i++) {
	cimodelflux[i]= (float *)calloc(modelres+1, sizeof(float));
	log_cimodelflux[i] = (float *)calloc(modelres+1, sizeof(float));
      }

      tmp_offerrorsplus_CI = (float *)calloc(regnumber_CI+1, sizeof(float));
      tmp_offerrorsminus_CI = (float *)calloc(regnumber_CI+1, sizeof(float));
		
      tmp_onerrorsplus_CI = (float *)calloc(regnumber_CI+1, sizeof(float));
      tmp_onerrorsminus_CI = (float *)calloc(regnumber_CI+1, sizeof(float));
		
      imgnum_CI = (int *)calloc(regnumber_CI+1, sizeof(int));

      numupperlimits = (int *)calloc(regnumber_CI+1, sizeof(int));
      upperlimits = (int **)calloc(regnumber_CI+1, sizeof(int *));

      // Make sure we pass at least an array of 0s for each region
      for (i=1; i<=regnumber_CI; i++) {
	upperlimits[i] = (int *)calloc(1, sizeof(int));
      }

      // The order here is reversed compared to normal JP/KP fitting. It is sent to a tmp_array for plotting later anyway and is a major pain otherwise. May require seperate spetral index etc. modules. The variables regflux_CI and fluxerror_CI require a dummy array at the end for legacy purposes.

      regflux_CI = (float ***)calloc(regnumber_CI+1, sizeof(float **));
      frequency_CI = (float **)calloc(regnumber_CI+1, sizeof(float *));
      fluxerror_CI = (float ***)calloc(regnumber_CI+1, sizeof(float **));
      fluxorder_CI = (int **)calloc(regnumber_CI+1, sizeof(int *));
  
      printf("Reading data file...\n");

      if ((cidatafile = fopen(cidatafilename, "r")) != NULL) {

	for (i=1; i<=regnumber_CI; i++) { // Doing the nested loops this way around is much easier to handle

	  CI_linenum = 1;
  
	  while( fgets(CI_line_storage, sizeof(CI_line_storage), cidatafile) != NULL )  {

	    if (sscanf(CI_line_storage,"%100[^,],%f,%f,%f,%d", tmp_index, &tmp_frequency, &tmp_flux, &tmp_error, &tmp_upperlimit) != 5) {
	      printf("Missing entry or incorrect data format found on line %d. Skipping this entry...\n", CI_linenum);
	    }
	    else {

	      if (strcasecmp(tmp_index, CI_index[i]) == 0) {

		// Once we have a match process it

		// Order has been reversed. See previous memory allocation for details. Regflux_CI and fluxerror_CI also need a dummy region array at the end for legacy purposes
		regflux_CI[i] = (float **)realloc(regflux_CI[i], (imgnum_CI[i]+1) * sizeof(float *));
		fluxerror_CI[i] = (float **)realloc(fluxerror_CI[i], (imgnum_CI[i]+1) * sizeof(float *));
		frequency_CI[i] = (float *)realloc(frequency_CI[i], (imgnum_CI[i]+1) * sizeof(float));
	        fluxorder_CI[i] = (int *)realloc(fluxorder_CI[i], (imgnum_CI[i]+1) * sizeof(int));

		regflux_CI[i][imgnum_CI[i]] = (float *)calloc(2, sizeof(float)); // Dummy array
		fluxerror_CI[i][imgnum_CI[i]] = (float *)calloc(2, sizeof(float)); // Dummy array
	      
		if (tmp_upperlimit == 1) {
		  numupperlimits[i]++;
		  upperlimits[i] = (int *)realloc(upperlimits[i], numupperlimits[i] * sizeof(int));
		  upperlimits[i][numupperlimits[i]-1] = imgnum_CI[i];
		}
		
		frequency_CI[i][imgnum_CI[i]] = tmp_frequency;
		regflux_CI[i][imgnum_CI[i]][1] = tmp_flux;
		fluxerror_CI[i][imgnum_CI[i]][1] = tmp_error;

		imgnum_CI[i]++;
	      }
	    }
	    CI_linenum++;
	  }
	  rewind(cidatafile);
	
	  if (imgnum_CI[i] == 0) {
	    printf("*** Warning: No matching data found for header index %s. Skipping this entry... ***\n", CI_index[i]);
	  }
	  else if (imgnum_CI[i] < minnummaps) {
	    printf("*** Warning: Only %d data point(s) found for header index %s (minimum of %d required). Skipping this entry... ***\n", imgnum_CI[i], CI_index[i], minnummaps);
	  }

	}

	fclose(cidatafile);
      }
      else {
	fprintf(stderr,"*** Error: Unable to locate or access the given data file. Please check the file exists and the correct permission are set. ***\n\n");
        break;
      }
    
      int rejected, success, setvalid;
      float curminfreq, curmaxfreq;

      setvalid = 0;
      rejected = 0;
      success = 0;

      // Set the flux order arrayfor each dataset
      for (i=1; i <= regnumber_CI; i++) {
	
	if (imgnum_CI[i] >= minnummaps) {

	  setvalid++; // Indicates how many header/data pairs are valid for fitting used in e.g. closing the plotting window

	  for (j=0; j < imgnum_CI[i]; j++) {

	    if (j == 0) {
	      curminfreq = 1e-32;
	    }
	    else {
	      curminfreq = frequency_CI[i][fluxorder_CI[i][j-1]];
	    }
	    curmaxfreq = 1e+32;

	    for (a=0; a < imgnum_CI[i]; a++) {

	      if (a==0) {
		fluxorder_CI[i][j] = a;
	      }
	      if ((frequency_CI[i][a] > curminfreq) && (frequency_CI[i][a] < curmaxfreq)) {
		curmaxfreq = frequency_CI[i][a];
		fluxorder_CI[i][j] = a;
	      }
	    }
	  }
	}
      }

    
      if (model == 5) {
	if (minmyears < 1) {
	  printf("\nThe on time cannot be exactly 0 for a CI off model as the source must be on for at least some period to be visible. Setting the minimum search range to 0.01 Myr.\n\n");
	  tmp_minmyears = 0;
	}
	else {
	  tmp_minmyears = minmyears;
	}

	if (myears < 1) {
	  printf("\nThe maximum on time search range cannot be less that 1 for a CI off model as the source must be on for at least some period to be visible. Setting the maximum on time to 1 Myr.\n\n");
	  tmp_myears = 1;
	}
	else {
	  tmp_myears = myears;
	}
      }
      else if (model == 6) {
	if (minmyears < 1) {
	  printf("\nThe on time cannot be exactly 0 for a CI model as the source must be on for at least some period to be visible. Setting the minimum search range to 0.01 Myr.\n\n");
	  tmp_minmyears = 0;
	}
	else {
	  tmp_minmyears = minmyears;
	}
      
	if (myears < 1) {
	  printf("\nThe maximum on time search range cannot be less that 1 for a CI model as the source must be on for at least some period to be visible. Setting the maximum on time to 1 Myr.\n\n");
	  tmp_myears = 1;
	}
	else {
	  tmp_myears = myears;
	}
      }
      else {
	printf("\n *** Error: Unable to determine model type. Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org *** \n\n");
	break;
      }
    
      int regcounter=0;
      int hundreds=0;
      float tmp_ymin, tmp_ymax;
      int skipcounter;


      printf("Calculating model fits (this may take some time)\n Currently processing region: 1 of %d valid regions / targets...\n", setvalid);

      skipcounter=1;

      for (i=1; i<=regnumber_CI; i++) {

	if (imgnum_CI[i] >= minnummaps) {

	  // If we are on a skip round and not exporting, increment i by 1 and try again
	  if  ( (skipcounter != skip) && (export != 1) ) {
	    skipcounter++;
	    continue;
	  }

	  tmp_ymin=1e32;
	  tmp_ymax=1e-32;

	  regcounter++;

	  if (regcounter >= 100) {
	    hundreds++;
	    printf("%d00...\n", hundreds);
	    regcounter = 0;
	  }

	  for (j=0; j<imgnum_CI[i]; j++) {
	    frequency_CI[i][j] *= (1+CI_redshift[i]);
	  }

	  // Check by return if the fit is valid

	  int suppress; // Suppresses being kicked out if we hit the minimum fitting flux. Always 1 or 0.
	  int tmp_suppress; // If > 1, something was suppressed
	
	  suppress = 1;
	  tmp_suppress = suppress; // The relation between suppress and tmp_suppress is currently redundant, but it makes it easier if we want to use this elsewhere later on as we dont change suppress in the global scope


	  // This is for Leah's modelling. Make usable by the general public later.
	  int fixedage = 0; // Switch on/off
	  float tmp_fixedage = -1.0; // Value placeholder

	  //yyyyyyy
	  // This is for Leah's modelling. Make usable by the general public later.
	  if (fixedage == 1) {
	    char ciagefilename[1024];
	    FILE *ciagefile;
	    int foundentry = 0;

	    sprintf(ciagefilename, "nopeaked_sample_highz_ages.dat");


	    if ((ciagefile = fopen(ciagefilename, "r")) != NULL) {

	      CI_linenum = 1;
  
	      while( fgets(CI_line_storage, sizeof(CI_line_storage), ciagefile) != NULL )  {

		if (sscanf(CI_line_storage,"%100[^,],%f", tmp_index, &tmp_fixedage) != 2) {
		  printf("Missing entry or incorrect age format found on line %d. Skipping this entry...\n", CI_linenum);
		}
		else {//yyyyy
		  if (strcasecmp(tmp_index, CI_index[i]) == 0) {
		    //printf("index: %s fixedage: %.4f\n",CI_index[i], tmp_fixedage);
		    foundentry = 1;
		    break;
		  }
		}
		CI_linenum++;
	      }

	      rewind(ciagefile);
	
	      if (foundentry != 1) {
		printf("*** Warning: No matching age data found for header index %s. Skipping this entry... ***\n", CI_index[i]);
		fclose(ciagefile);
		continue;
	      }
	
	      fclose(ciagefile);
	    }
	    else {
	      fprintf(stderr,"*** Error: Unable to locate or access the given age file. Please check the file exists and the correct permission are set. ***\n\n");
	      break;
	    }

	  }
	
	  printf("Fitting for source with index %s\n", CI_index[i]);

	  if(ciageingmodels(regflux_CI[i], imgnum_CI[i], fluxorder_CI[i], frequency_CI[i], i, cichisquared, cibestage, cibestoff, fluxerror_CI[i], cimodelalpha, cimodelflux, usr_gmin, usr_gmax, tmp_myears, tmp_minmyears, minoff, maxoff, ageresolution, levels, cibestnorm, modelres, printresults, CI_inject[i], CI_bfield[i], model, printreject, CI_redshift[i], tmp_offerrorsplus_CI, tmp_offerrorsminus_CI, tmp_onerrorsplus_CI, tmp_onerrorsminus_CI, numupperlimits[i], upperlimits[i], cibcmb, ciageingb, cibreakon, cibreakoff, ciconflvl, &tmp_suppress, CI_index[i], export, suppresscdf, extrapolatemodel, extrapolationfrequency, tmp_fixedage) != 5) {

	    modelfreq_CI = (float *)calloc(modelres+1, sizeof(float));
	    fluxfreq_CI = (float *)calloc(imgnum_CI[i], sizeof(float));
	    fluxerrorplus_CI = (float *)calloc(imgnum_CI[i], sizeof(float));
	    fluxerrorminus_CI = (float *)calloc(imgnum_CI[i], sizeof(float));
	    tmp_regflux_CI = (float *)calloc(imgnum_CI[i], sizeof(float));
	    sprintf(citarget, "%s", CI_index[i]);
	
	    for (j=0; j<=modelres; j++) {

	      if ( (extrapolatemodel == 1) && (extrapolationfrequency > frequency_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]])) {
		modelfreq_CI[j] = log10(frequency_CI[i][fluxorder_CI[i][0]] + ( ( (extrapolationfrequency - frequency_CI[i][fluxorder_CI[i][0]]) / modelres ) * j ) );
	      }
	      else if ( (extrapolatemodel == 1) && (extrapolationfrequency < frequency_CI[i][fluxorder_CI[i][0]]) ){	
		modelfreq_CI[j] = log10(extrapolationfrequency + ( ( (frequency_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]] - extrapolationfrequency) / modelres ) * j ) );
	      }
	      else {
		modelfreq_CI[j] = log10(frequency_CI[i][fluxorder_CI[i][0]] + ( ( (frequency_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]] - frequency_CI[i][fluxorder_CI[i][0]]) / modelres ) * j ) );
	      }

	      //modelfreq_CI[j] = log10(frequency_CI[i][fluxorder_CI[i][0]] + ( ( (frequency_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]] - frequency_CI[i][fluxorder_CI[i][0]]) / modelres ) * j ) );
	      log_cimodelflux[i][j] = log10(cimodelflux[i][j]);
	    }
	
	    for (j=0; j<imgnum_CI[i]; j++) {
	      fluxfreq_CI[j] = log10(frequency_CI[i][j]);
	      fluxerrorplus_CI[j] = log10(regflux_CI[i][j][1] + fluxerror_CI[i][j][1]);
	      fluxerrorminus_CI[j] = log10(regflux_CI[i][j][1] - fluxerror_CI[i][j][1]);
	      tmp_regflux_CI[j] = log10(regflux_CI[i][j][1]);
	    }

	    // If exporting, output open a new session for each with a different name
	    if (export == 1) {

	      time( &currenttime );
	      time_struct = localtime ( &currenttime );

	      strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	      if (model == 5) {
		sprintf(output,"%s/%s_CIOFF_%s.%s/%s", imageloc, citarget, timebuff, imagetype, imagetype);
	      }
	      else  if (model == 6) {
		sprintf(output,"%s/%s_CI_%s.%s/%s", imageloc, citarget, timebuff, imagetype, imagetype);
	      }
	      else {
		sprintf(output,"%s/%s_UnknownModel_%s.%s/%s", imageloc, settarget[currentset], timebuff, imagetype, imagetype);
	      }
	      // printf("Exporting image as %s\n", output);
	    }
	    else {
	      sprintf(output,"/xs");
	    }

	    cpgopen(output);

	    if (model == 5) { // CI Off Model
	      freeparams = 3;
	    }
	    else if (model == 6) { // CI Model
	      freeparams = 2;
	    }
	    else {
	      fprintf(stderr,"\nError: Model number is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");
	      break;
	    }
	    
	    if (usereduced == 1) {
	      tmp_chisquared = cichisquared[1] / (imgnum_CI[i]-numupperlimits[i]-freeparams);
	    }
	    else {
	      tmp_chisquared = cichisquared[1];
	    }
	  
	    if ((regflux_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]][1] - fluxerror_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]][1]) <= cimodelflux[i][modelres]) {

	      if ( (numupperlimits[i] > 0) && (upperlimits[i][numupperlimits[i]-1] == fluxorder_CI[i][imgnum_CI[i]-1]) ) { // If the last frequency is an upper limit 
		tmp_ymin = regflux_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]][1] - (regflux_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]][1]*0.3);
	      }
	      else {
		tmp_ymin = regflux_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]][1] - fluxerror_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]][1];

		/*if (tmp_ymin <= 0.0) {
		  fluxerror_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]][1] += ( (fluxerror_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]][1] - regflux_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]][1]) + 1e-32); //Stop us ever going negative
		}*/
	      }
	    }
	    else {
	      tmp_ymin = cimodelflux[i][modelres];
	    }

	    if ((regflux_CI[i][fluxorder_CI[i][0]][1] + fluxerror_CI[i][1][1]) > cimodelflux[i][0]) {
	      tmp_ymax = regflux_CI[i][0][1] + fluxerror_CI[i][0][1];
	    }
	    else {
	      tmp_ymax = cimodelflux[i][0];
	    }

	    // Catch a negative flux
	    /*if (tmp_ymin <= 0.0) {
	      tmp_ymin = 1e-32;
	      }*/

	    // Creating the plots

	    // Print out if we have had to suppress anything
	    if ( (tmp_suppress > 1) && (export !=1) ) {
	      	  printf("*** Warning: Could not fit over the entire requested search range (model flux effectively 0 or unphysically high). If you are fitting over a very large search range this should be fine, but the results should be checked. ***\n");
	    }
	
	    // Reverse the regflux_CI etc nested array order here, but only for the region we are currently interested in. This is for backwards compatibility reasons (the alternative is a nightmare)

	    if ( (extrapolatemodel == 1) && (extrapolationfrequency > frequency_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]]) ) {
	      plotmodelvsflux(log_cimodelflux[i], tmp_regflux_CI, fluxerrorplus_CI, fluxerrorminus_CI, symbol, modelres, imgnum_CI[i], fluxfreq_CI,  modelfreq_CI, autoscalex, autoscaley, userminx, usermaxx, userminy, usermaxy, frequency_CI[i][fluxorder_CI[i][0]], extrapolationfrequency, tmp_ymin, tmp_ymax, citarget, titles, labels, axismarks, CI_inject[i], CI_bfield[i], cibestage[1], cibestoff[1], tmp_chisquared, paramlabels, model, usereduced, numupperlimits[i], upperlimits[i], largetxt);
	      }
	      else if ( (extrapolatemodel == 1) && (extrapolationfrequency < frequency_CI[i][fluxorder_CI[i][0]]) ){	
		plotmodelvsflux(log_cimodelflux[i], tmp_regflux_CI, fluxerrorplus_CI, fluxerrorminus_CI, symbol, modelres, imgnum_CI[i], fluxfreq_CI,  modelfreq_CI, autoscalex, autoscaley, userminx, usermaxx, userminy, usermaxy, extrapolationfrequency, frequency_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]], tmp_ymin, tmp_ymax, citarget, titles, labels, axismarks, CI_inject[i], CI_bfield[i], cibestage[1], cibestoff[1], tmp_chisquared, paramlabels, model, usereduced, numupperlimits[i], upperlimits[i], largetxt);
	      }
	      else {
		plotmodelvsflux(log_cimodelflux[i], tmp_regflux_CI, fluxerrorplus_CI, fluxerrorminus_CI, symbol, modelres, imgnum_CI[i], fluxfreq_CI,  modelfreq_CI, autoscalex, autoscaley, userminx, usermaxx, userminy, usermaxy, frequency_CI[i][fluxorder_CI[i][0]], frequency_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]], tmp_ymin, tmp_ymax, citarget, titles, labels, axismarks, CI_inject[i], CI_bfield[i], cibestage[1], cibestoff[1], tmp_chisquared, paramlabels, model, usereduced, numupperlimits[i], upperlimits[i], largetxt);
	      }


	    // Exporting the data
	    if (export == 1) {

	      if ( (filestream = fopen(dataoutput,"a") ) == NULL ) {
		printf("\n*** Error: Cannot open specified file to output data (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in dataloc ***\n\n", dataoutput);
		break;
	      }
	      else {

		// Check if we have suppressed anything and put in in a more intuative 0/1 format for output
		if (tmp_suppress > 1) {
		  tmp_suppress = 1;
		}
		else {
		  tmp_suppress = 0;
		}

		// At the moment, most of the data is in array [1] which is then overwritten for each iteration. This is rather wasteful as the array sizes are of order [regnum+1]. This does allow us to store the data later if we wish though. The true region number is passed across to ciageing but requires a rework of that module if thats what we want to do. Altervatively, use a tmp_ and pass it once we are back here. Best way forward is still TBD, but this works for now as we are not memory intensive at this point anyway.
		fprintf(filestream, "%s,%.4f,%.4e,%.4f,%d,%d,%d,%.4e,%.4e,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4e,%.4f,%d,%.4e,%.4e,%.2e,%.2e,%d", CI_index[i], CI_redshift[i], CI_bfield[i], CI_inject[i], imgnum_CI[i], numupperlimits[i], model, cibcmb[i], ciageingb[i], (cibestage[1] + cibestoff[1]), sqrt(pow(tmp_onerrorsplus_CI[1],2)+pow(tmp_offerrorsplus_CI[1],2)), sqrt(pow(tmp_onerrorsminus_CI[1],2)+pow(tmp_offerrorsminus_CI[1],2)), cibestage[1], tmp_onerrorsplus_CI[1], tmp_onerrorsminus_CI[1], cibestoff[1], tmp_offerrorsplus_CI[1], tmp_offerrorsminus_CI[1], cichisquared[1], (cichisquared[1] / (imgnum_CI[i]-numupperlimits[i]-freeparams)), (imgnum_CI[i] - numupperlimits[i] - freeparams), cibestnorm[1], cibreakon[i], cibreakoff[i], ciconflvl[i], tmp_suppress);

		fprintf(filestream, "\n");
		fclose(filestream);
	      }
	    
	      // If we are exporting the full data set with model fluxes
	      if (exportcompression == 0) {

		if ( (filestream = fopen(dataoutput_ext,"a") ) == NULL ) {

		  printf("\n*** Error: Cannot open specified file to output the model data (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in dataloc ***\n\n", dataoutput);
		  break;
		}
	      
		else {

		  for (j=0; j <= modelres; j++) {

		    if ( (extrapolatemodel == 1) && (extrapolationfrequency > frequency_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]])){
		    fprintf(filestream, "%s,%.4e,%.4e", CI_index[i], (frequency_CI[i][fluxorder_CI[i][0]] + ( ( (extrapolationfrequency - frequency_CI[i][fluxorder_CI[i][0]]) / modelres ) * j)), cimodelflux[i][j]);
		    }
		    else if ( (extrapolatemodel == 1) && (extrapolationfrequency < frequency_CI[i][fluxorder_CI[i][0]]) ){
		      fprintf(filestream, "%s,%.4e,%.4e", CI_index[i], (extrapolationfrequency + ( ( (frequency_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]] - extrapolationfrequency) / modelres ) * j)), cimodelflux[i][j]);
		    }
		    else {
		      fprintf(filestream, "%s,%.4e,%.4e", CI_index[i], (frequency_CI[i][fluxorder_CI[i][0]] + ( ( (frequency_CI[i][fluxorder_CI[i][imgnum_CI[i]-1]] - frequency_CI[i][fluxorder_CI[i][0]]) / modelres ) * j)), cimodelflux[i][j]);
		    }

		    //fprintf(filestream, "%s,%.4e,%.4e", CI_index[i], (frequency_CI[i][fluxorder_CI[i][0]] + ( ( (frequency_CI[fluxorder_CI[i][imgnum_CI[i]-1]] - frequency_CI[fluxorder_CI[i][0]]) / modelres ) * j)), cimodelflux[i][j]);
		    fprintf(filestream, "\n");
		  }
		    fclose(filestream);
		}
	      }
	    }

	    if ( (skipcounter == skip) && (export != 1) && (i != regnumber_CI) ) {
	      printf("Model fitting for source with index %s complete\n\n", CI_index[i]);
	      printf("Press enter to continue...\n");
	      getchar();
	      skipcounter=0;
	      cpgclos();
	    }
	    else if (i == regnumber_CI) {
	      skipcounter=0;
	      cpgclos();
	    }
	    else if (export == 1) {
	      cpgclos();
	    }
	    skipcounter++;
	    success++;
	  }
	  else {
	    rejected++;
	  }
	}
      }

      if (setvalid != 0) {
	free(fluxfreq_CI);

        if (skip > 1) { // Warn if we have skipped things
	  printf("\n*** Warning: Skip is currently set to %d. Fitting has not be attempted on all possible data. ***\n", skip);
	}

	printf("\nModel fitting completed sucessfully (Attempted to fit %d of %d records. Succeeded: %d Failed: %d)\n\n",  (success+rejected), regnumber_CI, success, rejected);
      }
      else {
	printf("\n*** Error: No valid header / data pairs found. Unable to perform model fitting. ***\n\n");
      }

      // Free up the memory
      free(numupperlimits);
      free(tmp_onerrorsplus_CI);
      free(tmp_onerrorsminus_CI); 
      free(tmp_offerrorsplus_CI);
      free(tmp_offerrorsminus_CI); 
      free(modelfreq_CI);
      free(regflux_CI);
      free(frequency_CI);
      free(cichisquared);
      free(cibestage);
      free(cibestoff);
      free(cimodelalpha);
      free(fluxerror_CI);
      free(log_cimodelflux);
      free(cimodelflux);
      free(cibestnorm);
      free(fluxorder_CI);
      free(fluxerrorplus_CI);
      free(fluxerrorminus_CI);
      free(tmp_regflux_CI);
      free(upperlimits);
      free(CI_redshift);
      free(cibcmb);
      free(ciageingb);
      free(cibreakon);
      free(cibreakoff);
      free(ciconflvl);

      break;
      

    case 120: // Set minimum on age for the CI models (must be > 0)

      validentry = 0;
      tmp_minmodelmyears = minmodelmyears;

      printf("Minimum example model age/on time currently set to %d Myr (DEFAULT: %d Myr)\n", minmodelmyears, DEFAULTMINMODELMYEARS);

      while (validentry == 0) {

	cmdbuffer = readline("Enter the minimum age/on time in megayears to plot when outputting models (e.g. plotcimodel, plotjpmodel) (integer >= 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	minmodelmyears = strtol(cmdbuffer, &endptr, 10);

	if ( (minmodelmyears == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (minmodelmyears < 0) {
	    printf("Escaping command...\n");
	    minmodelmyears = tmp_minmodelmyears;
	    break;
	  }
	  else {
	    printf("\nThe minimum on time / age for the example models is now set to %d Myrs\n\n", minmodelmyears);
	    validentry = 1;
	  }
	}
      }

      break;


    case 121: // Set minimum off time for the CI models (must be > 0)

      validentry = 0;
      tmp_modelminoff = modelminoff;

      printf("Minimum example model off time currently set to %d Myr (DEFAULT: %d Myr)\n", modelminoff, DEFAULTMODELMINOFF);


      while (validentry == 0) {

	cmdbuffer = readline("Enter the minimum off time in megayears to plot when outputting example CI off models (e.g. plotcioff) (integer >= 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	modelminoff = strtol(cmdbuffer, &endptr, 10);

	if ( (modelminoff == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (modelminoff < 0) {
	    printf("Escaping command...\n");
	    modelminoff = tmp_modelminoff;
	    break;
	  }
	  else {
	    printf("\nThe minimum off time for the example CI model is now set to %d Myrs\n\n", modelminoff);
	    validentry = 1;
	  }
	}
      }

      break;


    case 122: // Set maximum age off for the CI models (must be > 0)

      validentry = 0;
      tmp_modelmaxoff = modelmaxoff;

      printf("Maximum example model off time currently set to %d Myr (DEFAULT: %d Myr)\n", modelmaxoff, DEFAULTMODELMAXOFF);

      while (validentry == 0) {

	cmdbuffer = readline("Enter the maximum off time in megayears to plot when outputting example CI off models (e.g. plotcioff) (integer > 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	modelmaxoff = strtol(cmdbuffer, &endptr, 10);

	if ( (modelmaxoff == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (modelmaxoff < 0) {
	    printf("Escaping command...\n");
	    modelmaxoff = tmp_modelmaxoff;
	    break;
	  }
	  else if (modelmaxoff < 1) {
	    printf("\nMaximum off time cannot be less that 1 Myr. Please use the plotcimodel command for a spectrum with no off component.\n\n");
	    continue;
	  }
	  else {
	    printf("\nThe maximum off time for the example CI model is now set to %d Myrs\n\n", modelmaxoff);
	    validentry = 1;
	  }
	}
      }

      break;


    case 123: // Switch between varying off and on ages in example CI models

      if (varyoffage == 1) {
	varyoffage = 0;
	printf("CI off example models will now vary the on time between ciminmodelmyears and modelmyears (currently: %d to %d Myrs).\n", minmodelmyears, myears);
      }
      else if (varyoffage == 0) {
	varyoffage = 1;
	printf("CI off example models will now vary the off time between modelminoff and modelmaxoff (currently: %d to %d Myrs).\n", modelminoff, modelmaxoff);
      }
      else {
	printf("\nThe value of export is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", export);
      }

      break;


    case 124: // Set the redshift to use when plotting ageing models

      validentry = 0;
      tmp_exampleredshift = exampleredshift;

      printf("Redshift for example models currently set to %.3f (DEFAULT: %.3f)\n", exampleredshift, DEFAULTEXAMPLEREDSHIFT);

      while (validentry == 0) {

	cmdbuffer = readline("Enter the redshift to use when plotting ageing models (e.g. plotcioff, plotjpmodel etc.). (integer >= 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	exampleredshift = strtod(cmdbuffer, &endptr);

	if ( (exampleredshift == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (exampleredshift < 0) {
	    printf("Escaping command...\n");
	    exampleredshift = tmp_exampleredshift;
	    break;
	  }
	  else {
	    printf("\nThe redshift to be used in example spectral ageing model plots is now set to %.3f\n\n", exampleredshift);
	    validentry = 1;
	  }
	}
      }

      break;


    case 125: // Set minimum off time for the CI models (must be >= 0)

      validentry = 0;
      tmp_minoff = minoff;

      printf("Minimum off time for CI models currently set to %d Myr (DEFAULT: %d Myr)\n", minoff, DEFAULTMINOFF);

      while (validentry == 0) {

	cmdbuffer = readline("Enter the minimum off time in megayears to use when fitting CI off models (e.g. fitcioff). (integer >= 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	minoff = strtol(cmdbuffer, &endptr, 10);

	if ( (minoff == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (minoff < 0) {
	    printf("Escaping command...\n");
	    minoff = tmp_minoff;
	    break;
	  }
	  else {
	    printf("\nThe minimum off time to be used when fitting CI models is now set to %d Myrs\n\n", minoff);
	    validentry = 1;
	  }
	}
      }

      break;


    case 126: // Set maximum age off for the CI s (must be > 0)

      validentry = 0;
      tmp_maxoff = maxoff;

      printf("Maximum off time for CI models currently set to %d (DEFAULT: %d)\n", maxoff, DEFAULTMAXOFF);

      while (validentry == 0) {

	cmdbuffer = readline("Enter the maximum off time in megayears to use when fitting CI off models (e.g. fitcioff). (integer > 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	maxoff = strtol(cmdbuffer, &endptr, 10);

	if ( (maxoff == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (maxoff < 0) {
	    printf("Escaping command...\n");
	    maxoff = tmp_maxoff;
	    break;
	  }
	  else if (maxoff < 1) {
	    printf("\nMaximum off time cannot be less that 1 Myr. Please use the fitcimodel command a spectrum with no off component.\n\n");
	    continue;
	  }
	  else {
	    printf("\nThe maximum off time to be used when fitting CI off models is now set to %d Myrs\n\n", maxoff);
	    validentry = 1;
	  }
	}
      }

      break;


    case 127: // Set minimum on time / age for all ageing models

      validentry = 0;
      tmp_minmyears = minmyears;

      printf("Minimum age currently set to %d Myr (DEFAULT: %d Myr)\n", minmyears, DEFAULTMINMYEARS);

      while (validentry == 0) {

	cmdbuffer = readline("Enter the minimum age in megayears to attempt when fitting spectral ageing models (e.g. fitjpmodel, fitcioff) (integer >= 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	minmyears = strtol(cmdbuffer, &endptr, 10);

	if ( (minmyears == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (minmyears < 0) {
	    printf("Escaping command...\n");
	    minmyears = tmp_minmyears;
	    break;
	  }
	  else {
	    printf("\nThe minimum age to attempt when fitting models is now set to %d Myrs\n\n", minmyears);
	    validentry = 1;
	  }
	}
      }

      break;


    case 128: // Export model data for an arbitary normalisation

      validentry = 0;
      exactage = 0;
      escape = 0;

      if (export == 1) {

	// Open the directory that us currently set and check it exists
	dir = opendir(dataloc);

	if (dir != NULL) {
	  closedir (dir); // If we find it move on without comment
	}
	else {
	  printf("\n*** Error: Cannot open specified directory to output data (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in dataloc ***\n\n", dataloc);
	  break;
	}
      }
      

      while (validentry == 0) {

	cmdbuffer = readline("Output a range of ages between minmodelmyears and modelmyears (0) or an exact age (1)? (integer, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	exactage = strtol(cmdbuffer, &endptr, 10);

	if ( (exactage == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (exactage < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if (exactage > 1) {
	    printf("\nInvalid selection. Must be either 0 (range) or 1 (exact), please try again...\n\n");
	  }
	  else {
	    validentry = 1;
	  }
	}
      }

      if (escape == 1) {
	break;
      }

      validentry = 0;
      escape = 0;

      if (exactage == 1) {
	while (validentry == 0) {

	  if ( (model <= 4) || (model == 6) ) {
	    cmdbuffer = readline("Please enter an age to be output in Myr? (float, -1 to escape): ");
	  }
	  else if (model == 5) {
	    cmdbuffer = readline("Please enter the on time to be output in Myr? (float, -1 to escape): ");
	  }
	  else {
	    cmdbuffer = readline("Please enter an age to be output in Myr? (float, -1 to escape): ");
	  }

	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  exactmyears = strtod(cmdbuffer, &endptr);

	  if ( (exactmyears == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }
	  else {
	    if (exactmyears < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else {
	      validentry = 1;
	    }
	  }
	}

	if (model == 5) {

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Please enter the off time to be output in Myr? (float, -1 to escape): ");

	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    exactoff = strtod(cmdbuffer, &endptr);

	    if ( (exactoff == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }
	    else {
	      if (exactoff < 0) {
		printf("Escaping command...\n");
		escape = 1;
		break;
	      }
	      else {
		validentry = 1;
	      }
	    }
	  }
	}
      }
      else {
	exactmyears = 0.0;
      }

      if (escape == 1) {
	break;
      }

      printf("\nCreating example model data between %.2e Hz and %.2e Hz\n\n", minmodelfreq, maxmodelfreq);

      // Setup the output name and path

      time( &currenttime );
      time_struct = localtime ( &currenttime );
      
      strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

      if (model == 1) {      
	sprintf(output,"%s/JP_ModelData_%s.txt", dataloc, timebuff);
      }
      if (model == 2) {      
      sprintf(output,"%s/KP_ModelData_%s.txt", dataloc, timebuff);
      }
      if (model == 3) {      
      sprintf(output,"%s/Tribble_ModelData_%s.txt", dataloc, timebuff);
      }
      if (model == 4) {      
      sprintf(output,"%s/KP_TribbleModelData_%s.txt", dataloc, timebuff);
      }
      if (model == 5) {      
      sprintf(output,"%s/CIOff_ModelData_%s.txt", dataloc, timebuff);
      }
      if (model == 6) {      
      sprintf(output,"%s/CI_ModelData_%s.txt", dataloc, timebuff);
      }

      printf("Exporting as %s\n", output);

      if (model <= 4) {      
	if (modeldata(minmodelfreq, maxmodelfreq, inject, fieldstrength, model, usr_gmin, usr_gmax, minmodelmyears, modelmyears, output, exampleredshift, exactage, exactmyears, dataintervals) == 0) {
	  printf("\nModel data sucessfully exported\n\n");
	}
      }
      else if ( (model == 5) || (model == 6) ) {      
	if (cimodeldata(minmodelfreq, maxmodelfreq, inject, fieldstrength, model, usr_gmin, usr_gmax, minmodelmyears, modelmyears, output, exampleredshift, modelminoff, modelmaxoff, skip, varyoffage, exactage, exactmyears, exactoff, dataintervals) == 0) {
	  printf("\nModel data sucessfully exported\n\n");
	}
      }
      else {
	fprintf(stderr,"\nError: Model number is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");
	break;
      }

      model = 1; // Revert back to model 1 for consistency

      break;


    case 129: // Switch between exporting maps as fits and standard images

      if (exportfits == 1) {
	exportfits = 0;
	printf("Maps will now be exported in standard image format as set by imagetype (currenty: %s)\n", imagetype);
      }
      else if (exportfits == 0) {
	exportfits = 1;
	printf("Maps will now be exported in FITS format\n");
      }
      else {
	printf("\nThe value of export is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", export);
      }

      break;


  case 130: // Set the number of intervals used when determining the spectrum of example ageing data sets 

      validentry = 0;
      tmp_dataintervals = dataintervals;

      printf("The nunber of intervals for example data sets is currently set to %d (DEFAULT: %d)\n", dataintervals, DEFAULTDATAINTERVALS);

      while (validentry == 0) {

	cmdbuffer = readline("Enter the number of data intervals to use when determinine example ageing model data (e.g. jpdata). Total output is N + 1 datapoints (integer >= 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	dataintervals= strtol(cmdbuffer, &endptr, 10);

	if ( (dataintervals == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (dataintervals < 0) {
	    printf("Escaping command...\n");
	    dataintervals = tmp_dataintervals;
	    break;
	  }
	  else {
	    printf("\nThe number of data intervals to use when determinine example ageing model data is now set to %d\n\n", dataintervals);
	    validentry = 1;
	  }
	}
      }

      break;


    case 131: // Set the target name for a data set

      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");
	// Check at least 1 dataset is loaded
	if (numdatasets > 0) {

	  // Loop through each of the sets
	  for (i=0;i<numdatasets;i++) {

	    list(i, imgnum[i], setname[i], setreg[i], setbg[i]);
	  }
	}
	else {
	  fprintf(stderr," Error: No datasets have yet been loaded!\n");
	}

	printf("\n========================================================================\n");


	tmp_currentset = currentset;

	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the number of the data set you wish to change the target/data set name for (-1 to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    currentset = tmp_currentset;
	    escape = 1;
	    break;
	  }
	  else if (currentset >= numdatasets) { // Check the set exists
	    printf("\nThis value is greater than the number of data sets that exists! Please try again...\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}
      
	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}


	validentry = 0;
	strcpy(tmp_maptarget, settarget[currentset]);


	while (validentry == 0) {

	  cmdbuffer = readline("Please enter a new name for the target/data set (esc to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  strcpy(settarget[currentset], cmdbuffer);

	  if ( (settarget[currentset] == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }
	  else {
	    if ( (strstr(settarget[currentset], "esc") != NULL) && (strlen(settarget[currentset]) == 3) ) {
	      strcpy(settarget[currentset], tmp_maptarget);
	      printf("Escaping command...\n");
	      break;
	    }
	    else {
	      printf("\nThe name of data set %d is now set to %s\n\n", currentset, settarget[currentset]);
	      validentry = 1;
	    }
	  }
	}
      }
      else {
	fprintf(stderr," Error: No data sets have yet been loaded!\n");
      }


      break;


    case 132: // Set the values for RMS noise manually

      changearange = 0;
      doallmaps = 0;

      if (numdatasets > 0) {

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to set the RMS noise (-1 to escape, integers only, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}
      
	// All maps or just one?
	validentry = 0;
     
	if (currentset == 888) {
	  printf("================================================================================================\n\n");
	  printf(" Dataset |  Index  |  Frequency (obs.)  |  Frequency (Rest)  |  RMS (Jy/Beam)  |  RMS (Jy/Pixel)  \n\n");

	  // Loop through each of the sets
	  for (i=0; i<numdatasets; i++) {
	    for (a=0; a<imgnum[i]; a++) {

	      printf("   %d          %d         %.2e Hz          %.2e Hz          %.2e          %.2e\n", i, fluxorder[i][a], frequencyobs[i][fluxorder[i][a]], frequency[i][fluxorder[i][a]], rms[i][fluxorder[i][a]]*beamarea[i], rms[i][fluxorder[i][a]]);

	    }
	  }
	  printf("\n================================================================================================\n");

	}

	else {

	  printf("=====================================================================================\n\n");
	  printf(" Index  |  Frequency (Obs.)  |  Frequency (Rest)  |  RMS (Jy/Beam)  |  RMS (Jy/Pixel)  \n\n");

	  // Loop through each of the sets
	  for (a=0; a<imgnum[currentset]; a++) {

	    printf("   %d         %.2e Hz          %.2e Hz          %.2e          %.2e\n", fluxorder[currentset][a], frequencyobs[currentset][fluxorder[currentset][a]], frequency[currentset][fluxorder[currentset][a]], rms[currentset][fluxorder[currentset][a]]*beamarea[currentset], rms[currentset][fluxorder[currentset][a]]);

	  }

	  printf("\n=====================================================================================\n");
    
	}
    

	while (validentry == 0) {

	  cmdbuffer = readline("Which map(s) would you like to change the RMS noise for? (-1 to escape, 888 all maps, 777 to select a frequency range): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  mapselect = strtol(cmdbuffer, &endptr, 10);

	  if ( (mapselect == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (mapselect < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }

	  else if ( (mapselect >= imgnum[currentset]) && (currentset != 888) && (mapselect != 888) && (mapselect != 777)) { // Check the map exists
	    printf("\nInvalid selection, please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset == 888) && (mapselect != 888) && (mapselect != 777) ){ // Check the map exists
	    count = 0;
	    for (i=0; i<numdatasets; i++) {
	      if (mapselect >= imgnum[i]) {
		printf("\n*** Warning: Data set %d does not contain an image with index %d. RMS will not be updated for this data set... ***\n", i, mapselect);
	      }
	      else {
		count++;
	      }
	    }
	    if (count == 0) {
	      printf("\nError: No loaded data set contains a map with index %d, please try again...\n\n", mapselect);
	      continue;
	    }
	    else {
	      validentry = 1;
	    }
	  }
	  else if (mapselect == 888) {
	    doallmaps = 1;
	    validentry = 1;
	  }
	  else if (mapselect == 777) {
	    changearange = 1;
	    doallmaps = 1;
	    validentry = 1;
	  }
	  else {
	    validentry = 1;
	  }
	}
      
	if (escape == 1) {
	  break;
	}

	// Get the frequency range to change if required

	if (mapselect == 777) {

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Greater than (0) or less than (1) a frequency? (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    gtlt = strtol(cmdbuffer, &endptr, 10);

	    if ( (gtlt == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (gtlt < 0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }

	    else if (gtlt > 1) { // Check the entry is valid
	      printf("\nInvalid selection, please try again...\n\n");
	      continue;
	    }

	    else {
	      validentry = 1;
	    }
	  }
      
	  if (escape == 1) {
	    break;
	  }


	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Enter the (observed) frequency above or below which the RMS noise should be applied in Hz (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    errorchange = strtod(cmdbuffer, &endptr);

	    if ( (errorchange == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (errorchange < 0.0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else {
	      validentry = 1;
	    }
	  }
      
	  if (escape == 1) {
	    break;
	  }

	  validentry = 0;

	  while (validentry == 0) {
    
	    cmdbuffer = readline("Enter the RMS noise for this selection in units of Jy/Beam (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    flaterror = strtod(cmdbuffer, &endptr);

	    if ( (flaterror == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (flaterror < 0.0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else {
	      validentry = 1;
	    }
	  }
      
	  if (escape == 1) {
	    break;
	  }

	}

	else {

	  validentry = 0;

	  while (validentry == 0) {

	    cmdbuffer = readline("Enter the RMS noise to apply in units of Jy/Beam (-1 to escape): ");
	    if (cmdbuffer[0] != 0) {
	      add_history(cmdbuffer);
	    }

	    chomp(cmdbuffer);

	    flaterror = strtod(cmdbuffer, &endptr);

	    if ( (flaterror == 0) && (cmdbuffer == endptr) ) {
	      printf("\nInvalid input, please try again...\n\n");
	      continue;
	    }

	    if (flaterror < 0.0) {
	      printf("Escaping command...\n");
	      escape = 1;
	      break;
	    }
	    else {
	      validentry = 1;
	    }
	  }
      
	  if (escape == 1) {
	    break;
	  }

	}

	// Setup the loops
	if (currentset == 888) {
	  looplimit = numdatasets;
	}
	else {
	  looplimit = 1;
	}


	for (i=0; i<looplimit; i++){

	  if (looplimit > 1) {
	    currentset = i;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = i;
	  }


	  if (doallmaps == 1) {
	    looplimit2 = imgnum[currentset];
	  }
	  else {
	    looplimit2 = 1;
	  }


	  for (j=0; j<looplimit2; j++){

	    if (looplimit2 > 1) {
	      mapselect  = j;
	    }
	    else if ( (looplimit2 == 1) && ( (mapselect == 888) || (mapselect == 777) ) ) {
	      mapselect  = j;
	    }

	    // Convert to Jy/Px and update the RMS noise values

	    // Do the changes if a range is requested
	    if (changearange == 1) {

	      if ( (gtlt == 0) && (frequencyobs[currentset][fluxorder[currentset][mapselect]] > errorchange) ) {

		rms[currentset][fluxorder[currentset][mapselect]] = flaterror/beamarea[currentset];
	      }
	      else if ( (gtlt == 1) && (frequencyobs[currentset][fluxorder[currentset][mapselect]] < errorchange) ) {

		rms[currentset][fluxorder[currentset][mapselect]] = flaterror/beamarea[currentset];
	      }
	    }

	    else {
	      rms[currentset][mapselect] = flaterror/beamarea[currentset];
	    }

	    //printf("Flux calibration error for dataset %d at %.2e Hz now set to %.2f%%\n", currentset, frequency[currentset][mapselect], flaterror*100);

	  }

	}

	printf("\nRMS noise now set. Please re-run the appropriate region selection command to update the region errors.\n\n");

      }
      else {
	fprintf(stderr,"\n Error: No datasets have yet been loaded!\n\n");
      }
    
      break;


    case 133: // Switch between exporting maps as fits and standard images

      if (suppresscdf == 1) {
	suppresscdf = 0;
	printf("Chi-squared confidence tables will now be shown.\n");
      }
      else if (suppresscdf == 0) {
	suppresscdf = 1;
	printf("The automatic output of chi-squared confidence tables will now be suppressed.\n");
      }
      else {
	printf("\nThe value of suppresscdf is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", suppresscdf);
      }

      break;



    case 134: // Plot the best fitting injection index by region

      if (numdatasets > 0) {

	if (export == 1) {

	  // Open the directory that is currently set and check it exists
	  dir = opendir(imageloc);

	  if (dir != NULL) {
	    closedir (dir); // If we find it move on without comment
	  }
	  else {
	    printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	    break;
	  }
	}

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to map the best fitting injection index by region (-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (mininjectset[currentset][1] < 1) && (mininjectset[currentset][2] < 1) && (mininjectset[currentset][3] < 1) && (mininjectset[currentset][4] < 1) ) {
	    printf("\nfindinject has not yet been run for this dataset. Please try again...\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	// Ask which model to map and check its validity

	validentry = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Which model would you like to map? (-1 to escape, 1 for JP model, 2 for KP model, 3 Tribble (JP), 4 Tribble (KP)): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  model = strtol(cmdbuffer, &endptr, 10);

	  if ( (model == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (model < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (model > NUMBEROFMODELS) ||  (model == 0) ) { // Check the set exists
	    printf("\nInvalid model number. Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (model == 1) && (mininjectset[currentset][1] < 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: JP injection index values have not yet been calculated for this data set! Please run fitjpmodel first.\n\n");
	    continue;
	  }

	  else if ( (currentset != 888) && (model == 2) && (mininjectset[currentset][2] < 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: KP injection index values have not yet been calculated for this data set! Please run fitkpmodel first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (model == 3) && (mininjectset[currentset][3] < 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Tribble (JP) injection index values have not yet been calculated for this data set! Please run fitjptribble first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (model == 4) && (mininjectset[currentset][4] < 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Tribble (KP) injection index values have not yet been calculated for this data set! Please run fitkptribble first.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}
	
	if ( (contours == 1) && (currentset != 888) ) {

	  if ( (exportfits != 1)  || (export != 1) ) {

	    validentry = 0;

	    printf("============================================\n\n");
	    printf(" Index  |  Frequency (Obs.)  |  Frequency (Rest)  \n\n");

	    // Loop through each of the sets
	    for (a=0; a<imgnum[currentset]; a++) {

	      printf("   %d       %.2e Hz             %.2e Hz\n", fluxorder[currentset][a], frequencyobs[currentset][fluxorder[currentset][a]], frequency[currentset][fluxorder[currentset][a]]);

	    }

	    printf("\n============================================\n");
    

    
	    while (validentry == 0) {

	      cmdbuffer = readline("Select a map frequency for the contours (-1 to escape): ");
	      if (cmdbuffer[0] != 0) {
		add_history(cmdbuffer);
	      }

	      chomp(cmdbuffer);

	      currentmap = strtol(cmdbuffer, &endptr, 10);

	      if ( (currentmap == 0) && (cmdbuffer == endptr) ) {
		printf("\nInvalid input, please try again...\n\n");
		continue;
	      }

	      if (currentmap < 0) {
		printf("Escaping command...\n");
		escape = 1;
		break;
	      }   
	      else if (currentmap >= imgnum[currentset]) { // Check the set exists
		printf("Invalid selection, please try again...\n");
		continue;
	      }
	      else {
		validentry = 1;
	      }
	    }

	  }
	}
	
	else if ( (contours == 1) && (currentset == 888) ) { // Multimaps always uses map 0 as reference. This can easily be changed if needed
	  currentmap = 0;
	}
	else { //If we are not using contours, make sure its still a sensible value
	  currentmap = 0;
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if ( (currentset == 888) && (export == 0) ) {
	  looplimit = numdatasets;
	  extwincontrol = 1;
	  sprintf(output,"/xs");
	  cpgopen(output);
	}
	else if ( (currentset == 888) && (export == 1) ) {
	  looplimit = numdatasets;
	  extwincontrol = 0;
	}
	else {
	  looplimit = 1;
	  extwincontrol = 0;
	}
      

	for (j=0; j<looplimit; j++) {

	  if (looplimit > 1) {
	    currentset = j;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = j;
	  }

	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets
	  if (regionsset[currentset] != 1) {
	    printf("\nThe regions for data set %d have not been set! Skipping this set...\n", j);
	    continue;	  
	  }

	  if ( (model == 1) && (mininjectset[currentset][1] < 1) ) {
	    printf("\nThe JP injection index values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  else if ( (model == 2) && (mininjectset[currentset][2] < 1) ) {
	    printf("\nThe KP injection index values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  else if ( (model == 3) && (mininjectset[currentset][3] < 1) ) {
	    printf("\nThe Tribble (JP) injection index values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  else if ( (model == 4) && (mininjectset[currentset][4] < 1) ) {
	    printf("\nThe Tribble (KP) injection index values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  time( &currenttime );
	  time_struct = localtime ( &currenttime );
	    
	  strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	  
	  if ( (export == 1) && (exportfits != 1) ) {

	    if (model == 1) {
	      sprintf(modelbuff,"JP");
	    }
	    else if (model == 2) {
	      sprintf(modelbuff,"KP");
	    }
	    else if (model == 3) {
	      sprintf(modelbuff,"TribbleJP");
	    }
	    else if (model == 4) {
	      sprintf(modelbuff,"TribbleKP");
	    }
	    else {
	      sprintf(modelbuff,"Unknown");
	    }


	    if (contours == 1) {
	      sprintf(output,"%s/%s_InjectByRegion_Contours_%s_%s.%s/%s", imageloc, settarget[currentset], modelbuff, timebuff, imagetype, imagetype);
	      }
	    else {
	      sprintf(output,"%s/%s_InjectByRegion_%s_%s.%s/%s", imageloc, settarget[currentset], modelbuff, timebuff, imagetype, imagetype);
	    }

	    printf("\nExporting as %s\n", output);
	  }

	  else if ( ( export == 1 ) && (exportfits == 1) ) {
	    if (model == 1) {
	      sprintf(output,"%s/%s_InjectByRegion_JP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else if (model == 2) {
	      sprintf(output,"%s/%s_InjectByRegion_KP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else if (model == 3) {
	      sprintf(output,"%s/%s_InjectByRegion_TribbleJP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else if (model == 4) {
	      sprintf(output,"%s/%s_InjectByRegion_TribbleKP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else {
	      sprintf(output,"%s/%s_InjectByRegion_UnknownModel_%s.fits", imageloc, settarget[currentset], timebuff);
	    }

	    printf("\nExporting as %s\n", output);
	  }
	  else {
	    sprintf(output,"/xs");
	  }
	
	  if (model <= 4) {
	    injectbyregion(inj_regbestinject[currentset][model], regionarray[currentset], flux[currentset][currentmap], xdim[currentset], ydim[currentset], settarget[currentset], border[currentset], zoom, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, firstcontour, contourlogs, output, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, xshift, yshift, export, exportfits, largetxt);
	  }
	  else {
	    fprintf(stderr,"\nError: Model number is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");
	    break;
	  }

	}

	if (extwincontrol == 1) {
	  cpgclos();
	  extwincontrol = 0;
	}

      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;



    case 135: // Plot the best fitting injection index by region

      if (numdatasets > 0) {

	if (export == 1) {

	  // Open the directory that is currently set and check it exists
	  dir = opendir(imageloc);

	  if (dir != NULL) {
	    closedir (dir); // If we find it move on without comment
	  }
	  else {
	    printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	    break;
	  }
	}

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

	// Loop through each of the sets
	for (i=0;i<numdatasets;i++) {

	  list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

	}

	printf("\n========================================================================\n");
    
    
	// Get the reference data set
	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the data set for which to map the chi-squared values of the best fitting injection index by region (-1 to escape, 888 for all): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (regionsset[currentset] != 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Regions have not yet been set for this data set! Please run setregions first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (mininjectset[currentset][1] < 1) && (mininjectset[currentset][2] < 1) && (mininjectset[currentset][3] < 1) && (mininjectset[currentset][4] < 1) ) {
	    printf("\nfindinject has not yet been run for this dataset. Please try again...\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	// Ask which model to map and check its validity

	validentry = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Which model would you like to map? (-1 to escape, 1 for JP model, 2 for KP model, 3 Tribble (JP), 4 Tribble (KP)): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  model = strtol(cmdbuffer, &endptr, 10);

	  if ( (model == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (model < 0) {
	    printf("Escaping command...\n");
	    escape = 1;
	    break;
	  }
	  else if ( (model > NUMBEROFMODELS) ||  (model == 0) ) { // Check the set exists
	    printf("\nInvalid model number. Please try again...\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (model == 1) && (mininjectset[currentset][1] < 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: JP injection index values have not yet been calculated for this data set! Please run fitjpmodel first.\n\n");
	    continue;
	  }

	  else if ( (currentset != 888) && (model == 2) && (mininjectset[currentset][2] < 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: KP injection index values have not yet been calculated for this data set! Please run fitkpmodel first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (model == 3) && (mininjectset[currentset][3] < 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Tribble (JP) injection index values have not yet been calculated for this data set! Please run fitjptribble first.\n\n");
	    continue;
	  }
	  else if ( (currentset != 888) && (model == 4) && (mininjectset[currentset][4] < 1) ) { // Check the regions have been set
	    fprintf(stderr,"\nError: Tribble (KP) injection index values have not yet been calculated for this data set! Please run fitkptribble first.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}
	
	if ( (contours == 1) && (currentset != 888) ) {

	  if ( (exportfits != 1)  || (export != 1) ) {

	    validentry = 0;

	    printf("============================================\n\n");
	    printf(" Index  |  Frequency (Obs.)  |  Frequency (Rest)  \n\n");

	    // Loop through each of the sets
	    for (a=0; a<imgnum[currentset]; a++) {

	      printf("   %d       %.2e Hz             %.2e Hz\n", fluxorder[currentset][a], frequencyobs[currentset][fluxorder[currentset][a]], frequency[currentset][fluxorder[currentset][a]]);

	    }

	    printf("\n============================================\n");
    

    
	    while (validentry == 0) {

	      cmdbuffer = readline("Select a map frequency for the contours (-1 to escape): ");
	      if (cmdbuffer[0] != 0) {
		add_history(cmdbuffer);
	      }

	      chomp(cmdbuffer);

	      currentmap = strtol(cmdbuffer, &endptr, 10);

	      if ( (currentmap == 0) && (cmdbuffer == endptr) ) {
		printf("\nInvalid input, please try again...\n\n");
		continue;
	      }

	      if (currentmap < 0) {
		printf("Escaping command...\n");
		escape = 1;
		break;
	      }   
	      else if (currentmap >= imgnum[currentset]) { // Check the set exists
		printf("Invalid selection, please try again...\n");
		continue;
	      }
	      else {
		validentry = 1;
	      }
	    }

	  }
	}
	
	else if ( (contours == 1) && (currentset == 888) ) { // Multimaps always uses map 0 as reference. This can easily be changed if needed
	  currentmap = 0;
	}
	else { //If we are not using contours, make sure its still a sensible value
	  currentmap = 0;
	}

	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if ( (currentset == 888) && (export == 0) ) {
	  looplimit = numdatasets;
	  extwincontrol = 1;
	  sprintf(output,"/xs");
	  cpgopen(output);
	}
	else if ( (currentset == 888) && (export == 1) ) {
	  looplimit = numdatasets;
	  extwincontrol = 0;
	}
	else {
	  looplimit = 1;
	  extwincontrol = 0;
	}
      

	for (j=0; j<looplimit; j++) {

	  if (looplimit > 1) {
	    currentset = j;
	  }
	  else if ( (looplimit == 1) && (currentset == 888) ) {
	    currentset = j;
	  }

	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets
	  if (regionsset[currentset] != 1) {
	    printf("\nThe regions for data set %d have not been set! Skipping this set...\n", j);
	    continue;	  
	  }

	  if ( (model == 1) && (mininjectset[currentset][1] < 1) ) {
	    printf("\nThe JP injection index values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  else if ( (model == 2) && (mininjectset[currentset][2] < 1) ) {
	    printf("\nThe KP injection index values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  else if ( (model == 3) && (mininjectset[currentset][3] < 1) ) {
	    printf("\nThe Tribble (JP) injection index values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  else if ( (model == 4) && (mininjectset[currentset][4] < 1) ) {
	    printf("\nThe Tribble (KP) injection index values for data set %d have not been calculated! Skipping this set...\n", j);
	    continue;	  
	  }

	  time( &currenttime );
	  time_struct = localtime ( &currenttime );
	    
	  strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	  if ( (export == 1) && (exportfits != 1) ) {

	    if (model == 1) {
	      sprintf(modelbuff,"JP");
	    }
	    else if (model == 2) {
	      sprintf(modelbuff,"KP");
	    }
	    else if (model == 3) {
	      sprintf(modelbuff,"TribbleJP");
	    }
	    else if (model == 4) {
	      sprintf(modelbuff,"TribbleKP");
	    }
	    else {
	      sprintf(modelbuff,"Unknown");
	    }

	    if (contours == 1) {
	      sprintf(output,"%s/%s_InjectChiSquaredByRegion_Contours_%s_%s.%s/%s", imageloc, settarget[currentset], modelbuff, timebuff, imagetype, imagetype);
	    }
	    else {
	      sprintf(output,"%s/%s_InjectChiSquaredByRegion_%s_%s.%s/%s", imageloc, settarget[currentset], modelbuff, timebuff, imagetype, imagetype);
	    }

	    
	    printf("\nExporting as %s\n", output);
	  }

	  else if ( ( export == 1 ) && (exportfits == 1) ) {
	    if (model == 1) {
	      sprintf(output,"%s/%s_InjectChiSquaredByRegion_JP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else if (model == 2) {
	      sprintf(output,"%s/%s_InjectChiSquaredByRegion_KP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else if (model == 3) {
	      sprintf(output,"%s/%s_InjectChiSquaredByRegion_TribbleJP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else if (model == 4) {
	      sprintf(output,"%s/%s_InjectChiSquaredByRegion_TribbleKP_%s.fits", imageloc, settarget[currentset], timebuff);
	    }
	    else {
	      sprintf(output,"%s/%s_InjectChiSquaredByRegion_UnknownModel_%s.fits", imageloc, settarget[currentset], timebuff);
	    }

	    printf("\nExporting as %s\n", output);
	  }
	  else {
	    sprintf(output,"/xs");
	  }
	
	  dof = (imgnum[currentset] - 2);

	  if (model <= 4) {
	    injectchibyregion(inj_regchisquared[currentset][model], regionarray[currentset], flux[currentset][currentmap], xdim[currentset], ydim[currentset], settarget[currentset], border[currentset], zoom, titles, labels, axismarks, extwincontrol, smoothmaps, contours, contourlevels, contourcolour, contourlinestyle, firstcontour, contourlogs, output, model, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, xshift, yshift, export, exportfits, usereduced, dof, largetxt);
	  }
	  else {
	    fprintf(stderr,"\nError: Model number is unknown please try restarting BRATS. If this problem persists please contact Jeremy.Harwood@physics.org\nExiting...\n\n");
	    break;
	  }

	}

	if (extwincontrol == 1) {
	  cpgclos();
	  extwincontrol = 0;
	}

      }
      else {
	fprintf(stderr," Error: No datasets have yet been loaded!\n");
      }

      break;


  case 136: // Extend model data beyond the observed frequnecies


      validentry = 0;
      validentry2 = 0;
      tmp_extrapolatemodel = extrapolatemodel;

      if (extrapolatemodel == 1) {
	printf("Extending of model data is currently set to ON\n");
      }
      else {
	printf("Extending of model data is currently set to OFF\n");
      }

      while (validentry == 0) {

      if (extrapolatemodel == 1) {
	cmdbuffer = readline("Would you like the extending of model data to remain on? (No [0], Yes [1], -1 to escape): ");
      }
      else {
	cmdbuffer = readline("Would you like to turn on the extending of model data? (No [0], Yes [1], -1 to escape): ");
      }

	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	extrapolatemodel = strtol(cmdbuffer, &endptr, 10);

	if ( (extrapolatemodel == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}	
	else {
	  if (extrapolatemodel < 0) {
	    printf("Escaping command...\n");
	    extrapolatemodel = tmp_extrapolatemodel;
	    validentry2 = 1;
	    break;
	  }
	  else if (extrapolatemodel > 1 ) {
	  printf("\nEntry must be either 0 (off) or 1 (on), please try again...\n\n");
	  continue;
	  }
	  else {
	    validentry = 1;
	  }
	}
      }

      if (validentry2 == 1) {
	break;
      }

      if (extrapolatemodel == 0) {
	printf("\nExtending of model data is now set to OFF\n\n");
	break;
      }


      validentry = 0;
      tmp_extrapolationfrequency = extrapolationfrequency;

      printf("Models are currently extended to %.3e Hz (DEFAULT: %.3e Hz)\n", extrapolationfrequency, DEFAULTEXTRAPOLATIONFREQUENCY);

      while (validentry == 0) {

	cmdbuffer = readline("Enter the frequency which the model should be extended to in Hz (float > 0, -1 to escape): ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	extrapolationfrequency= strtod(cmdbuffer, &endptr);

	if ( (extrapolationfrequency == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (extrapolationfrequency <= 0) {
	    printf("Escaping command...\n");
	    extrapolatemodel = tmp_extrapolatemodel; // Also revert back to the origional on/off status
	    extrapolationfrequency = tmp_extrapolationfrequency;
	    break;
	  }
	  else {
	    printf("\nModel data will now be extended to %.3e Hz\n\n", extrapolationfrequency);
	    validentry = 1;
	  }
	}
      }

      break;

      
    case 137: // Switch between standard and large text formats

      if (largetxt == 1) {
	largetxt = 0;
	printf("Text in map and plots will now use the standard size format\n");
      }
      else if (largetxt == 0) {
	largetxt = 1;
	printf("Text in maps and plots will now use the large size format\n");
      }
      else {
	printf("\nThe value of export is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", largetxt);
      }

      break;

            
    case 138: // Plot spectral index error maps

      if (numdatasets > 0) {

	if (export == 1) {

	  // Open the directory that is currently set and check it exists
	  dir = opendir(imageloc);

	  if (dir != NULL) {
	    closedir (dir); // If we find it move on without comment
	  }
	  else {
	    printf("\n*** Error: Cannot open specified directory to output images (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in imageloc ***\n\n", imageloc);
	    break;
	  }
	
	}

	// List the available datasets

	printf("========================================================================\n\n");
	printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");
	// Check at least 1 dataset is loaded
	if (numdatasets > 0) {

	  // Loop through each of the sets
	  for (i=0;i<numdatasets;i++) {

	    list(i, imgnum[i], setname[i], setreg[i], setbg[i]);
	  }
	}
	else {
	  fprintf(stderr," Error: No datasets have yet been loaded!\n");
	}

	printf("\n========================================================================\n");


	tmp_currentset = currentset;

	validentry = 0;
	escape = 0;

	while (validentry == 0) {

	  cmdbuffer = readline("Enter the number of the data set to determine the spectral index map (-1 to escape): ");
	  if (cmdbuffer[0] != 0) {
	    add_history(cmdbuffer);
	  }

	  chomp(cmdbuffer);

	  currentset = strtol(cmdbuffer, &endptr, 10);

	  if ( (currentset == 0) && (cmdbuffer == endptr) ) {
	    printf("\nInvalid input, please try again...\n\n");
	    continue;
	  }

	  if (currentset < 0) {
	    printf("Escaping command...\n");
	    currentset = tmp_currentset;
	    escape = 1;
	    break;
	  }
	  else if ( (currentset >= numdatasets) && (currentset != 888) ) { // Check the set exists
	    printf("\nThis value is greater than the number of data sets that exists! Please try again...\n\n");
	    continue;
	  }
	  else if ( (specindexset[currentset] != 1) && (currentset != 888) ) { // Check the spectral index has been set
	    fprintf(stderr,"\nError: Spectral indices have not yet been fitted to this data set! Please run the specindex command first.\n\n");
	    continue;
	  }
	  else if ( (specindextype[currentset] != 2) && (currentset != 888) ) { // Check the spectral index has been set
	    fprintf(stderr,"\nError: Invalid fitting method used. Weighted (GSL) least squares fitting must be used to enable this function.\n\n");
	    continue;
	  }
	  else {
	    validentry = 1;
	  }
	}
      
	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if (contours == 1) {

	  if ( (exportfits != 1)  || (export != 1) ) {

	    validentry = 0;

	    printf("============================================\n\n");
	    printf(" Index  |  Frequency (Obs.)  |  Frequency (Rest)  \n\n");

	    // Loop through each of the sets
	    for (a=0; a<imgnum[currentset]; a++) {

	      printf("   %d       %.2e Hz             %.2e Hz\n", fluxorder[currentset][a], frequencyobs[currentset][fluxorder[currentset][a]], frequency[currentset][fluxorder[currentset][a]]);

	    }

	    printf("\n============================================\n");
    

	    while (validentry == 0) {

	      cmdbuffer = readline("Select a maps frequency for the contours (-1 to escape): ");
	      if (cmdbuffer[0] != 0) {
		add_history(cmdbuffer);
	      }

	      chomp(cmdbuffer);

	      currentmap = strtol(cmdbuffer, &endptr, 10);

	      if ( (currentmap == 0) && (cmdbuffer == endptr) ) {
		printf("\nInvalid input, please try again...\n\n");
		continue;
	      }

	      if (currentmap < 0) {
		printf("Escaping command...\n");
		escape = 1;
		break;
	      }
	
	      else if (currentmap >= imgnum[currentset]) { // Check the set exists
		printf("\nInvalid selection, please choose again...\n\n");
		continue;
	      }
	      else {
		validentry = 1;
	      }

	    }
	  }
	}
	
	else if ( (contours == 1) && (currentset == 888) ) { // Multimaps always uses map 0 as reference. This can easily be changed if needed
	  currentmap = 0;
	}
	else {
	  currentmap = 0;
	}
      
	// Break the second loop if they want to escape
	if (escape == 1) {
	  break;
	}

	if ( (currentset == 888) && (export == 0) ) {
	  looplimit = numdatasets;
	  extwincontrol = 1;
	  sprintf(output,"/xs");
	  cpgopen(output);
	}
	else if ( (currentset == 888) && (export == 1) ) {
	  looplimit = numdatasets;
	  extwincontrol = 0;
	}
	else {
	  looplimit = 1;
	  extwincontrol = 0;
	}
	
	for (j=0; j<looplimit; j++){

	  if ( (looplimit > 1) || ( (looplimit == 1) && (currentset == 888) ) ) {
	    currentset = j;
	  }
	  
	  // This will only execute in the case of 888, a previous (similar) statement catches single data sets
	  if (specindexset[currentset] != 1) {
	    printf("\n*** Warning: The spectral indices for data set %d have not been set! Skipping this set... ***\n\n", j);
	    continue;	  
	  }
	  else if (specindextype[currentset] != 2) { // Check the set exists
	    fprintf(stderr,"\n*** Warning: Invalid fitting method used for dataset %d (weighted (GSL) least squares fitting only). Skipping this set... ***\n\n", j);
	    continue;
	  }

	  if (export == 1) {

	    extwincontrol = 1; // We always want external control if exporting

	    time( &currenttime );
	    time_struct = localtime ( &currenttime );

	    strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

	    if ( ( export == 1 ) && (exportfits != 1) ) {

	      if (contours == 1) {
		sprintf(output,"%s/%s_SpectralIndexErrors_Contours_%s.%s/%s", imageloc, settarget[currentset], timebuff, imagetype, imagetype);
	      }
	      else {
		sprintf(output,"%s/%s_SpectralIndexErrors_%s.%s/%s", imageloc, settarget[currentset], timebuff, imagetype, imagetype);
	      }
	      printf("\nExporting as %s\n", output);
	      cpgopen(output);
	    }
	    else if ( ( export == 1 ) && (exportfits == 1) ) {
	      sprintf(output,"%s/%s_SpectralIndexErrors_%s.fits", imageloc, settarget[currentset], timebuff);
	      printf("\nExporting as %s\n", output);
	    }
	  }
	  else {
	    sprintf(output,"/xs");
	  }

	  sprintf(maptitle,"     Spectral Index Error Map of %s (%.2e to %.2e Hz)", settarget[currentset],frequencyobs[currentset][fluxorder[currentset][0]], frequencyobs[currentset][fluxorder[currentset][imgnum[currentset]-1]]);

	  specindexerrormap(specindexerror, regnumber[currentset], xdim, ydim, settarget[currentset], minsi, maxsi, regionarray, currentset, maptitle, zoom, border[currentset], titles, labels, axismarks, symbol, smoothmaps, flux[currentset][currentmap], extwincontrol, contours, contourlevels, contourcolour, contourlinestyle, firstcontour, contourlogs, output, ra[currentset][posmap], dec[currentset][posmap], delt1[currentset][posmap], delt2[currentset][posmap], crpix1[currentset][posmap], crpix2[currentset][posmap], crota1[currentset][posmap], crota2[currentset][posmap], eq[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, wedgemultiplier, xshift, yshift, frequency[currentset][fluxorder[currentset][0]], frequency[currentset][fluxorder[currentset][imgnum[currentset]-1]], export, exportfits, largetxt);

	  if ( (export == 1) && (exportfits != 1) && (extwincontrol == 1) ){
	    cpgclos();
	  }
	}

	if (extwincontrol == 1) {
	  extwincontrol = 0;
	}
    
      }
      else {
	fprintf(stderr," Error: No data sets have yet been loaded!\n");
      }
      
      break;


    case 139: // Select the colour scheme for mapping

      // List the available choices
      printf("======================================================\n");
      printf("  Selection  |  Data type                             \n");
      printf("------------------------------------------------------\n");
      printf("      0      |  Rainbow (default)                     \n");
      printf("      1      |  Heat                                  \n");
      printf("      2      |  Viridis                               \n");
      printf("      3      |  Red                                   \n");
      printf("      4      |  Green                                 \n");
      printf("      5      |  Blue                                  \n");
      printf("      6      |  Grey scale                            \n");
      printf("      -1     |  Escape command                        \n");
      printf("======================================================\n");


      validentry = 0;
      tmp_pltcolour = pltcolour;

      printf("Colour scheme currently set to %d\n", pltcolour);

      while (validentry == 0) {

	cmdbuffer = readline("Select a colour scheme to use when mapping: ");
	if (cmdbuffer[0] != 0) {
	  add_history(cmdbuffer);
	}

	chomp(cmdbuffer);

	pltcolour = strtol(cmdbuffer, &endptr, 10);

	if ( (pltcolour == 0) && (cmdbuffer == endptr) ) {
	  printf("\nInvalid input, please try again...\n\n");
	  continue;
	}
	else {
	  if (pltcolour < 0) {
	    printf("Escaping command...\n");
	    pltcolour = tmp_pltcolour;
	    break;
	  }
	  else if (pltcolour > 5) {
	    printf("\nInvalid entry, please try again...\n\n");
	    continue;
	  }
	  else {
	    printf("\nThe colour scheme for mapping is now set to %d\n\n", pltcolour);
	    validentry = 1;
	  }
	}
      }

      
      break;

    case 140: // Select if the colour scheme should be inverted or not

      if (invertcolours == 1) {
	invertcolours = 0;
	printf("The colour scheme for mapping is now set to its standard orientation\n");
      }
      else if (invertcolours == 0) {
	invertcolours = 1;
	printf("The colour scheme for mapping is now inverted\n");
      }
      else {
	printf("\nThe value of invertcolours is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", invertcolours);
      }
      
      break;


    case 141: // Set what type of fitting should be done for the spectral index

      if (forceerrortype == 1) {
	forceerrortype = 0;
	printf("Error type forcing is now tuned off. The error type will no longer be forced to standard deviation for spectral index fitting with 2 data point.\n");
      }
      else if (forceerrortype == 0) {
	forceerrortype = 1;
	printf("Error type forcing is now turned on. The error type will now be forced to standard deviation for spectral index fitting with 2 data point.\n");
      }
      else {
	printf("\nThe value of forceerrortype is set to an undefined case of %d! Please try restarting BRATS. If this problem continues, please contact Jeremy.Harwood@physics.org\n", forceerrortype);
      }
      
      break;


    case 999:
      fprintf(stderr,"Unknown command - type help for a list of commands\n");
      break;

    case '?':
      fprintf(stderr,"Unknown command - type help for a list of commands\n");
      break;

    default:
      fprintf(stderr,"Unknown command - type help for a list of commands\n");
      break;
   
    }

  }


  // Cleanup here (This is mainly to keep valgrind happy and debugging clean!)

  // This is is only required for developement
  /*
    if (jpmodelmemoryset == 1) {

    for (a=0; a<MAXDATASETS; a++) {

    // If the memory has already been set for this dataset, free it up so we can make some mew ones!
    if (jpfitmemset[a] == 1) {

    free(jpchisquared[a]);
    free(jpbestage[a]);
    free(jpmodelalpha[a]);
    free(jpbestnorm[a]);
    free(jpmodelflux[a]);
    }
    }

    free(jpchisquared);
    free(jpbestage);
    free(jpmodelalpha);
    free(jpbestnorm);
    free(jpmodelflux);
    }

    if (jptribmodelmemoryset == 1) {

    for (a=0; a<MAXDATASETS; a++) {

    // If the memory has already been set for this dataset, free it up so we can make some mew ones!
    if (jptribfitmemset[a] == 1) {

    free(jptribchisquared[a]);
    free(jptribbestage[a]);
    free(jptribmodelalpha[a]);
    free(jptribbestnorm[a]);
    free(jptribmodelflux[a]);
    }
    }

    free(jptribchisquared);
    free(jptribbestage);
    free(jptribmodelalpha);
    free(jptribbestnorm);
    free(jptribmodelflux);
    }


    if (kpmodelmemoryset == 1) {

    for (a=0; a<MAXDATASETS; a++) {

    // If the memory has already been set for this dataset, free it up so we can make some mew ones!
    if (kpfitmemset[a] == 1) {

    free(kpchisquared[a]);
    free(kpbestage[a]);
    free(kpmodelalpha[a]);
    free(kpbestnorm[a]);
    free(kpmodelflux[a]);
    }
    }

    free(kpchisquared);
    free(kpbestage);
    free(kpmodelalpha);
    free(kpbestnorm);
    free(kpmodelflux);
    }



    if (kptribmodelmemoryset == 1) {

    for (a=0; a<MAXDATASETS; a++) {

    // If the memory has already been set for this dataset, free it up so we can make some mew ones!
    if (kptribfitmemset[a] == 1) {

    free(kptribchisquared[a]);
    free(kptribbestage[a]);
    free(kptribmodelalpha[a]);
    free(kptribbestnorm[a]);
    free(kptribmodelflux[a]);
    }
    }

    free(kptribchisquared);
    free(kptribbestage);
    free(kptribmodelalpha);
    free(kptribbestnorm);
    free(kptribmodelflux);
    }


    if (curvearrayset == 1) {

    for (a=0; a<MAXDATASETS; a++) {

    if (curvatureset[a] == 1) {

    for(i=0;i<=regnumber[a];i++) {

    free(curvearray[a][i]);
    }

    free(curvearray[a]);
    }
    }

    free(curvearray);
    free(mincurve);
    free(maxcurve);
    }


    if (regmemallocated == 1) {
  
    for (a=0; a<MAXDATASETS; a++) {

    for (i=0; i<xdim[a]; i++) {

    free(regionarray[a][i]);
    }

    free(regionarray[a]);
    free(regionlocx[a]);
    free(regionlocy[a]);
    }
  
    for (a=0; a<MAXDATASETS; a++) {
    for (i=0; i<imgnum[a]; i++) {

    free(fluxerror[a][i]);
    free(regflux[a][i]);
    }
    free(fluxerror[a]);
    free(regflux[a]);
    }

    free(regmaxflux);
    free(regminflux);
    free(regflux);
    free(regionarray);
    free(regnumber);
    free(fluxerror);
    free(regionlocx);
    free(regionlocy);
    free(regionsize);
    free(setaveraged);
    }

    if (specallocated == 1) {

    for (i=0; i<MAXDATASETS; i++) {
    free(alpha[i]);
    }

    free(alpha);

    free(minsi);
    free(maxsi);
    }

    if (arraysallocated >= 1) {

    for (a=0; a<MAXDATASETS; a++) {
    for (i=0; i<imgnum[a]; i++) {
    for (j=0; j<xdim[a]; j++) {

    free(flux[a][i][j]);
    }

    free(flux[a][i]);
    }
    }     

    for (a=0; a<MAXDATASETS; a++) {
    free(maxflux[a]);
    free(minflux[a]);  
    free(frequency[a]);
    free(rms[a]);
    free(fluxorder[a]);
    free(flux[a]);
    free(ra[a]);
    free(dec[a]);
    }

    free(xdim);
    free(ydim);
    free(bmaj);
    free(bmin);
    free(bpa);
    free(beamarea);
    free(imgnum);
    free(rms);
    free(flux);
    free(fluxorder);
    free(frequency);
    free(ra);
    free(dec);
    free(maxflux);
    free(minflux);
    }

    for (a=0; a<numdatasets; a++) {
    free(settarget[a]);
    }

    // Free the persistant variable memory
    free(regionsset);
    free(settarget);
    free(curvatureset);
    free(specindexset);
    free(jpfitmemset);
    free(kpfitmemset);
    free(jptribfitmemset);
    free(kptribfitmemset);
    free(jpinject);
    free(kpinject);
    free(bgbuff);

  */

  printf("Exiting...\n");

  return 0;
}




