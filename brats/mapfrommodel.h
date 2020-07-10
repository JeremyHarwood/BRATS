/*
  Create a fits image based on a given model and dataset.

  Requires the file "chomp.h" and <funtools.h> plus the standard headers to have been previously declared

*/

//Definitions
#ifndef MAXCMDLENGTH
#define MAXCMDLENGTH 1024
#endif

#ifndef MAXNUMFILES
#define MAXNUMFILES 50
#endif

#ifndef GFACTOR
#define GFACTOR (2.0*sqrt(2.0*log(2.0)))
#endif


// Declare the function prototypes

int syscom(const char[]);

int findmodelvalues(double usr_gmin, double usr_gmax, double fieldstrength, int model, float *reconmapflux, int xdim, int ydim, float usr_frequency, float *bestnorm, float *bestage, int **regionarray, float inject, float beamarea, float redshift);

int mapfrommodel(char *imageloc, char *dirname, int xdim, int ydim, char *maptarget, double usr_gmin, double usr_gmax, double fieldstrength, int model, float *bestnorm, float *bestage, int **regionarray, float inject, float beamarea, float redshift) { 


  // Declare the variables for this task

  char **fitsfiles, modelname[32], *endptr, *cmdbuffer, datebuff[32], timebuff[32], historytext[MAXCMDLENGTH], *units;

  double ra, dec, delt1, delt2, datamax, datamin, eq, crpix1, crota1, crpix2, crota2, crpix4, crota4, crval4, delt4;

  int required, i, j, got;
  float *reconmapflux;
  char *fileloc, *reconnloc, imagename[256];

  //FITSHead head;
  Fun *fitsconn, *modelconn;
  DIR *direct;

  float usr_frequency = 1e9;

  struct dirent *entry;

  // Set up variables and structs for times and dates
  time_t currenttime;
  struct tm * time_struct;


  // Get the directory name and format it
  required = 1;

  while (required == 1) {

    // Get the scaling factor
    cmdbuffer = readline("Please enter a frequency (in Hz) for which to construct the image (-1 to escape): ");
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    chomp(cmdbuffer);

    usr_frequency = strtod(cmdbuffer, &endptr);

    if ( (usr_frequency == 0) && (cmdbuffer == endptr) ) {
      printf("\nInvalid input, please try again...\n\n");
      continue;
    }
    else if (usr_frequency < 0.0) {
      printf("Escaping command...\n\n");
      return(100);
    }
    else {
      printf("Combined map set to a frequency of %.2e Hz\n", usr_frequency);
      required = 0;
    }
  }


  /* Check all the files and directories exists then put them in the correct format  */

  // Allocating memory

  fitsfiles = (char **)calloc(MAXNUMFILES, sizeof(char *)); 
  
  //Open the directory
  direct = opendir(dirname);

  //Set the map counter to 0
  i = 0;

  //Check the directory exists
  if (direct != NULL) {

    // While there are maps to read...
    while ((entry = readdir(direct))) {
      if (strlen(entry->d_name) > MAXCMDLENGTH+1) {
	fprintf(stderr,"\n*** Error: File path must be less than %d characters long! ***\n\n", MAXCMDLENGTH);
	return 2;
      }
      //Do nothing if we are reading . or ..
      else if ((strcmp(entry->d_name, ".") == 0) || (strcmp(entry->d_name, "..") == 0)) {
	;
      }

      // Get the fits file data...
      else {
	fitsfiles[i] = (char *)calloc((MAXCMDLENGTH+1),sizeof(char));
	strcpy(fitsfiles[i], entry->d_name);

	//Increase the map counter
	i++;
      }
    }

    closedir (direct);
  }

  //If we can't find or access the directory, kick back to the command prompt
  else {
    fprintf(stderr,"*** Error: Unable to locate FITS directory ***\n\n");
    return 404;
  }

  /* Go off and get the data */

  printf("Determining the reconstructed map fluxes...\n");

  reconmapflux = (float *)calloc(xdim*ydim, sizeof(float));

  findmodelvalues(usr_gmin, usr_gmax, fieldstrength, model, reconmapflux, xdim, ydim, usr_frequency, bestnorm, bestage, regionarray, inject, beamarea, redshift);

  // Find the min max values
  datamax = -1e32;
  datamin = 1e32;

  for(i=0; i<xdim; i++) {
    for (j=0; j<ydim; j++){

      if (reconmapflux[i*xdim+j] > datamax) {
	datamax = reconmapflux[i*xdim+j];
      }
      if (reconmapflux[i*xdim+j] < datamin) {
	datamin = reconmapflux[i*xdim+j];
      }
    }
  }

  fileloc = (char *)calloc(MAXCMDLENGTH,sizeof(char));
  reconnloc = (char *)calloc(MAXCMDLENGTH,sizeof(char));

  printf("Loading (real) map parameters...\n");

   // Create the map

  if (model == 1) {
    sprintf(modelname,"JP");
  }
  else if (model == 2) {
    sprintf(modelname,"KP");
  }
  else if (model == 3) {
    sprintf(modelname,"JPTRIB");
  }
  else if (model == 4) {
    sprintf(modelname,"KPTRIB");
  }
  else {
    sprintf(modelname,"Unknown");
  }


  //Get the current date and time
  time( &currenttime );
  time_struct = localtime ( &currenttime );

  strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

  sprintf(imagename,"%s_%.2e_%s_recon_%s.fits", maptarget, usr_frequency, modelname, timebuff);
  sprintf(reconnloc,"./%s/%s", imageloc, imagename);

   // Create the file and write the header (currently AIPS format)
    
    modelconn = FunOpen(reconnloc, "w", NULL);

    // Check we can access the file
    if (modelconn == NULL){
      fprintf(stderr,"\n*** Error: Unable to open %s, please check the file exists and access permissions are correct. ***\n\n", fitsfiles[i]);
      return 404;
    }
    
    // Get the header of the first map loaded
    sprintf(fileloc,"./%s/%s", dirname, fitsfiles[0]);
    fitsconn = FunOpen(fileloc, "r", NULL);

    // Check we can access the file
    if (fitsconn == NULL){
      fprintf(stderr,"\n*** Error: Unable to open %s, please check the file exists and access permissions are correct. ***\n\n", fitsfiles[i]);
      return 404;
    }

    //FunInfoGet(fitsconn, FUN_HEADER, &head, 0);

    // Set this as the header for the new map
    //FunInfoPut(modelconn, FUN_HEADER, &head, 0);

    //Reformat the datetime to that required for the header
    strftime(datebuff,32,"%Y-%m-%d", time_struct);

    //Retrieve the remaining header info needed
    eq = FunParamGetd(fitsconn, "EQUINOX", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain the EQUINOX of %s! ***\n\n", fitsfiles[0] );
      return 4;
    }
    ra = FunParamGetd(fitsconn, "OBSRA", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain the RA of %s! ***\n\n", fitsfiles[0] );
      return 4;
    }
    dec = FunParamGetd(fitsconn, "OBSDEC", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain the Dec of %s! ***\n\n", fitsfiles[0] );
      return 4;  
    }

    crpix1 = FunParamGetd(fitsconn, "CRPIX1", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain CRPIX1 for %s! ***\n\n", fitsfiles[0] );
      return 4;  
    }
    crota1 = FunParamGetd(fitsconn, "CROTA1", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain CROTA1 for %s! ***\n\n", fitsfiles[0] );
      return 4;  
    }
    crpix2 = FunParamGetd(fitsconn, "CRPIX2", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain CRPIX2 for %s! ***\n\n", fitsfiles[0] );
      return 4;  
    }
    crota2 = FunParamGetd(fitsconn, "CROTA2", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain CROTA2 for %s! ***\n\n", fitsfiles[0] );
      return 4;  
    }
    crval4 = FunParamGetd(fitsconn, "CRVAL4", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain CRVAL4 for %s! ***\n\n", fitsfiles[0] );
      return 4;  
    }
    crpix4 = FunParamGetd(fitsconn, "CRPIX4", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain CRPIX4 for %s! ***\n\n", fitsfiles[0] );
      return 4;  
    }
    crota4 = FunParamGetd(fitsconn, "CROTA4", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain CROTA4 for %s! ***\n\n", fitsfiles[0] );
      return 4;  
    }
    delt4 = FunParamGetd(fitsconn, "CDELT4", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain CDELT4 for %s! ***\n\n", fitsfiles[0] );
      return 4;  
    }
    units = FunParamGets(fitsconn, "BUNIT", 0, NULL, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain CDELT4 for %s! ***\n\n", fitsfiles[0] );
      return 4;  
    }
    else if (strcasecmp("JY/BEAM", units)) {
      fprintf(stderr,"\n*** Error: The units used in the map %s (as specified in the header) are not JY/BEAM! ***\n\n", fitsfiles[i] );
      return 4;
    }
    delt1=FunParamGetd(fitsconn, "CDELT1", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to determine the pixel size! ***\n\n");
      return 4;
    }
    delt2=FunParamGetd(fitsconn, "CDELT2", 0, 0.0, &got);
    if (got) {   
      if (delt2!=-delt1) {
	fprintf(stderr,"\n*** Error: Pixels aren't square! ***\n\n");
	return 4;
      }
    } else {
      fprintf(stderr,"\n*** Error: Unable to determine the second pixel size! ***\n\n");
      return 4;
    }
      
    //delt1*=3600.0;

    // Write the header
    FunParamPutb(modelconn, "SIMPLE", 0, 1, NULL, 1);
    FunParamPuti(modelconn, "BITPIX", 0, -32, NULL, 1);
    FunParamPuti(modelconn, "NAXIS", 0, 4, "", 1);
    FunParamPuti(modelconn, "NAXIS1", 0, xdim, "X dimension (px)", 1);
    FunParamPuti(modelconn, "NAXIS2", 0, ydim, "Y dimension (px)", 1);
    FunParamPuti(modelconn, "NAXIS3", 0, 1, NULL, 1);
    FunParamPuti(modelconn, "NAXIS4", 0, 1, NULL, 1);
    FunParamPutb(modelconn, "EXTEND", 0, 1, "Tables following main image", 1);
    FunParamPutb(modelconn, "BLOCKED", 0, 1, "Tape may be blocked", 1);
    FunParamPuts(modelconn, "OBJECT", 0,  maptarget, "Source name", 1);
    FunParamPuts(modelconn, "TELESCOPE", 0, "BRATS", "Broadband Radio Astronomy Tools", 1);
    FunParamPuts(modelconn, "INSTRUME", 0, "", NULL, 1);
    FunParamPuts(modelconn, "OBSERVER", 0, "BRATS", "Broadband Radio Astronomy Tools", 1);
    //FunParamPuts(modelconn,"DATE-OBS", 0, "", "Obs start date YYYY-MM-DD", 1);
    FunParamPuts(modelconn, "DATE-MAP", 0, datebuff, "Last processing date YYYY-MM-DD", 1);
    FunParamPutd(modelconn, "BSCALE", 0, 1.0, 11, "REAL = TAPE * BSCALE + BZERO", 1);
    FunParamPutd(modelconn, "BZERO", 0, 0.0, 9, NULL, 1);
    FunParamPuts(modelconn, "BUNIT", 0,  units, "Units of flux", 1);
    FunParamPutd(modelconn, "EQUINOX", 0, eq, 11, "Epoch of RA DEC", 1);
    FunParamPutd(modelconn, "OBSRA", 0, ra, 11, "Antenna pointing RA", 1);
    FunParamPutd(modelconn, "OBSDEC", 0, dec, 11, "Antenna pointing DEC", 1);
    FunParamPutd(modelconn, "DATAMAX", 0, datamax, 11, "Maximum pixel value", 1);
    FunParamPutd(modelconn, "DATAMIN", 0, datamin, 11, "Minimum pixel value", 1);
    FunParamPuts(modelconn, "CTYPE1", 0, "RA---SIN", NULL, 1);
    FunParamPutd(modelconn, "CRVAL1", 0, ra, 11, "RA", 1);
    FunParamPutd(modelconn, "CDELT1", 0, delt1, 11, NULL, 1);
    FunParamPutd(modelconn, "CRPIX1", 0, crpix1, 11, NULL, 1);
    FunParamPutd(modelconn, "CROTA1", 0, crota1, 11, NULL, 1);
    FunParamPuts(modelconn, "CTYPE2", 0, "DEC--SIN", NULL, 1);
    FunParamPutd(modelconn, "CRVAL2", 0, dec, 11, "Dec", 1);
    FunParamPutd(modelconn, "CDELT2", 0, delt2, 11, NULL, 1);
    FunParamPutd(modelconn, "CRPIX2", 0, crpix2, 11, NULL, 1);
    FunParamPutd(modelconn, "CROTA2", 0, crota2, 11, NULL, 1);
    FunParamPuts(modelconn, "CTYPE3", 0, "FREQ", NULL, 1);
    FunParamPutd(modelconn, "CRVAL3", 0, (double)usr_frequency, 11, "Frequency", 1);
    FunParamPutd(modelconn, "CDELT3", 0, 0.0, 11, NULL, 1);
    FunParamPutd(modelconn, "CRPIX3", 0, 0.0, 11, NULL, 1);
    FunParamPutd(modelconn, "CROTA3", 0, 0.0, 11, NULL, 1);
    FunParamPuts(modelconn, "CTYPE4", 0, "STOKES", NULL, 1);
    FunParamPutd(modelconn, "CRVAL4", 0, crval4, 11, "Stokes", 1);
    FunParamPutd(modelconn, "CDELT4", 0, delt4, 11, NULL, 1);
    FunParamPutd(modelconn, "CRPIX4", 0, crpix4, 11, NULL, 1);
    FunParamPutd(modelconn, "CROTA4", 0, crota4, 11, NULL, 1);
    FunParamPuts(modelconn, "ORIGIN", 0, "BRATS", "Broadband Radio Astronomy Tools", 1);


    strcpy(historytext, "--------------------------------------------------------------------");
    FunParamPuts(modelconn, "HISTORY", 0, historytext, NULL, 1);
    strcpy(historytext, "BRATS  Created by the Broadband Radio Astronomy Tools");
    FunParamPuts(modelconn, "HISTORY", 0, historytext, NULL, 1);
    strcpy(historytext, "BRATS  http://www.askanastronomer.co.uk/brats");
    FunParamPuts(modelconn, "HISTORY", 0, historytext, NULL, 1);
    strcpy(historytext, "BRATS  Extrapolated/interpolated image made using the ");
    strcat(historytext, modelname);
    strcat(historytext, " model of spectral ageing.");
    FunParamPuts(modelconn, "HISTORY", 0, historytext, NULL, 1);


    //FunParamPutd(modelconn, "CRVAL3", 0,  (double)usr_frequency, 11, "Map Frequency", 1);

    // Output the image
    FunImagePut(modelconn, reconmapflux, xdim, ydim, -32, NULL);
    FunClose(modelconn);
    FunClose(fitsconn);


    printf("\nReconstructed map has been exported to %s\n\n", reconnloc);


  //Free up the memory
  for (i=0; i<MAXNUMFILES; i++) {
    free(fitsfiles[i]);
  }

  free(fileloc);
  free(fitsfiles);
  free(reconmapflux);

  return 0;

}

