/*
  Loads a series of maps and combine them in the image plane.

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

int diffmap(int casadata, char *imageloc, int ***regionarray, int *regionsset, int *reg_xdim, int *reg_ydim, int numdatasets, int *imgnum, char setname[][MAXCMDLENGTH], char setreg[][MAXCMDLENGTH], char setbg[][MAXCMDLENGTH]) { 


  // Declare the variables for this task

  int required, i, j, a, got, m, *xdim, *ydim, currentset, tmp_casadata;
  float *bmaj, *bmin, *bpa, *beamarea, **imgbuff, *subtractedmap, ra, dec, *frequency;
  double delt1, delt2, datamax, datamin, eq, crpix1, crota1, crpix2, crota2, crpix4, crota4, crval4, delt4;
  char *fileloc, **units, **maptarget, *diffloc, imagename[256], **fitsfiles, *endptr, *cmdbuffer, datebuff[32], timebuff[32], historytext[MAXCMDLENGTH], dbltochar[MAXCMDLENGTH];

  // dirname[MAXCMDLENGTH], 

  FITSHead head;
  FITSCard card;
  Fun *fitsconn, *funcomb;

  // Set up variables and structs for times and dates
  time_t currenttime;
  struct tm * time_struct;

  int validentry = 0;
  int useregions = 0;

  char *str_casadata, str_compare_aips[8], str_compare_casa[8];
  strcpy(str_compare_aips, "AIPS");
  strcpy(str_compare_casa, "CASA");

  fitsfiles = (char **)calloc(2, sizeof(char *)); 
  fitsfiles[0] = (char *)calloc((MAXCMDLENGTH+1),sizeof(char));
  fitsfiles[1] = (char *)calloc((MAXCMDLENGTH+1),sizeof(char));

  //Get the file names and format them

  required = 1;

  while (required == 1) {

    //Get the background file name and format it
    cmdbuffer = readline("Map to be subtracted from: ");
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    chomp(cmdbuffer);

    strcpy(fitsfiles[0], cmdbuffer);

    if ( (fitsfiles[0] == 0) && (cmdbuffer == endptr) ) {
      printf("\nInvalid input, please try again...\n\n");
      continue;
    }

    if (strlen(fitsfiles[0]) >= MAXCMDLENGTH-1) {
      fprintf(stderr,"\n*** Error: File path must be less than %d characters long! ***\n\n", MAXCMDLENGTH);
      continue;
    }
    else if (strlen(fitsfiles[0]) == 0) {
      printf("\nPlease enter a file name.\n\n");
      continue;
    }
    else if ((strstr(fitsfiles[0], "esc") != NULL) && (strlen(fitsfiles[0]) == 3)) {
      printf("Escaping command...\n");
      return 100;
    }
    else if ((strstr(fitsfiles[0], "ls") != NULL) && (strlen(fitsfiles[0]) == 2)) {
      syscom("ls");
      continue;
    }
    else if( access(fitsfiles[0], R_OK ) == -1 ) {
      printf("\nSpecified file either does not exist or has incorrect permissions set (unable to read). Please try again.\n\n");
      continue;
    }
    else {
      required = 0;
    }
  }

  required = 1;

  while (required == 1) {

    //Get the background file name and format it
    cmdbuffer = readline("Map to subtract: ");
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    chomp(cmdbuffer);

    strcpy(fitsfiles[1], cmdbuffer);

    if ( (fitsfiles[1] == 0) && (cmdbuffer == endptr) ) {
      printf("\nInvalid input, please try again...\n\n");
      continue;
    }

    if (strlen(fitsfiles[1]) >= MAXCMDLENGTH-1) {
      fprintf(stderr,"\n*** Error: Folder path must be less than %d characters long! ***\n\n", MAXCMDLENGTH);
      continue;
    }
    else if (strlen(fitsfiles[1]) == 0) {
      printf("\nPlease enter a file name.\n\n");
      continue;
    }
    else if ((strstr(fitsfiles[1], "esc") != NULL) && (strlen(fitsfiles[1]) == 3)) {
      printf("Escaping command...\n");
      return 100;
    }
    else if ((strstr(fitsfiles[1], "ls") != NULL) && (strlen(fitsfiles[1]) == 2)) {
      syscom("ls");
    }
    else if( access(fitsfiles[1], R_OK ) == -1 ) {
      printf("\nSpecified file either does not exist or has incorrect permissions set (unable to read). Please try again.\n\n");
      continue;
    }
    else {
      required = 0;
    }
  }

  validentry = 0;

  while (validentry == 0) {

    cmdbuffer = readline("Would you like to restrict the subtraction to specific regions (0 flux elsewhere)? (No [0], Yes [1], -1 to escape): ");
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    chomp(cmdbuffer);

    useregions = strtol(cmdbuffer, &endptr, 10);

    if ( (useregions == 0) && (cmdbuffer == endptr) ) {
      printf("\nInvalid input, please try again...\n\n");
      continue;
    }
    else if (useregions > 1) {
      printf("\nYou must select either no (0) or yes (1). Please try again.\n\n");
      continue;
    }
    else if (useregions < 0) {
      printf("Escaping command...\n");
      return 100;
    }
    else {
      validentry = 1;
    }
  }

  if (useregions == 1) {

    if (numdatasets > 0) {

      // List the available datasets

      printf("========================================================================\n\n");
      printf(" Index  |    Directory    |    Background    |    Region   \n\n");

      // Loop through each of the sets
      for (i=0;i<numdatasets;i++) {

	list(i, imgnum[i], setname[i], setreg[i], setbg[i]);

      }

      printf("\n========================================================================\n");


      validentry = 0;

      while (validentry == 0) {

	cmdbuffer = readline("Enter the data set which contains the regions to be used (-1 to escape, integers only): ");
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
	  return 100;
	}

	else if (currentset >= numdatasets) { // Check the set exists
	  fprintf(stderr,"\nError: This dataset does not exist! Please try again...\n\n");
	  continue;
	}
	else if (regionsset[currentset] != 1) {
	  printf("\nRegions have not yet been created for this dataset. Please run the setregions command and try again, or choose a different dataset...\n\n");
	  continue;
	}
	else {
	  validentry = 1;
	}
      }

    }

    else {
      useregions = 0;
      fprintf(stderr,"\n*** Error: No datasets have yet been loaded! Proceeding with full map subtraction... *** \n");
    }
  }


  // Put some space between the input and the output
  printf("\n");
  /* Check all the files and directories exists then put them in the correct format  */

  // Allocating memory

  xdim = (int *)calloc(2, sizeof(int));
  ydim = (int *)calloc(2, sizeof(int));
  

  /* Start looping through the maps and getting the data */

  fileloc = (char *)calloc(MAXCMDLENGTH,sizeof(char));
  diffloc = (char *)calloc(MAXCMDLENGTH,sizeof(char));
  units = (char **)calloc(2,sizeof(char *));
  maptarget = (char **)calloc(2,sizeof(char *));
  bmaj = (float *)calloc(2,sizeof(float));
  bmin = (float *)calloc(2,sizeof(float));
  bpa = (float *)calloc(2,sizeof(float));
  beamarea = (float *)calloc(2,sizeof(float));


  frequency = (float *)calloc(2,sizeof(float));
  imgbuff = (float **)calloc(2,sizeof(float *));

  for(i=0; i < 2; i++) {
  
    //Make the full path name for each map
    sprintf(fileloc,"%s", fitsfiles[i]);
    
    //Open the connection
    fitsconn = FunOpen(fileloc, "r", NULL);
 
    //Check we can access the files
    if (fitsconn == NULL){
      fprintf(stderr,"\n*** Error: Unable to open %s, please check the file exists and access permissions are correct. ***\n\n", fitsfiles[i]);
      return 404;
    }

    if( !(imgbuff[i] = FunImageGet(fitsconn, NULL, "bitpix=-32")) ) {
      fprintf(stderr, "\n*** Error: Unable to use ImageGet on %s, please check the file is in the correct format. ***\n\n", fitsfiles[i]);
      return 404;
    }

    /* Get map details and check that they are suitable for use */
  
    FunInfoGet(fitsconn, FUN_SECT_DIM1, &xdim[i], FUN_SECT_DIM2, &ydim[i], 0);
    if ((i > 0) && ((xdim[0] != xdim[i]) || (ydim[0] != ydim[i]))) {
      fprintf(stderr,"\n*** Error: Dimensions of %s do not match those of %s, the map sizes must be the same! ***\n\n", fitsfiles[i], fitsfiles[0] );
      return 4;
    }

    maptarget[i] = FunParamGets(fitsconn, "OBJECT", 0, NULL, &got);
 
    if (strcasecmp(maptarget[0], maptarget[i])) {
      fprintf(stderr,"\n*** Warning: The target of %s (as specified in the header) does not match that of %s! This may simply be due to different naming conventions between maps, but should be checked. ***\n\n", fitsfiles[i], fitsfiles[0]);
    }
 
    units[i] = FunParamGets(fitsconn, "BUNIT", 0, NULL, &got);
    if (strcasecmp("JY/BEAM", units[i])) {
      fprintf(stderr,"\n*** Error: The units used in the map %s (as specified in the header) are not JY/BEAM! ***\n\n", fitsfiles[i] );
      return 4;
    }

    frequency[i] = FunParamGetd(fitsconn, "CRVAL3", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain the frequency of %s! ***\n\n", fitsfiles[i] );
      return 4;  
    }
  
    delt1=FunParamGetd(fitsconn, "CDELT1", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to determine the pixel size! ***\n\n");
      return 4;
    }

    /*if (got) {
      cellsize = 3600 * delt1;
      } else {
      fprintf(stderr,"\n*** Error: Unable to determine the pixel size! ***\n\n");
      return 4;
      }*/
      
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

    FunInfoGet(fitsconn, FUN_HEADER, &head, 0);
      
    // Automatically determine the header type

    str_casadata=FunParamGets(fitsconn, "ORIGIN", 0, NULL, &got);

    if ( strncmp(str_casadata, str_compare_aips, 4) == 0 ) {
      tmp_casadata = 0;
    }
    else if ( strncmp(str_casadata, str_compare_casa, 4) == 0) {
      tmp_casadata = 1;
    }
    else {
      printf("Unable to automatically determine map origin. Assuming the manually set header type\n");
      tmp_casadata = casadata;
    }

    if(tmp_casadata == 0) {

      // This deals with AIPS's crazy header encoding for beam parameters
      for(m=1; m<=head->ncard; m++){
	card = ft_cardnth(head,m);
	if( card && card->c && *card->c){
	  if (!strncmp(card->c, "HISTORY AIPS   CLEAN BMAJ",25)) {
	    sscanf(card->c+26,"%f",&bmaj[i]);
	    sscanf(card->c+44,"%f",&bmin[i]);
	    sscanf(card->c+61,"%f",&bpa[i]);
	    bmaj[i]*=3600; // Convert to arcseconds
	    bmin[i]*=3600;
	    bpa[i] += 90.0; // Convert from AIPS orientation to (0 is Y axis) to plot format (0 is X axis)
	    //printf("BMAJ: %.3f BMIN %.3f BPA %.2f \n", bmaj[i], bmin[i], bpa[i]);
	  }
	}
      }
    }
    else {
      bmaj[i] = FunParamGetd(fitsconn, "BMAJ", 0, 0.0, &got);
      bmin[i] = FunParamGetd(fitsconn, "BMIN", 0, 0.0, &got);
      bmin[i] = FunParamGetd(fitsconn, "BPA", 0, 0.0, &got); // Check this for CASA!
      bmaj[i]*=3600;
      bmin[i]*=3600;
      bpa[i] += 90.0; // Check this for CASA
      //printf("BMAJ: %.3f BMIN %.3f \n", bmaj[i], bmin[i]);
    }

    if ( (fabs(bmaj[i]) - fabs(bmaj[0]) ) > (0.01*fabs(bmaj[0]) ) || (fabs(bmin[i]) - fabs(bmin[0]) ) > (0.01*fabs(bmin[0]) ) ) {
      fprintf(stderr,"\n*** Error: Beam parameters for %s and %s don't match! ***\n\n", fitsfiles[0], fitsfiles[i]);
      return 4;
    }

    if ( (fabs(bpa[i]) - fabs(bpa[0]) ) > (0.01 * fabs(bpa[0]) ) ) {
      fprintf(stderr,"\n*** Error: Beam position angle for %s and %s don't match! ***\n\n", fitsfiles[0], fitsfiles[i]);
      return 4;
    }

    // Parse the headers for encoding of the CLEAN beam size
    if (bmaj[i]>0.0) {
      beamarea[i]=2.0*PI*(bmaj[i]*bmin[i])/(GFACTOR*GFACTOR*delt1*delt1*3600*3600);
    } else {
      fprintf(stderr,"\n*** Failed to find CLEAN beam parameters -- is this a radio map? ***\n\n");
      return 4;
    }

    FunClose(fitsconn);
  }


  // Output what we found (More to do)

  printf("========================================================================\n");

  printf("Pixel size is %.2e arcsec\n", delt1);
  printf("Beam area is %.2f pixels\n", beamarea[0]);
  printf("Beam major axis %.2f arcsec, minor axis %.2f arcsec\n", bmaj[0], bmin[0]);
  printf("The map dimensions are %d by %d\n", xdim[0], ydim[0]);

  printf("========================================================================\n\n");


  printf("Loading map data...\n");

  subtractedmap = (float *)calloc(xdim[0]*ydim[0], sizeof(float));

  datamax = -1e32;
  datamin = 1e32;

  // Loop through by map and find the source flux
  if (useregions == 0) {
    for (a=0; a < 2; a++) {

      // Loop through by pixel
      for (i=0; i<xdim[0]; i++){
	for (j=0; j<ydim[0]; j++){

	  if (a == 0) {
	    subtractedmap[i*xdim[0]+j] = imgbuff[a][i*xdim[0]+j];
	  }
	  else {
	    subtractedmap[i*xdim[0]+j] -= imgbuff[a][i*xdim[0]+j];

	    //Set the min/max values
	    if (subtractedmap[i*xdim[0]+j] > datamax) {
	      datamax = subtractedmap[i*xdim[0]+j];
	    }

	    if (subtractedmap[i*xdim[0]+j] < datamin) {
	      datamin = subtractedmap[i*xdim[0]+j];
	    }

	  }
	}
      }

      // Clean up the image buffer as we go
      free(imgbuff[a]);
    }
  }

  else {

    for (a=0; a < 2; a++) {

      // Loop through by pixel
      for (i=0; i<xdim[0]; i++){
	for (j=0; j<ydim[0]; j++){
	  if ( (i < reg_xdim[currentset]) && (j < reg_ydim[currentset]) ) {
	    if (regionarray[currentset][i][j] > 0) {
	      if (a == 0) {
		subtractedmap[i*xdim[0]+j] = imgbuff[a][i*xdim[0]+j];
	      }
	      else {
		subtractedmap[i*xdim[0]+j] -= imgbuff[a][i*xdim[0]+j];
	      }
	    }
	  }
	}
      }

      // Clean up the image buffer as we go
      free(imgbuff[a]);
    }
  }

  // Create the map

  //Get the current date and time
  time( &currenttime );
  time_struct = localtime ( &currenttime );

  strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);
    
  sprintf(imagename,"%s_%.2e_subtracted_%s.fits",  maptarget[0], frequency[0], timebuff);
  sprintf(diffloc,"%s/%s", imageloc, imagename);

  // Create the file and write the header

  funcomb = FunOpen(diffloc, "w", NULL);

  // Check we can access the file
  if (funcomb == NULL){
    fprintf(stderr,"\n*** Error: Unable to open %s, please check the file exists and access permissions are correct. ***\n\n", fitsfiles[i]);
    return 404;
  }
 

  //Reformat the datetime to that required for the header
  strftime(datebuff,32,"%Y-%m-%d", time_struct);

  //Retrieve the remaining header info needed
  sprintf(fileloc,"%s", fitsfiles[0]);
  fitsconn = FunOpen(fileloc, "r", NULL);

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

  //Set the new header

  FunParamPutb(funcomb, "SIMPLE", 0, 1, NULL, 1);
  FunParamPuti(funcomb, "BITPIX", 0, -32, NULL, 1);
  FunParamPuti(funcomb, "NAXIS", 0, 4, "", 1);
  FunParamPuti(funcomb, "NAXIS1", 0, xdim[0], "X dimension (px)", 1);
  FunParamPuti(funcomb, "NAXIS2", 0, ydim[0], "Y dimension (px)", 1);
  FunParamPuti(funcomb, "NAXIS3", 0, 1, NULL, 1);
  FunParamPuti(funcomb, "NAXIS4", 0, 1, NULL, 1);
  FunParamPutb(funcomb, "EXTEND", 0, 1, "Tables following main image", 1);
  FunParamPutb(funcomb, "BLOCKED", 0, 1, "Tape may be blocked", 1);
  FunParamPuts(funcomb, "OBJECT", 0,  maptarget[0], "Source name", 1);
  FunParamPuts(funcomb, "TELESCOPE", 0, "BRATS", "Broadband Radio Astronomy Tools", 1);
  FunParamPuts(funcomb, "INSTRUME", 0, "", NULL, 1);
  FunParamPuts(funcomb, "OBSERVER", 0, "BRATS", "Broadband Radio Astronomy Tools", 1);
  //FunParamPuts(funcomb,"DATE-OBS", 0, "", "Obs start date YYYY-MM-DD", 1);
  FunParamPuts(funcomb, "DATE-MAP", 0, datebuff, "Last processing date YYYY-MM-DD", 1);
  FunParamPutd(funcomb, "BSCALE", 0, 1.0, 11, "REAL = TAPE * BSCALE + BZERO", 1);
  FunParamPutd(funcomb, "BZERO", 0, 0.0, 9, NULL, 1);
  FunParamPuts(funcomb, "BUNIT", 0,  units[0], "Units of flux", 1);
  FunParamPutd(funcomb, "EQUINOX", 0, eq, 11, "Epoch of RA DEC", 1);
  FunParamPutd(funcomb, "OBSRA", 0, ra, 11, "Antenna pointing RA", 1);
  FunParamPutd(funcomb, "OBSDEC", 0, dec, 11, "Antenna pointing DEC", 1);
  FunParamPutd(funcomb, "DATAMAX", 0, datamax, 11, "Maximum pixel value", 1);
  FunParamPutd(funcomb, "DATAMIN", 0, datamin, 11, "Minimum pixel value", 1);
  FunParamPuts(funcomb, "CTYPE1", 0, "RA---SIN", NULL, 1);
  FunParamPutd(funcomb, "CRVAL1", 0, ra, 11, "RA", 1);
  FunParamPutd(funcomb, "CDELT1", 0, delt1, 11, NULL, 1);
  FunParamPutd(funcomb, "CRPIX1", 0, crpix1, 11, NULL, 1);
  FunParamPutd(funcomb, "CROTA1", 0, crota1, 11, NULL, 1);
  FunParamPuts(funcomb, "CTYPE2", 0, "DEC--SIN", NULL, 1);
  FunParamPutd(funcomb, "CRVAL2", 0, dec, 11, "Dec", 1);
  FunParamPutd(funcomb, "CDELT2", 0, delt2, 11, NULL, 1);
  FunParamPutd(funcomb, "CRPIX2", 0, crpix2, 11, NULL, 1);
  FunParamPutd(funcomb, "CROTA2", 0, crota2, 11, NULL, 1);
  FunParamPuts(funcomb, "CTYPE3", 0, "FREQ", NULL, 1);
  FunParamPutd(funcomb, "CRVAL3", 0, 0.0, 11, "Frequency", 1);
  FunParamPutd(funcomb, "CDELT3", 0, 0.0, 11, NULL, 1);
  FunParamPutd(funcomb, "CRPIX3", 0, 0.0, 11, NULL, 1);
  FunParamPutd(funcomb, "CROTA3", 0, 0.0, 11, NULL, 1);
  FunParamPuts(funcomb, "CTYPE4", 0, "STOKES", NULL, 1);
  FunParamPutd(funcomb, "CRVAL4", 0, crval4, 11, "Stokes", 1);
  FunParamPutd(funcomb, "CDELT4", 0, delt4, 11, NULL, 1);
  FunParamPutd(funcomb, "CRPIX4", 0, crpix4, 11, NULL, 1);
  FunParamPutd(funcomb, "CROTA4", 0, crota4, 11, NULL, 1);
  FunParamPuts(funcomb, "ORIGIN", 0, "BRATS", "Broadband Radio Astronomy Tools", 1);


  strcpy(historytext, "--------------------------------------------------------------------");
  FunParamPuts(funcomb, "HISTORY", 0, historytext, NULL, 1);
  strcpy(historytext, "BRATS  Created by the Broadband Radio Astronomy Tools");
  FunParamPuts(funcomb, "HISTORY", 0, historytext, NULL, 1);
  strcpy(historytext, "BRATS  http://www.askanastronomer.co.uk/brats");
  FunParamPuts(funcomb, "HISTORY", 0, historytext, NULL, 1);
  strcpy(historytext, "BRATS  Difference map between ");
  sprintf(dbltochar, "%.4e", frequency[0]);
  strcat(historytext, dbltochar);
  strcat(historytext, " and ");
  sprintf(dbltochar, "%.4e", frequency[1]);
  strcat(historytext, dbltochar);
  strcat(historytext, " Hz");
  FunParamPuts(funcomb, "HISTORY", 0, historytext, NULL, 1);


  // Output the image
  FunImagePut(funcomb, subtractedmap, xdim[0], ydim[0], -32, NULL);
  FunClose(funcomb);
  //FunClose(fitsconn);
    
  printf("\nDifference map has been exported to %s\n\n", diffloc);

  //Free up the memory
  for (i=0; i<2; i++) {
    free(fitsfiles[i]);
  }

  free(maptarget);
  free(fileloc);
  free(units);
  free(bmaj);
  free(bmin);
  free(beamarea);
  free(imgbuff);
  free(xdim);
  free(ydim);
  free(fitsfiles);

  return 0;

}
