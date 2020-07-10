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

int mergeimages(int casadata, char *imageloc) { 


  // Declare the variables for this task

  register char **fitsfiles;

  int required, i, j, a, got, m, *xdim, *ydim, imgnum, tmp_casadata;
  float *bmaj, *bmin, *bpa, *beamarea, **imgbuff, scalingfactor, *combinedmap, *ra, *dec, *frequency, cellsize;
  double delt, delt2, rapix, decpix;
  char *fileloc, **units, **maptarget, dirname[MAXCMDLENGTH], *combloc, imagename[256], *endptr, *cmdbuffer;

  FITSHead head;
  FITSCard card;
  Fun *fitsconn, *funcomb;
  DIR *direct;

  float combfreq = 1e9;

  char *str_casadata, str_compare_aips[8], str_compare_casa[8];
  strcpy(str_compare_aips, "AIPS");
  strcpy(str_compare_casa, "CASA");

  struct dirent *entry;

  //Get the directory name and format it
 
  required = 1;

  while (required == 1) {

    strcpy(dirname, "");

    cmdbuffer = readline("Enter directory name containing the maps to be merged: ");
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    chomp(cmdbuffer);

    strcpy(dirname, cmdbuffer);

    if ( (dirname == 0) && (cmdbuffer == endptr) ) {
      printf("\nInvalid input, please try again...\n\n");
      continue;
    }

    // Check a directory name has been enetered isn't over the maximum character limit
    if (strlen(dirname) >= MAXCMDLENGTH-1) {
      fprintf(stderr,"\n*** Error: Folder path must be less than %d characters long! ***\n\n", MAXCMDLENGTH);
      return(2);
    }
    else if (strlen(dirname) == 0) {
      continue;
    }
    else if ((strstr(dirname, "esc") != NULL) && (strlen(dirname) == 3)) {
      printf("Escaping command...\n");
      return 100;
    }
    else if ((strstr(dirname, "ls") != NULL) && (strlen(dirname) == 2)) {
      syscom("ls");
    }
    else {
      required = 0;
    }
  }

  required = 1;

  while (required == 1) {

    // Get the scaling factor

    cmdbuffer = readline("Please enter a frequency (in Hz) for the combined image (-1 to escape): ");
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    chomp(cmdbuffer);

    combfreq = strtod(cmdbuffer, &endptr);

    if ( (combfreq == 0) && (cmdbuffer == endptr) ) {
      printf("\nInvalid input, please try again...\n\n");
      continue;
    }
    else if (combfreq < 0.0) {
      printf("Escaping command...\n\n");
      return(100);
    }
    else {
      printf("Combined map set to a frequency of %.2e Hz\n", combfreq);
      required = 0;
    }
  }

  required = 1;

  while (required == 1) {

    // Get the scaling factor
    cmdbuffer = readline("Please enter a value for the power law scaling (-1 to escape): ");
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    chomp(cmdbuffer);

    scalingfactor = strtod(cmdbuffer, &endptr);

    if ( (scalingfactor == 0) && (cmdbuffer == endptr) ) {
      printf("\nInvalid input, please try again...\n\n");
      continue;
    }
    else if (scalingfactor < 0.0) {
      printf("Escaping command...\n\n");
      return(100);
    }
    else {
      printf("Scaling factor set to %.2f\n", scalingfactor);
      required = 0;
    }
  }


  // Put some space between the input and the output
  printf("\n");
  /* Check all the files and directories exists then put them in the correct format  */

  // Setting image counter and allocating memory
  imgnum = 0;

  xdim = (int *)calloc(MAXNUMFILES, sizeof(int));
  ydim = (int *)calloc(MAXNUMFILES, sizeof(int));
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
	imgnum++;

	// Check we havent exceeeded the maxmim number of files we are allowed 
	if (imgnum >= MAXNUMFILES) {
	  fprintf(stderr,"\n*** Error: A maximum of %d files can currently be processed ***\n\n", MAXNUMFILES);
	  return 3;
	}
      }
    }

    closedir (direct);

  }

  //If we can't find or access the directory, kick back to the command prompt
  else {
    fprintf(stderr,"*** Error: Unable to locate FITS directory ***\n\n");
    return 404;
  }

  /* Start looping through the maps and getting the data */

  fileloc = (char *)calloc(MAXCMDLENGTH,sizeof(char));
  combloc = (char *)calloc(MAXCMDLENGTH,sizeof(char));
  units = (char **)calloc(imgnum,sizeof(char *));
  maptarget = (char **)calloc(imgnum,sizeof(char *));
  bmaj = (float *)calloc(imgnum,sizeof(float));
  bmin = (float *)calloc(imgnum,sizeof(float));
  bpa = (float *)calloc(imgnum,sizeof(float));
  beamarea = (float *)calloc(imgnum,sizeof(float));
  ra = (float *)calloc(imgnum,sizeof(float));
  dec = (float *)calloc(imgnum,sizeof(float));


  frequency = (float *)calloc(imgnum,sizeof(float));
  imgbuff = (float **)calloc(imgnum,sizeof(float *));

  for(i=0; i < imgnum; i++) {
  
    //Make the full path name for each map
    sprintf(fileloc,"./%s/%s", dirname, fitsfiles[i]);
    
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
 
    if (strcmp(maptarget[0], maptarget[i])) {
      fprintf(stderr,"\n*** Warning: The target of %s (as specified in the header) does not match that of %s! This may simply be due to different naming conventions between maps, but should be checked. ***\n\n", fitsfiles[i], fitsfiles[0]);
    }
 
    units[i] = FunParamGets(fitsconn, "BUNIT", 0, NULL, &got);
    if (strcmp("JY/BEAM", units[i])) {
      fprintf(stderr,"\n*** Error: The units used in the map %s (as specified in the header) are not JY/BEAM! ***\n\n", fitsfiles[i] );
      return 4;
    }

    frequency[i] = FunParamGetd(fitsconn, "CRVAL3", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain the frequency of %s! ***\n\n", fitsfiles[i] );
      return 4;  
    }
   
    if (i > 0) {
      for (a=0; a < i; a++) {
	if (frequency[a] == frequency[i]) {
	  fprintf(stderr,"\n*** Error: The frequency of the maps %s and %s is the same. This would lead to some odd results! ***\n\n", fitsfiles[a], fitsfiles[i]);
	  return 4; 
	}
      }
    }
  
    delt=FunParamGetd(fitsconn, "CDELT2", 0, 0.0, &got);
      
    if (got) {
      cellsize = 3600 * delt;
    } else {
      fprintf(stderr,"\n*** Error: Unable to determine the pixel size! ***\n\n");
      return 4;
    }
      
    delt2=FunParamGetd(fitsconn, "CDELT1", 0, 0.0, &got);
      
    if (got) {
    
      if (delt2!=-delt) {
	fprintf(stderr,"\n*** Error: Pixels aren't square! ***\n\n");
	return 4;
      }
    } else {
      fprintf(stderr,"\n*** Error: Unable to determine the second pixel size! ***\n\n");
      return 4;
    }
      
    delt*=3600.0;

    FunInfoGet(fitsconn, FUN_HEADER, &head, 0);
      
    // At least all params below need adjusting to automatically accept to CASA data

    ra[i] = FunParamGetd(fitsconn, "OBSRA", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain the RA of %s! ***\n\n", fitsfiles[i] );
      return 4;  
    }

    dec[i] = FunParamGetd(fitsconn, "OBSDEC", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain the Dec of %s! ***\n\n", fitsfiles[i] );
      return 4;  
    }

    // Adjust the RA and DEC so they are on the centre pixel

    rapix = FunParamGetd(fitsconn, "CRPIX1", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain the pixel number for the RA of %s! ***\n\n", fitsfiles[i] );
      return 4;  
    }

    ra[i] += ( ( (xdim[i]/2) - rapix) * (cellsize / 3600) );


   decpix = FunParamGetd(fitsconn, "CRPIX2", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain the pixel number for the DEC of %s! ***\n\n", fitsfiles[i] );
      return 4;
    }

    dec[i] += ( ( (xdim[i]/2) - decpix) * ( cellsize / 3600) );


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
      bmin[i] = FunParamGetd(fitsconn, "BPA", 0, 0.0, &got);
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
      beamarea[i]=2.0*PI*(bmaj[i]*bmin[i])/(GFACTOR*GFACTOR*delt*delt);
    } else {
      fprintf(stderr,"\n*** Failed to find CLEAN beam parameters -- is this a radio map? ***\n\n");
      return 4;
    }

    FunClose(fitsconn);
  }


  // Output what we found (More to do)

  printf("========================================================================\n");

  printf("There were %d maps found in the directory %s\n", imgnum, dirname);
  printf("Pixel size is %.2e arcsec\n", delt);
  printf("Beam area is %.2f pixels\n", beamarea[0]);
  printf("Beam major axis %.2f arcsec, minor axis %.2f arcsec\n", bmaj[0], bmin[0]);
  printf("The map dimensions are %d by %d\n", xdim[0], ydim[0]);

  printf("========================================================================\n\n");


  printf("Loading map data...\n");

  combinedmap = (float *)calloc(xdim[0]*ydim[0], sizeof(float));


  // Loop through by map and find the source flux

   for (a=0; a < imgnum; a++) {

    // Loop through by pixel
    for (i=(0); i<xdim[0]; i++){
      for (j=(0); j<ydim[0]; j++){

	combinedmap[i*xdim[0]+j] += (imgbuff[a][i*xdim[0]+j] * pow((frequency[a]/combfreq), scalingfactor) );
      }
    }

    // Clean up the image buffer as we go
    free(imgbuff[a]);
   }


   // Average the flux
   for (i=(0); i<xdim[0]; i++){
      for (j=(0); j<ydim[0]; j++){
	combinedmap[i*xdim[0]+j] /= imgnum;
      }
    }

   // Create the map

   sprintf(imagename,"%s_%.2e_combined.fits",  maptarget[0], combfreq);
   sprintf(combloc,"./%s/%s", imageloc, imagename);

   // Create the file and write the header (currently AIPS format)
    
    funcomb = FunOpen(combloc, "w", NULL);

    // Check we can access the file
    if (fitsconn == NULL){
      fprintf(stderr,"\n*** Error: Unable to open %s, please check the file exists and access permissions are correct. ***\n\n", fitsfiles[i]);
      return 404;
    }
    
    // Get the header of the first map loaded
    sprintf(fileloc,"./%s/%s", dirname, fitsfiles[0]);
    fitsconn = FunOpen(fileloc, "r", NULL);
    FunInfoGet(fitsconn, FUN_HEADER, &head, 0);

    // Set this as the header for the new map
    FunInfoPut(funcomb, FUN_HEADER, &head, 0);

    // Set any parameters which have been changed.
    /*
    FunParamPuts(funcomb, "OBJECT", 0,  maptarget[0], "Target Source", 1);
    FunParamPuts(funcomb, "BUNIT", 0,  units[0], "Units", 1);
    FunParamPutd(funcomb, "CRVAL3", 0,  (double)combfreq, 11, "Map Frequency", 1);
    FunParamPutd(funcomb, "CDELT1", 0,  (double)delt2, 9, "", 1);
    FunParamPutd(funcomb, "CDELT2", 0,  (double)delt, 9, "", 1);
    */

    FunParamPutd(funcomb, "CRVAL3", 0,  (double)combfreq, 11, "Map Frequency", 1);



    // Output the image
    FunImagePut(funcomb, combinedmap, xdim[0], ydim[0], -32, NULL);
    FunClose(funcomb);
    FunClose(fitsconn);
    
    printf("\nCombined image has been exported to %s\n\n", fileloc);

  //Free up the memory
  for (i=0; i<MAXNUMFILES; i++) {
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
