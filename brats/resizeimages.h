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

#include <time.h>


// Declare the function prototypes

int syscom(const char[]);

int resizeimages(char *imageloc) { 

  // Declare the variables for this task

  register char **fitsfiles;

  int required, i, a, x, y, got, got2, got3, *xdim, *ydim, imgnum, newxdim, newydim, xshift, yshift, freq_ok;
  float **imgbuff, *frequency;
  double newcentrex, newcentrey, centrex, centrey;
  char *fileloc, **maptarget, dirname[MAXCMDLENGTH], *newloc, imagename[256], timebuff[32], *endptr, *cmdbuffer, historytext[MAXCMDLENGTH], inttochar[MAXCMDLENGTH], *tmp_valcheck;

  FITSHead head;
  Fun *fitsconn, *funnew;
  DIR *direct;

  struct dirent *entry;

  // Set up variables and structs for times and dates

  time_t currenttime;
  struct tm * time_struct;

  //Get the directory name and format it
 
  required = 1;

  while (required == 1) {

    strcpy(dirname, "");

    cmdbuffer = readline("Enter directory name containing the maps to be resized: ");
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
    cmdbuffer = readline("Please enter the size in pixels for the new X-axis (-1 to escape): ");
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    chomp(cmdbuffer);

    newxdim = strtol(cmdbuffer, &endptr, 10);

    if ( (newxdim == 0) && (cmdbuffer == endptr) ) {
      printf("\nInvalid input, please try again...\n\n");
      continue;
    }
    else if (newxdim < 0.0) {
      printf("Escaping command...\n\n");
      return(100);
    }
    else if ( (newxdim < 10.0) && (newxdim < 10.0) ) {
      printf("\nThis is too small to be sensible. Please try again...\n\n");
      continue;
    }
    else {
      required = 0;
    }
  }


  required = 1;

  while (required == 1) {

    // Get the scaling factor
    cmdbuffer = readline("Please enter the size in pixels for the new Y-axis (-1 to escape): ");
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    chomp(cmdbuffer);

    newydim= strtol(cmdbuffer, &endptr, 10);

    if ( (newydim == 0) && (cmdbuffer == endptr) ) {
      printf("\nInvalid input, please try again...\n\n");
      continue;
    }
    else if (newydim < 0.0) {
      printf("Escaping command...\n\n");
      return(100);
    }
    else if ( (newydim < 10.0) && (newydim < 10.0) ) {
      printf("\nThis is too small to be sensible. Please try again...\n\n");
      continue;
    }
    else {
      required = 0;
    }
  }

  printf("The dimensions of the new map will be %d x %d pixels\n", newxdim, newydim);


  // Put some space between the input and the output
  printf("\n");

  printf("Loading map data...\n");

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
    fprintf(stderr,"\n*** Error: Unable to locate FITS directory ***\n\n");
    return 404;
  }

  /* Start looping through the maps and getting the data */

  fileloc = (char *)calloc(MAXCMDLENGTH,sizeof(char));
  newloc = (char *)calloc(MAXCMDLENGTH,sizeof(char));
  maptarget = (char **)calloc(imgnum,sizeof(char *));
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

    // Get the current map dimensions
  
    FunInfoGet(fitsconn, FUN_SECT_DIM1, &xdim[i], FUN_SECT_DIM2, &ydim[i], 0);
    printf("Map %d currently has the dimensions %d x %d\n", i, xdim[i], ydim[i]);

    maptarget[i] = FunParamGets(fitsconn, "OBJECT", 0, NULL, &got);
    if (!got) {
      fprintf(stderr,"\n*** Warning: Unable to obtain the source name of %s! Setting as UNKNOWN ***\n\n", fitsfiles[i] );
      maptarget[i] = (char *)calloc(10,sizeof(char));
      strcpy(maptarget[i], "UNKNOWN");
    }

    /* frequency[i] = FunParamGetd(fitsconn, "CRVAL3", 0, 0.0, &got); */
    /* if (!got) { */
    /*   fprintf(stderr,"\n*** Error: Unable to obtain the frequency of %s! ***\n\n", fitsfiles[i] ); */
    /*   return 4;   */
    /* } */

    freq_ok = 0;
    frequency[i] = FunParamGetd(fitsconn, "REFFREQ", 0, 0.0, &got);

    if (!got) {

      tmp_valcheck = FunParamGets(fitsconn, "CTYPE3", 0, NULL, &got2);

      if (!got2) { // Check we actually have something

	frequency[i] = atof(FunParamGets(fitsconn, "FREQ", 0, NULL, &got3)); // Used by GLEAM, is a string so must be converted

	if (!got3) { // If thats not there....
	
	  printf("I've ran out of ideas where the observed frequency could be located. Please contact the developer to have your header style added...\n");
	  fprintf(stderr,"\n*** Error: Unable to obtain the frequency of %s! ***\n\n", fitsfiles[i] );
	  return 4;
	}
	else {
	  freq_ok = 1;
	}
      }
      else {
	chomp(tmp_valcheck);

	if (strcasecmp("FREQ", tmp_valcheck) == 0) {
	  frequency[i] = FunParamGetd(fitsconn, "CRVAL3", 0, 0.0, &got);
	
	  if (!got) {

	    printf("I've managed to find the frequency under CTYPE3, but the accompanying CRVAL3 appears to be missing. Unable to proceed. Please add the values to the header or contact the developer to have your header style added...\n");
	    fprintf(stderr,"\n*** Error: Unable to obtain the frequency of %s! ***\n\n", fitsfiles[i] );
	    return 4;

	  }
	  else {
	    freq_ok = 1;
	  }
	}

	if (freq_ok == 0) {
	  printf("I've ran out of ideas where the observed frequency could be located. Please contact the developer to have your header style added...\n");
	  fprintf(stderr,"\n*** Error: Unable to obtain the frequency of %s! ***\n\n", fitsfiles[i] );
	  return 4;  
	}
      }
    }





    centrex = FunParamGetd(fitsconn, "CRPIX1", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain the centre pixel (x-axis) of %s! ***\n\n", fitsfiles[i] );
      return 4;  
    }

    centrey = FunParamGetd(fitsconn, "CRPIX2", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain the centre pixel (y-axis) of %s! ***\n\n", fitsfiles[i] );
      return 4;  
    }

    FunClose(fitsconn);
  }


  // Loop through by mapz and copy over the flux to the new location

  printf("Creating the resized maps...\n");

  for (a=0; a < imgnum; a++) {

    float *newmap, *newmapP1;
    newmap = (float *)calloc(newxdim*newydim, sizeof(float));

    // Set the shift needed to put the data in the middle of the new map size
    xshift=(newxdim-xdim[a])/2;
    yshift=(newydim-ydim[a])/2;

    // Loop through by pixel

    // If there is no change, skip resizing
    if ( (newxdim == xdim[a]) && (newydim == ydim[a]) ){
      printf("The new dimensions match those of this original. Skipping map %d...\n", a);
      continue;
    }

    // Else if we are contracting in the both directions...
    else if ( (newxdim < xdim[a]) && (newydim < ydim[a]) ) {

      // Do a 2 step process, contract in one direction then in the other. This allows for ANY contacted dimensions to be used, rather than just squares
      newmapP1 = (float *)calloc(xdim[a]*newydim, sizeof(float));

      // Fixing X to the original size and contacting in Y
      for (x=0; x<xdim[a]; x++) {
	for (y=0; y<newydim; y++) {
	  newmapP1[(y*xdim[a])+x] = imgbuff[a][(abs(yshift)*xdim[a])+(y*xdim[a])+x];
	}
      }

      // Now contacting the new array in X with a fixed Y...
      for (x=0; x<newxdim; x++){
	for (y=0; y<newydim; y++){
	  newmap[(y*newxdim)+x] = newmapP1[(x-xshift)+(xdim[a]*y)];
	}
      }
    }

    // Else if we are expanding in the both directions...
    else if ( (newxdim > xdim[a]) && (newydim > ydim[a]) ) {

      newmapP1 = (float *)calloc(xdim[a]*newydim, sizeof(float));

      // Fixing X to the original size and expanding in Y
      for (x=0; x<xdim[a]; x++) {
	for (y=0; y<ydim[a]; y++) {
	  //newmapP1[yshift+(y*xdim[a])+x] = imgbuff[a][(x*xdim[a])+y];
	  newmapP1[(yshift*xdim[a])+(y*xdim[a])+x] = imgbuff[a][(y*xdim[a])+x];
	}
      }

      // Now expanding the new array in X with a fixed Y...
      for (x=0; x<xdim[a]; x++){
	for (y=0; y<newydim; y++){
	  newmap[xshift+(y*newxdim)+x] = newmapP1[x+(xdim[a]*y)];
	}
      }
    }


    // Else if we are contracting in the x and expanding (or fixed) in the y directions...
    else if ( (newxdim < xdim[a]) && (newydim >= ydim[a]) ){
      for (x=0; x<newxdim; x++){
	for (y=yshift; y<(yshift+ydim[a]); y++){
	  newmap[(newxdim*yshift)+((y-yshift)*newxdim)+x] = imgbuff[a][abs(xshift)+((y-yshift)*xdim[a])+x];
	}
      }
    }

    // Else if we are contracting in the y and expanding (or fixed) in the x directions...
    else if ( (newxdim >= xdim[a]) && (newydim < ydim[a]) ){
      for (x=xshift; x<(xshift+xdim[a]); x++) {
	for (y=0; y<newydim; y++) {
	  newmap[(xshift+(y*newxdim)+(x-xshift))] = imgbuff[a][(abs(yshift)*xdim[a])+(y*xdim[a])+(x-xshift)];
	}
      }
    }

    // Create the new map

    // Make a time stamp and create the image name
    time( &currenttime );
    time_struct = localtime ( &currenttime );

    strftime(timebuff,32,"%d%m%y%H%M%S", time_struct);

    sprintf(imagename,"%s_%d_%dby%d_%s_resized.fits", maptarget[a], a, newxdim, newydim, timebuff);
    sprintf(newloc,"%s/%s", imageloc, imagename);

    // Create the file and write the header (currently AIPS format)
    funnew = FunOpen(newloc, "w", NULL);

    // Check we can access the file
    if (funnew == NULL){
      fprintf(stderr,"\n*** Error: Unable to open %s, please check the file exists and access permissions are correct. ***\n\n", newloc);
      return 404;
    }

    // Get the header of the first map loaded
    sprintf(fileloc,"./%s/%s", dirname, fitsfiles[a]);
    fitsconn = FunOpen(fileloc, "r", NULL);
    FunInfoGet(fitsconn, FUN_HEADER, &head, 0);

    // Set this as the header for the new map
    FunInfoPut(funnew, FUN_HEADER, &head, 0);

    newcentrex = centrex + xshift;
    newcentrey = centrey + yshift;

    //printf("xshift: %d yshift: %d centrex: %.4e centrey: %.4e newcentrex: %.4e newcentrey: %.4e\n", xshift, yshift, centrex, centrey,newcentrex ,newcentrey );

    // Set the new dimensions
    FunParamPuti(funnew, "NAXIS1", 0, newxdim, NULL, 1);
    FunParamPuti(funnew, "NAXIS2", 0, newydim, NULL, 1);
    FunParamPutd(funnew, "CRPIX1", 0, newcentrex, 11, NULL, 1);
    FunParamPutd(funnew, "CRPIX2", 0, newcentrey, 11, NULL, 1);

    // Add a line to the history saying what we have done
    strcpy(historytext, "BRATS  Image resized using the Broadband Radio Astronomy Tools");
    FunParamPuts(funnew, "HISTORY", 0, historytext, NULL, 1);
    strcpy(historytext, "BRATS  http://www.askanastronomer.co.uk/brats");
    FunParamPuts(funnew, "HISTORY", 0, historytext, NULL, 1);
    strcpy(historytext, "BRATS  Image resized to ");
    sprintf(inttochar, "%d", newxdim);
    strcat(historytext, inttochar);
    strcat(historytext, " x ");
    sprintf(inttochar, "%d", newydim);
    strcat(historytext, inttochar);
    strcat(historytext, " px");

    FunParamPuts(funnew, "HISTORY", 0, historytext, NULL, 1);


    // Output the image
    FunImagePut(funnew, newmap, newxdim, newydim, -32, NULL);
    FunClose(funnew);
    FunClose(fitsconn);
    
    printf("Image %d has been resized and exported to %s\n", a, newloc);

    // Clean up the image buffer as we go
    free(imgbuff[a]);
    free(newmap);
  }


  //Free up the memory
  for (i=0; i<MAXNUMFILES; i++) {
    free(fitsfiles[i]);
  }

  free(maptarget);
  free(fileloc);
  free(imgbuff);
  free(xdim);
  free(ydim);
  free(fitsfiles);

  return 0;

}
