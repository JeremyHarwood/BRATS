/*

  Loads a series of maps. Must be run first before most other tasks.

  Requires the file "chomp.h", "fluxmap.h" and <funtools.h> plus the standard headers to have been previously declared

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

int fluxmap(float ****flux, float **rms, int *xdim, int *ydim, int mapnumber, float **maxflux, float **minflux, int currentset, int border, float zoom, float sigmalevel, char fluxmaptitle[128], int titles, int labels, int axismarks, int extwincontrol, int smoothmaps, int fluxlogs, int fluxcut, float usercut, float wedgemultiplier, char *output, float ra, float dec, float cellsize, int scaletype, float bmaj, float bmin, float bpa, int beamposx, int beamposy, int plotbeam, int xshift, int yshift, int largetxt, int export);

int load(char dirname[], int currentset, int numdatasets, int *imgnum, int *xdim, int *ydim, float *bmaj, float *bmin, float *bpa, float *beamarea, float **frequency, float **frequencyobs, char regname[], char bgname[], float ****flux, float **rms, float sigmalevel, int **fluxorder,float zoom, int *border, int titles, int labels, int axismarks, char tmp_maptarget[], float **yaxismax, float **yaxismin, int casadata, int fluxlogs, int fluxcut, float usercut, float wedgemultiplier, float **ra, float **dec, float **delt1, float **delt2, float **crpix1, float **crpix2, float **crota1, float **crota2, float **eq, float *cellsize, float scaletype, int posmap, int beamposx, int beamposy, int plotbeam, float ***bgbuff, int xshift, int yshift, float *redshift, int largetxt, float frac_border) {


  //Declare the variables for this task

  register char **fitsfiles, **bgfiles;

  float celltol = 0.000000000000001; // This needs to be made variable

  int usereg, required, i, j, a, got, got2, got3, m, *tmp_xdim, *tmp_ydim, count, tmp_casadata, freq_ok, istart, iend;
  float *tmp_bmaj, *tmp_bmin, *tmp_bpa, *tmp_beamarea, bgflux, **imgbuff, curminfreq, curmaxfreq, tmp_redshift, avgdec, tmp_delt, zoomedx;
  char *fileloc, **units, *bgloc, **maptarget, fluxmaptitle[128], output[16], **tmp_valcheck;

  int suppress = 0; // Suppress repeating messages

  float degtorad = 3.14159/180;
  //float wholecirclesecs= 86400.0;

  FITSHead head;
  FITSCard card;
  Fun *fitsconn, *bgconn;
  DIR *direct;

  char *str_casadata, str_compare_aips[8], str_compare_casa[8], *cmdbuffer, *endptr;
  strcpy(str_compare_aips, "AIPS");
  strcpy(str_compare_casa, "CASA");

  struct dirent *entry;

  // Auto complete paths and file names
  //rl_bind_key('\t', rl_complete);

  //Get the directory name and format it 
  required = 1;

  while (required == 1) {

    strcpy(dirname, "");

    cmdbuffer = readline("Enter directory name: ");
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    strcpy(dirname, cmdbuffer);

    chomp(dirname);

    // Check a directory name has been enetered isn't over the maximum character limit
    if (strlen(dirname) >= MAXCMDLENGTH-1) {
      fprintf(stderr,"\n*** Error: Folder path must be less than %d characters long! ***\n\n", MAXCMDLENGTH);
      return(2);
    }
    else if (strlen(dirname) == 0) {
      printf("Please enter a directory name\n");
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

    //Get the background file name and format it

    cmdbuffer = readline("Background file path: ");
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    strcpy(bgname, cmdbuffer);

    chomp(bgname);

    if (strlen(bgname) >= MAXCMDLENGTH-1) {
      fprintf(stderr,"\n*** Error: Folder path must be less than %d characters long! ***\n\n", MAXCMDLENGTH);
      return 2;
    }
    else if (strlen(bgname) == 0) {
      printf("Please enter a directory name\n");
    }
    else if ((strstr(bgname, "esc") != NULL) && (strlen(bgname) == 3)) {
      printf("Escaping command...\n");
      return 100;
    }
    else if ((strstr(bgname, "ls") != NULL) && (strlen(bgname) == 2)) {
      syscom("ls");
    }
    else {
      required = 0;
    }
  }

  required = 1;

  while (required == 1) {

    //Get the region file name and format it
    cmdbuffer = readline("Region file path (Default: full map): ");
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    strcpy(regname, cmdbuffer);

    chomp(regname);

    if (strlen(regname) >= MAXCMDLENGTH-1) {
      fprintf(stderr,"\n*** Error: File path must be less than %d characters long! ***\n\n", MAXCMDLENGTH);
      return 2;
    }
    else if ((strstr(regname, "esc") != NULL) && (strlen(regname) == 3)) {
      printf("Escaping command...\n");
      return 100;
    }
    else if ((strstr(regname, "ls") != NULL) && (strlen(regname) == 2)) {
      syscom("ls");
    }
    else {
      //Set if a region file is being used or not
      if (strlen(regname) == 0) {
	usereg = 0;
	required = 0;
      }
      else {
	usereg = 1;
	required = 0;
      }
    }
  }

  required = 1;

  
  while (required == 1) {

    cmdbuffer = readline("Redshift of the target source (-1 to exit): ");
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    chomp(cmdbuffer);

    tmp_redshift = strtod(cmdbuffer, &endptr);

    //  printf("Redshift of the target source (-1 to exit): ");

    if ( (tmp_redshift == 0) && (cmdbuffer == endptr) ) {
		 
      printf("\nInvalid input, please try again...\n\n");
      continue;

      //if (scanf("%f", &tmp_redshift) == 0) {
      //printf("Invalid input, escaping command to protect the existing data...\n");
      //return 100;
    }

    if (tmp_redshift < 0) {
      printf("Escaping command...\n");
      return 100;
    }

    else {
      redshift[currentset]=tmp_redshift;
      required = 0;
    }
  }

  // Put some space between the input and the output
  printf("\n");
  /* Check all the files and directories exists then put them in the correct format  */

  // Setting image counter and allocating memory
  imgnum[currentset] = 0;

  tmp_xdim = (int *)calloc(MAXNUMFILES, sizeof(int));
  tmp_ydim = (int *)calloc(MAXNUMFILES, sizeof(int));
  fitsfiles = (char **)calloc(MAXNUMFILES, sizeof(char *));
  bgfiles = (char **)calloc(MAXNUMFILES, sizeof(char *));
 
  
  //Open the directory
  direct = opendir(dirname);

  //Set the map counter to 0
  i = 0;

  //Check the directory exists
  if (direct != NULL) {

    //While there are maps to read...
    while ((entry = readdir(direct))) {
      if (strlen(entry->d_name) > MAXCMDLENGTH+1) {
	fprintf(stderr,"\n*** Error: File path must be less than %d characters long! ***\n\n", MAXCMDLENGTH);
	return 2;
      }
      //Do nothing if we are reading . or ..
      else if ((strcmp(entry->d_name, ".") == 0) || (strcmp(entry->d_name, "..") == 0)) {
	;
      }

      //Append the region and background file name to the map name...
      else {
	fitsfiles[i] = (char *)calloc((MAXCMDLENGTH+1),sizeof(char));
	strcpy(fitsfiles[i], entry->d_name);
	if (usereg != 0){
	  if( access(regname, F_OK) != -1 ) {
	    if( access(bgname, F_OK) != -1 ) {
	      bgfiles[i] = (char *)calloc((MAXCMDLENGTH+1),sizeof(char));
	      sprintf(bgfiles[i],"%s[@%s]", fitsfiles[i], bgname);     
	    }
	    else {
	      fprintf(stderr,"\n*** Error: Cannot open the background region file %s. ***\n\n", bgname);
	      return 404; 
	    }
	    // This must be second so as to pass the file name without the region appended!
	    sprintf(fitsfiles[i],"%s[@%s]", fitsfiles[i], regname);
	  }
	  else {
	    fprintf(stderr,"\n*** Error: Cannot open the region file %s. ***\n\n", regname);
	    return 404; 
	  }
	}

	// ... or if we are doing full map processing, just setup the appended background file name

	else {
	  if( access(bgname, F_OK) != -1 ) {
	    bgfiles[i] = (char *)calloc((MAXCMDLENGTH+1),sizeof(char));
	    sprintf(bgfiles[i],"%s[@%s]", fitsfiles[i], bgname);  
	  }   
	  else {
	    fprintf(stderr,"\n*** Error: Cannot open the background region file %s. ***\n\n", bgname);
	    return 404; 
	  }
	}

	//Increase the map counter
	i++;
	imgnum[currentset]++;

	// Check we havent exceeeded the maxmim number of files we are allowed 
	if (imgnum[currentset] >= MAXNUMFILES) {
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
  units = (char **)calloc(imgnum[currentset],sizeof(char *));
  tmp_valcheck = (char **)calloc(imgnum[currentset],sizeof(char *));
  maptarget = (char **)calloc(imgnum[currentset],sizeof(char *));
  tmp_bmaj = (float *)calloc(imgnum[currentset],sizeof(float));
  tmp_bmin = (float *)calloc(imgnum[currentset],sizeof(float));
  tmp_bpa = (float *)calloc(imgnum[currentset],sizeof(float));
  tmp_beamarea = (float *)calloc(imgnum[currentset],sizeof(float));
  bgloc = (char *)calloc(MAXCMDLENGTH,sizeof(char));
  ra[currentset] = (float *)calloc(imgnum[currentset],sizeof(float));
  dec[currentset] = (float *)calloc(imgnum[currentset],sizeof(float));
  delt1[currentset] = (float *)calloc(imgnum[currentset],sizeof(float));
  delt2[currentset] = (float *)calloc(imgnum[currentset],sizeof(float));
  crpix1[currentset] = (float *)calloc(imgnum[currentset],sizeof(float));
  crpix2[currentset] = (float *)calloc(imgnum[currentset],sizeof(float));
  crota1[currentset] = (float *)calloc(imgnum[currentset],sizeof(float));
  crota2[currentset] = (float *)calloc(imgnum[currentset],sizeof(float));
  eq[currentset] = (float *)calloc(imgnum[currentset],sizeof(float));

  frequency[currentset] = (float *)calloc(imgnum[currentset],sizeof(float));
  frequencyobs[currentset] = (float *)calloc(imgnum[currentset],sizeof(float));
  imgbuff = (float **)calloc(imgnum[currentset],sizeof(float *));
  bgbuff[currentset] =(float **)calloc(imgnum[currentset],sizeof(float *));

  for(i=0; i < imgnum[currentset]; i++) {
  
    //Make the full path name for each map
    sprintf(fileloc,"./%s/%s", dirname, fitsfiles[i]);
    sprintf(bgloc,"./%s/%s", dirname, bgfiles[i]);
    
    //Open the connection
    fitsconn = FunOpen(fileloc, "r", NULL);
 
    //Check we can access the files
    if (fitsconn == NULL){
      fprintf(stderr,"\n*** Error: Unable to open %s, please check the file exists and access permissions are correct. ***\n\n", fitsfiles[i]);
      return 404;
    }

    /* Check that the images are square */
    FunInfoGet(fitsconn, FUN_SECT_DIM1, &tmp_xdim[i], FUN_SECT_DIM2, &tmp_ydim[i], 0);
    if (tmp_xdim[i] != tmp_ydim[i]) {
      fprintf(stderr,"\n*** Error: The image %s is not square (this is currently a requirement). This can be easily corrected by using the resizeimages command to expand all images to the largest dimension on both sides. ***\n\n", fitsfiles[i]);
      return 4;
    }

    if( !(imgbuff[i] = FunImageGet(fitsconn, NULL, "bitpix=-32")) ) {
      fprintf(stderr, "\n*** Error: Unable to use ImageGet on %s, please check the file is in the correct format. ***\n\n", fitsfiles[i]);
      return 404;
    }

    /* Get map details and check that they are suitable for use */
  
    FunInfoGet(fitsconn, FUN_SECT_DIM1, &tmp_xdim[i], FUN_SECT_DIM2, &tmp_ydim[i], 0);
    if ((i > 0) && ((tmp_xdim[0] != tmp_xdim[i]) || (tmp_ydim[0] != tmp_ydim[i]))) {
      fprintf(stderr,"\n*** Error: Dimensions of %s do not match those of %s, the map sizes must be the same! ***\n\n", fitsfiles[i], fitsfiles[0] );
      return 4;
    }


    maptarget[i] = FunParamGets(fitsconn, "OBJECT", 0, NULL, &got);

    if (!got) {
      maptarget[i] = (char *)calloc(18,sizeof(char));
      strcpy(maptarget[i], "Unknown Target");
    }

    if (suppress == 0) {
      if (strcasecmp(maptarget[0], maptarget[i])) {
	fprintf(stderr,"*** Warning: The target of %s (as specified in the header) does not match that of %s! This may simply be due to different naming conventions between maps, but should be checked. ***\n\nThis message will be shown only once. Any future instances will be suppressed.\n\n", fitsfiles[i], fitsfiles[0]);
	suppress = 1;
      }
    }
    
    units[i] = FunParamGets(fitsconn, "BUNIT", 0, NULL, &got);
    if ( (strcasecmp("JY/BEAM", units[i])) != 0 ) {
      fprintf(stderr,"\n*** Error: The units used in the map %s (as specified in the header) are not JY/BEAM! ***\n\n", fitsfiles[i] );
      return 4;
    }
    
    freq_ok = 0;
    frequencyobs[currentset][i] = FunParamGetd(fitsconn, "REFFREQ", 0, 0.0, &got);

    
    if (!got) {
      frequencyobs[currentset][i] = FunParamGetd(fitsconn, "FREQ", 0, 0.0, &got); // For GLEAM
    }
      
    if (!got) {
      tmp_valcheck[i] = FunParamGets(fitsconn, "CTYPE3", 0, NULL, &got2);

      if (!got2) { // Check we actually have something

	tmp_valcheck[i] = FunParamGets(fitsconn, "VOBS", 0, NULL, &got3); // For Miriad
	  
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
	chomp(tmp_valcheck[i]);
	if ( (strcasecmp("FREQ", tmp_valcheck[i]) == 0) || (strcasecmp("FREQ-LSR", tmp_valcheck[i]) == 0) ) { // LSR is Miriad
	  frequencyobs[currentset][i] = FunParamGetd(fitsconn, "CRVAL3", 0, 0.0, &got);

	  if (!got) {
	    printf("I've managed to find the frequency under CTYPE3, but the accompanying CRVAL3 appears to be missing. Unable to proceed. Please add the values to the header or contact the developer to have your header style added...\n");
	    fprintf(stderr,"\n*** Error: Unable to obtain the frequency of %s! ***\n\n", fitsfiles[i] );
	    return 4;

	  }
	  else {
	    freq_ok = 1;
	  }
	}
	/*if (freq_ok == 0) { // This probably isnt valid...
	  frequencyobs[currentset][i] = FunParamGetd(fitsconn, "ALTRVAL", 0, 0.0, &got);
	  
	  if (got) {
	  freq_ok = 1;
	  }
	  }*/
	if (freq_ok == 0) {
	  printf("I've ran out of ideas where the observed frequency could be located. Please contact the developer to have your header style added...\n");
	  fprintf(stderr,"\n*** Error: Unable to obtain the frequency of %s! ***\n\n", fitsfiles[i] );
	  return 4;  
	}
      }
    }

    frequency[currentset][i] = frequencyobs[currentset][i] * (redshift[currentset] + 1);

    /*
      frequency[currentset][i] = FunParamGetd(fitsconn, "CRVAL3", 0, 0.0, &got2);
      if (!got2) {
      fprintf(stderr,"\n*** Error: Unable to obtain the frequency of %s! ***\n\n", fitsfiles[i] );
      return 4;  
      }
    */

    if (i > 0) {
      for (a=0; a < i; a++) {
	if (frequencyobs[currentset][a] == frequencyobs[currentset][i]) {
	  fprintf(stderr,"\n*** Error: The frequency of the maps %s and %s is the same. This would lead to some odd results! ***\n\n", fitsfiles[a], fitsfiles[i]);
	  return 4; 
	}
      }
    }

    delt2[currentset][i]=FunParamGetd(fitsconn, "CDELT2", 0, 0.0, &got);

    if (!got) {

      delt2[currentset][i]=FunParamGetd(fitsconn, "CD2_2", 0, 0.0, &got2); // GLEAM uses this

      if (!got2) {
	fprintf(stderr,"\n*** Error: Unable to determine the pixel size! (Failed on map frequency %.3e, DELT2) ***\n\n", frequencyobs[currentset][i]);
	return 4;
      }
    }

    cellsize[currentset] = 3600 * delt2[currentset][i]; // Convert

   
    delt1[currentset][i]=FunParamGetd(fitsconn, "CDELT1", 0, 0.0, &got);

    if (!got) {

      delt1[currentset][i]=FunParamGetd(fitsconn, "CD1_1", 0, 0.0, &got2); // GLEAM uses this

      if (!got2) {
	fprintf(stderr,"\n*** Error: Unable to determine the pixel size! (Failed on map frequency %.3e, DELT1) ***\n\n", frequencyobs[currentset][i]);
	return 4;
      }
    }

    if ( ((abs(delt2[currentset][i]) - abs(delt1[currentset][i])) > celltol) || ((abs(delt1[currentset][i]) - abs(delt2[currentset][i])) > celltol) ) {
      fprintf(stderr,"\n*** Error: Pixels aren't square! (Failed on map frequency %.3e) ***\n\n", frequencyobs[currentset][i]);
      return 4;
    }

 
    tmp_delt = delt2[currentset][i] * 3600.0;

    FunInfoGet(fitsconn, FUN_HEADER, &head, 0);
      
    // At least all params below need adjusting to automatically accept to CASA data

    ra[currentset][i] = FunParamGetd(fitsconn, "CRVAL1", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain the RA of %s! ***\n\n", fitsfiles[i] );
      return 4;  
    }

    dec[currentset][i] = FunParamGetd(fitsconn, "CRVAL2", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain the Dec of %s! ***\n\n", fitsfiles[i] );
      return 4;  
    }

    // Adjust the RA and DEC so they are on the centre pixel

    crpix2[currentset][i] = FunParamGetd(fitsconn, "CRPIX2", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain the pixel number for the DEC of %s! ***\n\n", fitsfiles[i] );
      return 4;
    }

    avgdec = (dec[currentset][i] + (dec[currentset][i]+(((tmp_xdim[0]/2) - crpix2[currentset][i]) * ( cellsize[currentset] / 3600)) ) ) / 2 ;

    dec[currentset][i] += ( ( (tmp_xdim[0]/2) - crpix2[currentset][i]) * ( cellsize[currentset] / 3600) );

    //printf("Dec: %.4e\n",dec[currentset][i]);


    crpix1[currentset][i] = FunParamGetd(fitsconn, "CRPIX1", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"\n*** Error: Unable to obtain the pixel number for the RA of %s! ***\n\n", fitsfiles[i] );
      return 4;  
    }

    crota1[currentset][i] = FunParamGetd(fitsconn, "CROTA1", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"*** Warning: Unable to obtain the rotation value for the RA of %s! This is normal for some header types. Assuming a value of 0.0... ***\n", fitsfiles[i] );
      crota1[currentset][i] = 0.0;
      //return 4;
    }
    
    crota2[currentset][i] = FunParamGetd(fitsconn, "CROTA2", 0, 0.0, &got);
    if (!got) {
      fprintf(stderr,"*** Warning: Unable to obtain the rotation value for the DEC of %s! This is normal for some header types. Assuming a value of 0.0... ***\n", fitsfiles[i] );
      crota2[currentset][i] = 0.0;
      //return 4;
    }
  
    eq[currentset][i] = FunParamGetd(fitsconn, "EQUINOX", 0, 0.0, &got);
    if (!got) {
      
      eq[currentset][i] = FunParamGetd(fitsconn, "EPOCH", 0, 0.0, &got2);
	  
      if (!got2) {
	fprintf(stderr,"\n*** Error: Unable to obtain the EQUINOX/EPOCH of %s! ***\n\n", fitsfiles[i] );
	//eq[currentset][i] = 2000;
	return 4;
      }
    }


    float adjcell;

    //printf("Cellsize: %.4f\n", cellsize[currentset]);

    adjcell = (cellsize[currentset]/3600)*cos(degtorad * avgdec);

    ra[currentset][i] += ((crpix1[currentset][i] - (tmp_xdim[0]/2))*adjcell);

      /*
    printf("adjcell: %.4f\n", adjcell);
    
    conv_ra = 240*(360 + ra[currentset][i]); // arcsec

    conv_ra +=  ((wholecirclesecs*cos(degtorad * avgdec) ) / 360) * ((rapix - (tmp_xdim[0]/2)) * (cellsize[currentset]/3600));

    ra[currentset][i] = (conv_ra/240)-360;
      */    

    /*
    float delta_ra, sep;

    conv_dec = degtorad * avgdec; // Convert to radians for the trig functions

    delta_ra = ((tmp_xdim[0]/2)-rapix) * (cellsize[currentset]/3600);

    sep =  delta_ra * (cos(conv_dec));

    ra[currentset][i] += sep;
    */

    //conv_dec = dec[currentset][i] * 3600; // arcsec

    // Convert the RA and Dec to seconds
    //conv_ra = 240*(360 + ra[currentset][i]); 
    
    // Figure out the RA of the centre pixel for the given dec
    //conv_ra += ((wholecirclesecs*cos(degtorad * dec[currentset][i]) ) / 360) * ((rapix - (tmp_xdim[i]/2)) * (cellsize[currentset]/3600));


    // Convert back to decimal format
    //ra[currentset][i] = (conv_ra/240)-360;

    //ra[currentset][i] += ( ( (tmp_xdim[i]/2) - rapix) * (cellsize[currentset] / 3600) );


    // Automatically determine the header type

    str_casadata=FunParamGets(fitsconn, "ORIGIN", 0, NULL, &got);

    //printf("str_casadata: %s\n", str_casadata);

    if (got) {
      if ( strncasecmp(str_casadata, str_compare_aips, 4) == 0 ) {
	tmp_casadata = 0;
      }
      else if ( strncasecmp(str_casadata, str_compare_casa, 4) == 0) {
	tmp_casadata = 1;
      }
      else {
	printf("Unable to automatically determine map origin. Assuming the manually set header type\n\n");
	tmp_casadata = casadata;
      }
    }
    else {
      printf("Unable to automatically determine map origin. Assuming the manually set header type\n\n");
      tmp_casadata = casadata;
    }

    if(tmp_casadata == 0) {

      // This deals with AIPS's crazy header encoding for beam parameters
      for(m=1; m<=head->ncard; m++){
	card = ft_cardnth(head,m);
	if( card && card->c && *card->c){
	  if (!strncmp(card->c, "HISTORY AIPS   CLEAN BMAJ",25)) {
	    sscanf(card->c+26,"%f",&tmp_bmaj[i]);
	    sscanf(card->c+44,"%f",&tmp_bmin[i]);
	    sscanf(card->c+61,"%f",&tmp_bpa[i]);
	    tmp_bmaj[i]*=3600; // Convert to arcseconds
	    tmp_bmin[i]*=3600;
	    tmp_bpa[i] += 90.0; // Convert from AIPS orientation to (0 is Y axis) to plot format (0 is X axis)
	    //printf("BMAJ: %.3f BMIN %.3f BPA %.2f \n", tmp_bmaj[i], tmp_bmin[i], tmp_bpa[i]);
	  }
	}
      }
    }
    else {
      tmp_bmaj[i] = FunParamGetd(fitsconn, "BMAJ", 0, 0.0, &got);
      tmp_bmin[i] = FunParamGetd(fitsconn, "BMIN", 0, 0.0, &got);
      tmp_bpa[i] = FunParamGetd(fitsconn, "BPA", 0, 0.0, &got);
      tmp_bmaj[i]*=3600;
      tmp_bmin[i]*=3600;
      tmp_bpa[i] += 90.0; // Check this for CASA
      //printf("BMAJ: %.3f BMIN %.3f \n", tmp_bmaj[i], tmp_bmin[i]);
    }
    
    if ( (fabs(tmp_bmaj[i]) - fabs(tmp_bmaj[0]) ) > (0.01*fabs(tmp_bmaj[0]) ) || (fabs(tmp_bmin[i]) - fabs(tmp_bmin[0]) ) > (0.01*fabs(tmp_bmin[0]) ) ) {
      fprintf(stderr,"\n*** Error: Beam parameters for %s and %s don't match! ***\n\n", fitsfiles[0], fitsfiles[i]);
      return 4;
    }

    // Only check BPA is there is a >1 per cent fifference between bmin and bmaj. This is to account for the fact CASA sometimes randomises the BPA in the header for a circular beam, even if set to 0.
    if ( (tmp_bmaj[i] - tmp_bmin[i]) > (0.01 * tmp_bmin[i]) ) {
      printf("Beam is not circular (> 1%% difference). Checking the BPA matches to within 1 degree.\n");
      // Allow a 1 degree difference in BPA
      if(fabs(tmp_bpa[i]) > fabs(tmp_bpa[0])) {
	if ( (fabs(tmp_bpa[i]) - fabs(tmp_bpa[0]) ) > 1 ) {
	  fprintf(stderr,"\n*** Error: Beam position angle for %s and %s don't match! ***\n\n", fitsfiles[0], fitsfiles[i]);
	  return 4;
	}
      }
      else {
	  if ( (fabs(tmp_bpa[0]) - fabs(tmp_bpa[i]) ) > 1 ) {
	    fprintf(stderr,"\n*** Error: Beam position angle for %s and %s don't match! ***\n\n", fitsfiles[0], fitsfiles[i]);
	    return 4;
	  }
	}
    }
    
    // Parse the headers for encoding of the CLEAN beam size
    if (tmp_bmaj[i]>0.0) {
      tmp_beamarea[i]=2.0*PI*(tmp_bmaj[i]*tmp_bmin[i])/(GFACTOR*GFACTOR*tmp_delt*tmp_delt);
    } else {
      fprintf(stderr,"\n*** Failed to find CLEAN beam parameters -- Please check you have the correct header type selected (e.g. casadata) and the folder only contains radio maps? ***\n\n");
      return 4;
    }
    
    FunClose(fitsconn);

    // Get the info needed later for the background of this map

    bgconn = FunOpen(bgloc, "r", NULL);
   
    if (bgconn == NULL){
      fprintf(stderr,"\nError: Unable to open %s, please check the file exists and access permissions are correct.\nExiting...\n\n", fitsfiles[i]);
      exit(2);
    }

    if( !(bgbuff[currentset][i] = FunImageGet(bgconn, NULL, "bitpix=-32")) ) {
      fprintf(stderr, "Error: Unable to use ImageGet on %s, please check the file is in the correct format.\nExiting...\n\n", fitsfiles[i]);
    }

    FunClose(bgconn);

  }

  // Set the global dimensions to those found locally 
  xdim[currentset] = tmp_xdim[0];
  ydim[currentset] = tmp_ydim[0];
  bmaj[currentset] = tmp_bmaj[0];
  bmin[currentset] = tmp_bmin[0];
  beamarea[currentset] = tmp_beamarea[0];

  // Output what we found (More to do)

  printf("========================================================================\n");

  if (numdatasets == 0) {
    printf("Set %d was added or changed. There is now %d dataset in total\n", currentset, (numdatasets+1));
  }
  else {
    printf("Set %d was added or changed. There are now %d datasets in total\n", currentset, (numdatasets+1));
  }

  printf("There were %d maps found in the directory %s\n", imgnum[currentset], dirname);
  printf("Pixel size is %.2e arcsec\n", cellsize[currentset]);
  printf("Beam area is %.2f pixels\n", beamarea[currentset]);
  printf("Beam major axis %.2f arcsec, minor axis %.2f arcsec\n", bmaj[currentset],bmin[currentset]);
  printf("The map dimensions are %d by %d\n", xdim[currentset], ydim[currentset]);

  printf("========================================================================\n\n");

  
  printf("Loading map data...\n");

  // Allocate memory to the overall flux, rms arrays and the min max flux values

  flux[currentset] = (float ***)calloc(imgnum[currentset],sizeof(float **)); // Maps
  yaxismax[currentset] = (float *)calloc(imgnum[currentset],sizeof(float));
  yaxismin[currentset] = (float *)calloc(imgnum[currentset],sizeof(float));

  for (a=0; a<imgnum[currentset]; a++) {

    flux[currentset][a] = (float **)calloc(xdim[currentset],sizeof(float *));

    for (i=0; i<xdim[currentset]; i++) {

      flux[currentset][a][i] = (float *)calloc(ydim[currentset],sizeof(float));
    }
  }

  rms[currentset] = (float *)calloc(imgnum[currentset],sizeof(float));

  //xaxismax = -1e+20;
  //xaxismin = 1e+20;

  float averagebackground;
  float rmsnumerator;

  // Find a RMS value for each map

  printf("\nRMS Noise (from background region):\n");

  for (a=0; a < imgnum[currentset]; a++) {
    bgflux = 0;
    count = 0;
    averagebackground = 0;
    rmsnumerator = 0;

    // Calculate the average flux
    for (i=0; i<xdim[currentset]; i++){
      for (j=0; j<ydim[currentset]; j++){
	bgflux += (bgbuff[currentset][a][i*ydim[currentset]+j] / beamarea[currentset]); // In per pixel format

	// If the flux is in the region (allowing for negative values) then increase the count
	if ( ( bgbuff[currentset][a][i*xdim[currentset]+j] < -1e-20) || ( bgbuff[currentset][a][i*xdim[currentset]+j] > 1e-20) ) {
	  count++;
	}
      }
    }

    //printf("Count: %d \n", count);

    averagebackground = bgflux / count;

    // Find the rms
    for (i=0; i<xdim[currentset]; i++){
      for (j=0; j<ydim[currentset]; j++){
	// Ensure we are in the region
	if ( ( bgbuff[currentset][a][i*xdim[currentset]+j] < -1e-20) || ( bgbuff[currentset][a][i*xdim[currentset]+j] > 1e-20) ) {
	  rmsnumerator += pow((bgbuff[currentset][a][i*xdim[currentset]+j] / beamarea[currentset]) - averagebackground, 2);

	}
      }
    }

    rms[currentset][a] =  sqrt(rmsnumerator / count); // This gives RMS in a per pixel format
    
    // Check we dont have a NaN or infinite value (usually caused by the region being outside the edge of the map)
    if ( isnan(rms[currentset][a]) || isinf(rms[currentset][a]) ) {
	fprintf(stderr,"\n*** Error: Unable to detemrine the RMS noise for the map at %.2eHz as it is returning a value of %.4e. Background region may be beyond the extent of the image. ***\n\n", frequencyobs[currentset][a], rms[currentset][a]);
	return 404; 
	  }

        printf("%.2e Hz (Observed), %.2e Hz (Rest) = %.4e Jy/Pixel  (%.4e Jy/Beam)\n", frequencyobs[currentset][a], frequency[currentset][a], rms[currentset][a], rms[currentset][a]*beamarea[currentset]);
  }

  // Create some space so it is clearer
  printf("\n");

  // Loop through by map and find the source flux

  for (a=0; a < imgnum[currentset]; a++) {

    yaxismax[currentset][a] = -1e+20;
    yaxismin[currentset][a] = 1e+20;
 
    // Loop through by pixel
    for (i=(0); i<xdim[currentset]; i++){
 
      for (j=(0); j<ydim[currentset]; j++){
	    
	if ((imgbuff[a][i*ydim[currentset]+j] > (sigmalevel * rms[currentset][a]) ) && (imgbuff[a][i*ydim[currentset]+j] >1e-11) ) { // Quick hack to fix Tom's maps

	    // if (imgbuff[a][i*xdim[currentset]+j] > (sigmalevel * rms[currentset][a]) ) {

	  flux[currentset][a][i][j] = imgbuff[a][i*ydim[currentset]+j] / beamarea[currentset];


	  if (flux[currentset][a][i][j] > yaxismax[currentset][a]) {
	    yaxismax[currentset][a] = flux[currentset][a][i][j];
	  }

	  if (flux[currentset][a][i][j] < yaxismin[currentset][a]) {
	    yaxismin[currentset][a] = flux[currentset][a][i][j];
	  }
	}
	else {
	  flux[currentset][a][i][j] = 0.0;
	}
      }
    }

    // Clean up the image buffer as we go
    free(imgbuff[a]);
  }

  // Sort the order of the maps in terms of frequency

  fluxorder[currentset] =  (int *)calloc((imgnum[currentset]+1),sizeof(int));

  for (i=0; i < imgnum[currentset]; i++) {
    if (i == 0) {
      curminfreq = 1e-32;
    }
    else {
      curminfreq = frequency[currentset][fluxorder[currentset][i-1]];
    }
    curmaxfreq = 1e+32;

    for (a=0; a < imgnum[currentset]; a++) {
      if (a==0) {
	
	fluxorder[currentset][i] = a;
      }
      if ((frequency[currentset][a] > curminfreq) && (frequency[currentset][a] < curmaxfreq)) {
	curmaxfreq = frequency[currentset][a];
	fluxorder[currentset][i] = a;
      }
    }
  }


  // Set the fractional border
  zoomedx =  xdim[currentset] / zoom;
  istart = xshift + ( (xdim[currentset] - (int)zoomedx) / 2);
  iend = xshift + xdim[currentset] - ( (xdim[currentset] - (int)zoomedx) / 2);
  border[currentset] = ceil((iend-istart) * frac_border); // Always round up

  sprintf(output,"/xs");

  sprintf(fluxmaptitle,"Flux of %s at %.2e Hz (%.2f Sigma Detection, Raw Data)", maptarget[fluxorder[currentset][0]], frequencyobs[currentset][fluxorder[currentset][0]], sigmalevel);

  fluxmap(flux, rms, xdim, ydim, fluxorder[currentset][0], yaxismax, yaxismin, currentset, border[currentset], zoom, sigmalevel, fluxmaptitle, titles, labels, axismarks, 0, 1, fluxlogs, fluxcut, usercut, wedgemultiplier, output, ra[currentset][posmap], dec[currentset][posmap], cellsize[currentset], scaletype, bmaj[currentset], bmin[currentset], bpa[currentset], beamposx, beamposy, plotbeam, xshift, yshift, largetxt, 0);

  //Free up the memory
  for (i=0; i<MAXNUMFILES; i++) {
    free(fitsfiles[i]);
    free(bgfiles[i]);
  }

  //Copy the target over to the transfering variable
  strcpy(tmp_maptarget, maptarget[0]);

  free(maptarget);
  free(fileloc);
  free(units);
  free(tmp_valcheck);
  free(tmp_bmaj);
  free(tmp_bmin);
  free(tmp_beamarea);
  free(bgloc);
  free(imgbuff);
  free(tmp_xdim);
  free(tmp_ydim);
  free(fitsfiles);
  free(bgfiles);

  return 0;

}
