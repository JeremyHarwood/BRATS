/*
  Functions for importing a full dataset in .brats format 
*/

int importcoords(char *importfilename, int currentset, float **ra, float **dec, int imgnum) {

  char line_storage[1024], buffer[1024];
  int line_num = 1, foundra = 0, founddec = 0, raline, decline;	

  FILE *importfile = fopen(importfilename, "r");

  ra[currentset] = (float *)calloc(imgnum,sizeof(float));
  dec[currentset] = (float *)calloc(imgnum,sizeof(float));

  printf("Reading the coordinates...\n");

  if (importfile != NULL) {

    while( fgets(line_storage, sizeof(line_storage), importfile) != NULL )  {

      sscanf(line_storage,"%s", buffer);

      if (foundra == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - raline) ) ) {
	  printf("\n*** ERROR: Finished reading the RA before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - raline) ) ) {
	  printf("RA imported successfully\n");
	  foundra = 2;
	}
	else {
	  sscanf(buffer,"%f", &ra[currentset][(line_num - raline) - 1]);
	  // printf("RA of map %d is %f (taken from the string %s)\n", (line_num - raline) - 1,  ra[(line_num - raline) - 1], buffer);
	}
      }

      if (founddec == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - decline) ) ) {
	  printf("\n*** ERROR: Finished reading the DEC before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - decline) ) ) {
	  printf("DEC imported successfully\n");
	  founddec = 2;
	}
	else {
	  sscanf(buffer,"%f", &dec[currentset][(line_num - decline) - 1]);
	  // printf("DEC of map %d is %f (taken from the string %s)\n", (line_num - decline) - 1,  dec[(line_num - decline) - 1], buffer);
	}
      }

      // If we are done, go back to the control script so we don't waste time or hit the EOF
      if ( (foundra == 2) && (founddec == 2) ) {
	return 0;
      }

      if (strcmp(buffer,"###RA###") == 0)  {
	printf("RA found on line %d\n", line_num);
	foundra = 1;
	raline = line_num;
      }

      if (strcmp(buffer,"###DEC###") == 0)  {
	printf("DEC found on line %d\n", line_num);
	founddec = 1;
	decline = line_num;
      }

      if (strcmp(buffer,"###ENDOFIMPORTFILE###") == 0)  {
	printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EFI1001) ***\n\n");
	fclose(importfile);
	return 404;
      }

      line_num++;
    }

    fclose(importfile);

    if ( (foundra == 0) || (founddec == 0) ) {
      printf("\n*** ERROR: Unable to find RA and / or DEC fields! The data may be corrupt ***\n\n");
      return 404;
    }

  }
  else {
	fprintf(stderr,"*** Error: Unable to locate or access the .brats file requested ***\n\n");
	return 404;
      }

  return 0;

}


int importfrequencies(char *importfilename, int currentset, float **frequency, float **frequencyobs, int **fluxorder, int imgnum) {

  char line_storage[1024], buffer[1024];
  int line_num = 1, foundfreq = 0, foundfreqobs = 0, foundfluxorder = 0, freqline, freqobsline, fluxorderline;	

  FILE *importfile = fopen(importfilename, "r");

  frequency[currentset] = (float *)calloc(imgnum,sizeof(float));
  frequencyobs[currentset] = (float *)calloc(imgnum,sizeof(float));
  fluxorder[currentset] = (int *)calloc(imgnum,sizeof(int));


  printf("Reading the map frequencies...\n");

  if (importfile != NULL) {

    while( fgets(line_storage, sizeof(line_storage), importfile) != NULL )  {

      sscanf(line_storage,"%s", buffer);

      if (foundfreq == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - freqline) ) ) {
	  printf("\n*** ERROR: Finished reading FREQUENCY before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - freqline) ) ) {
	  printf("FREQUENCY imported successfully\n");
	  foundfreq = 2;
	}
	else {
	  sscanf(buffer,"%f", &frequency[currentset][(line_num - freqline) - 1]);
	  // printf("FREQUENCY of map %d is %f (taken from the string %s)\n", (line_num - freqline) - 1,  frequency[(line_num - freqline) - 1], buffer);
	}
      }

      if (foundfreqobs == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - freqobsline) ) ) {
	  printf("\n*** ERROR: Finished reading the FREQUENCYOBS before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - freqobsline) ) ) {
	  printf("FREQUENCYOBS imported successfully\n");
	  foundfreqobs = 2;
	}
	else {
	  sscanf(buffer,"%f", &frequencyobs[currentset][(line_num - freqobsline) - 1]);
	  // printf("FREQUENCYOBS of map %d is %f (taken from the string %s)\n", (line_num - freqobsline) - 1,  frequencyobs[currentset][(line_num - freqobsline) - 1], buffer);
	}
      }


      if (foundfluxorder == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - fluxorderline) ) ) {
	  printf("\n*** ERROR: Finished reading the FLUXORDER before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - fluxorderline) ) ) {
	  printf("FLUXORDER imported successfully\n");
	  foundfluxorder = 2;
	}
	else {
	  sscanf(buffer,"%d", &fluxorder[currentset][(line_num - fluxorderline) - 1]);
	  // printf("FLUXORDER of map %d is %d (taken from the string %s)\n", (line_num - fluxorderline) - 1,  fluxorder[currentset][(line_num - fluxorderline) - 1], buffer);
	}
      }


      // If we are done, go back to the control script so we don't waste time or hit the EOF
      if ( (foundfreq == 2) && (foundfreqobs == 2) && (foundfluxorder == 2) ) {
	return 0;
      }


      if (strcmp(buffer,"###FREQUENCY###") == 0)  {
	printf("FREQUENCY found on line %d\n", line_num);
	foundfreq = 1;
	freqline = line_num;
      }

      if (strcmp(buffer,"###FREQUENCYOBS###") == 0)  {
	printf("FREQUENCYOBS found on line %d\n", line_num);
	foundfreqobs = 1;
	freqobsline = line_num;
      }

      if (strcmp(buffer,"###FLUXORDER###") == 0)  {
	printf("FLUXORDER found on line %d\n", line_num);
	foundfluxorder = 1;
	fluxorderline = line_num;
      }
      if (strcmp(buffer,"###ENDOFIMPORTFILE###") == 0)  {
	printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EFI1002) ***\n\n");
	fclose(importfile);
	return 404;
      }

      line_num++;
    }

    fclose(importfile);

    if ( (foundfreq == 0) || (foundfreqobs == 0) || (foundfluxorder == 0)) {
      printf("\n*** ERROR: Unable to find FREQUENCY and / or FREQUENCYOBS fields! The data may be corrupt ***\n\n");
      return 404;
    }

  }
  else {
	fprintf(stderr,"*** Error: Unable to locate or access the .brats file requested ***\n\n");
	return 404;
      }

  return 0;

}



int importrawfluxdata(char *importfilename, int currentset, int imgnum, float **rms, float **fluxcalerror, float **minflux, float **maxflux, float ****flux, float ***bgbuff, int xdim, int ydim) {

  char line_storage[1024], buffer[1024], readstr1[100], readstr2[100], readstr3[100], readstr4[100];
  int line_num = 1, foundrms = 0, foundfluxcalerror = 0, foundminflux = 0, foundmaxflux = 0, foundflux = 0, foundbgbuff = 0, rmsline, fluxcalerrorline, minfluxline, maxfluxline, tmp_a, tmp_x, tmp_y, a, i;

  FILE *importfile = fopen(importfilename, "r");


  // Setup the required memory
  bgbuff[currentset] =(float **)calloc(imgnum,sizeof(float *));
  flux[currentset] = (float ***)calloc(imgnum,sizeof(float **)); // Maps
  minflux[currentset] = (float *)calloc(imgnum,sizeof(float));
  maxflux[currentset] = (float *)calloc(imgnum,sizeof(float));
  rms[currentset] = (float *)calloc(imgnum,sizeof(float));


  for (a=0; a<imgnum; a++) {

    flux[currentset][a] = (float **)calloc(xdim,sizeof(float *));
    bgbuff[currentset][a] =(float *)calloc((xdim*ydim),sizeof(float));

    for (i=0; i<xdim; i++) {
      
      flux[currentset][a][i] = (float *)calloc(ydim,sizeof(float));
    }
  }


  printf("Reading the raw flux data...\n");

  if (importfile != NULL) {

    while( fgets(line_storage, sizeof(line_storage), importfile) != NULL )  {

      sscanf(line_storage,"%s", buffer);

      if (foundrms == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - rmsline) ) ) {
	  printf("\n*** ERROR: Finished reading RMS before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - rmsline) ) ) {
	  printf("RMS imported successfully\n");
	  foundrms = 2;
	}
	else {
	  sscanf(buffer,"%f", &rms[currentset][(line_num - rmsline) - 1]);
	  // printf("RMS of map %d is %f (taken from the string %s)\n", (line_num - rmsline) - 1,  rms[currentset][(line_num - rmsline) - 1], buffer);
	}
      }

      if (foundfluxcalerror == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - fluxcalerrorline) ) ) {
	  printf("\n*** ERROR: Finished reading the FLUXCALERROR before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - fluxcalerrorline) ) ) {
	  printf("FLUXCALERROR imported successfully\n");
	  foundfluxcalerror = 2;
	}
	else {
	  sscanf(buffer,"%f", &fluxcalerror[currentset][(line_num - fluxcalerrorline) - 1]);
	  // printf("FLUXCALERROR of map %d is %f (taken from the string %s)\n", (line_num - fluxcalerrorline) - 1,  fluxcalerror[currentset][(line_num - fluxcalerrorline) - 1], buffer);
	}
      }


      if (foundminflux == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - minfluxline) ) ) {
	  printf("\n*** ERROR: Finished reading the MINFLUX before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - minfluxline) ) ) {
	  printf("MINFLUX imported successfully\n");
	  foundminflux = 2;
	}
	else {
	  sscanf(buffer,"%f", &minflux[currentset][(line_num - minfluxline) - 1]);
	  // printf("MINFLUX of map %d is %e (taken from the string %s)\n", (line_num - minfluxline) - 1,  minflux[currentset][(line_num - minfluxline) - 1], buffer);
	}
      }


      if (foundmaxflux == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - maxfluxline) ) ) {
	  printf("\n*** ERROR: Finished reading the MAXFLUX before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - maxfluxline) ) ) {
	  printf("MAXFLUX imported successfully\n");
	  foundmaxflux = 2;
	}
	else {
	  sscanf(buffer,"%f", &maxflux[currentset][(line_num - maxfluxline) - 1]);
	  // printf("MAXFLUX of map %d is %e (taken from the string %s)\n", (line_num - maxfluxline) - 1,  maxflux[currentset][(line_num - maxfluxline) - 1], buffer);
	}
      }


      if (foundflux == 1) {
	if (strcmp(buffer,"###END###") == 0) {
	printf("FLUX imported successfully\n");
	foundflux = 2;
	}
	else {
	  // Get the map number and XY coordinates
	  sscanf(buffer,"%100[^,],%100[^,],%100[^,]", readstr1, readstr2, readstr3);
	  tmp_a = atoi(readstr1);
	  tmp_x = atoi(readstr2);
	  tmp_y = atoi(readstr3);
	  // printf("buffer: %s tmp_a: %d tmp_x: %d tmp_y: %d \n", buffer, tmp_a, tmp_x, tmp_y);

	  // Retrieve the flux and add it to the array
	  sscanf(buffer,"%*100[^,],%*100[^,],%*100[^,],%s", readstr4);

	  flux[currentset][tmp_a][tmp_x][tmp_y] = atof(readstr4);

	  //printf("FLUX of map %d at X:%d Y:%d is %e (taken from the string %s)\n", tmp_a, tmp_x, tmp_y, flux[currentset][tmp_a][tmp_x][tmp_y], buffer);

	}
      }


      if (foundbgbuff == 1) {
	if (strcmp(buffer,"###END###") == 0) {
	  printf("BGBUFF imported successfully\n");
	  foundbgbuff = 2;
	}
	else {
	  // Get the map number and XY coordinates
	  sscanf(buffer,"%100[^,],%100[^,],%100[^,]", readstr1, readstr2, readstr3);
	  tmp_a = atoi(readstr1);
	  tmp_x = atoi(readstr2);
	  tmp_y = atoi(readstr3);

	  // Retrieve the bgbuff and add it to the array
	  sscanf(buffer,"%*100[^,],%*100[^,],%*100[^,],%s", readstr4);

	  bgbuff[currentset][tmp_a][tmp_x*xdim+tmp_y] = atof(readstr4);

	  // printf("BGBUFF of map %d at X:%d Y:%d is %e (taken from the string %s)\n", tmp_a, tmp_x, tmp_y, bgbuff[currentset][tmp_a][tmp_x*xdim+tmp_y], buffer);
	}
      }


      // If we are done, go back to the control script so we don't waste time or hit the EOF
      if ( (foundrms == 2) && (foundfluxcalerror == 2) && (foundminflux == 2) && (foundmaxflux == 2) && (foundflux == 2) && (foundbgbuff == 2) ) {
	return 0;
      }

      if (strcmp(buffer,"###RMS###") == 0)  {
	printf("RMS found on line %d\n", line_num);
	foundrms = 1;
	rmsline = line_num;
      }

      if (strcmp(buffer,"###FLUXCALERROR###") == 0)  {
	printf("FLUXORDER found on line %d\n", line_num);
	foundfluxcalerror = 1;
	fluxcalerrorline = line_num;
      }

      if (strcmp(buffer,"###MINFLUX###") == 0)  {
	printf("MINFLUX found on line %d\n", line_num);
	foundminflux = 1;
	minfluxline = line_num;
      }

      if (strcmp(buffer,"###MAXFLUX###") == 0)  {
	printf("MAXFLUX found on line %d\n", line_num);
	foundmaxflux = 1;
	maxfluxline = line_num;
      }

      if (strcmp(buffer,"###FLUX###") == 0)  {
	printf("FLUX found on line %d\n", line_num);
	foundflux = 1;
      }
      
      if (strcmp(buffer,"###BGBUFF###") == 0)  {
	printf("BGBUFF found on line %d\n", line_num);
	foundbgbuff = 1;
      }
      
      if (strcmp(buffer,"###ENDOFIMPORTFILE###") == 0)  {
	printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EFI1003) ***\n\n");
	fclose(importfile);
	return 404;
      }

      line_num++;
    }

    fclose(importfile);

    if ( (foundrms == 0) || (foundfluxcalerror == 0) || (foundminflux == 0) || (foundmaxflux == 0) || (foundflux == 0) || (foundbgbuff == 0) ) {
      printf("\n*** ERROR: Unable to find all of the required raw flux fields! The data may be corrupt ***\n\n");
      return 404;
    }

  }
  else {
	fprintf(stderr,"*** Error: Unable to locate or access the .brats file requested ***\n\n");
	return 404;
      }

  return 0;

}


int importnumberofregions(char *importfilename, int currentset, int *regnumber) {

  char line_storage[1024], buffer[1024];
  int line_num = 1, foundregnumber = 0, regnumberline;

  FILE *importfile = fopen(importfilename, "r");

  if (importfile != NULL) {

    while( fgets(line_storage, sizeof(line_storage), importfile) != NULL )  {

      sscanf(line_storage,"%s", buffer);

      if (foundregnumber == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((line_num - regnumberline) != 2) ) {
	  printf("\n*** ERROR: Finished reading REGNUMBER before the expected number of entries had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((line_num - regnumberline) == 2) ) {
	  printf("REGNUMBER imported successfully\n");
	  foundregnumber = 2;
	  fclose(importfile);
	  return 0; // Kick out of the loop start importing the data
	}
	else {
	  sscanf(buffer,"%d", &regnumber[currentset]);
	  printf("REGNUMBER is %d (taken from the string %s)\n", regnumber[currentset], buffer);
	}
      }
  
      if (strcmp(buffer,"###REGNUMBER###") == 0)  {
	printf("REGNUMBER found on line %d\n", line_num);
	foundregnumber = 1;
	regnumberline = line_num;
      }
    
        
      if (strcmp(buffer,"###ENDOFIMPORTFILE###") == 0)  {
	printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EFI1004) ***\n\n");
	fclose(importfile);
	return 404;
      }
      
      line_num++;
    }

    if (foundregnumber == 0) {
      printf("\n*** ERROR: Unable to find the number of regions set! The data may be corrupt ***\n\n");
      return 404;
    }

  }
  else {
    fprintf(stderr,"*** Error: Unable to locate or access the .brats file requested ***\n\n");
    return 404;
  }
 

  return 0;
}

int importregiondata(char *importfilename, int importcompression, int currentset, int imgnum, int xdim, int ydim, int ***regionarray, float ***regflux, int *regnumber, float ***fluxerror, int *setaveraged, float **regionlocx, float **regionlocy, int **regionsize, float **regmaxflux, float **regminflux) {

  char line_storage[1024], buffer[1024], readstr1[100], readstr2[100], readstr3[100];
  int line_num = 1, foundregionarray = 0, foundregflux = 0, foundaveraged = 0, foundregionsize = 0, foundregminmax = 0, foundfluxerror = 0,regionarrayline, regfluxline, averagedline, regionsizeline, regminmaxline, fluxerrorline, tmp_x, tmp_y, a, i, j,regionarrayblanked = 0, totalregionpixels=0;

  FILE *importfile = fopen(importfilename, "r");


  // Setup the required memory
  regflux[currentset] = (float **)calloc(imgnum, sizeof(float *));
  regionarray[currentset] = (int **)calloc(xdim, sizeof(int *));
  fluxerror[currentset] = (float **)calloc(imgnum, sizeof(float *));
  regionlocx[currentset] = (float *)calloc((xdim*ydim)+1, sizeof(float));
  regionlocy[currentset] = (float *)calloc((ydim*xdim)+1, sizeof(float));
  regionsize[currentset] = (int *)calloc((xdim*ydim)+1, sizeof(int));
  regmaxflux[currentset] = (float *)calloc(imgnum, sizeof(float));
  regminflux[currentset] = (float *)calloc(imgnum, sizeof(float));

  for (a=0; a<imgnum; a++) {
    regflux[currentset][a] = (float *)calloc(xdim*ydim, sizeof(float));
    fluxerror[currentset][a] = (float *)calloc(xdim*ydim, sizeof(float));
  }

  for (i=0; i<xdim; i++) {
    regionarray[currentset][i] = (int *)calloc(ydim, sizeof(int));
  }

  printf("Reading the region data...\n");

  if (importfile != NULL) {

    while( fgets(line_storage, sizeof(line_storage), importfile) != NULL )  {

      sscanf(line_storage,"%s", buffer);


      if (foundregflux == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && (((imgnum*regnumber[currentset])+1) != (line_num - regfluxline) ) ) {
	  printf("\n*** ERROR: Finished reading the REGFLUX before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && (((imgnum*regnumber[currentset])+1) == (line_num - regfluxline) ) ) {
	  printf("REGFLUX imported successfully\n");
	  foundregflux = 2;
	}
	else {

	  sscanf(buffer,"%100[^,],%s", readstr1, readstr2);
	    
	  tmp_x = atoi(readstr1);

	  regflux[currentset][tmp_x][( (line_num - regfluxline) - (tmp_x*regnumber[currentset]) )] = atof(readstr2);

	  // printf("REGFLUX of map %d, region %d is %f (taken from the string %s)\n", tmp_x, ( (line_num - regfluxline) - (tmp_x*regnumber[currentset]) ), regflux[currentset][tmp_x][( (line_num - regfluxline) - (tmp_x*regnumber[currentset]) )], buffer);
	}
      }


      if (foundregminmax == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - regminmaxline) ) ) {
	  printf("\n*** ERROR: Finished reading the REGMINMAXFLUX before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	  else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - regminmaxline) ) ) {
	    printf("REGMINMAXFLUX imported successfully\n");
	    foundregminmax = 2;
	  }
	  else {
	    // Get the values from the compressed format and add them to the arrays 
	    sscanf(buffer,"%100[^,],%s", readstr1, readstr2);
	    
	    regminflux[currentset][(line_num - regminmaxline)-1] = atof(readstr1);
	    regmaxflux[currentset][(line_num - regminmaxline)-1] = atof(readstr2);

	    // printf("REGMINMAXFLUX for map %d is min: %e max %e (taken from the string %s)\n",(line_num - regminmaxline)-1, regminflux[currentset][(line_num - regminmaxline)-1], regmaxflux[currentset][(line_num - regminmaxline)-1], buffer);
	    
	  }
	
      }


      if (foundregionsize == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber[currentset]+1) != (line_num - regionsizeline) ) ) {
	  printf("\n*** ERROR: Finished reading the REGSIZE before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber[currentset]+1) == (line_num - regionsizeline) ) ) {
	  printf("%d regions over %d pixels\n", regnumber[currentset], totalregionpixels);
	  printf("REGSIZE imported successfully\n");
	  foundregionsize = 2;
	}
	else {
	  sscanf(buffer,"%d", &regionsize[currentset][(line_num - regionsizeline)]);
	  totalregionpixels+=regionsize[currentset][(line_num - regionsizeline)];
	  // printf("REGSIZE of region %d is %d (taken from the string %s)\n", (line_num - regionsizeline), regionsize[currentset][(line_num - regionsizeline)], buffer);
	}
      }


      if (foundaveraged == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((line_num - averagedline) != 2) ) {
	  printf("\n*** ERROR: Finished reading SETAVERAGED before the expected number of entries had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((line_num - averagedline) == 2) ) {
	  printf("SETAVERAGED imported successfully\n");
	  foundaveraged = 2;
	}
	else {
	  sscanf(buffer,"%d", &setaveraged[currentset]);
	  // printf("SETAVERAGED is %d (taken from the string %s)\n", setaveraged[currentset], buffer);
	}
      }
  

      if (foundregionarray == 1) {
	
	// If the data is in compact format, set everything to -1 (below flux limit) before we begin
	if ( (regionarrayblanked == 0) && (importcompression = 1) ) {
	  
	  for (i=0; i<xdim; i++) {
	    for (j=0; j<ydim; j++) {
	      regionarray[currentset][i][j] = -1;
	    }
	  }
	  regionarrayblanked = 1;
	}

	if (importcompression == 1) {

	  if ( (strcmp(buffer,"###END###") == 0) && ((totalregionpixels+1) != (line_num - regionarrayline) ) ) {
	    printf("\n*** ERROR: Finished reading the REGIONARRAY before the expected number of regions had been reached (compact format)! The data may be corrupt ***\n\n");
	    fclose(importfile);
	    return 404;
	  }
	  else if ( (strcmp(buffer,"###END###") == 0) && ((totalregionpixels+1) == (line_num - regionarrayline) ) ) {
	  printf("REGIONARRAY imported successfully\n");
	  foundregionarray = 2;
	  }
	  else {
	    // Get the values from the compressed format and add them to the arrays 
	    sscanf(buffer,"%100[^,],%100[^,],%s", readstr1, readstr2, readstr3);

	    tmp_x = atoi(readstr1);
	    tmp_y = atoi(readstr2);

	    regionarray[currentset][tmp_x][tmp_y] = atoi(readstr3);

	    regionlocx[currentset][regionarray[currentset][tmp_x][tmp_y]] = tmp_x;
	    regionlocy[currentset][regionarray[currentset][tmp_x][tmp_y]] = tmp_y;


	    // printf("REGIONARRAY at X:%d Y:%d is %d (regionloc x: %f y: %f (taken from the string %s)\n", tmp_x, tmp_y, regionarray[currentset][tmp_x][tmp_y], regionlocx[currentset][regionarray[currentset][tmp_x][tmp_y]], regionlocy[currentset][regionarray[currentset][tmp_x][tmp_y]], buffer);
	    
	  }
	}
	else if (importcompression == 2) {
	  
	  if ( (strcmp(buffer,"###END###") == 0) && ( ((xdim*ydim)+1) != (line_num - regionarrayline) ) ) {
	    printf("\n*** ERROR: Finished reading the REGIONARRAY before the expected number of regions had been reached (full format)! The data may be corrupt ***\n\n");
	    fclose(importfile);
	    return 404;
	  }
	  else if ( (strcmp(buffer,"###END###") == 0) && ( ((xdim*ydim)+1) == (line_num - regionarrayline) ) ) {
	    printf("REGIONARRAY imported successfully\n");
	    foundregionarray = 2;
	  }
	  else {
	    // If we are uncompressed, loop through the entire map size and add the values to the arrays
	    sscanf(buffer,"%100[^,],%100[^,],%s", readstr1, readstr2, readstr3);

	    tmp_x = atoi(readstr1);
	    tmp_y = atoi(readstr2);
	    
	    regionarray[currentset][i][j] = atoi(readstr3);
		
	    // If the values if >0, we have a valid region so add it to the regionloc arrays
	    if (regionarray[currentset][i][j] > 0) {
	      regionlocx[currentset][regionarray[currentset][tmp_x][tmp_y]] = tmp_x;
	      regionlocy[currentset][regionarray[currentset][tmp_x][tmp_y]] = tmp_y;
	    }
	    
	  }
	}
	else {
	  printf("\n*** ERROR: Unknown compression type! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
      }


      if (foundfluxerror == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && (((imgnum*regnumber[currentset])+1) != (line_num - fluxerrorline) ) ) {
	  printf("\n*** ERROR: Finished reading the FLUXERROR before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && (((imgnum*regnumber[currentset])+1) == (line_num - fluxerrorline) ) ) {
	  printf("FLUXERROR imported successfully\n");
	  foundfluxerror = 2;
	}
	else {

	  sscanf(buffer,"%100[^,],%100[^,],%s", readstr1, readstr2, readstr3);
	    
	  tmp_x = atoi(readstr1);
	  tmp_y = atoi(readstr2);

	  fluxerror[currentset][tmp_x][( (line_num - fluxerrorline) - (tmp_x*regnumber[currentset]) )] = atof(readstr3);

	  // printf("FLUXERROR of map %d, region %d is %e (taken from the string %s)\n", tmp_x, ( (line_num - fluxerrorline) - (tmp_x*regnumber[currentset]) ), fluxerror[currentset][tmp_x][( (line_num - fluxerrorline) - (tmp_x*regnumber[currentset]) )], buffer);
	}
      }



      // If we are done, go back to the control script so we don't waste time or hit the EOF
      if ( (foundregionarray == 2) && (foundregflux == 2) && (foundaveraged == 2) && (foundregionsize == 2) && (foundregminmax == 2) && (foundregflux == 2) && (foundfluxerror == 2) ) {
	return 0;
      }


      if (strcmp(buffer,"###REGIONARRAY###") == 0)  {
	printf("REGIONARRAY found on line %d\n", line_num);
	foundregionarray = 1;
	regionarrayline = line_num;
      }

      if (strcmp(buffer,"###SETAVERAGED###") == 0)  {
	printf("SETAVERAGED found on line %d\n", line_num);
	foundaveraged = 1;
	averagedline = line_num;
      }

      if (strcmp(buffer,"###REGSIZE###") == 0)  {
	printf("REGSIZE found on line %d\n", line_num);
	foundregionsize = 1;
	regionsizeline = line_num;
      }
      
      if (strcmp(buffer,"###REGMINMAXFLUX###") == 0)  {
	printf("REGMINMAXFLUX found on line %d\n", line_num);
	foundregminmax = 1;
	regminmaxline = line_num;
      }

      if (strcmp(buffer,"###REGFLUX###") == 0)  {
	printf("REGFLUX found on line %d\n", line_num);
	foundregflux = 1;
	regfluxline = line_num;
      }

      if (strcmp(buffer,"###FLUXERROR###") == 0)  {
	printf("FLUXERROR found on line %d\n", line_num);
	foundfluxerror = 1;
	fluxerrorline = line_num;
      }
      
      if (strcmp(buffer,"###ENDOFIMPORTFILE###") == 0)  {
	printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EFI1005) ***\n\n");
	fclose(importfile);
	return 404;
      }

      line_num++;
    }
  
    fclose(importfile);

    if ( (foundregionarray == 0) || (foundregflux == 0) || (foundaveraged == 0) || (foundregionsize == 0) || (foundregminmax == 0) || (foundfluxerror == 0) ) {
      printf("\n*** ERROR: Unable to find all of the required region fields! The data may be corrupt ***\n\n");
      return 404;
    }
    
  }
  else {
    fprintf(stderr,"*** Error: Unable to locate or access the .brats file requested ***\n\n");
    return 404;
  }

  return 0;

}




int importjpmodel(char *importfilename, int currentset, int imgnum, int xdim, int ydim, int regnumber, float **jpchisquared, float **jpbestage, float **jpmodelalpha, float ***jpmodelflux, float **jpbestnorm, float **jpageerrorsplus, float **jpageerrorsminus, int setmodelres) {

  char line_storage[1024], buffer[1024], readstr1[100], readstr2[100], readstr3[100];
  int line_num = 1, foundbestage = 0, foundbestnorm = 0, foundchisquared = 0, foundmodelalpha = 0, foundmodelflux = 0, foundageerrorsplus = 0, foundageerrorsminus = 0, bestageline, bestnormline, chisquaredline, modelalphaline, modelfluxline, ageerrorsplusline, ageerrorsminusline, tmp_regno, tmp_modresno;

  FILE *importfile = fopen(importfilename, "r");

  printf("Reading the JP model fitting data...\n");

  if (importfile != NULL) {

    while( fgets(line_storage, sizeof(line_storage), importfile) != NULL )  {

      sscanf(line_storage,"%s", buffer);


      if (foundbestage == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - bestageline) ) ) {
	  printf("\n*** ERROR: Finished reading the JPBESTAGE before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - bestageline) ) ) {
	  printf("JPBESTAGE imported successfully\n");
	  foundbestage = 2;
	}
	else {
	  sscanf(buffer,"%f", &jpbestage[currentset][(line_num - bestageline)]);
	  // printf("JPBESTAGE of region %d is %f (taken from the string %s)\n", (line_num - bestageline), jpbestage[currentset][(line_num - bestageline)], buffer);
	}

      }
    

      if (foundbestnorm == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - bestnormline) ) ) {
	  printf("\n*** ERROR: Finished reading the JPBESTNORM before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - bestnormline) ) ) {
	  printf("JPBESTNORM imported successfully\n");
	  foundbestnorm = 2;
	}
	else {
	  sscanf(buffer,"%f", &jpbestnorm[currentset][(line_num - bestnormline)]);
	  // printf("JPBESTNORM of region %d is %.6e (taken from the string %s)\n", (line_num - bestnormline), jpbestnorm[currentset][(line_num - bestnormline)], buffer);
	}

      }
    

      if (foundchisquared == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - chisquaredline) ) ) {
	  printf("\n*** ERROR: Finished reading the JPCHISQUARED before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - chisquaredline) ) ) {
	  printf("JPCHISQUARED imported successfully\n");
	  foundchisquared = 2;
	}
	else {
	  sscanf(buffer,"%f", &jpchisquared[currentset][(line_num - chisquaredline)]);
	  // printf("JPCHISQUARED of region %d is %.6e (taken from the string %s)\n", (line_num - chisquaredline), jpchisquared[currentset][(line_num - chisquaredline)], buffer);
	}

      }


      if (foundmodelalpha == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - modelalphaline) ) ) {
	  printf("\n*** ERROR: Finished reading the JPMODELALPHA before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - modelalphaline) ) ) {
	  printf("JPMODELALPHA imported successfully\n");
	  foundmodelalpha = 2;
	}
	else {
	  sscanf(buffer,"%f", &jpmodelalpha[currentset][(line_num - modelalphaline)]);
	  // printf("JPMODELALPHA of region %d is %.6e (taken from the string %s)\n", (line_num - modelalphaline), jpmodelalpha[currentset][(line_num - modelalphaline)], buffer);
	}

      }


      if (foundageerrorsplus == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - ageerrorsplusline) ) ) {
	  printf("\n*** ERROR: Finished reading the JPAGEERRORSPLUS before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - ageerrorsplusline) ) ) {
	  printf("JPAGEERRORSPLUS imported successfully\n");
	  foundageerrorsplus = 2;
	}
	else {
	  sscanf(buffer,"%f", &jpageerrorsplus[currentset][(line_num - ageerrorsplusline)]);
	  // printf("JPAGEERRORSPLUS of region %d is %.6e (taken from the string %s)\n", (line_num - ageerrorsplusline), jpageerrorsplus[currentset][(line_num - ageerrorsplusline)], buffer);
	}

      }


      if (foundageerrorsminus == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - ageerrorsminusline) ) ) {
	  printf("\n*** ERROR: Finished reading the JPAGEERRORSMINUS before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - ageerrorsminusline) ) ) {
	  printf("JPAGEERRORSMINUS imported successfully\n");
	  foundageerrorsminus = 2;
	}
	else {
	  sscanf(buffer,"%f", &jpageerrorsminus[currentset][(line_num - ageerrorsminusline)]);
	  // printf("JPAGEERRORSMINUS of region %d is %.6e (taken from the string %s)\n", (line_num - ageerrorsminusline), jpageerrorsminus[currentset][(line_num - ageerrorsminusline)], buffer);
	}

      }


      if (foundmodelflux == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ( ((regnumber*setmodelres)+regnumber+1) != (line_num - modelfluxline) ) ) {
	  printf("\n*** ERROR: Finished reading the MODELFLUX before the expected number of regions had been reached (full format)! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ( ((regnumber*setmodelres)+regnumber+1) == (line_num - modelfluxline) ) ) {
	  printf("MODELFLUX imported successfully\n");
	  foundmodelflux = 2;
	}
	else {

	  sscanf(buffer,"%100[^,],%100[^,],%s", readstr1, readstr2, readstr3);
	  
	  tmp_regno = atoi(readstr1);
	  tmp_modresno = atoi(readstr2);

	  // (line_num - modelfluxline) - ((tmp_regno-1)*setmodelres)

	    
	  jpmodelflux[currentset][tmp_regno][tmp_modresno] = atof(readstr3);
	  
	  // printf("MODELFLUX of region %d is %.6e for res point %d (taken from the string %s)\n", tmp_regno,  jpmodelflux[currentset][tmp_regno][tmp_modresno], tmp_modresno, buffer);

	}

      }


      // If we are done, go back to the control script so we don't waste time or hit the EOF
      if ( (foundbestage == 2) && (foundbestnorm == 2) && (foundchisquared == 2) && (foundmodelalpha == 2) && (foundmodelflux == 2) && (foundageerrorsplus == 2) && (foundageerrorsminus == 2) ) {
	return 0;
      }


      if (strcmp(buffer,"###JPBESTAGE###") == 0)  {
	printf("JPBESTAGE found on line %d\n", line_num);
	foundbestage = 1;
	bestageline = line_num;
      }

      if (strcmp(buffer,"###JPBESTNORM###") == 0)  {
	printf("JPBESTNORM found on line %d\n", line_num);
	foundbestnorm = 1;
	bestnormline = line_num;
      }

      if (strcmp(buffer,"###JPCHISQUARED###") == 0)  {
	printf("JPCHISQUARED found on line %d\n", line_num);
	foundchisquared = 1;
	chisquaredline = line_num;
      }

      if (strcmp(buffer,"###JPMODELALPHA###") == 0)  {
	printf("JPMODELALPHA found on line %d\n", line_num);
	foundmodelalpha = 1;
	modelalphaline = line_num;
      }
      
      if (strcmp(buffer,"###JPMODELFLUX###") == 0)  {
	printf("JPMODELFLUX found on line %d\n", line_num);
	foundmodelflux = 1;
	modelfluxline = line_num;
      }

      if (strcmp(buffer,"###JPAGEERRORSPLUS###") == 0)  {
	printf("JPAGEERRORSPLUS found on line %d\n", line_num);
	foundageerrorsplus = 1;
	ageerrorsplusline = line_num;
      }

      if (strcmp(buffer,"###JPAGEERRORSMINUS###") == 0)  {
	printf("JPAGEERRORSMINUS found on line %d\n", line_num);
	foundageerrorsminus = 1;
	ageerrorsminusline = line_num;
      }

      
      if (strcmp(buffer,"###ENDOFIMPORTFILE###") == 0)  {
	printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EFI1012) ***\n\n");
	fclose(importfile);
	return 404;
      }

      line_num++;
    }
  
    fclose(importfile);



    if ( (foundbestage == 0) || (foundbestnorm == 0) || (foundchisquared == 0) || (foundmodelalpha == 0) || (foundmodelflux == 0) || (foundageerrorsplus == 0) || (foundageerrorsminus == 0) ) {
      printf("\n*** ERROR: Unable to find all of the required region fields! The data may be corrupt ***\n\n");
      return 404;
    }
    

  }
  else {
    fprintf(stderr,"*** Error: Unable to locate or access the .brats file requested ***\n\n");
    return 404;
  }

  return 0;

}




int importkpmodel(char *importfilename, int currentset, int imgnum, int xdim, int ydim, int regnumber, float **kpchisquared, float **kpbestage, float **kpmodelalpha, float ***kpmodelflux, float **kpbestnorm, float **kpageerrorsplus, float **kpageerrorsminus, int setmodelres) {

  char line_storage[1024], buffer[1024], readstr1[100], readstr2[100], readstr3[100];
  int line_num = 1, foundbestage = 0, foundbestnorm = 0, foundchisquared = 0, foundmodelalpha = 0, foundmodelflux = 0, foundageerrorsplus = 0, foundageerrorsminus = 0, bestageline, bestnormline, chisquaredline, modelalphaline, modelfluxline, ageerrorsplusline, ageerrorsminusline, tmp_regno, tmp_modresno;

  FILE *importfile = fopen(importfilename, "r");

  printf("Reading the KP model fitting data...\n");

  if (importfile != NULL) {

    while( fgets(line_storage, sizeof(line_storage), importfile) != NULL )  {

      sscanf(line_storage,"%s", buffer);

      if (foundbestage == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - bestageline) ) ) {
	  printf("\n*** ERROR: Finished reading the KPBESTAGE before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - bestageline) ) ) {
	  printf("KPBESTAGE imported successfully\n");
	  foundbestage = 2;
	}
	else {
	  sscanf(buffer,"%f", &kpbestage[currentset][(line_num - bestageline)]);
	  // printf("KPBESTAGE of region %d is %f (taken from the string %s)\n", (line_num - bestageline), kpbestage[currentset][(line_num - bestageline)], buffer);
	}

      }
    
      if (foundbestnorm == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - bestnormline) ) ) {
	  printf("\n*** ERROR: Finished reading the KPBESTNORM before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - bestnormline) ) ) {
	  printf("KPBESTNORM imported successfully\n");
	  foundbestnorm = 2;
	}
	else {
	  sscanf(buffer,"%f", &kpbestnorm[currentset][(line_num - bestnormline)]);
	  // printf("KPBESTNORM of region %d is %.6e (taken from the string %s)\n", (line_num - bestnormline), kpbestnorm[currentset][(line_num - bestnormline)], buffer);
	}
      }
    

      if (foundchisquared == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - chisquaredline) ) ) {
	  printf("\n*** ERROR: Finished reading the KPCHISQUARED before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - chisquaredline) ) ) {
	  printf("KPCHISQUARED imported successfully\n");
	  foundchisquared = 2;
	}
	else {
	  sscanf(buffer,"%f", &kpchisquared[currentset][(line_num - chisquaredline)]);
	  // printf("KPCHISQUARED of region %d is %.6e (taken from the string %s)\n", (line_num - chisquaredline), kpchisquared[currentset][(line_num - chisquaredline)], buffer);
	}

      }


      if (foundmodelalpha == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - modelalphaline) ) ) {
	  printf("\n*** ERROR: Finished reading the KPMODELALPHA before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - modelalphaline) ) ) {
	  printf("KPMODELALPHA imported successfully\n");
	  foundmodelalpha = 2;
	}
	else {
	  sscanf(buffer,"%f", &kpmodelalpha[currentset][(line_num - modelalphaline)]);
	  // printf("KPMODELALPHA of region %d is %.6e (taken from the string %s)\n", (line_num - modelalphaline), kpmodelalpha[currentset][(line_num - modelalphaline)], buffer);
	}
      }


      if (foundageerrorsplus == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - ageerrorsplusline) ) ) {
	  printf("\n*** ERROR: Finished reading the KPAGEERRORSPLUS before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - ageerrorsplusline) ) ) {
	  printf("KPAGEERRORSPLUS imported successfully\n");
	  foundageerrorsplus = 2;
	}
	else {
	  sscanf(buffer,"%f", &kpageerrorsplus[currentset][(line_num - ageerrorsplusline)]);
	  // printf("KPAGEERRORSPLUS of region %d is %.6e (taken from the string %s)\n", (line_num - ageerrorsplusline), kpageerrorsplus[currentset][(line_num - ageerrorsplusline)], buffer);
	}
      }


      if (foundageerrorsminus == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - ageerrorsminusline) ) ) {
	  printf("\n*** ERROR: Finished reading the KPAGEERRORSMINUS before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - ageerrorsminusline) ) ) {
	  printf("KPAGEERRORSMINUS imported successfully\n");
	  foundageerrorsminus = 2;
	}
	else {
	  sscanf(buffer,"%f", &kpageerrorsminus[currentset][(line_num - ageerrorsminusline)]);
	  // printf("KPAGEERRORSMINUS of region %d is %.6e (taken from the string %s)\n", (line_num - ageerrorsminusline), kpageerrorsminus[currentset][(line_num - ageerrorsminusline)], buffer);
	}
      }


      if (foundmodelflux == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ( ((regnumber*setmodelres)+regnumber+1) != (line_num - modelfluxline) ) ) {
	  printf("\n*** ERROR: Finished reading the MODELFLUX before the expected number of regions had been reached (full format)! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ( ((regnumber*setmodelres)+regnumber+1) == (line_num - modelfluxline) ) ) {
	  printf("MODELFLUX imported successfully\n");
	  foundmodelflux = 2;
	}
	else {

	  sscanf(buffer,"%100[^,],%100[^,],%s", readstr1, readstr2, readstr3);
	  
	  tmp_regno = atoi(readstr1);
	  tmp_modresno = atoi(readstr2);

	  // (line_num - modelfluxline) - ((tmp_regno-1)*setmodelres)

	  kpmodelflux[currentset][tmp_regno][tmp_modresno] = atof(readstr3);
	  
	  // printf("MODELFLUX of region %d is %.6e for res point %d (taken from the string %s)\n", tmp_regno,  kpmodelflux[currentset][tmp_regno][tmp_modresno], tmp_modresno, buffer);

	}

      }


      // If we are done, go back to the control script so we don't waste time or hit the EOF
      if ( (foundbestage == 2) && (foundbestnorm == 2) && (foundchisquared == 2) && (foundmodelalpha == 2) && (foundmodelflux == 2) && (foundageerrorsplus == 2) && (foundageerrorsminus == 2) ) {
	return 0;
      }


      if (strcmp(buffer,"###KPBESTAGE###") == 0)  {
	printf("KPBESTAGE found on line %d\n", line_num);
	foundbestage = 1;
	bestageline = line_num;
      }

      if (strcmp(buffer,"###KPBESTNORM###") == 0)  {
	printf("KPBESTNORM found on line %d\n", line_num);
	foundbestnorm = 1;
	bestnormline = line_num;
      }

      if (strcmp(buffer,"###KPCHISQUARED###") == 0)  {
	printf("KPCHISQUARED found on line %d\n", line_num);
	foundchisquared = 1;
	chisquaredline = line_num;
      }

      if (strcmp(buffer,"###KPMODELALPHA###") == 0)  {
	printf("KPMODELALPHA found on line %d\n", line_num);
	foundmodelalpha = 1;
	modelalphaline = line_num;
      }
      
      if (strcmp(buffer,"###KPMODELFLUX###") == 0)  {
	printf("KPMODELFLUX found on line %d\n", line_num);
	foundmodelflux = 1;
	modelfluxline = line_num;
      }

      if (strcmp(buffer,"###KPAGEERRORSPLUS###") == 0)  {
	printf("KPAGEERRORSPLUS found on line %d\n", line_num);
	foundageerrorsplus = 1;
	ageerrorsplusline = line_num;
      }

      if (strcmp(buffer,"###KPAGEERRORSMINUS###") == 0)  {
	printf("KPAGEERRORSMINUS found on line %d\n", line_num);
	foundageerrorsminus = 1;
	ageerrorsminusline = line_num;
      }

      
      if (strcmp(buffer,"###ENDOFIMPORTFILE###") == 0)  {
	printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EFI1013) ***\n\n");
	fclose(importfile);
	return 404;
      }

      line_num++;
    }
  
    fclose(importfile);



    if ( (foundbestage == 0) || (foundbestnorm == 0) || (foundchisquared == 0) || (foundmodelalpha == 0) || (foundmodelflux == 0) || (foundageerrorsplus == 0) || (foundageerrorsminus == 0) ) {
      printf("\n*** ERROR: Unable to find all of the required region fields! The data may be corrupt ***\n\n");
      return 404;
    }
    

  }
  else {
    fprintf(stderr,"*** Error: Unable to locate or access the .brats file requested ***\n\n");
    return 404;
  }

  return 0;

}







int importjptribmodel(char *importfilename, int currentset, int imgnum, int xdim, int ydim, int regnumber, float **jptribchisquared, float **jptribbestage, float **jptribmodelalpha, float ***jptribmodelflux, float **jptribbestnorm, float **jptribbleageerrorsplus, float **jptribbleageerrorsminus, int setmodelres) {

  char line_storage[1024], buffer[1024], readstr1[100], readstr2[100], readstr3[100];
  int line_num = 1, foundbestage = 0, foundbestnorm = 0, foundchisquared = 0, foundmodelalpha = 0, foundmodelflux = 0, foundageerrorsplus = 0, foundageerrorsminus = 0, bestageline, bestnormline, chisquaredline, modelalphaline, modelfluxline, ageerrorsplusline, ageerrorsminusline, tmp_regno, tmp_modresno;

  FILE *importfile = fopen(importfilename, "r");

  printf("Reading the JPTRIB model fitting data...\n");

  if (importfile != NULL) {

    while( fgets(line_storage, sizeof(line_storage), importfile) != NULL )  {

      sscanf(line_storage,"%s", buffer);


      if (foundbestage == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - bestageline) ) ) {
	  printf("\n*** ERROR: Finished reading the JPTRIBBESTAGE before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - bestageline) ) ) {
	  printf("JPTRIBBESTAGE imported successfully\n");
	  foundbestage = 2;
	}
	else {
	  sscanf(buffer,"%f", &jptribbestage[currentset][(line_num - bestageline)]);
	  // printf("JPTRIBBESTAGE of region %d is %f (taken from the string %s)\n", (line_num - bestageline), jptribbestage[currentset][(line_num - bestageline)], buffer);
	}

      }
    

      if (foundbestnorm == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - bestnormline) ) ) {
	  printf("\n*** ERROR: Finished reading the JPTRIBBESTNORM before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - bestnormline) ) ) {
	  printf("JPTRIBBESTNORM imported successfully\n");
	  foundbestnorm = 2;
	}
	else {
	  sscanf(buffer,"%f", &jptribbestnorm[currentset][(line_num - bestnormline)]);
	  // printf("JPTRIBBESTNORM of region %d is %.6e (taken from the string %s)\n", (line_num - bestnormline), jptribbestnorm[currentset][(line_num - bestnormline)], buffer);
	}

      }
    

      if (foundchisquared == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - chisquaredline) ) ) {
	  printf("\n*** ERROR: Finished reading the JPTRIBCHISQUARED before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - chisquaredline) ) ) {
	  printf("JPTRIBCHISQUARED imported successfully\n");
	  foundchisquared = 2;
	}
	else {
	  sscanf(buffer,"%f", &jptribchisquared[currentset][(line_num - chisquaredline)]);
	  // printf("JPTRIBCHISQUARED of region %d is %.6e (taken from the string %s)\n", (line_num - chisquaredline), jptribchisquared[currentset][(line_num - chisquaredline)], buffer);
	}

      }


      if (foundmodelalpha == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - modelalphaline) ) ) {
	  printf("\n*** ERROR: Finished reading the JPTRIBMODELALPHA before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - modelalphaline) ) ) {
	  printf("JPTRIBMODELALPHA imported successfully\n");
	  foundmodelalpha = 2;
	}
	else {
	  sscanf(buffer,"%f", &jptribmodelalpha[currentset][(line_num - modelalphaline)]);
	  // printf("JPTRIBMODELALPHA of region %d is %.6e (taken from the string %s)\n", (line_num - modelalphaline), jptribmodelalpha[currentset][(line_num - modelalphaline)], buffer);
	}

      }


      if (foundageerrorsplus == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - ageerrorsplusline) ) ) {
	  printf("\n*** ERROR: Finished reading the JPTRIBBLEAGEERRORSPLUS before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - ageerrorsplusline) ) ) {
	  printf("JPTRIBBLEAGEERRORSPLUS imported successfully\n");
	  foundageerrorsplus = 2;
	}
	else {
	  sscanf(buffer,"%f", &jptribbleageerrorsplus[currentset][(line_num - ageerrorsplusline)]);
	  // printf("JPTRIBBLEAGEERRORSPLUS of region %d is %.6e (taken from the string %s)\n", (line_num - ageerrorsplusline), jptribbleageerrorsplus[currentset][(line_num - ageerrorsplusline)], buffer);
	}

      }


      if (foundageerrorsminus == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - ageerrorsminusline) ) ) {
	  printf("\n*** ERROR: Finished reading the JPTRIBBLEAGEERRORSMINUS before the expected number of regions had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - ageerrorsminusline) ) ) {
	  printf("JPTRIBBLEAGEERRORSMINUS imported successfully\n");
	  foundageerrorsminus = 2;
	}
	else {
	  sscanf(buffer,"%f", &jptribbleageerrorsminus[currentset][(line_num - ageerrorsminusline)]);
	  // printf("JPTRIBBLEAGEERRORSMINUS of region %d is %.6e (taken from the string %s)\n", (line_num - ageerrorsminusline), jptribbleageerrorsminus[currentset][(line_num - ageerrorsminusline)], buffer);
	}

      }


      if (foundmodelflux == 1) {

	if ( (strcmp(buffer,"###END###") == 0) && ( ((regnumber*setmodelres)+regnumber+1) != (line_num - modelfluxline) ) ) {
	  printf("\n*** ERROR: Finished reading the MODELFLUX before the expected number of regions had been reached (full format)! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ( ((regnumber*setmodelres)+regnumber+1) == (line_num - modelfluxline) ) ) {
	  printf("MODELFLUX imported successfully\n");
	  foundmodelflux = 2;
	}
	else {

	  sscanf(buffer,"%100[^,],%100[^,],%s", readstr1, readstr2, readstr3);
	  
	  tmp_regno = atoi(readstr1);
	  tmp_modresno = atoi(readstr2);

	  // (line_num - modelfluxline) - ((tmp_regno-1)*setmodelres)

	    
	  jptribmodelflux[currentset][tmp_regno][tmp_modresno] = atof(readstr3);
	  
	  // printf("MODELFLUX of region %d is %.6e for res point %d (taken from the string %s)\n", tmp_regno,  jptribmodelflux[currentset][tmp_regno][tmp_modresno], tmp_modresno, buffer);

	}

      }


      // If we are done, go back to the control script so we don't waste time or hit the EOF
      if ( (foundbestage == 2) && (foundbestnorm == 2) && (foundchisquared == 2) && (foundmodelalpha == 2) && (foundmodelflux == 2) && (foundageerrorsplus == 2) && (foundageerrorsminus == 2) ) {
	return 0;
      }


      if (strcmp(buffer,"###JPTRIBBESTAGE###") == 0)  {
	printf("JPTRIBBESTAGE found on line %d\n", line_num);
	foundbestage = 1;
	bestageline = line_num;
      }

      if (strcmp(buffer,"###JPTRIBBESTNORM###") == 0)  {
	printf("JPTRIBBESTNORM found on line %d\n", line_num);
	foundbestnorm = 1;
	bestnormline = line_num;
      }

      if (strcmp(buffer,"###JPTRIBCHISQUARED###") == 0)  {
	printf("JPTRIBCHISQUARED found on line %d\n", line_num);
	foundchisquared = 1;
	chisquaredline = line_num;
      }

      if (strcmp(buffer,"###JPTRIBMODELALPHA###") == 0)  {
	printf("JPTRIBMODELALPHA found on line %d\n", line_num);
	foundmodelalpha = 1;
	modelalphaline = line_num;
      }
      
      if (strcmp(buffer,"###JPTRIBMODELFLUX###") == 0)  {
	printf("JPTRIBMODELFLUX found on line %d\n", line_num);
	foundmodelflux = 1;
	modelfluxline = line_num;
      }

      if (strcmp(buffer,"###JPTRIBBLEAGEERRORSPLUS###") == 0)  {
	printf("JPTRIBBLEAGEERRORSPLUS found on line %d\n", line_num);
	foundageerrorsplus = 1;
	ageerrorsplusline = line_num;
      }

      if (strcmp(buffer,"###JPTRIBBLEAGEERRORSMINUS###") == 0)  {
	printf("JPTRIBBLEAGEERRORSMINUS found on line %d\n", line_num);
	foundageerrorsminus = 1;
	ageerrorsminusline = line_num;
      }

      
      if (strcmp(buffer,"###ENDOFIMPORTFILE###") == 0)  {
	printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EFI1014) ***\n\n");
	fclose(importfile);
	return 404;
      }

      line_num++;
    }
  
    fclose(importfile);



    if ( (foundbestage == 0) || (foundbestnorm == 0) || (foundchisquared == 0) || (foundmodelalpha == 0) || (foundmodelflux == 0) || (foundageerrorsplus == 0) || (foundageerrorsminus == 0) ) {
      printf("\n*** ERROR: Unable to find all of the required region fields! The data may be corrupt ***\n\n");
      return 404;
    }
    

  }
  else {
    fprintf(stderr,"*** Error: Unable to locate or access the .brats file requested ***\n\n");
    return 404;
  }

  return 0;

}



int importmapprops(char *importfilename, int currentset, float **eq, float **delt1, float **delt2, float **crpix1, float **crpix2, float **crota1, float **crota2, int imgnum) {

  char line_storage[1024], buffer[1024];
  int line_num = 1, foundeq = 0, founddelt1 = 0, founddelt2 = 0, foundcrpix1 = 0, foundcrpix2 = 0, foundcrota1 = 0, foundcrota2 = 0, eqline, delt1line, delt2line, crpix1line, crpix2line, crota1line, crota2line;	

  FILE *importfile = fopen(importfilename, "r");

  eq[currentset] = (float *)calloc(imgnum,sizeof(float));
  delt1[currentset] = (float *)calloc(imgnum,sizeof(float));
  delt2[currentset] = (float *)calloc(imgnum,sizeof(float));
  crpix1[currentset] = (float *)calloc(imgnum,sizeof(float));
  crpix2[currentset] = (float *)calloc(imgnum,sizeof(float));
  crota1[currentset] = (float *)calloc(imgnum,sizeof(float));
  crota2[currentset] = (float *)calloc(imgnum,sizeof(float));

  printf("Reading the coordinates...\n");

  if (importfile != NULL) {

    while( fgets(line_storage, sizeof(line_storage), importfile) != NULL )  {

      sscanf(line_storage,"%s", buffer);

      if (foundeq == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - eqline) ) ) {
	  printf("\n*** ERROR: Finished reading the EQ before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - eqline) ) ) {
	  printf("EQ imported successfully\n");
	  foundeq = 2;
	}
	else {
	  sscanf(buffer,"%f", &eq[currentset][(line_num - eqline) - 1]);
	}
      }

      if (founddelt1 == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - delt1line) ) ) {
	  printf("\n*** ERROR: Finished reading the DELT1 before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - delt1line) ) ) {
	  printf("DELT1 imported successfully\n");
	  founddelt1 = 2;
	}
	else {
	  sscanf(buffer,"%f", &delt1[currentset][(line_num - delt1line) - 1]);
	}
      }

    if (founddelt2 == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - delt2line) ) ) {
	  printf("\n*** ERROR: Finished reading the DELT2 before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - delt2line) ) ) {
	  printf("DELT2 imported successfully\n");
	  founddelt2 = 2;
	}
	else {
	  sscanf(buffer,"%f", &delt2[currentset][(line_num - delt2line) - 1]);
	}
      }


      if (foundcrpix1 == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - crpix1line) ) ) {
	  printf("\n*** ERROR: Finished reading the CRPIX1 before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - crpix1line) ) ) {
	  printf("CRPIX1 imported successfully\n");
	  foundcrpix1 = 2;
	}
	else {
	  sscanf(buffer,"%f", &crpix1[currentset][(line_num - crpix1line) - 1]);
	}
      }

    if (foundcrpix2 == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - crpix2line) ) ) {
	  printf("\n*** ERROR: Finished reading the CRPIX2 before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - crpix2line) ) ) {
	  printf("CRPIX2 imported successfully\n");
	  foundcrpix2 = 2;
	}
	else {
	  sscanf(buffer,"%f", &crpix2[currentset][(line_num - crpix2line) - 1]);
	}
      }


      if (foundcrota1 == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - crota1line) ) ) {
	  printf("\n*** ERROR: Finished reading the CROTA1 before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - crota1line) ) ) {
	  printf("CROTA1 imported successfully\n");
	  foundcrota1 = 2;
	}
	else {
	  sscanf(buffer,"%f", &crota1[currentset][(line_num - crota1line) - 1]);
	}
      }

    if (foundcrota2 == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) != (line_num - crota2line) ) ) {
	  printf("\n*** ERROR: Finished reading the CROTA2 before the expected number of maps had been reached! The data may be corrupt ***\n\n");
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((imgnum+1) == (line_num - crota2line) ) ) {
	  printf("CROTA2 imported successfully\n");
	  foundcrota2 = 2;
	}
	else {
	  sscanf(buffer,"%f", &crota2[currentset][(line_num - crota2line) - 1]);
	}
      }


      // If we are done, go back to the control script so we don't waste time or hit the EOF
      if ( (foundeq == 2) && (founddelt1 == 2) && (founddelt2 == 2) && (foundcrpix1 == 2) && (foundcrpix2 == 2) && (foundcrota1 == 2) && (foundcrota2 == 2) ) {
	return 0;
      }

      if (strcmp(buffer,"###EQ###") == 0)  {
	printf("EQ found on line %d\n", line_num);
	foundeq = 1;
	eqline = line_num;
      }

      if (strcmp(buffer,"###DELT1###") == 0)  {
	printf("DELT1 found on line %d\n", line_num);
	founddelt1 = 1;
	delt1line = line_num;
      }

      if (strcmp(buffer,"###DELT2###") == 0)  {
	printf("DELT2 found on line %d\n", line_num);
	founddelt2 = 1;
	delt2line = line_num;
      }

      if (strcmp(buffer,"###CRPIX1###") == 0)  {
	printf("CRPIX1 found on line %d\n", line_num);
	foundcrpix1 = 1;
	crpix1line = line_num;
      }

      if (strcmp(buffer,"###CRPIX2###") == 0)  {
	printf("CRPIX2 found on line %d\n", line_num);
	foundcrpix2 = 1;
	crpix2line = line_num;
      }

      if (strcmp(buffer,"###CROTA1###") == 0)  {
	printf("CROTA1 found on line %d\n", line_num);
	foundcrota1 = 1;
	crota1line = line_num;
      }

      if (strcmp(buffer,"###CROTA2###") == 0)  {
	printf("CROTA2 found on line %d\n", line_num);
	foundcrota2 = 1;
	crota2line = line_num;
      }

      if (strcmp(buffer,"###ENDOFIMPORTFILE###") == 0)  {
	printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EFI1006) ***\n\n");
	fclose(importfile);
	return 404;
      }

      line_num++;
    }

    fclose(importfile);

    if ( (foundeq == 0) || (founddelt1 == 0) || (founddelt2 == 0) || (foundcrpix1 == 0) || (foundcrpix2 == 0) || (foundcrota1 == 0) || (foundcrota2 == 0) ) {
      printf("\n*** ERROR: Unable to find all of the map property fields! The data may be corrupt ***\n\n");
      return 404;
    }

  }
  else {
	fprintf(stderr,"*** Error: Unable to locate or access the .brats file requested ***\n\n");
	return 404;
      }

  return 0;
}


int convertmapprops(char *importfilename, int currentset, int exportversion, float **eq, float **delt1, float **delt2, float **crpix1, float **crpix2, float **crota1, float **crota2, int imgnum, float cellsize, int xdim, int ydim) {

  char  *cmdbuffer, *endptr;
  int i, validentry, eqselect = -1;	

  eq[currentset] = (float *)calloc(imgnum,sizeof(float));
  delt1[currentset] = (float *)calloc(imgnum,sizeof(float));
  delt2[currentset] = (float *)calloc(imgnum,sizeof(float));
  crpix1[currentset] = (float *)calloc(imgnum,sizeof(float));
  crpix2[currentset] = (float *)calloc(imgnum,sizeof(float));
  crota1[currentset] = (float *)calloc(imgnum,sizeof(float));
  crota2[currentset] = (float *)calloc(imgnum,sizeof(float));

  validentry = 0;

  while (validentry == 0) {

    cmdbuffer = readline("Unable to automatically determine the epoch of the coordinate system. Is the data in B1950 (0) or J2000 (1)?: ");
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    chomp(cmdbuffer);

    eqselect = strtol(cmdbuffer, &endptr, 10);

    if ( (eqselect == 0) && (cmdbuffer == endptr) ) {
      printf("\nInvalid input, please try again...\n\n");
      continue;
    }
    else if (eqselect < 0) {
      printf("Escaping command...\n");
      return 404;
    }
    else if ( (eqselect != 0) && (eqselect != 1) ) {
      printf("\nInvalid selection, please try again...\n\n");
      continue;
    }
    else if (eqselect == 0) {
      for(i=0; i<imgnum; i++) {
	eq[currentset][i] = 1950;
      }
      validentry = 1;
    }
    else if (eqselect == 1) {
      for(i=0; i<imgnum; i++) {
	eq[currentset][i] = 2000;
      }
      validentry = 1;
    }
  }

  printf("Automatically determining DELT1 and DELT2 from the reference map cell size...\n");

  for(i=0; i<imgnum; i++) {
    delt1[currentset][i] = cellsize / 3600;
    delt2[currentset][i] = cellsize / 3600;
  }

 printf("No information available for CRPIX1 and CRPIX2. Assuming they are the standard midpoints...\n");

  for(i=0; i<imgnum; i++) {
    crpix1[currentset][i] = (ydim / 2);
    crpix2[currentset][i] = (xdim / 2) + 1;
  }

  printf("No rotational information available. Assuming CROTA1 and CROTA2 are 0...\n");

  for(i=0; i<imgnum; i++) {
    crota1[currentset][i] = 0.0;
    crota2[currentset][i] = 0.0;
  }

  return 0;
}


int importfindinject(char *importfilename, float *inj_sumchisquared, float *inj_regbestinject, float *inj_regchisquared, int mininjectset, float *mininjectstore, float *maxinjectstore, int model, int regnumber, int exportversion) {

  char line_storage[1024], buffer[1024], str_model[32], comp_store[64], comp_injsum[64], comp_regchiinject[64], comp_regbestinject[64];
  int line_num = 1, foundstore = 0, foundinjsum = 0, foundregbestinject = 0, foundregchiinject = 0, storeline, injsumline, regbestinjectline, regchiinjectline, expectedheaders = 3, tmp_injectint; 	

  FILE *importfile;

  // Set up the comparison strings for later
  sprintf(str_model, "%d", model);

  strcpy(comp_store, "###MININJECTHEADMOD");
  strcat(comp_store, str_model);
  strcat(comp_store, "###");

  strcpy(comp_injsum, "###INJ_SUMCHISQUAREDMOD");
  strcat(comp_injsum, str_model);
  strcat(comp_injsum, "###");

  strcpy(comp_regbestinject, "###INJ_REGBESTINJECTMOD");
  strcat(comp_regbestinject, str_model);
  strcat(comp_regbestinject, "###");

  strcpy(comp_regchiinject, "###INJ_REGCHISQUAREDMOD");
  strcat(comp_regchiinject, str_model);
  strcat(comp_regchiinject, "###");


  printf("Reading the findinject data...\n");

  if ( (importfile = fopen(importfilename,"r") ) != NULL ) {

    while( fgets(line_storage, sizeof(line_storage), importfile) != NULL )  {

      sscanf(line_storage,"%s", buffer);

      if (foundstore == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((2) != (line_num - storeline) ) ) {
	  printf("\n*** ERROR: Finished reading MININJECTHEADMOD%d before the expected number of maps had been reached! The data may be corrupt ***\n\n", model);
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((2) == (line_num - storeline) ) ) {
	  printf("MININJECTHEADMOD%d imported successfully\n", model);
	  foundstore = 2;
	}
	else {
	  if (sscanf(line_storage,"%d,%f,%f", &tmp_injectint, mininjectstore, maxinjectstore) != expectedheaders) {
	    printf("*** Error: Found an unexpected number of header entries for MININJECTHEADMOD%d. The data may be corrupt ***\n", model);
	    fclose(importfile);
	    return 404;
	  }
	  else if (tmp_injectint != mininjectset) {
	    printf("*** Error: The number of injection intervals in the main header does not match that of MININJECTHEADMOD%d. The data may be corrupt ***\n", model);
	    fclose(importfile);
	    return 404;
	  }

	}
      }

      if (foundinjsum == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((mininjectset+2) != (line_num - injsumline) ) ) { // +2 as the number of datapoints is mininjectset +1
	  printf("\n*** ERROR: Finished reading the INJ_SUMCHISQUAREDMOD%d before the expected number of intervals had been reached! The data may be corrupt ***\n\n", model);
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((mininjectset+2) == (line_num - injsumline) ) ) {
	  printf("INJ_SUMCHISQUAREDMOD%d imported successfully\n", model);
	  foundinjsum = 2;
	}
	else {
	  sscanf(buffer,"%f", &inj_sumchisquared[(line_num - injsumline) - 1]);
	  printf("%f\n", inj_sumchisquared[(line_num - injsumline) - 1]);
	}
      }

      if (foundregbestinject == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - regbestinjectline) ) ) {
	  printf("\n*** ERROR: Finished reading the INJ_REGBESTINJECTMOD%d before the expected number of regions had been reached! The data may be corrupt ***\n\n", model);
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - regbestinjectline) ) ) {
	  printf("INJ_REGBESTINJECTMOD%d imported successfully\n", model);
	  foundregbestinject = 2;
	}
	else {
	  sscanf(buffer,"%f", &inj_regbestinject[(line_num - regbestinjectline)]); // No -1 as we are not indexed from 0
	}
      }

      if (foundregchiinject == 1) {
	if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) != (line_num - regchiinjectline) ) ) {
	  printf("\n*** ERROR: Finished reading the INJ_REGCHISQUAREDMOD%d before the expected number of regions had been reached! The data may be corrupt ***\n\n", model);
	  fclose(importfile);
	  return 404;
	}
	else if ( (strcmp(buffer,"###END###") == 0) && ((regnumber+1) == (line_num - regchiinjectline) ) ) {
	  printf("INJ_REGCHISQUAREDMOD%d imported successfully\n", model);
	  foundregchiinject = 2;
	}
	else {
	  sscanf(buffer,"%f", &inj_regchisquared[(line_num - regchiinjectline)]); // No -1 as we are not indexed from 0
	}
      }

      // If we are done, go back to the control script so we don't waste time or hit the EOF
      if (exportversion == 1) {
	if ( (foundstore == 2) && (foundinjsum == 2) ) {
	  return 0;
	}
      }
      else {
	if ( (foundstore == 2) && (foundinjsum == 2) && (foundregbestinject == 2) && (foundregchiinject == 2) ) {	  
	  return 0;
	}
      }

      if (strcmp(buffer, comp_store) == 0)  {
	printf("MININJECTHEADMOD%d found on line %d\n", model, line_num);
	foundstore = 1;
	storeline = line_num;
      }

      if (strcmp(buffer, comp_injsum) == 0)  {
	printf("SUMCHISQUAREDMOD%d found on line %d\n", model, line_num);
	foundinjsum = 1;
	injsumline = line_num;
      }

      if (strcmp(buffer, comp_regbestinject) == 0)  {
	printf("REGBESTINJECTMOD%d found on line %d\n", model, line_num);
	foundregbestinject = 1;
	regbestinjectline = line_num;
      }

      if (strcmp(buffer, comp_regchiinject) == 0)  {
	printf("REGCHISQUAREDMOD%d found on line %d\n", model, line_num);
	foundregchiinject = 1;
	regchiinjectline = line_num;
      }

      //printf("%d %d %d %d\n", foundstore, foundinjsum, foundregbestinject, foundregchiinject);
      
      if (strcmp(buffer,"###ENDOFIMPORTFILE###") == 0)  {
	printf("\n*** ERROR: Reached the end of the file without finding all of the required fields! The data may be corrupt (EFI1007) ***\n\n");
	fclose(importfile);
	return 404;
      }

      line_num++;
    }

    fclose(importfile);


    if (exportversion == 1) {
      if ( (foundstore == 0) || (foundinjsum == 0) ) {
	printf("\n*** ERROR: Unable to find all of the required the findinject fields! The data may be corrupt ***\n\n");
	return 404;
      }
    }
    else {
      if ( (foundstore == 0) || (foundinjsum == 0) || (foundregbestinject == 0) || (foundregchiinject == 0) ) {
	printf("\n*** ERROR: Unable to find all of the required the findinject fields! The data may be corrupt ***\n\n");
	return 404;
      }
    }

  }
  else {
    fprintf(stderr,"*** Error: Unable to locate or access the .brats file requested ***\n\n");
    return 404;
  }
  
  return 0;
}
