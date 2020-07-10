/*

  Scales the flux of a dataset or subset maps within a dataset

  Requires the file "chomp.h" and <funtools.h> plus the standard headers to have been previously declared

*/

// Declare the function prototypes

int syscom(const char[]);

int scalefluxes(int currentset, int mapopt, float minfreq, float maxfreq, int singlemap, int scaleopt, float factorfreq1, float factorfreq2, float factor1, float factor2, float singlefactor, float ****flux, int xdim, int ydim, int *fluxorder, int imgnum, float *frequency, float ***bgbuff, float **rms, float beamarea) { 


  // Declare the variables for this task

  int i, j, a, startmap, finishmap, count;
  float tmp_minfreq, tmp_maxfreq, tmp_scalefactor, averagebackground, rmsnumerator, bgflux;

  int escape = 0;

  
  //Get the directory name and format it

  printf("\nScaling data set %d\n", currentset);

  if (mapopt == 1) {

    printf("Scaling map %d by a factor of %.4f\n", singlemap, singlefactor);

    tmp_scalefactor = singlefactor;

    for (i=0; i<xdim; i++) {
      for (j=0; j<ydim; j++){
	flux[currentset][singlemap][i][j] *= tmp_scalefactor;
      }
    }
      
    // Reset the temp RMS values for each map level pass
    bgflux = 0.0;
    count = 0;
    averagebackground = 0.0;
    rmsnumerator = 0.0;

    // Calculate the new background flux
    for (i=0; i<xdim; i++){
      for (j=0; j<ydim; j++){
	bgbuff[currentset][singlemap][i*xdim+j] *= tmp_scalefactor; // Rescale the buffer
	bgflux += (bgbuff[currentset][singlemap][i*xdim+j] / beamarea);
	if ( ( bgbuff[currentset][singlemap][i*xdim+j] < -1e-20) || ( bgbuff[currentset][singlemap][i*xdim+j] > 1e-20) ) {
	  count++;
	}
      }
    }

    //printf("Count: %d \n", count);

    averagebackground = bgflux / count;

    // Find the new RMS
    for (i=0; i<xdim; i++){
      for (j=0; j<ydim; j++){
	// Ensure we are in the region
	if ( (bgbuff[currentset][singlemap][i*xdim+j] < -1e-20) || (bgbuff[currentset][singlemap][i*xdim+j] > 1e-20) ) {
	  rmsnumerator += pow((bgbuff[currentset][singlemap][i*xdim+j] / beamarea) - averagebackground, 2);
	}
      }
    }

    rms[currentset][singlemap] =  sqrt(rmsnumerator / count); // This gives RMS in a per pixel format
    //printf("New RMS: %.2e Hz = %.4e Jy/Pixel  (%.4e Jy/Beam)\n", frequency[a], rms[currentset][a], rms[currentset][a]*beamarea);

    printf("\nRMS is now:\n");
    printf("%.2e Hz = %.4e Jy/Pixel  (%.4e Jy/Beam)\n", frequency[singlemap], rms[currentset][singlemap], rms[currentset][singlemap]*beamarea);
 
  }



  else if ( (mapopt == 2) || (mapopt == 3) ) {

    // Get the start, finish and number of maps to be processed if we are not doing them all
    if (mapopt == 3) {
      startmap = 0;
      finishmap = imgnum;
    }
    else if (mapopt == 2) {

      startmap = -1;

      tmp_minfreq = 1e27;
      tmp_maxfreq = -1e27;
    
      for (a=0; a<imgnum; a++){

	if ( (frequency[fluxorder[a]] >= minfreq) && (frequency[fluxorder[a]] <= maxfreq) ) {
	  if (frequency[fluxorder[a]] < tmp_minfreq) {
	    startmap = a; // How far have we moved along the flux order rows?
	    tmp_minfreq = frequency[fluxorder[a]];
	    //printf("startmap: %d\n", startmap);
	  }
	  if (frequency[fluxorder[a]] > tmp_maxfreq) {
	    finishmap = a+1; // Add one to correct for the imgnum indexing convention
	    tmp_maxfreq = frequency[fluxorder[a]];
	    //printf("finishmap: %d\n", finishmap);
	  }
	}
      }

      //Escape here if no maps in freq range
      if (startmap == -1) {
	printf("\n*** ERROR: No maps in frequency range. Exiting...\n\n");
	escape = 1;
	return 404;
      }
    }

    // Break the second loop if there is an error
    if (escape == 1) {
      return 404;
    }


      for (a=startmap; a<finishmap; a++){

	if (scaleopt == 1) {

	  printf("Scaling map %d by a factor of %.4f\n", fluxorder[a], singlefactor);

	  tmp_scalefactor = singlefactor;

	  for (i=0; i<xdim; i++){
	    for (j=0; j<ydim; j++){
	      flux[currentset][fluxorder[a]][i][j] *= singlefactor; // Get the scaled flux
	    }
	  }
	}
      
	else if (scaleopt == 2) {

	  tmp_scalefactor = ( factor1 + (factor2 - factor1) * ((frequency[fluxorder[a]] - factorfreq1) / (factorfreq2-factorfreq1)) );
	
	  printf("Scaling map %d by a factor of %.4f (linear interp.)\n", fluxorder[a], tmp_scalefactor);

	  //printf("factor1:%.4f factor2:%.4f frequency[fluxorder[a]]:%.3e minfreq:%.4e maxfreq:%.4e\n", factor1, factor2, frequency[fluxorder[a]], factorfreq1, factorfreq2);

	  for (i=0; i<xdim; i++){
	    for (j=0; j<ydim; j++){
	      flux[currentset][fluxorder[a]][i][j] *= tmp_scalefactor; // Get the scaled flux
	    }
	  }
	}

      
	// Reset the temp RMS values for each map level pass
	bgflux = 0.0;
	count = 0;
	averagebackground = 0.0;
	rmsnumerator = 0.0;

	// Calculate the new background flux
	for (i=0; i<xdim; i++){
	  for (j=0; j<ydim; j++){
	    bgbuff[currentset][fluxorder[a]][i*xdim+j] *= tmp_scalefactor; // Rescale the buffer
	    bgflux += (bgbuff[currentset][fluxorder[a]][i*xdim+j] / beamarea);
	    if ( (bgbuff[currentset][fluxorder[a]][i*xdim+j] < -1e-20) || (bgbuff[currentset][fluxorder[a]][i*xdim+j] > 1e-20) ) {
	      count++;
	    }
	  }
	}

	averagebackground = bgflux / count;

	//printf("Count: %d \n", count);
    

	// Find the new RMS
	for (i=0; i<xdim; i++){
	  for (j=0; j<ydim; j++){
	    // Ensure we are in the region
	    if ( (bgbuff[currentset][fluxorder[a]][i*xdim+j] < -1e-20) || (bgbuff[currentset][fluxorder[a]][i*xdim+j] > 1e-20) ) {
	      rmsnumerator += pow( (bgbuff[currentset][fluxorder[a]][i*xdim+j] / beamarea) - averagebackground, 2);
	    }
	  }
	}

	rms[currentset][fluxorder[a]] =  sqrt(rmsnumerator / count); // This gives RMS in a per pixel format
	//printf("New RMS: %.2e Hz = %.4e Jy/Pixel  (%.4e Jy/Beam)\n", frequency[a], rms[currentset][a], rms[currentset][a]*beamarea);
	
      }

      printf("\nRMS values are now:\n");

      for (a=0; a < imgnum; a++) {
	printf("%.2e Hz = %.4e Jy/Pixel  (%.4e Jy/Beam)\n", frequency[fluxorder[a]], rms[currentset][fluxorder[a]], rms[currentset][fluxorder[a]]*beamarea);
      }

  }


    printf("\n *** WARNING: Scaling is applied to the raw flux tables. Region selection should be re-run before any other commands if this is to be applied in a meaningful way e.g. setregions ***\n\n");


  return 0;

}
