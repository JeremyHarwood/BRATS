/* Set the percentage flux error based on instrument */


int fluxerrors(float *fluxcalerror, float *frequency, int imgnum, int instrument, float flaterror) {

  int i;

  if (instrument == 0) { // Flat percentage
    printf("Flux calibration for all maps set to a flat %.3f%% error\n\n", flaterror *100);

    for(i=0; i<imgnum; i++) {
      fluxcalerror[i] = flaterror;
    }
  }

  else if (instrument == 1) { // EVLA

    printf("Flux calibration errors set to standard EVLA values (2%% at <12 GHz, 5%% at > 12 GHz)\n\n");

    for(i=0; i<imgnum; i++) {
      if (frequency[i] > 1.2e+10) {
	fluxcalerror[i] = 0.05;
      }
      else {
	fluxcalerror[i] = 0.02;
      }
    }
  }

  else {

    printf("Unknown telescope. Defaulting to flat calibration errors...\n\n");

    for(i=0; i<imgnum; i++) {
      fluxcalerror[i] = flaterror;
    }
  }


  return 0;

}
