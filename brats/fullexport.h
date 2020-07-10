/*

  Exports the full dataset in .brats format 

*/

// Declare the function prototypes

#include <time.h>

int fullexport(char *filename, int exportcompression, char *settarget, int imgnum, float redshift, int xdim, int ydim, float bmaj, float bmin, float bpa, float beamarea, float *ra, float *dec, float *eq, float *delt1, float *delt2, float *crpix1, float *crpix2, float *crota1, float *crota2, float cellsize, char *setname, char *setbg, char *setreg, int posmap, char *dataloc, int regionsset, int jpmodelset, int kpmodelset, int jptribmodelset, int *mininjectset, float *frequency, float *frequencyobs, int *fluxorder, float *rms, float ***flux, float *minflux, float *maxflux, float **bgbuff, float *fluxcalerror, int setaveraged, int regnumber, float *regionlocx, float *regionlocy, float *regminflux, float *regmaxflux, float **fluxerror, int **regionarray, float **regflux, float *mininjectstore, float *maxinjectstore, float **inj_sumchisquared, int *setmodelres, float jpfield, float jpinject, float *jpbestage, float *jpbestnorm, float *jpchisquared, float *jpmodelalpha, float **jpmodelflux, float *jpageerrorsplus, float *jpageerrorsminus, float kpfield, float kpinject, float *kpbestage, float *kpbestnorm, float *kpchisquared, float *kpmodelalpha, float **kpmodelflux, float *kpageerrorsplus, float *kpageerrorsminus, float jptribfield, float jptribinject, float *jptribbestage, float *jptribbestnorm, float *jptribchisquared, float *jptribmodelalpha, float **jptribmodelflux, float *jptribbleageerrorsplus, float *jptribbleageerrorsminus, int *regionsize, float **inj_regbestinject, float **inj_regchisquared) { 


  // Declare the variables for this task

  int i, j, a;

  FILE *filestream;

  // Open up a new file and start exporting

  printf("Exporting data, please wait... ");

  if ( (filestream = fopen(filename,"w") ) == NULL ) {

    printf("\n*** Error: Cannot open specified directory to output data (%s)! Please check the directory exists, the correct permission are set and the location is correctly set in dataloc ***\n\n", dataloc);
    return 404;
  }

  // If there is no target in the fits header, set a dummy one up (a blank target sometimes causes a seg fault on import).
  if(settarget[0] == '\0') {
    strcpy(settarget, "NoNameFoundInFITSHeader");
  }

  // Write the prelims
  fprintf(filestream, "######## EXPORTED BRATS DATA FOR THE SOURCE %s ########", settarget);
  fprintf(filestream, "\n\n");

  //Write thea user readable header
  fprintf(filestream, "######## USER READABLE HEADER ########\n\n");
  fprintf(filestream, "FORMAT: BRATSV2\n");
  if (exportcompression == 1) {
    fprintf(filestream, "EXPORTTYPE: Compact\n");
  }
  else {
    fprintf(filestream, "EXPORTTYPE: Full\n");
  }
  fprintf(filestream, "TARGET: %s\n", settarget);
  fprintf(filestream, "REDSHIFT: %.6f\n", redshift);
  fprintf(filestream, "RA: %.12e\n", ra[posmap]);
  fprintf(filestream, "DEC: %.12e\n", dec[posmap]);
  fprintf(filestream, "EQ: %.12e\n", eq[posmap]);
  fprintf(filestream, "IMAGES: %d\n", imgnum);
  fprintf(filestream, "XDIM: %d\n", xdim);
  fprintf(filestream, "YDIM: %d\n", ydim);
  fprintf(filestream, "BMAJ: %.6f\n", bmaj);
  fprintf(filestream, "BMIN: %.6f\n", bmin);
  fprintf(filestream, "BPA: %.4f\n", bpa);
  fprintf(filestream, "CRPIX1: %.4f\n", crpix1[posmap]);
  fprintf(filestream, "CRPIX2: %.4f\n", crpix2[posmap]);
  fprintf(filestream, "CROTA1: %.4f\n", crota1[posmap]);
  fprintf(filestream, "CROTA2: %.4f\n", crota2[posmap]);
  fprintf(filestream, "BEAMAREA: %.6f\n", beamarea);
  fprintf(filestream, "CELLSIZE: %.6f\n", cellsize);
  fprintf(filestream, "SETLOC: %s\n", setname);
  fprintf(filestream, "REG: %s\n", setreg);
  fprintf(filestream, "BGREG: %s\n", setbg);
  fprintf(filestream, "\n");
  fprintf(filestream, "#################################\n\n");


  // Start writing what the software will read
  fprintf(filestream, "######## DATA ########\n\n");
  fprintf(filestream, "!!!!!!!! Warning: Any adjustment to the format or content of this file beyond this point may cause importing to fail !!!!!!!!");
  fprintf(filestream, "\n\n");
  
  // Write the general header info
  fprintf(filestream, "###GENHEAD###\n");
  fprintf(filestream, "BRATSV2\n");
  if (exportcompression == 1) {
    fprintf(filestream, "Compact\n");
  }
  else {
    fprintf(filestream, "Full\n");
  }
  fprintf(filestream, "1HEAD1,%s,%.6f,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%s,%s\n", settarget, redshift, imgnum, xdim, ydim, bmaj, bmin, bpa, beamarea, cellsize, setname, setreg, setbg);
  fprintf(filestream, "###END###\n");

  // Write the RA and DEC of each map
  fprintf(filestream, "###RA###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%.12e\n", ra[i]);
  }
  fprintf(filestream, "###END###\n");

  fprintf(filestream, "###DEC###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%.12e\n", dec[i]);
  }
  fprintf(filestream, "###END###\n");

  // Write the EQ of each map
  fprintf(filestream, "###EQ###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%.12e\n", eq[i]);
  }
  fprintf(filestream, "###END###\n");

 // Write the rotation of each map
  fprintf(filestream, "###CROTA1###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%.12e\n", crota1[i]);
  }
  fprintf(filestream, "###END###\n");

  fprintf(filestream, "###CROTA2###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%.12e\n", crota2[i]);
  }
  fprintf(filestream, "###END###\n");

 // Write the pixel size of each map
  fprintf(filestream, "###DELT1###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%.12e\n", delt1[i]);
  }
  fprintf(filestream, "###END###\n");

  fprintf(filestream, "###DELT2###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%.12e\n", delt2[i]);
  }
  fprintf(filestream, "###END###\n");

 // Write the coordinate reference pixel of each map
  fprintf(filestream, "###CRPIX1###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%.12e\n", crpix1[i]);
  }
  fprintf(filestream, "###END###\n");

  fprintf(filestream, "###CRPIX2###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%.12e\n", crpix2[i]);
  }
  fprintf(filestream, "###END###\n");

  // Write which functions have actually been exported and what memory was allocated at time of export
  fprintf(filestream, "###SETHEAD###\n");
  fprintf(filestream, "%d\n%d\n%d\n%d\n%d\n%d\n%d\n", regionsset, jpmodelset, kpmodelset, jptribmodelset, mininjectset[1], mininjectset[2], mininjectset[3]);
  fprintf(filestream, "###END###\n");
  
  // Write the rest and observed frequencies
  fprintf(filestream, "###FREQUENCY###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%.12e\n", frequency[i]);
  }
  fprintf(filestream, "###END###\n");
  
  fprintf(filestream, "###FREQUENCYOBS###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%.12e\n", frequencyobs[i]);
  }
  fprintf(filestream, "###END###\n");

  // Write indexing of the maps in flux order
  fprintf(filestream, "###FLUXORDER###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%d\n", fluxorder[i]);
  }
  fprintf(filestream, "###END###\n");

  // Write the RMS of each map
  fprintf(filestream, "###RMS###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%.12e\n", rms[i]);
  }
  fprintf(filestream, "###END###\n");

  // Write the flux cal errror of each map
  fprintf(filestream, "###FLUXCALERROR###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%.6f\n", fluxcalerror[i]);
  }
  fprintf(filestream, "###END###\n");

  // Write the min and max flux of each map
  fprintf(filestream, "###MINFLUX###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%.12e\n", minflux[i]);
  }
  fprintf(filestream, "###END###\n");

    fprintf(filestream, "###MAXFLUX###\n");
  for(i=0; i<imgnum; i++) {
    fprintf(filestream, "%.12e\n", maxflux[i]);
  }
  fprintf(filestream, "###END###\n");

  // Write the rawflux of each map
  fprintf(filestream, "###FLUX###\n");
  for(a=0; a<imgnum; a++) {
    for(i=0; i<xdim; i++) {
      for(j=0; j<ydim; j++) {
	if (flux[a][i][j] > 0.0000000) {
	  fprintf(filestream, "%d,%d,%d,%.12e\n", a, i, j, flux[a][i][j]);
	}
      }
    }
  }
  fprintf(filestream, "###END###\n");

  // Write the background buffer array for each map
  fprintf(filestream, "###BGBUFF###\n");
  for(a=0; a<imgnum; a++) {
    for(i=0; i<xdim; i++) {
      for(j=0; j<ydim; j++) {
	if (bgbuff[a][i*xdim+j] > 0.0000000) {
	  fprintf(filestream, "%d,%d,%d,%.12e\n", a, i, j, bgbuff[a][i*xdim+j]);
	}
      }
    }
  }
  fprintf(filestream, "###END###\n");

  // Write region related data (any check for anyting that requires regions) if they have been set
  if (regionsset == 1) {

    // Write if the regions are averaged or not
    fprintf(filestream, "###SETAVERAGED###\n");
    fprintf(filestream, "%d\n", setaveraged);
    fprintf(filestream, "###END###\n");

    // Write the number of regions
    fprintf(filestream, "###REGNUMBER###\n");
    fprintf(filestream, "%d\n", regnumber);
    fprintf(filestream, "###END###\n");

    // Write the region sizes
    fprintf(filestream, "###REGSIZE###\n");
    for(i=1; i<=regnumber; i++) {
	  fprintf(filestream, "%d\n", regionsize[i]);
    }
    fprintf(filestream, "###END###\n");

    // Write the regionloc array
    fprintf(filestream, "###REGLOC###\n");
    for(i=1; i<=regnumber; i++) {
      fprintf(filestream, "%.6f,%.6f\n", regionlocx[i], regionlocy[i]);
    }
    fprintf(filestream, "###END###\n");

    // Write the regionlocmin and max flux per region per map
    fprintf(filestream, "###REGMINMAXFLUX###\n");
    for(i=0; i<imgnum; i++) {
      fprintf(filestream, "%.12e,%.12e\n", regminflux[i], regmaxflux[i]);
    }
    fprintf(filestream, "###END###\n");

    // Write the fluxerror array
    fprintf(filestream, "###FLUXERROR###\n");
    for(a=0; a<imgnum; a++) {
      for(j=1; j<=regnumber; j++) {
	fprintf(filestream, "%d,%d,%.12e\n", a, j, fluxerror[a][j]);
      }
    }
    fprintf(filestream, "###END###\n");

    // Write the region array
    fprintf(filestream, "###REGIONARRAY###\n");
    if (exportcompression == 1) {
      for(i=0; i<xdim; i++) {
	for(j=0; j<ydim; j++) {
	  if (regionarray[i][j] > 0) {
	    fprintf(filestream, "%d,%d,%d\n", i, j, regionarray[i][j]);
	  }
	}
      }
    }
    else {
      for(i=0; i<xdim; i++) {
	for(j=0; j<ydim; j++) {
	    fprintf(filestream, "%d,%d,%d\n", i, j, regionarray[i][j]);
	}
      }
    }
    fprintf(filestream, "###END###\n");

    // Write the region fluxes
    fprintf(filestream, "###REGFLUX###\n");
    for(i=0; i<imgnum; i++) {
      for(j=1; j<=regnumber; j++) {
	fprintf(filestream, "%d,%.12e\n", i, regflux[i][j]);
      }
    }
    fprintf(filestream, "###END###\n");

    // Output the injection index minimization arrays if they have been run (note that mininjectset also acts as the injectinterval storage parameters)
    // yyyyyyyy
    if (mininjectset[1] > 0) {
      fprintf(filestream, "###MININJECTHEADMOD1###\n");
      fprintf(filestream, "%d,%.6f,%.6f\n", mininjectset[1], mininjectstore[1], maxinjectstore[1]);
      fprintf(filestream, "###END###\n");

      fprintf(filestream, "###INJ_SUMCHISQUAREDMOD1###\n");
      for(i=0; i<=mininjectset[1]; i++) {
	fprintf(filestream, "%.12e\n", inj_sumchisquared[1][i]);
      }
      fprintf(filestream, "###END###\n");

      fprintf(filestream, "###INJ_REGBESTINJECTMOD1###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", inj_regbestinject[1][i]);
      }
      fprintf(filestream, "###END###\n");
      
      fprintf(filestream, "###INJ_REGCHISQUAREDMOD1###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", inj_regchisquared[1][i]);
      }
      fprintf(filestream, "###END###\n");
    }



    if (mininjectset[2] > 0) {
      fprintf(filestream, "###MININJECTHEADMOD2###\n");
      fprintf(filestream, "%d,%.6f,%.6f\n", mininjectset[2], mininjectstore[2], maxinjectstore[2]);
      fprintf(filestream, "###END###\n");

      fprintf(filestream, "###INJ_SUMCHISQUAREDMOD2###\n");
      for(i=0; i<=mininjectset[2]; i++) {
	fprintf(filestream, "%.12e\n", inj_sumchisquared[2][i]);
      }
      fprintf(filestream, "###END###\n");

      fprintf(filestream, "###INJ_REGBESTINJECTMOD2###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", inj_regbestinject[2][i]);
      }
      fprintf(filestream, "###END###\n");
      
      fprintf(filestream, "###INJ_REGCHISQUAREDMOD2###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", inj_regchisquared[2][i]);
      }
      fprintf(filestream, "###END###\n");
    }

    if (mininjectset[3] > 0) {
      fprintf(filestream, "###MININJECTHEADMOD3###\n");
      fprintf(filestream, "%d,%.6f,%.6f\n", mininjectset[3], mininjectstore[3], maxinjectstore[3]);
      fprintf(filestream, "###END###\n");

      fprintf(filestream, "###INJ_SUMCHISQUAREDMOD3###\n");
      for(i=0; i<=mininjectset[3]; i++) {
	fprintf(filestream, "%.12e\n", inj_sumchisquared[3][i]);
      }
      fprintf(filestream, "###END###\n");

      fprintf(filestream, "###INJ_REGBESTINJECTMOD3###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", inj_regbestinject[3][i]);
      }
      fprintf(filestream, "###END###\n");
      
      fprintf(filestream, "###INJ_REGCHISQUAREDMOD3###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", inj_regchisquared[3][i]);
      }
      fprintf(filestream, "###END###\n");
    }
  
    // Output the JP model arrays if it has been fitted
    if (jpmodelset == 1) {

      // Write the header
      fprintf(filestream, "###JPFITHEADER###\n");
      fprintf(filestream, "%.6e\n", jpfield);
      fprintf(filestream, "%.6e\n", jpinject);
      fprintf(filestream, "%d\n", setmodelres[1]);
      fprintf(filestream, "###END###\n");

      // Write the best ages
      fprintf(filestream, "###JPBESTAGE###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", jpbestage[i]);
      }
      fprintf(filestream, "###END###\n");

      // Write the best normalizations
      fprintf(filestream, "###JPBESTNORM###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", jpbestnorm[i]);
      }
      fprintf(filestream, "###END###\n");

      // Write the chi squared values
      fprintf(filestream, "###JPCHISQUARED###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", jpchisquared[i]);
      }
      fprintf(filestream, "###END###\n");

      // Write the model alpha
      fprintf(filestream, "###JPMODELALPHA###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", jpmodelalpha[i]);
      }
      fprintf(filestream, "###END###\n");

      // Write the model flux
      fprintf(filestream, "###JPMODELFLUX###\n");
      for(i=1; i<=regnumber; i++) {
	for(j=0; j<=setmodelres[1]; j++) {
	  fprintf(filestream, "%d,%d,%.12e\n", i, j, jpmodelflux[i][j]);
	}
      }
      fprintf(filestream, "###END###\n");

      // Write the model age errors
      fprintf(filestream, "###JPAGEERRORSPLUS###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", jpageerrorsplus[i]);
      }
      fprintf(filestream, "###END###\n");

      fprintf(filestream, "###JPAGEERRORSMINUS###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", jpageerrorsminus[i]);
      }
      fprintf(filestream, "###END###\n");

    }

   
    
    // Output the KP model arrays if it has been fitted
    if (kpmodelset == 1) {

      // Write the header
      fprintf(filestream, "###KPFITHEADER###\n");
      fprintf(filestream, "%.6e\n", kpfield);
      fprintf(filestream, "%.6e\n", kpinject);
      fprintf(filestream, "%d\n", setmodelres[2]);
      fprintf(filestream, "###END###\n");

      // Write the best ages
      fprintf(filestream, "###KPBESTAGE###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", kpbestage[i]);
      }
      fprintf(filestream, "###END###\n");

      // Write the best normalizations
      fprintf(filestream, "###KPBESTNORM###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", kpbestnorm[i]);
      }
      fprintf(filestream, "###END###\n");

      // Write the chi squared values
      fprintf(filestream, "###KPCHISQUARED###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", kpchisquared[i]);
      }
      fprintf(filestream, "###END###\n");

      // Write the model alpha
      fprintf(filestream, "###KPMODELALPHA###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", kpmodelalpha[i]);
      }
      fprintf(filestream, "###END###\n");

      // Write the model flux
      fprintf(filestream, "###KPMODELFLUX###\n");
      for(i=1; i<=regnumber; i++) {
	for(j=0; j<=setmodelres[2]; j++) {
	  fprintf(filestream, "%d,%d,%.12e\n", i, j, kpmodelflux[i][j]);
	}
      }
      fprintf(filestream, "###END###\n");

      // Write the model age errors
      fprintf(filestream, "###KPAGEERRORSPLUS###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", kpageerrorsplus[i]);
      }
      fprintf(filestream, "###END###\n");

      fprintf(filestream, "###KPAGEERRORSMINUS###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", kpageerrorsminus[i]);
      }
      fprintf(filestream, "###END###\n");

    }


    // Output the JPTRIB model arrays if it has been fitted
    if (jptribmodelset == 1) {

      // Write the header
      fprintf(filestream, "###JPTRIBFITHEADER###\n");
      fprintf(filestream, "%.6e\n", jptribfield);
      fprintf(filestream, "%.6e\n", jptribinject);
      fprintf(filestream, "%d\n", setmodelres[3]);
      fprintf(filestream, "###END###\n");

      // Write the best ages
      fprintf(filestream, "###JPTRIBBESTAGE###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", jptribbestage[i]);
      }
      fprintf(filestream, "###END###\n");

      // Write the best normalizations
      fprintf(filestream, "###JPTRIBBESTNORM###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", jptribbestnorm[i]);
      }
      fprintf(filestream, "###END###\n");

      // Write the chi squared values
      fprintf(filestream, "###JPTRIBCHISQUARED###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", jptribchisquared[i]);
      }
      fprintf(filestream, "###END###\n");

      // Write the model alpha
      fprintf(filestream, "###JPTRIBMODELALPHA###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", jptribmodelalpha[i]);
      }
      fprintf(filestream, "###END###\n");

      // Write the model flux
      fprintf(filestream, "###JPTRIBMODELFLUX###\n");
      for(i=1; i<=regnumber; i++) {
	for(j=0; j<=setmodelres[3]; j++) {
	  fprintf(filestream, "%d,%d,%.12e\n", i, j, jptribmodelflux[i][j]);
	}
      }
      fprintf(filestream, "###END###\n");

      // Write the model age errors
      fprintf(filestream, "###JPTRIBBLEAGEERRORSPLUS###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", jptribbleageerrorsplus[i]);
      }
      fprintf(filestream, "###END###\n");

      fprintf(filestream, "###JPTRIBBLEAGEERRORSMINUS###\n");
      for(i=1; i<=regnumber; i++) {
	fprintf(filestream, "%.12e\n", jptribbleageerrorsminus[i]);
      }
      fprintf(filestream, "###END###\n");

    }

  }
  

  fprintf(filestream, "###ENDOFIMPORTFILE###\n");

  fclose(filestream);

  printf("\nData exported to %s\n\n", filename);

  return 0;

}
