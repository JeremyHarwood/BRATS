/*
List of commands available for use and their descriptions
*/

#include <stdio.h>

int help() {

  // setup the variable and memory

  int help_mainloop = 1;
  int tmp_model = 0;
  int help_cmdnum = -1;
  char help_cmdtxt[1024], *tmp_helpchar, *cmdbuffer;

  printf("========================================================================\n");
  printf("                           MAIN COMMAND LIST                            \n");
  printf("------------------------------------------------------------------------\n");
  printf("         load          |       setregions       |        fluxmap        \n");
  printf("       specindex       |     setsingleregion    |       findinject      \n");
  printf("      plotspecindex    |      fixedregions      |       scaleflux       \n");
  printf("       fitjpmodel      |       fitkpmodel       |      fitjptribble     \n");
  printf("       specagemap      |      chisquaredmap     |        errormap       \n");
  printf("       plotjpmodel     |       plotkpmodel      |      plotjptribble    \n");
  printf("      combineimages    |        diffmap         |       resizeimage     \n");
  printf("      fitintegrated    |      mapfrommodel      |        fitpoly        \n");
  printf("      plotmodelobs     |        plotflux        |     specindexerrors   \n");
  printf("      plotcurvespec    |        contours        |       colourcolour    \n");
  printf("       fullexport      |       fullimport       |       exportdata      \n");
  printf("        export         |        rmsnoise        |     specchisquared    \n");
  printf("      conflevels       |       fitcimodel       |        fitcioff       \n");
  printf("     injectionmap      |     injectionchimap    |       plotcimodel     \n");
  printf("       plotcioff       |                        |                       \n");
  printf("       (basics)        |       (plotting)       |       (modeldata)     \n");

 
  /*
    printf("\n");
    printf("                       |                        |                       \n");
  */
  
  printf("\nWelcome to the BRATS help terminal. Type 'return' at any time to return to the main command prompt.\n\n");

  printf("Typing any of the main commands listed above into the help terminal will give a brief description of the task, along with any sub-commands and parameters which control it's usage.\n\n");

  printf("Those shown in brackets represent a help section which is not in itself a command, but groups together similar sub-commands used to control certain functions within BRATS e.g. typing 'plotting' will describe all of the commands which control how BRATS creates plots, maps and images.\n\n");

  printf("To view this list again type 'tasks' into the help terminal. For additional support, please visit the BRATS website.\n\n");

  ///*
  
  while (help_mainloop == 1) {

    // Set the command number to 0
    help_cmdnum = 0;
   
    // Print the command promt and capture any input
    cmdbuffer = readline("BRATS_HELP: ");

    // If it's not blank, add to the history
    if (cmdbuffer[0] != 0) {
      add_history(cmdbuffer);
    }

    strcpy(help_cmdtxt, cmdbuffer);

    // Format out any return trailing white space
    chomp(help_cmdtxt);


    for (tmp_helpchar = help_cmdtxt; *tmp_helpchar != '\0'; ++tmp_helpchar) {
      *tmp_helpchar = tolower(*tmp_helpchar);
    }

    if (strlen(help_cmdtxt) == 0) {
      help_cmdnum = 0;
    }
    else if ((strstr(help_cmdtxt, "return") != NULL) && (strlen(help_cmdtxt) == 6)) {
      help_cmdnum = 1;
    }
    else if ( ((strstr(help_cmdtxt, "exit") != NULL) || (strstr(help_cmdtxt, "quit") != NULL)) && (strlen(help_cmdtxt) == 4)) {
      help_cmdnum = 2;
    }
    else if ((strstr(help_cmdtxt, "tasks") != NULL) && (strlen(help_cmdtxt) == 5)) {
      help_cmdnum = 3;
    }
    else if ((strstr(help_cmdtxt, "basics") != NULL) && (strlen(help_cmdtxt) == 6)) {
      help_cmdnum = 4;
    }
    else if ((strstr(help_cmdtxt, "load") != NULL) && (strlen(help_cmdtxt) == 4)) {
      help_cmdnum = 5;
    }
    else if ((strstr(help_cmdtxt, "setregions") != NULL) && (strlen(help_cmdtxt) == 10)) {
      help_cmdnum = 6;
    }
    else if ((strstr(help_cmdtxt, "specindex") != NULL) && (strlen(help_cmdtxt) == 9)) {
      help_cmdnum = 7;
    }
    else if ((strstr(help_cmdtxt, "fitjpmodel") != NULL) && (strlen(help_cmdtxt) == 10)) {
      tmp_model=1;
      help_cmdnum = 8;
    }
    else if ((strstr(help_cmdtxt, "fitkpmodel") != NULL) && (strlen(help_cmdtxt) == 10)) {
      tmp_model=2;
      help_cmdnum = 8;
    }
    else if ((strstr(help_cmdtxt, "fitjptribble") != NULL) && (strlen(help_cmdtxt) == 12)) {
      tmp_model=3;
      help_cmdnum = 8;
    }
    else if ((strstr(help_cmdtxt, "fitintegrated") != NULL) && (strlen(help_cmdtxt) == 13)) {
      help_cmdnum = 9;
    }
    else if ((strstr(help_cmdtxt, "mapfrommodel") != NULL) && (strlen(help_cmdtxt) == 12)) {
      help_cmdnum = 10;
    }
    else if ((strstr(help_cmdtxt, "diffmap") != NULL) && (strlen(help_cmdtxt) == 7)) {
      help_cmdnum = 11;
    }
    else if ((strstr(help_cmdtxt, "findinject") != NULL) && (strlen(help_cmdtxt) == 10)) {
      help_cmdnum = 12;
    }
    else if ((strstr(help_cmdtxt, "specagemap") != NULL) && (strlen(help_cmdtxt) == 10)) {
      help_cmdnum = 13;
    }
    else if ((strstr(help_cmdtxt, "chisquaredmap") != NULL) && (strlen(help_cmdtxt) == 13)) {
      help_cmdnum = 14;
    }
    else if ((strstr(help_cmdtxt, "errormap") != NULL) && (strlen(help_cmdtxt) == 8)) {
      help_cmdnum = 15;
    }
    else if ((strstr(help_cmdtxt, "plotmodelobs") != NULL) && (strlen(help_cmdtxt) == 12)) {
      help_cmdnum = 16;
    }
    else if ((strstr(help_cmdtxt, "plotflux") != NULL) && (strlen(help_cmdtxt) == 8)) {
      help_cmdnum = 17;
    }
    else if ((strstr(help_cmdtxt, "plotspecindex") != NULL) && (strlen(help_cmdtxt) == 13)) {
      help_cmdnum = 18;
    }
    else if ((strstr(help_cmdtxt, "plotcurvespec") != NULL) && (strlen(help_cmdtxt) == 13)) {
      help_cmdnum = 19;
    }
    else if ((strstr(help_cmdtxt, "fluxmap") != NULL) && (strlen(help_cmdtxt) == 7)) {
      help_cmdnum = 20;
    }
    else if ((strstr(help_cmdtxt, "plotting") != NULL) && (strlen(help_cmdtxt) == 8)) {
      help_cmdnum = 21;
    }
    else if ((strstr(help_cmdtxt, "contours") != NULL) && (strlen(help_cmdtxt) == 8)) {
      help_cmdnum = 22;
    }
    else if ((strstr(help_cmdtxt, "export") != NULL) && (strlen(help_cmdtxt) == 6)) {
      help_cmdnum = 23;
    }
    else if ((strstr(help_cmdtxt, "exportdata") != NULL) && (strlen(help_cmdtxt) == 10)) {
      help_cmdnum = 24;
    }
    else if ((strstr(help_cmdtxt, "fullexport") != NULL) && (strlen(help_cmdtxt) == 10)) {
      help_cmdnum = 25;
    }
    else if ((strstr(help_cmdtxt, "fullimport") != NULL) && (strlen(help_cmdtxt) == 10)) {
      help_cmdnum = 26;
    }
    else if ((strstr(help_cmdtxt, "combineimages") != NULL) && (strlen(help_cmdtxt) == 13)) {
      help_cmdnum = 27;
    }
    else if ((strstr(help_cmdtxt, "fixedregions") != NULL) && (strlen(help_cmdtxt) == 12)) {
      help_cmdnum = 28;
    }
    else if ((strstr(help_cmdtxt, "colourcolour") != NULL) && (strlen(help_cmdtxt) == 12)) {
      help_cmdnum = 29;
    }
    else if ((strstr(help_cmdtxt, "colorcolor") != NULL) && (strlen(help_cmdtxt) == 12)) {
      help_cmdnum = 29;
    }
    else if ((strstr(help_cmdtxt, "fitpoly") != NULL) && (strlen(help_cmdtxt) == 7)) {
      help_cmdnum = 30;
    }
    else if ((strstr(help_cmdtxt, "plotjpmodel") != NULL) && (strlen(help_cmdtxt) == 11)) {
      help_cmdnum = 31;
    }
    else if ((strstr(help_cmdtxt, "plotkpmodel") != NULL) && (strlen(help_cmdtxt) == 11)) {
      help_cmdnum = 32;
    }
    else if ((strstr(help_cmdtxt, "plotjptribble") != NULL) && (strlen(help_cmdtxt) == 13)) {
      help_cmdnum = 33;
    }
    else if ((strstr(help_cmdtxt, "resizeimage") != NULL) && (strlen(help_cmdtxt) == 11)) {
      help_cmdnum = 34;
    }
    else if ((strstr(help_cmdtxt, "scaleflux") != NULL) && (strlen(help_cmdtxt) == 9)) {
      help_cmdnum = 35;
    }
    else if ((strstr(help_cmdtxt, "setsingleregion") != NULL) && (strlen(help_cmdtxt) == 15)) {
      help_cmdnum = 36;
    }
    else if ((strstr(help_cmdtxt, "specchisquared") != NULL) && (strlen(help_cmdtxt) == 14)) {
      help_cmdnum = 37;
    }
    else if ((strstr(help_cmdtxt, "rmsnoise") != NULL) && (strlen(help_cmdtxt) == 8)) {
      help_cmdnum = 38;
    }
    else if ((strstr(help_cmdtxt, "modeldata") != NULL) && (strlen(help_cmdtxt) == 9)) {
      help_cmdnum = 39;
    }
    else if ((strstr(help_cmdtxt, "conflevels") != NULL) && (strlen(help_cmdtxt) == 10)) {
      help_cmdnum = 40;
    }
    else if ((strstr(help_cmdtxt, "fitcimodel") != NULL) && (strlen(help_cmdtxt) == 10)) {
      tmp_model = 4;
      help_cmdnum = 41;
    }
    else if ((strstr(help_cmdtxt, "fitcioff") != NULL) && (strlen(help_cmdtxt) == 8)) {
      tmp_model = 5;
      help_cmdnum = 41;
    }
    else if ((strstr(help_cmdtxt, "injectionmap") != NULL) && (strlen(help_cmdtxt) == 12)) {
      help_cmdnum = 42;
    }
    else if ((strstr(help_cmdtxt, "injectionchimap") != NULL) && (strlen(help_cmdtxt) == 15)) {
      help_cmdnum = 43;
    }
    else if ((strstr(help_cmdtxt, "specindexerrors") != NULL) && (strlen(help_cmdtxt) == 15)) {
      help_cmdnum = 44;
    }
    else if ((strstr(help_cmdtxt, "plotcimodel") != NULL) && (strlen(help_cmdtxt) == 11)) {
      help_cmdnum = 45;
    }
    else if ((strstr(help_cmdtxt, "plotcioff") != NULL) && (strlen(help_cmdtxt) == 9)) {
      help_cmdnum = 46;
    }
    else {
      help_cmdnum = 999;
    }


    switch (help_cmdnum) {

    case 0:
      ;
      break;


    case 1:
      return 0;
      break;


    case 2:
      printf("Please return to the main commmand prompt using the 'return' command before exiting BRATS\n");
      break;


    case 3:

      printf("\n");
      printf("========================================================================\n");
      printf("                           MAIN COMMAND LIST                            \n");
      printf("------------------------------------------------------------------------\n");
      printf("         load          |       setregions       |        fluxmap        \n");
      printf("       specindex       |     setsingleregion    |       findinject      \n");
      printf("      plotspecindex    |      fixedregions      |       scaleflux       \n");
      printf("       fitjpmodel      |       fitkpmodel       |      fitjptribble     \n");
      printf("       specagemap      |      chisquaredmap     |        errormap       \n");
      printf("       plotjpmodel     |       plotkpmodel      |      plotjptribble    \n");
      printf("      combineimages    |        diffmap         |       resizeimage     \n");
      printf("      fitintegrated    |      mapfrommodel      |        fitpoly        \n");
      printf("      plotmodelobs     |        plotflux        |     specindexerrors   \n");
      printf("      plotcurvespec    |        contours        |       colourcolour    \n");
      printf("       fullexport      |       fullimport       |       exportdata      \n");
      printf("        export         |        rmsnoise        |     specchisquared    \n");
      printf("      conflevels       |       fitcimodel       |        fitcioff       \n");
      printf("     injectionmap      |     injectionchimap    |       plotcimodel     \n");
      printf("       plotcioff       |                        |                       \n");
      printf("       (basics)        |       (plotting)       |       (modeldata)     \n");
      printf("\n");

      break;


    case 4:

      printf("\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("    brats -version     | (When starting BRATS) Display the version information.                        \n");
      printf("    brats -help        | (When starting BRATS) Display the basics statup help.                         \n");
      printf("    quit               | Quit BRATS and return to the system terminal.                                 \n");
      printf("    esc                | Escape to the command screen from (most) tasks requiring a text based input.  \n");
      printf("    -1                 | Escape to the command screen from (most) tasks requiring a numeric input.     \n");
      printf("    ls                 | List the contents of the current directory. Can also used where tasks require \n");
      printf("                       | a text based input e.g. a directory location.                                 \n");
      printf("    shell              | Enter shell command prompt mode which allows the use of multi word system     \n");
      printf("                       | commands e.g. mkdir mynewimagesfolder                                         \n");
      printf("    list               | Output to the terminal all of the currently loaded datasets.                  \n");
      printf("    props              | Output to terminal the basic properties of a dataset.                         \n");
      printf("    zoom               | Set the zoom level to use when outputting maps etc.                           \n");
      printf("    shiftimage         | Shift the centre of the output image by a give number of pixels in the X and  \n");
      printf("                       | Y direction e.g. for sources not at the pointing centre.                      \n");
      printf("    targetname         | Rename the target source/region of a given data set.                          \n");
      printf("\n");


      printf("\n");

      break;


    case 5:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        load           | Loads a new dataset into BRATS .                                              \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       sigma           | Change the initial source detection level. Calculated using the RMS of the    \n");
      printf("                       | background region with no on-source multiplier applied (DEFAULT = 5).         \n");
      printf("       casadata        | Select which reduction type should be defaulted to if load cannot automatica- \n");
      printf("                       | lly determine the header type (DEFAULT = CASA).                               \n");
      printf("       telescope       | Change the default flux calibration errors set when loading a map to a given  \n");
      printf("                       | telescope, or, a fixed user value (DEFAULT = JVLA).                           \n");
      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("     directoryname     | Location of directory containing images to be loaded.                         \n");
      printf("     backgroundregion  | Location of the DS9 region file defining the area in which the off-source     \n");
      printf("                       | background noise should be calculated.                                        \n");
      printf("     sourceregion      | (Optional, but recommended) Location of the DS9 region file encompassing the  \n");
      printf("                       | target source.                                                                \n");
      printf("     redshift          | Redshift of the target source. Required for spectral age fitting. For non-    \n");
      printf("                       | spectral age analysis where the redshift is not known, any reasonable dummy   \n");
      printf("                       | value can be entered.                                                         \n");

      printf("\n");

      break;


    case 6:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       setregions      | Determines and sets adaptive regions for a given data set.                    \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      signaltonoise    | Sets the signal-to-noise level that a region must reach. A value of 1 creates \n");
      printf("                       | regions on a pixel-by-pixel basis (DEFAULT = 1).                              \n");
      printf("      sigma            | Change the initial source detection level. Calculated using the RMS of the    \n");
      printf("                       | background region with no on-source multiplier applied (DEFAULT = 5).         \n");
      printf("      searcharea       | Sets the maximum area a region can be before it is flagged (DEFAULT = 1).     \n");
      printf("      fluxcalerror     | Sets the flux calibration error for a single or range of maps and datasets.   \n");
      printf("                       | Automatically determined for a (limited) number or radio telescopes.          \n");
      printf("      viewerrors       | Displays the currently set values for the flux calibration errors.            \n");
      printf("      onsource         | Change the value by which the RMS is multiplied for the on-source noise       \n");
      printf("                       | (DEFAULT = 3).                                                                \n");
      printf("      hotpixels        | Sets the limit at which pixels will be considered hot and cold as a decimal   \n");
      printf("                       | fraction of the surrounding pixels flux (DEFAULT = 0.2).                      \n");
      printf("      maptomap         | Set the limit at which map-to-map variations will cause flagging as a decimal \n");
      printf("                       | of the flux (DEFAULT = -1 [OFF]).                                             \n");

      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       dataset         | Select for which data set the regions should be determined.                   \n");

      printf("\n");

      break;


    case 7:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       specindex       | Calculate and map the spectral indices by region of a given data set.         \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("   specindexcalctype   | Set the method used when fitting spectral indices to a data set (DEFAULT = 2  \n");
      printf("                       | [weighted GSL least squares fit]).                                            \n");
      printf("     forceerrortype    | Forces the standard deviation method to be used for errors for data sets      \n");
      printf("                       | containing only 2 frequencies (DEFAULT = OFF).                                \n");
      printf("      printindex       | Turn on and off the output of spectral index values to the terminal window as \n");
      printf("                       | they are calculated (DEFAULT = OFF).                                          \n");
      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the spectral indices should be determined.          \n");

      printf("\n");
      
      break;


    case 8:
      if (tmp_model == 1) {
	printf("\n");
	printf("Main Command:\n");
	printf("=======================================================================================================\n");
	printf("        Command        | Description                                                                   \n");
	printf("-------------------------------------------------------------------------------------------------------\n");
	printf("       fitjpmodel      | Fit the JP model of spectral ageing to a given data set.                      \n");

      }
      else if (tmp_model == 2) {
	printf("\n");
	printf("Main Command:\n");
	printf("=======================================================================================================\n");
	printf("        Command        | Description                                                                   \n");
	printf("-------------------------------------------------------------------------------------------------------\n");
	printf("       fitkpmodel      | Fit the KP model of spectral ageing to a given data set.                      \n");
      }
      else if (tmp_model == 3) {
	printf("\n");
	printf("Main Command:\n");
	printf("=======================================================================================================\n");
	printf("        Command        | Description                                                                   \n");
	printf("-------------------------------------------------------------------------------------------------------\n");
	printf("      fitjptribble     | Fit the Tribble model of spectral ageing to a given data set.                 \n");
      }
      else {
	printf("***WARNING*** Unable to determine which model is help requires! Please contact brats@askanastronomer if this problem persists\n");
	printf("Displaying common commands and parameters for all single injection models instead... \n");
      }

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       bfield          | Set the magnetic field strength to use in model fitting (DEFAULT =  1e-9 T).  \n");
      printf("    injectionindex     | Set the injection index to use in model fitting (DEFAULT = 0.6)               \n");
      printf("      modelpower       | Alternative to ‘injectionindex’. Sets the injection index to use in model     \n");
      printf("                       | fitting based on the model power, where α inj = (P − 1)/2 (DEFAULT = 2.2).    \n");
      printf("       myears          | Set the maximum spectral age in megayears which for the fitting should be     \n");
      printf("                       | attempted (DEFAULT = 10).                                                     \n");
      printf("      minmyears        | Set the minimum spectral age in megayears which for the fitting should be     \n");
      printf("                       | attempted (DEFAULT = 0).                                                      \n");
      printf("       ageres          | Change the resolution of the ages attempted when fitting models e.g. 10 gives \n");
      printf("                       | a step size 0.1 for the first level of fitting between 1 and 10 myears        \n");
      printf("                       | (DEFAULT = 10).                                                               \n");
      printf("       levels          | The number of (age) levels deep to go when model fitting. 3 is recommended    \n");
      printf("                       | for testing and initial results, 5 for final results. Note setting this to    \n");
      printf("                       | high values drastically increases computation time! (DEFAULT = 3).            \n");
      printf("        gmin           | Set the minimum value of gamma for model fitting (DEFAULT = 10).              \n");
      printf("        gmax           | Set the maximum value of gamma for model fitting (DEFAULT = 1000000)          \n");
      printf("     printresults      | Turn on and off the printing of individual results when model fitting. The    \n");
      printf("                       | The final summary will still be printed on completion of a model fit (DEFAULT \n");
      printf("                       | = OFF).                                                                       \n");
      printf("     printreject       | Turn on and off the output whether it should be stated that a model should be \n");
      printf("                       | rejected if the average chi-squared is above the 90 per cent confidence inte- \n");
      printf("                       | rval. This only affects the print out to the terminal in the summary, not the \n");
      printf("                       | fitting itself. Note if there are a number of regions where you expect a bad  \n");
      printf("                       | fit e.g. from dynamic range issues, the average may not be a suitable value   \n");
      printf("                       | to use (DEFAULT = OFF).                                                       \n");
      printf("     suppressconf      | Stops the chi-squared confidence level tables from being automatically output.\n");
      printf("                       | This is useful for avoiding errors where the number of DoF of a data set is   \n");
      printf("                       | larger than can be handled by the CDF function (see known issues).            \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the spectral age model should be fitted.            \n");

      printf("\n");

      tmp_model = 0;

      break;
      

    case 9:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("     fitintegrated     | Fit a model to integrated flux values. Note integrated fitting values are     \n");
      printf("                       | output in text format only (Section 3.4.2 of the BRATS manual for details).   \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       bfield          | Set the magnetic field strength to use in model fitting (DEFAULT =  1e-9 T).  \n");
      printf("    injectionindex     | Set the injection index to use in model fitting (DEFAULT = 0.6)               \n");
      printf("      modelpower       | Alternative to ‘injectionindex’. Sets the injection index to use in model     \n");
      printf("                       | fitting based on the model power, where α inj = (P − 1)/2 (DEFAULT = 2.2).    \n");
      printf("       myears          | Set the maximum spectral age in megayears which for the  fitting should be    \n");
      printf("                       | attempted (DEFAULT = 10).                                                     \n");
      printf("       ageres          | Change the resolution of the ages attempted when fitting models e.g. 10 gives \n");
      printf("                       | a step size 0.1 for the first level of fitting between 1 and 10 myears        \n");
      printf("                       | (DEFAULT = 10).                                                               \n");
      printf("       levels          | The number of (age) levels deep to go when model fitting. 3 is recommended    \n");
      printf("                       | for testing and initial results, 5 for final results. Note setting this to    \n");
      printf("                       | high values drastically increases computation time! (DEFAULT = 3).            \n");
      printf("        gmin           | Set the minimum value of gamma for model fitting (DEFAULT = 10).              \n");
      printf("        gmax           | Set the maximum value of gamma for model fitting (DEFAULT = 1000000)          \n");
      printf("     printresults      | Turn on and off the printing of individual results when model fitting. The    \n");
      printf("                       | The final summary will still be printed on completion of a model fit (DEFAULT \n");
      printf("                       | = OFF).                                                                       \n");
      printf("      printreject      | Turn on and off the output whether it should be stated that a model should be \n");
      printf("                       | rejected if the average chi-squared is above the 90 per cent confidence inte- \n");
      printf("                       | rval. This only affects the print out to the terminal in the summary, not the \n");
      printf("                       | fitting itself. Note if there are a number of regions where you expect a bad  \n");
      printf("                       | fit e.g. from dynamic range issues, the average may not be a suitable value   \n");
      printf("                       | to use (DEFAULT = OFF).                                                       \n");
 
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       model           | Select which model to be fitted.                                              \n");
      printf("       redshift        | Enter the redshift of the source. Note that for standard non-integrated fitt- \n");
      printf("                       | ting, redshift is set during the ‘load’ command.                              \n");
      printf("  numberoffrequencies  | Number of data points to be fitted to.                                        \n");
      printf("       frequency       | Enter the frequency of the nth data point.                                    \n");
      printf("       flux            | Enter the flux of the nth data point.                                         \n");
      printf("       error           | Enter the error of the nth data point.                                        \n");

      printf("\n");

      break;


    case 10:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      mapfrommodel     | Creates a simulated radio map in FITS format at any user defined frequency    \n");
      printf("                       | given any spectral ageing model. Model fitting of a real data set must first  \n");
      printf("                       | be performed to obtain the required spectral ages and normalisations.         \n");


      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        None           |                                                                               \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select which data set should provide the required model fitting data.         \n");
      printf("        frequency      | Set to what frequency the spectral ageing model should be extrapolated /      \n");
      printf("                       | iterpolated.                                                                  \n");
      printf("        name           | Set a name from the resulting reconstructed image file.                       \n");

      printf("\n");

      break;


    case 11:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        diffmap        | Subtracts the flux of one FITS file from another and outputs a third FITS     \n");
      printf("                       | image. Subtraction can either be full map or confined to the regions of a     \n");
      printf("                       | data set.                                                                     \n");


      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        None           |                                                                               \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       sourcemap       | Image which should be subtracted from.                                        \n");
      printf("     subtractionmap    | Image to subtract.                                                            \n");
      printf("    restricttoregions  | Sets whether to restrict the subtraction to the regions of a given data set.  \n");

      printf("\n");

      break;


    case 12:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       findinject      | Run a model fit (or multiple model fits) with varying injection indices in an \n");
      printf("                       | attempt to find the best fitting injection index. Warning: This function can  \n");
      printf("                       | take a VERY long time to run, especially if multiple models over multiple     \n");
      printf("                       | datasets are being attempted.                                                \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       mininject       | Set the minimum injection index be tested with findinject (DEFAULT = 0.5).    \n");
      printf("       maxinject       | Set the maximum injection index be tested with findinject (DEFAULT = 1.0).    \n");
      printf("    injectintervals    | Set the number of intervals between mininject and maxinject to be tested by   \n");
      printf("                       | findinject.  n + 1 intevals are tested to find the best fit e.g. with         \n");
      printf("                       | mininject 0.5 and and maxinject 1.0, n = 10 will result in steps of 0.5       \n");
      printf("                       | (DEFAULT = 10).                                                               \n");
      printf("       bfield          | Set the magnetic field strength to use in model fitting (DEFAULT = 1e-9 T).   \n");
      printf("       modelpower      | Alternative to ‘injectionindex’. Sets the injection index to use in model     \n");
      printf("                       | fitting based on the model power, where α inj = (P − 1)/2 (DEFAULT = 2.2).    \n");
      printf("       myears          | Set the maximum spectral age in megayears for which the fitting should be     \n");
      printf("                       | attempted (DEFAULT = 10).                                                     \n");
      printf("       ageres          | Change the resolution of the ages attempted when fitting models e.g. 10 gives \n");
      printf("                       | a step size 0.1 for the first level of fitting between 1 and 10 myears        \n");
      printf("                       | (DEFAULT = 10).                                                               \n");
      printf("       levels          | The number of (age) levels deep to go when model fitting. 3 is recommended    \n");
      printf("                       | for testing and initial results, 5 for final results. Note setting this to    \n");
      printf("                       | high values drastically increases computation time! (DEFAULT = 3).            \n");
      printf("        gmin           | Set the minimum value of gamma for model fitting (DEFAULT = 10).              \n");
      printf("        gmax           | Set the maximum value of gamma for model fitting (DEFAULT = 1000000)          \n");

     
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the injection index should be determined.           \n");
      printf("      ageingmodel      | Select which spectral ageing model should be fitted.                          \n");


      printf("\n");


      break;


    case 13:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      specagemap       | Plot a map of the spectral age as a function of position.                     \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("     alwayszeroage     | Turn on and off fixing the minimum age range to zero for spectral ageing      \n");
      printf("                       | maps. Note: this only adjusts the wedge and colour range, not the data        \n");
      printf("                       | values (DEFAULT = ON).                                                        \n");
      printf("       contours        | Turn on and off the overlay of flux contours.                                 \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the model fitting results should be plotted.        \n");
      printf("         model         | Select for which model the results should be plotted.                         \n");

      printf("\n");

      break;


    case 14:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("     chisquaredmap     | Plot a map of the chi-squared values as a function of position.               \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       contours        | Turn on and off the overlay of flux contours.                                 \n");
      printf("     suppressconf      | Stops the chi-squared confidence level tables from being automatically output.\n");
      printf("                       | This is useful for avoiding errors where the number of DoF of a data set is   \n");
      printf("                       | larger than can be handled by the CDF function (see known issues).            \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the model fitting results should be plotted.        \n");
      printf("         model         | Select for which model the results should be plotted.                         \n");

      printf("\n");

      break;


    case 15:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        errormap       | Plot a map of spectral age error as a function of position.               \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       contours        | Turn on and off the overlay of flux contours.                                 \n");


      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the model fitting results should be plotted.        \n");
      printf("         model         | Select for which model the results should be plotted.                         \n");
      printf("       errortype       | Select whether positive or negative errors should beplotted (positive errors  \n");
      printf("                       | recommended).                                                                 \n");

      printf("\n");

      break;


    case 16:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      plotmodelobs     | Plot the flux values of both the model and observed data as a function of     \n");
      printf("                       | frequency.                                                                    \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      (plotting)       | Standard global plotting parameters (see the 'plotting' section of help).     \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the results should be plotted.                      \n");
      printf("        model          | Select for which model the results should be plotted.                         \n");

      printf("\n");

      break;


    case 17:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        plotflux       | Plot region flux and a function of frequency.                                 \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      (plotting)       | Standard global plotting parameters (see the 'plotting' section of help).     \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the results should be plotted.                      \n");

      printf("\n");

      break;


    case 18:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("     plotspecindex     | Plot the spectral index fits for each region.                                 \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      (plotting)       | Standard global plotting parameters (see the 'plotting' section of help).     \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the results should be plotted.                      \n");

      printf("\n");

      break;


    case 19:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      plotcurvespec    | Plot curvature (fitpoly) against spectral index for a given dataset.          \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      (plotting)       | Standard global plotting parameters (see the 'plotting' section of help).     \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the results should be plotted.                      \n");

      printf("\n");

      break;


    case 20:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        fluxmap        | Create a standard flux map by plotting the flux of a source as a function of  \n");
      printf("                       | position. Note only flux above the original flux cut made when loading the    \n");
      printf("                       | data will be plotted.                                                         \n");


      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      (plotting)       | Standard global plotting parameters (see the 'plotting' section of help).     \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       dataset         | Select for which data set the flux should be plotted.                         \n");
      printf("       frequency       | Select which frequency should be plotted.                                     \n");
      printf("       rawvsregions    | Select if the raw data or adaptive regions should be used when creating the   \n");
      printf("                       | flux maps.                                                                    \n");
      printf("       average         | (adaptive regions only) Select if the total or average flux of each regions   \n");
      printf("                       | should be mapped.                                                             \n");

      printf("\n");

      break;


    case 21:

      printf("\n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       border          | Set the border size as a percentage of the total size for mapping             \n");
      printf("                       | (DEFAULT = 3.0).                                                              \n");
      printf("       titles          | Turn on and off titles in maps and plots.                                     \n");
      printf("       labels          | Turn on and off axis labels in maps and plots.                                \n");
      printf("       axismarks       | Turn on and off axis markings in maps and plots                               \n");
      printf("       autoscale       | Turn on and off autoscaling of all plot axis, and set or change user defined  \n");
      printf("                       | values.                                                                       \n");
      printf("       uselogs         | Turn on and off the use of logarithmic values for plots.                      \n");
      printf("       skip            | Set a value of how many regions to skip when plotting (i.e. only plot every   \n");
      printf("                       | nth region).                                                                  \n");
      printf("       symbol          | Change the symbol which is used to mark values when creating plots.           \n");
      printf("       zoom            | Set the zoom level to use when outputting maps etc.                           \n");
      printf("       shiftimage      | Shift the centre of the output image by a give number of pixels in the X and  \n");
      printf("                       | Y direction e.g. for sources not at the pointing centre that require the zoom \n");
      printf("                       | function to be used. NOTE: Currently not applied to the setregions command.   \n");
      printf("     wedgemultiplier   | Multiplier to apply to the maximim value for a given dataset that defines the \n");
      printf("                       | upper limit of the wedge colour. This is useful in avoiding real values from  \n");
      printf("                       | taking the background colour (DEFAULT = 1.15).                                \n");
      printf("       scaletype       | Set the scale to be used when mapping. Options are: Pixels, Arcseconds (from  \n");
      printf("                       | centre of map), Arcseconds in degrees, minutes, seconds format or WCS.        \n");
      printf("                       | (DEFAULT =  WCS).                                                             \n");
      printf("       posmap          | Set which map should be used as reference for co-ordinates when plotting maps \n");
      printf("                       | maps in WCS (DEFAULT = 0).                                                    \n");
      printf("       beam            | Turn on and off the plotting of the beam to maps (DEFAULT = ON).              \n");
      printf("       beampos         | Set the position to plot the beam when making maps (DEFAULT X = 95, Y = 5).  \n");
      printf("       usereduced      | Turn on and off the use of reduced chi-squared values for plotting and        \n");
      printf("                       | mapping (DEFAULT = ON).                                                       \n");
      printf("       largetext       | Turn on and off the use of large text format for plotting and mapping         \n");
      printf("                       | (DEFAULT = OFF).                                                              \n");
      printf("       colourscheme    | (Mapping only) Change the colour scheme used for mapping (DEFAULT = RAINBOW). \n");
      
      printf("\n");

      break;


    case 22:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        contours       | Turn on and off the printing of contours (DEFAULT = OFF).                     \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("     contourlevels     | Set the number of contours to plot between min and max values (DEFAULT = 8).  \n");
      printf("     contourcolour     | Set the colour of the contour lines (DEFAULT = BLACK).                        \n");
      printf("     contourlinestyle  | Set the style of contour lines (DEFAULT = SOLID).                             \n");
      printf("     contourlogs       | Set whether to use linear or logarithmic spacings between contours (DEFAULT = \n");
      printf("                       | LOG).                                                                         \n");
      printf("     firstcontour      | The level of the first contour to plot (above the minimum) (DEFAULT= 1).      \n");
      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       frequency       | Set what frequency map should be used to determine the contours to be         \n");
      printf("                       | overlaid. An additional option will automatically appear for this selection   \n");
      printf("                       | to be made in each main command where contours are relevant e.g. the          \n");
      printf("                       | ‘specagemap’ command.                                                         \n");

      printf("\n");

      break;


    case 23:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        export         | Turn on and off the exporting of data (DEFAULT = OFF).                        \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       imageloc        | Set the location for images to be exported (DEFAULT = ./images).              \n");
      printf("       imagetype       | Set the format of exported images (DEFAULT = PNG).                            \n");
      printf("       exportasfits    | Exports maps in FITS format. Overrides imagetype. (DEFAULT = OFF).            \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("         None          |                                                                               \n");

      printf("\n");

      break;


    case 24:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       exportdata      | Exports numerical data such as spectral ages, errors and injection indices    \n");
      printf("                       | to a text file.                                                               \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       dataloc         | Set the location for data to be exported (DEFAULT = ./data).                  \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       dataset         | Select which data set for which the data should be exported.                  \n");
      printf("   datatobeexported    | Select which data should be exported.                                         \n");
      printf("       filename        | Enter a name to be appended to the end of the standard naming format.         \n");

      printf("\n");

      break;


    case 25:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       fullexport      | Export a full dataset, including any model fitting, in .brats format.         \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("   exportcompression   | Set exporting of a full dataset to either full or compressed.  Note that      \n");
      printf("                       | compressed is recommended in nearly all circumstances. (DEFAULT = Compressed).\n");
      printf("       dataloc         | Set the location for data to be exported (DEFAULT = ./data).                  \n");
      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       dataset         | Select which data set should be exported.                                     \n");

      printf("\n");

      break;


    case 26:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       fullimport      | Import a full dataset from a .brats format file.                              \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        None           |                                                                               \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      importfile       | Select the .brats file from which the data should be imported.                \n");

      printf("\n");

      break;


    case 27:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("     combineimages     | Combines a set of radio maps to a common (user defined) frequency in the      \n");
      printf("                       | image plane.                                                                  \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        None           |                                                                               \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       imagefolder     | Select a folder containing the images to be combined.                         \n");
      printf("       frequency       | Set the frequency at which the combined map should be made.                   \n");
      printf("     powerlawscaling   | Set the power law that the images should be scaled by.                        \n");

      printf("\n");

      break;


    case 28:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      fixedregions     | Applies the regions of one data set to another using pixel coordinates.       \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        None           |                                                                               \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       regionsfrom     | Select which data set’s regions should be used as reference.                  \n");
      printf("       regionsto       | Select which data set the regions should be applied to.                       \n");

      printf("\n");

      break;
      

    case 29:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      colourcolour     | Create a colour-colour plot of a given data set.                              \n");
      printf("      colorcolor       | Alternative US spelling for the colourcolour command.                         \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        None           |                                                                               \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       dataset         | Select which data set a colour-colour plot should be produced for.            \n");
      printf("    lowestfrequency    | Select a map to be used as the lowest frequency data point.                   \n");
      printf("     midfrequency      | Select a map to be used as the mid frequency data point.                      \n");
      printf("    highfrequency      | Select a map to be used as the high frequency data point.                     \n");

      printf("\n");

      break;


    case 30:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        fitpoly        | Fit an nth order polynomial to the regions of a given data set.               \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       setpoly         | Define the order of polynomial to use for fitting (DEFAULT = 2).              \n");
      printf("       curvecon        | Switch the sign convention used when fitting curves (DEFAULT = OFF).          \n");
      printf("       plotres         | Set the number of points to use when plotting the curve (DEFAULT = 1000,      \n");
      printf("                       | MINIMUM 10).                                                                  \n");

      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       dataset         | Select which data set the polynomial fit should be performed on.              \n");

      printf("\n");


      break;


    case 31:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       plotjpmodel     | Plot an example JP model between ‘minmodelfreq’ and ‘maxmodelfreq’ for an     \n");
      printf("                       | arbitrary normalisation.                                                      \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("     modelredshift     | Set the redshift for example spectral ageing models (DEFAULT = 0.2).          \n");
      printf("     minmodelfreq      | Change the minimum frequency in Hz for example spectral ageing models         \n");
      printf("                       | (DEFAULT = 1e7 Hz).                                                           \n");
      printf("     maxmodelfreq      | Change the minimum frequency in Hz for example spectral ageing models         \n");
      printf("                       | (DEFAULT = 1e13 Hz).                                                          \n");
      printf("     modelmyears       | Set the maximum age to plot when outputting models. Model will output at      \n");
      printf("                       | intervals of 1 megayear (DEFAULT = 10 Myr).                                   \n");
      printf("    minmodelmyears     | Set the minimum age to plot when outputting models (DEFAULT = 0 Myr).         \n");
      printf("      modelres         | Set the number of frequency data points to plot for each age (DEFAULT = 100). \n");
      printf("       bfield          | Set the magnetic field strength to use in model fitting (DEFAULT =  1e-9 T).  \n");
      printf("    injectionindex     | Set the injection index to use in model fitting (DEFAULT = 0.6)               \n");
      printf("      modelpower       | Alternative to ‘injectionindex’. Sets the injection index to use in model     \n");
      printf("                       | fitting based on the model power, where α inj = (P − 1)/2 (DEFAULT = 2.2).    \n");
      printf("        gmin           | Set the minimum value of gamma for model fitting (DEFAULT = 10).              \n");
      printf("        gmax           | Set the maximum value of gamma for model fitting (DEFAULT = 1000000)          \n");
      printf("        skip           | Set the value of how many ages to skip when plotting e.g. a value of 10 plots \n");
      printf("                       | models between minmodelmyears and modelmyears at 10 Myr intervals.            \n");

      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("         None          |                                                                               \n");

      printf("\n");

      break;


    case 32:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      plotkpmodel      | Plot an example KP model between ‘minmodelfreq’ and ‘maxmodelfreq’ for an     \n");
      printf("                       | arbitrary normalisation.                                                      \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("     modelredshift     | Set the redshift for example spectral ageing models (DEFAULT = 0.2).          \n");
      printf("     minmodelfreq      | Change the minimum frequency in Hz for example spectral ageing models         \n");
      printf("                       | (DEFAULT = 1e7 Hz).                                                           \n");
      printf("     maxmodelfreq      | Change the minimum frequency in Hz for example spectral ageing models         \n");
      printf("                       | (DEFAULT = 1e13 Hz).                                                          \n");
      printf("     modelmyears       | Set the maximum age to plot when outputting models. Model will output at      \n");
      printf("                       | intervals of 1 megayear (DEFAULT = 10 Myr).                                   \n");
      printf("    minmodelmyears     | Set the minimum age to plot when outputting models (DEFAULT = 0 Myr).         \n");
      printf("      modelres         | Set the number of frequency data points to plot for each age (DEFAULT = 100). \n");
      printf("       bfield          | Set the magnetic field strength to use in model fitting (DEFAULT =  1e-9 T).  \n");
      printf("    injectionindex     | Set the injection index to use in model fitting (DEFAULT = 0.6)               \n");
      printf("      modelpower       | Alternative to ‘injectionindex’. Sets the injection index to use in model     \n");
      printf("                       | fitting based on the model power, where α inj = (P − 1)/2 (DEFAULT = 2.2).    \n");
      printf("        gmin           | Set the minimum value of gamma for model fitting (DEFAULT = 10).              \n");
      printf("        gmax           | Set the maximum value of gamma for model fitting (DEFAULT = 1000000)          \n");
      printf("        skip           | Set the value of how many ages to skip when plotting e.g. a value of 10 plots \n");
      printf("                       | models between minmodelmyears and modelmyears at 10 Myr intervals.            \n");
      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("         None          |                                                                               \n");

      printf("\n");

      break;


    case 33:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("     plotjptribble     | Plot an example Tribble  model between ‘minmodelfreq’ and ‘maxmodelfreq’ for  \n");
      printf("                       | an arbitrary normalisation.                                                   \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("     modelredshift     | Set the redshift for example spectral ageing models (DEFAULT = 0.2).          \n");
      printf("     minmodelfreq      | Change the minimum frequency in Hz for example spectral ageing models         \n");
      printf("                       | (DEFAULT = 1e7 Hz).                                                           \n");
      printf("     maxmodelfreq      | Change the minimum frequency in Hz for example spectral ageing models         \n");
      printf("                       | (DEFAULT = 1e13 Hz).                                                          \n");
      printf("     modelmyears       | Set the maximum age to plot when outputting models. Model will output at      \n");
      printf("                       | intervals of 1 megayear (DEFAULT = 10 Myr).                                   \n");
      printf("    minmodelmyears     | Set the minimum age to plot when outputting models (DEFAULT = 0 Myr).         \n");
      printf("      modelres         | Set the number of frequency data points to plot for each age (DEFAULT = 100). \n");
      printf("       bfield          | Set the magnetic field strength to use in model fitting (DEFAULT =  1e-9 T).  \n");
      printf("    injectionindex     | Set the injection index to use in model fitting (DEFAULT = 0.6)               \n");
      printf("      modelpower       | Alternative to ‘injectionindex’. Sets the injection index to use in model     \n");
      printf("                       | fitting based on the model power, where α inj = (P − 1)/2 (DEFAULT = 2.2).    \n");
      printf("        gmin           | Set the minimum value of gamma for model fitting (DEFAULT = 10).              \n");
      printf("        gmax           | Set the maximum value of gamma for model fitting (DEFAULT = 1000000)          \n");
      printf("        skip           | Set the value of how many ages to skip when plotting e.g. a value of 10 plots \n");
      printf("                       | models between minmodelmyears and modelmyears at 10 Myr intervals.            \n");
      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("         None          |                                                                               \n");

      printf("\n");

      break;


    case 34:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      resizeimage      | Resize a series of FITS images to a new, user defined size.                   \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       imageloc        | Set the location for the resized images to be exported (DEFAULT = ./images).  \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       imagefolder     | Select which folder contains the images to be resized.                        \n");
      printf("       X size          | Size of the new image’s X axis (in pixels).                                   \n");
      printf("       Y size          | Size of the new image’s Y axis (in pixels).                                   \n");

      printf("\n");

      break;


    case 35:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       scaleflux       | Scales the raw flux values by a given scaling factor. Either single scaling   \n");
      printf("                       | value or interpolation / extrapolation of two values can be used for either   \n");
      printf("                       | all or a subsection of maps. Region selection should be re-run after the      \n");
      printf("                       | scaling has been applied.                                                     \n");


      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        None          |                                                                                \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the flux should be scaled.                          \n");
      printf("     mapstobescaled    | Select which maps within the data set should be scaled.                       \n");
      printf("      minfrequency     | (Frequency range only) Select the minimum frequency range to be scaled.       \n");
      printf("      maxfrequency     | (Frequency range only) Select the maximum frequency range to be scaled.       \n");
      printf("       scaletype       | Select the type of scaling to be used.                                        \n");
      printf("     firstfrequency    | (Interpolation / extrapolation only) Enter the frequency of the first scaling \n");
      printf("                       | factor to be used.                                                            \n");
      printf("   firstscalingvalue   | (Interpolation / extrapolation only) Enter the value of first scaling factor  \n");
      printf("                       | to be used.                                                                   \n");
      printf("     secondfrequency   | (Interpolation / extrapolation only) Enter the frequency of the second scal-  \n");
      printf("                       | ing factor to be used.                                                        \n");
      printf("   secondscalingvalue  | (Interpolation / extrapolation only) Enter the value of second scaling factor \n");
      printf("                       | factor to be used.                                                            \n");
      printf("   singlescalingvalue  | (Single value only) Enter the value of the scaling factor to be used.         \n");

      printf("\n");

      break;


    case 36:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("     setsingleregion   | Applies a single region to a selected data set.                               \n");


      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      signaltonoise    | Sets the signal-to-noise level that thr region must reach.                     \n");
      printf("      sigma            | Change the initial source detection level. Calculated using the RMS of the    \n");
      printf("                       | background region with no on-source multiplier applied (DEFAULT = 5).         \n");
      printf("      fluxcalerror     | Sets the flux calibration error for a single or range of maps and datasets.   \n");
      printf("                       | Automatically determined for a (limited) number or radio telescopes.          \n");
      printf("      viewerrors       | Displays the currently set values for the flux calibration errors.            \n");
      printf("      onsource         | Change the value by which the RMS is multiplied for the on-source noise       \n");
      printf("                       | (DEFAULT = 3).                                                                \n");

      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       dataset         | Select for which data set the region should be determined.                       \n");

      printf("\n");

      break;


    case 37:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("     specchisquared    | (Weighted GSL method only) Map the chi-squared values for the spectral index  \n");
      printf("                       | fitting as a function of position.                                            \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      (plotting)       | Standard global plotting parameters (see the 'plotting' section of help).     \n");
      printf("     suppressconf      | Stops the chi-squared confidence level tables from being automatically output.\n");
      printf("                       | This is useful for avoiding errors where the number of DoF of a data set is   \n");
      printf("                       | larger than can be handled by the CDF function (see known issues).            \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the results should be plotted.                      \n");

      printf("\n");

      break;


    case 38:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       rmsnoise        | Manually set the rms noise of a dataset. Can be applied to a single map, all  \n");
      printf("                       | maps, or to a specified frequency range.                                      \n");


      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        None          |                                                                                \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data the rms noise should be applied.                        \n");
      printf("     mapstobescaled    | Select which maps within the data set should be scaled.                       \n");
      printf("      greaterorless    | (Frequency range only) Select whether greater or less than a given frequnecy. \n");
      printf("       frequency       | (Frequency range only) Enter the frequency above/below which should be scaled.\n");
      printf("          rms          | Enter the new rms noise to be applied in Jy/Beam.                             \n");

      printf("\n");

      break;


    case 39:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        jpdata         | Exports example JP model data to a comma delimited textfile.                  \n");
      printf("        kpdata         | Exports example KP model data to a comma delimited textfile.                  \n");
      printf("        tribbledata    | Exports example Tribble model data to a comma delimited textfile.             \n");
      printf("        cidata         | Exports example CI model data to a comma delimited textfile.                  \n");
      printf("        cioffdata      | Exports example CI off model data to a comma delimited textfile.              \n");


      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       dataloc         | Set the location for data to be exported (DEFAULT = ./data).                  \n");
      printf("     modelredshift     | Set the redshift for example spectral ageing models (DEFAULT = 0.2).          \n");
      printf("     minmodelfreq      | Change the minimum frequency in Hz for example spectral ageing models         \n");
      printf("                       | (DEFAULT = 1e7 Hz).                                                           \n");
      printf("     maxmodelfreq      | Change the minimum frequency in Hz for example spectral ageing models         \n");
      printf("                       | (DEFAULT = 1e13 Hz).                                                          \n");
      printf("     modelmyears       | Set the maximum age to plot when outputting models. Model will output at      \n");
      printf("                       | intervals of 1 megayear (DEFAULT = 10 Myr).                                   \n");
      printf("    minmodelmyears     | Set the minimum age to plot when outputting models (DEFAULT = 0 Myr).         \n");
      printf("    dataintervals      | Set the number of frequency data points to export for each age (DEFAULT = 100)\n");
      printf("       bfield          | Set the magnetic field strength to use in model fitting (DEFAULT =  1e-9 T).  \n");
      printf("    injectionindex     | Set the injection index to use in model fitting (DEFAULT = 0.6)               \n");
      printf("      modelpower       | Alternative to ‘injectionindex’. Sets the injection index to use in model     \n");
      printf("                       | fitting based on the model power, where α inj = (P − 1)/2 (DEFAULT = 2.2).    \n");
      printf("        gmin           | Set the minimum value of gamma for model fitting (DEFAULT = 10).              \n");
      printf("        gmax           | Set the maximum value of gamma for model fitting (DEFAULT = 1000000)          \n");
      printf("        gmax           | Set the maximum value of gamma for model fitting (DEFAULT = 1000000)          \n");
      printf("     minmodeloff       | (CI off only) Set the minimum off time to plot when outputting models. Model  \n");
      printf("                       | will output at  intervals of 1 megayear (DEFAULT = 0 Myr).                    \n");
      printf("     maxmodeloff       | (CI off only) Set the maximum off time to output (DEFAULT = 20).              \n");
      printf("     varyoffage        | (CI off only) Set whether the off age should be varied (DEFAULT = ON).        \n");


      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("    exactageselect     | Select whether to output a range of ages between minmodelmyears and           \n");
      printf("                       | modelmyears or an exact age.                                                  \n");
      printf("       exactage        | (Exact age only) Set the age in Myr to be output.                             \n");

      printf("\n");

      break;

    case 40:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       conflevels      | Display a table of confidence levels (both standard and reduced) for a given  \n");
      printf("                       | dataset or user defined number of degrees of freedom.                         \n");


      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        None           |                                                                                \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the confidence levels should be calculated.         \n");
      printf("          dof          | (Manual entry only) Set the number of degrees of freedom that should be used  \n");
      printf("                       | when calculating the confidence levels.                                       \n");

      printf("\n");

      break;

    case 41:
      if (tmp_model == 4) {
	printf("\n");
	printf("Main Command:\n");
	printf("=======================================================================================================\n");
	printf("        Command        | Description                                                                   \n");
	printf("-------------------------------------------------------------------------------------------------------\n");
	printf("       fitcimodel      | Fit the CI model of spectral ageing to a given data set.                      \n");
	
	printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
	printf("=======================================================================================================\n");
	printf("       Command         | Description                                                                   \n");
	printf("-------------------------------------------------------------------------------------------------------\n");
	printf("       bfield          | Set the magnetic field strength to use in model fitting (DEFAULT =  1e-9 T).  \n");
	printf("    injectionindex     | Set the injection index to use in model fitting (DEFAULT = 0.6)               \n");
	printf("      modelpower       | Alternative to ‘injectionindex’. Sets the injection index to use in model     \n");
	printf("                       | fitting based on the model power, where α inj = (P − 1)/2 (DEFAULT = 2.2).    \n");
	printf("       myears          | Set the maximum spectral age (on time) in megayears for which the fitting     \n");
	printf("                       | should be attempted (DEFAULT = 10).                                           \n");
	printf("      minmyears        | Set the minimum spectral age (on time) in megayears for which the fitting     \n");
	printf("                       | should be attempted (DEFAULT = 0).                                            \n");
	printf("       ageres          | Change the resolution of the ages attempted when fitting models e.g. 10 gives \n");
	printf("                       | a step size 0.1 for the first level of fitting between 1 and 10 myears        \n");
	printf("                       | (DEFAULT = 10).                                                               \n");
	printf("       levels          | The number of (age) levels deep to go when model fitting. 3 is recommended    \n");
	printf("                       | for testing and initial results, 5 for final results. Note setting this to    \n");
	printf("                       | high values drastically increases computation time! (DEFAULT = 3).            \n");
	printf("        gmin           | Set the minimum value of gamma for model fitting (DEFAULT = 10).              \n");
	printf("        gmax           | Set the maximum value of gamma for model fitting (DEFAULT = 1000000)          \n");
	printf("     printresults      | Turn on and off the printing of individual results when model fitting. The    \n");
	printf("                       | The final summary will still be printed on completion of a model fit (DEFAULT \n");
	printf("                       | = OFF).                                                                       \n");
	printf("     printreject       | Turn on and off the output whether it should be stated that a model should be \n");
	printf("                       | rejected if the average chi-squared is above the 90 per cent confidence inte- \n");
	printf("                       | rval. This only affects the print out to the terminal in the summary, not the \n");
	printf("                       | fitting itself. Note if there are a number of regions where you expect a bad  \n");
	printf("                       | fit e.g. from dynamic range issues, the average may not be a suitable value   \n");
	printf("                       | to use (DEFAULT = OFF).                                                       \n");
	printf("     suppressconf      | Stops the chi-squared confidence level tables from being automatically output.\n");
	printf("                       | This is useful for avoiding errors where the number of DoF of a data set is   \n");
	printf("                       | larger than can be handled by the CDF function (see known issues).            \n");
	printf("     extendmodel       | Extends the exported model data beyond the observed values to a user defined  \n");
	printf("                       | frequency (DEFAULT = OFF).                                                    \n");

      }
      else if (tmp_model == 5) {
	printf("\n");
	printf("Main Command:\n");
	printf("=======================================================================================================\n");
	printf("        Command        | Description                                                                   \n");
	printf("-------------------------------------------------------------------------------------------------------\n");
	printf("        fitcioff       | Fit the CI off (KGJP) model of spectral ageing to a given data set.           \n");

	printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
	printf("=======================================================================================================\n");
	printf("       Command         | Description                                                                   \n");
	printf("-------------------------------------------------------------------------------------------------------\n");
	printf("       bfield          | Set the magnetic field strength to use in model fitting (DEFAULT =  1e-9 T).  \n");
	printf("    injectionindex     | Set the injection index to use in model fitting (DEFAULT = 0.6)               \n");
	printf("      modelpower       | Alternative to ‘injectionindex’. Sets the injection index to use in model     \n");
	printf("                       | fitting based on the model power, where α inj = (P − 1)/2 (DEFAULT = 2.2).    \n");
	printf("       myears          | Set the maximum spectral age (on time) in megayears for which the fitting     \n");
	printf("                       | should be attempted (DEFAULT = 10).                                           \n");
	printf("      minmyears        | Set the minimum spectral age (on time) in megayears for which the fitting     \n");
	printf("                       | should be attempted (DEFAULT = 0).                                            \n");
	printf("       ageres          | Change the resolution of the ages attempted when fitting models e.g. 10 gives \n");
	printf("                       | a step size 0.1 for the first level of fitting between 1 and 10 myears        \n");
	printf("                       | (DEFAULT = 10).                                                               \n");
	printf("       minoff          | Set the minimum off time in megayears for which the fitting should be         \n");
	printf("                       | attempted (DEFAULT = 0).                                                      \n");
	printf("       maxoff          | Set the maximum off time in megayears for which the fitting should be         \n");
	printf("                       | attempted (DEFAULT = 20).                                                     \n");
	printf("       levels          | The number of (age) levels deep to go when model fitting. 3 is recommended    \n");
	printf("                       | for testing and initial results, 5 for final results. Note setting this to    \n");
	printf("                       | high values drastically increases computation time! (DEFAULT = 3).            \n");
	printf("        gmin           | Set the minimum value of gamma for model fitting (DEFAULT = 10).              \n");
	printf("        gmax           | Set the maximum value of gamma for model fitting (DEFAULT = 1000000)          \n");
	printf("     printresults      | Turn on and off the printing of individual results when model fitting. The    \n");
	printf("                       | The final summary will still be printed on completion of a model fit (DEFAULT \n");
	printf("                       | = OFF).                                                                       \n");
	printf("     printreject       | Turn on and off the output whether it should be stated that a model should be \n");
	printf("                       | rejected if the average chi-squared is above the 90 per cent confidence inte- \n");
	printf("                       | rval. This only affects the print out to the terminal in the summary, not the \n");
	printf("                       | fitting itself. Note if there are a number of regions where you expect a bad  \n");
	printf("                       | fit e.g. from dynamic range issues, the average may not be a suitable value   \n");
	printf("                       | to use (DEFAULT = OFF).                                                       \n");
	printf("     suppressconf      | Stops the chi-squared confidence level tables from being automatically output.\n");
	printf("                       | This is useful for avoiding errors where the number of DoF of a data set is   \n");
	printf("                       | larger than can be handled by the CDF function (see known issues).            \n");
	printf("     extendmodel       | Extends the exported model data beyond the observed values to a user defined  \n");
	printf("                       | frequency (DEFAULT = OFF).                                                    \n");
      }

      else {
	printf("***WARNING*** Unable to determine which model is help requires! Please contact brats@askanastronomer if this problem persists\n");
	printf("Displaying common commands and parameters for all single injection models instead... \n");
      

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       bfield          | Set the magnetic field strength to use in model fitting (DEFAULT =  1e-9 T).  \n");
      printf("    injectionindex     | Set the injection index to use in model fitting (DEFAULT = 0.6)               \n");
      printf("      modelpower       | Alternative to ‘injectionindex’. Sets the injection index to use in model     \n");
      printf("                       | fitting based on the model power, where α inj = (P − 1)/2 (DEFAULT = 2.2).    \n");
      printf("       myears          | Set the maximum spectral age (on time) in megayears for which the fitting     \n");
      printf("                       | should be attempted (DEFAULT = 10).                                           \n");
      printf("      minmyears        | Set the minimum spectral age (on time) in megayears for which the fitting     \n");
      printf("                       | should be attempted (DEFAULT = 0).                                            \n");
      printf("       ageres          | Change the resolution of the ages attempted when fitting models e.g. 10 gives \n");
      printf("                       | a step size 0.1 for the first level of fitting between 1 and 10 myears        \n");
      printf("                       | (DEFAULT = 10).                                                               \n");
      printf("       minoff          | Set the minimum off time in megayears for which the fitting should be         \n");
      printf("                       | attempted (DEFAULT = 0).                                                      \n");
      printf("       maxoff          | Set the maximum off time in megayears for which the fitting should be         \n");
      printf("                       | attempted (DEFAULT = 20).                                                     \n");
      printf("       levels          | The number of (age) levels deep to go when model fitting. 3 is recommended    \n");
      printf("                       | for testing and initial results, 5 for final results. Note setting this to    \n");
      printf("                       | high values drastically increases computation time! (DEFAULT = 3).            \n");
      printf("        gmin           | Set the minimum value of gamma for model fitting (DEFAULT = 10).              \n");
      printf("        gmax           | Set the maximum value of gamma for model fitting (DEFAULT = 1000000)          \n");
      printf("     printresults      | Turn on and off the printing of individual results when model fitting. The    \n");
      printf("                       | The final summary will still be printed on completion of a model fit (DEFAULT \n");
      printf("                       | = OFF).                                                                       \n");
      printf("     printreject       | Turn on and off the output whether it should be stated that a model should be \n");
      printf("                       | rejected if the average chi-squared is above the 90 per cent confidence inte- \n");
      printf("                       | rval. This only affects the print out to the terminal in the summary, not the \n");
      printf("                       | fitting itself. Note if there are a number of regions where you expect a bad  \n");
      printf("                       | fit e.g. from dynamic range issues, the average may not be a suitable value   \n");
      printf("                       | to use (DEFAULT = OFF).                                                       \n");
      printf("     suppressconf      | Stops the chi-squared confidence level tables from being automatically output.\n");
      printf("                       | This is useful for avoiding errors where the number of DoF of a data set is   \n");
      printf("                       | larger than can be handled by the CDF function (see known issues).            \n");
      }
      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       headerloc       | Enter the location of the file containing the header information.             \n");
      printf("        dataloc        | Enter the location of the file containing the data to be fitted.              \n");

      printf("\n");

      tmp_model = 0;

      break;

    case 42:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      injectionmap     | Plot a map of the best fitting injection index for each region.               \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       contours        | Turn on and off the overlay of flux contours.                                 \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the model fitting results should be plotted.        \n");
      printf("         model         | Select for which model the results should be plotted.                             \n");

      printf("\n");

      break;


    case 43:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("    injectionchimap    | Plot a map of the chi-squared for the best fitting injection index of each    \n");
      printf("                       | region.                                                                       \n");


      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       contours        | Turn on and off the overlay of flux contours.                                 \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the model fitting results should be plotted.        \n");
      printf("         model         | Select for which model the results should be plotted.                         \n");

      printf("\n");

      break;


    case 44:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("    specindexerrors    | (Weighted GSL method only) Map the error values for the spectral index fitting\n");
      printf("                       | as a function of position.                                                    \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("      (plotting)       | Standard global plotting parameters (see the 'plotting' section of help).     \n");
      
 
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("        dataset        | Select for which data set the results should be plotted.                      \n");

      printf("\n");

      break;

      
    case 45:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       plotcimodel     | Plot an example JP model between ‘minmodelfreq’ and ‘maxmodelfreq’ for an     \n");
      printf("                       | arbitrary normalisation.                                                      \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("     modelredshift     | Set the redshift for example spectral ageing models (DEFAULT = 0.2).          \n");
      printf("     minmodelfreq      | Change the minimum frequency in Hz for example spectral ageing models         \n");
      printf("                       | (DEFAULT = 1e7 Hz).                                                           \n");
      printf("     maxmodelfreq      | Change the minimum frequency in Hz for example spectral ageing models         \n");
      printf("                       | (DEFAULT = 1e13 Hz).                                                          \n");
      printf("     modelmyears       | Set the maximum age to plot when outputting models. Model will output at      \n");
      printf("                       | intervals of 1 megayear (DEFAULT = 10 Myr).                                   \n");
      printf("    minmodelmyears     | Set the minimum age to plot when outputting models (DEFAULT = 0 Myr).         \n");
      printf("      modelres         | Set the number of frequency data points to plot for each age (DEFAULT = 100). \n");
      printf("       bfield          | Set the magnetic field strength to use in model fitting (DEFAULT =  1e-9 T).  \n");
      printf("    injectionindex     | Set the injection index to use in model fitting (DEFAULT = 0.6)               \n");
      printf("      modelpower       | Alternative to ‘injectionindex’. Sets the injection index to use in model     \n");
      printf("                       | fitting based on the model power, where α inj = (P − 1)/2 (DEFAULT = 2.2).    \n");
      printf("        gmin           | Set the minimum value of gamma for model fitting (DEFAULT = 10).              \n");
      printf("        gmax           | Set the maximum value of gamma for model fitting (DEFAULT = 1000000)          \n");
      printf("        skip           | Set the value of how many ages to skip when plotting e.g. a value of 10 plots \n");
      printf("                       | models between minmodelmyears and modelmyears at 10 Myr intervals.            \n");


      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("         None          |                                                                               \n");

      printf("\n");

      break;

      
    case 46:

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("       plotcioff       | Plot an example JP model between ‘minmodelfreq’ and ‘maxmodelfreq’ for an     \n");
      printf("                       | arbitrary normalisation.                                                      \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("     modelredshift     | Set the redshift for example spectral ageing models (DEFAULT = 0.2).          \n");
      printf("     minmodelfreq      | Change the minimum frequency in Hz for example spectral ageing models         \n");
      printf("                       | (DEFAULT = 1e7 Hz).                                                           \n");
      printf("     maxmodelfreq      | Change the minimum frequency in Hz for example spectral ageing models         \n");
      printf("                       | (DEFAULT = 1e13 Hz).                                                          \n");
      printf("     modelmyears       | Set the maximum age to plot when outputting models. Model will output at      \n");
      printf("                       | intervals of 1 megayear (DEFAULT = 10 Myr).                                   \n");
      printf("    minmodelmyears     | Set the minimum age to plot when outputting models (DEFAULT = 0 Myr).         \n");
      printf("     minmodeloff       | Set the minimum off age to plot when outputting models (DEFAULT = 0 Myr).     \n");
      printf("     varyoffage        | Set whether the off age should be varied for example models (DEFAULT = YES).  \n");
      printf("      modelres         | Set the number of frequency data points to plot for each age (DEFAULT = 100). \n");
      printf("       bfield          | Set the magnetic field strength to use in model fitting (DEFAULT =  1e-9 T).  \n");
      printf("    injectionindex     | Set the injection index to use in model fitting (DEFAULT = 0.6)               \n");
      printf("      modelpower       | Alternative to ‘injectionindex’. Sets the injection index to use in model     \n");
      printf("                       | fitting based on the model power, where α inj = (P − 1)/2 (DEFAULT = 2.2).    \n");
      printf("        gmin           | Set the minimum value of gamma for model fitting (DEFAULT = 10).              \n");
      printf("        gmax           | Set the maximum value of gamma for model fitting (DEFAULT = 1000000)          \n");
      printf("        skip           | Set the value of how many ages to skip when plotting e.g. a value of 10 plots \n");
      printf("                       | models between minmodelmyears and modelmyears at 10 Myr intervals.            \n");


      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("         None          |                                                                               \n");

      printf("\n");

      break;
      

    case 999:
      fprintf(stderr,"Unknown command - please try again (type 'tasks' for a list of possible commands or 'return' to go back to the main BRATS prompt) \n");
      break;

    case '?':
      fprintf(stderr,"Unknown command - please try again\n");
      break;

    default:
      fprintf(stderr,"Unknown command - please try again\n");
      break;
    
    }
  }

 

  // Help template
      /*

      printf("\n");
      printf("Main Command:\n");
      printf("=======================================================================================================\n");
      printf("        Command        | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("                  |                                                \n");

      printf("\nCommands to set global parameters prior to running the main task (run from the main command prompt): \n");
      printf("=======================================================================================================\n");
      printf("       Command         | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("                  |     \n");

      
      printf("\nParameters set within the main task (requested when the main command is run):\n");
      printf("=======================================================================================================\n");
      printf("       Parameter       | Description                                                                   \n");
      printf("-------------------------------------------------------------------------------------------------------\n");
      printf("          |                          \n");

      // printf("         |              \n");

      printf("\n");

      */


    // */


  return 0;

}
