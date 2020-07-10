/* 
   A function to print to screen the currently loaded datasets and their properties.


   What needs to be passed:

   int indexnumber - Index number of the set
   int imgnum - Number of maps in the set
   char dirname[] - Directory name
   char bgname[] - Background region name
   char regname[] - Region name

   Basic usage:

   #include "list.h"

   ...

 printf("========================================================================\n\n");
    printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

    if (numdatasets > 0) {
      for (i=0;i<numdatasets;i++) {

	list(i, imgnum[i], setname[i], setreg[i], setbg[i]);
      }
    }
    else {
      fprintf(stderr," Error: No datasets have yet been loaded!\n");
    }

    printf("\n========================================================================\n");



   If you have any question please email Jeremy.Harwood@physics.org and I will try and answer them as soon as I can.

*/

#include <stdio.h>
#include <stdlib.h>

int list(const int indexnumber, const int imgnum, const char dirname[], const char regname[], const char bgname[]) {

  printf("   %d       %d         %s           %s       %s \n", indexnumber, imgnum, dirname, bgname, regname);

  return 0;

}

/* 
   A function to print to screen the currently loaded datasets and their properties.


   What needs to be passed:

   int indexnumber - Index number of the set
   int imgnum - Number of maps in the set
   char dirname[] - Directory name
   char bgname[] - Background region name
   char regname[] - Region name
   int xdim - Size of the map in the  X axis
   int ydim - Size of the map in the  Y axis
   float bmaj - Major beam axis size
   float bmin- Minor beam axis size
   float beamarea - Area of the beam

   Basic usage:

   #include "list.h"

   ...

   if (numdatasets > 0) {

	printf("Enter dataset number (-1 to escape): ");
	scanf("%d", &setnumber);
 
	if (setnumber < 0) {
	  printf("Escaping command...\n");
	  break;
	}
	else if (setnumber >= numdatasets) {

	  fprintf(stderr,"\nError: This dataset does not exist! Use list to see those currently loaded\n\n");
	}

	else {

	  i = setnumber; // Keeps the command a bit cleaner
	  listprops(i, imgnum[i], setname[i], setreg[i], setbg[i], xdim[i], ydim[i], bmaj[i], bmin[i], beamarea[i]);
	}
      }
      else {
	fprintf(stderr,"\nError: No datasets have yet been loaded!\n\n");
      }

   If you have any question please email Jeremy.Harwood@physics.org and I will try and answer them as soon as I can.

*/


int listprops(const int indexnumber, const int imgnum, const char dirname[], const char regname[], const char bgname[], const int xdim, const int ydim, const float bmaj, const float bmin, const float beamarea) {

  printf("========================================================================\n\n");

  printf(" Index  | Maps |    Directory    |    Background    |    Region   \n\n");

  printf("   %d       %d         %s           %s       %s \n\n", indexnumber, imgnum, dirname, bgname, regname);

  printf(" Map Dimensions  |     BMaj     |     BMin     |    Beam Area     \n\n");
 
  printf("    %d x %d          %.2f           %.2f            %.2f \n", xdim, ydim, bmaj, bmin, beamarea);

  printf("\n========================================================================\n\n");

  return 0;

}
