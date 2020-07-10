/* 
   Simple routine to run a system command. Simply pass the string to be run by system() such as 'ls -r'.
*/

#include <stdio.h>
#include <stdlib.h>

int syscom(const char command[]) {

  //Checking if we can reach the shell
  if (!system(NULL)) {
    fprintf(stderr,"\n*** Error: Unable to call the processor! Returning to the main promt ***\n\n");
    return 100;
  }

  //Send the command
  system(command);

  return 0;

}
