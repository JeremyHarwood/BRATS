/* 
A function for removing newline characters (\n) from strings, similar to the Perl chomp command.
It is not extensively tested as supplied as is, so you should ensure it meets your needs before using it.

What needs to be passed:

char extstr - A string of the text to be formatted

Basic usage:

#include "chomp.c"

...

chomp(anystring);

Where anystring is a previously defined string. 

If you have any question please email Jeremy.Harwood@physics.org and I will try and answer them as soon as I can.

*/


#include <stdio.h>
#include <string.h>


int chomp(char extstr[]) {

  int len;
  char *eos;

  // Removing newline from string
  //printf("extstr: %s\n", extstr);
  len = strlen(extstr);
  //printf("len: %d\n", len);


  if (len != 0) { // If it's not completely blank

    if(extstr[len-1] == '\n') { 
      extstr[len-1] = 0; }
  
    // Removing trailing whitespace

    eos = extstr + len - 1;

    while(eos > extstr && isspace(*eos)) {
      eos--;
    }

    *(eos+1) = 0;

  }

  return 0;

}
