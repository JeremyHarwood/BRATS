/* 
   Custom functions used by readline for tab auto complete
*/



char **commandcompletion(const char *commandtext, int start, int end) {
  
  rl_completion_append_character = '\0';
  return rl_completion_matches(commandtext, commandgenerator);
  
}


char *commandgenerator(const char *commandtext, int state) {
  
  static int list_index, len;
  char *command;

  if (!state) {
    list_index = 0;
    len = strlen(commandtext);
  }

  while ((command = commandarray[list_index++])) {
    if (strncmp(command, commandtext, len) == 0) {
      return strdup(command);
    }
  }

  return NULL;

}
