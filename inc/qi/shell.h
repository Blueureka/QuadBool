#ifndef _SHELL_HH_
#define _SHELL_HH_

/** rpl */
#include <libqi/rpl/bigint_matrix.h>
#include <libqi/rpl/math_vector.h>

/** QI */
#include <libqi/qi.h>
#include <qi/settings.h>

typedef struct _SHELL_VAR_ {

  char name[16];
  bigint_matrix value;

}ShellVar;

#ifndef MAX_VARIABLES
#define MAX_VARIABLES 1024
#endif

#define MAX_COMMAND_SIZE	10000


void print_shell_help (void);
short index_of(char *varname);
void set_value_of(char *varname, char *value);
void noisify (bigint_matrix *q1, bigint_matrix *q2);

/** Main */
void shell_main (void);

#endif
