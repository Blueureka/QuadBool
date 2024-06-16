#ifndef _PARSE_ARGS_H_
#define _PARSE_ARGS_H_

#include <qi/settings.h>
#include <string>

using std::string;

void print_command_help (string progname);

void parseArgs (int argc, char **argv);

/** Check whether the libqi library has been
 *  compiled with the right options.
 *  Uses the flag "LIBQI_TESTING_COMPLIANT"
 *  set at configuration time */
void is_testing_compliant(void);

#endif
