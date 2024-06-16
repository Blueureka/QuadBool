#include <qi/settings.h>

//#ifdef HAVE_CONFIG_H
//#include "config.h"
//#endif

/** Optimization option (flag) - CPU expensive */
char opt_optimize = false;

/** Requested size of transformed input quadrics, obtained
  by adding random noise to them */
bigint noise = 0;
  
/** Format of output's input quadrics : default = false */
bool affineQuadrics = false;

/** Format of output intersection (affine or projective) : default = false */
bool affineParametrizations = false;

/** Omit components which aren't in R^3 : default = false */
bool omitComponentsNotInRealSpace = false;

/** Show the Euclidean type of the input quadrics */
bool showEuclideanType = true;

/** Show cut params : default = false */
bool showCutParams = false;

/** Verbosity level: default is 2 (exhaustive) - see QIOutput.h for a list of
    choices */
VerbosityLevel verbosityLevel = VERBOSITY_EXHAUSTIVE;

/** Quiet mode: no shell messages nor prompt */
bool quiet = false;

/** Output in HTML format */
bool web = false;

/** Prevents the IO layer from producing multiline output */
bool no_multiline = false;
