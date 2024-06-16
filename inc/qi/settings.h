/** Defines a list of settings used by the "qi" program. */

#ifndef _SETTINGS_H_
#define _SETTINGS_H_

#include <libqi/rpl/bigint.h>
#include <libqi/qi_io.h>

using namespace rpl;

/** Optimization option (flag) - CPU expensive */
extern char   opt_optimize;

/** Requested size of transformed input quadrics, obtained
    by adding random noise to them. */
extern bigint noise;

/** Format of output's input quadrics : default = false */
extern bool affineQuadrics;

/** Format of output intersection (affine or projective) : default = false */
extern bool affineParametrizations;

/** Omit components which aren't in R^3 : default = false */
extern bool omitComponentsNotInRealSpace;

/** Show the Euclidean type of the input quadrics */
extern bool showEuclideanType;

/** Show cut params : default = false */
extern bool showCutParams;

/** Verbosity level */
extern VerbosityLevel verbosityLevel;

/** Quiet mode: no shell messages nor prompt */
extern bool quiet;

/** Output in HTML format */
extern bool web;

/** Prevents the IO layer from producing multiline output */
extern bool no_multiline;

#endif

