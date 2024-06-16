#ifndef _qi_io_h_
#define _qi_io_h_

/** Include this file if you want to use the
    QI's Input/Output layer, to parse input
    from the user, and pretty display outputs.
    The IO package acts as an interface between
    the QI's calculus kernel and the end user.
*/

#include <libqi/io/QIAnsicolors.h>
#include <libqi/io/QIParser.h>
#include <libqi/io/QIWriter.h>
#include <libqi/io/QIOutputter.h>
#include <libqi/io/QIConsoleWriter.h>
#include <libqi/io/QIHTMLWriter.h>

/** Macro stringification's macro */
#define xstringify(s)	stringify(s)
#define stringify(s)	#s

#endif
