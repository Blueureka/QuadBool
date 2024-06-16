#ifndef _QIWriter_h_
#define _QIWriter_h_

#include <libqi/io/QIOutput.h>
#include <string>
#include <iostream>

using namespace std;

/** This class is in charge of formatting
    the results of the intersection in a
    clear manner for the end user.
    It needs to be subclassed by specific
    "writers", able to produce a document
    in a particular format like html or
    plain text.
    See: 
	"QIConsoleWriter.h" (plain text formatter)
	"QIHTMLWriter.h"    (HTML formatter)

    A verbosity level has been defined, to
    specify how much information will be added
    to the result itself, to make the output
    more understandable and friendly:

	   0:	Brute results
	   1:	Results with labels (what they correspond to)
(default)* 2:	Complete page with titles and explications
		
    See "qi_io.h" for the declaration of verbosity levels.

*/


		
class QIWriter {
  
 protected:
  /** Information to be printed out */
  QIOutputInfo *_information;
  
  typedef enum _SPECIAL_CHAR_ {
    SYM_REAL         = 0,
    SYM_COMPLEX      = 1,
    SYM_UNION        = 2,
    SYM_INFINITE     = 3,
    SYM_NEWLINE      = 4
  }SpecialChar;

  VerbosityLevel _verbosityLevel;
  string	 _noiseValue;
  
  /** Frees memory pointed to by "_information" */
  void _cleanup (void);

  
  
 public:

  inline QIWriter() { 
	  _verbosityLevel = VERBOSITY_EXHAUSTIVE;
 	  _noiseValue	  = "";
  }

  inline virtual ~QIWriter() {}

  /** Keep just a reference on a QIOutputInfo structure.
      No hard copy is performed, and the structure is
      automatically freed after the "write" method had
      finished its job. */
  inline void setOutputInformation (QIOutputInfo *info) { _information = info; }

  /** Calls the following format* methods to build
      a structured output of the intersection result. */
  void write (void);
  
  
  /** Sets the verbosity level */
  inline void setVerbosityLevel (VerbosityLevel level) { _verbosityLevel = level; }

  /** Specifies that input quadrics have been perturbated by
   *  a noise */
  inline void setNoiseValue (string value)		{ _noiseValue = string(value); }
  
  /** Adjust some variables to adapt comments in the output */
  void formatComments ();
 
  /** Formatting methods. Depending on the output
      format, the implementation differs. Have a look
      to subclasses for examples of implementations. */
  virtual string formatBeginDoc     (void) = 0;               /** opens a document */
  virtual string formatDefault      (string rawText) = 0;
  virtual string formatSection      (string sectionName) = 0;
  virtual string formatSubsection   (string subsectionName) = 0;
  virtual string formatMath         (string expr) = 0;
  virtual string formatInfo         (string information) = 0;
  virtual string formatExponent     (string exponent) = 0;
  virtual string formatSpecial      (SpecialChar id) = 0;
  virtual string formatEndDoc       (void) = 0;              /** closes a document */

};


#endif

