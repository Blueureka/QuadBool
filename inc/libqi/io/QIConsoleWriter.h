#ifndef _QIConsoleWriter_h_
#define _QIConsoleWriter_h_

#include "QIWriter.h"

extern const char *SPECIAL_CHARS_CONSOLE [];
   
class QIConsoleWriter : public QIWriter {

 public:
  inline QIConsoleWriter()  {}
  inline virtual ~QIConsoleWriter() {}
 

  string formatBeginDoc     (void);
  string formatDefault      (string rawText);
  string formatSection      (string sectionName);
  string formatSubsection   (string subsectionName);
  string formatMath         (string expr);
  string formatInfo         (string information);
  string formatExponent     (string exponent);
  string formatSpecial      (SpecialChar id);
  string formatEndDoc       (void);

};

#endif
