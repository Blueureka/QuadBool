#ifndef _QIHTMLWriter_h_
#define _QIHTMLWriter_h_

#include <libqi/io/QIWriter.h>

static const char *SPECIAL_CHARS_HTML [] = {

    "&real;",
    "C",
    "&cup;",
    "&infin;",
    "</br>"
 
};

class QIHTMLWriter : public QIWriter {

 public:
  inline QIHTMLWriter() {}
  inline virtual ~QIHTMLWriter(){}

  
  inline string formatBeginDoc     (void)
    {
      
      string result;
      
      result += "<html>\n";
      result += " <head>\n";
      result += "  <link rel=\"stylesheet\" type=\"text/css\" href=\"style_output.css\">\n";
      result += " <body>";

      return result;

    }

  /** A slight modification to convert "\n"'s into HTML <br/> */
  inline string formatDefault      (string rawText)           
    { 
	    
	stringstream ss;
	unsigned short k;

	for (k = 0; k < rawText.length(); k++) {

		if (rawText[k] == '\n') {
			ss << "<br/>";
		}
		else
		 ss << rawText[k];

	}
		
	return ss.str();

    }

  inline string formatSection      (string sectionName)
    { return "<div class=\"section\">" + formatDefault(sectionName) + "</div>\n"; }

  inline string formatSubsection   (string subsectionName)
    { return "<div class=\"subsection\">" + formatDefault(subsectionName) + "</div>\n"; }

  inline string formatMath         (string expr)
    { return "<b class=\"math\">" + formatDefault(expr) + "</b>"; }

  inline string formatInfo         (string information)
    { return "<b class=\"info\">" + formatDefault(information) + "</b>"; }

  inline string formatExponent     (string exponent)
    { return "<sup>" + formatDefault(exponent) + "</sup>"; }
  
  inline string formatSpecial      (SpecialChar id) 
    { return SPECIAL_CHARS_HTML[id]; }

  inline string formatEndDoc     (void)
    {

      string result;
      
      result += " </body>\n";
      result += "</html>";

      return result;

    }
  

};

#endif
