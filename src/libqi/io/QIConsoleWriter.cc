#include <libqi/io/QIConsoleWriter.h>
#include <libqi/io/QIAnsicolors.h>

#define SECTION_LINE    "======================================================================"
#define SUBSECTION_LINE "----------------------------------------------------------------------"

const char *SPECIAL_CHARS_CONSOLE[] = {
    "|R",
    "|C",
    "U",
    "inf",
    "\n"
};

string QIConsoleWriter::formatBeginDoc     (void)                 { /** TESTING ONLY: return "\n--- BEGIN ---\n"; */ return "";}
string QIConsoleWriter::formatDefault      (string rawText)       { return rawText; }
string QIConsoleWriter::formatMath         (string expr)          { return FG_BOLD + expr + FG_NORMAL; }
string QIConsoleWriter::formatInfo         (string information)   { return FG_FACE3 + information + FG_NORMAL; }
string QIConsoleWriter::formatExponent     (string exponent)      { return string("^"+exponent); }
string QIConsoleWriter::formatSpecial      (SpecialChar id)       { return SPECIAL_CHARS_CONSOLE[id]; }
string QIConsoleWriter::formatEndDoc       (void)                 { /** TESTING ONLY: return "\n--- END ---\n"; */ return "";}


/** Example of section:

    ==============================================
                    INPUT QUADRICS               
    ==============================================

*/
string QIConsoleWriter::formatSection      (string sectionName){
  
  string output;
  short n_spaces = (string(SECTION_LINE).length() - sectionName.length()) / 2;

  output += "\n";
  output += SECTION_LINE;
  output += "\n";

  while (n_spaces -- > 0) output += " ";

  output += sectionName;
  output += "\n";

  output += SECTION_LINE;
  output += "\n\n";

  return output;
  
}

/** Example of subsection:
 * 		----------------------------------
                     Smooth quartic, branch 1
		----------------------------------
 *
*/
string QIConsoleWriter::formatSubsection   (string subsectionName) {

  string output;
  short n_spaces = (string(SUBSECTION_LINE).length() - subsectionName.length()) / 2;

  output += "\n";
  output += SUBSECTION_LINE;
  output += "\n";

  while (n_spaces -- > 0) output += " ";

  output += subsectionName;
  output += "\n";

  output += SUBSECTION_LINE;
  output += "\n\n";

  return output;


}

