#ifndef _QIParser_h_
#define _QIParser_h_

#include <libqi/rpl/bigint_matrix.h>
#include <libqi/rpl/math_vector.h>
#include <string>
#include <sstream>

#include <libqi/io/QIAnsicolors.h>

/** --- QI parser ---

    In QI, input data are polynomials.
    This data is provided as plain text,
    and may conform to different formats.
    For instance, the 3 following descriptions
    are 3 different representations of the
    same polynomial:

    "-6*x^2+6*x*y+y^2-z"
    "-6*x^2+6*x*y+y^2-z*w"
    "[-6 6 0 0 1 0 0 0 -1 0]"

    Each of these representations must be parsed,
    and transformed into a data structure
    understandable by the calculus kernel.

    QIParser is a parser that takes a string
    as input, recognizes its format, parses it, 
    then creates the associated data structure.

    The internal data structure used by the kernel
    is a matricial representation of the coefficients.
    However, it's handy to use a vectorial representation
    first, as converting it to a matrix is trivial.
*/


using namespace std;
using namespace rpl;

class QIParser {

private:
  
  /** The input string */
  string                 _quadricDesc;
  unsigned short         _string_length;
  
  /** Vectorial representation of the
      polynom specified by the input
      string. */
  math_vector<bigint>    _vectorialDesc;

  /** Matricial representation of the
      polynom specified by the input
      string. */
  bigint_matrix          _matricialDesc;

  /** Vectorial parsing function.
      Example of input: "[-6 6 0 0 1 0 0 0 -1 0]" */
  void _parse_vect (void) noexcept(false);
  
  /** Polynomial parsing function
      Examples of inputs: "-6*x^2+6*x*y+y^2-z"
                          "-6*x^2+6*x*y+y^2-z*w" */ 
  void _parse_poly (void) noexcept(false);

  /** *****************************************  */
  /** Polynomial parser's specific declarations  */
  /** Avoids to declare a list of macros and get */
  /** confused with the source code              */

#include <libqi/io/QIPolynomialParser.h>

  /** ***************************************** */

public:
  

  QIParser ();
  ~QIParser ();

  /** Stores the input string for further
      use by the parser. */
  void setQuadricDesc (string quadricDesc);
  
  /** Returns the original description. */
  string getQuadricDesc (void);

  /** Returns the vectorial description of
      the polynomial, after parsing. */
  math_vector<bigint>   getVectorialDesc (void);

  /** Returns the matricial description of
      the polynomial, after parsing.
      This format is used by the kernel. */
  bigint_matrix         getMatricialDesc (void);


  /** Call this function to launch the parser. */
  void parse (void) noexcept(false);
  

};


/** Non-object functions */

/** Converts a list of coefficients [x^2, xy, xz, xw, y^2, yz, yw, z^2, zw, w^2]
    into its corresponding matrix form */
bigint_matrix vector2quadric (const math_vector <bigint> &qvec);

#endif

