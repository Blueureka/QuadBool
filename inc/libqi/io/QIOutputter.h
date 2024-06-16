#ifndef _QIOutputter_h_
#define _QIOutputter_h_

#include <libqi/io/QIOutput.h>
#include <libqi/kernel/QIHompoly.h>
#include <libqi/kernel/QIQsicStruct.h>
#include <libqi/kernel/QIElem.h> // For conversion routines
#include <libqi/kernel/QIBench.h>
#include <libqi/io/QIComponentLabels.h>

#include <libqi/rpl/bigint_matrix.h>

using namespace std;
using namespace QI;

/** This class takes the resulting structures
    from the kernel, after calculating the
    intersection.It transforms data into
    a structure of strings (cf. QIOutput.h), 
    used by "QIWriter" which format properly the
    output to the end user.

	$$ JC
	TODO: add an option to show "cut parameters" 

    Note: THIS CLASS DOESN'T SUPPORT TEMPLATES, AS IT
	BRINGS CONFUSION. I'VE CHOOSEN A TEMPORARY
	SOLUTION WHICH IS TO INSTANCIATE DIRECTLY
	WITH "bigint"

*/

/** Threshold to determine how to display the
 *  parametrization: either as a vector, or as
 *  a system.
 *  If one of the four equations X(u,v), Y(u,v), Z(u,v), W(u,v)
 *  has a length greater than the following value, then
 *  the result is displayed as a system. */
#define MAX_PARAMS_LENGTH	20
class QIOutputter {

	private:
		
		VerbosityLevel _verbosityLevel;
		
		/** Resulting output */
		QIOutputInfo *_outputInfo;
		
		/** Options for customizing the output */
		bool	_affineInputQuadrics;
		bool	_affineParametrizations;
		bool	_omitImaginaryParametrizations;
		bool	_multilineEnabled; /** display as a system if needed (see above) */
		bool	_showInputEuclideanType;
		bool    _showCutParams;

		/** Support for special LaTeX notations */
  		bool _useLaTeX;	
		
		/** Utility printing function: recup from the original code inside qi_print.h,
		    with some slight modifications to take the previous options in account. */
		/**inline*/void print_an(const hom_polynomial<bigint> &p_c, const bigint &a, bool &flag, ostream &s);
		/**inline*/void print_hp(const hom_polynomial<bigint> &p1, const hom_polynomial<bigint> &p2,
				const hom_polynomial<bigint> &p3, const hom_polynomial<bigint> &p4,
				const bigint &a, const hom_polynomial<bigint> &D1,
				const hom_polynomial<bigint> &D2, ostream &s);
		/**inline*/void print_cp(const curve_param<bigint> &c1, const curve_param<bigint> &c2,
				const curve_param<bigint> &c3, const curve_param<bigint> &c4,
				const bigint &a, const hom_polynomial<bigint> &D1,
				const hom_polynomial<bigint> &D2, ostream &s);
		
		void print_delta(const bigint &a, const hom_polynomial<bigint> &D1, const hom_polynomial<bigint> &D2, ostream &s);
	
				/** Specific output functions (internal) */
		void _output_component (component <bigint> &c, const bigint_matrix &q1, const bigint_matrix &q2, unsigned short indice);
		
		/** For surfaces */
		string _output_surface (component <bigint> &c);
		
		/** For lines with constraints */
		string _output_lineWithConstraints (component <bigint> &c, const bigint_matrix &q1, const bigint_matrix &q2);
		
		/** For smooth quartics */
		string _output_smoothQuartic (component <bigint> &c, const bigint_matrix &q1, const bigint_matrix &q2);
		
		/** For every other cases */
		string _output_otherCases (component <bigint> &c, const bigint_matrix &q1, const bigint_matrix &q2);
		
	public:
		
		QIOutputter  ();	 /** _output is initialized here, in the constructor */
		inline ~QIOutputter () { /** Don't free _output, it'll be done by QIWriter
					  	which uses it. */ }
		
		
		/** Changes the default options */
		inline void setAffineInputQuadrics    	     (void) 	{ _affineInputQuadrics = true; 		 }
		inline void setAffineParametrizations 	     (void)	{ _affineParametrizations = true; 	 }
		inline void setOmitImaginaryParametrizations (void)	{ _omitImaginaryParametrizations = true; }
		inline void disableMultiline 	 	     (void)	{ _multilineEnabled = false; 		 }		
		inline void showInputEuclideanType 	     (void)	{ _showInputEuclideanType = true;	 }
	
		inline void showCutParams 	     (void)	{ _showCutParams = true;	 }
	
		/** Sets the verbosity level */
		inline void setVerbosityLevel (VerbosityLevel level) { _verbosityLevel = level; }
		
		/** Enables LaTeX notation for sqrt(delta) */
		inline void useLaTeX (bool b)	{ _useLaTeX = b; }
		
		/** Launchs the transformation process. */
		void output (quad_inter<bigint> &quad, const bigint_matrix & q1, const bigint_matrix & q2);
		
		std::vector<std::string> combine_vectors(const quad_inter <bigint> &quad);
		
		/** Gets the output results */
		inline QIOutputInfo *getOutput (void) { return _outputInfo; }
		
		/** $$ JC Added: converts the new type codes (integers)
			in strings */
		void stringify_type(quad_inter<bigint> *quad);
		

};


#endif
