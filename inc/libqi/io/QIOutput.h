/** Data structures storing the result
    of an intersection process, to
    be interpreted differently according to
    the output media (html, console, etc.)
*/
#ifndef _QIOutput_h_
#define _QIOutput_h_

#include <string>

using namespace std;

typedef short VerbosityLevel;

#define	VERBOSITY_PARAM_ONLY	0	/** Only prints parametrizations and Deltas */
#define VERBOSITY_BRUTE		1	/** Brute results */
#define	VERBOSITY_LABELS	2	/** Result with a label */
#define	VERBOSITY_EXHAUSTIVE	3	/** Exhaustive results */


typedef struct _QI_OUTPUT_PARAMETRIZATION_ {

	string label;        /** Ex: "Parametrization of smooth quartic, branch 1" */
	string param;        /** Parametrization [x, y, z, w] */
  /** string x_param;       Parametrization consists in a 4 elements vector:
  string y_param;        [x, y, z, w]                                     
  string z_param;       This parametrization can be in projective or      
  string w_param;       affine.					     */
  	string delta_param;  /** Parametrization of delta.			     */
  	string optimality;   /** "optimal" or "near-optimal"                       */

}QIOutputParametrization;

typedef struct _QI_BENCHMARK_INFO_ {

	string cpu_used_sec;		/** CPU Time used, in seconds */
	
}QIBenchmarkInfo;

typedef struct _QI_OUTPUT_INFO_ {
	
	string interTypeComplexProj; /** Type of intersection in the complex projective space     */
  	string interTypeRealProj;    /** Type of intersection in the real projective space        */
  	string interTypeRealAffine;  /** Type of intersection in the real affine space            */
  	string q1;	             /** Equation of the first quadric (in affine or projective)  */
  	string q2;    	       	     /** Equation of the second quadric (in affine or projective) */
  	string q1_euclideanType;     /** Affine type of the first quadric                         */
  	string q2_euclideanType;     /** Affine type of the second quadric                        */
  	QIOutputParametrization parametrizations[4]; /** List of parametrizations (one for each algebraic
							   component of the intersection) */
	short n_components;		/** Number of components <=> number of parametrizations to
						display */

	/** Additional information needed to customize
	 *  some comments in exhaustive mode (verbose level 3) */
	bool affineParametrizations;

	QIBenchmarkInfo		benchInfo; /** Metrics concerning the calculus itself (CPU used) */

	string surfacePrametrizations[4];/*pencil parametrization*/
}QIOutputInfo;


#endif
