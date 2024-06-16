#include <libqi/io/QIWriter.h>

/** A useful macro to avoid polluting
    the source code with "<<" and "cout"
    everywhere ... */
#define PRINT(str) cout<<str

static string _param 		= "";
static string _paramSet 	= ""; 
static string _paramSetDesc	= "";
static string _monomial		= "";

void QIWriter::formatComments () {

 if (_information->affineParametrizations) {
 	_param		= "u";	
 	_paramSet 	= "the closure of " + formatSpecial(SYM_REAL);
 	_paramSetDesc   = formatSpecial(SYM_REAL) + " " + formatSpecial(SYM_UNION) + " {" + formatSpecial(SYM_INFINITE) + "}";
 	_monomial	= "u" + formatExponent ("k");
 }
 else {
	_param		= "u, v";	
 	_paramSet 	= "real projective space";
 	_paramSetDesc   = "P" + formatExponent("1") + "(" + formatSpecial(SYM_REAL)+")"; 
	_monomial	= "u" + formatExponent ("p") + "*" + "v" + formatExponent ("q");
 }

}

/** To change the layout of the output,
    just change the contents of the following
    method. */
void QIWriter::write () {

  short k;

  /** Adjust some variables to adapt comments in the output */
  formatComments ();
  
  PRINT( formatBeginDoc() );
  
  /** Print benchmark information */
  if ( _verbosityLevel >= VERBOSITY_LABELS) {
	PRINT ( formatDefault("CPU Time: ") + formatMath(_information->benchInfo.cpu_used_sec) +formatDefault (" ms") + formatSpecial(SYM_NEWLINE));
  }
  
  if ( _verbosityLevel >= VERBOSITY_EXHAUSTIVE ) {

   PRINT( formatSection("Input quadrics") );

   if ( _noiseValue.length() > 0 )
	   PRINT ( formatDefault("Noise: ") + formatMath (_noiseValue) + formatSpecial (SYM_NEWLINE) + formatSpecial(SYM_NEWLINE) );

   PRINT( formatDefault("Quadric 1: ") );
   if ( _information->q1_euclideanType.length() > 0)
	   PRINT ( formatInfo(_information->q1_euclideanType) );
   PRINT( formatSpecial(SYM_NEWLINE) );
   PRINT( formatMath (_information->q1) );

   PRINT ( formatSpecial(SYM_NEWLINE) + formatSpecial(SYM_NEWLINE) );
   
   PRINT( formatDefault("Quadric 2: ") );
   if ( _information->q2_euclideanType.length() > 0)
	   PRINT ( formatInfo(_information->q2_euclideanType) );
   PRINT( formatSpecial(SYM_NEWLINE) );
   PRINT( formatMath (_information->q2) );
   PRINT( formatSpecial(SYM_NEWLINE) );
   
   PRINT( formatSection("Type of the intersection") );
  }
  
    
  if ( _verbosityLevel >= VERBOSITY_LABELS ) {
   PRINT( formatDefault("Type in real projective space ") );
   PRINT( formatMath("P"+ formatExponent("3") + "(" + formatSpecial(SYM_REAL) +"): ") );
  }
  
  if ( _verbosityLevel >= VERBOSITY_BRUTE)
   PRINT ( formatInfo(_information->interTypeRealProj) + formatSpecial(SYM_NEWLINE) );
  
  if ( _verbosityLevel >= VERBOSITY_LABELS ) {
   PRINT( formatDefault("Type in complex projective space ") );
   PRINT( formatMath("P"+ formatExponent("3") + "(" + formatSpecial(SYM_COMPLEX) +"): ") );
  }
  
  if ( _verbosityLevel >= VERBOSITY_BRUTE)
   PRINT( formatInfo(_information->interTypeComplexProj) + formatSpecial(SYM_NEWLINE) );

  /** Not yet implemented */
  /**PRINT( formatDefault("Type in real affine space ") );
  PRINT( formatMath(formatSpecial(SYM_REAL)+ formatExponent("3") + ": ") + formatInfo(_information->interTypeRealAffine) + formatSpecial(SYM_NEWLINE));*/

  if (_information->n_components == 0) {
	  if ( _verbosityLevel >= VERBOSITY_LABELS ) {
	   PRINT( formatSpecial(SYM_NEWLINE) );
	   PRINT( formatMath("** There's no component belonging to the real affine space. **") );
	   PRINT( formatSpecial(SYM_NEWLINE) );
	   PRINT( formatSpecial(SYM_NEWLINE) );
	  }
	  PRINT( formatEndDoc() );
	  _cleanup();
	  return;
  }
  
  if ( _verbosityLevel >= VERBOSITY_EXHAUSTIVE ) {
   PRINT( formatSection("Parametrization of the intersection") );

   PRINT( formatDefault("Parametrization of each component of the intersection in ") + formatMath(formatSpecial(SYM_REAL) + formatExponent("3")) +
		  formatDefault(" in homogeneous coordinates ") + formatMath("[x("+_param+"); y("+_param+"); z("+_param+"); w("+_param+")]") + formatDefault(", ")  +
		  formatDefault("where (") + formatMath(_param) + formatDefault(") is the parameter in ") + formatMath(_paramSet + " ( " +_paramSetDesc+ " )") +
		  formatSpecial(SYM_NEWLINE) );

   PRINT( formatDefault("The corresponding affine parametrization is: ") + formatMath("[x("+_param+")/w("+_param+"); y("+_param+")/w("+_param+"); z("+_param+")/w("+_param+")]") +
		  formatSpecial(SYM_NEWLINE) );
  }
  
  k = 0;
  while ( k < _information->n_components ) {
  	
	if ( _verbosityLevel >= VERBOSITY_EXHAUSTIVE ) {
	 PRINT( formatSubsection("["+_information->parametrizations[k].label+"]") );
	}
	else {
	 if (_verbosityLevel >= VERBOSITY_BRUTE)
	  PRINT( formatSpecial(SYM_NEWLINE) );
	 PRINT( formatSpecial(SYM_NEWLINE) );	
	}
	
	if ( _verbosityLevel >= VERBOSITY_BRUTE)
	 PRINT( formatDefault("Parametrization is ") + formatInfo(_information->parametrizations[k].optimality) );
	
	if ( _verbosityLevel >= VERBOSITY_EXHAUSTIVE ) {
	 if ( _information->parametrizations[k].optimality == "OPTIMAL" ) {
		PRINT( formatDefault(": the number of square roots in the coefficients of the ") + formatMath(_monomial)+
				formatDefault(" is minimal.") );
	 }else if ( _information->parametrizations[k].optimality == "NEAR-OPTIMAL" ) {
		PRINT (formatDefault(": there might be one extra square root in the coefficients of the ") + formatMath(_monomial) );
	 }
	}
	
	if ( _verbosityLevel >= VERBOSITY_BRUTE)
	 PRINT ( formatSpecial(SYM_NEWLINE) );
	PRINT ( formatSpecial(SYM_NEWLINE) );
	
	PRINT ( formatMath( _information->parametrizations[k].param) );
	
	k++;
  }

  PRINT( formatEndDoc() );

  /** Do some cleanup */
  _cleanup();
  
}

void QIWriter::_cleanup () {PRINT( formatSpecial(SYM_NEWLINE) );
	if (_information != NULL)
		delete _information;
}


/** The following code is an example:

void QIWriter::write () {

  PRINT( formatBeginDoc() );
  
  PRINT( formatSection("Input quadrics") );
  PRINT( formatDefault("Q1: ") + formatMath("x^2 + y^2 - z^2 - 1") + formatSpecial(SYM_NEWLINE));
  PRINT( formatDefault("Q2: ") + formatMath("x*y - 2*z") );

  PRINT( formatSection("Type of the intersection") );
  PRINT( formatDefault("Type in SYM_COMPLEX projective space ") );
  PRINT( formatMath("P"+ formatExponent("3") + "(" + formatSpecial(SYM_COMPLEX) +"): ") + formatInfo("UNSET") + formatSpecial(SYM_NEWLINE));
  PRINT( formatDefault("Type in SYM_REAL projective space ") );
  PRINT( formatMath("P"+ formatExponent("3") + "(" + formatSpecial(SYM_REAL) +"): ") + formatInfo("UNSET") + formatSpecial(SYM_NEWLINE));
  PRINT( formatDefault("Type in SYM_REAL affine space ") );
  PRINT( formatMath(formatSpecial(SYM_REAL)+ formatExponent("3") + ": ") + formatInfo("UNSET") + formatSpecial(SYM_NEWLINE));

  PRINT( formatSection("Parametrization of the intersection") );
  PRINT( formatDefault("Parametrization of each component of the intersection in ") + formatMath(formatSpecial(SYM_REAL) + formatExponent("3")) +
	   formatDefault(" in homogeneous coordinates ") + formatMath("[x(u), y(u), z(u), w(u)]") + 
	   formatDefault(", where ") + formatMath("u") + formatDefault(" is the parameter in the closure of ") +
	   formatMath(formatSpecial(SYM_REAL) + " (" + formatSpecial(SYM_REAL) + " " + formatSpecial(SYM_UNION) + " {" + formatSpecial(INF) + "})") +
	   formatSpecial(SYM_NEWLINE) );

  PRINT( formatDefault("The corresponding affine parametrization is: ") + formatMath("[x(u)/w(u), y(u)/w(u), z(u)/w(u)]") +
	   formatSpecial(SYM_NEWLINE) );
	   
  PRINT( formatDefault("The parametrizations are ") + formatInfo("UNKNOWN:") + 
	   formatDefault(" the number of square roots in the coefficients of the ") + formatMath("u" + formatExponent("k")) +
	   formatDefault(" is minimal.") );

  PRINT( formatSubsection("Smooth quartic, branch 1") );
  PRINT( formatMath("[2*sqrt(SYM_DELTA), u^3 - 4*u, u*sqrt(SYM_DELTA), u^2 - 4]") + formatSpecial(SYM_NEWLINE));
  PRINT( formatMath(" Delta = u^4 - 5*u^2 + 4") );


  PRINT( formatSubsection("Smooth quartic, branch 2") );
  PRINT( formatMath("[- 2*sqrt(SYM_DELTA), u^3 - 4*u, - u*sqrt(SYM_DELTA), u^2 - 4]") + formatSpecial(SYM_NEWLINE));
  PRINT( formatMath(" Delta = u^4 - 5*u^2 + 4") );

  PRINT( formatEndDoc() );


}*/
