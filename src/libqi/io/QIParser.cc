#include <libqi/rpl/rpl.h>
#include <libqi/io/QIParser.h>
#include <boost/typeof/typeof.hpp>
QIParser::QIParser () {}

QIParser::~QIParser () {}

/** Removes spaces, tabs and newlines inside a string */
void trim (string &s) {
  size_t i;

  while ( 
	  (i = s.find (" ", 0)) != string::npos   || 
	  (i = s.find ("\t", 0)) != string::npos  ||
	  (i = s.find ("\n", 0)) != string::npos
	)
    s.erase (i, 1);
}

void QIParser::setQuadricDesc (string quadricDesc) {

  _quadricDesc = quadricDesc;
  _string_length = _quadricDesc.length();

}

/** Returns the original description. */
string QIParser::getQuadricDesc (void) { return _quadricDesc; }

/** Returns the vectorial description of
    the polynomial, after parsing. */
math_vector<bigint>   QIParser::getVectorialDesc (void) { return _vectorialDesc; }

/** Returns the matricial description of
    the polynomial, after parsing.
    This one is used by the kernel. */
bigint_matrix         QIParser::getMatricialDesc (void) { return _matricialDesc; }

/** Parser selector */
void QIParser::parse (void) noexcept(false) {

  /** Found a vectorial format. */
  if (_quadricDesc.find  ("[", 0) != string::npos &&
      _quadricDesc.find  ("]", 0) != string::npos)
    {
      try {
	/** Do NOT trim the string because
	    spaces are used as coefficient
	    delimiters */
	_parse_vect ();
      } catch (string e) { throw e; }
    }
  else {
    /** Attempts to parse the string as
	a polynomial format. */
    try {
      trim(_quadricDesc);
      _parse_poly ();
    }catch (string e) { throw e; }

  }

  _matricialDesc =  vector2quadric (_vectorialDesc);

}

void QIParser::_parse_poly (void) noexcept(false) {

  CURRENT_SYMBOL = 0;
  _vectorialDesc = rpl::math_vector<bigint>(10, 10);
  
  /** The "IS_A" macro is defined inside "QIPolynomialParser.h"
      Its purpose is to check whether the input string is
      a valid polynomial. */
  try { IS_A(poly); } catch (string e) { throw FG_FACE2 + e + FG_NORMAL; }
    
  _matricialDesc =  vector2quadric (_vectorialDesc);
  
}


void QIParser::_parse_vect (void) noexcept(false) {

  stringstream ss(stringstream::in|stringstream::out);

  ss << _quadricDesc;
  _vectorialDesc.resize(10);
  int nbr; 
  char c;
  ss >>c; // pour enlever ce @#<>~ de charactere "[" au debut!
  for(int i = 0; i < 10 ; ++i)  { 
  ss >> nbr;
  _vectorialDesc[i] = bigint(nbr); 
  }

  if (_vectorialDesc.size() != 10) {
    throw string (FG_FACE2"Bad vector size. Expected 10 elements."FG_NORMAL);
  }

}



bigint_matrix vector2quadric (const math_vector <bigint> &qvec)
{

  bigint_matrix qmat = bigint_matrix(4,4);

  if ((qvec[1].is_even()) && (qvec[2].is_even()) && (qvec[3].is_even()) &&
      (qvec[5].is_even()) && (qvec[6].is_even()) && (qvec[8].is_even()))
    {
      qmat.sto(0,0,qvec[0]);
      qmat.sto(0,1,qvec[1]/2);
      qmat.sto(1,0,qvec[1]/2);
      qmat.sto(0,2,qvec[2]/2);
      qmat.sto(2,0,qvec[2]/2);
      qmat.sto(0,3,qvec[3]/2);
      qmat.sto(3,0,qvec[3]/2);
      qmat.sto(1,1,qvec[4]);  
      qmat.sto(1,2,qvec[5]/2);  
      qmat.sto(2,1,qvec[5]/2);  
      qmat.sto(1,3,qvec[6]/2);    
      qmat.sto(3,1,qvec[6]/2);    
      qmat.sto(2,2,qvec[7]);   
      qmat.sto(2,3,qvec[8]/2);     
      qmat.sto(3,2,qvec[8]/2);     
      qmat.sto(3,3,qvec[9]);     
    }
  else
    {
      qmat.sto(0,0,2*qvec[0]);
      qmat.sto(0,1,qvec[1]);
      qmat.sto(1,0,qvec[1]);
      qmat.sto(0,2,qvec[2]);
      qmat.sto(2,0,qvec[2]);
      qmat.sto(0,3,qvec[3]);
      qmat.sto(3,0,qvec[3]);
      qmat.sto(1,1,2*qvec[4]);  
      qmat.sto(1,2,qvec[5]);  
      qmat.sto(2,1,qvec[5]);  
      qmat.sto(1,3,qvec[6]);    
      qmat.sto(3,1,qvec[6]);    
      qmat.sto(2,2,2*qvec[7]);   
      qmat.sto(2,3,qvec[8]);     
      qmat.sto(3,2,qvec[8]);     
      qmat.sto(3,3,2*qvec[9]);     
    }

  return qmat;

}

/** Polynomial parser specific code. */
#include "QIPolynomialParser.cc"
