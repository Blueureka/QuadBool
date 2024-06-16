#include <libqi/io/QIOutputter.h>
#include <libqi/kernel/QICheck.h>
#include <sstream>

QIOutputter::QIOutputter () {
	
	_outputInfo = new QIOutputInfo;
	
	_affineInputQuadrics		 = false;
	_affineParametrizations		 = true;
	_omitImaginaryParametrizations   = false;
	_verbosityLevel = VERBOSITY_EXHAUSTIVE;
	_multilineEnabled		 = true;
	_showInputEuclideanType		 = false;
	_showCutParams  		 = false;

	_outputInfo->q1_euclideanType	 = "";
	_outputInfo->q2_euclideanType	 = "";

	_useLaTeX			 = false;
}

void QIOutputter::print_an(const hom_polynomial <bigint> &p, const bigint &a, bool &flag, ostream &s)
{
	/** If the product is zero, skip printing */
	if ( (p.is_zero()) || (a == 0) ) return;
	
	/** If the square roots equals 1, we can skip it in the
	    product. */
	if (a == 1) 
	{
		if (flag)
		{
			s << " ";
			if (p[xdegree(p)] > 0)
				s << "+ ";
		}
		if (_affineParametrizations)
			p.print_verbose(s,'u','1');
		else
			p.print_verbose (s, 'u', 'v');
	}
	/** Else, we need to print it */
	else
	{
		if (xdegree(p) + ydegree(p) == p.degree()) // a != 1 and p is a monomial
		{
			if (flag)
			{
				s << " ";
				if (p[xdegree(p)] > 0)
					s << "+ ";
			}

			if ((p.degree() == 0) && (p[0] == -1)) // constant polynomial p = -1
				s << "- ";
			else if ((p.degree() != 0) || (p[0] != 1)) {
		// $$ JC
			if (_affineParametrizations)
				p.print_verbose(s,'u','1');
			else
				p.print_verbose (s, 'u', 'v');
		     }		
	      // else p = 1
		}
		else  // a != 1 and p is not a monomial
		{
			if (flag)
				s << " + ";
			s << "(";
	      		if (_affineParametrizations)
				p.print_verbose(s,'u','1');
			else
				p.print_verbose (s, 'u', 'v');
			s << ")";
		}
	  
		if ((p.degree() == 0) && (abs(p[0]) == 1)) // p = +/- 1
		{
			if ( _useLaTeX)
				s << "\\sqrt{" << a << "}";
			else
				s << "sqrt(" << a << ")";
		}else{ 
			if ( _useLaTeX)
				s << "*\\sqrt{" << a << "}";
			else
				s << "*sqrt(" << a << ")";
		}
	}

	flag = 1;
}

void QIOutputter::print_hp(const hom_polynomial <bigint> &p1, const hom_polynomial <bigint> &p2,
					 const hom_polynomial <bigint> &p3, const hom_polynomial <bigint> &p4,
					 const bigint &a, const hom_polynomial <bigint> &D1,
					 const hom_polynomial <bigint> &D2, ostream &s) {
						 
  bool flag = 0;
  bool flag2 = 0;
  bool flag3 = 0;
  
  print_an(p1,(bigint)1,flag,s);
  print_an(p2,a,flag,s); 
  
  //  s << p1 << "  " << p2 << "  " << p3 << "   " << p4 << " " << a << "  " << D1 << "  " << D2 << endl;
  
  if (((!D1.is_zero()) || ((!D2.is_zero()) && (a != 0))) // D1+sqrt(a)*D2 != 0 
      && ((!p3.is_zero()) || ((!p4.is_zero()) && (a != 0)))) // p3 + p4*sqrt(a) != 0 
    {
      if ((xdegree(p3)+ydegree(p3) == p3.degree()) && ((p4.is_zero()) || (a == 0)))
	{
	  if ((p3.degree() == 0) && (abs(p3[0]) == 1))
	    {
	      if (p3[0] == -1)
		{
		  s << " - ";
		  flag = 1;
		}
	      else if (flag)
		s << " + ";
	      
	      flag3 = 1;
	    }
	  else
	    print_an(p3,(bigint)1,flag,s);
	}
      else if ((xdegree(p4)+ydegree(p4) == p4.degree()) && (p3.is_zero()))
	print_an(p4,(bigint)a,flag,s);
      else 
	{
	  if (flag)
	    s << " + ";
	  s << "(";
	  print_an(p3,(bigint)1,flag2,s);
	  print_an(p4,a,flag2,s);
	  s << ")";
	}
      
      if (!flag3)
	s << "*";
      
      if ((D1.degree()>0) ||(D2.degree()>0)) {
	if (_useLaTeX)
		s << "\\sqrt{\\delta}";
	else
      	      	s << "sqrt(Delta)";
      }
      else
	{
	  
	  if ( _useLaTeX)
		  s << "\\sqrt{";
	  else
	  	  s << "sqrt(";
	  
	  flag2 = 0;
	  print_an(D1,(bigint)1,flag2,s);
	  print_an(D2,a,flag2,s);
	  s << ")";
	}
      flag = 1;
    }  
  
  if (!flag)
    s << 0;
  return;						 

}

void QIOutputter::print_cp(const curve_param <bigint> &c1, const curve_param <bigint> &c2,
			 const curve_param <bigint> &c3, const curve_param <bigint> &c4,
			 const bigint &a, const hom_polynomial <bigint> &D1,
			 const hom_polynomial <bigint> &D2, ostream &s)
{
  	// Print the curve param
	/** Prints as a row if the output isn't too large,
	 *  as a system in case it doesn't fit */
	stringstream *params = new stringstream[c1.capacity()];
	bool outofbounds = false;
	
	for (rpl_size_t i = 0; i < c1.capacity(); i++)
	{
		print_hp(c1[i],c2[i],c3[i],c4[i],a,D1,D2,params[i]);
		
		if (params[i].str().length() > MAX_PARAMS_LENGTH)
			outofbounds = true;
	}

	/** Row format */
	if ( ! outofbounds || ! _multilineEnabled ) {

		s << "[";

		 for (rpl_size_t i = 0; i < c1.capacity(); i++){
			if ( i > 0 ) s << ", ";
		 	s << params[i].str();
		 }
			
		s << "]" << endl;
		
	}
	/** System format */
	else {

		string parameters = (_affineParametrizations) ? ("u") : ("u, v");
		string variables[] = {"x", "y", "z", "w"};

		 for (rpl_size_t i = 0; i < c1.capacity(); i++)
			s << variables[i % 4] << "(" << parameters << ") = " << params[i].str() << endl;
	}

	if ((D1.degree() > 0) || (D2.degree() > 0))
		print_delta(a, D1, D2, s);

	delete [] params;
	
	return;
}

void QIOutputter::print_delta(const bigint &a, const hom_polynomial <bigint> &D1, const hom_polynomial <bigint> &D2, ostream &s) {
	 	
	bool flag2 = 0;
	s << endl << "Delta = ";
	print_an(D1,(bigint)1,flag2,s);
	print_an(D2,a,flag2,s);
	s << endl;
}

// Display the type of the intersection
void QIOutputter::stringify_type(quad_inter <bigint> *quad) {
	
	unsigned int ct = quad->ctype, rt = quad->rtype;
	
	if (ct == 1)
	{
		quad->sc = "smooth quartic";

		if (rt == 1)
			quad->sr = "empty set";
		else if (rt == 2)
			quad->sr = "smooth quartic, two finite components";
		else if (rt == 3)
			quad->sr = "smooth quartic, one finite component";
		else if (rt == 4)
			quad->sr = "smooth quartic, two infinite components";
	}
	else if (ct == 2)    
	{
		quad->sc = "nodal quartic";

		if (rt == 1)
			quad->sr = "point";
		else if (rt == 2)
			quad->sr = "nodal quartic with isolated node";
		else if (rt == 3)
			quad->sr = "nodal quartic, affinely finite";
		else if (rt == 4)
			quad->sr = "nodal quartic, affinely infinite";
	}
	else if (ct == 3)    
	{
		quad->sc = "two secant conics";

		if (rt == 1)
			quad->sr = "empty set";
		else if (rt == 2)
			quad->sr = "two points";
		else if (rt == 3)
			quad->sr = "two non-secant conics";
		else if (rt == 4)
			quad->sr = "two secant conics, affinely finite";
		else if (rt == 5)
			quad->sr = "one conic";
		else if (rt == 6)
			quad->sr = "two secant conics, affinely infinite";
	}
	else if (ct == 4)    
	{
		quad->sc = "cuspidal quartic";
      
		if (rt == 1)
			quad->sr = quad->sc;
	}
	else if (ct == 5)    
	{
		quad->sc = "two tangent conics";

		if (rt == 1)
			quad->sr = "double point";
		else if (rt == 2)
			quad->sr = quad->sc;
	}
	else if (ct == 6)    
	{
		quad->sc = "double conic";

		if (rt == 1)
			quad->sr = "empty set";
		else if (rt == 2)
			quad->sr = quad->sc;
	}
	else if (ct == 7)    
	{
		quad->sc = "cubic and tangent line";

		if (rt == 1)
			quad->sr = quad->sc;
	}
	else if (ct == 8)    
	{
		quad->sc = "conic and two lines crossing";

		if (rt == 1)
			quad->sr = "conic";
		else if (rt == 2)
			quad->sr = quad->sc;
	}
	else if (ct == 9)    
	{
		quad->sc = "two simple lines and a double line";

		if (rt == 1)
			quad->sr = "double line";
		else if (rt == 2)
			quad->sr = quad->sc;
	}
	else if (ct == 10)    
	{
		quad->sc = "two secant double lines";

		if (rt == 1)
			quad->sr = "point";
		else if (rt == 2)
			quad->sr = quad->sc;
	}
	else if (ct == 11)    
	{
		quad->sc = "smooth quadric";

		if (rt == 1)
			quad->sr = "empty set";
		else if (rt == 2)
			quad->sr = "quadric of inertia (3,1)";
		else if (rt == 3)
			quad->sr = "quadric of inertia (2,2)";
	}
	else if (ct == 12)    
	{
		quad->sc = "cubic and secant line";

		if (rt == 1)
			quad->sr = quad->sc;
		else if (rt == 2)
			quad->sr = "cubic and non-secant line";
	}
	else if (ct == 13)    
	{
		quad->sc = "conic and two lines not crossing";

		if (rt == 1)
			quad->sr = "point";
		else if (rt == 2)
			quad->sr = "conic and point";
		else if (rt == 3)
			quad->sr = quad->sc;
	}
	else if (ct == 14)    
	{
		quad->sc = "skew quadrilateral";
      
		if (rt == 1)
			quad->sr = "empty set";
		else if (rt == 2)
			quad->sr = quad->sc;
		else if (rt == 3)
			quad->sr = "two points";
		else if (rt == 4)
			quad->sr = "two lines";
	}
	else if (ct == 15)    
	{
		quad->sc = "conic and double line";

		if (rt == 1)
			quad->sr = quad->sc;
	}
	else if (ct == 16)    
	{
		quad->sc = "four concurrent lines";
      
		if (rt == 1)
			quad->sr = "point";
		else if (rt == 2)
			quad->sr = "two concurrent lines";
		else if (rt == 3)
			quad->sr = quad->sc;
	}
	else if (ct == 17)    
	{
		quad->sc = "two concurrent lines and a double line";

		if (rt == 1)
			quad->sr = "double line";
		else if (rt == 2)
			quad->sr = quad->sc;
	}
	else if (ct == 18)    
	{
		quad->sc = "line and triple line";
      
		if (rt == 1)
			quad->sr = quad->sc;
	}
	else if (ct == 19)    
	{
		quad->sc = "two concurrent double lines";

		if (rt == 1)
			quad->sr = "point";
		else if (rt == 2)
			quad->sr = quad->sc;
	}
	else if (ct == 20)    
	{
		quad->sc = "quadruple line";

		if (rt == 1)
			quad->sr = quad->sc;
	}
	else if (ct == 21)    
	{
		quad->sc = "projective cone";
      
		if (rt == 1)
			quad->sr = "quadric of inertia (3,0)";
		else if (rt == 2)
			quad->sr = "quadric of inertia (2,1)";
	}
	else if (ct == 22)    
	{
		quad->sc = "line and plane";

		if (rt == 1)
			quad->sr = quad->sc;
	}
	else if (ct == 23)    
	{
		quad->sc = "quadruple line";

		if (rt == 1)
			quad->sr = quad->sc;
		else if (rt == 2)
			quad->sr = quad->sc;
	}
	else if (ct == 24)    
	{
		quad->sc = "plane";

		if (rt == 1)
			quad->sr = quad->sc;
	}
	else if (ct == 25)    
	{
		quad->sc = "pair of planes";

		if (rt == 1)
			quad->sr = "quadric of inertia (2,0)";
		else if (rt == 2)
			quad->sr = "quadric of inertia (1,1)";
	}
	else if (ct == 26)    
	{
		quad->sc = "double plane";

		if (rt == 1)
			quad->sr = quad->sc;
	}
	else if (ct == 27)    
	{
		quad->sc = "universe";

		if (rt == 1)
			quad->sr = quad->sc;
	}
}

string QIOutputter::_output_surface (component <bigint> &comp) {
	
	stringstream ss(stringstream::in|stringstream::out);
		
	/** Universe */
	if (comp.type == INTER_TYPE_UNIVERSE) {
		ss << "[s, t, u, v]" << endl;
		return ss.str();
		}
	else {
		/** Rational plane */
		if (comp.type == INTER_TYPE_PLANE)
			print_plane ( kernel( trans(comp.m)), ss);
		else
			print_quadric (comp.m, ss);
			
		ss << " = 0" << endl;
		return ss.str();
			
	}
}

string QIOutputter::_output_lineWithConstraints (component <bigint> &comp, const bigint_matrix &q1, const bigint_matrix &q2) {
	
	stringstream ss(stringstream::in|stringstream::out);
	
	curve_param <bigint> c1 = comp.c[0], c2 = comp.c[1];
	hom_polynomial <bigint> p = comp.d[0];

	  // The lines
	surface_param <bigint> lines(4),l_tmp(4);
	hom_polynomial <bigint> u_pol,v_pol;

	  // The first side
	u_pol.assign_x();
	multiply(lines,c1,u_pol);

	  // The other side
	v_pol.assign_y();
	multiply(l_tmp,c2,v_pol);

	  // The sum
	add(lines,l_tmp,lines);

	ss << lines << endl;

	if ( _verbosityLevel >= VERBOSITY_BRUTE) {
	 ss << endl << "constraint polynomial: ";
	}
	 /** $$ JC: I don't know if I should replace 't' with '1'
		when affine output is requested ... */
	 p.print_verbose(ss,'s','t');
	 ss << " = 0" << endl;	

#ifndef NDEBUG
	/**  ss << ">> cut parameter: (u, v) = [";

	  hom_polynomial <bigint> zp;
	  bigint z;

	  // Well, this is just [1 0] but let's pretend we don't know...
	  print_hp((hom_polynomial <bigint>)comp.p[0].u[0],zp,zp,zp,z,zp,zp,ss);
	  ss << ", ";
	  print_hp((hom_polynomial <bigint>)comp.p[0].v,zp,zp,zp,z,zp,zp,ss);
	  ss << "], component cut: " << comp.p[0].id << endl;*/

		ss << endl;
	  check_param(q1,q2,lines,p,ss);
#endif
	
	return ss.str();
}

string QIOutputter::_output_smoothQuartic (component <bigint> &comp, const bigint_matrix &q1, const bigint_matrix &q2) {
	
	stringstream ss(stringstream::in|stringstream::out);
	curve_param <bigint> c1(4), c2(4), c3(4), c4(4);
	bigint D;
	hom_polynomial <bigint> D1,D2;

	// Component is
	//	 c = c1 + sqrt(D). c2 + eps. sqrt(Delta). (c3 + sqrt(D). c4)
  // where Delta = D1 + sqrt(D). D2
	
	c1 = comp.c[0];
	
	if (comp.nb_cp == 2)
		{
			c3 = comp.c[1];
			D = (bigint)1;
			D1 = comp.d[0];
		}
	else // comp.nb_cp == 4
		{
			c2 = comp.c[1];
			c3 = comp.c[2];
		  c4 = comp.c[3];
			D = comp.d[0][0];
			D1 = comp.d[1];
			D2 = comp.d[2];
		}
	
	print_cp(c1,c2,c3,c4,D,D1,D2,ss);

#ifndef NDEBUG
	     /** bigfloat::set_precision(5);

	      // Size of input
	      bigfloat size = size_of_input(q1,q2);

	      ss << ">> size of input: " << size << ", height of Delta: " << 
			      height_of_polynomial(comp.d[0],size) << endl;  */

	      ss << endl;
				check_param(q1,q2,c1,c2,c3,c4,D,D1,D2,ss);
#endif

	return ss.str();
}

string QIOutputter::_output_otherCases (component <bigint> &comp, const bigint_matrix &q1, const bigint_matrix &q2) {
	
	stringstream ss(stringstream::in|stringstream::out);
	curve_param <bigint> c1(4),c2(4),c3(4),c4(4);
	bigint D;
	hom_polynomial <bigint> D1,D2;

	if (comp.nb_cp >= 1)
		c1 = comp.c[0];
	if (comp.nb_cp >= 2)
	{
		c2 = comp.c[1];
		D = comp.d[0][0];
	}
	if (comp.nb_cp >= 3)
	{
		c3 = comp.c[2];
		D1 = comp.d[1];
	}
	if (comp.nb_cp >= 4)
	{
		c4 = comp.c[3];
		D2 = comp.d[2];
	}

	print_cp(c1,c2,c3,c4,D,D1,D2,ss);

/**
	if (_showCutParams)
	  ss << endl << "cutpar: " << comp.nb_cut << endl;
	
	for (unsigned int j = 0; j < comp.nb_cut; j++)
	{
		curve_param <bigint> c(4);

		bigint u0, u1, u2, u3, v;
		bigint D;
		hom_polynomial <bigint> D1, D2;
		curve_param <bigint> cut1(2), cut2(2), cut3(2), cut4(2);

		u0 = comp.p[j].u[0];
		v = comp.p[j].v;

		cut1[0] = u0;
		cut1[1] = v;

		if (comp.p[j].l == 2)
		{
			u1 = comp.p[j].u[1];

			cut2[0] = u1;

			D = comp.p[j].D[0];
		}
		else if (comp.p[j].l == 4)
		{
			u1 = comp.p[j].u[1];
			u2 = comp.p[j].u[2];
			u3 = comp.p[j].u[3];

			cut2[0] = u1;
			cut3[0] = u2;
			cut4[0] = u3;

			D = comp.p[j].D[0];
			D1 = (hom_polynomial <bigint>)comp.p[j].D[1];
		}

#ifndef NDEBUG
		  ss << ">> cut parameter: (u, v) = ";
		  print_cp(cut1,cut2,cut3,cut4,D,D1,D2,ss);
		  ss << ", component cut: " << comp.p[j].id << endl; 
#endif
	}
*/

#ifndef NDEBUG
	    /**  bigfloat::set_precision(5);

	      // Size of input
	      bigfloat size = size_of_input(q1,q2);

	      if (comp.type == 1)
	      {
		      bigfloat m1 = height_of_polynomial(D1,size), m2 = height_of_polynomial(D2,size);

		      if (m1 > m2)
			      ss << ">> size of input: " << size << ", height of Delta: " << m1 << endl;
		      else
			      ss << ">> size of input: " << size << ", height of Delta: " << m2 << endl;
	      }
	      else
		      ss << ">> size of input: " << size << ", height of output: " << 
				      height_of_output(c1,size) << endl;  */

	      ss << endl;
	      check_param(q1,q2,c1,c2,c3,c4,D,D1[0],D2[0],ss);  
#endif

	return ss.str();
}

void QIOutputter::_output_component (component <bigint> &comp, const bigint_matrix &q1, const bigint_matrix &q2, unsigned short indice) {

	QIOutputParametrization *param;
	 
	param = &(_outputInfo->parametrizations[indice]);
	param->label = InterComponentLabels[comp.type];
	param->optimality = ( (comp.optiflag) ? ("OPTIMAL") : ("NEAR-OPTIMAL") );
	
	if (comp.type >= INTER_TYPE_SMOOTH_QUADRIC)		param->param = _output_surface (comp);
	else if (comp.type == INTER_TYPE_LINES_WITH_CONSTRAINT) 	param->param = _output_lineWithConstraints (comp,q1,q2);
	else {
		string multiplicity = "";
		if (comp.mult == 2)
			multiplicity = "double ";
		else if (comp.mult == 3)
			multiplicity = "triple ";
		else if (comp.mult == 4)
			multiplicity = "quadruple ";
		
		param->label = multiplicity + param->label;	
		
		if ( (comp.type == INTER_TYPE_SMOOTH_QUARTIC_BRANCH_1 || comp.type == INTER_TYPE_SMOOTH_QUARTIC_BRANCH_2) ) 
			param->param = _output_smoothQuartic (comp,q1,q2);
		else 
			param->param = _output_otherCases (comp,q1,q2);
	}
}

std::vector<std::string> QIOutputter::combine_vectors(const quad_inter <bigint> &quad)
{
	std::vector<std::string> vectorA_str=quad.s1.surface_paramToString();
	std::vector<std::string> vectorB_str=quad.s2.surface_paramToString();
	bigint D=quad.numberOfSqrt;

	for(int i=0;i<vectorA_str.size();i++)
    {
        vectorA_str[i]+="+sqrt("+D.val.get_str()+")*("+vectorB_str[i]+")";
    }
    return vectorA_str;
}


void QIOutputter::output (quad_inter <bigint> &quad, const bigint_matrix &q1, const bigint_matrix &q2) {


	unsigned short k = 0;
 	stringstream	ss(stringstream::in|stringstream::out);
	
	_outputInfo->q1 = quad2string (q1, !_affineInputQuadrics);
	_outputInfo->q2 = quad2string (q2, !_affineInputQuadrics);
	stringify_type (&quad);
	_outputInfo->interTypeComplexProj = quad.sc;
	_outputInfo->interTypeRealProj    = quad.sr;
	
	_outputInfo->n_components = 0;

	/** We need to forward this information to customize some comments
	 *  produced by "QIWriter" */
	_outputInfo->affineParametrizations = _affineParametrizations;
	
	if ( _showInputEuclideanType ) {
		InputEuclideanType euclidQ1 = getEuclideanType(q1);
		InputEuclideanType euclidQ2 = getEuclideanType(q2);
		
		if (euclidQ1 == 9999) 
			 _outputInfo->q1_euclideanType = "undefined";
		else
  			 _outputInfo->q1_euclideanType = InputEuclideanLabels[euclidQ1];
		if (euclidQ2 == 9999) 
			 _outputInfo->q2_euclideanType = "undefined";
		else
			 _outputInfo->q2_euclideanType = InputEuclideanLabels[euclidQ2];
	}

	// _outputInfo.interTypeRealAffine= quad.sra; */
	
	for (k = 0 ; k < quad.nb_cc; k ++) {
		
		if (_omitImaginaryParametrizations) {
			/**Check if the component is in |R^3 or not */
			if ( ! quad.cc[k].isInRealAffineSpace() ) continue;
		}			
		
		_output_component (quad.cc[k], q1, q2, _outputInfo->n_components);
		_outputInfo->n_components ++;
	}

	std::vector<std::string> surfacePrametrizations=combine_vectors(quad);
	_outputInfo->surfacePrametrizations[0]=surfacePrametrizations[0];
	_outputInfo->surfacePrametrizations[1]=surfacePrametrizations[1];
	_outputInfo->surfacePrametrizations[2]=surfacePrametrizations[2];
	_outputInfo->surfacePrametrizations[3]=surfacePrametrizations[3];

	/** Outputs benchmark information */
	ss << QI_CPU_TIME_MS;
	_outputInfo->benchInfo.cpu_used_sec = ss.str();
}

