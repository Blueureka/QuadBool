// Structures for storing qsics, components, ...

#ifndef _qi_qsic_struct_h_
#define _qi_qsic_struct_h_

/** QI */
#include <libqi/kernel/QIHompoly.h>
#include <libqi/kernel/QIParamStruct.h>
#include <libqi/kernel/QIComponentTypes.h>

using namespace std;
using namespace rpl;

// Enter namespace QI
namespace QI {

// Three structures here: quad_inter, component and cut_param

template <class T> class quad_inter;
template <class T> class component;
template <class T> class cut_param;

// quad_inter
// The structure used to store an intersection between two quadrics: nothing
// more than a counter and a pointer on a structure of type component

// Possible real pencil type
//  (1,1): empty set
//  (1,2): smooth quartic, two finite components
//  (1,3): smooth quartic, one finite component
//  (1,4): smooth quartic, two infinite components
//  (2,1): point
//  (2,2): nodal quartic with isolated node
//  (2,3): nodal quartic, affinely finite
//  (2,4): nodal quartic, affinely infinite
//  (3,1): empty set
//  (3,2): two points
//  (3,3): two non-secant conics
//  (3,4): two secant conics, affinely finite
//  (3,5): one conic
//  (3,6): two secant conics, affinely infinite
//  (4,1): cuspidal quartic
//  (5,1): double point
//  (5,2): two tangent conics
//  (6,1): empty set
//  (6,2): double conic
//  (7,1): cubic and tangent line
//  (8,1): conic 
//  (8,2): conic and two lines crossing
//  (9,1): double line
//  (9,2): two simple lines and a double line
// (10,1): point
// (10,2): two secant double lines
// (11,1): empty set
// (11,2): quadric of inertia (3,1)
// (11,3): quadric of inertia (2,2)
// (12,1): cubic and secant line
// (12,2): cubic and non-secant line
// (13,1): point
// (13,2): conic and point
// (13,3): conic and two lines not crossing
// (14,1): empty set
// (14,2): skew quadrilateral
// (14,3): two points
// (14,4): two lines
// (15,1): conic and double line
// (16,1): point
// (16,2): two concurrent lines
// (16,3): four concurrent lines
// (17,1): double line
// (17,2): two concurrent lines and a double line
// (18,1): line and triple line
// (19,1): point
// (19,2): two concurrent double lines
// (20,1): quadruple line
// (21,1): quadric of inertia (3,0)
// (21,2): quadric of inertia (2,1)
// (22,1): line and plane
// (23,1): quadruple line
// (23,2): quadruple line
// (24,1): plane
// (25,1): quadric of inertia (2,0)
// (25,2): quadric of inertia (1,1)
// (26,1): double plane
// (27,1): universe

template <class T>
class quad_inter
  {
  public:
    unsigned int nb_cc;

    component <T> *cc;
	
    surface_param <T> s1, s2; // For smooth quartic case

	bigint numberOfSqrt;
    
    unsigned int ctype, rtype; // Type of intersection: ctype (complex), rtype (real)
    /** $$ JC Added: for backcompatibility*/
    /** String versions of each type */
    string sc;
    string sr;
    
    // Constructors
    quad_inter()
      {
				this->create();
      }

    // Create a quad_inter with nb components
    quad_inter(const unsigned int nb)
      {
				this->nb_cc = nb;

				this->cc = new component <T>[nb];
      }

    // Copy constructor
    quad_inter(const quad_inter <T> &a)
      {
				this->create();

				this->assign(a);
      }

    // Destructor
    ~quad_inter()
      {
				this->destroy();
      }

    // Copy assignment
    inline quad_inter <T> & operator = (const quad_inter <T> &a)
      {
				if (this != &a) // Beware of self-assignment
				  {
				    this->destroy();
			
				    this->create();
			
				    this->assign(a);
				  }
			
				return *this;
      }

    // Set type of intersection
    inline void set_type(const unsigned int ct, const unsigned int rt)
      {
				this->ctype = ct;
				this->rtype = rt;
      }

    // Create component cc[id] and conjugate cc[id+1]
    inline void create_two_components(const unsigned int id,
				      const unsigned int tp, const curve_param <T> &c1,
				      const curve_param <T> &c2, const T &D)
      {
				this->cc[id].create_component(tp,c1,c2,D);
				this->cc[id+1].create_component(tp,c1,-c2,D);
      }

    // Create component cc[id] and conjugate cc[id+1]
    inline void create_two_components(const unsigned int id,
				      const unsigned int tp, const curve_param <T> &c1,
				      const curve_param <T> &c2, const curve_param <T> &c3,
				      const curve_param <T> &c4, 
				      const T &D1, const T &D2, const T &D3)
      {
				this->cc[id].create_component(tp,c1,c2,c3,c4,D1,D2,D3);
				this->cc[id+1].create_component(tp,c1,c2,-c3,-c4,D1,D2,D3);
      }

    // Create four components by changing signs
    inline void create_four_components(const unsigned int tp, const curve_param <T> &c1,
				       const curve_param <T> &c2, const curve_param <T> &c3,
				       const T &D1, const T &D2)
      {
				this->cc[0].create_component(tp,c1,c2,c3,D1,D2);
				this->cc[1].create_component(tp,c1,c2,-c3,D1,D2);
				this->cc[2].create_component(tp,c1,-c2,c3,D1,D2);
				this->cc[3].create_component(tp,c1,-c2,-c3,D1,D2);
      }

    // Create four components by changing signs
    inline void create_four_components(const unsigned int tp, const curve_param <T> &c1,
				       const curve_param <T> &c2, const curve_param <T> &c3,
				       const curve_param <T> &c4, const T &D1, 
				       const T &D2)
      {
				this->cc[0].create_component(tp,c1,c2,c3,c4,D1,D2,0);
				this->cc[1].create_component(tp,c1,c2,-c3,-c4,D1,D2,0);
				this->cc[2].create_component(tp,c1,-c2,c3,-c4,D1,D2,0);
				this->cc[3].create_component(tp,c1,-c2,-c3,c4,D1,D2,0);
      }

    // Set the optiflag to the value of the parameter for all components of the intersection
    inline void set_optiflag(const bool &flag)
      {
				for (unsigned int i = 0; i < nb_cc; i++)
	  			this->cc[i].optiflag = flag;
      }

    // Assignment
    inline void assign(const quad_inter <T> &a)
      {
				this->nb_cc = a.nb_cc;
			
				this->ctype = a.ctype;
				this->rtype = a.rtype;
			
				this->s1 = a.s1;
				this->s2 = a.s2;
			
				if (this->nb_cc > 0)
				  {
				    this->cc = new component <T>[this->nb_cc];
			
				    for (unsigned int i = 0; i < this->nb_cc; i++)
				      this->cc[i] = a.cc[i];
				  }
				this->numberOfSqrt=a.numberOfSqrt;
      }

  private:
    inline void create()
      {
				this->nb_cc = 0;
				this->cc = NULL;
				this->s1.assign(surface_param <T>(0));
				this->s2.assign(surface_param <T>(0));
				this->ctype = INTER_TYPE_UNDEFINED;
				this->rtype = INTER_TYPE_UNDEFINED;
				/** $$ JC Added */
				this->sc = "Undefined";
				this->sr = "Undefined";
      }

    inline void destroy()
      {
				if (this->cc != NULL)
	  			delete[] this->cc;
      }
  };

//////////////////////////////////////////////////////////////////////////////////////////////
// Now to component class

// The structure used to store each component of the intersection curve
// In the worst case, the format of the component is:
//
//   c[0] + sqrt(d[0])*c[1] + (c[3] + sqrt(d[0])*c[4])*sqrt(d[1] + sqrt(d[0])*d[2]),
//
// where the cp's are curve_params, d[0] is a bigint, and d[1], d[2] are hom_polys
//
// - type: type of real component (1: smooth quartic, 2: nodal quartic, 
//     3: cuspidal quartic, 4: cubic, 5: conic, 6: line, 7: lines
//     with constraint (four concurrent lines), 8: point, 9: smooth
//     quadric, 10: cone, 11: pair of planes, 12: plane, 13: universe
// - nb_cp: number of cp's in the representation of the parameterization
// - cp's: the cp[] in the representation above
// - d's: the d[] in the representation above 
// - m's: matrices when the output is a surface
// - h: used for bihomogeneous equation when smooth quartic

// For the fucking four concurrent lines with a constraint polynomial, I have
// decided not to create a surface_param but to enter the two endpoints: the first 
// is a point and the second is a curve_param of degree 2

template <class T>
class component
  {
  public:
    unsigned int type;
    bool optiflag; // true if optimal, false if quasi-optimal
    unsigned int mult; // multiplicity of intersection, default is 1

    unsigned int nb_cp;
    unsigned int nb_cut;
    curve_param <T> *c;
    hom_polynomial <T> *d;
    math_matrix <T> m;
    cut_param <T> *p;

    // Constructor
    component()
      {
				this->create();
      }

    // Copy constructor
    component(const component <T> &a)
      {
				this->create();

				this->assign(a);
      }

    // Destructor
    ~component()
      {
				this->destroy();
      }

    // Copy assignment
    component <T> & operator = (const component <T> &a)
      {
				if (this != &a) // Beware of self-assignment
				  {
				    // Clean old stuff
				    this->destroy();
			
				    this->create();
			
				    this->assign(a);
				  }
			
				return *this;
      }

    // Assignment
    inline void assign(const component <T> &a)
      {
				this->nb_cp = a.nb_cp;
				this->nb_cut = a.nb_cut;
				this->type = a.type;
				this->optiflag = a.optiflag;
				this->mult = a.mult;
			
				if ((this->type >= INTER_TYPE_SMOOTH_QUARTIC_BRANCH_1) && (this->type <= INTER_TYPE_POINT)) 
				  // curve component
				  {
				    if (a.c != NULL)
				      {
					this->c = new curve_param <T>[this->nb_cp];
				    
					for (unsigned int i = 0; i < this->nb_cp; i++)
					  this->c[i] = a.c[i];
				    
					if (this->nb_cp > 1)
					  {
					    this->d = new hom_polynomial <T>[this->nb_cp-1];
					
					    for (unsigned int i = 0; i < this->nb_cp-1; i++)
					      this->d[i] = a.d[i];		    
					  }
				      }
				    
				    if (this->nb_cut > 0)
				      {
								this->p = new cut_param <T>[this->nb_cut];
					
								for (unsigned int i = 0; i < this->nb_cut; i++)
					  			this->p[i] = a.p[i];		 
				      }
				  }
				else if ((this->type >= INTER_TYPE_SMOOTH_QUADRIC) && (this->type <= INTER_TYPE_PLANE)) // surface component
				  this->m = a.m;
      }

    // Add a cut parameter [uval0, vval0]
    inline void create_cut_parameter(const int id0, const T &uval0, const T &vval0)
      {
				this->nb_cut = 1;
			
				this->p = new cut_param <T>[1];
			
				this->p[0].create_cut_parameter(id0,uval0,vval0);
      }

    // Add a cut parameter [val0[0], val0[1]]
    inline void create_cut_parameter(const int id0, const math_vector <T> &val0)
      {
				this->create_cut_parameter(id0,val0[0],val0[1]);
      }

    // Add a cut parameter [val[0]+sqrt(D1)*val[1], val[2]] 
    inline void create_cut_parameter(const int id0, 
				     const math_vector <T> &val, const T &D1)
      {
				this->nb_cut = 1;
			
				this->p = new cut_param <T>[1];
			
				if (val[1].is_zero())
				  this->p[0].create_cut_parameter(id0,val[0],val[2]);
				else
				  this->p[0].create_cut_parameter(id0,val[0],val[1],val[2],D1);
      }

    // Add cut parameters [uval0, vval0] and [uval1, vval1]
    inline void create_cut_parameters(const int id0, const T &uval0, const T &vval0, 
				      const int id1, const T &uval1, const T &vval1)
      {
				this->nb_cut = 2;
			
				this->p = new cut_param <T>[2];
			
				this->p[0].create_cut_parameter(id0,uval0,vval0);
				this->p[1].create_cut_parameter(id1,uval1,vval1);
      }

    // Add cut parameters [val0[0], val0[1]] and [val1[0], val1[1]]
    inline void create_cut_parameters(const int id0, const math_vector <T> &val0, 
				      const int id1, const math_vector <T> &val1)
      {
				this->create_cut_parameters(id0,val0[0],val0[1],id1,val1[0],val1[1]);
      }

    // Add cut parameters [val[0], val[1]] and [val[2], val[3]]
    inline void create_cut_parameters(const int id0, const int id1, 
				      const math_vector <T> &val)
      {
				this->create_cut_parameters(id0,val[0],val[1],id1,val[2],val[3]);
      }

    // Add cut parameters [uval0_0+sqrt(D1)*uval1_0, vval_0] and
    //                    [uval0_1+sqrt(D1)*uval1_1, vval_1]
    inline void create_cut_parameters(const int id0,
				      const T &uval0_0, const T &uval1_0, const T &vval_0,
				      const int id1,
				      const T &uval0_1, const T &uval1_1, const T &vval_1, 
				      const T &D1)
      {
				this->nb_cut = 2;
			
				this->p = new cut_param <T>[2];
			
				if (uval1_0.is_zero())
				  this->p[0].create_cut_parameter(id0,uval0_0,vval_0);
				else
				  this->p[0].create_cut_parameter(id0,uval0_0,uval1_0,vval_0,D1);
			
				if (uval1_1.is_zero())
				  this->p[1].create_cut_parameter(id1,uval0_1,vval_1);
				else
				  this->p[1].create_cut_parameter(id1,uval0_1,uval1_1,vval_1,D1);
      }

    // Add cut parameters [val[0]+sqrt(D1)*val[1], val[2]] and
    //                    [val[3]+sqrt(D1)*val[4], val[5]]
    inline void create_cut_parameters(const int id0, const int id1,
				      const math_vector <T> &val, const T &D1)
      {
				this->create_cut_parameters(id0,val[0],val[1],val[2],id1,val[3],val[4],val[5],D1);
      }







    inline void create_cut_parameters(const int id0,
				      const T &uval0_0, const T &uval1_0, const T &uval2_0,
 				      const T &uval3_0, const T &vval_0, 
				      const int id1,
 				      const T &uval0_1, const T &uval1_1, const T &uval2_1, 
 				      const T &uval3_1, const T &vval_1, 
 				      const T &D1, const T &D2) 
       { 
			 	this->nb_cut = 2; 
			 	
			 	this->p = new cut_param <T>[2];
			 	
			 	
			 	
			 	
			 	// Particular cases to handle
			 	
			 	this->p[0].create_cut_parameter(id0,uval0_0,uval1_0,uval2_0,uval3_0,vval_0,D1,D2);
			 	this->p[1].create_cut_parameter(id1,uval0_1,uval1_1,uval2_1,uval3_1,vval_1,D1,D2);


       }

    // Universe...
    inline void create_universe()
      { 
				this->type = INTER_TYPE_UNIVERSE;
				this->nb_cp = 1;
      }

    // Component made of a surface
    inline void create_component(const unsigned int tp, const math_matrix <T> &m1)
      {
				this->type = tp;
				this->nb_cp = 1;

				this->m = m1;
      }

    // Component of the form c1
    inline void create_component(const unsigned int tp, const curve_param <T> &c1)
      {
				this->type = tp;
				this->nb_cp = 1;
			
				this->c = new curve_param <T>[1];
			
				this->c[0] = c1;
      }
	
    // Component of the form c1 + sqrt(D)*c2
    inline void create_component(const unsigned int tp, const curve_param <T> &c1,
				 const curve_param <T> &c2, const T &D)
      {
				this->type = tp;
				this->nb_cp = 2;
			
				this->c = new curve_param <T>[2];
			
				this->c[0] = c1;
				this->c[1] = c2;
			
				this->d = new hom_polynomial <T>[1];
			
				this->d[0] = D;
      }

    // Component of the form c1 + sqrt(D)*c2, D is polynomial
    // or lines with constraints: u*c1+v*c2(s,t) subject to D(s,t) = 0
    inline void create_component(const unsigned int tp, const curve_param <T> &c1,
				 const curve_param <T> &c2, const hom_polynomial <T> &D)
      {
				this->type = tp;
				this->nb_cp = 2;
			
				this->c = new curve_param <T>[2];
			
				this->c[0] = c1;
				this->c[1] = c2;
			
				this->d = new hom_polynomial <T>[1];
			
				this->d[0] = D;
      }

    // Component of the form c1 + sqrt(D1)*c2 + sqrt(D2)*c3
    inline void create_component(const unsigned int tp, const curve_param <T> &c1,
				 const curve_param <T> &c2, const curve_param <T> &c3,
				 const T &D1, const T &D2)
      {
				this->type = tp;
				this->nb_cp = 3;
			
				this->c = new curve_param <T>[3];
			
				this->c[0] = c1;
				this->c[1] = c2;
				this->c[2] = c3;
			
				this->d = new hom_polynomial <T>[2];
			
				this->d[0] = D1;
				this->d[1] = D2;
      }

    // Component of the form c1 + sqrt(D1)*c2 + (c3 + sqrt(D1)*c4)*sqrt(D2 +
    // D3*sqrt(D1)), D2 and D3 are polynomials
    inline void create_component(const unsigned int tp, const curve_param <T> &c1,
				 const curve_param <T> &c2, const curve_param <T> &c3,
				 const curve_param <T> &c4, const T &D1, 
				 const hom_polynomial <T> &D2, const hom_polynomial <T> &D3)
      {
				this->type = tp;
				this->nb_cp = 4;
			
				this->c = new curve_param <T>[4];
			
				this->c[0] = c1;
				this->c[1] = c2;
				this->c[2] = c3;
				this->c[3] = c4;
			
				this->d = new hom_polynomial <T>[3];
			
				this->d[0] = D1;
				this->d[1] = D2;
				this->d[2] = D3;
      }

    // Component of the form c1 + sqrt(D1)*c2 + (c3 + sqrt(D1)*c4)*sqrt(D2 + D3*sqrt(D1))
    inline void create_component(const unsigned int tp, const curve_param <T> &c1,
				 const curve_param <T> &c2, const curve_param <T> &c3,
				 const curve_param <T> &c4, 
				 const T &D1, const T &D2, const T &D3)
      {
				this->create_component(tp,c1,c2,c3,c4,D1,hom_polynomial <T>(D2),
			       hom_polynomial <T>(D3));
      }

      /** $$ JC: checks whether the component is in R^3 */
     bool isInRealAffineSpace (void) {
	
	      unsigned short k;
	
	      for (k = 0; k < nb_cp; k++) {
		      /** $$ JC : c[k][3] is not a pointer so cannot be NULL... */
		      /*if ( c == NULL || c[k][3] == NULL) continue;*/
		      if (c == NULL) continue;
		      if (! c[k][3].is_zero() ) return true;	
	      }
	
	      return false;
     }
      
  private:
    inline void destroy()
      {
				if (this->c != NULL)
				  delete[] this->c;
			
				if (this->d != NULL)
				  delete[] this->d;
			
				if (this->p != NULL)
				  delete[] this->p; 
      }

    inline void create()
      {
				this->optiflag = false;
				this->type = 0;
				this->mult = 1;
			
				this->nb_cp = 0;
				this->nb_cut = 0;
				this->c = NULL;
				this->d = NULL;
				this->m = math_matrix<T>(0,0);
				this->p = NULL;
      }
  };

//////////////////////////////////////////////////////////////////////////////////////////////

// Cut params can be (u0, v0), (u0+sqrt(D1)*u1, v0), (u0+sqrt(D1)*u1+sqrt(D2)*u2, v0) ...
// l is the number of ui's
// id is the index of a component on which the cut param falls

template <class T>
class cut_param
  {
  public:
    unsigned int l;
    int id; // -1 means I don't know

    T *u;
    T v;
    T *D;

    // Constructor
    cut_param()
      {
				this->create();
      }

    // Copy constructor
    cut_param(const cut_param <T> &a)
     	{
				this->create();

				this->assign(a);
      }

    // Destructor
    ~cut_param()
      {
				this->destroy();
      }

    // Copy assignment
    cut_param <T> & operator = (const cut_param <T> &a)
      {
				if (this != &a) // Beware of self-assignment
				  {
				    // Clean old stuff
				    this->destroy();
			
				    this->create();
			
				    this->assign(a);
				  }
			
				return *this;
      }

    inline void assign(const cut_param <T> &a)
      {
				this->l = a.l;
				this->id = a.id;
				this->v = a.v;
			
				if (this->l > 0)
				  {
				    this->u = new T[this->l];
			
				    for (unsigned int i = 0; i < this->l; i++)
				      this->u[i] = a.u[i];
				  }
			
				if (this->l > 1)
				  {
				    unsigned int nbD;
				    
				    if (this->l == 2)
				      nbD = 1;
				    else // l == 4
				      nbD = 2;
				    
				    this->D = new T[nbD];
				    
				    for (unsigned int i = 0; i < nbD; i++)
				      this->D[i].assign(a.D[i]);
				  }
      }

    // A cut parameter (u0, v0), id0 is the index of the other component that is cut
    inline void create_cut_parameter(const int id0, const T &u0, const T &v0)
      {
				this->l = 1;
			
				this->id = id0;
			
			  	this->u = new T[1];
			 	this->u[0] = u0;
				
				this->v = v0;
      }

    // A cut parameter (u0+sqrt(D1)*u1, v0)
    inline void create_cut_parameter(const int id0,
				     const T &u0, const T &u1, const T &v0, const T &D1)
      {
				this->l = 2;
			
				this->id = id0;
			
			 	this->u = new T[2];
			
				this->D = new T[1];
			
			 	this->u[0] = u0;
			 	this->u[1] = u1;
			
				this->D[0] = D1;
				
				this->v = v0;
      }

    // A cut parameter (u0+sqrt(D1)*u1+sqrt(D2)*u2+sqrt(D1*D2)*u3, v0)
    inline void create_cut_parameter(const int id0, 
				     const T &u0, const T &u1, const T &u2, const T &u3,
				     const T &v0, const T &D1, const T &D2)
      {
				this->l = 4;
			
				this->id = id0;
			
				this->u = new T[4];
			
				this->D = new T[2];
			
				this->u[0] = u0;
				this->u[1] = u1;
				this->u[2] = u2;
				this->u[3] = u3;
			
				this->D[0] = D1;
				this->D[1] = D2;
				
				this->v = v0;
      }

  private:
    inline void create()
      {
				this->l = 0;
				this->id = 0;
			      
				this->u = 0;
				this->v = 0;
				this->D = 0;
      }

    inline void destroy()
      {
				if (this->u != NULL)
				  delete[] this->u;
			
				if (this->D != NULL)
				  delete[] this->D;
      }
  };

} // end of namespace QI

#endif
