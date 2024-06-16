// Elementary quadrics calculations
#include <libqi/rpl/rpl.h>
#include <libqi/kernel/QIElem.h>
#include <libqi/kernel/QINumber.h>

// Enter namespace QI
namespace QI {

// Find inertia of quadric using Descartes' rule
math_vector <int> inertia(const bigint_matrix &q)
{
	math_vector <int> inert(2,2);
	rpl_size_t d = q.get_no_of_rows();

	bigint_matrix id(d,d);
	id.diag(1,0);

	hom_polynomial <bigint> pol;

	if (d == 4)
		pol = det_pencil(q,id);
	else if (d == 3)
		pol = det_pencil3(q,id);
	else
		pol = det_pencil2(q,id);

	rpl_size_t eigen_p = descartes(pol);

	// If the number of positive eigenvalues is the size of the matrix, no need to 
	// go further
	if (eigen_p == d)
		{
			inert[0] = d;
			inert[1] = 0;
		}
	else
		{
			// Compute p(-x)
			for (rpl_size_t i = 0; i <= d; i = i+2)
				pol[i] = -pol[i];

			rpl_size_t eigen_m = descartes(pol);

			if (eigen_m > eigen_p)
				{
					inert[0] = eigen_m;
					inert[1] = eigen_p;
				}
			else
				{
		 			inert[0] = eigen_p;
		 			inert[1] = eigen_m;
				}
		}

	return(inert);
}

// Signed inertia
math_vector <int> signed_inertia(const bigint_matrix &q)
{
	math_vector <int> inert(2,2);
	rpl_size_t d = q.get_no_of_rows();

	bigint_matrix id(d,d);
	id.diag(1,0);

	hom_polynomial <bigint> pol;

	if (d == 4)
		pol = det_pencil(q,id);
	else if (d == 3)
		pol = det_pencil3(q,id);
	else
		pol = det_pencil2(q,id);

	inert[0] = descartes(pol);

	// If the number of positive eigenvalues is the size of the matrix, no need to 
	// go further
	if (inert[0] == d)
		inert[1] = 0;
	else
		{
			// Compute p(-x)
			for (rpl_size_t i = 0; i <= d; i = i+2)
				pol[i] = -pol[i];

			inert[1] = descartes(pol);
		}

	return(inert);
}

// Find inertia of quadric using Descartes' rule when we know its rank is r
math_vector <int> inertia_known_rank(const bigint_matrix &q, 
																			const rpl_size_t &r)
{
	math_vector <int> inert(2,2);
	rpl_size_t d = q.get_no_of_rows();

	bigint_matrix id(d,d);
	id.diag(1,0);

	hom_polynomial <bigint> pol;

	if (d == 4)
		pol = det_pencil(q,id);
	else if (d == 3)
		pol = det_pencil3(q,id);
	else
		pol = det_pencil2(q,id);

	rpl_size_t eigen_p = descartes(pol);

	if (2*eigen_p > r)
		{
			inert[0] = eigen_p;
			inert[1] = r-eigen_p;
		}
	else
		{
			inert[0] = r-eigen_p;
			inert[1] = eigen_p;
		}

	return(inert);
}

// Find the inertia of a non-singular non-rational conic Q whose 3x3 matric is 
// q1+sqrt(D).q2. Being non-singular, the inertia is [3,0] or [2,1]. 
// Output 0 if inertia [3,0] and 1 otherwise 
bool inertia_non_rational_conic(const bigint_matrix &q1, const bigint_matrix &q2, 
																const bigint &D)
{
	bool tmp;

	// The characteric polynomial det(q1+sqrt(D).q2-x.I)= 
	//	P(x) = -x^3 + c2*x^2 +c1*x +c0 
	//	P(-x) = x^3 + c2*x^2 -c1*x +c0 
	// The inertia of Q is [3,0] iff P has 0 positive roots or 0 negative roots 
	// (since all 3 roots are real Descartes gives the exact number of positive roots)
	// tha is by Descarte rule iff
	//	c2<=0, c1<=0, c0<=0 or 
	//	c2>=0, c1<=0, c0>=0. 
	// Let ci = ai +bi*sqrt(D)
	// Recall that sign(a+b*sqrt(D)) =	sign(bigint a, bigint b, bigint D)
	
	// Computed with Maple c0, c1 c2 are as follows.
	bigint a1 = (-q2.member(0,0)*q2.member(1,1)+q2.member(0,2)*q2.member(0,2)
					-q2.member(0,0)*q2.member(2,2)-q2.member(1,1)*q2.member(2,2)
					+q2.member(0,1)*q2.member(0,1)+q2.member(1,2)*q2.member(1,2))*D
					-q1.member(0,0)*q1.member(2,2)-q1.member(0,0)*q1.member(1,1)
					+q1.member(0,2)*q1.member(0,2)-q1.member(1,1)*q1.member(2,2)
					+q1.member(0,1)*q1.member(0,1)+q1.member(1,2)*q1.member(1,2);

	bigint b1 = -q2.member(0,0)*q1.member(2,2)+2*q1.member(0,2)*q2.member(0,2)
					-q1.member(1,1)*q2.member(2,2)-q2.member(0,0)*q1.member(1,1)
					+2*q1.member(1,2)*q2.member(1,2)-q1.member(0,0)*q2.member(1,1)
					-q2.member(1,1)*q1.member(2,2)-q1.member(0,0)*q2.member(2,2)
					+2*q1.member(0,1)*q2.member(0,1);

	if (sign(a1, b1, D)>0)
		tmp = 1; // inertia [2,1]
	else // c1 <= 0
		{
			bigint a2 = q1.member(2,2)+q1.member(1,1)+q1.member(0,0);
			bigint b2 = q2.member(2,2)+q2.member(1,1)+q2.member(0,0);

			bigint a0 = (-q2.member(0,1)*q2.member(0,1)*q1.member(2,2)
			 		-2*q2.member(0,0)*q1.member(1,2)*q2.member(1,2)
			 		-2*q1.member(0,2)*q2.member(0,2)*q2.member(1,1)
			 		-q2.member(0,2)*q2.member(0,2)*q1.member(1,1)
			 		+q2.member(0,0)*q1.member(1,1)*q2.member(2,2)
			 		+2*q1.member(0,1)*q2.member(0,2)*q2.member(1,2)
			 		-q1.member(0,0)*q2.member(1,2)*q2.member(1,2)
			 		+2*q2.member(0,1)*q2.member(0,2)*q1.member(1,2)
			 		+q2.member(0,0)*q2.member(1,1)*q1.member(2,2)
			 		-2*q1.member(0,1)*q2.member(0,1)*q2.member(2,2)
			 		+2*q2.member(0,1)*q1.member(0,2)*q2.member(1,2)
			 		+q1.member(0,0)*q2.member(1,1)*q2.member(2,2))*D
					-q1.member(0,2)*q1.member(0,2)*q1.member(1,1)
					-q1.member(0,1)*q1.member(0,1)*q1.member(2,2)
					-q1.member(0,0)*q1.member(1,2)*q1.member(1,2)
					+q1.member(0,0)*q1.member(1,1)*q1.member(2,2)
					+2*q1.member(0,1)*q1.member(0,2)*q1.member(1,2);

			bigint b0 = (-q2.member(0,1)*q2.member(0,1)*q2.member(2,2)
					+q2.member(0,0)*q2.member(1,1)*q2.member(2,2)
					-q2.member(0,0)*q2.member(1,2)*q2.member(1,2)
					-q2.member(0,2)*q2.member(0,2)*q2.member(1,1)
					+2*q2.member(0,1)*q2.member(0,2)*q2.member(1,2))*D
					+(-q2.member(0,0)*q1.member(1,2)*q1.member(1,2)
					+q2.member(0,0)*q1.member(1,1)*q1.member(2,2)
					+q1.member(0,0)*q1.member(1,1)*q2.member(2,2)
					+2*q1.member(0,1)*q2.member(0,2)*q1.member(1,2)
					+2*q2.member(0,1)*q1.member(0,2)*q1.member(1,2)
					+q1.member(0,0)*q2.member(1,1)*q1.member(2,2)
					+2*q1.member(0,1)*q1.member(0,2)*q2.member(1,2)
					-2*q1.member(0,0)*q1.member(1,2)*q2.member(1,2)
					-q1.member(0,1)*q1.member(0,1)*q2.member(2,2)
					-2*q1.member(0,2)*q2.member(0,2)*q1.member(1,1)
					-2*q1.member(0,1)*q2.member(0,1)*q1.member(2,2)
					-q1.member(0,2)*q1.member(0,2)*q2.member(1,1));

			int sign_c2 = sign(a2, b2, D);
			int sign_c0 = sign(a0, b0, D);
			if (((sign_c2 >= 0) && (sign_c0 >= 0))	|| ((sign_c2 <= 0) && (sign_c0 <= 0)))
				tmp = 0; // inertia [3,0]
			else
				tmp = 1; // inertia [2,1]
		}

	return tmp;
}

// Descartes: count number of sign changes in sequence
rpl_size_t descartes(const polynomial <bigint> &p)
{
	rpl_size_t s_old = p[0].sign();
	rpl_size_t s_changes = 0;

	for (rpl_size_t i = 1; i <= p.degree(); i++)
		if (s_old != 0)
			{
				rpl_size_t s_cur = p[i].sign();
				if (((s_old == 1) && (s_cur == -1)) || ((s_old == -1) && (s_cur == 1)))
		 			{
			 			s_changes++;
			 			s_old = s_cur;
		 			}
			}
		else
			s_old = p[i].sign();
		
	return(s_changes);
}

// Compute the determinantal equation of the pencil generated by (q1,q2), 4x4 matrices
hom_polynomial <bigint> det_pencil(const bigint_matrix &q1, const bigint_matrix &q2)
{
	hom_polynomial <bigint> pol;
	pol.set_degree(4);
	pol[4] = det(q1);
	pol[0] = det(q2);
	
	bigint_matrix qi,qj,qk;
	
	for (rpl_size_t i = 0; i < 4; i++)
		{
			qi = q1;
			column(qi, i) = column(q2, i);
			pol[3] = pol[3]+det(qi);
			for (rpl_size_t j = i+1; j < 4; j++)
				if (i != j)
		 			{
			 			qj = qi;
			 			column(qj, j) = column(q2, j);
			 			pol[2] = pol[2]+det(qj);
			 			for (rpl_size_t k = j+1; k < 4; k++)
				 			if ((k != j) && (k != i))
								{
									qk = qj;
									column(qk, k) = column(q2, k); 
									pol[1] = pol[1]+det(qk);
								}
		 			}
		}
	
	return(pol);
}

// Compute the determinantal equation of the pencil generated by (q1,q2), 3x3 matrices
hom_polynomial <bigint> det_pencil3(const bigint_matrix &q1, const bigint_matrix &q2)
{
	hom_polynomial <bigint> pol;
	pol.set_degree(3);
	pol[3] = det(q1);
	pol[0] = det(q2);

	bigint_matrix qi,qj;
	
	for (rpl_size_t i = 0; i < 3; i++)
		{
			qi = q1;
			column(qi, i) = column(q2, i);
			pol[2] = pol[2]+det(qi);
			for (rpl_size_t j = i+1; j < 3; j++)
				if (i != j)
		 			{
			 			qj = qi;
			 			column(qj, j) = column(q2, j);
			 			pol[1] = pol[1]+det(qj);
		 			}
		}

	return(pol);
}

// Compute the determinantal equation of the pencil generated by (q1,q2), 2x2 matrices
hom_polynomial <bigint> det_pencil2(const bigint_matrix &q1, const bigint_matrix &q2)
{
	hom_polynomial <bigint> pol;
	pol.set_degree(2);
	pol[2] = det(q1);
	pol[0] = det(q2);

	bigint_matrix qi;
	
	for (rpl_size_t i = 0; i < 2; i++)
		{
			qi = q1;
			column(qi, i) = column(q2, i); 
			pol[1] = pol[1]+det(qi);
		}

	return(pol);
}

// Nice priting of plane equation
void print_plane(const bigint_matrix &m, ostream &s)
{
	const char * labels[] = {"x","y","z","w"};
	bool flag = 0;

	for (size_t i = 0; i < m.get_no_of_rows(); i++)
		{
			if (m(i,0) != 0)
				{
					if (flag)
						s << " ";
					if (m(i,0) < 0)
						s << "- ";
					else if (flag)
						s << "+ ";
					if (abs(m(i,0)) != 1)
						s << abs(m(i,0)) << "*";
					s << labels[i];
			
					flag = 1;
				}
		}
}			

// Nice printing of equation of quadric
void print_quadric(const bigint_matrix &m, ostream &s)
{
	const char * labels[] = {"x^2","x*y","x*z","x*w","y^2","y*z","y*w","z^2","z*w","w^2"};
	rpl_size_t i_r[] = {0,0,0,0,1,1,1,2,2,3};
	rpl_size_t i_c[] = {0,1,2,3,1,2,3,2,3,3};
	bool flag = 0, mult;
	bigint co;

	if ((m(0,0).is_even()) && (m(1,1).is_even()) && (m(2,2).is_even()) &&
			(m(3,3).is_even()))
		mult = 1;
	else
		mult = 0;

	for (rpl_size_t i = 0; i < 10; i++)
		{
			if ((i == 0) || (i == 4) || (i == 7) || (i == 9))
				co = m(i_r[i],i_c[i]);
			else
				co = 2*m(i_r[i],i_c[i]);
			if (mult)
				co = co/2;
			if (co != 0)
				{
		 			if (flag)
			 			s << " ";
		 			if (co < 0)
			 			s << "- ";
		 			else if (flag)
			 			s << "+ ";
		 			if (abs(co) != 1)
			 			s << abs(co) << "*";
		 			s << labels[i];

		 			flag = 1;
				}
		}

	if (!flag)
		s << 0;
}

// $$ JC Added
string quad2string(const bigint_matrix quad, bool homogeneous) {
	
	math_vector <bigint> vect = quadtovec(quad);
	
	const char * hom_labels[] = {"x^2","x*y","x*z","x*w","y^2","y*z","y*w","z^2","z*w","w^2"};
	const char * aff_labels[] = {"x^2","x*y","x*z","x","y^2","y*z","y","z^2","z"};
	const char **labels = ( (homogeneous) ? (hom_labels) : (aff_labels) );
	
	short k;
	bool firstCoef = true;
	stringstream result(stringstream::in|stringstream::out);
	bigint coef;
	
	for ( k = 0; k < 10; k ++) {
		
		if (vect[k] == 0) continue;
			
		if ( !firstCoef && vect[k] > 0 ) 
			result << "+";
		else if (vect[k] < 0)			
			result << "-";
		
		coef = abs(vect[k]);
		
		if ( coef != 1 || ( k > 8 && !homogeneous) )
			result << coef;
		
		if (homogeneous || (k <= 8)) {
			if (coef != 1) result << "*";
				result << labels[k];
		}
		
		firstCoef = false;
	}
	
	return (result.str().length() > 0) ? result.str() : "0";
}

/** Returns the Euclidean type of the matrix */
InputEuclideanType getEuclideanType (const bigint_matrix &mat) {

	math_vector<int> S_inertia = inertia (mat);	
	
	bigint_matrix Su(3, 3);
	bigint_matrix B(3, 1);
	bigint_matrix C(1, 3);
	bigint_matrix D(1, 1);

	mat.split_t (Su, B, C, D);
	
	math_vector<int> Su_inertia = inertia (Su);

	/** Seems that C++ misunderstand the meaning of
	 *	 memorization classes ... */
	//unsigned char code;
	int code;
	
	/** According to both inertias (S and Su), we can deduce
	 *	 the Euclidean type of the quadric.
	 *	 Instead of doing a bunch of if...else, let's combine
	 *	 these two inertias to form a unique value and identify
	 *	 the Euclidean type. (hashing of the two inertia) */
	code = (S_inertia[0] << 3) + (S_inertia[1] << 2) + (Su_inertia[0] << 1) + Su_inertia[1];

	switch (code) {
		case 34: return EUCL_INPUT_TYPE_ELLIPSOID;
		case 33: return EUCL_INPUT_TYPE_HYPERBOLOID_OF_TWO_SHEETS;
		case 32: return EUCL_INPUT_TYPE_ELLIPTIC_PARABOLOID;
		case 30: return EUCL_INPUT_TYPE_POINT;
		case 29: return EUCL_INPUT_TYPE_HYPERBOLOID_OF_ONE_SHEET;
		case 27: return EUCL_INPUT_TYPE_HYPERBOLIC_PARABOLOID;
		case 25: return EUCL_INPUT_TYPE_CONE;
		case 24: return EUCL_INPUT_TYPE_ELLIPTIC_CYLINDER;
		case 23: return EUCL_INPUT_TYPE_HYPERBOLIC_CYLINDER;
		case 22: return EUCL_INPUT_TYPE_PARABOLIC_CYLINDER;
		case 20: return EUCL_INPUT_TYPE_LINE;
		case 15: return EUCL_INPUT_TYPE_INTERSECTING_PLANES;
		case 14: return EUCL_INPUT_TYPE_PARALLEL_PLANES;
		case 12: return EUCL_INPUT_TYPE_SIMPLE_PLANE;
		case 10: return EUCL_INPUT_TYPE_DOUBLE_PLANE;
		/** We discriminate the following cases of empty sets to
		 *	give more information about what happens at infinity
		 *	to the user */
		case 38: return	EUCL_INPUT_TYPE_EMPTY;
		case 28: return EUCL_INPUT_TYPE_EMPTY_ONE_POINT_AT_INF;
		case 18: return EUCL_INPUT_TYPE_EMPTY_ONE_LINE_AT_INF;
		case 8:	 return EUCL_INPUT_TYPE_EMPTY_DOUBLE_PLANE_AT_INF;
		case 0:	 return EUCL_INPUT_TYPE_UNIVERSE;

		default: return EUCL_INPUT_TYPE_UNDEFINED;
	}

	return EUCL_INPUT_TYPE_UNDEFINED;
}

// Convert list of coeffs (x^2, xy, xz, xw, y^2, yz, yw, z^2, zw, w^2) to matrix form
bigint_matrix vectoquad(const math_vector <bigint> &qvec)
{
	bigint_matrix qmat(4,4);

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

// Convert matrix to list of coeffs
math_vector <bigint> quadtovec(const bigint_matrix &q)
{
	math_vector <bigint> vec(10,10);

	vec[0] = q(0,0);
	vec[1] = 2*q(0,1);
	vec[2] = 2*q(0,2);
	vec[3] = 2*q(0,3);
	vec[4] = q(1,1);
	vec[5] = 2*q(1,2);
	vec[6] = 2*q(1,3);
	vec[7] = q(2,2);
	vec[8] = 2*q(2,3);
	vec[9] = q(3,3);

	return vec;
}

// Kill row k and column l of matrix q
bigint_matrix mat_minor(const bigint_matrix &q, const rpl_size_t &k, 
												const rpl_size_t &l)
{
	rpl_size_t s = q.get_no_of_rows()-1;
	bigint_matrix q_minor(s,s);

	rpl_size_t inc_i = 0;
	rpl_size_t inc_j;
	
	for (rpl_size_t i = 0; i < s; i++)
		{
			if (i == k)
				inc_i++;
			inc_j = 0;
			for (rpl_size_t j = 0; j < s; j++)
				{
					if (j >= l)
						 inc_j = 1;
					q_minor.sto(i,j,q(i+inc_i,j+inc_j));
				}
		}

	return q_minor;
}

// Cofactor matrix
bigint_matrix cofactors_mat(const bigint_matrix &q)
{
	rpl_size_t d = q.get_no_of_rows();

	bigint_matrix a(d,d);

	for (rpl_size_t i = 0; i < d; i++)
		for (rpl_size_t j = 0; j < d; j++)
			if ( ! ((i+j) % 2) )
				a.sto(i,j,det(mat_minor(q,i,j)));
			else
				a.sto(i,j,-det(mat_minor(q,i,j)));

	return a;
}

// Compute the singular locus of matrix q - Output is a matrix of vectors forming
// a basis of singular locus
bigint_matrix singular(const bigint_matrix &q)
{
	return kernel(q);
}

// Intersect two linear spaces
bigint_matrix linear_intersect(const bigint_matrix &m1, const bigint_matrix &m2)
{
	// Horizontal concatenation of matrices
	rpl_size_t col_size1 = m1.get_no_of_columns();
	rpl_size_t col_sum = col_size1+m2.get_no_of_columns();
	rpl_size_t row_size = m1.get_no_of_rows();
	bigint_matrix m(row_size,col_sum);
	m.compose_h(m1,m2);

	// Compute the relations between the columns of m
	bigint_matrix result = kernel(m);

	// Zeros everywhere: kernel has dimension -1
	if (result.is_column_zero(0))
		result.set_no_of_columns(0);
	else
		{
			// Take only the first set of rows
			result.resize(col_size1,result.get_no_of_columns());

		 	// Multiply by m1
			result = prod(m1, result);
		}

	return result;
}

// Compute a dxd projective transformation sending infinity (0 ... 0 1) to the point
bigint_matrix send_to_infinity(const math_vector <bigint> &point)
{
	rpl_size_t d = point.capacity();

	bigint_matrix proj_mat(d,d);
	
	// Put point as first column of transformation matrix
	column(proj_mat, 0) = point;

	bigint_matrix ker = kernel(trans(proj_mat));

	// Now put point as last column
	proj_mat.swap_columns(0,d-1);

	// Complete matrix with vectors of kernel
	for (size_t i = 0; i < ker.get_no_of_columns(); i++)
		column(proj_mat, i) = column(ker, i);

	return proj_mat;
}

// Compute a 4x4 projective transformation sending the line point1-point2 to z = w = 0
bigint_matrix send_to_zw(const math_vector <bigint> &point1, 
			 const math_vector <bigint> &point2)
{
	bigint_matrix proj_mat(4,4);

	// Put point1 as first column of transformation matrix
	column(proj_mat, 0) = point1;

	// Put point2 as second column of transformation matrix
	column(proj_mat, 1) = point2;

	bigint_matrix ker = kernel(trans(proj_mat));

	// Now put points back in last and next to last position
	proj_mat.swap_columns(0,2);
	proj_mat.swap_columns(1,3);	 

	// Complete matrix with vectors of kernel
	for (size_t i = 0; i < ker.get_no_of_columns(); i++) {
		assert( ker.size1() == 4);
		column(proj_mat, i)	 = column(ker, i);
	}
	
	return proj_mat;
}

// Nice printout of non-rational roots
// projective (a,b+c sqrt{d}) if real
// projective (a,b+c i sqrt{-d}) if complex
void print_root(const math_vector <bigint> &root, ostream &s, const int sig)
{
	if (root[3].is_gt_zero()) // Root is real
		{
			s << "[ " << root[0] << " ";
			if (!root[1].is_zero())
				{
					s << root[1];
					if (sig*root[2] > 0)
						s << "+";
				}
			if (sig*root[2] != 1) 
				{
					if (sig*root[2] == -1)
		 				s << "-";
					else
		 				s << sig*root[2] << "*";
			 	}
			s << "sqrt(" << root[3] << ") ]";
		}
	else // Root is complex
		{
			s << "[ " << root[0] << " ";
			if (!root[1].is_zero())
				{
		 			s << root[1];
		 			if (sig*root[2] > 0)
			 			s << "+";
				}
			if (sig*root[2] != 1) 
				{
					if (sig*root[2] == -1)
		 				s << "-";
					else
		 				s << sig*root[2] << "*";
			 	}
			if (root[3] != -1)
				s << "sqrt(" << -root[3] << ")*";
			s << "i ]";
		}
}

// Compare two (projective) points
bool are_equal(const math_vector <bigint> &v, const math_vector <bigint> &w)
{
	rpl_size_t i;

	// Find first non-zero element
	for (i = 0; i < v.capacity(); i++)
		if (!v[i].is_zero())
			break;

	return (static_cast<math_vector<bigint> >(v*w[i])) == (static_cast<math_vector<bigint> >(w*v[i]));
}

// Tell if k1+sqrt(d)*k2 = (a+sqrt(d)*b)*(k1'+sqrt(d)*k2') if d > 0
// Tell if k1+i*sqrt(d)*k2 = (a+i*sqrt(d)*b)*(k1'+i*sqrt(d)*k2') if d < 0
bool are_equal(const bigint_matrix &k1, const bigint_matrix &k2, 
								const bigint_matrix &k1p, const bigint_matrix &k2p, const bigint &d)
{
	rpl_size_t r = k1.get_no_of_rows();
	rpl_size_t c = k1.get_no_of_columns();

	bigint_matrix comp_tmp1(r,2*c),comp_tmp2(r,2*c),comp(2*r,3*c);
	comp_tmp1.compose_h(k1,d*k2);
	comp_tmp2.compose_h(k2,k1);
	comp.compose_t(comp_tmp1,k1p,comp_tmp2,k2p);
	
	if (rpl::rank(comp) == 3)
		return 0;
	else
		return 1;
}

// Output a quadric of the pencil through the given rational point
bigint_matrix pencil_quadric_through_ratpoint(const bigint_matrix &q1, 
																							const bigint_matrix &q2,
																							const math_vector <bigint> &pt,
																							math_vector <bigint> &l)
{
	rpl_size_t d = q1.get_no_of_rows();
	math_vector <bigint> tmp1(d,d),tmp2(d,d);
	multiply(tmp1,q1,pt);
	multiply(tmp2,q2,pt);

	// The linear equation is l[1]*\lambda - l[0]*\mu = 0
	// Dot products
	multiply(l[1],pt,tmp1);
	multiply(l[0],pt,tmp2);
	l[0].negate();

	// If both l[0] and l[1] are zero: the point is on all quadrics - take q1
	if ((l[0].is_zero()) && (l[1].is_zero()))
		l[0].assign_one();

	optimize(l);

	return l[0]*q1+l[1]*q2;
}

// Sign of the algebraic number a+b*sqrt(c)
int sign(const bigint &a, const bigint &b, const bigint &c)
{
	if (a == 0)
		return b.sign();
	else if ((a > 0) && (b > 0))
		return 1;
	else if ((a < 0) && (b < 0))
		return -1;
	else
		return a.sign()*(a*a-c*b*b).sign();
}

// Replace the ith column by cj times the jth column plus ci times the ith column
void op_col(bigint_matrix &q, const rpl_size_t &i, const rpl_size_t &j, 
						const bigint &ci = 1, const bigint &cj = 1)
{
	math_vector <bigint> tmp = ci* column(q,i)+ cj*column(q,j);
		
	column(q, i) = tmp;
}

// Replace the ith row by cj times the jth row plus ci times the ith row
void op_row(bigint_matrix &q, const rpl_size_t &i, const rpl_size_t &j, 
			 			const bigint &ci = 1, const bigint &cj = 1)
{
	math_vector <bigint> tmp = ci*(math_vector <bigint>)q.get_row_vector(i)+
															cj*(math_vector <bigint>)q.get_row_vector(j);
		
	row(q, i)= tmp;
}

// Gauss reduction of quadratic form
// can is canonical form, tm^T*q*tm = can
// This procedure is not nice and clean... but we don't need something nice and clean
void gauss(const bigint_matrix &q, bigint_matrix &tm, bigint_matrix &can)
{
	rpl_size_t d = q.get_no_of_columns();

	can = q;

	tm.set_no_of_columns(d);
	tm.set_no_of_rows(d);
	tm.diag(1,0);
	rpl_size_t k;
	bigint tmpi = 0;
	bigint tmpj;

	for (rpl_size_t j = 0; j < d; j++)
		{
			k = j+1;
			// Obtain a pivot != 0
			while ((can(j,j) == 0) && (k < d))
				{
		 			op_col(can,j,k);
		 			op_row(can,j,k);
		 			op_row(tm,j,k);
		 			k++;
				}

			if (can(j,j) == 0)
				break;

			// Cancel the rest of the row and column
			for (rpl_size_t i = j+1; i < d; i++)
				{
					if (tmpi == 0)
					  tmpi = can(j,j);
					else
					  tmpi = can(j,j);
			
					tmpj = -can(i,j);
			
					op_col(can,i,j,tmpi,tmpj);
					op_row(can,i,j,tmpi,tmpj);
					op_row(tm,i,j,tmpi,tmpj);
				}
		}

	math_vector <int> in_can = signed_inertia(can);

	// Sort diagonal elements
	for (rpl_size_t i = 0; i < d-1; i++)
		for (rpl_size_t j = i+1; j < d; j++)
			if (((can(i,i) == 0) && (can(j,j) == 0)) ||
		 		((can(i,i) < 0) && (can(j,j) > 0) && (in_can[0] >= in_can[1])) ||
		 		((can(i,i) > 0) && (can(j,j) < -1) && (in_can[1] > in_can[0])))
				{
					can.swap_columns(i,j);
					can.swap_rows(i,j);
					tm.swap_rows(i,j);
				}

	if (can(0,0) < 0)
		can = -can;

	tm = trans(tm);
}

// Compute a random matrix of bigints
bigint_matrix rand_mat(const bigint &size)
{
	bigint_matrix r_mat(4,4);

	bigint sign,tmp;

	for (rpl_size_t i = 0; i < 4; i++)
		for (rpl_size_t j = i; j < 4; j++)
			{
				tmp.randomize(size+1);
				sign.randomize(2);
				if (sign == 0)
					sign = -1;

				r_mat.sto(i,j,sign*tmp);
				r_mat.sto(j,i,r_mat(i,j));
			}

	return r_mat;
}

// Compute a random bigint (including negative numbers)
bigint rand_bigint(const bigint &size)
{
	bigint sign,tmp;
	
	tmp.randomize(size+1);
	sign.randomize(2);

	if (sign == 0)
		sign = -1;

	return sign*tmp;
}	 

// Solve the equation q*v = l*dir in integers (l is a constant)
// q is the matrix of a cone
math_vector <bigint> solve_proj(const bigint_matrix &q, const math_vector <bigint> &dir, 
																const math_vector <bigint> &sing)
{
	bigint_matrix tr = send_to_infinity(sing);

	bigint_matrix q3 = prod( base_matrix<bigint>(prod(trans(tr),q)),tr);
	q3.resize(3,3);

	math_vector <bigint> dir3 = prod(trans(tr), dir);
	dir3.set_size(3);

	optimize(dir3);

	math_vector <bigint> so = prod(cofactors_mat(q3), dir3);
	so.set_size(4);
	so = prod(tr, so);

	optimize(so);

	return so;
}

} // end of namespace QI

