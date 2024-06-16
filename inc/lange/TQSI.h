#pragma once

#include <libqi/qi.h>
#include<Eigen/Eigen>
#include <unsupported/Eigen/Polynomials>

class CSGtree;
#include "CSG/csg.h"
#include<lange/Global.h>


namespace TQSI
{
    struct Factor
    {
        GiNaC::ex expression;
        double exponent;
        Factor(GiNaC::ex e,double ee) : expression(e),exponent(ee){}
    };
    
    struct Root
    {
        double parameter;
        //multiplicity now not need, now just find root and get point, may be alg2 need
        double multiplicity;
        Root(double p,double m): parameter(p),multiplicity(m){}
    };

    struct seq
    {
        std::vector<TQSI::intPoint> points;
        void print()
        {
            for(auto& point:points)
            {
                point.print();
            }   
        }
    };

    // injective map ϕ
    struct injectiveMap
    {
        std::unordered_map<int,int> map;
    };
};


void calculateFraction(double x, int iteration, int& numerator, int& denominator);

// 连分数展开
void continuedFractionExpansion(double x, int maxIterations, int& numerator, int& denominator,double tolerance=1e-8);

GiNaC::ex get_exponent(const GiNaC::ex &expression);

std::vector<TQSI::Factor> split_expression(const GiNaC::ex &expression);

bool solvePolynomialEquation(const TQSI::Factor& factor,const GiNaC::symbol& x,std::vector<TQSI::Root>& roots );

GiNaC::ex polynomialEquationRationalize(const GiNaC::ex &expression,const GiNaC::symbol& x);


std::vector<TQSI::Root> findTQSI(const GiNaC::ex &expression,const GiNaC::symbol& u=GlobalVariable::getU());

inline GiNaC::numeric bigintToNumeric(rpl::bigint a);

bool isPointOnSurface(const TQSI::point& p,const bigint_matrix &Qu);

GiNaC::ex XTAX(const GiNaC::ex& Xc,const bigint_matrix& A);

// GiNaC::ex findIntersectionPointSmoothQuartic(curve_equation& intersectCurve,const bigint_matrix &Qu);
GiNaC::ex findIntersectionPointSmoothQuartic(unsigned int curveIdx,const bigint_matrix &Qu);

TQSI::seq intersectionTriple(curve_equation& intersectCurve,const bigint_matrix &Qu);
TQSI::seq intersectionTriple(const std::pair<unsigned int,unsigned int>& beginEnd,const bigint_matrix &Qu);



TQSI::injectiveMap matchSeq(TQSI::seq& seq1,TQSI::seq& seq2);

void matchVertex(curve_equation& intersectCurve,TQSI::seq& seq);

void matchTwoSeq(TQSI::seq &seq1, TQSI::seq &seq2, double tolerance);
void matchThreeSeq(TQSI::seq& seq1,TQSI::seq& seq2,TQSI::seq& seq3);

void curveAddSeqToPointList(curve_equation& intersectCurve,const TQSI::seq& seq);

void divideCurve(curve_equation& intersectCurve,CSGtree* tree);
void divideCurve(curve_equation& intersectCurve,CSGtree* tree,TQSI::seq& seq);



