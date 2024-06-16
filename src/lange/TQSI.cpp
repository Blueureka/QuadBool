#include <lange/TQSI.h>

void calculateFraction(double x, int iteration, int &numerator, int &denominator)
{
    if (iteration == 0 || x == 0)
    {
        numerator = 1;
        denominator = round(1.0 / x);
        return;
    }
    double partialDenominator = 1.0 / x;
    // 计算下一个部分的系数
    int a = static_cast<int>(partialDenominator);

    // 计算剩余部分的倒数
    double newFractionPart = partialDenominator - a;

    calculateFraction(newFractionPart, iteration - 1, numerator, denominator);

    // 更新分子和分母
    int newNumerator = denominator;
    int newDenominator = a * denominator + numerator;
    // 将结果更新到传入的变量中
    numerator = newNumerator;
    denominator = newDenominator;
}

void continuedFractionExpansion(double x, int maxIterations, int &numerator, int &denominator, double tolerance)
{
    int a0 = static_cast<int>(x); // 整数部分
    double fractionPart = x - a0; // 剩余部分
    if (fractionPart <= tolerance)
    {
        numerator = a0;
        denominator = 1;
        return;
    }
    numerator = 1, denominator = 1;
    calculateFraction(fractionPart, maxIterations - 1, numerator, denominator);
    // 加上整数部分
    numerator += a0 * denominator;
}

GiNaC::ex get_exponent(const GiNaC::ex &expression)
{
    // 检查表达式是否是幂运算
    if (GiNaC::is_a<GiNaC::power>(expression))
    {
        const GiNaC::power &power_expression = GiNaC::ex_to<GiNaC::power>(expression);
        return power_expression.op(1);
    }
    else
    {
        // 如果不是幂运算，则默认指数为1
        return 1;
    }
}

std::vector<TQSI::Factor> split_expression(const GiNaC::ex &expression)
{
    std::vector<TQSI::Factor> ans;
    std::vector<GiNaC::ex> factors;
    // 检查是否为乘积
    if (GiNaC::is_a<GiNaC::mul>(expression))
    {
        // 遍历所有因子
        for (const auto &factor : expression)
        {
            // factors.emplace_back(expression.op(0),GiNaC::ex_to<GiNaC::numeric>(expression.op(1)).to_double());
            factors.push_back(factor);
            // cout<<factor<<endl;
            //  cout<<factor.op(0)<<" "<<factor.op(1)<<endl;
        }
    }
    else
    {
        // 如果不是乘积，直接将整个表达式视为单个因子
        // factors.emplace_back(expression.op(0),GiNaC::ex_to<GiNaC::numeric>(expression.op(1)).to_double());
        factors.push_back(expression);
    }

    for (const auto &factor : factors)
    {
        // ans.emplace_back(factor.op(0),GiNaC::ex_to<GiNaC::numeric>(factor.op(1)).to_double());
        // cout<<factor<<endl;
        auto ee = get_exponent(factor);
        if (ee == 1)
        {
            ans.emplace_back(factor, 1);
        }
        else
        {
            ans.emplace_back(factor.op(0), GiNaC::ex_to<GiNaC::numeric>(factor.op(1)).to_double());
        }
    }
    return ans;
}

bool solvePolynomialEquation(const TQSI::Factor &factor, const GiNaC::symbol &x, std::vector<TQSI::Root> &roots)
{
    GiNaC::ex expression = factor.expression;

    int degree = expression.degree(x);

    if (degree == 0)
        return false;

    static Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    // Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    Eigen::VectorXd coefficients(degree + 1);
    for (int i = 0; i <= degree; i++)
    {
        auto co = GiNaC::ex_to<GiNaC::numeric>(expression.coeff(x, i)).to_double();
        coefficients(i) = co;
    }
    // cout<<"coefficients= "<<coefficients<<endl;
    solver.compute(coefficients);
    const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType &gen = solver.roots();

    bool haveImageRoot = false;
    for (int i = 0; i < gen.size(); ++i)
    {
        // todo：可能需要加上容差
        if (gen(i).imag() == 0)
        {
            roots.emplace_back(gen(i).real(), factor.exponent);
        }
        else // 虚数根
        {
            // ±∞ 到时候改定义
            //  roots.emplace_back(-1000,factor.exponent);
            //  roots.emplace_back(1000,factor.exponent);
            haveImageRoot = true;
        }
    }
    return haveImageRoot;
}

GiNaC::ex polynomialEquationRationalize(const GiNaC::ex &expression, const GiNaC::symbol &x)
{
    PRINT_VAR(expression);
    GiNaC::ex ans;
    int degree = expression.degree(x);
    int cntIter=15;
    for (int i = 0; i <= degree; i++)
    {
        auto coeff = GiNaC::ex_to<GiNaC::numeric>(expression.coeff(x, i).evalf());
        // cout<<"i= "<<i<<endl;
        // cout<<"coeff= "<<coeff<<endl;
        if (!coeff.is_rational())
        {
            double coeffDouble = coeff.to_double();
            double rund = round(coeffDouble);
            GiNaC::numeric nr(rund);
            int numerator, denominator;
            if ((GiNaC::abs(nr - coeff) - GiNaC::numeric(1e-8)).is_negative())
            {
                // cout<<"rund= "<<rund<<endl;
                continuedFractionExpansion(rund, cntIter, numerator, denominator);
            }
            else
            {
                continuedFractionExpansion(coeffDouble, cntIter, numerator, denominator);
            }

            GiNaC::numeric num(numerator);
            GiNaC::numeric denom(denominator);
            GiNaC::numeric fac = num / denom;
            ans += fac * pow(x, i);
        }
    }
    return ans;
}

std::vector<TQSI::Root> findTQSI(const GiNaC::ex &expression, const GiNaC::symbol &u)
{
    // PRINT_VAR(expression);
    // auto expression_expand=expression.expand();
    // PRINT_VAR(expression_expand);

    auto rationalize = polynomialEquationRationalize(expression.expand(), u);
    cout << rationalize << endl;

    GiNaC::ex BiVarPol = GiNaC::sqrfree(rationalize);
    cout << "BiVarPol= " << BiVarPol << endl;

    std::vector<TQSI::Factor> factors = split_expression(BiVarPol);

    cout << "Factors:" << endl;
    for (const auto &factor : factors)
    {
        cout << factor.expression << " exponent= " << factor.exponent << endl;
    }

    // cout << "solve" << endl;
    std::vector<TQSI::Root> roots;
    bool haveInfRoot = false;
    for (auto &factor : factors)
    {
        if (solvePolynomialEquation(factor, u, roots))
        {
            haveInfRoot = true;
        }
    }
    if (haveInfRoot)
    {
        // roots.emplace_back(-1000,1);
        roots.emplace_back(INF_PARAMETER, 1);
    }

    // for (int i = 0; i < roots.size(); ++i)
    // {
    //     std::cout << roots[i].parameter << " MUL= " << roots[i].multiplicity << std::endl;
    // }
    std::cout << "findTQSI finish" << std::endl;

    return roots;
}

inline GiNaC::numeric bigintToNumeric(rpl::bigint a)
{
    // may improve performance
    return a.val.fits_sint_p() ? GiNaC::numeric(a.val.get_si()) : GiNaC::numeric(a.val.get_str().c_str());
}

TQSI::seq intersectionTriple(const bigint_matrix &q1, const bigint_matrix &q2, const bigint_matrix &q3)
{
    TQSI::seq ans;
    return ans;
}

bool isPointOnSurface(const TQSI::point& p,const bigint_matrix &Qu)
{
    // std::cout<<"isPointOnSurface : begin"<<std::endl;
    if(p.x==INFINITY)return false;

    // std::cout<<"isPointOnSurface : before h"<<std::endl;

    GiNaC::numeric h;
    h += p.x * p.x * bigintToNumeric(Qu(0, 0));     // xx
    h += 2 * p.x * p.y * bigintToNumeric(Qu(0, 1)); // xy
    h += 2 * p.x * p.z * bigintToNumeric(Qu(0, 2)); // xz
    h += 2 * p.x  * bigintToNumeric(Qu(0, 3)); // x
    h += p.y * p.y * bigintToNumeric(Qu(1, 1));     // yy
    h += 2 * p.y * p.z * bigintToNumeric(Qu(1, 2)); // yz
    h += 2 * p.y  * bigintToNumeric(Qu(1, 3)); // y
    h += p.z * p.z * bigintToNumeric(Qu(2, 2));     // zz
    h += 2 * p.z  * bigintToNumeric(Qu(2, 3)); // z
    h +=  bigintToNumeric(Qu(3, 3));     // w

    // std::cout<<"isPointOnSurface : h done"<<std::endl;

    return h*h<=1e-8;
}

GiNaC::ex XTAX(const GiNaC::ex* Xc,const bigint_matrix& A)
{
    GiNaC::ex h;
    h += Xc[0] * Xc[0] * bigintToNumeric(A(0, 0));     // xx
    h += 2 * Xc[0] * Xc[1] * bigintToNumeric(A(0, 1)); // xy
    h += 2 * Xc[0] * Xc[2] * bigintToNumeric(A(0, 2)); // xz
    h += 2 * Xc[0] * Xc[3] * bigintToNumeric(A(0, 3)); // x
    h += Xc[1] * Xc[1] * bigintToNumeric(A(1, 1));     // yy
    h += 2 * Xc[1] * Xc[2] * bigintToNumeric(A(1, 2)); // yz
    h += 2 * Xc[1] * Xc[3] * bigintToNumeric(A(1, 3)); // y
    h += Xc[2] * Xc[2] * bigintToNumeric(A(2, 2));     // zz
    h += 2 * Xc[2] * Xc[3] * bigintToNumeric(A(2, 3)); // z
    h += Xc[3] * Xc[3] * bigintToNumeric(A(3, 3));     // w
    return h;
}

GiNaC::ex findIntersectionPointSmoothQuartic(unsigned int curveIdx,const bigint_matrix &Qu)
{
    auto v=GlobalVariable::getV();
    auto s=GlobalVariable::getS();
    auto t=GlobalVariable::getT();

    std::cout<<"findIntersectionPointSmoothQuartic: begin()"<<std::endl;
    PRINT_VAR(curveIdx);
    PRINT_VAR(GlobalVariable::getAllCurveEquation()->size());
    auto intersectCurve=(*GlobalVariable::getAllCurveEquation())[curveIdx];
    PRINT_VAR(intersectCurve.first_index);
    auto q1=(*GlobalVariable::getQMatrix())[intersectCurve.first_index];
    std::cout<<"findIntersectionPointSmoothQuartic: q1 done"<<std::endl;

    auto f=XTAX(intersectCurve.surfacePrametrizations,q1).expand();
    auto g=XTAX(intersectCurve.surfacePrametrizations,Qu).expand();

    std::cout<<"findIntersectionPointSmoothQuartic: f、g done"<<std::endl;

    std::vector<GiNaC::ex> ai;
    std::vector<GiNaC::ex> bi;
    std::vector<GiNaC::ex> a;
    std::vector<GiNaC::ex> b;
    ai.push_back(f.coeff(t,2));
    ai.push_back(f.coeff(s,1).coeff(t,1));
    ai.push_back(f.coeff(s,2));
    bi.push_back(g.coeff(t,2));
    bi.push_back(g.coeff(s,1).coeff(t,1));
    bi.push_back(g.coeff(s,2));
    
    for(int e=0;e<3;e++)
    {
        a.push_back(ai[e].subs(v==1));
        b.push_back(bi[e].subs(v==1));
        // std::cout<<"ai"<<e<<"= "<<ai[e].expand()<<std::endl;
        // std::cout<<"bi"<<e<<"= "<<bi[e].expand()<<std::endl;
    }

    GiNaC::ex deltaU=a[1]*a[1]-4*a[0]*a[2];
    GiNaC::ex s01=a[0]*b[1]-a[1]*b[0];
    GiNaC::ex s02=a[0]*b[2]-a[2]*b[0];
    GiNaC::ex s12=a[1]*b[2]-a[2]*b[1];
    GiNaC::ex resU=s02*s02-s01*s12;

    return resU;
}

TQSI::seq intersectionTriple(curve_equation& intersectCurve, const bigint_matrix &Qu)
{
    TQSI::seq ans;
    // std::cout<<Qu<<std::endl;
    GiNaC::ex h;

    if(intersectCurve.have_delta)
    {
        std::cout<<"intersectionTriple findIntersectionPointSmoothQuartic"<<std::endl;
        // h= findIntersectionPointSmoothQuartic(intersectCurve,Qu);
        std::cout<<"intersectionTriple findIntersectionPointSmoothQuartic done"<<std::endl;
    }
    else
    {
        std::vector<GiNaC::ex> Xc;
        Xc.push_back(intersectCurve.x_equation);
        Xc.push_back(intersectCurve.y_equation);
        Xc.push_back(intersectCurve.z_equation);
        Xc.push_back(intersectCurve.w_equation);
        
        h += Xc[0] * Xc[0] * bigintToNumeric(Qu(0, 0));     // xx
        h += 2 * Xc[0] * Xc[1] * bigintToNumeric(Qu(0, 1)); // xy
        h += 2 * Xc[0] * Xc[2] * bigintToNumeric(Qu(0, 2)); // xz
        h += 2 * Xc[0] * Xc[3] * bigintToNumeric(Qu(0, 3)); // x
        h += Xc[1] * Xc[1] * bigintToNumeric(Qu(1, 1));     // yy
        h += 2 * Xc[1] * Xc[2] * bigintToNumeric(Qu(1, 2)); // yz
        h += 2 * Xc[1] * Xc[3] * bigintToNumeric(Qu(1, 3)); // y
        h += Xc[2] * Xc[2] * bigintToNumeric(Qu(2, 2));     // zz
        h += 2 * Xc[2] * Xc[3] * bigintToNumeric(Qu(2, 3)); // z
        h += Xc[3] * Xc[3] * bigintToNumeric(Qu(3, 3));     // w

        std::cout << "h(ξ)= " << h << std::endl;
        if(intersectCurve.have_delta)
        {
            h=h.subs(intersectCurve.Delta==intersectCurve.delta_equation);
            std::cout << "h(ξ) with delta = " << h << std::endl;
        }
    }

    
    // if {h(ξ) ≡ 0} report C and return.
    if (h.is_zero())
    {
        return ans;
    }
    std::vector<TQSI::Root> roots;
    roots = std::move(findTQSI(h, intersectCurve.u));

    std::cout << "print roots" << std::endl;
    for (int i = 0; i < roots.size(); ++i)
    {
        std::cout << roots[i].parameter << " multiplicity = " << roots[i].multiplicity << std::endl;
        auto parameter = roots[i].parameter;
        if (parameter != INF_PARAMETER)
        {
            auto point=intersectCurve.getPoint(parameter, 1, 1);
            if(!isPointOnSurface(point,Qu))
            {

            }
            ans.points.emplace_back(point, parameter);
        }
        else
        {
            if (intersectCurve.isLoop)
            {
                if(isPointOnSurface(intersectCurve.positiveInf,Qu))
                ans.points.emplace_back((intersectCurve.positiveInf + intersectCurve.negativeInf) / 2, -parameter);
            }
            else
            {
                if(isPointOnSurface(intersectCurve.negativeInf,Qu))
                ans.points.emplace_back(intersectCurve.negativeInf, parameter);
                if(isPointOnSurface(intersectCurve.positiveInf,Qu))
                ans.points.emplace_back(intersectCurve.positiveInf, parameter);
            }
        }
    }
    std::sort(ans.points.begin(), ans.points.end(), [&](TQSI::intPoint a, TQSI::intPoint b)
              { return a.parameter < b.parameter; });

    // for(auto s:ans.points)
    // {
    //     cout<<"parameter= "<<s.parameter<<endl;
    // }

    return ans;
}

TQSI::seq intersectionTriple(const std::pair<unsigned int,unsigned int>& beginEnd,const bigint_matrix &Qu)
{   
    //smooth quartic
    TQSI::seq ans;
    // std::cout<<Qu<<std::endl;

    //!TODO: AABB 

    auto no_op_deleter = [](curve_equation*) {};
    GiNaC::ex h;
    // std::shared_ptr<curve_equation> branch1SPtr(&(*GlobalVariable::getAllCurveEquation())[beginEnd.first],no_op_deleter);
    // std::shared_ptr<curve_equation> branch2SPtr;
    curve_equation* branch1SPtr=&(*GlobalVariable::getAllCurveEquation())[beginEnd.first];
    curve_equation* branch2SPtr;

    cout<<"all curve size = "<<GlobalVariable::getAllCurveEquation()->size()<<endl;
    PRINT_VAR(beginEnd.first);
    PRINT_VAR(beginEnd.second);

    std::cout<<"intersectionTriple init done"<<std::endl;

    if(beginEnd.first+1!=beginEnd.second)
    {
        std::cout<<"intersectionTriple findIntersectionPointSmoothQuartic"<<std::endl;
        // branch2SPtr.reset(&(*GlobalVariable::getAllCurveEquation())[beginEnd.second],no_op_deleter);
        branch2SPtr=&(*GlobalVariable::getAllCurveEquation())[beginEnd.second];
        // PRINT_VAR(beginEnd.first);
        // h= findIntersectionPointSmoothQuartic((*GlobalVariable::getAllCurveEquation())[beginEnd.first],Qu);
        h= findIntersectionPointSmoothQuartic(beginEnd.first,Qu);
        std::cout<<"intersectionTriple findIntersectionPointSmoothQuartic done"<<std::endl;
    }
    else
    {
        std::vector<GiNaC::ex> Xc;
        if(branch1SPtr==nullptr)cout<<"branch1 is null"<<endl;


        Xc.push_back(branch1SPtr->x_equation);
        Xc.push_back(branch1SPtr->y_equation);
        Xc.push_back(branch1SPtr->z_equation);
        Xc.push_back(branch1SPtr->w_equation);
        
        h += Xc[0] * Xc[0] * bigintToNumeric(Qu(0, 0));     // xx
        h += 2 * Xc[0] * Xc[1] * bigintToNumeric(Qu(0, 1)); // xy
        h += 2 * Xc[0] * Xc[2] * bigintToNumeric(Qu(0, 2)); // xz
        h += 2 * Xc[0] * Xc[3] * bigintToNumeric(Qu(0, 3)); // x
        h += Xc[1] * Xc[1] * bigintToNumeric(Qu(1, 1));     // yy
        h += 2 * Xc[1] * Xc[2] * bigintToNumeric(Qu(1, 2)); // yz
        h += 2 * Xc[1] * Xc[3] * bigintToNumeric(Qu(1, 3)); // y
        h += Xc[2] * Xc[2] * bigintToNumeric(Qu(2, 2));     // zz
        h += 2 * Xc[2] * Xc[3] * bigintToNumeric(Qu(2, 3)); // z
        h += Xc[3] * Xc[3] * bigintToNumeric(Qu(3, 3));     // w

        std::cout << "h(ξ)= " << h << std::endl;
        if(branch1SPtr->have_delta)
        {
            h=h.subs(branch1SPtr->Delta==branch1SPtr->delta_equation);
            std::cout << "h(ξ) with delta = " << h << std::endl;
        }
    }

    std::cout<<"intersectionTriple h done"<<std::endl;

    if (h.is_zero())
    {
        return ans;
    }
    std::vector<TQSI::Root> roots;
    roots = std::move(findTQSI(h));

    std::cout<<"intersectionTriple find root done"<<std::endl;
    
    // if(beginEnd.first!=beginEnd.second)
    // {
    //     h= findIntersectionPointSmoothQuartic((*GlobalVariable::getAllCurveEquation())[beginEnd.first],Qu);
    // }
    // else
    // {

    // }

    std::cout << "print roots" << std::endl;
    for (int i = 0; i < roots.size(); ++i)
    {
        std::cout << roots[i].parameter << " multiplicity = " << roots[i].multiplicity << std::endl;
        auto parameter = roots[i].parameter;
        if (parameter != INF_PARAMETER)
        {
            auto point=branch1SPtr->getPoint(parameter, 1, 1);
            if(point.x==INFINITY)
            {
                point=branch1SPtr->getPoint(-parameter, 1, 1);
                if(point.x==INFINITY)continue;//image root
                else
                {
                    if(isPointOnSurface(point,Qu))
                    {
                        branch1SPtr->pointList.emplace_back(point,-parameter);
                    }
                    else if(beginEnd.first!=beginEnd.second)
                    {
                        auto point2=branch2SPtr->getPoint(-parameter, 1, 1);
                        if(isPointOnSurface(point2,Qu))
                        {
                            branch2SPtr->pointList.emplace_back(point2,-parameter);
                        }
                    }
                    
                }
            }
            else
            {
                if(isPointOnSurface(point,Qu))
                {
                    branch1SPtr->pointList.emplace_back(point,parameter);
                }
                else if(beginEnd.first+1!=beginEnd.second)
                {
                    auto point2=branch2SPtr->getPoint(parameter, 1, 1);
                    if(isPointOnSurface(point2,Qu))
                    {
                        branch2SPtr->pointList.emplace_back(point2,parameter);
                    }
                }
            }
        }
        else
        {
            
            if(isPointOnSurface(branch1SPtr->negativeInf,Qu))
            ans.points.emplace_back(branch1SPtr->negativeInf, parameter);
            if(isPointOnSurface(branch1SPtr->positiveInf,Qu))
            ans.points.emplace_back(branch1SPtr->positiveInf, parameter);
            if(beginEnd.first+1!=beginEnd.second)
            {
                if(isPointOnSurface(branch2SPtr->negativeInf,Qu))
                ans.points.emplace_back(branch2SPtr->negativeInf, parameter);
                if(isPointOnSurface(branch2SPtr->positiveInf,Qu))
                ans.points.emplace_back(branch2SPtr->positiveInf, parameter);
            }
        }
    }
    std::sort(ans.points.begin(), ans.points.end(), [&](TQSI::intPoint a, TQSI::intPoint b)
              { return a.parameter < b.parameter; });

    // for(auto s:ans.points)
    // {
    //     cout<<"parameter= "<<s.parameter<<endl;
    // }

    std::cout<<"intersectionTriple : done"<<std::endl;
    return ans;
}


TQSI::injectiveMap matchSeq(TQSI::seq seq1, TQSI::seq seq2)
{
    TQSI::injectiveMap ans;
    return ans;
}

void matchTwoSeq(TQSI::seq &seq1, TQSI::seq &seq2, double tolerance)
{
    int idxA=0;
    int k=0;
    while(idxA<seq1.points.size()&&k<1000)
    {
        int cnt=0;
        std::vector<unsigned int>intersect;
        for(unsigned int i=0;i<seq2.points.size();i++)
        {
            if(seq1.points[idxA].point.isEqual(seq2.points[i].point,tolerance))
            {
                intersect.push_back(i);
            }
        }

        if(intersect.size()==1)
        {
            auto idxB=intersect[0];
            if(seq1.points[idxA].parameter==INF_PARAMETER||seq1.points[idxA].parameter==-INF_PARAMETER)
            {
                if(seq2.points[idxB].parameter!=INF_PARAMETER&&seq2.points[idxB].parameter!=-INF_PARAMETER)
                {
                    seq1.points[idxA].point=seq2.points[idxB].point;
                }
            }
            else if(seq2.points[idxB].parameter==INF_PARAMETER||seq2.points[idxB].parameter==-INF_PARAMETER)
            {
               if(seq1.points[idxA].parameter!=INF_PARAMETER&&seq1.points[idxA].parameter!=-INF_PARAMETER)
               {
                    seq2.points[idxB].point=seq1.points[idxA].point;
               }
            }

            idxA++;
        }
        else if(intersect.empty())
        {
            tolerance*=10;
        }
        else
        {
            tolerance/=10;
        }
        k++;
    }
}

//取参数为非无穷的点
void matchThreeSeq(TQSI::seq &seq1, TQSI::seq &seq2, TQSI::seq &seq3)
{
    double tolerance = 1e-4;

    matchTwoSeq(seq1,seq2,tolerance);
    matchTwoSeq(seq1,seq3,tolerance);
    matchTwoSeq(seq2,seq3,tolerance);
}

void curveAddSeqToPointList(curve_equation& intersectCurve,const TQSI::seq& seq)
{
    for(auto& intPoint:seq.points)
    {
        intersectCurve.pointList.push_back(intPoint);
    }
}

void divideCurve(curve_equation& intersectCurve,CSGtree* tree)
{
    std::cout<<"divideCurve"<<std::endl;

    std::vector<double> parameters;
    bool addClose=false;
    auto& seq=intersectCurve.pointList;
    if(seq.empty())
    {
        std::cout<<"seq is empty"<<std::endl;
        if(intersectCurve.isAllInterval)
        {
            //deal ring edge   
        }
        return;
    }
    
    intersectCurve.isAllInterval=false;

    for(auto& s:seq)
    {
        parameters.push_back(s.parameter);
    }
    if(intersectCurve.isLoop)
    {
        if(seq[0].parameter==-INF_PARAMETER)
        {
            parameters.push_back(INF_PARAMETER);
            addClose=true;
        }
    }
    int numOfPtsOnCurve = parameters.size();
    bool isCurSegValid = false;
    bool startPushedOverlap = false;

    // cout<<"collect parameter"<<endl;
    // for(auto para:parameters)
    // {
    //     cout<<para<<endl;
    // }

    for (int i = 0; i < numOfPtsOnCurve - 1; i++)
    {
        bool isPrevSegValid = isCurSegValid;
        TQSI::point midPoint=intersectCurve.getPoint((parameters[i]+parameters[i+1])/2,1,1);
        isCurSegValid=tree->isPointOn(midPoint);
        if(isCurSegValid)
        {
            intersectCurve.validSegment.emplace_back(parameters[i],parameters[i+1]);
        }

        if (i == 0)
        {
            startPushedOverlap = isCurSegValid;
        }
        else
        {
            // Todo: there is problem when here several overlap segments for surface-curve
            if (!isPrevSegValid && !isCurSegValid)
            {
                intersectCurve.validSegment.emplace_back(parameters[i],parameters[i]);
            }
        }
    }
    bool endPushedOverlap = isCurSegValid;

    if(addClose)
    {
        if (!startPushedOverlap && !endPushedOverlap)
        {
            // (i == 0) and (i == numOfPtsOnCurve - 1) is surface curve overlap start(end) point
            bool isValid=tree->isPointOn((intersectCurve.negativeInf+intersectCurve.positiveInf)/2,1e-4);
            // cout<<"sdf= "<<tree->get_sdf(intersectCurve.negativeInf);
            // cout<<"sdf= "<<tree->get_sdf((intersectCurve.negativeInf+intersectCurve.positiveInf)/2);
            if (isValid)
            {
                intersectCurve.validSegment.emplace_back(INF_PARAMETER,INF_PARAMETER);
            }
        }
    }
    else
    {
        if (!startPushedOverlap)
        {
            // i == 0 is surface curve overlap start point
            bool isValid=tree->isPointOn(intersectCurve.getPoint(parameters[0],1,1));
            if (isValid)
            {
                intersectCurve.validSegment.emplace_back(parameters[0],parameters[0]);
            }
        }
        if (!endPushedOverlap)
        {
            // i == numOfPtsOnCurve - 1 is surface curve overlap end point
            bool isValid=tree->isPointOn(intersectCurve.getPoint(parameters.back(),1,1));
            if (isValid)
            {
                intersectCurve.validSegment.emplace_back(parameters.back(),parameters.back());
            }
        }
    }

    for(auto segment:intersectCurve.validSegment)
    {
        cout<<"startParameter= "<<segment.startParameter<<" endParameter= "<<segment.endParameter<<endl;
    }

    std::cout<<"divideCurve finish"<<std::endl;
    //interPoints list to segment
    intersectCurve.pointList.clear();

    //according segment save result to pool
    for(auto segment:intersectCurve.validSegment)
    {

    }

}

void divideCurve(curve_equation& intersectCurve,CSGtree* tree,TQSI::seq& seq)
{
    std::cout<<"divideCurve"<<std::endl;

    std::vector<double> parameters;
    bool addClose=false;
    for(auto& s:seq.points)
    {
        parameters.push_back(s.parameter);
    }
    if(intersectCurve.isLoop)
    {
        if(seq.points[0].parameter==-INF_PARAMETER)
        {
            parameters.push_back(INF_PARAMETER);
            addClose=true;
        }
    }
    int numOfPtsOnCurve = parameters.size();
    bool isCurSegValid = false;
    bool startPushedOverlap = false;

    // cout<<"collect parameter"<<endl;
    // for(auto para:parameters)
    // {
    //     cout<<para<<endl;
    // }

    for (int i = 0; i < numOfPtsOnCurve - 1; i++)
    {
        bool isPrevSegValid = isCurSegValid;
        TQSI::point midPoint=intersectCurve.getPoint((parameters[i]+parameters[i+1])/2,1,1);
        isCurSegValid=tree->isPointOn(midPoint);
        if(isCurSegValid)
        {
            intersectCurve.validSegment.emplace_back(parameters[i],parameters[i+1]);
        }

        if (i == 0)
        {
            startPushedOverlap = isCurSegValid;
        }
        else
        {
            // Todo: there is problem when here several overlap segments for surface-curve
            if (!isPrevSegValid && !isCurSegValid)
            {
                intersectCurve.validSegment.emplace_back(parameters[i],parameters[i]);
                // m_result->m_intPoints.push_back(faceCurvePts[i]);
            }
        }
    }
    bool endPushedOverlap = isCurSegValid;

    if(addClose)
    {
        if (!startPushedOverlap && !endPushedOverlap)
        {
            // (i == 0) and (i == numOfPtsOnCurve - 1) is surface curve overlap start(end) point
            bool isValid=tree->isPointOn(intersectCurve.negativeInf);
            if (isValid)
            {
                intersectCurve.validSegment.emplace_back(INF_PARAMETER,INF_PARAMETER);
            }
        }
        // else
        // {
        //     faceCurvePts[0].m_intPoint->m_prevRelation = faceCurvePts[numOfPtsOnCurve - 1].m_intPoint->m_prevRelation;
        // }
    }
    else
    {
        if (!startPushedOverlap)
        {
            // i == 0 is surface curve overlap start point
            bool isValid=tree->isPointOn(intersectCurve.getPoint(parameters[0],1,1));
            if (isValid)
            {
                // filterRedundantPoint(m_result->m_intPoints, faceCurvePts[0]);
                intersectCurve.validSegment.emplace_back(parameters[0],parameters[0]);
            }
        }
        if (!endPushedOverlap)
        {
            // i == numOfPtsOnCurve - 1 is surface curve overlap end point
            bool isValid=tree->isPointOn(intersectCurve.getPoint(parameters.back(),1,1));
            if (isValid)
            {
                intersectCurve.validSegment.emplace_back(parameters.back(),parameters.back());
            }
        }
    }

    for(auto segment:intersectCurve.validSegment)
    {
        cout<<"startParameter= "<<segment.startParameter<<" endParameter= "<<segment.endParameter<<endl;
    }

    std::cout<<"divideCurve finish"<<std::endl;
}







