#ifndef _qi_h_
#define _qi_h_

#define PRINT_VAR(value) std::cout << #value << " = " << (value) << std::endl;
#define INF_PARAMETER 100000

/** Top level header to be included in the case
    you intend to use both the calculus kernel
    and the IO kit for parsing user input and
    pretty display results. */

#include <glm/glm.hpp>
#include <ginac/ginac.h>
#include <igl/marching_cubes.h>
#include <igl/voxel_grid.h>
#include <string>
#include <vector>

#include <unsupported/Eigen/Polynomials>

#include "qi_io.h"
#include "qi_kernel.h"

#include"lange/Global.h"
#include "lange/AABB.h"
#include"lange/Global.h"

struct vertex_equation;
struct curve_equation;
struct edge_equation;



struct GinacParams
{
    GiNaC::symbol u;
    GiNaC::symbol Delta;
    GiNaC::parser *reader;
    GiNaC::symtab *table;
    GinacParams(){};
    GinacParams(GiNaC::symbol _u, GiNaC::symbol _Delta, GiNaC::symtab *_table) : u(_u), Delta(_Delta), table(_table){};
};
namespace TQSI
{
    struct point
    {
        double x = 0;
        double y = 0;
        double z = 0;
        point()
        {
            x = 0;
            y = 0;
            z = 0;
        }
        point(double _x, double _y, double _z)
            : x(_x), y(_y), z(_z){};
        // bool isEqual(point a,point b,double tol=1e-8)
        // {

        // }
        double sqrDistance(const point &ptOther) const
        {
            double aD = 0;
            double aDD = 0;

            aDD = ptOther.x;
            aDD -= (*this).x;
            aDD *= aDD;
            aD += aDD;
            aDD = ptOther.y;
            aDD -= (*this).y;
            aDD *= aDD;
            aD += aDD;
            aDD = ptOther.z;
            aDD -= (*this).z;
            aDD *= aDD;
            aD += aDD;
            return aD;
        }
        double &operator[](int n)
        {
            assert(n >= 0 && n <= 2);
            return (&x)[n];
        }

        const double &operator[](int n) const
        {
            assert(n >= 0 && n <= 2);
            return (&x)[n];
        }
        friend std::ostream &operator<<(std::ostream &lhs, const point rhs)
        {
            lhs << "[" << rhs.x << "," << rhs.y << "," << rhs.z << "]";
            return lhs;
        }
        point operator+(const point &vec)
        {
            return point(this->x + vec.x, this->y + vec.y, this->z + vec.z);
        }
        point operator-(const point &vec)
        {
            return point(this->x - vec.x, this->y - vec.y, this->z - vec.z);
        }
        point operator/(const int a)
        {
            return point(this->x / a, this->y / a, this->z / a);
        }

        point operator/(const double a)
        {
            return point(this->x / a, this->y / a, this->z / a);
        }
        bool isEqual(const point &ptOther, const double tol = 1e-8 /*= PSGMConstants::LengthTolerance*/) const
        {
            return sqrDistance(ptOther) <= tol * tol;
        }
    };

    struct intPoint
    {
        TQSI::point point;
        double parameter;
        intPoint(double _x, double _y, double _z, double _p)
            : point(_x, _y, _z), parameter(_p){};
        intPoint(TQSI::point _point, double _p) : point(_point), parameter(_p){};

        bool operator==(const intPoint &other)
        {
            return parameter == other.parameter;
        }
        void print()
        {
            cout << "x= " << point.x << " y= " << point.y << " z= " << point.z << " parameter= " << parameter << endl;
        }
    };

    struct interval
    {
        double startParameter, endParameter;
        interval(double _s, double _e) : startParameter(_s), endParameter(_e){};

        bool contains(const double t, const double eps = 1e-8)
        {
            if (t >= startParameter - eps && t <= endParameter + eps)
            {
                return true;
            }
            return false;
        }
    };
}

struct edge_equation
{
    std::shared_ptr<curve_equation>  curveSPtr;
    std::vector<TQSI::interval> validSegment;
};

struct vertex_equation
{
    double x, y, z;
    unsigned int first_index, second_index, third_index;
    std::vector<curve_equation *> components;
    std::vector<double> parameters;
    std::vector<vertex_equation *> nextVertex;
    std::vector<vertex_equation *> preVertex;
    vertex_equation(double _x, double _y, double _z, unsigned int fi, unsigned int si, unsigned int ti)
        : x(_x), y(_y), z(_z), first_index(fi), second_index(si), third_index(ti){};
    vertex_equation(double _x, double _y, double _z): x(_x), y(_y), z(_z){};
    vertex_equation(Eigen::MatrixXd point)
    {
        x = point(0);
        y = point(1);
        z = point(2);
    }
    vertex_equation(TQSI::point point)
    {
        x = point.x;
        y = point.y;
        z = point.z;
    }
};

struct curve_equation
{
    GiNaC::ex x_equation, y_equation, z_equation, w_equation, delta_equation;
    GiNaC::ex x_diff, y_diff, z_diff, w_diff; // now dont need diff information
    GiNaC::ex surfacePrametrizations[4];
    GiNaC::symbol u, Delta;
    GiNaC::symbol v, s,t;
    bool have_delta = false;
    unsigned int first_index, second_index;
    std::vector<vertex_equation *> verticesList;
    std::vector<TQSI::intPoint> pointList;
    std::vector<TQSI::intPoint> endPointList;
    bool isLoop=false;
    bool isAllInterval = true;
    unsigned int visualizaitonPoint_cnt = 100;
    std::vector<TQSI::interval> validSegment;
    TQSI::point negativeInf, positiveInf;

    curve_equation(std::vector<GiNaC::ex> equ,
                   GiNaC::symbol v,
                   GiNaC::symbol d,
                   unsigned int fi,
                   unsigned int si,
                   std::vector<GiNaC::ex> _surfacePrametrizations)
    {
        first_index = fi;
        second_index = si;
        x_equation = equ[0];
        y_equation = equ[1];
        z_equation = equ[2];
        w_equation = equ[3];

        // (弃用)直接符号处理可能会导致速度下降，但精度会上升
        // 后续需要单独的x,y,z,w equation
        // x_equation /= w_equation;
        // y_equation /= w_equation;
        // z_equation /= w_equation;

        u = v;
        Delta = d;
        if (equ.size() == 5)
        {
            delta_equation = equ[4];
            have_delta = true;
            isAllInterval=false;

            surfacePrametrizations[0]=_surfacePrametrizations[0];
            surfacePrametrizations[1]=_surfacePrametrizations[1];
            surfacePrametrizations[2]=_surfacePrametrizations[2];
            surfacePrametrizations[3]=_surfacePrametrizations[3];

            int degree = delta_equation.degree(u);
            static Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
            // Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
            Eigen::VectorXd coefficients(degree + 1);
            for (int i = 0; i <= degree; i++)
            {
                auto co = GiNaC::ex_to<GiNaC::numeric>(delta_equation.coeff(u, i).evalf()).to_double();
                coefficients(i) = co;
            }
            // cout<<"coefficients= "<<coefficients<<endl;
            solver.compute(coefficients);
            const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType &gen = solver.roots();

            std::vector<double> roots;
            roots.push_back(-INF_PARAMETER);
            for(int i=0;i<gen.size();i++)
            {
                if(!gen(i).imag())
                {
                    roots.push_back(gen(i).real());
                }
            }
            sort(roots.begin(),roots.end());
            roots.push_back(INF_PARAMETER);
            for (size_t i = 1; i < roots.size(); ++i) 
            {
                double left = roots[i-1];
                double right = roots[i];
                double testPoint;

                // 选择测试点
                if (left == -INF_PARAMETER) {
                    testPoint = right - 1; // 选择接近无穷小的点
                } else if (right == INF_PARAMETER) {
                    testPoint = left + 1; // 选择接近无穷大的点
                } else {
                    testPoint = (left + right) / 2; // 选择区间中点
                }

                // 计算测试点的函数值
                if(GiNaC::ex_to<GiNaC::numeric>(delta_equation.subs(u==testPoint).evalf()).is_positive())
                {
                    validSegment.emplace_back(left,right);
                    double wLeft = GiNaC::ex_to<GiNaC::numeric>(w_equation.subs(GiNaC::lst{u == left, Delta == 0}).evalf()).to_double();
                    if (std::abs(wLeft) >= 1e-8)
                    {
                        GiNaC::ex xdivw = x_equation / w_equation;
                        GiNaC::ex ydivw = y_equation / w_equation;
                        GiNaC::ex zdivw = z_equation / w_equation;
                        double xLeft = GiNaC::ex_to<GiNaC::numeric>(xdivw.subs(GiNaC::lst{u == left, Delta == 0}).evalf()).to_double();
                        double yLeft = GiNaC::ex_to<GiNaC::numeric>(ydivw.subs(GiNaC::lst{u == left, Delta == 0}).evalf()).to_double();
                        double zLeft = GiNaC::ex_to<GiNaC::numeric>(zdivw.subs(GiNaC::lst{u == left, Delta == 0}).evalf()).to_double();
                        endPointList.emplace_back(xLeft,yLeft,zLeft,left);
                    }
                    else
                    {
                        endPointList.emplace_back(INFINITY,INFINITY,INFINITY,left);
                    }

                    double wRight = GiNaC::ex_to<GiNaC::numeric>(w_equation.subs(GiNaC::lst{u == right, Delta == 0}).evalf()).to_double();
                    if (std::abs(wRight) >= 1e-8)
                    {
                        GiNaC::ex xdivw = x_equation / w_equation;
                        GiNaC::ex ydivw = y_equation / w_equation;
                        GiNaC::ex zdivw = z_equation / w_equation;
                        double xRight = GiNaC::ex_to<GiNaC::numeric>(xdivw.subs(GiNaC::lst{u == right, Delta == 0}).evalf()).to_double();
                        double yRight = GiNaC::ex_to<GiNaC::numeric>(ydivw.subs(GiNaC::lst{u == right, Delta == 0}).evalf()).to_double();
                        double zRight = GiNaC::ex_to<GiNaC::numeric>(zdivw.subs(GiNaC::lst{u == right, Delta == 0}).evalf()).to_double();
                        endPointList.emplace_back(xRight,yRight,zRight,right);
                    }
                    else
                    {
                        endPointList.emplace_back(INFINITY,INFINITY,INFINITY,right);
                    }
                    // 打印区间和符号
                    std::cout << "Interval (" << left << ", " << right << "): f(x) is " << sign << std::endl;
                }
                
            }

        }

        // this->negativeInf=this->getPoint(-1000,1,1);
        // this->positiveInf=this->getPoint(1000,1,1);
        negativeInf = getPoint(-INF_PARAMETER, 1, 0);
        positiveInf = getPoint(INF_PARAMETER, 1, 0);

        if (negativeInf.isEqual(positiveInf, 1e-4))
        {
            isLoop = true;
        }

        GiNaC::ex dx;
        GiNaC::ex dy;
        GiNaC::ex dz;
        GiNaC::ex dw;
        if (have_delta)
        {
            dx = x_equation.subs(Delta == delta_equation);
            dy = y_equation.subs(Delta == delta_equation);
            dz = z_equation.subs(Delta == delta_equation);
            dw = w_equation.subs(Delta == delta_equation);
        }
        else
        {
            dx = x_equation;
            dy = y_equation;
            dz = z_equation;
            dw = w_equation;
        }

        // cout << "==================\n";
        // cout << dx << endl;
        // cout << dx.subs(u == u_equation) << endl;
        // cout << dx.subs(u == u_equation).diff(u) << endl;

        
        x_diff = dx.diff(u);
        y_diff = dy.diff(u);
        z_diff = dz.diff(u);
        w_diff = dw.diff(u);
    }

    curve_equation(std::vector<GiNaC::ex> equ,
                   unsigned int fi,
                   unsigned int si,
                   std::vector<GiNaC::ex> _surfacePrametrizations):curve_equation(equ,GlobalVariable::getU(),GlobalVariable::getDelta(),fi,si,_surfacePrametrizations)
    {
        // curve_equation(equ,GlobalVariable::getU(),GlobalVariable::getDelta(),fi,si,_surfacePrametrizations);
    }               

    glm::dvec3 getPoint(double t)
    {
        double j = tan(t / 2);
        double x, y, z, w;
        GiNaC::ex xdivw = x_equation / w_equation;
        GiNaC::ex ydivw = y_equation / w_equation;
        GiNaC::ex zdivw = z_equation / w_equation;

        if (have_delta)
        {
            double delta = GiNaC::ex_to<GiNaC::numeric>(delta_equation.subs(u == t).evalf()).to_double();
            w = GiNaC::ex_to<GiNaC::numeric>(w_equation.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
            if (std::abs(w) <= 1e-8)
            {
                return {0, 0, 0};
            }
            x = GiNaC::ex_to<GiNaC::numeric>(xdivw.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
            y = GiNaC::ex_to<GiNaC::numeric>(ydivw.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
            z = GiNaC::ex_to<GiNaC::numeric>(zdivw.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
        }
        else
        {
            w = GiNaC::ex_to<GiNaC::numeric>(w_equation.subs(u == t).evalf()).to_double();
            if (std::abs(w) <= 1e-8)
            {
                return {0, 0, 0};
            }
            x = GiNaC::ex_to<GiNaC::numeric>(xdivw.subs(u == t).evalf()).to_double();
            y = GiNaC::ex_to<GiNaC::numeric>(ydivw.subs(u == t).evalf()).to_double();
            z = GiNaC::ex_to<GiNaC::numeric>(zdivw.subs(u == t).evalf()).to_double();
        }
        return {x, y, z};
        // return {x / w, y / w, z / w};
    }

    glm::dvec3 getPoint(double t, bool a)
    {
        // double j = tan(t / 2);
        double x, y, z, w;
        GiNaC::ex xdivw = x_equation / w_equation;
        GiNaC::ex ydivw = y_equation / w_equation;
        GiNaC::ex zdivw = z_equation / w_equation;
        if (have_delta)
        {
            double delta = GiNaC::ex_to<GiNaC::numeric>(delta_equation.subs(u == t).evalf()).to_double();
            w = GiNaC::ex_to<GiNaC::numeric>(w_equation.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
            if (std::abs(w) <= 1e-8)
            {
                return {0, 0, 0};
            }
            x = GiNaC::ex_to<GiNaC::numeric>(xdivw.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
            y = GiNaC::ex_to<GiNaC::numeric>(ydivw.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
            z = GiNaC::ex_to<GiNaC::numeric>(zdivw.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
        }
        else
        {
            w = GiNaC::ex_to<GiNaC::numeric>(w_equation.subs(u == t).evalf()).to_double();
            if (std::abs(w) <= 1e-8)
            {
                return {0, 0, 0};
            }
            x = GiNaC::ex_to<GiNaC::numeric>(xdivw.subs(u == t).evalf()).to_double();
            y = GiNaC::ex_to<GiNaC::numeric>(ydivw.subs(u == t).evalf()).to_double();
            z = GiNaC::ex_to<GiNaC::numeric>(zdivw.subs(u == t).evalf()).to_double();
        }
        return {x, y, z};
        // return {x / w, y / w, z / w};
    }

    TQSI::point getPoint(double t, bool a, bool b)
    {
        if (b)
        {
            if (t == INF_PARAMETER)
                return positiveInf;
            if (t == -INF_PARAMETER)
                return negativeInf;
            for(const auto& endPoint:endPointList)
            {
                if(t==endPoint.parameter)
                {
                    return endPoint.point;
                }
            }
        }
        
        bool isValid=false;
        for(const auto& segment:validSegment)
        {
            if(segment.startParameter<=t&&t<=segment.endParameter)
            {
                isValid=true;
                break;
            }
        }
        if(!isValid)
        {   
            return {INFINITY, INFINITY, INFINITY};
        }

        // double j = tan(t / 2);
        double x, y, z, w;
        GiNaC::ex xdivw = x_equation / w_equation;
        GiNaC::ex ydivw = y_equation / w_equation;
        GiNaC::ex zdivw = z_equation / w_equation;
        if (have_delta)
        {
            double delta = GiNaC::ex_to<GiNaC::numeric>(delta_equation.subs(u == t).evalf()).to_double();
            w = GiNaC::ex_to<GiNaC::numeric>(w_equation.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
            if (std::abs(w) <= 1e-8)
            {
                return {INFINITY, INFINITY, INFINITY};
            }
            x = GiNaC::ex_to<GiNaC::numeric>(xdivw.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
            y = GiNaC::ex_to<GiNaC::numeric>(ydivw.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
            z = GiNaC::ex_to<GiNaC::numeric>(zdivw.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
        }
        else
        {
            w = GiNaC::ex_to<GiNaC::numeric>(w_equation.subs(u == t).evalf()).to_double();
            if (std::abs(w) <= 1e-8)
            {
                return {0, 0, 0};
            }
            x = GiNaC::ex_to<GiNaC::numeric>(xdivw.subs(u == t).evalf()).to_double();
            y = GiNaC::ex_to<GiNaC::numeric>(ydivw.subs(u == t).evalf()).to_double();
            z = GiNaC::ex_to<GiNaC::numeric>(zdivw.subs(u == t).evalf()).to_double();
        }
        return {x, y, z};
        // return {x / w, y / w, z / w};
    }

    glm::dvec4 get4DPoint(double t)
    {
        double j = tan(t / 2);
        double x, y, z, w;
        GiNaC::ex xdivw = x_equation / w_equation;
        GiNaC::ex ydivw = y_equation / w_equation;
        GiNaC::ex zdivw = z_equation / w_equation;
        if (have_delta)
        {
            double delta = GiNaC::ex_to<GiNaC::numeric>(delta_equation.subs(u == t).evalf()).to_double();
            w = GiNaC::ex_to<GiNaC::numeric>(w_equation.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
            if (std::abs(w) <= 1e-8)
            {
                return {0, 0, 0, 0};
            }
            x = GiNaC::ex_to<GiNaC::numeric>(xdivw.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
            y = GiNaC::ex_to<GiNaC::numeric>(ydivw.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
            z = GiNaC::ex_to<GiNaC::numeric>(zdivw.subs(GiNaC::lst{u == t, Delta == delta}).evalf()).to_double();
        }
        else
        {
            w = GiNaC::ex_to<GiNaC::numeric>(w_equation.subs(u == t).evalf()).to_double();
            if (std::abs(w) <= 1e-8)
            {
                return {0, 0, 0, 0};
            }
            x = GiNaC::ex_to<GiNaC::numeric>(xdivw.subs(u == t).evalf()).to_double();
            y = GiNaC::ex_to<GiNaC::numeric>(ydivw.subs(u == t).evalf()).to_double();
            z = GiNaC::ex_to<GiNaC::numeric>(zdivw.subs(u == t).evalf()).to_double();
        }
        return {x, y, z, w};
    }

    glm::dvec4 get4DPoint(double t, bool a)
    {

        double x, y, z, w;
        if (have_delta)
        {
            double delta = GiNaC::ex_to<GiNaC::numeric>(delta_equation.subs(u == t)).to_double();
            x = GiNaC::ex_to<GiNaC::numeric>(x_equation.subs(GiNaC::lst{u == t, Delta == delta})).to_double();
            y = GiNaC::ex_to<GiNaC::numeric>(y_equation.subs(GiNaC::lst{u == t, Delta == delta})).to_double();
            z = GiNaC::ex_to<GiNaC::numeric>(z_equation.subs(GiNaC::lst{u == t, Delta == delta})).to_double();
            w = GiNaC::ex_to<GiNaC::numeric>(w_equation.subs(GiNaC::lst{u == t, Delta == delta})).to_double();
        }
        else
        {
            x = GiNaC::ex_to<GiNaC::numeric>(x_equation.subs(u == t)).to_double();
            y = GiNaC::ex_to<GiNaC::numeric>(y_equation.subs(u == t)).to_double();
            z = GiNaC::ex_to<GiNaC::numeric>(z_equation.subs(u == t)).to_double();
            w = GiNaC::ex_to<GiNaC::numeric>(w_equation.subs(u == t)).to_double();
        }
        return {x, y, z, w};
    }

    glm::dvec4 getDerivatives(double t)
    {
        // GiNaC::symbol u("u");
        // GiNaC::symbol y("y");
        // GiNaC::ex x = pow(u, 2) + sqrt(y);
        // GiNaC::ex y_expr = pow(u, 3);

        // // 将y表达式代入x方程中
        // GiNaC::ex x_with_y = x.subs(y == y_expr);

        // // 输出结果
        // std::cout << "x(u) with y(u) substituted: " << x_with_y << std::endl;

        // std::cout<<x_equation<<std::endl;
        // std::cout<<x_equation.diff(u)<<std::endl;
        // std::cout<<x_equation.diff(u).subs(u==t)<<std::endl;

        // if(have_delta)
        // {
        //   dx=x_equation.subs(Delta==delta_equation);
        //   dy=y_equation.subs(Delta==delta_equation);
        //   dz=z_equation.subs(Delta==delta_equation);
        //   dw=w_equation.subs(Delta==delta_equation);
        // }
        // else
        // {
        //   dx=x_equation;
        //   dy=y_equation;
        //   dz=z_equation;
        //   dw=w_equation;
        // }

        // double j = tan(t / 2);
        //  std::cout << x_equation << std::endl;
        //  std::cout << x_diff << std::endl;
        //  std::cout << "j= " << j << std::endl;
        //  std::cout << x_diff.subs(u == j) << std::endl;
        //  std::cout<<dy<<std::endl;
        //  std::cout<<dy.diff(u).subs(u==t)<<std::endl;
        //  std::cout<<dz<<std::endl;
        //  std::cout<<dz.diff(u).subs(u==t)<<std::endl;
        //  std::cout<<dw<<std::endl;
        //  std::cout<<dw.diff(u).subs(u==t)<<std::endl;

        glm::dvec4 ans(GiNaC::ex_to<GiNaC::numeric>(x_diff.subs(u == t)).to_double(),
                       GiNaC::ex_to<GiNaC::numeric>(y_diff.subs(u == t)).to_double(),
                       GiNaC::ex_to<GiNaC::numeric>(z_diff.subs(u == t)).to_double(),
                       GiNaC::ex_to<GiNaC::numeric>(w_diff.subs(u == t)).to_double());
        return ans;
    }

    glm::dvec4 getSecondDerivatives(double t)
    {
        // double j = tan(t / 2);
        // std::cout << "second\n";
        // std::cout << x_diff << std::endl;
        // std::cout << x_diff.diff(u) << std::endl;
        // std::cout << x_diff.diff(u).subs(u == t) << std::endl;
        glm::dvec4 ans(GiNaC::ex_to<GiNaC::numeric>(x_diff.diff(u).subs(u == t)).to_double(),
                       GiNaC::ex_to<GiNaC::numeric>(y_diff.diff(u).subs(u == t)).to_double(),
                       GiNaC::ex_to<GiNaC::numeric>(z_diff.diff(u).subs(u == t)).to_double(),
                       GiNaC::ex_to<GiNaC::numeric>(w_diff.diff(u).subs(u == t)).to_double());
        return ans;
    }

    void print()
    {
        PRINT_VAR(first_index);
        PRINT_VAR(second_index);
        PRINT_VAR(validSegment.size());
        for (int i = 0; i < validSegment.size(); i++)
        {
            cout << "segment " << i << " start= " << validSegment[i].startParameter << " end= " << validSegment[i].endParameter << endl;
        }
        // PRINT_VAR(x_equation);
        // PRINT_VAR(y_equation);
        // PRINT_VAR(z_equation);
        // PRINT_VAR(w_equation);
        // std::cout << x_equation << std::endl;
        // std::cout << y_equation << std::endl;
        // std::cout << z_equation << std::endl;
        // std::cout << w_equation << std::endl;
        // std::cout << x_diff << std::endl;
        // std::cout << y_diff << std::endl;
        // std::cout << z_diff << std::endl;
        // std::cout << w_diff << std::endl;
        // std::cout << x_diff.diff(u) << std::endl;
        // std::cout << y_diff.diff(u) << std::endl;
        // std::cout << z_diff.diff(u) << std::endl;
        // std::cout << w_diff.diff(u) << std::endl;
    }
};

// class langeBVH : public PSGMIntBVH
// {
// public:
//     curve_equation* PolygonCurve;
//     void buildCurveBVH(int depth, int index, double t)
//     {
//         PSGMBBox box;
//         PSGMIntBVHNode node;
//         double resolution = 1 << (depth - 1);

//         double tstep = 2 * M_PI / resolution;

//         node.m_depth = depth;
//         node.m_param.emplace_back(t, t + tstep);

//         if (depth == getHeight())
//         {
//             //box.expand(m_pCurve->getPoint(t));
//             //box.expand(m_pCurve->getPoint(t + tstep));

//             box.expand(PolygonCurve->getPoint(t));
//             box.expand(PolygonCurve->getPoint(t + tstep));

//             node.m_firstChild = -1;
//             node.m_box = box;
//             m_data[index] = node;
//             return;
//         }
//         else
//         {
//             int firstIndex = 2 * index + 1;
//             buildCurveBVH(depth + 1, firstIndex, t);
//             buildCurveBVH(depth + 1, firstIndex + 1, t + tstep / 2.0);

//             box.merge(m_data[firstIndex].m_box);
//             box.merge(m_data[firstIndex + 1].m_box);

//             node.m_firstChild = firstIndex;
//             node.m_box = box;
//             m_data[index] = node;
//         }
//     }

//     /**
//      * @brief       Construct a BVH tree for a NurbsCurve
//      * @param[in]   bvhHeight   the height of the bvh
//      * @param[in]   pCur  the object that the bvh contains
//     */
//     langeBVH(int bvhHeight, curve_equation* pCur)
//     {
//         m_bvhHeight = bvhHeight;
//         PolygonCurve = pCur;
//         m_dataSize = (1 << m_bvhHeight) - 1;
//         m_data.resize(getSize());
//         m_type = PSGMIntBVHType::Curve;
//     }
// };

bigint_matrix vector2quadric(const Eigen::Matrix<float, 10, 1> qvec);

/**@brief 把传入的浮点数列向量全部同乘一个系数，使其转化成整数列向量
 * @param[in] qvec 传入的列向量
 * @param[in] mask_length 掩码长度，即丢弃浮点位数后mask_length位
 */
void vcf2vcbigint(Eigen::Matrix<float, 10, 1> &qvec, short mask_length);

mpz_class bigint_gcd(mpz_class a, mpz_class b);

bigint_matrix vcf2gmp2bigintmatrix(Eigen::Matrix<float, 10, 1> &qvec);

void trim(string &s);

int trim_deng(string &s);

string delay_no_doubt(string &s);

std::vector<string> spilt4ginac(string &str, string &delta);

std::vector<GiNaC::ex> str2ginacex(std::vector<string> spilt);

double quadric_sdf_function(double x, double y, double z, Eigen::Matrix<float, 10, 1> v);

void handleCurveComponent(const QIOutputInfo *outputInfo, const GinacParams &ginacParams,
                          std::vector<curve_equation> &all_curve_equation,
                          std::vector<vertex_equation> &all_vertex_equation, const unsigned int first_index, const unsigned int second_index);

void handleCurveComponent(const QIOutputInfo *outputInfo,
                          const unsigned int first_index, const unsigned int second_index);

std::pair<unsigned int,unsigned int> intersectionDouble(const std::vector<bigint_matrix> &Q, std::map<pair<unsigned int,unsigned int>,pair<unsigned int,unsigned int>> &firstSecondToBeginEnd,const unsigned int first_index, const unsigned int second_index, const GinacParams &ginacParams,
                        std::vector<curve_equation> &all_curve_equation,
                        std::vector<vertex_equation> &all_vertex_equation);

std::pair<unsigned int,unsigned int> intersectionDouble(const std::vector<bigint_matrix> &Q, std::map<pair<unsigned int,unsigned int>,pair<unsigned int,unsigned int>> &firstSecondToBeginEnd,const unsigned int first_index, const unsigned int second_index);






#endif
