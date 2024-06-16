#pragma once

#include "geometry/PSGMLCS.h"
#include "lange/TQSI.h"
#include "lange/AABB.h"

double op_intersection(double d1, double d2);
double op_union(double d1, double d2);
double op_substraction(double d1, double d2);

class CSGtree;



class primitive
{
    unsigned int id;
    const char *name;//maybe use std::string

public:
    virtual double sdf(double x, double y, double z) = 0;
};

enum TREEtype
{
    PRIMITIVE = 0,
    OPERATION,
};
enum CSGoperation
{
    UNION = 0,
    INTERSECTION,
    SUBTRACTION,
    REVERSESUBTRACTION
};
class CSGtree
{
    CSGtree *left, *right;
    unsigned short type;
    unsigned short operation;
    primitive *p;
    double minX, minY, minZ, maxX, maxY, maxZ;

public:
    CSGtree()
    {
        left = NULL;
        right = NULL;
        type = PRIMITIVE;
    }
    CSGtree(const CSGtree& other)
    {
        left=other.left;
        right=other.right;
        type=other.type;
        operation=other.operation;
        p=other.p;
    }
    CSGtree(primitive *_p)
        : p(_p)
    {
        left = NULL;
        right = NULL;
        type = PRIMITIVE;
    }
    CSGtree(int op, primitive *A, primitive *B)
        : operation(op), type(OPERATION)
    {
        CSGtree *l = new CSGtree(A);
        CSGtree *r = new CSGtree(B);
        left = l;
        right = r;
    }
    CSGtree(int op, CSGtree *A, CSGtree *B)
        : operation(op), left(A), right(B), type(OPERATION)
    {
    }
    CSGtree(int op, primitive *A, CSGtree *B)
        : operation(op), right(B), type(OPERATION)
    {
        CSGtree *l = new CSGtree(A);
        left = l;
    }
    CSGtree(int op, CSGtree *A, primitive *B)
        : operation(op), left(A), type(OPERATION)
    {
        CSGtree *r = new CSGtree(B);
        right = r;
    }
    
    unsigned short getType()
    {
        return this->type;
    }

    void print()
    {
        std::cout<<((int)type==0?"PRIMITIVE":"OPERATION")<<std::endl;
        std::cout<<((int)operation==0?"UNION":"no UNION")<<std::endl;
        std::cout<<(left==nullptr?"nullptr":"have left")<<std::endl;
        std::cout<<(right==nullptr?"nullptr":"have right")<<std::endl;
        std::cout<<(p==nullptr?"nullptr":"have p")<<std::endl;
    }

    double get_sdf(double x, double y, double z);
    double get_sdf(TQSI::point p);

    bool isPointInside(TQSI::point p,double tolerance=1e-8);
    bool isPointInside(double x,double y,double z,double tolerance=1e-8);

    bool isPointOn(TQSI::point p,double tolerance=1e-8);
    bool isPointOn(double x,double y,double z,double tolerance=1e-8);

    bool isPointOutside(TQSI::point p,double tolerance=1e-8);
    bool isPointOutside(double x,double y,double z,double tolerance=1e-8);
};

class sphere_primitive : public primitive
{
    double centerX, centerY, centerZ, radius;

public:
    sphere_primitive(double cx, double cy, double cz, double r)
        : centerX(cx), centerY(cy), centerZ(cz), radius(r)
    {
    }
    sphere_primitive(PSGMLCS lcs, double r);
    double sdf(double x, double y, double z) override;
};

class box_primitive : public primitive
{
    double conerX, conerY, conerZ;
    double heighX, heighY, heighZ;

public:
    box_primitive(double cx, double cy, double cz,
                  double hx, double hy, double hz)
        : conerX(cx), conerY(cy), conerZ(cz), heighX(hx), heighY(hy), heighZ(hz)
    {
    }
    double sdf(double x, double y, double z) override;
};

class cylinder_primitive : public primitive
{
    double centerX, centerY, centerZ, radius;
    double heighX, heighY, heighZ;

public:
    cylinder_primitive(double cx, double cy, double cz,
                       double hx, double hy, double hz, double r)
        : centerX(cx), centerY(cy), centerZ(cz), heighX(hx), heighY(hy), heighZ(hz), radius(r)
    {
    }
    cylinder_primitive(double cx, double cy, double cz,
                       double hz, double r)
        : centerX(cx), centerY(cy), centerZ(cz), heighZ(hz), radius(r)
    {
    }
    double sdf(double x, double y, double z) override;
};

class cone_primitive : public primitive
{
    double centerX, centerY, centerZ, radius;
    double heighX, heighY, heighZ;

public:
    cone_primitive(double cx, double cy, double cz,
                   double hx, double hy, double hz, double r)
        : centerX(cx), centerY(cy), centerZ(cz), heighX(hx), heighY(hy), heighZ(hz), radius(r)
    {
    }
    cone_primitive(double cx, double cy, double cz,
                   double hz, double r)
        : centerX(cx), centerY(cy), centerZ(cz), heighZ(hz), radius(r)
    {
    }
    double sdf(double x, double y, double z) override;
};
