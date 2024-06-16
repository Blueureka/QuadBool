#include <CSG/csg.h>

double op_intersection(double d1, double d2)
{
    return d1 > d2 ? d1 : d2;
}

double op_union(double d1, double d2)
{
    return d1 < d2 ? d1 : d2;
}

double op_substraction(double d1, double d2)
{
    return d1 > -d2 ? d1 : -d2;
}

bool CSGtree::isPointInside(TQSI::point p,double tolerance)
{
    return this->get_sdf(p.x,p.y,p.z) <= -tolerance;
}
bool CSGtree::isPointInside(double x,double y,double z,double tolerance)
{
    return this->get_sdf(x,y,z) <= -tolerance;
}

bool CSGtree::isPointOn(TQSI::point p,double tolerance)
{
    double sdf=this->get_sdf(p.x,p.y,p.z);
    return sdf>-tolerance && sdf<tolerance;
}
bool CSGtree::isPointOn(double x,double y,double z,double tolerance)
{
    double sdf=this->get_sdf(x,y,z);
    return sdf>-tolerance && sdf<tolerance;
}

bool CSGtree::isPointOutside(TQSI::point p,double tolerance)
{
    return this->get_sdf(p.x,p.y,p.z) >= -tolerance;
}
bool CSGtree::isPointOutside(double x,double y,double z,double tolerance)
{
    return this->get_sdf(x,y,z) >= -tolerance;
}


double sphere_primitive::sdf(double x, double y, double z)
{
    return (x - this->centerX) * (x - this->centerX) + (y - this->centerY) * (y - this->centerY) + (z - this->centerZ) * (z - this->centerZ) - this->radius * this->radius;
}

sphere_primitive::sphere_primitive(PSGMLCS lcs, double r)
    : radius(r)
{
    auto p = lcs.getOrigin();
    centerX = p[0];
    centerY = p[1];
    centerZ = p[2];
}

double box_primitive::sdf(double x, double y, double z)
{
    double bottom = conerZ - z;
    double up = z - conerZ - heighZ;
    double left = conerY - y;
    double right = y - conerY - heighY;
    double front = conerX - x;
    double back = x - conerX - heighX;
    return op_intersection(op_intersection(op_intersection(bottom, up), op_intersection(left, right)), op_intersection(front, back));
}

double cylinder_primitive::sdf(double x, double y, double z)
{
    double d1 = (x - this->centerX) * (x - this->centerX) + (y - this->centerY) * (y - this->centerY) - this->radius * this->radius;
    double d2 = centerZ - z;
    double d3 = z - centerZ - heighZ;
    return op_intersection(op_intersection(d2, d3), d1);
}

double cone_primitive::sdf(double x, double y, double z)
{
    double d1 = z + (1.0 / this->radius) * sqrt((x - this->centerX) * (x - this->centerX) + (y - this->centerY) * (y - this->centerY)) - this->centerZ;
    double d2 = this->centerZ - heighZ - z;

    return op_intersection(d2, d1);
}

double CSGtree::get_sdf(double x, double y, double z)
{
    if (this->type == PRIMITIVE)
    {
        return this->p->sdf(x, y, z);
    }
    double l = 0, r = 0;
    l = this->left->get_sdf(x, y, z);
    r = this->right->get_sdf(x, y, z);
    if (operation == UNION)
    {
        return op_union(l, r);
    }
    else if (operation == INTERSECTION)
    {
        return op_intersection(l, r);
    }
    else if (operation == SUBTRACTION)
    {
        return op_substraction(l, r);
    }
}

double CSGtree::get_sdf(TQSI::point p)
{
    return get_sdf(p.x,p.y,p.z);
}