#pragma once

#include <igl/opengl/glfw/Viewer.h>

#include <libqi/qi.h>
#include <igl/marching_cubes.h>
#include <igl/voxel_grid.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

#include "CSG/csg.h"

struct RGB 
{
    int r;
    int g;
    int b;
};

namespace TQSI
{

class Viewer
{
public:
    Viewer()=delete;

    Viewer(igl::opengl::glfw::Viewer* v);

    void forDeBug();

    void addSurface(const Eigen::Matrix<float, 10, 1>& surface);
    void addSurface(const std::vector<Eigen::Matrix<float, 10, 1>>& surfaces);

    void addCurve(const curve_equation& curve);
    void addCurve(const std::vector<curve_equation>& curves);

    void addPoint(const Eigen::MatrixXd& point);
    void addPoint(const std::vector<Eigen::MatrixXd>& points);

    void addCSG(const CSGtree csg);
    void addCSG(const std::vector<CSGtree>& csgs);

    void drawSurface();
    void saveMeshToOff(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const std::string& filename);

    void sampleCurvePoint(const curve_equation& eqution,Eigen::MatrixXd& mat,Eigen::MatrixXd& color,const double step);
    void sampleCurvePoint(const curve_equation& eqution,Eigen::MatrixXd& mat,const double step);

    bool isParameterInValidSegement(const curve_equation& eqution,const double parameter);
    void sampleCurvePointWithOneInterval(const curve_equation& eqution,Eigen::MatrixXd& mat,const double step);
    void sampleCurvePointWithIntervals(curve_equation& eqution,std::vector<Eigen::MatrixXd>& mat,const double point_cnt);

    void drawCurve();
    void drawCurveByPoint();
    
    void drawPoint();

    void drawCSG();

    void visualizate();

    

    void setVisualizateSurface(bool b);
    void setVisualizateCurve(bool b);
    void setVisualizatePoint(bool b);
    void setVisualizateCSG(bool b);

    

private:
    igl::opengl::glfw::Viewer* viewerPointer;
    std::vector<Eigen::Matrix<float, 10, 1>> surfaceList;
    std::vector<curve_equation> curveList;
    std::vector<Eigen::MatrixXd> pointList;
    std::vector<CSGtree> CSGList;
    bool visualizateSurface=false;
    bool visualizateCurve=false;
    bool visualizatePoint=false;
    bool visualizateCSG=false;
};

}


//色相、饱和度、亮度
RGB HSVtoRGB(float H, float S, float V);

std::vector<RGB> generateColors(int n);