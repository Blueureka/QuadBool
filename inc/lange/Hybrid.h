#pragma once

#include "PSUtils.h"

#include <CSG/csg.h>
#include <unordered_map>
#include <libqi/qi.h>

#include<lange/Global.h>
#include "lange/AABB.h"

//TODO::make cone


struct PSGMShape
{
    PSGMTopoBody *body{nullptr};
    psgmstd::vector<PSGMTopoLump *> lumps{};
    psgmstd::vector<PSGMTopoShell *> shells{};
    psgmstd::vector<PSGMTopoFace *> faces{};
    psgmstd::vector<PSGMTopoLoop *> loops{};
    psgmstd::vector<PSGMTopoEdge *> edges{};
    psgmstd::vector<PSGMTopoCoedge *> coedges{};
    psgmstd::vector<PSGMTopoVertex *> vertices{};
};

class hybridPresentation
{
public:
    PSGMShape shape;
    // CSGtree *tree{nullptr};
    primitive* tree{nullptr};
    // std::unordered_map<PSGMTopoFace *, Eigen::Matrix<float, 10, 1>> *facePtr2impliceFunction{nullptr};
    // std::unordered_map<Eigen::Matrix<float, 10, 1>, PSGMTopoFace *> *impliceFunction2facePtr{nullptr};
    std::vector<PSGMTopoFace *> *facesPointer{nullptr};
    std::vector<Eigen::Matrix<float, 10, 1>> *functionsPointer{nullptr};
    std::map<PSGMTopoFace *,unsigned int> *faceQMatrixIndex{nullptr};

    std::vector<PSGMTopoEdge *> *edgesPointer{nullptr};
    std::map<PSGMTopoEdge *,unsigned int> *edgeFuctionIndex{nullptr};
    AABB box;

    hybridPresentation()
    {
        facesPointer = new std::vector<PSGMTopoFace *>;
        functionsPointer = new std::vector<Eigen::Matrix<float, 10, 1>>;
        faceQMatrixIndex = new std::map<PSGMTopoFace *,unsigned int>;
        edgeFuctionIndex =new std::map<PSGMTopoEdge *,unsigned int>;
    }
    unsigned int getQMatriIndex(PSGMTopoFace* face)
    {
        return (*faceQMatrixIndex)[face];
    }
protected:
    PSGMLCS lcs;
};

class box_hybrid : public hybridPresentation
{
    double conerX, conerY, conerZ;
    double heighX, heighY, heighZ;
public:
    inline void makeBox(PSGMISessionContext *context, double cx, double cy, double cz,
                        double hx, double hy, double hz, PSGMShape &block);

    box_hybrid(PSGMISessionContext *context, double cx, double cy, double cz,
               double hx, double hy, double hz);
    
};

class cylinder_hybrid : public hybridPresentation
{
    double conerX, conerY, conerZ;
    double heighZ;
    double radius;
    // PSGMLCS lcs;
public:
    inline void makeCylinder(PSGMISessionContext *context, double cx, double cy, double cz,
                             double hz, double radius, PSGMShape &cylinder);

    cylinder_hybrid(PSGMISessionContext *context, double cx, double cy, double cz,
                    double hz, double r);
    
    inline void makeCylinder(PSGMISessionContext *context, PSGMLCS lcs,
                             double height, double radius, PSGMShape &cylinder);
    cylinder_hybrid(PSGMISessionContext *context, PSGMLCS lcs,
                    double hz, double r);
};

class sphere_hybrid : public hybridPresentation
{
    double conerX, conerY, conerZ;
    double radius;
    // PSGMLCS lcs;
public:
    inline void makeSphere(PSGMISessionContext *context, double cx, double cy, double cz,
                           double radius, PSGMShape &sphere);

    sphere_hybrid(PSGMISessionContext *context, double cx, double cy, double cz,
                  double r);
};

inline PSGMTopoVertex *addVertex(PSGMShape &shape, const PSGMPoint &point)
{
    auto vertex = shape.body->createVertex(point);
    shape.vertices.emplace_back(vertex);
    return vertex;
}

inline PSGMTopoEdge *addEdge(PSGMShape &shape, PSGMCurveSPtr curve)
{
    auto entCurve = convertCurveToEntity(shape.body->getContext(), curve);
    auto edge = shape.body->createEdge(entCurve);
    shape.edges.emplace_back(edge);
    return edge;
}

inline PSGMTopoEdge *
addEdge(PSGMShape &shape, PSGMCurveSPtr curve, PSGMTopoVertex *bgn, PSGMTopoVertex *end, bool sense = true)
{
    auto entCurve = convertCurveToEntity(shape.body->getContext(), curve);
    auto edge = shape.body->createEdge(bgn, end, entCurve, sense);
    shape.edges.emplace_back(edge);
    return edge;
}

inline PSGMTopoCoedge *addCoedgeToLoop(PSGMShape &shape, PSGMTopoLoop *loop, PSGMTopoEdge *edge, bool sense = true)
{
    auto coedge = edge->createCoedge(sense);
    shape.coedges.emplace_back(coedge);
    loop->appendCoedge(coedge);
    return coedge;
}

inline PSGMEntCurve *createLine(PSGMISessionContext *context, PSGMPoint const &pt, PSGMDirection const &dir)
{
    PSGMCurveLineSPtr const &line{std::shared_ptr<PSGMCurveLine>(PSGMCurveLine::createByDir(pt, dir))};
    return PSGMEntCurveLine::create(context, line);
}

inline PSGMEntCurve *createLine(PSGMISessionContext *context, PSGMPoint const &pt1, PSGMPoint const &pt2)
{
    PSGMCurveLineSPtr const &line{std::shared_ptr<PSGMCurveLine>(PSGMCurveLine::createByPoints(pt1, pt2))};
    return PSGMEntCurveLine::create(context, line);
}

inline PSGMEntSurface *createPlane(PSGMISessionContext *context, PSGMPoint const &pt, PSGMDirection const &dir)
{
    PSGMSurfPlaneSPtr const &plane{std::shared_ptr<PSGMSurfPlane>(PSGMSurfPlane::create(pt, dir))};
    return PSGMEntSurfPlane::create(context, plane);
}

inline PSGMEntCurve *createCircle(PSGMISessionContext *context, const PSGMLCS &lcs, const double radius)
{
    return PSGMEntCurveCircle::create(context, PSGMCurveCircle::createByRadius(lcs, radius));
}

inline PSGMEntCurve *
createEllipse(PSGMISessionContext *context, const PSGMLCS &lcs, const double majorRadius, const double minorRadius)
{
    return PSGMEntCurveEllipse::create(context, PSGMCurveEllipse::createByRadius(lcs, majorRadius, minorRadius));
}
