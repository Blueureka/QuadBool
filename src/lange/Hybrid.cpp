#include <lange/Hybrid.h>

inline void box_hybrid::makeBox(PSGMISessionContext *context, double cx, double cy, double cz,
                                double hx, double hy, double hz, PSGMShape &block)
{
    block.body = PSGMTopoBody::create(context, PSGMBodyType::Solid);
    PSGMTopoLump *lump{block.body->createLump(PSGMLumpType::Solid)};
    block.lumps.push_back(lump);

    PSGMTopoShell *shell{lump->createShell()};
    block.shells.push_back(shell);

    PSGMPoint minPt{cx, cy, cz};
    PSGMPoint maxPt{cx + hx, cy + hy, cz + hz};
    psgmstd::vector<PSGMPoint> pts{};
    pts.emplace_back(minPt.x(), minPt.y(), minPt.z());
    pts.emplace_back(maxPt.x(), minPt.y(), minPt.z());
    pts.emplace_back(maxPt.x(), maxPt.y(), minPt.z());
    pts.emplace_back(minPt.x(), maxPt.y(), minPt.z());
    pts.emplace_back(minPt.x(), minPt.y(), maxPt.z());
    pts.emplace_back(maxPt.x(), minPt.y(), maxPt.z());
    pts.emplace_back(maxPt.x(), maxPt.y(), maxPt.z());
    pts.emplace_back(minPt.x(), maxPt.y(), maxPt.z());

    size_t id{0};
    for (id = 0; id < 8; ++id)
    {
        block.vertices.emplace_back(block.body->createVertex(pts[id]));
    }

    psgmstd::vector<PSGMEntSurface *> surfs{};
    surfs.emplace_back(createPlane(context, minPt, PSGMDirection{0.0, 0.0, -1.0}));
    surfs.emplace_back(createPlane(context, maxPt, PSGMDirection{0.0, 0.0, 1.0}));
    surfs.emplace_back(createPlane(context, minPt, PSGMDirection{0.0, -1.0, 0.0}));
    surfs.emplace_back(createPlane(context, maxPt, PSGMDirection{1.0, 0.0, 0.0}));
    surfs.emplace_back(createPlane(context, maxPt, PSGMDirection{0.0, 1.0, 0.0}));
    surfs.emplace_back(createPlane(context, minPt, PSGMDirection{-1.0, 0.0, 0.0}));

    for (id = 0; id < 6; ++id)
    {
        block.faces.emplace_back(shell->createFace(surfs[id]));
    }

    // add map
    Eigen::Matrix<float, 10, 1> tempFunction;
    PSGMDirection direction;
    auto createMapVector = [&](PSGMDirection direction, PSGMPoint point, Eigen::Matrix<float, 10, 1> &tempFunction)
    {
        double dx = direction.x();
        double dy = direction.y();
        double dz = direction.z();
        double x = point.x();
        double y = point.y();
        double z = point.z();
        tempFunction << 0, 0, 0, dx, 0, 0, dy, 0, dz, -(x * dx + y * dy + z * dz);
    };

    auto v=GlobalVariable::getImplicitMatrixVector();

    direction.set(0.0, 0.0, -1.0);
    createMapVector(direction, minPt, tempFunction);
    // (*facePtr2impliceFunction)[block.faces[0]] = tempFunction;
    // (*impliceFunction2facePtr)[tempFunction] = block.faces[0];
    (*facesPointer).emplace_back(block.faces[0]);
    (*functionsPointer).emplace_back(tempFunction);
    v->push_back(tempFunction);
    

    direction.set(0.0, 0.0, 1.0);
    createMapVector(direction, maxPt, tempFunction);
    // (*facePtr2impliceFunction)[block.faces[1]] = tempFunction;
    // (*impliceFunction2facePtr)[tempFunction] = block.faces[1];
    // facesPointer.emplace_back(block.faces[1]);
    // functionsPointer.emplace_back(tempFunction);
    (*facesPointer).emplace_back(block.faces[1]);
    (*functionsPointer).emplace_back(tempFunction);
    v->push_back(tempFunction);

    direction.set(0.0, -1.0, 0.0);
    createMapVector(direction, minPt, tempFunction);
    // (*facePtr2impliceFunction)[block.faces[2]] = tempFunction;
    // (*impliceFunction2facePtr)[tempFunction] = block.faces[2];
    // facesPointer.emplace_back(block.faces[2]);
    // functionsPointer.emplace_back(tempFunction);
    (*facesPointer).emplace_back(block.faces[2]);
    (*functionsPointer).emplace_back(tempFunction);
    v->push_back(tempFunction);

    direction.set(1.0, 0.0, 0.0);
    createMapVector(direction, maxPt, tempFunction);
    // (*facePtr2impliceFunction)[block.faces[3]] = tempFunction;
    // (*impliceFunction2facePtr)[tempFunction] = block.faces[3];
    // facesPointer.emplace_back(block.faces[3]);
    // functionsPointer.emplace_back(tempFunction);
    (*facesPointer).emplace_back(block.faces[3]);
    (*functionsPointer).emplace_back(tempFunction);
    v->push_back(tempFunction);

    direction.set(0.0, 1.0, 0.0);
    createMapVector(direction, maxPt, tempFunction);
    // (*facePtr2impliceFunction)[block.faces[4]] = tempFunction;
    // (*impliceFunction2facePtr)[tempFunction] = block.faces[4];
    // facesPointer.emplace_back(block.faces[4]);
    // functionsPointer.emplace_back(tempFunction);
    (*facesPointer).emplace_back(block.faces[4]);
    (*functionsPointer).emplace_back(tempFunction);
    v->push_back(tempFunction);

    direction.set(-1.0, 0.0, 0.0);
    createMapVector(direction, minPt, tempFunction);
    // (*facePtr2impliceFunction)[block.faces[5]] = tempFunction;
    // (*impliceFunction2facePtr)[tempFunction] = block.faces[5];
    // facesPointer.emplace_back(block.faces[5]);
    // functionsPointer.emplace_back(tempFunction);
    (*facesPointer).emplace_back(block.faces[5]);
    (*functionsPointer).emplace_back(tempFunction);
    v->push_back(tempFunction);

    psgmstd::vector<PSGMEntCurve *> lines{};
    lines.emplace_back(createLine(context, pts[0], PSGMDirection{1.0, 0.0, 0.0}));
    lines.emplace_back(createLine(context, pts[1], PSGMDirection{0.0, 1.0, 0.0}));
    lines.emplace_back(createLine(context, pts[2], PSGMDirection{-1.0, 0.0, 0.0}));
    lines.emplace_back(createLine(context, pts[3], PSGMDirection{0.0, -1.0, 0.0}));
    lines.emplace_back(createLine(context, pts[4], PSGMDirection{1.0, 0.0, 0.0}));
    lines.emplace_back(createLine(context, pts[5], PSGMDirection{0.0, 1.0, 0.0}));
    lines.emplace_back(createLine(context, pts[6], PSGMDirection{-1.0, 0.0, 0.0}));
    lines.emplace_back(createLine(context, pts[7], PSGMDirection{0.0, -1.0, 0.0}));
    lines.emplace_back(createLine(context, pts[0], PSGMDirection{0.0, 0.0, 1.0}));
    lines.emplace_back(createLine(context, pts[1], PSGMDirection{0.0, 0.0, 1.0}));
    lines.emplace_back(createLine(context, pts[2], PSGMDirection{0.0, 0.0, 1.0}));
    lines.emplace_back(createLine(context, pts[3], PSGMDirection{0.0, 0.0, 1.0}));

    block.edges.emplace_back(block.body->createEdge(block.vertices[0], block.vertices[1], lines[0]));
    block.edges.emplace_back(block.body->createEdge(block.vertices[1], block.vertices[2], lines[1]));
    block.edges.emplace_back(block.body->createEdge(block.vertices[2], block.vertices[3], lines[2]));
    block.edges.emplace_back(block.body->createEdge(block.vertices[3], block.vertices[0], lines[3]));
    block.edges.emplace_back(block.body->createEdge(block.vertices[4], block.vertices[5], lines[4]));
    block.edges.emplace_back(block.body->createEdge(block.vertices[5], block.vertices[6], lines[5]));
    block.edges.emplace_back(block.body->createEdge(block.vertices[6], block.vertices[7], lines[6]));
    block.edges.emplace_back(block.body->createEdge(block.vertices[7], block.vertices[4], lines[7]));
    block.edges.emplace_back(block.body->createEdge(block.vertices[0], block.vertices[4], lines[8]));
    block.edges.emplace_back(block.body->createEdge(block.vertices[1], block.vertices[5], lines[9]));
    block.edges.emplace_back(block.body->createEdge(block.vertices[2], block.vertices[6], lines[10]));
    block.edges.emplace_back(block.body->createEdge(block.vertices[3], block.vertices[7], lines[11]));

    for (id = 0; id < 6; ++id)
    {
        block.loops.emplace_back(block.faces[id]->createLoop());
    }

    for (id = 0; id < 12; ++id)
    {
        block.coedges.emplace_back(block.edges[id]->createCoedge(true));
        block.coedges.emplace_back(block.edges[id]->createCoedge(false));
    }

    block.loops[0]->appendCoedge(block.coedges[1]);
    block.loops[0]->appendCoedge(block.coedges[7]);
    block.loops[0]->appendCoedge(block.coedges[5]);
    block.loops[0]->appendCoedge(block.coedges[3]);

    block.loops[1]->appendCoedge(block.coedges[8]);
    block.loops[1]->appendCoedge(block.coedges[10]);
    block.loops[1]->appendCoedge(block.coedges[12]);
    block.loops[1]->appendCoedge(block.coedges[14]);

    block.loops[2]->appendCoedge(block.coedges[0]);
    block.loops[2]->appendCoedge(block.coedges[18]);
    block.loops[2]->appendCoedge(block.coedges[9]);
    block.loops[2]->appendCoedge(block.coedges[17]);

    block.loops[3]->appendCoedge(block.coedges[2]);
    block.loops[3]->appendCoedge(block.coedges[20]);
    block.loops[3]->appendCoedge(block.coedges[11]);
    block.loops[3]->appendCoedge(block.coedges[19]);

    block.loops[4]->appendCoedge(block.coedges[4]);
    block.loops[4]->appendCoedge(block.coedges[22]);
    block.loops[4]->appendCoedge(block.coedges[13]);
    block.loops[4]->appendCoedge(block.coedges[21]);

    block.loops[5]->appendCoedge(block.coedges[6]);
    block.loops[5]->appendCoedge(block.coedges[16]);
    block.loops[5]->appendCoedge(block.coedges[15]);
    block.loops[5]->appendCoedge(block.coedges[23]);
}

box_hybrid::box_hybrid(PSGMISessionContext *context, double cx, double cy, double cz,
                       double hx, double hy, double hz)
    : conerX(cx), conerY(cy), conerZ(cz), heighX(hx), heighY(hy), heighZ(hz)
{
    box_primitive *newBox = new box_primitive(conerX, conerY, conerZ, heighX, heighY, heighZ);
    // tree = new CSGtree(newBox);
    tree=(primitive*)newBox;
    facesPointer = new std::vector<PSGMTopoFace *>;
    functionsPointer = new std::vector<Eigen::Matrix<float, 10, 1>>;
    makeBox(context, conerX, conerY, conerZ, heighX, heighY, heighZ, shape);
}

void cylinder_hybrid::makeCylinder(PSGMISessionContext *context, double cx, double cy, double cz,
                                   double hz, double radius, PSGMShape &cylinder)
{
    double height = hz;

    cylinder.body = PSGMTopoBody::create(context, PSGMBodyType::Solid);
    PSGMTopoLump *lump{cylinder.body->createLump(PSGMLumpType::Solid)};
    cylinder.lumps.push_back(lump);

    PSGMTopoShell *shell{lump->createShell()};
    cylinder.shells.push_back(shell);

    const PSGMLCS lcs(PSGMPoint{cx, cy, cz}, PSGMDirection{1.0, 0, 0}, PSGMDirection{0, 1.0, 0});
    PSGMLCS const lcs1{lcs.getOrigin(), lcs.getXDir(), -lcs.getYDir()};
    PSGMLCS const lcs2{lcs.getOrigin() + lcs.getZDir() * height, lcs.getXDir(), lcs.getYDir()};

    PSGMSurfPlaneSPtr const &plane1{std::shared_ptr<PSGMSurfPlane>(PSGMSurfPlane::create(lcs1))};
    PSGMSurfPlaneSPtr const &plane2{std::shared_ptr<PSGMSurfPlane>(PSGMSurfPlane::create(lcs2))};
    PSGMSurfCylinderSPtr cylnd{nullptr};
    PSGMSurfCylinder::create(lcs, radius, &cylnd);

    PSGMCurveCircleSPtr const &circ1{std::shared_ptr<PSGMCurveCircle>(PSGMCurveCircle::createByRadius(lcs1, radius))};
    PSGMCurveCircleSPtr const &circ2{std::shared_ptr<PSGMCurveCircle>(PSGMCurveCircle::createByRadius(lcs2, radius))};

    PSGMEntSurface *entPlane1{PSGMEntSurfPlane::create(context, plane1)};
    PSGMEntSurface *entPlane2{PSGMEntSurfPlane::create(context, plane2)};
    PSGMEntSurface *entCylinder{PSGMEntSurfCylinder::create(context, cylnd)};

    PSGMEntCurve *circle1{PSGMEntCurveCircle::create(context, circ1)};
    PSGMEntCurve *circle2{PSGMEntCurveCircle::create(context, circ2)};

    cylinder.faces.emplace_back(cylinder.shells[0]->createFace(entPlane1));
    cylinder.faces.emplace_back(cylinder.shells[0]->createFace(entPlane2));
    cylinder.faces.emplace_back(cylinder.shells[0]->createFace(entCylinder));

    cylinder.edges.emplace_back(cylinder.body->createEdge(circle1));
    cylinder.edges.emplace_back(cylinder.body->createEdge(circle2));

    auto v=GlobalVariable::getImplicitMatrixVector();
    auto Q=GlobalVariable::getQMatrix();
    auto qIndex=Q->size();
    auto firstSecondToBeginEnd=GlobalVariable::getFirstSecondToBeginEnd();
    auto qMatrixIndexToPSGMsurfaceSPtrSPtr=GlobalVariable::getQMatrixIndexToPSGMsurfaceSPtrSPtr();

    // add map
    Eigen::Matrix<float, 10, 1> tempFunction;
    PSGMDirection direction;
    auto createMapVector = [&](PSGMLCS lcs, Eigen::Matrix<float, 10, 1> &tempFunction)
    {
        PSGMPoint point = lcs.getOrigin();
        PSGMDirection direction = lcs.getZDir();
        double dx = direction.x();
        double dy = direction.y();
        double dz = direction.z();
        double x = point.x();
        double y = point.y();
        double z = point.z();
        tempFunction << 0, 0, 0, dx, 0, 0, dy, 0, dz, -(x * dx + y * dy + z * dz);
    };

    createMapVector(lcs1, tempFunction);
    facesPointer->emplace_back(cylinder.faces[0]);
    functionsPointer->emplace_back(tempFunction);
    v->push_back(tempFunction);
    (*faceQMatrixIndex)[cylinder.faces[0]]=Q->size();
    Q->push_back(vcf2gmp2bigintmatrix(tempFunction));
    qMatrixIndexToPSGMsurfaceSPtrSPtr->push_back(cylinder.faces[0]->getSurface());

    createMapVector(lcs2, tempFunction);
    (*facesPointer).emplace_back(cylinder.faces[1]);
    (*functionsPointer).emplace_back(tempFunction);
    v->push_back(tempFunction);
    (*faceQMatrixIndex)[cylinder.faces[1]]=Q->size();
    Q->push_back(vcf2gmp2bigintmatrix(tempFunction));
    qMatrixIndexToPSGMsurfaceSPtrSPtr->push_back(cylinder.faces[1]->getSurface());

    {
        PSGMPoint point = lcs.getOrigin();
        PSGMDirection direction = lcs.getZDir();
        double dx = direction.x();
        double dy = direction.y();
        double dz = direction.z();
        double x = point.x();
        double y = point.y();
        double z = point.z();
        tempFunction << 1, 0, 0, 0, 1, 0, 0, 0, 0, -(radius * radius);

        (*facesPointer).emplace_back(cylinder.faces[2]);
        (*functionsPointer).emplace_back(tempFunction);
        v->push_back(tempFunction);
        (*faceQMatrixIndex)[cylinder.faces[2]]=Q->size();
        Q->push_back(vcf2gmp2bigintmatrix(tempFunction));
        qMatrixIndexToPSGMsurfaceSPtrSPtr->push_back(cylinder.faces[2]->getSurface());
    }

    (*edgeFuctionIndex)[cylinder.edges[0]]=GlobalVariable::getAllCurveEquation()->size();
    intersectionDouble(*Q,*firstSecondToBeginEnd,qIndex,qIndex+2,*GlobalVariable::getGinacParamsSPtr(),*GlobalVariable::getAllCurveEquation(),*GlobalVariable::getAllVertexEquation());

    (*edgeFuctionIndex)[cylinder.edges[1]]=GlobalVariable::getAllCurveEquation()->size();
    intersectionDouble(*Q,*firstSecondToBeginEnd,qIndex+1,qIndex+2,*GlobalVariable::getGinacParamsSPtr(),*GlobalVariable::getAllCurveEquation(),*GlobalVariable::getAllVertexEquation());
    

    cylinder.loops.emplace_back(cylinder.faces[0]->createLoop(PSGMLoopType::Outer));
    cylinder.loops.emplace_back(cylinder.faces[1]->createLoop(PSGMLoopType::Outer));
    cylinder.loops.emplace_back(cylinder.faces[2]->createLoop(PSGMLoopType::Winding));
    cylinder.loops.emplace_back(cylinder.faces[2]->createLoop(PSGMLoopType::Winding));

    cylinder.loops[0]->appendCoedge(cylinder.edges[0]->createCoedge(true));
    cylinder.loops[1]->appendCoedge(cylinder.edges[1]->createCoedge(true));

    cylinder.loops[2]->appendCoedge(cylinder.edges[0]->createCoedge(false));
    cylinder.loops[3]->appendCoedge(cylinder.edges[1]->createCoedge(false));
}

cylinder_hybrid::cylinder_hybrid(PSGMISessionContext *context, double cx, double cy, double cz,
                                 double hz, double r)
    : conerX(cx), conerY(cy), conerZ(cz), heighZ(hz), radius(r)
{
    cylinder_primitive *newCylinder = new cylinder_primitive(conerX, conerY, conerZ, heighZ, radius);
    // tree = new CSGtree(newCylinder);
    tree=(primitive*)newCylinder;
    // facesPointer = new std::vector<PSGMTopoFace *>;
    // functionsPointer = new std::vector<Eigen::Matrix<float, 10, 1>>;
    makeCylinder(context, conerX, conerY, conerZ, heighZ, radius, shape);
}


inline void cylinder_hybrid::makeCylinder(PSGMISessionContext *context, PSGMLCS lcs,
                             double height, double radius, PSGMShape &cylinder)
{
    cylinder.body = PSGMTopoBody::create(context, PSGMBodyType::Solid);
    PSGMTopoLump* lump{cylinder.body->createLump(PSGMLumpType::Solid)};
    cylinder.lumps.push_back(lump);

    PSGMTopoShell* shell{lump->createShell()};
    cylinder.shells.push_back(shell);

    PSGMLCS const lcs1{lcs.getOrigin(), lcs.getXDir(), -lcs.getYDir()};
    PSGMLCS const lcs2{lcs.getOrigin() + lcs.getZDir() * height, lcs.getXDir(), lcs.getYDir()};

    PSGMSurfPlaneSPtr const& plane1{std::shared_ptr<PSGMSurfPlane>(PSGMSurfPlane::create(lcs1))};
    PSGMSurfPlaneSPtr const& plane2{std::shared_ptr<PSGMSurfPlane>(PSGMSurfPlane::create(lcs2))};
    PSGMSurfCylinderSPtr cylnd{nullptr};
    PSGMSurfCylinder::create(lcs, radius, &cylnd);

    PSGMCurveCircleSPtr const& circ1{std::shared_ptr<PSGMCurveCircle>(PSGMCurveCircle::createByRadius(lcs1, radius))};
    PSGMCurveCircleSPtr const& circ2{std::shared_ptr<PSGMCurveCircle>(PSGMCurveCircle::createByRadius(lcs2, radius))};

    PSGMEntSurface* entPlane1{PSGMEntSurfPlane::create(context, plane1)};
    PSGMEntSurface* entPlane2{PSGMEntSurfPlane::create(context, plane2)};
    PSGMEntSurface* entCylinder{PSGMEntSurfCylinder::create(context, cylnd)};

    PSGMEntCurve* circle1{PSGMEntCurveCircle::create(context, circ1)};
    PSGMEntCurve* circle2{PSGMEntCurveCircle::create(context, circ2)};

    cylinder.faces.emplace_back(cylinder.shells[0]->createFace(entPlane1));
    cylinder.faces.emplace_back(cylinder.shells[0]->createFace(entPlane2));
    cylinder.faces.emplace_back(cylinder.shells[0]->createFace(entCylinder));

    cylinder.edges.emplace_back(cylinder.body->createEdge(circle1));
    cylinder.edges.emplace_back(cylinder.body->createEdge(circle2));

    auto v=GlobalVariable::getImplicitMatrixVector();
    auto Q=GlobalVariable::getQMatrix();
    auto QAABBSPtr=GlobalVariable::getQAABBSPtr();
    auto qIndex=Q->size();
    auto firstSecondToBeginEnd=GlobalVariable::getFirstSecondToBeginEnd();
    auto qMatrixIndexToPSGMsurfaceSPtrSPtr=GlobalVariable::getQMatrixIndexToPSGMsurfaceSPtrSPtr();

     // add map
    Eigen::Matrix<float, 10, 1> tempFunction;
    PSGMDirection direction;

    auto createPlaneVector = [&](PSGMLCS lcs, Eigen::Matrix<float, 10, 1> &tempFunction)
    {
        PSGMPoint point = lcs.getOrigin();
        PSGMDirection direction = lcs.getZDir();
        double dx = direction.x();
        double dy = direction.y();
        double dz = direction.z();
        double x = point.x();
        double y = point.y();
        double z = point.z();
        tempFunction << 0, 0, 0, dx, 0, 0, dy, 0, dz, -(x * dx + y * dy + z * dz);
    };

    auto createCylinderVector = [&](PSGMLCS lcs, Eigen::Matrix<float, 10, 1> &tempFunction)
    {
        PSGMPoint point = lcs.getOrigin();
        PSGMDirection direction = lcs.getZDir();
        double a = direction.x();
        double b = direction.y();
        double c = direction.z();
        double x0 = point.x();
        double y0 = point.y();
        double z0 = point.z();
        double aa=a*a;
        double bb=b*b;
        double cc=c*c;
        double ab=a*b;
        double ac=a*c;
        double bc=b*c;
        tempFunction<<1-a*a,-a*b,-a*c,-x0*(1-a*a)+a*b*y0+a*c*z0,1-b*b,-b*c,-y0*(1-b*b)+a*b*x0+b*c*z0,1-c*c,-z0*(1-c*c)+a*c*x0+b*c*y0,x0*x0*(1-a*a)+y0*y0*(1-b*b)+z0*z0*(1-c*c)-2*x0*y0;    
    };

    createPlaneVector(lcs1, tempFunction);
    facesPointer->emplace_back(cylinder.faces[0]);
    functionsPointer->emplace_back(tempFunction);
    v->push_back(tempFunction);
    (*faceQMatrixIndex)[cylinder.faces[0]]=Q->size();
    Q->push_back(vcf2gmp2bigintmatrix(tempFunction));
    qMatrixIndexToPSGMsurfaceSPtrSPtr->push_back(cylinder.faces[0]->getSurface());

    createPlaneVector(lcs2, tempFunction);
    (*facesPointer).emplace_back(cylinder.faces[1]);
    (*functionsPointer).emplace_back(tempFunction);
    v->push_back(tempFunction);
    (*faceQMatrixIndex)[cylinder.faces[1]]=Q->size();
    Q->push_back(vcf2gmp2bigintmatrix(tempFunction));
    qMatrixIndexToPSGMsurfaceSPtrSPtr->push_back(cylinder.faces[1]->getSurface());

    createCylinderVector(lcs,tempFunction);
    (*facesPointer).emplace_back(cylinder.faces[2]);
    (*functionsPointer).emplace_back(tempFunction);
    v->push_back(tempFunction);
    (*faceQMatrixIndex)[cylinder.faces[2]]=Q->size();
    Q->push_back(vcf2gmp2bigintmatrix(tempFunction));
    qMatrixIndexToPSGMsurfaceSPtrSPtr->push_back(cylinder.faces[2]->getSurface());

    //surface aabb
    if(lcs.getZDir().isParallel({0,0,1}))
    {
        auto bottomOrgin=lcs1.getOrigin();
        QAABBSPtr->emplace_back(-INFINITY,-INFINITY,bottomOrgin[2],INFINITY,INFINITY,bottomOrgin[2],1e-8);

        auto upOrgin=lcs2.getOrigin();
        QAABBSPtr->emplace_back(-INFINITY,-INFINITY,upOrgin[2],INFINITY,INFINITY,upOrgin[2],1e-8);

        auto cylinderOrgin=lcs.getOrigin();
        QAABBSPtr->emplace_back(cylinderOrgin[0]-radius,cylinderOrgin[1]-radius,-INFINITY,cylinderOrgin[0]+radius,cylinderOrgin[1]+radius,INFINITY);
    }
    else
    {
        QAABBSPtr->emplace_back(-INFINITY,-INFINITY,-INFINITY,INFINITY,INFINITY,INFINITY,1e-8);
        QAABBSPtr->emplace_back(-INFINITY,-INFINITY,-INFINITY,INFINITY,INFINITY,INFINITY,1e-8);
        QAABBSPtr->emplace_back(-INFINITY,-INFINITY,-INFINITY,INFINITY,INFINITY,INFINITY);
    }

    //!TODO:body AABB

    (*edgeFuctionIndex)[cylinder.edges[0]]=GlobalVariable::getAllCurveEquation()->size();
    intersectionDouble(*Q,*firstSecondToBeginEnd,qIndex,qIndex+2,*GlobalVariable::getGinacParamsSPtr(),*GlobalVariable::getAllCurveEquation(),*GlobalVariable::getAllVertexEquation());

    (*edgeFuctionIndex)[cylinder.edges[1]]=GlobalVariable::getAllCurveEquation()->size();
    intersectionDouble(*Q,*firstSecondToBeginEnd,qIndex+1,qIndex+2,*GlobalVariable::getGinacParamsSPtr(),*GlobalVariable::getAllCurveEquation(),*GlobalVariable::getAllVertexEquation());

    cylinder.loops.emplace_back(cylinder.faces[0]->createLoop(PSGMLoopType::Outer));
    cylinder.loops.emplace_back(cylinder.faces[1]->createLoop(PSGMLoopType::Outer));
    cylinder.loops.emplace_back(cylinder.faces[2]->createLoop(PSGMLoopType::Winding));
    cylinder.loops.emplace_back(cylinder.faces[2]->createLoop(PSGMLoopType::Winding));

    cylinder.loops[0]->appendCoedge(cylinder.edges[0]->createCoedge(true));
    cylinder.loops[1]->appendCoedge(cylinder.edges[1]->createCoedge(true));

    cylinder.loops[2]->appendCoedge(cylinder.edges[0]->createCoedge(false));
    cylinder.loops[3]->appendCoedge(cylinder.edges[1]->createCoedge(false));
}

cylinder_hybrid::cylinder_hybrid(PSGMISessionContext *context, PSGMLCS lcs,
                                 double hz, double r)
{
    cylinder_primitive *newCylinder = new cylinder_primitive(conerX, conerY, conerZ, heighZ, radius);
    // tree = new CSGtree(newCylinder);
    tree=(primitive*)newCylinder;
    // facesPointer = new std::vector<PSGMTopoFace *>;
    // functionsPointer = new std::vector<Eigen::Matrix<float, 10, 1>>;
    makeCylinder(context, lcs, heighZ, radius, shape);
}


inline void sphere_hybrid::makeSphere(PSGMISessionContext *context, double cx, double cy, double cz,
                                      double radius, PSGMShape &sphere)
{
    sphere.body = PSGMTopoBody::create(context, PSGMBodyType::Solid);
    PSGMTopoLump *lump{sphere.body->createLump(PSGMLumpType::Solid)};
    sphere.lumps.push_back(lump);

    PSGMTopoShell *shell{lump->createShell()};
    sphere.shells.push_back(shell);

    const PSGMLCS lcs(PSGMPoint{cx, cy, cz}, PSGMDirection{1.0, 0, 0}, PSGMDirection{0, 1.0, 0});
    PSGMSurfSphereSPtr surfSphere{nullptr};
    PSGMSurfSphere::create(lcs, radius, &surfSphere);
    PSGMEntSurface *entSphere{PSGMEntSurfSphere::create(context, surfSphere)};
    sphere.faces.emplace_back(sphere.shells[0]->createFace(entSphere));

    // add map
    Eigen::Matrix<float, 10, 1> tempFunction;

    PSGMPoint point = lcs.getOrigin();
    PSGMDirection direction = lcs.getZDir();
    double dx = direction.x();
    double dy = direction.y();
    double dz = direction.z();
    double x = point.x();
    double y = point.y();
    double z = point.z();
    tempFunction << 1, 0, 0, -2 * x, 1, 0, -2 * y, 1, -2 * z, x * x + y * y + z * z - radius * radius;

    (*facesPointer).emplace_back(sphere.faces[0]);
    (*functionsPointer).emplace_back(tempFunction);
    GlobalVariable::getImplicitMatrixVector()->push_back(tempFunction);
    (*faceQMatrixIndex)[sphere.faces[0]]=GlobalVariable::getQMatrix()->size();
    GlobalVariable::getQMatrix()->push_back(vcf2gmp2bigintmatrix(tempFunction));
    // (*GlobalVariable::getQMatrixIndexToPSGMsurfaceSPtrSPtr())[GlobalVariable::getQMatrix()->size()]=sphere.faces[0];
    GlobalVariable::getQMatrixIndexToPSGMsurfaceSPtrSPtr()->push_back(sphere.faces[0]->getSurface());
    GlobalVariable::getQAABBSPtr()->emplace_back(x-radius,y-radius,z-radius,x+radius,y+radius,z+radius);
}

sphere_hybrid::sphere_hybrid(PSGMISessionContext *context, double cx, double cy, double cz,
                             double r)
    : conerX(cx), conerY(cy), conerZ(cz), radius(r)
{
    sphere_primitive *newSphere = new sphere_primitive(conerX, conerY, conerZ, radius);
    // tree = new CSGtree(newSphere);
    // tree=dynamic_pointer_cast<primitive>(newSphere);
    tree=(primitive*)newSphere;
    facesPointer = new std::vector<PSGMTopoFace *>;
    functionsPointer = new std::vector<Eigen::Matrix<float, 10, 1>>;
    makeSphere(context, conerX, conerY, conerZ, radius, shape);
}