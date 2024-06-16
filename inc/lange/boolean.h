#pragma once

#include <lange/Hybrid.h>
#include <CSG/csg.h>
#include <libqi/qi.h>
#include <lange/TQSI.h>

struct HEPair
{
    int vno;
    PSGMTopoCoedge* coedge;
    bool used;
    HEPair(PSGMTopoCoedge* _coedge):vno(-1),coedge(_coedge),used(false){}
};

struct FaceEdgePool
{
    //temp id, use index directly
    // unsigned int faceId;
    //may change
    std::shared_ptr<std::vector<HEPair>> same=std::make_shared<std::vector<HEPair>>();
    std::shared_ptr<std::vector<HEPair>> diff=std::make_shared<std::vector<HEPair>>();
};

struct faceInfo
{
    // unsigned int tempId;

    PSGMTopoFace* face;
    unsigned int qMatrixSurfaceId;

    std::shared_ptr<FaceEdgePool> pool=std::make_shared<FaceEdgePool>();
    unsigned int owner;
    faceInfo(PSGMTopoFace* _face,unsigned int _qMatrixSurfaceId):face(_face),qMatrixSurfaceId(_qMatrixSurfaceId)
    {
        // pool=new FaceEdgePool();
    }
};

class BooleanOperation
{
public:
    BooleanOperation(hybridPresentation* _bodyA,hybridPresentation* _bodyB,CSGoperation _operation,PSGMISessionContext *_context):bodyA(_bodyA),bodyB(_bodyB),operation(_operation),context(_context)
    {
        //init
        for(const auto& face:bodyA->shape.faces)
        {
            unsigned int qMatrixIndex=(*bodyA->faceQMatrixIndex)[face];
            // (*QMatrixIndexToFacePointerSPtr)[qMatrixIndex]=face;
            // (*faceToQMatrixIndexSPtr)[face]=qMatrixIndex;

            // (*faceToTempIndexSPtr)[face]=faceList->size();
            (*QMatrixIndexToTempIndexSPtr)[qMatrixIndex]=faceList->size();
            faceList->emplace_back(face,qMatrixIndex);
        }

        for(const auto& face:bodyB->shape.faces)
        {
            unsigned int qMatrixIndex=(*bodyB->faceQMatrixIndex)[face];
            // (*QMatrixIndexToFacePointerSPtr)[qMatrixIndex]=face;
            // (*faceToQMatrixIndexSPtr)[face]=qMatrixIndex;

            // (*faceToTempIndexSPtr)[face]=faceList->size();
            (*QMatrixIndexToTempIndexSPtr)[qMatrixIndex]=faceList->size();
            faceList->emplace_back(face,qMatrixIndex);
        }

        // for(const auto& edge:bodyA->shape.edges)
        // {
        //     unsigned int curveEquationIndex=(*bodyA->edgeFuctionIndex)[edge];
        //     (*CurveIndexToEdgePointerSPtr)[curveEquationIndex]=edge;
        //     (*edgeToCurveIndexSPtr)[edge]=curveEquationIndex;
        // }
        // for(const auto& edge:bodyB->shape.edges)
        // {
        //     unsigned int curveEquationIndex=(*bodyB->edgeFuctionIndex)[edge];
        //     (*CurveIndexToEdgePointerSPtr)[curveEquationIndex]=edge;
        //     (*edgeToCurveIndexSPtr)[edge]=curveEquationIndex;
        // }

        result->shape.body = PSGMTopoBody::create(context, PSGMBodyType::Solid);
        PSGMTopoLump *lump{result->shape.body->createLump(PSGMLumpType::Solid)};
        result->shape.lumps.push_back(lump);

        PSGMTopoShell *shell{lump->createShell()};
        result->shape.shells.push_back(shell);
    }

    void handleRingEdge(curve_equation& curve,CSGtree* tree);
    
    void divideCurveAndMakeEdge(curve_equation& intersectCurve,CSGtree* tree);

    void allBoolean(CSGtree finalTree);

    int findFirst(int faceIndex,bool same);
    int findNext(int faceIndex,int edgeIndex,bool same);

    std::shared_ptr<std::vector<faceInfo>> getFaceList()
    {
        return this->faceList;
    }

private:
    hybridPresentation* bodyA;
    hybridPresentation* bodyB;
    CSGoperation operation;
    PSGMISessionContext *context;

    std::shared_ptr<hybridPresentation> result=std::make_shared<hybridPresentation>();

    std::shared_ptr<std::vector<faceInfo>> faceList=std::make_shared<std::vector<faceInfo>>();

    std::shared_ptr<std::map<PSGMTopoFace*, unsigned int>> faceToTempIndexSPtr;
    std::shared_ptr<std::map<unsigned int,PSGMTopoFace*>> QMatrixIndexToFacePointerSPtr;
    std::shared_ptr<std::map<unsigned int,unsigned int>> QMatrixIndexToTempIndexSPtr=std::make_shared<std::map<unsigned int,unsigned int>>();
    std::shared_ptr<std::map<PSGMTopoFace*, unsigned int>> faceToQMatrixIndexSPtr;
    
    std::shared_ptr<std::map<unsigned int,PSGMTopoEdge*>> CurveIndexToEdgePointerSPtr;
    std::shared_ptr<std::map<PSGMTopoEdge*,unsigned int>> edgeToCurveIndexSPtr;
    std::vector<unsigned int> tempCurveIndex;
};



void interpolateBSpline(curve_equation& curve,double step,double startParameter,double endParameter,hybridPresentation& result,PSGMISessionContext* context);

void interpolateIntersctionCurve();



