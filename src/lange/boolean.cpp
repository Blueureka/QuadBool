#include <lange/boolean.h>


void interpolateBSpline(curve_equation& curve,double step,double startParameter,double endParameter,hybridPresentation& result,PSGMISessionContext* context)
{
    psgmstd::vector<PSGMPoint> points;
    {
        auto p=curve.getPoint(tan(startParameter / 2),1,1);
        points.emplace_back(p.x,p.y,p.z);
        PSGMPoint psgmPoint{p.x,p.y,p.z};
        result.shape.vertices.emplace_back(result.shape.body->createVertex(psgmPoint));
    }
    for (double i = startParameter+step; i < endParameter; i += step)
    {
        double j = tan(i / 2);
        auto p=curve.getPoint(j,1,1);
        points.emplace_back(p.x,p.y,p.z);
    }
    {
        auto p=curve.getPoint(tan(endParameter / 2),1,1);
        points.emplace_back(p.x,p.y,p.z);
        PSGMPoint psgmPoint{p.x,p.y,p.z};
        result.shape.vertices.emplace_back(result.shape.body->createVertex(psgmPoint));
    }
    PSGMCurveBSplineSPtr bspSPtr
    {
        PSGMCurveBSpline::createByPoints(points, curve.isLoop, 3)
    };
    PSGMEntCurve* bsp1{PSGMEntCurveBSpline::create(context,bspSPtr)};
    result.shape.edges.emplace_back(result.shape.body->createEdge(result.shape.vertices[result.shape.vertices.size()-2], result.shape.vertices[result.shape.vertices.size()-1], bsp1));
}


void BooleanOperation::handleRingEdge(curve_equation& curve,CSGtree* tree)
{
    auto infPoint=curve.negativeInf;
    if(!tree->isPointOn(infPoint))
    {
        return;
    }
    psgmstd::vector<PSGMPoint> points;
    points.emplace_back(infPoint.x,infPoint.y,infPoint.z);
    double step= (2*M_PI)/curve.visualizaitonPoint_cnt;
    for (double i = -M_PI+step; i < M_PI; i += step)
    {
        double j = tan(i / 2);
        auto p=curve.getPoint(j,1,1);
        points.emplace_back(p.x,p.y,p.z);
    }
    PSGMCurveBSplineSPtr bspSPtr
    {
        PSGMCurveBSpline::createByPoints(points, true, 3)
    };
    PSGMEntCurve* bsp1{PSGMEntCurveBSpline::create(context,bspSPtr)};
    result->shape.edges.emplace_back(result->shape.body->createEdge(bsp1));
    
    auto qMatrixIndexToPSGMsurfaceSPtrSPtr=GlobalVariable::getQMatrixIndexToPSGMsurfaceSPtrSPtr();
    
    auto surface1SPtr=(*qMatrixIndexToPSGMsurfaceSPtrSPtr)[curve.first_index];
    auto surface2SPtr=(*qMatrixIndexToPSGMsurfaceSPtrSPtr)[curve.second_index];
    
    
    auto p0=curve.getPoint(0,1,1);
    auto p1=curve.getPoint(1,1,1);
    
    PSGMPoint pp0{p0.x,p0.y,p0.z};
    PSGMPoint pp1{p1.x,p1.y,p1.z};
    auto vec=pp1-pp0;

    auto norm1=surface1SPtr->getNormal(surface1SPtr->getParamAt(pp0));
    auto norm2=surface2SPtr->getNormal(surface2SPtr->getParamAt(pp0));

    bool sense1 = norm1.cross(norm2).dot(vec) > 0;
    bool sense2 = !sense1;
    
    //UNION
    unsigned int tempId1=(*QMatrixIndexToTempIndexSPtr)[curve.first_index];
    unsigned int tempId2=(*QMatrixIndexToTempIndexSPtr)[curve.second_index];
    switch (operation)
    {
    case UNION:
    {
        (*faceList)[tempId2].pool->same->emplace_back(result->shape.edges.back()->createCoedge(!sense2));
        (*faceList)[tempId1].pool->same->emplace_back(result->shape.edges.back()->createCoedge(!sense1));
        break;
    }
    case INTERSECTION:
    {
        (*faceList)[tempId1].pool->same->emplace_back(result->shape.edges.back()->createCoedge(sense1));
        (*faceList)[tempId2].pool->same->emplace_back(result->shape.edges.back()->createCoedge(sense2));
        break;
    }
    case SUBTRACTION:
    {
        (*faceList)[tempId1].pool->same->emplace_back(result->shape.edges.back()->createCoedge(sense1));
        (*faceList)[tempId2].pool->same->emplace_back(result->shape.edges.back()->createCoedge(sense2));
        break;
    }
    default:
        break;
    }    
}

void BooleanOperation::divideCurveAndMakeEdge(curve_equation& intersectCurve,CSGtree* tree)
{
    std::cout<<"divideCurveAndMakeEdge"<<std::endl;

    std::vector<double> parameters;
    bool addClose=false;
    auto& seq=intersectCurve.pointList;
    if(seq.empty())
    {
        std::cout<<"seq is empty"<<std::endl;
        if(intersectCurve.isAllInterval)
        {   
            handleRingEdge(intersectCurve,tree);
        }
        return;
    }
    
    edge_equation edge;
    // edge.curveSPtr.reset(&intersectCurve);

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
            edge.validSegment.emplace_back(parameters[i],parameters[i+1]);
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
                edge.validSegment.emplace_back(parameters[i],parameters[i]);
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
                edge.validSegment.emplace_back(INF_PARAMETER,INF_PARAMETER);
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
                edge.validSegment.emplace_back(parameters[0],parameters[0]);
            }
        }
        if (!endPushedOverlap)
        {
            // i == numOfPtsOnCurve - 1 is surface curve overlap end point
            bool isValid=tree->isPointOn(intersectCurve.getPoint(parameters.back(),1,1));
            if (isValid)
            {
                intersectCurve.validSegment.emplace_back(parameters.back(),parameters.back());
                edge.validSegment.emplace_back(parameters.back(),parameters.back());
            }
        }
    }

    for(auto segment:edge.validSegment)
    {
        std::cout<<"startParameter= "<<segment.startParameter<<" endParameter= "<<segment.endParameter<<endl;
    }

    std::cout<<"divideCurve finish"<<std::endl;
    //interPoints list to segment
    intersectCurve.pointList.clear();


    //according segment save result to pool
    // for(auto segment:intersectCurve.validSegment)
    // {
    //     unsigned int face1id;
    //     unsigned int face2id;
    //     //face1id= intersectCurve.first_index (surface id ) to face id 
    //     //face2id= intersectCurve.first_index (surface id ) to face id
    //     bool someCondition;
    //     if(someCondition)
    //     {
    //         PSGMTopoCoedge* coedge;
    //         faceEdgePools[face1id].same.push_back(coedge);
    //     }
    // }
    for(auto segment:edge.validSegment)
    {
        psgmstd::vector<PSGMPoint> points;
        auto startPoint=intersectCurve.getPoint(segment.startParameter,1,1);
        auto endPoint=intersectCurve.getPoint(segment.endParameter,1,1);
        result->shape.vertices.emplace_back(result->shape.body->createVertex({startPoint.x,startPoint.y,startPoint.z}));
        result->shape.vertices.emplace_back(result->shape.body->createVertex({endPoint.x,endPoint.y,endPoint.z}));

        points.emplace_back(startPoint.x,startPoint.y,startPoint.z);
        double step= (segment.endParameter-segment.startParameter)/intersectCurve.visualizaitonPoint_cnt;
        for (double i = segment.startParameter+step; i < segment.endParameter; i += step)
        {
            auto p=intersectCurve.getPoint(i,1,1);
            points.emplace_back(p.x,p.y,p.z);
        }
        points.emplace_back(endPoint.x,endPoint.y,endPoint.z);
        

        // for(auto p:points)
        // {
        //     cout<<"{"<<p[0]<<","<<p[1]<<","<<p[2]<<"},";
        // }cout<<endl;

        PSGMCurveBSplineSPtr bspSPtr
        {
            PSGMCurveBSpline::createByPoints(points, false, 3)
        };
        PSGMEntCurve* bsp1{PSGMEntCurveBSpline::create(context,bspSPtr)};
        auto verticesSize=result->shape.vertices.size();
        result->shape.edges.emplace_back(result->shape.body->createEdge(result->shape.vertices[verticesSize-2],result->shape.vertices[verticesSize-1],bsp1));

        auto qMatrixIndexToPSGMsurfaceSPtrSPtr=GlobalVariable::getQMatrixIndexToPSGMsurfaceSPtrSPtr();
        auto surface1SPtr=(*qMatrixIndexToPSGMsurfaceSPtrSPtr)[intersectCurve.first_index];
        auto surface2SPtr=(*qMatrixIndexToPSGMsurfaceSPtrSPtr)[intersectCurve.second_index];
        
        
        // auto p0=intersectCurve.getPoint(segment.startParameter,1,1);
        // auto p1=intersectCurve.getPoint((segment.startParameter+segment.endParameter)/2,1,1);
        
        PSGMPoint pp0=bspSPtr->getPoint(0);
        PSGMPoint pp1=bspSPtr->getPoint(0.01);
        auto vec=pp1-pp0;

        auto norm1=surface1SPtr->getNormal(surface1SPtr->getParamAt(pp0));
        auto norm2=surface2SPtr->getNormal(surface2SPtr->getParamAt(pp0));

        bool sense1 = norm1.cross(norm2).dot(vec) > 0;
        bool sense2 = !sense1;

        unsigned int tempId1=(*QMatrixIndexToTempIndexSPtr)[intersectCurve.first_index];
        unsigned int tempId2=(*QMatrixIndexToTempIndexSPtr)[intersectCurve.second_index];
        switch (operation)
        {
        case UNION:
        {
            (*faceList)[tempId2].pool->same->emplace_back(result->shape.edges.back()->createCoedge(!sense2));
            (*faceList)[tempId1].pool->same->emplace_back(result->shape.edges.back()->createCoedge(!sense1));
            break;
        }
        case INTERSECTION:
        {
            (*faceList)[tempId1].pool->same->emplace_back(result->shape.edges.back()->createCoedge(sense1));
            (*faceList)[tempId2].pool->same->emplace_back(result->shape.edges.back()->createCoedge(sense2));
            break;
        }
        case SUBTRACTION:
        {
            (*faceList)[tempId1].pool->same->emplace_back(result->shape.edges.back()->createCoedge(!sense1));
            (*faceList)[tempId2].pool->same->emplace_back(result->shape.edges.back()->createCoedge(!sense1));
            break;
        }
        default:
            break;
        }    

    }

}

int BooleanOperation::findFirst(int faceIndex,bool same)
{
    auto poolSPtr=same?((*faceList)[faceIndex].pool->same):((*faceList)[faceIndex].pool->diff);
    for(int i=0;i<poolSPtr->size();i++)
    {
        if(!(*poolSPtr)[i].used)return i;
    }
    return -1;
}

int BooleanOperation::findNext(int faceIndex,int edgeIndex,bool same)
{
    int ans=0;
    auto poolSPtr=same?((*faceList)[faceIndex].pool->same):((*faceList)[faceIndex].pool->diff);

    auto curEdge=(*poolSPtr)[edgeIndex].coedge;
    auto endVertex=curEdge->getEndVertex();
    if(endVertex)
    {
        auto endPoint=endVertex->getEntPoint()->getPoint();
        for(;ans<poolSPtr->size();ans++)
        {
            if((*poolSPtr)[ans].used)continue;
            auto edge=(*poolSPtr)[ans].coedge;
            auto stratVertex=edge->getStartVertex();
            if(stratVertex)
            {
                auto startPoint=stratVertex->getEntPoint()->getPoint();
                if(startPoint.isEqual(endPoint))return ans;
            }
            else
            {
                continue;
            }
        }
    }
    else//nullptr
    {
        return -1;
    }
    return -1;
}

void BooleanOperation::allBoolean(CSGtree finalTree)
{
    std::cout<<"BooleanOperation::allBoolean : begin"<<std::endl;
    
    // if(operation==CSGoperation::INTERSECTION)
    // {
    //     if(bodyA->box.contains(bodyB->box))
    //     {
    //         result.reset(bodyB);
    //         return;
    //     }
    //     else if(bodyB->box.contains(bodyA->box))
    //     {
    //         result.reset(bodyA);
    //         return;
    //     }
    //     else if(!bodyA->box.intersects(bodyB->box))
    //     {
    //         return;
    //     }
    // }
    // else if(operation==CSGoperation::UNION)
    // {
    //     if(bodyA->box.contains(bodyB->box))
    //     {
    //         result.reset(bodyA);
    //         return;
    //     }
    //     else if(bodyB->box.contains(bodyA->box))
    //     {
    //         result.reset(bodyB);
    //         return;
    //     }
    //     else if(!bodyA->box.intersects(bodyB->box))
    //     {
    //         //!!TODO:
    //         //result=a+b
    //         return;
    //     }
    // }
    // else if(operation==CSGoperation::SUBTRACTION)
    // {
    //     if(bodyB->box.contains(bodyA->box))
    //     {
    //         //!TODO:
    //         //result = null
    //         return;
    //     }
    //     else if(!bodyA->box.intersects(bodyB->box))
    //     {
    //         result.reset(bodyA);
    //         return;
    //     }
    // }

    //init
    //may be class member

    //find surface
    std::vector<unsigned int> tempQMatrixIndex;
    //find face
    std::vector<PSGMTopoFace *> tempQMatrixIndexToFaceSPtr;
    //find curve
    std::vector<unsigned int> tempCurveIndex;
    
    //!!TODO:!!
    //may change
    for(auto faceInfo:(*faceList))
    {
        tempQMatrixIndex.push_back(faceInfo.qMatrixSurfaceId);
    }

    //!!TODO:!!
    //intersectionDouble some input parameters become SPtr
    //or some input ignore, 在函数里通过GlobalVariable获得
    //

    //global Variable
    auto Q=(*GlobalVariable::getQMatrix());
    auto firstSecondToBeginEnd=(*GlobalVariable::getFirstSecondToBeginEnd());
    auto ginacParams=(*GlobalVariable::getGinacParamsSPtr());
    auto all_curve_equationSPtr=GlobalVariable::getAllCurveEquation();
    auto all_vertex_equation=(*GlobalVariable::getAllVertexEquation());

    quad_inter<bigint> qi_ic ;
    QIOutputter outputter;
    outputter.setOmitImaginaryParametrizations();
    QIOutputInfo *aasdadadadas;

    PRINT_VAR(tempQMatrixIndex.size());
    //find all intersection
    for(int i=0;i<tempQMatrixIndex.size();i++)
    {
        unsigned int first_index=tempQMatrixIndex[i];
        PRINT_VAR(first_index);
        for(int j=i+1;j<tempQMatrixIndex.size();j++)
        {
            unsigned int second_index=tempQMatrixIndex[j];
            PRINT_VAR(second_index);

            // pair<unsigned int,unsigned int> beginEnd=intersectionDouble(Q,firstSecondToBeginEnd,first_index,second_index,ginacParams,all_curve_equation,all_vertex_equation);
            pair<unsigned int,unsigned int> beginEnd=intersectionDouble(Q,firstSecondToBeginEnd,first_index,second_index);
            std::cout<<"{begin,end}={ "<<beginEnd.first<<", "<<beginEnd.second<<"}"<<std::endl;

            if(beginEnd.first==beginEnd.second)
            {
                
                continue;
            }
            for(int begin=beginEnd.first;begin<beginEnd.second;begin++){tempCurveIndex.push_back(begin);}
            //TSQI 
            //get three QuadSurface intersection point 
            for(int k=j+1;k<tempQMatrixIndex.size();k++)
            {
                unsigned int third_index=tempQMatrixIndex[k];
                PRINT_VAR(third_index);

                TQSI::seq seqST= intersectionTriple(beginEnd,Q[third_index]);
                // for(int begin=beginEnd.first;begin<beginEnd.second;begin++)
                // {
                    
                //     TQSI::seq seqST= intersectionTriple(all_curve_equation[begin],Q[third_index]);
                //     curveAddSeqToPointList(all_curve_equation[begin],seqST);
                // }

                // divideCurve(all_curve_equation[beginEnd.first],&finalTree,seqST);

                // TQSI::seq seqSU= intersectionTriple(Q[first_index], Q[third_index],Q[second_index]);

                // matchSeq(seqST,seqSU);

                pair<unsigned int,unsigned int> firstThird={first_index,third_index};
                if(firstSecondToBeginEnd.find(firstThird)!=firstSecondToBeginEnd.end())
                {
                    beginEnd=firstSecondToBeginEnd[firstThird];
                }
                else
                {
                    // beginEnd.first=all_curve_equation.size();
                    
                    // qi_ic = intersection(Q[first_index], Q[third_index], true, std::cout);
                    // outputter.output(qi_ic, Q[first_index], Q[third_index]);
                    // aasdadadadas = outputter.getOutput();
                    // handleCurveComponent(aasdadadadas,ginacParams,all_curve_equation,all_vertex_equation,first_index,third_index);

                    // beginEnd.second=all_curve_equation.size();
                    // firstSecondToBeginEnd[firstThird]=beginEnd;

                    pair<unsigned int,unsigned int> beginEnd=intersectionDouble(Q,firstSecondToBeginEnd,first_index,third_index);
                    if(beginEnd.first==beginEnd.second)
                    {
                        continue;
                    }
                }
                std::cout<<"{begin,end}={ "<<beginEnd.first<<", "<<beginEnd.second<<"}"<<std::endl;

                TQSI::seq seqSU= intersectionTriple(beginEnd,Q[second_index]);
                for(int begin=beginEnd.first;begin<beginEnd.second;begin++)
                {
                    tempCurveIndex.push_back(begin);
                    // TQSI::seq seqSU= intersectionTriple(all_curve_equation[begin],Q[second_index]);
                    // curveAddSeqToPointList(all_curve_equation[begin],seqSU);
                }
                // divideCurve(all_curve_equation[beginEnd.first],&finalTree,seqSU);

                // TQSI::seq seqTU= intersectionTriple(Q[second_index], Q[third_index],Q[first_index]);

                // matchSeq(seqST,seqTU);

                pair<unsigned int,unsigned int> secondThird={second_index,third_index};
                if(firstSecondToBeginEnd.find(secondThird)!=firstSecondToBeginEnd.end())
                {
                    beginEnd=firstSecondToBeginEnd[secondThird];
                }
                else
                {
                    // beginEnd.first=all_curve_equation.size();
                    
                    // cout<<"qq1= "<<Q[second_index]<<endl;
                    // cout<<"qq2= "<<Q[third_index]<<endl;

                    // qi_ic = intersection(Q[second_index], Q[third_index], true, std::cout);
                    // outputter.output(qi_ic, Q[second_index], Q[third_index]);
                    // aasdadadadas = outputter.getOutput();
                    // handleCurveComponent(aasdadadadas,ginacParams,all_curve_equation,all_vertex_equation,second_index,third_index);

                    // beginEnd.second=all_curve_equation.size();
                    // firstSecondToBeginEnd[secondThird]=beginEnd;

                    pair<unsigned int,unsigned int> beginEnd=intersectionDouble(Q,firstSecondToBeginEnd,second_index,third_index);
                    if(beginEnd.first==beginEnd.second)
                    {
                        continue;
                    }
                }
                std::cout<<"{begin,end}={ "<<beginEnd.first<<", "<<beginEnd.second<<"}"<<std::endl;
                
                TQSI::seq seqTU= intersectionTriple(beginEnd,Q[first_index]);
                for(int begin=beginEnd.first;begin<beginEnd.second;begin++)
                {
                    tempCurveIndex.push_back(begin);
                    // TQSI::seq seqTU= intersectionTriple(all_curve_equation[begin],Q[first_index]);
                    // curveAddSeqToPointList(all_curve_equation[begin],seqTU);
                }
            }
        }
    }
    std::cout<<"tempCurveIndex.size()= "<<tempCurveIndex.size()<<std::endl;

    //after interction, obtain all curve with point_list
    //divide curve and result save in curve.segment. 
    //according segment to make edge,add to pool
    std::cout<<"all_curve_equation.size = "<<all_curve_equationSPtr->size()<<endl;
    std::sort(tempCurveIndex.begin(),tempCurveIndex.end());
    
    for(auto aa:tempCurveIndex)
    {
        cout<<aa<<endl;
    }
    
    tempCurveIndex.erase(std::unique(tempCurveIndex.begin(),tempCurveIndex.end()),tempCurveIndex.end());
    std::cout<<"tempCurveIndex.size()= "<<tempCurveIndex.size()<<std::endl;

    for(int i=0;i<tempCurveIndex.size();i++)
    {
        auto& curve=(*all_curve_equationSPtr)[i];
        for(const auto& segment:curve.validSegment)
        {
            if(segment.startParameter==segment.endParameter)
            {
                curve.pointList.emplace_back(curve.getPoint(segment.startParameter,1,1),segment.startParameter);
            }
            else
            {
                curve.pointList.emplace_back(curve.getPoint(segment.startParameter,1,1),segment.startParameter);
                curve.pointList.emplace_back(curve.getPoint(segment.endParameter,1,1),segment.endParameter);
            }
        }
        std::sort(curve.pointList.begin(),curve.pointList.end(),[](TQSI::intPoint a,TQSI::intPoint b)
        {
            return a.parameter<b.parameter;
        });
        curve.pointList.erase(std::unique(curve.pointList.begin(), curve.pointList.end(), [](TQSI::intPoint a, TQSI::intPoint b)
                                { return a.parameter == b.parameter; }),
                    curve.pointList.end());
        
        // divideCurve(curve,&finalTree);
        divideCurveAndMakeEdge(curve,&finalTree);
    }


    for(int i=0;i<faceList->size();i++)
    {
        auto faceInfo=(*faceList)[i];
        if(faceInfo.pool->diff->size()||faceInfo.pool->same->size())
        {
            
            result->shape.faces.emplace_back(result->shape.shells[0]->createFace(faceInfo.face->getEntSurface()));
            auto face=result->shape.faces.back();
            
            for(int sd = 0; sd < 2; sd++)
            {
                bool same= sd==0? true:false;
                for(;;)
                {
                    int first=findFirst(i,same);
                    if(first==-1)break;

                    auto poolSPtr=same?((*faceList)[i].pool->same):((*faceList)[i].pool->diff);
                    result->shape.loops.emplace_back(face->createLoop(PSGMLoopType::Unclear));
                    auto loop=result->shape.loops.back();
                    loop->appendCoedge((*poolSPtr)[first].coedge);
                    (*poolSPtr)[first].used=true;

                    int he_init = first, he_prev = first;

                    bool bad_loop = false;
                    for(;;)
                    {
                        int nextIndex=findNext(i,first,same);
                        if(nextIndex==-1)
                        {
                            bad_loop=true;
                            break;
                        }

                        loop->appendCoedge((*poolSPtr)[nextIndex].coedge);
                        (*poolSPtr)[nextIndex].used=true;

                    }
                }
            }
        }
    }

    

}



