#include <iostream>
#include <string>
#include <igl/opengl/glfw/Viewer.h>
#include <libqi/qi.h>
#include <cmath>
// using namespace GiNaC;

#include <lange/Hybrid.h>
#include <lange/TQSI.h>
#include <lange/Visualization.h>
#include <lange/boolean.h>

#include <TEST/testcase.h>

#include <lange/Global.h>

igl::opengl::glfw::Viewer viewer;
TQSI::Viewer Viewer(&viewer);

std::string Guioutput = "";
std::string Debugoutput = "";
std::string statusstr = "";
bool isTracing = false;
bool generate = false;

class Timer
{
    std::string name;
    std::chrono::_V2::system_clock::time_point start, end;

public:
    Timer(const char *n) : name(n)
    {
        start = std::chrono::system_clock::now();
    }
    Timer(std::string n) : name(n)
    {
        start = std::chrono::system_clock::now();
    }
    ~Timer()
    {
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << name << elapsed.count() * 1000 << "ms\n";
    }
};



#include "visualization/tessellate/PSGMTessCacheManager.h"

#include "transaction/PSGMISession.h"
#include "cache/PSGMBodyCacheManager.h"
struct PSGMContext
{
    PSGMApiDeltaCallbacks m_callbacks;
    PSGMSessionStartOption m_option{1, "apirecords.json", 0};
    PSGMISessionSPtr m_session;
    PSGMISessionContext* m_context{nullptr};
    PSGMAlgPrimitiveModelSPtr m_model;
    PSGMIEntManagerSPtr m_entMgr;
    PSGMIEngine* m_engine{nullptr};
    PSGMTag m_sysColor{0};
    PSGMTag m_sysColor2{0};
    // PSGMSessionStartOption m_option{0, "apirecords.json", 0};
    PSGMTopologySPtr m_topology;
    // PSGMTag m_sysColor{0};
    // PSGMTag m_sysColor2{0};
    PSGMContext()
    {
        // m_topology = std::make_shared<PSGMAlgTopology>(m_context);
        PSGMApiEngineOption option;
        m_engine = startEngine(option);

        PSGMAlgorithm::init();
        PSGMTessCacheManager::init();
        PSGMBodyCacheManager::init();

        m_session = PSGMISession::create();
        m_session->start(m_option);
        m_context = m_session->getContext().get();
        m_model = std::make_shared<PSGMAlgPrimitiveModel>(m_context);
        m_entMgr = m_context->getEntManager();
    }
} PSGMContext;

#define printOnce(X) std::cout<<X<<std::endl;
#define PRINT_VAR(value) std::cout << #value << " = " << (value) << std::endl;




int main()
{
    // GiNaC::symbol u("u");
    // GiNaC::symbol v("v");
    // GiNaC::symbol Delta("Delta");
    // GiNaC::symtab table;
    // table["u"] = u;
    // table["v"] = v;
    // table["Delta"] = Delta;

    GiNaC::symbol u=GlobalVariable::getU();
    GiNaC::symbol v=GlobalVariable::getV();
    GiNaC::symbol Delta=GlobalVariable::getDelta();
    GiNaC::symtab table=(*GlobalVariable::getTableSPtr());
    GiNaC::parser reader(table);
    GinacParams ginacParams=GinacParams(u,Delta,&table);

    auto ginacParamsSPtr=GlobalVariable::getGinacParamsSPtr();
    GlobalVariable::getGinacParamsSPtr().reset(&ginacParams);

    std::vector<curve_equation> all_curve_equation=(*GlobalVariable::getAllCurveEquation());
    std::vector<vertex_equation> all_vertex_equation=(*GlobalVariable::getAllVertexEquation());
    auto qMatrixIndexToPSGMsurfaceSPtrSPtr=GlobalVariable::getQMatrixIndexToPSGMsurfaceSPtrSPtr();



    // {Eigen::Matrix<float, 10, 1> vv1 = Eigen::Matrix<float, 10, 1>::Zero();
    // Eigen::Matrix<float, 10, 1> vv2 = Eigen::Matrix<float, 10, 1>::Zero();
    // vv1<<1,0,0,-2,1,0,0,1,0,0;
    // vv2<<1,0,0,0,1,0,0,0,0,-1; 
    // // vv1<<-1,0,0,0,0,0,2,-1,0;
    // // vv2<<-3,0,0,0,1,0,0,-1,0;
    // // vv1<<1,0,0,0,-1,0,4,1,0;
    // // vv2<<-3,0,0,0,1,0,0,1,0,0;
    // bigint_matrix qq1=vcf2gmp2bigintmatrix(vv1); bigint_matrix qq2=vcf2gmp2bigintmatrix(vv2);
    // cout<<"qq1= "<<qq1<<endl;
    // cout<<"qq2= "<<qq2<<endl;
    // quad_inter<bigint> qi_ic = intersection(qq2, qq1, true, std::cout);
    // QIOutputter outputter;
    // outputter.setOmitImaginaryParametrizations();
    //              outputter.output(qi_ic, qq1,qq2);
                 
                 
    // auto aasdadadadas = outputter.getOutput();
                
    //             handleCurveComponent(aasdadadadas,ginacParams,all_curve_equation,all_vertex_equation,1,2);
    //             Viewer.addCurve(all_curve_equation);
    // Viewer.setVisualizateCurve(true);
    // Viewer.addSurface(vv2);
    // Viewer.addSurface(vv1);
    // Viewer.setVisualizateSurface(true);
    // // Viewer.setVisualizatePoint(true);
    // // Viewer.setVisualizateCSG(true);
    // }
    // std::cout<<"aaaaa"<<std::endl;
    //  Viewer.forDeBug();std::cout<<"bbbbbbb"<<std::endl;

    // return 0;


    cylinder_hybrid cylinder10(PSGMContext.m_context, 0, 0, -2, 4, 1);
    // sphere_hybrid cylinder10(PSGMContext.m_context, -0.5, 0, 0, 1);
    cout<<"make cylinder done"<<endl;
    sphere_hybrid sphere10(PSGMContext.m_context, 1, 0, 0, 1);
    cout<<"make sphere done"<<endl;
    // CSGtree ccc=CSGtree(CSGoperation::SUBTRACTION,cylinder10.tree,sphere10.tree);
    CSGtree ccc=CSGtree(CSGoperation::SUBTRACTION,sphere10.tree,cylinder10.tree);
    BooleanOperation bl(&sphere10,&cylinder10,CSGoperation::SUBTRACTION,PSGMContext.m_context);
    // BooleanOperation bl(&cylinder10,&sphere10,CSGoperation::SUBTRACTION,PSGMContext.m_context);
    bl.allBoolean(ccc);
    std::cout<<"bool done"<<std::endl;

    auto faceList=bl.getFaceList();
    for(auto faceInfo:(*faceList))
    {
        std::cout<<"same.size()"<<faceInfo.pool->same->size()<<std::endl;
        std::cout<<"diff.size()"<<faceInfo.pool->diff->size()<<std::endl;
    }

    // PSGMSessionSaveOption SaveOption;
    // SaveOption.m_transmitFormat = PSGMTransmitFormatType::Step;
    // PSGMContext.m_session->save("abcd123",SaveOption);
    // cout<<"save done"<<endl;
    PS_CLOUD_VIEW_MODEL(PSGMContext.m_session);
    return 0;
    
    Viewer.addCurve((*GlobalVariable::getAllCurveEquation()));
    Viewer.addCSG(ccc);
    Viewer.setVisualizateCurve(true);
    // Viewer.setVisualizatePoint(true);
    Viewer.setVisualizateCSG(true);
    Viewer.visualizate();
    return 0;



    // auto tt=GlobalVariable::getTest();
    // std::cout<<*tt<<std::endl;
    // *tt=4;
    // std::cout<<*GlobalVariable::getTest()<<std::endl;
    // return 0;
    Eigen::Matrix<float, 10, 1> v1 = Eigen::Matrix<float, 10, 1>::Zero();
    Eigen::Matrix<float, 10, 1> v2 = Eigen::Matrix<float, 10, 1>::Zero();

    std::vector<hybridPresentation*> hybridPresentations;

    
    

    

    CSGtree finalTree;

    box_hybrid box(PSGMContext.m_context, 0, 0, 0, 2, 2, 1);
    cylinder_hybrid cylinder(PSGMContext.m_context, 1, 1, -1, 4, 0.5);
    sphere_hybrid sphere1(PSGMContext.m_context, 0, 0, 0, 1);
    sphere_hybrid sphere2(PSGMContext.m_context, 1, 0, 0, 1);
    goto sixball;

{
    Eigen::Matrix<float, 10, 1> xEqu0 = Eigen::Matrix<float, 10, 1>::Zero();
    xEqu0<<0,0,0,1,0,0,0,0,0,0;
    v1 = (*sphere1.functionsPointer)[0];
    auto q1 = vcf2gmp2bigintmatrix(v1);
    v2 = (*sphere2.functionsPointer)[0];
    auto q2 = vcf2gmp2bigintmatrix(v2);
    quad_inter<bigint> qi_ic = intersection(q1, q2, true, std::cout);
    QIOutputter outputter;
    QIWriter *writer = NULL;
    writer = new QIConsoleWriter();
    outputter.output(qi_ic, q1, q2);
    auto aasdadadadas = outputter.getOutput();
    std::cout << aasdadadadas->q1 << endl;
    std::cout << aasdadadadas->q2 << endl;
    std::cout << aasdadadadas->n_components << endl; // 相交分支数
    auto qweqwe = aasdadadadas->n_components;

    std::cout << "eaqution after intersection" << std::endl;
    for (int k = 0; k < qweqwe; k++)
    {
        std::cout << "{" + aasdadadadas->parametrizations[k].label + "}" << endl;
        string delta = "";
        auto ajskdhasjkdhakj = aasdadadadas->parametrizations[k].param;
        std::cout << ajskdhasjkdhakj << std::endl;
        auto spilt = spilt4ginac(ajskdhasjkdhakj, delta);

        GiNaC::parser reader(table);
        GiNaC::ex temp;
        std::vector<GiNaC::ex> bgelaqiwro;

        std::cout << "eaqution in ginac" << std::endl;
        for (int i = 0; i < spilt.size(); i++)
        {
            temp = reader(spilt[i]);
            std::cout << temp << " ";
            bgelaqiwro.push_back(temp);
        }
        std::cout << std::endl;

        //alg2(bgelaqiwro,u,vcf2gmp2bigintmatrix(xEqu0));
        //alg2(bgelaqiwro,u,q1);

        // all_curve_equation.emplace_back(bgelaqiwro, u, Delta, 1, 1);
    }
    std::cout << std::endl;
    return 0;
}

{
    Eigen::Matrix<float, 10, 1> xEqu0 = Eigen::Matrix<float, 10, 1>::Zero();
    xEqu0<<0,0,0,1,0,0,0,0,0,0;
    v1 = (*box.functionsPointer)[0];
    auto q1 = vcf2gmp2bigintmatrix(v1);
    v2 = (*cylinder.functionsPointer)[2];
    auto q2 = vcf2gmp2bigintmatrix(v2);
    quad_inter<bigint> qi_ic = intersection(q1, q2, true, std::cout);
    QIOutputter outputter;
    QIWriter *writer = NULL;
    writer = new QIConsoleWriter();
    outputter.output(qi_ic, q1, q2);
    auto aasdadadadas = outputter.getOutput();
    std::cout << aasdadadadas->q1 << endl;
    std::cout << aasdadadadas->q2 << endl;
    std::cout << aasdadadadas->n_components << endl; // 相交分支数
    auto qweqwe = aasdadadadas->n_components;

    std::cout << "eaqution after intersection" << std::endl;
    for (int k = 1; k < qweqwe; k++)
    {
        std::cout << "{" + aasdadadadas->parametrizations[k].label + "}" << endl;
        string delta = "";
        auto ajskdhasjkdhakj = aasdadadadas->parametrizations[k].param;
        std::cout << ajskdhasjkdhakj << std::endl;
        auto spilt = spilt4ginac(ajskdhasjkdhakj, delta);

        GiNaC::parser reader(table);
        GiNaC::ex temp;
        std::vector<GiNaC::ex> bgelaqiwro;

        std::cout << "eaqution in ginac" << std::endl;
        for (int i = 0; i < spilt.size(); i++)
        {
            temp = reader(spilt[i]);
            std::cout << temp << " ";
            bgelaqiwro.push_back(temp);
        }
        std::cout << std::endl;

        //alg2(bgelaqiwro,u,vcf2gmp2bigintmatrix(xEqu0));
        //alg2(bgelaqiwro,u,q1);

        // all_curve_equation.emplace_back(bgelaqiwro, u, Delta, 1, 1);
       
    }
    std::cout << std::endl;
    
}


    std::cout << "in\n";
    // for (int i = 0; i < (*sphere1.functionsPointer).size(); i++)
    // {
    //     v1 = (*sphere1.functionsPointer)[i];
    //     auto q1 = vcf2gmp2bigintmatrix(v1);
    //     for (int j = 0; j < (*sphere2.functionsPointer).size(); j++)
    //     {
    //         v2 = (*sphere2.functionsPointer)[j];
    //         auto q2 = vcf2gmp2bigintmatrix(v2);
    //         std::cout << q1 << std::endl;
    //         std::cout << q1(0, 0) << std::endl;
    //         auto qrqrqrqqqq = q1(0, 0);
    //         quad_inter<bigint> qi_ic = intersection(q1, q2, true, std::cout);
    //         QIOutputter outputter;
    //         QIWriter *writer = NULL;
    //         writer = new QIConsoleWriter();
    //         outputter.output(qi_ic, q1, q2);
    //         auto aasdadadadas = outputter.getOutput();
    //         std::cout << aasdadadadas->q1 << endl;
    //         std::cout << aasdadadadas->q2 << endl;
    //         std::cout << aasdadadadas->n_components << endl; // 相交分支数
    //         auto qweqwe = aasdadadadas->n_components;

    //         std::cout << "eaqution after intersection" << std::endl;
    //         for (int k = 0; k < qweqwe; k++)
    //         {
    //             std::cout << "{" + aasdadadadas->parametrizations[k].label + "}" << endl;
    //             string delta = "";
    //             auto ajskdhasjkdhakj = aasdadadadas->parametrizations[k].param;
    //             std::cout << ajskdhasjkdhakj << std::endl;
    //             auto spilt = spilt4ginac(ajskdhasjkdhakj, delta);

    //             GiNaC::parser reader(table);
    //             GiNaC::ex temp;
    //             std::vector<GiNaC::ex> bgelaqiwro;

    //             std::cout << "eaqution in ginac" << std::endl;
    //             for (int i = 0; i < spilt.size(); i++)
    //             {
    //                 temp = reader(spilt[i]);
    //                 std::cout << temp << " ";
    //                 bgelaqiwro.push_back(temp);
    //             }
    //             std::cout << std::endl;

    //             all_curve_equation.emplace_back(bgelaqiwro, u, Delta, 1, 1);
    //         }
    //         std::cout << std::endl;
    //     }
    // }



    // for (int i = 0; i < (*box.functionsPointer).size(); i++)
    // {
    //     v1 = (*box.functionsPointer)[i];
    //     auto q1 = vcf2gmp2bigintmatrix(v1);
    //     for (int j = 0; j < (*cylinder.functionsPointer).size(); j++)
    //     {
    //         v2 = (*cylinder.functionsPointer)[j];
    //         auto q2 = vcf2gmp2bigintmatrix(v2);
    //         quad_inter<bigint> qi_ic = intersection(q1, q2, true, std::cout);
    //         QIOutputter outputter;
    //         QIWriter *writer = NULL;
    //         writer = new QIConsoleWriter();
    //         outputter.output(qi_ic, q1, q2);
    //         auto aasdadadadas = outputter.getOutput();
    //         std::cout << aasdadadadas->q1 << endl;
    //         std::cout << aasdadadadas->q2 << endl;
    //         std::cout << aasdadadadas->n_components << endl; // 相交分支数
    //         auto qweqwe = aasdadadadas->n_components;
    //         for (int k = 0; k < qweqwe; k++)
    //         {
    //             std::cout << "{" + aasdadadadas->parametrizations[k].label + "}" << endl;
    //             string delta = "";
    //             auto ajskdhasjkdhakj = aasdadadadas->parametrizations[k].param;
    //             std::cout << ajskdhasjkdhakj << std::endl;
    //             // auto spilt = spilt4ginac(ajskdhasjkdhakj, delta);

    //             // GiNaC::parser reader(table);
    //             // GiNaC::ex temp;
    //             // std::vector<GiNaC::ex> bgelaqiwro;

    //             // for (int i = 0; i < spilt.size(); i++)
    //             // {
    //             //     temp = reader(spilt[i]);
    //             //     std::cout << temp << " ";
    //             //     bgelaqiwro.push_back(temp);
    //             // }
    //             // std::cout << std::endl;

    //             // all_curve_equation.emplace_back(bgelaqiwro, u, Delta, 1, 1);
    //         }
    //         std::cout << std::endl;
    //     }
    // }

    // std::cout<<"box=\n";
    // for (int i = 0; i < (*box.functionsPointer).size(); i++)
    // {
    //     std::cout << (*box.functionsPointer)[i]<<std::endl;
    // }
    // std::cout<<"cylinder=\n";
    // for (int i = 0; i < (*cylinder.functionsPointer).size(); i++)
    // {
    //     std::cout << (*cylinder.functionsPointer)[i]<<std::endl;
    // }
    // PS_CLOUD_VIEW_MODEL(PSGMContext.m_session);
    std::cout << "out\n";

    std::cout << "out\n";
    return 0;


sixball:
    bigint_matrix q1;
    bigint_matrix q2;
    bigint_matrix q3; // test

    // Eigen::Matrix<float, 10, 1> v1 = Eigen::Matrix<float, 10, 1>::Zero();
    // Eigen::Matrix<float, 10, 1> v2 = Eigen::Matrix<float, 10, 1>::Zero();
    Eigen::Matrix<float, 10, 1> v3 = Eigen::Matrix<float, 10, 1>::Zero();
    Eigen::Matrix<float, 10, 1> v4 = Eigen::Matrix<float, 10, 1>::Zero();
    Eigen::Matrix<float, 10, 1> v5 = Eigen::Matrix<float, 10, 1>::Zero();
    Eigen::Matrix<float, 10, 1> v6 = Eigen::Matrix<float, 10, 1>::Zero();

    std::vector<Eigen::Matrix<float, 10, 1>> vectorV;

    // v1<<1 ,0, 0, 0, 1, 0, 0, 1, 0, -1;
    // v2<<2 ,0, 0, 0, 4, 0, 0, 0, 0, -1;

    // no intersect
    // v1<<1,0,0,0,1,0,0,1,0,-1;
    // v2<<2,0,0,0,4,0,0,3,0,-1;

    // easy case intersect with point
    //  v1<<-1,0,0,0,0,0,2,-1,0,0;
    //  v2<<-3,0,0,0,1,0,0,-1,0,0;

    // case
    // v1<<0,0.02,0,0,0,0,0,0,0.02,0;
    // v2<<-2,0,0,0,2,0,0,-4,2,4;

    v1 << -1, 0, 0, 0, 0, 0, 2, -1, 0, 0;
    v2 << -3, 0, 0, 0, 1, 0, 0, -1, 0, 0;

    v1 << 1, 0, 0, 0, 1, 0, 0, 1, 0, -1;
    v2 << 2, 0, 0, 0, 2, 0, 0, 0, 0, -1;

    v1 << 1, 0, 0, 0, -1, 0, 4, 1, 0, 0;
    v2 << -3, 0, 0, 0, 1, 0, 0, 1, 0, 0;

    v1 << 1, 0, 0, 0, 1, 0, 0, 1, 0, -1;
    v2 << 2, 0, 0, 0, 4, 0, 0, 0, 0, -1;

    // type3 output case
    //  v1<<1,0,0,0,-1,0,4,1,0,0;
    //  v2<<-3,0,0,0,1,0,0,1,0,0;

    // type1 output case
    //  v1<<1,0,0,0,1,0,0,1,0,-1;
    //  v2<<2,0,0,0,4,0,0,0,0,-1 ;

    v3 << 1, 2, 0, 0, 1, 0, 0, 1, 0, -1;
    // v3<<0, 2, 0, 0, 0, 0, 0, 0, 2, 0;

    // six sphere
    v1 << 4, 0, 0, 4, 4, 0, 0, 4, 0, -3;
    v2 << 4, 0, 0, -4, 4, 0, 0, 4, 0, -3;
    v3 << 4, 0, 0, 0, 4, 0, 4, 4, 0, -3;
    v4 << 4, 0, 0, 0, 4, 0, -4, 4, 0, -3;
    v5 << 4, 0, 0, 0, 4, 0, 0, 4, 4, -3;
    v6 << 4, 0, 0, 0, 4, 0, 0, 4, -4, -3;

    // v1 << 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
    // v2 << 0, 0, 0, 1, 0, 0, 0, 0, 0, 1;
    // v3 << 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;

//    v1<<1,0,0,0,1,0,0,1,0,-1;
// v2<<2,0,0,0,4,0,0,0,0,-1;
    // v1<<1,0,0,0,1,0,0,1,0,-1;



    vectorV.push_back(v1);
    vectorV.push_back(v2);
    vectorV.push_back(v3);
    vectorV.push_back(v4);
    vectorV.push_back(v5);
    vectorV.push_back(v6);
    // vectorV.resize(3);
    // vectorV.resize(2);
    // vectorV.resize(1);
    vectorV.resize(4);
    vectorV.resize(0);

    //load test case
    InitParams initParams(PSGMContext.m_context, &hybridPresentations, &vectorV,&finalTree,{});
    initCase(initParams,"两个一样大的sphere");
    finalTree=*(initParams.finalTree);

    // std::cout<<"qMatrixIndexToPSGMsurfaceSPtrSPtr.size = "<<qMatrixIndexToPSGMsurfaceSPtrSPtr.size()<<std::endl;
    // std::cout<<"hybridPresentations.size = "<<hybridPresentations.size()<<std::endl;
    // std::cout<<"vectorV.size = "<<vectorV.size()<<std::endl;
    // std::cout<<"tree.type = "<<finalTree.getType()<<std::endl;
    // return 0;

    Viewer.addSurface(vectorV);

    //add CSG
    // sphere_primitive s1(-0.5,0,0,1);
    // sphere_primitive s2(0.5,0,0,1);
    // sphere_primitive s3(0,-0.5,0,1);
    // sphere_primitive s4(0,0.5,0,1);
    // CSGtree ss1(&s1);
    // CSGtree ss2(&s2);
    // CSGtree ss3(&s3);
    // CSGtree ss4(&s4);
    // CSGtree s1s2(CSGoperation::INTERSECTION,&ss1,&ss2);
    // CSGtree s1s2s3(CSGoperation::INTERSECTION,&s1s2,&ss3);
    // CSGtree s1s2s3s4(CSGoperation::INTERSECTION,&s1s2s3,&ss4);

    // finalTree(&s1s2s3s4);
    // finalTree=s1s2s3s4;
    // finalTree.print();
    // Viewer.addCSG(s1s2s3s4);
    Viewer.addCSG(finalTree);

    std::vector<bigint_matrix> Q;
    for (int i = 0; i < vectorV.size(); i++)
    {
        bigint_matrix q;
        q = vcf2gmp2bigintmatrix(vectorV[i]);
        Q.push_back(q);
    }

    int point_cnt = 100;
    double step = 2 * M_PI / point_cnt;
    float pointSize = 20;
    viewer.data().point_size=pointSize;
    
    // line_width is NOT SUPPORTED on Mac OS and Windows
    viewer.data().line_width=50.0f;
    
    int pqrowrqr = -2;

    // draw ssi_curve
    int cnt_renderCurve=3;
    int now_renderCurve=0;

    std::map<pair<unsigned int,unsigned int>,pair<unsigned int,unsigned int>> firstSecondToBeginEnd;
    quad_inter<bigint> qi_ic;
    QIOutputter outputter;
    outputter.setOmitImaginaryParametrizations();
    QIWriter *writer = NULL;
    writer = new QIConsoleWriter();
    QIOutputInfo* aasdadadadas;

    
    
    
    Timer ECQJ("ECQJ");
    for (unsigned int first_index = 0; first_index < Q.size(); first_index++)
    {
        for (unsigned int second_index = first_index + 1; second_index < Q.size(); second_index++)
        {
            // std::cout<<"first_index= "<<first_index+1<<" second_index= "<<second_index+1<<endl;
            pair<unsigned int,unsigned int> beginEnd=intersectionDouble(Q,firstSecondToBeginEnd,first_index,second_index,ginacParams,all_curve_equation,all_vertex_equation);
            // pair<unsigned int,unsigned int> beginEnd;
            // pair<unsigned int,unsigned int> firstSecond={first_index,second_index};
            // if(firstSecondToBeginEnd.find(firstSecond)!=firstSecondToBeginEnd.end())
            // {
            //     beginEnd=firstSecondToBeginEnd[firstSecond];
            // }
            // else
            // {
            //     beginEnd.first=all_curve_equation.size();

            //     qi_ic = intersection(Q[first_index], Q[second_index], true, std::cout);
            //     outputter.output(qi_ic, Q[first_index], Q[second_index]);

            //     // writer->setOutputInformation(outputter.getOutput());
            //     // writer->setVerbosityLevel(3);
            //     // writer->write ();/** Produces the result to /dev/stdout by default */

            //     // QIOutputParametrization parametrizations[4]; /** List of parametrizations (one for each algebraic
            //     // 						   component of the intersection) */

            //     // aasdadadadas = outputter.getOutput();
            //     // handleCurveComponent(aasdadadadas,ginacParams,all_curve_equation,all_vertex_equation,first_index,second_index);
            //     handleCurveComponent(outputter.getOutput(),ginacParams,all_curve_equation,all_vertex_equation,first_index,second_index);
                
            //     beginEnd.second=all_curve_equation.size();
            //     firstSecondToBeginEnd[firstSecond]=beginEnd;
            // }
            
            //TSQI 
            //get three QuadSurface intersection point 
            for (unsigned int third_index = second_index + 1; third_index < Q.size(); third_index++)
            {
                // std::cout<<"third_index= "<<third_index+1<<endl;
                // std::cout<<"begin= "<<beginEnd.first<<" End= "<<beginEnd.second<<endl;
                // TQSI::seq seqST= intersectionTriple(Q[first_index], Q[second_index],Q[third_index]);

                for(int begin=beginEnd.first;begin<beginEnd.second;begin++)
                {
                    TQSI::seq seqST= intersectionTriple(all_curve_equation[begin],Q[third_index]);
                    curveAddSeqToPointList(all_curve_equation[begin],seqST);
                }
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
                    beginEnd.first=all_curve_equation.size();

                    qi_ic = intersection(Q[first_index], Q[third_index], true, std::cout);
                    outputter.output(qi_ic, Q[first_index], Q[third_index]);
                    aasdadadadas = outputter.getOutput();
                    handleCurveComponent(aasdadadadas,ginacParams,all_curve_equation,all_vertex_equation,first_index,third_index);

                    beginEnd.second=all_curve_equation.size();
                    firstSecondToBeginEnd[firstThird]=beginEnd;
                }

                for(int begin=beginEnd.first;begin<beginEnd.second;begin++)
                {
                    TQSI::seq seqSU= intersectionTriple(all_curve_equation[begin],Q[second_index]);
                    curveAddSeqToPointList(all_curve_equation[begin],seqSU);
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
                    beginEnd.first=all_curve_equation.size();
                    
                    qi_ic = intersection(Q[second_index], Q[third_index], true, std::cout);
                    outputter.output(qi_ic, Q[second_index], Q[third_index]);
                    aasdadadadas = outputter.getOutput();
                    handleCurveComponent(aasdadadadas,ginacParams,all_curve_equation,all_vertex_equation,second_index,third_index);

                    beginEnd.second=all_curve_equation.size();
                    firstSecondToBeginEnd[secondThird]=beginEnd;
                }
                
                for(int begin=beginEnd.first;begin<beginEnd.second;begin++)
                {
                    TQSI::seq seqTU= intersectionTriple(all_curve_equation[begin],Q[first_index]);
                    curveAddSeqToPointList(all_curve_equation[begin],seqTU);
                }
                // divideCurve(all_curve_equation[beginEnd.first],&finalTree,seqTU);

                //取代无穷参数的点
                //有可能出现两个分支，就会有两个seq，到时需要match的话需要修改
                // matchThreeSeq(seqST,seqSU,seqTU);
                
                // cout<<"cout seq:"<<endl;
                // seqST.print();
                // seqSU.print();
                // seqTU.print();
                //deal one of seq
                //deal seqST to get vertex
            }

        }
    }
    ECQJ.~Timer();
    delete writer;

    cout<<"all_curve_equation.size = "<<all_curve_equation.size()<<endl;
    for(auto& curve:all_curve_equation)
    {
        sort(curve.pointList.begin(),curve.pointList.end(),[](TQSI::intPoint a,TQSI::intPoint b)
        {
            return a.parameter<b.parameter;
        });
        curve.pointList.erase(std::unique(curve.pointList.begin(), curve.pointList.end(), [](TQSI::intPoint a, TQSI::intPoint b)
                                { return a.parameter == b.parameter; }),
                    curve.pointList.end());
        
        divideCurve(curve,&finalTree);
        // for(auto point:curve.pointList)
        // {
        //     point.print();
        // }
        // cout<<"done"<<endl;
    }

    Viewer.addCurve(all_curve_equation);
    // Viewer.addCurve(all_curve_equation[4]);

    for(auto& curve:all_curve_equation)
    {
        //Otherwise, use the endpoint directly
        if(curve.isLoop&&curve.isAllInterval)
        {
            // handleRingEdge();
        }
        else
        {
            
        }
        auto p0=curve.getPoint(0,1,1);
        auto p1=curve.getPoint(1,1,1);
        PSGMPoint pp0{p0.x,p0.y,p0.z};
        PSGMPoint pp1{p1.x,p1.y,p1.z};
        auto vec=pp1-pp0;
        // auto vec=p1-p0;
        auto fi=curve.first_index;
        auto si=curve.second_index;

        auto face1SPtr=(*qMatrixIndexToPSGMsurfaceSPtrSPtr)[fi];
        auto face2SPtr=(*qMatrixIndexToPSGMsurfaceSPtrSPtr)[si];
        // auto surface1=face1SPtr->getEntSurface()->getSurface();
        // auto surface2=face2SPtr->getEntSurface()->getSurface();
        auto norm1=face1SPtr->getNormal(face1SPtr->getParamAt(pp0));
        auto norm2=face2SPtr->getNormal(face2SPtr->getParamAt(pp0));
        // auto norm1=surface1->getNormal(surface1->getParamAt(pp0));
        // auto norm2=surface2->getNormal(surface2->getParamAt(pp1));

        bool sense1 = norm1.cross(norm2).dot(vec) > 0;
        bool sense2 = !sense1;
        std::cout << sense1 << std::endl;
        std::cout << sense2 << std::endl;
    }

    // for(int i=0;i<all_curve_equation.size();i++)
    // {
    //     cout<<"curve"<<i<<" = "<<endl;
    //     all_curve_equation[i].print();
    // }

    // Eigen::MatrixXd mat;
    // sampleCurvePoint(all_curve_equation[0],mat,step);
    // viewer.data().add_points(mat, Eigen::RowVector3d(0.0, 0.0, 1.0));

    Eigen::MatrixXd point0(1, 3);
    Eigen::MatrixXd point1(1, 3);

    auto p0=all_curve_equation[0].getPoint(0,1);
    auto p1=all_curve_equation[0].getPoint(200,1);
    // auto p2=all_curve_equation[1].getPoint(0.33333,1);
    // auto p3=all_curve_equation[1].getPoint(200,1);
    // auto p4=all_curve_equation[1].getPoint(1000,1);
   
    
    auto neg100=all_curve_equation[0].getPoint(-100,1);
    auto pos100=all_curve_equation[0].getPoint(100,1);
    


    // auto negative=s1s2s3.get_sdf(neg100[0],neg100[1],neg100[2]);
    // auto positive=s1s2s3.get_sdf(pos100[0],pos100[1],pos100[2]);
    // PRINT_VAR(negative);PRINT_VAR(positive);
    
    point0<<p0[0],p0[1],p0[2];
    point1<<p1[0],p1[1],p1[2];


    Viewer.addPoint(point0);
    Viewer.addPoint(point1);

    // viewer.data().add_points(point0, Eigen::RowVector3d(0.0, 1.0, 0));
    // viewer.data().add_points(point1, Eigen::RowVector3d(0.0, 0.5, 0.5));
    
    // ttpoint<<0,-0.005,-0.866;
    // viewer.data().add_points(ttpoint, Eigen::RowVector3d(0.0, 1.0, 0));

    
    // Viewer.setVisualizateSurface(true);
    Viewer.setVisualizateCurve(true);
    Viewer.setVisualizatePoint(true);
    Viewer.setVisualizateCSG(true);
    Viewer.visualizate();
    // Viewer.drawCurveByPoint();
    
    // viewer.launch();

    // PS_CLOUD_VIEW_MODEL(PSGMContext.m_session);


    return 0;
}