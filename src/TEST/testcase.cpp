#include <TEST/testcase.h>

using json = nlohmann::json;

int parserJson(InitParams& initParams,const std::string &filePath, const std::string &caseName)
{
    // 读取 JSON 文件
    ifstream ifs(filePath);
    if (!ifs.is_open())
    {
        std::cerr << "Failed to open the file." << std::endl;
        return 1;
    }
    else
    {
        std::cout << "Successfully open the file." << std::endl;
    }

    // 从文件流读取 JSON 数据
    json j;
    ifs >> j;
    // string caseName = "两个一样大的sphere";

    // 检查是否解析成功
    if (j.is_null())
    {
        cerr << "Failed to parse JSON." << endl;
        return 1;
    }
    else if ((j.contains(caseName) && j[caseName].is_object()))
    {
        std::cout << "have this case" << std::endl;
    }
    else
    {
        cerr << "no this case" << endl;
        return 1;
    }

    auto context=initParams.context;
    auto hybridPresentations=initParams.hybridPresentations;
    auto v=initParams.v;
    auto vindexToFace=initParams.vindexToFace;

    // auto finalTree=initParams.finalTree;
    // std::vector<CSGtree*> CSGvector;
    std::vector<primitive*> CSGvector;
    auto op=CSGoperation::UNION;

    // very important!!!
    j = j[caseName];

    // 从 JSON 中提取数据
    if (j.is_object())
    {
        if (j.contains("operation") && j["operation"].is_string())
        {
            std::string operation=j["operation"];
            if(operation=="UNION")
            {
                op=CSGoperation::UNION;
            }
            else if(operation=="INTERSECTION")
            {
                op=CSGoperation::INTERSECTION;
            }

        }
        // 提取body
        if (j.contains("bodys") && j["bodys"].is_array())
        {
            auto bodys = j["bodys"];
            cout << "bodys:" << endl;
            for (auto &body : bodys)
            {
                cout << "  " << body << endl;
                if (body.contains("type") && body["type"].is_string())
                {
                    auto type = body["type"];
                    if (type == "sphere")
                    {
                        double x = body["x"];
                        double y = body["y"];
                        double z = body["z"];
                        auto radius= body["radius"];
                        sphere_hybrid sphere1(context, x, y, z, radius);
                        hybridPresentation *ap = (hybridPresentation *)&sphere1;
                        hybridPresentations->push_back(ap);
                        // std::cout<<"in function hybridPresentations.size= "<<hybridPresentations->size()<<std::endl;
                        auto implictFunctionsVector = *sphere1.functionsPointer;
                        auto faceVector=*sphere1.facesPointer;

                        for(int i=0;i<implictFunctionsVector.size();i++)
                        {
                            unsigned int vindex=v->size();
                            v->push_back(implictFunctionsVector[i]);
                           (*(vindexToFace))[vindex]=faceVector[i];
                        }
                        auto caga = sphere1.tree;
                        CSGvector.push_back(caga);
                        // finalTree=CSGtree(*caga);
                        // std::cout << x << std::endl;
                    }
                }
            }
        }
        
    }

    primitive* tempPrim=CSGvector[0];
    // CSGtree tempCSG=CSGtree(*(CSGvector[0]));
    // tempCSG.print();
    CSGtree tempCSG;
    tempCSG=CSGtree(op,tempPrim,CSGvector[1]);
    // for(int i=2;i<CSGvector.size();i++)
    // {
    //     auto newCSG=CSGtree(op,&tempCSG,CSGvector[i]);
    //     tempCSG=newCSG;
    //     // tempCSG.print();
    // }
    initParams.finalTree=&tempCSG;
    return 0;
}

void initCase(InitParams& initParams,const std::string &caseName)
{
    
    if (parserJson(initParams,"../src/TEST/example.json", "四个一样大的sphere"))
    {
        // 抛异常
        return;
    }
    return;

    PSGMISessionContext *context = initParams.context;
    std::vector<hybridPresentation *> hybridPresentations = *(initParams.hybridPresentations);
    std::vector<Eigen::Matrix<float, 10, 1>> v = *(initParams.v);
    CSGtree* finalTree = initParams.finalTree;
    v.resize(0);


    return ;
    // two sphere union

    sphere_hybrid sphere1(context, -0.5, 0, 0, 1);
    sphere_hybrid sphere2(context, 0.5, 0, 0, 1);

    hybridPresentation *ap = (hybridPresentation *)&sphere1;
    hybridPresentation *bp = (hybridPresentation *)&sphere2;
    hybridPresentations.push_back(ap);
    hybridPresentations.push_back(bp);

    auto a = *sphere1.functionsPointer;
    for (auto &vv : a)
    {
        v.push_back(vv);
    }

    auto b = *sphere2.functionsPointer;
    for (auto &vv : b)
    {
        v.push_back(vv);
    }

    auto caga = sphere1.tree;
    auto cagb = sphere2.tree;
    finalTree = new CSGtree(CSGoperation::UNION, caga, cagb);
}