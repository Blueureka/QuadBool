#pragma once

#include <vector>
#include <Eigen/Eigen>
#include <string>

#include <lange/Hybrid.h>
#include <nlohmann/json.hpp>

struct InitParams
{
    PSGMISessionContext *context;
    std::vector<hybridPresentation *> *hybridPresentations;
    std::vector<Eigen::Matrix<float, 10, 1>> *v;
    CSGtree* finalTree;
    std::map<unsigned int,PSGMTopoFace*>* vindexToFace;
    InitParams(PSGMISessionContext *_context, std::vector<hybridPresentation *> *_hybridPresentations, std::vector<Eigen::Matrix<float, 10, 1>> *_v, CSGtree* _finalTree,std::map<unsigned int,PSGMTopoFace*>* _vindexToFace) : context(_context), hybridPresentations(_hybridPresentations), v(_v), finalTree(_finalTree),vindexToFace(_vindexToFace)
    {
    }
};

int parserJson(InitParams &initParams, const std::string &filePath, const std::string &caseName);

void initCase(PSGMISessionContext *context, std::vector<hybridPresentation *> &hybridPresentations, std::vector<Eigen::Matrix<float, 10, 1>> &v, CSGtree &finalTree);
void initCase(InitParams &initParams, const std::string &caseName);