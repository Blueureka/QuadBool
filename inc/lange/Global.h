#pragma once

#include <memory>
#include <vector>
#include <map>

#include <ginac/ginac.h>
#include<lange/AABB.h>


#include "PSUtils.h"

struct vertex_equation;
struct curve_equation;
struct GinacParams;

class GlobalVariable
{
public:
    // 获取单例实例
    static GlobalVariable &getInstance()
    {
        static GlobalVariable instance; // 在第一次调用时创建唯一的实例
        return instance;
    }

    // 防止拷贝和赋值
    GlobalVariable(const GlobalVariable &) = delete;
    void operator=(const GlobalVariable &) = delete;

    static std::shared_ptr<std::vector<curve_equation>> &getAllCurveEquation()
    {
        return getInstance().all_curve_equation;
    }
    static std::shared_ptr<std::vector<vertex_equation>> &getAllVertexEquation()
    {
        return getInstance().all_vertex_equation;
    }
    static std::shared_ptr<std::vector<PSGMSurfaceSPtr>> &getQMatrixIndexToPSGMsurfaceSPtrSPtr()
    {
        return getInstance().qMatrixIndexToPSGMsurfaceSPtrSPtr;
    }
    static std::shared_ptr<std::vector<Eigen::Matrix<float, 10, 1>>> &getImplicitMatrixVector()
    {
        return getInstance().implicitMatrixVector;
    }
    static std::shared_ptr<std::vector<bigint_matrix>> &getQMatrix()
    {
        return getInstance().QMatrix;
    }
    static std::shared_ptr<std::map<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, unsigned int>>> &getFirstSecondToBeginEnd()
    {
        return getInstance().firstSecondToBeginEnd;
    }
    static std::shared_ptr<GinacParams> &getGinacParamsSPtr()
    {
        return getInstance().ginacParamsSPtr;
    }
    static GiNaC::symbol& getU()
    {
        return getInstance().u;
    }
    static GiNaC::symbol& getV()
    {
        return getInstance().v;
    }
    static GiNaC::symbol& getS()
    {
        return getInstance().s;
    }
    static GiNaC::symbol& getT()
    {
        return getInstance().t;
    }
    static GiNaC::symbol& getDelta()
    {
        return getInstance().Delta;
    }

    static std::shared_ptr<GiNaC::symtab>& getTableSPtr()
    {
        return getInstance().tableSPtr;
    }

    static std::shared_ptr<std::vector<AABB>>& getQAABBSPtr()
    {
        return getInstance().QAABBSPtr;
    }

    static std::shared_ptr<int> getTest()
    {
        return getInstance().test;
    }

private:
    // 私有构造函数，防止外部直接实例化
    // 默认构造函数初始化
    GlobalVariable() : all_curve_equation(std::make_shared<std::vector<curve_equation>>()), all_vertex_equation(std::make_shared<std::vector<vertex_equation>>()), qMatrixIndexToPSGMsurfaceSPtrSPtr(std::make_shared<std::vector<PSGMSurfaceSPtr>>()), implicitMatrixVector(std::make_shared<std::vector<Eigen::Matrix<float, 10, 1>>>()),QMatrix(std::make_shared<std::vector<bigint_matrix>>()),firstSecondToBeginEnd(std::make_shared<std::map<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, unsigned int>>>() ),ginacParamsSPtr(std::make_shared<GinacParams>()),test(std::make_shared<int>(0)),u("u"),v("v"),s("s"),t("t"),Delta("Delta") 
    {
        (*tableSPtr)["u"] = u;
        (*tableSPtr)["v"] = v;
        (*tableSPtr)["s"] = s;
        (*tableSPtr)["t"] = t;
        (*tableSPtr)["Delta"] = Delta;
    }

    // 其他成员函数和成员变量...
    std::shared_ptr<std::vector<curve_equation>> all_curve_equation;
    std::shared_ptr<std::vector<vertex_equation>> all_vertex_equation;
    
    //find one of corresponding PSGMsurfaceSPtr,and use PS's function
    std::shared_ptr<std::vector<PSGMSurfaceSPtr>> qMatrixIndexToPSGMsurfaceSPtrSPtr;
    std::shared_ptr<std::vector<Eigen::Matrix<float, 10, 1>>> implicitMatrixVector;

    std::shared_ptr<std::vector<bigint_matrix>> QMatrix;
    std::shared_ptr<std::vector<AABB>> QAABBSPtr=std::make_shared<std::vector<AABB>>();
    std::shared_ptr<std::map<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, unsigned int>>> firstSecondToBeginEnd;
    std::shared_ptr<GinacParams> ginacParamsSPtr;

    GiNaC::symbol u;
    GiNaC::symbol v;
    GiNaC::symbol s;
    GiNaC::symbol t;
    GiNaC::symbol Delta;
    std::shared_ptr<GiNaC::symtab> tableSPtr=std::make_shared<GiNaC::symtab>();

    std::shared_ptr<int> test;
};

// 全局单例实例声明
// extern GlobalVariable& globalVariables;