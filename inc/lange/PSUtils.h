#pragma once

// include all head file which i need

#include "algorithm/PSGMAlgTopologyTypes.h"
#include "algorithm/PSGMAlgorithm.h"
#include "algorithm/primitive/PSGMAlgPrimitiveModel.h"
#include "api/PSGMApi.h"

#include "model/topology/PSGMTopoBody.h"
#include "model/topology/PSGMTopoCoedge.h"
#include "model/topology/PSGMTopoEdge.h"
#include "model/topology/PSGMTopoFace.h"
#include "model/topology/PSGMTopoLoop.h"
#include "model/topology/PSGMTopoLump.h"
#include "model/topology/PSGMTopoShell.h"
#include "model/topology/PSGMTopoVertex.h"

#include "geometry/PSGMLCS.h"
#include "geometry/curve/PSGMCurveCircle.h"
#include "geometry/curve/PSGMCurveLine.h"
#include "geometry/surface/PSGMSurfCone.h"
#include "geometry/surface/PSGMSurfCylinder.h"
#include "geometry/surface/PSGMSurfPlane.h"
#include "geometry/surface/PSGMSurfSphere.h"
#include "geometry/surface/PSGMSurfTorus.h"

#include "model/geometry/curve/PSGMEntCurve.h"
#include "model/geometry/curve/PSGMEntCurveBSpline.h"
#include "model/geometry/curve/PSGMEntCurveCircle.h"
#include "model/geometry/curve/PSGMEntCurveEllipse.h"
#include "model/geometry/curve/PSGMEntCurveIntersect.h"
#include "model/geometry/curve/PSGMEntCurvePolyline.h"

#include "transaction/PSGMISession.h"

#include "geometry/PSGMLCS.h"
#include "geometry/PSGMPoint2D.h"
#include "geometry/curve/PSGMCurveCircle.h"
#include "geometry/curve/PSGMCurveLine.h"
#include "geometry/surface/PSGMSurfCone.h"
#include "geometry/surface/PSGMSurfCylinder.h"
#include "geometry/surface/PSGMSurfPlane.h"
#include "geometry/surface/PSGMSurfSphere.h"
#include "geometry/surface/PSGMSurfTorus.h"

#include "model/geometry/curve/PSGMEntCurve.h"
#include "model/geometry/curve/PSGMEntCurveBSpline.h"
#include "model/geometry/curve/PSGMEntCurveCircle.h"
#include "model/geometry/curve/PSGMEntCurveEllipse.h"
#include "model/geometry/curve/PSGMEntCurveIntersect.h"
#include "model/geometry/curve/PSGMEntCurvePolyline.h"

#include "model/geometry/PSGMEntGeometryConvertor.h"
#include "model/geometry/curve/PSGMEntCurveLine.h"
#include "model/geometry/surface/PSGMEntSurfBSpline.h"
#include "model/geometry/surface/PSGMEntSurfBlend.h"
#include "model/geometry/surface/PSGMEntSurfCone.h"
#include "model/geometry/surface/PSGMEntSurfCylinder.h"
#include "model/geometry/surface/PSGMEntSurfPlane.h"
#include "model/geometry/surface/PSGMEntSurfSphere.h"
#include "model/geometry/surface/PSGMEntSurfSpin.h"
#include "model/geometry/surface/PSGMEntSurfTorus.h"
#include "model/geometry/surface/PSGMEntSurface.h"

#include "model/topology/PSGMTopoBody.h"
#include "model/topology/PSGMTopoCoedge.h"
#include "model/topology/PSGMTopoEdge.h"
#include "model/topology/PSGMTopoFace.h"
#include "model/topology/PSGMTopoLoop.h"
#include "model/topology/PSGMTopoLump.h"
#include "model/topology/PSGMTopoShell.h"
#include "model/topology/PSGMTopoVertex.h"

/**
 * @file       PSGMCloudUtils.h
 * @brief      cloud utils for model viewing
 * @details
 * @author     qiuxinyun(qiuxinyun@poissonsoft.com)
 * @date       2023-09-08
 * @version    V1.0
 * @copyright  Copyright (c) 2022-2025  Shen Zhen Poisson Software Technology CO.,LTD
 */

#pragma once

#include <chrono>
#include <ctime>
#include <filesystem>
#include <sstream>
#ifdef _WIN32
#include <windows.h>
#else
#include <cstdlib>
#endif

#include "api/PSGMApi.h"
#include "api/PSGMApiSession.h"
#include "apiProxy/PSGMApiProxy.h"
#include "debug/PSGMIApiDebug.h"
#include "persistence/external/json/PSGMBrepJsonWriter.h"
#include "transaction/PSGMISession.h"

#ifdef ENABLE_PS_CLOUD_VIEW
#define PS_CLOUD_VIEW_MODEL(session)                                                                        \
    {                                                                                                       \
        std::string testCaseName = "ljl";                                                                   \
        testCaseName.erase(std::remove(testCaseName.begin(), testCaseName.end(), '/'), testCaseName.end()); \
        ModelViewUtils::view(session, testCaseName);                                                        \
    }
#else
#define PS_CLOUD_VIEW_MODEL(session)
#endif

#define PS_CLOUD_MODEL_SAVE_PATH "/3rdparty/ps-kernel/gm/vis"
#define PS_CLOUD_NGINX_START_SCRIPT_PATH "/gm/test/nginx/start.bat"
#define INDENT_SPACES(n) std::string(n, ' ')

// Callback function for renderTopoFacet
class GoCallbacks : public PSGMIApiGraphicOutputCallbacks
{
    struct Point
    {
        float m_xyz[3]{0.0, 0.0, 0.0};
    };

    struct Face
    {
        int m_faceTag{0};
        std::vector<Point> m_triangularPoints;
    };

    struct Line
    {
        int m_curveTag{0};
        std::vector<Point> m_points;
    };

    struct Body
    {
        int m_bodyTag{0};
        std::vector<Face> m_faces;
        std::vector<Line> m_curves;
    };

public:
    GoCallbacks()
    {
    }

    virtual void segmentError(const PSGMApiGOSegmentType& segmentType,
                              const int32_t tagCount,
                              const int32_t* tags,
                              const PSGMApiGOErrorType& errorType,
                              PSGMApiGOStatus* iFail) override
    {
        switch (errorType)
        {
        case PSGMApiGOErrorType::Unspecified:
            m_error = -1;
            break;
        case PSGMApiGOErrorType::Rubber:
            m_error = -2;
            break;
        default:
            m_error = -3;
            break;
        }

        *iFail = PSGMApiGOStatus::Continue;
    }

    virtual void segmentFacet(const int32_t tagCount,
                              const int32_t* tags,
                              const int32_t geomCount,
                              const double* geoms,
                              const int32_t occurenceCount,
                              const PSGMApiGOFacetInfo& facetInfo,
                              PSGMApiGOStatus* iFail) override
    {
        Point p;

        for (int i = 0; i < 3; ++i)
        {
            int posOffset = 3 * i;
            p.m_xyz[0] = float(geoms[posOffset]);
            p.m_xyz[1] = float(geoms[posOffset + 1]);
            p.m_xyz[2] = float(geoms[posOffset + 2]);
            m_currFace.m_triangularPoints.push_back(p);
        }

        *iFail = PSGMApiGOStatus::Continue;
    }

    virtual void segmentLine(const PSGMApiGOSegmentType& segmentType,
                             const int32_t tagCount,
                             const int32_t* tags,
                             const int32_t geomCount,
                             const double* geoms,
                             const int32_t occurenceCount,
                             const PSGMApiGOLineInfo& lineInfo,
                             PSGMApiGOStatus* iFail) override
    {
        // sizeof geoms = 3 * geomCount, this is compatibal with parasolid
        for (int i = 0; i < geomCount; ++i)
        {
            Point p;
            int posOffset = 3 * i;
            p.m_xyz[0] = double(geoms[posOffset]);
            p.m_xyz[1] = double(geoms[posOffset + 1]);
            p.m_xyz[2] = double(geoms[posOffset + 2]);
            m_currCurve.m_points.push_back(p);
        }

        *iFail = PSGMApiGOStatus::Continue;
        return;
    }

    virtual void openSegment(const PSGMApiGOSegmentType& segmentType,
                             const int32_t tagCount,
                             const int32_t* tags,
                             const int32_t geomCount,
                             const double* geoms,
                             const int32_t lineTypeCount,
                             const int32_t* lineTypes,
                             PSGMApiGOStatus* iFail) override
    {
        ++m_deep;
        m_processing = true;

        onOpenClose(segmentType, tags);
        *iFail = PSGMApiGOStatus::Continue;
    }

    virtual void closeSegment(const PSGMApiGOSegmentType& segmentType,
                              const int32_t tagCount,
                              const int32_t* tags,
                              const int32_t geomCount,
                              const double* geoms,
                              const int32_t lineTypeCount,
                              const int32_t* lineTypes,
                              PSGMApiGOStatus* iFail) override
    {
        --m_deep;
        onOpenClose(segmentType, tags);
        *iFail = PSGMApiGOStatus::Continue;
    }

    /**
     * @brief Wait to process all segment
     * @return int
     */
    int proc()
    {
        int count = 0;
        while (count < 3000)
        {
            if (m_processing && m_deep == 0)
            {
                // process over
                break;
            }
            else
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }

            count++;
        }
        return 0;
    }

    /**
     * @brief Get serialized json data
     * @return string
     */
    std::string getFacetJson()
    {
        bool remain = false;
        if (m_faces.size() > 0)
        {
            m_currBody.m_faces.clear();
            std::copy(m_faces.begin(), m_faces.end(), std::back_inserter(m_currBody.m_faces));
            m_faces.clear();
            remain = true;
        }

        if (m_curves.size() > 0)
        {
            m_currBody.m_curves.clear();
            std::copy(m_curves.begin(), m_curves.end(), std::back_inserter(m_currBody.m_curves));
            m_curves.clear();
            remain = true;
        }

        if (remain)
        {
            m_bodies.push_back(m_currBody);

            m_currBody.m_faces.clear();
            m_currBody.m_curves.clear();
            m_currBody.m_bodyTag = 0;
        }

        std::stringstream json;
        json << std::setprecision(15);
        json << "{\n";
        json << INDENT_SPACES(2) << "\"bodies\": [\n";

        for (int i = 0; i < m_bodies.size(); i++)
        {
            const Body& b = m_bodies[i];
            json << INDENT_SPACES(4) << "{\n";
            json << INDENT_SPACES(6) << "\"bodyTag\": { \"m_PSGMApiTag\": " << b.m_bodyTag << " },\n";
            json << INDENT_SPACES(6) << "\"faces\": [\n";

            for (int j = 0; j < b.m_faces.size(); j++)
            {
                const Face& f = b.m_faces[j];
                json << INDENT_SPACES(8) << "{\n";
                json << INDENT_SPACES(10) << "\"faceTag\": { \"m_PSGMApiTag\": " << f.m_faceTag << " },\n";
                json << INDENT_SPACES(10) << "\"triangularPoints\": [\n";

                for (int k = 0; k < f.m_triangularPoints.size(); k++)
                {
                    const Point& p = f.m_triangularPoints[k];
                    json << INDENT_SPACES(12) << "{\n";
                    json << INDENT_SPACES(14) << "\"xyz\": [" << p.m_xyz[0] << ", " << p.m_xyz[1] << ", " << p.m_xyz[2]
                         << "]\n";
                    json << INDENT_SPACES(12) << "}";

                    if (k != f.m_triangularPoints.size() - 1)
                    {
                        json << ",";
                    }

                    json << "\n";
                }

                json << INDENT_SPACES(10) << "]\n";
                json << INDENT_SPACES(8) << "}";

                if (j != b.m_faces.size() - 1)
                {
                    json << ",";
                }

                json << "\n";
            }

            json << INDENT_SPACES(6) << "\n],\n";
            json << INDENT_SPACES(6) << "\"curves\": [\n";
            for (int j = 0; j < b.m_curves.size(); j++)
            {
                const Line& f = b.m_curves[j];
                json << INDENT_SPACES(8) << "{\n";
                json << INDENT_SPACES(10) << "\"curveTag\": { \"m_PSGMApiTag\": " << f.m_curveTag << " },\n";
                json << INDENT_SPACES(10) << "\"points\": [\n";

                for (int k = 0; k < f.m_points.size(); k++)
                {
                    const Point& p = f.m_points[k];
                    json << INDENT_SPACES(12) << "{\n";
                    json << INDENT_SPACES(14) << "\"xyz\": [" << p.m_xyz[0] << ", " << p.m_xyz[1] << ", " << p.m_xyz[2]
                         << "]\n";
                    json << INDENT_SPACES(12) << "}";

                    if (k != f.m_points.size() - 1)
                    {
                        json << ",";
                    }

                    json << "\n";
                }

                json << INDENT_SPACES(10) << "]\n";
                json << INDENT_SPACES(8) << "}";

                if (j != b.m_curves.size() - 1)
                {
                    json << ",";
                }

                json << "\n";
            }

            json << INDENT_SPACES(6) << "]\n";
            json << INDENT_SPACES(4) << "}";

            if (i != m_bodies.size() - 1)
            {
                json << ",";
            }

            json << "\n";
        }

        json << INDENT_SPACES(2) << "]\n";
        json << "}\n";

        return json.str();
    }

    int getError() const
    {
        return m_error;
    }

    const std::vector<Body>& getResult()
    {
        return m_bodies;
    }

private:
    /**
     * @brief Process GOSegment
     * @param[in] segmentType GOSegment type
     * @param[in] tags tags of entities
     * @return
     */
    void onOpenClose(const PSGMApiGOSegmentType& segmentType, const int32_t* tags)
    {
        // std::cout << "segmentType:" << (int)segmentType << " tag:" << *tags << std::endl;
        if (segmentType == PSGMApiGOSegmentType::Face)
        {
            if (m_currFace.m_faceTag > 0 && m_currFace.m_triangularPoints.size() > 0)
            {
                m_faces.push_back(m_currFace);
            }

            m_currFace.m_triangularPoints.clear();
            m_currFace.m_faceTag = *tags;
        }
        else if (segmentType == PSGMApiGOSegmentType::Body)
        {
            if (m_currBody.m_bodyTag > 0)
            {
                if (m_faces.size() > 0)
                {
                    m_currBody.m_faces.clear();
                    std::copy(m_faces.begin(), m_faces.end(), std::back_inserter(m_currBody.m_faces));
                    m_faces.clear();
                }

                if (m_curves.size() > 0)
                {
                    m_currBody.m_curves.clear();
                    std::copy(m_curves.begin(), m_curves.end(), std::back_inserter(m_currBody.m_curves));
                    m_curves.clear();
                }

                m_bodies.push_back(m_currBody);
            }
            m_currBody.m_faces.clear();
            m_currBody.m_curves.clear();
            m_currBody.m_bodyTag = *tags;
        }
        else if (segmentType == PSGMApiGOSegmentType::Edge)
        {
            m_currCurve.m_curveTag = *tags;

            if (m_currCurve.m_curveTag > 0 && m_currCurve.m_points.size() > 0)
            {
                Line item;
                item.m_curveTag = m_currCurve.m_curveTag;
                std::copy(m_currCurve.m_points.begin(), m_currCurve.m_points.end(), std::back_inserter(item.m_points));
                m_curves.push_back(item);
            }

            m_currCurve.m_points.clear();
            // m_currCurve.m_curveTag = *tags;
        }
    }

private:
    std::vector<Face> m_faces;
    std::vector<Body> m_bodies;
    std::vector<Line> m_curves;
    Line m_currCurve;
    Face m_currFace;
    Body m_currBody;

    // open add 1 close sub 1
    int m_deep{0};
    // open set to true
    bool m_processing{false};
    int m_error{0};
};

class ModelViewUtils
{
public:
    /**
     * @brief Open a web page to view model with default browser, only support windows platform
     * @tparam T psgmstd::shared_ptr<PSGMISession> or PSGMIApiSession*
     * @param  session from Testcase
     * @param  testCaseName the name of Testcase
     * @return int
     */
    template <typename T>
    static int view(T session, const std::string& testCaseName)
    {
        std::cout << "View model in testcase: " << testCaseName << std::endl;

        static_assert(
            std::is_same<T, psgmstd::shared_ptr<PSGMISession>>::value || std::is_same<T, PSGMIApiSession*>::value,
            "Invalid session type. Only psgmstd::shared_ptr<PSGMISession> and PSGMIApiSession* are supported.");

        auto now = std::chrono::system_clock::now();
        auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count();

        std::string kernelPath = PS_KERNEL_PATH;
        std::cout << "kernelPath: " << PS_KERNEL_PATH << std::endl;

        std::string startTime = getCurrentTime(timestamp);
        int status = 1;
        std::string filePath = kernelPath + PS_CLOUD_MODEL_SAVE_PATH + "/" + testCaseName + "_" + startTime;

        std::cout << "Begin to save BrepJson File" << std::endl;
        int ret = writeBrepJsonToFile(session, filePath);

        if (ret != 0)
        {
            std::cout << "Failed to save BrepJson File" << std::endl;
            status = 0;
        }

        if (ret == 0)
        {
            std::cout << "Begin to save Facet File" << std::endl;
            ret = writeFacetToFile(session, filePath);

            if (ret != 0)
            {
                std::cout << "Failed to save Facet File" << std::endl;
                status = 0;
            }
        }

        std::string urlParam = testCaseName + "_" + startTime;

        if (status == 1)
        {
            std::cout << "Model View testCase: "
                      << "\033[34m" << urlParam << "\033[0m" << std::endl;
        }

#ifdef _WIN32
        now = std::chrono::system_clock::now();
        auto timestamp2 = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count();

        std::stringstream command;
        command << kernelPath << PS_CLOUD_NGINX_START_SCRIPT_PATH << " " << status << " " << urlParam << " "
                << testCaseName << " " << timestamp << " " << timestamp2;
        system(command.str().c_str());
#endif
        return 0;
    }

private:
    /**
     * @brief Write string data to a file
     * @param[in] data string
     * @param[in] savePath save path of the file
     * @return int
     */
    static int writeStringToFile(const std::string& data, const std::string& savePath)
    {
        std::ofstream file(savePath);

        if (!file.is_open())
        {
            std::cout << "Failed to Open file: " << savePath << std::endl;
            return -1;
        }

        file << data;

        if (!file.good())
        {
            std::cout << "Failed to Write file: " << savePath << std::endl;
            return -1;
        }

        file.close();

        return 0;
    }

    /**
     * @brief get current time for save files
     * @param[in] timestamp long long
     * @return string YearMouthDayHourMinSec
     */
    static std::string getCurrentTime(long long timestamp)
    {
        std::time_t time = timestamp / 1000;
        std::tm tm;

#ifdef _WIN32
        localtime_s(&tm, &time);
#else
        localtime_r(&time, &tm);
#endif

        std::stringstream ss;
        ss << std::put_time(&tm, "%Y%m%d%H%M%S");

        return ss.str();
    }

    /**
     * @brief Write brep json to a file
     * @param[in] session used session from testcase
     * @param[in] savePath save path of the file
     * @return int
     */
    static int writeBrepJsonToFile(PSGMIApiSession* session, const std::string& savePath)
    {
        static auto debug = PSGMIApiDebug::create();
        char* json = nullptr;
        PSGMApiErrorCode errCode = debug->getBrepJson(session, 0, &json);

        if (errCode != PSGMApiErrorCode::NoError)
        {
            std::cout << "Failed to getBrepJson, PSGMApiErrorCode: " << (int)errCode << std::endl;
            PSGMMemory::free(json);
            return -1;
        }

        auto ret = writeStringToFile(json, savePath + "_brep.json");
        PSGMMemory::free(json);
        return ret;
    }

    /**
     * @brief Write brep json to a file
     * @param[in] session used session from testcase
     * @param[in] savePath save path of the file
     * @return int
     */
    static int writeBrepJsonToFile(psgmstd::shared_ptr<PSGMISession> session, const std::string& savePath)
    {
        auto json = PSGMBrepJsonWriter::getJsonString(session->getContext().get(), false);
        if (json.empty())
        {
            return -1;
        }

        int ret = writeStringToFile(json.c_str(), savePath + "_brep.json");
        return ret;
    }

    /**
     * @brief Write facet json to a file
     * @param[in] visualization used to renderTopoFacet
     * @param[in] bodyTags tags of all bodies
     * @param[in] savePath save path of the file
     * @return int
     */
    static int writeFacetJsonToFile(PSGMIApiVisualization* visualization,
                                    const PSGMTags& bodyTags,
                                    const PSGMTags& geomTags,
                                    const std::string& savePath)
    {
        GoCallbacks callback;
        visualization->registerGraphicOutputCallbacks(&callback);

        PSGMApiErrorCode errCode = PSGMApiErrorCode::NoError;

        if (bodyTags.size() > 0)
        {
            std::cout << "Begin to renderTopoFacet" << std::endl;
            PSGMApiRenderTopoFacetOption options;
            errCode = visualization->renderTopoFacet(
                static_cast<int>(bodyTags.size()), bodyTags.data(), {}, PSGM_API_NULL_TAG, options);
        }

        if (errCode != PSGMApiErrorCode::NoError)
        {
            std::cout << "Failed to renderTopoFacet, PSGMApiErrorCode: " << (int)errCode << std::endl;
            return -1;
        }

        if (geomTags.size() > 0)
        {
            std::cout << "Begin to renderGeometry" << std::endl;
            PSGMApiRenderGeomOption geomOption;
            errCode = visualization->renderGeometry(static_cast<int>(geomTags.size()), geomTags.data(), {}, geomOption);
        }

        if (errCode != PSGMApiErrorCode::NoError)
        {
            std::cout << "Failed to renderGeometry, PSGMApiErrorCode: " << (int)errCode << std::endl;
            // return -1;
        }

        callback.proc();

        if (callback.getError() != 0)
        {
            std::cout << "Failed to renderTopoFacet, PSGMApiGOStatus from GoCallback: " << callback.getError()
                      << std::endl;
            return -1;
        }

        std::string facetJson = callback.getFacetJson();
        auto ret = writeStringToFile(facetJson, savePath + "_facet.json");

        if (ret != 0)
        {
            std::cout << "retfacetJson:" << ret << std::endl;
            return ret;
        }

        return ret;
    }

    /**
     * @brief Write topofacet data to a file
     * @param[in] session used session from testcase
     * @param[in] savePath save path of the file
     * @return int
     */
    static int writeFacetToFile(PSGMIApiSession* session, const std::string& savePath)
    {
        int partitionCount{0};
        PSGMApiTag* partitionTags{nullptr};

        PSGMApiErrorCode errCode = session->getPartitions(&partitionCount, &partitionTags);

        if (errCode != PSGMApiErrorCode::NoError)
        {
            PSGMMemory::free(partitionTags);
            std::cout << "Failed to getPartitions, PSGMApiErrorCode: " << (int)errCode << std::endl;
            return -1;
        }

        psgmstd::vector<PSGMApiTag> bodyTags;
        PSGMApiTag* partitionBodyTags{nullptr};

        psgmstd::vector<PSGMApiTag> geomTags;
        PSGMApiTag* partitionGeomTags{nullptr};

        for (int i = 0; i < partitionCount; i++)
        {
            PSGMApiTag partitionTag = partitionTags[i];

            int partitionBodyCount{0};
            errCode = session->partitionGetBodies(partitionTag, &partitionBodyCount, &partitionBodyTags);
            if (errCode != PSGMApiErrorCode::NoError)
            {
                PSGMMemory::free(partitionTags);
                PSGMMemory::free(partitionBodyTags);
                PSGMMemory::free(partitionGeomTags);
                std::cout << "Failed to partitionGetBodies, PSGMApiErrorCode: " << (int)errCode << std::endl;
                return -1;
            }

            int partitionGeomCount{0};
            errCode = session->partitionGetGeometries(partitionTag, &partitionGeomCount, &partitionGeomTags);
            if (errCode != PSGMApiErrorCode::NoError)
            {
                PSGMMemory::free(partitionTags);
                PSGMMemory::free(partitionBodyTags);
                PSGMMemory::free(partitionGeomTags);
                std::cout << "Failed to partitionGetGeometries, PSGMApiErrorCode: " << (int)errCode << std::endl;
                return -1;
            }

            for (int j = 0; j < partitionBodyCount; j++)
            {
                bodyTags.emplace_back(partitionBodyTags[j]);
            }

            for (int j = 0; j < partitionGeomCount; j++)
            {
                geomTags.emplace_back(partitionGeomTags[j]);
            }
        }
        PSGMMemory::free(partitionTags);
        PSGMMemory::free(partitionBodyTags);
        PSGMMemory::free(partitionGeomTags);

        PSGMIApiVisualization* visualization = session->getApiVisualization();
        return writeFacetJsonToFile(visualization, bodyTags, geomTags, savePath);
    }

    /**
     * @brief Write topofacet data to a file
     * @param[in] session used session from testcase
     * @param[in] savePath save path of the file
     * @return int
     */
    static int writeFacetToFile(psgmstd::shared_ptr<PSGMISession> session, const std::string& savePath)
    {
        PSGMTags bodyTags;
        PSGMTags geomTags;
        const auto& partitions = session->getPartitions();

        for (const auto& partition : partitions)
        {
            const auto& tags = partition->getWorld()->getBodies();
            bodyTags.insert(bodyTags.end(), tags.begin(), tags.end());

            const auto geoms = partition->getWorld()->getOrphanGeometries();
            geomTags.insert(geomTags.end(), geoms.begin(), geoms.end());
        }

        PSGMIApiVisualizationSPtr visualization = PSGMApiVisualization::createInstance(session->getContext().get(), 1);
        return writeFacetJsonToFile(visualization.get(), bodyTags, geomTags, savePath);
    }
};