#include <lange/Visualization.h>

RGB HSVtoRGB(float H, float S, float V)
{
    float C = S * V;
    float X = C * (1 - std::abs(fmod(H / 60.0, 2) - 1));
    float m = V - C;

    float r, g, b;

    if (H >= 0 && H < 60)
    {
        r = C, g = X, b = 0;
    }
    else if (H >= 60 && H < 120)
    {
        r = X, g = C, b = 0;
    }
    else if (H >= 120 && H < 180)
    {
        r = 0, g = C, b = X;
    }
    else if (H >= 180 && H < 240)
    {
        r = 0, g = X, b = C;
    }
    else if (H >= 240 && H < 300)
    {
        r = X, g = 0, b = C;
    }
    else
    {
        r = C, g = 0, b = X;
    }

    RGB rgb;
    rgb.r = (r + m) * 255;
    rgb.g = (g + m) * 255;
    rgb.b = (b + m) * 255;

    return rgb;
}

std::vector<RGB> generateColors(int n)
{
    std::vector<RGB> colors;
    for (int i = 0; i < n; i++)
    {
        float H = static_cast<float>(i) / static_cast<float>(n) * 360.0f;
        float S = 1.0f;
        float V = 1.0f;
        colors.push_back(HSVtoRGB(H, S, V));
    }
    return colors;
}

TQSI::Viewer::Viewer(igl::opengl::glfw::Viewer *v)
    : viewerPointer(v)
{
}

void TQSI::Viewer::addSurface(const Eigen::Matrix<float, 10, 1> &surface)
{
    this->surfaceList.push_back(surface);
}

void TQSI::Viewer::addSurface(const std::vector<Eigen::Matrix<float, 10, 1>> &surfaces)
{
    for (auto &surface : surfaces)
    {
        this->surfaceList.push_back(surface);
    }
}

void TQSI::Viewer::addCurve(const curve_equation &curve)
{
    this->curveList.push_back(curve);
}

void TQSI::Viewer::addCurve(const std::vector<curve_equation> &curves)
{
    for (auto &curve : curves)
    {
        this->curveList.push_back(curve);
    }
}

void TQSI::Viewer::addPoint(const Eigen::MatrixXd &point)
{
    this->pointList.push_back(point);
}

void TQSI::Viewer::addPoint(const std::vector<Eigen::MatrixXd> &points)
{
    for (auto &point : points)
    {
        this->pointList.push_back(point);
    }
}

void TQSI::Viewer::addCSG(const CSGtree csg)
{
    CSGList.push_back(csg);
}

void TQSI::Viewer::addCSG(const std::vector<CSGtree> &csgs)
{
    for (auto &csg : csgs)
    {
        CSGList.push_back(csg);
    }
}

void saveMeshToOff(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const std::string &filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }

    // Write header
    file << "OFF\n";
    file << V.rows() << " " << F.rows() << " 0\n";

    // Write vertices
    for (int i = 0; i < V.rows(); ++i)
    {
        file << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
    }

    // Write faces
    for (int i = 0; i < F.rows(); ++i)
    {
        file << "3 " << F(i, 0) << " " << F(i, 1) << " " << F(i, 2) << "\n";
    }

    file.close();
    std::cout << "Mesh saved to " << filename << std::endl;
}

void TQSI::Viewer::drawSurface()
{
    std::vector<Eigen::MatrixXd> S;

    Eigen::MatrixXd GV;
    Eigen::MatrixXd BBOX(2, 3);
    BBOX << -1, -1, -1, 1, 1, 1;
    BBOX *= 5;
    Eigen::RowVector3i res;
    igl::voxel_grid(BBOX, 0, 128, 0, GV, res);

    for (int i = 0; i < surfaceList.size(); i++)
    {
        Eigen::MatrixXd temp(GV.rows(), 1);
        for (int j = 0; j < GV.rows(); j++)
        {
            temp(j) = quadric_sdf_function(GV(j, 0), GV(j, 1), GV(j, 2), surfaceList[i]);
        }
        S.push_back(temp);
    }

    // marching_cubes
    std::vector<Eigen::MatrixXi> F;
    std::vector<Eigen::MatrixXd> V;

    for (int i = 0; i < S.size(); i++)
    {
        Eigen::MatrixXd tempv;
        Eigen::MatrixXi tempf;
        igl::marching_cubes(S[i], GV, res(0), res(1), res(2), 0, tempv, tempf);
        F.push_back(tempf);
        V.push_back(tempv);
    }

    Eigen::MatrixXd mesh_V;
    Eigen::MatrixXi mesh_F;

    for (int i = 0; i < F.size(); ++i)
    {
        mesh_F.conservativeResize(mesh_F.rows() + F[i].rows(), F[i].cols());
        mesh_F.bottomRows(F[i].rows()) = F[i].array() + mesh_V.rows();
        mesh_V.conservativeResize(mesh_V.rows() + V[i].rows(), V[i].cols());
        mesh_V.bottomRows(V[i].rows()) = V[i];
    }

    int x = mesh_F.rows();
    Eigen::MatrixXd Colors(x, 3);
    Eigen::RowVector3d color_red = Eigen::RowVector3d(1, 0, 0);
    Eigen::RowVector3d color_green = Eigen::RowVector3d(0, 1, 0);
    Eigen::RowVector3d color1 = Eigen::RowVector3d(1, 0, 0);
    Eigen::RowVector3d color2 = Eigen::RowVector3d(0, 1, 0);
    Colors << color1.replicate(F[0].rows(), 1),
        color2.replicate(F[1].rows(), 1);
    // Colors << color_red.replicate(x, 1);

    this->viewerPointer->data().set_mesh(mesh_V, mesh_F);
    this->viewerPointer->data().set_colors(Colors);
    // saveMeshToOff(mesh_V, mesh_F, "sphere.off");
}

void TQSI::Viewer::sampleCurvePoint(const curve_equation &eqution, Eigen::MatrixXd &mat, Eigen::MatrixXd &color, const double step)
{
    // 假设我们有一个 std::vector 存储 Eigen::Vector3d 对象
    std::vector<Eigen::Vector3d> vec;
    std::vector<Eigen::Vector3d> colorVec;

    auto u = eqution.u;
    auto Delta = eqution.Delta;
    for (double i = -M_PI; i < M_PI; i += step)
    {
        double j = tan(i / 2);
        double tx = 0.0;
        double ty = 0.0;
        double tz = 0.0;
        double tw = 0.0;
        if (eqution.have_delta)
        {
            double delta = GiNaC::ex_to<GiNaC::numeric>(eqution.delta_equation.subs(u == j)).to_double();
            tx = GiNaC::ex_to<GiNaC::numeric>(eqution.x_equation.subs(GiNaC::lst{u == j, Delta == delta})).to_double();
            ty = GiNaC::ex_to<GiNaC::numeric>(eqution.y_equation.subs(GiNaC::lst{u == j, Delta == delta})).to_double();
            tz = GiNaC::ex_to<GiNaC::numeric>(eqution.z_equation.subs(GiNaC::lst{u == j, Delta == delta})).to_double();
            tw = GiNaC::ex_to<GiNaC::numeric>(eqution.w_equation.subs(GiNaC::lst{u == j, Delta == delta})).to_double();
        }
        else
        {
            tx = GiNaC::ex_to<GiNaC::numeric>(eqution.x_equation.subs(u == j)).to_double();
            ty = GiNaC::ex_to<GiNaC::numeric>(eqution.y_equation.subs(u == j)).to_double();
            tz = GiNaC::ex_to<GiNaC::numeric>(eqution.z_equation.subs(u == j)).to_double();
            tw = GiNaC::ex_to<GiNaC::numeric>(eqution.w_equation.subs(u == j)).to_double();
        }
        if (tw == 0)
        {
            tx = 0.0;
            ty = 0.0;
            tz = 0.0;
        }
        else
        {
            tx /= tw;
            ty /= tw;
            tz /= tw;
        }

        vec.push_back(Eigen::Vector3d(tx, ty, tz));
        colorVec.push_back(Eigen::Vector3d(0.0, 0.0, 1.0));
    }

    // 转换 std::vector<Eigen::Vector3d> 为 Eigen::MatrixXd
    mat.resize(vec.size(), 3);   // 初始化大小为 vec.size() x 3 的矩阵
    color.resize(vec.size(), 3); // 初始化大小为 vec.size() x 3 的矩阵
    for (size_t i = 0; i < vec.size(); ++i)
    {
        mat.row(i) = vec[i];
        color.row(i) = colorVec[i];
    }
    // std::cout << "Matrix converted from std::vector:\n" << mat << std::endl;
}

void TQSI::Viewer::sampleCurvePoint(const curve_equation &eqution, Eigen::MatrixXd &mat, const double step)
{
    // 假设我们有一个 std::vector 存储 Eigen::Vector3d 对象
    std::vector<Eigen::Vector3d> vec;

    auto u = eqution.u;
    auto Delta = eqution.Delta;
    for (double i = -M_PI; i < M_PI; i += step)
    {
        double j = tan(i / 2);
        double tx = 0.0;
        double ty = 0.0;
        double tz = 0.0;
        double tw = 0.0;
        if (eqution.have_delta)
        {
            double delta = GiNaC::ex_to<GiNaC::numeric>(eqution.delta_equation.subs(u == j)).to_double();
            tx = GiNaC::ex_to<GiNaC::numeric>(eqution.x_equation.subs(GiNaC::lst{u == j, Delta == delta})).to_double();
            ty = GiNaC::ex_to<GiNaC::numeric>(eqution.y_equation.subs(GiNaC::lst{u == j, Delta == delta})).to_double();
            tz = GiNaC::ex_to<GiNaC::numeric>(eqution.z_equation.subs(GiNaC::lst{u == j, Delta == delta})).to_double();
            tw = GiNaC::ex_to<GiNaC::numeric>(eqution.w_equation.subs(GiNaC::lst{u == j, Delta == delta})).to_double();
        }
        else
        {
            tx = GiNaC::ex_to<GiNaC::numeric>(eqution.x_equation.subs(u == j)).to_double();
            ty = GiNaC::ex_to<GiNaC::numeric>(eqution.y_equation.subs(u == j)).to_double();
            tz = GiNaC::ex_to<GiNaC::numeric>(eqution.z_equation.subs(u == j)).to_double();
            tw = GiNaC::ex_to<GiNaC::numeric>(eqution.w_equation.subs(u == j)).to_double();
        }
        if (tw == 0)
        {
            tx = 0.0;
            ty = 0.0;
            tz = 0.0;
        }
        else
        {
            tx /= tw;
            ty /= tw;
            tz /= tw;
        }

        vec.push_back(Eigen::Vector3d(tx, ty, tz));
    }

    // 转换 std::vector<Eigen::Vector3d> 为 Eigen::MatrixXd
    mat.resize(vec.size(), 3); // 初始化大小为 vec.size() x 3 的矩阵
    for (size_t i = 0; i < vec.size(); ++i)
    {
        mat.row(i) = vec[i];
    }
    // std::cout << "Matrix converted from std::vector:\n" << mat << std::endl;
}

bool TQSI::Viewer::isParameterInValidSegement(const curve_equation &eqution, const double parameter)
{
    for (auto segment : eqution.validSegment)
    {
        if (segment.contains(parameter))
        {
            return true;
        }
    }
    return false;
}

void TQSI::Viewer::sampleCurvePointWithOneInterval(const curve_equation &eqution, Eigen::MatrixXd &mat, const double step)
{
    std::vector<Eigen::Vector3d> vec;

    auto u = eqution.u;
    auto Delta = eqution.Delta;
    // std::vector<double> poiii={-31.8205, -15.8945, -10.5789, -7.91582, -6.31375, -5.24218, -4.59512, -3.89182, -3.4761, -3.17805, -2.94444, -2.75154, -2.58669, -2.44235, -2.31363, -2.19722, -2.09074, -1.99243, -1.90096, -1.81529, -1.7346, -1.65823, -1.58563, -1.51635, -1.45001, -1.38629, -1.32493, -1.26567, -1.20831, -1.15268, -1.09861, -1.04597, -0.994623, -0.944462, -0.895384, -0.847298, -0.800119, -0.753772, -0.708185, -0.663294, -0.619039, -0.575364, -0.532217, -0.489548, -0.447312, -0.405465, -0.363965, -0.322773, -0.281851, -0.241162, -0.200671, -0.160343, -0.120144, -0.0800427, -0.0400053, 8.88178e-16, 0.0400053, 0.0800427, 0.120144, 0.160343, 0.200671, 0.241162, 0.281851, 0.322773, 0.363965, 0.405465, 0.447312, 0.489548, 0.532217, 0.575364, 0.619039, 0.663294, 0.708185, 0.753772, 0.800119, 0.847298, 0.895384, 0.944462, 0.994623, 1.04597, 1.09861, 1.15268, 1.20831, 1.26567, 1.32493, 1.38629, 1.45001, 1.51635, 1.58563, 1.65823, 1.7346, 1.81529, 1.90096, 1.99243, 2.09074, 2.19722, 2.31363, 2.44235, 2.58669, 2.75154, 2.94444, 3.17805, 3.4761, 3.89182, 4.59512,  5.24218,  6.31375,  7.91582,  10.5789,  15.8945,  31.8205}; 
    std::vector<double> poiii={-31.8205, -15.8945, -6.91582,-4.59512, -3.89182, -3.4761, -3.17805, -2.94444, -2.75154, -2.58669, -2.44235, -2.31363, -2.19722, -2.09074, -1.99243, -1.90096, -1.81529, -1.7346, -1.65823, -1.58563, -1.51635, -1.45001, -1.38629, -1.32493, -1.26567, -1.20831, -1.15268, -1.09861, -1.04597, -0.994623, -0.944462, -0.895384, -0.847298, -0.800119, -0.753772, -0.708185, -0.663294, -0.619039, -0.575364, -0.532217, -0.489548, -0.447312, -0.405465, -0.363965, -0.322773, -0.281851, -0.241162, -0.200671, -0.160343, -0.120144, -0.0800427, -0.0400053, 8.88178e-16, 0.0400053, 0.0800427, 0.120144, 0.160343, 0.200671, 0.241162, 0.281851, 0.322773, 0.363965, 0.405465, 0.447312, 0.489548, 0.532217, 0.575364, 0.619039, 0.663294, 0.708185, 0.753772, 0.800119, 0.847298, 0.895384, 0.944462, 0.994623, 1.04597, 1.09861, 1.15268, 1.20831, 1.26567, 1.32493, 1.38629, 1.45001, 1.51635, 1.58563, 1.65823, 1.7346, 1.81529, 1.90096, 1.99243, 2.09074, 2.19722, 2.31363, 2.44235, 2.58669, 2.75154, 2.94444, 3.17805, 3.4761, 3.89182, 4.59512,      6.91582,  15.8945,  31.8205,};
    std::vector<double> poiii2={-31.8205,-15.8945,-6.91582,-4.59512,-3.89182,-3.4761,-3.17805,-2.84799,-2.51452,-2.25542,-2.04158,-1.81695,-1.58674,-1.38708,-1.20889,-1.0464,-0.871816,-0.686072,-0.489981,-0.282084,-0.080107,0,0.080107,0.282084,0.489981,0.686072,0.871816,1.0464,1.20889,1.38708,1.58674,1.81695,2.04158,2.25542,2.51452,2.84799,3.17805,3.4761,3.89182,4.59512,6.91582,15.8945,31.8205};
    PRINT_VAR(poiii.size());
    std::vector<double> linshi;
    int cnt=0;
    // for (double i = -M_PI; i < M_PI; i += step)
    // for (double i = step; i < 1; i += step)
    // for (double i = -30; i <= 30; i += step)
    // for(auto i:poiii)
    // for(int i=0;i<poiii.size();i++)
    for(int i=0;i<poiii2.size();i++)
    {
        // double j = tan(i / 2);
        // double j=std::log(i/(1-i));
        // double j=i;
        // double j=poiii[i];
        double j=poiii2[i];
        // if(abs(j)<=0.75)
        // {
        //     j+=poiii[++i];
        //     j+=poiii[++i];
        //     j+=poiii[++i];
        //     j+=poiii[++i];
        //     j/=5;
        // }
        // else if(abs(j)<=1)
        // {
        //     j+=poiii[++i];
        //     j+=poiii[++i];
        //     j+=poiii[++i];
        //     j/=4;
        // }
        // else if(abs(j)<=2)
        // {
        //     // j+=0.125*step;
        //     j+=poiii[++i];
        //     j+=poiii[++i];
        //     j/=3;
        // }
        // else if(abs(j)<=3)
        // {
        //     j+=poiii[++i];
        //     j/=2;
        // }
        std::cout<<cnt++<<"-th ="<<j<<endl;
        
        // if(eqution.isAllInterval)
        // {

        // }
        // else if (!isParameterInValidSegement(eqution, j))
        //     continue;
        // if(!(j<eqution.startParameter||j>eqution.endParameter))continue;
        double tx = 0.0;
        double ty = 0.0;
        double tz = 0.0;
        double tw = 0.0;
        if (eqution.have_delta)
        {
            double delta = GiNaC::ex_to<GiNaC::numeric>(eqution.delta_equation.subs(u == j).evalf()).to_double();
            tx = GiNaC::ex_to<GiNaC::numeric>(eqution.x_equation.subs(GiNaC::lst{u == j, Delta == delta}).evalf()).to_double();
            ty = GiNaC::ex_to<GiNaC::numeric>(eqution.y_equation.subs(GiNaC::lst{u == j, Delta == delta}).evalf()).to_double();
            tz = GiNaC::ex_to<GiNaC::numeric>(eqution.z_equation.subs(GiNaC::lst{u == j, Delta == delta}).evalf()).to_double();
            tw = GiNaC::ex_to<GiNaC::numeric>(eqution.w_equation.subs(GiNaC::lst{u == j, Delta == delta}).evalf()).to_double();
        }
        else
        {
            tx = GiNaC::ex_to<GiNaC::numeric>(eqution.x_equation.subs(u == j).evalf()).to_double();
            ty = GiNaC::ex_to<GiNaC::numeric>(eqution.y_equation.subs(u == j).evalf()).to_double();
            tz = GiNaC::ex_to<GiNaC::numeric>(eqution.z_equation.subs(u == j).evalf()).to_double();
            tw = GiNaC::ex_to<GiNaC::numeric>(eqution.w_equation.subs(u == j).evalf()).to_double();
        }
        if (tw == 0)
        {
            tx = 0.0;
            ty = 0.0;
            tz = 0.0;
        }
        else
        {
            tx /= tw;
            ty /= tw;
            tz /= tw;
        }

        vec.push_back(Eigen::Vector3d(tx, ty, tz));
    }

    for(int i=linshi.size()-1;i>=0;i--)cout<<linshi[i]<<",";
    // 转换 std::vector<Eigen::Vector3d> 为 Eigen::MatrixXd
    mat.resize(vec.size(), 3); // 初始化大小为 vec.size() x 3 的矩阵
    for (size_t i = 0; i < vec.size(); ++i)
    {
        mat.row(i) = vec[i];
        // cout << "{" << vec[i][0] << ", " << vec[i][1] << ",  " << vec[i][2] << "}" << ",";
    }
}

void TQSI::Viewer::sampleCurvePointWithIntervals(curve_equation &eqution, std::vector<Eigen::MatrixXd> &mat, const double point_cnt)
{
    if(!eqution.validSegment.size())return;
    
    // double i;
    // double curSegment=0;
    // for (i = -M_PI; i < M_PI; i += step)
    // {
    //     double j = tan(i / 2);
    //     if(eqution.validSegment[0].contains(j))break;
    // }
    // std::vector<Eigen::Vector3d> vec;
    
    // while(i<M_PI)
    // {
    //     double j = tan(i / 2);
    //     if(curSegment==eqution.validSegment.size())break;

    //     if(!eqution.validSegment[curSegment].contains(j))
    //     {
    //         Eigen::MatrixXd temp;
    //         temp.resize(vec.size()+1, 3);

    //         auto sp=eqution.validSegment[curSegment].startParameter;
    //         auto startPoint=-INF_PARAMETER?(eqution.negativeInf):(eqution.getPoint(sp,1,1));
    //         temp.row(0)<<startPoint.x,startPoint.y,startPoint.z;
            
    //         for (size_t i = 0; i < vec.size(); ++i)
    //         {
    //             temp.row(i+1) = vec[i];
    //         }
    //         // auto ep=eqution.validSegment[curSegment].endParameter;
    //         // auto endPoint=INF_PARAMETER?(eqution.positiveInf):(eqution.getPoint(ep,1,1));
    //         curSegment++;
    //         // temp.row(vec.size()+1)<<endPoint.x,endPoint.y,endPoint.z;
    //         mat.push_back(temp);
    //         vec.clear();
    //     }
    //     else
    //     {
    //         TQSI::point samplePoint = eqution.getPoint(j, true, true);
    //         vec.push_back(Eigen::Vector3d(samplePoint.x, samplePoint.y, samplePoint.z));
    //         i+=step;
    //     }
    // }

    for (auto &segment : eqution.validSegment)
    {
        // cout<<"segment.startParameter= "<<segment.startParameter<<",segment.endParameter= "<<segment.endParameter<<endl;
        std::vector<Eigen::Vector3d> vec;
        Eigen::MatrixXd temp;
        double step=(segment.endParameter-segment.startParameter)/point_cnt;
        for (double i = segment.startParameter; i < segment.endParameter; i += step)
        {
            TQSI::point samplePoint = eqution.getPoint(i, true, true);
            vec.push_back(Eigen::Vector3d(samplePoint.x, samplePoint.y, samplePoint.z));
        }
        TQSI::point samplePoint = eqution.getPoint(segment.endParameter, true, true);
        vec.push_back(Eigen::Vector3d(samplePoint.x, samplePoint.y, samplePoint.z));
        temp.resize(vec.size(), 3);
        for (size_t i = 0; i < vec.size(); ++i)
        {
            temp.row(i) = vec[i];
        }
        mat.push_back(temp);
    }
}

void TQSI::Viewer::drawCurve()
{
    std::vector<Eigen::MatrixXd> mat;
    auto colors=generateColors(6);
    int curColor=1;

    for (auto &curve : curveList)
    {
        int point_cnt = curve.visualizaitonPoint_cnt;
        double step = 2 * M_PI / point_cnt;
        // sampleCurvePoint(curve,mat,step);
        
        // sampleCurvePointWithIntervals(curve, mat, step);
        // for (auto &m : mat)
        // {
        //     for (int i = 0; i < m.rows() - 1; i++)
        //     {
        //         viewerPointer->data().add_edges(
        //             m.row(i),
        //             m.row(i + 1),
        //             // Eigen::RowVector3d(0, 0, 1)
        //             Eigen::RowVector3d(colors[curColor].r,colors[curColor].g,colors[curColor].b)
        //             );
        //     }
        // }

        mat.push_back(Eigen::MatrixXd());
        sampleCurvePointWithOneInterval(curve,mat[0],step);
        for(int i=0;i<mat[0].rows()-1;i++)
        {
            viewerPointer->data().add_edges
            (
            mat[0].row(i),
            mat[0].row(i+1),
            // Eigen::RowVector3d(0,0,1)
            Eigen::RowVector3d(colors[curColor].r,colors[curColor].g,colors[curColor].b)
            );
        }

        curColor++;

        // if(curve.isLoop)
        // {
        //     viewerPointer->data().add_edges
        //     (
        //     mat.row(mat.rows()-1),
        //     mat.row(0),
        //     Eigen::RowVector3d(0,0,1)
        //     );
        // }
    }
}

void TQSI::Viewer::drawCurveByPoint()
{
    Eigen::MatrixXd mat;
    pointList.clear();
    for (auto &curve : curveList)
    {
        cout << "new curve" << endl;
        int point_cnt = curve.visualizaitonPoint_cnt;
        double step = 2 * M_PI / point_cnt;
        sampleCurvePoint(curve, mat, step);
        int n = mat.rows();
        for (int i = 0; i < n; i++)
        {
            pointList.push_back(mat.row(i));
            // viewerPointer->data().add_points(mat.row(i), Eigen::RowVector3d(0.0, 1.0*i/n, 0));
            viewerPointer->data().add_points(mat.row(i), Eigen::RowVector3d(0.0, 0, 1.0));
        }
        auto p0 = pointList[0];
        auto p1 = pointList.back();
        auto p2 = p1 - p0;
        auto normm = p2.squaredNorm();
        auto tan1 = pointList[1] - p0;
        auto tan2 = pointList[pointList.size() - 2] - p1;
        double dianji = tan1(0) * tan2(0) + tan1(1) * tan2(1) + tan1(2) * tan2(2);
        PRINT_VAR(p0);
        PRINT_VAR(p1);
        PRINT_VAR(p2);
        PRINT_VAR(normm);
        PRINT_VAR(tan1);
        PRINT_VAR(tan2);
        PRINT_VAR(dianji);
        pointList.clear();
    }
}

void TQSI::Viewer::drawPoint()
{
    for (auto &point : pointList)
    {
        viewerPointer->data().add_points(point, Eigen::RowVector3d(0.0, 1.0, 0.0));
    }
}

void TQSI::Viewer::drawCSG()
{
    std::vector<Eigen::MatrixXd> S;

    Eigen::MatrixXd GV;
    Eigen::MatrixXd BBOX(2, 3);
    BBOX << -1, -1, -1, 1, 1, 1;
    BBOX *= 5;
    Eigen::RowVector3i res;
    igl::voxel_grid(BBOX, 0, 128, 0, GV, res);

    for (int i = 0; i < CSGList.size(); i++)
    {
        Eigen::MatrixXd temp(GV.rows(), 1);
        for (int j = 0; j < GV.rows(); j++)
        {
            temp(j) = CSGList[i].get_sdf(GV(j, 0), GV(j, 1), GV(j, 2));
        }
        S.push_back(temp);
    }

    // marching_cubes
    std::vector<Eigen::MatrixXi> F;
    std::vector<Eigen::MatrixXd> V;

    for (int i = 0; i < S.size(); i++)
    {
        Eigen::MatrixXd tempv;
        Eigen::MatrixXi tempf;
        igl::marching_cubes(S[i], GV, res(0), res(1), res(2), 0, tempv, tempf);
        F.push_back(tempf);
        V.push_back(tempv);
    }

    Eigen::MatrixXd mesh_V;
    Eigen::MatrixXi mesh_F;

    for (int i = 0; i < F.size(); ++i)
    {
        mesh_F.conservativeResize(mesh_F.rows() + F[i].rows(), F[i].cols());
        mesh_F.bottomRows(F[i].rows()) = F[i].array() + mesh_V.rows();
        mesh_V.conservativeResize(mesh_V.rows() + V[i].rows(), V[i].cols());
        mesh_V.bottomRows(V[i].rows()) = V[i];
    }

    int x = mesh_F.rows();
    Eigen::MatrixXd Colors(x, 3);
    Eigen::RowVector3d color_red = Eigen::RowVector3d(1, 0, 0);
    Eigen::RowVector3d color_green = Eigen::RowVector3d(0, 1, 0);
    // Eigen::RowVector3d color1 = Eigen::RowVector3d(1, 0, 0);
    // Eigen::RowVector3d color2 = Eigen::RowVector3d(0, 1, 0);
    // Colors << color1.replicate(F[0].rows(), 1),
    //     color2.replicate(F[1].rows(), 1);
    Colors << color_red.replicate(x, 1);

    this->viewerPointer->data().set_mesh(mesh_V, mesh_F);
    this->viewerPointer->data().set_colors(Colors);
    // saveMeshToOff(mesh_V, mesh_F, "sphere.off");
}

void TQSI::Viewer::visualizate()
{
    viewerPointer->data().double_sided = true;
    viewerPointer->data().show_lines = false;
    int colorCount = 0;
    if (visualizateSurface)
        colorCount += surfaceList.size();
    if (visualizateCurve)
        colorCount += curveList.size();
    if (visualizatePoint)
        colorCount += pointList.size();
    if (visualizateCSG)
        colorCount += CSGList.size();

    if (!colorCount)
        return;

    auto colors = generateColors(colorCount);

    if (visualizateSurface)
    {
        drawSurface();
    }

    if (visualizateCurve)
    {
        // curveList[0].startParameter=0;
        // curveList[0].endParameter=1000;
        // curveList[1].startParameter=0.333333;
        // curveList[1].endParameter=1000;
        // curveList[2].startParameter=0.333333;
        // curveList[2].endParameter=1000;
        drawCurve();
    }

    if (visualizatePoint)
    {
        drawPoint();
    }

    if (visualizateCSG)
    {
        drawCSG();
    }

    viewerPointer->launch();
}

void TQSI::Viewer::setVisualizateSurface(bool b)
{
    visualizateSurface = b;
}

void TQSI::Viewer::setVisualizateCurve(bool b)
{
    visualizateCurve = b;
}

void TQSI::Viewer::setVisualizatePoint(bool b)
{
    visualizatePoint = b;
}

void TQSI::Viewer::setVisualizateCSG(bool b)
{
    visualizateCSG = b;
}

void TQSI::Viewer::forDeBug()
{
    viewerPointer->data().double_sided = true;
    viewerPointer->data().show_lines = false;
    // viewerPointer->data().show_custom_labels = true;
    drawSurface();
    std::vector<Eigen::MatrixXd> mat;
    auto colors=generateColors(6);
    int curColor=1;
    int aa=0;

    for (auto &curve : curveList)
    {
        // if(!aa)
        // {
        //     aa++;
        //     continue;
        // }
        int point_cnt = curve.visualizaitonPoint_cnt;
        double step = 2 * M_PI / point_cnt;

        mat.push_back(Eigen::MatrixXd());
        
        sampleCurvePointWithIntervals(curve,mat,point_cnt);
        for(int i=0;i<mat.size();i++)
        {
            // viewerPointer->data().add_points(mat[i],Eigen::RowVector3d(colors[curColor].r,colors[curColor].g,colors[curColor].b));  
            for(int j=0;j<mat[i].rows()-1;j++)
            {
                viewerPointer->data().add_edges
                (
                mat[i].row(j),
                mat[i].row(j+1),
                // Eigen::RowVector3d(0,0,1)
                Eigen::RowVector3d(colors[curColor].r,colors[curColor].g,colors[curColor].b)
                );  
            }
        }

        // step=1.0 / 100;
        // #include <chrono>
        // auto start = std::chrono::high_resolution_clock::now();
        // sampleCurvePointWithOneInterval(curve,mat[0],step);
        //  // 获取程序结束运行时的时间点
        // auto end = std::chrono::high_resolution_clock::now();

        // // 计算运行时间，单位为毫秒
        // std::chrono::duration<double, std::milli> elapsed = end - start;

        // // 打印运行时间
        // std::cout << "Elapsed time: " << elapsed.count() << " ms" << std::endl;

        // for(int i=0;i<mat[0].rows();i++)
        // {
        //     viewerPointer->data().add_points(mat[0],Eigen::RowVector3d(colors[curColor].r,colors[curColor].g,colors[curColor].b)); 
        //     // viewerPointer->data().add_edges
        //     // (
        //     // mat[0].row(i),
        //     // mat[0].row(i+1),
        //     // // Eigen::RowVector3d(0,0,1)
        //     // Eigen::RowVector3d(colors[curColor].r,colors[curColor].g,colors[curColor].b)
        //     // );
        //     // viewerPointer->data().add_label(mat[0].row(i),to_string(i));
        // }
        
        curColor++;

        // if(curve.isLoop)
        // {
        //     viewerPointer->data().add_edges
        //     (
        //     mat.row(mat.rows()-1),
        //     mat.row(0),
        //     Eigen::RowVector3d(0,0,1)
        //     );
        // }
    }
    viewerPointer->launch();
}







