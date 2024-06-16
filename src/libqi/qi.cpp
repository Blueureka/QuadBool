#include <libqi/qi.h>

bigint_matrix vector2quadric(const Eigen::Matrix<float, 10, 1> qvec)
{

    bigint_matrix qmat = bigint_matrix(4, 4);

    if ((long(qvec(1)) % 2 == 0) && (long(qvec(2)) % 2 == 0) && (long(qvec(3)) % 2 == 0) && (long(qvec(5)) % 2 == 0) &&
        (long(qvec(6)) % 2 == 0) && (long(qvec(8)) % 2 == 0))
    {
        qmat.sto(0, 0, bigint(qvec(0)));
        qmat.sto(0, 1, bigint(qvec(1)) / 2);
        qmat.sto(1, 0, bigint(qvec(1)) / 2);
        qmat.sto(0, 2, bigint(qvec(2)) / 2);
        qmat.sto(2, 0, bigint(qvec(2)) / 2);
        qmat.sto(0, 3, bigint(qvec(3)) / 2);
        qmat.sto(3, 0, bigint(qvec(3)) / 2);
        qmat.sto(1, 1, bigint(qvec(4)));
        qmat.sto(1, 2, bigint(qvec(5)) / 2);
        qmat.sto(2, 1, bigint(qvec(5)) / 2);
        qmat.sto(1, 3, bigint(qvec(6)) / 2);
        qmat.sto(3, 1, bigint(qvec(6)) / 2);
        qmat.sto(2, 2, bigint(qvec(7)));
        qmat.sto(2, 3, bigint(qvec(8)) / 2);
        qmat.sto(3, 2, bigint(qvec(8)) / 2);
        qmat.sto(3, 3, bigint(qvec(9)));
    }
    else
    {
        qmat.sto(0, 0, 2 * bigint(qvec(0)));
        qmat.sto(0, 1, bigint(qvec(1)));
        qmat.sto(1, 0, bigint(qvec(1)));
        qmat.sto(0, 2, bigint(qvec(2)));
        qmat.sto(2, 0, bigint(qvec(2)));
        qmat.sto(0, 3, bigint(qvec(3)));
        qmat.sto(3, 0, bigint(qvec(3)));
        qmat.sto(1, 1, 2 * bigint(qvec(4)));
        qmat.sto(1, 2, bigint(qvec(5)));
        qmat.sto(2, 1, bigint(qvec(5)));
        qmat.sto(1, 3, bigint(qvec(6)));
        qmat.sto(3, 1, bigint(qvec(6)));
        qmat.sto(2, 2, 2 * bigint(qvec(7)));
        qmat.sto(2, 3, bigint(qvec(8)));
        qmat.sto(3, 2, bigint(qvec(8)));
        qmat.sto(3, 3, 2 * bigint(qvec(9)));
    }

    return qmat;
}

void vcf2vcbigint(Eigen::Matrix<float, 10, 1> &qvec, short mask_length)
{
    // int mask = 0xFFFFF000;
    unsigned int mask = 0xFFFFFFFF;
    mask_length = (mask_length < 22 ? mask_length : 22);

    mask = mask << mask_length;
    // cout<<"mask="<<mask<<endl;

    // unsigned char* ptr = reinterpret_cast<unsigned char*>(&mask);
    //   for (int bnvnv = sizeof(int) - 1; bnvnv >= 0; bnvnv--)
    //     {
    //       for (int j = 7; j >= 0; j--)
    //       {
    //         std::cout << ((ptr[bnvnv] >> j) & 1);
    //       }
    //     }
    //     cout<<endl;

    unsigned int exp_mask = 0x7F800000;
    int exp_mask_rshift = 0;

    int min_exp = 200;

    float farr[10];
    for (int i = 0; i < 10; i++)
    {
        farr[i] = qvec(i);
    }

    for (int i = 0; i < 10; i++)
    {
        // void* fadd = &qvec(i);

        // void* fadd = &farr[i];
        // int *ifa = (int*)fadd;

        // 后续根据精度需求调整mask长度，以及容忍左移的次数
        unsigned int x_bits = *reinterpret_cast<unsigned int *>(&farr[i]); // 将浮点数x转换为unsigned int型整数
        unsigned int y_bits = mask & x_bits;
        unsigned int exp_bits = (x_bits & 0x7F800000) >> 23; // 获取浮点数的指数部分
        int exp = static_cast<int>(exp_bits) - 127;          // 计算浮点数的指数值
        float lkjmcnv = *reinterpret_cast<float *>(&y_bits);
        farr[i] = lkjmcnv;
        min_exp = min(min_exp, exp);

        // get exp
        //  int *fexp=(int*)fadd;
        //  (*fexp)=(*fexp)&exp_mask;
        //  max_exp=max((((*fexp)>>23)-1),max_exp);

        // qvec(i)=lkjmcnv;

        // print float in binary style
        //  unsigned char* ptr = reinterpret_cast<unsigned char*>(&farr[i]);
        //  for (int bnvnv = sizeof(float) - 1; bnvnv >= 0; bnvnv--)
        //  {
        //    for (int j = 7; j >= 0; j--)
        //    {
        //      std::cout << ((ptr[bnvnv] >> j) & 1);
        //    }
        //  }
        //  cout<<endl;
    }

    for (int i = 0; i < 10; i++)
    {
        qvec(i) = farr[i];
    }

    // cout<<qvec<<endl;

    // exp_mask_rshift=min(127-max_exp,max_rshift);
    // for(int i=0;i<exp_mask_rshift;i++)
    // {
    //       for(int j=0;j<10;j++)
    //       {
    //             if((qvec[j]) - (int)(qvec[j]))
    //             {
    //                   qvec=2*qvec;
    //                   break;
    //             }
    //       }
    // }
    int max_rshift = 23 - mask_length;
    int cnt_shift = max_rshift - min_exp;
    for (int i = 0; i < cnt_shift; i++)
    {
        bool all_int = true;
        for (int j = 0; j < 10; j++)
        {
            float int_part, frac_part;
            frac_part = std::modf(qvec[j], &int_part); // 获取浮点数的整数部分和小数部分
            if (frac_part != 0)
            {
                all_int = false;
                break;
            }
        }
        if (!all_int)
        {
            qvec = 2 * qvec;
        }
        else
        {
            break;
        }
    }
}

mpz_class bigint_gcd(mpz_class a, mpz_class b)
{
    if (b != 0)
    {
        return bigint_gcd(b, a % b);
    }
    else
        return a;
}

bigint_matrix vcf2gmp2bigintmatrix(Eigen::Matrix<float, 10, 1> &qvec)
{
    int min_zhishu = 256;
    int cur_zhishu = 0;
    for (int i = 0; i < 10; i++)
    {
        frexp(qvec(i), &cur_zhishu);
        min_zhishu = min(min_zhishu, cur_zhishu);
    }
    int cnt_shift = 23;
    if (min_zhishu < 0)
    {
        cnt_shift += (-min_zhishu);
    }

    bigint_matrix qmat = bigint_matrix(4, 4);

    math_vector<bigint> bigint_qvec;
    bigint_qvec.set_size(10);
    for (int i = 0; i < 10; i++)
    {
        mpf_class linshi = qvec(i);
        linshi = linshi << cnt_shift;
        bigint linshiint;
        linshiint.val = (mpz_class)linshi;
        bigint_qvec[i].assign(linshiint);
    }

    // gcd here
    mpz_class fenmu = bigint_gcd(bigint_qvec[0].val, bigint_qvec[1].val);
    for (int i = 2; i < 10; i++)
    {
        /* code */
        fenmu = bigint_gcd(fenmu, bigint_qvec[i].val);
    }
    if (fenmu < 0)
        fenmu = -fenmu;

    for (int i = 0; i < 10; i++)
    {
        // bigint_qvec[i].val = bigint_qvec[i].val / fenmu 
        bigint_qvec[i].val = bigint_qvec[i].val / fenmu *2;
        // cout<<i<<"="<<bigint_qvec[i].val<<endl;
    }
    // gcd done

    if ((bigint_qvec[1].is_even()) && (bigint_qvec[2].is_even()) && (bigint_qvec[3].is_even()) &&
        (bigint_qvec[5].is_even()) && (bigint_qvec[6].is_even()) && (bigint_qvec[8].is_even()))
    {
        qmat.sto(0, 0, bigint_qvec[0]);
        qmat.sto(0, 1, bigint_qvec[1] / 2);
        qmat.sto(1, 0, bigint_qvec[1] / 2);
        qmat.sto(0, 2, bigint_qvec[2] / 2);
        qmat.sto(2, 0, bigint_qvec[2] / 2);
        qmat.sto(0, 3, bigint_qvec[3] / 2);
        qmat.sto(3, 0, bigint_qvec[3] / 2);
        qmat.sto(1, 1, bigint_qvec[4]);
        qmat.sto(1, 2, bigint_qvec[5] / 2);
        qmat.sto(2, 1, bigint_qvec[5] / 2);
        qmat.sto(1, 3, bigint_qvec[6] / 2);
        qmat.sto(3, 1, bigint_qvec[6] / 2);
        qmat.sto(2, 2, bigint_qvec[7]);
        qmat.sto(2, 3, bigint_qvec[8] / 2);
        qmat.sto(3, 2, bigint_qvec[8] / 2);
        qmat.sto(3, 3, bigint_qvec[9]);
    }
    else
    {
        qmat.sto(0, 0, 2 * bigint_qvec[0]);
        qmat.sto(0, 1, bigint_qvec[1]);
        qmat.sto(1, 0, bigint_qvec[1]);
        qmat.sto(0, 2, bigint_qvec[2]);
        qmat.sto(2, 0, bigint_qvec[2]);
        qmat.sto(0, 3, bigint_qvec[3]);
        qmat.sto(3, 0, bigint_qvec[3]);
        qmat.sto(1, 1, 2 * bigint_qvec[4]);
        qmat.sto(1, 2, bigint_qvec[5]);
        qmat.sto(2, 1, bigint_qvec[5]);
        qmat.sto(1, 3, bigint_qvec[6]);
        qmat.sto(3, 1, bigint_qvec[6]);
        qmat.sto(2, 2, 2 * bigint_qvec[7]);
        qmat.sto(2, 3, bigint_qvec[8]);
        qmat.sto(3, 2, bigint_qvec[8]);
        qmat.sto(3, 3, 2 * bigint_qvec[9]);
    }
    return qmat;
}

void trim(string &s)
{
    size_t i;

    while ((i = s.find(" ", 0)) != string::npos || (i = s.find("\t", 0)) != string::npos ||
           (i = s.find("\n", 0)) != string::npos)
        s.erase(i, 1);
}

int trim_deng(string &s)
{
    auto i = s.find("=", 0);
    if (i != string::npos)
    {
        s.erase(0, i + 1);
    }
    if (s[i - 1] == 'x')
    {
        return 0;
    }
    else if (s[i - 1] == 'y')
    {
        return 1;
    }
    else if (s[i - 1] == 'z')
    {
        return 2;
    }
    else if (s[i - 1] == 'w')
    {
        return 3;
    }
    else if (s[i - 1] == 'a')
    {
        return 4;
    }
}

string delay_no_doubt(string &s)
{
    string ans = "";

    int cur = 0;
    while (cur < s.size())
    {
        if (s[cur] == 'u')
        {
            if (s[cur + 1] == '^')
            {
                ans = ans + "pow(u," + s[cur + 2] + ')';
                cur += 3;
            }
            else
            {
                ans += s[cur++];
            }
        }
        else if (s[cur] == 's') // sqrt
        {
            cur += 5;
            ans += "sqrt(";
            string temp = "";
            // sqrt(-2 + 2 * sqrt(2))
            int cnt_left = 1;

            while (cnt_left)
            {
                if (s[cur] == '(')
                {
                    cnt_left++;
                    temp += s[cur++];
                }
                else if (s[cur] == ')')
                {
                    cnt_left--;
                    if (cnt_left)
                        temp += s[cur++];
                }
                else
                {
                    temp += s[cur++];
                }
            }
            // ans = ans + delay_no_doubt(temp) + ",0.5";
            ans = ans + delay_no_doubt(temp);
        }
        else
        {
            ans += s[cur++];
        }
    }

    return ans;
}

std::vector<string> spilt4ginac(string &str, string &delta)
{
    // spilt
    std::vector<string> spi;
    string temp = "";

    if (str[0] == '[') //[x,y,z,w]
    {
        trim(str);
        int pos;
        auto k = str.size();
        if ((pos = str.find("Delta=", 0)) != string::npos) // have delta
        {
            for (; pos < k; pos++)
            {
                delta += str[pos];
            }
        }
        // spilt
        for (int i = 1; i < k; i++)
        {
            if (str[i] == ',' || str[i] == ']')
            {
                if (temp.size())
                {
                    // temp += '+';
                    spi.push_back(temp);
                    temp = "";
                }
            }
            else
            {
                temp += str[i];
            }
        }
        if (delta.size())
        {
            // maybe delta need erase "Delta="
            trim_deng(delta);
            spi.push_back(delta);
        }
    }
    else
    {
        auto k = str.size();
        for (int i = 0; i < k; i++)
        {
            if (str[i] == '\n')
            {
                if (temp.size())
                {
                    // temp += '+';
                    spi.push_back(temp);
                    temp = "";
                }
            }
            else
            {
                temp += str[i];
            }
        }
        // trim
        for (int i = 0; i < spi.size(); i++)
        {
            trim(spi[i]);
            trim_deng(spi[i]);
        }
    }

    for (int i = 0; i < spi.size(); i++)
    {
        temp = delay_no_doubt(spi[i]);
        spi[i] = temp;
    }

    return spi;
}

std::vector<GiNaC::ex> str2ginacex(std::vector<string> spilt)
{
    std::vector<GiNaC::ex> ans;

    GiNaC::symbol u("u");
    GiNaC::symtab table;
    table["u"] = u;
    GiNaC::parser reader(table);
    GiNaC::ex temp;

    for (int i = 0; i < spilt.size(); i++)
    {
        temp = reader(spilt[i]);
        ans.push_back(temp);
    }

    return ans;
}

double quadric_sdf_function(double x, double y, double z, Eigen::Matrix<float, 10, 1> v)
{
    double a = v(0);
    double b = v(1);
    double c = v(2);
    double d = v(3);
    double e = v(4);
    double f = v(5);
    double g = v(6);
    double h = v(7);
    double i = v(8);
    double j = v(9);
    // return -(a * x * x + b * y * y + c * z * z + d * x * y + e * x * z + f * y * z + g * x + h * y + i * z + j);
    return (a * x * x + e * y * y + h * z * z + b * x * y + c * x * z + f * y * z + d * x + g * y + i * z + j);
}

void handleCurveComponent(const QIOutputInfo *outputInfo, const GinacParams &ginacParams,
                          std::vector<curve_equation> &all_curve_equation,
                          std::vector<vertex_equation> &all_vertex_equation, const unsigned int first_index, const unsigned int second_index)
{
    // cnt of collect intersection curve
    int pqrowrqr = -100;

    // std::cout << outputInfo->q1 << endl;
    // std::cout << outputInfo->q2 << endl;

    // std::cout << outputInfo->n_components << endl; // 相交分支数
    auto qweqwe = outputInfo->n_components;

    while (qweqwe--)
    {

        if (pqrowrqr > 1)
            break;

        std::cout << "{" + outputInfo->parametrizations[qweqwe].label + "}" << endl;
        if (outputInfo->parametrizations[qweqwe].label == "point")
        {
            Eigen::MatrixXd ppp(1, 3);
            // ppp=extract point coordinate()
            ppp << 0, 0, 0;
            all_vertex_equation.emplace_back(ppp);
            // Viewer.addPoint(ppp);

            // point
            std::cout << "pointpointpointpoint" << endl;

            continue;
        }
        else
        {
        }
        std::cout << "[\n" + outputInfo->parametrizations[qweqwe].param + "]" << endl;

        auto ajskdhasjkdhakj = outputInfo->parametrizations[qweqwe].param;
        

        string delta = "";
        auto spilt = spilt4ginac(ajskdhasjkdhakj, delta);

        GiNaC::ex temp;
        std::vector<GiNaC::ex> bgelaqiwro;
        std::vector<GiNaC::ex> surfacePrametrizations;
        // std::cout << "eaqution in ginac" << std::endl;
        GiNaC::parser reader((*GlobalVariable::getTableSPtr()));
        for (int i = 0; i < spilt.size(); i++)
        {
            // PRINT_VAR(spilt[i]);
            temp = reader(spilt[i]);
            bgelaqiwro.push_back(temp);
            // std::cout << temp << std::endl;
        }
        for(int i=0;i<4;i++)
        {
            temp = reader(outputInfo->surfacePrametrizations[i]);
            surfacePrametrizations.emplace_back(temp);
        }
        all_curve_equation.emplace_back(bgelaqiwro, ginacParams.u, ginacParams.Delta, first_index, second_index,surfacePrametrizations);

        pqrowrqr++;
    }
}

void handleCurveComponent(const QIOutputInfo *outputInfo,
                           const unsigned int first_index, const unsigned int second_index)
{
    // cnt of collect intersection curve
    int pqrowrqr = -100;

    // std::cout << outputInfo->q1 << endl;
    // std::cout << outputInfo->q2 << endl;

    // std::cout << outputInfo->n_components << endl; // 相交分支数
    auto qweqwe = outputInfo->n_components;

    while (qweqwe--)
    {

        if (pqrowrqr > 1)
            break;

        std::cout << "{" + outputInfo->parametrizations[qweqwe].label + "}" << endl;
        if (outputInfo->parametrizations[qweqwe].label == "point")
        {
            Eigen::MatrixXd ppp(1, 3);
            // ppp=extract point coordinate()
            ppp << 0, 0, 0;
            // all_vertex_equation.emplace_back(ppp);
            GlobalVariable::getAllVertexEquation()->emplace_back(ppp);
            // Viewer.addPoint(ppp);

            // point
            std::cout << "pointpointpointpoint" << endl;

            continue;
        }
        else
        {
        }
        std::cout << "[\n" + outputInfo->parametrizations[qweqwe].param + "]" << endl;

        auto ajskdhasjkdhakj = outputInfo->parametrizations[qweqwe].param;
        

        string delta = "";
        auto spilt = spilt4ginac(ajskdhasjkdhakj, delta);

        GiNaC::ex temp;
        std::vector<GiNaC::ex> bgelaqiwro;
        std::vector<GiNaC::ex> surfacePrametrizations;
        // std::cout << "eaqution in ginac" << std::endl;
        GiNaC::parser reader((*GlobalVariable::getTableSPtr()));
        for (int i = 0; i < spilt.size(); i++)
        {
            // PRINT_VAR(spilt[i]);
            temp = reader(spilt[i]);
            bgelaqiwro.push_back(temp);
            // std::cout << temp << std::endl;
        }
        for(int i=0;i<4;i++)
        {
            temp = reader(outputInfo->surfacePrametrizations[i]);
            surfacePrametrizations.emplace_back(temp);
        }
        GlobalVariable::getAllCurveEquation()->emplace_back(bgelaqiwro, first_index, second_index,surfacePrametrizations);
        // all_curve_equation.emplace_back(bgelaqiwro, first_index, second_index,surfacePrametrizations);

        pqrowrqr++;
    }
}

std::pair<unsigned int,unsigned int> intersectionDouble(const std::vector<bigint_matrix> &Q, std::map<pair<unsigned int,unsigned int>,pair<unsigned int,unsigned int>> &firstSecondToBeginEnd,const unsigned int first_index, const unsigned int second_index, const GinacParams &ginacParams,
                        std::vector<curve_equation> &all_curve_equation,
                        std::vector<vertex_equation> &all_vertex_equation)
{
    pair<unsigned int,unsigned int> beginEnd;
    pair<unsigned int,unsigned int> firstSecond={first_index,second_index};
    if(firstSecondToBeginEnd.find(firstSecond)!=firstSecondToBeginEnd.end())
    {
        std::cout<<"already intersect"<<std::endl;
        beginEnd=firstSecondToBeginEnd[firstSecond];
    }
    else
    {
        std::cout<<"no intersect"<<std::endl;
        beginEnd.first=all_curve_equation.size();
        
        PRINT_VAR(Q[first_index]);
        PRINT_VAR(Q[second_index]);

        quad_inter<bigint> qi_ic = intersection(Q[first_index], Q[second_index], true, std::cout);
        QIOutputter outputter;
        outputter.setOmitImaginaryParametrizations();
        outputter.output(qi_ic, Q[first_index], Q[second_index]);
        // writer->setOutputInformation(outputter.getOutput());
        // writer->setVerbosityLevel(3);
        // writer->write ();/** Produces the result to /dev/stdout by default */

        // QIOutputParametrization parametrizations[4]; /** List of parametrizations (one for each algebraic
        // 						   component of the intersection) */

        // aasdadadadas = outputter.getOutput();
        handleCurveComponent(outputter.getOutput(),ginacParams,all_curve_equation,all_vertex_equation,first_index,second_index);
        
        beginEnd.second=all_curve_equation.size();
        firstSecondToBeginEnd[firstSecond]=beginEnd;
    }
    return beginEnd;
}

std::pair<unsigned int,unsigned int> intersectionDouble(const std::vector<bigint_matrix> &Q, std::map<pair<unsigned int,unsigned int>,pair<unsigned int,unsigned int>> &firstSecondToBeginEnd,const unsigned int first_index, const unsigned int second_index)
{
    std::cout<<"intersectionDouble:begin"<<std::endl;
    auto QAABB=*GlobalVariable::getQAABBSPtr();
    if(!QAABB[first_index].intersects(QAABB[second_index]))
    {
        std::cout<<"intersectionDouble:two box no intersect"<<std::endl;
        return{-1,-1};
    }
    
    pair<unsigned int,unsigned int> beginEnd;
    pair<unsigned int,unsigned int> firstSecond={first_index,second_index};

    std::cout<<"intersectionDouble:all curve size= "<<GlobalVariable::getAllCurveEquation()->size()<<std::endl;

    if(firstSecondToBeginEnd.find(firstSecond)!=firstSecondToBeginEnd.end())
    {
        std::cout<<"intersectionDouble:already intersect"<<std::endl;
        beginEnd=firstSecondToBeginEnd[firstSecond];
    }
    else
    {
        beginEnd.first=GlobalVariable::getAllCurveEquation()->size();

        quad_inter<bigint> qi_ic = intersection(Q[first_index], Q[second_index], true, std::cout);
        QIOutputter outputter;
        outputter.setOmitImaginaryParametrizations();
        outputter.output(qi_ic, Q[first_index], Q[second_index]);
        // writer->setOutputInformation(outputter.getOutput());
        // writer->setVerbosityLevel(3);
        // writer->write ();/** Produces the result to /dev/stdout by default */

        // QIOutputParametrization parametrizations[4]; /** List of parametrizations (one for each algebraic
        // 						   component of the intersection) */

        // aasdadadadas = outputter.getOutput();
        handleCurveComponent(outputter.getOutput(),first_index,second_index);
        
        beginEnd.second=GlobalVariable::getAllCurveEquation()->size();
        firstSecondToBeginEnd[firstSecond]=beginEnd;
    }
    std::cout<<"intersectionDouble:all curve size= "<<GlobalVariable::getAllCurveEquation()->size()<<std::endl;
    std::cout<<"intersectionDouble:done "<<std::endl;
    return beginEnd;
}

