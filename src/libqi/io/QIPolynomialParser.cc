
#include <libqi/io/QIPolynomialParser.h>

/** ************************************** */
/** Part of the QIParser.cc implementation */
/** ************************************** */

/** Flags involved in a boolean expression */
/** Non-terminal symbols's flags */
decl_is_a(term);
decl_is_a(mono);
decl_is_a(var);
decl_is_a(coef);
decl_is_a(sign);
decl_is_a(eol);

/** Terminal symbol's flags */
decl_is_a(mult);
decl_is_a(power);
decl_is_a(two);

const short m_hash[17] = {-1, -1, 9, 8, 7, 6, 5, -1, 4, 3, 2, -1, 1, -1, -1, -1, 0};

PARSE_IMPL(poly)
{

    /** Initializes internal parser data */

    // init ();
    reset_code;
    m_coefficient = 1;
    m_sign = 1;

    /** First sign is not mandatory, so the
        parsing exception is trapped */
    try
    {
        IS_A(sign);
    }
    catch (...)
    {
    }
    try
    {
        IS_A(term);
    }
    catch (string e)
    {
        throw e;
    }

    for (;;)
    {

        try
        {

            is_a_sign = false;
            is_a_term = false;
            try
            {
                IS_A(sign);
                is_a_sign = true;
            }
            catch (string e)
            {
                throw e;
            }
            try
            {
                IS_A(term);
                is_a_term = true;
            }
            catch (string e)
            {
                throw e;
            }
        }
        catch (string e)
        {

            if (!is_a_sign)
                break;
            if (is_a_sign && !is_a_term)
            { /*cout << "FALSE" ;*/
                throw e;
            }
        }
    }

    try
    {
        IS_A(eol);
    }
    catch (string e)
    { /*cout << "FALSE" << endl ;*/
        throw e;
    }

    /** Finalizes the last coefficient (stores it) */
    /*  next_coeff();*/
}

PARSE_IMPL(term)
{

    string error;

    is_a_coef = false;
    is_a_mono = false;

    /** Note: goto is usually a bad way of programming, but there are special
        cases like parsing where it's useful and quick. */
    try
    {
        IS_A(coef);
        is_a_coef = true;
        goto next;
    }
    catch (string e)
    {
        error = e;
    }
    try
    {
        IS_A(mono);
        is_a_mono = true;
    }
    catch (string e)
    {
        error = e;
    }

next:
    if (!is_a_coef && !is_a_mono)
        throw error;

    for (;;)
    {

        try
        {

            is_a_mult = false;
            is_a_mono = false;
            try
            {
                PARSE_TERMINAL('*');
                is_a_mult = true;
            }
            catch (string e)
            {
                throw e + ": expected symbol was: '*'";
            }
            try
            {
                IS_A(mono);
                is_a_mono = true;
            }
            catch (string e)
            {
                throw e;
            }
        }
        catch (string e)
        {

            if (!is_a_mult)
                break;
            if (is_a_mult && !is_a_mono)
            { /*cout << "FALSE" << endl;*/
                throw e;
            }
        }
    }

    next_coeff();
}

PARSE_IMPL(mono)
{

    string error;

    try
    {
        IS_A(var);
    }
    catch (string e)
    { /*cout << "FALSE" << endl;*/
        throw e;
    }

    is_a_power = false;
    is_a_two = false;
    try
    {
        PARSE_TERMINAL('^');
        is_a_power = true;
    }
    catch (string e)
    {
        error = e + ": expected symbol was: '^'";
    }
    try
    {
        PARSE_TERMINAL('2');
        raise_var;
        is_a_two = true;
    }
    catch (string e)
    {
        error = e + ": expected symbol was: '2'";
    }

    if (is_a_power && !is_a_two)
    { /*cout << "FALSE"  << endl;*/
        throw error;
    }
}

PARSE_IMPL(var)
{

    try
    {
        PARSE_TERMINAL('w');
        raise_w;
        return;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('W');
        raise_w;
        return;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('x');
        raise_x;
        return;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('X');
        raise_x;
        return;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('y');
        raise_y;
        return;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('Y');
        raise_y;
        return;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('z');
        raise_z;
        return;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('Z');
        raise_z;
        return;
    }
    catch (string e)
    { /*cout << "FALSE" ;*/
        throw e + ": expected symbol was: [w|W|x|X|y|Y|z|Z]";
    }
}

PARSE_IMPL(coef)
{

    /** See remark inside PARSE(term) above about the usage of goto */

    /** At least one digit */
    try
    {
        PARSE_TERMINAL('0');
        m_coefficient = 0;
        goto nextdigit;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('1');
        m_coefficient = 1;
        goto nextdigit;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('2');
        m_coefficient = 2;
        goto nextdigit;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('3');
        m_coefficient = 3;
        goto nextdigit;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('4');
        m_coefficient = 4;
        goto nextdigit;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('5');
        m_coefficient = 5;
        goto nextdigit;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('6');
        m_coefficient = 6;
        goto nextdigit;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('7');
        m_coefficient = 7;
        goto nextdigit;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('8');
        m_coefficient = 8;
        goto nextdigit;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('9');
        m_coefficient = 9;
        goto nextdigit;
    }
    catch (string e)
    { /*cout << "FALSE" << endl ;*/
        throw e + ": expected symbol was: [0-9]";
    }

nextdigit:
    /** Some extra digits */
    for (;;)
    {

        try
        {
            PARSE_TERMINAL('0');
            m_coefficient = (10 * m_coefficient) + 0;
            continue;
        }
        catch (...)
        {
        }
        try
        {
            PARSE_TERMINAL('1');
            m_coefficient = (10 * m_coefficient) + 1;
            continue;
        }
        catch (...)
        {
        }
        try
        {
            PARSE_TERMINAL('2');
            m_coefficient = (10 * m_coefficient) + 2;
            continue;
        }
        catch (...)
        {
        }
        try
        {
            PARSE_TERMINAL('3');
            m_coefficient = (10 * m_coefficient) + 3;
            continue;
        }
        catch (...)
        {
        }
        try
        {
            PARSE_TERMINAL('4');
            m_coefficient = (10 * m_coefficient) + 4;
            continue;
        }
        catch (...)
        {
        }
        try
        {
            PARSE_TERMINAL('5');
            m_coefficient = (10 * m_coefficient) + 5;
            continue;
        }
        catch (...)
        {
        }
        try
        {
            PARSE_TERMINAL('6');
            m_coefficient = (10 * m_coefficient) + 6;
            continue;
        }
        catch (...)
        {
        }
        try
        {
            PARSE_TERMINAL('7');
            m_coefficient = (10 * m_coefficient) + 7;
            continue;
        }
        catch (...)
        {
        }
        try
        {
            PARSE_TERMINAL('8');
            m_coefficient = (10 * m_coefficient) + 8;
            continue;
        }
        catch (...)
        {
        }
        try
        {
            PARSE_TERMINAL('9');
            m_coefficient = (10 * m_coefficient) + 9;
            continue;
        }
        catch (...)
        {
        }

        /** If no more digits are found, return. */
        return;
    }
}

PARSE_IMPL(sign)
{

    /** If we encounter a "+" or a "-", the "sign" variable contains
        the sign of the expression, "1" for a positive value, "-1" otherwise.
        Default value is "1" */
    try
    {
        PARSE_TERMINAL('+'); /*next_coeff();*/
        return;
    }
    catch (...)
    {
    }
    try
    {
        PARSE_TERMINAL('-'); /*next_coeff();*/
        m_sign = -1;
        return;
    }
    catch (string e)
    { /*cout << "FALSE" << endl ;*/
        throw e + ": expected symbol was: [+|-]";
    }
}

PARSE_IMPL(eol)
{

    try
    {
        PARSE_TERMINAL('\0');
        return;
    }
    catch (string e)
    { /*cout << "FALSE" << endl ;*/
        throw e + ": expected symbol was end of string";
    }
}
