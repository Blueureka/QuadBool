/** ********************************** */
/** Part of the QIParser.h declaration */
/** ********************************** */

#ifndef _qi_polynomial_parser_h_
#define _qi_polynomial_parser_h_

#include <libqi/io/QIParser.h>
#include <libqi/rpl/bigint.h>

unsigned short _current_symbol;

/** I. */
/** *************************** */
/** Input stream symbol reading */
/** *************************** */
#define CURRENT_SYMBOL _current_symbol
#define NEXT_SYMBOL CURRENT_SYMBOL++
#define PREV_SYMBOL --CURRENT_SYMBOL

#define GET_CURRENT_SYMBOL() ((_current_symbol < _string_length) ? (_quadricDesc.c_str()[CURRENT_SYMBOL]) : '\0')
#define GET_NEXT_SYMBOL() ((_current_symbol < _string_length) ? (_quadricDesc.c_str()[NEXT_SYMBOL]) : '\0')
#define GET_PREV_SYMBOL() ((_current_symbol > 0) ? (_quadricDesc.c_str()[PREV_SYMBOL]) : _quadricDesc[0])

/** II. */
/** ***************************************************** */
/** Monomial parsing: see the file "monomial_parsing.txt" */
/** ***************************************************** */
short m_code[4];   /** Power of X,Y,Z,W */
short m_sum;       /** Value of X+Y+Z+W */
short current_var; /** Last encountered variable */

/** Current monomial's sign */
short m_sign;

/** Current monomial's coefficient value */
rpl::bigint m_coefficient;

/** Associate an indice inside m_code for each variable X,Y,Z,W */
#define VAR_X 3
#define VAR_Y 2
#define VAR_Z 1
#define VAR_W 0

/** Resets the code */
#define reset_code m_sum = m_code[0] = m_code[1] = m_code[2] = m_code[3] = 0

/** Update each variable's power value as it is encountered */
#define raise_var ++(m_code[current_var])
#define raise_x          \
    current_var = VAR_X; \
    raise_var
#define raise_y          \
    current_var = VAR_Y; \
    raise_var
#define raise_z          \
    current_var = VAR_Z; \
    raise_var
#define raise_w          \
    current_var = VAR_W; \
    raise_var

/** Sums each digit of the code */
#define sum_code m_sum = m_code[0] + m_code[1] + m_code[2] + m_code[3]

/** Homogenises the code */
#define hom_code m_code[0] += (2 - m_sum)

/** A general purpose macro that checks the code
    and homogenises it if necessary */
#define check_code                                 \
    sum_code;                                      \
    if (m_sum > 2)                                 \
    {                                              \
        throw "Invalid monomial: " + decode();     \
    }                                              \
    else if (m_sum < 2)                            \
    {                                              \
        if (m_code[0] != 0)                        \
        {                                          \
            throw "Invalid monomial: " + decode(); \
        }                                          \
        else                                       \
            hom_code;                              \
    }

/** Does W+2Z+4Y+8X  */
#define hash_code m_code[0] + (m_code[1] << 1) + (m_code[2] << 2) + (m_code[3] << 3)

/** Returns the indice inside the vector */
#define indice_vector m_hash[hash_code]

// #define init()         \
//     reset_code;        \
//     m_coefficient = 1; \
//     m_sign = 1;\
//     printf("init() called at %s:%d\n", __FILE__, __LINE__);

/** When a "+" or a "-" is encountered, we store the
    previously parsed monomial (if exists), then
    we reinit all to prepare for the next monomial, via
    the "init" macro */
#define next_coeff()                                            \
    try                                                         \
    {                                                           \
        check_code;                                             \
    }                                                           \
    catch (string e)                                            \
    {                                                           \
        throw e;                                                \
    }                                                           \
    if (indice_vector >= 0 && indice_vector < 10)               \
        _vectorialDesc[indice_vector] = m_sign * m_coefficient; \
    reset_code;                                                 \
    m_coefficient = 1;                                          \
    m_sign = 1;

/** III. */
/** ****************************************** */
/** Templates defining parsing functions for   */
/** each non terminal of the grammar.          */
/** ****************************************** */

/** Handy macros to define prototypes and skeletons */

/** Header */
#define PARSE_DECL(name) void parse_##name(void) noexcept(false)

/** Declaration (used in the QIParser's implementation file */
#define PARSE_IMPL(name) void QIParser::parse_##name(void) noexcept(false)

/** Call to a non terminal checking function */
#define IS_A(name) parse_##name();

/** Booleans used in logical expression */
#define decl_is_a(name) bool is_a_##name = false

/** Automatically declares a prototype for each
    non terminal symbol of the grammar. */
PARSE_DECL(poly);
PARSE_DECL(term);
PARSE_DECL(mono);
PARSE_DECL(var);
PARSE_DECL(coef);
PARSE_DECL(sign);
PARSE_DECL(eol);

/** IV. */
/** **************** */
/** Terminal parsing */
/** **************** */
/** Please note that if the read character doesn't match, we *don't*
 *  advance to the next character, i.e. "NEXT_SYMBOL" is skipped by
 *  the "throw" statement. */
#define PARSE_TERMINAL(character)                                                                     \
    do                                                                                                \
    {                                                                                                 \
        if (!(GET_CURRENT_SYMBOL() == character))                                                     \
        {                                                                                             \
            stringstream s;                                                                           \
            s << "Parse error: '" << GET_CURRENT_SYMBOL() << "' at position: " << CURRENT_SYMBOL + 1; \
            throw s.str();                                                                            \
        }                                                                                             \
        NEXT_SYMBOL;                                                                                  \
    } while (0)

/** V. */
/** *************************************** */
/** Debugging function                      */
/** Prints the current monomial as a string */
/** *************************************** */
inline string decode()
{
    stringstream ss;

    if (m_code[3] > 0)
    {
        ss << "x";
        if (m_code[3] >= 2)
            ss << "^" << m_code[3];
        ss << " ";
    }
    if (m_code[2] > 0)
    {
        ss << "y";
        if (m_code[2] >= 2)
            ss << "^" << m_code[2];
        ss << " ";
    }
    if (m_code[1] > 0)
    {
        ss << "z";
        if (m_code[1] >= 2)
            ss << "^" << m_code[1];
        ss << " ";
    }
    if (m_code[0] > 0)
    {
        ss << "w";
        if (m_code[0] >= 2)
            ss << "^" << m_code[0];
    }
    return ss.str();
}

#endif

// #undef init