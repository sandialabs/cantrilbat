/**
 * @file tok_input_util.h
 *
 * $Author: hkmoffa $
 * $Revision: 508 $
 * $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */


#ifndef TOK_INPUT_UTIL_H
#define TOK_INPUT_UTIL_H

#ifndef TRUE
#  define TRUE  1
#  define FALSE 0
#endif

#ifndef BOOLEAN
# define BOOLEAN int
#endif

#include <vector>
#include <cstdio>

namespace TKInput
{

#ifndef MAX_INPUT_STR_LN
#define MAX_INPUT_STR_LN 16100
#endif
#ifndef MAX_TOKEN_STR_LN
#define MAX_TOKEN_STR_LN 255
#endif
#ifndef MAXTOKENS
#define MAXTOKENS   255
#endif
#ifndef MAXTOKENSBIG
#define MAXTOKENSBIG   25500
#endif

//==================================================================================================================================
//! TOKEN structure
/*!
 *  This structure is used to parse strings. The original string is
 *  tokenized into a set of tokens via separation wrt white space.
 *  Both the tokens and the original string are storred within the structure
 */
class TOKEN
{
public:
    //! Default constructor
    TOKEN();

    //! Regular constructor for initializing a %TOKEN structure
    /*!
     *   @param[in] str                 String input for the token
     *   @param[in] nstd_delims         Nonstandard delimiters for the tokenization. 
     *                                  This defaults to 0, which uses the white space delimiters
     *                                  Usually, it the standard white space characters (i.e., isspace()), " \t\n\f\r\v"
     */
    TOKEN(const char* str, const char* nstd_delims = 0);

    //! Copy constructor
    /*!
     *  @param[in] right   Object to be copied
     */
    TOKEN(const TOKEN& right);

    //!  Assignment operator
    /*!
     *   @param[in]  right   Object to be copied
     *
     *   @return     Returns the current copied object
     */
    TOKEN& operator=(const TOKEN& right);

    //!   Destructor
    ~TOKEN();

    //! Original string
    char* orig_str;

    //! Modified string. Start of white space is modified by putting in a null character to create separate strings
    char* tok_str;

    //! Vector of pointers to the start of tokens in the string
    std::vector <char*> tok_ptrV;

    //!  Number of tokens in the structure
    int ntokes;

    //! If we are to tokenize the character string with nonstandard DELIMS, store the string here
    char* nstd_delims_;
};


#define NO_DEFAULT_INT     -68361
#define NO_DEFAULT_BOOLEAN NO_DEFAULT_INT
#define NO_DEFAULT_DOUBLE  -1.241E+11
#define NO_DEFAULT_STR    "NO_DEFAULT"


/**************************************************************************/
/*           Prototypes for routines tok_input_util.c                     */
/**************************************************************************/

extern void    set_tok_input_print_flag(int);
extern BOOLEAN get_next_keyLine(FILE*, TOKEN*, TOKEN*);
extern int     tok_to_int(const TOKEN*, const int, const int,
                          const int,    BOOLEAN*);
extern int     str_to_int(const char*,    const int,    const int,
                          const int,    BOOLEAN*);
extern double  tok_to_double(const TOKEN*, const double, const double,
                             const double, BOOLEAN*);
extern double  str_to_double(const char*,    const double, const double,
                             const double, BOOLEAN*);
extern BOOLEAN tok_to_boolean(const TOKEN*, const int, BOOLEAN*);
extern BOOLEAN str_to_boolean(const char*, const int, BOOLEAN*);
extern char*   tok_to_string(const TOKEN*,  const int,  const int,
                             const char*, BOOLEAN*);
extern char*   str_to_string(const char*, const char*, BOOLEAN*);
extern int    scan_for_int(FILE*, const char*, const int, const int);
extern double scan_for_double(FILE*,  const char*, const double,
                              const double);
extern char*  scan_for_string(FILE*,  const char*, const int, const int);
extern BOOLEAN scan_for_boolean(FILE*, const char*);
extern int    scan_for_line(FILE*, const char*, char [], const char,
                            const int, bool errorPrinting = true);
extern int  read_line(FILE*, char [], const int);
extern int  read_string(FILE*, char [], const char);
extern int strip(char []);
extern void lower_case(char []);
extern char* TokToStrng(const TOKEN*);
extern int  stokenize(char*, const char*, char* [], const int);

extern int  stokenizeV(char*, const char*, std::vector<char*>& tok_ptr, const int);

extern BOOLEAN strmatch(const char*, const char*);
extern BOOLEAN strstrmatch(const char*, const char*);
extern BOOLEAN strtokmatch(const TOKEN*, const char*);
extern BOOLEAN toktokmatch(const TOKEN*, const TOKEN*);
extern BOOLEAN toktokincluded(const TOKEN* tok1, const TOKEN* tok2);
extern void fillTokStruct(TOKEN*, const char*);



extern void copyTokStruct(TOKEN*, const TOKEN*);

//!   Finds a match of one string against a list of strings.
/*! 
 *   Returns the position that the first match occurred. If no match occurred, returns -1.
 *   The comparisons ignore differences in white space.
 *
 *     @param[in]   str1        C  String to be matched
 *     @param[in]   list        C Pointer vector to C Strings to be matched
 *     @param[in]   numList     Numver of strings in the list vector
 *
 *     @return      Returns the position that the first match occurs, or -1 if no match occurs
 */
extern int in_char_list(const char* const str1, const char** const list, int numList);

//! Make a copy of a string returning the copy
/*!
 *  Mallocs a new character string and copies the old string to it. Note, this routine always mallocs something
 *
 *    NOTE: Memory leak may result if the calling program doesn't free
 *          the malloced space
 *
 *  @param[in] cstr   C string to be copied
 *  @return           Returns a malloced copy of a string
 */
char* copy_string(const char* const cstr);

//! Change the token by taking eliminating the iword'th token from the token structure and original string
/*!
 *    @param[in]     iword      Token number to be eliminated from the token structure
 *    @param[in,out] tok        Token object to be manipulated
 */
void strip_item_from_token(int iword, TOKEN* tok);

}
/***************************************************************************************************************************************/
#endif /* END OF TOK_INPUT_UTIL_H */
/***************************************************************************************************************************************/

