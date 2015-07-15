/**
 * @file tok_input_util.h
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
//!  Define true as 1
#define TRUE  1
//!  Define false as 0
#define FALSE 0
#endif

#ifndef BOOLEAN
//!  Define the BOOLEAN macro as an int
# define BOOLEAN int
#endif


#include <vector>
#include <cstdio>
//----------------------------------------------------------------------------------------------------------------------------------
//!  Namespace for the manipulation of TOKEN structures and for the reading and writing of ascii input files based on 
//!  reading each line and tokenizing their input.
/*!
 *  This is largely a C routine based namespace. It is efficient and fast.  It relies on the C standard library routines
 *  stokenize() to manipulate and break up character strings.
 */
namespace TKInput
{
//==================================================================================================================================
#ifndef MAX_INPUT_STR_LN
//!        This is the maximum line length in characters
#define MAX_INPUT_STR_LN 16100
#endif
#ifndef MAXTOKENS
//!        This is the maximum number of tokens in a TOKEN structure
#define MAXTOKENS   255
#endif
#ifndef MAXTOKENSBIG
//!        This is the maximum number of tokens in a large string that is tokenized
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
//==================================================================================================================================
//!  Special int value that denotes that there is no default value that can be supplied
#define NO_DEFAULT_INT     -68361

//!  Special int value that denotes that there is no default boolean that can be supplied
#define NO_DEFAULT_BOOLEAN NO_DEFAULT_INT

//!  Special double value that denotes that there is no default double value that can be supplied
#define NO_DEFAULT_DOUBLE  -1.241E+11

//!  Special character string that denotes that there is no default character string value that can be supplied
#define NO_DEFAULT_STR    "NO_DEFAULT"

/**********************************************************************************************************************************/
/*                                     Prototypes for routines tok_input_util.cpp                                                 */
/**********************************************************************************************************************************/

//! This sets the static variable TKInput::PrintInputFile, which determines whether the input file
//! is copied to the stdout during the reading of the input deck
/*!
 *    @param[in]  print_flag    Value for the variable, 1 means to print, and 0 means to be silent.
 */
void set_tok_input_print_flag(int print_flag);

//!  This routine reads the input file to obtain the next line of uncommented data. The results are returned in two
//! TOKEN structures. keyLineTok contains the key Line (everything before the first equals sign).
//!   keyArgTok contains everything after the equals sign.
//!  Note - Either keyLineTok or keyArgTok may be the null token (but not both)
/*!
 *
 *  @param[in]    ifp                       FILE pointer to read lines from
 *  @param[out]   keyLineTok                Pointer to the TOKEN structure for everything before the equals sign
 *  @param[out]   keyArgTok                 Pointer to the TOKEN structure for everything after the equals sign
 *
 *   The definition of a token structure, given in .h file, is as follows:
 *
 *      struct TOKEN {
 *        char  orig_str[MAX_INPUT_STR_LN + 1];
 *        char  tok_str[MAX_INPUT_STR_LN + 1];
 *        char *tok_ptr[MAXTOKENS];
 *        int   ntokes;
 *      };
 *   
 *     orig_str Contains the original string, unmodified.
 *     tok_str  Contains a modified version of the string,
 *              whose positions
 *              are pointed to by tok_ptr[i] values. It is usually not
 *              referenced directly.
 *     tok_ptr[i] Contains the i_th token of the original string. This
 *              is a stripped character string. 0 <=i <= i-1
 *     ntokes   Number of tokens in the string.
 *
 *
 *   Comments are denoted by either '!' or '#'.
 *   Everything after the comment character on a
 *   line is stripped first. The comment character can occur
 *   anywhere in the line.
 *
 *   Arguments to the keyLine are denoted by everything after a
 *   '=' character on the line.
 *
 *   Example:
 *   ---------------------------------------------------------------
 *   ! Jack and Jill went up the hill to fetch a pale of water
 *   ! Jack by nimble, Jack be swift; Jack jump over the candle stick
 *
 *   The meaning of life is =    36.243     24  136 Not even ! close
 *   -----------------------------------------------------------------
 *
 *   Then, the routine would return (amongst other things):
 *     keyLineTok->orig_str = "The meaning of life is"
 *     keyArgTok->orig_str = "36.243     24   36 Not even"
 *     keyArgTok->ntokes = 5
 *     keyArgTok->tok_ptr[0] = "36.243"
 *     keyArgTok->tok_ptr[1] = "24"
 *     keyArgTok->tok_ptr[2] = "136"
 *     keyArgTok->tok_ptr[3] = "Not"
 *     keyArgTok->tok_ptr[4] = "Even"
 *
 *   @return             The function returns TRUE if there is a next line to process.
 *                       It returns false if an EOF is encountered.
 */
BOOLEAN get_next_keyLine(FILE* ifp, TOKEN* keyLineTok, TOKEN* keyArgTok);

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


//!      Interprets the arguments in a TOKEN structure as a string.
//!      It mallocs new space for the string, and returns the pointer to it.
/*!
 *      The number of tokens in the string is checked before returning.
 *      The value must be between the maxTok and minTok; it can equal the
 *      max or min value.
 *
 *      A default may be specified on this routine. If the TOKEN hasn't been filled yet,
 *      the value of the default is returned. The absence of a
 *      default may also be specified by setting devault_value to
 *      NO_DEFAULT_INT.
 *
 *      If there is an error, *error is set to TRUE. *error isn't touched
 *      if there isn't an error.
 *
 *  @param[in]    tokPtr       Pointer to the input TOKEN structure
 *  @param[in]    maxTok       Maximum number of tokens
 *  @param[in]    minTok       Minimum number of tokens
 *  @param[in]    defaultVal   Input default value of the c string. If there is no default use the
 *                             value NO_DEFAULT_STR
 *  @param[out]   error        Set to True if there is an error condition on the card
 *
 *  @return                    Returns the value of the string as a c string, i.e., null terminated. This string
 *                             is always malloced, and therefore always has to be freed by the calling program.
 */
char* tok_to_string(const TOKEN* tokPtr, const int maxTok, const int minTok, const char* defaultVal, BOOLEAN* error);

//! Interprets the argument string as a string, then mallocs new space for the string, and returns the pointer to it
/*!
 *      Certain ascii strings are checked for first (case is insignificant):
 *
 *                 String              Retn_Value
 *                 ---------        --------------
 *                 default            defaultVal
 *
 *      A default may be specified on the command line. The absence of a default may also be specified by setting
 *      devault_value to NO_DEFAULT_INT.
 *
 *      If there is an error, *error is set to TRUE. *error isn't touched if there isn't an error.
 *
 *  @param[in]    str          Input c string to be interpreted
 *  @param[in]    defaultVal   Input default value of the c string
 *  @param[out]   error        True if there is an error condition on the card
 *
 *  @return                    Returns the value of the string as a c string, i.e., null terminated. This string
 *                             is always malloced.
 */
char* str_to_string(const char* str, const char* defaultVal, BOOLEAN* error);

extern int    scan_for_int(FILE*, const char*, const int, const int);
extern double scan_for_double(FILE*,  const char*, const double,
                              const double);
extern char*  scan_for_string(FILE*,  const char*, const int, const int);
extern BOOLEAN scan_for_boolean(FILE*, const char*);
extern int    scan_for_line(FILE*, const char*, char [], const char,
                            const int, bool errorPrinting = true);
extern int  read_line(FILE*, char [], const int);
extern int  read_string(FILE*, char [], const char);

//!  This routine strips off blanks and tabs (only leading and trailing characters) 
/*!
 *   Comments are excluded -> All instances of the comment character, '!', are replaced by '\0' thereby terminating the string.
 *
 *    @param[in,out]  str              On output 'str' contains the same characters as on
 *                                     input except the leading and trailing white space and comments have been removed.
 *
 *    @return                          Returns the number of characters remaining in the string (excluding the null character).
 */
int strip(char str[]);

//!  Translates a string delimited by a NULL character to lower case.
/*!
 *  There is no error checking in this version.  Relies on stlib function, tolower.
 *
 *   @param[in,out]      str            On input it contains the string. On output it contains the converted string.
 */
void lower_case(char str[]);

//!  Mallocs a new character string and copies the tokens character string to it, appending all tokens together
//!  into a single string separated by a single space character.
/*!
 *      @param[in]   keyptr             Pointer to the TOKEN structure to be copied
 *
 *      @return                         Returns the pointer to the newly malloced string. The new string should be freed 
 *                                      when no longer needed.
 */
char* TokToStrng(const TOKEN* keyptr);

//! This function will break up a string into its respective "tokens".
//!  It returns the number of tokens found. See the strtok(3) man page.
/*!
 *
 *  @param[in]      string                    String to be tokenized.  Note, that the string is
 *                                            changed by this procedure. Null characters are put between each symbol.
 *
 *  @param[in]      delimiters                String containing a list of delimiters. The example below covers 'white space'
 *                                            e.g., char *delimiters  = " \t\n"
 *  @param[in]      max_tokens                Maximum number of tokens to be found
 *
 *  @param[out]    tok_ptr                    Vector of pointers to char*. This must already have been malloced before entrance
 *                                            to the routine, with length max_tokens.
 *
 *  @return                                   Returns the number of tokens in the string
 */
int stokenize(char* string, const char* delimiters, char* tok_ptr[], const int max_tokens);

//! This function will break up a string into its respective "tokens".
//!  It returns the number of tokens found. See the strtok(3) man page.
/*!
 *
 *  @param[in]      string                    String to be tokenized.  Note, that the string is
 *                                            changed by this procedure. Null characters are put between each symbol.
 *
 *  @param[in]      delimiters                String containing a list of delimiters. The example below covers 'white space'
 *                                            e.g., char *delimiters  = " \t\n"
 *  @param[in]      max_tokens                Maximum number of tokens to be found
 *
 *  @param[out]     tok_ptrV                  Vector of pointers to strings, that contain the input string's tokens
 *
 *  @return                                   Returns the number of tokens in the string
 */
int stokenizeV(char* string, const char* delimiters, std::vector<char*>& tok_ptrV, const int max_tokens);

//!   This routine checks whether one string is the same as another.
//!   Upper case is transformed into lower case before the comparison is done.
//!   Thus, case doesn't matter in the comparison. However, white space does matter in this comparison.
/*!
 *    @param[in]        s1                Const pointer to null terminated C string
 *    @param[in]        s2                Const pointer to null terminated C string
 *
 *    @return                             If there is a match it returns true. If they aren't, it returns false.
 */
BOOLEAN strmatch(const char* s1, const char* s2);

//! This routine checks whether two strings are the same modulo differences in their white space
/*!
 *    @param[in]        s1                Const pointer to null terminated C string
 *    @param[in]        s2                Const pointer to null terminated C string
 *
 *    @return                             If there is a match it returns true. If they aren't, it returns false.
 */
BOOLEAN strstrmatch(const char* s1, const char* s2);

//!  THIS routine checks whether a string matches the string contained in the tokens of a keyLineStr.
//!  White space and case are ignored.
/*!
 *    @param[in]        keyptr            Const pointer to first TOKEN structure
 *    @param[in]        s2                Const pointer to null terminated C string
 *
 *    @return   If there is a match it returns true. If they aren't, it returns false.
 */
BOOLEAN strtokmatch(const TOKEN* keyptr, const char* s2);

//! This routine checks whether two %TOKEN structures contain the same data up to differences in white space.
//!  Case is ignored as well, as strmatch is called.
/*!
 *  @param[in]           keyptr1                   Const pointer to first TOKEN structure
 *  @param[in]           keyptr2                   Const pointer to second TOKEN structure
 *
 *  @return                                        Returns true if the two TOKENs match, false otherwise
 */
BOOLEAN toktokmatch(const TOKEN* keyptr1, const TOKEN* keyptr2);

//!  This routine checks whether the first TOKEN structure is included in the second TOKEN structure.
//!  Case and whitespace are ignored, as strmatch is called.
/*!
 *    If the first TOKEN is empty, this routine returns true.
 *
 *  @param[in]           keyptr1                   Const pointer to first TOKEN structure
 *  @param[in]           keyptr2                   Const pointer to second TOKEN structure
 *
 *  @return                                        Returns true if the first TOKEN is included in the second TOKEN pointer.
 */
BOOLEAN toktokincluded(const TOKEN* keyptr1, const TOKEN* keyptr2);

//! Fill in a %TOKEN structure with a string. Use the definition of white space
//! at the start of the file to tokenize the string, storring it in the TOKEN structure.
/*!
 *   @param[in]       s2          Null terminated C string that will be copied into the %TOKEN structure
 *   @param[out]      keyptr1     Pointer to the TOKEN structure that will be filled by the character string.
 */
void fillTokStruct(TOKEN* keyptr1, const char* s2);

//!  Copies the information stored in TOKEN keyptr2 into TOKEN keyptr1
/*!
 *    @param[in]     keyptr2         Pointer to the TOKEN object to be copied
 *    @param[out]    keyptr1         Pointer to the TOKEN object to be copied into.
 */
void copyTokStruct(TOKEN* keyptr1, const TOKEN* keyptr2);

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

