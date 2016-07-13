/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 * $RCSfile: ctmlDiff.cpp,v $
 * $Author: hkmoffa $
 * $Date: 2013-03-26 10:44:21 -0600 (Tue, 26 Mar 2013) $
 * $Revision: 571 $
 * $Name:  $
 *====================================================================*/
/*
 *
 *  xmlSolnDiff File1.xml File2.xml
 *
 *  Compares the variable values in two Cantera solution xml
 *  files.
 *  The comparison is done using a weighted norm basis.
 *
 *  The two files should be basically equal. However, File1.xml is
 *  taken as the reference file, that has precedence, when there is
 *  something to be decided upon.
 *
 *  Arguments:
 *   -h = prints this usage information
 *
 *  Shell Return Values
 *    1 = Comparison was successful
 *    0 = One or more nodal values failed the comparison
 *   -1 = Apples to oranges, the files can not even be compared against
 *        one another.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#if defined(__CYGWIN__)
#include <getopt.h>
#endif

#include "tok_input_util.h"

using namespace std;
using #ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif;
using namespace TKInput;

#ifdef WIN32
#pragma warning(disable:4996)
#endif
#ifndef MAX
#  define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif
#ifndef MIN
#  define MIN(x,y)    (( (x) < (y) ) ? (x) : (y))
#endif
#ifndef SAFE_DELETE
#define SAFE_DELETE(a) if (a) { delete (a) ; a = 0 ; }
#endif

int Debug_Flag = 1;
double grtol = 1.0E-3;
double gatol = 1.0E-9;
std::map<std::string, double> VarRtol;
std::map<std::string, double> VarAtol;

/*
 * Possible entries in a value list
 */
enum valueType_enum {
    UNIDENTIFIED_VT,
    FLOAT_VT,
    FLOATLIST_CSV_VT,
    INTEGER_VT,
    INTEGERLIST_CSV_VT,
    STRING_VT,
    STRING_CSV_VT,
    PAIR_STRING_FLOAT_CSV_VT,
    MATRIX_STRING_STRING_FLOAT_CSV_VT
};

//======================================================================================================================
#ifdef WINMSVC
/*
 * Windows doesn't have getopt(). This is an incomplete version that
 * does enough to handle required functionality.
 */
int optind = -1;
char* optarg = 0;

int getopt(int argc, char** argv, const char*)
{
    static int currArg = 1;
    static int currOptInd = 1;
    string tok;
    static int charPos = 0;
    int rc = -1;
    if (currArg >= argc) {
        optarg = 0;
        return -rc;
    }
    tok = string(argv[currArg]);
    currOptInd = currArg+1;
    if (currOptInd > argc - 1) {
        currOptInd = -1;
        optarg = 0;
    } else {
        optarg = argv[currArg+1];
    }
    size_t len = strlen(tok.c_str());
    if (charPos == 0) {
        bool found = false;
        do {
            tok = string(argv[currArg]);
            len = strlen(tok.c_str());
            if (len > 1 && tok[0] == '-') {
                found = true;
                charPos = 1;
                if (len > 2 && tok[1] == '-') {
                    charPos = 2;
                }
            } else {
                if (optind == -1) {
                    optind = currArg;
                }
            }
            if (!found) {
                if (currArg < (argc-1)) {
                    currArg++;
                } else {
                    optarg = 0;
                    return -1;
                }
            }
        } while (!found);
    }

    rc = tok[charPos];
    if (charPos < static_cast<int>(len - 1)) {
        charPos++;
    } else {
        charPos = 0;
    }
    return rc;
}

#endif
//======================================================================================================================
double rtolVar(std::string var)
{
    std::map<std::string, double>::const_iterator i = VarRtol.find(var);
    if (i != VarRtol.end()) {
        return i->second;
    } else {
        VarRtol[var] = grtol;
    }
    return VarRtol[var];
}
//======================================================================================================================
double atolVar(std::string var)
{
    std::map<std::string, double>::const_iterator i = VarAtol.find(var);
    if (i != VarAtol.end()) {
        return i->second;
    } else {
        VarAtol[var] = grtol;
    }
    return VarAtol[var];
}
//======================================================================================================================
int diff_double(double d1, double d2, double rtol, double atol)

/*
 * Compares 2 doubles. If they are not within tolerance, then this
 * function returns true.
 */
{
    if (fabs(d1 - d2) > (atol + rtol * 0.5 * (fabs(d1) + fabs(d2)))) {
        return 1;
    }
    return 0;
}
//======================================================================================================================
static int diff_double_slope(double d1, double d2, double rtol,
                             double atol, double xtol, double slope1, double slope2)

/*
 * Compares 2 doubles. If they are not within tolerance, then this
 * function returns true.
 */
{
    double atol2 = xtol * (fabs(slope1) + fabs(slope2));
    if (fabs(d1 - d2) > (atol + atol2 + rtol * 0.5 * (fabs(d1) + fabs(d2)))) {
        return 1;
    }
    return 0;
}
//======================================================================================================================
static double calc_rdiff(double d1, double d2, double rtol, double atol)
{
    double rhs, lhs;
    rhs = fabs(d1 - d2);
    lhs = atol + rtol * 0.5 * (fabs(d1) + fabs(d2));
    return (rhs / lhs);
}
//======================================================================================================================
static double get_atol(std::vector<double>& values, const int nvals, const double atol)
{
    double sum = 0.0, retn;
    if (nvals <= 0) {
        return gatol;
    }
    for (int i = 0; i < nvals; i++) {
        retn = values[i];
        sum += retn * retn;
    }
    sum /= nvals;
    retn = sqrt(sum);
    return ((retn + 1.0) * atol);
}
//======================================================================================================================
bool compareFieldVectorFlts(std::vector<double>& comp1,
                            std::vector<double>& comp2,
                            std::vector<double>& xpos,
                            std::string varName)
{

    double rtol = rtolVar(varName);
    double atol = atolVar(varName);

    double max_diff = 0.0;
    int ndiff = 0;
    double rel_diff = 0.0;
    double xatol, slope1, slope2;

    int nDataRows1 = comp1.size();
    int nDataRows2 = comp2.size();
    int nDataRowsMIN = MIN(nDataRows1, nDataRows2);
    int nDataRowsMAX = MAX(nDataRows1, nDataRows2);
    // Adjust for large values
    double atol_j = get_atol(comp1, nDataRows1, atol);
    atol_j = MIN(atol_j, get_atol(comp2, nDataRows2, atol));
    int jmax = 0;

    for (int j = 0; j < nDataRowsMIN; j++) {

        slope1 = 0.0;
        slope2 = 0.0;

        if (j == 0) {
            slope1 = (comp1[j + 1] - comp1[j]) / (xpos[j + 1] - xpos[j]);
            slope2 = (comp2[j + 1] - comp2[j]) / (xpos[j + 1] - xpos[j]);
            xatol = fabs(grtol * (xpos[1] - xpos[0]));
        } else if (j == (nDataRowsMIN - 1)) {
            slope1 = (comp1[j] - comp1[j - 1]) / (xpos[j] - xpos[j - 1]);
            slope2 = (comp2[j] - comp2[j - 1]) / (xpos[j] - xpos[j - 1]);
            xatol = fabs(grtol * (xpos[j] - xpos[j - 1]));
        } else {
            slope1 = (comp1[j + 1] - comp1[j - 1]) / (xpos[j + 1] - xpos[j - 1]);
            slope2 = (comp2[j + 1] - comp2[j - 1]) / (xpos[j + 1] - xpos[j - 1]);
            xatol = fabs(grtol * 0.5 * (xpos[j + 1] - xpos[j - 1]));
        }
        bool notOK = diff_double_slope(comp1[j], comp2[j],
                                       rtol,
                                       atol_j, xatol, slope1, slope2);
        if (notOK) {
            ndiff++;
            rel_diff = calc_rdiff((double) comp1[j], (double) comp2[j], rtol, atol_j);
            if (rel_diff > max_diff) {
                jmax = j;
                max_diff = rel_diff;
            }
            if (ndiff < 10 || (jmax == j)) {
                printf("\tfield variable %s at Node  %d ", varName.c_str(), j);
                printf(" differs: %g %g\n", comp1[j], comp2[j]);
            }
        }
    }

    if (nDataRowsMIN != nDataRowsMAX) {
        ndiff += nDataRowsMAX - nDataRowsMIN;
        if (ndiff < 10) {
            if (nDataRows1 > nDataRows2) {
                for (int j = nDataRowsMIN; j < nDataRowsMAX; j++) {
                    printf("\tColumn variable %s at data row %d ", varName.c_str(), j + 1);
                    printf(" differ: %g      NA\n", comp1[j]);
                }
            } else {
                for (int j = nDataRowsMIN; j < nDataRowsMAX; j++) {
                    printf("\tColumn variable %s at data row %d ", varName.c_str(), j + 1);
                    printf(" differ: NA     %g \n", comp2[j]);
                }
            }
        }
    }
    if (ndiff == 0) {
        printf("Field Variable %s passed test on %d Nodes\n", varName.c_str(), nDataRowsMIN);
    } else {
        printf("Field Variable %s failed test on %d Nodes out of %d Nodes\n", varName.c_str(), ndiff, nDataRowsMIN);
    }
    return (ndiff == 0);
}
//======================================================================================================================
bool compareVectorFlts(std::vector<double>& comp1,
                       std::vector<double>& comp2,
                       std::string varName)
{

    double rtol = rtolVar(varName);
    double atol = atolVar(varName);

    double max_diff = 0.0;
    int ndiff = 0;
    double rel_diff = 0.0;

    int nDataRows1 = comp1.size();
    int nDataRows2 = comp2.size();
    int nDataRowsMIN = MIN(nDataRows1, nDataRows2);
    int nDataRowsMAX = MAX(nDataRows1, nDataRows2);
    // Adjust for large values
    double atol_j = get_atol(comp1, nDataRows1, atol);
    atol_j = MIN(atol_j, get_atol(comp2, nDataRows2, atol));
    int jmax = 0;

    for (int j = 0; j < nDataRowsMIN; j++) {

        bool notOK = diff_double(comp1[j], comp2[j], rtol, atol_j);
        if (notOK) {
            ndiff++;
            rel_diff = calc_rdiff((double) comp1[j], (double) comp2[j], rtol, atol_j);
            if (rel_diff > max_diff) {
                jmax = j;
                max_diff = rel_diff;
            }
            if (ndiff < 10 || (jmax == j)) {
                printf("\tVariable %s at Node  %d ", varName.c_str(), j);
                printf(" differs: %g %g\n", comp1[j], comp2[j]);
            }
        }
    }

    if (nDataRowsMIN != nDataRowsMAX) {
        ndiff += nDataRowsMAX - nDataRowsMIN;
        if (ndiff < 10) {
            if (nDataRows1 > nDataRows2) {
                for (int j = nDataRowsMIN; j < nDataRowsMAX; j++) {
                    printf("\tvariable %s at data row %d ", varName.c_str(), j + 1);
                    printf(" differ: %g      NA\n", comp1[j]);
                }
            } else {
                for (int j = nDataRowsMIN; j < nDataRowsMAX; j++) {
                    printf("\tvariable %s at data row %d ", varName.c_str(), j + 1);
                    printf(" differ: NA     %g \n", comp2[j]);
                }
            }
        }
    }
    if (ndiff == 0) {
        printf("Variable %s passed test on %d Nodes\n", varName.c_str(), nDataRowsMIN);
    } else {
        printf("Variable %s failed test on %d Nodes out of %d Nodes\n", varName.c_str(), ndiff, nDataRowsMIN);
    }
    return (ndiff == 0);
}

//======================================================================================================================
XML_Node* getSimul(XML_Node* xmlTop, std::string id_tag)
{
    XML_Node* xctml = xmlTop;
    if (xctml->name() != "ctml") {
        xctml = xmlTop->findByName("ctml");
        if (!xctml) {
            throw CanteraError("countSimulations", "can't find ctml node");
        }
    }
    XML_Node* node = xctml->findNameID("simulation", id_tag);
    return node;
}

//======================================================================================================================
int countSimulations(XML_Node* xmlTop)
{
    XML_Node* xctml = xmlTop;
    if (xctml->name() != "ctml") {
        xctml = xmlTop->findByName("ctml");
        if (!xctml) {
            throw CanteraError("countSimulations", "can't find ctml node");
        }
    }
    std::vector<XML_Node*> ccc;
    xctml->getChildren("simulation", ccc);
    int sz = ccc.size();
    return sz;
}
//======================================================================================================================
XML_Node* getSimulNum(XML_Node* xmlTop, int num)
{
    XML_Node* xctml = xmlTop;
    if (xctml->name() != "ctml") {
        xctml = xmlTop->findByName("ctml");
        if (!xctml) {
            throw CanteraError("countSimulations", "can't find ctml node");
        }
    }
    std::vector<XML_Node*> ccc;
    xctml->getChildren("simulation", ccc);
    int sz = ccc.size();
    if (num < 0 || num >= sz) {
        throw CanteraError("getSimulNum", "out of bounds");
    }
    return ccc[num];
}
//======================================================================================================================
XML_Node* readXML(std::string inputFile)
{

    if (inputFile.size() == 0) {
        throw CanteraError("constructXML", "input file is null");
    }
    string path = findInputFile(inputFile);
    std::ifstream fin(path.c_str());
    if (!fin) {
        throw CanteraError("HMWSoln:constructPhaseFile", "could not open "
                           + path + " for reading.");
    }
    /*
     * The phase object automatically constructs an XML object.
     * Use this object to store information.
     */
    XML_Node* fxml = new XML_Node();
    fxml->build(fin);
    return fxml;
}

//! Define a spaces char vector, which includes everything returned by isspace()
const char* CHAR_SPACES = " \t\f\n\r\v";
const char* CHAR_DIGITS = "0123456789";
//======================================================================================================================
//! Returns the length of the next integer, not caring about what comes after the integer
/*!
 *   integer [spaces] [+|-]<unit>[GARBAGE]
 * @param  cstr
 * @return  returns the length of the integer including the
 */
int IntLen(const char* cstr)
{

    int k, n = 0;
    if (cstr) {
        n = strspn(cstr, CHAR_SPACES);
        cstr += n;

        if (*cstr == '-' || *cstr == '+') {
            ++cstr;
            ++n;
        }
        k = strspn(cstr, CHAR_DIGITS);
        n = k ? n + k : 0;
    }
    return n;
}
//========================================================================================================================
int isInteger(const char* estr)
{
    int retn = 0;
    char* cstr = (char*) malloc(strlen(estr) + 1);
    strcpy(cstr, estr);
    int ll = stripLTWScstring(cstr);

    int k, n = 0;
    if (cstr) {
        n = strspn(cstr, CHAR_SPACES);
        cstr += n;
        if (*cstr == '-' || *cstr == '+') {
            ++cstr;
            ++n;
        }
        k = strspn(cstr, CHAR_DIGITS);
        n = k ? n + k : 0;
    }
    if (n == ll) {
        retn = 1;
    }
    free(cstr);
    return retn;
}
//=====================================================================================================================

int isFloat(const char* estr)
{
    int retn = 0;
    char* cstr = (char*) malloc(strlen(estr) + 1);
    strcpy(cstr, estr);
    int ll = stripLTWScstring(cstr);
    int numDot = 0;
    int numExp = 0;
    char ch;
    int istart = 0;
    ch = estr[0];
    if (ch == '+' || ch == '-') {
        istart = 1;
        ch = estr[1];
    }
    if (!isdigit(ch)) {
        goto doneA;
    }
    for (int i = istart; i < ll; i++) {
        ch = estr[i];
        if (isdigit(ch)) {
        } else if (ch == '.') {
            numDot++;
            if (numDot > 1) {
                goto doneA;
            }
        } else if (ch == 'e' || ch == 'E' || ch == 'd' || ch == 'D') {
            numExp++;
            if (numExp > 1) {
                goto doneA;
            }
            ch = estr[i + 1];
            if (ch == '+' || ch == '-') {
                i++;
                ch = estr[i + 1];
            }
            if (!isdigit(ch)) {
                goto doneA;
            }
        } else {
            goto doneA;
        }
    }
    retn = 1;
doneA:
    free(cstr);
    return retn;
}
//======================================================================================================================
//! Identifies a string as a single string that can be a valid Cantera name
/*!
 *   Returns true if the c string is a valid single Cantera string name. A valid name doesn't contain any
 *   internal spaces, must start with an alphanumerical character, and doesn't contain any weird punctuation
 *
 * @param   estr  c string must not contain any internal spaces
 * @return  1 for yes and 0 for no
 */
int isStringName(const char* estr)
{
    int retn = 0;
    char ch;
    char* cstr = (char*) malloc(strlen(estr) + 1);
    strcpy(cstr, estr);
    int ll = stripLTWScstring(cstr);
    if (ll == 0) {
        goto doneA;
    }
    ch = cstr[0];
    if (!isalpha(ch)) {
        goto doneA;
    }
    for (int i = 1; i < ll; i++) {
        ch = estr[i];
        if (isalnum(ch)) {
        } else if (ispunct(ch)) {
            if (ch == '\'' || ch == '\"') {
                goto doneA;
            }
            if (ch == '\\') {
                goto doneA;
            }
            if (ch == ':' || ch == '`') {
                goto doneA;
            }
        } else {
            goto doneA;
        }
    }
    retn = 1;
doneA:
    free(cstr);
    return retn;
}
//======================================================================================================================
//! Returns a boolean indicating whether a c string has a decimal point or not
/*!
 *
 * @param c1    Character to be queried
 * @return      Returns true or false boolean
 */
int countDecimals(const char* c1)
{
    int index = 0;
    char* c2;
    if (c1 == NULL) {
        return 0;
    }
    do {
        c2 = strchr(c1, '.');
        if (c2) {
            index++;
            c1 = c2 + sizeof(char);
            if (!c1) {
                return index;
            }
        }
    } while (c2);
    return index;
}
//======================================================================================================================
//! Returns a boolean indicating whether a c string has a decimal point or not
/*!
 *
 * @param c1    Character to be queried
 * @return      Returns true or false boolean
 */
int countColons(const char* c1)
{
    int index = 0;
    char* c2 = 0;
    if (c1 == NULL) {
        return 0;
    }
    do {
        c2 = strchr(c1, ':');
        if (c2) {
            index++;
            c1 = c2 + sizeof(char);
            if (!c1) {
                return index;
            }
        }
    } while (c2);
    return index;
}
//======================================================================================================================
//! Returns a boolean indicating how many commas are in the c string
/*!
 *
 * @param c1    Character to be queried
 * @return      Returns true or false boolean
 */
int countCommas(const char* c1)
{
    int index = 0;
    char* c2 = 0;
    if (c1 == NULL) {
        return 0;
    }
    do {
        c2 = strchr(c1, ',');
        if (c2) {
            index++;
            c1 = c2 + sizeof(char);
            if (!c1) {
                return index;
            }
        }
    } while (c2);
    return index;
}
//======================================================================================================================
//!  Identify the contents of a field
/*!
 *
 * @param valueFieldc
 * @return
 */
enum valueType_enum identifyType(const std::string valueFieldc)
{
    valueType_enum retnT = STRING_VT;
    std::string valueField = stripws(valueField);
    int nS = valueField.size();
    int numColons;
    const char* delimColons = ": \t\n\f\r\v";
    const char* delimCommas = ", \t\n\f\r\v";

    std::string c1 = valueField;
    /*
     *  Count tokens
     */

    TKInput::TOKEN* t1 = new TKInput::TOKEN(c1.c_str());

    int ntokesWS = t1->ntokes;

    int numCommas = countCommas(c1.c_str());

    /*
     *  Handle the case where there is one token in the value field
     */
    if (ntokesWS == 1 && numCommas == 0) {
        numColons = countColons(t1->tok_ptrV[0]);
        if (numColons == 1) {
            TKInput::TOKEN* t2 = new TKInput::TOKEN(c1.c_str(), delimColons);
            if (!isStringName(t2->tok_ptrV[0])) {
                printf("identifyType: expected stringName in first colon field: \"%s\"", c1.c_str());
                return UNIDENTIFIED_VT;
            }
            if (!isFloat(t2->tok_ptrV[2])) {
                printf("identifyType: expected float in seond colon field: \"%s\"", c1.c_str());
                return UNIDENTIFIED_VT;
            }
            retnT = PAIR_STRING_FLOAT_CSV_VT;
            return retnT;
        }
        if (numColons == 2) {
            TKInput::TOKEN* t2 = new TKInput::TOKEN(c1.c_str(), delimColons);
            if (t2->ntokes != 3) {
                printf("identifyType: Unexpected value field: \"%s\"", c1.c_str());
                return UNIDENTIFIED_VT;
            }
            if (!isStringName(t2->tok_ptrV[0])) {
                printf("identifyType: expected stringName in first colon field: \"%s\"", c1.c_str());
                return UNIDENTIFIED_VT;
            }
            if (!isStringName(t2->tok_ptrV[1])) {
                printf("identifyType: expected stringName in second colon field: \"%s\"", c1.c_str());
                return UNIDENTIFIED_VT;
            }
            if (!isFloat(t2->tok_ptrV[2])) {
                printf("identifyType: expected float in third colon field: \"%s\"", c1.c_str());
                return UNIDENTIFIED_VT;
            }
            retnT = MATRIX_STRING_STRING_FLOAT_CSV_VT;
            return retnT;
        }
        if (isInteger(t1->tok_ptrV[0])) {
            return INTEGER_VT;
        }
        if (isFloat(t1->tok_ptrV[0])) {
            return FLOAT_VT;
        }
        if (isStringName(t1->tok_ptrV[0])) {
            return STRING_VT;
        }
    }

    TKInput::TOKEN* tc = new TKInput::TOKEN(c1.c_str(), delimCommas);





// int stokenize(c1.c_str(), const char *delimiters, char *tok_ptr[],
//                    const int max_tokens)

    /*
     * Count commas
     */

    /*
     * Count colons
     */

    return retnT;
}

//======================================================================================================================
class printRecord
{
public:
    printRecord(const XML_Node* node1, const XML_Node* node2, int etype = -1) :
        node1p_(node1),
        node2p_(node2),
        errorMessage_(""),
        eType_(etype),
        str1E(""),
        str2E("") {

    }

    void print() {
        string sbuf;
        printf("XML Nodes differ in their name");
        XML_Node* pp = node1p_->parent();
        sbuf = pp->name();
        printf(" first node with parent named %s", sbuf.c_str());
        //node1p_->write(std::cout, 2, 0);

    }

    const XML_Node* node1p_;
    const XML_Node* node2p_;
    std::string errorMessage_;

//! Error types
    /*!
     *  1  Name of the XML element's don't match
     *  2  Attribute values are different
     */
    int eType_;
    std::string str1E;
    std::string str2E;

};
//======================================================================================================================
std::vector<printRecord*> global_pr;
//======================================================================================================================

// Compare two nodes recursively
/*
 *
 */
int compareNodes(const XML_Node& node_1, const XML_Node& node_2)
{
    string sbuf;
    int nErrors = 0;
    /*
     *  Compare the name first
     */
    if (node_1.name() != node_2.name()) {
        printRecord* rec = new printRecord(&node_1, &node_2, 1);
        global_pr.push_back(rec);
        return 1;
    }
    const std::map<std::string, std::string>& attr1 = node_1.attribsConst();
    const std::map<std::string, std::string>& attr2 = node_2.attribsConst();
    int nAttr1 = attr1.size();
    int nAttr2 = attr2.size();
    int* nCov1 = new int[nAttr1];
    int* nCov2 = new int[nAttr2];

    std::map<std::string, std::string>::const_iterator it1;
    std::map<std::string, std::string>::const_iterator it2;
    int ii1 = 0;
    int ii2 = 0;
    for (it1 = attr1.begin(); it1 != attr1.end(); it1++) {
        bool iFound = false;
        ii2 = 0;
        for (it2 = attr2.begin(); it2 != attr2.end(); it2++) {
            if (nCov2[ii2] == 0) {
                if ((*it1).first == (*it2).first) {
                    iFound = true;
                    if ((*it1).second != (*it2).second) {
                        printRecord* rec = new printRecord(&node_1, &node_2, 2);
                        global_pr.push_back(rec);
                        rec->str1E = (*it1).first + " = " + (*it1).second;
                        rec->str2E = (*it2).first + " = " + (*it2).second;
                        nErrors++;
                    }
                    nCov2[ii2]++;
                    break;
                }
                ii2++;
            }
        }
        if (!iFound) {
            printRecord* rec = new printRecord(&node_1, &node_2, 2);
            global_pr.push_back(rec);
            rec->str1E = (*it1).first + " = " + (*it1).second;
            rec->str2E = (*it1).first + " = \"\"";
            nErrors++;
        }
    }
    ii1 = 0;
    ii2 = 0;
    for (it2 = attr2.begin(); it2 != attr2.end(); it2++) {
        bool iFound = false;
        ii1 = 0;
        for (it1 = attr1.begin(); it1 != attr1.end(); it1++) {
            if (nCov1[ii1] == 0) {
                if ((*it1).first == (*it2).first) {
                    iFound = true;
                    if ((*it1).second != (*it2).second) {
                        printRecord* rec = new printRecord(&node_1, &node_2, 2);
                        global_pr.push_back(rec);
                        rec->str1E = (*it1).first + " = " + (*it1).second;
                        rec->str2E = (*it2).first + " = " + (*it2).second;
                        nErrors++;
                    }
                    nCov1[ii1]++;
                    break;
                }
                ii1++;
            }
        }
        if (!iFound) {
            printRecord* rec = new printRecord(&node_1, &node_2, 2);
            global_pr.push_back(rec);
            rec->str1E = (*it2).first + " = \"\"";
            rec->str2E = (*it2).first + " = " + (*it2).second;
            nErrors++;
        }
    }

    /*
     *  Check the value
     */

    /*
     * Check the children
     */

    return nErrors;
}

//======================================================================================================================
static void print_usage()
{
    printf("\t\n");
    printf("  xmlSolnDiff [-h] [-a atol] [-r rtol] File1.xml File2.xml\n");
    printf("\t\n");
    printf("\tCompares two solution files in Cantera Solution XML format.\n");
    printf("\tThe comparison is done using a weighted norm basis.\n");
    printf("\t\n");
    printf("\tThe two files should be basically equal. However, File1.xml is\n");
    printf("\ttaken as the reference file that has precedence, when there is\n");
    printf("\tsomething to be decided upon.\n");
    printf("\t\n");
    printf("\t Arguments:\n");
    printf("\t  -h      = Usage info\n");
    printf("\t  -a atol = Set absolute tolerance parameter - default = 1.0E-9\n");
    printf("\t  -r rtol = Set relative tolerance parameter - default = 1.0E-3\n");
    printf("\t\n");
    printf("\t Shell Return Values:\n");
    printf("\t   1 = Comparison was successful\n");
    printf("\t   0 = One or more nodal values failed the comparison\n");
    printf("\t  -1 = Apples to oranges, the files can not even be compared against\n");
    printf("\t       one another.\n");
    printf("\t\n");
}
//====================================================================================================================
int main(int argc, char* argv[])
{
    int opt_let;
    std::string sbuf;
    int testFailed = 1;
    char* fileName1 = NULL, *fileName2 = NULL;
    XML_Node* topNode_1 = 0;
    XML_Node* topNode_2 = 0;
    XML_Node* ctmlNode_1 = 0;
    XML_Node* ctmlNode_2 = 0;
    int testPassed = 1;
    double atol_arg = 0.0, rtol_arg = 0.0;
    int id = 0;
    int id2 = 0;
    char* ggg = 0;
    char* rrr = 0;
    /*
     * Interpret command line arguments
     */
    /* Loop over each command line option */
    while ((opt_let = getopt(argc, argv, "ha:r:")) != EOF) {

        /* case over the option letter */
        switch (opt_let) {

        case 'h':
            /* Usage info was requested */
            print_usage();
            exit(0);

        case 'a':
            /* atol parameter */

            ggg = optarg;
            //printf("a = %s\n", ggg);
            id = sscanf(ggg, "%lg", &atol_arg);
            if (id != 1) {
                printf(" atol param bad: %s\n", ggg);
                exit(-1);
            }
            gatol = atol_arg;
            break;

        case 'r':
            /* rtol parameter */

            rrr = optarg;
            //printf("r = %s\n", ggg);
            id2 = sscanf(rrr, "%lg", &rtol_arg);
            if (id2 != 1) {
                printf(" rtol param bad: %s\n", rrr);
                exit(-1);
            }
            grtol = rtol_arg;
            break;

        default:
            /* Default case. Error on unknown argument. */
            printf("default called opt_let = %c\n", opt_let);
            fprintf(stderr, "ERROR in command line usuage:\n");
            print_usage();
            return 0;
        } /* End "switch(opt_let)" */

    } /* End "while((opt_let=getopt(argc, argv, "i")) != EOF)" */

    if (optind != argc - 2) {
        print_usage();
        exit(-1);
    } else {
        fileName1 = argv[argc - 2];
        fileName2 = argv[argc - 1];
    }

    /*
     *      Print Out Header
     */
    printf("\n");
    printf("----------------------------------------------------------------------\n");
    printf("ctmlDiff:    ctml Soln File comparison utility program for Cantera files\n");
    printf("             Version $Revision: 571 $\n");
    printf("             Harry K. Moffat Div. 1516 Sandia National Labs\n");
    printf("           \n");
    printf("             First  Cantera XML file = %s\n", fileName1);
    printf("             Second Cantera XML file = %s\n", fileName2);
    printf("\n");
    printf("             Absolute tol = %g\n", gatol);
    printf("             Relative tol = %g\n", grtol);
    printf("----------------------------------------------------------------------\n");
    printf("\n");

    /*
     *  Open up the two xml Files #1 and #2
     */

    if (!(topNode_1 = readXML(fileName1))) {
        fprintf(stderr, "Error opening up file1, %s\n", fileName1);
        exit(-1);
    }
    if (!(topNode_2 = readXML(fileName2))) {
        fprintf(stderr, "Error opening up file2, %s\n", fileName2);
        exit(-1);
    }
    if (topNode_1->name() != "ctml") {
        ctmlNode_1 = topNode_1->findByName("ctml");
        if (!ctmlNode_1) {
            sbuf = topNode_1->name();
            fprintf(stderr, "ctmlDiff:  %s with topname of  %s didn't have a ctml element\n", fileName1, sbuf.c_str());
            exit(-1);
        }
    } else {
        ctmlNode_1 = topNode_1;
    }
    if (topNode_2->name() != "ctml") {
        ctmlNode_2 = topNode_2->findByName("ctml");
        if (!ctmlNode_2) {
            sbuf = topNode_2->name();
            fprintf(stderr, "ctmlDiff:  %s with topname of  %s didn't have a ctml element\n", fileName2, sbuf.c_str());
            exit(-1);
        }
    } else {
        ctmlNode_2 = topNode_2;
    }

    int iTotal = compareNodes(*ctmlNode_1, *ctmlNode_2);

    if (iTotal) {
        printf(" One or more XML elements differ. Test failed\n");
        testFailed = 1;
    } else {
        printf("Test passed\n");
        testFailed = 0;
    }

    SAFE_DELETE(topNode_1);
    SAFE_DELETE(topNode_2);

    return testFailed;
}
//======================================================================================================================
