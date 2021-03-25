/*
 *
 *  xmlSolnDiff File1.xml File2.xml
 *
 *  Compares the variable values in two Zuzax solution xml
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

#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <fstream>
#include <unistd.h>
#include <iostream>


#include "zuzax/base/xml.h"
#include "zuzax/base/ctml.h"
#include "zuzax/base/ctexceptions.h"
#include "zuzax/base/zzcompare.h"
#include "zuzax/base/VarType.h"

#include "zuzax/numerics/SolnOutput.h"
#include "zuzax/numerics/SolnState.h"
#include "zuzax/numerics/TimeSpec.h"
#include "zuzax/numerics/GridSpec.h"
#include "zuzax/base/Distrib1D.h"

#include <cstdio>
#include <cstring>


using namespace std;
using namespace Zuzax;
//==================================================================================================================================
//! Read an XML file into a XML_Node Tree structure
/*!
 *  @param[in]               inputFile           Name of the input file
 */
static XML_Node* readXML(std::string inputFile)
{

    if (inputFile.size() == 0) {
        throw ZuzaxError("readXML()",  "input file is null");
    }
    std::string path = findInputFile(inputFile);
    std::ifstream fin(path.c_str());
    if (!fin) {
        throw ZuzaxError("readXML()","could not open " +path+" for reading.");
    }
    XML_Node* fxml = new XML_Node();
    fxml->build(fin);
    return fxml;
}

//==================================================================================================================================
static XML_Node* getPhaseXML(XML_Node* xmlTop, int num = 0)
{
    XML_Node* xctml = xmlTop;
    if (xctml->name() != "ctml") {
        xctml = xmlTop->findByName("ctml");
        if (!xctml) {
            throw ZuzaxError("getSimulNum","Can't find ctml node");
        }
    }

    XML_Node* pp;
    if (xctml->hasChild("phase")) {
        pp = & xctml->child("phase");
    } else {
        throw ZuzaxError("readXML()","XML file doesn't have a phase in it");
    }
    return pp;
}
//==================================================================================================================================
static XML_Node* getThermoXML(XML_Node* xmlPhase, int num = 0)
{
    XML_Node* tt = nullptr;
    if (xmlPhase->hasChild("thermo")) {
        tt = & xmlPhase->child("thermo");
    } else {
        throw ZuzaxError("readXML()","XML file doesn't have a phase in it");
    }
    return tt;
}
//==================================================================================================================================
static XML_Node* getACXML(XML_Node* ttPhase, int num = 0)
{
    XML_Node* acNode = nullptr;
    if (ttPhase->hasChild("activityCoefficients")) {
        acNode = & ttPhase->child("activityCoefficients");
    } else {
        throw ZuzaxError("getACXML()","XML file doesn't have a phase in it");
    }
    return acNode;
}
//==================================================================================================================================
static int getLambdaParamXML(XML_Node* acNode, std::vector<XML_Node*>& ccc)
{
    acNode->getChildren("lambdaNeutral", ccc);
    int sz = ccc.size();
    return sz;
}

//==================================================================================================================================
static void printUsage()
{
    std::cout << "usage: xmlPitzerAC_Lambda [-h] [-help_cmdfile] [fileName]  "
              <<  std::endl;
    std::cout << "    -h           help" << endl;
    std::cout << "    fileName:   XML File to parse for lambdaNeutral blocks" << endl;
    std::cout << std::endl;
}
//==================================================================================================================================
int main(int argc, char* argv[])
{
    char buf[100000];
    std::string  fileName1 = "HMW_CuNaOCl_full.xml";
    XML_Node* xmlTop = nullptr;
    std::string optArgSS, arglc;
    std::vector<XML_Node*> ccc;
    std::string src;
    fileName1 = "HMW_CuNaOCl_full.xml";
    std::string fileNameArg = "";

    if (argc > 1) {
        std::string tok;
        for (int j = 1; j < argc; j++) {
            tok = std::string(argv[j]);
            if (tok[0] == '-') {
                int nopt = static_cast<int>(tok.size());
                for (int n = 1; n < nopt; n++) {
                    if (!strcmp(tok.c_str() + 1, "help_cmdfile")) {
                    } else if (tok[n] == 'h') {
                        printUsage();
                        exit(1);
                    }
                }
            } else {
                if (fileNameArg != "") {
                    printUsage();
                    exit(1);
                }
                fileNameArg = tok;
            }
        }
    }
    if (fileNameArg != "") {
        fileName1 = fileNameArg;
    }

    /*
     * Interpret command line arguments
     */
    if (!(xmlTop = readXML(fileName1))) {
        fprintf(stderr,"Error opening up file1, %s\n", fileName1.c_str());
        exit(-1);
    }


    Zuzax::XML_Node* xmlPhase = getPhaseXML(xmlTop);

    Zuzax::XML_Node* ttPhase = getThermoXML(xmlPhase);

    Zuzax::XML_Node* acNode = getACXML(ttPhase);

    size_t numParamBlocks;

    FILE* fcsv = fopen("LambdaParam.csv", "w");

    // Buffer the csv line into buf;
    for (size_t ii = 0; ii < 1; ++ii) {
        std::string itype = "lambdaNeutral";
        if (ii == 0) {
            numParamBlocks = getLambdaParamXML(acNode, ccc);
            printf("Got  %d binary lambda parameter points\n", (int) numParamBlocks);
        }
        for (size_t iP = 0; iP < numParamBlocks; iP++) {
            XML_Node* bspNode = ccc[iP];
            buf[0] = '\0';
            int slen = 0;

            // neutral , ion or neutral  species name

            std::string spName1, spName2;
            if (ii == 0) {
                spName1 = bspNode->attrib("species1");
                spName2 = bspNode->attrib("species2");
            }
            if (spName1 == "") {
                throw ZuzaxError("xml", "blank species name");
            }
            if (spName2 == "") {
                throw ZuzaxError("xml", "blank species name");
            }
            sprintf(buf + slen, "%s , ", spName1.c_str());
            slen = strlen(buf);

            sprintf(buf + slen, "%s , ", spName2.c_str());
            slen = strlen(buf);


            std::vector<double> lambdaB;
            lambdaB.resize(0);
            int num;
            // check for all lower case on theta!
            if (bspNode->hasChildOnlyOne("lambda")) {
                num = ztml::getFloatArray(*bspNode, lambdaB, "", "", "lambda");
                if ((int) lambdaB.size() != num) {
                    printf("Lambda block for  %s  %s  return and size differs: %d %d\n", spName1.c_str(), spName2.c_str(), num,
                           (int) lambdaB.size());
                }
            } else {
                printf("Lambda block for  %s  %s  not there\n", spName1.c_str(), spName2.c_str());
                exit(-1) ;
            }
            if (lambdaB.size() != 5) {
                printf("Lambda block for  %s  %s  only has %d entries\n", spName1.c_str(), spName2.c_str(), num);
                exit(-1) ;
            }
            // In everything that we have been doing, the quadratic term has been zero.
            if (lambdaB[2] != 0.0) {
                printf("lambda[2] != 0.0: %g\n", lambdaB[2]);
                exit(-1);
            }
            sprintf(buf + slen, "% .2f , %.2f , %.2f , %.2f , ", lambdaB[0], lambdaB[1], lambdaB[3], lambdaB[4]);
            slen = strlen(buf);


            const XML_Node* xsource = bspNode->findByName("source");
            src = "";
            if (xsource) {
                src = xsource->value();
            }
            sprintf(buf + slen, " \"%s\" ",  src.c_str());
            slen = strlen(buf);

            fprintf(fcsv,"%s\n", buf);
        }
    }

    fclose(fcsv);

    return 0;
}

//==================================================================================================================================
