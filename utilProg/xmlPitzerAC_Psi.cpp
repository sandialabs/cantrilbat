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

#include <stdio.h>
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
static int getPsiCommonCationParamXML(XML_Node* acNode, std::vector<XML_Node*>& ccc)
{
    acNode->getChildren("psiCommonCation", ccc);
    int sz = ccc.size();
    return sz;
}
//==================================================================================================================================
static int getPsiCommonAnionParamXML(XML_Node* acNode, std::vector<XML_Node*>& ccc)
{
    acNode->getChildren("psiCommonAnion", ccc);
    int sz = ccc.size();
    return sz;
}
//==================================================================================================================================
static void printUsage()
{
    std::cout << "usage: xmlPitzerAC_Psi [-h] [-help_cmdfile] [fileName]  "
              <<  std::endl;
    std::cout << "    -h           help" << endl;
    std::cout << "    fileName:   XML File to parse for psi blocks" << endl;
    std::cout << std::endl;
}
//==================================================================================================================================
int main(int argc, char* argv[])
{
    char buf[100000];
    XML_Node* xmlTop = nullptr;
    std::string optArgSS, arglc;
    std::vector<XML_Node*> ccc;
    std::string src;
    std::string fileName1 = "HMW_CuNaOCl_full.xml";
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

    FILE* fcsv = fopen("PsiParam.csv", "w");

    // Buffer the csv line into buf;
    for (size_t ii = 0; ii < 2; ++ii) {
        std::string itype = "psiCommonCation";
        if (ii == 0) {
            numParamBlocks = getPsiCommonCationParamXML(acNode, ccc);
            printf("Got  %d ternary  psi common cation parameter points\n", (int) numParamBlocks);
        } else if (ii == 1) {
            itype = "thetaAnion";
            numParamBlocks = getPsiCommonAnionParamXML(acNode, ccc);
            printf("Got  %d ternary psi common anion parameter points\n", (int) numParamBlocks);
        }
        for (size_t iP = 0; iP < numParamBlocks; iP++) {
            XML_Node* bspNode = ccc[iP];
            buf[0] = '\0';
            int slen = 0;

            // cation species name
            std::string spName1, spName2, spName3;
            if (ii == 0) {
                spName1 = bspNode->attrib("cation");
                spName2 = bspNode->attrib("anion1");
                spName3 = bspNode->attrib("anion2");
            } else if (ii == 1) {
                spName1 = bspNode->attrib("anion");
                spName2 = bspNode->attrib("cation1");
                spName3 = bspNode->attrib("cation2");
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

            sprintf(buf + slen, "%s , ", spName3.c_str());
            slen = strlen(buf);

            std::vector<double> psiB;
            psiB.resize(0);
            int num;
            // check for all lower case on theta!
            if (bspNode->hasChildOnlyOne("psi")) {
                num = ztml::getFloatArray(*bspNode, psiB, "", "", "psi");
                if ((int) psiB.size() != num) {
                    printf("Psi block for  %s  %s %s  return and size differs: %d %d\n", spName1.c_str(), spName2.c_str(),
                             spName3.c_str(), num, (int) psiB.size());

                }
            } else {
                printf("Psi block for  %s  %s %s not there\n", spName1.c_str(), spName2.c_str(), spName3.c_str());
                exit(-1) ;
            }
            if (psiB.size() != 5) {
                printf("Psi block for  %s  %s %s only has %d entries\n", spName1.c_str(), spName2.c_str(), spName3.c_str(), num);
                exit(-1) ;
            }
            if (psiB[2] != 0.0) {
                printf("Psi[2] != 0.0: %g\n", psiB[2]);
                exit(-1);
            }
            sprintf(buf + slen, "% .2f , %.2f , %.2f , %.2f , ", psiB[0], psiB[1], psiB[3], psiB[4]);
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
