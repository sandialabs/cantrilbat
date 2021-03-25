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
    if (xmlPhase->hasChild("thermo") ) {
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
static int getBinarySaltParamXML(XML_Node* acNode, std::vector<XML_Node*>& ccc)
{
    acNode->getChildren("binarySaltParameters", ccc);
    int sz = ccc.size();
    return sz;
}

//==================================================================================================================================
int main(int argc, char* argv[])
{
    char buf[100000];
    const char*  fileName1=NULL;
    XML_Node* xmlTop = nullptr;
    std::string optArgSS, arglc;
    std::vector<XML_Node*> ccc;
    std::string src;
    fileName1 = "HMW_CuNaOCl_full.xml";
    double alpha1 = 2.0, alpha2 = 12.0;
    /*
     * Interpret command line arguments
     */
    if (!(xmlTop = readXML(fileName1))) {
        fprintf(stderr,"Error opening up file1, %s\n", fileName1);
        exit(-1);
    }


    Zuzax::XML_Node* xmlPhase = getPhaseXML(xmlTop);

    Zuzax::XML_Node* ttPhase = getThermoXML(xmlPhase);

    Zuzax::XML_Node* acNode = getACXML(ttPhase);

    size_t numParamBlocks = getBinarySaltParamXML(acNode, ccc);

    FILE *fcsv = fopen("BinarySaltParam.csv", "w");    
    printf("Got  %d binary salt parameter points\n", (int) numParamBlocks);

    // Buffer the csv line into buf;
    for (size_t iP = 0; iP < numParamBlocks; iP++) {
        XML_Node* bspNode = ccc[iP];
        buf[0] = '\0';
        int slen = 0;

        // cation species name
        std::string spCatName = bspNode->attrib("cation");
        if (spCatName == "") {
            throw ZuzaxError("xml", "blank species name");
        }
        sprintf(buf + slen, "%s , ", spCatName.c_str());
        slen = strlen(buf);

        // anion species name
        std::string spAnionName = bspNode->attrib("anion");
        if (spAnionName == "") {
            throw ZuzaxError("xml", "blank species name");
        }
        sprintf(buf + slen, "%s , ", spAnionName.c_str());
        slen = strlen(buf);

       alpha1 = 2.0;
       if (bspNode->hasChild("Alpha1")) {
           alpha1 = ztml::getFloat(*bspNode, "Alpha1");
       }

       alpha2 = 12.0;
       if (bspNode->hasChild("Alpha2")) {
           alpha2 = ztml::getFloat(*bspNode, "Alpha2");
       }
       sprintf(buf + slen, "% .2f , %.2f , ", alpha1 , alpha2);
       slen = strlen(buf);

       std::vector<double> beta0;
       int num;
       if (bspNode->hasChildOnlyOne("beta0")) {
          num = ztml::getFloatArray(*bspNode, beta0, "", "", "beta0");
       }
       if (beta0[2] != 0.0) {
           printf("beta0[2] != 0.0: %g\n", beta0[2]);
           exit(-1);
       }
       sprintf(buf + slen, "% .2f , %.2f , %.2f , %.2f , ", beta0[0], beta0[1], beta0[3], beta0[4]);
       slen = strlen(buf);

       std::vector<double> beta1;
       num = 5;
       if (bspNode->hasChildOnlyOne("beta1")) {
          num = ztml::getFloatArray(*bspNode, beta1, "", "", "beta1");
       }
       if (beta1[2] != 0.0) {
           printf("beta1[2] != 0.0: %g\n", beta1[2]);
           exit(-1);
       }
       sprintf(buf + slen, "% .2f , %.2f , %.2f , %.2f , ", beta1[0], beta1[1], beta1[3], beta1[4]);
       slen = strlen(buf);

       std::vector<double> beta2;
       num = 5;
       if (bspNode->hasChildOnlyOne("beta2")) {
          num = ztml::getFloatArray(*bspNode, beta2, "", "", "beta2");
       }
       if (beta2[2] != 0.0) {
           printf("beta2[2] != 0.0: %g\n", beta2[2]);
           exit(-1);
       }
       sprintf(buf + slen, "% .2f , %.2f ,  %.2f , %.2f , ", beta2[0], beta2[1], beta2[3], beta2[4]);
       slen = strlen(buf);

       std::vector<double> Cphi;
       num = 5;
       if (bspNode->hasChildOnlyOne("Cphi")) {
          num = ztml::getFloatArray(*bspNode, Cphi, "", "", "Cphi");
       }
       if (Cphi[2] != 0.0) {
           printf("Cphi[2] != 0.0: %g\n", Cphi[2]);
           exit(-1);
       }
       sprintf(buf + slen, "% .2f , %.2f , %.2f , %.2f , ", Cphi[0], Cphi[1], Cphi[3], Cphi[4]);
       slen = strlen(buf);


       const XML_Node* xsource = bspNode->findByName("source");
       src = "";
       if (xsource) {
          src = xsource->value();
       }
       sprintf(buf + slen, " \"%s\" ",  src.c_str() );
       slen = strlen(buf);
        
       fprintf(fcsv,"%s\n", buf); 

    }

    fclose (fcsv);

    return 0;
}
//==================================================================================================================================
