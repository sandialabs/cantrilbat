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
//! Returns the num'th SolnOutput record from the XML tree.
/*!
 *  @param[in]               xmlTop              XML Tree including a "ctml" XMLNode with multiple SolnOutput records
 *  @param[in]               num                 Number of the record
 */
static XML_Node* getSpeciesData(XML_Node* xmlTop, int num = 0)
{
    XML_Node* xctml = xmlTop;
    if (xctml->name() != "ctml") {
        xctml = xmlTop->findByName("ctml");
        if (!xctml) {
            throw ZuzaxError("getSimulNum","Can't find ctml node");
        }
    }
    std::vector<XML_Node*> ccc;
    xctml->getChildren("speciesData", ccc);
    int sz = ccc.size();
    if (num < 0 || num >= sz) {
        throw ZuzaxError("getSpeciesData", "out of bounds, requested %d record, only have %d records", num, sz);
    }
    return ccc[num];
}
//==================================================================================================================================
static size_t getSpeciesXML(XML_Node* xmlTop, std::vector<XML_Node*>& ccc)
{
    XML_Node* xctml = xmlTop;
    xctml->getChildren("species", ccc);
    size_t sz = ccc.size();
    return sz;
}

//==================================================================================================================================
int main(int argc, char* argv[])
{
    char buf[100000];
    const char*  fileName1=NULL;
    XML_Node* fp1 = nullptr;
    std::string optArgSS, arglc;
    std::vector<XML_Node*> ccc;
    fileName1 = "HMW_CuNaOCl_full.xml";
    double DG0, DH0, S0;
    double a1, a2, a3, a4, c1, c2, omega, charge;
    /*
     * Interpret command line arguments
     */
    if (!(fp1 = readXML(fileName1))) {
        fprintf(stderr,"Error opening up file1, %s\n", fileName1);
        exit(-1);
    }



    Zuzax::XML_Node* sd0 = getSpeciesData(fp1, 0);

    size_t numSp = getSpeciesXML(sd0, ccc);

    FILE *fcsv = fopen("species.csv", "w");    
    printf("Got  %d species points\n", (int) numSp);

    // Buffer the csv line into buf;
    for (size_t iS = 0; iS < numSp; iS++) {
        XML_Node* spNode = ccc[iS];
        buf[0] = '\0';
        int slen = 0;

        // species name
        std::string spName = spNode->attrib("name");
        if (spName == "") {
            throw ZuzaxError("xml", "blank species name");
        }
        sprintf(buf + slen, "%s , ", spName.c_str());
        slen = strlen(buf);

        XML_Node* thNode = & spNode->child("thermo");
        if (! thNode->hasChildOnlyOne("HKFT")) {
            continue;
        }
        XML_Node* hh = & thNode->child("HKFT");

        if (hh->hasChildOnlyOne("DG0_f_Pr_Tr")) {
           DG0 = ztml::getFloat(*hh, "DG0_f_Pr_Tr");
        }

        if (hh->hasChildOnlyOne("DH0_f_Pr_Tr")) {
            DH0 = ztml::getFloat(*hh, "DH0_f_Pr_Tr");
        }

        if (hh->hasChildOnlyOne("S0_Pr_Tr")) {
          S0 = ztml::getFloat(*hh, "S0_Pr_Tr");
        }

        sprintf(buf + slen, "% .2f , %.2f , %.2f ,", DG0 , DH0, S0);
        slen = strlen(buf);

        XML_Node* ss = & spNode->child("standardState");

        a1 = ztml::getFloat(*ss, "a1");
        a2 = ztml::getFloat(*ss, "a2");
        a3 = ztml::getFloat(*ss, "a3");
        a4 = ztml::getFloat(*ss, "a4");
        c1 = ztml::getFloat(*ss, "c1");
        c2 = ztml::getFloat(*ss, "c2");
        omega = ztml::getFloat(*ss, "omega_Pr_Tr");

        sprintf(buf + slen, "% .2f , %.2f , %.2f , %.2f , %.2f , %.2f , %.2f ,", 
                a1, a2, a3, a4, c1, c2, omega);
        slen = strlen(buf);

        charge =  ztml::getFloat(*spNode, "charge");
        std::string src;

        const XML_Node* xsource = spNode->findByName("source");
        src = "";
        if (xsource) {
           src = xsource->value();
        }


        sprintf(buf + slen, "% .2f , \"%s\" ", charge, "ref" );
        slen = strlen(buf);
        
        
        fprintf(fcsv,"%s\n", buf); 

    }

    fclose (fcsv);

    return 0;
}
//==================================================================================================================================
