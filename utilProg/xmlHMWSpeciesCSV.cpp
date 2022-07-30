/*
 *
 *  xmlHMWSpeciesCSV File1.xml
 *
 *  This utility program will print out the entries from all thermo XML blocks that
 *  have the HMW standard state formulation, in a comma separated variable format.
 *  All entries that make up the standard state formulation are printed out.
 *
 *  Arguments:
 *   -h = prints this usage information
 *
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

#include "zuzax/thermo/HMWSoln.h"

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


//===================================================================================================================================
// Gather the entropy of the elements of a species at 298 K. This is useful for going back and forth from the
// gibbs free energy of formation and the absolute gibbs free energy in NIST format.
double entropyElem298(Zuzax::ThermoPhase* g_ptr, size_t k)
{
    double se;
    double stotal = 0.0;
    for (size_t m = 0; m < g_ptr->nElements(); m++) {
        double na = g_ptr->nAtoms(k, m);
        if (na != 0.0) {
            se = g_ptr->entropyElement298(m, true);
            if (se == ENTROPY298_UNKNOWN) {
                return ENTROPY298_UNKNOWN;
            }
            stotal += se * na;
        }
    }
    double ch = g_ptr->charge(k);
    if (ch != 0.0) {
        double entCh = g_ptr->entropyCharge(ch);
        stotal += entCh;
    }

    return stotal;
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
static void printUsage()
{
    std::cout << "usage: xmlHMWSpeciesCVS [-h] [-u units] [-help_cmdfile] [fileName.xml]  "
              <<  std::endl;
    std::cout << "    -h           help" << endl;
    std::cout << "    -u   kJ or cals" << endl;
    std::cout << "    fileName.xml: XML File to parse for Standard State Specification block" << endl;
    std::cout << std::endl;
    std::cout << "        Prints out a CSV file consisting of the entries for the specification of the standard state\n" <<
              std::endl;
    std::cout << "        The file is named fileName_species.csv\n" << std::endl;
}

//==================================================================================================================================
int main(int argc, char* argv[])
{
    char buf[100000];
    XML_Node* fp1 = nullptr;
    std::string optArgSS, arglc;
    std::vector<XML_Node*> ccc;
    std::string fileName1 = "HMW_CuNaOCl_full.xml";
    std::string fileNameArg = "";

    int units = 0;
    double DG0, DH0, S0;
    double a1, a2, a3, a4, c1, c2, omega, charge;
    if (argc > 1) {
        std::string tok, tok2;
        for (int j = 1; j < argc; j++) {
            tok = std::string(argv[j]);
            if (tok[0] == '-') {
                int nopt = static_cast<int>(tok.size());
                for (int n = 1; n < nopt; n++) {
                    if (!strcmp(tok.c_str() + 1, "help_cmdfile")) {
                    } else if (tok[n] == 'h') {
                        printUsage();
                        exit(1);
                    } else if (tok[n] == 'u') {
                        j++;
                        tok2 = std::string(argv[j]);
                        if (!strcmp(tok2.c_str(), "kJ")) {
                            units = 1;
                        } else if (!strcmp(tok2.c_str(), "cals")) {
                            units = 0;
                        }
                        else {
                            printUsage();
                            exit(1);
                        }
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
    if (fileNameArg != "") {
        fileName1 = fileNameArg;
    }
    if (!(fp1 = readXML(fileName1))) {
        fprintf(stderr,"Error opening up file1, %s\n", fileName1.c_str());
        exit(-1);
    }

    std::string fileNameOut = fileName1;
    size_t sz = fileNameOut.size();
    fileNameOut = fileNameOut.substr(0,sz-4);
    fileNameOut += "_species.csv";


    Zuzax::XML_Node* sd0 = getSpeciesData(fp1, 0);

    size_t numSp = getSpeciesXML(sd0, ccc);

    FILE* fcsv = fopen(fileNameOut.c_str(), "w");
    printf("Got  %d species points\n", (int) numSp);

    std::string ssH = "    SpeciesName   ,";
    ssH += "      deltaG0 ,      deltaH0 ,    S0     ,";
    ssH += "      deltaS0 ,";
    ssH += "       a1  ,        a2 ,        a3 ,        a4 ,        c1 ,        c2  ,      omega,";
    ssH += " charge,  REFERENCE_FIELD";

    fprintf(fcsv,"%s\n", ssH.c_str());

    std::string ssU = "                  ,";
    if (units == 0) {
    ssU += "     cal/gmol ,  cal/gmol    , cal/gmolK ,";
    ssU += "   cal/gmolK , ";
    } else {
    ssU += "      kJ/gmol ,   kJ/gmol    ,  J/gmolK  ,";
    ssU += "   J/gmolK   , ";
    }
    ssU += "cal/gmol/bar, cal/gmol , cal K/gmol/bar, cal K/gmol, cal/gmol/K, cal K/gmol, cal/gmol,";
    ssU += "     ,               ";
    fprintf(fcsv,"%s\n", ssU.c_str());

    HMWSoln hmw(fileName1);
    double elemEntrop;

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
        size_t ik = hmw.speciesIndex(spName);

        if (ik != npos) {
            elemEntrop = entropyElem298(&hmw, ik);
        } else {
            elemEntrop = 0.0;
        }

        sprintf(buf + slen, " %16.16s , ", spName.c_str());
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
        if (units == 1) {
           DG0 *= 4.184 / 1.0E3;
           DH0 *= 4.184 / 1.0E3;
           S0 *= 4.184;
        }

        sprintf(buf + slen, "% 12.2f , % 12.2f , %12.3f ,", DG0 , DH0, S0);
        slen = strlen(buf);

        if (units == 1) {
        // Translate from J/kmol/K to J/gmol/K
          elemEntrop /= 1000;
        } else if (units == 0) {
        // Translate from J/kmol/K to cals/gmol/K
          elemEntrop *= 1.0 /( 1000 * 4.184);
        }
        double deltaS0 = S0 - elemEntrop;
        sprintf(buf + slen, "%12.3f ,", deltaS0);
        slen = strlen(buf);

        XML_Node* ss = & spNode->child("standardState");

        a1 = ztml::getFloat(*ss, "a1");
        a2 = ztml::getFloat(*ss, "a2");
        a3 = ztml::getFloat(*ss, "a3");
        a4 = ztml::getFloat(*ss, "a4");
        c1 = ztml::getFloat(*ss, "c1");
        c2 = ztml::getFloat(*ss, "c2");
        omega = ztml::getFloat(*ss, "omega_Pr_Tr");

        sprintf(buf + slen, "% 9.2f , %9.2f , %9.2f , %9.2f , %9.2f , % 10.2f , %9.2f ,",
                a1, a2, a3, a4, c1, c2, omega);
        slen = strlen(buf);

        charge =  ztml::getFloat(*spNode, "charge");
        std::string src;

        const XML_Node* xsource = spNode->findByName("source");
        src = "";
        if (xsource) {
            src = xsource->value();
        }


        sprintf(buf + slen, "% 6.2f , \"%s\" ", charge, "ref");
        slen = strlen(buf);


        fprintf(fcsv,"%s\n", buf);

    }

    fclose(fcsv);

    return 0;
}
//==================================================================================================================================
