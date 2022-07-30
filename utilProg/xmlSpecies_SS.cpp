/*
 *
 *          xmlHMWSpeciesCVS [-h] [-help_cmdfile] [fileName]  "
 *
 *    Creates an CSV file consisting of the key entries for the specification of the standard state
 *    The analysis file is named speciesSS.csv 
 *
 *  Arguments:
 *   -h = prints this usage information
 *    fileName   selects that XML file for analysis
 *
 *  Shell Return Values
 *    0 = Analysis was successful
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

#include "zuzax/thermo.h"

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
static void printUsage()
{
    std::cout << "usage: xmlHMWSpeciesCVS [-h] [-help_cmdfile] [fileName.xml]  "
              <<  std::endl;
    std::cout << "    -h           help" << endl;
    std::cout << "    fileName:   XML File to parse for Standard State Specification block" << endl;
    std::cout << std::endl;
    std::cout << "        Creates out a CSV file consisting of the key entries for the specification of the standard state\n"; 
    std::cout << "        The file is named fileName_speciesSS.csv \n" << std::endl;
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
    fileNameOut += "_speciesSS.csv";
    

    thermo_t_double& tp = *newPhase(fileName1);
    size_t nsp = tp.nSpecies();
    std::vector<double> xmol(tp.nSpecies());
    tp.getMoleFractions(xmol.data());
    double T298 = 298.15;
    tp.setState_TPX(T298, OneBar, xmol.data());

    std::vector<double>  h_ref(nsp);
    tp.getEnthalpy_RT_ref(h_ref.data());

    std::vector<double> s_ref(nsp), g_ref(nsp);
    tp.getEntropy_R_ref(s_ref.data());
    tp.getGibbs_RT_ref(g_ref.data());

    std::vector<double> u_ref(nsp), cp_ref(nsp);

    tp.getIntEnergy_RT(u_ref.data());

    tp.getCp_R_ref(cp_ref.data());

     std::vector<double> mv_ref(nsp);
    tp.getStandardVolumes_ref(mv_ref.data());

    double RT = GasConstant * 298.15;
    for (size_t k = 0; k < nsp; k++) {
        h_ref[k] *= RT;
        g_ref[k] *= RT;
        s_ref[k] *= GasConstant; 
        cp_ref[k] *= GasConstant; 
        u_ref[k] *= RT;
        if (fabs(h_ref[k]) < 1.0E-12) {
            h_ref[k] = 0.0;
        }
        if (fabs(g_ref[k]) < 1.0E-12) {
            g_ref[k] = 0.0;
        }
        if (fabs(s_ref[k]) < 1.0E-12) {
            s_ref[k] = 0.0;
        }
        if (fabs(u_ref[k]) < 1.0E-12) {
            u_ref[k] = 0.0;
        }
        if (fabs(cp_ref[k]) < 1.0E-12) {
            cp_ref[k] = 0.0;
        }
        if (fabs(mv_ref[k]) < 1.0E-12) {
            mv_ref[k] = 0.0;
        }
    }

    std::vector<double> hf_298(nsp);
    for (size_t k = 0; k < nsp; k++) {
        hf_298[k] = tp.Hf298SS(k);
        if (fabs(hf_298[k]) < 1.0E-12) {
           hf_298[k] = 0.0;
        }
    }
    std::vector<double> gf_298(nsp);
    for (size_t k = 0; k < nsp; k++) {
        gf_298[k] = tp.Gf298Ref(k);
        if (fabs(gf_298[k]) < 1.0E-12) {
           gf_298[k] = 0.0;
        }
    }

    FILE *fcsv = fopen(fileNameOut.c_str(), "w");    
    printf("Printed out  %d species 298.15K reference state values at 1 bar\n", (int) nsp);

    fprintf(fcsv, "    SpeciesName  ,      HF_298  ,      S_ref  ,      G_ref  , delta_Gf_298,     Cp_ref  , molarVol_ref\n");
    fprintf(fcsv, "                 ,     kJ/gmol  ,    J/gmol/K ,    kJ/gmol  ,     kJ/gmol ,   J/gmol/K  ,    cm3/gmol \n");

    // Buffer the csv line into buf;
    for (size_t k = 0; k < nsp; k++) {
        buf[0] = '\0';
        int slen = 0;

        // species name
        std::string spName = tp.speciesName(k);
        sprintf(buf + slen, "%16.16s , ", spName.c_str());
        slen = strlen(buf);

        sprintf(buf + slen, " % 12.5E, % 12.5E,", 1.0E-6*hf_298[k], 1.0E-3*s_ref[k]);
        slen = strlen(buf);

        sprintf(buf + slen, " % 12.5E, % 12.5E, % 12.5E,", 1.0E-6*g_ref[k], 1.0E-6*gf_298[k], 1.0E-3*cp_ref[k]);
        slen = strlen(buf);

        sprintf(buf + slen, " % 12.5E, ", 1.0E3*mv_ref[k]);
        slen = strlen(buf);
        
        fprintf(fcsv,"%s\n", buf); 

    }

    fclose (fcsv);

    return 0;
}
//==================================================================================================================================
