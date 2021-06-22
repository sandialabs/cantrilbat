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

#include "zuzax/thermo/HMWSoln.h"
#include "zuzax/thermo/MetalSHEelectrons.h"
#include "zuzax/thermo/StoichSubstance.h"


#include "zuzax/kinetics/InterfaceKinetics.h"
#include "zuzax/kinetics/ElectrodeKinetics.h"
#include "zuzax/kinetics/RxnMolChange.h"
#include "zuzax/kinetics.h"

#include "zuzax/equil/vcs_MultiPhaseEquil.h"

bool doAnode = false;

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

InterfaceKinetics* gKinetics = nullptr;
HMWSoln *brine = nullptr;
StoichSubstance* cuMetal = nullptr;
MetalSHEelectrons* eeSHE = nullptr;
StoichSubstance* Cu2O_S = nullptr;

//==================================================================================================================================
// handle to initialize all of the phases
static void initializeInterfaceRxn( std::vector<thermo_t_double*>& thAnodeVec)
{
    ElectrodeKinetics* cathodeKin = nullptr;
    thAnodeVec.clear();
    brine = new HMWSoln("HMW_CuNaOCl_full.xml");

    brine->setState_TPM(298.15, OneBar, "Cu++:0.00001 Cu+:0.00001 Na+:4.0 Cl-:4.00003, O2(aq):1.0E-5 H+:1.0E-8 OH-:1.0E-8");
 
    brine->setElectricPotential(-0.40);

    cuMetal = new StoichSubstance("CuSolid.xml");
    eeSHE = new MetalSHEelectrons("metalSHEelectrons.xml");
    Cu2O_S = new StoichSubstance("Cu2O_solid.xml");


    thAnodeVec.push_back(brine);
    thAnodeVec.push_back(cuMetal);
    thAnodeVec.push_back(eeSHE);
    thAnodeVec.push_back(Cu2O_S);

    if (doAnode) {
        cathodeKin = (ElectrodeKinetics*) newInterfaceKineticsMgrFromFile(thAnodeVec, "Cu_electrode.xml");
    } else {
        cathodeKin = 
          (ElectrodeKinetics*) newInterfaceKineticsMgrFromFile(thAnodeVec, "OH_cathode_damk.xml");
    }

    gKinetics = cathodeKin;

}
//==================================================================================================================================
// handle to set the bath conditions for the analysis
static void setBathConditions(std::vector<thermo_t_double*>& thermo_vec, InterfaceKinetics* gKinetics)
{

      gKinetics->setStateKin_TP(298.15, OneBar);
      brine->addMolalitiesByName(" Cu+:5.5E-2 CuCl2-:5.5E-2");
      brine->addMolalitiesByName(" Cu+:1.0 OH-:1.0");
      brine->addMolalitiesByName(" Cu++:1.0 OH-:2.0");
      int nsp = brine->nSpecies();
      if (nsp != 21) {
        printf("error\n");
        exit(-1);
      }
      double Phi_lyte;

      // dump at time = 40053 seconds
      double x[30];
      if (doAnode) {
      x[0] = 8.3221E-01 ;
      x[1] = 2.1241E-03 ;
      x[2] = 4.4922E-13 ;
      x[3] = 1.0156E-01;
      x[4] =  8.5796E-06;
      x[5] =  5.0077E-11 ;
      x[6] =  1.0069E-25;
      x[7] =  1.2043E-05;
      x[8] =   1.6896E-10;
      x[9] = 4.5062E-13;
      x[10] = 4.2135E-07;
      x[11] = 5.5681E-05 ;
      x[12] = 2.6480E-02;
      x[13] =  1.2366E-15;
      x[14] =  2.0402E-08;
      x[15] =   1.5479E-17 ;
      x[16] =  3.5418E-02 ;
      x[17] =   2.1250E-03 ;
      x[18] =   2.4134E-10;
      x[19] = 6.0858E-10;
      x[20] =  4.4503E-07 ;
      Phi_lyte = -0.291692;
      } else {
         Phi_lyte = -0.400308;
         x[0] = 8.9927E-1;
         x[1] = 3.3606E-2 ;
         x[2] = 2.0593E-16 ;
         x[3] = 5.0367E-2;
         x[4] = 1.6757E-2 ;
         x[5] = 1.7992E-20 ;
         x[6] = 8.7839E-25;
         x[7] = 1.7887E-12 ;
         x[8] = 1.6985E-17 ;
         x[9] = 1.7584E-17;
         x[10] = 5.5571E-11;
         x[11] = 2.6981E-11;
         x[12] = 9.5548E-8;
         x[13] = 6.8062E-19;
         x[14] = 1.8432E-9;
         x[15] = 6.3853E-20;
         x[16] = 1.9143E-6;
         x[17] = 1.8162E-7 ;
         x[18] = 2.9582E-17;
         x[19] = 7.2343E-8;
         x[20] = 1.5158E-6;
      }
      brine->setMoleFractions(x);

      double molal[100];

      brine->setElectricPotential(Phi_lyte);

      // Equilibrate the brine system. -> note this doesn't equilibrate Cu++ and Cu- wrt electrode
      //vcs_equilibrate(*brine, "TP");

      brine->getMolalities(molal);

      printf(" setBathConditions: New Molalities Vector: \n\n");
      for (size_t k = 0; k < brine->nSpecies() ; ++k) {
         printf("\t %3d , %24s , %12.4E \n", (int) k, brine->speciesName(k).c_str(), molal[k]);
      }
      printf("\n\n\n");

}
double af, bf, Ea, io_af, io_bf, io_Ea;
//==================================================================================================================================
static void calculateBVForm( size_t iR  )
{

    // get forward rate constant

    

     // convert

     //gKinetics->convertExchangeCurrentDensityFormulation(doublevalue* const kfwd);


     gKinetics->updateROP();

     gKinetics->getBVFormKFrwd(iR, af, bf, Ea, io_af, io_bf, io_Ea);

    // First see if it is a BV form

    


}
//==================================================================================================================================
// handle to set the bath conditions for the analysis
static void  cleanup(std::vector<thermo_t_double*>& thermo_vec, InterfaceKinetics* gKinetics)
     // Ea is in temperature
{
   delete gKinetics;
   delete brine;
   delete cuMetal;
   delete eeSHE;
   delete Cu2O_S;

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
    /*
     * Interpret command line arguments
     */
    if (!(xmlTop = readXML(fileName1))) {
        fprintf(stderr,"Error opening up file1, %s\n", fileName1);
        exit(-1);
    }

    std::vector<thermo_t_double*> thermo_vec;

    initializeInterfaceRxn(thermo_vec);

    size_t nr = gKinetics->nReactions();
    printf("Number of interfacial reactions = %d \n", (int) nr);

    ElectrodeKinetics* eK = dynamic_cast<ElectrodeKinetics*>(gKinetics);

    if (!eK) {
        printf("Interface isn't an electrode interface\n");
        return -1;
    }

    FILE *fcsv = fopen("InterfaceReactions.csv", "w");    

    setBathConditions(thermo_vec, gKinetics);

    // Note that the ReactionData class is not retained after instantiation of the kinetics class.
    // Since we are doing calculations with the InterfaceKinetics class, we need to not use the ReactionData
    // class.

    double nStoichElectrons, OCV, io, overPotential, beta, resistance, damk;
    double deltaG[40], deltaSSG[40];
    // Buffer the csv line into buf;

    printf("Species Names\n");
    for( size_t kKin = 0; kKin < gKinetics->nKinSpecies(); kKin++) {
        printf("\t %3d  %s\n", (int) kKin, gKinetics->kineticsSpeciesName(kKin).c_str());
    }

    fprintf(fcsv, " EqnString , ");
    fprintf(fcsv, " chargeTransfer? , ");
    fprintf(fcsv, " ImposedVoltage , ");
    fprintf(fcsv, " Ess , ");
    fprintf(fcsv, " E , ");
    fprintf(fcsv, " OCV , ");
    fprintf(fcsv, " overPotential , ");
    fprintf(fcsv, " io , ");
    fprintf(fcsv, " beta , ");
    fprintf(fcsv, " nStoichElectrons , ");
    fprintf(fcsv, " af ,   bf ,  Ea_kjgm,   io_af,   io_bf, io_Ea_kjgm , damk");

    double TT = gKinetics->reactionPhaseThermo().temperature();

    fprintf(fcsv, "\n");

    for (size_t iR = 0; iR < nr; iR++) {
        buf[0] = '\0';
        int slen = 0;
        RxnMolChange* rmc = new RxnMolChange(gKinetics, iR);
        int ct = rmc->m_ChargeTransferInRxn;

        size_t iMetal = eK->metalPhaseIndex();

        size_t iSoln = eK->solnPhaseIndex();


        double phi0Metal = (eK->thermo(iMetal)).electricPotential();
        double phi0Soln = (eK->thermo(iSoln)).electricPotential();
        double V0 = phi0Metal - phi0Soln;

        eK->getExchangeCurrentDensityFormulationD(iR, nStoichElectrons, OCV, io, overPotential, beta, resistance, damk);
        eK->getDeltaGibbs(DATA_PTR(deltaG));
        eK->getDeltaSSGibbs(DATA_PTR(deltaSSG));
        double deltaGrxn = deltaG[iR];
        double deltaSSGrxn = deltaSSG[iR];
        double Erxn =  deltaGrxn/Faraday/nStoichElectrons;
        double ESSrxn =  deltaSSGrxn/Faraday/nStoichElectrons;

        double ii = eK->calcCurrentDensity(overPotential, nStoichElectrons,io, beta, TT, 0.0, damk);

        printf(" calculated current density = %g A/m2\n", ii);


        calculateBVForm(iR);
    
        // Convert to kJ/gm for printing out (and only printing out)
        double Ea_kjgm = Ea * 1.0E-6;
        double io_Ea_kjgm = io_Ea * 1.0E-6;

        // Calculate the equation string.
        const std::string& rS = gKinetics->reactionString(iR);
        sprintf(buf + slen, " %s , " , rS.c_str());
        slen = strlen(buf);

        sprintf(buf + slen, " %d , " , ct );
        slen = strlen(buf);

        sprintf(buf + slen, " %g , " , V0);
        slen = strlen(buf);

        sprintf(buf + slen, " %g , " , ESSrxn);
        slen = strlen(buf);

        sprintf(buf + slen, " %g , " , Erxn);
        slen = strlen(buf);

        sprintf(buf + slen, " %g , " , OCV);
        slen = strlen(buf);

        sprintf(buf + slen, " %g , " , overPotential);
        slen = strlen(buf);

        sprintf(buf + slen, " %g , " , io);
        slen = strlen(buf);

        sprintf(buf + slen, " %g , " , beta);
        slen = strlen(buf);

        sprintf(buf + slen, " %g , " , nStoichElectrons);
        slen = strlen(buf);

        sprintf(buf + slen, " %g , %g , %g , %g , %g , %g , " , af, bf, Ea_kjgm, io_af, io_bf, io_Ea_kjgm);
        slen = strlen(buf);

        sprintf(buf + slen, " %g , " , damk);
        slen = strlen(buf);

        fprintf(fcsv,"%s\n", buf); 

        // New section to

        printf(" Reaction %d\n", (int) iR); 


        double kf = af * pow(TT,bf) * exp(-Ea / (GasConstant * TT));
        printf("  kf = %g \n", kf);

        double io_kf =  io_af * pow(TT,io_bf) * exp(-io_Ea / (GasConstant * TT)); 
        printf("  io_kf = %g \n", io_kf);

        double actMolal[30];
        brine->getActivities(actMolal);

        double deltaGSS[5];
        gKinetics->getDeltaSSGibbs(deltaGSS);

        if (iR == 0 && doAnode) {

           double deltaG0 = deltaGSS[0];
           double prodReactantCs = 1.0;
           double iok =  kf * Faraday * exp(beta *  deltaG0 / (GasConstant * TT) ) * prodReactantCs;

           printf(" these should be equal: %g %g \n", iok, io_kf);
 
           size_t iCupp = brine->speciesIndex("Cu++");
           double pc = 1.0 * (1.0 - beta);           
           double aM_Cupp = actMolal[iCupp];
           printf("     aM_Cupp = %g\n", aM_Cupp);
           double tmp = pow(aM_Cupp , pc);
           double ii_oc = io_kf * nStoichElectrons * tmp;

           printf("   These are two ways to calculate the exchange current density. They should be equal:\n");
           printf("     resulting ioc_alt = %g\n", ii_oc);
           printf("               io   = %g\n", io);

        }

    }

    fclose (fcsv);

    cleanup(thermo_vec, gKinetics);

    return 0;
}
//==================================================================================================================================
