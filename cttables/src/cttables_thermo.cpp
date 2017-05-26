/**
 *  @file cttables_thermo.cpp
 *
 *
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "cttables.h"

#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/IdealSolidSolnPhase.h"

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

#include <cstdio>

using std::endl;
using std::cout;

//===================================================================================================================================
// Gather the entropy of the elements of a species at 298 K. This is useful for going back and forth from the
// gibbs free energy of formation and the absolute gibbs free energy in NIST format.
//
double entropyElem298(ZZCantera::ThermoPhase* g_ptr, size_t k)
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
    return stotal;
}
//===================================================================================================================================
void printThermoCoeffSpecies(ThermoPhase* g_ptr, int k)
{

    int j, nint, nzones = 1, i;
    std::string sName = g_ptr->speciesName(k);
    double* cptr;
    int index;
    size_t kindex;
    int type;
    PDSS* pdss_ptr = nullptr;
    // Get the species thermo manager
    SpeciesThermo& sp = g_ptr->speciesThermo();
    VPStandardStateTP* vpss = nullptr;
    PDSS_enumType ptype;
    double DH0_tr_pr;
    double S0_tr_pr;
    double Mu0_tr_pr;
    double DG0_tr_pr;
    double dg_consistent, as, bs, cs;

    // Get the report Type
    int rt = sp.reportType(k);
    double c[200];

    // Get the coefficients of the parameterization
    int rt2;
    double minTemp, maxTemp, refPressure;
    sp.reportParams(k, rt2, c, minTemp, maxTemp, refPressure);
    if (rt != rt2) {
        printf("report types don't match ERROR\n");
        exit(-1);
    }
    double stotal = entropyElem298(g_ptr, k);

    // switch on the report type, printing out a formatted representation of the polynomials
    switch (rt) {
    case NASA2:
        dnt(1); printf("NASA Polynomial format: 2 zones\n");
        dnt(2); printf("Low  temperature polynomials: %g < T < %g: \n", minTemp, c[0]);
        dnt(2); printf("%12.5g %12.5g %12.5g %12.5g\n", c[3], c[4],  c[5],  c[6]);
        dnt(2); printf("%12.5g %12.5g %12.5g\n", c[7], c[1], c[2]);
        dnt(2); printf("High temperature polynomials: %g < T < %g: \n", c[0], maxTemp);
        dnt(2); printf("%16.10g %16.10g %16.10g %16.10g\n", c[10], c[11], c[12],  c[13]);
        dnt(2); printf("%16.10g %16.10g %16.10g\n", c[14], c[8], c[9]);
        break;

    case SHOMATE2:
        dnt(1); printf("SHOMATE Polynomial format: 2 zones\n");
        dnt(2); printf("Low  temperature polynomials: %g < T < %g: \n", minTemp, c[0]);
        dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n", c[1],  c[2], c[3], c[4]);
        dnt(2); printf("%17.11g %17.11g %17.11g\n", c[5],  c[6], c[7]);
        dnt(2); printf("High temperature polynomials: %g < T < %g: \n", c[0], maxTemp);
        dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n", c[8],  c[9], c[10], c[11]);
        dnt(2); printf("%17.11g %17.11g %17.11g\n", c[12], c[13], c[14]);
        break;

    case NASA1:
        dnt(1); printf("NASA Polynomial format: 1 zones\n");
        dnt(2); printf("temperature polynomials: %g < T < %g: \n", minTemp, maxTemp);
        dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n", c[1], c[2],  c[3],  c[4]);
        dnt(2); printf("%17.11g %17.11g %17.11g\n", c[5], c[6], c[7]);
        break;

    case SHOMATE1:
        dnt(1); printf("SHOMATE Polynomial format: 1 zones\n");
        dnt(2); printf("temperature polynomials: %g < T < %g: \n", minTemp, maxTemp);
        dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n", c[0],  c[1], c[2], c[3]);
        dnt(2); printf("%17.11g %17.11g %17.11g\n", c[4],  c[5], c[6]);
        break;

    case CONSTANT_CP:
    case SIMPLE:
        dnt(1); printf("CONSTANT_CP format:\n");
        dnt(2); printf(" Valid Range: %g < T < %g: \n", minTemp, maxTemp);
        dnt(2); printf(" at T = %5g K, H  = %17.11g J kmol-1\n", c[0], c[1]);
        dnt(2); printf("                 S  = %17.11g J kmol-1 K-1\n", c[2]);
        dnt(2); printf("                 Cp = %17.11g J kmol-1 K-1\n", c[3]);
        break;

    case MU0_INTERP:
        dnt(1); printf("MU0_POLY format:\n");
        nint = (int) c[0];
        dnt(2); printf(" Valid Range: %g < T < %g: \n", minTemp, maxTemp);
        dnt(2); printf(" T        mu_0\n");
        j = 2;
        for (i =0 ; i < nint; i++) {
            dnt(2); printf(" %g   %17.11g\n", c[j], c[j+1]);
            j += 2;
        }
        break;

    case NASA9:
        nzones = 1;
        dnt(1); printf("NASA9 Polynomial format: %d zones\n", nzones);
        cptr = c;
        dnt(2); printf("%2d:  Temperature polynomials: %g < T < %g: \n", 0, minTemp, maxTemp);
        dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n", cptr[3], cptr[4],  cptr[5],  cptr[6]);
        dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n", cptr[7], cptr[8], cptr[9], cptr[10]);
        dnt(2); printf("%17.11g\n", cptr[11]);
        break;

    case NASA9MULTITEMP:
        nzones = (int) c[0];
        dnt(1);
        printf("NASA9 Polynomial format: %d zones\n", nzones);
        index = 1;
        for (i = 0; i < nzones; i++) {
            minTemp = c[index];
            maxTemp = c[index+1];
            cptr = c+index+2;
            dnt(2); printf("%2d:  Temperature polynomials: %g < T < %g: \n", i, minTemp, maxTemp);
            dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n", cptr[0], cptr[1],  cptr[2],  cptr[3]);
            dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n", cptr[4], cptr[5], cptr[6], cptr[7]);
            dnt(2); printf("%17.11g\n", cptr[8]);
            index += 11;
        }
        break;

    case MINEQ3:
        nzones = 1;
        DH0_tr_pr = c[7];
        S0_tr_pr = c[8];
        Mu0_tr_pr = c[9];
        DG0_tr_pr = c[10];
        as = c[0] / 4.184;
        bs = c[1] / (4.184 * 1000.);
        cs = c[4] * 1.0E6 / 4.184;

        if (stotal == ENTROPY298_UNKNOWN) {
            dg_consistent = ENTROPY298_UNKNOWN;
        } else {
            dg_consistent =  Mu0_tr_pr + 298.15 * (stotal);
        }
        dnt(1); printf("MinEQ3 format:  (a varient of Shomate1 format) \n");
        dnt(2); printf("temperature polynomials (Shomate Form): %g < T < %g: \n", minTemp, maxTemp);
        dnt(2); printf("%17.11g %17.11g %17.11g %17.11g\n", c[0],  c[1], c[2], c[3]);
        dnt(2); printf("%17.11g %17.11g %17.11g\n",  c[4],  c[5], c[6]);
        dnt(2); printf("  Delta G0_Tr_Pr = %16.9E cal/gmol\n", DG0_tr_pr / (4.184 * 1.0E3));
        dnt(2); printf("                 = %16.6g kJ /gmol\n", DG0_tr_pr/1.0E6);
        dnt(2); printf("  Delta H0_Tr_Pr = %16.9E cal/gmol\n", DH0_tr_pr / (4.184 * 1.0E3));
        dnt(2); printf("        S0_Tr_Pr = %16.9g cal/gmol/K\n", S0_tr_pr / (4.184 * 1.0E3));
        dnt(2); printf("                 = %16.6g  J /gmol/K\n", S0_tr_pr / 1.0E3);
        dnt(2); printf("       mu0_Tr_Pr = %16.6g kJ /gmol\n",   Mu0_tr_pr / 1.0E6);
        if (dg_consistent == ENTROPY298_UNKNOWN) {
            dnt(2); printf(" Delta G0_consis = [UNAVAILABE BECAUSE ENTROPY298 NOT INPUT]\n");
        } else {
            dnt(2); printf(" Delta G0_consis = %16.6g kJ /gmol\n", dg_consistent / 1.0E6);
        }
        dnt(2); printf("               a = %16.9g cal/gkmol/K\n", as);
        dnt(2); printf("               b = %16.9g cal/kmol/K2\n", bs);
        dnt(2); printf("               c = %16.9g cal-K/gmol\n", cs);
        break;

    case PDSS_TYPE:
        dnt(1); printf("Presure Dependent Standard State form\n");
        vpss = dynamic_cast<VPStandardStateTP*>(g_ptr);
        pdss_ptr = vpss->providePDSS(k);
        ptype = pdss_ptr->reportPDSSType();

        double minTemp, maxTemp, refPressure;
        pdss_ptr->reportParams(kindex, type, &c[0], minTemp, maxTemp, refPressure);

        dnt(2); printf("PDSS Species index    = %d\n", (int) kindex);
        dnt(2); printf("PDSS min Temperature  = %g\n", minTemp);
        dnt(2); printf("PDSS max Temperature  = %g\n", maxTemp);
        dnt(2); printf("PDSS reference pres   = %g\n", refPressure);

        switch (ptype) {
        case cPDSS_IDEALGAS:
            dnt(2); printf("Ideal Gas PDSS thermo type\n");
            break;
        case  cPDSS_CONSTVOL:
            dnt(2); printf("ConstVol PDSS thermo type\n");
            break;
        case cPDSS_MOLAL_CONSTVOL:
            dnt(2); printf("Molal ConstVol PDSS thermo type\n");
            break;
        case cPDSS_WATER:
            dnt(2); printf("Water PDSS thermo type\n");
            break;
        case cPDSS_MOLAL_HKFT:
            dnt(2); printf("PDSS_MOLAL_HKFT PDSS thermo type\n");
            dnt(3); printf("deltaG_formation_tr_pr = %g\n", c[0]);
            dnt(3); printf("deltaH_formation_tr_pr = %g\n", c[1]);
            dnt(3); printf("Mu0_tr_pr              = %g\n", c[2]);
            dnt(3); printf("Entrop_tr_pr           = %g\n", c[3]);
            dnt(3); printf("a1                     = %g\n", c[4]);
            dnt(3); printf("a2                     = %g\n", c[5]);
            dnt(3); printf("a3                     = %g\n", c[6]);
            dnt(3); printf("a4                     = %g\n", c[7]);
            dnt(3); printf("c1                     = %g\n", c[8]);
            dnt(3); printf("c2                     = %g\n", c[9]);
            dnt(3); printf("omega_pr_tr            = %g\n", c[10]);

            break;
        case cPDSS_UNDEF:
            dnt(2); printf("Undefined PDSS thermo type\n");
            break;
        default:
            printf("unknown species PDSS thermo type %d\n", (int) ptype);
        }

        break;
    default:
        printf("unknown species reference thermo type %d\n", rt);
    }
}
//==================================================================================================================================
/*
 *  This routine adds the molar volume information onto a species table.
 *  This is not part of the reference state.  This is part of the description of the phase.
 */
void printVolSpecies(ThermoPhase* const g_ptr, size_t k)
{
    std::string sName = g_ptr->speciesName(k);
    int eos = g_ptr->eosType();
    //int dim = g_ptr->nDim();
    size_t nS = g_ptr->nSpecies();
    std::vector<double> ssVol(nS, 0.0);
    g_ptr->getStandardVolumes(DATA_PTR(ssVol));
    std::vector<double> pmVol(nS, 0.0);
    g_ptr->getPartialMolarVolumes(DATA_PTR(pmVol));
    double currTemp = g_ptr->temperature();
    double currPres = g_ptr->pressure();

    dnt(1); printf("Partial    Molar Volume at (T,P) of (%g,%g) = %g m%1d/kmol\n", currTemp, currPres, pmVol[k], 3);
    dnt(1); printf("StandState Molar Volume at (T,P) of (%g,%g) = %g m%1d/kmol\n", currTemp, currPres, ssVol[k], 3);

    /*
     * Loop over the known types of eos's. Print out what's available.
     */
    if (eos == cIdealSolidSolnPhase0  || eos == cIdealSolidSolnPhase1 || eos == cIdealSolidSolnPhase2) {

        IdealSolidSolnPhase* issp = static_cast<IdealSolidSolnPhase*>(g_ptr);
        double mV = issp->speciesMolarVolume(k);
        dnt(1); cout << "species Molar Volume = " << mV << " m**3 / kmol\n";

        switch (eos) {
        case cIdealSolidSolnPhase0:
            dnt(1); cout << "standard Concentration value =  unity\n";
            break;
        case cIdealSolidSolnPhase1:
            dnt(1); cout << "standard Concentration value =  1/(molar volume) units of kmol/m3\n";
            break;
        case cIdealSolidSolnPhase2:
            dnt(1); cout << "standard Concentration value =  1/(solvent volume) units of kmol/m3\n";
            break;
        }

        dnt(1); cout << "Pres dependence -> (U, S, V, and Cp are independent)\n";
        dnt(1); cout << "                -> (H and G dependent): h(T) = u(T) + pv" << endl;

    } else if (eos == cStoichSubstance) {

        //StoichSubstance* issp = static_cast<StoichSubstance*>(g_ptr);
        double dens = g_ptr->density();
        double mw = g_ptr->molecularWeight(0);
        double mV = mw / (dens);
        dnt(1); cout << "StoichSubstance constant density = " << dens << " kg m-3 \n";
        dnt(1); cout << "species Molar Volume = " << mV << " m**3 / kmol <- constant\n";
        dnt(1); cout << "Press dependence -> (U, S, V, and Cp are independent)\n";
        dnt(1); cout << "                 -> (H and G dependent): h(T) = u(T) + pv\n";
        dnt(1); cout << "standard Concentration value =  unity" << endl;

    } else if (eos == cIdealGas) {

        dnt(1); cout << "Volume Properties are from ideal gas law\n";
        dnt(1); cout << "Standard concentration value =  P / (R T) with units of kmol/m3" << endl;

    } else if (eos == cDebyeHuckel0  || eos == cDebyeHuckel1  || eos == cDebyeHuckel2) {

        double mV = pmVol[k];
        dnt(1); cout << "species partial Molar Volume = " << mV << " m**3 / kmol\n";
        switch (eos) {
        case cDebyeHuckel0:
            dnt(1); cout << "standard Concentration value =  unity\n";
            break;
        case cDebyeHuckel1:
            dnt(1); cout << "standard Concentration value =  1/(molar volume) units of kmol/m3\n";
            break;
        case cDebyeHuckel2:
            dnt(1); cout << "standard Concentration value =  1/(solvent volume) units of kmol/m3\n";
            break;
        }
        dnt(1); cout << "Pres dependence -> (U, S, V, and Cp are independent)\n";
        dnt(1); cout << "                -> (H and G dependent): h(T) = u(T) + pv" << endl;

    } else if (eos == cHMWSoln0  || eos == cHMWSoln1  || eos == cHMWSoln2) {

        double mV = pmVol[k];
        dnt(1); cout << "species partial Molar Volume = " << mV << " m**3 / kmol\n";
        switch (eos) {
        case cHMWSoln0:
            dnt(1); cout << "standard Concentration value =  unity" << endl;
            break;
        case  cHMWSoln1 :
            dnt(1); cout << "standard Concentration value =  1/(molar volume) units of kmol/m3" << endl;
            break;
        case cHMWSoln2:
            dnt(1); cout << "standard Concentration value =  1/(solvent volume) units of kmol/m3" << endl;
            break;
        }

    } else if (eos == cSurf) {

        dnt(1); cout << "Volume Properties are from fixed surface site equation\n";
        SurfPhase* isurf = static_cast<SurfPhase*>(g_ptr);
        double sited = isurf->siteDensity();
        sited  *= Avogadro / 1.0E4;
        dnt(1); cout << "Surface Site Density = " << sited << " # per cm2 " << endl;

    } else {

        dnt(1); cout << "Molar Volume formulation for eos type " << eos << " is not known to cttables" << endl;

    }
}
//==================================================================================================================================

