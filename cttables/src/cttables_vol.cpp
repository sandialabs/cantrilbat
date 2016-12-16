/**
 *  @file cttables_vol.cpp
 *
 */

#include "cttables.h"

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/mix_defs.h"
#include "cantera/thermo/IdealSolidSolnPhase.h"

#include <iostream>

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

using  std::cout;
using  std::endl;

//==================================================================================================================================
/*
 *  This routine adds onto a species table the molar volume information
 *  supplied for the species. This is not part of the reference state.  This is part of the description of the phase.
 */
void printVolSpecies(ThermoPhase* g_ptr, size_t k)
{
    std::string sName = g_ptr->speciesName(k);
    int eos = g_ptr->eosType();
    int dim = g_ptr->nDim();
    size_t nS = g_ptr->nSpecies();
    std::vector<double> ssVol(nS, 0.0);
    g_ptr->getStandardVolumes(DATA_PTR(ssVol));
    std::vector<double> pmVol(nS, 0.0);
    g_ptr->getPartialMolarVolumes(DATA_PTR(pmVol));
    double currTemp = g_ptr->temperature();
    double currPres = g_ptr->pressure();
    dnt(1);
    printf("Partial    Molar Volume at (T,P) of (%g,%g) = %g m%1d/kmol\n", currTemp, currPres, pmVol[k], dim);
    dnt(1);
    printf("StandState Molar Volume at (T,P) of (%g,%g) = %g m%1d/kmol\n", currTemp, currPres, ssVol[k], dim);

    /*
     * Loop over the known types of eos's. Print out what's available.
     */
    if (eos == cIdealSolidSolnPhase0  || eos == cIdealSolidSolnPhase1 || eos == cIdealSolidSolnPhase2) {

        IdealSolidSolnPhase* issp = static_cast<IdealSolidSolnPhase*>(g_ptr);
        double mV = issp->speciesMolarVolume(k);
        dnt(1);
        cout << "species Molar Volume = " << mV << " m**3 / kmol\n";
        switch (eos) {
        case cIdealSolidSolnPhase0:
            dnt(1);
            cout << "standard Concentration value =  unity\n";
            break;
        case cIdealSolidSolnPhase1:
            dnt(1);
            cout << "standard Concentration value =  1/(molar volume) units of kmol/m3\n";
            break;
        case cIdealSolidSolnPhase2:
            dnt(1);
            cout << "standard Concentration value =  1/(solvent volume) units of kmol/m3\n";
            break;
        }
        dnt(1);
        cout << "Pres dependence -> (U, S, V, and Cp are independent)\n";
        dnt(1);
        cout << "                -> (H and G dependent): h(T) = u(T) + pv" << endl;
    } else if (eos == cStoichSubstance) {
        //StoichSubstance* issp = static_cast<StoichSubstance*>(g_ptr);
        double dens = g_ptr->density();
        double mw = g_ptr->molecularWeight(0);
        double mV = mw / (dens);
        dnt(1);
        cout << "StoichSubstance constant density = " << dens << " kg m-3 \n";
        dnt(1);
        cout << "species Molar Volume = " << mV << " m**3 / kmol <- constant\n";
        dnt(1);
        cout << "Press dependence -> (U, S, V, and Cp are independent)\n";
        dnt(1);
        cout << "                 -> (H and G dependent): h(T) = u(T) + pv\n";
        dnt(1);
        cout << "standard Concentration value =  unity" << endl;
    } else if (eos == cIdealGas) {
        dnt(1);
        cout << "Volume Properties are from ideal gas law\n";
        dnt(1);
        cout << "Standard concentration value =  P / (R T) with units of kmol/m3" << endl;

    } else if (eos == cDebyeHuckel0  ||
               eos == cDebyeHuckel1  ||
               eos == cDebyeHuckel2) {

        double mV = pmVol[k];
        dnt(1);
        cout << "species partial Molar Volume = " << mV << " m**3 / kmol\n";
        switch (eos) {
        case cDebyeHuckel0:
            dnt(1);
            cout << "standard Concentration value =  unity\n";
            break;
        case cDebyeHuckel1:
            dnt(1);
            cout << "standard Concentration value =  1/(molar volume) units of kmol/m3\n";
            break;
        case cDebyeHuckel2:
            dnt(1);
            cout << "standard Concentration value =  1/(solvent volume) units of kmol/m3\n";
            break;
        }
        dnt(1);
        cout << "Pres dependence -> (U, S, V, and Cp are independent)\n";
        dnt(1);
        cout << "                -> (H and G dependent): h(T) = u(T) + pv" << endl;
    } else if (eos == cHMWSoln0  ||
               eos == cHMWSoln1  ||
               eos == cHMWSoln2) {
        double mV = pmVol[k];
        dnt(1);
        cout << "species partial Molar Volume = " << mV << " m**3 / kmol\n";
        switch (eos) {
        case cHMWSoln0:
            dnt(1);
            cout << "standard Concentration value =  unity" << endl;
            break;
        case  cHMWSoln1 :
            dnt(1);
            cout << "standard Concentration value =  1/(molar volume) units of kmol/m3" << endl;
            break;
        case cHMWSoln2:
            dnt(1);
            cout << "standard Concentration value =  1/(solvent volume) units of kmol/m3" << endl;
            break;
        }
    } else if (eos == cSurf) {
        dnt(1);
        cout << "Volume Properties are from fixed surface site equation\n";
        SurfPhase* isurf = static_cast<SurfPhase*>(g_ptr);
        double sited = isurf->siteDensity();
        sited  *= Avogadro / 1.0E4;
        dnt(1);
        cout << "Surface Site Density = " << sited << " # per cm2 " << endl;
    } else {
        dnt(1);
        cout << "Molar Volume formulation for eos type " << eos << " is not known to cttables" << endl;
    }
}
//==================================================================================================================================

