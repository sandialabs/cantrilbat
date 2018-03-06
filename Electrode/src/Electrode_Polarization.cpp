/**
 *  @file Electrode_Polarization.cpp
 */

/*
 * Copywrite 2018 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "Electrode_Polarization.h"
#include "Electrode_Exception.h"

#include "Electrode.h"

#include <cmath>

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
PolarizationSurfRxnResults::PolarizationSurfRxnResults(int electrodeDomainNumber, int electrodeCellNumber, 
                                                       size_t surfIndex, size_t rxnIndex) :
    electrodeDomainNumber_(electrodeDomainNumber),
    electrodeCellNumber_(electrodeCellNumber),
    ee(nullptr),
    isurf_(surfIndex),
    iRxnIndex_(rxnIndex)
{
}
//==================================================================================================================================
void PolarizationSurfRxnResults::addSubStep(struct PolarizationSurfRxnResults& sub)
{
    // We assume that the reaction #, the surface id, the cell # and the domain # are all the same
   
    // Develop a weighting of currents.
    double bef = fabs(icurrSurf) / (fabs(icurrSurf) + fabs(sub.icurrSurf) + 1.0E-100);
    double aft = 1.0 - bef;

    // Add the currents.
    icurrSurf += sub.icurrSurf;

    if (ocvSurfRxn != sub.ocvSurfRxn) {
        ocvSurfRxn = sub.ocvSurfRxn;
    }

    bool same = true;
    if (VoltageElectrode  != sub.VoltageElectrode) {
        VoltageElectrode = sub.VoltageElectrode;
        same = false;
    }

    // Right now it seems the best idea is to replace the polarization list with last polarization list ?!?
    /*
     *  We may want to try an averaging procedure based on currents
     */
    if (same) {
        for (size_t i = 0 ; i < voltsPol_list.size(); ++i) {
            struct VoltPolPhenom& vpp = voltsPol_list[i];
            const struct VoltPolPhenom& vpps = sub.voltsPol_list[i];
            vpp.voltageDrop = bef * vpp.voltageDrop + aft * vpps.voltageDrop;
            if (vpp.ipolType != vpps. ipolType) {
                throw Electrode_Error("PolarizationSurfRxnResults::addSubStep()", "assumptions wrong");
            }
        }
    } else {
        voltsPol_list = sub.voltsPol_list;
    }

}
//==================================================================================================================================
// add contribution for solid electronic conduction through the electrode's solid
void PolarizationSurfRxnResults::addSolidPol(double phiCurrentCollector, int region, bool dischargeDir)
{
    double voltsS;
    bool anodeDischargeDir = true;
    if (region == 0) {
       if (!dischargeDir) anodeDischargeDir = false;
    } else if (region == 2) {
       if (dischargeDir) anodeDischargeDir = false;
    }
    if (anodeDischargeDir) {
        voltsS  = phiCurrentCollector - phiMetal;
    } else {
        voltsS  = phiMetal - phiCurrentCollector;
    }
    VoltPolPhenom ess(ELECTRICAL_CONDUCTION_LOSS_PL, region, voltsS);
    bool found = false;
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == ELECTRICAL_CONDUCTION_LOSS_PL) {
            found = true;
            vp = ess;
        }
    }
    if (! found) {
        voltsPol_list.push_back(ess);
    }
 
    // Adjust the total voltage from cathode to anode    
    if (dischargeDir) {
        VoltageTotal -= voltsS;
    } else {
        VoltageTotal += voltsS;
    }
}
//==================================================================================================================================
void PolarizationSurfRxnResults::addSolidCCPol(double phiTerminal, double phiCurrentCollector, int region, bool dischargeDir)
{
    double voltsS;
    bool anodeDischargeDir = true;
    if (region == 0) {
       if (!dischargeDir) anodeDischargeDir = false;
    } else if (region == 2) {
       if (dischargeDir) anodeDischargeDir = false;
    }
    if (anodeDischargeDir) {
        voltsS  = phiTerminal - phiCurrentCollector;
    } else {
        voltsS  = phiCurrentCollector - phiTerminal;
    }
    VoltPolPhenom ess(ELECT_COND_COLLECTORS_LOSS_PL, region, voltsS);
    bool found = false;
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == ELECT_COND_COLLECTORS_LOSS_PL) {
            found = true;
            vp = ess;
        }
    }
    if (! found) {
        voltsPol_list.push_back(ess);
    }
 
    // Adjust the total voltage from cathode to anode    
    if (dischargeDir) {
        VoltageTotal -= voltsS;
    } else {
        VoltageTotal += voltsS;
    }

}
//==================================================================================================================================
void PolarizationSurfRxnResults::addLyteCondPol(double phiLyteElectrode, double phiLyteBoundary, int region, bool dischargeDir)
{
    double voltsS;
    bool anodeDischargeDir = true;
    if (region == 0) {
       if (!dischargeDir) anodeDischargeDir = false;
    } else if (region == 2) {
       if (dischargeDir) anodeDischargeDir = false;
    }
    if (anodeDischargeDir) {
        voltsS  = phiLyteElectrode - phiLyteBoundary;
    } else {
        voltsS  = phiLyteBoundary - phiLyteElectrode;
    }
    VoltPolPhenom ess(VOLT_LOSS_LYTE_PL, region, voltsS);
    bool found = false;
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == VOLT_LOSS_LYTE_PL) {
            found = true;
            vp = ess;
        }
    }
    if (! found) {
        voltsPol_list.push_back(ess);
    }

    // Adjust the total voltage from cathode to anode
    if (dischargeDir) {
        VoltageTotal -= voltsS;
    } else {
        VoltageTotal += voltsS;
    }

}
//==================================================================================================================================
void PolarizationSurfRxnResults::addLyteConcPol(double* state_Lyte_Electrode, double* state_Lyte_SeparatorMid,
                                                int region, bool dischargeDir)
{
    bool anodeDischargeDir = true;
    if (region == 0) {
       if (!dischargeDir) anodeDischargeDir = false;
    } else if (region == 2) {
       if (dischargeDir) anodeDischargeDir = false;
    }
    double signADD = 1.0;
    if (!anodeDischargeDir) {
        signADD = -1.0;
    }



    // Fetch the ReactingSurDomain object for the current
    ReactingSurDomain* rsd = ee->reactingSurface(isurf_);

    // Fetch the number of stoichiometric electrons for the current reaction
    doublevalue nStoich = rsd->nStoichElectrons(iRxnIndex_);

    // extract temperatures
    double T_Electrode =  state_Lyte_Electrode[0];
    double T_Separator =  state_Lyte_SeparatorMid[0];

    size_t numKinSpecies = rsd->nKinSpecies();
    std::vector<double> netStoichVec(numKinSpecies);

    // Get the phase number of the electrolyte within the Electrode object
    size_t lytePN = ee->solnPhaseIndex();
    size_t nspLyte = ee->numSolnPhaseSpecies();

    // First calculate the activities at the electrode object

    /*  get the ThermoPhase for Lyte */
    thermo_t& tpLyte = ee->thermo(lytePN);

    /*  Set the state to that at the electrode object */
    tpLyte.restoreState(3+nspLyte, state_Lyte_Electrode);

    std::vector<double> actElectrode(nspLyte);
    tpLyte.getActivities(actElectrode.data());

    tpLyte.restoreState(3+nspLyte, state_Lyte_SeparatorMid);

    std::vector<double> actSeparator(nspLyte);
    tpLyte.getActivities(actSeparator.data());

    // Get the vector of stoichiometric coefficients

    rsd->getNetStoichCoeffVector(iRxnIndex_, netStoichVec.data());

    size_t lyte_KinP = rsd->solnPhaseIndex();
    size_t kinStart = rsd->kineticsSpeciesIndex(lyte_KinP, 0);

    double pLoss = 0.0;
    for (size_t k = 0; k < nspLyte; ++k) {
        double sc = netStoichVec[kinStart + k];
        if (sc != 0.0) {
            double deltaV = GasConstant * T_Electrode * log(actElectrode[k]) - GasConstant * T_Separator * log(actSeparator[k]);
            deltaV *= signADD * sc /(Faraday * nStoich);
            pLoss += deltaV;
        }
    }

    // Keep track of the adjustments in the effective OCV.
    ocvSurfRxnAdj += pLoss;

    VoltPolPhenom ess(CONC_LOSS_LYTE_PL, region, pLoss);
    bool found = false;
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == CONC_LOSS_LYTE_PL) {
            found = true;
            vp = ess;
        }
    }
    if (! found) {
        voltsPol_list.push_back(ess);
    }



}
//==================================================================================================================================
void PolarizationSurfRxnResults::addSolidElectrodeConcPol(double* mf_OuterSurf_Electrode, double* mf_Average_Electrode,
                                                          int region, bool dischargeDir)
{

    Electrode_Types_Enum eT = ee->electrodeType();
    if (eT != SIMPLE_DIFF_ET) {
        return;
    }

    double signADD = 1.0;
    if (!dischargeDir) {
        signADD = -1.0;
    }

    // Fetch the ReactingSurDomain object for the current
    //ReactingSurDomain* rsd = ee->reactingSurface(isurf_);

    // Fetch the number of stoichiometric electrons for the current reaction
    //doublevalue nStoich = rsd->nStoichElectrons(iRxnIndex_);
    //size_t numKinSpecies = rsd->nKinSpecies();

    // Get the phase number of the electrolyte within the Electrode object
    //size_t lytePN = ee->solnPhaseIndex();
    //size_t nspLyte = ee->numSolnPhaseSpecies();
    //size_t lyte_KinP = rsd->solnPhaseIndex();
    //size_t kinStart = rsd->kineticsSpeciesIndex(lyte_KinP, 0);

    double ocv_mixAvg = ee->openCircuitVoltage_MixtureAveraged(isurf_);

    double ocv_diff = ee->openCircuitVoltage(isurf_);

    double pLoss = 0.0;
    double deltaV = ocv_diff - ocv_mixAvg;
    deltaV *= signADD;
    pLoss += deltaV;

    // Keep track of the adjustments in the effective OCV.
    ocvSurfRxnAdj += deltaV;

    VoltPolPhenom ess(SOLID_DIFF_CONC_LOSS_PL, region, pLoss);
    bool found = false;
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == SOLID_DIFF_CONC_LOSS_PL) {
            found = true;
            vp = ess;
        }
    }
    if (! found) {
        voltsPol_list.push_back(ess);
    }

}
//==================================================================================================================================
} 
//----------------------------------------------------------------------------------------------------------------------------------
