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

#include <cstdio>

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
    isurf(surfIndex),
    iRxnIndex(rxnIndex)
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
void PolarizationSurfRxnResults::addSolidPol(double phiCurrentCollector, int region)
{
    double voltsS;
    if (region == 0) {
        voltsS  = phiCurrentCollector - phiMetal;
    } else if (region == 2) {
        voltsS  = phiMetal - phiCurrentCollector;
    }
    VoltPolPhenom ess(ELECTRICAL_CONDUCTION_LOSS_PL, region, voltsS) ;
    voltsPol_list.push_back(ess);
 
    // Adjust the total voltage from cathode to anode    
    VoltageTotal -= voltsS;
}
//==================================================================================================================================
} 
//----------------------------------------------------------------------------------------------------------------------------------
