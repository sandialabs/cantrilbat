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
PolarizationSurfRxnResults::PolarizationSurfRxnResults(int electrodeDomainNumber, int electrodeCellNumber,  Electrode* ee,
                                                       size_t surfIndex, size_t rxnIndex) :
    electrodeDomainNumber_(electrodeDomainNumber),
    electrodeCellNumber_(electrodeCellNumber),
    ee_(ee),
    iSurf_(surfIndex),
    iRxn_(rxnIndex),
    electronProd_(0.0),
    deltaTime_ (0.0)
{
    ocvSurfRxn = ee_->openCircuitVoltageRxn(iSurf_, iRxn_);
    ocvSurf = ocvSurfRxn;
    VoltageTotal = ocvSurfRxn;
}

//==================================================================================================================================
void PolarizationSurfRxnResults::addSubStep(struct PolarizationSurfRxnResults& sub)
{
    // We assume that the reaction #, the surface id, the cell # and the domain # are all the same
   
    // Develop a weighting of currents.
    double bef = fabs(electronProd_) / (fabs(electronProd_) + fabs(sub.electronProd_) + 1.0E-100);
    double aft = 1.0 - bef;

    // Add the currents. -> currents are deprecated
    icurrSurf_ += sub.icurrSurf_;
    electronProd_ += sub.electronProd_;
    deltaTime_ += sub.deltaTime_;

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
void PolarizationSurfRxnResults::subtractSubStep(struct PolarizationSurfRxnResults& sub)
{
    // We assume that the reaction #, the surface id, the cell # and the domain # are all the same
   

    // subtract the charges 
    icurrSurf_ -= sub.icurrSurf_;
    electronProd_ -= sub.electronProd_;
    deltaTime_ -= sub.deltaTime_;

    if (fabs(VoltageElectrode  - sub.VoltageElectrode) > 0.00001) {
        throw Electrode_Error("PolarizationSurfRxnResults::subtractSubStep()", "inconsistent voltages: %g %g\n",
                             VoltageElectrode, sub.VoltageElectrode);
    }
}
//==================================================================================================================================
void PolarizationSurfRxnResults::addOverPotentialPol(double overpotential, double nStoichElectrons,  int region, bool dischargeDir)
{
    // this has to be first
    double voltsS;
    bool anodeDischargeDir = true;
    if (region == 0) {
       if (!dischargeDir) anodeDischargeDir = false;
    } else if (region == 2) {
       if (dischargeDir) anodeDischargeDir = false;
    }
    if (anodeDischargeDir) {
        voltsS  = nStoichElectrons * overpotential;
    } else {
        voltsS  = - nStoichElectrons * overpotential;
    }





}
//==================================================================================================================================
// Add contribution for solid electronic conduction through the electrode's solid
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
    if (region == 0) {
        phi_anode_point_ = phiCurrentCollector;
    } else  if (region == 2) {
        phi_cathode_point_ =  phiCurrentCollector;
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
void PolarizationSurfRxnResults::addLyteCondPol_Sep(double phiLyteBoundary, double phiLyte_Spoint, int region, bool dischargeDir)
{
    double voltsS;
    bool anodeDischargeDir = true;
    if (region == 0) {
       if (!dischargeDir) anodeDischargeDir = false;
    } else if (region == 2) {
       if (dischargeDir) anodeDischargeDir = false;
    }
    if (anodeDischargeDir) {
        voltsS  = phiLyteBoundary - phiLyte_Spoint;
    } else {
        voltsS  = phiLyte_Spoint - phiLyteBoundary;
    }
    VoltPolPhenom ess(VOLT_LOSS_LYTE_SEP_PL, region, voltsS);
    bool found = false;
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == VOLT_LOSS_LYTE_SEP_PL) {
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
void PolarizationSurfRxnResults::addLyteConcPol(double* state_Lyte_Electrode, double* state_Lyte_SeparatorBdry,
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
    ReactingSurDomain* rsd = ee_->reactingSurface(iSurf_);

    // Fetch the number of stoichiometric electrons for the current reaction
    doublevalue nStoich = rsd->nStoichElectrons(iRxn_);

    // extract temperatures
    double T_Electrode =  state_Lyte_Electrode[0];
    double T_Separator =  state_Lyte_SeparatorBdry[0];

    size_t numKinSpecies = rsd->nKinSpecies();
    std::vector<double> netStoichVec(numKinSpecies);

    // Get the phase number of the electrolyte within the Electrode object
    size_t lytePN = ee_->solnPhaseIndex();
    size_t nspLyte = ee_->numSolnPhaseSpecies();

    // First calculate the activities at the electrode object

    /*  get the ThermoPhase for Lyte */
    thermo_t& tpLyte = ee_->thermo(lytePN);

    /*  Set the state to that at the electrode object */
    tpLyte.restoreState(3+nspLyte, state_Lyte_Electrode);

    std::vector<double> actElectrode(nspLyte);
    tpLyte.getActivities(actElectrode.data());

    tpLyte.restoreState(3+nspLyte, state_Lyte_SeparatorBdry);

    std::vector<double> actSeparator(nspLyte);
    tpLyte.getActivities(actSeparator.data());

    // Get the vector of stoichiometric coefficients

    rsd->getNetStoichCoeffVector(iRxn_, netStoichVec.data());

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
 
    // Adjust the total voltage from cathode to anode    
    if (dischargeDir) {
        VoltageTotal -= pLoss;
    } else {
        VoltageTotal += pLoss;
    }
}
//==================================================================================================================================
void PolarizationSurfRxnResults::addLyteConcPol_Sep(double* state_Lyte_SepBdry, double* state_Lyte_SepMid, 
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
    ReactingSurDomain* rsd = ee_->reactingSurface(iSurf_);

    // Fetch the number of stoichiometric electrons for the current reaction
    doublevalue nStoich = rsd->nStoichElectrons(iRxn_);

    // extract temperatures
    double T_SepBdry = state_Lyte_SepBdry[0];
    double T_SepMid  = state_Lyte_SepMid[0];

    size_t numKinSpecies = rsd->nKinSpecies();
    std::vector<double> netStoichVec(numKinSpecies);

    // Get the phase number of the electrolyte within the Electrode object
    size_t lytePN = ee_->solnPhaseIndex();
    size_t nspLyte = ee_->numSolnPhaseSpecies();

    // First calculate the activities at the electrode object

    /*  get the ThermoPhase for Lyte */
    thermo_t& tpLyte = ee_->thermo(lytePN);

    /*  Set the state to that at the electrode object */
    tpLyte.restoreState(3+nspLyte, state_Lyte_SepBdry);

    std::vector<double> actSepBdry(nspLyte);
    tpLyte.getActivities(actSepBdry.data());

    tpLyte.restoreState(3+nspLyte, state_Lyte_SepMid);
    std::vector<double> actSepMid(nspLyte);
    tpLyte.getActivities(actSepMid.data());

    // Get the vector of stoichiometric coefficients

    rsd->getNetStoichCoeffVector(iRxn_, netStoichVec.data());

    size_t lyte_KinP = rsd->solnPhaseIndex();
    size_t kinStart = rsd->kineticsSpeciesIndex(lyte_KinP, 0);

    double pLoss = 0.0;
    for (size_t k = 0; k < nspLyte; ++k) {
        double sc = netStoichVec[kinStart + k];
        if (sc != 0.0) {
            double deltaV = GasConstant * T_SepBdry * log(actSepBdry[k]) - GasConstant * T_SepMid * log(actSepMid[k]);
            deltaV *= signADD * sc /(Faraday * nStoich);
            pLoss += deltaV;
        }
    }

    // Keep track of the adjustments in the effective OCV.
    ocvSurfRxnAdj += pLoss;

    VoltPolPhenom ess(CONC_LOSS_LYTE_SEP_PL, region, pLoss);
    bool found = false;
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == CONC_LOSS_LYTE_SEP_PL) {
            found = true;
            vp = ess;
        }
    }
    if (! found) {
        voltsPol_list.push_back(ess);
    }
 
    // Adjust the total voltage from cathode to anode    
    if (dischargeDir) {
        VoltageTotal -= pLoss;
    } else {
        VoltageTotal += pLoss;
    }

}
//==================================================================================================================================
/*
 *  Add polarization due to the solid phase conduction
 */
void PolarizationSurfRxnResults::addSolidElectrodeConcPol(int region, bool dischargeDir)
{

    Electrode_Types_Enum eT = ee_->electrodeType();
    if (eT != SIMPLE_DIFF_ET) {
        return;
    }

    double signADD = 1.0;
    if (!dischargeDir) {
        signADD = -1.0;
    }
    double signAC = 1.0;
    if (region == 0) {
        signAC = -1.0;
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

    double ocv_mixAvg = ee_->openCircuitVoltage_MixtureAveraged(iSurf_);

    double ocv_diff = ee_->openCircuitVoltage(iSurf_);

    double deltaV = ocv_diff - ocv_mixAvg;
    double pLoss = deltaV * signAC * signADD;
    if (pLoss > 0.0 && dischargeDir) {
        throw Electrode_Error("PolarizationSurfRxnResults::addSolidElectrodeConcPol()",
                              "Negative polarization contribution, investigate %g %g %g", 
                              pLoss, ocv_diff, ocv_mixAvg);
    }
    if (pLoss < 0.0 && dischargeDir == false) {
        throw Electrode_Error("PolarizationSurfRxnResults::addSolidElectrodeConcPol()",
                              "Negative polarization contribution, investigate %g %g %g", 
                              pLoss, ocv_diff, ocv_mixAvg);
    }

    // Keep track of the adjustments in the effective OCV.
    ocvSurfRxnAdj += deltaV;

    VoltPolPhenom ess(SOLID_DIFF_CONC_LOSS_PL, region, pLoss);
    bool found = false;
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == SOLID_DIFF_CONC_LOSS_PL) {
            found = true;
            vp = ess;
            break;
        }
    }
    if (! found) {
        voltsPol_list.push_back(ess);
    }
    
    // adjust the storred OCV -> we have to adjust the storred OCV because it was previously
    // attributed to the OCV, but in fact it was diffusion resistance.
    found = false;
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == SURFACE_OCV_PL) {
            found = true;
            if (region == 0) {
               deltaV = -deltaV;
            }
            vp.voltageDrop += deltaV;
            break;
        }
    }
    if (!found) {
        throw Electrode_Error("PolarizationSurfRxnResults::addSolidElectrodeConcPol()", "logic error");
    }

    // Adjust the total voltage from cathode to anode    
    if (dischargeDir) {
        VoltageTotal -= pLoss;
    } else {
        VoltageTotal += pLoss;
    }
}
//==================================================================================================================================
bool PolarizationSurfRxnResults::checkConsistency(const double gvoltageTotal)
{
    bool res = true;
    if (fabs( gvoltageTotal - VoltageTotal) > 1.0E-5) {
        printf("PolarizationSurfRxnResults:: Error record has total voltage drop of %g instread of %g",
               VoltageTotal,  gvoltageTotal);
        res = false;
    }
    /*
     *  Sum up all of the records and see that it is equal to VoltageTotal
     */

    double vv = 0.0;
    for (VoltPolPhenom& vp : voltsPol_list) {
        vv += vp.voltageDrop;
    }

    if (fabs( vv  - VoltageTotal) > 1.0E-5) {
        printf("PolarizationSurfRxnResults:: Error record has total voltage drop of %g but contributions sum to %g",
               VoltageTotal,  vv);
        res = false;
    }

    return res; 
}
//==================================================================================================================================
std::string polString(enum Polarization_Loss_Enum plr)
{
   std::string s = "";
   switch (plr) {
   case UNKNOWN_PL:
       s = "unknown";
       break;
   case SURFACE_OCV_PL:
      s = "Surf_OCV"; 
      break;
    case SURFACE_OVERPOTENTIAL_PL:
      s = "surf_overpotential";
      break;
    case CONC_LOSS_LYTE_PL:
      s = "ConcLoss_Lyte";
      break;
    case VOLT_LOSS_LYTE_SEP_PL:
      s = "VoltLoss_Lyte_Sep";
      break;
   case  CONC_LOSS_LYTE_SEP_PL:
      s = "ConcLoss_Lyte_Sep";
      break;
    case VOLT_LOSS_LYTE_PL:
      s = "VoltLoss_Lyte";
      break;
    case CONC_LOSS_BL_LYTE_PL:
      s = "ConcLoss_LyteBL"; 
      break;
    case RESISTANCE_FILM_PL:
      s = "Resistance_Film";
      break;
   case  ELECTRICAL_CONDUCTION_LOSS_PL:
      s = "ElectCond_Electrode";
      break;
   case  ELECT_COND_COLLECTORS_LOSS_PL:
      s = "ElectCond_Collectors";
      break;
   case  SOLID_DIFF_CONC_LOSS_PL:
      s = "ConcLoss SolidDiff";
      break;
   case  CONC_LOSS_SOLID_FILM_PL:
      s = "ConcLoss SEI";
      break;
   case  VOLT_LOSS_SOLID_FILM_PL:
      s = "VoltLoss SEI";
      break;
    //! Other loss mechanisms not yet identified
   case   OTHER_PL:
      s = "other";
      break;
   default:
      throw ZuzaxError("Electrode::pString()", "Unknown poltype %d", plr);
      break;
   }

   return s;
}
//==================================================================================================================================
double totalElectronSource(const std::vector<struct PolarizationSurfRxnResults>& polarSrc_list)
{
#ifdef DEBUG_MODE
    double dt = 0.0;
#endif
    
    double totalSrc = 0.0;
    for (size_t i = 0; i < polarSrc_list.size(); ++i) {
        const struct PolarizationSurfRxnResults& ipol = polarSrc_list[i];
        totalSrc += ipol.electronProd_;
#ifdef DEBUG_MODE
        if (i == 0) {
            dt = ipol.deltaTime_;
        }  else {
           if (fabs( dt - ipol.deltaTime_) > std::max(0.001 * dt , 1.0E-12)) {
               printf("Caution delta times aren't the same: %g %g\n", dt, ipol.deltaTime_);
           }
        }
#endif
    }
    return totalSrc;
}
//==================================================================================================================================
} 
//----------------------------------------------------------------------------------------------------------------------------------
