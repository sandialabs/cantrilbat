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
    dischargeDir_ (true),
    electrodeDomainNumber_(electrodeDomainNumber),
    electrodeCellNumber_(electrodeCellNumber),
    ee_(ee),
    iSurf_(surfIndex),
    iRxn_(rxnIndex),
    electronProd_(0.0),
    deltaTime_ (0.0),
    ocvSurf(0.0),
    ocvSurfRxnAdj(0.0),
    VoltageElectrode_(0.0),
    phiMetal_(0.0),
    phi_lyteAtElectrode(0.0),
    phi_lyteSpoint(0.0),
    phi_anode_point_(0.0),
    phi_cathode_point_(0.0),
    VoltageTotal(0.0)
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
    electronProd_ += sub.electronProd_;
    deltaTime_ += sub.deltaTime_;

    if (ocvSurfRxn != sub.ocvSurfRxn) {
        ocvSurfRxn = sub.ocvSurfRxn;
    }

    bool same = true;
    if (VoltageElectrode_ != sub.VoltageElectrode_) {
        VoltageElectrode_ = sub.VoltageElectrode_;
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
    electronProd_ -= sub.electronProd_;
    deltaTime_ -= sub.deltaTime_;

    if (fabs(VoltageElectrode_  - sub.VoltageElectrode_) > 0.00001) {
        throw Electrode_Error("PolarizationSurfRxnResults::subtractSubStep()", "inconsistent voltages: %g %g\n",
                             VoltageElectrode_, sub.VoltageElectrode_);
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
    VoltPolPhenom ess(SURFACE_OVERPOTENTIAL_PL, region, voltsS);

}
//==================================================================================================================================
// Add contribution for solid electronic conduction through the electrode's solid
void PolarizationSurfRxnResults::addSolidPol(double phiCurrentCollector, int region, bool dischargeDir)
{
    double voltsS;
    bool found = false;
    double vOld = 0.0;

    if (region == 0) {
        voltsS = phiMetal_ - phiCurrentCollector;
    } else if (region == 2) {
        voltsS = phiCurrentCollector - phiMetal_;
    }

    VoltPolPhenom ess(ELECTRICAL_CONDUCTION_LOSS_PL, region, voltsS);
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == ELECTRICAL_CONDUCTION_LOSS_PL) {
            found = true;
            vOld = vp.voltageDrop;
            vp = ess;
            break;
        }
    }
    if (! found) {
        voltsPol_list.push_back(ess);
    }
 
    // Adjust the total voltage from cathode to anode, sin
    if (dischargeDir) {
        VoltageTotal += voltsS - vOld;
    } else {
        VoltageTotal += voltsS - vOld;
    }
    if (region == 0) {
        phi_anode_point_ = phiCurrentCollector;
    }
    if (region == 2) {
        phi_cathode_point_ = phiCurrentCollector;
    }
} 
//==================================================================================================================================
//  do this after the addSolidPol() calc.
void PolarizationSurfRxnResults::addSolidCCPol(double phiTerminal, double phiCurrentCollector, int region, bool dischargeDir)
{
    double voltsS;
    bool found = false;
    double vOld = 0.0;

    if (region == 0) {
        // usually negative
        voltsS  = phiCurrentCollector - phiTerminal;
    } else if (region == 2){
        // usually negative for discharge dir
        voltsS  = phiTerminal - phiCurrentCollector;
    } else {
        throw Electrode_Error("PolarizationSurfRxnResults::addSolidCCPol", "unknown region");
    }

    VoltPolPhenom ess(ELECT_COND_COLLECTORS_LOSS_PL, region, voltsS);
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == ELECT_COND_COLLECTORS_LOSS_PL) {
            found = true;
            vOld = vp.voltageDrop;
            vp = ess;
            break;
        }
    }
    if (! found) {
        voltsPol_list.push_back(ess);
    }
 
    // Adjust the total voltage from cathode to anode    
    if (dischargeDir) {
        VoltageTotal += voltsS - vOld;
    } else {
        VoltageTotal += voltsS - vOld;
    }
    if (region == 0) {
        phi_anode_point_ = phiTerminal;
    }
    if (region == 2) {
        phi_cathode_point_ = phiTerminal;
    }

}
//==================================================================================================================================
// Treats to boundary of electrode - separator region
void PolarizationSurfRxnResults::addLyteCondPol(double phiLyteElectrode, double phiLyteBoundary, int region, bool dischargeDir)
{
    double voltsS;
    double vOld = 0.0;
    if (phiLyteElectrode != phi_lyteAtElectrode) {
        throw Electrode_Error("PolarizationSurfRxnResults::addLyteCondPol()", "phiLyteElectrode !=  phi_lyteAtElectrode %g %g\n",
                             phiLyteElectrode,  phi_lyteAtElectrode);
    }
    if (region == 0) {
        voltsS  = phiLyteBoundary - phiLyteElectrode;
    } else if (region == 2) {
        voltsS  = phiLyteElectrode - phiLyteBoundary;
    }
    VoltPolPhenom ess(VOLT_LOSS_LYTE_PL, region, voltsS);
    bool found = false;
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == VOLT_LOSS_LYTE_PL) {
            found = true;
            vOld = vp.voltageDrop;
            vp = ess;
            break;
        }
    }
    if (! found) {
        voltsPol_list.push_back(ess);
    }

    // Adjust the total voltage from cathode to anode
    if (dischargeDir) {
        VoltageTotal += voltsS - vOld;
    } else {
        VoltageTotal += voltsS - vOld;
    }
    if (region == 0) {
        phi_cathode_point_ = phiLyteBoundary;
    }
    if (region == 2) {
        phi_anode_point_ = phiLyteBoundary;
    }

}
//==================================================================================================================================
void PolarizationSurfRxnResults::addLyteCondPol_Sep(double phiLyteBoundary, double phiLyte_Spoint, int region, bool dischargeDir)
{
    double voltsS;
    double vOld = 0.0;
    bool found = false;

    if (region == 0) {
        voltsS  = phiLyte_Spoint - phiLyteBoundary;
    } else if (region == 2) {
        voltsS  = phiLyteBoundary - phiLyte_Spoint;
    }

    VoltPolPhenom ess(VOLT_LOSS_LYTE_SEP_PL, region, voltsS);
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == VOLT_LOSS_LYTE_SEP_PL) {
            found = true;
            vOld = vp.voltageDrop;
            vp = ess;
            break;
        }
    }
    if (! found) {
        voltsPol_list.push_back(ess);
    }

    // Adjust the total voltage from cathode to anode
    if (dischargeDir) {
        VoltageTotal += voltsS - vOld;
    } else {
        VoltageTotal += voltsS - vOld;
    }
    if (region == 0) {
        phi_cathode_point_ = phiLyte_Spoint;
    }
    if (region == 2) {
        phi_anode_point_ = phiLyte_Spoint;
    }

}
//==================================================================================================================================
void PolarizationSurfRxnResults::addLyteConcPol(const double* state_Lyte_Electrode, const double* state_Lyte_SeparatorBdry,
                                                int region, bool dischargeDir)
{
    double signADD = 1.0;
    if (region == 2) {
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
            double deltaV = - GasConstant * T_Electrode * log(actElectrode[k]) + GasConstant * T_Separator * log(actSeparator[k]);
            deltaV *= signADD * sc /(Faraday * nStoich);
            pLoss += deltaV;
        }
    }
    if (dischargeDir) {
        if (pLoss > 0.0)  {
            throw Electrode_Error("PolarizationSurfRxnResults::addLyteConcPol()", 
                                  "Possible sign error on discargeDir ploss = %g\n", pLoss);
        }
    }

    VoltPolPhenom ess(CONC_LOSS_LYTE_PL, region, pLoss);
    bool found = false;
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == CONC_LOSS_LYTE_PL) {
            found = true;
            vp = ess;
            break;
        }
    }
    if (! found) {
        voltsPol_list.push_back(ess);
    }

    // Adjust the storred OCV -> we have to adjust the storred OCV because it was previously
    // attributed to the OCV, but in fact it was diffusion resistance.
    found = false;
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == SURFACE_OCV_PL) {
            found = true;
            vp.voltageDrop -= pLoss;
            break;
        }
    }
    if (!found) {
        throw Electrode_Error("PolarizationSurfRxnResults::addSolidElectrodeConcPol()", "logic error");
    }
    // Keep track of the adjustments in the effective OCV.
    ocvSurfRxnAdj -= pLoss;

}
//==================================================================================================================================
void PolarizationSurfRxnResults::addLyteConcPol_Sep(const double* state_Lyte_SepBdry, const double* state_Lyte_SepMid, 
                                                    int region, bool dischargeDir)
{
    double signADD = 1.0;
    if (region == 2) {
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
            double deltaV = - GasConstant * T_SepBdry * log(actSepBdry[k]) + GasConstant * T_SepMid * log(actSepMid[k]);
            deltaV *= signADD * sc /(Faraday * nStoich);
            pLoss += deltaV;
        }
    }

    VoltPolPhenom ess(CONC_LOSS_LYTE_SEP_PL, region, pLoss);
    bool found = false;
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == CONC_LOSS_LYTE_SEP_PL) {
            found = true;
            vp = ess;
            break;
        }
    }
    if (! found) {
        voltsPol_list.push_back(ess);
    }
    if (dischargeDir) {
        if (pLoss > 0.0)  {
            throw Electrode_Error("PolarizationSurfRxnResults::addLyteConcPol()", 
                                  "Possible sign error on discargeDir ploss = %g\n", pLoss);
        }
    }

    // Adjust the storred OCV -> we have to adjust the storred OCV because it was previously
    // attributed to the OCV, but in fact it was electrolyte diffusion resistance.
    found = false;
    for (VoltPolPhenom& vp : voltsPol_list) {
        if (vp.ipolType == SURFACE_OCV_PL) {
            found = true;
            vp.voltageDrop -= pLoss;
            break;
        }
    }
    if (!found) {
        throw Electrode_Error("PolarizationSurfRxnResults::addSolidElectrodeConcPol()", "logic error");
    }
    // Keep track of the adjustments in the effective OCV.
    ocvSurfRxnAdj += pLoss;

}
//==================================================================================================================================
/*
 *  Add polarization due to the solid phase conduction - HKM checked
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
    ocvSurfRxnAdj -= deltaV;

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
            vp.voltageDrop -= deltaV;
            break;
        }
    }
    if (!found) {
        throw Electrode_Error("PolarizationSurfRxnResults::addSolidElectrodeConcPol()", "logic error");
    }

    // Do not adjust the total voltage from cathode to anode, because we adjusted the OCV record

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
/*
 *  Bin according to surface and rxn number
 */
void agglomerate_init(std::vector<struct PolarizationSurfRxnResults>& polarSrc_agglom,
                      const std::vector<struct PolarizationSurfRxnResults>& polarSrc_list)
{
    bool found = false;
    for (size_t i = 0; i < polarSrc_list.size(); ++i) {
        const struct PolarizationSurfRxnResults& ipol = polarSrc_list[i];
        found = false; 
        for (size_t j = 0; j < polarSrc_agglom.size(); ++j) {
            struct PolarizationSurfRxnResults& jpol = polarSrc_agglom[j];
            if ((jpol.iSurf_ == ipol.iSurf_) && (jpol.iRxn_ == ipol.iRxn_)) {
                found = true;
                jpol.electronProd_ += ipol.electronProd_;
                break;
            }
        }
        if (!found) {
            polarSrc_agglom.push_back(ipol);
            struct PolarizationSurfRxnResults& jpol = polarSrc_agglom.back();
            jpol.electrodeCellNumber_ = -1;
            jpol.voltsPol_list.clear();
        }
    }
}
//==================================================================================================================================
void agglomerate_IV_add(struct PolarizationSurfRxnResults& apol,
                        const std::vector<struct PolarizationSurfRxnResults>& polarSrc_list)
{
    double ivprod;
    bool found;
    for (size_t i = 0; i < polarSrc_list.size(); ++i) {
        const struct PolarizationSurfRxnResults& ipol = polarSrc_list[i];
        if ((apol.iSurf_ == ipol.iSurf_) && (apol.iRxn_ == ipol.iRxn_)) {
             // Add the electron production found
             apol.electronProd_ += ipol.electronProd_;

             // For each record  VoltPolPhenom record in ipol
             for (const VoltPolPhenom& vp : ipol.voltsPol_list) {
                 ivprod = vp.voltageDrop * ipol.electronProd_;
                 // search for a match
                 found = false;
                 for (VoltPolPhenom& avp : apol.voltsPol_list) {
                     if ((vp.ipolType == avp.ipolType) && (vp.regionID == avp.regionID)) {
                         found = true;
                         avp.voltageDrop += ivprod;              
                         break;
                     }
                 }
                 if (!found) {
                     apol.voltsPol_list.push_back(vp);
                     VoltPolPhenom& avp = apol.voltsPol_list.back();
                     avp.voltageDrop = ivprod;
                 }
             }
             break;
        }
    }
}
//==================================================================================================================================
} 
//----------------------------------------------------------------------------------------------------------------------------------
