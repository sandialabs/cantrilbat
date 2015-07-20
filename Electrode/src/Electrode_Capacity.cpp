/*
 * $Id: Electrode_Capacity.cpp 590 2013-05-06 17:47:49Z hkmoffa $
 */

#include "Electrode.h"

#include "cantera/base/ctexceptions.h"

#include <cstdio>

namespace Cantera
{


// -----------------------------------------------------------------------------------------------------------------
// --------------------------  CAPACITY SETUP AND OUTPUT -----------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------


// ------------------------------------------------------------------------------------------------
// --------------------------- CAPACITY CALCULATION OUTPUT  ---------------------------------------
// ------------------------------------------------------------------------------------------------
//====================================================================================================================
// Returns the total capacity of the electrode in Amp seconds
/*
 *  Returns the capacity of the electrode in Amps seconds.
 *  This is the same as the number of coulombs that can be delivered at any voltage.
 *  Note, this number differs from the capacity of electrodes that is usually quoted for
 *  a battery. That number depends on the rate of discharge and also depends on the
 *  specification of a cutoff voltage. Here, we dispense with both of these specifications.
 *  So, it should be considered a theoretical capacity at zero current and minimal cutoff voltage.
 *  It will also include all plateaus that are defined by the electrode object.
 *
 *  @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
 *
 * @return returns the theoretical capacity of the electrode in Amp sec = coulombs.
 */
double Electrode::capacity(int platNum) const
{
    if (electrodeChemistryModelType_ == 2) {
        setCapacityCoeff_FeS2();
    }
    double capZeroDoD = 0.0;
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        if (iph == (size_t) solnPhase_ || iph == (size_t) metalPhase_) {
            continue;
        }
        int kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        int nspPhase = tp.nSpecies();
        for (int k = 0; k < nspPhase; k++) {
            double ll = spMoles_final_[kStart + k];
            capZeroDoD += ll * capacityZeroDoDSpeciesCoeff_[kStart + k];
        }
    }
    double tmp = capZeroDoD * Faraday;
    return tmp;
}

//====================================================================================================================
// Initial capacity of the elctrode in Amp seconds
/*
 *  This is the initial capacity of the electrode before any degradation occurs.
 *
 *  @param platNum Plateau number. Default is -1 which treats all plateaus as a single entity.
 */
double Electrode::capacityInitial(int platNum) const
{
    return capacityInitialZeroDod_;
}
//====================================================================================================================
// Amount of charge that the electrode that has available to be discharged
/*
 *   We report the number in terms of Amp seconds = coulombs
 */
double Electrode::capacityLeft(int platNum, double voltsMax, double voltsMin) const
{
    if (electrodeChemistryModelType_ == 2) {
        setCapacityCoeff_FeS2();
    }
    double capLeft = 0.0;
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        if (iph == (size_t) solnPhase_ || iph == (size_t) metalPhase_) {
            continue;
        }
        int kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        int nspPhase = tp.nSpecies();
        for (int k = 0; k < nspPhase; k++) {
            double ll = spMoles_final_[kStart + k];
            capLeft += ll * capacityLeftSpeciesCoeff_[kStart + k];
        }
    }
    return capLeft * Faraday;
}
//====================================================================================================================
// Report the current depth of discharge in amps sec = coulumbs
/*
 * Report the current depth of discharge. This is roughly equal to the total
 * number of electrons that has been theoretically discharged from a fully charged state compared
 * to the number of electrons that can theoretically be discharged.
 *
 * Usually this is reported as a function of the discharge rate and there is a
 * cutoff voltage at which the electron counting is turned off. Neither of these
 * concepts is employed here.
 *
 *  @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
 *
 * @return  returns the depth of discharge in amp seconds = coulombs
 */
double Electrode::depthOfDischarge(int platNum) const
{
    double capLeft = capacityLeft(platNum);
    double capZeroDod = capacity(platNum);
    double dod = capZeroDod - capLeft;
    return (dod);
}
//====================================================================================================================
//  Check on the accounting of the capacity, done at t_final
/*
 *  This check should be as extensive as possible. It should be more extensive in the child routines. The
 *  parent routines does a cursory accounting.
 *  
 *   Virtual from Electrode
 *
 *  @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
 *
 *  @return  true if everything checks out. False otherwise
 */
bool Electrode::checkCapacityBalances_final(int platNum) const
{
#ifdef DEBUG_MASSLOSS
    static int firstTime = 1;
#endif
    // The following is a default treatment 
    // It's not very 
    double capLeft = capacityLeft(platNum);

    double capDischarged = capacityDischarged(platNum);

    double cap = capacity(platNum);
    double capS = depthOfDischargeStarting(platNum);

    double capOrig = capacityInitial(platNum);
    double rel = capOrig - cap;
    double denom = std::max(capOrig, rel);
    denom = std::max(denom, 1.0E-50);
    if (abs(rel) > 1.0E-12 * denom) {
	printf("Electrode::checkCapacityBalances_final ERROR: init cap % 19.12E  not same as curr cap % 19.12E (%d)\n",
	       capOrig, cap, platNum);
#ifdef DEBUG_MODE
	throw CanteraError(" Electrode::checkCapacityBalances_final ERROR", "Capacity loss mechanism unaccounted for");
#endif
	return true;
    }

    rel = capOrig - capLeft - capDischarged - capS;
    denom = std::max(capOrig, capDischarged);
    denom = std::max(denom, 1.0E-50);
#ifdef DEBUG_MASSLOSS
    string fff = "capCheck_" + int2str(electrodeCellNumber_) + ".txt";
    FILE *fp = fopen(fff.c_str(), "a");
    if (firstTime) {
	fprintf(fp,"                  \"cap\"               \"capLeft\"             \" capDischarged\"            \"spMoles4\"   "
		"\"capLeft + capDisc\"    \"RelativeOff\"\n");
	firstTime = 0;
    }
    fprintf(fp,"(BeforeB)  % 21.15E % 21.15E  % 21.15E  % 21.15E % 21.15E % 21.15E\n", cap/Faraday, capLeft/Faraday, 
	    capDischarged/Faraday, capS/Faraday,
	    capLeft/Faraday + capDischarged/Faraday,  rel/Faraday);
    fclose(fp);
#endif
    if (abs(rel) > 1.0E-12 * denom) {
        if (printLvl_ > 7) {
	  printf("Electrode::checkCapacityBalances_final  cap not accounted for by % 19.12E coulombs (%d)\n", rel, platNum);
	  printf("\t\t\t  capOrig       = % 21.14E coulombs  \n", capOrig);
	  printf("\t\t\t -capAccounted  = % 21.14E\n", capLeft + capDischarged + capS);
  	  printf("\t\t\t- capleft       = % 21.14E\n", capLeft);
	  printf("\t\t\t- capDischarged = % 21.14E\n", capDischarged);
	  printf("\t\t\t- DoD_starting  = % 21.14E\n", capS);

	  printf("\t\t\t= Missing       = % 21.14E\n", rel  );
	  printf("\t\t\t                          =  % 21.15E moles\n", rel/Faraday);
        }
#ifdef DEBUG_MODE
	//throw CanteraError(" Electrode::checkCapacityBalances_final ERROR", "Capacity unaccounted for");
#endif
	return true;
    }
    return false;
}
//====================================================================================================================
// Experimental routine to enforce a balance on the electrode capacity given previous balance information
/* 
 *  This routine equalizes the capacity 
 */
void Electrode::fixCapacityBalances_final()
{
    double capLeft = capacityLeft();
    double capDischarged = capacityDischarged();
    double cap = capacity();
    double capS = depthOfDischargeStarting();
    double capOrig = capacityInitial();
    double capLost = capOrig - cap;
    // rel is the missing capacity that is numerically unaccounted for
    double rel = capOrig - capLeft - capDischarged - capS - capLost;
    // relMole is the moles that are missing.
    double relMole = rel / Faraday;
    int iKadd = -1;
    int iKsub = -1;
    if (relMole != 0.0) {

	double coeffA = 0.0;
	double coeffB = 0.0;
	for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
	    if (iph == (size_t) solnPhase_ || iph == (size_t)  metalPhase_) {
		continue;
	    }
	    int kStart = m_PhaseSpeciesStartIndex[iph];
	    ThermoPhase& tp = thermo(iph);
	    int nspPhase = tp.nSpecies();
	    for (int k = 0; k < nspPhase; k++) {
		if (capacityLeftSpeciesCoeff_[kStart + k] > 0.0 && spMoles_final_[kStart + k] > 0.0) {
		    iKadd = k + kStart;
		    coeffA = capacityLeftSpeciesCoeff_[kStart + k];
		}
		if (capacityLeftSpeciesCoeff_[kStart + k] == 0.0 && spMoles_final_[kStart + k] > 0.0) {
		    if (capacityZeroDoDSpeciesCoeff_[kStart + k] > 0.0) {
			iKsub = k + kStart;
			coeffB = capacityZeroDoDSpeciesCoeff_[kStart + k];
		    }
		}
	    }
	}

	if (iKadd >= 0) {
	    spMoles_final_[iKadd] -= relMole / coeffA;
	}
	if (iKsub >= 0) {
	    spMoles_final_[iKsub] += relMole / coeffB;
	}
    }
}
//====================================================================================================================
// Report the current depth of discharge in percent
/*
 * Report the current depth of discharge. This is roughly equal to the total
 * number of electrons that has been theoretically discharged from a fully charged state compared
 * to the number of electrons that can theoretically be discharged.
 *
 * Usually this is reported as a function of the discharge rate and there is a
 * cutoff voltage at which the electron counting is turned off. Neither of these
 * concepts is employed here.
 *
 *  @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
 *
 * @return  returns the depth of discharge in percent
 */
double Electrode::depthOfDischargeFraction(int platNum) const
{
    double capLeft = capacityLeft(platNum);
    double capZeroDod = capacity(platNum);
    double dod = capZeroDod - capLeft;
    return (dod/capZeroDod);
}
//====================================================================================================================
// Amount of charge that the electrode has discharged up to this point
/*
 *   We report the number in terms of Amp seconds = coulombs
 *
 *   @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
 */
double Electrode::capacityDischarged(int platNum) const
{
    //  We use the variable electronKmolDischargedToDate_ to keep track of what has gone before.
    double tmp = electronKmolDischargedToDate_;
    if (pendingIntegratedStep_) {
        tmp += spMoleIntegratedSourceTerm_[kElectron_];
    }
    double sgn = 1.0;
    if (electrodeCapacityType_ == CAPACITY_CATHODE_ECT) {
       sgn = -1.0;
    }
    return sgn * tmp * Cantera::Faraday;
}
//======================================================================================================================
double Electrode::depthOfDischargeStarting(int platNum) const
{
    return depthOfDischargeStarting_;
}
//======================================================================================================================
double Electrode::calcRelativeExtentRxn_final() const
{
    return depthOfDischargeFraction(-1);
}
//======================================================================================================================
// Reset the counters that keep track of the amount of discharge to date
void Electrode::resetCapacityDischargedToDate() 
{
    if (pendingIntegratedStep_) {
       throw CanteraError("Electrode::resetCapacityDischargedToDate() ERROR",
                          "called during a pending integration step");
    }
    electronKmolDischargedToDate_ = 0.0;
    depthOfDischargeStarting_ = depthOfDischarge();
}
//======================================================================================================================
Electrode_Capacity_Type_Enum Electrode::capacityType() const
{
    return electrodeCapacityType_;
}
//======================================================================================================================
void  Electrode::setCapacityType(Electrode_Capacity_Type_Enum electrodeCapacityType)
{
    electrodeCapacityType_ = electrodeCapacityType;
}
//====================================================================================================================
// Set parameters that tell the object how to calculate the capacity of the electrode
/*
 * @param sName          Name of the species that contains the capacity
 * @param coeffLeft      Coefficient describing how many electrons are left
 * @param coeffZeroDoD   Coefficient describing how many electrons are originally there
 */
void Electrode::setCapacityCalcParams(std::string sName, double coeffLeft, double coeffZeroDoD)
{
    bool found = false;
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        if (iph == (size_t) solnPhase_ || iph == (size_t) metalPhase_) {
            continue;
        }
        int kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        int nspPhase = tp.nSpecies();
        for (int k = 0; k < nspPhase; k++) {
            int iGlobSpeciesIndex = kStart + k;
            std::string sss = speciesName(iGlobSpeciesIndex);
            if (sName == sss) {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = coeffLeft;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = coeffZeroDoD;
                found = true;
                break;
            }
        }
    }
    if (!found) {
        throw CanteraError("Electrode::setCapacityCalcParams", "species not found");
    }
}
//====================================================================================================================
//   Set the capacity coefficients from the input file
/*
 *  Here we set :
 *
 *  These are set from the input file.
 *
 *      capacityLeftSpeciesCoeff_[iGlobSpeciesIndex
 *	    capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex]
 *
 *  It is a fatal error to call this member function if the coefficients have
 *  not been set from the input file.
 *
 *  @param ei   Input Key file object
 */
void Electrode::setCapacityCoeffFromInput(const ELECTRODE_KEY_INPUT* const ei)
{
    double** CapLeftPhase = ei->m_BG->CapLeftCoeffPhases;
    double** CapZeroDoDPhase = ei->m_BG->CapZeroDoDCoeffPhases;
    double* capLeft = 0;
    double* capZero = 0;
    bool foundNonZeroLeft = false;
    bool foundNonZeroDoD = false;
    if (!CapLeftPhase) {
	throw CanteraError("Electrode::setCapacityCoeffFromInput()", "Capacity Left coefficients not set in the input file");
    }
    if (!CapZeroDoDPhase) {
	throw CanteraError("Electrode::setCapacityCoeffFromInput()", "Capacity Zero DoD coefficients not set in the input file");
    }
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        if (iph == (size_t) metalPhase_) {
            continue;
        }
        int kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        int nspPhase = tp.nSpecies();
	capLeft = CapLeftPhase[iph];
	capZero = CapZeroDoDPhase[iph];
        for (int k = 0; k < nspPhase; k++) {
            int iGlobSpeciesIndex = kStart + k;
            std::string sss = speciesName(iGlobSpeciesIndex);
	    capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] =  capLeft[k];
	    capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = capZero[k];
	    if (capLeft[k] != 0.0) {
		foundNonZeroLeft = true;
		if (capZero[k] == 0.0) {
		    throw Electrode_Error("Electrode::setCapacityCoeffFromInput(", "Capacity Zero is zero but capacity Left is nonzero");
		}
	    }
            
	    if (capZero[k] != 0.0) {
		foundNonZeroDoD = true;
	    }
        }
    }
    if (!foundNonZeroLeft) {
        throw CanteraError("Electrode::setCapacityCoeffFromInput()", "Didn't find any nonzero capacity Left coefficients");
    }   
    if (!foundNonZeroDoD) {
        throw CanteraError("Electrode::setCapacityCoeffFromInput()", "Didn't find any nonzero Dod capacity coefficients");
    }  
}
//====================================================================================================================
void Electrode::setCapacityCoeff_LiSi() const
{
    double DoDElectronTmp;
    double DoDMoleTmp;
    double tmp;

    DoDElectronTmp = 11./3.;
    DoDMoleTmp = 4./3.;

    DoDElectronTmp += DoDMoleTmp * 13./7.;
    DoDMoleTmp *= 3./7.;

    DoDElectronTmp += DoDMoleTmp * 12.;

    double  capacityLeftSpeciesCoeffLi13Si4 =  DoDElectronTmp;

    for (int iph = 0; iph < NumVolPhases_; iph++) {
        if (iph == solnPhase_ || iph == metalPhase_) {
            continue;
        }

        int kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        int nspPhase = tp.nSpecies();
        std::string pname = PhaseNames_[iph];
        bool found = false;
        for (int k = 0; k < nspPhase; k++) {
            int iGlobSpeciesIndex = kStart + k;
            std::string sss = speciesName(iGlobSpeciesIndex);
            found = false;

            if (sss == "Li13Si4(S)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] =  capacityLeftSpeciesCoeffLi13Si4;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = capacityLeftSpeciesCoeffLi13Si4;

                capacityLeftSpeciesCoeffPlat_(0, iGlobSpeciesIndex) = 11./3.;
                capacityZeroDoDSpeciesCoeffPlat_(0, iGlobSpeciesIndex) = 11./3.;

                tmp = 4./3. * 13. / 7.;
                capacityLeftSpeciesCoeffPlat_(1,iGlobSpeciesIndex) = tmp;
                capacityZeroDoDSpeciesCoeffPlat_(1,iGlobSpeciesIndex) = tmp;

                tmp = 4./3. * 3. / 7. * 12.;
                capacityLeftSpeciesCoeffPlat_(2, iGlobSpeciesIndex) = tmp;
                capacityZeroDoDSpeciesCoeffPlat_(2, iGlobSpeciesIndex) = tmp;


                found = true;
            }


            if (sss == "Li7Si3(S)") {
                DoDElectronTmp = 13./7.;
                DoDMoleTmp = 3./7.;
                DoDElectronTmp += DoDMoleTmp * 12.;
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = DoDElectronTmp;
                DoDMoleTmp = 3./4.;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = DoDMoleTmp * capacityLeftSpeciesCoeffLi13Si4;

                capacityLeftSpeciesCoeffPlat_(0, iGlobSpeciesIndex) = 0.0;
                capacityZeroDoDSpeciesCoeffPlat_(0, iGlobSpeciesIndex) = 3./4. * 11. / 3.;

                capacityLeftSpeciesCoeffPlat_(1,iGlobSpeciesIndex) = 13./7.;
                capacityZeroDoDSpeciesCoeffPlat_(1,iGlobSpeciesIndex) = 13./7.;

                tmp = 3. / 7. * 12.;
                capacityLeftSpeciesCoeffPlat_(2,iGlobSpeciesIndex) = tmp;
                capacityZeroDoDSpeciesCoeffPlat_(2,iGlobSpeciesIndex) = tmp;

                found = true;
            }

            if (sss == "Li12Si7(S)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = 12.0;
                DoDMoleTmp = 7./3.;
                DoDMoleTmp *= 3./4.;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] =  DoDMoleTmp * capacityLeftSpeciesCoeffLi13Si4;

                capacityLeftSpeciesCoeffPlat_(0, iGlobSpeciesIndex) = 0.0;
                capacityZeroDoDSpeciesCoeffPlat_(0, iGlobSpeciesIndex) = 3./4. * 7./3. * 11. / 3.;

                capacityLeftSpeciesCoeffPlat_(1, iGlobSpeciesIndex) = 0.0;
                capacityZeroDoDSpeciesCoeffPlat_(1, iGlobSpeciesIndex) = 7./3. * 13. / 7.;

                capacityLeftSpeciesCoeffPlat_(2, iGlobSpeciesIndex) = 12.;
                capacityZeroDoDSpeciesCoeffPlat_(2, iGlobSpeciesIndex) = 12.;


                found = true;
            }

            if (sss == "Si(S)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = 0.0;
                DoDMoleTmp = 1./7.;
                DoDMoleTmp *= 7./3.;
                DoDMoleTmp *= 3./4.;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] =  DoDMoleTmp * capacityLeftSpeciesCoeffLi13Si4;


                capacityLeftSpeciesCoeffPlat_(0, iGlobSpeciesIndex) = 0.0;
                capacityZeroDoDSpeciesCoeffPlat_(0, iGlobSpeciesIndex) = 1./7. * 3./4. * 7./3. * 11. / 3.;

                capacityLeftSpeciesCoeffPlat_(1, iGlobSpeciesIndex) = 0.0;
                capacityZeroDoDSpeciesCoeffPlat_(1, iGlobSpeciesIndex) = 1./ 7. * 7./3. * 13. / 7.;

                capacityLeftSpeciesCoeffPlat_(2, iGlobSpeciesIndex) = 0.0;
                capacityZeroDoDSpeciesCoeffPlat_(2, iGlobSpeciesIndex) = 1./7. * 12.;

                found = true;
            }

            if (!found) {
                throw CanteraError(":setCapacityCoeff_LiSi()", "unknown species: " + sss);
            }

        }
    }
}
//====================================================================================================================
//
void Electrode::setCapacityCoeff_LiSi_Li() const
{
    double DoDElectronTmp;
    double DoDMoleTmp;

    DoDElectronTmp = 11./3.;
    DoDMoleTmp = 4./3.;


    double  capacityLeftSpeciesCoeffLi13Si4 =  DoDElectronTmp;

    for (int iph = 0; iph < NumVolPhases_; iph++) {
        if (iph == solnPhase_ || iph == metalPhase_) {
            continue;
        }

        int kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        int nspPhase = tp.nSpecies();
        std::string pname = PhaseNames_[iph];
        bool found = false;
        for (int k = 0; k < nspPhase; k++) {
            int iGlobSpeciesIndex = kStart + k;
            std::string sss = speciesName(iGlobSpeciesIndex);
            found = false;

            if (sss == "Li13Si4(S)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] =  capacityLeftSpeciesCoeffLi13Si4;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = capacityLeftSpeciesCoeffLi13Si4;
                found = true;
            }

            if (sss == "Li(i)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = 1.0;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = 0.0;
                found = true;
            }

            if (sss == "V(i)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = 0.0;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = 0.0;
                found = true;
            }

            if (sss == "Li7Si3(S)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = 0.0;
                DoDMoleTmp = 3./4.;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] =  DoDMoleTmp * capacityLeftSpeciesCoeffLi13Si4;
                found = true;
            }

            if (!found) {
                throw CanteraError(":setCapacityCoeff_LiSi_Li()", "unknown species: " + sss);
            }

        }
    }
}
//====================================================================================================================
// Set the Capacity coefficients for the FeS2 cathode system
/*
 *  @todo Move or get rid of
 */
void Electrode::setCapacityCoeff_FeS2() const
{
    double DoDElectronTmp;
    double DoDMoleTmp;

    DoDElectronTmp = 3./2.;
    DoDMoleTmp = 1./2.;

    DoDElectronTmp += DoDMoleTmp * 1.0;
    DoDMoleTmp *= 2.0;

    DoDElectronTmp += DoDMoleTmp * 2.0;

    double  capacityLeftSpeciesCoeffFeS2 =  DoDElectronTmp;

    for (int iph = 0; iph < NumVolPhases_; iph++) {
        if (iph == solnPhase_ || iph == metalPhase_) {
            continue;
        }

        int kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        int nspPhase = tp.nSpecies();
        std::string pname = PhaseNames_[iph];
        bool found = false;
        for (int k = 0; k < nspPhase; k++) {
            int iGlobSpeciesIndex = kStart + k;
            std::string sss = speciesName(iGlobSpeciesIndex);
            found = false;

            if (sss == "FeS2(S)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] =  capacityLeftSpeciesCoeffFeS2;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = capacityLeftSpeciesCoeffFeS2;
                found = true;
            }

            if (sss == "Li3Fe2S4(S)") {
                DoDElectronTmp = 1.0;
                DoDMoleTmp = 2.0;
                DoDElectronTmp += DoDMoleTmp * 2.0;
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = DoDElectronTmp;
                DoDMoleTmp = 2.0;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = DoDMoleTmp * capacityLeftSpeciesCoeffFeS2;
                found = true;
            }


            int i_FeS_B = globalSpeciesIndex("FeS(B)");
            int i_VaS_B = globalSpeciesIndex("VaS(B)");
            int i_Li73_A = globalSpeciesIndex("Li[7/3]Fe[2/3]S2(A)");

            /*
             *   Method for converting Li[7/3]Fe[2/3]S2(A) and FeS phase into DoD.
             *
             *   First the tie-line approach indicates that the coefficients are actually functions of the
             *   relative moles of the Li[7/3]Fe[2/3]S2(A) FeS_B, and VaS_B.
             *
             *   We will use the three reactions below to convert Li[7/3]Fe[2/3]S2(A) FeS_B, and VaS_B.
             *   into Li and Li2FeS2(A) and Fe(S). The extents of reactions, sigma, are used to solve
             *   for this system
             *
             *   6  Li[7/3]Fe[2/3]S2(A) + 4 FeS(B) + 2 Li+ + 2 e- = 8 Li2FeS2(A)       sigma_1
             *
             *   0.5 FeS(B) + 0.5 VaS(B) +  Li+ +  e- = 0.5 Li2FeS2(A)                 sigma_2
             *
             *   FeS(B)  [=]   Fe(S)  +  VaS(B)                                        sigma_3
             *
             *   sigma_1 = m_Li73 / 6
             *
             *   Left over FeS from the first reaction:
             *
             *       m1_FeS = m_FeS - 4 sigma_1
             *
             *       sigma_3 doesn't actually matter. Sigma2 may be negative or positive. It doesn't matter.
             *
             *    sigma2 = m1_FeS + m_VaS
             *
             *     Amount of electrons + Li+ consumed to get to Li2FeS2(A)
             *
             *       #electrons = 2 sigma_1 + sigma_2 = m_Li73/3 + m1_FeS + m_VaS = m_Li73/3 + m_FeS - 2/3 m_Li73 + m_VaS =
             *                  = m_FeS + m_VaS - 1/3 m_Li73
             *
             *       #Li2FeS2(A) = 8 sigma_1 + 0.5 sigma_2
             *                   = 4/3  m_Li73 + 0.5 m1_FeS +  0.5 m_VaS
             *                   = 4/3  m_Li73 + 0.5 m_FeS - 2 sigma_1 +  0.5 m_VaS
             *                   = 4/3  m_Li73 + 0.5 m_FeS - 1/3 m_Li73 +  0.5 m_VaS
             *                   =      m_Li73 + 0.5 m_FeS +  0.5 m_VaS
             *
             */


            if (sss == "Li[7/3]Fe[2/3]S2(A)") {
                DoDElectronTmp = -1./3.;
                DoDMoleTmp = 1.0;
                DoDElectronTmp += DoDMoleTmp * 2.0;
                capacityLeftSpeciesCoeff_[i_Li73_A] = DoDElectronTmp;

                DoDElectronTmp = 1.0;
                DoDMoleTmp = 0.5;
                DoDElectronTmp += DoDMoleTmp * 2.0;
                capacityLeftSpeciesCoeff_[i_FeS_B] = DoDElectronTmp;

                DoDElectronTmp = 1.0;
                DoDMoleTmp = 0.5;
                DoDElectronTmp += DoDMoleTmp * 2.0;
                capacityLeftSpeciesCoeff_[i_VaS_B] = DoDElectronTmp;

                DoDMoleTmp = 1.0;
                DoDMoleTmp *= 1./2.;
                DoDMoleTmp *= 2.0;
                capacityZeroDoDSpeciesCoeff_[i_Li73_A] =  DoDMoleTmp * capacityLeftSpeciesCoeffFeS2;

                DoDMoleTmp = 0.5;
                DoDMoleTmp *= 1./2.;
                DoDMoleTmp *= 2.0;
                capacityZeroDoDSpeciesCoeff_[i_FeS_B] =  DoDMoleTmp * capacityLeftSpeciesCoeffFeS2;

                DoDMoleTmp = 0.5;
                DoDMoleTmp *= 1./2.;
                DoDMoleTmp *= 2.0;
                capacityZeroDoDSpeciesCoeff_[i_VaS_B] =  DoDMoleTmp * capacityLeftSpeciesCoeffFeS2;
                found = true;
            }

            if (sss == "FeS(B)" || sss == "VaS(B)") {
                found = true;
            }

            if (sss == "Li2FeS2(A)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = 2.0;
                DoDMoleTmp = 1./2.;
                DoDMoleTmp *= 2.0;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] =  DoDMoleTmp * capacityLeftSpeciesCoeffFeS2;
                found = true;
            }

            if (sss == "Li2S(S)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = 0.0;
                DoDMoleTmp = 1./2.;
                DoDMoleTmp *= 1./2.;
                DoDMoleTmp *= 2.;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] =  DoDMoleTmp * capacityLeftSpeciesCoeffFeS2;
                found = true;
            }

            if (sss == "Fe(S)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = 0.0;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = 0.0;
                found = true;
            }

            if (!found) {
                throw CanteraError(":setCapacityCoeff_FeS2()", "unknown species: " + sss);
            }

        }
    }
}
//====================================================================================================================
// Set the Capacity coefficients for the modified FeS2 cathode system
/*
 *  @todo Move or get rid of
 */
void Electrode::setCapacityCoeff_FeS2_Combo() const
{
    double DoDElectronTmp;
    double DoDMoleTmp;

    DoDElectronTmp = 3./2.;
    DoDMoleTmp = 1./2.;

    DoDElectronTmp += DoDMoleTmp * 1.0;
    DoDMoleTmp *= 2.0;

    DoDElectronTmp += DoDMoleTmp * 2.0;

    double  capacityLeftSpeciesCoeffFeS2 =  DoDElectronTmp;

    for (int iph = 0; iph < NumVolPhases_; iph++) {
        if (iph == solnPhase_ || iph == metalPhase_) {
            continue;
        }

        int kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        int nspPhase = tp.nSpecies();
        std::string pname = PhaseNames_[iph];
        bool found = false;
        for (int k = 0; k < nspPhase; k++) {
            int iGlobSpeciesIndex = kStart + k;
            std::string sss = speciesName(iGlobSpeciesIndex);
            found = false;

            if (sss == "FeS2(S)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] =  capacityLeftSpeciesCoeffFeS2;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = capacityLeftSpeciesCoeffFeS2;
                found = true;
            }

            if (sss == "Li3Fe2S4(S)") {
                DoDElectronTmp = 1.0;
                DoDMoleTmp = 2.0;
                DoDElectronTmp += DoDMoleTmp * 2.0;
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = DoDElectronTmp;
                DoDMoleTmp = 2.0;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = DoDMoleTmp * capacityLeftSpeciesCoeffFeS2;
                found = true;
            }




            if (sss == "LitFe1S2(S)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = 2.38;
                DoDMoleTmp = 1./2.;
                DoDMoleTmp *= 2.0;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] =  DoDMoleTmp * capacityLeftSpeciesCoeffFeS2;
                found = true;
            }

            if (sss == "Li2Fe1S2(S)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = 2.0;
                DoDMoleTmp = 1./2.;
                DoDMoleTmp *= 2.0;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] =  DoDMoleTmp * capacityLeftSpeciesCoeffFeS2;
                found = true;
            }

            if (sss == "Li2S(S)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = 0.0;
                DoDMoleTmp = 1./2.;
                DoDMoleTmp *= 1./2.;
                DoDMoleTmp *= 2.;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] =  DoDMoleTmp * capacityLeftSpeciesCoeffFeS2;
                found = true;
            }

            if (sss == "Fe(S)") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = 0.0;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = 0.0;
                found = true;
            }

            if (!found) {
                throw CanteraError(":setCapacityCoeff_FeS2_Combo()", "unknown species: " + sss);
            }

        }
    }
}
//====================================================================================================================
//  Capacity Parameters for the MCMB anode
/*
 *  Here we supply parameters for the arrays
 *          capacityLeftSpeciesCoeff_[]
 *          capacityZeroDoDSpeciesCoeff_[]
 *  This is an intercalated electrode, named MCMB_Interstitials_anode, with the following species
 *         Li_C6-bulk
 *          V_C6-bulk
 *
 */
void Electrode::setCapacityCoeff_MCMB() const
{
    double  capacityLeftSpeciesCoeff=  1.0;

    for (int iph = 0; iph < NumVolPhases_; iph++) {
        if (iph == solnPhase_ || iph == metalPhase_) {
            continue;
        }

        int kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        int nspPhase = tp.nSpecies();
        std::string pname = PhaseNames_[iph];
        bool found = false;
        for (int k = 0; k < nspPhase; k++) {
            int iGlobSpeciesIndex = kStart + k;
            std::string sss = speciesName(iGlobSpeciesIndex);
            found = false;

            if (sss == "Li_C6-bulk") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex]    = capacityLeftSpeciesCoeff;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = capacityLeftSpeciesCoeff;
                found = true;
            }

            if (sss == "V_C6-bulk") {
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex]    = 0.0;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = capacityLeftSpeciesCoeff;
                found = true;
            }

            if (!found) {
                printf(":setCapacityCoeff_MCMB(): unknown species: %s  in phase  %s Setting Capacity to zero\n",
                       sss.c_str(), pname.c_str());
                capacityLeftSpeciesCoeff_[iGlobSpeciesIndex]    = 0.0;
                capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = 0.0;
            }
        }
    }
}
//====================================================================================================================
} // End of namespace Cantera
//======================================================================================================================
