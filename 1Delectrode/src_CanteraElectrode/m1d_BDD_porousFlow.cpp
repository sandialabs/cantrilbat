/**
 * @file m1d_BDT_porAnode_LiKCl.cpp
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_defs.h"
#include "m1d_BDD_porousFlow.h"
#include "m1d_porousFlow_dom1D.h"
#include "m1d_ProblemStatementCell.h"
#include "m1d_CanteraElectrodeGlobals.h" 
#include "m1d_BatteryResidEval.h"

#include "m1d_DomainLayout.h"

using namespace std;
using namespace Cantera;

namespace m1d
{

//====================================================================================================================
BDD_porousFlow::BDD_porousFlow(DomainLayout *dl_ptr,
			       std::string domainName) :
    BulkDomainDescription(dl_ptr, domainName),
    ionicLiquid_(0),
    trans_(0),
    solidSkeleton_(0),
    nSpeciesElectrolyte_(npos),
    iMFS_index_(npos),
    iCN_index_(npos)
{
    IsAlgebraic_NE.resize(7,0);
    IsArithmeticScaled_NE.resize(7,0);
}
//=====================================================================================================================
BDD_porousFlow::BDD_porousFlow(const BDD_porousFlow &r) :
    BulkDomainDescription(r.DL_ptr_), 
    ionicLiquid_(0), 
    trans_(0),
    solidSkeleton_(0),
    nSpeciesElectrolyte_(npos),
    iMFS_index_(npos),
    iCN_index_(npos)
{
    *this = r;
}
//=====================================================================================================================
BDD_porousFlow::~BDD_porousFlow()
{
  /*
   * Delete objects that we own
   */
  safeDelete(ionicLiquid_);

  safeDelete(trans_);

  safeDelete(solidSkeleton_);
 
}
//=====================================================================================================================
BDD_porousFlow &
BDD_porousFlow::operator=(const BDD_porousFlow &r)
{
  if (this == &r) {
    return *this;
  }
  BulkDomainDescription::operator=(r);

  safeDelete(ionicLiquid_);
  if (r.ionicLiquid_) {
      ionicLiquid_ = (r.ionicLiquid_)->duplMyselfAsThermoPhase();
  }
  safeDelete(trans_);
  if (r.trans_) {
       trans_ = (r.trans_)->duplMyselfAsTransport();
  }
  safeDelete(solidSkeleton_);
  if (r.solidSkeleton_) {
       solidSkeleton_ = (r.solidSkeleton_)->duplMyselfAsThermoPhase();
  }

  nSpeciesElectrolyte_  = r.nSpeciesElectrolyte_;
  iMFS_index_           = r.iMFS_index_;
  iCN_index_            = r.iCN_index_;

  return *this;
}
//=====================================================================================================================
void
BDD_porousFlow::ReadModelDescriptions()
{
    /*
     * Store a copy of the electrolyte ThermoPhase object
     */
    int iph = (PSCinput_ptr->PhaseList_)->globalPhaseIndex(PSCinput_ptr->electrolytePhase_);
    if (iph < 0) {
      throw CanteraError("BDD_porousFlow::ReadModelDescriptions()",
                         "Can't find the phase in the phase list: " + PSCinput_ptr->electrolytePhase_);
    }
    ThermoPhase* tmpPhase = & (PSCinput_ptr->PhaseList_)->thermo(iph);
    ionicLiquid_ = tmpPhase->duplMyselfAsThermoPhase();
    nSpeciesElectrolyte_ =  ionicLiquid_->nSpecies();
    
    iph = (PSCinput_ptr->PhaseList_)->globalPhaseIndex(PSCinput_ptr->separatorPhase_);
    if (iph < 0) {
      throw CanteraError("BDD_porousFlow::ReadModelDescriptions()",
                         "Can't find the phase in the phase list: " + PSCinput_ptr->separatorPhase_);
    }
    tmpPhase = & (PSCinput_ptr->PhaseList_)->thermo(iph);
    if (!tmpPhase) {
        throw CanteraError("BDD_porousFlow::ReadModelDescriptions()",
                           "Can't find the ThermoPhase in the phase list: " + PSCinput_ptr->separatorPhase_);
    }
    solidSkeleton_ = tmpPhase->duplMyselfAsThermoPhase();
}
//=====================================================================================================================
//  Make list of the equations and variables
/*
 *  We also set the ordering here.
 */
void
BDD_porousFlow::SetEquationsVariablesList()
{
    int eqnIndex = 0;
    //
    //  Get the Problem object. This will have information about what type of equations are to be solved
    //
    ProblemResidEval *pb = DL_ptr_->problemResid_;
    BatteryResidEval* batRes = static_cast<BatteryResidEval*>(pb);

    if (!pb) {
	throw m1d_Error("BDD_porousFlow::setEquationsVariablesList()",
			"ProblemResidEval not set yet");
    }
    if (!batRes) {
	throw m1d_Error("BDD_porousFlow::setEquationsVariablesList()",
			"problem isn't a battery problem");
    }
    /*
     *  Clear the list. We may have set the list previously
     */
    EquationNameList.clear();
    VariableNameList.clear();

    if (batRes->hasPressureEquation_ == 0) {
	// Continuity is used to solve for bulk velocity
	// Note that this is a single phase continuity so phase change will result in a source term
	//         Equation 0: = Continuity         variable 0 = Axial Velocity
	
	EquationNameList.push_back(EqnType(Continuity, 0, "Continuity: Bulk Velocity"));
	VariableNameList.push_back(VarType(Velocity_Axial, 0, 0));
	IsAlgebraic_NE[eqnIndex] = 1;
	IsArithmeticScaled_NE[eqnIndex] = 1;
	eqnIndex++;
    } else {
	// Continuity is used to solve for pressure
	//
	//         Equation 0: = Continuity         variable 0 = Pressure
	//         Equation 1: = Darcy's equation   variable 1 = Axial Velocity
	EquationNameList.push_back(EqnType(Continuity, 0, "Continuity: Pressure"));
	VariableNameList.push_back(VarType(Pressure_Axial, 0, 0));
	IsAlgebraic_NE[eqnIndex] = 0;          // Pressure has a time derivative!
	IsArithmeticScaled_NE[eqnIndex] = 1;   // liquid Pressure may not be positive -> think cavitation. 
	eqnIndex++;

	EquationNameList.push_back(EqnType(Momentum_Axial, 0, "Electrolyte Darcy Velocity"));
	VariableNameList.push_back(VarType(Velocity_Axial, 0, 0));
	IsAlgebraic_NE[eqnIndex] = 1;          // Velocity doesn't have a time derivative -> It is instantaneous
	                                       //  This means it's part of the DAE problem
	IsArithmeticScaled_NE[eqnIndex] = 1;   // velocity may be pos or neg 
	eqnIndex++;
    }
    
    // List of species in the electrolyte
    const std::vector<std::string> & namesSp = ionicLiquid_->speciesNames();
    nSpeciesElectrolyte_ =  ionicLiquid_->nSpecies();

    /*
     *  Loop over the species in the electrolyte phase. Each gets its own equation.
     *  Here, we hard code the mole fraction sum equation to the solvent, ECDMC, and we
     *  hardcode the charge conservation equation to PF6m. All other species are assigned
     *  the species conservation equation.
     */
    for (size_t k = 0; k <  nSpeciesElectrolyte_; k++) {
	if (namesSp[k] == "ECDMC") {
	    iMFS_index_ = k;
	    VariableNameList.push_back(VarType(MoleFraction_Species, k, (namesSp[k]).c_str()));
	    EquationNameList.push_back(EqnType(MoleFraction_Summation, 0));
	    IsAlgebraic_NE[1 + k] = 2;
	} else if (namesSp[k] == "PF6-") {
	    iCN_index_ = k;
	    VariableNameList.push_back(VarType(MoleFraction_Species, k, (namesSp[k]).c_str()));
	    EquationNameList.push_back(EqnType(ChargeNeutrality_Summation, 0));
	    IsAlgebraic_NE[1 + k] = 2;
	} else {
	    VariableNameList.push_back(VarType(MoleFraction_Species, k, (namesSp[k]).c_str()));
	    EquationNameList.push_back(EqnType(Species_Conservation, k, (namesSp[k]).c_str()));
	    IsAlgebraic_NE[1 + k] = 0;
	}
	eqnIndex++;
    }

    if (iMFS_index_ == npos) {
	iMFS_index_ = 0;
    }
    if (iCN_index_ == npos) {
	for (size_t k = 0; k <  nSpeciesElectrolyte_; k++) {
	    if (k != iMFS_index_) {
		double chg = ionicLiquid_->charge(k);
		if (chg < 0.0) {
		    iCN_index_ = k;
		    break;
		}
	    }
	}
	if (iCN_index_ == npos) {
	    throw CanteraError("sep", "no negative charge species");
	}
    }
    
    //   Current conservation is used to solve for electrostatic potential
    //           Equation 4: Current Conservation - Electrolyte   Variable 4: Volts_Electrolyte
    EquationNameList.push_back(EqnType(Current_Conservation, 0, "Current Conservation"));
    VariableNameList.push_back(VarType(Voltage, 0, "Electrolyte"));
    IsAlgebraic_NE[eqnIndex] = 1;
    IsArithmeticScaled_NE[eqnIndex] = 1;
    eqnIndex++;
    
    // Enthalpy conservation is used to solve for the temperature
    if (pb->energyEquationProbType_ == 3) {
	/*
	 * For equation type 3 push an enthalpy conservation equation to the rear.
	 */
	EquationNameList.push_back(EqnType(Enthalpy_Conservation, 0, "Enthalpy Conservation"));
	VariableNameList.push_back(VarType(Temperature, 0, 0));
	IsAlgebraic_NE[eqnIndex] = 0;
	IsArithmeticScaled_NE[eqnIndex] = 0;
	eqnIndex++;

    } else if (pb->energyEquationProbType_ == 4) {
        /*
         * For equation type 4 push an enthalpy conservation equation to the rear.
         */
        EquationNameList.push_back(EqnType(Thermal_Conservation, 0, "Thermal Conservation"));
        VariableNameList.push_back(VarType(Temperature, 0, 0));
        IsAlgebraic_NE[eqnIndex] = 0;
        IsArithmeticScaled_NE[eqnIndex] = 0;
        eqnIndex++;

    } else if (pb->energyEquationProbType_ == 2) {
        /*
         * For equation type 2 push a dirichilet equation to the rear
         */
        EquationNameList.push_back(EqnType(Thermal_Dirichilet, 0, "Thermal Dirichilet"));
        VariableNameList.push_back(VarType(Temperature, 0, 0));
        IsAlgebraic_NE[eqnIndex] = 0;
        IsArithmeticScaled_NE[eqnIndex] = 0;
        eqnIndex++;

    }    
    // CBL add the mechanical model solution vars here. 
#ifdef MECH_MODEL
    if (pb->solidMechanicsProbType_ == 1) {
        EquationNameList.push_back(EqnType(Mechanical_Model_Axial, 0, "Mech Strain"));
        VariableNameList.push_back(VarType(Solid_Stress_Axial, 0, 0));
        IsAlgebraic_NE[eqnIndex] = 0;
        IsArithmeticScaled_NE[eqnIndex] = 0;
        eqnIndex++;
    }
#endif
}
//=====================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual efficiently
 */
BulkDomain1D *
BDD_porousFlow::mallocDomain1D()
{
  BulkDomainPtr_ = new porousFlow_dom1D(*this);
  return BulkDomainPtr_;
}
//=====================================================================================================================
void
BDD_porousFlow::DetermineConstitutiveModels()
{
  if (!trans_) {
     delete trans_;
  }
  /*
   *  Create and Store a pointer to the Transport Manager
   */
  trans_ = Cantera::newTransportMgr("Simple", ionicLiquid_, 1);
}
//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================
