/**
 * @file m1d_BDT_porAnode_LiKCl.cpp
 */
/*
 *  $Id: m1d_BDD_porousElectrode.cpp 552 2013-03-01 21:25:03Z hkmoffa $
 */
#include "m1d_defs.h"
#include "m1d_BDD_porousElectrode.h"
#include "m1d_porousElectrode_dom1D.h"
#include "m1d_ProblemStatementCell.h"
#include "m1d_CanteraElectrodeGlobals.h" 
#include "Electrode_Factory.h"

#include "m1d_DomainLayout.h"

#include "m1d_BatteryResidEval.h"

#include "m1d_exception.h"

using namespace std;

namespace m1d
{
//====================================================================================================================
BDD_porousElectrode::BDD_porousElectrode(DomainLayout *dl_ptr, int electrodeType,
					 std::string domainName) :
  BDD_porousFlow(dl_ptr, domainName),
  Electrode_(0),
  metalPhase_(0),
  electrodeType_(electrodeType)
{
  IsAlgebraic_NE.resize(7,0);
  IsArithmeticScaled_NE.resize(7,0);
  /*
   * Store a copy of the electrolyte ThermoPhase object
   */

}
//=====================================================================================================================
BDD_porousElectrode::BDD_porousElectrode(const BDD_porousElectrode &r) :
    BDD_porousFlow(r), 
    Electrode_(0),
    metalPhase_(0)
{
  *this = r;
}
//=====================================================================================================================
BDD_porousElectrode::~BDD_porousElectrode()
{
  /*
   * Delete objects that we own
   */
  safeDelete(Electrode_);
  safeDelete(metalPhase_);
}
//=====================================================================================================================
BDD_porousElectrode &
BDD_porousElectrode::operator=(const BDD_porousElectrode &r)
{
  if (this == &r) {
    return *this;
  }
  BDD_porousFlow::operator=(r);

  safeDelete(Electrode_);
  // ok this is wrong and needs to be changed
  Electrode_ = r.Electrode_->duplMyselfAsElectrode();

  electrodeType_ = r.electrodeType_;

  safeDelete(metalPhase_);
  metalPhase_ = r.metalPhase_->duplMyselfAsThermoPhase();

  return *this;
}
//=====================================================================================================================
void
BDD_porousElectrode::ReadModelDescriptions()
{
    //
    // pick up the default setup of ionicLiquid_, solidSkeleton_.
    // Many cases probably don't need to read solidSkeleton_
    //
    BDD_porousFlow::ReadModelDescriptions();


    /*
     *  Find the hook into the input for the cathode electrode model
     */
    ProblemStatementCell* psc_ptr = PSCinput_ptr;

    ELECTRODE_KEY_INPUT* eki = psc_ptr->anode_input_;
    if (electrodeType_ == 1) {
	eki = psc_ptr->cathode_input_;
    } 
    if (electrodeType_ > 1) {
	throw  m1d_Error("BDD_porousElectrode::ReadModelDescription()",
                         "unhandled");
    }
    /*
     *  Use the ElectrodeModelName value as input to the electrode factory to create the electrode
     */
    Electrode_  = newElectrodeObject(eki->electrodeModelName);
    if (!Electrode_) {
        throw  m1d_Error("BDD_porousElectrode::ReadModelDescriptions()",
                         "newElectrodeObject failed");
    }
    ELECTRODE_KEY_INPUT* eki_new = newElectrodeKeyInputObject(eki->electrodeModelName);
    string commandFile = eki->commandFile_;
    BEInput::BlockEntry* cfA = new BEInput::BlockEntry("command_file");
    /*
     *  Parse the complete child input file
     */
    int retn = eki_new->electrode_input_child(commandFile, cfA);
    if (retn == -1) {
        throw  m1d_Error("BDD_porousElectrode::instantiateElectrodeCells()",
                         "Electrode input child method failed");
    }
    /*
     * Switch the pointers around so that the child input file is returned.
     * Delete the original pointer.
     */
    delete eki;
    if (electrodeType_ == 0) {
	PSCinput_ptr->anode_input_ = eki_new;
    } else if (electrodeType_ == 1) {
	PSCinput_ptr->cathode_input_ = eki_new;
    }

    /*
     *   Initialize the electrode model using the input from the ELECTRODE_KEY_INPUT object
     */
    if (electrodeType_ == 0) { 
	retn = Electrode_->electrode_model_create(PSCinput_ptr->anode_input_);
	if (retn == -1) {
	    throw CanteraError("BDD_porousElectrode::ReadModelDescriptions()",
			       "Error initializing the anode electrode object");
	}
	retn = Electrode_->setInitialConditions(PSCinput_ptr->anode_input_);
	if (retn == -1) {
	    throw CanteraError("BDD_porousElectrode::ReadModelDescriptions()",
			       "Electrode::setInitialConditions() for anode failed");
	}
    } else if (electrodeType_ == 1) {
	retn = Electrode_->electrode_model_create(PSCinput_ptr->cathode_input_);
	if (retn == -1) {
	    throw CanteraError("BDD_porousElectrode::ReadModelDescriptions()",
			       "Error initializing the cathode electrode object");
	}
	retn = Electrode_->setInitialConditions(PSCinput_ptr->cathode_input_);
	if (retn == -1) {
	    throw CanteraError("BDD_porousElectrode::ReadModelDescriptions()",
			       "Electrode::setInitialConditions() for cathode failed");
	}
    }

    int metalIndex = Electrode_->metalPhaseIndex();
    ThermoPhase*mP = & Electrode_->thermo(metalIndex);
    metalPhase_ = mP->duplMyselfAsThermoPhase();
    string ss = metalPhase_->id();
    //printf("name of metal phase = %s\n", ss.c_str());

    delete cfA;
}
//=====================================================================================================================
//  Make list of the equations and variables
/*
 *  We also set the ordering here.
 */
void
BDD_porousElectrode::SetEquationsVariablesList()
{
    //
    //  Get the Problem object. This will have information about what type of equations are to be solved
    //
    ProblemResidEval *pb = DL_ptr_->problemResid_;
    BatteryResidEval *batres = dynamic_cast<BatteryResidEval*>(pb);
    if (!batres) {
	printf("contention not true\n");
	exit(-1);
    }

    int eqnIndex = 0;
    /*
     *  Create a vector of Equation Names
     *  This is the main place to specify the ordering of the equations within the code
     */
    EquationNameList.clear();
    VariableNameList.clear();


    // Continuity is used to solve for bulk velocity
    // Note that this is a single phase continuity so phase change will result in a source term
    //         Equation 0: = Continuity         variable 0 = Axial Velocity

    EquationNameList.push_back(EqnType(Continuity, 0, "Continuity: Bulk Velocity"));
    VariableNameList.push_back(VarType(Velocity_Axial, 0, 0));
    IsAlgebraic_NE[eqnIndex] = 1;
    IsArithmeticScaled_NE[eqnIndex] = 1;
    eqnIndex++;

    // List of species in the electrolyte
    const std::vector<std::string> & namesSp = ionicLiquid_->speciesNames();
    nSpeciesElectrolyte_ = ionicLiquid_->nSpecies();

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

    if (electrodeType_ == 0) {
	// Current conservation is used to solve for electrostatic potential
	//  Equation 5: Current Conservation - Cathode   Variable 5: Volts_cathode
	EquationNameList.push_back(EqnType(Current_Conservation, 1, "Anode Current Conservation"));
	VariableNameList.push_back(VarType(Voltage, 1, "AnodeVoltage"));
	IsAlgebraic_NE[eqnIndex] = 1;
	IsArithmeticScaled_NE[eqnIndex] = 1;
	eqnIndex++;

    } else if (electrodeType_ == 1) {
	// Current conservation is used to solve for electrostatic potential
	//  Equation 5: Current Conservation - Cathode   Variable 5: Volts_cathode
	EquationNameList.push_back(EqnType(Current_Conservation, 2, "Cathode Current Conservation"));
	VariableNameList.push_back(VarType(Voltage, 2, "CathodeVoltage"));
	IsAlgebraic_NE[eqnIndex] = 1;
	IsArithmeticScaled_NE[eqnIndex] = 1;
	eqnIndex++;
    } else {
	exit(-1);
    }
    //
    //  Enthalpy conservation is used to solve for the temperature
    //           EquationNameList.push_back(EqnType(Enthalpy_conservation, 0, "Enthalpy Conservation"));
    //           This is not an algebraic equation. It has a time derivative.
    if (PSCinput_ptr->Energy_equation_prob_type_ == 3) {
	EqnType ht = EqnType(Enthalpy_Conservation, 0, "Enthalpy Conservation");
	EquationNameList.push_back(ht);
	VariableNameList.push_back(VarType(Temperature, 0, 0));
	IsAlgebraic_NE[eqnIndex] = 0;
	IsArithmeticScaled_NE[eqnIndex] = 0;
	eqnIndex++;
    }
    //
    //  Thermal conservation is used to solve for the temperature
    //           
    if (PSCinput_ptr->Energy_equation_prob_type_ == 4) {
	EqnType ht = EqnType(Thermal_Conservation, 0, "Thermal_Conservation");
	EquationNameList.push_back(ht);
	VariableNameList.push_back(VarType(Temperature, 0, 0));
	IsAlgebraic_NE[eqnIndex] = 0;
	IsArithmeticScaled_NE[eqnIndex] = 0;
	eqnIndex++;
    }
    //
    //  Thermal dirichilet is used to solve for the temperature
    //       
    if (PSCinput_ptr->Energy_equation_prob_type_ == 2) {
        /*
         * For equation type 2 push a dirichilet equation to the rear
         */
        EquationNameList.push_back(EqnType(Thermal_Dirichilet, 0, "Thermal Dirichilet"));
        VariableNameList.push_back(VarType(Temperature, 0, 0));
        IsAlgebraic_NE[eqnIndex] = 0;
        IsArithmeticScaled_NE[eqnIndex] = 0;
        eqnIndex++;
    } 
    //
    //  Mechanical strain equation
    //
#ifdef MECH_MODEL
    if (batres->solidMechanicsProbType_ == 1) {
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
BDD_porousElectrode::mallocDomain1D()
{
    BulkDomainPtr_ = new porousElectrode_dom1D(*this);
    return BulkDomainPtr_;
}
//=====================================================================================================================
void
BDD_porousElectrode::DetermineConstitutiveModels()
{
    BDD_porousFlow::DetermineConstitutiveModels();
}
//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================
