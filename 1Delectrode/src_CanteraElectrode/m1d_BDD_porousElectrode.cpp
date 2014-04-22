/**
 * @file m1d_BDT_porAnode_LiKCl.cpp
 */

/*
 *  $Id: m1d_BDD_porousElectrode.cpp 552 2013-03-01 21:25:03Z hkmoffa $
 */

#include "m1d_BDD_porousElectrode.h"
#include "m1d_porousElectrode_dom1D.h"
//#include "Electrode_Factory.h"
#include "m1d_ProblemStatementCell.h"
extern m1d::ProblemStatementCell PSinput;

using namespace std;
using namespace Cantera;

namespace m1d
{

//====================================================================================================================
  BDD_porousElectrode::BDD_porousElectrode(DomainLayout *dl_ptr,
					   std::string domainName, 
					   ELECTRODE_KEY_INPUT *electrode_input ) :
  BulkDomainDescription(dl_ptr, domainName),
  ionicLiquid_(0),
  trans_(0), 
  Electrode_(0)
{
  int eqnIndex = 0;
  IsAlgebraic_NE.resize(7,0);
  IsArithmeticScaled_NE.resize(7,0);
  /*
   * Store a copy of the electrolyte ThermoPhase object
   */
  int iph = (PSinput.PhaseList_)->globalPhaseIndex(PSinput.electrolytePhase_);
  if (iph < 0) {
    throw CanteraError("BDD_porousElectrode::BDD_porousElectrode()",
                       "Can't find the phase in the phase list: " + PSinput.electrolytePhase_);
  }
  ThermoPhase* tmpPhase = & (PSinput.PhaseList_)->thermo(iph);
  ionicLiquid_ = tmpPhase->duplMyselfAsThermoPhase();

  /*
   *  Create and Store a pointer to the Transport Manager
   */
  trans_ = Cantera::newTransportMgr(ionicLiquid_, 1);

  /*
   *  Find the hook into the input for the electrode model
   */
  ProblemStatementCell *psc_ptr = &PSinput;
  ELECTRODE_KEY_INPUT *ai = psc_ptr->anode_input_;

  /*
   *  Use the ElectrodeModelName value as input to the electrode factory to create the electrode
   */
  Electrode_  = newElectrodeObject(electrode_input->electrodeModelName);
  if (!Electrode_) {
    throw  m1d_Error("BDD_porousElectrode::instantiateElectrodeCells()",
		     "Electrode factory method failed");
  }
  /*
   *   Initialize the electrode model using the input from the ELECTRODE_KEY_INPUT object
   */
  retn = Electrode_->electrode_model_create(electrode_input_);
  if (retn == -1) {
    throw CanteraError("BDD_porousElectrode::BDD_porousElectrode()",
		       "Error initializing the anode electrode object");
  }
  retn = Electrode_->setInitialConditions(PSinput.anode_input_);
  if (retn == -1) {
    throw CanteraError("BDD_porousElectrode::BDD_porousElectrode()",
		       "Electrode::setInitialConditions() failed");
  }

  delete cfA;

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
  int nsp = ionicLiquid_->nSpecies();

  /*
   *  Loop over the species in the electrolyte phase. Each gets its own equation.
   *  Here, we hard code the mole fraction sum equation to the solvent, ECDMC, and we
   *  hardcode the charge conservation equation to PF6m. All other species are assigned
   *  the species conservation equation.
   */
  int iMFS = -1;
  int iCN = -1;
  for (int k = 0; k < nsp; k++) {
     if (namesSp[k] == "ECDMC") {
       iMFS = k;
       VariableNameList.push_back(VarType(MoleFraction_Species, 0, (namesSp[k]).c_str()));
       EquationNameList.push_back(EqnType(MoleFraction_Summation, 0));
       IsAlgebraic_NE[1 + k] = 2;
     } else if (namesSp[k] == "PF6-") {
       iCN = k;
       VariableNameList.push_back(VarType(MoleFraction_Species, 0, (namesSp[k]).c_str()));
       EquationNameList.push_back(EqnType(ChargeNeutrality_Summation, 0));
       IsAlgebraic_NE[1 + k] = 2;
     } else {
       VariableNameList.push_back(VarType(MoleFraction_Species, 0, (namesSp[k]).c_str()));
       EquationNameList.push_back(EqnType(Species_Conservation, 0, (namesSp[k]).c_str()));
       IsAlgebraic_NE[1 + k] = 0;
     }
     eqnIndex++;
  }
  if (iMFS < 0) {
    throw CanteraError("sep", "no ECDMC");
  }
  if (iCN < 0) {
    throw CanteraError("sep", "no PF6-");
  }



  //   Current conservation is used to solve for electrostatic potential
  //           Equation 4: Current Conservation - Electrolyte   Variable 4: Volts_Electrolyte
  EquationNameList.push_back(EqnType(Current_Conservation, 0, "Current Conservation"));
  VariableNameList.push_back(VarType(Voltage, 0, "Electrolyte"));
  IsAlgebraic_NE[eqnIndex] = 1;
  IsArithmeticScaled_NE[eqnIndex] = 1;
  eqnIndex++;

  // Current conservation is used to solve for electrostatic potential
  //  Equation 5: Current Conservation - Cathode   Variable 5: Volts_cathode

  EquationNameList.push_back(EqnType(Current_Conservation, 1, "Anode Current Conservation"));
  VariableNameList.push_back(VarType(Voltage, 1, "AnodeVoltage"));
  IsAlgebraic_NE[eqnIndex] = 1;
  IsArithmeticScaled_NE[eqnIndex] = 1;
  eqnIndex++;

  // Enthalpy conservation is used to solve for the temperature
  // EquationNameList.push_back(EqnType(Enthalpy_conservation, 0, "Enthalpy Conservation"));
}
//=====================================================================================================================
BDD_porousElectrode::BDD_porousElectrode(const BDD_porousElectrode &r) :
  BulkDomainDescription(r.DL_ptr_), ionicLiquid_(0), trans_(0), Electrode_(0)
{
  *this = r;
}
//=====================================================================================================================
BDD_porousElectrode::~BDD_porousElectrode()
{
  /*
   * Delete objects that we own
   */
  safeDelete(ionicLiquid_);
  safeDelete(trans_);
  safeDelete(Electrode_);
}
//=====================================================================================================================
BDD_porousElectrode &
BDD_porousElectrode::operator=(const BDD_porousElectrode &r)
{
  if (this == &r) {
    return *this;
  }
  BulkDomainDescription::operator=(r);

  delete ionicLiquid_;
  ionicLiquid_ = (r.ionicLiquid_)->duplMyselfAsThermoPhase();

  delete trans_;
  trans_ = Cantera::newTransportMgr("Simple", ionicLiquid_, 1);

  delete Electrode_;
  // ok this is wrong and needs to be changed
  Electrode_ = r.Electrode_->duplMyselfAsElectrode();

  exit(-1);

  return *this;
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
  BulkDomainPtr_ = new porousLiIon_Anode_dom1D(*this);
  return BulkDomainPtr_;
}
//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================
