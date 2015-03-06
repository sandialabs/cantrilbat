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

#include "m1d_exception.h"

using namespace std;

namespace m1d
{

//====================================================================================================================
BDD_porousElectrode::BDD_porousElectrode(DomainLayout *dl_ptr, int electrodeType,
					 std::string domainName) :
  BDD_porousFlow(dl_ptr, domainName),
  ionicLiquid_(0),
  trans_(0), 
  Electrode_(0),
  electrodeType_(electrodeType)
{
  IsAlgebraic_NE.resize(7,0);
  IsArithmeticScaled_NE.resize(7,0);
  /*
   * Store a copy of the electrolyte ThermoPhase object
   */
  int iph = (PSCinput_ptr->PhaseList_)->globalPhaseIndex(PSCinput_ptr->electrolytePhase_);
  if (iph < 0) {
    throw CanteraError("BDD_porousElectrode::BDD_porousElectrode()",
                       "Can't find the phase in the phase list: " + PSCinput_ptr->electrolytePhase_);
  }
  ThermoPhase* tmpPhase = & (PSCinput_ptr->PhaseList_)->thermo(iph);
  ionicLiquid_ = tmpPhase->duplMyselfAsThermoPhase();

  /*
   *  Create and Store a pointer to the Transport Manager
   */
  trans_ = Cantera::newTransportMgr("Simple", ionicLiquid_, 1);


}
//=====================================================================================================================
BDD_porousElectrode::BDD_porousElectrode(const BDD_porousElectrode &r) :
  BDD_porousFlow(r.DL_ptr_), ionicLiquid_(0), trans_(0), Electrode_(0)
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

  electrodeType_ = r.electrodeType_;

  exit(-1);

  return *this;
}
//=====================================================================================================================
//  Make list of the equations and variables
/*
 *  We also set the ordering here.
 */
void
BDD_porousElectrode::SetEquationsVariablesList()
{
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
       VariableNameList.push_back(VarType(MoleFraction_Species, k, (namesSp[k]).c_str()));
       EquationNameList.push_back(EqnType(MoleFraction_Summation, 0));
       IsAlgebraic_NE[1 + k] = 2;
     } else if (namesSp[k] == "PF6-") {
       iCN = k;
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
  EquationNameList.push_back(EqnType(Current_Conservation, 2, "Catudne Current Conservation"));
  VariableNameList.push_back(VarType(Voltage, 2, "CathodeVoltage"));
  IsAlgebraic_NE[eqnIndex] = 1;
  IsArithmeticScaled_NE[eqnIndex] = 1;
  eqnIndex++;
  } else {
    exit(-1);
  }


  // Enthalpy conservation is used to solve for the temperature
  // EquationNameList.push_back(EqnType(Enthalpy_conservation, 0, "Enthalpy Conservation"));
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
} /* End of Namespace */
//=====================================================================================================================
