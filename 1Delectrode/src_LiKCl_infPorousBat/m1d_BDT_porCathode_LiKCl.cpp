/**
 * @file m1d_BDT_porCathode_LiKCl.cpp
 */

/*
 *  $Id: m1d_BDT_porCathode_LiKCl.cpp 598 2013-05-15 15:22:09Z hkmoffa $
 */

#include "m1d_BDT_porCathode_LiKCl.h"
#include "m1d_porousLiKCl_FeS2Cathode_dom1D.h"

#include "Electrode_input.h"
#include "Electrode_InfCapacity.h"
#include "Electrode_SimplePhaseChangeDiffusion.h"
#include "Electrode_Factory.h"

#include "m1d_ProblemStatementCell.h"
extern m1d::ProblemStatementCell PSinput;

using namespace std;
using namespace Cantera;

namespace m1d
{

//=====================================================================================================================
BDT_porCathode_LiKCl::BDT_porCathode_LiKCl(DomainLayout *dl_ptr) :
  BulkDomainDescription(dl_ptr), ionicLiquid_(0), trans_(0), m_position(1), Electrode_(0)
{ 
  int eqnIndex = 0;
  IsAlgebraic_NE.resize(7,0);  
  IsArithmeticScaled_NE.resize(7,0);
  //! initialize the ions from liquid
  //  ionicLiquid_ = new Cantera::IonsFromNeutralVPSSTP("LiKCl_recipmoltenSalt_trans.xml");
  // ionicLiquid_ = new Cantera::IonsFromNeutralVPSSTP( PSinput.electrolyteFile_ );
  int iph = (PSinput.PhaseList_)->globalPhaseIndex(PSinput.electrolytePhase_);
  if (iph < 0) {
    throw CanteraError("BDT_porCathode_LiKCl::BDT_porCathode_LiKCl",
                       "Can't find the phase in the phase list: " + PSinput.electrolytePhase_);
  }
  ThermoPhase* tmpPhase = & (PSinput.PhaseList_)->thermo(iph);
  ionicLiquid_ = dynamic_cast<Cantera::IonsFromNeutralVPSSTP *>( tmpPhase->duplMyselfAsThermoPhase() );

  trans_ = Cantera::newTransportMgr("Liquid", ionicLiquid_, 1);

  /*
   *   Initialize the electrode model
   */
  ProblemStatementCell *psc_ptr = &PSinput;
  ELECTRODE_KEY_INPUT *ci = psc_ptr->cathode_input_;
  Electrode_  = newElectrodeObject(ci->electrodeModelName);
  if (!Electrode_) {
      throw  m1d_Error("BDT_porCathode_LiKCl::BDT_porCathode_LiKCl()",
		       "Electrode factory method failed");
  }
  ELECTRODE_KEY_INPUT *ci_new = newElectrodeKeyInputObject(ci->electrodeModelName);  
  string commandFile = ci->commandFile_;
  BEInput::BlockEntry *cfC = new BEInput::BlockEntry("command_file");

  /*
   *  Parse the complete child input file
   */
  int retn = ci_new->electrode_input_child(commandFile, cfC);
  if (retn == -1) {
    throw  m1d_Error("BDT_porCathode_LiKCl::BDT_porCathode_LiKCl()",
                     "Electrode input child method failed");
  }
  /*
   * Switch the pointers around so that the child input file is returned.
   * Delete the original pointer.
   */
  delete ci;
  psc_ptr->cathode_input_ = ci_new;

  retn = Electrode_->electrode_model_create(PSinput.cathode_input_);
  if (retn == -1) {
    throw  m1d_Error("BDT_porCathode_LiKCl::BDT_porCathode_LiKCl()", 
                     "Electrode model create method failed");
  }
  retn = Electrode_->setInitialConditions(PSinput.cathode_input_);
  if (retn == -1) {
    throw  m1d_Error("BDT_porCathode_LiKCl::BDT_porCathode_LiKCl()", 
                     "Electrode::setInitialConditions method failed");
  }

  delete cfC;

  EquationNameList.clear();

  // Continuity is used to solve for bulk velocity
  // Note that this is a single phase continuity 
  // so phase change will result in a source term
  //         Equation 0: = Continuity         variable 0 = Axial Velocity

  EquationNameList.push_back(EqnType(Continuity, 0, "Continuity: Bulk Velocity"));
  VariableNameList.push_back(VarType(Velocity_Axial, 0, 0));
  IsAlgebraic_NE[0] = 1;
  IsArithmeticScaled_NE[0] = 1;
  eqnIndex++;

  // Equation 1: Species Conservation Li+        Variable 1: Li+ Mole Fraction
  // Equation 2: Species Conservation K+         Variable 2: K+  Mole Fraction
  // Equation 3: Species Conservation Cl-        Variable 3: Cl- Mole Fraction
  //list of species equations
  std::vector<std::string> namesSp;
  //  namesSp.resize(mp->nSpecies());
  //for (int i = 0 ; i < namesSp.size(); i++) { 
  // namesSp.push_back(mp->speciesName(i));
  //}
  /*
   *  Hard Code Names for the moment.
   */
  namesSp.push_back("Li+");
  namesSp.push_back("K+");
  namesSp.push_back("Cl-");
  VariableNameList.push_back(VarType(MoleFraction_Species, 0, (namesSp[0]).c_str()));
  VariableNameList.push_back(VarType(MoleFraction_Species, 1, (namesSp[1]).c_str()));
  VariableNameList.push_back(VarType(MoleFraction_Species, 2, (namesSp[2]).c_str()));

  // Species conservation is used to solve for mole fractions

  EquationNameList.push_back(EqnType(Species_Conservation, 0, (namesSp[0]).c_str()));
  EquationNameList.push_back(EqnType(MoleFraction_Summation, 0));
  IsAlgebraic_NE[2] = 2;
  EquationNameList.push_back(EqnType(ChargeNeutrality_Summation, 0));
  IsAlgebraic_NE[3] = 2;
  eqnIndex++;
  eqnIndex++;
  eqnIndex++;

  //Current conservation is used to solve for electrostatic potential
  //  Equation 4: Current Conservation - Electrolyte   Variable 4: Volts_Electrolyte

  EquationNameList.push_back(EqnType(Current_Conservation, 0, "Current Conservation"));
  VariableNameList.push_back(VarType(Voltage, 0, "Electrolyte"));
  IsAlgebraic_NE[4] = 1;
  IsArithmeticScaled_NE[eqnIndex] = 1;
  eqnIndex++;

  //Current conservation is used to solve for electrostatic potential
  //  Equation 5: Current Conservation - Cathode   Variable 5: CathodeVoltage

  EqnType cc = EqnType(Current_Conservation, 2, "Cathode Current Conservation");
  EquationNameList.push_back(cc);
  VariableNameList.push_back(VarType(Voltage, 2, "CathodeVoltage"));
  IsAlgebraic_NE[5] = 1;
  IsArithmeticScaled_NE[eqnIndex] = 1;
  eqnIndex++;

  //Enthalpy conservation is used to solve for the temperature

  // EquationNameList.push_back(EqnType(Enthalpy_conservation, 0, "Enthalpy Conservation"));

}
//=====================================================================================================================
BDT_porCathode_LiKCl::BDT_porCathode_LiKCl(const BDT_porCathode_LiKCl &r) :
  BulkDomainDescription(r.DL_ptr_), ionicLiquid_(0), trans_(0), m_position(1), Electrode_(0)
{
  *this = r;
}
//=====================================================================================================================
BDT_porCathode_LiKCl::~BDT_porCathode_LiKCl()
{
  /*
   * Delete objects that we own
   */
  safeDelete(ionicLiquid_);
  safeDelete(trans_);
  safeDelete(Electrode_);
}
//=====================================================================================================================
BDT_porCathode_LiKCl &
BDT_porCathode_LiKCl::operator=(const BDT_porCathode_LiKCl &r)
{
  if (this == &r) {
    return *this;
  }

  BulkDomainDescription::operator=(r);

  delete ionicLiquid_;
  ionicLiquid_ = new Cantera::IonsFromNeutralVPSSTP(*(r.ionicLiquid_));

  delete trans_;
  trans_ = Cantera::newTransportMgr("Liquid", ionicLiquid_, 1);

  m_position = r.m_position;

  delete Electrode_;
  Electrode_ = r.Electrode_->duplMyselfAsElectrode();

  return *this;
}
//=====================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
BulkDomain1D *
BDT_porCathode_LiKCl::mallocDomain1D()
{
  BulkDomainPtr_ = new porousLiKCl_FeS2Cathode_dom1D(*this);
  return BulkDomainPtr_;
}
//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================
