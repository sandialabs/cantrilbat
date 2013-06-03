/**
 * @file m1d_BDT_porousLiKCl.cpp
 */

/*
 *  $Id: m1d_BDT_porousLiKCl.cpp 359 2012-08-16 23:34:13Z hkmoffa $
 */

#include "m1d_BDT_porousLiKCl.h"

#include "m1d_porousLiKCl_dom1D.h"

#include "m1d_ProblemStatementCell.h"
extern m1d::ProblemStatementCell PSinput;

using namespace Cantera;

namespace m1d
{
//======================================================================================================================
//======================================================================================================================
//======================================================================================================================
BDT_porousLiKCl::BDT_porousLiKCl(DomainLayout *dl_ptr) :
  BulkDomainDescription(dl_ptr), ionicLiquid_(0), trans_(0)
{

  IsAlgebraic_NE.resize(6,0);
  IsArithmeticScaled_NE.resize(6, 0);
  //  ionicLiquid_ = new Cantera::IonsFromNeutralVPSSTP("LiKCl_recipmoltenSalt_trans.xml");
  //ionicLiquid_ = new Cantera::IonsFromNeutralVPSSTP( PSinput.electrolyteFile_ );
  int iph = (PSinput.PhaseList_)->globalPhaseIndex(PSinput.electrolytePhase_);
  if (iph < 0) {
    throw CanteraError("BDT_porousLiKCl::BDT_porousLiKCl()", 
                       "Can't find the phase named " + PSinput.electrolytePhase_);
  }
  ThermoPhase* tmpPhase = & (PSinput.PhaseList_)->thermo(iph);
  ionicLiquid_ = dynamic_cast<Cantera::IonsFromNeutralVPSSTP *>( tmpPhase->duplMyselfAsThermoPhase() );

  trans_ = Cantera::newTransportMgr("Liquid", ionicLiquid_, 1);

  EquationNameList.clear();

  // Continuity is used to solve for bulk velocity
  // Note that this is a single phase continuity 
  // so phase change will result in a source term
  EquationNameList.push_back(EqnType(Continuity, 0, "Continuity: Bulk Velocity"));
  IsAlgebraic_NE[0] = 1;
  // the velocity can cross the origin with impunity
  IsArithmeticScaled_NE[0] = 1;

  VariableNameList.push_back(VarType(Velocity_Axial, 0, 0));
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

  //Current conservation is used to solve for electrostatic potential
  EquationNameList.push_back(EqnType(Current_Conservation, 0, "Current Conservation"));
  VariableNameList.push_back(VarType(Voltage, 0, "Electrolyte"));
  IsAlgebraic_NE[4] = 1;
  IsArithmeticScaled_NE[4] = 1;

  // Enthalpy conservation is used to solve for the temperature
  // EquationNameList.push_back(EqnType(Enthalpy_conservation, 0, "Enthalpy Conservation"));
}
//======================================================================================================================
BDT_porousLiKCl::BDT_porousLiKCl(const BDT_porousLiKCl &r) :
  BulkDomainDescription(r.DL_ptr_), ionicLiquid_(0), trans_(0)
{
  *this = r;
}
//======================================================================================================================
BDT_porousLiKCl::~BDT_porousLiKCl()
{
  /*
   * Delete the objects that we own
   */
  delete ionicLiquid_;
  delete trans_;
}
//======================================================================================================================
BDT_porousLiKCl &
BDT_porousLiKCl::operator=(const BDT_porousLiKCl &r)
{
  if (this == &r) {
    return *this;
  }

  BulkDomainDescription::operator=(r);

  EquationID = r.EquationID;

  delete ionicLiquid_;
  ionicLiquid_ = new Cantera::IonsFromNeutralVPSSTP(*(r.ionicLiquid_));

  delete trans_;
  trans_ = Cantera::newTransportMgr("Liquid", ionicLiquid_, 1);

  return *this;
}
//======================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
BulkDomain1D *
BDT_porousLiKCl::mallocDomain1D()
{
  BulkDomainPtr_ = new porousLiKCl_dom1D(*this);
  return BulkDomainPtr_;
}
//======================================================================================================================
} /* End of Namespace */
//======================================================================================================================

