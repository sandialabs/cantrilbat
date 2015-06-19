/**
 * @file m1d_BDT_porousLiKCl.cpp
 */

/*
 *  $Id: m1d_BDT_porousLiKCl.cpp 361 2012-08-21 00:39:02Z hkmoffa $
 */

#include "m1d_BDT_porousLiKCl.h"

#include "m1d_porousLiKCl_infPlate_dom1D.h"

#include "m1d_ProblemStatementCell.h"
extern m1d::ProblemStatementCell PSinput;

using namespace Cantera;

namespace m1d
{

//====================================================================
//====================================================================
//====================================================================
BDT_porousLiKCl::BDT_porousLiKCl(DomainLayout *dl_ptr) :
    BDD_porousFlow(dl_ptr), 
    ionicLiquidIFN_(0)
{
  int eqnIndex = 0;
  IsAlgebraic_NE.resize(7,0);
  IsArithmeticScaled_NE.resize(7,0);

  //  ionicLiquid_ = new Cantera::IonsFromNeutralVPSSTP("LiKCl_recipmoltenSalt_trans.xml");
  int iph = (PSinput.PhaseList_)->globalPhaseIndex(PSinput.electrolytePhase_);
  ThermoPhase* tmpPhase = &(PSinput.PhaseList_)->thermo(iph);
  ionicLiquidIFN_ = dynamic_cast<Cantera::IonsFromNeutralVPSSTP *>( tmpPhase->duplMyselfAsThermoPhase() );
  
  if (trans_) {
    delete trans_;
  }
  trans_ = Cantera::newTransportMgr("Liquid", ionicLiquidIFN_, 1);

  EquationNameList.clear();

  // Continuity is used to solve for bulk velocity
  // Note that this is a single phase continuity 
  // so phase change will result in a source term
  EquationNameList.push_back(EqnType(Continuity, 0, "Continuity: Bulk Velocity"));

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
  eqnIndex++;

  EquationNameList.push_back(EqnType(MoleFraction_Summation, 0));
  eqnIndex++;

  EquationNameList.push_back(EqnType(ChargeNeutrality_Summation, 0));
  eqnIndex++;

  //Current conservation is used to solve for electrostatic potential
  EquationNameList.push_back(EqnType(Current_Conservation, 0, "Current Conservation"));
  VariableNameList.push_back(VarType(Voltage, 0, "Electrolyte"));
  IsAlgebraic_NE[eqnIndex] = 1;
  IsArithmeticScaled_NE[eqnIndex] = 1;
  eqnIndex++;

  //Enthalpy conservation is used to solve for the temperature
  // EquationNameList.push_back(EqnType(Enthalpy_conservation, 0, "Enthalpy Conservation"));

}
//==================================================================
BDT_porousLiKCl::BDT_porousLiKCl(const BDT_porousLiKCl &r) :
    BDD_porousFlow(r),
    ionicLiquidIFN_(0)
{
  *this = r;
}
//==================================================================
BDT_porousLiKCl::~BDT_porousLiKCl()
{
  /*
   * Delete objects that we own
   */
  delete ionicLiquidIFN_; 
}
//=====================================================================================================================
//  Make list of the equations and variables
/*
 *  We also set the ordering here.
 */
void
BDT_porousLiKCl::SetEquationsVariablesList()
{
    int eqnIndex = 0;
    EquationNameList.clear();
    VariableNameList.clear();

    // Continuity is used to solve for bulk velocity
    // Note that this is a single phase continuity 
    // so phase change will result in a source term
    EquationNameList.push_back(EqnType(Continuity, 0, "Continuity: Bulk Velocity"));
    
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
    eqnIndex++;
    
    EquationNameList.push_back(EqnType(MoleFraction_Summation, 0));
    eqnIndex++;
    
    EquationNameList.push_back(EqnType(ChargeNeutrality_Summation, 0));
    eqnIndex++;
    
    //Current conservation is used to solve for electrostatic potential
    EquationNameList.push_back(EqnType(Current_Conservation, 0, "Current Conservation"));
    VariableNameList.push_back(VarType(Voltage, 0, "Electrolyte"));
    IsAlgebraic_NE[eqnIndex] = 1;
    IsArithmeticScaled_NE[eqnIndex] = 1;
    eqnIndex++;
    
    //Enthalpy conservation is used to solve for the temperature
    // EquationNameList.push_back(EqnType(Enthalpy_conservation, 0, "Enthalpy Conservation"));
}
//==================================================================
BDT_porousLiKCl &
BDT_porousLiKCl::operator=(const BDT_porousLiKCl &r)
{
  if (this == &r) {
    return *this;
  }

  BulkDomainDescription::operator=(r);

  EquationID = r.EquationID;
  
  delete ionicLiquidIFN_;
  ionicLiquidIFN_ = new Cantera::IonsFromNeutralVPSSTP(*(r.ionicLiquidIFN_));

  setupTransport();

  return *this;
}
//=====================================================================================================================
//  Make list of the equations and variables
/*
 *  We also set the ordering here.
 */
void
BDT_porousLiKCl::setupTransport()
{
    delete trans_;
    trans_ = Cantera::newTransportMgr("Liquid", ionicLiquidIFN_, 1);
}
//==================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
BulkDomain1D *
BDT_porousLiKCl::mallocDomain1D()
{
  BulkDomainPtr_ = new porousLiKCl_infPlate_dom1D(*this);
  return BulkDomainPtr_;
}
//=====================================================================================================================
void
BDT_porousLiKCl::DetermineConstitutiveModels()
{
    setupTransport();
}
//==================================================================
} /* End of Namespace */
//==================================================================

