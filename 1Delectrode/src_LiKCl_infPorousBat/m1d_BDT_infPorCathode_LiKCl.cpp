/**
 * @file m1d_BDT_infPorCathode_LiKCl.cpp
 */

#include "m1d_BDT_infPorCathode_LiKCl.h"
#include "m1d_infPorousLiKCl_FeS2Cathode_dom1D.h"
#include "m1d_exception.h"
#include "m1d_defs.h"

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
BDT_infPorCathode_LiKCl::BDT_infPorCathode_LiKCl(DomainLayout *dl_ptr) :
    BDD_porousElectrode(dl_ptr, 1), 
    ionicLiquidIFN_(0), 
    m_position(1)
{ 
  IsAlgebraic_NE.resize(7,0);  
  IsArithmeticScaled_NE.resize(7,0);




}
//=====================================================================================================================
BDT_infPorCathode_LiKCl::BDT_infPorCathode_LiKCl(const BDT_infPorCathode_LiKCl &r) :
    BDD_porousElectrode(r.DL_ptr_, 1 ), 
    ionicLiquidIFN_(0), 
    m_position(1)
{
  *this = r;
}
//=====================================================================================================================
BDT_infPorCathode_LiKCl::~BDT_infPorCathode_LiKCl()
{
  /*
   * Delete objects that we own
   */
  ionicLiquidIFN_ = 0;
  //safeDelete(Electrode_);
}
//=====================================================================================================================
BDT_infPorCathode_LiKCl &
BDT_infPorCathode_LiKCl::operator=(const BDT_infPorCathode_LiKCl &r)
{
  if (this == &r) {
    return *this;
  }

  BDD_porousElectrode::operator=(r);
  ionicLiquidIFN_ = (IonsFromNeutralVPSSTP *) ionicLiquid_;

  m_position = r.m_position;

  return *this;

}
//=====================================================================================================================
void
BDT_infPorCathode_LiKCl::ReadModelDescriptions()
{
    BDD_porousElectrode::ReadModelDescriptions();

    ionicLiquidIFN_ = dynamic_cast<Cantera::IonsFromNeutralVPSSTP *>( ionicLiquid_ );

   
}
//=====================================================================================================================
void
BDT_infPorCathode_LiKCl::SetEquationsVariablesList()
{
    int eqnIndex = 0;
    EquationNameList.clear();
    VariableNameList.clear();
 

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
    
    // Current conservation is used to solve for electrostatic potential
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
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
BulkDomain1D *
BDT_infPorCathode_LiKCl::mallocDomain1D()
{
  BulkDomainPtr_ = new infPorousLiKCl_FeS2Cathode_dom1D(this);
  return BulkDomainPtr_;
}

//=====================================================================================================================
void
BDT_infPorCathode_LiKCl::DetermineConstitutiveModels()
{
    if (!trans_) {
	delete trans_;
    }
    /*
     *  Create and Store a pointer to the Transport Manager
     */
    trans_ = Cantera::newTransportMgr("Liquid", ionicLiquidIFN_, 1);
}
//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================
