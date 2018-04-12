/**
 * @file m1d_BDD_porousLiKCl.cpp
 */



#include "m1d_BDD_porousLiKCl.h"

#include "m1d_porousLiKCl_dom1D.h"

#include "m1d_ProblemStatementCell.h"
#include "m1d_CanteraElectrodeGlobals.h"
#include "m1d_exception.h"
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

extern m1d::ProblemStatementCell PSinput;


namespace m1d
{
//======================================================================================================================
//======================================================================================================================
//======================================================================================================================
BDD_porousLiKCl::BDD_porousLiKCl(DomainLayout *dl_ptr, std::string domainName) :
    BDD_porousFlow(dl_ptr, domainName),
    ionicLiquidIFN_(0)
{

  IsAlgebraic_NE.resize(6,0);
  IsArithmeticScaled_NE.resize(6, 0);
}
//======================================================================================================================
BDD_porousLiKCl::BDD_porousLiKCl(const BDD_porousLiKCl &r) :
    BDD_porousFlow(r.DL_ptr_), 
    ionicLiquidIFN_(0)
{
  *this = r;
}
//======================================================================================================================
BDD_porousLiKCl::~BDD_porousLiKCl()
{
  /*
   * Delete the objects that we own
   */
  //delete ionicLiquidIFN_;
}
//======================================================================================================================
BDD_porousLiKCl &
BDD_porousLiKCl::operator=(const BDD_porousLiKCl &r)
{
  if (this == &r) {
    return *this;
  }

  BDD_porousFlow::operator=(r);

  EquationID = r.EquationID;

  // delete ionicLiquidIFN_;
  //ionicLiquidIFN_ = new ZZCantera::IonsFromNeutralVPSSTP(*(r.ionicLiquidIFN_));
  ionicLiquidIFN_ = (IonsFromNeutralVPSSTP*) ionicLiquid_;

  // delete trans_;
  //trans_ = ZZCantera::newTransportMgr("Liquid", ionicLiquid_, 1);

  return *this;
}
//=====================================================================================================================
void
BDD_porousLiKCl::ReadModelDescriptions()
{

    BDD_porousFlow::ReadModelDescriptions();

 
     ionicLiquidIFN_ = dynamic_cast<ZZCantera::IonsFromNeutralVPSSTP *>( ionicLiquid_ );
     if (!ionicLiquidIFN_) {
	 throw m1d_Error("BDD_porousLiKCl::ReadModelDescriptions()", 
			 "ionicLiquidIFN_  failed on dynamic cast");
     }
}
//======================================================================================================================
void
BDD_porousLiKCl::SetEquationsVariablesList()
{
  /*
   *  Create a vector of Equation Names
   *  This is the main place to specify the ordering of the equations within the code
   */
  EquationNameList.clear();
  VariableNameList.clear();
  IsAlgebraic_NE.resize(6,0);
  IsArithmeticScaled_NE.resize(6, 0);

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
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
BulkDomain1D *
BDD_porousLiKCl::mallocDomain1D()
{
  BulkDomainPtr_ = new porousLiKCl_dom1D(this);
  return BulkDomainPtr_;
}
//=====================================================================================================================
void
BDD_porousLiKCl::DetermineConstitutiveModels()
{
    if (!trans_) {
	delete trans_;
    }
    /*
     *  Create and Store a pointer to the Transport Manager
     */
    trans_ = ZZCantera::newTransportMgr("Liquid", ionicLiquidIFN_, 1);
}
//======================================================================================================================
} /* End of Namespace */
//======================================================================================================================

