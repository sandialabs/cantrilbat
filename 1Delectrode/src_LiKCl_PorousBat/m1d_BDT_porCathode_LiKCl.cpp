/**
 * @file m1d_BDT_porCathode_LiKCl.cpp
 */

/*
 *  $Id: m1d_BDT_porCathode_LiKCl.cpp 598 2013-05-15 15:22:09Z hkmoffa $
 */

#include "m1d_BDT_porCathode_LiKCl.h"
#include "m1d_porousLiKCl_FeS2Cathode_dom1D.h"
#include "m1d_exception.h"
#include "m1d_defs.h"

#include "Electrode_input.h"
#include "Electrode_InfCapacity.h"
#include "Electrode_SimplePhaseChangeDiffusion.h"
#include "Electrode_Factory.h"

#include "m1d_ProblemStatementCell.h"
extern m1d::ProblemStatementCell PSinput;

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

namespace m1d
{

//=====================================================================================================================
BDT_porCathode_LiKCl::BDT_porCathode_LiKCl(DomainLayout *dl_ptr) :
    BDD_porousElectrode(dl_ptr, 1, ""), 
    ionicLiquidIFN_(0), 
    m_position(1)
{
  IsAlgebraic_NE.resize(7,0);
  IsArithmeticScaled_NE.resize(6, 0);
}
//=====================================================================================================================
BDT_porCathode_LiKCl::BDT_porCathode_LiKCl(const BDT_porCathode_LiKCl &r) :
    BDD_porousElectrode(r), 
    ionicLiquidIFN_(0),
    m_position(1)
{
    *this = r;
}
//=====================================================================================================================
BDT_porCathode_LiKCl::~BDT_porCathode_LiKCl()
{
  /*
   * Delete objects that we own
   */

}
//=====================================================================================================================
BDT_porCathode_LiKCl &
BDT_porCathode_LiKCl::operator=(const BDT_porCathode_LiKCl &r)
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
BDT_porCathode_LiKCl::ReadModelDescriptions()
{

    int iph = (PSinput.PhaseList_)->globalPhaseIndex(PSinput.electrolytePhase_);
    if (iph < 0) {
	throw CanteraError("BDT_porCathode_LiKCl::BDT_porAnode_LiKCl",
			   "Can't find the phase in the phase list: " + PSinput.electrolytePhase_);
    }
    ionicLiquid_ = & (PSinput.PhaseList_)->thermo(iph);
    ionicLiquidIFN_ = dynamic_cast<ZZCantera::IonsFromNeutralVPSSTP *>( ionicLiquid_->duplMyselfAsThermoPhase() );
    ionicLiquid_ = (ThermoPhase*) ionicLiquidIFN_;

    ELECTRODE_KEY_INPUT *ci = PSCinput_ptr->cathode_input_;
    Electrode_ = newElectrodeObject(ci->electrodeModelName);
    if (!Electrode_) {
	throw  m1d_Error("BDT_porCathode_LiKCl::ReadModelDescriptions()", "Electrode factory method failed");
    }
    ELECTRODE_KEY_INPUT *ci_new = newElectrodeKeyInputObject(ci->electrodeModelName);  
    string commandFile = ci->commandFile_;
    BEInput::BlockEntry *cfA = new BEInput::BlockEntry("command_file");
    
    /*
     *  Parse the complete child input file
     */
    int retn = ci_new->electrode_input_child(commandFile, cfA);
    if (retn == -1) {
	throw  m1d_Error("BDT_porCathode_LiKCl::ReadModelDescriptions()",
			 "Electrode input child method failed");
    }
    /*
     * Switch the pointers around so that the child input file is returned.
     * Delete the original pointer.
     */
    delete ci;
    PSCinput_ptr->cathode_input_ = ci_new;
  
    retn = Electrode_->electrode_model_create(PSCinput_ptr->cathode_input_);
    if (retn == -1) {
	throw  m1d_Error("BDT_porCathode_LiKCl::ReadModelDescriptions()", 
			 "Electrode model create method failed");
    }
    retn = Electrode_->setInitialConditions(PSCinput_ptr->anode_input_);
    if (retn == -1) {
	throw  m1d_Error("BDT_porCathode_LiKCl::ReadModelDescriptions()", 
			 "setInitialConditions method failed");
    }
    
    delete cfA;
}
//=====================================================================================================================
//  Make list of the equations and variables
/*
 *  We also set the ordering here.
 */
void
BDT_porCathode_LiKCl::SetEquationsVariablesList()
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
    const std::vector<std::string> & namesSp = ionicLiquidIFN_->speciesNames();
    //int nsp = ionicLiquidIFN_->nSpecies();

    if (namesSp[0] != "Li+") {
        exit(-1);
    }
    if (namesSp[1] != "K+") {
        exit(-1);
    }
    if (namesSp[2] != "Cl-") {
        exit(-1);
    }
    VariableNameList.push_back(VarType(MoleFraction_Species, 0, (namesSp[0]).c_str()));
    VariableNameList.push_back(VarType(MoleFraction_Species, 1, (namesSp[1]).c_str()));
    VariableNameList.push_back(VarType(MoleFraction_Species, 2, (namesSp[2]).c_str()));

    // Species conservation is used to solve for mole fractions
    // Equation 1: Species Conservation Li+        Variable 1: Li+ Mole Fraction
    // Equation 2: Species Conservation K+         Variable 2: K+  Mole Fraction
    // Equation 3: Species Conservation Cl-        Variable 3: Cl- Mole Fraction

    EquationNameList.push_back(EqnType(Species_Conservation, 0, (namesSp[0]).c_str()));
    EquationNameList.push_back(EqnType(MoleFraction_Summation, 0));
    IsAlgebraic_NE[2] = 2;
    EquationNameList.push_back(EqnType(ChargeNeutrality_Summation, 0));
    IsAlgebraic_NE[3] = 2;
    eqnIndex++;
    eqnIndex++;
    eqnIndex++;

    //   Current conservation is used to solve for electrostatic potential
    //           Equation 4: Current Conservation - Electrolyte   Variable 4: Volts_Electrolyte
    EquationNameList.push_back(EqnType(Current_Conservation, 0, "Current Conservation"));
    VariableNameList.push_back(VarType(Voltage, 0, "Electrolyte"));
    IsAlgebraic_NE[eqnIndex] = 1;
    IsArithmeticScaled_NE[eqnIndex] = 1;
    eqnIndex++;

    // Current conservation is used to solve for electrostatic potential
    //  Equation 5: Current Conservation - Cathode   Variable 5: Volts_cathode
    EquationNameList.push_back(EqnType(Current_Conservation, 2, "Cathode Current Conservation"));
    VariableNameList.push_back(VarType(Voltage, 2, "CathodeVoltage"));
    IsAlgebraic_NE[eqnIndex] = 1;
    IsArithmeticScaled_NE[eqnIndex] = 1;
    eqnIndex++;

    // Enthalpy conservation is used to solve for the temperature
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
BDT_porCathode_LiKCl::mallocDomain1D()
{
    BulkDomainPtr_ = new porousLiKCl_FeS2Cathode_dom1D(this);
    return BulkDomainPtr_;
}
//=====================================================================================================================
void
BDT_porCathode_LiKCl::DetermineConstitutiveModels()
{
    trans_ = ZZCantera::newTransportMgr("Liquid", ionicLiquidIFN_, 1);
}
//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================
