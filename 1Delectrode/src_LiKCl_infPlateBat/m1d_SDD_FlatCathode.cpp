/**
 * @file m1d_SurfDomainTypes.cpp
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "m1d_SDD_FlatCathode.h"
#include "m1d_SurDomain_FlatFeS2Cathode.h"
#include "m1d_defs.h"

#include "m1d_ProblemStatementCell.h"

#include "cantera/base/ctexceptions.h"

#include "Electrode.h"
#include "Electrode_Factory.h"
#include <cantera/transport.h>      // transport properties
#include <cantera/thermo.h>      // transport properties
#include <cantera/thermo/IonsFromNeutralVPSSTP.h>  // ion properties



#include  <string>

extern m1d::ProblemStatementCell  PSinput;

#ifdef useZuzaxNamespace
#ifndef ZZCantera
#define ZZCantera Zuzax
#endif
#else
#ifndef ZZCantera
#define ZZCantera Cantera 
#endif
#endif
using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif


//=====================================================================================================================
namespace m1d
{
//=====================================================================================================================
SDD_FlatCathode::SDD_FlatCathode(DomainLayout *dl_ptr, int pos) :
    SDD_Mixed(dl_ptr), 
    m_position(pos),
    ElectrodeC_(nullptr), 
    voltageVarBCType_(0), 
    icurrCathodeSpecified_(0.0)
{
/*
  ElectrodeC_ = new ZZCantera::Electrode_SuccessiveSubstitution();

  ZZCantera::ELECTRODE_KEY_INPUT *electrodeC_input = new ZZCantera::ELECTRODE_KEY_INPUT();

  std::string commandFileC = "cathode.inp";
   // Initialize a block input structure for the command file
  BEInput::BlockEntry *cfC = new BEInput::BlockEntry("command_file");

   // Go get the problem description from the input file
  electrodeC_input->printLvl_ = 5;
  int retn = electrodeC_input->electrode_input(commandFileC, cfC);

  if (retn == -1) {
    printf("exiting with error\n");
    exit(-1);
  }

  retn = ElectrodeC_->electrode_model_create(electrodeC_input);
  if (retn == -1) {
    printf("exiting with error\n");
    exit(-1);
  }

  retn = ElectrodeC_->setInitialConditions(electrodeC_input);
  if (retn == -1) {
    printf("exiting with error\n");
    exit(-1);
  }

  voltageVarBCType_ = PSinput.cathodeBCType_;
  icurrCathodeSpecified_ = - PSinput.icurrDischargeSpecified_ * 1.0E4;
*/

    size_t iph = (PSinput.PhaseList_)->globalPhaseIndex(PSinput.electrolytePhase_);
    if (iph == npos) {
        //throw ZuzaxError("SDD_FlatCathode::SDD_FlatCathode", "Can't find the phase in the phase list: " + PSinput.electrolytePhase_);
        throw Zuzax::ZuzaxError("SDD_FlatCathode::SDD_FlatCathode", "Can't find the phase in the phase list: " );
    }
    thermo_t_double* ionicLiquid_ = & (PSinput.PhaseList_)->thermo(iph);
    ionicLiquidIFN_ = dynamic_cast<ZZCantera::IonsFromNeutralVPSSTP *>( ionicLiquid_->duplMyselfAsThermoPhase() );
    ionicLiquid_ = (ThermoPhase*) ionicLiquidIFN_;

    ELECTRODE_KEY_INPUT *ci = PSinput.cathode_input_;
    ElectrodeC_ = newElectrodeObject(ci->electrodeModelName);
    if (!ElectrodeC_) {
        throw  m1d_Error("SDD_FlatCathode::SDD_FlatCathode", "Electrode factory method failed");
    }
    ELECTRODE_KEY_INPUT *ci_new = newElectrodeKeyInputObject(ci->electrodeModelName);
    std::string commandFile = ci->commandFile_;
    BEInput::BlockEntry *cfC = new BEInput::BlockEntry("command_file");

    
    //  Parse the complete child input file
    int retn = ci_new->electrode_input_child(commandFile, cfC);
    if (retn == -1) {
        throw  m1d_Error("SDD_FlatCathode::SDD_FlatCathode", "Electrode input child method failed");
    }
    /*
     * Switch the pointers around so that the child input file is returned.
     * Delete the original pointer.
     */
    delete ci;
    PSinput.cathode_input_ = ci_new;

    retn = ElectrodeC_->electrode_model_create(PSinput.cathode_input_);
    if (retn == -1) {
        throw  m1d_Error("SDD_FlatCathode::SDD_FlatCathode", "Electrode model create method failed");
    }
    retn = ElectrodeC_->setInitialConditions(PSinput.cathode_input_);
    if (retn == -1) {
        throw  m1d_Error("SDD_FlatCathode::SDD_FlatCathode", "setInitialConditions method failed");
    }



  voltageVarBCType_ = PSinput.cathodeBCType_;
  icurrCathodeSpecified_ = - PSinput.icurrDischargeSpecified_ * 1.0E4;

  safeDelete(cfC);
}
//=====================================================================================================================
SDD_FlatCathode::SDD_FlatCathode(const SDD_FlatCathode &r) :
    SDD_Mixed(r.DL_ptr_), 
    m_position(0), ElectrodeC_(0), voltageVarBCType_(0), icurrCathodeSpecified_(0.0)
{
  *this = r;
}
//=====================================================================================================================
SDD_FlatCathode::~SDD_FlatCathode()
{
  safeDelete(ElectrodeC_);
}
//=====================================================================================================================
SDD_FlatCathode &
SDD_FlatCathode::operator=(const SDD_FlatCathode &r)
{
  if (this == &r) {
    return *this;
  }

  SDD_Mixed::operator=(r);

  m_position = r.m_position;

  delete ElectrodeC_;
  ElectrodeC_ = (r.ElectrodeC_)->duplMyselfAsElectrode();

  voltageVarBCType_ = r.voltageVarBCType_;
  icurrCathodeSpecified_ = r.icurrCathodeSpecified_;

  return *this;
}
//======================================================================================================================
// Determine the list of Equations and Variables
/*
 *  This routine is responsible for setting the variables:
 *    - VariableNameList
 *    - EquationNameList
 */
void
SDD_FlatCathode::SetEquationsVariablesList() {

    EquationNameList.clear();
    VariableNameList.clear();  
    /*
     *  Add an equation for this surface domain
     *    For the cathode we will install a boundary condition on it of either
     *    a constant voltage or a constant current.
     */
    if (voltageVarBCType_ == 0) {
	EquationNameList.push_back(EqnType(Voltage_Specification, 2, "Cathode Voltage Specification"));
    } else {
	EquationNameList.push_back(EqnType(Current_Specification, 2, "Cathode Current Conservation"));
    }
    VariableNameList.push_back(VarType(Voltage, 2, "CathodeVoltage"));
}
//=====================================================================================================================
// Set the equation description
/*
 *  This routine is responsible for setting the variables:
 *    - NumEquationsPerNode
 *    - VariableNameList
 *    - EquationNameList
 *    - EquationIndexStart_EqName
 */
void
SDD_FlatCathode::SetEquationDescription()
{

  /*
   * Set the policy for connecting bulk domains
   * This really isn't set yet.
   */
  setRLMapping(1);
  /*
   * Fill in the rest of the information
   */
  SurfDomainDescription::SetEquationDescription();

  /*
   *  If we are just fixing the voltage at the cathode, we can set the plain Dirichlet condition here.
   *  If we are setting the current, we will add in a residual equation by hand in residEval().
   */
  if (voltageVarBCType_ == 0) {
    EqnType e1 = EquationNameList[0];
    VarType v1 = VariableNameList[0];
    addDirichletCondition(e1, v1, PSinput.CathodeVoltageSpecified_);
  }
}
//=====================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
SurDomain1D *
SDD_FlatCathode::mallocDomain1D()
{
  SurDomain_FlatFeS2Cathode * s1d = new SurDomain_FlatFeS2Cathode(*this, 1);
  return s1d;
}

//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================

