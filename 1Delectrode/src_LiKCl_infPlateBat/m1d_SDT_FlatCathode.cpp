/**
 * @file m1d_SurfDomainTypes.cpp
 */

/*
 *  $Id: m1d_SDT_FlatCathode.cpp 550 2013-03-01 20:53:19Z hkmoffa $
 */

#include "m1d_SDT_FlatCathode.h"
#include "m1d_SurDomain_FlatFeS2Cathode.h"

#include "m1d_ProblemStatementCell.h"

#include "Electrode_input.h"
#include "Electrode.h"
#include "Electrode_SuccessiveSubstitution.h"

#include "BlockEntry.h"

#include  <string>

extern m1d::ProblemStatementCell  PSinput;

//=====================================================================================================================
namespace m1d
{
//=====================================================================================================================
SDT_FlatCathode::SDT_FlatCathode(DomainLayout *dl_ptr, int pos) :
  SDT_Mixed(dl_ptr), m_position(pos), ElectrodeC_(0), voltageVarBCType_(0), icurrCathodeSpecified_(0.0)
{
  ElectrodeC_ = new Cantera::Electrode_SuccessiveSubstitution();

  ELECTRODE_KEY_INPUT *electrodeC_input = new ELECTRODE_KEY_INPUT();

  std::string commandFileC = "cathode.inp";
  /**
   * Initialize a block input structure for the command file
   */
  BEInput::BlockEntry *cfC = new BEInput::BlockEntry("command_file");

  /*
   * Go get the problem description from the input file
   */
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

  safeDelete(cfC);
  safeDelete(electrodeC_input);
}
//=====================================================================================================================
SDT_FlatCathode::SDT_FlatCathode(const SDT_FlatCathode &r) :
  SDT_Mixed(r.DL_ptr_), m_position(0), ElectrodeC_(0), voltageVarBCType_(0), icurrCathodeSpecified_(0.0)
{
  *this = r;
}
//=====================================================================================================================
SDT_FlatCathode::~SDT_FlatCathode()
{
  safeDelete(ElectrodeC_);
}
//=====================================================================================================================
SDT_FlatCathode &
SDT_FlatCathode::operator=(const SDT_FlatCathode &r)
{
  if (this == &r) {
    return *this;
  }

  SDT_Mixed::operator=(r);

  m_position = r.m_position;

  delete ElectrodeC_;
  ElectrodeC_ = new Cantera::Electrode_SuccessiveSubstitution((const Cantera::Electrode_SuccessiveSubstitution&)*r.ElectrodeC_);

  voltageVarBCType_ = r.voltageVarBCType_;
  icurrCathodeSpecified_ = r.icurrCathodeSpecified_;

  return *this;
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
SDT_FlatCathode::SetEquationDescription()
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
SDT_FlatCathode::mallocDomain1D()
{
  SurDomain_FlatFeS2Cathode * s1d = new SurDomain_FlatFeS2Cathode(*this, 1);
  return s1d;
}

//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================

