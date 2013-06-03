/**
 * @file m1d_SDT_FlatAnode.cpp
 */

/*
 *  $Id: m1d_SDT_FlatAnode.cpp 550 2013-03-01 20:53:19Z hkmoffa $
 */

#include "m1d_SDT_FlatAnode.h"
#include "m1d_SurDomain_FlatLiSiAnode.h"

#include "Electrode_input.h"
#include "Electrode_SuccessiveSubstitution.h"

#include "BlockEntry.h"

#include  <string>

using namespace std;
using namespace Cantera;


//=====================================================================================================================
namespace m1d
{
//=====================================================================================================================
SDT_FlatAnode::SDT_FlatAnode(DomainLayout *dl_ptr, int pos) :
  SDT_Mixed(dl_ptr), m_position(pos), ElectrodeA_(0)
{
  ElectrodeA_ = new Cantera::Electrode_SuccessiveSubstitution();

  ELECTRODE_KEY_INPUT *electrodeA_input = new ELECTRODE_KEY_INPUT();

  std::string commandFileA = "anode.inp";
  /**
   * Initialize a block input structure for the command file
   */
  BEInput::BlockEntry *cfA = new BEInput::BlockEntry("command_file");

  /*
   * Go get the problem description from the input file
   */
  electrodeA_input->printLvl_ = 5;
  int retn = electrodeA_input->electrode_input(commandFileA, cfA);

  if (retn == -1) {
    printf("exiting with error\n");
    exit(-1);
  }

  retn = ElectrodeA_->electrode_model_create(electrodeA_input);
  if (retn == -1) {
    printf("exiting with error\n");
    exit(-1);
  }

  retn = ElectrodeA_->setInitialConditions(electrodeA_input);
  if (retn == -1) {
    printf("exiting with error\n");
    exit(-1);
  }


  /*
   *  Add an equation for this surface domain
   *    For the anode we will install a boundary condition on it of a constant voltage
   *    This is the voltage datum for the system.
   */
  EquationNameList.push_back(EqnType(Current_Conservation, 1, "Anode Current Conservation"));
  VariableNameList.push_back(VarType(Voltage, 1, "AnodeVoltage"));

  safeDelete(cfA);
  safeDelete(electrodeA_input);
}
//=====================================================================================================================
SDT_FlatAnode::SDT_FlatAnode(const SDT_FlatAnode &r) :
  SDT_Mixed(r.DL_ptr_), m_position(0)
{
  *this = r;
}
//=====================================================================================================================
SDT_FlatAnode::~SDT_FlatAnode()
{
  safeDelete(ElectrodeA_);
}
//=====================================================================================================================
SDT_FlatAnode &
SDT_FlatAnode::operator=(const SDT_FlatAnode &r)
{
  if (this == &r) {
    return *this;
  }

  SDT_Mixed::operator=(r);
  m_position = r.m_position;


  delete ElectrodeA_;
  ElectrodeA_ = new Cantera::Electrode_SuccessiveSubstitution((const Cantera::Electrode_SuccessiveSubstitution &)*(r.ElectrodeA_));

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
SDT_FlatAnode::SetEquationDescription()
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
   *  Add the dirichlet condition
   */
  EqnType e1 = EquationNameList[0];
  VarType v1 = VariableNameList[0];
  addDirichletCondition(e1, v1, 0.0);
}
//=====================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
SurDomain1D *
SDT_FlatAnode::mallocDomain1D()
{
  SurDomain_FlatLiSiAnode * s1d = new SurDomain_FlatLiSiAnode(*this, 1);
  return s1d;
}

//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================

