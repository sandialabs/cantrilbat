/**
 * @file m1d_SurfDomainTypes.cpp
 */

/*
 *  $Id: m1d_SDT_CathodeCollector.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "m1d_SDT_CathodeCollector.h"
#include "m1d_SurDomain_CathodeCollector.h"

#include "m1d_ProblemStatementCell.h"

#include "Electrode_input.h"
#include "Electrode.h"

#include "BlockEntry.h"

#include  <string>

extern m1d::ProblemStatementCell PSinput;

//=====================================================================================================================
namespace m1d
{
//=====================================================================================================================
SDT_CathodeCollector::SDT_CathodeCollector(DomainLayout *dl_ptr, int pos, const char *domainName) :
  SDT_Mixed(dl_ptr,domainName), m_position(pos), voltageVarBCType_(0), icurrCathodeSpecified_(0.0), voltageCathodeSpecified_(1.9)
{

  voltageVarBCType_ = PSinput.cathodeBCType_;
  icurrCathodeSpecified_ = PSinput.icurrDischargeSpecified_;

  /*
   *  Add an equation for this surface domain
   *    For the cathode we will install a boundary condition on it of either
   *    a constant voltage or a constant current.
   */

}
//=====================================================================================================================
SDT_CathodeCollector::SDT_CathodeCollector(const SDT_CathodeCollector &r) :
  SDT_Mixed(r.DL_ptr_), m_position(0), voltageVarBCType_(0), icurrCathodeSpecified_(0.0), voltageCathodeSpecified_(1.9)
{
  *this = r;
}
//=====================================================================================================================
SDT_CathodeCollector::~SDT_CathodeCollector()
{
}
//=====================================================================================================================
SDT_CathodeCollector &
SDT_CathodeCollector::operator=(const SDT_CathodeCollector &r)
{
  if (this == &r) {
    return *this;
  }

  SDT_Mixed::operator=(r);

  m_position = r.m_position;

  voltageVarBCType_ = r.voltageVarBCType_;
  icurrCathodeSpecified_ = r.icurrCathodeSpecified_;
  voltageCathodeSpecified_ = r.voltageCathodeSpecified_;

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
SDT_CathodeCollector::SetEquationDescription()
{
  /*
   * Set the policy for connecting bulk domains
   * This really isn't set yet.
   */
  setRLMapping(0);
  /*
   * Fill in the rest of the information
   */
  SurfDomainDescription::SetEquationDescription();

  /*
   *  If we are just fixing the voltage at the cathode, we can set the plain Dirichlet condition here.
   *  If we are setting the current, we will add in a residual equation by hand in residEval().
   */
  voltageCathodeSpecified_ = PSinput.CathodeVoltageSpecified_;
  double (*timeDepFunction)(double) = PSinput.TimeDepFunction_;
  BoundaryCondition * BC_timeDep = PSinput.BC_TimeDep_;
  EqnType e1(Current_Conservation, 2, "Cathode Current Conservation");
  VarType v1(Voltage, 2, "CathodeVoltage");
  switch (voltageVarBCType_)
    {
    case 0:
      addDirichletCondition(e1, v1, voltageCathodeSpecified_ );
      break;
    case 1:
      addFluxCondition(e1, v1, icurrCathodeSpecified_ );
      break;
    case 2:
      addDirichletCondition(e1, v1, voltageCathodeSpecified_, timeDepFunction );
      break;
    case 3:
      addFluxCondition(e1, v1, icurrCathodeSpecified_, timeDepFunction );
      break;
    case 4:
    case 6:
    case 8:
      addDirichletCondition(e1, v1, voltageVarBCType_, BC_timeDep );
      break;
    case 5:
    case 7:
    case 9:
      addFluxCondition(e1, v1, voltageVarBCType_, BC_timeDep);
      break;
    default:
      throw m1d_Error("SDT_CathodeCollector::SetEquationDescription", 
		      "voltageVarBCType not 0, 1, 2 for Dirichlet, Neumann, or Diriclet with sin oscillation" );
    }

  /*
   *  All of the other boundary conditions default to zero flux at the interface
   *      This includes:
   *
   *           flux of Li+
   *           flux of K+
   *           flux of Cl-
   *           flux of current
   *
   *  Because of the staggered grid, a dirichlet condition must be set on the axial velocity. Or else the
   *  last axial velocity unknown is not represented in the solution vector, leading to a singular matrix.
   *
   */
  EqnType e2(Continuity, 0, "Continuity: Bulk Velocity");
  VarType v2(Velocity_Axial, 0, "Axial_Velocity");
  addDirichletCondition(e2, v2, 0.0);
}
//=====================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
SurDomain1D *
SDT_CathodeCollector::mallocDomain1D()
{
  SurDomain_CathodeCollector * s1d = new SurDomain_CathodeCollector(*this, 1);
  return s1d;
}

//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================

