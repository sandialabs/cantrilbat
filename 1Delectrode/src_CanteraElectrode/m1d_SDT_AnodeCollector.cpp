/**
 * @file m1d_SDT_AnodeCollector.cpp
 */

/*
 *  $Id: m1d_SDT_AnodeCollector.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "m1d_SDT_AnodeCollector.h"
#include "m1d_SurDomain_AnodeCollector.h"

#include "Electrode_input.h"
#include "Electrode.h"

#include "BlockEntry.h"

#include  <string>

using namespace std;
using namespace Cantera;

//=====================================================================================================================
namespace m1d
{
//=====================================================================================================================
  SDT_AnodeCollector::SDT_AnodeCollector(DomainLayout *dl_ptr, int pos, const char *domainName) :
    SDT_Mixed(dl_ptr,domainName), m_position(pos)
{
  /*
   * For this implementation, we are not adding a separate equation for this surface domain.
   * The current equation for the electrode phase from the bulk domain is used as the current in the anode
   * collector.
   *
   * This object will set an anode voltage of zero on that equation.
   */
}
//=====================================================================================================================
SDT_AnodeCollector::SDT_AnodeCollector(const SDT_AnodeCollector &r) :
  SDT_Mixed(r.DL_ptr_), m_position(0)
{
  *this = r;
}
//=====================================================================================================================
SDT_AnodeCollector::~SDT_AnodeCollector()
{
}
//=====================================================================================================================
SDT_AnodeCollector &
SDT_AnodeCollector::operator=(const SDT_AnodeCollector &r)
{
  if (this == &r) {
    return *this;
  }

  SDT_Mixed::operator=(r);
  m_position = r.m_position;

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
SDT_AnodeCollector::SetEquationDescription()
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
    *  Add the Dirichlet condition onto the bulk domain equation
    */
   EqnType e1(Current_Conservation, 1, "Anode Current Conservation");
   VarType v1(Voltage, 1, "AnodeVoltage");
   addDirichletCondition(e1, v1, 0.0);
   /*
    *  All of the other boundary conditions default to zero flux at the interface
    *      This includes:
    *
    *           axial velocity
    *           flux of Li+
    *           flux of K+
    *           flux of Cl-
    *           flux of current
    */

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
   VarType v2(Velocity_Axial, 0, "Axial velocity");
   addFluxCondition(e2, v2, 0.0);
}
//=====================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
SurDomain1D *
SDT_AnodeCollector::mallocDomain1D()
{
  SurDomain_AnodeCollector * s1d = new SurDomain_AnodeCollector(*this, 1);
  return s1d;
}

//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================

