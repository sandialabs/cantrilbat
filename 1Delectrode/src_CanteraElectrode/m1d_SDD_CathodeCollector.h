/**
 * @file m1d_CathodeCollector.h
 */

#ifndef M1D_SDD_CATHODECOLLECTOR_H_
#define M1D_SDD_CATHODECOLLECTOR_H_

#include "m1d_SDD_Mixed.h"

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
class ELECTRODE_MODEL;
}
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! This class specifies that all equations are handled
//! by a simple Dirichlet condition or simple flux conditions of the third type
/*!
 *
 */
class SDD_CathodeCollector : public SDD_Mixed
{
public:

  //! Constructor
  /*!
   *   We construct the object but don't actually specify any Dirichlet conditions.
   *   Later we can add dirichlet conditions into the object.
   *
   * In the constructor, we have typically been laying out what the unknowns are
   * and what the equations are, that are solved within the domain.
   *
   *
   * @param dl_ptr  Domain Layout object that owns this description.
   */
  SDD_CathodeCollector(DomainLayout *dl_ptr, int position, const char *domainName = "");

  //! Destructor
  virtual
  ~SDD_CathodeCollector();

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  SDD_CathodeCollector(const SDD_CathodeCollector &r);

  //! Assignment operator
  /*!
   * @param r    Object to be copied
   * @return     Returns a changeable reference to the current object
   */
  SDD_CathodeCollector &
  operator=(const SDD_CathodeCollector &r);

  //! Set the equation description
  /*!
   *  This routine is responsible for setting the variables:
   *    - NumEquationsPerNode
   *    - VariableNameList
   *    - EquationNameList
   *    - EquationIndexStart_EqName
   */
  virtual void
  SetEquationDescription();

  //! Malloc and Return the object that will calculate the residual efficiently
  /*!
   *
   * @return  Returns a pointer to the object that will calculate the residual
   *          efficiently
   */
  virtual SurDomain1D *
  mallocDomain1D();

  // --------------------------------------------------------------------------------------------------------------
  //                                   DATA
  // --------------------------------------------------------------------------------------------------------------

  //! top or bottom of the domain
  /*!
   *   0 - top, right
   *   1 - bottom, left
   */
  int m_position;

  //! Type of the boundary condition specified on the cathode
  /*!
   *   0 specify the voltage
   *   1 specify the current
   *
   *  10 Set a Robin boundary condition ---   current = R_anodeCC (v_acc - vanode)
   */
  int voltageVarBCType_;

  //! Specified current in the cathode
  /*!
   *  This is actually the current from the cathode into the electrolyte.
   *  Therefore, during a normal discharge operation of the battery, this will be a
   *  negative quantity.
   *
   *   Note, this is only relevant when voltageVarBCType_ = 1, 3, 5, 7, 9
   */
  double icurrCathodeSpecified_;

  //! Specified voltage at the cathode and sometimes at the cathode current collector
  double voltageCathodeSpecified_;

  //!  Thickness of the cathode current collector
  double cathodeCCThickness_;

  //! Extra resistance attached to the entire battery that is still included in heat buildup
  /*!
   *    Note this is not cross section based. This is for the entire battery. Therefore, the
   *    cross section must be multiplied in.
   *
   *    Units: ohms
   */
  double extraResistanceCathode_;

  //! Load resistance attached to the entire battery 
  /*!
   *    Note this is not cross section based. This is for the entire battery. Therefore, the
   *    cross section must be multiplied in
   *
   *    units ohms
   */
  double ResistanceLoad_;
  double VoltageLoad_;


  //!  Type of the temperature boundary condition
  /*!
   *  There are three possibilities
   *     0 set the cathode temperature to a constant
   *     1 set the flux to a constant
   *    10 Set a Robin boundary condition - heatflux = h ( T - T_cath)
   */
  int cathodeTempBCType_;
  double cathodeTempCollector_;
  double cathodeHeatTranCoeff_;


  //! Make the SurDomain1D class a friend so that it can access all of the stuff in this class
  friend class SurDomain_CathodeCollector;
};

}

#endif /*  */
