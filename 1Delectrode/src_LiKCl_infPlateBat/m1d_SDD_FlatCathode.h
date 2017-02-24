/**
 * @file m1d_SurfDomainTypes.h
 */

#ifndef M1D_SDD_FLATCATHODE_H_
#define M1D_SDD_FLATCATHODE_H_

#include "m1d_SDD_Mixed.h"

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
class Electrode;
}

namespace m1d
{

//! This class specifies that all equations are handled
//! by a simple Dirichlet condition
/*!
 *
 */
class SDD_FlatCathode : public SDD_Mixed
{
public:

  //! Constructor
  /*!
   *   We construct the object but don't actually specify any Dirichlet conditions.
   *   Later we can add dirichlet conditions into the object.
   *      *
   * In the constructor, we have typically been laying out what the unknowns are
   * and what the equations are, that are solved within the domain.
   *
   *
   * @param dl_ptr  Domain Layout object that owns this description.
   */
  SDD_FlatCathode(DomainLayout *dl_ptr, int position);

  //! Destructor
  virtual
  ~SDD_FlatCathode();

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  SDD_FlatCathode(const SDD_FlatCathode &r);

  //! Assignment operator
  /*!
   * @param r    Object to be copied
   * @return     Returns a changeable reference to the current object
   */
  SDD_FlatCathode &
  operator=(const SDD_FlatCathode &r);


  //! Set the equation and variables list
  /*!
   *  This routine is responsible for setting the variables:
   *    - VariableNameList
   *    - EquationNameList
   */
  virtual void SetEquationsVariablesList();

  //! Set the equation description
  /*!
   *  This routine is responsible for setting the variables:
   *    - NumEquationsPerNode
   *    - EquationIndexStart_EqName
   */
  virtual void SetEquationDescription();

  //! Malloc and Return the object that will calculate the residual efficiently
  /*!
   *
   * @return  Returns a pointer to the object that will calculate the residual
   *          efficiently
   */
  virtual SurDomain1D *
  mallocDomain1D();

  //! top or bottom of the domain
  /*!
   *   0 - top, right
   *   1 - bottom, left
   */
  int m_position;

  //! Pointer to the electrode object
  /*!
   * We own the electrode object.
   */
  ZZCantera::Electrode *ElectrodeC_;

  //! Type of the boundary condition specified on the cathode
  /*!
   *   0 specify the voltage
   *   1 specify the current
   */
  int voltageVarBCType_;

  //! Specified current in the cathode
  /*!
   *  This is actually the current from the cathode into the electrolyte.
   *  Therefore, during a normal discharge operation of the battery, this will be a
   *  negative quantity.
   *
   *   Note, this is only relevant when voltageVarBCType_ = 1
   */
  double icurrCathodeSpecified_;

  //! Make the SurDomain1D class a friend so that it can access all of the stuff in this class
  friend class SurDomain_FlatFeS2Cathode;
};

}

#endif /*  */
