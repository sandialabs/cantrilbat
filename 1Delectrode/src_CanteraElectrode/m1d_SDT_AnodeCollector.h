/**
 * @file m1d_SurfDomainTypes.h
 */
/*
 * $Id: m1d_SDT_AnodeCollector.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef M1D_SDT_ANODECOLLECTOR_H_
#define M1D_SDT_ANODECOLLECTOR_H_

#include "m1d_SurfDomainTypes.h"

namespace Cantera
{
class ELECTRODE_MODEL;
}
namespace m1d
{

//! This class specifies that all equations are handled
//! by a simple Dirichlet condition
/*!
 *
 */
class SDT_AnodeCollector : public SDT_Mixed
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
   * @param dl_ptr  Domain Layout object that owns this description.
   */
  SDT_AnodeCollector(DomainLayout *dl_ptr, int position, const char *domainName = "");

  //! Destructor
  virtual ~SDT_AnodeCollector();

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  SDT_AnodeCollector(const SDT_AnodeCollector &r);

  //! Assignment operator
  /*!
   * @param r    Object to be copied
   * @return     Returns a changeable reference to the current object
   */
  SDT_AnodeCollector &
  operator=(const SDT_AnodeCollector &r);

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

  //! top or bottom of the domain
  /*!
   *   0 - top, right
   *   1 - bottom, left
   */
  int m_position;

  //! Make the SurDomain1D class a friend so that it can access all of the stuff in this class
  friend class SurDomain_AnodeCollector;
};

}

#endif /*  */
