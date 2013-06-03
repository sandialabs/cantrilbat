/**
 * @file m1d_SurfDomainTypes.h
 */
/*
 * $Id: m1d_SDT_Mixed_Cu2S.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef M1D_SDT_MIXED_CU2S_TOP_H_
#define M1D_SDT_MIXED_CU2S_TOP_H_

#include "m1d_SurfDomainTypes.h"

namespace m1d
{

//! This class specifies that all equations are handled
//! by a simple Dirichlet condition
/*!
 *
 */
class SDT_Mixed_Cu2S : public SDT_Mixed {
public:

  //! Constructor
  /*!
   *   We construct the object but don't actually specify any Dirichlet conditions.
   *   Later we can add dirichlet conditions into the object.
   *
   * @param dl_ptr  Domain Layout object that owns this description.
   */
  SDT_Mixed_Cu2S(DomainLayout *dl_ptr, int position);

  //! Destructor
  virtual
  ~SDT_Mixed_Cu2S();

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  SDT_Mixed_Cu2S(const SDT_Mixed_Cu2S &r);

  //! Assignment operator
  /*!
   * @param r    Object to be copied
   * @return     Returns a changeable reference to the current object
   */
  SDT_Mixed_Cu2S &
  operator=(const SDT_Mixed_Cu2S &r);


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
   *   0 - top
   *   1 - bottom
   */
  int m_position;
};


}

#endif /*  */
