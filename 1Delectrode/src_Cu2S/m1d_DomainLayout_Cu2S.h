/**
 * @file m1d_DomainLayout.h  This is a top level object that describes
 *         the domain layout of the problem.
 */

/*
 *  $Id: m1d_DomainLayout_Cu2S.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
#ifndef M1D_DOMAINLAYOUT_CU2S_H
#define M1D_DOMAINLAYOUT_CU2S_H

#include "m1d_DomainLayout.h"

namespace m1d
{



//! DomainLayout is a light class that describes the overall
//! disposition of the domains in the problem
/*!
 * 
 */
class DomainLayout_Cu2S : public DomainLayout
{
public:
  //! Constructor
  DomainLayout_Cu2S();

  //! Constructor - hard coded problem
  DomainLayout_Cu2S(int probNum, ProblemStatement *psInput_ptr);

  //! Copy constructor
  /*!
   * @param r  Object to be copied
   */
  DomainLayout_Cu2S(const DomainLayout_Cu2S &r);

  //! Destructor
  virtual
  ~DomainLayout_Cu2S();

  //! Assignment operator
  /*!
   * @param r  Object to be copied.
   * @return   returns a changeable reference to the current object
   */
  DomainLayout_Cu2S &
  operator=(const DomainLayout_Cu2S &r);

  //! Allocate the domain structure
  virtual void
  malloc_domains();


  //==================================================================================
  //! Hard coded problem number
  int ProbNum_;

};
//=====================================================================================
} // end namespace m1d
//=====================================================================================
#endif
