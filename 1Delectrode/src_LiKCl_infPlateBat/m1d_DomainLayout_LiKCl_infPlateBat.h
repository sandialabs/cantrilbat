/**
 * @file m1d_DomainLayout_LiKCl_infPlateBat.h  
 *         This is a top level object that describes
 *         the domain layout of the problem.
 */

/*
 *  $Id: m1d_DomainLayout_LiKCl_infPlateBat.h 506 2013-01-07 22:43:59Z hkmoffa $
 */
#ifndef M1D_DOMAINLAYOUT_LIKCL_INFPLATEBAT_H
#define M1D_DOMAINLAYOUT_LIKCL_INFPLATEBAT_H

//#include "zuzax/transport.h"      // transport properties
//#include "zuzax/thermo.h"      // transport properties


#include "m1d_DomainLayout.h"

namespace m1d
{
  class ProblemStatementCell;
  //! DomainLayout is a light class that describes the overall
  //! disposition of the domains in the problem
  /*!
   *   Most of the information is located in the base class. The base class
   *   maintains a listing of the ordering of the bulk and surface domains.
   *   The class is responsible for mallocing the domain des
   */
  class DomainLayout_LiKCl_infPlateBat : public DomainLayout
  {
  public:
    //! Constructor
    DomainLayout_LiKCl_infPlateBat(ProblemStatement *psInput_ptr);

    //! Constructor - hard coded problem
    DomainLayout_LiKCl_infPlateBat(int probNum, ProblemStatement *psInput_ptr);


    //! Copy constructor
    /*!
     * @param r  Object to be copied
     */
    DomainLayout_LiKCl_infPlateBat(const DomainLayout_LiKCl_infPlateBat&r);

    //! Destructor
    virtual
    ~DomainLayout_LiKCl_infPlateBat();

    //! Assignment operator
    /*!
     * @param r  Object to be copied.
     * @return   returns a changeable reference to the current object
     */
    DomainLayout_LiKCl_infPlateBat &
    operator=(const DomainLayout_LiKCl_infPlateBat &r);

    //! Allocate the domain structure
    virtual void
    malloc_domains();


    //==================================================================================
    //! Hard coded problem number
    int ProbNum_;

    ProblemStatementCell * pscInput_ptr_;


  protected:
  
  };
  //=====================================================================================
} // end namespace m1d
//=====================================================================================
#endif
