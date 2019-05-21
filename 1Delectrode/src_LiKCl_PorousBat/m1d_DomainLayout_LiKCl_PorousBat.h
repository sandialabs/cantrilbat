/**
 * @file m1d_DomainLayout_LiKCl_PorousBat.h  
 *         This is a top level object that describes
 *         the domain layout of the problem.
 */

/*
 *  $Id: m1d_DomainLayout_LiKCl_PorousBat.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
#ifndef M1D_DOMAINLAYOUT_LIKCL_POROUSBAT_H
#define M1D_DOMAINLAYOUT_LIKCL_POROUSBAT_H

//#include "zuzax/Cantera.h"
//#include "zuzax/transport.h"      // transport properties
//#include "zuzax/thermo.h"      // transport properties


#include "m1d_DomainLayout.h"

namespace m1d
{
  class  ProblemStatementCell;
  //! DomainLayout is a light class that describes the overall
  //! disposition of the domains in the problem
  /*!
   *   Most of the information is located in the base class. The base class
   *   maintains a listing of the ordering of the bulk and surface domains.
   *   The class is responsible for mallocing the domain des
   */
  class DomainLayout_LiKCl_infPorousBat : public DomainLayout
  {
  public:
    //! Constructor
    DomainLayout_LiKCl_infPorousBat(ProblemStatement *psInput_ptr);

    //! Constructor - hard coded problem
    DomainLayout_LiKCl_infPorousBat(int probNum, ProblemStatement *psInput_ptr);


    //! Copy constructor
    /*!
     * @param r  Object to be copied
     */
    DomainLayout_LiKCl_infPorousBat(const DomainLayout_LiKCl_infPorousBat&r);

    //! Destructor
    virtual
    ~DomainLayout_LiKCl_infPorousBat();

    //! Assignment operator
    /*!
     * @param r  Object to be copied.
     * @return   returns a changeable reference to the current object
     */
    DomainLayout_LiKCl_infPorousBat &
    operator=(const DomainLayout_LiKCl_infPorousBat &r);

    //! Allocate the domain structure
    virtual void
    malloc_domains();


    //==================================================================================
    //! Hard coded problem number
    int ProbNum_;

    ProblemStatementCell *pscInput_ptr_;


  protected:
  
  };
  //=====================================================================================
} // end namespace m1d
//=====================================================================================
#endif
