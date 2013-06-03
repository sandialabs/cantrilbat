/**
 * @file m1d_DomainLayout_LiIon_PorousBat.h  
 *         This is a top level object that describes
 *         the domain layout of the problem.
 */

/*
 *  $Id: m1d_DomainLayout_LiIon_PorousBat.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
#ifndef M1D_DOMAINLAYOUT_LIION_POROUSBAT_H
#define M1D_DOMAINLAYOUT_LIION_POROUSBAT_H

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
  class DomainLayout_LiIon_PorousBat : public DomainLayout
  {
  public:
    //! Constructor
    DomainLayout_LiIon_PorousBat(ProblemStatement *psInput_ptr);

    //! Constructor - hard coded problem
    DomainLayout_LiIon_PorousBat(int probNum, ProblemStatement *psInput_ptr);


    //! Copy constructor
    /*!
     * @param r  Object to be copied
     */
    DomainLayout_LiIon_PorousBat(const DomainLayout_LiIon_PorousBat &r);

    //! Destructor
    virtual
    ~DomainLayout_LiIon_PorousBat();

    //! Assignment operator
    /*!
     * @param r  Object to be copied.
     * @return   returns a changeable reference to the current object
     */
    DomainLayout_LiIon_PorousBat &
    operator=(const DomainLayout_LiIon_PorousBat &r);

    //! Allocate the domain structure
    virtual void
    malloc_domains();

    // ---------------------------------------------------------------------------------------------
    //                                 DATA
    // ---------------------------------------------------------------------------------------------
   
    //! Hard coded problem number
    int ProbNum_;

    ProblemStatementCell *pscInput_ptr_;

  };
  //=====================================================================================
} // end namespace m1d
//=====================================================================================
#endif
