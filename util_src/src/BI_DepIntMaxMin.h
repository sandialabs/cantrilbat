/**
 * @file BI_DepIntMaxMin.h
 *   Declarations for the BI_DepIntMaxMin class
 */
/*
 * $Author: hkmoffa $
 * $Revision: 5 $
 * $Date: 2012-02-23 14:34:18 -0700 (Thu, 23 Feb 2012) $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */
#ifndef BI_DEPINTMAXMIN_H
#define BI_DEPINTMAXMIN_H


#include "BI_Dependency.h"
#include "BaseEntry.h"

namespace BEInput {

  //! Conditional Dependency check that queries the target BaseEntry expecting
  //! an int response. The condition is that the int response must be
  //! between a range of values
  /*!
   *  This is a derived dependency.  The key here is that this dependency
   *  requires that the target BaseEntry be able to return an int.
   *
   * BIDT_type
   *    This is the exact definition of what it means for the
   *    dependency to be satisfied. There are two choices
   *    that make sense for this class:
   *        -  BIDT_ONEINT
   *        -  BIDT_INTMAXMIN
   *
   *    The ResultType for the dependency can be almost any type.
   */
  class BI_DepIntMaxMin : public BI_Dependency {
  public:

    //! constructor
    /*!
     * @param be   This is the BaseEntry object that the dependent BaseEntry
     *             is dependent on or is a prerequisite for.
     *
     * @param maxTargetI  max int value for a positive response
     *
     * @param minTargetI  min int value for a positive response
     *
     * @param BIDT_type  This is the exact definition of what it means for the
     *                   dependency to be satisfied. There are two choices
     *                   that make sense here:
     *                   -   BIDT_ONEINT
     *                   -   BIDT_INTMAXMIN
     *
     * @param BIDRT_type   This is the result -> Basically what happens if the
     *                     dependency is met or not met. The default is to 
     *                     throw an error condition if the dependency is not
     *                     met and to do nothing if the dependency is met.
     */
    BI_DepIntMaxMin(BaseEntry *be, BIDT_TYPE BIDT_type,
		    int maxTargetI , int minTargetI, 
		    BIDRT_TYPE BIDRT_type = BIDRT_PT_ERROR);

    //! Copy Constuctor
    /*!
     * @param b object to be copied
     */
    BI_DepIntMaxMin(const BI_DepIntMaxMin& b);

    //! Assignment operator
    /*!
     * @param b object to be copied
     */
    BI_DepIntMaxMin& operator=(const BI_DepIntMaxMin&b);

    //! destructor
    virtual ~BI_DepIntMaxMin();

    //! duplicator
    virtual BI_Dependency* duplicateMyself() const;

    //! Check to see if the dependency is satisfied
    /*!
     *  This requires that the target dependency return an int.
     */
    virtual bool checkDependencySatisfied() const;

    //! Check to see if the dependency is satisfied
    /*!
     *  This returns the int as well.
     *
     * @param returnInt Value of the returned int
     *
     * @return
     *  returns a boolean indicating whether the dependency has
     *  been satisified.
     */
    virtual bool checkDepOneInt(int &returnInt) const;

  protected:

    //! max value to trigger the dependency check
    int m_maxValue;

    //! min value to trigger the dependency check
    int m_minValue;

    //! holds the int variable returned from the
    //! dependency query.
    mutable int m_RetnOneInt; 
  };

}

#endif
