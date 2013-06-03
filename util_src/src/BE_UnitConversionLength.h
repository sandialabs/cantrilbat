/**
 * @file BE_UnitConversion.ch
 *
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


#ifndef BE_UNITCONVERSIONLENGTH_H
#define BE_UNITCONVERSIONLENGTH_H

#include "BE_UnitConversion.h"

namespace BEInput {

  class BE_UnitConversionLength : public BE_UnitConversion {
  public:
   /*
     * Constructor is defined to be private so that we can't make multiple
     * copies of this.
     */
    BE_UnitConversionLength();
 
   //! Constructor is defined to be protected so that we can't make multiple
    //! copies of this, yet we can inherit from this object.
    /*!
     * In the constructor, we define all of the conversion amounts.
     */
    BE_UnitConversionLength(const BE_UnitConversionLength &right);

    //! Assignment operator
    /*!
     *  @param right Object to be copied
     */
    BE_UnitConversionLength& operator=(const BE_UnitConversionLength& right);

    //! duplicator
    virtual BE_UnitConversion * duplMyselfAsUnitConversion() const;
    /**
     * Destructor for the object. -> Note: we don't destroy the static
     * object here, because that would create an inifinite loop if the
     * destructor is called for the static object.
     */
    virtual ~BE_UnitConversionLength();

    /**
     * Return a string containing a short description of what the 
     * units converter does.
     */
    std::string returnUsage() const;
  };

}

#endif
