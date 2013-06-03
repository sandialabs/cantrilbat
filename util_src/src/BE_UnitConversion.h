/**
 * @file BE_UnitConversion.h
 *  Header for the base object that handles units conversion
 *  (see \ref blockentryModule and class 
 *  \link BEInput::BE_UnitConversion BE_UnitConversion\endlink).
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

#ifndef BE_UNITCONVERSION_H
#define BE_UNITCONVERSION_H

#include <string>
#include <map>

namespace BEInput {

  /**
   * @defgroup unitconversion Automatic Unit Conversion
   *
   * These classes conduct automatic unit conversion for
   * double-like data within the input system.
   * 
   * All units are internally converted into the MKS system.
   *
   * 
   *
   *
   */
  class BE_UnitConversion {
  public:
    //! Constructor is defined to be protected so that we can't make multiple
    //! copies of this, yet we can inherit from this object.
    /*!
     * In the constructor, we define all of the conversion amounts.
     */
    BE_UnitConversion();

    //! Constructor is defined to be protected so that we can't make multiple
    //! copies of this, yet we can inherit from this object.
    /*!
     * In the constructor, we define all of the conversion amounts.
     */
    BE_UnitConversion(const BE_UnitConversion &right);

    //! Assignment operator
    /*!
     *  @param right Object to be copied
     */
    BE_UnitConversion& operator=(const BE_UnitConversion& right);

    virtual BE_UnitConversion * duplMyselfAsUnitConversion() const;
    /**
     * This static function creates a one and only one instantiation
     * of this object.
     */
    // static BE_UnitConversion* units();

    /**
     * Use this function to delete the static pointer.
     * Note you can not use the destructor to delete or an
     * infinite loop results when the destructor for the static
     * object is called.
     */
    //virtual void deleteStatic();

    /**
     * Destructor for the object. -> Note: we don't destroy the static
     * object here, because that would create an inifinite loop if the
     * destructor is called for the static object.
     */
    virtual ~BE_UnitConversion();
   
     //! The function converts a string expression to 
     //! a conversion value.
     /*!
      * The following syntax is used
      * / -> stands for division
      * * -> stands for multiplication
      * [1-9] -> stands for the powers of the previous string
      *          expression.
      *
      * @param unitString String input that needs to be converted
      *                   into a multiplier
      */
    virtual double toSI(std::string unitString) const;

    /**
     * Return a string containing  a short description of what the
     *  units converter does.
     */
    virtual std::string returnUsage() const;

  protected:
    
    //! Map from the string to the multiplicative conversion amount.
    mutable std::map<std::string, double> m_u;
   
 
  };
}
#endif
