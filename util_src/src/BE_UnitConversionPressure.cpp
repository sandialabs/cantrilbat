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

#include "BE_UnitConversionPressure.h"

namespace BEInput {


  BE_UnitConversionPressure::BE_UnitConversionPressure(const BE_UnitConversionPressure &right) :
    BE_UnitConversion(right)
  {
  }

  BE_UnitConversionPressure& 
  BE_UnitConversionPressure::operator=(const BE_UnitConversionPressure& right) {
    if (this != &right) {
      BE_UnitConversion::operator=(right);
    }
    return *this;
  }

  BE_UnitConversion * BE_UnitConversionPressure::duplMyselfAsUnitConversion() const {
    BE_UnitConversionPressure *bec = new BE_UnitConversionPressure(*this);
    return (BE_UnitConversion *) bec;
  }

 

  /*********************************************************************
   *
   * returnUsage():
   *
   *  return a short string describing the conversion usage.
   */
  std::string BE_UnitConversionPressure::returnUsage() const {
    std::string hhh;
    hhh = "Pressure Units (Pa, torr, atm)";
    hhh += " (Default = Pa)";
    return hhh;
  }


  /*********************************************************************
   *
   * BE_UnitConversionPressure():
   *
   * Constructor for the object. -> 
   *  Here, we just call the base class constructor. The base map
   *  is fine for our purposes.
   */
  BE_UnitConversionPressure::BE_UnitConversionPressure() :
    BE_UnitConversion()
  {
  }

  /*********************************************************************
   *
   * ~BE_UnitConversionPressure():
   *
   * Destructor for the object. -> Note: we don't destroy the static
   * object here, because that would create an inifinite loop if the
   * destructor is called for the static object.
   */
  BE_UnitConversionPressure::~BE_UnitConversionPressure()
  {
  }
  /********************************************************************/
}
