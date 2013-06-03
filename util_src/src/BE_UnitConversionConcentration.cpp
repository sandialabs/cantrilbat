/**
 * @file BE_UnitConversionConcentration.cpp
 *
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


#include "BE_UnitConversionConcentration.h"

namespace BEInput {


  BE_UnitConversionConcentration::BE_UnitConversionConcentration(const BE_UnitConversionConcentration &right) :
    BE_UnitConversion(right)
  {
  }

  BE_UnitConversionConcentration& 
  BE_UnitConversionConcentration::operator=(const BE_UnitConversionConcentration& right) {
    if (this != &right) {
      BE_UnitConversion::operator=(right);
    }
    return *this;
  }

  BE_UnitConversion * BE_UnitConversionConcentration::duplMyselfAsUnitConversion() const {
    BE_UnitConversionConcentration *bec = new BE_UnitConversionConcentration(*this);
    return (BE_UnitConversion *) bec;
  }

 

  /*********************************************************************
   *
   * returnUsage():
   *
   *  return a short string describing the conversion usage.
   */
  std::string BE_UnitConversionConcentration::returnUsage() const {
    std::string hhh;
    hhh = "Concentration Units (Pa, torr, atm)";
    hhh += " (Default = Pa)";
    return hhh;
  }


  /*********************************************************************
   *
   * BE_UnitConversionConcentration():
   *
   * Constructor for the object. -> 
   *  Here, we just call the base class constructor. The base map
   *  is fine for our purposes.
   */
  BE_UnitConversionConcentration::BE_UnitConversionConcentration()  
  {
    //length
    m_u["m"]        = 1.0;
    m_u["cm"]       = 0.01;
    m_u["km"]       = 1.0e3;
    m_u["mm"]       = 1.0e-3;
    m_u["micron"]   = 1.0e-6;
    m_u["nm"]       = 1.0e-9;
    m_u["A"]        = 1.0e-10;
    m_u["inches"]   = 0.0254;
    m_u["ft"]       = 0.3048;
    const double Avogadro =  6.022136736e26;

    // quantity
    m_u["mol"]      = 1.0e-3;
    m_u["gmol"]     = 1.0e-3;
    m_u["mole"]     = 1.0e-3;
    m_u["kmol"]     = 1.0;
    m_u["molec"]    = 1.0/Avogadro;
  }

  /*********************************************************************
   *
   * ~BE_UnitConversionConcentration():
   *
   * Destructor for the object. -> Note: we don't destroy the static
   * object here, because that would create an inifinite loop if the
   * destructor is called for the static object.
   */
  BE_UnitConversionConcentration::~BE_UnitConversionConcentration()
  {
  }

}
