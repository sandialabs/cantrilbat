/**
 * @file BE_UnitConversionLength.cpp
 *   Definitions for the object that handles unit conversions based on length or powers of the length
 *  ((see \ref blockentryModule and class \link BEInput::BE_UnitConversionLength BE_UnitConversionLength\endlink).
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "BE_UnitConversionLength.h"

//--------------------------------------------------------------------------------------------------------------------------------
namespace BEInput
{
//================================================================================================================================
 /*
  * BE_UnitConversionLength():
  *
  * Constructor for the object. ->
  *  Here, we just call the base class constructor. The base map
  *  is fine for our purposes.
  */
BE_UnitConversionLength::BE_UnitConversionLength() :
    BE_UnitConversion()
{
    m_u.clear();
    // length
    m_u["m"]        = 1.0;
    m_u["cm"]       = 0.01;
    m_u["km"]       = 1.0e3;
    m_u["mm"]       = 1.0e-3;
    m_u["micron"]   = 1.0e-6;
    m_u["microns"]  = 1.0e-6;
    m_u["nm"]       = 1.0e-9;
    m_u["A"]        = 1.0e-10;
    m_u["inches"]   = 0.0254;
    m_u["ft"]       = 0.3048;
}
//================================================================================================================================
BE_UnitConversionLength::BE_UnitConversionLength(const BE_UnitConversionLength& right) :
    BE_UnitConversion(right)
{
}
//================================================================================================================================
BE_UnitConversionLength& BE_UnitConversionLength::operator=(const BE_UnitConversionLength& right)
{
    if (this != &right) {
        BE_UnitConversion::operator=(right);
    }
    return *this;
}
//================================================================================================================================
BE_UnitConversion* BE_UnitConversionLength::duplMyselfAsUnitConversion() const
{
    BE_UnitConversionLength* bec = new BE_UnitConversionLength(*this);
    return (BE_UnitConversion*) bec;
}
//================================================================================================================================
/*
 *
 * returnUsage():
 *
 *  return a short string describing the conversion usage.
 */
std::string BE_UnitConversionLength::returnUsage() const
{
    std::string hhh;
    hhh = "Length Units (m cm mm um nm inches ft)";
    hhh += " (Default = m)";
    return hhh;
}
//================================================================================================================================
BE_UnitConversionLength::~BE_UnitConversionLength()
{
}
//================================================================================================================================
}
//--------------------------------------------------------------------------------------------------------------------------------
