/**
 * @file BE_UnitConversionConcentration.cpp
 *
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */


#include "BE_UnitConversionEnergy.h"

namespace BEInput
{


BE_UnitConversionEnergy::BE_UnitConversionEnergy(const BE_UnitConversionEnergy& right) :
    BE_UnitConversion(right)
{
}

BE_UnitConversionEnergy&
BE_UnitConversionEnergy::operator=(const BE_UnitConversionEnergy& right)
{
    if (this != &right) {
        BE_UnitConversion::operator=(right);
    }
    return *this;
}

BE_UnitConversion* BE_UnitConversionEnergy::duplMyselfAsUnitConversion() const
{
    BE_UnitConversionEnergy* bec = new BE_UnitConversionEnergy(*this);
    return (BE_UnitConversion*) bec;
}



/*********************************************************************
 *
 * returnUsage():
 *
 *  return a short string describing the conversion usage.
 */
std::string BE_UnitConversionEnergy::returnUsage() const
{
    std::string hhh;
    hhh = "Energy Units (J/kmol kJ/kmol kJ/gmol kcals/gmol ...)";
    hhh += " (Default = J/kmol)";
    return hhh;
}


/*********************************************************************
 *
 * toSI
 */
double BE_UnitConversionEnergy::toSI(std::string unitString) const
{
    /*
     * right now this is a shell. However, we will add
     * checks to make sure that the string is in fact a well
     * formed expression representing energy units.
     */
    double f = BE_UnitConversion::toSI(unitString);
    return f;
}

/*********************************************************************
 *
 * BE_UnitConversionEnergy():
 *
 * Constructor for the object. ->
 *  Here, we just call the base class constructor. The base map
 *  is fine for our purposes.
 */
BE_UnitConversionEnergy::BE_UnitConversionEnergy()  :
    BE_UnitConversion()
{

}

/*********************************************************************
 *
 * ~BE_UnitConversionEnergy():
 *
 * Destructor for the object. -> Note: we don't destroy the static
 * object here, because that would create an inifinite loop if the
 * destructor is called for the static object.
 */
BE_UnitConversionEnergy::~BE_UnitConversionEnergy()
{
}
/********************************************************************/
}
