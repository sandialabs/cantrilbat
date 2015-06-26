/**
 * @file BE_UnitConversionPressure.h
 *
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef BE_UNITCONVERSIONPRESSURE_H
#define BE_UNITCONVERSIONPRESSURE_H

#include "BE_UnitConversion.h"

namespace BEInput
{

//!  Unit conversion operator for quantities that are of type pressure
/*!
 *    The member function toSi() converts pressure unit strings to a conversion factor for pascals.
 */
class BE_UnitConversionPressure : public BE_UnitConversion
{
public:

    //! Constructor
    BE_UnitConversionPressure();

    //! Constructor is defined to be protected so that we can't make multiple
    //! copies of this, yet we can inherit from this object.
    /*!
     *     In the constructor, we define all of the conversion amounts.
     *
     *    @param[in]  right             Object to be copied
     */
    BE_UnitConversionPressure(const BE_UnitConversionPressure& right);

    //! Assignment operator
    /*!
     *  @param right Object to be copied
     *
     *  @return                   Returns a reference to the current object
     */
    BE_UnitConversionPressure& operator=(const BE_UnitConversionPressure& right);

    //! Duplicator
    /*!
     *  @return                   Returns a pointer to a copy of the current object
     */
    virtual BE_UnitConversion* duplMyselfAsUnitConversion() const;

    /**
     * Destructor for the object. -> Note: we don't destroy the static
     * object here, because that would create an inifinite loop if the
     * destructor is called for the static object.
     */
    virtual ~BE_UnitConversionPressure();

    //! Return a string containing a short description of what the units converter does.
    /*!
     *   @return                  Returns a description as a string.
     */
    std::string returnUsage() const;

private:


};

}

#endif
