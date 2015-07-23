/**
 * @file BE_UnitConversionLength.h
 *   Declarations for the object that handles unit conversions based on length or multiples of the length
 *   ((see \ref blockentryModule and class \link BEInput::BE_UnitConversionLength BE_UnitConversionLength\endlink).
 *
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
//----------------------------------------------------------------------------------------------------------------------------------
namespace BEInput
{
//==================================================================================================================================
//! Unit conversion object for quantities that have units of length, returning values in meters
/*!
 *    This may include powers of the length.
 */
class BE_UnitConversionLength : public BE_UnitConversion
{
public:

    //! Default constructor
    BE_UnitConversionLength();

    //! Copy Constructor
    /*!
     * In the constructor, we define all of the conversion amounts.
     *
     *  @param[in]   right        Object to be copied
     */
    BE_UnitConversionLength(const BE_UnitConversionLength& right);

    //! Assignment operator
    /*!
     *  @param[in] right          Object to be copied
     *
     *  @return          Returns a reference to the current object
     */
    BE_UnitConversionLength& operator=(const BE_UnitConversionLength& right);

    //! Duplicator
    /*!
     *  @return          Returns a duplicate of the current object
     */
    virtual BE_UnitConversion* duplMyselfAsUnitConversion() const;

    //! Destructor for the object
    /*!
     *      Note: we don't destroy the static object here, because that would create an inifinite loop if the
     *            destructor is called for the static object.
     */
    virtual ~BE_UnitConversionLength();
   
    //!Return a string containing a short description of what the
    /*!
     *  @return             Returns a string containing the description
     */
    std::string returnUsage() const;
};
//==============================================================================================================================
}
//------------------------------------------------------------------------------------------------------------------------------
#endif
