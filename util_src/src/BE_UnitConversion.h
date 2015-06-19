/**
 * @file BE_UnitConversion.h
 *  Header for the base object that handles units conversion
 *  (see \ref blockentryModule and class
 *  \link BEInput::BE_UnitConversion BE_UnitConversion\endlink).
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

namespace BEInput
{

/**
 * @defgroup unitconversion Automatic Unit Conversion
 *
 * These classes conduct automatic unit conversion for
 * double-like data within the input system.
 *
 * All units are internally converted into the MKS system.
 *
 *     The main member function is toSI(). Given a string representing units,
 *    this function returns a mulitplicative constant that can be used to change
 *    the value given in those units into the MKS unit system.
 *
 *
 */

//!  Base class for the unit conversion system
/*!
 *    The main member function is toSI(). Given a string representing units,
 *    this function returns a mulitplicative constant that can be used to change
 *    the value given in those units into the MKS unit system.
 */
class BE_UnitConversion
{
public:

    //! Default constructor
    /*!
     * In the constructor, we define all of the conversion amounts.
     */
    BE_UnitConversion();

    //! Copy constructor
    /*!
     * In the constructor, we define all of the conversion amounts.
     *
     *  @param[in] right Object to be copied
     */
    BE_UnitConversion(const BE_UnitConversion& right);

    //! Assignment operator
    /*!
     *  @param[in] right Object to be copied
     */
    BE_UnitConversion& operator=(const BE_UnitConversion& right);

    //! Duplicator function
    virtual BE_UnitConversion* duplMyselfAsUnitConversion() const;

    //! Destructor for the object.
    virtual ~BE_UnitConversion();

    //! The function converts a string expression into a multiplicative conversion value
    /*!
     * The following syntax is used
     * / -> stands for division
     * * -> stands for multiplication
     * [1-9] -> stands for the powers of the previous string
     *          expression.
     *
     * @param unitString String input that needs to be converted
     *                   into a multiplier
     *
     * @return           double, the multiplicative value to turn the units into the MKS system
     */
    virtual double toSI(std::string unitString) const;

   
    //! Return a string containing  a short description of what the units converter does
    /*!
     * @return             returns a string
     */
    virtual std::string returnUsage() const;

protected:

    //! Map from the string to the multiplicative conversion amount.
    mutable std::map<std::string, double> m_u;
};
}
#endif
