/**
 * @file BE_UnitConversionEnergy.h
 *      Declarations for the class  BE_UnitConversionEnergy 
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef BE_UNITCONVERSIONENERGY_H
#define BE_UNITCONVERSIONENERGY_H

#include "BE_UnitConversion.h"
//----------------------------------------------------------------------------------------------------------------------------------
namespace BEInput
{
//==================================================================================================================================
//! Class to group units of type molar energy
/*!
 *  @ingroup unitconversion
 *
 *       This class is used to convert units of type molar energy (kg m2 / s2 / kmol ).
 *                  
 */
class BE_UnitConversionEnergy : public BE_UnitConversion
{
public:

    //! Constructor
    BE_UnitConversionEnergy();

    //! Copy constructor
    /*!
     * In the constructor, we define all of the conversion amounts.
     *
     * @param[in]  right       Object to be copied
     */
    BE_UnitConversionEnergy(const BE_UnitConversionEnergy& right);

    //! Assignment operator
    /*!
     *  @param    right        Object to be copied
     *
     *  @return                Returns a reference to the current object
     */
    BE_UnitConversionEnergy& operator=(const BE_UnitConversionEnergy& right);

    //! Duplicator
    virtual BE_UnitConversion* duplMyselfAsUnitConversion() const;

    //! Destructor
    virtual ~BE_UnitConversionEnergy();
    
    //! Return a string containing a short description of what the units converter does.
    /*!
     * @return    Returns a string containing the usage
     */
    std::string returnUsage() const;

    //! Does the actual conversion to SI units based on the unitString value
    //! The function converts a string expression into a multiplicative conversion value
    /*!
     * The following syntax is used
     * / -> stands for division
     * * -> stands for multiplication
     * [1-9] -> stands for the powers of the previous string
     *          expression.
     *
     * @param unitString String input that needs to be converted into a multiplier
     *
     * @return           double, the multiplicative value to turn the units into the MKS system
     */
    virtual double toSI(std::string unitString) const;

};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
