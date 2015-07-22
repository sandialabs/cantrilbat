/**
 * @file BE_UnitConversionConcentration.h
 *   Declarations for the UnitConversion class that handles concentrations 
 *   (see \ref blockentryModule and class \link BEInput::BE_UnitConversionConcentration BE_UnitConversionConcentration\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef BE_UNITCONVERSIONCONCENTRATION_H
#define BE_UNITCONVERSIONCONCENTRATION_H

#include "BE_UnitConversion.h"
//----------------------------------------------------------------------------------------------------------------------------------
namespace BEInput
{

//==================================================================================================================================
//! Class to group units of type concentration
/*!
 *  @ingroup unitconversion
 *
 *       This class is used to convert units of the concentration (kmol / m3).
 *                  
 */
class BE_UnitConversionConcentration : public BE_UnitConversion
{
public:

    //! Default constructor
    BE_UnitConversionConcentration();

    //! Copy constructor
    /*!
     * In the constructor, we define all of the conversion amounts.
     *
     * @param[in]   right             
     */
    BE_UnitConversionConcentration(const BE_UnitConversionConcentration& right);

    //! Assignment operator
    /*!
     *  @param[in] right              Object to be copied
     *
     * @return                        returns a reference to the current object
     */
    BE_UnitConversionConcentration& operator=(const BE_UnitConversionConcentration& right);

    //! Duplicator
    /*!
     *   @return                      Returns a pointer to the duplicate of the current object
     */
    virtual BE_UnitConversion* duplMyselfAsUnitConversion() const;

    //! Destructor
    virtual ~BE_UnitConversionConcentration();
    
     //! Return a string containing a short description of what the units converter does.
     /*!
      * @return                   Returns a descriptive string
      */
    std::string returnUsage() const;

private:
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
