/*
 * $Id: Electrode_MP_RxnExtent.h 571 2013-03-26 16:44:21Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_MP_RXNEXTENT_FES2_H
#define _ELECTRODE_MP_RXNEXTENT_FES2_H


#include "Electrode_MP_RxnExtent.h"

namespace Cantera
{



class Electrode_MP_RxnExtent_FeS2 : public Cantera::Electrode_MP_RxnExtent
{
public:

    //! Constructor
    Electrode_MP_RxnExtent_FeS2();

    //! Destructor
    virtual ~Electrode_MP_RxnExtent_FeS2();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    Electrode_MP_RxnExtent_FeS2(const Electrode_MP_RxnExtent_FeS2& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    Electrode_MP_RxnExtent_FeS2& operator=(const Electrode_MP_RxnExtent_FeS2& right);


    //! Return the type of electrode
    /*!
     *  Returns the enum type of the electrode. This is used in the factory routine.
     *
     *  @return Returns an enum type, called   Electrode_Types_Enum
     */
    virtual Electrode_Types_Enum electrodeType() const;

  
    //! Returns the equilibrium standard state open circuit voltage for the current conditions (virtual)
    /*!
     *  Returns the standard state open circuit voltage within the current region at the relative
     *  extent of reaction. The region parameter takes precedence and no checking is done.
     *
     *  (virtual function of Electrode_MP_RxnExtent.h)
     *
     *   @param relativeExtentRxn      Relative extent of reaction
     *   @param xRegion                Integer indicating the region of the extent of reaction
     *
     *   @return Returns the standard state voltage.
     */
    virtual double openCircuitVoltageSS_Region_NoCheck(double relativeExtentRxn, int xRegion) const;
};

// ------------------------- End Namespace Cantera -------------------------------------
}


#endif
/*****************************************************************************/
                                                                                
