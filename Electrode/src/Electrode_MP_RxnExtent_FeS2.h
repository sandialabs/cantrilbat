/**
 * @file Electrode_MP_RxnExtent.h
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

//----------------------------------------------------------------------------------------------------------------------------------
namespace Zuzax
{
//==================================================================================================================================
//! Reaction extent model constructed for the FeS2 thermal battery system
/*!
 *
 */
class Electrode_MP_RxnExtent_FeS2 : public Zuzax::Electrode_MP_RxnExtent
{
public:

    //! Constructor
    Electrode_MP_RxnExtent_FeS2();

    //! Destructor
    virtual ~Electrode_MP_RxnExtent_FeS2();

    //! Copy Constructor
    /*!
     * @param[in]            right               Object to be copied
     */
    Electrode_MP_RxnExtent_FeS2(const Electrode_MP_RxnExtent_FeS2& right);

    //! Assignment operator
    /*!
     *  @param[in]           right               object to be copied
     *
     *  @return                                  Returns a reference to the current object
     */
    Electrode_MP_RxnExtent_FeS2& operator=(const Electrode_MP_RxnExtent_FeS2& right);

    //! Return the type of electrode
    /*!
     *  Returns the enum type of the electrode. This is used in the factory routine.
     *
     *  @return                                  Returns an enum type, called   Electrode_Types_Enum
     */
    virtual Electrode_Types_Enum electrodeType() const override;
  
    //! Returns the equilibrium standard state open circuit voltage for the current conditions 
    /*!
     *  Returns the standard state open circuit voltage within the current region at the relative
     *  extent of reaction. The region parameter takes precedence and no checking is done.
     *
     *  (virtual function of Electrode_MP_RxnExtent.h)
     *
     *   @param[in]          relativeExtentRxn   Relative extent of reaction
     *   @param[in]          xRegion             Integer indicating the region of the extent of reaction
     *
     *   @return                                Returns the standard state voltage.
     */
    virtual double openCircuitVoltageSS_Region_NoCheck(double relativeExtentRxn, int xRegion) const override;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

