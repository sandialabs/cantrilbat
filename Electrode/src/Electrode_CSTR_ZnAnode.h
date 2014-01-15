/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_CSTR_ZNANODE_H
#define _ELECTRODE_CSTR_ZNANODE_H

#include "Electrode_CSTR.h"

class ELECTRODE_KEY_INPUT;

namespace Cantera
{

//! Child CSTR class for Zn anode in alkaline batteries.
class Electrode_CSTR_ZnAnode : public Electrode_CSTR
{
public:
    //! Constructor
    Electrode_CSTR_ZnAnode();

    //! Destructor
    virtual ~Electrode_CSTR_ZnAnode();

    //! Copy Constructor
    Electrode_CSTR_ZnAnode(const Electrode_CSTR_ZnAnode& right);

    //! Assignment operator
    Electrode_CSTR_ZnAnode& operator=(const Electrode_CSTR_ZnAnode& right);

    //! Return the type of electrode
    virtual Electrode_Types_Enum electrodeType() const;

    //!  Setup the electrode
    virtual int electrode_model_create(ELECTRODE_KEY_INPUT* ei);


    //! Calculate the relative extent of reaction from the current state of the object
    /*!
     *  Calculate the relative extent of reaction from the final state, spmoles_final.
     *  This is a virtual function because there is no way to do this except by knowing about
     *  the system.
     *  The relative extent of reaction is a dimensionless number that varies. It doesn't
     *  always vary between 0 and 1. Sometimes there are Li's that can be reacted or sites
     *  that can't be filled with Li....
     *
     *  @return returns the relative extent of reaction (dimensionless).
     */
    virtual double calcRelativeExtentRxn_final() const;

    //! Set the final state of the electrode using the relExtentRxn
    /*!
     *  This sets the state of the system, i.e., spmoles_final_[] for the solid phase
     *  components of the electrode using a single number.
     *
     *  It is virtual because there is no way to do this except by knowing about the system
     *
     *  It must be the case that  calcRelativeExtentRxn_final() and etStateFinal_fromRelativeExtentRxn()
     *  are inverses of one another. Note, this means that if the state of the system has more than one rank,
     *  then the other ranks are unperturbed by the round trip.
     *
     *  @param relExtentRxn  input of the relative extent of reaction
     */
    virtual void setStateFinal_fromRelativeExtentRxn(double relExtentRxn);


private:
    int zn_phase_index;
    int zn_sp_index;
};

}

#endif
