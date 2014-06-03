/*
 * $Id: Electrode_CSTR_LiCoO2Cathode.h 571 2013-03-26 16:44:21Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_CSTR_LICOO2CATHODE_H
#define _ELECTRODE_CSTR_LICOO2CATHODE_H


#include "Electrode_CSTR.h"


class ELECTRODE_KEY_INPUT;

namespace Cantera
{



//! Electrode_CSTR class is an electrode that models a particle as a CSTR
/*!
 *  The base class is close to being a CSTR. This class will have to ensure
 *  that there are no plateaus that create complicated morphologies. There
 *  will just be an interfacial reaction on the exterior of the particle.
 */
class Electrode_CSTR_LiCoO2Cathode : public Electrode_CSTR
{
public:
    //! Constructor
    Electrode_CSTR_LiCoO2Cathode();

    //! Destructor
    virtual ~Electrode_CSTR_LiCoO2Cathode();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    Electrode_CSTR_LiCoO2Cathode(const Electrode_CSTR_LiCoO2Cathode& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    Electrode_CSTR_LiCoO2Cathode& operator=(const Electrode_CSTR_LiCoO2Cathode& right);

    //! Return the type of electrode
    /*!
     *  Returns the enum type of the electrode. This is used in the factory routine.
     *
     *  @return Returns an enum type, called   Electrode_Types_Enum
     */
    virtual Electrode_Types_Enum electrodeType() const;


    void setCapacityCoeff_LiCoO2();

    //!  Setup the electrode
    /*!
     * @param ei    ELECTRODE_KEY_INPUT pointer object
     */
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
    virtual void setState_relativeExtentRxn(double relativeExtentRxn);


    int Global_LiCoO2_Model_;


private:
    int ig_SolidLi_;
    int ig_SolidV_;

    //! global index of the LiCoO2 phase in the phaselist
    int ip_LiCoO2_;
};

}



#endif
/*****************************************************************************/
