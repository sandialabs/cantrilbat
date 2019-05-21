/**
 *  @file Electrode_CSTR.h
 *     Headers for the declarations of the Electrode_CSTR class, used to model 
 *     Electrode processes in particles with no transport limits hard-coded for a MCMB Anode
 *     (see \ref electrode_mgr and class \link Zuzax::Electrode_CSTR Electrode_CSTR\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_CSTR_MCMBANODE_H
#define _ELECTRODE_CSTR_MCMBANODE_H

#include "Electrode_CSTR.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace Zuzax
{
//==================================================================================================================================
//! Electrode_CSTR class is an electrode that models a particle as a CSTR
/*!
 *  The base class is close to being a CSTR. This class will have to ensure
 *  that there are no plateaus that create complicated morphologies. There
 *  will just be an interfacial reaction on the exterior of the particle.
 */
class Electrode_CSTR_MCMBAnode : public Electrode_CSTR
{
public:

    //! Constructor
    Electrode_CSTR_MCMBAnode();

    //! Destructor
    virtual ~Electrode_CSTR_MCMBAnode();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    Electrode_CSTR_MCMBAnode(const Electrode_CSTR_MCMBAnode& right);

    //! Assignment operator
    /*!
     *  @param[in]           right               object to be copied
     *
     *  @return                                  Returns a reference to the current object
     */
    Electrode_CSTR_MCMBAnode& operator=(const Electrode_CSTR_MCMBAnode& right);
 
    //! Duplicator function
    /*!
     *  Duplicate the current Electrode object, returning a base Electrode pointer
     *
     *  @return                                   Returns a duplicate of the current object as a base class pointer
     */
    virtual Electrode* duplMyselfAsElectrode() const override;

    //! Return the type of electrode
    /*!
     *  Returns the enum type of the electrode. This is used in the factory routine.
     *
     *  @return Returns an enum type, called   Electrode_Types_Enum
     */
    virtual Electrode_Types_Enum electrodeType() const override;

    //!  Setup the electrode
    /*!
     * @param[in]            ei                  ELECTRODE_KEY_INPUT pointer object
     *
     * @return                                   Returns 0 if successful, -1 if not
     */
    virtual int electrode_model_create(ELECTRODE_KEY_INPUT* ei) override;

    //! Calculate the relative extent of reaction from the current state of the object
    /*!
     *  (virtual from Electrode.h)
     *
     *  Calculate the relative extent of reaction from the final state, spmoles_final.
     *  This is a virtual function because there is no way to do this except by knowing about
     *  the system.
     *  The relative extent of reaction is a dimensionless number that varies. It doesn't
     *  always vary between 0 and 1. Sometimes there are Li's that can be reacted or sites
     *  that can't be filled with Li....
     *
     *  @return                                  returns the relative extent of reaction (dimensionless).
     */
    virtual double calcRelativeExtentRxn_final() const override;

    //! Set the final state of the electrode using the relExtentRxn
    /*!
     *  This sets the state of the system, i.e., spmoles_final_[] for the solid phase
     *  components of the electrode using a single number.
     *
     *  It is virtual because there is no way to do this except by knowing about the system
     *
     *  It must be the case that  calcRelativeExtentRxn_final() and setState_relativeExtentRxn()
     *  are inverses of one another. Note, this means that if the state of the system has more than one rank,
     *  then the other ranks are unperturbed by the round trip.
     *
     *  @param[in]           relExtentRxn        input of the relative extent of reaction
     */
    virtual void setState_relativeExtentRxn(double relExtentRxn) override;

private:
    //! global index of the lithium containing species 
    size_t ig_SolidLi_;

    //! globa index of the vacancy-type species
    size_t ig_SolidV_;

    //! global index of the MCMB phase in the phaselist
    size_t ip_MCMB_;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

