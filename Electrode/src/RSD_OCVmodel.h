/**
 *  @file RSD_OCVmodel.h 
 *    Declarations for the base object that contains the relative extent and  OCV calculation for models
 *    that override the normal OCV calculation based on species thermodynamics 
 *    (see \ref ExtendedPhaseGroups and class \link Cantera::RSD_OCVmodel RSD_OCVmodel\endlink).
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef RSD_OCVMODEL_H
#define RSD_OCVMODEL_H

#include "Electrode_defs.h"

#include <string>
#include <vector>
#include <map>

//! Model for OCV Override functions

namespace Cantera 
{

class ThermoPhase;
//==================================================================================================================================

//! OCV Model ID for an Anode that can be used to set an arbitrary OCV using a member function.
//!                                     
#define  OCVAnode_CONSTANT               100

//! OCV Model ID for MCMB 2528 graphite measured by Chris Bogatu 2000, 
//!                                     Telcordia and PolyStor materials.
//!                                     Modified May 2003 to match data from Joongpyo Shim for 0.01 < x < 0.99
//!
#define OCVAnode_MCMB2528                101

//! OCV Model ID for  MCMB 2510 carbon (Bellcore)
//!                                     Telcordia and PolyStor materials.
//!                                     Modified May 2003 to match data from Joongpyo Shim for 0.01 < x < 0.99
//!
#define OCVAnode_MCMB2528_dualfoil       102

//! OCV Model ID for CoO2 (Cobalt dioxide)
//!                                      Measured by Oscar Garcia 2001 using Quallion electrodes for
//!                                      0.5 < y < 0.99.  Fit revised by Karen Thomas in May 2003 to
//!                                      match Doyle's fit for y < 0.4 and Garcia's data at larger y.
//!                                      Valid for 0 < y < 0.99. Note that capacity fade is found to
//!                                      occur experimentally if y goes below 0.5; this is not included in the model.
//!
#define OCVCathode_CoO2_dualfoil         201

//!  Create a 1 to 1 map between an int OCV model number and a string OCV model name. We have global functions to go between
//!  the two values.
/*!
 *   We create a map between an int description and a string. This is kept as a global variable with the namespace.
 *   See the functions defined in Electrode_Factory.h for a description of these variables and global functions.
 *
 *   @param[in]            smap                Mapping between an int, defined above with prefixes OCV and a 
 *                                             string identifying the model when needed for printouts and input
 */
extern void createOCVmodel_map(std::map<int, std::string>& smap);

//==================================================================================================================================
//!  Base class for specification of models for the open circuit voltage override
//!  calculation
/*!
 *       These models will calculate the open circuit voltage in a different way
 *
 *
 *   Temperature and pressure are specified by the thermophase. They are not set within this interface.
 *   The determination of the relative extent of reaction must be specified by this interface. However, the relative extent of
 *   reaction can't be changed by this interface.
 *
 *   What is the definition of the open circuit voltage for this class? This is not necessarily straightforward.
 *   It involves the combination of standard state gibbs free energies and chemical potentials for the species.
 *   We use the stndard state G for the solution phase species always, but leave the electrode phases as regular G values.
 *   There is no other way to define this, as the concentration of the solution affects the open circuit voltage, but can't
 *   be a part of the calculation. This is the same definition and is equivalent to one where we assume that the other electrode
 *   is the Lithium reference electrode. 
 *   (so we could substitute Li metal thermophase for this treatment. We may in the future).
 *
 *   Take for example an anode reaction.
 *
 *      Li-Theta -> Theta + Li+ + e-
 *
 *   DeltaG0 = u0_Theta + u0_Lip + u0_em  - u0_Li-Theta
 *   DeltaG  = u_Theta + u_Lip + u_em  - u_Li-Theta
 *
 *   The open circuit voltage is defined as the following statement.
 *
 *      (1) (F) (OCV) = u_Theta + + u0_Lip + u0_em - u_Li-Theta
 *
 *   
 *
 *   Relative Extent of reaction variable
 *   -------------------------------------
 *      This denotes the relative extent of reaction of the solid phase.  The model will not be restricted
 *      from this variable going strictly between 0 and 1. It may go over a more limited range, or a slightly
 *      larger range. However, it is unitless with an order of magnitude of 1. 
 *
 *      The usual setup will be to specify the mole fraction of a particular species as being the relative
 *      depth of discharge.
 *
 *   Model names and modelID
 *   -------------------------------------
 *
 *      We maintain modelIDs and modelName strings. There is a 1 to 1 mapping between the two.
 *      To add a new model, these mappings need to be updated with the new values.
 *
 *       
 */
class RSD_OCVmodel
{
  public:

    //! Constructor, with a modelID as a parameter
    /*!
     *   @param        modelId                  Parameter that specifies the model
     *                   
     */
    RSD_OCVmodel(int modelId = -1);

    //! Copy constructor
    /*!
     *   @param        right                    Object to be copied
     */
    RSD_OCVmodel(const RSD_OCVmodel& right);

    //! Destructor
    virtual ~RSD_OCVmodel();

    //! Assignment operator
    /*!
     *  @param[in]     right                     Object to be copied
     *
     *  @return                                  Returns a reference to the current object
     */
    RSD_OCVmodel& operator=(const RSD_OCVmodel& right);

    //! Duplicator function for this class
    /*!
     *  @param[in]    solidPhase                 Pointer to the ThermoPhase class to use within the duplicate object.
     *                                           The default is NULL, which indicates that the current pointer to the existing
     *                                           ThermoPhase object will be used in the duplicate object.
     *
     *  @return                                  Returns a duplication of the current state as a pointer to the base class
     */
    virtual RSD_OCVmodel* duplMyselfAsOCVmodel(ThermoPhase* solidPhase = 0) const;

    //!  Assign the shallow pointer ThermoPhase object for the solidPhase
    /*!
     *   @param[in]    solidPhase                Pointer to the ThermoPhase class to use within this object
     */
    void assignShallowPointers(ThermoPhase* solidPhase);

    //!  Returns the pointer to the assigned solidPhase model
    /*!
     *   @return                                 Returns a const pointer to the solid phase object
     */
    const ThermoPhase* solidPhasePtr() const;

    //! Set how the relative extent calculation is carried out
    /*!
     *      The usual setup will be to specify the mole fraction of a particular species as being the relative
     *      depth of discharge.
     *
     *   @param        tp                       Pointer to the ThermoPhase class
     *   @param        kspec                    Species whose mole fraction will be assigned as the rel extent
     *   @param        dvec                     double vector for later expansion
     *   @param        ivec                     int vector for later expansion
     */
    virtual void setup_RelExtent(ThermoPhase *tp, size_t kspec, double *dvec = 0, int *ivec = 0);

    //!  Return the open circuit voltage given the relative extent of reaction
    /*!
     *   @returns                               Returns the open circuit voltage at the current relative extent of reaction
     */
    virtual double OCV_value() const;

    //!  Return the derivative of the open circuit voltage wrt the relative extent of reaction
    /*!
     *   @returns                               Return the derivative of the open circuit voltage wrt the relative extent of
     *                                          reaction
     */
    virtual double OCV_dvaldExtent() const;

    //!  Return the derivative of the open circuit voltage wrt the relative extent of reaction
    /*!
     *   @returns                               Return the derivative of the open circuit voltage wrt the relative 
     *                                          extent of reaction
     */
    virtual double OCV_dvaldT() const;

    //! Return the model name
    /*!
     *  @return                                 Returns the string model name
     */
    virtual std::string modelName() const;

    //! Return the relative extent of the reaction
    /*!
     *  One characteristic of this object is that the object defines what the relative extent variable is.
     *  This is discussed in the documentation for the class.
     *
     *  @return                                 Returns the value of the relative extent
     */
    virtual double RelExtent() const;

protected:
    //!  Calculate the relative extent of reaction as defined by this class object or child objects
    /*!
     *   Calculate the relative extent of reaction. This internal routine stores the result in relExtent_.
     */
    virtual void calcRelExtent() const;

    // -----------------------------------------------------------------------------------------------------------------------------
    // ------------------------------------------------      DATA       ------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------------------
protected:
    //!  Model id
    int modelID_;

    //!  Model name
    std::string modelName_;

    //!  Pointer to the Underlying %ThermoPhase model for the solid phase of the electrode.
    ThermoPhase* solidPhaseModel_;

    //!  Particular species within the ThermoPhase which indicates the relative depth of discharge
    /*!
     *   The usual setup will be to specify the mole fraction of a particular species as being the determiner of the relative
     *   depth of discharge. This works for the normal case of a two species solid phase model. However, it may not work
     *   for more complicated cases, in which case this species index becomes undefined or redefined.
     */
    size_t kSpecies_DoD_;

    //!  Value of the relative extent
    /*!
     *   For calculations related to OCV calculations, we keep an intermediate value of the relative extent.
     */
    mutable double relExtent_;

    //!  Vector of mole fractions from solidPhaseModel_
    /*!
     *   For computations related to the calculation of relative extents and OCV values, 
     *   we keep a copy of the solid phase mole fractions here.
     */
    mutable std::vector<double> xMF_;

    //! Vector available for use by models
    std::vector<double> dvec_;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
