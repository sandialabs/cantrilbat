/**
 *  @file RSD_OCVmodel.h
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

// MCMB 2528 graphite measured by Chris Bogatu 2000, 
//           Telcordia and PolyStor materials.
//    Modified May 2003 to match data from Joongpyo Shim
//    for 0.01 < x < 0.99

#define OCVAnode_MCMB2528                101


#define OCVCathode_MCMB2528              301

//!  create a map
/*!
 *
 */
extern void createOCVmodel_map(std::map<int, std::string>& smap);

//==================================================================================
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
 *   Relative Extent of reaction variable
 *   -------------------------------------
 *      This denotes the relative extent of reaction of the solid phase.  The model will not be restricted
 *      from this variable going strictly between 0 and 1. It may go over a more limited range, or a slightly
 *      larger range. However, it is unitless with an order of magnitude of 1. 
 *
 *      The usual setup will be to specify the mole fraction of aa particular species as being the relative
 *      depth of discharge.
 *
 *       
 */
class RSD_OCVmodel
{
  public:

    //! Constructor, with a modelID as a parameter
    /*!
     *   @param modelId  Parameter that specifies the model
     *                   
     */
    RSD_OCVmodel(int modelId = -1);

    //! Copy constructor
    /*!
     *   @param right   Object to be copied
     */
    RSD_OCVmodel(const RSD_OCVmodel& right);

    //! destructor
    virtual ~RSD_OCVmodel();

    //! Assignment operator
    /*!
     *   @param right   Object to be copied
     *
     *  @return returns a reference to the current object
     */
    RSD_OCVmodel& operator=(const  RSD_OCVmodel& right);

    //! Duplicator function for this class
    /*!
     *  @return Returns a duplication of the current state as a pointer to the base class
     */
    virtual RSD_OCVmodel* duplMyselfAsOCVmodel(ThermoPhase* solidPhase = 0) const;

    void assignShallowPointers(ThermoPhase* solidPhase);

    //! Set how the relative extent calculation is carried out
    /*!
     *      The usual setup will be to specify the mole fraction of a particular species as being the relative
     *      depth of discharge.
     *
     *   @param tp                Pointer to the ThermoPhase class
     *   @param kspec             species whose mole fraction will be assigned as the rel extent
     *   @param dvec              double vector for later expansion
     *   @param ivec              int vector for later expansion
     */
    virtual void setup_RelExtent(ThermoPhase *tp, size_t kspec, double *dvec = 0, int *ivec = 0);

    //!  Return the open circuit voltage given the relative extent of reaction
    /*!
     *   @returns Returns the open circuit voltage at the current relative extent of reaction
     */
    virtual double OCV_value() const;

    //!  Return the derivative of the open circuit voltage wrt the relative extent of reaction
    /*!
     *   @returns Return the derivative of the open circuit voltage wrt the relative extent of reaction
     */
    virtual double OCV_dvaldExtent() const;

    //!  Return the derivative of the open circuit voltage wrt the relative extent of reaction
    /*!
     *   @returns Return the derivative of the open circuit voltage wrt the relative extent of reaction
     */
    virtual double OCV_dvaldT() const;

    //! Return the model name
    virtual std::string modelName() const;

    //! Return the relative extent
    virtual double RelExtent() const;

protected:
    //!  Calculate the relative extent of reaction
    /*!
     *  Internal routine that stores the result in relExtent_;
     */
    virtual void calcRelExtent() const;

    // - ----------------------------------------------------------------------------------------------------------------
    // --------------------------------------------      DATA       -----------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------
    //!  Model id
    int modelID_;

    //! model name
    std::string modelName_;
public:
    //! underlying ThermoPhase model
    ThermoPhase* solidPhaseModel_;
protected:
    //! Particular species within the ThermoPhase which indicates the relative depth of discharge
    /*!
     *    The usual setup will be to specify the mole fraction of aa particular species as being the relative
     *    depth of discharge
     */
    size_t kSpecies_DoD_;

    //! Value of the relative extent
    mutable double relExtent_;

    //! Vector of mole fractions from solidPhaseModel_
    mutable std::vector<double> xMF_;
    

};


}
#endif
