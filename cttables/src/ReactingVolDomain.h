/**
 * @file ReactingVolDomain.h
 *
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef REACTINGVOLDOMAIN_H
#define REACTINGVOLDOMAIN_H

#include "zuzax/kinetics/GasKinetics.h"
#include "zuzax/kinetics/InterfaceKinetics.h"
#include "zuzax/transport.h"

#include <string>
#include <iostream>

//----------------------------------------------------------------------------------------------------------------------------------
namespace Zuzax
{

// Forward references
class Transport;
class PhaseList;

//==================================================================================================================================
//! Utility class that wraps PhaseList and adds homogeneous or heterogeneous kinetics to it
/*!
 *  
 */
class ReactingVolDomain
{
public:
    //! Base constructor
    ReactingVolDomain();

    //! Virtual destructor
    virtual ~ReactingVolDomain();

    //! Import the ReactingVolDomain from a single XML source
    /*!
     *  @param[in]           root                XML root of the file from which we will gather all phases
     */    
    bool importFromXML(XML_Node& root);

    //! Import and setup Surface Kinetics
    /*!
     *  @param[in]           pl                  Pointer to the input PhaseList
     *  @param[in]           iphSurKin           Surface phase index that owns the kinetics
     *
     *  @return                                  Returns true if InterfaceKinetics is setup and ready to go
     */
    bool importSurKinFromPL(PhaseList* pl, size_t iphSurKin);

    //! Import and setup Volume  Kinetics, but not both
    /*!
     *  @param[in]           pl                  Pointer to the input PhaseList
     *  @param[in]           iphVolKin           volume phase index that owns the kinetics
     *
     *  @return                                  Returns true if Kinetics is setup and ready to go
     */
    bool importVolKinFromPL(PhaseList* pl, size_t iphVolKin);

    //! Print out the information in a ReactingVolumeDomain
    friend std::ostream& operator<<(std::ostream& s, ReactingVolDomain& vd);

    //! Return the ThermoPhase Reference for the nth phase
    ThermoPhase& thermo(int n=0)
    {
        return *(tpList[n]);
    }

    //! This kinetics operator is associated with just one homogeneous phase, associated with tpList[0] phase
    /*!
     * This object owns the Kinetics object
     */
    Kinetics* m_kinetics;

    //! This kinetics operator is associated with multiple homogeneous and surface phases.
    /*!
     * This object owns the Kinetics object
     */
    InterfaceKinetics* m_InterfaceKinetics;

    //! Pointer to the transport mechanism
    Transport* m_transport;

    //! Number of phases in the kinetics manager
    size_t m_NumKinPhases;

    //! Number of phases in the PhaseList
    size_t m_NumPLPhases;

    //! Number of species in the kinetics manager
    size_t m_NumKinSpecies;

    //! List of ThermoPhase pointers   
    std::vector<ThermoPhase*>tpList;

    //! Mapping between the phase order in the InterfaceKinetics or Kinetics object and the overall phase order
    //! in the PhaseList object
    /*!
     *  Note in the PhaseList object, surface phases are listed after the volume phases.
     *  Length is the number of phases in the InterfaceKinetics object.
     *  Value is the id of the phase in the PhaseList object.
     *
     *  kinOrder[kph] = iph;
     *        kph = phase index in the InterfaceKinetics object
     *        iph = phase index in the PhaseList object
     */
    std::vector<size_t> kinOrder;

    //! Vector of the indecises of each phase in the ElectrodeKinetics_PL object
    //! given the index within the PhaseList object
    /*!
     *       jph = PLtoKinPhaseIndex_[iph];
     *
     *          iph refers to the index of the phase in the PhaseList object 
     *          jph refers to the index of the phase in the heterogeneous kinetics object
     *
     *  Length = number of phases in the PhaseList
     *
     *  A value of -1 or npos in this slot means that the phase doesn't participate in the
     *  current ElectrodeKinetics_intoPL object
     */
    std::vector<size_t> PLtoKinPhaseIndex_;

    //! Vector of the indexes of each species in the ElectrodeKinetics object
    //! given the index within the PhaseList object
    /*!
     *       jsp = PLtoKinSpeciesIndex_[isp];
     *
     *          isp refers to the index of the species in the PhaseList object
     *          jsp refers to the index of the species in the heterogeneous kinetics object
     *
     *  Length = number of species in the PhaseList object
     *
     *  A value of npos in this slot means that the species doesn't participate in the
     *  current ElectrodeKinetics_intoPL object
     */
    std::vector<size_t> PLtoKinSpeciesIndex_;

    //! Index mapping kinetics species index to the PhaseList species index.
    /*!
     *   Length is the number of species in the kinetics species list
     *   Length:   m_NumKinSpecies;
     */
    std::vector<size_t> KintoPLSpeciesIndex_;

    //! Phase Index within the PhaseList which is the owning phase of the kinetics object
    size_t m_iphGlobKin;

    //! String representing the transport model
    std::string transportModel;

    //! Do surface kinetics
    bool m_DoSurfKinetics;

    //! do homogeneous kinetics
    bool m_DoHomogKinetics;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
