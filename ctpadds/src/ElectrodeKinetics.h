/**
 * @file ElectrodeKinetics.h
 *
 */


/*
 *  $Id: ElectrodeKinetics.h 508 2013-01-07 22:54:04Z hkmoffa $
 */

/*
 * Copywrite 2007 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef CT_ELECTRODEKINETICS_H
#define CT_ELECTRODEKINETICS_H

#include <fstream>
#include <math.h>
#include <map>
#include <stdlib.h>

#include "cantera/thermo/mix_defs.h"
#include "cantera/kinetics/Kinetics.h"

#include "cantera/base/utilities.h"
#include "cantera/kinetics/RateCoeffMgr.h"
#include "cantera/kinetics/ReactionStoichMgr.h"
#include "cantera/kinetics/InterfaceKinetics.h"

namespace Cantera {

  // forward references

  class ReactionData;
  class InterfaceKineticsData;
  class ThermoPhase;
  class SurfPhase;
  class ImplicitSurfChem;

  /*!
   * Integer ID, add to mix_defs.h
   */
  const int cElectrodeKinetics = 89; 

    
  //!  A kinetics manager for heterogeneous reaction mechanisms between
  //! surfaces and water involving charge transfer
  /*!
   *  The reactions are assumed to occur at a 2D interface between two
   *  3D phases.
   */
  class ElectrodeKinetics : public InterfaceKinetics {

  public:

    //! Constructor
    /*!
     * Took out the option to initialize with a ThermoPhase Object.
     */
    ElectrodeKinetics();

    //! Destructor.
    virtual ~ElectrodeKinetics();

    //! Return the ID type of the kinetics mechanism
    virtual int ID() const { return cElectrodeKinetics; }

    //! Return the type of the kinetics mechanism
    virtual int type() const { return cElectrodeKinetics; }

    //!  Add a phase to the kinetics manager object.
    /*!
     * This must be done before the function init() is called or
     * before any reactions are input.
     * The following fields are updated:
     *  m_start -> vector of integers, containing the
     *             starting position of the species for
     *             each phase in the kinetics mechanism.
     *  m_surfphase -> index of the surface phase.
     *  m_thermo -> vector of pointers to ThermoPhase phases
     *              that participate in the kinetics
     *              mechanism.
     *  m_phaseindex -> map containing the std::string id of each
     *              ThermoPhase phase as a key and the
     *              index of the phase within the kinetics
     *              manager object as the value.
     *
     * @param thermo    Reference to the ThermoPhase to be added.
     */
    virtual void addPhase(thermo_t& thermo);

    ///
    ///  @name Reaction Rates Of Progress
    ///
    //@{

  
 

    //@}
    /**
     * @name Species Production Rates
     */
    //@{


    //@}
    /**
     * @name Reaction Mechanism Informational Query Routines
     */
    //@{


    //@}
    /**
     * @name Reaction Mechanism Construction
     */
    //@{

   
    //! Import a reaction mechanism for a phase or an interface. 
    //! This is modified for ElectrodeKinetics only pending a generalized treatment.
    /*!
     *
     * @param phase This is an xml node containing a description
     *              of a phase. Within the phase is a XML element
     *              called reactionArray containing the location
     *              of the description of the reactions that make
     *              up the kinetics object. 
     *              Also within the phase is an XML element called
     *              phaseArray containing a listing of other phases
     *              that participate in the kinetics mechanism.
     *
     * @param th    This is a list of ThermoPhase pointers which must
     *              include all of
     *              the phases that participate in the kinetics
     *              operator. All of the phases must have already
     *              been initialized and formed within Cantera.
     *              However, their pointers should not have been
     *              added to the Kinetics object; this addition
     *              is carried out here. Additional phases may
     *              be include; these have no effect.
     *
     * @return
     *         Return true if successful, false otherwise.
     */
    bool importKinetics(const XML_Node& phase, std::vector<ThermoPhase*> th);

    //! Read and Install Reactions into Kinetics Object
    /*!
     *  Take information from the XML tree, p, about reactions
     *  and install them into the kinetics object, kin. 
     *  default_phase is the default phase to assume when
     *  looking up species.
     *
     *  At this point, p usually refers to the phase xml element.
     *  One of the children of this element is reactionArray,
     *  the element which determines where in the xml file to
     *  look up the reaction rate data pertaining to the phase.
     *
     *  @param p   XML_Node phase element
     *  @param owning_phase String name of the owning phase
     *  @param check_for_duplicates  Boolean indicating whether
     *        an operation should be done to check for duplicates.
     *
     *  @return 
     *    If reaction instantiation goes correctly, return true.
     *    If there is a problem, return false.
     */
    bool installReactionArrays(const XML_Node& p,  
			       std::string default_phase,
			       bool check_for_duplicates);

    /**
     * Prepare the class for the addition of reactions. This function
     * must be called after instantiation of the class, but before
     * any reactions are actually added to the mechanism.
     * This function calculates m_kk the number of species in all
     * phases participating in the reaction mechanism. We don't know
     * m_kk previously, before all phases have been added.
     */
    virtual void init();
  
    
    //! Add a single reaction to the mechanism. 
    /*!
     *   This routine
     *   must be called after init() and before finalize().
     *   This function branches on the types of reactions allowed
     *   by the interfaceKinetics manager in order to install
     *   the reaction correctly in the manager.
     *   The manager allows the following reaction types
     *         - Elementary
     *         - Surface
     *         - Global  
     *   There is no difference between elementary and surface 
     *   reactions.
     *
     * @param r Reference to the ReactionData data type
     */
    virtual void addReaction(ReactionData& r);
    
    //! Finish adding reactions and prepare for use.
    /*!
     * This function
     * must be called after all reactions are entered into the mechanism
     * and before the mechanism is used to calculate reaction rates.
     */
    virtual void finalize();

    virtual bool ready() const;



  
  protected:

  
    void addElementaryReaction(ReactionData& r);
    void addGlobalReaction(const ReactionData& r);
    void installReagents(const ReactionData& r);


    //! Apply modifications for the forward reaction rate for interfacial charge transfer reactions
    /*!
     * For reactions that transfer charge across a potential difference,
     * the activation energies are modified by the potential difference.
     * (see, for example, ...). This method applies this correction.
     *
     * @param kfwd  Vector of forward reaction rate constants on which to have
     *              the voltage correction applied
     */
    void applyVoltageKfwdCorrection(doublereal* const kfwd);

  };
}

#endif
