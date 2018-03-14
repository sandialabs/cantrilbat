/**
 * @file ReactingSurDomain.h
 *  Declarations for the ElectrodeKinetics object that does handles interactions with the PhaseList object
 *  (see class \link Zuzax::ReactingSurDomain ReactingSurDomain\endlink).
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef REACTINGSURDOMAIN_H
#define REACTINGSURDOMAIN_H

#include "cantera/kinetics/ElectrodeKinetics.h"
#include "cantera/multiphase/PhaseList.h"

#include "RSD_OCVmodel.h"

#include <string>
#include <vector>

class RxnMolChange;
//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{

// Forward declartion for structure. This is fully defined in Electrode_input.h
struct OCV_Override_input;

//==================================================================================================================================
//!  ReactingSurDomain is a class of reaction that combines the PhaseList information with the Zuzax ElectrodeKinetics class
/*!
 *      This is an inheritance class for the ElectrodeKinetics class of Zuzax.
 *      Essentially this class takes care of the bookkeeping between Zuzax and the PhaseList class
 *
 *      The class also implements the OCV override.
 *
 *      Usually the Interface class calculates temporary vectors over and over again for thermodynamic information.
 *      Here we keep vectors for phase gibbs free energies, enthalpies and entropies. Then we override the 
 *      these entries for the particular species that is designitated to receive the OCV override information.
 *
 *      The standard state thermodynamic functions for the delta of reactions are not overridden. They 
 *      still refer to the unchanged thermodynamics values.
 *
 *      todo: Working on splitting this class into 2. The first one will be located in Zuzax/kinetics: ElectrodeKinetics_into_PL
 *
 *      todo: This class is missing some obvious member functions, that are probably carried out manually within the Electrode
 *            object
 */
class ReactingSurDomain : public ZZCantera::ElectrodeKinetics
{
public:

    //! Default constructor
    ReactingSurDomain();

    //! Constructor based on a PhaseList object
    /*!
     *  @param[in]         pl                Pointer to the phase list object 
     *  @param[in]         iskin             integer number for the  interfacial kinetics object within the PhaseList (which may contain
     *                                       multiple surfaces each with its own associated interfacial kinetics object, 
     *                                       that this object will be inherited from. 
     */
    ReactingSurDomain(ZZCantera::PhaseList* pl, size_t iskin);

    //! Copy Constructor for the %Kinetics object.
    /*!
     *  Currently, this is not fully implemented. If called it will  throw an exception.
     *
     *  @param[in]       right                Reference to %Kinetics object to be copied into the current one.
     */
    ReactingSurDomain(const ReactingSurDomain& right);

    //! Assignment operator
    /*!
     *
     *  @param[in]       right                Reference to %Kinetics object to be copied into the current one.
     *
     *  @return                               Returns a reference to the current object
     */
    ReactingSurDomain& operator=(const ReactingSurDomain& right);

    //! Default destructor
    virtual ~ReactingSurDomain();

    //!  Duplication routine for objects which inherit from Kinetics
    /*!
     *   This virtual routine can be used to duplicate %Kinetics objects
     *   inherited from %Kinetics even if the application only has a pointer to %Kinetics to work with.
     *
     *   These routines are basically wrappers around the derived copy constructor. Here we provide an
     *   opportunity to associate the new Kinetics duplicate with a different set of ThermoPhase objects.
     *   This means that we provide a vector of ThermoPhase pointers as input to this routine.
     *
     *   @param[in]     tpVector              Vector of shallow pointers to ThermoPhase objects. This is the
     *                                        m_thermo vector within this object.
     *
     *   @return                              Returns a pointer to the duplicate object.
     */
    virtual Kinetics* duplMyselfAsKinetics(const std::vector<thermo_t*>& tpVector) const override;

    //!  Imports and initializes the surface kinetics for a kinetics object given a PhaseList object that is already set up
    /*!
     *   This routine initializes a ReactingSurDomain object and its underlying interfacial kinetics
     *   given a PhaseList object. It uses the XML data written into the PhaseList object to go find
     *   and initialize the surface kinetics. On input the surface phase containing the surface kinetics
     *   must be identified via the parameter list. It then reads the XML data for the surface phase.
     *   It then finds the ThermoPhase objects within the PhaseList associated with the surface kinetics.
     *
     *   The routine finds all of the ThermoPhase objects that are part of the Kinetics object by
     *   querying the XML node, phaseArray. It then finds these ThermoPhase objects within the PhaseList
     *   objects to see if the kinetics object can be successfully formulated.
     *   It then creates the kinetics object by calling Zuzax's importPhase() routine.
     *   
     *   Then, other initializations are carried out such as formulating the kinetics species list
     *   and its indexing into the Phaselist species list.
     *
     *   Right now it is an error for there not to be a kinetics mechanism associated with the argument surface phase.
     *
     *   @param[in]     pl                         Fully formed pointer to the PhaseList object that will be associated
     *                                             with this object.
     *
     *   @param[in]     iphSurKin                  The surface phase index of the surface phase within the PhaseList that has the 
     *                                             surface kinetics associated with it. 
     *
     *   @return                                   Returns true upon proper instanteation of the kinetics. Returns false
     *                                             if there was a problem.
     */
    bool importFromPL(ZZCantera::PhaseList* const pl, size_t iphSurKin);

    //! Reassign the internal PhaseList point to a new PhaseList object
    /*!
     *  This is used in the copy constructor.
     *  Right now, we reassign the PhaseList without much checking. However, we should add more checking.
     *
     *  @param[in]           pl_ptr              New PhaseList pointer. 
     */
    void reassignPhaseList(PhaseList* pl_ptr);

    //! Redo the initialization
    /*!
     *  This redoes the index initializations between the kinetics object and the PhaseList object.
     *  This must be called if the PhaseList class has added or subtracted phases,
     *  or if a ThermoPhase class has added or subtracted species.
     */
    virtual void reinitializeIndexing();

    //! Returns a reference to the PhaseList object
    /*!
     *  @return                                  Returns a reference to the PhaseList object
     */
    const PhaseList& pl() {
       return *m_pl;
    }

    //! Routine to be called after all species and phases have been defined for the object
    /*!
     *  This initializes based on number of species
     */
    virtual void init() override;

    //! Routine to be called after all reactions have been defined for the object.
    /*!
     *  This routine initializes vectors based on the number of reactions.
     */
    virtual void finalize() override;

    //
    //------------------------------------------------------------------------------------------------------------------------------
    //! @name               Conversion of Phase and Species Indecise Between PhaseList and Kinetics Objects
    //------------------------------------------------------------------------------------------------------------------------------
    //@{

    //! Returns the global phase index of the phase in the PhaseList given the index in the Kinetics object
    /*!
     *  @param[in]           iphKin              Phase index within the kinetics object
     *
     *  @return                                  Returns the global phase index in the PhaseList
     */
    size_t globalPhaseIndex_fromKP(size_t iphKin) const;

    //! Returns the global phase index within the PhaseList corresponding to the phase which owns the kinetic object
    /*!
     *  @return                                  Returns the index within PhaseList that owns the Kinetics object
     */
    size_t globalReactionPhaseIndex() const {
        return m_iphGlobKin;
    }

    //! Returns the kinetics phase index of the phase given the global index in the PhaseList
    /*!
     *  If not present returns npos
     *  @param[in]           iphGlob             Global phase index within the Phaselist object
     *
     *  @return                                  Returns the kinetics phase index or npos if not present
     */
    size_t kineticsPhaseIndex_fromPLP(size_t iphGlob) const;

    //! Returns the offset from the start of the PhaseList species vector for the argument Kinetics phase index
    /*!
     *  @param[in]           iphKin              Phase index within the Kinetics object
     *
     *  @return                                  Returns the offset within a PhaseList species vector from the start of the species
     *                                           belonging to the Kinetics phase listed in the parameter
     */
    size_t globalSpeciesStart_fromKP(size_t iphKin) const;

    //! Returns the offset from the start of the Kinetics object species vector for the species in the argument Phaselist phase index
    /*!
     *  Phases which aren't in the Kinetics object will return an npos value.
     *
     *  @param[in]           iphGlob             Phase index within the PhaseList object
     *
     *  @return                                  Returns the offset within a Kinetics species vector from the start of the species
     *                                           belonging to the PhaseList phase listed in the parameter. 
     *                                           If that phase doesn't exist, this returns npos.
     */
    size_t kineticsSpeciesStart_fromPLP(size_t iphGlob) const;

    //@}

    //! Returns a reference to the calculated production rates of species from this interfacial Kinetics class
    /*!
     *  This routine calls thet getNetProductionRate() function and then returns a reference to the result.
     *
     *  @return                               Vector of length m_kk containing the species net
     *                                        production rates (kmol s-1 m-2)
     */
    const std::vector<doublevalue>& calcNetSurfaceProductionRateDensities();

    //! Returns a reference to the calculated limited production rates of species from this interfacial Kinetics class
    /*!
     *   This routine first calls the limitROP() function to limit the rates near a phase boundary
     *   This routine calls thet getNetProductionRate function and then returns a reference to the result.
     *
     *  nMoles represents a scaled total moles system, where the sum of the total moles of all species in all phases has been scaled to 
     *  a value roughly equal to one before the calculation has started. 
     *  Therefore, absolute cutoffs of phase mole numbers can be compared to the value of 1 kmol.
     *
     *  @param[in]          nMoles            nMoles is a vector of species mole numbers for all species in the PhaseList
     *                                                 Length: global number of species in PhaseList
     *                                                 Units: kmol (scaled)
     *
     *  @return                               Vector of length m_kk containing the species net
     *                                        production rates (kmol s-1 m-2)
     */
    const std::vector<doublevalue>& calcNetLimitedSurfaceProductionRateDensities(const doublevalue* const nMoles);

    //!  Apply smooth limiter to ROP
    /*!
     *   Given the distance towards getting rid of a phase, this routine will limit the rate of progress
     *   vectors within the Electrode object.
     * 
     *   @param[in]          nMoles              Scaled vector of species mole numbers for all species in the PhaseList.
     *                                           The scale requires that the total moles of solid species be roughly
     *                                           equal to 1 kmol.
     *                                             Length: number of species in PhaseList
     *                                             Units: kmol
     *
     *         Seems to basically satisfy the needs of turning interfacial kinetics into homogeneous kinetics.
     *         There is a basic exponential decay algorithm created out of the interfacial kinetics, which is
     *         not exponential decay.
     *
     *  HKM -> several problems. One is that the surface area doesn't enter into this limiter. Therefore, the
     *         limiter doesn't have the correct units, and will fail as surface area scales to different values
     *         than what it was set up to be. That way you could also talk about the algorithm setting a minimum "time
     *         constant" for phase removal.
     *
     *         Second is that rate turns off for phase_moles = 1.0E-20. There's no reason for that as this is straight
     *         exponential decay. The 1.0E-20 value has to be reconciled with other small numbers in the algorithm
     *         to see if there isn't a conflict.
     *
     *         Third, the algorithm only works for phase_moles > 0.0. Phase pop problems typically have times when
     *         phase_moles < 0.0. What happens in these cases. Also non-max principle algorithms will have problems.
     *
     *         Fourth, algorithm should maybe not be based on phase_moles, but individual moles. One really wants
     *         all moles to stay positive, not just the phase moles.
     */
    void limitROP(const doublevalue* const nMoles);

    //! Returns a constant reference to the vector of reaction rates of progress
    /*!
     *  This is just a wrapper around updateROP() and then it returns the internal vector m_ropnet;
     *
     *  @return                                 Returns a const reference to  the vector of reaction rates of progress.
     *                                          The units are kmol m-2 s-1.
     */
    const std::vector<doublevalue>& calcNetSurfaceROP();

    //! Returns a reference to the calculated creation rates of species
    /*!
     *   This routine calls thet getCreationRate function
     *   and then returns a reference to the result.
     *
     * @return Vector of length m_kk containing the species creation rates
     */
    const std::vector<doublevalue>& calcSurfaceCreationRateDensities();

    //! Returns a reference to the calculated destruction rates of species
    /*!
     *   This routine calls thet getDestructionRate function
     *   and then returns a reference to the result.
     *
     *  @return                                  Vector of length m_kk containing the species destruction rates
     */
    const std::vector<doublevalue>& calcSurfaceDestructionRateDensities();

    //!  Update the standard state chemical potentials and species equilibrium constant entries
    /*!
     *  Virtual because it is overwritten when dealing with experimental open circuit voltage overrides
     */
    virtual void updateMu0() override;

    //! Get the net current for the set of reactions on this surface in amps m-2.
    /*!
     *  Get the net current.  We use the value of kElectronIndex_ to identify the index for the electron species.
     *  Then the net electron production is used for the net current across the interface. Note a positive electron
     *  production means that positive charge is being put from the solid into the solution. We identify this
     *  as a positive current. Therefore, under normal operations an anode has positive current and a cathode 
     *  has a negative current.
     * 
     *  @param[out]      currentDensityRxn          Returns a vector containing the net current from all
     *                                              reactions in the mechanism. On input this must be of length
     *                                              nReactions(). On output the units are amps m-2.
     *
     *  @return                                     Returns the net current in amps m-2 from all reactions
     *                                              at the current conditions.
     */
    double getCurrentDensityRxn(doublevalue* const currentDensityRxn = 0);

    //! Get the net current density for the set of reactions on this surface in amps m-2.
    /*!
     *  Get the net current.  We use the value of kElectronIndex_ to identify the index for the electron species.
     *  Then the net electron production is used for the net current across the interface. Note a positive electron
     *  production means that positive charge is being put from the solid into the solution. We identify this
     *  as a positive current. Therefore, under normal operations an anode has positive current and a cathode 
     *  has a negative current.
     * 
     *  This function calls limitROP() to limit the ROP near end of phases.
     *
     *  nMoles represents a scaled total moles system, where the sum of the total moles of all species in all phases has been scaled to 
     *  a value roughly equal to one. Therefore, absolute cutoffs of phase mole numbers can be compared to the value of 1 kmol
     *  for understanding what cutoffs values mean. 
     * 
     *   @param[in]          nMoles                  nMoles is a vector of species mole numbers for all species in the PhaseList
     *                                                 Length: global number of species in PhaseList
     *                                                 Units: kmol (scaled)
     *
     *   @return                                     Returns the net limited current in amps m-2 from all reactions
     *                                               at the current conditions.
     */
    double getLimitedCurrentDensityRxn(const doublevalue* const nMoles);

    //! Add an open circuit voltage override feature to the current reacting surface
    /*!
     *     The open circuit voltage will replace the deltaG calculation for this reaction.
     *     It will do this by replacing the thermo functions for one of the species in the mechanism.
     *
     *      @param ocv_ptr    The open circuit voltage feature is described by this pointer.
     *
     */
    void addOCVoverride(OCV_Override_input *ocv_ptr);

    //! Calculate the effective chemical potential of the replaced species
    /*!
     *  We calculate the effects of the OCV override on the storred thermodynamics of the species
     *  We modify the following storage variables
     *
     *       After Transformation            Before Transformation
     *          m_mu
     *          m_mu0 
     *          m_GibbsOCV_rspec             
     *
     *                                       deltaGRxnOCV_Before_
     */
    void deriveEffectiveChemPot();

    //! Calculate the effective thermodynamic variables of all of the species due to the OCV override
    /*!
     *  We calculate the effects of the OCV override on the storred thermodynamics of the species
     *  We modify the following storage variables
     *
     *       After Transformation            Before Transformation
     *          m_mu
     *          m_mu0 
     *          m_GibbsOCV_rspec           
     *          m_Entropies_rspec            m_Entropies_Before_rspec
     *          m_Enthalpies_rspec           m_Enthalpies_Before_rspec
     *
     *                                       deltaGRxn_Before_
     *                                       deltaHRxn_Before_
     *                                       deltaSRxn_Before_
     *
     *                                       deltaGRxnOCV_Before_
     *     
     */
    void deriveEffectiveThermo();

    //!  Get the vector of deltaG values for all reactions defined in the kinetics object
    /*!
     *   (Virtual from Kinetics)
     *
     *   This routine provides an override to the normal calculation of deltaG, when the thermodynamics
     *   is modified by a specification of the open circuit potential.
     *
     *   @param[out]        deltaG      Vector of deltaG values. Must be at least of length equal to the number of reactions.
     */
    virtual void getDeltaGibbs(doublevalue* const deltaG) override;

    //!  Get the vector of deltaG values for all reactions defined in the kinetics object
    /*!
     *
     *   This routine provides the normal calculation of deltaG in all cases. It ignores the
     *   possible presence of an override.
     *
     *   @param[out]        deltaG      Vector of deltaG values. Must be at least of length equal to the number of reactions.
     */
    void getDeltaGibbs_Before(doublevalue* const deltaG = 0);

    //! Return the vector of values for the electrochemical free energy change of reaction.
    /*!
     *  (virtual from Kinetics)
     *  These values depend upon the concentration of the solution and the  voltage of the phases
     *  involved with each reaction that has a charged species as a participant.
     *
     *  units = J kmol-1
     *
     *  @param[out]          deltaM              Vector of Reaction delta electrochemical potentials
     *                                           If 0, this updates the internally stored values only, in m_deltaM
     */
    virtual void getDeltaElectrochemPotentials(doublevalue* const deltaM) override;

    //!  Get the vector of deltaH values for all reactions defined in the kinetics object
    /*!
     *   (Virtual from Kinetics.h)
     *
     *   This routine provides an override to the normal calculation of deltaH, when the thermodynamics
     *   is modified by a specification of the open circuit potential and the derivative of the
     *   OCV wrt temperature.
     *
     *   @param[out]        deltaH      Vector of deltaH values. Must be at least of length equal
     *                                  to the number of reactions
     */
    virtual void getDeltaEnthalpy(doublevalue* const deltaH) override;

    //!  Get the vector of deltaH values for all reactions defined in the kinetics object
    /*!
     *   This routine provides the normal calculation of deltaH, before any override modification.
     *
     *   @param[out]        deltaH      Vector of deltaH values. Must be at least of length equal
     *                                  to the number of reactions, m_ii
     */
    void getDeltaEnthalpy_Before(doublevalue* const deltaH);

    //! This gets the deltaG for each reaction in the mechanism, but using the standard state
    //! chemical potential for the electrolyte.
    /*!
     *          @param  deltaG_special  DeltaG for each reaction using standard state chemical
     *                                  potentials for the electrolyte. 
     *                     length = nReactions(), J/kmol
     */
    void getDeltaGibbs_electrolyteSS(doublevalue* const deltaG_special);

    //!  Get the vector of deltaS values for all reactions defined in the kinetics object
    /*!
     *   (Virtual from Kinetics.h)
     *
     *   This routine provides an override to the normal calculation of deltaS, when the thermodynamics
     *   is modified by a specification of the open circuit potential and the derivative of the
     *   OCV wrt temperature.
     *
     *   @param[out]        deltaS      Vector of deltaS values. Must be at least of length equal
     *                                  to the number of reactions. Units are J kmol-1 K-1.
     */
    virtual void getDeltaEntropy(doublevalue* const deltaS) override;

    //!  Get the vector of deltaS values for all reactions defined in the kinetics object before OCV override
    /*!
     *   This routine does a normal calculation of deltaS, before any override of thermodynamics is carried out.
     *
     *   @param[out]        deltaS      Vector of deltaS values. Must be at least of length equal
     *                                  to the number of reactions. Units are J kmol-1 K-1.
     */
    void getDeltaEntropy_Before(doublevalue* const deltaS);

    //!  Return the vector of values for the reaction standard state gibbs free energy change.  These values don't depend upon
    //!  the concentration of the solution.
    /*!   
     *  (virtual from Kinetics.h)
     *  units = J kmol-1
     *
     * @param[out]      deltaG            Output vector of ss deltaG's for reactions Length: m_ii.
     */
    virtual void getDeltaSSGibbs(doublevalue* const deltaG) override;

    //!  Return the vector of values for the change in the standard
    //! state enthalpies of reaction.  These values don't depend
    //! upon the concentration of the solution.
    /*!
     *  (virtual from Kinetics.h)
     *  units = J kmol-1
     *
     * @param[out]     deltaH              Output vector of ss deltaH's for reactions Length: m_ii.
     */
    virtual void getDeltaSSEnthalpy(doublevalue* const deltaH) override;
    
    //! Return the vector of values for the change in the standard state entropies for each reaction.  These values don't
    //! depend upon the concentration of the solution.
    /*! 
     *  (virtual from Kinetics.h)
     *  units = J kmol-1 Kelvin-1
     *
     *     @param[out]    deltaS           Output vector of ss deltaS's for reactions Length: m_ii.
     */
    virtual void getDeltaSSEntropy(doublevalue* const deltaS) override;

    //!  Get the OCV thermodynamic functions offsets for the species that is replaced when carrying out
    //!  an OCV override step
    /*!
     *      @param[out]   deltaG_species   Change in the value of chemical potential
     *      @param[out]   deltaH_species   Change in the value of the enthalpy
     *      @param[out]   deltaS_species   Change in the value of the entropy
     */
    void getOCVThermoOffsets_ReplacedSpecies(double& deltaG_species, double& deltaH_species, double& deltaS_species);

public:
    //! Declare a printing routine as a friend to this class
    /*!
     *   @param[in]        s              Reference to the ostream that will be used for the printing
     *   @param[in]        rsd            Reference to the ReactingSurDomain whose values will be printed
     *
     *   @return                          Returns a reference to the input ostream, as required for chaining these commands together.
     */
    friend std::ostream& operator<<(std::ostream& s, ReactingSurDomain& rsd);

    //
    //   -----------------------------------   DATA --------------------------------------------------------------------------
    //

    //! Mapping between the phase order in the InterfaceKinetics object and the overall order in the PhaseList object
    /*!
     *  Note in the phase list object, surface phases are listed after volume phases, and edge phases are listed last.
     *  Length is the number of phases in the InterfaceKinetics object. Value is the index of the phase in the PhaseList object.
     *
     *  KinToPL_PhaseIndex_[kph] = iph;
     *        kph = phase index in the InterfaceKinetics object
     *        iph = phase index in the PhaseList object
     *
     *  Length:   Number of phases in the Kinetics object 
     *  Indexing: iphKin =  index of the kinetics phase within the Kinetics object
     */
    std::vector<size_t> KinToPL_PhaseIndex_;

    //! Vector of the indexes of each phase in the ReactionSurfaceDomain object
    //! given the index within the PhaseList object
    /*!
     *       jph = PLToKin_PhaseIndex_[iph];
     *
     *          iph refers to the index of the phase in the Electrode_Model object 
     *          jph refers to the index of the phase in the heterogeneous kinetics object
     *
     *
     *  A value of -1 or npos in this slot means that the phase doesn't participate in the
     *  current ReactingSurDomain object
     *
     *  Length:   m_NumTotPhases =  number of phases in the PhaseList object 
     *  Indexing: iphGlob = global phase indexing within the PhaseList object
     */
    std::vector<size_t> PLToKin_PhaseIndex_;

    //! Index mapping kinetics species start index to the PhaseList species start index.
    /*!
     *  Value is the offset index for that kinetics phase within the species PhaseList vector
     *
     *  Length:   Number of phases in the Kinetics object 
     *  Indexing: iphKin =  index of the kinetics phase within the Kinetics object
     */
    std::vector<size_t> KinToPL_SpeciesStartIndex_;

    //! Index mapping phaselist species start index to the kinetics species start index.
    /*!
     *  Value is the offset index for that Phaselist phase within the species Kinetics vector
     *
     *  Phases which are missing in the kinetics object have a npos value in their entry
     *
     *  Length:   m_NumTotPhases =  number of phases in the PhaseList object 
     *  Indexing: iphGlob = global phase indexing within the PhaseList object
     */
    std::vector<size_t> PLToKin_SpeciesStartIndex_;

    //! Vector of the indexes of each species in the ReactionSurDomain object given the index within the PhaseList object
    /*!
     *       ksp = PLToKin_SpeciesIndex_[kGlob];
     *
     *          kGlob refers to the index of the global species in the PhaseList object
     *          ksp refers to the index of the species in the heterogeneous kinetics object
     *
     *
     *  A value of npos in this slot means that the species doesn't participate in the
     *  current ReactingSurDomain kinetics object
     *
     *  Length:   m_pl->m_NumTotSpecies =  number of species in the PhaseList object
     *  Indexing: kGlob  = global species indexing within the PhaseList object
     */
    std::vector<size_t> PLToKin_SpeciesIndex_;

    //! Index mapping kinetics species index to the PhaseList global species index.
    /*!
     *   kGlob = KinToPl_SpeciesIndex_[ksp];
     *
     *   Length:    m_NumKinSpecies = number of species in the kinetics species list
     *   Indexing:  kKin = kinetic species indexing
     */
    std::vector<size_t> KinToPL_SpeciesIndex_;

    //! Global phase index of the phase in the PhaseList object that has the kinetics
    //! object for this reacting surface.
    size_t m_iphGlobKin;

    //! List of id()s that constitute the ThermoPhases needed for the kinetics object
    /*!
     *   Currently this is needed to fix up the shallow pointer copy operation. 
     *   The list is needed to carry out shallow pointer assignments when an Electrode object is copied.
     */
    std::vector<std::string> tpList_IDs_;

    //! If there is a surface kinetics mechanism associated with this object, this is true. 
    /*!
     *   The index of the surface phase is kept in the variable, iphaseKin_ 
     */
    bool m_DoSurfKinetics;

    //! Vector that will expose the species production rates for this kinetics object
    /*!
     *  Length:   m_NumKinSpecies = Number of species in the kinetics vector
     *  Units:    kmol m-2 s-1.
     *  Indexing: kKin = Kinetics species indexing
     */
    std::vector<double> speciesProductionRates_;

    //! Vector that will expose the species creation rates for this kinetics object
    /*!
     *  Length:   m_NumKinSpecies = Number of species in the kinetics vector
     *  Units:    kmol m-2 s-1.
     *  Indexing: kKin = Kinetics species indexing
     */
    std::vector<double> speciesCreationRates_;

    //! Vector that will expose the species destruction rates for this kinetics object
    /*!
     *  Length:   m_NumKinSpecies = Number of species in the kinetics vector
     *  Units:    kmol m-2 s-1.
     *  Indexing: kKin = Kinetics species indexing
     */
    std::vector<double> speciesDestructionRates_;

    //! Internal Vector of deltaG of reaction for all reactions
    /*!
     *  This is used to store the DeltaG of reaction before modification due to OCV override.
     *  The length is equal to nreactions(), m_ii. The units are Joules kmol-1.
     *
     *  Length:   m_ii = Number of reactions
     *  Units:    Joules kmol-1
     *  Indexing: Reaction number
     */
    std::vector<double> deltaGRxn_Before_;

    //! Internal Vector of deltaH of reaction for all reactions
    /*!
     *  This is used to store the DeltaH of reaction before modification due to OCV override.
     *
     *  Length:   m_ii = Number of reactions
     *  Units:    Joules kmol-1
     *  Indexing: Reaction number
     */
    std::vector<double> deltaHRxn_Before_;

    //! Internal Vector of deltaS of reaction for all reactions
    /*!
     *  This is used to store the DeltaS of reaction before modification due to OCV override.
     *
     *  Length:   m_ii = Number of reactions
     *  Units:    Joules kmol-1
     *  Indexing: Reaction number
     */
    std::vector<double> deltaSRxn_Before_;

protected:
    //! Pointer to the PhaseList object that contains the ThermoPhase objects.
    /*!
     *  This object doesn't own this. However, it uses this heavily. It is a shallow pointer.
     */
    ZZCantera::PhaseList* m_pl;

    //! This is true if the PhaseList and the KinSpecies list are exactly the same thing
    /*!
     *  This has advantages in terms of the speed for posting results for source terms
     */
    bool m_OneToOne;

    //! This is true if the KinSpecies list is contiguous within the PhaseList list
    /*!
     *  This will be true if the phases within the Kinetics object are ordered in the same as in 
     *  the PhaseList object. There may be additional phases in the PhaseList object in the beginning or the end of the object.
     *  However, the kinetic species vector may be considered to be a subset of the PhaseList vector.
     */
    bool m_IsContiguous;

public:

    // -----------------------------------------------------------------------------------------------------------------------------
    //           DATA for the limiting ROP model
    // -----------------------------------------------------------------------------------------------------------------------------

    //! Vector of Limited Rates of Progress of the reactions
    /*!
     *  Length is the number of reactions, n_ii. The units are kmol m-2 s-1.
     */
    std::vector<double> limitedROP_;

    //! Vector of Limited Rates of Progress of the species production rates
    /*!
     *  Length is the number of species. The units are kmol m-2 s-1.
     */
    std::vector<double> limitedNetProductionRates_;

    // -----------------------------------------------------------------------------------------------------------------------------
    //           DATA FOR THE OCV override mode
    // -----------------------------------------------------------------------------------------------------------------------------

    //!  Pointer to an OCV_Override_input object which is used to override the thermodynamics of an electrode
    //!  reaction given input
    /*!
     *      The default value is NULL, indicating that there isn't any override.
     */
    OCV_Override_input *ocv_ptr_;

    //! Pointer to the OCV model that may be used to override the thermodynamics for the OCV for this interface
    /*!
     *  If there isn't an override, this is set to zero.
     */
    RSD_OCVmodel* OCVmodel_;

    //!  Kinetic species index for the species whose thermodynamics representation is replaced by an experimentally
    //!  determined open circuit voltage expression. If there is none, then this value is npos.
    size_t kReplacedSpeciesRS_;

    //!  Vector of the chemical potential for all species in all phases that participate in the reaction mechanism, modified for OCVoverride
    /*!
     *   We keep a vector of Gibbs energy here over all kinetics species. This vector gets modified from its
     *   strictly thermodynamic origin when an OCV override is done, so that the measured OCVoverride may be matched.
     *   This is a specified combination of regular chemical potentials and standard state chemical potential depending 
     *   on the OCV_Format value.
     *
     *   Length:   m_NumKinSpecies
     *   Units:    Joules/kmol
     *   Indexing: kKin = kinetic species index
     */
    std::vector<double> m_GibbsOCV_rspec;

    //! Internal Vector of deltaG of reaction for all reactions (without the electrolyte mixing entropy term)
    /*!
     *  This is used to store the DeltaG of reaction before modification due to OCV override. Note this isn't
     *  deltaG, necessarily because it is missing the electorlyte mixing entropy Term. ther exact formulation depends
     *  on OCV format.
     *
     *  Length:   m_ii = Number of reactions
     *  Units:    Joules kmol-1
     *  Indexing: Reaction number
     */
    std::vector<double> deltaGRxnOCV_Before_;

    //!  Vector of the enthalpies for all species in all phases that participate in the reaction mechanism, modified for OCVoverride
    /*!
     *   We keep a vector of enthalpies here over all reaction species. This vector gets modified from its
     *   strictly thermodynamic origin when an OCV override is done.
     *
     *   Length: m_NumKinSpecies
     *   Units:  Joules/kmol
     *   Indexing: kinetic species index
     */
    std::vector<double> m_Enthalpies_rspec;

    //!  Vector of the enthalpies for all species in all phases that participate in the reaction mechanism
    /*!
     *   We keep a vector of enthalpies here over all reaction species. This vector does not get modified from its
     *   strictly thermodynamic origin when an OCV override is done.
     *
     *   Length: m_NumKinSpecies
     *   Units:  Joules/kmol
     *   Indexing: kinetic species index
     */
    std::vector<double> m_Enthalpies_Before_rspec;

    //!  Vector of the entropies for all species in all phases that participate in the reaction mechanism, modified for OCVoverride
    /*!
     *   We keep a vector of entropies here over all kinetics species. This vector gets modified from its
     *   strictly thermodynamic origin when an OCV override is done.
     *
     *   Length: m_NumKinSpecies
     *   Units:  Joules/kmol/K
     *   Indexing: kinetic species index
     */
    std::vector<double> m_Entropies_rspec;

    //!  Vector of the entropies for all species in all phases that participate in the reaction mechanism
    /*!
     *   We keep a vector of enthalpies here over all reaction species. This vector does not get modified from its
     *   strictly thermodynamic origin when an OCV override is done.
     *
     *   Length: m_NumKinSpecies
     *   Units:  Joules/kmol/K
     *   Indexing: kinetic species index
     */
    std::vector<double> m_Entropies_Before_rspec;

    //!  Vector of the chemical potential for all species in all phases that participate in the reaction mechanism
    /*!
     *   We keep a vector of gibbs energy here over all reaction species. This vector does not get modified from its
     *   strictly thermodynamic origin when an OCV override is done.
     *
     *   Length:   m_NumKinSpecies
     *   Units:    Joules/kmol
     *   Indexing: kKin = kinetic species index
     */
    std::vector<double> m_Gibbs_Before_rspec;

    //!  Value of deltaG, i.e., the change in the chemical potential, for the particular species that gets modified  
    //!  due to the OCV override 
    /*!
     *   This is the delta value of G at the T and P of solution to generate the correct G value as measured
     *   from data.  This deltaG value is added to the G value of the particular species.
     *
     *   Kinetics species index kReplacedSpeciesRS_
     *   Units:  Joules/kmol
     */
    double deltaG_species_;

    //!  Value of deltaS, i.e., the change in the entropy, for the particular species that gets modified  
    //!  due to the OCV override 
    /*!
     *   This is the delta value of S at the T and P of solution to generate the correct dGdT value as measured
     *   from data.  This deltaS value is added to the S value of the particular species.
     *
     *   Kinetics species index kReplacedSpeciesRS_
     *   Units:  Joules/kmol/K
     */
    double deltaS_species_;

    //!  Value of deltaH, i.e., the change in the enthalpy, for the particular species that gets modified  
    //!  due to the OCV override 
    /*!
     *   This is the value of deltaH at the T and P of solution to generate the correct deltaG value to match
     *   experimental data. This deltaH value is added to the H value of the particular species.
     *
     *   Kinetics species index kReplacedSpeciesRS_
     *   Units:  Joules/kmol
     */
    double deltaH_species_;

protected:

    friend class RSD_OCVmodel;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
