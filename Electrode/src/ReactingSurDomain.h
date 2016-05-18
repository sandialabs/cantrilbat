/**
 * @file ReactingSurDomain.h
 *  Declarations for the ElectrodeKinetics object that does handles interactions with the PhaseList object
 *  (see \ref ExtendedPhaseGroups and class \link Cantera::ReactingSurDomain ReactingSurDomain\endlink).

 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef REACTINGSURDOMAIN_H
#define REACTINGSURDOMAIN_H

#include "cantera/kinetics/ElectrodeKinetics.h"

#include "RSD_OCVmodel.h"
#include "PhaseList.h"

#include <string>
#include <vector>

class RxnMolChange;

namespace Cantera
{

// Forward declartion for structure. This is fully defined in Electrode_input.h
struct OCV_Override_input;

//!  ReactingSurDomain is a class of reaction that combines the PhaseList information with the Cantera ElectrodeKinetics class
/*!
 *      This is an inheritance class for the ElectrodeKinetics class of Cantera.
 *      Essentially this class takes care of the bookkeeping between Cantera and the PhaseList class
 *
 *      The class also implements the OCV override.
 *
 *      Usually the Interface class calculates temporary vectors over and over again for thermodynamic information.
 *      Here we keep vectors for phase gibbs free energies, enthalpies and entropies. Then we override the 
 *      these entries for the particular species that is designitated to receive the OCV override information.
 *
 *      The standard state thermodynamic functions for the delta of reactions are not overridden. They 
 *      still refer to the unchanged thermodynamics values.
 */
class ReactingSurDomain : public Cantera::ElectrodeKinetics
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
    ReactingSurDomain(Cantera::PhaseList* pl, size_t iskin);

    //! Copy Constructor for the %Kinetics object.
    /*!
     *  Currently, this is not fully implemented. If called it will  throw an exception.
     *
     *  @param           right                Reference to %Kinetics object to be copied into the current one.
     */
    ReactingSurDomain(const ReactingSurDomain& right);

    //! Assignment operator
    /*!
     *
     *  @param           right                Reference to %Kinetics object to be copied into the current one.
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
     *   @param         tpVector              Vector of shallow pointers to ThermoPhase objects. This is the
     *                                        m_thermo vector within this object.
     *
     *   @return                              Returns a pointer to the duplicate object.
     */
    virtual Kinetics* duplMyselfAsKinetics(const std::vector<thermo_t*>& tpVector) const;

    //!  Import all the phases from a PhaseList and initialize the surface kinetics for this object
    /*!
     *   This routine initializes a ReactingSurDomain object and its underlying interfacial kinetics
     *   given a PhaseList object. It uses the XML data written into the PhaseList object to go find
     *   and initialize the surface kinetics. On input the surface phase containing the kinetics
     *   must be identified via the parameter list. It then reads the XML data for the surface
     *   phase. It then finds the ThermoPhase objects within the PhaseList associated with the
     *   surface kinetics.
     *
     *   The routine finds all of the ThermoPhase objects that are part of the Kinetics object by
     *   querying the XML node phaseArray. It then finds these ThermoPhase objects within the PhaseList
     *   objects to see if the kinetics object can be successfully formulated.
     *   It then creates the kinetics object by calling Cantera's importPhase() routine.
     *   
     *   Then other initializations are carried out such as formulating the kinetics species list
     *   and its indexing into the Phaselist species list.
     *
     *   Right now it is an error for there not to be a kinetics mechanism.
     *
     *   @param[in]     pl                         Fully formed pointer to the PhaseList object that will be associated
     *                                             with this object.
     *
     *   @param[in]     iskin                      The index of the surface phase within the PhaseList that has the 
     *                                             surface kinetics associated with it.
     *
     *   @return                                   Returns true upon proper instanteation of the kinetics. Returns false
     *                                             if there was a problem.
     */
    bool importFromPL(Cantera::PhaseList* const pl, size_t iskin);

    //! Routine to be called after all species and phases have been defined for the object
    /*!
     *  This initializes based on number of species
     */
    virtual void init();

    //! Routine to be called after all reactions have been defined for the object.
    /*!
     *    This routine initializes vectors based on the number of reactions.
     */
    virtual void finalize();

    //! Returns a reference to the calculated production rates of species from this interfacial Kinetics class
    /*!
     *   This routine calls thet getNetProductionRate function and then returns a reference to the result.
     *
     *  @return                               Vector of length m_kk containing the species net
     *                                        production rates (kmol s-1 m-2)
     */
    const std::vector<double>& calcNetSurfaceProductionRateDensities();

    //! Returns a reference to the calculated limited production rates of species from this interfacial Kinetics class
    /*!
     *   This routine first calls the limitROP() function to limit the rates near a phase boundary
     *   This routine calls thet getNetProductionRate function and then returns a reference to the result.
     *
     *  @param[in]          n                 n is a vector of species mole numbers for all species in the PhaseList
     *
     *  @return                               Vector of length m_kk containing the species net
     *                                        production rates (kmol s-1 m-2)
     */
    const std::vector<double>& calcNetLimitedSurfaceProductionRateDensities(const double* n);

    //!  Apply smooth limiter to ROP
    /*!
     *   Given the distance towards getting rid of a phase, this routine will limit the rate of progress
     *   vectors within the Electrode object.
     * 
     *   @param[in]          n                 n is a vector of species mole numbers for all species in the PhaseList
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
    void limitROP(const double* const n);

    //! Returns a constant reference to the vector of reaction rates of progress
    /*!
     *  This is just a wrapper around updateROP() and then it returns the internal vector m_ropnet;
     *
     *  @return                                 Returns a const reference to  the vector of reaction rates of progress.
     *                                          The units are kmol m-2 s-1.
     */
    const std::vector<double>& calcNetSurfaceROP();

    //! Returns a reference to the calculated creation rates of species
    /*!
     *   This routine calls thet getCreationRate function
     *   and then returns a reference to the result.
     *
     * @return Vector of length m_kk containing the species creation rates
     */
    const std::vector<double>& calcSurfaceCreationRateDensities();

    //! Returns a reference to the calculated destruction rates of species
    /*!
     *   This routine calls thet getDestructionRate function
     *   and then returns a reference to the result.
     *
     * @return Vector of length m_kk containing the species destruction rates
     */
    const std::vector<double>& calcSurfaceDestructionRateDensities();

    //!  Update the standard state chemical potentials and species equilibrium constant entries
    /*!
     *  Virtual because it is overwritten when dealing with experimental open circuit voltage overrides
     */
    virtual void updateMu0();

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
    double getCurrentDensityRxn(double * const currentDensityRxn = 0);

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
     *   @param[in]          n                       n is a vector of species mole numbers for all species in the PhaseList
     *
     *   @return                                     Returns the net limited current in amps m-2 from all reactions
     *                                               at the current conditions.
     */
    double getLimitedCurrentDensityRxn(const double* n);

#ifdef DONOTREMOVE
    //! Get the exchange current density formulation for the current reaction rate
    /*!
     *  Returns the parameterization for a single reaction that defines the exchange
     *  current density formulation for that reaction. Note the exchange current density depends
     *  on the activity concentrations of the reactants and the products.
     *
     *         i  = io [ exp( (Betaf  Nu  nSt  F)/ RT) - exp(- ((1 - Betaf) Nu nSt F)/ RT) ]
     *
     *              where Nu = E - OCV is the overpotential
     *
     *              and E = Phi(Metal) - Phi(Soln)
     *
     *   where the current density (amps m-2) is defined as
     *
     *         i = - z_e- * F * d[e-]/dt
     *
     *    and
     *
     *         d[e-]/dt = (n) ROP_net
     *
     *   A note about signs
     *
     *    i is defined here as being positive if there is a net current density going from the metal into solution.
     *    This in turn means that if there is a positive generation of electrons, then the current
     *    is positive.
     *
     *    In this formulation the stoichiometric number of electrons can be positive or negative
     *    It will be negative if the electrons appears on the reactant side of the equation.
     *
     *    In this formulation, io can be positive or negative. It's sign will be determined by the
     *    sign of the nSt, the stoichiometric electrons.
     *
     *   \todo          Understand if this can be a const function
     *   \todo          Generalize this to the case where activities aren't included in io.
     *
     *   @param[in]    irxn                  Reaction id
     *   @param[out]   nStoich               Number of stoichiometric electrons transferred (can be negative)
     *   @param[out]   OCV                   Open circuit voltage (volts)
     *   @param[out]   io                    Exchange Current density value (can be negative) units coulombs / sec / m^2
     *   @param[out]   nu                    Overpotential for the reaction (can be positive or negative)
     *   @param[out]   beta                  Symmetry factor
     *   @param[out]   resist_ptr            Resistivity pointer
     *
     *   @return                             Returns the current density for the reaction (amps m-2)
     */
    double getExchangeCurrentDensityFormulation(size_t irxn, doublereal* nStoich, doublereal* OCV,
                                                doublereal* io, doublereal* nu, doublereal *beta, doublereal* resist_ptr);
#endif

#ifdef DONOTREMOVE
    //! Utility routine to calculate the current density given the parameters for
    //! an exchange current density formulation of the reaction rate
    /*!
     *  @param[in]  nu                       Overpotential for the reaction (can be positive or negative) (volts)
     *  @param[in]  nStoich                  Number of stoichiometric electrons transfered
     *  @param[in]  io                       Exchange Current density value (can be negative) units coulombs / sec / m^2
     *  @param[in]  beta                     Symmetry factor (unitless)
     *  @param[in]  temp                     temperature (Kelvin)
     *
     *  @return                              Returns the current density (amps m-2)
     */
    double calcCurrentDensity(double nu, double nStoich, double io, double beta, double temp) const;
#endif


 

    //!  Identify the metal phase and the electrons species
    /*!
     *   We fill in the internal variables, metalPhaseRS_ and kElectronRS_ here
     */
    void identifyMetalPhase();

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
    void deriveEffectiveChemPot();

    //! Calculate the effectit thermodynamic variables of all of the species due to the OCV override
    /*!
     *  We calculate the effects of the OCV override on the storred thermodynamics of the species
     *
     *   We modify the following storage variables
     *
     *     
     */
    void deriveEffectiveThermo();

    //!  Get the vector of deltaG values for all reactions defined in the kinetics object
    /*!
     *   (Virtual from Kinetics.h)
     *
     *   This routine provides an override to the normal calculation of deltaG, when the thermodynamics
     *   is modified by a specification of the open circuit potential.
     *
     *   @param[out]        deltaG      Vector of deltaG values. Must be at least of length equal to the number of reactions.
     */
    virtual void getDeltaGibbs(doublereal* deltaG);

    //!  Get the vector of deltaG values for all reactions defined in the kinetics object
    /*!
     *
     *   This routine provides the normal calculation of deltaG in all cases. It ignores the
     *   possible presence of an override.
     *
     *   @param[out]        deltaG      Vector of deltaG values. Must be at least of length equal to the number of reactions.
     */
    void getDeltaGibbs_Before(doublereal* const deltaG = 0);

    //! Return the vector of values for the electrochemical free energy change of reaction.
    /*!
     * These values depend upon the concentration of the solution and the  voltage of the phases
     * involved with each reaction that has a charged species as a participant.
     *
     *  units = J kmol-1
     *
     *  @param[out]        deltaM       Output vector of  deltaM's for all of the reactions. The length is m_ii.
     */
    virtual void getDeltaElectrochemPotentials(doublereal* deltaM);

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
    virtual void getDeltaEnthalpy(doublereal* deltaH);

    //!  Get the vector of deltaH values for all reactions defined in the kinetics object
    /*!
     *   This routine provides the normal calculation of deltaH, before any override modification.
     *
     *   @param[out]        deltaH      Vector of deltaH values. Must be at least of length equal
     *                                  to the number of reactions, m_ii
     */
    void getDeltaEnthalpy_Before(doublereal* const deltaH);

    //! This gets the deltaG for each reaction in the mechanism, but using the standard state
    //! chemical potential for the electrolyte.
    /*!
     *          @param  deltaG_special  DeltaG for each reaction using standard state chemical
     *                                  potentials for the electrolyte. 
     *                     length = nReactions(), J/kmol
     */
    void getDeltaGibbs_electrolyteSS(doublereal* deltaG_special);

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
    virtual void getDeltaEntropy(doublereal* deltaS);

    //!  Get the vector of deltaS values for all reactions defined in the kinetics object before OCV override
    /*!
     *   This routine does a normal calculation of deltaS, before any override of thermodynamics is carried out.
     *
     *   @param[out]        deltaS      Vector of deltaS values. Must be at least of length equal
     *                                  to the number of reactions. Units are J kmol-1 K-1.
     */
    void getDeltaEntropy_Before(doublereal* const deltaS);

    //!  Return the vector of values for the reaction standard state gibbs free energy change.  These values don't depend upon
    //!  the concentration of the solution.
    /*!   
     *  (virtual from Kinetics.h)
     *  units = J kmol-1
     *
     * @param[out]      deltaG            Output vector of ss deltaG's for reactions Length: m_ii.
     */
    virtual void getDeltaSSGibbs(doublereal* deltaG);

    //!  Return the vector of values for the change in the standard
    //! state enthalpies of reaction.  These values don't depend
    //! upon the concentration of the solution.
    /*!
     *  (virtual from Kinetics.h)
     *  units = J kmol-1
     *
     * @param[out]     deltaH              Output vector of ss deltaH's for reactions Length: m_ii.
     */
    virtual void getDeltaSSEnthalpy(doublereal* deltaH);

    
    //!  Return the vector of values for the change in the standard
    //!  state entropies for each reaction.  These values don't
    //! depend upon the concentration of the solution.
    /*! 
     *  (virtual from Kinetics.h)
     *  units = J kmol-1 Kelvin-1
     *
     *     @param[out]    deltaS           Output vector of ss deltaS's for reactions Length: m_ii.
     */
    virtual void getDeltaSSEntropy(doublereal* deltaS);

    //!   Sets the temperature and pressure for all phases that are part of the reacting surface
    /*!
     *     This calls the underlying %ThermoPhase routines for all phases that are part of the surface
     *
     *      @param[in]     temp            Temperature (Kelvin)
     *      @param[in]     pres            Pressure    (Pascal)
     */
    void setState_TP(double temp, double pres);

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

    //! Number of phases within the class
    size_t numPhases_;

    //!  Vector of pointers to xml trees
    std::vector<XML_Node*> xmlList;

    //! Mapping between the phase order in the interface object
    //! and the overall order in the PhaseList object
    /*!
     *  Note in the phase list object, surface phases are listed last.
     *  Length is the number of phases in the interface kinetics object
     *  value is the id of the phase in the PhaseList object.
     *
     *   kinOrder[kph] = iph;
     *        kph = phase index in the interface kinetics object
     *        iph = phase index in the PhaseList object
     */
    std::vector<size_t> kinOrder;

    //! Vector of the indexes of each phase in the ReactionSurfaceDomain object
    //! given the index within the PhaseList object
    /*!
     *       jph = PLtoKinPhaseIndex_[iph];
     *
     *          iph refers to the index of the phase in the Electrode_Model object 
     *          jph refers to the index of the phase in the heterogeneous kinetics object
     *
     *  Length = number of phases in the PhaseList
     *
     *  A value of -1 or npos in this slot means that the phase doesn't participate in the
     *  current ReactingSurDomain object
     */
    //std::vector<int> PLtoKinPhaseIndex_;
    std::vector<size_t> PLtoKinPhaseIndex_;

    //! Vector of the indexes of each species in the ReactionSurfaceDomain object
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
     *  current ReactingSurDomain object
     */
    std::vector<size_t> PLtoKinSpeciesIndex_;

    //! Index mapping kinetics species index to the PhaseList species index.
    /*!
     *   Length is the number of species in the kinetics species list
     */
    std::vector<size_t> KintoPLSpeciesIndex_;

    //! Global phase Index of the phase in the PhaseList object that has the kinetics
    //! object for this reacting surface.
    size_t iphaseKin_;

    //! List of names that constitute the ThermoPhases needed for the kinetics object
    /*!
     *    Currently this is needed to fix up the shallow pointer copy operation. 
     *    The list is needed to carry out shallow pointer assignments when an Electrode object is copied.
     */
    std::vector<std::string> tpList_IDs_;

    //! Temp vector that may be eliminated in the future.
    std::vector<int> tplRead;

    //! If there is a surface kinetics mechanism associated with this object, this is true. 
    /*!
     *   The index of the surface phase is kept in the variable, iphaseKin_ 
     */
    bool m_DoSurfKinetics;

    //! Vector that will expose the species production rates for this kinetics object
    /*!
     *  Length is the number of species in the kinetics vector. The units are kmol m-2 s-1.
     */
    std::vector<double> speciesProductionRates_;

    //! Vector of Limited Rates of Progress of the reactions
    /*!
     *  Length is the number of reactions, n_ii. The units are kmol m-2 s-1.
     */
    std::vector<double> limitedROP_;

    //! Vector that will expose the species creation rates for this kinetics object
    /*!
     *  Length is the number of species in the kinetics vector
     */
    std::vector<double> speciesCreationRates_;

    //! Vector that will expose the species destruction rates for this kinetics object
    /*!
     *  Length is the number of species in the kinetics vector
     */
    std::vector<double> speciesDestructionRates_;

    //! Internal Vector of deltaG of reaction for all reactions
    /*!
     *  This is used to store the DeltaG of reaction before modification due to OCV override.
     *  The length is equal to nreactions(), m_ii. The units are Joules kmol-1.
     */
    std::vector<double> deltaGRxn_Before_;

    //! Internal Vector of deltaG of reaction for all reactions (without the electrolyte mixing entropy term)
    /*!
     *  This is used to store the DeltaG of reaction before modification due to OCV override.
     *  The length is equal to nreactions(), m_ii. The units are Joules kmol-1.
     */
    std::vector<double> deltaGRxnOCV_Before_;

    //! Internal Vector of deltaH of reaction for all reactions
    /*!
     *  This is used to store the DeltaH of reaction before modification due to OCV override.
     *  The length is equal to nreactions(), m_ii. The units are Joules kmol-1.
     */
    std::vector<double> deltaHRxn_Before_;

    //! Internal Vector of deltaS of reaction for all reactions
    /*!
     *  This is used to store the DeltaS of reaction before modification due to OCV override.
     *  The length is equal to nreactions(), m_ii. The units are Joules kmol-1 K-1.
     */
    std::vector<double> deltaSRxn_Before_;

    //! Pointer to the PhaseList object that contains the ThermoPhase objects.
    /*!
     *  This object doesn't own this. However, it uses this heavily. It is a shallow pointer.
     */
    Cantera::PhaseList* m_pl;

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

    //!  Vector of the enthalpies for all species in all phases that participate in the reaction mechanism, modified for OCVoverride
    /*!
     *   We keep a vector of enthalpies here over all reaction species. This vector gets modified from its
     *   strictly thermodynamic origin when an OCV override is done.
     */
    std::vector<double> m_Enthalpies_rspec;

    //!  Vector of the enthalpies for all species in all phases that participate in the reaction mechanism
    /*!
     *   We keep a vector of enthalpies here over all reaction species. This vector does not get modified from its
     *   strictly thermodynamic origin when an OCV override is done.
     */
    std::vector<double> m_Enthalpies_Before_rspec;

    //!  Vector of the entropies for all species in all phases that participate in the reaction mechanism, modified for OCVoverride
    /*!
     *   We keep a vector of enthalpies here over all reaction species. This vector gets modified from its
     *   strictly thermodynamic origin when an OCV override is done.
     */
    std::vector<double> m_Entropies_rspec;

    //!  Vector of the entropies for all species in all phases that participate in the reaction mechanism
    /*!
     *   We keep a vector of enthalpies here over all reaction species. This vector does not get modified from its
     *   strictly thermodynamic origin when an OCV override is done.
     */
    std::vector<double> m_Entropies_Before_rspec;

    //!  Vector of the chemical potential for all species in all phases that participate in the reaction mechanism, modified for OCVoverride
    /*!
     *   We keep a vector of gibbs energy here over all reaction species. This vector gets modified from its
     *   strictly thermodynamic origin when an OCV override is done.
     */
    std::vector<double> m_GibbsOCV_rspec;

    //!  Vector of the chemical potential for all species in all phases that participate in the reaction mechanism
    /*!
     *   We keep a vector of gibbs energy here over all reaction species. This vector does not get modified from its
     *   strictly thermodynamic origin when an OCV override is done.
     */
    std::vector<double> m_Gibbs_Before_rspec;

    //!  Value of deltaG, i.e., the change in the chemical potential, for the particular species that gets modified  
    //!  due to the OCV override 
    double deltaG_species_;

    //!  Value of deltaS, i.e., the change in the entropy, for the particular species that gets modified  
    //!  due to the OCV override 
    double deltaS_species_;

    //!  Value of deltaH, i.e., the change in the enthalpy, for the particular species that gets modified  
    //!  due to the OCV override 
    double deltaH_species_;


protected:

    friend class RSD_OCVmodel;

};

}
#endif
