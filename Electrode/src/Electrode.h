/*
 * $Id: Electrode.h 604 2013-05-24 16:27:35Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_H
#define _ELECTRODE_H


#include "cantera/equilibrium.h"
#include "cantera/numerics/ResidEval.h"
#include "cantera/numerics/ResidJacEval.h"

#include "tok_input_util.h"

#include "PhaseList.h"
#include "ReactingSurDomain.h"

#include "Electrode_defs.h"
#include "ExtraGlobalRxn.h"
#include "EState.h"

#include "BlockEntry.h"
#include "Electrode_Exception.h"
#include "Electrode_input.h"

#include <string>
#include <vector>

class ELECTRODE_KEY_INPUT;
class EGRInput;
namespace Cantera
{

#ifndef SAFE_DELETE
#define SAFE_DELETE(a) if (a) { delete (a); a = 0; }
#endif


//!  Enum Type identifying the models for the Electrodes
/*!
 *  Note that this Enum class may be extended in other contexts.
 *  Therefore, loops over this enum should use a case statement with a default that either throws
 *  a catchable signal or returns a similar type of error message.
 */
enum Electrode_Types_Enum {
    UNKNOWN_ET = -1,
    BASE_TYPE_ET = 0,
    INF_CAPACITY_ET,
    MP_RXNEXTENT_ET,
    MULTIPLATEAU_NODIFF_ET,
    //! Diffusion
    SIMPLE_DIFF_ET,
    SIMPLE_PHASE_CHANGE_DIFFUSION_ET,
    CSTR_ET,
    CSTR_ZN_ANODE_ET,
    CSTR_MCMB_ANODE_ET,
    CSTR_LICO2_CATHODE_ET,
    SUCCESSIVE_SUBSTITUTION_ET,
    MP_RXNEXTENT_FES2_ET,
    MP_RXNEXTENT_LISI_ET,


    //! network of one or more radial diffusion regions with surrounding interfacial kinetics regions
    RADIAL_DIFF_REGIONS_ET,

    DIFF_TALE_ET
        
};

//!  The SOURCES enum lists the source terms that are supplied by the electrode object to the calling
//!  routine.  This is currently only used in the Jacobian wrapper. However, it may be used
//!  more widely in the future later.
enum SOURCES
{
    CURRENT_SOURCE,
    ELECTROLYTE_PHASE_SOURCE,
    ENTHALPY_SOURCE,
    SPECIES_SOURCE,
    MAX_SOURCE
};

//!  The DOFS enum lists the independent degrees of freedom that the electrode object can
//!  handle. This is currently only used in the Jacobian wrapper. However, it may be used
//!  more widely in the future later.
enum DOFS
{
    SOLID_VOLTAGE,
    LIQUID_VOLTAGE,
    TEMPERATURE,
    PRESSURE,
    SPECIES,
    MAX_DOF
};

//! Interpolation of exterior fields
enum Electrode_Exterior_Field_Interpolation_Scheme_Enum {

    //! Field interpolation scheme based on assuming final_final values persist throughout the time step
    T_FINAL_CONST_FIS = 0,

    //! Field interpolation scheme based on assuming that external fields linearly vary
    //! from the init_init_ values to the final_final values.
    LINEAR_INTERP_FIS
};


//! Base structure for storing the external state of the electrode at a particular time
/*!
 *  There will be external states for the t_init_init and t_final_final times.
 *  there will be an interpolation strategy between the two.
 */
struct Electrode_ExternalField_State {

    double phiMetal_;

    double phiElectrode_;

    std::vector<double> electrolyteMoleFractions_;

    double Temperature_;

    double Pressure_;
};

class Electrode_Equilibrium;

//! NEW CONCEPT
//!  Integration routine needs to know the context in which it is being called in order
//!  to make decisions about how it will be treating the time stepping.
/*!
 *  This enum provides the variations of the context with which it is being called.
 */
enum Subgrid_Integration_RunType_Enum {

    //! Time integration stategy is developed for the first time along the interval and storred
    //! for later use.
    /*!
     *  This is used in the calculation of numerical deltas of the Electrode objects, where
     *  we want to minimize the noise in the calculation due to time stepping.
     *   -> this has been the default
     */
    BASE_TIMEINTEGRATION_SIR= 0,

    //! Time integration stategy is reevalulated after the first time along the interval and storred
    //! for later use.
    /*!
     *  This is used in the calculation of numerical deltas of the Electrode objects, where
     *  we want to minimize the noise in the calculation due to time stepping. Here,
     *  we use the previous history of the time step over the interval to smooth out the calculation
     *  so as to equalize the errors incurred between steps.
     *
     *   (HKM -> delay implementation or delete)
     */
    BASE_REEVALUATION_TIMEINTEGRATION_SIR,

    //! Time integration strategy is fixed at a set number of subcycles from the input file
    /*!
     *   Special cases can change the actual number of subcylces used to be greater than
     *   the requested number
     */
    FIXEDSUBCYCLENUMBER_TIMEINTEGRATION_SIR,

    //! Reuse the time integration stategy using a delta of the value of the field variables
    /*!
     *  This is used in the  calculation of numerical deltas of the Electrode objects, where
     *  we want to minimize the noise in the calculation due to time stepping.
     */
    FVDELTA_TIMEINTEGRATION_SIR
};


//! Electrode class is the base class used to model porous electrodes
/*!
 * Complete problem statement
 *
 *  This class also serves as the placeholder for the state of the electrode
 *
 *  The basic structure is the PhaseList. The volume phases are listed first.
 *  then the surface phases are listed. All vectors over phases in this
 *  object follow that pattern. Species list in this object also follow
 *  that pattern, with all species grouped by phase first.
 *
 *  It may be noted that this class uses extensive variables to describe the
 *  volume of material, the number of moles of solid phase and electrolyte
 *  species, and the surface area of each reacting surface on the particle.
 *
 *  In this base class, there is only one reacting surface. The surface area is constant.
 *  The activity of the electrode may change due
 *
 *  The bulk phases are considered to be uniform. They may change in mole numbers
 *
 *  The Morphology consists of a bunch of sphere which touch each other in close
 *  packed form.
 *
 *  The surface area may be calculated from the above morphology.
 *
 *  Surface areas are separated into two categories: external and internal. External
 *  surface areas are defined to be surfaces involving the surrounding electrolyte.
 *  Internal surfaces are defined to be surface areas between internal solid regions.
 *
 *  The electrical conductivity is a linear combination of volumes of the bulk phase
 *  electrical conductivity.
 *
 *
 *  There is a convention for the sign of the current.  The current is positive for
 *  current going into the electrode  and then into the electrolyte. Thus, under
 *  normal battery operation where the anode is negative and the cathode is positive,
 *  the current is positive going into the anode and negative going into the cathode.
 *
 *
 *  More detailed models will be created in child objects.
 *
 *  Calculation of source terms within the object
 * ---------------------------------------------------------
 *  The electrode contains an idea of morphology. The morphology may change as a function
 *  of time. There may be internal variables within the electrode object that are not
 *  presented to the exterior of the object. Therefore, all source term calculations must
 *  actually be couched in terms of integrations in time. There is an initial state, which
 *  is kept within the electrode object. The integration occurs over the time step to the
 *  final state. The final state is kept within the object.
 *
 *  The integration depends on the variables in the final state which are shared within the
 *  code surrounding the electrode object. For example, we will take a backwards euler approximatino
 *  The external variables of the electrode are shared with the surrounding code. The backwards
 *  euler approximation assumes that these external variables are constant during the
 *  time integration and equal to the values at the end of the integration. For stiff
 *  problems, a jacobian is needed to be calculated by the electrode object consisting of the
 *  the dependence of the source term with respect to these external variables.
 *  Therefore, the integration has to be carried out for a matrix of conditions, consisting
 *  of variations of the solution with respect to these external variables.
 *  Note that this jacobian may only be approximate as we are only really concerned with relaxing
 *  the system with the jacobian.
 *
 *  The time step control which controls the length of the integration depends upon an
 *  estimate of the time step error. This is obtained by taking the difference between
 *  the predicted and the final solution at the current time step. The predicted solution
 *  is obtained by integrating the Electrode equations using the external variables obtained
 *  at the previous time step. for stiff systems this may lead to instabilities in the time
 *  stepping without the implicit be step. however, one can get the time step error this way.
 *
 *
 *
 *    this will include a distributed solid domain.
 *    This will include multiple materials that create multiple plateaus.
 *
 *      Classifications of member functions
 * ----------------------------------------------------
 *
 *    Each function should be categorized by how they treat virtual functions and whether they can be called during the
 *    creation phase
 *
 *       serial virtual - post creation       Virtual functions that call child member functions one or more times.
 *                                            These are restricted to post-creation phase.
 *       virtual onion Out                    Virtual or non-virtual functions that don't call child members.
 *                                            carried out from parent to child
 *       virtual onion In                     Virtual or non-virtual functions that don't call child members.
 *                                            carried out from child to parent
 *       virtual overloaded                   Virtual functions that complete eliminate the functionality
 *                                            of parent member functions.
 *
 *    Each function should be categorized by which state they interact with. Most functions will change
 *    the final state. Some functions can only be called in between integration steps because they assume
 *    that all states are the same.
 *            -> during integration, affects final state or final_final state (post creation)
 *            -> after integration, affects final state  (post creation)
 *            -> between integration steps (after call to reset_condition or before any integration steps)
 *                     affects all states equally.
 *
 *    Origin of virtual
 *            -> State the base class that first assigns the member function os virtual
 *
 */
class Electrode : public Cantera::PhaseList
{
public:

    //! Static function to initialize static variables that are read from the environment.
    /*!
     *     Environmental Variable                  Static Variable
     *     ELECTRODE_TURN_OFF_PC_PRINTING                     s_printLvl_PREDICTOR_CORRECTOR
     */
    static void readEnvironmentalVariables();

    // -----------------------------------------------------------------------------------------------------------------
    // ------------------------------------- BASIC SETUP ROUTINES  -----------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------

    //! Constructor
    Electrode();

    //! Destructor
    virtual ~Electrode();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    Electrode(const Electrode& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    Electrode& operator=(const Electrode& right);

    //! Duplicator function
    /*!
     *  Duplicate the current electrode object, returning a base electrode pointer
     */
    virtual Electrode* duplMyselfAsElectrode() const;

    //! Set the electrode ID information
    void setID(int domainNum, int cellNum);

    //! Return the type of electrode
    /*!
     *  Returns the enum type of the electrode. This is used in the factory routine.
     *  Note, this does not indicate whether the electrode is being used as
     *  an anode or a cathode.
     *
     *  @return Returns an enum type, called   Electrode_Types_Enum
     */
    virtual Electrode_Types_Enum electrodeType() const;


    //! Add additional Keylines for child electrode objects, and then read them in
    /*!
     *   (virtual function from Electrode)
     *   (overload virtual function - replaces parent members)
     *
     *   This function will replace the ELECTRODE_KEY_INPUT structure with an expanded
     *   child member structure containing the extra information.
     *
     *   If the command file has been read before, it will then reparse the command file
     *   storring the new information in the ELECTRODE_KEY_INPUT structure.
     *
     *    @param ei  Handle to the ELECTRODE_KEY_INPUT base pointer. This handle may change
     *               as the child class of  ELECTRODE_KEY_INPUT gets malloced.
     *
     *    @return  0 successful but no change in ei
     *             1 Successful and ei has changed
     *            -1 unsuccessful fatal error of some kind.
     */
    virtual int electrode_input_child(ELECTRODE_KEY_INPUT** ei);

    //!  Setup the electrode for first time use
    /*!
     * (virtual from Electrode  - onion Out)
     *
     *    This is one of the most important routines. It sets up the electrode's internal structures
     *    After the call to this routine, the electrode should be internally ready to be integrated
     *    and reacted.
     *    It takes its input from an ELECTRODE_KEY_INPUT object which specifies the setup of the electrode
     *    object and the initial state of that object.
     *    The routine works like an onion Out initialization. The parent object is initialized before the
     *    child. This means the child object first calls the parent, before it does its own initializations.
     *
     *    There are some virtual member functions that won't work until this routine is called.
     *    That's because the data structures won't be set up for base and child Electrode objects until 
     *    this is called.
     *
     *  @param ei   BASE ELECTRODE_KEY_INPUT pointer object. Note, it must have the correct child class
     *              for the child electrode object.
     *
     *  @return  Returns zero if successful, and -1 if not successful.
     */
    virtual int electrode_model_create(ELECTRODE_KEY_INPUT* ei);

    //!  Set the electrode initial conditions from the input file.
    /*!
     *   (virtual from Electrode)
     *   (This is a serial virtual function or an overload function)
     *
     *    This is one of the most important routines. It sets up the initial conditions of the electrode
     *    from the input file. The electrode itself has been set up from a call to electrode_model_create().
     *    After the call to this routine, the electrode should be internally ready to be integrated and reacted.
     *    It takes its input from an ELECTRODE_KEY_INPUT object which specifies the setup of the electrode
     *    object and the initial state of that object.
     *
     *    The routine works like an onion initialization. The parent object is initialized before the
     *    child. This means the child object first calls the parent, before it does its own initializations.
     *
     * @param ei    ELECTRODE_KEY_INPUT pointer object
     *
     *  @return  Returns zero if successful, and -1 if not successful.
     */
    virtual int setInitialConditions(ELECTRODE_KEY_INPUT* ei);

    //! Create an object that saves the electrode state and can print out an XML solution to file
    /*!
     *  The pointer to the malloced object is saved in the internal variable eState_final_ .
     *  Because there is an object, the state of the electrode will be saved at each step.
     *
     *  @param ei Electrode Key Input
     *
     *  @return  Returns zero if successful, and -1 if not successful.
     */
    virtual int electrode_stateSave_create();

    //! Set the sizes of the electrode from the input parameters
    /*!
     *  We resize all of the information within the electrode from the input parameters
     *
     * @param electrodeArea       Area of the electrode
     * @param electrodeThickness  Width of the electrode
     * @param porosity            Volume of the electrolyte phase
     */
    virtual void setElectrodeSizeParams(doublereal electrodeArea, doublereal electrodeThickness,
                                        doublereal porosity);

protected:
    //! Resize the solid phase and electrolyte mole numbers within the object
    /*!
     *  (virtual from Electrode)
     *  This routine uses particleDiameter_ , particleNumberToFollow_, and porosity_ to recalculate
     *  all the mole numbers in the electrode. This is done by rescaling all of the numbers.
     *  At the end of the process, the total volume of the electrode object is
     *
     *    grossVol = SolidVol() / ( 1.0 - porosity_)
     *
     *  where the SolidVol() is equal to
     *
     *   SolidVol() =  particleNumberToFollow_  Pi *  particleDiameter_**3 / 6.0;
     *
     */
    virtual void resizeMoleNumbersToGeometry();

public:
    //! Resize the electrolyte mole numbers so that it consistent with the porosity
    /*!
     *  This routine resizes the electrolyte mole numbers within the object to be consistent
     *  with the porosity. This has the effect of changing the total volume of the
     *  object.
     *
     *  @param porosityReset   If positive, this resets the porosity within the electrode
     *                         to the given value.
     */
// Uncalled
    void resizeSolutionNumbersToPorosity(doublereal porosityReset = -1.0);

    //! resize the surface areas according to the input geometry.
    /*!
     *  We resize the surface areas of the Reacting Surfaces to a value which is
     *  equal to the particule exterior surface area multiplied by the number of particles.
     *
     *  Note we can't do anything better in this routine because we do not know the specifics
     *  of the particles.
     */
//Can protect
    void resizeSurfaceAreasToGeometry();

    // We change the particleNumbertoFollow_ field to comply with the number of moles and the particle diameter
//Can protect
    void resizeParticleNumbersToMoleNumbers();


    //! Returns a pointer to the current outer reacting surface object
    /*!
     *  We currently assume that the first Active reacting surface is the outer surface
     *
     *  @return Returns a pointer to first active reacting surface
     */
    Cantera::ReactingSurDomain* currOuterReactingSurface();

    //! Returns a reacting surface object
    /*!
     *  @param iSurf  returns a pointer to that reacting surface
     */
    Cantera::ReactingSurDomain* reactingSurface(int iSurf);

    //! Set the temperature and pressure of the electrode
    /*!
     *   Depends on Serial Virtual Functions that must be used after electrode_model_create()
     *
     * @param temperature    Temperature (Kelvin)
     * @param pressure       Pressure (Pa)
     */
    void setState_TP(double temperature, double pressure);

    //! return the current temperature
    /*!
     *   Kelvin
     */
    double temperature() const;

    //! Return the current pressure
    /*!
     *  units = Pa
     */
    double pressure() const;

    //! Return the current porosity
    double porosity() const;

    //
    // --------------------------------------  QUERY VOLUMES -----------------------------------------------------------
    //
    //!    Return the total volume of solid material
    /*!
     *  (virtual from Electrode.h)
     *
     *       This is the main routine for calculating the
     *       We leave out the solnPhase_ volume from the calculation
     *       units = m**3
     */
    virtual double SolidVol() const;

    //! Return the total volume in the electrode
    /*!
     *  (virtual from Electrode.h)
     *
     *   This returns the value in the _final_ state
     *
     *   @return total volume (m**3)
     */
//Can protect?
    virtual double TotalVol(bool ignoreErrors = false) const;

    //!  Returns the total moles in the electrode phases of the electrode
    /*!
     *  @return total moles (kmol)
     */
//Can protect
    virtual  double SolidTotalMoles() const;

    //! Return a vector of the phase volumes for all phases in the electrode
    /*!
     *  Note the vector is over surface phases as well. Currently, all surface phases have zero
     *  volume.
     *
     *  length = m_NumTotPhases
     *  units = m**3
     */
// Deprecate
    virtual void getPhaseVol(double* const phaseVols) const;
    //
    // --------------------------------------  QUERY HEAT CAPACITY  -----------------------------------------------------
    //

    //!  Returns the total Heat Capacity of the Material in the Solid Electrode at constant volume
    /*!
     *  This is an extensive quantity.
     *
     *  @return Joule K-1
     */
    virtual double SolidHeatCapacityCV() const;

    //-------------------------------------------------------------------------------------------------------------------
    // --------------------------------------------- SURFACE AREAS ------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------
    //! Return the number of surfaces in the Electrode object
    /*!
     *  The number of surfaces may differ from the number of surface kinetics objects now.
     *  All references to surfaces within the Electrode loop over the number of surfaces,
     *  with an array determining there is an Electrode kinetics object located at that surface.
     */
    int nSurfaces() const;

    //! Return the number of reactions in the a reacting surface domain
    int nReactions(int isk) const ;

    //! Return a vector of surface phase areas.
    /*!
     *  Returns a vector of surface areas for surfaces defined in the electrode object.
     *
     * @param surfArea  vector of surface areas.
     *                  units = m2
     *                  length = number of surfaces defined in the Electrode
     */
    void getSurfaceAreas(double* const surfArea) const;

    //-------------------------------------------------------------------------------------------------------------------
    // ----------------------------------- QUERY AND SET VOLTAGES -------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------
    //! This sets the metal and solution voltage
    /*!
     *  This sets the metal and the electrolyte voltages
     *  This assumes that there are only two voltages in the system.
     *   The voltage of the interface is defined as VOLTS = phiMetal - phiElectrolyte
     *
     *  @param[in]  phiMetal     Potential of the metal
     *  @param[in]  phElectrolyte Potential of the electrolyte
     */
    void setVoltages(const double phiMetal, const double phiElectrolyte);

    //! This returns the voltage of the electrode at the final state conditions
    /*!
     *  The voltage of the electrode is defined as the potential of the metal minus
     *  the potential of the electrolyte.
     *
     *  @return  returns voltage in volts.
     */
    double voltage() const;

    //! Return the voltage of a phase
    /*!
     * @param iph  Phase id
     *
     * @return Returns the voltage in volts
     */
    double phaseVoltage(int iph) const;

    //! Set the voltage of a phase
    /*!
     * @param iph  Phase id
     * @param volts volts
     */
    void setPhaseVoltage(int iph, double volts);

    //-------------------------------------------------------------------------------------------------------------------
    // -------------------------------- QUERY AND SET MOLE NUMBERS-------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------

    //! Set the mole numbers in the electrolyte phase
    /*!
     *  We set the mole numbers of the electrolyte phase separately from
     *  the rest of the phases.
     *
     *  We always make sure that mole numbers are positive by clipping. We
     *  always make sure that mole fractions sum to one.
     *
     *  If we are not following mole numbers in the electrolyte, we set the
     *  total moles to the internal constant, electrolytePseudoMoles_, while
     *  using this vector to set the mole fractions.
     *
     * @param electrolyteMoleNum vector of mole numbers of the species in the
     *                     electrolyte phase.
     *                     units = kmol
     *                     size = number of species in the electrolyte phase
     * @param setInitial   Boolean indicating that we should set the initial values of the
     *                     electrolyte mole numbers as well. We also set init_init and 
     *                     final_final values too. 
     */
    virtual void setElectrolyteMoleNumbers(const double* const electrolyteMoleNum, bool setInitial);

    //! Set the mole numbers in a single phase
    /*!
     *  We set the mole numbers of a single phase separately from the rest of the phases.
     *
     *  We always make sure that mole numbers are positive by clipping. We
     *  always make sure that mole fractions sum to one.
     *
     *  If we are not following mole numbers in the electrode, we set the
     *  total moles to the internal constant, electrolytePseudoMoles_, while
     *  using this vector to set the mole fractions.
     *
     * @param iph     Phase id.
     * @param moleNum vector of mole numbers of the species in the
     *                     electrolyte phase.
     *                     units = kmol
     *                     size = number of species in the electrolyte phase
     */
    virtual void setPhaseMoleNumbers(int iph, const double* const moleNum);



    //! Update all mole numbers in the object from the mole numbers in the spMoles_final_[] vector
    /*!
     *  We use the field spMoles_final_[] to set the field phaseMoles_final_[].
     *
     *  We set the mole numbers of a single phase separately from the rest of the phases.
     *
     *  We do not clip the mole numbers to be positive. We allow negative mole numbers.
     *
     *  We make sure that the mole fractions sum to one.
     *
     *   The following fields in this object are set:
     *
     *            spMf_final_[]
     *            VolPM_[]
     *            spElectroChemPot_[]
     *
     *            phaseMoles_final_[iph]
     *            phaseVoltages_[iph]
     *            phaseMolarVolumes_[iph]
     *
     *  If we are not following  the mole numbers in the electrode, we set the
     *  total moles to the internal constant, electrolytePseudoMoles_, while
     *  using this vector to set the mole fractions, using the ThermoPhase object
     *  to get the mole fractions.
     *
     * @param iph     Phase id.
     */
//Can protect
    void updatePhaseNumbers(int iph);


    // -----------------------------------------------------------------------------------------------------------------
    // ------------------------------------ CALCULATE INSTANTANEOUS RATES ----------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------
    //     -- These are now called "ProductionRate"  or SourceTerms values

    //! Calculate the instantaneous time derivative of the species vector as determined by all source terms
    /*!
     *  This is the rate of change in the moles of species defined in the electrode  at t_final.
     *  This calculation does not necessarily use an interval of time to calculate anything.
     *
     *  @param spMoleDot   The end result in terms of the rate of change in moles of species in the
     *                     electrode. (kmol s-1)
     */
    virtual void speciesProductionRates(doublereal* const spMoleDot);

    //!  Returns the current and the net production rates of all species in the electrode object
    //!  at the current conditions from one surface kinetics object
    /*!
     * @param isk   Surface index to get the net production rates from
     * @param net   Species net production rates [kmol/m^2/s]. Return the species
     *
     *   @return   Returns the current columb sec-1 m-2
     */
    doublereal getNetSurfaceProductionRatesCurrent(const int isk, doublereal* const net) const;

    //! Get the net production rates of all species in the electrode object
    //! at the current conditions from one surface kinetics object
    /*!
     * @param isk   Surface index to get the net production rates from
     * @param net   Species net production rates [kmol/m^2/s]. Return the species
     *
     */
//Can protect
    void getNetSurfaceProductionRates(const int isk, doublereal* const net) const;

    //!  Returns the current and the net production rates of the phases in kg/m2/s from a single surface
    /*!
     *  Returns the net production rates of all phases from reactions on a single surface
     *
     *  @param isk Surface ID to get the fluxes from.
     *  @param phaseMassFlux  Returns the mass fluxes of the phases
     */
// Deprecate
    void getPhaseMassFlux(const int isk, doublereal* const phaseMassFlux);


    //!  Returns the current and the net production rates of the phases in kmol/m2/s from a single surface
    /*!
     *  Returns the net production rates of all phases from reactions on a single surface
     *
     *  @param isk Surface ID to get the fluxes from.
     *  @param phaseMassFlux  Returns the mass fluxes of the phases
     */
// Deprecate
    void getPhaseMoleFlux(const int isk, doublereal* const phaseMoleFlux);

    void getPhaseProductionRates(const double * const speciesProductionRates, doublereal* const phaseMoleFlux) const;


    //! Overpotential term for the heat generation from a single surface
    /*!
     *   @param irxn Surface index 
     */
    virtual double thermalEnergySourceTerm_overpotential(int isk);

    //! Reversible Entropy term leading to  heat generation
    /*!
     *    @param isk   
     */
    virtual double thermalEnergySourceTerm_reversibleEntropy(size_t isk);

    //! Reversible Entropy term leading to  heat generation
    /*!
     *
     */
    virtual double thermalEnergySourceTerm_EnthalpyFormulation(size_t isk);

    // -----------------------------------------------------------------------------------------------------------------
    // ------------------------------------ CARRY OUT INTEGRATION OF EQUATIONS -----------------------------------------
    // -----------------------------------------------------------------------------------------------------------------

    //! Get the final state time
    /*!
     *  @return  Returns the time used in the _final_ state.
     */
    double getFinalTime() const;

    //! Returns the initial global time
    /*!
     *  @return Returns the initial time
     */
    double timeInitInit() const;
  
    //! Returns the final global time
    /*!
     *  @return Returns the final time
     */
    double timeFinalFinal() const;

    //! The internal state of the electrode must be kept for the initial and final times of an integration step.
    /*!
     *  This function advances the initial state to the final state that was calculated
     *  in the last integration step.
     *
     * @param Tinitial   This is the new initial time. This time is compared against the "old"
     *                   final time, to see if there is any problem.
     */
    virtual void resetStartingCondition(double Tinitial, bool doTestsAlways = false);

    //! Revert the object's conditions to the initial conditions
    /*!
     *  If there was a solution in t_final, this is wiped out and replaced with the solution at t_init.
     *  We get rid of the pending flag here as well.
     *
     * @param revertToInitInit revert to the t_init_init solution
     *                         Defaults to true
     */
    virtual void revertToInitialTime(bool revertToInitInit = true);

    //!  Calculate the change in the state of the system when integrating from T_initial_initial_
    //!  to t_final_final_
    /*!
     *  All information is kept internal within this routine. This may be done continuously
     *  and the solution is not updated.
     *
     *  Note the tolerance parameters refere to the nonlinear solves within the calculation
     *  They do not refer to time step parameters.
     *
     *  @param deltaT        DeltaT for the integration step.
     *  @param GlobalRtolSrcTerm Relative tolerance for the source term vector calcualted from
     *                       this routine.
     *                       Defaults to 1.0E-3
     *  @param fieldInterpolationType Type of interpolation of field variables defaults to T_FINAL_CONST_FIS,
     *  @param subIntegrationType     Type of subintegration. Defaults to BASE_TIMEINTEGRATION_SIR.
     *                                In this integration, the program determines its own strategy
     *                                for the time step.
     *
     *  @return Returns the number of subcycle steps it took to complete the full step.
     *          Failures to complete the integration due to time truncation error issues return a -1.
     *          Failures due to invalid function calculation attempts return a -2.
     *          Failures due to invalid arguments return a -3.
     */
    virtual int integrate(double deltaT, double  GlobalRtolSrcTerm = 1.0E-3,
                          Electrode_Exterior_Field_Interpolation_Scheme_Enum fieldInterpolationType = T_FINAL_CONST_FIS,
                          Subgrid_Integration_RunType_Enum subIntegrationType = BASE_TIMEINTEGRATION_SIR) = 0;


    //!  Residual calculation for the solution of the Nonlinear integration problem
    /*!
     *    Given tfinal and delta_t, and given y and ydot which are estimates of the solution
     *    and solution derivative at tfinal, this function calculates the residual equations.
     *    It is the residual function used in the nonlinear solver that relaxes the equations
     *    at each time step.
     *
     *    This is typically called from evalResidNJ(), which is called directly from the
     *    nonlinear solver. However, we expose this routine so that the residual can be queried
     *    given all of the inputs.
     *
     * @param[in]  tfinal      Time                    (input)
     * @param[in]  delta_t     The current value of the time step (input)
     * @param[in]  y           Solution vector (input, do not modify)
     * @param[in]  ydot        Rate of change of solution vector. (input, do not modify)
     * @param[out] resid       Value of the residual that is computed (output)
     * @param[in]  evalType    Type of the residual being computed (defaults to Base_ResidEval)
     * @param[in]  id_x        Index of the variable that is being numerically differenced to find
     *                         the jacobian (defaults to -1, which indicates that no variable is being
     *                         differenced or that the residual doesn't take this issue into account)
     * @param[in]  delta_x     Value of the delta used in the numerical differencing
     *
     * @return                Returns an integer that gets fed back through evalResidNJ() to the
     *                        nonlinear solver. Anything other than a 1 causes an immediate failure
     *                        of the nonlinear solver to occur.
     */
    virtual int integrateResid(const doublereal tfinal, const doublereal delta_t,
                               const doublereal* const y, const doublereal* const ydot,
                               doublereal* const resid,
                               const ResidEval_Type_Enum evalType, const int id_x, const doublereal delta_x);

    //!  Calculate the change in the state of the system when integrating from Tinitial to Tfinal
    //!  at constant current, finding the required voltage to produce that current in an
    //!  average sense.
    /*!
     *  This function finds the voltage that is necessary to produce an average current over
     *  an interval of time, deltaT.
     *
     *  There are instances where the returned deltaT may be different than the final deltaT.
     *  The algorithm reduces deltaT if there is a problem with the original deltaT,
     *  in order to find the requested current at a shorter time step.
     *
     *  This routine is designed to return an integration result in all cases. That's right.
     *  However, the returned current and the returned deltaT must
     *  all be checked to see if they are the requested values. Frequently the returned voltage
     *  will return at the requested max and min values, indicating that we are at a boundary
     *  or that we have run out of electrons in the electrode.
     *
     *  @param current    Current in amps. On return it returns the current actually obtained.
     *
     *  @param deltaT     DeltaT for the integration step. On output it contains the actual
     *                    deltaT taken. This value may be different than the input deltaT.
     *                    The algorithm reduces deltaT if there is a problem with the original
     *                    deltaT, in order to find the requested current at a shorter time
     *                    step.
     *
     *  @param phiMax     Maximum value of the voltage. Defaults to  100.
     *                    This is an input. The voltage must be bounded for the rootfinder to
     *                    work efficiently and to avoid overflows, roundoff errors due to
     *                    large reaction rate constants.
     *
     *  @param phiMin     Minimum value of the voltage. Defaults to -100.
     *                    This is an input.
     *
     *  @param maxIntegrationSteps max number of integration steps. defaults to 5000
     *
     *  @return Return the voltage used to obtain the current.
     */
    virtual double integrateConstantCurrent(doublereal& current, doublereal& deltaT,
                                            double phiMax = 100., double phiMin = -100.,
                                            int maxIntegrationSteps = 5000);

    //! Set the deltaT used for the subcycle step
    /*!
     *  The default setting for this is infinite. If it's infinite then deltaTSubcycle is set to the
     *  deltaT of the integration step. However, there are many times where a smaller natural delta T is
     *  required to solve the equation system
     */
// Deprecate? Is there any reason not to have child objects handle this as part of their adaptive time stepping?
    virtual void setDeltaTSubcycle(doublereal deltaTSubcycle);

    //! Set the maximum value of the subcycle delta T
    /*!
     *  @param deltaTSubcycle  Maximum value of the subcycle time step
     */
// Deprecate
    void setDeltaTSubcycleMax(doublereal deltaTSubcycle);

    //!  Calculate the time derivatives of the mole numbers at the current subcycle step
    /*!
     *   This may be used as a virtual function
     */
// Deprecate
    virtual void calculateTimeDerivatives(doublereal deltaTsubcycle);

    //------------------------------------------------------------------------------------------------------------------
    // ----------------------------------------  ERROR ANALYSIS OF THE INTEGRATION -------------------------------------
    //------------------------------------------------------------------------------------------------------------------

    //! Report the number of state variables and their relative integration errors during the
    //! current global integration step
    /*!
     *   Note rtol doesn't factor into this immediately. Therefore, a value or 1E-3
     *                                  would mean the error in the value is 1 part in 1000.
     *
     *  @param[out] numSV               Returns the number of state variables
     *  @param[out] errorVector         Returns a vector of errors in the state variables for the global step
     *                                  Note rtol doesn't factor into this immediately. Therefore, a value or 1E-3
     *                                  would mean the error in the value is 1 part in 1000.
     *  @return     Returns the largest value of the errors in the errorVector.
     */
    double reportStateVariableIntegrationError(int& numSV, double* const errorVector) const;

    //------------------------------------------------------------------------------------------------------------------
    // ---------------------------- INTEGRATED SOURCE TERM QUERIES -----------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------

    //! Report the integrated source term for the electrode over an interval in time.
    /*!
     *  This is the net change in the moles of species defined in the electrode over that
     *  interval of time. The conditions at the end of the interval are used to carry
     *  out the integrations.
     *
     *  @param spMoleDelta The end result in terms of the change in moles of species in the
     *                     electrode.
     *
     *  @return Tfinal    Final time to integrate to.
     *
     */
    virtual double integratedSpeciesSourceTerm(doublereal* const spMoleDelta);

    //! Report the enthalpy source term for the electrode over an interval in time
    /*!
     *   FIX!!!
     *  HKM -> This formula is in error for standard reactions for multispecies phases,
     *         and it is wholly inadequate for electrode reactions. 
     *  Sum over phases ( enthalpy phase * (phaseMoles_final_ - phaseMoles_init_init_) )
     *  This should only be called after integrate() has finished running.
     */
    virtual double integratedEnthalpySourceTerm();

  

    //! Energy released 
    /*!
     *     Energy released within the electrode during a local time step
     *
     *   @param energy released (joules)
     */
    virtual double thermalEnergySourceTerm_EnthalpyFormulation_SingleStep();

    virtual double thermalEnergySourceTerm_ReversibleEntropy_SingleStep();
    virtual double thermalEnergySourceTerm_Overpotential_SingleStep();

    //! Get the integrated source term values for one of a set of sources
    /*!
     *     @param sourceType   The enum source term value. Species indecises are 
     *                         designated by indexing on top of the base SPECIES_SOURCE
     *
     *     @return Returns the source term
     */
    virtual double getIntegratedSourceTerm(SOURCES sourceType);

    //!  Returns the net production rates of all species in the electrode object
    //!  over the last integration step
    /*!
     *  We calculate a rate here by taking the total production amounts and then
     *  divide by the time step.
     *
     *   @param net   Species net production rates [kmol/s].
     *   @return   Returns the current, amps = columb sec-1
     */
    doublereal getIntegratedProductionRatesCurrent(doublereal* const net) const;

    //!  Returns the net current in the electrode object
    //!  at the current conditions over the current last local time step
    /*!
     *   Note we must have integrated a local time step previously.
     *
     *   @return   Returns the current columb sec-1 = amps
     */
//Can protect
    doublereal integratedLocalCurrent() const;

    //!  Returns the net production rates of all species in the electrode object
    //!  over the last integration step
    /*!
     *  We calculate a rate here by taking the total production amounds and then
     *  divide by the time step.
     *
     *   @param net   Species net production rates [kmol/s].
     */
//Can protect
    void getIntegratedProductionRates(doublereal* const net) const;

  
    //!  Returns the net current in the electrode object
    //!  at the current conditions over the current global time step
    /*!
     *   Note we must have integrated a global time step previously.
     *
     *   @return   Returns the current columb sec-1 = amps
     */
//Can protect
    doublereal integratedCurrent() const;

    //! Returns the integrated moles transfered for each phase in the electrode object
    //! over the time step
    /*!
     *    @param  phaseMolesTransfered vector of moles transfered (length = number of total
     *            phases in the electrode object)
     *            units = kmol
     */
    virtual void getIntegratedPhaseMoleTransfer(doublereal* const phaseMolesTransfered);

    //! Returns the integrated thermal energy source term (Joules)
    /*!
     *    Returns the heat release that has occurred during the global time step. 
     *
     *  @param Returns the heat release (joules)
     */
    virtual double getIntegratedThermalEnergySourceTerm();

    //! The thermal energy source term can be broken into two parts. This part is the irreversible
    //! heat generation term due to the non-zero overpotential 
    /*!
     *   @return returns an 
     */
    virtual double getIntegratedThermalEnergySourceTerm_overpotential();

    //! The thermal energy source term can be broken into two parts. This part is the reversible
    //! heat generation term due to the entropy change of reaction
    /*!
     *   @return returns a double
     */
    virtual double getIntegratedThermalEnergySourceTerm_reversibleEntropy();


    // -----------------------------------------------------------------------------------------------------------------
    // ---------------------------- SOLUTION OF NONLINEAR TIME DEPENDENT SYSTEM  ---------------------------------------
    // -----------------------------------------------------------------------------------------------------------------

    //! Class used in the nonlinear solve of the integration equations
    /*!
     *   This is a child of the ResidJacEval class
     */
// Deprecate
    class integrate_ResidJacEval : public Cantera::ResidJacEval
    {
    public:
        integrate_ResidJacEval(Electrode* ee);

        //! Evaluate the residual function
        /*!
         * @param t             Time                    (input)
         * @param delta_t       The current value of the time step (input)
         * @param y             Solution vector (input, do not modify)
         * @param ydot          Rate of change of solution vector. (input, do not modify)
         * @param resid         Value of the residual that is computed (output)
         * @param evalType      Type of the residual being computed (defaults to Base_ResidEval)
         * @param id_x          Index of the variable that is being numerically differenced to find
         *                      the jacobian (defaults to -1, which indicates that no variable is being
         *                      differenced or that the residual doesn't take this issue into account)
         * @param delta_x       Value of the delta used in the numerical differencing
         */
        int  evalResidNJ(const doublereal t, const doublereal delta_t,
                         const doublereal* const y,
                         const doublereal* const ydot,
                         doublereal* const resid,
                         const ResidEval_Type_Enum evalType = Base_ResidEval,
                         const int id_x = -1,
                         const doublereal delta_x = 0.0);

        int  getInitialConditions(const doublereal t0, doublereal* const y,
                                  doublereal* const ydot);


        //! Return the number of equations in the equation system
        int nEquations() const;

        //! Pointer that is used as a self-reference
        Electrode* ee_;
    };
    // -----------------------------------------------------------------------------------------------------------------

// Deprecate
    int phasePopResid(int iphaseTarget, const double* const Xf_phase,
                      double deltaTsubcycle, double* const resid);


// Deprecate
    int phasePop(int iphaseTarget, double* const Xmf_stable, double deltaTsubcycle);

    //------------------------------------------------------------------------------------------------------------------
    // -------------------------------  SetState Functions -------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------

    //! Set the time
    /*!
     *   Set the time
     */
    void setTime(double time);


  private:
    //! Set the internal initial intermediate and initial global state from the internal final state
    /*!
     *  (non-virtual function)  -> function should onionize in-first.
     *
     *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
     *  routine as well.
     *
     * @param setInitInit   Boolean indicating whether you should set the init_init state as 
     */
    void setInitStateFromFinal_Oin(bool setInitInit = false);
  public:
    //! Set the internal initial intermediate and initial global state from the internal final state
    /*!
     *  (virtual function)
     *
     *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
     *  routine as well.
     *
     * @param setInitInit   Boolean indicating whether you should set the init_init state as well
     */
//Can protect
    virtual void setInitStateFromFinal(bool setInitInit = false);


    //! Set the internal initial intermediate and initial global state from the internal final_final state
    /*!
     *  (virtual function)
     *
     *  Set the intial  and init_int state and the final_final from the final state. We also can set the init_init state from this
     *  routine as well.
     *
     */
//Can protect
    virtual void setInitInitStateFromFinalFinal();

    //! Set the internal final intermediate and from the internal init state
    /*!
     *  (non-virtual function)  -> function should onionize in-first.
     *
     *  Set the final state from the init state. This is commonly called during a failed time step
     */
//Can protect
    void setFinalStateFromInit_Oin();

    //! Set the internal final intermediate state from the internal init state
    /*!
     *  (virtual function from Electrode)
     *
     *  Set the final state from the init state. This is commonly called during a failed time step
     */
//Can protect
    virtual void setFinalStateFromInit();

    //! Set the internal initial intermediate from the internal initial global state
    /*!
     *  Set the intial state from the init init state. We also can set the final state from this
     *  routine as well.
     *
     *  The final_final is not touched.
     *
     * @param setFinal   Boolean indicating whether you should set the final as well
     */
//Can protect
    virtual void setInitStateFromInitInit(bool setFinal = false);

    //! Set the internal final global state from the internal final intermediate state
    /*!
     *  (virtual function from Electrode)
     *
     *  Set the final_final state from the final state. This is commonly called at the end of successful base integration
     */
//Can protect
    virtual void setFinalFinalStateFromFinal();

    //! Set the internal final global state from the internal final intermediate state
    /*!
     *  (non-virtual function from Electrode)
     *
     *  Set the final_final state from the final state. This is commonly called at the end of successful base integration.
     *  This is an onionize in function
     */
//Can protect
    void setFinalFinalStateFromFinal_Oin();

    //! Calculates the change in the surface area of all external and internal interfaces within the electrode
    /*!
     *  (virtual)
     *  variables to be potentially altered
     *   surfaceAreaRS_[];
     *   isExternalSurface[]
     *   numExternalInterfacialSurfaces_;
     */
//Can protect
// deprecate
    virtual double calcSurfaceAreaChange(double deltaT);


    //! Recalulate the electrolyte mole number that fills up the designated space in the electrode
    //! and rebaseline the electrolyte mole numbers
    /*!
     *  We calculate the volume taken up by the electrolyte.
     *  Then we calculate the mole number for that volume
     */
    virtual double updateElectrolytePseudoMoles();

//Can protect, Should always appear externally as if this is off, can keep these if it's important to keep it on
// temporarily during the solve
    virtual void turnOffFollowElectrolyteMoles();
    virtual void turnOnFollowElectrolyteMoles();

    //! Get the ID of the external surface, if there is one
    const std::vector<bool>& getExternalSurfaceBooleans() const;


    //! Return the number of moles in a phase
    /*!
     * @param iph  Phase id

     * @return Returns the number of moles in kmol
     */
// I think a lot of the following functions are only used for I/O purposes, if so we should probably protect
// them and do all I/O through the friend class
    double phaseMoles(int iph) const;

    //! Returns the number of moles of an element
    /*!
     *  @param  ie   the index of the element in the global list
     *
     * @return Returns the number of moles in kmol
     */
    double elementMoles(int ie) const;


    //! Returns the number of moles of an element
    /*!
     *  @param  eName String Name of the element - two characters with the first capitalized
     *
     * @return Returns the number of moles in kmol
     */
    double elementMoles(std::string eName) const;

    //! Returns the number of moles of an element not including the electrolyte
    /*!
     *  @param  ie   the index of the element in the global list
     *
     * @return Returns the number of moles in kmol
     */
    double elementSolidMoles(int ie) const;


    //! Returns the number of moles of an element not including the electrolyte
    /*!
     *  @param  eName String Name of the element - two characters with the first capitalized
     *
     * @return Returns the number of moles in kmol
     */
    double elementSolidMoles(std::string eName) const;

    //! Returns the electrochemical potential of a single species
    /*!
     *     @param[in] Value of the global species index of the species
     *     @return Returns the electrochemical potential (J / kmol)
     */
    double speciesElectrochemPotential(int globalSpeciesIndex) const;

    double speciesChemPotential(int globalSpeciesIndex) const;

    //! Get mole fractions of all species in the phase object
    /*!
     *   @param x  Vector of mole fractions. Index is the same as PhaseList's global species index
     *             Length = number of species in PhaseList
     */
    void getMoleFractions(doublereal* const x) const;

    //! Return the mole fraction of a single species
    /*!
     *  @param globalSpeciesIndex  PhaseList's global species index
     */
    double moleFraction(int globalSpeciesIndex) const;

    //! Get mole numbers of all species in the phase object
    /*!
     *   @param x  Vector of mole numbers. Index is the same as PhaseList's global species index
     *             Length = number of species in PhaseList
     *             units are kmol
     */
    void getMoleNumSpecies(doublereal* const n) const;
    const std::vector<doublereal> & getMoleNumSpecies() const;

    //! Get mole numbers of all phases in the phase object
    /*!
     *   @param x  Vector of mole numbers. Index is the same as PhaseList  index
     *             Length = number of phases in PhaseList
     *             units are kmol
     */
    void getMoleNumPhases(doublereal* const np) const;

    //! Return the mole number of a single species
    /*!
     *  @param globalSpeciesIndex  PhaseList's global species index
     *
     *  @return Returns the mole number. Units are kmol.
     */
    double moleNumSpecies(int globalSpeciesIndex) const;


    //! Returns the index of a phase in the ReactionSurfaceDomain object
    //! given the index of that phase in the PhaseList object
    /*!
     * @param isk  Surface phase index, used to look up the ReactingSurfaceDomain object
     * @param PLph index of the phase in the PhaseList object, which is also the
     *             Electrode_Model object.
     *
     *  @return  Returns the index of the phase in the current ReactingSurDomain
     *           object. A value of -1 in this slot means that the phase doesn't
     *           participate in the  current ReactingSurDomain object
     */
// Deprecate
    int ReactingSurfacePhaseIndex(int isk, int PLph) const;


    //! Take the state (i.e., the final state) within the Electrode object and push it down
    //! to the ThermoPhase objects and other variables that are part of the Electrode object
    /*!
     *  (virtual function from Electrode - serial virtual)
     *  This virtual function is written so that child routines are called by parent routines.
     *  It may be the case that child routines will surplant parent routines, in order to create efficiency.
     *
     *  We take the values of spMoles_final_[], and the number of particles, and their size,
     *  which are the default specification of the state variables within the Electrode object,
     *  and propagate them down to the ThermoPhase objects in the electrode.
     *  We also calculate the volumetric properties of the Electrode, the phase moles,
     *  and the mole fractions, and the external radius of the particle
     *
     *  All of these properties are defined for the _final_ state.
     *
     *  Thus, this is the main routine that reconciles all of the state information within the object.
     *  At the end of this routine, all aspects of the final state are consistent with each other.
     *
     *  prerequisites: The object must have been already been created.
     *
     *
     *  Summary: State Variables
     *              spMoles_final_[kGlobal]  : kGlobal in solid phase species
     *              particleNumberToFollow_
     *
     *           Independent Variables
     *              spMoles_final_[kGlobal]  :   kGlobal in electrolyte phase
     *              phaseVoltages_[iph] :        iph in solid phase electrodes
     *              Temperature
     *              Pressure
     *
     *   Dependent StateVariables
     *              phaseMoles_final_[iph]       iph in solid phase electrode
     *              spMf_final_[kGlobal]         kGlobal in solid phase species
     *              VolPM_[kGlobal]              kGlobal in solid phase species
     *              spElectroChemPot_[kGlobal]   kGlobal in solid phase species
     *              ThermoPhase[iph]             iph in solid phase electrode
     *        	    phaseMolarVolumes_[iph]      iph in solid phase electrode
     *              ElectrodeSolidVolume_        (this is called from SolidVolume();
     *              Radius_exterior_final_       (this is calculated from ElectrodeSolidVolume_)
     *
     */
//Should probably protect
    virtual void updateState();

    //! Take the state (i.e., the final state) within the Electrode object and push it down
    //! to the ThermoPhase objects and other variables that are part of the Electrode object
    /*!
     *  This function is written so that child routines are not called by parent routines.
     *  It may be the case that child routines will surplant parent routines, in order to create efficiency.
     *
     *  We take the values of spMoles_final_[], and the number of particles, and their size,
     *  which are the default specification of the state variables within the Electrode object,
     *  and propagate them down to the ThermoPhase objects in the electrode.
     *  We also calculate the volumetric properties of the Electrode, the phase moles,
     *  and the mole fractions, and the external radius of the particle
     *
     *  All of these properties are defined for the _final_ state.
     *
     *  Thus, this is the main routine that reconciles all of the state information within the object.
     *  At the end of this routine, all aspects of the final state are consistent with each other.
     *
     *  prerequisites: The object must have been already been created.
     *
     *
     *  Summary: State Variables
     *              spMoles_final_[kGlobal]  : kGlobal in solid phase species
     *              particleNumberToFollow_
     *
     *           Independent Variables
     *              spMoles_final_[kGlobal]  :   kGlobal in electrolyte phase
     *              phaseVoltages_[iph] :        iph in solid phase electrodes
     *              Temperature
     *              Pressure
     *
     *   Dependent StateVariables
     *              phaseMoles_final_[iph]       iph in solid phase electrode
     *              spMf_final_[kGlobal]         kGlobal in solid phase species
     *              VolPM_[kGlobal]              kGlobal in solid phase species
     *              spElectroChemPot_[kGlobal]   kGlobal in solid phase species
     *              ThermoPhase[iph]             iph in solid phase electrode
     *        	    phaseMolarVolumes_[iph]      iph in solid phase electrode
     *              ElectrodeSolidVolume_        (this is called from SolidVolume();
     *              Radius_exterior_final_       (this is calculated from ElectrodeSolidVolume_)
     *
     */
//Can protect
    void updateState_OnionOut();

    //!  Recalculate the surface areas of the surfaces for the final state
    /*!
     *    (virtual function from Electrode)
     *
     *    Here, we don't know anything about the morphology of the particle. Models that do know will override this function.
     *    Hwere we assume that the surface area is equal to the exterior surface area multiplied by the numbers of particles.
     *    We also respect keeping the surface area as zero, if it is initially set to zero.
     *
     *    Dependent StateVariables Used
     *         Radius_exterior_final;
     *         particleNumberToFollow_
     *
     *    Dependent StateVariables Calculated
     *          surfaceAreaRS_final_[]
     */
//Can protect
    virtual void updateSurfaceAreas();

    //! This is used to set the phase information that is implicit but not set by a restart or an initialization
    /*!
     *  (virtual function from Electrode)
     *
     *  This is called immediately after the restart file's contents are loaded into the electrode object.
     *  We then call this function to calculate the internal flags. Then, we call updateState() to make sure
     *  all information about the state is self-consistent.
     *
     *  Extra information that may be needed in advance of a successful updateState() call that specifies all of the
     *  information in the state
     *
     *  @param flagErrors If true any changes in the current flags caused by a mismatch between the state
     *                    and the values of the flags will cause an error exit.
     *
     *  @return  Returns whether there has been any discovered errors (currently ignored)
     */
//Can protect
    virtual bool stateToPhaseFlagsReconciliation(bool flagErrors);

    //!  Reactant stoichiometric coefficient
    /*!
     *   Get the reactant stoichiometric coefficient for the kth global species
     *   in the ith reaction of the reacting surface domain with index isk.
     */
//Can protect
    double reactantStoichCoeff(const int isk, int kGlobal, int i) const;

    //! Product stoichiometric coefficient
    /*!
     * Get the product stoichiometric coefficient for the kth global species
     * in the ith reaction of the reacting surface domain with index isk.
     */
//Can protect
    double productStoichCoeff(const int isk, int kGlobal, int i) const;

    //! Returns the current potential drop across the electrode
    /*!
     *   This is equal to phiMetal - phiSoln
     */
    double potentialDrop() const ;

    //! Set the phase existence flag in the electrode kinetics object so that kinetics
    //! are calculated correctly
    /*!
     *    Flags are set in the kinetics object to tell the kinetics object which phases
     *    have zero moles.  The zero mole indicator is taken from phaseMoles_final_[]. Therefore,
     *    the final state is queried.
     *    There is a special case. Phases that are in justBornPhase_[] vector are allowed to
     *    be set to exist even if their phase moles are zero.
     *
     * @param assumeStableSingleSpeciesPhases Assume that single phases are stable. This
     *                         allows their production rates to be calculated
     */
    virtual void setPhaseExistenceForReactingSurfaces(bool assumeStableSingleSpeciesPhases = false);

    //! Print conditions of the electrode for the current integration step to stdout
    /*!
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values
     */
    virtual void printElectrode(int pSrc = 1, bool subTimeStep = false);

    //! Print Capacity Utilization conditions of the electrode for the current integration step to stdout
    /*!
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values
     */
    virtual void printElectrodeCapacityInfo(int pSrc, bool subTimeStep);

    //! Print a table of values for each phase in the problem  for the current integration step to stdout
    /*!
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values
     */
    virtual void printElectrodePhaseList(int pSrc = 1, bool subTimeStep = false);

    //! Print condition of a phase in the electrode
    /*!
     *  @param iPhase        Print the phase
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values
     */
    virtual void printElectrodePhase(int iPhase, int pSrc = 1,  bool subTimeStep = false);


    //! Return the number of extra print tables
// Deprecate
    virtual int getNumPrintTables() const;

    //! Get the values that are printed in tables for the 1D code.
    /*!
     *   @param itable    table id
     *   @param colNames   string names of the header (length is the length of the column)
     *   @param colValues    Value of the columns (length is the length of the column)
     */
// Deprecate? this doesn't appear to be used
    virtual void getPrintTable(int itable, std::vector<std::string>& colNames,
                               std::vector<double>& colValues) const;

    //! Toggle switch for Printing for the predictor corrector
    /*!
     *  To get printing from the predictor-corrector printing routines, both the general
     *  printing level must be high enough and this static routine should be turned on.
     *  Also the ifdef DEBUG_PREDICTOR should be turned on if printing is desired at a
     *  low level of printLvl_
     *
     *   0 No printing 
     *   1 predictorCorrect printing is turned on (default)
     */
    static int s_printLvl_PREDICTOR_CORRECTOR;

    /********************************************************************************************************************
     *  OPEN CIRCUIT VOLTAGE
     *******************************************************************************************************************/

    //! Returns the standard state E0 for the electrode based on a single reaction (virtual)
    /*!
     *  When there is more than a single reaction,
     *  pick open circuit potential for reaction that is
     *  closest to equilibrium given the cell voltage since this one
     *  is the one for which open circuit is most relevant.
     *  Note that it will be possible for the standard state OCV
     *  to be computed for a different reaction relative to the
     *  method openCircuitVoltage(isk) that computes OCV from
     *  the current concentrations because it looks for the reaction
     *  closest to equilibrium given the current cell voltage.
     *
     *
     *  (virtual function from Electrode.h)
     *
     * @param isk  Reacting surface domain id
     *  @param iReaction   Explicit index of the reaction. If -1, then it attempts
     *                     to pick the reaction that best represents the open circuit potential.
     */
    virtual double openCircuitVoltageSSRxn(int isk, int iReaction = -1) const;

    //! Returns the equilibrium OCV for the selected ReactingSurfaceDomain, current conditions
    //! based on a single reaction (virtual)
    /*!
     *  When there is more than a single reaction,
     *  pick open circuit potential for a reaction that is
     *  closest to equilibrium given the cell voltage, since this one
     *  is the one for which open circuit is most relevant.
     *
     *   @param isk         Reacting surface domain id
     *   @param iReaction   Explicit index of the reaction. If -1, then it attempts
     *                      to pick the reaction that best represents the open circuit potential.
     */
    virtual double openCircuitVoltageRxn(int isk, int iReaction = -1) const;

    //! Returns the equilibrium OCV for the selected ReactingSurfaceDomain and current conditions (virtual)
    /*!
     *  This routine uses a root finder to find the voltage at which there
     *  is zero net electron production.  It leaves the object unchanged. However, it
     *  does change the voltage of the phases during the calculation, so this is a non const function.
     *
     * @param isk  Reacting surface domain id
     */
    virtual double openCircuitVoltage(int isk);

    //! Returns the vector of OCV's for the selected ReactingSurfaceDomain and current conditions.
    /*
     * @param isk  Reacting surface domain id
     */
    void  getOpenCircuitVoltages(int isk, double* ocv) const;

    //! Returns the exchange current density for a given surface reaction in A/m^2
    /*!
     *  The parameterization for a single reaction defines the exchange
     *  current density formulation for that reaction. Note the exchange current density depends
     *  on the activity concentrations of the reactants and the products.
     *
     *         i  = io [ exp( (Betaf  Nu  nSt  F)/ RT) - exp(- ((1 - Betaf) Nu nSt F)/ RT) ]
     *
     *              where Nu = E - OCV is the overpotential
     *
     *              and E = Phi(Metal) - Phi(Soln)
     *
     *   where the current is defined as
     *
     *         i = - z_e- * F * d[e-]/dt
     *
     *    and
     *
     *         d[e-]/dt = (n) ROP_net
     *
     *   A note about signs
     *
     *    i is defined here as being positive if there is a net current going from the metal into solution.
     *    This in turn means that if there is a positive generation of electrons, then the current
     *    is positive.
     *
     *    In this formulation the stoichiometric number of electrons can be positive or negative
     *    It will be negative if the electrons appears on the reactant side of the equation.
     *
     *    In this formulation, io can be positive or negative. It's sign will be determined by the
     *    sign of the nSt, the stoichiometric electrons.
     *
     *    This routine returns io in the above discusion.
     *
     *  @param isk   Reacting surface index.
     *  @param irxn  Reaction number on the reacting surface
     *
     *  @return  returns the exchange current in units of amps / m2
     */
// Deprecate
    double  getExchangeCurrentDensity(int isk, int irxn) const;

    //! Returns the overpotential for the current conditions
    /*!
     *  The overpotential is the current voltage minus the open circuit voltage.
     *  This routine calculates the open circuit voltage for the surface via a rootfinder.
     *
     *   @param isk   reacting surface number
     */
    double overpotential(int isk);

    //! Returns the overpotential for the current conditions based on a particular reaction
    /*!
     *  The overpotential is the current voltage minus the open circuit voltage.
     *  This routine calculates the open circuit voltage for the surface by calling a
     *  openCircuitVoltageRxn() routine.
     *
     *   @param isk   reacting surface number
     *   @param irxn  Reaction number
     */
    double overpotentialRxn(int isk, int irxn = -1);

    //! Return the kinetics species index of the electron for the surface phase, isph
    //! Return the kinetics species index of the electron for the surface phase, isph
    /*!
     *  @param isph   Surface phase index
     */
    int kKinSpecElectron(const int isph) const;

    //!  Return the global species index of the electron in the Electrode object
    int kSpecElectron() const;

    //! Return the index of the phase corresponding to the currently active metal
    /*!
     *  Note this may change during the simulation
     */
    int metalPhaseIndex() const;

    //! return the index of the phase corresponding to the soln
    int solnPhaseIndex() const;

    virtual int numSolnPhaseSpecies() const;

    // -----------------------------------------------------------------------------------------------------------------
    // --------------------------- CAPACITY CALCULATION OUTPUT  --------------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------

    //  Capacity versus DepthOfDischarge
    //   The difference occurs when there is a loss of capacity due to irreversible changes in the
    //   electrode.
    //   Capacity calculations are current at all times.
    //   Depth of discharge calculations are counters for electrons gained/lost in the electrode.
    //   Therefore, during a cycle the DoD counter may not go back to zero. The Capacity counter
    //   will register zero.
    //   So far it doesn't matter as we don't have any loss mechanisms.

    //! Returns the total capacity of the electrode in Amp seconds
    /*!
     *  Returns the capacity of the electrode in Amps seconds.
     *  This is the same as the number of coulombs that can be delivered at any voltage.
     *  Note, this number differs from the capacity of electrodes that is usually quoted for
     *  a battery. That number depends on the rate of discharge and also depends on the
     *  specification of a cutoff voltage. Here, we dispense with both of these specifications.
     *  So, it should be considered a theoretical capacity at zero current and minimal cutoff voltage
     *  considering the current state of the battery. The initial theoretical capacity given
     *  ideal conditions is given by capacityInitial().
     *
     *  It will also include all plateaus that are defined by the electrode object.
     *
     *  This capacity may change as degradation mechanisms cause the electrode to lose capability.
     *  Therefore, the capacity will be a function of time.
     *  At all times the following relation holds:
     *
     *  capacity() = capacityDischarged() + capacityLeft().
     *
     *  @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     *                   If positive or zero, each plateau is treated as a separate entity.
     *
     *  @return returns the theoretical capacity of the electrode in Amp seconds = coulombs.
     */
    virtual double capacity(int platNum = -1) const;

    //! Initial capacity of the electrode in Amp seconds
    /*!
     *  This is the initial capacity of the electrode before any degradation occurs.
     *
     *  At all times the following relation holds:
     *
     *  capacity() = capacityDischarged() + capacityLeft().
     *
     *  @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     *                   If positive or zero, each plateau is treated as a separate entity.
     */
    virtual double capacityInitial(int platNum = -1) const;

    //! Amount of charge that the electrode has discharged up to this point (coulombs)
    /*!
     *   We report the number in terms of Amp seconds = coulombs.
     *   Note the capacity discharged for cathodes will be defined as the negative of the electron
     *   source term, as this refers to the forward discharge of a cathode.
     *   This definition is necessary for the following to be true.
     *
     *         capacity() = capacityDischarged() + capacityLeft()
     *
     *   Note, the current is defined as the amount
     *   of positive charge that goes from the solid into the electrolyte.
     *
     *   @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     *                   If positive or zero, each plateau is treated as a separate entity.
     */
    virtual double capacityDischarged(int platNum = -1) const;

    //! Reset the counters that keep track of the amount of discharge to date
    virtual void resetCapacityDischargedToDate();

    //! Amount of charge that the electrode that has available to be discharged
    /*!
     *  We report the number in terms of Amp seconds = coulombs. This accounts for loss mechanisms.
     *
     *        At all times the following relation holds:
     *
     *  capacity() = capacityDischarged() + capacityLeft() + capacityStarting().
     *
     *   @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     *                   If positive or zero, each plateau is treated as a separate entity.
     */
    virtual double capacityLeft(int platNum = -1, double voltsMax = 50.0, double voltsMin = -50.0) const;


    //! Initial starting depth of discharge in coulombs
    /*!
     * 
     *
     *   @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     *                   If positive or zero, each plateau is treated as a separate entity.
     */
    virtual double depthOfDischargeStarting(int platNum = -1) const;

    //! Report the current depth of discharge in Amp seconds
    /*!
     * Report the current depth of discharge. This is roughly equal to the total
     * number of electrons that has been theoretically discharged from a fully charged state.
     *  For multiple cycles, this becomes the true electron counter for the electrode.
     *
     * Usually this is reported as a function of the discharge rate and there is a
     * cutoff voltage at which the electron counting is turned off. Neither of these
     * concepts is employed here.
     *
     *  @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     *
     *  @return  returns the depth of discharge in Amp seconds
     */
    double depthOfDischarge(int platNum = -1) const;

    //! Report the current depth of discharge as a fraction of the total capacity
    /*!
     * Report the current depth of discharge. This is roughly equal to the total
     * number of electrons that has been theoretically discharged from a fully charged state compared
     * to the number of electrons that can theoretically be discharged.
     *
     * Usually this is reported as a function of the discharge rate and there is a
     * cutoff voltage at which the electron counting is turned off. Neither of these
     * concepts is employed here.
     *
     *  @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     *
     *  @return  returns the depth of discharge as a fraction of the total possible discharged
     */
    double depthOfDischargeFraction(int platNum = -1) const;

    //! Report the current depth of discharge divided by the molar size of the electrode
    /*!
     * Report the current depth of discharge divided by the molar size of the electrode.
     *
     * Usually this is reported as a function of the discharge rate and there is a
     * cutoff voltage at which the electron counting is turned off. Neither of these
     * concepts is employed here.
     *
     *  @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     *
     *  @return  returns the depth of discharge in percent
     */
    double depthOfDischargePerMole(int platNum = -1) const;

    bool checkCapacityBalances_final(int platNum = -1) const;

    //! Experimental routine to enforce a balance on the electrode capacity given previous balance information
    /*! 
     *  This routine equalizes the capacity 
     */
    virtual void fixCapacityBalances_final();


    // -----------------------------------------------------------------------------------------------------------------
    // --------------------------  CAPACITY CALCULATION SETUP ---- -----------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------
    //       Definitions are in the file Electrode_Capacity.cpp


    //! Returns the capacity type of the electrode
    /*!
     *  @return returns whether the electrode capacity calculation is designated a
     *                  anode or a cathode type
     */
    Electrode_Capacity_Type_Enum  capacityType() const;

    //! Sets the capacity type of the electrode
    /*!
     *  @param electrodeCapacityType  Sets the electrode capacity calculation is designated a
     *                                anode or a cathode type
     */
    void setCapacityType(Electrode_Capacity_Type_Enum electrodeCapacityType);

    //! Set the relative current capacity discharged per mole
    /*!
     * (virtual function from Electrode)
     * ( this is a serial virtual function - not to be used in creation)
     *
     * This is roughly equal to the total number of electrons that has been discharged
     * from a fully charged state divided by the total moles of solid species in the electrode
     *
     *  @param relDischargedperMole      Relative value of the discharge per mole. Always goes between 0 and number of electrons
     *                                   per active mole, num
     *                                  0 means that the electrode is fully charged, num means that it is fully discharged.
     *
     *  @param platNum            Plateau number. Default is -1 which treats all plateaus as a single entity and
     *                            the relative discharged as a single combined fraction. If platNum is
     *                            >= 0, then the discharge is relative to the current plateau.
     */
    virtual void setRelativeCapacityDischargedPerMole(double relDischargedPerMole, int platNum = -1);

    //! Set parameters that tell the object how to calculate the capacity of the electrode
    /*!
     * @param sName          Name of the species that contains the capacity
     * @param coeffLeft      Coefficient describing how many electrons are left
     * @param coeffZeroDoD   Coefficient describing how many electrons are originally there
     */
    void setCapacityCalcParams(std::string sName, double coeffLeft, double coeffZeroDoD);

    //! Set the Capacity coefficients for the LiSi anode system
    /*!
     *  @todo Move or get rid of
     */
    void setCapacityCoeff_LiSi() const;

    //! Set the Capacity coefficients for the FeS2 cathode system
    /*!
     *  @todo Move or get rid of
     */
    void setCapacityCoeff_FeS2() const;


    //! Set the Capacity coefficients for the LiSi anode system
    /*!
     *  @todo Move or get rid of
     */
    void setCapacityCoeff_LiSi_Li() const;


    //! Set the Capacity coefficients for the modified FeS2 cathode system
    /*!
     *  @todo Move or get rid of
     */
    void setCapacityCoeff_FeS2_Combo() const;

    //!  Capacity Parameters for the MCMB anode
    /*!
     *  Here we supply parameters for the arrays
     *          capacityLeftSpeciesCoeff_[]
     *          capacityZeroDoDSpeciesCoeff_[]
     *  This is an intercalated electrode, named MCMB_Interstitials_anode, with the following species
     *         Li_C6-bulk
     *          V_C6-bulk
     *
     *    @todo Move or get rid of
     */
    void setCapacityCoeff_MCMB() const;

    //! Set the capacity coefficients from the input file
    /*!
     *  Here we set :
     *
     *  These are set from the input file.
     *
     *      capacityLeftSpeciesCoeff_[iGlobSpeciesIndex
     *	    capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex]
     *
     *  It is a fatal error to call this member function if the coefficients have
     *  not been set from the input file.
     *
     *  @param ei   Input Key file object
     */
    void setCapacityCoeffFromInput(const ELECTRODE_KEY_INPUT* const ei);


    // --------------            EXTRA GLOBAL RXN PATHWAYS  -----------------------------

    //! Add a global reaction object to the internal list
    /*!
     *  These global reaction pathways are made up of a linear combination of
     *  existing reactions.
     *
     *     @param egr_ptr       Input struct describing the reaction
     */
    void addExtraGlobalRxn(EGRInput* egr_ptr);

    //! Get the RxnMolChange pointer object for a single extra global reaciton
    /*!
     *  int iegr  Index of the extra global reaction
     */
    RxnMolChange*   rxnMolChangesEGR(int iegr);

    //!  return the extra global rxn pathway
    /*!
     *  @param iegr  index of the pathway
     */
    Cantera::ExtraGlobalRxn* extraGlobalRxnPathway(int iegr);


    int processExtraGlobalRxnPathways();

    int numExtraGlobalRxnPathways() const;

// Deprecate
    void updatePhaseNumbersTmp(std::vector<doublereal>& spMoles_tmp,
                               std::vector<doublereal>& phaseMoles_tmp, std::vector<doublereal>& spMf_tmp);

    //-----------------------------------------------------------------------------------------------------------------
    // -----------------------  STATE and PRINTING FUNCTIONS ----------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------

    //! Determines the level of printing for each step.
    /*!
     *   0 -> absolutely nothing is printed for a single call to integrate.
     *   1 -> One line summary per integrate call
     *   2 -> short description, points of interest: Table of nonlinear solve - one line per iteration
     *   3 -> Table is included -> More printing per nonlinear iteration (default) that occurs during the table
     *   4 -> Summaries of the nonlinear solve iteration as they are occurring -> table no longer printed
     *   5 -> Algorithm information on the nonlinear iterates are printed out
     *   6 -> Additional info on the nonlinear iterates are printed out
     *   7 -> Additional info on the linear solve is printed out.
     *   8 -> Info on a per iterate of the linear solve is printed out.
     */
    void setPrintLevel(int printLvl);


    //! Write out CSV tabular data on the integrations
    /*!
     *  This routine write out a line of tabular data for the current conditions.
     *  The idea is to print out tabular data about each intermediate run and about each
     *  global run.
     *
     *  Whether this is called from the routine integrate() or not is determined by the variable
     *  printCSVLvl_.  
     * 
     *  This routine is used extensively in the tests, as it prints out a line for each time step.
     *
     *  @param itype Type of the data
     *            - 0      Initialization information
     *            - 1      Normal intermediate information
     *            - 2      Normal intermediate information at the end of a global time step
     *            - -1     Failed intermediate calculation, a failure from a nonlinear solver step
     *            - -2     Failed calculation from the predictor step - not necessarily significant.
     */
    void writeCSVData(int itype);

    //! Make a XML data tree on the integrations
    /*!
     * Create an XML data tree containing the timeIncrement XML nodes containing the state of the electrodes,
     * and the time history of the integration.
     *
     *  @param itypeXML Type of the data
     *            - 0      Initialization information
     *            - 1      Information for the current global step that will be used in the time advance
     *            - 2      Information for the current global step as it stands now
     *            - 3      Write the timeState t_final information only
     *            - 4      Write timeIncrement information out as if intermediate steps were global steps.
     */
    Cantera::XML_Node* makeXMLTimeIncrementDataTree(int itypeXML);

    //!  Create a timeIncrement XML element to store the results for intermediate steps of the solver.
    /*!
     *    Creates the following XML tree structure.
     *    If the boolean addInitState is false the timeState t_init is not written.
     *
     *     <timeIncrement   type="intermediate">
     *        <timeState type="t_init">
     *              ....xmlStateData_init_
     *        </timeState>
     *        <timeState type="t_intermediate">
     *              ....xmlStateData_final_
     *        </timeState>
     *     </timeIncrement>
     */
    void makeXML_TI_intermediate(bool addInitState = true);

    //! Adds to a timeIncrement XML element to store the results for intermediate or global-final steps of the solver.
    /*!
     *    Adds a timeState record to the following XML tree structure.
     *
     *     <timeIncrement   type="intermediate">
     *        <timeState type="t_init">
     *              ....xmlStateData_init_
     *        </timeState>
     *        <timeState type="t_intermediate">
     *              ....xmlStateData_final_
     *        </timeState>
     *        <timeState type="t_inal">
     *              ....xmlStateData_final_
     *        </timeState>
     *     </timeIncrement>
     *
     *    This XML Tree is storred in the variable  xmlTimeIncrementData_
     *
     *  @param notDone   Boolean that if true sets the type attribute to t_intermediate
     *                   If false, the type attribute is set to t_final
     */
    void addtoXML_TI_final(bool notDone);

    //! Creates a timeIncrement XML element to store the results for global steps of the Electrode solver
    /*!
     *   Creates a XML Tree structure of the following form. What this function does is to delete the old record
     *   and start a new record. Later calls to addtoXML_TI_final() adds records to the XML tree.
     *
     *     <timeIncrement   type="global">
     *        <timeState type="t_init">
     *              ....xmlStateData_init_
     *        </timeState>
     *        <timeState type="t_intermediate">
     *              ....xmlStateData_final_
     *        </timeState>
     *        <timeState type="t_inal">
     *              ....xmlStateData_final_
     *        </timeState>
     *     </timeIncrement>
     *
     *    This XML Tree is storred in the variable  xmlTimeIncrementData_
     *
     *
     *  @param addInitState   Boolean that if true adds the initial state to the tree.
     *                        The default is true.
     */
    void startXML_TI_final(bool addInitState = true);

    //! Specifies the amount of output that the Electrode object writes to its XML solution file
    /*!
     *    The level is given by the following table
     *            - 0     Zero output (default)
     *            - 1     Global time end
     *            - 2     global time increment
     *            - 3     intermediate steps
     *
     *    @param level Specifies the level of output
     *    @baseName    Basename of the output file. The domain and cell number are tacked on.
     */
    void specifySolutionFileLevel(int level, const char* baseName);

    //!  Wrap the Time increment XML element with a solution XML element and then write it out to
    //!  an output file
    /*!
     *       We assume that the XML file has the following topology
     *
     *
     *         <ctml>
     *           <electrodeOutput index = 1>
     *             <timeIncrement index = 1>
     *
     *             </timeIncrement index>
     *           </electrodeOutput>
     *         </ctml>
     *
     *  This routine adds another <timeIncrement> XML element to the end of the file. It
     *  takes care to first eliminate any existing to backspace over the last
     *  </electrodeOutput> and </ctml> entries before writing the new <timeIncrement> XML
     *  element.
     */
    void writeSolutionTimeIncrement();

    //!  Write the state of the Electrode object at the t_final time out to the output XML_Node
    /*!
     *  This routine is used by the 1Dsolvers to write a restart file for the electrode object out
     *  to an XML file. The following XML records are written out by this routine, usually into 
     *  a surrounding domain XML_Node
     *
     *    <domain id="BulkDomain1D_0" numVariables="6" points="10" type="bulk">
     *       <TimeState Cell Number="0" Domain="0" type="t_final">   <------------------------ Routines writes this out
     *          <time> 1e-08 </time>
     *          <electrodeState>
     *              . . .
     *          </electrodeState>
     *       </TimeState>
     *       . . .
     *    </domain>
     *
     *   The electrodeState record is written out by the child Electrode objects from saved data.
     *   Usually, all objects within the domain write their own records.
     *
     *   @param bb  Output XML mode
     *
     *   @return Returns whether the step was successful or not. Note, the electrode object doesn't
     *           necessarily create the XML information to be saved, in which case this routine
     *           will return false.
     */
    bool writeTimeStateFinal_toXML(XML_Node&  bb);

    //!  Select the globa time step increment record by the consequatively numbered record index number
    /*!
     *    @param   xSoln               Solution file for the electrode object
     *    @param   globalTimeStepNum   Time step number to select
     *
     *  @TODO shouldn't this be a static routine.
     */
    XML_Node* selectGlobalTimeStepIncrement(XML_Node* xSoln, int globalTimeStepNum);

    //! Given a Time increment record this routine loads the saved solution for t_final into the electrode object
    /*!
     *
     *     This routine reads an XML tree layout, such as the example below. It is given the globalTimeStep XML element
     *     pointer. It then reads the timeState  final record and final time, and initializes the current object with
     *     the state contained in the record.
     *
     *            <globalTimeStep index = 312>
     *               <timeIncrement type="global">
     *                  <timeState type="t_init">
     *                     ....xmlStateData_init_
     *                  </timeState>
     *                  <timeState type="t_intermediate">
     *                      ....xmlStateData_final_
     *                 </timeState>
     *                 <timeState type="t_final">
     *                    <time>   3.45 <time>             <-----------  time used.
     *                    <electrodeState>                 <-----------  Record read
     *                    </electrodeState>
     *                 </timeState>
     *               </timeIncrement>
     *             </globalTimeStep>
     *
     *
     *  @param xGSTI    Global time step increment record
     */
    void loadGlobalTimeStepTFinalState(XML_Node* xGTSI);


    double loadTimeStateFinal(XML_Node& xTimeStateFinal);

    //------------------------------------------------------------------------------------------------------------------

protected:

    //! Capacity type of the electrode
    /*!
     *  Either CAPACITY_ANODE_ECT or CAPACITY_CATHODE_ECT.
     *
     *  This determines what it direction it means to have positive capacity and
     *  the direction of the state of charge as well.
     */
    Electrode_Capacity_Type_Enum  electrodeCapacityType_;

    //! true if there is a pending integration step
    /*!
     *  We keep track of whether there is a pending integration step with this variable
     */
    int pendingIntegratedStep_;

    //! Multiphase object
    /*!
     *  This includes all the phases which are equilibrated during the
     *  initialization of the program
     */
// Deprecate?
    //Cantera::MultiPhase* m_mp;
    //Cantera::MultiPhase* m_mpTemp;


    //! Integer representing the Problem type.
    /*!
     *  The identity of  what is held constant. Currently,
     *   T and P are held constant, and this input is ignored
     *   0 : we are tracking mole numbers of electrolyte
     *   1 - We are not tracking mole numbers of electrolyte
     *       Basically, mole number tracking stops at the electrolyte
     *       electrode interface.
     */
// Deprecate
    int prob_type;

    //! Number of surfaces in the problem
    /*!
     *   The number of surfaces is initially identified with the number of surface phases
     *   However, this can change. Each surface will have a surface area associated with it
     *   and it may have an active Kinetics object associated with it.
     */
    int numSurfaces_;

    //!  Integer representing whehter we keep a conserved quantity for the electrolyte phase mole numbers or not 
    /*!
     *  0:  Constant electrolyte composition but varying amount of electrolyte,
     *      specified porosity
     *      During a subintegration step the electrode can change volume by changing the solidVolume.
     *      The electrolyte moles are not tracked so that the porosity of the electrode is kept constant
     *      or at a specified functional value. Electrolyte flows into or out of the domain
     *      in order to maintain the specified porosity. 
     *      (default implementation)
     *
     *  1:  Tracked electrolyte composition and amount of electrolyte, 
     *      calculated porosity
     *      We are tracking the mole numbers of electrolyte. The extra flux of
     *      of ions needed to maintain neutrality is preserved as well. We must maintain electroneutrality
     *      within the phase at all times.
     *      The mole numbers of electrolyte may vary if there are reactions creating
     *      neutrals or compensating ions. The volume of electrolyte, electrode both
     *      varies, so the total volume of the electrode varies. The porosity changes
     *      according to the relative amount of solid to electrolyte
     *
     *      Note constant volume or constant pressure are additional constraints that must
     *      be specified elsewhere to specify the solid mechanics conditions. This constraint has
     *      to do with whether the electrolyte flows in and out of the electrolyte domain to 
     *      satisfy the porosity constraints. The constant volume or pressure constraint is
     *      a separate degree of freedom.
     */
    int followElectrolyteMoles_;

    //! If this is negative, then we are tracking the mole numbers of electrolyte species
    /*!
     *  If this is positive, then we multiple this number by the
     *  mole fraction vector to determine the number of moles of
     *  electrolyte species.
     *  @deprecate This is confusing. We don't need it. Just use spMoleNumber(solnPhase).
     *              HKM - concur 2/6/14
     */
    double electrolytePseudoMoles_;

    //! Type of interpolation of the external fields within the object
    /*!
     *  The default is to use Backwards Euler interpolation using the t_final field values.
     *  This is called T_FINAL_CONST_FIS
     *
     *  -> To be expanded later into a workable concept
     */
    enum Electrode_Exterior_Field_Interpolation_Scheme_Enum externalFieldInterpolationType_;

    //! Initial value of the time, at the start of the current global time step
    double t_init_init_;

    //! Final value of the time, at the end of the current global time step
    double t_final_final_;

    //! Initial value of the time, at the start of the current intermediate time step
    double tinit_;

    //! Final value of the time, at the end of the current intermediate time step
    double tfinal_;

    //! Maximum Value of deltaT used in the subcycle integration at the start of each integration
    /*!
     *  The default for this is inf, so that the deltaTsubcycle_ value is chosen to be equal to deltaT.
     */
    double deltaTsubcycleMax_;

    //! Value of deltaT used in the sybcycle integration at the start of the current interval
    /*!
     *   This object is meant to do the same integration over the same time interval over and over again as part of a
     *   jacobian evaluation. In order to do this effectively, the time step history has to be matched as best
     *   as it can be. Therefore, we must start the integration with the same time step value. This
     *   is the storred value of the initial time step.
     */
    double deltaTsubcycle_init_init_;

    //! Current value of deltaT for the subcycle of the time integration
    double deltaTsubcycle_;

    //! Value of deltaT to be used for the next subcycle
    /*!
     *  The default for this is inf, so that the deltaTsubcycle chosen is equal to deltaT.
     */
    double deltaTsubcycleNext_;

    //! Value of deltaT used in the sybcycle integration at the start of the next global time step interval
    /*!
     *   We calculate the subcycle deltaT value to be used in next global time step interval
     */
    double deltaTsubcycle_init_next_;



public:
    //! Flag for the choice for the initial subcycle time step for a new global time interval
    /*!
     *   0  Use the value of the time step chosen by the last subcycle of the previous global time step (default)
     *   1  Use the value of the time step chosen by the first subcycle of the previous global time step
     *   2  The subgrid integration is set to 1/10 of the global time step
     */
    int choiceDeltaTsubcycle_init_;

    //! Number of subcyles taken on the last
    int numIntegrationSubCycles_final_final_;

    //! Boolean indicating whether we should be doing thermal property calculations during updateState()
    //! calculations.
    /*!
     *   This is public and can be changed externally
     */
    bool doThermalPropertyCalculations_;

protected:

    //! Temperature of the electrode (Kelvin)
    double temperature_;

    //! Pressure of the electrode (MKS - Pascal)
    double pressure_;

    //! Total volume of the electrode material
    /*!
     *  This is equal to the sum of VolPM[k] * spMole[k] for the solid phase
     *  components of the electrode. The volume of the electrolyte is not included.
     */
    double ElectrodeSolidVolume_;

    //! Vector of molar volumes of each of the phases in the mechanism
    /*!
     *  This has units of m**3 / kmol.
     *  It has length equal to the total number of phases (vol and surf) in the mechanism.
     *  Here we assign 0 to the molar volume for surface phases. This is the actual
     *  model used by Cantera.
     */
    std::vector<double> phaseMolarVolumes_;

    //! Vector of molar areas of each of the surface phases in the mechanism
    /*!
     *  This has units of m**2 / kmol.
     *  It has length equal to the total number of surface phases in the mechanism.
     *  We actually get this number by calling molarVolume() from the Cantera
     *  Surface ThermoPhase object.
     */
    std::vector<double> sphaseMolarAreas_;

    //! Partial molar volumes of all of the species species
    /*!
     *  Length = global number of species in PhaseList species vector
     *  units = m**3 kmol-1
     */
// Is this a duplicate of phaseMolarVolumes_?
    std::vector<double> VolPM_;

    //! Partial molar Heat Capacity at constant volume of all of the species
    /*!
     *  Length = global number of species in PhaseList species vector
     *  Units = Joules / Kelvin
     */
    mutable std::vector<double> CvPM_;

    //! Number of moles of each species in each phase at the end of each
    //! subcycle of the integration step
    /*!
     *  INDEPENDENT STATE VARIABLE FOR THE ELECTRODE PROBLEM
     *   Number of moles of each species in each phase at the end of each
     *   subcycle of the integration step.
     *
     *    length = Number of species in PhaseList
     *    indexing of PhaseList
     */
    std::vector<double> spMoles_final_;

    //! Number of moles of each species in each phase at the end of each
    //! outer loop of the integration step
    /*!
     *  INDEPENDENT STATE VARIABLE FOR THE ELECTRODE PROBLEM
     *   Number of moles of each species in each phase at the end of each
     *   subcycle of the integration step.
     *
     *    length = Number of species in PhaseList
     *    indexing of PhaseList
     */
    std::vector<double> spMoles_final_final_;

    //! Number of moles of each species in each phase at the start of each
    //! subcycle of the integration step.
    /*!
     *  INDEPENDENT STATE VARIABLE FOR THE ELECTRODE PROBLEM
     *   Number of moles of each species in each phase at the start of each
     *   subcycle of the integration step.
     *
     *    length = Number of species in PhaseList
     *    indexing of PhaseList
     */
    std::vector<double> spMoles_init_;

    //! Number of moles of each species in each phase at the start of the integration step
    /*!
     *  This will only be updated when we are moving onto the next outside time step.
     *  this is true of all entities entitled "_init_init_"
     *
     *    length = Number of species in PhaseList
     *    indexing of PhaseList
     */
    std::vector<double> spMoles_init_init_;

    //! Time derivative of the species mole numbers
    /*!
     *  This is the time derivative of the subcycle time step.
     *  It's calculated from spMoles_final_[] and spMoels_init_[].
     *  It may not be kept current.
     */
// Deprecate, appears unused
    std::vector<double> spMoles_dot_;

    //! Predicted Mole vector
    std::vector<double> spMoles_predict_;

//TODO: Do we need to track both spMoles and spMf?
    //! Mole fraction of each species in each phase at the start of the time step
    /*!
     *  DEPENDENT STATE VARIABLE for the electrode.
     *  This variable is kept synched with the spMoles_init_ variable at all times.
     *
     *    length = Number of species in PhaseList
     *    indexing of PhaseList
     */
    std::vector<double> spMf_init_;

    //! Mole fraction of each species in each phase at the start of the global time step
    /*!
     *  DEPENDENT STATE VARIABLE for the electrode.
     *  This variable is kept synched with the spMoles_init_init_ variable at all times.
     *
     *    length = Number of species in PhaseList
     *    indexing of PhaseList
     */
    std::vector<double> spMf_init_init_;

    //! Mole fraction of each species in each phase
    /*!
     *  DEPENDENT STATE VARIABLE for the electrode.
     *  This variable is kept synched with the spMoles_final variable at all times.
     *
     *    length = Number of species in PhaseList
     *    indexing of PhaseList
     */
    std::vector<double> spMf_final_;


    //! Mole fraction of each species in each phase at the start of the global time step
    /*!
     *  DEPENDENT STATE VARIABLE for the electrode.
     *  This variable is kept synched with the spMoles_init_init_ variable at all times.
     *
     *    length = Number of species in PhaseList
     *    indexing of PhaseList
     */
    std::vector<double> spMf_final_final_;


    //! Vector of species Electrochemical potentials
    /*!
     *    length = Number of species in PhaseList
     *    indexing of PhaseList
     */
// Deprecate? How often is this called that it is worth storing in memory rather than computing if needed?
    mutable std::vector<double> spElectroChemPot_;

    //! Vector of phase voltages
    /*!
     *  length = Number of total phases in PhaseList
     *  indexing of PhaseList
     */
    std::vector<double> phaseVoltages_;

    //! Vector of ReactingSurface objects in the problem
    /*!
     *  This vector has a length equal to the number of surfaces in the problem.
     *  If a surface doesn't have kinetics associated with it, the position is set to null.
     *
     *     Length = number of surfaces that may be present: numSurfaces_
     */
    std::vector<ReactingSurDomain*> RSD_List_;

    //! Vector of the number of reactions in a ReactingSurface object in the problem
    /*!
     *  This vector has length number of surfaces in the problem.
     *  If a surface doesn't have kinetics associated with it, the position is set to 0.
     *
     *     Length = number of surfaces that may be present: numSurfaces_
     */
    std::vector<int> numRxns_;

    //! Vector indicating that a surface is currently kinetically active
    /*!
     *   To be true, the surface must have a kinetics object and the surface area must have a
     *   nonzero positive value. Mole numbers for one side of the interfacial reaction
     *   should also be present.
     *
     *   This is used to trigger the calculation of the rates of progress of reactions that
     *   are on the surface.
     *
     *     Length = number of surfaces that may be present: numSurfaces_
     */
    std::vector<int> ActiveKineticsSurf_;

    //! Number of moles in each phase
    /*!
     *  This is kept synched with spMoles[] at all times
     */
    std::vector<double> phaseMoles_init_;

    //! Number of moles in each phase at the start of the integration time period
    /*!
     *  This is kept synched with spMoles_init_init[] at all times
     */
    std::vector<double> phaseMoles_init_init_;

    //! Number of moles in each phase
    /*!
     *  This is kept synched with spMoles_final[] at all times
     */
    std::vector<double> phaseMoles_final_;

    //! Number of moles in each phase at the end of the global step
    /*!
     *  This is kept synched with spMoles_final[] at all times
     */
    std::vector<double> phaseMoles_final_final_;

    //! Time derivative of the phase moles
    std::vector<double> phaseMoles_dot_;

    //! Vector of phase numbers for phases that are just coming into existence
    /*!
     *   This is a vector that is pertinent for intermediate time steps.
     *   A phase will be just Born when it is coming into existence at that
     *   particular intermediate time step.
     *
     *   length = number of phases being born.
     */
    std::vector<int> justBornPhase_;

    //! Vector of phase numbers for phases that are just dying
    /*!
     *   This is a vector that is pertinent for intermediate time steps.
     *   A phase will be just Died when it goes out of existence at that
     *   particular intermediate time step.
     *
     *  Note this vector may have an effect on the algorithm.
     *
     *  Length = number of total phases
     */
    std::vector<int> justDiedPhase_;

    // ---------------   SURFACE AREAS -----------------------------------------

    //! Number of external interfacial areas
    int numExternalInterfacialSurfaces_;

    //! True if the surface is an external surface, false otherwise
    /*!
     *  Dimensioned numSurfaces_
     */
    std::vector<bool> isExternalSurface_;

    //! Vector of the surface area for each Reacting Surface
    //! in the electrode.
    /*!
     *  Each surface is assumed to have a surface area. The surface area is a
     *  constitutive function  of the composition (i.e., State of Charge) of the electrode.
     *  In most constitutive models only one external reacting surface will be present at
     *  any one time.  This vector is over internal and external surfaces.
     *
     *  length = number of surfaces that may be present: numSurfaces_
     *  units m**2
     *
     *  This is the initial value at each time step
     */
    std::vector<double> surfaceAreaRS_init_;

    //! Vector of the surface area for each Reacting Surface in the electrode.
    /*!
     *  Each surface is assumed to have a surface area. The surface area is a
     *  constitutive function  of the composition (i.e., State of Charge) of the electrode.
     *  In most constitutive models only one external reacting surface will be present at
     *  any one time.  This vector is over internal and external surfaces.
     *
     *  length = number of surfaces that may be present: numSurfaces_
     *  units m**2
     *
     *  This is the final value at each time step
     */
    std::vector<double> surfaceAreaRS_final_;

    //! Vector of the surface area for each Reacting Surface in the electrode.
    /*!
     *  Each surface is assumed to have a surface area. The surface area is a
     *  constitutive function  of the composition (i.e., State of Charge) of the electrode.
     *  In most constitutive models only one external reacting surface will be present at
     *  any one time.  This vector is over internal and external surfaces.
     *
     *  length = number of surfaces that may be present: numSurfaces_
     *  units m**2
     *
     *  This is the initial value at the start of the global time step
     */
    std::vector<double> surfaceAreaRS_init_init_;

    //! Vector of the surface area for each Reacting Surface in the electrode.
    /*!
     *  Each surface is assumed to have a surface area. The surface area is a
     *  constitutive function  of the composition (i.e., State of Charge) of the electrode.
     *  In most constitutive models only one external reacting surface will be present at
     *  any one time.  This vector is over internal and external surfaces.
     *
     *  length = number of surfaces that may be present: numSurfaces_
     *  units m**2
     *
     *  This is the final global value at each time step
     */
    std::vector<double>  surfaceAreaRS_final_final_;

    //! Internal array containing the net production rate per area for each reacting surface
    /*!
     *   This is in the PhaseList indexing scheme
     *
     *  The column is the reacting surface index, while the row is the species PhaseList index
     *
     *  The dimensioning is m_NumTotSpecies by numSurfaces_.
     */
    Array2D spNetProdPerArea_List_;

    //! Internal vector containing the net production rate over the global integration period
    //! of all species in the electrode object.
    /*!
     *  This vector will have nonzero entries for the electrolyte phase and the electron
     *  phase even if the moles for these phases are not allowed to change within the object
     *  itself.
     *
     *   units kmols
     */
    std::vector<double> spMoleIntegratedSourceTerm_;

    //! Internal vector containing the net production rate over last intermediate step
    //! of all species in the electrode object.
    /*!
     *  This vector will have nonzero entries for the electrolyte phase and the electron
     *  phase even if the moles for these phases are not allowed to change within the object
     *  itself.
     *
     *   units kmols
     */
    std::vector<double> spMoleIntegratedSourceTermLast_;

    
    std::vector<double> enthalpyMolar_init_init_;
    std::vector<double> enthalpyMolar_init_;
    std::vector<double> enthalpyMolar_final_;
    std::vector<double> enthalpyMolar_final_final_;


    std::vector<double> entropyMolar_init_init_;
    std::vector<double> entropyMolar_init_;
    std::vector<double> entropyMolar_final_;
    std::vector<double> entropyMolar_final_final_;

    std::vector<double> chempotMolar_init_init_;
    std::vector<double> chempotMolar_init_;
    std::vector<double> chempotMolar_final_;
    std::vector<double> chempotMolar_final_final_;

    double integratedThermalEnergySourceTerm_;
    double integratedThermalEnergySourceTermLast_;

    double integratedThermalEnergySourceTerm_overpotential_;
    double integratedThermalEnergySourceTerm_overpotential_Last_;
    double integratedThermalEnergySourceTerm_reversibleEntropy_;
    double integratedThermalEnergySourceTerm_reversibleEntropy_Last_;

    //! Name of the electrode to be used in printouts
    std::string electrodeName_;

    //! Number of extra global reaction pathways specified
    int numExtraGlobalRxns;

    std::vector<EGRInput*> m_EGRList;

public:
    //! Pointer vector of ExtraGlobalRxn objects
    /*!
     *       Extra global reactions are described by this object
     */
    std::vector<Cantera::ExtraGlobalRxn*> m_egr;

    //! Pointer vector of RxmMolChange objects
    /*!
     *      Extra information for each reaction are described by this object
     */
    std::vector<RxnMolChange*> m_rmcEGR;

protected:
    //! Phase ID of the metal
    /*!
     *  this is the phase where the electron species exists.
     *      This object doesn't have to have a positive mole number to be active
     */
    int metalPhase_;

    //! Phase ID of the electrolyte solution
    /*!
     *  This is the phase where the product ions exist
     */
    int solnPhase_;

    //! Global index within this object of the electron species
    int kElectron_;

    //! Global index within each of the Reacting surface/ interfacial kinetics object
    //! for the electron species
    /*!
     *  Length is equal to the number of surfaces
     *
     *  A value of -1 indicates that either there is no surface kinetics object or there
     *  is no electron that participates in the surface reactions
     */
    std::vector<int> kKinSpecElectron_sph_;

    //! This is the Phi_metal - Phi_soln.
    //! In other words, the voltage drop across the interface.
    double deltaVoltage_;

    //! Amount of electrons that have left the electrode
    /*!
     *  Note this can be a negative quantity for cathode operations in batteries.
     *  Units for this internal storage are kmol.
     */
    double electronKmolDischargedToDate_;

    //! Discharge Capacity Coefficient
    /*!
     *  The capacity of an anode electrode is equal to the number of electrons that can
     *  be released. The capacity of a cathode is equal to the number of electrons that
     *  can be accumulated.
     *
     *  We determine the capacity left in an electrode by multiplying the number of moles
     *  of each species by this capacityLeftSpeciesCoeff_[] value to get the number
     *  of electrons that can be discharged.
     *
     *  We determine capacityLeftSpeciesCoeff_[] in the instantiation phase. Basically
     *  we hard-code the evaluation process based on species names.
     *
     *  The depth of discharge will be equal to the amount of electrons that
     *  have been discharged divided by the number of electrons that can be released.
     *
     *  For the capacity in multiple plateau situations, we assume that the
     *  capacity is equal to the number of electrons that can be extracted from all plateaus
     *  defined in the problem regardless of the E_0 voltage.
     */
    mutable std::vector<double> capacityLeftSpeciesCoeff_;


    //! Array of discharge Capacity coefficients
    mutable Cantera::Array2D capacityLeftSpeciesCoeffPlat_;

    //! Capacity Coefficient for Calculation of theoreteical zero DoD capacity
    /*!
     *  The capacity of an anode electrode is equal to the number of electrons that can
     *  be released. The capacity of a cathode is equal to the number of electrons that
     *  can be accumulated.
     *
     *  We determine the capacity of an electrode by multiplying the number of moles
     *  of each species by this capacityZeroDoDSpeciesCoeff_[] value to get the number
     *  of electrons that can be discharged.
     *
     *  We determine capacityZeroDoDSpeciesCoeff_[] in the instantiation phase. Basically
     *  we hard-code the evaluation process based on species names.
     *
     *  The depth of discharge will be equal to the amount of electrons that
     *  have been discharged divided by the number of electrons that can be released.
     *
     *  For the capacity in multiple plateau situations, we assume that the
     *  capacity is equal to the number of electrons that can be extracted from all plateaus
     *  defined in the problem regardless of the E_0 voltage.
     */
    mutable std::vector<double> capacityZeroDoDSpeciesCoeff_;

    //! Array of capacity coefficients
    mutable Cantera::Array2D capacityZeroDoDSpeciesCoeffPlat_;

    //! Initial value of the capacity of the Electrode
    /*!
     *  This is the initial capacity of the electrode before any degradation
     *
     *  units = kmol
     */
    double capacityInitialZeroDod_;

    //! Starting depth of discharge in coulombs
    /*!
     *   Initial value of the depth of discharge
     */
    double depthOfDischargeStarting_;

    //! Current going through the electrode
    /*!
     * Note the following convention:
     *      The current is positive for current going into the electrode
     *      and then into the electrolyte. Thus, under normal battery operation
     *      where the anode is negative and the cathode is positive, the
     *      current is positive going into the anode and negative going
     *      into the cathode.
     */
    double Icurrent_;

    //! Value of deltaG for each reaction in the current reacting domain (can also be used for species storage)
    /*!
     *  This is dimensioned to be the max of the max number of reactions in any reacting domain
     *  and the number of species.
     */
    mutable std::vector<double> deltaG_;

    //! Particle diameter that is input from the input file
    double inputParticleDiameter_;

    //! Number of Particles to follow
    /*!
     *   All the extrinsic properties of the object are multiplied by this value.
     *   This is the number of particles that are in the electrode object
     */
    double particleNumberToFollow_;

    //!  Radius of the exterior of the particle
    /*!
     *  This is at the start of the global time step
     */
    double Radius_exterior_init_init_;

    //!  Radius of the exterior of the particle
    /*!
     *  This is at the start of the subgrid time step
     */
    double Radius_exterior_init_;

    //!  Radius of the exterior of the particle
    /*!
     *  This is at the end of the subgrid time step
     */
    double Radius_exterior_final_;

    //!  Radius of the exterior of the particle
    /*!
     *  This is at the end of the global time step
     */
    double Radius_exterior_final_final_;

    //! Porosity of the electrode
    /*!
     *    This is the percentage of space that is occupied by the electrolyte. The electrode consists of
     *    spherical particles within this electrolyte
     */
    double porosity_;

    //! Molar value of atol
    double molarAtol_;

protected:
    XML_Node* xmlTimeIncrementData_;

    XML_Node* xmlTimeIncrementIntermediateData_;

    //! Storage of the  external state of the system in terms of an XML data structure
    XML_Node* xmlExternalData_init_init_;
    XML_Node* xmlExternalData_init_;
    XML_Node* xmlExternalData_final_;
    XML_Node* xmlExternalData_final_final_;

    //! Storage of the state of the system in terms of an XML data structure
    XML_Node* xmlStateData_init_init_;
    XML_Node* xmlStateData_init_;
    XML_Node* xmlStateData_final_;
    XML_Node* xmlStateData_final_final_;


    //!  Pointer to the object that is in charge of formulating the saved state and writing
    //!  that state out to an XML object
    EState* eState_final_;

    //! Base name of the solution file, currently defaults to "soln"
    std::string baseNameSoln_;

public:

    //! Electrode Chemistry model type
    /*!
     *   Electrode Chemistry model type is used to set the capacity coefficients for specific chemical mechanisms.
     *   This is used as an expedient to set up some problems. Names of phases are used to trigger the usage
     *   of this. The default is not to use this expedient.
     *
     *      0  Undetermined chemistry model type
     *      1  LiSi anode with an initial composition of pure Li13Si4(S) and a final state of pure Si(s).
     *      2  FeS2 anode with an initial composition of pure FeS2(S) and a final state of Li2S(s) and Fe(s).
     *      3  LiSi anode with one plateau and interstitial diffusion model
     *      4  FeS2 cathode with a simplified thermo
     */
    int electrodeChemistryModelType_;

    //! Domain number of the electrode object
    /*!
     *  This number refers to the domain number within 1DElectrode.
     */
    int electrodeDomainNumber_;

    //! Cell number within the domain
    int electrodeCellNumber_;

    //! Number of integrations (successful or failed) carried out by this object
    int counterNumberIntegrations_;

    //! Number of subintegrations (successful or failed) carried out by this object
    int counterNumberSubIntegrations_;

    //! Amount of printing to be carried out by the object
    /*!
     *   0 -> absolutely nothing is printed for a single call to integrate.
     *   1 -> One line summary per integrate call
     *   2 -> short description, points of interest: Table of nonlinear solve - one line per iteration
     *   3 -> Table is included -> More printing per integrate call
     *   4 -> Table is included -> More printing per integrate call, source terms included
     *   5 -> Algorithm information on the nonlinear iterates are printed out
     *   6 -> Additional info on the nonlinear iterates are printed out
     *   7 -> Additional info on the linear solve is printed out.
     *   8 -> Info on a per iterate of the linear solve is printed out.
     */
    int printLvl_;

    //! Specifies the amount of output that the Electrode object writes to its solution file
    /*!
     *    The level is given by the following table
     *            - 0     Zero output (default)
     *            - 1     Global time end
     *            - 2     global time increments
     *            - 3     intermediate steps
     *            - 4     intermediate Steps written out as global steps
     *            - 13    Jacobian steps written out as intermediate steps, with multiple steps
     *                    per global time steps
     */
    int printXMLLvl_;

    //! Print level of the CSV solution.
    /*!
     *  If this is turned on to greater than zero, than a CSV file is printed with the values at each time step
     *  printed on a single line.
     */
    int printCSVLvl_;

    //! Detailed Nonlinear Residual printouts
    /*!
     *  Select a particular level of detailed printouts from the Nonlinear Residual layer.
     *  This circumscribes the printLvl_ level above.
     */
    int detailedResidPrintFlag_;

    //! Boolean that turns on and off Nonlinear Residual layer printing
    bool enableExtraPrinting_;


private:
    // -----------------------------------------------------------------------------------------------------------------
    //! Class that calculates the residual for a phase Pop evaluation
    /*!
     *
     */
    class phasePop_Resid : public Cantera::ResidEval
    {
    public:
        phasePop_Resid(Electrode* ee, int iphaseTarget, double* const Xmf_stable,
                       double deltaTsubcycle);
        int evalSS(const doublereal t, const doublereal* const y,
                   doublereal* const r);

        int getInitialConditions(const doublereal t0, doublereal* const y,
                                 doublereal* const ydot);


        //! Return the number of equations in the equation system
        int nEquations() const;

        Electrode* ee_;
        int iphaseTarget_;
        double*   Xmf_stable_;
        double deltaTsubcycle_;

    public:

    };
    // -----------------------------------------------------------------------------------------------------------------



    friend class Cantera::EState;
    friend void Cantera::EState::copyElectrode_intoState(const Cantera::Electrode* const e);

    friend class Cantera::Electrode_Equilibrium;
};

}

int electrode_input(ELECTRODE_KEY_INPUT* input,  std::string commandFile,
                    BEInput::BlockEntry* cf);



#endif
/**********************************************************************************************************************/
