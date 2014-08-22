/*
 * $Id: InterfacialMassTransfer.h 507 2013-01-07 22:48:29Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _INTERFACIALMASSTRANSFER_H
#define _INTERFACIALMASSTRANSFER_H


#include "cantera/equilibrium.h"
#include "cantera/numerics/ResidEval.h"
#include "cantera/numerics/ResidJacEval.h"

#include "tok_input_util.h"

#include "PhaseList.h"
#include "ReactingSurDomain.h"
#include "ExternalField.h"

#include "InterfacialMassTransfer_input.h"

#include "mdp_allo.h"
#include <string>
#include <vector>
/*
 *-----------------------------------------------------------------------------
 *
 * Include file containing constant declarations for inputs to 
 * mpequil
 *
 *-----------------------------------------------------------------------------
 */


namespace Cantera {
  class IMT_KEY_INPUT;
}
class SubIntegrationHistory;

namespace Cantera {

  //!  Enum Type identifying the models for the Electrodes
  /*!
   *  
   */
  enum IMT_Types_Enum {
    UNKNOWN_IMT = -1,
    BASE_TYPE_IMT = 0,
    PSEUDOSS_NOSURF_IMT,
    KIN_1TO1DISTRIB_IMT,
    PSEUDOSS_NOSURF_DIFFBL_IMT
  };


  class TimeStepHistory
  {
  public:
    TimeStepHistory() :
      t_init_(0.0),
      t_final_(0.0),
      numNonLinSolves_(0),
      timeTypeSoln_(0),
      solnNorm_(0.0),
      wdotNorm_(0.0)
    {
    }

    TimeStepHistory(const TimeStepHistory &right) :
      t_init_(right.t_init_),
      t_final_(right.t_final_),
      numNonLinSolves_(right.numNonLinSolves_),
      timeTypeSoln_(right.timeTypeSoln_),
      solnNorm_(right.solnNorm_),
      wdotNorm_(right.wdotNorm_)
    {
    }

    double t_init_;
    double t_final_;
    int numNonLinSolves_;
    int timeTypeSoln_;
    double solnNorm_;
    double wdotNorm_;
  };

  class SubIntegrationHistory 
  {
  public:
    SubIntegrationHistory();

    SubIntegrationHistory(const SubIntegrationHistory &right);
 
   ~SubIntegrationHistory();

    SubIntegrationHistory &operator=(const SubIntegrationHistory &right);

    int nTimeSteps_;

    std::vector<TimeStepHistory> TimeStepList_;
    
    double GsolnNorm_;
    double GwdotNorm_;

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
   */
  class InterfacialMassTransfer : public Cantera::PhaseList {
  public:
    
    // ---------------------------------------------------------------------------------------------
    // ----------------------- BASIC SETUP ROUTINES  -----------------------------------------------
    // ---------------------------------------------------------------------------------------------

    //! Constructor
    InterfacialMassTransfer();

    //! Destructor
    virtual ~InterfacialMassTransfer();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    InterfacialMassTransfer(const InterfacialMassTransfer &right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    InterfacialMassTransfer & operator=(const InterfacialMassTransfer &right);

    //! Return the type of interfacial mass transport object 
    /*!
     *  Returns the enum type of the object. This is used in the factory routine.
     *
     *  @return Returns an enum type, called   IMT_Types_Enum
     */
    virtual IMT_Types_Enum imtType() const;

    //! Set the interface ID information
    void setID(int domainNum, int cellNum);

    //!  Setup the interface object
    /*!
     * @param ei    IMT_KEY_INPUT pointer object
     */
    virtual int model_create(IMT_KEY_INPUT *ei);

    //!  Set the initial conditions from the input file.
    /*!   
     *   (virtual from InterfacialMassTransfer)
     *   (This is a serial virtual function or an overload function)
     *
     *    This is one of the most important routines. It sets up the initial conditions of the interface
     *    from the input file. The interface object itself has been set up from a call to model_create().
     *    After the call to this routine, the interface should be internally ready to be integrated and reacted. 
     *    It takes its input from an IMT_KEY_INPUT object which specifies the setup of the interface
     *    object and the initial state of that object.
     *    
     *    The routine works like an onion initialization. The parent object is initialized before the 
     *    child. This means the child object first calls the parent, before it does its own initializations.
     * 
     * @param ei    IMT_KEY_INPUT pointer object
     *  
     *  @return  Returns zero if successful, and -1 if not successful.
     */
    virtual int setInitialConditions(IMT_KEY_INPUT *ei);

    // ---------------------------------------------------------------------------------------------
    // ----------------------- SPECIFY AND OBTAIN PROBLEM PARAMETERS -------------------------------
    // ---------------------------------------------------------------------------------------------

    // ------------------------------ OBTAIN STATIC PROBLEM INFORMATION ----------------------------

    //! Returns a pointer to the current outer reacting surface object
    /*!
     *
     */
    Cantera::ReactingSurDomain * currOuterReactingSurface();

    //! Returns a reacting surface object
    /*!
     *  @param iSurf  returns a pointer to that reacting surface
     */
    Cantera::ReactingSurDomain * reactingSurface(int iSurf);

  
    //! Set the phase existence flag in the electrode kinetics object so that kinetics 
    //! are calculated correctly
    /*!
     *    Flags are set in the kinetics object to tell the kinetics object which phases
     *    have zero moles.  The zero mole indicator is taken from phaseMoles_final_[]. Therefore,
     *    the final state is queried.
     *    There is a special case. Phases that are in justBornPhase_[] vector are allowed to
     *    be set to exist even if their phase moles are zero.
     *
     * @param doInactive   Boolean indicating whether to do inactive surface reactive surface domains 
     * @param assumeStableSingleSpeciesPhases Assume that single phases are stable. This
     *                         allows their production rates to be calculated
     */
    virtual void setPhaseExistenceForReactingSurfaces(bool doInactive, 
						      bool assumeStableSingleSpeciesPhases = false);

 

    //! Returns the index of a phase in the ReactionSurfaceDomain object
    //! given the index of that phase in the PhaseList object
    /*!
     * @param isk  Surface phase index, used to look up the ReactingSurfaceDomain object
     * @param PLph index of the phase in the PhaseList object, which is also the
     *             InterfacialMassTransfer_Model object.
     *
     *  @return  Returns the index of the phase in the current ReactingSurDomain
     *           object. A value of -1 in this slot means that the phase doesn't
     *           participate in the  current ReactingSurDomain object
     */
    int ReactingSurfacePhaseIndex(int isk, int PLph) const;

    //! Reactant stoichiometric coefficient
    /*!
     * Get the reactant stoichiometric coefficient for the kth global species
     * in the ith reaction of the reacting surface domain with index isk.
     */
    double reactantStoichCoeff(const int isk, int kGlobal, int i);

    //! Product stoichiometric coefficient
    /*!
     * Get the product stoichiometric coefficient for the kth global species
     * in the ith reaction of the reacting surface domain with index isk.
     */
    double productStoichCoeff(const int isk, int kGlobal, int i);

    //! Specify the external fields are discretized with respect to the time coordinate
    /*!
     *       0   Behavior within the global step is akin to backwards Euler. A step jump is 
     *           assumed to the global values at the end of the global time step even for
     *           intermediate times
     *       1   Behaviow within the global step is treated as a linear function between the 
     *           beginning values and the end values. 
     *
     *  @param externFieldTimeBehaviorType  Parameter describing the behavior.
     */
    void specifyExternalFieldTimeBehavior(EF_FieldTimeBehavior_Enum externFieldTimeBehaviorType);

    //! Report  how the external fields are discretized with respect to the time coordinate
    /*!
     *       0   Behavior within the global step is akin to backwards Euler. A step jump is 
     *           assumed to the global values at the end of the global time step even for
     *           intermediate times
     *       1   Behaviow within the global step is treated as a linear function between the 
     *           beginning values and the end values. 
     *
     *  @return  Returns a parameter describing the behavior.
     */
    EF_FieldTimeBehavior_Enum reportExternalFieldTimeBehavior() const;

    // ------------------------------ SPECIFY BASIC THERMO CONDITIONS  ------------------------------------

    //! Set the temperature and pressure of the electrode
    /*!
     * @param temperature    Temperature (Kelvin)
     * @param pressure       Pressure (Pa)
     */
    void setState_TP(double temperature, double presA, double presB = -1.0);

     
    //! Set the global final_final time
    /*!
     *  When we do this we are setting a pending state flag
     *
     *  @param t_final_final  Final time of the global step
     *  @param setFinal       set the final state as well as the final_final
     *                           Defaults to true.
     */
    void setStateFF_time(double t_final_final, bool setFinal = true);

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

   //! Set the time
    /*!
     *   Set the time
     */
    virtual void setTime(double time);

    // ------------------------------ SPECIFY PROBLEM PARAMETERS ------------------------------------

    //! Set the mole numbers and internal state of the solnA phase at the final time of the
    //! global step.
    /*!
     *  We set the mole numbers of the solnA phase separately from 
     *  the rest of the phases.
     *
     *  We always make sure that mole numbers are positive by clipping. We
     *  always make sure that mole fractions sum to one.
     *
     *  If we are not following mole numbers in the electrode, we set the
     *  total moles to the internal constant, electrolytePseudoMoles_, while
     *  using this vector to set the mole fractions.
     *
     * @param solnAMoleNum vector of mole numbers of the species in the
     *                     electrolyte phase. 
     *                     units = kmol
     *                     size = number of species in the electrolyte phase 
     *
     * @param tempA        Temperature of phase A
     *
     * @param presA        Pressure of phase A
     *
     * @param extFieldAValues External field values of A, the first one is usually the voltage.
     *
     * @param setInitial   Boolean indicating that we should set the initial values of the
     *                     solnA instead of the final values.(default false)
     */
    void setSolnA_BoundaryConcentrations(const double * const solnAMoleNum, double tempA, double presA,
					 const double * extFieldAValues = 0,
					 bool setInitial = false);

    //! Set the mole numbers and internal state of the solnB phase at the final time of the
    //! global step.
    /*!
     *  We set the mole numbers of the solnB phase separately from 
     *  the rest of the phases.
     *
     *  We always make sure that mole numbers are positive by clipping. We
     *  always make sure that mole fractions sum to one.
     *
     *  If we are not following mole numbers in the electrode, we set the
     *  total moles to the internal constant, electrolytePseudoMoles_, while
     *  using this vector to set the mole fractions.
     *
     * @param solnBMoleNum vector of mole numbers of the species in the
     *                     electrolyte phase. 
     *                     units = kmol
     *
     *                     size = number of species in the electrolyte phase 
     *
     * @param tempB        Temperature of phase B
     *
     * @param presB        Pressure of phase B
     *
     * @param extFieldBValues External field values of B, the first one is usually the voltage.
     *
     * @param setInitial   Boolean indicating that we should set the initial values of the
     *                     solnB instead of the final values (default false)
     */
    void setSolnB_BoundaryConcentrations(const double * const solnBMoleNum, double tempB, double presB,
					 const double * extFieldBValues = 0,
					 bool setInitial = false);

    //! Set the surface areas of surfaces within the model at the final time of the global
    //! step.
    /*!
     *   It is up to the object to determine what to do with this information, if anything.
     *
     *  @param surfaceAreasRS   Vector of surface areas (length equal to number of surfaces)
     *
     *  @param setInitial   Boolean indicating that we should set the initial values
     *                      instead of the final values (default false).
     */
    void setSurfaceAreas(const double * const surfaceAreasRS, bool setInitial = false);

    //! This sets the metal and solution voltage
    /*!
     *    This sets the metal and the electrolyte voltages
     *  This assumes that there are only two voltages in the system.
     */
    void setVoltages(const double phiMetal, const double phiElectrolyte);

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
    virtual void setPhaseMoleNumbers(int iph, const double * const moleNum);

    // ---------------------------------------------------------------------------------------------
    // ----------------------------- CARRY OUT INTEGRATIONS -----------------------------------------
    // ---------------------------------------------------------------------------------------------

    //! The internal state e must be kept for the initial and final times of an integration step.
    /*!
     *  This function advances the initial state to the final state that was calculated
     *  in the last integration step. 
     *
     * @param Tinitial   This is the new initial time. This time is compared against the "old"
     *                   final time, to see if there is any problem.
     */
    virtual void resetStartingCondition(double Tinitial);
    
   
    //!  Calculate the change in the state of the system when integrating from Tinitial to Tfinal
    /*!
     *  All information is kept internal within this routine. This may be done continuously
     *  and the solution is not updated. 
     *
     *  Note the tolerance parameters refere to the nonlinear solves within the calculation
     *  They do not refer to time step parameters.
     *
     *  @param deltaT              DeltaT for the integration step.
     *  @param GlobalRtolSrcTerm   Relative tolerance for nonlinear solves within the calculation
     *                             Defaults to 1.0E-3
     *  @param GlobalAtolSrcTerm   Absolute tolerance for nonlinear solves within the calculation
     *                             Defaults to 1.0E-12
     *
     *  @param fieldInterpolationType
     *
     *  @param subIntegrationType
     *  @param sih
     *
     *  @return Returns the number of subcycle steps it took to complete the full step.
     */
    virtual int integrate(double deltaT, double GlobalRtolSrcTerm = 1.0E-3, double GlobalAtolSrcTerm = 1.0E-12,
                          int fieldInterpolationType = 0, int subIntegrationType = 0, SubIntegrationHistory * sih = 0);

    //! Set the deltaT used for the subcycle step
    /*!
     *  The default setting for this is infinite. If it's infinite then deltaTSubcycle is set to the
     *  deltaT of the integration step. However, there are many times where a smaller natural delta T is
     *  required to solve the equation system
     */
    virtual void setDeltaTSubcycle(doublereal deltaTSubcycle);

    //! Set the maximum value of the subcycle delta T
    /*!
     *  @param deltaTSubcycle  Maximum value of the subcycle time step
     */
    void setDeltaTSubcycleMax(doublereal deltaTSubcycle);

    //!  Calculate the time derivatives of the mole numbers at the current subcycle step
    /*!
     *   This may be used as a virtual function
     */
    virtual void calculateTimeDerivatives(doublereal deltaTsubcycle);

    //! This function calculates the integrated source term given a final_final state value 
    /*!
     *  Thus, this is a consistency calculation.
     */
    void calcIntegratedSourceTerm();

    //!  Residual calculation for the solution of the Nonlinear integration problem
    /*!
     *  This is an intermediate integration step from tinit_ to tfinal_
     *
     * @param tfinal        Time                    (input) 
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
    virtual int integrateResid(const doublereal tfinal, const doublereal delta_t,
			       const doublereal * const y, const doublereal * const ydot,
			       doublereal * const resid,
			       const ResidEval_Type_Enum evalType, const int id_x, const doublereal delta_x);

 
    // ----------------------------------------------------------------------------------------------
    // ----------------------------- GET CONDITIONS OUT --------------------------------------------
    // ----------------------------------------------------------------------------------------------
    //
    //       (unless specified this is always at the final conditions and time
    //

    // ----------------------------- GET INSTANTANEOUS SOURCE TERMS --------------------------------

    //! Get the net production rates of all species in the object
    //! at the current final conditions from all surface kinetics objects
    /*!
     * @param isk   Surface index to get the net production rates from
     * @param net   Species net production rates [kmol/s]. Return the species
     */
    void getNetProductionRates(doublereal* const net) const;

    //! Get the net production rates of all species in the object
    //! at the current final conditions from one surface kinetics object
    /*!
     * @param isk   Surface index to get the net production rates from
     * @param net   Species net production rates [kmol/s]. Return the species
     */
    void getNetProductionRatesRSD(const int isk, doublereal* const net) const;

    //!  Returns the current and the net production rates of the phases in kg/m2/s from a single surface
    //!  at the current final conditions and t_final
    /*!
     *  Returns the net production rates of all phases from reactions on a single surface
     *  
     *  @param isk Surface ID to get the fluxes from.      
     *  @param phaseMassFlux  Returns the mass fluxes of the phases
     */
    void getPhaseMassFlux(doublereal* const phaseMassFlux) const;

    //!  Returns the current and the net production rates of the phases in kmol/m2/s from a single surface
    //!  at the current final conditions and t_final
    /*!
     *  Returns the net production rates of all phases from reactions on a single surface
     *  
     *  @param isk Surface ID to get the fluxes from.      
     *  @param phaseMassFlux  Returns the mass fluxes of the phases
     */
    void getPhaseMoleFlux(const int isk, doublereal* const phaseMoleFlux) const;

    //!  Returns the computed mass flux of species into the A phase at final conditions
    /*!
     *  This is also the creation rate of phase A at the final conditions.
     *
     *      @return Returns mass flux  of species into phase A (kg s-1)
     */
    double getPhaseAMassSourceTerm() const;

    //!  Returns the computed mass flux of species into the Bphase at final conditions
    /*!
     *  This is also the creation rate of phase B at the final conditions.
     *
     *      @return Returns mass flux  of species into phase B (kg s-1)
     */
    double getPhaseBMassSourceTerm() const;

    //! Returns the computed Stefan velocity for phase A in the object
    //! at the current final conditions. 
    /*!
     *  Multiply by the surface area to get the stefan volume production rate at the
     *  current final conditions.
     * 
     *    @return returns stefan velociy created in phase A (m s-1)
     */
    virtual double StefanVelocityPhaseA() const;

    //! Returns the integrated volume creation rate for phase B in the object
    //! over the current global time step
    /*!
     *  This can be turned into a velocity by dividing by the surface area within the
     *  object.
     *    @return returns volume created in phase B (m3 s-1)
     */
    virtual double StefanVelocityPhaseB() const;

    // ---------------------------- GET INTEGRATED SOURCE TERMS -------------------------------------

    //!  Returns the net production rates of all species in the object
    //!  over the last global integration step
    /*!
     *  We calculate a rate here by taking the total production amounds and then divide by the time step.
     *
     *   @param net   Species net production rates [kmol/s].
     */
    void getIntegratedProductionRates(doublereal* const net) const;

    //! Report the integrated source term for the electrode over an interval in time.
    /*!
     *  This is the net change in the moles of species defined in the electrode over that
     *  interval of time.
     *  
     *  @param spMoleDelta The end result in terms of the change in moles of species in the
     *                     electrode.
     *
     *  @return Tfinal    Final time to integrate to.
     */
    double integratedSourceTerm(doublereal* const spMoleDelta) const;

    //! Returns the integrated moles created for each phase in the object
    //! over the current global time step
    /*!
     *    @param  phaseMolesTransfered vector of moles transfered (length = number of total 
     *            phases in the object)
     *            units = kmol
     */
    virtual void getIntegratedPhaseMoleSourceTerm(doublereal* const phaseMolesCreated) const;

    //! Returns the integrated mass created for a particular phase in the object
    //! over the current global time step
    /*!
     *    @param  phaseMolesTransfered vector of moles transfered (length = number of total 
     *            phases in the object)
     *            units = kmol
     */
    virtual double getIntegratedPhaseMassSourceTerm(int iph) const;

    //! Returns the integrated mass created for phase A in the object
    //! over the current global time step
    /*!
     */
    virtual double getIntegratedPhaseAMassSourceTerm() const;
 
    //! Returns the integrated mass created for phase B in the object
    //! over the current global time step
    /*!
     *    @return returns mass created in phase B
     */
    virtual double getIntegratedPhaseBMassSourceTerm() const;

    //! Returns the integrated volume creation rate for phase A in the object
    //! over the current global time step
    /*!
     *  This can be turned into a velocity by dividing by the surface area within the
     *  object.
     *    @return returns volume created in phase A (m3 s-1)
     */
    virtual double integratedStefanVolumeRatePhaseA () const;

    //! Returns the integrated volume creation rate for phase B in the object
    //! over the current global time step
    /*!
     *  This can be turned into a velocity by dividing by the surface area within the
     *  object.
     *    @return returns volume created in phase B (m3 s-1)
     */
    virtual double integratedStefanVolumeRatePhaseB () const;

    // --------------------------- GET MOLE NUMBERS ------------------------------------------------

    //! Return the number of moles in a phase
    /*!
     * @param iph  Phase id

     * @return Returns the number of moles in kmol
     */
    double phaseMoles(int iph) const;

    //! Returns the number of moles of a n elemenet
    /*!
     *  @param  ie   the index of the element
     *
     * @return Returns the number of moles in kmol
     */
    double elementMoles(int ie) const;

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

    //! Get the total number of moles in the system
    /*!
     * @return returns the total number of moles in the system
     */
    double totalMoles() const;

    // --------------------------- GET THERMO AND VOLTAGE  ------------------------------------------------

    //! Get chemical potential of a species
    /*!
     *  @param globalSpeciesIndex   Index of a species
     *
     *  @return  returns the chemical potential (J/kmol)
     */
    double speciesChemPotential(int globalSpeciesIndex) const;

    //! Returns the current potential drop across the electrode
    /*!
     *   This is equal to phiMetal - phiSoln
     */
    double potentialDrop() const ;

    //! Return the voltage of a phase
    /*!
     * @param iph  Phase id
     *
     * @return Returns the voltage in volts
     */
    double phaseVoltage(int iph) const;

    // ------------------------- GET VOLUMES -----------------------------------------------------------

    //! Return the total volume of solid material in the object
    /*!
     *  We leave out the solnPhase_ volume
     *  units = m**3
     */
    virtual double SolidVol() const;

    //! Return the total volume in the InterfacialMassTransfer object
    /*!
     *   @return total volume (m**3)
     */
    double TotalVol() const;

    //! Return a vector of the phase volumes for all phases in the electrode
    /*!
     *  Note the vector is over surface phases as well. Currently, all surface phases have zero
     *  volume.
     *
     *  length = m_NumtotPhases
     *  units = m**3
     */
    void getPhaseVol(double * const phaseVols) const;

    // ---------------------- GET SURFACE AREAS -------------------------------------------------------

    //! Return a vector of surface phase areas.
    /*!
     *  Returns a vector of surface areas for surfaces defined in the electrode object.
     *
     * @param surfArea  vector of surface areas.
     *                  units = m2
     *                  length = number of surfaces defined in the PhaseList
     */
    void getSurfaceAreas(double * const surfArea) const;

    //! Sets the vector of surface phase areas.
    /*!
     * Overwrites the surface area values witin the object
     *
     * @param surfArea  vector of surface areas.
     *                  units = m2
     *                  length = number of surfaces defined in the PhaseList
     */
    void setSurfaceAreas(const double *const surfArea);

    //! Calculates the change in the surface area of all external and internal interfaces 
    /*!
     *  (virtual)
     *  variables to be potentially altered
     *   surfaceAreaRS_[];
     *   isExternalSurface[]
     *   numExternalInterfacialSurfaces_;
     */
    virtual double calcSurfaceAreaChange(double deltaT);

  
    //! Get the ID of the external surface, if there is one
    const std::vector<bool> & getExternalSurfaceBooleans() const;

    // --------------------------- INTERNAL UPDATE FUNCTIONS --------------------------------------

    //! Take the state (i.e., the final state) within the InterfacialMassTransfer_Model and push it down
    //! to the ThermoPhase Objects
    /*!
     *  We take the values of spMoles_final_[] and propagate them down to the ThermoPhase
     *  objects
     */
    virtual void updateState();

    //! Update the velocities within the object at the final conditions.
    /*!
     *  (virtual from InterfacialMassTransfer)
     * 
     *  In order to do this, we need the instantaenous source terms and the current value of the
     *  surface area over the local time step. Therefore, this routine must be called 
     *  after DspMoles_final_ has been calculated. 
     */
    virtual void updateVelocities();

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
    virtual void updatePhaseNumbers(int iph);

    // ---------------------------------------------------------------------------------------------
    // -------------------------------  SetState Functions -------------------------------------------------------
    // ---------------------------------------------------------------------------------------------

    //! Set the internal initial intermediate and initial global state from the internal final state
    /*!
     *  (non-virtual function)  -> function should onionize in-first.
     *
     *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
     *  routine as well.
     *
     * @param setInitInit   Boolean indicating whether you should set the init_init state as well
     */
    void setInitStateFromFinal_Oin(bool setInitInit = false);

    //! Set the internal initial intermediate and initial global state from the internal final state
    /*!
     *  (virtual function) 
     *
     *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
     *  routine as well.
     *
     * @param setInitInit   Boolean indicating whether you should set the init_init state as well
     */
    virtual void setInitStateFromFinal(bool setInitInit = false);

    //! Set the internal final intermediate and from the internal init state
    /*!
     *  (non-virtual function)  -> function should onionize in-first.
     *
     *  Set the final state from the init state. This is commonly called during a failed time step
     */
    void setFinalStateFromInit_Oin();


    //! Set the internal final intermediate state from the internal init state
    /*!
     *  (virtual function from Electrode) 
     *
     *  Set the final state from the init state. This is commonly called during a failed time step
     */
    virtual void setFinalStateFromInit();

    //! Set the internal initial intermediatefrom the internal initial global state
    /*!
     *  Set the intial state from the init init state. We also can set the final state from this
     *  routine as well.
     *
     *  The final_final is not touched.
     *
     * @param setFinal   Boolean indicating whether you should set the final as well
     */
    virtual void setInitStateFromInitInit(bool setFinal = false);

    //! Set the internal final global state from the internal final intermediate state
    /*!
     *  (non-virtual function from Electrode) 
     *
     *  Set the final_final state from the final state. This is commonly called at the end of successful base integration.
     *  This is an onionize in function
     */
    void setFinalFinalStateFromFinal_Oin();

    //! Set the internal final global state from the internal final intermediate state
    /*!
     *  (virtual function from Electrode) 
     *
     *  Set the final_final state from the final state. This is commonly called at the end of successful base integration
     */
    virtual void setFinalFinalStateFromFinal();

    // ---------------------------------------------------------------------------------------------
    // ------------------------------ PRINT ROUTINES -----------------------------------------------
    // ---------------------------------------------------------------------------------------------

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

    //! Print conditions of the electrode for the current integration step to stdout
    /*!
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values 
     */
    virtual void printInterfacialMassTransfer(int pSrc = 1, bool subTimeStep = false);

    //! Print a table of values for each phase in the problem  for the current integration step to stdout
    /*!
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values 
     */
    virtual void printInterfacialMassTransferPhaseList(int pSrc = 1, bool subTimeStep = false);

    //! Print condition of a phase in the electrode
    /*!
     *  @param iPhase        Print the phase
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values 
     */
    virtual void printInterfacialMassTransferPhase(int iPhase, int pSrc = 1,  bool subTimeStep = false);

 
    //! Write out CSV tabular data on the integrations
    /*!
     *  The idea is to print out tabular data about each intermediate run and about each
     *  global run
     *
     *  @param itype Type of the data
     *            - 0      Initialization information
     *            - 1      Normal intermediate information
     *            - 2      Normal intermediate information at the end of a global time step
     *            - -1     Failed intermediate calculation, a failure from a nonlinear solver step
     *            - -2     Failed calculation from the predictor step - not necessarily significant.
     */
    void writeCSVData(int itype);

    // ---------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------

  protected:

    //! true if there is a pending integration step
    /*!
     *  We keep track of whether there is a pending integration step with this variable
     */
    int pendingIntegratedStep_;

    //! Include A and B mole numbers in the solution vector
    /*!
     *     0 (default) mole number for A and B are not part of the solution vector
     *     1           Mole numbers for A and B are part of the solution vector.
     *                 Additional boundary conditions for the solution A and B sides
     *                 must be introduced.
     *  (one instance of A and B moles being part of the solution vector has been created)
     */
    int includeABMolesInSolnVector_;

    //! Integer representing the Problem type.
    /*!
     *  The identity of  what is held constant. Currently, 
     *   T and P are held constant, and this input is ignored 
     *   0 : we are tracking mole numbers of electrolyte
     *   1 - We are not tracking mole numbers of electrolyte
     *       Basically, mole number tracking stops at the electrolyte
     *       electrode interface.
     */
    int prob_type;

    //! Mapping between the phases in the MultiPhase object
    //! and the phases in this object/the PhaseList
    /*!
     *  Length = number of phases in the m_mp object
     *  iph = PhaseIndex_mp[impPhase]
     *     impPhase is the phase index in the mp object
     *     iph is the phase in this object (and the PhaseList object).
     */
    std::vector<int> PhaseIndex_mp_;

    //! InterfacialMassTransfer Type
    /*!
     */
    IMT_Types_Enum interfaceType_;

    //!  Specify the external fields are discretized with respect to the time coordinate
    /*!
     *       0   Behavior within the global step is akin to backwards Euler. A step jump is 
     *           assumed to the global values at the end of the global time step even for
     *           intermediate times.
     *       1   Behaviow within the global step is treated as a linear function between the 
     *           beginning values and the end values. 
     *
     *  The default is 0
     */
    EF_FieldTimeBehavior_Enum  externFieldTimeBehaviorType_;

    //! Initial value of the time, at the start of the current global time step
    double t_init_init_;

    //! Final value of the time, at the end of the current global time step
    double t_final_final_;

    //! Initial value of the time, at the start of the current intermediate time step
    double t_init_;

    //! Final value of the time, at the end of the current intermediate time step
    double t_final_;

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

 
    //! Flag for the choice for the initial subcycle time step for a new global time interval
    /*!
     *   0  Use the value of the time step chosen by the last subcycle of the previous global time step (default)
     *   1  Use the value of the time step chosen by the first subcycle of the previous global time step
     *   2  The subgrid integration is set to 1/10 of the global time step
     */
    int choiceDeltaTsubcycle_init_;

  protected:

    // ------------------------------ UNKNOWNS REPRESENTING THE INTERFACE PROBLEM --------------------------------------

    //! Temperature of Interface. Assuming it's isothermal here
    double Temp_;

    //! Pressure (Pascal) in Phase A at start of global step
    double Pres_A_Interface_init_init_;

    //! Pressure (Pascal) in Phase A at start of subgrid step
    double Pres_A_Interface_init_;

    //! Pressure (Pascal) in Phase A at end  of subgrid step
    double Pres_A_Interface_final_;

    //! Pressure (Pascal) in Phase A at start of global step
    double Pres_A_Interface_final_final_;


    //! Pressure (Pascal) in Phase A at start of global step
    double Pres_B_Interface_init_init_;

    //! Pressure (Pascal) in Phase A at start of subgrid step
    double Pres_B_Interface_init_;

    //! Pressure (Pascal) in Phase A at end  of subgrid step
    double Pres_B_Interface_final_;

    //! Pressure (Pascal) in Phase A at start of global step
    double Pres_B_Interface_final_final_;


    //! Interface mass-averaged Velocity at start of global step
    double Velo_S_Interface_init_init_;

    //! Interface  mass-averaged Velocity at start of subgrid step
    double Velo_S_Interface_init_;

    //! Interface mass-averaged  Velocity at end of subgrid step
    double Velo_S_Interface_final_;

    //! Interface mass-averaged Velocity at end of global step
    double Velo_S_Interface_final_final_;

    //! Normal mass-averaged Velocity of phase A at start of global step
    double Velo_A_Interface_init_init_;

    //! Normal mass-averaged Velocity of phase A  at start of subgrid step
    double Velo_A_Interface_init_;

    //! Normal mass-averaged Velocity of phase A  at end of subgrid step
    double Velo_A_Interface_final_;

    //! Normal mass-averaged Velocity of phase A at end of global step
    double Velo_A_Interface_final_final_;

    //! Normal mass-averaged Velocity of phase B at start of global step
    double Velo_B_Interface_init_init_;

    //! Normal mass-averaged Velocity of phase B  at start of subgrid step
    double Velo_B_Interface_init_;

    //! Normal mass-averaged Velocity of phase B  at end of subgrid step
    double Velo_B_Interface_final_;

    //! Normal mass-averaged Velocity of phase B at end of global step
    double Velo_B_Interface_final_final_;

    //! Reference frame for the velocity for the calculation
    /*!
     *   Everything is migrating with a certain velocity. This is subtracted from all normal velocity expressions.
     */
    double Velo_ReferenceFrame_;

    double BLThickness_A_init_;
    double BLThickness_A_final_;
    double BLThickness_A_init_init_;
    double BLThickness_A_final_final_;

    double BLThickness_B_init_;
    double BLThickness_B_final_;
    double BLThickness_B_init_init_;
    double BLThickness_B_final_final_;

    //! Number of moles of each species in each phase at the start of the integration step
    /*!
     *  This will only be updated when we are moving onto the next outside time step.
     *  this is true of all entities entitled "_init_init_"
     *
     *    length = Number of species in PhaseList
     *    indexing of PhaseList
     */
    std::vector<double> spMoles_init_init_;

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

    //! Number of moles of each species in each phase at the end of each subcycle of the integration step
    /*!
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

    std::vector<double> spMf_init_init_;

    //! Mole fraction of each species in each phase at the start of the time step
    /*!
     *  DEPENDENT STATE VARIABLE for the electrode.
     *  This variable is kept synched with the spMoles_init_ variable at all times.
     *
     *    length = Number of species in PhaseList
     *    indexing of PhaseList
     */
    std::vector<double> spMf_init_;
  

    //! Mole fraction of each species in each phase
    /*!
     *  DEPENDENT STATE VARIABLE for the electrode.
     *  This variable is kept synched with the spMoles_final variable at all times.
     *
     *    length = Number of species in PhaseList
     *    indexing of PhaseList
     */
    std::vector<double> spMf_final_;

    //! Time derivative of the species mole numbers
    /*!
     *  This is the time derivative of the subcycle time step.
     *  It's calculated from spMoles_final_[] and spMoels_init_[].
     *  It may not be kept current.
     */
    std::vector<double> spMoles_dot_;

    //! Predicted 
    std::vector<double> spMoles_predict_;
 
    //! Vector of molar volumes of each of the phases in the mechanism
    /*!
     *  This has units of m**3 / kmol. 
     *  It has length equal to the total number of phases (vol and surf) in the mechanism
     */
    std::vector<double> phaseMolarVolumes_;


    double Pres_solnA_BC_init_init_;
    double Pres_solnA_BC_init_;
    double Pres_solnA_BC_final_;
    double Pres_solnA_BC_final_final_;



    double Pres_solnB_BC_init_init_;
    double Pres_solnB_BC_init_;
    double Pres_solnB_BC_final_;
    double Pres_solnB_BC_final_final_;


    //!  Mole fraction of solnA at the boundary at the start of global time step
    std::vector<doublereal> spMF_solnA_BC_init_init_;
       //!  Mole fraction of solnA at the boundary at the end of global time step
    std::vector<doublereal> spMF_solnA_BC_init_;

   //!  Mole fraction of solnA at the boundary at the end of global time step
    std::vector<doublereal> spMF_solnA_BC_final_;

    //!  Mole fraction of solnA at the boundary at the end of global time step
    std::vector<doublereal> spMF_solnA_BC_final_final_;

    //!  Mole fraction of solnB at the boundary at the start of global time step
    std::vector<doublereal> spMF_solnB_BC_init_init_;
    
    //!  Mole fraction of solnB at the boundary at the end of subgrid iteration step
    std::vector<doublereal> spMF_solnB_BC_init_;

   //!  Mole fraction of solnB at the boundary at the end of subgrid iteration step
    std::vector<doublereal> spMF_solnB_BC_final_;

    //!  Mole fraction of solnB at the boundary at the end of global time step
    std::vector<doublereal> spMF_solnB_BC_final_final_;

 

    //! Partial molar volumes of all of the species species
    /*!
     *  Length = global number of species in PhaseList species vector
     *  units = m**3 kmol-1
     */
    std::vector<double> VolPM_; 
 

    //! Vector of species chemical potentials
    /*!
     *    length = Number of species in PhaseList
     *    indexing of PhaseList
     */
    mutable std::vector<double> spChemPot_;

 

    //! Vector of ReactingSurface objects in the problem
    /*!
     *  This vector has length number of surfaces in the problem
     *  If a surface doesn't have kinetics associated with it, the position is set to null.
     */
    std::vector<ReactingSurDomain *> RSD_List_;

    //! Vector of the number of reactions in a ReactingSurface object in the problem
    /*!
     *  This vector has length number of surfaces in the problem
     *  If a surface doesn't have kinetics associated with it, the position is set to null.
     */
    std::vector<int> numRxns_;

    //! Vector indicating that a surface is currently kinetically active
    /*!
     *   To be true, the surface must have a kinetics object and the surface area must have a 
     *   nonzero positive value. Mole numbers for one side of the interfacial reaction
     *   should also be present.
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

    //! Time derivative of the phase moles
    std::vector<double> phaseMoles_dot_;
   

    // ---------------   SURFACE AREAS -----------------------------------------
  
    //! Number of surfaces in the problem
    /*!
     *   The number of surfaces is initially identified with the number of surface phases
     *   However, this can change. Each surface will have a surface area associated with it
     *   and it may have an active Kinetics object associated with it.
     *   
     *   Surfaces which don't have a surface phase associated with it are added onto the end 
     *   of the domain. They also aren't associated with a surface species of any kind. And,
     *   currently, they don't exist.
     */
    int numSurfaces_;

    //! Vector of the surface area for each Reacting Surface
    //! in the electrode. 
    /*!
     *  Each surface is assumed to have a surface area. The surface area is a 
     *  constitutive function  of the composition (i.e., State of Charge) of the electrode.
     *  In most constitutive models only one external reacting surface will be present at
     *  any one time.  This vector is over internal and external surfaces.
     *
     *  length = number of surfaces that may be present
     *  units m**2
     *
     *  This is the initial value at each time step
     */
    std::vector<double> surfaceAreaRS_init_init_;

    //! Vector of the surface area for each Reacting Surface
    //! in the electrode. 
    /*!
     *  Each surface is assumed to have a surface area. The surface area is a 
     *  constitutive function  of the composition (i.e., State of Charge) of the electrode.
     *  In most constitutive models only one external reacting surface will be present at
     *  any one time.  This vector is over internal and external surfaces.
     *
     *  length = number of external surfaces that may be present (numSurfaces_)
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
     *  length = number of external surfaces that may be present
     *  units m**2
     *
     *  This is the final value at each time step
     */
    std::vector<double> surfaceAreaRS_final_;

    //! Internal array containing the net production rate per area for each reacting surface
    /*!
     *   This is in the PhaseList indexing scheme
     *
     *  The column is the reacting surface index, while the row is the species PhaseList index
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

    //! Internal vector containing the net production rate of all species in the interface object.
    /*!
     *  This vector will have nonzero entries for the solnA and solnB phases
     *  even if the moles for these phases are not allowed to change within the object
     *  itself.
     * 
     *   units kmols / sec
     */
    std::vector<double> DspMoles_final_;

  public:
    //! Title of the electrode
    std::string Title_;

  
    //! Phase ID of the solution on the rhs
    /*!
     * 
     */
    int solnAPhase_;

    //! Number of species in the A phase
    int nSpeciesA_;

    //! Phase ID of the solution on the lhs
    /*!
     * 
     */
    int solnBPhase_; 

    //! Number of species in the B phase
    int nSpeciesB_;


    //! global index within this object for the electron species
    int kElectron_; 

  protected:
    //! Global index within each of the Reacting surface/ interfacial kinetics object
    //! for the electron species
    /*!
     *  Length is equal to the number of surface phases
     *
     *  A value of -1 indicates that either there is no surface kinetics object or there
     *  is no electron that participates in the surface reactions
     */
    std::vector<int> kKinSpecElectron_sph_;

    // Estimate of the open circuit voltage at the present time. This is the deltaVoltage_
    // defined just above for which the net current transfer is zero. 
    /*!
     * Obviously this depends upon the concentrations.
     */
    double openCircuitVoltageEst_;


    //! value of deltaG
    mutable std::vector<double> deltaG_;


    //! Molar value of atol
    double molarAtol_;

  public:


    //! Domain number of the electrode object
    /*!
     *  This number refers to the domain number within 1DInterfacialMassTransfer.
     */
    int DomainNumber_;

    //! Cell number within the domain
    int CellNumber_;

    //! Number of integrations (successful or failed) carried out by this object
    int counterNumberIntegrations_;
   
    //! Number of subintegrations (successful or failed) carried out by this object
    int counterNumberSubIntegrations_;
   
    //! Number of subintegrations carried out by this object on the last subintegration 
    int counterNumberLastSubIntegrations_;

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

    //! Detailed Nonlinear Residual printouts
    /*!
     *  Select a particular level of detailed printouts from the Nonlinear Residual layer.
     *  This circumscribes the printLvl_ level above.
     */
    int detailedResidPrintFlag_;

    //! Boolean that turns on and off Nonlinear Residual layer printing
    bool enableExtraPrinting_;


  };

  int imt_model_input(Cantera::IMT_KEY_INPUT *input,  std::string commandFile, BEInput::BlockEntry *cf);

}



#endif 
/*****************************************************************************/
