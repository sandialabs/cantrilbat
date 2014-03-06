/*
 * $Id: Electrode_MP_RxnExtent.h 604 2013-05-24 16:27:35Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_MP_RXNEXTENT_H
#define _ELECTRODE_MP_RXNEXTENT_H




#include "tok_input_util.h"

#include "Electrode.h"
#include "Electrode_Integrator.h"

#include "cantera/numerics/NonlinearSolver.h"

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


class ELECTRODE_KEY_INPUT;
class EGRInput;
namespace Cantera
{

class ELECTRODE_MP_RxnExtent_KEY_INPUT : public ELECTRODE_KEY_INPUT
{
public:

    //! Constructor
    ELECTRODE_MP_RxnExtent_KEY_INPUT(int printLvl = 0);

    //! Destructor
    virtual ~ELECTRODE_MP_RxnExtent_KEY_INPUT();

    ELECTRODE_MP_RxnExtent_KEY_INPUT(const ELECTRODE_MP_RxnExtent_KEY_INPUT& right);

    ELECTRODE_MP_RxnExtent_KEY_INPUT& operator=(const ELECTRODE_MP_RxnExtent_KEY_INPUT& right);

    //!  First pass through the child setup system
    /*!
     *    Typically we will fill in all vectors that depend on the value of numRegions_ in this
     *    pass.
     *  @param cf    Pointer to the BlockEntry record
     */
    void setup_input_child1(BEInput::BlockEntry* cf);

    //!  Second pass through the child setup system
    /*!
     *    Typically we will fill in all vectors that depend on the value of numRegions_ in this
     *    pass.
     *  @param cf    Pointer to the BlockEntry record
     */
    void setup_input_child2(BEInput::BlockEntry* cf);

    //! Number of plateau regions in the model
    /*!
     *  The number of regions refers to the continuity of the electrode open circuit model.
     *  Within each region the ss open circuit model is continuous. Either the voltage is
     *  constant or it is a linear function of the extent of reaction variable.
     */
    int numRegions_;

    //! Solid state diffusion model identification
    /*!
     *   0 = There is no diffusion model (default)
     *   1 = There is solid state diffusion limitations. The model assumes pseudo-steady state
     *       The limiting reacting surface location is determined by the  locationOfReactingSurface_ model.
     *       The details of the model is given in the writeup.
     */
    int solidDiffusionModel_;

    //! Location of the reacting surface
    /*!
     *   0 = Surface area is fixed at the initial conditions. (default)
     *   1 = Surface is located internal to the particle. Surface area of the reacting surface is calculated as moving
     *       front dependent on the extent of reaction during the current plateau
     *   2 = Surface is located at the exterior of the particle, which may be growing.
     */
    int locationOfReactingSurface_;


    //! Diffusion coefficient in the outer region of a two region spherical model
    /*!
     *  We identify the solids by Zones, which are related to regions. Generally, when the extent of reaction is
     *  in the zeroeth region, we are using diffusionCoeffRegions_[1], i.e., the first zone, because that is
     *  zone 1 in the model. Zone 0 would be diffusion in the inner core, whose value we never
     *  actually need. It's actually conceptually simpler this way.
     *
     *  Length = numRegions_ + 1
     *  Units = m**2 s-1
     */

    std::vector<double> diffusionCoeffRegions_;


    //! Molar volume of each of the region ends
    /*!
     *     Length = numRegions_ + 1
     *     units = m**3 kmol-1
     */
    std::vector<doublereal> molarVolumeRegionBoundaries_;

    std::vector<doublereal> rxnPerturbRegions_ ;

};

//!  Newman Model
/*!
 *
 *
 *   Treatment of surfaces within the model
 *
 *   Formally there are now two surfaces within the model. The first surface is the
 *   assigned to the inner radius. The second surface is the external surface.
 *
 *   The Reacting surface domain is assigned to either of the surfaces according to whether the user
 *   wants to treat the rate limiting step as the inner reaction or the outer reaction.
 *   This choice is currently hardwared to the outer via setting locationOfReactingSurface_ = 1.
 *
 *   The inner radius is calculated from a calculation involving the relative distance through
 *   the current region. This relative value is equated with the volume of the total plateau
 *   that is current calculated, and then translated into a relative radius of the external
 *   radius value.
 *
 */
class Electrode_MP_RxnExtent : public Cantera::Electrode_Integrator
{
public:

    //! Constructor
    Electrode_MP_RxnExtent();

    //! Destructor
    virtual ~Electrode_MP_RxnExtent();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    Electrode_MP_RxnExtent(const Electrode_MP_RxnExtent& right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    Electrode_MP_RxnExtent& operator=(const Electrode_MP_RxnExtent& right);

    //! Return the type of electrode
    /*!
     *  Returns the enum type of the electrode. This is used in the factory routine.
     *
     *  @return Returns an enum type, called   Electrode_Types_Enum
     */
    virtual Electrode_Types_Enum electrodeType() const;

    //! Set the electrode ID information
    void setID(int domainNum, int cellNum);

    //! Change the ELECTRODE_KEY_INPUT to the correct type, parsing the command file again, and then returning the
    //!  the struct in its place
    virtual int electrode_input_child(ELECTRODE_KEY_INPUT** ei);

    //!  Setup the electrode
    /*!
     * @param ei    ELECTRODE_KEY_INPUT pointer object
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



    virtual int electrode_stateSave_create();

    //! Set the sizes of the electrode from the input parameters
    /*!
     *  We resize all of the information within the electrode from the input parameters
     *
     * @param electrodeArea   Area of the electrode
     * @param electrodeThickness  Width of the electrode
     * @param porosity        Volume of the electrolyte phase
     */
    virtual void setElectrodeSizeParams(doublereal electrodeArea, doublereal electrodeThickness, doublereal porosity);

    //! Resize the solid phase and electrolyte mole numbers within the object
    /*!
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

    //! Find the initial conditions for the OCV
    void developBaseE0();

    //!  Change the Heat of formation of compound B in the reaction in order to generate a given open circuit
    //!  voltage at the current temperature
    /*
     *  @param  Value of the standard state OCV that is desired
     *  @param  doPrint  if nonzero, then detailed printing is done from this routine
     */
    void changeToE0(double E0, int doPrint = 0);

    //! Calculate the mole numbers of the solid phases given the relative extent of reaction
    /*!
     *
     */
    void relativeExtentRxnToMoles_final();

    //! Calculate the mole numbers of the solid phases given the relative extent of reaction
    /*!
     *  Calculate the relative extent of reaction given the mole numbers of the solid phases.
     *  This is a reversal of the previous operation. The two must be inverses of each other.
     */
    void molesToRelativeExtentRxn_final();

    //! Set the temperature and pressure of the electrode
    /*!
     * @param temperature    Temperature (Kelvin)
     * @param pressure       Pressure (Pa)
     */
    void setState_TP(double temperature, double pressure);

    // -------------------------  VOLUMES -----------------------------------------------------------

    //!    Return the total volume of solid material
    /*!
     *  (virtual from Electrode.h)
     *
     *   This returns the value in the _final_ state
     *
     *       This is the main routine for calculating the volume of electrode material.
     *       Here we redo the model. The total moles of A + B are used. However,
     *       the molar volume is obtained from the relative extent of reaction within the current plateau
     *       and the mole volumes of the two end points of the current plateau.
     *
     *       We leave out the solnPhase_ volume from the calculation
     *       units = m**3
     */
    virtual double SolidVol() const;

    //! Return the total volume in the electrode
    /*!
     *  (virtual from Electrode.h)
     *
     *   This returns the value in the _final_ state
     *   We have to redo these because the molar volume approximation is redone in these routines.
     *
     *   @return total volume (m**3)
     */
    virtual double TotalVol(bool ignoreErrors = false) const;

    //!  Returns the total moles in the electrode phases of the electrode
    /*!
     *  @return total moles (kmol)
     */
    virtual  double SolidTotalMoles() const;

    //! Return a vector of the phase volumes for all phases in the electrode
    /*!
     *  Note the vector is over surface phases as well. Currently, all surface phases have zero
     *  volume.
     *
     *  length = m_NumTotPhases
     *  units = m**3
     */
    virtual void getPhaseVol(double* const phaseVols) const;

    // ---------------------- SURFACE AREAS -------------------------------------------------------


    //! Take the state (i.e., the final state) within the Electrode_Model and push it down
    //! to the ThermoPhase Objects
    /*!
     *  We take the values of spMoles_final_[] and propagate them down to the ThermoPhase
     *  objects in the electrode.
     *  Here we also take the value of RelativeExtentRxn_final_ and modify the open circuit potential
     *  and other parameters based on its value.
     */
    virtual void updateState();

    //!  Recalculate the surface areas of the surfaces for the final state
    /*!
     *    (virtual function from Electrode)
     *
     *    We used the internal variable locationOfReactingSurface_ to determine the behavior.
     *    A value of zero indicates that the surface 0 follows the reaction front as it goes from outer to inner as
     *    a function of the % though the plateau.
     *    A value of locationOfReactingSurface_ = 1 indicates that the surface 0 follows the exterior surface of the particle
     *
     *    We also assume that the surface area is equal to the particle surface area multiplied by the numbers of particles.
     *
     *
     *    Dependent StateVariables Used
     *         Radius_exterior_final;
     *         particleNumberToFollow_
     *
     *    Dependent StateVariables Calculated
     *          surfaceAreaRS_final_[]
     */
    virtual void updateSurfaceAreas();

    //!  Extract the ROP of the two reaction fronts from Cantera within this routine
    /*!
     *  In this routine we calculate the rates of progress from the two surfaces
     *  The vectors are filled in:
     *
     *        ROP_outer_[jRxn]
     *        ROP_inner_[jRxn]
     *        spNetProdPerArea_List_[isk][kIndexKin]
     *        justBornPhase_[jph]
     */
    void extractInfo();

    double modifyROPForDiffusion();

    //! Extract various quantities from the lumped reaction for use in calculating the diffusion approximation
    void calculatekABForDa();

    //!   Calcualte the Damkoeler number assuming the reaction occurs on the outer surface, and we have interstial diffusion
    //!   of a neutral diffusing through a region to a new reaction front, whose reactions are fast so that they are
    //!   in equilibrium
    /*!
     *  We also can do checks here to make sure that the approximation makes sense. In particular we should
     *  make sure that the mole fraction is between 0 and 1 for the interstitial diffuser.
     */
    void calculateDaOuter();

    void calculateDaInner();

    //! Update the molar production rates for all annular regions
    /*!
     *  Requires that SolidInnerKSpecies_IS_, SolidInnerKSpeciesStoichCoeff_IS_, and the like
     *  and ROP_inner_ and ROP_outer_ are already known
     */
    void updateMoleRatesFinal();

    //! Calculate the production rate of species in the electrode at the final time of the time step
    /*!
     *   We take the ROP_inner_[] and ROP_outer_[] rates of progress, combine them with the surface area calculation,
     *   and the stoichiometric coefficients to calculate the DspMoles_final_[], which is production rate for
     *   all species in the electrode.
     */
    void updateSpeciesMoleChangeFinal();

    //! returns the old region

    int changeRegion(int newRegion);


    // ---------------------------------------------------------------------------------------------
    // ----------------------------- CARRY OUT INTEGRATIONS -----------------------------------------
    // ---------------------------------------------------------------------------------------------

    //! The internal state of the electrode must be kept for the initial and
    //! final times of an integration step.
    /*!
     *  This function advances the initial state to the final state that was calculated
     *  in the last integration step. If the initial time is input, then the code doesn't advance
     *  or change anything.
     *
     * @param Tinitial   This is the New initial time. This time is compared against the "old"
     *                   final time, to see if there is any problem.
     */
    void  resetStartingCondition(double Tinitial, bool doTestsAlways = false);



    //! Set the base tolerances for the nonlinear solver within the integrator
    /*!
     *   The tolerances are based on controlling the integrated electron source term
     *   for the electrode over the integration interval.  The integrated source term
     *   has units of kmol.
     *
     *   Because the electron is only one molar quantity within a bunch of molar quantities,
     *   this requirement will entail that we control the source terms of all species within the
     *   electrode to the tolerance requirements of the electron source term.
     *
     *   @param rtolResid  Relative tolerance allowed for the electron source term over the interval.
     *                     This is a unitless quantity
     *   @param atolResid  Absolute tolerance cutoff for the electron source term over the interval.
     *                     Below this value we do not care about the results.
     *                     atol has units of kmol.
     */
    virtual void setNLSGlobalSrcTermTolerances(double rtolResid);

    //! Set the Residual absolute error tolerances
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *   Set the absolute error tolerances fror the nonlinear solvers. This is called at the top
     *   of the integrator() routine.
     *
     *   Calculates atolNLS_[]
     *   Calculates atolResidNLS_[]
     */
    virtual void setResidAtolNLS();


    //! Predict the solution
    /*!
     * Ok at this point we have a time step deltalimiTsubcycle_
     * and initial conditions consisting of phaseMoles_init_ and spMF_init_.
     * We now calculate predicted solution components from these conditions.
     *
     * @return   Returns the success of the operation
     *                 1  A predicted solution is achieved
     *                 2  A predicted solution with a multispecies phase pop is acheived
     *                 0  A predicted solution is not achieved, but go ahead anyway
     *                -1  The predictor suggests that the time step be reduced and a retry occur.
     */
    virtual int predictSoln();

    //! Predict the solution using a special limiting behavior
    /*!
     * @return   Returns the success of the operation
     *                 1  A predicted solution is achieved
     *                 2  A predicted solution with a multispecies phase pop is acheived
     *                 0  A predicted solution is not achieved, but go ahead anyway
     *                -1  The predictor suggests that the time step be reduced and a retry occur.
     */
    int  predictSoln_SpeciaType1();

    // ---------------------------------------------------------------------------------------------
    // ---------------------------- SOLUTION OF NONLINEAR TIME DEPENDENT SYSTEM  --------------------
    // ---------------------------------------------------------------------------------------------


    //! Return the number of equations in the equation system
    virtual  int nEquations() const;

    //! Unpack the soln vector
    /*!
     *  This function unpacks the solution vector into deltaTsubcycleCalc_  and   RelativeExtentRxn_final_
     */
    void unpackNonlinSolnVector(const double* const y);


    // Main internal routine to calculate the rate constant
    /*
     *  This routine calculates the functional at the current stepsize, deltaTsubcycle_.
     *  A new stepsize, deltaTsubcycleCalc_, is calculated within this routine for changes
     *  in topology.
     *
     *  This routine calcules yval_retn, which is the calculated value of the residual for the
     *  nonlinear function
     *
     *   resid[i] = y[i] - yval_retn[i]
     *
     *   resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;
     *   resid[1] = RxnExtent_final_ - ( RxnExtent_init_ + Icurr * sa * (t_final_ - t_init_) * DiffResistance)
     *
     *  The formulation of the solution vector is as follows. The solution vector will consist of the following form
     *
     *     y =   deltaTsubcycleCalc_
     *          RxnExtent_final_
     *
     *  @param resid    value of the residual
     *
     *
     *  @return  1 Means a good calculation that produces a valid result
     *           0 Bad calculation that means that the current nonlinear iteration should be terminated
     */
    int calcResid(doublereal* const resid, const ResidEval_Type_Enum evalType);

    //!  Gather the predicted solution values and the predicted integrated source terms
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *  Both the predicted solution values and the predicted integrated source terms are used
     *  in the time step control
     */
    virtual void gatherIntegratedSrcPrediction();

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



    //!  Calculate the integrated source terms and do other items now that we have a completed time step
    /*!
     *  Calculate source terms on completion of a step. At this point we have solved the nonlinear problem
     *  for the current step, and we are calculating post-processed quantities like source terms.
     */
    virtual void calcSrcTermsOnCompletedStep();


    //! Accumulate src terms and other results from the local step into the global holding bins.
    /*!
     *  Accumulate source terms on completion of a step. At this point we have solved the nonlinear problem
     *  for the current step and we have satisfied all accuracy requirements.
     *  The step is good. We now accumulate the results before going on to a new local step.
     */
    virtual void accumulateSrcTermsOnCompletedStep(bool remove = false);


    double l0normM(const std::vector<double>& v1, const std::vector<double>& v2, int num,
                   const std::vector<double>& atolVec, const double rtol) const;

    //!
    void checkRegion(int regionID) const;

    //! Pack the nonlinear solver proplem
    /*!
     *  formulate the nonlinear solver problem to be solved.
     *     Fields to be filled in
     *             yvalNLS_
     *             ylowNLS_
     *             yhighNLS_
     *             deltaBoundsMagnitudesNLS_
     */
    virtual void initialPackSolver_nonlinFunction();


    //!  Calculate the norm of the difference between the predicted answer and the final converged answer
    //!  for the current time step
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *   The norm calculated by this routine is used to determine whether the time step is accurate enough.
     *
     *  @return    Returns the norm of the difference. Normally this is the L2 norm of the difference
     */
    virtual double predictorCorrectorWeightedSolnNorm(const std::vector<double>& yvalNLS_);

    //! Calculate the vector of predicted errors in the source terms that this integrator is responsible for
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *    In the base implementation we assume that the there are just one source term, the electron
     *    source term.
     *    However, this will be wrong in almost all cases.
     *    The number of source terms is unrelated to the number of unknowns in the nonlinear problem.
     *    Source terms will have units associated with them.
     *    For example the integrated source term for electrons will have units of kmol
     */
    virtual void predictorCorrectorGlobalSrcTermErrorVector();

    //!  Calculate the norm of the errors in the global source terms
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *   This routine make use of the source term error vector along with rtols and atols for the
     *   individual source terms to calculated a normalized error measure. This is the single number
     *   that the integration routine will try to control as it calculates a time stepping strategy.
     *
     *   @return  Returns a single nondimensional number representing the normalized error
     *            for the calculation of the source term
     */
    virtual double predictorCorrectorGlobalSrcTermErrorNorm();

    virtual void predictorCorrectorPrint(const std::vector<double>& yval,
                                         double pnormSrc, double pnormSoln) const;

    //! Check to see that the preceding step is a successful one
    /*!
     *  (virtual from Electrode_Integrator)
     *
     *   We check to see if the preceding step is a successful one.
     *
     *  @return Returns a bool true if the step is acceptable, and false if it is unacceptable.
     */
    virtual bool  checkSubIntegrationStepAcceptable() const;

    //! Possibly change the solution due to phase births and deaths.
    /*!
     *   (virtual from Electrode_Integrator)
     *
     *  @return  Returns true if the solution step is bad. It returns false if there is not a problem.
     */
    virtual bool changeSolnForBirthDeaths();


    // -----------------------------------------------------------------------------------------------------------------




    //! Print conditions of the electrode for the current integration step to stdout
    /*!
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values
     */
    virtual void printElectrode(int pSrc = 1, bool subTimeStep = false);


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

    //! Returns the molar volume of the Electrode given the extent of reaction
    /*!
     *  The molar volume model for this Electrode is different than other Electrode objects.
     *  It depends solely on the relativeExtentRxn variable, as it must
     *  Within each region the molar volume is calculated as a linear interpolation of the
     *  values on the region boundaries, which are supplied from the input file.
     *  We neglect the temperature dependence here.
     *
     *  The total volume of solid material is then given by the following relation:
     *
     *           SolidVol() =  molarVol *  spMoles_FeS2_Normalization_
     *
     *  Note that spMoles_FeS2_Normalization_ stays constant throughout the calculation once
     *  the size of the electrode is established.
     */
    double molarVolume_relExtentRxn(double relativeExtentRxn) const;

    /*************************************************************************************************************************
     *  OPEN CIRCUIT VOLTAGE
     *************************************************************************************************************************/

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
     *  (virtual function from Electrode.h)
     *
     *   @param isk         Reacting surface domain id
     *   @param iReaction   Explicit index of the reaction. If -1, then it attempts
     *                      to pick the reaction that best represents the open circuit potential.
     */
    virtual double openCircuitVoltageSSRxn(int isk, int iReaction = -1) const;
    

    //! Calculate the open circuit voltage at the final state
    double openCircuitVoltageSS_final() const;

    //! Calculate the inner radius at the final state
    /*!
     *  An inner radius is needed for diffusion approximations. Here we calculate the inner radius
     *  by the extent of reaction carried out in the current region. Everything is relative
     *  to the exterior radius.
     *
     *  When the extent of reaction in the current plateau is zero, we assume that the inner
     *  radius is equal to the external radius (mult by (1 -small)). When the extent of reaction is
     *  nearly equal to the end of the region, we assume the internal radius is nearly zero
     *     (radius_ext * (small))
     *
     *  Here small is defined as 1.0E-8 (will experiment with that number.
     */
    double calculateRadiusInner(double relativeExtentRxn) const;

    //! Returns the equilibrium OCV for the current conditions (virtual)
    /*!
     *  The final conditions are assumed
     *
     * @param isk  Reacting surface domain id
     *
     * @return Returns the voltage
     */
    virtual double openCircuitVoltage(int isk);

    virtual double openCircuitVoltageRxn(int isk, int iReaction = -1) const;

    //! Returns the equilibrium  open circuit voltage for the current conditions 
    /*!
     *  Returns the open circuit voltage within the current region at the relative
     *  extent of reaction. The function  openCircuitVoltageSS_Region()
     *  is called to calculate the standard state OCV.
     *  Then, the activity coefficient of the lithium in the electrolyte is subtracted off of this result: 
     *
     *   Volts = - DeltaG / F
     *
     *   DeltaG = DeltaG_SS - (muLiP - muLiPSS)
     *
     *   We are assuming that the following reaction is used
     *
     *                           Li+    +   Electron-  +  A  ->   B
     *
     *   Then, the difference between the OCV and the SS OCV is the activity of the Li+ with a stoichiometric
     *   coefficient of -1 applied.
     *
     *   @todo   Understand how this changes when an anode is modeled.
     *
     *
     *   @param relativeExtentRxn      Relative extent of reaction
     *   @param xRegion                Integer indicating the region of the extent of reaction
     *
     *   @return returns the standard state voltage.
     */
    double openCircuitVoltage_Region(double relativeExtentRxn, int xRegion) const;

    //! Returns the equilibrium standard state open circuit voltage for the current conditions 
    /*!
     *  Returns the standard state open circuit voltage within the current region at the relative
     *  extent of reaction. Several checks are carried out. Then, the virtual function  openCircuitVoltageSS_Region_NoCheck()
     *  is called. 
     *
     *   @param relativeExtentRxn      Relative extent of reaction
     *   @param xRegion                Integer indicating the region of the extent of reaction
     *
     *   @return returns the standard state voltage.
     */
    double openCircuitVoltageSS_Region(double relativeExtentRxn, int xRegion) const;

    //! Returns the equilibrium standard state open circuit voltage for the current conditions (virtual)
    /*!
     *  Returns the standard state open circuit voltage within the current region at the relative
     *  extent of reaction. The region parameter takes precedence and no checking is done.
     *
     *
     *                           Volts = - DeltaG_SS / F
     *
     *   We are assuming that the following reaction is used
     *
     *                           Li+    +   Electron-  +  A  ->   B
     *
     *   The standard state deltaG, deltaGSS, is defined as 
     *
     *                           deltaG_SS = umSS_B - (umSS_A + umSS_Li+ + umSS_Electron-) 
     *
     *   Then, the difference between the OCV and the SS OCV is the activity of the Li+ with a stoichiometric
     *   coefficient of -1 applied.
     *
     *   @TODO   Understand how this changes when an anode is modeled.
     *
     *  (virtual function of Electrode_MP_RxnExtent.h)
     *
     *   @param relativeExtentRxn      Relative extent of reaction
     *   @param xRegion                Integer indicating the region of the extent of reaction
     *
     *   @return Returns the standard state voltage.
     */
    virtual double openCircuitVoltageSS_Region_NoCheck(double relativeExtentRxn, int xRegion) const;

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
     *  It will also include all plateaus that are defined by the electrode object.
     *
     *  @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     *
     *  @return returns the theoretical capacity of the electrode in Amp seconds = coulombs.
     */
    virtual double capacity(int platNum = -1) const;

    //! Initial capacity of the electrode in Amp seconds
    /*!
     *  This is the initial capacity of the electrode before any degradation occurs.
     *
     *  @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     */
    virtual double capacityInitial(int platNum = -1) const;

    //! Amount of charge that the electrode that has available to be discharged
    /*!
     *   We report the number in terms of Amp seconds = coulombs
     *
     *   @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
     */
    virtual double capacityLeft(int platNum = -1, double voltsMax = 50.0, double voltsMin = -50.0) const;

    //! Set the relative current capacity discharged
    /*!
     *  (virtual function from Electrode)  This routine is called from the electrode level to set the
     *        state of the Electrode to a relative state of charge.
     *
     *
     * This is roughly equal to the total number of electrons that has been discharged
     * from a fully charged state divided by the total moles of solid species in the electrode
     *
     *  @param relDischargedperMole      Relative value of the discharge per mole. Always goes between 0 and number of electrons
     *                                   per active mole, num
     *                                  0 means that the electrode is fully charged, num means that it is fully discharged.
     *
     *
     *  @param platNum           Plateau number. Default is -1 which treats all plateaus as a single entity and
     *                            the relative discharged as a single combined fraction. If platNum is
     *                            >= 0, then the discharge is relative to the current plateau.
     */
    virtual void setRelativeCapacityDischargedPerMole(double relDischargedPerMole, int platNum = -1);


    // -------------------------------  SetState Functions -------------------------------------------------------


    //!   Set the current capacity discharged in amp seconds
    /*!
     * This is roughly equal to the total number of electrons that has been discharged
     * from a fully charged state.
     *
     *  @param  relativeExtentRxn  Relative extent of reaction variable (input)
     */
    void setState_relativeExtentRxn(double relativeExtentRxn);

    //! Set the internal initial intermediate and initial global state from the internal final state
    /*!
     *  (non - virtual function -  onionize in-first)
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
     *  (virtual function from Electrode)
     *
     *  Set the final state from the init state. This is commonly called during a failed time step
     *
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

    //! Set the internal initial intermediate and initial global state from the internal final_final state
    /*!
     *  (virtual function)
     *
     *  Set the intial  and init_int state and the final_final from the final state. We also can set the init_init state from this
     *  routine as well.
     *
     */
    virtual void setInitInitStateFromFinalFinal();

    //! Set the internal final global state from the internal final intermediate state
    /*!
     *  (virtual function from Electrode)
     *
     *  Set the final_final state from the final state. This is commonly called at the end of successful base integration
     */
    virtual void setFinalFinalStateFromFinal();

    //! This is used to set Phase information that is implicit but not set by a restart or an initialization
    /*!
     *  (virtual function from Electrode)
     *
     *  Extra informat that may be needed in advance of a successful updateState() call that specifies all of the
     *  information in the state
     *
     *  @param flagErrors If true any changes in the current flags caused by a mismatch between the state
     *                    and the values of the flags will cause an error exit.
     */
    virtual bool stateToPhaseFlagsReconciliation(bool flagErrors);



    //! Set voltage vs extent of reaction for FeS2
    /*!
     *  hard coded for FeS2 electrode, this sets paramaters such as RegionBoundaries_ExtentRxn_
     *
     */
    int setVoltageVsExtent_FeS2();

    //! Find region
    int findRegion(double RelativeExtentRxn) const;

    //! Return the relative extent of reaction
    double relativeExtentRxn(double time) const;


    //! Return the number of extra print tables
    virtual int getNumPrintTables() const;

    //! Get the values that are printed in tables for the 1D code.
    /*!
     *   @param itable    table id
     *   @param colNames   string names of the header (length is the length of the column)
     *   @param colValues    Value of the columns (length is the length of the column)
     */
    virtual void getPrintTable(int itable, std::vector<std::string>& colNames,
                               std::vector<double>& colValues) const;


    //! Possibly change the solution due to phase births and deaths after phase has been accepted.
    /*!
     *   (virtual from Electrode_Integrator)
     *
     *  This routine is carried out after the step is deemed a success. Massaging of the solution
     *  must be carried out within strict tolerances.
     */
    virtual void manageBirthDeathSuccessfulStep();


    //! Error check on the routine step
    /*!
     *    (virtual from Electrode_Integrator)
     *
     *   Error checks go here. All errors are fatal exits.
     */
    virtual void check_final_state();

    //! Check the nonlinear residual equations for completeness and the ability to be solved
    /*!
     *   @return  0 Everything is good
     *           -1 residual isn't good. We need to cut the time step and retry again.
     */
    virtual int check_nonlinResidConditions();

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
     *
     * @return
     */
    virtual int evalResidNJ(const doublereal t, const doublereal delta_t,
                            const doublereal* const y,
                            const doublereal* const ydot,
                            doublereal* const resid,
                            const ResidEval_Type_Enum evalType = Base_ResidEval,
                            const int id_x = -1,
                            const doublereal delta_x = 0.0);

    //! Fill in the initial conditions
    /*!
     * Values for both the solution and the value of ydot may be provided.
     *
     * @param t0            Time                    (input)
     * @param y             Solution vector (output)
     * @param ydot          Rate of change of solution vector. (output)
     */
    virtual int getInitialConditions(const doublereal t0, doublereal* const y,
                                     doublereal* const ydot);




    //!  Return a vector of delta y's for calculation of the numerical Jacobian
    /*!
     *   There is a default algorithm provided.
     *
     *        delta_y[i] = atol[i] + 1.0E-6 ysoln[i]
     *        delta_y[i] = atol[i] + MAX(1.0E-6 ysoln[i] * 0.01 * solnWeights[i])
     *
     * @param t             Time                    (input)
     * @param y             Solution vector (input, do not modify)
     * @param ydot          Rate of change of solution vector. (input, do not modify)
     * @param delta_y       Value of the delta to be used in calculating the numerical jacobian
     * @param solnWeights   Value of the solution weights that are used in determining convergence (default = 0)
     *
     * @return Returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *            0  Means an unsuccessful operation
     */
    virtual int calcDeltaSolnVariables(const doublereal t, const doublereal* const ySoln,
                                       const doublereal* const ySolnDot, doublereal* const deltaYSoln,
                                       const doublereal* const solnWeights);

    //! Apply a filtering process to the step
    /*!
     *  @param timeCurrent    Current value of the time
     *  @param ybase          current value of the solution
     *  @param step0          Value of the step in the solution vector that will be filtered.
     *                        The filter is applied to the step values.
     *
     *  @return Returns the norm of the value of the amount filtered
     */
    virtual doublereal  filterNewStep(const doublereal timeCurrent,
                                      const doublereal* const ybase,
                                      doublereal* const step0);





protected:


    //! Major internal variable is the Extent
    int numRegions_;

    //! Initial relative extent at the start of the current step
    /*!
     *  units = unitless
     */
    double RelativeExtentRxn_init_;

    //! Initial relative extent at the start of the global step
    double RelativeExtentRxn_init_init_;

    //! final extent of reaction at the end of the current step
    /*!
     *  A note about what relative extent is.
     *
     *    Relative extent is the amount of reaction per mole of reactant. It extends across
     *    multiple plateaus. For graphite, it would correspond to the amount of Lithium per 6
     *    moles of carbon.
     *    units = unitless
     */
    double RelativeExtentRxn_final_;

    //! final extent of reaction at the end of the global step
    double RelativeExtentRxn_final_final_;

    //! temporary variable containing current estimate of the final extent of reaction
    double  RelativeExtentRxn_tmp_;

    //! Initial region for the extent of reaction
    int xRegion_init_;

    //! Initial region for the extent of reaction of the global step
    int xRegion_init_init_;

    //! Final region for the extent of reaction
    int xRegion_final_;

    //! Final region for the extent of reaction of the global step
    int xRegion_final_final_;

    //! The time step variable, this is an unknown in the solution vector
    //   double deltaTsubcycleCalc_;


    std::vector<double> RegionBoundaries_ExtentRxn_;

    //! Molar source rate for the extent of reaction in the electrode object
    //! for the final time during a time step
    /*!
     *  units kmol s-1
     */
    double SrcDot_ExtentRxn_final_;

    //! Source rate for the relative extent of reaction in the electrode object
    //! for the final time during a time step
    /*!
     *  units = unitless
     */
    double SrcDot_RelativeExtentRxn_final_;

    double deltaTdeath_;

    //! Vector containing the rates of progress of reactions defined on the surfaces
    /*!
     *   This concept only works when there is one surface that is reacting
     *   units =   kmol m-2 s-1
     *
     *  Length is over the number of reactions that are defined on that surface, numRxns_[isurf];
     */
    std::vector<double> ROP_;

    //! Molar source rate for the species vector of all species in the electrode object
    //! for the final time during a time step
    /*!
     *  units kmol s-1
     */
    std::vector<double> DspMoles_final_;

    std::vector<int> phaseIndexSolidPhases_;
    std::vector<int> numSpecInSolidPhases_;

    ThermoPhase* Li_liq_;

    //!  This is the heat of formation of B such that the OCV at 723.15 K for the first transition
    //!  of the FeS2 electrode is 2.05565 volts vs. liquid lithium.
    double Hf_B_base_;

    double Hf_B_current_;

    //! Standard State open circuit voltage value
    double volts_OCV_SS_final_;

    //! Current open circuit voltage value
    double volts_OCV_final_;

    //! Normalizing quantity for the depth of discharge
    double spMoles_FeS2_Normalization_;

    //! This integer describes if the system is current on a Region boundary at the start of a subgrid
    //! integration step
    /*!
     *  If it isn't, we set it to -1. If it is, we set it to the region boundary index
     */
    int onRegionBoundary_init_;

    //! This integer describes if the system is current on a Region boundary at the end of a subgrid
    //! integration step
    /*!
     *  If it isn't, we set it to -1. If it is, we set it to the region boundary index
     */
    int onRegionBoundary_final_;

    //!  This Boolean is true when we are at a plateau boundary and the voltage is inbetween the
    //!  top value and the bottom value.
    int goNowhere_;

    //! Molar volume of the electrode at the final state
    /*!
     *  units m3 kmol-1
     */
    doublereal molarVolume_final_;

    //! Molar volume of the electrode at the final global step state
    /*!
     *  units m3 kmol-1
     */
    doublereal molarVolume_final_final_;

    //! Molar volume of the electrode at the initial subgrid step state
    /*!
     *  units m3 kmol-1
     */
    doublereal molarVolume_init_;

    //! Molar volume of the electrode at the initial global step state
    /*!
     *  units m3 kmol-1
     */
    doublereal molarVolume_init_init_;

    //! Internal Radius at the final state (m)
    /*!
     *   In order to do diffusion we need an inner radius
     *   This is defined in terms of the Degree that we have progressed along the extone of reaction variable
     *   between region boundaries.
     */
    doublereal Radius_internal_final_;

    //! Internal radius at the final final state (m)
    doublereal Radius_internal_final_final_;

    //! Internal radius at the initial state
    doublereal Radius_internal_init_;

    //! Internal radius at the initial initial state
    doublereal Radius_internal_init_init_;

    //! Location of the reacting surface
    /*!
     *   0 = Surface area is fixed at the initial conditions. (default)
     *   1 = Surface is located internal to the particle. Surface area of the reacting surface is calculated as moving
     *       front dependent on the extent of reaction during the current plateau
     *   2 = Surface is located at the exterior of the particle, which may be growing.
     */
    int locationOfReactingSurface_;

    //!  Type of modification of the rate of progress due to solid state diffusion
    /*!
     *   0 = none
     *   1 = Inner surface Damkoeler number
     *   2 = Outer surface Damkoeler number
     */
    int ROPModificationType_;

    //! Index of the one reacting surface in the surface list
    /*!
     *  0 interior surface
     *  1 exterior surface
     */
    int indexOfReactingSurface_;

    //! Solid state diffusion model identification
    /*!
     *   0 = There is no diffusion model (default)
     *   1 = There is solid state diffusion limitations. The model assumes pseudo-steady state
     *       The limiting reacting surface location is determined by the  locationOfReactingSurface_ model.
     *       The details of the model is given in the writeup.
     */
    int solidDiffusionModel_;

    //! Diffusion coefficient in the outer region of a two region spherical model
    /*!
     *  We identify the solids by Zones, which are related to regions. Generally, when the extent of reaction is
     *  in the zeroeth region, we are using diffusionCoeffRegions_[1], i.e., the first zone, because that is
     *  zone 1 in the model. Zone 0 would be diffusion in the inner core, whose value we never
     *  actually need. It's actually conceptually simpler this way.
     *
     *  Length = numRegions_ + 1
     *  Units = m**2 s-1
     */
public:
    std::vector<double> diffusionCoeffRegions_;

    //! This is the data field for the calculation of molar volumes
    /*!
     *  Length = numRegions_ + 1
     *  Units = m**3 / kmol
     */
    std::vector<double> molarVolumeRegions_;

    //! Reaction multiplier for the regions
    std::vector<double> rxnPerturbRegions_;

protected:
    //! Equilibrium activities of interstitials of Zone A with Zone A-1
    /*!
     *  We identify the solids by Zones, which are related to regions. Generally, when the extent of reaction is
     *  in the zeroeth region, we are using diffusionCoeffRegions_[1], i.e., the first zone, because that is
     *  zone 1 in the model. Zone 0 would be diffusion in the inner core, whose value we never
     *  actually need. It's actually conceptually simpler this way.
     *
     *  Length = numRegions_ + 1
     *  Units = unitless
     */
    std::vector<double> actEquilInterstitialsRegions_;

    int kf_id_;
    int kf_dir_;

    //
    //        Diffusion model parameters -> The nomenclature is from the memo
    //

    //! Rate constant used in the calculation of the outer  Da #. It must have units of kmol m-2 s-1.
    /*!
     *  This is the rate constant on the solid side of the outer boundary. The reaction is equal to
     *
     *    ROP_f_ext =   kfOuter_  * ca_Li_Interstitial
     */
    double kfExt_;

    //! Rate constant used in the calculation of the outer  Da #. It must have units of kmol m-2 s-1.
    /*!
     *  This is the rate constant on the electrolyte side of the outer boundary. The reaction is equal to
     *
     *    ROP_r_ext =   krOuter_  * ca_Li+ * ca_e-
     */
    double krExt_;

    //! Rate constant used in the calculation of the inn Da #. It must have units of kmol m-2 s-1.
    /*!
     *  This is the rate constant on the solid outer side of the inner boundary. The reaction is equal to
     *
     *    ROP_f_inner =   kfInner_  * ca_Li_Interstitial
     */
    double kfInner_;

    //! Rate constant used in the calculation of the inner Da #. It must have units of kmol m-2 s-1.
    /*!
     *  This is the rate constant on the inner side of the internal boundary. The reaction is equal to
     *
     *      ROP_r_inner = krInner_
     */
    double krInner_;

    //! Damkoeler number calculated assuming the reaction occurs on the outer surface, and we have interstitial diffusion
    //! of a neutral diffusing through a region to a new reaction front, whose reactions are fast so that they are
    //! in equilibrium.
    double DaOuter_;

    //! Damkoeler number divided by r_in calculated assuming the reaction occurs on the outer surface, and we have interstitial diffusion
    //! of a neutral diffusing through a region to a new reaction front, whose reactions are fast so that they are
    //! in equilibrium.
    /*!
     *   This number is equal to
     *
     *        DaOuter_Bar_ = DaOuter_ / Radius_inner_final_
     */
    double DaOuter_Bar_;
    double radiusSmallBound_;

    //! Damkoeler number calculated assuming the reaction occurs on the inner surface, and we have interstitial diffusion
    //! of a neutral diffusing through a region to a new reaction front, whose reactions are slow so that they are
    //! not in equilibrium. Reaction on outer boundary is assumed to be fast compared to inner reaction and to diffusion.
    double DaInner_;

    //! Value of the mole fraction of the lithium interstitial at the inner surface
    /*!
     *   dimensionless
     */
    double Lin_;

    //! Value of the mole fraction of the lithium interstitial at the outer surface
    /*!
     *   dimensionless
     */
    double Lout_;

    //! Forward reaction rate constant for the Lumped reaction
    double ROP_AB_;
    double kfAB_;
    double krAB_;
    double ca_Lip_;
    double Eocv_;

    double betaF_AB_;
    double betaR_AB_;

    //! This is 0 or 1
    /*!
     *       0 means that it is normal
     *       1 means that Da_outer goes to inf as r_in goes to zero.
     */
    int limitingEquationBehavior_;

    //!
    /*!
     *  See notes
     *     1
     *     2
     */
    int innerDaTreatmentType_;



    //! Absolute tolerance for the integrated global src term vectors
    //  std::vector<double> atol_IntegratedSrc_global_;
public:

    //! Hook for the nonlinear solver
    //   Cantera::NonlinearSolver *pSolve_;

    //  calcRate_ResidJacEval  *pSolve_Res_;


    // -----------------------------------------------------------------------------------------------------------------

    friend class EState;

};

}


#endif
/*****************************************************************************/
