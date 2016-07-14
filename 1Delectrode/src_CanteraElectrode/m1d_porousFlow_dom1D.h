/**
 * @file m1d_porousFlow_dom1D.h
 */


#ifndef M1D_POROUSFLOW_DOM1D_H_
#define M1D_POROUSFLOW_DOM1D_H_

#include "m1d_defs.h"
#include "m1d_BulkDomain1D.h"
#include "m1d_BDD_porousFlow.h"

#include "m1d_valCellTmps_porousFlow.h"

#include <cantera/transport.h>

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
class Transport;
}

namespace m1d
{
class LocalNodeIndices;
class cellTmps;
class ExtraPhase;

//======================================================================================================================
//! This is derived class  provides the function
//! evaluation for a porous electrode.
/*!
 * The porous electrolyte domain is characterized by a
 * current conservation equation and several species
 * conservation equations describing the electrolyte.
 * A porosity/tortuosity is also associated with the domain.
 *
 * There is a 1 to 1 mapping between the local control volume indexing and
 * the Local Consecutive Ordering indexing
 *
 * There is a 1 to 1 mapping between the global control volume indexing
 * and the Global node number indexing that is given by a single offset.
 */
class porousFlow_dom1D : public BulkDomain1D
{

public:

    //! Constructor
    /*!
     * @param bdd   Contains the bulk domain description.
     */
    porousFlow_dom1D(m1d::BDD_porousFlow* bdd_pf_ptr);

    //! Copy constructor
    /*!
     * @param r      Object to be copied into the current object
     */
    porousFlow_dom1D(const porousFlow_dom1D& r);

    //! Destructor
    virtual  ~porousFlow_dom1D();

    //! Assignment operator
    /*!
     * @param r      Object to be copied into the current object
     * @return       Returns a changeable reference to the current object
     */
    porousFlow_dom1D& operator=(const porousFlow_dom1D& r);

    //! Prepare all of the indices for fast calculation of the residual
    /*!
     *  Ok, at this point, we will have figured out the number of equations
     *  to be calculated at each node point. The object NodalVars will have
     *  been fully formed.
     *
     *  We use this to figure out what local node numbers/ cell numbers are
     *  needed and to set up indices for their efficient calling.
     *
     *  Child objects of this one will normally call this routine in a
     *  recursive fashion.
     */
    virtual void
    domain_prep(LocalNodeIndices* li_ptr);

    //! Function that gets called at end the start of every time step
    /*!
     *  This function provides a hook for a residual that gets called whenever a
     *  time step has been accepted and we are about to move on to the next time step.
     *  The call is made with the current time as the time
     *  that is accepted. The old time may be obtained from t and rdelta_t_accepted.
     *
     *  After this call interrogation of the previous time step's results will not be
     *  valid.
     *
     *   @param  doTimeDependentResid  This is true if we are solving a time dependent
     *                                 problem.
     *   @param  soln_ptr              Solution value at the current time
     *   @param  solnDot_ptr           derivative of the solution at the current time.
     *   @param  solnOld_ptr           Solution value at the old time step, n-1
     *   @param  t                     current time to be accepted, n
     *   @param  t_old                 previous time step value
     */
    virtual void
    advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector* soln_ptr,
                        const Epetra_Vector* solnDot_ptr, const Epetra_Vector* solnOld_ptr,
                        const double t, const double t_old);

    virtual void
    residEval_PreCalc(const bool doTimeDependentResid,
                      const Epetra_Vector* soln_ptr,
                      const Epetra_Vector* solnDot_ptr,
                      const Epetra_Vector* solnOld_ptr,
                      const double t,
                      const double rdelta_t,
                      const ResidEval_Type_Enum residType,
                      const Solve_Type_Enum solveType);

    //!  Setup shop at a particular nodal point in the domain, calculating intermediate quantites
    //!  and updating Cantera's objects
    /*!
     *  All member data with the suffix, _Curr_, are updated by this function.
     *
     * @param nv         NodalVars structure for the current node
     * @param soln_Curr  Current value of the solution vector
     */
    virtual void SetupThermoShop1(const NodalVars* const nv, const doublereal* const soln_Curr);


    //! Function updates the ThermoPhase object for the electrolyte given the solution vector
    /*!
     *  This function calculates the values at the cell center
     *
     * @param nv Nodal Values for the current node
     * @param solnElectrolyte
     */
    virtual void
    updateElectrolyte(const NodalVars* const nv, const doublereal* const solnElectrolyte);

    virtual void
    getVoltages(const NodalVars* const nv, const double* const solnElectrolyte);

    virtual void
    getMFElectrolyte_soln(const NodalVars* const nv, const double* const solnElectrolyte);

    //! Calculate a thermodynamically acceptable mole fraction vector from the current raw mole
    //! fraction vector
    /*!
     *  The following conditions are enforced:
     *
     *     charge neutral solution
     *     no negative mole fractions
     *     sums to 1
     *
     *  Note, this routine is really meant for phases with charged species in them. It works
     *  for neutral phases, but is inefficient.
     */
    void calcMFElectrolyte_Thermo(const double* const mf, double* const mf_Thermo) const;

    virtual double heatSourceLastStep() const;

    virtual double heatSourceAccumulated() const;

    virtual void heatSourceZeroAccumulated() const;

    //! Set up tmps for quick calculation of residuals
    void residSetupTmps();

    //! Calculates and returns an estimate of the effective areal resistance of the layer
    /*!
     *   (virtual in porous_flow_dom1D)
     *
     *   resistance = ((potCathodic - potAnodic) - voltOCV) / current
     *
     *  @param potAnodic potential in the anodic direction. If the anode, this returns the potential of the
     *                   solid in the anode next to the anode current collector.
     *  @param potCathodic potential in the cathode direction. If the anode, this returns the potential of the
     *                   electrolyte in the anode next to the separator.
     *  @param voltOCV  OCV calculated in a quick manner.
     *  @param current   Returns the total current going through the domain ( amps m-2)
     *
     *  @return returns the effective resistance of the layer (ohm m2)
     */
    virtual double effResistanceLayer(double& potAnodic, double&  potCathodic, double& voltOCV, double& current);

    //! Generate the initial conditions
    /*!
     *   The basic algorithm is to loop over the volume domains.
     *   Then, we loop over the surface domains
     *
     * @param doTimeDependentResid    Boolean indicating whether we should
     *                                formulate the time dependent residual
     * @param soln                    Solution vector. This is the input to
     *                                the residual calculation.
     * @param solnDot                 Solution vector. This is the input to
     *                                the residual calculation.
     * @param t                       Time
     * @param delta_t                 delta_t for the initial time step
     */
    virtual void
    initialConditions(const bool doTimeDependentResid, Epetra_Vector* soln, Epetra_Vector* solnDot,
                      const double t, const double delta_t);

    //! Calculate the thermal conductivity of the porous matrix at the current cell.
    virtual double thermalCondCalc_PorMatrix();

    //! Volume fraction of other components than the Electrode object and the electrolyte phase
    /*!
     *   Calculates the volume fraction at a given cell.
     *
     *   @param[in]       iCell       Cell number
     *
     *   @return                      Returns the volume fraction 
     */
    virtual double volumeFractionOther(size_t iCell);

    virtual double calcPorosity(size_t iCell);


    // -----------------------------------------------------------------------------------------------------------------------------
    // -----------------------------------------   DATA   --------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------------------

    //! Pointer to the BDD object that is most derived
    BDD_porousFlow* BDT_ptr_;
   
    // ------------------  Options for Processing ------------------------------------------------------------------------



    // ------------------- Thermodynamics quantities on the domain -------------------------------------------------------

    //! Molar Heat Capacity of the electrolyte phase located in all of the cells
    /*!
     *   Vector of Molar heat capacity const press (iCell)
     *   Units of Joules/(K)
     */
    std::vector<doublereal> CpMolar_lyte_Cell_;

    std::vector<doublereal> CpMolar_solid_Cell_;
    std::vector<doublereal> CpMolar_total_Cell_;

    //!  Partial molar Enthalpy  of the electrolyte species located in all of the cells
    /*!
     *   Vector of partial molar enthalpy  (KRSpecies, iCell)
     *   Units of Joules/(kmol)
     */
    std::vector<doublereal> EnthalpyPM_lyte_Cell_;

    //!  Molar Enthalpy of the electrolyte phase located in all of the cells
    /*!
     *   Vector of molar enthalpy  (iCell)
     *   Units of Joules/(kmol)
     */
    std::vector<doublereal> EnthalpyMolar_lyte_Cell_;

    //
    // ------------------- Porosity of the Domain -----------------------------------------------------------------------
    //

protected:
    //! Volume Fraction of the electrolyte within each control volume
    /*!
     * (change to CV)
     *  Length is number of cells on the processor.
     */
    std::vector<double> porosity_Cell_;

    //! Volume Fraction of the electrolyte within the cell at the previous time step
    /*!
     *  Length is number of cells on the processor.
     */
    std::vector<double> porosity_Cell_old_;

    //! Temperature within each cell at the previous time step
    /*!
     *  Length is number of cells on the processor.
     */
    std::vector<double> Temp_Cell_old_;

    //! Mesh Values
    

    //! Control Volume thickenss for all cells
    /*!
     *  Length is number of cells on the processor
     *  units = m
     */
    std::vector<double> xdelCell_Cell_;


    //!  Number of extra condensed phases
    /*!
     *   The first phase refers to the ThermoPhase solidSkeleton. Other phases are problem dependent.
     *   For electrode phases, the solidSkeleton phase is normally the binder.
     */
    size_t numExtraCondensedPhases_;

    //!  Volume fraction of Extra Condensed phases in each phase.
    /*!
     *  The first phase refers to the ThermoPhase solidSkeleton. Other phases are problem dependent.
     *
     *      volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell_ + jPhase]
     */
    std::vector<double> volumeFraction_Phases_Cell_;

    //!  Volume fraction of Extra Condensed phases in each phase at the last time step
    /*!
     *  The first phase refers to the ThermoPhase solidSkeleton. Other phases are problem dependent.
     *
     *      volumeFraction_Phases_Cell_old_[numExtraCondensedPhases_ * iCell_ + jPhase]
     */
    std::vector<double> volumeFraction_Phases_Cell_old_;

    //!  Mole number of each extra phase (kmol / m2)
    /*!
     *  This is a per cross-sectional area quantity. 
     *  The first phase refers to the ThermoPhase solidSkeleton. Other phases are problem dependent.
     *
     *      moleNumber_Phases_Cell_[numExtraCondensedPhases_ * iCell_ + jPhase]
     */
    std::vector<double> moleNumber_Phases_Cell_;

    //!  Mole number of each extra phase at the old time step (kmol / m2)
    /*!
     *  This is a per cross-sectional area quantity. 
     *  The first phase refers to the ThermoPhase solidSkeleton. Other phases are problem dependent.
     *
     *      moleNumber_Phases_Cell_[numExtraCondensedPhases_ * iCell_ + jPhase]
     */
    std::vector<double> moleNumber_Phases_Cell_old_;

    //
    // --------------------- Cell Storage -------------------------------------------------------------------------------
    //
    //!  Cell storage -> storage of cell related quantities at the current cell, cIndex_cc_

    //! Cell index number
    int cIndex_cc_;

    //
    // ------------------- Locally derived quantities that are valid at the point of current interest --------------------
    //
    //                         ( these are intermediate values and all have the suffix _Curr_ )
    //                         ( these all refer to the new time value)
    //                         ( these are all calculated by the routine SetupShopThermo1() or SetupShopThermo2()
    //                         ( these are either node quantities or boundary quantities )

    //! Temperature at the current point (Kelvin)
    double temp_Curr_;

    //! 1d stress at the current point pa
    //    double mm_stress_Curr_;

    //! Local value of the pressure (Pascal)
    double pres_Curr_;

    //!  Total concentration of the electrolyte at the current position (kmol m-3)
    double concTot_Curr_;

    //!  Current value of the electrolyte voltage (volts)
    double phiElectrolyte_Curr_;

    //! Current porosity
    double porosity_Curr_;

    //! Thermal Conductivity at the current point
    /*!
     *  units = Watts / (m K) = Joule / (s m K) = Volt Amp / (m K)
     */
    double thermalCond_Curr_;

    //! Heat flux at the current point
    /*!
     *  Units = Watts / m2
     */
    double heatFlux_Curr_;

    //! Heat flux of the Enhanced EnthalpyPhi at the current point
    /*!
     *  Units = Watts / m2
     */
    double jFlux_EnthalpyPhi_Curr_;

    //! Vector of temporary indexing quantities for each cell
    /*!
     * These are calculated once at the start of the program
     */
    std::vector<cellTmps> cellTmpsVect_Cell_;

    //! Current value of the Electrolyte mole fraction vector - noncropped
    std::vector<double> mfElectrolyte_Soln_Curr_;

    //! Current value of the Electrolyte mole fraction vector - cropped to always be positive
    /*!
     *  We also ensure that the charge neutrality constraint is satisfied.
     *  We ensure that the mole fractions sum to one.
     *  Thermo is not defined if any of the three are violated.
     */
    std::vector<double> mfElectrolyte_Thermo_Curr_;

    //!  Partial molar Enthalpy of the electrolyte species at the Curr point with Phi enhancement
    /*!
     *   Vector of partial molar enthalpy  (KRSpecies)
     *   Units of Joules/(kmol)
     */
    std::vector<doublereal> EnthalpyPM_lyte_Curr_;
    std::vector<doublereal> EnthalpyPhiPM_lyte_Curr_;

    //! Value of the molar Enthalpy of the electrolyte at the current location
    /*!
     *  Units are Joules/kmol
     */
    double EnthalpyMolar_lyte_Curr_;
 
    //! Source of heat during the current time step
    /*!
     *  qSource is the heat generated in each cell per time for the current time step.
     *  Therefore, the source term is integrated in the axial direction and it is
     *  integrated wrt to the time interval.
     *
     *  Length = number of local cells
     *  units = Joules / m2
     */
    mutable std::vector<double> qSource_Cell_curr_;

    //! Source of heat during the global interval
    /*!
     *  qSource is the heat generated in each cell per time and accumulated over
     *  the time steps that make up the current interval.
     *  Therefore, the source term is integrated in the axial direction and it is
     *  integrated wrt to the time interval.
     *
     *  Length = number of local cells
     *  units = Joules / m2
     */
    mutable std::vector<double> qSource_Cell_accumul_;

    //! Source of joule heating in the electrolyte during the current time interval for the current cell
    /*!
     *  This is the heat generated due to joule heating in the electrolyte in each cell per time for the current time step
     *  Therefore, the source term is integrated in the axial direction and it is
     *  integrated wrt to the time interval.
     *
     *  Length = number of local cells
     *  units = Joules / m2
     */
    std::vector<double> jouleHeat_lyte_Cell_curr_;

    //! Source of joule heating in the solid during the current time interval for the current cell
    /*!
     *  This is the heat generated due to joule heating in the electrolyte in each cell per time for the current time step
     *  Therefore, the source term is integrated in the axial direction and it is
     *  integrated wrt to the time interval.
     *
     *  Length = number of local cells
     *  units = Joules / m2
     */
    std::vector<double> jouleHeat_solid_Cell_curr_;

    //! Source of heat release in the electrode during the current time interval for the current cell
    /*!
     *  This is the heat generated due to diffusion and reaction in the electrode
     *  in each cell per time for the current time step
     *  Therefore, the source term is integrated in the axial direction and it is
     *  integrated wrt to the time interval.
     *
     *  Length = number of local cells
     *  units = Joules / m2
     */
    std::vector<double> electrodeHeat_Cell_curr_;

    //! Source of overpotential heat release in the electrode during the current time interval for the current cell
    /*!
     *  This is the heat generated due to irreversible reaction and diffusion processes in the electrode
     *  in each cell per time for the current time step
     *  Therefore, the source term is integrated in the axial direction and it is
     *  integrated wrt to the time interval.
     *
     *  Length = number of local cells
     *  units = Joules / m2
     */
    std::vector<double> overPotentialHeat_Cell_curr_;

    //! Source of reversible entropic heat release in the electrode during the current time interval for the current cell
    /*!
     *  This is the heat generated due to reversible reaction and diffusion processes in the electrode
     *  in each cell per time for the current time step
     *  Therefore, the source term is integrated in the axial direction and it is
     *  integrated wrt to the time interval.
     *
     *  Length = number of local cells
     *  units = Joules / m2
     */
    std::vector<double> deltaSHeat_Cell_curr_;

    //!  Anodic electric potential of the domain for the current location
    /*!
     *  This is the potential that is closest to the anode side of the domain.
     *  If we are in the anode, this is the potential of the elctrode near the anode
     *  current collector.
     *
     *  Quantity is evaluated in post processor
     */
    double potentialAnodic_;

    //!  Cathodic electric potential of the domain
    /*!
     *  This is the potential that is closest to the cathode side of the domain.
     *
     *  Quantity is evaluated in post processor
     */
    double potentialCathodic_;

    //!  Total enthalpy within each cell at the current conditions (Joules/m2)
    /*!
     *   This only depends on the current conditions of temperatue, pressure, mole numbers of species, and volume
     *   fractions. It is an extensive quantity.
     *
     *   This quantity is malloced always. But it is only calculated when doEnthalpyEquation_ is true.
     *
     *   Length = total number of cells.
     */
    std::vector<double> nEnthalpy_New_Cell_;

    //!  Total enthalpy within each cell at the previous time step (Joules/m2)
    /*!
     *   This only depends on the previous time step conditions of temperatue, pressure, mole numbers of species, and 
     *   volume fractions. It is an extensive quantity.
     *
     *   This quantity is malloced always. But it is only calculated when doEnthalpyEquation_ is true.
     *
     *   Length = total number of cells.
     */
    std::vector<double> nEnthalpy_Old_Cell_;

    //! Vector of average thermal conductivities of the matrix comprising each cell
    std::vector<double> thermalCond_Cell_;

    //! Vector of temporary solution variables
    std::vector<valCellTmps> valCellTmpsVect_Cell_;

    //! Velocity basis of the transport equations
    ZZCantera::VelocityBasis ivb_;

    //! Pointer to the thermo object for the electrolyte
    /*!
     *   We do not own this object
     */
    ZZCantera::ThermoPhase* ionicLiquid_;

    //! Pointer to the transport object for the electrolyte
    ZZCantera::Transport* trans_;

    //! Pointer to the solid skeleton
    ZZCantera::ThermoPhase* solidSkeleton_;

    //! Vector of extra phases
    std::vector<ExtraPhase*> ExtraPhaseList_;

    //! Porosity Equation Problem Type
    /*! 
     *   List of available options : 
     *    0 ->  Constant                   
     *    1 ->  Calculated Out Of Equation System                    
     *    2 ->  Calculated in equation system    
     *    3 ->  Calculated in equation System as  part Of Mechanics system
     *          Stress equation causes change in strains causing change of geometry
     */
    int Porosity_prob_type_;

};
//======================================================================================================================
}
//======================================================================================================================
#endif
