/**
 * @file m1d_porousFlow_dom1D.h Base class for calculating residuals within a porous flow domain
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef M1D_POROUSFLOW_DOM1D_H_
#define M1D_POROUSFLOW_DOM1D_H_

#include "m1d_defs.h"
#include "m1d_BulkDomain1D.h"
#include "m1d_BDD_porousFlow.h"

#include "m1d_valCellTmps_porousFlow.h"

#include "zuzax/transport.h"

namespace Zuzax
{
class Transport;
}
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
class LocalNodeIndices;
class cellTmps;
class ExtraPhase;

//==================================================================================================================================
//! This is derived class provides the function evaluation for a porous electrode.
/*!
 *  The porous electrolyte domain is characterized by a current conservation equation and several species
 *  conservation equations describing the electrolyte. A porosity/tortuosity is also associated with the domain.
 *
 *  There is a 1 to 1 mapping between the local control volume indexing and the Local Consecutive Ordering indexing
 *
 *  There is a 1 to 1 mapping between the global control volume indexing and the Global node number indexing
 *  that is given by a single offset.
 */
class porousFlow_dom1D : public BulkDomain1D
{

public:
    //! Constructor
    /*!
     *  @param[in]           bdd_pf_ptr          Contains the bulk domain description.
     */
    porousFlow_dom1D(m1d::BDD_porousFlow* bdd_pf_ptr);

    //! Copy constructor
    /*!
     *  @param[in]           r                   Object to be copied into the current object
     */
    porousFlow_dom1D(const porousFlow_dom1D& r);

    //! Destructor
    virtual  ~porousFlow_dom1D();

    //! Assignment operator
    /*!
     *  @param[in]           r                   Object to be copied into the current object
     *
     *  @return                                  Returns a changeable reference to the current object
     */
    porousFlow_dom1D& operator=(const porousFlow_dom1D& r);

    //! Prepare all of the indices for fast calculation of the residual
    /*!
     *  (virtual from Domain1D)
     *
     *  Ok, at this point, we will have figured out the number of equations
     *  to be calculated at each node point. The object NodalVars will have  been fully formed.
     *
     *  We use this to figure out what local node numbers/ cell numbers are
     *  needed and to set up indices for their efficient calling.
     *
     *  Child objects of this one will normally call this routine in a recursive fashion.
     *
     *  @param[in]           li_ptr              Pointer to the LocalNodeIndices object
     */
    virtual void domain_prep(LocalNodeIndices* li_ptr) override;

    //! Function that gets called at end the start of every time step
    /*!
     *  (virtual from domain1d)
     *
     *  This function provides a hook for a residual that gets called whenever a
     *  time step has been accepted and we are about to move on to the next time step.
     *  The call is made with the current time as the time that is accepted.
     *  The old time may be obtained from t and rdelta_t_accepted.
     *
     *  After this call interrogation of the previous time step's results will not be valid.
     *
     *  Note, when t is equal to t_old, soln_ptr should equal solnOld_ptr values. However,
     *  solnDot_ptr values may not be zero.
     *
     *
     *  @param[in]           doTimeDependentResid  This is true if we are solving a time dependent problem.
     *  @param[in]           soln_ptr            Solution value at the current time
     *  @param[in]           solnDot_ptr         derivative of the solution at the current time.
     *  @param[in]           solnOld_ptr         Solution value at the old time step, n-1
     *  @param[in]           t                   current time to be accepted, n
     *  @param[in]           t_old               previous time step value
     */
    virtual void
    advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector* const soln_ptr,
                        const Epetra_Vector* solnDot_ptr, const Epetra_Vector* const solnOld_ptr,
                        const double t, const double t_old) override;

    //! Utility function to calculate quantities before the main residual routine.
    /*!
     *  (virtual from Domain1D)
     *
     *  This is used for a loop over nodes in the domain. All calculated quantities must be internally storred within the
     *  domain structure. Currently this is called during the residual evalultion of the problem.
     *
     *  @param[in]             doTimeDependentResid  Boolean indicating whether the time dependent residual is requested
     *  @param[in]             soln_ptr            Solution vector at which the residual should be evaluated
     *  @param[in]             solnDot_ptr         Solution dot vector at which the residual should be evaluated.
     *  @param[in]             solnOld_ptr         Pointer to the solution vector at the old time step
     *  @param[in]             t                   time
     *  @param[in]             rdelta_t            inverse of delta_t
     *  @param[in]             residType           Type of evaluation of the residual. Uses the ResidEval_Type_Enum type.
     *                                             Defaults to Base_ResidEval
     *  @param[in]             solveType           Type of solution Type. Uses the Solve_Type_Enum  type.
     *                                             Defaults to  TimeDependentAccurate_Solve
     */
    virtual void residEval_PreCalc(const bool doTimeDependentResid, const Epetra_Vector* const soln_ptr,
                      const Epetra_Vector* const solnDot_ptr, const Epetra_Vector* const solnOld_ptr,
                      const double t, const double rdelta_t,
                      const Zuzax::ResidEval_Type residType, const Zuzax::Solve_Type solveType) override;

    //! Get the state vector at a particular cell in the domain
    /*!
     *  The state here is defined as the following vector:
     *      -  0         temperature
     *      -  1         pressure
     *      -  2         mole fractions of the electrolyte species
     *      -  2 + nsp   electric potential of electrolyte
     *
     *  @param[in]           iCell               Cell number
     *  @param[in]           soln                Reference to the Ghosted Epetra_Vector containing the solution
     *
     *  @param[out]          state_Lyte          Reference to the return vector which will contain the state 
     */
    void getState_Lyte_atCell(size_t iCell, const Epetra_Vector_Ghosted& soln, std::vector<double>& state_Lyte);


    //! Get the state vector at a particular node in the domain
    /*!
     *  The state here is defined as the following vector:
     *      -  0         temperature
     *      -  1         pressure
     *      -  2         mole fractions of the electrolyte species
     *      -  2 + nsp   electric potential of electrolyte
     *
     *  @param[in]           nv                  NodalVars structure for the current node
     *  @param[in]           solnNode_Curr       Pointer to the start solution vector at the current node
     *
     *  @param[out]          state_Lyte          Reference to the return vector which will contain the state 
     */
    void getState_Lyte(const NodalVars* const nv, const double* const solnNode_Curr, std::vector<double>& state_Lyte);

    //! Setup shop at a particular nodal point in the domain, calculating intermediate quantites and updating Zuzax's objects
    /*!
     *  (virtual porousFlow_dom1D)
     *
     *  All member data with the suffix, _Curr_, are updated by this function.
     *
     *  Calculated quantities:
     *             porosity_Cell_[cIndex_cc_]
     *             porosity_Curr_ 
     *  Called functions:
     *             updateElectrolyte()
     *
     *  @param[in]           nv                  NodalVars structure for the current node
     *  @param[in]           solnNode_Curr       Pointer to the start solution vector at the current node
     */
    virtual void SetupThermoShop1(const NodalVars* const nv, const double* const solnNode_Curr);

    //! Function updates the ThermoPhase object for the electrolyte given the solution vector
    /*!
     *  (virtual porousFlow_dom1D)
     *
     *  This function calculates the values at the node
     *
     *  Calculated quantities:
     *              temp_Curr_
     *              pres_Curr_
     *              concTot_Curr_
     *  Called functions:
     *              getMFElectrolyte_soln()
     *              getVoltages()
     *  Sets the ThermoPhase state of the electrolyte, ionicLiquid_, to the node value.
     *
     *  @param[in]           nv                  Nodal Values for the current node
     *  @param[in]           solnNode_Curr       Pointer to the start solution vector at the current node
     */
    virtual void updateElectrolyte(const NodalVars* const nv, const double* const solnNode_Curr);

    //! Function extracts the voltages
    /*!
     *  (virtual porousFlow_dom1D)
     *
     *  This function extracts the voltage values at the node.
     *
     *  Calculated quantities:
     *              phiElectrolyte_Curr_
     *
     *  @param[in]           nv                  Nodal Values for the current node
     *  @param[in]           solnNode_Curr       Pointer to the start solution vector at the current node
     */
    virtual void getVoltages(const NodalVars* const nv, const double* const solnNode_Curr);

    //! Function extracts the mole fractions of the electrolyte species
    /*!
     *  (virtual porousFlow_dom1D)
     *
     *  This function extracts the mole fractions of the electrolyte species
     *
     *  Calculated quantities:
     *              mfElectrolyte_Soln_Curr_[k]
     *              mfElectrolyte_Thermo_Curr_[k]
     *
     *  @param[in]           nv                  Nodal Values for the current node
     *  @param[in]           solnNode_Curr       Pointer to the start solution vector at the current node
     */
    virtual void getMFElectrolyte_soln(const NodalVars* const nv, const double* const solnNode_Curr);

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
     *
     *  @param[in]           mf                  Input mole fractions
     *                                              Length: m_kk = number of species in electrolyte phase)
     *  @param[out]          mf_Thermo           Output mole fractions
     *                                              Length: m_kk = number of species in electrolyte phase)
     */
    void calcMFElectrolyte_Thermo(const double* const mf, double* const mf_Thermo) const;

    //! Calculate the heat source
    /*!
     *  (virtual from porousFlow_dom1D)
     *
     *  Sums up all of the quantities in qSource_Cell_curr_[iCell] and returns the result.
     *
     *  @return                                  Returns the sum of the heat sources within the domain
     *                                              Units: Joules / m2
     */
    virtual double heatSourceLastStep() const;

    //! Calculate the accumulated heat source from all previous steps (Joules /m2)
    /*!
     *  (virtual from porousFlow_dom1D)
     *   
     *  Sums up all of the quantities in qSource_Cell_accumul_[iCell] and returns the result.
     *
     *  @return                                  Returns the sum of the heat sources within the domain
     *                                             Units: Joules / m2
     */
    virtual double heatSourceAccumulated() const;

    //! Zero the accumulated heat source terms
    /*!
     *  (virtual from porousFlow_dom1D)
     */
    virtual void heatSourceZeroAccumulated() const;

    //! Set up tmps for quick calculation of residuals
    void residSetupTmps();

    //! Calculates and returns an estimate of the effective areal resistance of the layer
    /*!
     *  (virtual in porous_flow_dom1D)
     *
     *   resistance = ((potCathodic - potAnodic) - voltOCV) / current
     *
     *  @param[in]            potAnodic          potential in the anodic direction. If the anode, this returns the potential of the
     *                                           solid in the anode next to the anode current collector.
     *  @param[in]            potCathodic        potential in the cathode direction. If the anode, this returns the potential of the
     *                                           electrolyte in the anode next to the separator.
     *  @param[in]            voltOCV            OCV calculated in a quick manner.
     *  @param[in]            current            Returns the total current going through the domain ( amps m-2)
     *
     *  @return                                  returns the effective resistance of the layer (ohm m2)
     */
    virtual double effResistanceLayer(double& potAnodic, double&  potCathodic, double& voltOCV, double& current);

    //! Generate the initial conditions
    /*!
     *  (virtual from Domain1D)
     *
     *  The basic algorithm is to loop over the volume domains. Then, we loop over the surface domains
     *
     *  @param[in]        doTimeDependentResid   Boolean indicating whether we should
     *                                           formulate the time dependent residual
     *  @param[in]        soln                   Solution vector. This is the input to the residual calculation.
     *  @param[in]        solnDot                Solution vector. This is the input to the residual calculation.
     *  @param[in]        t                      Time
     *  @param[in]        delta_t                Time step to be used in the initial time step
     */
    virtual void initialConditions(const bool doTimeDependentResid, Epetra_Vector* soln, Epetra_Vector* solnDot,
                                   const double t, const double delta_t) override;

    //! Calculate the thermal conductivity of the porous matrix at the current cell.
    /*!
     *  (virtual in porous_flow_dom1D)
     *
     *  Currently it is hardwired to return 0.5 watts m-1 K-1, a number characteristic of water.
     *
     *  @return                                   Returns the thermal conductivity of the effective porous media
     *                                              Units = watts m-1 K-1
     */
    virtual double thermalCondCalc_PorMatrix();

    //! Volume fraction of other components than the Electrode object and the electrolyte phase
    /*!
     *  (virtual in porous_flow_dom1D)
     *
     *  Calculates the volume fraction of phases other than the electrolyte and the electrode at a given cell
     *  Sums up the phase volumeFraction_Phases_Cell_[] function.
     *
     *  @param[in]           iCell               Cell number
     *
     *  @return                                  Returns the volume fraction 
     */
    virtual double volumeFractionOther(size_t iCell);

    //! Calculate the porosity
    /*!
     *  (virtual in porous_flow_dom1D)
     *
     *  @param[in]           Cell               Cell number
     *
     *  @return                                 Returns the porosity
     */
    virtual double calcPorosity(size_t iCell);

    // -----------------------------------------------------------------------------------------------------------------------------
    // -----------------------------------------   DATA   --------------------------------------------------------------------------
    // -----------------------------------------------------------------------------------------------------------------------------

    //! Pointer to the BDD object that is most derived
    BDD_porousFlow* BDD_ptr_;
   
    // ------------------  Options for Processing ------------------------------------------------------------------------



    // ------------------- Thermodynamics quantities on the domain -------------------------------------------------------

    //! Molar Heat Capacity of the electrolyte phase located in all of the cells
    /*!
     *   Vector of Molar heat capacity const press (iCell)
     *   Units of Joules/(K)
     */
    std::vector<double> CpMolar_lyte_Cell_;

    std::vector<double> CpMolar_solid_Cell_;
    std::vector<double> CpMolar_total_Cell_;

    //!  Partial molar Enthalpy  of the electrolyte species located in all of the cells
    /*!
     *   Vector of partial molar enthalpy  (KRSpecies, iCell)
     *   Units of Joules/(kmol)
     */
    std::vector<double> EnthalpyPM_lyte_Cell_;

    //!  Molar Enthalpy of the electrolyte phase located in all of the cells
    /*!
     *   Vector of molar enthalpy  (iCell)
     *   Units of Joules/(kmol)
     */
    std::vector<double> EnthalpyMolar_lyte_Cell_;

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
     *  We also ensure that the charge neutrality constraint is satisfied.  We ensure that the mole fractions sum to one.
     *  Thermo is not defined if any of the three are violated.
     */
    std::vector<double> mfElectrolyte_Thermo_Curr_;

    //!  Partial molar Enthalpy of the electrolyte species at the Curr point with Phi enhancement
    /*!
     *   Vector of partial molar enthalpy  (KRSpecies)
     * 
     *   Units:   Joules / kmol
     */
    std::vector<double> EnthalpyPM_lyte_Curr_;
    std::vector<double> EnthalpyPhiPM_lyte_Curr_;

    //! Value of the molar Enthalpy of the electrolyte at the current location
    /*!
     *  Units:    Joules / kmol
     */
    double EnthalpyMolar_lyte_Curr_;
 
    //! Source of heat during the current time step
    /*!
     *  qSource is the heat generated in each cell per time for the current time step.
     *  Therefore, the source term is integrated in the axial direction and it is integrated wrt to the time interval.
     *
     *  Length:   number of local cells
     *  Units:    Joules / m2
     */
    mutable std::vector<double> qSource_Cell_curr_;

    //! Source of heat during the global interval
    /*!
     *  qSource is the heat generated in each cell per time and accumulated over
     *  the time steps that make up the current interval.
     *  Therefore, the source term is integrated in the axial direction and it is integrated wrt to the time interval.
     *
     *  Length:   number of local cells
     *  Units:    Joules / m2
     */
    mutable std::vector<double> qSource_Cell_accumul_;

    //! Source of joule heating in the electrolyte during the current time interval for the current cell
    /*!
     *  This is the heat generated due to joule heating in the electrolyte in each cell per time for the current time step
     *  Therefore, the source term is integrated in the axial direction and it is
     *  integrated wrt to the time interval.
     *
     *  Length:  number of local cells
     *  Units:   Joules / m2
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
    Zuzax::VelocityBasis ivb_;

    //! Pointer to the thermo object for the electrolyte
    /*!
     *   We do not own this object
     */
    Zuzax::ThermoPhase* ionicLiquid_;

    //! Pointer to the transport object for the electrolyte
    Zuzax::Transport* trans_;

    //! Pointer to the solid skeleton
    Zuzax::ThermoPhase* solidSkeleton_;

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
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
