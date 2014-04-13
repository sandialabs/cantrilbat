/*
 * m1d_porousLiIon_Cathode_dom1D.h
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 */

#ifndef M1D_POROUSLIION_CATHODE_DOM1D_H_
#define M1D_POROUSLIION_CATHODE_DOM1D_H_

#include "m1d_porousElectrode_dom1D.h"

namespace Cantera
{
class Electrode;
class Transport;
}

namespace m1d
{
class LocalNodeIndices;


// --------------------------------------------------------------------------------------------------

//! Intermediate bookkeeping information for nodes based on loops over cells
class NodeTmps
{
public:
 
    //!  Pointer to the NodalVars struct 
    NodalVars *nv;

    //!  Offset of the nodal variables from the start of the solution vector
    size_t index_EqnStart;

    //! Offset of axial
    size_t Offset_Velocity_Axial;

    //!  Offset of variables wrt the start of the nodal solution vector.
    size_t Offset_Voltage;

    //!  Offset of variables wrt the start of the nodal solution vector.
    size_t Offset_MoleFraction_Species;

    //!  Offset of Residual for the current conservation equation
    size_t RO_Current_Conservation;
    size_t RO_Electrolyte_Continuity;
    size_t RO_Species_Eqn_Offset;
    size_t RO_MFSum_offset;
    size_t RO_ChargeBal_offset;


};

// --------------------------------------------------------------------------------------------------

//! Intermediate bookkeeping information for loops over cells
/*!
 *    If we are on a left boundary, there will be no Left Node. Instead, the left and center node
 *    are really collocated.  In that case NodeLeft_ will be a duplicate of NodeCenter_.
 *    And, nodeLeft member value will be set to zero.
 *
 *    An analogous treatment of right boundaries where there is no right node is also done. 
 */
class cellTmps
{
public:
    NodalVars* nvLeft_;
    NodalVars* nvCent_;
    NodalVars* nvRight_;

    NodeTmps NodeTmpsLeft_;
    NodeTmps NodeTmpsCenter_;
    NodeTmps NodeTmpsRight_;

    double xdelL_; // Distance from the center node to the left node
    double xdelR_; // Distance from the center node to the right node
    double xdelCell_; // cell width - right boundary minus the left boundary.
    double xCellBoundaryL_; //cell boundary left
    double xCellBoundaryR_; //cell boundary right
};

//=====================================================================================================================
//! This is derived class provides the function
//! evaluation for a porous electrolyte cathode .
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
class porousLiIon_Cathode_dom1D : public porousElectrode_dom1D
{

public:

    //! Constructor
    /*!
     * @param bdd   Contains the bulk domain description.
     */
    porousLiIon_Cathode_dom1D(m1d::BulkDomainDescription& bdd);

    //! Copy constructor
    /*!
     * @param r      Object to be copied into the current object
     */
    porousLiIon_Cathode_dom1D(const porousLiIon_Cathode_dom1D& r);

    //! Destructor
    virtual
    ~porousLiIon_Cathode_dom1D();

    //! Assignment operator
    /*!
     * @param r      Object to be copied into the current object
     * @return       Returns a changeable reference to the current object
     */
    porousLiIon_Cathode_dom1D&
    operator=(const porousLiIon_Cathode_dom1D& r);

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


    //! Set up tmps for quick calculation of residuals
    void 
    residSetupTmps();

    //! Revert the domain object's conditions to the conditions at the start of the global time step
    /*!
     *  If there was a solution in t_final, this is wiped out and replaced with the solution at t_init_init.
     *
     *  Virtual from m1d_domain.h
     */
    virtual void
    revertToInitialGlobalTime();


    //! Basic function to calculate the residual for the domain.
    /*!
     *  All residual terms are written with the following sign convention
     *  based on keeping the time derivative term positive.
     *
     *       res = dcdt - dc2 /dx2 - src = 0
     *
     * @param res  Output vector containing the residual
     * @param doTimeDependentResid  boolean indicating whether the time
     *                         dependent residual is requested
     * @param soln_ptr     solution vector at which the residual should be
     *                     evaluated
     * @param solnDot_ptr  solution dot vector at which the residual should
     *                     be evaluated.
     *  @param t           time
     *  @param rdelta_t    inverse of delta_t
     */
    virtual void
    residEval(Epetra_Vector& res,
              const bool doTimeDependentResid,
              const Epetra_Vector* soln_ptr,
              const Epetra_Vector* solnDot_ptr,
              const Epetra_Vector* solnOld_ptr,
              const double t,
              const double rdelta_t,
              const ResidEval_Type_Enum residType = Base_ResidEval,
              const Solve_Type_Enum solveType = TimeDependentAccurate_Solve);

  virtual void
  eval_PostSoln(
            const bool doTimeDependentResid,
            const Epetra_Vector *soln_ptr,
            const Epetra_Vector *solnDot_ptr,
            const Epetra_Vector *solnOld_ptr,
            const double t,
            const double rdelta_t);

    //! Utility function to calculate quantities before the main residual routine.
    /*!
     *  This is used for a loop over nodes. All calculated quantities must be internally storred.
     *
     *  Currently this is called during the residual evalultion of the problem.
     *
     * @param res  Output vector containing the residual
     * @param doTimeDependentResid  boolean indicating whether the time
     *                         dependent residual is requested
     * @param soln_ptr     solution vector at which the residual should be
     *                     evaluated
     * @param solnDot_ptr  solution dot vector at which the residual should
     *                     be evaluated.
     * @param solnOld_ptr  Pointer to the solution vector at the old time step
     *  @param t           time
     *  @param rdelta_t    inverse of delta_t
     */
    virtual void
    residEval_PreCalc(const bool doTimeDependentResid,
                      const Epetra_Vector* soln_ptr,
                      const Epetra_Vector* solnDot_ptr,
                      const Epetra_Vector* solnOld_ptr,
                      const double t,
                      const double rdelta_t,
                      const ResidEval_Type_Enum residType,
                      const Solve_Type_Enum solveType = TimeDependentAccurate_Solve);

    //!  Calculate the electrode reaction rates and store it in internal variables
    /*!
     *  @return Returns the maximum normal electrode subintegration steps
     */
    int
    calcElectrode();

    //!  Setup shop at a particular nodal point in the domain, calculating intermediate quantites
    //!  and updating Cantera's objects
    /*!
     *  All member data with the suffix, _Curr_, are updated by this function.
     *
     * @param solnElectrolyte_Curr  Current value of the solution vector
     */
    void
    SetupThermoShop1(const NodalVars* const nv, const doublereal* const solnElectrolyte_Curr);

    //!  Setup shop at a particular point in the domain, calculating intermediate quantites
    //!  and updating Cantera's objects
    /*!
     *  This routine will set up shop at a point intermediate to a left and a right point
     *  All member data with the suffix, _Curr_, are updated by this function.
     *
     * @param solnElectrolyte_CurrL  Current value of the solution vector at the left side
     * @param solnElectrolyte_CurrR  Current value of the solution vector at the right side
     * @param type                  Type of call
     *                              0 - at the left cell boundary
     *                              1 - at the right cell boundary
     */
    void
    SetupThermoShop2(const NodalVars* const nvL, const doublereal* const solnElectrolyte_CurrL,
                     const NodalVars* const nvR, const doublereal* const solnElectrolyte_CurrR,
                     int type);

    //! Calculate gradients and fluxes at the current point
    /*!
     * @param xdel            size of the cell
     * @param type            Type of call
     *                              0 - at the left cell boundary
     *                              1 - at the right cell boundary
     */
    void
    SetupTranShop(const double xdel, const int type);

    //! Function updates the ThermoPhase object for the electrolyte
    //! given the solution vector
    /*!
     *
     * @param solnElectrolyte
     */
    void
    updateElectrolyte(const NodalVars* const nv, const doublereal* const solnElectrolyte);

    //! Functions updates the Electrode object from the current values that are stored within the object
    void
    updateElectrode();

    //! Retrieves the voltages from the solution vector and puts them into local storage
    /*!
     * @param solnElectrolyte start of the solution vector at the current node
     */
    void
    getVoltages(const NodalVars* const nv, const double* const solnElectrolyte);

    void
    getMFElectrolyte_soln(const NodalVars* const nv, const double* const solnElectrolyte);

    //! Base Class for reading the solution from the saved file
    /*!
     *  This class assumes that the XML_Node is <domain> in the example below.
     *
     *  <simulation id="0">
     *    <time type="time" units="s" vtype="float"> 0.000000000000000E+00 </time>
     *    <delta_t type="time" units="s" vtype="float"> 1.000000000000000E-08 </delta_t>
     *    <StepNumber type="time" vtype="integer"> 0 </StepNumber>
     *    <domain id="BulkDomain1D_0" numVariables="6" points="10" type="bulk">
     *      <grid_data>
     *        <floatArray size="10" title="X0" type="length" units="m">
     *          0.000000000000000E+00,   8.748888888888889E-05,   1.749777777777778E-04,
     *          2.624666666666667E-04,   3.499555555555555E-04,   4.374444444444444E-04,
     *          5.249333333333334E-04,   6.124222222222222E-04,   6.999111111111111E-04,
     *          7.873999999999999E-04
     *        </floatArray>
     *     </domain>
     *  </simulation>
     *
     * @param domainNode           Reference to the XML_Node, named domain, to read the solution from
     * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
     * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
     */
    virtual void
    readDomain(const Cantera::XML_Node& domainNode,
               Epetra_Vector* const soln_GlAll_ptr,
               Epetra_Vector* const solnDot_GlAll_ptr);

    //! Base class for saving the solution on the domain in an xml node.
    /*!
     *
     * @param oNode                Reference to the XML_Node
     * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
     * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
     * @param t                    time
     *
     * @param duplicateOnAllProcs  If this is true, all processors will include
     *                             the same XML_Node information as proc 0. If
     *                             false, the xml_node info will only exist on proc 0.
     */
    virtual void
    saveDomain(Cantera::XML_Node& oNode,
               const Epetra_Vector* soln_GlAll_ptr,
               const Epetra_Vector* solnDot_GlAll_ptr,
               const double t,
               bool duplicateOnAllProcs = false);

    // Method for writing the header for the surface domain to a tecplot file.
    /*
     * Only proc0 will write tecplot files.
     */
    void writeSolutionTecplotHeader();

    // Method for writing the solution on the surface domain to a tecplot file.
    /*
     * Only proc0 will write tecplot files.
     *
     * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
     * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
     * @param t                    time
     *
     */
    void writeSolutionTecplot(const Epetra_Vector* soln_GlAll_ptr,
                              const Epetra_Vector* solnDot_GlAll_ptr,
                              const double t);

    //! Base class for writing the solution on the domain to a logfile.
    /*!
     *
     * @param soln_GlALL_ptr       Pointer to the Global-All solution vector
     * @param solnDot_GlALL_ptr    Pointer to the Global-All solution dot vector
     * @param soln_ptr             Pointer to the solution vector
     * @param solnDot_ptr          Pointer to the time derivative of the solution vector
     * @param solnOld_ptr          Pointer to the solution vector at the old time step
     * @param residInternal _ptr   Pointer to the current value of the residual just calculated
     *                             by a special call to the residEval()
     * @param t                    time
     * @param rdelta_t             The inverse of the value of delta_t
     * @param indentSpaces         Indentation that all output should have as a starter
     * @param duplicateOnAllProcs  If this is true, all processors will include
     *                             the same log information as proc 0. If
     *                             false, the loginfo will only exist on proc 0.
     */
    virtual void
    showSolution(const Epetra_Vector* soln_GlAll_ptr,
                 const Epetra_Vector* solnDot_GlAll_ptr,
                 const Epetra_Vector* soln_ptr,
                 const Epetra_Vector* solnDot_ptr,
                 const Epetra_Vector* solnOld_ptr,
                 const Epetra_Vector_Owned* residInternal_ptr,
                 const double t,
                 const double rdelta_t,
                 int indentSpaces,
                 bool duplicateOnAllProcs = false);

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
    initialConditions(const bool doTimeDependentResid,
                      Epetra_Vector* soln,
                      Epetra_Vector* solnDot,
                      const double t,
                      const double delta_t);

    //!  Fill the vector atolVector with the values from the DomainDescription for abs tol
    /*!
     * @param atolDefault             Default atol value
     * @param soln                    Solution vector. This is a constant
     *                                the residual calculation.
     * @param atolVector              (OUTPUT) Reference for the atol vector to fill up
     */
    virtual void setAtolVector(double atolDefault, const Epetra_Vector_Ghosted& soln,
                               Epetra_Vector_Ghosted& atolVector,
                               const Epetra_Vector_Ghosted* const atolV = 0);

    //!  Fill the vector atolVector with the values from the DomainDescription for abs tol
    /*!
     * @param atolDefault             Default atol value
     * @param soln                    Solution vector. This is a constant
     *                                the residual calculation.
     * @param atolVector              (OUTPUT) Reference for the atol vector to fill up
     */
    virtual void setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted& soln,
                                       const Epetra_Vector_Ghosted& solnDot,
                                       Epetra_Vector_Ghosted& atolVector_DAEInit,
                                       const Epetra_Vector_Ghosted* const atolV = 0);

    //! Evaluates the atol vector used in the delta damping process.
    /*!
     *   @param relcoeff     Relative constant to multiply all terms by
     *   @param soln         current solution vector.
     *   @param atolDeltaDamping      If non-zero, this copies the vector into the object as input
     *                      The default is zero.
     */
    virtual void
    setAtolDeltaDamping(double atolDefault, double relcoeff,
                        const Epetra_Vector_Ghosted& soln,
                        Epetra_Vector_Ghosted& atolDeltaDamping,
                        const Epetra_Vector_Ghosted* const atolV = 0);


//! Evaluates the atol vector used in the delta damping process for the DAE problem
    /*!
     *   @param relcoeff     Relative constant to multiply all terms by
     *   @param soln         current solution vector.
     *   @param solnDot      Current solutionDot vector.
     *   @param atolDeltaDamping       If non-zero, this copies the vector into the object as input
     *                       The default is zero.
     */
    virtual void
    setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff,
                                const Epetra_Vector_Ghosted& soln,
                                const Epetra_Vector_Ghosted& solnDot,
                                Epetra_Vector_Ghosted& atolDeltaDamping,
                                const Epetra_Vector_Ghosted* const atolV = 0);

    //! Set up the vector of electrode cells within the object
    void instantiateElectrodeCells();

    /**
     * Method to check for precipitation of the salts.
     * Returns index of offending cation or -1 if no precipitation
     */
    int checkPrecipitation();

    // -----------------------------------------------------------------------------------------------
    //                                 DATA
    // -----------------------------------------------------------------------------------------------

protected:

    //! Pointer to the thermo object for the molten salt
    /*!
     *   We do not own this object
     */
    Cantera::ThermoPhase* ionicLiquid_;

    //! Pointer to the transport object for the molten salt
    /*!
     * We do not own this object
     */
    Cantera::Transport* trans_;

    //! Number of phases solved
    int nph_;

    //! Number of species solved for within the domain
    int nsp_;


    //! Total concentration of the electolyte
    /*!
     *  Units = kg / m3
     */
    double concTot_cent_;

    //! Total concentration of the electolyte at the previous time step
    /*!
     *  Units = kg / m3
     */
    double concTot_cent_old_;


    //Need to define a void fraction variable here.
    //The void fraction does not change with time
    //for the electrolyte, but it will for the electrode.
    //Since I'm not sure how I want to define this porosity
    //at this instant, I will set it to a constant for now.

    // ------------------------------------------------------------------------

    //! Surface Area of the electrolyte - electrode interface within the cell
    /*!
     *  This is now a constant.
     *  Length is number of cells on the processor.
     *  units = m2 / m2. (surface area per cross-sectional area)
     */
    std::vector<double> surfaceArea_Cell_;

    //! Electrode Current per surface area of the electrode
    /*!
     *  The surface area is specifically defined as the external surface of the electrode
     *  Length is number of cells on the processor.
     *  units = amps / m2
     */
    std::vector<double> icurrInterfacePerSurfaceArea_Cell_;

    //! delta X for the current cell
    /*!
     *  Length is number of cells on the processor
     *  units = m
     */
    std::vector<double> xdelCell_Cell_;

    //! Total concentration of the electrolyte at cell centers
    /*!
     *  Length is number of cells on the processor.
     */
    std::vector<double> concTot_Cell_;

    //!Total concentration of the electrolyte at cell centers
    /*!
     *  Length is number of cells on the processor.
     */
    std::vector<double> concTot_Cell_old_;

    //! Electrode Cross-Sectional Area
    /*!
     *  The electrode objects are extrinsic objects that have a cross-sectional area associated
     *  with them. This is the value of that area. The rest of the code is based on a per m2 basis.
     *  Therefore, we need to
     */
    double electrodeCrossSectionalArea_;

    //! Electrode Cell data for the anode cells
    /*!
     *  Length is the number of owned Cells on the processor
     *   -> not sure what owned here means. This will have to be figured out in the future.
     *   Note, in contrast to the other quantities in the calculation, the electrode object is extrinsic
     *   it has a volume associated with it.  The volume is equal to the control volume thickness
     *   multiplied by the cross-sectional area.
     */
    std::vector<Cantera::Electrode*> Electrode_Cell_;

    //!  Capacity discharged by the particular electrode cell
    /*!
     *   Units:  amps * sec / m2  = coulumbs / m2
     *
     *   We calculate this quantity by taking the capacityDischarged() from the electrode object
     *   and dividing by the cross sectional area of the extrinsic electrode object
     */
    std::vector<double> capacityDischarged_Cell_;

    //!  Depth of Discharge of this particular electrode cell
    /*!
     *   Units:  amps * sec / m2  = coulumbs / m2
     *
     *   We calculate this quantity by taking the depthOfDischarge() from the electrode object
     *   and dividing by the cross sectional area of the extrinsic electrode object
     */
    std::vector<double> depthOfDischarge_Cell_;

    //!  Capacity left in this particular electrode cell
    /*!
     *   Units:  amps * sec / m2  = coulumbs / m2
     *
     *   We calculate this quantity by taking the capacityLeft() from the electrode object
     *   and dividing by the cross sectional area of the extrinsic electrode object
     */
    std::vector<double> capacityLeft_Cell_;

    //!  Capacity of this particular electrode cell if it were at zero depth of discharge
    /*!
     *   Units:  amps * sec / m2  = coulumbs / m2
     *
     *   We calculate this quantity by taking the capacityZeroDoD() from the electrode object
     *   and dividing by the cross sectional area of the extrinsic electrode object
     */
    std::vector<double> capacityZeroDoD_Cell_;

    // ------------------------------------------------------------------------
    //!  Cell storage -> storage of cell related quantities

    //! Cell index number
    int cIndex_cc_;

    //! Axial velocity - left cell boundary
    double Fleft_cc_;
    //! Axial Velocity - right cell boundary
    double Fright_cc_;

    //! Electrostatic potential - Left cell
    double Vleft_cc_;
    //! Electrostatic potential - center cell
    double Vcent_cc_;
    //! Electrostatic potential - right cell
    double Vright_cc_;

    //! Electrostatic potential in the electrode - Left cell
    double VElectrodeLeft_cc_;
    //! Electrostatic potential in the electrode  - center cell
    double VElectrodeCent_cc_;
    //! Electrostatic potential  in the electrode - right cell
    double VElectrodeRight_cc_;

    double t_final_;
    double t_init_;

    //! Mole fraction of electrolyte species in the left cell
    /*!
     * Length = number of electrolyte species = 3
     */
    std::vector<double> Xleft_cc_;

    //! Mole fraction of electrolyte species in the center cell (i.e., the current cell)
    /*!
     * Length = number of electrolyte species = 3
     */
    std::vector<double> Xcent_cc_;

    //! Mole fraction of electrolyte species in the right cell
    /*!
     * Length = number of electrolyte species = 3
     */
    std::vector<double> Xright_cc_;

    //! Charge of the species in the electrolyte
    /*!
     * Length = number of electrolyte species = 3
     */
    std::vector<double> spCharge_;

    // -----------------------------------------------------------------------
    //!  Current Thermo value of quantities at the current point


    std::vector<double> mfElectrolyte_Soln_Curr_;

    std::vector<double> mfElectrolyte_Thermo_Curr_;

    Cantera::Array2D mfElectrolyte_Soln_Cell_old_;

    //! Current value of the cathode voltage
    double phiElectrode_Curr_;

    //! Electrical conductivity of the electrode
    /*!
     *   units are S m-1
     *   The default is 1 E6
     */
    double conductivityElectrode_;

    // --------------------------------------------------------------------------
    //!  Current transport values of quantities at the current point

    //! Gradient of the temperature
    double gradT_trCurr_;

    //! Gradient of the Electric Potential
    double gradV_trCurr_;

    //! Gradient of the potential in the electrode phase
    double gradVElectrode_trCurr_;

    // Gradient of the Mole fraction
    std::vector<double> gradX_trCurr_;

    std::vector<double> Vdiff_trCurr_;

    //! Diffusive flux of species in the electrolyte
    /*!
     *
     */
    std::vector<double> jFlux_trCurr_;

    //! current flow in the electrode due to conduction
    double icurrElectrode_trCurr_;

    //! Number of species in the electrode object
    int nSpeciesElectrode_;

    //! Number of interphase surfaces or plateaus in the electrode object
    int nSurfsElectrode_;

    //! Electrode production rate delta for all species in the Electrode object at all cells
    /*!
     *   Units are kmols
     *   Length is the number of species in the electrode object multiplied by the number of cells
     *   Outer loop is over cells:
     *     electrodeSpeciesMoleDelta_Cell_[iCell * nSpeciesElectrode_ + k]
     */
    std::vector<double> electrodeSpeciesMoleDelta_Cell_;

    //! Value of the interface current from the electrode in the cell
    /*!
     *   Units are amps / m2 and the area is the cross sectional area of the electrode
     */
    std::vector<double> icurrInterface_Cell_;

    //! Phase mole transfers from the electrode reactions
    /*!
     *  Vector of mole transfers for each phase in the electrode
     *  Units = kmol
     */
    std::vector<double> phaseMoleTransfer_;

    //! soln Phase mole fluxes from the electrode reactions
    /*!
     *  Units = kmol /m2 /sec
     */
    std::vector<double> solnMoleFluxInterface_Cell_;

    //! Electrode Current at the cell boundaries: left boundary
    /*!
     *  Length is the number of cells
     */
    std::vector<double> icurrElectrode_CBL_;

    //! Electrode Current at the cell boundaries: right boundary
    /*!
     *  Length is the number of cells
     */
    std::vector<double> icurrElectrode_CBR_;

    //! Electrolyte Current at the cell boundaries - left
    std::vector<double> icurrElectrolyte_CBL_;

    //! Electrolyte Current at the cell boundaries - right
    std::vector<double> icurrElectrolyte_CBR_;

    //! Electrostatic potential difference between electrolyte and metal phases
    std::vector<double> deltaV_Cell_;

    //! open circuit potential for each surface (plateau)
    std::vector<double> Ess_Cell_;
    //!overpotential for each surface (plateau)
    std::vector<double> overpotential_Cell_;
    std::vector<double> icurrRxn_Cell_;
    std::vector<double> LiFlux_Cell_;

    // --------------------------------------------------------------------------

    //! species index of the solvent species in the electrolyte phase - 0
    /*!
     *   This will get the sum of mole fractions equals 1 equation applied to it.
     */
    int iECDMC_;

    //! species index of the Li+ species in the electrolyte phase - 1
    int iLip_;

    //! species index of the Li+ species in the electrolyte phase - 2
    /*!
     *   This will get the charge balance equation applied to it
     */
    int iPF6m_;

    std::vector<double> solnTemp;

    //! Vector of temporary indexing quantities for each cell
    /*!
     * These are calculated once at the start of the program
     */
    std::vector<cellTmps> cellTmpsVect_Cell_;

    //! Velocity basis of the transport equations
    Cantera::VelocityBasis ivb_;

private:
    void
    err(const char* msg);

public:

};
//=====================================================================================================================
}
//=====================================================================================================================
#endif // M1D_POROUSLIKCL_FES2CATHODE_DOM1D_H_
