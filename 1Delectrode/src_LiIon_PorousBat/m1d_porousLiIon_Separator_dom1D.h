/*
 * m1d_porousLiKCl_dom1D.h
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 */

#ifndef M1D_POROUSLIION_SEPARATOR_DOM1D_H_
#define M1D_POROUSLIION_SEPARATOR_DOM1D_H_

#include "m1d_porousFlow_dom1D.h"

#include <cantera/transport.h>      // transport properties

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

// --------------------------------------------------------------------------------------------------

//!  This is derived class  provides the function evaluation for a porous electrolyte bulk domain.
/*!
 *  The porous electrolyte domain is characterized by a
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
class porousLiIon_Separator_dom1D : public porousFlow_dom1D
{

public:

    //! Constructor
    /*!
     * @param bdd   Contains the bulk domain description.
     */
    porousLiIon_Separator_dom1D(m1d::BulkDomainDescription& bdd);

    //! Copy constructor
    /*!
     * @param r      Object to be copied into the current object
     */
    porousLiIon_Separator_dom1D(const porousLiIon_Separator_dom1D& r);

    //! Destructor
    virtual
    ~porousLiIon_Separator_dom1D();

    //! Assignment operator
    /*!
     * @param r      Object to be copied into the current object
     * @return       Returns a changeable reference to the current object
     */
    porousLiIon_Separator_dom1D&
    operator=(const  porousLiIon_Separator_dom1D& r);

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

    void
    getVoltages(const NodalVars* const nv, const double* const solnElectrolyte);

    void
    getMFElectrolyte_soln(const NodalVars* const nv, const double* const solnElectrolyte);

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

    //! number of species solved
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

    //! Local value of the temperature
    //double temp_Curr_;

    //! Local value of the pressure
    // double pres_Curr_;

    //!  Current value of the voltage
    //double phiElectrolyte_Curr_;

    //!  Total concentration
    double concTot_Curr_;

    //! Current porosity
    double porosity_Curr_;

    // --------------------------------------------------------------------------
    //!  Current transport values of quantities at the current point

    //! Gradient of the temperature
    double gradT_trCurr_;

    // Gradient of the Electric Potential
    double gradV_trCurr_;

    // Gradient of the Mole fraction
    std::vector<double> gradX_trCurr_;

    std::vector<double> Vdiff_trCurr_;

    //! diffusive flux
    /*!
     *
     */
    std::vector<double> jFlux_trCurr_;

    // --------------------------------------------------------------------------

    //! Electrolyte Current at the cell boundaries - left
    std::vector<double> icurrElectrolyte_CBL_;

    //! Electrolyte Current at the cell boundaries - right
    std::vector<double> icurrElectrolyte_CBR_;

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

}

#endif // M1D_POROUSLIKCL_DOM1D_H_
