/**
 *  @file m1d_BatteryResidEval.h
 *  File that contains the description of a single m1d problem that
 *  will be solved.
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef M1D_BATTERYRESIDEVAL_H
#define M1D_BATTERYRESIDEVAL_H

#include "m1d_ProblemResidEval.h"
#include "m1d_BoundaryCondition.h"
#include "m1d_SurDomain1D.h"

#include "cantera/base/Array.h"
#include <vector>

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d {

//==================================================================================================================================
//! Class that specializes the global heat balance for batteries
/*!
 *  This is filled up on a domain basis.
 *  It inherits from globalHeatBalVals to add terms due to batteries.
 *
 *  This heat balance class is used for both volume domains and surface domains.
 */
class globalHeatBalValsBat : public globalHeatBalVals
{
public:
    //! Default constructor
    globalHeatBalValsBat() :
    globalHeatBalVals(),
      totalEnthalpyInit(0.0),
      totalEnthalpyFinal(0.0),
      currentRight(0.0),
      currentLeft(0.0),
      JHelecRight(0.0),
      JHelecLeft(0.0),
      phiSolid(0.0),
      phiLyte(0.0),
      jouleHeat_lyte(0.0),
      enthalpyIVfluxRight(0.0),
      enthalpyIVfluxLeft(0.0),
      sourceTermExtra(0.0)
      {
	  
      }
    
    //! Destructor
    virtual ~globalHeatBalValsBat()
    {
    }


    //! Zero all of the values within the structure
    virtual void zero()
    {
	globalHeatBalVals::zero();
	totalEnthalpyInit = 0.0;
	totalEnthalpyFinal = 0.0;
	currentRight = 0.0;
	currentLeft = 0.0;
	JHelecRight = 0.0;
	JHelecLeft = 0.0;
	phiSolid = 0.0;
	phiLyte = 0.0;
        jouleHeat_lyte = 0.0;
        enthalpyIVfluxRight = 0.0;
        enthalpyIVfluxLeft = 0.0;
        sourceTermExtra = 0.0;
	size_t nsp = species_Lyte_New_Total.size();
	for (size_t k = 0; k < nsp; ++k) {
	    species_Lyte_New_Total[k] = 0.0;
	    species_Lyte_Old_Total[k] = 0.0;
	    species_Lyte_Src_Total[k] = 0.0;
	    species_convRight[k] = 0.0;
	    species_convLeft[k] = 0.0;
	    species_jFluxRight[k] = 0.0;
	    species_jFluxLeft[k] = 0.0;
	}
        elem_Lyte_New.resize(10, 0.0);
        elem_Solid_New.resize(10, 0.0);
        elem_Lyte_Old.resize(10, 0.0);
        elem_Solid_Old.resize(10, 0.0);
	for (size_t k = 0; k < 10; ++k) {
            elem_Lyte_New[k] = 0.0;
            elem_Solid_New[k] = 0.0;
            elem_Lyte_Old[k] = 0.0;
            elem_Solid_Old[k] = 0.0;
	}

    }

    //! Resize all of the variables associated with the Electrolyte according to the argument value
    /*!
     *  @param[in]           nsp                 Number of species in the electrolyte
     */
    virtual void sizeLyte(size_t nsp)
    {
	species_Lyte_New_Total.resize(nsp, 0.0);
	species_Lyte_Old_Total.resize(nsp, 0.0);
	species_Lyte_Src_Total.resize(nsp, 0.0);
	species_convRight.resize(nsp, 0.0);
	species_convLeft.resize(nsp, 0.0);
	species_jFluxRight.resize(nsp, 0.0);
	species_jFluxLeft.resize(nsp, 0.0);

        elem_Lyte_New.resize(10, 0.0);
        elem_Solid_New.resize(10, 0.0);
        elem_Lyte_Old.resize(10, 0.0);
        elem_Solid_Old.resize(10, 0.0);
    }

    // ---------------------------------------------- D A T A -------------------------------------------------------------

    //! Extrinsic total initial enthalpy - initial
    /*!
     *  Units:  Joules
     */
    double totalEnthalpyInit;

    //! Extrinsic total initial enthalpy - final
    /*!
     *  Units:  Joules
     */
    double totalEnthalpyFinal;

    //! Current leaving the domain on the right.
    /*!
     *  This is filled in only for the cathode domain ATM
     *  Units: amps / m2  (area is the cross-sectional aria of the battery)
     */
    double currentRight;

    //! Current leaving the domain on the left.
    /*!
     *  This is filled in only for the anode domain ATM
     *  Units: amps / m2  (area is the cross-sectional aria of the battery)
     */
    double currentLeft;

    //! Value of the flux of electrons multiplied by the enthalpyPhi of the electrons at the rhs of domain
    /*!
     *  EnthalpyPhi = H_elec(T,P) + z_elec F Phi  -> units:  Joule / kmol
     *  flux = - icurr / F                        -> units:  kmol / sec /m2
     *
     *  Used:   For volume domains
     *  Units:  Joule / sec / m2
     */
    double JHelecRight;

    //! Value of the flux of electrons multiplied by the enthalpyPhi of the electrons at the lhs of domain
    /*!
     *  EnthalpyPhi = H_elec(T,P) + z_elec F Phi  -> units:  Joule / kmol
     *  flux = - icurr / F                        -> units:  kmol / sec /m2
     *
     *  Used:   For volume domains
     *  Units:  Joule / sec / m2
     */
    double JHelecLeft;

    //! Value of the electric potential in the solid electrode
    /*!
     *  Used:   For surface domains
     *  Units:  volts
     */
    double phiSolid;

    //! Value of the electric potential in the electrolyte
    /*!
     *  Used:   For surface domains
     *  Units:  volts
     */
    double phiLyte;

    //! Value of the joule heat term for the domain
    //! Source of joule heating in the electrolyte during the current time step for the domain
    /*!
     *  This is the heat generated due to joule heating in the electrolyte in in the domain per time for the current time step
     *  Therefore, the source term is integrated in the axial direction and it is integrated wrt to the time interval.
     *
     *  Units:   Joules / m2
     */
    double jouleHeat_lyte;

    //! Value of heat flux carried by the current-voltage combination out the RHS of domain, i.e., the cathode current collector
    /*!
     *  This is essentially i * V
     *
     *  Units = (amps / m2 ) ( V ) = (amps / m2 ) ( Joule / coul ) = Joule / sec / m2
     */ 
    double enthalpyIVfluxRight;

    //! Value of heat flux carried by the current-voltage combination out the LHS of domain, i.e., the anode current collector
    /*!
     *  This is essentially i * V
     *
     *  Units = (amps / m2 ) ( V ) = (amps / m2 ) ( Joule / coul ) = Joule / sec / m2
     */ 
    double enthalpyIVfluxLeft;

    double sourceTermExtra;

    std::vector<double> elem_Lyte_New;
    std::vector<double> elem_Solid_New;
    std::vector<double> elem_Lyte_Old;
    std::vector<double> elem_Solid_Old;

    std::vector<double> species_Lyte_New_Total;
    std::vector<double> species_Lyte_Old_Total;
    std::vector<double> species_Lyte_Src_Total;
    std::vector<double> species_convRight;
    std::vector<double> species_convLeft;
    std::vector<double> species_jFluxRight;
    std::vector<double> species_jFluxLeft;

};
//==================================================================================================================================
//!  Residual for 1D cell battery evaluations
/*!
 *   We add in a lot of specific functionality for batteries here.
 *   We tailor the output to print out battery related information.
 */
class BatteryResidEval: public ProblemResidEval
{

public:
    //! Default constructor
    /*!
     *
     * @param atol   Absolute tolerance for calculation
     */
    BatteryResidEval(double atol = 1.0e-13);

    //! Destructor
    virtual ~BatteryResidEval();

    //! Default copy constructor
    /*!
     *
     * @param r  Object to be copied
     */
    BatteryResidEval(const BatteryResidEval &r);

    //! Assignment operator
    /*!
     *
     * @param r  Object to be copied
     * @return   Returns a copy of the current problem
     */
    BatteryResidEval& operator=(const BatteryResidEval&r);

    void residSetupTmps();

    //! Set the underlying state of the system from the solution vector
    /*!
     *   Note this is an important routine for the speed of the solution.
     *   It would be great if we could supply just exactly what is changing here.
     *   This routine is always called at the beginning of the residual evaluation process.
     *
     *   This is a natural place to put any precalculations of nodal quantities that
     *   may be needed by the residual before its calculation.
     *
     *   Also, this routine is called with rdelta_t = 0. This implies that a step isn't being taken. However, the
     *   the initial conditions must be propagated.
     *
     * @param doTimeDependentResid
     * @param soln
     * @param solnDot
     * @param t                    This is not necessarily equal to t_old + delta_t, but can be anywhere between the two values
     *                   
     * @param delta_t If zero then delta_t equals 0.
     */
    virtual void
    setStateFromSolution(const bool doTimeDependentResid, const Epetra_Vector_Ghosted *soln, const Epetra_Vector_Ghosted *solnDot,
                         const double t, const double delta_t, const double t_old);

    //! Calculate the initial conditions
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  This calls the parent class initialConditions method to loop over the volume and surface domains.
     *  Then the method tried to better estimate the electrolyte potential.
     *
     *  @param[in]           doTimeDependentResid  Boolean indicating whether we should formulate the time dependent residual
     *  @param[in]           soln                Solution vector. This is the input to  the residual calculation.
     *  @param[in]           solnDot             Solution vector. This is the input to the residual calculation.
     *  @param[in]           t                   Time
     *  @param[in]           delta_t             delta_t for the initial time step
     */
    virtual void
    initialConditions(const bool doTimeDependentResid, Epetra_Vector_Ghosted *soln, Epetra_Vector_Ghosted *solnDot, double &t,
                      double &delta_t, double& delta_t_np1) override;

    //! Make an attempt to improve the initial guess for the electrolyte voltage  based on the boundary conditions
    /*!
     *  @param[in, out]      soln                Pointer to the solution vector. This is the input to the residual calculation.
     */
    void improveInitialConditions(Epetra_Vector_Ghosted * const soln);

    //! Calculate a residual vector, given a solution vector at a new time and delta t 
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  This is the main routine to calculate a solution vector.
     *  The basic algorithm is to loop over the volume domains. Then, we loop over the surface domains.
     *
     *  @param[out]          res                 Residual output
     *  @param[in]           doTimeDependentResid   Boolean indicating whether the time dependent residual is requested
     *  @param[in]           soln                Pointer to the solution vector. This is the input to the residual calculation.
     *  @param[in]           solnDot             Pointer to the solution Dot vector. This is the input to the residual calculation.
     *  @param[in]           t                   current time
     *  @param[in]           rdelta_t            delta t inverse
     *  @param[in]           residType           Residual type
     *                                                Defaults to a a base residual evaluation
     *  @param[in]           solveType           Solve type . 
     *                                                Defaults to a time dependent solve
     */
    virtual void
    residEval(Epetra_Vector* const &  res,
	      const bool doTimeDependentResid,
	      const Epetra_Vector_Ghosted *soln,
	      const Epetra_Vector_Ghosted *solnDot,
	      const double t,
	      const double rdelta_t,
	      const Zuzax::ResidEval_Type residType = Zuzax::ResidEval_Type::Base_ResidEval,
	      const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve);

    //! Write the solution to either the screen or to a log file
    /*!
     *  This is a general output utility to Cantera's logfile.
     *  It's not hooked into the IO algorithm at all. It should be conditionally called depending on the whims of the user.
     *
     *  @param[in]           ievent              Type of the event. The following form is used:
     *                                             - 0 Initial conditions
     *                                             - 1 Completion of a successful intermediate step.
     *                                             - 2 Final successful conditions.
     *                                             - 3 Intermediate nonlinear step
     *                                             - -1 unsuccessful step
     *  @param[in]           doTimeDependentResid   Do the time dependent residual calculation
     *  @param[in]           t                   Current time
     *  @param[in]           delta_t             delta t
     *  @param[in]           y_n                 Current value of the solution vector
     *  @param[in]           ydot_n              Current value of the derivative of the solution vector
     *  @param[in]           solveType           Solve type
     *  @param[in]           delta_t_np1         delta_t for the next step
     */
    virtual void
    showProblemSolution(const int ievent,
			bool doTimeDependentResid,
			const double t,
			const double delta_t,
			const Epetra_Vector_Owned &y_n,
			const Epetra_Vector_Owned * const ydot_n,
			const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve,
			const double delta_t_np1 = 0.0) override;
    

    //! Write out to a file or to standard output the current solution
    /*!
     *   These functions are affected by the print controls of the nonlinear solver
     *   and the time stepper.
     *
     *      ievent is a description of the event that caused this
     *      function to be called.
     *
     *      @param ievent  Event that's causing this routine to be called.
     *                     =  0 Initial conditions for a calculation
     *                     =  1 Completion of a successful intermediate time step.
     *                     =  2 Completion of a successful Final time or final calculation.
     *                     =  3 Completion of a successful Intermediate nonlinear step
     *                     = -1 unsuccessful time step that converged, but failed otherwise
     *                     = -2 unsuccessful nonlinear step.
     *      @param doTimeDependentResid  true if solving a time dependent problem
     *      @param time_current      Current time
     *      @param delta_t_n         Current value of delta_t
     *      @param istep             Current step number
     *      @param soln_n               Current value of the solution vector
     *      @param solnDot_n_ptr            Current value of the time deriv of the solution vector
     *      @param delta_t_np1       Suggested next delta t value (defaults to 0.0 if there isn't a
     *                               good guess for the next delta_t).
     */
    virtual void
    writeSolution(const int ievent, const bool doTimeDependentResid, const double time_current, const double delta_t_n,
                  const int istep, const Epetra_Vector_Ghosted &soln_n, const Epetra_Vector_Ghosted * const solnDot_n_ptr,
                  const Zuzax::Solve_Type solveType = Zuzax::Solve_Type::TimeDependentAccurate_Solve, 
	          const double delta_t_np1 = 0.0) override;

    //! This function may be used to create output at various points in the execution of an application.
    /*!
     *  These functions are not affected by the print controls of the nonlinear solver  and the time stepper.
     *
     *  ievent is a description of the event that caused this function to be called.
     *
     *  @param[in]           ievent              Event that's causing this routine to be called.
     *                                              =  0 Initial conditions for a calculation
     *                                              =  1 Completion of a successful intermediate time step.
     *                                              =  2 Completion of a successful Final time or final calculation.
     *                                              =  3 Completion of a successful Intermediate nonlinear step
     *                                              =  4 Write out current and voltage to timeDep_IV.dat
     *                                              = -1 unsuccessful time step that converged, but failed otherwise
     *                                              = -2 unsuccessful nonlinear step.
     *
     *  @param[in]           time_current        Current time
     *  @param[in]           delta_t_n           Current value of delta_t
     *  @param[in]           istep               Current step number
     *  @param[in]           y_n                 Current value of the solution vector
     *  @param[in]           ydot_n_ptr          Current value of the time deriv of the solution vector
     */
    virtual void
    user_out(const int ievent, const double time_current, const double delta_t_n, const int istep, const Epetra_Vector_Ghosted &y_n,
             const Epetra_Vector_Ghosted* const ydot_n_ptr);

    //! Write a tecplot file consisting of the voltage and current as a function of time
    /*!
     *      @param ievent  Event that's causing this routine to be called.
     *                     =  0 Initial conditions for a calculation
     *                     =  1 Completion of a successful intermediate time step.
     *                     =  2 Completion of a successful Final time or final calculation.
     *                     =  3 Completion of a successful Intermediate nonlinear step
     *                     =  4 Write out current and voltage to timeDep_IV.dat
     *                     = -1 unsuccessful time step that converged, but failed otherwise
     *                     = -2 unsuccessful nonlinear step.
     *      @param doTimeDependentResid  true if solving a time dependent problem
     *      @param time_current      Current time
     *      @param delta_t_n         Current value of delta_t
     *      @param istep             Current step number
     *      @param y_n               Current value of the solution vector
     *      @param ydot_n_ptr        Current value of the time deriv of the solution vector
     */
    virtual void 
    write_IV(const int ievent, const bool doTimeDependentResid, const double time_current, 
             const double delta_t_n, int istep, const Epetra_Vector_Ghosted &y_n, 
             const Epetra_Vector_Ghosted * const ydot_n_ptr);

    //! Write a global tecplot file that includes variables which span all of the bulk domains
    /*!
     *   (Note, this may be replaced with a tecplot substitute that we are working on
     *
     *   @param ievent
     *
     */
    virtual void
    writeGlobalTecplot(const int ievent,
		       const bool doTimeDependentResid,
		       const double time_current,
		       const double delta_t_n,
		       int istep,
		       const Epetra_Vector_Ghosted &y_n,
		       const Epetra_Vector_Ghosted * const ydot_n_ptr,
		       const Zuzax::Solve_Type solveType, 
		       const double delta_t_np1);

    virtual void
    writeGlobalTecplotHeader(const int ievent,
			     const bool doTimeDependentResid,
			     const double time_current,
			     const double delta_t_n,
			     int istep,
			     const Epetra_Vector_Ghosted &y_n,
			     const Epetra_Vector_Ghosted * const ydot_n_ptr,
			     const Zuzax::Solve_Type solveType, 
			     const double delta_t_np1);
    
    //! Evaluate a supplemental set of equations that are not part of the solution vector, but are considered
    //! to be time dependent, at the end of every successful time step
    /*!
     *  (virtual from ProblemResidEval) 
     *
     *  Equations in this system are evaluated using the time discretization scheme of the nonlinear solver.
     *  It can be used to track supplemental quantities in the calculation, especially if they need to be
     *  integrated in time.
     *
     *  An example of this may be total flux quantites that are dumped into a phase.
     *
     *  This routine is called at the beginning of the time stepping, in order to set up any initializations,
     *  and it is called after every successful time step, once and only once.
     *
     *  This is used to calculate the heat source terms in a battery.
     *
     *  @param[in]           ifunc               Identity of the type of call to the routine:
     *                                            - 0 Initial call to the function, done whenever the time stepper is entered
     *                                            - 1 Called after every successful time step.
     *  @param[in]           t                   Current time
     *  @param[in]           deltaT              Current value of deltaT
     *  @param[in]           soln                Current value of the solution vectors
     *  @param[in]           solnDot_ptr         Current value of time derivative of the solution vectors.
     */
    virtual void
    evalTimeTrackingEqns(const int ifunc, const double t, const double deltaT, const Epetra_Vector_Ghosted& soln,
			 const Epetra_Vector_Ghosted* const solnDot_ptr) override;

    //! Carries out a balance on EnthalpyPhi over the domains covered by the simulation
    /*!
     *  We define EnthalpyPhi as the enthalpy function with additions from the potential energy due to charged ions in an
     *  electric potential.
     *
     *  This prints out the total EnthalpyPhi balance over all domains. We can make this more general
     *  in the future. But for now, this works as a proof of strict conservation of EnthalpyPhi up to round-off-error.
     *
     *  @param[in]           ifunc               Identity of the type of call to the routine:
     *                                            - 0 Initial call to the function, done whenever the time stepper is entered
     *                                            - 1 Called after every successful time step.
     *  @param[in]           t                   Current time
     *  @param[in]           deltaT              Current value of deltaT
     *  @param[in]           soln                Current value of the solution vectors
     *  @param[in]           solnDot_ptr         Current value of time derivative of the solution vectors.
     */
    void doHeatAnalysis(const int ifunc, const double t, const double deltaT, const Epetra_Vector_Ghosted& soln,
	                const Epetra_Vector_Ghosted* const solnDot_ptr);

    //! Do an analysis of the polarization losses at the current global time step
    /*!
     *  (virtual from BatteryResidEval) 
     *
     *  Part of the analysis is that we have to find the global extent of reaction of the anode and cathode in order to
     *  find a global open circuit potential.
     *
     *  @param[in]           ifunc               0 Initial call to the function, done whenever the time stepper is entered
     *                                           1 Called after every successful time step.
     *  @param[in]           t                   Current time
     *  @param[in]           deltaT              Current value of deltaT
     *  @param[in]           soln                Current value of the solution vectors
     *  @param[in]           solnDot_ptr         Current value of time derivative of the solution vectors.
     */
    virtual void 
    doPolarizationAnalysis(const int ifunc, const double t, const double deltaT, const Epetra_Vector_Ghosted& soln,
                           const Epetra_Vector_Ghosted* const solnDot_ptr);

    //! Carries out a species element balance over all of the domains
    /*!
     *  This prints out the total element balance of Lithium over all domains, atm. We can make this more general
     *  in the future. But for now, this works as a proof of strict conservation up to round-off-error.
     *
     *  @param[in]           ifunc               0 Initial call to the function, done whenever the time stepper is entered
     *                                           1 Called after every successful time step.
     *  @param[in]           t                   Current time
     *  @param[in]           deltaT              Current value of deltaT
     *  @param[in]           y                   Current value of the solution vectors
     *  @param[in]           solnDot_ptr         Current value of time derivative of the solution vectors.
     */
    void doSpeciesAnalysis(const int ifunc, const double t, const double deltaT, const Epetra_Vector_Ghosted& y,
                           const Epetra_Vector_Ghosted* const solnDot_ptr);
    
    //! Set a solution parameter 
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  @param[in]           paramName           String identifying the parameter to be set
     *  @param[in]           paramVal            Single double value of the parameter to be set
     *
     *  @return                                  Returns the number of parameters set by the command, zero or a positive number.
     *                                           Returns a negative number if the parameter is unknown
     */
    virtual int
    setSolutionParam(std::string paramName, double paramVal) override;

    //! Get a solution parameter 
    /*!
     *  (virtual from ProblemResidEval)
     *
     *  @param[in]           paramName           String identifying the parameter to be retrieved
     *  @param[in]           paramVal            Vector of parameters returned
     *
     *  @return                                  Returns the number of parameters returned.
     */
    virtual int
    getSolutionParam(std::string paramName, double* const paramVal) override;

    //! Report on the boundary condition applied to the cathode voltage equation
    /*!
     *
     *  @param[in]  time        Current time for evaluating time dependent BC
     *  @param[out] BC_Type     Type of the boundary condition
     *  @param[out] value       Value of the dirichlet condition or flux - default 0.0
     *  @param[out] BC_TimeDep  BoundaryCondition Pointers for time dependent BC for BC_Tppe = 3,4,5
     *                            (default 0)
     *  @param[out] TimeDep     Function pointer to a function that returns a double given a single parameter (the time).
     *                           Defaults to a NULL pointer.
     */
    void
    reportCathodeVoltageBC(double time, int &BC_Type, double &value, BoundaryCondition * &BC_TimeDep,
                           TimeDepFunctionPtr &TimeDep) const;

    //! Change the boundary condition applied to the cathode voltage equation
    /*!
     *  @param[in]           BC_Type             Type of the boundary condition
     *  @param[in]           value               Value of the Dirichlet condition or flux - default 0.0
     *  @param[in]           BC_TimeDep          BoundaryCondition Pointers for time dependent BC for BC_Type = 3,4,5
     *                                             Defaults to a NULL pointer.
     *  @param[in]           TimeDep             Function pointer to a function that returns a double given a single parameter (the time).
     *                                             Defaults to a NULL pointer.
     */
    void
    changeCathodeVoltageBC(int BC_Type, double value, BoundaryCondition * BC_TimeDep = 0, TimeDepFunctionPtr TimeDep = 0);

    //! Report the cathode voltage
    /*!
     *  (virtual from BatteryResidEval)
     *  Note the anode voltage is set to zero. So this function will return the voltage of the battery.
     *
     *  @return                                  Returns the cathode voltage (volts)
     */
    virtual double reportCathodeVoltage() const;

    //! Report the cathode current assuming an area of the battery.
    /*!
     *  (virtual from BatteryResidEval)
     *  Note the anode current should be exactly equal to the cathode current
     *  The current will be positive for the normal operation of the battery.
     *
     *  @return                                  Returns the cathode current (amp)
     */
    virtual double reportCathodeCurrent() const;

    //! Get the max value of the sub grid time step number from the last residual calculation
    /*!
     *  @return                                  Returns the max subgrid time step number from all of the electrodes. Special steps
     *                                           aren't counted in this number.
     */
    int getMaxSubGridTimeSteps() const;

    double heatSourceLastStep() const;

    double heatSourceAccumulated() const;

    void heatSourceZeroAccumulated() const;

    //! Here we gather statistics about the current state of the battery by summing up all of the 
    //! individual cells for the anode and cathode
    /*!
     *   The following variables are calculated by this routine
     */
    void gatherCapacityStatistics();

    // ------------------------------------------------- MEMBER DATA ----------------------------------------------------------------

    //! Boolean indicating whether to calculate Heat Source Time tracking terms and output file
    int doHeatSourceTracking_;

    //! Tracking of the polarization resistance in the domain
    int doPolarizationAnalysis_;

    //! Boolean indicating whether to calculate the Electrical resistance tracking terms and output file
    int doResistanceTracking_;

    //! Type of anode
    /*!
     *  Determines whether the anode is distributed or is part of a SurDomain1D.
     *  Right now, this defaults to 0, but can be set to 1 because it's public.
     *
     *       0 Anode is distributed as a domain1D
     *       1 Anode is a SurDomain_Electrode type
     */
    int anodeType_;

    //! Type of the cathode
    /*!
     *  Determines whether the cathode is distributed or is part of a SurDomain1D
     *  Right now, this defaults to 0, but can be set to 1 because it's public.
     *
     *       0    Cathode is distributed as a domain1D
     *       1    Cathode is a SurDomain_Electrode type
     */
    int cathodeType_;

 protected:
    //! 
    int maxSubGridTimeSteps_;

    //! Heat generation term for the battery sandwidge at the time step t_n
    /*!
     *   units are Joules per time per m2.
     *   It's the heat generation per unit area of the battery at the t_n time step.
     */
    double QdotPerArea_n_;

    //! Heat generation term for the battery sandwidge at the time step t_n
    /*!
     *   units are Joules per time per m2.
     *   It's the heat generation per unit area of the battery at the t_n time step.
     */
    double QdotPerArea_nm1_;

    //! Heat generation term for the battery anode at the time step t_n
    /*!
     *   units are Joules per time per m2.
     *   It's the heat generation per unit area of the battery at the t_n time step.
     */
    double QdotAnodePerArea_n_;

    //! Heat generation term for the battery separator at the time step t_n
    /*!
     *   units are Joules per time per m2.
     *   It's the heat generation per unit area of the battery at the t_n time step.
     */
    double QdotSeparatorPerArea_n_;

    //! Heat generation term for the battery cathode at the time step t_n
    /*!
     *   units are Joules per time per m2.
     *   It's the heat generation per unit area of the battery at the t_n time step.
     */
    double QdotCathodePerArea_n_;

    //! Initial Capacity of the anode per Area
    /*!
     *    units Amps sec m-2
     */
    double capacityAnodePA_;

    //! Initial Capacity of the cathode per Area
    /*!
     *    units Amps sec m-2
     */
    double capacityCathodePA_;

    //! Remaining Capacity of the anode per Area
    /*!
     *    units Amps sec m-2
     */
    double capacityLeftAnodePA_;

    //! Remaining Capacity of the cathode per Area
    /*!
     *    units Amps sec m-2
     */
    double capacityLeftCathodePA_;
    
    //! Total Capacity discharged from the anode per Area
    /*!
     *    units Amps sec m-2
     */
    double capacityDischargedAnodePA_;

    //! Total Capacity discharged from the cathode per Area
    /*!
     *    units Amps sec m-2
     */
    double capacityDischargedCathodePA_;

    //! Current depth of discharge in Amp seconds per cross-sectional area for the anode
    /*!
     *  Report the current depth of discharge. This is roughly equal to the total
     *  number of electrons that has been theoretically discharged from a fully charged state.
     *  For multiple cycles, this becomes the true electron counter for the electrode.
     */
    double dodAnodePA_;

    //! Current depth of discharge in Amp seconds per cross-sectional area for the cathode
    /*!
     *  Report the current depth of discharge. This is roughly equal to the total
     *  number of electrons that has been theoretically discharged from a fully charged state.
     *  For multiple cycles, this becomes the true electron counter for the electrode.
     */
    double dodCathodePA_;

    //! Total effective capacity of the battery
    /*!
     *   This is the minimum of the anode and cathode capacities.
     *    units Amps sec m-2
     */
    double capacityPAEff_;

    //! Total effective capacity left of the battery
    /*!
     *   This is the minimum of the anode and cathode capacities left.
     *    units Amps sec m-2
     */
    double capacityLeftPAEff_;

    //! Time left to discharge the battery fully.
    /*!
     *   timeLeft = capacityLeft / current
     *
     *  Units:   seconds
     */
    double timeLeft_;

    //!  Current value of the Crate.
    /*!
     *  A Crate would mean that the battery would discharge from initial
     *  capacity to empty capacity in 1 hour. If the current I_1 were
     *  defined as the current with a crate of 1, then
     * 
     *      Crate = current / I_1
     *
     *  units = unitless
     */
    double Crate_current_;

    //! Quick determination of the open circuit voltage of the anode
    double ocvAnode_; 

    //! Quick determination of the open circuit voltage of the cathode
    double ocvCathode_;

public:
    //!  If true, the residual equations have a pressure equation and we have Darcy's equation for velocity
    /*!
     *  If true there the total continuity equation is associated with the pressure equation.
     *  And, the axial velocity is associated with the Darcy expresion
     *
     *  If not true the the axial velocity is associated with the pressure equation
     */
    int hasPressureEquation_;

    //! If true we have a separate equation for the gas in the container.
    int hasGasResevoir_;

    //! Anode Bulk Domain Index
    size_t Anode_BDI_;

    //! Separator Bulk Domain Index
    size_t Separator_BDI_;

    //! Cathode Bulk Domain Index
    size_t Cathode_BDI_;

};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

