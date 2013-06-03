/*
 * $Id: InterfacialMassTransfer_1to1Distrib.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _INTERFACIALMASSTRANSFER_PSEUDOSS_H
#define _INTERFACIALMASSTRANSFER_PSEUDOSS_H


#include "InterfacialMassTransfer_Integrator.h"

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


  //! Class that assumes you have a pseudo steady state treatment of the integrated equations
  /*!
   *  This approximation means that the surface relaxation times are fast compared to the
   *  volumetric relaxation times. It also means that there are no state variables in the
   *  object that have to be tracked as a function of time that change over time.
   *
   *  These assumptions then mean that integrations in time are greatly simplified. Therefore,
   *  an assumption that pseudo steady state within that routine can be used. 
   */
  class InterfacialMassTransfer_PseudoSS : public InterfacialMassTransfer_Integrator {
  public:
 
    // ---------------------------------------------------------------------------------------------
    // ----------------------- BASIC SETUP ROUTINES  -----------------------------------------------
    // ---------------------------------------------------------------------------------------------

    //! Constructor
    InterfacialMassTransfer_PseudoSS();

    //! Destructor
    virtual ~InterfacialMassTransfer_PseudoSS();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    InterfacialMassTransfer_PseudoSS(const InterfacialMassTransfer_PseudoSS &right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    InterfacialMassTransfer_PseudoSS & operator=(const InterfacialMassTransfer_PseudoSS &right);

    //!  Setup the electrode
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


    // ------------------------------ SPECIFY PROBLEM PARAMETERS ------------------------------------

    // ---------------------------------------------------------------------------------------------
    // ----------------------------- CARRY OUT INTEGRATIONS -----------------------------------------
    // ---------------------------------------------------------------------------------------------

    //! Report the integrated source term for the interface over an interval in time.
    /*!
     *  This is the net change in the moles of species defined in the electrode over that
     *  interval of time. The conditions at the end of the interval are used to carry
     *  out the integrations.
     *  
     *  @param spMoleDelta The end result in terms of the change in moles of species in the
     *                     interfacial object
     *
     *  @return Tfinal    Final time to integrate to.
     *                       
     */
    double integratedSourceTerm(doublereal* const spMoleDelta);

    //! Calculate the integrated source term for the electrode over an interval in time.
    /*!
     *  This is the net change in the moles of species defined in the electrode over that
     *  interval of time. The conditions at the beginning of the interval are used to carry
     *  out the integrations. This may be used at the start of the time step, since no
     *  unknown conditions are required.
     *
     *  @param deltaT     time to integrate
     *  @param spMoleDelta The end result in terms of the change in moles of species in the
     *                     electrode.
     *
     *  @return Tfinal    Final time to integrate to.
     */
    double integrateAndPredictSourceTerm(doublereal deltaT, doublereal* const spMoleDelta);

  

    //! The internal state of the electrode must be kept for the initial and 
    //! final times of an integration step.
    /*
     *  This function advances the initial state to the final state that was calculated
     *  in the last integration step. If the initial time is input, then the code doesn't advance
     *  or change anything.
     *
     * @param Tinitial   This is the New initial time. This time is compared against the "old"
     *                   final time, to see if there is any problem.
     */
    void resetStartingCondition(double Tinitial);

    // --------------------------- GET MOLE NUMBERS ------------------------------------------------



    // --------------------------- GET THERMO AND VOLTAGE  ------------------------------------------------

 

    // -------------------------  GET VOLUMES -----------------------------------------------------------

    // ---------------------- GET SURFACE AREAS -------------------------------------------------------



    // --------------------------- INTERNAL UPDATE FUNCTIONS --------------------------------------

    //! Take the state (i.e., the final state) within the InterfacialMassTransfer_Model and push it down
    //! to the ThermoPhase Objects
    /*!
     *  We take the values of spMoles_final_[] and propagate them down to the ThermoPhase
     *  objects
     */
    virtual void updateState();

    //! Carry out the pseudo steady state calculation
    /*!
     *
     */
    virtual  void solvePseudoSteadyStateProblem(int ifuncOverride, doublereal timeScaleOverride);


    // ---------------------------------------------------------------------------------------------
    // -------------------------------  SetState Functions -------------------------------------------------------
    // ---------------------------------------------------------------------------------------------

    // ---------------------------------------------------------------------------------------------
    // ------------------------------ PRINT ROUTINES -----------------------------------------------
    // ---------------------------------------------------------------------------------------------
 
   

    void printInterfacialMassTransfer_List(int pSrc, bool subTimeStep);
 
    //! Print conditions of the electrode for the current integration step to stdout
    /*!
     *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
     *                       to the final_final time.
     *                       The default is to print out the source terms
     *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
     *                       time step. The default is to print out the global values 
     */
    void printInterfacialMassTransfer(int pSrc, bool subTimeStep);


    void printInterfacialMassTransfer_Phase(int iph, int pSrc, bool subTimeStep);

    //! Write out CSV tabular data on the integrations
    /*!
     *  The idea is to print out tabular data about each intermediate run and about each
     *  global run
     *
     *  @param itype Type of the data
     *            - 0      Initialization information
     *            - 1      Normal intermediate information
     *            - 2      Normal global information at the end of a global time step
     *            - -1     Failed intermediate calculation, a failure from a nonlinear solver step
     *            - -2     Failed calculation from the predictor step - not necessarily significant.
     */
    void writeCSVData(int itype); 

    //--------------------------------------------------------------------------------------

  protected:

    //! surface chemistry problem
    Cantera::ImplicitSurfChem *isc_prob;

    //! Number of unknowns in the surface problem
    int neq_isc_prob_;


  };

}


#endif 
/*****************************************************************************/
