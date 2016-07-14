/*
 * $Id: InterfacialMassTransfer_1to1Distrib.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _IMTPSS_NOSURF_H
#define _IMTPSS_NOSURF_H


#include "InterfacialMassTransfer_PseudoSS.h"

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



#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif 
{



  class imtPSS_NoSurf: public InterfacialMassTransfer_PseudoSS {
  public:
        
    // ---------------------------------------------------------------------------------------------
    // ----------------------- BASIC SETUP ROUTINES  -----------------------------------------------
    // ---------------------------------------------------------------------------------------------

    //! Constructor
    imtPSS_NoSurf();

    //! Destructor
    virtual ~imtPSS_NoSurf();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    imtPSS_NoSurf(const imtPSS_NoSurf &right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    imtPSS_NoSurf & operator=(const imtPSS_NoSurf &right);

    //! Return the type of interfacial mass transport object 
    /*!
     *  Returns the enum type of the object. This is used in the factory routine.
     *
     *  @return Returns an enum type, called   IMT_Types_Enum
     */
    virtual IMT_Types_Enum imtType() const;


    //! Set the electrode ID information
    void setID(int domainNum, int cellNum);

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


    // ------------------------------ SPECIFY BASIC THERMO CONDITIONS  ------------------------------------
  
    // ------------------------------ SPECIFY PROBLEM PARAMETERS ------------------------------------

    // ---------------------------------------------------------------------------------------------
    // ----------------------------- CARRY OUT INTEGRATIONS -----------------------------------------
    // ---------------------------------------------------------------------------------------------

    //! Predict the solution
    /*!
     *  (virtual from InterfacialMassTransfer_Integrator)
     *
     * Ok at this point we have a time step 
     * and initial conditions consisting of phaseMoles_init_ and spMF_init_.
     * We now calculate predicted solution components from these conditions.
     * Additionally, we evaluate whether any multispecies phases are going to pop into existence,
     * adding in potential seed values, or die, creating a different equation set and time step equation.
     * We set the phaseExistence flags in the kinetics solver to reflect phase pops.
     *
     * @return   Returns the success of the operation
     *                 1  A predicted solution is achieved 
     *                 2  A predicted solution with a multispecies phase pop is acheived
     *                 0  A predicted solution is not achieved, but go ahead anyway
     *                -1  The predictor suggests that the time step be reduced and a retry occur.
     */
    virtual int predictSoln();
 
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
    void  resetStartingCondition(double Tinitial);


    //! Pack the nonlinear solver problem
    /*!
     *  formulate the nonlinear solver problem to be solved.
     *     Fields to be filled in
     *             yvalNLS_
     *             ylowNLS_
     *             yhighNLS_
     *             atolNLS_
     *             deltaBoundsMagnitudesNLS_     
     */
    virtual void initialPackSolver_nonlinFunction();

    //! Unpack the soln vector
    /*!
     *  (virtual from InterfacialMassTransfer_Integrator)
     *
     *  This function unpacks the solution vector into  phaseMoles_final_,  spMoles_final_, and spMf_final_[]
     */
    virtual void unpackNonlinSolnVector(const double * const y);

    //! calculate the residual
    /*!
     *   
     *  (virtual fucntion from InterfacialMassTransfer_Integrator)
     *
     *  @return Returns 1 if everything is ok.
     *          Returns 0 if the current conditions can not be calculated.
     */
    virtual int calcResid(doublereal * const resid, const ResidEval_Type_Enum evalType);

    // ----------------------------------------------------------------------------------------------
    // ----------------------------- GET CONDITIONS OUT --------------------------------------------
    // ----------------------------------------------------------------------------------------------
   //
    //       (unless specified this is always at the final conditions and time
    //

    // ----------------------------- GET INSTANTANEOUS SOURCE TERMS --------------------------------


    // ---------------------------- GET INTEGRATED SOURCE TERMS -------------------------------------

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

    // --------------------------- GET MOLE NUMBERS ------------------------------------------------


    // --------------------------- GET THERMO AND VOLTAGE  ------------------------------------------------

    // ------------------------- GET VOLUMES -----------------------------------------------------------
    
 

    // ---------------------- GET SURFACE AREAS -------------------------------------------------------





    // --------------------------- INTERNAL UPDATE FUNCTIONS --------------------------------------


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


    // ---------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------
  protected:


    // ------------------------------ UNKNOWNS REPRESENTING THE INTERFACE PROBLEM --------------------------------------



  };

}


#endif 
/*****************************************************************************/
