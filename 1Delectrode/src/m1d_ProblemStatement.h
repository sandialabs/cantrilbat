/**
 *  @file m1d_ProblemStatement.h   Class that stores the user input from an input file.
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _M1D_PROBLEMSTATEMENT_H
#define _M1D_PROBLEMSTATEMENT_H

#include "m1d_ProblemResidEval.h"

#include "BlockEntryGlobal.h"

//----------------------------------------------------------------------------------------------------------------------------------
// Forward entry for BlockEntry class
namespace BEInput
{
class BlockEntry;
}

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! Storage for Command file input
/*!
 *  Complete problem statement is storred here.
 *  This class is virtual as child classes will contain more user input. All data here is public to facilitate processing.
 *
 *  All member functions are virtual. They will be overridden by child classes. However, the child classes
 *  are expected to call the base classes, so that cards and capabilities are added on in an onion-like fashion.
 *
 *  This is the current command file specification of the problem statement.
 */
class ProblemStatement
{
public:
    //! Constructor
    ProblemStatement();

    //! Copy constructor
    /*!
     *  @param[in]       right               Object to be copied
     */
    ProblemStatement(const ProblemStatement& right);

    //! Destructor
    virtual ~ProblemStatement();

    //! Assignment operator
    /*!
     *  @param[in]       right               Object to be copied
     *  @return                              Returns a reference to the current object
     */
    ProblemStatement& operator=(const ProblemStatement& right);

    //! Set-up the input deck cards using a first pass
    /*!
     *  
     *  @param[in,out]       cf                  Pointer to the the block of the input file description
     *                                           BlockEntry is a card and block nested-list of permissible values
     */
    virtual void setup_input_pass1(BEInput::BlockEntry* cf);

    //! Set-up the input deck cards using a second pass
    /*!
     *  
     *  @param[in,out]       cf                  Pointer to the the block of the input file description
     *                                           BlockEntry is a card and block nested-list of permissible values
     */
    virtual void setup_input_pass2(BEInput::BlockEntry* cf);

    //! Set-up the input deck cards using a third pass
    /*!
     *  
     *  @param[in,out]       cf                  Pointer to the the block of the input file description
     *                                           BlockEntry is a card and block nested-list of permissible values
     */
    virtual void setup_input_pass3(BEInput::BlockEntry* cf);

    //! Process each of the passes
    /*!
     *  @param[in]           cf                  Pointer to the the block of the input file description
     *                                           BlockEntry is a card and block nested-list of permissible values
     *
     *  @param[in]           fileName            Name of the input file
     *  @param[in]           printFlag           Amount of printing. 
     *                                                - 0 = no output
     *                                                - 1 = error message output
     *                                                - 2 = output
     *
     *  @return                                  Return true if everything is ok
     */
    virtual bool process_input(BEInput::BlockEntry* cf, std::string fileName, int printFlag);

    //! Set up the input
    /*!
     *  Initial parsing to be done here. Domain layout is specified by the end of this level.
     *
     *  @param[in]           commandFile         String name for the input file to be parsed
     *  @return                                  Return 0 if everything is ok
     */
    virtual int parse_input_1(std::string commandFile);

    //! Set up the input and parse for the second pass
    /*!
     *  This is where we specify the number of nodes in each bulk region.
     *
     *  @return                                  Return 0 if everything is ok
     */
    virtual int parse_input_2();

    //! Set up the input and parse for the third time
    /*!
     *  Further processing may be necessary.
     *
     *  @return                                  Return 0 if everything is ok
     */
    virtual int parse_input_3();

    //! Carry out other preparation steps
    /*!
     *
     */
    virtual void InitForInput();

    //! Do any post processing required.
    /*!
     *  This might include unit conversions, opening files, etc.
     *
     *  @param[in]           cf                  Pointer to the the block of the input file description
     *                                           BlockEntry is a card and block nested-list of permissible values
     */
    virtual void post_process_input(BEInput::BlockEntry* cf);

    //============================================================================================================
    //!        DATA FROM INPUT FILES
    //============================================================================================================

    //! BlockEntry structure for the command file.
    BEInput::BlockEntry* cf_;

    //! Name of the command file
    std::string commandFile_;

    //! Title of the simulation
    std::string Title_;

    //! Default temperature in (Kelvin)
    /*!
     *  defaults to 298.15 K
     */
    double TemperatureReference_;

    //! Default Pressure in (Pascal)
    /*!
     *  defaults to OneAtm
     */
    double PressureReference_;

    //! Integer representing the Problem type.
    /*!
     *   The identity of  what is held constant. Currently,
     *   T and P are held constant, and this input is ignored
     */
    int prob_type_;

    //! Write a start and end file record
    int writeStartEndFile_;

    //! Integer representing the energy equation problem type
    /*!
         *  0 -> isothermal               Don't solve an energy equation (default)
         *  1 -> Fixed Temperature Profile Don't solve an energy equation
         *  2 -> Dirichlet Equation       Solve a Dirichlet equation for temperature.
         *                                This is a way to do the fixed system while keeping the
         *                                matrix structure the same.
         *  3 -> Enthalpy Equation        Solve a full enthalpy equation for the temperature
         *  4 -> Temperature Equation     Solve a Cp dT/dt formulation for the temperature
         */
    int Energy_equation_prob_type_;

    //! Integer representing the solid mechanics problem type
    /*!
     *  0 -> none                     Don't solve an stress-strain relationship for mesh motion
     *  1 -> LinearElastic            Solve for mesh motion using a global simple stress-strain relationship
     */
    int Solid_Mechanics_prob_type_;

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

    //! Level of solution printing done to stdout
    /*!
     *   0 -> Don't print anything
     *   1 -> Print only about significant issues going on
     *   2 -> Print status information at regular intervals.
     *   3 -> Print ShowSolution at regular intervals
     *   4 -> Print ShowSolution at all successful time steps
     *   5 -> Print additional information at each time step
     *   6 -> Print some information about each electrode object at each time step
     *   7 -> Print a lot of information about each electrode object at each time step
     */
    int SolutionBehavior_printLvl_;

    //! Level of printing from the Time Stepper
    /*!
     *     The environmental variable PRE_printFlagEnv must be set to greater than 0 as well
     *
     *   0 -> absolutely nothing is printed for a single time step.
     *   1 -> One line summary per time step
     *   2 -> short description, points of interest
     *   3 -> More printed per time step -> major algorithm issues are displayed
     *   4 -> Additional time step error control information is printed out
     *        One line summary of the nonlinear solve
     *   5 -> Summaries of the nonlinear solve iterates are printed out
     *   6 -> Algorithm information on the nonlinear iterates are printed out
     *   7 -> Additional info on the nonlinear iterates are printed out
     *   8 -> Additional info on the linear solve is printed out.
     *   9 -> Info on a per iterate of the linear solve is printed out.
     */
    int TimeStepper_printLvl_;

    //! Log Level for the nonlinear solver
    //!  Set the level of printing that occurs during the nonlinear solve
    /*!
     *  -1 -> Takes the value that the time stepper printLvl wants
     *   0 -> absolutely nothing is printed for a single time step.
     *   1 -> One line summary per time step
     *   2 -> short description, points of interest
     *   3 -> More printed per time step -> major algorithm issues are displayed
     *   4 -> Additional time step error control information is printed out
     *        One line summary of the nonlinear solve
     *   5 -> Summaries of the nonlinear solve iterates are printed out
     *   6 -> Algorithm information on the nonlinear iterates are printed out
     *   7 -> Additional info on the nonlinear iterates are printed out
     *   8 -> Additional info on the linear solve is printed out.
     *   9 -> Info on a per iterate of the linear solve is printed out.
     */
    int NonlinSolver_printLvl_;

    //! Level of residual information printing done to stdout
    /*!
     *   0 -> Don't print anything
     *   1 -> Print only about significant issues going on
     *   2 -> Print status information at regular intervals.
     *   3 -> Print ShowResidual at regular intervals
     *   4 -> Print ShowResidual at all successful time steps
     *   5 -> Print additional information when ShowSolution is called.
     *   6 -> Print additional information when any residual is called.
     *   7 -> Print a lot of information about each when ShowSolution is called.
     *   8 -> Print a lot of information when any base or show residual is called
     *   9 -> Print a lot of information when any residual is called
     */
    int Residual_printLvl_;

    //!  Record number for the restart capability
    /*!
     *  Defaults to -1, in which case there is no restart.
     *  If the value is greater than the number of records, then the last record is used.
     */
    int restartRecordNumber_;

    //!  String name of the restart solution file
    /*!
     *  Defaults to the name solutionRestart.xml
     */
    std::string restartFileName_;

    //! Pointer to the solver input record block
    RecordTree_base* I_LinearSolverBlock;

    //! Start time for simulation
    double startTime_;

    //! End time for simulation
    double endTime_;

    //! initial time step for simulation
    double initialTimeStep_;

    //! Maximum value of the time step
    double MaxTimeStep_;

    //! Minimum value of the time step
    double MinTimeStep_;

    //! Time step for algebraic discontinuities. Below this time step algebraic unknowns start
    //! to matter less to the time step truncation error.
    double TimeStep_AUD_;

    //! Number of initial time steps to take where the time truncation error tolerances are not checked. Instead
    //! the delta T is uniform
    int  m_numInitialConstantDeltaTSteps;

    //! Relative tolerance for time integration
    double absTol_;

    //! Relative tolerance for time integration
    double relTol_;

    //! Initial number of cells in the domain
    int initDefaultNumCVsPerDomain_;

    //! Maximum Number of time steps to be taken
    int maxNumTimeSteps_;

    //! The type of coordinate system that is used
    /*!
     *  There are two that are envisioned: Rectinear_Coordinates and Cylindrical_Coordinates
     */
    CoordinateSystem_Type_Enum coordinateSystemType_;

    //! Cross sectional area, if in cartesian coordinates
    /*!
     *  The overwhelming output from the program is on a per-crosssectional area basis.
     *  However, there are some times when the cross-section is needed. This is the place where it is supplied.
     *
     *    units m2
     */
    double crossSectionalArea_;

    //! Cylinder Length, if in cylindrical coordinates
    /*!
     *  The overwhelming output from the program is on a per-crosssectional area basis
     *  However, there are some times when the cross-section is needed. This is the place
     *  where it is supplied.  We assume 2 pi radians always, i.e., a full radius
     *
     *    units m
     */
    double cylinderLength_;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

