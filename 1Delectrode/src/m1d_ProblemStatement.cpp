/*
 * $Id: m1d_ProblemStatement.cpp 567 2013-03-21 23:03:11Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"

#include "m1d_globals.h"
#include "m1d_ProblemStatement.h"
#include "m1d_EpetraJac.h"

#include "Epetra_Comm.h"

#include "BlockEntryGlobal.h"
#include "mdp_allo.h"

using namespace std;
using namespace BEInput;
using namespace TKInput;
using namespace mdpUtil;

namespace m1d
{
bool PrintInputFormat = false;

  ProblemStatement * PSinput_ptr = 0;
//=====================================================================================================================
/*
 *
 *  MPEQUIL_INPUT: constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.  
 */
  ProblemStatement::ProblemStatement() : 
    cf_(0),
    commandFile_(""),
    prob_type(0),
    Energy_equation_prob_type_(0),
    SolutionBehavior_printLvl_(4),
    TimeStepper_printLvl_(1),
    NonlinSolver_printLvl_(-1),
    Residual_printLvl_(0),
    I_LinearSolverBlock(0),
    startTime_(0.0),
    endTime_(3000.0),
    initialTimeStep_(1e-8),
    MaxTimeStep_(1.0E6),
    MinTimeStep_(1.0E-15),
    absTol_(1e-8),
    relTol_(1e-3),
    initDefaultNumCVsPerDomain_(10),
    maxNumTimeSteps_(1000)
  {
  }
  //=====================================================================================================================
  ProblemStatement::ProblemStatement(const ProblemStatement &right) : 
    cf_(0),
    commandFile_(""),
    prob_type(0),
    Energy_equation_prob_type_(0),
    SolutionBehavior_printLvl_(4),
    TimeStepper_printLvl_(1),
    NonlinSolver_printLvl_(-1),
    Residual_printLvl_(0),
    I_LinearSolverBlock(0),
    startTime_(0.0),
    endTime_(3000.0),
    initialTimeStep_(1e-8),
    MaxTimeStep_(1.0E6),
    MinTimeStep_(1.0E-15),
    absTol_(1e-8),
    relTol_(1e-3),
    initDefaultNumCVsPerDomain_(10),
    maxNumTimeSteps_(1000)
  {
    *this = right;
  }
  //=====================================================================================================================
  ProblemStatement::~ProblemStatement()
  {
    if (cf_) {
      delete cf_;
      cf_ = 0;
    }
    delete I_LinearSolverBlock; 
  }
 //=====================================================================================================================  
  ProblemStatement& ProblemStatement::operator=(const ProblemStatement &right) {

    if (this == &right) return *this;

    if (cf_) {
      delete (cf_);
    }
    cf_ =  (right.cf_)->duplMyselfAsBlockEntry();

    commandFile_               = right.commandFile_;
    prob_type                  = right.prob_type;
    Energy_equation_prob_type_ = right.Energy_equation_prob_type_;
    SolutionBehavior_printLvl_ = right.SolutionBehavior_printLvl_;
    TimeStepper_printLvl_      = right.TimeStepper_printLvl_;
    NonlinSolver_printLvl_     = right.NonlinSolver_printLvl_;
    Residual_printLvl_         = right.Residual_printLvl_;
    I_LinearSolverBlock        = right.I_LinearSolverBlock;
    startTime_                 = right.startTime_;
    endTime_                   = right.endTime_;
    initialTimeStep_           = right.initialTimeStep_;
    MaxTimeStep_               = right.MaxTimeStep_;
    MinTimeStep_               = right.MinTimeStep_;
    absTol_                    = right.absTol_;
    relTol_                    = right.relTol_;
    initDefaultNumCVsPerDomain_= right.initDefaultNumCVsPerDomain_;
    maxNumTimeSteps_           = right.maxNumTimeSteps_;

    return *this;
  }

//=====================================================================================================================
/*
 * Do any post processing required.
 * This might include unit conversions, opening files, etc.
 */
void 
ProblemStatement::post_process_input()
{
}
  //====================================================================================================================
  //! other preparation steps
void
ProblemStatement::InitForInput()
{
  /*
   * nothing to be done at this level
   */
}
//=====================================================================================================================
void
ProblemStatement::setup_input_pass1(BlockEntry *cf)
{
  /*
   * Nothing to do yet
   */
   BaseEntry::set_SkipUnknownEntries(true);

   /* --------------------------------------------------------------------
    * 
    * Energy Equation Problem Type
    *
    *   0  = None (default)
    *   1  = Fixed
    *   2  = Dirichlet Equation
    *   3  = Enthalpy Equation
    *   4  = Temperature Equation
    */
   const char *energyEqList[5] = {"None", "Fixed", "Dirichlet Equation", "Enthalpy Equation", "Temperature Equation"};

   LE_PickList *lepkm = new LE_PickList("Energy Equation Problem Type", &Energy_equation_prob_type_,
                                        energyEqList, 5, 0, "Energy_equation_prob_type_");
   lepkm->set_default(0);
   cf->addLineEntry(lepkm);
}
//=====================================================================================================================
void
ProblemStatement::setup_input_pass2(BlockEntry *cf)
{

 /* ------------------------------------------------------------------
   * Line Input For the problem type
   *
   */
  int reqd = 0;
  LE_OneInt *i2 = new LE_OneInt("Number of CVs Per Domain", &(initDefaultNumCVsPerDomain_),
			        reqd, "numCVs");
  i2->set_default(10);
  i2->set_limits(1000, 1);
  cf->addLineEntry(i2);
 
}
//=====================================================================================================================
void
ProblemStatement::setup_input_pass3(BlockEntry *cf)
{
  int reqd = 0;
  /* ---------------------------------------------------------------
   *
   */
  LE_OneStr *s1 = new LE_OneStr("Title", &(Title), 100000, 1, 0, "Title");
  cf->addLineEntry(s1);

  /* ------------------------------------------------------------------
   * Line Input For the problem type
   *
   */
  reqd = 1;
  LE_OneInt *i2 = new LE_OneInt("Problem Type", &(prob_type), reqd, "prob_type");
  i2->set_default(0);
  i2->set_limits(1000, 0);
  cf->addLineEntry(i2);
  /*
   */

  /* ------------------------------------------------------------------
   * Line Input For Simulation Start Time
   *
   */
  reqd = 0;
  LE_OneDbl *d1 = new LE_OneDbl("Start Time", &(startTime_), reqd, "startTime_");
  d1->set_default(0.0);
  cf->addLineEntry(d1);

  /* ------------------------------------------------------------------
   * Line Input For Simulation End Time
   *
   */
  reqd = 0;
  LE_OneDbl *d2 = new LE_OneDbl("End Time", &(endTime_), reqd, "endTime_");
  d2->set_default(3000.0);
  cf->addLineEntry(d2);
  /*
   */

  /* ------------------------------------------------------------------
   * Line Input For Initial Time Step
   *
   */
  reqd = 0;
  LE_OneDbl *d3 = new LE_OneDbl("Initial Time Step", &(initialTimeStep_), reqd, "initialTimeStep_");
  d3->set_default(1e-8);
  cf->addLineEntry(d3);
  /*
   */

 /* ------------------------------------------------------------------
   * Line Input For Maximum value of the Time Step
   *
   */
  reqd = 0;
  LE_OneDbl *dmax = new LE_OneDbl("Maximum Time Step", &(MaxTimeStep_), reqd, "MaxTimeStep_");
  dmax->set_default(1e6);
  cf->addLineEntry(dmax);
  /*
   */

  /* ------------------------------------------------------------------
   * Line Input For Minimum value of the Time Step
   *
   */
  reqd = 0;
  LE_OneDbl *dmin = new LE_OneDbl("Minimum Time Step", &(MinTimeStep_), reqd, "MinTimeStep_");
  dmin->set_default(1e-15);
  cf->addLineEntry(dmin);
  /*
   */

  /* ------------------------------------------------------------------
   * Line Input For Absolute Tolerance
   *
   */
  reqd = 0;
  LE_OneDbl *d4 = new LE_OneDbl("Absolute Tolerance", &(absTol_), reqd, "absTol_");
  d4->set_default(1e-8);
  cf->addLineEntry(d4);
  /*
   */

  /* ------------------------------------------------------------------
   * Line Input For Relative Tolerance
   *
   */
  reqd = 0;
  LE_OneDbl *d5 = new LE_OneDbl("Relative Tolerance", &(relTol_), reqd, "relTol_");
  d5->set_default(1e-3);
  cf->addLineEntry(d5);

  /* --------------------------------------------------------------------------------------------------------------------------
   *   Set the level of printing that occurs for each time step
   *
   * Level of solution printing done to stdout
   *
   *   0 -> Don't print anything
   *   1 -> Print only about significant issues going on
   *   2 -> Print status information at regular intervals.
   *   3 -> Print ShowSolution at regular intervals
   *   4 -> Print ShowSolution at all successful time steps
   *   5 -> Print additional information at each time step
   *   6 -> Print some information about each electrode object at each time step
   *   7 -> Print a lot of information about each electrode object at each time step

   */
  reqd = 0;
  LE_OneInt *sb1 = new LE_OneInt("Solution Behavior Print Level", &(SolutionBehavior_printLvl_), reqd, 
				 "SolutionBehavior_printLvl_");
  sb1->set_default(4);
  sb1->set_limits(20, 0);
  cf->addLineEntry(sb1);

  /* -------------------------------------------------------------------------------------------------------------------
   *   Set the level of printing that occurs for each time step
   *
   *  Time Stepper Print Level - int       (optional)
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
  reqd = 0;
  LE_OneInt *ts1 = new LE_OneInt("Time Stepper Print Level", &(TimeStepper_printLvl_), reqd, "TimeStepper_printLvl_");
  ts1->set_default(1);
  ts1->set_limits(20, 0);
  cf->addLineEntry(ts1);

  /* -------------------------------------------------------------------------------------------------------------------
   * Line Input For the nonlinear solver log level
   *
   *  Nonlinear Solver Log Level - int       (optional)
   *
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
  reqd = 0;
  LE_OneInt *i3 = new LE_OneInt("Nonlinear Solver Log Level", &(NonlinSolver_printLvl_), reqd, "NonlinSolver_printLvl");
  i3->set_default(-1);
  i3->set_limits(20, -1);
  cf->addLineEntry(i3);


  /* -------------------------------------------------------------------------------------------------------------------
   * Level of residual information printing done to stdout
   *
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
  reqd = 0;
  LE_OneInt *irp = new LE_OneInt("Residual Print Level", &(Residual_printLvl_), reqd, "Residual_printLvl");
  irp->set_default(0);
  irp->set_limits(20, 0);
  cf->addLineEntry(irp);


  /* -------------------------------------------------------------------------------------------------------------------
   * Line Input for the maximum number of time steps
   *
   *  Maximum Number of Time Steps - int       (optional)
   *
   *  Default = 1000
   */
  reqd = 0;
  LE_OneInt *mts1 = new LE_OneInt("Maximum Number of Time Steps", &(maxNumTimeSteps_), reqd, "maxNumTimeSteps_");
  mts1->set_default(1000);
  mts1->set_limits(1000000000, 0);
  cf->addLineEntry(mts1);

  /* -------------------------------------------------------------------------------------------------------------------
   *   Setup a block for the linear solver 
   */

  I_LinearSolverBlock = EpetraJac::setupMDinput_pass1(cf);


  BaseEntry::set_SkipUnknownEntries(false);
}
//======================================================================================================================
/*
 *  Read the input file
 *
 *  printFlag 0 = no output
 *            1 = error message output
 *            2 = output
 */
bool
ProblemStatement::process_input(BEInput::BlockEntry *cf, std::string fileName, int printFlag)
{
  int my_pid = Comm_ptr->MyPID();
  static int pass = 0;
  pass++;
  cf->ZeroLineCount();
  const TOKEN tok_in;
  TOKEN tok_out;
  FILE *ifp = fopen(fileName.c_str(), "r");
  if (!ifp) {
    if (printFlag) {
      cout << "ERROR can't open file " << fileName << endl;
    }
    return false;
  }
  if (!my_pid && printFlag > 1) {
    printf("==========================================================\n");
    printf(" STARTING PROCESSING COMMAND FILE %s, PASS # %d\n", fileName.c_str(), pass);
    printf("==========================================================\n");
  }
  try {
    /*
     * Call the block read function at the main level
     */
    cf->ZeroLineCount();
    cf->read_block(ifp, &tok_out, &tok_in, 0);
  }
  catch (BI_InputError &bi) {
    /*
     * This catches error messages
     */
    cout << bi.errorMessage() << endl;
    return false;
  }
  fclose(ifp);
  if (!my_pid && printFlag > 1) {
    printf("=========================================================\n");
    printf(" FINISHED PROCESSING COMMAND FILE %s, PASS # %d\n", fileName.c_str(), pass);
    printf("=========================================================\n");
  }
  return true;
}
//======================================================================================================================
static void
print_char(const char letter, const int num)
{
  int i;
  for (i = 0; i < num; i++)
    printf("%c", letter);
}
//=====================================================================================================================
/*
 *
 * mpequil_input():
 *
 *    Input for vcs_Cantera. This routine will combine a text input file
 *    with a cantera xml or cti file to create an equilibrium problem
 *    to solve.
 */
int
ProblemStatement::parse_input_1(std::string commandFile)
{
  int retn = 0;
  PSinput_ptr = this;
  commandFile_ = commandFile;
  int my_pid = Comm_ptr->MyPID();

  int printBIProclevel = 9;
  set_tok_input_print_flag(0);
  if (!my_pid) {
    BlockEntry::set_printProcessedLine(true);
  } else {
    BlockEntry::set_printProcessedLine(false);
  }

  if (!my_pid && printBIProclevel) {
    printf("\n");
    print_char('=', 80);
    printf("\n");
    print_char('=', 20);
    printf(" m1d_input: START OF PROBLEM STATEMENT ");
    print_char('=', 21);
    printf("\n");
    print_char('=', 80);
    printf("\n\n");
  }

  /**
   * Initialize a block input structure for the command file
   */
  if (cf_) {
    delete cf_;
  }
  cf_ = new BlockEntry("command_file");

  /**
   * Setup and process the input deck for first time.
   * -> Might have to print out the input and quit as well.
   */
  setup_input_pass1(cf_);
  if (PrintInputFormat) {
     printf("\n");
     printf("Command Line Format for Pass 1:\n\n");
     cf_->print_usage(1);
   }

  bool ok = process_input(cf_, commandFile, printBIProclevel);
  if (!ok) {
    return -1;
  }


  return retn;
}
  //====================================================================================================================
  // Set up the input and parse for the second pass
  /*
   *  This is where we specify the number of nodes in each bulk region. 
   *
   * @param commandFile String name for the input file to be parsed
   * @return  Return 0 if everything is ok
   */
  int ProblemStatement::parse_input_2()
  {
    int retn  = 0;
    int printBIProclevel = 9;
    /*
     * Setup and process the input deck for second time.
     * -> Might have to print out the input and quit as well.
     */
    setup_input_pass2(cf_);
    if (PrintInputFormat) {
      printf("\n");
      printf("Command Line Format for Pass 2:\n\n");
      cf_->print_usage(1);
    }
    
    bool ok = process_input(cf_, commandFile_, printBIProclevel);
    if (!ok) {
      return -1;
    }

    /*
     * Setup internally for next pass through the input file.
     */
    InitForInput();
    
    return retn;
  }
 //====================================================================================================================
  // Set up the input and parse for the second pass
  /*
   *  This is where we specify the number of nodes in each bulk region. 
   *
   * @param commandFile String name for the input file to be parsed
   * @return  Return 0 if everything is ok
   */
  int ProblemStatement::parse_input_3()
  {
   int retn  = 0; 
   int printBIProclevel = 9;
   int my_pid = Comm_ptr->MyPID();

   /*
    * Setup and process the input deck for second time
    * -> Might have to print out the input and quit as well.
    */
   setup_input_pass3(cf_);
   if (PrintInputFormat) {
     printf("\n");
     printf("Command Line Format for Pass 3:\n\n");
     cf_->print_usage(1);
   }
   /*
    * Process the first pass of the input file ->
    *   We are just after the information needed to initialize
    *   the Cantera structures and size the problem
    */
   bool ok = process_input(cf_, commandFile_, printBIProclevel);
   if (!ok) {
     return -1;
   } 


  /*
   *          Printout the species information: PhaseID's and mole numbers
   */
  if (!my_pid && printBIProclevel) {
    printf("\n");
    print_char('-', 80);
    printf("\n");

    print_char('=', 20);
    printf(" m1d_input: END OF PROBLEM STATEMENT ");
    print_char('=', 23);
    printf("\n");
    print_char('=', 80);
    printf("\n\n");
  }



   post_process_input();


   return retn;
  }
//=====================================================================================================================

}
//=====================================================================================================================
