/**
 * @file tddiff.cpp
 *
 */

/*
 *  $Id: LiKCl_infPorousBat.cpp 504 2013-01-07 22:32:48Z hkmoffa $
 */


#include <cantera/transport.h>      // transport properties
#include <cantera/thermo.h>      // transport properties
#include <cantera/thermo/IonsFromNeutralVPSSTP.h>  // ion properties
#include <cantera/thermo/StoichSubstance.h>  // separator

#include "Ifpack.h"
#include "AztecOO.h"

#include "m1d_defs.h"

#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "m1d_defs.h"
#include "m1d_Comm.h"
#include "m1d_EpetraExtras.h"
#include "m1d_globals.h"
#include "m1d_BulkDomain1D.h"

#include "m1d_DomainLayout.h"
#include "m1d_SolGlobalNonlinear.h"
#include "m1d_SolNonlinear_CurrentSolve.h"
#include "m1d_ProblemStatement.h"
#include "m1d_DomainLayout_LiKCl_infPorousBat.h"

#include "m1d_ProblemStatementCell.h"
#include "m1d_CanteraElectrodeGlobals.h"

#include "m1d_BEulerInt.h"

#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"


#include "LiKCl_infPorousBat.h"
#include "m1d_solvers.h"

#include "m1d_GlobalIndices.h"
#include "m1d_LocalNodeIndices.h"
#include "m1d_exception.h"
#include "m1d_VBRIndices.h"


//=====================================================================================
/*
 * HKM
 *
 *     Formulation of the structure of the problem.
 *
 *
 *    We will create a system for solving one dimensional problems based
 *    on using Trilinos and Cantera.
 *
 *    The matrix stencil for the calculations will be based on  a set of
 *    of domains that are connected with each other by tie-regions.
 *    The tie-regions may have a set of equations associated with
 *    the tie region. This set of equations will have dependencies
 *    on the neighboring nodes of the bounding domains.
 *
 *
 *
 *    |  Tie        |  Domain  | Tie      |  Domain  |  Tie     |
 *    |  Region     |  Zero    | Region   |  One     |  Region  |
 *    |  0          |          | 1        |          |  2       |
 *    |             |          |          |          |          |
 *    |  Ne^T_0     |  Ne^D_0  | Ne^T_1   |  Ne^D_1  |  Ne_T_2  |
 *
 *
 *
 *
 *   The system will be able to do solved using a MP formulation. However,
 *   we will not build into the formulation any kind of extensibility. All
 *   nodes will know everything about the problem. All nodes will have all
 *   data. MP will be achieved by breaking the nodes of the problem up
 *   into sections. The residual for each section will belong to one processor
 *   The matrix row for that section will belong to one processor.
 *   Matrix solves will occur via a MP method.
 *
 *
 *
 *  For our first step we will assume one domain, with no tie regions !
 */

using namespace std;
using namespace m1d;
using namespace beuler;

m1d::BatteryResidEval *PS_ptr = 0;

m1d::ProblemStatementCell PSinput;

int flagPrecipitation = -1;

void
printUsage()
{
  cout << "usage: LiKCl_infPlateBat [-h] [-help_cmdfile] [-d #] LiKCl_infPlateBat.inp" << endl;
  cout << "    -h           help" << endl;
  cout << "    -d           #   : level of debug printing" << endl;
  cout << "  separator.inp  #   : command file" << endl;
  cout << endl;
}

//=====================================================================================

int
main(int argc, char **argv)
{

  ofstream oftmp;
  m1d::stream0 w0;
#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc, &argv);
#endif
  try {
#ifdef HAVE_MPI
    // Create the One New communications object for global communications
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif
    Comm_ptr = &Comm;
    /*
     * Process the command line arguments
     */
    int retn = 0;
    std::string commandFile = "LiKCl_infPorousBat.inp";
    if (argc > 1) {
      std::string tok;
      for (int j = 1; j < argc; j++) {
        tok = string(argv[j]);
        if (tok[0] == '-') {
          int nopt = static_cast<int> (tok.size());
          for (int n = 1; n < nopt; n++) {
            if (!strcmp(tok.c_str() + 1, "help_cmdfile")) {
              m1d::PrintInputFormat = true;
              break;
            } else if (tok[n] == 'h') {
              printUsage();
              exit(1);
            } else if (tok[n] == 'd') {
              if (j < (argc - 1)) {
                string tokla = string(argv[j + 1]);
                if (strlen(tokla.c_str()) > 0) {
                  n = nopt - 1;
                  j += 1;
                }
              }
            } else {
              printUsage();
              exit(1);
            }
          }
        } else if (commandFile == "" || commandFile == "LiKCl_infPorousBat.inp") {
          commandFile = tok;
        } else {
          printUsage();
          exit(1);
        }
      }
    }

    /*
     * Go get the problem description from the input file
     */
    PSCinput_ptr = &PSinput;
    retn = PSinput.parse_input_1(commandFile);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }
    retn = PSinput.parse_input_2();
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }
    retn = PSinput.parse_input_3();
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

    PSinput.readAnodeInputFile();
    PSinput.readCathodeInputFile();

    m1d::BatteryResidEval *ps = new BatteryResidEval(1.0E-13);
    PS_ptr = ps;

    /*
     *  Initialize the domain structure for the problem
     * Some hard coded things here until we get the parser
     * integrated.
     */

    //not used?
    /* map boundary conditions to their string label */
    /* since I haven't yet invented boundary conditions, 
     * we'll just let them be ints 
     */
    map<string, int> boundaryMap;
    boundaryMap["AnodeCollectorPlate"] = 0;
    boundaryMap["CathodeCollectorPlate"] = 1;

    /*
     *  In this step we assign all of the domains in the problem.
     *  We assign the number of nodes per domain.
     *  In setting up the domains, we assign the number of unknowns
     *  per node within the domain.
     */
    DomainLayout *dl = new DomainLayout_LiKCl_infPorousBat(&PSinput);

    /*
     *  Here, we assign the domainlayout to the problem solver
     *  object.
     */
    ps->specifyProblem(dl, &PSinput);

    /*
     *  Because we know the number of nodes and the number of 
     *  unknowns per node, we can now lay out the global indexing
     *  for the problem.
     */
    ps->generateGlobalIndices();

    // (Now that MPI is initialized, we can access/use MPI_COMM_WORLD.)
    m1d::GlobalIndices &GI = *ps->GI_ptr_;
#ifdef HAVE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &(GI.NumProc));
    MPI_Comm_rank(MPI_COMM_WORLD, &(GI.MyProcID));
#else
    GI.NumProc = 1;
    GI.MyProcID = 0;
#endif

    /*
     *  Generate and fill up the local node vectors on this processor
     *    These indices have to do with accessing the local nodes on the processor
     *    At the end of this routine, the object LocalNodeIncides is fully formed
     *    on this processor.
     */
    ps->generateLocalIndices();

    /*
     * Malloc and start to prep the Domains that will be used for the efficient
     * calculation of the residual
     */
    generateDomain1D(ps);

    if (m1d::PrintInputFormat) {
      exit(0);
    }

    /*
     *  Set the reference temperature for all domains
     *  Do this after we have malloced the Domain1D structures
     */
    dl->set_TP_Reference(PSinput.TemperatureReference_, PSinput.PressureReference_);

    // bool copyMode = false;
    // ps->LRN_VBR_ptr_ = new m1d::LocalRowNodeVBRIndices(Comm_ptr, copyMode, ps->GI_ptr_, ps->LI_ptr_);

    // Declare a matrix and two vectors, to form a linear system: A * v = b.

    // The following will cause an error in the importer, because of an internal error in epetra
    //   -> const ElementSize is not handled correctly.
    // this will cause an error in the importer -> Epetra_Vector *v = new Epetra_Vector(*GI.GbEqnstoOwnedLcEqnsMap);
    Epetra_Vector *v = new Epetra_Vector(*((ps->LI_ptr_)->GbBlockNodeEqnstoLcBlockNodeEqnsColMap));
    Epetra_Vector *b = new   Epetra_Vector(*((ps->LI_ptr_)->GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap));
    Epetra_Vector *res = new Epetra_Vector(*((ps->LI_ptr_)->GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap));
    Epetra_Vector *soln = new Epetra_Vector(*((ps->LI_ptr_)->GbBlockNodeEqnstoLcBlockNodeEqnsColMap));

    v->PutScalar(0.0);

    b->PutScalar(0.0);
    soln->PutScalar(0.0);

    ps->domain_prep();

    EpetraJac *jac = new EpetraJac(*ps);
    jac->allocateMatrix();
#ifdef DEBUG_MATRIX_STRUCTURE
    // this block always produces different results on every call. Therefore, it's not appropriate
    // for a test suite.
    print0_sync_start(true, w0, Comm);
    ostringstream ssSave;
    jac->queryMatrixStructure(ssSave);
    w0 << ssSave.str();
    w0 << endl;
    ssprint0(w0);
    print0_sync_end(true, w0, Comm);
#endif
    ps->residSetupTmps();

    double t_init = 0.0;
    double delta_t = 1.0E-8;
    double delta_t_np1;
    ps->initialConditions(false, soln, 0, t_init, delta_t, delta_t_np1);
    ps->residEval(res, false, soln, 0, 0.0, 0.0);

    print0_epMultiVector(*res, "Residual Value");

    print0_epMultiVector(*soln, "Solution Values");

    BEulerInt t1;
#ifdef DO_INIT_CALC
    char *resp_str = getenv("1DBat_DODAEINIT");
    if (resp_str) {
      if (resp_str[0] == 'y') {
        t1.m_doSpecialStartCalc = 1;
      }
    }
#endif

    t1.initializePRE(*ps);
    t1.determineInitialConditions(PSinput.startTime_, delta_t);

    if (PSinput.rootFinderForConstantCurrent_ > 0) {
       SolGlobalNonlinear *sn = new SolNonlinear_CurrentSolve();
       t1.specifyNonLinearSolver(sn);
    }

    Epetra_Vector_Ghosted &solnInt = t1.solnVector();

    // Ok after the initialization of the solution vector, redo the calculation
    ps->residEval(res, false, &solnInt, 0, 0.0, 0.0);

    print0_epMultiVector(*res, "Residual Value");

    print0_epMultiVector(solnInt, "Solution Values");

    t1.setProblemType(BEULER_JAC_ANAL);

    //Epetra_Vector *av = new Epetra_Vector(*(t1.m_ownedMap));
    //av->PutScalar(PSinput.absTol_);
    //t1.setTolerancesEpetra(PSinput.relTol_, *av);
    //delete av;
    //av = 0;
    
    const Epetra_Vector_Owned &  abstol = ps->atolVector();
    t1.setTolerancesEpetra(PSinput.relTol_, abstol);

    /***************************************************************************/
    // debugging section to do a matrix formation and a solve

    ps->residEval(res, false, &solnInt, 0, 0.0, 0.0);
    jac->matrixEval(false, &solnInt, 0, 0.0, 0.0, Zuzax::Solve_Type::SteadyState_Solve);

    print0_epVbrMatrix(*(jac->A_));

    for (int i = 0; i < ps->LI_ptr_->NumLcOwnedEqns; i++) {
      (*b)[i] = -(*res)[i];
    }

    // Call a function to solve the linear system using an iterative method from AztecOO.
    double linNorm;
    int linearIts;
    jac->solve(b, v, linearIts, linNorm, true);

    print0_epBlockMap((*soln).Map());
    print0_epBlockMap((*v).Map());
    /***************************************************************************/

    t1.setInitialTimeStep( PSinput.initialTimeStep_ );
    int numInitialConstantDeltaTSteps = 3;
    t1.setNumInitialConstantDeltaTSteps(numInitialConstantDeltaTSteps);
    int printNumToTout = 30;
    int printSolnInterval = 1;
    int printSolnFirstSteps = 30;
    int dumpJacobians = 0;
    t1.setPrintSolnOptions(printSolnInterval, printNumToTout, printSolnFirstSteps, dumpJacobians);
    int min_newt_its = 2;
    int matrixConditioning = 0;
    bool colScaling = 1;
    bool rowScaling = 1;
    t1.setNonLinOptions(min_newt_its, matrixConditioning, colScaling, rowScaling);
    t1.setPrintFlag(PSinput.TimeStepper_printLvl_, PSinput.NonlinSolver_printLvl_);
    t1.setMaxNumTimeSteps(PSinput.maxNumTimeSteps_);

    /*
     * Start the simulation
     */
    //double tstart = second();
    set<double> stepTimes;
    stepTimes.insert( PSinput.endTime_ );
    set<double>::iterator step;
    
    //cycle through step times and add to the set of stepTimes
    if ( PSinput.BC_TimeDep_ ) {
      BoundaryCondition *BC = PSinput.BC_TimeDep_ ;
      double nextTime = BC->nextStep();
      //cycle through step times and add to the set of stepTimes
      do { 
	if ( nextTime < PSinput.endTime_ ) {
	  stepTimes.insert( nextTime );
	}
      } while ( ( nextTime = BC->nextStep() ) > 0 );  
      BC->resetSteps();
    }

    //double TInit = PSinput.startTime_ ;
    double TFinal = PSinput.startTime_ ;
    double Tstop;

    for ( step = stepTimes.begin(); step != stepTimes.end(); step++ )  {

      if (PSinput.initialTimeStep_ > 0.0) {
          t1.setInitialTimeStep(PSinput.initialTimeStep_);
      } else {
          t1.setInitialTimeStep(delta_t_np1);
      }
      
      fprintf(stderr, "BOUNDARY CONDITION time step until %f\n", *step );

      TFinal = *step;
      
      // For step-change BC, the predictor corrector requires 
      // values at both ends of the interval.  One of these 
      // will be on the other end of the step.  We overcome 
      // this by reducing the interval by 1% of the initialTimeStep_
      // We should actually base this on a roundoff error of the TFinal.
      if ( PSinput.cathodeBCType_  == 6 
	   || PSinput.cathodeBCType_  == 7 ) {
	double smallTimeShift = std::max( 1e-4 * PSinput.initialTimeStep_, 1e-10 * TFinal );
	//double smallTimeShift = 1e-10 * TFinal;
	Tstop = TFinal - smallTimeShift;
	Tstop = t1.integratePRE( Tstop );
      } 
      else {
	Tstop = t1.integratePRE( TFinal );
      }
      
    }

    //double tend = second();
    // fprintf(stderr, "Total time = %g seconds\n", tend - tstart);
    //printf("Total time = %g seconds\n", tend - tstart);


    // safeDelete(wAll);

    /*
     * Cleanup
     */
    ZZCantera::appdelete();
    safeDelete(jac);
    safeDelete(ps);

    safeDelete(b);
    safeDelete(v);
    safeDelete(soln);
    safeDelete(res);

#undef PRECIPITATE
#ifdef  PRECIPITATE
    if ( PSinput.useDakota_ ) {
      std::ofstream dfp;
      dfp.open( PSinput.toDakotaFileName_.c_str() , std::ios_base::out );
      dfp << flagPrecipitation << std::endl;
      dfp.close();
    }
#endif // PRECIPITATE
#undef PRECIPITATE

    ZZCantera::appdelete();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
  }
  /*
   *  If we have thrown a CanteraError of some sort, we catch it here and print out
   *  an informative error message.
   */
  catch (ZZCantera::CanteraError) {
    ZZCantera::showErrors();
    return -1;
  }
  /*
   *  Epetra (and trilinos) throw errors of type int!!!!
   *  This is kind of unusual. However, we can catch these errors
   *  and print out something useful, rather than letting them go uncaught
   *  to the operating system.
   */
  catch (int eNum) {
    cerr << "main:: Caught Trilinos integer error: " << eNum << "\n\n" << endl;
    return eNum;
  }
  catch (m1d::m1d_Error &mE) {
    cerr << "caught an error\n" << endl;
    m1d::showErrors(std::cerr);
    return -1;
  }
  return 0;
}

//==============================================================================
void
generateDomain1D(BatteryResidEval *ps)
{
  (ps->DL_ptr_)->generateDomain1D(ps);
}
//==================================================================================

//==================================================================================
//==================================================================================


