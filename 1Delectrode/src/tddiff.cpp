/**
 * @file tddiff.cpp
 *
 */

/*
 *  $Id: tddiff.cpp 560 2013-03-06 23:50:04Z hkmoffa $
 */


#include <AztecOO.h>
#include <Ifpack.h>
#include "m1d_defs.h"
#include "mdp_allo.h"

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
#include "m1d_SolNonlinear.h"
#include "m1d_ProblemStatement.h"

#include "BEulerInt.h"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_VbrMatrix.h>
#include "m1d_EpetraJac.h"
#include <Teuchos_ParameterList.hpp>

#include <vector>
#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

#include "tddiff.h"
#include "m1d_solvers.h"

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

#include "m1d_GlobalIndices.h"
#include "m1d_LocalNodeIndices.h"
#include "m1d_exception.h"
#include "m1d_VBRIndices.h"
#include "m1d_ProblemResidEval.h"

using namespace std;
using namespace m1d;
using namespace beuler;

m1d::ProblemResidEval *PS_ptr = 0;

ProblemStatement PSinput;


//m1d::ProblemStatement PSinput;

void
printUsage()
{
  cout << "usage: mpequil [-h] [-help_cmdfile] [-d #] [mpequil.inp]" << endl;
  cout << "    -h           help" << endl;
  cout << "    -d           #   : level of debug printing" << endl;
  cout << "  mpequil.inp    : command file" << endl;
  cout << "                     : (if missing, assume mpequil.inp)" << endl;
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
    std::string commandFile = "tddiff.inp";
    if (argc > 1) {
      string tok;
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
        } else if (commandFile == "" || commandFile == "tddiff.inp") {
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

    m1d::ProblemResidEval *ps = new ProblemResidEval(1.0E-13);
    PS_ptr = ps;

    /*
     *  Initialize the domain structure for the problem
     */

    ps->specifyProblem(2, &PSinput);

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

    BEulerInt t1;

    t1.m_doSpecialStartCalc = 1;
    char *resp_str = getenv("1DBat_DODAEINIT");
    if (resp_str) {
      if (resp_str[0] == 'n') {
        t1.m_doSpecialStartCalc = 0;
      }
    }


    b->PutScalar(0.0);
    soln->PutScalar(0.0);

    ps->domain_prep();
    ps->setAtolVector(1.0E-8, *soln);
   
    EpetraJac *jac = new EpetraJac(*ps);
    jac->allocateMatrix();

#ifdef DEBUG_MATRIX_STRUCTURE
    print0_sync_start(true, w0, Comm);
    ostringstream ssSave;
    jac->queryMatrixStructure(ssSave);
    w0 << ssSave.str();
    w0 << endl;
    ssprint0(w0);
    print0_sync_end(true, w0, Comm);
#endif

    ps->residEval(res, false, soln, 0, 0.0, 0.0);

    print0_epMultiVector(*res, "Residual Value");

    t1.setInitialTimeStep(1.0E-8);
    t1.initializePRE(*ps);

    int min_newt_its = 2;
    int matrixConditioning = 0;
    bool colScaling = 1;
    bool rowScaling = 1;
    t1.setNonLinOptions(min_newt_its, matrixConditioning, colScaling, rowScaling);

    t1.setProblemType(BEULER_JAC_ANAL);



    int printSolnSteps = 3;
    int printSolnInterval = 1;
    int printSolnFirstSteps = 4;
    int dumpJacobians = 1;
    t1.setPrintSolnOptions(printSolnSteps, printSolnInterval, printSolnFirstSteps, dumpJacobians);
 
    int printFlag;
#ifdef DEBUG_HKM
   // printFlag = 10;
    //t1.setPrintFlag(printFlag);
    //SolNonlinear::s_print_NumJac = true;
#endif

    t1.determineInitialConditions(0.0);

    printFlag = 6;
    t1.setPrintFlag(printFlag);

    /*
     *  Ok, if we have turned on the consistent initial conditions capability, the initial
     *  time-dependent residuals should be zero, now. This has been checked out in practice.
     *  Let's print it out to make sure.
     */
    Epetra_Vector_Ghosted & solnBE    = t1.solnVector();
    Epetra_Vector_Ghosted & solnDotBE = t1.solnDotVector();

    ps->setAtolVector_DAEInit(1.0E-4, solnBE, solnDotBE);


    ps->residEval(res, true, &solnBE, &solnDotBE, 0.0, 0.0);
    print0_epMultiVector(*res, "Residual Value");


    const Epetra_Vector_Owned &abstol = ps->atolVector();
    t1.setTolerancesEpetra(PSinput.relTol_, abstol);

    int numInitialConstantDeltaTSteps = 3;
    t1.setNumInitialConstantDeltaTSteps(numInitialConstantDeltaTSteps);


    printSolnSteps = 3;
    printSolnInterval = 1;
    printSolnFirstSteps = 4;
    dumpJacobians = 0;
    t1.setPrintSolnOptions(printSolnSteps, printSolnInterval, printSolnFirstSteps, dumpJacobians);

  
    t1.setMaxNumTimeSteps(PSinput.maxNumTimeSteps_);

    /*
     * Start the simulation
     */
    //double tstart = second();
    double TFinal = 100;
    t1.integratePRE(TFinal);


    /*
     * Cleanup
     */
    safeDelete(jac);
    safeDelete(ps);

    // safeDelete(snp);
    safeDelete(b);
    safeDelete(v);
    safeDelete(soln);
    safeDelete(res);
    Cantera::appdelete();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
  }
  /*
   *  Epetra (and trilinos) throw errors of type int!!!!
   *  This is kind of unusual. However, we can catch these errors
   *  and print out something useful, rather than letting them go uncaught
   *  to the operating system.
   */
  catch (int eNum) {
    cerr << "main:: Caught Trilinos integer error: " << eNum << "\n\n" << endl;
  }
  catch (m1d::m1d_Error &mE) {
    cerr << "caught an error\n" << endl;
    m1d::showErrors(std::cerr);
  }

}

//==============================================================================
void
generateDomain1D(ProblemResidEval *ps)
{
  (ps->DL_ptr_)->generateDomain1D(ps);
}
//==================================================================================

//==================================================================================
//==================================================================================


