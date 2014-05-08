/**
 * @file Multi1DDomain.cpp
 *
 */

/*
 *  $Id: Multi1DDomain.cpp 504 2013-01-07 22:32:48Z hkmoffa $
 */

#include "m1d_defs.h"

#include "../../config.h"

#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "Multi1DDomain.h"
#include "m1d_EpetraExtras.h"
#include "m1d_globals.h"
#include "m1d_DomainLayout.h"
#include "m1d_SolNonlinear.h"
#include "m1d_ProblemStatement.h"
#include "m1d_GlobalIndices.h"

int
alloc_double_matrix_3d(double ****double_matrix_3d, int L, int M, int N);

void
ex_write_output_file(M1D_MPI_Comm mpi_comm, Epetra_CrsMatrix*& A, Epetra_Vector*& v, Epetra_Vector*& b);

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

m1d::ProblemResidEval *PS_ptr = 0;

m1d::ProblemStatement PSinput;

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
  int linearIts;
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
    std::string commandFile = "Multi1DDomain.inp";
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
        } else if (commandFile == "" || commandFile == "Multi1DDomain.inp") {
          commandFile = tok;
        } else {
          printUsage();
          exit(1);
        }
      }
    }

    /*
     * Go get the problem description from the input file
     *   Read what the problem is
     */
    retn = PSinput.parse_input_1(commandFile);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }
 

    m1d::ProblemResidEval *ps = new ProblemResidEval(1.0E-13);
    PS_ptr = ps;

    /*
     *  Second call to the input file, read the number of nodes.
     *   read whether this is a restart file or not.
     */
    retn = PSinput.parse_input_2();
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }


    /*
     *  Initialize the domain structure for the problem, specify the number of nodes
     */

    ps->specifyProblem(1, &PSinput);

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

    /*
     *  Read the rest of the input
     */
    retn = PSinput.parse_input_3();
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

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

    b->PutScalar(0.0);
    soln->PutScalar(0.0);

    ps->domain_prep();

    EpetraJac *jac = new EpetraJac(*ps);
    jac->allocateMatrix();
    jac->process_BEinput(PSinput.I_LinearSolverBlock);

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

    SolNonlinear *snp = new SolNonlinear();

    // Setup the problem for solution.
    Solve_Type_Enum stype = SteadyState_Solve;
    snp->setup_problem(stype, soln, soln, 0.0, *ps, *jac);

    int num_newt_its = 0;
    int num_linear_solves = 0;
    int num_backtracks = 0;

    snp->setRowScaling(true);
    snp->setColScaling(true, 2);
    snp->setPrintFlag(PSinput.NonlinSolver_printLvl_);
    int rtn = snp->solve_nonlinear_problem(stype, soln, 0, 0.0, 0.0, num_newt_its, num_linear_solves, num_backtracks);
    if (!rtn) {
      printf("exception code from nonlinear solve\n");
    }

    std::string fn = "solution";

    Epetra_DataAccess eee = View;
    double *V;
    soln->ExtractView(&V);
    Epetra_Vector *soln_owned = new Epetra_Vector(eee, 
						  *((ps->LI_ptr_)->GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap),
						  V);

    ps->saveSolutionEnd(0, fn, *soln_owned, 0, 0.0, 0.0);

    ps->residEval(res, false, soln, 0, 0.0, 0.0);
    jac->matrixEval(false, soln, 0, 0.0, 0.0, SteadyState_Solve);

    print0_epVbrMatrix(*(jac->A_));

    for (int i = 0; i < ps->LI_ptr_->NumLcOwnedEqns; i++) {
      (*b)[i] = -(*res)[i];
    }

    // Call a function to solve the linear system using an iterative method from AztecOO.
    double linNorm;
    jac->solve(b, v, linearIts, linNorm, true);

    print0_epBlockMap((*soln).Map());
    print0_epBlockMap((*v).Map());
    //b = v;

    /*
     * Update the non-ghost nodes
     */
    for (int i = 0; i < ps->LI_ptr_->NumLcOwnedEqns; i++) {
      (*soln)[i] += (*v)[i];
    }

    print0_epMultiVector(*soln);

    // Write the solution vector to an output mesh.
    // The code in ex_write_output_file() does a gather and single-file write,
    //   so no need to serialize the write operations.
    //  ex_write_output_file(MPI_COMM_WORLD, A, v, b);

    print0_sync_start(true, w0, Comm);
    w0.print0("proc %d writes an int %d\n", GI.MyProcID, 10);
    w0 << "Solution vector proc (" << GI.MyProcID << ")" << std::endl;

    ssprint0(w0);

    //Print0_mv(w0, *v);
    w0 << *soln;

    print0_sync_end(true, w0, Comm);

    ps->residEval(res, false, soln, 0, 0.0, 0.0);

    print0_epMultiVector(*res, "Residual Value Again");

    printOn0(*soln, &Comm);

    /*
     * Let's test out the distributed to GlobalAll communications routine
     */
    Epetra_Vector *wAll = gatherOnAll(*soln, &Comm);

    print0_sync_start(true, w0, Comm);
    Print0_epMultiVector(w0, *wAll);
    print0_sync_end(true, w0, Comm);

    safeDelete(wAll);

    /*
     * Cleanup
     */
    safeDelete(soln_owned);
    safeDelete(jac);
    safeDelete(ps);

    safeDelete(snp);
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


