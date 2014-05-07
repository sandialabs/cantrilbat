/**
 * @file m1d_solvers.cpp
 *
 */

/*
 * Copywrite 2013 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include <Ifpack.h>
#include <AztecOO.h>
#include "Amesos.h"

#include "Epetra_LinearProblem.h"
#ifdef HAVE_MPI
#include "Epetra_Comm.h"
#endif
#include <Epetra_Map.h>
#include <Epetra_VbrMatrix.h>

#include "m1d_defs.h"
#include "m1d_solvers.h"
#include "m1d_GlobalIndices.h"
#include "m1d_LocalNodeIndices.h"
#include "m1d_VBRIndices.h"
#include <Epetra_VbrRowMatrix.h>
#include "m1d_EpetraJac.h"


using namespace std;

namespace m1d
{

int PrintLevel_LinearSolver = 0;

//=====================================================================================
void
fill_matrix(Epetra_VbrMatrix*& A,
            Epetra_Vector*& v,
            Epetra_Vector*& b,
            m1d::GlobalIndices *gi_ptr,
            m1d::LocalNodeIndices *li_ptr,
            m1d::LocalRowNodeVBRIndices *lrn_vbr_ptr)
{
  bool copyMode = lrn_vbr_ptr->CopyMode;
  int numLcRowNodes = lrn_vbr_ptr->NumLcRowNodes;

  for (int iBlockRow = 0; iBlockRow < numLcRowNodes; iBlockRow++) {
    int GbNode = li_ptr->IndexGbNode_LcNode[iBlockRow];

    int numBlocks = lrn_vbr_ptr->NumColBlocks_LcRNodes[iBlockRow];

    int *bIndexes = lrn_vbr_ptr->IndexLcNodeColBlock_LcRNodes[iBlockRow];
    //int *colSize_row = m1d::LRN_VBR_ptr->ColSizeColBlock_LcRNodes[iBlockRow];
    Epetra_SerialDenseMatrix **rowBlock = lrn_vbr_ptr->BlockMatrices[iBlockRow];
    double value = 0.0;

    /*
     * Set up when to overwrite with dirichlet conditions
     */
    bool doDir = false;
    if (GbNode == 0 || GbNode == gi_ptr->NumGbNodes - 1) {
      doDir = true;
    }

    if (copyMode) {
      A->BeginReplaceMyValues(iBlockRow, numBlocks, bIndexes);
    }

    for (int jColIndex = 0; jColIndex < numBlocks; jColIndex++) {
      int blockColIndex = bIndexes[jColIndex];

      Epetra_SerialDenseMatrix *rowColBlock = rowBlock[jColIndex];

      if (blockColIndex == li_ptr->IDLeftLcNode_LcNode[iBlockRow]) {
        if (doDir) {
          value = 0.0;
        }
        else {
          value = 1.0;
        }

        (*rowColBlock)[0][0] = value;
        if (GbNode != 0) {
          if (copyMode) {
            A->SubmitBlockEntry(*rowColBlock);
          }
        }
      }
      else if (blockColIndex == iBlockRow) {
        // Block Diagonal
        if (doDir) {
          value = 1.0;
        }
        else {
          value = -2.0;
        }
        rowColBlock = rowBlock[jColIndex];
        (*rowColBlock)[0][0] = value;
        if (copyMode) {
          A->SubmitBlockEntry(*rowColBlock);
        }
      }
      else if (blockColIndex == li_ptr->IDRightLcNode_LcNode[iBlockRow]) {
        // Right
        if (doDir) {
          value = 0.0;
        }
        else {
          value = 1.0;
        }
        if (GbNode != gi_ptr->NumGbNodes - 1) {
          rowColBlock = rowBlock[jColIndex];
          (*rowColBlock)[0][0] = value;
          if (copyMode) {
            A->SubmitBlockEntry(*rowColBlock);
          }
        }
      }
      else {
        printf("confused\n");
        exit(-1);
      }
    }

    if (copyMode) {
      A->EndSubmitEntries();
    }
    Epetra_Vector &vv = *v;
    Epetra_Vector &bb = *b;
    value = 0.0;
    vv.ReplaceMyValues(1, &value, &iBlockRow);
    // vv[iBlockRow] = 0.0;
    bb[iBlockRow] = 0.0;

    if (GbNode == 0) {
      bb[iBlockRow] = 1.0;
    }
  }

}

//=====================================================================================
int
solve_aztecoo(Epetra_VbrMatrix* A, Epetra_Vector* v, Epetra_Vector* b, EpetraJac *jac)
{

#ifdef HAVE_MPI
  const Epetra_Comm &cc = A->Comm();
  int mypid = cc.MyPID();
#endif
  // First, let's set up an Ifpack preconditioner to use with AztecOO.
  Ifpack factory;
  //GlobalIndices &GI = *GI_ptr;
  // Using preconditioner_type("ILU") seems to cause a hang
  //   after the call to P->SetParameters(parameters) below.
  //std::string preconditioner_type("ILUT");
  //std::string preconditioner_type("ILU");
  std::string preconditioner_type("Amesos");
  //int overlap = 1;

  //printf("CPU%03i:   Creating Ifpack preconditioner...\n",  GI.my_procID);

  //Ifpack_Preconditioner* P = factory.Create(preconditioner_type, A, overlap);
  Ifpack_Preconditioner* P = factory.Create(preconditioner_type, A);

  // We'll use a Teuchos::ParameterList for passing control parameters to the preconditioner.
  Teuchos::ParameterList parameters;

  //printf( "CPU%03i:   Setting Teuchos parameters...",  GI.my_procID);

  parameters.set("fact: drop tolerance", 1E-6);
  parameters.set("fact: level-of-fill", 3);

  // Now give the parameters to the preconditioner.
  //printf( "CPU%03i:   Passing Teuchos parameters to preconditioner...\n",  GI.my_procID);
  P->SetParameters(parameters);

  // Now have the preconditioner go ahead and compute the factorization
  //printf("CPU%03i:   Preconditioning A matrix...\n",  GI.my_procID);

  P->Initialize();
  P->Compute();

  // Next, create the AztecOO solver.
  AztecOO azoo(A, v, b);

  // Give our preconditioner to AztecOO.
  if (finite(P->Condest()) && fabs(P->Condest()) < 1E6) {
    //printf("CPU%03i:   Passing precondioner %E to AztecOO...\n",  GI.my_procID, P->Condest());

    azoo.SetPrecOperator(P);

  }
  else {
#ifdef HAVE_MPI
    printf("CPU%03i:   Ifpack computed non-finite preconditioner - not skipping, but  maybe we should have.\n",  mypid);
#else
    printf("CPU0:   Ifpack computed non-finite preconditioner - not skipping, but  maybe we should have.\n");
#endif

    azoo.SetPrecOperator(P);
  }

  int max_num_iterations = 10000;
  double tolerance = 1.0E-3;

  // Direct AztecOO to use the GMRES method, with a restart value of 40.
  //printf("CPU%03i:   Setting Aztec options...\n", GI.my_procID);

  azoo.SetAztecOption(AZ_solver, AZ_gmres);
  //	azoo.SetAztecOption(AZ_solver, AZ_bicgstab);
  azoo.SetAztecOption(AZ_kspace, 40);
  // Do we want AztecOO to print out residual norms during the solve?
  //   (If not, use 'AZ_none' instead of '1'.)
  // azoo.SetAztecOption(AZ_output, 1);
  azoo.SetAztecOption(AZ_output, AZ_warnings);

  // Now solve the linear-system.
  //printf("CPU%03i:   Invoking AztecOO Iterate(%i, %.3E)...",
  //	 GI.my_procID, max_num_iterations, tolerance);

  azoo.Iterate(max_num_iterations, tolerance);
  //printf("done after %i iterations.", azoo.NumIters());


  /** Print out the solution vector. **/
  //std::cout << "Solution vector: " << std::endl << *v << std::endl;


  // Destroy the preconditioner that we had the factory create above.
  delete P;

  return 0;
}
//=====================================================================================
int
solve_aztecoo(Epetra_VbrRowMatrix* Arow, Epetra_Vector* v, Epetra_Vector* b, EpetraJac *jac)
{

#ifdef HAVE_MPI
  const Epetra_Comm &cc = Arow->Comm();
  int mypid = cc.MyPID();
#endif
  // First, let's set up an Ifpack preconditioner to use with AztecOO.
  Ifpack factory;
  //GlobalIndices &GI = *GI_ptr;
  // Using preconditioner_type("ILU") seems to cause a hang
  //   after the call to P->SetParameters(parameters) below.
  //std::string preconditioner_type("ILUT");
  //std::string preconditioner_type("ILU");
  std::string preconditioner_type("Amesos");
  //int overlap = 1;

  //printf("CPU%03i:   Creating Ifpack preconditioner...\n",  GI.my_procID);

  //Ifpack_Preconditioner* P = factory.Create(preconditioner_type, A, overlap);
  Ifpack_Preconditioner* P = factory.Create(preconditioner_type, Arow);

  // We'll use a Teuchos::ParameterList for passing control parameters to the preconditioner.
  Teuchos::ParameterList parameters;

  //printf( "CPU%03i:   Setting Teuchos parameters...",  GI.my_procID);

  parameters.set("fact: drop tolerance", 1E-6);
  parameters.set("fact: level-of-fill", 3);

  // Now give the parameters to the preconditioner.
  //printf( "CPU%03i:   Passing Teuchos parameters to preconditioner...\n",  GI.my_procID);
  P->SetParameters(parameters);

  // Now have the preconditioner go ahead and compute the factorization
  //printf("CPU%03i:   Preconditioning A matrix...\n",  GI.my_procID);

  P->Initialize();
  P->Compute();

  // Next, create the AztecOO solver.
  AztecOO azoo(Arow, v, b);

  // Give our preconditioner to AztecOO.
  if (finite(P->Condest()) && fabs(P->Condest()) < 1E6) {
    //printf("CPU%03i:   Passing precondioner %E to AztecOO...\n",  GI.my_procID, P->Condest());

    azoo.SetPrecOperator(P);

  }
  else {
#ifdef HAVE_MPI
    printf("CPU%03i:   Ifpack computed non-finite preconditioner - not skipping, but  maybe we should have.\n",  mypid);
#else
    printf("CPU0:   Ifpack computed non-finite preconditioner - not skipping, but  maybe we should have.\n");
#endif

    azoo.SetPrecOperator(P);
  }

  int max_num_iterations = 10000;
  double tolerance = 1.0E-3;

  // Direct AztecOO to use the GMRES method, with a restart value of 40.
  //printf("CPU%03i:   Setting Aztec options...\n", GI.my_procID);

  azoo.SetAztecOption(AZ_solver, AZ_gmres);
  //	azoo.SetAztecOption(AZ_solver, AZ_bicgstab);
  azoo.SetAztecOption(AZ_kspace, 40);
  // Do we want AztecOO to print out residual norms during the solve?
  //   (If not, use 'AZ_none' instead of '1'.)
  // azoo.SetAztecOption(AZ_output, 1);
  azoo.SetAztecOption(AZ_output, AZ_warnings);

  // Now solve the linear-system.
  //printf("CPU%03i:   Invoking AztecOO Iterate(%i, %.3E)...",
  //	 GI.my_procID, max_num_iterations, tolerance);

  azoo.Iterate(max_num_iterations, tolerance);
  //printf("done after %i iterations.", azoo.NumIters());


  /** Print out the solution vector. **/
  //std::cout << "Solution vector: " << std::endl << *v << std::endl;


  // Destroy the preconditioner that we had the factory create above.
  delete P;

  return 0;
}
//==================================================================================
int
solve_amesos(Epetra_VbrMatrix* A, Epetra_Vector* x, Epetra_Vector* b, EpetraJac *jac)
{
  // Creates an epetra linear problem, contaning matrix
  // A, solution x and rhs b.
  Epetra_LinearProblem Problem(A, x, b);

   // Initializes the Amesos solver. This is the base class for
   // Amesos. It is a pure virtual class (hence objects of this
   // class cannot be allocated, and can exist only as pointers
   // or references).
   //
   Amesos_BaseSolver* Solver = 0;

   // Initializes the Factory. Factory is a function class (a
   // class that contains methods only, no data). Factory
   // will be used to create Amesos_BaseSolver derived objects.
   //
   Amesos Factory;

   // Specifies the solver. String ``SolverType'' can assume one
   // of the following values:
   // - Lapack
   // - Klu
   // - Umfpack
   // - Pardiso
   // - Taucs
   // - Superlu
   // - Superludist
   // - Mumps
   // - Dscpack
   //
   // HKM checked these out  - as working.
   //std::string SolverType = "Klu";
 //  std::string SolverType = "Umfpack";
   std::string SolverType = jac->directSolverName_;
   //std::string SolverType = "Superludist";
   Solver = Factory.Create(SolverType, Problem);
   // Factory.Create() returns 0 if the requested solver
   // is not available
   if (Solver == 0) {
     std::cerr << "Selected solver (" << SolverType << ") is not available" << std::endl;
     std::cerr << "     available solvers: Umfpack(default) Klu, Superludist" << std::endl;
     return(-1);
   }

   // Calling solve to compute the solution. This calls the symbolic
   // factorization and the numeric factorization.
   Solver->Solve();
   // Print out solver timings and get timings in parameter list.
   //Solver->PrintStatus();
   //Solver->PrintTiming();

   // =========================================== //
    // E N D   O F   T H E   A M E S O S   P A R T //
    // =========================================== //

    // Computes ||Ax - b|| //

    double residual = 0.0;
    if ((PrintLevel_LinearSolver) > 1) {

      Epetra_Vector Ax(b->Map());
      A->Multiply(false, *x, Ax);
      Ax.Update(1.0, *b, -1.0);
      Ax.Norm2(&residual);

#ifdef HAVE_MPI
      const Epetra_Comm &cc = x->Comm();
      int mypid = cc.MyPID();
      if (mypid == 0) {
        std::cout << "Successful Amesos solution: ||b - A * x||_2 = " << residual << std::endl;
      }
#else
        std::cout << "Successful Amesos solution: ||b - A * x||_2 = " << residual << std::endl;
#endif
    }

    // delete Solver. Do this before calling MPI_Finalize() because
    // MPI calls can occur.
    delete Solver;

    if (residual > 1e-5) {
      return(1);
    }
    return 0;
}
//==================================================================================
int
solve_amesos(Epetra_VbrRowMatrix* Arow, Epetra_Vector* x, Epetra_Vector* b, EpetraJac *jac)
{

  /*
   * 
   */

  double **xPoint;
  x->ExtractView(&xPoint);
  //Epetra_Map *acol = new Epetra_Map(Arow->RowMatrixColMap()); 

  double **bPoint;
  b->ExtractView(&bPoint);
  const Epetra_Map& arow = Arow->RowMatrixRowMap();
  Epetra_MultiVector* x_row = new Epetra_MultiVector(View, arow, xPoint, 1);

  Epetra_MultiVector* b_row = new Epetra_MultiVector(View, arow, bPoint, 1);




  // Creates an epetra linear problem, contaning matrix
  // A, solution x and rhs b.
  Epetra_LinearProblem Problem(Arow, x_row, b_row);

  /*
   *   CheckInput mainly checks that the x vector has the same Map as OperatorDomain
   */
  int retn2 = Problem.CheckInput();
  if (retn2) {
     std::cerr << "Epetra_LinearProblem CheckInput routine returned an error" << retn2 << ". This is fatal" << std::endl;
     return(-1);
  }

   // Initializes the Amesos solver. This is the base class for
   // Amesos. It is a pure virtual class (hence objects of this
   // class cannot be allocated, and can exist only as pointers
   // or references).
   //
   Amesos_BaseSolver* Solver = 0;

   // Initializes the Factory. Factory is a function class (a
   // class that contains methods only, no data). Factory
   // will be used to create Amesos_BaseSolver derived objects.
   //
   Amesos Factory;

   // Specifies the solver. String ``SolverType'' can assume one
   // of the following values:
   // - Lapack
   // - Klu
   // - Umfpack
   // - Pardiso
   // - Taucs
   // - Superlu
   // - Superludist
   // - Mumps
   // - Dscpack
   //
   // HKM checked these out  - as working.
   //std::string solverType = "Klu";
   std::string solverType = jac->directSolverName_;
   //std::string SolverType = "Umfpack";
   //std::string SolverType = "Superludist";
   Solver = Factory.Create(solverType, Problem);
   // Factory.Create() returns 0 if the requested solver
   // is not available
   if (Solver == 0) {
     std::cerr << "Selected solver (" << solverType << ") is not available" << std::endl;
     std::cerr << "     available solvers: Umfpack(default) Klu, Superludist" << std::endl;
     return(-1);
   }

   // Calling solve to compute the solution. This calls the symbolic
   // factorization and the numeric factorization.
   int retn = Solver->Solve();
   if (retn != 0) {
     std::cerr << "Solver returned error " << retn << "\n" << std::endl;
     return(-1);
   }
   // Print out solver timings and get timings in parameter list.
   //Solver->PrintStatus();
   //Solver->PrintTiming();

   Teuchos::ParameterList TimingsList;
   Solver->GetTiming( TimingsList );

    // =========================================== //
    // E N D   O F   T H E   A M E S O S   P A R T //
    // =========================================== //

    // Computes ||Ax - b|| //

    double residual = 0.0;
    if ((PrintLevel_LinearSolver) > 1) {

      Epetra_Vector Ax(b->Map());
      Arow->Multiply(false, *x, Ax);
      Ax.Update(1.0, *b, -1.0);
      Ax.Norm2(&residual);

#ifdef HAVE_MPI
      const Epetra_Comm &cc = x->Comm();
      int mypid = cc.MyPID();
      if (mypid == 0) {
        std::cout << "Successful Amesos solution: ||b - A * x||_2 = " << residual << std::endl;
      }
#else
        std::cout << "Successful Amesos solution: ||b - A * x||_2 = " << residual << std::endl;
#endif
    }

    // delete Solver. Do this before calling MPI_Finalize() because
    // MPI calls can occur.
    delete Solver;

    delete x_row;
    delete b_row;

    if (residual > 1e-5) {
      return(1);
    }
    return 0;
}
//==================================================================================
}
//==================================================================================


