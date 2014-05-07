/**
 * @file LiKCl_PorousBat.h
 *
 */
/*
 * Copywrite 2013 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "Epetra_ConfigDefs.h"

#include "m1d_Comm.h"
#include "m1d_globals.h"
#include "m1d_BatteryResidEval.h"

//! Forward declarations of the classes within the m1d namespacer
namespace m1d
{
class ProblemResidEval;
class GlobalIndices;
class LocalRowNodeVBRIndices;
class ProblemStatementCell;
}

//! Pointer to the Structure which is responsible for forming the residual
/*!
 *   There is one and only one structure. However, this may change in the future
 */
extern m1d::BatteryResidEval* PS_ptr;

//! Global Problem input structure
/*!
 *   This contains the input data for the problem.
 *   We've made it a global structure, as there is one and only one instance of the structure
 */
extern m1d::ProblemStatementCell PSinput;


Epetra_VbrMatrix*
alloc_VbrMatrix(M1D_MPI_Comm mpi_comm);


int
alloc_double_matrix_3d(double**** double_matrix_3d, int L, int M, int N);

void
ex_write_output_file(M1D_MPI_Comm mpi_comm,
                     Epetra_CrsMatrix*& A,
                     Epetra_Vector*& v,
                     Epetra_Vector*& b);

//!  Generate Domain Objects
/*!
 *
 */
void
generateDomain1D(m1d::BatteryResidEval* ps);

