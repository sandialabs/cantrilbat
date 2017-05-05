/**
 * @file m1d_solvers.h
 *
 */

/*
 * Copywrite 2013 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef M1D_SOLVERS_H
#define M1D_SOLVERS_H

#include "m1d_defs.h"

#include <Epetra_Vector.h>
#include <Epetra_VbrMatrix.h>
#include <Epetra_VbrRowMatrix.h>
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
class LocalRowNodeVBRIndices;
class GlobalIndices;
class EpetraJac;

//! Construct the matrix
/*!
 *
 * @param A  Matrix to be solved for
 * @param v
 * @param b  rhs of the matrix problem
 */
void
fill_matrix(Epetra_VbrMatrix*& A,
            Epetra_Vector*& v,
            Epetra_Vector*& b,
            m1d::GlobalIndices *gi_ptr,
            m1d::LocalRowNodeVBRIndices *lrn_vbr_ptr);

//! Print level for the linear solver
/*!
 *   0 absolutely no printouts ever
 *   1 Print out only if a fatal error condition occurs
 *   2 One line of printouts can occur
 *   3 5 lines per call can occur
 *   4 20 lines per call can occur 
 *   5 One line per iteration can occur
 */
extern int PrintLevel_LinearSolver;

int solve_aztecoo(Epetra_VbrMatrix* const A, Epetra_Vector* const v, Epetra_Vector* const b, EpetraJac* const jac);

int solve_aztecoo(Epetra_VbrRowMatrix* const A, Epetra_Vector* const v, Epetra_Vector* const b,  EpetraJac* const jac);

//! Solve the system
/*!
 *   Note, I've tried making b const, and it doesn't work. The format is to assume that b may change during the solution process.
 *   The actually underlying software assumes that the matrix is an Epetra_RowMatrix and b is an Epetra_MultiVector. This should be
 *   a priority to switch to that format. (and figure out what the reasons are for the differentiation)
 */
int solve_amesos(Epetra_VbrMatrix* const A, Epetra_Vector* const v, Epetra_Vector* b, EpetraJac* const jac);

int solve_amesos(Epetra_VbrRowMatrix* const A, Epetra_Vector* const v, Epetra_Vector* const b,  EpetraJac* const jac);

}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
