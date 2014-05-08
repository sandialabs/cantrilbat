/**
 * @file tddiff.h
 *
 */

/*
 *  $Id: tddiff.h 504 2013-01-07 22:32:48Z hkmoffa $
 */



#include "Epetra_ConfigDefs.h"

#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "m1d_Comm.h"
#include "m1d_globals.h"



namespace m1d
{
class ProblemResidEval;
class GlobalIndices;
class LocalRowNodeVBRIndices;
}

extern m1d::ProblemResidEval *PS_ptr;

Epetra_VbrMatrix* alloc_VbrMatrix(M1D_MPI_Comm mpi_comm);

int
alloc_double_matrix_3d(double ****double_matrix_3d, int L, int M, int N);

void
ex_write_output_file(M1D_MPI_Comm mpi_comm,
                     Epetra_CrsMatrix*& A,
                     Epetra_Vector*& v,
                     Epetra_Vector*& b);



//!  Generate Domain Objects
/*!
 *s
 */
void
generateDomain1D(m1d::ProblemResidEval *ps);

