/**
 * @file LiKCl_PorousBat.h
 *
 */

/*
 *  $Id: LiIon_PorousBat.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "Epetra_ConfigDefs.h"

#include "m1d_Comm.h"
#include "m1d_globals.h"
#include "m1d_BatteryResidEval.h"

#include "m1d_ProblemStatementCell.h"


namespace m1d
{
class ProblemResidEval;
class GlobalIndices;
class LocalRowNodeVBRIndices;
}

extern m1d::BatteryResidEval* PS_ptr;

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
 *s
 */
void
generateDomain1D(m1d::BatteryResidEval* ps);

