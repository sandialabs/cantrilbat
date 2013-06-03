/**
 * @file tddiff.h
 *
 */

/*
 *  $Id: Cu2S_growth.h 506 2013-01-07 22:43:59Z hkmoffa $
 */

#include "cantera/base/mdp_allo.h"

#define HAVE_MPI

#include "Epetra_ConfigDefs.h"

#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "m1d_Comm.h"
#include "m1d_EpetraExtras.h"
#include "m1d_globals.h"

#include "m1d_DomainLayout.h"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_VbrMatrix.h>
#include <Teuchos_ParameterList.hpp>
#include <Ifpack.h>
#include <AztecOO.h>

#include <vector>
#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

namespace m1d
{
class ProblemResidEval;
class GlobalIndices;
class LocalRowNodeVBRIndices;
}

extern m1d::ProblemResidEval *PS_ptr;

Epetra_VbrMatrix*
alloc_VbrMatrix(MPI_Comm mpi_comm);


int
alloc_double_matrix_3d(double ****double_matrix_3d, int L, int M, int N);

void
ex_write_output_file(MPI_Comm mpi_comm,
                     Epetra_CrsMatrix*& A,
                     Epetra_Vector*& v,
                     Epetra_Vector*& b);



//!  Generate Domain Objects
/*!
 *s
 */
void
generateDomain1D(m1d::ProblemResidEval *ps);

