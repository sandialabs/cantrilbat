/**
 * @file LiKCl_PorousBat.h
 *
 */

/*
 *  $Id: LiKCl_PorousBat.h 5 2012-02-23 21:34:18Z hkmoffa $
 */


#include "config.h"
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
#include "m1d_BatteryResidEval.h"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
//#include "Teuchos_ParameterList.hpp"
#include "Ifpack.h"
#include "AztecOO.h"

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

extern m1d::BatteryResidEval *PS_ptr;

Epetra_VbrMatrix*
alloc_VbrMatrix(M1D_MPI_Comm mpi_comm);


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
generateDomain1D(m1d::BatteryResidEval *ps);

