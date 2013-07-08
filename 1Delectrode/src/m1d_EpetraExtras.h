/**
 * @file m1d_EpectraExtras.h
 *
 */

/*
 *   $Id: m1d_EpetraExtras.h 504 2013-01-07 22:32:48Z hkmoffa $
 */

#ifndef _M1D_EPECTRAEXTRAS_H
#define _M1D_EPECTRAEXTRAS_H

#include "Epetra_ConfigDefs.h"

#ifdef HAVE_MPI
#include <mpi.h>
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_VbrMatrix.h>
#include <Epetra_Import.h>

namespace m1d
{
//! Take a distributed vector and gather it all on processor 0
/*!
 *   This returns an Epetra_Object whose data is all located on processor 0. It mallocs
 *   the Epetra_Vector and therefore, it must be freed by the calling program.
 *
 *  @param distribV     Epetra vector that is distributed on all processors. It may be ghosted or
 *                      unghosted. It does not matter.
 *  @param comm_ptr     Pointer to the Epetra_comm object
 *
 *  @return Returns a pointer to the malloced vector
 */
Epetra_Vector *
gatherOn0(Epetra_Vector &distribV, Epetra_Comm *comm_ptr = 0);

void
printOn0(Epetra_Vector &distribV, Epetra_Comm *acomm_ptr = 0);

//! Take a distributed vector and gather it all on all processors
/*!
 *
 * @param comm         Comm object
 * @param distribV     Epetra vector that is distributed on all processors
 */
Epetra_Vector *
gatherOnAll(const Epetra_Vector &distribV, Epetra_Comm *comm_ptr = 0);

//! 
void
gather_nodeV_OnAll(Epetra_Vector & global_node_V, const Epetra_Vector &distrib_node_V, Epetra_Comm *acomm_ptr=0);

void
gather_nodeIntV_OnAll(Epetra_IntVector & global_node_IV, const Epetra_IntVector &distrib_node_IV, Epetra_Comm *acomm_ptr=0);

//! Create a new view of the original vector using potentially a different map
Epetra_Vector * new_EpetraVectorView(const Epetra_Vector& orig, const Epetra_BlockMap& nmap);



}
#endif
