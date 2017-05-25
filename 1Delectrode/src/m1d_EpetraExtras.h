/**
 * @file m1d_EpetraExtras.h Utility routines that carry out distributed Epetra operations
 *
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
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! Take a distributed vector and gather it all on processor 0
/*!
 *  This returns an Epetra_Object whose data is all located on processor 0. It mallocs
 *  the Epetra_Vector and therefore, it must be freed by the calling program.
 *
 *  @param[in]               distribV            Epetra vector that is distributed on all processors. It may be ghosted or
 *                                               unghosted. It does not matter.
 *  @param[in]               comm_ptr            Pointer to the Epetra_comm object
 *
 *  @return                                      Returns a pointer to the malloced vector
 */
Epetra_Vector* gatherOn0(const Epetra_Vector& distribV, const Epetra_Comm* const comm_ptr = nullptr);

//==================================================================================================================================
//! Print a distributed vector from processor 0
/*!
 *  We call gatherOn0 to collect the vector on processor 0. Then we just do a printf
 *
 *  @param[in]               distribV            Epetra vector that is distributed on all processors. It may be ghosted or
 *                                               unghosted. It does not matter.
 *  @param[in]               comm_ptr            Pointer to the Epetra_comm object
 *
 */
void printOn0(const Epetra_Vector& distribV, const Epetra_Comm* const comm_ptr = nullptr);

//==================================================================================================================================
//! Take a distributed vector and gather it all on all of the processors
/*!
 *  This is done for the nodal positions, for example.
 *
 *  This has been shown to only work if the element size is equal to one.
 *  In order to fix this we need to first do a gatherOnAll for the element sizes across all processors!
 *  This is not impossible, but delayed implementation.
 *
 *  @param[in]               distribV            Epetra vector that is distributed on all processors. It may be ghosted or
 *                                               unghosted. It does not matter. Only the owned values are used during the gather
 *                                               from each processor.
 *
 *  @param[in]               comm_ptr            Pointer to the Epetra_comm object
 *
 *  @return                                      Returns a pointer to the malloced Epetra_Vector
 */
Epetra_Vector* gatherOnAll(const Epetra_Vector& distribV, const Epetra_Comm* const comm_ptr = nullptr);

//==================================================================================================================================
//! Gather all of a distributed vector of element size one onto an existing Global_All_ALl vector
/*!
 *  @param[out]              global_node_V       Epetra_Vector that is an all-all vector (all indecises on all processors)
 *
 *  @param[in]               distrib_node_V      Epetra vector that is distributed on all processors. It may be ghosted or
 *                                               unghosted. It does not matter. Only the owned values are used during the gather
 *                                               from each processor.
 *
 *  @param[in]               comm_ptr            Pointer to the Epetra_comm object. Defaults to nullptr
 */
void gather_nodeV_OnAll(Epetra_Vector& global_node_V, const Epetra_Vector& distrib_node_V,
                        const Epetra_Comm* const comm_ptr=nullptr);

//==================================================================================================================================
//! Gather all of a distributed vector of ints of element size one onto an existing Global_All_ALl int vector
/*!
 *  @param[out]              global_node_IV      Epetra_IntVector that is an all-all vector (all indecises on all processors)
 *
 *  @param[in]               distrib_node_IV     Epetra IntVector that is distributed on all processors. It may be ghosted or
 *                                               unghosted. It does not matter. Only the owned values are used during the gather
 *                                               from each processor.
 *
 *  @param[in]               comm_ptr            Pointer to the Epetra_comm object. Defaults to nullptr
 */
void gather_nodeIntV_OnAll(Epetra_IntVector& global_node_IV, const Epetra_IntVector& distrib_node_IV,
                           const Epetra_Comm* const comm_ptr=nullptr);

//==================================================================================================================================
//! Create a new view of the original vector using potentially a different map
/*!
 *  @param[in]               orig                Original Epetra_Vector
 *
 *  @param[in]               nmap                New map of the vector
 *
 *  @return                                      Returns a malloced Epetra_Vector containing a new view of the original information
 */
Epetra_Vector* new_EpetraVectorView(const Epetra_Vector& orig, const Epetra_BlockMap& nmap);

//==================================================================================================================================
//! Scatter to All from Processor zero.
/*!
 *  @param[in]               globalSoln          Epetra_Vector that is a  global_all_all vector
 *
 *  @param[out]              distribV            Distributed vector
 *
 *  @param[in]               comm_ptr            Pointer to the Epetra_comm object. Defaults to nullptr
 */
void scatterToAllFrom0(const Epetra_Vector& globalSoln, Epetra_Vector& distribV, const Epetra_Comm* const comm_ptr = nullptr);

//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
