/**
 * @file m1d_EpectraExtras.cpp
 *
 */

/*
 *  $Id: m1d_EpetraExtras.cpp 560 2013-03-06 23:50:04Z hkmoffa $
 */
#include "m1d_EpetraExtras.h"

#include "m1d_globals.h"
#include "m1d_exception.h"

#include <stdio.h>
#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "m1d_globals.h"
#include "m1d_EpetraExtras.h"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_VbrMatrix.h>
#include <Epetra_Import.h>

namespace m1d
{
//=====================================================================================================================
Epetra_Comm *Comm_ptr = 0;

//=====================================================================================================================
Epetra_Vector *
gatherOn0(Epetra_Vector &distribV, Epetra_Comm *acomm_ptr)
{
  Epetra_Comm *comm_ptr = acomm_ptr;
  if (acomm_ptr == 0) {
    if (Comm_ptr) {
      comm_ptr = Comm_ptr;
    } else {
      throw m1d_Error("gatherOnO", "no comm pointer");
    }
  }
  int procID = comm_ptr->MyPID();

  // Elements are the number of block rows in the grid
  // HKM -> Check for numEqs > 1 compatibility
  int numMyElements = 0;
  int numGlobalElements = distribV.GlobalLength();
  if (procID == 0) {
    numMyElements = numGlobalElements;
  }

  /*
   * Create a map with all of the unknowns on processor 0
   */
  Epetra_Map *e0_map = new Epetra_Map(numGlobalElements, numMyElements, 0,
                                      *comm_ptr);

  /*
   *  Extract the map for the distributed vector from the object
   */
  const Epetra_BlockMap &distribMap = distribV.Map();

  /*
   *  Create an import object that describes the communication to bring
   *  all of the unknowns from the distributed object onto processor 0
   */
  Epetra_Import *import0 = new Epetra_Import(*e0_map, distribMap);

  /*
   *  Malloc and create the Epetra object to hold all of the data on proc 0
   */
  Epetra_Vector *e0 = new Epetra_Vector(*e0_map, true);

  /*
   *  Bring all of the data onto processor zero and store it in e0.
   */
  e0->Import(distribV, *import0, Insert, 0);

  delete import0;
  delete e0_map;

  return e0;
}
//=====================================================================================================================
void
printOn0(Epetra_Vector &distribV, Epetra_Comm *acomm_ptr)
{
  Epetra_Comm *comm_ptr = acomm_ptr;
  if (acomm_ptr == 0) {
    if (Comm_ptr) {
      comm_ptr = Comm_ptr;
    } else {
      throw m1d_Error("gatherOnO", "no comm pointer");
    }
  }
  int procID = comm_ptr->MyPID();
  Epetra_Vector *e0 = gatherOn0(distribV, acomm_ptr);
  int numGlobalElements = distribV.GlobalLength();
  if (procID == 0) {
    for (int i = 0; i < numGlobalElements; i++) {
      printf(" i = %d, distribV = %g\n", i, (*e0)[i]);
    }
  }

  delete e0;
}
//=====================================================================================================================
/*
 * This has been shown to only work if the element size is equal to one.
 * In order to fix this we need to first do a gatherOnAll for the element sizes across all processors!
 * This is not impossible, but delayed implementation.
 */
Epetra_Vector *
gatherOnAll(const Epetra_Vector &distribV, Epetra_Comm *acomm_ptr)
{
  Epetra_Comm *comm_ptr = acomm_ptr;
  if (acomm_ptr == 0) {
    if (Comm_ptr) {
      comm_ptr = Comm_ptr;
    } else {
      throw m1d_Error("gatherOnO", "no comm pointer");
    }
  }
  //int procID = comm_ptr->MyPID();
  /*
   * Assign all elements to all processors
   */
  int numMyElements = distribV.GlobalLength();
  int numGlobalElements = numMyElements;
  /*
   * Create a map with all of the unknowns on all processors
   */
  Epetra_Map *eAll_map = new Epetra_Map(numGlobalElements, numMyElements, 0,
                                      *comm_ptr);
  // HKM needs fixing.
  //new Epetra_BlockMap(-1, NumOwnedLcNodes, DATA_PTR(IndexGbNode_LcNode),
  //      DATA_PTR(NumEqns_LcNode), 0, *Comm_ptr_);
  //}
  /*
   *  Extract the map on this processor for the distributed vector from the object
   */
  const Epetra_BlockMap &distribMap = distribV.Map();
  /*
   *  Create an import object that describes the communication to bring
   *  all of the unknowns from the distributed object onto all of the processors
   */
  Epetra_Import *importAll = new Epetra_Import(*eAll_map, distribMap);
  /*
   *  Malloc and create the Epetra object to hold all of the data on all procs
   */
  Epetra_Vector *eAll = new Epetra_Vector(*eAll_map, true);

  /*
   *  Bring all of the data onto all processors and store it in e0.
   */
  eAll->Import(distribV, *importAll, Insert, 0);

  delete importAll;
  delete eAll_map;
  return eAll;
}
//=====================================================================================================================
/*
 * We assume here that global_node_V and distrib_node_V are internally consistent.
 */
void
gather_nodeV_OnAll(Epetra_Vector & global_node_V, const Epetra_Vector &distrib_node_V, Epetra_Comm *acomm_ptr)
{
  Epetra_Comm *comm_ptr = acomm_ptr;
  if (acomm_ptr == 0) {
    if (Comm_ptr) {
      comm_ptr = Comm_ptr;
    } else {
      throw m1d_Error("gather_nodeV_OnAll", "no comm pointer");
    }
  }
  /*
   * Assign all elements to all processors
   */
  int numMyElements = distrib_node_V.GlobalLength();
  //int numGlobalElements = numMyElements;
  if (global_node_V.MyLength() != numMyElements) {
    printf("(global length = %d) != (distrib length = %d) \n", global_node_V.MyLength(), numMyElements);
  }
  // AssertTrace(global_node_V.MyLength() == numMyElements);

  int ng = global_node_V.GlobalLength();
  int nl = global_node_V.MyLength();
  AssertTrace(ng == nl);
  /*
   * Create a map with all of the unknowns on all processors
   */
  //Epetra_Map *eAll_map = new Epetra_Map(numGlobalElements, numMyElements, 0,
    //                                  *comm_ptr);
  /*
   *  Extract the map on this processor for the distributed vector from the object
   */
  const Epetra_BlockMap &distribMap = distrib_node_V.Map();
  /*
   *  Create an import object that describes the communication to bring
   *  all of the unknowns from the distributed object onto all of the processors
   */
  Epetra_Import *importAll = new Epetra_Import(global_node_V.Map(), distribMap);

  /*
   *  Malloc and create the Epetra object to hold all of the data on all procs
   */
  //Epetra_Vector *eAll = new Epetra_Vector(*eAll_map, true);

  /*
   *  Bring all of the data onto all processors and store it in e0.
   */
  global_node_V.Import(distrib_node_V, *importAll, Insert, 0);

  delete importAll;
 // delete eAll_map;

  return;
}
//=====================================================================================================================
/*
 * We assume here that global_node_IV and distrib_node_IV are internally consistent.
 */
void
gather_nodeIntV_OnAll(Epetra_IntVector & global_node_IV, const Epetra_IntVector &distrib_node_IV, Epetra_Comm *acomm_ptr)
{
  Epetra_Comm *comm_ptr = acomm_ptr;
  if (acomm_ptr == 0) {
    if (Comm_ptr) {
      comm_ptr = Comm_ptr;
    } else {
      throw m1d_Error("gather_nodeV_OnAll", "no comm pointer");
    }
  }
  /*
   * Assign all elements to all processors
   */ 
  int numMyElements = distrib_node_IV.GlobalLength();
  //int numGlobalElements = numMyElements;
  if (global_node_IV.MyLength() != numMyElements) {
    printf("(global length = %d) != (distrib length = %d) \n", global_node_IV.MyLength(), numMyElements);
  }
  AssertTrace(global_node_IV.MyLength() == numMyElements);
  
  int ng = global_node_IV.GlobalLength();
  int nl = global_node_IV.MyLength();
  AssertTrace(ng == nl);
  /*
   * Create a map with all of the unknowns on all processors
   */
  //Epetra_Map *eAll_map = new Epetra_Map(numGlobalElements, numMyElements, 0,
    //                                  *comm_ptr);
  /*
   *  Extract the map on this processor for the distributed vector from the object
   */
  const Epetra_BlockMap &distribMap = distrib_node_IV.Map();
  /*
   *  Create an import object that describes the communication to bring
   *  all of the unknowns from the distributed object onto all of the processors
   */
  Epetra_Import *importAll = new Epetra_Import(global_node_IV.Map(), distribMap);
  
  /*
   *  Malloc and create the Epetra object to hold all of the data on all procs
   */
  //Epetra_Vector *eAll = new Epetra_Vector(*eAll_map, true);
  
  /*
   *  Bring all of the data onto all processors and store it in e0.
   */
  global_node_IV.Import(distrib_node_IV, *importAll, Insert, 0);
  
  delete importAll;
 // delete eAll_map;

  return;
}

//=====================================================================================================================
}
//=====================================================================================================================
