/**
 * @file m1d_GlobalIndices
 *  Global structure containing global information about macrodomains, nodes, and equations.
 */

/*
 *  $Id: m1d_GlobalIndices.h 361 2012-08-21 00:39:02Z hkmoffa $
 */

#ifndef M1D_GLOBALINDICES_H
#define M1D_GLOBALINDICES_H

#include <vector>

#include "m1d_defs.h"

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"

class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Comm;

namespace m1d
{

class NodalVars;
class DomainLayout;

//! Global indices class is the same on all processors
/*!
 *  This is a global structure, containing global information.
 *  All information in this structure, except for my_procID, numLnNodes,
 *  and numLnEqns, is the same on all processors.
 *
 *  How the global indexing works
 *  --------------------------------------------------------
 *
 *  Global indexing allows for different numbers of equations per node. Simply,
 *  all equations at a single node are contiguous. Global nodes are numbered
 *  sequentially from left to right in the domain.
 *  If there are global equations, they are located at a bottom (i.e. right) 
 *  "virtual" global node of the domain.
 *
 *   GbEqn = IndexStartGbEqns_GbNode(GbNode) + local eqn #
 *
 *  When there is only one processor in the problem, the local node numbering
 *  will be exactly the same as the global node numbering.
 *
 */
class GlobalIndices {
public:

  //! Constructor
  GlobalIndices(Epetra_Comm *comm_ptr);

  //! Destructor
  ~GlobalIndices();

  //! Copy constructor
  /*!
   * @param r Object to be copied
   */
  GlobalIndices(const GlobalIndices &r);

  //! Assignment operator
  /*!
   *  @param r   Object to be copied
   */
  GlobalIndices &
  operator=(const GlobalIndices &r);

  //! Initialize the number of global nodes
  /*!
   * Size all arrays based on the total number of global nodes
   */
  void
  init(DomainLayout *dl_ptr);

  //! Calculate the total number of unknowns per node
  int
  discoverNumEqnsPerNode();

  void
  procDivide();

  void
  initNodeMaps();

  void
  initBlockNodeMaps(int *numEqns_LcNode);

  //  Epetra_Vector *
  //createDistribEqnVector();

  void
  InitMesh();

  //! This utility function will return the global node number
  //! given the global equation number
  /*!
   * @param   rowEqnNum returns the node equation number
   * @return  Returns the global node number. Will return -1 if there
   *          is a problem.
   */
  int
  GbEqnToGbNode(const int GbEqnNum, int & rowEqnNum) const;

  //! Take a distributed Epetra_Vector vector and make it into a globally-all distributed
  //! vector of the node positions, and then feed it to the nodal Values object
  /*!
   *
   * @param Xpos_LcNode_p
   */
  void
  updateGlobalPositions(Epetra_Vector *Xpos_LcNode_p);

  //! Communications object
  Epetra_Comm *Comm_ptr_;

  //! Number of processors
  int NumProc;

  //! My procID
  int MyProcID;

  //! Number of global nodes
  int NumGbNodes;

  //! Number of global equations
  int NumGbEqns;

  //! Index of the starting global node number owned by each processor.
  /*!
   * Length = number of processors.
   */
  std::vector<int> IndexStartGbNode_Proc;

  //! Index of the starting owned global equation number owned by each processor
  /*!
   * Length = number of processors.
   */
  std::vector<int> IndexStartGbEqns_Proc;

  //! Number of local nodes owned by each processor
  /*!
   * Length = number of processors.
   */
  std::vector<int> NumOwnedLcNodes_Proc;

  //! Number of Local Equations owned by each processor
  /*!
   * Length = number of processors.
   */
  std::vector<int> NumOwnedLcEqns_Proc;

  //! Number of equations at each global node
  /*!
   *  Length  = number of global nodes
   */
  std::vector<int> NumEqns_GbNode;

  //! The starting global equation index at each global node
  /*!
   *   Length = number of global nodes.
   */
  std::vector<int> IndexStartGbEqns_GbNode;

  //! Number of local nodes owned by this processor
  int NumOwnedLcNodes;

  //! Number of local equations owned by this processor
  int NumOwnedLcEqns;

  //! Epetra_Map object containing the mapping of the equations to processors.
  /*!
   *  This describes the distribution of the equations onto the processors.
   *  Input to create this map involves specifying the number of local equations
   *  that each processor owns. This map then describes the global indexing that
   *  results.
   *
   *  It sets up the local processor indexing vs the global indexing
   *  for the solution vector.
   *
   *  This is a distributed object. It's used in the creation of all distributed
   *  vectors of the equation system.
   *
   * HKM -> WARNING THIS IS NOT A BLOCK MAP!!!
   */
  //  Epetra_Map *GbEqnstoOwnedLcEqnsMap;

  //! Epetra_Map object containing the mapping of the nodes to processors.
  /*!
   *  This describes the distribution of the nodes onto the processors.
   *  Input to create this map involves specifying the number of local equations
   *  that each processor owns. This map then describes the global indexing that
   *  results.
   *
   *  It sets up the local processor indexing vs the global indexing
   *  for the node vector.
   *
   *  This is a distributed object. It's used in the creation of all distributed
   *  vectors over global nodes within the module.
   */
  Epetra_Map *GbNodetoOwnedLcNodeMap;

  //! Epetra_Map object containing the mapping of the block node eqns to processors.
  /*!
   *  This describes the distribution of the block node eqns onto the processors.
   *
   *  Input to create this map involves specifying the number of local equations
   *  that each processor owns. This map then describes the global indexing that
   *  results.
   *
   *  It sets up the local processor indexing vs the global indexing
   *  for the node vector.
   *
   *  This is a distributed object. It's used in the creation of all distributed
   *  vectors over global nodes within the module.
   */
  Epetra_BlockMap *GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap;

  //! Global Vector of nodal variable descriptions.
  /*!
   * NodalVars object contains a complete description of the variables at
   * a global node. 
   * length = number of global nodes.
   */
  std::vector<NodalVars *> NodalVars_GbNode;

  //! Positions at the global nodes
  /*!
   *   We keep a global vector of these on each node
   */
  Epetra_Vector *XNodePos_GbNode;

  //! Epetra_Map object containing the mapping of all equations as if they were on all processors.
  /*!
   *  This describes the result of a gather_to_all operation involving equations
   */
  Epetra_BlockMap *GbEqnstoAllMap;

  //! Vector containing the entire solution on all processors.
  Epetra_Vector *SolnAll;

  Epetra_IntVector *SolnIntAll;

  //! Vector containing the entire solution on all processors.
  Epetra_Vector *SolnDotAll;

  //! Domain Layout object
  /*!
   *
   */
  DomainLayout *DL_ptr_;
};

//=====================================================================================
//! Each processor contains a pointer to one global instance of this class.
//extern GlobalIndices *GI_ptr;
//=====================================================================================

}

#endif
