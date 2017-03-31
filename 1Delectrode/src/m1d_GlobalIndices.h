/**
 * @file m1d_GlobalIndices.h
 *   Definitions for the global structure containing global information about macrodomains, nodes, and equations.
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef M1D_GLOBALINDICES_H
#define M1D_GLOBALINDICES_H

#include <vector>

#include "m1d_defs.h"

#include "Epetra_IntVector.h"

class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Comm;
class Epetra_Vector;

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{

class NodalVars;
class DomainLayout;
//==================================================================================================================================
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
 *  Note on Ordinals
 *  ----------------------------------
 *
 *   Epetra either uses "int" or "long long" as its ordinals. It does not used "unsigned long = size_t"
 *   The choice is between int and long log. I chose int.
 *   Therefore all calls to Epetra should use int ordinal types.
 *
 */
class GlobalIndices {
public:

  //! Constructor
  /*!
   *  @param[in]             comm_ptr            Epetra_Comm pointer
   */
  GlobalIndices(Epetra_Comm *comm_ptr);

  //! Destructor
  ~GlobalIndices();

  //! Copy constructor
  /*!
   *  @param[in]             r                   Object to be copied
   */
  GlobalIndices(const GlobalIndices &r);

  //! Assignment operator
  /*!
   *  @param[in]             r                   Object to be copied
   *
   *  @return                                    Returns a reference to the global object
   */
  GlobalIndices& operator=(const GlobalIndices &r);

  //! Initialize the number of global nodes by reading in the domain layout object
  /*!
   *  Size all arrays based on the total number of global nodes
   *
   *  Creates:
   *     NumGbNodes          Number of global nodes
   *     NumEqns_GbNode      global vector of number of equations on each node
   *     NodalVars_GbNode    global Vector of NodeVars objects for each global node
   *     XNodePos_GbNode     Distributed Epetra vector of doubles with each global index corresponding to the local index on 
   *                         each node.
   *
   *  @param[in]             dl_ptr              Pointer to the domain layout object. 
   *                                             The total number of global nodes has already been calculated here.
   */
  void init(DomainLayout *dl_ptr);

  //! Calculate the total number of unknowns per node
  /*!
   *  We find out which domains are located at each node, and we find out what equations are located 
   *  at each node.  We calculate a couple of key quantities here, and store them in this global structure.
   *
   *    NumEqns_GbNode[iGbNode]                Number of equations at the global node iGbNode
   *    IndexStartGbEqns_GbNode[iGbNode]       Starting index for the equations at iGbNode in the vector of global equations
   *    NumGbEqns                              Total number of global equations
   *
   *   We store the starting global equation index back into the NodeVars object
   *
   *  @return                                    Returns the total number of global equations in the problem
   */
  size_t discoverNumEqnsPerNode();

  //! Divide up the problem amongst the processors
  /*!
   *  Here we divide the nodes and equations amongst the processors.
   *  We use a simple linear division
   *  The following variables are formed here:
   *
   *   IndexStartGbNode_Proc[NumProc]       Index of the first global node on the ith processor
   *   IndexStartGbEqns_Proc[NumProc]       Index of the first global equation on the ith processor
   *   NumOwnedLcNodes_Proc[NumProc]        Number of nodes owned by the ith processor
   *   NumOwnedLcEqn[NumProc]               Number of owned equation by the ith processor
   * 
   */
  void procDivide();

  //! Formulate an Epetra map that describes the layout of owned nodes indecises on the processor
  /*!
   *  The following epetra map is formed:
   *
   *   GbNodetoOwnedLcNodeMap              Global length = NumGbNodes
   *                                       Local length = NumOwnedLcNodes_Proc[i]
   */
  void initNodeMaps();

  //! Initialize BlockMaps that will describe the equation layouts on the local processors
  /*!
   *   The following BlockMaps are initialized
   *      GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap  BlockMap of distributed node indecises with the number of equations per node
   *                                                   assigned to each block
   *
   *      GbEqnstoAllMap                               BlockMap with the full global block equation structure replicated on each processor
   *
   *   GlobalAll solution vectors are allocated on each processor
   *       SolnAll, SolnDotAll,  SolnIntAll 
   *
   */ 
  void initBlockNodeMaps();

  //! Initialize the position of the nodes of the mesh
  /*!
   *  This is done by calling the domain layout object
   */
  void InitMesh();

  //! This utility function will return the global node number given the global equation number
  /*!
   * @param   rowEqnNum returns the node equation number
   * @return  Returns the global node number. Will return -1 if there
   *          is a problem.
   */
  int GbEqnToGbNode(const int GbEqnNum, int & rowEqnNum) const;

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

  //! Number of global nodes as a size_t
  size_t NumGbNodes_s;

  //! Number of global equations
  int NumGbEqns;

  //! Number of global equations as a size_t
  int NumGbEqns_s;

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
   *  Length: number of processors.
   */
  std::vector<int> NumOwnedLcNodes_Proc;

  //! Total Number of Local Equations owned by each processor
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
   *
   * length = number of global nodes. Note this is repeated for all processors.
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
//==================================================================================================================================
//=====================================================================================
//! Each processor contains a pointer to one global instance of this class.
//extern GlobalIndices *GI_ptr;
//=====================================================================================

}

//----------------------------------------------------------------------------------------------------------------------------------
#endif
