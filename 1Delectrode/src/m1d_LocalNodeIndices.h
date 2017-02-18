/**
 * @file m1d_LocalNodeIndices.h
 *
 */

#ifndef M1D_LOCALNODEINDICES_H
#define M1D_LOCALNODEINDICES_H

#include "m1d_defs.h"

#include <vector>

class Epetra_Comm;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_VbrMatrix;
class Epetra_MapColoring;
class Epetra_IntVector;
class Epetra_Vector;
class Epetra_Import;

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{

class NodalVars;
class DomainLayout;
class GlobalIndices;

//==================================================================================================================================
//! These indices have to do with accessing the local nodes on the processor
/*!
 *   All numbers are in local row node format. Local row node format is defined as the following. All owned nodes
 *   come first. They are ordered in terms of increasing global node number.
 *   Then the right ghost node is listed.  Then the left ghost node is listed
 *   Then, the "globally-all-connected node is listed, if available. All of this structure spans domains
 *   starting from left to right. Owned nodes are defined irrespective of the domain structure.
 *
 *   For example if there are 5 processors, numbered from 0 to 4, with
 *   10 nodes per processor, with one globally-all-connected" node. Then,
 *   the following would be the layout for proc 1:
 *
 *      Local Row Node Index      GbNodeIndex         LcNodeIsExt
 *       0                            10                false
 *       . . .                        12-18             false
 *       NumOwnedLRNodes - 1=9        19                false
 *       10                           20                true
 *       11                            9                true
 *       12                           50                true   (global node)
 *
 *   There are 13 local nodes in all. Three of them are externals.
 *
 *    For last processor on the mesh the globally-all-connected node is slightly different.
 *    Here is the layout. The last processor owns the globally-all-connected node
 *    so its index occurs before
 *
 *      Local Row Node Index      GbNodeIndex         LcNodeIsExt
 *       0                            40                false
 *       . . .                        42-48             false
 *       NumOwnedLRNodes - 1=9        49                false
 *       10                           50                false   (global node)
 *       11                           39                true
 *
 *   The key variable for this structure is NumLcNodes, the number of local nodes defined for this processor.
 *   Another key variable is NumOwnedLcNodes. This is the number of owned local nodes defined on this processor.
 *   The difference is that NumLcNodes contains ghost nodes which this processor knows about, but doesn't own.
 *
 *   Also, MyProcID in this structure is the processor number that owns this part of the mesh. In a given problem,
 *   there are numProcs number of these structures in the problem.
 *
 */
class LocalNodeIndices {
public:

  //! Constructor
  LocalNodeIndices(Epetra_Comm *comm_ptr, GlobalIndices *gi_ptr);

  //! Destructor
  ~LocalNodeIndices();

  //! Copy constructor
  /*!
   * @param r  Object to be copied
   */
  LocalNodeIndices(const LocalNodeIndices &r);

  //! Assignment operator
  /*!
   * @param r  Object to be copied
   * @return   Return the current object
   */
  LocalNodeIndices & operator=(const LocalNodeIndices &r);

  //! Initialize the sizes and the contents of the member arrays
  //! and determine the nodes needed by the current processor to
  //! do its calculations
  /*!
   *  In this routine, we figure out the external nodes that are
   *  needed to handle the calculations on the node. We then create
   *  the Epetra maps that store that information.
   *
   *   This also determines the matrix stencil.
   */
  void determineLcNodeMaps(DomainLayout *dl_ptr);

  //! Initialize the LcEqn arrays, including the color array
  /*!
   *  In this routine, we figure out the color array for equations.
   */
  void determineLcEqnMaps();

  //! Initialize the Epetra_Map objects associated with specifying
  //! the local nodes on the processor, including whether they are owned or not
  /*!
   *    GbNodetoLcNodeColMap   a map of the local nodes and ghost nodes. The global ids are the global node numbers.
   *                           This map will included ghosted nodes.
   *                            This map is not a 1 to 1 map. 
   *   
   *    GbNodetoOwnedLcNodeMap A map of the local nodes. The global ids are the global node numbers.
   *                           This map will not include ghosted nodes.
   *                           This map is a 1 to 1 map. 
   */
  void initLcNodeMaps();

  //! Initialize the Epetra_Map objects associated with specifying
  //! the local block equations on the processor, including whether they are owned or not
  /*!
   *   GbBlockNodeEqnstoLcBlockNodeEqnsColMap
   *                           Create a block map of the local eqns and ghost eqns.
   *                           The global element ids are the global node numbers.
   *                           The number of rows in the block are the number of equations defined at each node.
   *
   *   GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap
   *                           Create a map of the local equations only. The global ids are the global node numbers.
   *                           This map will not include ghosted nodes, and therefore will be 1 to 1.
   */
  void initLcBlockNodeMaps();

  //! Construct a coloring map for the LcNodes on this processor.
  /*!
   *   This map includes the external nodes that are defined on the processor
   *   as well.
   *
   *   This coloring map will be used to calculate the jacobian. It's based on
   *   an assumption about the matrix stencil. Every node that is adjacent to
   *   another node or is involved with its residual evaluation receives
   *   a different color. Right now we assume the matrix stencil is strictly
   *   an adjacent node operation.   If this assumption is not
   *   correct, this is the location to change the algorithm.
   *
   *   This routine may be called after the LcNodes map is created.
   */
  void makeNodeColors();

  //! Construct a coloring map for the locally defined equations on this processor.
  /*!
   *   Equations which are owned or not owned
   *
   *   This coloring map will be used to calculate the jacobian. It's based on
   *   an assumption about the matrix stencil. Every node that is adjacent to
   *   another node or is involved with its residual evaluation receives
   *   a different color. Right now we assume the matrix stencil is strictly
   *   an adjacent node operation.   If this assumption is not
   *   correct, this is the location to change the algorithm.
   *
   *   This routine may be called after the LcNodes map is created.
   */
  void makeEqnColors();

  //! Global node to local node mapping
  /*!
   * Given a global node, this function will return the local node value.
   * If the global node is not on this processor, then this function returns -1.
   * The local node may or may not be owned by this processor.
   *
   * @param gbNode  global node
   * @return returns the local node value
   */
  int GbNodeToLcNode(const int gbNode) const;

  //! Global eqn to local eqn mapping
  /*!
   * Given a global eqn, this function will return the local eqn value.
   * If the global eqn is not on this processor, then this function returns -1.
   * The local eqn may or may not be owned by this processor.
   *
   * @param gbEqn  global node
   * @return returns the local eqn value
   */
  int GbEqnToLcEqn(const int gbEqn) const;

  void InitializeLocalNodePositions();

  //! Update the nodal variables structure
  /*!
   *   This routine will update the position and local node number
   */
  void UpdateNodalVarsPositions();

  //!  Extract the positions from the solution vector and propagate them into all other structures
  /*!
   *  The following member data are updated:
   *      this->Xpos_LcNode_p[];
   *      nv->XNodePos
   *
   * @param soln  Solution vector which contains displacements as one of the solution components
   */
  void ExtractPositionsFromSolution(const Epetra_Vector * const soln_p);

  //!  We set initial conditions here that make sense from a global perspective.
  //!  This should be done as a starting point. If there are better answers, it should be overridden.
  /*!
   *   Current activities:
   *
   *     Set Displacement_Axial unknowns to 0.0
   *     Assert that nodal value of x0NodePos is equal to the local vector of *Xpos_LcNode_p
   */
  void setInitialConditions(const bool doTimeDependentResid, Epetra_Vector *soln, Epetra_Vector *solnDot, 
                            const double t, const double delta_t);

  //! Generate the nodal variables structure
  /*!
   *   This routine will update the pointer to the NodalVars structure
   *   and it will update the positions
   */
  void GenerateNodalVars();

  //! Generate the importers
  /*!
   *   Data members calculated
   *
   *        Importer_GhostEqns 
   *        Importer_NodalValues
   */
  void makeImporters();

  //! Calculate the number of equations at each node
  /*!
   *  The number of equations per node and the total # of equations is determined here.
   *  The values are calculated from the NodalVars Data structure.
   *
   *  Data members calculated:
   *        NumEqns_LcNode[]
   *        IndexLcEqns_LcNode[]
   *        NumLcOwnedEqns
   *        NumLcEqns
   */
  int UpdateEqnCount();

  //! Generate Equation mapping vectors
  /*!
   *  The values are calculated from the NodalVars Data structure.
   *
   *  This must be called after the total number of local equations on a processor, NumLcEqns, is found out.
   *  Data members calculated:
   *        IndexGbEqns_LcEqns[]
   */
  void generateEqnMapping();

  void updateGhostEqns(Epetra_Vector * const solnV, const Epetra_Vector * const srcV);
  /* ---------------------------------------------------------- */

  //! local copy of the Epetra_Comm ptr object
  Epetra_Comm *Comm_ptr_;

  //! My processor ID
  int MyProcID;

  //! This is the number of local nodes, including external nodes
  int NumLcNodes;

  //! Number of locally owned nodes
  int NumOwnedLcNodes;

  //! Number of external nodes
  int NumExtNodes;

  //! Number of local row nodes in the matrix
  /*!
   *  This is usually equal to NumOwnedLcNodes
   */
  int NumLcRowNodes;

  //! Number of equations on this processor whether they are
  //! owned by this processor or not.
  int NumLcEqns;

  //! Number of owned equations on this processor
  int NumLcOwnedEqns;

  //! Identity of the globally-connected node
  /*!
   *   This node is connected to all other nodes. It is the location
   *   where global equations should be located. It does not have a spatial
   *   position. If the value is -1, there is no globally connected node.
   *   If present this global node is the last node defined in the problem.
   *   The last processor will own this node.
   */
  int GCNIndexGbNode;

  //! Identity of the globally-connected node in local node number
  /*!
   *   This node is connected to all other nodes. It is the location
   *   where global equations are located. It does not have a spatial
   *   position. If the value is -1, there is no globally connected node.
   */
  int GCNIndexLcNode;

  //! Global index of each of the local nodes
  /*!
   *  Length = NumLcNodes
   */
  std::vector<int> IndexGbNode_LcNode;

  //! Number of equations at each of the local nodes
  /*!
   *  Length = NumLcNodes
   */
  std::vector<int> NumEqns_LcNode;

  //! Starting local equation index for each local node
  /*!
   *  Provides a map between the local node index and the starting position
   *  for the equations corresponding to that node in the local solutions
   *  vector.
   *  
   *         Length = NumLcNodes
   */
  std::vector<int> IndexLcEqns_LcNode;

  //! Mapping between the local equation index and the global equation index
  /*!
   * Length = NumLcEqns
   */
  std::vector<int> IndexGbEqns_LcEqns;

  //! Vector indicating whether the current local Node is owned by
  //! the current processor
  /*!
   *  Boolean indicating whether the current local node is not
   *  owned by this processor.
   *
   *  Length = NumLcNodes
   */
  std::vector<bool> IsExternal_LcNode;

  //! Vector containing the ID of the node to the left of the
  //! current local node
  std::vector<int> IDLeftLcNode_LcNode;

  //! Vector containing the ID of the node to the right of the
  //! current local node
  /*!
   * *
   *  Length = NumLcNodes
   */
  std::vector<int> IDRightLcNode_LcNode;

  //! Vector containing the Global ID of the node to the left of the
  //! current local node
  /*!
   *
   *  Length = NumLcNodes
   */
  std::vector<int> IDLeftGbNode_LcNode;

  //! Vector containing the Global ID of the node to the right of the
  //! current local node
  /*
   *  Length = NumLcNodes
   */
  std::vector<int> IDRightGbNode_LcNode;

  //! This is the node with the highest axial z position in the processor
  int RightLcNode;

  //! This is the node with the lowest axial z position in the processor
  int LeftLcNode;

  //! Epetra_Map object containing the mapping of the nodes corresponding
  //! to column row nodes needed by each processor
  /*!
   *  This describes the distribution of the column nodes
   *  needed by this processor for each processor. This is not a
   *  1 to 1 map.
   *
   *  This is a distributed object. It's used to set up the columns of the
   *  matrix owned by this processor. It's also used to describe the
   *  all of the local nodes defined on a processor, whether they are owned
   *  or not
   *
   *  This is called the overlap map for the nodes. Might switch
   *  to this nomenclature in the future.
   */
  Epetra_Map *GbNodetoLcNodeColMap;

 //! Epetra_Map object containing the mapping of the nodes corresponding
 //! to row nodes needed by each processor
 /*!
  *  This describes the distribution of the owned nodes
  *  needed by this processor. This is a 1 to 1 map.
  *
  *  This is a distributed object. It's used to set up the rows of the
  *  node-only matrix owned by this processor. It's also used to describe the
  *  all of the ownedlocal nodes defined on a processor.
  */
 Epetra_Map *GbNodetoOwnedLcNodeMap;

  //! Epetra_Map object containing the mapping of the block node
  //! equations corresponding  to all column rows equations needed by this processor
  //! whether they are owned or are external
  /*!
   *  This describes the distribution of the block node equations corresponding to
   *  column nodes needed by this processor for each processor. This is not a
   *  1 to 1 map.
   *
   *  This is a distributed object. It's used to set up the columns of the
   *  matrix owned by this processor.
   *
   *  This is called the overlap map for the solution vector. Might switch
   *  to this nomenclature in the future.
   *
   *  This is the main map for setting up ghosted solution vectors within the code.
   *
   *  It called the Epetra_Vector_Ghosted type within m1d.
   */
  Epetra_BlockMap *GbBlockNodeEqnstoLcBlockNodeEqnsColMap;

  //! Epetra_Map object containing the mapping of the block node
  //! equations corresponding  to column rows equations owned by this processor
  /*!
   *  This describes the distribution of the block node equations corresponding to
   *  row nodes needed by this processor for each processor. This is a
   *  1 to 1 map with the solution vector. It is used to residual vectors and some
   *  types of solution vectors such as delta variable solution vectors.
   *
   *  This is a distributed object. It's used to set up the row of the
   *  matrix owned by this processor.
   *
   *  This is called the domain map for the solution vector. Might switch
   *  to this nomenclature in the future.
   *
   *  It called the Epetra_Vector_Owned type within m1d.
   */
  Epetra_BlockMap *GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap;

  //! Specifies the communications pattern for  Imports/updates of
  //! ghost unknowns on processors.
  /*!
   *
   */
  Epetra_Import *Importer_GhostEqns;

  //! Specifies the communications pattern for  Imports/updates of
  //! ghost unknowns on processors.
  /*!
   *
   */
  Epetra_Import *Importer_NodalValues;

  //! Number of Colors 
  int NumNodeColors;

  //! Color map for nodes, both external and internal on the processor
  /*!
   *   This coloring map will be used to calculate the jacobian.
   */
  Epetra_MapColoring *NodeColorMap;

  //! Number of equations colors
  int NumEqnColors;

  //! Color map for equations, both external and internal on the processor
  /*!
   *   This coloring map will be used to calculate the Jacobian. Note, This is currently
   *   broken, because the coloring scheme only works on elements within the Epetra
   *   object, and not on points. DON'T USE -> JUST KEEPING IT AROUND SO WE CAN
   *   WORK IT OUT WITH THE DEVELOPERS.
   */
  Epetra_MapColoring *EqnColorMap;

  //! Color map for equations, both external and internal on the processor
  /*!
   *   This coloring map will be used to calculate the Jacobian. Note, This is the fixup
   *   vector that we have implemented because Epetra_MapColoring is currently broken.
   */
  Epetra_IntVector *EqnColors;

  //! Initial Position of the nodes
  /*!
   *   Initial position of the nodes. There may be displacements that alter the
   *   position
   */
  std::vector<double> X0pos_LcNode;

  //! Actual Position of the nodes
  /*!
   *   Initial position of the nodes. There may be displacements that alter the
   *   position. Length is the number of nodes defined on the processor
   */
  Epetra_Vector *Xpos_LcNode_p;

  Epetra_Vector *Xpos_LcOwnedNode_p;

  //! Vector of nodal variables descriptors.
  /*!
   * NodalVars object contains a complete description of the variables at
   * a local node. These are shallow copies. The actual objects are held by
   * the GlobalIndices object.
   *
   * length = number of local nodes identified on the current processor.
   */
  std::vector<NodalVars *> NodalVars_LcNode;

  //! Pointer to the global index object
  GlobalIndices *GI_ptr_;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
