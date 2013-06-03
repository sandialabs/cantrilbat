/**
 * @file m1d_NodalVars.h
 *
 */

/*
 *  $Id: m1d_NodalVars.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
#ifndef M1D_NODALVARS_H
#define M1D_NODALVARS_H

#include "m1d_Eqn_Names.h"
#include "m1d_EqnVarTypes.h"

#include <vector>
#include <map>

namespace m1d
{

class DomainLayout;

//! Class that describes the environment at one node in the domain
/*!
 *
 *   The ordering of the equation at a single node is as follows.
 *
 *      domain_Left equations
 *
 *      domain_Right equations (if not mapped into domain_left equations)
 *
 *      domain_other equations*     (* not yet implemented)
 *  
 *      surfaceDomain_0 equations
 *      surfaceDomain_1 equations
 *       . . .
 *      surfaceDomain_N-1 equations
 *
 *   The surfaceDomain_0 object determines what equations from the domain_Right
 *   bulk domain gets mapped into the domain_Left equations.
 *   
 *   When a bulk equation gets mapped from one bulk domain to another, it means
 *   that there will be created a single conservation equation for the unknown
 *   striding both domains, and the variable is continuous across the interface
 *   as well.
 */
class NodalVars {
public:

  //! Default constructor
  NodalVars(DomainLayout *dl_ptr);

  //! Constructor
  NodalVars(const int gbnode, DomainLayout *dl_ptr);

  //! Copy Constructor
  /*!
   *   @param r  object to be copied.
   */
  NodalVars(const NodalVars &r);

  //! Destructor
  ~NodalVars();

  //! Assignment Operator
  /*!
   *   @param r  object to be copied.
   */
  NodalVars &
  operator=(const NodalVars &r);

  //! This routine will search the global information to discover what
  //! domains are located at this node.
  /*!
   *
   */
  void
  DiscoverDomainsAtThisNode();

  //! This routine discovers and orders the equation unknowns at a node
  /*!
   * This is the ONE PLACE IN THE CODE that discovers and orders the
   * equations at a node.
   */
  void
  GenerateEqnOrder();

  //! Generate a name for the kth variable at the node
  /*!
   *
   * @param k  input the index of the variable on the local node
   * @return  Returns a string representing the variable.
   */
  std::string
  VariableName(int k);

  //! Change the node position
  /*!
   * @param xNodePos  new value of the node position
   */
  void
  changeNodePosition(double xNodePos);

  //! Set up the initial node positions
  /*!
   * @param xNodePos
   * @param x0NodePos
   * @param xFracNodePos
   */
  void
  setupInitialNodePosition(double x0NodePos, double xFracNodePos);

  //! Returns the node position
  double
  xNodePos() const;

  //! Returns the initial node position
  double
  x0NodePos() const;

  //! Return the fraction node position from the left boundary
  double
  xFracNodePos() const;

  //! What Global Node am I
  int GbNode;

  //! Number of equations located at this node
  int NumEquations;

  //! Number of bulk domains located at the node
  int NumBulkDomains;

  //! Number of surface domains located at the node
  int NumSurfDomains;

  //! Starting global index for equations located at this node
  /*!
   * This is in terms of the global equation index, GbEqnIndex,
   * and is independent of processor number
   */
  int EqnStart_GbEqnIndex;

  //! Listing of the bulk domains at this node
  /*!
   *  Length = number of domains, NumBulkDomains
   */
  std::vector<int> BulkDomainIndex_BDN;

  //! Offset index for the equations corresponding
  //! to each domain located at this node
  /*!
   *  Length = number of domains
   */
  std::vector<int> OffsetIndex_BulkDomainEqnStart_BDN;


  //! Map from the bulk domain id to the order index at the current node
  /*!
   *
   */
  std::map<int, int> BulkDomainIndex_fromID;


  //! Listing of the surface domains at this node
  /*!
   *  Length = number of domains, NumDomains
   */
  std::vector<int> SurfDomainIndex_SDN;

  //! Offset index for the equations corresponding
  //! to each domain located at this node
  /*!
   *  Length = number of domains
   */
  std::vector<int> OffsetIndex_SurfDomainEqnStart_SDN;


  //! Map from the surface domain id to the order index at the current node
  /*!
   *
   */
  std::map<int, int> SurfDomainIndex_fromID;


  //! Listing of the Variable Type for each dof at a node
  /*!
   *   This is a listing of the variables types, with the index
   *   being the equation number
   *
   *  Length = number of variables at the node
   */
  std::vector<VarType> VariableNameList_EqnNum;

  //! Listing of the variable subtype
  std::vector<VAR_TYPE_SUBNUM> VariableSubType_EqnNum;

  //! Listing of the Equation Types for each dof at a node
  /*!
   *   This is a listing of the variables types, with the index
   *   being the equation number
   *
   *  Length = number of variables at the node
   */
  std::vector<EqnType> EquationNameList_EqnNum;

  //! Listing of the variable subtype
  std::vector<VAR_TYPE_SUBNUM> EquationSubType_EqnNum;

  //! Map between the variable type and offset of the variable in the unknowns for the node
  /*!
   *   Note dangerous, because the unknowns are not contiguous wrt variable.
   *   Replace with the preferred treatment where the unknown offsets are located in m1d_DomainDescription
   *   They are contiguous wrt domains.
   *       Offset_VarType[MoleFraction_Species] is the offset for the mole fraction variables.
   *
   *      @deprecated
   */
  std::map<VAR_TYPE, int> Offset_VarType;

protected:
  //! Current Spatial position of the node
  /*!
   *
   */
  double XNodePos;

  //! Initial Spatial position of the node
  /*!
   *
   */
  double X0NodePos;

  //! Fraction of the nodal position from the left to the right
  double XFracNodePos;

public:
  //! Right cell boundary position
  //double XcellBoundRight;

  //! Left cell boundary position;
  //double XcellBoundLeft;

  DomainLayout *DL_ptr_;

};

}

#endif
