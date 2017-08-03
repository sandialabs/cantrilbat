/*
 * m1d_SimpleTDDiff_dom1D.cpp
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 */

//! This is a heavyweight base class that provides the function
//! evaluation for a single bulk domain.
#include "m1d_SimpleTDDiff_dom1D.h"

#include "m1d_NodalVars.h"
#include "m1d_LocalNodeIndices.h"

#include "m1d_GlobalIndices.h"
#include "m1d_exception.h"

#include "cantera/base/ctml.h"

#include "stdio.h"
#include "stdlib.h"

using namespace std;
using namespace m1d;

namespace m1d
{

//==============================================================================
SimpleTDDiff_dom1D::SimpleTDDiff_dom1D(BulkDomainDescription* bdd_ptr) :
  BulkDomain1D(bdd_ptr)
{

}
//==============================================================================
SimpleTDDiff_dom1D::SimpleTDDiff_dom1D(const SimpleTDDiff_dom1D &r) :
  BulkDomain1D(r.BDD_ptr_)
{
  SimpleTDDiff_dom1D::operator=(r);
}
//==============================================================================
SimpleTDDiff_dom1D::~SimpleTDDiff_dom1D() {
}
//==============================================================================
SimpleTDDiff_dom1D &
SimpleTDDiff_dom1D::operator=(const SimpleTDDiff_dom1D &r)
{
  if (this == &r) {
    return *this;
  }
  // Call the parent assignment operator
  BulkDomain1D::operator=(r);

  return *this;
}
//==============================================================================
// Prepare all of the indices for fast calculation of the residual
/*
 *  Ok, at this point, we will have figured out the number of equations
 *  to be calculated at each node point. The object NodalVars will have
 *  been fully formed.
 *
 *  We use this to figure out what local node numbers/ cell numbers are
 *  needed and to set up indices for their efficient calling.
 *
 *  Child objects of this one will normally call this routine in a
 *  recursive fashion.
 */
void
SimpleTDDiff_dom1D::domain_prep(LocalNodeIndices *li_ptr)
{
  /*
   * First call the parent domain prep to get the node information
   */
  BulkDomain1D::domain_prep(li_ptr);
}
//==============================================================================
// Basic function to calculate the residual for the domain.
/*
 *
 */
void
SimpleTDDiff_dom1D::residEval(Epetra_Vector &res,
                              const bool doTimeDependentResid,
                              const Epetra_Vector *soln_ptr,
                              const Epetra_Vector *solnDot_ptr,
                              const Epetra_Vector *solnOld_ptr,
                              const double t,
                              const double rdelta_t,
                              const ResidEval_Type_Enum residType,
			      const Zuzax::Solve_Type solveType)
{
  residType_Curr_ = residType;
  const double surfArea = 1.0;
  int index_RightLcNode;
  int index_LeftLcNode;
  int index_CentLcNode;
  // NodalVars *nodeRight = 0;
  NodalVars *nodeLeft = 0;
  NodalVars *nodeCent = 0;
  NodalVars *nodeRight = 0;
  double xdelL; // Distance of the cell from the center node to the left cell boundary
  double xdelR; // Distance of the cell from the center node to the right cell boundary
  double xdelCell; // Distance of the cell - right boundary minus the left boundary.
  double Cright;
  double Cleft;
  double Ccent;
  double fluxR = 0.0;
  double DiffCoeff = 1.0;
  const Epetra_Vector &soln = *soln_ptr;

  int indexLeft_EqnStart;
  int indexCent_EqnStart;
  int indexRight_EqnStart;

  /*
   * Special section to do the left boundary
   */

  double fluxL = 0.0;
  fluxR = fluxL;

  /*
   * Special section if we own the left node of the domain. If we do
   * it will be cell 0
   */
  if (IOwnLeft) {
    DiffFluxLeftBound_LastResid_NE[0] = fluxL;
  }

  bool doLeftFluxCalc = true;
  /*
   *  Loop over the number of Cells in this domain on this processor
   *  This loop is done from left to right.
   */
  for (int iCell = 0; iCell < NumLcCells; iCell++) {

    /*
     * get the index for the left node
     */
    index_LeftLcNode = Index_LeftLcNode_LCO[iCell];
    if (index_LeftLcNode < 0) {
      nodeLeft = 0;
      indexLeft_EqnStart = -1;

    }
    else {
      // get the node structure for the left node
      nodeLeft = LI_ptr_->NodalVars_LcNode[index_LeftLcNode];
      /*
       * Get the starting index of the equation for the left node
       */
      indexLeft_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_LeftLcNode]
          + nodeLeft->OffsetIndex_BulkDomainEqnStart_BDN[0];
    }
    /*
     * If we are past the first cell, then we have already done the calculation
     * for this flux at the right cell edge of the previous cell
     */
    if (iCell > 0) {
      doLeftFluxCalc = false;
    }

    /*
     * Get the index for the center node
     */
    index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
        + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];

    /*
     * Get the indexes for the right node
     */
    index_RightLcNode = Index_RightLcNode_LCO[iCell];
    if (index_RightLcNode < 0) {
      nodeRight = 0;
    }
    else {
      nodeRight = LI_ptr_->NodalVars_LcNode[index_RightLcNode];
      indexRight_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_RightLcNode]
          + nodeRight->OffsetIndex_BulkDomainEqnStart_BDN[0];
    }
    /*
     * --------------------------------------------------------------------------
     * Add in the fluxes from the left boundary
     * --------------------------------------------------------------------------
     */
    /*
     * Calculate the distance between the left and center node points
     */
    if (nodeLeft) {
      xdelL = nodeCent->xNodePos() - nodeLeft->xNodePos();
    }
    else {
      xdelL = 0.0;
    }

    if (doLeftFluxCalc) {
      if (nodeLeft == 0) {
        /*
         *  We are here if we are at the left node boundary and we need a flux
         *  condition. The default now is to set the flux to zero. We could
         *  put in a more sophisticated treatment
         */
        fluxR = 0.0;
        fluxL = 0.0;
      }
      else {

        /*
         *  Establish the environment at the left cell boundary
         */
        // Calculate temperature

        // Calculate pressure

        // Calculate concentrations

        // Set the Cantera state

        // Calculate the density of the fluid

        /*
         * Calculate the transport coefficients at the left cell boundary
         */
        DiffCoeff = 1.0;
        /*
         *  Add in the fluxes from the left cell boundary
         */
        Cleft = soln[indexLeft_EqnStart];

        Ccent = soln[indexCent_EqnStart];
        /*
         * Calculate the flux at the right boundary by transferring from right to left
         */
        fluxL = -DiffCoeff * surfArea * (Ccent - Cleft) / (xdelL);
      }
    }
    else {
      /*
       * Copy the fluxes from the storred right side
       */
      fluxL = fluxR;
    }

    /*
     * Add the source term into the residual
     */
    res[indexCent_EqnStart] -= fluxL;

    /*
     * --------------------------------------------------------------------------
     *  Add in the source terms at the current cell center
     * --------------------------------------------------------------------------
     */
    /*
     *  Establish the environment at the center of the cell
     */
    if (nodeRight == 0) {
      xdelR = 0.0;
    }
    else {
      xdelR = nodeRight->xNodePos() - nodeCent->xNodePos();
    }
    xdelCell = xdelL + xdelR;

    // Calculate temperature

    // Calculate pressure

    // Calculate concentrations

    // Set the Cantera state

    // Calculate the density of the fluid

    /*
     * Supplementary case where we are adding on time dependent terms to the residual
     */
    if (doTimeDependentResid) {
      double CDotcent = (*solnDot_ptr)[indexCent_EqnStart];
      res[indexCent_EqnStart] += xdelCell * surfArea * CDotcent;
    }

    /*
     * --------------------------------------------------------------------------
     * Add in the fluxes from the right boundary
     * --------------------------------------------------------------------------
     */

    /*
     * Calculate the flux distance
     */
    if (nodeRight == 0) {
      /*
       *  We are here if we are at the right node boundary and we need a flux
       *  condition. The default now is to set the flux to zero. We could
       *  put in a more sophisticated treatment
       */
      AssertTrace(iCell == NumLcCells-1);
      fluxR = 0.0;
    }
    else {
      /*
       *  Establish the environment at the right cell boundary
       */
      // Calculate temperature

      // Calculate pressure

      // Calculate concentrations

      // Set the Cantera state

      // Calculate the density of the fluid

      /*
       * Calculate the transport coefficients at the right cell boundary
       */
      DiffCoeff = 1.0;
      /*
       *  Add in the fluxes from the left cell boundary
       */
      Cright = soln[indexRight_EqnStart];

      Ccent = soln[indexCent_EqnStart];
      /*
       * Establish the environment at the right node
       *
       */
      // Calculate temperature

      // Calculate pressure

      // Calculate concentrations

      // Set the Cantera state

      // Calculate the density of the fluid

      /*
       * Calculate the fluxes
       */
      fluxR = -DiffCoeff * surfArea * (Cright - Ccent) / (xdelR);

    }
    /*
     * Add the flux term into the residual
     */
    res[indexCent_EqnStart] += fluxR;

  }

  /*
   * Special section to do the right boundary
   */
  /*
   * Special section if we own the left node of the domain. If we do
   * it will be last cell
   */
  if (IOwnRight) {
    DiffFluxRightBound_LastResid_NE[0] = fluxR;
  }

}
//=====================================================================================================================
// saving the solution on the domain in an xml node.
/*
 *
 * @param oNode                Reference to the XML_Node
 * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
 * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
 * @param t                    time
 *
 * @param duplicateOnAllProcs  If this is true, all processors will include
 *                             the same XML_Node information as proc 0. If
 *                             false, the xml_node info will only exist on proc 0.
 */
void
SimpleTDDiff_dom1D::saveDomain(ZZCantera::XML_Node& oNode,
                             const Epetra_Vector *soln_GLALL_ptr,
                             const Epetra_Vector *solnDot_GLALL_ptr,
                             const double t,
                             bool duplicateOnAllProcs)
{
  // get the NodeVars object pertaining to this global node
  GlobalIndices *gi = LI_ptr_->GI_ptr_;

  // Add a child for this domain
  ZZCantera::XML_Node& bdom = oNode.addChild("domain");

  // Number of equations per node
  int numEquationsPerNode = BDD_ptr_->NumEquationsPerNode;

  // Vector containing the variable names as they appear in the solution vector
  std::vector<VarType> &variableNameList = BDD_ptr_->VariableNameList;

  //! First global node of this bulk domain
  int firstGbNode = BDD_ptr_->FirstGbNode;

  //! Last Global node of this bulk domain
  int lastGbNode = BDD_ptr_->LastGbNode;
  int numNodes = lastGbNode - firstGbNode + 1;

  bdom.addAttribute("id", id());
  bdom.addAttribute("points", numNodes);
  bdom.addAttribute("type", "bulk");
  bdom.addAttribute("numVariables", numEquationsPerNode);

  // Dump out the coordinates
  ZZCantera::XML_Node& gv = bdom.addChild("grid_data");

  std::vector<double> varContig(numNodes);

  int i = 0;
  for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
    NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
    varContig[i] = nv->x0NodePos();
  }
  ZZctml::addNamedFloatArray(gv, "X0", varContig.size(), &(varContig[0]), "m", "length");

  for (int iVar = 0; iVar < numEquationsPerNode; iVar++) {
    VarType vt = variableNameList[iVar];
    i = 0;
    std::string nmm = vt.VariableName(200);
    for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
      int istart = nv->EqnStart_GbEqnIndex;
      varContig[i] = (*soln_GLALL_ptr)[istart+ibulk+iVar];
    }
    ZZctml::addNamedFloatArray(gv, nmm, varContig.size(), &(varContig[0]), "kmol/m3", "concentration");

  }
}
//=================================================================================
void
SimpleTDDiff_dom1D::err(const char *msg)
{
  printf("SimpleTDDiff_dom1D: function not implemented: %s\n", msg);
  exit(-1);
}
//=================================================================================
//=================================================================================
}
//=================================================================================


