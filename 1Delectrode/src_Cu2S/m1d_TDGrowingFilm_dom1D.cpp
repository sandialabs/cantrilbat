/*
 * m1d_TDGrowingFilm_dom1D.cpp
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 */

//! This is a heavyweight base class that provides the function
//! evaluation for a single bulk domain.
#include "m1d_TDGrowingFilm_dom1D.h"

#include "m1d_NodalVars.h"

#include "m1d_DomainLayout.h"

using namespace std;
using namespace m1d;

namespace m1d
{

//==============================================================================
TDGrowingFilm_dom1D::TDGrowingFilm_dom1D(BulkDomainDescription* bdd_ptr) :
  BulkDomain1D(bdd_ptr)
{

}
//==============================================================================
TDGrowingFilm_dom1D::TDGrowingFilm_dom1D(const TDGrowingFilm_dom1D &r) :
  BulkDomain1D(r.BDD_ptr_)
{
  TDGrowingFilm_dom1D::operator=(r);
}
//==============================================================================
TDGrowingFilm_dom1D::~TDGrowingFilm_dom1D()
{
}
//==============================================================================
TDGrowingFilm_dom1D &
TDGrowingFilm_dom1D::operator=(const TDGrowingFilm_dom1D &r)
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
TDGrowingFilm_dom1D::domain_prep(LocalNodeIndices *li_ptr)
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
TDGrowingFilm_dom1D::residEval(Epetra_Vector &res, const bool doTimeDependentResid,
                               const Epetra_Vector* const soln_ptr,
                               const Epetra_Vector* const solnDot_ptr,
                               const Epetra_Vector* const solnOld_ptr,
                               const double t,
                               const double rdelta_t,
                               const Zuzax::ResidEval_Type residType,
			       const Zuzax::Solve_Type solveType)
{
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
  double xCellBoundaryL;
  double xCellBoundaryR;
  double Cright;
  double Cleft;
  double Ccent;
  double fluxR = 0.0;
  double DiffCoeff = 1.0;
  const Epetra_Vector &soln = *soln_ptr;
  bool doMeshConTerm = true;

  int indexLeft_EqnStart;
  int indexCent_EqnStart;
  int indexRight_EqnStart;

  // Will have to generalize this indexing capability
  int indexMeshDisplacement = 0;
  int indexConc0 = 1;
  DomainLayout *dl = BDD_ptr_->DL_ptr_;

  incrementCounters(residType);

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

  /*
   * Calculate the cell velocities
   */
  std::vector<double> cellboundVelocL(NumLcCells, 0.0);
  std::vector<double> cellboundVelocR(NumLcCells, 0.0);
  std::vector<double> cellboundL(NumLcCells, 0.0);
  std::vector<double> cellboundR(NumLcCells, 0.0);
  std::vector<double> cellboundL_old(NumLcCells, 0.0);
  std::vector<double> cellboundR_old(NumLcCells, 0.0);
  double nodeVelocC = 0.0, nodePosC, nodePosC_old;
  double nodeVelocL, nodePosL, nodePosL_old;
  double nodeVelocR, nodePosR, nodePosR_old;
  for (int iCell = 0; iCell < NumLcCells; iCell++) {

    index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
        + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
    if (solnDot_ptr) {
      nodeVelocC = (*solnDot_ptr)[indexCent_EqnStart + indexMeshDisplacement];
    }
    nodePosC = (*soln_ptr)[indexCent_EqnStart + indexMeshDisplacement] + nodeCent->x0NodePos();
    AssertTrace(nodePosC == nodeCent->xNodePos());
    if (solnOld_ptr) {
      nodePosC_old = (*solnOld_ptr)[indexCent_EqnStart + indexMeshDisplacement] + nodeCent->x0NodePos();
    } else {
      nodePosC_old = nodePosC;
    }

    index_LeftLcNode = Index_LeftLcNode_LCO[iCell];
    if (index_LeftLcNode >= 0) {
      nodeLeft = LI_ptr_->NodalVars_LcNode[index_LeftLcNode];
      indexLeft_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_LeftLcNode]
          + nodeLeft->OffsetIndex_BulkDomainEqnStart_BDN[0];
      if (solnDot_ptr) {
        nodeVelocL = (*solnDot_ptr)[indexLeft_EqnStart + indexMeshDisplacement];
      }
      nodePosL = (*soln_ptr)[indexLeft_EqnStart + indexMeshDisplacement] + nodeLeft->x0NodePos();
      //AssertTrace(nodePosL == nodeLeft->xNodePos());
      if (nodePosL != nodeLeft->xNodePos()) {
	  printf("we are here %g %g \n", nodePosL, nodeLeft->xNodePos());
      }
      if (solnOld_ptr) {
        nodePosL_old = (*solnOld_ptr)[indexLeft_EqnStart + indexMeshDisplacement] + nodeLeft->x0NodePos();
      } else {
        nodePosL_old = nodePosL;
      }
    } else {
      nodeVelocL = nodeVelocC;
      nodePosL = nodePosC;
      nodePosL_old = nodePosC_old;
    }

    index_RightLcNode = Index_RightLcNode_LCO[iCell];
    if (index_RightLcNode >= 0) {
      nodeRight = LI_ptr_->NodalVars_LcNode[index_RightLcNode];
      indexRight_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_RightLcNode]
          + nodeRight->OffsetIndex_BulkDomainEqnStart_BDN[0];
      if (solnDot_ptr) {
        nodeVelocR = (*solnDot_ptr)[indexRight_EqnStart + indexMeshDisplacement];
      }
      nodePosR = (*soln_ptr)[indexRight_EqnStart + indexMeshDisplacement] + nodeRight->x0NodePos();
      AssertTrace(nodePosR == nodeRight->xNodePos());
      if (solnOld_ptr) {
        nodePosR_old = (*solnOld_ptr)[indexRight_EqnStart + indexMeshDisplacement]+ nodeRight->x0NodePos();
      } else {
        nodePosR_old = nodePosR;
      }
    } else {
      nodeVelocR = nodeVelocC;
      nodePosR = nodePosC;
      nodePosR_old = nodePosC_old;
    }

    cellboundL[iCell] = 0.5 * (nodePosL + nodePosC);
    cellboundL_old[iCell] = 0.5 * (nodePosL_old + nodePosC_old);

    cellboundR[iCell] = 0.5 * (nodePosR + nodePosC);
    cellboundR_old[iCell] = 0.5 * (nodePosR_old + nodePosC_old);
    if (solnDot_ptr) {
      cellboundVelocL[iCell] = 0.5 * (nodeVelocL + nodeVelocC);
      cellboundVelocR[iCell] = 0.5 * (nodeVelocR + nodeVelocC);
    }
  }

  bool doLeftFluxCalc = true;
  /*
   *  Loop over the number of Cells in this domain on this processor
   *  This loop is done from left to right.
   */
  for (int iCell = 0; iCell < NumLcCells; iCell++) {

#ifdef DEBUG_HKM
    if (counterResBaseCalcs_ > 125 && residType == Zuzax::ResidEval_Type::Base_ResidEval) {
      if (iCell == NumLcCells - 1) {
        // printf("we are here\n");
      }
    }
#endif

    /*
     * Get the index for the center node
     */
    index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
        + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
    /*
     * get the index for the left node
     */
    index_LeftLcNode = Index_LeftLcNode_LCO[iCell];
    if (index_LeftLcNode < 0) {
      nodeLeft = 0;
      indexLeft_EqnStart = -1;
    } else {
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
     * Get the indexes for the right node
     */
    index_RightLcNode = Index_RightLcNode_LCO[iCell];
    if (index_RightLcNode < 0) {
      nodeRight = 0;
      indexRight_EqnStart = -1;
    } else {
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
      xCellBoundaryL = 0.5 * (nodeLeft->xNodePos() + nodeCent->xNodePos());
    } else {
      xdelL = 0.0;
      xCellBoundaryL = 0.0;
    }

    Ccent = soln[indexCent_EqnStart + indexConc0];
    if (indexLeft_EqnStart >= 0) {
      Cleft = soln[indexLeft_EqnStart + indexConc0];
    } else {
      Cleft = Ccent;
    }
    if (indexRight_EqnStart >= 0) {
      Cright = soln[indexRight_EqnStart + indexConc0];
    } else {
      Cright = Ccent;
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
      } else {

        /*
         *  Establish the environment at the left cell boundary
         */
        // Calculate temperature

        // Calculate pressure

        // Calculate concentrations

        // Set the Zuzax state

        // Calculate the density of the fluid

        /*
         * Calculate the transport coefficients at the left cell boundary
         */
        DiffCoeff = 1.0;

        /*
         * Calculate the flux at the right boundary by transferring from right to left
         */
        fluxL = -DiffCoeff * surfArea * (Ccent - Cleft) / (xdelL);
      }
    } else {
      /*
       * Copy the fluxes from the storred right side
       */
      fluxL = fluxR;
    }

    /*
     * Add the left flux term into the residual
     */
    res[indexCent_EqnStart + indexConc0] -= fluxL;

    /*
     * Convective flux - mesh movement with NO material movement
     */
    if (doMeshConTerm) {
      if (cellboundVelocL[iCell] >= 0.0) {
        double flux = (cellboundL[iCell] - cellboundL_old[iCell]) * Ccent * rdelta_t;
        res[indexCent_EqnStart + indexConc0] += flux;
      } else {
        double flux = (cellboundL[iCell] - cellboundL_old[iCell]) * Cleft * rdelta_t;
        res[indexCent_EqnStart + indexConc0] += flux;
      }
    }
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
      xCellBoundaryR = nodeCent->xNodePos();
    } else {
      xdelR = nodeRight->xNodePos() - nodeCent->xNodePos();
      xCellBoundaryR = 0.5 * (nodeRight->xNodePos() + nodeCent->xNodePos());
    }
    xdelCell = xCellBoundaryR - xCellBoundaryL;

    // Calculate temperature

    // Calculate pressure

    // Calculate concentrations

    // Set the Zuzax state

    // Calculate the density of the fluid

    /*
     * Supplementary case where we are adding on time dependent terms to the residual
     */
    if (doTimeDependentResid) {
      double newStuff = Ccent * xdelCell;
      double Ccent_old = (*solnOld_ptr)[indexCent_EqnStart + indexConc0];
      double xdelCell_old = cellboundR_old[iCell] - cellboundL_old[iCell];
      double oldStuff = Ccent_old * xdelCell_old;
      res[indexCent_EqnStart + indexConc0] += surfArea * (newStuff - oldStuff) * rdelta_t;
    }

    // Do the mesh equation for the node in the middle of the cell
    //          res[indexCent_EqnStart+indexMeshDisplacment] =m_nodePos[i] - m_xbot0
    //                        - m_fracNodePos[i] / m_fracNodePos[i+ 1] * (m_nodePos[i + 1] - m_xbot0);
    //
    if (nodeRight != 0) {
      res[indexCent_EqnStart + indexMeshDisplacement] = (nodeCent->xNodePos() - dl->XLoc_LeftBoundary)
          - (nodeCent->xFracNodePos() / nodeRight->xFracNodePos()) * (nodeRight->xNodePos() - dl->XLoc_LeftBoundary);
    } else {
      /*
       *  This is a dummy statement that fills out the specification of the residual. However,
       *  It's needed for a stable zero-state for the system.
       */
      res[indexCent_EqnStart + indexMeshDisplacement] = nodeCent->xNodePos() - nodeCent->x0NodePos();
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
    } else {
      /*
       *  Establish the environment at the right cell boundary
       */
      // Calculate temperature

      // Calculate pressure

      // Calculate concentrations

      // Set the Zuzax state

      // Calculate the density of the fluid

      /*
       * Calculate the transport coefficients at the right cell boundary
       */
      DiffCoeff = 1.0;
      /*
       *  Add in the fluxes from the left cell boundary
       */
      // Cright = soln[indexRight_EqnStart + indexConc0];

      // Ccent = soln[indexCent_EqnStart + indexConc0];
      /*
       * Establish the environment at the right node
       *
       */
      // Calculate temperature

      // Calculate pressure

      // Calculate concentrations

      // Set the Zuzax state

      // Calculate the density of the fluid

      /*
       * Calculate the fluxes
       */
      fluxR = -DiffCoeff * surfArea * (Cright - Ccent) / (xdelR);

    }
    /*
     * Add the flux term into the residual
     */
    res[indexCent_EqnStart + indexConc0] += fluxR;

    /*
     * Convective flux - mesh movement with NO material movement
     */
    if (nodeRight) {
      if (doMeshConTerm) {
        if (cellboundVelocR[iCell] >= 0.0) {
          double flux = (cellboundR[iCell] - cellboundR_old[iCell]) * Cright * rdelta_t;
          res[indexCent_EqnStart + indexConc0] -= flux;
        } else {
          double flux = (cellboundR[iCell] - cellboundR_old[iCell]) * Ccent * rdelta_t;
          res[indexCent_EqnStart + indexConc0] -= flux;
        }
      }
    }
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
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void TDGrowingFilm_dom1D::setAtolVector(double atolDefault, const Epetra_Vector_Ghosted & soln, 
					Epetra_Vector_Ghosted & atolVector,
					const Epetra_Vector_Ghosted * const atolV)
{
  double xmin = 1E300;
  double xmax = -1.0E300;
  /*
   * Offsets for the variable unknowns in the solution vector for the bulk domain
   */
  int iVar_MeshDisplacement_BD = 0;
  int iVar_Conc_BD = 1;
  
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    NodalVars * nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    int indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
      + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
    double nodePosC = soln[indexCent_EqnStart_BD + iVar_MeshDisplacement_BD] + nodeCent->x0NodePos();
    if ( nodePosC <xmin) {
      xmin = nodePosC;
    }
    if (nodePosC > xmax) {
      xmax = nodePosC;
    }

  }
  double deltaX = xmax - xmin;


  for (int iCell = 0; iCell < NumLcCells; iCell++) {

    /*
     *    Get the index for the center node 
     */
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    
    /*
     *   Get the pointer to the NodalVars object for the center node
     */
    NodalVars * nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

    /*
     *  Index of the first equation in the bulk domain of center node
     */
    int indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
      + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
  
  
    atolVector[indexCent_EqnStart_BD + iVar_Conc_BD] = 1.0E-14;
    atolVector[indexCent_EqnStart_BD + iVar_MeshDisplacement_BD] = 1.0 * deltaX;
  }
}
//=====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void TDGrowingFilm_dom1D::setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted & soln, 
						const Epetra_Vector_Ghosted & solnDot, 
						Epetra_Vector_Ghosted & atolVector,
						const Epetra_Vector_Ghosted * const atolV)
{
  double xmin = 1E300;
  double xmax = -1.0E300;
  /*
   * Offsets for the variable unknowns in the solution vector for the bulk domain
   */
  int iVar_MeshDisplacement_BD = 0;
  int iVar_Conc_BD = 1;
  
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    NodalVars * nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    int indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
      + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
    double nodePosC = soln[indexCent_EqnStart_BD + iVar_MeshDisplacement_BD] + nodeCent->x0NodePos();
    if ( nodePosC <xmin) {
      xmin = nodePosC;
    }
    if (nodePosC > xmax) {
      xmax = nodePosC;
    }

  }
  double deltaX = xmax - xmin;


  for (int iCell = 0; iCell < NumLcCells; iCell++) {

    /*
     *    Get the index for the center node 
     */
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    
    /*
     *   Get the pointer to the NodalVars object for the center node
     */
    NodalVars * nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

    /*
     *  Index of the first equation in the bulk domain of center node
     */
    int indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
      + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
  
  
    atolVector[indexCent_EqnStart_BD + iVar_Conc_BD] = 1.0E-14;
    atolVector[indexCent_EqnStart_BD + iVar_MeshDisplacement_BD] = 1.0 * deltaX;
  }
}
//=====================================================================================================================
void
TDGrowingFilm_dom1D::setAtolDeltaDamping(double atolDefault, double relcoeff, 
					 const Epetra_Vector_Ghosted & soln, 
					 Epetra_Vector_Ghosted & atolDeltaDamping,
					 const Epetra_Vector_Ghosted * const atolV)
{
  double xmin = 1E300;
  double xmax = -1.0E300;
 
  /*
   * Offsets for the variable unknowns in the solution vector for the bulk domain
   */
  int iVar_MeshDisplacement_BD = 0;
  int iVar_Conc_BD = 1;
  
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    NodalVars * nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    int indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
      + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
    double nodePosC = soln[indexCent_EqnStart_BD + iVar_MeshDisplacement_BD] + nodeCent->x0NodePos();
    if ( nodePosC <xmin) {
      xmin = nodePosC;
    }
    if (nodePosC > xmax) {
      xmax = nodePosC;
    }

  }
  double deltaX = xmax - xmin;
  for (int iCell = 0; iCell < NumLcCells; iCell++) {


    /*
     *  ---------------- Get the index for the center node ---------------------------------
     */
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    /*
     *   Get the pointer to the NodalVars object for the center node
     */
    NodalVars * nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

    /*
     *  Index of the first equation in the bulk domain of center node
     */
    int indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
        + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
  

   

    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Conc_BD] = 1.0E-8 * relcoeff;
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_MeshDisplacement_BD] = 1.0 * deltaX * relcoeff;
  }
}
//=====================================================================================================================
void
TDGrowingFilm_dom1D::setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff, 
						 const Epetra_Vector_Ghosted & soln, 
						 const Epetra_Vector_Ghosted & solnDot, 
						 Epetra_Vector_Ghosted & atolDeltaDamping,
						 const Epetra_Vector_Ghosted * const atolV)
{
  double xmin = 1E300;
  double xmax = -1.0E300;
 
  /*
   * Offsets for the variable unknowns in the solution vector for the bulk domain
   */
  int iVar_MeshDisplacement_BD = 0;
  int iVar_Conc_BD = 1;
  
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    NodalVars * nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    int indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
      + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
    double nodePosC = soln[indexCent_EqnStart_BD + iVar_MeshDisplacement_BD] + nodeCent->x0NodePos();
    if ( nodePosC <xmin) {
      xmin = nodePosC;
    }
    if (nodePosC > xmax) {
      xmax = nodePosC;
    }

  }
  double deltaX = xmax - xmin;
  for (int iCell = 0; iCell < NumLcCells; iCell++) {


    /*
     *  ---------------- Get the index for the center node ---------------------------------
     */
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    /*
     *   Get the pointer to the NodalVars object for the center node
     */
    NodalVars * nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

    /*
     *  Index of the first equation in the bulk domain of center node
     */
    int indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
        + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];

    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Conc_BD] = 1.0E-8 * relcoeff;
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_MeshDisplacement_BD] = 1.0 * deltaX * relcoeff;
  }
}
//=================================================================================
void
TDGrowingFilm_dom1D::err(const char *msg)
{
  printf("TDGrowingFilm_dom1D: function not implemented: %s\n", msg);
  exit(-1);
}
//=================================================================================
//=================================================================================
}
//=================================================================================


