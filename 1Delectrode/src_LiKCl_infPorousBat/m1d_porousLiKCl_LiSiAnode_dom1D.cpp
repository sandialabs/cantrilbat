/*
 * m1d_porousLiKCl_LiSiAnode_dom1D.cpp
 *
 
 */

//  This is a heavyweight base class that provides the function
//evaluation for a single bulk domain.
#include "m1d_porousLiKCl_LiSiAnode_dom1D.h"
#include "m1d_BDT_porAnode_LiKCl.h"

#include "m1d_NodalVars.h"
#include "m1d_LocalNodeIndices.h"

#include "m1d_GlobalIndices.h"
#include "m1d_exception.h"
#include "m1d_DomainLayout.h"
#include "m1d_Comm.h"

#include "Epetra_Comm.h"

#include "cantera/base/ctml.h"
#include "cantera/transport/Tortuosity.h"

//next two lines added for salt precipitation
#include "cantera/thermo.h"
#include "cantera/thermo/MargulesVPSSTP.h"
extern int flagPrecipitation;

#include "stdio.h"
#include "stdlib.h"

using namespace std;
using namespace Cantera;

#include "m1d_ProblemStatementCell.h"
extern m1d::ProblemStatementCell PSinput;
//=====================================================================================================================
namespace m1d
{

//=====================================================================================================================
porousLiKCl_LiSiAnode_dom1D::porousLiKCl_LiSiAnode_dom1D(BulkDomainDescription & bdd) :
      porousElectrode_dom1D(bdd), 
      ionicLiquid_(0), trans_(0), Electrode_(0), nph_(0), nsp_(0), concTot_cent_(0.0),
  concTot_cent_old_(0.0), porosity_Cell_(0), porosity_Cell_old_(0), surfaceAreaDensity_Cell_(0), icurrInterfacePerSurfaceArea_Cell_(0),
  xdelCell_Cell_(0),
      concTot_Cell_(0), concTot_Cell_old_(0), cIndex_cc_(0), Fleft_cc_(0.0), Fright_cc_(0.0), Vleft_cc_(0.0),
      Vcent_cc_(0.0), Vright_cc_(0.0), VElectrodeLeft_cc_(0.0), VElectrodeCent_cc_(0.0), VElectrodeRight_cc_(0.0),
      Xleft_cc_(0), Xcent_cc_(0), Xright_cc_(0), spCharge_(0), mfElectrolyte_Soln_Curr_(0),
      mfElectrolyte_Thermo_Curr_(0), temp_Curr_(TemperatureReference_), pres_Curr_(0.0), phiElectrolyte_Curr_(-10000.),
      phiElectrode_Curr_(0.0), concTot_Curr_(0.0), porosity_Curr_(0.0), conductivityElectrode_(1.0E6), gradT_trCurr_(0.0),
      gradV_trCurr_(0.0), gradVElectrode_trCurr_(0.0), gradX_trCurr_(0), Vdiff_trCurr_(0), jFlux_trCurr_(0),
      icurrElectrode_trCurr_(0.0), electrodeSpeciesProdRates_(0), icurrInterface_Curr_(0.0), phaseMoleFlux_(0),
      solnMoleFluxInterface_Curr_(0.0), icurrElectrode_CBL_(0), icurrElectrode_CBR_(0), icurrElectrolyte_CBL_(0),
      icurrElectrolyte_CBR_(0), deltaV_Cell_(0), Ess_Cell_(0), overpotential_Cell_(0), icurrRxn_Cell_(0),  LiFlux_Cell_(0), solnTemp(0),
      ivb_(VB_MOLEAVG)
{
  BDT_porAnode_LiKCl *fa = dynamic_cast<BDT_porAnode_LiKCl *> (&bdd);
  if (!fa) {
    throw m1d_Error("confused", "confused");
  }
  /*
   * This is a shallow pointer copy. The BDT object owns the ionicLiquid_ object
   */
  ionicLiquid_ = fa->ionicLiquid_;
  /* 
   *  This is a shallow pointer copy. The BDT object owns the transport object
   */
  trans_ = fa->trans_;
  /*
   *  This is a shallow pointer copy. The BDT object owns the Electrode object
   */
  Electrode_ = fa->Electrode_;
  Electrode_->setID(0, 0);
  nsp_ = 3;
  nph_ = 1;
}
//=====================================================================================================================
porousLiKCl_LiSiAnode_dom1D::porousLiKCl_LiSiAnode_dom1D(const porousLiKCl_LiSiAnode_dom1D &r) :
  porousElectrode_dom1D(r.BDD_), ionicLiquid_(0), trans_(0), Electrode_(0), nph_(0), nsp_(0), concTot_cent_(0.0),
  concTot_cent_old_(0.0), porosity_Cell_(0), porosity_Cell_old_(0), surfaceAreaDensity_Cell_(0), icurrInterfacePerSurfaceArea_Cell_(0), 
  xdelCell_Cell_(0),
      concTot_Cell_(0), concTot_Cell_old_(0), cIndex_cc_(0), Fleft_cc_(0.0), Fright_cc_(0.0), Vleft_cc_(0.0),
      Vcent_cc_(0.0), Vright_cc_(0.0), VElectrodeLeft_cc_(0.0), VElectrodeCent_cc_(0.0), VElectrodeRight_cc_(0.0),
      Xleft_cc_(0), Xcent_cc_(0), Xright_cc_(0), spCharge_(0), mfElectrolyte_Soln_Curr_(0),
      mfElectrolyte_Thermo_Curr_(0), temp_Curr_(TemperatureReference_), pres_Curr_(0.0), phiElectrolyte_Curr_(-10000.),
      phiElectrode_Curr_(0.0), concTot_Curr_(0.0), porosity_Curr_(0.0), conductivityElectrode_(1.0E6), gradT_trCurr_(0.0),
      gradV_trCurr_(0.0), gradVElectrode_trCurr_(0.0), gradX_trCurr_(0), Vdiff_trCurr_(0), jFlux_trCurr_(0),
      icurrElectrode_trCurr_(0.0), electrodeSpeciesProdRates_(0), icurrInterface_Curr_(0.0), phaseMoleFlux_(0),
      solnMoleFluxInterface_Curr_(0.0), icurrElectrode_CBL_(0), icurrElectrode_CBR_(0), icurrElectrolyte_CBL_(0),
      icurrElectrolyte_CBR_(0), deltaV_Cell_(0), Ess_Cell_(0), overpotential_Cell_(0), icurrRxn_Cell_(0),  LiFlux_Cell_(0), solnTemp(0),
      ivb_(VB_MOLEAVG)
{
  porousLiKCl_LiSiAnode_dom1D::operator=(r);
}
//=====================================================================================================================
porousLiKCl_LiSiAnode_dom1D::~porousLiKCl_LiSiAnode_dom1D()
{
}
//=====================================================================================================================
porousLiKCl_LiSiAnode_dom1D &
porousLiKCl_LiSiAnode_dom1D::operator=(const porousLiKCl_LiSiAnode_dom1D &r)
{
  if (this == &r) {
    return *this;
  }
  // Call the parent assignment operator
  BulkDomain1D::operator=(r);

  ionicLiquid_ = r.ionicLiquid_;
  trans_ = r.trans_;
  Electrode_ = r.Electrode_;

  nph_ = r.nph_;
  nsp_ = r.nsp_;
  concTot_cent_ = r.concTot_cent_;
  concTot_cent_old_ = r.concTot_cent_old_;
  porosity_Cell_ = r.porosity_Cell_;
  porosity_Cell_old_ = r.porosity_Cell_old_;
  surfaceAreaDensity_Cell_ = r.surfaceAreaDensity_Cell_;
  icurrInterfacePerSurfaceArea_Cell_ = r.icurrInterfacePerSurfaceArea_Cell_;
  xdelCell_Cell_ = r.xdelCell_Cell_;
  concTot_Cell_ = r.concTot_Cell_;
  concTot_Cell_old_ = r.concTot_Cell_old_;
  capacityDischarged_Cell_ = r.capacityDischarged_Cell_;
  depthOfDischarge_Cell_ = r.depthOfDischarge_Cell_;
  capacityLeft_Cell_ = r.capacityLeft_Cell_;
  capacityZeroDoD_Cell_ = r.capacityZeroDoD_Cell_;
  cIndex_cc_ = r.cIndex_cc_;
  Fleft_cc_ = r.Fleft_cc_;
  Fright_cc_ = r.Fright_cc_;
  Vleft_cc_ = r.Vleft_cc_;
  Vcent_cc_ = r.Vcent_cc_;
  Vright_cc_ = r.Vright_cc_;
  VElectrodeLeft_cc_ = r.VElectrodeLeft_cc_;
  VElectrodeCent_cc_ = r.VElectrodeCent_cc_;
  VElectrodeRight_cc_ = r.VElectrodeRight_cc_;
  Xleft_cc_ = r.Xleft_cc_;
  Xcent_cc_ = r.Xcent_cc_;
  Xright_cc_ = r.Xright_cc_;
  spCharge_ = r.spCharge_;
  mfElectrolyte_Soln_Curr_ = r.mfElectrolyte_Soln_Curr_;
  mfElectrolyte_Thermo_Curr_ = r.mfElectrolyte_Thermo_Curr_;
  temp_Curr_ = r.temp_Curr_;
  pres_Curr_ = r.pres_Curr_;
  phiElectrolyte_Curr_ = r.phiElectrolyte_Curr_;
  phiElectrode_Curr_ = r.phiElectrode_Curr_;
  concTot_Curr_ = r.concTot_Curr_;
  porosity_Curr_ = r.porosity_Curr_;
  conductivityElectrode_ = r.conductivityElectrode_;
  gradT_trCurr_ = r.gradT_trCurr_;
  gradV_trCurr_ = r.gradV_trCurr_;
  gradVElectrode_trCurr_ = r.gradVElectrode_trCurr_;
  gradX_trCurr_ = r.gradX_trCurr_;
  Vdiff_trCurr_ = r.Vdiff_trCurr_;
  jFlux_trCurr_ = r.jFlux_trCurr_;
  icurrElectrode_trCurr_ = r.icurrElectrode_trCurr_;
  electrodeSpeciesProdRates_ = r.electrodeSpeciesProdRates_;
  icurrInterface_Curr_ = r.icurrInterface_Curr_;
  phaseMoleFlux_ = r.phaseMoleFlux_;
  solnMoleFluxInterface_Curr_ = r.solnMoleFluxInterface_Curr_;
  icurrElectrode_CBL_ = r.icurrElectrode_CBL_;
  icurrElectrode_CBR_ = r.icurrElectrode_CBR_;
  icurrElectrolyte_CBL_ = icurrElectrolyte_CBL_;
  icurrElectrolyte_CBR_ = r.icurrElectrolyte_CBR_;
  deltaV_Cell_ = r.deltaV_Cell_;
  Ess_Cell_ = r.Ess_Cell_;
  overpotential_Cell_ = r.overpotential_Cell_;
  icurrRxn_Cell_ = r.icurrRxn_Cell_;
  LiFlux_Cell_ = r.LiFlux_Cell_;
  solnTemp = r.solnTemp;
  ivb_ = r.ivb_;

  return *this;
}
//=====================================================================================================================
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
porousLiKCl_LiSiAnode_dom1D::domain_prep(LocalNodeIndices *li_ptr)
{
  /*
   * First call the parent domain prep to get the node information
   */
  BulkDomain1D::domain_prep(li_ptr);
  
  /*
   * Porous electrode domain prep
   */
  porosity_Cell_.resize(NumLcCells, 0.5);
  porosity_Cell_old_.resize(NumLcCells, 0.5);
  surfaceAreaDensity_Cell_.resize(NumLcCells, 1.0E5);
  icurrInterfacePerSurfaceArea_Cell_.resize(NumLcCells, 0.0);
  xdelCell_Cell_.resize(NumLcCells, 0.0);
  concTot_Cell_.resize(NumLcCells, 0.0);
  concTot_Cell_old_.resize(NumLcCells, 0.0);
  capacityDischarged_Cell_.resize(NumLcCells, 0.0);
  depthOfDischarge_Cell_.resize(NumLcCells, 0.0);
  capacityLeft_Cell_.resize(NumLcCells, 0.0);
  capacityZeroDoD_Cell_.resize(NumLcCells, 0.0);

  Xleft_cc_.resize(3, 0.0);
  Xcent_cc_.resize(3, 0.0);
  Xright_cc_.resize(3, 0.0);
  spCharge_.resize(3, 0.0);
  spCharge_[0] = 1.0;
  spCharge_[1] = 1.0;
  spCharge_[2] = -1.0;
  
  mfElectrolyte_Soln_Curr_.resize(3, 0.0);
  mfElectrolyte_Thermo_Curr_.resize(3, 0.0);
  
  gradX_trCurr_.resize(3, 0.0);
  Vdiff_trCurr_.resize(3, 0.0);
  jFlux_trCurr_.resize(3, 0.0);
  electrodeSpeciesProdRates_.resize(30, 0.0);
  phaseMoleFlux_.resize(10, 0.0);

  solnTemp.resize(10, 0.0);

  for (int i = 0; i < NumLcCells; i++) {
    porosity_Cell_[i] = 0.2;  //hardwired to initial anode porosity for test design
  }

  icurrElectrode_CBL_.resize(NumLcCells, 0.0);
  icurrElectrode_CBR_.resize(NumLcCells, 0.0);
  
  icurrElectrolyte_CBL_.resize(NumLcCells, 0.0);
  icurrElectrolyte_CBR_.resize(NumLcCells, 0.0);
  
  deltaV_Cell_ .resize(NumLcCells, 0.0);
  Ess_Cell_.resize(NumLcCells, 0.0);
  overpotential_Cell_.resize(NumLcCells, 0.0);
  icurrRxn_Cell_.resize(NumLcCells, 0.0);
  LiFlux_Cell_.resize(NumLcCells, 0.0);

  /*
   *  Set the velocity basis of the transport object. We are using
   *  mole-averaged velocities as the basis.
   */
  trans_->setVelocityBasis(ivb_);
}
//=====================================================================================================================
// Basic function to calculate the residual for the domain.
/*
 *  This class is used just for volumetric domains with an electrolyte.
 *
 *  All residual terms are written with the following sign convention
 *  based on keeping the time derivative term positive.
 *
 *       res = dcdt - dc2 /dx2 - src = 0
 *
 * @param res  Output vector containing the residual
 * @param doTimeDependentResid  boolean indicating whether the time
 *                         dependent residual is requested
 * @param soln_ptr     solution vector at which the residual should be
 *                     evaluated
 * @param solnDot_ptr  solution dot vector at which the residual should
 *                     be evaluated.
 *  @param t           time
 *  @param rdelta_t    inverse of delta_t
 *
 */
void
porousLiKCl_LiSiAnode_dom1D::residEval(Epetra_Vector &res,
                                       const bool doTimeDependentResid,
                                       const Epetra_Vector *soln_ptr,
                                       const Epetra_Vector *solnDot_ptr,
                                       const Epetra_Vector *solnOld_ptr,
                                       const double t,
                                       const double rdelta_t,
                                       const ResidEval_Type_Enum residType,
				       const Solve_Type_Enum solveType)
{
  int index_RightLcNode;
  int index_LeftLcNode;
  int index_CentLcNode;

  NodalVars *nodeLeft = 0;
  NodalVars *nodeCent = 0;
  NodalVars *nodeRight = 0;

  double xdelL; // Distance from the center node to the left node
  double xdelR; // Distance from the center node to the right node
  double xdelCell; // cell width - right boundary minus the left boundary.
  double xCellBoundaryL; //cell boundary left
  double xCellBoundaryR; //cell boundary right

  //  Electrolyte mole fluxes - this is c V dot n at the boundaries of the cells
  double fluxFright = 0.;
  double fluxFleft;

  // Flux of current in the electrolyte
  //double fluxVright;
  //double fluxVleft;

  // Flux of current in the electrode phase
  double fluxVElectrodeRight = 0.;
  double fluxVElectrodeLeft;

  //mole fraction fluxes
  std::vector<double> fluxXright(nsp_, 0.0);
  std::vector<double> fluxXleft(nsp_, 0.0);

  double fluxL = 0.0;
  ;
  // double fluxR;

  const Epetra_Vector &soln = *soln_ptr;
  const Epetra_Vector &solnOld = *solnOld_ptr;

  /*
   * Index of the first equation at the left node corresponding to the first bulk domain, which is the electrolyte
   */
  int indexLeft_EqnStart_BD;
  /*
   * Index of the first equation at the center node corresponding to the first bulk domain, which is the electrolyte
   */
  int indexCent_EqnStart_BD;
  /*
   * Index of the first equation at the right node corresponding to the first bulk domain, which is the electrolyte
   */
  int indexRight_EqnStart_BD;

  /*
   *   Find the species index for the first species in the electrode object pertaining to the electrolyte
   */
  int sf = Electrode_->solnPhaseIndex();
  int speciesIndex0 = Electrode_->getGlobalSpeciesIndex(sf, 0);

  /*
   * offset of the electolyte solution unknowns at the current node
   */
  index_CentLcNode = Index_DiagLcNode_LCO[0];
  nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

  /*
   *  Offsets for the equation unknowns in the residual vector for the electrolyte domain
   */
  int EQ_Current_offset_BD = BDD_.EquationIndexStart_EqName[Current_Conservation];
  int EQ_Current_offset_ED = EQ_Current_offset_BD + 1;
  int EQ_TCont_offset_BD = BDD_.EquationIndexStart_EqName[Continuity];
  int EQ_Species_offset_BD = BDD_.EquationIndexStart_EqName[Species_Conservation];
  int EQ_MFSum_offset_BD = BDD_.EquationIndexStart_EqName[MoleFraction_Summation];
  int EQ_ChargeBal_offset_BD = BDD_.EquationIndexStart_EqName[ChargeNeutrality_Summation];

  /*
   * Offsets for the variable unknowns in the solution vector for the electrolyte domain
   */
  int iVAR_Vaxial_BD = BDD_.VariableIndexStart_VarName[Velocity_Axial];
  int iVar_Species_BD = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
  int iVar_Voltage_BD = BDD_.VariableIndexStart_VarName[Voltage];

  incrementCounters(residType);
  Fright_cc_ = 0.0;

  /*
   *  -------------------- Special section to do the left boundary -------------------------------
   */

  /*
   * Special section if we own the left node of the domain. If we do
   * it will be cell 0
   */
  if (IOwnLeft) {
    DiffFluxLeftBound_LastResid_NE[0] = fluxL;
  }
  /*
   *  Boolean flag that specifies whether the left flux should be recalculated.
   *  Usually you don't need to calculate the left flux because it can be copied from the right flux
   *  from the previous cell.
   */
  bool doLeftFluxCalc = true;
#ifdef DEBUG_HKM
  if (residType == Base_ResidEval) {
    if (doTimeDependentResid) {
      // printf("we are here\n");
    }
  }
#endif
  /*
   *  ------------------------------ LOOP OVER CELL -------------------------------------------------
   *  Loop over the number of Cells in this domain on this processor
   *  This loop is done from left to right.
   */
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    cIndex_cc_ = iCell;

#ifdef DEBUG_HKM
    if (counterResBaseCalcs_ > 125 && residType == Base_ResidEval) {
      if (iCell == NumLcCells - 1) {
        // printf("we are here\n");
      }
    }
#endif

    /*
     *  ---------------- Get the index for the center node ---------------------------------
     */
    index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    /*
     *   Get the pointer to the NodalVars object for the center node
     */
    nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

    /*
     *  Index of the first equation in the bulk domain of center node
     */
    indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
        + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
    /*
     *  ------------------- Get the index for the left node -----------------------------
     *    There may not be a left node if we are on the left boundary. In that case
     *    set the pointer to zero and the index to -1. Hopefully, we will get a segfault on an error.
     */
    index_LeftLcNode = Index_LeftLcNode_LCO[iCell];
    if (index_LeftLcNode < 0) {
      /*
       *  We assign node object to zero.
       */
      nodeLeft = 0;
      /*
       *  If there is no left node, we assign the left solution index to the center solution index
       */
      indexLeft_EqnStart_BD = indexCent_EqnStart_BD;
    } else {
      // get the node structure for the left node
      nodeLeft = LI_ptr_->NodalVars_LcNode[index_LeftLcNode];
      //index of first equation in the electrolyte of the left node
      indexLeft_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_LeftLcNode]
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
     * ------------------------ Get the indexes for the right node ------------------------------------
     */
    index_RightLcNode = Index_RightLcNode_LCO[iCell];
    if (index_RightLcNode < 0) {
      nodeRight = 0;
      /*
       *  If there is no right node, we assign the right solution index to the center solution index
       */
      indexRight_EqnStart_BD = indexCent_EqnStart_BD;
    } else {
      //NodalVars
      nodeRight = LI_ptr_->NodalVars_LcNode[index_RightLcNode];
      //index of first equation of right node
      indexRight_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_RightLcNode]
          + nodeRight->OffsetIndex_BulkDomainEqnStart_BDN[0];
    }

    /*
     * --------------------------- CALCULATE POSITION AND DELTA_X Variables -----------------------------
     * Calculate the distance between the left and center node points
     */
    if (nodeLeft) {
      xdelL = nodeCent->xNodePos() - nodeLeft->xNodePos();
      xCellBoundaryL = 0.5 * (nodeLeft->xNodePos() + nodeCent->xNodePos());
    } else {
      xdelL = 0.0;
      xCellBoundaryL = nodeCent->xNodePos();
    }
    /*
     * Calculate the distance between the right and center node points
     */
    if (nodeRight == 0) {
      xdelR = 0.0;
      xCellBoundaryR = nodeCent->xNodePos();
    } else {
      xdelR = nodeRight->xNodePos() - nodeCent->xNodePos();
      xCellBoundaryR = 0.5 * (nodeRight->xNodePos() + nodeCent->xNodePos());
    }
    /*
     * Calculate the cell width
     */
    xdelCell = xCellBoundaryR - xCellBoundaryL;
    xdelCell_Cell_[iCell] =  xdelCell;

    /*
     * --------------------------- DO PRE-SETUPSHOP RASTER OVER LEFT,CENTER,RIGHT -----------------------------
     * Calculate the distance between the left and center node points
     */
    /*
     * Get current velocity, mole fraction, temperature, potential
     * from the solution
     */

    for (int k = 0; k < nsp_; k++) {
      Xcent_cc_[k] = soln[indexCent_EqnStart_BD + iVar_Species_BD + k];
    }
    Vcent_cc_ = soln[indexCent_EqnStart_BD + iVar_Voltage_BD];
    VElectrodeCent_cc_ = soln[indexCent_EqnStart_BD + iVar_Voltage_BD + 1];

    if (nodeLeft != 0) {
      /*
       * Find the velocity located at the left cell boundary. 
       * The left cell boundary velocity is stored at the previous (left)
       * cell index as per our conventions.
       */
      Fleft_cc_ = soln[indexLeft_EqnStart_BD + iVAR_Vaxial_BD];
      for (int k = 0; k < nsp_; k++) {
        Xleft_cc_[k] = soln[indexLeft_EqnStart_BD + iVar_Species_BD + k];
      }
      Vleft_cc_ = soln[indexLeft_EqnStart_BD + iVar_Voltage_BD];
      VElectrodeLeft_cc_ = soln[indexLeft_EqnStart_BD + iVar_Voltage_BD + 1];
    } else {
      /*
       * We are here when we are at the left most part of the boundary. Then, there is no
       * left node. Or, the left node is actually the center node. The boundary conditions
       * on the fluxes are set by boundary condition routines.
       */
      /*
       *   We set the left flux to zero in this case. It may or may not be overridden by
       *   a boundary condition in the adjoining surface domain
       */
      Fleft_cc_ = 0.0;

      for (int k = 0; k < nsp_; k++) {
        Xleft_cc_[k] = Xcent_cc_[k];
      }
      Vleft_cc_ = Vcent_cc_;
      VElectrodeLeft_cc_ = VElectrodeCent_cc_;
    }
    /*
     * The right velocity is stored at the current cell index.
     */
    Fright_cc_ = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];

    if (nodeRight != 0) {

      for (int k = 0; k < nsp_; k++) {
        Xright_cc_[k] = soln[indexRight_EqnStart_BD + iVar_Species_BD + k];
      }
      Vright_cc_ = soln[indexRight_EqnStart_BD + iVar_Voltage_BD];
      VElectrodeRight_cc_ = soln[indexRight_EqnStart_BD + iVar_Voltage_BD + 1];
    } else {

      for (int k = 0; k < nsp_; k++) {
        Xright_cc_[k] = Xcent_cc_[k];
      }
      Vright_cc_ = Vcent_cc_;
      VElectrodeRight_cc_ = VElectrodeCent_cc_;
    }

    /*
     * ------------------- CALCULATE FLUXES AT THE LEFT BOUNDARY -------------------------------
     *
     */
    if (doLeftFluxCalc) {
      if (nodeLeft == 0) {
        /*
         *  We are here if we are at the left node boundary and we
         *  need a flux condition. The default now is to
         *  set the flux to zero. We could put in a more
         *  sophisticated treatment
         */
        fluxFleft = 0.0;
        //fluxVleft = 0.0;
        icurrElectrolyte_CBL_[iCell] = 0.0;
        fluxVElectrodeLeft = 0.0;
        icurrElectrode_CBL_[iCell] = 0.0;
        for (int k = 0; k < nsp_; k++) {
          fluxXleft[k] = 0.0;
        }
      } else {

        /*
         *  Establish the environment at the left cell boundary
         */
        SetupThermoShop2(&(soln[indexLeft_EqnStart_BD]), &(soln[indexCent_EqnStart_BD]), 0);

        SetupTranShop(xdelL, 0);

        /*
         * Calculate the flux at the left boundary for each equation
         */
        fluxFleft = Fleft_cc_ * concTot_Curr_ * porosity_Curr_;

        /*
         * Calculate the flux of species and the flux of charge
         *   - the flux of charge must agree with the flux of species
         */
        //fluxVleft = 0.0;
        icurrElectrolyte_CBL_[iCell] = 0.0;
        fluxVElectrodeLeft = icurrElectrode_trCurr_;
        icurrElectrode_CBL_[iCell] = icurrElectrode_trCurr_;
        for (int k = 0; k < nsp_; k++) {
          fluxXleft[k] = jFlux_trCurr_[k] * porosity_Curr_;
          //fluxVleft += fluxXleft[k] * spCharge_[k];
          icurrElectrolyte_CBL_[iCell] += fluxXleft[k] * spCharge_[k];
          if (Fleft_cc_ > 0.0) {
            fluxXleft[k] += Fleft_cc_ * Xleft_cc_[k] * concTot_Curr_ * porosity_Curr_;
          } else {
            fluxXleft[k] += Fleft_cc_ * Xcent_cc_[k] * concTot_Curr_ * porosity_Curr_;
          }
        }
        icurrElectrolyte_CBL_[iCell] *= (Cantera::Faraday);
      }
    } else {  // !doLeftFluxCalc
      /*
       * Copy the fluxes from the stored right side
       */
      fluxFleft = fluxFright;
      icurrElectrolyte_CBL_[iCell] = icurrElectrolyte_CBR_[iCell - 1];
      fluxVElectrodeLeft = fluxVElectrodeRight;
      icurrElectrode_CBL_[iCell] = icurrElectrode_CBR_[iCell - 1];
      for (int k = 0; k < nsp_; k++) {
        fluxXleft[k] = fluxXright[k];
      }
    }
    /*
     * ------------------- CALCULATE FLUXES AT THE RIGHT BOUNDARY -------------------------------
     *
     */
    if (nodeRight == 0) {
      /*
       *  We are here if we are at the right node boundary and we need a flux
       *  condition. The default now is to set the flux to zero. We could
       *  put in a more sophisticated treatment
       */
      AssertTrace(iCell == NumLcCells-1);
      // fluxFright = 0.0;
      Fright_cc_ = 0.0;
      SetupThermoShop1(&(soln[indexCent_EqnStart_BD]), 0);
      fluxFright = Fright_cc_ * concTot_Curr_ * porosity_Curr_;
      //fluxVright = 0.0;
      icurrElectrolyte_CBR_[iCell] = 0.0;
      fluxVElectrodeRight = 0.0;
      icurrElectrode_CBR_[iCell] = 0.0;
      for (int k = 0; k < nsp_; k++) {
        fluxXright[k] = 0.0;
      }
    } else {
      /*
       *  Establish the environment at the right cell boundary
       */
      SetupThermoShop2(&(soln[indexCent_EqnStart_BD]), &(soln[indexRight_EqnStart_BD]), 1);

      SetupTranShop(xdelR, 1);

      /*
       * Calculate the flux at the right boundary for each equation
       * This is equal to
       *       Conc * Vaxial * phi
       */
      fluxFright = Fright_cc_ * concTot_Curr_ * porosity_Curr_;

      /*
       * Calculate the flux of species and the flux of charge
       *   - the flux of charge must agree with the flux of species
       */
      //fluxVright = 0.0;
      icurrElectrolyte_CBR_[iCell] = 0.0;
      fluxVElectrodeRight = icurrElectrode_trCurr_;
      icurrElectrode_CBR_[iCell] = icurrElectrode_trCurr_;
      for (int k = 0; k < nsp_; k++) {
        fluxXright[k] = jFlux_trCurr_[k];
        // fluxVright += fluxXright[k] * spCharge_[k];
        icurrElectrolyte_CBR_[iCell] += fluxXright[k] * spCharge_[k];
        if (Fright_cc_ > 0.0) {
          //fluxXright[k] += Fright_cc_ * Xcent_cc_[k] * concTot_Curr_ * porosity_Curr_;
          fluxXright[k] += Fright_cc_ * mfElectrolyte_Thermo_Curr_[k] * concTot_Curr_ * porosity_Curr_;
        } else {
          fluxXright[k] += Fright_cc_ * mfElectrolyte_Thermo_Curr_[k] * concTot_Curr_ * porosity_Curr_;
        }
      }
      icurrElectrolyte_CBR_[iCell] *= (Cantera::Faraday);
    }

#ifdef DEBUG_CELLRES_HKM
    if (doTimeDependentResid) {

      if (residType == Base_ResidEval) {
        printf(" Cell = %d, Totalflux_K+ = %10.3e,  Totalflux_Cl- = %10.3e \n", iCell, fluxXright[1], fluxXright[2]);
        printf("           Vmolal = %10.3e, jd_Li+ = %10.3e  jd_K+ = %10.3e jd_Cl- = %10.3e\n", Fright_cc_,
               jFlux_trCurr_[0], jFlux_trCurr_[1], jFlux_trCurr_[2]);
        printf("           Vmolal = %10.3e, vd_Li+ = %10.3e  vd_K+ = %10.3e vd_Cl- = %10.3e\n", Fright_cc_,
               Vdiff_trCurr_[0], Vdiff_trCurr_[1], Vdiff_trCurr_[2]);
      }
    }
#endif

    /*
     * --------------- ADD FLUX TERMS INTO THE RESIDUALS --------------------------------------
     */
#ifdef DEBUG_HKM
    if (iCell == 9 && residType == Base_ResidEval) {
      if (iCell == NumLcCells - 1) {
       // printf("we are here fluxFLeft = %g\n", fluxFleft);
      }
    }
#endif

    /*
     *  Total continuity equation - fluxFright and fluxFleft represent the total mole
     *                              fluxes coming and going from the cell.
     *                    R =   d rho d t + del dot (rho V) = 0
     */
    res[indexCent_EqnStart_BD + EQ_TCont_offset_BD] += (fluxFright - fluxFleft);

    /*
     * Species continuity Equation - 2
     */
    for (int k = 0; k < nsp_ - 2; k++) {
      res[indexCent_EqnStart_BD + EQ_Species_offset_BD + k] += (fluxXright[k] - fluxXleft[k]);
    }

    /*
     *   Current conservation equation
     */
    // res[indexCent_EqnStart_BD + EQ_Current_offset_BD] += (fluxVright - fluxVleft);
    res[indexCent_EqnStart_BD + EQ_Current_offset_BD] += icurrElectrolyte_CBR_[iCell] - icurrElectrolyte_CBL_[iCell];

    /*
     *   Current conservation equation - electrode
     */
    res[indexCent_EqnStart_BD + EQ_Current_offset_ED] += (fluxVElectrodeRight - fluxVElectrodeLeft);

    /*
     *   ------------------- ADD SOURCE TERMS TO THE CURRENT CELL CENTER --------------------------------------
     */

    SetupThermoShop1(&(soln[indexCent_EqnStart_BD]), 0);

    /*
     *  Calculate the electrode reactions
     */
    calcElectrode();

    /*
     * Species 0 Conservation equation
     *    Source terms for the species production rate of Li+.
     */
    res[indexCent_EqnStart_BD + EQ_Species_offset_BD] -= electrodeSpeciesProdRates_[speciesIndex0]
        * surfaceAreaDensity_Cell_[iCell] * xdelCell;

    /*
     * Mole fraction summation equation
     */
    res[indexCent_EqnStart_BD + EQ_MFSum_offset_BD] = 1.0;
    for (int k = 0; k < nsp_; k++) {
      res[indexCent_EqnStart_BD + EQ_MFSum_offset_BD] -= Xcent_cc_[k];
    }

    /*
     * Electroneutrality equation
     */
    res[indexCent_EqnStart_BD + EQ_ChargeBal_offset_BD] = 0.0;
    for (int k = 0; k < nsp_; k++) {
      res[indexCent_EqnStart_BD + EQ_ChargeBal_offset_BD] += Xcent_cc_[k] * spCharge_[k];
    }

    /*
     *  Total continuity equation for moles in the electrolyte phase of the cell
     *     Add in the molar flux from the electrode into the electrolyte phase
     *     We are assuming for the current problem that the volumes stay constant
     */
    res[indexCent_EqnStart_BD + EQ_TCont_offset_BD] -= solnMoleFluxInterface_Curr_ * surfaceAreaDensity_Cell_[iCell]
        * xdelCell;

    /*
     *   Current conservation equation
     *      These are written as a source term. icurrInterface_Curr_ will be positive for the anode
     *      where current flows into the electrolyte and will be negative for the cathode where current
     *      flows out of the cathode
     */
    res[indexCent_EqnStart_BD + EQ_Current_offset_BD] -= icurrInterface_Curr_ * surfaceAreaDensity_Cell_[iCell]
        * xdelCell;

    /*
     *   Current conservation equation for the current in the electrode material
     *      These are written as a sink term They will be exactly opposite to the electrolyte current
     *      source terms
     */
    res[indexCent_EqnStart_BD + EQ_Current_offset_ED] += icurrInterface_Curr_ * surfaceAreaDensity_Cell_[iCell]
        * xdelCell;

    /*
     * Special section if we own the left node of the domain. If we do
     * it will be cell 0. Store currents for later use.
     * These are the correct currents that work for the global balances
     */
    if (IOwnLeft && iCell == 0) {
      if (residType == Base_ResidEval) {
        DiffFluxLeftBound_LastResid_NE[EQ_Current_offset_BD] = res[indexCent_EqnStart_BD + EQ_Current_offset_BD];
        DiffFluxLeftBound_LastResid_NE[EQ_Current_offset_ED] = res[indexCent_EqnStart_BD + EQ_Current_offset_ED];
      }
      if (residType == Base_ShowSolution || residType == Base_ResidEval) {
        icurrElectrode_CBL_[iCell] += +icurrInterface_Curr_ * surfaceAreaDensity_Cell_[iCell] * xdelCell
            + icurrElectrode_CBR_[iCell];

      }
    }
    /*
     * Special section if we own the right node of the domain. If we do
     * it will be cell 0. Store currents for later use.
     * These are the correct currents that work for the global balances
     */
    if (IOwnRight && iCell == (NumLcCells - 1)) {
      if (residType == Base_ResidEval) {
        DiffFluxRightBound_LastResid_NE[EQ_Current_offset_BD] = res[indexCent_EqnStart_BD + EQ_Current_offset_BD];
        DiffFluxRightBound_LastResid_NE[EQ_Current_offset_ED] = res[indexCent_EqnStart_BD + EQ_Current_offset_ED];
      }
      if (residType == Base_ShowSolution || residType == Base_ResidEval) {
        icurrElectrolyte_CBR_[iCell] = +icurrInterface_Curr_ * surfaceAreaDensity_Cell_[iCell] * xdelCell
            + icurrElectrolyte_CBL_[iCell];

      }
    }
    if (residType == Base_ShowSolution) {
      deltaV_Cell_[iCell] = Electrode_->potentialDrop();
      Ess_Cell_[iCell] = Electrode_->openCircuitVoltage(0);
      overpotential_Cell_[iCell] = Electrode_->overpotential(0);
      icurrRxn_Cell_[iCell] = icurrInterface_Curr_ * surfaceAreaDensity_Cell_[iCell] * xdelCell;
      LiFlux_Cell_[iCell] = jFlux_trCurr_[0];
    }

    /*
     *   ------------------ ADD IN THE TIME DEPENDENT TERMS ----------------------------------------------------
     */
    if (doTimeDependentResid) {

#ifdef DEBUG_CELLRES_HKM
      // double concTot_CurrNew = concTot_Curr_;
      // double mfNew = mfElectrolyte_Soln_Curr_[0];
      if (residType == Base_ResidEval) {
        printf(" Cell = %d, Totalflux_Li+_r = %10.3e,  = %10.3e, Totalflux_Li+_l ", iCell, fluxXright[0], fluxXleft[0]);
      }
#endif
      double newStuffTC = concTot_Curr_ * porosity_Curr_ * xdelCell;
      double newStuffSpecies0 = Xcent_cc_[0] * newStuffTC;

      /*
       *   .................... Calculate quantities needed at the previous time step
       */
      /*
       * Setup shop with the old time step
       */
      SetupThermoShop1(&(solnOld[indexCent_EqnStart_BD]), 0);

      double oldStuffTC = concTot_Curr_ * porosity_Curr_ * xdelCell;
      double oldStuffSpecies0 = mfElectrolyte_Soln_Curr_[0] * oldStuffTC;
      double tmp = (newStuffSpecies0 - oldStuffSpecies0) * rdelta_t;

#ifdef DEBUG_CELLRES_HKM
      if (residType == Base_ResidEval) {
        printf(" deltaT term = %10.3e BulkSum = %10.3e\n", tmp, tmp + (fluxXright[0] - fluxXleft[0]));
      }
#endif
      /*
       *   .................... Add these terms in the residual
       */
      /*
       *  Add in the time term for species 0
       */
      res[indexCent_EqnStart_BD + EQ_Species_offset_BD + 0] += tmp;

      /*
       *   Add in the time term for the total continuity equation
       *         note: the current problem will have this term equally zero always.
       *               However, we put it in here for the next problem.
       */
      res[indexCent_EqnStart_BD + EQ_TCont_offset_BD] += (newStuffTC - oldStuffTC) * rdelta_t;

      /*
       *   .................... Go back to setting up shop at the current time
       */
      SetupThermoShop1(&(soln[indexCent_EqnStart_BD]), 0);
    }

  }

}
//=====================================================================================================================
void
porousLiKCl_LiSiAnode_dom1D::calcElectrode()
{
#undef PRECIPITATE
#ifdef PRECIPITATE
  //print something if one of the salts is precipitating
  int iPrec = checkPrecipitation();
  flagPrecipitation = ( flagPrecipitation > iPrec ) ? flagPrecipitation : iPrec;
  if ( iPrec >= 0 ) {
    //    icurrInterface_Curr_ = 0.0;
    //    solnMoleFluxInterface_Curr_ = 0.0;
    //    return;
  }
#endif

  /*
   * Calculate the rates of production of all species in the Electrode
   * and determine the current
   */
  icurrInterface_Curr_ = Electrode_->getNetSurfaceProductionRatesCurrent(0, &electrodeSpeciesProdRates_[0]);
  //  Let's store this for future use in showSolution() to debug the electrode
  icurrInterfacePerSurfaceArea_Cell_[cIndex_cc_] = icurrInterface_Curr_;

  /*
   * Get the phase mole flux
   */
  Electrode_->getPhaseMoleFlux(0, &phaseMoleFlux_[0]);
  int sf = Electrode_->solnPhaseIndex();
  solnMoleFluxInterface_Curr_ = phaseMoleFlux_[sf];

  if (residType_Curr_ == Base_ShowSolution) {
    capacityDischarged_Cell_[cIndex_cc_] = Electrode_->capacityDischarged();
    depthOfDischarge_Cell_[cIndex_cc_] = Electrode_->depthOfDischarge();
    capacityLeft_Cell_[cIndex_cc_] = Electrode_->capacityLeft();
    capacityZeroDoD_Cell_[cIndex_cc_]= Electrode_->capacity();
  }
}
//=====================================================================================================================
void
porousLiKCl_LiSiAnode_dom1D::SetupThermoShop1(const doublereal * const solnElectrolyte_Curr, int type)
{
  if (type == 0) {
    porosity_Curr_ = porosity_Cell_[cIndex_cc_];
  }
  updateElectrolyte(solnElectrolyte_Curr);
  updateElectrode();
}
//=====================================================================================================================
void
porousLiKCl_LiSiAnode_dom1D::SetupThermoShop2(const doublereal * const solnElectrolyte_CurrL,
                                              const doublereal * const solnElectrolyte_CurrR,
                                              int type)
{
  for (int i = 0; i < BDD_.NumEquationsPerNode; i++) {
    solnTemp[i] = 0.5 * (solnElectrolyte_CurrL[i] + solnElectrolyte_CurrR[i]);
  }
  if (type == 0) {
    porosity_Curr_ = 0.5 * (porosity_Cell_[cIndex_cc_ - 1] + porosity_Cell_[cIndex_cc_]);
  } else {
    porosity_Curr_ = 0.5 * (porosity_Cell_[cIndex_cc_ + 1] + porosity_Cell_[cIndex_cc_]);
  }
  updateElectrolyte(&solnTemp[0]);
  updateElectrode();
}
//=====================================================================================================================
// Function updates the ThermoPhase object for the electrolyte
// given the solution vector
/*
 *   Routine will update the molten salt ThermoPhase object with the current state of the electrolyte
 *
 * @param solnElectrolyte
 */
void
porousLiKCl_LiSiAnode_dom1D::updateElectrolyte(const doublereal * const solnElectrolyte_Curr)
{
  /*
   * Get the temperature: Check to see if the temperature is in the solution vector.
   *   If it is not, then use the reference temperature
   */
  temp_Curr_ = TemperatureReference_;
  int iTemp = BDD_.VariableIndexStart_VarName[Temperature];
  if (iTemp >= 0) {
    temp_Curr_ = solnElectrolyte_Curr[iTemp];
  }
  /*
   * Get the pressure
   */
  pres_Curr_ = PressureReference_;

  getMFElectrolyte_soln(solnElectrolyte_Curr);
  getVoltages(solnElectrolyte_Curr);

  ionicLiquid_->setState_TPX(temp_Curr_, pres_Curr_, &mfElectrolyte_Thermo_Curr_[0]);

  ionicLiquid_->setElectricPotential(phiElectrolyte_Curr_);

  // Calculate the total concentration of the electrolyte kmol m-3.
  concTot_Curr_ = ionicLiquid_->molarDensity();

}
//=====================================================================================================================
void
porousLiKCl_LiSiAnode_dom1D::updateElectrode()
{
  /*
   * set the properties in the Electrode object
   *  -> temperature and pressure
   *  -> voltages of the phases
   */
  Electrode_->setState_TP(temp_Curr_, pres_Curr_);
  Electrode_->setVoltages(phiElectrode_Curr_, phiElectrolyte_Curr_);
  Electrode_->setElectrolyteMoleNumbers(&(mfElectrolyte_Thermo_Curr_[0]), true);

  /*
   * Set the internal objects within the electrode
   */
  Electrode_->updateState();
}
//=====================================================================================================================
// Retrieves the voltages from the solution vector and puts them into local storage
/*
 * @param solnElectrolyte start of the solution vector at the current node
 */
void
porousLiKCl_LiSiAnode_dom1D::getVoltages(const double * const solnElectrolyte)
{
  int indexVS = BDD_.VariableIndexStart_VarName[Voltage];
  phiElectrolyte_Curr_ = solnElectrolyte[indexVS];
  phiElectrode_Curr_ = solnElectrolyte[indexVS + 1];
}
//=====================================================================================================================
void
porousLiKCl_LiSiAnode_dom1D::getMFElectrolyte_soln(const double * const solnElectrolyte_Curr)
{
  int indexMF = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
  mfElectrolyte_Soln_Curr_[0] = solnElectrolyte_Curr[indexMF];
  mfElectrolyte_Soln_Curr_[1] = solnElectrolyte_Curr[indexMF + 1];
  mfElectrolyte_Soln_Curr_[2] = solnElectrolyte_Curr[indexMF + 2];
  double mf0 = MAX(mfElectrolyte_Soln_Curr_[0], 0.0);
  double mf1 = MAX(mfElectrolyte_Soln_Curr_[1], 0.0);
  double tmp = mf0 + mf1;

  mfElectrolyte_Thermo_Curr_[0] = (mf0) * 0.5 / tmp;
  mfElectrolyte_Thermo_Curr_[1] = (mf1) * 0.5 / tmp;
  mfElectrolyte_Thermo_Curr_[2] = 0.5;
}
//=====================================================================================================================
void
porousLiKCl_LiSiAnode_dom1D::SetupTranShop(const double xdel, const int type)
{

  /*
   * Determine diffusion velocities
   */

  //set gradients
  gradT_trCurr_ = 0.0;

  if (type == 0) {
    // Left boundary
    gradV_trCurr_ = (Vcent_cc_ - Vleft_cc_) / xdel;
    gradVElectrode_trCurr_ = (VElectrodeCent_cc_ - VElectrodeLeft_cc_) / xdel;

  } else {
    // Right boundary
    gradV_trCurr_ = (Vright_cc_ - Vcent_cc_) / xdel;
    gradVElectrode_trCurr_ = (VElectrodeRight_cc_ - VElectrodeCent_cc_) / xdel;
  }

  if (type == 0) {
    for (int k = 0; k < nsp_; k++) {
      gradX_trCurr_[k] = (Xcent_cc_[k] - Xleft_cc_[k]) / xdel;
    }
  } else {
    for (int k = 0; k < nsp_; k++) {
      gradX_trCurr_[k] = (Xright_cc_[k] - Xcent_cc_[k]) / xdel;
    }
  }

  trans_->getSpeciesVdiffES(1, &gradT_trCurr_, nsp_, &gradX_trCurr_[0], nsp_, &gradV_trCurr_, &Vdiff_trCurr_[0]);

  //  Correct species diffusivities according to the tortuosity
  double bruggemannExp = 2.5;
  Tortuosity tort(bruggemannExp );
  for (int k = 0; k < nsp_; k++) {
    Vdiff_trCurr_[k] *= tort.tortuosityFactor( porosity_Curr_ );
  }
  //Convert from diffusion velocity to diffusion flux
  for (int k = 0; k < nsp_; k++) {
    jFlux_trCurr_[k] = mfElectrolyte_Soln_Curr_[k] * concTot_Curr_ * porosity_Curr_ * Vdiff_trCurr_[k];
  }

  icurrElectrode_trCurr_ = -conductivityElectrode_ * gradVElectrode_trCurr_;
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
porousLiKCl_LiSiAnode_dom1D::saveDomain(Cantera::XML_Node& oNode,
                                        const Epetra_Vector *soln_GLALL_ptr,
                                        const Epetra_Vector *solnDot_GLALL_ptr,
                                        const double t,
                                        bool duplicateOnAllProcs)
{
  // get the NodeVars object pertaining to this global node
  GlobalIndices *gi = LI_ptr_->GI_ptr_;

  // Add a child for this domain
  Cantera::XML_Node& bdom = oNode.addChild("domain");

  // Number of equations per node
  int numEquationsPerNode = BDD_.NumEquationsPerNode;

  // Vector containing the variable names as they appear in the solution vector
  std::vector<VarType> &variableNameList = BDD_.VariableNameList;

  //! First global node of this bulk domain
  int firstGbNode = BDD_.FirstGbNode;

  //! Last Global node of this bulk domain
  int lastGbNode = BDD_.LastGbNode;
  int numNodes = lastGbNode - firstGbNode + 1;

  bdom.addAttribute("id", id());
  bdom.addAttribute("points", numNodes);
  bdom.addAttribute("type", "bulk");
  bdom.addAttribute("numVariables", numEquationsPerNode);

  // Dump out the coordinates
  Cantera::XML_Node& gv = bdom.addChild("grid_data");

  std::vector<double> varContig(numNodes);

  int i = 0;
  for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
    NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
    varContig[i] = nv->x0NodePos();
  }
  ctml::addNamedFloatArray(gv, "X0", varContig.size(), &(varContig[0]), "m", "length");

  for (int iVar = 0; iVar < numEquationsPerNode; iVar++) {
    VarType vt = variableNameList[iVar];
    i = 0;
    std::string nmm = vt.VariableName(200);
    for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
      int istart = nv->EqnStart_GbEqnIndex;
      varContig[i] = (*soln_GLALL_ptr)[istart + ibulk + iVar];
    }
    ctml::addNamedFloatArray(gv, nmm, varContig.size(), &(varContig[0]), "kmol/m3", "concentration");

  }
}
//=====================================================================================================================
static void
drawline(int sp, int ll)
{
  for (int i = 0; i < sp; i++) {
    Cantera::writelog(" ");
  }
  for (int i = 0; i < ll; i++) {
    Cantera::writelog("-");
  }
  Cantera::writelog("\n");
}
//=====================================================================================================================
static void
drawline0(int sp, int ll)
{
  for (int i = 0; i < sp; i++) {
    sprint0(" ");
  }
  for (int i = 0; i < ll; i++) {
    sprint0("-");
  }
  sprint0("\n");
}

//=====================================================================================================================
// Method for writing the header for the surface domain to a tecplot file.
void
porousLiKCl_LiSiAnode_dom1D::writeSolutionTecplotHeader()
{
  int mypid = LI_ptr_->Comm_ptr_->MyPID();
  bool doWrite = !mypid ; //only proc 0 should write

  if (doWrite) {
    
    //open tecplot file
    FILE* ofp;
    string sss = id();
    char filename[20];
    sprintf(filename,"%s%s",sss.c_str(),".dat");
    ofp = fopen( filename, "w");

    //write title and variable list
    fprintf( ofp, "TITLE = \"Solution on Domain %s\"\n",sss.c_str() );

    // Number of equations per node
    int numVar = BDD_.NumEquationsPerNode;
    // Vector containing the variable names as they appear in the solution vector
    std::vector<VarType> &variableNameList = BDD_.VariableNameList;
    //! First global node of this bulk domain

    fprintf( ofp, "VARIABLES = ");
    fprintf( ofp, "\"x [m]\"  \n" );

    for (int k = 0; k < numVar; k++) {
      VarType &vt = variableNameList[k];
      string name = vt.VariableName(15);
      fprintf( ofp, "\"%s\" \t", name.c_str() );
    }
    fprintf( ofp, "\n" );

    /////////////////////////////////////////
    //// end BulkDomain1D section
    /////////////////////////////////////////

    fprintf(ofp, "\"Potential Drop Electrode-Electrolyte [V]\" \t" );
    fprintf(ofp, "\"Equilibrium Potential Drop Electrode-Electrolyte [V]\" \t" );
    fprintf(ofp, "\"Overpotential [V]\" \t" );
    fprintf(ofp, "\"Current Source [A/m^2]\" \t" );
    fprintf(ofp, "\"Li+ flux [kmol/m^2/s]\" \t" );
    fprintf(ofp, "\"Electrolyte Current [A/m^2]\" \t" );
    fprintf(ofp, "\"Electrode Current [A/m^2]\" \t" );
    fprintf(ofp, "\n");
    fclose(ofp);
  }
}
//=====================================================================================================================
// Method for writing the solution on the surface domain to a tecplot file.
/*
 *
 * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
 * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
 * @param t                    time
 *
 * @param duplicateOnAllProcs  If this is true, all processors will include
 *                             the same log information as proc 0. If
 *                             false, the loginfo will only exist on proc 0.
 */
void
porousLiKCl_LiSiAnode_dom1D::writeSolutionTecplot(const Epetra_Vector *soln_GlAll_ptr,
						const Epetra_Vector *solnDot_GlAll_ptr,
						const double t )
{
  int mypid = LI_ptr_->Comm_ptr_->MyPID();
  bool doWrite = !mypid ; //only proc 0 should write

  if (doWrite) {
    
    // get the NodeVars object pertaining to this global node
    GlobalIndices *gi = LI_ptr_->GI_ptr_;
    // Number of equations per node
    int numVar = BDD_.NumEquationsPerNode;
    //! First global node of this bulk domain
    int firstGbNode = BDD_.FirstGbNode;
    //! Last Global node of this bulk domain
    int lastGbNode = BDD_.LastGbNode;
    int numNodes = lastGbNode - firstGbNode + 1;
    
    //open tecplot file
    FILE* ofp;
    string sss = id();
    char filename[20];
    sprintf(filename,"%s%s",sss.c_str(),".dat");
    ofp = fopen( filename, "a");
    
    fprintf(ofp, "ZONE T = \"t = %g [s]\" I = %d SOLUTIONTIME = %19.13E\n", t, numNodes, t);

    for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++) {
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];

      //x-position
      fprintf( ofp, "%g \t", nv->xNodePos() );

      for (int iVar = 0; iVar < numVar; iVar++) {

	//other variables
	int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
	int istart = nv->EqnStart_GbEqnIndex;
	fprintf( ofp, "%g \t", (*soln_GlAll_ptr)[istart + ibulk + iVar] );	
      }
      fprintf( ofp, "\n");

    /////////////////////////////////////////
    //// end BulkDomain1D section
    /////////////////////////////////////////

      int iCell = iGbNode - firstGbNode;
      // surface reaction data
      fprintf(ofp, "%g \t", deltaV_Cell_[iCell]);
      fprintf(ofp, "%g \t", Ess_Cell_[iCell]);
      fprintf(ofp, "%g \t", overpotential_Cell_[iCell]); 
      fprintf(ofp, "%g \t", icurrRxn_Cell_[iCell]); 
      fprintf(ofp, "%g \t", LiFlux_Cell_[iCell]); 

      // current flux -- average of left and right fluxes
      fprintf(ofp, "%g \t", 0.5 * ( icurrElectrolyte_CBL_[iCell] + icurrElectrolyte_CBR_[iCell] ) );
      fprintf(ofp, "%g \t", 0.5 * ( icurrElectrode_CBL_[iCell] + icurrElectrode_CBR_[iCell] ) );
      fprintf(ofp, "\n");
    }

    fclose(ofp);

#undef DAKOTAOUT
#ifdef  DAKOTAOUT
    double firstOutputTime = 0.1;
    static double outputTime = 0.1;
    if ( t > outputTime ) {
      std::ofstream dfp;
      if ( outputTime > firstOutputTime )
	dfp.open( "results.out", std::ios_base::app );
      else
	dfp.open( "results.out", std::ios_base::out );
      int indexVS = BDD_.VariableIndexStart_VarName[Voltage];
      double potentialE_left = (*soln_GlAll_ptr)[indexVS];
      double lithium_left = (*soln_GlAll_ptr)[1];
      
      NodalVars *nv = gi->NodalVars_GbNode[lastGbNode-1];
      int istart = nv->EqnStart_GbEqnIndex;
      double potentialE_right = (*soln_GlAll_ptr)[istart+indexVS];
      double lithium_right = (*soln_GlAll_ptr)[istart+1];
      double ohmicLoss =  potentialE_left - potentialE_right;
      double moleFracDelta = lithium_left - lithium_right;

      // uncomment next line to get ohmic loss in electrolyte
      //dfp << ohmicLoss << "  t" << outputTime << std::endl;

      // uncomment next line to get change in litium mole fraction
      //dfp << moleFracDelta << "  t" << outputTime << std::endl;
      
      double Li_eutectic = 0.585/2.0;
      double Rgas_local = 8.3144;
      double Faraday_local = 9.6485e4;
      double concentrationPotential = Rgas_local * TemperatureReference_ / Faraday_local
	* log( lithium_left / lithium_right ) ;

      //uncomment following line to get concentration potential loss
      dfp << concentrationPotential << "  t" << outputTime << std::endl;

      //uncomment following line to get combined ohmic and concentration losses
      outputTime *= 10.0;
    }
      
#endif
#undef DAKOTAOUT
  }
}



//=====================================================================================================================
// Base class for writing the solution on the domain to a logfile.
/*
 *
 * @param soln_GlALL_ptr       Pointer to the Global-All solution vector
 * @param solnDot_GlALL_ptr    Pointer to the Global-All solution dot vector
 * @param soln_ptr             Pointer to the solution vector
 * @param solnDot_ptr          Pointer to the time derivative of the solution vector
 * @param solnOld_ptr          Pointer to the solution vector at the old time step
 * @param residInternal _ptr   Pointer to the current value of the residual just calculated
 *                             by a special call to the residEval()
 * @param t                    time
 * @param rdelta_t             The inverse of the value of delta_t
 * @param indentSpaces         Indentation that all output should have as a starter
 * @param duplicateOnAllProcs  If this is true, all processors will include
 *                             the same log information as proc 0. If
 *                             false, the loginfo will only exist on proc 0.
 */
void
porousLiKCl_LiSiAnode_dom1D::showSolution(const Epetra_Vector *soln_GlAll_ptr,
                                          const Epetra_Vector *solnDot_GlAll_ptr,
                                          const Epetra_Vector *soln_ptr,
                                          const Epetra_Vector *solnDot_ptr,
                                          const Epetra_Vector *solnOld_ptr,
                                          const Epetra_Vector_Owned *residInternal_ptr,
                                          const double t,
                                          const double rdelta_t,
                                          int indentSpaces,
                                          bool duplicateOnAllProcs)
{

  // nn is the number of block rows in the printout
  int nn = NumDomainEqns / 5;
  int mypid = LI_ptr_->Comm_ptr_->MyPID();
  bool do0Write = (!mypid || duplicateOnAllProcs);

  std::vector<VarType> &variableNameList = BDD_.VariableNameList;
  int iBlock;
  int iGbNode;
  int n;
  int nPoints = BDD_.LastGbNode - BDD_.FirstGbNode + 1;

  std::string indent = "";
  for (int i = 0; i < indentSpaces; i++) {
    indent += " ";
  }
  const char *ind = indent.c_str();
  doublereal v;
  GlobalIndices *gi = LI_ptr_->GI_ptr_;
  // Number of points in each vector
  string sss = id();
  stream0 ss;

  if (do0Write) {
    drawline0(indentSpaces, 80);
    ss.print0("%s  Solution on Bulk Domain %12s : Number of variables = %d\n", ind, sss.c_str(), NumDomainEqns);
    ss.print0("%s                                         : Number of Nodes = %d\n", ind, nPoints);
    ss.print0("%s                                         : Beginning pos %g\n", ind, BDD_.Xpos_start);
    ss.print0("%s                                         : Ending    pos %g\n", ind, BDD_.Xpos_end);
  }
  if (do0Write) {
    for (iBlock = 0; iBlock < nn; iBlock++) {
      drawline0(indentSpaces, 80);
      ss.print0("%s        z   ", ind);
      for (n = 0; n < 5; n++) {
        int ivar = iBlock * 5 + n;
        VarType vt = variableNameList[ivar];
        string name = vt.VariableName(15);
        ss.print0(" %15s", name.c_str());
      }
      ss.print0("\n");
      drawline(indentSpaces, 80);

      for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {

        NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
        doublereal x = nv->xNodePos();
        ss.print0("\n%s    %-10.4E ", ind, x);
        int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
        int istart = nv->EqnStart_GbEqnIndex;
        for (n = 0; n < 5; n++) {
          int ivar = iBlock * 5 + n;
          VarType vt = variableNameList[ivar];
          v = (*soln_GlAll_ptr)[istart + ibulk + iBlock * 5 + n];
          ss.print0(" %-10.4E ", v);
        }
      }
      ss.print0("\n");
    }

    int nrem = NumDomainEqns - 5 * nn;
    if (nrem > 0) {
      drawline(indentSpaces, 80);
      ss.print0("%s        z   ", ind);
      for (n = 0; n < nrem; n++) {
        int ivar = nn * 5 + n;
        VarType vt = variableNameList[ivar];
        string name = vt.VariableName(15);
        ss.print0(" %15s", name.c_str());
      }
      ss.print0("\n");
      drawline0(indentSpaces, 80);

      for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
        NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
        doublereal x = nv->xNodePos();
        ss.print0("%s    %-10.4E ", ind, x);
        int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
        int istart = nv->EqnStart_GbEqnIndex;
        for (n = 0; n < nrem; n++) {
          int ivar = iBlock * 5 + n;
          VarType vt = variableNameList[ivar];
          v = (*soln_GlAll_ptr)[istart + ibulk + nn * 5 + n];
          ss.print0(" %-10.4E ", v);
        }
        ss.print0("\n");
      }
    }
    drawline(indentSpaces, 80);
  }
  if (do0Write) {
    // ----------------------------------------------------
    // --             PRINT REACTION RATES within the cell --
    // ----------------------------------------------------
    ss.print0("%s        z       Delta_V        Ess    Overpotential icurrRxnCell", ind);
    ss.print0("\n");
    drawline0(indentSpaces, 80);
  }
  doublereal x;
  int iCell;
  for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
      iCell = iGbNode - BDD_.FirstGbNode;
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      x = nv->xNodePos();
      ss.print0("%s    %-10.4E ", ind, x);
      ss.print0("%11.4E ", deltaV_Cell_[iCell]);
      ss.print0("%11.4E ", Ess_Cell_[iCell]);
      ss.print0("%11.4E ", overpotential_Cell_[iCell]);
      ss.print0("%11.4E ", icurrRxn_Cell_[iCell]);
      ss.print0("\n");
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
  }

  if (do0Write) {
    // -----------------------------------------------------------------------------------------------------------------
    // --             PRINT VOLUME DETAILS --
    // -----------------------------------------------------------------------------------------------------------------
    ss.print0("\n");
    drawline0(indentSpaces, 80);
    ss.print0("%s        z      Porosity", ind);
    ss.print0("\n");
    drawline0(indentSpaces, 80);
  }
  for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
      iCell = iGbNode - BDD_.FirstGbNode;
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      x = nv->xNodePos();
      ss.print0("%s    %-10.4E ", ind, x);
      ss.print0("%11.4E ", porosity_Cell_[iCell]);
      ss.print0("\n");
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
  }

  if (do0Write) { 
    ss.print0("\n");
    drawline0(indentSpaces, 80);
    // -----------------------------------------------------------------------------------------------------------------
    // --             PRINT SURFACE REACTION DETAILS ABOUT EACH CELL --
    // -----------------------------------------------------------------------------------------------------------------
    ss.print0("%s        z      SurfaceArea  currPerSA   DeltaZ", ind);
    ss.print0("\n");
    drawline0(indentSpaces, 80);
  }
  for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
      iCell = iGbNode - BDD_.FirstGbNode;
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      x = nv->xNodePos();
      ss.print0("%s    %-10.4E ", ind, x);
      ss.print0("%11.4E ", surfaceAreaDensity_Cell_[iCell]);
      ss.print0("%11.4E ", icurrInterfacePerSurfaceArea_Cell_[iCell]);
      ss.print0("%11.4E ", xdelCell_Cell_[iCell]);
      ss.print0("\n");
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
  }


  if (do0Write) {
    // -----------------------------------------------------------------------------------------------------------------
    // --             PRINT DEPTH OF DISCHARGE VALUES FOR EACH CELL --
    // -----------------------------------------------------------------------------------------------------------------
    ss.print0("\n");
    drawline0(indentSpaces, 80);
    ss.print0("%s        z    capDischarged DepthDischge  capLeft  capZeroDOD  ", ind);
    ss.print0("\n");
    drawline0(indentSpaces, 80);
  }
  for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
      iCell = iGbNode - BDD_.FirstGbNode;
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      x = nv->xNodePos();
      ss.print0("%s    %-10.4E ", ind, x);
      ss.print0("%11.4E ", capacityDischarged_Cell_[iCell]);
      ss.print0("%11.4E ", depthOfDischarge_Cell_[iCell]);
      ss.print0("%11.4E ", capacityLeft_Cell_[iCell]);
      ss.print0("%11.4E ", capacityZeroDoD_Cell_[iCell]);
  
      ss.print0("\n");
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
  }
  

  if (do0Write) {
    // ----------------------------------------------------
    // --             PRINT FLUXES AT THE CELL BOUNDARIES --
    // ----------------------------------------------------
    ss.print0("\n");
    drawline0(indentSpaces, 80);
    ss.print0("%s    CellBound    z     IcurrElectrolyte IcurrElectrode ", ind);
    ss.print0("\n");
    drawline(indentSpaces, 80);
  }

  NodalVars *nvl;
  NodalVars *nvr;

  for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {

      ss.print0("%s    ", ind);
      if (iGbNode == BDD_.FirstGbNode) {

        iCell = 0;
        nvr = gi->NodalVars_GbNode[BDD_.FirstGbNode];
        x = nvr->xNodePos();
        ss.print0("Lft-0     %11.4E ", x);
        ss.print0("%11.4E ", icurrElectrolyte_CBL_[0]);
        ss.print0("%11.4E \n", icurrElectrode_CBL_[0]);
        ss.print0("%s    ", ind);
        iCell = iGbNode - BDD_.FirstGbNode;
        nvl = gi->NodalVars_GbNode[iGbNode];
        nvr = gi->NodalVars_GbNode[iGbNode + 1];
        x = 0.5 * (nvl->xNodePos() + nvr->xNodePos());
        ss.print0("%3d-%-3d   %11.4E ", iCell, iCell + 1, x);
        ss.print0("%11.4E ", icurrElectrolyte_CBR_[iCell]);
        ss.print0("%11.4E ", icurrElectrode_CBR_[iCell]);
      } else if (iGbNode == BDD_.LastGbNode) {
        iCell = BDD_.LastGbNode - BDD_.FirstGbNode;
        nvr = gi->NodalVars_GbNode[BDD_.LastGbNode];
        x = nvr->xNodePos();
        ss.print0("%3d-Rgt   %11.4E ", iCell, x);
        ss.print0("%11.4E ", icurrElectrolyte_CBR_[iCell]);
        ss.print0("%11.4E ", icurrElectrode_CBR_[iCell]);
      } else {
        iCell = iGbNode - BDD_.FirstGbNode;
        nvl = gi->NodalVars_GbNode[iGbNode];
        nvr = gi->NodalVars_GbNode[iGbNode + 1];
        x = 0.5 * (nvl->xNodePos() + nvr->xNodePos());
        ss.print0("%3d-%-3d   %11.4E ", iCell, iCell + 1, x);
        ss.print0("%11.4E ", icurrElectrolyte_CBR_[iCell]);
        ss.print0("%11.4E ", icurrElectrode_CBR_[iCell]);
      }
      ss.print0("\n");
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
  }

}
//=====================================================================================================================
// Generate the initial conditions
/*
 *   The basic algorithm is to loop over the volume domains.
 *   Then, we loop over the surface domains
 *
 * @param doTimeDependentResid    Boolean indicating whether we should
 *                                formulate the time dependent residual
 * @param soln                    Solution vector. This is the input to
 *                                the residual calculation.
 * @param solnDot                 Solution vector. This is the input to
 *                                the residual calculation.
 * @param t                       Time
 * @param delta_t                 delta_t for the initial time step
 */
void
porousLiKCl_LiSiAnode_dom1D::initialConditions(const bool doTimeDependentResid,
                                               Epetra_Vector *soln_ptr,
                                               Epetra_Vector *solnDot,
                                               const double t,
                                               const double delta_t)
{
  Epetra_Vector &soln = *soln_ptr;
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    cIndex_cc_ = iCell;

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
    /*
     * Offsets for the variable unknowns in the solution vector for the electrolyte domain
     */
    int iVAR_Vaxial_BD = BDD_.VariableIndexStart_VarName[Velocity_Axial];
    int iVar_Species_BD = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
    int iVar_Voltage_BD = BDD_.VariableIndexStart_VarName[Voltage];
    int iVar_Voltage_ED = iVar_Voltage_BD + 1;

    soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = 0.0;

    /*
     * Get initial mole fractions from PSinput
     */
    int iLiCl = PSinput.PhaseList_->globalSpeciesIndex ("LiCl(L)");
    int iKCl = PSinput.PhaseList_->globalSpeciesIndex ("KCl(L)");
    int iLi = PSinput.PhaseList_->globalSpeciesIndex ("Li+");
    int iK = PSinput.PhaseList_->globalSpeciesIndex ("K+");

    soln[indexCent_EqnStart_BD + iVar_Species_BD + 0] = 
      PSinput.electrolyteMoleFracs_[iLiCl] / 2.0 +
      PSinput.electrolyteMoleFracs_[iLi];
    soln[indexCent_EqnStart_BD + iVar_Species_BD + 1] = 
      PSinput.electrolyteMoleFracs_[iKCl] / 2.0 +
      PSinput.electrolyteMoleFracs_[iK];
    soln[indexCent_EqnStart_BD + iVar_Species_BD + 2] = 0.5;

    soln[indexCent_EqnStart_BD + iVar_Voltage_BD] = -0.07;
    soln[indexCent_EqnStart_BD + iVar_Voltage_ED] = 0.0;

  }
}
//=====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void
porousLiKCl_LiSiAnode_dom1D::setAtolVector(double atolDefault, const Epetra_Vector_Ghosted & soln, 
					   Epetra_Vector_Ghosted & atolVector,
					   const Epetra_Vector_Ghosted * const atolV)
{
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    cIndex_cc_ = iCell;

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
    /*
     * Offsets for the variable unknowns in the solution vector for the electrolyte domain
     */
    int iVAR_Vaxial_BD = BDD_.VariableIndexStart_VarName[Velocity_Axial];
    int iVar_Species_BD = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
    int iVar_Voltage_BD = BDD_.VariableIndexStart_VarName[Voltage];
    int iVar_Voltage_ED = iVar_Voltage_BD + 1;

    /*
     * Set the atol value for the axial velocity
     *   arithmetically scaled -> so this is a characteristic value
     */
    double vax = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];
    atolVector[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = MAX(1.0E-4, 1.0E-1 * vax);

    /*
     * Set atol values for the species mole fractions
     */
    double val = atolDefault * 1.0E-2;
    if (val < 1.0E-12) {
      val = 1.0E-12;
    }
    for (int k = 0; k < nsp_; k++) {
      atolVector[indexCent_EqnStart_BD + iVar_Species_BD + k] = val;
    }

    /*
     * Set the atol value for the electrolyte voltage
     *      arithmetically scaled.-> so this is a characteristic value
     *         1 kcal gmol-1 = 0.05 volts
     */
    atolVector[indexCent_EqnStart_BD + iVar_Voltage_BD] = 0.05;
    atolVector[indexCent_EqnStart_BD + iVar_Voltage_ED] = 0.05;
  }
}
//=====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void
porousLiKCl_LiSiAnode_dom1D::setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted & soln, 
						   const Epetra_Vector_Ghosted & solnDot,
						   Epetra_Vector_Ghosted & atolVector,
						   const Epetra_Vector_Ghosted * const atolV )
{
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    cIndex_cc_ = iCell;

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
    /*
     * Offsets for the variable unknowns in the solution vector for the electrolyte domain
     */
    int iVAR_Vaxial_BD = BDD_.VariableIndexStart_VarName[Velocity_Axial];
    int iVar_Species_BD = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
    int iVar_Voltage_BD = BDD_.VariableIndexStart_VarName[Voltage];
    int iVar_Voltage_ED = iVar_Voltage_BD + 1;

    /*
     * Set the atol value for the axial velocity
     *   arithmetically scaled -> so this is a characteristic value
     */
    double vax = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];
    atolVector[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = MAX(1.0E-4, 1.0E-1 * vax);

    /*
     * Set atol values for the species mole fractions
     */
    double val = atolDefault * 1.0E6;
    if (val < 1.0E-5) {
      val = 1.0E-5;
    }
    for (int k = 0; k < nsp_; k++) {
      atolVector[indexCent_EqnStart_BD + iVar_Species_BD + k] = val;
    }

    /*
     * Set the atol value for the electrolyte voltage
     *      arithmetically scaled.-> so this is a characteristic value
     *         1 kcal gmol-1 = 0.05 volts
     */
    atolVector[indexCent_EqnStart_BD + iVar_Voltage_BD] = 0.05;
    atolVector[indexCent_EqnStart_BD + iVar_Voltage_ED] = 0.05;
  }
}
//=====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void
porousLiKCl_LiSiAnode_dom1D::setAtolDeltaDamping(double atolDefault, double relcoeff, 
						 const Epetra_Vector_Ghosted & soln, 
						 Epetra_Vector_Ghosted & atolDeltaDamping,
						 const Epetra_Vector_Ghosted * const atolV)
{
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    cIndex_cc_ = iCell;

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
    /*
     * Offsets for the variable unknowns in the solution vector for the electrolyte domain
     */
    int iVAR_Vaxial_BD = BDD_.VariableIndexStart_VarName[Velocity_Axial];
    int iVar_Species_BD = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
    int iVar_Voltage_BD = BDD_.VariableIndexStart_VarName[Voltage];
    int iVar_Voltage_ED = iVar_Voltage_BD + 1;

    /*
     * Set the atol value for the axial velocity
     *   arithmetically scaled -> so this is a characteristic value
     */
    double vax = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];
    atolDeltaDamping[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = MAX(1.0E-4, 1.0E-1 * vax) * relcoeff;

    /*
     * Set atol values for the species mole fractions
     */
    double val = atolDefault * 1.0E6;
    if (val < 1.0E-4) {
      val = 1.0E-4;
    }
    for (int k = 0; k < nsp_; k++) {
      atolDeltaDamping[indexCent_EqnStart_BD + iVar_Species_BD + k] = val * relcoeff;
    }

    /*
     * Set the atol value for the electrolyte voltage
     *      arithmetically scaled.-> so this is a characteristic value
     *         1 kcal gmol-1 = 0.05 volts
     */
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Voltage_BD] = 0.05 * relcoeff;
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Voltage_ED] = 0.05 * relcoeff;
  }
}
//=====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void
porousLiKCl_LiSiAnode_dom1D::setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff, 
							 const Epetra_Vector_Ghosted & soln, 
							 const Epetra_Vector_Ghosted & solnDot, 
							 Epetra_Vector_Ghosted & atolDeltaDamping,
							 const Epetra_Vector_Ghosted * const atolV)
{
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    cIndex_cc_ = iCell;

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
    /*
     * Offsets for the variable unknowns in the solution vector for the electrolyte domain
     */
    int iVAR_Vaxial_BD = BDD_.VariableIndexStart_VarName[Velocity_Axial];
    int iVar_Species_BD = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
    int iVar_Voltage_BD = BDD_.VariableIndexStart_VarName[Voltage];
    int iVar_Voltage_ED = iVar_Voltage_BD + 1;

    /*
     * Set the atol value for the axial velocity
     *   arithmetically scaled -> so this is a characteristic value
     */
    double vax = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];
    atolDeltaDamping[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = MAX(1.0E-4, 1.0E-1 * vax) * relcoeff;

    /*
     * Set atol values for the species mole fractions
     */
    double val = atolDefault * 1.0E6;
    if (val < 1.0E-1) {
      val = 1.0E-1;
    }
    for (int k = 0; k < nsp_; k++) {
      atolDeltaDamping[indexCent_EqnStart_BD + iVar_Species_BD + k] = val * relcoeff;
    }

    /*
     * Set the atol value for the electrolyte voltage
     *      arithmetically scaled.-> so this is a characteristic value
     *         1 kcal gmol-1 = 0.05 volts
     */
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Voltage_BD] = 0.05 * relcoeff;
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Voltage_ED] = 0.05 * relcoeff;
  }
}

//=====================================================================================================================
/**
 * Method to check for precipitation of the salts.  
 * Returns index of offending cation or -1 if no precipitation
 */
int 
porousLiKCl_LiSiAnode_dom1D::checkPrecipitation(  ) {

  // molten salt phase
  string id_salt = "LiKCl_Margules";
  int iph = (PSinput.PhaseList_)->globalPhaseIndex(id_salt);
  if (iph < 0) {
    throw CanteraError("porousLiKCl_LiSiAnode_dom1D::checkPrecipitation()",
                       "Can't find the phase in the phase list: " + id_salt);
  }
  ThermoPhase* tmpPhase = & (PSinput.PhaseList_)->thermo(iph);
  
  MargulesVPSSTP *salt ;
  salt = dynamic_cast<Cantera::MargulesVPSSTP *>( tmpPhase->duplMyselfAsThermoPhase() );
  
  int iKCl_l = salt->speciesIndex("KCl(L)");
  int iLiCl_l = salt->speciesIndex("LiCl(L)");
  int iK_ion = ionicLiquid_->speciesIndex("K+");
  int iLi_ion = ionicLiquid_->speciesIndex("Li+");

  //solid LiCl phase
  id_salt = "LiCl(S)";
  iph = (PSinput.PhaseList_)->globalPhaseIndex(id_salt);
  if (iph < 0) {
    throw CanteraError("porousLiKCl_LiSiAnode_dom1D::checkPrecipitation()",
                       "Can't find the phase in the phase list: " + id_salt);
  }
  tmpPhase = & (PSinput.PhaseList_)->thermo(iph);
  Cantera::ThermoPhase *LiCl_solid = tmpPhase->duplMyselfAsThermoPhase() ;

  //solid KCl phase
  id_salt = "KCl(S)";
  iph = (PSinput.PhaseList_)->globalPhaseIndex(id_salt);
  if (iph < 0) {
    throw CanteraError("porousLiKCl_LiSiAnode_dom1D::checkPrecipitation()",
                       "Can't find the phase in the phase list: " + id_salt);
  }
  tmpPhase = & (PSinput.PhaseList_)->thermo(iph);
  Cantera::ThermoPhase *KCl_solid = tmpPhase->duplMyselfAsThermoPhase() ;



  //set current states
  //mole fraction of the electrolyte ions are held in
  // mfElectrolyte_Thermo_Curr_[3]
  double x[2];
  x[iKCl_l] = 2.0 * mfElectrolyte_Thermo_Curr_[iK_ion];
  x[iLiCl_l] = 2.0 * mfElectrolyte_Thermo_Curr_[iLi_ion];
  salt->setState_TPX(temp_Curr_, pres_Curr_, x);
  
  /*  
  string f_licl = "LiCl_solid.xml";
  string id = "LiCl(S)";
  Cantera::ThermoPhase *LiCl_solid = Cantera::newPhase(f_licl, id);
  */
  LiCl_solid->setState_TP(temp_Curr_, pres_Curr_);
  
  /*  
  string f_kcl = "KCl_solid.xml";
  id = "KCl(S)";
  Cantera::ThermoPhase *KCl_solid = Cantera::newPhase(f_kcl, id);
  */
  KCl_solid->setState_TP(temp_Curr_, pres_Curr_);
  
  //get chemical potentials
  double mu[3];
  //molten salt
  salt->getChemPotentials(mu);
  double mu_LiCl_liq = mu[iLiCl_l];
  double mu_KCl_liq = mu[iKCl_l];
  //solid LiCl
  LiCl_solid->getChemPotentials(mu);
  double mu_LiCl_solid = mu[0];
  //solid KCl  
  KCl_solid->getChemPotentials(mu);
  double mu_KCl_solid = mu[0];
  
  if ( mu_KCl_solid < mu_KCl_liq ) {
    std::cout << "KCl is precipitating at mole fraction " 
	      << 2.0 * mfElectrolyte_Thermo_Curr_[iK_ion] << std::endl;
    return iK_ion;
  }
  if ( mu_LiCl_solid < mu_LiCl_liq ) {
    std::cout << "LiCl is precipitating at mole fraction " 
	      << 2.0 * mfElectrolyte_Thermo_Curr_[iLi_ion] << std::endl;
    return iLi_ion;
  }
  
  //no precipitation
  return -1;
  
}

//=====================================================================================================================
void
porousLiKCl_LiSiAnode_dom1D::err(const char *msg)
{
  printf("porousLiKCl_LiSiAnode_dom1D: function not implemented: %s\n", msg);
  exit(-1);
}
//=====================================================================================================================
} //namespace m1d
//=====================================================================================================================


