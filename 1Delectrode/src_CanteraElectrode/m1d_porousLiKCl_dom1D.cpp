/*
 * m1d_porousLiKCl_dom1D.cpp
 *
 */

/*
 *  $Id: m1d_porousLiKCl_dom1D.cpp 504 2013-01-07 22:32:48Z hkmoffa $
 */

//  This is a heavyweight base class that provides the function
//  evaluation for a single bulk domain.
#include "m1d_porousLiKCl_dom1D.h"
#include "m1d_BDT_porousLiKCl.h"

#include "m1d_NodalVars.h"
#include "m1d_LocalNodeIndices.h"

#include "m1d_GlobalIndices.h"
#include "m1d_exception.h"
#include "m1d_DomainLayout.h"
#include "m1d_Comm.h"

#include "cantera/base/ctml.h"
#include "cantera/transport/Tortuosity.h"

//next two lines added for salt precipitation
#include "cantera/thermo.h"
#include "cantera/thermo/MargulesVPSSTP.h"
#include "m1d_ProblemStatementCell.h"
#include "m1d_CanteraElectrodeGlobals.h"

using namespace std;
using namespace Cantera;


namespace m1d
{

//=====================================================================================================================
porousLiKCl_dom1D::porousLiKCl_dom1D(BulkDomainDescription & bdd) :
  porousFlow_dom1D(bdd),
  ionicLiquid_(0), trans_(0), nph_(0), nsp_(0), concTot_cent_(0.0), concTot_cent_old_(0.0),
      concTot_Cell_(0), concTot_Cell_old_(0), cIndex_cc_(0), Fleft_cc_(0.0),
      Fright_cc_(0.0), Vleft_cc_(0.0), Vcent_cc_(0.0), Vright_cc_(0.0), Xleft_cc_(0), Xcent_cc_(0), Xright_cc_(0),
      spCharge_(0), 
  mfElectrolyte_Soln_Curr_(0), 
  mfElectrolyte_Thermo_Curr_(0), 
  mfElectrolyte_Soln_Cell_old_(0,0),
  mfElectrolyte_SolnDot_Curr_(0),
  pmVolElectrolyte_Curr_(0),
  concTotDot_Curr_(0.0),
      gradT_trCurr_(0.0), gradV_trCurr_(0), gradX_trCurr_(0), Vdiff_trCurr_(0), jFlux_trCurr_(0),
      icurrElectrolyte_CBL_(0), icurrElectrolyte_CBR_(0), solnTemp(0), ivb_(VB_MOLEAVG)
{

  BDT_porousLiKCl *fa = dynamic_cast<BDT_porousLiKCl *> (&bdd);
  if (!fa) {
    throw m1d_Error("confused", "confused");
  }
  ionicLiquid_ = fa->ionicLiquid_;
  trans_ = fa->trans_;
  nsp_ = 3;
  nph_ = 1;

}
//=====================================================================================================================
porousLiKCl_dom1D::porousLiKCl_dom1D(const porousLiKCl_dom1D &r) :
  porousFlow_dom1D(r.BDD_),
  ionicLiquid_(0), trans_(0), nph_(0), nsp_(0), concTot_cent_(0.0), concTot_cent_old_(0.0),
      concTot_Cell_(0), concTot_Cell_old_(0), cIndex_cc_(0), Fleft_cc_(0.0),
      Fright_cc_(0.0), Vleft_cc_(0.0), Vcent_cc_(0.0), Vright_cc_(0.0), Xleft_cc_(0), Xcent_cc_(0), Xright_cc_(0),
      spCharge_(0),
  mfElectrolyte_Soln_Curr_(0), 
  mfElectrolyte_Thermo_Curr_(0),
  mfElectrolyte_Soln_Cell_old_(0,0),
  mfElectrolyte_SolnDot_Curr_(0),
  pmVolElectrolyte_Curr_(0),
  concTotDot_Curr_(0.0),
      gradT_trCurr_(0.0), gradV_trCurr_(0), gradX_trCurr_(0), Vdiff_trCurr_(0), jFlux_trCurr_(0),
      icurrElectrolyte_CBL_(0), icurrElectrolyte_CBR_(0), solnTemp(0), ivb_(VB_MOLEAVG)
{
  porousLiKCl_dom1D::operator=(r);
}
//=====================================================================================================================
porousLiKCl_dom1D::~porousLiKCl_dom1D()
{
}
//=====================================================================================================================
porousLiKCl_dom1D &
porousLiKCl_dom1D::operator=(const porousLiKCl_dom1D &r)
{
  if (this == &r) {
    return *this;
  }
  // Call the parent assignment operator
  porousFlow_dom1D::operator=(r);

  ionicLiquid_ = r.ionicLiquid_;
  trans_ = r.trans_;
  nph_ = r.nph_;
  nsp_ = r.nsp_;
  concTot_cent_ = r.concTot_cent_;
  concTot_cent_old_ = r.concTot_cent_old_;
  concTot_Cell_ = r.concTot_Cell_;
  concTot_Cell_old_ = r.concTot_Cell_old_;
  cIndex_cc_ = r.cIndex_cc_;
  Fleft_cc_ = r.Fleft_cc_;
  Fright_cc_ = r.Fright_cc_;
  Vleft_cc_ = r.Vleft_cc_;
  Vcent_cc_ = r.Vcent_cc_;
  Vright_cc_ = r.Vright_cc_;
  Xleft_cc_ = r.Xleft_cc_;
  Xcent_cc_ = r.Xcent_cc_;
  Xright_cc_ = r.Xright_cc_;
  spCharge_ = r.spCharge_;
  mfElectrolyte_Soln_Curr_ = r.mfElectrolyte_Soln_Curr_;
  mfElectrolyte_Thermo_Curr_ = r.mfElectrolyte_Thermo_Curr_;
  mfElectrolyte_Soln_Cell_old_ = r.mfElectrolyte_Soln_Cell_old_;
  mfElectrolyte_SolnDot_Curr_ = r.mfElectrolyte_SolnDot_Curr_;
  pmVolElectrolyte_Curr_ = r.pmVolElectrolyte_Curr_;
  concTotDot_Curr_ = r.concTotDot_Curr_;
  gradT_trCurr_ = r.gradT_trCurr_;
  gradV_trCurr_ = r.gradV_trCurr_;
  gradX_trCurr_ = r.gradX_trCurr_;
  Vdiff_trCurr_ = r.Vdiff_trCurr_;
  jFlux_trCurr_ = r.jFlux_trCurr_;
  icurrElectrolyte_CBL_ = r.icurrElectrolyte_CBL_;
  icurrElectrolyte_CBR_ = r.icurrElectrolyte_CBR_;
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
porousLiKCl_dom1D::domain_prep(LocalNodeIndices *li_ptr)
{
  /*
   * First call the parent domain prep to get the node information
   */
  porousFlow_dom1D::domain_prep(li_ptr);

  /*
   * Figure out what the mass of the separator is
   * and then figure out its volume fraction to  
   * determine the cell porosity.
   *
   * We should read in the MgO.xml file to get the MgO density
   */

  int iph = (PSCinput_ptr->PhaseList_)->globalPhaseIndex(PSCinput_ptr->separatorPhase_);
  ThermoPhase* tmpPhase = & (PSCinput_ptr->PhaseList_)->thermo(iph);
  StoichSubstance* inert = dynamic_cast<Cantera::StoichSubstance *>( tmpPhase );

  double volumeSeparator = PSCinput_ptr->separatorArea_ * PSCinput_ptr->separatorThickness_;
  double volumeInert = PSCinput_ptr->separatorMass_ / inert->density() ;
  double porosity = 1.0 - volumeInert / volumeSeparator;

  std::cout << "Separator volume is " << volumeSeparator << " m^3 with "
            << volumeInert << " m^3 inert and porosity " << porosity <<  std::endl;

  for (int i = 0; i < NumLcCells; i++) {
    porosity_Cell_[i] = porosity;
    porosity_Cell_old_[i] = porosity;
  }

  /*
   * Porous electrode domain prep
   */
  concTot_Cell_.resize(NumLcCells, 0.0);
  concTot_Cell_old_.resize(NumLcCells, 0.0);

  Xleft_cc_.resize(3, 0.0);
  Xcent_cc_.resize(3, 0.0);
  Xright_cc_.resize(3, 0.0);
  spCharge_.resize(3, 0.0);
  spCharge_[0] = 1.0;
  spCharge_[1] = 1.0;
  spCharge_[2] = -1.0;

  mfElectrolyte_Soln_Curr_.resize(3, 0.0);
  mfElectrolyte_Thermo_Curr_.resize(3, 0.0);
  mfElectrolyte_Soln_Cell_old_.resize(3, NumLcCells, 0.0);
  mfElectrolyte_SolnDot_Curr_.resize(3, 0.0);
  pmVolElectrolyte_Curr_.resize(3, 0.0);

  gradX_trCurr_.resize(3, 0.0);
  Vdiff_trCurr_.resize(3, 0.0);
  jFlux_trCurr_.resize(3, 0.0);

  solnTemp.resize(10, 0.0);

  icurrElectrolyte_CBL_.resize(NumLcCells, 0.0);
  icurrElectrolyte_CBR_.resize(NumLcCells, 0.0);

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
porousLiKCl_dom1D::residEval(Epetra_Vector &res,
                             const bool doTimeDependentResid,
                             const Epetra_Vector *soln_ptr,
                             const Epetra_Vector *solnDot_ptr,
                             const Epetra_Vector *solnOld_ptr,
                             const double t,
                             const double rdelta_t,
                             ResidEval_Type_Enum residType,
			     const Solve_Type_Enum solveType)
{
  residType_Curr_ = residType;
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

  //  Electrolyte mass fluxes - this is rho V dot n at the boundaries of the cells
  double fluxFright = 0.;
  double fluxFleft;

  //mole fraction fluxes
  std::vector<double> fluxXright(nsp_, 0.0);
  std::vector<double> fluxXleft(nsp_, 0.0);

  double fluxL = 0.0;
  double fluxR = 0.0;

  const Epetra_Vector &soln = *soln_ptr;
  const doublereal * solnDot_Cellptr = 0;
  //const Epetra_Vector &solnDot = *solnDot_ptr;
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
   * offset of the electolyte solution unknowns at the current node
   */
  index_CentLcNode = Index_DiagLcNode_LCO[0];
  nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

  /*
   *  Offsets for the equation unknowns in the residual vector for the electrolyte domain
   */
  int EQ_Current_offset_BD = BDD_.EquationIndexStart_EqName[Current_Conservation];
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
  /*
   *  ------------------------------ LOOP OVER CELL -------------------------------------------------
   *  Loop over the number of Cells in this domain on this processor
   *  This loop is done from left to right.
   */
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    cIndex_cc_ = iCell;

#ifdef DEBUG_HKM_NOT
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
    } else {

      for (int k = 0; k < nsp_; k++) {
        Xright_cc_[k] = Xcent_cc_[k];
      }
      Vright_cc_ = Vcent_cc_;
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
        icurrElectrolyte_CBL_[iCell] = 0.0;
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
        fluxFleft = Fleft_cc_ * concTot_Curr_;

        /*
         * Calculate the flux of species and the flux of charge
         *   - the flux of charge must agree with the flux of species
         */
        icurrElectrolyte_CBL_[iCell] = 0.0;
        for (int k = 0; k < nsp_; k++) {
          fluxXleft[k] = jFlux_trCurr_[k];
          icurrElectrolyte_CBL_[iCell] += fluxXleft[k] * spCharge_[k];
          if (Fleft_cc_ > 0.0) {
            fluxXleft[k] += Fleft_cc_ * Xleft_cc_[k] * concTot_Curr_;
          } else {
            fluxXleft[k] += Fleft_cc_ * Xcent_cc_[k] * concTot_Curr_;
          }
        }
        icurrElectrolyte_CBL_[iCell] *= (Cantera::Faraday);
      }
    } else {
      /*
       * Copy the fluxes from the stored right side
       */
      fluxFleft = fluxFright;
      icurrElectrolyte_CBL_[iCell] = icurrElectrolyte_CBR_[iCell - 1];
      for (int k = 0; k < nsp_; k++) {
        fluxXleft[k] = fluxXright[k];
      }
    }
    /*
     * ------------------- CALCULATE FLUXES AT THE RIGHT BOUNDARY -------------------------------
     *
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
      // fluxFright = 0.0;
      if (solnDot_ptr) {
	solnDot_Cellptr =  &((*solnDot_ptr)[indexCent_EqnStart_BD]);
      } else {
	solnDot_Cellptr = 0;
      }
      SetupThermoShop1(&(soln[indexCent_EqnStart_BD]), solnDot_Cellptr, 0);
      Fright_cc_ = 0.0;
      fluxFright = Fright_cc_ * concTot_Curr_;
      icurrElectrolyte_CBR_[iCell] = 0.0;
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
      fluxFright = Fright_cc_ * concTot_Curr_;

      /*
       * Calculate the flux of species and the flux of charge
       *   - the flux of charge must agree with the flux of species
       */
      icurrElectrolyte_CBR_[iCell] = 0.0;
      for (int k = 0; k < nsp_; k++) {
        fluxXright[k] = jFlux_trCurr_[k];
        icurrElectrolyte_CBR_[iCell] += fluxXright[k] * spCharge_[k];
        if (Fright_cc_ > 0.0) {
          fluxXright[k] += Fright_cc_ * mfElectrolyte_Thermo_Curr_[k] * concTot_Curr_;
        } else {
          fluxXright[k] += Fright_cc_ * mfElectrolyte_Thermo_Curr_[k] * concTot_Curr_;
        }
      }
      icurrElectrolyte_CBR_[iCell] *= (Cantera::Faraday);
    }

#ifdef DEBUG_HKM_NOT
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
     * Add the flux terms into the residual
     */

    /*
     *  Total continuity equation - fluxFright and fluxFleft represent the total mass
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
     * Mole fraction summation equation
     */

    /*
     * Electroneutrality equation
     */

    /*
     *   Current conservation equation
     */
    res[indexCent_EqnStart_BD + EQ_Current_offset_BD] += (icurrElectrolyte_CBR_[iCell] - icurrElectrolyte_CBL_[iCell]);

    /*
     * --------------------------------------------------------------------------
     *  Add in the source terms at the current cell center
     * --------------------------------------------------------------------------
     */
    if (solnDot_ptr) {
      solnDot_Cellptr =  &((*solnDot_ptr)[indexCent_EqnStart_BD]);
    } else {
      solnDot_Cellptr = 0;
    }
    SetupThermoShop1(&(soln[indexCent_EqnStart_BD]), solnDot_Cellptr, 0);

#undef PRECIPITATE
#ifdef PRECIPITATE
    //print something if the salt should be precipitating
    checkPrecipitation();
#endif //PRECIPITATE
    /*
     *  Total continuity equation -
     */


    /*
     * Mole fraction summation equation
     *       For the DAE problem we specify the time derivative formulation. 
     *       For the regular problem, we specify the mole fraction formulation.
     *       Both are equivalent - > set m_isAlgebraic = 2
     */
    if (solveType == DAESystemInitial_Solve) { 
      res[indexCent_EqnStart_BD + EQ_MFSum_offset_BD] = 0.0;
      for (int k = 0; k < nsp_; k++) {
	res[indexCent_EqnStart_BD + EQ_MFSum_offset_BD] += mfElectrolyte_SolnDot_Curr_[k]; 
      }
    } else {
      res[indexCent_EqnStart_BD + EQ_MFSum_offset_BD] = 1.0;
      for (int k = 0; k < nsp_; k++) {
	res[indexCent_EqnStart_BD + EQ_MFSum_offset_BD] -= Xcent_cc_[k];
      }
    }

    /*
     * Electroneutrality equation
     *       For the DAE problem we specify the time derivative formulation. 
     *       For the regular problem, we specify the mole fraction formulation.
     *       Both are equivalent - > set m_isAlgebraic = 2
     */
    if (solveType == DAESystemInitial_Solve) {
      res[indexCent_EqnStart_BD + EQ_ChargeBal_offset_BD] = 0.0;
      for (int k = 0; k < nsp_; k++) {
	res[indexCent_EqnStart_BD + EQ_ChargeBal_offset_BD] += mfElectrolyte_SolnDot_Curr_[k] * spCharge_[k];
      }
    } else {
      res[indexCent_EqnStart_BD + EQ_ChargeBal_offset_BD] = 0.0;
      for (int k = 0; k < nsp_; k++) {
	res[indexCent_EqnStart_BD + EQ_ChargeBal_offset_BD] += Xcent_cc_[k] * spCharge_[k];
      }
    }

    /*
     *   Current conservation equation
     */

    /*
     * Special section if we own the left node of the domain. If we do
     * it will be cell 0. Store currents for later use.
     * These are the correct currents that work for the global balances
     */
    if (IOwnLeft && iCell == 0) {
      if (residType == Base_ShowSolution) {
        icurrElectrolyte_CBL_[iCell] =  icurrElectrolyte_CBR_[iCell];
      }
    }
    /*
     * Special section if we own the right node of the domain. If we do
     * it will be cell 0. Store currents for later use.
     * These are the correct currents that work for the global balances
     */
    if (IOwnRight && iCell == (NumLcCells - 1)) {
      if (residType == Base_ShowSolution) {
        icurrElectrolyte_CBR_[iCell] =  icurrElectrolyte_CBL_[iCell];
      }
    }
    /*
     *   ------------------ ADD IN THE TIME DEPENDENT TERMS ----------------------------------------------------
     */
    if (doTimeDependentResid) {

#ifdef DEBUG_HKM_NOT
      if (residType == Base_ResidEval) {
        printf(" Cell = %d, Totalflux_Li+_r = %10.3e,  = %10.3e, Totalflux_Li+_l ", iCell, fluxXright[0], fluxXleft[0]);
      }
#endif
      double newStuffTC = concTot_Curr_ * porosity_Curr_ * xdelCell;
      double newStuffSpecies0 = Xcent_cc_[0] * newStuffTC;

      /*
       * Setup shop with the old time step
       */
      SetupThermoShop1(&(solnOld[indexCent_EqnStart_BD]), 0,  0);

      double oldStuffTC = concTot_Curr_ * porosity_Curr_ * xdelCell;
      double oldStuffSpecies0 = mfElectrolyte_Soln_Curr_[0] * oldStuffTC;
      double tmp = (newStuffSpecies0 - oldStuffSpecies0) * rdelta_t;
      // double delta_concTot_CurrNew = concTot_CurrNew - concTot_Curr_;
      //  double delta_mf = mfNew - mfElectrolyte_Soln_Curr_[0];
#ifdef DEBUG_HKM_NOT
      if (residType == Base_ResidEval) {
        printf(" deltaT term = %10.3e BulkSum = %10.3e\n", tmp, tmp + (fluxXright[0] - fluxXleft[0]));
      }
#endif
      // HKM WARNING -> THIS IS DIFFERENT and will have to be looked at. I've down away with SetupThermoShop1 on the
      //                old solution in the other objects.
      /*
       *   .................... Add these terms in the residual
       */
      /*
       *  For this calculation we do the time derivative another way
       */
      if (solveType == DAESystemInitial_Solve) {
	
	double tmp2;
	/*
	 * Species continuity Equation - 2
	 */
	for (int k = 0; k < nsp_ - 2; k++) {
	  tmp2 = xdelCell * (concTot_Curr_    * porosity_Curr_ * mfElectrolyte_SolnDot_Curr_[k] + 
			     concTotDot_Curr_ * porosity_Curr_ * Xcent_cc_[k] );
	  res[indexCent_EqnStart_BD + EQ_Species_offset_BD + k] += tmp2;
	}
	
      } else {
	res[indexCent_EqnStart_BD + EQ_Species_offset_BD + 0] += tmp;
      }
      /*
       *   Add in the time term for the total continuity equation
       *         note: the current problem will have this term equally zero always.
       *               However, we put it in here for the next problem.
       */
      res[indexCent_EqnStart_BD + EQ_TCont_offset_BD] += (newStuffTC - oldStuffTC) * rdelta_t;
      /*
       *   .................... Go back to setting up shop at the current time
       */
      if (solnDot_ptr) {
	solnDot_Cellptr =  &((*solnDot_ptr)[indexCent_EqnStart_BD]);
      } else {
	solnDot_Cellptr = 0;
      }
      SetupThermoShop1(&(soln[indexCent_EqnStart_BD]), solnDot_Cellptr, 0);
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
void
porousLiKCl_dom1D::SetupThermoShop1(const doublereal * const solnElectrolyte_Curr, 
				    const doublereal * const solnDotElectrolyte_Curr, int type)
{
  updateElectrolyte(solnElectrolyte_Curr, solnDotElectrolyte_Curr);
  if (type == 0) {
    porosity_Curr_ = porosity_Cell_[cIndex_cc_];
  }
}
//=====================================================================================================================
void
porousLiKCl_dom1D::SetupThermoShop2(const doublereal * const solnElectrolyte_CurrL,
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
  updateElectrolyte(&solnTemp[0], 0);
}
  //=====================================================================================================================
  // Function updates the ThermoPhase object for the electrolyte given the solution vector and solution dot vector
  /*
   *   Routine will update the molten salt ThermoPhase object with the current state of the electrolyte.
   *   It also fills up the _Curr_ fields in this object refering to the Electrolyte properties.
   *
   *  @param solnElectrolyte      Vector of the solution at the current cell and bulk domain
   *  @param solnDotElectrolyte   Vector of the solution dot at the current cell and bulk domain
   */
  void
  porousLiKCl_dom1D::updateElectrolyte(const doublereal * const solnElectrolyte_Curr, 
				       const doublereal * const solnDotElectrolyte_Curr)
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
    /*
     * Fill up _Curr_ mf vectors for the electrolyte
     */
    getMFElectrolyte_soln(solnElectrolyte_Curr, solnDotElectrolyte_Curr);
    /*
     * Fill up _Curr_ Voltage values for the electrolyte
     */
    getVoltages(solnElectrolyte_Curr, 0);
    /*
     *  Set the ThermoPhase object with the current solution values
     */
    ionicLiquid_->setState_TPX(temp_Curr_, pres_Curr_, &mfElectrolyte_Thermo_Curr_[0]);
    ionicLiquid_->setElectricPotential(phiElectrolyte_Curr_);
    /*
     * Get the partial molar volumes
     */
    ionicLiquid_->getPartialMolarVolumes(DATA_PTR(pmVolElectrolyte_Curr_));

    // Calculate the total concentration of the electrolyte kmol m-3.
    concTot_Curr_ = ionicLiquid_->molarDensity();


    // Calculate the time derivative of the concentration of the electrolyte
    concTotDot_Curr_ = 0.0;
    for (int k = 0; k < 3; k++) {
      concTotDot_Curr_ += pmVolElectrolyte_Curr_[k] * mfElectrolyte_SolnDot_Curr_[k];
    }
    concTotDot_Curr_ *= (- concTot_Curr_ * concTot_Curr_);

  }
//=====================================================================================================================
void
porousLiKCl_dom1D::getVoltages(const double * const solnElectrolyte_Curr, const double * const solnSolid_Curr)
{
  int indexVS = BDD_.VariableIndexStart_VarName[Voltage];
  phiElectrolyte_Curr_ = solnElectrolyte_Curr[indexVS];
}
//=====================================================================================================================
void
porousLiKCl_dom1D::getMFElectrolyte_soln(const double * const solnElectrolyte_Curr, 
					 const double * const solnDotElectrolyte_Curr)
{
  int indexMF = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
  mfElectrolyte_Soln_Curr_[0] = solnElectrolyte_Curr[indexMF];
  mfElectrolyte_Soln_Curr_[1] = solnElectrolyte_Curr[indexMF + 1];
  mfElectrolyte_Soln_Curr_[2] = solnElectrolyte_Curr[indexMF + 2];
  double mf0 = std::max(mfElectrolyte_Soln_Curr_[0], 0.0);
  double mf1 = std::max(mfElectrolyte_Soln_Curr_[1], 0.0);
  double tmp = mf0 + mf1;

  mfElectrolyte_Thermo_Curr_[0] = (mf0) * 0.5 / tmp;
  mfElectrolyte_Thermo_Curr_[1] = (mf1) * 0.5 / tmp;
  mfElectrolyte_Thermo_Curr_[2] = 0.5;

  if (solnDotElectrolyte_Curr) {
    mfElectrolyte_SolnDot_Curr_[0] = solnDotElectrolyte_Curr[indexMF];
    mfElectrolyte_SolnDot_Curr_[1] = solnDotElectrolyte_Curr[indexMF + 1];
    mfElectrolyte_SolnDot_Curr_[2] = solnDotElectrolyte_Curr[indexMF + 2];
  }
}
//=====================================================================================================================
void
porousLiKCl_dom1D::SetupTranShop(const double xdel, const int type)
{

  /*
   * Determine diffusion velocities
   */

  //set gradients
  gradT_trCurr_ = 0.0;

  if (type == 0) {
    // Left boundary
    gradV_trCurr_ = (Vcent_cc_ - Vleft_cc_) / xdel;

  } else {
    // Right boundary
    gradV_trCurr_ = (Vright_cc_ - Vcent_cc_) / xdel;
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
    jFlux_trCurr_[k] = mfElectrolyte_Soln_Curr_[k] * concTot_Curr_ * Vdiff_trCurr_[k];
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
porousLiKCl_dom1D::showSolution(const Epetra_Vector *soln_GlAll_ptr,
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
  //bool doWrite = (NumOwnedNodes > 0) || duplicateOnAllProcs;
  std::vector<VarType> &variableNameList = BDD_.VariableNameList;
  int iBlock;
  int iGbNode;
  int n;
  int nPoints = BDD_.LastGbNode - BDD_.FirstGbNode + 1;

  char buf[100];
  string indent = "";
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
    drawline(indentSpaces, 80);
    ss.print0("%s  Solution on Bulk Domain %12s : Number of variables = %d\n", ind, sss.c_str(), NumDomainEqns);
    ss.print0("%s                                         : Number of Nodes = %d\n", ind, nPoints);
    ss.print0("%s                                         : Beginning pos %g\n", ind, BDD_.Xpos_start);
    ss.print0("%s                                         : Ending    pos %g\n", ind, BDD_.Xpos_end);
  }
  if (do0Write) {
    for (iBlock = 0; iBlock < nn; iBlock++) {
      drawline(indentSpaces, 80);
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
      Cantera::writelog(buf);
      for (n = 0; n < nrem; n++) {
        int ivar = nn * 5 + n;
        VarType vt = variableNameList[ivar];
        string name = vt.VariableName(15);
        ss.print0(" %15s", name.c_str());
        Cantera::writelog(buf);
      }
      ss.print0("\n");
      drawline(indentSpaces, 80);

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
    // ----------------------------------------------------
    // --             PRINT FLUXES AT THE CELL BOUNDARIES --
    // ----------------------------------------------------
    ss.print0("%s    CellBound    z     IcurrElectrolyte ", ind);
    ss.print0("\n");
    drawline(indentSpaces, 80);
  }
  NodalVars *nvl;
  NodalVars *nvr;
  doublereal x;
  int iCell;
  for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
      ss.print0("%s    ", ind);
      if (iGbNode == BDD_.FirstGbNode) {

        iCell = 0;
        nvr = gi->NodalVars_GbNode[BDD_.FirstGbNode];
        x = nvr->xNodePos();
        ss.print0("Lft-0     %11.4E ", x);
        ss.print0("%11.4E ", icurrElectrolyte_CBL_[iCell]);
      } else if (iGbNode == BDD_.LastGbNode) {
        iCell = BDD_.LastGbNode - BDD_.FirstGbNode;
        nvr = gi->NodalVars_GbNode[BDD_.LastGbNode];
        x = nvr->xNodePos();
        ss.print0("%3d-Rgt   %11.4E ", iCell, x);
        ss.print0("%11.4E ", icurrElectrolyte_CBR_[iCell]);
      } else {
        iCell = iGbNode - BDD_.FirstGbNode;
        nvl = gi->NodalVars_GbNode[iGbNode];
        nvr = gi->NodalVars_GbNode[iGbNode + 1];
        x = 0.5 * (nvl->xNodePos() + nvr->xNodePos());
        ss.print0("%3d-%-3d   %11.4E ", iCell, iCell + 1, x);
        ss.print0("%11.4E ", icurrElectrolyte_CBR_[iCell]);
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
porousLiKCl_dom1D::initialConditions(const bool doTimeDependentResid,
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

    soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = 0.0;

    /*
     * Get initial mole fractions from PSinput
     */
    int iLiCl = PSCinput_ptr->PhaseList_->globalSpeciesIndex ("LiCl(L)");
    int iKCl =  PSCinput_ptr->PhaseList_->globalSpeciesIndex ("KCl(L)");
    int iLi =   PSCinput_ptr->PhaseList_->globalSpeciesIndex ("Li+");
    int iK =    PSCinput_ptr->PhaseList_->globalSpeciesIndex ("K+");

    soln[indexCent_EqnStart_BD + iVar_Species_BD + 0] = 
      PSCinput_ptr->electrolyteMoleFracs_[iLiCl] / 2.0 + PSCinput_ptr->electrolyteMoleFracs_[iLi];
    soln[indexCent_EqnStart_BD + iVar_Species_BD + 1] = 
      PSCinput_ptr->electrolyteMoleFracs_[iKCl] / 2.0 + PSCinput_ptr->electrolyteMoleFracs_[iK];
    soln[indexCent_EqnStart_BD + iVar_Species_BD + 2] = 0.5;

    soln[indexCent_EqnStart_BD + iVar_Voltage_BD] = -0.075;  
  
    // Let's compute a characteristic electrical conductivity based on these initial conditions.
    //We only need to do this for one control volume    
    if (iCell == 0 ) {
      //Fill thermodynamic state
      SetupThermoShop1(&(soln[indexCent_EqnStart_BD]), 0, 0);
      //This computes conductivity based on unit potential gradient and no species gradients
      electrolyteConduct_ = trans_->getElectricConduct();
      fprintf(stderr, "The separator conductivity is %g S/m before porosity\n",electrolyteConduct_ );
      //  Correct species diffusivities according to the tortuosity
      // Porosity is known from porosity_Cell_ set in domain_prep
      double bruggemannExp = 2.5;
      Tortuosity tort(bruggemannExp );
      electrolyteConduct_ *= tort.tortuosityFactor( porosity_Cell_[iCell] );
      fprintf(stderr, "and %g S/m after porosity\n",electrolyteConduct_ );
    }  //end if (iCell==0)
    
  }  
}
//=====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 * @param atolVector Reference for the atol vector to fill up
 */
void porousLiKCl_dom1D::setAtolVector(double atolDefault, const Epetra_Vector_Ghosted & soln, 
				      Epetra_Vector_Ghosted & atolVector,
				      const Epetra_Vector_Ghosted * const atolV)
{
 
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
    int iVAR_Vaxial_BD  = BDD_.VariableIndexStart_VarName[Velocity_Axial];
    int iVar_Species_BD = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
    int iVar_Voltage_BD = BDD_.VariableIndexStart_VarName[Voltage];

    /*
     * Set the atol value for the axial velocity
     *   arithmetically scaled.
     */
    double vax = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];
    atolVector[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = std::max(1.0E-4, 1.0E-1 * vax);

    /*
     * Set atol values for the species mole fractions
     */
    double val = atolDefault * 1.0E-2;
    if (val < 1.0E-12) {
      val = 1.0E-12;
    }
    atolVector[indexCent_EqnStart_BD + iVar_Species_BD + 0] = val;
    atolVector[indexCent_EqnStart_BD + iVar_Species_BD + 1] = val;
    atolVector[indexCent_EqnStart_BD + iVar_Species_BD + 2] = val;

    /*
     * Set the atol value for the electrolyte voltage
     *      arithmetically scaled.
     */
    atolVector[indexCent_EqnStart_BD + iVar_Voltage_BD] = 0.05;
  }  
}
//=====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 * @param atolVector Reference for the atol vector to fill up
 */
void porousLiKCl_dom1D::setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted & soln, 
					      const Epetra_Vector_Ghosted & solnDot,
					      Epetra_Vector_Ghosted & atolVector,
					      const Epetra_Vector_Ghosted * const atolV)
{
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

    /*
     * Set the atol value for the axial velocity
     *   arithmetically scaled.
     */
    double vax = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];
    atolVector[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = std::max(1.0E-4, 1.0E-1 * vax);

    /*
     * Set atol values for the species mole fractions
     */
    double val = atolDefault * 1.0E-2;
    if (val < 1.0E-6) {
      val = 1.0E-6;
    }
    atolVector[indexCent_EqnStart_BD + iVar_Species_BD + 0] = val;
    atolVector[indexCent_EqnStart_BD + iVar_Species_BD + 1] = val;
    atolVector[indexCent_EqnStart_BD + iVar_Species_BD + 2] = val;

    /*
     * Set the atol value for the electrolyte voltage
     *      arithmetically scaled.
     */
    atolVector[indexCent_EqnStart_BD + iVar_Voltage_BD] = 0.05;
  }  
}
//=====================================================================================================================
void
porousLiKCl_dom1D::setAtolDeltaDamping(double atolDefault, double relcoeff, 
				       const Epetra_Vector_Ghosted & soln, 
				       Epetra_Vector_Ghosted & atolDeltaDamping,
				       const Epetra_Vector_Ghosted * const atolV)
{

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

    /*
     * Set the atol value for the axial velocity
     *   arithmetically scaled.
     */
    if (iVAR_Vaxial_BD >= 0) {
      //   double vax = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];
      atolDeltaDamping[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = 1.0E-4 * relcoeff;
    }

    /*
     * Set atol values for the species mole fractions
     */
    double val = atolDefault * 1.0E3;
    if (val < 1.0E-4) {
      val = 1.0E-4;
    }
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Species_BD + 0] = val * relcoeff;
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Species_BD + 1] = val * relcoeff;
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Species_BD + 2] = val * relcoeff;
    
    /*
     * Set the atol value for the electrolyte voltage
     *      arithmetically scaled.
     */
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Voltage_BD] = 0.05 * relcoeff;
  }
}
//=====================================================================================================================
void
porousLiKCl_dom1D::setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff, 
					       const Epetra_Vector_Ghosted & soln, 
					       const Epetra_Vector_Ghosted & solnDot, 
					       Epetra_Vector_Ghosted & atolDeltaDamping,
					       const Epetra_Vector_Ghosted * const atolV)
{

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

    /*
     * Set the atol value for the axial velocity
     *   arithmetically scaled.
     */
    if (iVAR_Vaxial_BD >= 0) {
      //   double vax = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];
      atolDeltaDamping[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = 1.0E-2 * relcoeff;
    }

    /*
     * Set atol values for the species mole fractions
     */
    double val = atolDefault * 1.0E7;
    if (val < 1.0E-1) {
      val = 1.0E-1;
    }
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Species_BD + 0] = val * relcoeff;
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Species_BD + 1] = val * relcoeff;
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Species_BD + 2] = val * relcoeff;
    
    /*
     * Set the atol value for the electrolyte voltage
     *      arithmetically scaled.
     */
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Voltage_BD] = 0.05 * relcoeff;
  }
}

//=====================================================================================================================
void
porousLiKCl_dom1D::err(const char *msg)
{
  printf("porousLiKCl_dom1D: function not implemented: %s\n", msg);
  exit(-1);
}
//=====================================================================================================================
/**
 * Return an estimate for the electrolyte conductivity
 */
int 
porousLiKCl_dom1D::getConductivity(  ) {
  return electrolyteConduct_;
}
//=====================================================================================================================
/**
 * Method to check for precipitation of the salts.  
 * Returns index of offending cation or -1 if no precipitation
 */
int 
porousLiKCl_dom1D::checkPrecipitation(  ) {

  //mole fraction of the electrolyte ions are held in
  // mfElectrolyte_Thermo_Curr_[3]

  // molten salt phase
  string id_salt = "LiKCl_Margules";
  int iph = (PSCinput_ptr->PhaseList_)->globalPhaseIndex(id_salt);
  if (iph < 0) {
    throw CanteraError("porousLiKCl_LiSiAnode_dom1D::checkPrecipitation()",
                       "Can't find the phase in the phase list: " + id_salt);
  }
  ThermoPhase* tmpPhase = & (PSCinput_ptr->PhaseList_)->thermo(iph);
  
  MargulesVPSSTP *salt ;
  salt = dynamic_cast<Cantera::MargulesVPSSTP *>( tmpPhase->duplMyselfAsThermoPhase() );
  
  int iKCl_l = salt->speciesIndex("KCl(L)");
  int iLiCl_l = salt->speciesIndex("LiCl(L)");
  int iK_ion = ionicLiquid_->speciesIndex("K+");
  int iLi_ion = ionicLiquid_->speciesIndex("Li+");

  //solid LiCl phase
  id_salt = "LiCl(S)";
  iph = (PSCinput_ptr->PhaseList_)->globalPhaseIndex(id_salt);
  if (iph < 0) {
    throw CanteraError("porousLiKCl_LiSiAnode_dom1D::checkPrecipitation()",
                       "Can't find the phase in the phase list: " + id_salt);
  }
  tmpPhase = & (PSCinput_ptr->PhaseList_)->thermo(iph);
  Cantera::ThermoPhase *LiCl_solid = tmpPhase->duplMyselfAsThermoPhase() ;

  //solid KCl phase
  id_salt = "KCl(S)";
  iph = (PSCinput_ptr->PhaseList_)->globalPhaseIndex(id_salt);
  if (iph < 0) {
    throw CanteraError("porousLiKCl_LiSiAnode_dom1D::checkPrecipitation()",
                       "Can't find the phase in the phase list: " + id_salt);
  }
  tmpPhase = & (PSCinput_ptr->PhaseList_)->thermo(iph);
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
}
//=====================================================================================================================
