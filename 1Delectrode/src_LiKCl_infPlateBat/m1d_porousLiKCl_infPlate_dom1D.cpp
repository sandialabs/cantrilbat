/*
 * m1d_porousLiKCl_infPlate_dom1D.cpp
 *
 *  $Id: m1d_porousLiKCl_dom1D.cpp 540 2013-02-27 22:18:26Z hkmoffa $ 
 */

//  This is a heavyweight base class that provides the function
//evaluation for a single bulk domain.
#include "m1d_porousLiKCl_infPlate_dom1D.h"
#include "m1d_BDT_porousLiKCl.h"

#include "m1d_NodalVars.h"
#include "m1d_LocalNodeIndices.h"

#include "m1d_GlobalIndices.h"
#include "m1d_exception.h"
#include "m1d_DomainLayout.h"

#include "cantera/base/ctml.h"

#include "stdio.h"
#include "stdlib.h"

#include "m1d_ProblemStatementCell.h"
#include "m1d_porousFlow_dom1D.h"

extern m1d::ProblemStatementCell PSinput;

using namespace std;
using namespace Cantera;

namespace m1d
{

//==============================================================================
porousLiKCl_infPlate_dom1D::porousLiKCl_infPlate_dom1D(BDT_porousLiKCl& bdd) :
  porousFlow_dom1D(bdd), 
  temp_Curr_(TemperatureReference_), 
  ivb_(VB_MOLEAVG)
{

  BDT_porousLiKCl *fa = dynamic_cast<BDT_porousLiKCl *> (&bdd);
  if (!fa) {
    throw m1d_Error("confused", "confused");
  }
  ionicLiquid_ = fa->ionicLiquidIFN_;
  trans_ = fa->trans_;
  nsp_ = 3;
  nph_ = 1;

}
//==============================================================================
porousLiKCl_infPlate_dom1D::porousLiKCl_infPlate_dom1D(const porousLiKCl_infPlate_dom1D &r) :
  porousFlow_dom1D((BDT_porousLiKCl&) r.BDD_), 
  temp_Curr_(TemperatureReference_), 
  ivb_(VB_MOLEAVG)
{
  porousLiKCl_infPlate_dom1D::operator=(r);
}
//==============================================================================
porousLiKCl_infPlate_dom1D::~porousLiKCl_infPlate_dom1D()
{
}
//==============================================================================
porousLiKCl_infPlate_dom1D &
porousLiKCl_infPlate_dom1D::operator=(const porousLiKCl_infPlate_dom1D &r)
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

  porosity_Cell_ = r.porosity_Cell_;
  temp_Curr_ = r.temp_Curr_;
  ivb_ = r.ivb_;

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
porousLiKCl_infPlate_dom1D::domain_prep(LocalNodeIndices *li_ptr)
{
  /*
   * First call the parent domain prep to get the node information
   */
  BulkDomain1D::domain_prep(li_ptr);

  /*
   * Figure out what the mass of the separator is
   * and then figure out its volume fraction to  
   * determine the cell porosity.
   *
   * We should read in the MgO.xml file to get the MgO density
   */
  //  StoichSubstanceDef inert( PSinput.separatorXMLFile_, 
  //			    PSinput.separatorPhase_ );
  int iph = (PSinput.PhaseList_)->globalPhaseIndex(PSinput.separatorPhase_);
  ThermoPhase *tmpPhase = &(PSinput.PhaseList_)->thermo(iph);
 // ThermoPhase* tmpPhase = (PSinput.PhaseList_)->getPhase( PSinput.separatorPhase_ );
  StoichSubstance* inert = dynamic_cast<Cantera::StoichSubstance *>(tmpPhase);

  //need to convert inputs from cgs to SI
  double volumeSeparator = 
    PSinput.separatorArea_ * PSinput.separatorThickness_; 
  double volumeInert = 
    PSinput.separatorMass_ / inert->density() ;
  double porosity = 1.0 - volumeInert / volumeSeparator;

  std::cout << "Separator volume is " 
	    << volumeSeparator << " m^3 with "
	    << volumeInert << " m^3 inert and porosity "
	    << porosity 
	    <<  std::endl;

  /*
   * Porous electrode domain prep
   */
  porosity_Cell_.resize(NumLcCells, porosity);
  porosity_Cell_old_.resize(NumLcCells, porosity);
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

  gradX_trCurr_.resize(3, 0.0);
  Vdiff_trCurr_.resize(3, 0.0);
  jFlux_trCurr_.resize(3, 0.0);

  solnTemp.resize(10, 0.0);

  for (int i = 0; i < NumLcCells; i++) {
    porosity_Cell_[i] = porosity;
  }

  /*
   *  Set the velocity basis of the transport object. We are using
   *  mole-averaged velocities as the basis.
   */
  trans_->setVelocityBasis(ivb_);
}
//==============================================================================
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
porousLiKCl_infPlate_dom1D::residEval(Epetra_Vector &res,
                             const bool doTimeDependentResid,
                             const Epetra_Vector *soln_ptr,
                             const Epetra_Vector *solnDot_ptr,
                             const Epetra_Vector *solnOld_ptr,
                             const double t,
                             const double rdelta_t,
                             const ResidEval_Type_Enum residType,
			     const Solve_Type_Enum solveType )
{
  residType_Curr_ = residType;
  const double surfArea = 1.0;

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

  //  Electrostatic fluxes (current)
  double fluxVright = 0.;
  double fluxVleft;
  //mole fraction fluxes
  std::vector<double> fluxXright(nsp_, 0.0);
  std::vector<double> fluxXleft(nsp_, 0.0);

  double fluxL = 0.0;
  double fluxR = 0.0;

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
        fluxVleft = 0.0;
        for (int k = 0; k < nsp_; k++) {
          fluxXleft[k] = 0.0;
        }
      } else {

        /*
         *  Establish the environment at the left cell boundary
         */
	  SetupThermoShop2(nodeLeft, &(soln[indexLeft_EqnStart_BD]), nodeCent, &(soln[indexCent_EqnStart_BD]), 0);

        SetupTranShop(xdelL, 0);

        /*
         * Calculate the flux at the left boundary for each equation
         */
        fluxFleft = Fleft_cc_ * concTot_Curr_ * porosity_Curr_;

        /*
         * Calculate the flux of species and the flux of charge
         *   - the flux of charge must agree with the flux of species
         */
        fluxVleft = 0.0;
        for (int k = 0; k < nsp_; k++) {
          fluxXleft[k] = jFlux_trCurr_[k] * porosity_Curr_;
          fluxVleft += fluxXleft[k] * spCharge_[k];
          if (Fleft_cc_ > 0.0) {
            fluxXleft[k] += Fleft_cc_ * Xleft_cc_[k] * concTot_Curr_ * porosity_Curr_;
          } else {
            fluxXleft[k] += Fleft_cc_ * Xcent_cc_[k] * concTot_Curr_ * porosity_Curr_;
          }
        }
        fluxVleft *= (Cantera::Faraday);
      }
    } else {
      /*
       * Copy the fluxes from the stored right side
       */
      fluxFleft = fluxFright;
      fluxVleft = fluxVright;
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
      SetupThermoShop1Old(&(soln[indexCent_EqnStart_BD]), 0);
      fluxFright = Fright_cc_ * concTot_Curr_ * porosity_Curr_;
      fluxVright = 0.0;
      for (int k = 0; k < nsp_; k++) {
        fluxXright[k] = 0.0;
      }
    } else {
      /*
       *  Establish the environment at the right cell boundary
       */
	SetupThermoShop2(nodeCent, &(soln[indexCent_EqnStart_BD]), nodeRight, &(soln[indexRight_EqnStart_BD]), 1);

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
      fluxVright = 0.0;
      for (int k = 0; k < nsp_; k++) {
        fluxXright[k] = jFlux_trCurr_[k];
        fluxVright += fluxXright[k] * spCharge_[k];
        if (Fright_cc_ > 0.0) {
          //fluxXright[k] += Fright_cc_ * Xcent_cc_[k] * concTot_Curr_ * porosity_Curr_;
	  fluxXright[k] += Fright_cc_ * mfElectrolyte_Thermo_Curr_[k] * concTot_Curr_ * porosity_Curr_;
        } else {
          fluxXright[k] += Fright_cc_ * mfElectrolyte_Thermo_Curr_[k]  * concTot_Curr_ * porosity_Curr_;
        }
      }
      fluxVright *= (Cantera::Faraday);
    }



#ifdef DEBUG_HKM
    /*
    if (doTimeDependentResid) {
    
      if (residType == Base_ResidEval) {
        printf(" Cell = %d, Totalflux_K+ = %10.3e,  Totalflux_Cl- = %10.3e \n", iCell,  fluxXright[1], fluxXright[2]);
	    printf("           Vmolal = %10.3e, jd_Li+ = %10.3e  jd_K+ = %10.3e jd_Cl- = %10.3e\n",
	       Fright_cc_ ,  jFlux_trCurr_[0], jFlux_trCurr_[1], jFlux_trCurr_[2]);
	    printf("           Vmolal = %10.3e, vd_Li+ = %10.3e  vd_K+ = %10.3e vd_Cl- = %10.3e\n",
	       Fright_cc_ , Vdiff_trCurr_[0], Vdiff_trCurr_[1],  Vdiff_trCurr_[2]);
      }
    }
    */
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
    res[indexCent_EqnStart_BD + EQ_Current_offset_BD] += (fluxVright - fluxVleft);

    /*
     * --------------------------------------------------------------------------
     *  Add in the source terms at the current cell center
     * --------------------------------------------------------------------------
     */

    SetupThermoShop1Old(&(soln[indexCent_EqnStart_BD]), 0);
    /*
     *  Total continuity equation -
     */

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
     *   Current conservation equation
     */

    /*
     * Supplementary case where we are adding on time dependent terms to the residual
     */
    if (doTimeDependentResid) {

#ifdef DEBUG_HKM
      //double concTot_CurrNew = concTot_Curr_;
      //double mfNew = mfElectrolyte_Soln_Curr_[0];
      if (residType == Base_ResidEval) {
        //printf(" Cell = %d, Totalflux_Li+_r = %10.3e,  = %10.3e, Totalflux_Li+_l ", iCell, fluxXright[0], fluxXleft[0]);
      }
#endif
      double newStuff = Xcent_cc_[0] * concTot_Curr_ * porosity_Curr_ * xdelCell;

      /*
       * Setup shop with the old time step
       */
      SetupThermoShop1Old(&(solnOld[indexCent_EqnStart_BD]), 0);

      double oldStuff = mfElectrolyte_Soln_Curr_[0] * concTot_Curr_ * porosity_Curr_ * xdelCell;
      double tmp = surfArea * (newStuff - oldStuff) * rdelta_t;
     // double delta_concTot_CurrNew = concTot_CurrNew - concTot_Curr_;
     // double delta_mf = mfNew - mfElectrolyte_Soln_Curr_[0];
#ifdef DEBUG_HKM
      if (residType == Base_ResidEval) {
       // printf(" deltaT term = %10.3e BulkSum = %10.3e\n", tmp, tmp + (fluxXright[0] - fluxXleft[0]));
      }
#endif
      res[indexCent_EqnStart_BD + EQ_Species_offset_BD + 0] += tmp;

      SetupThermoShop1Old(&(soln[indexCent_EqnStart_BD]), 0);
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
porousLiKCl_infPlate_dom1D::SetupThermoShop1Old(const doublereal * const solnElectrolyte_Curr, int type)
{
    //updateElectrolyteOld(solnElectrolyte_Curr);
    if (type == 0) {
	porosity_Curr_ = porosity_Cell_[cIndex_cc_];
    }
}
//=====================================================================================================================
//  Setup shop at a particular point in the domain, calculating intermediate quantites
//  and updating Cantera's objects
/*
 *  All member data with the suffix, _Curr_, are updated by this function.
 *
 * @param solnElectrolyte_Curr  Current value of the solution vector
 * @param type                  Type of call
 *                              0 - at the current cell center
 */
void
porousLiKCl_infPlate_dom1D::SetupThermoShop1(const NodalVars* const nv, const doublereal* const solnElectrolyte_Curr)
{
    porosity_Curr_ = porosity_Cell_[cIndex_cc_];
    updateElectrolyte(nv, solnElectrolyte_Curr);
}

//===================================================================================================================
void
porousLiKCl_infPlate_dom1D::SetupThermoShop2(const NodalVars* const nvL, const doublereal* const solnElectrolyte_CurrL,
                                              const NodalVars* const nvR, const doublereal* const solnElectrolyte_CurrR,
                                              int type)
{

// Needs major work
   // for (int i = 0; i < BDD_.NumEquationsPerNode; i++) {
   //     solnTemp[i] = 0.5 * (solnElectrolyte_CurrL[i] + solnElectrolyte_CurrR[i]);
   // }

    double tempL = getPointTemperature(nvL, solnElectrolyte_CurrL);
    double tempR = getPointTemperature(nvR, solnElectrolyte_CurrR);
    temp_Curr_ = 0.5 * (tempL + tempR);
    /*
     * Get the pressure
     */
    pres_Curr_ = PressureReference_;

    size_t indexMFL = nvL->indexBulkDomainVar0(MoleFraction_Species);
    size_t indexMFR = nvR->indexBulkDomainVar0(MoleFraction_Species);


    mfElectrolyte_Soln_Curr_[0] = 0.5 * (solnElectrolyte_CurrL[indexMFL] +solnElectrolyte_CurrR[indexMFR]);
    mfElectrolyte_Soln_Curr_[1] = 0.5 * (solnElectrolyte_CurrL[indexMFL+1] +solnElectrolyte_CurrR[indexMFR+1]);
    mfElectrolyte_Soln_Curr_[2] = 0.5 * (solnElectrolyte_CurrL[indexMFL+2] +solnElectrolyte_CurrR[indexMFR+2]);
    double mf0 = std::max(mfElectrolyte_Soln_Curr_[0], 0.0);
    double mf1b = std::max(mfElectrolyte_Soln_Curr_[1], 0.0);
    double mf2b = std::max(mfElectrolyte_Soln_Curr_[2], 0.0);
    double mf1 = mf1b;
    double mf2 = mf2b;
    if (mf1b != mf2b) {
        mf1 = 0.5 * (mf1b + mf2b);
        mf2 = 0.5 * (mf1b + mf2b);
    }
    double tmp = mf0 + mf1 + mf2;

    mfElectrolyte_Thermo_Curr_[0] = mf0 / tmp;
    mfElectrolyte_Thermo_Curr_[1] = mf1 / tmp;
    mfElectrolyte_Thermo_Curr_[2] = mf2 / tmp;
 size_t indexVS = nvL->indexBulkDomainVar0(Voltage);
    double phiElectrolyteL = solnElectrolyte_CurrL[indexVS];
    indexVS = nvR->indexBulkDomainVar0(Voltage);
    double phiElectrolyteR = solnElectrolyte_CurrR[indexVS];
    phiElectrolyte_Curr_ = 0.5 * (phiElectrolyteL + phiElectrolyteR);

    if (type == 0) {
        porosity_Curr_ = 0.5 * (porosity_Cell_[cIndex_cc_ - 1] + porosity_Cell_[cIndex_cc_]);
    } else {
        porosity_Curr_ = 0.5 * (porosity_Cell_[cIndex_cc_ + 1] + porosity_Cell_[cIndex_cc_]);
    }
    /*
     *  Set the ThermoPhase states
     */
    ionicLiquid_->setState_TPX(temp_Curr_, pres_Curr_, &mfElectrolyte_Thermo_Curr_[0]);
    ionicLiquid_->setElectricPotential(phiElectrolyte_Curr_);
    //
    // Calculate the total concentration of the electrolyte kmol m-3 and store into concTot_Curr_
    //
    concTot_Curr_ = ionicLiquid_->molarDensity();

    //
    // Calculate the matrix thermal conductivity from a series resistance on the two sides
    //
    if (type == 0) {
        tmp = 1.0 / thermalCond_Cell_[cIndex_cc_ - 1] + 1.0 / thermalCond_Cell_[cIndex_cc_];
    } else {
        tmp = 1.0 / thermalCond_Cell_[cIndex_cc_ + 1] + 1.0 / thermalCond_Cell_[cIndex_cc_];
    }
    thermalCond_Curr_ = 2.0 / tmp;


}
//=====================================================================================================================
// Function updates the ThermoPhase object for the electrolyte
// given the solution vector
/*
 *   Routine will update the molten salt ThermoPhase object with the current state of the electrolyte
 *
 * @param solnElectrolyte
 */
/*
void
porousLiKCl_infPlate_dom1D::updateElectrolyteOld(const doublereal * const solnElectrolyte_Curr)
{
 
  temp_Curr_ = TemperatureReference_;
  int iTemp = BDD_.VariableIndexStart_VarName[Temperature];
  if (iTemp >= 0) {
    temp_Curr_ = solnElectrolyte_Curr[iTemp];
  }

  pres_Curr_ = PressureReference_;

  //getMFElectrolyte_solnOld(solnElectrolyte_Curr);
  //getVoltagesOld(solnElectrolyte_Curr, 0);

  ionicLiquid_->setState_TPX(temp_Curr_, pres_Curr_, &mfElectrolyte_Thermo_Curr_[0]);

  ionicLiquid_->setElectricPotential(phiElectrolyte_Curr_);

  // Calculate the total concentration of the electrolyte kmol m-3.
  concTot_Curr_ = ionicLiquid_->molarDensity();

}
*/
//=====================================================================================================================
// Function updates the ThermoPhase object for the electrolyte given the solution vector
void
porousLiKCl_infPlate_dom1D::updateElectrolyte(const NodalVars* const nv, const doublereal* const solnElectrolyte_Curr)
{
    /*
     * Get the temperature: Check to see if the temperature is in the solution vector.
     *   If it is not, then use the reference temperature
     */
    temp_Curr_ = getPointTemperature(nv, solnElectrolyte_Curr);
    /*
     * Get the pressure
     */
    pres_Curr_ = PressureReference_;
    /*
     *  Assemble electrolyte mole fractions into mfElectrolyte_Thermo_Curr_[]
     */
    getMFElectrolyte_soln(nv, solnElectrolyte_Curr);
    /*
     *  assemble electrolyte potential into phiElectrolyte_Curr_
     */
    getVoltages(nv, solnElectrolyte_Curr);
    /*
     *  Set the ThermoPhase states
     */
    ionicLiquid_->setState_TPX(temp_Curr_, pres_Curr_, &mfElectrolyte_Thermo_Curr_[0]);
    ionicLiquid_->setElectricPotential(phiElectrolyte_Curr_);
    //
    // Calculate the total concentration of the electrolyte kmol m-3 and store into concTot_Curr_
    //
    concTot_Curr_ = ionicLiquid_->molarDensity();
}

//=====================================================================================================================
/*
void
porousLiKCl_infPlate_dom1D::getVoltagesOld(const double * const solnElectrolyte_Curr, const double * const solnSolid_Curr)
{
  int indexVS = BDD_.VariableIndexStart_VarName[Voltage];
  phiElectrolyte_Curr_ = solnElectrolyte_Curr[indexVS];

  //int indexVE = SDD_.VariableIndexStart_VarName[Voltage];
  //phiAnode_ = solnSolid[indexVE];
}
*/
//=====================================================================================================================
// Retrieves the voltages from the solution vector and puts them into local storage
/*
 * @param solnElectrolyte start of the solution vector at the current node
 */
void
porousLiKCl_infPlate_dom1D::getVoltages(const NodalVars* const nv, const double* const solnElectrolyte)
{  
    size_t indexVS = nv->indexBulkDomainVar0(Voltage);
    phiElectrolyte_Curr_ = solnElectrolyte[indexVS];
    // phiElectrode_Curr_ = solnElectrolyte[indexVS + 1];
}
//=====================================================================================================================
/*
void
porousLiKCl_infPlate_dom1D::getMFElectrolyte_solnOld(const double * const solnElectrolyte_Curr)
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
}
*/
//=====================================================================================================================
void
porousLiKCl_infPlate_dom1D::getMFElectrolyte_soln(const NodalVars* const nv, const double* const solnElectrolyte_Curr)
{
    // We are assuming that mf species start with subindex 0 here
    size_t indexMF = nv->indexBulkDomainVar0(MoleFraction_Species);
    mfElectrolyte_Soln_Curr_[0] = solnElectrolyte_Curr[indexMF];
    mfElectrolyte_Soln_Curr_[1] = solnElectrolyte_Curr[indexMF + 1];
    mfElectrolyte_Soln_Curr_[2] = solnElectrolyte_Curr[indexMF + 2];
    double mf0 = std::max(mfElectrolyte_Soln_Curr_[0], 0.0);
    double mf1 = std::max(mfElectrolyte_Soln_Curr_[1], 0.0);
    double tmp = mf0 + mf1;
    
    mfElectrolyte_Thermo_Curr_[0] = (mf0) * 0.5 / tmp;
    mfElectrolyte_Thermo_Curr_[1] = (mf1) * 0.5 / tmp;
    mfElectrolyte_Thermo_Curr_[2] = 0.5;
}
//=====================================================================================================================
void
porousLiKCl_infPlate_dom1D::SetupTranShop(const double xdel, const int type)
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

  for (int k = 0; k < nsp_; k++) {
    jFlux_trCurr_[k] = mfElectrolyte_Soln_Curr_[k] * concTot_Curr_ * pow( porosity_Curr_, 0.5 ) * Vdiff_trCurr_[k];
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
porousLiKCl_infPlate_dom1D::initialConditions(const bool doTimeDependentResid, Epetra_Vector *soln_ptr,
                                     Epetra_Vector *solnDot, const double t,
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

    //Eutectic composition
    double xLiCl_eutectic = 0.585;
    soln[indexCent_EqnStart_BD + iVar_Species_BD + 0] = 0.5*(xLiCl_eutectic);
    soln[indexCent_EqnStart_BD + iVar_Species_BD + 1] = 0.5*(1-xLiCl_eutectic);
    soln[indexCent_EqnStart_BD + iVar_Species_BD + 2] = 0.5;

    soln[indexCent_EqnStart_BD + iVar_Voltage_BD] = -0.07;

  }
}
//=====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void porousLiKCl_infPlate_dom1D::setAtolVector(double atolDefault, const Epetra_Vector_Ghosted & soln, 
				      Epetra_Vector_Ghosted & atolVector, const Epetra_Vector_Ghosted * const atolV)
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

    /*
     * Set the atol value for the axial velocity
     *   arithmetically scaled -> so this is a characteristic value
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
    for (int k = 0; k < nsp_; k++) {
      atolVector[indexCent_EqnStart_BD + iVar_Species_BD + k] = val;
    }

    /*
     * Set the atol value for the electrolyte voltage
     *      arithmetically scaled.-> so this is a characteristic value
     *         1 kcal gmol-1 = 0.05 volts
     */
    atolVector[indexCent_EqnStart_BD + iVar_Voltage_BD] = 0.05; 
  }
}
//=====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 * @param atolVector Reference for the atol vector to fill up
 */
void porousLiKCl_infPlate_dom1D::setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted & soln, 
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
    if (val < 1.0E-8) {
      val = 1.0E-8;
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
porousLiKCl_infPlate_dom1D::setAtolDeltaDamping(double atolDefault, double relcoeff, 
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
    if (val < 1.0E-6) {
      val = 1.0E-6;
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
porousLiKCl_infPlate_dom1D::setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff, 
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
porousLiKCl_infPlate_dom1D::err(const char *msg)
{
  printf("porousLiKCl_infPlate_dom1D: function not implemented: %s\n", msg);
  exit(-1);
}
//=====================================================================================================================
//=====================================================================================================================
}
//=====================================================================================================================


