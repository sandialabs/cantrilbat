/**
 * @file m1d_porousLiIon_Separator_dom1D.cpp
 */

/*
 * Copywrite 2013 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "LiIon_PorousBat.h"

#include "m1d_porousLiIon_Separator_dom1D.h"
#include "m1d_BDD_porSeparator_LiIon.h"
#include "m1d_NodalVars.h"
#include "m1d_GlobalIndices.h"
#include "m1d_ProblemStatementCell.h"
#include "m1d_SurfDomainDescription.h"
#include "m1d_CanteraElectrodeGlobals.h"

#include "cantera/transport/Tortuosity.h"

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

namespace m1d
{
//==================================================================================================================================
static void
drawline(int sp, int ll)
{
    for (int i = 0; i < sp; i++) {
        ZZCantera::writelog(" ");
    }
    for (int i = 0; i < ll; i++) {
        ZZCantera::writelog("-");
    }
    ZZCantera::writelog("\n");
}
//==================================================================================================================================
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
//==================================================================================================================================
porousLiIon_Separator_dom1D::porousLiIon_Separator_dom1D(BDD_porSeparator_LiIon* bdd_sep_ptr) :
    porousFlow_dom1D(bdd_sep_ptr),
    BDT_sep_ptr_(bdd_sep_ptr),
    nph_(0), nsp_(0),
    concTot_Cell_(0), 
    concTot_Cell_old_(0), 
    Fleft_cc_(0.0),
    Fright_cc_(0.0), Vleft_cc_(0.0), Vcent_cc_(0.0), Vright_cc_(0.0), 
    t_final_(0.0), t_init_(0.0), Xleft_cc_(0), Xcent_cc_(0), Xright_cc_(0),
    spCharge_(0),
    gradT_trCurr_(0.0), gradV_trCurr_(0), gradX_trCurr_(0), Vdiff_trCurr_(0), jFlux_trCurr_(0),
    icurrElectrolyte_CBL_(0), 
    icurrElectrolyte_CBR_(0),
    iECDMC_(-1),
    iLip_(-1),
    iPF6m_(-1),
    solnTemp(0)
{

    //BDT_ptr_ = dynamic_cast<BDT_porSeparator_LiIon*>(&BDD_);
    if (!BDT_sep_ptr_) {
        throw m1d_Error("porousLiIon_Separator_dom1D", "confused");
    }
    nsp_ = ionicLiquid_->nSpecies();
    nph_ = 1;


    iECDMC_ = ionicLiquid_->speciesIndex("ECDMC");
    if (iECDMC_ < 0) {
        throw CanteraError("confused", "confused");
    }
    iLip_ = ionicLiquid_->speciesIndex("Li+");
    if (iLip_ < 0) {
        throw CanteraError("confused", "confused");
    }
    iPF6m_ = ionicLiquid_->speciesIndex("PF6-");
    if (iPF6m_ < 0) {
        throw CanteraError("confused", "confused");
    }

}
//==================================================================================================================================
porousLiIon_Separator_dom1D::porousLiIon_Separator_dom1D(const porousLiIon_Separator_dom1D& r) :
    porousFlow_dom1D(r.BDT_sep_ptr_),
    BDT_sep_ptr_(r.BDT_sep_ptr_),
    nph_(0), nsp_(0),
    concTot_Cell_(0), concTot_Cell_old_(0),
    Fleft_cc_(0.0),
    Fright_cc_(0.0), Vleft_cc_(0.0), Vcent_cc_(0.0), Vright_cc_(0.0), 
    t_final_(0.0), t_init_(0.0), Xleft_cc_(0), Xcent_cc_(0), Xright_cc_(0),
    spCharge_(0),
    gradT_trCurr_(0.0), gradV_trCurr_(0), gradX_trCurr_(0), Vdiff_trCurr_(0), jFlux_trCurr_(0),
    icurrElectrolyte_CBL_(0), icurrElectrolyte_CBR_(0),
    iECDMC_(-1),
    iLip_(-1),
    iPF6m_(-1),
    solnTemp(0)
{
    porousLiIon_Separator_dom1D::operator=(r);
}
//==================================================================================================================================
porousLiIon_Separator_dom1D::~porousLiIon_Separator_dom1D()
{
}
//==================================================================================================================================
porousLiIon_Separator_dom1D&
porousLiIon_Separator_dom1D::operator=(const porousLiIon_Separator_dom1D& r)
{
    if (this == &r) {
        return *this;
    }
    // Call the parent assignment operator
    porousFlow_dom1D::operator=(r);

    BDT_sep_ptr_                          = r.BDT_sep_ptr_;
    ionicLiquid_                          = r.ionicLiquid_;
    nph_                                  = r.nph_;
    nsp_                                  = r.nsp_;
    concTot_Cell_                         = r.concTot_Cell_;
    concTot_Cell_old_                     = r.concTot_Cell_old_;

    Fleft_cc_                             = r.Fleft_cc_;
    Fright_cc_                            = r.Fright_cc_;
    Vleft_cc_                             = r.Vleft_cc_;
    Vcent_cc_                             = r.Vcent_cc_;
    Vright_cc_                            = r.Vright_cc_;
    t_final_                              = r.t_final_;
    t_init_                               = r.t_init_;
    Xleft_cc_                             = r.Xleft_cc_;
    Xcent_cc_                             = r.Xcent_cc_;
    Xright_cc_                            = r.Xright_cc_;
    spCharge_                             = r.spCharge_;
    mfElectrolyte_Soln_Cell_old_          = r.mfElectrolyte_Soln_Cell_old_;
    gradT_trCurr_                         = r.gradT_trCurr_;
    gradV_trCurr_                         = r.gradV_trCurr_;
    gradX_trCurr_                         = r.gradX_trCurr_;
    Vdiff_trCurr_                         = r.Vdiff_trCurr_;
    jFlux_trCurr_                         = r.jFlux_trCurr_;
    icurrElectrolyte_CBL_                 = r.icurrElectrolyte_CBL_;
    icurrElectrolyte_CBR_                 = r.icurrElectrolyte_CBR_;
    iECDMC_                               = r.iECDMC_;
    iLip_                                 = r.iLip_;
    iPF6m_                                = r.iPF6m_;
    solnTemp                              = r.solnTemp;

    return *this;
}
//==================================================================================================================================
// Prepare all of the indices for fast calculation of the residual
/*
 *  Ok, at this point, we will have figured out the number of equations
 *  to be calculated at each node point. The object NodalVars will havex1
 *  been fully formed.
 *
 *  We use this to figure out what local node numbers/ cell numbers are
 *  needed and to set up indices for their efficient calling.
 *
 *  Child objects of this one will normally call this routine in a
 *  recursive fashion.
 */
void
porousLiIon_Separator_dom1D::domain_prep(LocalNodeIndices* li_ptr)
{
    /*
     * First call the parent domain prep to get the node information
     */
    porousFlow_dom1D::domain_prep(li_ptr);

    double domainThickness = BDT_ptr_->Xpos_end - BDT_ptr_->Xpos_start;
    double porosity = -1.0;
    double grossVolumeSeparator = PSCinput_ptr->separatorArea_ * domainThickness;
    double volumeInert = 0.0;
    double volumeFractionInert = 0.0;
    if (solidSkeleton_) {
        if (PSCinput_ptr->separatorMass_ > 0.0) {
            volumeInert = PSCinput_ptr->separatorMass_ / solidSkeleton_->density() ;
	    volumeFractionInert = volumeInert / grossVolumeSeparator;
        } else if (PSCinput_ptr->separatorSolid_vf_ > 0.0) {
	    volumeFractionInert = PSCinput_ptr->separatorSolid_vf_;
        }
        if (volumeFractionInert >= 1.0) {
            throw m1d_Error("porousLiIon_Separator_dom1D::domain_prep()",
                            "Volume fraction of separator solid more than 1.0: " + fp2str(volumeFractionInert));
        }
    }
    //
    // If there is a solidSkeleton ThermoPhase, then identify that with the first volume fraction of the extra condensed phases.
    // We'll keep the mole number and volume fraction in the extra phases lists.
    //
   
    if (numExtraCondensedPhases_ > 0) {
	for (size_t iCell = 0; iCell < (size_t) NumLcCells; ++iCell) {
            int  index_CentLcNode = Index_DiagLcNode_LCO[iCell];
            NodalVars* nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
            NodalVars* nodeLeft = nodeCent;
            NodalVars* nodeRight = nodeCent;
            if (iCell !=0) {
                   int index_LeftLcNode = Index_LeftLcNode_LCO[iCell];
               nodeLeft =  LI_ptr_->NodalVars_LcNode[ index_LeftLcNode];
            }
            if ((int) iCell !=NumLcCells - 1) {
                   int index_RightLcNode = Index_RightLcNode_LCO[iCell];
               nodeRight = LI_ptr_->NodalVars_LcNode[index_RightLcNode];
            }
            double cellThickness = 0.5*(nodeRight->xNodePos() - nodeLeft->xNodePos());
	    porosity = 1.0;
	    for (size_t k = 0; k < ExtraPhaseList_.size(); ++k) {
		ExtraPhase* ep = ExtraPhaseList_[k];
		ThermoPhase* tp = ep->tp_ptr;
		tp->setState_TP(temp_Curr_, pres_Curr_);
		double mvp = tp->molarVolume();
		volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + k] = ep->volFraction;
		volumeFraction_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + k] = ep->volFraction;
		moleNumber_Phases_Cell_[numExtraCondensedPhases_ * iCell + k] = ep->volFraction * cellThickness * mvp ;
		moleNumber_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + k] = ep->volFraction * cellThickness * mvp;
		porosity -= ep->volFraction;
	    }
	    porosity_Cell_[iCell] = porosity;
	    porosity_Cell_old_[iCell] = porosity;
	}
    } else {
	for (size_t iCell = 0; iCell < (size_t) NumLcCells; ++iCell) {
	    porosity_Cell_[iCell] = 1.0;
	    porosity_Cell_old_[iCell] = 1.0;
	}
    }

    /*
     * Porous electrode domain prep
     */
    concTot_Cell_.resize(NumLcCells, 0.0);
    concTot_Cell_old_.resize(NumLcCells, 0.0);

    Xleft_cc_.resize(nsp_, 0.0);
    Xcent_cc_.resize(nsp_, 0.0);
    Xright_cc_.resize(nsp_, 0.0);

    spCharge_.resize(nsp_, 0.0);
    for (int k = 0; k < nsp_; k++) {
        spCharge_[k] = ionicLiquid_->charge(k);
    }

    gradX_trCurr_.resize(nsp_, 0.0);
    Vdiff_trCurr_.resize(nsp_, 0.0);
    jFlux_trCurr_.resize(nsp_, 0.0);

    solnTemp.resize(10, 0.0);

    icurrElectrolyte_CBL_.resize(NumLcCells, 0.0);
    icurrElectrolyte_CBR_.resize(NumLcCells, 0.0);

    mfElectrolyte_Soln_Cell_old_.resize(nsp_, NumLcCells, 0.0);

    /*
     *  Set the velocity basis of the transport object. We are using
     *  mole-averaged velocities as the basis.
     */
    trans_->setVelocityBasis(ivb_);
#ifdef TRACK_LOCATION
    std::cout << "Leaving porousLiIon_Separator_dom1D::domain_prep"<<std::endl;
#endif

}
//==================================================================================================================================
// Function that gets called at end the start of every time step
/*
 *  This function provides a hook for a residual that gets called whenever a
 *  time step has been accepted and we are about to move on to the next time step.
 *  The call is made with the current time as the time
 *  that is accepted. The old time may be obtained from t and rdelta_t_accepted.
 *
 *  After this call interrogation of the previous time step's results will not be
 *  valid.
 *
 *   @param  doTimeDependentResid  This is true if we are solving a time dependent
 *                                 problem.
 *   @param  soln_ptr              Solution value at the current time
 *   @param  solnDot_ptr           derivative of the solution at the current time.
 *   @param  solnOld_ptr           Solution value at the old time step, n-1
 *   @param  t                     current time to be accepted, n
 *   @param  t_old                 previous time step value
 */
void
porousLiIon_Separator_dom1D::advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector* soln_ptr,
                                                 const Epetra_Vector* solnDot_ptr, const Epetra_Vector* solnOld_ptr,
                                                 const double t, const double t_old)
{
    const Epetra_Vector& soln = *soln_ptr;
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell;
        /*
         *  ---------------- Get the index for the center node ---------------------------------
         */
        int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
        /*
         *   Get the pointer to the NodalVars object for the center node
         */
        NodalVars* nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
        /*
         *  Index of the first equation in the bulk domain of center node
         */
        int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];

        SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));

        concTot_Cell_old_[iCell] = concTot_Curr_;
        porosity_Cell_old_[iCell] = porosity_Cell_[iCell] ;
	Temp_Cell_old_[iCell] = temp_Curr_;

        double* mfElectrolyte_Soln_old = mfElectrolyte_Soln_Cell_old_.ptrColumn(iCell);
        for (size_t k = 0; k < (size_t) nsp_; ++k) {
            mfElectrolyte_Soln_old[k] = mfElectrolyte_Soln_Curr_[k];
        }

        if (energyEquationProbType_ == 3) { 
	        double volCellNew = xdelCell_Cell_[iCell];
		/*
                double solidMolarEnthalpyNew = solidSkeleton_->enthalpy_mole();
                double solidConcNew =  1.0 / solidSkeleton_->molarVolume();
                double volSolidNew = (1.0 - porosity_Curr_) * volCellNew;
                double solidEnthalpyNew = solidMolarEnthalpyNew * solidConcNew *  volSolidNew;
		*/
                double lyteMolarEnthalpyNew = ionicLiquid_->enthalpy_mole();
                double volLyteNew = porosity_Curr_ * volCellNew;
                double lyteEnthalpyNew =  lyteMolarEnthalpyNew * concTot_Curr_ * volLyteNew;
		double nEnthalpyInertNew = 0.0;
		for (size_t jPhase = 0; jPhase < numExtraCondensedPhases_; ++jPhase) {
		    ExtraPhase* ep = ExtraPhaseList_[jPhase];
		    ThermoPhase* tp = ep->tp_ptr;
		    tp->setState_TP(temp_Curr_, pres_Curr_);
		    double mEnth = tp->enthalpy_mole();
		    double mn = moleNumber_Phases_Cell_[numExtraCondensedPhases_ * iCell + jPhase];
		    nEnthalpyInertNew += mEnth * mn;
		}
		
                double nEnthalpy_New  = lyteEnthalpyNew + nEnthalpyInertNew;

		if (! checkDblAgree( nEnthalpy_New, nEnthalpy_New_Cell_[iCell] ) ) {
                    throw m1d_Error("porousLiIon_Separator_dom1D::advanceTimeBaseline", 
                                    "Disagreement on new enthalpy calc");
                }
		
		nEnthalpy_Old_Cell_[iCell] = nEnthalpy_New; 
        }
        if (numExtraCondensedPhases_ > 0) {
	for (size_t jPhase = 0; jPhase < numExtraCondensedPhases_; jPhase++) {
	    moleNumber_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + jPhase] =  
	      moleNumber_Phases_Cell_[numExtraCondensedPhases_ * iCell + jPhase];
	    volumeFraction_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + jPhase] =  
	      volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + jPhase];  
	}
        }
    }
}
//==================================================================================================================================
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
porousLiIon_Separator_dom1D::residEval(Epetra_Vector& res,
                                       const bool doTimeDependentResid,
                                       const Epetra_Vector* soln_ptr,
                                       const Epetra_Vector* solnDot_ptr,
                                       const Epetra_Vector* solnOld_ptr,
                                       const double t,
                                       const double rdelta_t,
                                       const Zuzax::ResidEval_Type residType,
                                       const Zuzax::Solve_Type solveType)
{
  static int tmpsSetup = 0;
  if (!tmpsSetup) {
    residSetupTmps();
    tmpsSetup = 1;
  }
  residType_Curr_ = residType;
  // convert to static_cast!!
  //BDT_porSeparator_LiIon* fa = dynamic_cast<BDT_porSeparator_LiIon*>(&BDD_);
  t_final_ = t;
  if (rdelta_t < 1.0E-200) {
    t_init_ = t;
  } else {
    t_init_ = t - 1.0/rdelta_t;
  }

  NodalVars* nodeLeft = 0;
  NodalVars* nodeCent = 0;
  NodalVars* nodeRight = 0;

  double xdelL; // Distance from the center node to the left node
  double xdelR; // Distance from the center node to the right node
  double xdelCell; // cell width - right boundary minus the left boundary.

  //  Electrolyte mass fluxes - this is rho V dot n at the boundaries of the cells
  double moleFluxRight = 0.0;
  double moleFluxLeft = 0.0;
  //  Thermal fluxes
  double fluxTright = 0.0;
  double fluxTleft = 0.0;
  double fluxL_JHPhi = 0.0;
  double fluxR_JHPhi = 0.0;
  double enthConvRight = 0.0;
  double enthConvLeft = 0.0;
  double oldStuffTC = 0.0;
  double newStuffTC = 0.0;
  double* mf_old = 0;

  //mole fraction fluxes
  std::vector<double> fluxXright(nsp_, 0.0);
  std::vector<double> fluxXleft(nsp_, 0.0);
 
  std::vector<double> resLocal_Species(nsp_, 0.0);

  double res_Cont_0 = 0.0;
 
  const Epetra_Vector& soln = *soln_ptr;
  const Epetra_Vector& solnOld = *solnOld_ptr;

#ifdef MECH_MODEL
  // scratch pad to calculate the expansion of the matrix
  std::vector<double> xratio(NumLcCells,0.0); 
  
  // scratch pad to calculate the new positions due to the strain of the elements
  std::vector<double> new_node_pos(NumLcCells,0.0);

  // for the seperator, we use the following values:
  // These numbers are only 'good' to ~8-20% strain. 
  // v poisson's ratio             0.5 @28.5C 0.46 @80C
  // K  Bulk Modulus           615 Mpa @ 28.5C,  481 MPa @ 55C, 327MPa @ 80C, all saturated with liquid. 
  // youngs mod = 3K(1-2v) 
  // sher mod G = 3*K(1-2v)/2(1+v)
  // extracted from graphs on page 59-60 of http://dc.uwm.edu/cgi/viewcontent.cgi?article=1039&context=etd
  //Martinsen, Michael James, "Material Behavior Characterization of a Thin Film Polymer Used in Lithium-Ion Batteries" (2012).Theses and Dissertations.Paper 36
  // thermal expansion coefficient is 35-110 e-6 depending on what kind of mechanical process the seperator polyethylene has be subject to. Until better info exists, pick 50e-6 
  double Thermal_Expansion = 50e-6;
  double poisson = 0.45;
  double Eyoung = -1;

  int pMECH_MODEL_extra = 1;
#endif
  /*
   * Index of the first equation at the left node corresponding to the first bulk domain, which is the electrolyte
   */
  size_t indexLeft_EqnStart;
  /*
   * Index of the first equation at the center node corresponding to the first bulk domain, which is the electrolyte
   */
  size_t indexCent_EqnStart;
  /*
   * Index of the first equation at the right node corresponding to the first bulk domain, which is the electrolyte
   */
  size_t indexRight_EqnStart;

  incrementCounters(residType);
  Fright_cc_ = 0.0;
  /*
   *  -------------------- Special section to do the left boundary -------------------------------
   */

  
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

    cellTmps& cTmps      = cellTmpsVect_Cell_[iCell];
    ///valCellTmps& valTmps = valCellTmpsVect_Cell_[iCell];

    NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
    NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;
    NodeTmps& nodeTmpsRight  = cTmps.NodeTmpsRight_;

#ifdef DEBUG_HKM_NOT
    if (counterResBaseCalcs_ > 125 && residType == Zuzax::ResidEval_Type::Base_ResidEval) {
      if (iCell == NumLcCells - 1) {
	// printf("we are here\n");
      }
    }
#endif
    

    /*
     *  ---------------- Get the index for the center node ---------------------------------
     *   Get the pointer to the NodalVars object for the center node
     *   Index of the first equation in the bulk domain of center node
     */
    nodeCent = cTmps.nvCent_;
    indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;

    /*
     *  ------------------- Get the index for the left node -----------------------------
     *    There may not be a left node if we are on the left boundary. In that case
     *    set the pointer to zero and the index to -1.
     *    The solution index is set to the center solution index in that case as well.
     */
    nodeLeft = cTmps.nvLeft_;
    indexLeft_EqnStart = nodeTmpsLeft.index_EqnStart;
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
    nodeRight = cTmps.nvRight_;
    indexRight_EqnStart = nodeTmpsRight.index_EqnStart;

    /*
     * --------------------------- CALCULATE POSITION AND DELTA_X Variables -----------------------------
     * Calculate the distance between the left and center node points
     */
    xdelL = cTmps.xdelL_;
    /*
     * Calculate the distance between the right and center node points
     */
    xdelR = cTmps.xdelR_;
    /*
     * Calculate the cell width
     */
    xdelCell = cTmps.xdelCell_; 
    xdelCell_Cell_[iCell] = xdelCell;

    /*
     * --------------------------- DO PRE-SETUPSHOP RASTER OVER LEFT,CENTER,RIGHT -----------------------------
     * Calculate the distance between the left and center node points
     */
    /*
     * Get current velocity, mole fraction, temperature, potential
     * from the solution
     */

    for (int k = 0; k < nsp_; k++) {
      Xcent_cc_[k] = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_MoleFraction_Species + k];
    }
    Vcent_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Voltage];

    if (nodeLeft != 0) {
      /*
       * Find the velocity located at the left cell boundary.
       * The left cell boundary velocity is stored at the previous (left)
       * cell index as per our conventions.
       */
      Fleft_cc_ = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Velocity_Axial];
      for (int k = 0; k < nsp_; k++) {
	Xleft_cc_[k] = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_MoleFraction_Species + k];
      }
      Vleft_cc_ = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Voltage];
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
    Fright_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Velocity_Axial];

    if (nodeRight != 0) {

      for (int k = 0; k < nsp_; k++) {
	Xright_cc_[k] = soln[indexRight_EqnStart + nodeTmpsRight.Offset_MoleFraction_Species + k];
      }
      Vright_cc_ = soln[indexRight_EqnStart + nodeTmpsRight.Offset_Voltage];
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
    // We've turned this off, because it is not necessary
    if (doLeftFluxCalc) {
      if (nodeLeft == 0) {
	AssertTrace(iCell == 0);
	/*
	 *  We are here if we are at the left node boundary and we
	 *  need a flux condition. The default now is to
	 *  set the flux to zero. This is good if we expect the negate/equalize the fluxes at a domain boundary,
	 *  since both sides are set to zero. We could put in a more
	 *  sophisticated treatment.
	 *
	 *  The other case is when we are specifying a flux boundary condition. In that case,
	 *  we set the flux here to zero. Then, in the surface domain we add the specification of
	 *  the flux to the residual of the node at the boundary.
	 */
	moleFluxLeft = 0.0;
	fluxTleft = 0.0;
	fluxL_JHPhi = 0.0;
	enthConvLeft = 0.0;
	icurrElectrolyte_CBL_[iCell] = 0.0;
	for (int k = 0; k < nsp_; k++) {
	  fluxXleft[k] = 0.0;
	}
      } else {

	/*
	 *  Establish the environment at the left cell boundary
	 */
	SetupThermoShop2(nodeLeft, &(soln[indexLeft_EqnStart]), nodeCent, &(soln[indexCent_EqnStart]), 0);

	SetupTranShop(xdelL, 0);

	/*
	 * Calculate the flux at the left boundary for each equation
	 */
	moleFluxLeft = Fleft_cc_ * concTot_Curr_;
	fluxTleft = heatFlux_Curr_;
	fluxL_JHPhi = jFlux_EnthalpyPhi_Curr_;
	enthConvLeft = moleFluxLeft * EnthalpyMolar_lyte_Curr_;

	/*
	 * Calculate the flux of species and the flux of charge
	 *   - the flux of charge must agree with the flux of species
	 *   - At the left boundary, we've set Fleft_cc_ to zero. This is necessary if the domain isn't the left-most one.
	 *     In that case the cell is shared between half-cells on either side of the domain. So, we will want to 
	 *     not add any convection term here.
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
	icurrElectrolyte_CBL_[iCell] *= (ZZCantera::Faraday);
      }
    } else {
      /*
       * Copy the fluxes from the stored right side of the previous cell to the left side of the current cell
       */
      moleFluxLeft = moleFluxRight;
      fluxTleft = fluxTright;
      fluxL_JHPhi = fluxR_JHPhi;
      enthConvLeft = enthConvRight;
      icurrElectrolyte_CBL_[iCell] = icurrElectrolyte_CBR_[iCell - 1];
      for (int k = 0; k < nsp_; k++) {
	fluxXleft[k] = fluxXright[k];
      }
    }
    /*
     * ------------------- CALCULATE FLUXES AT THE RIGHT BOUNDARY -------------------------------
     *
     */
    if (nodeRight == 0) {
      AssertTrace(iCell == NumLcCells-1);
      SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
      SetupThermoShop1Extra(nodeCent, &(soln[indexCent_EqnStart]));
      /*
       *  We are here if we are at the right node boundary and we need a flux
       *  condition. The default now is to set the flux to zero. 
       *  This is good if we expect the negate/equalize the fluxes at a domain boundary,
       *  since both sides are set to zero.
       *
       *  The other case is when we are specifying a flux boundary condition. In that case,
       *  we set the flux here to zero. Then, in the surface domain we add the specification of
       *  the flux to the residual of the node at the boundary.
       *
       *  Fright_cc_ is set to zero, because we are at an internal boundary. Convection terms are set to zero at these 
       *  types of internal boundaries.
       */
      Fright_cc_ = 0.0;
      fluxTright = 0.0;
      fluxR_JHPhi = 0.0;
      enthConvRight = 0.0;
      moleFluxRight = Fright_cc_ * concTot_Curr_;
      icurrElectrolyte_CBR_[iCell] = 0.0;
      for (int k = 0; k < nsp_; k++) {
	fluxXright[k] = 0.0;
      }
    } else {
      /*
       *  Establish the environment at a normal right cell boundary
       */
      SetupThermoShop2(nodeCent, &(soln[indexCent_EqnStart]), nodeRight, &(soln[indexRight_EqnStart]), 1);

      SetupTranShop(xdelR, 1);

      /*
       * Calculate the flux at the right boundary for each equation
       * This is equal to
       *     fluxFright =  Conc * Vaxial * phi
       */
      moleFluxRight = Fright_cc_ * concTot_Curr_;
      /*
       * Calculate the heat flux - all of the types
       */
      fluxTright = heatFlux_Curr_;
      fluxR_JHPhi = jFlux_EnthalpyPhi_Curr_;
      enthConvRight = moleFluxRight * EnthalpyMolar_lyte_Curr_;
           
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
      icurrElectrolyte_CBR_[iCell] *= (ZZCantera::Faraday);
    }

#ifdef DEBUG_HKM_NOT
    if (doTimeDependentResid) {

      if (residType == Zuzax::ResidEval_Type::Base_ResidEval) {
	printf(" Cell = %d, Totalflux_K+ = %10.3e,  Totalflux_Cl- = %10.3e \n", iCell, fluxXright[1], fluxXright[2]);
	printf("           Vmolal = %10.3e, jd_Li+ = %10.3e  jd_K+ = %10.3e jd_Cl- = %10.3e\n", Fright_cc_,
	       jFlux_trCurr_[0], jFlux_trCurr_[1], jFlux_trCurr_[2]);
	printf("           Vmolal = %10.3e, vd_Li+ = %10.3e  vd_K+ = %10.3e vd_Cl- = %10.3e\n", Fright_cc_,
	       Vdiff_trCurr_[0], Vdiff_trCurr_[1], Vdiff_trCurr_[2]);
      }
    }
#endif

    /*
     * -------------------- ADD THE FLUX TERMS INTO THE RESIDUAL EQUATIONS -----------------------------------------------
     */

#ifdef DEBUG_RESID
    double residBefore = 0.0;
    if (IOwnLeft && iCell == 0) {
      if (residType == Zuzax::ResidEval_Type::Base_ShowSolution) {
	residBefore =  res[indexCent_EqnStart + nodeTmpsCent.RO_Electrolyte_Continuity];
	double tmp3 = 0.0;
	double sum5 = residBefore + fluxFright - tmp3;
	printf("ResidSep Cell = %d, ResBefore = -fluxleft = %10.3e, FluxRight = %10.3e,"
	       "Prod = %10.3e total = %10.3e \n", iCell,
	       residBefore, fluxFright, tmp3, sum5);

      }
    }
#endif
    /*
     *  Total continuity equation - fluxFright and fluxFleft represent the total mass
     *                              fluxes coming and going from the cell.
     *                    R =   d rho d t + del dot (rho V) = 0
     */
    if (ivb_ == VB_MOLEAVG) {
      res[indexCent_EqnStart + nodeTmpsCenter.RO_Electrolyte_Continuity] += (moleFluxRight - moleFluxLeft);

	   
    } else {
      exit(-1); 
    }

    /*
     * Species continuity Equation
     */
    for (int k = 0; k < nsp_; k++) {
      if ((k != (int) BDT_sep_ptr_->iMFS_index_) && (k != iPF6m_)) {
	res[indexCent_EqnStart + nodeTmpsCenter.RO_Species_Eqn_Offset + k] += (fluxXright[k] - fluxXleft[k]);
      }
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
    res[indexCent_EqnStart + nodeTmpsCenter.RO_Current_Conservation] +=
      (icurrElectrolyte_CBR_[iCell] - icurrElectrolyte_CBL_[iCell]);
    /*
     *  Energy Equation
     */
    if (PS_ptr->energyEquationProbType_ == 3) {
      AssertTrace(nodeTmpsCenter.RO_Enthalpy_Conservation != npos);
      res[indexCent_EqnStart + nodeTmpsCenter.RO_Enthalpy_Conservation] += (fluxTright - fluxTleft);
      res[indexCent_EqnStart + nodeTmpsCenter.RO_Enthalpy_Conservation] += (fluxR_JHPhi - fluxL_JHPhi);
      res[indexCent_EqnStart + nodeTmpsCenter.RO_Enthalpy_Conservation] += (enthConvRight - enthConvLeft);
    }
    /*
     * --------------------------------------------------------------------------
     *             Add in the source terms at the current cell center
     * --------------------------------------------------------------------------
     */
    SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));

    /*
     *  Total continuity equation -
     */

    /*
     * Mole fraction summation equation
     */
    res[indexCent_EqnStart + nodeTmpsCenter.RO_MFSum_offset] = 1.0;
    for (int k = 0; k < nsp_; k++) {
      res[indexCent_EqnStart + nodeTmpsCenter.RO_MFSum_offset] -= Xcent_cc_[k];
    }

    /*
     * Electroneutrality equation
     */
    res[indexCent_EqnStart + nodeTmpsCenter.RO_ChargeBal_offset] = 0.0;
    for (int k = 0; k < nsp_; k++) {
      res[indexCent_EqnStart + nodeTmpsCenter.RO_ChargeBal_offset] += Xcent_cc_[k] * spCharge_[k];
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
      if (residType == Zuzax::ResidEval_Type::Base_ShowSolution) {
	icurrElectrolyte_CBL_[iCell] =  icurrElectrolyte_CBR_[iCell];
      }
    }
    /*
     * Special section if we own the right node of the domain. If we do
     * it will be cell 0. Store currents for later use.
     * These are the correct currents that work for the global balances
     */
    if (IOwnRight && iCell == (NumLcCells - 1)) {
      if (residType == Zuzax::ResidEval_Type::Base_ShowSolution) {
	icurrElectrolyte_CBR_[iCell] =  icurrElectrolyte_CBL_[iCell];
      }
    }
    /*
     *   ------------------ ADD IN THE TIME DEPENDENT TERMS ----------------------------------------------------
     */
    if (doTimeDependentResid) {

#ifdef DEBUG_HKM_NOT
      if (residType == Zuzax::ResidEval_Type::Base_ResidEval) {
	printf(" Cell = %d, Totalflux_Li+_r = %10.3e,  = %10.3e, Totalflux_Li+_l ", iCell, fluxXright[0], fluxXleft[0]);
      }
#endif
      newStuffTC = concTot_Curr_ * porosity_Curr_ * xdelCell;
      double newStuffSpecies0 = Xcent_cc_[iLip_] * newStuffTC;
      /*
       *   .................... Calculate quantities needed at the previous time step
       */
      /*
       * Setup shop with the old time step
       */
      SetupThermoShop1(nodeCent, &(solnOld[indexCent_EqnStart]));

      mf_old = mfElectrolyte_Soln_Cell_old_.ptrColumn(iCell);
      oldStuffTC = concTot_Cell_old_[iCell] * porosity_Cell_old_[iCell] * xdelCell;

      // oldStuffTC = concTot_Curr_ * porosity_Curr_ * xdelCell;
      double oldStuffSpecies0 = mfElectrolyte_Soln_Curr_[iLip_] * oldStuffTC;

      //AssertTrace(mf_old[iLip_] ==  mfElectrolyte_Soln_Curr_[iLip_]); 
      //(HKM Note, this uncovered a print error in time stepper reentry
      //if (mf_old[iLip_] !=  mfElectrolyte_Soln_Curr_[iLip_]) {
      //printf("We are here\n");
      // }
      // double tmp = (newStuffSpecies0 - oldStuffSpecies0) * rdelta_t;
    

      /*
       *   .................... Add these terms in the residual
       */
      /*
       *  Add in the time term for species iLip
       */

      for (int k = 0; k < nsp_; k++) {
	if (k != iECDMC_ && k != iPF6m_) {
	  double tmp = (newStuffSpecies0 - oldStuffSpecies0) * rdelta_t;
	  res[indexCent_EqnStart + nodeTmpsCenter.RO_Species_Eqn_Offset + k] += tmp;
	}
      }

      /*
       *   Add in the time term for the total continuity equation
       *         note: the current problem will have this term equally zero always.
       *               However, we put it in here for the next problem and for consistency of the time terms
       */
      res[indexCent_EqnStart + nodeTmpsCenter.RO_Electrolyte_Continuity] += (newStuffTC - oldStuffTC) * rdelta_t;

      if (iCell == 0) {
	res_Cont_0 =  (newStuffTC - oldStuffTC) * rdelta_t + (moleFluxRight - moleFluxLeft);
      }



      if  (energyEquationProbType_ == 3) {
	//
	// Do old time enthalpy calculation -> will be replaced
	//   (volcell might have changed -> this will work if it hasn't)
	  /*
	double volCellOld = xdelCell;
	double solidMolarEnthalpyOld = solidSkeleton_->enthalpy_mole();
	double solidConcOld =  1.0 / solidSkeleton_->molarVolume();
	double volSolidOld = (1.0 - porosity_Curr_) * volCellOld;
	double solidEnthalpyOld = solidMolarEnthalpyOld * solidConcOld *  volSolidOld;

	double lyteMolarEnthalpyOld = ionicLiquid_->enthalpy_mole();
	double volLyteOld = porosity_Curr_ * volCellOld;
	double lyteEnthalpyOld =  lyteMolarEnthalpyOld * concTot_Curr_ * volLyteOld;
	
	double nEnthalpy_Old =  solidEnthalpyOld + lyteEnthalpyOld;


	if (! checkDblAgree( nEnthalpy_Old, nEnthalpy_Old_Cell_[iCell] ) ) {
	  throw m1d_Error("", "Disagreement on old enthalpy calc");
	}
	  */

	SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));

	double lyteMolarEnthalpyNew = ionicLiquid_->enthalpy_mole();
	double volLyteNew = porosity_Curr_ * xdelCell;
	double lyteEnthalpyNew =  lyteMolarEnthalpyNew * concTot_Curr_ * volLyteNew;
	// 
	// Calculate the enthalpy of other phases that may be present
	//
	double nEnthalpyInertNew = 0.0;
	for (size_t jPhase = 0; jPhase < numExtraCondensedPhases_; ++jPhase) {
	    ExtraPhase* ep = ExtraPhaseList_[jPhase];
	    ThermoPhase* tp = ep->tp_ptr;
	    tp->setState_TP(temp_Curr_, pres_Curr_);
	    double mEnth = tp->enthalpy_mole();
	    double mn = moleNumber_Phases_Cell_[numExtraCondensedPhases_ * iCell + jPhase];
	    nEnthalpyInertNew += mEnth * mn;
	}
	
	nEnthalpy_New_Cell_[iCell] = lyteEnthalpyNew + nEnthalpyInertNew;

	double tmp = (nEnthalpy_New_Cell_[iCell] - nEnthalpy_Old_Cell_[iCell]) * rdelta_t;
	res[indexCent_EqnStart + nodeTmpsCenter.RO_Enthalpy_Conservation] += tmp;
      }

      // ----------------------  Calculations that are based on having a complete local residual -----------------------------

      /*
       *   .................... Go back to setting up shop at the current time
       */
      SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
    }
    //
    //  Section to find the Axial velocity at the left domain boundary, but only for show and only when the 
    //          residual is already solved.
    //
    if (residType == Zuzax::ResidEval_Type::Base_ShowSolution) {
      if (IOwnLeft) {
	if (iCell == 0) {
	  res_Cont_0 = (newStuffTC - oldStuffTC) * rdelta_t + (moleFluxRight );
	  moleFluxLeft = res_Cont_0;
	  Fleft_cc_ = moleFluxLeft / concTot_Curr_;
	  DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity] = Fleft_cc_;
	  //
	  //  Section to find the Diffusive flux of species
	  //        fluxXleft[k] is set to zero. Also some species are replaced by constraints.
	  //        Here, we solve for everything, revealing the hidden and satisfied missing species continuity equations
	  for (size_t k = 0; k < (size_t) nsp_; k++) {
	    resLocal_Species[k] = (fluxXright[k] - fluxXleft[k]);
	    double newStuffSpecies0 = Xcent_cc_[k] * newStuffTC;
	    double oldStuffSpecies0 = mf_old[k] * oldStuffTC;
	    resLocal_Species[k] += (newStuffSpecies0 - oldStuffSpecies0) * rdelta_t;
	    //
	    // Save the local domain residual to a vector for later validation studies
	    //
	    DomainResidVectorLeftBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k] = resLocal_Species[k];

	    TotalFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k] = - resLocal_Species[k];
	    DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k] =  
	      TotalFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k]
	      -  Fleft_cc_ * mfElectrolyte_Soln_Curr_[k] * concTot_Curr_;
	  }
	  if  (energyEquationProbType_ == 3) {
	    double resLocal_Enth = (fluxTright - fluxTleft);
	    resLocal_Enth += (fluxR_JHPhi - fluxL_JHPhi);
	    resLocal_Enth += (enthConvRight - enthConvLeft);
	    resLocal_Enth += (nEnthalpy_New_Cell_[iCell] - nEnthalpy_Old_Cell_[iCell]) * rdelta_t;
	    DomainResidVectorLeftBound_LastResid_NE[nodeTmpsCenter.RO_Enthalpy_Conservation] = resLocal_Enth;
	  }
	}
      }
      //
      //  Section to find the Axial velocity at the right domain boundary
      //
      if (IOwnRight) {
	if (iCell == NumLcCells - 1) {
	  SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
	  SetupThermoShop1Extra(nodeCent, &(soln[indexCent_EqnStart]));
	  //
	  // Calculate the axial velocity that is needed to zero out the residual in this domain only
	  //
	  res_Cont_0 = (newStuffTC - oldStuffTC) * rdelta_t + ( - moleFluxLeft );
	  moleFluxRight = - res_Cont_0;
	  Fright_cc_ = moleFluxRight / concTot_Curr_;
	  DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity] = Fright_cc_;
	  //
	  //  Section to find the Diffusive flux of species
	  //        fluxXright[k] is set to zero. Also some species are replaced by constraints.
	  //        Here, we solve for everything, revealing the hidden and satisfied missing species continuity equations
	  //
	  for (size_t k = 0; k < (size_t) nsp_; k++) {
	    resLocal_Species[k] = (fluxXright[k] - fluxXleft[k]);
	    double newStuffSpecies0 =  Xcent_cc_[k] * newStuffTC;
	    double oldStuffSpecies0 =  mf_old[k] * oldStuffTC;
	    resLocal_Species[k] += (newStuffSpecies0 - oldStuffSpecies0) * rdelta_t;
	    //
	    // Save the local domain residual to a vector for later validation studies
	    //
	    DomainResidVectorRightBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k] = resLocal_Species[k];

	    TotalFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k] = - resLocal_Species[k];
	    fluxXright[k] = - resLocal_Species[k];
	    DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k] =  
	      TotalFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k]
	      -  Fright_cc_ * mfElectrolyte_Soln_Curr_[k] * concTot_Curr_;
	  }

	  if  (energyEquationProbType_ == 3) {
	      double resLocal_Enth = (fluxTright - fluxTleft);
	      resLocal_Enth += (fluxR_JHPhi - fluxL_JHPhi);
	      resLocal_Enth += (enthConvRight - enthConvLeft);
	      resLocal_Enth += (nEnthalpy_New_Cell_[iCell] - nEnthalpy_Old_Cell_[iCell]) * rdelta_t;
	      DomainResidVectorRightBound_LastResid_NE[nodeTmpsCenter.RO_Enthalpy_Conservation] = resLocal_Enth;
	  }
	}
      }
	       
    }
#ifdef MECH_MODEL
    if (solidMechanicsProbType_ > 0) {
      valCellTmps& valTmps = valCellTmpsVect_Cell_[iCell];
      // need the average temp to get the correct K
      cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
      NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
      NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;	
      NodalVars* nodeCent = cTmps.nvCent_;
      NodalVars* nodeLeft = cTmps.nvLeft_;
      NodalVars* nodeRight = cTmps.nvRight_;

      indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;
      indexLeft_EqnStart = nodeTmpsLeft.index_EqnStart;
      indexRight_EqnStart = nodeTmpsRight.index_EqnStart;
      nodeTmpsCenter.Offset_Displacement_Axial   = nodeCent->indexBulkDomainVar0((size_t) Displacement_Axial);
   
      double ReferenceTemperature =    valTmps.Temperature.center;
      
      double BulkMod=-1;
      if(ReferenceTemperature < 28.5+273.15 ) BulkMod = 615;
      else if(ReferenceTemperature > 80+273.15 ) BulkMod = 327;
      else if( ReferenceTemperature < 55+273.15) {
	double low = ReferenceTemperature - (28.5+273.14);
	BulkMod= 615 - ((615-481)*low/(55-28.5));
      }
      else if (ReferenceTemperature >= 55+273.15) {
	double low = ReferenceTemperature - (55+273.15);
	BulkMod = 481 - (481-327)*low/(80-55);
      }

      BulkMod*=1.0e6; // in MPa
      Eyoung =       3*BulkMod*(1.0 - 2*poisson);

      xratio[iCell] = 1.0;

      // CBL:: There is no change in the solid volume of the separator other than temperature effects???? \todo
      if ( (Domain1D::ChemEx | Domain1D::TempEx | Domain1D::All)  & solidMechanicsProbType_) {
	if(iCell <  NumLcCells-1) {
	  xratio[iCell] = 1.0;
	  nodeTmpsRight.Offset_Displacement_Axial    = nodeRight->indexBulkDomainVar0((size_t) Displacement_Axial);
	  xratio[iCell] =  (Thermal_Expansion+1.0)*(valTmps.Temperature.center+valTmps.Temperature.right)*0.5/ TemperatureReference_;
	  //*****************
	  //*****************
	  //***** NOTE *****
	  // Do _not_ use xNodePos() in any calculation of the displacement, as  
	  // xNodePos() = x0NodePos()+soln[displacement]
	  //  as a part of the set up for each iteration. 
	  //*****************
	  //*****************

	   xratio[iCell] *= 
	     ( nodeRight->x0NodePos()  - nodeCent->x0NodePos()) /
	     (( nodeRight->x0NodePos()+soln[indexRight_EqnStart + nodeTmpsRight.Offset_Displacement_Axial ]) - (nodeCent->xNodePos() +soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial ] ))  ;
	}
	else
	  xratio[iCell] = -9e9; // should never be used. 
      }

      // the divergence of the pressure == - the trace of the STRESS tensor No factor of 3 as we are 1 d

      // HOWEVER my understanding is that the pressure variable is the fluid pressure, not the solid matrix pressure. 
      // SO the below result may need to be modified by the non-solid-volume ratio. 
      // Also!!!
      // With the assumption that the time step is much much greater than the sound speed across the whole battery layer;
      // then this is quasi static: the """correct""" way to calculate the pressure induced strain would be to average the pressure
      // across the whole battery, and adjust the strain in each cell so that average would be equal to the right boundary pressure. 
      // While the implicit solution generated by calculating the residuals as below will give the same result, the convergance
      // of the solution may be greatly lenghtened by not explicitly calculating the quasi-static solution. 

      size_t iVar_Pressure = nodeCent->indexBulkDomainVar0((size_t) Pressure_Axial);
      double pressure_STRESS = 0;
      if (iVar_Pressure != npos) {
	if(nodeRight && ! nodeLeft)
	  pressure_STRESS = -(1.0/BulkMod)*(soln[indexRight_EqnStart + iVar_Pressure]-soln[indexCent_EqnStart + iVar_Pressure]);
	else if(nodeLeft && ! nodeRight)
	  pressure_STRESS = -(1.0/BulkMod)*(soln[indexLeft_EqnStart + iVar_Pressure]-soln[indexCent_EqnStart + iVar_Pressure]);
	else if (nodeLeft && nodeRight) {
	  pressure_STRESS = -(1.0/BulkMod)*(
					    (soln[indexLeft_EqnStart + iVar_Pressure]-soln[indexCent_EqnStart + iVar_Pressure])+
					    (soln[indexRight_EqnStart + iVar_Pressure]-soln[indexCent_EqnStart + iVar_Pressure])
					    );	    	    

	}
      }

      double pressure_strain = pressure_STRESS/Eyoung;
      if ( (Domain1D::FluidPr | Domain1D::All) & solidMechanicsProbType_) {
	xratio[iCell]*= (1.0+pressure_strain);
      }
      // // now do the Solid Stess calculation
      // nodeTmpsCenter.Offset_Solid_Stress_Axial = nodeCent->indexBulkDomainVar0((size_t) Solid_Stress_Axial);
      // if(nodeLeft) 
      // 	nodeTmpsLeft.Offset_Solid_Stress_Axial   = nodeLeft->indexBulkDomainVar0((size_t) Solid_Stress_Axial);
      // else
      // 	nodeTmpsLeft.Offset_Solid_Stress_Axial   = nodeTmpsCenter.Offset_Solid_Stress_Axial;
      // if(nodeRight)
      // 	nodeTmpsRight.Offset_Solid_Stress_Axial   = nodeRight->indexBulkDomainVar0((size_t) Solid_Stress_Axial);
      // else
      // 	nodeTmpsRight.Offset_Solid_Stress_Axial   = nodeTmpsCenter.Offset_Solid_Stress_Axial;
	
      // double left_matrix_stress = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Solid_Stress_Axial] ;
      // double center_matrix_stress = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Solid_Stress_Axial] ;
      // double right_matrix_stress = soln[indexRight_EqnStart + nodeTmpsRight.Offset_Solid_Stress_Axial] ;
      // double matrix_pressure_left = - (left_matrix_stress- center_matrix_stress);
      // double matrix_pressure_right= - (center_matrix_stress - right_matrix_stress);
      // double matrix_LP_center = matrix_pressure_left - matrix_pressure_right;
      // xratio[iCell] *= (1+ matrix_LP_center/Eyoung);
	
      // avg_delta_matrix_pressure += matrix_pressure_left;
    }
    //    avg_delta_matrix_pressure /= NumLcCells;
#endif	
  } // end of iCell loop

#ifdef MECH_MODEL
  if (solidMechanicsProbType_ > 0 ) {
      new_node_pos.resize(NumLcCells,0.0);
      // node that this is offset by one node. 
      // The left most node is the right most of the anode, and has already been calculated. 
      for (int iCell = 1; iCell < NumLcCells; iCell++) {
	  cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
	  NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
	  NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;
	  NodalVars* nodeCent   = cTmps.nvCent_;
	  NodalVars* nodeLeft   = cTmps.nvLeft_;
	  nodeTmpsCenter.Offset_Displacement_Axial   = nodeCent->indexBulkDomainVar0((size_t) Displacement_Axial);
	  if (nodeLeft) {
	      nodeTmpsLeft.Offset_Displacement_Axial   = nodeLeft->indexBulkDomainVar0((size_t) Displacement_Axial);
	  } else {
	      nodeTmpsLeft.Offset_Displacement_Axial = -1;
	  }
	  indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;
	  indexLeft_EqnStart = nodeTmpsLeft.index_EqnStart;

	  // NOTE new_node_pos[0] must be set, using the information from the Anode
	  // this is a bit tortured, but on the anode, the right most node residual is set to

	  // res[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial ] = 
	  //   soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial]
	  //   -(new_node_pos[iCell]  - nodeCent->x0NodePos());

	  // so to extract the calculated new position, new_node_pos, reverse the summation:=?
	  if (iCell == 1) 
	  {
	    new_node_pos[0] =  res[indexLeft_EqnStart + nodeTmpsLeft.Offset_Displacement_Axial ] + 
	      - soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Displacement_Axial] - nodeLeft->x0NodePos();
	    new_node_pos[0] *=-1;
	    if(pMECH_MODEL_extra) std::cout << " separator new_node_pos[0] = "<< new_node_pos[0]<<" vs old pos " <<nodeLeft->x0NodePos()<< std::endl;
	    
	  }
	  double delta_0 =  (nodeCent->x0NodePos() + soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial]) 
	    - (nodeLeft->x0NodePos() + soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Displacement_Axial]);
	  double new_delta = delta_0 *  xratio[iCell-1]; // 
	  new_node_pos[iCell] = new_node_pos[iCell-1] + new_delta;

	  // // stress
	  // double left_matrix_stress = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Solid_Stress_Axial] ;
	  // double center_matrix_stress = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Solid_Stress_Axial] ;
	  // double lc_pressure = -(left_matrix_stress-center_matrix_stress);
	  // res[indexCent_EqnStart + nodeTmpsCenter.Offset_Solid_Stress_Axial] = left_matrix_stress + (avg_delta_matrix_pressure-lc_pressure); 
      }
      for (int iCell = 1; iCell < NumLcCells; iCell++) {
	// res[ left ] has already been set by the Anode code; hence start at iCell==1
	  cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
	  NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
	  NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;
	  NodalVars* nodeLeft  = cTmps.nvLeft_;
	  NodalVars* nodeCent  = cTmps.nvCent_;
	  indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;
	  indexLeft_EqnStart = nodeTmpsLeft.index_EqnStart;
	  nodeTmpsCenter.Offset_Displacement_Axial   = nodeCent->indexBulkDomainVar0((size_t) Displacement_Axial);
	  if (nodeLeft) {
	      nodeTmpsLeft.Offset_Displacement_Axial   = nodeLeft->indexBulkDomainVar0((size_t) Displacement_Axial);
	  } else {
	      nodeTmpsLeft.Offset_Displacement_Axial = -1;
	  }

	  res[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial ] =
	    soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial ] -
	    (new_node_pos[iCell] - nodeCent->x0NodePos() );
	  
	  if (pMECH_MODEL_extra) {
	    if(iCell==1)  std::cout << " Sep::residEval sol,res,r-1.0,x-x0:: ";
	    std::cout << " ( "<<soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial ]<<" "
		      << res[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial ] <<" "
		      << xratio[iCell]-1.0<<" "
		      <<(new_node_pos[iCell]  - nodeCent->x0NodePos())<<" ) ";
	    if(iCell == NumLcCells-1) std::cout<<endl;
	  }
      }
  }
#endif
}
//==================================================================================================================================
void
porousLiIon_Separator_dom1D::residEval_PreCalc(const bool doTimeDependentResid,
					       const Epetra_Vector* soln_ptr,
					       const Epetra_Vector* solnDot_ptr,
					       const Epetra_Vector* solnOld_ptr,
					       const double t,
					       const double rdelta_t,
					       const Zuzax::ResidEval_Type residType,
					       const Zuzax::Solve_Type solveType)
{
    static int tmpsSetup = 0;
    if (!tmpsSetup) {
        residSetupTmps();
        tmpsSetup = 1;
    }
    // Special Section to determine where to get temperature and pressure
    cellTmps& cTmps          = cellTmpsVect_Cell_[0];
    NodalVars* nodeCent = cTmps.nvCent_;
    bool haveTemp = true;
    size_t iVar_Temperature = nodeCent->indexBulkDomainVar0((size_t) Temperature);
    if (iVar_Temperature == npos) {
        haveTemp = false;
    }

    residType_Curr_ = residType;
    const Epetra_Vector& soln = *soln_ptr;

    t_final_ = t;
    if (rdelta_t < 1.0E-200) {
        t_init_ = t;
    } else {
        // We want an infinitly small time step 
        if (solveType == Zuzax::Solve_Type::TimeDependentInitial) {
           t_init_ = t;
        } else {
           t_init_ = t - 1.0/rdelta_t;
        }
    }

    /*
     *  ------------------------------ LOOP OVER CELL -------------------------------------------------
     *  Loop over the number of Cells in this domain on this processor
     *  This loop is done from left to right.
     */
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell;

        cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
	valCellTmps& valTmps = valCellTmpsVect_Cell_[iCell];

	NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
        NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;
        NodeTmps& nodeTmpsRight  = cTmps.NodeTmpsRight_;

	NodalVars* nodeCent = cTmps.nvCent_;

	if (haveTemp) {
            AssertTrace( nodeTmpsLeft.Offset_Temperature != npos);
            valTmps.Temperature.left   = soln[nodeTmpsLeft.index_EqnStart   + nodeTmpsLeft.Offset_Temperature];
            valTmps.Temperature.center = soln[nodeTmpsCenter.index_EqnStart + nodeTmpsCenter.Offset_Temperature];
            valTmps.Temperature.right  = soln[nodeTmpsRight.index_EqnStart  + nodeTmpsRight.Offset_Temperature];
        } else {
            AssertTrace( nodeTmpsLeft.Offset_Temperature == npos);
            valTmps.Temperature.left =  TemperatureReference_;
            valTmps.Temperature.center =  TemperatureReference_;
            valTmps.Temperature.right =  TemperatureReference_;
        }

	/*
         *   Get the pointer to the NodalVars object for the center node
         */
        int indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;
        /*
         * Calculate the cell width
         */
	xdelCell_Cell_[iCell] = cTmps.xCellBoundaryR_ - cTmps.xCellBoundaryL_;

        for (int k = 0; k < nsp_; k++) {
            Xcent_cc_[k] = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_MoleFraction_Species + k];
        }
        Vcent_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Voltage];

        //
        //  Alter the porosity here for some problem types. Other problem types involve calculating an equation 
        //      -> we do this here when the porosity is not a formal variable.
        //
        if (porosityEquationProbType_  & Porosity_EqnType_Status::CalculatedOutOfEqnSystem) {
            if (porosityEquationProbType_  & Porosity_EqnType_Status::PartOfMechanics) {
               // do something different
            } else {
                porosity_Cell_[iCell] = calcPorosity(iCell);
            }
        }
      
        /*
         * Setup the thermo
         */
        SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
       
	//
	//  Get the total lyte, total solid and total cell heat capacity of the current cell
	//  Store them in vectors for later use.
	//
	if (energyEquationProbType_) {
	    CpMolar_total_Cell_[iCell] = getCellHeatCapacity(nodeCent, &(soln[indexCent_EqnStart]));
	    //
	    // Calculate the thermal conductivity of the porous matrix if we are calculating the energy equation
	    //
	    thermalCond_Cell_[iCell] = thermalCondCalc_PorMatrix();
	}
    } // end of icell loop
}
//==================================================================================================================================
void
porousLiIon_Separator_dom1D::eval_PostSoln(
    const bool doTimeDependentResid,
    const Epetra_Vector *soln_ptr,
    const Epetra_Vector *solnDot_ptr,
    const Epetra_Vector *solnOld_ptr,
    const double t,
    const double rdelta_t)
{
    NodalVars* nodeCent = 0;
    NodalVars* nodeLeft = 0;
    NodalVars* nodeRight = 0;
    int indexCent_EqnStart, indexLeft_EqnStart, indexRight_EqnStart;
    const Epetra_Vector& soln = *soln_ptr;
    double xdelL; // Distance from the center node to the left node
    double xdelR; // Distance from the center node to the right node

    double deltaT = 0.0;
    if (rdelta_t > 0.0) {
        deltaT = 1.0 / rdelta_t;
    }

    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell;

        cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
        NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
        NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;
        NodeTmps& nodeTmpsRight  = cTmps.NodeTmpsRight_;

	qSource_Cell_curr_[iCell] = 0.0;
        jouleHeat_lyte_Cell_curr_[iCell] = 0.0;
        jouleHeat_solid_Cell_curr_[iCell] = 0.0;

	/*
         *  ---------------- Get the index for the center node ---------------------------------
         *   Get the pointer to the NodalVars object for the center node
	 *   Index of the first equation in the bulk domain of center node
         */
	nodeCent = cTmps.nvCent_;
	indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;
        for (int k = 0; k < nsp_; k++) {
            Xcent_cc_[k] = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_MoleFraction_Species + k];
        }
        Vcent_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Voltage];
	if (iCell == 0) {
	    potentialAnodic_ = Vcent_cc_;
	}
	if (iCell == NumLcCells-1) {
	    potentialCathodic_ = Vcent_cc_;
	}

        /*
         *  ------------------- Get the index for the left node -----------------------------
         *    There may not be a left node if we are on the left boundary. In that case
         *    set the pointer to zero and the index to -1.
	 *    The solution index is set to the center solution index in that case as well.
         */
	nodeLeft = cTmps.nvLeft_;
	indexLeft_EqnStart = nodeTmpsLeft.index_EqnStart;    
	for (int k = 0; k < nsp_; k++) {
            Xleft_cc_[k] = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_MoleFraction_Species + k];
        }
        Vleft_cc_ = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Voltage];

        /*
         * If we are past the first cell, then we have already done the calculation
         * for this flux at the right cell edge of the previous cell
         */
    
	xdelL = cTmps.xdelL_;
        /*
         * Calculate the distance between the right and center node points
         */
        xdelR = cTmps.xdelR_;
        /*
         * Calculate the cell width
         */
      

        /*
         * ------------------------ Get the indexes for the right node ------------------------------------
         */
	nodeRight = cTmps.nvRight_;
	indexRight_EqnStart = nodeTmpsRight.index_EqnStart;
	for (int k = 0; k < nsp_; k++) {
            Xright_cc_[k] = soln[indexRight_EqnStart + nodeTmpsRight.Offset_MoleFraction_Species + k];
        }
        Vright_cc_ = soln[indexRight_EqnStart + nodeTmpsRight.Offset_Voltage];


	if (nodeLeft != 0) {
	    /*
	     *  Establish the environment at the left cell boundary
	     */
	    SetupThermoShop2(nodeLeft, &(soln[indexLeft_EqnStart]), nodeCent, &(soln[indexCent_EqnStart]), 0);
	    
	    SetupTranShop(xdelL, 0);
	    /*
	     * Calculate the flux at the left boundary for each equation
	     */
	    gradV_trCurr_ = (Vcent_cc_ - Vleft_cc_) / xdelL;

	    /*
	     * Calculate the flux of species and the flux of charge
	     *   - the flux of charge must agree with the flux of species
	     */
	    icurrElectrolyte_CBL_[iCell] = 0.0;
	    for (int k = 0; k < nsp_; k++) {
		icurrElectrolyte_CBL_[iCell] += jFlux_trCurr_[k] * spCharge_[k];
	    }
	    icurrElectrolyte_CBL_[iCell] *= (ZZCantera::Faraday);
            if (nodeRight == 0) {
                 icurrElectrolyte_CBR_[iCell] = icurrElectrolyte_CBL_[iCell];
            }

	    qSource_Cell_curr_[iCell]       += - gradV_trCurr_ * icurrElectrolyte_CBL_[iCell] * xdelL * 0.5 * deltaT;
	    jouleHeat_lyte_Cell_curr_[iCell]+= - gradV_trCurr_ * icurrElectrolyte_CBL_[iCell] * xdelL * 0.5 * deltaT;
	}

	if (nodeRight != 0) {
	    /*
             *  Establish the environment at the right cell boundary
             */
            SetupThermoShop2(nodeCent, &(soln[indexCent_EqnStart]), nodeRight, &(soln[indexRight_EqnStart]), 1);

            SetupTranShop(xdelR, 1);

            /*
             * Calculate the flux at the right boundary for each equation
             * This is equal to
             *       Conc * Vaxial * phi
             */
	    gradV_trCurr_ = (Vright_cc_ - Vcent_cc_) / xdelR;

            /*
             * Calculate the flux of species and the flux of charge
             *   - the flux of charge must agree with the flux of species
             */
            icurrElectrolyte_CBR_[iCell] = 0.0;
            for (int k = 0; k < nsp_; k++) {
                icurrElectrolyte_CBR_[iCell] += jFlux_trCurr_[k]* spCharge_[k];
            }
            icurrElectrolyte_CBR_[iCell] *= (ZZCantera::Faraday);
            if (nodeLeft == 0) {
                 icurrElectrolyte_CBL_[iCell] = icurrElectrolyte_CBR_[iCell];
            }

	    qSource_Cell_curr_[iCell]       += - gradV_trCurr_ * icurrElectrolyte_CBR_[iCell] * xdelR * 0.5 * deltaT;
	    jouleHeat_lyte_Cell_curr_[iCell]+= - gradV_trCurr_ * icurrElectrolyte_CBR_[iCell] * xdelR * 0.5 * deltaT;
  
	}
	qSource_Cell_accumul_[iCell] += qSource_Cell_curr_[iCell];
    }
}
//==================================================================================================================================
void
porousLiIon_Separator_dom1D::eval_HeatBalance(const int ifunc,
					      const double t,
					      const double deltaT,
					      const Epetra_Vector *soln_ptr,
					      const Epetra_Vector *solnDot_ptr,
					      const Epetra_Vector *solnOld_ptr,
					      class globalHeatBalVals& dVals)
{
    NodalVars* nodeCent = 0;
    NodalVars* nodeLeft = 0;
    NodalVars* nodeRight = 0;
    globalHeatBalValsBat* dValsB_ptr = dynamic_cast<globalHeatBalValsBat*>(&dVals);
    int indexCent_EqnStart, indexLeft_EqnStart, indexRight_EqnStart;
    const Epetra_Vector& soln = *soln_ptr;
    double xdelL; // Distance from the center node to the left node
    double xdelR; // Distance from the center node to the right node
 
#ifdef DEBUG_PRINT_CELL_TABLES
    int doPrint = 1;
#else
    int doPrint = 0;
#endif
    int doTimes = 2;
    int nColsTable = 173;

    //
    // Find the pointer for the left Anode
    SurfDomainDescription *leftS = BDD_ptr_->LeftSurf;
    BulkDomainDescription *leftDD = leftS->LeftBulk;
    BulkDomain1D *leftD_anode = leftDD->BulkDomainPtr_;
    //
    // Find the pointer for the left Cathode
    SurfDomainDescription *rightS = BDD_ptr_->RightSurf;
    BulkDomainDescription *rightDD = rightS->RightBulk;
    BulkDomain1D *rightD_cathode = rightDD->BulkDomainPtr_;

    //double phiIcurrL = 0.0;
    //double phiIcurrR = 0.0;
    double moleFluxLeft = 0.0;
    double moleFluxRight = 0.0;
    double residAdd = 0.0;

 
    for (int itimes = 0; itimes < doTimes; itimes++) {
	if (doPrint) {
	    drawline(1, nColsTable);
	    if (itimes == 0) {
		printf("           CELL BALANCES FOR ENTHALPY:    SEPARATOR\n");
		printf("Cell|   NewnEnth       OldnEnth        deltanEnth | ");
		printf("    fluxTLeft     FluxTRight | ");
		printf("    fluxL_JHPhi  fluxR_JHPhi | ");
		printf("  enthConvLeft enthConvRight | ");
		printf("    Resid_Add     Residual  |");
		printf("\n");
	    }
	    if (itimes == 1) {
		printf("\n\n                Analysys of Heat Source Terms:   SEPARATOR \n");
		printf("      ");
		printf("                JOULE HEATING  |   Conduction   |        REACTIONS\n");
		printf("Cell| ");
		printf("  lyteJHSource   solidJHSource |   HeatSource   | RxnHeatSrcTerm");
		printf("\n");
	    }
	}
	double fluxTright = 0.0, fluxTleft = 0.0, fluxR_JHPhi = 0.0, fluxL_JHPhi = 0.0, enthConvRight = 0.0,
	  enthConvLeft = 0.0, resid;
	double jouleHeat_lyte_total = 0.0;
	dVals.oldNEnthalpy = 0.0;
	dVals.newNEnthalpy = 0.0;
	dVals.totalHeatCapacity = 0.0;
	for (int iCell = 0; iCell < NumLcCells; iCell++) {
	    cIndex_cc_ = iCell;

	    cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
	    NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
	    NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;
	    NodeTmps& nodeTmpsRight  = cTmps.NodeTmpsRight_;

	    /*
	     *  ---------------- Get the index for the center node ---------------------------------
	     *   Get the pointer to the NodalVars object for the center node
	     *   Index of the first equation in the bulk domain of center node
	     */
	    nodeCent = cTmps.nvCent_;
	    indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;
	    for (int k = 0; k < nsp_; k++) {
		Xcent_cc_[k] = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_MoleFraction_Species + k];
	    }
	    Vcent_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Voltage];
	    if (iCell == 0) {
		potentialAnodic_ = Vcent_cc_;
	    }
	    if (iCell == NumLcCells-1) {
		potentialCathodic_ = Vcent_cc_;
	    }

	    /*
	     *  ------------------- Get the index for the left node -----------------------------
	     *    There may not be a left node if we are on the left boundary. In that case
	     *    set the pointer to zero and the index to -1.
	     *    The solution index is set to the center solution index in that case as well.
	     */
	    nodeLeft = cTmps.nvLeft_;
	    indexLeft_EqnStart = nodeTmpsLeft.index_EqnStart;    
	    for (int k = 0; k < nsp_; k++) {
		Xleft_cc_[k] = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_MoleFraction_Species + k];
	    }
	    Vleft_cc_ = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Voltage];

	    /*
	     * If we are past the first cell, then we have already done the calculation
	     * for this flux at the right cell edge of the previous cell
	     */
    
	    xdelL = cTmps.xdelL_;
	    /*
	     * Calculate the distance between the right and center node points
	     */
	    xdelR = cTmps.xdelR_;
	    /*
	     * Calculate the cell width
	     */
	    //double xdelCell = cTmps.xdelCell_;
	    /*
	     * ------------------------ Get the indexes for the right node ------------------------------------
	     */
	    nodeRight = cTmps.nvRight_;
	    indexRight_EqnStart = nodeTmpsRight.index_EqnStart;
	    for (int k = 0; k < nsp_; k++) {
		Xright_cc_[k] = soln[indexRight_EqnStart + nodeTmpsRight.Offset_MoleFraction_Species + k];
	    }
	    Vright_cc_ = soln[indexRight_EqnStart + nodeTmpsRight.Offset_Voltage];


	    if (nodeLeft != 0) {
		/*
		 *  Establish the environment at the left cell boundary
		 */
		SetupThermoShop2(nodeLeft, &(soln[indexLeft_EqnStart]), nodeCent, &(soln[indexCent_EqnStart]), 0);
	    
		SetupTranShop(xdelL, 0);
		/*
		 * Calculate the flux at the left boundary for each equation
		 */
		gradV_trCurr_ = (Vcent_cc_ - Vleft_cc_) / xdelL;
		fluxTleft = heatFlux_Curr_;
		fluxL_JHPhi = jFlux_EnthalpyPhi_Curr_;

		Fleft_cc_ =soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Velocity_Axial];
		moleFluxLeft = Fleft_cc_ * concTot_Curr_;
		enthConvLeft = moleFluxLeft * EnthalpyMolar_lyte_Curr_;

		/*
		 * Calculate the flux of species and the flux of charge
		 *   - the flux of charge must agree with the flux of species
		 */
		icurrElectrolyte_CBL_[iCell] = 0.0;
		for (int k = 0; k < nsp_; k++) {
		    icurrElectrolyte_CBL_[iCell] += jFlux_trCurr_[k] * spCharge_[k];
		}
		icurrElectrolyte_CBL_[iCell] *= (ZZCantera::Faraday);
		if (nodeRight == 0) {
		    icurrElectrolyte_CBR_[iCell] = icurrElectrolyte_CBL_[iCell];
		}

		//phiIcurrL = icurrElectrolyte_CBL_[iCell] * phiElectrolyte_Curr_;

	
	    } else {
		//
		//  Setup shop at the left node
		//
		SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));

		//phiIcurrL = icurrElectrolyte_CBL_[iCell] * phiElectrolyte_Curr_;

		Fleft_cc_ = DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity];
		moleFluxLeft = Fleft_cc_ * concTot_Curr_;

		fluxTleft = 0.0;
		SetupThermoShop1Extra(nodeCent, &(soln[indexCent_EqnStart]));

		for (int k = 0; k < nsp_; k++) {
		    jFlux_trCurr_[k] = DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k];
		}
		jFlux_EnthalpyPhi_Curr_ = 0.0;
		for (size_t k = 0; k < (size_t) nsp_; ++k) {
		    jFlux_EnthalpyPhi_Curr_ += jFlux_trCurr_[k] * EnthalpyPhiPM_lyte_Curr_[k];
		}
		//fluxL_JHPhi = jFlux_EnthalpyPhi_Curr_;
		fluxL_JHPhi = 0.0;

		//enthConvLeft = moleFluxLeft * EnthalpyMolar_lyte_Curr_;
		enthConvLeft = 0.0;
	    }

	    if (nodeRight != 0) {
		/*
		 *  Establish the environment at the right cell boundary
		 */
		SetupThermoShop2(nodeCent, &(soln[indexCent_EqnStart]), nodeRight, &(soln[indexRight_EqnStart]), 1);

		SetupTranShop(xdelR, 1);

		/*
		 * Calculate the flux at the right boundary for each equation
		 * This is equal to
		 *       Conc * Vaxial * phi
		 */
		gradV_trCurr_ = (Vright_cc_ - Vcent_cc_) / xdelR;
		fluxTright = heatFlux_Curr_;
		fluxR_JHPhi = jFlux_EnthalpyPhi_Curr_;
		Fright_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Velocity_Axial];
		moleFluxRight = Fright_cc_ * concTot_Curr_;
		enthConvRight = moleFluxRight * EnthalpyMolar_lyte_Curr_;

		/*
		 * Calculate the flux of species and the flux of charge
		 *   - the flux of charge must agree with the flux of species
		 */
		icurrElectrolyte_CBR_[iCell] = 0.0;
		for (int k = 0; k < nsp_; k++) {
		    icurrElectrolyte_CBR_[iCell] += jFlux_trCurr_[k] * spCharge_[k];
		}
		icurrElectrolyte_CBR_[iCell] *= (ZZCantera::Faraday);
		if (nodeLeft == 0) {
		    icurrElectrolyte_CBL_[iCell] = icurrElectrolyte_CBR_[iCell];
		}

		//phiIcurrR = icurrElectrolyte_CBR_[iCell] * phiElectrolyte_Curr_;
	    } else {

		//
		//  Setup shop at the right node
		//
		SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));

		//phiIcurrR = icurrElectrolyte_CBR_[iCell] * phiElectrolyte_Curr_;

		Fright_cc_ = DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity];
		moleFluxRight = Fright_cc_ * concTot_Curr_;

		SetupThermoShop1Extra(nodeCent, &(soln[indexCent_EqnStart]));
		//enthConvRight = moleFluxRight * EnthalpyMolar_lyte_Curr_;
		enthConvRight = 0.0;
		for (int k = 0; k < nsp_; k++) {
		    jFlux_trCurr_[k] = DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k];
		}
		jFlux_EnthalpyPhi_Curr_ = 0.0;
		for (size_t k = 0; k < (size_t) nsp_; ++k) {
		    jFlux_EnthalpyPhi_Curr_ += jFlux_trCurr_[k] * EnthalpyPhiPM_lyte_Curr_[k];
		}
		// fluxR_JHPhi = jFlux_EnthalpyPhi_Curr_;
		fluxR_JHPhi = 0.0;

		//fluxTright = DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Enthalpy_Conservation];
		fluxTright = 0.0;
	    }

	
	    if (doPrint) { 
		double deltanEnth = nEnthalpy_New_Cell_[iCell] - nEnthalpy_Old_Cell_[iCell];
		if (itimes == 0) {
		    resid = 0.0;
		   
		    resid = deltanEnth + deltaT *( fluxTright - fluxTleft);
		    resid += deltaT * (fluxR_JHPhi - fluxL_JHPhi + enthConvRight - enthConvLeft);

		    residAdd = 0.0;
		    if (iCell == 0) {
			residAdd = deltaT * 
			  leftD_anode->DomainResidVectorRightBound_LastResid_NE[nodeTmpsCenter.RO_Enthalpy_Conservation];
		    }
		    if (iCell == NumLcCells - 1) {
			residAdd = deltaT *
			  rightD_cathode->DomainResidVectorLeftBound_LastResid_NE[nodeTmpsCenter.RO_Enthalpy_Conservation];
		    }
		    resid += residAdd;

		    printf("%3d |  % 12.6E  % 12.6E  % 12.5E |", iCell, nEnthalpy_New_Cell_[iCell], 
			   nEnthalpy_Old_Cell_[iCell], deltanEnth);
		    printf("   % 12.5E  % 12.5E |",   deltaT *fluxTleft,  deltaT *fluxTright);
		    printf("   % 12.5E  % 12.5E |",   deltaT *fluxL_JHPhi,  deltaT *fluxR_JHPhi);
		    printf("   % 12.5E  % 12.5E |",   deltaT *enthConvLeft,  deltaT *enthConvRight);
		    printf("  % 12.5E  % 12.5E |", residAdd, resid);
		    printf("  \n");
		}
		if (itimes == 1) {
	
		    printf("%3d |  % 12.6E  % 12.6E  | ", iCell, jouleHeat_lyte_Cell_curr_[iCell], 0.0);
		    //double fluxL_JH = fluxL_JHPhi - phiIcurrL;
		    //double fluxR_JH = fluxR_JHPhi - phiIcurrR;
		    //double Einflux = deltaT * (enthConvLeft - enthConvRight + fluxL_JH -  fluxR_JH);
		    double rxnEnt = 0.0;
		    double HeatIn = deltaT * (fluxTleft - fluxTright);
		    printf(" % 12.6E | % 12.6E  ", 
			   HeatIn, rxnEnt);
		    printf("\n");
		}
	    }


	    dValsB_ptr->jouleHeat_lyte = jouleHeat_lyte_total;
	    dVals.totalHeatCapacity +=CpMolar_total_Cell_[iCell];
	    //
	    //  Count up the total old and new cell enthalpies
	    //
	    dVals.oldNEnthalpy += nEnthalpy_Old_Cell_[iCell];
	    dVals.newNEnthalpy += nEnthalpy_New_Cell_[iCell];
	}
        if (doPrint) {
	    drawline(1, nColsTable);
        }
    }
}
//==================================================================================================================================
// Evaluate species and element balances
void
porousLiIon_Separator_dom1D::eval_SpeciesElemBalance(const int ifunc,
						     const double t,
						     const double deltaT,
						     const Epetra_Vector *soln_ptr,
						     const Epetra_Vector *solnDot_ptr,
						     const Epetra_Vector *solnOld_ptr,
						     class globalHeatBalVals& dVals)
{
    NodalVars* nodeCent = 0;
    NodalVars* nodeLeft = 0;
    NodalVars* nodeRight = 0;
    globalHeatBalValsBat* dValsB_ptr = dynamic_cast<globalHeatBalValsBat*>(&dVals);
    int indexCent_EqnStart, indexLeft_EqnStart, indexRight_EqnStart;
    const Epetra_Vector& soln = *soln_ptr;
    double xdelL; // Distance from the center node to the left node
    double xdelR; // Distance from the center node to the right node
    double xdelCell;
#ifdef DEBUG_PRINT_CELL_TABLES
    int doPrint = 1;
#else
    int doPrint = 0;
#endif

    // Initially, we'll limit it to Li transport, PF6- bal, and solvent balance
    int doTimes = 1;
    int nColsTable = 175;
  
    //
    // Find the pointer for the left Anode
    SurfDomainDescription *leftS = BDD_ptr_->LeftSurf;
    BulkDomainDescription *leftDD = leftS->LeftBulk;
    BulkDomain1D *leftD_anode = leftDD->BulkDomainPtr_;
    //
    // Find the pointer for the left Cathode
    SurfDomainDescription *rightS = BDD_ptr_->RightSurf;
    BulkDomainDescription *rightDD = rightS->RightBulk;
    BulkDomain1D *rightD_cathode = rightDD->BulkDomainPtr_;

    double moleFluxLeft = 0.0;
    double moleFluxRight = 0.0;
    double residAdd = 0.0;

    ZZCantera::Array2D elem_Lyte_New_Cell(nsp_, NumLcCells, 0.0);
    ZZCantera::Array2D elem_Lyte_Old_Cell(nsp_, NumLcCells, 0.0);

    std::vector<double>& elemLi_Lyte_New = dValsB_ptr->elem_Lyte_New;
    std::vector<double>& elemLi_Lyte_Old = dValsB_ptr->elem_Lyte_Old;

    std::vector<double> icurrCBL_Cell(NumLcCells, 0.0);
    std::vector<double> icurrCBR_Cell(NumLcCells, 0.0);

    std::vector<double> convRight(nsp_, 0.0);
    std::vector<double> convLeft(nsp_, 0.0);
    std::vector<double> jFluxRight(nsp_, 0.0);
    std::vector<double> jFluxLeft(nsp_, 0.0);
    std::vector<double> species_Lyte_New_Curr(nsp_, 0.0);
    std::vector<double> species_Lyte_Old_Curr(nsp_, 0.0);

    std::vector<double>& species_Lyte_New_Total = dValsB_ptr->species_Lyte_New_Total;
    std::vector<double>& species_Lyte_Old_Total = dValsB_ptr->species_Lyte_Old_Total;   
    std::vector<double> species_convRight  = dValsB_ptr->species_convRight;
    std::vector<double> species_convLeft   = dValsB_ptr->species_convLeft;
    std::vector<double> species_jFluxRight = dValsB_ptr->species_jFluxRight;
    std::vector<double> species_jFluxLeft  = dValsB_ptr->species_jFluxLeft;

   
    ZZCantera::Array2D species_Lyte_New_Cell(nsp_, NumLcCells, 0.0);
    ZZCantera::Array2D species_Lyte_Old_Cell(nsp_, NumLcCells, 0.0);
    ZZCantera::Array2D res_Species(nsp_,NumLcCells, 0.0);
  
    for (int itimes = 0; itimes < doTimes; itimes++) {
	if (doPrint) {
	    if (itimes == 0) {
		drawline(1, nColsTable);
		printf(" \n\n            SEPARATOR::          Species Cell balances \n");
		printf("Cell Spec|     NewnBal       OldnBal     deltanSpecies | ");
		printf("  - SrcTerm   | ");
		printf("    fluxLeft       FluxRight | ");
		printf("    convLeft    convRight | ");
		printf("  Residual-additions  Residual  |  Current_CBL  Current CBR ");
		printf("\n");
		drawline(1, nColsTable);
	    }
	  
	}
	double resid;

	for (int iCell = 0; iCell < NumLcCells; iCell++) {
	    cIndex_cc_ = iCell;

	    cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
	    NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
	    NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;
	    NodeTmps& nodeTmpsRight  = cTmps.NodeTmpsRight_;
	    double icurrCBL = 0.0;
	    double icurrCBR = 0.0;

	    /*
	     *  ---------------- Get the index for the center node ---------------------------------
	     *   Get the pointer to the NodalVars object for the center node
	     *   Index of the first equation in the bulk domain of center node
	     */
	    nodeCent = cTmps.nvCent_;
	    indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;
	    for (int k = 0; k < nsp_; k++) {
		Xcent_cc_[k] = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_MoleFraction_Species + k];
	    } 
	    Vcent_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Voltage];
	  
	    /*
	     *  ------------------- Get the index for the left node -----------------------------
	     *    There may not be a left node if we are on the left boundary. In that case
	     *    set the pointer to zero and the index to -1.
	     *    The solution index is set to the center solution index in that case as well.
	     */
	    nodeLeft = cTmps.nvLeft_;
	    indexLeft_EqnStart = nodeTmpsLeft.index_EqnStart;    
	    for (int k = 0; k < nsp_; k++) {
		Xleft_cc_[k] = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_MoleFraction_Species + k];
	    }
	    Vleft_cc_ = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Voltage];


	    xdelL = cTmps.xdelL_;
	    xdelR = cTmps.xdelR_;
	    xdelCell = cTmps.xdelCell_;


	    double* mfElectrolyte_Soln_old = mfElectrolyte_Soln_Cell_old_.ptrColumn(iCell);

	    /*
	     * ------------------------ Get the indexes for the right node ------------------------------------
	     */
	    nodeRight = cTmps.nvRight_;
	    indexRight_EqnStart = nodeTmpsRight.index_EqnStart;
	    for (size_t k = 0; k < (size_t) nsp_; k++) {
		Xright_cc_[k] = soln[indexRight_EqnStart + nodeTmpsRight.Offset_MoleFraction_Species + k];
	    }
	    Vright_cc_ = soln[indexRight_EqnStart + nodeTmpsRight.Offset_Voltage];

	    if (nodeLeft != 0) {
		/*
		 *  Establish the environment at the left cell boundary
		 */
		SetupThermoShop2(nodeLeft, &(soln[indexLeft_EqnStart]), nodeCent, &(soln[indexCent_EqnStart]), 0);	    
		SetupTranShop(xdelL, 0);
	
		Fleft_cc_ =soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Velocity_Axial];
		moleFluxLeft = Fleft_cc_ * concTot_Curr_;
		for (size_t k = 0; k < (size_t) nsp_; k++) {
		    convLeft[k] =  moleFluxLeft *  mfElectrolyte_Soln_Curr_[k];
		    jFluxLeft[k] = jFlux_trCurr_[k];
		    icurrCBL += jFluxLeft[k] *spCharge_[k] * Faraday;
		}
	    } else {
		//
		//  Setup shop at the left node
		//
		SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
		SetupThermoShop1Extra(nodeCent, &(soln[indexCent_EqnStart]));

		Fleft_cc_ = DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity];
		moleFluxLeft = Fleft_cc_ * concTot_Curr_;
		for (size_t k = 0; k <  (size_t) nsp_; k++) {
		    convLeft[k] = 0.0;
		}
		for (size_t k = 0; k < (size_t) nsp_; k++) {
		    jFlux_trCurr_[k] = - DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k];
		    jFluxLeft[k] = 0.0;
		    icurrCBL += jFlux_trCurr_[k] *spCharge_[k] * Faraday;
		}
	    }

	    if (nodeRight != 0) {
		/*
		 *  Establish the environment at the right cell boundary
		 */
		SetupThermoShop2(nodeCent, &(soln[indexCent_EqnStart]), nodeRight, &(soln[indexRight_EqnStart]), 1);
		SetupTranShop(xdelR, 1);

		Fright_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Velocity_Axial];
		moleFluxRight = Fright_cc_ * concTot_Curr_;
		for (size_t k = 0; k <  (size_t) nsp_; k++) {
		    convRight[k] = moleFluxRight * mfElectrolyte_Soln_Curr_[k];
		    jFluxRight[k] = jFlux_trCurr_[k];
		    icurrCBR += jFluxRight[k] *spCharge_[k] * Faraday;
		}
	    } else {
		//
		//  Setup shop at the right node
		//
		SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));	
		SetupThermoShop1Extra(nodeCent, &(soln[indexCent_EqnStart]));

		Fright_cc_ = DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity];
		moleFluxRight = Fright_cc_ * concTot_Curr_;
		for (size_t k = 0; k <  (size_t) nsp_; k++) {
		    convRight[k] = 0.0;
		}
		for (size_t k = 0; k <  (size_t) nsp_; k++) {
		    jFlux_trCurr_[k] = DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k];
		    jFluxRight[k] = 0.0;
		    icurrCBR += jFlux_trCurr_[k] * spCharge_[k] * Faraday;
		}
	    }

	    //
	    // Return to the center node
	    //
	    SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
	    for (size_t k = 0; k < (size_t) nsp_; k++) {
		elem_Lyte_New_Cell(k, iCell) = porosity_Curr_ * concTot_Curr_ * xdelCell * mfElectrolyte_Soln_Curr_[k];
		elem_Lyte_Old_Cell(k, iCell) = (porosity_Cell_old_[iCell] * concTot_Cell_old_[iCell] *
						xdelCell  * mfElectrolyte_Soln_old[k]);
	    
		elemLi_Lyte_New[k] += elem_Lyte_New_Cell(k, iCell);
		elemLi_Lyte_Old[k] += elem_Lyte_Old_Cell(k, iCell);
	    }


	    for (size_t k = 0; k < (size_t) nsp_; k++) {
		species_Lyte_New_Curr[k] = porosity_Curr_ * concTot_Curr_* xdelCell  * mfElectrolyte_Soln_Curr_[k];
		species_Lyte_Old_Curr[k] = (porosity_Cell_old_[iCell] * concTot_Cell_old_[iCell]* xdelCell  * 
					    mfElectrolyte_Soln_old[k]);
		species_Lyte_New_Total[k] += species_Lyte_New_Curr[k];
		species_Lyte_Old_Total[k] += species_Lyte_Old_Curr[k];
		if (iCell == 0) {
		    species_convLeft[k] = convLeft[k];
		    species_jFluxLeft[k] =  jFluxLeft[k];	 
		}
		if (iCell == NumLcCells-1) {
		    species_convRight[k]  = convRight[k];
		    species_jFluxRight[k] = jFluxRight[k];
		}
	    }

	    if (doPrint) {
		if (itimes == 0) {
		    resid = 0.0;
		    for (size_t k = 0 ; k <  (size_t) nsp_; k++) {
			double deltanSp = species_Lyte_New_Curr[k] - species_Lyte_Old_Curr[k];
			double convSpRight = convRight[k];
			double convSpLeft = convLeft[k];
			double jFluxSpRight =  jFluxRight[k];
			double jFluxSpLeft =  jFluxLeft[k];
			resid = deltanSp;
			resid += deltaT * ( convSpRight - convSpLeft );
			resid += deltaT * ( jFluxSpRight - jFluxSpLeft );

			residAdd = 0.0;
			if (iCell == 0) {
			    residAdd = deltaT * 
			      leftD_anode->DomainResidVectorRightBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k];
			}
			if (iCell == NumLcCells - 1) {
			    residAdd = deltaT *
			      rightD_cathode->DomainResidVectorLeftBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k];
			}
			resid += residAdd;
			printf("%3d, %3d |  % 12.6E  % 12.6E  % 12.5E |", iCell, (int) k, 
			       species_Lyte_New_Total[k], species_Lyte_Old_Total[k], deltanSp);
			printf("  % 12.5E  % 12.5E |",  - deltaT * jFluxSpLeft,  deltaT * jFluxSpRight);
			printf("  % 12.5E  % 12.5E |",  - deltaT * convSpLeft,  deltaT * convSpRight);
			printf("  % 12.5E  % 12.5E |", residAdd,  resid);
			printf("  % 12.5E  % 12.5E |", icurrCBL, icurrCBR);
			printf("  \n");
		    }
		}
		drawline(1, nColsTable);
		printf("\n");
	    }
	}
    }
}
//==================================================================================================================================
void
porousLiIon_Separator_dom1D::SetupThermoShop1(const NodalVars* const nv, const double* const solnElectrolyte_Curr)
{
    updateElectrolyte(nv, solnElectrolyte_Curr);
  
    solidSkeleton_->setState_TP(temp_Curr_, pres_Curr_);

    porosity_Curr_ = porosity_Cell_[cIndex_cc_];
    
}
//==================================================================================================================================
void
porousLiIon_Separator_dom1D::SetupThermoShop1Extra(const NodalVars* const nv, 
						   const double* const solnElectrolyte_Curr)
{
    //
    // Calculate the EnthalpyPhi values at the CV interface and store these in  EnthalpyPhiPM_lyte_Curr_[]
    //
    ionicLiquid_->getPartialMolarEnthalpies(&EnthalpyPM_lyte_Curr_[0]);
    EnthalpyMolar_lyte_Curr_ = 0.0;
    for (size_t k = 0; k < BDT_ptr_->nSpeciesElectrolyte_; ++k) {
        double z = ionicLiquid_->charge(k);
        EnthalpyPhiPM_lyte_Curr_[k] = EnthalpyPM_lyte_Curr_[k] + Faraday * z * phiElectrolyte_Curr_;
	EnthalpyMolar_lyte_Curr_ += mfElectrolyte_Soln_Curr_[k] * EnthalpyPM_lyte_Curr_[k];
    }
}
//==================================================================================================================================
void
porousLiIon_Separator_dom1D::SetupThermoShop2(const NodalVars* const nvL, const double* const solnElectrolyte_CurrL,
                                              const NodalVars* const nvR, const double* const solnElectrolyte_CurrR,
                                              int type)
{
   // Needs major work
   // for (int i = 0; i < BDD_ptr_->NumEquationsPerNode; i++) {
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

    for (size_t k = 0; k < BDT_ptr_->nSpeciesElectrolyte_; ++k) {
        mfElectrolyte_Soln_Curr_[k] = 0.5 * (solnElectrolyte_CurrL[indexMFL+k] +solnElectrolyte_CurrR[indexMFR+k]);
    }

    calcMFElectrolyte_Thermo(&mfElectrolyte_Soln_Curr_[0], &mfElectrolyte_Thermo_Curr_[0]); 

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
    // Calculate the EnthalpyPhi values at the CV interface and store these in  EnthalpyPhiPM_lyte_Curr_[]
    ionicLiquid_->getPartialMolarEnthalpies(&EnthalpyPM_lyte_Curr_[0]);
    EnthalpyMolar_lyte_Curr_ = 0.0;
    for (size_t k = 0; k < BDT_ptr_->nSpeciesElectrolyte_; ++k) {
        double z = ionicLiquid_->charge(k);
        EnthalpyPhiPM_lyte_Curr_[k] = EnthalpyPM_lyte_Curr_[k] + Faraday * z * phiElectrolyte_Curr_;
	EnthalpyMolar_lyte_Curr_ += mfElectrolyte_Soln_Curr_[k] * EnthalpyPM_lyte_Curr_[k];
    }

    //
    // Calculate the matrix thermal conductivity from a series resistance on the two sides
    //
    double tmp;
    if (type == 0) {
        tmp = 1.0 / thermalCond_Cell_[cIndex_cc_ - 1] + 1.0 / thermalCond_Cell_[cIndex_cc_];
    } else {
        tmp = 1.0 / thermalCond_Cell_[cIndex_cc_ + 1] + 1.0 / thermalCond_Cell_[cIndex_cc_];
    }
    thermalCond_Curr_ = 2.0 / tmp;
}
//==================================================================================================================================
// Function updates the ThermoPhase object for the electrolyte given the solution vector
void
porousLiIon_Separator_dom1D::updateElectrolyte(const NodalVars* const nv, const double* const solnElectrolyte_Curr)
{
    /*
     *  Get the temperature: Check to see if the temperature is in the solution vector.
     *  If it is not, then use the reference temperature
     */
    temp_Curr_ = getPointTemperature(nv, solnElectrolyte_Curr);  
    /*
     *  Get the pressure
     */
    pres_Curr_ = PressureReference_;
    /*
     *  Assemble electrolyte mole fractions into mfElectrolyte_Thermo_Curr_[]
     */
    getMFElectrolyte_soln(nv, solnElectrolyte_Curr);
    /*
     *  Assemble electrolyte potential into phiElectrolyte_Curr_
     */
    getVoltages(nv, solnElectrolyte_Curr);
    /*
     *  Set the ThermoPhase states
     */
    ionicLiquid_->setState_TPX(temp_Curr_, pres_Curr_, &mfElectrolyte_Thermo_Curr_[0]);
    ionicLiquid_->setElectricPotential(phiElectrolyte_Curr_);
    //
    //  Calculate the total concentration of the electrolyte kmol m-3 and store into concTot_Curr_
    //
    concTot_Curr_ = ionicLiquid_->molarDensity();
}
//=====================================================================================================================
void
porousLiIon_Separator_dom1D::getVoltages(const NodalVars* const nv, const double* const solnElectrolyte_Curr)
{
    size_t indexVS = nv->indexBulkDomainVar0(Voltage);
    phiElectrolyte_Curr_ = solnElectrolyte_Curr[indexVS];
}
//=====================================================================================================================
/*
 *  We assume that setupthermoshop1 has been called.
 *  This calculates the Heat capacity per cross-sectional area (Joules/K m2)
 *
 */
double
porousLiIon_Separator_dom1D::getCellHeatCapacity(const NodalVars* const nv, const double* const solnElectrolyte_Curr)
{
    double CpMolar = ionicLiquid_->cp_mole();
    double lyteVol = porosity_Curr_ * xdelCell_Cell_[cIndex_cc_];
    CpMolar_lyte_Cell_[cIndex_cc_] = lyteVol * concTot_Curr_ * CpMolar;

    double solidVol = (1.0 - porosity_Curr_) * xdelCell_Cell_[cIndex_cc_];
    double cpMolarSolid = solidSkeleton_->cp_mole();
    double concSolid = solidSkeleton_->molarDensity();
    CpMolar_solid_Cell_[cIndex_cc_] = solidVol * concSolid * cpMolarSolid;
    double cptotal = CpMolar_lyte_Cell_[cIndex_cc_] + CpMolar_solid_Cell_[cIndex_cc_];
    return cptotal;
}
//==================================================================================================================================
void
porousLiIon_Separator_dom1D::SetupTranShop(const double xdel, const int type)
{
    /*
     * Determine diffusion velocities
     */

    //set gradients
    gradT_trCurr_ = 0.0;
    valCellTmps& valTmps = valCellTmpsVect_Cell_[cIndex_cc_];


    //
    // Calculate the Temperature derivative
    //
    if (type == 0) {
        // Left boundary
        gradV_trCurr_ = (Vcent_cc_ - Vleft_cc_) / xdel;
        
        gradT_trCurr_ = (valTmps.Temperature.center - valTmps.Temperature.left) / xdel;

    } else {
        // Right boundary
        gradV_trCurr_ = (Vright_cc_ - Vcent_cc_) / xdel;
        gradT_trCurr_ = (valTmps.Temperature.right - valTmps.Temperature.center) / xdel;
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
    Tortuosity tort(bruggemannExp);
    for (int k = 0; k < nsp_; k++) {
        Vdiff_trCurr_[k] *= tort.tortuosityFactor(porosity_Curr_);
    }
    //  Convert from diffusion velocity to diffusion flux
    for (int k = 0; k < nsp_; k++) {
        jFlux_trCurr_[k] = mfElectrolyte_Soln_Curr_[k] * concTot_Curr_ * Vdiff_trCurr_[k];
    }

    jFlux_EnthalpyPhi_Curr_ = 0.0;
    for (size_t k = 0; k < (size_t) nsp_; ++k) {
        jFlux_EnthalpyPhi_Curr_ += jFlux_trCurr_[k] * EnthalpyPhiPM_lyte_Curr_[k];
    }

    //
    // Calculate the heat flux
    //
    heatFlux_Curr_ = - thermalCond_Curr_ * gradT_trCurr_;
}

//==================================================================================================================================
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
porousLiIon_Separator_dom1D::saveDomain(ZZCantera::XML_Node& oNode,
                                    const Epetra_Vector* const soln_GLALL_ptr,
                                    const Epetra_Vector* const solnDot_GLALL_ptr,
                                    const double t,
                                    bool duplicateOnAllProcs)
{
    // get the NodeVars object pertaining to this global node
    GlobalIndices* gi = LI_ptr_->GI_ptr_;

    // Add an XML child for this domain. We will add the data to this child
    ZZCantera::XML_Node& bdom = oNode.addChild("domain");

    // Number of equations per node
    int numEquationsPerNode = BDD_ptr_->NumEquationsPerNode;

    // Vector containing the variable names as they appear in the solution vector
    std::vector<VarType>& variableNameList = BDD_ptr_->VariableNameList;

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
        NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
        varContig[i] = nv->x0NodePos();
    }
    ZZctml::addNamedFloatArray(gv, "X0", varContig.size(), &(varContig[0]), "m", "length");

    for (int iVar = 0; iVar < numEquationsPerNode; iVar++) {
        VarType vt = variableNameList[iVar];
        i = 0;
        std::string nmm = vt.VariableName(200);
        for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
            NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
            int istart = nv->EqnStart_GbEqnIndex;
	    size_t offset = nv->indexBulkDomainVar(vt.VariableType, vt.VariableSubType);
	    if (offset == npos) {
		throw m1d_Error("porousLiIon_Separator_dom1D::saveDomain()", "cant find a variable");
	    }
            varContig[i] = (*soln_GLALL_ptr)[istart + offset];
        }
        ZZctml::addNamedFloatArray(gv, nmm, varContig.size(), &(varContig[0]), "kmol/m3", "concentration");
    }

    if (PS_ptr->doHeatSourceTracking_) {
       std::string nmm = "qSource_Cell_curr_";
       ZZctml::addNamedFloatArray(gv, nmm, numNodes, &(qSource_Cell_curr_[0]), "Joule/s/m2", "");
    }

}
//
//==================================================================================================================================
// Method for writing the header for the bulk domain to a tecplot file.
void
porousLiIon_Separator_dom1D::writeSolutionTecplotHeader()
{
    int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = !mypid ; //only proc 0 should write

    if (doWrite) {
    
	// open tecplot file
	string sss = id();
	char filename[50];
	sprintf(filename,"%s%s",sss.c_str(),".dat");
	FILE* ofp = fopen(filename, "w");

	//write title and variable list
	fprintf(ofp, "TITLE = \"Solution on Domain %s\"\n",sss.c_str() );

	// Number of equations per node
	size_t numVar = BDD_ptr_->NumEquationsPerNode;
	// Vector containing the variable names as they appear in the solution vector
	std::vector<VarType> &variableNameList = BDD_ptr_->VariableNameList;
	//! First global node of this bulk domain

	fprintf(ofp, "VARIABLES = ");
	fprintf(ofp, "\"x [m]\"  \n" );

	for (size_t k = 0; k < numVar; k++) {
	    VarType &vt = variableNameList[k];
	    std::string name = vt.VariableName(15);
	    fprintf(ofp, "\"%s\" \t", name.c_str() );
	}

        // print porosity, specific surface area, thickness for each control volume
	fprintf(ofp, "\"Porosity []\" \t");
	fprintf(ofp, "\"Control volume thickness [m]\" \t");

        //
        // volume fractions of control volume
        //
        for (size_t i = 0; i < numExtraCondensedPhases_; i++) {
            ExtraPhase* ep = ExtraPhaseList_[i];
            fprintf(ofp, "\"volFrac %s []\" \t", ep->phaseName.c_str());
        }

	//print thermal source terms
	// check dimensions!!
	fprintf(ofp, "\"qHeat_accum [J/m3]\" \t");
	fprintf(ofp, "\"qHeat_step [W/m3]\" \t");
	fprintf(ofp, "\"Joule_Lyte [W/m3]\" \t");
	fprintf(ofp, "\"Joule_Solid [W/m3]\" \t");
    
	fprintf(ofp, "\n" );
	fclose(ofp);
    }
}
//==================================================================================================================================
// Method for writing the solution on the bulk domain to a tecplot file.
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
porousLiIon_Separator_dom1D::writeSolutionTecplot(const Epetra_Vector *soln_GlAll_ptr, const Epetra_Vector *solnDot_GlAll_ptr,
			           const double t)
{
    int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = !mypid ; //only proc 0 should write
    int rsize = 0;
    if (doWrite) {

 
	double deltaTime = t_final_ - t_init_;
    
	// get the NodeVars object pertaining to this global node
	GlobalIndices *gi = LI_ptr_->GI_ptr_;
	// Number of equations per node
	int numVar = BDD_ptr_->NumEquationsPerNode;
	//! First global node of this bulk domain
	int firstGbNode = BDD_ptr_->FirstGbNode;
	//! Last Global node of this bulk domain
	int lastGbNode = BDD_ptr_->LastGbNode;
	int numNodes = lastGbNode - firstGbNode + 1;
	std::vector<double> vars(numNodes);
    
	//open tecplot file
	FILE* ofp;
	string sss = id();
	char filename[20];
	sprintf(filename,"%s%s",sss.c_str(),".dat");
	ofp = fopen(filename, "a");
    
	/*
	 *  Write out the Heading for the solution at the current time. It's put in a ZONE structure
	 *  with T being the heading and SOLUTIONTIME being the value of the time
	 */
  
#define NEW_TECPLOT
#ifndef NEW_TECPLOT
    fprintf(ofp, "ZONE T = \"t = %11.6E [s]\", I = %d, SOLUTIONTIME = %12.6E\n", t, numNodes, t);
  
    // ----------------------------------------------------------------------------------------------------------------------
    for (size_t iGbNode = (size_t) firstGbNode; iGbNode <= (size_t) lastGbNode; iGbNode++) {
	NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
	//x-position
	fprintf(ofp, "%g \t", nv->xNodePos());
	for (size_t iVar = 0; iVar < (size_t) numVar; iVar++) {
	    //other variables
	    int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
	    int istart = nv->EqnStart_GbEqnIndex;
	    fprintf(ofp, "%g \t", (*soln_GlAll_ptr)[istart + ibulk + iVar]);	
	}
	fprintf(ofp, "\n");
	// print porosity, specific surface area, thickness for each control volume
	fprintf(ofp, "%g \t", porosity_Cell_[iCell]);
	fprintf(ofp, "%g \t", xdelCell_Cell_[iCell]);
	//
	// volume fractions of control volume
	//
	for (size_t i = 0; i < numExtnumExtraCondensedPhases_; i++) {
	    fprintf(ofp, "%g \t", volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + i]);
	}
	// print thermal source terms
	int iCell = iGbNode - firstGbNode;
	fprintf(ofp, "%g \t", qSource_Cell_accumul_[iCell] / xdelCell_Cell_[iCell] );
	if ( deltaTime > 1e-80 ) {
	    fprintf(ofp, "%g \t", qSource_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime );
	    fprintf(ofp, "%g \t", jouleHeat_lyte_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime );
	    fprintf(ofp, "%g \t", jouleHeat_solid_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime );
	} else {
	    fprintf(ofp, "0.0 \t 0.0 \t 0.0 \t " );
	}
	fprintf(ofp, "\n");
    }
#else
    fprintf(ofp, "ZONE T = \"SEP , %10.4E\", I = %d, SOLUTIONTIME = %12.6E\n", t, numNodes, t);
    fprintf(ofp, "ZONETYPE = ORDERED\n");
    fprintf(ofp, "DATAPACKING = BLOCK\n");
    fprintf(ofp, "STRANDID = 2\n");
    // ----------------------------------------------------------------------------------------------------------------------
    // Print the x field
    //
   
    for (size_t iGbNode = (size_t) firstGbNode, i = 0; iGbNode <= (size_t) lastGbNode; ++iGbNode, ++i) {
	NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
	vars[i] = nv->xNodePos();
    }
    fwriteTecplotVector(ofp, vars, 13, 10);
    //
    // Print all of the solution variables
    //x
    for (size_t iVar = 0; iVar < (size_t) numVar; iVar++) {
	const VarType& vt = BDD_ptr_->VariableNameList[iVar];
	for (size_t iGbNode = (size_t) firstGbNode, i = 0; iGbNode <= (size_t) lastGbNode; ++iGbNode, ++i) {
	    NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
	    size_t istart = nv->EqnStart_GbEqnIndex;
	    size_t ioffset = nv->indexBulkDomainVar(vt);
	    if (ioffset == npos) {
		throw m1d_Error("", "index prob");
	    }
	    vars[i] = (*soln_GlAll_ptr)[istart + ioffset];
	}
	fwriteTecplotVector(ofp, vars, 13, 10);
    }
    //
    // Porosity of control volume
    //
    for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	vars[iCell] = porosity_Cell_[iCell];
    }
    fwriteTecplotVector(ofp, vars, 13);
    //
    // Thickness of each control volume units = m
    //
    for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	vars[iCell] = xdelCell_Cell_[iCell];
    }
    fwriteTecplotVector(ofp, vars, 13);
    //
    // volume fractions of control volume
    //
    for (size_t i = 0; i < numExtraCondensedPhases_; i++) {
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + i];
	}
	fwriteTecplotVector(ofp, vars, 13);
    }

    //
    // Print the qSource_Cell_accumul_ field
    //
    for (size_t iGbNode = (size_t) firstGbNode, i = 0; iGbNode <= (size_t) lastGbNode; ++iGbNode, ++i) {
	int iCell = iGbNode - firstGbNode;
	double val = qSource_Cell_accumul_[iCell] / xdelCell_Cell_[iCell];
	if (deltaTime < 1e-80) {
	    val = 0.0;
	}
	vars[i] = val;
    }
    fwriteTecplotVector(ofp, vars, 13, 10);
 
    //
    // Print the qSource_Cell_curr_ field
    //
    for (size_t iGbNode = (size_t) firstGbNode, i = 0; iGbNode <= (size_t) lastGbNode; ++iGbNode, ++i) {
	int iCell = iGbNode - firstGbNode;
	double val = qSource_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime;
	if ( deltaTime < 1e-80 ) {
	    val = 0.0;
	}
	vars[i] = val;
    }
    fwriteTecplotVector(ofp, vars, 13, 10);

    //
    // Print the qSource_Cell_curr_ field jouleHeat_lyte_Cell_curr_ field
    //
    for (size_t iGbNode = (size_t) firstGbNode, i = 0; iGbNode <= (size_t) lastGbNode; ++iGbNode, ++i) {
	int iCell = iGbNode - firstGbNode;
	double val = jouleHeat_lyte_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime;
	if (deltaTime < 1e-80) {
	    val = 0.0;
	}
	vars[i] = val;
    }
    fwriteTecplotVector(ofp, vars, 13, 10);

    //
    // Print  the jouleHeat_solid_Cell_curr_ field
    //
    for (size_t iGbNode = (size_t) firstGbNode; iGbNode <= (size_t) lastGbNode; iGbNode++) {
	int iCell = iGbNode - firstGbNode;
	double val = jouleHeat_solid_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime;
	if ( deltaTime < 1e-80 ) {
	    val = 0.0;
	}
	fprintf(ofp, "%20.13E ", val);
	if (++rsize >= 10) {
	    fprintf(ofp, "\n");
	    rsize = 0;
        }
    }
    if (rsize != 0) {
	fprintf(ofp, "\n");
    }

#endif

    fclose(ofp);
  }
}
//==================================================================================================================================
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
porousLiIon_Separator_dom1D::showSolution(const Epetra_Vector* soln_GlAll_ptr,
                                          const Epetra_Vector* solnDot_GlAll_ptr,
                                          const Epetra_Vector* soln_ptr,
                                          const Epetra_Vector* solnDot_ptr,
                                          const Epetra_Vector* solnOld_ptr,
                                          const Epetra_Vector_Owned* residInternal_ptr,
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
    std::vector<VarType>& variableNameList = BDD_ptr_->VariableNameList;
    int iBlock;
    int iGbNode;
    int n;
    int nPoints = BDD_ptr_->LastGbNode - BDD_ptr_->FirstGbNode + 1;

    char buf[100];
    string indent = "";
    for (int i = 0; i < indentSpaces; i++) {
        indent += " ";
    }
    const char* ind = indent.c_str();
    double v;
    GlobalIndices* gi = LI_ptr_->GI_ptr_;
    // Number of points in each vector
    string sss = id();
    stream0 ss;
    double oldVaxial = 0.0;
    double newVaxial = 0.0;


    if (do0Write) {
        drawline(indentSpaces, 80);
        ss.print0("%s  Solution on Bulk Domain %12s : Number of variables = %d\n", ind, sss.c_str(), NumDomainEqns);
        ss.print0("%s                                         : Number of Nodes = %d\n", ind, nPoints);
        ss.print0("%s                                         : Beginning pos %g\n", ind, BDD_ptr_->Xpos_start);
        ss.print0("%s                                         : Ending    pos %g\n", ind, BDD_ptr_->Xpos_end);
#ifdef MECH_MODEL
	int firstGbNode = BDD_ptr_->FirstGbNode;
	int lastGbNode = BDD_ptr_->LastGbNode;
	if (solidMechanicsProbType_ > 0) {
	    for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++) {
		NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
		double diff = nv->xNodePos() -nv->x0NodePos();
		ss.print0("%s                                         : Node[ %d ] displacement %g\n", ind, iGbNode,diff);
	    }
	}
#endif 

    }
    if (do0Write) {
        for (iBlock = 0; iBlock < nn; iBlock++) {
            drawline(indentSpaces, 80);
            ss.print0("%s   iGbNode  z   ", ind);
            for (n = 0; n < 5; n++) {
                int ivar = iBlock * 5 + n;
                VarType vt = variableNameList[ivar];
                string name = vt.VariableName(15);
                ss.print0(" %15s", name.c_str());
            }
            ss.print0("\n");
            drawline(indentSpaces, 80);

            for (iGbNode = BDD_ptr_->FirstGbNode; iGbNode <= BDD_ptr_->LastGbNode; iGbNode++) {
                NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
                double x = nv->xNodePos();
                ss.print0("\n%s   %4d  % -10.4E ", ind, iGbNode, x);
                int istart = nv->EqnStart_GbEqnIndex;
                for (n = 0; n < 5; n++) {
                    int ivar = iBlock * 5 + n;
                    VarType vt = variableNameList[ivar];
		    size_t offset = nv->indexBulkDomainVar(vt.VariableType, vt.VariableSubType);
		    if (offset == npos) {
			throw m1d_Error("porousLiIon_Separator_dom1D::showSolution()", "cant find a variable");
		    }
                    v = (*soln_GlAll_ptr)[istart + offset];
		    /*
                     *  Special case the axial velocity because it's not located at the nodes.
                     */
                    if (vt.VariableType == Velocity_Axial) {
                        newVaxial = v;
                        if (iGbNode == BDD_ptr_->FirstGbNode) {
 	                    cellTmps& cTmps = cellTmpsVect_Cell_[0];
	                    NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
                            // we've applied a boundary condition here!
                            v = DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity];
                        } else if (iGbNode == BDD_ptr_->LastGbNode) {
			    cellTmps& cTmps = cellTmpsVect_Cell_[NumLcCells-1];
	                    NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
                            // we've applied a boundary condition here!
                            v = DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity];
			    newVaxial = v;
			} else {
                            v = 0.5 * (oldVaxial + newVaxial);
                        }
                        oldVaxial = newVaxial;
                    }


                    ss.print0(" % -11.4E ", v);
                }
            }
            ss.print0("\n");
        }

        int nrem = NumDomainEqns - 5 * nn;
        if (nrem > 0) {
            drawline(indentSpaces, 80);
            ss.print0("%s    iGbNode  z   ", ind);
            ZZCantera::writelog(buf);
            for (n = 0; n < nrem; n++) {
                int ivar = nn * 5 + n;
                VarType vt = variableNameList[ivar];
                string name = vt.VariableName(15);
                ss.print0(" %15s", name.c_str());
                ZZCantera::writelog(buf);
            }
            ss.print0("\n");
            drawline(indentSpaces, 80);

            for (iGbNode = BDD_ptr_->FirstGbNode; iGbNode <= BDD_ptr_->LastGbNode; iGbNode++) {
                NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
                double x = nv->xNodePos();
                ss.print0("%s   %4d  % -10.4E ", ind, iGbNode, x);
                int istart = nv->EqnStart_GbEqnIndex;
                for (n = 0; n < nrem; n++) {
                    int ivar = iBlock * 5 + n;
                    VarType vt = variableNameList[ivar];
		    size_t offset = nv->indexBulkDomainVar(vt.VariableType, vt.VariableSubType);
		    if (offset == npos) {
			throw m1d_Error("porousLiIon_Separator_dom1D::showSolution()", "cant find a variable");
		    }
                    v = (*soln_GlAll_ptr)[istart + offset];
		    if (vt.VariableType == Velocity_Axial) {
                        newVaxial = v;
                        if (iGbNode == BDD_ptr_->FirstGbNode) {
 	                    cellTmps& cTmps = cellTmpsVect_Cell_[0];
	                    NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
                            // we've applied a boundary condition here!
                            v = DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity];
                        } else if (iGbNode == BDD_ptr_->LastGbNode) {
			    cellTmps& cTmps = cellTmpsVect_Cell_[0];
	                    NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
                            // we've applied a boundary condition here!
                            v = DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity];
			    newVaxial = v;
			} else {
                            v = 0.5 * (oldVaxial + newVaxial);
                        }
                        oldVaxial = newVaxial;
                    }
                    ss.print0(" % -11.4E ", v);
                }
                ss.print0("\n");
            }
        }
	print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
	if (do0Write) {
	    // -----------------------------------------------------------------------------------------------------------------
	    // --             print porosity within the cell --
	    // -----------------------------------------------------------------------------------------------------------------
	    ss.print0("\n");
	    drawline0(indentSpaces, 80);
	    ss.print0("%s  iGbNode     z      Porosity  ", ind);
	    for (size_t jPhase = 0; jPhase < numExtraCondensedPhases_; ++jPhase) {
		ExtraPhase* ep = ExtraPhaseList_[jPhase];
		ss.print0("%-11.11s ", ep->phaseName.c_str());
	    }
	    ss.print0("\n");
	    drawline0(indentSpaces, 80);
	}
	print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
	for (iGbNode = BDD_ptr_->FirstGbNode; iGbNode <= BDD_ptr_->LastGbNode; iGbNode++) {
	    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
	    if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
		int iCell = iGbNode - BDD_ptr_->FirstGbNode;
		NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
		double x = nv->xNodePos();
		ss.print0("%s %4d  %-10.4E ", ind, iGbNode, x);
		ss.print0("% -11.4E ", porosity_Cell_[iCell]);
		for (size_t jPhase = 0; jPhase < numExtraCondensedPhases_; ++jPhase) {
		    ss.print0("% -11.4E ", volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + jPhase]);
		}
		ss.print0("\n");
	    }
	    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
	}

#ifdef MECH_MODEL
    // if (do0Write) {
    //   drawline0(indentSpaces, 80);
    //   ss.print0("%s    CellBound    SolidStressAxial ", ind);
    //   ss.print0("\n");
    //   drawline(indentSpaces, 80);
    //   const Epetra_Vector& soln = *soln_ptr;
     
    //   for (int iCell = 1; iCell < NumLcCells; iCell++) {
    //   	int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    //   	int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
    //   	cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
    //   	NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
    //   	indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;
    //   	double Solid_Stress_Axial = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Solid_Stress_Axial];
    //   	ss.print0("%s %d-%d   % -11.4E ",ind,iCell-1,iCell, Solid_Stress_Axial);
    //   	ss.print0("\n");
    //   }
      
    // }
#endif 
        drawline(indentSpaces, 80);
        // ----------------------------------------------------
        // --             PRINT FLUXES AT THE CELL BOUNDARIES --
        // ----------------------------------------------------
        ss.print0("%s    CellBound    z     IcurrElectrolyte ", ind);
        ss.print0("\n");
        drawline(indentSpaces, 80);
    }
    NodalVars* nvl;
    NodalVars* nvr;
    double x;
    int iCell;
    for (iGbNode = BDD_ptr_->FirstGbNode; iGbNode <= BDD_ptr_->LastGbNode; iGbNode++) {
        print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
        if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
            ss.print0("%s    ", ind);
            if (iGbNode == BDD_ptr_->FirstGbNode) {
                iCell = 0;
                nvr = gi->NodalVars_GbNode[BDD_ptr_->FirstGbNode];
                x = nvr->xNodePos();
                ss.print0("Lft-0     % -10.4E ", x);
                ss.print0("% -11.4E ", icurrElectrolyte_CBL_[iCell]);
            } else if (iGbNode == BDD_ptr_->LastGbNode) {
                iCell = BDD_ptr_->LastGbNode - BDD_ptr_->FirstGbNode;
                nvr = gi->NodalVars_GbNode[BDD_ptr_->LastGbNode];
                x = nvr->xNodePos();
                ss.print0("%3d-Rgt   % -10.4E ", iCell, x);
                ss.print0("% -11.4E ", icurrElectrolyte_CBR_[iCell]);
            } else {
                iCell = iGbNode - BDD_ptr_->FirstGbNode;
                nvl = gi->NodalVars_GbNode[iGbNode];
                nvr = gi->NodalVars_GbNode[iGbNode + 1];
                x = 0.5 * (nvl->xNodePos() + nvr->xNodePos());
                ss.print0("%3d-%-3d   % -10.4E ", iCell, iCell + 1, x);
                ss.print0("% -11.4E ", icurrElectrolyte_CBR_[iCell]);
            }
            ss.print0("\n");
        }
        print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
    }

    if (PS_ptr->doHeatSourceTracking_) {
        if (do0Write) {
            // ----------------------------------------------------
            // --         PRINT HEAT GENERATION TERMS for DOMAIN --
            // ----------------------------------------------------

            ss.print0("\n");
            drawline0(indentSpaces, 80);
            ss.print0("%s   iGbNode  z     qHeat_step    qHeat_accum Joule_Lyte Joule_Solid ", ind);
            ss.print0("\n");
            drawline0(indentSpaces, 80);
        }
        //NodalVars* nvl;
        //NodalVars* nvr;
        for (iGbNode = BDD_ptr_->FirstGbNode; iGbNode <= BDD_ptr_->LastGbNode; iGbNode++) {
            print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
            if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
                iCell = iGbNode - BDD_ptr_->FirstGbNode;
                NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
                x = nv->xNodePos();
                ss.print0("%s %4d  % -10.4E ", ind, iGbNode, x);
                //Electrode* ee = Electrode_Cell_[iCell];
                ss.print0("% -11.4E ", qSource_Cell_curr_[iCell]);
                ss.print0("% -11.4E ",  qSource_Cell_accumul_[iCell]);
                ss.print0("% -11.4E ",  jouleHeat_lyte_Cell_curr_[iCell]);
                ss.print0("% -11.4E ",  jouleHeat_solid_Cell_curr_[iCell]);
                ss.print0("\n");
            }
            print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
        }
    }
}
//==================================================================================================================================
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
porousLiIon_Separator_dom1D::initialConditions(const bool doTimeDependentResid,
                                               Epetra_Vector* soln_ptr,
                                               Epetra_Vector* solnDot,
                                               const double t,
                                               const double delta_t)
{
    Epetra_Vector& soln = *soln_ptr;


    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell;

        /*
         *  ---------------- Get the index for the center node ---------------------------------
         */
        int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
        /*
         *   Get the pointer to the NodalVars object for the center node
         */
        NodalVars* nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
	NodalVars* nodeLeft = nodeCent;
	NodalVars* nodeRight = nodeCent;
	if (iCell !=0) {
	    int index_LeftLcNode = Index_LeftLcNode_LCO[iCell];
	    nodeLeft =  LI_ptr_->NodalVars_LcNode[ index_LeftLcNode];
	}
	if ((int) iCell !=NumLcCells - 1) {
	    int index_RightLcNode = Index_RightLcNode_LCO[iCell];
	    nodeRight = LI_ptr_->NodalVars_LcNode[index_RightLcNode];
	}
	double cellThickness = 0.5*(nodeRight->xNodePos() - nodeLeft->xNodePos());

        /*
         *  Index of the first equation in the bulk domain of center node
         */
        int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
        /*
         * Offsets for the variable unknowns in the solution vector for the electrolyte domain
         */
	size_t iVAR_Vaxial  = nodeCent->indexBulkDomainVar0((size_t) Velocity_Axial);
        size_t iVar_Species = nodeCent->indexBulkDomainVar0((size_t) MoleFraction_Species);
	size_t iVar_Temperature = nodeCent->indexBulkDomainVar0((size_t) Temperature);
	size_t iVar_Pressure = nodeCent->indexBulkDomainVar0((size_t) Pressure_Axial);
	size_t iVar_Voltage = nodeCent->indexBulkDomainVar0((size_t) Voltage);
	//
	// Find the start of the solution at the current node
	//
	//const double *solnCentStart = &(soln[indexCent_EqnStart]);
	//
	// Set Vaxial to zero
	//
	AssertThrow( iVAR_Vaxial != npos, "iVar_Vaxial isn't found");
        soln[indexCent_EqnStart + iVAR_Vaxial] = 0.0;
	//
	// Set the temperature if it is part of the solution vector
	//	
	temp_Curr_ = PSinput.TemperatureReference_;
        if (iVar_Temperature != npos) {
            soln[indexCent_EqnStart + iVar_Temperature] = PSinput.TemperatureReference_;
        }
	Temp_Cell_old_[iCell] = temp_Curr_;

	//
	// Set the pressure if it is part of the solution vector
	//
	pres_Curr_ = PSinput.PressureReference_;
        if (iVar_Pressure != npos) {
            soln[indexCent_EqnStart + iVar_Pressure] = PSinput.PressureReference_;
        }
        /*
         * Get initial mole fractions from PSinput
         */
	for (size_t k = 0; k < BDT_ptr_->nSpeciesElectrolyte_; ++k) {
	    soln[indexCent_EqnStart + iVar_Species + k] = PSinput.electrolyteMoleFracs_[k];
	}
	//
	// Set some arbitrary values for the one voltage
	//
        soln[indexCent_EqnStart + iVar_Voltage] = -0.07;

        cIndex_cc_ = iCell;

	//
	// Porosity Set up
	// 
	double porosity = 1.0;
	for (size_t k = 0; k < ExtraPhaseList_.size(); ++k) {
	    ExtraPhase* ep = ExtraPhaseList_[k];
	    ThermoPhase* tp = ep->tp_ptr;
	    tp->setState_TP(temp_Curr_, pres_Curr_);
	    double mvp = tp->molarVolume();
	    volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + k] = ep->volFraction;
	    volumeFraction_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + k] = ep->volFraction;
	    moleNumber_Phases_Cell_[numExtraCondensedPhases_ * iCell + k] = ep->volFraction * cellThickness * mvp;
	    moleNumber_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + k] = ep->volFraction * cellThickness * mvp;
	    porosity -= ep->volFraction;
	}
	porosity_Cell_[iCell] = porosity;
	porosity_Cell_old_[iCell] = porosity;

        SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));

        concTot_Cell_old_[iCell] = concTot_Curr_;     

        double* mfElectrolyte_Soln_old = mfElectrolyte_Soln_Cell_old_.ptrColumn(iCell);
        for (size_t k = 0; k < (size_t) nsp_; ++k) {
            mfElectrolyte_Soln_old[k] = mfElectrolyte_Soln_Curr_[k];
        }

        if (energyEquationProbType_ == 3) { 
	    double volCellNew = xdelCell_Cell_[iCell];
	    double lyteMolarEnthalpyNew = ionicLiquid_->enthalpy_mole();
	    double volLyteNew = porosity_Curr_ * volCellNew;
	    double lyteEnthalpyNew =  lyteMolarEnthalpyNew * concTot_Curr_ * volLyteNew;
	    double nEnthalpyInertNew = 0.0;
            for (size_t jPhase = 0; jPhase < numExtraCondensedPhases_; ++jPhase) {
                ExtraPhase* ep = ExtraPhaseList_[jPhase];
                ThermoPhase* tp = ep->tp_ptr;
                tp->setState_TP(temp_Curr_, pres_Curr_);
                double mEnth = tp->enthalpy_mole();
                double mn = moleNumber_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + jPhase];
                nEnthalpyInertNew += mEnth * mn;
            }

	    double nEnthalpy_New  = lyteEnthalpyNew + nEnthalpyInertNew;	    
	    nEnthalpy_Old_Cell_[iCell] = nEnthalpy_New; 
	    nEnthalpy_New_Cell_[iCell] = nEnthalpy_New; 
        }

    }
}
//=====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void porousLiIon_Separator_dom1D::setAtolVector(double atolDefault, const Epetra_Vector_Ghosted& soln,
                                                Epetra_Vector_Ghosted& atolVector,
                                                const Epetra_Vector_Ghosted* const atolV)
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
        NodalVars* nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
        /*
         *  Index of the first equation in the bulk domain of center node
         */
        int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
        /*
         * Offsets for the variable unknowns in the solution vector for the electrolyte domain
         */
	size_t iVAR_Vaxial  = nodeCent->indexBulkDomainVar0((size_t) Velocity_Axial);
        size_t iVar_Species = nodeCent->indexBulkDomainVar0((size_t) MoleFraction_Species);
	size_t iVar_Voltage = nodeCent->indexBulkDomainVar0((size_t) Voltage);
        /*
         * Set the atol value for the axial velocity
         *   arithmetically scaled -> so this is a characteristic value
         */
        double vax = soln[indexCent_EqnStart + iVAR_Vaxial];
        atolVector[indexCent_EqnStart + iVAR_Vaxial] = std::max(1.0E-4, 1.0E-1 * vax);
        /*
         * Set atol values for the species mole fractions
         */
        double val = atolDefault * 1.0E-2;
        if (val < 1.0E-12) {
            val = 1.0E-12;
        }
        for (int k = 0; k < nsp_; k++) {
            atolVector[indexCent_EqnStart + iVar_Species + k] = val;
        }
        /*
         * Set the atol value for the electrolyte voltage
         *      arithmetically scaled.-> so this is a characteristic value
         *         1 kcal gmol-1 = 0.05 volts
         */
        atolVector[indexCent_EqnStart + iVar_Voltage] = 1.0E-6;
        /*
         * Set the tolerance on the temperature
         */
        size_t iVar_Temperature = nodeCent->indexBulkDomainVar0((size_t) Temperature);
        if (iVar_Temperature != npos) {
            atolVector[indexCent_EqnStart + iVar_Temperature] = 1.0E-7;
        }
    }
}
//==================================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void porousLiIon_Separator_dom1D::setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted& soln,
                                                        const Epetra_Vector_Ghosted& solnDot,
                                                        Epetra_Vector_Ghosted& atolVector,
                                                        const Epetra_Vector_Ghosted* const atolV)
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
        NodalVars* nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
        /*
         *  Index of the first equation in the bulk domain of center node
         */
        int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
        /*
         * Offsets for the variable unknowns in the solution vector for the electrolyte domain
         */
        int iVAR_Vaxial  = nodeCent->indexBulkDomainVar0((size_t) Velocity_Axial);
        int iVar_Species = nodeCent->indexBulkDomainVar0((size_t) MoleFraction_Species);
	int iVar_Voltage = nodeCent->indexBulkDomainVar0((size_t) Voltage);
        /*
         * Set the atol value for the axial velocity
         *   arithmetically scaled -> so this is a characteristic value
         */
        double vax = soln[indexCent_EqnStart + iVAR_Vaxial];
        atolVector[indexCent_EqnStart + iVAR_Vaxial] = std::max(1.0E-4, 1.0E-1 * vax);

        /*
         * Set atol values for the species mole fractions time derivatives
         */
        double val = atolDefault * 1.0E-1;
        if (val < 1.0E-6) {
            val = 1.0E-6;
        }
        for (int k = 0; k < nsp_; k++) {
            atolVector[indexCent_EqnStart + iVar_Species + k] = val;
        }

        /*
         * Set the atol value for the electrolyte voltage
         *      arithmetically scaled.-> so this is a characteristic value
         *         1 kcal gmol-1 = 0.05 volts
         */
        atolVector[indexCent_EqnStart + iVar_Voltage] = 1.0E-6;
        /*
         * Set the tolerance on the temperature
         */
        size_t iVar_Temperature = nodeCent->indexBulkDomainVar0((size_t) Temperature);
        if (iVar_Temperature != npos) {
            atolVector[indexCent_EqnStart + iVar_Temperature] = 1.0E-7;
        }
    }
}
//==================================================================================================================================
void
porousLiIon_Separator_dom1D::setAtolDeltaDamping(double atolDefault, double relcoeff,
                                                 const Epetra_Vector_Ghosted& soln,
                                                 Epetra_Vector_Ghosted& atolDeltaDamping,
                                                 const Epetra_Vector_Ghosted* const atolV)
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
        NodalVars* nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
        /*
         *  Index of the first equation in the bulk domain of center node
         */
        int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
        /*
         * Offsets for the variable unknowns in the solution vector for the electrolyte domain
         */
	size_t iVAR_Vaxial  = nodeCent->indexBulkDomainVar0((size_t) Velocity_Axial);
        size_t iVar_Species = nodeCent->indexBulkDomainVar0((size_t) MoleFraction_Species);
	size_t iVar_Voltage = nodeCent->indexBulkDomainVar0((size_t) Voltage);
        /*
         * Set the atol value for the axial velocity
         *   arithmetically scaled -> so this is a characteristic value
         */
        atolDeltaDamping[indexCent_EqnStart + iVAR_Vaxial] = 1.0E-4 * relcoeff;
        /*
         * Set atol values for the species mole fractions
         */
        double val = atolDefault * 1.0E3;
        if (val < 1.0E-4) {
            val = 1.0E-4;
        }
        for (int k = 0; k < nsp_; k++) {
            atolDeltaDamping[indexCent_EqnStart + iVar_Species + k] = val * relcoeff;
        }
        /*
         * Set the atol value for the electrolyte voltage
         *      arithmetically scaled.-> so this is a characteristic value
         *         1 kcal gmol-1 = 0.05 volts
         */
        atolDeltaDamping[indexCent_EqnStart + iVar_Voltage] = 0.05 * relcoeff;
        /*
         * Set the tolerance on the temperature
         */
        size_t iVar_Temperature = nodeCent->indexBulkDomainVar0((size_t) Temperature);
        if (iVar_Temperature != npos) {
            atolDeltaDamping[indexCent_EqnStart + iVar_Temperature] = 5.0 * relcoeff;
        }

	/*
         * Set the tolerance on the displacment
	 *   -> turning it off for the moment.
         */
        size_t iVar_disp = nodeCent->indexBulkDomainVar0((size_t) Displacement_Axial);
        if (iVar_disp != npos) {
            atolDeltaDamping[indexCent_EqnStart + iVar_disp] = 5.0 * relcoeff;
        }

	/*
         * Set the tolerance on the Velocity_Radial
         */
        size_t iVar_Vrad = nodeCent->indexBulkDomainVar0((size_t) Velocity_Radial);
        if (iVar_Vrad != npos) {
            atolDeltaDamping[indexCent_EqnStart + iVar_Vrad] = 5.0 * relcoeff;
        }



    }
}
//==================================================================================================================================
void
porousLiIon_Separator_dom1D::setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff,
                                                         const Epetra_Vector_Ghosted& soln,
                                                         const Epetra_Vector_Ghosted& solnDot,
                                                         Epetra_Vector_Ghosted& atolDeltaDamping,
                                                         const Epetra_Vector_Ghosted* const atolV)
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
        NodalVars* nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
        /*
         *  Index of the first equation in the bulk domain of center node
         */
        int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
        /*
         * Offsets for the variable unknowns in the solution vector for the electrolyte domain
         */
	int iVAR_Vaxial  = nodeCent->indexBulkDomainVar0((size_t) Velocity_Axial);
        int iVar_Species = nodeCent->indexBulkDomainVar0((size_t) MoleFraction_Species);
	int iVar_Voltage = nodeCent->indexBulkDomainVar0((size_t) Voltage);
        /*
         * Set the atol value for the axial velocity
         *   arithmetically scaled -> so this is a characteristic value
         */
        atolDeltaDamping[indexCent_EqnStart + iVAR_Vaxial] = 1.0E-4 * relcoeff;
        /*
         * Set atol values for the species mole fractions
         */
        double val = atolDefault * 1.0E5;
        if (val < 1.0E-1) {
            val = 1.0E-1;
        }
        for (int k = 0; k < nsp_; k++) {
            atolDeltaDamping[indexCent_EqnStart + iVar_Species + k] = val * relcoeff;
        }
        /*
         * Set the atol value for the electrolyte voltage
         *      arithmetically scaled.-> so this is a characteristic value
         *         1 kcal gmol-1 = 0.05 volts
         */
        atolDeltaDamping[indexCent_EqnStart + iVar_Voltage] = 0.05 * relcoeff;
        /*
         * Set the tolerance on the temperature
         */
        size_t iVar_Temperature = nodeCent->indexBulkDomainVar0((size_t) Temperature);
        if (iVar_Temperature != npos) {
            atolDeltaDamping[indexCent_EqnStart + iVar_Temperature] = 5.0;
        }
    }
}
//==================================================================================================================================
void
porousLiIon_Separator_dom1D::err(const char* msg)
{
    printf("porousLiIon_Separator_dom1D: function not implemented: %s\n", msg);
    exit(-1);
}
//==================================================================================================================================
/*
 * Method to check for precipitation of the salts.
 * Returns index of offending cation or -1 if no precipitation
 */
int
porousLiIon_Separator_dom1D::checkPrecipitation()
{

    //no precipitation
    return -1;

}
//==================================================================================================================================
//  WARNING -> will fail for multiprocessors, becauase the accumulation of information within eval_PostProc will fail.
double porousLiIon_Separator_dom1D::effResistanceLayer(double &potAnodic, double &potCathodic, 
                                                       double &voltOCV, double &current)
{
    static double resistance = 0.0;
    potAnodic = potentialAnodic_;
    potCathodic = potentialCathodic_;
    current = 0.5 * (icurrElectrolyte_CBL_[0] + icurrElectrolyte_CBR_[0]);
    voltOCV = 0.0;
    //
    //  Calculate the effective resistance of the separator layer by taking the potential drop across the domain
    //  and dividing it by the current.
    //
    resistance = 0.0;
    if (fabs(current) > 1.0E-200) {
	resistance = (potAnodic - potCathodic) / current;
    }
    return resistance;
}
//==================================================================================================================================
}
//==================================================================================================================================


