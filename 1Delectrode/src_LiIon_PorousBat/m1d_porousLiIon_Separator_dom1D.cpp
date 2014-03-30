/**
 * @file m1d_porousLiIon_Separator_dom1D.cpp
 */

/*
 *   $Id: m1d_porousLiIon_Separator_dom1D.cpp 506 2013-01-07 22:43:59Z hkmoffa $
 */

#include "m1d_porousLiIon_Separator_dom1D.h"
#include "m1d_BDT_porSeparator_LiIon.h"
#include "m1d_NodalVars.h"
#include "m1d_GlobalIndices.h"
#include "m1d_Comm.h"

#include "cantera/transport/Tortuosity.h"
#include "cantera/thermo/StoichSubstance.h"

#include "m1d_ProblemStatementCell.h"
extern m1d::ProblemStatementCell PSinput;

#include "stdio.h"
#include "stdlib.h"

using namespace std;
using namespace Cantera;

namespace m1d
{

//=====================================================================================================================
porousLiIon_Separator_dom1D::porousLiIon_Separator_dom1D(BulkDomainDescription& bdd) :
    porousFlow_dom1D(bdd),
    ionicLiquid_(0), trans_(0), nph_(0), nsp_(0), concTot_cent_(0.0), concTot_cent_old_(0.0),
    concTot_Cell_(0), concTot_Cell_old_(0), cIndex_cc_(0), Fleft_cc_(0.0),
    Fright_cc_(0.0), Vleft_cc_(0.0), Vcent_cc_(0.0), Vright_cc_(0.0), Xleft_cc_(0), Xcent_cc_(0), Xright_cc_(0),
    spCharge_(0), mfElectrolyte_Soln_Curr_(0), mfElectrolyte_Thermo_Curr_(0), temp_Curr_(TemperatureReference_),
    pres_Curr_(PressureReference_), phiElectrolyte_Curr_(0.0), concTot_Curr_(0.0), porosity_Curr_(0.0),
    gradT_trCurr_(0.0), gradV_trCurr_(0), gradX_trCurr_(0), Vdiff_trCurr_(0), jFlux_trCurr_(0),
    icurrElectrolyte_CBL_(0), icurrElectrolyte_CBR_(0),
    iECDMC_(-1),
    iLip_(-1),
    iPF6m_(-1),
    solnTemp(0), ivb_(VB_MOLEAVG)
{

    BDT_porSeparator_LiIon* fa = dynamic_cast<BDT_porSeparator_LiIon*>(&bdd);
    if (!fa) {
        throw m1d_Error("porousLiIon_Separator_dom1D", "confused");
    }
    ionicLiquid_ = fa->ionicLiquid_;
    trans_ = fa->trans_;
    nsp_ = 3;
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
//=====================================================================================================================
porousLiIon_Separator_dom1D::porousLiIon_Separator_dom1D(const porousLiIon_Separator_dom1D& r) :
    porousFlow_dom1D(r.BDD_),
    ionicLiquid_(0), trans_(0), nph_(0), nsp_(0), concTot_cent_(0.0), concTot_cent_old_(0.0),
    concTot_Cell_(0), concTot_Cell_old_(0), cIndex_cc_(0), Fleft_cc_(0.0),
    Fright_cc_(0.0), Vleft_cc_(0.0), Vcent_cc_(0.0), Vright_cc_(0.0), Xleft_cc_(0), Xcent_cc_(0), Xright_cc_(0),
    spCharge_(0), mfElectrolyte_Soln_Curr_(0), mfElectrolyte_Thermo_Curr_(0), temp_Curr_(TemperatureReference_),
    pres_Curr_(PressureReference_), phiElectrolyte_Curr_(0.0), concTot_Curr_(0.0), porosity_Curr_(0.0),
    gradT_trCurr_(0.0), gradV_trCurr_(0), gradX_trCurr_(0), Vdiff_trCurr_(0), jFlux_trCurr_(0),
    icurrElectrolyte_CBL_(0), icurrElectrolyte_CBR_(0),
    iECDMC_(-1),
    iLip_(-1),
    iPF6m_(-1),
    solnTemp(0),
    ivb_(VB_MOLEAVG)
{
    porousLiIon_Separator_dom1D::operator=(r);
}
//=====================================================================================================================
porousLiIon_Separator_dom1D::~porousLiIon_Separator_dom1D()
{
}
//=====================================================================================================================
porousLiIon_Separator_dom1D&
porousLiIon_Separator_dom1D::operator=(const porousLiIon_Separator_dom1D& r)
{
    if (this == &r) {
        return *this;
    }
    // Call the parent assignment operator
    porousFlow_dom1D::operator=(r);

    ionicLiquid_                          = r.ionicLiquid_;
    trans_                                = r.trans_;
    nph_                                  = r.nph_;
    nsp_                                  = r.nsp_;
    concTot_cent_                         = r.concTot_cent_;
    concTot_cent_old_                     = r.concTot_cent_old_;
    concTot_Cell_                         = r.concTot_Cell_;
    concTot_Cell_old_                     = r.concTot_Cell_old_;

    cIndex_cc_                            = r.cIndex_cc_;
    Fleft_cc_                             = r.Fleft_cc_;
    Fright_cc_                            = r.Fright_cc_;
    Vleft_cc_                             = r.Vleft_cc_;
    Vcent_cc_                             = r.Vcent_cc_;
    Vright_cc_                            = r.Vright_cc_;
    Xleft_cc_                             = r.Xleft_cc_;
    Xcent_cc_                             = r.Xcent_cc_;
    Xright_cc_                            = r.Xright_cc_;
    spCharge_                             = r.spCharge_;
    mfElectrolyte_Soln_Curr_              = r.mfElectrolyte_Soln_Curr_;
    mfElectrolyte_Thermo_Curr_            = r.mfElectrolyte_Thermo_Curr_;
    mfElectrolyte_Soln_Cell_old_          = r.mfElectrolyte_Soln_Cell_old_;
    temp_Curr_                            = r.temp_Curr_;
    pres_Curr_                            = r.pres_Curr_;
    phiElectrolyte_Curr_                  = r.phiElectrolyte_Curr_;
    concTot_Curr_                         = r.concTot_Curr_;
    porosity_Curr_                        = r.porosity_Curr_;
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
    ivb_                                  = r.ivb_;

    return *this;
}
//=====================================================================================================================
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
    ThermoPhase* tmpPhase = & (PSinput.PhaseList_)->thermo(iph);
    StoichSubstance* inert = dynamic_cast<Cantera::StoichSubstance*>(tmpPhase);

    //need to convert inputs from cgs to SI
    double volumeSeparator =
        PSinput.separatorArea_ * PSinput.separatorThickness_;
    double volumeInert =
        PSinput.separatorMass_ / inert->density() ;
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

    Xleft_cc_.resize(nsp_, 0.0);
    Xcent_cc_.resize(nsp_, 0.0);
    Xright_cc_.resize(nsp_, 0.0);

    spCharge_.resize(nsp_, 0.0);
    for (int k = 0; k < nsp_; k++) {
        spCharge_[k] = ionicLiquid_->charge(k);
    }


    mfElectrolyte_Soln_Curr_.resize(nsp_, 0.0);
    mfElectrolyte_Thermo_Curr_.resize(nsp_, 0.0);

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
}
//====================================================================================================================
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
        int indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode] +
                                    nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];

        SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart_BD]), 0);

        concTot_Cell_old_[iCell] = concTot_Curr_;
        porosity_Cell_old_[iCell] = porosity_Curr_;

        double* mfElectrolyte_Soln_old = mfElectrolyte_Soln_Cell_old_.ptrColumn(iCell);
        mfElectrolyte_Soln_old[0] = mfElectrolyte_Soln_Curr_[0];
        mfElectrolyte_Soln_old[1] = mfElectrolyte_Soln_Curr_[1];
        mfElectrolyte_Soln_old[2] = mfElectrolyte_Soln_Curr_[2];

    }
}
//=====================================================================================================================
void
porousLiIon_Separator_dom1D:: residSetupTmps()
{
    size_t index_CentLcNode;

    NodalVars *nodeCent = 0;
    NodalVars *nodeLeft = 0;
    NodalVars *nodeRight = 0;

    size_t  indexCent_EqnStart;
    size_t  indexLeft_EqnStart;
    size_t  indexRight_EqnStart;

    int  index_LeftLcNode;
    int  index_RightLcNode;

    cellTmpsVect_Cell_.resize(NumLcCells);
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell;

	cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
	NodeTmps& nodeTmpsCenter = cTmps.NodeCenter_;
	NodeTmps& nodeTmpsLeft   = cTmps.NodeLeft_;
	NodeTmps& nodeTmpsRight  = cTmps.NodeRight_;

	/*
         *  ---------------- Get the index for the center node ---------------------------------
         */
        index_CentLcNode = Index_DiagLcNode_LCO[iCell];

        /*
         *   Get the pointer to the NodalVars object for the center node
         */
        nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
	nodeTmpsCenter.nv = nodeCent;
	cTmps.nodeCent = nodeCent; 
        /*
         *  Index of the first equation in the bulk domain of center node
         */
        indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
	nodeTmpsCenter.index_EqnStart = indexCent_EqnStart;
	
	/*
	 * Offsets for the variable unknowns in the solution vector for the electrolyte domain
	 */
	nodeTmpsCenter.Offset_Voltage = nodeCent->indexBulkDomainVar0((size_t)Voltage);



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
            indexLeft_EqnStart = indexCent_EqnStart;
	    nodeTmpsLeft.index_EqnStart = indexLeft_EqnStart;

	    nodeTmpsLeft.Offset_Voltage = nodeTmpsCenter.Offset_Voltage;

        } else {
            // get the node structure for the left node
            nodeLeft = LI_ptr_->NodalVars_LcNode[index_LeftLcNode];
            //index of first equation in the electrolyte of the left node
            indexLeft_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_LeftLcNode];
	    nodeTmpsLeft.index_EqnStart = indexLeft_EqnStart;
	    /*
	     *
	     */
	 
	    nodeTmpsLeft.Offset_Voltage = nodeLeft->indexBulkDomainVar0((size_t)Voltage);

        }
	cTmps.nodeLeft = nodeLeft;
	

	/*
         * ------------------------ Get the indexes for the right node ------------------------------------
         */
        index_RightLcNode = Index_RightLcNode_LCO[iCell];
        if (index_RightLcNode < 0) {
            nodeRight = 0;
            /*
             *  If there is no right node, we assign the right solution index to the center solution index
             */
            indexRight_EqnStart = indexCent_EqnStart;
	    nodeTmpsRight.index_EqnStart = indexRight_EqnStart;

	    nodeTmpsRight.Offset_Voltage = nodeTmpsCenter.Offset_Voltage;

        } else {
            //NodalVars
            nodeRight = LI_ptr_->NodalVars_LcNode[index_RightLcNode];
            //index of first equation of right node
            indexRight_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_RightLcNode];
	    nodeTmpsRight.index_EqnStart = indexRight_EqnStart;

	    nodeTmpsRight.Offset_Voltage = nodeRight->indexBulkDomainVar0((size_t)Voltage);
        }
	cTmps.nodeRight = nodeRight; 
    }
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
porousLiIon_Separator_dom1D::residEval(Epetra_Vector& res,
                                       const bool doTimeDependentResid,
                                       const Epetra_Vector* soln_ptr,
                                       const Epetra_Vector* solnDot_ptr,
                                       const Epetra_Vector* solnOld_ptr,
                                       const double t,
                                       const double rdelta_t,
                                       const ResidEval_Type_Enum residType,
                                       const Solve_Type_Enum solveType)
{
    static int tmpsSetup = 0;
    if (!tmpsSetup) {
	residSetupTmps();
	tmpsSetup = 1;
    }
    residType_Curr_ = residType;
    int index_RightLcNode;
    int index_LeftLcNode;
    int index_CentLcNode;

    NodalVars* nodeLeft = 0;
    NodalVars* nodeCent = 0;
    NodalVars* nodeRight = 0;

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

    const Epetra_Vector& soln = *soln_ptr;
    const Epetra_Vector& solnOld = *solnOld_ptr;

    /*
     * Index of the first equation at the left node corresponding to the first bulk domain, which is the electrolyte
     */
    int indexLeft_EqnStart;
    /*
     * Index of the first equation at the center node corresponding to the first bulk domain, which is the electrolyte
     */
    int indexCent_EqnStart;
    /*
     * Index of the first equation at the right node corresponding to the first bulk domain, which is the electrolyte
     */
    int indexRight_EqnStart;

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
        DiffFluxLeftBound_LastResid_NE[EQ_TCont_offset_BD] = fluxL;
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

	cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
	NodeTmps& nodeTmpsCenter = cTmps.NodeCenter_;
	NodeTmps& nodeTmpsLeft   = cTmps.NodeLeft_;
	NodeTmps& nodeTmpsRight  = cTmps.NodeRight_;


#ifdef DEBUG_HKM_NOT
        if (counterResBaseCalcs_ > 125 && residType == Base_ResidEval) {
            if (iCell == NumLcCells - 1) {
                // printf("we are here\n");
            }
        }
#endif
        /*
         *  Placeholder for debugging
         */
        if (IOwnLeft && iCell == 0) {
            if (residType == Base_ShowSolution) {
                cIndex_cc_ = iCell;
            }
        }
        if (IOwnLeft && iCell == 9) {
            if (residType == Base_ShowSolution) {
                cIndex_cc_ = iCell;
            }
        }

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
        indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
	/*
	 * Offsets for the variable unknowns in the solution vector for the electrolyte domain
	 */
	
//	size_t iVar_Voltage_Cent = nodeCent->indexBulkDomainVar0((size_t)Voltage);
	//size_t iVar_Voltage_Left = iVar_Voltage_Cent;
//	size_t iVar_Voltage_Right = iVar_Voltage_Cent;

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
            indexLeft_EqnStart = indexCent_EqnStart;
        } else {
            // get the node structure for the left node
            nodeLeft = LI_ptr_->NodalVars_LcNode[index_LeftLcNode];
            //index of first equation in the electrolyte of the left node
            indexLeft_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_LeftLcNode];
	    /*
	     *
	     */
	    //iVar_Voltage_Left = nodeLeft->indexBulkDomainVar0((size_t)Voltage);
        }
	//nodeLeft = cellTmps.NodeLeft_;
	AssertTrace(cTmps.nodeLeft == nodeLeft);
	AssertTrace(indexLeft_EqnStart == (int) nodeTmpsLeft.index_EqnStart);
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
            indexRight_EqnStart = indexCent_EqnStart;
        } else {
            //NodalVars
            nodeRight = LI_ptr_->NodalVars_LcNode[index_RightLcNode];
            //index of first equation of right node
            indexRight_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_RightLcNode];
        }
	AssertTrace(cTmps.nodeRight == nodeRight);
	AssertTrace(indexRight_EqnStart == (int) nodeTmpsRight.index_EqnStart);

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
            Xcent_cc_[k] = soln[indexCent_EqnStart + iVar_Species_BD + k];
        }
        Vcent_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Voltage];

        if (nodeLeft != 0) {
            /*
             * Find the velocity located at the left cell boundary.
             * The left cell boundary velocity is stored at the previous (left)
             * cell index as per our conventions.
             */
            Fleft_cc_ = soln[indexLeft_EqnStart + iVAR_Vaxial_BD];
            for (int k = 0; k < nsp_; k++) {
                Xleft_cc_[k] = soln[indexLeft_EqnStart + iVar_Species_BD + k];
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
        Fright_cc_ = soln[indexCent_EqnStart + iVAR_Vaxial_BD];

        if (nodeRight != 0) {

            for (int k = 0; k < nsp_; k++) {
                Xright_cc_[k] = soln[indexRight_EqnStart + iVar_Species_BD + k];
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
                SetupThermoShop2(nodeLeft, &(soln[indexLeft_EqnStart]), nodeCent, &(soln[indexCent_EqnStart]), 0);

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
            SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]), 0);
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
            SetupThermoShop2(nodeCent, &(soln[indexCent_EqnStart]), nodeRight, &(soln[indexRight_EqnStart]), 1);

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

#ifdef DEBUG_RESID
        double residBefore = 0.0;
        if (IOwnLeft && iCell == 0) {
            if (residType == Base_ShowSolution) {
                residBefore =  res[indexCent_EqnStart_BD + EQ_TCont_offset_BD];
                double tmp3 = 0.0;
                double sum5 = residBefore + fluxFright - tmp3;
                printf("ResidSep Cell = %d, ResBefore = -fluxleft = %10.3e, FluxRight = %10.3e, Prod = %10.3e total = %10.3e \n", iCell,
                       residBefore, fluxFright, tmp3, sum5);

            }
        }
#endif
        /*
         *  Total continuity equation - fluxFright and fluxFleft represent the total mass
         *                              fluxes coming and going from the cell.
         *                    R =   d rho d t + del dot (rho V) = 0
         */
        res[indexCent_EqnStart + EQ_TCont_offset_BD] += (fluxFright - fluxFleft);

        /*
         * Species continuity Equation
         */
        for (int k = 0; k < nsp_; k++) {
            if ((k != iECDMC_) && (k != iPF6m_)) {
                res[indexCent_EqnStart + EQ_Species_offset_BD + k] += (fluxXright[k] - fluxXleft[k]);
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
        res[indexCent_EqnStart + EQ_Current_offset_BD] += (icurrElectrolyte_CBR_[iCell] - icurrElectrolyte_CBL_[iCell]);

        /*
         * --------------------------------------------------------------------------
         *  Add in the source terms at the current cell center
         * --------------------------------------------------------------------------
         */

        SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]), 0);


        /*
         *  Total continuity equation -
         */

        /*
         * Mole fraction summation equation
         */
        res[indexCent_EqnStart + EQ_MFSum_offset_BD] = 1.0;
        for (int k = 0; k < nsp_; k++) {
            res[indexCent_EqnStart + EQ_MFSum_offset_BD] -= Xcent_cc_[k];
        }

        /*
         * Electroneutrality equation
         */
        res[indexCent_EqnStart + EQ_ChargeBal_offset_BD] = 0.0;
        for (int k = 0; k < nsp_; k++) {
            res[indexCent_EqnStart + EQ_ChargeBal_offset_BD] += Xcent_cc_[k] * spCharge_[k];
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
            double newStuffSpecies0 = Xcent_cc_[iLip_] * newStuffTC;

            /*
             * Setup shop with the old time step
             */
            SetupThermoShop1(nodeCent, &(solnOld[indexCent_EqnStart]), 0);

            double oldStuffTC = concTot_Cell_old_[iCell] * porosity_Cell_old_[iCell] * xdelCell;
            oldStuffTC = concTot_Curr_ * porosity_Curr_ * xdelCell;
            double oldStuffSpecies0 = mfElectrolyte_Soln_Curr_[iLip_] * oldStuffTC;
            double tmp = (newStuffSpecies0 - oldStuffSpecies0) * rdelta_t;
            // double delta_concTot_CurrNew = concTot_CurrNew - concTot_Curr_;
            //  double delta_mf = mfNew - mfElectrolyte_Soln_Curr_[0];
#ifdef DEBUG_HKM_NOT
            if (residType == Base_ResidEval) {
                printf(" deltaT term = %10.3e BulkSum = %10.3e\n", tmp, tmp + (fluxXright[0] - fluxXleft[0]));
            }
#endif

            /*
             *   .................... Add these terms in the residual
             */
            /*
             *  Add in the time term for species iLip
             */

            for (int k = 0; k < nsp_; k++) {
                if (k != iECDMC_ && k != iPF6m_) {
                    res[indexCent_EqnStart + EQ_Species_offset_BD + k] += tmp;
                }
            }

            /*
             *   Add in the time term for the total continuity equation
             *         note: the current problem will have this term equally zero always.
             *               However, we put it in here for the next problem.
             */
            //  res[indexCent_EqnStart_BD + EQ_TCont_offset_BD] += (newStuffTC - oldStuffTC) * rdelta_t;
            /*
             *   .................... Go back to setting up shop at the current time
             */
            SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]), 0);
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
        DiffFluxRightBound_LastResid_NE[EQ_TCont_offset_BD] = fluxR;
    }

}
//=====================================================================================================================
void
porousLiIon_Separator_dom1D::SetupThermoShop1(const NodalVars* const nv, const doublereal* const solnElectrolyte_Curr, int type)
{
    updateElectrolyte(nv, solnElectrolyte_Curr);
    if (type == 0) {
        porosity_Curr_ = porosity_Cell_[cIndex_cc_];
    }
}
//=====================================================================================================================
void
porousLiIon_Separator_dom1D::SetupThermoShop2(const NodalVars* const nvL, const doublereal* const solnElectrolyte_CurrL,
                                              const NodalVars* const nvR, const doublereal* const solnElectrolyte_CurrR,
                                              int type)
{
   // Needs major work
    for (int i = 0; i < BDD_.NumEquationsPerNode; i++) {
        solnTemp[i] = 0.5 * (solnElectrolyte_CurrL[i] + solnElectrolyte_CurrR[i]);
    }
    if (type == 0) {
        porosity_Curr_ = 0.5 * (porosity_Cell_[cIndex_cc_ - 1] + porosity_Cell_[cIndex_cc_]);
    } else {
        porosity_Curr_ = 0.5 * (porosity_Cell_[cIndex_cc_ + 1] + porosity_Cell_[cIndex_cc_]);
    }
    updateElectrolyte(nvR, &solnTemp[0]);
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
porousLiIon_Separator_dom1D::updateElectrolyte(const NodalVars* const nv, const doublereal* const solnElectrolyte_Curr)
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

    getMFElectrolyte_soln(nv, solnElectrolyte_Curr);
    getVoltages(nv, solnElectrolyte_Curr);

    ionicLiquid_->setState_TPX(temp_Curr_, pres_Curr_, &mfElectrolyte_Thermo_Curr_[0]);

    ionicLiquid_->setElectricPotential(phiElectrolyte_Curr_);

    // Calculate the total concentration of the electrolyte kmol m-3.
    concTot_Curr_ = ionicLiquid_->molarDensity();

}
//=====================================================================================================================
void
porousLiIon_Separator_dom1D::getVoltages(const NodalVars* const nv, const double* const solnElectrolyte_Curr)
{
    int indexVS = BDD_.VariableIndexStart_VarName[Voltage];
    int newIndexVS = nv->indexBulkDomainVar(Voltage, 0);
    AssertTrace(indexVS == newIndexVS);

    phiElectrolyte_Curr_ = solnElectrolyte_Curr[indexVS];
}
//=====================================================================================================================
void
porousLiIon_Separator_dom1D::getMFElectrolyte_soln(const NodalVars* const nv, const double* const solnElectrolyte_Curr)
{
    int indexMF = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
    // We are assuming that mf species start with subindex 0 here
    int newIndexMF = nv->indexBulkDomainVar(MoleFraction_Species, 0);
    AssertTrace(indexMF == newIndexMF);

    mfElectrolyte_Soln_Curr_[0] = solnElectrolyte_Curr[indexMF];
    mfElectrolyte_Soln_Curr_[1] = solnElectrolyte_Curr[indexMF + 1];
    mfElectrolyte_Soln_Curr_[2] = solnElectrolyte_Curr[indexMF + 2];
    double mf0 = MAX(mfElectrolyte_Soln_Curr_[0], 0.0);
    double mf1b = MAX(mfElectrolyte_Soln_Curr_[1], 0.0);
    double mf2b = MAX(mfElectrolyte_Soln_Curr_[2], 0.0);
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
}
//=====================================================================================================================
void
porousLiIon_Separator_dom1D::SetupTranShop(const double xdel, const int type)
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
    Tortuosity tort(bruggemannExp);
    for (int k = 0; k < nsp_; k++) {
        Vdiff_trCurr_[k] *= tort.tortuosityFactor(porosity_Curr_);
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
    std::vector<VarType>& variableNameList = BDD_.VariableNameList;
    int iBlock;
    int iGbNode;
    int n;
    int nPoints = BDD_.LastGbNode - BDD_.FirstGbNode + 1;

    char buf[100];
    string indent = "";
    for (int i = 0; i < indentSpaces; i++) {
        indent += " ";
    }
    const char* ind = indent.c_str();
    doublereal v;
    GlobalIndices* gi = LI_ptr_->GI_ptr_;
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
                NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
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
                NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
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
    NodalVars* nvl;
    NodalVars* nvr;
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
        int igECDMC  = PSinput.PhaseList_->globalSpeciesIndex("ECDMC");
        if (igECDMC < 0) {
            throw CanteraError("confused", "confused");
        }
        int igLip = PSinput.PhaseList_->globalSpeciesIndex("Li+");
        if (igLip < 0) {
            throw CanteraError("confused", "confused");
        }
        int igPF6m = PSinput.PhaseList_->globalSpeciesIndex("PF6-");
        if (igPF6m < 0) {
            throw CanteraError("confused", "confused");
        }

        soln[indexCent_EqnStart_BD + iVar_Species_BD + iECDMC_] = PSinput.electrolyteMoleFracs_[igECDMC];
        soln[indexCent_EqnStart_BD + iVar_Species_BD + iLip_] = PSinput.electrolyteMoleFracs_[igLip];
        soln[indexCent_EqnStart_BD + iVar_Species_BD + iPF6m_] = PSinput.electrolyteMoleFracs_[igPF6m];


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
    }
}
//=====================================================================================================================
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
        atolVector[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = MAX(1.0E-4, 1.0E-1 * vax);

        /*
         * Set atol values for the species mole fractions time derivatives
         */
        double val = atolDefault * 1.0E-1;
        if (val < 1.0E-6) {
            val = 1.0E-6;
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
//======================================================================================================================
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
        //double vax = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];
        atolDeltaDamping[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = 1.0E-4 * relcoeff;

        /*
         * Set atol values for the species mole fractions
         */
        double val = atolDefault * 1.0E3;
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
    }
}
//======================================================================================================================
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
        //double vax = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];
        atolDeltaDamping[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = 1.0E-4 * relcoeff;

        /*
         * Set atol values for the species mole fractions
         */
        double val = atolDefault * 1.0E5;
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
    }
}
//=====================================================================================================================
void
porousLiIon_Separator_dom1D::err(const char* msg)
{
    printf("porousLiIon_Separator_dom1D: function not implemented: %s\n", msg);
    exit(-1);
}
//=====================================================================================================================
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

//=====================================================================================================================
}
//=====================================================================================================================


