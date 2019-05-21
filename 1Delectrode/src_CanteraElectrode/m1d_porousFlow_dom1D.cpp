/**
 * @file m1d_porousFlow_dom1D.cpp  Definitions for the porous flow residual calculator
 */
/*
 * Copywrite 2013 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_porousFlow_dom1D.h"

#include "m1d_NodalVars.h"
#include "m1d_cellTmps_PorousFlow.h"

#include "m1d_CanteraElectrodeGlobals.h"
#include "m1d_ProblemStatementCell.h"
#include "m1d_exception.h"


using namespace std;
using namespace Zuzax;

//==================================================================================================================================
//! Global Problem input structure
/*!
 *   This contains the input data for the problem.
 *   We've made it a global structure, as there is one and only one instance of the structure
 */
extern m1d::ProblemStatementCell PSinput;

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{

//==================================================================================================================================
porousFlow_dom1D::porousFlow_dom1D(BDD_porousFlow* bdd_pf_ptr) :
    BulkDomain1D(bdd_pf_ptr),
    BDD_ptr_(bdd_pf_ptr),
    porosity_Cell_(0),
    porosity_Cell_old_(0),
    Temp_Cell_old_(0),
    xdelCell_Cell_(0),
    numExtraCondensedPhases_(0),
    volumeFraction_Phases_Cell_(0),
    volumeFraction_Phases_Cell_old_(0),
    moleNumber_Phases_Cell_(0),
    moleNumber_Phases_Cell_old_(0),
    cIndex_cc_(-1),
    temp_Curr_(TemperatureReference_),
    //    mm_stress_Curr_(SolidStressAxialRef_),
    pres_Curr_(PressureReference_),
    concTot_Curr_(0.0),
    phiElectrolyte_Curr_(0.0),
    porosity_Curr_(0.0),
    thermalCond_Curr_(0.0),
    heatFlux_Curr_(0.0), 
    jFlux_EnthalpyPhi_Curr_(0.0),
    EnthalpyMolar_lyte_Curr_(0.0),
    ivb_(VB_MOLEAVG),
    ionicLiquid_(0),
    trans_(0),
    solidSkeleton_(0),
    Porosity_prob_type_(0)
{
    //BDD_ptr_ = static_cast<BDD_porousFlow*>(&BDD_);
    ionicLiquid_ = BDD_ptr_->ionicLiquid_;
    trans_ = BDD_ptr_->trans_;
    solidSkeleton_ = BDD_ptr_->solidSkeleton_;

    energyEquationProbType_ = PSinput.Energy_equation_prob_type_;
    solidMechanicsProbType_ = PSinput.Solid_Mechanics_prob_type_;

    size_t sz = BDD_ptr_->ExtraPhaseList_.size();
    ExtraPhaseList_.resize(sz);
    for (size_t i = 0; i < ExtraPhaseList_.size(); ++i) {
        ExtraPhaseList_[i] = new ExtraPhase(*(BDD_ptr_->ExtraPhaseList_[i]));
    }
    Porosity_prob_type_        = BDD_ptr_->Porosity_prob_type_;
    porosityEquationProbType_  = BDD_ptr_->porosityEquationProbType_;
    crossSectionalArea_ = PSinput.crossSectionalArea_;
}
//=====================================================================================================================
porousFlow_dom1D::porousFlow_dom1D(const porousFlow_dom1D &r) :
    BulkDomain1D(r.BDD_ptr_),
    BDD_ptr_(r.BDD_ptr_),
    porosity_Cell_(0),
    porosity_Cell_old_(0),
    Temp_Cell_old_(0),
    xdelCell_Cell_(0),
    numExtraCondensedPhases_(0),
    volumeFraction_Phases_Cell_(0),
    volumeFraction_Phases_Cell_old_(0),
    moleNumber_Phases_Cell_(0),
    moleNumber_Phases_Cell_old_(0),
    cIndex_cc_(-1),
    temp_Curr_(TemperatureReference_),
    //    mm_stress_Curr_(SolidStressAxialRef_),
    pres_Curr_(PressureReference_),
    concTot_Curr_(0.0),
    phiElectrolyte_Curr_(0.0),
    porosity_Curr_(0.0),
    thermalCond_Curr_(0.0),
    heatFlux_Curr_(0.0),
    jFlux_EnthalpyPhi_Curr_(0.0),
    EnthalpyMolar_lyte_Curr_(0.0),
    ivb_(VB_MOLEAVG),
    ionicLiquid_(0),
    trans_(0),
    solidSkeleton_(0),
    ExtraPhaseList_(0),
    Porosity_prob_type_(0)
{
    //BDD_ptr_ = static_cast<BDD_porousFlow*>(&BDD_);
    porousFlow_dom1D::operator=(r);
}
//=====================================================================================================================
porousFlow_dom1D::~porousFlow_dom1D()
{
}
//=====================================================================================================================
porousFlow_dom1D& porousFlow_dom1D::operator=(const porousFlow_dom1D &r)
{
    if (this == &r) {
      return *this;
    }
    // Call the parent assignment operator
    BulkDomain1D::operator=(r);

    BDD_ptr_                  = r.BDD_ptr_;
    CpMolar_lyte_Cell_        = r.CpMolar_lyte_Cell_;
    CpMolar_solid_Cell_       = r.CpMolar_solid_Cell_;
    CpMolar_total_Cell_       = r.CpMolar_total_Cell_;
    EnthalpyPM_lyte_Cell_     = r.EnthalpyPM_lyte_Cell_;
    porosity_Cell_            = r.porosity_Cell_;
    porosity_Cell_old_        = r.porosity_Cell_old_;
    Temp_Cell_old_            = r.Temp_Cell_old_;
    xdelCell_Cell_            = r.xdelCell_Cell_;
    numExtraCondensedPhases_  = r.numExtraCondensedPhases_;
    volumeFraction_Phases_Cell_=r.volumeFraction_Phases_Cell_;
    volumeFraction_Phases_Cell_old_=r.volumeFraction_Phases_Cell_old_;
    moleNumber_Phases_Cell_   =r.moleNumber_Phases_Cell_;
    moleNumber_Phases_Cell_old_=r.moleNumber_Phases_Cell_old_;
    cIndex_cc_                = r.cIndex_cc_;
    temp_Curr_                = r.temp_Curr_;
    //    mm_stress_Curr_           = r.mm_stress_Curr_;
    pres_Curr_                = r.pres_Curr_;
    concTot_Curr_             = r.concTot_Curr_;
    phiElectrolyte_Curr_      = r.phiElectrolyte_Curr_;
    porosity_Curr_            = r.porosity_Curr_;
    thermalCond_Curr_         = r.thermalCond_Curr_;
    heatFlux_Curr_            = r.heatFlux_Curr_;
    jFlux_EnthalpyPhi_Curr_   = r.jFlux_EnthalpyPhi_Curr_;

    cellTmpsVect_Cell_        = r.cellTmpsVect_Cell_;
    mfElectrolyte_Soln_Curr_   = r.mfElectrolyte_Soln_Curr_;
    mfElectrolyte_Thermo_Curr_ = r.mfElectrolyte_Thermo_Curr_;
    EnthalpyPM_lyte_Curr_      = r.EnthalpyPM_lyte_Curr_;
    EnthalpyPhiPM_lyte_Curr_   = r.EnthalpyPhiPM_lyte_Curr_;
    EnthalpyMolar_lyte_Curr_   = r.EnthalpyMolar_lyte_Curr_;

    qSource_Cell_curr_        = r.qSource_Cell_curr_;
    qSource_Cell_accumul_     = r.qSource_Cell_accumul_;
    jouleHeat_lyte_Cell_curr_ = r.jouleHeat_lyte_Cell_curr_;
    jouleHeat_solid_Cell_curr_ = r.jouleHeat_solid_Cell_curr_;
    electrodeHeat_Cell_curr_ = r.electrodeHeat_Cell_curr_; 
    overPotentialHeat_Cell_curr_ = r.overPotentialHeat_Cell_curr_;
    deltaSHeat_Cell_curr_      = r.deltaSHeat_Cell_curr_;
    potentialAnodic_           = r.potentialAnodic_;
    potentialCathodic_         = r.potentialCathodic_;
    nEnthalpy_New_Cell_        = r.nEnthalpy_New_Cell_;
    nEnthalpy_Old_Cell_        = r.nEnthalpy_Old_Cell_;

    thermalCond_Cell_          = r.thermalCond_Cell_;
    valCellTmpsVect_Cell_      = r.valCellTmpsVect_Cell_;
    ivb_                       = r.ivb_;
  
    ionicLiquid_               = r.ionicLiquid_->duplMyselfAsThermoPhase();
    trans_                     = r.trans_->duplMyselfAsTransport();
    solidSkeleton_             = r.solidSkeleton_->duplMyselfAsThermoPhase();

    for (size_t i = 0; i < ExtraPhaseList_.size(); ++i) {
       delete ExtraPhaseList_[i];
    }
    ExtraPhaseList_            = r.ExtraPhaseList_;
    for (size_t i = 0; i < ExtraPhaseList_.size(); ++i) {
        ExtraPhaseList_[i]     = new ExtraPhase(*(r.ExtraPhaseList_[i]));
    }
    Porosity_prob_type_        = r.Porosity_prob_type_;

    return *this;
}
//=====================================================================================================================
void porousFlow_dom1D::domain_prep(LocalNodeIndices *li_ptr)
{
    /*
     * First call the parent domain prep to get the node information
     */
    BulkDomain1D::domain_prep(li_ptr);

    double domainThickness = BDD_ptr_->Xpos_end - BDD_ptr_->Xpos_start;
    double porosity = -1.0;
    double volumeSeparator = PSCinput_ptr->separatorArea_ * domainThickness;
    double volumeInert = 0.0;
    double volumeFractionInert = 0.0;
    //double mv = 0.0;
    if (solidSkeleton_) {
        if (PSCinput_ptr->separatorMass_ > 0.0) {
            volumeInert = PSCinput_ptr->separatorMass_ / solidSkeleton_->density() ;
	    volumeFractionInert = volumeInert / volumeSeparator;
        } else if (PSCinput_ptr->separatorSolid_vf_ > 0.0) {
	    volumeFractionInert = PSCinput_ptr->separatorSolid_vf_;
        }
        if (volumeFractionInert >= 1.0) {
            throw m1d_Error("porousFlow_dom1D::domain_prep()", 
                            "Input volume fraction of separator solid is greater than one: " + fp2str(volumeFractionInert));
        }
	ExtraPhase* ep = ExtraPhaseList_[0];
	if (!doubleEqual(ep->volFraction, volumeFractionInert)) {
	    throw m1d_Error("porousFlow_dom1D::domain_prep", "solidSkel");
	}
	ep->volFraction = volumeFractionInert;
    }
    size_t offS = 0;
    //
    // If there is a solidSkeleton ThermoPhase, then identify that with the first volume fraction of the extra condensed phases.
    // We'll keep the mole number and volume fraction in the extra phases lists.
    //
    numExtraCondensedPhases_ = ExtraPhaseList_.size();

    porosity_Cell_.resize(NumLcCells, porosity);
    porosity_Cell_old_.resize(NumLcCells, porosity);
    Temp_Cell_old_.resize(NumLcCells, TemperatureReference_);
    xdelCell_Cell_.resize(NumLcCells, -1.0);
    volumeFraction_Phases_Cell_.resize(NumLcCells*numExtraCondensedPhases_, 0.0);
    volumeFraction_Phases_Cell_old_.resize(NumLcCells*numExtraCondensedPhases_, 0.0);
    moleNumber_Phases_Cell_.resize(NumLcCells*numExtraCondensedPhases_, 0.0);
    moleNumber_Phases_Cell_old_.resize(NumLcCells*numExtraCondensedPhases_, 0.0);

    double cellThickness = 0.0;
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
            cellThickness = 0.5*(nodeRight->xNodePos() - nodeLeft->xNodePos());
            xdelCell_Cell_[iCell] = cellThickness;
      
	    porosity = 1.0;
	    /*
	    volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell] = volumeFractionInert;
	    volumeFraction_Phases_Cell_old_[numExtraCondensedPhases_ * iCell] = volumeFractionInert;
	    moleNumber_Phases_Cell_[numExtraCondensedPhases_ * iCell] = volumeFractionInert * cellThickness * mv;
	    moleNumber_Phases_Cell_old_[numExtraCondensedPhases_ * iCell] = volumeFractionInert * cellThickness * mv;
	    */
	    for (size_t k = 0; k < ExtraPhaseList_.size(); ++k) {
		ExtraPhase* ep = ExtraPhaseList_[k];
		ThermoPhase* tp = ep->tp_ptr;
		tp->setState_TP(temp_Curr_, pres_Curr_);
		double mvp = tp->molarVolume();
		volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + offS + k] = ep->volFraction;
		volumeFraction_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + offS + k] = ep->volFraction;
		moleNumber_Phases_Cell_[numExtraCondensedPhases_ * iCell + offS + k] = ep->volFraction * cellThickness * mvp;
		moleNumber_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + offS + k] = ep->volFraction * cellThickness * mvp;
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

    cellTmpsVect_Cell_.resize(NumLcCells);
    qSource_Cell_curr_.resize(NumLcCells, 0.0);
    qSource_Cell_accumul_.resize(NumLcCells, 0.0);
    jouleHeat_lyte_Cell_curr_.resize(NumLcCells, 0.0);
    jouleHeat_solid_Cell_curr_.resize(NumLcCells, 0.0);
    electrodeHeat_Cell_curr_.resize(NumLcCells, 0.0);
    overPotentialHeat_Cell_curr_.resize(NumLcCells, 0.0);
    deltaSHeat_Cell_curr_.resize(NumLcCells, 0.0);
    nEnthalpy_New_Cell_.resize(NumLcCells, 0.0);
    nEnthalpy_Old_Cell_.resize(NumLcCells, 0.0);

    valCellTmpsVect_Cell_.resize(NumLcCells);

    thermalCond_Cell_.resize(NumLcCells);

    size_t nsp = ionicLiquid_->nSpecies();

   
    mfElectrolyte_Soln_Curr_.resize(nsp, 0.0);
    mfElectrolyte_Thermo_Curr_.resize(nsp, 0.0);
    EnthalpyPM_lyte_Curr_.resize(nsp, 0.0);
    EnthalpyPhiPM_lyte_Curr_.resize(nsp, 0.0);
 
    if  (energyEquationProbType_) {
	CpMolar_lyte_Cell_.resize(NumLcCells, 0.0);
	CpMolar_solid_Cell_.resize(NumLcCells, 0.0);
	CpMolar_total_Cell_.resize(NumLcCells, 0.0);
	EnthalpyPM_lyte_Cell_.resize(NumLcCells * nsp, 0.0);
	EnthalpyMolar_lyte_Cell_.resize(NumLcCells, 0.0);
    }
}
//==================================================================================================================================
double porousFlow_dom1D::heatSourceLastStep() const
{
    double q = 0.0;
    for (size_t iCell = 0; iCell < (size_t) NumLcCells; ++iCell) {
	q +=  qSource_Cell_curr_[iCell];
    }
    return q;
}
//==================================================================================================================================
double porousFlow_dom1D::heatSourceAccumulated() const
{
    double q = 0.0;
    for (size_t iCell = 0; iCell < (size_t) NumLcCells; ++iCell) {
	q += qSource_Cell_accumul_[iCell];
    }
    return q;
}
//==================================================================================================================================
void porousFlow_dom1D::heatSourceZeroAccumulated() const
{
    for (size_t iCell = 0; iCell < (size_t) NumLcCells; ++iCell) {
	qSource_Cell_accumul_[iCell] = 0.0;
    }
}
//==================================================================================================================================
void porousFlow_dom1D::residSetupTmps()
{
    size_t index_CentLcNode;

    NodalVars *nodeCent = nullptr;
    NodalVars *nodeLeft = nullptr;
    NodalVars *nodeRight = nullptr;

    size_t indexCent_EqnStart;
    size_t indexLeft_EqnStart;
    size_t indexRight_EqnStart;

    int  index_LeftLcNode;
    int  index_RightLcNode;

    cellTmpsVect_Cell_.resize(NumLcCells);
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell;

        cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
        NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
        NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;
        NodeTmps& nodeTmpsRight  = cTmps.NodeTmpsRight_;
        /*
         *  ---------------- Get the index for the center node ---------------------------------
         */
        index_CentLcNode = Index_DiagLcNode_LCO[iCell];

        /*
         *   Get the pointer to the NodalVars object for the center node
         */
        nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
        nodeTmpsCenter.nv = nodeCent;
        cTmps.nvCent_ = nodeCent;
        /*
         *  Index of the first equation in the bulk domain of center node
         */
        indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
        nodeTmpsCenter.index_EqnStart = indexCent_EqnStart;

        /*
         * Offsets for the variable unknowns in the solution vector for the electrolyte domain
         */
        nodeTmpsCenter.Offset_Voltage              = nodeCent->indexBulkDomainVar0((size_t) Voltage);
        nodeTmpsCenter.Offset_MoleFraction_Species = nodeCent->indexBulkDomainVar0((size_t) MoleFraction_Species);
	nodeTmpsCenter.Offset_Displacement_Axial   = nodeCent->indexBulkDomainVar0((size_t) Displacement_Axial);
        nodeTmpsCenter.Offset_Velocity_Axial       = nodeCent->indexBulkDomainVar0((size_t) Velocity_Axial);
        nodeTmpsCenter.Offset_Temperature          = nodeCent->indexBulkDomainVar0((size_t) Temperature);
        nodeTmpsCenter.Offset_Pressure             = nodeCent->indexBulkDomainVar0((size_t) Pressure_Axial);
#ifdef MECH_MODEL
	nodeTmpsCenter.Offset_Solid_Stress_Axial   = nodeCent->indexBulkDomainVar0((size_t) Solid_Stress_Axial);
#endif
        nodeTmpsCenter.RO_Current_Conservation     = nodeCent->indexBulkDomainEqn0((size_t) Current_Conservation);
        nodeTmpsCenter.RO_Electrolyte_Continuity   = nodeCent->indexBulkDomainEqn0((size_t) Continuity);
        nodeTmpsCenter.RO_Species_Eqn_Offset       = nodeCent->indexBulkDomainEqn0((size_t) Species_Eqn_Offset);
        nodeTmpsCenter.RO_MFSum_offset             = nodeCent->indexBulkDomainEqn0((size_t) MoleFraction_Summation);
        nodeTmpsCenter.RO_ChargeBal_offset         = nodeCent->indexBulkDomainEqn0((size_t) ChargeNeutrality_Summation);
        nodeTmpsCenter.RO_Enthalpy_Conservation    = nodeCent->indexBulkDomainEqn0((size_t) Enthalpy_Conservation);

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
            nodeLeft = nullptr;
            /*
             *  If there is no left node, we assign the left solution index to the center solution index
             */
            indexLeft_EqnStart = indexCent_EqnStart;
            nodeTmpsLeft.index_EqnStart = indexLeft_EqnStart;

            nodeTmpsLeft.Offset_Voltage              = nodeTmpsCenter.Offset_Voltage;
            nodeTmpsLeft.Offset_MoleFraction_Species = nodeTmpsCenter.Offset_MoleFraction_Species;
            nodeTmpsLeft.Offset_Displacement_Axial   = nodeTmpsCenter.Offset_Displacement_Axial;
            nodeTmpsLeft.Offset_Velocity_Axial       = nodeTmpsCenter.Offset_Velocity_Axial;
            nodeTmpsLeft.Offset_Temperature          = nodeTmpsCenter.Offset_Temperature;
            nodeTmpsLeft.Offset_Pressure             = nodeTmpsCenter.Offset_Pressure;
#ifdef MECH_MODEL
	    nodeTmpsLeft.Offset_Solid_Stress_Axial   = nodeTmpsCenter.Offset_Solid_Stress_Axial;
#endif
            nodeTmpsLeft.RO_Current_Conservation     = nodeTmpsCenter.RO_Current_Conservation;
            nodeTmpsLeft.RO_Electrolyte_Continuity   = nodeTmpsCenter.RO_Electrolyte_Continuity;
            nodeTmpsLeft.RO_Species_Eqn_Offset       = nodeTmpsCenter.RO_Species_Eqn_Offset;
            nodeTmpsLeft.RO_MFSum_offset             = nodeTmpsCenter.RO_MFSum_offset;
            nodeTmpsLeft.RO_ChargeBal_offset         = nodeTmpsCenter.RO_ChargeBal_offset;
            nodeTmpsLeft.RO_Enthalpy_Conservation    = nodeTmpsCenter.RO_Enthalpy_Conservation;
        } else {
            // get the node structure for the left node
            nodeLeft = LI_ptr_->NodalVars_LcNode[index_LeftLcNode];
            //index of first equation in the electrolyte of the left node
            indexLeft_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_LeftLcNode];
            nodeTmpsLeft.index_EqnStart = indexLeft_EqnStart;

            nodeTmpsLeft.Offset_Voltage              = nodeLeft->indexBulkDomainVar0((size_t) Voltage);
            nodeTmpsLeft.Offset_MoleFraction_Species = nodeLeft->indexBulkDomainVar0((size_t) MoleFraction_Species);
            nodeTmpsLeft.Offset_Displacement_Axial   = nodeLeft->indexBulkDomainVar0((size_t) Displacement_Axial);
            nodeTmpsLeft.Offset_Velocity_Axial       = nodeLeft->indexBulkDomainVar0((size_t) Velocity_Axial);
            nodeTmpsLeft.Offset_Temperature          = nodeLeft->indexBulkDomainVar0((size_t) Temperature);
            nodeTmpsLeft.Offset_Pressure             = nodeLeft->indexBulkDomainVar0((size_t) Pressure_Axial);
#ifdef MECH_MODEL
	    nodeTmpsLeft.Offset_Solid_Stress_Axial   = nodeLeft->indexBulkDomainVar0((size_t) Solid_Stress_Axial);
#endif
            nodeTmpsLeft.RO_Current_Conservation     = nodeLeft->indexBulkDomainEqn0((size_t) Current_Conservation);
            nodeTmpsLeft.RO_Electrolyte_Continuity   = nodeLeft->indexBulkDomainEqn0((size_t) Continuity);
            nodeTmpsLeft.RO_Species_Eqn_Offset       = nodeLeft->indexBulkDomainEqn0((size_t) Species_Eqn_Offset);
            nodeTmpsLeft.RO_MFSum_offset             = nodeLeft->indexBulkDomainEqn0((size_t) MoleFraction_Summation);
            nodeTmpsLeft.RO_ChargeBal_offset         = nodeLeft->indexBulkDomainEqn0((size_t) ChargeNeutrality_Summation);
            nodeTmpsLeft.RO_Enthalpy_Conservation    = nodeLeft->indexBulkDomainEqn0((size_t) Enthalpy_Conservation);
        }
        cTmps.nvLeft_ = nodeLeft;
      /*
         * ------------------------ Get the indexes for the right node ------------------------------------
         */
        index_RightLcNode = Index_RightLcNode_LCO[iCell];
        if (index_RightLcNode < 0) {
            nodeRight = nullptr;
            /*
             *  If there is no right node, we assign the right solution index to the center solution index
             */
            indexRight_EqnStart = indexCent_EqnStart;
            nodeTmpsRight.index_EqnStart = indexRight_EqnStart;

            nodeTmpsRight.Offset_Voltage              = nodeTmpsCenter.Offset_Voltage;
            nodeTmpsRight.Offset_MoleFraction_Species = nodeTmpsCenter.Offset_MoleFraction_Species;
            nodeTmpsRight.Offset_Displacement_Axial   = nodeTmpsCenter.Offset_Displacement_Axial;
            nodeTmpsRight.Offset_Velocity_Axial       = nodeTmpsCenter.Offset_Velocity_Axial;
            nodeTmpsRight.Offset_Temperature          = nodeTmpsCenter.Offset_Temperature;
            nodeTmpsRight.Offset_Pressure             = nodeTmpsCenter.Offset_Pressure;
#ifdef MECH_MODEL
	    nodeTmpsRight.Offset_Solid_Stress_Axial   = nodeTmpsCenter.Offset_Solid_Stress_Axial;
#endif
            nodeTmpsRight.RO_Current_Conservation     = nodeTmpsCenter.RO_Current_Conservation;
            nodeTmpsRight.RO_Electrolyte_Continuity   = nodeTmpsCenter.RO_Electrolyte_Continuity;
            nodeTmpsRight.RO_Species_Eqn_Offset       = nodeTmpsCenter.RO_Species_Eqn_Offset;
            nodeTmpsRight.RO_MFSum_offset             = nodeTmpsCenter.RO_MFSum_offset;
            nodeTmpsRight.RO_ChargeBal_offset         = nodeTmpsCenter.RO_ChargeBal_offset;
            nodeTmpsRight.RO_Enthalpy_Conservation    = nodeTmpsCenter.RO_Enthalpy_Conservation;
        } else {
            //NodalVars
            nodeRight = LI_ptr_->NodalVars_LcNode[index_RightLcNode];
            //index of first equation of right node
            indexRight_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_RightLcNode];
            nodeTmpsRight.index_EqnStart = indexRight_EqnStart;

            nodeTmpsRight.Offset_Voltage              = nodeRight->indexBulkDomainVar0((size_t) Voltage);
            nodeTmpsRight.Offset_MoleFraction_Species = nodeRight->indexBulkDomainVar0((size_t) MoleFraction_Species);
            nodeTmpsRight.Offset_Displacement_Axial   = nodeRight->indexBulkDomainVar0((size_t) Displacement_Axial);
            nodeTmpsRight.Offset_Velocity_Axial       = nodeRight->indexBulkDomainVar0((size_t) Velocity_Axial);
            nodeTmpsRight.Offset_Temperature          = nodeRight->indexBulkDomainVar0((size_t) Temperature);
            nodeTmpsRight.Offset_Pressure             = nodeRight->indexBulkDomainVar0((size_t) Pressure_Axial);
#ifdef MECH_MODEL
	    nodeTmpsRight.Offset_Solid_Stress_Axial   = nodeRight->indexBulkDomainVar0((size_t) Solid_Stress_Axial);
#endif
            nodeTmpsRight.RO_Current_Conservation     = nodeRight->indexBulkDomainEqn0((size_t) Current_Conservation);
            nodeTmpsRight.RO_Electrolyte_Continuity   = nodeRight->indexBulkDomainEqn0((size_t) Continuity);
            nodeTmpsRight.RO_Species_Eqn_Offset       = nodeRight->indexBulkDomainEqn0((size_t) Species_Eqn_Offset);
            nodeTmpsRight.RO_MFSum_offset             = nodeRight->indexBulkDomainEqn0((size_t) MoleFraction_Summation);
            nodeTmpsRight.RO_ChargeBal_offset         = nodeRight->indexBulkDomainEqn0((size_t) ChargeNeutrality_Summation);
            nodeTmpsRight.RO_Enthalpy_Conservation    = nodeRight->indexBulkDomainEqn0((size_t) Enthalpy_Conservation);
        }
        cTmps.nvRight_ = nodeRight;    
        /*
         * --------------------------- CALCULATE POSITION AND DELTA_X Variables -----------------------------
         * Calculate the distance between the left and center node points
         */
        if (nodeLeft) {
            cTmps.xdelL_ = nodeCent->xNodePos() - nodeLeft->xNodePos();
            cTmps.xCellBoundaryL_ = 0.5 * (nodeLeft->xNodePos() + nodeCent->xNodePos());
        } else {
            cTmps.xdelL_ = 0.0;
            cTmps.xCellBoundaryL_ = nodeCent->xNodePos();
        }
        /*
         * Calculate the distance between the right and center node points
         */
        if (nodeRight == 0) {
            cTmps.xdelR_ = 0.0;
            cTmps.xCellBoundaryR_ = nodeCent->xNodePos();
        } else {
            cTmps.xdelR_ = nodeRight->xNodePos() - nodeCent->xNodePos();
            cTmps.xCellBoundaryR_ = 0.5 * (nodeRight->xNodePos() + nodeCent->xNodePos());
        }
        /*
         * Calculate the cell width
         */
        cTmps.xdelCell_ = cTmps.xCellBoundaryR_ - cTmps.xCellBoundaryL_;

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
 *
 */
void
porousFlow_dom1D::initialConditions(const bool doTimeDependentResid,
				    Epetra_Vector* soln_ptr,
				    Epetra_Vector* solnDot,
				    const double t,
				    const double delta_t)
{
    Epetra_Vector& soln = *soln_ptr;

    int index_CentLcNode;
    NodalVars* nodeCent = 0;
    int indexCent_EqnStart;
    PhaseList* pl_ptr = PSCinput_ptr->PhaseList_;

    int iph = pl_ptr->globalPhaseIndex(PSCinput_ptr->electrolytePhase_);
    if (iph < 0) {
	throw ZuzaxError("BDD_porousElectrode::BDD_porousElectrode()",
			   "Can't find the phase in the phase list: " + PSCinput_ptr->electrolytePhase_);
    }
    //ThermoPhase* tmpPhase = & (PSCinput_ptr->PhaseList_)->thermo(iph);
    //int nSp = tmpPhase->nSpecies();
#ifdef TRACK_LOCATION
    std::cout << " porousFlow_dom1D::initialConditions after pl_ptr->globalPhaseIndex(PSCinput_ptr->electrolytePhase_) "<<std::endl;
#endif
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell;
    
        index_CentLcNode = Index_DiagLcNode_LCO[iCell];
        // pointer to the NodalVars object for the center node
        nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
        // Index of the first equation in the bulk domain of center node
        indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];

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
         * Offsets for the variable unknowns in the solution vector for the electrolyte domain
         */
        size_t iVAR_Vaxial  = nodeCent->indexBulkDomainVar0((size_t) Velocity_Axial);
        //size_t iVar_Species = nodeCent->indexBulkDomainVar0((size_t) MoleFraction_Species);
        size_t iVar_Voltage = nodeCent->indexBulkDomainVar0((size_t) Voltage);
        size_t iVar_Temperature = nodeCent->indexBulkDomainVar0((size_t) Temperature);
        size_t iVar_Pressure = nodeCent->indexBulkDomainVar0((size_t) Pressure_Axial);
        size_t iVar_Voltage_ED = iVar_Voltage + 1;
        //
        // Find the start of the solution at the current node
        //
        //const double *solnCentStart = &(soln[indexCent_EqnStart]);

        soln[indexCent_EqnStart + iVAR_Vaxial] = 0.0;
        //
        // Set the temperature if it is part of the solution vector
        //      
        temp_Curr_ = PSinput.TemperatureReference_;
        if (iVar_Temperature != npos) {
            soln[indexCent_EqnStart + iVar_Temperature] = PSinput.TemperatureReference_;
        }
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
        int igECDMC  = PSinput.PhaseList_->globalSpeciesIndex("ECDMC");
        if (igECDMC < 0) {
            throw ZuzaxError("confused", "confused");
        }
        int igLip = PSinput.PhaseList_->globalSpeciesIndex("Li+");
        if (igLip < 0) {
            throw ZuzaxError("confused", "confused");
        }
        int igPF6m = PSinput.PhaseList_->globalSpeciesIndex("PF6-");
        if (igPF6m < 0) {
            throw ZuzaxError("confused", "confused");
        }


	//for (size_t i = 0; i < KRSpecies) 

	// soln[indexCent_EqnStart + iVar_Species + iECDMC_] = PSinput.electrolyteMoleFracs_[igECDMC];
	// soln[indexCent_EqnStart + iVar_Species + iLip_  ] = PSinput.electrolyteMoleFracs_[igLip];
        //soln[indexCent_EqnStart + iVar_Species + iPF6m_ ] = PSinput.electrolyteMoleFracs_[igPF6m];
        soln[indexCent_EqnStart + iVar_Voltage] = -0.07;

        //double icurr = PSinput.icurrDischargeSpecified_;
        double volt = PSinput.CathodeVoltageSpecified_;
        soln[indexCent_EqnStart + iVar_Voltage_ED] = volt;
        //
        //  Fill in phiElectroyte_Curr_ 
        //
	// getVoltages(nodeCent, solnCentStart);
        //
        //  fill in mfElectrolyte_Soln_Curr[]  mfElectrolyte_Thermo_Curr_[]
        //
	//   getMFElectrolyte_soln(nodeCent, solnCentStart);
        //
        //
        // update porosity as computed from electrode input
        //
	//  porosity_Cell_[iCell] = porosity;

	//
	// Porosity Set up
	// 
	double volumeSeparator = PSCinput_ptr->separatorArea_ * PSCinput_ptr->separatorThickness_;
	double volumeInert = 0.0;
	double volumeFractionInert = 0.0;
        int offS = 0;
	//double mv = 0.0;
	double porosity = 1.0;	
	//double domainThickness = BDD_ptr_->Xpos_end - BDD_ptr_->Xpos_start;
        if (solidSkeleton_) {
	    offS = 1;
	    solidSkeleton_->setState_TP(temp_Curr_, pres_Curr_);
	    if (PSCinput_ptr->separatorMass_ > 0.0 ) {
		volumeInert = PSCinput_ptr->separatorMass_ / solidSkeleton_->density() ;
		volumeFractionInert = volumeInert / volumeSeparator;
	    } else if (PSCinput_ptr->separatorSolid_vf_ > 0.0) {
		volumeFractionInert = PSCinput_ptr->separatorSolid_vf_;
	    } else {
		throw m1d_Error("initial_conditions", "volumeFractionInert not set");
	    }
            ExtraPhase* ep = ExtraPhaseList_[0];
            ep->volFraction = volumeFractionInert;
	}

	for (size_t k = 0; k < ExtraPhaseList_.size(); ++k) {
	    ExtraPhase* ep = ExtraPhaseList_[k];
	    ThermoPhase* tp = ep->tp_ptr;
	    tp->setState_TP(temp_Curr_, pres_Curr_);
	    double mvp = tp->molarVolume();
	    volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + offS + k] = ep->volFraction;
	    volumeFraction_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + offS + k] = ep->volFraction;
	    moleNumber_Phases_Cell_[numExtraCondensedPhases_ * iCell + offS + k] = ep->volFraction * cellThickness  * mvp;
	    moleNumber_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + offS + k] = ep->volFraction * cellThickness  * mvp;
	    porosity -= ep->volFraction;
	}
	porosity_Cell_[iCell] = porosity;
	porosity_Cell_old_[iCell] = porosity;

    } // iCell
}
//==================================================================================================================================
void
porousFlow_dom1D::advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector* const soln_ptr,
                                      const Epetra_Vector* solnDot_ptr, const Epetra_Vector* const solnOld_ptr,
                                      const double t, const double t_old)
{
}
//==================================================================================================================================
void
porousFlow_dom1D::residEval_PreCalc(const bool doTimeDependentResid,
                                             const Epetra_Vector* const soln_ptr,
                                             const Epetra_Vector* const solnDot_ptr,
                                             const Epetra_Vector* const solnOld_ptr,
                                             const double t,
                                             const double rdelta_t,
                                             const Zuzax::ResidEval_Type residType,
                                             const Zuzax::Solve_Type solveType)

{

    // iCell == 0 special left section

    const Epetra_Vector& soln = *soln_ptr;
 
    // Special Section to determine where to get temperature and pressure
    cellTmps& cTmps          = cellTmpsVect_Cell_[0];
    NodalVars* nodeCent = cTmps.nvCent_;
    bool haveTemp = true;
    size_t iVar_Temperature = nodeCent->indexBulkDomainVar0((size_t) Temperature);
    if (iVar_Temperature == npos) {
	haveTemp = false;
    }
  
    // size_t iVar_Pressure_Axial = nodeCent->indexBulkDomainVar0((size_t) Pressure_Axial);
    /*
    bool havePres = true;
    if (iVar_Pressure_Axial == npos) {
	havePres = false;
    }
    */
    
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell;

        cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];

	valCellTmps& valTmps = valCellTmpsVect_Cell_[iCell];

        NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
	NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;
        NodeTmps& nodeTmpsRight  = cTmps.NodeTmpsRight_;

	nodeCent = cTmps.nvCent_;

	// This is the right and left control volume boundary
	valTmps.AxialVeloc.left   = soln[nodeTmpsLeft.index_EqnStart   +  nodeTmpsLeft.Offset_Velocity_Axial];
	valTmps.AxialVeloc.center = soln[nodeTmpsCenter.index_EqnStart + nodeTmpsCenter.Offset_Velocity_Axial];

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


        // Calculate the thermal conductivity of the porous matrix if we are calculating the energy equation

	thermalCond_Cell_[iCell] = thermalCondCalc_PorMatrix();
   
    }

}
//===================================================================================================================================
void porousFlow_dom1D::getState_Lyte_atCell(size_t iCell, const Epetra_Vector_Ghosted& soln, std::vector<double>& state_Lyte)
{
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    NodalVars* nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
    getState_Lyte(nodeCent, &(soln[indexCent_EqnStart]), state_Lyte);
}
//===================================================================================================================================
void porousFlow_dom1D::getState_Lyte(const NodalVars* const nv, const double* const solnNode_Curr, std::vector<double>& state_Lyte)
{
    SetupThermoShop1(nv, solnNode_Curr);
    size_t nspLyte = ionicLiquid_->nSpecies();
    state_Lyte.resize(nspLyte + 3);
    state_Lyte[0] = temp_Curr_;
    state_Lyte[1] = pres_Curr_;
    for (size_t k = 0; k < nspLyte; ++k) {
        state_Lyte[2 + k] = mfElectrolyte_Thermo_Curr_[k];
    }
    state_Lyte[2+nspLyte] = phiElectrolyte_Curr_;
}
//===================================================================================================================================
/*
 *  Setup shop at a particular point in the domain, calculating intermediate quantites and updating Zuzax's objects
 *  
 *  All member data with the suffix, _Curr_, are updated by this function.
 */
void porousFlow_dom1D::SetupThermoShop1(const NodalVars* const nv, const double* const solnNode_Curr)
{
    if (porosityEquationProbType_  &  Porosity_EqnType_Status::CalculatedOutOfEqnSystem) {
	double vfo = volumeFractionOther(cIndex_cc_);
	porosity_Cell_[cIndex_cc_] = 1.0 - vfo;
    }
    porosity_Curr_ = porosity_Cell_[cIndex_cc_];
    updateElectrolyte(nv, solnNode_Curr);
}
//===================================================================================================================================
void
porousFlow_dom1D::updateElectrolyte(const NodalVars* const nv, const double* const solnNode_Curr)
{
    /*
     * Get the temperature: Check to see if the temperature is in the solution vector.
     *   If it is not, then use the reference temperature
     */
    temp_Curr_ = getPointTemperature(nv, solnNode_Curr);
    /*
     * Get the pressure
     */
    pres_Curr_ = getPointPressure(nv, solnNode_Curr);
    /*
     *  Assemble electrolyte mole fractions into mfElectrolyte_Thermo_Curr_[]
     */
    getMFElectrolyte_soln(nv, solnNode_Curr);
    /*
     *  assemble electrolyte potential into phiElectrolyte_Curr_
     */
    getVoltages(nv, solnNode_Curr);
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
void
porousFlow_dom1D::getVoltages(const NodalVars* const nv, const double* const solnNode_Curr)
{
    size_t indexVS = nv->indexBulkDomainVar0(Voltage);
    phiElectrolyte_Curr_ = solnNode_Curr[indexVS];
}
//=====================================================================================================================
void
porousFlow_dom1D::getMFElectrolyte_soln(const NodalVars* const nv, const double* const solnNode_Curr)
{
    size_t indexMF = nv->indexBulkDomainVar0(MoleFraction_Species);
    for (size_t k = 0; k < BDD_ptr_->nSpeciesElectrolyte_; ++k) {
	mfElectrolyte_Soln_Curr_[k] = solnNode_Curr[indexMF + k];
    }
    calcMFElectrolyte_Thermo(&mfElectrolyte_Soln_Curr_[0], &mfElectrolyte_Thermo_Curr_[0]);
}
//=====================================================================================================================
/*
 *  This complicated logic assures us that  mf_Thermo_Curr[] is
 *    * charge neutral
 *    * no negative mole fractions
 *    * sums to 1
 */
void porousFlow_dom1D::calcMFElectrolyte_Thermo(const double* const mf, double* const mf_Thermo) const
{
    double sum, z;  
    double posSum = 0.0;
    double negSum = 0.0;
    for (size_t k = 0; k < BDD_ptr_->nSpeciesElectrolyte_; ++k) {
	mf_Thermo[k] = std::max(mf[k], 0.0);
        z = ionicLiquid_->charge(k);
	if (z > 0.0) {
	    posSum += mf_Thermo[k] * z;
	} else if (z < 0.0) {
	    negSum -= mf_Thermo[k] * z;
	}
    }
    if (posSum > 0.0 || negSum > 0.0) {	
	if (negSum <= 0.0) {
	    for (size_t k = 0; k < BDD_ptr_->nSpeciesElectrolyte_; ++k) {
		if (z > 0.0) {
		    mf_Thermo[k] = 0.0;
		}
	    }
	} else if (posSum <= 0.0) {
	    for (size_t k = 0; k < BDD_ptr_->nSpeciesElectrolyte_; ++k) {
		if (z < 0.0) {
		    mf_Thermo[k] = 0.0;
		}
	    }
	} else {
	    if (posSum != negSum) {
		sum = 0.5 * (posSum + negSum);
		double pratio = sum / posSum;
		double nratio = sum / negSum;
		for (size_t k = 0; k < BDD_ptr_->nSpeciesElectrolyte_; ++k) {
		    z = ionicLiquid_->charge(k);
		    if (z > 0.0) {
			mf_Thermo[k] *= pratio;
		    } else if (z < 0.0) {
			mf_Thermo[k] *= nratio;
		    }
		}
	    } 
	}
    }
    sum = 0.0;
    for (size_t k = 0; k < BDD_ptr_->nSpeciesElectrolyte_; ++k) {
	sum += mf_Thermo[k];
    }
    if (sum != 1.0) {
	for (size_t k = 0; k < BDD_ptr_->nSpeciesElectrolyte_; ++k) {
	    mf_Thermo[k] /= sum;
	}
    }
}
//==================================================================================================================================
double porousFlow_dom1D::effResistanceLayer(double &potAnodic, double &potCathodic, double &voltOCV, double &current)
{
    voltOCV=0.0;
    return 0.0;
}
//==================================================================================================================================
// Calculate the thermal conductivity of the porous matrix at the current cell.
double
porousFlow_dom1D::thermalCondCalc_PorMatrix()
{
    // Brief research
    //   Thermal conductivity of polycrystalline graphite = 80 Watts m-1 K-1
    //   Thermal conductivity of water = 0.58 Watts m-1 K-1
    //   Thermal conductivity of paper = 0.05 Watts m-1 K-1
    //   Thermal conductivity of ethylene glycol = 0.25 Watts m-1 K-1
    //   Thermal conductivity of glass = 1.05 Watts m-1 K-1
    //   Thermal conductivity of copper = 400. Watts m-1 K-1
    //   Thermal conductivity of aluminum = 200. Watts m-1 K-1
    //                        of olive oil = 0.16 Watts m-1 K-1
    //                        of gravel = 0.7  Watts m-1 K-1
    // temp value
    //  Judging from above, we are not more that a factor of 2-3 off here. This may not be necessary to go into much
    //  more detail.
    //
    return 0.5;
}
//==================================================================================================================================
double porousFlow_dom1D::volumeFractionOther(size_t iCell)
{
    double vf = 0.0;
    for (size_t jPhase = 0; jPhase < numExtraCondensedPhases_; ++jPhase) {
	vf += volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + jPhase];
    }
    return vf;
}
//==================================================================================================================================
//
// Calculate the porosity of a single cell
//     This routine can handle thermal expansion of the stoichiometric phases
//      Uses:
//             temp_Curr_          Needs the current temperature and pressuer
//             pres_Curr_
//             cTmps.xdelCell_     Needs geometry of cell
//             crossSectionalArea_
//             moleNumber_Phases_Cell_[]  Needs to know current moles of Skeletal and Other Phases
//
double porousFlow_dom1D::calcPorosity(size_t iCell) 
{
   cellTmps& cTmps = cellTmpsVect_Cell_[iCell];
   double xdelCell = cTmps.xdelCell_;
   double p = 1.0;
   for (size_t jPhase = 0; jPhase < numExtraCondensedPhases_; ++jPhase) {
        ExtraPhase* ep = ExtraPhaseList_[jPhase];
	ThermoPhase* tp = ep->tp_ptr;
	tp->setState_TP(temp_Curr_, pres_Curr_);
	double mv = tp->molarVolume();
        double vf = moleNumber_Phases_Cell_[numExtraCondensedPhases_ * iCell + jPhase] / (mv * xdelCell);
	volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + jPhase] = vf;
        p -= vf;
   }
   return p;
}
//==================================================================================================================================
} //namespace m1d
//==================================================================================================================================
