/**
 * @file m1d_porousFlow_dom1D.cpp
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

//! Global Problem input structure
/*!
 *   This contains the input data for the problem.
 *   We've made it a global structure, as there is one and only one instance of the structure
 */
extern m1d::ProblemStatementCell PSinput;
#include "m1d_exception.h"


using namespace std;

namespace m1d
{
 
//=====================================================================================================================
porousFlow_dom1D::porousFlow_dom1D(BDD_porousFlow &bdd) :
    BulkDomain1D(bdd),
    BDT_ptr_(0),
    porosity_Cell_(0),
    porosity_Cell_old_(0),
    Temp_Cell_old_(0),
    cIndex_cc_(-1),
    temp_Curr_(TemperatureReference_),
#ifdef MECH_MODEL
    mm_stress_Curr_(SolidStressAxialRef_),
#endif
    pres_Curr_(PressureReference_),
    concTot_Curr_(0.0),
    phiElectrolyte_Curr_(0.0),
    porosity_Curr_(0.0),
    thermalCond_Curr_(0.0),
    heatFlux_Curr_(0.0), 
    jFlux_EnthalpyPhi_Curr_(0.0),
    EnthalpyMolar_lyte_Curr_(0.0),
    ivb_(VB_MOLEAVG)
{
    BDT_ptr_ = static_cast<BDD_porousFlow*>(&BDD_);
    ionicLiquid_ = BDT_ptr_->ionicLiquid_;
    trans_ = BDT_ptr_->trans_;
    solidSkeleton_ = BDT_ptr_->solidSkeleton_;

    energyEquationProbType_ = PSinput.Energy_equation_prob_type_;
    solidMechanicsProbType_ = PSinput.Solid_Mechanics_prob_type_;
}
//=====================================================================================================================
  porousFlow_dom1D::porousFlow_dom1D(const porousFlow_dom1D &r) :
      BulkDomain1D(r.BDD_),
      BDT_ptr_(0),
      porosity_Cell_(0),
      porosity_Cell_old_(0),
      Temp_Cell_old_(0),
      cIndex_cc_(-1),
      temp_Curr_(TemperatureReference_),
#ifdef MECH_MODEL
      mm_stress_Curr_(SolidStressAxialRef_),
#endif
      pres_Curr_(PressureReference_),
      concTot_Curr_(0.0),
      phiElectrolyte_Curr_(0.0),
      porosity_Curr_(0.0),
      thermalCond_Curr_(0.0),
      heatFlux_Curr_(0.0),
      jFlux_EnthalpyPhi_Curr_(0.0),
      EnthalpyMolar_lyte_Curr_(0.0),
      ivb_(VB_MOLEAVG)
  {
      BDT_ptr_ = static_cast<BDD_porousFlow*>(&BDD_);
      porousFlow_dom1D::operator=(r);
  }
//=====================================================================================================================
  porousFlow_dom1D::~porousFlow_dom1D()
  {
  }
  //=====================================================================================================================
  porousFlow_dom1D &
  porousFlow_dom1D::operator=(const porousFlow_dom1D &r)
  {
    if (this == &r) {
      return *this;
    }
    // Call the parent assignment operator
    BulkDomain1D::operator=(r);

    BDT_ptr_                  = r.BDT_ptr_;
    CpMolar_lyte_Cell_        = r.CpMolar_lyte_Cell_;
    CpMolar_solid_Cell_       = r.CpMolar_solid_Cell_;
    CpMolar_total_Cell_       = r.CpMolar_total_Cell_;
    EnthalpyPM_lyte_Cell_     = r.EnthalpyPM_lyte_Cell_;
    porosity_Cell_            = r.porosity_Cell_;
    porosity_Cell_old_        = r.porosity_Cell_old_;
    Temp_Cell_old_            = r.Temp_Cell_old_;
    cIndex_cc_                = r.cIndex_cc_;
    temp_Curr_                = r.temp_Curr_;
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
  
    ionicLiquid_               = r.ionicLiquid_;
    trans_                     = r.trans_;
    solidSkeleton_             = r.solidSkeleton_;

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
  porousFlow_dom1D::domain_prep(LocalNodeIndices *li_ptr)
  {
    /*
     * First call the parent domain prep to get the node information
     */
    BulkDomain1D::domain_prep(li_ptr);

    double porosity = -1.0;

    porosity_Cell_.resize(NumLcCells, porosity);
    porosity_Cell_old_.resize(NumLcCells, porosity);
    Temp_Cell_old_.resize(NumLcCells, TemperatureReference_);

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
//=====================================================================================================================
double porousFlow_dom1D::heatSourceLastStep() const
{
    double q = 0.0;
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
	q +=  qSource_Cell_curr_[iCell];
    }
    return q;
}
//=====================================================================================================================
double porousFlow_dom1D::heatSourceAccumulated() const
{
    double q = 0.0;
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
	q += qSource_Cell_accumul_[iCell];
    }
    return q;
}
//=====================================================================================================================
void porousFlow_dom1D::heatSourceZeroAccumulated() const
{
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
	qSource_Cell_accumul_[iCell] = 0.0;
    }
}
//=====================================================================================================================
void
porousFlow_dom1D::residSetupTmps()
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
        nodeTmpsCenter.Offset_Velocity_Axial       = nodeCent->indexBulkDomainVar0((size_t) Velocity_Axial);
        nodeTmpsCenter.Offset_Temperature          = nodeCent->indexBulkDomainVar0((size_t) Temperature);

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
            nodeLeft = 0;
            /*
             *  If there is no left node, we assign the left solution index to the center solution index
             */
            indexLeft_EqnStart = indexCent_EqnStart;
            nodeTmpsLeft.index_EqnStart = indexLeft_EqnStart;

            nodeTmpsLeft.Offset_Voltage              = nodeTmpsCenter.Offset_Voltage;
            nodeTmpsLeft.Offset_MoleFraction_Species = nodeTmpsCenter.Offset_MoleFraction_Species;
            nodeTmpsLeft.Offset_Velocity_Axial       = nodeTmpsCenter.Offset_Velocity_Axial;
            nodeTmpsLeft.Offset_Temperature          = nodeTmpsCenter.Offset_Temperature;
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
            nodeTmpsLeft.Offset_Velocity_Axial       = nodeLeft->indexBulkDomainVar0((size_t) Velocity_Axial);
            nodeTmpsLeft.Offset_Temperature          = nodeLeft->indexBulkDomainVar0((size_t) Temperature);
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
            nodeRight = 0;
            /*
             *  If there is no right node, we assign the right solution index to the center solution index
             */
            indexRight_EqnStart = indexCent_EqnStart;
            nodeTmpsRight.index_EqnStart = indexRight_EqnStart;

            nodeTmpsRight.Offset_Voltage              = nodeTmpsCenter.Offset_Voltage;
            nodeTmpsRight.Offset_MoleFraction_Species = nodeTmpsCenter.Offset_MoleFraction_Species;
            nodeTmpsRight.Offset_Velocity_Axial       = nodeTmpsCenter.Offset_Velocity_Axial;
            nodeTmpsRight.Offset_Temperature          = nodeTmpsCenter.Offset_Temperature;
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
            nodeTmpsRight.Offset_Velocity_Axial       = nodeRight->indexBulkDomainVar0((size_t) Velocity_Axial);
            nodeTmpsRight.Offset_Temperature          = nodeRight->indexBulkDomainVar0((size_t) Temperature);
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
	throw CanteraError("BDD_porousElectrode::BDD_porousElectrode()",
			   "Can't find the phase in the phase list: " + PSCinput_ptr->electrolytePhase_);
    }
    //ThermoPhase* tmpPhase = & (PSCinput_ptr->PhaseList_)->thermo(iph);
    //int nSp = tmpPhase->nSpecies();

    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell;
    
        index_CentLcNode = Index_DiagLcNode_LCO[iCell];
        // pointer to the NodalVars object for the center node
        nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
        // Index of the first equation in the bulk domain of center node
        indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];

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


	//for (size_t i = 0; i < KRSpecies) 

	// soln[indexCent_EqnStart + iVar_Species + iECDMC_] = PSinput.electrolyteMoleFracs_[igECDMC];
	// soln[indexCent_EqnStart + iVar_Species + iLip_  ] = PSinput.electrolyteMoleFracs_[igLip];
        //soln[indexCent_EqnStart + iVar_Species + iPF6m_ ] = PSinput.electrolyteMoleFracs_[igPF6m];
        soln[indexCent_EqnStart + iVar_Voltage] = -0.07;

        //double icurr = PSinput.icurrDischargeSpecified_;
        double volt = PSinput.CathodeVoltageSpecified_;
        soln[indexCent_EqnStart + iVar_Voltage_ED] = volt;
        //
        //  Fill in phiElectroyte_Curr_ and phiElectrode_Curr_
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


    }
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
porousFlow_dom1D::advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector* soln_ptr,
                                      const Epetra_Vector* solnDot_ptr, const Epetra_Vector* solnOld_ptr,
                                      const double t, const double t_old)
{
}
//========================================================================================================================
void
porousFlow_dom1D::residEval_PreCalc(const bool doTimeDependentResid,
                                             const Epetra_Vector* soln_ptr,
                                             const Epetra_Vector* solnDot_ptr,
                                             const Epetra_Vector* solnOld_ptr,
                                             const double t,
                                             const double rdelta_t,
                                             const ResidEval_Type_Enum residType,
                                             const Solve_Type_Enum solveType)

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

//=====================================================================================================================
//  Setup shop at a particular point in the domain, calculating intermediate quantites
//  and updating Cantera's objects
/*
 *  
 *  All member data with the suffix, _Curr_, are updated by this function.
 *
 * @param soln_Curr  Current value of the solution vector at the current node
 */
void
porousFlow_dom1D::SetupThermoShop1(const NodalVars* const nv, const doublereal* const soln_Curr)
{
    porosity_Curr_ = porosity_Cell_[cIndex_cc_];
    updateElectrolyte(nv, soln_Curr);
}
//=====================================================================================================================
void
porousFlow_dom1D::updateElectrolyte(const NodalVars* const nv, const doublereal* const solnElectrolyte_Curr)
{
    /*
     * Get the temperature: Check to see if the temperature is in the solution vector.
     *   If it is not, then use the reference temperature
     */
    temp_Curr_ = getPointTemperature(nv, solnElectrolyte_Curr);  
    /*
     * Get the pressure
     */
    pres_Curr_ = getPointPressure(nv, solnElectrolyte_Curr);
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
void
porousFlow_dom1D::getVoltages(const NodalVars* const nv, const double* const solnElectrolyte_Curr)
{
    size_t indexVS = nv->indexBulkDomainVar0(Voltage);
    phiElectrolyte_Curr_ = solnElectrolyte_Curr[indexVS];
}
//=====================================================================================================================
void
porousFlow_dom1D::getMFElectrolyte_soln(const NodalVars* const nv, const double* const solnElectrolyte_Curr)
{
    size_t indexMF = nv->indexBulkDomainVar0(MoleFraction_Species);
    for (size_t k = 0; k < BDT_ptr_->nSpeciesElectrolyte_; ++k) {
	mfElectrolyte_Soln_Curr_[k] = solnElectrolyte_Curr[indexMF + k];
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
    for (size_t k = 0; k < BDT_ptr_->nSpeciesElectrolyte_; ++k) {
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
	    for (size_t k = 0; k < BDT_ptr_->nSpeciesElectrolyte_; ++k) {
		if (z > 0.0) {
		    mf_Thermo[k] = 0.0;
		}
	    }
	} else if (posSum <= 0.0) {
	    for (size_t k = 0; k < BDT_ptr_->nSpeciesElectrolyte_; ++k) {
		if (z < 0.0) {
		    mf_Thermo[k] = 0.0;
		}
	    }
	} else {
	    if (posSum != negSum) {
		sum = 0.5 * (posSum + negSum);
		double pratio = sum / posSum;
		double nratio = sum / negSum;
		for (size_t k = 0; k < BDT_ptr_->nSpeciesElectrolyte_; ++k) {
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
    for (size_t k = 0; k < BDT_ptr_->nSpeciesElectrolyte_; ++k) {
	sum += mf_Thermo[k];
    }
    if (sum != 1.0) {
	for (size_t k = 0; k < BDT_ptr_->nSpeciesElectrolyte_; ++k) {
	    mf_Thermo[k] /= sum;
	}
    }
}
//=====================================================================================================================
double porousFlow_dom1D::effResistanceLayer(double &potAnodic, double &potCathodic, double &voltOCV, double &current)
{
    voltOCV=0.0;
    return 0.0;
}
//=====================================================================================================================
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
     //                        of gravel = 0.7  .7  .7  .7  .7  .7  .7  Watts m-1 K-1
     // temp value
     //  Judging from above, we are not more that a factor of 2-3 off here. This may not be necessary to go into much
     //  more detail.
     //
     return 0.5;
 }
//=====================================================================================================================
} //namespace m1d
//=====================================================================================================================
