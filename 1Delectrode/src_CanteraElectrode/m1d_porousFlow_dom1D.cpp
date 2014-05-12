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

using namespace std;

namespace m1d
{
 
  //=====================================================================================================================
  porousFlow_dom1D::porousFlow_dom1D(BulkDomainDescription & bdd) :
    BulkDomain1D(bdd),
    doEnthalpyEquation_(0),
    porosity_Cell_(0),
    porosity_Cell_old_(0),
    temp_Curr_(TemperatureReference_),
    pres_Curr_(PressureReference_),
    concTot_Curr_(0.0),
    phiElectrolyte_Curr_(0.0),
    porosity_Curr_(0.0)
  {
  }
  //=====================================================================================================================
  porousFlow_dom1D::porousFlow_dom1D(const porousFlow_dom1D &r) :
    BulkDomain1D(r.BDD_),
    doEnthalpyEquation_(0),
    porosity_Cell_(0),
    porosity_Cell_old_(0),
    temp_Curr_(TemperatureReference_),
    pres_Curr_(PressureReference_),
    concTot_Curr_(0.0),
    phiElectrolyte_Curr_(0.0),
    porosity_Curr_(0.0)
  {
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

    doEnthalpyEquation_ = r.doEnthalpyEquation_;
    porosity_Cell_ = r.porosity_Cell_;
    porosity_Cell_old_ = r.porosity_Cell_old_;
    temp_Curr_ = r.temp_Curr_;
    concTot_Curr_ = r.concTot_Curr_;
    pres_Curr_ = r.pres_Curr_;
    phiElectrolyte_Curr_ = r.phiElectrolyte_Curr_;
    porosity_Curr_ = r.porosity_Curr_;
    qSource_Cell_curr_ = r.qSource_Cell_curr_;
    qSource_Cell_accumul_ = r.qSource_Cell_accumul_;
    jouleHeat_lyte_Cell_curr_ = r.jouleHeat_lyte_Cell_curr_;
    jouleHeat_solid_Cell_curr_ = r.jouleHeat_solid_Cell_curr_;
    nEnthalpy_New_Cell_        = r.nEnthalpy_New_Cell_;
    nEnthalpy_Old_Cell_        = r.nEnthalpy_Old_Cell_;

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

    qSource_Cell_curr_.resize(NumLcCells, 0.0);
    qSource_Cell_accumul_.resize(NumLcCells, 0.0);
    jouleHeat_lyte_Cell_curr_.resize(NumLcCells, 0.0);
    jouleHeat_solid_Cell_curr_.resize(NumLcCells, 0.0);
    electrodeHeat_Cell_curr_.resize(NumLcCells, 0.0);
    overPotentialHeat_Cell_curr_.resize(NumLcCells, 0.0);
    deltaSHeat_Cell_curr_.resize(NumLcCells, 0.0);
    nEnthalpy_New_Cell_.resize(NumLcCells, 0.0);
    nEnthalpy_Old_Cell_.resize(NumLcCells, 0.0);

    cellTmpsVect_Cell_.resize(NumLcCells);
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
	q +=  qSource_Cell_accumul_[iCell];
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
//=====================================================================================================================
} //namespace m1d
//=====================================================================================================================
