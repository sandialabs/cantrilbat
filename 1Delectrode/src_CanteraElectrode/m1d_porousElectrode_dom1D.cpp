/**
 * @file m1d_porousElectrode_dom1D.cpp
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */


#include "m1d_porousElectrode_dom1D.h"
#include "m1d_cellTmps_PorousFlow.h"
#include "m1d_valCellTmps_porousFlow.h"
#include "m1d_ProblemStatementCell.h"
#include "Electrode.h"
#include "m1d_exception.h"
#include "m1d_globals.h"

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

namespace m1d
{
 
//=====================================================================================================================
porousElectrode_dom1D::porousElectrode_dom1D(BDD_porousElectrode* bdd_pe_ptr) :
    porousFlow_dom1D(bdd_pe_ptr),
    BDD_PE_ptr_(bdd_pe_ptr),
    Electrode_Cell_(0),
    maxElectrodeSubIntegrationSteps_(0),
    surfaceArea_Cell_(0),
    nEnthalpy_Electrode_New_Cell_(0),
    nEnthalpy_Electrode_Old_Cell_(0),
    nVol_zeroStress_Electrode_Cell_(0),
    jFlux_EnthalpyPhi_metal_trCurr_(0.0),
    metalPhase_(0)
{
    // assign the BDD_porousElectrode object pointer associated with this child object to a dedicated variable
    //BDD_PE_ptr_ = static_cast<BDD_porousElectrode*>(&BDD_);

    // Assign the metal phase pointer 
    metalPhase_ = BDD_PE_ptr_->metalPhase_;
    numElectrodeSubCycles_Cell_.resize(30);  // Problematic Statement -> why?
    // Debugging
    wBufff_[0] = 0;
    wBufff_[1] = 0;
    wBufff_[2] = 0;
    wBufff_[3] = 0;
    wBufff_[4] = 0;
}
//=====================================================================================================================
porousElectrode_dom1D::porousElectrode_dom1D(const porousElectrode_dom1D &r) :
    porousFlow_dom1D(r.BDD_PE_ptr_),
    BDD_PE_ptr_(r.BDD_PE_ptr_),
    Electrode_Cell_(0),
    maxElectrodeSubIntegrationSteps_(0),
    surfaceArea_Cell_(0),
    nEnthalpy_Electrode_New_Cell_(0),
    nEnthalpy_Electrode_Old_Cell_(0),
    nVol_zeroStress_Electrode_Cell_(0),
    jFlux_EnthalpyPhi_metal_trCurr_(0.0),
    metalPhase_(0)
{
    operator=(r);
}
//=====================================================================================================================
porousElectrode_dom1D::~porousElectrode_dom1D()
{
    for (size_t iCell = 0; iCell < Electrode_Cell_.size(); iCell++) {
	delete Electrode_Cell_[iCell];
	Electrode_Cell_[iCell] = 0;
    }
}
//=====================================================================================================================
porousElectrode_dom1D &
porousElectrode_dom1D::operator=(const porousElectrode_dom1D &r)
{
    if (this == &r) {
	return *this;
    }
    // Call the parent assignment operator
    porousFlow_dom1D::operator=(r);
    
    BDD_PE_ptr_ = r.BDD_PE_ptr_;
    // first do a shallow pointer copy
    Electrode_Cell_ = r.Electrode_Cell_;
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        // Dupl routine  is not implemented totally
        Electrode_Cell_[iCell] = (r.Electrode_Cell_[iCell])->duplMyselfAsElectrode();
    }
    
    maxElectrodeSubIntegrationSteps_ = r.maxElectrodeSubIntegrationSteps_;
    surfaceArea_Cell_                = r.surfaceArea_Cell_;
    nEnthalpy_Electrode_New_Cell_    = r.nEnthalpy_Electrode_New_Cell_;
    nEnthalpy_Electrode_Old_Cell_    = r.nEnthalpy_Electrode_Old_Cell_;
    nVol_zeroStress_Electrode_Cell_  = r. nVol_zeroStress_Electrode_Cell_;
    nVol_zeroStress_Electrode_Old_Cell_ = r.nVol_zeroStress_Electrode_Old_Cell_;
    jFlux_EnthalpyPhi_metal_trCurr_  = r.jFlux_EnthalpyPhi_metal_trCurr_;
    EnthalpyPhiPM_metal_Curr_        = r.EnthalpyPhiPM_metal_Curr_;
    metalPhase_                      = r.metalPhase_;
    elem_Solid_Old_Cell_             = r.elem_Solid_Old_Cell_;
    numElectrodeSubCycles_Cell_      = r.numElectrodeSubCycles_Cell_;
    
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
porousElectrode_dom1D::domain_prep(LocalNodeIndices *li_ptr)
{
    /*
     * First call the parent domain prep to get the node information
     * Also all arrays defined by parent objects are sized appropriately.
     * Also, an initial attempt is made to calculate the porosity and the volumeFraction_Phases is made using the reference
     * temperature and pressure to calculate molar volumes.
     */
    numElectrodeSubCycles_Cell_.resize(30);  // Problematic Statement -> why?
    porousFlow_dom1D::domain_prep(li_ptr);
    /*
     *  Size the arrays that are defined in this object
     */
    Electrode_Cell_.resize(NumLcCells, NULL);
    surfaceArea_Cell_.resize(NumLcCells, 0.0);
    nEnthalpy_Electrode_New_Cell_.resize(NumLcCells, 0.0);
    nEnthalpy_Electrode_Old_Cell_.resize(NumLcCells, 0.0);
    nVol_zeroStress_Electrode_Cell_.resize(NumLcCells, 0.0);
    nVol_zeroStress_Electrode_Old_Cell_.resize(NumLcCells, 0.0);
    EnthalpyPhiPM_metal_Curr_.resize(1, 0.0);
 
    //BDD_porousElectrode* bdde = static_cast<BDD_porousElectrode*>(&BDD_);
    Electrode* ee =  BDD_PE_ptr_->Electrode_;
    size_t neSolid = ee->nElements();
    size_t nl =  NumLcCells;
    elem_Solid_Old_Cell_.resize(neSolid , nl, 0.0);
    //int nn = NumLcCells;
    //numElectrodeSubCycles_Cell_.resize(nn, 0);

    /*
     *  Note, in child objects we would call instantiate_ElectrodeCells(). However, we don't have enough
     *  information to create electrode objects here.
     */
}
//====================================================================================================================
//  An electrode object must be created and initialized for every cell in the domain
/*
 * Create electrode objects for every cell. Correct the volume and number of moles of 
 * active material within each of these electrode objects to correspond to the discretized volume.
 */
void
porousElectrode_dom1D::instantiateElectrodeCells() 
{
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
porousElectrode_dom1D::advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector* soln_ptr,
                                             const Epetra_Vector* solnDot_ptr, const Epetra_Vector* solnOld_ptr,
                                             const double t, const double t_old)
{

    porousFlow_dom1D::advanceTimeBaseline(doTimeDependentResid, soln_ptr, solnDot_ptr, solnOld_ptr, t, t_old);

    const Epetra_Vector& soln = *soln_ptr;
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell; 

	cellTmps& cTmps      = cellTmpsVect_Cell_[iCell];

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

        porosity_Cell_old_[iCell] = porosity_Curr_;
        Temp_Cell_old_[iCell] = temp_Curr_;

	/*
        double* mfElectrolyte_Soln_old = mfElectrolyte_Soln_Cell_old_.ptrColumn(iCell);
        for (size_t k = 0; k < (size_t) nsp_; ++k) {
            mfElectrolyte_Soln_old[k] = mfElectrolyte_Soln_Curr_[k];
        }
	*/
	/*
         * Tell the electrode object to accept the current step and prep for the next step.
         *
         * We might at this point do a final integration to make sure we nailed the conditions of the last step.
         * However, we will hold off at implementing this right now
         */
        Electrode* ee = Electrode_Cell_[iCell];
        ee->resetStartingCondition(t);

        size_t neSolid = ee->nElements();
        for (size_t elem = 0; elem < neSolid; ++elem) {
            elem_Solid_Old_Cell_(elem,iCell) = ee->elementSolidMoles("Li");
        }
        //
        //  this is needed for a proper startup 
        //
        ee->updateState();
        //
        //  this is needed for a proper startup - sync initinit with final
        //
        ee->setInitStateFromFinal(true);

        if (energyEquationProbType_ == 3) {
            double volCellNew = cTmps.xdelCell_;
            // double volElectrodeCell = solidVolCell / crossSectionalArea_;
            double solidEnthalpy = ee->SolidEnthalpy() / crossSectionalArea_;
            double solidEnthalpyNew = solidEnthalpy;

            double lyteMolarEnthalpyNew = ionicLiquid_->enthalpy_mole();
            double volLyteNew = porosity_Curr_ * volCellNew;
            double lyteEnthalpyNew =  lyteMolarEnthalpyNew * concTot_Curr_ * volLyteNew;

            double nEnthalpy_New  = solidEnthalpyNew + lyteEnthalpyNew;

            if (! checkDblAgree( nEnthalpy_New, nEnthalpy_New_Cell_[iCell] ) ) {
                throw m1d_Error("", "Disagreement on new enthalpy calc");
            }

            nEnthalpy_Old_Cell_[iCell] = nEnthalpy_New_Cell_[iCell];
            nEnthalpy_Electrode_Old_Cell_[iCell] = nEnthalpy_Electrode_New_Cell_[iCell];
        }

	nVol_zeroStress_Electrode_Old_Cell_[iCell] = nVol_zeroStress_Electrode_Old_Cell_[iCell];

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
//====================================================================================================================
int porousElectrode_dom1D::getMaxSubGridTimeSteps() const
{
    return maxElectrodeSubIntegrationSteps_;
}
//=====================================================================================================================
double porousElectrode_dom1D::capacityPA(int platNum) const
{
    return 0.0;
}
//=====================================================================================================================
double porousElectrode_dom1D::capacityDischargedPA(int platNum) const
{
    return 0.0;
}
//=====================================================================================================================
double porousElectrode_dom1D::capacityLeftPA(int platNum, double voltsMax, double voltsMin) const
{
    return 0.0;
}
//=====================================================================================================================
double porousElectrode_dom1D::depthOfDischargePA(int platNum) const
{
    return 0.0;
}
//=====================================================================================================================
double porousElectrode_dom1D::depthOfDischargeStartingPA(int platNum) const
{
    return 0.0;
}
//=====================================================================================================================
void porousElectrode_dom1D::resetCapacityDischargedToDate() 
{
}
//=====================================================================================================================
// Return a value for the open circuit potential without doing a formally correct calculation
/*
 *  Currently this is defined as the open circuit potential on the outside electrode.
 *
 *   @return return the open circuit potential 
 */
double porousElectrode_dom1D::openCircuitPotentialQuick() const
{
    return 0.0;
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
double porousElectrode_dom1D::calcPorosity(size_t iCell) 
{
    cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
    double xdelCell = cTmps.xdelCell_;
    double volCell = xdelCell * crossSectionalArea_;
    size_t offS = 0;
    double mv, volS;
    double vf = 0.0;
    double vfE =  nVol_zeroStress_Electrode_Cell_[iCell] / volCell;
    double p = 1.0 - vfE;
    if (solidSkeleton_) {
        offS = 1;
      	solidSkeleton_->setState_TP(temp_Curr_, pres_Curr_);
        mv = solidSkeleton_->molarVolume();
        volS = moleNumber_Phases_Cell_[numExtraCondensedPhases_ * iCell] / mv;
        vf = volumeFraction_Phases_Cell_[iCell*numExtraCondensedPhases_] = volS / volCell;
    }
    
    for (size_t jPhase = 0; jPhase < numExtraCondensedPhases_; ++jPhase) {
        ExtraPhase* ep = ExtraPhaseList_[jPhase];
	ThermoPhase* tp = ep->tp_ptr;
	tp->setState_TP(temp_Curr_, pres_Curr_);
	mv = tp->molarVolume();
        volS = moleNumber_Phases_Cell_[numExtraCondensedPhases_ * iCell + offS + jPhase] / mv;
	volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + offS + jPhase] = volS / xdelCell;
        vf += volS / xdelCell;
    }
    p -= vf;
   return p; 
}
//==================================================================================================================================
void porousElectrode_dom1D::doPolarizationAdditions(double phiCurrentCollector, int region)
{
    bool dischargeDir = true;
    for (int iCell = 0; iCell < NumLcCells; ++iCell) {
         Electrode* ee = Electrode_Cell_[iCell];
         if (ee->doPolarizationAnalysis_) {
            (void) ee->polarizationAnalysisSurf(ee->polarSrc_list_Last_);
            for (size_t n = 0; n < ee->polarSrc_list_Last_.size(); ++n) {
                PolarizationSurfRxnResults& psr = ee->polarSrc_list_Last_[n];
                // Add contribution for addition of electrode's solid-phase conduction 
                psr.addSolidPol(phiCurrentCollector, region, dischargeDir);
            }

         }
    }

}
//==================================================================================================================================
} 
//----------------------------------------------------------------------------------------------------------------------------------

