/**
 * @file m1d_porousElectrode_dom1D.cpp
 */

/*
 *   $Id: m1d_porousElectrode_dom1D.cpp 564 2013-03-08 23:35:51Z hkmoffa $
 */

#include "m1d_porousElectrode_dom1D.h"

#include "Electrode.h"

using namespace std;
using namespace Cantera;

namespace m1d
{
 
//=====================================================================================================================
porousElectrode_dom1D::porousElectrode_dom1D(BDD_porousElectrode& bdd) :
    porousFlow_dom1D(bdd),
    BDD_PE_ptr_(0),
    Electrode_Cell_(0),
    maxElectrodeSubIntegrationSteps_(0),
    surfaceArea_Cell_(0),
    nEnthalpy_Electrode_New_Cell_(0),
    nEnthalpy_Electrode_Old_Cell_(0),
    jFlux_EnthalpyPhi_metal_trCurr_(0.0),
    EnthalpyPhiPM_metal_Curr_(0),
    elem_Solid_Old_Cell_()
{
    BDD_PE_ptr_ = static_cast<BDD_porousElectrode*>(&BDD_);

    BDD_porousElectrode* bdd_pe_ptr = &bdd;
    
    metalPhase_ = BDD_PE_ptr_->metalPhase_;
}
//=====================================================================================================================
porousElectrode_dom1D::porousElectrode_dom1D(const porousElectrode_dom1D &r) :
    porousFlow_dom1D((BDD_porousElectrode&)r.BDD_),
    BDD_PE_ptr_(0),
    Electrode_Cell_(0),
    maxElectrodeSubIntegrationSteps_(0),
    surfaceArea_Cell_(0),
    nEnthalpy_Electrode_New_Cell_(0),
    nEnthalpy_Electrode_Old_Cell_(0),
    jFlux_EnthalpyPhi_metal_trCurr_(0.0),
    EnthalpyPhiPM_metal_Curr_(0),
    elem_Solid_Old_Cell_()
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
    jFlux_EnthalpyPhi_metal_trCurr_  = r.jFlux_EnthalpyPhi_metal_trCurr_;
    EnthalpyPhiPM_metal_Curr_        = r.EnthalpyPhiPM_metal_Curr_;
    metalPhase_                      = r.metalPhase_;
    elem_Solid_Old_Cell_             = r.elem_Solid_Old_Cell_;
    
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
     */
    porousFlow_dom1D::domain_prep(li_ptr);
    
    Electrode_Cell_.resize(NumLcCells, 0);
    surfaceArea_Cell_.resize(NumLcCells, 0.0);
    nEnthalpy_Electrode_New_Cell_.resize(NumLcCells, 0.0);
    nEnthalpy_Electrode_Old_Cell_.resize(NumLcCells, 0.0);
    EnthalpyPhiPM_metal_Curr_.resize(1, 0.0);
 
    BDD_porousElectrode* bdde = static_cast<BDD_porousElectrode*>(&BDD_);
    Electrode*ee = bdde->Electrode_;
    int neSolid = ee->nElements();
    elem_Solid_Old_Cell_.resize(neSolid , NumLcCells, 0.0);
}
//====================================================================================================================
//  An electrode object must be created and initialized for every cell in the domain
/*
 * Create electrode objects for every cell.  
 * Correct the volume and number of moles of 
 * active material within each of these electrode 
 * objects to correspond to the discretized volume.
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
//=====================================================================================================================
} //namespace m1d
//=====================================================================================================================

