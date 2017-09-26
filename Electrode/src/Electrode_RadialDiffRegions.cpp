/*
 * $Id: Electrode_RadialDiffRegions.cpp 298 2012-08-08 20:15:48Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"

#include "Electrode_RadialDiffRegions.h"
#include "cantera/integrators.h"

#include "BlockEntryGlobal.h"

using namespace std;
using namespace BEInput;
using namespace TKInput;

#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif

#include <set>
//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//====================================================================================================================
//====================================================================================================================
ELECTRODE_RadialDiffRegions_KEY_INPUT::ELECTRODE_RadialDiffRegions_KEY_INPUT(int printLvl) :
    ELECTRODE_KEY_INPUT(printLvl),
    numRegions_(1),
    numRegionsEntered_(0),
    solidDiffusionModel_(0),
    diffusiveFluxModel_(0),
    numRadialCellsRegions_(0),
    diffusionCoeffRegions_(0),
    rxnPerturbRegions_(0)
{
    numRadialCellsRegions_.resize(1, 5);
}
//====================================================================================================================
ELECTRODE_RadialDiffRegions_KEY_INPUT::~ELECTRODE_RadialDiffRegions_KEY_INPUT()
{
}
//====================================================================================================================
ELECTRODE_RadialDiffRegions_KEY_INPUT::
ELECTRODE_RadialDiffRegions_KEY_INPUT(const ELECTRODE_RadialDiffRegions_KEY_INPUT& right) :
    ELECTRODE_KEY_INPUT(right),
    numRegions_(right.numRegions_),
    numRegionsEntered_(right.numRegionsEntered_),
    solidDiffusionModel_(right.solidDiffusionModel_),
    diffusiveFluxModel_(right.diffusiveFluxModel_)
{
    diffusionCoeffRegions_          = right.diffusionCoeffRegions_;
    numRadialCellsRegions_          = right.numRadialCellsRegions_;
    rxnPerturbRegions_              = right.rxnPerturbRegions_;
    rregions_                       = right.rregions_;
}
//====================================================================================================================
ELECTRODE_RadialDiffRegions_KEY_INPUT&
ELECTRODE_RadialDiffRegions_KEY_INPUT::operator=(const ELECTRODE_RadialDiffRegions_KEY_INPUT& right)
{
    if (this == &right) {
        return *this;
    }
    ELECTRODE_KEY_INPUT::operator=(right);

    numRegions_                     = right.numRegions_;
    numRegionsEntered_              = right.numRegionsEntered_;
    solidDiffusionModel_            = right.solidDiffusionModel_; 
    diffusiveFluxModel_             = right.diffusiveFluxModel_;
    numRadialCellsRegions_          = right.numRadialCellsRegions_;
    diffusionCoeffRegions_          = right.diffusionCoeffRegions_;
    rxnPerturbRegions_              = right.rxnPerturbRegions_;
    rregions_                       = right.rregions_;

    return *this;
}
//====================================================================================================================
void ELECTRODE_RadialDiffRegions_KEY_INPUT::setup_input_child1(BEInput::BlockEntry* cf)
{
    /*  -------------------------------------------------------------------------------------------------------
     * Obtain the number of regions
     */
    LE_OneInt* s1 = new LE_OneInt("Number of Regions", &(numRegions_), 1, "numRegions");
    s1->set_default(1);
    cf->addLineEntry(s1);
    BaseEntry::set_SkipUnknownEntries(3);

    /*  -------------------------------------------------------------------------------------------------------
     * Obtain the diffusion model
     *   - Don't eliminate, but populate thinking about the case of LiFePO4
     */
    LE_OneInt* sdm = new LE_OneInt("Solid Diffusion Model", &(solidDiffusionModel_), 0, "solidDiffusionModel");
    sdm->set_default(0);
    cf->addLineEntry(sdm);
    BaseEntry::set_SkipUnknownEntries(3);

    /*  -------------------------------------------------------------------------------------------------------
     *   Diffusive flux model
     *            
     */
    LE_OneInt* dfm = new LE_OneInt("Diffusive Flux model", &(diffusiveFluxModel_), 0, "diffusiveFluxModel");
    dfm->set_default(0);
    cf->addLineEntry(dfm);

    // Not finished yet
    BaseEntry::set_SkipUnknownEntries(3);
}
//====================================================================================================================
void ELECTRODE_RadialDiffRegions_KEY_INPUT::setup_input_child2(BEInput::BlockEntry* cf)
{
    /*
     * Obtain the number of regions
     */
    int reqd = 0;
    if (solidDiffusionModel_) {
        reqd = 1;
    }

    /*
     *  Resize the output regions
     */
    if (numRegions_ != 1) {
       numRadialCellsRegions_.resize(numRegions_, 5);
    } 
    diffusionCoeffRegions_.resize(numRegions_, 1.0E-12);
    /*
     * Put a potential loop over regions here
     */ 
    // - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    BE_MultiBlockNested* be_rdr = new BE_MultiBlockNested("Radial Diffusion Region", &numRegionsEntered_, 1, cf);
    cf->addSubBlock(be_rdr);

    // - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    LE_OneInt* s1 = new LE_OneInt("Index of the Region", 0, 0, "indexRegion");
    s1->set_default(0);
    be_rdr->addLineEntry(s1);

    // - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    int iCell = 0;
    LE_OneInt* nRm = new LE_OneInt("Number of Cells in Region", &(numRadialCellsRegions_[iCell]), 0, "numRadialCellsRegions");
    nRm->set_default(5);
    be_rdr->addLineEntry(nRm);

    // - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    LE_MultiCStr* pid = new LE_MultiCStr("Phase Names within Distributed region", 0, 10, 1, 1, "phaseNames");
    be_rdr->addLineEntry(pid);

    /*
     * Make the diffusion coefficients optional unless the solid diffusion model is turned on
     */
    reqd = 0;
    if (solidDiffusionModel_) {
        reqd = numRegions_  + 1;
    }
    // - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    reqd = 0;
    LE_OneDbl* dDDC = new LE_OneDbl("Default Diffusion Coefficient", &(diffusionCoeffRegions_[iCell]), reqd,
				    "defaultDiffusionCoefficient");
    dDDC->set_default(1.0E-12);
    be_rdr->addLineEntry(dDDC);

    // - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    LE_StdVecDblVarLength* rr1 = new LE_StdVecDblVarLength("Reaction Rate Constant Perturbations for Regions",
            &(rxnPerturbRegions_),
            numRegions_, reqd, "rxnPerturbRegions_");
    rr1->set_default(1.0);
    rr1->set_limits(1.0E30, 1.0E-30);
    be_rdr->addLineEntry(rr1);

    BaseEntry::set_SkipUnknownEntries(0);
}
//====================================================================================================================
void ELECTRODE_RadialDiffRegions_KEY_INPUT::setup_input_child3(BEInput::BlockEntry* cf)
{
}
//======================================================================================================================
void ELECTRODE_RadialDiffRegions_KEY_INPUT::post_input_child2(BEInput::BlockEntry* cf)
{
    cf->print_usage(0);
    rregions_.clear();
    /*
     *  Create an ELECTRODE_RadialRegion_KEY_INPUT() object and save it in a vector
     *  We will use it to store information
     */
    rregions_.push_back( ELECTRODE_RadialRegion_KEY_INPUT() );
    /*
     *  Search for the radial diffusion region block in the block entries
     */
    const BEInput::BlockEntry* be = cf->searchBlockEntry("Radial Diffusion Region", false);
    const BEInput::BlockEntry* be_cand;
    /*
     *  Collect a set of block entries for the Radial Diffusion Rgions
     */
    std::set<const BlockEntry*> cc = cf->collectBlockEntries("Radial Diffusion Region", false);

    std::set<const BlockEntry*>::iterator cc_ptr;
    for (cc_ptr = cc.begin(); cc_ptr != cc.end(); cc_ptr++) {
      be_cand = *cc_ptr;
      int numT = be_cand->get_NumTimesProcessed(); 
      if (numT > 0) {
          be = *cc_ptr;
          // be_cand->print_usage(0);
      }
    }

    const BEInput::BE_MultiBlockNested* be_rdr = dynamic_cast<const BE_MultiBlockNested*>(be);
    BEInput::LineEntry *le = be_rdr->searchLineEntry("Number of Cells in Region");
    BEInput::LE_OneInt* le_int = dynamic_cast<LE_OneInt*>(le);

    BEInput::LineEntry *lepn = be_rdr->searchLineEntry("Phase Names within Distributed Region");
    BEInput::LE_MultiCStr* lepn_int = dynamic_cast<LE_MultiCStr*>(lepn);
   
    int iCell = 0;
    int numCells = le_int->currentTypedValue();
    numRadialCellsRegions_[iCell] = numCells;

    const char** curPN = lepn_int->currentTypedValue();
    int numPS = lepn_int->get_NumTimesProcessed();
    ELECTRODE_RadialRegion_KEY_INPUT& r0 = rregions_[0];

    for (int i = 0; i < numPS; i++) {
       const char * cP = curPN[i];
       string sP(cP);
       int k = m_pl->globalPhaseIndex(sP);
       if (k < 0) {
	   throw Electrode_Error("ELECTRODE_RadialDiffRegions_KEY_INPUT::post_input_child2() ERROR",
				 "Could not find Phase Name within Distributed region, " + sP + 
				 ", in the list of available phases");
       }

       (r0.phaseIndeciseKRsolidPhases_).push_back(k);
    }

    le = be_rdr->searchLineEntry("Default Diffusion Coefficient");
    BEInput::LE_OneDbl* d_lddc = dynamic_cast<LE_OneDbl*>(le);
    if (d_lddc->get_NumTimesProcessed() > 0) {
	r0.defaultDiffusionCoeff_ = d_lddc->currentTypedValue();
    }

}
//======================================================================================================================
void ELECTRODE_RadialDiffRegions_KEY_INPUT::post_input_child3(BEInput::BlockEntry* cf)
{
}
//======================================================================================================================
Electrode_RadialDiffRegions::Electrode_RadialDiffRegions() :
    Electrode_Integrator(),
    numRadialRegions_(-1),
    RadialRegionList_(0),
    SurfaceRegionList_(0)
{

}
//======================================================================================================================
Electrode_RadialDiffRegions::Electrode_RadialDiffRegions(const Electrode_RadialDiffRegions& right) :
    Electrode_Integrator(),
    numRadialRegions_(-1),
    RadialRegionList_(0),
    SurfaceRegionList_(0)
{
    /*
     * Call the assignment operator.
     */
    operator=(right);
}
//======================================================================================================================
Electrode_RadialDiffRegions&
Electrode_RadialDiffRegions::operator=(const Electrode_RadialDiffRegions& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    Electrode_Integrator::operator=(right);

    for (int i = 0; i < numRadialRegions_; i++) {
	delete RadialRegionList_[i];
	delete SurfaceRegionList_[i];
    }
    numRadialRegions_   = right.numRadialRegions_;
    RadialRegionList_.resize(numRadialRegions_);
    SurfaceRegionList_.resize(numRadialRegions_);
    for (int i = 0; i < numRadialRegions_; i++) {
	RadialRegionList_[i] = new Electrode_RadialRegion(*(right.RadialRegionList_[i]));
	SurfaceRegionList_[i] = new Electrode_SurfaceRegion(*(right.SurfaceRegionList_[i]));
    }
  
    MolarVolume_Ref_               = right.MolarVolume_Ref_;
    rnodePos_final_final_          = right.rnodePos_final_final_;
    rnodePos_final_                = right.rnodePos_final_;
    rnodePos_init_                 = right.rnodePos_init_;
    rnodePos_init_init_            = right.rnodePos_init_init_;
    rRefPos_final_                 = right.rRefPos_final_;
    rRefPos_init_                  = right.rRefPos_init_;
    rRefPos_init_init_             = right.rRefPos_init_init_;
    fracNodePos_                   = right.fracNodePos_;

    return *this;
}
//======================================================================================================================
Electrode_RadialDiffRegions::~Electrode_RadialDiffRegions()
{
   for (int i = 0; i < numRadialRegions_; i++) {
	delete RadialRegionList_[i];
	delete SurfaceRegionList_[i];
    }
}
//======================================================================================================================
Electrode* Electrode_RadialDiffRegions::duplMyselfAsElectrode() const
{
    Electrode_RadialDiffRegions* dd = new Electrode_RadialDiffRegions(*this);
    return dd;
}
//======================================================================================================================
Electrode_Types_Enum Electrode_RadialDiffRegions::electrodeType() const
{
    return RADIAL_DIFF_REGIONS_ET;
}
//======================================================================================================================
int Electrode_RadialDiffRegions::electrode_model_create(ELECTRODE_KEY_INPUT* eibase)
{
    Electrode_Integrator::electrode_model_create(eibase);

    /*
     *  Downcast the Key input to make sure we are being fed the correct child object
     */
    ELECTRODE_RadialDiffRegions_KEY_INPUT* ei = dynamic_cast<ELECTRODE_RadialDiffRegions_KEY_INPUT*>(eibase);
    if (!ei) {
        throw Electrode_Error(" Electrode_RadiaDiffRegions::electrode_model_create()",
                              " Expecting a child ELECTRODE_KEY_INPUT object and didn't get it");
    }

    /*
     * Find the number of radial diffusion regions
     */
    numRadialRegions_ = ei->numRegions_;
    AssertTrace(numRadialRegions_ == 1);
    RadialRegionList_.clear();
    SurfaceRegionList_.clear();
    for (int nn = 0 ; nn <  numRadialRegions_; nn++) {
	Electrode_RadialRegion* rr = new Electrode_RadialRegion();
	RadialRegionList_.push_back(rr);

	Electrode_SurfaceRegion* sr = new Electrode_SurfaceRegion();
	SurfaceRegionList_.push_back(sr);
    }
 

    for (int nn = 0 ; nn <  numRadialRegions_; nn++) {
        Electrode_RadialRegion* rr =  RadialRegionList_[nn];
        rr->electrode_model_create(eibase);

        Electrode_SurfaceRegion* sr = SurfaceRegionList_[nn];
        sr->electrode_model_create(eibase);
    }

    /*
     * Initialize the arrays in this object now that we know the number of equations
     */
    init_sizes();

    /*
     *  Initialize the species
     */
    initializeAsEvenDistribution();

    return 0;
}
//=====================================================================================================================
int
Electrode_RadialDiffRegions::setInitialConditions(ELECTRODE_KEY_INPUT* eibase)
{

    /*
     *  Downcast the Key input to make sure we are being fed the correct child object
     */
    ELECTRODE_RadialDiffRegions_KEY_INPUT* ei = dynamic_cast<ELECTRODE_RadialDiffRegions_KEY_INPUT*>(eibase);
    if (!ei) {
        throw Electrode_Error(" Electrode_RadialDiffRegions::electrode_model_create()",
                              " Expecting a child ELECTRODE_KEY_INPUT object and didn't get it");
    }

    Electrode::setInitialConditions(ei);


    Electrode::setInitStateFromFinal(true);

    setPhaseExistenceForReactingSurfaces();
    // printf("spMoles_final_[2] = %g\n", spMoles_final_[2]);

    neq_ = nEquations();

    /*
     *  Call a routine to
     */
    create_solvers();

    return 0;
}
//====================================================================================================================
void
Electrode_RadialDiffRegions::init_sizes()
{
    neq_ = nEquations();
    initSizes();
}
//====================================================================================================================
void
Electrode_RadialDiffRegions::initializeAsEvenDistribution()
{

}
//==================================================================================================================================
bool Electrode_RadialDiffRegions::resetStartingCondition(double Tinitial, bool doResetAlways)
{
    bool resetToInitInit = Electrode_Integrator::resetStartingCondition(Tinitial, doResetAlways);
    return resetToInitInit;
}
//==================================================================================================================================
void Electrode_RadialDiffRegions::updateState()
{
}
//==================================================================================================================================
/*
 * There is a small dependence on mf_external and mf_internal exhibited by this function
 */
void  Electrode_RadialDiffRegions::extractInfoJustBorn(std::vector<size_t>& justBornMultiSpecies)
{
    updateState();
}
//====================================================================================================================
void Electrode_RadialDiffRegions::printElectrode(int pSrc, bool subTimeStep)
{
    double* netROP = new double[m_NumTotSpecies];
    double egv = TotalVol();
    printf("   ===============================================================\n");
    if (subTimeStep) {
        printf("      Electrode at intermediate-step time final = %g\n", tfinal_);
        printf("                   intermediate-step time init  = %g\n", tinit_);
    } else {
        printf("      Electrode at time final = %g\n", t_final_final_);
        printf("                   time init  = %g\n", t_init_init_);
    }
    printf("   ===============================================================\n");
    printf("          Number of external surfaces = %d\n", (int) numExternalInterfacialSurfaces_);
    printf("          Solid Volume = %10.3E\n", ElectrodeSolidVolume_);
    printf("          Total Volume = %10.3E\n", egv);
    printf("          Temperature = %g\n", temperature_);
    printf("          Pressure = %g\n", pressure_);
    double capacd = capacityDischarged();
    printf("          Capacity Discharged = %g coulombs = %g Ah\n", capacd, capacd / 3600.);
    printf("\n");
    printf("          followElectrolyteMoles = %d\n", followElectrolyteMoles_);
    printf("          ElectrolytePseudoMoles = %g\n",  electrolytePseudoMoles_);

    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        printElectrodePhase(iph, pSrc);
        printf("     ===============================================================\n");
    }
    delete [] netROP;
}
//===================================================================================================================
void Electrode_RadialDiffRegions::printElectrodePhase(size_t iph, int pSrc, bool subTimeStep)
{
    size_t isph = npos;
    double* netROP = new double[m_NumTotSpecies];
    ThermoPhase& tp = thermo(iph);
    std::string pname = tp.id();
    size_t istart = m_PhaseSpeciesStartIndex[iph];
    size_t nsp = tp.nSpecies();
    printf("     ===============================================================\n");
    printf("          Phase %d %s \n", static_cast<int>(iph), pname.c_str());
    printf("                Total moles = %g\n", phaseMoles_final_[iph]);
    if (iph == metalPhase_) {
        double deltaT = t_final_final_ - t_init_init_;
        if (subTimeStep) {
            deltaT = tfinal_ - tinit_;
        }
        if (deltaT > 1.0E-200) {
            double amps = spMoleIntegratedSourceTerm_[istart] / deltaT * Faraday;
            if (subTimeStep) {
                amps = spMoleIntegratedSourceTermLast_[istart] / deltaT * Faraday;
            }
            printf("                Current = %g amps \n", amps);
        } else {
            printf("                Current = NA amps \n");
        }
    }
    if (iph == metalPhase_ || iph == solnPhase_) {
        printf("                  Voltage = %g\n", tp.electricPotential());
    }
    if (iph >= m_NumVolPhases) {
        isph = iph - m_NumVolPhases;
        printf("                surface area (final) = %11.5E m2\n",  surfaceAreaRS_final_[isph]);
        printf("                surface area (init)  = %11.5E m2\n",  surfaceAreaRS_init_[isph]);
        int ddd =  isExternalSurface_[isph];
        printf("                IsExternalSurface = %d\n", ddd);
        double oc = openCircuitVoltage(isph);
        if (oc != 0.0) {
            printf("                 Open Circuit Voltage = %g\n", oc);
        }
    }
    printf("\n");
    printf("                Name               MoleFrac_final  kMoles_final kMoles_init SrcTermLastStep(kMoles)\n");
    for (size_t k = 0; k < nsp; k++) {
        std::string sname = tp.speciesName(k);
        if (pSrc) {
            if (subTimeStep) {
                printf("                %-22s %10.3E %10.3E   %10.3E  %10.3E\n", sname.c_str(), spMf_final_[istart + k],
                       spMoles_final_[istart + k], spMoles_init_[istart + k],
                       spMoleIntegratedSourceTermLast_[istart + k]);
            } else {
                printf("                %-22s %10.3E %10.3E   %10.3E  %10.3E\n", sname.c_str(), spMf_final_[istart + k],
                       spMoles_final_[istart + k], spMoles_init_init_[istart + k],
                       spMoleIntegratedSourceTerm_[istart + k]);
            }
        } else {
            if (subTimeStep) {
                printf("                %-22s %10.3E %10.3E   %10.3E\n", sname.c_str(), spMf_final_[istart + k],
                       spMoles_final_[istart + k],   spMoles_init_[istart + k]);
            } else {
                printf("                %-22s %10.3E %10.3E   %10.3E\n", sname.c_str(), spMf_final_[istart + k],
                       spMoles_final_[istart + k],   spMoles_init_init_[istart + k]);
            }
        }
    }
    if (iph >= m_NumVolPhases) {
        const std::vector<double>& rsSpeciesProductionRates = RSD_List_[isph]->calcNetSurfaceProductionRateDensities();
        RSD_List_[isph]->getNetRatesOfProgress(netROP);

        double* spNetProdPerArea = (double*) spNetProdPerArea_List_.ptrColumn(isph);
        std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
        size_t nphRS = RSD_List_[isph]->nPhases();
        size_t kIndexKin = 0;
        for (size_t kph = 0; kph < nphRS; kph++) {
            size_t jph = RSD_List_[isph]->kinOrder[kph];
            size_t istart = m_PhaseSpeciesStartIndex[jph];
            size_t nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
            for (size_t k = 0; k < nsp; k++) {
                spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
                kIndexKin++;
            }
        }
        printf("\n");
        printf("                           spName                  Source (kmol/m2/s) \n");
        for (size_t k = 0; k <  m_NumTotSpecies; k++) {
            std::string ss = speciesName(k);
            printf("                           %-22s %10.3E\n", ss.c_str(), spNetProdPerArea[k]);
        }
    }
    delete [] netROP;
}
//==================================================================================================================================
} // End of #ifdef useZuzaxNamespace
//----------------------------------------------------------------------------------------------------------------------------------
