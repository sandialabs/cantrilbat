/*
 * $Id: m1d_ProblemStatementCell.cpp 564 2013-03-08 23:35:51Z hkmoffa $
 */

#include "m1d_globals.h"
#include "m1d_ProblemStatementCell.h"
#include "m1d_CanteraElectrodeGlobals.h"
#include "m1d_exception.h"

#include "BlockEntryGlobal.h"

#include "Electrode_Factory.h"
#include "importPL.h"
#include "mdp_allo.h"

#define USE_DAKOTA
#ifdef USE_DAKOTA
#include "../../Experiment/src/exp_DakotaInterface.h"
#endif
using namespace std;
using namespace BEInput;
using namespace TKInput;

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

namespace m1d
{
//==================================================================================================================================
//!  Global storage for this global pointer
/*!
 *   This pointer needs to be properly taken care of in the calling program
 */
ProblemStatementCell* PSCinput_ptr = 0;

//==================================================================================================================================
ExtraPhase::ExtraPhase() :
    phaseName(""),
    canteraFileName(""),
    regions("all"),
    volFraction(0.0),
    tp_ptr(0)
{    
    for (size_t i = 0; i < 10; ++i) {
        bregionID[i] = npos;
        sregionID[i] = npos;
    }
}
//==================================================================================================================================
ExtraPhase::ExtraPhase(const ExtraPhase& b) :
    phaseName(b.phaseName),
    canteraFileName(b.canteraFileName),
    regions(b.regions),
    volFraction(b.volFraction),
    tp_ptr(b.tp_ptr)
{    
    for (size_t i = 0; i < 10; ++i) {
        bregionID[i] = b.bregionID[i];
        sregionID[i] = b.sregionID[i];
    }
    if (b.tp_ptr) {
        tp_ptr = (b.tp_ptr)->duplMyselfAsThermoPhase();
    }
}
//==================================================================================================================================
ProblemStatementCell::ProblemStatementCell() :
  ProblemStatement(), 
  NumberCanteraFiles(0), CanteraFileNames(0),
  anodeBCType_(0), 
  cathodeBCType_(0), 
  cathodeBCFile_("BC.xml"),
  icurrDischargeSpecified_(1.0),
  CathodeVoltageSpecified_(1.9),
  rootFinderForConstantCurrent_(0),
  anodeFile_("anode.inp"), cathodeFile_("cathode.inp"),
  electrolytePhase_("LiKCl_electrolyte"), separatorPhase_("MgO"),
  separatorMass_(0.0), separatorArea_(-1.0),
  separatorThickness_(70.0E-6), 
  separatorDiameter_(-1.0),
  separatorSolid_vf_(-1.0),
  conductivityAnode_(1.0E6),
  conductivityCathode_(100.),
  anodeCCThickness_(0.0),
  cathodeCCThickness_(0.0),
  extraCathodeResistance_(0.0),
  ResistanceLoad_(0.0),
  VoltageLoad_(0.0),
  useDakota_(false),
  maxSubgridTimeStepsPerGlobalTimeStep_(100),
  fromDakotaFileName_("params.in"), toDakotaFileName_("results.out"),
  nTotPhases_(0), nTotSpecies_(0), nTotElements_(0),
  SpeciesNames_(0),
  PhaseNames_(0), 
  ElementNames_(0),
  initDefaultNumCVsAnode_(10),
  initDefaultNumCVsCathode_(10),
  initDefaultNumCVsSeparator_(10),
  doHeatSourceTracking_(0),
  doResistanceTracking_(0),
  anodeTempBCType_(-1),
  anodeTempRef_(298.15),
  anodeHeatTranCoeff_(1000.),
  cathodeTempBCType_(-1),
  cathodeTempRef_(298.15),
  cathodeHeatTranCoeff_(1000.),
  Pressure_formulation_prob_type_(0),
  artificialCompressibilityInvAtm_(0.0),
  numExtraPhases_(0),
  ExtraPhaseList_(0)
{
    PhaseList_ = new ZZCantera::PhaseList();

    struct ExtraPhase* ep = new ExtraPhase();
    ExtraPhaseList_.push_back(ep); 
}
//==================================================================================================================================
ProblemStatementCell::~ProblemStatementCell()
{
    if (CanteraFileNames) {
        for (int i = 0; CanteraFileNames[i] != 0; i++) {
            mdpUtil::mdp_safe_free((void **) &(CanteraFileNames[i]));
        }
        mdpUtil::mdp_safe_free((void **) &(CanteraFileNames));
    }

    if (PhaseList_){
        delete PhaseList_;
        PhaseList_ = 0;
    }

    if (anode_input_) {
        delete anode_input_;
        anode_input_ = 0;
    }

    if (cathode_input_) {
        delete cathode_input_;
        cathode_input_ = 0;
    }
    if (BC_TimeDep_) {
        delete BC_TimeDep_;
        BC_TimeDep_ = 0;
    }

    for (size_t i = 0; i < ExtraPhaseList_.size(); ++i) {
        delete ExtraPhaseList_[i];
    }
}
//==================================================================================================================================
void
ProblemStatementCell::setup_input_pass1(BlockEntry *cf)
{
  PSCinput_ptr = this;
  ProblemStatement::setup_input_pass1(cf);
  /*
   * Obtain the number of cantera files to be read
   */
  LE_OneInt *s1 = new LE_OneInt("Number of Cantera Files", &(NumberCanteraFiles), 0, "NumCanteraFiles");
  s1->set_default(0);
  cf->addLineEntry(s1);
  BaseEntry::set_SkipUnknownEntries(3);

  /* -------------------------------------------------------------------------
   *  
   *   Cross Sectional Area
   *
   *      Change this card to being mandatory if the coordinate system is rectilinear 
   */
  LineEntry* cs1 = cf->searchLineEntry("Cross Sectional Area");
  LineEntry* cpl = cf->searchLineEntry("Coordinate System");
  BI_DepIntMaxMin* dmm2 = new BI_DepIntMaxMin(cpl, BIDT_INTMAXMIN, 0, 0, BIDRT_ONENUMTR);
  cs1->declareDependency(dmm2);

  /* -------------------------------------------------------------------------
   *
   *     Heat Source Tracking
   *
   *   0 = Do not do any heat source tracking
   *   1 = calculate heat sources, print them out, and have them available for external communication
   */
  int reqd = 0;
  LE_OneInt *i2 = new LE_OneInt("Heat Source Tracking", &(doHeatSourceTracking_), reqd, "doHeatSourceTracking");
  i2->set_default(0);
  i2->set_limits(9, 0);
  cf->addLineEntry(i2);

  /* -------------------------------------------------------------------------
   *
   *     Resistance Tracking
   *
   *   0 = Do not do any resistance tracking
   *   1 = calculate resistances for layers, print them out, and have them available for external communication
   */
  LE_OneInt *ir = new LE_OneInt("Resistance Tracking", &(doResistanceTracking_), reqd, "doResistanceTracking");
  ir->set_default(0);
  ir->set_limits(1, 0);
  cf->addLineEntry(ir);

  /* --------------------------------------------------------------------
   * 
   * Pressure Formulation Problem Type
   *
   *   0  = None (default)
   *   1  = Darcy Flow
   *   2  = Darcy Flow with Gas Reservoir
   *   3  = Partial Saturation
   */
  const char *pressureFormList[4] = {"None", "Darcy Flow", "Darcy Flow with Gas Reservoir", "Partial Saturation"};
  LE_PickList *lepkp = new LE_PickList("Pressure Formulation Problem Type", &Pressure_formulation_prob_type_,
				       pressureFormList, 4, 0, "Pressure_formulation_prob_type_");
  lepkp->set_default(0);
  cf->addLineEntry(lepkp);

}
//===================================================================================================================================
void
ProblemStatementCell::setup_input_pass2(BlockEntry *cf)
{
    ProblemStatement::setup_input_pass2(cf);

    LineEntry *sle1 = 0;
    /*
     *  Get the input deck for
     *  Cantera description of the model.
     */
    LE_MultiCStr *s1 = new LE_MultiCStr("Cantera File Name", &CanteraFileNames, 1, 1, 0, "CanteraFileNames");
    s1->set_default("gas.cti");
    /*
     * Set up a dependency on the input from the Number of cantera
     * Files card
     */
    sle1 = cf->searchLineEntry("Number of Cantera Files");
    int numF = 1;
    bool okbefore = sle1->ansDepCheckOneInt(numF);
    if (okbefore) {
	// printf("Found it before\n");
	s1->set_NumTimesRequired(numF);
    } else {
	printf("Num Lines not in input deck -> no dependency\n");
	// Note -> this is not right -> should be one or more dependencies.
    }
    cf->addLineEntry(s1);

    /* ------------------------------------------------------------------------
     *  Option to specify if you want to use the root finder for constant current applications
     *        0  - no
     *        1  - yes
     *        2  - Yes, but only when in trouble (under construction)
     */
    int reqd = 0;
    LE_OneInt *iRoot = new LE_OneInt("Root Finder for Constant Current", &(rootFinderForConstantCurrent_), reqd,
				     "rootFinderForConstantCurrent");
    iRoot->set_default(0);
    iRoot->set_limits(2, 0);
    cf->addLineEntry(iRoot);

    /* -------------------------------------------------------------------------
     *
     * Anode Temp BC Type - int (not required)
     *    -1 - not set
     *     0 - Specify a fixed temperature at the anode
     *     1 - Specify a fixed heat flux through the battery
     *     2 - Time dependent temp at the anode
     *    10 - Specify heat transfer coeff temperature bc.
     */
    reqd = 0;
    LE_OneInt *at = new LE_OneInt("Anode TemperatureBC Type", &(anodeTempBCType_), reqd, "anodeTemp_bc_type");
    at->set_default(-1);
    at->set_limits(10, 0);
    cf->addLineEntry(at);

    /* -------------------------------------------------------------------------
     *  Anode Collector Temperature
     */
    reqd = 0;
    LE_OneDbl *datemp = new LE_OneDbl("Anode Collector Temperature", &(anodeTempRef_), reqd, "anodeTempCollector");
    datemp->set_default(298.15);
    cf->addLineEntry(datemp);

    /* -------------------------------------------------------------------------
     *  Anode Heat Transfer Coeff
     */
    reqd = 0;
    LE_OneDbl *dacoeff = new LE_OneDbl("Anode Heat Transfer Coefficient", &(anodeHeatTranCoeff_), reqd, "anodeHeatTranCoeff");
    dacoeff->set_default(1000.);
    cf->addLineEntry(dacoeff);

    /* -------------------------------------------------------------------------
     *
     * Cathode Temp BC Type - int (not required)
     *    -1 - not set
     *     0 - Specify a fixed temperature at the anode
     *     1 - Specify a fixed heat flux through the battery
     *     2 - Time dependent temp at the cathode
     *    10 - Specify heat transfer coeff temperature bc.
     */
    reqd = 0;
    LE_OneInt *ct = new LE_OneInt("Cathode TemperatureBC Type", &(cathodeTempBCType_), reqd, "cathodeTemp_bc_type");
    ct->set_default(-1);
    ct->set_limits(10, 0);
    cf->addLineEntry(ct);

    /* -------------------------------------------------------------------------
     *  Cathode Collector Temperature
     */
    reqd = 0;
    LE_OneDbl *dctemp = new LE_OneDbl("Cathode Collector Temperature", &(cathodeTempRef_), reqd, "cathodeTempCollector");
    dctemp->set_default(298.15);
    cf->addLineEntry(dctemp);

    /* -------------------------------------------------------------------------
     *  Cathode Heat Transfer Coeff
     */
    reqd = 0;
    LE_OneDbl *dccoeff = new LE_OneDbl("Cathode Heat Transfer Coefficient", &(cathodeHeatTranCoeff_), reqd, "cathodeHeatTranCoeff_");
    dccoeff->set_default(1000.);
    cf->addLineEntry(dccoeff);

    /* ------------------------------------------------------------------------------------------------------------
     *  Extra Phase
     */
    int numPhasesRequired = 0;
    BE_MultiBlockNested* be_exPh = new BE_MultiBlockNested("Extra Phase", &numExtraPhases_, numPhasesRequired, cf);
    cf->addSubBlock(be_exPh);
    ExtraPhase* ep_ptr = ExtraPhaseList_[0];
    /* ------------------------------------------------------------------------------------------------------------
     * Extra Phases: Phase Name
     *          name of the phases
     */
    LE_OneStr* ppn = new LE_OneStr("Phase Name", &(ep_ptr->phaseName), 10, 1, 1, "epPhaseName");
    be_exPh->addLineEntry(ppn);

    /* ------------------------------------------------------------------------------------------------------------
     * Extra Phases: Phase Name
     *          name of the phases
     */
    LE_OneStr* pcfn = new LE_OneStr("Cantera File Name ", &(ep_ptr->canteraFileName), 10, 1, 1, "epCanteraFileName");
    be_exPh->addLineEntry(pcfn);

    /* ------------------------------------------------------------------------------------------------------------
     * Extra Phases: Regions
     *          name of the regions
     *             anode, cathode, separator, all
     *             default = "all"
     */
    LE_OneStr* pcreg = new LE_OneStr("Region", &(ep_ptr->regions), 10, 1, 1, "epCanteraFileName");
    pcreg->set_default("all");
    be_exPh->addLineEntry(pcreg);

    /* ------------------------------------------------------------------------------------------------------------
     * Extra Phases: Regions
     *          name of the regions
     *             anode, cathode, separator, all
     *             default = "all"
     */
    LE_OneDbl* pvf = new LE_OneDbl("Volume Fraction", &(ep_ptr->volFraction), 1, "epVolumeFraction");
    pvf->set_default(0.0);
    be_exPh->addLineEntry(pvf);

}
//===================================================================================================================================
void
ProblemStatementCell::setup_input_pass3(BlockEntry *cf)
{

  //Need to call parent class to setup for parent class data members
  ProblemStatement::setup_input_pass3(cf);

  int reqd = 0;


  /* -------------------------------------------------------------------------
   *
   * Cathode BC Type - int (required)
   *     0 - Specify a fixed voltage at the cathode
   *     1 - Specify a fixed discharge current through the battery
   *     2 - Time dependent voltage at the cathode
   *     3 - Time dependent current at the cathode
   *     4 - specify time dependent voltage BoundaryCondition BCconstant
   *     5 - specify time dependent current BoundaryCondition BCconstant
   *     6 - specify time dependent voltage BoundaryCondition BCsteptable
   *     7 - specify time dependent current BoundaryCondition BCsteptable
   *     8 - specify time dependent voltage BoundaryCondition BClineartable
   *     9 - specify time dependent current BoundaryCondition BClineartable
   *    10 - Extra resistance or closed-loop boundary condition
   *    11 - External Load
   */
  reqd = 1;
  LE_OneInt *i2 = new LE_OneInt("Cathode BC Type", &(cathodeBCType_), reqd, "cathode_bc_type");
  i2->set_default(0);
  i2->set_limits(11, 0);
  cf->addLineEntry(i2);

  /* ------------------------------------------------------------------------
   * Name of XML formatted input file for the cathode boundary condition
   */
  reqd = 0;
  LE_OneStr *s7 = new LE_OneStr("Input Cathode BC File",
  				 &cathodeBCFile_, 1, 1, reqd,
  				 "cathodeBCFile" );
  s7->set_default("BC.xml");
  cf->addLineEntry(s7);

  /* -------------------------------------------------------------------------
   *
   * Discharge Current - double [Amps / cm2] (no default)
   * (conditionally required if Problem Type = 1)
   *
   * Specify the current in amps per cm2 for the battery.
   * This is only needed when the Problem Type is 1
   */
  BE_UnitConversion *ucCurr = new BE_UnitConversion();

  LE_OneDblUnits *d1 = new LE_OneDblUnits("Discharge Current", &(icurrDischargeSpecified_), 0, "icurrDischargeSpecified", ucCurr);
  d1->set_default(1.0);
  d1->set_limits(1.0E6, 0.0);
  cf->addLineEntry(d1);

  /* -------------------------------------------------------------------------
   *
   * Cathode Voltage - double [volts] (no default) (required)
   *
   * Specify the initial cathode voltage. Note the anode voltage is
   * always set to 0.0.
   * For a problem type of 0, this is the Dirichlet condition on the
   * cathode voltage. For a problem type of 1, this value is used to
   * set the initial condition of the cathode voltage.
   */
  BE_UnitConversion *ucVolt = new BE_UnitConversion();
  LE_OneDblUnits *d2 = new LE_OneDblUnits("Cathode Voltage", &(CathodeVoltageSpecified_), 0, "voltageSpecified", ucVolt);
  d2->set_default(1.9);
  d2->set_limits(6.7, -1.0);
  cf->addLineEntry(d2);

  /* ------------------------------------------------------------------------
   * Name of Electrolyte phase
   * This phase should be found in one of the listed Cantera files.
   */
  reqd = 0;
  LE_OneStr *s1 = new LE_OneStr("Electrolyte Phase",
				 &electrolytePhase_, 1, 1, reqd,
				 "electrolytePhase" );
  s1->set_default("LiKCl_electrolyte");
  cf->addLineEntry(s1);

  /* ------------------------------------------------------------------------
   * Electrolyte composition in mole fraction units
   */
  electrolyteMoleFracs_ = mdpUtil::mdp_alloc_dbl_1( nTotSpecies_, 0.0);

  reqd = 0;
  int construct = 1;
  BE_MoleComp *mf1 = new BE_MoleComp("Electrolyte Mole Fractions",
				      &electrolyteMoleFracs_, reqd,
				      SpeciesNames_, nTotSpecies_,
				      construct, "electrolyteMoleFracs_",
				      cf );
  mf1->generateDefLE();
  cf->addSubBlock(mf1);

  for ( int i = 0; i < nTotSpecies_; i++ )
    std::cout << "Species name " << i
	      << " is " << SpeciesNames_[i] << std::endl;

   /* ------------------------------------------------------------------------
   * Name of Anode Input file
   */
  reqd = 0;
  LE_OneStr *s2 = new LE_OneStr("Input Anode File",
				 &anodeFile_, 1, 1, reqd,
				 "anodeFile" );
  s2->set_default("anode.inp");
  cf->addLineEntry(s2);

  /* ------------------------------------------------------------------------
   * Name of Cathode Input file
   */
  reqd = 0;
  LE_OneStr *s3 = new LE_OneStr("Input Cathode File",
				 &cathodeFile_, 1, 1, reqd,
				 "cathodeFile" );
  s3->set_default("cathode.inp");
  cf->addLineEntry(s3);

  /* ------------------------------------------------------------------------
   * Name of Separator Phase
   */
  reqd = 0;
  LE_OneStr *s4 = new LE_OneStr("Separator Phase",
				 &separatorPhase_, 1, 1, reqd,
				 "separatorPhase" );
  s4->set_default("MgO");
  cf->addLineEntry(s4);

  /* -------------------------------------------------------------------------
   * Separator Mass
   * This phase should be found in one of the listed Cantera files.
   */
  BE_UnitConversion *ucMass = new BE_UnitConversion();
  reqd = 0;
  LE_OneDblUnits *d3 = new LE_OneDblUnits("Separator Mass", &(separatorMass_), reqd, "separatorMass", ucMass);

  d3->set_default(0.0);
  cf->addLineEntry(d3);


  /* -------------------------------------------------------------------------
   * Separator Area
   */
  BE_UnitConversion *ucL1 = new BE_UnitConversionLength();
  reqd = 0;
  LE_OneDblUnits *d4 = new LE_OneDblUnits("Separator Area", &(separatorArea_), reqd, "separatorArea", ucL1);

  d4->set_default(-1.0);
  cf->addLineEntry(d4);

  /* -------------------------------------------------------------------------
   * Separator Thickness
   */
  BE_UnitConversion *ucL2 = new BE_UnitConversionLength();
  reqd = 0;
  LE_OneDblUnits *d5 = new LE_OneDblUnits("Separator Thickness", &(separatorThickness_), reqd, "separatorThickness", ucL2);
  d5->set_default(70.0E-6);
  cf->addLineEntry(d5);

  /* -------------------------------------------------------------------------
   * Separator Diameter
   */
  BE_UnitConversion *ucL3 = new BE_UnitConversionLength();
  reqd = 0;
  LE_OneDblUnits *d6 = new LE_OneDblUnits("Separator Diameter", &(separatorDiameter_), reqd, "separatorDiameter", ucL3);
  d6->set_default(-1.0);
  cf->addLineEntry(d6);

  // If we specify the area we cannot specify the diameter
  BI_Dependency * dep_sepDia_sepArea = new BI_Dependency(d4, BIDT_ENTRYPROCESSED, BIDRT_ANTITHETICAL_ERROR);
  d6->declareDependency(dep_sepDia_sepArea);

  /* -------------------------------------------------------------------------
   * Separator Solid Skeleton Volume Fraction
   * This is one way to specify the porosity - we make it the default unless Separator Mass is specified.
   */
  reqd = 1;
  LE_OneDbl *dvf = new LE_OneDbl("Separator Solid Skeleton Volume Fraction", &(separatorSolid_vf_), reqd, "separatorSoild_vf");
  dvf->set_default(-1.0);
  cf->addLineEntry(dvf);

  // If we specify the separator mass we cannot specify the volume fraction
  BI_Dependency* dep_sepMass_sepVF = new BI_Dependency(d3, BIDT_ENTRYPROCESSED, BIDRT_ANTITHETICAL_ERROR);
  dvf->declareDependency(dep_sepMass_sepVF);

  // If we specify the separator vf we cannot specify the mass
  BI_Dependency* dep_sepVF_sepMass = new BI_Dependency(dvf, BIDT_ENTRYPROCESSED, BIDRT_ANTITHETICAL_ERROR);
  d3->declareDependency(dep_sepVF_sepMass);

  // If we specify the separator mass, then this card becomes optional
  BI_Dependency* dep_sepMassOpt_sepVF = new BI_Dependency(d3, BIDT_ENTRYPROCESSED, BIDRT_ZERONUMTIMESREQUIRED);
  dvf->declareDependency(dep_sepMassOpt_sepVF);

  /* ------------------------------------------------------------------------------------------------------------------
   *  Anode Current Collector thickness
   */
  reqd = 0;
  BE_UnitConversion *ucL4 = new BE_UnitConversionLength();
  LE_OneDblUnits *dacc = new LE_OneDblUnits("Anode Current Collector Thickness", &(anodeCCThickness_), reqd,
                                            "anodeCCThickness", ucL4);
  dacc->set_default(0.0);
  cf->addLineEntry(dacc);

  /* ------------------------------------------------------------------------------------------------------------------
   *  Cathode Current Collector thickness
   */
  reqd = 0;
  BE_UnitConversion *ucL5 = new BE_UnitConversionLength();
  LE_OneDblUnits *dccc = new LE_OneDblUnits("Cathode Current Collector Thickness", &(cathodeCCThickness_), reqd,
                                            "cathodeCCThickness", ucL5);
  dccc->set_default(0.0);
  cf->addLineEntry(dccc);

  /* ------------------------------------------------------------------------------------------------------------------
   * Electrical Conductivity of the Anode Electrode =
   */
  reqd = 0;
  LE_OneDbl *deca = new LE_OneDbl("Electrical Conductivity of the Anode Electrode", &(conductivityAnode_), reqd,
                                            "conductivityAnode_");
  deca->set_default(1.0E6);
  cf->addLineEntry(deca);

  /* ------------------------------------------------------------------------------------------------------------------
   * Electrical Conductivity of the Cathode Electrode =
   */
  reqd = 0;
  LE_OneDbl *decc = new LE_OneDbl("Electrical Conductivity of the Cathode Electrode", &(conductivityCathode_), reqd,
                                            "conductivityCathode_");
  decc->set_default(1.0E2);
  cf->addLineEntry(decc);


  /* ------------------------------------------------------------------------------------------------------------------
   *  Extra Resistance in Series with Cathode
   */
  reqd = 0;
  LE_OneDblUnits *derc = new LE_OneDblUnits("Extra Resistance in Series with Cathode", 
					    &(extraCathodeResistance_), reqd, "extraCathodeResistance");
  derc->set_default(0.0);
  cf->addLineEntry(derc);

  /* ------------------------------------------------------------------------------------------------------------------
   *  Extra Resistance characterized by load
   */
  reqd = 0;
  LE_OneDblUnits *derl = new LE_OneDblUnits("Load Resistance",
                                            &(ResistanceLoad_), reqd, "ResistanceLoad");
  derl->set_default(0.0);
  cf->addLineEntry(derl);

  /* ------------------------------------------------------------------------------------------------------------------
   *  Extra Resistance characterized by load
   */
  reqd = 0;
  LE_OneDblUnits *devl = new LE_OneDblUnits("Load Voltage",
                                            &(VoltageLoad_), reqd, "VoltageLoad");
  devl->set_default(0.0);
  cf->addLineEntry(devl);


  /* -------------------------------------------------------------------------------------------------------------------
   * Number of control volumes in anode
   */
  LE_OneInt *iCVA = new LE_OneInt("Number of CVs in Anode Domain", &initDefaultNumCVsAnode_);
  iCVA->set_default(10);
  iCVA->set_limits(1000, 1);
  cf->addLineEntry(iCVA);
 
  /* -------------------------------------------------------------------------------------------------------------------
   * Number of control volumes in anode
   */
  LE_OneInt *iCVC = new LE_OneInt("Number of CVs in Cathode Domain", &initDefaultNumCVsCathode_);
  iCVC->set_default(10);
  iCVC->set_limits(1000, 1);
  cf->addLineEntry(iCVC);

  /* -------------------------------------------------------------------------------------------------------------------
   * Number of control volumes in anode
   */
  LE_OneInt *iCVS = new LE_OneInt("Number of CVs in Separator Domain", &initDefaultNumCVsSeparator_);
  iCVS->set_default(10);
  iCVS->set_limits(1000, 1);
  cf->addLineEntry(iCVS);

  /* -------------------------------------------------------------------------------------------------------------------
   * Maximum Number of subcycle steps per global steps
   *     - optional
   *     default = 100;
   */ 
  reqd = 0;
  LE_OneInt *iMAX_SG = new LE_OneInt("Maximum Number of Subcycle Steps per Global Steps", &maxSubgridTimeStepsPerGlobalTimeStep_,
				     reqd, "maxSubgridTimeStepsPerGlobalTimeStep");
  iMAX_SG->set_default(100);
  iMAX_SG->set_limits(10000, 1);
  cf->addLineEntry(iMAX_SG);


  /* ------------------------------------------------------------------------------------------------------------------
   *  Artificial Compressibility
   */
  reqd = 0;
  LE_OneDbl *dac = new LE_OneDbl("Artificial Compressibility", &(artificialCompressibilityInvAtm_), reqd,
                                            "artficialCompressibility");
  dac->set_default(0.0);
  cf->addLineEntry(dac);



  /************************************ DAKOTA INTERFACE ********************************************/
  /* -------------------------------------------------------------------------
   * Use Dakota interface = TRUE
   */
  reqd = 0;
  LE_OneBool *b1 = new LE_OneBool("Use Dakota Interface", &(useDakota_), reqd, "useDakota");
  b1->set_default(false);
  cf->addLineEntry(b1);

  /* ------------------------------------------------------------------------
   * Name of file name to receive parameters from Dakota
   * Dakota Parameters File Name = params.in
   */
  reqd = 0;
  LE_OneStr *s5 = new LE_OneStr("Dakota Parameters File Name",
				 &fromDakotaFileName_, 1, 1, reqd,
				 "fromDakotaFileName" );
  s5->set_default("params.in");
  cf->addLineEntry(s5);

  /* ------------------------------------------------------------------------
   * Name of file name to send results to Dakota
   * Dakota Results File Name = results.out
   */
  reqd = 0;
  LE_OneStr *s6 = new LE_OneStr("Dakota Results File Name",
				 &toDakotaFileName_, 1, 1, reqd,
				 "toDakotaFileName" );
  s6->set_default("results.out");
  cf->addLineEntry(s6);

  BaseEntry::set_SkipUnknownEntries(0);
}
//===================================================================================================================================
/**
 * Do any post processing required.
 * This might include unit conversions, opening files, etc.
 */
void
ProblemStatementCell::post_process_input(BEInput::BlockEntry *cf)
{
    /*
     * Setup Dakota interface if required
     * This parses the file fromDakotaFileName_
     * and then looks for variable names in that which
     * match the member data in this class.
     * Note that you need to go in and manually align variable
     * names in fromDakotaFileName_ with the
     *    if ( di.hasParam("variableName" )
     * statements below.
     */
#ifdef USE_DAKOTA
    if ( useDakota_ ) {
	DakotaInterface di( fromDakotaFileName_, toDakotaFileName_ );

	if ( di.hasParam( "temperature" ) ) {
	    TemperatureReference_ = di.value( "temperature" );
	}
	if ( di.hasParam( "current" ) ) {
	    icurrDischargeSpecified_ = di.value( "current" );
	}
	if ( di.hasParam( "separatorThickness" ) ) {
	    separatorThickness_ = di.value( "separatorThickness" );
	}
    }
#endif
    /*
     * Conversions to SI units
     * NOT NEEDED now that we use unit converters above
     */
    /*
      separatorMass_ *= 1e-3;          // [g] to [kg]
      separatorThickness_ *= 1e-2;     // [cm] to [m]
      separatorArea_ *= 1e-4;          // [cm^2] to [m^2]
      icurrDischargeSpecified_ *= 1e4; // [A/cm^2] to [A/m^2]
    */

 

    // Next check to see if we have area or diameter for each layer
    if (!(separatorArea_ > 0.0)) {
	std::cout << "Warning::ProblemStatementCell() : separator area or diameter not specified" << std::endl;
    } 

    if (anodeCCThickness_ > 0.0) {
	if (anodeBCType_ == 0) {
	    anodeBCType_ = 10;
	}
    }

    if (cathodeCCThickness_ > 0.0 ||  extraCathodeResistance_ > 0.0) {
	if (cathodeBCType_ == 0) {
	    cathodeBCType_ = 10;
	}
    }

    /*
     *  Search for the radial diffusion region block in the block entries
     */
    const BEInput::BlockEntry* be = cf->searchBlockEntry("Extra Phase", false);
  
    const BEInput::BlockEntry* be_cand;
    /*
     *  Collect a set of block entries for the Extra Phase
     */
    std::set<const BlockEntry*> cc = cf->collectBlockEntries("Extra Phase", false);

    std::set<const BlockEntry*>::iterator cc_ptr;

    ExtraPhaseList_.resize(numExtraPhases_);
  
    for (size_t k = 0; k < (size_t) numExtraPhases_; ++k) {
	if (!ExtraPhaseList_[k]) {
	    ExtraPhaseList_[k] = new ExtraPhase();
	}
	ExtraPhase* ep = ExtraPhaseList_[k];
	bool found = false;
	for (cc_ptr = cc.begin(); cc_ptr != cc.end(); cc_ptr++) {
	    be_cand = *cc_ptr;
	    int numT = be_cand->get_NumTimesProcessed();
	    if (numT > 0) {
		be = *cc_ptr;
		int ii = be->multiContribIndex();
		if ((int) k == ii) {
		    found = true;
		    BEInput::LineEntry *le = be->searchLineEntry("Phase Name");
		    BEInput::LE_OneStr* le_str = dynamic_cast<LE_OneStr*>(le);
		    ep->phaseName = le_str->currentTypedValue();

		    le = be->searchLineEntry("Cantera File Name");
		    le_str = dynamic_cast<LE_OneStr*>(le);
		    ep->canteraFileName = le_str->currentTypedValue();

		    le = be->searchLineEntry("Region");
		    le_str = dynamic_cast<LE_OneStr*>(le);
		    ep->regions = le_str->currentTypedValue();
		    string ss = lowercase( ep->regions);

		    le = be->searchLineEntry("Volume Fraction");
		    BEInput::LE_OneDbl*  le_dbl = dynamic_cast<LE_OneDbl*>(le);
		    ep->volFraction = le_dbl->currentTypedValue();

		    std::vector<std::string> v;
		    tokenizeString(ss, v);
		    for (size_t i = 0; i < v.size(); ++i) {
			if (v[i] == "all") {
			    ep->bregionID[0] = 0;
			    ep->bregionID[1] = 1;
			    ep->bregionID[2] = 2;
			} else if (v[i] == "anode") {
			    ep->bregionID[0] = 0;
			} else if (v[i] == "separator") {
			    ep->bregionID[0] = 1;
			} else if (v[i] == "cathode") {
			    ep->bregionID[0] = 2;
			} else {
			    throw m1d_Error("ProblemStatementCell::post_process_input()", "unknown region: " + ep->regions);
			}
		    }
	
		}
	    }
	    if (found) {
		break;
	    }
	}
	if (!found) {
	    throw m1d_Error("post_process", "didn't find the extra phase " + int2str(k));
	}
    }


    /**
     * If we are using time dependent boundary conditions,
     * read in appropriate XML files and generate BoundaryConditions.
     * Note that an alternate means of generating these objects is with an XML node
     */
    if (cathodeBCType_ == 6 || cathodeBCType_ == 7 ) {
	BC_TimeDep_ = new BCsteptable( cathodeBCFile_ );
    }
    if ( cathodeBCType_ == 8 || cathodeBCType_ == 9 ) {
	BC_TimeDep_ = new BClineartable( cathodeBCFile_ );
    }
  
    if ( cathodeBCType_ == 4 ) {
	BC_TimeDep_ = new BCconstant( CathodeVoltageSpecified_, 
				      "Cathode Voltage", "s", "V" ) ;
	BC_TimeDep_->setLowerLimit( startTime_ );
	BC_TimeDep_->setUpperLimit( endTime_ );
    }
    if ( cathodeBCType_ == 5 ) {
	BC_TimeDep_ = new BCconstant( icurrDischargeSpecified_, 
				      "Cathode Current", "s", "Amp/m2" ) ;
	BC_TimeDep_->setLowerLimit( startTime_ );
	BC_TimeDep_->setUpperLimit( endTime_ );
    }
}
//===================================================================================================================================
void ProblemStatementCell::readAnodeInputFile(Electrode_Factory *f )
{
  /**
   * set up anode ELECTRODE_KEY_INPUT based on anodeFile_
   */
  ELECTRODE_KEY_INPUT *ei_tmp = new ELECTRODE_KEY_INPUT();
  /**
   * Initialize a block input structure for the command file
   */
  BEInput::BlockEntry *cfA = new BEInput::BlockEntry("command_file");
  /*
   * Go get the anode description from the anode input file
   */
  // ei_tmp->printLvl_ = 3;
  int retn = ei_tmp->electrode_input(anodeFile_, cfA);
  if (retn == -1) {
    throw  m1d_Error("ProblemStatementCell::post_process_input()",
                     "Electrode failed first creation");
  }
  /*
   * Given the electrodeModelName from first parse, go through and
   * recreate this object.  TODO: create it properly the first time...
   */
  anode_input_ = newElectrodeKeyInputObject(ei_tmp->electrodeModelName, f);
  // anode_input_->printLvl_ = 3;
  /*
   *  Parse the complete child input file
   */
  retn = anode_input_->electrode_input_child(anodeFile_, cfA);
  if (retn == -1) {
    throw  m1d_Error("ProblemStatementCell::post_process_input()",
                     "Electrode input child method failed");
  }
  delete ei_tmp;
  delete cfA;

  if (!(anode_input_->electrodeGrossThickness > 0.0 ) ) {
    std::cout << "Warning::ProblemStatementCell() : anode thickness not specified" << std::endl;
  }
  if (anode_input_->electrodeGrossDiameter > 0.0 ) {
      anode_input_->electrodeGrossArea = 0.25 * Pi
					 * anode_input_->electrodeGrossDiameter
					 * anode_input_->electrodeGrossDiameter;
  }
  if ( !( anode_input_->electrodeGrossArea > 0.0 ) ) {
      std::cout << "Warning::ProblemStatementCell() : anode area or diameter not specified" << std::endl;
  }

}
//===================================================================================================================================
void ProblemStatementCell::readCathodeInputFile(Electrode_Factory *f )
{
  /**
   * set up cathode ELECTRODE_KEY_INPUT based on cathodeFile_
   */
  ELECTRODE_KEY_INPUT *ei_tmp = new ELECTRODE_KEY_INPUT();
  /**
   * Initialize a block input structure for the command file
   */
  BEInput::BlockEntry *cfA = new BEInput::BlockEntry("command_file");
  /*
   * Go get the cathode description from the cathode input file
   */
  // ei_tmp->printLvl_ = 3;
  int retn = ei_tmp->electrode_input(cathodeFile_, cfA);
  if (retn == -1) {
    throw  m1d_Error("ProblemStatementCell::post_process_input()",
                     "Electrode failed first creation");
  }
  /*
   * Given the electrodeModelName from first parse, go through and
   * recreate this object.  TODO: create it properly the first time...
   */
  cathode_input_ = newElectrodeKeyInputObject(ei_tmp->electrodeModelName, f);
  // cathode_input_->printLvl_ = 3;
  /*
   *  Parse the complete child input file
   */
  retn = cathode_input_->electrode_input_child(cathodeFile_, cfA);
  if (retn == -1) {
    throw  m1d_Error("ProblemStatementCell::post_process_input()",
                     "Electrode input child method failed");
  }
  delete ei_tmp;
  delete cfA;

  if (!(cathode_input_->electrodeGrossThickness > 0.0 ) ) {
    std::cout << "Warning::ProblemStatementCell() : cathode thickness not specified" << std::endl;
  }
  if (cathode_input_->electrodeGrossDiameter > 0.0 ) {
      cathode_input_->electrodeGrossArea = 0.25 * Pi * cathode_input_->electrodeGrossDiameter * cathode_input_->electrodeGrossDiameter;
  }
  if ( !( cathode_input_->electrodeGrossArea > 0.0 ) ) {
      std::cout << "Warning::ProblemStatementCell() : cathode area or diameter not specified" << std::endl;
  }

}
//===================================================================================================================================
bool ProblemStatementCell::AnodeCathodeCompatibility()
{
    double electrodeGrossAreaA;
    //see if we know the electrode area
    if (anode_input_->electrodeGrossArea > 0.0) {
        electrodeGrossAreaA = anode_input_->electrodeGrossArea;
    } else if (anode_input_->electrodeGrossDiameter > 0.0) {
        electrodeGrossAreaA = Pi * 0.25 * anode_input_->electrodeGrossDiameter *
                              anode_input_->electrodeGrossDiameter;
    }

   double electrodeGrossAreaC;
    //see if we know the electrode area
    if (cathode_input_->electrodeGrossArea > 0.0) {
        electrodeGrossAreaC = cathode_input_->electrodeGrossArea;
    } else if (anode_input_->electrodeGrossDiameter > 0.0) {
        electrodeGrossAreaC = Pi * 0.25 * cathode_input_->electrodeGrossDiameter *
                              cathode_input_->electrodeGrossDiameter;
    }

    if(! doubleEqual(electrodeGrossAreaA, electrodeGrossAreaC, 1.0E-13, 10)) {
      throw m1d_Error("problemStatementCell::AnodeCathodeCompatibility()",
                      " areas differ: " + ZZCantera::fp2str(electrodeGrossAreaA) + " " + ZZCantera::fp2str(electrodeGrossAreaC));
    } else {
       cathode_input_->electrodeGrossArea = electrodeGrossAreaA;
    } 
                      
    return true;
}
//===================================================================================================================================
/**
 *  This processes the phases in the Cantera input files,
 * fills the PhaseList_ object and other auxiliary data like
 * the numbers and names of species, elements and phases.
 */
void
ProblemStatementCell::InitForInput()
{
  ProblemStatement::InitForInput();

  /* Count number of Cantera files found and make sure we
   * found the number expected
   */
  int ifiles = 0;
  for (; CanteraFileNames[ifiles] != 0; ifiles++) {
  }
  if (ifiles != NumberCanteraFiles) {
    printf("Number of requested files differ\n");
    exit(-1);
  }

  /*
   * Read in all of the phase specifications from the cantera
   * input files into the PhaseList structure.
   */
  ZZCantera::PhaseList *pl = PhaseList_;
  std::string fn;
  bool surNotFound = true;
  for (int i = 0; i < NumberCanteraFiles; i++) {
    fn = CanteraFileNames[i];
    importAllCTMLIntoPhaseList(pl, fn);
    if (surNotFound && (pl->nSurPhases() > 0)) {
      surNotFound = false;
      //   pl->CanteraFNSurface = fn;
    }
  }

  //Count number of phases, species, elements
  nTotPhases_  = pl->nPhases();
  nTotSpecies_ = pl->nSpecies();
  nTotElements_ = pl->nElements();

  SpeciesNames_ = mdpUtil::mdp_alloc_VecFixedStrings(nTotSpecies_, MPEQUIL_MAX_NAME_LEN_P1);
  PhaseNames_ = mdpUtil::mdp_alloc_VecFixedStrings(nTotPhases_, MPEQUIL_MAX_NAME_LEN_P1);
  ElementNames_ = mdpUtil::mdp_alloc_VecFixedStrings(nTotElements_, MPEQUIL_MAX_NAME_LEN_P1);

  /*
   * Generate list of Phases and Species for reference in parsing
   */
  int kT = 0;
  for (int iphase = 0; iphase < nTotPhases_; iphase++) {
    ZZCantera::ThermoPhase *tPhase = &(pl->thermo(iphase));

    std::string id = tPhase->id();
    strncpy(PhaseNames_[iphase], id.c_str(), MPEQUIL_MAX_NAME_LEN);
    int nspecies = tPhase->nSpecies();
    for (int k = 0; k < nspecies; k++) {
      string sname = tPhase->speciesName(k);
      strncpy(SpeciesNames_[kT], sname.c_str(), MPEQUIL_MAX_NAME_LEN);
      kT++;
    }
  }

  const ZZCantera::Elements *eObj = pl->getGlobalElements();

  for (int e = 0; e < nTotElements_; e++) {
    string eName = eObj->elementName(e);
    strncpy(ElementNames_[e], eName.c_str(), MPEQUIL_MAX_NAME_LEN);
  }


}


} //namespace m1d
//===================================================================================================================================
