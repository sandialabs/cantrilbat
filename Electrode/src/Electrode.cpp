/**
 *  @file Electrode.cpp
 *     Definitions of member functions for the base Electrode class, used to model Electrode processes
 *     (see \ref electrode_mgr and class \link Zuzax::Electrode Electrode\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "Electrode.h"
#include "Electrode_Factory.h"
#include "Electrode_FuncCurrent.h"
#include "Electrode_Exception.h"
#include "Electrode_input.h"

#include "cantera/numerics/RootFind.h"
#include "cantera/numerics/solveProb.h"
#include "cantera/kinetics/ExtraGlobalRxn.h"
#include "cantera/base/vec_functions.h"
#include "ApplBase_print.h"


#ifndef SAFE_DELETE
//! Delete a malloced quantity and set the pointer to 0
#define SAFE_DELETE(x)  if ((x)) { delete (x) ; x = nullptr ; }
#endif

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif 
{
//======================================================================================================================
// By default predictor_corrector printing is turned on, at least to the printLvl_ level.
int Electrode::s_printLvl_PREDICTOR_CORRECTOR = 1;

// by default special debug printing is turned off
int Electrode::s_printLvl_DEBUG_SPECIAL = 0;

//==================================================================================================================================
//
//  We read the first character to determine boolean values
//
void Electrode::readEnvironmentalVariables() {
     char *PC_PRINTING = getenv("ELECTRODE_TURN_OFF_PC_PRINTING");
     if (PC_PRINTING) {
        //printf ("TURN_OFF_PC_PRINTING = %s\n", PC_PRINTING);
        if (PC_PRINTING[0] != 'f' && PC_PRINTING[0] != 'F' && PC_PRINTING[0] != '0') {
           Electrode::s_printLvl_PREDICTOR_CORRECTOR=0;
        }
     } 

     char *DS = getenv("ELECTRODE_DEBUG_SPECIAL");
     if (DS) {
        if (DS[0] == 'f' || DS[0] == 'F' || DS[0] == '0') {
           Electrode::s_printLvl_DEBUG_SPECIAL = 0;
        }
        if (DS[0] == 't' || DS[0] == 'T' || DS[0] == '1') {
           Electrode::s_printLvl_DEBUG_SPECIAL = 1;
        }
     }
}
//======================================================================================================================
Electrode::Electrode() :
                PhaseList(),
                electrodeCapacityType_(CAPACITY_ANODE_ECT),
                pendingIntegratedStep_(0),

                prob_type(TP),
                numSurfaces_(0),
                followElectrolyteMoles_(0),
                electrolytePseudoMoles_(1.0),
                externalFieldInterpolationType_(T_FINAL_CONST_FIS),
                t_init_init_(0.0),
                t_final_final_(0.0),
                tinit_(0.0),
                tfinal_(0.0),
                deltaTsubcycleMax_(1.0E300),
                deltaTsubcycle_init_init_(1.0E300),
                deltaTsubcycleNext_(1.0E300),
                deltaTsubcycle_init_next_(1.0E300),
                choiceDeltaTsubcycle_init_(0),
                numIntegrationSubCycles_final_final_(1),
                doThermalPropertyCalculations_(false),
                temperature_(298.15),
                pressure_(1.0E5),
                ElectrodeSolidVolume_(0.0),
                phaseMolarVolumes_(0),
                sphaseMolarAreas_(0),
                VolPM_(0),
                CvPM_(0),
                spMoles_final_(0),
                spMoles_final_final_(0),
                spMoles_init_(0),
                spMoles_init_init_(0),
                spMoles_dot_(0),
                spMoles_predict_(0),
                spMf_init_(0),
                spMf_init_init_(0),
                spMf_final_(0),
                spElectroChemPot_(0),
                phaseVoltages_(0),
                RSD_List_(0),
                phaseMoles_init_(0),
                phaseMoles_init_init_(0),
                phaseMoles_final_(0),
                phaseMoles_final_final_(0),
                phaseMoles_dot_(0),
                numExternalInterfacialSurfaces_(0),
                isExternalSurface_(0),
                surfaceAreaRS_init_(0),
                surfaceAreaRS_final_(0),
                surfaceAreaRS_init_init_(0),
                surfaceAreaRS_final_final_(0),
                spNetProdPerArea_List_(0, 0),
		enthalpyMolar_init_init_(0),
		enthalpyMolar_init_(0),
		enthalpyMolar_final_(0),
		enthalpyMolar_final_final_(0),
		entropyMolar_init_init_(0),
		entropyMolar_init_(0),
		entropyMolar_final_(0),
		entropyMolar_final_final_(0),
		chempotMolar_init_init_(0),
		chempotMolar_init_(0),
		chempotMolar_final_(0),
		chempotMolar_final_final_(0),
                integratedThermalEnergySourceTerm_(0.0),
                integratedThermalEnergySourceTermLast_(0.0),
		integratedThermalEnergySourceTerm_overpotential_(0.0),
		integratedThermalEnergySourceTerm_overpotential_Last_(0.0),
		integratedThermalEnergySourceTerm_reversibleEntropy_(0.0),
		integratedThermalEnergySourceTerm_reversibleEntropy_Last_(0.0),
                electrodeName_("Electrode"),
                numExtraGlobalRxns(0),
                m_egr(0),
                m_rmcEGR(0),
                OCVoverride_ptrList_(0),
                metalPhase_(npos),
                solnPhase_(npos),
                kElectron_(npos),
                deltaVoltage_(0.0),
                electronKmolDischargedToDate_(0.0),
                capacityLeftSpeciesCoeff_(0),
                capacityZeroDoDSpeciesCoeff_(0),
		capacityInitialZeroDod_(0.0),
		depthOfDischargeStarting_(0.0),
                Icurrent_(0.0),
                deltaG_(0),
                inputParticleDiameter_(1.0E-6),
                particleNumberToFollow_(-1.0),
                Radius_exterior_init_init_(5.0E-7),
                Radius_exterior_init_(5.0E-7),
                Radius_exterior_final_(5.0E-7),
                Radius_exterior_final_final_(5.0E-7),
                porosity_(0.0),
                molarAtol_(1.0E-16),
                xmlTimeIncrementData_(nullptr),
                xmlTimeIncrementIntermediateData_(nullptr),
                xmlExternalData_init_init_(nullptr),
                xmlExternalData_init_(nullptr),
                xmlExternalData_final_(nullptr),
                xmlExternalData_final_final_(nullptr),
                xmlStateData_init_init_(nullptr),
                xmlStateData_init_(nullptr),
                xmlStateData_final_(nullptr),
                xmlStateData_final_final_(nullptr),
                eState_save_(nullptr),
                baseNameSoln_("soln"),
                electrodeChemistryModelType_(0),
                electrodeDomainNumber_(0),
                electrodeCellNumber_(0),
                counterNumberIntegrations_(0),
                counterNumberSubIntegrations_(0),
                globalTimeStepNumber_(1),
                writeRestartFileOnSuccessfulStep_(0),
                printLvl_(4),
                printXMLLvl_(0),
                printCSVLvl_(0),
                detailedResidPrintFlag_(0),
                enableExtraPrinting_(false)
{
}
//======================================================================================================================
Electrode::Electrode(const Electrode& right) :
    PhaseList(),
    electrodeCapacityType_(CAPACITY_ANODE_ECT),
    pendingIntegratedStep_(0),
                prob_type(TP),
                numSurfaces_(0),
                followElectrolyteMoles_(0),
                electrolytePseudoMoles_(1.0),
                externalFieldInterpolationType_(T_FINAL_CONST_FIS),
                t_init_init_(0.0),
                t_final_final_(0.0),
                tinit_(0.0),
                tfinal_(0.0),
                deltaTsubcycleMax_(1.0E300),
                deltaTsubcycle_init_init_(1.0E300),
                deltaTsubcycleNext_(1.0E300),
                deltaTsubcycle_init_next_(1.0E300),
                choiceDeltaTsubcycle_init_(0),
                numIntegrationSubCycles_final_final_(1),
                doThermalPropertyCalculations_(false),
                temperature_(298.15),
                pressure_(1.0E5),
                ElectrodeSolidVolume_(0.0),
                phaseMolarVolumes_(0),
                sphaseMolarAreas_(0),
                VolPM_(0),
                CvPM_(0),
                spMoles_final_(0),
                spMoles_final_final_(0),
                spMoles_init_(0),
                spMoles_init_init_(0),
                spMoles_dot_(0),
                spMoles_predict_(0),
                spMf_init_(0),
                spMf_init_init_(0),
                spMf_final_(0),
                spElectroChemPot_(0),
                phaseVoltages_(0),
                RSD_List_(0),
                phaseMoles_init_(0),
                phaseMoles_init_init_(0),
                phaseMoles_final_(0),
                phaseMoles_final_final_(0),
                phaseMoles_dot_(0),
                numExternalInterfacialSurfaces_(0),
                isExternalSurface_(0),
                surfaceAreaRS_init_(0),
                surfaceAreaRS_final_(0),
                surfaceAreaRS_init_init_(0),
                surfaceAreaRS_final_final_(0),
                spNetProdPerArea_List_(0, 0),
		enthalpyMolar_init_init_(0),
		enthalpyMolar_init_(0),
		enthalpyMolar_final_(0),
                enthalpyMolar_final_final_(0),
                entropyMolar_init_init_(0),
                entropyMolar_init_(0),
                entropyMolar_final_(0),
                entropyMolar_final_final_(0),
                chempotMolar_init_init_(0),
                chempotMolar_init_(0),
                chempotMolar_final_(0),
                chempotMolar_final_final_(0),
                integratedThermalEnergySourceTerm_(0.0),
                integratedThermalEnergySourceTermLast_(0.0),
                integratedThermalEnergySourceTerm_overpotential_(0.0),
		integratedThermalEnergySourceTerm_overpotential_Last_(0.0),
		integratedThermalEnergySourceTerm_reversibleEntropy_(0.0),
		integratedThermalEnergySourceTerm_reversibleEntropy_Last_(0.0),
                electrodeName_("Electrode"),
                numExtraGlobalRxns(0),
                m_egr(0),
                m_rmcEGR(0),
                OCVoverride_ptrList_(0),
                metalPhase_(npos),
                solnPhase_(npos),
                kElectron_(npos),
                deltaVoltage_(0.0),
                electronKmolDischargedToDate_(0.0),
                capacityLeftSpeciesCoeff_(0),
                capacityZeroDoDSpeciesCoeff_(0),
                capacityInitialZeroDod_(0.0),
		depthOfDischargeStarting_(0.0),
                Icurrent_(0.0),
                deltaG_(0),
                inputParticleDiameter_(1.0E-6),
                particleNumberToFollow_(-1.0),
                Radius_exterior_init_init_(5.0E-7),
                Radius_exterior_init_(5.0E-7),
                Radius_exterior_final_(5.0E-7),
                Radius_exterior_final_final_(5.0E-7),
                porosity_(0.0),
                molarAtol_(1.0E-16),
                xmlTimeIncrementData_(nullptr),
                xmlTimeIncrementIntermediateData_(nullptr),
                xmlExternalData_init_init_(0),
                xmlExternalData_init_(0),
                xmlExternalData_final_(0),
                xmlExternalData_final_final_(0),
                xmlStateData_init_init_(0),
                xmlStateData_init_(0),
                xmlStateData_final_(0),
                xmlStateData_final_final_(0),
                eState_save_(nullptr),
                baseNameSoln_("soln"),
                electrodeChemistryModelType_(0),
                electrodeDomainNumber_(0),
                electrodeCellNumber_(0),
                counterNumberIntegrations_(0),
                counterNumberSubIntegrations_(0),
                globalTimeStepNumber_(1),
                writeRestartFileOnSuccessfulStep_(0),
                printLvl_(4),
                printXMLLvl_(0),
                printCSVLvl_(0),
                detailedResidPrintFlag_(0),
                enableExtraPrinting_(false)
{
    operator=(right);
}
//======================================================================================================================
Electrode& Electrode::operator=(const Electrode& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    PhaseList::operator=(right);
    /*
     * We do a straight assignment operator on all of the
     * data. The vectors are copied.
     */
    electrodeCapacityType_ = right.electrodeCapacityType_;
    pendingIntegratedStep_ = right.pendingIntegratedStep_;
    prob_type = right.prob_type;
    numSurfaces_ = right.numSurfaces_;
    followElectrolyteMoles_ = right.followElectrolyteMoles_;
    electrolytePseudoMoles_ = right.electrolytePseudoMoles_;
    externalFieldInterpolationType_ = right.externalFieldInterpolationType_;
    t_init_init_ = right.t_init_init_;
    t_final_final_ = right.t_final_final_;

    /*
     * Copy over the ReactingSurDomain list
     * The electrode object owns the ReactingSurDomain. 
     */
    for (size_t i = 0; i < RSD_List_.size(); ++i) {
        delete RSD_List_[i];
    }
    RSD_List_ = right.RSD_List_;  
   
    for (size_t i = 0; i < right.numSurfaces_; i++) {
        bool idHit = false;
        //if (right.m_rSurDomain == right.RSD_List_[i]) {
        //	idHit = true;
        //}
        if (right.RSD_List_[i]) {
            //
            //  Call the copy constructor for the copying operation.
            //
            RSD_List_[i] = new ReactingSurDomain(*(right.RSD_List_[i]));
            //
            // NOTE: we must fix up shallow pointers in the new ReactingSurDomain object
            //       so that everything points into this object.
	    //
	    RSD_List_[i]->m_pl = this;
	    std::vector<ThermoPhase*> tpList(0);
	    for (size_t iph = 0; iph < RSD_List_[i]->nPhases(); ++iph) {
		std::string ss = RSD_List_[i]->tpList_IDs_[iph];
		bool notFound = true;
		for (size_t jph = 0; jph < nPhases(); jph++) {
		    ThermoPhase *tp = & thermo(jph);
		    if (tp->id() == ss) {
			notFound = false;
			tpList.push_back(tp);
			break;
		    }
		}
		if (notFound) {
		    throw Electrode_Error("Electrode assignment operator", "could not find id for phase " + ss);
		}
	    }
	    RSD_List_[i]->assignShallowPointers(tpList);

            if (RSD_List_[i]->OCVmodel_) {
              RSD_OCVmodel* mm = RSD_List_[i]->OCVmodel_;
              std::string ss = (mm->solidPhasePtr())->id();
	      bool notFound = true;
              for (size_t jph = 0; jph <  nPhases(); jph++)  {
                    ThermoPhase *tp = & thermo(jph);
                    if (tp->id() == ss) {
                        notFound = false;
                        mm->assignShallowPointers(tp);
                        break;
                    }
                }
                if (notFound) {
                    throw Electrode_Error("Electrode assignment operator", "could not find id for phase " + ss);
                }

            }

            if (idHit) {
                // m_rSurDomain = RSD_List_[i];
            }
        }
    }
    tinit_ = right.tinit_;
    tfinal_ = right.tfinal_;
    deltaTsubcycleMax_ = right.deltaTsubcycleMax_;
    deltaTsubcycle_init_init_ = right.deltaTsubcycle_init_init_;
    deltaTsubcycle_ = right.deltaTsubcycle_;
    deltaTsubcycleNext_ = right.deltaTsubcycleNext_;
    deltaTsubcycle_init_next_ = right.deltaTsubcycle_init_next_;
    choiceDeltaTsubcycle_init_ = right.choiceDeltaTsubcycle_init_;
    numIntegrationSubCycles_final_final_ = right.numIntegrationSubCycles_final_final_;
    doThermalPropertyCalculations_ = right.doThermalPropertyCalculations_;
    temperature_ = right.temperature_;
    pressure_ = right.pressure_;
    ElectrodeSolidVolume_ = right.ElectrodeSolidVolume_;
    phaseMolarVolumes_ = right.phaseMolarVolumes_;
    sphaseMolarAreas_ = right.sphaseMolarAreas_;
    VolPM_ = right.VolPM_;
    CvPM_ = right.CvPM_;
    spMoles_final_ = right.spMoles_final_;
    spMoles_final_final_ = right.spMoles_final_final_;
    spMoles_init_ = right.spMoles_init_;
    spMoles_init_init_ = right.spMoles_init_init_;
    spMoles_dot_ = right.spMoles_dot_;
    spMoles_predict_ = right.spMoles_predict_;
    spMf_init_ = right.spMf_init_;
    spMf_init_init_ = right.spMf_init_init_;
    spMf_final_ = right.spMf_final_;
    spMf_final_final_ = right.spMf_final_final_;
    spMf_final_dot = right.spMf_final_dot;
    spElectroChemPot_ = right.spElectroChemPot_;
    phaseVoltages_ = right.phaseVoltages_;

    //RSD_LIST_

    numRxns_ = right.numRxns_;
    ActiveKineticsSurf_ = right.ActiveKineticsSurf_;
    phaseMoles_init_ = right.phaseMoles_init_;
    phaseMoles_init_init_ = right.phaseMoles_init_init_;
    phaseMoles_final_ = right.phaseMoles_final_;
    phaseMoles_final_final_ = right.phaseMoles_final_final_;
    phaseMoles_dot_ = right.phaseMoles_dot_;
    justBornPhase_ = right.justBornPhase_;
    justDiedPhase_ = right.justDiedPhase_;
    numExternalInterfacialSurfaces_ = right.numExternalInterfacialSurfaces_;
    isExternalSurface_ = right.isExternalSurface_;
    surfaceAreaRS_init_ = right.surfaceAreaRS_init_;
    surfaceAreaRS_final_ = right.surfaceAreaRS_final_;
    surfaceAreaRS_init_init_ = right.surfaceAreaRS_init_init_;
    surfaceAreaRS_final_final_ = right.surfaceAreaRS_final_final_;
    spNetProdPerArea_List_ = right.spNetProdPerArea_List_;
    spMoleIntegratedSourceTerm_ = right.spMoleIntegratedSourceTerm_;
    spMoleIntegratedSourceTermLast_ = right.spMoleIntegratedSourceTermLast_;
    enthalpyMolar_init_init_ = right.enthalpyMolar_init_init_;
    enthalpyMolar_init_ = right.enthalpyMolar_init_;
    enthalpyMolar_final_ = right.enthalpyMolar_final_;
    enthalpyMolar_final_final_ = right.enthalpyMolar_final_final_;
    entropyMolar_init_init_ = right.entropyMolar_init_init_;
    entropyMolar_init_ = right.entropyMolar_init_;
    entropyMolar_final_ = right.entropyMolar_final_;
    entropyMolar_final_final_ = right.entropyMolar_final_final_;
    chempotMolar_init_init_ = right.chempotMolar_init_init_;
    chempotMolar_init_ = right.chempotMolar_init_;
    chempotMolar_final_ = right.chempotMolar_final_;
    chempotMolar_final_final_ = right.chempotMolar_final_final_;
    integratedThermalEnergySourceTerm_ = right.integratedThermalEnergySourceTerm_;
    integratedThermalEnergySourceTermLast_ = right.integratedThermalEnergySourceTermLast_;
    integratedThermalEnergySourceTerm_overpotential_ = right.integratedThermalEnergySourceTerm_overpotential_;
    integratedThermalEnergySourceTerm_overpotential_Last_ = right.integratedThermalEnergySourceTerm_overpotential_Last_;
    integratedThermalEnergySourceTerm_reversibleEntropy_ = right.integratedThermalEnergySourceTerm_reversibleEntropy_;
    integratedThermalEnergySourceTerm_reversibleEntropy_Last_ = right.integratedThermalEnergySourceTerm_reversibleEntropy_Last_;
    electrodeName_ = right.electrodeName_;
    numExtraGlobalRxns = right.numExtraGlobalRxns;

    ZZCantera::deepStdVectorPointerCopy<ExtraGlobalRxn>(right.m_egr, m_egr);
    ZZCantera::deepStdVectorPointerCopy<RxnMolChange>(right.m_rmcEGR, m_rmcEGR);
    ZZCantera::deepStdVectorPointerCopy<OCV_Override_input>(right.OCVoverride_ptrList_, OCVoverride_ptrList_);

    metalPhase_ = right.metalPhase_;
    solnPhase_ = right.solnPhase_;
    kElectron_ = right.kElectron_;
    kKinSpecElectron_sph_ = right.kKinSpecElectron_sph_;
    deltaVoltage_ = right.deltaVoltage_;
    electronKmolDischargedToDate_ = right.electronKmolDischargedToDate_;
    capacityLeftSpeciesCoeff_ = right.capacityLeftSpeciesCoeff_;
    capacityZeroDoDSpeciesCoeff_ = right.capacityZeroDoDSpeciesCoeff_;
    capacityInitialZeroDod_ = right.capacityInitialZeroDod_;
    depthOfDischargeStarting_ = right.depthOfDischargeStarting_;
    Icurrent_ = right.Icurrent_;
    deltaG_ = right.deltaG_;
    inputParticleDiameter_ = right.inputParticleDiameter_;
    particleNumberToFollow_ = right.particleNumberToFollow_;
    Radius_exterior_init_init_ = right.Radius_exterior_init_init_;
    Radius_exterior_init_ = right.Radius_exterior_init_;
    Radius_exterior_final_ = right.Radius_exterior_final_;
    Radius_exterior_final_final_ = right.Radius_exterior_final_final_;
    porosity_ = right.porosity_;
    molarAtol_ = right.molarAtol_;

    SAFE_DELETE(xmlTimeIncrementData_);
    if (right.xmlTimeIncrementData_) {
        xmlTimeIncrementData_ = new XML_Node(*right.xmlTimeIncrementData_);
   }

    SAFE_DELETE(xmlTimeIncrementIntermediateData_);
    if (right.xmlTimeIncrementIntermediateData_) {
        xmlTimeIncrementIntermediateData_ = new XML_Node(*right.xmlTimeIncrementIntermediateData_);
    }

    SAFE_DELETE(xmlExternalData_init_init_);
    if (right.xmlExternalData_init_init_) {
        xmlExternalData_init_init_ = new XML_Node(*right.xmlExternalData_init_init_);
    }

    SAFE_DELETE(xmlExternalData_init_);
    if (right.xmlExternalData_init_) {
        xmlExternalData_init_ = new XML_Node(*right.xmlExternalData_init_);
    }

    SAFE_DELETE(xmlExternalData_final_);
    if (right.xmlExternalData_final_) {
        xmlExternalData_final_ = new XML_Node(*right.xmlExternalData_final_);
    }

    SAFE_DELETE(xmlExternalData_final_final_);
    if (right.xmlExternalData_final_final_) {
        xmlExternalData_final_final_ = new XML_Node(*right.xmlExternalData_final_final_);
    }

    SAFE_DELETE(xmlStateData_init_init_);
    if (right.xmlStateData_init_init_) {
        xmlStateData_init_init_ = new XML_Node(*right.xmlStateData_init_init_);
    }

    SAFE_DELETE(xmlStateData_init_);
    if (right.xmlStateData_init_) {
        xmlStateData_init_ = new XML_Node(*right.xmlStateData_init_);
    }

    SAFE_DELETE(xmlStateData_final_);
    if (right.xmlStateData_final_) {
        xmlStateData_final_ = new XML_Node(*right.xmlStateData_final_);
    }

    SAFE_DELETE(xmlStateData_final_final_);
    if (right.xmlStateData_final_final_) {
        xmlStateData_final_final_ = new XML_Node(*right.xmlStateData_final_final_);
    }

    SAFE_DELETE(eState_save_);
    if (right.eState_save_) {
        eState_save_ = right.eState_save_->duplMyselfAsEState(this);
    }

    baseNameSoln_ = right.baseNameSoln_;

    electrodeChemistryModelType_ = right.electrodeChemistryModelType_;
    electrodeDomainNumber_ = right.electrodeDomainNumber_;
    electrodeCellNumber_ = right.electrodeCellNumber_;
    counterNumberIntegrations_ = 0;
    counterNumberSubIntegrations_ = 0;
    globalTimeStepNumber_ = right.globalTimeStepNumber_;
    writeRestartFileOnSuccessfulStep_  = right.writeRestartFileOnSuccessfulStep_;
    printLvl_ = right.printLvl_;
    printXMLLvl_ = right.printXMLLvl_;
    printCSVLvl_ = right.printCSVLvl_;
    detailedResidPrintFlag_ = right.detailedResidPrintFlag_;
    enableExtraPrinting_ = right.enableExtraPrinting_;

    return *this;
}
//======================================================================================================================
Electrode::~Electrode()
{
    if (RSD_List_.size() > 0) {
    for (size_t i = 0; i < (size_t) numSurfaces_; i++) {
       
        if (RSD_List_[i]) {
            delete RSD_List_[i];
        }
    }
    }
    //m_rSurDomain = 0;

    size_t is = m_egr.size();
    for (size_t i = 0; i < is; i++) {
        if (m_egr[i]) {
            delete m_egr[i];
            m_egr[i] = 0;
        }
    }
    is = m_rmcEGR.size();
    for (size_t i = 0; i < is; i++) {
        if (m_rmcEGR[i]) {
            delete m_rmcEGR[i];
            m_rmcEGR[i] = 0;
        }
    }
    for (size_t i = 0; i < OCVoverride_ptrList_.size(); ++i) {
         delete OCVoverride_ptrList_[i];
    } 

    SAFE_DELETE(xmlTimeIncrementData_);
    SAFE_DELETE(xmlTimeIncrementIntermediateData_);
    SAFE_DELETE(xmlExternalData_init_init_);
    SAFE_DELETE(xmlExternalData_init_);
    SAFE_DELETE(xmlExternalData_final_);
    SAFE_DELETE(xmlExternalData_final_final_);

    SAFE_DELETE(xmlStateData_init_init_);
    SAFE_DELETE(xmlStateData_init_);
    SAFE_DELETE(xmlStateData_final_);
    SAFE_DELETE(xmlStateData_final_final_);

    SAFE_DELETE(eState_save_);
}
//======================================================================================================================
// Duplicator function
/*
 *  Duplicate the current electrode object, returning a base electrode pointer
 *
 *  (Virtual function from Electrode.h)
 */
Electrode* Electrode::duplMyselfAsElectrode() const
{
    return 0;
}
//======================================================================================================================
// Set the electrode ID information
void Electrode::setID(int domainNum, int cellNum)
{
    electrodeDomainNumber_ = domainNum;
    electrodeCellNumber_ = cellNum;
}
//======================================================================================================================
Electrode_Types_Enum Electrode::electrodeType() const
{
    return BASE_TYPE_ET;
}
//======================================================================================================================
// Add additional Keylines for child electrode objects, and then read them in
/*
 *   (virtual function from Electrode)
 *   (overload virtual function - replaces parent members)
 *
 *   This function will replace the ELECTRODE_KEY_INPUT structure with an expanded
 *   child member structure containing the extra information.
 *
 *   If the command file has been read before, it will then reparse the command file
 *   storring the new information in the ELECTRODE_KEY_INPUT structure.
 *
 *    @param ei  Handle to the ELECTRODE_KEY_INPUT base pointer. This handle may change
 *               as the child class of  ELECTRODE_KEY_INPUT gets malloced.
 *
 *    @return  0 successful but no change in ei
 *             1 Successful and ei has changed
 *            -1 unsuccessful fatal error of some kind.
 */
int Electrode::electrode_input_child(ELECTRODE_KEY_INPUT** ei)
{
    return 0;
}
//==================================================================================================================================
//! Throw an error
/*!
 *  @param[in]               pos                 Position within the phaselist names
 *  @param[in]               actual              actual phase name found
 *  @param[in]               expected            expected phase name
 */
static void ErrorModelType(int pos, std::string actual, std::string expected)
{
    throw Electrode_Error("Electrode::electrode_model_create() model id",
                       "At pos " + int2str(pos) + ", expected phase " + expected + " but got phase " + actual);
}
//==================================================================================================================================
//! check match against expected Phase ID's
/*!
 *  Throw an error if the expected phase name isn't where it's expected.
 *  We are hard-coding expected phase names in the xml files here.
 *
 *  @param[in]               pn                  list of phase names
 *  @param[in]               pos                 position in the vector
 *  @param[in]               expected            Expected phase name
 */
static void pmatch(std::vector<std::string>& pn, int pos, std::string expected)
{
    if (pos <= ((int) pn.size() - 1)) {
        if (pn[pos] != expected) {
            ErrorModelType(pos, pn[pos], expected);
        }
    }
}
//==================================================================================================================================
int Electrode::electrode_model_create(ELECTRODE_KEY_INPUT* ei)
{
    Electrode_Types_Enum ieos = electrodeType();
    Electrode_Types_Enum seos = string_to_Electrode_Types_Enum(ei->electrodeModelName);

    if (ieos != seos) {
        throw Electrode_Error(
                " Electrode::electrode_model_create()",
                "Electrode Object Type, " + Electrode_Types_Enum_to_string(ieos) + ", is different than requested type, "
                        + ei->electrodeModelName);
    }
    /*
     *  Read environmental variables and initialize statics
     */
    Electrode::readEnvironmentalVariables();

    // use the assignment operator to transfer for now
    // May get more sophisticated
    PhaseList::operator=(*(ei->m_pl));

    // resize m_NumTotPhases = total number of phase vectors
    phaseMoles_init_.resize(m_NumTotPhases, 0.0);
    phaseMoles_init_init_.resize(m_NumTotPhases, 0.0);
    phaseMoles_final_.resize(m_NumTotPhases, 0.0);
    phaseMoles_final_final_.resize(m_NumTotPhases, 0.0);
    phaseMoles_dot_.resize(m_NumTotPhases, 0.0);
    phaseVoltages_.resize(m_NumTotPhases, 0.0);
    phaseMolarVolumes_.resize(m_NumTotPhases, 0.0);
    justBornPhase_.resize(m_NumTotPhases, 0);
    justDiedPhase_.resize(m_NumTotPhases, 0);

    // resize volume phase vectors

    // Identify the number of surfaces with the number of surface phases for now
    numSurfaces_ = m_NumSurPhases;
    // resize surface phase vectors
    isExternalSurface_.resize(numSurfaces_, true);
    surfaceAreaRS_init_.resize(numSurfaces_, 0.0);
    surfaceAreaRS_final_.resize(numSurfaces_, 0.0);
    surfaceAreaRS_init_init_.resize(numSurfaces_, 0.0);
    surfaceAreaRS_final_final_.resize(numSurfaces_, 0.0);
    RSD_List_.resize(numSurfaces_, 0);
    numRxns_.resize(numSurfaces_, static_cast<size_t>(0));
    ActiveKineticsSurf_.resize(numSurfaces_, 0);
    sphaseMolarAreas_.resize(numSurfaces_, 0.0);

    // Resize species vectors
    spMoles_init_.resize(m_NumTotSpecies, 0.0);
    spMoles_final_.resize(m_NumTotSpecies, 0.0);
    spMoles_final_final_.resize(m_NumTotSpecies, 0.0);
    spMoles_init_init_.resize(m_NumTotSpecies, 0.0);
    spMoles_dot_.resize(m_NumTotSpecies, 0.0);
    spMoles_predict_.resize(m_NumTotSpecies, 0.0);
    spMf_final_.resize(m_NumTotSpecies, 0.0);
    spMf_init_.resize(m_NumTotSpecies, 0.0);
    spMf_init_init_.resize(m_NumTotSpecies, 0.0);
    spMf_final_final_.resize(m_NumTotSpecies, 0.0);

    enthalpyMolar_init_init_.resize(m_NumTotSpecies, 0.0);
    enthalpyMolar_init_.resize(m_NumTotSpecies, 0.0);
    enthalpyMolar_final_.resize(m_NumTotSpecies, 0.0);
    enthalpyMolar_final_final_.resize(m_NumTotSpecies, 0.0);
    entropyMolar_init_init_.resize(m_NumTotSpecies, 0.0);
    entropyMolar_init_.resize(m_NumTotSpecies, 0.0);
    entropyMolar_final_.resize(m_NumTotSpecies, 0.0);
    entropyMolar_final_final_.resize(m_NumTotSpecies, 0.0);
    chempotMolar_init_init_.resize(m_NumTotSpecies, 0.0);
    chempotMolar_init_.resize(m_NumTotSpecies, 0.0);
    chempotMolar_final_.resize(m_NumTotSpecies, 0.0);
    chempotMolar_final_final_.resize(m_NumTotSpecies, 0.0);

    VolPM_.resize(m_NumTotSpecies, 0.0);
    CvPM_.resize(m_NumTotSpecies, 0.0);
    spElectroChemPot_.resize(m_NumTotSpecies, 0.0);

    capacityLeftSpeciesCoeff_.resize(m_NumTotSpecies, 0.0);
    capacityZeroDoDSpeciesCoeff_.resize(m_NumTotSpecies, 0.0);

    capacityLeftSpeciesCoeffPlat_.resize(5, m_NumTotSpecies);
    capacityZeroDoDSpeciesCoeffPlat_.resize(5, m_NumTotSpecies);

    /*
     * OK, Find the first kinetics object
     */
    size_t isphfound = npos;
    for (size_t isph = 0; isph < m_NumSurPhases; isph++) {
        if (SurPhaseHasKinetics_[isph]) {
            isphfound = isph;
            break;
        }
    }
    if (isphfound == npos) {
        for (size_t i = 0; i < m_NumVolPhases; i++) {
            if (VolPhaseHasKinetics_[i]) {
                break;
            }
        }
    }

    /*
     * Assign a reacting surface to each surface in the PhaseList object
     */
    for (size_t isurf = 0; isurf < numSurfaces_; isurf++) {
        // Right now there is a correspondence between isph and isurf for the initial read in.
        // However this may break in the future. Therefore, I've put in this logic for future expansion
        int isph = isurf;
        if (SurPhaseHasKinetics_[isph]) {
            ReactingSurDomain* rsd = new ReactingSurDomain();
            int ok = rsd->importFromPL(this, isph);
            if (!ok) {
                throw Electrode_Error("cttables main:", "rSurDomain returned an error");
            }
            //
            // Check to see if we have entered an OCVoverride for this species. If we have then
            // modify the reacting surface
            //
	    ZZCantera::OCV_Override_input *ocv_ptr = ei->OCVoverride_ptrList[isurf];
            if (ocv_ptr->numTimes > 0) {

               rsd->addOCVoverride(ocv_ptr);
            }

            // if  OCVoverride_ptrList
            // We carry a list of pointers to ReactingSurDomain
            RSD_List_[isurf] = rsd;
            numRxns_[isurf] = rsd->nReactions();
            ActiveKineticsSurf_[isurf] = 1;
        }
    }
    ZZCantera::deepStdVectorPointerCopy<OCV_Override_input>(ei->OCVoverride_ptrList, OCVoverride_ptrList_);
    //OCVoverride_ptrList_ = ei->OCVoverride_ptrList;

    /*
     * Calculate the number of external interfacial surfaces
     */
    numExternalInterfacialSurfaces_ = 1;

    /*
     * Resize the species production vector for the electrode
     *  -> per area and total net
     */
    spMoleIntegratedSourceTerm_.resize(m_NumTotSpecies, 0.0);
    spMoleIntegratedSourceTermLast_.resize(m_NumTotSpecies, 0.0);
    spNetProdPerArea_List_.resize(m_NumTotSpecies, numSurfaces_, 0.0);
    integratedThermalEnergySourceTerm_ = 0.0;
    integratedThermalEnergySourceTermLast_ = 0.0;
    integratedThermalEnergySourceTerm_overpotential_ = 0.0;
    integratedThermalEnergySourceTerm_overpotential_Last_ = 0.0;
    integratedThermalEnergySourceTerm_reversibleEntropy_ = 0.0;
    integratedThermalEnergySourceTerm_reversibleEntropy_Last_ = 0.0;

    /*
     * Load up the temperature and pressure
     */
    temperature_ = ei->Temperature;
    pressure_ = ei->Pressure;
    
    /*
     *  Loop Over all phases in the PhaseList, adding these
     *  formally to the Electrode object.
     */
    size_t nspecies = 0;
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
	// Get the pointer to the ThermoPhase object
        ThermoPhase* tphase = &(thermo(iph));
        size_t nSpecies = tphase->nSpecies();

        // Find the name of the input block
        std::string phaseBath = "Bath Specification for Phase ";
        std::string phaseNm = tphase->name();
        phaseBath += phaseNm;

        /*
         * Search the ReactingSurDomain to see if the current phase is in
         * the object
         */
        std::string pname = tphase->id();

        size_t kstart = nspecies;
        nspecies += nSpecies;

        /*
         *  We copy the mole numbers and mole fraction information from the bath gas section of the input file
         *  here. Note, the total quantitites of electrode material and electrolyte material may be overwritten
         *  later.
	 *  We only really need to populate the _final_ state. We will set all states equal to each other at the end 
	 *  of this process.
         */
        double sum = 0.0;
        for (size_t k = 0; k < nSpecies; k++) {
            spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k] = ei->MoleNumber[m_PhaseSpeciesStartIndex[iph] + k];
            spMoles_init_[m_PhaseSpeciesStartIndex[iph] + k] = spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k];
            sum += ei->MoleNumber[m_PhaseSpeciesStartIndex[iph] + k];
        }

        if (sum > 1.0E-300) {
            for (size_t k = 0; k < nSpecies; k++) {
                spMf_final_[kstart + k] = spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k] / sum;
                spMf_init_[kstart + k] = spMf_final_[kstart + k];
                spMf_init_init_[kstart + k] = spMf_final_[kstart + k];
                spMf_final_final_[kstart + k] = spMf_final_[kstart + k];
            }
	    /*
	     *  If the phase formally exists, we set the mole fraction vector from the usual calculation
	     *  and post it down to the ThermoPhase object
	     */
            tphase->setMoleFractions(&(spMf_final_[kstart]));
        } else {
	    /*
	     *  If the phase doesn't formally exist, we get the mole fraction vector from the
	     *  ThermoPhase object itself. This is important as the ThermoPhase object has an
	     *  idea of what a reasonable initial state would be, while we do not have any idea.
	     */
            tphase->getMoleFractions(&(spMf_final_[kstart]));
            tphase->getMoleFractions(&(spMf_init_[kstart]));
        }
	/*
	 *  set the voltage from the thermophase object 
	 */
        phaseVoltages_[iph] = tphase->electricPotential();
    }

    /* set the absolute tolerance to 1e-5 the total number of moles */
    // HEWSON -- 12/1/10 -- this particular tolerance value
    // and all magic tolerance values throughout are subject to revision
    double tMoles = 0.0;
    for (size_t k = 0; k < m_NumTotSpecies; k++) {
        tMoles += spMoles_final_[k];
    }
    molarAtol_ = tMoles * 1.0E-5;

    /*
     * Now that spMoles_final_ is sized appropriately, we can call updatePhaseNumbers()
     */
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        Electrode::updateState_Phase(iph);
    }

    electrodeName_ = ei->electrodeName;

    if (ei->electrodeCapacityType == 0) {
        electrodeCapacityType_ = CAPACITY_ANODE_ECT;
    } else {
        electrodeCapacityType_ = CAPACITY_CATHODE_ECT;
    }

    //
    // Copy over the input parameters for extra global reactions
    //
    numExtraGlobalRxns = ei->numExtraGlobalRxns;

    for (size_t iph = 0; iph < m_NumVolPhases; iph++) {
        phaseMoles_init_[iph] = phaseMoles_final_[iph];
        phaseMoles_init_init_[iph] = phaseMoles_final_[iph];
    }

    /*
     * Identify the phases which will ID the voltage drops
     *  metalPhase_ = Phase where the electron exists.
     */
    metalPhase_ = npos;
    solnPhase_ = npos;
    kElectron_ = npos;
    //   InterfaceKinetics *iK = m_rSurDomain;
    ReactingSurDomain* rsd = RSD_List_[0];
    std::vector<RxnMolChange*>& rmcV = rsd->rmcVector;
    size_t nr = rmcV.size();
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        ThermoPhase* tp = &(thermo(iph));
        size_t nSpecies = tp->nSpecies();
        size_t nElements = tp->nElements();
        size_t eElectron = tp->elementIndex("E");
        if (eElectron != npos) {
            for (size_t k = 0; k < nSpecies; k++) {
                if (tp->nAtoms(k, eElectron) == 1) {
                    int ifound = 1;
                    for (size_t e = 0; e < nElements; e++) {
                        if (tp->nAtoms(k, e) != 0.0) {
                            if (e != eElectron) {
                                ifound = 0;
                            }
                        }
                    }
                    if (ifound == 1) {
                        metalPhase_ = iph;
                        kElectron_ = m_PhaseSpeciesStartIndex[iph] + k;
                    }
                }
            }
        }
        if (iph != metalPhase_) {
            size_t jph = rsd->PLtoKinPhaseIndex_[iph];
            if (jph != npos) {
                for (size_t i = 0; i < nr; i++) {
                    RxnMolChange* rmc = rmcV[i];
                    if (rmc->m_phaseChargeChange[jph] != 0) {
                        if (rmc->m_phaseDims[jph] == 3) {
                            solnPhase_ = iph;
                            break;
                        }
                    }
                }
            }
        }
    }
    if (solnPhase_ == npos) {
        throw Electrode_Error("Electrode::electrode_model_create()", "Couldn't find soln phase");
    }
    if (metalPhase_ == npos) {
        throw Electrode_Error("Electrode::electrode_model_create()", "Couldn't find metal phase");
    }

    kKinSpecElectron_sph_.resize(numSurfaces_, npos);
    for (size_t isurf = 0; isurf < numSurfaces_; isurf++) {
        kKinSpecElectron_sph_[isurf] = npos;
        ReactingSurDomain* rsd = RSD_List_[isurf];
        if (rsd) {
            std::vector<RxnMolChange*>& rmcV = rsd->rmcVector;
            if (rsd) {
                InterfaceKinetics* iK = rsd;
                for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
                    ThermoPhase* tp = &(thermo(iph));
                    size_t nSpecies = tp->nSpecies();
                    size_t nElements = tp->nElements();
                    size_t eElectron = tp->elementIndex("E");
                    if (eElectron != npos) {
                        for (size_t k = 0; k < nSpecies; k++) {
                            if (tp->nAtoms(k, eElectron) == 1) {
                                int ifound = 1;
                                for (size_t e = 0; e < nElements; e++) {
                                    if (tp->nAtoms(k, e) != 0.0) {
                                        if (e != eElectron) {
                                            ifound = 0;
                                        }
                                    }
                                }
                                if (ifound == 1) {
                                    metalPhase_ = iph;
                                    kElectron_ = m_PhaseSpeciesStartIndex[iph] + k;
                                    kKinSpecElectron_sph_[isurf] = iK->kineticsSpeciesIndex(k, iph);
                                }
                            }
                        }
                    }
                    if (iph != metalPhase_) {
                        size_t jph = rsd->PLtoKinPhaseIndex_[iph];
                        if (jph != npos) {
                            for (size_t i = 0; i < nr; i++) {
                                RxnMolChange* rmc = rmcV[i];
                                if (rmc->m_phaseChargeChange[jph] != 0) {
                                    if (rmc->m_phaseDims[jph] == 3) {
                                        solnPhase_ = iph;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Identify the electrode voltage and calculate it
    deltaVoltage_ = phaseVoltages_[metalPhase_] - phaseVoltages_[solnPhase_];

    /*
     *  Find the maximum number of reactions in any reacting domain and the max of the number of species in the 
     *  Electrode object.
     */
    size_t mR = spMoles_final_.size();
    for (size_t isk = 0; isk <  RSD_List_.size(); isk++) {
        ReactingSurDomain* rsd = RSD_List_[isk];
        if (rsd) {
            size_t nr = rsd->nReactions();
            if (mR < nr) {
                mR = nr;
            }
        }
    }
    deltaG_.resize(mR);

    //
    // Process Extra global reactions
    //
    processExtraGlobalRxnPathways(ei->m_EGRList);

    /*
     *  Set up the particleDiameter_ and particleNumberToFollow_ fields.
     *   If any geometry information is supplied in the input file, then the geometry is
     *   used to rescale the input mole numbers of species. If no geometry information
     *   is input, the the partcleNumbertoFollow_ variable is calculated to be consistent
     *   with the input mole numbers.
     */
    // particle diameter is a required parameter
    if (ei->particleDiameter <= 0.0) {
        throw Electrode_Error("Electrode::electrode_model_create()", "Particle diameter is <= 0");
    }
    inputParticleDiameter_ = ei->particleDiameter;
    Radius_exterior_final_ = inputParticleDiameter_ / 2.0;
    Radius_exterior_init_ = Radius_exterior_final_;
    Radius_exterior_init_init_ = Radius_exterior_final_;
    Radius_exterior_final_final_ = Radius_exterior_final_;
    double partVol = Pi * inputParticleDiameter_ * inputParticleDiameter_ * inputParticleDiameter_ / 6.0;
    bool rewriteMoleNumbers = false;
    /*
     *  Gross volume is the total superficial volume of the electrode. It is the sum of the
     *  solid volume and the liquid electrylte volume.
     *    - Compute the gross volume if this information has been supplied through the input file
     */
    double grossVolume = 0.0;
    if (ei->electrodeGrossArea > 0.0 && ei->electrodeGrossThickness > 0.0) {
        grossVolume = ei->electrodeGrossArea * ei->electrodeGrossThickness;
    } else if (ei->electrodeGrossDiameter > 0.0 && ei->electrodeGrossThickness > 0.0) {
        grossVolume = Pi * 0.25 * ei->electrodeGrossDiameter * ei->electrodeGrossDiameter * ei->electrodeGrossThickness;
    }

    /* Consider four cases
     * 1) specify the particle number and porosity --> calculate total superficial volume, renormalize solid mole numbers
     * 2) Specify the Total superficial volume of electrode and porosity --> calculate solid volume, renormalize solid mole numbers
     * 3) specify the mass/moles electrodes and Total superficial volume --> calculate porosity
     * 4) specify the mass/moles electrodes and porosity --> calculate total superficial volume
     */
    if (ei->particleNumberToFollow > 0.0) {
        particleNumberToFollow_ = ei->particleNumberToFollow;
        rewriteMoleNumbers = true;
        ElectrodeSolidVolume_ = particleNumberToFollow_ * partVol;
        if (ei->porosity > 0.0) {
            porosity_ = ei->porosity;
            grossVolume = ElectrodeSolidVolume_ / (1.0 - porosity_);
        } else {
            if (grossVolume > 0.0) {
                porosity_ = (grossVolume - ElectrodeSolidVolume_) / grossVolume;
            } else {
                throw Electrode_Error(" Electrode::electrode_model_create() ",
                                   " Need to either specify the total superficial volume or the porosity");
            }
            porosity_ = 0.3;
        }
        printf("Electrode::electrode_model_create \n\tWARNING: reset moles/masses "
               "based on \'Particle Number to Follow\' keyword.\n");
    } else if ((ei->porosity > 0.0) && (grossVolume > 0.0)) {
        porosity_ = ei->porosity;
        ElectrodeSolidVolume_ = grossVolume * (1.0 - porosity_);
        particleNumberToFollow_ = ElectrodeSolidVolume_ / partVol;
        rewriteMoleNumbers = true;
    } else if (Electrode::SolidVol() > 0.0) {
        if (grossVolume > 0.0) {
            ElectrodeSolidVolume_ = Electrode::SolidVol();
            porosity_ = (grossVolume - ElectrodeSolidVolume_) / grossVolume;
            if (porosity_ < 0.0) {
                throw Electrode_Error(" Electrode::electrode_model_create() ",
                                   " porosity_ is negative " + fp2str(porosity_));
            }
            particleNumberToFollow_ = ElectrodeSolidVolume_ / partVol;
        } else if (ei->porosity > 0.0) {
            ElectrodeSolidVolume_ = Electrode::SolidVol();
            particleNumberToFollow_ = ElectrodeSolidVolume_ / partVol;
            grossVolume = ElectrodeSolidVolume_ / (1.0 - ei->porosity);
        } else {
            throw Electrode_Error(
                    " Electrode::electrode_model_create() ",
                    " no specification of cell volume through \'Electrode Porosity\' "
                            "or \'Electrode Thickness\', etc., keywords. grossVol gives " + fp2str(grossVolume));
        }
    } else {
        throw Electrode_Error(
                " Electrode::electrode_model_create() ",
                " no specification of solid volume through \'Particle Number to Follow\', "
                        " Porosity and Superficial Volume, "
                        "\'Phase Moles\', or \'Phase Mass\' keyword.  SolidVol() gives "
                        + fp2str(Electrode::SolidVol()));
    }

    /*
     *  Calculate the total number of moles of species in the electrode that are consistent
     *  with the input parameters, when the input parameters have specified the total solid volume of the electrode.
     */
    if (rewriteMoleNumbers) {
        Electrode::resizeMoleNumbersToGeometry();
    }

    phaseMoles_init_ = phaseMoles_final_;
    phaseMoles_init_init_ = phaseMoles_final_;
    phaseMoles_final_final_ = phaseMoles_final_;

    /*
     *  Resize the value of the surface areas to ones appropriate for the input geometry.
     */
    Electrode::resizeSurfaceAreasToGeometry();

    Electrode::updateState_OnionOut();

    spMoles_init_init_ = spMoles_final_;
    spMoles_final_final_ = spMoles_final_;

    ElectrodeBath* BG = &(ei->m_BG);

    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        if (iph == solnPhase_ || iph == metalPhase_) {
            continue;
        }
        size_t kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        size_t nspPhase = tp.nSpecies();

        double* vv = BG->CapZeroDoDCoeffPhases[iph];
        double* vl = BG->CapLeftCoeffPhases[iph];
        for (size_t k = 0; k < nspPhase; k++) {
            size_t iGlobSpeciesIndex = kStart + k;
            capacityLeftSpeciesCoeff_[iGlobSpeciesIndex] = vl[k];
            capacityZeroDoDSpeciesCoeff_[iGlobSpeciesIndex] = vv[k];

        }

        //ei->SpeciesNames+kstart,
        //nSpecies
    }

    /*
     *  Process the Electrode Capacity Type information from the input file.
     */
    if (ei->electrodeCapacityType == 0) {
      electrodeCapacityType_ = CAPACITY_ANODE_ECT;
    } else if (ei->electrodeCapacityType == 1) {
      electrodeCapacityType_ = CAPACITY_CATHODE_ECT;
    }

    if (ei->methodCapacityCalc == 2 || ei->methodCapacityCalc == 3) {
        setCapacityCoeffFromInput(ei);
    } else {

    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {

        if (PhaseNames_[iph] == "LiFeS_X_combo") {
            electrodeChemistryModelType_ = 4;
            break;
        }
        if (PhaseNames_[iph] == "FeS2(S)") {
            electrodeChemistryModelType_ = 2;
        }
        if (PhaseNames_[iph] == "Li13Si4(S)") {
            if (electrodeChemistryModelType_ != 3) {
                electrodeChemistryModelType_ = 1;
            }
        }
        if (PhaseNames_[iph] == "Li7Si3_Interstitial") {
            electrodeChemistryModelType_ = 3;
        }
        if (PhaseNames_[iph] == "MCMB_Interstitials_anode") {
            electrodeChemistryModelType_ = 5;
        }
    }

    /*
     *  Now that we have identified the model type, make sure our list of phases is what
     *  we expect from the model type.
     */
    size_t nRSD = RSD_List_.size();


    if (electrodeChemistryModelType_ == 1) {
        pmatch(PhaseNames_, 2, "Li13Si4(S)");
        pmatch(PhaseNames_, 3, "Li7Si3(S)");
        if (nRSD > 1) {
            pmatch(PhaseNames_, 4, "Li12Si7(S)");
        }
        if (nRSD > 2) {
            pmatch(PhaseNames_, 5, "Si(S)");
        }
    }
    if (electrodeChemistryModelType_ == 2) {
        pmatch(PhaseNames_, 2, "FeS2(S)");
        pmatch(PhaseNames_, 3, "Li3Fe2S4(S)");
        if (nRSD > 1) {
            pmatch(PhaseNames_, 4, "Li[2+x]Fe[1-x]S2(S)");
            pmatch(PhaseNames_, 5, "Fe[1-x]S(S)");
        }
        if (nRSD > 2) {
            pmatch(PhaseNames_, 6, "Li2S(S)");
            pmatch(PhaseNames_, 7, "Fe(S)");
        }
    }
    if (electrodeChemistryModelType_ == 4) {
        pmatch(PhaseNames_, 2, "FeS2(S)");
        pmatch(PhaseNames_, 3, "Li3Fe2S4(S)");
        if (nRSD > 1) {
            pmatch(PhaseNames_, 4, "LiFeS_X_combo");
        }
        if (nRSD > 2) {
            pmatch(PhaseNames_, 5, "Li2S(S)");
            pmatch(PhaseNames_, 6, "Fe(S)");
        }
    }

    if (electrodeChemistryModelType_ == 1) {
        setCapacityCoeff_LiSi();
    } else if (electrodeChemistryModelType_ == 2) {
        setCapacityCoeff_FeS2();
    } else if (electrodeChemistryModelType_ == 3) {
        setCapacityCoeff_LiSi_Li();
    } else if (electrodeChemistryModelType_ == 4) {
        setCapacityCoeff_FeS2_Combo();
    } else if (electrodeChemistryModelType_ == 5) {
        setCapacityCoeff_MCMB();
    } else {
	Electrode_Types_Enum  etype = electrodeType();
	if (etype != MP_RXNEXTENT_ET && etype != MP_RXNEXTENT_FES2_ET &&
            etype != CSTR_ZN_ANODE_ET &&
            etype != CSTR_MCMB_ANODE_ET &&
            etype != CSTR_LICO2_CATHODE_ET &&
            etype != MP_RXNEXTENT_LISI_ET ) {
	    setCapacityCoeffFromInput(ei);
	}
    }
    }

    if (ei->RelativeCapacityDischargedPerMole != -1) {

        // setRelativeCapacityDischarged(ei->RelativeCapacityDischarged);
    }

    // Set the initial zeroDoD capacity by calling the capacity function at the beginning.
    capacityInitialZeroDod_ = capacity();
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];

    return 0;
}
//==================================================================================================================================
int Electrode::setInitialConditions(ELECTRODE_KEY_INPUT* ei)
{
    size_t nspecies = 0;
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {

        ThermoPhase* tphase = &(thermo(iph));
        size_t nSpecies = tphase->nSpecies();
        size_t kstart = nspecies;
        nspecies += nSpecies;

        /*
         *  We copy the mole numbers and mole fraction information from the bath gas section of the input file
         *  here. Note, the total quantitites of electrode material and electrolyte material may be overwritten
         *  later.
         */
        double sum = 0.0;
        for (size_t k = 0; k < nSpecies; k++) {
            spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k] = ei->MoleNumber[m_PhaseSpeciesStartIndex[iph] + k];
            sum += ei->MoleNumber[m_PhaseSpeciesStartIndex[iph] + k];
        }

        if (sum > 0.0) {
            for (size_t k = 0; k < nSpecies; k++) {
                spMf_final_[kstart + k] = spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k] / sum;
            }
            tphase->setMoleFractions(&(spMf_final_[kstart]));
        } else {
            tphase->getMoleFractions(&(spMf_final_[kstart]));
        }
        phaseVoltages_[iph] = tphase->electricPotential();
    }

    deltaVoltage_ = phaseVoltages_[metalPhase_] - phaseVoltages_[solnPhase_];

    /* set the absolute tolerance to 1e-5 the total number of moles */
    double tMoles = 0.0;
    for (size_t k = 0; k < m_NumTotSpecies; k++) {
        tMoles += spMoles_final_[k];
    }
    molarAtol_ = tMoles * 1.0E-5;

    /*
     * Now that spMoles_final_ is sized appropriately, we can call updateState_Phase()
     */
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        Electrode::updateState_Phase(iph);
    }

    /*
     *  Set up the particleDiameter_ and particleNumberToFollow_ fields.
     *   If any geometry information is supplied in the input file, then the geometry is
     *   used to rescale the input mole numbers of species. If no geometry information
     *   is input, the the partcleNumbertoFollow_ variable is calculated to be consistent
     *   with the input mole numbers.
     */
    // particle diameter is a required parameter
    if (ei->particleDiameter <= 0.0) {
        throw Electrode_Error("Electrode::electrode_model_create()", "Particle diameter is <= 0");
    }
    inputParticleDiameter_ = ei->particleDiameter;
    Radius_exterior_final_ = inputParticleDiameter_ / 2.0;

    double partVol = Pi * inputParticleDiameter_ * inputParticleDiameter_ * inputParticleDiameter_ / 6.0;
    bool rewriteMoleNumbers = false;
    /*
     *  Gross volume is the total superficial volume of the electrode. It is the sum of the
     *  solid volume and the liquid electrylte volume.
     *    - Compute the gross volume if this information has been supplied through the input file
     */
    double grossVolume = 0.0;
    if (ei->electrodeGrossArea > 0.0 && ei->electrodeGrossThickness > 0.0) {
        grossVolume = ei->electrodeGrossArea * ei->electrodeGrossThickness;
    } else if (ei->electrodeGrossDiameter > 0.0 && ei->electrodeGrossThickness > 0.0) {
        grossVolume = Pi * 0.25 * ei->electrodeGrossDiameter * ei->electrodeGrossDiameter * ei->electrodeGrossThickness;
    }

    /* Consider four cases
     * 1) Specify the particle number and porosity --> calculate total superficial volume, renormalize solid mole numbers
     * 2) Specify the Total superficial volume of electrode and porosity --> calculate solid volume, renormalize solid mole numbers
     * 3) Specify the mass/moles electrodes and Total superficial volume --> calculate porosity
     * 4) Specify the mass/moles electrodes and porosity --> calculate total superficial volume
     */
    if (ei->particleNumberToFollow > 0.0) {
        particleNumberToFollow_ = ei->particleNumberToFollow;
        rewriteMoleNumbers = true;
        ElectrodeSolidVolume_ = particleNumberToFollow_ * partVol;
        if (ei->porosity > 0.0) {
            porosity_ = ei->porosity;
            grossVolume = ElectrodeSolidVolume_ / (1.0 - porosity_);
        } else {
            if (grossVolume > 0.0) {
                porosity_ = (grossVolume - ElectrodeSolidVolume_) / grossVolume;
            } else {
                throw Electrode_Error(" Electrode::electrode_model_create() ",
                                   " Need to either specify the total superficial volume or the porosity");
            }
            porosity_ = 0.3;
        }
        printf("Electrode::electrode_model_create \n\tWARNING: reset moles/masses "
               "based on \'Particle Number to Follow\' keyword.\n");
    } else if ((ei->porosity > 0.0) && (grossVolume > 0.0)) {
        porosity_ = ei->porosity;
        ElectrodeSolidVolume_ = grossVolume * (1.0 - porosity_);
        particleNumberToFollow_ = ElectrodeSolidVolume_ / partVol;
        rewriteMoleNumbers = true;
    } else if (SolidVol() > 0.0) {
        if (grossVolume > 0.0) {
            ElectrodeSolidVolume_ = SolidVol();
            porosity_ = (grossVolume - ElectrodeSolidVolume_) / grossVolume;
            if (porosity_ < 0.0) {
                throw Electrode_Error(" Electrode::electrode_model_create() ",
                                   " porosity_ is negative " + fp2str(porosity_));
            }
            particleNumberToFollow_ = ElectrodeSolidVolume_ / partVol;
        } else if (ei->porosity > 0.0) {
            ElectrodeSolidVolume_ = SolidVol();
            particleNumberToFollow_ = ElectrodeSolidVolume_ / partVol;
            grossVolume = ElectrodeSolidVolume_ / (1.0 - ei->porosity);
        } else {
            throw Electrode_Error(
                    " Electrode::setInitialConditions() ",
                    " no specification of cell volume through \'Electrode Porosity\' "
                            "or \'Electrode Thickness\', etc., keywords. grossVol gives " + fp2str(grossVolume));
        }
    } else {
        throw Electrode_Error(
                " Electrode::setInitialConditions() ",
                " no specification of solid volume through \'Particle Number to Follow\', "
                        " Porosity and Superficial Volume, "
                        "\'Phase Moles\', or \'Phase Mass\' keyword.  SolidVol() gives " + fp2str(SolidVol()));
    }

    /*
     *  Calculate the total number of moles of species in the electrode that are consistent
     *  with the input parameters, when the input parameters have specified the total solid volume of the electrode.
     */
    if (rewriteMoleNumbers) {
        resizeMoleNumbersToGeometry();
    }

    /*
     *  Resize the value of the surface areas to ones appropriate for the input geometry.
     */
    resizeSurfaceAreasToGeometry();

    // Set the initial ZeroDoD capacity by calling the capacity() here.
    capacityInitialZeroDod_ = capacity();
    resetCapacityDischargedToDate();

    Electrode::updateState();
    Electrode::updateSurfaceAreas();

    /*
     *  New addition -> we should be able to put this here, but let's run test suite.
     */
    if (ei->RelativeCapacityDischargedPerMole != -1) {
        setRelativeCapacityDischargedPerMole(ei->RelativeCapacityDischargedPerMole);
	Electrode::updateState();
	Electrode::updateSurfaceAreas();
    }

    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];
    /*
     *  Now transfer the final state to the init and init init states.
     */
    Electrode::setInitStateFromFinal_Oin(true);

    return 0;
}
//==================================================================================================================================
int Electrode::electrode_stateSave_create(bool force)
{
   if (eState_save_) {
        if (!force) {
             return 0;
        }
        delete eState_save_;
    }
    eState_save_ = new EState();
    int rr = eState_save_->initialize(this);
    if (rr >= 0) {
        rr = 0;
    }
    if (!xmlStateData_final_) {
        eState_save_->copyElectrode_intoState(this);
        xmlStateData_final_ = eState_save_->write_electrodeState_ToXML();
    }
    return rr;
}
//==================================================================================================================================
// Set the sizes of the electrode from the input parameters
// This may be called after the model is created.
// Onion out routine (so we don't allow non-onion out virtual functions to call child objects)
/*
 *  We resize all of the information within the electrode from the input parameters
 *    electrodeArea   Area of the electrode
 *    electrodeThickness  Width of the electrode
 *    porosity        Volume of the electrolyte phase and other non-electrode volume phases
 */
void Electrode::setElectrodeSizeParams(double electrodeArea, double electrodeThickness, double porosity)
{
    if (pendingIntegratedStep_) {
        throw Electrode_Error("Electrode::setElectrodeSizeParams()", "called during a pending integration step");
    }
    porosity_ = porosity;
    double grossVolume = electrodeArea * electrodeThickness;
    double partVol = 4.0 * Pi * Radius_exterior_final_ * Radius_exterior_final_ * Radius_exterior_final_ / 3.0;
    if (grossVolume > 0.0) {
        ElectrodeSolidVolume_ = (1.0 - porosity_) * grossVolume;
        particleNumberToFollow_ = ElectrodeSolidVolume_ / partVol;
        // virtual function
        resizeMoleNumbersToGeometry();
        // virtual function
        resizeSurfaceAreasToGeometry();
    } else {
        throw Electrode_Error("Electrode::setElectrodeSizeParams()", "current volume is zero or less");
    }
    // reset the initial ZeroDoD capacity by calling the capacity() here.
    capacityInitialZeroDod_ = Electrode::capacity();
    Electrode::resetCapacityDischargedToDate();

    updateState();
    updateSurfaceAreas();

    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];
    /*
     *  Now transfer the final state to the init and init init states.
     */
    setInitStateFromFinal_Oin(true);
}
//==================================================================================================================================
// Resize the solid phase and electrolyte mole numbers within the object
/*
 *
 *  (Initialization safe, because this function doesn't call child routines)
 *
 *  This routine uses particleDiameter_ , particleNumberToFollow_, and porosity_ to recalculate
 *  all the mole numbers in the electrode. This is done by rescaling all of the numbers.
 *  At the end of the process, the total volume of the electrode object is
 *
 *    grossVol = SolidVol() / ( 1.0 - porosity_)
 *
 *  where the SolidVol() is equal to
 *
 *   SolidVol() =  particleNumberToFollow_  Pi *  particleDiameter_**3 / 6.0;
 *
 */
void Electrode::resizeMoleNumbersToGeometry()
{
    double partVol = 4.0 * Pi * Radius_exterior_final_ * Radius_exterior_final_ * Radius_exterior_final_ / 3.0;
    double targetSolidVol = partVol * particleNumberToFollow_;
    double currentSolidVol = Electrode::SolidVol();
    if (currentSolidVol <= 0.0) {
        throw Electrode_Error("Electrode::resizeMoleNumbersToGeometry()", "current solid volume is zero or less");
    }
    double ratio = targetSolidVol / currentSolidVol;
    double tMoles = 0.0;
    for (size_t k = 0; k < m_NumTotSpecies; k++) {
        spMoles_final_[k] *= ratio;
        spMoles_init_[k] = spMoles_final_[k];
        spMoles_init_init_[k] = spMoles_final_[k];
        tMoles += spMoles_final_[k];
    }
    molarAtol_ = tMoles * 1.0E-5;
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        updateState_Phase(iph);
    }
    phaseMoles_init_ = phaseMoles_final_;
    phaseMoles_init_init_ = phaseMoles_final_;
    phaseMoles_final_final_ = phaseMoles_final_;

    /*
     * Resize the electrolyte so that the total volume of the electrolyte is consistent with the given
     * porosity, porosity_
     */
    // Electrode::resizeSolutionNumbersToPorosity();
    currentSolidVol = Electrode::SolidVol();
    double grossVol = currentSolidVol / (1.0 - porosity_);
    double targetSolnVol = grossVol - currentSolidVol;
    double totalVol = Electrode::TotalVol();
    double currentSolnVol = totalVol - currentSolidVol;
    ratio = targetSolnVol / currentSolnVol;

    size_t istart = m_PhaseSpeciesStartIndex[solnPhase_];
    size_t nspSoln = m_PhaseSpeciesStartIndex[solnPhase_ + 1] - m_PhaseSpeciesStartIndex[solnPhase_];
    for (size_t kk = 0; kk < nspSoln; kk++) {
        size_t k = istart + kk;
        spMoles_final_[k] *= ratio;
        spMoles_init_[k] = spMoles_final_[k];
        spMoles_init_init_[k] = spMoles_final_[k];
        spMoles_final_final_[k] = spMoles_final_[k];
    }
    updateState_Phase(solnPhase_);
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];

    currentSolidVol = Electrode::SolidVol();
    totalVol = Electrode::TotalVol();
    double calcPor = (totalVol - currentSolidVol) / totalVol;
    if (fabs(calcPor - porosity_) > 1.0E-6) {
        throw Electrode_Error("Electrode::resizeMoleNumbersToGeometry() Error",
                              "Couldn't set the porosity correctly: " + fp2str(calcPor) + " vs " + fp2str(porosity_));
    }
}
//====================================================================================================================
//   Resize the electrolyte mole numbers so that it consistent with the porosity
/*
 *  (serial virtual function -> calls children)
 *
 *  This routine resizes the electrolyte mole numbers within the object to be consistent
 *  with the porosity. This has the effect of changing the total volume of the
 *  object.
 *
 *  @param porosityReset   If positive, this resets the porosity within the electrode
 *                         to the given value.
 */
/*
void Electrode::resizeSolutionNumbersToPorosity(double porosityReset)
{
    if (porosityReset > 0.0) {
        porosity_ = porosityReset;
    }
    double currentSolidVol = SolidVol();
    double grossVol = currentSolidVol / (1.0 - porosity_);
    double targetSolnVol = grossVol - currentSolidVol;
    double currentSolnVol = TotalVol() - currentSolidVol;
    double ratio = targetSolnVol / currentSolnVol;

    int istart = m_PhaseSpeciesStartIndex[solnPhase_];
    int nspSoln = m_PhaseSpeciesStartIndex[solnPhase_ + 1] - m_PhaseSpeciesStartIndex[solnPhase_];
    for (int kk = 0; kk < nspSoln; kk++) {
        int k = istart + kk;
        spMoles_final_[k] *= ratio;
        spMoles_init_[k] = spMoles_final_[k];
        spMoles_init_init_[k] = spMoles_final_[k];
    }
    updateState_Phase(solnPhase_);

}
*/
//====================================================================================================================
//    Resize the surface areas according to the input geometry.
/*
 *  We resize the surface areas of the Reacting Surfaces to a value which is
 *  equal to the particle exterior surface area multiplied by the number of particles.
 *
 *  Note we can't do anything better in this routine because we do not know the specifics
 *  of the particles.
 */
void Electrode::resizeSurfaceAreasToGeometry()
{
    double totalSurfArea = 4.0 * Pi * Radius_exterior_final_ * Radius_exterior_final_ * particleNumberToFollow_;

    double total = 0.0;
    double maxVal = -1.0;
    for (size_t i = 0; i < numSurfaces_; i++) {
        if (maxVal < surfaceAreaRS_final_[i]) {
            maxVal = surfaceAreaRS_final_[i];
        }
        total += surfaceAreaRS_final_[i];
    }

    if (total == 0.0) {
        for (size_t i = 0; i < numSurfaces_; i++) {
            surfaceAreaRS_final_[i] = totalSurfArea;
            surfaceAreaRS_init_[i] = totalSurfArea;
            surfaceAreaRS_init_init_[i] = totalSurfArea;
            surfaceAreaRS_final_final_[i] = totalSurfArea;

            if (RSD_List_.size() > 0) {
                if (RSD_List_[i]) {
                    ActiveKineticsSurf_[i] = 1;
                } else {
                    ActiveKineticsSurf_[i] = 0;
                }
            }
        }
    } else {
        double ratio = totalSurfArea / maxVal;
        for (size_t i = 0; i < numSurfaces_; i++) {
            surfaceAreaRS_final_[i] *= ratio;
            surfaceAreaRS_init_[i] = surfaceAreaRS_final_[i];
            surfaceAreaRS_init_init_[i] = surfaceAreaRS_final_[i];
            surfaceAreaRS_final_final_[i] = surfaceAreaRS_final_[i];
            ActiveKineticsSurf_[i] = 0;
            if (surfaceAreaRS_final_[i] > 0.0) {
                if (RSD_List_.size() > 0) {
                    if (RSD_List_[i]) {
                        ActiveKineticsSurf_[i] = 1;
                    }
                }
            }
        }
    }

}
//====================================================================================================================
// We change the particleNumbertoFollow_ field to comply with the number of moles and the particle diameter
void Electrode::resizeParticleNumbersToMoleNumbers()
{
    double partVol = 4.0 * Pi * Radius_exterior_final_ * Radius_exterior_final_ * Radius_exterior_final_ / 3.0;
    double currentSolidVol = SolidVol();
    particleNumberToFollow_ = currentSolidVol / partVol;
    resizeSurfaceAreasToGeometry();
}
//====================================================================================================================
ZZCantera::ReactingSurDomain* Electrode::currOuterReactingSurface()
{
    for (size_t isk = 0; isk < numSurfaces_; isk++) {
        if (ActiveKineticsSurf_[isk]) {
            return RSD_List_[isk];
        }
    }
    return 0;
}
//==================================================================================================================================
ZZCantera::ReactingSurDomain* Electrode::reactingSurface(size_t iSurf)
{
    return RSD_List_[iSurf];
}
//==================================================================================================================================
// Set the mole numbers in the electrolyte phase
/*
 *  We set the mole numbers of the electrolyte phase separately from
 *  the rest of the phases.
 *
 *  We always make sure that mole numbers are positive by clipping. We
 *  always make sure that mole fractions sum to one.
 *
 *  If we are not following mole numbers in the electrode, we set the
 *  total moles to the internal constant, electrolytePseudoMoles_, while
 *  using this vector to set the mole fractions.
 *
 * @param electrolyteMoleNum vector of mole numbers of the species in the
 *                     electrolyte phase.
 *                     units = kmol
 *                     size = number of species in the electrolyte phase
 *
 * @param setInitial   Boolean indicating that we should set the initial values of the
 *                     electrolyte mole numbers as well. 
 */
void Electrode::setElectrolyteMoleNumbers(const double* const electrolyteMoleNum, bool setInitial)
{
    size_t istart = m_PhaseSpeciesStartIndex[solnPhase_];
    ThermoPhase& tp = thermo(solnPhase_);
    size_t nsp = tp.nSpecies();
    AssertTrace(nsp == (m_PhaseSpeciesStartIndex[solnPhase_+1] - m_PhaseSpeciesStartIndex[solnPhase_]));
    double tmp = 0.0;
    for (size_t k = 0; k < nsp; k++) {
        spMoles_final_[istart + k] = std::max(electrolyteMoleNum[k], 0.0);
        tmp += spMoles_final_[istart + k];
    }
    phaseMoles_final_[solnPhase_] = tmp;

    if (tmp > 1.0E-200) {
        for (size_t k = 0; k < nsp; k++) {
            spMf_final_[istart + k] = spMoles_final_[istart + k] / tmp;
        }
        /*
         *  @TODO: take this section out. If we are setting the moles at the interface,
         *         then we should set the moles. See the explanation in the Electrode.h file.
         */
        if (!followElectrolyteMoles_) {
            for (size_t k = 0; k < nsp; k++) {
                spMoles_final_[istart + k] *= electrolytePseudoMoles_ / tmp;
            }
            phaseMoles_final_[solnPhase_] = electrolytePseudoMoles_;
        }
    }
    /*
     *  Update the internal state 
     */
    updateState();

    /*
     *  Set other states retroactively
     */
    if (setInitial) {
	if (pendingIntegratedStep_) {
	    throw Electrode_Error("Electrode::setElectrolyteMoleNumbers ERROR",
			       "Trying to set initial electroltye mole numbers");
	}
        for (size_t k = 0; k < nsp; k++) {
            spMoles_init_[istart + k] = spMoles_final_[istart + k];
            spMoles_init_init_[istart + k] = spMoles_final_[istart + k];
        }
        phaseMoles_init_[solnPhase_] = phaseMoles_final_[solnPhase_];
        phaseMoles_init_init_[solnPhase_] = phaseMoles_final_[solnPhase_];
        phaseMoles_final_final_[solnPhase_] = phaseMoles_final_[solnPhase_];
    

	/*
	 *  Set the init and the init_init state and the final_final state from the final state
	 */
	setInitStateFromFinal();
	setFinalFinalStateFromFinal();
    }
}
//================================================================================================
// Set the mole numbers in a single phase
/*
 *  We set the mole numbers of a single phase separately from
 *  the rest of the phases.
 *
 *  We always make sure that mole numbers are positive by clipping. We
 *  always make sure that mole fractions sum to one.
 *
 *  If we are not following mole numbers in the electrode, we set the
 *  total moles to the internal constant, electrolytePseudoMoles_, while
 *  using this vector to set the mole fractions.
 *
 * @param iph     Phase id.
 * @param moleNum vector of mole numbers of the species in the
 *                     electrolyte phase.
 *                     units = kmol
 *                     size = number of species in the electrolyte phase
 */
void Electrode::setPhaseMoleNumbers(size_t iph, const double* const moleNum)
{
    if (pendingIntegratedStep_) {
        if (iph != solnPhase_) {
            throw Electrode_Error("Electrode::setPhaseMoleNumbers", "called when there is a pending step");
        }
    }
    size_t istart = m_PhaseSpeciesStartIndex[iph];
    size_t nsp = m_PhaseSpeciesStartIndex[iph + 1] - istart;
    for (size_t k = 0; k < nsp; k++) {
        spMoles_final_[istart + k] = std::max(moleNum[k], 0.0);
        spMoles_init_[istart + k] = spMoles_final_[istart + k];
        spMoles_init_init_[istart + k] = spMoles_final_[istart + k];
    }
    updateState_Phase(iph);

    spMoles_init_init_ = spMoles_final_;
    spMoles_init_ = spMoles_final_;
    spMoles_final_final_ = spMoles_final_;
    phaseMoles_init_init_ = phaseMoles_final_;
    phaseMoles_init_ = phaseMoles_final_;
    phaseMoles_final_final_ = phaseMoles_final_;

    /*
     *  Change the particleNumbersToFollow_ field
     */
    resizeParticleNumbersToMoleNumbers();

    /*
     *  Call the virtual function updateState. This may fill in further state information that
     *  is not known to this function. 
     */
    updateState();

    /*
     *  Set the init and the init_init state and the final_final state from the final state
     */
    setInitStateFromFinal(true);
    setFinalFinalStateFromFinal();
}
//==================================================================================================================
//! Set the current state of the electrode object based on a relative extent of reaction
/*!
 *   The relative extent of reaction is a dimensionless number on the order of one
 *   that represents the state of the electrode. A value of zero represents the
 *   fully charged state, while a value of one (or equivalent) represents a fully
 *   discharged state. 
 *
 *  @param  relativeExtentRxn  Relative extent of reaction variable (input)
 */
void Electrode::setState_relativeExtentRxn(double relativeExtentRxn)
{
    throw Electrode_Error("Electrode::setState_relativeExtentRxn",
	                  "Base class called but not implemented");
}

//==================================================================================================================================
double Electrode::relativeExtentRxn(double time) const
{
    double cap, capLeft, rel = -1;;
    bool ok = false;
    if (fabs(time - tfinal_) < 1.0E-50) {
        ok = true;
        cap = capacityInitial();
        capLeft = capacityLeft();
        rel = (cap - capLeft) / cap;
        return rel;
    }
    if (fabs(time - tinit_) < 1.0E-50) {
        ok = true;
    }
    if (fabs(time - t_init_init_) < 1.0E-50) {
        ok = true;
    }
    if (fabs(time - t_final_final_) < 1.0E-50) {
        ok = true;
    }
    if (!ok) {
        throw Electrode_Error("Electrode::relativeExtentRxn(double time)", "unknown time: " + fp2str(time));
    }
    throw Electrode_Error("Electrode::relativeExtentRxn(double time)", "Base class called but not implemented");
    return 0.0;
}
//==================================================================================================================================


// -----------------------------------------------------------------------------------------------------------------
// ------------------------------------ CARRY OUT INTEGRATION OF EQUATIONS -----------------------------------------
// -----------------------------------------------------------------------------------------------------------------

//==================================================================================================================================
void Electrode::setTime(double time)
{
    // if we are in a pending step, this is an error.
    if (pendingIntegratedStep_) {
        throw Electrode_Error("Electrode::setTime()", "called when there is a pending step");
    }
    t_init_init_ = time;
    tfinal_ = t_init_init_;
    tinit_ = t_init_init_;
    t_final_final_ = t_init_init_;
}
//==================================================================================================================================
// Sets the state of the Electrode object given an EState object
/*
 *   (virtual function)
 *   This sets all of the states within the object to the same state.
 *   It is an error to call this function during a pending step where there can be a difference between t_init and t_final.
 *
 *   @param[in]  es          const reference to the EState object.  Must be the correct EState object for the
 *                           current Electrode object, or else it will throw an error. However, there is an option to 
 *                           read EState objects with less information. 
 */
void Electrode::setState_EState(const EState& es)
{
    Electrode::setState_EStateBase(es);
}
//==================================================================================================================================
void Electrode::setState_EStateBase(const EState& es)
{
    if (pendingIntegratedStep_) {
        throw Electrode_Error("Electrode::setState_EState()", "called when there is a pending step");
    }
    spMoles_final_                     = es.spMoles_;
    phaseVoltages_                     = es.phaseVoltages_;
    temperature_                       = es.temperature_;
    pressure_                          = es.pressure_;
    electrodeChemistryModelType_       = es.electrodeChemistryModelType_;
    electrodeDomainNumber_             = es.electrodeDomainNumber_;
    electrodeCellNumber_               = es.electrodeCellNumber_;
    particleNumberToFollow_            = es.particleNumberToFollow_;
    ElectrodeSolidVolume_              = es.electrodeSolidVolume_;
    // es.grossVolume_ 
    Radius_exterior_final_             = es.radiusExterior_;
    surfaceAreaRS_final_               = es.surfaceAreaRS_;
    // electrodeMoles_                    = es.electrodeMoles_;
    electrodeCapacityType_             = es.electrodeCapacityType_;
    //capacityLeft_                      = es.capacityLeft_; (no explicit storage of this object in Electrode object)
    capacityInitialZeroDod_            = es.capacityInitial_;

    // depthOfDischarge_                  = es.depthOfDischarge_; (no explicit storage of this object in Electrode object)
    depthOfDischargeStarting_          = es.depthOfDischargeStarting_;
    // relativeElectronsDischargedPerMole_
    // relativeDepthOfDischarge_          = es.relativeDepthOfDischarge_ 
    // capacityDischargedToDate_
    electronKmolDischargedToDate_      = es.electronKmolDischargedToDate_;
  
    deltaTsubcycle_init_next_          = es.deltaTsubcycle_init_next_;

    // Calculate the porosity from the EState gross volume and solid volume
    double grossVol = es.grossVolume_;
    double currentSolidVol = es.electrodeSolidVolume_;
    porosity_ = currentSolidVol / grossVol;

    // Update the _final_ state locally. Calculate phase moles and phase mole fractions
    Electrode::updateState();

    // Set the _init_ and _init_init_ states from the final state
    Electrode::setInitStateFromFinal(true);
    // Set the _final_final_ state 
    Electrode::setFinalFinalStateFromFinal();
}
//==================================================================================================================================
// Update all mole numbers in the object from the mole numbers in the spMoles_final_[] vector
/*
 *  We use the field spMoles_final_[] to set the field phaseMoles_final_[].
 *  We set the mole numbers of a single phase separately from the rest of the phases.
 *  We do not clip the mole numbers to be positive. We allow negative mole numbers.
 *  We make sure that the mole fractions sum to one.
 *  The following fields are state variables
 *
 *  Summary: State Variables
 *              spMoles_final_[kGlob]  :  iph entries in Global Species PhaseList vector
 *
 *              phaseMoles_final_[iph]  : When the phase moles are negative
 *              spMf_final_[kGlobal]      iph entries in Global Species PhaseList vector
 *
 *           Independent Variables
 *              spMoles_final_[kGlobal]  :    kGlobal in electrolyte phase
 *              phaseVoltages_[iph] :         iph in solid phase electrodes
 *              Temperature
 *              Pressure
 *
 *   The following member data fields are set by this function:
 *
 *            spElectroChemPot_[]            iph entries in Global Species PhaseList vector
 *            phaseMoles_final_[iph]         iph entries in Global Phase PhaseList vector
 *            spMf_final_[]                  iph entries in Global Species PhaseList vector
 *            VolPM_[]                       iph entries in Global Species PhaseList vector
 *            spElectroChemPot_[]            iph entries in Global Species PhaseList vector
 *            ThermoPhase[iph]               setState using TPX and electricPotential
 *            phaseMolarVolumes_[]           iph entries in Global Species PhaseList vector 
 *
 *  If we are not following  the mole numbers in the electrode, we set the
 *  total moles to the internal constant, electrolytePseudoMoles_, while
 *  using this vector to set the mole fractions, using the ThermoPhase object
 *  to get the mole fractions.
 *
 * A Note about negative Phase Numbers.
 * -----------------------------------------------
 *     The algorithms shouldn't rely on this routine to handle negative mole numbers or negative phase numbers.
 *  The algorithms should already have arranged this issue before calling this routine.
 *  However, this routine does have a method for handling this. For positive mole numbers, the spMoles vector
 *  is treated as the state vector. For negative mole numbers,
 */
void Electrode::updateState_Phase(size_t iph)
{
    AssertTrace(iph != npos);
    size_t istart = m_PhaseSpeciesStartIndex[iph];
    ThermoPhase& tp = thermo(iph);
    size_t nsp = m_PhaseSpeciesStartIndex[iph + 1] - istart;
    double tmp = 0.0;
    for (size_t k = 0; k < nsp; k++) {
        tmp += spMoles_final_[istart + k];
    }
    phaseMoles_final_[iph] = tmp;
    if (tmp > 1.0E-200) {
        for (size_t k = 0; k < nsp; k++) {
            spMf_final_[istart + k] = spMoles_final_[istart + k] / tmp;
        }
#ifdef OLD_FOLLOW
        if (iph == solnPhase_) {
            if (!followElectrolyteMoles_) {
                ThermoPhase& tp = thermo(solnPhase_);
                tp.getMoleFractions(&spMf_final_[istart]);
                for (size_t k = 0; k < nsp; k++) {
                    spMoles_final_[istart + k] = electrolytePseudoMoles_ * spMf_final_[istart + k];
                }
                phaseMoles_final_[iph] = electrolytePseudoMoles_;
            }
        }
#endif
    } else if (tmp < -1.0E-200) {
        /*
         *  We are here when the phase moles are less than zero. This is a permissible condition when we
         *  are calculating a phase death, as we solve for the time when the phase moles is exactly zero.
         *  Letting the phase moles go negative during the solve is a necessary requirement
         */
        /*
         * For the current method the spMf_final_[] contains the valid information.
         */
        for (size_t k = 0; k < nsp; k++) {
            tmp = spMf_final_[istart + k];
            if (tmp < 0.0 || tmp > 1.0) {
                throw Electrode_Error("Electrode::updatePhaseNumbers()",
                                   "Mole fractions out of bounds:" + int2str(k) + " " + fp2str(tmp));
            }
            // TODO: this may be wrong. Experiment with leaving this alone
            spMoles_final_[istart + k] = 0.0;
        }
    } else {
        // We are here when the mole numbers of the phase are zero. In this case, we respect the
	// mole fraction vector in spMf_final_[]. There are cases where we are solving for the
	// mole fraction vector when the phase is zero, and we are just trying to get a 1 on the
	// diagonal. In this case we need to respect  spMF_final_[].
	// However we check for out of bounds values.
        for (size_t k = 0; k < nsp; k++) {
            tmp = spMf_final_[istart + k];
            if (tmp < 0.0 || tmp > 1.0) {
                throw Electrode_Error("Electrode::updatePhaseNumbers()",
                                   "Mole fractions out of bounds:" + int2str(k) + " " + fp2str(tmp));
            }
            spMoles_final_[istart + k] = 0.0;
            phaseMoles_final_[iph] = 0.0;
        }
    }

    // Here we set the state within the phase object
    tp.setState_TPX(temperature_, pressure_, &spMf_final_[istart]);
    tp.setElectricPotential(phaseVoltages_[iph]);
    tp.getPartialMolarVolumes(&(VolPM_[istart]));
    tp.getElectrochemPotentials(&(spElectroChemPot_[istart]));
    if (iph < m_NumVolPhases) {
        phaseMolarVolumes_[iph] = tp.molarVolume();
    } else {
        // compatibilty
        //phaseMolarVolumes_[iph] = tp.molarVolume();
        phaseMolarVolumes_[iph] = 0.0;
        size_t isurf = iph - m_NumVolPhases;
        sphaseMolarAreas_[isurf] = tp.molarArea();
    }
    // Right now we use the Cp calculation from Cantera until we expand Cantera to calculate Cv
    if (doThermalPropertyCalculations_) {
       tp.getPartialMolarCp(&(CvPM_[istart]));
       tp.getPartialMolarEnthalpies(&(enthalpyMolar_final_[istart]));
       tp.getPartialMolarEntropies(&(entropyMolar_final_[istart]));
       tp.getChemPotentials(&(chempotMolar_final_[istart]));
    }
}
//==================================================================================================================================
/*
void Electrode::updatePhaseNumbersTmp(const std::vector<double>& spMoles_tmp, std::vector<double>& phaseMoles_tmp,
        std::vector<double>& spMf_tmp) const
{
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        size_t istart = m_PhaseSpeciesStartIndex[iph];
        // ThermoPhase &tp = thermo(solnPhase_);
        size_t nsp = m_PhaseSpeciesStartIndex[iph + 1] - istart;
        double tmp = 0.0;
        for (size_t k = 0; k < nsp; k++) {
            tmp += spMoles_tmp[istart + k];
        }
        phaseMoles_tmp[iph] = tmp;
        if (tmp > 1.0E-200) {
            for (size_t k = 0; k < nsp; k++) {
                spMf_tmp[istart + k] = spMoles_tmp[istart + k] / tmp;
            }
        } else {
            // We are here when the mole numbers of the phase are zero. In this case, we still need
            // a valid mole fraction vector. This is kept in the ThermoPhase object.
            ThermoPhase& tp = thermo(iph);
            tp.getMoleFractions(&(spMf_tmp[istart]));
        }
    }
}
*/
//==================================================================================================================================
/*
 *  Calculates the change in the surface area of all external and internal interfaces within the electrode
 *  variables to be potentially altered
 *   surfaceAreaRS_[];
 *   isExternalSurface[]
 *   numExternalInterfacialSurfaces_;
 */
double Electrode::calcSurfaceAreaChange(double deltaT)
{
    double sa_final = surfaceAreaRS_init_[0];
    return sa_final;
}
//==================================================================================================================================
double Electrode::updateElectrolytePseudoMoles()
{
    /*
     * Call updateState to update spMf_final_ and phaseMoles_final_
     *   We also save the value of   electrolytePseudoMoles_  here.
     */
    if (pendingIntegratedStep_) {
        throw Electrode_Error(" Electrode::updateElectrolytePseudoMoles()",
                           " call is illegal when there is a pending step");
    }
    int origFEM = followElectrolyteMoles_;
    followElectrolyteMoles_ = 1;
    updateState();
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];
    double currentSolidVol = SolidVol();
    double totalVol = TotalVol();
    double currentSolnVol = totalVol - currentSolidVol;
    phaseMoles_final_[solnPhase_] = currentSolnVol / phaseMolarVolumes_[solnPhase_];
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];
    followElectrolyteMoles_ = 0;
    updateState();
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];

    followElectrolyteMoles_ = 0;
    followElectrolyteMoles_ = origFEM;
    return electrolytePseudoMoles_;
}
//==================================================================================================================================
void Electrode::turnOffFollowElectrolyteMoles()
{
    followElectrolyteMoles_ = 0;
    if (phaseMoles_final_[solnPhase_] <= 0.0) {
        throw Electrode_Error("Electrode::turnOffFollowElectrolyteMoles()", "electrolyte moles set negative or zero");
    }
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];
    updateState_Phase(solnPhase_);
}
//==================================================================================================================================
void Electrode::turnOnFollowElectrolyteMoles()
{
    if (followElectrolyteMoles_ == 0) {
        phaseMoles_final_[solnPhase_] = electrolytePseudoMoles_;
    }

    followElectrolyteMoles_ = 1;
    updateState_Phase(solnPhase_);

    if (phaseMoles_final_[solnPhase_] <= 0.0) {
        if (electrolytePseudoMoles_ > 0.0) {
            followElectrolyteMoles_ = 0;
            updateState_Phase(solnPhase_);
            followElectrolyteMoles_ = 1;
        }
        if (phaseMoles_final_[solnPhase_] <= 0.0) {
            throw Electrode_Error("Electrode::turnOnFollowElectrolyteMoles()", "electrolyte moles set negative or zero");
        }
    }
}
//==================================================================================================================================
const std::vector<bool>& Electrode::getExternalSurfaceBooleans() const
{
    return isExternalSurface_;
}
//==================================================================================================================================
double Electrode::phaseMoles(size_t iph) const
{
    return phaseMoles_final_[iph];
}
//==================================================================================================================================
double Electrode::elementMoles(size_t ieG) const
{
    const Elements* elemList = globalElements();
    std::string eN = elemList->elementName(ieG);
    return elementMoles(eN);
}
//==================================================================================================================================
double Electrode::elementMoles(std::string eN) const
{
    double sum = 0.0;
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        const ThermoPhase* tp_ptr = &thermo(iph);
        size_t ie = tp_ptr->elementIndex(eN);
        if (ie != npos) {
            size_t nsp = thermo(iph).nSpecies();
            size_t kStart = m_PhaseSpeciesStartIndex[iph];
            for (size_t ik = 0; ik < nsp; ik++) {
                double na = thermo(iph).nAtoms(ik, ie);
                sum += na * spMoles_final_[kStart + ik];
            }
        }
    }
    return sum;
}
//==================================================================================================================================
double Electrode::elementSolidMoles(size_t ieG) const
{
    const Elements* elemList = globalElements();
    std::string eN = elemList->elementName(ieG);
    return elementMoles(eN);
}
//==================================================================================================================================
double Electrode::elementSolidMoles(std::string eN) const
{
    double sum = 0.0;
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        if (iph != solnPhase_) {
            const ThermoPhase* tp_ptr = &thermo(iph);
            size_t ie = tp_ptr->elementIndex(eN);
            if (ie != npos) {
                size_t nsp = thermo(iph).nSpecies();
                size_t kStart = m_PhaseSpeciesStartIndex[iph];
                for (size_t ik = 0; ik < nsp; ik++) {
                    double na = thermo(iph).nAtoms(ik, ie);
                    sum += na * spMoles_final_[kStart + ik];
                }
            }
        }
    }
    return sum;
}
//================================================================================================
double Electrode::speciesElectrochemPotential(size_t iGlobalSpIndex) const
{
    // DANGER -> Check to see this is updated

    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        ThermoPhase* tphase = &(thermo(iph));
        size_t kStart = m_PhaseSpeciesStartIndex[iph];
        tphase->getElectrochemPotentials(&(spElectroChemPot_[kStart]));
    }
    return spElectroChemPot_[iGlobalSpIndex];
}
//================================================================================================
double Electrode::speciesChemPotential(size_t iGlobalSpIndex) const
{
    // DANGER -> Check to see this is updated

    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        ThermoPhase* tphase = &(thermo(iph));
        size_t kStart = m_PhaseSpeciesStartIndex[iph];
        tphase->getChemPotentials(&(spElectroChemPot_[kStart]));
    }
    double tmp = spElectroChemPot_[iGlobalSpIndex];
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        ThermoPhase* tphase = &(thermo(iph));
        size_t kStart = m_PhaseSpeciesStartIndex[iph];
        tphase->getElectrochemPotentials(&(spElectroChemPot_[kStart]));
    }
    return tmp;
}
//================================================================================================
void Electrode::getMoleFractions(double* const x) const
{
    std::copy( spMf_final_.begin(), spMf_final_.end(), x );
}
//================================================================================================
void Electrode::getMoleNumSpecies(double* const n) const
{
    std::copy( spMoles_final_.begin(), spMoles_final_.end(), n );
}
//================================================================================================
const std::vector<double> & Electrode::getMoleNumSpecies() const
{
  return spMoles_final_;
}
//================================================================================================
void Electrode::getMoleNumPhases(double* const np) const
{
    std::copy( phaseMoles_final_.begin(), phaseMoles_final_.end(), np);
}
//================================================================================================
double Electrode::moleFraction(size_t globalSpeciesIndex) const
{
    return spMf_final_[globalSpeciesIndex];
}
//==================================================================================================================================
double Electrode::moleNumSpecies(size_t globalSpeciesIndex) const
{
    return spMoles_final_[globalSpeciesIndex];
}
//==================================================================================================================================
// Returns the index of a phase in the ReactionSurfaceDomain object
// given the index of that phase in the PhaseList object
/*
 * @param PLph index of the phase in the PhaseList object, which is also the
 *             Electrode_Model object.
 *
 *  @return  Returns the index of the phase in the current ReactingSurDomain
 *           object. A value of -1 in this slot means that the phase doesn't
 *           participate in the  current ReactingSurDomain object
 */
size_t Electrode::ReactingSurfacePhaseIndex(size_t isk, size_t PLph) const
{
    ReactingSurDomain* rsd = RSD_List_[isk];
    if (!rsd) {
        throw Electrode_Error("ReactingSurfacePhaseIndex", "Reacting surface not found");
    }
    return (int) rsd->PLtoKinPhaseIndex_[PLph];
}
//==================================================================================================================================
void Electrode::setState_TP(double temperature, double pressure)
{
    temperature_ = temperature;
    pressure_ = pressure;
    updateState();
}
//==================================================================================================================================
double Electrode::temperature() const
{
    return temperature_;
}
//==================================================================================================================================
double Electrode::pressure() const
{
    return pressure_;
}
//==================================================================================================================================
double Electrode::porosity() const
{
    return porosity_;
}
//==================================================================================================================================
double Electrode::SolidVol() const
{
    double vol = 0.0;
    for (size_t iph = 0; iph < m_NumVolPhases; iph++) {
        size_t kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        size_t nspPhase = tp.nSpecies();
        if (iph != solnPhase_) {
            for (size_t k = 0; k < nspPhase; k++) {
                vol += spMoles_final_[kStart + k] * VolPM_[kStart + k];
            }
        }
    }
    return vol;
}
//==================================================================================================================================
double Electrode::SolidTotalMoles() const
{
    double tot = 0.0;
    for (size_t iph = 0; iph < m_NumVolPhases; iph++) {
        size_t kStart = m_PhaseSpeciesStartIndex[iph];
        //ThermoPhase& tp = thermo(iph);
        //size_t nspPhase = tp.nSpecies();
        size_t nspPhase = NumSpeciesList_[iph];
        if (iph != solnPhase_ && (iph != metalPhase_)) {
            for (size_t k = 0; k < nspPhase; k++) {
                tot += spMoles_final_[kStart + k];
            }
        }
    }
    return tot;
}
//==================================================================================================================================
double Electrode::TotalVol(bool ignoreErrors) const
{
    double vol = 0.0;

    for (size_t iph = 0; iph < m_NumVolPhases; iph++) {
        double psum = 0.0;
        size_t kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        size_t nspPhase = tp.nSpecies();
        for (size_t k = 0; k < nspPhase; k++) {
            vol += spMoles_final_[kStart + k] * VolPM_[kStart + k];
            psum += spMoles_final_[kStart + k] * VolPM_[kStart + k];
        }
        double mv = tp.molarVolume();
        // compatibility statement until it is changed
        if (iph >= m_NumVolPhases) {
            mv = 0.0; 
        }
        double palt = mv * phaseMoles_final_[iph];
        if (!ignoreErrors) {
            if (palt < -1.0E-15) {
                throw Electrode_Error(" Electrode::TotalVol() ", " phase volume is negative " + fp2str(palt));
            }
            if (psum < -1.0E-15) {
                throw Electrode_Error(" Electrode::TotalVol() ", " phase volume is negative " + fp2str(psum));
            }
        }
        double denom = palt + psum + 1.0E-9;
        if (tp.eosType() != cLattice) {
            if (!ignoreErrors && 0) {
                if (fabs((palt - psum) / denom) > 1.0E-4) {
                    throw Electrode_Error(" Electrode::TotalVol() ",
                                       " internal inconsistency " + fp2str(palt) + " " + fp2str(psum));
                }
            }
        }
    }
    return vol;
}
//==================================================================================================================================
void Electrode::getPhaseVol(double* const phaseVols) const
{
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        phaseVols[iph] = 0.0;
        size_t kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        size_t nspPhase = tp.nSpecies();
        for (size_t k = 0; k < nspPhase; k++) {
            phaseVols[iph] += spMoles_final_[kStart + k] * VolPM_[kStart + k];
        }
        if (iph >= m_NumVolPhases) {
            phaseVols[iph] = 0.0;
        }
    }
}
//====================================================================================================================
double Electrode::SolidHeatCapacityCV() const
{
    //
    //   Calculate CvPM_ on the fly if it isn't part of the normal calculation.
    //         -> debatable about how to do this most efficiently.
    //
    if (!doThermalPropertyCalculations_) {
	for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
	    int istart = m_PhaseSpeciesStartIndex[iph];
	    ThermoPhase& tp = thermo(iph);
	    tp.getPartialMolarCp(&(CvPM_[istart]));
	}
    }

    double heatCapacity = 0.0;
    for (size_t iph = 0; iph < m_NumVolPhases; iph++) {
        size_t kStart = m_PhaseSpeciesStartIndex[iph];
        ThermoPhase& tp = thermo(iph);
        size_t nspPhase = tp.nSpecies();
        if (iph != solnPhase_) {
            for (size_t k = 0; k < nspPhase; k++) {
               heatCapacity += spMoles_final_[kStart + k] * CvPM_[kStart + k];
            }
        }
    }
    return heatCapacity;
}
//====================================================================================================================
double Electrode::SolidEnthalpy() const
{
    //
    //   Calculate C_ on the fly if it isn't part of the normal calculation.
    //         -> debatable about how to do this most efficiently.
    //
    if (!doThermalPropertyCalculations_) {
	for (size_t iph = 0; iph < (size_t) m_NumTotPhases; iph++) {
	    int istart = m_PhaseSpeciesStartIndex[iph];
	    ThermoPhase& tp = thermo(iph);
	    tp.getPartialMolarEnthalpies(&(enthalpyMolar_final_[istart]));
	}
    }

    double enthalpy = 0.0;
    for (size_t iph = 0; iph < m_NumVolPhases; iph++) {
        int kStart = m_PhaseSpeciesStartIndex[iph];
        if (iph != solnPhase_) {
            for (size_t k = 0; k < thermo(iph).nSpecies(); k++) {
               enthalpy += spMoles_final_[kStart + k] * enthalpyMolar_final_[kStart + k];
            }
        }
    }
    return enthalpy;
}
//====================================================================================================================
size_t Electrode::nSurfaces() const
{
    return numSurfaces_;
}
//====================================================================================================================
size_t Electrode::nReactions(size_t isk) const
{
    ReactingSurDomain* rsd = RSD_List_[isk];
    if (!rsd) {
        return (size_t) 0;
    } else {
        return rsd->rmcVector.size();
    }
}
//====================================================================================================================
void Electrode::getSurfaceAreas(double* const surfArea) const
{
    for (size_t i = 0; i < numSurfaces_; i++) {
        surfArea[i] = surfaceAreaRS_final_[i];
    }
}
//====================================================================================================================
// ----------------------------------- QUERY AND SET VOLTAGES --------------------------------------------------------
//====================================================================================================================
void Electrode::setVoltages(const double phiMetal, const double phiElectrolyte)
{
    phaseVoltages_[metalPhase_] = phiMetal;
    phaseVoltages_[solnPhase_] = phiElectrolyte;
    deltaVoltage_ = phaseVoltages_[metalPhase_] - phaseVoltages_[solnPhase_];
    updateState();
}
//====================================================================================================================
double Electrode::voltage() const
{
  return deltaVoltage_;
}
//====================================================================================================================
double Electrode::phaseElectricPotential(size_t iph) const
{
    return phaseVoltages_[iph];
}
//====================================================================================================================
void Electrode::setPhaseElectricPotential(size_t iph, double phi)
{
    phaseVoltages_[iph] = phi;
    updateState();
}
//====================================================================================================================
//   Take the state (i.e., the final state) within the Electrode object and push it down
//   to the ThermoPhase objects and other variables that are part of the Electrode object
/*
 *  (virtual function from Electrode)
 *  This virtual function should be written so that child routines are not called by parent routines.
 *  It may be the case that child routines will surplant parent routines, in order to create efficiency.
 *
 *  We take the values of spMoles_final_[], and the number of particles, and their size,
 *  which are the default specification of the state variables within the Electrode object,
 *  and propagate them down to the ThermoPhase objects in the electrode.
 *  We also calculate the volumetric properties of the Electrode, the phase moles,
 *  and the mole fractions, and the external radius of the particle
 *
 *  All of these properties are defined for the _final_ state.
 *
 *  Thus, this is the main routine that reconciles all of the state information within the object.
 *  At the end of this routine, all aspects of the final state are consistent with each other.
 *
 *  prerequisites: The object must have been already been created.
 *
 *
 *  Summary: State Variables
 *              spMoles_final_[kGlobal]  : kGlobal in solid phase species
 *              particleNumberToFollow_
 *
 *           Independent Variables
 *              spMoles_final_[kGlobal]  :   kGlobal in electrolyte phase
 *              phaseVoltages_[iph] :        iph in solid phase electrodes
 *              Temperature
 *              Pressure
 *
 *   Dependent StateVariables
 *              phaseMoles_final_[iph]       iph in solid phase electrode
 *              spMf_final_[kGlobal]         kGlobal in solid phase species
 *              VolPM_[kGlobal]              kGlobal in solid phase species
 *              spElectroChemPot_[kGlobal]   kGlobal in solid phase species
 *              ThermoPhase[iph]             iph in solid phase electrode
 *        	    phaseMolarVolumes_[iph]      iph in solid phase electrode
 *              ElectrodeSolidVolume_        (this is called from SolidVolume();
 *              Radius_exterior_final_       (this is calculated from ElectrodeSolidVolume_)
 *
 */
void Electrode::updateState()
{
    //   Loop over all phases in the object, volume and surface phases
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        updateState_Phase(iph);
    }
    deltaVoltage_ = phaseVoltages_[metalPhase_] - phaseVoltages_[solnPhase_];
    
    //    Calculate the volume of the electrode phase. This is the main routine to do this.
    ElectrodeSolidVolume_ = SolidVol();

    //  Calculate the external radius of the particle (in meters) assuming that all particles are the same
    //  size and are spherical
    double volParticle = ElectrodeSolidVolume_ / particleNumberToFollow_;
    Radius_exterior_final_ = pow(volParticle * 3.0 / (4.0 * Pi), 0.3333333333333333);
}
//====================================================================================================================
//
/*
 *  Note this version doesn't call child member functions
 */
void Electrode::updateState_OnionOut()
{
    /*
     * However, I want to make sure mole fractions are
     * consistent with final moles.
     */
    for (size_t i = 0; i < m_NumTotPhases; i++) {
        Electrode::updateState_Phase(i);
    }
    /*
     * Loop over all phases in the object
     */
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        ThermoPhase* tphase = &(thermo(iph));
        std::string pName = tphase->id();
        size_t kStart = m_PhaseSpeciesStartIndex[iph];
        /*
         * A direct call to tphase to set the state is unnecessary,
         * because m_mp will do it anyway.
         */
        tphase->setState_TPX(temperature_, pressure_, &spMf_final_[kStart]);
        tphase->setElectricPotential(phaseVoltages_[iph]);

        /*
         * Ok, we have set the state. Now upload the Vol and mu's.
         */
        tphase->getPartialMolarVolumes(&(VolPM_[kStart]));
        tphase->getElectrochemPotentials(&(spElectroChemPot_[kStart]));

        if (iph < m_NumVolPhases) {
            phaseMolarVolumes_[iph] = tphase->molarVolume();
        } else {
            phaseMolarVolumes_[iph] = 0.0;
        }

    }
    /*
     * Calculate the volume of the electrode phase. This is the main routine to do this.
     */
    ElectrodeSolidVolume_ = Electrode::SolidVol();

    double vol = ElectrodeSolidVolume_ / particleNumberToFollow_;

    /*
     *  Calculate the external radius of the particle (in meters) assuming that all particles are the same
     *  size and are spherical
     */
    Radius_exterior_final_ = pow(vol * 3.0 / (4.0 * Pi), 0.3333333333333333);
}
//====================================================================================================================
//  Recalcualte the surface areas of the surfaces for the final state
/*
 * (virtual function from Electrode)
 *
 *    Here, we don't know anything about the morphology of the particle. Models that do know will override this
 *    function.
 *    Hwere we assume that the surface area is equal to the exterior surface area multiplied by the number
 *    of particles.
 *    We also respect keeping the surface area as zero, if it is initially set to zero.
 */
void Electrode::updateSurfaceAreas()
{
    double totalSurfArea = 4.0 * Pi * Radius_exterior_final_ * Radius_exterior_final_ * particleNumberToFollow_;
    for (size_t i = 0; i < numSurfaces_; i++) {
        if (surfaceAreaRS_final_[i] != 0.0) {
            surfaceAreaRS_final_[i] = totalSurfArea;
        }
    }
    /*
     * Here we assume that the surface just creates species enough to fill the available surface sites
     */
    int i = 0;
    for (size_t iph = m_NumVolPhases; iph < m_NumTotPhases; iph++, i++) {
        ThermoPhase* tphase = &(thermo(iph));
        sphaseMolarAreas_[i] = tphase->molarArea();
        int kstart = m_PhaseSpeciesStartIndex[iph];
        int nsp = m_PhaseSpeciesStartIndex[iph + 1] - kstart;
        phaseMoles_final_[iph] = surfaceAreaRS_final_[i] / sphaseMolarAreas_[i];
        for (int k = 0; k < nsp; k++) {
            spMoles_final_[k + kstart] = spMf_final_[k] * phaseMoles_final_[iph];
        }
    }
}
//====================================================================================================================
//   This is used to set the phase information that is implicit but not set by a restart or an initialization
/*
 *  (virtual function from Electrode)
 *
 *  This is called immediately after the restart file's contents are loaded into the electrode object.
 *  We then call this function to calculate the internal flags. Then, we call updateState() to make sure
 *  all information about the state is self-consistent.
 *
 *  @param flagErrors If true any changes in the current flags caused by a mismatch between the state
 *                    and the values of the flags will cause an error exit.
 */
bool Electrode::stateToPhaseFlagsReconciliation(bool flagErrors)
{
    /*
     *  We probably need to set the phase existance for Reacting Surfaces.
     */
    setPhaseExistenceForReactingSurfaces(false);
    return false;
}
//====================================================================================================================
//   Reactant stoichiometric coefficient
/*
 * Get the reactant stoichiometric coefficient for the kth global species
 * in the ith reaction of the reacting surface domain with index isk.
 */
double Electrode::reactantStoichCoeff(const size_t isk, size_t kGlobal, size_t i) const
{
    ReactingSurDomain* rsd = RSD_List_[isk];
    size_t krsd = rsd->PLtoKinSpeciesIndex_[kGlobal];
    if (krsd == npos) {
        return 0.0;
    }
    double rst = rsd->reactantStoichCoeff(krsd, i);
    return rst;
}
//====================================================================================================================
//   Reactant stoichiometric coefficient
/*
 * Get the reactant stoichiometric coefficient for the kth global species
 * in the ith reaction of the reacting surface domain with index isk.
 */
double Electrode::productStoichCoeff(const size_t isk, size_t kGlobal, size_t i) const
{
    ReactingSurDomain* rsd = RSD_List_[isk];
    size_t krsd = rsd->PLtoKinSpeciesIndex_[kGlobal];
    if (krsd == npos) {
        return 0.0;
    }
    double rst = rsd->productStoichCoeff(krsd, i);
    return rst;
}
//====================================================================================================================
// Get the net production rates of all species in the electrode object at the current conditions
/*
 *  This routine assumes that the underlying objects have been updated
 */
void Electrode::getNetSurfaceProductionRates(const size_t isk, double* const net) const
{
    std::fill_n(net, m_NumTotSpecies, 0.);
    /*
     *  This routine basically translates between species lists for the reacting surface
     *  domain and the Electrode.
     *  Later, when we have more than one reacting surface domain in the electrode object,
     *  this will do a lot more
     */

    /*
     *  For each Reacting surface
     */
    if (ActiveKineticsSurf_[isk]) {
        /*
         *  Get the species production rates for the reacting surface
         */
        const std::vector<double>& rsSpeciesProductionRates = RSD_List_[isk]->calcNetSurfaceProductionRateDensities();

        /*
         *  loop over the phases in the reacting surface
         */
        size_t nphRS = RSD_List_[isk]->nPhases();
        size_t kIndexKin = 0;
        for (size_t kph = 0; kph < nphRS; kph++) {
            size_t jph = RSD_List_[isk]->kinOrder[kph];
            size_t istart = m_PhaseSpeciesStartIndex[jph];
            size_t nsp = m_PhaseSpeciesStartIndex[jph + 1] - istart;
            for (size_t k = 0; k < nsp; k++) {
                net[istart + k] += rsSpeciesProductionRates[kIndexKin];
                kIndexKin++;
            }
        }
    }
}
//==================================================================================================================================
void Electrode::getIntegratedSpeciesProductionRates(double* const net) const
{
    if (pendingIntegratedStep_ != 1) {
        throw Electrode_Error("Electrode::integratedSpeciesProductionRates()",
                              "No integration step has been carried out");
    }
    double invDelT = 1.0;
    if (t_final_final_ > t_init_init_) {
        invDelT = 1.0 / (t_final_final_ - t_init_init_);
    }
    for (size_t k = 0; k < m_NumTotSpecies; k++) {
        net[k] = invDelT * spMoleIntegratedSourceTerm_[k];
    }
}
//==================================================================================================================================
//  Returns the current and the net production rates per unit surface area of surface isk
/*
 *    columb sec-1 m-2
 */
double Electrode::getNetSurfaceProductionRatesCurrent(const size_t isk, double* const net) const
{
    getNetSurfaceProductionRates(isk, net);
    // kmol sec-1 m-2
    double Eprod = net[kElectron_];
    // coulomb / kmol
    return Eprod * Faraday;
}
//==================================================================================================================================
//  Return the global species index of the electron in the Electrode object
size_t Electrode::kSpecElectron() const
{
    return kElectron_;
}
//==================================================================================================================================
double Electrode::getIntegratedProductionRatesCurrent(double* const net) const
{
    if (pendingIntegratedStep_ != 1) {
        double* netOne = &deltaG_[0];
        zeroD(m_NumTotSpecies, net);
        for (size_t isk = 0; isk < numSurfaces_; isk++) {
            if (ActiveKineticsSurf_[isk]) {
                getNetSurfaceProductionRates(isk, netOne);
                for (size_t k = 0; k < m_NumTotSpecies; k++) {
                    net[k] += netOne[k] * surfaceAreaRS_final_[isk];
                }
            }
        }
        double Eprod = net[kElectron_];
        return Eprod * Faraday;
    }
    getIntegratedSpeciesProductionRates(net);
    // kmol sec-1
    double Eprod = net[kElectron_];
    // coulomb / kmol
    return Eprod * Faraday;
}
//==================================================================================================================================
double Electrode::integratedCurrent() const
{
    if (pendingIntegratedStep_ != 1) {
        std::vector<double> net(m_NumTotSpecies, 0.0);
        double* netOne = &deltaG_[0];
        for (size_t isk = 0; isk < numSurfaces_; isk++) {
            if (ActiveKineticsSurf_[isk]) {
                getNetSurfaceProductionRates(isk, netOne);
                for (size_t k = 0; k < m_NumTotSpecies; k++) {
                    net[k] += netOne[k] * surfaceAreaRS_final_[isk];
                }
            }
        }
        double Eprod = net[kElectron_];
        // coulomb / kmol
        return Eprod * Faraday;
    }
    double invDelT = 1.0;
    if (t_final_final_ > t_init_init_) {
        invDelT = 1.0 / (t_final_final_ - t_init_init_);
    }
    // kmol sec-1
    double Eprod = invDelT * spMoleIntegratedSourceTerm_[kElectron_];
    // coulomb / kmol
    return Eprod * Faraday;
}
//====================================================================================================================
double Electrode::integratedLocalCurrent() const
{
    if (pendingIntegratedStep_ != 1) {
        std::vector<double> net(m_NumTotSpecies, 0.0);
        double* netOne = &deltaG_[0];
        for (size_t isk = 0; isk < numSurfaces_; isk++) {
            if (ActiveKineticsSurf_[isk]) {
                getNetSurfaceProductionRates(isk, netOne);
                for (size_t k = 0; k < m_NumTotSpecies; k++) {
                    net[k] += netOne[k] * surfaceAreaRS_final_[isk];
                }
            }
        }
        double Eprod = net[kElectron_];
        // coulomb / kmol
        return Eprod * Faraday;
    }
    double invDelT = 1.0;
    if (t_final_final_ > t_init_init_) {
        invDelT = 1.0 / (tfinal_ - tinit_);
    }
    // kmol sec-1
    double Eprod = invDelT * spMoleIntegratedSourceTermLast_[kElectron_];
    // returns coulomb / sec
    return Eprod * Faraday;
}
//====================================================================================================================
//  Returns the current and the net production rates of the phases in kg/m2/s from a single surface
/*
 *  Returns the net production rates of all phases from reactions on a single surface
 *
 *  @param isk Surface ID to get the fluxes from.
 *  @param phaseMassFlux  Returns the mass fluxes of the phases
 */
void Electrode::getPhaseMassFlux(const size_t isk, double* const phaseMassFlux)
{
    std::fill_n(phaseMassFlux, m_NumTotPhases, 0.);

    if (ActiveKineticsSurf_[isk]) {
        const std::vector<double>& rsSpeciesProductionRates = RSD_List_[isk]->calcNetSurfaceProductionRateDensities();

        size_t nphRS = RSD_List_[isk]->nPhases();
        size_t kIndexKin = 0;
        for (size_t kph = 0; kph < nphRS; kph++) {
            size_t jph = RSD_List_[isk]->kinOrder[kph];
            ThermoPhase& tp = thermo(jph);
            size_t istart = m_PhaseSpeciesStartIndex[jph];
            size_t nsp = m_PhaseSpeciesStartIndex[jph + 1] - istart;
            for (size_t k = 0; k < nsp; k++) {
                double net = rsSpeciesProductionRates[kIndexKin];
                double mw = tp.molecularWeight(k);
                phaseMassFlux[jph] += net * mw;
                kIndexKin++;
            }
        }
    }
}
//================================================================================================
//  Returns the current and the net production rates of the phases in kmol/m2/s from a single surface
/*
 *  Returns the net production rates of all phases from reactions on a single surface
 *
 *  @param isk Surface ID to get the fluxes from.
 *  @param phaseMassFlux  Returns the mass fluxes of the phases
 */
void Electrode::getPhaseMoleFlux(const size_t isk, double* const phaseMoleFlux)
{
    std::fill_n(phaseMoleFlux, m_NumTotPhases, 0.);
    if (ActiveKineticsSurf_[isk]) {
        const std::vector<double>& rsSpeciesProductionRates = RSD_List_[isk]->calcNetSurfaceProductionRateDensities();
        size_t nphRS = RSD_List_[isk]->nPhases();
        size_t kIndexKin = 0;
        for (size_t RSph = 0; RSph < nphRS; RSph++) {
            size_t PLph = RSD_List_[isk]->kinOrder[RSph];
            size_t istart = m_PhaseSpeciesStartIndex[PLph];
            size_t nsp = m_PhaseSpeciesStartIndex[PLph + 1] - istart;
            for (size_t k = 0; k < nsp; k++) {
                double net = rsSpeciesProductionRates[kIndexKin];
                phaseMoleFlux[PLph] += net;
                kIndexKin++;
            }
        }
    }
}
//================================================================================================
void Electrode::getPhaseProductionRates(const double* const speciesProductionRates,
					double* const phaseProductionRates) const
{
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        phaseProductionRates[iph] = 0.0;
        ThermoPhase& tp = thermo(iph);
        //string pname = tp.id();
        size_t istart = m_PhaseSpeciesStartIndex[iph];
        size_t nsp = tp.nSpecies();
        for (size_t ik = 0; ik < nsp; ik++) {
            size_t k = istart + ik;
            phaseProductionRates[iph] += speciesProductionRates[k];
        } 
    }
}
//================================================================================================
void Electrode::getIntegratedPhaseMoleTransfer(double* const phaseMolesTransfered)
{
    if (!pendingIntegratedStep_) {
        throw Electrode_Error(" Electrode::getIntegratedPhaseMoleTransfer", "no pending integration step");
    }
    double sum = 0.0;
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        phaseMolesTransfered[iph] = 0.0;
        ThermoPhase& tp = thermo(iph);
        //string pname = tp.id();
        size_t istart = m_PhaseSpeciesStartIndex[iph];
        size_t nsp = tp.nSpecies();
        for (size_t ik = 0; ik < nsp; ik++) {
            size_t k = istart + ik;
            phaseMolesTransfered[iph] += spMoleIntegratedSourceTerm_[k];
        }
        sum += fabs(phaseMolesTransfered[iph]);
    }
    /*
     *  TO DO : We can put in a check here on phaseMoles_Final vs phaseMoles_init_
     */
    /* THIS DOESNT WORK FOR ELECTRON PHASE and perhaps toher phases */
    /*
     double deltaPhaseMolesIntegrated = 0.0;

     sum += 1.0E-14;
     for (int iph = 0; iph < m_NumTotPhases; iph++) {
     deltaPhaseMolesIntegrated = phaseMoles_final_[iph] - phaseMoles_init_init_[iph];
     if (fabs(deltaPhaseMolesIntegrated - phaseMolesTransfered[iph])/ sum > 1.0E-7) {
     throw Electrode_Error(" Electrode::getIntegratedPhaseMoleTransfer",
     " mismatched phase deltas");
     }
     }
     */
}
//=================================================================================================
// Returns the integrated thermal energy source term (Joules)
/*
 *    Returns the heat release that has occurred during the global time step. 
 *
 *  @param Returns the heat release (joules)
 */
double Electrode::getIntegratedThermalEnergySourceTerm()
{
  return integratedThermalEnergySourceTerm_;
}
//=================================================================================================
//  Remember this is in joules
double Electrode::getIntegratedThermalEnergySourceTerm_overpotential()
{
  return integratedThermalEnergySourceTerm_overpotential_;
}
//=================================================================================================
//  Remember this is in joules
double Electrode::getIntegratedThermalEnergySourceTerm_reversibleEntropy()
{
  return integratedThermalEnergySourceTerm_reversibleEntropy_;
}
//=================================================================================================
// Returns the standard state OCV for the selected ReactingSurfaceDomain.
/*
 *  When there is more than a single reaction,
 *  pick open circuit potential for reaction that is
 *  closest to equilibrium given the cell voltage since this one
 *  is the one for which open circuit is most relevant.
 *  Note that it will be possible for the standard state OCV
 *  to be computed for a different reaction relative to the
 *  method openCircuitVoltage(isk) that computes OCV from
 *  the current concentrations because it looks for the reaction
 *  closest to equilibrium given the current cell voltage.
 */
double Electrode::openCircuitVoltageSSRxn(size_t isk, size_t iReaction) const
{
    ReactingSurDomain* rsd = RSD_List_[isk];
    if (!rsd) {
        return 0.0;
    }
    size_t rxnIndex = 0;
    if (iReaction != npos) {
        rxnIndex = iReaction;
    }
    rsd->getDeltaSSGibbs(DATA_PTR(deltaG_));
    size_t length = rsd->rmcVector.size();
    double ERxn = 0.0;

    size_t metalPhaseRS = rsd->PLtoKinPhaseIndex_[metalPhase_];
    if (metalPhaseRS != npos) {
        RxnMolChange* rmc = rsd->rmcVector[rxnIndex];
        /*
         *  Compute open circuit potential from first reaction
         */
        double nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseRS];
        if (nStoichElectrons != 0.0) {
            ERxn = deltaG_[rxnIndex] / Faraday / nStoichElectrons;
        }
        if (iReaction != npos) {
            return ERxn;
        }

        /*
         *  Compute open circuit potential from other reactions.
         *  Pick open circuit potential for reaction that is
         *  closest to equilibrium given cell voltage since this one
         *  is the one for which open circuit is most relevant.
         */
        if (length > 1) {
            for (size_t i = 1; i < length; i++) {
                rmc = rsd->rmcVector[i];
                nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseRS];
                if (nStoichElectrons != 0.0) {
                    if (fabs(deltaVoltage_ - ERxn) > fabs(deltaVoltage_ - deltaG_[i] / Faraday / nStoichElectrons)) {
                        rxnIndex = i;
                        ERxn = deltaG_[i] / Faraday / nStoichElectrons;
                    }
                }
            }
        }
        return ERxn;
    }
    return 0.0;
}
//==================================================================================================================================
// Returns the equilibrium OCV for the selected ReactingSurfaceDomain and current conditions.
/*
 *  When there is more than a single reaction,
 *  pick open circuit potential for reaction that is
 *  closest to equilibrium given the cell voltage since this one
 *  is the one for which open circuit is most relevant.
 */
double Electrode::openCircuitVoltageRxn(size_t isk, size_t iReaction, bool comparedToReferenceElectrode) const
{
    ReactingSurDomain* rsd = RSD_List_[isk];
    if (!rsd) {
        return 0.0;
    }
    size_t rxnIndex = 0;
    if (iReaction != npos) {
        rxnIndex = iReaction;
    }

    if (comparedToReferenceElectrode) {
	rsd->getDeltaGibbs_electrolyteSS(DATA_PTR(deltaG_));
    } else {
	rsd->getDeltaGibbs(DATA_PTR(deltaG_));
    }

    size_t nR = rsd->nReactions();
    double ERxn = 0.0; // store open circuit voltage here

    size_t metalPhaseRS = rsd->PLtoKinPhaseIndex_[metalPhase_];
    if (metalPhaseRS != npos) {

        RxnMolChange* rmc = rsd->rmcVector[rxnIndex];
        /*
         *  Compute open circuit potential from first reaction
         */
        double nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseRS];
        if (nStoichElectrons != 0.0) {
            ERxn = deltaG_[rxnIndex] / Faraday / nStoichElectrons;
        }
  
        if (iReaction != npos) {
            return ERxn;
        }
        /*
         *  Compute open circuit potential from other reactions.  Pick open circuit potential for reaction that is
         *  closest to equilibrium given cell voltage since this one is the one for which open circuit is most relevant.
         */
        if (nR > 1) {
            for (size_t i = 1; i < nR; i++) {
                rmc = rsd->rmcVector[i];
                nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseRS];
                if (nStoichElectrons != 0.0) {
                    if (fabs(deltaVoltage_ - ERxn) > fabs(deltaVoltage_ - deltaG_[i] / Faraday / nStoichElectrons)) {
                        rxnIndex = i;
                        ERxn = deltaG_[i] / Faraday / nStoichElectrons;
                    }
                }
            }
        }
        return ERxn;
    }
    return 0.0;
}
//==================================================================================================================================
// Returns the equilibrium OCV for all reactions on the selected ReactingSurfaceDomain at the current conditions.
/*
 *  When there is more than a single reaction,
 *  there may be many open circuit voltages, one for each for reaction.
 */
void Electrode::getOpenCircuitVoltages(size_t isk, double* Erxn, bool comparedToReferenceElectrode) const
{
    ReactingSurDomain* rsd = RSD_List_[isk];
    if (!rsd) {
        Erxn[0] = 0.0;
        return;
    }
    /*
     *  get the reaction Gibbs free energy for all reactions defined on the reacting surface 
     */
    if (comparedToReferenceElectrode) {
	rsd->getDeltaGibbs_electrolyteSS(DATA_PTR(deltaG_));
    } else {
	rsd->getDeltaGibbs(DATA_PTR(deltaG_));
    }
    // find number of reactions
    size_t nr = rsd->rmcVector.size();
    size_t metalPhaseRS = rsd->PLtoKinPhaseIndex_[metalPhase_];
    if (metalPhaseRS != npos) {
        RxnMolChange* rmc;
        double nStoichElectrons;
        for (size_t i = 0; i < nr; i++) {
            rmc = rsd->rmcVector[i];
            // 
            //  nStoichElectrons are the net number of product electrons
            //
            nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseRS];
            if (nStoichElectrons != 0.0) {
                Erxn[i] = deltaG_[i] / Faraday / nStoichElectrons;
            } else {
                Erxn[i] = 0.0;
            }
        }

    } else {
        // no metalPhase_
        for (size_t i = 0; i < nr; i++) {
            Erxn[0] = 0.0;
        }
    }
    return;
}
//======================================================================================================
double Electrode::openCircuitVoltage_MixtureAveraged(size_t isk,  bool comparedToReferenceElectrode)
{
    return Electrode::openCircuitVoltage(isk, comparedToReferenceElectrode);
}
//==================================================================================================================================
// A calculation of the open circuit voltage of an Electrode object
/*
 *  This routine uses a root finder to find the voltage at which there is zero net electron production.  
 *  It leaves the object unchanged. However, it does change the voltage of the phases during the calculation,
 *  so this is a nonconst function.
 */
double Electrode::openCircuitVoltage(size_t isk, bool comparedToReferenceElectrode)
{
#ifdef DEBUG_OCV
    static int s_numIters = 0;
    s_numIters++;
#endif
    int printDebug = 30;
    static int oIts = 0;
    double nStoichElectrons = 0.0;
    double deltaT = 0.005;
    size_t rxnIndex = 0;
    ReactingSurDomain* rsd = RSD_List_[isk];
    if (!rsd) {
        return 0.0;
    }

  startOver: ;
    /*
     *  Get the current value of deltaG_[] for reactions defined on the current reacting surface
     */
    if (comparedToReferenceElectrode) {
	rsd->getDeltaGibbs_electrolyteSS(DATA_PTR(deltaG_));
    } else {
	rsd->getDeltaGibbs(DATA_PTR(deltaG_));
    }
    if (printDebug < 0) {
        printf("we are here, oIts = %d\n", oIts);
    }
    RxnMolChange* rmc = 0;
    size_t nR = rsd->nReactions();
    // If we don't have any reactions we can't do the calculation
    if (nR == 0) {
        return 0.0;
    }
    size_t nP = rsd->nPhases();
    double ERxn = 0.0; // store open circuit voltage here

    std::vector<int> phaseExistsInit(nP, 1);
    std::vector<int> phaseStabInit(nP, 1);

    for (size_t iph = 0; iph < nP; iph++) {
        phaseExistsInit[iph] = rsd->phaseStability(iph);
        phaseStabInit[iph] = rsd->phaseStability(iph);
        if (printLvl_ > printDebug) {
            ThermoPhase& tp = rsd->thermo(iph);
            std::string pname = tp.name();
            printf(" %20s %d  %d \n", pname.c_str(), phaseExistsInit[iph], phaseStabInit[iph]);
        }
    }

    // If we don't have a metal phase, we don't know how to do the calculation (lazy I think)
    size_t metalPhaseRS = rsd->PLtoKinPhaseIndex_[metalPhase_];
    if (metalPhaseRS == npos) {
        return 0.0;
    }
    double phiMin = +100.;
    double phiMax = -100.;
    int nT = 0;
    // Counter for the number of reactions producing or consuming electrons
    int numERxns = 0;
    std::vector<int> reactionTurnedOn(nR, 1);
    bool thereIsOneTurnedOn = false;
    //
    //  Our first check is to see whether there is one reaction occurring between phases which are allowed to exist
    //
    for (size_t i = 0; i < nR; i++) {
        size_t rxnIndex = i;
        rmc = rsd->rmcVector[rxnIndex];
        /*
         *  Compute open circuit potential for the current reaction reaction
         */
        nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseRS];
        if (nStoichElectrons != 0.0) {
	    numERxns++;
            for (size_t jph = 0; jph < nP; jph++) {
		//  \todo check to see whether this needs deltaG additions to see what direction it is going.
		// If this phase participates in the reaction mechanism ...
                if (rmc->m_phaseReactantMoles[jph] > 0.0 || rmc->m_phaseProductMoles[jph] > 0.0) {
		    // If this phase doesn't exist or if this phase is not stable then discount this reaction.
		    // \todo Check that this logic is right.
                    if (phaseExistsInit[jph] == 0 || phaseStabInit[jph] == 0) {
                        reactionTurnedOn[i] = 0;
                    }
                }
            }
        } else {
            reactionTurnedOn[i] = 0;
        }
	if (reactionTurnedOn[i]) {
	    nT++;
            ERxn = deltaG_[rxnIndex] / Faraday / nStoichElectrons;
            phiMax = std::max(phiMax, ERxn);
            phiMin = std::min(phiMin, ERxn);
	}

        if (printLvl_ > printDebug) {
            std::string rstring = rsd->reactionString(i);
            printf("%d %100.100s %g  \n", static_cast<int>(i), rstring.c_str(), ERxn);
        }
        if (reactionTurnedOn[i]) {
            thereIsOneTurnedOn = true;
        }
    }

    if (printLvl_ > printDebug) {
        printf("    deltaVoltage_ = %g\n", deltaVoltage_);
    }
    /*
     *  If there is only one participating electron reaction then the
     *  open circuit potential is given by the Nernst equation for that reaction.
     *  We just calculated it above. Therefore, we can do a quick return.
     */
    if (nT == 1) {
        return ERxn;
    }
    // 
    // Early exist if there aren't any reactions with electrons, we return 0.0
    //
    if (numERxns == 0) {
	return ERxn;
    }
    //
    //  ----------------- When we are here, we have a more complicated case ----------------------------------
    //
    if (comparedToReferenceElectrode) {
     	rsd->getDeltaGibbs(DATA_PTR(deltaG_));
    } 
    double phiMetalInit = phaseVoltages_[metalPhase_];
    deltaVoltage_ = phaseVoltages_[metalPhase_] - phaseVoltages_[solnPhase_];
    //
    //  Make copies of the phaseExistsInit and phaseStabInit vectors. We will turn phases on now
    //
    std::vector<int> phaseExists(phaseExistsInit);
    std::vector<int> phaseStab(phaseStabInit);

#ifdef DEBUG_OCV
    printf("s_numIters = %d\n", s_numIters);
    if (s_numIters == 389) {
	printf("we are here\n");
    }
#endif

    double ERxnBest = 1.0E300;
    double ERxnMin = 0.0;;
    double ERxnMax = 0.0;
    numERxns = 0;
    size_t iR = npos;
    for (size_t i = 0; i < nR; i++) {
        bool electReact = false;
        bool electProd = false;
        /*
         *  Get the extra information for the reaction
         */
        rmc = rsd->rmcVector[i];
        /*
         *  First figure out if electron is a reactant or a product
         */
        nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseRS];
        if (nStoichElectrons > 0) {
            electReact = true; 
        }
        if (nStoichElectrons < 0) {
            electProd = true;
        }
        if (nStoichElectrons != 0.0) {
	    numERxns++;
	    iR = i;
	    //
	    //  Calculate the open circuit voltage for the current reaction
	    //
            ERxn = deltaG_[i] / Faraday / nStoichElectrons;
	    if (ERxnBest == 1.0E300) {
		ERxnMin = ERxn;
		ERxnMax = ERxn;
	    }
	    if (ERxn < ERxnMin) {
		ERxnMin = ERxn;
	    }
	    if (ERxn > ERxnMax) {
		ERxnMax = ERxn;
	    }
	    //
	    //  Pick out the OCV closest to the current detla voltage
            if (fabs(deltaVoltage_ - ERxnBest) > fabs(deltaVoltage_ - ERxn)) {
                ERxnBest = ERxn;
            }
            if (!thereIsOneTurnedOn) {
                if (electProd) {
                    if (ERxn > deltaVoltage_) {
                        for (size_t jph = 0; jph < nP; jph++) {
                            if (rmc->m_phaseReactantMoles[jph] > 0.0 || rmc->m_phaseProductMoles[jph] > 0.0) {
                                phaseExists[jph] = 1;
                                phaseStab[jph] = 1;
                            }
                        }
                    }
                }
                if (electReact) {
                    if (ERxn <= deltaVoltage_) {
                        for (size_t jph = 0; jph < nP; jph++) {
                            if (rmc->m_phaseReactantMoles[jph] > 0.0 || rmc->m_phaseProductMoles[jph] > 0.0) {
                                phaseExists[jph] = 1;
                                phaseStab[jph] = 1;
                            }
                        }
                    }
                }
            }
        }
    }
    if (numERxns == 1) {
	if (comparedToReferenceElectrode) {
	    rsd->getDeltaGibbs_electrolyteSS(DATA_PTR(deltaG_)); 
	    ERxn = deltaG_[iR] / Faraday / nStoichElectrons;
	}
	return ERxn;
    } else {
	if (comparedToReferenceElectrode) {
#ifdef DEBUG_OCV
	    printf("openCircuitWarning() warning: Calculating OCV for multiple electron reaction not carried out correctly yet\n");
#endif
	}
    }
    
    if (!thereIsOneTurnedOn) {
        rmc = rsd->rmcVector[rxnIndex];
        for (size_t jph = 0; jph < nP; jph++) {
            if (rmc->m_phaseReactantMoles[jph] > 0.0 || rmc->m_phaseProductMoles[jph] > 0.0) {
                phaseExists[jph] = 1;
                phaseStab[jph] = 1;
            }
        }
    }

    for (size_t iph = 0; iph < nP; iph++) {
        rsd->setPhaseExistence(iph, phaseExists[iph]);
        rsd->setPhaseStability(iph, phaseStab[iph]);
        if (printLvl_ > printDebug) {
            ThermoPhase& tp = rsd->thermo(iph);
            std::string pname = tp.name();
            printf(" %20s %d  %d \n", pname.c_str(), phaseExists[iph], phaseStab[iph]);
        }
    }

    /*
     *  bump the values
     */
    phiMin = phiMin - fabs(0.001 + fabs((phiMax - phiMin) * 0.01));
    phiMax = phiMax + fabs(0.001 + fabs((phiMax - phiMin) * 0.01));

    RSD_ElectronProduction ep(rsd, metalPhase_, deltaT);

    if (printLvl_ > printDebug) {
        ep.printLvl_ = 15;
    }

    RootFind rf(&ep);

    ERxn = openCircuitVoltageSSRxn(isk);
    ERxn = ERxnBest;
    double phiMetalRxn = ERxn    + phaseVoltages_[solnPhase_];
    double phiMetalMax = ERxnMax + phaseVoltages_[solnPhase_] + 0.01;
    double phiMetalMin = ERxnMin + phaseVoltages_[solnPhase_] - 0.01;
#ifdef DEBUG_OCV
    rf.setPrintLvl(15);
#else
    if (printLvl_ > printDebug) {
        rf.setPrintLvl(15);
    }
#endif
    rf.setTol(1.0E-5, 1.0E-10);
    rf.setFuncIsGenerallyIncreasing(true);
    rf.setDeltaX(0.01);

    double funcTargetValue = 0.0;

    int retn = rf.solve(phiMetalMin, phiMetalMax, 100, funcTargetValue, &phiMetalRxn);

    if (retn != 0) {
        printf(" Electrode::openCircuitVoltage: rootSolver returned an error = %d\n", retn);
        printf("            phiMetalMin = %g, phiMetalMax = %g, phiMetalRxn = %g, funcTargetValue = %g\n", phiMetalMin,
               phiMetalMax, phiMetalRxn, funcTargetValue);
        printf("            Redoing calculating with printing turned on, then terminating, oIts = %d\n", oIts);
#ifdef DEBUG_OCV
        printf("            OCV iteration# = %d\n", s_numIters);
#endif
        if (printDebug >= 0) {
            printDebug = -1;
            goto startOver;
        }
        exit(-1);
    }
    phaseVoltages_[metalPhase_] = phiMetalInit;
    ERxn = phiMetalRxn - phaseVoltages_[solnPhase_];

    for (size_t iph = 0; iph < rmc->m_nPhases; iph++) {
        rsd->setPhaseExistence(iph, phaseExistsInit[iph]);
        rsd->setPhaseStability(iph, phaseStabInit[iph]);
        if (printLvl_ > printDebug) {
            ThermoPhase& tp = rsd->thermo(iph);
            std::string pname = tp.name();
            printf(" %20s %d  %d \n", pname.c_str(), phaseExistsInit[iph], phaseStabInit[iph]);
        }
    }

    updateState();

    oIts++;

    return ERxn;
}
//====================================================================================================================
// Returns the overpotential for the current conditions
/*
 *  The overpotential is the current voltage minus the open circuit voltage.
 */
double Electrode::overpotential(size_t isk)
{
    double Erxn = openCircuitVoltage(isk);
    return (deltaVoltage_ - Erxn);
}
//====================================================================================================================
double Electrode::overpotentialRxn(size_t isk, size_t irxn)
{
    double Erxn = openCircuitVoltageRxn(isk, irxn);
    return (deltaVoltage_ - Erxn);
}
//==================================================================================================================================
double Electrode::openCircuitVoltageTotal(bool comparedToReferenceElectrode)
{
    if (numSurfaces_ == 1) {
        return openCircuitVoltage(0, comparedToReferenceElectrode);
    }
    int numA = 0;
    double OCV = 0.0;
    for (size_t iSurf = 0; iSurf < numSurfaces_; iSurf++) {
        if (ActiveKineticsSurf_[iSurf]) {
            numA++;
            OCV = openCircuitVoltage(0, comparedToReferenceElectrode);
        }
    } 
    // prelim implementation -> Doesn't work in complicated situations
    return OCV;
}
//==================================================================================================================================
double Electrode::openCircuitVoltageTotal_MixtureAveraged(bool comparedToReferenceElectrode)
{
    if (numSurfaces_ == 1) {
        return openCircuitVoltage_MixtureAveraged(0, comparedToReferenceElectrode);
    }
    int numA = 0;
    double OCV = 0.0;
    for (size_t iSurf = 0; iSurf < numSurfaces_; iSurf++) {
        if (ActiveKineticsSurf_[iSurf]) {
            numA++;
            OCV = openCircuitVoltage(0, comparedToReferenceElectrode);
        }
    } 
    // prelim implementation -> Doesn't work in complicated situations
    return OCV;
}
//==================================================================================================================================
size_t Electrode::kKinSpecElectron(size_t isurf) const
{
    return kKinSpecElectron_sph_[isurf];
}
//==================================================================================================================================
size_t Electrode::metalPhaseIndex() const
{
    return metalPhase_;
}
//==================================================================================================================================
size_t Electrode::solnPhaseIndex() const
{
    return solnPhase_;
}
//==================================================================================================================================
size_t Electrode::numSolnPhaseSpecies() const
{
    return phasePtr(phase_name(solnPhase_).c_str())->nSpecies();
}
//==================================================================================================================================
double Electrode::polarizationAnalysisSurf(std::vector<PolarizationSurfRxnResults>& psr_list)
{
  // Don't compare to reference electrode
  const bool comparedToReferenceElectrode = false;

  psr_list.clear();
 
  // Whatever path we take, the volts taken must be the same.
  double volts = voltage();

  // Create a PolarizationSurfResults record for each active surface in the problem
  /*
   *  We will only look at the exit conditions. Then, we'll normalize by the total current over the interval
   *  The base class doesn't have many polarization modes accessible
   */
  for (size_t iSurf = 0; iSurf < numSurfaces_; iSurf++) {
        if (ActiveKineticsSurf_[iSurf]) {
            ReactingSurDomain* rsd = RSD_List_[iSurf];
             
            // Calculate the non mixture averaged OCV of the surface. This is the OCV at which no electrons are produced 
            // given the current state of the bulk phases next to it and given the current surface adsorbate state.
            double ocvSurf = openCircuitVoltage(iSurf, comparedToReferenceElectrode);

            size_t nr = rsd->nReactions();

            bool eok;
            doublevalue nStoich, OCV, io, overPotential, beta, resistance;
            size_t numErxn = 0;
            for (size_t iRxn = 0; iRxn < nr; iRxn++) {
                eok = rsd->getExchangeCurrentDensityFormulation(iRxn, nStoich, OCV, io, overPotential, beta, resistance);
                double ocvSurfRxn = OCV; 
                double ocvSurfRxn_local = ocvSurfRxn;
                if (eok) {
                    if (nStoich != 0.0) {
                        double icurrPerArea = rsd->calcCurrentDensity(overPotential, nStoich, OCV, beta, temperature_, resistance);

                        // Create a record
                        psr_list.emplace_back(electrodeDomainNumber_, electrodeCellNumber_, iSurf, iRxn);
                        PolarizationSurfRxnResults& psr = psr_list.back();
                        psr.icurrSurf = 0.0;
                        numErxn++;
                        if (resistance != 0.0) {
                           VoltPolPhenom vpp;
               
                           double voltsRes = psr.icurrSurf * resistance;
                           vpp.ipolType = RESISTANCE_FILM_PL;
                           vpp.voltageDrop = voltsRes;
                           psr.voltsPol_list.push_back(vpp);
                           ocvSurfRxn_local -= voltsRes;
                           overPotential -= voltsRes;
                        }

                        VoltPolPhenom vpp;
                        vpp.ipolType = SURFACE_OVERPOTENTIAL_PL;
                        vpp.voltageDrop = overPotential; 
                        psr.voltsPol_list.push_back(vpp);
                        
                        psr.VoltageElectrode = volts;
                        psr.VoltageTotal = volts;
                        psr.ocvSurf = ocvSurf;
                        psr.ocvSurfRxn = ocvSurfRxn_local;
                        psr.icurrSurf = 0.5 * (surfaceAreaRS_final_final_[iSurf] + surfaceAreaRS_init_init_[iSurf]) * icurrPerArea;
                    }
                }
             }
             if (numErxn > 1) {
                 throw Electrode_Error("", "complicated situation not handled yet.");
             }
             
        }
    }

    // Now normalize it by the actual value of the total integrated current.
    /*
     *  Below we do a lot of cheating. We want the total integrated current through all channels allowed by the Electrode
     *  object. Currently, we are not storring this information in all its details. Therefore, we assume that the 
     *  conditions at the end of the global time step can be used throughout the global time step, but using a 
     *  constant factor to get the integrated current correct throughout the time-step.
     */
    double intCurrentTotal = integratedCurrent();
    if (fabs(intCurrentTotal) < 1.0e-200) {
        for (size_t iRec = 0; iRec < psr_list.size(); iRec++) {
            PolarizationSurfRxnResults& psr = psr_list[iRec];
            psr.icurrSurf = 0.0;
        }
    } else {
        if (psr_list.size() > 0) {
            double icurrSum = 0.0;
            size_t iRecMax = 0;
            double iRmax = 0.0;
            for (size_t iRec = 0; iRec < psr_list.size(); iRec++) {
                PolarizationSurfRxnResults& psr = psr_list[iRec];
                icurrSum += psr.icurrSurf;
                if (fabs(psr.icurrSurf) > iRmax) {
                    iRecMax = iRec;
                    iRmax = fabs(psr.icurrSurf);
                }
            }
            if (fabs(icurrSum) < 1.0E-100) {
               PolarizationSurfRxnResults& psr = psr_list[iRecMax];
               psr.icurrSurf += (intCurrentTotal - icurrSum);
            } else {
                double factor = intCurrentTotal / icurrSum;
                for (size_t iRec = 0; iRec < psr_list.size(); iRec++) {
                    PolarizationSurfRxnResults& psr = psr_list[iRec];
                    psr.icurrSurf *= factor;
                }
            }
        } else {
            throw Electrode_Error("Electrode::polarizationAnalysSurf()", "Current is zero but we didn't pick up any modes");
        }
    }

    return intCurrentTotal;
}
//==================================================================================================================================
// Return the number of extra print tables
int Electrode::getNumPrintTables() const
{
    return 0;
}
//==================================================================================================================================
//! Get the values that are printed in tables for the 1D code.
void Electrode::getPrintTable(int itable, std::vector<std::string>& colNames, std::vector<double>& colValues) const
{
}
//==================================================================================================================================
// Set the current depth of discharge
/*
 * This is roughly equal to the total number of electrons that has been discharged from a fully charged state
 * divided by the total moles
 *
 * @return  returns the depth of discharge in amps sec
 * @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
 */
void Electrode::setRelativeCapacityDischargedPerMole(double relDischargedPerMole, int platNum)
{
    if (platNum != -1) {
        throw Electrode_Error(" Electrode::setRelativeCapacityDischargedPerMole()", "platNum not -1");
    }
    if (pendingIntegratedStep_) {
        throw Electrode_Error(" Electrode::setRelativeCapacityDischargedPerMole()", "pending step");
    }
    throw Electrode_Error(" Electrode::setRelativeCapacityDischargedPerMole()", "not implemented");
}
//==================================================================================================================================
RxnMolChange* Electrode::rxnMolChangesEGR(size_t iegr)
{
    return m_rmcEGR[iegr];
}
//==================================================================================================================================
void Electrode::addExtraGlobalRxn(const EGRInput& egri)
{
    size_t rsdI = egri.m_RSD_index;
    InterfaceKinetics* iKA = RSD_List_[rsdI];
    if (!iKA) {
       throw Electrode_Error("Electrode::addExtraGlobalRxn()", "No kinetics object for reacting surface domain " + int2str(rsdI));
    }
    size_t nReactionsA = iKA->nReactions();

    ZZCantera::ExtraGlobalRxn* egr = new ExtraGlobalRxn(*iKA);
    double* RxnVector = new double[nReactionsA];
    for (size_t i = 0; i < nReactionsA; i++) {
        RxnVector[i] = 0.0;
    }

    for (size_t iErxn = 0; iErxn < (size_t) egri.m_numElemReactions; iErxn++) {
        ERSSpec* ers_ptr = egri.m_ERSList[iErxn];
        RxnVector[ers_ptr->m_reactionIndex] = ers_ptr->m_reactionMultiplier;
    }
    egr->setupElemRxnVector(RxnVector, egri.m_SS_KinSpeciesKindex);

    RxnMolChange* rmcEGR = new RxnMolChange(iKA, egr);

    m_egr.push_back(egr);
    m_rmcEGR.push_back(rmcEGR);

    delete[] RxnVector;
}
//==================================================================================================================================
size_t Electrode::processExtraGlobalRxnPathways(const std::vector<EGRInput*>& EGRList)
{
    for (size_t i = 0; i < numExtraGlobalRxns; i++) {
        addExtraGlobalRxn(*(EGRList[i]));
    }
    return numExtraGlobalRxns;
}
//==================================================================================================================================
ZZCantera::ExtraGlobalRxn* Electrode::extraGlobalRxnPathway(size_t iegr)
{
    return m_egr[iegr];
}
//==================================================================================================================================
size_t Electrode::numExtraGlobalRxnPathways() const
{
    return numExtraGlobalRxns;
}
//==================================================================================================================================
Electrode::phasePop_Resid::phasePop_Resid(Electrode* ee, size_t iphaseTarget, double* const Xmf_stable,
        double deltaTsubcycle) :
    ResidEval(),
    ee_(ee),
    iphaseTarget_(iphaseTarget),
    Xmf_stable_(Xmf_stable),
    deltaTsubcycle_(deltaTsubcycle)
{
}
//==================================================================================================================================
int Electrode::phasePop_Resid::evalResidSS(const double t, const double* const y, double* const r)
{
    int retn = ee_->phasePopResid(iphaseTarget_, y, deltaTsubcycle_, r);
    return retn;
}
//==================================================================================================================================
int Electrode::phasePop_Resid::getInitialConditions(const double t0, double* const y, double* const ydot)
{
    size_t ne = nEquations();
    for (size_t k = 0; k < ne; k++) {
        y[k] = Xmf_stable_[k];
    }
    return 1;
}
//===================================================================================================================
int Electrode::phasePop_Resid::nEquations() const
{
    ThermoPhase* tp = &(ee_->thermo(iphaseTarget_));
    int nsp = tp->nSpecies();
    return nsp;
}
//===================================================================================================================
int Electrode::phasePopResid(size_t iphaseTarget, const double* const Xf_phase, double deltaTsubcycle, double* const resid)
{
    size_t k;
    std::vector<double> spMoles_tmp(m_NumTotSpecies, 0.0);

    ThermoPhase* tptarget = PhaseList_[iphaseTarget];
    size_t nspPhase = tptarget->nSpecies();
    if (nspPhase == 1) {
        resid[0] = Xf_phase[0] - 1.0;
        return 0;
    }

    size_t kstartTarget = m_PhaseSpeciesStartIndex[iphaseTarget];
    /*
     *   Set the internal objects to the correct conditions
     *    -> This will be the final conditions.
     *    -> Because we do this within an iterator loop, we are essentially trying to emulate
     *       an implicit algorithm
     */
    updateState();
    tptarget->setState_PX(pressure_, (double*) Xf_phase);

    /*
     * Loop over surface phases, filling in the phase existence fields within the
     * kinetics operator
     */
    for (size_t isk = 0; isk < numSurfaces_; isk++) {
        /*
         *  Loop over phases, figuring out which phases have zero moles.
         *  Volume phases exist if the initial or final mole numbers are greater than zero
         *  Surface phases exist if the initial or final surface areas are greater than zero.
         */
        if (ActiveKineticsSurf_[isk]) {
            ReactingSurDomain* rsd = RSD_List_[isk];
            size_t nph = rsd->nPhases();
            for (size_t jph = 0; jph < nph; jph++) {
                size_t iph = rsd->kinOrder[jph];
                if (iph == metalPhase_) {
                    continue;
                }
                double mm = phaseMoles_init_[iph];
                double mmf = phaseMoles_final_[iph];
                ThermoPhase& tp = thermo(iph);
                size_t nsp = tp.nSpecies();
                if (iph >= m_NumVolPhases) {
                    // we are in a surface phase
                    size_t isur = iph - m_NumVolPhases;
                    double sa_init = surfaceAreaRS_init_[isur];
                    double sa_final = surfaceAreaRS_final_[isur];
                    if (sa_init > 0.0 || sa_final > 0.0) {
                        rsd->setPhaseExistence(jph, true);
                    } else {
                        rsd->setPhaseExistence(jph, false);
                    }
                } else {
                    if (mm <= 0.0 && mmf <= 0.0) {
                        rsd->setPhaseExistence(jph, false);
                        if (nsp == 1 || nsp == 2) {
                            rsd->setPhaseStability(jph, true);
                        }
                    } else {
                        rsd->setPhaseExistence(jph, true);
                    }
                }
                if (iph == iphaseTarget) {
                    rsd->setPhaseExistence(jph, true);
                }
            }
        }
    }
    /*
     *  This routine basically translates between species lists for the reacting surface
     *  domain and the Electrode.
     *  Later, when we have more than one reacting surface domain in the electrode object,
     *  this will do a lot more
     */

    for (size_t isk = 0; isk < numSurfaces_; isk++) {
        // Loop over phases, figuring out which phases have zero moles.

        if (ActiveKineticsSurf_[isk]) {

            /*
             *  For each Reacting surface
             *
             *  Get the species production rates for the reacting surface
             */
            //    m_rSurDomain->getNetProductionRates(&RSSpeciesProductionRates_[0]);
            const std::vector<double>& rsSpeciesProductionRates = RSD_List_[isk]->calcNetSurfaceProductionRateDensities();

            double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(isk);
            /*
             *  loop over the phases in the reacting surface
             *  Get the net production vector
             */
            std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
            size_t nphRS = RSD_List_[isk]->nPhases();
            size_t jph, kph;
            size_t kIndexKin = 0;
            for (kph = 0; kph < nphRS; kph++) {
                jph = RSD_List_[isk]->kinOrder[kph];
                size_t istart = m_PhaseSpeciesStartIndex[jph];
                size_t nsp = m_PhaseSpeciesStartIndex[jph + 1] - istart;
                for (size_t k = 0; k < nsp; k++) {
                    spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
                    if (rsSpeciesProductionRates[kIndexKin] > 0.0) {
                        if (phaseMoles_init_[jph] <= 0.0) {
                            if (nsp > 1) {
                                size_t bornMultiSpecies = jph;
                                if (iphaseTarget != bornMultiSpecies) {
                                    throw Electrode_Error("", "two multispecies phases");
                                }
                            }
                        }
                    }
                    kIndexKin++;
                }
            }
        }
    }
    /*
     *  Find the initial surface area to use
     */
    double sa_init = surfaceAreaRS_init_[0];

    /*
     *  Find the final surface area to use
     */
    double sa_final = sa_init;

    /*
     *  Calculate the change in the moles of all of the species
     */

    for (size_t k = 0; k < m_NumTotSpecies; k++) {
        spMoles_tmp[k] = spMoles_init_[k];
        for (size_t isk = 0; isk < numSurfaces_; isk++) {
            if (ActiveKineticsSurf_[isk]) {
                sa_init = surfaceAreaRS_init_[isk];
                sa_final = surfaceAreaRS_final_[isk];
                double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(isk);
                spMoles_tmp[k] += 0.5 * deltaTsubcycle * (sa_init + sa_final) * spNetProdPerArea[k];
            }
        }
    }

    /*
     * Calculate the total moles
     */
    double ptotal = 0.0;
    for (size_t kp = 0; kp < nspPhase; kp++) {
        k = kp + kstartTarget;
        ptotal += spMoles_tmp[k];
    }
    if (ptotal <= 0.0) {
        return -1;
    }
    size_t kmax = 0;
    double xmax = Xf_phase[0];
    double xsum = 0.0;
    for (size_t kp = 0; kp < nspPhase; kp++) {
        k = kp + kstartTarget;
        resid[kp] = Xf_phase[kp] - spMoles_tmp[k] / ptotal;
        if (Xf_phase[kp] > xmax) {
            kmax = kp;
            xmax = Xf_phase[kp];
        }
        xsum += Xf_phase[kp];
    }
    resid[kmax] = 1.0 - xsum;

    return 0;
}
//===================================================================================================================
int Electrode::phasePop(size_t iphaseTarget, double* const Xmf_stable, double deltaTsubcycle)
{
    size_t k;
    std::vector<double> Xf_phase(m_NumTotSpecies, 0.0);
    int retn = 1;
    ThermoPhase* tptarget = PhaseList_[iphaseTarget];
    size_t nspPhase = tptarget->nSpecies();
    if (nspPhase == 1) {
        Xmf_stable[0] = 1.0;
        return 0;
    }
    for (k = 0; k < nspPhase; k++) {
        Xf_phase[k] = Xmf_stable[k];
    }
    tptarget->setMoleFractions(DATA_PTR(Xf_phase));

    size_t kstartTarget = m_PhaseSpeciesStartIndex[iphaseTarget];

    /*
     * Check starting conditions
     */
    for (size_t kp = 0; kp < nspPhase; kp++) {
        k = kp + kstartTarget;
        if (spMoles_init_[k] != 0.0) {
            throw Electrode_Error(" Electrode::phasePop", "spMoles_init_[k] != 0.0");
        }
        if (spMoles_final_[k] != 0.0) {
            throw Electrode_Error(" Electrode::phasePop", "spMoles_final_[k] != 0.0");
        }
    }

    /*
     * Load with an initial guess
     */
    Electrode::phasePop_Resid* pSolve_Res = new phasePop_Resid(this, iphaseTarget, Xmf_stable, deltaTsubcycle);

    ZZCantera::solveProb* pSolve = new solveProb(pSolve_Res);
    pSolve->m_ioflag = 10;
    pSolve->setAtolConst(1.0E-11);
    retn = pSolve->solve(solveProb::RESIDUAL, 1.0E-6, 1.0E-3);
    if (retn == 0) {
        pSolve->reportState(Xmf_stable);
        vector_fp resid(nspPhase);
        phasePopResid(iphaseTarget, Xmf_stable, deltaTsubcycle, DATA_PTR(resid));
    }

    delete pSolve;
    delete pSolve_Res;
    return retn;
}
//====================================================================================================================
int Electrode::numSubcycles() const
{
    return numIntegrationSubCycles_final_final_;
}
//====================================================================================================================
double Electrode::reportStateVariableIntegrationError(int& numSV, double* const errorVector) const
{
    numSV = 0;
    throw Electrode_Error("Electrode::reportStateVariableIntegrationError()", "Base Class Called");
    return 0.0;
}
//====================================================================================================================
double Electrode::reportTimeLimit(int allowedSubSteps, double allowedErrorStateVariables, double allowedSourceTermError)
{
    throw Electrode_Error("Electrode::reportTimeLimit()", "Base Class Called");
    return 0.0;
}
//==================================================================================================================================
bool Electrode::compareLocalInterval(const Electrode* const eGuest, int nDigits)
{
    //bool doComplicated = false;
    double atol = molarAtol_;
    double a1 = 8.0, a2 = 8.0;
    
    if (m_NumTotSpecies != eGuest->m_NumTotSpecies) {
	//doComplicated = true;
    }

    int res = compareOtherPL(eGuest);
    if (res != 0) {
	printf("Electrode::compareLocalInterval:  Underlying PhaseLists of two electrode objects are different, res = %d\n", res);
        return false;
    }

    bool goodComp = esmodel::doubleEqual(a1, a2, atol, nDigits);

    printf("Electrode::compareLocalInterval NOT CODED YET\n");
    return goodComp;
}
//==================================================================================================================================
Electrode::integrate_ResidJacEval::integrate_ResidJacEval(Electrode* ee) :
    ResidJacEval(),
    ee_(ee)
{
}
//  -----------------------------------------------------------------------------------------------------------------
int Electrode::integrate_ResidJacEval::evalResidNJ(double t, const double deltaT, const double* const y,
        const double* const ydot, double* const resid, const ResidEval_Type evalType, int id_x,
        double delta_x)
{
    int retn = ee_->integrateResid(t, deltaT, y, ydot, resid, evalType, id_x, delta_x);
    return retn;
}
//  -----------------------------------------------------------------------------------------------------------------
int Electrode::integrate_ResidJacEval::getInitialConditionsWithDot(const double t0, double* const y, double* const ydot)
{
    int ne = nEquations();
    for (int k = 0; k < ne; k++) {
        y[k] = 0.0;
    }
    return 1;
}
//  -----------------------------------------------------------------------------------------------------------------
int Electrode::integrate_ResidJacEval::nEquations() const
{
    return neq_;
}
//==================================================================================================================================
//  Residual calculation for the solution of the Nonlinear integration problem
/*
 *
 * @param t             Time                    (input)
 * @param delta_t       The current value of the time step (input)
 * @param y             Solution vector (input, do not modify)
 * @param ydot          Rate of change of solution vector. (input, do not modify)
 * @param resid         Value of the residual that is computed (output)
 * @param evalType      Type of the residual being computed (defaults to Base_ResidEval)
 * @param id_x          Index of the variable that is being numerically differenced to find
 *                      the jacobian (defaults to -1, which indicates that no variable is being
 *                      differenced or that the residual doesn't take this issue into account)
 * @param delta_x       Value of the delta used in the numerical differencing
 */
int Electrode::integrateResid(const double tfinal, const double deltaTsubcycle, const double* const y,
        const double* const ydot, double* const resid, const ResidEval_Type evalType, const int id_x,
        const double delta_x)
{
    throw Electrode_Error(" Electrode::integrateResid()", "Base class called");
    return 0;
}
//====================================================================================================================
// The internal state of the electrode must be kept for the initial and final times of an integration step.
/*
 *  This function advances the initial state to the final state that was calculated  in the last integration step. 
 *
 *  If the initial time is input, then the code doesn't advance
 *  or change anything. It reverts the code to the t_init_init condition wiping clean variables that are filled
 *  with the integration.
 *
 *  We also assume that the final state is equal to the final_final state. If this is not the case then this is an error.
 */
void Electrode::resetStartingCondition(double Tinitial, bool doResetAlways)
{
    bool resetToInitInit = false;
    /*
     * If this routine is called with Tinitial = t_init_init_, then we should return without doing anything
     * We have already advanced the time step to the new time.
     */
    double tbase = std::max(t_init_init_, 1.0E-50);
    if (fabs(Tinitial - t_init_init_) < (1.0E-13 * tbase) && !doResetAlways) {
        resetToInitInit = true;
    }

    if (!resetToInitInit) {
       if (Tinitial > t_init_init_) {
           if (writeRestartFileOnSuccessfulStep_) {
                if (writeRestartFileOnSuccessfulStep_ > 1) {
                    std::string extraN = "_" + int2str(globalTimeStepNumber_);
                    writeRestartFile(globalTimeStepNumber_, 1, extraN);
                } else {
                    writeRestartFile(globalTimeStepNumber_);
                }
           }
           globalTimeStepNumber_++;
       }
    }
    /*
     *  The final_final time must be equal to the new Tinitial time if we are not resetting to the initial condition
     */ 
    tbase = std::max(Tinitial, tbase);
    tbase = std::max(tbase, t_final_final_);
    if (!resetToInitInit) {
	if (fabs(Tinitial - t_final_final_) > (1.0E-13 * tbase)) {
	    throw Electrode_Error("Electrode::resetStartingCondition()",
			       "Tinitial " + fp2str(Tinitial) + " is not compatible with t_final_final_ " + fp2str(t_final_final_));
	}
    }
    /*
     *  Make sure that tfinal and tfinal_final_ are the same. This is a comfort condition
     */
    if (fabs(tfinal_ - t_final_final_) > (1.0E-13 * tbase)) {
        throw Electrode_Error("Electrode::resetStartingCondition()",
                           "tfinal_ " + fp2str(tfinal_) + " is not equal to t_final_final_ " + fp2str(t_final_final_) + 
                           " This condition is needed for some transfers.");
    }
    /*
     * Set the new time to the new value, Tinitial, which is also equal to tfinal_ and t_final_final_
     */
    if (resetToInitInit) {
	tinit_ = t_init_init_;
        tfinal_ = t_init_init_;
        t_final_final_ = t_init_init_;
    } else {
	t_init_init_ = Tinitial;
	tinit_ = Tinitial;
        tfinal_ = Tinitial;
        t_final_final_ = Tinitial;
    }
    /*
     *  Here is where we store the electrons discharged
     */
    if (!resetToInitInit) {
	if (pendingIntegratedStep_ == 1) {
	    electronKmolDischargedToDate_ += spMoleIntegratedSourceTerm_[kElectron_];
	}
    }
    pendingIntegratedStep_ = 0;
    /*
     *  Below is close to a  redo of Electrode::setInitInitStateFromFinalFinal()
     *  Not sure if I should combine the two treatments.
     */
    //
    // Major change: do a full state change function here eventually
    //
    if (resetToInitInit) {
        //setInitStateFromInitInit();
    } else {
	Electrode::setInitStateFromFinal(true);
    }

    zeroD(m_NumTotSpecies, spMoleIntegratedSourceTerm_.data());
    zeroD(m_NumTotSpecies, spMoleIntegratedSourceTermLast_.data());
    integratedThermalEnergySourceTerm_ = 0.0;
    integratedThermalEnergySourceTermLast_ = 0.0;
    integratedThermalEnergySourceTerm_overpotential_ = 0.0;
    integratedThermalEnergySourceTerm_overpotential_Last_ = 0.0;
    integratedThermalEnergySourceTerm_reversibleEntropy_ = 0.0;
    integratedThermalEnergySourceTerm_reversibleEntropy_Last_ = 0.0;

    /*
     *  Change the initial subcycle time delta here. Note, we should not change it during the integration steps
     *  because we want jacobian calculations to mainly use the same time step history, so that the answers are
     *  comparible irrespective of the time step truncation error.
     */;
    if (!resetToInitInit) {
	if (deltaTsubcycle_init_next_ < 1.0E299) {
	    deltaTsubcycle_init_init_ = deltaTsubcycle_init_next_;
	}
	deltaTsubcycle_init_next_ = 1.0E300;
    }
}
//==================================================================================================================================
int Electrode::globalTimeStepNumber() const
{
    return globalTimeStepNumber_;
}
//==================================================================================================================================
void Electrode::setGlobalTimeStepNumber(int timestep)
{
    globalTimeStepNumber_ = timestep;
}
//==================================================================================================================================
double Electrode::timeInit() const
{
    return tinit_;
}
//==================================================================================================================================
double Electrode::timeFinal() const
{
    return tfinal_;
}
//==================================================================================================================================
double Electrode::timeInitInit() const
{
    return t_init_init_;
}
//==================================================================================================================================
double Electrode::timeFinalFinal() const
{
    return t_final_final_;
}
//==================================================================================================================================
// Set the internal initial intermediate and initial global state from the internal final state
/*
 *  (non-virtual function)  -> function should onionize in-first.
 *
 *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
 *  routine as well.
 *
 * @param setInitInit   Boolean indicating whether you should set the init_init state as well
 */
void Electrode::setInitStateFromFinal_Oin(bool setInitInit)
{
    // reset surface quantities
    for (size_t i = 0; i < (size_t) numSurfaces_; i++) {
        surfaceAreaRS_init_[i] = surfaceAreaRS_final_[i];
    }
    if (setInitInit) {
        for (size_t i = 0; i < (size_t) numSurfaces_; i++) {
            surfaceAreaRS_init_init_[i] = surfaceAreaRS_final_[i];
        }
    }

    // Reset total species quantities
    for (size_t k = 0; k < m_NumTotSpecies; k++) {
        spMoles_init_[k] = spMoles_final_[k];
        spMf_init_[k] = spMf_final_[k];
        spMoles_final_final_[k] = spMoles_final_[k];
        enthalpyMolar_init_[k] = enthalpyMolar_final_[k];
	entropyMolar_init_[k] = entropyMolar_final_[k];
	chempotMolar_init_[k] = chempotMolar_final_[k];
        if (setInitInit) {
            spMf_init_init_[k] = spMf_final_[k];
            spMoles_init_init_[k] = spMoles_final_[k];
            enthalpyMolar_init_init_[k] = enthalpyMolar_final_[k];
	    entropyMolar_init_init_[k] = entropyMolar_final_[k];
	    chempotMolar_init_init_[k] = chempotMolar_final_[k];
        }
    }

    // Reset the total phase moles quantities
    for (size_t i = 0; i < m_NumTotPhases; i++) {
        phaseMoles_init_[i] = phaseMoles_final_[i];
        if (setInitInit) {
            phaseMoles_init_init_[i] = phaseMoles_final_[i];
        }
    }

    // Reset the particle size
    if (setInitInit) {
#ifdef DEBUG_MODE
        if (tfinal_ != t_init_init_) {
        if (pendingIntegratedStep_ ) {
            printf("ElectrodesetInitStateFromFinal(): t_final_ = %20.13E, t_init_init_ = %20.13E\n",
                        tfinal_ , t_init_init_);
        }
        }
        if (pendingIntegratedStep_ ) {
        if (tfinal_ != t_init_init_) {
            printf("ElectrodesetInitStateFromFinal(): t_final_ = %20.13E, t_init_init_ = %20.13E\n",
                        tfinal_ , t_init_init_);
        }
        }
#endif
        if (fabs(tfinal_ - t_init_init_) > 1.0E-40) {
            if (pendingIntegratedStep_) {
                throw Electrode_Error("Electrode::setInitStateFromFinal(true)",
                                   "Function called to overwrite init_init during a pending step");
            }
        }
        Radius_exterior_init_init_ = Radius_exterior_final_;
    }
    Radius_exterior_init_ = Radius_exterior_final_;
    Radius_exterior_final_final_ = Radius_exterior_final_;

    if (xmlStateData_final_) {
        SAFE_DELETE(xmlStateData_init_);
        xmlStateData_init_ = new XML_Node(*xmlStateData_final_);
        SAFE_DELETE(xmlStateData_final_final_);
        xmlStateData_final_final_ = new XML_Node(*xmlStateData_final_);
        if (setInitInit) {
            SAFE_DELETE(xmlStateData_init_init_);
            xmlStateData_init_init_ = new XML_Node(*xmlStateData_final_);
        }
    }

    if (xmlExternalData_final_) {
        SAFE_DELETE(xmlExternalData_init_);
        xmlStateData_init_ = new XML_Node(*xmlExternalData_final_);
        SAFE_DELETE(xmlExternalData_final_final_);
        xmlExternalData_final_final_ = new XML_Node(*xmlExternalData_final_);
        if (setInitInit) {
            SAFE_DELETE(xmlExternalData_init_init_);
            xmlExternalData_init_init_ = new XML_Node(*xmlExternalData_final_);
        }
    }
    // Set the times to the final time
    tinit_ = tfinal_;
    //t_final_final_ = tfinal_;
    if (setInitInit) {
        t_init_init_ = tfinal_;
    }
}
//====================================================================================================================
// Set the internal initial intermediate and initial global state from the internal final state
/*
 *  (virtual function)
 *
 *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
 *  routine as well.
 *
 * @param setInitInit   Boolean indicating whether you should set the init_init state as well
 */
void Electrode::setInitStateFromFinal(bool setInitInit)
{
    setInitStateFromFinal_Oin(setInitInit);
}
//====================================================================================================================
//! Set the internal initial intermediate and initial global state from the internal final_final state
/*!
 *  (virtual function)
 *
 *  Set the intial  and init_int state and the final_final from the final state. We also can set the init_init state from this
 *  routine as well.
 *
 */
void Electrode::setInitInitStateFromFinalFinal()
{
    // reset surface quantities
    for (size_t i = 0; i < (size_t) numSurfaces_; i++) {
        surfaceAreaRS_init_[i] = surfaceAreaRS_final_final_[i];
        surfaceAreaRS_init_init_[i] = surfaceAreaRS_final_final_[i];
    }

    // Reset total species quantities
    for (size_t k = 0; k < m_NumTotSpecies; k++) {
        spMoles_init_[k] = spMoles_final_final_[k];
        spMf_init_[k] = spMf_final_final_[k];
        spMf_init_init_[k] = spMf_final_final_[k];
        spMoles_init_init_[k] = spMoles_final_final_[k];
        enthalpyMolar_init_[k] = enthalpyMolar_final_final_[k];
        enthalpyMolar_init_init_[k] = enthalpyMolar_final_final_[k];
	entropyMolar_init_[k] = entropyMolar_final_final_[k];
	entropyMolar_init_init_[k] = entropyMolar_final_final_[k];
	chempotMolar_init_[k] = chempotMolar_final_final_[k];
	chempotMolar_init_init_[k] = chempotMolar_final_final_[k];
    }

    // Reset the total phase moles quantities
    for (size_t i = 0; i < m_NumTotPhases; i++) {
        phaseMoles_init_[i] = phaseMoles_final_final_[i];
        phaseMoles_init_init_[i] = phaseMoles_final_final_[i];
    }

    if (pendingIntegratedStep_) {
        throw Electrode_Error("Electrode::setInitInitStateFromFinalFinal(true)",
                           "Function called to overwrite init_init during a pending step");
    }

    Radius_exterior_init_init_ = Radius_exterior_final_final_;
    Radius_exterior_init_ = Radius_exterior_final_final_;

    if (xmlStateData_final_final_) {
        SAFE_DELETE(xmlStateData_init_);
        xmlStateData_init_ = new XML_Node(*xmlStateData_final_final_);
        SAFE_DELETE(xmlStateData_init_init_);
        xmlStateData_init_init_ = new XML_Node(*xmlStateData_final_final_);
    }

    if (xmlExternalData_final_final_) {
        SAFE_DELETE(xmlExternalData_init_);
        xmlStateData_init_ = new XML_Node(*xmlExternalData_final_final_);
        SAFE_DELETE(xmlExternalData_init_init_);
        xmlExternalData_init_init_ = new XML_Node(*xmlExternalData_final_final_);
    }

    tinit_ = t_final_final_;
    //t_final_final_ = tfinal_;

    t_init_init_ = t_final_final_;

}
//====================================================================================================================
// Revert the object's conditions to the initial conditions
/*!
 *  If there was a solution in t_final, this is wiped out and replaced with the solution at t_init.
 *  We get rid of the pending flag here as well.
 *
 * @param revertToInitInit revert to the t_init_init solution
 *                         Defaults to true.
 */
void Electrode::revertToInitialTime(bool revertToInitInit)
{
    if (revertToInitInit) {
       setInitStateFromInitInit(true);
       // need this to reset t_final_final_
       setFinalFinalStateFromFinal();
       pendingIntegratedStep_ = 0;
    } else {
       setFinalStateFromInit();
       setFinalFinalStateFromFinal();
    }
}
//====================================================================================================================
// Set the internal final intermediate from the internal init state
/*
 *  (nonvirtual function onionized from outer to inner)
 *
 *  Set the final state from the init state.
 *
 */
void Electrode::setFinalStateFromInit_Oin()
{
    // reset surface quantities
    for (size_t i = 0; i < (size_t) numSurfaces_; i++) {
        surfaceAreaRS_final_[i] = surfaceAreaRS_init_[i];
    }
    // Reset total species quantities
    for (size_t i = 0; i < m_NumTotSpecies; i++) {
        spMoles_final_[i] = spMoles_init_[i];
        spMf_final_[i] = spMf_init_[i];
        enthalpyMolar_final_[i] = enthalpyMolar_init_[i]; 
	entropyMolar_final_[i] = entropyMolar_init_[i]; 
	chempotMolar_final_[i] = chempotMolar_init_[i]; 
    }
    // Reset the total phase moles quantities
    for (size_t i = 0; i < m_NumTotPhases; i++) {
        phaseMoles_final_[i] = phaseMoles_init_[i];
    }
    // Reset the particle size
    Radius_exterior_final_ = Radius_exterior_init_;

    if (xmlStateData_init_) {
        SAFE_DELETE(xmlStateData_final_);
        xmlStateData_final_ = new XML_Node(*xmlStateData_init_);
    }
    if (xmlExternalData_init_) {
        SAFE_DELETE(xmlExternalData_final_);
        xmlExternalData_final_ = new XML_Node(*xmlExternalData_init_);
    }
    tfinal_ = tinit_;
}
//====================================================================================================================
void Electrode::setFinalStateFromInit()
{
    setFinalStateFromInit_Oin();
}
//====================================================================================================================
//    Set the internal initial intermediatefrom the internal initial global state
/*
 *  Set the intial state from the init init state. We also can set the final state from this
 *  routine as well.
 *
 *  The final_final is not touched.
 *
 * @param setFinal   Boolean indicating whether you should set the final as well
 */
void Electrode::setInitStateFromInitInit(bool setFinal)
{

    // reset surface quantities
    for (size_t i = 0; i < numSurfaces_; i++) {
        surfaceAreaRS_init_[i] = surfaceAreaRS_init_init_[i];
    }
    // Reset total species quantities
    if (setFinal) {
        for (size_t k = 0; k < m_NumTotSpecies; k++) {
            spMoles_init_[k] = spMoles_init_init_[k];
            spMoles_final_[k] = spMoles_init_init_[k];
            spMf_init_[k] = spMf_init_init_[k];
            spMf_final_[k] = spMf_init_init_[k];
            enthalpyMolar_init_[k] = enthalpyMolar_init_init_[k]; 
            enthalpyMolar_final_[k] = enthalpyMolar_init_init_[k]; 
	    entropyMolar_init_[k] = entropyMolar_init_init_[k];
	    entropyMolar_final_[k] = entropyMolar_init_init_[k];
	    chempotMolar_init_[k] = chempotMolar_init_init_[k];
	    chempotMolar_final_[k] = chempotMolar_init_init_[k];
        }
    } else {
        for (size_t k = 0; k < m_NumTotSpecies; k++) {
            spMoles_init_[k] = spMoles_init_init_[k];
            spMf_init_[k] = spMf_init_init_[k];
            enthalpyMolar_init_[k] = enthalpyMolar_init_init_[k]; 
	    entropyMolar_init_[k] = entropyMolar_init_init_[k];
	    chempotMolar_init_[k] = chempotMolar_init_init_[k];
        }

    }
    // Reset the total phase moles quantities
    for (size_t i = 0; i < m_NumTotPhases; i++) {
        phaseMoles_init_[i] = phaseMoles_init_init_[i];
    }
    if (setFinal) {
        for (size_t i = 0; i < m_NumTotPhases; i++) {
            phaseMoles_final_[i] = phaseMoles_init_init_[i];
        }
    }
    // Reset the particle size
    Radius_exterior_init_ = Radius_exterior_init_init_;
    if (setFinal) {
        Radius_exterior_final_ = Radius_exterior_init_init_;
    }
    // Copy XML data structures
    if (xmlStateData_init_init_) {
        SAFE_DELETE(xmlStateData_init_);
        xmlStateData_init_ = new XML_Node(*xmlStateData_init_init_);
        if (setFinal) {
            SAFE_DELETE(xmlStateData_final_);
            xmlStateData_final_ = new XML_Node(*xmlStateData_init_init_);
        }
    }
    if (xmlExternalData_init_init_) {
        SAFE_DELETE(xmlExternalData_init_);
        xmlStateData_init_ = new XML_Node(*xmlExternalData_init_init_);
        if (setFinal) {
            SAFE_DELETE(xmlExternalData_final_);
            xmlExternalData_final_ = new XML_Node(*xmlExternalData_init_init_);
        }
    }
    // reset the time
    tinit_ = t_init_init_;
    if (setFinal) {
        tfinal_ = t_init_init_;
    }
}
//====================================================================================================================
// Set the internal final global state from the internal final intermediate state
/*
 *  (virtual function from Electrode)
 *
 *  Set the final_final state from the final state. This is commonly called at the end of successful base integration
 */
void Electrode::setFinalFinalStateFromFinal()
{
    Electrode::setFinalFinalStateFromFinal_Oin();
}
//====================================================================================================================
//   Set the internal final global state from the internal final intermediate state
/*
 *  (non-virtual function from Electrode)
 *
 *  Set the final_final state from the final state. This is commonly called at the end of successful base integration.
 *  This is an onionize in function
 */
void Electrode::setFinalFinalStateFromFinal_Oin()
{
    for (size_t i = 0; i < numSurfaces_; i++) {
        surfaceAreaRS_final_final_[i] = surfaceAreaRS_final_[i];
    }

    for (size_t i = 0; i < m_NumTotPhases; i++) {
        phaseMoles_final_final_[i] = phaseMoles_final_[i];
    }

    for (size_t k = 0; k < m_NumTotSpecies; k++) {
        spMoles_final_final_[k] = spMoles_final_[k];
        spMf_final_final_[k] = spMf_final_[k];
    }

    Radius_exterior_final_final_ = Radius_exterior_final_;

    if (xmlStateData_final_) {
        SAFE_DELETE(xmlStateData_final_final_);
        xmlStateData_final_final_ = new XML_Node(*xmlStateData_final_);
    }

    if (xmlExternalData_final_) {
        SAFE_DELETE(xmlExternalData_final_final_);
        xmlExternalData_final_final_ = new XML_Node(*xmlExternalData_final_);
    }
    // We may not have reached the final time
    // 
    t_final_final_ = tfinal_;
}
//====================================================================================================================
//  Calculate the change in the state of the system when integrating from Tinitial to Tfinal
//  at constant current
/*
 *
 *  @param current    Current in amps. On return it returns the current actually obtained
 *  @param deltaT     DeltaT for the integration step.
 *  @param phiMax     Maximum value of the voltage. Defaults to  100.
 *  @param phiMin     Minimum value of the voltage. Defaults to -100.
 *  @param maxIntegrationSteps  defaults to 5000
 *
 *  @return Return the voltage used to obtain the current
 */
double Electrode::integrateConstantCurrent(double& current, double& deltaT, double phiMax,
                                           double phiMin, double rtolInt, int maxIntegrationSteps)
{
    int status;
    double currentNeeded = current;
    Electrode_ECurr ec(this, deltaT);
    ec.rtolInt_ = rtolInt;
    //printf("current needed = %g\n", currentNeeded);

    /*
     *  We take our cues about the best voltages from the electrode object itself
     */
    double phiM = phaseVoltages_[metalPhase_];
    double phimax = phiMax;
    double phimin = phiMin;

    if (phiMax == 100. && phiMin == -100.) {
        phimax = phiM + 2.0;
        phimin = phiM - 2.0;
    }
    if (phiMax < phiM) {
	throw Electrode_Error("Electrode::integrateConstantCurrent()",
			   "phiMax , " + fp2str(phiMax) + ", is less than starting phi, " + fp2str(phiM));
    }
    if (phiMin > phiM) {
	throw Electrode_Error("Electrode::integrateConstantCurrent()",
			   "phiMin , " + fp2str(phiMin) + ", is greater than starting phi, " + fp2str(phiM));
    }

    RootFind rf(&ec);
    rf.setTol(1.0E-5, 1.0E-10);
    rf.setFuncIsGenerallyIncreasing(true);
    rf.setDeltaX(0.01);
    double xbest = phiM;
    int oldP = printLvl_;
    double deltaT_curr = deltaT;
    // Turn down printing for lower levels
    printLvl_ = std::max(0, printLvl_ - 3);
    int printLvlLocal = oldP;
    if (printLvlLocal > 0) {
        rf.setPrintLvl(printLvlLocal);
        ec.printLvl_ = printLvlLocal;
    } else {
        rf.setPrintLvl(0);
        ec.printLvl_ = 0;
    }
    Electrode_Integrator* eei = dynamic_cast<Electrode_Integrator*>(this);
  
    int numSteps = integrate(deltaT, rtolInt,  T_FINAL_CONST_FIS, BASE_TIMEINTEGRATION_SIR);
    if (numSteps > maxIntegrationSteps) {
	if (eei) {
	    SubIntegrationHistory& sih = eei->timeHistory();
	    if (printLvlLocal > 2) {
		sih.print(5);
	    } 
	    int nn = 3 * maxIntegrationSteps / 4 + 1;
	    TimeStepHistory& tsh = sih.TimeStepList_[nn];
	    deltaT_curr = tsh.t_final_calc_;
	} else {
	    deltaT_curr = 0.5 * deltaT *  maxIntegrationSteps /  numSteps;
	}
    } else {
	if (printLvlLocal > 1) {
	    if (eei) {
		SubIntegrationHistory& sih = eei->timeHistory(); 
		sih.print(5);
	    }
	}
    }


    /*
     *  We put the rootfinder in a loop here. Sometimes the rootfinder will fail. The algorithm is to cut
     *  the time in half and retry the problem again.
     */
    for (int iTrial = 0; iTrial < 5; iTrial++) {
        ec.set_deltaT(deltaT_curr);
        double currentObtained = currentNeeded;
        xbest = phiM;
        status = rf.solve(phimin, phimax, 100, currentObtained, &xbest);

        /*
         *  Status = 0 means success. Return with a success indicator.
         */
        if (status == 0) {
            if (printLvl_ > 1) {
                printf("Electrode::integrateConstantCurrent(): Volts (%g amps) = %g at deltaT = %g\n", currentObtained,
                       xbest, deltaT_curr);
            }
	    numSteps = ec.nIntegrationSteps();
	    if (eei) {
		SubIntegrationHistory& sih = eei->timeHistory(); 
		if (printLvl_ > 2) {
		    sih.print(5);
		}
	    
		if (numSteps > maxIntegrationSteps) {
		    int nn = 3 * maxIntegrationSteps / 4 + 1;
		    TimeStepHistory& tsh = sih.TimeStepList_[nn];
		    deltaT_curr = tsh.t_final_calc_ - t_init_init_;
		    deltaT = deltaT_curr;
		    current = currentObtained;
		    continue;
		}
	    }
            deltaT = deltaT_curr;
            current = currentObtained;
            printLvl_ = oldP;
            return xbest;
        }
        /*
         *  Any other status indicator, we print out a warning message
         */
        if (printLvlLocal > 1) {
            printf("Electrode::integrateConstantCurrent(): bad status = %d Volts (%g amps) = %g at deltaT = %g\n",
                   status, currentObtained, xbest, deltaT_curr);
        }
        /*
         *  Deal with the returned voltage being at the top limit voltage
         */
        if (xbest >= phiMax) {
            /*
             *  If we have tried the rootfinder 3 times with decreasing time steps and we are at the max voltage,
             *  then give up early. We are probably at an end-of-electron condition. With this condition, we need to
             *  exit gracefully and handle the situation at a higher level.
             */
            if (iTrial >= 2) {
                if (printLvlLocal > 1) {
                    printf("                            volts at top = %g : returning with voltMax solution\n", phiMax);
                }
                current = currentObtained;
                deltaT = deltaT_curr;
                printLvl_ = oldP;
                return xbest;
            } else {
                if (printLvlLocal > 1) {
		    if (fabs(currentObtained ) < 1.0E-10) {
			printf("                            volts at top = %g probably because no electrons left, curr = %g\n", phiMax, currentObtained);
		    } else {
			printf("                            volts at top = %g\n", phiMax);
		    }
                }
		if (fabs(currentObtained ) < 1.0E-10) {
		    current = currentObtained;
		    deltaT = deltaT_curr;
		    printLvl_ = oldP;
		    return xbest;
		}
            }
        }
        /*
         *  Deal with the returned voltage being at the bottom limit voltage
         */
        if (xbest <= phiMin) {
            /*
             *  If we have tried the rootfinder 3 times with decreasing time steps and we are at the min voltage,
             *  then give up early. We are probably at an end-of-electron condition. With this condition, we need to
             *  exit gracefully and handle the situation at a higher level.
             */
            if (iTrial >= 2) {
                if (printLvlLocal > 1) {
                    printf("                            volts at bottom = %g: returning with voltMin solution\n",
                           phiMin);
                }
                current = currentObtained;
                deltaT = deltaT_curr;
                printLvl_ = oldP;
                return xbest;
            } else {
		if (printLvlLocal > 1) {
		    if (fabs(currentObtained ) < 1.0E-6) {
			printf("                            volts at bottom = %g probably because no electrons left, curr = %g\n", phiMin, currentObtained);
		    } else {
			printf("                            volts at bottom = %g, current = %g\n", phiMin,  currentObtained);
		    }
                }
		if (fabs(currentObtained ) < 1.0E-6) {
		    current = currentObtained;
		    deltaT = deltaT_curr;
		    printLvl_ = oldP;
		    return xbest;
		}
            }
        }
        /*
         *  We are here for ealy cases of max or min voltages, where we are running out of electrons. We are
         *  also here when the rootfinder fails for indeterminate reasons. For this latter case, we will try
         *  reducing the time step to see if this helps solve the system. In some cases, this has actually worked,
         *  while in other cases, it has not worked.
         */
        /*
         *  Fill in the return variables, current and deltaT, with the last quantities actually obtained from
         *  the rootfinder. Then reduce the time step by a factor of two and resolve the system.
         *  Maybe we can find the requested current at a smaller time step increment.  This often
         *  happens when we are running out of electrons in the electrode.
         */
        current = currentObtained;
        deltaT = deltaT_curr;
        deltaT_curr *= 0.5;
    }
    /*
     *  Return with failed indicator values. We return with the last voltage, xbest, the last current, current,
     *  and the last deltaT, deltaT, actually returned from the rootfinder.
     */
    if (printLvl_ > 1) {
        printf("                            Returning with failed rootfinder values\n");
    }
    printLvl_ = oldP;
    return xbest;
}
//==================================================================================================================================
/*
 *  Report the integrated source term for the electrode over an interval in time.
 *  This is the net change in the moles of species defined in the electrode over that interval of time.
 *  The number of subcycles is returned, as well.
 */
size_t Electrode::integratedSpeciesSourceTerm(double* const spMoleDelta)
{
#ifdef DEBUG_MODE
    if (tfinal_ == tinit_) {
        throw Electrode_Error(" Electrode::integratedSpeciesSourceTerm()", "tfinal == tinit");
    }
#endif
    Zuzax::copyD(m_NumTotSpecies, spMoleDelta, spMoleIntegratedSourceTerm_.data());
    return numIntegrationSubCycles_final_final_;
}
//==================================================================================================================================
// Calculate the instantaneous time derivative of the species vector as determined by all source terms
/*
 *  This is the rate of change in the moles of species defined in the electrode  at t_final.
 *  This calculation does not necessarily use an interval of time to calculate anything.
 *
 *  @param spMoleDot   The end result in terms of the rate of change in moles of species in the
 *                     electrode. phaseList format. (kmol s-1)
 */
void Electrode::speciesProductionRates(double* const spMoleDot)
{
    std::fill_n(spMoleDot, m_NumTotSpecies, 0.);
    //
    //  For non-pending we calculate the instantaneous value
    //
    if (pendingIntegratedStep_ == 1) {
	printf("WARNING: speciesProductionRate called with pendingIntegratedStep_ = 1\n");
    }
    //
    // Look over active kinetics surfaces
    //
    for (size_t isk = 0; isk < numSurfaces_; isk++) {
	if (ActiveKineticsSurf_[isk]) {
	    /*
	     *  For each Reacting surface
	     *      (  m_rSurDomain->getNetProductionRates(&RSSpeciesProductionRates_[0]);
	     *  Get the species production rates for the reacting surface
	     */
	    // TODO: Check this logic for end of region conditions and goNowhere issues
	    const std::vector<double>& rsSpeciesProductionRates = RSD_List_[isk]->calcNetSurfaceProductionRateDensities();
	    
	    /*
	     *  loop over the phases in the reacting surface
	     *  Get the net production vector
	     */
	    size_t nphRS = RSD_List_[isk]->nPhases();
	    int kIndexKin = 0;
	    for (size_t kph = 0; kph < nphRS; kph++) {
		size_t jph = RSD_List_[isk]->kinOrder[kph];
		int istart = m_PhaseSpeciesStartIndex[jph];
		int nsp = m_PhaseSpeciesStartIndex[jph + 1] - istart;
		for (int k = 0; k < nsp; k++) {
		    spMoleDot[istart + k] += rsSpeciesProductionRates[kIndexKin] * surfaceAreaRS_final_[isk];
		    kIndexKin++;
		}
	    }
	}
    }
}
//====================================================================================================================
//! Report the enthalpy source term for the electrode over an interval in time
/*!
 * FIX!
 *  Sum over phases ( enthalpy phase * (phaseMoles_final_ - phaseMoles_init_init_) )
 *  This should only be called after integrate() has finished running.
 */
double Electrode::integratedEnthalpySourceTerm()
{
    double energySource = 0.0;
    for (size_t iph = 0; iph < m_NumVolPhases; iph++) {
        ThermoPhase& tp = thermo(iph);
        //int nspPhase = tp.nSpecies();
	// HKM -> this formula is in error. Enthalpy per mole could have changed over the interval.
        energySource += tp.enthalpy_mole() * (phaseMoles_final_[iph] - phaseMoles_init_init_[iph]);
    }
    return energySource;
}
//==================================================================================================================================
// Overpotential term for the heat generation
double Electrode::thermalEnergySourceTerm_overpotential(size_t isk)
{   
    double nstoich, ocv, io, nu, beta, resist;
    double iCurr;
    double q = 0.0;
    if (RSD_List_[isk]) {
         ReactingSurDomain* rsd = RSD_List_[isk];
         double sa = surfaceAreaRS_final_[isk];
         size_t nr = rsd->nReactions();
         for (size_t irxn = 0; irxn < nr; irxn++) {
             double overpotential = overpotentialRxn(isk, (int) irxn);
      
             bool okk = rsd->getExchangeCurrentDensityFormulation(irxn, nstoich, ocv, io, nu, beta, resist);
             if (okk) {
                iCurr = rsd->calcCurrentDensity(nu, nstoich, io, beta, temperature_, resist);
             } else {
                iCurr = 0.0;
             }

	     if (nstoich != 0.0) {
		 q += sa * iCurr * overpotential;
	     }
         }
    }
    return q;
}
//===================================================================================================================================
// Reversible Enthalpy term leading to heat generation
double Electrode::thermalEnergySourceTerm_reversibleEntropy(size_t isk)
{
    double nstoich, ocv, io, nu, beta, resist;
    double iCurr;
    double q = 0.0;
    static std::vector<double> s_deltaS;
    static std::vector<double> iCurrDens;
    if (RSD_List_[isk]) {
         ReactingSurDomain* rsd = RSD_List_[isk];
         double sa = surfaceAreaRS_final_[isk];
         size_t nr = rsd->nReactions();
	 size_t ss = std::max(s_deltaS.size(), nr);
	 s_deltaS.resize(ss, 0.0);
	 iCurrDens.resize(ss, 0.0);
	 rsd->getDeltaEntropy(&(s_deltaS[0]));
         iCurr = rsd->getCurrentDensityRxn(&(iCurrDens[0]));
	 double tt = temperature_;

         for (size_t irxn = 0; irxn < nr; irxn++) {
             bool okk = rsd->getExchangeCurrentDensityFormulation(irxn, nstoich, ocv, io, nu, beta, resist);
             if (okk) {
                iCurr = rsd->calcCurrentDensity(nu, nstoich, io, beta, temperature_, resist);
             } else {
                iCurr = 0.0;
             }

	     if (nstoich != 0.0) {
		 q -= sa * iCurr * tt * s_deltaS[irxn] / Faraday;
	     } else {

	     }

         }

    }
    return q; 
}
//====================================================================================================================
// Reversible Entropy term leading to  heat generation
/*
 *
 */
double Electrode::thermalEnergySourceTerm_EnthalpyFormulation(size_t isk)
{
    double nstoich, ocv, io, nu, beta, resist;
    double iCurr;
    double q = 0.0;
    static std::vector<double> s_deltaH;
    static std::vector<double> iCurrDens;
   
    if (RSD_List_[isk]) {
         ReactingSurDomain* rsd = RSD_List_[isk];
         double sa = surfaceAreaRS_final_[isk];
         size_t nr = rsd->nReactions();
	 size_t ss = std::max(s_deltaH.size(), nr);
	 s_deltaH.resize(ss, 0.0);
	 iCurrDens.resize(ss, 0.0);
	 rsd->getDeltaEnthalpy(&(s_deltaH[0]));
         iCurr = rsd->getCurrentDensityRxn(&(iCurrDens[0]));
	 const std::vector<double>& ROP = rsd->calcNetSurfaceROP();
#ifdef DEBUG_THERMAL
         for (int k = 0; k < m_NumTotSpecies; k++) {
            printf("h %d = %20.9E\n", k, enthalpyMolar_final_[k]);
         }
#endif
         // HKM -> what happens here when not in exchange current format?
         for (size_t irxn = 0; irxn < nr; irxn++) {
             bool okk = rsd->getExchangeCurrentDensityFormulation(irxn, nstoich, ocv, io, nu, beta, resist);
             if (okk) {
                iCurr = rsd->calcCurrentDensity(nu, nstoich, io, beta, temperature_, resist);
             } else {
                iCurr = 0.0;
             }
#ifdef DEBUG_THERMAL
	     double deltaM = - iCurr * sa / Faraday;
             printf ("delta moles = %g\n", deltaM * deltaTsubcycle_);
             printf ("delta moles rate = %g\n", deltaM);
             printf(" deltaH = %20.9E\n", s_deltaH[irxn]);
             double term0 =  deltaM * enthalpyMolar_final_[0];
             printf(" term0 =  %g \n", deltaM * enthalpyMolar_final_[0]);
             printf(" term0_int =  %g \n", deltaM * enthalpyMolar_final_[0] * deltaTsubcycle_);
             double term4 = - deltaM * enthalpyMolar_final_[4];
             printf(" term4 =  %g \n", - deltaM * enthalpyMolar_final_[4]);
             printf(" term4_int =  %g \n", - deltaM * enthalpyMolar_final_[4] * deltaTsubcycle_);

             double volt0_term = - deltaM * Faraday * phaseVoltages_[0];
             printf(" volt0 term = %g \n",  - deltaM * Faraday * phaseVoltages_[0]);
             printf(" volt0 term_int = %g \n",  - deltaM * Faraday * phaseVoltages_[0] * deltaTsubcycle_);
#endif
	     if (nstoich != 0.0) {
		 q -= sa * iCurr * s_deltaH[irxn] / (Faraday * nstoich);
		 q += sa * iCurr * deltaVoltage_;
	     } else {
		 q -= sa * ROP[irxn] * s_deltaH[irxn];
	     }
#ifdef DEBUG_THERMAL
             printf("  term0 + term4 + volt0 term = %g\n", term0 + term4 + volt0_term );
#endif
         }
    }
    return q; 
}
//==============================================================================================================
double Electrode::thermalEnergySourceTerm_EnthalpyFormulation_SingleStep_Old()
{
    double q = 0.0;
    double phiPhase = 0.0;
     
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        ThermoPhase& tp = thermo(iph);
        phiPhase = phaseVoltages_[iph];

        size_t istart = m_PhaseSpeciesStartIndex[iph];
        size_t nsp = tp.nSpecies();
        if (iph == metalPhase_ || iph == solnPhase_) {
	    for (size_t ik = 0; ik < nsp; ik++) {
		size_t k = istart + ik;
		double cc = tp.charge(ik);
		double deltaNH = enthalpyMolar_final_[k] * spMoleIntegratedSourceTermLast_[k];
#ifdef DEBUG_THERMAL
		double deltaN = spMoles_final_[k] - spMoles_init_[k]; 
		printf("deltaNH_%d = %g,   deltaN = %g\n", k, -deltaNH, deltaN);
#endif
		q -= deltaNH;
		if (cc != 0.0) {
		    q -= cc * Faraday * phiPhase * spMoleIntegratedSourceTermLast_[k];
#ifdef DEBUG_THERMAL
		    printf("deltV_%d = %g\n", k, - cc * Faraday * phiPhase * spMoleIntegratedSourceTermLast_[k]);
#endif
		}
	    }
        } else {
	    for (size_t ik = 0; ik < nsp; ik++) {
		size_t k = istart + ik;  
		double cc = tp.charge(ik);
		double deltaNH = enthalpyMolar_final_[k] * spMoles_final_[k] - enthalpyMolar_init_[k] * spMoles_init_[k];
#ifdef DEBUG_THERMAL
		double deltaN = spMoles_final_[k] - spMoles_init_[k];
		printf("deltaNH_%d = %g,   deltaN = %g\n", k, -deltaNH, deltaN);
#endif
		q -= deltaNH;
		if (cc != 0.0) {
		    q -= cc * Faraday * phiPhase * spMoleIntegratedSourceTermLast_[k];
#ifdef DEBUG_THERMAL
		    printf("deltV_%d = %g\n", k, - cc * Faraday * phiPhase * spMoleIntegratedSourceTermLast_[k]);
#endif
		}
	    }
        }
    }
    return q;
}
//==============================================================================================================
/*
 *   sdot h_k delta_t
 */
double Electrode::thermalEnergySourceTerm_EnthalpyFormulation_SingleStep()
{
    double q = 0.0;     
    for (size_t iph = 0; iph < (size_t) m_NumTotPhases; iph++) {
        ThermoPhase& tp = thermo(iph);
        double phiPhase = phaseVoltages_[iph];

        size_t istart = m_PhaseSpeciesStartIndex[iph];
        size_t nsp = tp.nSpecies();
	for (size_t ik = 0; ik < nsp; ik++) {
	    size_t k = istart + ik;
	    double cc = tp.charge(ik);
#ifdef DEBUG_THERMAL
	    double deltaNH = enthalpyMolar_final_[k] * spMoleIntegratedSourceTermLast_[k];
	    double deltaN = spMoles_final_[k] - spMoles_init_[k]; 
	    printf("deltaNH_%d = %g,   deltaN = %g\n", k, -deltaNH, deltaN);
#endif
	    q -= enthalpyMolar_final_[k] * spMoleIntegratedSourceTermLast_[k];
	    if (cc != 0.0) {
		q -= cc * Faraday * phiPhase * spMoleIntegratedSourceTermLast_[k];
#ifdef DEBUG_THERMAL
		printf("deltV_%d = %g\n", k, - cc * Faraday * phiPhase * spMoleIntegratedSourceTermLast_[k]);
#endif
	    }
	}
    }
    return q;
}
//==============================================================================================================
// This is called at the end of a single step.
/*
 *    we assume here that the change in moles can be represented by the set of ROP of the surface reactions
 *    on all of the surfaces.
 */
double Electrode::thermalEnergySourceTerm_ReversibleEntropy_SingleStep()
{
    double q_alt = 0.0;
    for (size_t iph = 0; iph < (size_t) m_NumTotPhases; iph++) {
        ThermoPhase& tp = thermo(iph);
        size_t istart = m_PhaseSpeciesStartIndex[iph];
        size_t nsp = tp.nSpecies();
	for (size_t ik = 0; ik < nsp; ik++) {
	    size_t  k = istart + ik;
	    double deltaNS = temperature_ * entropyMolar_final_[k] * spMoleIntegratedSourceTermLast_[k];
#ifdef DEBUG_THERMAL
	    double deltaN = spMoleIntegratedSourceTermLast_[k]; 
	    printf("deltaNS_%d = %g,   deltaN = %g\n", k, -deltaNS, deltaN);
#endif
	    q_alt -= deltaNS;
	}
    }
    return q_alt;
}
//==============================================================================================================
// This is called at the end of a single step.
/*
 *    we assume here that the change in moles can be represented by the set of ROP of the surface reactions
 *    on all of the surfaces.
 *  
 */
double Electrode::thermalEnergySourceTerm_ReversibleEntropy_SingleStep_Old()
{
    //
    //  Ok we have calculated q one way, calc q another way so that we can have options.
    // 
    double q_alt = 0.0;
    double tt = temperature_;
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        ThermoPhase& tp = thermo(iph);
        size_t istart = m_PhaseSpeciesStartIndex[iph];
        size_t nsp = tp.nSpecies();
        if (iph == metalPhase_ || iph == solnPhase_) {
	    for (size_t ik = 0; ik < nsp; ik++) {
		size_t k = istart + ik;
		double deltaNS = tt * entropyMolar_final_[k] * spMoleIntegratedSourceTermLast_[k];
#ifdef DEBUG_THERMAL
		double deltaN = spMoles_final_[k] - spMoles_init_[k]; 
		printf("deltaNS_%d = %g,   deltaN = %g\n", k, -deltaNS, deltaN);
#endif
		q_alt -= deltaNS;
	    }
        } else {
	    for (size_t ik = 0; ik < nsp; ik++) {
		size_t k = istart + ik;
		double deltaNS = tt * (entropyMolar_final_[k] * spMoles_final_[k] - entropyMolar_init_[k] * spMoles_init_[k]);
#ifdef DEBUG_THERMAL
		double deltaN = spMoles_final_[k] - spMoles_init_[k];
		printf("deltaNS_%d = %g,   deltaN = %g\n", k, -deltaNS, deltaN);
#endif
		q_alt -= deltaNS;
	    }
        }
    }
    return q_alt;
}
//==============================================================================================================
// This is called at the end of a single step.
/*
 *    we assume here that the change in moles can be represented by the set of ROP of the surface reactions
 *    on all of the surfaces.
 */
double Electrode::thermalEnergySourceTerm_Overpotential_SingleStep_Old()
{
    double q_alt = 0.0;
    double phiPhase = 0.0;    
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        ThermoPhase& tp = thermo(iph);
        phiPhase = phaseVoltages_[iph];

        size_t istart = m_PhaseSpeciesStartIndex[iph];
        size_t nsp = tp.nSpecies();
        if (iph == metalPhase_ || iph == solnPhase_) {
	    for (size_t ik = 0; ik < nsp; ik++) {
		size_t k = istart + ik;
		double cc = tp.charge(ik);
		double deltaNG = chempotMolar_final_[k] * spMoleIntegratedSourceTermLast_[k];
#ifdef DEBUG_THERMAL
		double deltaN = spMoles_final_[k] - spMoles_init_[k]; 
		printf("deltaNG_%d = %g,   deltaN = %g\n", k, -deltaNG, deltaN);
#endif
		q_alt -= deltaNG;
		if (cc != 0.0) {
		    q_alt -= cc * Faraday * phiPhase * spMoleIntegratedSourceTermLast_[k];
#ifdef DEBUG_THERMAL
		    printf("deltV_%d = %g\n", k, - cc * Faraday * phiPhase * spMoleIntegratedSourceTermLast_[k]);
#endif
		}
	    }
        } else {
	    for (size_t ik = 0; ik < nsp; ik++) {
		size_t k = istart + ik;  
		double cc = tp.charge(ik);
		double deltaNG = chempotMolar_final_[k] * spMoles_final_[k] - chempotMolar_init_[k] * spMoles_init_[k];
#ifdef DEBUG_THERMAL
		double deltaN = spMoles_final_[k] - spMoles_init_[k];
		printf("deltaNG_%d = %g,   deltaN = %g\n", k, -deltaNG, deltaN);
#endif
		q_alt -= deltaNG;
		if (cc != 0.0) {
		    q_alt -= cc * Faraday * phiPhase * spMoleIntegratedSourceTermLast_[k];
#ifdef DEBUG_THERMAL
		    printf("deltV_%d = %g\n", k, - cc * Faraday * phiPhase * spMoleIntegratedSourceTermLast_[k]);
#endif
		}
	    }
        }
    }
    return q_alt;
}
//==================================================================================================================================
// This is called at the end of a single step.
//
double Electrode::thermalEnergySourceTerm_Overpotential_SingleStep()
{
    double q_alt = 0.0;
    for (size_t iph = 0; iph < (size_t) m_NumTotPhases; iph++) {
        ThermoPhase& tp = thermo(iph);
        double phiPhase = phaseVoltages_[iph];
        size_t istart = m_PhaseSpeciesStartIndex[iph];
        size_t nsp = tp.nSpecies();
	for (size_t ik = 0; ik < nsp; ik++) {
	    size_t k = istart + ik;
	    double cc = tp.charge(ik);
	    double deltaNG = chempotMolar_final_[k] * spMoleIntegratedSourceTermLast_[k];
#ifdef DEBUG_THERMAL
	    double deltaN = spMoles_final_[k] - spMoles_init_[k]; 
	    printf("deltaNG_%d = %g,   deltaN = %g\n", k, -deltaNG, deltaN);
#endif
	    q_alt -= deltaNG;
	    if (cc != 0.0) {
		q_alt -= cc * Faraday * phiPhase * spMoleIntegratedSourceTermLast_[k];
#ifdef DEBUG_THERMAL
		printf("deltV_%d = %g\n", k, - cc * Faraday * phiPhase * spMoleIntegratedSourceTermLast_[k]);
#endif
	    }
	}
    }
    return q_alt;
}
//==================================================================================================================================
double Electrode::getIntegratedSourceTerm(SOURCES sourceType)
{
  double result = 0.0;

  size_t species_index = 0;
  size_t solnSpeciesStart = m_PhaseSpeciesStartIndex[solnPhase_];
  size_t nsp = thermo(solnPhase_).nSpecies();
  if (sourceType >= VOLUME_SOURCE) {

  }
  else if (sourceType > SPECIES_SOURCE) {
    species_index = sourceType - SPECIES_SOURCE;
    AssertThrow(species_index < nsp, "Electrode::getSourceTerm");
    sourceType = SPECIES_SOURCE;
  }

  switch (sourceType) {
  case CURRENT_SOURCE:
    result = spMoleIntegratedSourceTerm_[kElectron_];
    break;
  case ELECTROLYTE_PHASE_SOURCE:
    for (size_t ik = 0; ik < nsp; ik++) {
      size_t k = solnSpeciesStart + ik;
      result += spMoleIntegratedSourceTerm_[k];
    }
    break;
  case ENTHALPY_SOURCE:
    result = integratedEnthalpySourceTerm();
    break;
  case SPECIES_SOURCE:
    result = spMoleIntegratedSourceTerm_[solnSpeciesStart + species_index];
    break;
  case VOLUME_SOURCE:
    // For the volume source, we report the extrinsic volume of the solid portion of the electrode as it exists right now
    //  units = m^3
    result = SolidVol();
    break;
  case MAX_SOURCE:
  default:
    result = 0.0;
    break;
  }
  return result;
}
//====================================================================================================================
// Set the phase existence flag in the electrode kinetics object so that kinetics
// are calculated correctly
/*
 *    Flags are set in the kinetics object to tell the kinetics object which phases
 *    have zero moles.  The zero mole indicator is taken from phaseMoles_final_[]. Therefore,
 *    the final state is queried.
 *    There is a special case. Phases that have a 1 in the justBornPhase_[] vector are allowed to
 *    be set to exist even if their phase moles are zero.
 *
 * @param assumeStableSingleSpeciesPhases Assume that single phases are stable. This
 *                         allows their production rates to be calculated
 */
void Electrode::setPhaseExistenceForReactingSurfaces(bool assumeStableSingleSpeciesPhases)
{
    for (size_t isk = 0; isk < numSurfaces_; isk++) {
        /*
         *  Loop over phases, figuring out which phases have zero moles.
         *  Volume phases exist if the initial or final mole numbers are greater than zero
         *  Surface phases exist if the initial or final surface areas are greater than zero.
         */
        if (RSD_List_[isk]) {
            ReactingSurDomain* rsd = RSD_List_[isk];
            size_t nph = rsd->nPhases();
            for (size_t jph = 0; jph < nph; jph++) {
                size_t iph = rsd->kinOrder[jph];
                if (iph == metalPhase_) {
                    continue;
                }
                double mm = phaseMoles_init_[iph];
                double mmf = phaseMoles_final_[iph];
                ThermoPhase& tp = thermo(iph);
                size_t nsp = tp.nSpecies();
                if (iph >= m_NumVolPhases) {
                    // we are in a surface phase
                    size_t isur = iph - m_NumVolPhases;
                    double sa_init = surfaceAreaRS_init_[isur];
                    double sa_final = surfaceAreaRS_final_[isur];
                    if (sa_init > 0.0 || sa_final > 0.0) {
                        rsd->setPhaseExistence(jph, true);
                    } else {
                        rsd->setPhaseExistence(jph, false);
                    }
                } else {
                    if (mm <= 0.0 && mmf <= 0.0) {
                        rsd->setPhaseExistence(jph, false);
                        if (nsp == 1 && assumeStableSingleSpeciesPhases) {
                            rsd->setPhaseStability(jph, true);
                        }
                    } else {
                        rsd->setPhaseExistence(jph, true);
                    }
                }

                if (justBornPhase_[iph]) {
                    rsd->setPhaseExistence(jph, true);
                }
            }
        }
    }
}
//====================================================================================================================
void Electrode::printElectrodePhaseList(int pSrc, bool subTimeStep)
{
    printf("\n");
    printf("         PName                   MoleNum      molarVol     Volume       FractVol     Voltage   \n");
    printf("     ============================================================================================\n");
    double egv = TotalVol();
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        std::string pname = PhaseNames_[iph];
        ThermoPhase& tp = thermo(iph);
        double mv = tp.molarVolume();
        if (iph >= m_NumVolPhases) {
            mv = 0.0;
        }
        double pv = mv * phaseMoles_final_[iph];
        printf("     ");
        ca_ab::pr_sf_lj(pname, 24, 1);
        printf(" %12.3E", phaseMoles_final_[iph]);
        printf(" %12.3E", mv);
        printf(" %12.3E", pv);
        printf(" %12.3E", pv / egv);
        if (iph == metalPhase_ || iph == solnPhase_) {
            printf(" %12.5E", tp.electricPotential());
        } else {
            printf("             ");
        }
        printf("\n");
    }
    printf("     ============================================================================================\n");
}
//===================================================================================================================================
// Print conditions of the electrode for the current integration step to stdout
/*
 *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
 *                       to the final_final time.
 *                       The default is to print out the source terms
 *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
 *                       time step. The default is to print out the global values
 */
void Electrode::printElectrode(int pSrc, bool subTimeStep)
{
    size_t iph;
    double egv = TotalVol();
    double tsm = SolidTotalMoles();
    printf("   ==============================================================================================\n");
    if (subTimeStep) {
        printf("      Electrode at intermediate-step time final = %12.5E\n", tfinal_);
        printf("                   intermediate-step time init  = %12.5E                deltaT = %g\n", 
	       tinit_, deltaTsubcycle_);
        printf("                ChemModel Type = %3d , DomainNumber = %2d , CellNumber = %2d , SubIntegrationCounter = %d\n",
               electrodeChemistryModelType_, electrodeDomainNumber_, electrodeCellNumber_, counterNumberSubIntegrations_);
    } else {
        printf("      Electrode at time final_final = %12.5E\n", t_final_final_);
        printf("                   time init_init   = %12.5E                      deltaTglobal = %g\n",
	       t_init_init_, t_final_final_ - t_init_init_);
        printf("                ChemModel Type = %3d , DomainNumber = %2d , CellNumber = %2d , IntegrationCounter = %d\n",
               electrodeChemistryModelType_, electrodeDomainNumber_, electrodeCellNumber_, counterNumberIntegrations_);
	printf("                numIntegratioSubCycles = %d, SubIntegrationCounter = %d\n",
	       static_cast<int>(numIntegrationSubCycles_final_final_), counterNumberSubIntegrations_);
    }
    printf("   ==============================================================================================\n");
    printf("          Voltage (phiMetal - phiElectrolyte) = %12.5E volts\n", deltaVoltage_);
    if (subTimeStep) {
        double curr = integratedLocalCurrent();
        printf("          Current = %12.5E amps\n", curr);
    } else {
        double curr = integratedCurrent();
        printf("          Current = %12.5E amps\n", curr);
    }
    printf("          Number of external surfaces = %d\n", (int) numExternalInterfacialSurfaces_);
    printf("          Solid Volume = %11.4E m**3\n", ElectrodeSolidVolume_);
    printf("          Total Volume = %11.4E m**3\n", egv);
    if (egv > 0.0) {
        printf("          Porosity     = %10.3E\n", (egv - ElectrodeSolidVolume_) / egv);
    } else {
        printf("          Porosity     =       NA\n");
    }
    printf("          Temperature = %g K\n", temperature_);
    printf("          Pressure = %g Pa\n", pressure_);
    printf("          Total Solid Moles = %11.4E kmol\n", tsm);
    if (tsm > 0.0) {
        printf("          Molar Volume of Solid = %11.4E cm3 gmol-1\n", ElectrodeSolidVolume_ / tsm * 1.0E3);
    } else {
    }
    printf("          Particle Number to Follow = %11.4E\n", particleNumberToFollow_);

    printf("          followElectrolyteMoles = %d\n", followElectrolyteMoles_);
    printf("          ElectrolytePseudoMoles = %g\n", electrolytePseudoMoles_);
    printf("\n");

    printElectrodeCapacityInfo(pSrc, subTimeStep);
    printf("\n");
    printElectrodePhaseList(pSrc, subTimeStep);

    size_t m = m_NumTotPhases;
    if (numSurfaces_ > m_NumSurPhases) {
        m = numSurfaces_ + m_NumVolPhases;
    }
    for (iph = 0; iph < m; iph++) {
        printElectrodePhase(iph, pSrc, subTimeStep);
    }
    printf("   ==============================================================================================\n");
}
//===================================================================================================================
void Electrode::printElectrodeCapacityInfo(int pSrc, bool subTimeStep)
{
    double capacd = capacityDischarged();
    printf("          Capacity Discharged Since Start = %12.6g coulombs = %12.6g Ah\n", capacd, capacd / 3600.);
    double dod = depthOfDischarge();
    double rdod = depthOfDischargeFraction();
    printf("          Depth of Discharge (Current)    = %12.6g coulombs    Relative Depth of Discharge = %-12.6g\n",
           dod, rdod);
    double capLeft = capacityLeft();
    double capZero = capacity();
    printf("          Capacity Left                   = %12.6g coulombs    Capacity at Zero DOD = %-12.6g coulombs\n",
           capLeft, capZero);
}
//===================================================================================================================
void Electrode::printElectrodePhase(size_t iph, int pSrc, bool subTimeStep)
{
    size_t isurf = npos;
    double* netROP = new double[m_NumTotSpecies];
    ThermoPhase& tp = thermo(iph);
    std::string pname = tp.name();
    size_t istart = m_PhaseSpeciesStartIndex[iph];
    size_t nsp = tp.nSpecies();
    if (printLvl_ <= 1) {
        return;
    }
    printf("     ============================================================================================\n");
    printf("          PHASE %d %s \n", static_cast<int>(iph), pname.c_str());
    printf("                Total Moles  = %11.5E kmol\n", phaseMoles_final_[iph]);
    double mv = tp.molarVolume();
    double ma = 0.0;
    if (iph >= m_NumVolPhases) {
        mv = 0.0;
        ma = tp.molarArea();
        printf("                Molar Volume = %11.5E cm3 gmol-1\n", mv * 1.0E3);
        printf("                Molar Area   = %11.5E cm2 gmol-1\n", ma * 10.);
    } else {
        printf("                Molar Volume = %11.5E cm3 gmol-1\n", mv * 1.0E3);
    }
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
        printf("                Voltage = %g volts\n", tp.electricPotential());
    }
    /*
     * Do specific surface phase printouts
     */
    double radius;
    if (iph >= m_NumVolPhases) {
        isurf = iph - m_NumVolPhases;
        radius = sqrt(surfaceAreaRS_final_[isurf] / (4.0 * Pi * particleNumberToFollow_));
        radius *= 1.0E6;
        printf("                surface area (final) = %11.5E m**2      Radius(final) = %11.5E um\n",
               surfaceAreaRS_final_[isurf], radius);
        if (subTimeStep) {
            radius = sqrt(surfaceAreaRS_init_[isurf] / (4.0 * Pi * particleNumberToFollow_));
            radius *= 1.0E6;
            printf("                surface area (init)  = %11.5E m**2      Radius(init)  = %11.5E um\n",
                   surfaceAreaRS_init_[isurf], radius);
        } else {
            radius = sqrt(surfaceAreaRS_init_init_[isurf] / (4.0 * Pi * particleNumberToFollow_));
            radius *= 1.0E6;
            printf("                surface area (init)  = %11.5E m**2      Radius(init)  = %11.5E um\n",
                   surfaceAreaRS_init_init_[isurf], radius);
        }
        int ttt = isExternalSurface_[isurf];
        printf("                IsExternalSurface = %d\n", ttt);
        ttt = ActiveKineticsSurf_[isurf];
        printf("                ActiveKineticsSurf = %d\n", ttt);
        //double oc = openCircuitVoltage(isurf);
        double oc = openCircuitVoltage(isurf);
        if (oc != 0.0) {
            printf("               Open Circuit Voltage = %g Volts\n", oc);
        }
    }
    if (printLvl_ >= 3) {
        printf("\n");
        printf("                Name              MoleFrac_final kMoles_final kMoles_init SrcTermIntegrated(kmol)\n");
        for (size_t k = 0; k < nsp; k++) {
            std::string sname = tp.speciesName(k);
            if (pSrc) {
                if (subTimeStep) {
                    printf("                %-22s %10.3E %10.3E   %10.3E  %10.3E\n", sname.c_str(),
                           spMf_final_[istart + k], spMoles_final_[istart + k], spMoles_init_[istart + k],
                           spMoleIntegratedSourceTermLast_[istart + k]);
                } else {
                    printf("                %-22s %10.3E %10.3E   %10.3E  %10.3E\n", sname.c_str(),
                           spMf_final_[istart + k], spMoles_final_[istart + k], spMoles_init_init_[istart + k],
                           spMoleIntegratedSourceTerm_[istart + k]);
                }
            } else {
                if (subTimeStep) {
                    printf("                %-22s %10.3E %10.3E   %10.3E\n", sname.c_str(), spMf_final_[istart + k],
                           spMoles_final_[istart + k], spMoles_init_[istart + k]);
                } else {
                    printf("                %-22s %10.3E %10.3E   %10.3E\n", sname.c_str(), spMf_final_[istart + k],
                           spMoles_final_[istart + k], spMoles_init_init_[istart + k]);
                }
            }
        }
    }
    if (printLvl_ >= 4) {
        if (iph >= m_NumVolPhases) {
            const std::vector<double>& rsSpeciesProductionRates = RSD_List_[isurf]->calcNetSurfaceProductionRateDensities();
            RSD_List_[isurf]->getNetRatesOfProgress(netROP);

            double* spNetProdPerArea = (double*) spNetProdPerArea_List_.ptrColumn(isurf);
            std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
            size_t nphRS = RSD_List_[isurf]->nPhases();
            size_t kIndexKin = 0;
            for (size_t kph = 0; kph < nphRS; kph++) {
                size_t jph = RSD_List_[isurf]->kinOrder[kph];
                size_t istart = m_PhaseSpeciesStartIndex[jph];
                size_t nsp = m_PhaseSpeciesStartIndex[jph + 1] - istart;
                for (size_t k = 0; k < nsp; k++) {
                    spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
                    kIndexKin++;
                }
            }
            printf("\n");
            printf("                           spName                  SourceRateLastStep (kmol/m2/s) \n");
            for (size_t k = 0; k < m_NumTotSpecies; k++) {
                std::string ss = speciesName(k);
                printf("                           %-22s %10.3E\n", ss.c_str(), spNetProdPerArea[k]);
            }
        }
    }
    printf("     ============================================================================================\n");
    delete[] netROP;
}
//====================================================================================================================
// Determines the level of printing for each step.
/*
 *   0 -> absolutely nothing is printed for a single call to integrate.
 *   1 -> One line summary per integrate call
 *   2 -> short description, points of interest: Table of nonlinear solve - one line per iteration
 *   3 -> Table is included -> More printing per nonlinear iteration (default) that occurs during the table
 *   4 -> Summaries of the nonlinear solve iteration as they are occurring -> table no longer printed
 *   5 -> Algorithm information on the nonlinear iterates are printed out
 *   6 -> Additional info on the nonlinear iterates are printed out
 *   7 -> Additional info on the linear solve is printed out.
 *   8 -> Info on a per iterate of the linear solve is printed out.
 */
void Electrode::setPrintLevel(int printLvl)
{
    printLvl_ = printLvl;
}
//====================================================================================================================
//   Set the deltaT used for the subcycle step
/*
 *  The default setting for this is infinite. If it's infinite then deltaTSubcycle is set to the
 *  deltaT of the integration step. However, there are many times where a smaller natural delta T is
 *  required to solve the equation system
 */
void Electrode::setDeltaTSubcycle(double deltaTsubcycle)
{
    deltaTsubcycleNext_ = deltaTsubcycle;
    deltaTsubcycle_init_next_ = deltaTsubcycle;
    deltaTsubcycle_init_init_ = deltaTsubcycle;
}
//====================================================================================================================
double Electrode::deltaTsubcycle_init_init() const
{
    return deltaTsubcycle_init_init_;
}
//====================================================================================================================
void Electrode::setDeltaTSubcycleMax(double deltaTsubcycle)
{
    deltaTsubcycleMax_ = deltaTsubcycle;
}
//====================================================================================================================
// Write out CSV tabular data on the integrations
/*
 *  The idea is to print out tabular data about each intermediate run and about each
 *  global run
 *
 *  @param itype Type of the data
 *            - 0      Initialization information - overwrite file
 *            - 1      Normal intermediate information
 *            - 2      Normal global information at the end of a global time step
 *            - -1     Failed intermediate calculation, a failure from a nonlinear solver step
 *            - -2     Failed calculation from the predictor step - not necessarily significant.
 */
void Electrode::writeCSVData(int itype)
{
    //if (printLvl_ < 2) return;
    //int k;
    int kstart;
    static std::string globOutputName = "";
    static std::string intOutputName = "";
    static FILE* fpG = 0;
    static FILE* fpI = 0;
    static int runNumber = -1;
    bool ignoreErrors = false;
    size_t isk = 0;  // reacting surface assumption
    if (itype < 0) {
        ignoreErrors = true;
    }

    if (itype == 0 || fpI == 0) {
        runNumber++;
        if (runNumber == 0) {
            intOutputName = (electrodeName_ + "_intResults_" + int2str(electrodeDomainNumber_) + "_"
                    + int2str(electrodeCellNumber_) + ".csv");
        } else {
            intOutputName = (electrodeName_ + "_intResults_" + int2str(electrodeDomainNumber_) + "_"
                    + int2str(electrodeCellNumber_) + "_" + int2str(runNumber) + ".csv");
        }
        fpI = fopen(intOutputName.c_str(), "w");

        fprintf(fpI, "         Tinit ,        Tfinal ,");

        fprintf(fpI, "      Volts_Soln ,  Volts_Electrode ,");
        fprintf(fpI, " OpenCircVolt ,      OCV_avg ,");
        fprintf(fpI, " OCV_RefElec  , OCV_avg_refE ,");

        fprintf(fpI, "      Current ,");

        fprintf(fpI, " RelExtent    ,");
        fprintf(fpI, " DoD_Frac     ,");
        fprintf(fpI, "  CapDischarged ,");
        fprintf(fpI, "  CapacityLeft ,");
        fprintf(fpI, "  Capacity ,");

        for (size_t k = 0; k < m_NumTotSpecies; k++) {
            std::string sss = speciesName(k);
            fprintf(fpI, " MN_%-20.20s,", sss.c_str());
        }

        for (size_t k = 0; k < m_NumTotSpecies; k++) {
            std::string sss = speciesName(k);
            fprintf(fpI, " SRC_%-20.20s,", sss.c_str());
        }

        fprintf(fpI, " iType,");

        fprintf(fpI, " SolidVol,");
        fprintf(fpI, " ElectrolyteVol,");
        fprintf(fpI, " TotalVol,");
        fprintf(fpI, " Porosity");

        fprintf(fpI, "\n");
        fclose(fpI);
    }

    if (itype == 0 || fpG == 0) {
        if (runNumber == 0) {
            globOutputName = (electrodeName_ + "_globalResults_" + int2str(electrodeDomainNumber_) + "_"
                    + int2str(electrodeCellNumber_) + ".csv");
        } else {
            intOutputName = (electrodeName_ + "_globalResults_" + int2str(electrodeDomainNumber_) + "_"
                    + int2str(electrodeCellNumber_) + "_" + int2str(runNumber) + ".csv");
        }
        fpG = fopen(globOutputName.c_str(), "w");
        fprintf(fpG, "         Tinit ,        Tfinal ,");

        fprintf(fpG, "      Volts_Soln ,  Volts_Electrode ,");
        fprintf(fpG, " OpenCircVolt ,      OCV_avg ,");
        fprintf(fpG, " OCV_RefElec  , OCV_avg_refE ,");

        fprintf(fpG, "      Current ,");

        fprintf(fpG, " RelExtent    ,");
        fprintf(fpG, " DoD_Frac     ,");
        fprintf(fpG, "  CapDischarged ,");
        fprintf(fpG, "  CapacityLeft ,");
        fprintf(fpG, "  Capacity ,");

        for (size_t k = 0; k < m_NumTotSpecies; k++) {
            std::string sss = speciesName(k);
            fprintf(fpG, " MN_%-20.20s,", sss.c_str());
        }

        for (size_t k = 0; k < m_NumTotSpecies; k++) {
            std::string sss = speciesName(k);
            fprintf(fpG, " SRC_%-20.20s,", sss.c_str());
        }

        fprintf(fpG, " iType,");

        fprintf(fpG, " SolidVol,");
        fprintf(fpG, " ElectrolyteVol,");
        fprintf(fpG, " TotalVol,");
        fprintf(fpG, " Porosity");
        fprintf(fpG, "\n");
        fclose(fpG);
    }

    // Provide a hook here for the debugger
    if (itype == -1) {
        kstart = m_PhaseSpeciesStartIndex[metalPhase_];
    }

    if (itype == 1 || itype == 2) {
        fpI = fopen(intOutputName.c_str(), "a");
        fprintf(fpI, "  %12.5E ,  %12.5E ,", tinit_, tfinal_);

        double phiS = phaseVoltages_[solnPhase_];
        double phiM = phaseVoltages_[metalPhase_];
        fprintf(fpI, "    %12.5E ,     %12.5E ,", phiS, phiM);

        double ocv_avg = openCircuitVoltage_MixtureAveraged(isk, false);
        double ocv = openCircuitVoltage(isk, false);
        fprintf(fpI, " %12.5E , %12.5E , ", ocv, ocv_avg);

        double ocv_avg_refE = openCircuitVoltage_MixtureAveraged(isk, true);
        double ocv_refE = openCircuitVoltage(isk, true);
        fprintf(fpI, " %12.5E , %12.5E , ", ocv_refE, ocv_avg_refE);

        kstart = m_PhaseSpeciesStartIndex[metalPhase_];
        double deltaT = tfinal_ - tinit_;
        double amps = 0.0;
        if (deltaT > 1.0E-200) {
            amps = spMoleIntegratedSourceTermLast_[kstart] / deltaT * Faraday;
        }
        fprintf(fpI, " %12.5E ,", amps);

        double relE = calcRelativeExtentRxn_final();
        fprintf(fpI, " %12.5E ,", relE);

        double dodF = depthOfDischargeFraction(-1);
        fprintf(fpI, " %12.5E ,", dodF);

        double cap = capacityDischarged();
        fprintf(fpI, "   %12.5E ,", cap);

        cap = capacityLeft();
        fprintf(fpI, "   %12.5E ,", cap);

        cap = capacity();
        fprintf(fpI, "   %12.5E ,", cap);

        for (size_t k = 0; k < m_NumTotSpecies; k++) {
            fprintf(fpI, " %12.5E           ,", spMoles_final_[k]);
        }

        for (size_t k = 0; k < m_NumTotSpecies; k++) {
            fprintf(fpI, " %12.5E            ,", spMoleIntegratedSourceTermLast_[k]);
        }

        fprintf(fpI, "   %3d ,", itype);

        double solidVol = SolidVol();
        double grossVol = TotalVol(ignoreErrors);
        double lyteVol = grossVol - solidVol;
        double porosity = lyteVol / grossVol;
        //  This will fail until we figure out how we want to handle this
        // if (fabs(porosity - porosity_) > 1.0E-6) {
        //	throw Electrode_Error("writeCSVData() ERROR", "porosity isn't consistent");
        //}
        fprintf(fpI, " %12.5E ,", solidVol);
        fprintf(fpI, " %12.5E ,", lyteVol);
        fprintf(fpI, " %12.5E ,", grossVol);
        fprintf(fpI, " %12.5E ", porosity);
        fprintf(fpI, " \n");

        fclose(fpI);
    }

    if (itype == 2) {
        fpG = fopen(globOutputName.c_str(), "a");
        fprintf(fpG, "  %12.5E ,  %12.5E ,", t_init_init_, t_final_final_);

        double phiS = phaseVoltages_[solnPhase_];
        double phiM = phaseVoltages_[metalPhase_];
        fprintf(fpG, "    %12.5E ,     %12.5E ,", phiS, phiM);

        double ocv_avg = openCircuitVoltage_MixtureAveraged(isk);
        double ocv = openCircuitVoltage(isk);
        fprintf(fpG, " %12.5E , %12.5E , ", ocv, ocv_avg);

        double ocv_avg_refE = openCircuitVoltage_MixtureAveraged(isk, true);
        double ocv_refE = openCircuitVoltage(isk, true);
        fprintf(fpG, " %12.5E , %12.5E , ", ocv_refE, ocv_avg_refE);

        kstart = m_PhaseSpeciesStartIndex[metalPhase_];
        double deltaT = t_final_final_ - t_init_init_;
        double amps = 0.0;
        if (deltaT > 1.0E-200) {
            amps = spMoleIntegratedSourceTerm_[kstart] / deltaT * Faraday;
        }
        fprintf(fpG, " %12.5E ,", amps);

        double relE = calcRelativeExtentRxn_final();
        fprintf(fpG, " %12.5E ,", relE);

        double dodF = depthOfDischargeFraction(-1);
        fprintf(fpG, " %12.5E ,", dodF);

        double cap = capacityDischarged();
        fprintf(fpG, "   %12.5E ,", cap);

        cap = capacityLeft();
        fprintf(fpG, "   %12.5E ,", cap);

        cap = capacity();
        fprintf(fpG, "   %12.5E ,", cap);

        for (size_t k = 0; k < m_NumTotSpecies; k++) {
            fprintf(fpG, " %12.5E           ,", spMoles_final_[k]);
        }

        for (size_t k = 0; k < m_NumTotSpecies; k++) {
            fprintf(fpG, " %12.5E            ,", spMoleIntegratedSourceTerm_[k]);
        }

        fprintf(fpG, "   %3d,", itype);
        double solidVol = SolidVol();
        double grossVol = TotalVol();
        double lyteVol = grossVol - solidVol;
        double porosity = lyteVol / grossVol;

        fprintf(fpG, " %12.5E ,", solidVol);
        fprintf(fpG, " %12.5E ,", lyteVol);
        fprintf(fpG, " %12.5E ,", grossVol);
        fprintf(fpG, " %12.5E ", porosity);
        fprintf(fpG, " \n");
        fclose(fpG);
    }
}
//==================================================================================================================================
}// End of namespace
//----------------------------------------------------------------------------------------------------------------------------------
