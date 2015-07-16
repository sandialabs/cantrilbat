/*
 * $Id: EState_RadialDiffusion.cpp 571 2013-03-26 16:44:21Z hkmoffa $
 */


#include "cantera/base/ctml.h"

#include "EState_RadialDistrib.h"

#include "Electrode_SimpleDiff.h"
#include "Electrode_DiffTALE.h"
#include "EState_XML.h"

#include <cstdio>
#include <cstdlib>
using namespace Cantera;
using namespace std;

namespace Cantera
{
//======================================================================================================================
/*
 * EState constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.
 */
EState_RadialDistrib::EState_RadialDistrib() :
    EState(),
    numRCells_(0),
    numKRSpecies_(0),
    numSPhases_(0),
    rnodePos_(0),
    cellBoundR_(0),
    concTot_SPhase_Cell_(0),
    concKRSpecies_Cell_(0),
    spMoles_KRsolid_Cell_(0),
    onRegionBoundary_(-1)
{
    EST_fileToBeWritten_ = EST_RADIALDISTRIB;
    electrodeTypeString_ = "SimpleDiff";
}
//======================================================================================================================
/*
 * EState constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.
 */
EState_RadialDistrib::EState_RadialDistrib(std::string electrodeType) :
    EState(),
    numRCells_(0),
    numKRSpecies_(0),
    numSPhases_(0),
    rnodePos_(0),
    cellBoundR_(0),
    concTot_SPhase_Cell_(0),
    concKRSpecies_Cell_(0),
    spMoles_KRsolid_Cell_(0),
    onRegionBoundary_(-1)
{
    EST_fileToBeWritten_ = EST_RADIALDISTRIB;
    if (electrodeType == "SimpleDiff") {
        electrodeTypeString_ = "SimpleDiff";
    } else if  (electrodeType == "SimpleDiff") {
        electrodeTypeString_ = "DiffTALE";
    } else {
         throw Electrode_Error("EState_RadialDistrib::EState_RadialDistrib()", "wrong electrode type");  
    }
}
//======================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
EState_RadialDistrib::EState_RadialDistrib(const EState_RadialDistrib& right) :
    EState(),
    numRCells_(0),
    numKRSpecies_(0),
    numSPhases_(0),
    rnodePos_(0),
    cellBoundR_(0),
    concTot_SPhase_Cell_(0),
    concKRSpecies_Cell_(0),
    spMoles_KRsolid_Cell_(0),
    onRegionBoundary_(-1)
{
    EST_fileToBeWritten_ = EST_MULTIPLATEAU;

    /*
     * Call the assignment operator.
     */
    EState_RadialDistrib::operator=(right);
}
//======================================================================================================================
// Assignment operator
/*
 *  @param right object to be copied
 */
EState_RadialDistrib& EState_RadialDistrib::operator=(const EState_RadialDistrib& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    EState::operator=(right);

    numRCells_                             = right.numRCells_;
    numKRSpecies_                          = right.numKRSpecies_;
    numSPhases_                            = right.numSPhases_;
    rnodePos_                              = right.rnodePos_;
    cellBoundR_                            = right.cellBoundR_;
    concTot_SPhase_Cell_                   = right.concTot_SPhase_Cell_;
    concKRSpecies_Cell_                    = right.concKRSpecies_Cell_;
    spMoles_KRsolid_Cell_                  = right.spMoles_KRsolid_Cell_;
    onRegionBoundary_                      = right.onRegionBoundary_;
    /*
     * Return the reference to the current object
     */
    return *this;
}
//======================================================================================================================
EState_RadialDistrib::~EState_RadialDistrib()
{
}
//======================================================================================================================
// Duplicator function for this class
/*
 *  @return Returns a duplication of the current state as a pointer to the base class
 */
EState* EState_RadialDistrib::duplMyselfAsEState() const
{
    EState_RadialDistrib* es = new EState_RadialDistrib(*this);
    return dynamic_cast<EState*>(es);
}
//======================================================================================================================
int EState_RadialDistrib::initialize(const Cantera::Electrode* const ebase)
{
    EState::initialize(ebase);

    const Electrode_SimpleDiff* const e = dynamic_cast<const Electrode_SimpleDiff* const>(ebase);
    if (e) {
	numRCells_                   = e->numRCells_;
	numKRSpecies_                = e->numKRSpecies_;
	numSPhases_                  = e->numSPhases_;
	rnodePos_                    = e->rnodePos_final_;
	cellBoundR_                  = e->cellBoundR_final_;
	concTot_SPhase_Cell_         = e->concTot_SPhase_Cell_final_;
	concKRSpecies_Cell_          = e->concKRSpecies_Cell_final_;
	spMoles_KRsolid_Cell_        = e->spMoles_KRsolid_Cell_final_;
    } else {
        const Electrode_DiffTALE* const edt = dynamic_cast<const Electrode_DiffTALE* const>(ebase);
	if (edt) {
	    numRCells_                   = edt->numRCells_;
	    numKRSpecies_                = edt->numKRSpecies_;
	    numSPhases_                  = edt->numSPhases_;
	    rnodePos_                    = edt->rnodePos_final_;
	    cellBoundR_                  = edt->cellBoundR_final_;
	    concTot_SPhase_Cell_         = edt->concTot_SPhase_Cell_final_;
	    concKRSpecies_Cell_          = edt->concKRSpecies_Cell_final_;
	    spMoles_KRsolid_Cell_        = edt->spMoles_KRsolid_Cell_final_;
	} else {
	    throw Electrode_Error("EState_RadialDistrib::initialize()",
                                  "Need Electrode_SimpleDiff or Electrode_DiffTALE pointer type");
        }
    }
    return 1;
}
//======================================================================================================================
// Write the ElectrodeState to an XML_Node tree
/*
 *  @return pointer to the XML_Node tree
 */
XML_Node*   EState_RadialDistrib::writeIdentificationToXML() const
{
    XML_Node* x = new XML_Node("ElectrodeIdentification");

    ctml::addNamedString(*x, "electrodeTypeString", electrodeTypeString_);
    ctml::addInteger(*x, "EState_Type",         EST_fileToBeWritten_);
    ctml::addNamedString(*x, "EState_Type_String", esmodel::EState_Type_Enum_to_string(EST_fileToBeWritten_));
    ctml::addInteger(*x, "fileVersionNumber",  EST_Version_lastFileRead_);
    ctml::addInteger(*x, "electrodeModelType",  electrodeChemistryModelType_);
    ctml::addInteger(*x, "electrodeDomainNumber",  electrodeDomainNumber_);
    ctml::addInteger(*x, "electrodeCellNumber",  electrodeCellNumber_);
   if (electrodeCapacityType_ == CAPACITY_ANODE_ECT) {
	ctml::addNamedString(*x, "electrodeCapacityType", "Capacity_Anode");
    } else if (electrodeCapacityType_ == CAPACITY_CATHODE_ECT) {
	ctml::addNamedString(*x, "electrodeCapacityType", "Capacity_Cathode");
    } else {
	ctml::addNamedString(*x, "electrodeCapacityType", "Capacity_Other");
    }

    return x;
}
//======================================================================================================================
// Write the ElectrodeState to an XML_Node tree
/*
 *  @return pointer to the XML_Node tree
 */
XML_Node* EState_RadialDistrib::write_electrodeState_ToXML() const
{
    XML_Node* x = EState::write_electrodeState_ToXML();

    ctml::addInteger(*x, "numRCells", numRCells_);
    ctml::addInteger(*x, "numKRSpecies", numKRSpecies_);
    ctml::addInteger(*x, "numSPhases", numSPhases_);

    ctml::addNamedFloatArray(*x, "rnodePos", numRCells_, DATA_PTR(rnodePos_), "m");
    ctml::addNamedFloatArray(*x, "cellBoundR", numRCells_, DATA_PTR(cellBoundR_), "m");
    ctml::addNamedFloatArray(*x, "concTot_SPhase_Cell", numRCells_ * numSPhases_, DATA_PTR(concTot_SPhase_Cell_), "kmol/m3");
    ctml::addNamedFloatArray(*x, "concKRSpecies_Cell", numRCells_ * numKRSpecies_, DATA_PTR(concKRSpecies_Cell_), "kmol/m3");
    ctml::addNamedFloatArray(*x, "spMoles_KRsolid_Cell", numRCells_ * numKRSpecies_, DATA_PTR(spMoles_KRsolid_Cell_), "kmol");

    return x;
}
//======================================================================================================================
//  Read the state from the XML_Node  given by the argument
/*
 *  @return pointer to the XML_Node tree
 */
void EState_RadialDistrib::readStateFromXML(const XML_Node& xmlEState)
{
    EState::readStateFromXML(xmlEState);

    numRCells_ = ctml::getInteger(xmlEState, "numRCells");
    numKRSpecies_ = ctml::getInteger(xmlEState, "numKRSpecies");
    numSPhases_ = ctml::getInteger(xmlEState, "numSPhases");

    ctml::getFloatArray(xmlEState, rnodePos_, true, "", "rnodePos");
    ctml::getFloatArray(xmlEState, cellBoundR_, true, "", "cellBoundR"); 

    ctml::getFloatArray(xmlEState, concTot_SPhase_Cell_, true, "", "concTot_SPhase_Cell");
    ctml::getFloatArray(xmlEState, concKRSpecies_Cell_, true, "", "concKRSpecies_Cell");
    ctml::getFloatArray(xmlEState, spMoles_KRsolid_Cell_, true, "", "spMoles_KRsolid_Cell");
}
//======================================================================================================================
void EState_RadialDistrib::copyElectrode_SimpleDiff_intoState(const Cantera::Electrode_SimpleDiff* const emp)
{
    EState::copyElectrode_intoState(emp);

    rnodePos_                    = emp->rnodePos_final_;
    cellBoundR_                  = emp->cellBoundR_final_;
    concTot_SPhase_Cell_         = emp->concTot_SPhase_Cell_final_;
    concKRSpecies_Cell_          = emp->concKRSpecies_Cell_final_;
    spMoles_KRsolid_Cell_        = emp->spMoles_KRsolid_Cell_final_;
}
//======================================================================================================================
void EState_RadialDistrib::copyElectrode_DiffTALE_intoState(const Cantera::Electrode_DiffTALE* const emp)
{
    EState::copyElectrode_intoState(emp);

    rnodePos_                    = emp->rnodePos_final_;
    cellBoundR_                  = emp->cellBoundR_final_;
    concTot_SPhase_Cell_         = emp->concTot_SPhase_Cell_final_;
    concKRSpecies_Cell_          = emp->concKRSpecies_Cell_final_;
    spMoles_KRsolid_Cell_        = emp->spMoles_KRsolid_Cell_final_;
}
//======================================================================================================================
// Set the State of this object from the state of the Electrode object
/*
 *  (virtual function)
 *
 *  virtual function, because the base class is called from general code, allowing the child classes to be invoked.
 *
 *  This function takes the electrode objects _final_ state and copies it into this object.
 *  This function must be carried out before the XML_Node tree is written.
 *
 *  @param e   Pointer to the Electrode object. Note, this class may use dynamic casting
 *             to choose a child object, and then may invoke an error if the match isn't
 *             correct.
 */
void EState_RadialDistrib::copyElectrode_intoState(const Cantera::Electrode* const e)
{
    const Cantera::Electrode_SimpleDiff* const emp = dynamic_cast<const Cantera::Electrode_SimpleDiff* const>(e);
    if (emp) {
        copyElectrode_SimpleDiff_intoState(emp);
    } else {
       const Cantera::Electrode_DiffTALE* const edt = dynamic_cast<const Cantera::Electrode_DiffTALE* const>(e);
       if (edt) {
          copyElectrode_DiffTALE_intoState(edt);
       } else {
          throw CanteraError("EState_RadialDistrib::copyElectrode_intoState()","bad cast");
       }
    }
}
//======================================================================================================================
//    Set the state of the Electrode from the state of this object
void EState_RadialDistrib::setStateElectrode_SimpleDiff_fromEState(Cantera::Electrode_SimpleDiff* const emp) const
{
    EState::copyEState_toElectrode(emp);

    emp->rnodePos_final_             = 	rnodePos_;    
    emp->cellBoundR_final_           = 	cellBoundR_;
    emp->concTot_SPhase_Cell_final_  = 	concTot_SPhase_Cell_;
    emp->concKRSpecies_Cell_final_   = 	concKRSpecies_Cell_;
    emp->spMoles_KRsolid_Cell_final_ = 	spMoles_KRsolid_Cell_;
    /*
     * Now we can do an update
     */
    emp->updateState();
    emp->stateToPhaseFlagsReconciliation(false);
    emp->setInitStateFromFinal(true);
}
//======================================================================================================================
//    Set the state of the Electrode from the state of this object
void EState_RadialDistrib::setStateElectrode_DiffTALE_fromEState(Cantera::Electrode_DiffTALE* const emp) const
{
    EState::copyEState_toElectrode(emp);

    emp->rnodePos_final_             = 	rnodePos_;    
    emp->cellBoundR_final_           = 	cellBoundR_;
    emp->rLatticeCBR_final_          =  cellBoundR_;
    emp->concTot_SPhase_Cell_final_  = 	concTot_SPhase_Cell_;
    emp->concKRSpecies_Cell_final_   = 	concKRSpecies_Cell_;
    emp->spMoles_KRsolid_Cell_final_ = 	spMoles_KRsolid_Cell_;
    /*
     * Now we can do an update
     */
    emp->updateState();
    emp->stateToPhaseFlagsReconciliation(false);
    emp->setInitStateFromFinal(true);
}
//======================================================================================================================
//    Set the state of the Electrode from the state of this object
void EState_RadialDistrib::setStateElectrode_fromEState(Cantera::Electrode* const e) const
{
    Electrode_SimpleDiff* emp = dynamic_cast<Electrode_SimpleDiff*>(e);
    if (emp) {
	setStateElectrode_SimpleDiff_fromEState(emp); 
    } else {
       Cantera::Electrode_DiffTALE* edt = dynamic_cast<Cantera::Electrode_DiffTALE*>(e);
       if (edt) {
	   setStateElectrode_DiffTALE_fromEState(edt);
       } else {
	   throw CanteraError("EState_RadialDistrib::setStateElectrode_fromEState()","bad cast");
       }
    }
}
//======================================================================================================================
bool EState_RadialDistrib::compareOtherState(const EState* const ESguest, double molarAtol, int nDigits, 
			               	     bool includeHist, int printLvl) const
{
     bool btotal = true;
     bool boolR = true;
     const EState_RadialDistrib* ESG_rad = dynamic_cast<const EState_RadialDistrib*>(ESguest);
     if (!ESG_rad) {
          throw Electrode_Error("EState_RadialDistrib::compareOtherState()",
                                "Guest state isn't an EState_RadialDistrib object. This hasn't been done yet");
     }

     printHead(-1);

     if (printLvl > 1) {
	 printf("\t\tEState_RadialDistrib::compareOtherState() Start comparison of types %s vs. %s\n",
		electrodeTypeString_.c_str(), ESguest->electrodeTypeString_.c_str());
     }

     /*
      *  Compare the 0D outputs first. They must agree first
      */
     bool sub_ok = EState::compareOtherState(ESguest, molarAtol, nDigits, includeHist, printLvl); 
     btotal = btotal && sub_ok; 


     double volAtol =  molarAtol / 55.;
     double radiusAtol = pow(volAtol, 0.3333);
     // double surfaceAtol = 12 * radiusAtol * radiusAtol;

     boolR = doubleVectorEqual(rnodePos_, ESG_rad->rnodePos_, radiusAtol, nDigits);
     if (!boolR) {
	 printVecDiff("rnodeR_", rnodePos_, ESG_rad->rnodePos_, printLvl);
     }
     btotal = boolR && btotal;

     boolR = doubleVectorEqual(cellBoundR_, ESG_rad->cellBoundR_, radiusAtol, nDigits);
     if (!boolR) {
	 printVecDiff("cellboundR_", cellBoundR_, ESG_rad->cellBoundR_, printLvl);
     }
     btotal = boolR && btotal;

     double concTolAtol = molarAtol * 55.;
     boolR = doubleVectorEqual(concTot_SPhase_Cell_, ESG_rad->concTot_SPhase_Cell_, concTolAtol, nDigits);
     if (!boolR) {
	 printVecDiff("concTot_SPhase_Cell_", concTot_SPhase_Cell_, ESG_rad->concTot_SPhase_Cell_, printLvl);
     }
     btotal = boolR && btotal;

     boolR = doubleVectorEqual(concKRSpecies_Cell_, ESG_rad->concKRSpecies_Cell_, concTolAtol, nDigits);
     if (!boolR) {
	 printVecDiff("concKRSpecies_Cell_", concKRSpecies_Cell_, ESG_rad->concKRSpecies_Cell_, printLvl);
     }
     btotal = boolR && btotal;

     boolR = doubleVectorEqual(spMoles_KRsolid_Cell_, ESG_rad->spMoles_KRsolid_Cell_, molarAtol, nDigits);
     if (!boolR) {
	 printVecDiff("spMoles_KRsolid_Cell_", spMoles_KRsolid_Cell_, ESG_rad->spMoles_KRsolid_Cell_, printLvl);
     }
     btotal = boolR && btotal;

     boolR = onRegionBoundary_ == ESG_rad->onRegionBoundary_;
     if (!boolR) {
	 printDiff("onRegionBoundary_", true, onRegionBoundary_, ESG_rad->onRegionBoundary_, printLvl);
     }
     btotal = boolR && btotal;

     if (!btotal) {
	 if (printLvl == 1) {
	     printf("\t\tEState_RadialDistrib:: ERROR:   States are not the same\n");
	 }
     } else {
	 if (printLvl >= 1) {
	     printf("\t\tEState_RadialDistrib:: SUCCESS: States are the same\n");
	 }
     }

     return btotal;
}




//====================================================================================================================
} // End of namespace Cantera
//======================================================================================================================
