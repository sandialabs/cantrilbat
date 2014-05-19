/**
 * @file ReactingVolDomain.cpp
 *
 * $Author: hkmoffa $
 * $Revision: 571 $
 * $Date: 2013-03-26 10:44:21 -0600 (Tue, 26 Mar 2013) $
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "ReactingSurDomain.h"
#include "PhaseList.h"
#include "RxnMolChange.h"
#include "cantera/kinetics.h"
#include "Electrode_input.h"


using namespace std;
namespace Cantera
{
//====================================================================================================================
/*
 * Constructor #1 for the IdealReactingGas object.
 *
 * infile = input file location
 * id = ?
 */
ReactingSurDomain::ReactingSurDomain() :
    InterfaceKinetics(),
    numPhases_(0),
    xmlList(0),
    kinOrder(0),
    PLtoKinPhaseIndex_(0),
    PLtoKinSpeciesIndex_(0),
    iphaseKin_(-1),
    tpList_IDs_(0),
    tplRead(0),
    m_DoSurfKinetics(false),
    speciesProductionRates_(0),
    speciesCreationRates_(0),
    speciesDestructionRates_(0),
    deltaGRxn_(0),
    m_pl(0),
    metalPhaseRS_(-1),
    kElectronRS_(-1),
    solnPhaseRS_(-1),
    ocv_ptr_(0),
    OCVmodel_(0)
{
}
//====================================================================================================================
// Copy Constructor for the %Kinetics object.
/*
 * Currently, this is not fully implemented. If called it will
 * throw an exception.
 */
ReactingSurDomain::ReactingSurDomain(const ReactingSurDomain& right) :
    InterfaceKinetics(),
    numPhases_(0),
    xmlList(0),
    kinOrder(0),
    PLtoKinPhaseIndex_(0),
    PLtoKinSpeciesIndex_(0),
    iphaseKin_(-1),
    tpList_IDs_(0),
    tplRead(0),
    m_DoSurfKinetics(false),
    speciesProductionRates_(0),
    speciesCreationRates_(0),
    speciesDestructionRates_(0),
    deltaGRxn_(0),
    m_pl(0),
    metalPhaseRS_(-1),
    kElectronRS_(-1),
    solnPhaseRS_(-1),
    ocv_ptr_(0),
    OCVmodel_(0)

{
    /*
     * Call the assignment operator
     */
    operator=(right);
}
//====================================================================================================================
// Assignment operator
/*
 *  This is NOT a virtual function.
 *
 * @param right    Reference to %Kinetics object to be copied into the
 *                 current one.
 */
ReactingSurDomain&  ReactingSurDomain::operator=(const ReactingSurDomain& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }
    //
    //  Beware: The copy operation within InterfaceKinetics leaves shallow pointers to the
    //          underlying ThermoPhase classes in place. They must be fixed up at the Electrode
    //          object level where we copy the ThermoPhase classes.
    //
    InterfaceKinetics::operator=(right);

    numPhases_      = right.numPhases_;
    // Shallow copy of xmlList pointers -> beware
    xmlList         = right.xmlList;
    kinOrder        = right.kinOrder;
    PLtoKinPhaseIndex_ = right.PLtoKinPhaseIndex_;
    PLtoKinSpeciesIndex_ = right.PLtoKinSpeciesIndex_;
    iphaseKin_      = right.iphaseKin_;
    tpList_IDs_     = right.tpList_IDs_;
    tplRead         = right.tplRead;
    m_DoSurfKinetics = right.m_DoSurfKinetics;
    speciesProductionRates_ = right.speciesProductionRates_;
    speciesCreationRates_ = right.speciesCreationRates_;
    speciesDestructionRates_ = right.speciesDestructionRates_;
    deltaGRxn_ = right.deltaGRxn_;

    //
    // Beware -  Shallow copy of m_pl pointer
    //
    m_pl            = right.m_pl;

    metalPhaseRS_ = right.metalPhaseRS_;
    kElectronRS_ = right.kElectronRS_;
    solnPhaseRS_ = right.solnPhaseRS_;

    for (int i = 0; i < (int)rmcVector.size(); i++) {
        delete rmcVector[i];
        rmcVector[i] = 0;
    }
    rmcVector.resize(m_ii,0);
    for (size_t i = 0; i < m_ii; i++) {
        if (right.rmcVector[i]) {
            rmcVector[i] = new RxnMolChange(*(right.rmcVector[i]));
        }
    }

    if (right.ocv_ptr_) {
	delete ocv_ptr_;
	ocv_ptr_ =  new OCV_Override_input(*right.ocv_ptr_);
    }

    if (right.OCVmodel_) {
        delete OCVmodel_;
        //
        //  BEWARE: contains pointers to ThermoPhase object that will be needed to be updated
        //          to the owning ThermoPhase object
        //
        OCVmodel_ = new RSD_OCVmodel(*right.OCVmodel_);
    }


    return *this;
}
//====================================================================================================================
/*
 * Destructor for the ReactingSurDomain object.
 *
 * We must decide whether this object owns its own xml tree
 * structure.
 */
ReactingSurDomain::~ReactingSurDomain()
{
    int nr = nReactions();
    for (int i = 0; i < nr; i++) {
        delete rmcVector[i];
    }
}
//====================================================================================================================
// Duplication routine for objects which inherit from
// Kinetics
/*
 *  This virtual routine can be used to duplicate %ReactingSurDomain objects
 *  inherited from %Kinetics even if the application only has
 *  a pointer to %Kinetics to work with.
 *
 *  These routines are basically wrappers around the derived copy
 *  constructor.
 *
 * @param  tpVector Vector of shallow pointers to ThermoPhase objects. this is the
 *                  m_thermo vector within this object
 */
Kinetics* ReactingSurDomain::duplMyselfAsKinetics(const std::vector<thermo_t*>& tpVector) const
{
    ReactingSurDomain* rsd = new ReactingSurDomain(*this);
    rsd->assignShallowPointers(tpVector);
    return dynamic_cast<Kinetics*>(rsd);
}
//====================================================================================================================
// Returns a reference to the calculated production rates of species
/*
 *   This routine calls thet getNetProductionRate function
 *   and then returns a reference to the result.
 *
 * @return Vector of length m_kk containing the species net (kmol s-1 m-2)
 *         production rates
 */
const std::vector<double>& ReactingSurDomain::calcNetSurfaceProductionRateDensities()
{
    getNetProductionRates(&speciesProductionRates_[0]);
    return speciesProductionRates_;
}
//====================================================================================================================
const std::vector<double>& ReactingSurDomain::calcNetSurfaceROP()
{
    updateROP();
    return m_ropnet;
}
//====================================================================================================================
//    Returns a reference to the calculated creation rates of species
/*
 *   This routine calls thet getCreationRate function
 *   and then returns a reference to the result.
 *
 * @return Vector of length m_kk containing the species creation rates (kmol s-1 m-2)
 */
const std::vector<double>& ReactingSurDomain::calcSurfaceCreationRateDensities()
{
    getCreationRates(&speciesCreationRates_[0]);
    return speciesCreationRates_;
}
//====================================================================================================================
// Returns a reference to the calculated destruction rates of species
/*
 *   This routine calls thet getDestructionRate function
 *   and then returns a reference to the result.
 *
 * @return Vector of length m_kk containing the species destruction rates (kmol s-1 m-2)
 */
const std::vector<double>& ReactingSurDomain::calcSurfaceDestructionRateDensities()
{
    getDestructionRates(&speciesDestructionRates_[0]);
    return speciesDestructionRates_;
}
//====================================================================================================================
double ReactingSurDomain::getCurrentDensityRxn(double *currentDensityRxn) 
{
    double netCurrentDensity = 0.0;
    double ps, rs;
    if (kElectronRS_< 0) {
	return netCurrentDensity;
    }
    size_t nr = nReactions();
    // update rates of progress -> puts this into m_ropnet[]
    updateROP();
    if (currentDensityRxn) {
	for (size_t irxn = 0; irxn < nr; irxn++) {
	    rs = m_rrxn[kElectronRS_][irxn];
	    ps = m_prxn[kElectronRS_][irxn];
	    double electronProd = (ps - rs) * m_ropnet[irxn];
	    currentDensityRxn[irxn] = Faraday * electronProd;
	    netCurrentDensity += currentDensityRxn[irxn];
	}
    } else {
	for (size_t irxn = 0; irxn < nr; irxn++) {
	    rs = m_rrxn[kElectronRS_][irxn];
	    ps = m_prxn[kElectronRS_][irxn];
	    double electronProd = (ps - rs) * m_ropnet[irxn];
	    netCurrentDensity +=  Faraday * electronProd;
	}
    }
    return netCurrentDensity;
}
//====================================================================================================================

double ReactingSurDomain::getExchangeCurrentDensityFormulation(int irxn,
        double* nStoich, doublereal* OCV, doublereal* io,
        doublereal* overPotential, doublereal *beta)
{
    // This will calculate the equilibrium constant
    updateROP();

    RxnMolChange*   rmc = rmcVector[irxn];
    // could also get this from reactant and product stoichiometry
    double nStoichElectrons = -  rmc->m_phaseChargeChange[metalPhaseRS_];
    *nStoich =  nStoichElectrons;
    *OCV = 0.0;

    int nr = nReactions();
    std::vector<double> deltaG(nr);
    getDeltaGibbs(DATA_PTR(deltaG));

    if (nStoichElectrons != 0.0) {
        *OCV = deltaG[irxn]/Faraday/ nStoichElectrons;
    }
    //PROBABLY DELETE THIS CALL SINCE IT IS CALLED BY updateROP()
    // we have a vector of standard concentrations calculated from the routine below
    //            m_StandardConc[ik]
    updateExchangeCurrentQuantities();


    //rkc is reciprocal equilibrium constant
#define OLDWAY
#ifdef OLDWAY
    const vector_fp& rf = m_rfn;
    const vector_fp& rkc= m_rkcn;
#else
    const vector_fp& rf = m_kdata->m_rfn;
    const vector_fp& rkc= m_kdata->m_rkcn;
#endif

    // start with the forward reaction rate
    double iO = rf[irxn] * Faraday * nStoichElectrons;

    if (m_beta[irxn] > 0.0) {
        iO *= pow(rkc[irxn], m_beta[irxn]);
    }
    double b = m_beta[irxn];
    *beta = m_beta[irxn];
    double omb = 1.0 - b;

    int m_nTotalSpecies = m_conc.size();

    for (int k = 0; k < m_nTotalSpecies; k++) {
        doublereal reactCoeff = reactantStoichCoeff(k, irxn);
        doublereal prodCoeff =  productStoichCoeff(k, irxn);

        if (reactCoeff != 0.0) {
            iO *= pow(m_conc[k], reactCoeff*omb);
            iO *= pow(m_StandardConc[k], reactCoeff*b);
        }
        if (prodCoeff != 0.0) {
            iO *= pow(m_conc[k], prodCoeff*b);
            iO /= pow(m_StandardConc[k], prodCoeff*omb);
        }
    }
    *io = iO;

    double phiMetal = thermo(metalPhaseRS_).electricPotential();
    double phiSoln = thermo(solnPhaseRS_).electricPotential();
    double E = phiMetal - phiSoln;
    *overPotential = E - *OCV;

    return calcCurrentDensity(*overPotential, *nStoich, *io, *beta, m_temp);
}
//====================================================================================================================
double ReactingSurDomain::calcCurrentDensity(double nu, double nStoich,  double io, double beta, double temp) const
{
     double exp1 = nu * Faraday * beta / (GasConstant * temp);
     double exp2 = -nu * Faraday * (1.0 - beta) / (GasConstant * temp);
     double val = io * (exp(exp1) - exp(exp2));
     return val;
}
//====================================================================================================================
/*
 *  This ostream function describes how to extend cout output
 *  functions to this object. The way this is done is to
 *  call the Cantera's report function, which takes a ThermoPhase
 *  object as its arugment.
 *  This is a "friend" function to the class IdealReactingGas.
 *  Both this function and report are in the Cantera namespace.
 *

 *  Note -> The output doesn't cover kinetics.
 */
std::ostream& operator<<(std::ostream& s,
                         ReactingSurDomain& mix)
{
    ThermoPhase* th;
    InterfaceKinetics* iK = &mix;
    for (int i = 0; i < mix.numPhases_; i++) {
        th = &(iK->thermo(i));
        std::string r = th->report(true);
        s << r;
    }
    return s;
}

//====================================================================================================================
//  Identify the metal phase and the electrons species
void ReactingSurDomain::identifyMetalPhase()
{
    metalPhaseRS_ = -1;
    kElectronRS_ = -1;
    int nr = nReactions();
    int np = nPhases();
    for (int iph = 0; iph < np; iph++) {
        ThermoPhase* tp = & (thermo(iph));
        int nSpecies = tp->nSpecies();
        int nElements = tp->nElements();
        int eElectron = tp->elementIndex("E");
        if (eElectron >= 0) {
            for (int k = 0; k < nSpecies; k++) {
                if (tp->nAtoms(k,eElectron) == 1) {
                    int ifound = 1;
                    for (int e = 0; e < nElements; e++) {
                        if (tp->nAtoms(k,e) != 0.0) {
                            if (e != eElectron) {
                                ifound = 0;
                            }
                        }
                    }
                    if (ifound == 1) {
                        metalPhaseRS_ = iph;
                        kElectronRS_ = m_start[iph] + k;
                    }
                }
            }
        }
        if (iph != metalPhaseRS_) {
            int jph = PLtoKinPhaseIndex_[iph];
            if (jph >= 0) {
                for (int i = 0; i < nr; i++) {
                    RxnMolChange* rmc = rmcVector[i];
                    if (rmc->m_phaseChargeChange[jph] != 0) {
                        if (rmc->m_phaseDims[jph] == 3) {
                            solnPhaseRS_ = iph;
                            break;
                        }
                    }
                }
            }
        }

    }
}
//====================================================================================================================
/*
 *      @param iskin   If this is zero or positive, we will formulate a kinetics mechanism using an 
 *                     interfacial kinetics object.
 */
bool ReactingSurDomain::
importFromPL(Cantera::PhaseList* pl, int iskin)
{
    try {
        int iph;
        //
        //  Store the PhaseList as a shallow pointer within the object
        //
        m_pl = pl;

        XML_Node* kinXMLPhase = 0;
        ThermoPhase* kinPhase = 0;

        if (iskin >= 0) {
            kinXMLPhase = pl->surPhaseXMLNode(iskin);
            kinPhase = &(pl->surPhase(iskin));
        } else {
            throw CanteraError("importFromPL", "iskin is neg");
        }

        int nPhasesFound = pl->nSurPhases() + pl->nVolPhases();
        /*
         * Resize the internal list of pointers and
         * get a pointer to the vacant ThermoPhase pointer
         */
        std::vector<ThermoPhase*> tpList;
        tpList.clear();
        tpList_IDs_.clear();
        tplRead.resize(nPhasesFound, 0);
        kinOrder.resize(nPhasesFound, -1);
        xmlList.clear();
        iphaseKin_ = -1;
        if (iskin >= 0) {
            iphaseKin_ = iskin + pl->nVolPhases();
            m_DoSurfKinetics = true;
        } 

        numPhases_ = 0;
        if (iphaseKin_ >= 0) {
            xmlList.push_back(kinXMLPhase);
            tpList.push_back(kinPhase);
            tpList_IDs_.push_back(kinPhase->id());
            tplRead[numPhases_] = 1;
            numPhases_++;
        }

        /*
         *  OK, we have settled in on the kinetics object that we will process.
         *  Now, go look at the phaseArray XML field to get a listing of the ThermoPhases
         *  involved with the kinetics object.
         */
        XML_Node* phaseArrayXML = 0;
        if (iphaseKin_ >= 0) {
            XML_Node* xmlPhase = kinXMLPhase;
            phaseArrayXML = xmlPhase->findNameID("phaseArray", "");
            if (phaseArrayXML) {
                vector<string> phase_ids;
                ctml::getStringArray(*phaseArrayXML, phase_ids);
                int npToFind = phase_ids.size();
                for (iph = 0; iph < npToFind; iph++) {
                    string phaseID = phase_ids[iph];
                    bool found = false;
                    for (int jph = 0; jph < pl->nSurPhases(); jph++) {
                        XML_Node* xmlPhase_j = pl->surPhaseXMLNode(jph);
                        string pname = xmlPhase_j->operator[]("id");
                        if (phaseID == pname) {
                            found = true;
                            xmlList.push_back(xmlPhase_j);
                            tpList.push_back(&(pl->surPhase(jph)));
                            tpList_IDs_.push_back(pl->surPhase(jph).id());
                            tplRead[jph] = 1;
                            numPhases_++;
                            break;
                        }
                    }
                    if (!found) {
                        for (int jph = 0; jph < pl->nVolPhases(); jph++) {
                            XML_Node* xmlPhase_j = pl->volPhaseXMLNode(jph);
                            string pname = xmlPhase_j->operator[]("id");
                            if (phaseID == pname) {
                                found = true;
                                xmlList.push_back(xmlPhase_j);
                                tpList.push_back(&(pl->volPhase(jph)));
                                tpList_IDs_.push_back(pl->volPhase(jph).id());
                                tplRead[jph] = 1;
                                numPhases_++;
                                break;
                            }
                        }
                    }

                    if (!found) {
                        throw CanteraError("ReactingSurDomain::importFromPL()",
                                           "Phase, requested in phaseArray, was not found: "
                                           + phaseID);
                    }
                }
            }
        } else {
            for (iph = 0; iph < pl->nSurPhases(); iph++) {
                xmlList.push_back(pl->surPhaseXMLNode(iph));
                tpList.push_back(&(pl->surPhase(iph)));
                tpList_IDs_.push_back(pl->surPhase(iph).id());
                tplRead[numPhases_] = 1;
                numPhases_++;
            }
            for (iph = 0; iph < pl->nVolPhases(); iph++) {
                xmlList.push_back(pl->volPhaseXMLNode(iph));
                tpList.push_back(&(pl->volPhase(iph)));
                tpList_IDs_.push_back(pl->volPhase(iph).id());
                tplRead[numPhases_] = 1;
                numPhases_++;
            }
        }

        /*
         * Fill in the kinetics object k, by querying the
         * const XML_Node tree located at xmlPhase. The source terms and
         * eventually the source term vector will be constructed
         * from the list of ThermoPhases in the vector, tpList
         */
        XML_Node* xmlPhase = pl->surPhaseXMLNode(iskin);
        bool ok = importKinetics(*xmlPhase, tpList, this);
        if (!ok) {
            throw CanteraError("ReactingSurDomain::importFromPL()", "importKinetics() returned an error");
        }

        /*
         *  Create a mapping between the ReactingSurfPhase to the PhaseList phase
         */
        int nKinPhases = nPhases();
        kinOrder.resize(nKinPhases, -1);
        PLtoKinPhaseIndex_.resize(pl->nPhases(), -1);
        PLtoKinSpeciesIndex_.resize(pl->nSpecies(), -1);
        for (int kph = 0; kph < nKinPhases; kph++) {
            ThermoPhase& tt = thermo(kph);
            string kname = tt.id();
            int jph = -1;
            for (int iph = 0; iph < pl->nPhases(); iph++) {
                ThermoPhase& pp = pl->thermo(iph);
                std::string iname = pp.id();
                if (iname == kname) {
                    jph = iph;

                    break;
                }
            }
            if (jph == -1) {
                throw CanteraError("ReactingSurDomain::importFromPL()", "phase not found");
            }
            kinOrder[kph] = jph;
            PLtoKinPhaseIndex_[jph] = kph;

            int PLkstart = pl->getGlobalSpeciesIndex(jph, 0);
            int nspPhase = tt.nSpecies();
            for (int k = 0; k < nspPhase; k++) {
                if (PLtoKinSpeciesIndex_[k + PLkstart] != -1) {
                    throw CanteraError("ReactingSurDomain::importFromPL()",
                                       "Indexing error found while initializing  PLtoKinSpeciesIndex_");
                }
                PLtoKinSpeciesIndex_[k + PLkstart] = m_start[kph] + k;
            }
        }

        /*
         * Resize the arrays based on species number
         */
        speciesProductionRates_.resize(m_kk, 0.0);
        speciesCreationRates_.resize(m_kk, 0.0);
        speciesDestructionRates_.resize(m_kk, 0.0);

        /*
         * Resize the arrays based on the number of reactions
         */ 
        rmcVector.resize(m_ii, 0);
        for (size_t i = 0; i < m_ii; i++) {
            rmcVector[i] = new RxnMolChange(this, i);
        }
        deltaGRxn_.resize(m_ii, 0.0);

        //
        //  Identify the electron phase
        //
        identifyMetalPhase();

        return true;

    } catch (CanteraError) {
        showErrors(cout);
        throw CanteraError("ReactingSurDomain::importFromPL()",
                           "error encountered");
        return false;
    }
}

//====================================================================================================================
// An an override for the OCV
void ReactingSurDomain::addOCVoverride(OCV_Override_input *ocv_ptr)
{
    //
    // Save the pointer for the input information  (Question, should I make a deep copy?)
    //
    ocv_ptr_ = ocv_ptr;



}
//====================================================================================================================
void  ReactingSurDomain::deriveEffectiveChemPot()
{
    //
    //  Calculate the delta G and the open circuit voltage normally
    //
    //
    //
    getDeltaSSGibbs(DATA_PTR(deltaGRxn_));
  
    double phiRxnOrig = 0.0;

    //
    //  If we don't have an electrode reaction bail as being confused
    //
    if (metalPhaseRS_ < 0) {
	throw CanteraError("", "shouldn't be here");
    }
    //
    //  If we don't have an open circuit potential override situation, bail
    if (!ocv_ptr_) {
	return;
    }

    //
    //  Figure out the reaction id for the override
    //
    int rxnID = ocv_ptr_->rxnID;
    //
    //  Get a pointer to the RxnMolChange struct, which contains more info about the reaction
    //
    RxnMolChange* rmc = rmcVector[rxnID];
    //
    //   Find the number of stoichiometric electrons in the reaction
    //
    double nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseRS_];
    //
    //  If the number of stoichiometric electrons is zero, we are in kind of a bind, as we shouldn't
    //  be specifying an open circuit voltage in the first place!
    //
   if (nStoichElectrons == 0.0) {
	throw CanteraError("", "shouldn't be here");
    }
    //
    //  Calculate the open circuit voltage from this relation that would occur if it was not being overwritten
    //
    phiRxnOrig = deltaGRxn_[rxnID] / Faraday / nStoichElectrons;
    
    //
    // In order to calculate the OCV, we need the relative extent of reaction value
    //


    //
    //  Now calculate the OCV value to be used from the fit
    //

 
  
    
    
}










//====================================================================================================================
} // End of namespace cantera
//======================================================================================================================

