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
//#include "Cantera/kernel/importKinetics.h"

//static int DebugPrinting = 1;


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
    m_transport(0),
    numPhases(0),
    xmlList(0),
    kinOrder(0),
    PLtoKinPhaseIndex_(0),
    PLtoKinSpeciesIndex_(0),
    iphaseKin(-1),
    tplRead(0),
    m_DoSurfKinetics(false),
    speciesProductionRates_(0),
    speciesCreationRates_(0),
    speciesDestructionRates_(0),
    m_pl(0),
    metalPhaseRS_(-1),
    kElectronRS_(-1),
    solnPhaseRS_(-1)
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
    m_transport(0),
    numPhases(0),
    xmlList(0),
    kinOrder(0),
    PLtoKinPhaseIndex_(0),
    PLtoKinSpeciesIndex_(0),
    iphaseKin(-1),
    tplRead(0),
    m_DoSurfKinetics(false),
    speciesProductionRates_(0),
    speciesCreationRates_(0),
    speciesDestructionRates_(0),
    m_pl(0),
    metalPhaseRS_(-1),
    kElectronRS_(-1),
    solnPhaseRS_(-1)
{
    /*
     * Call the assignment operator
     */
    *this = operator=(right);
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

    InterfaceKinetics::operator=(right);

    numPhases       = right.numPhases;
    // Shallow copy of xmlList pointers -> beware
    xmlList         = right.xmlList;
    kinOrder        = right.kinOrder;
    PLtoKinPhaseIndex_ = right.PLtoKinPhaseIndex_;
    PLtoKinSpeciesIndex_ = right.PLtoKinSpeciesIndex_;
    iphaseKin       = right.iphaseKin;
    tplRead         = right.tplRead;
    m_DoSurfKinetics = right.m_DoSurfKinetics;
    speciesProductionRates_ = right.speciesProductionRates_;
    speciesCreationRates_ = right.speciesCreationRates_;
    speciesDestructionRates_ = right.speciesDestructionRates_;

    // Shallow copy of m_pl -> beware
    m_pl            = m_pl;

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
 * @return Vector of length m_kk containing the species net
 *         production rates
 */
const std::vector<double>& ReactingSurDomain::calcNetProductionRates()
{
    getNetProductionRates(&speciesProductionRates_[0]);
    return speciesProductionRates_;
}
//====================================================================================================================
//    Returns a reference to the calculated creation rates of species
/*
 *   This routine calls thet getCreationRate function
 *   and then returns a reference to the result.
 *
 * @return Vector of length m_kk containing the species creation rates
 */
const std::vector<double>& ReactingSurDomain::calcCreationRates()
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
 * @return Vector of length m_kk containing the species destruction rates
 */
const std::vector<double>& ReactingSurDomain::calcDestructionRates()
{
    getDestructionRates(&speciesDestructionRates_[0]);
    return speciesDestructionRates_;
}

//====================================================================================================================

void ReactingSurDomain::getExchangeCurrentFormulation(int irxn,
        double* nStoich, doublereal* OCV, doublereal* io,
        doublereal* nu)
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
    // we have a vector of standard concentrations calculated from the routine below
    //            m_StandardConc[ik]
    getExchangeCurrentQuantities();



    const vector_fp& rf = m_kdata->m_rfn;
    const vector_fp& rkc= m_kdata->m_rkcn;


    // start with the forward reaction rate
    double iO = rf[irxn] * Faraday * nStoichElectrons;

    if (m_beta[irxn] > 0.0) {
        iO *= pow(rkc[irxn], m_beta[irxn]);
    }
    double b = m_beta[irxn];
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
    *nu = E - *OCV;
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
    for (int i = 0; i < mix.numPhases; i++) {
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
 */
bool ReactingSurDomain::
importFromPL(Cantera::PhaseList* pl, int ivkin, int iskin)
{
    try {
        int iph;
        m_pl = pl;

        XML_Node* kinXMLPhase = 0;
        ThermoPhase* kinPhase = 0;

        //XML_Node *vPhase = pl->SurPhaseXMLNodes[0];
        if (iskin >= 0) {
            kinXMLPhase = pl->surPhaseXMLNode(iskin);
            kinPhase = &(pl->surPhase(iskin));
        } else    if (ivkin >= 0) {
            kinXMLPhase = pl->surPhaseXMLNode(ivkin);
            kinPhase = &(pl->surPhase(ivkin));
        }
        //AssertThrow(vPhase, "vPhase must be defined");



        int nPhasesFound = pl->nSurPhases() + pl->nVolPhases();
        /*
         * Resize the internal list of pointers and
         * get a pointer to the vacant ThermoPhase pointer
         */
        std::vector<ThermoPhase*> tpList;
        tpList.clear();
        tplRead.resize(nPhasesFound, 0);
        kinOrder.resize(nPhasesFound, -1);
        xmlList.clear();
        iphaseKin = -1;
        if (iskin >= 0) {
            iphaseKin = iskin + pl->nVolPhases();
            m_DoSurfKinetics = true;
        } else if (ivkin >= 0) {
            throw CanteraError("", "err");
        }


        numPhases = 0;
        if (iphaseKin >= 0) {
            xmlList.push_back(kinXMLPhase);
            tpList.push_back(kinPhase);
            tplRead[numPhases] = 1;
            numPhases++;
        }

        /*
         *  OK, we have settled in on the kinetics object that we will process.
         *  Now, go look at the phaseArray XML field to get a listing of the ThermoPhases
         *  involved with the kinetics object.
         */
        XML_Node* phaseArrayXML = 0;
        if (iphaseKin >= 0) {
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
                            tplRead[jph] = 1;
                            numPhases++;
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
                                tplRead[jph] = 1;
                                numPhases++;
                                break;
                            }
                        }
                    }

                    if (!found) {
                        throw CanteraError("ReactingSurDomain::importFromPL",
                                           "Phase, requested in phaseArray, was not found: "
                                           + phaseID);
                    }
                }
            }
        } else {
            for (iph = 0; iph < pl->nSurPhases(); iph++) {
                xmlList.push_back(pl->surPhaseXMLNode(iph));
                tpList.push_back(&(pl->surPhase(iph)));
                tplRead[numPhases] = 1;
                numPhases++;
            }
            for (iph = 0; iph < pl->nVolPhases(); iph++) {
                xmlList.push_back(pl->volPhaseXMLNode(iph));
                tpList.push_back(&(pl->volPhase(iph)));
                tplRead[numPhases] = 1;
                numPhases++;
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
            throw CanteraError("", "err");
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
                string iname = pp.id();
                if (iname == kname) {
                    jph = iph;

                    break;
                }
            }
            if (jph == -1) {
                throw CanteraError("importFromPL", "not found");
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

        int nr = nReactions();
        rmcVector.resize(nr,0);
        for (int i = 0; i < nr; i++) {
            rmcVector[i] = new RxnMolChange(this, i);
        }

        identifyMetalPhase();


        return true;

    } catch (CanteraError) {
        showErrors(cout);
        throw CanteraError("ReactingSurDomain::importFromPL",
                           "error encountered");
        return false;
    }


}
//====================================================================================================================
} // End of namespace cantera
//======================================================================================================================

