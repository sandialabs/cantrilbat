/**
 * @file ReactingVolDomain.cpp
 *
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */


#include "ReactingSurDomain.h"

#include "cantera/multiphase/PhaseList.h"
#include "cantera/kinetics/RxnMolChange.h"
#include "cantera/kinetics.h"
#include "cantera/base/ctml.h"


#ifdef useZuzaxNamespace
#define ZZctml ztml
#else
#define ZZctml ctml
#endif

using namespace std;
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif 
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
    iphaseKin(npos),
    tplRead(0),
    m_DoSurfKinetics(false),
    speciesProductionRates_(0),
    speciesCreationRates_(0),
    speciesDestructionRates_(0),
    m_pl(0),
    m_ok(false), 
    m_XMLPhaseTree(0)
  {
  }
  //====================================================================================================================
  // Copy Constructor for the %Kinetics object.
  /*
   * Currently, this is not fully implemented. If called it will
   * throw an exception.
   */
   ReactingSurDomain::ReactingSurDomain(const ReactingSurDomain &right) :
    InterfaceKinetics(),
    m_transport(0),
    numPhases(0),
    xmlList(0),
    kinOrder(0),
    PLtoKinPhaseIndex_(0),
    PLtoKinSpeciesIndex_(0),
    iphaseKin(npos),
    tplRead(0),
    m_DoSurfKinetics(false),
    speciesProductionRates_(0),
    speciesCreationRates_(0),
    speciesDestructionRates_(0),
    m_pl(0),
    m_ok(false), 
    m_XMLPhaseTree(0)
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
  ReactingSurDomain&  ReactingSurDomain::operator=(const ReactingSurDomain &right)
  {
    /*
     * Check for self assignment.
     */
    if (this == &right) return *this;

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
    m_ok            = right.m_ok;
    // Shallow copy of m_XMLPhaseTree -> beware
    m_XMLPhaseTree  = m_XMLPhaseTree;
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
  ReactingSurDomain::~ReactingSurDomain() {
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
  Kinetics* ReactingSurDomain::duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const
  {
    ReactingSurDomain* rsd = new ReactingSurDomain(*this); 
    rsd->assignShallowPointers(tpVector);
    return dynamic_cast<Kinetics *>(rsd);  
  }
  //====================================================================================================================
  // Returns a reference to the calculated production rates of species
  /*
   *   This routine calls thet getNetProductionRate function
   *   and then returns a reference to the result.
   *
   * @return Vector of length m_NumKinSpecies containing the species net 
   *         production rates
   */
  const std::vector<double> & ReactingSurDomain::calcNetProductionRates() {
    getNetProductionRates(&speciesProductionRates_[0]);
    return speciesProductionRates_;
  }
  //====================================================================================================================
  //    Returns a reference to the calculated creation rates of species
  /*
   *   This routine calls thet getCreationRate function
   *   and then returns a reference to the result.
   *
   * @return Vector of length m_NumKinSpecies containing the species creation rates
   */
  const std::vector<double> & ReactingSurDomain::calcCreationRates() {
    getCreationRates(&speciesCreationRates_[0]);
    return speciesCreationRates_;
  }
  //====================================================================================================================
  // Returns a reference to the calculated destruction rates of species
  /*
   *   This routine calls thet getDestructionRate function
   *   and then returns a reference to the result.
   *
   * @return Vector of length m_NumKinSpecies containing the species destruction rates
   */
  const std::vector<double> & ReactingSurDomain::calcDestructionRates() {
    getDestructionRates(&speciesDestructionRates_[0]);
    return speciesDestructionRates_;
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
  std::ostream& operator<<(std::ostream& s, ReactingSurDomain& mix)
  {
    ThermoPhase *th;
    InterfaceKinetics *iK = &mix;
    for (size_t i = 0; i < mix.numPhases; i++) {
      th = &(iK->thermo(i));
      std::string r = th->report(true);
      s << r;
    }
    return s;
  }
  //====================================================================================================================
  /*
   */
  bool ReactingSurDomain::importFromPL(PhaseList *pl, int ivkin, int iskin)  
  {
    try {
      m_pl = pl;
      
      XML_Node *kinXMLPhase = 0;
      ThermoPhase *kinPhase = 0;

      if (iskin >= 0) {
	kinXMLPhase =&( pl->surPhaseXMLNode(iskin) );
	kinPhase = &(pl->surPhase(iskin));
      } else    if (ivkin >= 0) {
	kinXMLPhase = &( pl->volPhaseXMLNode(ivkin) );
	kinPhase = &(pl->volPhase(ivkin));
      }
      //AssertThrow(vPhase, "vPhase must be defined");



      int nPhasesFound = pl->nSurPhases() + pl->nVolPhases();
      /*
       * Resize the internal list of pointers and
       * get a pointer to the vacant ThermoPhase pointer
       */
      std::vector<ThermoPhase *> tpList;
      tpList.clear();
      tplRead.resize(nPhasesFound, 0);
      kinOrder.resize(nPhasesFound, npos);
      xmlList.clear();
      iphaseKin = npos; 
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
      XML_Node *phaseArrayXML = 0;
      if (iphaseKin >= 0) {
	XML_Node *xmlPhase = kinXMLPhase;
	phaseArrayXML = xmlPhase->findNameID("phaseArray", "");
	if (phaseArrayXML) {
	  std::vector<std::string> phase_ids;
	  ZZctml::getTokenArray(*phaseArrayXML, phase_ids);
	  size_t npToFind = phase_ids.size();
	  for (size_t iph = 0; iph < npToFind; iph++) {
	    std::string phaseID = phase_ids[iph];
	    bool found = false;
	    for (size_t jph = 0; jph < pl->nSurPhases(); jph++) {
	      XML_Node *xmlPhase_j = &( pl->surPhaseXMLNode(jph) );
	      std::string pname = xmlPhase_j->operator[]("id");
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
	      for (size_t jph = 0; jph < pl->nVolPhases(); jph++) {
		XML_Node *xmlPhase_j = &( pl->volPhaseXMLNode(jph) );
		std::string pname = xmlPhase_j->operator[]("id");
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
	for (size_t iph = 0; iph < pl->nSurPhases(); iph++) {
	  xmlList.push_back(&( pl->surPhaseXMLNode(iph) ));
	  tpList.push_back(&(pl->surPhase(iph)));
	  tplRead[numPhases] = 1;
	  numPhases++;
	}
	for (size_t iph = 0; iph < pl->nVolPhases(); iph++) {
	  xmlList.push_back(&( pl->volPhaseXMLNode(iph) ));
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
      XML_Node *xmlPhase = &( pl->surPhaseXMLNode(iskin) );
      bool ok = importKinetics(*xmlPhase, tpList, this);
      if (!ok) {
	throw CanteraError("", "err");
      }
 
      /*
       *  Create a mapping between the ReactingSurfPhase to the PhaseList phase
       */
      size_t nKinPhases = nPhases();
      kinOrder.resize(nKinPhases, npos);
      PLtoKinPhaseIndex_.resize(pl->nPhases(), npos);
      PLtoKinSpeciesIndex_.resize(pl->nSpecies(), npos);
      for (size_t kph = 0; kph < nKinPhases; kph++) {
	ThermoPhase &tt = thermo(kph);
	string kname = tt.id();
	size_t jph = npos;
	for (size_t iph = 0; iph < pl->nPhases(); iph++) {
	  ThermoPhase &pp = pl->thermo(iph);
	  string iname = pp.id();
	  if (iname == kname) {
	    jph = iph;

	    break;
	  }
	}
	if (jph == npos) {
	  throw CanteraError("importFromPL", "not found");
	}
	kinOrder[kph] = jph;
	PLtoKinPhaseIndex_[jph] = kph;

        size_t PLkstart = pl->globalSpeciesIndex(jph, 0);
        size_t nspPhase = tt.nSpecies(); 
	for (size_t k = 0; k < nspPhase; k++) {
          if (PLtoKinSpeciesIndex_[k + PLkstart] != npos) {
             throw CanteraError("ReactingSurDomain::importFromPL()",
                                "Indexing error found while initializing  PLtoKinSpeciesIndex_");
          }
          PLtoKinSpeciesIndex_[k + PLkstart] = m_start[kph] + k;
        }
      }
    
      /*
       * Resize the arrays based on species number
       */
      speciesProductionRates_.resize(m_NumKinSpecies, 0.0);
      speciesCreationRates_.resize(m_NumKinSpecies, 0.0);
      speciesDestructionRates_.resize(m_NumKinSpecies, 0.0);

      size_t nr = nReactions();
      rmcVector.resize(nr,0);
      for (size_t i = 0; i < nr; i++) {
	rmcVector[i] = new RxnMolChange(this, i);
      }
      return true;

    }
    catch (CanteraError) {
      showErrors(cout);
      throw CanteraError("ReactingSurDomain::importFromXML", "error encountered");
      return false;
    }
  }
  //====================================================================================================================
} // End of namespace cantera
//======================================================================================================================

