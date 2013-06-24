/**
 *  @file ElectrodeKinetics.cpp 
 *
 */
/*
 *  $Id: ElectrodeKinetics.cpp 508 2013-01-07 22:54:04Z hkmoffa $
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */



#include "ElectrodeKinetics.h"
#include "cantera/thermo/SurfPhase.h"

#include "cantera/kinetics/ReactionData.h"
#include "ReactionDataElectrode.h"
#include "cantera/kinetics/RateCoeffMgr.h"

#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/importKinetics.h"

#include <memory>

using namespace std;

namespace Cantera {

  //! these are all used to check for duplicate reactions
  class rxninfoEK {
  public:
    //! rdata
    std::vector< std::map<int, doublereal> > m_rdata;

    //! string name
    std::vector<std::string> m_eqn;

    //! string vector of ints
    std::vector<int> m_dup;

    //! string vector of ints
    std::vector<int>  m_nr;

    //! string vector of ints
    std::vector<int>  m_typ;

    //! vector of bools.
    std::vector<bool> m_rev;

    ~rxninfoEK() {
      m_eqn.clear();
      m_dup.clear();
      m_nr.clear();
      m_typ.clear();
      m_rdata.clear();
    }

    /**
   *  Install an individual reaction into a kinetics manager. The
   *  data for the reaction is in the xml_node r. In other words, r
   *  points directly to a ctml element named "reaction". i refers
   *  to the number id of the reaction in the kinetics object.
   * 
   * @param i Reaction number.
   * @param r XML_Node containing reaction data.
   * @param k Kinetics manager to which reaction will be added.
   * @param default_phase Default phase for locating a species
   * @param rule Rule for handling reactions with missing species 
   *             (skip or flag as error)
   * @param validate_rxn If true, check that this reaction is not a 
   *                     duplicate of one already entered, and check that the reaction 
   *                     balances.
   *
   * @ingroup kineticsmgr
   */
    bool installReaction(int i, const XML_Node& r, ElectrodeKinetics* k, 
			 std::string default_phase, const ReactionRules& rule, bool validate_rxn) ;
  };

  // Constructor
  /*
   * Took out the option to initialize with a ThermoPhase Object.
   */
  ElectrodeKinetics::ElectrodeKinetics() :
    InterfaceKinetics()
  {
  }
  
  // Destructor
  ElectrodeKinetics::~ElectrodeKinetics() {
  }
  


  /**
   * For reactions that transfer charge across a potential difference,
   * the activation energies are modified by the potential difference.
   * (see, for example, ...). This method applies this correction.
   */
  void ElectrodeKinetics::applyButlerVolmerCorrection(doublereal* kf) {
    size_t i;

    int n, nsp, k, ik=0;
    doublereal rt = GasConstant*thermo(0).temperature();
    doublereal rrt = 1.0/rt;
    int np = nPhases();

    // compute the electrical potential energy of each species
    for (n = 0; n < np; n++) {
      nsp = thermo(n).nSpecies();
      for (k = 0; k < nsp; k++) {
	m_pot[ik] = Faraday*thermo(n).charge(k)*m_phi[n];
	ik++;
      }
    }

    // compute the change in electrical potential energy for each
    // reaction. This will only be non-zero if a potential
    // difference is present.
    //fill(m_rwork.begin(), m_rwork.begin() + m_ii, 0.0);
    //m_reactantStoich.decrementReactions(m_pot.begin(), m_rwork.begin()); 
    //m_revProductStoich.incrementReactions(m_pot.begin(), m_rwork.begin());
    //m_irrevProductStoich.incrementReactions(m_pot.begin(), m_rwork.begin());
    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_pot), 
				 DATA_PTR(m_rwork));

    // modify the reaction rates. Only modify those with a
    // non-zero activation energy, and do not decrease the
    // activation energy below zero.
    doublereal ea, eamod;

    for (i = 0; i < m_ii; i++) {
      eamod = 0.5*m_rwork[i];
      if (eamod != 0.0 && m_E[i] != 0.0) {
	ea = GasConstant * m_E[i];
	if (eamod + ea < 0.0) {
	  eamod = -ea;
	  writelog("warning: modified E < 0.\n");
	}
	kf[i] *= exp(-eamod*rrt);


      }
    }
  }


  // Import a reaction mechanism for a phase or an interface. 
  // This is modified for ElectrodeKinetics only pending a generalized treatment.
  /*
   * Responsibilities:
   *     -    Figures out if check_for_duplicates should be attempted
   *     -    Processes the phaseArray XML element, and loads the
   *          Kinetics object with ThermoPhase pointers.
   *     -    Figures out the kinetics species vector 
   *     -    Calls the Kinetics virtual function init()
   *     -    Call  the installReactionArrays function
   *
   * @param phase This is an xml node containing a description
   *              of a phase. Within the phase is a XML element
   *              called reactionArray containing the location
   *              of the description of the reactions that make
   *              up the kinetics object. 
   *              Also within the phase is an XML element called
   *              phaseArray containing a listing of other phases
   *              that participate in the kinetics mechanism.
   *
   * @param th    This is a list of ThermoPhase pointers which must
   *              include all of
   *              the phases that participate in the kinetics
   *              operator. All of the phases must have already
   *              been initialized and formed within Cantera.
   *              However, their pointers should not have been
   *              added to the Kinetics object; this addition
   *              is carried out here. Additional phases may
   *              be include; these have no effect.
   *
   * @return
   *         Return true if successful, false otherwise.
   */
  bool ElectrodeKinetics::importKinetics(const XML_Node& phase,
					 std::vector<ThermoPhase*> th) {



    // This phase will be the owning phase for the kinetics operator
    // For interfaces, it is the surface phase between two volumes.
    // For homogeneous kinetics, it's the current volumetric phase.
    string owning_phase_name = phase["id"];

    bool check_for_duplicates = false;
    if (phase.parent()->hasChild("validate")) {
      const XML_Node& d = phase.parent()->child("validate");
      if (d["reactions"] == "yes") check_for_duplicates = true;
    }

    // if other phases are involved in the reaction mechanism,
    // they must be listed in a 'phaseArray' child
    // element. Homogeneous mechanisms do not need to include a
    // phaseArray element.

    vector<string> phase_ids;
    if (phase.hasChild("phaseArray")) {
      const XML_Node& pa = phase.child("phaseArray");
      ctml::getStringArray(pa, phase_ids);
    }
    phase_ids.push_back(owning_phase_name);
            
    int np = static_cast<int>(phase_ids.size());
    int nt = static_cast<int>(th.size());

    // for each referenced phase, attempt to find its id among those
    // phases specified. 
    bool phase_ok;

    string phase_id;
    string msg = "";
    for (int n = 0; n < np; n++) {
      phase_id = phase_ids[n];
      phase_ok = false;

      // loop over the supplied 'ThermoPhase' objects representing
      // phases, to find an object with the same id.
      for (int m = 0; m < nt; m++) {
	if (th[m]->id() == phase_id) {
	  phase_ok = true;

          // if no phase with this id has been added to
          // the kinetics manager yet, then add this one
          if (phaseIndex(phase_id) == npos) {
	    addPhase(*th[m]);
	  }
	}
	msg += " "+th[m]->id();
      }
      if (!phase_ok) {
	throw CanteraError("ElectrodeKinetics::importKinetics",
			   "phase "+phase_id+" not found. Supplied phases are:"+msg);
      }
    }

    // Allocates arrays, etc. Must be called after the phases have 
    // been added to 'kin', so that the number of species in each
    // phase is known.
    init();

    // Install the reactions.
    return installReactionArrays(phase, owning_phase_name, check_for_duplicates);
  }

  /*
   * Add a phase to the kinetics manager object. This must
   * be done before the function init() is called or
   * before any reactions are input.
   * The following fields are updated:
   *  m_start -> vector of integers, containing the
   *             starting position of the species for
   *             each phase in the kinetics mechanism.
   *  m_surfphase -> index of the surface phase.
   *  m_thermo -> vector of pointers to ThermoPhase phases
   *              that participate in the kinetics
   *              mechanism.
   *  m_phaseindex -> map containing the string id of each
   *              ThermoPhase phase as a key and the
   *              index of the phase within the kinetics
   *              manager object as the value.
   */
  void ElectrodeKinetics::addPhase(thermo_t& thermo) {

    // if not the first thermo object, set the start position
    // to that of the last object added + the number of its species
    if (m_thermo.size() > 0) {
      m_start.push_back(m_start.back()
			+ m_thermo.back()->nSpecies());
    }
    // otherwise start at 0
    else {
      m_start.push_back(0);
    }

    // the phase with lowest dimensionality is assumed to be the
    // phase/interface at which reactions take place
    if (thermo.nDim() <= m_mindim) {
      m_mindim = thermo.nDim();
      m_rxnphase = nPhases();
    }


    // there should only be one surface phase
    int ptype = -100;
    if (type() == cEdgeKinetics) {
        ptype = cEdge;
    } else if (type() == cInterfaceKinetics) {
        ptype = cSurf;
    } else if (type() == cElectrodeKinetics) {
        ptype = cSurf;
    }
    if (thermo.eosType() == ptype) {
        m_surfphase = nPhases();
        m_rxnphase = nPhases();
    }

    m_thermo.push_back(&thermo);
    m_phaseindex[m_thermo.back()->id()] = nPhases();
    //
    m_phaseExists.push_back(true);
  }

  // Read and Install Reactions into Kinetics Object
  /*
   *  Take information from the XML tree, p, about reactions
   *  and install them into the kinetics object, kin. 
   *  default_phase is the default phase to assume when
   *  looking up species.
   *
   *  At this point, p usually refers to the phase xml element.
   *  One of the children of this element is reactionArray,
   *  the element which determines where in the xml file to
   *  look up the reaction rate data pertaining to the phase.
   *
   *  On return, if reaction instantiation goes correctly, return true.
   *  If there is a problem, return false.
   */
  bool ElectrodeKinetics::installReactionArrays(const XML_Node& p,  
						std::string default_phase, 
						bool check_for_duplicates) {

    const std::auto_ptr< rxninfoEK > _rxns( new rxninfoEK );
  

    vector<XML_Node*> rarrays;
    int itot = 0;
    /*
     * Search the children of the phase element for the
     * xml element named reactionArray. If we can't find it,
     * then return signaling having not found any reactions.
     * Apparently, we allow multiple reactionArray elements here
     * Each one will be processed sequentially, with the
     * end result being purely additive.
     */
    p.getChildren("reactionArray",rarrays);
    int na = static_cast<int>(rarrays.size());
    if (na == 0) return false;
    for (int n = 0; n < na; n++) {
      /*
       * Go get a reference to the current xml element, 
       * reactionArray. We will process this element now.
       */
      const XML_Node& rxns = *rarrays[n];
      /*
       * The reactionArray element has an attribute called,
       * datasrc. The value of the attribute is the xml
       * element comprising the top of the
       * tree of reactions for the phase.
       * Find this datasrc element starting with the root
       * of the current xml node.
       */
      const XML_Node* rdata = get_XML_Node(rxns["datasrc"], &rxns.root());
      /*
       * If the reactionArray element has a child element named
       * "skip", and if the attribute of skip called "species" has
       * a value of "undeclared", we will set rxnrule = 1.
       * rxnrule is passed to the routine that parses each individual
       * reaction. I believe what this means is that the parser will
       * skip all reactions containing an undefined species without
       * throwing an error condition.
       */
      ReactionRules rxnrule;
      if (rxns.hasChild("skip")) {
	const XML_Node& sk = rxns.child("skip");
	string sskip = sk["species"];
	if (sskip == "undeclared") {
	  rxnrule.skipUndeclaredSpecies = true;
	}
      }
      int i, nrxns = 0;
      /*
       * Search for child elements called include. We only include
       * a reaction if it's tagged by one of the include fields.
       * Or, we include all reactions if there are no include fields.
       */
      vector<XML_Node*> incl;
      rxns.getChildren("include",incl);
      int ninc = static_cast<int>(incl.size());

      vector<XML_Node*> allrxns;
      rdata->getChildren("reaction",allrxns);
      nrxns = static_cast<int>(allrxns.size());
      // if no 'include' directive, then include all reactions
      if (ninc == 0) {
	for (i = 0; i < nrxns; i++) {
	  const XML_Node* r = allrxns[i];
	  if (r) {
	    if (_rxns->installReaction(itot, *r, this, 
				default_phase, rxnrule, check_for_duplicates)) ++itot;
	  }
	}
      }
      else {
	for (int nii = 0; nii < ninc; nii++) {
	  const XML_Node& ii = *incl[nii];
	  string imin = ii["min"];
	  string imax = ii["max"];

	  string::size_type iwild = string::npos;
	  if (imax == imin) {
	    iwild = imin.find("*");
	    if (iwild != string::npos) {
	      imin = imin.substr(0,iwild);
	      imax = imin;
	    }
	  }

	  for (i = 0; i < nrxns; i++) {
	    const XML_Node* r = allrxns[i];
	    string rxid;
	    if (r) {
	      rxid = (*r)["id"];
	      if (iwild != string::npos) {
		rxid = rxid.substr(0,iwild);
	      }
	      /*
	       * To decide whether the reaction is included or not
	       * we do a lexical min max and operation. This 
	       * sometimes has surprising results.
	       */
	      if ((rxid >= imin) && (rxid <= imax)) {
		if (_rxns->installReaction(itot, *r, this, 
				    default_phase, rxnrule, check_for_duplicates)) ++itot;
	      }
	    }
	  }
	}
      }
    }

    /*
     * Finalize the installation of the kinetics, now that we know
     * the true number of reactions in the mechanism, itot.
     */
    finalize();
  
    return true;
  }

  // Add a single reaction to the mechanism. 
  /*
   *   This routine
   *   must be called after init() and before finalize().
   *   This function branches on the types of reactions allowed
   *   by the interfaceKinetics manager in order to install
   *   the reaction correctly in the manager.
   *   The manager allows the following reaction types
   *         - Elementary
   *         - Surface
   *         - Global  
   *   There is no difference between elementary and surface 
   *   reactions.
   *
   * @param r Reference to the ReactionData data type
   */
  void ElectrodeKinetics::
    //Compile error if r is const, TODO: check if ::addReaction(r) should be able to take const
  addReaction(ReactionData& r) {
    InterfaceKinetics::addReaction(r);
  }


  void ElectrodeKinetics::
  addElementaryReaction(ReactionData& r) {
    InterfaceKinetics::addElementaryReaction(r);
  }

        
  void ElectrodeKinetics::installReagents(const ReactionData& r) {

    int n, ns, m; 
    doublereal nsFlt;
    /*
     * extend temporary storage by one for this rxn.
     */
    m_kdata->m_ropf.push_back(0.0);
    m_kdata->m_ropr.push_back(0.0);
    m_kdata->m_ropnet.push_back(0.0);
    m_kdata->m_rkcn.push_back(0.0);

    /*
     * Obtain the current reaction index for the reaction that we
     * are adding. The first reaction is labeled 0.
     */
    int rnum = nReactions();

    // vectors rk and pk are lists of species numbers, with
    // repeated entries for species with stoichiometric
    // coefficients > 1. This allows the reaction to be defined
    // with unity reaction order for each reactant, and so the
    // faster method 'multiply' can be used to compute the rate of
    // progress instead of 'power'.

    std::vector<size_t> rk;
    int nr = r.reactants.size();
    for (n = 0; n < nr; n++) {
      nsFlt = r.rstoich[n];
      ns = (int) nsFlt;
      if ((doublereal) ns != nsFlt) {
	if (ns < 1) ns = 1;
      }
      /*
       * Add to m_rrxn. m_rrxn is a vector of maps. m_rrxn has a length
       * equal to the total number of species for each species, there
       * exists a map, with the reaction number being the key, and the
       * reactant stoichiometric coefficient being the value.
       */
      m_rrxn[r.reactants[n]][rnum] = ns;
      for (m = 0; m < ns; m++) {
	rk.push_back(r.reactants[n]);
      }
    }
    /*
     * Now that we have rk[], we add it into the vector<vector_int> m_reactants
     * in the rnum index spot. Thus m_reactants[rnum] yields a vector
     * of reactants for the rnum'th reaction
     */
    m_reactants.push_back(rk);
        
    std::vector<size_t> pk;
    int np = r.products.size();
    for (n = 0; n < np; n++) {
      nsFlt = r.pstoich[n];
      ns = (int) nsFlt;
      if ((doublereal) ns != nsFlt) {
	if (ns < 1) ns = 1;
      }
      /*
       * Add to m_prxn. m_prxn is a vector of maps. m_prxn has a length
       * equal to the total number of species for each species, there
       * exists a map, with the reaction number being the key, and the
       * product stoichiometric coefficient being the value.
       */
      m_prxn[r.products[n]][rnum] = ns;
      for (m = 0; m < ns; m++) {
	pk.push_back(r.products[n]);
      }
    }
    /*
     * Now that we have pk[], we add it into the vector<vector_int> m_products
     * in the rnum index spot. Thus m_products[rnum] yields a vector
     * of products for the rnum'th reaction
     */
    m_products.push_back(pk);
    /*
     * Add this reaction to the stoichiometric coefficient manager. This
     * calculates rates of species production from reaction rates of 
     * progress.
     */
    m_rxnstoich.add(nReactions(), r);
    /*
     * register reaction in lists of reversible and irreversible rxns.
     */
    if (r.reversible) {
      m_revindex.push_back(nReactions());
      m_nrev++;
    } else {
      m_irrev.push_back( nReactions() );
      m_nirrev++;
    }        
  }


  /**
   * Prepare the class for the addition of reactions. This function
   * must be called after instantiation of the class, but before
   * any reactions are actually added to the mechanism.
   * This function calculates m_kk the number of species in all
   * phases participating in the reaction mechanism. We don't know
   * m_kk previously, before all phases have been added. 
   */
  void ElectrodeKinetics::init() {
    InterfaceKinetics::init();
  
  }

  /*
   * Finish adding reactions and prepare for use. This function
   * must be called after all reactions are entered into the mechanism
   * and before the mechanism is used to calculate reaction rates.
   *
   * Here, we resize work arrays based on the number of reactions,
   * since we don't know this number up to now.
   */
  void ElectrodeKinetics::finalize() {
    InterfaceKinetics::finalize();
  }


  bool ElectrodeKinetics::ready() const {
    return (m_finalized);
  }

  /**
   *  Install an individual reaction into a kinetics manager. The
   *  data for the reaction is in the xml_node r. In other words, r
   *  points directly to a ctml element named "reaction". i refers
   *  to the number id of the reaction in the kinetics object.
   * 
   * @param i Reaction number.
   * @param r XML_Node containing reaction data.
   * @param k Kinetics manager to which reaction will be added.
   * @param owning_phase Default phase for locating a species
   * @param rule Rule for handling reactions with missing species 
   *             (skip or flag as error)
   * @param validate_rxn If true, check that this reaction is not a 
   *                     duplicate of one already entered, and check that the reaction 
   *                     balances.
   *
   * @ingroup kineticsmgr
   */
  bool rxninfoEK::installReaction(int i, const XML_Node& reactionNode, ElectrodeKinetics* k, 
				  string owning_phase, const ReactionRules& rule,
				  bool validate_rxn) {

    ElectrodeKinetics& kin = *k;

    if (reactionNode.name() != "reaction") {
      throw CanteraError(" rxninfo::installReaction",
			 " expected xml node reaction, got " + reactionNode.name());
    }
    /*
     *  We use the ReactionData object to store initial values read
     *  in from the xml data. Then, when we have collected everything
     *  we add the reaction to the kinetics object, k, at the end
     * of the routine. (Someday this may be rewritten to skip building
     * the ReactionData object).
     */
    ReactionDataElectrode rdata;

    // Check to see if the reaction is specified to be a duplicate
    // of another reaction. It's an error if the reaction is a 
    // duplicate and this is not set.
    int dup = 0;
    if (reactionNode.hasAttrib("duplicate")) dup = 1;

    // Check to see if the reaction rate constant can be negative
    // It's an error if a negative rate constant is found and
    // this is not set.
    int negA = 0;
    if (reactionNode.hasAttrib("negative_A")) negA = 1;

    /*
     * This seemingly simple expression goes and finds the child element,
     * "equation". Then it treats all of the contents of the "equation"
     * as a string, and returns it the variable eqn. We post-process
     * the string to convert [ and ] characters into < and >, which 
     * cannot be stored in an XML file.
     * The string eqn is just used for IO purposes. It isn't parsed
     * for the identities of reactants or products.
     */
    string eqn = "<no equation>";
    if (reactionNode.hasChild("equation")) {
      eqn = reactionNode("equation");
    }
    int eqlen = static_cast<int>(eqn.size());
    int nn;
    for (nn = 0; nn < eqlen; nn++) {
      if (eqn[nn] == '[') eqn[nn] = '<';
      if (eqn[nn] == ']') eqn[nn] = '>';
    }

    /*
     * Get the reactants We store the id of products in rdata.reactants
     */
    bool ok = getReagents(reactionNode, kin, 1, owning_phase, rdata.reactants, 
			  rdata.rstoich, rdata.rorder, rule);

    /*
     * Get the products. We store the id of products in rdata.products
     */ 
    vector_fp dummy;
    ok = ok && getReagents(reactionNode, kin, -1, owning_phase, rdata.products, 
			   rdata.pstoich, dummy, rule);

    // if there was a problem getting either the reactants or the products, 
    // then abort.
    if (!ok) {
      return false;
    }

    // Check whether the reaction is specified to be reversible. Default is irreversible.
    // Check to see if the reaction is specified to be a duplicate
    // of another reaction. It's an error if the reaction is a 
    // duplicate and this is not set.
    rdata.reversible = false;
    string isrev = reactionNode["reversible"];
    if (isrev == "yes" || isrev == "true")
      rdata.reversible = true;

    /*
     * If reaction orders are specified, then this reaction
     * does not follow mass-action kinetics, and is not 
     * an elementary reaction. So check that it is not reversible,
     * since computing the reverse rate from thermochemistry only
     * works for elementary reactions. Set the type to global,
     * so that kinetics managers will know to process the reaction
     * orders.
     */
    if (reactionNode.hasChild("order")) {
      if (rdata.reversible == true) 
	throw CanteraError("installReaction",
			   "reaction orders may only be given for "
			   "irreversible reactions");
      rdata.global = true;
    }

    /*
     * Search the reaction element for the attribute "type".
     * If found, then branch on the type, to fill in appropriate
     * fields in rdata. 
     */ 
    rdata.reactionType = ELEMENTARY_RXN;
    string typ = reactionNode["type"];
    if (typ == "surface") {
      rdata.reactionType = SURFACE_RXN;
    }
    else if (typ == "edge") {
      rdata.reactionType = EDGE_RXN;
    }
    else if (typ != "") {
      throw CanteraError("installReaction",
			 "Unknown reaction type: " + typ);
    }
    /*
     * Look for undeclared duplicate reactions.
     */
    if (validate_rxn) {
      doublereal c = 0.0;

      map<int, doublereal> rxnstoich;
      rxnstoich.clear();
      int nr = rdata.reactants.size();
      for (nn = 0; nn < nr; nn++) {
	rxnstoich[-1 - rdata.reactants[nn]] -= rdata.rstoich[nn];
      }
      int np = rdata.products.size();
      for (nn = 0; nn < np; nn++) {
	rxnstoich[rdata.products[nn]+1] += rdata.pstoich[nn];
      }
      int nrxns = static_cast<int>(m_rdata.size());
      for (nn = 0; nn < nrxns; nn++) {
	if ((int(rdata.reactants.size()) == m_nr[nn]) 
	    && (rdata.reactionType == m_typ[nn])) {
	  c = isDuplicateReaction(rxnstoich, m_rdata[nn]);
	  if (c > 0.0 
	      || (c < 0.0 && rdata.reversible)
	      || (c < 0.0 && m_rev[nn])) {
	    if ((!dup || !m_dup[nn])) {
	      string msg = string("Undeclared duplicate reactions detected: \n")
		+"Reaction "+int2str(nn+1)+": "+m_eqn[nn]
		+"\nReaction "+int2str(i+1)+": "+eqn+"\n";
	      throw CanteraError("installReaction", msg);
	    }
	  }
	}
      }
      m_dup.push_back(dup);
      m_rev.push_back(rdata.reversible);
      m_eqn.push_back(eqn);
      m_nr.push_back(rdata.reactants.size());
      m_typ.push_back(rdata.reactionType);
      m_rdata.push_back(rxnstoich);
    }
        
    // Write information to the rdata object 
    rdata.equation = eqn;
    rdata.number = i;
    rdata.rxn_number = i;

    /*
     * Read the rate coefficient data from the XML file. Trigger an
     * exception for negative A unless specifically authorized.
     */
    rdata.getRateCoefficient(reactionNode.child("rateCoeff"), kin, negA);

    /*
     * Check to see that the elements balance in the reaction.
     * Throw an error if they don't
     */        
    if (validate_rxn) {
      checkRxnElementBalance(kin, rdata);
    }

    /*
     * Ok we have read everything in about the reaction. Add it
     * to the kinetics object by calling the Kinetics member function,
     * addReaction()
     */
    kin.addReaction(rdata);
    return true;
  }


}








