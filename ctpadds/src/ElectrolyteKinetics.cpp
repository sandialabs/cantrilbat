/**
 *  @file ElectrolyteKinetics.cpp 
 * Homogeneous kinetics in a liquid electrolyte phase.
 * 
 */
/*
 * $Id:
 */
/*
 * Copywrite 2007 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */


#include "ElectrolyteKinetics.h"
#include "cantera/thermo/IdealSolidSolnPhase.h"
#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/importKinetics.h"


using namespace std;

namespace Cantera {

  /*
   *
   * ElectrolyteKinetics():
   *
   * Constructor for ElectrolyteKinetics, starts with an empty reaction mechanism.
   * However, the thermo parameter is required. It makes no sense to 
   * construct
   * a reaction mechanism for a phase that doesn't have a ThermoPhase object
   * associated with it.
   */    
  ElectrolyteKinetics::
  ElectrolyteKinetics() :
    Kinetics(),
    m_kk(0),
    m_nirrev(0), 
    m_nrev(0),
    logC0AllTheSameConstant(false),
    logC0ProdVariable(true),
    m_finalized(false)
  {
    m_kdata = new ElectrolyteKineticsData;
  }

  /**************************************************************************
   *
   * ElectrolyteKinetics():
   *
   * Constructor for ElectrolyteKinetics, starts with an empty reaction mechanism.
   * However, the thermo parameter is required. It makes no sense to 
   * construct
   * a reaction mechanism for a phase that doesn't have a ThermoPhase object
   * associated with it.
   */    
  ElectrolyteKinetics::
  ElectrolyteKinetics(thermo_t* thermo_ptr) :
    Kinetics(),
    m_kk(0),
    m_nirrev(0), 
    m_nrev(0),
    logC0AllTheSameConstant(false),
    logC0ProdVariable(true),
    m_finalized(false)
  {
    m_kdata = new ElectrolyteKineticsData;

  
    if (!thermo_ptr) {
      throw CanteraError("ElectrolyteKinetics Constructor",
			 "Must supply a valid ThermoPhase object");
    }
    Kinetics::addPhase(*thermo_ptr);
    /**
     * query ThermoPhase object for the equation of state type.
     * Fill in the flags logC0ProdVariable and   logC0AllTheSameConstant 
     * based on its value.
     */
    int eosFlag = thermo().eosType();
    switch (eosFlag) {
    case cIdealSolidSolnPhase0:
      logC0ProdVariable = false;
      break;
    case cIdealSolidSolnPhase1:
    case cIdealSolidSolnPhase2:
      logC0ProdVariable = false;
      logC0AllTheSameConstant = true;
      break;
    }
  }

  /*************************************************************************
   *
   * ~ElectrolyteKinetics()
   *
   * Destructor for the ElectrolyteKinetics class
   */
  ElectrolyteKinetics::~ElectrolyteKinetics(){
    delete m_kdata;
  }

  /**************************************************************************
   *
   * getFwdRatesOfProgress()
   *
   * Returns a vector containing the forward rate of progress of the ith
   * reaction [kmol/m^3 s^1]. 
   *
   * @param netROP[i] On return it will contain the forward rate of progress
   *                  of the ith reaction [kmol/m^3 s^1]. Dimensioned
   *                  at least m_ii.
   */
  void ElectrolyteKinetics::getFwdRatesOfProgress(doublereal* fwdROP) { 
    updateROP(); 
    copy(m_kdata->m_ropf.begin(), m_kdata->m_ropf.end(), fwdROP);
  }

  /*************************************************************************
   *
   * getRevRatesOfProgress()
   *
   * Returns a vector containing the reverse rate of progress of the ith
   * reaction [kmol/m^3 s^1]. 
   *
   * @param netROP[i] On return it will contain the reverse rate of progress
   *                  of the ith reaction [kmol/m^3 s^1]. Dimensioned
   *                  at least m_ii.
   */
  void ElectrolyteKinetics::getRevRatesOfProgress(doublereal* revROP) { 
    updateROP(); 
    copy(m_kdata->m_ropr.begin(), m_kdata->m_ropr.end(), revROP);
  }

  /***************************************************************************
   *
   * getNetRatesOfProgress():
   *
   * Returns a vector containing the net rate of progress of the ith
   * reaction [kmol/m^3 s^1]. 
   *
   * @param netROP[i] On return it will contain the net rate of progress
   *                  of the ith reaction [kmol/m^3 s^1]. Dimensioned
   *                  at least m_ii.
   */
  void ElectrolyteKinetics::getNetRatesOfProgress(doublereal* netROP) { 
    updateROP(); 
    copy(m_kdata->m_ropnet.begin(), m_kdata->m_ropnet.end(), netROP);
  }

  /**************************************************************************
   *
   * getNetProductionRates():
   *
   * Species net production rates [kmol/m^3 s^1]. Return the species
   * net production rates (creation - destruction) in array
   * net, which must be dimensioned at least as large as the
   * total number of species.
   */
  void ElectrolyteKinetics::
  getNetProductionRates(doublereal* net) {
    /*
     * We do the work here. We calculate reactions rates of progress and
     * store them in the vector m_kdata->m_ropnet
     */
    updateROP();

    /*
     * Go call the stoichiometry managers to obtain the production
     * rates of the species.
     */   
    m_rxnstoich.getNetProductionRates(m_kk, 
				      &m_kdata->m_ropnet[0],
				      net);

    //	m_rxnstoich.getNetProductionRates(m_kk,
    //				  &m_kdata->m_ropnet[0],
    //				  net);
  }

  /**************************************************************************
   *
   * getCreationRates()
   *
   * Species net production rates [kmol/m^3]. Return the species
   * net production rates (creation - destruction) in array
   * net, which must be dimensioned at least as large as the
   * total number of species.
   */
  void ElectrolyteKinetics::getCreationRates(doublereal* cdot) {
    /*
     * Update the rates of progress of the reactions
     */
    updateROP();
    m_rxnstoich.getCreationRates(m_kk, &m_kdata->m_ropf[0], 
				 &m_kdata->m_ropr[0], cdot); 
  }

  /*************************************************************************
   *
   * getDestructionRates():
   *
   * Species destruction rates [kmol/m^3 s^1]. Return the species
   * destruction rates in array
   * ddot, which must be dimensioned at least as large as the
   * total number of species.
   */
  void ElectrolyteKinetics::getDestructionRates(doublereal* ddot) {
    /*
     * Update the rates of progress of the reactions
     */
    updateROP();
    m_rxnstoich.getDestructionRates(m_kk, &m_kdata->m_ropf[0], 
				    &m_kdata->m_ropr[0], ddot);
  }

  /************************************************************************
   *
   * update_T():
   *
   * Update temperature-dependent portions of reaction rates and
   * falloff functions.
   */
  void ElectrolyteKinetics::
  update_T() {
    _update_rates_T();
  }

  /************************************************************************
   *
   * update_C():
   *
   * Update concentration-dependent portions of reaction rates and
   * falloff functions.
   */
  void ElectrolyteKinetics::
  update_C() {
    _update_rates_C();
  }

  /**********************************************************************
   *
   * isReversible()
   *
   * This functions returns a boolean indicating whether a particular
   * reaction is reversible.
   *
   * @param i Reaction number
   */
  bool ElectrolyteKinetics::isReversible(size_t i) {
    if (find(m_revindex.begin(), m_revindex.end(), i) 
	< m_revindex.end()) return true;
    else return false;
  }

  /**********************************************************************
   *
   * _update_rates_T():
   *
   * Update the temperature dependent portions of the rate constants
   * Objects which are updated:
   *  m_rates - the rate constants. 
   *  m_kdata->m_rkc[] Inverse of the equilibrium constants
   */
  void ElectrolyteKinetics::
  _update_rates_T() {
    doublereal T = thermo().temperature();
    /*
     * Calculate logC0 if all the same.
     */
    if (logC0AllTheSameConstant) {
      m_kdata->m_logC0_scalar = log(thermo().standardConcentration(0));
    }

    doublereal logT = log(T);
    /**
     * Update the forward rate constants. We only update those
     * rate constants which have a temperature dependence in 
     * this step.
     */
    m_rates.update(T, logT, &m_kdata->m_rfn[0]);
    /**
     * Store the temperature at which the rate constants were
     * last evaluated.
     */
    m_kdata->m_temp = T;
    /*
     * Update the equilibrium constants for reversible reactions
     * only
     */
    updateKc();
    m_kdata->m_ROP_ok = false;
  };

  /***********************************************************************
   *
   * _update_rates_C():
   *
   * Update properties that depend on concentrations. Currently this
   * only involves the calculation of the generalized concentrations
   */         
  void ElectrolyteKinetics::
  _update_rates_C() {
    int np = nPhases();
    for (int n = 0; n < np; n++) {
      /*
       * We call the getActivityConcentrations function of each
       * ThermoPhase class that makes up this kinetics object to 
       * obtain the generalized concentrations for species within that 
       * class. This is collected in the vector m_conc. m_start[]
       * are integer indecises for that vector denoting the start of the
       * species for each phase.
       */
      thermo(n).getActivityConcentrations(&m_actConc[0] + m_start[n]);
    }
    m_kdata->m_ROP_ok = false;
  }

  /*********************************************************************
   *
   * updateKc():
   *
   * Update the inverse of the equilibrium constants in molar units. This is
   * carried out whenever _update_T() is called as the 
   * reverse rate constants are evaluated from the equalibrium
   * constants.
   * The inverse of the equilibrium constant is located in m_rkcn[].
   */
  void ElectrolyteKinetics::updateKc() {
    size_t i, n, k;

    /*
     *  Get the standard state chemical potentials of the species.
     *  This is the array of chemical potentials at unit activity 
     *  We define these here as the chemical potentials of the pure
     *  species at the temperature and pressure of the solution.
     */
    size_t nsp;
    size_t ik = 0;
    doublereal rt = GasConstant*thermo(0).temperature();
    doublereal rrt = 1.0/rt;
    size_t np = nPhases();
    for (n = 0; n < np; n++) {
      thermo(n).getStandardChemPotentials(&m_mu0[0] + m_start[n]);
      nsp = thermo(n).nSpecies();
      for (k = 0; k < nsp; k++) {
	m_mu0[ik] -= rt*thermo(n).logStandardConc(k);
	m_mu0[ik] += Faraday * m_phi[n] * thermo(n).charge(k);
	ik++;
      }
    }

    /*
     * Use the stoichiometric manager to find delta G^0 for each
     * reaction.
     * First, zero out m_rkc[]. Then fill it in for each reaction
     * whether reversible or not.
     */
    vector_fp& m_rkc = m_kdata->m_rkcn;
    m_rxnstoich.getRevReactionDelta(m_ii, &m_mu0[0], &m_rkc[0]);

    for (i = 0; i < m_nrev; i++) {
      size_t irxn = m_revindex[i];
      if (irxn < 0 || irxn >= nReactions()) {
	throw CanteraError("ElectrolyteKinetics",
			   "illegal value: irxn = "+int2str(irxn));
      }
      m_rkc[irxn] = exp(m_rkc[irxn]*rrt);
    }
	
    /*
     * Zero out the equilibrium constant for irreversible reactions.
     */
    for(i = 0; i != m_nirrev; ++i) {
      m_rkc[ m_irrevindex[i] ] = 0.0;
    }
  }

  /**************************************************************************
   *
   * getEqilibriumConstants():
   *
   * Get the equilibrium constants of all reactions, whether
   * reversible or not.
   *
   * @param kc[] - Vector of length m_ii which on return
   *               will contain the equilibrium constants
   *               for all reactions.
   */
  void ElectrolyteKinetics::
  getEquilibriumConstants(doublereal* kc) {
    int n, nsp, k, ik=0;
    doublereal rt = GasConstant*thermo(0).temperature();
    doublereal rrt = 1.0/rt;
    int np = nPhases();
    /*
     *  Get the standard state chemical potentials of the species.
     *  This is the array of chemical potentials at unit activity 
     *  We define these here as the chemical potentials of the pure
     *  species at the temperature and pressure of the solution.
     */
    for (n = 0; n < np; n++) {
      thermo(n).getStandardChemPotentials(&m_mu0[0] + m_start[n]);
      nsp = thermo(n).nSpecies();
      for (k = 0; k < nsp; k++) {
	m_mu0[ik] -= rt*thermo(n).logStandardConc(k);
	m_mu0[ik] += Faraday * m_phi[n] * thermo(n).charge(k);
	ik++;
      }
    }


    fill(kc, kc + m_ii, 0.0);

    //m_reactantStoich.decrementReactions(m_mu0.begin(), kc); 
    //m_revProductStoich.incrementReactions(m_mu0.begin(), kc);
    //m_irrevProductStoich.incrementReactions(m_mu0.begin(), kc);
    m_rxnstoich.getReactionDelta(m_ii, &m_mu0[0], kc);

    for (size_t i = 0; i < m_ii; i++) {
      kc[i] = exp(-kc[i]*rrt);
    }
  }

  /**************************************************************************
   *
   * getDeltaGibbs():
   *
   * Return the vector of values for the reaction gibbs free energy
   * change
   * These values depend upon the concentration
   * of the solution.
   *
   *  units = J kmol-1
   */
  void ElectrolyteKinetics::getDeltaGibbs(doublereal* deltaG) {
    /*
     * Get the chemical potentials of the species in the 
     * solid solution.Note this works as ElectrolyteKinetics is 
     * limited to one phase.
     */
    int np = nPhases();
    int n;
    for (n = 0; n < np; n++) {
      thermo(n).getChemPotentials(&m_grt[0] + m_start[n]);
    }
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, &m_grt[0], deltaG);
  }
    
  /**************************************************************************
   *
   * getDeltaEnthalpy():
   * 
   * Return the vector of values for the reactions change in
   * enthalpy.
   * These values depend upon the concentration
   * of the solution.
   *
   *  units = J kmol-1
   */
  void ElectrolyteKinetics::getDeltaEnthalpy(doublereal* deltaH) {
    /*
     * Get the partial molar enthalpy of all species in the 
     * ideal gas.
     */
    int np = nPhases();
    int n;
    for (n = 0; n < np; n++) {
      thermo(n).getPartialMolarEnthalpies(&m_grt[0] + m_start[n]);
    }
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, &m_grt[0], deltaH);
  }

  /************************************************************************
   *
   * getDeltaEntropy():
   *
   * Return the vector of values for the reactions change in
   * entropy.
   * These values depend upon the concentration
   * of the solution.
   *
   *  units = J kmol-1 Kelvin-1
   */
  void ElectrolyteKinetics::getDeltaEntropy( doublereal* deltaS) {
    /*
     * Get the partial molar entropy of all species in the
     * solid solution.
     */
    int np = nPhases();
    int n;
    for (n = 0; n < np; n++) {
      thermo(n).getPartialMolarEntropies(&m_grt[0] + m_start[n]);
    }
    /*
     * Use the stoichiometric manager to find deltaS for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, &m_grt[0], deltaS);
  }
    
  /***********************************************************************
   *
   * getDeltaSSGibbs():
   *
   * Return the vector of values for the reaction 
   * standard state gibbs free energy change.
   * These values don't depend upon the concentration
   * of the solution.
   *
   *  units = J kmol-1
   */
  void ElectrolyteKinetics::getDeltaSSGibbs(doublereal* deltaG) {
    /*
     *  Get the standard state chemical potentials of the species.
     *  This is the array of chemical potentials at unit activity 
     *  We define these here as the chemical potentials of the pure
     *  species at the temperature and pressure of the solution.
     */
    int np = nPhases();
    int n;
    for (n = 0; n < np; n++) {
      thermo(n).getStandardChemPotentials(&m_grt[0] + m_start[n]);
    }
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, &m_grt[0], deltaG);
  }

  /**********************************************************************
   *
   * getDeltaSSEnthalpy():
   *
   * Return the vector of values for the change in the
   * standard state enthalpies of reaction.
   * These values don't depend upon the concentration
   * of the solution.
   *
   *  units = J kmol-1
   */
  void ElectrolyteKinetics::getDeltaSSEnthalpy(doublereal* deltaH) {
    /*
     *  Get the standard state enthalpies of the species.
     *  This is the array of chemical potentials at unit activity 
     *  We define these here as the enthalpies of the pure
     *  species at the temperature and pressure of the solution.
     */
    size_t np = nPhases();
    for (size_t n = 0; n < np; n++) {
      thermo(n).getEnthalpy_RT(&m_grt[0] + m_start[n]);
    }
    doublereal RT = thermo().temperature() * GasConstant;
    for (size_t k = 0; k < m_kk; k++) {
      m_grt[k] *= RT;
    }
    /*
     * Use the stoichiometric manager to find deltaG for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, &m_grt[0], deltaH);
  }

  /*********************************************************************
   *
   * getDeltaSSEntropy():
   *
   * Return the vector of values for the change in the
   * standard state entropies for each reaction.
   * These values don't depend upon the concentration
   * of the solution.
   *
   *  units = J kmol-1 Kelvin-1
   */
  void ElectrolyteKinetics::getDeltaSSEntropy(doublereal* deltaS) {
    /*
     *  Get the standard state entropy of the species.
     *  We define these here as the entropies of the pure
     *  species at the temperature and pressure of the solution.
     */
    size_t np = nPhases();
    for (size_t n = 0; n < np; n++) {
      thermo(n).getEntropy_R(&m_grt[0] + m_start[n]);
    }
    doublereal R = GasConstant;
    for (size_t k = 0; k < m_kk; k++) {
      m_grt[k] *= R;
    }
    /*
     * Use the stoichiometric manager to find deltaS for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, &m_grt[0], deltaS);
  }

  void ElectrolyteKinetics::getActivationEnergies(doublereal *E) {
    copy(m_E.begin(), m_E.end(), E);
  }

  /*********************************************************************
   *
   * updateROP():
   *
   * Update the rate of progress for the reactions.
   * This key routine makes sure that the rate of progress vectors
   * located in the solid kinetics data class are up to date.
   */
  void ElectrolyteKinetics::updateROP() {

    _update_rates_T();
    _update_rates_C();

    /*
     * If the ROP's are up to date return without doing any
     * calculations -> The update rate coefficients above will
     * check to see if the temperature and concentrations have
     * changed. 
     */
    if (m_kdata->m_ROP_ok) return;

    const vector_fp& rf = m_kdata->m_rfn;
    const vector_fp& m_rkc = m_kdata->m_rkcn;
    vector_fp& ropf = m_kdata->m_ropf;
    vector_fp& ropr = m_kdata->m_ropr;
    vector_fp& ropnet = m_kdata->m_ropnet;

    /*
     * copy the forward rate coefficients, m_rfn, into forward
     * rate of progress vector, ropf, to start the process.
     */
    copy(rf.begin(), rf.end(), ropf.begin());

    /*
     * multiply the rate coefficient by the perturbation
     *  factors
     */
    multiply_each(ropf.begin(), ropf.end(), m_perturb.begin());
           
    /*
     * Copy the forward rate coefficients to the reverse rate of
     * progress vector, ropr
     */ 
    copy(ropf.begin(), ropf.end(), ropr.begin());
        
    // for reverse rates computed from thermochemistry, multiply
    // the forward rates copied into m_ropr by the reciprocals of
    // the equilibrium constants
    multiply_each(ropr.begin(), ropr.end(), m_rkc.begin());

    /*
     * Multiply the rate coefficients by the activity concentrations
     * to generate the forward rate of progress of the reaction
     */
    m_rxnstoich.multiplyReactants(&m_actConc[0], &ropf[0]);
  
    // for reversible reactions, multiply ropr by concentration
    // products
    m_rxnstoich.multiplyRevProducts(&m_actConc[0], &ropr[0]); 

    // do global reactions
    //m_globalReactantStoich.power(m_conc.begin(), ropf.begin());
  
    /**
     * Calculate the net rate of progress for a reaction from the
     * difference of the forward and the reverse.
     */
    for (size_t j = 0; j < m_ii; ++j) {
      ropnet[j] = ropf[j] - ropr[j];
    }

    /*
     * signal that the ROP's are up to date
     */
    m_kdata->m_ROP_ok = true;
  }

  /*********************************************************************
   *
   * getFwdRateConstants():
   *
   * Update the rate of progress for the reactions.
   * This key routine makes sure that the rate of progress vectors
   * located in the solid kinetics data class are up to date.
   */
  void ElectrolyteKinetics::
  getFwdRateConstants(doublereal *kfwd) {
    _update_rates_T();
    _update_rates_C();
    const vector_fp& rf = m_kdata->m_rfn;
    for (size_t i = 0; i < m_ii; i++) {
      kfwd[i] = rf[i] * m_perturb[i];
    }
  }

  /*********************************************************************
   *
   * getRevRateConstants():
   *
   * Return a vector of the reverse reaction rate constants
   *
   * Length is the number of reactions. units depends
   * on many issues. Note, this routine will return rate constants
   * for irreversible reactions if the default for
   * doIrreversible is overridden.
   */
  void ElectrolyteKinetics::
  getRevRateConstants(doublereal *krev, bool doIrreversible) {
    _update_rates_T();
    _update_rates_C();
    const vector_fp& rf = m_kdata->m_rfn;
    if (doIrreversible) {
      doublereal *tmpKc = &m_kdata->m_ropnet[0];
      getEquilibriumConstants(tmpKc);
      for (size_t i = 0; i < m_ii; i++) {
	krev[i] = rf[i] * m_perturb[i] / tmpKc[i];
      }
    } else {
      /*
       * m_rkc[] is zero for irreversibly reactions
       */
      const vector_fp& m_rkc = m_kdata->m_rkcn;
      for (size_t i = 0; i < m_ii; i++) {
	krev[i] = rf[i] * m_perturb[i] * m_rkc[i];
      }
    }
  }

  /*********************************************************************
   *
   * addReaction():
   *
   * Add a reaction to the internal reaction mechanism for this class.
   * Calls to this routine must be made after the call to init()
   * and before the call to finalize.
   */
  void ElectrolyteKinetics::
  addReaction(ReactionData& r) {
    Kinetics::addReaction(r);
  }

  /***********************************************************************
   *
   * addElementaryReaction():
   *
   * Add a simple Arrhenious reaction to the reaction mechanism
   * This routine installs an entry for the reaction in the rate coefficient
   * calculator object, m_rates.
   * It initializes and resizes m_rfn, and
   * adds an entry in the information map, m_index. 
   */
  void ElectrolyteKinetics::
  addElementaryReaction(const ReactionData& r) {
//    int iloc;

    // install rate coeff calculator
 //   iloc = m_rates.install(reactionNumber(),
//			   r.rateCoeffType, r.rateCoeffParameters.size(), 
//			   &r.rateCoeffParameters[0] );
    size_t iloc = m_rates.install(reactionNumber(), r);

    /*
     * add constant term to rate coeff value vector.
     * Note, for constant reaction rates coefficients, this is the
     * one and only place where the rate coefficient is installed
     * into the rate coefficient array.
     */
    m_kdata->m_rfn.push_back(r.rateCoeffParameters[0]);                

    /*
     * Fill in the m_index entry for this reaction. This is an
     * information map, keyed on the reaction number,
     * describing the reaction type (always ELEMENTARY_RXN for this
     * class) and the index of the Reaction Rate Coefficient Calculator. 
     */
    registerReaction(reactionNumber(), ELEMENTARY_RXN, iloc);
  }
   
  /***********************************************************************
   *
   * installReagents():
   *
   * Add the species reactants and products to the various
   * arrays.
   * The following arrays are initialized for the current
   * reaction:
   *  m_rrxn - Stoichiometric coefficients for reactants
   *  m_prxn - Stoichiometric coefficients for products.
   *  m_reactants - List of reactants (Kinetics Class)
   *  m_products  - List of products  (Kinetics Class)
   *  m_reactantStoich - Stoichiometric manager for reactants
   *  m_revProductStoich - Stoichiometric manager for products
   *                       of reversible reactions.
   *  m_irrevProductStoich Stoichiometric manager for products
   *                       of irreversible reactions.
   *  m_dn - Vector containing the net change in moles of the
   *         ith reaction
   *  m_revindex - Vector containing a listing of the indecices
   *               of all the reversible reactions.
   *  m_irrevindex -Vector containing a listing of the indecices
   *               of all the irreversible reactions. 
   */
  void ElectrolyteKinetics::installReagents(const ReactionData& r) {
    int n, m;
    double ns;
    /*
     * Extend the rates of progress vectors by one, as we are
     * adding a new reaction. The rops exist in the
     * ElectrolyteKineticsData class.
     */
    m_kdata->m_ropf.push_back(0.0);
    m_kdata->m_ropr.push_back(0.0);
    m_kdata->m_ropnet.push_back(0.0);
    m_kdata->m_rkcn.push_back(0.0);
    m_kdata->m_logProdC0.push_back(0.0);
    /*
     * Go get the total number of reactions in the mechanism
     * added to date (m_ii)
     */
    int rnum = reactionNumber();

    /*
     * FIRST, WE WILL TAKE CARE OF STORAGE FOR THE REACTANTS
     */
    /*
     * rk will be used as input to the stoichiometric
     * manager
     */
    vector<size_t> rk;
    /*
     * look up the number of reactants in the reaction, nr.
     * and then loop over them
     */
    int nr = r.reactants.size();
    for (n = 0; n < nr; n++) {
      /*
       * Look up the stoichiometric coefficient for the
       * species in the current reaction
       */
      ns = r.rstoich[n];
      /*
       * Add the current species to the stoichiometric
       * coefficient reactants. (Q: m_rrxn the right size?)
       */
      if (ns != 0.0) m_rrxn[r.reactants[n]][rnum] = ns;
      /*
       * Create a number integer vector entries equal to the
       * stoichiometric coefficient for the current reactant
       */
      for (m = 0; m < ns; m++) {
	rk.push_back(r.reactants[n]);
      }
    }
    /*
     * Add an entry for m_reactants, which is a member
     * of the Kinetics class.
     */
    m_reactants.push_back(rk);

    /*
     * NEXT, WE WILL TAKE CARE OF STORAGE FOR THE PRODUCTS
     */
    vector<size_t> pk;
    /*
     * look up the number of products in the reaction, np.
     * and then loop over them
     */	
    int np = r.products.size();
    for (n = 0; n < np; n++) {
      /*
       * Look up the stoichiometric coefficient for the
       * current species product in the current reaction
       */
      ns = r.pstoich[n];
      /*
       * Add the current species to the stoichiometric
       * coefficient products. (Q: m_prxn the right size?)
       */
      if (ns != 0) m_prxn[r.products[n]][rnum] = ns;
      /*
       * Create a number of integer vector entries equal to the
       * stoichiometric coefficient for the current product
       */
      for (m = 0; m < ns; m++) {
	pk.push_back(r.products[n]);
      }
    }
    /*
     * Add an entry for m_reactants, which is a member
     * of the Kinetics class.
     */
    m_products.push_back(pk);

    /*
     * Add this reaction to the stoichiometric coefficient manager. This
     * calculates rates of species production from reaction rates of 
     * progress.
     */
    m_rxnstoich.add( reactionNumber(), r);

    /*
     * m_dn - fill in the entry for the net change of moles
     *        due to this reaction. 
     */
    m_dn.push_back(pk.size() - rk.size());

    /*
     * Branch on whether this reaction is reversible
     * -> fill in entries in the stoichiometric manager
     *    and in the listing of reversable and irreversible
     *    reactions.
     */
    if (r.reversible) {
      m_revindex.push_back(reactionNumber());     
      m_nrev++;
    }
    else {
      m_irrev.push_back( reactionNumber() );    
      m_nirrev++;
    }
  }

  /*********************************************************************
   *
   * installGroups()
   *
   * Save the reaction and product groups, which are
   * part of the ReactionData class, in this class.
   * They aren't used for anything but reaction path
   * analysis, as far as I can figure out.
   */
  void ElectrolyteKinetics::installGroups(size_t irxn, 
					  const vector<grouplist_t>& r, const vector<grouplist_t>& p) {
    if (!r.empty()) {
      m_rgroups[reactionNumber()] = r;
      m_pgroups[reactionNumber()] = p;
    }
  }
    
  /*********************************************************************
   *
   * registerReaction():
   *
   * Creates an information map about the ith reaction.
   * For the ith reaction,
   * the reaction type and the iloc location in the m_rate array
   * are saved as a pair object in the internal variable m_index.
   */
  void ElectrolyteKinetics::registerReaction(size_t rxnNumber, 
					     int type, int loc) {
    m_index[rxnNumber] = pair<int, int>(type, loc);
  }

  /*********************************************************************
   *
   * init()
   *
   * Prepare the class for the addition of reactions. This function
   * must be called after instantiation of the class, but before
   * any reactions are actually added to the mechanism.
   *
   * This function currently resizes all arrays based on the number
   * of species in the phase.
   */
  void ElectrolyteKinetics::init() { 
    m_kk = thermo().nSpecies();
    if (m_kk <= 0) {
      throw CanteraError("ElectrolyteKinetics::init",
			 "m_kk is zero or less");
    }
    m_rrxn.resize(m_kk);
    m_prxn.resize(m_kk);
    m_actConc.resize(m_kk);
    m_grt.resize(m_kk);
    m_kdata->m_logC0_vector.resize(m_kk);
  }


 
  // Initialization of a ElectrolyteKinetics reaction mechanism
  // from an XML file.
  /*!
   *
   * This routine is a precursor to importMechanism(XML_Node)
   * routine, which does most of the work.
   *
   * @param infile XML file containing the description of the
   *        phase node
   *
   * @param id  Optional parameter identifying the name of the
   *            phase. If none is given, the first XML
   *            phase element will be checked and used if
   *            it matched the ThermoPhase object.
   *              defaults to "" in the .h file.
   */
  void ElectrolyteKinetics::
  importMechanism(string inputFile, string id) {
    string path = findInputFile(inputFile);
    ifstream fin(path.c_str());
    if (!fin) {
      throw CanteraError("ElectrolyteKinetics::importMechanism",
			 "could not open "
			 +path+" for reading.");
    }
    /*
     * Create a top lvl XML_Node and then fill it with the
     * contents of the file
     */
    XML_Node *fxml = new XML_Node();
    fxml->build(fin);
    /*
     * jump down the XML tree to the XML element named phase
     * with the correct id.
     */
    XML_Node *fxml_phase = findXMLPhase(fxml, id);
    /*
     * Call the import mechanism function now that we
     * at the correct spot.
     */
    importMechanism(*fxml_phase, id);
    /*
     * Tidy up
     */
    delete fxml;
  }

  /*********************************************************************
   *
   * importMechanism():
   * This routine will call init(), get all of the reactions
   * from the XML data file, and then call finalize.
   * Note, it is expected that the thermo has already been initialized
   * in the ThermoPhase object and the pointer set correctly in the
   * Kinetics object.
   *
   *         id defaults to "" in the .h file.
   */
  void ElectrolyteKinetics::
  importMechanism(XML_Node& phaseNode, string id) {
    string idp = phaseNode.id();
    if (id.size() > 0) {
      if (idp != id) {
	throw CanteraError("initThermo", 
			   "phasenode and Id are incompatible");
      }
    }

    /**
     * Check to see that the current XML phase is the same
     * phase that was used to initialize the ThermoPhase object
     * which was used to initialize the Kinetics object. If there
     * is an incompatibility, throw an exception.
     */
    int index = phaseIndex(idp);
    if (index != 0) {
      throw CanteraError("ElectrolyteKinetics::importMechanism",
			 "phaseNode set to incompatible phase");
    }

    // if other phases are involved in the reaction mechanism,
    // they must be listed in a 'phaseArray' child
    // element. Homogeneous mechanisms do not need to include a
    // phaseArray element.
    /*
     * This Kinetics Object doesn't handle other phases.
     * If there is a phaseArray object, throw an exception
     */
    if (phaseNode.hasChild("phaseArray")) {
      throw CanteraError("ElectrolyteKinetics::importMechanism", 
			 "phaseArray element not allowed" );
    }

    /*
     * Look for a child of the xml element phase called
     * "kinetics". It has an attribute name "model".
     * Store the value of that attribute in the variable kintype
     */
    string kintype = phaseNode.child("kinetics")["model"];
    if (kintype != "ElectrolyteKinetics") {
      throw CanteraError("ElectrolyteKinetics::importMechanism", 
			 "kinetics model must have model "
			 "attribute ElectrolyteKinetics");
	  
    }

    /*
     * Look for compatibility of submodel types between the thermo
     * and the kinetics.
     */
    string kinsubtype;
    if (phaseNode.hasChild("standardConc")) {
      kinsubtype = phaseNode.child("standardConc")["model"];
    } else {
      if (phaseNode.hasChild("kinetics")) {
	XML_Node& kinNode = phaseNode.child("kinetics");
	if (kinNode.hasChild("standardConc")) {
	  kinsubtype = kinNode.child("standardConc")["model"];
	} else if (kinNode.hasAttrib("submodel")) {
	  kinsubtype = kinNode["submodel"];
	}
      }
    }
    int eosT = thermo().eosType();
    switch (eosT) {
    case cIdealSolidSolnPhase0:
      if (kinsubtype != "unity") {
	throw CanteraError("ElectrolyteKinetics::importMechanism", 
			   "kinetics submodel must have submodel "
			   "attribute \"unity\"");
      }
      break;
    case cIdealSolidSolnPhase1:
      if (kinsubtype != "molar_volume") {
	throw CanteraError("ElectrolyteKinetics::importMechanism", 
			   "kinetics submodel must have submodel "
			   "attribute \"molar volume\"");
      }	    
      break;
    case cIdealSolidSolnPhase2:
      if (kinsubtype != "solvent_volume") {
	throw CanteraError("ElectrolyteKinetics::importMechanism", 
			   "kinetics submodel must have submodel "
			   "attribute \"solvent volume\"");
      }
      break;
    default:
      throw CanteraError("ElectrolyteKinetics::importMechanism", 
			 "Unknown standard concentration model");
      break;
    }

    /*
     * allocate arrays based on the number of species
     */
    init();

    /*
     * Install the reactions.
     * This function also calls finalize
     */
    bool ok = installReactionArrays(phaseNode, *this, idp);
    if (!ok) {
      throw CanteraError("ElectrolyteKinetics::importMechanism",
			 "installReactionArrays returned an error flag");
    }
  }

  /********************************************************************
   *
   * finalize():
   *
   * Finish adding reactions and prepare for use. This function
   * must be called after all reactions are entered into the mechanism
   * and before the mechanism is used to calculate reaction rates. 
   *
   * member data initialized here:
   *  m_rstoich[i][k] -
   *  m_rstoich[i][k] -
   *  m_kdata->m_logC0_vector
   *  m_kdata->m_logC0_scalar
   *  m_kdata->m_logProdC0
   */
  void ElectrolyteKinetics::finalize() {
    if (!m_finalized) {
      int j, nr, np;
 
      for (size_t i = 0; i < m_ii; i++) {
	nr = m_reactants[i].size();
	for (j = 0; j < nr; j++) {
	  m_rstoich[i][m_reactants[i][j]]++;
	}
	np = m_products[i].size();
	for (j = 0; j < np; j++) {
	  m_pstoich[i][m_products[i][j]]++;
	}
      }

      /*
       * Provide initial values for m_logC0_vector, the vector
       * of standard concentrations for the species. This may
       * or may not be overwritten as the conditions change
       * depending upon the formulation of the ThermoPhase
       * object.
       */
      vector_fp& logC0_vector = m_kdata->m_logC0_vector;
      for (size_t i = 0; i < m_kk; ++i) {
	logC0_vector[i] = thermo().logStandardConc(i);
      }
      m_kdata->m_logC0_scalar =  logC0_vector[0];


      m_finalized = true;
    }
  }
  /********************************************************************
   *
   * ready():
   *
   *  Indicate whether the reaction mechanism is ready for use.
   */
  bool ElectrolyteKinetics::ready() const {
    return (m_finalized);
  }

  /*********************************************************************/
}








