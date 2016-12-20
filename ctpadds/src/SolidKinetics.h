/**
 * @file SolidKinetics.h
 *
 */
/*
 * $Author: hkmoffa $
 * $Revision: 508 $
 * $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef CT_SOLIDKINETICS_H
#define CT_SOLIDKINETICS_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/Array.h"
#include "cantera/base/xml.h"
#include "cantera/kinetics/BulkKinetics.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/StoichManager.h"
#include "cantera/kinetics/RateCoeffMgr.h"
#include "cantera/kinetics/importKinetics.h"

#include <map>

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif 
{

  // forward references

  //const int cSolidKinetics = 999;

  /**
   * Holds mechanism-specific data.
   */
  class SolidKineticsData {
  public:
    SolidKineticsData() :
      m_logC0_scalar(0.0),
      m_ROP_ok(false), 
      m_temp(0.0)
    {
    }
    virtual ~SolidKineticsData(){}

    /**
     * When the standard concentrations are constants and equal
     * for all species, this entry will contain the log of the
     * common standard concentration.
     */
    doublevalue m_logC0_scalar;

    /**
     * Rate of progress vector for the forward direction of each
     * reaction, length = m_ii 
     */
    vector_fp m_ropf;

    /**
     * Rate of progress vector for the reverse direction of each
     * reaction, length = m_ii 
     */
    vector_fp m_ropr;

    /**
     * Rate of progress vector for net direction of each
     * reaction, length = m_ii 
     * m_ropnet = m_ropf - m_ropr
     */
    vector_fp m_ropnet;

    /**
     * Boolean indicating that the ROP vectors are up to date
     */
    bool m_ROP_ok;

    /**
     * Last temperature at which the ROP vectors were evalulated at
     */
    doublevalue m_temp;

    /**
     * Value of the forward rate constants at the current
     * temperature. length number of reactions, mii.
     * This vector is updated by _update_T().
     * This vector is initialized in addReaction() with the
     * preexponential part of the rate constant. For reactions which
     * don't have any temperature dependent parts, the entry is not
     * updated during any subsequent processing. This is suppose to
     * save computation by cutting down on the number of exponentials
     * evaluated, but at the cost of additional indirection.
     */
    vector_fp m_rfn;

    /**
     * Value of the inverse of the equlibrium constant for each reaction
     * length = m_ii, the total number of reactions.
     */
    vector_fp m_rkcn;

    /**
     * Value of the logarithm of the products of the powers of the
     * standard concentrations with respect to their stoichiometric 
     * coefficients.
     * length = m_ii, the total number of reactions.
     */
    vector_fp m_logProdC0;

    /*
     * Vector of standard concentrations under the current conditions.
     * length = m_kk, the total number of species
     */
    vector_fp m_logC0_vector;
  };

//=========================================================================================================
  /**
   * Kinetics manager for elementary solid-phase chemistry. This
   * kinetics manager implements standard mass-action reaction rate
   * expressions for condensed phases. It assumes that all
   * stoichiometric coefficients are integers.
   *
   *  NOTE: THIS IS NOT USED. This is just a holder or a stub for functionality that may be added in.
   *        I don't think there is anything different implemented here for solids than the standard
   *        Cantera treatment.
   */
  class SolidKinetics : public BulkKinetics {

  public:

   SolidKinetics();
    /**
     * @name Constructors and General Information about Mechanism
     */
    //@{
    /// Constructor.
   SolidKinetics(thermo_t* thermo);

   SolidKinetics(const SolidKinetics& right);


   SolidKinetics& operator=(const SolidKinetics& right);

    /// Destructor.
    virtual ~SolidKinetics();

    Kinetics* duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const;

    /**
     * Return the ID number for this class
     */
    virtual int ID() const { return cSolidKinetics; }

    //@}
    /**
     * @name Reaction Rates Of Progress
     */
    //@{

    /**
     * Returns a vector containing the forward rate of progress of the ith
     * reaction [kmol/m^3 s^1]. 
     *
     * @param netROP[i] On return it will contain the forward rate of progress
     *                  of the ith reaction [kmol/m^3 s^1]. Dimensioned
     *                  at least m_ii.
     */
    virtual void getFwdRatesOfProgress(doublevalue* fwdROP);

    /**
     * Returns a vector containing the reverse rate of progress of the ith
     * reaction [kmol/m^3 s^1]. 
     *
     * @param netROP[i] On return it will contain the reverse rate of progress
     *                  of the ith reaction [kmol/m^3 s^1]. Dimensioned
     *                  at least m_ii.
     */
    virtual void getRevRatesOfProgress(doublevalue* revROP);

    /**
     * Returns a vector containing the net rate of progress of the ith
     * reaction [kmol/m^3 s^1]. 
     *
     * @param netROP[i] On return it will contain the net rate of progress
     *                  of the ith reaction [kmol/m^3 s^1]. Dimensioned
     *                  at least m_ii.
     */
    virtual void getNetRatesOfProgress(doublevalue* netROP);

    /*
     * Get the equilibrium constants of all reactions, whether
     * reversible or not.
     *
     * @param kc[] - Vector of length m_ii which on return
     *               will contain the equilibrium constants
     *               for all reactions.
     */
    virtual void getEquilibriumConstants(doublevalue* kc);

    /**
     * Return the vector of values for the reaction gibbs free energy
     * change.
     * These values depend upon the concentration
     * of the solution.
     *
     *  units = J kmol-1
     */
    virtual void getDeltaGibbs(doublevalue* deltaG);

    /**
     * Return the vector of values for the reactions change in
     * enthalpy.
     * These values depend upon the concentration
     * of the solution.
     *
     *  units = J kmol-1
     */
    virtual void getDeltaEnthalpy(doublevalue* deltaH);

    /**
     * Return the vector of values for the reactions change in
     * entropy.
     * These values depend upon the concentration
     * of the solution.
     *
     *  units = J kmol-1 Kelvin-1
     */
    virtual void getDeltaEntropy(doublevalue* deltaS);

    /**
     * Return the vector of values for the reaction 
     * standard state gibbs free energy change.
     * These values don't depend upon the concentration
     * of the solution.
     *
     *  units = J kmol-1
     */
    virtual void getDeltaSSGibbs(doublevalue* deltaG);

    /**
     * Return the vector of values for the change in the
     * standard state enthalpies of reaction.
     * These values don't depend upon the concentration
     * of the solution.
     *
     *  units = J kmol-1
     */
    virtual void getDeltaSSEnthalpy(doublevalue* deltaH);

    /**
     * Return the vector of values for the change in the
     * standard state entropies for each reaction.
     * These values don't depend upon the concentration
     * of the solution.
     *
     *  units = J kmol-1 Kelvin-1
     */
    virtual void getDeltaSSEntropy(doublevalue* deltaS);

    //@}
    /**
     * @name Species Production Rates
     */
    //@{

    /**
     * Species net production rates [kmol/m^3 s^1]. Return the species
     * net production rates (creation - destruction) in array
     * net, which must be dimensioned at least as large as the
     * total number of species.
     */
    virtual void getNetProductionRates(doublevalue* net);

    /**
     * Species creation rates [kmol/m^3 s^1]. Return the species
     * creation rates in array
     * cdot, which must be dimensioned at least as large as the
     * total number of species.
     */
    virtual void getCreationRates(doublevalue* cdot);

    /**
     * Species destruction rates [kmol/m^3 s^1]. Return the species
     * destruction rates in array
     * ddot, which must be dimensioned at least as large as the
     * total number of species.
     */
    virtual void getDestructionRates(doublevalue* ddot);

    //@}
    /**
     * @name Reaction Mechanism Informational Query Routines
     */
    //@{



    /**
     *  Return the reactant stoichiometric coefficient for the
     *  kth species in the ith reaction
     */
    //virtual doublevalue reactantStoichCoeff(size_t k, size_t i) const {
     // return getValue(m_rrxn[k], i, 0.0);

      //return Kinetics::m_rrxn[k][i];
   // }

    /**
     *  Return the product stoichiometric coefficient for the
     *  kth species in the ith reaction
     */
   // virtual doublevalue productStoichCoeff(size_t k, size_t i) const {
   //   return m_prxn[k][i];
   // }

    /**
     *  Return the net stoichiometric coefficient for the
     *  kth species in the ith reaction
     */
    //virtual doublevalue netStoichCoeff(size_t k, size_t i) const {
     // return (m_prxn[k][i] - m_rrxn[k][i]);
    //}

    /**
     * Return the forward rate constants
     *
     * length is the number of reactions. units depends
     * on many issues.
     */
    virtual void getFwdRateConstants(doublevalue *kfwd);

    /**
     * Return the reverse rate constants.
     *
     * length is the number of reactions. units depends
     * on many issues. Note, this routine will return rate constants
     * for irreversible reactions if the default for
     * doIrreversible is overridden.
     */
    virtual void getRevRateConstants(doublevalue *krev, bool doIrreversible = false);

    /**
     * Prepare the class for the addition of reactions. This function
     * must be called after instantiation of the class, but before
     * any reactions are actually added to the mechanism.
     */
    virtual void init();

    /**
     * Initialization of a SolidKinetics reaction mechanism
     * from an XML file.
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
     */
    void importMechanism(std::string infile, std::string id="");

    /**
     * This routine will call init(), get all of the reactions
     * from the XML data file, and then call finalize
     */
    void importMechanism(XML_Node& phaseNode, std::string id="");

    /**
     *   Add a reaction to the mechanism.
     * 
     *  This routine is called from deep within importCTML.cpp routines
     *  that read the xml file. 
     */
    virtual void addReaction(ReactionData& r);

    /**
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
    virtual void finalize();

    /**
     *  Indicate whether the reaction mechanism is ready for use.
     */
    virtual bool ready() const;

 
    /**
     * Update concentration portion of reaction rate coefficients
     */
    virtual void update_C();

    /**
     * Update the rate of progress for the reactions.
     * This key routine makes sure that the rate of progress vectors
     * located in the solid kinetics data class are up to date.
     */
    void updateROP();

    /**
     * Return the reaction type of the ith reaction
     */
    // virtual int reactionType(size_t i) const {
    //  return m_index[(size_t) i].first;
    // }

    /**
     * Return the index number in the Rate coefficient
     * calculator for the ith reaction
     */
    // virtual size_t reactionCalculatorIndex(size_t i) const {
    //  return m_index[i].second;
    // }

    /** 
     * Function returns a string containing the expression for
     * the reactions (i.e., reactants <=> products)
     *
     * @param i Reaction number
     */
    virtual std::string reactionString(size_t i) const {
      return m_rxneqn[i];
    }

    /**
     * Return the reaction path groups for the ith reaction's
     * reactants
     */
    //  const std::vector<grouplist_t>& reactantGroups(size_t i)
    // { return m_rgroups[i]; }

    /**
     * Return the reaction path groups for the ith reaction's
     * products
     */
    // const std::vector<grouplist_t>& productGroups(size_t i)
    //{ return m_pgroups[i]; }

    /**
     * This functions returns a boolean indicating whether a particle
     * reaction is reversible.
     *
     * @param i Reaction number
     */
    virtual bool isReversible(size_t i);

    /**
     * Update the temperature dependent portions of the rate constants
     * Objects which are updated:
     *  m_rates - the rate constants. 
     *  m_kdata->m_rkc[] Inverse of the equilibrium constants
     */
    void _update_rates_T();

    /*
     * Update properties that depend on concentrations. Currently this
     * only involves the calculation of the generalized concentrations
     */         
    void _update_rates_C();

  /**
     * Save the reaction and product groups, which are
     * part of the ReactionData class, in this class.
     * They aren't used for anything but reaction path
     * analysis, as far as I can figure out.
     */
    //   void installGroups(size_t irxn, const std::vector<grouplist_t>& r,
//		       const std::vector<grouplist_t>& p);
	
    /*
     * Update the inverse of the equilibrium constants in molar units. 
     * This is
     * carried out whenever _update_T() is called as the
     * reverse rate constants are evaluated from the equilibrium
     * constants.
     * The inverse of the equilibrium constant is located 
     * internally in m_rkcn[].
     */
    void updateKc();

    /**
     * Creates an information map about the ith reaction.
     * For the ith reaction,
     * the reaction type and the iloc location in the m_rate array
     * are saved as a pair object in the internal variable m_index.
     */
    //  void registerReaction(size_t rxnNumber, int type, int loc);

 private:

    /**
     * Returns the total number of reactions in the mechanism
     */
    //  size_t reactionNumber(){ return m_ii;}

    /**
     * Add a simple Arrhenious reaction to the reaction mechanism
     * This routine installs an entry for the reaction in the rate coefficient
     * calculator object, m_rates.
     * It initializes and resizes m_rfn, initializes m_fwrdOrder, and
     * adds an entry in the information map, m_index. 
     */
    // void addElementaryReaction(const ReactionData& r);
        
    /**
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
    // void installReagents(const ReactionData& r); 
  
  

  protected:

    /**
     * Number of species in the phase.
     */
    // size_t                                 m_kk;

    /**
     * This objects contains all of the rate constants for all of the
     * reactions in the mechanism. The Rate1 template class
     * is located in RateCoeffMgr.h. The Arrhenius class which the Rate1
     * class uses as a template parameter is located in
     * RxnRate.h.
     */
    // Rate1<Arrhenius>                    m_rates;        
        
    /**
     * Information map about the ith reaction. For the ith reaction,
     * the reaction type and the iloc location in the m_rate array
     * are saved as a pair object.
     */
    // mutable std::map<size_t, std::pair<int, size_t> >   m_index;

    /**
     * Stoichiometric Coefficent Manager for the reactants.
     * All of the reactants are listed in this object. The object
     * is used to convert from reaction rates of progress to 
     * species production rates.
     */
    // StoichManagerN                      m_reactantStoich;
    /**
     * Stoichiometric Coefficent Manager for the products from
     * reversible reactions
     * All of the products from reversible reactions
     * are listed in this object. The object
     * is used to convert from reaction rates of progress to 
     * species production rates.
     */
    // StoichManagerN                      m_revProductStoich;
    /**
     * Stoichiometric Coefficent Manager for the products from
     * irreversible reactions.
     * All of the products from irreversible reactions
     * are listed in this object. The object
     * is used to convert from reaction rates of progress to 
     * species production rates.
     */
    // StoichManagerN                      m_irrevProductStoich;

    /**
     * Number of irreversible reactions in the mechanism
     */
    int m_nirrev;
    /**
     * Number of reversible reactions in the mechanism
     */
    int m_nrev;

    /**
     * Map containing reaction path information for
     * each reaction's reactants
     */
    //  std::map<size_t, std::vector<grouplist_t> > m_rgroups;

    /**
     * Map containing reaction path information for
     * each reaction's productss
     */
    //std::map<size_t, std::vector<grouplist_t> >      m_pgroups;

    /**
     * This vector of maps contains the stoichiometric
     * coefficient information for the reactants.
     * The vector has length number of species, m_kk,
     * in the phase.
     * Each vector entry contains a map. The map's
     * key, an int, is the reaction index.
     * The map's value, a double, is the stoichiometric coefficient
     * of that species as a reactant in that reaction.
     * Thus,
     *        m_rrxn[k][i] = ns
     *
     *  species k in reaction i has a reactant stoichiometric
     *  coefficient of ns.
     */
    //mutable std::vector<std::map<size_t, doublevalue> > m_rrxn;

    /**
     * This vector of maps contains the stoichiometric
     * coefficient information for the products.
     * The vector has length number of species, m_kk,
     * in the phase.
     * Each vector entry contains a map. The map's
     * key, an int, is the reaction index.
     * The map's value, a double, is the stoichiometric coefficient
     * of that species as a product in that reaction.
     * Thus,
     *        m_prxn[k][i] = ns
     *
     *  species k in reaction i has a product stoichiometric
     *  coefficient of ns.
     */
    //mutable std::vector<std::map<size_t, doublevalue> >     m_prxn;

    /*
     * Vector of length m_ii, containing the net change in
     * moles of the ith reaction
     */
    vector_int m_dn;
    /*
     * Vector of length mrev, containing a listing of the
     * reactions which are reversible
     */
    //vector_int m_revindex;
    /*
     * Vector of length mirrev, containing a listing of the
     * reactions which are irreversible
     */  
    //  std::vector<int> m_irrevindex;

    /**
     * Stoichiometric matrix in a reaction-first format.
     *
     *  m_rstoich[i][k] is the reactant stoichiometric
     *  coefficient for species k in reaction i.
     *
     */
    std::map<size_t, std::map<size_t, doublevalue> >  m_rstoich;

    /**
     * Stoichiometric matrix in a product-first format.
     *
     *  m_pstoich[i][k] is the product stoichiometric
     *  coefficient for species k in reaction i.
     */
    std::map<size_t, std::map<size_t, doublevalue> >  m_pstoich;

    /**
     * Vector of strings of length m_ii, the number of 
     * reactions, containing the
     * string expressions for each reaction
     * (e.g., reactants <=> product1 + product2)
     */
    std::vector<std::string> m_rxneqn;

    /**
     * Pointer to the class containing all of the temporary arrays
     * that are used during processing. This is malloced once for
     * each class instance.
     */
    SolidKineticsData* m_kdata;

    /**
     * Vector of the current activity concentrations.
     */
    vector_fp m_actConc;

    /**
     * Standard state chemical potentials of the species. This is
     * a temporary vector used in the evalulation of the equilibrium
     * constant.
     */
    vector_fp m_grt;

    /**
     * Boolean that is true if the standard concentration is a constant
     * and is the same for all species in the phase. The value is set up in
     * constructor based on the value of the eosType(). 
     * default value = false.
     */
    bool logC0AllTheSameConstant;

    /**
     * Boolean that is true if the standard concentration is a variable
     * of the temperature or pressure. The value is set up in the
     * constructor based on the value of the eosType(). If true, this
     * will trigger a recalculation of the standard concentrations in the
     * equilibrium constant routines and a reevaluation of the product
     * of the powers of the standard concentration wrt the stoichiometric
     * coefficient in the evaluation of Kc.
     * If false, the powers value is calculated once during initialization
     * of this object and reused multiple times.
     * default value = true;
     */
    bool logC0ProdVariable;

  private:

  
    std::vector<double> m_logProdC0;
    std::vector<double> m_logC0_vector;
    /**
     * boolean indicating that the class is ready for processing of 
     * reaction rates.
     */
    bool m_finalized;

    /**
     * When the standard concentrations are constants and equal
     * for all species, this entry will contain the log of the
     * common standard concentration.
     */
    doublevalue m_logC0_scalar;
  };
}

#endif
