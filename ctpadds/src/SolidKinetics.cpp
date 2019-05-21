/**
 *  @file SolidKinetics.cpp 
 *
 * Homogeneous kinetics in a condensed phase.
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

#include "SolidKinetics.h"
#include "zuzax/thermo/IdealSolidSolnPhase.h"
#include "zuzax/kinetics/ReactionData.h"

#include <iostream>
#include <fstream>

using namespace std;

namespace Zuzax
{
//===========================================================================================================
SolidKinetics::SolidKinetics() :
    BulkKinetics(),
    m_nirrev(0),
    m_nrev(0),
    m_kdata(0),
    logC0AllTheSameConstant(false),
    logC0ProdVariable(true),
    m_finalized(false),
    m_logC0_scalar(0.0)
{
    m_kdata = new SolidKineticsData;  
}
//=========================================================================================================
/*
 *
 * SolidKinetics():
 *
 * Constructor for SolidKinetics, starts with an empty reaction mechanism.
 * However, the thermo parameter is required. It makes no sense to 
 * construct
 * a reaction mechanism for a phase that doesn't have a ThermoPhase object
 * associated with it.
 */    
SolidKinetics::
SolidKinetics(thermo_t* thermo_ptr) :
    BulkKinetics(),
    m_nirrev(0),
    m_nrev(0),
    logC0AllTheSameConstant(false),
    logC0ProdVariable(true),
    m_finalized(false),
    m_logC0_scalar(0.0)
{
    m_kdata = new SolidKineticsData;

    if (!thermo_ptr) {
	throw ZuzaxError("SolidKinetics Constructor",
			   "Must supply a valid ThermoPhase object");
    }
    /*
     * Add the phase to the phase list kept in the Kinetics object.
     */
    addPhase(*thermo_ptr);
    /**
     * query ThermoPhase object for the equation of state type.
     * Fill in the flags logC0ProdVariable and   logC0AllTheSameConstant 
     * based on its value.
     */
    int eosFlag = thermo().eosType();
    switch (eosFlag) {
    case cIdealSolidSolnPhase0:
	logC0AllTheSameConstant = true;
	logC0ProdVariable = false;
	break;
    case cIdealSolidSolnPhase1:
	logC0ProdVariable = true;
	logC0AllTheSameConstant = false;
	break;
    case cIdealSolidSolnPhase2:
	logC0ProdVariable = false;
	logC0AllTheSameConstant = true;
	break;
    }
}

//====================================================================================================
SolidKinetics::SolidKinetics(const SolidKinetics& right) :
        BulkKinetics(),
        m_nirrev(0),
        m_nrev(0),
	logC0AllTheSameConstant(false),
	logC0ProdVariable(true),
	m_finalized(false),
        m_logC0_scalar(0.0)
{
    /*
     * Call the assignment operator
     */
    operator=(right);
}


SolidKinetics& SolidKinetics::operator=(const SolidKinetics& right) 
{
     /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    BulkKinetics::operator=(right);
    

    //  m_rates = right.m_rates;
    //m_index = right.m_index;
    // m_reactantStoich = right.m_reactantStoich;
    // m_revProductStoich = right.m_revProductStoich;
    // m_irrevProductStoich = right.m_irrevProductStoich;
    m_nirrev = right.m_nirrev;
    m_nrev = right. m_nrev;
    // m_rgroups = right. m_rgroups;
    //m_rrxn = right.m_rrxn;
    m_prxn = right.m_prxn;
    m_dn = right.m_dn;
    // m_revindex = right.m_revindex;
    // m_irrevindex = m_irrevindex;
    m_rstoich = right.m_rstoich;
    m_pstoich = right.m_pstoich;
    m_rxneqn = right.m_rxneqn;
    m_kdata = right.m_kdata;
    m_actConc = right.m_actConc;
    m_grt = right.m_grt;
    logC0AllTheSameConstant= right.logC0AllTheSameConstant;
    logC0ProdVariable	 = right.logC0ProdVariable;
    m_logProdC0        = right.m_logProdC0;
    m_logC0_vector = right.m_logC0_vector;
    m_finalized = right.m_finalized;
    m_logC0_scalar   = right.m_logC0_scalar;

    return *this;
}


/*************************************************************************
 *
 * ~SolidKinetics()
 *
 * Destructor for the SolidKinetics class
 */
SolidKinetics::~SolidKinetics() {
    delete m_kdata;
}

Kinetics* SolidKinetics::duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const
{
    SolidKinetics* gK = new SolidKinetics(*this);
    gK->assignShallowPointers(tpVector);
    return gK;
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
void SolidKinetics::getFwdRatesOfProgress(doublevalue* fwdROP) { 
    updateROP(); 
    copy(m_ropf.begin(), m_ropf.end(), fwdROP);
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
void SolidKinetics::getRevRatesOfProgress(doublevalue* revROP) { 
    updateROP(); 
    copy(m_ropr.begin(), m_ropr.end(), revROP);
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
void SolidKinetics::getNetRatesOfProgress(doublevalue* netROP) { 
    updateROP(); 
    copy(m_ropnet.begin(), m_ropnet.end(), netROP);
}

    /***************************************************************************
     *
     * getNetProductionRates():
     *
     * Species net production rates [kmol/m^3 s^1]. Return the species
     * net production rates (creation - destruction) in array
     * net, which must be dimensioned at least as large as the
     * total number of species.
     */
void SolidKinetics::
getNetProductionRates(doublevalue* net) {
	/*
	 * We do the work here. We calculate reactions rates of progress and
	 * store them in the vector m_kdata->m_ropnet
	 */
    updateROP();
	/*
	 * Zero out the return vector
	 */
    fill(net, net + m_NumKinSpecies, 0.0);
	/*
	 * Go call the stoichiometry managers to obtain the production
	 * rates of the species.
	 */
    m_revProductStoich.incrementSpecies(DATA_PTR(m_ropnet), net);
    m_irrevProductStoich.incrementSpecies(DATA_PTR(m_ropnet), net);
    m_reactantStoich.decrementSpecies(DATA_PTR(m_ropnet), net);
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
    void SolidKinetics::getCreationRates(doublevalue* cdot) {
	/*
	 * Update the rates of progress of the reactions
	 */
	updateROP();
	fill(cdot, cdot + m_NumKinSpecies, 0.0);
	m_revProductStoich.incrementSpecies(DATA_PTR(m_ropf), cdot);
	m_irrevProductStoich.incrementSpecies(DATA_PTR(m_ropf), cdot);
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
    void SolidKinetics::getDestructionRates(doublevalue* ddot) {
	/*
	 * Update the rates of progress of the reactions
	 */
	updateROP();
	/*
	 * zero the matrix
	 */
	std::fill(ddot, ddot + m_NumKinSpecies, 0.0);
	m_reactantStoich.incrementSpecies(DATA_PTR(m_ropf), ddot);
    }

    /************************************************************************
     *
     * update_C():
     *
     * Update concentration-dependent portions of reaction rates and
     * falloff functions.
     */
    void SolidKinetics::
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
    bool SolidKinetics::isReversible(size_t i) {
	if (find(m_revindex.begin(), m_revindex.end(), i) 
	    < m_revindex.end()) return true;
	else return false;
    }
//===============================================================================================
/**********************************************************************
 *
 * _update_rates_T():
 *
 * Update the temperature dependent portions of the rate constants
 * Objects which are updated:
 *  m_rates - the rate constants. 
 *  m_kdata->m_rkc[] Inverse of the equilibrium constants
 */
void SolidKinetics::
_update_rates_T() {
    doublevalue T = thermo().temperature();
    /*
     * Calculate logC0 if all the same.
     */
    if (logC0AllTheSameConstant) {
	m_logC0_scalar = log(thermo().standardConcentration(0));
    }
    
    doublevalue logT = log(T);
    /**
     * Update the forward rate constants. We only update those
     * rate constants which have a temperature dependence in 
     * this step.
     */
    // m_rates.update(T, logT, DATA_PTR(m_kdata->m_rfn));
    m_rates.update(T, logT, &m_rfn[0]);
    /*
     * Store the temperature at which the rate constants were
     * last evaluated.
     */
    m_temp = T;
    /*
     * Update the equilibrium constants for reversible reactions only
     */
    updateKc();
    
    m_ROP_ok = false;
}
//==================================================================================================================================
    /***********************************************************************
     *
     * _update_rates_C():
     *
     * Update properties that depend on concentrations. Currently this
     * only involves the calculation of the generalized concentrations
     */         
    void SolidKinetics::
    _update_rates_C() {
        thermo().getActivityConcentrations(DATA_PTR(m_actConc));
        m_ROP_ok = false;
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
    void SolidKinetics::updateKc() {
        int i, irxn;

	/*
	 *  Get the standard state chemical potentials of the species.
         *  This is the array of chemical potentials at unit activity 
	 *  We define these here as the chemical potentials of the pure
	 *  species at the temperature and pressure of the solution.
	 */
	thermo().getStandardChemPotentials(DATA_PTR(m_grt));
 
	/*
	 * Use the stoichiometric manager to find delta G^0 for each
	 * reaction.
	 * First, zero out m_rkc[]. Then fill it in for each reaction
	 * whether reversible or not.
	 */
	vector_fp& m_rkc = m_rkcn;
	fill(m_rkc.begin(), m_rkc.end(), 0.0);
        m_reactantStoich.decrementReactions(DATA_PTR(m_grt), DATA_PTR(m_rkc));
        m_revProductStoich.incrementReactions(DATA_PTR(m_grt),DATA_PTR(m_rkc));
	
	/*
	 * If required, recalculate the log of the products of the
	 * standard concentrations.
	 */
	vector_fp& logProdC0 = m_logProdC0;
	if (logC0ProdVariable) {
	    vector_fp& logC0_vector = m_logC0_vector;
	  for (size_t i = 0; i < m_NumKinSpecies; ++i) {
	    m_logC0_vector[i] = thermo().logStandardConc(i);
	  }
	  fill(logProdC0.begin(), logProdC0.end(), 0.0);
	  m_reactantStoich.decrementReactions(DATA_PTR(logC0_vector), 
					      DATA_PTR(logProdC0)); 
	  m_revProductStoich.incrementReactions(DATA_PTR(logC0_vector), 
						DATA_PTR(logProdC0));
	}
 
        doublevalue rrt = 1.0/(GasConstant * thermo().temperature());
	/*
	 * Branch on whether the standard concentration is the same 
	 * constant for all species or not.
	 */
	if (logC0AllTheSameConstant) {
	  doublevalue logc0 = m_logC0_scalar;
	  for (i = 0; i < m_nrev; i++) {
            irxn = m_revindex[i];
            m_rkc[irxn] = exp(m_rkc[irxn]*rrt - m_dn[irxn]*logc0);
	  }
	} else {
	  for (i = 0; i < m_nrev; i++) {
            irxn = m_revindex[i];
            m_rkc[irxn] = exp(m_rkc[irxn]*rrt - logProdC0[irxn]);
	  }
	}

	/*
	 * Zero out the equilibrium constant for irreversible reactions.
	 */
        for(i = 0; i != m_nirrev; ++i) {
            m_rkc[ m_irrev[i] ] = 0.0;
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
    void SolidKinetics::
    getEquilibriumConstants(doublevalue* kc) {
        size_t i;
	
	/*
	 *  Get the standard state chemical potentials of the species.
         *  This is the array of chemical potentials at unit activity 
	 *  We define these here as the chemical potentials of the pure
	 *  species at the temperature and pressure of the solution.
	 */
        thermo().getStandardChemPotentials(DATA_PTR(m_grt));

	/*
	 * Use the stoichiometric manager to find delta G^0 for each
	 * reaction.
	 * First, zero out m_rkc[]. Then fill it in for each reaction
	 * whether reversible or not.
	 */
	vector_fp& rkc = m_rkcn;
        fill(rkc.begin(), rkc.end(), 0.0);
        m_reactantStoich.decrementReactions(DATA_PTR(m_grt), DATA_PTR(rkc)); 
        m_revProductStoich.incrementReactions(DATA_PTR(m_grt),DATA_PTR(rkc));
        m_irrevProductStoich.incrementReactions(DATA_PTR(m_grt),DATA_PTR(rkc));
 	
	/*
	 * If required, recalculate the log of the products of the
	 * standard concentrations.
	 */
	vector_fp& logProdC0 = m_logProdC0;
	if (logC0ProdVariable) {
	  vector_fp& logC0_vector = m_logC0_vector;
	  for (size_t i = 0; i < m_NumKinSpecies; ++i) {
	    logC0_vector[i] = thermo().logStandardConc(i);
	  }
	  fill(logProdC0.begin(), logProdC0.end(), 0.0);
	  m_reactantStoich.decrementReactions(DATA_PTR(logC0_vector),
					      DATA_PTR(logProdC0)); 
	  m_revProductStoich.incrementReactions(DATA_PTR(logC0_vector), 
						DATA_PTR(logProdC0));
	  m_irrevProductStoich.incrementReactions(DATA_PTR(logC0_vector),
						  DATA_PTR(logProdC0));
	}

 	/*
	 * Branch on whether the standard concentration is the same 
	 * constant for all species or not.
	 */
	doublevalue rrt = 1.0/(GasConstant * thermo().temperature());
	if (logC0AllTheSameConstant) {
	  const doublevalue logC0 = m_logC0_scalar;
	  for (size_t i = 0; i < m_ii; i++) {
            kc[i] = exp(- rkc[i]*rrt + m_dn[i]*logC0);
	  }
	} else {
	  for (i = 0; i < m_ii; i++) {
            kc[i] = exp( - rkc[i]*rrt + logProdC0[i]);
	  }
	}
    }

    /***********************************************************************
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
    void SolidKinetics::getDeltaGibbs(doublevalue* deltaG) {
	/*
	 * Get the chemical potentials of the species in the 
	 * solid solution.Note this works as SolidKinetics is 
	 * limited to one phase.
	 */
	thermo().getChemPotentials(DATA_PTR(m_grt));

	/*
	 * Use the stoichiometric manager to find deltaG for each
	 * reaction.
	 * First, zero out deltaG. Then fill it in for each reaction
	 * whether reversible or not.
	 */
	
	for (size_t i = 0; i < m_ii; i++) {
	  deltaG[i] = 0.0;
	}
        m_reactantStoich.decrementReactions(DATA_PTR(m_grt), deltaG); 
        m_revProductStoich.incrementReactions(DATA_PTR(m_grt), deltaG);
        m_irrevProductStoich.incrementReactions(DATA_PTR(m_grt), deltaG);
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
    void SolidKinetics::getDeltaEnthalpy(doublevalue* deltaH) {
	/*
	 * Get the partial molar enthalpy of all species in the 
	 * solid solution. Note this works as SolidKinetics is 
	 * limited to one phase.
	 */
	thermo().getPartialMolarEnthalpies(DATA_PTR(m_grt));

	/*
	 * Use the stoichiometric manager to find deltaH for each
	 * reaction.
	 * First, zero out deltaH. Then fill it in for each reaction
	 * whether reversible or not.
	 */
	for (size_t i = 0; i < m_ii; i++) {
	  deltaH[i] = 0.0;
	}
        m_reactantStoich.decrementReactions(DATA_PTR(m_grt), deltaH); 
        m_revProductStoich.incrementReactions(DATA_PTR(m_grt), deltaH);
        m_irrevProductStoich.incrementReactions(DATA_PTR(m_grt), deltaH);
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
    void SolidKinetics::getDeltaEntropy( doublevalue* deltaS) {
	/*
	 * Get the partial molar entropy of all species in the
	 * solid solution. Note this works as SolidKinetics is 
	 * limited to one phase.
	 */
	thermo().getPartialMolarEntropies(DATA_PTR(m_grt));

	/*
	 * Use the stoichiometric manager to find deltaS for each
	 * reaction.
	 * First, zero out deltaS. Then fill it in for each reaction
	 * whether reversible or not.
	 */
	for (size_t i = 0; i < m_ii; i++) {
	  deltaS[i] = 0.0;
	}
        m_reactantStoich.decrementReactions(DATA_PTR(m_grt), deltaS); 
        m_revProductStoich.incrementReactions(DATA_PTR(m_grt), deltaS);
        m_irrevProductStoich.incrementReactions(DATA_PTR(m_grt), deltaS);
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
    void SolidKinetics::getDeltaSSGibbs(doublevalue* deltaG) {
	/*
	 *  Get the standard state chemical potentials of the species.
         *  This is the array of chemical potentials at unit activity 
	 *  We define these here as the chemical potentials of the pure
	 *  species at the temperature and pressure of the solution.
	 */
        thermo().getStandardChemPotentials(DATA_PTR(m_grt));
	/*
	 * Use the stoichiometric manager to find deltaG for each
	 * reaction.
	 * First, zero out deltaG. Then fill it in for each reaction
	 * whether reversible or not.
	 */
	for (size_t i = 0; i < m_ii; i++) {
	  deltaG[i] = 0.0;
	}
        m_reactantStoich.decrementReactions(DATA_PTR(m_grt), deltaG); 
        m_revProductStoich.incrementReactions(DATA_PTR(m_grt), deltaG);
        m_irrevProductStoich.incrementReactions(DATA_PTR(m_grt), deltaG);
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
    void SolidKinetics::getDeltaSSEnthalpy(doublevalue* deltaH) {
	/*
	 *  Get the standard state enthalpies of the species.
         *  This is the array of chemical potentials at unit activity 
	 *  We define these here as the enthalpies of the pure
	 *  species at the temperature and pressure of the solution.
	 */
	thermo().getEnthalpy_RT(DATA_PTR(m_grt));
	doublevalue RT = thermo().temperature() * GasConstant;
	for (size_t k = 0; k < m_NumKinSpecies; k++) {
	  m_grt[k] *= RT;
	}
	/*
	 * Use the stoichiometric manager to find deltaH for each
	 * reaction.
	 * First, zero out deltaH. Then fill it in for each reaction
	 * whether reversible or not.
	 */
	for (size_t i = 0; i < m_ii; i++) {
	  deltaH[i] = 0.0;
	}
        m_reactantStoich.decrementReactions(DATA_PTR(m_grt), deltaH); 
        m_revProductStoich.incrementReactions(DATA_PTR(m_grt), deltaH);
        m_irrevProductStoich.incrementReactions(DATA_PTR(m_grt), deltaH);
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
    void SolidKinetics::getDeltaSSEntropy(doublevalue* deltaS) {
	/*
	 *  Get the standard state entropy of the species.
	 *  We define these here as the entropies of the pure
	 *  species at the temperature and pressure of the solution.
	 */
	thermo().getEntropy_R(DATA_PTR(m_grt));
	doublevalue R = GasConstant;
	for (size_t k = 0; k < m_NumKinSpecies; k++) {
	  m_grt[k] *= R;
	}
	/*
	 * Use the stoichiometric manager to find deltaS for each
	 * reaction.
	 * First, zero out deltaS. Then fill it in for each reaction
	 * whether reversible or not.
	 */
	for (size_t i = 0; i < m_ii; i++) {
	  deltaS[i] = 0.0;
	}
        m_reactantStoich.decrementReactions(DATA_PTR(m_grt), deltaS); 
        m_revProductStoich.incrementReactions(DATA_PTR(m_grt), deltaS);
        m_irrevProductStoich.incrementReactions(DATA_PTR(m_grt), deltaS);
    }

    /*********************************************************************
     *
     * updateROP():
     *
     * Update the rate of progress for the reactions.
     * This key routine makes sure that the rate of progress vectors
     * located in the solid kinetics data class are up to date.
     */
void SolidKinetics::updateROP() {

        _update_rates_T();
        _update_rates_C();

	/*
	 * If the ROP's are up to date return without doing any
	 * calculations -> The update rate coefficients above will
	 * check to see if the temperature and concentrations have
	 * changed. 
	 */
        if (m_ROP_ok) return;

        /*
	 * copy the forward rate coefficients, m_rfn, into forward
	 * rate of progress vector, ropf, to start the process.
	 */
	const vector_fp& rf = m_rfn;
	vector_fp& ropf = m_ropf;
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
	vector_fp& ropr = m_ropr;
        copy(ropf.begin(), ropf.end(), ropr.begin());
        
        // for reverse rates computed from thermochemistry, multiply
        // the forward rates copied into m_ropr by the reciprocals of
        // the equilibrium constants
	const vector_fp& m_rkc = m_rkcn;
        multiply_each(ropr.begin(), ropr.end(), m_rkc.begin());

        /*
	 * Multiply the rate coefficients by the activity concentrations
	 * to generate the forward rate of progress of the reaction
	 */
        m_reactantStoich.multiply(DATA_PTR(m_actConc), DATA_PTR(ropf)); 

        // for reversible reactions, multiply ropr by the
	// activity concentration products
        m_revProductStoich.multiply(DATA_PTR(m_actConc), DATA_PTR(ropr));

	/**
	 * Calculate the net rate of progress for a reaction from the
	 * difference of the forward and the reverse.
	 */
	vector_fp& ropnet = m_ropnet;
	//printf("Size of ropnet = %d\n", ropnet.size());
	//printf("Size of ropf = %d\n", ropf.size());
	//printf("Size of ropr = %d\n", ropr.size());
        //printf("Size of m_ii = %d\n", m_ii);

        for (size_t j = 0; j < m_ii; ++j) {
            ropnet[j] = ropf[j] - ropr[j];
        }

	/*
	 * signal that the ROP's are up to date
	 */
        m_ROP_ok = true;
    }

    /*********************************************************************
     *
     * getFwdRateConstants():
     *
     * Update the rate of progress for the reactions.
     * This key routine makes sure that the rate of progress vectors
     * located in the solid kinetics data class are up to date.
     */
    void SolidKinetics::
    getFwdRateConstants(doublevalue *kfwd) {
        _update_rates_T();
	_update_rates_C();
	const vector_fp& rf = m_rfn;
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
    void SolidKinetics::
    getRevRateConstants(doublevalue *krev, bool doIrreversible) {
	_update_rates_T();
	_update_rates_C();
	const vector_fp& rf = m_rfn;
	if (doIrreversible) {
	  doublevalue *tmpKc = DATA_PTR(m_ropnet);
	  getEquilibriumConstants(tmpKc);
	  for (size_t i = 0; i < m_ii; i++) {
	    krev[i] = rf[i] * m_perturb[i] / tmpKc[i];
	  }
	} else {
	  /*
	   * m_rkc[] is zero for irreversibly reactions
	   */
	  const vector_fp& m_rkc = m_rkcn;
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
//============================================================================================
void SolidKinetics::
addReaction(ReactionData& r) {
    //int n, m;
    // double ns;
    if (r.reactionType == ELEMENTARY_RXN) {
        addElementaryReaction(r);
    } else {
	throw ZuzaxError("addReaction", "Unknown reaction type");
    }

    BulkKinetics::addReaction(r);

    /*
     * Check m_finalized
     */
    if (m_finalized) {
	throw ZuzaxError("addReaction():", 
			   "Mechanism already finalized");
    }
	
    /*
     * Check to see that the elements are balanced in the reaction
     */
    checkRxnElementBalance(*this, r);

  
    /*
     * Add the reactants and products for  m_ropnet;the current reaction
     * to the various stoichiometric coefficient arrays.
     */
    // installReagents(r);
    m_kdata->m_ropf.push_back(0.0);
    m_kdata->m_ropr.push_back(0.0);
    m_kdata->m_ropnet.push_back(0.0);
    m_kdata->m_rkcn.push_back(0.0);

    m_logProdC0.push_back(0.0);
    /*
     * Go get the total number of reactions in the mechanism
     * added to date (m_ii)
     */
    //int rnum = nReactions();

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
    // int nr = r.reactants.size();
    //  for (n = 0; n < nr; n++) {
	/*
	 * Look up the stoichiometric coefficient for the
	 * species in the current reaction
	 */
//	ns = r.rstoich[n];
	/*
	 * Add the current species to the stoichiometric
	 * coefficient reactants. (Q: m_rrxn the right size?)
	 */
//	if (ns != 0.0) m_rrxn[r.reactants[n]][rnum] += ns;
	/*
	 * Create a number integer vector entries equal to the
	 * stoichiometric coefficient for the current reactant
	 */
//	for (m = 0; m < ns; m++) {
    //rk.push_back(r.reactants[n]);
    //}
    //}
    /*
     * Add an entry for m_reactants, which is a member
     * of the Kinetics class.
     */
    // m_reactants.push_back(rk);
    /*
     * Add an entry for the reactants into the stoichiometric
     * manager.
     */
    // m_reactantStoich.add((size_t) nReactions(), rk);
    
    /*
     * NEXT, WE WILL TAKE CARE OF STORAGE FOR THE PRODUCTS
     */
    // std::vector<size_t> pk;
    /*
     * look up the number of products in the reaction, np.
     * and then loop over them
     */	
    // int np = r.products.size();
    // for (n = 0; n < np; n++) {
	/*
	 * Look up the stoichiometric coefficient for the
	 * current species product in the current reaction
	 */
//	ns = r.pstoich[n];
	/*
	 * Add the current species to the stoichiometric
	 * coefficient products. (Q: m_prxn the right size?)
	 */
//	if (ns != 0) m_prxn[r.products[n]][rnum] += ns;
	/*
	 * Create a number of integer vector entries equal to the
	 * stoichiometric coefficient for the current product
	 */
//	for (m = 0; m < ns; m++) {
    // pk.push_back(r.products[n]);
//	}
    // }
    /*
     * Add an entry for m_reactants, which is a member
     * of the Kinetics class.
     */
    //  m_products.push_back(pk);
    
    /*
     * m_dn - fill in the entry for the net change of moles
     *        due to this reaction. 
     */
    // m_dn.push_back(pk.size() - rk.size());
    
    /*
     * Branch on whether this reaction is reversible
     * -> fill in entries in the stoichiometric manager
     *    and in the listing of reversable and irreversible
     *    reactions.
     */
    
    if (r.reversible) {
//	m_revProductStoich.add(m_ii, pk);
//	m_revindex.push_back(m_ii);
	m_nrev++;
    }
    else {
//	m_irrevProductStoich.add(m_ii, pk);          
//	m_irrevindex.push_back(m_ii );
	m_nirrev++;
    }        
    
    /* 
     * Save the reaction and product groups, which are
     * part of the ReactionData class, in this class.
     * They aren't used for anything but reaction path
     * analysis.
     */
    //installGroups(m_ii, r.rgroups, r.pgroups);
    /*
     * Increase the internal number of reactions, m_ii, by one.
     * increase the size of m_perturb by one as well.
     */
    // incrementRxnCount();
    /*
     * Store the string expression for the current reaction in
     * the vector of strings.
     */
    //  m_rxneqn.push_back(r.equation);
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
//void SolidKinetics::
//addElementaryReaction(const ReactionData& r) {

//   m_rates.install(m_ii, r);

    // install rate coeff calculator
    // iloc = m_rates.install(reactionNumber(),
    //		       r.rateCoeffType, r.rateCoeffParameters.size(), 
    //		       DATA_PTR(r.rateCoeffParameters) );

    /*
     * add constant term to rate coeff value vector.
     * Note, for constant reaction rates coefficients, this is the
     * one and only place where the rate coefficient is installed
     * into the rate coefficient array.
     */
    //    m_kdata->m_rfn.push_back(r.rateCoeffParameters[0]);                

    /*
     * Fill in the m_index entry for this reaction. This is an
     * information map, keyed on the reaction number,
     * describing the reaction type (always ELEMENTARY_RXN for this
     * class) and the index of the Reaction Rate Coefficient Calculator. 
     */
    //     registerReaction(reactionNumber(), ELEMENTARY_RXN, iloc);


    // install rate coeff calculator
//  vector_fp rp = r.rateCoeffParameters;
    // store activation energy
    //m_E.push_back(r.rateCoeffParameters[2]);

//  if (r.beta > 0.0) {
        //m_has_electrochem_rxns = true;
        //m_beta.push_back(r.beta);
        //m_ctrxn.push_back(reactionNumber());
	//        m_ctrxn_ecdf.push_back(0);
//  }

    // add constant term to rate coeff value vector
//   m_kdata->m_rfn.push_back(r.rateCoeffParameters[0]);
    //registerReaction(reactionNumber(), ELEMENTARY_RXN, iloc);

//}
   
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


    /*********************************************************************
     *
     * installGroups()
     *
     * Save the reaction and product groups, which are
     * part of the ReactionData class, in this class.
     * They aren't used for anything but reaction path
     * analysis, as far as I can figure out.
     */
//void SolidKinetics::installGroups(size_t irxn, 
//				  const vector<grouplist_t>& r, const vector<grouplist_t>& p) {
//    if (!r.empty()) {
//	m_rgroups[m_ii] = r;
//	m_pgroups[m_ii] = p;
//    }
//}
    
    /*********************************************************************
     *
     * registerReaction():
     *
     * Creates an information map about the ith reaction.
     * For the ith reaction,
     * the reaction type and the iloc location in the m_rate array
     * are saved as a pair object in the internal variable m_index.
     */
    //void SolidKinetics::registerReaction(size_t rxnNumber, int type, int loc) {
//	m_index[rxnNumber] = pair<int, int>(type, loc);
//  }

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
    void SolidKinetics::init() { 
        BulkKinetics::init();
        m_NumKinSpecies = thermo().nSpecies();
	if (m_NumKinSpecies <= 0) {
	  throw ZuzaxError("SolidKinetics::init",
			     "m_NumKinSpecies is zero or less");
	}
        //m_rrxn.resize(m_NumKinSpecies);
        //m_prxn.resize(m_NumKinSpecies);
        m_actConc.resize(m_NumKinSpecies);
        m_grt.resize(m_NumKinSpecies);
	//m_kdata->m_logC0_vector.resize(m_NumKinSpecies);

	m_logC0_vector.resize(m_NumKinSpecies);

    }


    /*********************************************************************
     *
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
     *              defaults to "" in the .h file.
     */
    void SolidKinetics::
    importMechanism(std::string inputFile, std::string id) {
	string path = findInputFile(inputFile);
	std::ifstream fin(path.c_str());
	if (!fin) {
	  throw ZuzaxError("SolidKinetics::importMechanism",
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
    void SolidKinetics::
    importMechanism(XML_Node& phaseNode, string id) {
	string idp = phaseNode.id();
	if (id.size() > 0) {
	  if (idp != id) {
	    throw ZuzaxError("initThermo", 
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
	  throw ZuzaxError("SolidKinetics::importMechanism",
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
	  throw ZuzaxError("SolidKinetics::importMechanism", 
			     "phaseArray element not allowed" );
        }

	/*
	 * Look for a child of the xml element phase called
	 * "kinetics". It has an attribute name "model".
	 * Store the value of that attribute in the variable kintype
	 */
        string kintype = phaseNode.child("kinetics")["model"];
	if (kintype != "SolidKinetics") {
	  throw ZuzaxError("SolidKinetics::importMechanism", 
			     "kinetics model must have model "
			     "attribute SolidKinetics");
	  
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
	      throw ZuzaxError("SolidKinetics::importMechanism", 
				 "kinetics submodel must have submodel "
				 "attribute \"unity\"");
	    }
	    break;
        case cIdealSolidSolnPhase1:
	    if (kinsubtype != "molar_volume") {
	      throw ZuzaxError("SolidKinetics::importMechanism", 
				 "kinetics submodel must have submodel "
				 "attribute \"molar volume\"");
	    }	    
	    break;
	case cIdealSolidSolnPhase2:
	    if (kinsubtype != "solvent_volume") {
	      throw ZuzaxError("SolidKinetics::importMechanism", 
				 "kinetics submodel must have submodel "
				 "attribute \"solvent volume\"");
	    }
	    break;
	default:
	    throw ZuzaxError("SolidKinetics::importMechanism", 
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
	  throw ZuzaxError("SolidKinetics::importMechanism",
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
void SolidKinetics::finalize() {
    doublevalue rstoiV, pstoiV;
    if (!m_finalized) {
	size_t i;
        
	for (i = 0; i < m_ii; i++) {
	    std::map<size_t, doublevalue>& rstoichIrxn = m_rstoich[i];
	    std::map<size_t, doublevalue>& pstoichIrxn = m_pstoich[i];
	    for (size_t k = 0; k < m_NumKinSpecies; k++) {
                rstoiV = reactantStoichCoeff(k, i);
                if (rstoiV != 0.0) {
		    rstoichIrxn[k] = rstoiV;
                }
                pstoiV = productStoichCoeff(k, i);
                if (pstoiV != 0.0) {
		    pstoichIrxn[k] = pstoiV;
                }
	    }
	}
	
	/*
	 * Provide initial values for m_logC0_vector, the vector
	 * of standard concentrations for the species. This may
	 * or may not be overwritten as the conditions change
	 * depending upon the formulation of the ThermoPhase
	 * object.
	 */
	vector_fp& logC0_vector = m_logC0_vector;
	for (size_t i = 0; i < m_NumKinSpecies; ++i) {
	    logC0_vector[i] = thermo().logStandardConc(i);
	}
	m_logC0_scalar =  logC0_vector[0];
	
	/*
	 * Provide initial values for m_logProdC0, the vector 
	 * (over reactions) of the log of the products of 
	 * powers of the standard concentrations
	 * wrt the stoichiometric coefficients in the reaction.
	 * This may or may not be recalculated as conditions 
	 * change. However, an initial calculation is done here.
	 */
	vector_fp& logProdC0 = m_logProdC0;
	fill(logProdC0.begin(), logProdC0.end(), 0.0);
	m_reactantStoich.decrementReactions(DATA_PTR(logC0_vector), 
					    DATA_PTR(logProdC0)); 
	m_revProductStoich.incrementReactions(DATA_PTR(logC0_vector), 
					      DATA_PTR(logProdC0));
	m_finalized = true;
    }
}
    /********************************************************************
     *
     * ready():
     *
     *  Indicate whether the reaction mechanism is ready for use.
     */
bool SolidKinetics::ready() const {
    return (m_finalized);
}

    /*********************************************************************/
}








