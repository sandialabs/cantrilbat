/**
 * @file ReactingSurDomain.cpp
 *  Definitions for the ElectrodeKinetics object that does handles interactions with the PhaseList object
 *  (see \ref ExtendedPhaseGroups and class \link Cantera::ReactingSurDomain ReactingSurDomain\endlink).
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "ReactingSurDomain.h"
#include "PhaseList.h"
//#include "RxnMolChange.h"
#include "cantera/kinetics.h"
#include "Electrode_input.h"

#include "Electrode_Factory.h"

#include "mdp_allo.h"

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
    ElectrodeKinetics(),
    numPhases_(0),
    xmlList(0),
    kinOrder(0),
    PLtoKinPhaseIndex_(0),
    PLtoKinSpeciesIndex_(0),
    KintoPLSpeciesIndex_(0),
    iphaseKin_(-1),
    tpList_IDs_(0),
    tplRead(0),
    m_DoSurfKinetics(false),
    speciesProductionRates_(0),
    limitedROP_(0),
    speciesCreationRates_(0),
    speciesDestructionRates_(0),
    deltaGRxn_Before_(0),
    deltaGRxnOCV_Before_(0),
    deltaHRxn_Before_(0),
    deltaSRxn_Before_(0),
    m_pl(0),
    ocv_ptr_(0),
    OCVmodel_(0),
    kReplacedSpeciesRS_(-1),
    m_Enthalpies_rspec(0),
    m_Enthalpies_Before_rspec(0),
    m_Entropies_rspec(0),
    m_Entropies_Before_rspec(0),
    m_GibbsOCV_rspec(0),
    m_Gibbs_Before_rspec(0),
    deltaG_species_(0.0),
    deltaS_species_(0.0),
    deltaH_species_(0.0)
{
}
//====================================================================================================================
// Copy Constructor for the %Kinetics object.
/*
 * Currently, this is not fully implemented. If called it will
 * throw an exception.
 */
ReactingSurDomain::ReactingSurDomain(const ReactingSurDomain& right) :
    ElectrodeKinetics(),
    numPhases_(0),
    xmlList(0),
    kinOrder(0),
    PLtoKinPhaseIndex_(0),
    PLtoKinSpeciesIndex_(0),
    KintoPLSpeciesIndex_(0),
    iphaseKin_(-1),
    tpList_IDs_(0),
    tplRead(0),
    m_DoSurfKinetics(false),
    speciesProductionRates_(0),
    limitedROP_(0),
    speciesCreationRates_(0),
    speciesDestructionRates_(0),
    deltaGRxn_Before_(0),
    deltaGRxnOCV_Before_(0),
    deltaHRxn_Before_(0),
    deltaSRxn_Before_(0),
    m_pl(0),
    ocv_ptr_(0),
    OCVmodel_(0),
    kReplacedSpeciesRS_(-1),
    m_Enthalpies_rspec(0),
    m_Enthalpies_Before_rspec(0),
    m_Entropies_rspec(0),
    m_Entropies_Before_rspec(0),
    m_GibbsOCV_rspec(0),
    m_Gibbs_Before_rspec(0),
    deltaG_species_(0.0),
    deltaS_species_(0.0),
    deltaH_species_(0.0)
{
    /*
     * Call the assignment operator
     */
    operator=(right);
}
//==================================================================================================================================
ReactingSurDomain::ReactingSurDomain(Cantera::PhaseList* pl, int iskin) :
    ElectrodeKinetics(),
    numPhases_(0),
    kinOrder(pl->nPhases(), -1),
    PLtoKinPhaseIndex_(pl->nPhases(), -1),
    PLtoKinSpeciesIndex_(pl->nSpecies(), -1),
    KintoPLSpeciesIndex_(0),
    iphaseKin_(iskin + pl->nVolPhases()),
    tplRead(pl->nPhases(),0),
    m_DoSurfKinetics(true),
    speciesProductionRates_(0),
    limitedROP_(0),
    speciesCreationRates_(0),
    speciesDestructionRates_(0),
    deltaGRxn_Before_(0),
    deltaGRxnOCV_Before_(0),
    deltaHRxn_Before_(0),
    deltaSRxn_Before_(0),
    m_pl(pl),
    ocv_ptr_(0),
    OCVmodel_(0),
    kReplacedSpeciesRS_(-1),
    m_Enthalpies_rspec(0),
    m_Enthalpies_Before_rspec(0),
    m_Entropies_rspec(0),
    m_Entropies_Before_rspec(0),
    m_GibbsOCV_rspec(0),
    m_Gibbs_Before_rspec(0),
    deltaG_species_(0.0),
    deltaS_species_(0.0),
    deltaH_species_(0.0)
{
    bool ok = importFromPL(pl, iskin);
    if (!ok) {
	throw Electrode_Error("ReactingSurDomain::ReactingSurDomain()",
			      "import from Phase list failed");
    }
}
//==================================================================================================================================
// Assignment operator
/*
 *  This is NOT a virtual function.
 *
 * @param right    Reference to %Kinetics object to be copied into the
 *                 current one.
 */
ReactingSurDomain& ReactingSurDomain::operator=(const ReactingSurDomain& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }
    //
    //  Beware: The copy operation within ElectrodeKinetics leaves shallow pointers to the
    //          underlying ThermoPhase classes in place. They must be fixed up at the Electrode
    //          object level where we copy the ThermoPhase classes.
    //
    ElectrodeKinetics::operator=(right);

    numPhases_      = right.numPhases_;
    // Shallow copy of xmlList pointers -> beware
    xmlList         = right.xmlList;
    kinOrder        = right.kinOrder;
    PLtoKinPhaseIndex_ = right.PLtoKinPhaseIndex_;
    PLtoKinSpeciesIndex_ = right.PLtoKinSpeciesIndex_;
    KintoPLSpeciesIndex_ = right.KintoPLSpeciesIndex_;
    iphaseKin_      = right.iphaseKin_;
    tpList_IDs_     = right.tpList_IDs_;
    tplRead         = right.tplRead;
    m_DoSurfKinetics = right.m_DoSurfKinetics;
    speciesProductionRates_ = right.speciesProductionRates_;
    limitedROP_     = right.limitedROP_;
    speciesCreationRates_ = right.speciesCreationRates_;
    speciesDestructionRates_ = right.speciesDestructionRates_;

    deltaGRxn_Before_ = right.deltaGRxn_Before_;
    deltaGRxnOCV_Before_ = right.deltaGRxnOCV_Before_;
    deltaHRxn_Before_ = right.deltaHRxn_Before_;
    deltaSRxn_Before_ = right.deltaSRxn_Before_;
    //
    // Beware -  Shallow copy of m_pl pointer
    //
    m_pl            = right.m_pl;


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

    kReplacedSpeciesRS_ = right.kReplacedSpeciesRS_;

    m_Enthalpies_rspec = right.m_Enthalpies_rspec;
    m_Enthalpies_Before_rspec = right.m_Enthalpies_Before_rspec;
    m_Entropies_rspec = right.m_Entropies_rspec;
    m_Entropies_Before_rspec = right.m_Entropies_Before_rspec;
    m_GibbsOCV_rspec = right.m_GibbsOCV_rspec;
    m_Gibbs_Before_rspec = right.m_Gibbs_Before_rspec;

    deltaG_species_ = right.deltaG_species_;
    deltaS_species_ = right.deltaS_species_;
    deltaH_species_ = right.deltaH_species_;

    return *this;
}
//==================================================================================================================================
/*
 * Destructor for the ReactingSurDomain object.
 *
 * We must decide whether this object owns its own xml tree
 * structure.
 */
ReactingSurDomain::~ReactingSurDomain()
{
}
//==================================================================================================================================
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
//==================================================================================================================================
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
//==================================================================================================================================
void ReactingSurDomain::limitROP(const double* const n)
{
    // Apply smooth limiter to ROP
    /*
     *         Seems to basically satisfy the needs of turning interfacial kinetics into homogeneous kinetics.
     *         There is a basic exponential decay algorithm created out of the interfacial kinetics, which is
     *         not exponential decay.
     *
     *  HKM -> several problems. One is that the surface area doesn't enter into this limiter. Therefore, the
     *         limiter doesn't have the correct units, and will fail as surface area scales to different values
     *         than what it was set up to be. That way you could also talk about the algorithm setting a minimum "time
     *         constant" for phase removal.
     *
     *         Second is that rate turns off for phase_moles = 1.0E-20. There's no reason for that as this is straight
     *         exponential decay. The 1.0E-20 value has to be reconciled with other small numbers in the algorithm
     *         to see if there isn't a conflict.
     *
     *         Third, the algorithm only works for phase_moles > 0.0. Phase pop problems typically have times when
     *         phase_moles < 0.0. What happens in these cases. Also non-max principle algorithms will have problems.
     *
     *         Fourth, algorithm should maybe not be based on phase_moles, but individual moles. One really wants
     *         all moles to stay positive, not just the phase moles.
     */
    for (size_t j = 0; j < nReactions(); ++j) {
	double netRev = m_ropr[j] - m_ropf[j];
	double netFwd = -netRev;
	
	for (size_t p = 0; p < nPhases(); ++p) {
	    if (m_rxnPhaseIsProduct[j][p] && netRev > 0.0) {
		ThermoPhase* tp = m_thermo[p];
		int k = m_pl->getGlobalSpeciesIndex(tp);
		double phase_moles = 0.0;

		for (size_t kk = 0; kk < tp->nSpecies(); ++kk) {
		    phase_moles += n[k+kk];
		}
		
		  if (phase_moles < 1e-20) {
		      m_ropr[j] = m_ropf[j];
		  }
		  else if (phase_moles < 1e-7) {
		      m_ropr[j] = std::min(m_ropr[j], m_ropf[j]+tanh(1.0e9*phase_moles) * netRev);
		  }
	    }
	    else if (m_rxnPhaseIsReactant[j][p] && netFwd > 0.0) {
		ThermoPhase* tp = m_thermo[p];
		int k = m_pl->getGlobalSpeciesIndex(tp);
		double phase_moles = 0.0;

		for(size_t kk = 0; kk < tp->nSpecies(); ++kk) {
		    phase_moles += n[k+kk];
		}
		
		if (phase_moles < 1e-20) {
		    m_ropf[j] = m_ropr[j];
		}
		else if (phase_moles < 1e-7) {
		    m_ropf[j] = std::min(m_ropf[j], m_ropr[j] + tanh(1e9*phase_moles) * netFwd);
		}
	    }
	}
	
	m_ropnet[j] = m_ropf[j] - m_ropr[j];
    }
}
//==================================================================================================================================
const std::vector<double>& ReactingSurDomain::calcNetLimitedSurfaceProductionRateDensities(const double* n)
{
    updateROP();

    limitROP(n);

    getNetProductionRates(&limitedROP_[0]);

    return limitedROP_;
}
//==================================================================================================================================
const std::vector<double>& ReactingSurDomain::calcNetSurfaceROP()
{
    updateROP();
    return m_ropnet;
}
//==================================================================================================================================
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
//==================================================================================================================================
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
//==================================================================================================================================
// Note: signs have been checked to be correct in this routine.
double ReactingSurDomain::getCurrentDensityRxn(double * const currentDensityRxn) 
{
    double netCurrentDensity = 0.0;
    double ps, rs;
    if (kElectronIndex_< 0) {
	return netCurrentDensity;
    }
    size_t nr = nReactions();
    // update rates of progress -> puts this into m_ropnet[]
    updateROP();
    if (currentDensityRxn) {
	for (size_t irxn = 0; irxn < nr; irxn++) {
	    rs = m_rrxn[kElectronIndex_][irxn];
	    ps = m_prxn[kElectronIndex_][irxn];
	    double electronProd = (ps - rs) * m_ropnet[irxn];
	    currentDensityRxn[irxn] =  Faraday * electronProd;
	    netCurrentDensity += currentDensityRxn[irxn];
	}
    } else {
	for (size_t irxn = 0; irxn < nr; irxn++) {
	    rs = m_rrxn[kElectronIndex_][irxn];
	    ps = m_prxn[kElectronIndex_][irxn];
	    double electronProd = (ps - rs) * m_ropnet[irxn];
	    netCurrentDensity += Faraday * electronProd;
	}
    }
    return netCurrentDensity;
}
//==================================================================================================================================
double ReactingSurDomain::getLimitedCurrentDensityRxn(const double* n)
{
    double netCurrentDensity = 0.0;
    double ps, rs;
    if (kElectronIndex_ < 0) {
        return netCurrentDensity;
    }
    size_t nr = nReactions();
    // update rates of progress -> puts this into m_ropnet[]
    updateROP();
    limitROP(n);
    
    for (size_t irxn = 0; irxn < nr; irxn++) {
        rs = m_rrxn[kElectronIndex_][irxn];
        ps = m_prxn[kElectronIndex_][irxn];
        double electronProd = (ps - rs) * m_ropnet[irxn];
        netCurrentDensity += Faraday * electronProd;
    }   
    return netCurrentDensity;
}

//==================================================================================================================================
#ifdef DONOTREMOVE
double ReactingSurDomain::getExchangeCurrentDensityFormulation(int irxn,  doublereal* nStoich, doublereal* OCV, doublereal* io,
							       doublereal* overPotential, doublereal *beta, doublereal* resist_ptr)
{
    doublereal icurr = 0.0;
    *resist_ptr = 0.0;
  
    // This will calculate the equilibrium constant
    updateROP();
 
    double TT = m_surf->temperature();
    double rtdf = GasConstant * TT / Faraday;

    //
    //   Get the phase mole change structure, and then the number of stoichiometric electrons
    //
    AssertThrow(metalPhaseIndex_ != npos, "ReactingSurDomain::getExchangeCurrentDensityFormulation");
    RxnMolChange* rmc = rmcVector[irxn];
    double nStoichElectrons = - rmc->m_phaseChargeChange[metalPhaseIndex_];
    *nStoich = nStoichElectrons;
    //
    //  If the stoich electrons are zero, we can't proceed here (need to generalize this to general charge transfer case)
    //
    if (nStoichElectrons == 0.0) {
	*OCV = 0.0;
	*overPotential = 0.0;
        *io = 0.0;
	return 0.0;
    }
    //
    //  Calculate the Open circuit potential
    //
    getDeltaGibbs(0);
    *OCV = m_deltaG[irxn] / Faraday/ nStoichElectrons;
    //
    // Find the iBeta index and beta
    //
    size_t iBeta = npos;
    for (size_t iBetaT = 0; iBetaT < m_beta.size(); iBetaT++) {
	if (m_ctrxn[iBetaT] == (size_t) irxn) {
	    iBeta = iBetaT;
	    break;
	}
    }
    if (iBeta == npos) {
	throw CanteraError("ReactingSurDomain::getExchangeCurrentDensityFormulation()", " beta value not found");
    }
    *beta = m_beta[iBeta];
    //
    //  Find the resistance
    //
    double resist = m_ctrxn_resistivity_[iBeta];
    *resist_ptr = resist;
    //
    //   Calculate the voltage of the electrode (Note, can't really do this without the specification of the
    //   metal phase and the solution phase, as the sign of the voltage would be unspecified)
    //
    double voltage = m_phi[metalPhaseIndex_] - m_phi[solnPhaseIndex_];
    //
    //   Calculate the overpotential
    //
    double nu = voltage - *OCV;
    *overPotential = nu;
    int reactionType = m_rxntype[irxn];

    if (reactionType == BUTLERVOLMER_NOACTIVITYCOEFFS_RXN) {
	//
	// OK, the reaction rate constant contains the current density rate constant calculation
	// the rxnstoich calculation contained the dependence of the current density on the activity concentrations
	// We finish up with the ROP calculation
	//
	int iECDFormulation =  m_ctrxn_ecdf[iBeta];
	if (iECDFormulation == 0) {
	    throw CanteraError("ReactingSurDomain::getExchangeCurrentDensityFormulation()",
			       "Straight kfwrd with BUTLERVOLMER_NOACTIVITYCOEFFS_RXN not handled yet");
	}
	//
	//   Now calculate the exchange current density, io
	//
	//     Start with the exchange current reaction rate constant, which should
	//     be located in m_rfn[]. Multiply by stoichiometric electrons and perturbation factor
	//
	double ioc = m_rfn[irxn] * nStoichElectrons * m_perturb[irxn];
	//
	//   Now we need the mole fraction vector and we need the RxnOrders vector.
	//
	const RxnOrders* ro_fwd = m_ctrxn_ROPOrdersList_[iBeta];
	if (ro_fwd == 0) {
	    throw CanteraError("ReactingSurDomain::getExchangeCurrentDensityFormulation()",
			       "Forward orders pointer is zero ?!?");
	}
	double tmp = 1.0;
	double mfS = 0.0;
	const std::vector<size_t>& kinSpeciesIDs = ro_fwd->kinSpeciesIDs_;
	const std::vector<doublereal>& kinSpeciesOrders = ro_fwd->kinSpeciesOrders_;
	for (size_t j = 0; j < kinSpeciesIDs.size(); j++) {
	    size_t ks = kinSpeciesIDs[j];
	    thermo_t& th = speciesPhase(ks);
	    size_t n = speciesPhaseIndex(ks);
	    size_t klocal = ks - m_start[n];
	    mfS = th.moleFraction(klocal);
	    
	    double oo = kinSpeciesOrders[j];
	    tmp *= pow(mfS, oo);
	}
	ioc *= tmp;
	//
	//   Add in the film resistance here, later
	//
	double exp1 = nu * nStoichElectrons * (*beta) / rtdf;
	double exp2 = - nu * nStoichElectrons * (1.0 - (*beta)) / (rtdf);
        icurr = ioc * (exp(exp1) - exp(exp2));
	if (resist != 0.0) {
	    icurr = solveCurrentRes(nu, nStoichElectrons, ioc, (*beta), TT, resist, 0);
	}
	*io = ioc;

    } else if (reactionType == SURFACEAFFINITY_RXN) {
	
	size_t jjA = npos;
	for (size_t jj = 0; jj < m_numAffinityRxns; ++jj) {
	    //
	    // Get the list of data needed to calculate the affinity reaction
	    //
	    affinityRxnData& aJ = affinityRxnDataList_[jj];
	    size_t jrxn = aJ.rxn_id;
	    if (jrxn == (size_t) irxn) {
		jjA = jj;
	    }
	}
	affinityRxnData& aJ = affinityRxnDataList_[jjA];

	//PROBABLY DELETE THIS CALL SINCE IT IS CALLED BY updateROP()
	// we have a vector of standard concentrations calculated from the routine below
	//            m_StandardConc[ik]
	updateExchangeCurrentQuantities();
	//
	//   Now calculate the exchange current density, io
	//
	//   Start with the exchange current reaction rate constant, which should
	//   be located in m_rfn[]. Multiply by stoich electrons and perturbation factor
	//

	double iO = Faraday * m_rfn[irxn] * nStoichElectrons *  m_perturb[irxn];


	//
	// Undo the voltage correction term that was added in applyVoltageKfwdCorrection()
	//
	doublereal eamod = m_beta[iBeta] * deltaElectricEnergy_[irxn];
	doublereal rt = GasConstant * thermo(0).temperature();
	doublereal rrt = 1.0 / rt;
	iO /= exp(-eamod * rrt);

	double b = m_beta[iBeta];
	double omb = 1.0 - b;
        //
	// Apply the exp (beta DeltaG0) term
	//
	double mG0 =  m_deltaG0[irxn];
	if (m_beta[iBeta] > 0.0) {
	    double fac = exp(mG0 * b * rrt);
            iO *= fac;
	}


	//
	//  Eqn. (65) of writeup
	//
	for (size_t k = 0; k < m_kk; k++) {
	    doublereal reactCoeff = reactantStoichCoeff(k, irxn);
	    doublereal prodCoeff =  productStoichCoeff(k, irxn);

	    if (reactCoeff != 0.0) {
		iO *= pow(m_actConc[k],      reactCoeff*omb);
		iO /= pow(m_StandardConc[k], reactCoeff*omb);
	    }
	    if (prodCoeff != 0.0) {
		iO *= pow(m_actConc[k],      prodCoeff*b);
		iO /= pow(m_StandardConc[k], prodCoeff*b);
	    }
	}
	for (size_t k = 0; k < m_kk; k++) {
	    doublereal reactCoeff = reactantStoichCoeff(k, irxn);
	    if (reactCoeff != 0.0) {
		iO *= pow(m_actConc[k],      -reactCoeff);
		iO /= pow(m_StandardConc[k], -reactCoeff);
	    }
	}
	//const std::vector<size_t>& affinSpec_FRC = aJ.affinSpec_FRC;
	//const std::vector<doublereal>&  affinOrder_FRC = aJ.affinOrder_FRC;
	for (size_t kk = 0; kk < aJ.affinSpec_FRC.size(); ++kk) {
	    size_t kspec = aJ.affinSpec_FRC[kk];
	    double order = aJ.affinOrder_FRC[kk];
	    iO *= pow(m_actConc[kspec], order);
	}
	*io = iO;


    } else {


	//PROBABLY DELETE THIS CALL SINCE IT IS CALLED BY updateROP()
	// we have a vector of standard concentrations calculated from the routine below
	//            m_StandardConc[ik]
	updateExchangeCurrentQuantities();


	//rkc is reciprocal equilibrium constant

	const vector_fp& rf = m_rfn;
	const vector_fp& rkc= m_rkcn;

	// start with the forward reaction rate
	double iO = rf[irxn] * Faraday * nStoichElectrons;

	if (m_beta[irxn] > 0.0) {
	    iO *= pow(rkc[irxn], m_beta[irxn]);
	}
	double b = *beta;
	double omb = 1.0 - b;


	for (size_t k = 0; k < m_kk; k++) {
	    doublereal reactCoeff = reactantStoichCoeff(k, irxn);
	    doublereal prodCoeff =  productStoichCoeff(k, irxn);

	    if (reactCoeff != 0.0) {
		iO *= pow(m_actConc[k], reactCoeff*omb);
		iO *= pow(m_StandardConc[k], reactCoeff*b);
	    }
	    if (prodCoeff != 0.0) {
		iO *= pow(m_actConc[k], prodCoeff*b);
		iO /= pow(m_StandardConc[k], prodCoeff*omb);
	    }
	}
	*io = iO;

	double phiMetal = thermo(metalPhaseIndex_).electricPotential();
	double phiSoln = thermo(solnPhaseIndex_).electricPotential();
	double E = phiMetal - phiSoln;
	*overPotential = E - *OCV;

	icurr = calcCurrentDensity(*overPotential, *nStoich, *io, *beta, m_temp);
    }
    return icurr;
}
#endif
//====================================================================================================================
#ifdef DONOTREMOVE
double ReactingSurDomain::calcCurrentDensity(double nu, double nStoich, double io, double beta, double temp) const
{
     double exp1 = nu * nStoich * Faraday * beta / (GasConstant * temp);
     double exp2 = -nu * nStoich * Faraday * (1.0 - beta) / (GasConstant * temp);
     double val = io * (exp(exp1) - exp(exp2));
     return val;
}
#endif
//==================================================================================================================================
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
std::ostream& operator<<(std::ostream& s, ReactingSurDomain& rsd)
{
    ThermoPhase* th;
    ElectrodeKinetics* iK = &rsd;
    for (size_t i = 0; i < rsd.numPhases_; i++) {
        th = &(iK->thermo(i));
        std::string r = th->report(true);
        s << r;
    }
    return s;
}
//====================================================================================================================
//  Identify the metal phase and the electrons species
//  This can be taken out, because it's been moved to Cantera
//    -> Right now, change it into a double-check routine.
void ReactingSurDomain::identifyMetalPhase()
{
    metalPhaseIndex_ = -1;
    kElectronIndex_ = -1;
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
                        metalPhaseIndex_ = iph;
                        kElectronIndex_ = m_start[iph] + k;
                    }
                }
            }
        }
        if ((size_t) iph != metalPhaseIndex_) {

	    for (int i = 0; i < nr; i++) {
		RxnMolChange* rmc = rmcVector[i];
		if (rmc->m_phaseChargeChange[iph] != 0) {
		    if (rmc->m_phaseDims[iph] == 3) {
			solnPhaseIndex_ = iph;
			break;
		    }
		}
	    }
	}
    }

}
//==================================================================================================================================
void ReactingSurDomain::init()
{
    ElectrodeKinetics::init();
    m_Enthalpies_rspec.resize(m_kk, 0.0);
    m_Entropies_rspec.resize(m_kk, 0.0);
    m_GibbsOCV_rspec.resize(m_kk, 0.0);
    m_Enthalpies_Before_rspec.resize(m_kk, 0.0);
    m_Entropies_Before_rspec.resize(m_kk, 0.0);
    m_Gibbs_Before_rspec.resize(m_kk, 0.0);
    speciesProductionRates_.resize(m_kk, 0.0);
    speciesCreationRates_.resize(m_kk, 0.0);
    speciesDestructionRates_.resize(m_kk, 0.0);
    KintoPLSpeciesIndex_.resize(m_kk, npos);
}
//==================================================================================================================================
void ReactingSurDomain::finalize()
{
    ElectrodeKinetics::finalize();
    deltaGRxn_Before_.resize(m_ii, 0.0);
    deltaHRxn_Before_.resize(m_ii, 0.0);
    deltaSRxn_Before_.resize(m_ii, 0.0);
    deltaGRxnOCV_Before_.resize(m_ii, 0.0);
    limitedROP_.resize(m_ii, 0.0);
}
//==================================================================================================================================
/*
 *      @param iskin   If this is zero or positive, we will formulate a kinetics mechanism using an 
 *                     interfacial kinetics object.
 */
bool ReactingSurDomain::
importFromPL(Cantera::PhaseList* const pl, int iskin)
{
    try {
        int iph;
        //
        //  Store the PhaseList as a shallow pointer within the object
        //
        m_pl = pl;

        XML_Node* kinXMLPhase = 0;
        ThermoPhase* kinPhase = 0;

        if (iskin < 0 || iskin >= pl->nSurPhases()) {
           throw Electrode_Error("ReactingSurDomain::importFromPL()",
                                 "index of surface reaction not within bounds");
        }
        if (iskin >= 0) {
            kinXMLPhase = pl->surPhaseXMLNode(iskin);
            kinPhase = &(pl->surPhase(iskin));
        } 

        int nPhasesFound = pl->nSurPhases() + pl->nVolPhases();
        /*
         * Resize the internal list of pointers and get a pointer to the vacant ThermoPhase pointer
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
	KintoPLSpeciesIndex_.resize(m_kk, npos);
	//size_t kKinIndex = 0;
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
		KintoPLSpeciesIndex_[m_start[kph] + k] = k + PLkstart;
            }
        }

        /*
         * Resize the arrays based on kinetic species number
         */
        speciesProductionRates_.resize(m_kk, 0.0);
        speciesCreationRates_.resize(m_kk, 0.0);
        speciesDestructionRates_.resize(m_kk, 0.0);


	m_Enthalpies_rspec.resize(m_kk, 0.0);
	m_Entropies_rspec.resize(m_kk, 0.0);
	m_GibbsOCV_rspec.resize(m_kk, 0.0);

        /*
         * Resize the arrays based on the number of reactions
         */ 
        deltaGRxn_Before_.resize(m_ii, 0.0);
	deltaHRxn_Before_.resize(m_ii, 0.0);
	deltaSRxn_Before_.resize(m_ii, 0.0);
	limitedROP_.resize(m_ii, 0.0);

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
    //
    //  Go get the model from the factory routine
    //
    OCVmodel_ = newRSD_OCVmodel(ocv_ptr_->OCVModel);

    OCVmodel_->initialize(this, *ocv_ptr_);
   
    //
    // Now setup the internal structures within the model to calculate the relative extent
    // Find the phase id and phase name of the replaced global species. We will assume that it is also
    // the solid phase where we will get the relative extent.
    
    int phase_id = m_pl->getPhaseIndexFromGlobalSpeciesIndex(ocv_ptr_->replacedGlobalSpeciesID);
    string phaseName = m_pl->phaseName(phase_id);
    //
    //  Since the pointers must all be the same, we look up the ThermoPhase pointer in the phase list
    //  and send that to the RSD_OCVmodel object
    //
    ThermoPhase* tp = &(m_pl->thermo(phase_id));
    //
    //  We also find the species which we will use as the relative extent variable. Typically it is 
    //  not the same as the species we used as the replaced species, at least for anodes
    //
    size_t kspec = ocv_ptr_->MF_DoD_LocalSpeciesID;
    AssertThrow(kspec != npos, "ReactingSurDomain::addOCVoverride(): MF_DoD_LocalSpeciesID is not set");
    //
    //  Set up the relative extent capability within the OCVmodel. Note we may have to expand this
    //  model more in the future 
    //
    OCVmodel_->setup_RelExtent(tp, kspec);

    kReplacedSpeciesRS_  = PLtoKinSpeciesIndex_[ocv_ptr_->replacedGlobalSpeciesID];
 
}
//====================================================================================================================
void ReactingSurDomain::deriveEffectiveChemPot()
{
    /*
     *   Find and store the temperature
     */
    double TT = m_surf->temperature();
    //  Wastefull, but for now get a complete SSG and G vector.
    for (size_t n = 0; n < nPhases(); n++) {
	size_t nsp = thermo(n).nSpecies();
	size_t kinSpecOff = m_start[n];
	thermo(n).getChemPotentials(DATA_PTR(m_mu) + m_start[n]);
	thermo(n).getStandardChemPotentials(DATA_PTR(m_mu0) + m_start[n]);

        if (OCVmodel_->OCV_Format_ == 0) {
	    if (n == (size_t) solnPhaseIndex_) {
		mdpUtil::mdp_copy_dbl_1(DATA_PTR(m_GibbsOCV_rspec) + kinSpecOff, DATA_PTR(m_mu0) + kinSpecOff, nsp);
            }  else {
		mdpUtil::mdp_copy_dbl_1(DATA_PTR(m_GibbsOCV_rspec) + kinSpecOff, DATA_PTR(m_mu) + kinSpecOff, nsp);
	    }
        } else if (OCVmodel_->OCV_Format_ == 1) {
	    if (n == (size_t) solnPhaseIndex_) {
		mdpUtil::mdp_copy_dbl_1(DATA_PTR(m_GibbsOCV_rspec) + kinSpecOff, DATA_PTR(m_mu0) + kinSpecOff, nsp);
            }  else {
		mdpUtil::mdp_copy_dbl_1(DATA_PTR(m_GibbsOCV_rspec) + kinSpecOff, DATA_PTR(m_mu) + kinSpecOff, nsp);
	    }
        } else if (OCVmodel_->OCV_Format_ == 2) {
	    mdpUtil::mdp_copy_dbl_1(DATA_PTR(m_GibbsOCV_rspec) + kinSpecOff, DATA_PTR(m_mu) + kinSpecOff, nsp);
        } else if (OCVmodel_->OCV_Format_ == 3) {
	    mdpUtil::mdp_copy_dbl_1(DATA_PTR(m_GibbsOCV_rspec) + kinSpecOff, DATA_PTR(m_mu0) + kinSpecOff, nsp);
        } else {
	    throw CanteraError("", "not implemnented"); 
        }
    }
    //
    //  Calculate the delta G and the open circuit voltage normally
    //
    //
    //
    //
    //  Get the reaction delta based on the mixed G and G_SS values accummulated above
    //
    getReactionDelta(DATA_PTR(m_GibbsOCV_rspec), DATA_PTR(deltaGRxnOCV_Before_));

    double phiRxnOrig = 0.0;

    //
    //  If we don't have an electrode reaction bail as being confused
    //
    if (metalPhaseIndex_ < 0) {
	throw CanteraError("", "shouldn't be here");
    }
    //
    //  If we don't have an open circuit potential override situation, bail out of the routine
    if (!ocv_ptr_) {
	return;
    }

    //
    //  Figure out the reaction id for the override
    //
    int rxnID = ocv_ptr_->rxnID;
    int rxnID_deltaS = ocv_ptr_->rxnID_deltaS;
    //
    //  Get a pointer to the RxnMolChange struct, which contains more info about the reaction
    //
    RxnMolChange* rmc = rmcVector[rxnID];
    //
    //   Find the number of stoichiometric electrons in the reaction
    //
    double nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseIndex_];
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
    phiRxnOrig = deltaGRxnOCV_Before_[rxnID_deltaS] / Faraday / nStoichElectrons;
    //
    //    In order to calculate the OCV, we need the relative extent of reaction value. This is determined
    //    automatically within the OCVmodel object given that the ThermoPhase is current.
    //    We report the value here.
    //
    OCVmodel_->RelExtent();
    //
    //    Now calculate the OCV value to be used from the fit presumably from a fit to experiment.
    //
    double phiRxnExp = OCVmodel_->OCV_value(TT);
    //
    //    Calculate the deltaGibbs needed to turn phiRxnOrig into phiRxnExp
    //
    double deltaG_Exp_delta = (phiRxnExp - phiRxnOrig) * Faraday * nStoichElectrons;
    //
    //
    //   
    double fstoich =  reactantStoichCoeff( kReplacedSpeciesRS_,  ocv_ptr_->rxnID_deltaS);
    double rstoich =  productStoichCoeff( kReplacedSpeciesRS_, ocv_ptr_->rxnID_deltaS);
    double nstoich = rstoich - fstoich;
    deltaG_species_ = deltaG_Exp_delta / nstoich;
    //
    // Calculate new G value for replaced species that will yield the experimental open circuit voltage
    //
    m_GibbsOCV_rspec[kReplacedSpeciesRS_] += deltaG_species_;
    m_mu[kReplacedSpeciesRS_]  += deltaG_species_;
    m_mu0[kReplacedSpeciesRS_] += deltaG_species_;

    //
    // Now recalc deltaG and calc OCV   HKM This calculation has checked out for the MCMB model
    //   We should now be at the exp OCV
#ifdef DEBUG_NEW
    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_mu), DATA_PTR(deltaGRxnOCV_Before_));
    phiRxnOrig = deltaGRxnOCV_Before_[rxnID_deltaS] / Faraday / nStoichElectrons;
    //printf(" phiRxnOrig_halfcell  =  %g\n",  phiRxnOrig );
    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_grt), DATA_PTR(deltaGRxnOCV_Before_));
    phiRxnOrig = deltaGRxnOCV_Before_[rxnID_deltaS] / Faraday / nStoichElectrons;
    //printf(" phiRxnOrig_new  =  %g\n",  phiRxnOrig );
#endif
}
//==================================================================================================================================
void ReactingSurDomain::deriveEffectiveThermo()
{
    /*
     *   Find and store the temperature
     */
    double TT = m_surf->temperature();

    //  Wastefull, but for now get a complete SSG and G vector.
    // Also
    for (size_t n = 0; n < numPhases_; n++) {
	size_t nsp = thermo(n).nSpecies();
	size_t kinSpecOff = m_start[n];
	thermo(n).getChemPotentials(DATA_PTR(m_mu) + kinSpecOff);
	thermo(n).getStandardChemPotentials(DATA_PTR(m_mu0) + kinSpecOff);

	if (OCVmodel_->OCV_Format_ == 0) {
	    if (n == (size_t) solnPhaseIndex_) {
		mdpUtil::mdp_copy_dbl_1(DATA_PTR(m_GibbsOCV_rspec) + kinSpecOff, DATA_PTR(m_mu0) + kinSpecOff, nsp);
            }  else {
		mdpUtil::mdp_copy_dbl_1(DATA_PTR(m_GibbsOCV_rspec) + kinSpecOff, DATA_PTR(m_mu) + kinSpecOff, nsp);
	    }
        } else if (OCVmodel_->OCV_Format_ == 1) {
	    if (n == (size_t) solnPhaseIndex_) {
		mdpUtil::mdp_copy_dbl_1(DATA_PTR(m_GibbsOCV_rspec) + kinSpecOff, DATA_PTR(m_mu0) + kinSpecOff, nsp);
            }  else {
		mdpUtil::mdp_copy_dbl_1(DATA_PTR(m_GibbsOCV_rspec) + kinSpecOff, DATA_PTR(m_mu) + kinSpecOff, nsp);
	    }
        } else if (OCVmodel_->OCV_Format_ == 2) {
	    mdpUtil::mdp_copy_dbl_1(DATA_PTR(m_GibbsOCV_rspec) + kinSpecOff, DATA_PTR(m_mu) + kinSpecOff, nsp);
        } else if (OCVmodel_->OCV_Format_ == 3) {
	    mdpUtil::mdp_copy_dbl_1(DATA_PTR(m_GibbsOCV_rspec) + kinSpecOff, DATA_PTR(m_mu0) + kinSpecOff, nsp);
        } else {
	    throw CanteraError("deriveEffectiveThermo()", 
			       "OCV_Format model " + int2str(OCVmodel_->OCV_Format_) + " not implemented"); 
        }

    }
    // Calc Delta G's as before
    getDeltaGibbs_Before(DATA_PTR(deltaGRxn_Before_));
    getDeltaEnthalpy_Before(DATA_PTR(deltaHRxn_Before_));
    getDeltaEntropy_Before(DATA_PTR(deltaSRxn_Before_));

    
#ifdef DEBUG_MODE
    double dgc = deltaHRxn_Before_[1] - TT * deltaSRxn_Before_[1];
    bool de = esmodel::doubleEqual(dgc, deltaGRxn_Before_[1]);
    if (!de) {
	printf("error!\n");
	exit(-1);
    }
#endif
    //
    //
    //
    //  Get the reaction delta based on the mixed G and G_SS values accummulated above
    //
    getReactionDelta(DATA_PTR(m_GibbsOCV_rspec), DATA_PTR(deltaGRxnOCV_Before_));

    double phiRxnOrig = 0.0;

    //
    //  If we don't have an electrode reaction bail as being confused
    //
    if (metalPhaseIndex_ < 0) {
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
    int rxnID_deltaS =  ocv_ptr_->rxnID_deltaS;
    //
    //  Get a pointer to the RxnMolChange struct, which contains more info about the reaction
    //
    RxnMolChange* rmc = rmcVector[rxnID];
    //
    //   Find the number of stoichiometric electrons in the reaction
    //
    double nStoichElectrons = -rmc->m_phaseChargeChange[metalPhaseIndex_];
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
    phiRxnOrig = deltaGRxnOCV_Before_[rxnID_deltaS] / Faraday / nStoichElectrons;
    //
    //  Calculate the dOCV/dT from that relation that would occur if it was not being overwritten
    //  -> where is deltaSRxn_Before_[rxnID] calculated
    doublereal d_phiRxnOrig_dT = -deltaSRxn_Before_[rxnID_deltaS] / Faraday / nStoichElectrons;
    //
    //    In order to calculate the OCV, we need the relative extent of reaction value. This is determined
    //    automatically within the OCVmodel object given that the ThermoPhase is current.
    //    We report the value here.
    //
    (void) OCVmodel_->RelExtent();
    //
    //    Now calculate the OCV value to be used from the fit presumably from a fit to experiment.
    //
    double phiRxnExp = OCVmodel_->OCV_value(TT);
    //
    //   Now calculate the dOCV/dT value to be used presumably from a fit
    //
    double d_phiRxnExp_dT = OCVmodel_->OCV_dvaldT(TT);
    //
    //    Calculate the deltaGibbs needed to turn phiRxnOrig into phiRxnExp
    //
    //double deltaG_Exp = (phiRxnExp) * Faraday * nStoichElectrons;
    double deltaG_Exp_delta = (phiRxnExp - phiRxnOrig) * Faraday * nStoichElectrons;
    //
    //    Calculate the deltaEntropy needed to turn dphiRxnOrigdT into dphiRxnExpdT
    //
    //double deltaS_Exp = (- d_phiRxnExp_dT) * Faraday * nStoichElectrons;
    double deltaS_Exp_delta = (- d_phiRxnExp_dT + d_phiRxnOrig_dT) * Faraday * nStoichElectrons;
    //
    //
    //   
    double fstoich =  reactantStoichCoeff( kReplacedSpeciesRS_,  ocv_ptr_->rxnID_deltaS);
    double rstoich =  productStoichCoeff(kReplacedSpeciesRS_, ocv_ptr_->rxnID_deltaS);
    double nstoich = rstoich - fstoich;
    deltaG_species_ = deltaG_Exp_delta / nstoich;
    deltaS_species_ = deltaS_Exp_delta / nstoich;
    deltaH_species_ = deltaG_species_ + TT * deltaS_species_;
    
    //
    // calculate new G value for replaced species that will yield the experimental open circuit voltage
    //
    m_GibbsOCV_rspec[kReplacedSpeciesRS_] += deltaG_species_;
    m_mu[kReplacedSpeciesRS_] += deltaG_species_;
    m_mu0[kReplacedSpeciesRS_] += deltaG_species_;

    //
    // Calculate the replaced entropies and enthalpies
    //
    mdpUtil::mdp_copy_dbl_1(DATA_PTR(m_Entropies_rspec), DATA_PTR(m_Entropies_Before_rspec), m_kk);
    mdpUtil::mdp_copy_dbl_1(DATA_PTR(m_Enthalpies_rspec), DATA_PTR(m_Enthalpies_Before_rspec), m_kk);
    m_Entropies_rspec[kReplacedSpeciesRS_] += deltaS_species_;
    m_Enthalpies_rspec[kReplacedSpeciesRS_] += deltaH_species_;
    //
    //   Now recalc deltaG and calc OCV   HKM This calculation has checked out for the MCMB model
    //   We should now be at the exp OCV
#ifdef DEBUG_NEW
    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_mu), DATA_PTR(deltaGRxnOCV_Before_));
    phiRxnOrig = deltaGRxnOCV_Before_[rxnID_deltaS] / Faraday / nStoichElectrons;
    //printf(" phiRxnOrig_halfcell  =  %g\n",  phiRxnOrig );
    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(mGibbsOCV_rspec), DATA_PTR(deltaGRxnOCV_Before_));
    phiRxnOrig = deltaGRxnOCV_Before_[rxnID_deltaS] / Faraday / nStoichElectrons;
    //printf(" phiRxnOrig_new  =  %g\n",  phiRxnOrig );
#endif
}
//==================================================================================================================================
void ReactingSurDomain::updateMu0()
{
    /*
     * Get the vector of standard state electrochemical potentials for species in the Interfacial
     * kinetics object and store it in m_mu0[] and in m_mu0_Kc[]
     */
    size_t nsp, ik = 0;
    doublereal rt = GasConstant * thermo(0).temperature();
    size_t np = nPhases();
    for (size_t n = 0; n < np; n++) {
        thermo(n).getStandardChemPotentials(DATA_PTR(m_mu0) + m_start[n]);
    }
    //
    // This overwrite m_mu0 for the selected species so that we get the experimental OCV
    // HKM -> This is insufficient!! Need to overwrite the complete thermo functions for that
    //        species, so that deltaS deltaH and deltaCp are modified as well! 5/1/2015
    //
    if (ocv_ptr_) {
       deriveEffectiveChemPot();
    }

    for (size_t n = 0; n < np; n++) {
        nsp = thermo(n).nSpecies();
        for (size_t k = 0; k < nsp; k++) {
            m_mu0_Kc[ik] = m_mu0[ik] + Faraday * m_phi[n] * thermo(n).charge(k);
            m_mu0_Kc[ik] -= rt * thermo(n).logStandardConc(k);
            ik++;
        }
    }
}
//==================================================================================================================================
// Modification for OCV override when called for.
//         ( works for MCMB )
void ReactingSurDomain::getDeltaGibbs(doublereal* deltaG)
{
     /*
      * Get the chemical potentials of the species in the all of the phases used in the 
      * kinetics mechanism -> We store them in the vector m_mu[]. which later can get modified for OCV override
      * issues.
      */
     for (size_t n = 0; n < numPhases_; n++) {
         m_thermo[n]->getChemPotentials(DATA_PTR(m_mu) + m_start[n]);
     }

     //  If have an open circuit potential override situation, do extra work
     if (ocv_ptr_) {
	 //deriveEffectiveChemPot();
	 deriveEffectiveThermo();
     }

    //
    // Use the stoichiometric manager to find deltaG for each
    // reaction.
    //
    getReactionDelta(DATA_PTR(m_mu), DATA_PTR(m_deltaG));
    if (deltaG != 0 && (DATA_PTR(m_deltaG) != deltaG)) {
	for (size_t j = 0; j < m_ii; ++j) {
	    deltaG[j] = m_deltaG[j];
	}
    }
}
//==================================================================================================================================
void ReactingSurDomain::getDeltaGibbs_Before(doublereal* const deltaG)
{
    for (size_t n = 0; n < numPhases_; n++) {
         m_thermo[n]->getChemPotentials(DATA_PTR(m_Gibbs_Before_rspec) + m_start[n]);
    }
    if (deltaG) {
	getReactionDelta(DATA_PTR(m_Gibbs_Before_rspec), DATA_PTR(deltaG));
    }
}
//==================================================================================================================================
// Modification for OCV override when called for.
//         ( works for MCMB )
void ReactingSurDomain::getDeltaElectrochemPotentials(doublereal* deltaG)
{
     /*
      * Get the chemical potentials of the species in the all of the phases used in the 
      * kinetics mechanism -> We store them in the vector m_mu[]. which later can get modified for OCV override
      * issues.
      */
     for (size_t n = 0; n < numPhases_; n++) {
         m_thermo[n]->getChemPotentials(DATA_PTR(m_mu) + m_start[n]);
     }

     //  If have an open circuit potential override situation, do extra work
     if (ocv_ptr_) {
	 deriveEffectiveChemPot();
     }

     for (size_t n = 0; n < numPhases_; n++) {
	 ThermoPhase& tp = thermo(n);
	 double ve = Faraday * tp.electricPotential();
	 size_t nsp =  tp.nSpecies();
	 for (size_t k = 0; k < nsp; k++) {
	     m_grt[m_start[n] + k] = m_mu[m_start[n] + k] + ve * tp.charge(k);
        }
     }
    //
    // Use the stoichiometric manager to find deltaG for each
    // reaction.
    //
    getReactionDelta(DATA_PTR(m_grt), DATA_PTR(deltaG));
}
//==================================================================================================================================
// Modification for OCV override when called for.
void ReactingSurDomain::getDeltaEnthalpy(doublereal* deltaH)
{
    /*
     *  Get the partial molar enthalpy of all species
     */
    for (size_t n = 0; n < numPhases_; n++) {
        thermo(n).getPartialMolarEnthalpies(DATA_PTR(m_Enthalpies_rspec) + m_start[n]);
    }
    /*
     *  If have an open circuit potential override situation, do extra work involving fixing thermo
     *    That work gets posted to m_Entropies_rspec
     */
    if (ocv_ptr_) {
         deriveEffectiveThermo();
    }
    /*
     *  Use the stoichiometric manager to find deltaH for each reaction.
     *    (we do not sture the deltaH vector, as I can't think of a reason to do so).
     */
    if (deltaH) {
	getReactionDelta(DATA_PTR(m_Enthalpies_rspec), DATA_PTR(deltaH));
    }
}
//=======================================================================================================================
// Modification for OCV override when called for.
void ReactingSurDomain::getDeltaEnthalpy_Before(doublereal* const deltaH)
{
    for (size_t n = 0; n < numPhases_; n++) {
        thermo(n).getPartialMolarEnthalpies(DATA_PTR(m_Enthalpies_Before_rspec) + m_start[n]);
    }
    if (deltaH) {
	getReactionDelta(DATA_PTR(m_Enthalpies_Before_rspec), DATA_PTR(deltaH));
    }
}
//=======================================================================================================================
// Modification for OCV override when called for
void ReactingSurDomain::getDeltaEntropy(doublereal* deltaS)
{
    /*
     *   Get the partial molar entropy of all species in all phases of the kinetics object
     */
    for (size_t n = 0; n < numPhases_; n++) {
        thermo(n).getPartialMolarEntropies(DATA_PTR(m_Entropies_rspec) + m_start[n]);
    }
    /*
     *  If have an open circuit potential override situation, do extra work involving fixing thermo
     *    That work gets posted to m_Entropies_rspec
     */
    if (ocv_ptr_) {
         deriveEffectiveThermo();
    }
    /*
     *  Use the stoichiometric manager to find deltaS for each reaction.
     *    (we do not store the deltaS vector, as I can't think of a reason to do so).
     */
    getReactionDelta(DATA_PTR(m_Entropies_rspec), DATA_PTR(deltaS));
}
//=======================================================================================================================
// Modification for OCV override when called for
void ReactingSurDomain::getDeltaEntropy_Before(doublereal* const deltaS)
{
    for (size_t n = 0; n < numPhases_; n++) {
        thermo(n).getPartialMolarEntropies(DATA_PTR(m_Entropies_Before_rspec) + m_start[n]);
    }
    if (deltaS) {
	getReactionDelta(DATA_PTR(m_Entropies_Before_rspec), DATA_PTR(deltaS));
    }
}
//=======================================================================================================================
// This gets the deltaG for each reaction in the mechanism, but using the standard state
// chemical potential for the electrolyte.
/*
 *          @param  deltaG_special  DeltaG for each reaction using standard state chemical
 *                                  potentials for the electrolyte. 
 *                     length = nReactions(), J/kmol
 */
void ReactingSurDomain::getDeltaGibbs_electrolyteSS(doublereal* deltaG_special)
{
     /*
      * Use m_grt, because this vector we are creating is a mixed beast
      */
     for (size_t n = 0; n < numPhases_; n++) {
	 ThermoPhase* tp = m_thermo[n];
	 if (n == solnPhaseIndex_) {
	     tp->getStandardChemPotentials(DATA_PTR(m_grt) + m_start[n]);
	 } else {
	     tp->getChemPotentials(DATA_PTR(m_grt) + m_start[n]);
	 }
     }

     if (ocv_ptr_) {
	 deriveEffectiveChemPot();
     }
     /*
      * Use the stoichiometric manager to find deltaG for each
      * reaction.
      */
     getReactionDelta(DATA_PTR(m_grt), deltaG_special);
}
//==================================================================================================================================
// Modification for OCV override when called for
void ReactingSurDomain::getDeltaSSEnthalpy(doublereal* deltaH)
{
   /*
     *  Get the standard state entropy of the species.
     *  We define these here as the entropies of the pure
     *  species at the temperature and pressure of the solution.
     */
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getEnthalpy_RT(DATA_PTR(m_grt) + m_start[n]);
    }
    doublereal RT = GasConstant * m_surf->temperature();
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= RT;
    }
    /*
     *  If have an open circuit potential override situation, do extra work involving fixing thermo
     *  
     */
    if (ocv_ptr_) {
         deriveEffectiveThermo();
	 m_grt[kReplacedSpeciesRS_] += deltaH_species_;
    }

    /*
     * Use the stoichiometric manager to find deltaH for each
     * reaction.
     */
    getReactionDelta(DATA_PTR(m_grt), deltaH);  
}
//=======================================================================================================================
// Modification for OCV override when called for
void ReactingSurDomain::getDeltaSSGibbs(doublereal* deltaG)
{
    /*
     *  Get the standard state gibbs free energy of the species.
     */
    for (size_t n = 0; n < nPhases(); n++) {
	thermo(n).getStandardChemPotentials(DATA_PTR(m_mu0) + m_start[n]);
    }   
    /*
     *  If have an open circuit potential override situation, do extra work involving fixing thermo
     */
    if (ocv_ptr_) {
         deriveEffectiveChemPot();
    }
    /*
     * Use the stoichiometric manager to find deltaG for each reaction.
     */
    getReactionDelta(DATA_PTR(m_mu0), deltaG);  
}
//==================================================================================================================================
void ReactingSurDomain::getDeltaSSEntropy(doublereal* deltaS)
{
    /*
     *  Get the standard state entropy of the species.
     *  We define these here as the entropies of the pure species at the temperature and pressure of the solution.
     */
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getEntropy_R(DATA_PTR(m_grt) + m_start[n]);
    }
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= GasConstant;
    }
    /*
     *  If have an open circuit potential override situation, do extra work involving fixing the thermo to fit experiment
     */
    if (ocv_ptr_) {
         deriveEffectiveThermo();
	 m_grt[kReplacedSpeciesRS_] += deltaS_species_;
    }
    /*
     * Use the stoichiometric manager to find deltaS for each reaction.
     */
    getReactionDelta(DATA_PTR(m_grt), deltaS); 
}
//==================================================================================================================================
void ReactingSurDomain::getOCVThermoOffsets_ReplacedSpecies(double& deltaG_species, double& deltaH_species, double& deltaS_species)
{
    deltaG_species = 0.0;
    deltaH_species = 0.0;
    deltaS_species = 0.0;
    if (ocv_ptr_) {
	deriveEffectiveThermo();
	deltaG_species = deltaG_species_;
	deltaH_species = deltaH_species_;
	deltaS_species = deltaS_species_;
    }
}
//==================================================================================================================================
void ReactingSurDomain::setState_TP(double temp, double pres) 
{
    for (size_t n = 0; n < nPhases(); n++) {
	m_thermo[n]->setState_TP(temp, pres);
    }
}
//==================================================================================================================================
} // End of namespace cantera
//----------------------------------------------------------------------------------------------------------------------------------

