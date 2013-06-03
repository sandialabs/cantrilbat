

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
#include "ReactionDataElectrode.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/base/stringUtils.h"

using namespace std;
using namespace ctml;

namespace Cantera {

  typedef XML_Node                 node_t;
  typedef const vector<XML_Node*>  nodeset_t;

  ReactionDataElectrode::ReactionDataElectrode() :
     ReactionData()
  {

  }


  ReactionDataElectrode::~ReactionDataElectrode() {

  }

  /**
   * getArrhenius() parses the xml element called Arrhenius. 
   * The Arrhenius expression is
   * \f[        k =  A T^(b) exp (-E_a / RT). \f]
   */
  static void getArrhenius(const XML_Node& node, int& highlow, 
			   doublereal& A, doublereal& b, doublereal& E) {
        
    if (node["name"] == "k0") 
      highlow = 0;
    else highlow = 1;
    /*
     * We parse the children for the A, b, and E conponents.
     */
    A = getFloat(node, "A", "-");
    b = getFloat(node, "b");
    E = getFloat(node, "E", "actEnergy");
    E /= GasConstant;
  }           

  void ReactionDataElectrode::getCoverageDependence(const XML_Node& node,
						    ThermoPhase& surfphase) {
    vector<XML_Node*> n_cov;
    node.getChildren("coverage", n_cov);
    int k, nc = static_cast<int>(n_cov.size());
    doublereal e;
    string spname;
    if (nc > 0) {
      for (int n = 0; n < nc; n++) {
	const XML_Node& cnode = *n_cov[n];
	spname = cnode["species"];
	k = surfphase.speciesIndex(spname);
	cov.push_back(doublereal(k));
	cov.push_back(getFloat(cnode, "a"));
	cov.push_back(getFloat(cnode, "m"));
	e = getFloat(cnode, "e", "actEnergy");
	cov.push_back(e/GasConstant);
      }
    }
  }     

   /**
     * getStick() processes the XML element called Stick that specifies
     * the sticking coefficient reaction. This routine will 
     * translate the sticking coefficient value into a "normal"
     * rate constant for the surface reaction.
     *
     *  Output
     * -----------
     * Output is the normal Arrhenius expressions for a surface
     * reaction rate constant.
     * 
     *   A - units such that rate of rxn has kmol/m^2/s when
     *       A is multiplied by activity concentrations of 
     *       reactants in the normal manner.
     *   n - unitless
     *   E - Units 1/Kelvin
     */
  void ReactionDataElectrode::getStick(const XML_Node& node, Kinetics& kin,
				       doublereal& A, doublereal& b, doublereal& E) {
    int nr = reactants.size();
    int k, klocal, not_surf = 0;
    int np = 0;
    doublereal f = 1.0;
    doublereal n_order;
    /*
     * species is the name of the special reactant whose surface
     * flux rate will be calculated.
     *      isp = species # in the local phase
     *      ispKinetics = species # in the kinetics object
     *      ispPhaseIndex = phase # of the special species
     */
    string spname = node["species"];
    ThermoPhase& th = kin.speciesPhase(spname);
    int isp = th.speciesIndex(spname);
    int ispKinetics = kin.kineticsSpeciesIndex(spname);
    int ispPhaseIndex = kin.speciesPhaseIndex(ispKinetics);
  
    double ispMW = th.molecularWeights()[isp];
    double sc;

    // loop over the reactants
    for (int n = 0; n < nr; n++) {
      k = reactants[n];
      n_order = rorder[n];    // stoich coeff

      // get the phase species k belongs to
      np = kin.speciesPhaseIndex(k);
      const ThermoPhase& p = kin.thermo(np);

      // get the local index of species k in this phase
      klocal = p.speciesIndex(kin.kineticsSpeciesName(k));

      // if it is a surface species, divide f by the standard
      // concentration for this species, in order to convert
      // from concentration units used in the law of mass action
      // to coverages used in the sticking probability
      // expression
      if (p.eosType() == cSurf || p.eosType() == cEdge) {
	sc = p.standardConcentration(klocal);
	f /= pow(sc, n_order);
      }   
      // Otherwise:
      else {
	// We only allow one species to be in the phase
	// containing the special sticking coefficient
	// species.
	if (ispPhaseIndex == np) {
	  not_surf++;
	} 
	// Other bulk phase species on the other side
	// of ther interface are treated like surface
	// species.
	else {
	  sc = p.standardConcentration(klocal);
	  f /= pow(sc, n_order);
	}
      }
    }
    if (not_surf != 1) {
      throw CanteraError("getStick",
			 "reaction probabilities can only be used in "
			 "reactions with exactly 1 gas/liquid species.");
    }

    doublereal cbar = sqrt(8.0*GasConstant/(Pi*ispMW));
    A = 0.25 * getFloat(node, "A", "-") * cbar * f;
    b = getFloat(node, "b") + 0.5;
    E = getFloat(node, "E", "actEnergy");
    E /= GasConstant;
  }
     
  /*
   * Extract the rate coefficient for a reaction from the xml node, kf.
   * kf should point to a XML element named "rateCoeff".
   * rdata is the partially filled ReactionData object for the reaction.
   * This function will fill in more fields in the ReactionData object.
   */
  void ReactionDataElectrode::getRateCoefficient(const node_t& kf, Kinetics& kin, 
						 int negA) {

    int nc = kf.nChildren();
    nodeset_t& kf_children = kf.children();
    vector_fp clow(3,0.0), chigh(3,0.0);
    //        int nr = nReacMolecules(rdata);
    for (int m = 0; m < nc; m++) {
      const node_t& c = *kf_children[m];
      string nm = c.name();
      int highlow=0;

      if (nm == "Arrhenius") {
	vector_fp coeff(3);
	if (c["type"] == "stick") {
	  getStick(c, kin, coeff[0], coeff[1], coeff[2]);
	  chigh = coeff;
	}
	else {
	  getArrhenius(c, highlow, coeff[0], coeff[1], coeff[2]);
	  if (highlow == 1 || reactionType == THREE_BODY_RXN 
	      || reactionType == ELEMENTARY_RXN) 
	    chigh = coeff;
	  else clow = coeff;
	}
	if (reactionType == SURFACE_RXN) {
	  int spi = kin.surfacePhaseIndex();
	  AssertTrace(spi >= 0);
	  getCoverageDependence(c, kin.thermo(spi));
	}

	if (coeff[0] <= 0.0 && negA == 0) {
	  throw CanteraError("getRateCoefficient", 
			     "negative or zero A coefficient for reaction "+int2str(number));
	}
      }
      else if (nm == "electrochem") {
	beta = fpValue(c["beta"]);
      }
    }
    /*
     * Store the coefficients in the ReactionData object for return
     * from this function.
     */
    if (reactionType == CHEMACT_RXN) 
      rateCoeffParameters = clow;        
    else
      rateCoeffParameters = chigh;

    if (reactionType == FALLOFF_RXN)
      auxRateCoeffParameters = clow;
    else if (reactionType == CHEMACT_RXN) 
      auxRateCoeffParameters = chigh; 
  }


}
