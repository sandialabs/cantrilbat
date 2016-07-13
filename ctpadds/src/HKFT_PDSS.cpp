/*
 * $Id: HKFT_PDSS.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */
#include "ct_defs.h"
#include "xml.h"
#include "ctml.h"
#include "HKFT_PDSS.h"

#include "ThermoPhase.h"

using namespace std;

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif 
{
  /**
   * Basic list of constructors and duplicators
   */



  HKFT_PDSS::HKFT_PDSS(ThermoPhase *tp, int spindex) :
    PDSS(tp, spindex)
  {
  }


  HKFT_PDSS::HKFT_PDSS(ThermoPhase *tp, int spindex, string inputFile, string id) :
    PDSS(tp, spindex, inputFile, id)
  {
  }

  HKFT_PDSS::HKFT_PDSS(ThermoPhase *tp, int spindex, XML_Node& phaseRoot, string id) :
    PDSS(tp, spindex, phaseRoot, id)
  {
  }


  HKFT_PDSS::HKFT_PDSS(const HKFT_PDSS &b) :
    PDSS(b)
  {
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy construtor.
     */
    *this = b;
  }

  /**
   * Assignment operator
   */
  HKFT_PDSS& HKFT_PDSS::operator=(const HKFT_PDSS&b) {
    if (&b == this) return *this;
    m_tp = b.m_tp;
    m_spindex = b.m_spindex;
    m_temp = b.m_temp;
    m_dens = b.m_dens;
    m_mw = b.m_mw;
    return *this;
  }
  
  /**
   * Destructor for the HKFT_PDSS class
   */
  HKFT_PDSS::~HKFT_PDSS() { 
  }
  
  void HKFT_PDSS::constructHKFT_PDSS(ThermoPhase *tp, int spindex) {
    initThermo();
  }
  
   
  /**
   * constructHKFT_PDSSXML:
   *
   * Initialization of a Debye-Huckel phase using an
   * xml file.
   *
   * This routine is a precursor to initThermo(XML_Node*)
   * routine, which does most of the work.
   *
   * @param infile XML file containing the description of the
   *        phase
   *
   * @param id  Optional parameter identifying the name of the
   *            phase. If none is given, the first XML
   *            phase element will be used.
   */
  void HKFT_PDSS::constructHKFT_PDSSXML(ThermoPhase *tp, int spindex, 
					XML_Node& phaseNode, string id) {
    initThermo();
  }

   
  /**
   * constructHKFT_PDSSFile():
   *
   * Initialization of a Debye-Huckel phase using an
   * xml file.
   *
   * This routine is a precursor to initThermo(XML_Node*)
   * routine, which does most of the work.
   *
   * @param infile XML file containing the description of the
   *        phase
   *
   * @param id  Optional parameter identifying the name of the
   *            phase. If none is given, the first XML
   *            phase element will be used.
   */
  void HKFT_PDSS::constructHKFT_PDSSFile(ThermoPhase *tp, int spindex,
					 string inputFile, string id) {

    if (inputFile.size() == 0) {
      throw CanteraError("WaterTp::initThermo",
			 "input file is null");
    }
    string path = findInputFile(inputFile);
    ifstream fin(path.c_str());
    if (!fin) {
      throw CanteraError("HKFT_PDSS::initThermo","could not open "
			 +path+" for reading.");
    }
    /*
     * The phase object automatically constructs an XML object.
     * Use this object to store information.
     */

    XML_Node *fxml = new XML_Node();
    fxml->build(fin);
    XML_Node *fxml_phase = findXMLPhase(fxml, id);
    if (!fxml_phase) {
      throw CanteraError("HKFT_PDSS::initThermo",
			 "ERROR: Can not find phase named " +
			 id + " in file named " + inputFile);
    }	
    constructHKFT_PDSSXML(tp, spindex, *fxml_phase, id);
    delete fxml;
  }

  void HKFT_PDSS::
  initThermoXML(XML_Node& phaseNode, string id) {
    initThermo();
  }

  void HKFT_PDSS::initThermo() {
  }

  void HKFT_PDSS::
  setParametersFromXML(const XML_Node& eosdata) {
  }

  /**
   * Return the molar enthalpy in units of J kmol-1
   */
  doublereal HKFT_PDSS::
  enthalpy_mole() const {
    throw CanteraError("HKFT_PDSS::enthalpy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the internal energy in mks units of
   * J kmol-1 
   */
  doublereal HKFT_PDSS::
  intEnergy_mole() const {
    throw CanteraError("HKFT_PDSS::enthalpy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the entropy in mks units of 
   * J kmol-1 K-1
   */
  doublereal HKFT_PDSS::
  entropy_mole() const {

    throw CanteraError("HKFT_PDSS::entropy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the Gibbs free energy in mks units of
   * J kmol-1 K-1.
   */
  doublereal HKFT_PDSS::
  gibbs_mole() const {
    throw CanteraError("HKFT_PDSS::gibbs_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the constant pressure heat capacity
   * in mks units of J kmol-1 K-1
   */
  doublereal HKFT_PDSS::
  cp_mole() const {
    throw CanteraError("HKFT_PDSS::cp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the constant volume heat capacity
   * in mks units of J kmol-1 K-1
   */
  doublereal HKFT_PDSS::
  cv_mole() const {
    throw CanteraError("HKFT_PDSS::cv_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Return the difference in enthalpy between current p
   * and ref p0, in mks units of
   * in units of J kmol-1
   */
  doublereal HKFT_PDSS::
  enthalpyDelp_mole() const {
    throw CanteraError("HKFT_PDSS::enthalpy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate difference in the internal energy between current p
   * and ref p0, in mks units of
   * J kmol-1 
   */
  doublereal HKFT_PDSS::
  intEnergyDelp_mole() const {
    throw CanteraError("HKFT_PDSS::enthalpyDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   *  Return the difference in entropy between current p
   * and ref p0, in mks units of
   * J kmol-1 K-1
   */
  doublereal HKFT_PDSS::
  entropyDelp_mole() const {

    throw CanteraError("HKFT_PDSS::entropyDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the difference in Gibbs free energy between current p and
   * the ref p0, in mks units of
   * J kmol-1 K-1.
   */
  doublereal HKFT_PDSS::
  gibbsDelp_mole() const {
    throw CanteraError("HKFT_PDSS::gibbsDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the difference in the constant pressure heat capacity
   * between the current p and the ref p0,
   * in mks units of J kmol-1 K-1
   */
  doublereal HKFT_PDSS::
  cpDelp_mole() const {
    throw CanteraError("HKFT_PDSS::cpDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the difference in constant volume heat capacity
   * between the current p and the ref p0
   * in mks units of J kmol-1 K-1
   */
  doublereal HKFT_PDSS::
  cvDelp_mole() const {
    throw CanteraError("HKFT_PDSS::cvDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the pressure (Pascals), given the temperature and density
   *  Temperature: kelvin
   *  rho: density in kg m-3
   */
  doublereal HKFT_PDSS::
  pressure() const {
    throw CanteraError("HKFT_PDSS::pressure()", "unimplemented");
    return (0.0);
  }
        
  void HKFT_PDSS::
  setPressure(doublereal p) {
    throw CanteraError("HKFT_PDSS::pressure()", "unimplemented");
  }
 

  /// critical temperature 
  doublereal HKFT_PDSS::critTemperature() const { 
    throw CanteraError("HKFT_PDSS::critTemperature()", "unimplemented");
    return (0.0);
  }
        
  /// critical pressure
  doublereal HKFT_PDSS::critPressure() const {
    throw CanteraError("HKFT_PDSS::critPressure()", "unimplemented");
    return (0.0);
  }
        
  /// critical density
  doublereal HKFT_PDSS::critDensity() const {
    throw CanteraError("HKFT_PDSS::critDensity()", "unimplemented");
    return (0.0);
  }
        
  void HKFT_PDSS::setDensity(double dens) {
    m_dens = dens;
  }

  double HKFT_PDSS::density() const {
    return m_dens;
  }

  double HKFT_PDSS::temperature() const {
    return m_temp;
  }
 
  void HKFT_PDSS::setTemperature(double temp) {
    m_temp = temp;
  }

  doublereal HKFT_PDSS::molecularWeight() const {
    return m_mw;
  }
  void HKFT_PDSS::setMolecularWeight(double mw) {
    m_mw = mw;
  }

  void HKFT_PDSS::setState_TP(double temp, double pres) {
    throw CanteraError("HKFT_PDSS::setState_TP()", "unimplemented");
  }

  /// saturation pressure
  doublereal HKFT_PDSS::satPressure(doublereal t){
    throw CanteraError("HKFT_PDSS::satPressure()", "unimplemented");
    return (0.0);
  }
    

}
