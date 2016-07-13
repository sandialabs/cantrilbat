/**
 *  @file HKFT_PDSS.h
 *
 * Declares class PDSS pressure dependent standard state
 * for a single species
 */

/*  $Author: hkmoffa $
 *  $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
 *  $Revision: 508 $
 *
 * 
 */

#ifndef CT_HKFT_PDSS_H
#define CT_HKFT_PDSS_H
#include "ct_defs.h"

class XML_Node;
class ThermoPhase;

class WaterPropsIAPWS;
#include "cantera/thermo/PDSS.h"

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif 
{


  /**
   * Class for pressure dependent standard states corresponding to 
   * ionic solutes in electrolyte water.
   *
   * NOTE: This is largely not done or not complete. 
   *
   */
  class HKFT_PDSS : public PDSS {

  public:

    /**
     * Basic list of constructors and duplicators
     */
    HKFT_PDSS(ThermoPhase *tp, int spindex);
    HKFT_PDSS(const HKFT_PDSS &b);
    HKFT_PDSS& operator=(const HKFT_PDSS&b);
    HKFT_PDSS(ThermoPhase *tp, int spindex,
	      std::string inputFile, std::string id = "");
    HKFT_PDSS(ThermoPhase *tp, int spindex, 
	      XML_Node& phaseRef, std::string id = "");
    virtual ~HKFT_PDSS();
        
    /**
     *   
     * @name  Utilities  
     * @{
     */
    virtual int pdssType() const { return -1; }

    /**
     * @} 
     * @name  Molar Thermodynamic Properties of the Solution --------------
     * @{
     */
    virtual doublereal enthalpy_mole() const;
    virtual doublereal intEnergy_mole() const;
    virtual doublereal entropy_mole() const;
    virtual doublereal gibbs_mole() const;
    virtual doublereal cp_mole() const;
    virtual doublereal cv_mole() const;

    /*
     * Get the difference in the standard state thermodynamic properties
     * between the reference pressure, po, and the current pressure.
     */
    virtual doublereal enthalpyDelp_mole() const;
    virtual doublereal intEnergyDelp_mole() const;
    virtual doublereal entropyDelp_mole() const;
    virtual doublereal gibbsDelp_mole() const;
    virtual doublereal cpDelp_mole() const;
    virtual doublereal cvDelp_mole() const;

    //@}
    /// @name Mechanical Equation of State Properties ---------------------
    //@{

    virtual doublereal pressure() const;
    virtual void setPressure(doublereal p);

    //@}
    /// @name  Partial Molar Properties of the Solution -----------------
    //@{

    virtual void getChemPotentials(doublereal* mu) const {
      mu[0] = gibbs_mole();
    }

    //@}
    /// @name  Properties of the Standard State of the Species
    //          in the Solution --
    //@{
    

    /// critical temperature 
    virtual doublereal critTemperature() const;
 
    /// critical pressure
    virtual doublereal critPressure() const;
        
    /// critical density
    virtual doublereal critDensity() const;
        
    /// saturation temperature
    //virtual doublereal satTemperature(doublereal p) const;
        
    

    /// saturation pressure
    virtual doublereal satPressure(doublereal t);
    
    virtual void setDensity(double dens);
    double density() const;
    virtual void setTemperature(double temp);
    double temperature() const;
    virtual void setState_TP(double temp, double pres);

    doublereal molecularWeight() const;
    void setMolecularWeight(double mw);
    
    void constructHKFT_PDSS(ThermoPhase *tp, int spindex);
    void constructHKFT_PDSSFile(ThermoPhase *tp, int spindex, 
				   std::string inputFile, std::string id);
    void constructHKFT_PDSSXML(ThermoPhase *tp, int spindex, 
				  XML_Node& phaseNode, std::string id);
    virtual void initThermoXML(XML_Node& eosdata, std::string id);
    virtual void initThermo();
    virtual void setParametersFromXML(const XML_Node& eosdata);

  protected:


  };

}

#endif



