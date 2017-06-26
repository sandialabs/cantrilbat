
// Example of taking a Barin-Knack Table entry and creating a shomate polynomial for use in Zuzax

/*
 *  Theory
 * -------------------------
 *
 *  Barin expresses his units for heat capacity in cal/gmol/K
 *
 *  Cp (cal/gmol/K) = A + B 1.0E-3 T + C 1.0E5 T^-2 + D 1.0E-6 T^2
 *
 *
 *
 */
#include <stdio.h>
#include <cmath>
#include <cstdlib>

#include "cantera/base/ct_defs.h"

//==================================================================================================================================
//! print an indent
/*!
 *  @param[in]               n                   number of spaces to indent
 */
extern void ind(size_t n);

// -------------------------------------------------------
//! Heat Capacity in Barin format
/*!
 *  Hist inputs are in the form of 4 parameters listed as:
 *     A   B   C   D
 *
 *  The input here is the same as the table input
 */
struct Cp_Barin {
    double A ;
    double B ;
    double C ;
    double D ;

    Cp_Barin(double a, double b, double c, double d) :
        A(a),
        B(b),
        C(c),
        D(d)
    {
    }

    //! Calculate the value of Cp in J/gmol/K
    /*!
     *  @param[in]            T                   temperature
     *
     *  @return                                   returns the heat capacity in J/gmol/K
     */
    double value(double T)
    {
        double value_cals = A + B * 1.0E-3 * T + C * 1.0E5 / (T * T) + D * 1.0E-6 * T * T;
        double value_joules = value_cals * 4.184;
        return value_joules;
    }
};
//==================================================================================================================================
//! Thermodynamics polynomial for Shomate
/*!
 * Seven coefficients \f$(A,\dots,G)\f$ are used to represent
 * \f$ c_p^0(T)\f$, \f$ h^0(T)\f$, and \f$ s^0(T) \f$ as polynomials in the temperature, \f$ T \f$ :
 *
 * \f[
 *     \tilde{c}_p^0(T) = A + B t + C t^2 + D t^3 + \frac{E}{t^2}
 * \f]
 *
 * \f[
 *     \tilde{h}^0(T) = A t + \frac{B t^2}{2} + \frac{C t^3}{3} + \frac{D t^4}{4}  - \frac{E}{t}  + F.
 * \f]
 *
 * \f[
 *     \tilde{s}^0(T) = A\ln t + B t + \frac{C t^2}{2} + \frac{D t^3}{3} - \frac{E}{2t^2}  + G.
 * \f]
 *
 * In the above expressions, the thermodynamic polynomials are expressed
 * in dimensional units, but the temperature,\f$ t \f$, is divided by 1000. The
 * following dimensions are assumed in the above expressions:
 *
 *    - \f$ \tilde{c}_p^0(T)\f$ = Heat Capacity (J/gmol*K)
 *    - \f$ \tilde{h}^0(T) \f$ = standard Enthalpy (kJ/gmol)
 *    - \f$ \tilde{s}^0(T) \f$= standard Entropy (J/gmol*K)
 *    - \f$ t \f$= temperature (K) / 1000.
 *
 *  For more information about Shomate polynomials, see the NIST website,  http://webbook.nist.gov/
 *
*   @verbatim
 *    <Shomate Pref="1 atm" Tmax="   623.15" Tmin="   298.00">
         <floatArray size="7">
             A, B, C, D, E, F
          </floatArray>
       </Shomate>
      @endverbatim
 *
 *
 */
struct Thermo_Shomate {
    double a_s;
    double b_s;
    double c_s;
    double d_s;
    double e_s;
    double f_s;
    double g_s;
    double Hf298_s;
    double S298_s;

    //! Default constructor
    Thermo_Shomate():
        a_s(0.0), b_s(0.0), c_s(0.0), d_s(0.0), e_s(0.0), f_s(0.0), g_s(0.0) { }

    //! Constructor for Thermo_Shomate
    /*!
     *  @param[in]           a
     *  @param[in]           b
     *  @param[in]           c
     *  @param[in]           d
     *  @param[in]           e
     *  @param[in]           f
     *  @param[in]           g
     */
    Thermo_Shomate(double a, double b, double c, double d, double e, double f = 0.0, double g = 0.0) :
        a_s(a), b_s(b), c_s(c), d_s(d), e_s(e), f_s(f), g_s(g) { }

    //! Takes the heat capacity polynomials given in the Barin & Knacke book and fills in the shomate polynomial values.
    /*!
     *  @param[in]           cpb                 Reference to the Cp_Barin structure containing the polynomials
     */
    void convertCpBarinToShomate(const Cp_Barin& cpb);

    /*!
     *  @param[in]         Hf_Barin      heat of formation at 298.15 in kcal/gmol
     */
    void convert_Hf_Barin_ToShomate(double Hf_kcals);

    //! Set the S298 value
    /*!
     *  @param[in]           S298                Entropy at 298 (Joules / gmol /K)
     */
    void set_S298(double S298);

    //! Converts the Barin Entropy value (cals/gmol/K) at 298 K into a value for S298 that is used in the Shomate polynomial
    /*!
     *  This routine will calculate the value of g_s in the Shomate polynomial given the previously given heat capacity polynomial.
     *
     *  @param[in]           S298_cals          Value of the entropy at temperature 298.15 in cals/gmol/K
     */
    void convert_S298_Barin_ToShomate(double S298_cals);

    //! Converts the Barin delta Enthalpy value (kcals/mole) for a phase change into a value for Hf298 that is used in the Shomate polynomial
    /*!
     *  This routine will calculate the value of f_s in the Shomate polynomial given the previously given heat capacity polynomial
     *  and a previously computed complete Shomate polynomial for the lower temperature phase.
     *
     *  @param[in]           solid               Reference to the Thermo_Shomate solid object that is on the lower temperature side
     *                                           of the phase change.
     *  @param[in]           T                   Temperature in Kelvin for the phase change
     *  @param[in]           deltaH_kcals        Value of the delta enthalpy at temperature T in kcals/gmol
     */
    void convert_DeltaH_Barin_ToShomate(const Thermo_Shomate& solid, double T_melt, double deltaH_kcals);

    //! Converts the Barin delta Entropy value (kcals/mole) for a phase change into a value for S298 that is used in the Shomate polynomial
    /*!
     *  This routine will calculate the value of g_s in the Shomate polynomial given the previously given heat capacity polynomial
     *  and a previously computed complete Shomate polynomial for the lower temperature phase.
     *
     *  @param[in]           solid               Reference to the Thermo_Shomate solid object that is on the lower temperature side
     *                                           of the phase change.
     *  @param[in]           T                   Temperature in Kelvin for the phase change
     *  @param[in]           deltaS_cals         Value of the delta enthalpy at temperature T in cals/gmol/K
     */
    void convert_DeltaS_Barin_ToShomate(const Thermo_Shomate& solid, double T_melt, double deltaS_cals);

    //! Converts the Barin delta Enthalpy value (kcals/mole) and delta entropy value (cals/gmol/K for a phase change
    //! into a value for Hf298 and S298 that is used in the Shomate polynomial
    /*!
     *  This routine will calculate the value of f_s and g_s in the Shomate polynomial given the previously given heat capacity polynomial
     *  and a previously computed complete Shomate polynomial for the lower temperature phase.
     *  It also makes sure that the phase change occurs exactly at the phase change temperature by ensuring that
     *             deltaHs - T deltaS = 0
     *  deltas will be slightly adjusted to make this equation happen.
     *
     *  @param[in]           solid               Reference to the Thermo_Shomate solid object that is on the lower temperature side
     *                                           of the phase change.
     *  @param[in]           T                   Temperature in Kelvin for the phase change
     *  @param[in]           deltaH_kcals        Value of the delta enthalpy at temperature T in kcals/gmol
     */
    void convert_PhaseChange_Barin_ToShomate(const Thermo_Shomate& solid, double T_melt, double deltaH_kcals,
                                             double deltaS_cals);

    //! Converts the Barin Enthalpy value (kcals/mole) into a value for Hf298 that is used in the Shomate polynomial
    /*!
     *  This routine will calculate the value of f_s in the Shomate polynomial given the previously given heat capacity polynomial.
     *
     *  @param[in]           T                   Temperature in Kelvin for the input enthalpy value
     *  @param[in]           H_kcals             Value of the enthalpy at temperature T in kcals/gmol
     */
    void convert_H_Barin_ToShomate(double T, double H_kcals);

    //! Converts the Barin Entropy value (cals/gmol/K) into a value for S298 that is used in the Shomate polynomial
    /*!
     *  This routine will calculate the value of g_s in the Shomate polynomial given the previously given heat capacity polynomial.
     *
     *  @param[in]           T                   Temperature in Kelvin for the input enthalpy value
     *  @param[in]           S_cals              Value of the entropy at temperature T in cals/gmol/K
     */
    void convert_S_Barin_ToShomate(double T, double S_cals);

    //! Returns the enthalpy in units of kJ/gmol at a temperature
    /*!
     *  @param[in]           T                   Temperature (kelvin)
     *
     *  @return                                  Returns the enthalpy in units of kJ/gmol
     */
    double h(double T) const;

    //! Returns the entropy in units of J/gmol/K at a temperature
    /*!
     *  @param[in]           T                   Temperature (kelvin)
     *
     *  @return                                  Returns the entropy in units of J/gmol /K
     */
    double s(double T) const;

    //! Returns the constant pressure heat capacity in units of J/gmol/K at a temperature
    /*!
     *  @param[in]           T                   Temperature (kelvin)
     *
     *  @return                                  Returns the heat capacity in units of J/gmol /K
     */
    double cp(double T) const;

    //! Returns the gibbs free energy in units of kJ/gmol at a temperature
    /*!
     *  @param[in]           T                   Temperature (kelvin)
     *
     *  @return                                  Returns the gibbs free energy in units of kJ/gmol
     */
    double g(double T) const;

    //! Adjust the value of H given in kJ/gmol at a particular temperature
    /*!
     *  @param[in]          T                    Temperature (kelvin)
     *  @param[in]          Hinput               Value of the enthalpy to use at the temperature (kJ/gmol)
     *
     *  @return                                  Returns the net adjustment of H in kJ/gmol
     */
    double adjustH(double T, double Hinput);

    //! Prints out an XML thermo block for the parameterization
    /*!
     *  @param[in]           n                   Amount to indent each line
     *  @param[in]           tmax                Maximum temperature
     *  @param[in]           tmin                Minimum temperature
     */
    void printThermoBlock(int n, double tmax, double tmin = 250.);
};

//==================================================================================================================================
//! Thermodynamics polynomial for NASA
/*!
 *  Seven coefficients \f$(a_0,\dots,a_6)\f$ are used to represent \f$ c_p^0(T)\f$, \f$ h^0(T)\f$, and \f$ s^0(T) \f$
 *  as polynomials in \f$ T \f$ :
 *
 * \f[
 *          \frac{c_p(T)}{R} = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4
 * \f]
 * \f[
 *          \frac{h^0(T)}{RT} = a_0 + \frac{a_1}{2} T + \frac{a_2}{3} T^2 + \frac{a_3}{4} T^3 + \frac{a_4}{5} T^4  + \frac{a_5}{T}.
 * \f]
 * \f[
 *          \frac{s^0(T)}{R} = a_0\ln T + a_1 T + \frac{a_2}{2} T^2 + \frac{a_3}{3} T^3 + \frac{a_4}{4} T^4  + a_6.
 * \f]
 *
 *  Note assignment object and copy constructor cannot the default implementation, due to arrays. Need to fix.
 */
struct Thermo_NASA {

    //! Enthalpy at 298.15K and reference pressure (kJ/gmol)
    double Hf298_s;
    //! Entropy at 298.15K and reference pressure (J/gmol/K)
    double S298_s;
    //! Low temperature limit
    double tlow_;
    //! High temperature limit
    double thigh_;
    //! Reference pressure (bar)
    double pref_;
    //! Temporary temperature polynomials 
    mutable double tt[6];
    //! Coefficients of the representation. 
    /*!
     *  These are in the order that the NASA object within Zuzax uses. The order is changed
     *  from what is in the XML file.
     */
    double m_coeff[7];

    //! Default constructor
    Thermo_NASA(); 

    //! Main constructor for the NASA polynomial routine
    /*!
     *  @param[in]           tlow                Low temperature value of the polynomial
     *  @param[in]           thigh               high temperature value of the polynomial
     *  @param[in]           pref                Reference pressure for the representation
     *  @param[in]           XMLcoeffs           Vector of coefficients for the representation in  the order
     *                                           that they appear in the Zuzax XML files. Note they are transmuted
     *                                           into a different order within this object and within the NASA
     *                                           polynomial object within Zuzax.
     */
    Thermo_NASA(double tlow, double thigh, double pref, const double* XMLcoeffs);

private:
    //! Update the temperature polynomial
    /*!
     *  @param[in]          T                   Temperature in Kelvin
     */
    void updateTemperaturePoly(double T) const;

public:
    //! Returns the constant pressure heat capacity in units of J/gmol/K at the specified temperature
    /*!
     *  @param[in]           T                   Temperature (kelvin)
     *
     *  @return                                  Returns the heat capacity in units of J/gmol /K
     */
    double cp(double T) const;

    //! Returns the enthalpy in units of kJ/gmol at the specified temperature
    /*!
     *  @param[in]           T                   Temperature (kelvin)
     *
     *  @return                                  Returns the enthalpy in units of kJ/gmol
     */
    double h(double T) const;

    //! Returns the entropy in units of J/gmol/K at the specified temperature
    /*!
     *  @param[in]           T                   Temperature (kelvin)
     *
     *  @return                                  Returns the entropy in units of J/gmol /K
     */
    double s(double T) const;

    //! Returns the enthalpy difference at T from the 298.15K value
    /*!
     *  @param[in]           T                   Temperature
     *  @return                                  Returns the H-H298 in kJ/gmol
     */
    double deltaH(double T) const;

    //! Adjust the value of H given in kJ/gmol at a particular temperature
    /*!
     *  @param[in]          T                    Temperature (kelvin)
     *  @param[in]          Hinput               Value of the enthalpy to use at the temperature (kJ/gmol)
     *
     *  @return                                  Returns the net adjustment of H in kJ/gmol
     */
    double adjustH(double T, double Hinput);

    //! Adjust the value of S given in J/gmol/K at a particular temperature
    /*!
     *  @param[in]          T                    Temperature (kelvin)
     *  @param[in]          Sinput               Value of the entropy to use at the temperature (J/gmol/K)
     *
     *  @return                                  Returns the net adjustment of S in kJ/gmol/K
     */
    double adjustS(double T, double Sinput);

    //! Adjust the value of Cp given in J/gmol/K at a particular temperature
    /*!
     *  This is different than just adjusting the value of H or S. We change the coefficients for the value of H and S 
     *  after the adjust in Cp so that it gives the same value of H and S at T as the representation gave previously.
     *
     *  @param[in]          T                    Temperature to apply all of the adjustments (kelvin).
     *  @param[in]          Cpinput              Value of the heat capacity to use at the temperature (J/gmol/K)
     *
     *  @return                                  Returns the net adjustment of Cp in kJ/gmol/K
     */
    double adjustCp(double T, double Sinput);

    //! Print out a block of XML code representing the NASA representation
    /*! 
     *  @param[in]           n                   Minimum indentation of all of the lines
     *  @param[in]           p                   Precision in the writing of the block
     */
    void printNASABlock(int n, int p = 10);

    //! Print out a block of XML code representing the thermodynamic representation
    /*! 
     *  @param[in]           n                   Minimum indentation of all of the lines
     *  @param[in]           p                   Precision in the writing of the block
     */
    void printThermoBlock(int n, int p = 10);
};
//==================================================================================================================================
//! Estabilish a matching condition between two polynomials
/*!
 *     Note by design, this uses pointers. The temperature polynomials should be malloced and storred elsewhere.
 * 
 */
struct nasaMatch {
    double tempMatch;
    Thermo_NASA* nlow;
    Thermo_NASA* nhigh;

    //! Constructor for the nasaMatch structure
    /*!
     *  @param[in]           T                   Temperature of the matching condition
     *                                           This should also equal the high temp of the low polynomial, and the low temp of the
     *                                           high temperature polynomial
     *  @param[in]           nl                  Pointer to the low temperature NASA polynomial
     *  @param[in]           nh                  Pointer to the high temperature NASA polynomial
     *
     */
    nasaMatch(double T, Thermo_NASA* nl, Thermo_NASA* nh) :
        tempMatch(T),
        nlow(nl),
        nhigh(nh)
    {
        if (fabs(nlow->thigh_ - nhigh->tlow_) > 1.0E-10) {
            printf(" nasaMatch constructor: tlow and thigh are not set correctly: %25.15E %25.15E\n",
                     nlow->thigh_ , nhigh->tlow_);
            exit(-1);
        } 
    }

    //! Adjust the high temperature polynomial so that there is matching enthalpy condition at the  tempMatch temperature
    /*!
     *  @return                                  Returns the deltaH change in the high temperature value of H needed to
     *                                           accomplish the matching condition (kJ/gmol)
     */
    double adjustH_high();

    //! Adjust the high temperature polynomial so that there is matching entropy condition at the tempMatch temperature
    /*!
     *  @return                                  Returns the deltaS change in the high temperature value of S needed to
     *                                           accomplish the matching condition (J/gmol/K)
     */
    double adjustS_high();

    //! Adjust the high temperature polynomial so that there is matching Cp condition at the tempMatch temperature
    /*!
     *  @return                                  Returns the deltaCp change in the high temperature value of Cp needed to
     *                                           accomplish the matching condition (J/gmol/K)
     */
    double adjustCp_high();

    //! Adjust the low temperature polynomial so that there is matching enthalpy condition at the  tempMatch temperature
    /*!
     *  @return                                  Returns the deltaH change in the low temperature value of H needed to
     *                                           accomplish the matching condition (kJ/gmol)
     */
    double adjustH_low();

    //! Adjust the low temperature polynomial so that there is matching entropy condition at the tempMatch temperature
    /*!
     *  @return                                  Returns the deltaS change in the low temperature value of S needed to
     *                                           accomplish the matching condition (J/gmol/K)
     */
    double adjustS_low();

    //! Adjust the low temperature polynomial so that there is matching Cp condition at the tempMatch temperature
    /*!
     *  @return                                  Returns the deltaCp change in the low temperature value of Cp needed to
     *                                           accomplish the matching condition (J/gmol/K)
     */
    double adjustCp_low();

    //! Adjust the high and low temperature polynomial so that there is matching enthalpy condition at the  tempMatch temperature
    /*!
     *  This takes the high and low polynomial value of H, and uses the average value at tempMatch of H.
     *  Therefore both polynomials are changed.
     *
     *  @return                                  Returns the deltaH change in the high temperature value of H needed to
     *                                           accomplish the matching condition (kJ/gmol)
     */
    double adjustH_avg();

    //! Adjust the high temperature polynomial so that there is matching entropy condition at the tempMatch temperature
    /*!
     *  This takes the high and low polynomial value of S, and uses the average value at tempMatch of S.
     *  Therefore both polynomials are changed.
     *
     *  @return                                  Returns the deltaS change in the high temperature value of S needed to
     *                                           accomplish the matching condition (J/gmol/K)
     */
    double adjustS_avg();

    //! Adjust the high and low temperature polynomial so that there is matching Cp condition at the tempMatch temperature
    /*!
     *  This takes the high and low polynomial value of Cp, and uses the average value at tempMatch of Cp.
     *  Therefore both polynomials are changed.
     *
     *  @return                                  Returns the deltaCp change in the high temperature value of Cp needed to
     *                                           accomplish the matching condition (J/gmol/K)
     */
    double adjustCp_avg();
};
//==================================================================================================================================
//! Container for holding pointer to Multiple NASA polynomial interpolations representing the thermo of a single speices
/*!
 *   Designed for pointers
 *
 *  This could be much more sophisticated. However, it serves its purpose as is.
 */
class Regions_NASA {
public:
    //! Vector of thermo representations, must be ordered from low to high temperature regions
    std::vector<Thermo_NASA*> NASA_list;

    //! vector of matching temperature conditions
    std::vector<nasaMatch> match_list;

    //! Low value of the temperature
    double tlowReg_;

    //! High value of the temperatue
    double thighReg_;

    //! Added another temperature polynomial covering temperatures above the current region
    /*!
     *  Program checks that the temperature region is above the last region.
     *  Automatically adds a matchingNASA condition at the interface temperature.
     *
     *  @param[in]           tNasa               Pointer to the Thermo_NASA polynomial for the next higher temperature region
     */
    void addRegionPoly(Thermo_NASA* tNasa);

    //! Carry out matching conditions at the interface temperatures
    /*!
     *  Ensures continuity of H, S, and Cp at the interface temperatures.
     *  Goes from low to high temperature matching temperatures, changing the high temperature polynomial
     *
     *  @return                                   Returns the max enthalpy change
     */
    double adjust_high();

    //! Print out a block of XML code representing the thermodynamic representation
    /*! 
     *  @param[in]           n                   Minimum indentation of all of the lines
     *  @param[in]           p                   Precision in the writing of the block
     */
    void printThermoBlock(int n, int p);
};
//==================================================================================================================================

//==================================================================================================================================
//! Thermodynamics polynomial for NASA9
/*!
 *  Nine coefficients \f$(a_0,\dots,a_6)\f$ are used to represent \f$ c_p^0(T)\f$, \f$ h^0(T)\f$, and \f$ s^0(T) \f$
 *  as polynomials in \f$ T \f$ :
 *
 * \f[
 *     \frac{C_p^0(T)}{R} = a_0 T^{-2} + a_1 T^{-1} + a_2 + a_3 T + a_4 T^2 + a_5 T^3 + a_6 T^4
 * \f]
 *
 * \f[
 *     \frac{H^0(T)}{RT} = - a_0 T^{-2} + a_1 \frac{\ln T}{T} + a_2 + \frac{a_3}{2} T + \frac{a_4}{3} T^2  + \frac{a_5}{4} T^3 +
 *                          \frac{a_6}{5} T^4 + \frac{a_7}{T}
 * \f]
 *
 * \f[
 *    \frac{s^0(T)}{R} = - \frac{a_0}{2} T^{-2} - a_1 T^{-1} + a_2 \ln T + a_3 T + \frac{a_4}{2} T^2 + \frac{a_5}{3} T^3 
 *                       + \frac{a_6}{4} T^4 + a_8
 * \f]
 *
 *
 *  Note assignment object and copy constructor cannot the default implementation, due to arrays. Need to fix.
 */
struct Thermo_NASA9 {

    //! Enthalpy at 298.15K and reference pressure (kJ/gmol)
    double Hf298_s;
    //! Entropy at 298.15K and reference pressure (J/gmol/K)
    double S298_s;
    //! Low temperature limit
    double tlow_;
    //! High temperature limit
    double thigh_;
    //! Reference pressure (bar)
    double pref_;
    //! Temporary temperature polynomials 
    mutable double tt[7];
    //! Coefficients of the representation. 
    /*!
     *  These are in the order that the NASA object within Zuzax uses.
     */
    double m_coeff[9];

    //! Default constructor
    Thermo_NASA9(); 

    //! Main constructor for the NASA polynomial routine
    /*!
     *  @param[in]           tlow                Low temperature value of the polynomial
     *  @param[in]           thigh               high temperature value of the polynomial
     *  @param[in]           pref                Reference pressure for the representation
     *  @param[in]           XMLcoeffs           Vector of coefficients for the representation in  the order
     *                                           that they appear in the Zuzax XML files. Note they are transmuted
     *                                           into a different order within this object and within the NASA
     *                                           polynomial object within Zuzax.
     */
    Thermo_NASA9(double tlow, double thigh, double pref, const double* XMLcoeffs);

private:
    //! Update the temperature polynomial
    /*!
     *  @param[in]          T                   Temperature in Kelvin
     */
    void updateTemperaturePoly(double T) const;

public:
    //! Returns the constant pressure heat capacity in units of J/gmol/K at the specified temperature
    /*!
     *  @param[in]           T                   Temperature (kelvin)
     *
     *  @return                                  Returns the heat capacity in units of J/gmol /K
     */
    double cp(double T) const;

    //! Returns the enthalpy in units of kJ/gmol at the specified temperature
    /*!
     *  @param[in]           T                   Temperature (kelvin)
     *
     *  @return                                  Returns the enthalpy in units of kJ/gmol
     */
    double h(double T) const;

    //! Returns the entropy in units of J/gmol/K at the specified temperature
    /*!
     *  @param[in]           T                   Temperature (kelvin)
     *
     *  @return                                  Returns the entropy in units of J/gmol /K
     */
    double s(double T) const;

    //! Returns the enthalpy difference at T from the 298.15K value
    /*!
     *  @param[in]           T                   Temperature
     *  @return                                  Returns the H-H298 in kJ/gmol
     */
    double deltaH(double T) const;

    //! Adjust the value of H given in kJ/gmol at a particular temperature
    /*!
     *  @param[in]          T                    Temperature (kelvin)
     *  @param[in]          Hinput               Value of the enthalpy to use at the temperature (kJ/gmol)
     *
     *  @return                                  Returns the net adjustment of H in kJ/gmol
     */
    double adjustH(double T, double Hinput);

    //! Adjust the value of S given in J/gmol/K at a particular temperature
    /*!
     *  @param[in]          T                    Temperature (kelvin)
     *  @param[in]          Sinput               Value of the entropy to use at the temperature (J/gmol/K)
     *
     *  @return                                  Returns the net adjustment of S in kJ/gmol/K
     */
    double adjustS(double T, double Sinput);

    //! Adjust the value of Cp given in J/gmol/K at a particular temperature
    /*!
     *  This is different than just adjusting the value of H or S. We change the coefficients for the value of H and S 
     *  after the adjust in Cp so that it gives the same value of H and S at T as the representation gave previously.
     *
     *  @param[in]          T                    Temperature to apply all of the adjustments (kelvin).
     *  @param[in]          Cpinput              Value of the heat capacity to use at the temperature (J/gmol/K)
     *
     *  @return                                  Returns the net adjustment of Cp in kJ/gmol/K
     */
    double adjustCp(double T, double Sinput);

    //! Print out a block of XML code representing the NASA representation
    /*! 
     *  @param[in]           n                   Minimum indentation of all of the lines
     *  @param[in]           p                   Precision in the writing of the block
     */
    void printNASA9Block(int n, int p = 10);

    //! Print out a block of XML code representing the thermodynamic representation
    /*! 
     *  @param[in]           n                   Minimum indentation of all of the lines
     *  @param[in]           p                   Precision in the writing of the block
     */
    void printThermoBlock(int n, int p = 10);
};
//==================================================================================================================================

