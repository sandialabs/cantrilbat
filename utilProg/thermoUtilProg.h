
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
 */
struct Thermo_NASA {
    double Hf298_s;
    double S298_s;
    double tlow_;
    double thigh_;
    double pref_;
    mutable double tt[6];
    double m_coeffs[7];

    Thermo_NASA() :
        Hf298_s(0.0),
        S298_s(0.0),
        tlow_(5.0),
        thigh_(1000.),
        pref_(100000.)
    {
        for (size_t i = 0; i < 7; i++) {
            m_coeffs[i] = 0.0;
        }
    }

    Thermo_NASA(double tlow, double thigh, double pref, const double* XMLcoeffs) :
        Hf298_s(0.0),
        S298_s(0.0),
        tlow_(tlow),
        thigh_(thigh),
        pref_(pref)
    {
        m_coeffs[0] = XMLcoeffs[5];
        m_coeffs[1] = XMLcoeffs[6];
        for (size_t i = 0; i < 5; i++) {
            m_coeffs[2+i] = XMLcoeffs[i];
        }
        Hf298_s = h(298.15);
        S298_s = s(298.15);
    }

private:
    //! update the temperature polynomial
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
    double cp(double T) const
    {
        updateTemperaturePoly(T);
        double ct0 = m_coeffs[2];          // a0
        double ct1 = m_coeffs[3]*tt[0];    // a1 * T
        double ct2 = m_coeffs[4]*tt[1];    // a2 * T^2
        double ct3 = m_coeffs[5]*tt[2];    // a3 * T^3
        double ct4 = m_coeffs[6]*tt[3];    // a4 * T^4

        double cpR= ct0 + ct1 + ct2 + ct3 + ct4;
        return cpR * Zuzax::GasConstant * 1.0E-3;
    }

    //! Returns the enthalpy in units of kJ/gmol at the specified temperature
    /*!
     *  @param[in]           T                   Temperature (kelvin)
     *
     *  @return                                  Returns the enthalpy in units of kJ/gmol
     */
    double h(double T) const
    {
        updateTemperaturePoly(T);
        double ct0 = m_coeffs[2];          // a0
        double ct1 = m_coeffs[3]*tt[0];    // a1 * T
        double ct2 = m_coeffs[4]*tt[1];    // a2 * T^2
        double ct3 = m_coeffs[5]*tt[2];    // a3 * T^3
        double ct4 = m_coeffs[6]*tt[3];    // a4 * T^4

        double hRT = ct0 + 0.5*ct1 + ct2 / 3.0 + 0.25*ct3 + 0.2*ct4 + m_coeffs[0]*tt[4];     // last term is a5/T

        return hRT * T * Zuzax::GasConstant * 1.0E-6;
    }

    //! Returns the entropy in units of J/gmol/K at the specified temperature
    /*!
     *  @param[in]           T                   Temperature (kelvin)
     *
     *  @return                                  Returns the entropy in units of J/gmol /K
     */
    double s(double T) const
    {
        updateTemperaturePoly(T);
        double ct0 = m_coeffs[2];          // a0
        double ct1 = m_coeffs[3]*tt[0];    // a1 * T
        double ct2 = m_coeffs[4]*tt[1];    // a2 * T^2
        double ct3 = m_coeffs[5]*tt[2];    // a3 * T^3
        double ct4 = m_coeffs[6]*tt[3];    // a4 * T^4

        double sR = ct0*tt[5] + ct1 + 0.5*ct2 + ct3/3.0 +0.25*ct4 + m_coeffs[1];          // last term is a6

        return sR * Zuzax::GasConstant * 1.0E-3;
    }

    //! adjust the value of H given in kJ/gmol at a particular temperature
    /*!
     *   @param[in]
     */
    void adjustH(double T, double Hinput);

    void adjustS(double T, double Sinput);

    void printThermoBlock(int n);
};
//==================================================================================================================================

struct nasaMatch {
    double tempMatch;
    Thermo_NASA* nlow;
    Thermo_NASA* nhigh;

    nasaMatch(double T, Thermo_NASA* nl, Thermo_NASA* nh) :
        tempMatch(T),
        nlow(nl),
        nhigh(nh)
    {
    }
};

struct Regions_Nasa {
    std::vector<Thermo_NASA*> NASA_list;

    std::vector<nasaMatch> match_list;



};


// -------------------------------------------------------


