
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

    Thermo_Shomate():
         a_s(0.0), b_s(0.0), c_s(0.0), d_s(0.0), e_s(0.0), f_s(0.0), g_s(0.0) { }

    Thermo_Shomate(double a, double b, double c, double d, double e, double f = 0.0, double g = 0.0) :
         a_s(a), b_s(b), c_s(c), d_s(d), e_s(e), f_s(f), g_s(g) { }

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

    void convert_S298_Barin_ToShomate(double S298_cals);

    void convert_DeltaH_Barin_ToShomate(const Thermo_Shomate& solid, double T_melt, double deltaH_kcals);

    void convert_DeltaS_Barin_ToShomate(const Thermo_Shomate& solid, double T_melt, double deltaS_cals);

    void convert_PhaseChange_Barin_ToShomate(const Thermo_Shomate& solid, double T_melt, double deltaH_kcals,
                                             double deltaS_cals);

    void convert_H_Barin_ToShomate(double T, double H_kcals);

    void convert_S_Barin_ToShomate(double T, double S_cals);

    double h(double T) const;

    double s(double T) const;

    double cp(double T) const;

    double g(double T) const;

    void printThermoBlock( int n, double tmax, double tmin = 250.);
};
//==================================================================================================================================
//! Thermodynamics polynomial for NASA
/*!
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

      Thermo_NASA(double tlow, double thigh, double pref, const double* coeffs) :
          Hf298_s(0.0),
          S298_s(0.0),
          tlow_(tlow),
          pref_(pref)
     {
        for (size_t i = 0; i < 7; i++) {
            m_coeffs[i] = coeffs[i];
        }
        Hf298_s = h(298.15);
        S298_s = s(298.15);
     }

     void updateTemperaturePoly(double T) const {
        tt[0] = T;
        tt[1] = T * T;
        tt[2] = tt[1] * T;
        tt[3] = tt[2] * T;
        tt[4] = 1.0 / T;
        tt[5] = std::log(T);
    }

    double cp(double T) const {
        updateTemperaturePoly(T);
        double ct0 = m_coeffs[2];          // a0
        double ct1 = m_coeffs[3]*tt[0];    // a1 * T
        double ct2 = m_coeffs[4]*tt[1];    // a2 * T^2
        double ct3 = m_coeffs[5]*tt[2];    // a3 * T^3
        double ct4 = m_coeffs[6]*tt[3];    // a4 * T^4

        double cpR= ct0 + ct1 + ct2 + ct3 + ct4;
        return cpR * Zuzax::GasConstant;
    }

    double h(double T) const {
        updateTemperaturePoly(T);
        double ct0 = m_coeffs[2];          // a0
        double ct1 = m_coeffs[3]*tt[0];    // a1 * T
        double ct2 = m_coeffs[4]*tt[1];    // a2 * T^2
        double ct3 = m_coeffs[5]*tt[2];    // a3 * T^3
        double ct4 = m_coeffs[6]*tt[3];    // a4 * T^4

        double hRT = ct0 + 0.5*ct1 + ct2 / 3.0 + 0.25*ct3 + 0.2*ct4 + m_coeffs[0]*tt[4];     // last term is a5/T

        return hRT * T * Zuzax::GasConstant;
    }

    double s(double T) const {
        updateTemperaturePoly(T);
        double ct0 = m_coeffs[2];          // a0
        double ct1 = m_coeffs[3]*tt[0];    // a1 * T
        double ct2 = m_coeffs[4]*tt[1];    // a2 * T^2
        double ct3 = m_coeffs[5]*tt[2];    // a3 * T^3
        double ct4 = m_coeffs[6]*tt[3];    // a4 * T^4

        double sR = ct0*tt[5] + ct1 + 0.5*ct2 + ct3/3.0 +0.25*ct4 + m_coeffs[1];          // last term is a6
  
        return sR * Zuzax::GasConstant;
     }

     //! adjust the value of H given in kJ/gmol at a particular temperature
     /*!
      *   @param[in]
      */
     void adjustH(double T, double Hinput);

     void adjustS(double T, double Sinput);

    void printThermoBlock( int n );
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

struct Regions_Nasa
{
    std::vector<Thermo_NASA*> NASA_list;

    std::vector<nasaMatch> match_list;

    

};


// -------------------------------------------------------


