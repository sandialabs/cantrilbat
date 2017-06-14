
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
// -------------------------------------------------------


