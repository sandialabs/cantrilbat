
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
static void ind(size_t n)
{ 
     for (size_t i = 0; i < n; i++) {
          printf(" ");
     }
}
//==================================================================================================================================

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
    *
    */
   double value(double T) 
   {
      double value_cals = A + B * 1.0E-3 * T + C * 1.0E5 / (T * T) + D * 1.0E-6 * T * T;
      double value_joules = value_cals * 4.184;
      return value_joules;
   }
};
// -------------------------------------------------------

//==================================================================================================================================
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
void Thermo_Shomate::convertCpBarinToShomate(const Cp_Barin& cpb)
{
    a_s = cpb.A * 4.184;
    b_s = cpb.B * 4.184;
    c_s = cpb.D * 4.184;
    d_s = 0.0 ;
    e_s = cpb.C / 10.0 * 4.184;
}
//==================================================================================================================================
void Thermo_Shomate::convert_Hf_Barin_ToShomate (double Hf_kcals)
{
   //H= A t + \frac{B t^2}{2} + \frac{C t^3}{3} + \frac{D t^4}{4}  - \frac{E}{t}  + F.
  double t = 298.15 / 1000.;
  double HfmF =  a_s * t + b_s * t * t / 2.0 + c_s * t * t * t / 3.0 + d_s * t * t * t * t / 4.0  - e_s / t;
  Hf298_s = Hf_kcals  * 4.184;
  f_s = Hf298_s - HfmF;
}
//==================================================================================================================================
void Thermo_Shomate::set_S298(double S298)
{
   double t = 298.15 / 1000.;
   //  \tilde{s}^0(T) = A\ln t + B t + \frac{C t^2}{2} + \frac{D t^3}{3} - \frac{E}{2t^2}  + G.
   double SmG = a_s * log(t) +  b_s * t + c_s * t * t / 2.0 + d_s * t * t * t / 3.0  - e_s / (2.0 * t * t);
   S298_s = S298;
   g_s = S298_s - SmG;
}
//==================================================================================================================================
void Thermo_Shomate::convert_S298_Barin_ToShomate(double S298_cals)
{
   double t = 298.15 / 1000.;
   //  \tilde{s}^0(T) = A\ln t + B t + \frac{C t^2}{2} + \frac{D t^3}{3} - \frac{E}{2t^2}  + G.
   double SmG = a_s * log(t) +  b_s * t + c_s * t * t / 2.0 + d_s * t * t * t / 3.0  - e_s / (2.0 * t * t);
   S298_s = S298_cals * 4.184;
   g_s = S298_s - SmG;
}

//==================================================================================================================================
void Thermo_Shomate::convert_H_Barin_ToShomate(double T, double H_kcals)
{
   //H= A t + \frac{B t^2}{2} + \frac{C t^3}{3} + \frac{D t^4}{4}  - \frac{E}{t}  + F.

   double H_kJ = H_kcals * 4.184;

   double t = T / 1000.;
   double HfmF =  a_s * t + b_s * t * t / 2.0 + c_s * t * t * t / 3.0 + d_s * t * t * t * t / 4.0  - e_s / t;

   f_s = H_kJ - HfmF;
   Hf298_s = h(298.15);
}
//==================================================================================================================================
void Thermo_Shomate::convert_S_Barin_ToShomate(double T, double S_cals)
{
   //H= A t + \frac{B t^2}{2} + \frac{C t^3}{3} + \frac{D t^4}{4}  - \frac{E}{t}  + F.

   double S_joules = S_cals * 4.184;

   double t = T / 1000.;
   double SmG = a_s * log(t) +  b_s * t + c_s * t * t / 2.0 + d_s * t * t * t / 3.0  - e_s / (2.0 * t * t);

   g_s = S_joules - SmG;
   S298_s = s(298.15);
}
//==================================================================================================================================
void Thermo_Shomate::convert_DeltaH_Barin_ToShomate(const Thermo_Shomate& solid, double T_melt, double deltaH_kcals)
{
   //H= A t + \frac{B t^2}{2} + \frac{C t^3}{3} + \frac{D t^4}{4}  - \frac{E}{t}  + F.
   double H_tm_solid = solid.h(T_melt);

   double H_tm_liquid = H_tm_solid + deltaH_kcals * 4.184;

   double t = T_melt / 1000.;
   double HfmF =  a_s * t + b_s * t * t / 2.0 + c_s * t * t * t / 3.0 + d_s * t * t * t * t / 4.0  - e_s / t;

   f_s = H_tm_liquid - HfmF;
   Hf298_s = h(298.15);
}
//==================================================================================================================================
void Thermo_Shomate::convert_DeltaS_Barin_ToShomate(const Thermo_Shomate& solid, double T_melt, double deltaS_cals)
{
   //H= A t + \frac{B t^2}{2} + \frac{C t^3}{3} + \frac{D t^4}{4}  - \frac{E}{t}  + F.
   double t = T_melt / 1000.;
   double S_tm_solid = solid.s(T_melt);
   double S_tm_liquid = S_tm_solid + deltaS_cals * 4.184;

   double SmG = a_s * log(t) +  b_s * t + c_s * t * t / 2.0 + d_s * t * t * t / 3.0  - e_s / (2.0 * t * t);
   g_s = S_tm_liquid - SmG;

   S298_s = s(298.15);
}
//==================================================================================================================================
//! Do a phase change
/*!
 *
 */
void Thermo_Shomate::convert_PhaseChange_Barin_ToShomate(const Thermo_Shomate& solid, double T_melt, double deltaH_kcals,
                                                         double deltaS_cals)
{
   double deltaG_kcals = deltaH_kcals - T_melt * deltaS_cals / 1000.;
   double deltaS_calsm = deltaS_cals;
   if (deltaG_kcals > 1.0E-8) {
       if (deltaG_kcals > 0.1) {
           printf("deltaH and deltaS are not compatible\n");
           exit(-1);
       } else {
           deltaS_calsm = deltaH_kcals * 1000. / T_melt;
       }
       printf("deltaS changed from %g to %g to ensure deltaG is identically zero at T = %g\n", deltaS_cals, deltaS_calsm, T_melt);
   }
   convert_DeltaH_Barin_ToShomate(solid, T_melt, deltaH_kcals);
   convert_DeltaS_Barin_ToShomate(solid, T_melt, deltaS_calsm);
}
//==================================================================================================================================
// h in kJ/gmol/K
double Thermo_Shomate::h(double T) const
{
   double t = T / 1000.;
   double H =  a_s * t + b_s * t * t / 2.0 + c_s * t * t * t / 3.0 + d_s * t * t * t * t / 4.0  - e_s / t + f_s;
   return H;
}
//==================================================================================================================================
// s in j/gmol/K
double Thermo_Shomate::s(double T) const
{
   double t = T / 1000.;
   double S = a_s * log(t) +  b_s * t + c_s * t * t / 2.0 + d_s * t * t * t / 3.0  - e_s / (2.0 * t * t) + g_s;
   return S;
}
//==================================================================================================================================
// cp in j/gmol/K
double Thermo_Shomate::cp(double T) const
{
   double t = T / 1000.;
   double cp = a_s +  b_s * t + c_s * t * t  + d_s * t * t * t  + e_s / (t * t);
   return cp;
}
//==================================================================================================================================
// g in kJ/gmol/K
double Thermo_Shomate::g(double T) const
{
   double t = T / 1000.;
   double H = h(T);
   double S = s(T);
   double G = H - S * t;
   return G;
}
//==================================================================================================================================
void Thermo_Shomate::printThermoBlock(int n, double tmax, double tmin )
{
     ind(n); printf("<!-- thermo From Barin-Knacke Translation, Hf298 = %g kJ/gmol S298= %g J/gmol/K -->\n", Hf298_s, S298_s);
     ind(n); printf("<thermo>\n");
     ind(n); printf("  <Shomate Pref=\"1 bar\" Tmax=\"%g\" Tmin=\"%g\">\n", tmax, tmin);
     ind(n); printf("     <floatArray size=\"7\">\n");
     int m = n + 6;
     ind(m); printf("%-11.8g,\n", a_s);
     ind(m); printf("%-11.8g,\n", b_s);
     ind(m); printf("%-11.8g,\n", c_s);
     ind(m); printf("%-11.8g,\n", d_s);
     ind(m); printf("%-11.8g,\n", e_s);
     ind(m); printf("%-11.8g,\n", f_s);
     ind(m); printf("%-11.8g \n", g_s);
     ind(n); printf("    </floatArray>\n");
     ind(n); printf("  </Shomate>\n");
     ind(n); printf("</thermo>\n");
}
// -------------------------------------------------------
int main () {


  Cp_Barin Cp_Cgraph_SolA(0.026, 9.307, -0.354, -4.155);   //298 - 737
  Cp_Barin Cp_Cgraph_SolB(5.841, 0.104, -7.559, 0.0);

  Thermo_Shomate ts_Cgraph_SolA;

  ts_Cgraph_SolA.convertCpBarinToShomate(Cp_Cgraph_SolA);

  ts_Cgraph_SolA.convert_Hf_Barin_ToShomate(0.0);

  //ts_Cgraph_SolA.convert_S298_Barin_ToShomate(1.372);
  // Set it to the exact value that we are using for the element Ca entropy contribution
  ts_Cgraph_SolA.set_S298(5.740);


  double dens = 2.15;  // gm / cm3;
  double mw = 12.011; 
  double mv = mw / dens * 1.0E-3;
  
  printf("      <!-- species Cgraph_SolA   -->\n");
 
  printf("      <species name=\"Cgraph_SolA\">\n");
  printf("        <atomArray> C:1 </atomArray>\n");
  ts_Cgraph_SolA.printThermoBlock(8, 1100 , 200.);
  printf("        <standardState  model=\"constant_incompressible\">\n");
  printf("          <!-- molar volume from C(s) based on density of %g gm/cm3 --> \n", dens);
  printf("          <molarVolume units=\"m3/kmol\">  %-11.8g </molarVolume>\n", mv);
  printf("        </standardState>\n");
  printf("      </species>\n");
  printf("                \n");

  //---------------------------------------------------------------------------------------------------------------

  Thermo_Shomate ts_Cgraph_SolB;

  ts_Cgraph_SolB.convertCpBarinToShomate(Cp_Cgraph_SolB);
  ts_Cgraph_SolB.convert_PhaseChange_Barin_ToShomate(ts_Cgraph_SolA, 1100., 0.00, 0.0);

  printf("      <!-- species Cgraph_SolB  -->\n");

  printf("      <species name=\"Cgraph_SolB\">\n");
  printf("       <atomArray> C:1  </atomArray>\n");
  ts_Cgraph_SolB.printThermoBlock(8, 4073, 1100.);
  printf("       <standardState  model=\"constant_incompressible\">\n");
  printf("          <!-- molar volume from C(s) --> \n");
  printf("          <molarVolume units=\"m3/kmol\">  %-11.8g </molarVolume>\n", mv);
  printf("       </standardState>\n");
  printf("      </species>\n");
  printf("                \n");

  //---------------------------------------------------------------------------------------------------------------


  return 0;

}


