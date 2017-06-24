
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

#include "thermoUtilProg.h"

//==================================================================================================================================
//! print an indent 
/*!
 *  @param[in]               n                   number of spaces to indent
 */
void ind(size_t n)
{ 
     for (size_t i = 0; i < n; i++) {
          printf(" ");
     }
}
//==================================================================================================================================
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
double Thermo_Shomate::adjustH(double T, double Hinput)
{
   double hnow = h(T);
   double H_kJ = Hinput;
   double t = T / 1000.;
   double HfmF =  a_s * t + b_s * t * t / 2.0 + c_s * t * t * t / 3.0 + d_s * t * t * t * t / 4.0  - e_s / t;
   f_s = H_kJ - HfmF;
   Hf298_s = h(298.15);
   double hnew = h(T);
   if (fabs(hnew - Hinput) > 1.0E-10) {
       printf("Thermo_Shomate::adjustH() ERROR: There is failure somewhere\n");
       exit(-1);
   }
   return (hnew - hnow);
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
//==================================================================================================================================
//==================================================================================================================================
//==================================================================================================================================
Thermo_NASA::Thermo_NASA() :
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
//==================================================================================================================================
Thermo_NASA::Thermo_NASA(double tlow, double thigh, double pref, const double* XMLcoeffs) :
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
//==================================================================================================================================
void Thermo_NASA::updateTemperaturePoly(double T) const 
{
    tt[0] = T;
    tt[1] = T * T;
    tt[2] = tt[1] * T;
    tt[3] = tt[2] * T;
    tt[4] = 1.0 / T;
    tt[5] = std::log(T);
}
//==================================================================================================================================
double Thermo_NASA::cp(double T) const
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
//==================================================================================================================================
double Thermo_NASA::h(double T) const
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
//==================================================================================================================================
double Thermo_NASA::s(double T) const
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
//==================================================================================================================================
double Thermo_NASA::deltaH(double T) const
{
    return h(T) - Hf298_s;
}
//==================================================================================================================================
double Thermo_NASA::adjustH(double T, double Hinput)
{
    double hnow = h(T);
    double HinputRT = Hinput * 1.0E6 /( Zuzax::GasConstant * T);
    double hnowRT = hnow * 1.0E6 / ( Zuzax::GasConstant * T);

    if (fabs(hnowRT - HinputRT ) > 1.0E-14) {
         m_coeffs[0] += (HinputRT -  hnowRT ) * T;
    }
    double hnew = h(T);
    Hf298_s = h(298.15);
    if (fabs( Hf298_s ) < 1.0E-14) {
       Hf298_s = 0.0;
    }
    if (fabs(hnew - Hinput) > 1.0E-9) {
         printf("error\n");
         exit(-1);
    }
    return (hnew - hnow);
}
//==================================================================================================================================
double Thermo_NASA::adjustS(double T, double Sinput)
{
    double snow = s(T);
    double SinputR = Sinput * 1.0E3 /( Zuzax::GasConstant );
    double snowR = snow * 1.0E3 / ( Zuzax::GasConstant );

    if (fabs(snowR - SinputR) > 1.0E-14) {
         m_coeffs[1] += (SinputR - snowR);
    }
    double snew = s(T);
    S298_s = s(298.15);
    if (fabs(snew - Sinput) > 1.0E-9) {
         printf("error\n");
         exit(-1);
    }
    return (snew - snow);
}
//==================================================================================================================================
double Thermo_NASA::adjustCp(double T, double Cpinput)
{
    double snow = s(T);
    double hnow = h(T);
    double cpnow = cp(T);
    double delCoeff2 = (Cpinput - cpnow) / (Zuzax::GasConstant * 1.0E-3);
    m_coeffs[2] += delCoeff2;
    double cpnew = cp(T);
    (void) adjustH(T, hnow);
    (void) adjustS(T, snow);
    S298_s = s(298.15);
    Hf298_s = h(298.15);
    if (fabs(cpnew - Cpinput) > 1.0E-10) {
        printf("Thermo_NASA::adjustCp() ERROR: something didn't work\n");
        exit(-1);
    }
    return (cpnew - cpnow);
}
//==================================================================================================================================
void Thermo_NASA::printNASABlock(int n, int p)
{
     if (p < 5) {
         printf("p must be greater than 5\n");
         p = 5;
     } else if (p > 16) {
         printf("p must be loess than 17\n");
         p = 16;
     }
     char m_fmt[32];
     ind(n); printf("<NASA Tmax=\"%g\" Tmin=\"%g\" P0=\"100000\">\n", thigh_, tlow_);
     ind(n); printf("    <floatArray name=\"coeffs\" size=\"7\">\n");
     int wMin = p + 7;
     snprintf(m_fmt, 31, "%s -%d.%dE, ", "%", wMin, p);
     int m = n + 6;
     ind(m); printf(m_fmt, m_coeffs[2]);
     printf(m_fmt, m_coeffs[3]);
     printf(m_fmt, m_coeffs[4]);
     printf(m_fmt, m_coeffs[5]);
     printf("\n");
     ind(m); 
     printf(m_fmt, m_coeffs[6]);
     printf(m_fmt, m_coeffs[0]);
     snprintf(m_fmt, 31, "%s -%d.%dE ", "%", wMin, p);
     printf(m_fmt, m_coeffs[1]);
     ind(n); printf("</floatArray>\n");
     ind(n); printf("</NASA>\n");
}
//==================================================================================================================================
void Thermo_NASA::printThermoBlock(int n, int p)
{
     double val = Hf298_s;
     if (fabs(val) < 1.0E-5) {
         val = 0.0;
     }
     ind(n); printf("<!-- thermo , Hf298 = %g kJ/gmol S298= %g J/gmol/K -->\n", val, S298_s);
     ind(n); printf("<thermo>\n");
     printNASABlock(n+2, p); 
     ind(n); printf("</thermo>\n");
}
//==================================================================================================================================
//==================================================================================================================================
//==================================================================================================================================
double nasaMatch::adjustH_high()
{
    double hlow = nlow->h(tempMatch);
    double hhigh = nhigh->h(tempMatch);
    (void) nhigh->adjustH(tempMatch, hlow);
    return (hlow - hhigh);
}
//==================================================================================================================================
double nasaMatch::adjustS_high()
{
    double slow = nlow->s(tempMatch);
    double shigh = nhigh->s(tempMatch);
    (void) nhigh->adjustS(tempMatch, slow);
    return (slow - shigh);
}
//==================================================================================================================================
double nasaMatch::adjustCp_high()
{
    double Cplow = nlow->cp(tempMatch);
    double Cphigh = nhigh->cp(tempMatch);
    // This call will also ensure that H and S at tempMatch doesn't change
    (void) nhigh->adjustCp(tempMatch, Cplow);
    return (Cplow - Cphigh);
}
//==================================================================================================================================
double nasaMatch::adjustH_low()
{
    double hlow = nlow->h(tempMatch);
    double hhigh = nhigh->h(tempMatch);
    (void) nlow->adjustH(tempMatch, hhigh);
    return (hhigh - hlow);
}
//==================================================================================================================================
double nasaMatch::adjustS_low()
{
    double slow = nlow->s(tempMatch);
    double shigh = nhigh->s(tempMatch);
    (void) nlow->adjustS(tempMatch, shigh);
    return (shigh - slow);
}
//==================================================================================================================================
double nasaMatch::adjustCp_low()
{
    double Cplow = nlow->cp(tempMatch);
    double Cphigh = nhigh->cp(tempMatch);
    // This call will also ensure that H and S at tempMatch doesn't change
    (void) nlow->adjustCp(tempMatch, Cphigh);
    return (Cphigh - Cplow);
}
//==================================================================================================================================
double nasaMatch::adjustH_avg()
{
    double hlow = nlow->h(tempMatch);
    double hhigh = nhigh->h(tempMatch);
    double havg = 0.5 * (hlow - hhigh);
    (void) nlow->adjustH(tempMatch, havg);
    (void) nhigh->adjustH(tempMatch, havg);
    return (havg - hhigh);
}
//==================================================================================================================================
double nasaMatch::adjustS_avg()
{
    double slow = nlow->s(tempMatch);
    double shigh = nhigh->s(tempMatch);
    double savg = 0.5 * (slow + shigh);
    (void) nlow->adjustS(tempMatch, savg);
    (void) nhigh->adjustS(tempMatch, savg);
    return (savg - shigh);
}
//==================================================================================================================================
double nasaMatch::adjustCp_avg()
{
    double Cplow = nlow->cp(tempMatch);
    double Cphigh = nhigh->cp(tempMatch);
    double Cpnew = (Cplow + Cphigh) * 0.5;
    // This call will also ensure that H and S at tempMatch doesn't change
    (void) nlow->adjustCp(tempMatch, Cpnew);
    (void) nhigh->adjustCp(tempMatch, Cpnew);
    return (Cpnew - Cphigh);
}
//==================================================================================================================================
void Regions_NASA::addRegionPoly(Thermo_NASA* tNasa)
{
    if (NASA_list.size() == 0) {
        NASA_list.push_back(tNasa);
        tlowReg_ = tNasa->tlow_;
        thighReg_ = tNasa->thigh_;
    } else {
        Thermo_NASA* tNASA_ct = NASA_list.back();
        if (fabs(tNASA_ct->thigh_ - tNasa->tlow_) > 1.0E-8) {
          printf("Regions_NASA::addRegionPolyERROR: temperature regions aren't compatible");
          exit(-1);
        }
        NASA_list.push_back(tNasa);
        thighReg_ = tNasa->thigh_;
        nasaMatch newM(tNASA_ct->thigh_, tNASA_ct, tNasa);
        match_list.push_back(newM);
    }
}
//==================================================================================================================================
double Regions_NASA::adjust_high()
{
    double vH, vS, vCp;
    double deltaHmax = 0.0;
    int is = 0;
    for (nasaMatch& mm : match_list) {
        vH = mm.adjustH_high();
        vS = mm.adjustS_high();
        vCp = mm.adjustCp_high();
        double tt = mm.tempMatch;
        printf("Match %d at T = %g: Adjustment in regions: delta H = %g, deltaS = %g, deltaCp = %g \n", is, tt, vH, vS, vCp);
        is++;
        deltaHmax = std::max(vH, deltaHmax);
    }
    return deltaHmax;
}
//==================================================================================================================================
void Regions_NASA::printThermoBlock(int n, int p)
{

    Thermo_NASA* tNASA_ct = NASA_list[0];
    double hf298 = tNASA_ct->Hf298_s;
    double S298 = tNASA_ct->S298_s;
    double val = hf298; 
    if (fabs(val) < 1.0E-5) {
        val = 0.0;
    }
    ind(n); printf("<!-- thermo , Hf298 = %g kJ/gmol S298= %g J/gmol/K -->\n", val, S298);
    ind(n); printf("<thermo>\n");

    for (Thermo_NASA* tn : NASA_list) {
        tn->printNASABlock(n + 2, p);
    }
    ind(n); printf("</thermo>\n");
}
//=================================================================================================================================

