/*
 * $Id: Electrode_SimplePhaseChangeDiffusion.cpp 571 2013-03-26 16:44:21Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"

#include "Electrode_SimplePhaseChangeDiffusion.h"

using namespace std;
using namespace BEInput;
using namespace TKInput;

#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif

namespace Zuzax
{
//======================================================================================================================
Electrode_SimplePhaseChangeDiffusion::Electrode_SimplePhaseChangeDiffusion() :
    Electrode(),
    phaseIndexOuterSolidPhase_(-1),
    MD_OuterSolidPhase_(0.0),
    index_Li_int_(npos),
    pindex_Li7_int(npos),
    phaseIndexInnerSolidPhase_(-1),
    MD_InnerSolidPhase_(0.0),

    zeroedInnerRadius_(false) ,
    deltaTZeroed_(0.0),
    NoInnerSolid_(false),
    SolidInnerKSpecies_(-1),
    SolidInnerKSpeciesReacStoichCoeff_(0.0),
    SolidOuterKSpecies_(-1),
    SolidOuterKSpeciesProdStoichCoeff_(0.0),
    CAP_init_(0.0),
    CAP_final_(0.0)
{
}
//======================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
Electrode_SimplePhaseChangeDiffusion::Electrode_SimplePhaseChangeDiffusion(const Electrode_SimplePhaseChangeDiffusion& right) :
    Electrode_SimplePhaseChangeDiffusion()
{
    /*
     * Call the assignment operator.
     */
    *this = operator=(right);
}
//======================================================================================================================
// Assignment operator
/*
 *  @param right object to be copied
 */
Electrode_SimplePhaseChangeDiffusion&
Electrode_SimplePhaseChangeDiffusion::operator=(const Electrode_SimplePhaseChangeDiffusion& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    Electrode::operator=(right);


    zeroedInnerRadius_ = right.zeroedInnerRadius_;
    deltaTZeroed_ = right.deltaTZeroed_;
    /*
     * Return the reference to the current object
     */
    return *this;
}
//======================================================================================================================
Electrode_SimplePhaseChangeDiffusion::~Electrode_SimplePhaseChangeDiffusion()
{
}
//======================================================================================================================
Electrode* Electrode_SimplePhaseChangeDiffusion::duplMyselfAsElectrode() const
{
    Electrode_SimplePhaseChangeDiffusion* dd = new Electrode_SimplePhaseChangeDiffusion(*this);
    return dd;
}
//======================================================================================================================
Electrode_Types_Enum  Electrode_SimplePhaseChangeDiffusion::electrodeType() const
{
    return SIMPLE_PHASE_CHANGE_DIFFUSION_ET;
}
//====================================================================================================================
int
Electrode_SimplePhaseChangeDiffusion::electrode_model_create(ELECTRODE_KEY_INPUT* ei)
{

    Electrode::electrode_model_create(ei);


    surfIndexOuterSurface_ = 0;
    surfIndexInnerSurface_ = 1;

    index_Li_int_ = globalSpeciesIndex("Li(i)");

    pindex_Li7_int = globalPhaseIndex("Li7Si3_Interstitial");
    ThermoPhase* tp = & thermo(pindex_Li7_int);

    MD_int_ = tp->molarDensity();

    double Xmf[5];
    tp->getMoleFractions(Xmf);
    mf_internal_final_ = Xmf[0];
    mf_internal_init_  = Xmf[0];
    mf_external_final_ = Xmf[0];
    mf_external_init_  = Xmf[0];

    /*
     *  Identify the global index of the inner solid phase
     */
    size_t pindex_Li13_int = globalPhaseIndex("Li13Si4(S)");
    phaseIndexInnerSolidPhase_ = pindex_Li13_int;

    tp = & thermo(pindex_Li13_int);
    MD_InnerSolidPhase_ = tp->molarDensity();


    int pindex_Li7 = globalPhaseIndex("Li7Si3(S)");
    phaseIndexOuterSolidPhase_ = pindex_Li7;
    tp = & thermo(pindex_Li7);
    MD_OuterSolidPhase_ = tp->molarDensity();

    printf(" MD_OuterSolidPhase_ = %g\n", MD_OuterSolidPhase_);
    printf(" MD_int_    = %g\n", MD_int_);



    Radius_internal_final_ = inputParticleDiameter_/2.0;
    Radius_internal_init_  = inputParticleDiameter_/2.0;

    Radius_internal_init_ = Radius_exterior_init_ * (1.0 - 1.0E-8);
    Radius_internal_final_ =  Radius_internal_init_;


    DiffCoeff_ = 1.0E-9;


    double vol_inner = 4. * Pi / 3. *  Radius_internal_init_ *  Radius_internal_init_ *  Radius_internal_init_ *  particleNumberToFollow_;
    double num_inner =  MD_InnerSolidPhase_ * vol_inner;



    SolidInnerKSpecies_ = -1;
    //   ReactingSurDomain *rsd = RSD_List_[surfIndexInnerSurface_];
    for (size_t k = 0; k < m_NumTotSpecies; k++) {
        double  rgamma = reactantStoichCoeff(surfIndexInnerSurface_, k, 0);
        if (rgamma != 0.0) {
            if (SolidInnerKSpecies_!= -1) {
                throw ZuzaxError(" Electrode_SimplePhaseChangeDiffusion::electrode_model_create", "Undefined situation");
            }
            SolidInnerKSpecies_= k;
            SolidInnerKSpeciesReacStoichCoeff_ = rgamma;
        }
    }
    int kk =  m_PhaseSpeciesStartIndex[phaseIndexInnerSolidPhase_];
    if (SolidInnerKSpecies_ != kk) {
        throw ZuzaxError("Electrode_SimplePhaseChangeDiffusion::electrode_model_create", "confusion in SolidInnerKSpecies_");
    }


    for (size_t k = 0; k < m_NumTotSpecies; k++) {
        double  rgamma = productStoichCoeff(surfIndexInnerSurface_, k, 0);
        if (rgamma != 0.0) {
            int iph = globalPhaseIndexFromGlobalSpeciesIndex(k);
            ThermoPhase& tp = thermo(iph);
            if (tp.nSpecies() == 1) {
                SolidOuterKSpecies_= k;
                SolidOuterKSpeciesProdStoichCoeff_ = rgamma;
            }
            if (tp.nSpecies() == 1) {
                // SolidOuterKSpecies_= k;
                // SolidOuterKSpeciesReacStoichCoeff_ = rgamma;
            }
        }
    }

    spMoles_final_[SolidInnerKSpecies_] = num_inner;
    spMoles_init_[SolidInnerKSpecies_] = num_inner;

    double vol_outer =  4. * Pi / 3. * particleNumberToFollow_ *
                        (Radius_exterior_final_ * Radius_exterior_final_ * Radius_exterior_final_ -
                         Radius_internal_init_ * Radius_internal_init_ * Radius_internal_init_);
    double num_outer = vol_outer * MD_OuterSolidPhase_;

    spMoles_final_[SolidOuterKSpecies_] = num_outer;
    spMoles_init_[SolidOuterKSpecies_] = num_outer;

    spMoles_final_[SolidOuterKSpecies_+1] = num_outer * 10. * 0.01;
    spMoles_init_[SolidOuterKSpecies_+1] =   spMoles_final_[SolidOuterKSpecies_+1];

    spMoles_final_[SolidOuterKSpecies_+2] = num_outer * 10. * 0.99;
    spMoles_init_[SolidOuterKSpecies_ +2] =  spMoles_final_[SolidOuterKSpecies_+2];

    C_external_init_ = spMoles_init_[SolidOuterKSpecies_+1] / vol_outer;
    C_external_init_init_ =  C_external_init_;
    C_internal_init_ =  C_external_init_;
    C_internal_init_init_ = C_external_init_;
    C_external_final_ =  C_external_init_;
    C_internal_final_ =  C_external_init_;

    Radius_internal_init_init_ = Radius_internal_init_;
    Radius_internal_final_final_ = Radius_internal_init_;
    updateState();

    std::copy(spMoles_final_.begin(), spMoles_final_.end(), spMoles_init_.begin());
    std::copy(spMoles_final_.begin(), spMoles_final_.end(), spMoles_init_init_.begin());

    std::copy(phaseMoles_final_.begin(), phaseMoles_final_.end(), phaseMoles_init_.begin());
    std::copy(phaseMoles_final_.begin(), phaseMoles_final_.end(), phaseMoles_init_init_.begin());

    spMoleIntegratedSourceTermLast_.resize(m_NumTotSpecies, 0.0);

    return 0;
}
//====================================================================================================================
int
Electrode_SimplePhaseChangeDiffusion::setInitialConditions(ELECTRODE_KEY_INPUT* ei)
{

    Electrode::setInitialConditions(ei);



    index_Li_int_ = globalSpeciesIndex("Li(i)");

    pindex_Li7_int = globalPhaseIndex("Li7Si3_Interstitial");
    ThermoPhase* tp = & thermo(pindex_Li7_int);

    MD_int_ = tp->molarDensity();

    double Xmf[5];
    tp->getMoleFractions(Xmf);
    mf_internal_final_ = Xmf[0];
    mf_internal_init_  = Xmf[0];
    mf_external_final_ = Xmf[0];
    mf_external_init_  = Xmf[0];

    /*
     *  Identify the global index of the inner solid phase
     */
    int pindex_Li13_int = globalPhaseIndex("Li13Si4(S)");
    phaseIndexInnerSolidPhase_ = pindex_Li13_int;

    tp = & thermo(pindex_Li13_int);
    MD_InnerSolidPhase_ = tp->molarDensity();


    int pindex_Li7 = globalPhaseIndex("Li7Si3(S)");
    phaseIndexOuterSolidPhase_ = pindex_Li7;
    tp = & thermo(pindex_Li7);
    MD_OuterSolidPhase_ = tp->molarDensity();


    Radius_internal_final_ = inputParticleDiameter_/2.0;
    Radius_internal_init_  = inputParticleDiameter_/2.0;

    Radius_internal_init_ = Radius_exterior_init_ * (1.0 - 1.0E-8);
    Radius_internal_final_ =  Radius_internal_init_;


    DiffCoeff_ = 1.0E-9;


    double vol_inner = 4. * Pi / 3. *  Radius_internal_init_ *  Radius_internal_init_ *  Radius_internal_init_ *  particleNumberToFollow_;
    double num_inner =  MD_InnerSolidPhase_ * vol_inner;


    spMoles_final_[SolidInnerKSpecies_] = num_inner;
    spMoles_init_[SolidInnerKSpecies_] = num_inner;

    double vol_outer =  4. * Pi / 3. * particleNumberToFollow_ *
                        (Radius_exterior_final_ * Radius_exterior_final_ * Radius_exterior_final_ -
                         Radius_internal_init_ * Radius_internal_init_ * Radius_internal_init_);
    double num_outer = vol_outer * MD_OuterSolidPhase_;

    spMoles_final_[SolidOuterKSpecies_] = num_outer;
    spMoles_init_[SolidOuterKSpecies_] = num_outer;

    spMoles_final_[SolidOuterKSpecies_+1] = num_outer * 10. * 0.01;
    spMoles_init_[SolidOuterKSpecies_+1] =   spMoles_final_[SolidOuterKSpecies_+1];

    spMoles_final_[SolidOuterKSpecies_+2] = num_outer * 10. * 0.99;
    spMoles_init_[SolidOuterKSpecies_ +2] =  spMoles_final_[SolidOuterKSpecies_+2];

    C_external_init_ = spMoles_init_[SolidOuterKSpecies_+1] / vol_outer;
    C_external_init_init_ =  C_external_init_;
    C_internal_init_ =  C_external_init_;
    C_internal_init_init_ = C_external_init_;
    C_external_final_ =  C_external_init_;
    C_internal_final_ =  C_external_init_;

    Radius_internal_init_init_ = Radius_internal_init_;
    Radius_internal_final_final_ = Radius_internal_init_;
    updateState();

    std::copy(spMoles_final_.begin(), spMoles_final_.end(), spMoles_init_.begin());
    std::copy(spMoles_final_.begin(), spMoles_final_.end(), spMoles_init_init_.begin());

    std::copy(phaseMoles_final_.begin(), phaseMoles_final_.end(), phaseMoles_init_.begin());
    std::copy(phaseMoles_final_.begin(), phaseMoles_final_.end(), phaseMoles_init_init_.begin());



    return 0;
}
//====================================================================================================================
void  Electrode_SimplePhaseChangeDiffusion::check_initial_CAP()
{

    double ro2 =  Radius_exterior_init_ * Radius_exterior_init_;
    double ro3 = ro2 * Radius_exterior_init_;
    double ri2 = Radius_internal_init_ * Radius_internal_init_;
    double ri3 = ri2 * Radius_internal_init_;

    double cap = (Radius_exterior_init_ * C_external_init_ - Radius_internal_init_*C_internal_init_)/(Radius_exterior_init_ - Radius_internal_init_) * (4./3.*Pi*(ro3-ri3))
                 - 2.*Pi* Radius_exterior_init_*(C_external_init_ - C_internal_init_)/(Radius_exterior_init_ - Radius_internal_init_)*Radius_internal_init_*(ro2-ri2);

    cap *= particleNumberToFollow_;

    int kStart = m_PhaseSpeciesStartIndex[pindex_Li7_int];
    double spi = spMoles_init_[kStart];

    double denom = cap + spi + 1.0E4 * molarAtol_;

    //printf("cap_init = %11.5e , spi_init = %11.5e \n", cap, spi);
    if (fabs(cap-spi) / denom > 1.0E-4) {
        throw ZuzaxError("Electrode_SimplePhaseChangeDiffusion::check_initial_CAP() ",
                           "failed caps = " + fp2str(cap) + ", spi = " + fp2str(spi));
    }
}
//====================================================================================================================
void  Electrode_SimplePhaseChangeDiffusion::check_final_CAP()
{

    double ro2 =  Radius_exterior_final_ * Radius_exterior_final_;
    double ro3 = ro2 * Radius_exterior_final_;
    double ri2 = Radius_internal_final_ * Radius_internal_final_;
    double ri3 = ri2 * Radius_internal_final_;

    double cap = (Radius_exterior_final_*C_external_final_ - Radius_internal_final_*C_internal_final_)/(Radius_exterior_final_ - Radius_internal_final_) * (4./3.*Pi*(ro3-ri3))
                 - 2.*Pi*Radius_exterior_final_*(C_external_final_ - C_internal_final_)/(Radius_exterior_final_ - Radius_internal_final_)*Radius_internal_final_*(ro2-ri2);

    cap *= particleNumberToFollow_;

    size_t kStart = m_PhaseSpeciesStartIndex[pindex_Li7_int];
    double spi = spMoles_final_[kStart];

    double denom = cap + spi + 1.0E4 * molarAtol_;

    // printf("cap_final = %11.5e , spi_final = %11.5e \n", cap, spi);
    if (fabs(cap-spi) / denom > 1.0E-4) {
        if (spi == 0.0) {
        } else {
            throw ZuzaxError("Electrode_SimplePhaseChangeDiffusion::check_final_CAP() ",
                               "failed caps = " + fp2str(cap) + ", spi = " + fp2str(spi));
        }
    }
}
//====================================================================================================================
void Electrode_SimplePhaseChangeDiffusion::check_final_OuterVol()
{

    double ro2 =  Radius_exterior_final_ * Radius_exterior_final_;
    double ro3 = ro2 * Radius_exterior_final_;
    double ri2 = Radius_internal_final_ * Radius_internal_final_;
    double ri3 = ri2 * Radius_internal_final_;

    double vol = 4. * Pi / 3. * (ro3 - ri3);
    double cap = MD_OuterSolidPhase_ * vol * particleNumberToFollow_;

    int kStart = m_PhaseSpeciesStartIndex[phaseIndexOuterSolidPhase_];
    double spf = spMoles_final_[kStart];

    //    printf("SolidS2_final = %11.5e , spf_SolidS2_final = %11.5e \n", cap, spf);
    double denom = cap + spf + 1.0E4 * molarAtol_;
    if (fabs(cap-spf) / denom > 1.0E-4) {
        if (spf < 1.0E-20) {

        } else {
            throw ZuzaxError("Electrode_SimplePhaseChangeDiffusion::check_final_OuterVol() ",
                               "failed caps = " + fp2str(cap) + ", spi = " + fp2str(spf));
        }
    }
}
//==================================================================================================================================
bool Electrode_SimplePhaseChangeDiffusion::resetStartingCondition(double Tinitial, bool doTestsAlways)
{
    bool resetToInitInit = Electrode::resetStartingCondition(Tinitial);

    if (!resetToInitInit) {
    mf_internal_init_ =  mf_internal_final_;
    mf_external_init_ =  mf_external_final_;

    Radius_internal_init_init_ = Radius_internal_final_final_;
    Radius_exterior_init_init_ = Radius_exterior_final_final_;

    MN_internal_init_ = MN_internal_final_;

    zeroedInnerRadius_ = false;
    deltaTZeroed_ = 0.0;

    C_external_init_      = C_external_final_;
    C_external_init_init_ = C_external_init_;
    C_internal_init_      = C_internal_final_;
    C_internal_init_init_ = C_internal_final_;

    CAP_init_ = CAP_final_;
    }
    return resetToInitInit;
}
//==================================================================================================================================
void  Electrode_SimplePhaseChangeDiffusion::calcRate(double deltaTsubcycle)
{
    double  C_tmp, dadd_f, dadd_i;

    double ro2 =  Radius_exterior_final_ * Radius_exterior_final_;
    double ro3 = ro2 * Radius_exterior_final_;
    double ri2 = Radius_internal_final_ * Radius_internal_final_;
    double ri3 = ri2 * Radius_internal_final_;
    CAP_final_ = (Radius_exterior_final_ * C_external_final_ - Radius_internal_final_*C_internal_final_)/(Radius_exterior_final_ - Radius_internal_final_) * (4./3.*Pi*(ro3-ri3))
                 - 2.*Pi*Radius_exterior_final_*(C_external_final_ - C_internal_final_)/(Radius_exterior_final_ - Radius_internal_final_)*Radius_internal_final_*(ro2-ri2);
    deltaCAPdt_ = (CAP_final_ - CAP_init_)/deltaTZeroed_;
    double deltaCAPFdt_ = (CAP_final_)/deltaTZeroed_;
    double deltaCAPIdt_ = (CAP_init_)/deltaTZeroed_;




    if (C_external_final_ > 1.0E-200) {
        dadd_f = deltaCAPFdt_ / (4.*Pi*ro2* C_external_final_);
        dadd_i = deltaCAPIdt_ / (4.*Pi*ro2);
    } else {
        dadd_f = 0.0;
        dadd_i = 0.0;
    }

    if (Radius_internal_init_ <= 0.0) {
        if (Radius_internal_final_ <= 0.0) {
            NoInnerSolid_ = true;
        }
    }
    double C_eq_internal = k_f_internal_ / k_r_internal_;
    double C_eq_external = k_r_external_ * a_lip_ / k_f_external_;

    Da_s_ = k_r_internal_ *  Radius_internal_eff_ * (Radius_exterior_final_ - Radius_internal_eff_) /
            (Radius_exterior_final_ * DiffCoeff_);
#ifdef DEBUG_HKM
    checkFinite(Radius_exterior_final_);
#endif
    double r_ratio = Radius_internal_eff_ / Radius_exterior_final_;

    //   double denomN = 1.0 + Da_s_ + k_r_internal_ / k_f_external_ * r_ratio * r_ratio;
    double denomN = (1.0 + dadd_f/k_f_external_)*(1.0 + Da_s_) + k_r_internal_ / k_f_external_ * r_ratio * r_ratio;

    if (NoInnerSolid_) {
        C_external_final_ = (k_r_external_ * a_lip_ + dadd_i)/ (k_f_external_ + dadd_f);
        C_internal_final_ = C_external_final_;
        mf_external_final_ = C_external_final_ / MD_int_;
        mf_internal_final_ = C_internal_final_ / MD_int_;
    } else {

        Nflux_final_ = k_r_internal_ * (C_eq_internal*(1.0 + dadd_f/k_f_external_) - dadd_i/k_f_external_ - C_eq_external) / denomN;

        Nrate_final_ =  4. * Pi * Radius_exterior_final_ * Radius_exterior_final_ * Nflux_final_;

        //  denom = k_r_internal_*r_ratio*r_ratio + k_f_external_ * (1.0 + Da_s_);
        double denom = k_r_internal_*r_ratio*r_ratio + (k_f_external_ + dadd_f) * (1.0 + Da_s_);
        if (denom == 0.0) {
            printf("denom === 0.0\n");
        }

        C_tmp = (r_ratio * r_ratio * k_r_internal_ * C_eq_internal + (k_r_external_ * a_lip_ + dadd_i) * (1.0 + Da_s_)) / denom;
        if (C_tmp <= 0.0) {
            printf("Warning: C_external_final_ predicted to go negative %g\n", C_tmp);
        }
        if (C_tmp > 2.0 *  C_external_final_) {
            C_external_final_ = 2.0 * C_external_final_;
        } else {
            if (C_tmp <  C_external_final_ * 0.5) {
                C_external_final_ = 0.5 *  C_external_final_;
            } else {
                C_external_final_ = C_tmp;
            }
        }
#ifdef DEBUG_HKM
        checkFinite(C_external_final_);
#endif
        mf_external_final_ = C_external_final_ / MD_int_;

        C_internal_final_ = (Da_s_ * C_eq_internal + C_external_final_) / (1.0 + Da_s_);
        mf_internal_final_ = C_internal_final_ / MD_int_;
    }
    ROP_outer_ = k_f_external_ * C_external_final_ - k_r_external_ * a_lip_;
    ROP_rate_outer_ =  4. * Pi * Radius_exterior_final_ * Radius_exterior_final_ * ROP_outer_;

    if (NoInnerSolid_) {
        ROP_inner_ = 0.0;
    } else {
        ROP_inner_ = k_f_internal_ - k_r_internal_ * C_internal_final_;
    }
    ROP_rate_inner_ =  4. * Pi * Radius_internal_eff_ * Radius_internal_eff_ * ROP_inner_;
    // double diff_flux_outer =  DiffCoeff_ * (Radius_internal_eff_ / Radius_external_final_) *
    //                        (C_internal_final_ - C_external_final_) / (Radius_external_final_ -  Radius_internal_final_);

    //double diff_rate_outer = 4. * Pi * Radius_external_final_ * Radius_external_final_ * diff_flux_outer;


    MN_internal_init_ = 4. / 3. * Pi *  Radius_internal_init_ * Radius_internal_init_ *  Radius_internal_init_
                        * MD_InnerSolidPhase_ * particleNumberToFollow_;


    double delRad = - deltaTsubcycle  *  SolidInnerKSpeciesReacStoichCoeff_ * ROP_inner_ / MD_InnerSolidPhase_;

    Radius_internal_final_ =   Radius_internal_init_ + delRad;
    if (Radius_internal_init_ > 0.0) {
        if (Radius_internal_final_ <= 0.0) {
            deltaTZeroed_ =  Radius_internal_init_ * MD_InnerSolidPhase_ / (SolidInnerKSpeciesReacStoichCoeff_ * ROP_inner_);
            deltaTsubcycle = deltaTZeroed_;
            Radius_internal_final_ = 0.0;
            zeroedInnerRadius_ = true;
        } else {
            if (Radius_internal_final_ > Radius_exterior_final_) {
                Radius_internal_final_ = Radius_exterior_final_*(1.0-1.0E-6);
                delRad  = Radius_internal_final_- Radius_internal_init_;
                printf("Warning: R_i > R_o\n");
            }
        }
    } else {
        Radius_internal_final_ = 0.0;
    }


    if (Radius_internal_init_ > 0.0) {
        /*
         * Calculate the number of final moles of inner solid
         */
        MN_internal_final_ = 4. / 3. * Pi *  Radius_internal_final_ * Radius_internal_final_ *  Radius_internal_final_
                             * MD_InnerSolidPhase_ * particleNumberToFollow_;


        double delta_n = MN_internal_final_ -  MN_internal_init_;

        Radius_internal_eff_ =  0.5 * (Radius_internal_final_  + Radius_internal_init_);
        if (fabs(ROP_inner_) > 1.0E-30) {
            double r2 = delta_n / (ROP_inner_* (-SolidInnerKSpeciesReacStoichCoeff_) * deltaTsubcycle *  particleNumberToFollow_ * 4 * Pi);
            if (r2 > 0.0) {
                Radius_internal_eff_ = sqrt(r2);
            }
        }
#ifdef DEBUG_HKM
        checkFinite(Radius_internal_eff_);
#endif
    } else {
        Radius_internal_eff_  = 0.0;
    }
#ifdef DEBUG_HKM
    checkFinite(Radius_internal_eff_);
#endif
    /*
     * Calculate the change in the external radius due to the change in volume of the reaction
     */
    double deltaS2Vol = deltaTsubcycle * SolidOuterKSpeciesProdStoichCoeff_ * ROP_inner_ / MD_OuterSolidPhase_ * 4. * Pi * Radius_internal_eff_ * Radius_internal_eff_;
    double deltaro3 =  3./(4.*Pi) * deltaS2Vol + (Radius_internal_final_ * Radius_internal_final_ *  Radius_internal_final_ -
                       Radius_internal_init_  * Radius_internal_init_  *  Radius_internal_init_);
    ro3 = Radius_exterior_init_ * Radius_exterior_init_ * Radius_exterior_init_ + deltaro3;
    Radius_exterior_final_ = pow(ro3, 0.333333333333333);

    ro2 =  Radius_exterior_final_ * Radius_exterior_final_;
    ro3 = ro2 * Radius_exterior_final_;
    ri2 = Radius_internal_final_ * Radius_internal_final_;
    ri3 = ri2 * Radius_internal_final_;
    CAP_final_ = (Radius_exterior_final_*C_external_final_ - Radius_internal_final_*C_internal_final_)/(Radius_exterior_final_ - Radius_internal_final_) * (4./3.*Pi*(ro3-ri3))
                 - 2.*Pi*Radius_exterior_final_*(C_external_final_ - C_internal_final_)/(Radius_exterior_final_ - Radius_internal_final_)*Radius_internal_final_*(ro2-ri2);
#ifdef DEBUG_HKM
    checkFinite(CAP_final_);
#endif
    deltaCAPdt_ = (CAP_final_ - CAP_init_)/deltaTZeroed_;
}
//================================================================================================
/*
 * There is a small dependence on mf_external and mf_internal exhibited by this function
 */
void  Electrode_SimplePhaseChangeDiffusion::extractInfo(std::vector<size_t>& justBornMultiSpecies)
{
    size_t isk;
    size_t bornMultiSpecies = npos;
    double fwdROP[5], revROP[5], netROP[5];

    updateState();
    /*
     * Loop over surface phases, filling in the phase existence fields within the
     * kinetics operator
     */
    for (isk = 0; isk < numSurfaces_; ++isk) {
        /*
         *  Loop over phases, figuring out which phases have zero moles.
         *  Volume phases exist if the initial or final mole numbers are greater than zero
         *  Surface phases exist if the initial or final surface areas are greater than zero.
         */
        if (ActiveKineticsSurf_[isk]) {
            ReactingSurDomain* rsd = RSD_List_[isk];
            size_t nph = rsd->nPhases();
            for (size_t jph = 0; jph < nph; jph++) {
                size_t iph = rsd->globalPhaseIndex_fromKP(jph);
                if (iph == (size_t) metalPhase_) {
                    continue;
                }
                double mm = phaseMoles_init_[iph];
                double mmf = phaseMoles_final_[iph];
                if ((size_t) iph >=  m_NumVolPhases) {
                    // we are in a surface phase
                    int isur = iph -  m_NumVolPhases;
                    double sa_init = surfaceAreaRS_init_[isur];
                    double sa_final = surfaceAreaRS_final_[isur];
                    if (sa_init > 0.0 || sa_final > 0.0) {
                        rsd->setPhaseExistence(jph, true);
                    } else {
                        rsd->setPhaseExistence(jph, false);
                    }
                } else {
                    if (mm <= 0.0 && mmf <= 0.0) {
                        rsd->setPhaseExistence(jph, false);
                    } else {
                        rsd->setPhaseExistence(jph, true);
                    }
                }
                if (iph == bornMultiSpecies) {
                    rsd->setPhaseExistence(jph, true);
                }
                for (size_t iiph = 0; iiph < justBornMultiSpecies.size(); iiph++) {
                    if (iph == (size_t) justBornMultiSpecies[iiph]) {
                        rsd->setPhaseExistence(jph, true);
                    }
                }
            }
        }
    }


    for (isk = 0; isk < numSurfaces_; isk++) {
        // Loop over phases, figuring out which phases have zero moles.

        if (ActiveKineticsSurf_[isk]) {

            /*
             *  For each Reacting surface
             *
             *  Get the species production rates for the reacting surface
             */
            //    m_rSurDomain->getNetProductionRates(&RSSpeciesProductionRates_[0]);
            const vector<double>& rsSpeciesProductionRates = RSD_List_[isk]->veckin_NetProductionRates();
            RSD_List_[isk]->getNetRatesOfProgress(netROP);

            double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(isk);
            /*
             *  loop over the phases in the reacting surface
             *  Get the net production vector
             */
            std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
            size_t nphRS = RSD_List_[isk]->nPhases();
            size_t jph, kph;
            size_t kIndexKin = 0;
            for (kph = 0; kph < nphRS; ++kph) {
                jph = RSD_List_[isk]->globalPhaseIndex_fromKP(kph);
                size_t istart = m_PhaseSpeciesStartIndex[jph];
                size_t nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
                for (size_t k = 0; k < nsp; k++) {
                    spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
                    if (rsSpeciesProductionRates[kIndexKin] > 0.0) {
                        if ((phaseMoles_init_[jph] <= 0.0) && (jph != metalPhase_)) {
                            bool notFound = true;
                            for (size_t iiph = 0; iiph < justBornMultiSpecies.size(); iiph++) {
                                if (jph == justBornMultiSpecies[iiph]) {
                                    notFound = false;
                                }
                            }
                            if (notFound) {
                                if (nsp > 1) {
                                    bornMultiSpecies = jph;
                                } else {
                                    justBornMultiSpecies.push_back(jph);
                                }
                            }
                        }
                    }
                    kIndexKin++;
                }
            }
        }
    }

    updateState();
    size_t kStart = m_PhaseSpeciesStartIndex[pindex_Li7_int];
    ThermoPhase* tphase = & thermo(pindex_Li7_int);

    if (mf_external_final_ <= 1.0E-64) {
        printf("warning mf_external_final = %g\n", mf_external_final_);
        mf_external_final_ = 1.0E-64;
    }
    spMf_final_[kStart] = mf_external_final_;
    spMf_final_[kStart+1] = 1.0 - mf_external_final_;
    tphase->setState_TPX(temperature_, pressure_, &spMf_final_[kStart]);


    ReactingSurDomain* rsd_outer = RSD_List_[surfIndexOuterSurface_];
    rsd_outer->getFwdRatesOfProgress(fwdROP);
    rsd_outer->getRevRatesOfProgress(revROP);
    rsd_outer->getNetRatesOfProgress(netROP);
    C_external_final_ = MD_int_ * mf_external_final_;

    k_f_external_ = fwdROP[0] / C_external_final_;
    k_r_external_ = revROP[0];
    a_lip_ = 1.0;

    if (mf_internal_final_ <= 1.0E-64) {
        printf("warning mf_internal_final = %g\n", mf_internal_final_);
        mf_internal_final_ = 1.0E-64;
    }

    spMf_final_[kStart] = mf_internal_final_;
    spMf_final_[kStart+1] = 1.0 - mf_internal_final_;
    tphase->setState_TPX(temperature_, pressure_, &spMf_final_[kStart]);

    ReactingSurDomain* rsd_inner = RSD_List_[surfIndexInnerSurface_];
    rsd_inner->getFwdRatesOfProgress(fwdROP);
    rsd_inner->getRevRatesOfProgress(revROP);
    rsd_inner->getNetRatesOfProgress(netROP);
    C_internal_final_ = MD_int_ * mf_internal_final_;
    k_f_internal_ = fwdROP[0];
    k_r_internal_ = revROP[0]/ C_internal_final_;

}
//===============================================================================================================
void Electrode_SimplePhaseChangeDiffusion::updateState()
{
    /*
     * This may be redundant. However, I want to make sure mole fractions are
     * consistent with final moles.
     */
    for (size_t i = 0; i < m_NumTotPhases; i++) {
        updateState_Phase(i);
    }


    /*
     * Loop over all phases in the object
     */
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        ThermoPhase* tphase = &(thermo(iph));
        std::string pName = tphase->id();
        int kStart = m_PhaseSpeciesStartIndex[iph];
        /*
         * A direct call to tphase to set the state is unnecessary,
         * because m_mp will do it anyway.
         */
        tphase->setState_TPX(temperature_, pressure_, &spMf_final_[kStart]);
        tphase->setElectricPotential(phaseVoltages_[iph]);

        /*
         * Ok, we have set the state. Now upload the Vol and mu's.
         */
        tphase->getPartialMolarVolumes(& (VolPM_[kStart]));
        tphase->getElectrochemPotentials(& (spElectroChemPot_[kStart]));

        if (iph < m_NumVolPhases) {
            phaseMolarVolumes_[iph] = tphase->molarVolume();
        }

    }

    ElectrodeSolidVolume_ = SolidVol();

    //double vol =  ElectrodeSolidVolume_ / particleNumberToFollow_;

    // Radius_exterior_final_ = pow (vol * 3.0 / (4.0 * Pi), 0.3333333333333333);
}

//====================================================================================================================
//  Calculate the change in the state of the system when integrating from Tinitial to Tfinal
/*
 *  All information is kept internal within this routine. This may be done continuously
 *  and the solution is not updated.
 *
 *  Note the tolerance parameters refere to the nonlinear solves within the calculation
 *  They do not refer to time step parameters.
 *
 *
 *  @param deltaT        DeltaT for the integration step.
 *  @param GlobalRtolSrcTerm    Relative tolerance allowed for the electron source term over the interval.
 *                       This is a unitless quantity
 *                       Defaults to 1.0E-3
 *  @param fieldInterpolationType Type of interpolation of field variables defaults to 0
 *  @param subIntegrationType     Type of subintegration. defaults to 0
 *
 *  @return Returns the number of subcycle steps it took to complete the full step.
 */
int Electrode_SimplePhaseChangeDiffusion::integrate(double deltaT, double  GlobalRtolSrcTerm,
        Electrode_Exterior_Field_Interpolation_Scheme_Enum fieldInterpolationType,
        Subgrid_Integration_RunType_Enum subIntegrationType)
{

    counterNumberIntegrations_++;
    /*
     *   Set the internal objects to the correct conditions
     *    -> This will be the final conditions.
     */

    vector<double> Xf_tmp(m_NumTotSpecies, 0.0);
    vector<double> spMf_tmp(m_NumTotSpecies, 0.0);
    vector<double> spMoles_tmp(m_NumTotSpecies, 0.0);

    std::fill(spMoleIntegratedSourceTerm_.begin(), spMoleIntegratedSourceTerm_.end(), 0.);
    std::fill(spMoleIntegratedSourceTermLast_.begin(), spMoleIntegratedSourceTermLast_.end(), 0.);
    vector<double> deltaMoles(m_NumTotSpecies, 0.0);

#ifdef OLD_FOLLOW
    followElectrolyteMoles_ = 1;
    updateState();
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];
    followElectrolyteMoles_ = 0;
#else
    updateState();
    electrolytePseudoMoles_ = phaseMoles_final_[solnPhase_];
#endif
    int iterSubCycle = 0;
    bool notDone = true;

    double deltaTsubcycle = deltaT;
    double deltaTsubcycleNext =  deltaT;
    std::vector<size_t> justBornMultiSpecies(0);
    tinit_ = t_init_init_;
    tfinal_ = tinit_;
    t_final_final_ = t_init_init_ + deltaT;

    /*
     * We do a complete subcycle system here, even though the default is to use the same step as
     *  the main calculation. There will be times when we break up the step, and there may be
     *  problems where subcycling is warranted.
     */
    spMoles_init_  = spMoles_init_init_;
    // we can do better with this next line
    spMoles_final_ = spMoles_init_;

    Radius_internal_init_ =  Radius_internal_init_init_;
    Radius_exterior_init_ =  Radius_exterior_init_init_;

    check_initial_CAP();

    do {
        iterSubCycle++;
        counterNumberSubIntegrations_++;
        if (iterSubCycle > 1) {
            printf("we are here\n");
        }
        deltaTsubcycle = deltaTsubcycleNext;
        tinit_ = tfinal_;
        tfinal_ += deltaTsubcycleNext;
        if (tfinal_ > t_final_final_) {
            tfinal_ = t_final_final_;
            deltaTsubcycle = tfinal_ - tinit_;
        }
        deltaTZeroed_ = deltaTsubcycle;

        /*
         *   Update initial values of the state vector with the final values from the last subcycle
         *   iteration, if we are beyond the first subscycle
         */
        if (iterSubCycle > 1) {
            copy(spMoles_final_.begin(), spMoles_final_.end(), spMoles_init_.begin());
            copy(phaseMoles_final_.begin(), phaseMoles_final_.end(), phaseMoles_init_.begin());
            copy(surfaceAreaRS_final_.begin(), surfaceAreaRS_final_.end(), surfaceAreaRS_init_.begin());
            mf_internal_init_ = mf_internal_final_;
            Radius_internal_init_ = Radius_internal_final_;
            MN_internal_init_ = MN_internal_final_;
            Radius_exterior_init_ = Radius_exterior_final_;
            CAP_init_ = CAP_final_;
        }


        justBornMultiSpecies.clear();
        size_t bornMultiSpecies = npos;

        if (phaseMoles_init_[phaseIndexOuterSolidPhase_] <= 0.0) {
            justBornMultiSpecies.push_back(phaseIndexOuterSolidPhase_);
            bornMultiSpecies = pindex_Li7_int;
        }

restartStep:

        for (size_t isk = 0; isk < numSurfaces_; isk++) {
            // Loop over phases, figuring out which phases have zero moles.

            if (ActiveKineticsSurf_[isk]) {
                ReactingSurDomain* rsd = RSD_List_[isk];
                size_t nph = rsd->nPhases();
                for (size_t jph = 0; jph < nph; ++jph) {
                    size_t iph = rsd->globalPhaseIndex_fromKP(jph);
                    if (iph == metalPhase_) {
                        continue;
                    }
                    double mm = phaseMoles_init_[iph];
                    if (iph >= m_NumVolPhases) {
                        // we are in a surface phase
                        size_t isur = iph -  m_NumVolPhases;
                        double sa_init = surfaceAreaRS_init_[isur];
                        double sa_final = surfaceAreaRS_final_[isur];
                        if (sa_init > 0.0 || sa_final > 0.0) {
                            rsd->setPhaseExistence(jph, true);
                        } else {
                            rsd->setPhaseExistence(jph, false);
                        }
                    } else {
                        if (mm <= 0.0) {
                            rsd->setPhaseExistence(jph, false);
                        } else {
                            rsd->setPhaseExistence(jph, true);
                        }
                    }
                }
            }

        }


        /*
         *  This routine basically translates between species lists for the reacting surface
         *  domain and the Electrode.
         *  Later, when we have more than one reacting surface domain in the electrode object,
         *  this will do a lot more
         */

        for (size_t isk = 0; isk < numSurfaces_; isk++) {
            // Loop over phases, figuring out which phases have zero moles.

            if (ActiveKineticsSurf_[isk]) {

                /*
                 *  For each Reacting surface
                 *
                 *
                 *  Get the species production rates for the reacting surface
                 */
                //    m_rSurDomain->getNetProductionRates(&RSSpeciesProductionRates_[0]);
                const vector<double>& rsSpeciesProductionRates = RSD_List_[isk]->veckin_NetProductionRates();

                double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(isk);
                /*
                 *  loop over the phases in the reacting surface
                 *  Get the net production vector
                 */
                std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
                int nphRS = RSD_List_[isk]->nPhases();
                int jph, kph;
                int kIndexKin = 0;
                for (kph = 0; kph < nphRS; ++kph) {
                    jph = RSD_List_[isk]->globalPhaseIndex_fromKP(kph);
                    int istart = m_PhaseSpeciesStartIndex[jph];
                    int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
                    for (int k = 0; k < nsp; ++k) {
                        spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
                        kIndexKin++;
                    }
                }

            }
        }

        /*
         *  Find the initial surface area to use
         */
        //double sa_init = surfaceAreaRS_init_[0];

        /*
         *  Find the final surface area to use
         */
        double sa_final = calcSurfaceAreaChange(deltaTsubcycle);
        surfaceAreaRS_final_[0] = sa_final;
        /*
         * For phases which are just born, we need to start with a seed in order for the
         * algorithm to work.
         *  The seed needs to be set at a fraction of the initial mole number
         */
        if (bornMultiSpecies != npos) {
            ThermoPhase* tp = PhaseList_[bornMultiSpecies];
            tp->getMoleFractions(DATA_PTR(Xf_tmp));
            int retn = phasePop(bornMultiSpecies, DATA_PTR(Xf_tmp), deltaTsubcycle);
            if (retn == 0) {
                tp->setMoleFractions(DATA_PTR(Xf_tmp));
                size_t istart = m_PhaseSpeciesStartIndex[bornMultiSpecies];
                size_t nsp =  m_PhaseSpeciesStartIndex[bornMultiSpecies+1] - istart;
                for (size_t kp = 0; kp < nsp; kp++) {
                    size_t k = istart + kp;
                    spMf_tmp[k] = Xf_tmp[kp];
                }
                if (bornMultiSpecies == pindex_Li7_int) {
                    mf_internal_final_ = Xf_tmp[0];
                    mf_internal_init_  = Xf_tmp[0];
                    mf_external_final_ = Xf_tmp[0];
                    mf_external_init_  = Xf_tmp[0];
                }
            }
            justBornMultiSpecies.push_back(bornMultiSpecies);
            bornMultiSpecies = npos;
            goto   restartStep;
        }

        /*
         *  Note: I have found that it may be necessary to put extractInfo and calcRate
         *        in a loop to get convergence. It's due to the fact that the vacancy
         *        species creates inadequacies in the analytical solution. Successive
         *        substitution seems to work here.
         */
        Radius_internal_eff_ = Radius_internal_final_;
        extractInfo(justBornMultiSpecies);
        calcRate(deltaTsubcycle);
        double c_new = C_internal_final_;
        double c_old;
        double delta_c, convV;
        for (int iter = 0; iter < 50; iter++) {
            c_old = c_new;

            extractInfo(justBornMultiSpecies);

            calcRate(deltaTsubcycle);

            c_new = C_internal_final_;
            delta_c = c_new - c_old;
            convV = fabs(delta_c) / (c_old + c_new);
            if (convV < 1.0E-9) {
                break;
            }

        }
        if (convV > 1.0E-9) {
            printf("warning, not conv %g\n", convV);
        }


        // Calculate the net ROPs now that we have the concentrations

        /*
         * Set up the outer surface calculation
         */
        int kStart = m_PhaseSpeciesStartIndex[pindex_Li7_int];
        ThermoPhase* tphase = & thermo(pindex_Li7_int);
        double netROP[10];

        spMf_final_[kStart] = mf_external_final_;
        spMf_final_[kStart+1] = 1.0 - mf_external_final_;
        tphase->setState_TPX(temperature_, pressure_, &spMf_final_[kStart]);
        const vector<double>& rsSpeciesProductionRatesO = RSD_List_[surfIndexOuterSurface_]->veckin_NetProductionRates();
        RSD_List_[surfIndexOuterSurface_]->getNetRatesOfProgress(netROP);
        double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(surfIndexOuterSurface_);
        std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
        int nphRS = RSD_List_[surfIndexOuterSurface_]->nPhases();
        int jph, kph;
        int kIndexKin = 0;
        for (kph = 0; kph < nphRS; kph++) {
            jph = RSD_List_[surfIndexOuterSurface_]->globalPhaseIndex_fromKP(kph);
            int istart = m_PhaseSpeciesStartIndex[jph];
            int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
            for (int k = 0; k < nsp; k++) {
                spNetProdPerArea[istart + k] += rsSpeciesProductionRatesO[kIndexKin];
                kIndexKin++;
            }
        }
        /*
         * Set up the inner surface calculation
         */

        kStart = m_PhaseSpeciesStartIndex[pindex_Li7_int];
        tphase = & thermo(pindex_Li7_int);

        spMf_final_[kStart] = mf_internal_final_;
        spMf_final_[kStart+1] = 1.0 - mf_internal_final_;
        tphase->setState_TPX(temperature_, pressure_, &spMf_final_[kStart]);
        const vector<double>& rsSpeciesProductionRatesI = RSD_List_[surfIndexInnerSurface_]->veckin_NetProductionRates();
        RSD_List_[surfIndexInnerSurface_]->getNetRatesOfProgress(netROP);
        spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(surfIndexInnerSurface_);
        std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
        nphRS = RSD_List_[surfIndexInnerSurface_]->nPhases();
        kIndexKin = 0;
        for (kph = 0; kph < nphRS; kph++) {
            jph = RSD_List_[surfIndexInnerSurface_]->globalPhaseIndex_fromKP(kph);
            int istart = m_PhaseSpeciesStartIndex[jph];
            int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
            for (int k = 0; k < nsp; k++) {
                spNetProdPerArea[istart + k] += rsSpeciesProductionRatesI[kIndexKin];
                kIndexKin++;
            }
        }



        // ok we have the rate of a particle in Nrate_final

        MN_internal_init_ = 4. / 3. * Pi *  Radius_internal_init_ * Radius_internal_init_ *  Radius_internal_init_
                            * MD_InnerSolidPhase_ * particleNumberToFollow_;


        double delRad = - deltaTsubcycle  * SolidInnerKSpeciesReacStoichCoeff_ * ROP_inner_ / MD_InnerSolidPhase_;

        Radius_internal_final_ =   Radius_internal_init_ + delRad;
        if (Radius_internal_final_ < 0.0) {
            if (zeroedInnerRadius_) {
                Radius_internal_final_ = 0.0;
                tfinal_ -= deltaTsubcycle;
                deltaTsubcycle = deltaTZeroed_;
                tfinal_ += deltaTsubcycle;
            }
        }

        MN_internal_final_ = 4. / 3. * Pi *  Radius_internal_final_ * Radius_internal_final_ *  Radius_internal_final_
                             * MD_InnerSolidPhase_ * particleNumberToFollow_;

        surfaceAreaRS_final_[surfIndexOuterSurface_] =
            4. * Pi * Radius_exterior_final_ *  Radius_exterior_final_ * particleNumberToFollow_;
        surfaceAreaRS_final_[surfIndexInnerSurface_] =
            4. * Pi * Radius_internal_final_ *  Radius_internal_final_ * particleNumberToFollow_;

        /*
         * This is the effective Li consumption during the interval
         */

        double surAreaInternalEff = 4. * Pi *  Radius_internal_eff_ *  Radius_internal_eff_ * particleNumberToFollow_;



        std::fill(deltaMoles.begin(), deltaMoles.end(), 0.);
        for (size_t k = 0; k < m_NumTotSpecies; k++) {
            spMoles_tmp[k] = spMoles_init_[k];
            for (size_t isk = 0; isk < numSurfaces_; isk++) {
                if (ActiveKineticsSurf_[isk]) {
                    double sa = surAreaInternalEff;
                    if (isk == 0) {
                        sa = surfaceAreaRS_final_[surfIndexOuterSurface_];
                    }
                    double* spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(isk);
                    spMoles_tmp[k]  += deltaTsubcycle * (sa) * spNetProdPerArea[k];

                    deltaMoles[k] +=  deltaTsubcycle * (sa) * spNetProdPerArea[k];
                }
            }
        }

        for (size_t k = 0; k < m_NumTotSpecies; k++) {
            size_t iph = globalPhaseIndexFromGlobalSpeciesIndex(k);
            if (iph == metalPhase_) {
                continue;
            }
            /*
             * Negative moles may occur in Li(i) due to roundoff in the calculation. -> it really is trying
             *  to get it to zero.
             */
            if (spMoles_tmp[k] < 0.0) {
                printf("Warning sPMoles_tmp[%d, iph=%d] = %g\n", static_cast<int>(k), static_cast<int>(iph), spMoles_tmp[k]);
                spMoles_tmp[k] = 0.0;
            }
        }

        /*
         *  Accept the tmp vector and assign it to the final_ vector
         */
        for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
            ThermoPhase& tp = thermo(iph);
            string pname = tp.id();
            size_t istart = m_PhaseSpeciesStartIndex[iph];
            size_t nsp = tp.nSpecies();
            for (size_t ik = 0; ik < nsp; ik++) {
                size_t k = istart + ik;
                spMoleIntegratedSourceTerm_[k] += (spMoles_tmp[k] - spMoles_init_[k]);
                spMoleIntegratedSourceTermLast_[k] = (spMoles_tmp[k] - spMoles_init_[k]);
                if (iph == solnPhase_) {
                    continue;
                }
                if (iph == metalPhase_) {
                    continue;
                }
                spMoles_final_[k] = spMoles_tmp[k];
            }
        }

        if (tfinal_ >= (t_final_final_)) {
            notDone = false;
        }

        updateState();

        /*
         *  Do a check on the final molar amounts of the inner solid
         */

        if (fabs(MN_internal_final_ -   spMoles_final_[SolidInnerKSpecies_]) > molarAtol_) {
            throw ZuzaxError("Electrode_SimplePhaseChangeDiffusion::integrate()", " mb errr");
        } else {
            /*
             *  If a phase should be zeroed, make sure that phase is zeroed.
             */
            if (zeroedInnerRadius_ && (MN_internal_final_ <= 1.0E-200)) {
                spMoles_final_[SolidInnerKSpecies_] = 0.0;
            }
        }

        /*
         *  Do a balance check on the diffusing species
         */
        kStart      = m_PhaseSpeciesStartIndex[pindex_Li7_int];
        double deltaCap = spMoleIntegratedSourceTerm_[kStart];
        double fluxR = ROP_inner_ * 4. * Pi * Radius_internal_eff_ *  Radius_internal_eff_ * particleNumberToFollow_ * deltaTsubcycle;
        double fluxL = ROP_outer_ * 4. * Pi * Radius_exterior_final_ * Radius_exterior_final_ * particleNumberToFollow_ * deltaTsubcycle;
        double denomE = fabs(deltaCap) + fabs(fluxR) + fabs(fluxL) + 1.0E5*molarAtol_;
        double balLi = deltaCap - (fluxR - fluxL);
        if (denomE > 1.0E-200) {
            if (balLi/denomE > 1.0E-3) {
                throw ZuzaxError("", "balLi error");
            }
        }
        check_final_CAP();
        check_final_OuterVol();
    } while (notDone);

    /*
     *  Copy the results into the final holding pens.
     */
    spMoles_final_final_ = spMoles_final_;
    Radius_internal_final_final_ = Radius_internal_final_;
    Radius_exterior_final_final_ = Radius_exterior_final_;
    /*
     * indicate that we have a pending integrated step
     */
    pendingIntegratedStep_ = 1;
    return  iterSubCycle;
}
//====================================================================================================================
void Electrode_SimplePhaseChangeDiffusion::printElectrode(int pSrc, bool subTimeStep)
{
    double* netROP = new double[m_NumTotSpecies];
    double egv = TotalVol();
    printf("   ===============================================================\n");
    if (subTimeStep) {
        printf("      Electrode at intermediate-step time final = %g\n", tfinal_);
        printf("                   intermediate-step time init  = %g\n", tinit_);
    } else {
        printf("      Electrode at time final = %g\n", t_final_final_);
        printf("                   time init  = %g\n", t_init_init_);
    }
    printf("   ===============================================================\n");
    printf("          Number of external surfaces = %d\n", (int) numExternalInterfacialSurfaces_);
    printf("          Solid Volume = %10.3E\n", ElectrodeSolidVolume_);
    printf("          Total Volume = %10.3E\n", egv);
    printf("          Temperature = %g\n", temperature_);
    printf("          Pressure = %g\n", pressure_);
    double capacd = capacityDischarged();
    printf("          Capacity Discharged = %g coulombs = %g Ah\n", capacd, capacd / 3600.);
    printf("\n");
    printf("          followElectrolyteMoles = %d\n", followElectrolyteMoles_);
    printf("          ElectrolytePseudoMoles = %g\n",  electrolytePseudoMoles_);

    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
        printElectrodePhase(iph, pSrc);
        printf("     ===============================================================\n");
    }
    delete [] netROP;
}
//===================================================================================================================
void Electrode_SimplePhaseChangeDiffusion::printElectrodePhase(size_t iph, int pSrc, bool subTimeStep)
{
    size_t isph = npos;
    double* netROP = new double[m_NumTotSpecies];
    ThermoPhase& tp = thermo(iph);
    std::string pname = tp.name();
    size_t istart = m_PhaseSpeciesStartIndex[iph];
    size_t nsp = tp.nSpecies();
    printf("     ===============================================================\n");
    printf("          Phase %d %s \n", static_cast<int>(iph), pname.c_str());
    printf("                Total moles = %g\n", phaseMoles_final_[iph]);
    if ((size_t) iph == metalPhase_) {
        double deltaT = t_final_final_ - t_init_init_;
        if (subTimeStep) {
            deltaT = tfinal_ - tinit_;
        }
        if (deltaT > 1.0E-200) {
            double amps = spMoleIntegratedSourceTerm_[istart] / deltaT * Faraday;
            if (subTimeStep) {
                amps = spMoleIntegratedSourceTermLast_[istart] / deltaT * Faraday;
            }
            printf("                Current = %g amps \n", amps);
        } else {
            printf("                Current = NA amps \n");
        }
    }
    if (iph == metalPhase_ || iph == solnPhase_) {
        printf("                  Voltage = %g\n", tp.electricPotential());
    }
    if (iph >= m_NumVolPhases) {
        isph = iph - m_NumVolPhases;
        printf("                surface area (final) = %11.5E m2\n",  surfaceAreaRS_final_[isph]);
        printf("                surface area (init)  = %11.5E m2\n",  surfaceAreaRS_init_[isph]);
        int ddd =  isExternalSurface_[isph];
        printf("                IsExternalSurface = %d\n", ddd);
        double oc = openCircuitVoltage(isph);
        if (oc != 0.0) {
            printf("                 Open Circuit Voltage = %g\n", oc);
        }
    }
    printf("\n");
    printf("                Name               MoleFrac_final  kMoles_final kMoles_init SrcTermLastStep(kMoles)\n");
    for (size_t k = 0; k < nsp; k++) {
        std::string sname = tp.speciesName(k);
        if (pSrc) {
            if (subTimeStep) {
                printf("                %-22s %10.3E %10.3E   %10.3E  %10.3E\n", sname.c_str(), spMf_final_[istart + k],
                       spMoles_final_[istart + k], spMoles_init_[istart + k],
                       spMoleIntegratedSourceTermLast_[istart + k]);
            } else {
                printf("                %-22s %10.3E %10.3E   %10.3E  %10.3E\n", sname.c_str(), spMf_final_[istart + k],
                       spMoles_final_[istart + k], spMoles_init_init_[istart + k],
                       spMoleIntegratedSourceTerm_[istart + k]);
            }
        } else {
            if (subTimeStep) {
                printf("                %-22s %10.3E %10.3E   %10.3E\n", sname.c_str(), spMf_final_[istart + k],
                       spMoles_final_[istart + k],   spMoles_init_[istart + k]);
            } else {
                printf("                %-22s %10.3E %10.3E   %10.3E\n", sname.c_str(), spMf_final_[istart + k],
                       spMoles_final_[istart + k],   spMoles_init_init_[istart + k]);
            }
        }
    }
    if (iph >= m_NumVolPhases) {
        const vector<double>& rsSpeciesProductionRates = RSD_List_[isph]->veckin_NetProductionRates();
        RSD_List_[isph]->getNetRatesOfProgress(netROP);

        double* spNetProdPerArea = (double*) spNetProdPerArea_List_.ptrColumn(isph);
        std::fill_n(spNetProdPerArea, m_NumTotSpecies, 0.);
        size_t nphRS = RSD_List_[isph]->nPhases();
        size_t kIndexKin = 0;
        for (size_t kph = 0; kph < nphRS; ++kph) {
            size_t jph = RSD_List_[isph]->globalPhaseIndex_fromKP(kph);
            size_t istart = m_PhaseSpeciesStartIndex[jph];
            size_t nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
            for (size_t k = 0; k < nsp; ++k) {
                spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
                kIndexKin++;
            }
        }
        printf("\n");
        printf("                           spName                  Source (kmol/m2/s) \n");
        for (size_t k = 0; k <  m_NumTotSpecies; k++) {
            string ss = speciesName(k);
            printf("                           %-22s %10.3E\n", ss.c_str(), spNetProdPerArea[k]);
        }
    }
    delete [] netROP;

}
//==================================================================================================================================
} // End of namespace
//----------------------------------------------------------------------------------------------------------------------------------
