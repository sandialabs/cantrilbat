/**
 *  @file RSD_OCVmodel.cpp

 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */
#include "RSD_OCVmodel.h"
#include "Electrode_Factory.h"

#include "cantera/thermo/ThermoPhase.h"


#include "mdp_stringUtils.h"
#include "BI_InputError.h"
#include "tok_input_util.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>

//-----------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//=============================================================================================================================
//  Here we make a map of all ModelID's to string associations that we support
/*
 *
 */
void createOCVmodel_map(std::map<int, std::string>& smap)
{
   smap.clear();
   smap[OCVAnode_CONSTANT]          = "ANODE_CONSTANT";
   smap[OCVCathode_CONSTANT]        = "CATHODE_CONSTANT";
   smap[OCVAnode_MCMB2528]          = "MCMB2528";
   smap[OCVAnode_MCMB2528_dualfoil] = "MCMB2528_dualfoil";
   smap[OCVCathode_CoO2_dualfoil]   = "CoO2_dualfoil";
}
//==============================================================================================================================
//==============================================================================================================================
//  Constructor 
RSD_OCVmodel::RSD_OCVmodel(int modelID) :
    modelID_(modelID),
    rsd_ptr_(0),
    solidPhaseModel_(0),
    kSpecies_DoD_(npos),
    OCV_Format_(0),
    relExtent_(-1.0),
    xMF_(0),
    temperatureDerivType_(0),
    temperatureBase_(298.15),
    dvec_(0)
{
    //
    //  Match the modelID with a string name of the model
    //  It is not an error to not make a match. Inheritance may be used to provide matches.
    //
    modelName_ = modelID_to_stringName_RCD_OCVmodel(modelID);
}
//===============================================================================================================================
RSD_OCVmodel::RSD_OCVmodel(const RSD_OCVmodel& right) :
    modelID_ (right.modelID_),
    modelName_(right.modelName_),
    rsd_ptr_(right.rsd_ptr_),
    solidPhaseModel_(right.solidPhaseModel_),
    kSpecies_DoD_(right.kSpecies_DoD_),
    OCV_Format_(right.OCV_Format_),
    relExtent_ (right.relExtent_),
    xMF_(right.xMF_),
    temperatureDerivType_(right.temperatureDerivType_),
    temperatureBase_(right.temperatureBase_),
    dvec_(right.dvec_)
{
}
//===============================================================================================================================
RSD_OCVmodel::~RSD_OCVmodel()
{
}
//===============================================================================================================================
RSD_OCVmodel& RSD_OCVmodel::operator=(const RSD_OCVmodel& right)
{
    if (this == &right) {
       return *this;
    }
    modelID_ = right.modelID_;
    modelName_ = right.modelName_;
    //
    //  BEWARE:  This is a shallow copy
    //
    rsd_ptr_ = right.rsd_ptr_;
    //
    //  BEWARE:  This is a shallow copy. It will have to be fixed up outside
    //
    solidPhaseModel_ = right.solidPhaseModel_;
    kSpecies_DoD_ = right.kSpecies_DoD_;
    OCV_Format_   = right.OCV_Format_;
    relExtent_  = right.relExtent_;
    xMF_ = right.xMF_;
    temperatureDerivType_ = right.temperatureDerivType_;
    temperatureBase_ = right.temperatureBase_;
    dvec_ = right.dvec_;

    return *this;
}
//===============================================================================================================================
RSD_OCVmodel* RSD_OCVmodel::duplMyselfAsOCVmodel(ThermoPhase *solidPhase) const
{
    RSD_OCVmodel* pp = new RSD_OCVmodel(*this);
    if (solidPhase) {
	pp->assignShallowPointers(solidPhase);
    }
    return pp;    
}
//===============================================================================================================================
void RSD_OCVmodel::initialize(ReactingSurDomain * rsd_ptr, const OCV_Override_input& OCVinput)
{
    rsd_ptr_ = rsd_ptr;
    /*
    int numTimes;
    int surfacePhaseID;
    std::string surfacePhaseName;
    std::string OCVModel;
    std::string replacedSpeciesName;
   
    //! the global species id for the species whose thermo will be replaced
    int replacedGlobalSpeciesID;
    int replacedLocalSpeciesID;
    int replacedSpeciesPhaseID;
    //! OCV_Format defaults to 0
    int OCV_Format_;
    std::string DoDSurrogateSpeciesName;
    size_t MF_DoD_LocalSpeciesID;
    int rxnID;
    int temperatureDerivType;
    double temperatureBase;
    std::string OCVTempDerivModel;
    */
    // surfasePhaseID -> not relevant at this point as we are on a particular surface at this point
    // surfacePhaseName -> not used
    // OCVModel -> modelName_
    // numTimes -> not relevant at this point
    // replacedSpeciesName -> already processed to produce species index.
    if (! mdpUtil::LowerCaseStringEquals(modelName_, OCVinput.OCVModel)) {
	throw Electrode_Error("RSD_OCVmodel::initialize", "Model names aren't consistent: " + modelName_ + 
			      " vs " + OCVinput.OCVModel);
    }
    modelName_ = OCVinput.OCVModel;
    dvec_.resize(2);
    // replacedGlobalSpeciesID -> don't want to add PhaseList 
    // replacedLocalSpeciesID; -> don't want to put replaced species stuff here
    OCV_Format_ = OCVinput.OCV_Format_;

    temperatureDerivType_ =  OCVinput.temperatureDerivType;
    
    OCVTempDerivModel_ =  OCVinput.OCVTempDerivModel;
    TKInput::TOKEN toktd(OCVinput.OCVTempDerivModel);
    OCVTempDerivModel_ = toktd.tok_ptrV[0];
    if (toktd.ntokes > 1) {
      BOOLEAN err = 0;
      double val = TKInput::str_to_double(toktd.tok_ptrV[1], 1.0E300, -1.0E300, 0.0, &err);
      if (err) {
            throw BEInput::BI_InputError("TKInput::str_to_double()",
                                         "str_to_double interpretation failed: " + std::string(toktd.tok_ptrV[1]));
      }
      dvec_.resize(10);
      dvec_[1] =   val;
    }
    temperatureBase_ = OCVinput.temperatureBase;
    temperatureDerivModelType_ = stringName_RCD_OCVmodel_to_modelID(OCVTempDerivModel_);
    if (temperatureDerivModelType_ == -1) {
        throw Electrode_Error("RSD_OCVmodel::initialize()", "Unknown temperature derivative model: " + OCVTempDerivModel_);
    }
}
//===============================================================================================================================
 void  RSD_OCVmodel::assignShallowPointers(ThermoPhase* solidPhase)
 {
     solidPhaseModel_ = solidPhase;
 }
//=============================================================================================================================== 
void RSD_OCVmodel::setup_RelExtent(ThermoPhase *tp, size_t kspec, double *dvec, int *ivec)
{
    //
    //  Simplest treatment assigns the relative extent of reaction to a single mole fraction value in the phase
    //  But, other models may be used.
    //
    solidPhaseModel_ = tp;
    kSpecies_DoD_ = kspec;
    size_t nsp = solidPhaseModel_->nSpecies();
    xMF_.resize(nsp, 0.0);
    solidPhaseModel_->getMoleFractions(& xMF_[0]);

    if (modelID_ == OCVAnode_CONSTANT || modelID_ == OCVCathode_CONSTANT ) {
	dvec_.resize(2);
	dvec_[0] = dvec[0];
    }
}
//=============================================================================================================================== 
const ThermoPhase* RSD_OCVmodel::solidPhasePtr() const
{
    return solidPhaseModel_;
}
//=============================================================================================================================== 
void RSD_OCVmodel::calcRelExtent() const
{
    AssertThrow(kSpecies_DoD_ != npos, "RSD_OCVmodel: kSpecies_DoD_ is not assigned");
    // We do slightly more than is necessary so that models may have the full mole fraction vector
    solidPhaseModel_->getMoleFractions(&xMF_[0]);
    relExtent_ = xMF_[kSpecies_DoD_];
}
//=============================================================================================================================== 
double RSD_OCVmodel::OCV_value(double temp) const
{
    calcRelExtent();
    double volts ;
    double xLi = 1.0 - relExtent_;
    double xV;
    double dOCVdt = 0.0;
    if (temp != temperatureBase_) {
	dOCVdt = OCV_dvaldT(temp);
    }
    if (modelID_ == OCVAnode_CONSTANT) {
       volts = dvec_[0];

    } else if (modelID_ == OCVAnode_MCMB2528) {
                
       //  MCMB 2528 graphite measured by Chris Bogatu 2000, 
       //  Telcordia and PolyStor materials.
    
       volts = (  0.124   + 1.5 * exp(-150.0 * xLi) 
                   + 0.0351 * tanh( (xLi - 0.286) / 0.083)
                   - 0.0045 * tanh( (xLi - 0.90)  / 0.119)
                   - 0.035  * tanh( (xLi - 0.99)  / 0.05 )
                   - 0.0147 * tanh( (xLi - 0.50)  / 0.034)
                   - 0.102  * tanh( (xLi - 0.194) / 0.142)
                   - 0.022  * tanh( (xLi - 0.98 ) / 0.0164)
                   - 0.011  * tanh( (xLi - 0.124) / 0.0226)
                   + 0.0155 * tanh( (xLi - 0.105) / 0.029));

    } else if (modelID_ == OCVAnode_MCMB2528_dualfoil) { 

	/*
	 *   Line 4174 of dualfoil 5.1
	 
	 MCMB 2510 carbon (Bellcore)
	 c      c1=-0.160d0
	 c      c2=1.32d0
	 c      c3=-3.0d0
	 c      g0=c1+c2*expg(c3*xx(2+mpa,j)/ct1)
	 c      g0=g0+10.d0*expg(-2000.d0*xx(2+mpa,j)/ct1)
	 c      g1=c2*c3*expg(c3*xx(2+mpa,j)/ct1)/ct1
	 c      g1=g1-10.d0*2000.d0/ct1*expg(-2000.d0*xx(2+mpa,j)/ct1)
	 c     MCMB 2528 graphite measured by Chris Bogatu 2000, Telcordia and PolyStor materials
	 c     for 0.01 < x < 0.9
	 *
	 *   (note this agrees exactly with dualfoil)
	 */

	volts = ( 0.194 + 1.5 * exp(-120.0 * xLi)
		  + 0.0351 * tanh( (xLi - 0.286)   / 0.083)
		  - 0.0045 * tanh( (xLi - 0.849)   /0.119)
		  - 0.035 *  tanh( (xLi - 0.9233)  /0.05)
		  - 0.0147 * tanh( (xLi - 0.5)     /0.034)
		  - 0.102 *  tanh( (xLi - 0.194)   /0.142)
		  - 0.022 *  tanh( (xLi - 0.9)     /0.0164)
		  - 0.011 *  tanh( (xLi - 0.124)   /0.0226)
		  + 0.0155 * tanh( (xLi - 0.105)   /0.029)) ;

   } else if (modelID_ == OCVCathode_CoO2_dualfoil) {

        /*
         *   Line 4580 of dualfoil 5.1
	 * 
	 *       Measured by Oscar Garcia 2001 using Quallion electrodes for
	 *       0.5 < y < 0.99.  Fit revised by Karen Thomas in May 2003 to
         *       match Doyle's fit for y < 0.4 and Garcia's data at larger y.
	 *       Valid for 0 < y < 0.99. Note that capacity fade is found to
	 *       occur experimentally if y goes below 0.5; this is not included
	 *       in the model
	 */
       xV = relExtent_;
        volts = ( 2.16216 + 0.07645 * tanh(30.834 - 54.4806 * xV)
		  + 2.1581  * tanh( 52.294  - 50.294  * xV)
		  - 0.14169 * tanh( 11.0923 - 19.8543 * xV)
		  + 0.2051 *  tanh( 1.4684  - 5.4888  * xV)
		  + 0.2531 *  tanh( (-xV + 0.56478) / 0.1316)
		  - 0.02167 * tanh( ( xV - 0.525)   / 0.006)      );

    } else {
        printf("model not found\n");
        exit(-1);
    }
    /*
     *  interpolate assuming a linear relation. If we get something more complicated we'll have to revise
     */
    volts += dOCVdt * (temp - temperatureBase_);

    return volts;           
}
//===============================================================================================================================
double RSD_OCVmodel::OCV_dvaldExtent(double temp) const
{  
    calcRelExtent();
    throw Electrode_Error("", "not handled");
    return 0.0;
}
//===============================================================================================================================
double RSD_OCVmodel::OCV_dvaldT(double temp) const
{  
    double dOCVdt = 0.0;
    /*
     *  Calculate the relative extent of the model
     */
    calcRelExtent();

    double xLi = 1.0 - relExtent_;
   
    if (temperatureDerivType_ == 0) {
	/*
	 *  Temperature derivative is zero. -> this is a common choice to avoid egregiously incorrect behavior
	 *                                     and when there is a lack of information.
	 */
	return 0.0;
    } else if (temperatureDerivType_ == 1) {
	return 0.0;
    } else if (temperatureDerivType_ == 2) {
	if (temperatureDerivModelType_ == OCVAnode_CONSTANT) {
	    dOCVdt = dvec_[1];
	}

	else if (temperatureDerivModelType_ == OCVCathode_CONSTANT) {
	    dOCVdt = dvec_[1];
	}
	
	else if (temperatureDerivModelType_ == OCVAnode_MCMB2528_dualfoil) {
	    double xLi_2 = xLi * xLi;
	    double xLi_3 = xLi_2 * xLi;
	    double xLi_4 = xLi_3 * xLi;
	    double xLi_5 = xLi_4 * xLi;
	    double xLi_6 = xLi_5 * xLi;
	    double xLi_7 = xLi_6 * xLi;
	    double xLi_8 = xLi_7 * xLi;

	    double A = (0.00527             + 3.29927     * xLi   - 91.79326    * xLi_2 + 1004.91101  * xLi_3 - 5812.27813 * xLi_4 +
		        19329.75490 * xLi_5 - 37147.89470 * xLi_6 + 38379.18127 * xLi_7 - 16515.05308 * xLi_8);

	    double B = (1.0                - 48.09287     * xLi   + 1017.23480  * xLi_2 - 10481.80419  * xLi_3 + 59431.30001 * xLi_4 -
		        195881.64880 * xLi_5 + 374577.31520 * xLi_6 - 385821.16070 * xLi_7 + 165705.85970 * xLi_8);
	    
	    dOCVdt = 1.0E-3 * A / B;
	    return dOCVdt;
	}

	else {
	    throw Electrode_Error("RSD_OCVmodel::OCV_dvaldT()",
				  "Model ID, " + int2str(temperatureDerivModelType_) + ", doesn't have a temperature formulation: " 
				  + OCVTempDerivModel_);
	}
    } else {
	throw Electrode_Error("RSD_OCVmodel::OCV_dvaldT()",
			      "Unknown temperture deriv type: " + int2str(temperatureDerivType_));
    }

    return dOCVdt;
}
//===============================================================================================================================
std::string RSD_OCVmodel::modelName() const
{
    return modelName_;
}
//===============================================================================================================================
double  RSD_OCVmodel::RelExtent() const
{  
    calcRelExtent();
    return relExtent_;
}
//===============================================================================================================================
} // end of namespace
//-------------------------------------------------------------------------------------------------------------------------------

