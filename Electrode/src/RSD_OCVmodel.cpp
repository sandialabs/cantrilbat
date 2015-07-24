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

#include <cmath>
#include <cstdio>
#include <cstdlib>

//-----------------------------------------------------------------------------------------------------------------------------
namespace Cantera
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
   smap[OCVAnode_MCMB2528]          = "MCMB2528";
   smap[OCVAnode_MCMB2528_dualfoil] = "MCMB2528_dualfoil";
   smap[OCVCathode_CoO2_dualfoil]   = "CoO2_dualfoil";
}
//==============================================================================================================================
//==============================================================================================================================
//  Constructor 
RSD_OCVmodel::RSD_OCVmodel(int modelID) :
    modelID_(modelID),
    solidPhaseModel_(0),
    kSpecies_DoD_(npos),
    relExtent_(-1.0),
    xMF_(0),
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
    solidPhaseModel_(right.solidPhaseModel_),
    kSpecies_DoD_(right.kSpecies_DoD_),
    relExtent_ (right.relExtent_),
    xMF_(right.xMF_),
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
    //  BEWARE:  This is a shallow copy. It will have to be fixed up outside
    //
    solidPhaseModel_ = right.solidPhaseModel_;
    kSpecies_DoD_ = right.kSpecies_DoD_;
    relExtent_  = right.relExtent_;
    xMF_ = right.xMF_;
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

    if (modelID_ == OCVAnode_CONSTANT) {
           dvec_.resize(1);
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
double RSD_OCVmodel::OCV_value() const
{
    calcRelExtent();
    double volts ;
    double xLi = 1.0 - relExtent_;
    double xV;
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

        xV = 1.0 - relExtent_;
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
    return volts;           
}
//===============================================================================================================================
double RSD_OCVmodel::OCV_dvaldExtent() const
{  
    calcRelExtent();
    throw Electrode_Error("", "not handled");
    return 0.0;
}
//===============================================================================================================================
double RSD_OCVmodel::OCV_dvaldT() const
{  
    calcRelExtent();
    throw Electrode_Error("", "not handled");
    return 0.0;
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

