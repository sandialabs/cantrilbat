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

namespace Cantera
{
//=============================================================================================================================
//  Here we make a map of all ModelID's to string associations that we support
/*
 *
 */
void createOCVmodel_map(std::map<int, std::string>& smap)
{
   smap[OCVAnode_MCMB2528]         = "MCMB2528";
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
void RSD_OCVmodel::calcRelExtent() const
{
    // We do slightly more than is necessary so that models may have the full mole fraction vector
    solidPhaseModel_->getMoleFractions(&xMF_[0]);
    relExtent_ = xMF_[kSpecies_DoD_];
}
//=============================================================================================================================== 
double RSD_OCVmodel::OCV_value() const
{
    calcRelExtent();
    double volts ;
    if (modelID_ == OCVAnode_CONSTANT) {
       volts = dvec_[0];

    } else if (modelID_ == 101) {
                
       //  MCMB 2528 graphite measured by Chris Bogatu 2000, 
       //  Telcordia and PolyStor materials.
       double xLi = 1.0 - relExtent_;
       volts = (  0.124   + 1.5 * exp(-150.0 * xLi) 
                   + 0.0351 * tanh( (xLi - 0.286) / 0.083)
                   - 0.0045 * tanh( (xLi - 0.90)  / 0.119)
                   - 0.035  * tanh( (xLi - 0.99)  / 0.05 )
                   - 0.0147 * tanh( (xLi - 0.50)  / 0.034)
                   - 0.102  * tanh( (xLi - 0.194) / 0.142)
                   - 0.022  * tanh( (xLi - 0.98 ) / 0.0164)
                   - 0.011  * tanh( (xLi - 0.124) / 0.0226)
                   + 0.0155 * tanh( (xLi - 0.105) / 0.029));

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
//===============================================================================================================================

