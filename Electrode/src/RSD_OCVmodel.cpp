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

#include <cmath>
namespace Cantera
{
//=============================================================================================================
void createOCVmodel_map(std::map<int, std::string>& smap)
{
   smap[OCVAnode_MCMB2528]            = "MCMB2528";
}
//==============================================================================================================================
RSD_OCVmodel::RSD_OCVmodel(int modelID) :
 modelID_(modelID)
{
   modelName_ = modelID_to_stringName_RCD_OCVmodel(modelID);
}
//===============================================================================================================================
RSD_OCVmodel::RSD_OCVmodel(const RSD_OCVmodel& right) :
    modelID_ (right.modelID_),
    modelName_(right.modelName_)
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
    return *this;
}
//===============================================================================================================================
RSD_OCVmodel* RSD_OCVmodel::duplMyselfAsOCVmodel() const
{
    RSD_OCVmodel* pp = new RSD_OCVmodel(*this);
    return pp;    
}
//=============================================================================================================================== 
double RSD_OCVmodel::OCV_value(double relExtent)
{
    double volts ;
    if (modelID_ == 101) {
                
       //  MCMB 2528 graphite measured by Chris Bogatu 2000, 
       //  Telcordia and PolyStor materials.

       volts = (  0.124   + 1.5 * exp(-150.0 * relExtent) 
                   + 0.0351 * tanh( (relExtent - 0.286) / 0.083)
                   - 0.0045 * tanh( (relExtent - 0.90)  / 0.119)
                   - 0.035  * tanh( (relExtent - 0.99)  / 0.05 )
                   - 0.0147 * tanh( (relExtent - 0.50)  / 0.034)
                   - 0.102  * tanh( (relExtent - 0.194) / 0.142)
                   - 0.022  * tanh( (relExtent - 0.98 ) / 0.0164)
                   - 0.011  * tanh( (relExtent - 0.124) / 0.0226)
                   + 0.0155 * tanh( (relExtent - 0.105) / 0.029));

    } else {
        printf("model not found\n");
        exit(-1);
    }
    return volts;           
}
//===============================================================================================================================
double RSD_OCVmodel::OCV_dvaldExtent(double relExtent)
{
    return 0.0;
}
//===============================================================================================================================

double RSD_OCVmodel::OCV_dvaldT(double relExtent)
{
    return 0.0;
}
//===============================================================================================================================
std::string RSD_OCVmodel::modelName()
{
    return modelName_;
}
//===============================================================================================================================
} // end of namespace
//===============================================================================================================================

