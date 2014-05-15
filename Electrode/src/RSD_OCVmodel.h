/**
 *  @file RSD_OCVmodel.h
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */
#ifndef RSD_OCVMODEL_H
#define RSD_OCVMODEL_H

#include "Electrode_defs.h"

#include <string>
#include <vector>
#include <map>


//! Model for OCV Override functions

namespace Cantera 
{


// MCMB 2528 graphite measured by Chris Bogatu 2000, 
//           Telcordia and PolyStor materials.
//    Modified May 2003 to match data from Joongpyo Shim
//    for 0.01 < x < 0.99

#define OCVAnode_MCMB2528                101


#define OCVCathode_MCMB2528              301

//!  create a map
/*!
 *
 */
extern void createOCVmodel_map(std::map<int, std::string>& smap);

//==================================================================================

class RSD_OCVmodel
{
  public:
    RSD_OCVmodel(int modelId = -1);

    RSD_OCVmodel(const RSD_OCVmodel& right);

    virtual ~RSD_OCVmodel();

    RSD_OCVmodel& operator=(const  RSD_OCVmodel& right);

    //! Duplicator function for this class
    /*!
     *  @return Returns a duplication of the current state as a pointer to the base class
     */
    virtual RSD_OCVmodel* duplMyselfAsOCVmodel() const;

    
    virtual double OCV_value(double relExtent);

    virtual double OCV_dvaldExtent( double relExtent);

    virtual double OCV_dvaldT(double relExtent);

    virtual std::string modelName();

    int modelID_;
    std::string modelName_;

};


}
#endif
