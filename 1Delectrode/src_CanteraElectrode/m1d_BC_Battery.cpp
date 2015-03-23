/**
 * file m1d_BoundaryCondition.cpp
 *
 * Source code for class BoundaryCondition and subclasses.
 *
 * The BoundaryCondition class provides a single scalar dependent 
 * variable as a function a single independent variable.  This 
 * is suitable for either spatial or temporal boundary conditions.
 *
 * This file contains methods for subclass of the BoundaryCondition 
 * abstract class.  Methods for the base class are all in
 * exp_BoundaryCondition.h
 */

/*  $Author: hkmoffa $
 *  $Revision: 540 $
 *  $Date: 2013-02-27 15:18:26 -0700 (Wed, 27 Feb 2013) $
 *
 */

#include "cantera/base/ct_defs.h" 
#include "cantera/base/ctexceptions.h"
#include "m1d_BC_Battery.h" 
#include "m1d_materials.h"

namespace m1d {

//=====================================================================================================================

BC_anodeCC::BC_anodeCC(double thickness, double anodeCC_volts) :
        BoundaryCondition(),
        anodeCC_volts_(anodeCC_volts),
        thickness_(thickness)
{
}
//=====================================================================================================================
BC_anodeCC::~BC_anodeCC()
{
}
//=====================================================================================================================
BC_anodeCC::BC_anodeCC(const BC_anodeCC& right) :
   BoundaryCondition(right),
   anodeCC_volts_(right.anodeCC_volts_),
   thickness_(right.thickness_)
{
}
//=====================================================================================================================
BC_anodeCC& BC_anodeCC::operator=(const BC_anodeCC& right)
{
   if (&right == this) {
     return *this;
   }
   BoundaryCondition::operator=(right);
   anodeCC_volts_=right.anodeCC_volts_;
   thickness_=right.thickness_;
   return *this;
}
//=====================================================================================================================
double BC_anodeCC::valueAtTime(double time, double voltsAnode, int interval)
{
    double resistance = resistivity_copper(298.);
    double denom = resistance * thickness_;
    denom = std::max(denom, 1.0E-11);
    double val = (voltsAnode - anodeCC_volts_) / denom;
    return val; 
}
//=====================================================================================================================
//=====================================================================================================================

BC_cathodeCC::BC_cathodeCC(double thickness, double extraResistance, double electrodeCrossSectionalArea,
                           double cathodeCC_volts) :
        BoundaryCondition(),
        cathodeCC_volts_(cathodeCC_volts),
        thickness_(thickness),
        extraResistance_(extraResistance),
        electrodeCrossSectionalArea_(electrodeCrossSectionalArea)
{
}
//=====================================================================================================================
BC_cathodeCC::~BC_cathodeCC()
{
}
//=====================================================================================================================
BC_cathodeCC::BC_cathodeCC(const BC_cathodeCC& right) :
   BoundaryCondition(right),
   cathodeCC_volts_(right.cathodeCC_volts_),
   thickness_(right.thickness_),
   extraResistance_(right.extraResistance_),
   electrodeCrossSectionalArea_(right.electrodeCrossSectionalArea_)
{
}
//=====================================================================================================================
BC_cathodeCC& BC_cathodeCC::operator=(const BC_cathodeCC& right)
{
   if (&right == this) {
     return *this;
   }
   BoundaryCondition::operator=(right);
   cathodeCC_volts_=right.cathodeCC_volts_;
   thickness_=right.thickness_;
   extraResistance_ = right.extraResistance_;
   electrodeCrossSectionalArea_ = right.electrodeCrossSectionalArea_;

   return *this;
}
//=====================================================================================================================
double BC_cathodeCC::valueAtTime(double time, double voltsCathode, int interval)
{
    double resistivity = resistivity_aluminum(298.);
    double denom = resistivity * thickness_ + extraResistance_ * electrodeCrossSectionalArea_;
    denom = std::max(denom, 1.0E-11);
    double val = (voltsCathode - cathodeCC_volts_) / denom;
    return val;
}
//=====================================================================================================================
//=====================================================================================================================
BC_heatTransfer::BC_heatTransfer(double tranCoeff, double tempRef, double electrodeCrossSectionalArea) :
    BoundaryCondition(),
    tempRef_(tempRef),
    tranCoeff_(tranCoeff),
    electrodeCrossSectionalArea_(electrodeCrossSectionalArea)
{
}
//=====================================================================================================================
BC_heatTransfer::~BC_heatTransfer()
{
}
//=====================================================================================================================
BC_heatTransfer::BC_heatTransfer(const BC_heatTransfer& right) :
   BoundaryCondition(right),
   tempRef_(right.tempRef_),
   tranCoeff_(right.tranCoeff_),
   electrodeCrossSectionalArea_(right.electrodeCrossSectionalArea_)
{
}
//=====================================================================================================================
BC_heatTransfer& BC_heatTransfer::operator=(const BC_heatTransfer& right)
{
   if (&right == this) {
     return *this;
   }
   BoundaryCondition::operator=(right);
   tranCoeff_ = right.tranCoeff_;
   tempRef_   = right.tempRef_;
   electrodeCrossSectionalArea_ = right.electrodeCrossSectionalArea_;

   return *this;
}
//=====================================================================================================================
double BC_heatTransfer::valueAtTime(double time, double tempCollector, int interval)
{ 
    // This represents the net flux of heat out of the domain.
    double val = tranCoeff_ * electrodeCrossSectionalArea_ * (tempCollector - tempRef_);
    return val;
}
//=====================================================================================================================
}//namespace m1d
//=====================================================================================================================
