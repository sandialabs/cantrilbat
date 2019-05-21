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
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "zuzax/base/ct_defs.h" 
#include "zuzax/base/ctexceptions.h"
#include "m1d_BC_Battery.h" 
#include "m1d_materials.h"

//----------------------------------------------------------------------------------------------------------------------------------
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
double BC_anodeCC::valueAtTime(double time, double voltsAnode, int interval) const
{
    double resistance = resistivity_copper(298.);
    double denom = resistance * thickness_;
    denom = std::max(denom, 1.0E-11);
    double val = (voltsAnode - anodeCC_volts_) / denom;
    return val; 
}
//==================================================================================================================================
//==================================================================================================================================
//==================================================================================================================================

BC_cathodeCC::BC_cathodeCC(double thickness, double extraWireResistance, double electrodeWireCrossSectionalArea,
                           double cathodeCC_phi) :
        BoundaryCondition(),
        cathodeCC_phi_(cathodeCC_phi),
        thickness_(thickness),
        extraWireResistance_(extraWireResistance),
        electrodeWireCrossSectionalArea_(electrodeWireCrossSectionalArea)
{
}
//=====================================================================================================================
BC_cathodeCC::~BC_cathodeCC()
{
}
//=====================================================================================================================
BC_cathodeCC::BC_cathodeCC(const BC_cathodeCC& right) :
   BoundaryCondition(right),
   cathodeCC_phi_(right.cathodeCC_phi_),
   thickness_(right.thickness_),
   extraWireResistance_(right.extraWireResistance_),
   electrodeWireCrossSectionalArea_(right.electrodeWireCrossSectionalArea_)
{
}
//=====================================================================================================================
BC_cathodeCC& BC_cathodeCC::operator=(const BC_cathodeCC& right)
{
   if (&right == this) {
     return *this;
   }
   BoundaryCondition::operator=(right);
   cathodeCC_phi_=right.cathodeCC_phi_;
   thickness_=right.thickness_;
   extraWireResistance_ = right.extraWireResistance_;
   electrodeWireCrossSectionalArea_ = right.electrodeWireCrossSectionalArea_;

   return *this;
}
//=====================================================================================================================
double BC_cathodeCC::valueAtTime(double time, double voltsCathode, int interval) const
{
    double resistivity = resistivity_aluminum(298.);
    double denom = resistivity * thickness_ + extraWireResistance_ * electrodeWireCrossSectionalArea_;
    denom = std::max(denom, 1.0E-11);
    double val = (voltsCathode - cathodeCC_phi_) / denom;
    //  returns the current on a cross-sectional basis
    //  units = amps / m2
    return val;
}
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================

BC_cathodeCCLoad::BC_cathodeCCLoad(double thickness, double extraResistance, double electrodeCrossSectionalArea,
				   double cathodeCC_phi, double resistanceLoad, double voltageLoad) :
        BoundaryCondition(),
    //    cathodeCC_phi_(cathodeCC_phi),
        thickness_(thickness),
        extraResistance_(extraResistance),
        electrodeCrossSectionalArea_(electrodeCrossSectionalArea),
	resistanceLoad_(resistanceLoad),
	voltageLoad_(voltageLoad)
{
}
//=====================================================================================================================
BC_cathodeCCLoad::~BC_cathodeCCLoad()
{
}
//=====================================================================================================================
BC_cathodeCCLoad::BC_cathodeCCLoad(const BC_cathodeCCLoad& right) :
   BoundaryCondition(right),
   //cathodeCC_phi_(right.cathodeCC_phi_),
   thickness_(right.thickness_),
   extraResistance_(right.extraResistance_),
   electrodeCrossSectionalArea_(right.electrodeCrossSectionalArea_),
   resistanceLoad_(right.resistanceLoad_),
   voltageLoad_(right.voltageLoad_)
{
}
//=====================================================================================================================
BC_cathodeCCLoad& BC_cathodeCCLoad::operator=(const BC_cathodeCCLoad& right)
{
   if (&right == this) {
     return *this;
   }
   BoundaryCondition::operator=(right);
  // cathodeCC_phi_ = right.cathodeCC_phi_;
   thickness_ = right.thickness_;
   extraResistance_ = right.extraResistance_;
   electrodeCrossSectionalArea_ = right.electrodeCrossSectionalArea_;
   resistanceLoad_ = right.resistanceLoad_;
   voltageLoad_ = right.voltageLoad_;

   return *this;
}
//=====================================================================================================================
double BC_cathodeCCLoad::valueAtTime(double time, double phiCathode, int interval) const
{
    double resistivity = resistivity_aluminum(298.);
    double denom = resistivity * thickness_ + (extraResistance_ + resistanceLoad_) * electrodeCrossSectionalArea_;
    denom = std::max(denom, 1.0E-11);
    double val = (phiCathode - voltageLoad_) / denom;
    //  returns the current entering the battery at the cathode on a cross-sectional basis
    //  units = amps / m2
    return val;
}
//=====================================================================================================================
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
double BC_heatTransfer::valueAtTime(double time, double tempCollector, int interval) const
{ 
    // This represents the net flux of heat out of the domain.
    //  Units are Watts / m2 / K
    double val = tranCoeff_  * (tempCollector - tempRef_);
    return val;
}
//=====================================================================================================================
}//namespace m1d
//=====================================================================================================================
