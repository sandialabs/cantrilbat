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
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef CANTERA_APP
#define CANTERA_APP
#endif

#include "cantera/base/ct_defs.h" 
#include "cantera/base/ctexceptions.h"
#include "m1d_BoundaryCondition.h" 

#include <cstdio>

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d {

//==================================================================================================================================
BoundaryCondition::BoundaryCondition() :
    title_(""),
    lowerLim_(-1.0 * BigNumber),
    upperLim_(BigNumber),
    indepUnits_(""),
    depenUnits_(""),
    step_(0),
    stepMax_(0),
    hasExtendedDependentValue_(2),
    ifuncLowerLim_(0),
    depValLowerLim_(0.0)
{
}
//==================================================================================================================================
BoundaryCondition::~BoundaryCondition()
{
}
//==================================================================================================================================
BoundaryCondition::BoundaryCondition(const BoundaryCondition &right) :
    title_(""),
    lowerLim_(-1.0 * BigNumber),
    upperLim_(BigNumber),
    indepUnits_(""),
    depenUnits_(""),
    step_(0),
    stepMax_(0),
    hasExtendedDependentValue_(2),
    ifuncLowerLim_(0),
    depValLowerLim_(0.0)
{
    *this = right;
}
//==================================================================================================================================
BoundaryCondition& BoundaryCondition::operator=(const BoundaryCondition& right)
{
    if (&right != this) {
	return *this;
    }
    title_ = right.title_;
    lowerLim_ = right.lowerLim_;
    upperLim_ = right.upperLim_;
    indepUnits_ = right.indepUnits_;
    depenUnits_ = right.depenUnits_;
    step_ = right.step_;
    stepMax_ = right.stepMax_;
    hasExtendedDependentValue_ = right.hasExtendedDependentValue_;
    ifuncLowerLim_ = right.ifuncLowerLim_;
    depValLowerLim_ = right.depValLowerLim_;
    return *this;
}
//==================================================================================================================================
void BoundaryCondition::useXML(ZZCantera::XML_Node& bcNode)
{
    if (bcNode.name() != "BoundaryCondition") {
        m1d_Error("BoundaryCondition::useXML)" , "All BoundaryCondition XML nodes have the name, BoundaryCondition");
    }
    err("BoundaryCondition::useXML");
}
//==================================================================================================================================
double BoundaryCondition::value(double indVar, int interval) const
{
    return err("BoundaryCondition::value");
}
//==================================================================================================================================
double BoundaryCondition::valueAtTime(double time, double indVar, int interval) const
{
    return value(indVar, interval);
}
//==================================================================================================================================
double BoundaryCondition::valueAtTime_full(double time, const double* const solnVecNode, int interval) const
{
    return valueAtTime(time, solnVecNode[0], interval);
}
//==================================================================================================================================
double BoundaryCondition::nextStep() const
{
    return err("BoundaryCondition::nextStep");
}
//==================================================================================================================================
void BoundaryCondition::resetSteps()
{
    step_ = 0;
}
//===========================================================================================================
void BoundaryCondition::writeProfile() const
{
}
//===========================================================================================================
double BoundaryCondition::lowerLimit() const
{
    return lowerLim_;
}
//===========================================================================================================
double BoundaryCondition::upperLimit() const
{
    return upperLim_;
}
//=========================================================================================================== 
std::string BoundaryCondition::indepUnits() const 
{
    return indepUnits_;
}
//===========================================================================================================
std::string BoundaryCondition::depenUnits() const 
{
    return depenUnits_;
}
//===========================================================================================================
std::string BoundaryCondition::title() const
{
    return title_;
}
//==================================================================================================================================
void BoundaryCondition::setTitle(const std::string& name)
{
    title_ = name;
}
//==================================================================================================================================
void BoundaryCondition::setLowerLimitBoundsTreament(int ifunc, double depVal)
{
    ifuncLowerLim_ = ifunc;
    depValLowerLim_ = depVal;
}
//==================================================================================================================================
void BoundaryCondition::setLowerLimit(double indVal)
{
    lowerLim_ = indVal;
}
//==================================================================================================================================
void BoundaryCondition::setUpperLimit(double indVal)
{
    upperLim_ = indVal;
}
//==================================================================================================================================
void BoundaryCondition::setIndepUnits(const std::string& unitString)
{
    indepUnits_ = unitString;
}
//==================================================================================================================================
void BoundaryCondition::setDepenUnits(const std::string& unitString)
{
    depenUnits_ = unitString;
}
//===========================================================================================================
int BoundaryCondition::findStep(double indVar, int interval) const
{
    err("BoundaryCondition::findStep");
    return -1;
}
//===========================================================================================================
double BoundaryCondition::err(std::string msg) const
{
    throw m1d_Error("BoundaryCondition Base Class\n", "**** Method " + msg + " not implemented\n");
    return 0.0;
}
//==================================================================================================================================
//==================================================================================================================================
//==================================================================================================================================
BCconstant::BCconstant(double value, std::string titleName, std::string indepUnits, std::string depenUnits) :
    BoundaryCondition(),
    dependentVal_(value)
{
    setTitle(titleName);
    setIndepUnits(indepUnits);
    setDepenUnits(depenUnits);
    stepMax_ = 1;
}
//==================================================================================================================================
BCconstant::~BCconstant()
{
}
//==================================================================================================================================
double BCconstant::value(double indVar, int interval) const
{
    return dependentVal_;
}
//==================================================================================================================================
double BCconstant::nextStep() const
{
    if (step_ < stepMax_) {
	return upperLim_;
    } else {
	return -1;
    }
}
//==================================================================================================================================
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////         class BCsteptable        ///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//==================================================================================================================================
// Constructor taking two vectors of floats representing the independent variable and the dependent  variable
BCsteptable::BCsteptable(const ZZCantera::vector_fp& indValue, const ZZCantera::vector_fp& depValue, 
                         const ZZCantera::vector_fp& compareValue, const std::string& titleName, const std::string& indepUnits, 
                         const std::string& depenUnits) :
    BoundaryCondition(),
    indepVals_(indValue),
    depenVals_(depValue),
    compareVals_(compareValue)
{
    setTitle(titleName);
    setIndepUnits(indepUnits);
    setDepenUnits(depenUnits);

    if (indepVals_.size() != depenVals_.size()) {
        if (indepVals_.size() != depenVals_.size() + 1) {
            throw m1d_Error("BCsteptable constructor\n", "**** indepVals_ and depenVals_ unequal size\n");
        } else {
            hasExtendedDependentValue_ = 0;
            stepMax_ = depenVals_.size();
        }
    } else {
        hasExtendedDependentValue_ = 2;
        stepMax_ = depenVals_.size() - 1;
    }
    lowerLim_ = indepVals_[0];
    upperLim_ = indepVals_.back();
}
//=====================================================================================================================
// construct from filename
BCsteptable::BCsteptable(std::string filename) :
    BoundaryCondition()
{
    //convert input file name into XML node
    XML_Node* baseNode = get_XML_File(filename);
    //find the "boundaryCondition" node
    std::string targetName = "boundaryCondition";
    XML_Node* bcNode = get_XML_NameID(targetName, "", baseNode);
    //parse this node into class data structures
    useXML(*bcNode);
    close_XML_File(filename);
}
//=====================================================================================================================
//! construct from XMLnode
BCsteptable::BCsteptable(ZZCantera::XML_Node& baseNode) :
        BoundaryCondition()
{
    useXML(baseNode);
}
//=====================================================================================================================
BCsteptable::~BCsteptable()
{
}
//=====================================================================================================================
void BCsteptable::useXML(ZZCantera::XML_Node& bcNode)
{
    if (bcNode.name() != "BoundaryCondition") {
        m1d_Error("BoundaryCondition::useXML)" , "All BoundaryCondition XML nodes have the name, BoundaryCondition");
    }
    if (!bcNode.hasChild("independentVar"))
        throw m1d_Error("BCsteptable::useXML()", "no independentVar XML node.");

    if (!bcNode.hasChild("dependentVar"))
        throw m1d_Error("BCsteptable::useXML()", "no dependentVar XML node.");


    //get independentVar
    //XML_Node& indVarNode = bcNode.child("independentVar");
    ZZctml::getFloatArray(bcNode, indepVals_, true, "", "independentVar");


    //get dependentVar
    //XML_Node& depVarNode = bcNode.child("dependentVar");
    ZZctml::getFloatArray(bcNode, depenVals_, true, "", "dependentVar");
    //getNamedFloatArray(depVarNode, depenVals_, convert, depenUnits_);

    //get compareVar
    if (bcNode.hasChild("compareVar")) {
        //XML_Node& compVarNode = bcNode.child("compareVar");
        //getFloatArray(compVarNode, compareVals_, convert, compareUnits_);
        ZZctml::getFloatArray(bcNode, compareVals_, true, "", "compareVar");
    }

    /*
     *  Some of the data was written with equal numbers of independent and dependent vals.
     *  Need to massage the data so that the independent data has one more data point, because
     *  it is interval data.
     */
    if (indepVals_.size() == depenVals_.size() + 1) {
        stepMax_ = depenVals_.size();
        lowerLim_ = indepVals_[0];
        upperLim_ = indepVals_.back();
    } else if (indepVals_.size() == depenVals_.size()) {
        stepMax_ = indepVals_.size();
        indepVals_.push_back(0.0);
        for (int k = stepMax_; k != 0; --k) {
            indepVals_[k] = indepVals_[k-1];
        }
        indepVals_[0] = 0.0;
        stepMax_ = depenVals_.size();
        lowerLim_ = indepVals_[0];
        upperLim_ = indepVals_.back();
    } else {
        throw m1d_Error("BCsteptable::useXML()", "independent and dependent variable dimension mismatch");
    }
}
//==================================================================================================================================
double BCsteptable::value(double indVar, int interval) const
{
    step_ = findStep(indVar, interval);
    if (step_ < 0) {
        return depValLowerLim_;
    }
    return depenVals_[step_];
}
//==================================================================================================================================
double BCsteptable::nextStep() const
{
    if (step_ < (int) indepVals_.size()) {
        return indepVals_[step_+1];
    } else {
        return -1;
    }
}
//==================================================================================================================================
int BCsteptable::findStep(double indVar, int interval) const
{
    int step = -1;
    if (indepVals_.size() == 0) {
        return 0;
    }
    if (indVar < indepVals_[0]) {
        if (ifuncLowerLim_ == 0) {
            throw m1d_Error("BCsteptable::findStep()",
                            "Out of bounds error with step < 0\n\tProbably because indVar < indepVals_[0]");
        } else if (ifuncLowerLim_ == 1) {
            return 0;
        }
        return -1; 
    } else if (indVar == indepVals_[0]) {
        if (ifuncLowerLim_ != 0) {
            if (interval < 0) {
                return -1;
            }
        }
        return 0;
    }
    for (int i = 0; i < stepMax_; i++) {
        if (interval == i) {
            if (indVar <= indepVals_[i + 1]) {
                return i;
            }
        } else {
            if (indVar < indepVals_[i + 1]) {
                return i;
            }
        }
    }
    if (indVar == indepVals_[stepMax_]) {
        return stepMax_;
    }
    if (hasExtendedDependentValue_ == 2) {
       return stepMax_;
    } else if (hasExtendedDependentValue_ == 1) {
       return std::max(0, stepMax_- 1);
    }
    if (step < 0) {
        throw m1d_Error("BCsteptable::findStep()",
                        "Out of bounds error with step < 0\n\tProbably because indVar > indepVals_[ stepMax_ ] ");
    }
    return step;
}
//==================================================================================================================================
void BCsteptable::writeProfile() const
{
}
//==================================================================================================================================
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////               class BClineartable            ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * This subclass is designed to handle a table of dependent variable boundary conditions that are 
 * to be linearly interpolated between the values given.  For example, if the value pairs (0,0) 
 * and (1,1) are given, the value( 0.5 ) will  return 0.5.
 */
//==================================================================================================================================
// constructor taking vectors of floats
BClineartable::BClineartable(ZZCantera::vector_fp indValue, ZZCantera::vector_fp depValue, 
                             ZZCantera::vector_fp compareValue, std::string titleName,
                             std::string indepUnits, std::string depenUnits) :
    BoundaryCondition(),
    indepVals_(indValue),
    depenVals_(depValue),
    compareVals_(compareValue)
{
    setTitle(titleName);
    setIndepUnits(indepUnits);
    setDepenUnits(depenUnits);
    
    if (indepVals_.size() != depenVals_.size()) {
        throw m1d_Error("BClineartable constructor\n", "**** indepVals_ and depenVals_ unequal size\n");
	
    }
    stepMax_ = indepVals_.size();
    lowerLim_ = indepVals_[0];
    upperLim_ = indepVals_.back();
}
//=====================================================================================================================
BClineartable::BClineartable(std::string filename) :
        BoundaryCondition()
{
    //convert input file name into XML node
    XML_Node* baseNode = get_XML_File(filename);
    //find the "boundaryCondition" node
    std::string targetName = "boundaryCondition";
    XML_Node* bcNode = get_XML_NameID(targetName, "", baseNode);
    //parse this node into class data structures
    useXML(*bcNode);
    close_XML_File(filename);
}
//=====================================================================================================================
BClineartable::BClineartable(ZZCantera::XML_Node& baseNode) :
        BoundaryCondition()
{
    useXML(baseNode);
}
//=====================================================================================================================
BClineartable::~BClineartable()
{
}
//=====================================================================================================================
// fill independent and dependent values from XML_Node
void BClineartable::useXML(ZZCantera::XML_Node& bcNode)
{
    if (bcNode.name() != "BoundaryCondition") {
        m1d_Error("BoundaryCondition::useXML)" , "All BoundaryCondition XML nodes have the name, BoundaryCondition");
    }
    if (!bcNode.hasChild("independentVar")) {
        throw m1d_Error("BClineartable::useXML()", "no independentVar XML node.");
    }
    if (!bcNode.hasChild("dependentVar")) {
        throw m1d_Error("BClineartable::useXML()", "no dependentVar XML node.");
    }
    //bool convert = true;

    //get independentVar
    //XML_Node& indVarNode = bcNode.child("independentVar");
    //getFloatArray(indVarNode, indepVals_, convert, indepUnits_);
    ZZctml::getFloatArray(bcNode, indepVals_, true, "", "independentVar");

    //get dependentVar
    //XML_Node& depVarNode = bcNode.child("dependentVar");
    //getFloatArray(depVarNode, depenVals_, convert, depenUnits_);
    ZZctml::getFloatArray(bcNode, depenVals_, true, "", "dependentVar");

    //get compareVar
    if (bcNode.hasChild("compareVar")) {
        //XML_Node& compVarNode = bcNode.child("compareVar");
        //getFloatArray(compVarNode, compareVals_, convert, compareUnits_);
        ZZctml::getFloatArray(bcNode, compareVals_, true, "", "compareVar");
    }

    //err("BClineartable::useXML");

    if (indepVals_.size() == depenVals_.size()) {
        stepMax_ = depenVals_.size();
        lowerLim_ = indepVals_[0];
        upperLim_ = indepVals_.back();
    } else {
        throw m1d_Error("BClineartable::useXML()", "independent and dependent variable dimension mismatch");
    }
}
//=====================================================================================================================
double BClineartable::value(double indVar, int interval) const
{
    step_ = findStep(indVar, interval);
    if (step_ >= stepMax_ - 1) {
        return depenVals_[step_];
    } else {
        //
        // Do a linear interpolation within the step, step_ .  
        //
        double tmp = depenVals_[step_]
                + (indVar - indepVals_[step_]) * (depenVals_[step_ + 1] - depenVals_[step_])
                        / (indepVals_[step_ + 1] - indepVals_[step_]);
        return tmp;
    }
}
//=====================================================================================================================
double BClineartable::nextStep() const
{
    if (step_ < (int) indepVals_.size()) {
        return indepVals_[step_+1];
    } else {
        return -1;
    }
}
//=====================================================================================================================
int BClineartable::findStep(double indVar, int interval) const
{
    int step = -1;
    if (indVar < indepVals_[0]) {
        throw m1d_Error("BClineartable::findStep()",
                "Out of bounds error with step < 0\n\tProbably because indVar < indepVals_[0]");
    }
    for (int i = 0; i < stepMax_ - 1; i++) {
        if (interval == i) {
            if (indVar <= indepVals_[i + 1]) {
                step = i;
                break;
            }
        } else {
            if (indVar < indepVals_[i + 1]) {
                step = i;
                break;
            }
        }
    }
    if (indVar >= indepVals_[stepMax_ - 1]) {
        step = stepMax_ - 1;
    }
    if (step < 0) {
        throw m1d_Error("BClineartable::findStep()",
                "Out of bounds error with step < 0\n\tProbably because indVar > indepVals_[ stepMax_ ] ");
    }
    return step;
}
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
////////////////////////////////////////////////////////////
// class BCsinusoidal
////////////////////////////////////////////////////////////
//=====================================================================================================================
BCsinusoidal::BCsinusoidal(double baseDepValue, double oscAmplitude, double frequency, std::string titleName,
                           std::string indepUnits, std::string depenUnits) :
    BoundaryCondition(),
    baseDepValue_(baseDepValue),
    oscAmplitude_(oscAmplitude),
    frequency_(frequency)
{
    setTitle(titleName);
    setIndepUnits(indepUnits);
    setDepenUnits(depenUnits);

    stepMax_ = 1; // this will be updated in useXML
    lowerLim_ = 0.0;
    upperLim_ = 1.0; //this will be updated in useXML below
}
//=====================================================================================================================
// construct from filename
BCsinusoidal::BCsinusoidal(std::string filename) :
        BoundaryCondition()
{
    //convert input file name into XML node
    XML_Node* baseNode = get_XML_File(filename);
    //find the "sinusoidalFunction" node
    std::string targetName = "sinusoidalFunction";
    XML_Node* bcNode = get_XML_NameID(targetName, "", baseNode);
    //parse this node into class data structures
    useXML(*bcNode);
    close_XML_File(filename);
}
//=====================================================================================================================
BCsinusoidal::~BCsinusoidal()
{
}
//=====================================================================================================================
// construct from XMLnode
BCsinusoidal::BCsinusoidal(ZZCantera::XML_Node& baseNode) :
    BoundaryCondition()
{
    useXML(baseNode);
}
//=====================================================================================================================
// Fill independent and dependent values from XML_Node
void BCsinusoidal::useXML(ZZCantera::XML_Node& bcNode)
{
    if (bcNode.name() != "BoundaryCondition") {
        m1d_Error("BoundaryCondition::useXML)" , "All BoundaryCondition XML nodes have the name, BoundaryCondition");
    }
    if (!bcNode.hasChild("baseDependentValue"))
        throw m1d_Error("BCsinusoidal::useXML()", "no baseDependentValue XML node.");

    if (!bcNode.hasChild("oscillationAmplitude"))
        throw m1d_Error("BCsinusoidal::useXML()", "no oscillationAmplitude XML node.");

    if (!bcNode.hasChild("frequency"))
        throw m1d_Error("BCsinusoidal::useXML()", "no frequency XML node.");

    //get base amplitude
    baseDepValue_ = ZZctml::getFloat(bcNode, "baseDependentValue");

    //get oscillation amplitude
    oscAmplitude_ = ZZctml::getFloat(bcNode, "oscillationAmplitude");

    //get frequency of the oscillation
    frequency_ = ZZctml::getFloat(bcNode, "frequency");

    //set 100 steps per period
    deltaT_ = 1. / stepsPerPeriod_ / frequency_;
    lowerLim_ = 0.0;
    //set ten periods
    upperLim_ = periodsPerRun_ / frequency_;

    stepMax_ = periodsPerRun_ * stepsPerPeriod_ + 1;

}
//=====================================================================================================================
double BCsinusoidal::value(double indVar, int interval) const
{
    return baseDepValue_ + oscAmplitude_ * sin(2.0 * Pi * indVar * frequency_);
}
//=====================================================================================================================
// Return the next value for the independent variable at which the nature of the boundary condition changes.
/*
 * This is designed to guide grid generation and time stepping
 */
double BCsinusoidal::nextStep() const
{
    return deltaT_ * (step_+1);
}
//=====================================================================================================================
// Check to see which step we are at.
int BCsinusoidal::findStep(double indVar, int interval) const
{
    step_ = ceil(indVar / deltaT_);
    return step_;
}
//=====================================================================================================================
void BCsinusoidal::writeProfile() const
{
}
//=====================================================================================================================
void BCsinusoidal::setFrequency(double frequency)
{
    frequency_ = frequency;
    //update other values
    deltaT_ = 1. / stepsPerPeriod_ / frequency_;
    upperLim_ = periodsPerRun_ / frequency_;
}
//=====================================================================================================================
// specify the number of steps per period
void BCsinusoidal::setStepsPerPeriod(double stepsPerPeriod)
{
    stepsPerPeriod_ = stepsPerPeriod;
    //update other values
    deltaT_ = 1. / stepsPerPeriod_ / frequency_;
    stepMax_ = periodsPerRun_ * stepsPerPeriod_ + 1;
}
//=====================================================================================================================
// specify the number of periods to be run
void BCsinusoidal::setPeriodsPerRun(double periodsPerRun)
{
    periodsPerRun_ = periodsPerRun;
    //update other values
    stepMax_ = periodsPerRun_ * stepsPerPeriod_ + 1;
    upperLim_ = periodsPerRun_ / frequency_;
}
//=====================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
