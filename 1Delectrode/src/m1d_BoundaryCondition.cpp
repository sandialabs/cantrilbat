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
// Copyright 2010 Sandia National Laboratories
#ifndef CANTERA_APP
#define CANTERA_APP
#endif

#include "cantera/base/ct_defs.h" 
#include "cantera/base/ctexceptions.h"
#include "m1d_BoundaryCondition.h" 

#include <cstdio>

using namespace Cantera;

namespace m1d {

//============================================================================================================
BoundaryCondition::BoundaryCondition() :
    title_(""),
    lowerLim_(-1.0 * BigNumber),
    upperLim_(BigNumber),
    indepUnits_(""),
    depenUnits_(""),
    step_(0),
    stepMax_(0)
{
}
//============================================================================================================
BoundaryCondition::~BoundaryCondition()
{
}
//============================================================================================================
BoundaryCondition::BoundaryCondition(const BoundaryCondition &right) :
    title_(""),
    lowerLim_(-1.0 * BigNumber),
    upperLim_(BigNumber),
    indepUnits_(""),
    depenUnits_(""),
    step_(0),
    stepMax_(0)
{
    *this = right;
}
//============================================================================================================
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
    return *this;
}
//============================================================================================================
// Return the dependent variable value given the independent variable argument
/*
 *   @param indVar  Independent variable
 *   @param interval If greater than zero, then checking is done on the interval specified
 */
double BoundaryCondition::value(double indVar, int interval)
{
    return err("BoundaryCondition::value");
}
//============================================================================================================
double BoundaryCondition::valueAtTime(double time, double indVar, int interval)
{
    return value(indVar, interval);
}
//============================================================================================================
double BoundaryCondition::valueAtTime_full(double time, double* solnVecNode, int interval)
{
    return valueAtTime(time, solnVecNode[0], interval);
}
//===========================================================================================================
double BoundaryCondition::nextStep()
{
    return err("BoundaryCondition::nextStep");
}
//===========================================================================================================
void BoundaryCondition::resetSteps()
{
    step_ = 0;
}
//===========================================================================================================
void BoundaryCondition::writeProfile()
{
}
//===========================================================================================================
double BoundaryCondition::lowerLimit()
{
    return lowerLim_;
}
//===========================================================================================================
double BoundaryCondition::upperLimit()
{
    return upperLim_;
}
//=========================================================================================================== 
std::string BoundaryCondition::indepUnits()
{
    return indepUnits_;
}
//===========================================================================================================
std::string BoundaryCondition::depenUnits()
{
    return depenUnits_;
}
//===========================================================================================================
std::string BoundaryCondition::title()
{
    return title_;
}
//===========================================================================================================
void BoundaryCondition::setTitle(std::string name)
{
    title_ = name;
}
//===========================================================================================================
void BoundaryCondition::setLowerLimit(double indVal)
{
    lowerLim_ = indVal;
}
//===========================================================================================================
void BoundaryCondition::setUpperLimit(double indVal)
{
    upperLim_ = indVal;
}
//===========================================================================================================
void BoundaryCondition::setIndepUnits(std::string unitString)
{
    indepUnits_ = unitString;
}
//===========================================================================================================
void BoundaryCondition::setDepenUnits(std::string unitString)
{
    depenUnits_ = unitString;
}
//===========================================================================================================
int BoundaryCondition::findStep(double indVar, int interval)
{
    err("BoundaryCondition::findStep");
    return -1;
}
//===========================================================================================================
double BoundaryCondition::err(std::string msg) const
{
    throw CanteraError("BoundaryCondition Base Class\n", "**** Method " + msg + " not implemented\n");
    return 0.0;
}
//============================================================================================================
BCconstant::BCconstant(double value, std::string titleName, std::string indepUnits,
		       std::string depenUnits) :
    BoundaryCondition(),
    dependentVal_(value)
{
    setTitle(titleName);
    setIndepUnits(indepUnits);
    setDepenUnits(depenUnits);
        stepMax_ = 1;
}
//============================================================================================================
BCconstant::~BCconstant()
{
}
//============================================================================================================
double BCconstant::value(double indVar, int interval)
{
    return dependentVal_;
}
//============================================================================================================
double BCconstant::nextStep()
{
    if (step_ < stepMax_) {
	++step_;
	return upperLim_;
    } else {
	return -1;
    }
}
//============================================================================================================
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////         class BCsteptable        ///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * This subclass is designed to handle a table of 
 * dependent variable boundary conditions that are 
 * constant until the indicated value of the 
 * independent variable, at which point there is a
 * step change to a new value.  For example, if 
 * the (ind,dep) value pairs (0.5,0.0), (1.0,0.5)
 * and (2., 1.) are given, then value( 0.25 ) will 
 * return 0.0, value( 0.75 ) will return 0.5 and
 * value( 1.25 ) will return 1.0.
 */
//=====================================================================================================================
// Constructor taking two vectors of floats representing the independent variable and the dependent 
// variable
BCsteptable::BCsteptable(vector_fp indValue, vector_fp depValue, vector_fp compareValue, std::string titleName,
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
        throw CanteraError("BCsteptable constructor\n", "**** indepVals_ and depenVals_ unequal size\n");

    }
    stepMax_ = indepVals_.size();
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
BCsteptable::BCsteptable(XML_Node& baseNode) :
        BoundaryCondition()
{
    useXML(baseNode);
}
//=====================================================================================================================
BCsteptable::~BCsteptable()
{
}
//=====================================================================================================================
// fill independent and dependent values from XML_Node
void BCsteptable::useXML(XML_Node& bcNode)
{

    if (!bcNode.hasChild("independentVar"))
        throw CanteraError("BCsteptable::useXML()", "no independentVar XML node.");

    if (!bcNode.hasChild("dependentVar"))
        throw CanteraError("BCsteptable::useXML()", "no dependentVar XML node.");

    bool convert = true;

    //get independentVar
    XML_Node& indVarNode = bcNode.child("independentVar");
    ctml::getFloatArray(bcNode, indepVals_, true, "", "independentVar");


    //get dependentVar
    XML_Node& depVarNode = bcNode.child("dependentVar");
    ctml::getFloatArray(bcNode, depenVals_, true, "", "dependentVar");
    //getNamedFloatArray(depVarNode, depenVals_, convert, depenUnits_);

    //get compareVar
    if (bcNode.hasChild("compareVar")) {
        XML_Node& compVarNode = bcNode.child("compareVar");
        //getFloatArray(compVarNode, compareVals_, convert, compareUnits_);
        ctml::getFloatArray(bcNode, compareVals_, true, "", "compareVar");
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
        throw CanteraError("BCsteptable::useXML()",
                "independent and dependent variable dimension mismatch");
    }

}
//===========================================================================================================
// Return the dependent variable value given
// the independent variable argument
/*
 *   @param indVar   Independent variable
 *   @param interval If greater than zero, then checking is done on the interval specified
 *                   Also ties, i.e. numbers on the boundary go to the interval value.
 */
double BCsteptable::value(double indVar, int interval)
{
    step_ = findStep(indVar, interval);
    return depenVals_[step_];
}
//=====================================================================================================================
// return the next value for the independent variable at
// which the nature of the boundary condition changes.
/*
 * This is designed to guide grid generation and time stepping
 */
double BCsteptable::nextStep()
{
    if (step_ < (int) indepVals_.size()) {
        return indepVals_[step_++];
    } else {
        return -1;
    }
}
//=====================================================================================================================
// check to see which step we are at.
/*
 * Loop through the boundary condition independent variables.
 * The current step is the smallest value in indepVals_[ i ]
 * that indVar is less than.
 */
int BCsteptable::findStep(double indVar, int interval)
{
    int step = -1;
    if (indVar < indepVals_[0]) {
        throw CanteraError("BCsteptable::findStep()",
                "Out of bounds error with step < 0\n\tProbably because indVar < indepVals_[0]");
    }
    for (int i = 0; i < stepMax_; i++) {
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
    if (indVar == indepVals_[stepMax_]) {
        step = stepMax_;
    }
    if (step < 0) {
        throw CanteraError("BCsteptable::findStep()",
                "Out of bounds error with step < 0\n\tProbably because indVar > indepVals_[ stepMax_ ] ");
    }
#ifdef DEBUG_HKM
    if (step != 0) {
        printf("WE ARE HERE STEP = %d\n", step);
    }
#endif
    return step;
}
//=====================================================================================================================
// Write out the profile in tabular format.
void BCsteptable::writeProfile()
{

}
//=====================================================================================================================
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////               class BClineartable            /////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * This subclass is designed to handle a table of dependent variable boundary conditions that are 
 * to be linearly interpolated between the values given.  For example, if the value pairs (0,0) 
 * and (1,1) are given, the value( 0.5 ) will  return 0.5.
 */
//=====================================================================================================================
// constructor taking vectors of floats
BClineartable::BClineartable(vector_fp indValue, vector_fp depValue, vector_fp compareValue, std::string titleName,
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
        throw CanteraError("BClineartable constructor\n", "**** indepVals_ and depenVals_ unequal size\n");
	
    }
    stepMax_ = indepVals_.size();
    lowerLim_ = indepVals_[0];
    upperLim_ = indepVals_.back();
}
//=====================================================================================================================
//! construct from filename
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
// construct from XMLnode
BClineartable::BClineartable(XML_Node& baseNode) :
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
void BClineartable::useXML(XML_Node& bcNode)
{

    if (!bcNode.hasChild("independentVar")) {
        throw CanteraError("BClineartable::useXML()", "no independentVar XML node.");
    }
    if (!bcNode.hasChild("dependentVar")) {
        throw CanteraError("BClineartable::useXML()", "no dependentVar XML node.");
    }
    //bool convert = true;

    //get independentVar
    XML_Node& indVarNode = bcNode.child("independentVar");
    //getFloatArray(indVarNode, indepVals_, convert, indepUnits_);
    ctml::getFloatArray(bcNode, indepVals_, true, "", "independentVar");

    //get dependentVar
    XML_Node& depVarNode = bcNode.child("dependentVar");
    //getFloatArray(depVarNode, depenVals_, convert, depenUnits_);
    ctml::getFloatArray(bcNode, depenVals_, true, "", "dependentVar");

    //get compareVar
    if (bcNode.hasChild("compareVar")) {
        XML_Node& compVarNode = bcNode.child("compareVar");
        //getFloatArray(compVarNode, compareVals_, convert, compareUnits_);
        ctml::getFloatArray(bcNode, compareVals_, true, "", "compareVar");
    }

    //err("BClineartable::useXML");

    if (indepVals_.size() == depenVals_.size()) {
        stepMax_ = depenVals_.size();
        lowerLim_ = indepVals_[0];
        upperLim_ = indepVals_.back();
    } else {
        throw CanteraError("BClineartable::useXML()", "independent and dependent variable dimension mismatch");
    }
}
//=====================================================================================================================
// Return the dependent variable value given
// the independent variable argument
/*
 *   @param indVar   Independent variable
 *   @param interval If greater than zero, then checking is done on the interval specified
 *                   Also ties, i.e. numbers on the boundary go to the interval value.
 */
double BClineartable::value(double indVar, int interval)
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
// return the next value for the independent variable at
// which the nature of the boundary condition changes.
/*
 * This is designed to guide grid generation and time stepping
 */
double BClineartable::nextStep()
{
    if (step_ < (int) indepVals_.size()) {
        return indepVals_[step_++];
    } else {
        return -1;
    }
}
//=====================================================================================================================
// check to see which step we are at.
/*
 * We need to find the interval that we are in between
 * for interpolation purposes.
 */
int BClineartable::findStep(double indVar, int interval)
{
    int step = -1;
    if (indVar < indepVals_[0]) {
        throw CanteraError("BClineartable::findStep()",
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
        throw CanteraError("BClineartable::findStep()",
                "Out of bounds error with step < 0\n\tProbably because indVar > indepVals_[ stepMax_ ] ");
    }
#ifdef DEBUG_HKM
    if (step != 0) {
	printf("WE ARE HERE STEP = %d\n", step);
    }
#endif
    return step;
}
//=====================================================================================================================
////////////////////////////////////////////////////////////
// class BCsinusoidal
////////////////////////////////////////////////////////////

/**
 * This subclass is designed to provide a sinusoidally 
 * varying boundary condition.  A base value (voltage or current) 
 * is defined as well as an amplitude and frequency for the
 * oscillitory component.
 */

// constructor taking vectors of floats
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
BCsinusoidal::BCsinusoidal(XML_Node& baseNode) :
        BoundaryCondition()
{
    useXML(baseNode);
}
//=====================================================================================================================
// Fill independent and dependent values from XML_Node
void BCsinusoidal::useXML(XML_Node& bcNode)
{

    if (!bcNode.hasChild("baseDependentValue"))
        throw CanteraError("BCsinusoidal::useXML()", "no baseDependentValue XML node.");

    if (!bcNode.hasChild("oscillationAmplitude"))
        throw CanteraError("BCsinusoidal::useXML()", "no oscillationAmplitude XML node.");

    if (!bcNode.hasChild("frequency"))
        throw CanteraError("BCsinusoidal::useXML()", "no frequency XML node.");

    //get base amplitude
    baseDepValue_ = getFloat(bcNode, "baseDependentValue");

    //get oscillation amplitude
    oscAmplitude_ = getFloat(bcNode, "oscillationAmplitude");

    //get frequency
    frequency_ = getFloat(bcNode, "frequency");

    //err("BCsinusoidal::useXML");

    //set 100 steps per period
    deltaT_ = 1. / stepsPerPeriod_ / frequency_;
    lowerLim_ = 0.0;
    //set ten periods
    upperLim_ = periodsPerRun_ / frequency_;

    stepMax_ = periodsPerRun_ * stepsPerPeriod_ + 1;

}
//=====================================================================================================================
// Return the dependent variable value given
// the independent variable argument
double BCsinusoidal::value(double indVar, int interval)
{
    return baseDepValue_ + oscAmplitude_ * sin(indVar / (2 * Pi * frequency_));
}
//=====================================================================================================================
// Return the next value for the independent variable at which the nature of the boundary condition changes.
/*
 * This is designed to guide grid generation and time stepping
 */
double BCsinusoidal::nextStep()
{
    return deltaT_ * step_++;
}
//=====================================================================================================================
// Check to see which step we are at.
/*
 * Loop through the boundary condition independent variables.
 * The current step is the smallest value in indepVals_[ i ]
 * that indVar is less than.
 */
int BCsinusoidal::findStep(double indVar, int interval)
{
    step_ = ceil(indVar / deltaT_);
    return step_;
}
//=====================================================================================================================
// Write out the profile in tabular format.
void BCsinusoidal::writeProfile()
{

}
//=====================================================================================================================
// specify the frequency
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
}//namespace Cantera

//=====================================================================================================================
