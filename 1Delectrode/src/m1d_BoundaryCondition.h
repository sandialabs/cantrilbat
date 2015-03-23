/**
 * file m1d_BoundaryCondition.h
 * 
 * Header for class BoundaryCondition and subclasses.
 */

/*  $Author: hkmoffa $
 *  $Revision: 540 $
 *  $Date: 2013-02-27 15:18:26 -0700 (Wed, 27 Feb 2013) $
 *
 */
// Copyright 2010 Sandia National Laboratories
#ifndef M1D_BOUNDARYCONDITION
#define M1D_BOUNDARYCONDITION

#ifndef CANTERA_APP
#define CANTERA_APP
#endif

#include "cantera/base/ct_defs.h" 
#include "cantera/base/ctexceptions.h"

#include "cantera/base/xml.h"
#include "cantera/base/ctml.h"
using namespace ctml;
using namespace Cantera;

namespace m1d {

/** 
 * The BoundaryCondition class provides a single scalar dependent 
 * variable as a function a single independent variable.  This 
 * is suitable for either spatial or temporal boundary conditions.
 *
 * This is an abstract base class.  Subclasses provide specific
 * relationships between the dependent variable and independent variable.
 *
 * Boundary conditions representing fluxes are always added into the residual equations.
 * The residual equations always have their time derivative terms representing the accretion of 
 * the conserved quantity as positive. What this works out to is that the boundary conditions
 * representing the flux out of the domain is added to this residual.
 * This is true no matter if we are on the X = 0 side of the domain or the X = +TOP side of the domain.
 * 
 */
class BoundaryCondition
{

public:

    //!  Constructor.
    BoundaryCondition();

    //! Destructor.
    virtual ~BoundaryCondition();

    //!  Copy Constructor.
    BoundaryCondition(const BoundaryCondition &right);

    //! Assignment operator
    BoundaryCondition& operator=(const BoundaryCondition& right);

    //! Return the dependent variable value given
    //! the independent variable argument
    /*!
     *   @param indVar  Independentvariable
     *   @param interval If greater than zero, then checking is done on the interval specified
     *                   Also ties, i.e. numbers on the boundary go to the interval value.
     */
    virtual double value(double indVar, int interval = -1);

    //! Returns the dependent variable given the time and the independent variable
    /*!
     *   @param time   Value of the time
     *   @param indVar independent variable
     *   @param interval  If greater than zero, then checking is done on the interval specified
     *                    Also ties, i.e. numbers on the boundary go to the interval value.
     *
     *   @result   returns the value as a double.  
     */
    virtual double valueAtTime(double time, double indVar, int interval = -1);

    //! Returns the dependent variable given the time and the vector of solution values present at the local node
    /*!
     *   @param time   Value of the time
     *   @param solnVecNode  Vector of solution values at the current node
     *   @param interval  If greater than zero, then checking is done on the interval specified
     *                    Also ties, i.e. numbers on the boundary go to the interval value.
     *
     *   @result   returns the value as a double.  
     */
    virtual double valueAtTime_full(double time, double *solnVecNode, int interval = -1);

    //! Return the next value for the independent variable, i.e., time, at
    //! which the nature of the boundary condition changes.
    /*!
     *     This is designed to guide grid generation and time stepping
     */
    virtual double nextStep();

    //! Reset the step counter to zero
    void resetSteps();

    //! Write out the profile in tabular format.
    virtual void writeProfile();

    //! Lower limit of dependent variable for which BC applies
    double lowerLimit();

    //! Upper limit of dependent variable for which BC applies
    double upperLimit();

    //! The string defining the independent variable units
    std::string indepUnits();

    //! The string defining the dependent variable units
    std::string depenUnits();

    //! Return the title or name of boundary condition
    std::string title();

    //! Set title or name of boundary condition
    void setTitle(std::string name);

    //! Set lower limit of independent variable for which BC applies
    void setLowerLimit(double indVal);

    //! Set upper limit of independent variable for which BC applies
    void setUpperLimit(double indVal);

    //! set the string defining the independent variable units
    void setIndepUnits(std::string unitString);

    //! set the string defining the dependent variable units
    void setDepenUnits(std::string unitString);

protected:
    //! title or name of boundary condition
    std::string title_;

    //! lower limit of dependent variable for which BC applies
    double lowerLim_;

    //! upper limit of dependent variable for which BC applies
    double upperLim_;

    //! units string for the independent variable
    std::string indepUnits_;

    //! units string for the dependent variable
    std::string depenUnits_;

    //! current step in a sequence of values
    int step_;

    //! if step_ exceeds the number of steps that were input,
    //! then an out of bounds error should be generated
    int stepMax_;

    //! Calculate which interval of the independent variable we are in.
    /*!
     * If the independent variable argument exceeds the
     * current range, then increment the step counter and
     * check that this has not gone out of bounds.
     *
     *  @param interval   If positive, then ties goes to the value of the interval.
     *                    So, if indVar is at a boundary, then the interval chose is
     *                    equal to the value of the interval variable.
     */
    virtual int findStep(double indVar, int interval = -1);

    //! Error routine
    /*!
     * Throw an exception if an unimplemented method of this class is
     * invoked.
     *
     *  @param msg  Descriptive message string to add to the error report
     *
     *  @return  returns a double, though we will never get there
     */
    double err(std::string msg) const;
};

////////////////////////////////////////////////////////////
// class BCconstant
////////////////////////////////////////////////////////////

class BCconstant: public BoundaryCondition
{

public:

    BCconstant(double value = 0.0, std::string titleName = "BCtitle", std::string indepUnits = "unknownUnits",
               std::string depenUnits = "unknownUnits");

    virtual ~BCconstant();

    //! Return the dependent variable value given the independent variable argument
    /*!
     *   @param indVar  Independentvariable
     *   @param interval If greater than zero, then checking is done on the interval specified
     *                   Also ties, i.e. numbers on the boundary go to the interval value.
     */
    virtual double value(double indVar, int interval = -1);

    //! return the next value for the independent variable at
    //! which the nature of the boundary condition changes.
    /**
     * This is designed to guide grid generation and time stepping.
     * For this constant BC subclass, this provides no real information.
     */
    virtual double nextStep();

protected:

    //!      The fixed dependent variable
    double dependentVal_;
};

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
class BCsteptable: public BoundaryCondition
{
public:

    BCsteptable(vector_fp indValue, vector_fp depValue, vector_fp compareVals_, std::string titleName = "BCsteptable",
                std::string indepUnits = "unknownUnits", std::string depenUnits = "unknownUnits");

    //! construct from filename
    BCsteptable(std::string filename);

    //! construct from XMLnode
    BCsteptable(XML_Node& node);

    //! destructor
    virtual ~BCsteptable();

    //! fill independent and dependent values from XML_Node
    void useXML(XML_Node& node);

    //! Return the dependent variable value given
    //! the independent variable argument
    /*!
     *   @param indVar   Independent variable
     *   @param interval If greater than zero, then checking is done on the interval specified
     *                   Also ties, i.e. numbers on the boundary go to the interval value.
     */
    virtual double value(double indVar, int interval = -1);

    //! return the next value for the independent variable at
    //! which the nature of the boundary condition changes.
    /**
     * This is designed to guide grid generation and time stepping
     */
    virtual double nextStep();

    virtual void writeProfile();

protected:

    //! vector of independent variable values at which
    //! the dependent variable may change value
    vector_fp indepVals_;

    //! vector of dependent variable values appropriate
    //! for time/space after the corresponding indepVals_
    vector_fp depenVals_;

    //! vector of variable values for comparison purposes.
    //! For example, if current is input, these might be measured voltages
    vector_fp compareVals_;

    //! units string for a variable used for comparison purposes
    std::string compareUnits_;

    //! Calculate which interval of the independent variable we are in.
    /*!
     * If the independent variable argument exceeds the
     * current range, then increment the step counter and
     * check that this has not gone out of bounds.
     */
    virtual int findStep(double indVar, int interval = -1);

};

//! This subclass is designed to handle a table of  dependent variable boundary conditions that are
//! to be linearly interpolated between the values given.
/*!
 *   This implicitly implies that the function is continuous, and the derivative is piecewise continuous.
 *
 * For example, if the value pairs (0,0)
 * and (1,1) are given, the value( 0.5 ) will 
 * return 0.5.
 */
class BClineartable: public BoundaryCondition
{

public:

    BClineartable(vector_fp indValue, vector_fp depValue, vector_fp compareVals_, std::string titleName = "BClineartable",
                  std::string indepUnits = "unknownUnits", std::string depenUnits = "unknownUnits");

    //! construct from filename
    BClineartable(std::string filename);

    //! construct from XMLnode
    BClineartable(XML_Node& node);

    //! destructor
    virtual ~BClineartable();

    //! fill independent and dependent values from XML_Node
    void useXML(XML_Node& node);

    //! Return the dependent variable value given
    //! the independent variable argument
    /*!
     *   @param indVar   Independent variable
     *   @param interval If greater than zero, then checking is done on the interval specified.
     *                   Also ties, i.e. numbers on the boundary go to the interval value.
     */
    virtual double value(double indVar, int interval = -1);

    //! return the next value for the independent variable at
    //! which the nature of the boundary condition changes.
    /**
     * This is designed to guide grid generation and time stepping
     */
    virtual double nextStep();

protected:

    //! vector of indepedent variable values at which
    //! the dependent variable may change value
    vector_fp indepVals_;

    //! vector of depedent variable values appropriate
    //! for time/space after the corresponding indepVals_
    vector_fp depenVals_;

    //! vector of variable values for comparison purposes
    //! For example, if current is input, these might be measured voltages
    vector_fp compareVals_;

    //! units string for a variable used for comparison purposes
    std::string compareUnits_;

    //! Calculate which interval of the independent variable we are in.
    /*!
     * If the independent variable argument exceeds the
     * current range, then increment the step counter and
     * check that this has not gone out of bounds.
     *
     *  @param interval   If positive, then ties goes to the value of the interval.
     *                    So, if indVar is at a boundary, then the interval chose is
     *                    equal to the value of the interval variable.
     */
    virtual int findStep(double indVar, int interval = -1);

};

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
class BCsinusoidal: public BoundaryCondition
{

public:

    BCsinusoidal(double baseDepValue, double oscAmplitude, double frequency, std::string titleName = "BCsinusoidal",
                 std::string indepUnits = "unknownUnits", std::string depenUnits = "unknownUnits");

    //! construct from filename
    BCsinusoidal(std::string filename);

    //! Construct from XMLnode
    BCsinusoidal(XML_Node& node);

    //! Destructor
    virtual ~BCsinusoidal();

    //! fill independent and dependent values from XML_Node
    void useXML(XML_Node& node);

    //! Return the dependent variable value given
    //! the independent variable argument
    /*!
     *   @param indVar  Independent variable
     *   @param interval If greater than zero, then checking is done on the interval specified
     *                   Also ties, i.e. numbers on the boundary go to the interval value.
     */
    virtual double value(double indVar, int interval = -1);

    //! return the next value for the independent variable at
    //! which the nature of the boundary condition changes.
    /**
     * This is designed to guide grid generation and time stepping
     */
    virtual double nextStep();

    virtual void writeProfile();

    //! specify the number of steps per period
    void setStepsPerPeriod(double stepsPerPeriod);

    //! specify the number of periods to be run
    void setPeriodsPerRun(double periodsPerRun);

    //! specify the number of periods to be run
    void setFrequency(double frequency);

protected:

    //! Base value about which oscillations occur
    double baseDepValue_;

    //! Oscillation amplitude
    double oscAmplitude_;

    //! Frequency of oscillation
    double frequency_;

    //! Define a deltaT independent variable step size based on the frequency
    double deltaT_;

    //! specify the number of steps per period
    double stepsPerPeriod_;

    //! specify the number of periods to be run
    double periodsPerRun_;

    //! Calculate which interval of the independent variable we are in.
    /*!
     * If the independent variable argument exceeds the
     * current range, then increment the step counter and
     * check that this has not gone out of bounds.
     *
     *  @param interval   If positive, then ties goes to the value of the interval.
     *                    So, if indVar is at a boundary, then the interval chose is
     *                    equal to the value of the interval variable.
     */
    virtual int findStep(double indVar, int interval = -1);

};

} //namespace m1d

#endif // M1D_BOUNDARYCONDITION
