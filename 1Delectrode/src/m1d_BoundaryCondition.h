/**
 * file m1d_BoundaryCondition.h
 * 
 * Header for class BoundaryCondition and subclasses.
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef M1D_BOUNDARYCONDITION
#define M1D_BOUNDARYCONDITION

#include "m1d_defs.h"

#ifndef CANTERA_APP
#define CANTERA_APP
#endif

#include "cantera/base/ct_defs.h" 
#include "cantera/base/ctexceptions.h"

#include "cantera/base/xml.h"
#include "cantera/base/ctml.h"
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d {

//==================================================================================================================================
//! The BoundaryCondition class provides a single scalar dependent  variable as a function a single independent variable,
//! which is suitable for either spatial or temporal boundary conditions.
/*!
 *  This is an abstract base class.  Subclasses provide specific relationships between the dependent variable 
 *  and the independent variable.
 *
 *  This class can handle two situations. 
 *            One is where the independent variable is time, and we can have step jumps in values
 *            of the dependent variable at specified intervals in time.
 *            For this situation, use the member function value().
 *
 *
 *            Second is where the independent variable is another variable (mostly position, but could be something else)
 *            other than time. Or, the independent variable is the full solution vector at one node.
*             The relationship can have a time component as well, where the relationships
 *            change as a function of time, using intervals as well so that jumps can occur here as well.
 *            For this situation use the membmer function valueAtTime() and valueAtTime_full()
 *
 *
 *  Boundary conditions representing fluxes are always added into the residual equations.
 *  The residual equations always have their time derivative terms representing the accretion of 
 *  the conserved quantity as positive. What this works out to is that the boundary conditions
 *  representing the flux out of the domain is added to this residual.
 *  This is true no matter if we are on the X = 0 side of the domain or the X = +TOP side of the domain.
 *
 *  \todo Write the interval step logic into the base class. Make all classes conform to one set of logic.
 * 
 */
class BoundaryCondition
{
public:
    //! Default Constructor.
    BoundaryCondition();

    //! Destructor.
    virtual ~BoundaryCondition();

    //!  Copy Constructor.
    /*!
     *  @param[in]           right               Object to be copied
     */
    BoundaryCondition(const BoundaryCondition &right);

    //! Assignment operator
    /*!
     *  @param[in]           right               Object to be copied
     *
     *  @return                                  Reference to the current object
     */
    BoundaryCondition& operator=(const BoundaryCondition& right);

    //! Fill independent and dependent values from ZZCantera::XML_Node object named BoundaryCondition
    /*!
     *  XML_Node object contains all of the information for the boundary condition. The name of the node is "BoundaryCondition". 
     *
     *  (virtual from BoundaryCondition)
     *
     *  @param[in]           bcNode                Reference to an XML node where the information for the boundary condition is storred
     */
    virtual void useXML(ZZCantera::XML_Node& bcNode);

    //! Return the dependent variable value given the value of the independent variable argument
    /*!
     *  The independent variable is usually identified as the time. Note, this class can handle step jumps in the value of the 
     *  dependent variable, because there is a concept of a specified interval. This allows one to specify t = ti+ and t = ti- 
     *  values, where ti is a time where there is a step jump in the dependent variable value.
     *  (virtual from BoundaryCondition)
     *
     *  @param[in]           indVar              Independentvariable
     *  @param[in]           interval            If greater than zero, then checking is done on the interval specified.
     *                                           Also ties wrt to the independent variable are satisfied by goint to the interval
     *                                           specified by this value.
     *                                           Defaults to a value of -1.
     *
     *  @return                                  Returns the value at the value of the independent variable
     */
    virtual double value(double indVar, int interval = -1) const;

    //! Returns the dependent variable given the time and the independent variable
    /*!
     *  The independent variable, here, is defined as something other than the time. The boundary condition is now a function
     *  of this independent variable and the time.
     *  Note, this class can handle step jumps in the value of the 
     *  dependent variable, because there is a concept of a specified interval. This allows one to specify t = ti+ and t = ti- 
     *  values, where ti is a time where there is a step jump in the dependent variable value.
     *  (virtual from BoundaryCondition)

     *  @param[in]           time                Value of the time
     *  @param[in]           indVar              Independent variable
     *  @param[in]           interval            If greater than zero, then checking is done on the interval specified.
     *                                           Also ties wrt to the independent variable are satisfied by goint to the interval
     *                                           specified by this value.
     *                                           Defaults to a value of -1.
     *
     *  @result                                  Returns the value as a double.  
     */
    virtual double valueAtTime(double time, double indVar, int interval = -1) const;

    //! Returns the dependent variable given the time and the vector of solution values present at the local node
    /*!
     *  The independent variable, here, is defined as the full solution vector at a node. The boundary condition is now a function
     *  of these independent variables and the time.
     *  Note, this class can handle step jumps in the value of the 
     *  calculated dependent variable, because there is a concept of a specified interval. This allows one to specify t = ti+ and t = ti- 
     *  values, where ti is a time where there is a step jump in the dependent variable value.
     *  (virtual from BoundaryCondition)
     *
     *  @param[in]           time                Value of the time
     *  @param[in]           solnVecNode         Vector of solution values at the current node
     *  @param[in]           interval            If greater than zero, then checking is done on the interval specified.
     *                                           Also ties wrt to the independent variable are satisfied by goint to the interval
     *                                           specified by this value.
     *                                           Defaults to a value of -1.
     *
     *  @result                                  Returns the value as a double.  
     */
    virtual double valueAtTime_full(double time, const double* const solnVecNode, int interval = -1) const;

    //! Calculate which interval of the independent variable we are in.
    /*!
     *  If the independent variable argument exceeds the current interval range, then the interval is incremented to the next value.
     *
     *  If the independent value is below the lowest interval range, then the behavior is set by the setLowerLimitBoundsTreament()
     *  function. The default behavior is to throw an exception.
     *
     *  If the independent value is above the highest interval range, then the behavior is set by the bool hasExtendedDependentValue_.
     *  if  hasExtendedDependentValue_ is true, then The interval is set to stepMax_,
     *  and last value in the vector dependVals_[] is used to set the boundary condition.
     *  If hasExtenededDependentValue_ is false, then an error condition is thrown.
     *  stepMax_ is set to the number of intervals input, which defaults to zero initially. And, the value
     *  of hasExtendedDependentValue_ defaults to true.
     *  If the number of dependsVals_[] values is one less than the input indendVals_[], then  hasExtendedDependentValue_ 
     *  is set to false. 
     *
     *  @param[in]           indVar              Independent variable, usually the time.
     *  @param[in]           interval            If greater than zero, then checking is done on the interval specified.
     *                                           Also ties wrt to the independent variable are satisfied by going to the interval
     *                                           specified by this value.
     *                                           Defaults to a value of -1.
     *
     *  @return                                  Returns the interval within which the indVar lies.
     */
    virtual int findStep(double indVar, int interval = -1) const;

    //! Return the next value for the independent variable, i.e., time, at which the nature of the boundary condition changes.
    /*!
     *  This is designed to guide grid generation and time stepping.
     *  (virtual from BoundaryCondition)
     *
     *  @return                                  Returns the time for the next interval, where the bc changes its nature.
     */
    virtual double nextStep() const;

    //! Reset the step counter to zero
    void resetSteps();

    //! Write out the profile in tabular format.
    /*!
     *  This writes out the nature of the boundary condition to stdout.
     */
    virtual void writeProfile() const;

    //! Lower limit of the independent variable for which BC applies
    /*!
     *  @return                                  Returns the  Lower limit of the independent variable 
     */
    double lowerLimit() const;

    //! Upper limit of the inddependent variable for which BC applies
    /*!
     *  @return                                  Returns the upper limit of independent variable 
     */
    double upperLimit() const;

    //! The string defining the independent variable units
    /*!
     *  @return                                  Returns a string representing the independent variable units
     */
    std::string indepUnits() const;

    //! The string defining the dependent variable units
    /*!
     *  @return                                  Returns a string representing the dependent variable units
     */
    std::string depenUnits() const;

    //! Return the title or name of boundary condition
    /*!
     *  @return                                  Returns the title of the boundary condition
     */
    std::string title() const;

    //! Set title or name of boundary condition
    /*!
     *  @param[in]           name                Title of the boundary condition
     */
    void setTitle(const std::string &name);

    //! Set lower limit of independent variable for which BC applies
    /*!
     *  @param[in]           indVal              Set the lower limit for the independent variable
     */
    void setLowerLimit(double indVal);

    //! Determine what happens when the independent variable is below the lower limit
    /*!
     *  @param[in]           ifunc               Integer describing the treatment. The default treatment without calling this
     *                                           routine is 0.
     *                                             0:  Treat this as a fatal error and throw an error condition
     *                                             1:  Assume the same treatment as the zeroeth interval in the boundary condition.
     *                                                 No extra information is needed. functions of the independent variable are
     *                                                 continuous as the boundary condtion goes below the lower limit.
     *                                             2:  Another interval is assumed below the zeroeth interval with a double
     *                                                 value for the dependent value input in the parameter list.
     *
     *  @param[in]           depVal              Value of the dependent variable to use below the lower limit
     */
    void setLowerLimitBoundsTreament(int ifunc, double depVal = 0.0);

    //! Set upper limit of independent variable for which BC applies
    /*!
     *  @param[in]           indVal              Set the lower limit for the independent variable
     */
    void setUpperLimit(double indVal);

    //! set the string defining the independent variable units
    /*!
     *  @param[in]           unitString          Set the units string for the independent variable
     */
    void setIndepUnits(const std::string& unitString);

    //! set the string defining the dependent variable units
    /*!
     *  @param[in]           unitString          Set the units string for the dependent variable
     */
    void setDepenUnits(const std::string& unitString);

    //! Increment the time interval
    /*!
     *  @return                                  Returns true if there is a valid next time interval. Returns false if not.
     */
    bool incrStep();

    // --------------------------------------------------------- D A T A ----------------------------------------------
protected:
    //! title or name of boundary condition
    std::string title_;

    //! lower limit of the independent variable for which BC applies
    double lowerLim_;

    //! upper limit of the independent variable for which BC applies
    double upperLim_;

    //! units string for the independent variable
    std::string indepUnits_;

    //! units string for the dependent variable
    std::string depenUnits_;

    //! Vector of independent variable values at which the dependent variable may change value
    /*!
     *  Length: numIntervals+1
     */
    ZZCantera::vector_fp indepVals_;

    //! Vector of dependent variable values appropriate for time/space after the corresponding indepVals_
    /*!
     *  These are the values within the intervals. An optional last value is the dependent value
     *  for beyond the last interval.
     *
     *  Length: numIntervals or numIntervals+1
     */
    ZZCantera::vector_fp depenVals_;

    //! Vector of variable dependent values for comparison purposes.
    //! For example, if current is input, these might be measured voltages.
    /*!
     *  Length: numIntervals
     */ 
    ZZCantera::vector_fp compareVals_;

    //! Current step interval in a sequence of values
    /*!
     *  Starts at a value of zero. This represents the value of the current interval.
     */
    mutable int step_;

    //! Maximum value of the interval steps.
    /*!
     *  If step_ exceeds the number of steps that were input, then an out of bounds error should be generated.
     *  Starts at a value of zero.
     *  Value of the last interval. Defaults to 0, if no interval data is presented.
     */
    int stepMax_;

    //! If true then the depenVals_ vector has a last dependant value that is valid for independent values greater than the 
    //! last interval ending value.
    /*!
     *  Defaults to 2, as this is the normal case if no intervals are supplied.
     *                                             0:  Treat this as a fatal error and throw an error condition
     *                                             1:  Assume the same treatment as the last interval in the boundary condition.
     *                                                 No extra information is needed. functions of the independent variable are
     *                                                 continuous as the boundary condtion goes above the upper limit.
     *                                             2:  Another interval is assumed above the last interval with a double
     *                                                 value for the dependent value input in the parameter list.
     */
    int hasExtendedDependentValue_;

    //! Treatment of what to do when independent value is below the lower limit
    /*!
     *                                             0:  Treat this as a fatal error and throw an error condition
     *                                             1:  Assume the same treatment as the zeroeth interval in the boundary condition.
     *                                                 No extra information is needed. functions of the independent variable are
     *                                                 continuous as the boundary condtion goes below the lower limit.
     *                                             2:  Another interval is assumed below the zeroeth interval with a double
     *                                                 value for the dependent value input in the parameter list.
     */ 
    int ifuncLowerLim_;

    //! Value of the dependent value to use when below the lower limit
    double depValLowerLim_;

private:
    //! Error routine
    /*!
     *  Throw an exception if an unimplemented method of this class is invoked.
     *
     *  @param[in]           msg                 Descriptive message string to add to the error report
     *
     *  @return                                  returns a double, though we will never get there
     */
    double err(std::string msg) const;
};
//==================================================================================================================================
//! BCconstant is a simple class that specifies a constant value Dirichlet boundary condition
/*!
 *  No intervals are specified. Instead, the dependent variable is set to a constant
 */
class BCconstant: public BoundaryCondition
{
public:

    //! Default constructor for BCconstant
    /*!
     *  @param[in]           value               Value of the constant to set the dependent value to. Defaults to 0.0
     *  @param[in]           titleName           Name of the Boundary Condition. Defaults to "BCtitle".
     *  @param[in]           indepUnits          String containing units of independent variable. Defaults to "unknownUnits".
     *  @param[in]           depenUnits          String containing units of dependent variable. Defaults to "unknownUnits".
     */
    BCconstant(double value = 0.0, std::string titleName = "BCtitle", std::string indepUnits = "unknownUnits",
               std::string depenUnits = "unknownUnits");

    //! Virtual destructor
    virtual ~BCconstant();

    //! Return the dependent variable value given the independent variable argument
    /*!
     *   @param[in]          indVar              Independentvariable
     *   @param[in]          interval            If greater than zero, then checking is done on the interval specified
     *                                           Also ties, i.e. numbers on the boundary go to the interval value.
     *
     *  @return                                  Returns the value of the dependent variable
     */
    virtual double value(double indVar, int interval = -1) const override;

    //! Return the next value for the independent variable at which the nature of the boundary condition changes.
    /*!
     *  This is designed to guide grid generation and time stepping. For this constant BC subclass, 
     *  this provides no real information.
     *
     *  @return                                  Returns a max value, because there are no steps
     */
    virtual double nextStep() const override;

protected:

    //! The fixed dependent variable
    double dependentVal_;
};

//==================================================================================================================================
//! This subclass is designed to handle a table of dependent variable boundary conditions that are constant until the indicated 
//! value of the independent variable, at which point there is a step change to a new value of the dependent variable
/*!
 *  For example, if the (ind,dep) value pairs (0.5,0.0), (1.0,0.5) and (2., 1.) are given, then value( 0.25 ) will 
 *  return 0.0, value( 0.75 ) will return 0.5 and  value( 1.25 ) will return 1.0.
 *
 *  Usually, the size of the dep
 *  If the depValue vector size is 
 */
class BCsteptable: public BoundaryCondition
{
public:
    //! Default constructor for BCsteptable
    /*!
     *  @param[in]           indValue            Reference to a vector of independent values representing the boundaries of the
     *                                           intervals.
     *                                             Length: numIntervals+1
     *  @param[in]           depValue            Reference to a vector of dependent values representing the value of the dependent
     *                                           variable within the intervals. For a size of  (numIntervals+1), the last number
     *                                           represents the value of the dependent variable past the last interval value.
     *                                            Length: numIntervals or (numIntervals+1)
     *  @param[in]           compareVals_        Reference to a vector of data to compare to. (currently unused).
     *  @param[in]           titleName           Name of the Boundary Condition. Defaults to "BCsteptable".
     *  @param[in]           indepUnits          String containing units of independent variable. Defaults to "unknownUnits".
     *  @param[in]           depenUnits          String containing units of dependent variable. Defaults to "unknownUnits".
     */
    BCsteptable(const ZZCantera::vector_fp& indValue, const ZZCantera::vector_fp& depValue, const ZZCantera::vector_fp& compareVals_, 
                const std::string& titleName = "BCsteptable", const std::string& indepUnits = "unknownUnits", 
                const std::string& depenUnits = "unknownUnits");

    //! Construct the boundary condition from a filename
    /*!
     *  @param[in]           filename            The filename contains a complete description of the BCsteptable boundary condition
     *                                           input. The file is in XML format.
     */
    BCsteptable(std::string filename);

    //! Construct from XMLnode
    /*!
     *  @param[in]           node                Reference to an XML node where the information for the boundary condition is storred.
     */
    BCsteptable(ZZCantera::XML_Node& node);

    //! Destructor
    virtual ~BCsteptable();

    //! Fill independent and dependent values from ZZCantera::XML_Node object named BoundaryCondition
    /*!
     *  XML_Node object contains all of the information for the boundary condition. The name of the node is "BoundaryCondition". 
     *
     *  (virtual from BoundaryCondition)
     *
     *  @param[in]           bcNode                Reference to an XML node where the information for the boundary condition is storred
     */
    virtual void useXML(ZZCantera::XML_Node& bcNode) override;

    //! Return the dependent variable value given the independent variable argument
    /*!
     *  @param[in]           indVar              Independent variable
     *  @param[in]           interval            If greater than zero, then checking is done on the interval specified
     *                                           Also ties, i.e. numbers on the boundary go to the interval value.
     * 
     *  @return                                  Returns the value of the boundary condition
     */
    virtual double value(double indVar, int interval = -1) const override;

    //! Return the next value for the independent variable at which the nature of the boundary condition changes.
    /*!
     *  This is designed to guide grid generation and time stepping
     *
     *  @return                                  returns the value of the independent variable for the next big change
     */
    virtual double nextStep() const override;

protected:

    //! Units string for a variable used for comparison purposes
    std::string compareUnits_;

};
//==================================================================================================================================
//! This subclass is designed to handle a table of dependent variable boundary conditions that are
//! to be linearly interpolated between the values given.
/*!
 *  This implies that the function is continuous, and that the derivatives are piecewise continuous.
 *
 *  For example, if the value pairs (0,0) and (1,1) are given, the value( 0.5 ) will  return 0.5.
 */
class BClineartable: public BoundaryCondition
{
public:
    //! Default constructor for BClineartable
    /*!
     *  Function is assumed to be continuous, and its derivative is piecewise continuous.
     *
     *  @param[in]           indValue            Reference to a vector of independent values representing the boundaries of the
     *                                           intervals.
     *                                              Length: numIntervals+1
     *  @param[in]           depValue            Reference to a vector of dependent values representing the value of the dependent
     *                                           variables at the indValue[] points given in the previous parameter.
     *                                           represents the value of the dependent variable past the last interval value.
     *                                              Length: numIntervals+1
     *  @param[in]           compareVals_        Reference to a vector of data to compare to. (currently unused).
     *  @param[in]           titleName           Name of the Boundary Condition. Defaults to "BCsteptable".
     *  @param[in]           indepUnits          String containing units of independent variable. Defaults to "unknownUnits".
     *  @param[in]           depenUnits          String containing units of dependent variable. Defaults to "unknownUnits".
     */
    BClineartable(ZZCantera::vector_fp indValue, ZZCantera::vector_fp depValue, ZZCantera::vector_fp compareVals_, std::string titleName = "BClineartable",
                  std::string indepUnits = "unknownUnits", std::string depenUnits = "unknownUnits");

    //! Construct the BClineartable from information contained in an xml file
    /*!
     *  All of the parameters are read from an XML formatted file.
     *
     *  @param[in]           filename           Name of the file
     */
    BClineartable(std::string filename);

    //! Construct the boundary condition from an XMLnode object
    /*!
     *  All of the parameters are read from an XML formatted file.
     *
     *  @param[in]           node           Reference to an XML_Node object containing the boundary condition description
     */
    BClineartable(ZZCantera::XML_Node& node);

    //! Virtual destructor
    virtual ~BClineartable();

    //! Fill independent and dependent values from ZZCantera::XML_Node object named BoundaryCondition
    /*!
     *  XML_Node object contains all of the information for the boundary condition. The name of the node is "BoundaryCondition". 
     *
     *  (virtual from BoundaryCondition)
     *
     *  @param[in]           bcNode                Reference to an XML node where the information for the boundary condition is storred
     */
    virtual void useXML(ZZCantera::XML_Node& bcNode) override;

    //! Return the dependent variable value given
    //! the independent variable argument
    /*!
     *  @param[in]           indVar              Independent variable
     *  @param[in]           interval            If greater than zero, then checking is done on the interval specified.
     *                                           Also ties, i.e. numbers on the boundary go to the interval value.
     *
     *  @return                                  Returns the value of the dependent variable for the boundary condition
     */
    virtual double value(double indVar, int interval = -1) const override;

    //! return the next value for the independent variable at which the nature of the boundary condition changes.
    /*!
     *  This is designed to guide grid generation and time stepping
     *
     *  @return                                  Returns the next value of the independent variable for which the bc changes
     */
    virtual double nextStep() const override;

protected:

    //! units string for a variable used for comparison purposes
    std::string compareUnits_;
};

//==================================================================================================================================
//! This boundary condition imposes a sin() function on the value of the dependent variable that is function of the independent
//! variable
/*!
 *    The value of this boundary condition is given by the following functional form.
 *
 *      f(t) = baseDepValue + oscAmplitude * cos(2 Pi t frequency + phaseAngle) 
 */
class BCsinusoidal: public BoundaryCondition
{
public:

    //! Default constructor for BCsinusoidal
    /*!
     *  Function is assumed to be continuous, and its derivative is continuous as well.
     *
     *  @param[in]           baseDepValue        Base value of the dependent variable that is added onto the sine function.
     *                                           This is also the value of the function at T = 0;
     *
     *  @param[in]           oscAmplitude        Value of the amplitude of the dependent variable. This multiplies the sine() function.
     *
     *  @param[in]           frequency           Value of the frequency of the sine function. The period of the sine function is 
     *                                           equal to 1/frequency. 
     *
     *  @param[in]           phaseAngle          Value of the phase angle. Units: radians
     *
     *  @param[in]           titleName           Name of the Boundary Condition. Defaults to "BCsteptable".
     *  @param[in]           indepUnits          String containing units of independent variable. Defaults to "unknownUnits".
     *  @param[in]           depenUnits          String containing units of dependent variable. Defaults to "unknownUnits".
     */
    BCsinusoidal(double baseDepValue, double oscAmplitude, double frequency, double phaseAngle, std::string titleName = "BCsinusoidal",
                 std::string indepUnits = "unknownUnits", std::string depenUnits = "unknownUnits");

    //! Constructor from an XML file, named filename
    /*!
     *  @param[in]           filename            Filename to look up the boundary condition from
     */
    BCsinusoidal(std::string filename);

    //! Construct from XMLnode
    /*!
     *  @param[in]           node                Reference to an XML_Node containing all of the boundary condition information
     */
    BCsinusoidal(ZZCantera::XML_Node& node);

    //! Virtual Destructor
    virtual ~BCsinusoidal();

    //! Fill independent and dependent values from ZZCantera::XML_Node object named BoundaryCondition
    /*!
     *  XML_Node object contains all of the information for the boundary condition. The name of the node is "BoundaryCondition". 
     *
     *  (virtual from BoundaryCondition)
     *
     *  @param[in]           bcNode                Reference to an XML node where the information for the boundary condition is storred
     */
    void useXML(ZZCantera::XML_Node& bcNode) override;

    //! Return the dependent variable value given the independent variable argument
    /*!
     *  @param[in]           indVar              Independent variable
     *  @param[in]           interval            If greater than zero, then checking is done on the interval specified
     *                                           Also ties, i.e. numbers on the boundary go to the interval value.
     *
     *  @return                                  returns the value of the boundary condition
     */
    virtual double value(double indVar, int interval = -1) const override;

    //! return the next value for the independent variable at which the nature of the boundary condition changes.
    /*!
     *  This is designed to guide grid generation and time stepping
     *
     *  @return                                  Returns the next deltaT;
     */
    virtual double nextStep() const override;

    //! Specify the number of steps per period
    /*!
     *  @param[in]           stepsPerPeriod      Input the number of steps per period of the oscillation.
     *                                           This sets the maximum time step.
     */
    void setStepsPerPeriod(double stepsPerPeriod);

    //! specify the number of periods to be run
    /*!
     *  @param[in]           periodsPerRun       Number of periods to run the simulation. This sets the time of the run.
     */
    void setPeriodsPerRun(double periodsPerRun);

    //! specify the frequency of the oscillation
    /*!
     *  The frequency is the inverse of the period of the oscillation
     *
     *  @param[in]           frequency           Frequency of the oscillation
     *                                             Units: 1 / indValUnit
     */
    void setFrequency(double frequency);

    //! Specify the phase angle of the sinusoidal function
    /*!
     *  @param[in]           phaseAngle          Phase angle
     *                                             Units: radians
     */
    void setPhaseAngle(double phaseAngle);

protected:

    //! Base value about which oscillations occur
    double baseDepValue_;

    //! Oscillation amplitude
    double oscAmplitude_;

    //! Frequency of oscillation
    double frequency_;

    //! Phase Anglele
    double phaseAngle_;

    //! Define a deltaT independent variable step size based on the frequency
    double deltaT_;

    //! Specify the number of steps per period
    double stepsPerPeriod_;

    //! Specify the number of periods to be run in the simulation
    double periodsPerRun_;

    //! Calculate which interval of the independent variable we are in.
    /*!
     * If the independent variable argument exceeds the current range, then increment the step counter and
     * check that this has not gone out of bounds.
     *
     *  @param[in]           indVar              Value of the independent variable.
     * @param[in]            interval            If positive, then ties goes to the value of the interval.
     *                                           So, if indVar is at a boundary, then the interval chose is
     *                                           equal to the value of the interval variable.
     *
     *  @return                                  Returns the interval
     */
    virtual int findStep(double indVar, int interval = -1) const override;

};
//==================================================================================================================================
} 
//----------------------------------------------------------------------------------------------------------------------------------
#endif // M1D_BOUNDARYCONDITION
