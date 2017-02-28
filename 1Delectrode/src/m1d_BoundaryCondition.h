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
     * If the independent variable argument exceeds the current range, then increment the step counter and
     * check that this has not gone out of bounds.
     *
     *  @param[in]           indVar              Independent variable, usually the time.
     *  @param[in]           interval            If greater than zero, then checking is done on the interval specified.
     *                                           Also ties wrt to the independent variable are satisfied by goint to the interval
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

    //! Lower limit of dependent variable for which BC applies
    /*!
     *  @return                                  Returns the  Lower limit of dependent variable 
     */
    double lowerLimit() const;

    //! Upper limit of dependent variable for which BC applies
    /*!
     *  @return                                  Returns the  Upper limit of dependent variable 
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
     *  @return                                  Set the title string
     */
    void setTitle(const std::string &name);

    //! Set lower limit of independent variable for which BC applies
    /*!
     *  @param[in]           indVal              Set the lower limit for the independent variable
     */
    void setLowerLimit(double indVal);

    //! Set upper limit of independent variable for which BC applies
    /*!
     *  @param[in]           indVal              Set the lower limit for the independent variable
     */
    void setUpperLimit(double indVal);

    //! set the string defining the independent variable units
    /*!
     *  @param[in[           unitString          Set the units string for the independent variable
     */
    void setIndepUnits(const std::string& unitString);

    //! set the string defining the dependent variable units
    /*!
     *  @param[in[           unitString          Set the units string for the dependent variable
     */
    void setDepenUnits(const std::string& unitString);

    // --------------------------------------------------------- D A T A ----------------------------------------------
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

    //! Current step interval in a sequence of values
    /*!
     *  Starts at a value of zero
     */
    mutable int step_;

    //! Maximum value of the interval steps.
    /*!
     *  If step_ exceeds the number of steps that were input, then an out of bounds error should be generated
     *  Starts at a value of zero.
     */
    int stepMax_;

private:
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
//! value of the independent variable, at which point there is a step change to a new value. 
/*!
 *  For example, if the (ind,dep) value pairs (0.5,0.0), (1.0,0.5) and (2., 1.) are given, then value( 0.25 ) will 
 *  return 0.0, value( 0.75 ) will return 0.5 and
 *  value( 1.25 ) will return 1.0.
 */
class BCsteptable: public BoundaryCondition
{
public:

    //! Default constructor for BCsteptable
    /*!
     *  @param[in]           indValue            Reference to a vector of independent values representing the intervals
     *  @param[in]           depValue            Reference to a vector of dependent values representing the value of the dependent
     *                                           variable within the interval.
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
    BCsteptable(ZZCantera::XML_Node& node);

    //! Destructor
    virtual ~BCsteptable();

    //! Fill independent and dependent values from ZZCantera::XML_Node
    void useXML(ZZCantera::XML_Node& node);

    //! Return the dependent variable value given
    //! the independent variable argument
    /*!
     *   @param indVar   Independent variable
     *   @param interval If greater than zero, then checking is done on the interval specified
     *                   Also ties, i.e. numbers on the boundary go to the interval value.
     */
    virtual double value(double indVar, int interval = -1) const override;

    //! Return the next value for the independent variable at
    //! which the nature of the boundary condition changes.
    /**
     * This is designed to guide grid generation and time stepping
     */
    virtual double nextStep() const override;

    virtual void writeProfile() const override;

protected:

    //! vector of independent variable values at which
    //! the dependent variable may change value
    ZZCantera::vector_fp indepVals_;

    //! vector of dependent variable values appropriate
    //! for time/space after the corresponding indepVals_
    ZZCantera::vector_fp depenVals_;

    //! vector of variable values for comparison purposes.
    //! For example, if current is input, these might be measured voltages
    ZZCantera::vector_fp compareVals_;

    //! units string for a variable used for comparison purposes
    std::string compareUnits_;

    //! Calculate which interval of the independent variable we are in.
    /*!
     * If the independent variable argument exceeds the
     * current range, then increment the step counter and
     * check that this has not gone out of bounds.
     */
    virtual int findStep(double indVar, int interval = -1) const override;

};
//==================================================================================================================================
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

    BClineartable(ZZCantera::vector_fp indValue, ZZCantera::vector_fp depValue, ZZCantera::vector_fp compareVals_, std::string titleName = "BClineartable",
                  std::string indepUnits = "unknownUnits", std::string depenUnits = "unknownUnits");

    //! construct from filename
    BClineartable(std::string filename);

    //! construct from XMLnode
    BClineartable(ZZCantera::XML_Node& node);

    //! destructor
    virtual ~BClineartable();

    //! fill independent and dependent values from ZZCantera::XML_Node
    void useXML(ZZCantera::XML_Node& node);

    //! Return the dependent variable value given
    //! the independent variable argument
    /*!
     *   @param indVar   Independent variable
     *   @param interval If greater than zero, then checking is done on the interval specified.
     *                   Also ties, i.e. numbers on the boundary go to the interval value.
     */
    virtual double value(double indVar, int interval = -1) const override;

    //! return the next value for the independent variable at
    //! which the nature of the boundary condition changes.
    /**
     * This is designed to guide grid generation and time stepping
     */
    virtual double nextStep() const override;

protected:

    //! Vector of indepedent variable values at which the dependent variable may change value,
    //! have a discontinuity, or may change functional form.
    /*!
     *  Each interval, i, is defined as existing between indepVals_[i] and indepVals_[i+1].
     *  Length: numIntervals + 1
     */
    ZZCantera::vector_fp indepVals_;

    //! Vector of dependent variable values appropriate
    //! for time/space after the corresponding indepVals_
    /*!
     *  dependVals_[i] refers to the value that is appropriate for the ith interval
     *  dependVals_[numIntervals] refers to the value that is appropriate for the independent
     *    values beyond the last intervale.
     *  Length:  numIntervals + 1
     */
    ZZCantera::vector_fp depenVals_;

    //! Vector of variable values for comparison purposes
    //! For example, if current is input, these might be measured voltages
    ZZCantera::vector_fp compareVals_;

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
    virtual int findStep(double indVar, int interval = -1) const override;

};
//==================================================================================================================================
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
    BCsinusoidal(ZZCantera::XML_Node& node);

    //! Destructor
    virtual ~BCsinusoidal();

    //! fill independent and dependent values from ZZCantera::XML_Node
    void useXML(ZZCantera::XML_Node& node);

    //! Return the dependent variable value given
    //! the independent variable argument
    /*!
     *   @param indVar  Independent variable
     *   @param interval If greater than zero, then checking is done on the interval specified
     *                   Also ties, i.e. numbers on the boundary go to the interval value.
     */
    virtual double value(double indVar, int interval = -1) const override;

    //! return the next value for the independent variable at
    //! which the nature of the boundary condition changes.
    /**
     * This is designed to guide grid generation and time stepping
     */
    virtual double nextStep() const override;

    virtual void writeProfile() const override;

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
    virtual int findStep(double indVar, int interval = -1) const override;

};
//==================================================================================================================================
} 
//----------------------------------------------------------------------------------------------------------------------------------
#endif // M1D_BOUNDARYCONDITION
