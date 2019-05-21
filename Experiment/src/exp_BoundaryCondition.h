/**
 * file exp_BoundaryCondition.h
 * 
 * Header for class BoundaryCondition and subclasses.
 */

/*  $Author: hkmoffa $
 *  $Revision: 508 $
 *  $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
 *
 */
// Copyright 2010 Sandia National Laboratories

#ifndef EXP_BOUNDARYCONDITION
#define EXP_BOUNDARYCONDITION

#include "zuzax/base/ct_defs.h" 
#include "zuzax/base/ctexceptions.h"

#include "zuzax/base/xml.h"
#include "zuzax/base/ctml.h"
using namespace ztml;


namespace Zuzax {

/** 
 * The BoundaryCondition class provides a single scalar dependent 
 * variable as a function a single independent variable.  This 
 * is suitable for either spatial or temporal boundary conditions.
 *
 * This is an abstract base class.  Subclasses provide specific
 * relationships between the dependent variable and independent variable.
 */
class BoundaryCondition {

 public:

  //!  Constructor.
  BoundaryCondition(  ) :
    title_(""),
    lowerLim_( -1.0 * BigNumber ), 
    upperLim_( BigNumber ), 
    indepUnits_( "" ),
    depenUnits_( "" ),
    step_(0),
    stepMax_(0)
    { }

  //! Destructor.
  virtual ~BoundaryCondition(  ) { ; }

  //!  Copy Constructor.
  BoundaryCondition( const BoundaryCondition &right ) :
    title_(""),
    lowerLim_( -1.0 * BigNumber ), 
    upperLim_( BigNumber ), 
    indepUnits_( "" ),
    depenUnits_( "" ),
    step_(0),
    stepMax_(0)
    { *this = right; }

  //! Assignment operator
  BoundaryCondition&  operator=(const BoundaryCondition& right) {
    if (&right != this) { return *this; }
    title_           = right.title_;
    lowerLim_        = right.lowerLim_;
    upperLim_        = right.upperLim_;
    indepUnits_      = right.indepUnits_;
    depenUnits_      = right.depenUnits_;
    step_            = right.step_;
    stepMax_         = right.stepMax_;
    return *this;
  }

  //! return the dependent variable value given 
  //! the independent variable argument 
  virtual double value( double indVar )
  { return err("BoundaryCondition::value"); }

  //! return the next value for the independent variable at 
  //! which the nature of the boundary condition changes. 
  /** 
   * This is designed to guide grid generation and time stepping
   */
  virtual double nextStep( ) 
  { return err("BoundaryCondition::nextStep"); }

  //! Write out the profile in tabular format.
  virtual void writeProfile( ){};

  //!lower limit of dependent variable for which BC applies
  double lowerLimit() { return lowerLim_; }

  //!upper limit of dependent variable for which BC applies
  double upperLimit() { return upperLim_; }

  //! the string defining the independent variable units 
  std::string indepUnits() { return indepUnits_; }

  //! the string defining the dependent variable units 
  std::string depenUnits() { return depenUnits_; }

  //! return the title or name of boundary condition
  std::string title( ) { return title_ ; }

  //! set title or name of boundary condition
  void setTitle( std::string name ) { title_ = name; }

  //!set lower limit of dependent variable for which BC applies
  void setLowerLimit( double indVal ) { lowerLim_ = indVal;}

  //!set upper limit of dependent variable for which BC applies
  void setUpperLimit( double indVal ) { upperLim_ = indVal;}

  //! set the string defining the independent variable units 
  void setIndepUnits( std::string unitString ) { indepUnits_ = unitString; }

  //! set the string defining the dependent variable units 
  void setDepenUnits( std::string unitString ) { depenUnits_ = unitString; }

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

  //! check to see which step we are at.  
  /** 
   * If the independent variable argument exceeds the 
   * current range, then increment the step counter and
   * check that this has not gone out of bounds.
   */
  virtual void findStep( double indVar ) { 
     err("BoundaryCondition::findStep"); }
  
  //! Error routine
  /*!
   * Throw an exception if an unimplemented method of this class is
   * invoked. 
   *
   *  @param msg  Descriptive message string to add to the error report
   *
   *  @return  returns a double, though we will never get there
   */
  double err(std::string msg) const{
    throw ZuzaxError("BoundaryCondition Base Class\n",
		       "**** Method "+ msg +" not implemented\n");
    return 0.0;
  }


 private:


};

////////////////////////////////////////////////////////////
// class BCconstant
////////////////////////////////////////////////////////////

class BCconstant : public BoundaryCondition {

 public:

  BCconstant( double value = 0.0, 
	      std::string titleName = "BCtitle" , 
	      std::string indepUnits = "unknownUnits" , 
	      std::string depenUnits = "unknownUnits" ) :
    BoundaryCondition(),
    dependentVal_( value )
  {   
    setTitle( titleName );
    setIndepUnits( indepUnits );
    setDepenUnits( depenUnits );
  }

  ~BCconstant( ){ }

  //! return the dependent variable value given 
  //! the independent variable argument 
  virtual double value( double indVar )
  { return dependentVal_; }

  //! return the next value for the independent variable at 
  //! which the nature of the boundary condition changes. 
  /** 
   * This is designed to guide grid generation and time stepping.
   * For this constant BC subclass, this provides no real information.
   */
  virtual double nextStep( ) 
  { return upperLim_; }

 protected:
  
  //The fixed dependent variable
  double  dependentVal_;

};

/**
 * This subclass is designed to handle a table of 
 * dependent variable boundary conditions that are 
 * constant until the indicated value of the 
 * independent varaible, at which piont there is a
 * step change to a new value.  For example, if 
 * the (ind,dep) value pairs (0.5,0.0), (1.0,0.5)
 * and (2., 1.) are given, then value( 0.25 ) will 
 * return 0.0, value( 0.75 ) will return 0.5 and
 * value( 1.25 ) will return 1.0.
 */
class BCsteptable : public BoundaryCondition {

 public:

  BCsteptable( vector_fp indValue, vector_fp depValue, 
	       std::string titleName = "BCsteptable" , 
	       std::string indepUnits = "unknownUnits", 
	       std::string depenUnits = "unknownUnits" ); 


  //! construct from filename
  BCsteptable( std::string filename );

  //! construct from XMLnode
  BCsteptable( XML_Node& node );

  //! destructor
  virtual ~BCsteptable( ){ ; }

  //! fill independent and dependent values from XML_Node
  void useXML( XML_Node& node ) ;


  //! return the dependent variable value given 
  //! the independent variable argument 
  virtual double value( double indVar );

  //! return the next value for the independent variable at 
  //! which the nature of the boundary condition changes. 
  /** 
   * This is designed to guide grid generation and time stepping
   */
  virtual double nextStep( ) ;

  virtual void writeProfile( );

 protected:

  //!vector of indepedent variable values at which
  //! the dependent variable may change value
  vector_fp indepVals_;

  //!vector of depedent variable values appropriate 
  //! for time/space after the corresponding indepVals_ 
  vector_fp depenVals_;

  //! check to see which step we are at.  
  /** 
   * If the independent variable argument exceeds the 
   * current range, then increment the step counter and
   * check that this has not gone out of bounds.
   */
  virtual void findStep( double indVar ) ;


};

/**
 * This subclass is designed to handle a table of 
 * dependent variable boundary conditions that are 
 * to be linearly interpolated between the values
 * given.  For example, if the value pairs (0,0) 
 * and (1,1) are given, the value( 0.5 ) will 
 * return 0.5.
 */
class BClineartable : public BoundaryCondition {

 public:

  BClineartable( vector_fp indValue, vector_fp depValue, 
	       std::string titleName = "BClineartable" , 
	       std::string indepUnits = "unknownUnits", 
	       std::string depenUnits = "unknownUnits" ); 


  //! construct from filename
  BClineartable( std::string filename );

  //! construct from XMLnode
  BClineartable( XML_Node& node );

  //! destructor
  virtual ~BClineartable( ){ ; }

  //! fill independent and dependent values from XML_Node
  void useXML( XML_Node& node ) ;



  //! return the dependent variable value given 
  //! the independent variable argument 
  virtual double value( double indVar );


  //! return the next value for the independent variable at 
  //! which the nature of the boundary condition changes. 
  /** 
   * This is designed to guide grid generation and time stepping
   */
  virtual double nextStep( ) ;

 protected: 

  //!vector of indepedent variable values at which
  //! the dependent variable may change value
  vector_fp indepVals_;

  //!vector of depedent variable values appropriate 
  //! for time/space after the corresponding indepVals_ 
  vector_fp depenVals_;

  //! check to see which step we are at.  
  /** 
   * If the independent variable argument exceeds the 
   * current range, then increment the step counter and
   * check that this has not gone out of bounds.
   */
  virtual void findStep( double indVar ) ;

};


} 

#endif // EXP_BOUNDARYCONDITION
