/**
 * file exp_BoundaryCondition.cpp
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
 *  $Revision: 508 $
 *  $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
 *
 */
// Copyright 2010 Sandia National Laboratories

#ifndef CANTERA_APP
#define CANTERA_APP
#endif

#include "cantera/base/ct_defs.h" 
#include "cantera/base/ctexceptions.h"
#include "exp_BoundaryCondition.h" 

namespace Cantera {

  ////////////////////////////////////////////////////////////
  // class BCsteptable 
  ////////////////////////////////////////////////////////////

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

//! constructor taking vectors of floats
  BCsteptable::BCsteptable( vector_fp indValue, vector_fp depValue, 
			    std::string titleName , 
			    std::string indepUnits,
			    std::string depenUnits ) :
    BoundaryCondition(),
    indepVals_( indValue ),
    depenVals_( depValue )
  {
    setTitle( titleName );
    setIndepUnits( indepUnits );
    setDepenUnits( depenUnits );

    if ( indepVals_.size() != depenVals_.size() ) { 
      throw CanteraError("BCsteptable constructor\n",
			 "**** indepVals_ and depenVals_ unequal size\n");
      
    }
  }
   

  //! construct from filename
  BCsteptable::BCsteptable( std::string filename ) :
    BoundaryCondition() 
  {
    //convert input file name into XML node
    XML_Node* baseNode = get_XML_File( filename );
    //find the "boundaryCondition" node
    std::string targetName = "boundaryCondition";
    XML_Node* bcNode = get_XML_NameID( targetName, "", baseNode);
    //parse this node into class data structures
    useXML( *bcNode );
    close_XML_File( filename );
  }

  //! construct from XMLnode
  BCsteptable::BCsteptable( XML_Node& baseNode ) :
    BoundaryCondition() 
  {
    useXML( baseNode );
  }

  //! fill independent and dependent values from XML_Node
  void BCsteptable::useXML( XML_Node& bcNode ) {
    
    if ( !bcNode.hasChild("independentVar") ) 
      throw CanteraError("BCsteptable::useXML()",
			 "no independentVar XML node.");
    
    if ( !bcNode.hasChild("dependentVar") ) 
      throw CanteraError("BCsteptable::useXML()",
			 "no dependentVar XML node.");
    
    bool convert = true;

    //get independentVar
    XML_Node& indVarNode = bcNode.child("independentVar");
    getFloatArray( indVarNode, indepVals_, convert, indepUnits_ );

    //get dependentVar
    XML_Node& depVarNode = bcNode.child("dependentVar");
    getFloatArray( depVarNode, depenVals_, convert, depenUnits_ );


    //err("BCsteptable::useXML");




  }


  //! return the dependent variable value given 
  //! the independent variable argument 
  double BCsteptable::value( double indVar ) {
    findStep( indVar );
    return depenVals_[ step_ + 1 ];    
  }

  //! return the next value for the independent variable at 
  //! which the nature of the boundary condition changes. 
  /** 
   * This is designed to guide grid generation and time stepping
   */
  double BCsteptable::nextStep( ) {
    return indepVals_[ step_ + 1 ];
  }

  //! check to see which step we are at.  
  /** 
   * If the independent variable argument exceeds the 
   * current range, then increment the step counter and
   * check that this has not gone out of bounds.
   */
  void BCsteptable::findStep( double indVar ) { 

    if ( indVar > indepVals_[ step_ ] ) 
      { 
	step_++;
	if ( step_ > stepMax_ ) 
	  throw CanteraError("exp_BoundaryCondition::findStep",
			     "Out of bounds error with step_ > stepMax_");
	if ( indVar > indepVals_[ step_ ] ) 
	  throw CanteraError("exp_BoundaryCondition::findStep",
			     "BoundaryCondition independent variable step was greater than the BoundaryCondition resolution.");
      }
  }


  //! Write out the profile in tabular format.
  void BCsteptable::writeProfile( ) {
    
  }



  ////////////////////////////////////////////////////////////
  // class BClineartable
  ////////////////////////////////////////////////////////////
/**
 * This subclass is designed to handle a table of 
 * dependent variable boundary conditions that are 
 * to be linearly interpolated between the values
 * given.  For example, if the value pairs (0,0) 
 * and (1,1) are given, the value( 0.5 ) will 
 * return 0.5.
 */


//! constructor taking vectors of floats
  BClineartable::BClineartable( vector_fp indValue, vector_fp depValue, 
				std::string titleName, 
				std::string indepUnits, 
				std::string depenUnits ) :
    BoundaryCondition(),
    indepVals_( indValue ),
    depenVals_( depValue )
  {
    setTitle( titleName );
    setIndepUnits( indepUnits );
    setDepenUnits( depenUnits );

    if ( indepVals_.size() != depenVals_.size() ) { 
      throw CanteraError("BClineartable constructor\n",
			 "**** indepVals_ and depenVals_ unequal size\n");
      
    }
    stepMax_ = indepVals_.size() ;
  }
   

  //! construct from filename
  BClineartable::BClineartable( std::string filename ) :
    BoundaryCondition() 
  {
    //convert input file name into XML node
    XML_Node *baseNode;
    baseNode = get_XML_File( filename );
    //std::istream& infile( filename );
    //Cantera::XML_Node* baseNode;
    //baseNode.build( infile );
    useXML( *baseNode );
    close_XML_File( filename );
  }

  //! construct from XMLnode
  BClineartable::BClineartable( XML_Node& baseNode ) :
    BoundaryCondition() 
  {
    useXML( baseNode );
  }

  //! fill independent and dependent values from XML_Node
  void BClineartable::useXML( XML_Node& node ) {



    err("BClineartable::useXML");


    if ( indepVals_.size() != depenVals_.size() ) { 
      throw CanteraError("BClineartable::useXML\n",
			 "**** indepVals_ and depenVals_ unequal size\n");
      
    }
    stepMax_ = indepVals_.size() ;
  }


  //! return the dependent variable value given 
  //! the independent variable argument 
  double BClineartable::value( double indVar ) {
    findStep( indVar );
    return depenVals_[ step_ + 1 ];    
  }

  //! return the next value for the independent variable at 
  //! which the nature of the boundary condition changes. 
  /** 
   * This is designed to guide grid generation and time stepping
   */
  double BClineartable::nextStep( ) {
    return indepVals_[ step_ + 1 ];
  }

  //! check to see which step we are at.  
  /** 
   * If the independent variable argument exceeds the 
   * current range, then increment the step counter and
   * check that this has not gone out of bounds.
   */
  void BClineartable::findStep( double indVar ) { 

    if ( indVar > indepVals_[ step_ + 1 ] ) 
      { 
	step_++;
	if ( step_ > stepMax_ ) 
	  throw CanteraError("exp_BoundaryCondition::findStep",
			     "Out of bounds error with step_ > stepMax_");
	if ( indVar > indepVals_[ step_ + 1 ] ) 
	  throw CanteraError("exp_BoundaryCondition::findStep",
			     "BoundaryCondition independent variable step was greater than the BoundaryCondition resolution.");
      }
  }







}//namespace Cantera

