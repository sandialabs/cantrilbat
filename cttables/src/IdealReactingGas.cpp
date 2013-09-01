/**
 * @file IdealReactingGas.cpp
 *
 * $Author: hkmoffa $
 * $Revision: 497 $
 * $Date: 2013-01-07 14:17:04 -0700 (Mon, 07 Jan 2013) $
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "IdealReactingGas.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"

#include <cstdio>
#include <fstream>

using namespace std;
namespace Cantera {

    /**
     * Constructor #1 for the IdealReactingGas object.
     * 
     * infile = input file location
     * id = ?
     */
    IdealReactingGas::
    IdealReactingGas(string infile, string id) :
	m_ok(false),
	m_r(0) ,
	m_r_owned(true)
    {
	/*
	 * Call the global function in inputCTML.cpp that finds
	 * the input file. It looks in the current directory,
	 * plus a few other previously defined directories.
	 */
	string path = findInputFile(infile);
	ifstream fin(path.c_str());
	if (!fin) {
	  throw CanteraError("IdealReactingGas","could not open "
			     + path + " for reading.");
	}
	
	m_r = new XML_Node("-");
	m_r->build(fin);
	/*
	 * Call the function in importCTML.cpp that builds a model
	 * from an XML tree description of the mechanism.
	 */
	m_ok = buildSolutionFromXML(*m_r, id, "phase", this, this);
	if (!m_ok) {
	  throw CanteraError("IdealReactingGas",
			     "buildSolutionFromXML returned false");
	}
    }

    /**
     * Constructor #2 for the IdealReactingGas object.
     * 
     * root = Previously built xml tree structure (not owned
     *        by this object 
     * id = ?
     */
    IdealReactingGas::
    IdealReactingGas(XML_Node& root, string id) :
	m_ok(false),
	m_r(0),
	m_r_owned(false)
    {
	m_ok = buildSolutionFromXML(root, id, "phase", this, this);
	if (!m_ok) {
	  throw CanteraError("IdealReactingGas",
			     "buildSolutionFromXML returned false");
	}
    }
 
    /**
     * Destructor for the IdealReactingGas object.
     * 
     * We must decide whether this object owns its own xml tree
     * structure.
     */       
    IdealReactingGas::~IdealReactingGas() { 
	if (m_r_owned) {
	  delete m_r;
	}
    }
   
    /** 
     *  This ostream function describes how to extend cout output
     *  functions to this object. The way this is done is to
     *  call the Cantera's report function, which takes a ThermoPhase
     *  object as its arugment.
     *  This is a "friend" function to the class IdealReactingGas.
     *  Both this function and report are in the Cantera namespace.
     *
     *  Note -> The output doesn't cover kinetics.
     */
    std::ostream& operator<<(std::ostream& s, 
			     IdealReactingGas& mix) {
	std::string r = mix.report(true);
	s << r;
	return s;
    }
}

