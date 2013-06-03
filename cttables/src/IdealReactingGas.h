/**
 * @file IdealReactingGas.h
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

/*
 * The IdealReactingGas class combines the IdealGasPhase object with
 * the GasKinetics object. It's just about the same as 
 * IdealGasMix object which is defined in Cantera's include 
 * directory. I just thought the class should be pulled into 
 * the Particles src directory because it's an integral part of 
 * the definition of other important classes such as PartDiscGalerkin
 * class. (I renamed it to avoid conflicts as well). 
 */
#ifndef CP_IDEALREACTINGGAS
#define CP_IDEALREACTINGGAS

#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/GasKinetics.h"

#include <string>
#include <iostream>

namespace Cantera {

    /**
     * IdealReactingGas combines the IdealGasPhase object with
     * the GasKinetics object.
     */
    class IdealReactingGas : 
        public IdealGasPhase, public GasKinetics
    {
    public:
        IdealReactingGas(std::string infile, std::string id="");
        IdealReactingGas(XML_Node& root, std::string id);
        virtual ~IdealReactingGas();
        bool ready() const { return m_ok; }
        friend std::ostream& operator<<(std::ostream& s,
					IdealReactingGas& mix);
    protected:
        bool m_ok;
        XML_Node* m_r;
    private:
	bool m_r_owned;
    };
}
#endif
