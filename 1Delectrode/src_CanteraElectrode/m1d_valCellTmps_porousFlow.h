/*
 * m1d_porousLiKCl_dom1D.h
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 */

#ifndef M1D_VALTMPS_POROUSFLOW_H
#define M1D_VALTMPS_POROUSFLOW_H

#include "m1d_porousFlow_dom1D.h"

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
class Electrode;
class Transport;
}


namespace m1d
{

struct Value_CellRL {
    Value_CellRL(double val = 0.0);
    Value_CellRL(const Value_CellRL &r);
    Value_CellRL& operator=(const Value_CellRL &r);
    void fillNextRight(const Value_CellRL &valCellLeft, double newVal);
    double center;
    double right;
    double left;
};

// --------------------------------------------------------------------------------------------------

//! Intermediate bookkeeping information for loops over cells
/*!
 *    If we are on a left boundary, there will be no Left Node. Instead, the left and center node
 *    are really collocated.  In that case NodeLeft_ will be a duplicate of NodeCenter_.
 *    And, nodeLeft member value will be set to zero.
 *
 *    An analogous treatment of right boundaries where there is no right node is also done. 
 */
class valCellTmps
{
public:
    valCellTmps();

    valCellTmps(const valCellTmps &r);

    virtual ~valCellTmps();

    valCellTmps& operator=(const valCellTmps &r);

    Value_CellRL AxialVeloc;

    Value_CellRL Temperature;


};



}

#endif
