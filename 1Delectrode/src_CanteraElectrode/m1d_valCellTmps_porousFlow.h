/**
 *  @file  m1d_valCellTmps_porousFlow.h  Base class for the structure that holds the values of variables at the cell level.
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef M1D_VALCELLTMPS_POROUSFLOW_H
#define M1D_VALCELLTMPS_POROUSFLOW_H

#include "m1d_porousFlow_dom1D.h"

namespace Zuzax
{
class Electrode;
class Transport;
}

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! Structure containing the value of a single quantity at the cell, as well as its value in the right and left cells
/*!
 *  This is used for simplifying the calculation of cell level quantities.
 */
struct Value_CellRL 
{
    //! Default constructor
    /*!
     *  @param[in]           val                 Value of the quantity. Defaults to 0.
     */
    Value_CellRL(double val = 0.0);

    //! Copy constructor
    /*!
     *  @param[in]           r                   Object to be copied
     */
    Value_CellRL(const Value_CellRL &r);

    //! Assignment operator
    /*!
     *  @param[in]           r                   Object to be copied
     *  @return                                  Returns a current reference to the object
     */
    Value_CellRL& operator=(const Value_CellRL &r);

    //! Fill the values of the current structure given the values of the structure at the cell to the left
    /*!
     *  @param[in]           valCellLeft         corresponding structure from the left cell
     *  @param[in]           newRightVal         New value at the right cell to the current cell
     */
    void fillNextRight(const Value_CellRL &valCellLeft, double newRightVal);

    //! Value of the property at the cell
    double center;

    //! Value of the property at the right cell
    double right;

    //! Value of the property at the left cell
    double left;
};

//==================================================================================================================================
//! Intermediate bookkeeping information for loops over cells
/*!
 *  These values depend on the solution variables.
 *  They must be recalculated whenever the solution changes.
 *
 *  This class consists of public data only, which can be directly accessed.
 *
 *  If the left node doesn't exist then the left node value will be equal to the center value
 *  If the right node doesn't exist then the right node value will be equal to the center value
 */
class valCellTmps
{
public:

    //! Default constructor
    valCellTmps();

    //! Copy constructor
    /*!
     *  @param[in]           r                   Object to be copied
     */
    valCellTmps(const valCellTmps &r);

    //! Virtual default destructor
    virtual ~valCellTmps();

    //! Assignment operator
    /*!
     *  @param[in]           r                   Object to be copied
     *
     *  @return                                  Returns a reference to the object being copied
     */
    valCellTmps& operator=(const valCellTmps &r);

    //----------------------------------------------------- D A T A ----------------------------------------------------------------

    //!  Local value of the axial velocity (m/s) unknown at current cell and at RL cells
    Value_CellRL AxialVeloc;

    //!  Local value of the temperature (Kelvin) unknown at current cell and at RL cells
    Value_CellRL Temperature;

};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
