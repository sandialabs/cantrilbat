/*
 * m1d_porousLiKCl_dom1D.h
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 */

#ifndef M1D_CELLTMPS_POROUSFLOW_H
#define M1D_CELLTMPS_POROUSFLOW_H

#include "m1d_porousFlow_dom1D.h"

namespace Cantera
{
class Electrode;
class Transport;
}


namespace m1d
{

//! Intermediate bookkeeping information for nodes based on loops over cells
class NodeTmps
{
public:
    NodeTmps();

    NodeTmps(const NodeTmps &right);

    ~NodeTmps();

    NodeTmps &operator=(const NodeTmps &right);
    
    //!  Pointer to the NodalVars struct 
    NodalVars *nv;

    //!  Offset of the nodal variables from the start of the solution vector
    size_t index_EqnStart;

    //! Offset of axial displacement
    size_t Offset_Displacement_Axial;

    //! Offset of axial
    size_t Offset_Velocity_Axial; 

    //!  Offset of variables wrt the start of the nodal solution vector.
    size_t Offset_Voltage;

    //!  Offset of variables wrt the start of the nodal solution vector.
    size_t Offset_MoleFraction_Species;

    //! Offset of the temperature equation wrt the start of the nodal solution vector
    size_t Offset_Temperature;

    //! Offset of the pressure equation wrt the start of the nodal solution vector
    size_t Offset_Pressure;

#ifdef MECH_MODEL  
    size_t Offset_Solid_Stress_Axial;
#endif

    //!  Offset of Residual for the current conservation equation
    size_t RO_Current_Conservation;

    //!  Offset of Residual for the electrolyte total mass/mole continuity equation
    size_t RO_Electrolyte_Continuity;
    size_t RO_Species_Eqn_Offset;
    size_t RO_MFSum_offset;
    size_t RO_ChargeBal_offset;

    //!  Offset of Residual for the Enthalpy Conservation equation
    size_t RO_Enthalpy_Conservation;

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
class cellTmps
{
public:
    cellTmps();

    cellTmps(const cellTmps &r);

    virtual ~cellTmps();

    cellTmps& operator=(const cellTmps &r);

    NodalVars* nvLeft_;
    NodalVars* nvCent_;
    NodalVars* nvRight_;

    NodeTmps NodeTmpsLeft_;
    NodeTmps NodeTmpsCenter_;
    NodeTmps NodeTmpsRight_;

    double xdelL_; // Distance from the center node to the left node
    double xdelR_; // Distance from the center node to the right node
    double xdelCell_; // cell width - right boundary minus the left boundary.
    double xCellBoundaryL_; //cell boundary left
    double xCellBoundaryR_; //cell boundary right
};



}

#endif
