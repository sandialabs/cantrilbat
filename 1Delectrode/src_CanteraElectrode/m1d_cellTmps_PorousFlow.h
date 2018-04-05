/*
 * @file m1d_cellTmps_PorousFlow.h Definitions for a base class that handle loops over cells.
 *
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */


#ifndef M1D_CELLTMPS_POROUSFLOW_H
#define M1D_CELLTMPS_POROUSFLOW_H

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

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! Intermediate bookkeeping information for nodes based on loops over cells
/*!
 *  This is a non-virtual class, where all data is public.
 *  Once the mesh is created and not changed, the data in this object stays constant
 */
class NodeTmps
{
public:
    //! Default constructor
    /*!
     *  All offsets are initially set to npos  and remain at anpos if there are no equations at the node of that variable type
     */
    NodeTmps();

    //! Copy constructor
    /*!
     *  @param[in]           right               Object to be copied
     */
    NodeTmps(const NodeTmps &right);

    //! non-virtual Destructor
    ~NodeTmps();

    //! Assignment operator
    /*!
     *  @param[in]           right               Object to be copied
     *
     *  @return                                  returns a reference to the current object
     */
    NodeTmps &operator=(const NodeTmps &right);
    
    //!  Pointer to the NodalVars struct 
    NodalVars *nv;

    //!  Offset of the nodal variables from the start of the global solution vector
    size_t index_EqnStart;

    //! Offset of axial displacement variable wrt the start of the nodal solution vector
    size_t Offset_Displacement_Axial;

    //! Offset of axial velocity variables wrt the start of the nodal solution vector
    size_t Offset_Velocity_Axial; 

    //! Offset of voltage variables wrt the start of the nodal solution vector.
    size_t Offset_Voltage;

    //!  Offset of the mole fraction variables wrt the start of the nodal solution vector.
    size_t Offset_MoleFraction_Species;

    //! Offset of the temperature equation(s) wrt the start of the nodal solution vector
    size_t Offset_Temperature;

    //! Offset of the pressure equation wrt the start of the nodal solution vector
    size_t Offset_Pressure;

    //! Offset of the axial solid stress equation
    size_t Offset_Solid_Stress_Axial;

    //!  Offset of Residual for the current conservation equation
    size_t RO_Current_Conservation;

    //! Offset of Residual for the electrolyte total mass/mole continuity equation 
    //! from the start of the residual vector for that node
    size_t RO_Electrolyte_Continuity;

    //! Offset of the residual for the species equation from the start of the residual vector for that node
    size_t RO_Species_Eqn_Offset;

    //! Offset of the residual for the Mole fraction sum equals 1 residual vector for that node
    size_t RO_MFSum_offset;

    //! Offset of the charge balance residual equation from the start of the resid vector for that node
    size_t RO_ChargeBal_offset;

    //!  Offset of Residual for the Enthalpy Conservation equation
    size_t RO_Enthalpy_Conservation;

};

//==================================================================================================================================
//! Intermediate bookkeeping information for loops over cells
/*!
 *    If we are on a left boundary, there will be no left node. Instead, the left and center node
 *    are really collocated.  In that case NodeLeft_ will be a duplicate of NodeCenter_.
 *    And, nodeLeft member value will be set to zero.
 *
 *    An analogous treatment of right boundaries where there is no right node is also done. 
 */
class cellTmps
{
public:
    //! Constructor
    cellTmps();

    //! Copy constructor
    /*!
     *  @param[in]           r                   Object to be copied
     */
    cellTmps(const cellTmps &r);

    //! Virtual destructor
    virtual ~cellTmps();

    //! Assignment operator
    /*!
     *  @param[in]           r                   Object to be copied
     *
     *  @return                                  Returns a reference to the current object
     */
    cellTmps& operator=(const cellTmps &r);

    //! Pointer to the NodalVars structure to the left of the current cell
    /*!
     *  If there isn't one, then the value is nullptr
     */
    NodalVars* nvLeft_;

    //! Pointer to the NodalVars structure for the node in the center of the current node
    NodalVars* nvCent_;

    //! Pointer to the NodalVars structure to the left of the current cell
    /*!
     *  If there isn't one, then the value is nullptr
     */
    NodalVars* nvRight_;

    //! NodeTmps structure for the node to the left of the current cell
    NodeTmps NodeTmpsLeft_;

    //! NodeTmps structure for the node in the center of the current cell
    NodeTmps NodeTmpsCenter_;

    //! NodeTmps structure for the node to the right of the current cell
    NodeTmps NodeTmpsRight_;

    //! Distance from the center node to the left node
    double xdelL_; 

    //! Distance from the center node to the right node
    double xdelR_; 

    //! Cell width - right boundary minus the left boundary
    double xdelCell_; 

    //! Location of the left cell boundary
    double xCellBoundaryL_; 

    //! Location of the right cell boundary
    double xCellBoundaryR_; 
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------

#endif
