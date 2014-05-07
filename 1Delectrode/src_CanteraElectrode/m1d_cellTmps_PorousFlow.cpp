/*
 * m1d_porousLiKCl_dom1D.h
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 */

#include "m1d_cellTmps_PorousFlow.h"

namespace m1d
{
//=====================================================================================================================
NodeTmps::NodeTmps() :
    nv(0),
    index_EqnStart(npos),
    Offset_Velocity_Axial(npos),
    Offset_Voltage(npos),
    Offset_MoleFraction_Species(npos),
    Offset_Temperature(npos),
    RO_Current_Conservation(npos),
    RO_Electrolyte_Continuity(npos),
    RO_Species_Eqn_Offset(npos),
    RO_MFSum_offset(npos),
    RO_ChargeBal_offset(npos),
    RO_Enthalpy_Conservation(npos)
{
}
//=====================================================================================================================
NodeTmps::NodeTmps(const NodeTmps & r) :
    nv(r.nv),
    index_EqnStart(r.index_EqnStart),
    Offset_Velocity_Axial(r.Offset_Velocity_Axial),
    Offset_Voltage(r.Offset_Voltage),
    Offset_MoleFraction_Species(r.Offset_MoleFraction_Species),
    Offset_Temperature(r.Offset_Temperature),
    RO_Current_Conservation(r.RO_Current_Conservation),
    RO_Electrolyte_Continuity(r.RO_Electrolyte_Continuity),
    RO_Species_Eqn_Offset(r.RO_Species_Eqn_Offset),
    RO_MFSum_offset(r.RO_MFSum_offset),
    RO_ChargeBal_offset(r.RO_ChargeBal_offset),
    RO_Enthalpy_Conservation(r.RO_Enthalpy_Conservation)
{
}
//=====================================================================================================================
NodeTmps::~NodeTmps()
{
}
//=====================================================================================================================
NodeTmps& NodeTmps::operator=(const NodeTmps &r)
{
    if (this == &r) {
        return *this;
    }

    nv = r.nv;
    index_EqnStart = r.index_EqnStart;
    Offset_Velocity_Axial = r.Offset_Velocity_Axial;
    Offset_Voltage = r.Offset_Voltage;
    Offset_MoleFraction_Species = r.Offset_MoleFraction_Species;
    Offset_Temperature = r.Offset_Temperature;
    RO_Current_Conservation = r.RO_Current_Conservation;
    RO_Electrolyte_Continuity = r.RO_Electrolyte_Continuity;
    RO_Species_Eqn_Offset = r.RO_Species_Eqn_Offset;
    RO_MFSum_offset = r.RO_MFSum_offset;
    RO_ChargeBal_offset = r.RO_ChargeBal_offset;
    RO_Enthalpy_Conservation = r.RO_Enthalpy_Conservation;

    return *this;
}
//=====================================================================================================================
cellTmps::cellTmps() :
    nvLeft_(0),
    nvCent_(0),
    nvRight_(0),
    NodeTmpsLeft_(),
    NodeTmpsCenter_(),
    NodeTmpsRight_(),
    xdelL_(-1.0),
    xdelR_(-1.0),
    xdelCell_(-1.0),
    xCellBoundaryL_(-1.0),
    xCellBoundaryR_(-1.0)
{
}
//=====================================================================================================================
cellTmps::cellTmps(const cellTmps &r) :
    nvLeft_(r.nvLeft_),
    nvCent_(r.nvCent_),
    nvRight_(r.nvRight_),
    NodeTmpsLeft_(r.NodeTmpsLeft_),
    NodeTmpsCenter_(r.NodeTmpsCenter_),
    NodeTmpsRight_(r.NodeTmpsRight_),
    xdelL_(r.xdelL_),
    xdelR_(r.xdelR_),
    xdelCell_(r.xdelCell_),
    xCellBoundaryL_(r.xCellBoundaryL_),
    xCellBoundaryR_(r.xCellBoundaryR_)
{
}
//=====================================================================================================================
cellTmps::~cellTmps()
{
}
//=====================================================================================================================
cellTmps& cellTmps::operator=(const cellTmps &r)
{
    if (this == &r) {
        return *this;
    }

    nvLeft_ = r.nvLeft_;
    nvCent_  = r.nvCent_;
    nvRight_  = r.nvRight_;
    NodeTmpsLeft_ = r.NodeTmpsLeft_;
    NodeTmpsCenter_ = r.NodeTmpsCenter_;
    NodeTmpsRight_ = r.NodeTmpsRight_;
    xdelL_ = r.xdelL_;
    xdelR_ = r.xdelR_;
    xdelCell_ = r.xdelCell_;
    xCellBoundaryL_ = r.xCellBoundaryL_;
    xCellBoundaryR_ = r.xCellBoundaryR_;

    return *this;
}
//=====================================================================================================================
}

