/**
 *  @file Electrode_Polarization.h
 *     Headers for the declarations of the Electrode polarization classes, used to model 
 *     Electrode processes
 *     (see \ref electrode_mgr and class \link Zuzax::Electrode Electrode\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ELECTRODE_POLARIZATION_H
#define _ELECTRODE_POLARIZATION_H


#include "cantera/base/config.h"

#include "cantera/base/ct_defs.h"




//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{

class Electrode;
//==================================================================================================================================

//! The enum catalogs all of the types of polarization loss mechanisms that can occur within batteries that cantrilbat/aria 
//! can catalog
/*!
 *  
 */
enum Polarization_Loss_Enum {
    UNKNOWN_PL = -1,

    //! Surface OCV
    SURFACE_OCV_PL = 0,

    //! Voltage loss through a surface overpotential term
    SURFACE_OVERPOTENTIAL_PL,

    //! Voltage loss through a concentration gradient in the electrolyte
    CONC_LOSS_LYTE_PL,

    //! Voltage loss through voltage drop - electrolyte within the separator
    VOLT_LOSS_LYTE_SEP_PL,

    //! Voltage loss through a concentration gradient in the electrolyte within the separator
    CONC_LOSS_LYTE_SEP_PL,

    //! Voltage loss through voltage drop - electrolyte
    VOLT_LOSS_LYTE_PL,

    //! Concentration loss through a boundary layer near the surface of the electrolyte interface
    CONC_LOSS_BL_LYTE_PL,

    //! Resistance layer voltage loss at the surface of the electrolyte
    RESISTANCE_FILM_PL,

    //! Loss through conduction of electrons through the solid phase of the electrode
    ELECTRICAL_CONDUCTION_LOSS_PL,

    //! Loss through conduction of electrons through the Collectors
    ELECT_COND_COLLECTORS_LOSS_PL,

    //! Loss through solid phase diffusion within the electrode material
    SOLID_DIFF_CONC_LOSS_PL,

    //! Concentration polarization  from mixed ionic conduction through solid film, i.e., SEI.
    CONC_LOSS_SOLID_FILM_PL,

    //! voltage drop from mixed ionic conduction through solid film, i.e., SEI.
    VOLT_LOSS_SOLID_FILM_PL,

    //! Other loss mechanisms not yet identified
    OTHER_PL
};

//==================================================================================================================================
//! Structure to store the type and amount of voltage loss computed for a polarization loss
/*!
 *  The types are identified by the enum Polarization_Loss_Enum variable.
 *  The signs are always negative for a loss in power of the resulting electron produced by a battery in discharge, except for
 *  the OCV term. For the OCV term the cathode side is (phi_metal - phi_soln),  while the anode side is the negative of the
 *  cathode, (phi_soln - phi_metal).
 *
 *  Signs will be positive for charging polarization losses . 
 */
struct VoltPolPhenom
{
    //! Default constructor
    /*!
     *  @param[in]           pltype              Polarization type enum
     *  @param[in]           volts               volt drop
     *  @param[in]           region              Number of the region, to identify the region.
     */
    VoltPolPhenom(enum Polarization_Loss_Enum pltype, int region, double volts) :
        ipolType(pltype),
        regionID(region),
        voltageDrop(volts)
    {
    }

    //! Enum identifying the effect that we are cataloging
    enum Polarization_Loss_Enum ipolType;

    //! Region where the phenomena occurred
    /*!
     *  -1  -  unspecified
     *   0  -  anode
     *   1  -  separator
     *   2  -  cathode
     *   3  -  anode current collector
     *   4  -  cathode current collector
     */
    int regionID;

    //! Drop in voltage due to effect.
    /*!
     *  When indicating an effect, this will be the V_cathodeSide - V_anode_side. 
     *  Therefore, when you add up all of the effects from Cathode to anode, you will get V_cathode - V_anode.
     *
     *  Thus, all polarization effects will be a negative voltage for polarization losses in normal operation.
 
     *  When indicating OCV, this will be the voltage of the electrode (phiMetal - phiSoln) for the Cathode.
     *  It will be entered as - (phiMetal - phiSoln) for the anode.
     *  When the two half-cells are added together, OCV = V_Cathode - V_anode = voltageDrop_Cathode + voltageDrop_anode.
     *
     *  Therefore the following relation will hold for battery discharge
     *        phi_Cathode - phi_Anode = OCV - sum(polarizationLosses)  = sum (voltageDrop_i)
     *  The following will hold for battery charging direction
     *        phi_Cathode - phi_Anode = OCV + sum(polarizationLosses)  = sum (voltageDrop_i)
     * 
     */
    double voltageDrop;
};

//===================================================================================================================================
//! Structure for communicating polarization results
/*!
 *  Structure is for one reaction on one surface only. On many electrodes there will only be one surface active at a time.
 *  and there may be multiple reactions on each surface producing electrons. However, the only way to distinguish 
 *  the results is to differentiate by reaction. This is the only way to get the actual value of the overpotential for a reaction.
 *
 *  Production Development:
 *  I can imagine there will be loops on a particular surface where electrons are created and then destroyed by separate reactions.
 *  The net result will be a surface reactions which change speciation in the bulk or adsorbate phases.
 *  Pre-production of these loops to eliminate these net-zero electron production cases should occur. I can delay development
 *  of this code until the situation occurs in practice.
 */
struct PolarizationSurfRxnResults {

    //! Domain number of the electrode object
    /*!
     *  This number refers to the domain number within 1DElectrode.
     */
    int electrodeDomainNumber_ = 0;

    //! Cell number within the domain
    int electrodeCellNumber_ = 0;

    //! Shallow pointer to the Electrode object
    Electrode* ee_;

    //!  Index of the reacting surface within the Electrode that the summary is for
    size_t isurf_ = npos;

    //! Index of the surface reaction on that surface
    /*!
     *  Reaction index that is producing electrons
     */
    size_t iRxnIndex_ = npos;

    //! Total current through the surface that is using this particular reaction on this surface
    /*!
     *  Units: amps
     */
    double icurrSurf = 0.0;

    //! Vector of physical-based voltage losses
    std::vector<VoltPolPhenom> voltsPol_list;

    //! Value of the open circuit voltage for the surface reaction index
    double ocvSurfRxn = 0.0;

    //! Value of the open circuit voltage for the surface (assuming adsorbate coverage is at pseudo-equilibrium)
    double ocvSurf = 0.0;

    //! Adjustments to the ocvSurfRxn are made when we identify physical phenomena associated with the
    //! change in the OCV, such as electrolyte concentration polarization or solid state diffusion polarization.
    /*!
     *  We keep track of the adjusted OCV here.
     */ 
    double ocvSurfRxnAdj = 0.0;

    //! Value of the voltage for the electrode through which the electrons go through
    /*!
     *  This is phiMetal - phiSoln always no matter if it is the anode or the cathode
     */
    double VoltageElectrode = 0.0;

    //! Electric potential of the metal at the electrode
    /*!
     *  Need to know the absolute electric potentials in order to interpret the results.
     *
     *  phiSoln = phiMetal - VoltageElectrode
     */
    double phiMetal = 0.0;

    //! Value of the electric potential that is associated with the middle of the domain.
    /*!
     *  This is one end of the comparison voltage.
     */
    double phi_lyteSpoint = 0.0;

    //! Value for the total voltage drops accounted for by the voltsPol_List
    /*!
     *  We'll be post-processing this structure to add in voltage drops
     *  The voltages are processed from cathode to the anode.
     *  When fully formed the following condition holds:
     *         VoltageTotal = phi_Anode - phi_lyteSpoint
     */
    double VoltageTotal = 0.0;

    //! Default constructor
    /*!
     *  @param[in]           electrodeDomainNumber  Domain number
     *  @param[in]           electrodeCellNumber  cell number
     *  @param[in]           ee                  Pointer to the electrode object
     *  @param[in]           surfIndex           Surface index. Defaults to npos
     *  @param[in]           rxnIndex            Reaction index on the surface. Defaults to npos
     */
    PolarizationSurfRxnResults(int electrodeDomainNumber, int electrodeCellNumber, Electrode* ee, size_t surfIndex = npos, 
                               size_t rxnIndex = npos); 

    //! Add the polarization losses due to the electrical conduction through the Electrode's solid portion to the current collector
    /*!
     *  This gets added to the Electrode object's records.
     *
     *  @param[in]           phiCurrentCollector Electric potential of the current collector
     *  @param[in]           region              Region of the solid condition. Two possibilities
     *                                               - 0  anode
     *                                               - 2  cathode
     */
    void addSolidPol(double phiCurrentCollector, int region, bool dischargeDir);

    //! Add the polarization losses due to the electrical conduction through the Electrode Current Collector
    //! to the Electrode terminal
    /*!
     *  This gets added to the Electrode object's records.
     *
     *  @param[in]           phiTerminal         Electric potential at the terminal of the battery.
     *  @param[in]           phiCurrentCollector Electric potential of the current collector - electrode interface
     *  @param[in]           region              Region of the solid condition. Two possibilities
     *                                               - 0  anode
     *                                               - 2  cathode
     *  @param[in]           dischargeDir        Discharge direction:
     *                                               True if discharging battery
     *                                               False if charging the battery
     */
    void addSolidCCPol(double phiTerminal, double phiCurrentCollector, int region, bool dischargeDir);

    //! Add the polarization losses due to the ionic conduction through the electrolyte from the electrode
    //! to the edge of the electrode - separator material interface.
    /*!
     *  This gets added to the Electrode object's records. This assigns the Polarization loss type,
     *  VOLT_LOSS_LYTE_PL.
     *
     *  @param[in]           phiLyteElectrode    Electric potential of the electrolyte at the electrode position.
     *  @param[in]           phiLyteBoundary     Electric potential of the electrolyte at the electrode material - separator
     *                                           material boundary.
     *  @param[in]           region              Region of the solid condition. Two possibilities
     *                                               - 0  anode
     *                                               - 2  cathode
     *  @param[in]           discharngeDir       Discharge direction:
     *                                               True if discharging battery
     *                                               False if charging the battery
     */
    void addLyteCondPol(double phiLyteElectrode, double phiLyteBoundary, int region, bool dischargeDir);

    //! Add the polarization losses due to the ionic conduction through the electrolyte from the electrode
    //! to the edge of the electrode - separator material interface.
    /*!
     *  This gets added to the Electrode object's records. This assigns the Polarization loss type,
     *  VOLT_LOSS_LYTE_PL.
     *
     *  @param[in]           phiLyteElectrode    Electric potential of the electrolyte at the electrode position.
     *  @param[in]           phiLyteBoundary     Electric potential of the electrolyte at the electrode material - separator
     *                                           material boundary.
     *  @param[in]           phiLyte_Spoint      Electric potential of the electrolyte at middle of the separator
     *  @param[in]           region              Region of the solid condition. Two possibilities
     *                                               - 0  anode
     *                                               - 2  cathode
     *  @param[in]           discharngeDir       Discharge direction:
     *                                               True if discharging battery
     *                                               False if charging the battery
     */
    void addLyteCondPol_Sep(double phiLyteBoundary, double phiLyte_Spoint, int region, bool dischargeDir);

    //! Add the concentration polarization losses due to the ionic conduction through the electrolyte from the electrode
    //! to the edge of the electrode - separator material interface.
    /*!
     *  This gets added to the Electrode object's records. This assigns the Polarization loss type,
     *  CONC_LOSS_LYTE_PL.
     *
     *  @param[in]         state_Lyte_Electrode   Pointer to the state vector of the electrolyte at the Electrode object
     *  @param[in]     state_Lyte_SeparatorBdry   Pointer to the state vector of the electrolyte at the electrode 
     *                                           material - separator material boundary.
     *  @param[in]           region              Region of the solid condition. Two possibilities:
     *                                               - 0  anode
     *                                               - 2  cathode
     *  @param[in]           discharngeDir       Discharge direction:
     *                                               True if discharging battery
     *                                               False if charging the battery
     */
    void addLyteConcPol(double* state_Lyte_Electrode, double* state_Lyte_SeparatorBdry, int region, bool dischargeDir);

    //! Add the concentration polarization losses due to the ionic conduction through the electrolyte from the electrode
    //! to the edge of the electrode - separator material interface.
    /*!
     *  This gets added to the Electrode object's records. This assigns the Polarization loss type,
     *  CONC_LOSS_LYTE_PL.
     *
     *  @param[in]     state_Lyte_SeparatorBdry  Pointer to the state vector of the electrolyte at the electrode
     *                                            material - separator  material boundary.
     *  @param[in]     state_Lyte_SeparatorMid   Pointer to the state vector f the electrolyte at the
     *                                            separator midpoint 
     *  @param[in]           region              Region of the solid condition. Two possibilities
     *                                               - 0  anode
     *                                               - 2  cathode
     *  @param[in]           discharngeDir       Discharge direction:
     *                                               True if discharging battery
     *                                               False if charging the battery
     */
    void addLyteConcPol_Sep(double* mf_Lyte_SeparatorBdry, double* mf_Lyte_SeparatorMid, int region, bool dischargeDir);

    //! Add the concentration polarization losses due to the ionic conduction through the electrolyte from the electrode
    //! to the edge of the electrode - separator material interface.
    /*!
     *  This gets added to the Electrode object's records. This assigns the Polarization loss type,
     *  SOLID_DIFF_CONC_LOSS_PL.
     *
     *  @param[in]           mf_Lyte_Electrode   pointer to the mole fractions of the electrolyte at the electrode object
     *  @param[in]           mf_Lyte_Electrode   pointer to the mole fractions of the electrolyte at the electrode object
     *                                            e electrode material - separator
     *                                           material boundary.
     *  @param[in]           region              Region of the solid condition. Two possibilities
     *                                               - 0  anode
     *                                               - 2  cathode
     *  @param[in]           discharngeDir       Discharge direction:
     *                                               True if discharging battery
     *                                               False if charging the battery
     */
    void addSolidElectrodeConcPol(double* mf_OuterSurf_Electrode, double* mf_Average_Electrode, int region, bool dischargeDir);

    //! Add one PolarizationSurfRxnResults struct into another one
    /*!
     *  Here, I've made the current additive. This has to be checked for compatibility.
     *
     *  @param[in]           sub                 Struct to be added into the current struct
     */
    void addSubStep(struct PolarizationSurfRxnResults& sub);

    
    //! Once we have calculated all of the voltage losses, all records should have the same value of VoltageTotal
    /*!
     *  @param[in]            gvoltageTotal      Input voltage that this record will be checked against.
     *
     *  @return                                  true if the record is consistent
     */
    bool checkConsistency(const double gvoltageTotal);

};
//===================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

