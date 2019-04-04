/**
 *  @file EState_RadialDistrib.h  Extension of the save file state for radially distributed
 *                                electrode objects
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ESTATE_RADIALDISTRIB_H
#define _ESTATE_RADIALDISTRIB_H

#include "EState.h"
//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{

class Electrode_SimpleDiff;
class Electrode_DiffTALE;

//==================================================================================================================================
//! Child Class for the EState class for saving states associated with simple diffusion within a single domain
/*!
 *  This class assumes that some of the PhaseList phases are radially distributed. It then saves the concentration of those
 *  species in those phases as a function of the radial coordinate. The radial cell boundaries are also saved.
 *
 */
class EState_RadialDistrib : public EState
{

public:

    //! Default constructor for the base object
    /*!
     *  @param[in]           modelEState         string model type for                
     *  @param[in]           electrodeType       Type of the electrode
     */
    EState_RadialDistrib(const std::string& modelEState = "EState_RadialDistrib", 
                         const std::string& electrodeType = "SimpleDiff");

    //! Copy Constructor
    /*!
     * @param[in]            right               Object to be copied
     */
    EState_RadialDistrib(const EState_RadialDistrib& right);

    //! Destructor is virtual
    virtual ~EState_RadialDistrib();

    //! Assignment operator
    /*!
     *  @param[in]           right               object to be duplicated
     *
     *  @return                                  Returns a reference to the current object
     */
    EState_RadialDistrib& operator=(const EState_RadialDistrib& right);

    //! Duplicator function for this class
    /*!
     *  @param[in]           e                   Pointer to the Electrode object you want the new EState object to service
     *
     *  @return                                  Returns a duplication of the current state as a pointer to the base class
     */
    virtual EState* duplMyselfAsEState(Electrode* e = nullptr) const override;

    //! Initialize the object based on an electrode Base class
    /*!
     *  This call will initialize all of the arrays within this class.
     *  All of the species and phase identification information is created and the class is
     *  readied for use as a state maintainer.
     *
     *  @param[in]           e                   Pointer to the electrode base class
     *
     *  @return                                  Returns 0
     */
    virtual int initialize(const Electrode* const e) override;

    //! Create an indentification XML_Node element for this Electrode EState object
    /*!
     *  @return                                  Returns a malloced XML_Node tree.
     */
    XML_Node* writeIdentificationToXML() const;

    //! Write the ElectrodeState to an XML_Node tree
    /*!
     *  (virtual function)
     *
     *  virtual function, because the base class is called from general code, allowing the child classes to be invoked.
     *
     *  @return                                  pointer to the XML_Node tree
     */
    virtual XML_Node* write_electrodeState_ToXML() const override;

    //! Read the state from the XML_Node given by the argument
    /*!
     *  (Virtual from EState)
     *  This override contains information for the radially distributed part of the read.
     *
     *  @param[in]           xmlEState           electrodeState XML element to be read from.
     *                                           The electrodeState XML element contains the state of the electrode,
     *                                           at a particular time, though the time is not part of the record.
     */
    virtual void readStateFromXML(const XML_Node& xmlEState) override;

    //! Copy the current contents of a SimpleDiff Electrode into this state
    /*!
     *  @param[in]           emp                 Pointer to to a SimpleDiff electrode whose state will be copied into this
     *                                           container.
     *  @param[in]           doFinal             Copy the final quantities. Defaults to true
     *                                           If false it copies the init quantitites.
     */
    void copyElectrode_SimpleDiff_intoState(const ZZCantera::Electrode_SimpleDiff* const emp, bool doFinal = true);

    //! Copy the current contents of a SimpleDiff Electrode into this state
    /*!
     *  @param[in]           emp                 Pointer to to a DiffTALE electrode whose state will be copied into this
     *                                           container.
     *  @param[in]           doFinal             Copy the final quantities. Defaults to true
     *                                           If false it copies the init quantitites.
     */
    void copyElectrode_DiffTALE_intoState(const ZZCantera::Electrode_DiffTALE* const emp, bool doFinal = true);

    //! Set the State of this object from the state of the Electrode object
    /*!
     *  (virtual function)
     *
     *  virtual function, because the base class is called from general code, allowing
     *      the child classes to be invoked.
     *
     *  This function takes the electrode objects _final_ state and copies it into this object.
     *  This function must be carried out before the XML_Node tree is written.
     *
     *  @param e   Pointer to the Electrode object. Note, this class may use dynamic casting
     *             to choose a child object, and then may invoke an error if the match isn't
     *             correct.
     *  @param[in]           doFinal             Copy the final quantities. Defaults to true
     *                                           If false it copies the init quantitites.
     */
    virtual void copyElectrode_intoState(const ZZCantera::Electrode* const e, bool doFinal = true) override;

    //! Set the state of a SimpleDiff electrode object from the contents of this object
    /*!
     *  @param[in]           emp                 Pointer to to a SimpleDiff electrode whose state will be set to the state
     *                                           contained within this object.
     */
    void setStateElectrode_SimpleDiff_fromEState(ZZCantera::Electrode_SimpleDiff* const emp) const;

    //! Set the state of a DiffTALE electrode object from the contents of this object
    /*!
     *  @param[in]           emp                 Pointer to to a DiffTALE electrode whose state will be set to the state
     *                                           contained within this object.
     */
    void setStateElectrode_DiffTALE_fromEState(ZZCantera::Electrode_DiffTALE* const emp) const;

    //! Set the state of the Electrode from the state of this object
    /*!
     *  (virtual function)
     *
     *   Virtual function -> main way to get the process going and the child functions called.
     *
     *   This function is called from the Electrode object or a battery simulator to
     *   initialize the Electrode object from an input file.
     *
     *   This function will set the _final_ state of the Electrode object.
     *   This function will also set all of the internal degrees of freedom such that
     *   the next time step taken by the object will result in the same values and time step
     *   that would have occurred if there hadn't been a restart.
     *
     *  @param  e  Changeable pointer to the base class Electrode object. This function may
     *             do dynamic casting to get the correct child Electrode object.
     */
    virtual void setStateElectrode_fromEState(ZZCantera::Electrode* const e) const override;

    //!  Compare the current state of this object against another guest state to see if they are the same
    /*!
     *    We compare the state of the solution up to a certain number of digits.
     *
     *     @param[in]       ESguest          Guest state object to be compared against
     *     @param[in]       molarAtol        Absolute tolerance of the molar numbers in the state.
     *                                       Note from this value, we can get all other absolute tolerance inputs.
     *     @param[in]       nDigits          Number of digits to compare against
     *     @param[in]       includeHist      Include capacityDischarged and nextDeltaT variables in final bool comparison
     *     @param[in]       printLvl         print level of the routine
     *
     *     @return                           Returns true
     */
    virtual bool compareOtherState(const EState* const ESguest, double molarAtol, int nDigits, 
				   bool includeHist = false, int printLvl = 0) const override;

    /* --------------------------------------------------------------------------------------------------------------------------  */
    //                      D A T A 
    /* --------------------------------------------------------------------------------------------------------------------------  */
protected:

    //! Number of cells in the radial distribution
    int numRCells_;

    //! Define the number of species that are defined to have radially distributed distributions
    //! within the solid
    /*!
     *   Note, for this object there is only one radial distribution. This concept will be enhanced
     *   in later formulations.
     *
     *   There are a few arrays which have numKRSpecies_ as their inner loop dimension.
     */
    int numKRSpecies_;

    //! Number of phases which have radial distributions of their species
    int numSPhases_;

    //! Node position of the mesh - final
    std::vector<double> rnodePos_;  

    //! Vector of cell boundaries
    std::vector<double> cellBoundR_;
 
    //! Total concentration of each of the solid phases that are distributed - global init state
    /*!
     *  concTot_SPhase_Cell_init_init_[numSPhases_ * iCell + iJRPh]
     */
    std::vector<double> concTot_SPhase_Cell_;

    //! total concentration of the solid phases that are distributed - init state
    /*!
     *     concKRSpecies_Cell_init_init_[numKRSpecies_ * iCell + iKRSpecies]
     */
    std::vector<double> concKRSpecies_Cell_;

    //! Vector of solid species defined on the grid of the spherical particle
    /*!
     *
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies;
     *   init_init state value
     */
    std::vector<double> spMoles_KRsolid_Cell_;

    //! This integer describes if the system is current on a Region boundary at the end of a subgrid
    //! integration step
    /*!
     *   We define a region boundary here iff all cells are on that boundary
     *   If it isn't, we set it to -1. If it is, we set it to the region boundary index
     */
    int onRegionBoundary_;

    //! Statement that the Electrode class can access any information in this class
    /*!
     *  NOTE, I'm not sure that this direction of access is needed ATM.
     */
    //friend class ZZCantera::Electrode;
    friend class ZZCantera::Electrode_SimpleDiff;
    friend class ZZCantera::Electrode_DiffTALE;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
