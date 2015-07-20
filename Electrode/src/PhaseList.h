//
/**
 * @file PhaseList.h
 * This file contains definitions for the class PhaseList
 *  (see \ref ExtendedPhaseGroups and
 * class \link Cantera::PhaseList PhaseList\endlink).
 */


#ifndef CT_PHASELIST_H
#define CT_PHASELIST_H

#include "cantera/base/xml.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/Elements.h"

namespace Cantera
{

//!
//!   Class PhaseList is a list of Phases. It's used as a container class
//!
/*!
 *    PhaseList maintains a global species list, containing all
 *    of the species in all of the phases within PhaseList.
 *    The volume phases are first in this list, followed by
 *    the surface phases. The global species index is the
 *    index in this list.
 *
 *    Right now PhaseList doesn't contain any thermodynamic state information.
 *    Therefore, it doesn't include a temperature or pressure variable. It doesn't include
 *    mole numbers of species or a mole fraction vector.
 */
class PhaseList
{

public:
    /***********************************************************************/
    /*                   PUBLIC FUNCTIONS                                  */
    /***********************************************************************/
    /*
     * Constructors.
     *
     *   We default to setting up either no phases or just one
     *   single phase during the construction phase.
     *   We also default to not owning the phase pointers that
     *   this object keeps in its list.
     */
    PhaseList(ThermoPhase* const p = 0, XML_Node* vPhase = 0,
              double vol = 0.0,  bool ownership = true);

     //! Destructor():
    ~PhaseList();

    //! Copy constructor
    /*!
     * @param right PhaseList to copy
     */
    PhaseList(const PhaseList& right);

    //! Assignment operator
    /*!
     *  Does an appropriate copy operation
     *
     *  @param right PhaseList to copy
     *
     *  @return      Returns a reference to the current object
     */
    PhaseList& operator=(const PhaseList& right);


    //!   Add a volumetric phase to the list
    /*!
     *  @param vp             Previously initialized ThermoPhase object to be added to PhaseList
     *  @param vnode          XML_Node pointer to the phase XML Node for the current object
     *  @param canteraFile    File name of the file that contains the phase XML Node
     */
    void addVolPhase(Cantera::ThermoPhase* const vp, Cantera::XML_Node* vNode, std::string canteraFile = "");

    //! Add a volumetric phase
    void addVolPhase(std::string canteraFile);

    //! Add a surface phase to the list of surfaces within the object.
    /*!
     *  @param[in]   vp             Previously initialized ThermoPhase object to be added to PhaseList
     *  @param[in]   sNode          XML_Node pointer to the surface phase XML Node for the current object
     */
    void addSurPhase(Cantera::ThermoPhase* const vp, Cantera::XML_Node* sNode);

    //! Add a surface phase to the PhaseList object given an XML file name
    /*!
     *  Add a surface phase to the PhaseList object. This routine doesn't add the kinetics information to the
     *  PhaseList object.
     *
     *  This routine will pick the first phase listed in the file, open it up and initialize it. It will assume
     *  that it is a surface phase.
     *
     *      @param[in]  canteraFile       String containing the XML file description of a cantera surface phase. The phase
     *                                     may or may not have a kinetics object associated with it.
     */
    void addSurPhase(std::string canteraFile);
    /*
     * getVolPhaseIndex
     *     This routine returns the phase index of a phase. This
     *     number is the index value of the phase in the PhaseList
     *     object.
     */
    int getVolPhaseIndex(const ThermoPhase* const ttp) const;

    int getSurPhaseIndex(const ThermoPhase* const sp) const;

    //! Given a pointer to the ThermoPhase object this returns the global phase index of the phase
    /*!
     *       @param[in]    tp          Pointer to the ThermoPhase object
     * 
     *       @return                   Returns the global phase index.
     */
    int getGlobalPhaseIndex(const ThermoPhase* const tp) const;

    //! Get a pointer to the named phase from the list.
    /*!
     *   @param[in]        phaseName       CString containing the name of the phase
     *
     *   @return                           Returns the pointer to the phase. On errors, a NULL pointer is returned.
     */
    ThermoPhase* getPhase(const char* phaseName) const;

    //! Get the name of the phase
    /*!
     * @param globPhaseIndex     global phase Index of the volume or surface Phase.
     *
     * @return                   Returns the global phase name as a string
     */
    std::string phaseName(int globPhaseIndex) const;

    //! Return the global index of a phase given its name
    /*!
     *  @param   phaseName   Name of the phase
     *  @param   phaseIDAfter
     *
     *  @return  Returns the global index
     */
    int globalPhaseIndex(std::string phaseName, bool phaseIDAfter = false) const;

    //! Get the global species index given the thermophase and the local species index
    /*!
     *   @param ttp                        Pointer to the ThermoPhase index
     *   @param k                          Value of the local species index 
     *
     *   @return                           Returns the value of the global species index.
     */
    int getGlobalSpeciesIndex(const ThermoPhase* const ttp, int k = 0) const;

    //! Get the global species index for a volume or surface phase
    /*!
     * PhaseList maintains a global species list, containing all of the species in all of the phases within PhaseList.
     * The volume phases are first in this list, followed by the surface phases. The global species index is the
     * index in this list.
     *
     * @param[in]          globPhaseIndex        global phase Index of the volume or  surface Phase.
     * @param[in]          k                     local species index within the phase.  Default = 0.
     *
     * @return                               Returns the global species index within the PhaseList object. 
     *                                       If the parameters don't identify an index negative 1 is returned, or the an error is thrown.
     */
    int getGlobalSpeciesIndex(int globPhaseIndex, int k = 0) const;

    //! Return the global index of a species named 'name'
    /*!
     * @param speciesName String name of the species
     * @param phaseName   String name of the phase. This defaults to "", in which case the index
     *                    of the first match is returned
     *
     * @return  Returns the global index of the species.
     */
    int globalSpeciesIndex(const std::string speciesName, const std::string phaseName="") const;

    //! Get the global species index for a volume phase
    /*!
     * PhaseList maintains a global species list, containing all
     * of the species in all of the phases within PhaseList.
     * The volume phases are first in this list, followed by
     * the surface phases. The global species index is the
     * index in this list.
     *
     *
     * @param volPhaseIndex  volume Index of the volume Phase.
     * @param k              local species index within the phase.
     *                       Default = 0.
     *
     * @return
     *   Returns the global species index within the PhaseList object.
     */
    int getGlobalSpeciesIndexVolPhaseIndex(int volPhaseIndex, int k = 0) const;

    //! Get the global species index for a surface phase
    /*!
     * PhaseList maintains a global species list, containing all
     * of the species in all of the phases within PhaseList.
     * The volume phases are first in this list, followed by
     * the surface phases. The global species index is the
     * index in this list.
     *
     * @param surPhaseIndex  surface  Index of the surface Phase.
     * @param k              local species index within the phase.
     *                       Default = 0.
     *
     * @return
     *   Returns the global species index within the PhaseList object.
     */
    int getGlobalSpeciesIndexSurPhaseIndex(int surPhaseIndex, int k = 0) const;


    //! Get the global phase index from the global species index
    /*!
     *  @param   globalSpeciesIndex  species index
     *  @return  return the phase index
     */
    int getPhaseIndexFromGlobalSpeciesIndex(int globalSpeciesIndex) const;

    /*
     *
     */
    void
    getLocalIndecisesFromGlobalSpeciesIndex(int globalSpeciesIndex,
                                            int& phaseIndex,
                                            int& localSpeciesIndex) const;


    //! Compare against other phase lists, coming up with a number for the matching of phases and species
    /*!
     *  This is a messy routine where we try to determine whether two Phaselist's are logically the same, or 
     *  whether one PhaseList is a superset of the other PhaseList. 
     *
     *  @param[in]       plGuest            Pointer to a const PhaseList which will be compared.
     *
     *  @return             Returns an int describing whether the phases are the same:
     *        0  Phase lists are completely the same
     *        1  Phase lists are the same, but in a different phase order. Each ThermoPhase may have different
     *           numbers of species in the phase, but all phases map into another phase.
     *        2  Owning phase list is a superset of the other PhaseList. Each ThermoPhase of the guest PhaseList 
     *           maps into a unique ThermoPhase in the owning Phaselist. The number of species in that ThermoPhase may be different.
     *        3  Guest PhaseList is a superset of the owning PhaseList. Each ThermoPhase of the owing PhaseList
     *           maps into a unique ThermoPhase in the guest Phaselist. The number of species in that ThermoPhase may be different.
     *        4  PhaseLists are not compatible
     *  
     *
     */
    int compareOtherPL(const PhaseList* const plGuest) const;

    //! Return the reference to the %ThermoPhase of a single volume or surface phase
    /*!
     *  @param[in]        globalPhaseIndex  global phase index
     *
     *  @return                    Return a reference to the ThermoPhase object
     */
    ThermoPhase& thermo(int globalPhaseIndex) const;

    //! Return the common Elements object as a const pointer
    /*!
     *      @return              Return a pointer to a const Element object containing a description of the elements in the problem.
     */
    const Elements* getGlobalElements() const;

    //! Return the element name of the eth element
    /*!
     *  @param[in]       e                element index
     *  
     *  @return                           String Name of the element
     */
    std::string elementName(int e) const;

    //!       Index number of an element given its name
    /*!
     *    @param[in]           elementName              String name of the element
     *
     *    @return                                       returns the element index of the element.
     */
    size_t elementIndex(const std::string& elementName) const;

    //! Returns the number of elements in all of the phases, combined
    /*!
     * Note, we do not require all element objects have the same number of elements and the same ordering of elements.
     *
     *  @return   Returns the number of elements in all of the phases, combined.
     */
    int nElements() const;

    //! Return a pointer to the volume phase XML Node for a single volume phase
    /*!
     *    @param iVolIndex   Volume index of the volume phase
     *
     *    @return Returns a pointer to the volume phase XML Node
     */
    XML_Node* volPhaseXMLNode(int iVolIndex) const;

    //! Return a pointer to the surface phase XML Node for a single surface phase
    /*!
     * @param iSurIndex   Surface index of the surface phase
     *
     *   @return                  Return a pointer to the surface phase XML Node for a single surface phase.
     */
    XML_Node* surPhaseXMLNode(int iSurIndex) const;

    //! Returns a pointer to a single volume phase
    /*!
     *  @param iVolIndex   Volume index of the volume phase
     *
     *   @return                           Returns a reference to the %ThermoPhase object for the volume phase.
     */
    ThermoPhase& volPhase(int iVolIndex);

    //! Returns a reference to a single surface phase
    /*!
     *  @param[in]    iSurIndex            Surface index of the surface phase
     *
     *   @return                           Returns a reference to the %ThermoPhase object for the surface phase.
     */
    ThermoPhase& surPhase(int iSurIndex);

    //! Return the total number of phases
    /*!
     *   @return                           returns the total number of phases
     */
    int nPhases() const;

    //! Return the number of volume phases
    /*!
     *     @return                         Return the total number of volume phases
     */
    int nVolPhases() const;

    //! Return the number of surface phases
    /*!
     *     @return                         Return the total number of surface phases
     */
    int nSurPhases() const;

    //! Return the total number of volume species
    /*!
     *     @return                         Return the total number of volume species
     */
    int nVolSpecies() const;

    //! Return the total number of species in all volume and surface phases
    /*!
     *     @return                          Return the total number of species in all volume and surface phases
     */
    int nSpecies() const;

    //! Boolean indicating whether a volume phase has a kinetics object
    /*!
     *    @param[in]    iVolIndex          volume phase index 
     *
     *    @return                          Returns true if the selected volume phase has a kinetics object associated with it.
     */
    bool volPhaseHasKinetics(int iVolIndex) const ;

    //! Boolean indicating whether a surface phase has a surface kinetics object
    /*!
     *    @param[in]    iSurIndex          surface phase index 
     *
     *    @return                          Returns true if the selected surface phase has a kinetics object associated with it.
     */
    bool surPhaseHasKinetics(int iSurIndex) const ;

    //! Return the species name given the global species index
    /*!
     *  @param iGlobSpeciesIndex   global species index within the PhaseList object
     *
     *  @return returns the species name
     */
    std::string speciesName(int iGlobSpeciesIndex) const;

    //! Set the state of all the phases within the PhaseList to a given temperature and pressure
    /*!
     *  Calls the underlying function \link Cantera::ThermoPhase::setState_TP() setState_TP() \endlink
     *  for all phases in the group.
     *
     *   @param[in] temperature     Temperature in Kelvin
     *   @param[in] pressure        pressure in pascals
     */
    void setState_TP(doublereal temperature, doublereal pressure);


    /***********************************************************************/
    /*                BASIC INDEXING DATA                                  */
    /***********************************************************************/
protected:

    
    //! Total number of volume and surface phases defined in the object
    size_t m_NumTotPhases;

    //! Total number of species in all phases
    int m_NumTotSpecies;

    //!  Number of volume phases
    int NumVolPhases_;

    //!  Total number of volume phase species
    int m_totNumVolSpecies;

    //!  Vector of pointers to volume phases existing in the problem
    std::vector<ThermoPhase*> VolPhaseList;

    //! Vector of pointers to volume phase XML trees
    std::vector<XML_Node*> VolPhaseXMLNodes;

    //! Boolean vector indicating whether each volume phase has kinetics
    std::vector<int> VolPhaseHasKinetics;

    //! Number of surface phases
    int m_NumSurPhases;

    //! Total number of surface phase species.
    int m_totNumSurSpecies;

    //! Vector of surface phases existing in the problem
    std::vector<ThermoPhase*> SurPhaseList;

    //! Pointer to the XML nodes corresponding to each surface phase
    std::vector<XML_Node*> SurPhaseXMLNodes;

    //! Boolean vector indicating whether surface phase has kinetics
    std::vector<int> SurPhaseHasKinetics;

    //! Vector of phases in the problem
    std::vector<ThermoPhase*> PhaseList_;

    //! Vector of phase names in the problem
    std::vector<std::string> PhaseNames_;

    //! Number of elements in all of the phases.
    /*!
     * Note, we do not require
     * all element objects have the same number of elements and the
     * same ordering of elements.
     */
    int m_numElements;

    //! Index to the start of the species in a global species list
    /*!
     * length = m_NumTotPhases + 1
     */
    std::vector<int> m_PhaseSpeciesStartIndex;

    //! Name of the file containing the first surface phase encountered
    //std::string CanteraFNSurface;

    /*********************************************************************
     *       DERIVED DATA THAT DEPEND ON INDEPENDENT VARIABLES           *
     *********************************************************************/
private:


    /*********************************************************************
     *      INTERNAL DATA                                                *
     *********************************************************************/
    /*
     *  IOwnPhasePointers: boolean specifying whether this object
     *               owns its phase pointers.
     *
     *   Note: normally, the volume and surface domains own the
     *         phase pointers. However, in thermodynamics programs,
     *         where this object is the top level organizing object,
     *         this object will own phase pointers.
     *
     *     (ownership means that it is responsible for destroying
     *      the object).
     */
    bool IOwnPhasePointers;

    /*
     * Pointer to the element object
     */
    Elements* m_GlobalElementObj;
};

}

#endif
