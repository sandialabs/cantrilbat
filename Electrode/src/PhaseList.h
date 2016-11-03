/**
 * @file PhaseList.h
 * This file contains definitions for the class PhaseList (see \ref ExtendedPhaseGroups and
 * class \link Zuzax::PhaseList PhaseList\endlink).
 */

#ifndef CT_PHASELIST_H
#define CT_PHASELIST_H

#include "cantera/base/xml.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/Elements.h"
#include "cantera/base/config.h"
#include "cantera/base/Array.h"

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace 
namespace  Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
//!
//!   Class PhaseList is a list of Phases. It's used as a container class.
/*!
 *    %PhaseList maintains a global species list, containing all of the species in all of the phases within %PhaseList.
 *    The volume phases are first in this list, followed by  the surface phases. The global species index is the
 *    index in this list.  %PhaseList doesn't contain any thermodynamic state information.
 *    Therefore, it doesn't include a temperature or pressure variable. It doesn't include
 *    mole numbers of species or a mole fraction vector, either.
 *
 *    PhaseList class operates strictly on pointers of ThermoPhases. It maintains lists of pointers. Therefore, when
 *    the pointer to a ThermoPhase changes, the PhaseList must be reconstructed.
 *
 *    A note about the nomenclature. The following are official names for the indecises provided by this class:
 *
 *  <TABLE>
 *  <TR>
 *     <TD> Volume phase index number.</TD>
 *     <TD> Volumes phases have their own index number, ranging from 0 to the number of volume phases,  NumVolPhases_.</TD>
 *  </TR>
 *  <TR>
 *     <TD> Surface phase index number.</TD>
 *     <TD> Surface phases have their own index number, ranging from 0 to the number of surface phases,  NumSurPhases_. 
 *          Edge phases are included in the surface phase list. They are treated no differently than surface phases. 
 *          Therefore, 1 and 2 dimensional ThermoPhases are lumped together.</TD>
 *  </TR>
 *  <TR>
 *     <TD> Global phase index number.</TD>
 *     <TD> All phases are grouped first by volume and then by surface. This index is called the global phase
 *          index number. </TD>
 *  </TR>
 *  <TR>
 *     <TD> Local species index. </TD>
 *     <TD> The local species index is the species index of a species within its own %ThermoPhase. </TD>
 *  </TR>
 *  <TR>
 *     <TD> Global species index. </TD>
 *     <TD> %PhaseList maintains a vector of species located in all of the volume and surface phases within its
 *          structure. The global species index is the index into that vector. The species are contiguous first by their phase ids.
 *          The phases are organized according to the global phase index number, so that all of the volume phases
 *          are placed first in the phase list. </TD>
 *  </TR>
 *  </TABLE>
 *
 *  todo: add routines to take away phases.
 *        add routines to resize and copy old vectors to new vectors after phase addition or removal.
 *        Continue optimization.
 *        Change addPhase routines to return the global phase index.
 *        Change interface routine names one by one to make them conformal to a standard.
 *             ' get' name is used for return vectors in the parameter list.
 *             Determine whether to use globalSpeciesIndex() or speciesIndex() as a general naming convention.
 *
 */
class PhaseList
{
public:

    //! Constructor
    /*!
     *   We default to setting up either no phases or just one single phase during the construction phase.
     *   We also default to owning the phase pointers that this object keeps in its list.
     *
     *  @param[in]     Iownphases             Boolean indicating ownership of all of the ThermoPhase objects.
     */
    PhaseList(bool Iownphases = true);

    //! Destructor
    ~PhaseList();

    //! Copy constructor
    /*!
     * @param[in]     right                   Reference to the PhaseList to copy into this object
     */
    PhaseList(const PhaseList& right);

    //! Assignment operator
    /*!
     *  Does an appropriate copy operation
     *
     *  @param        right                    PhaseList to copy
     *
     *  @return                                Returns a reference to the current object
     */
    PhaseList& operator=(const PhaseList& right);

    //! Returns whether the PhaseList is ordered by dimensionality
    /*!
     *  @return                                  Returns whether the phase is ordered by dimensions
     */
    bool orderedByDims() const;
 
    //!   Add a volumetric phase to the list of phases held by this object
    /*!
     *   @param[in]          vp                  Previously initialized ThermoPhase object to be added to PhaseList
     *   @param[in]          vNode               XML_Node pointer to the phase XML Node for the current object
     *
     *   @param[in]          orderByDims         If true, the the volume phase will be placed at the end of the volume phases
     *                                           and before surface and edge phases. Defaults to true.
     */
    void addVolPhase(ThermoPhase* const vp, XML_Node* vNode = 0, bool orderByDims = true);

    //!  Add a phase to the list of phases held by this object
    /*!
     *  This routine will pick the first phase listed in the file, open it up and initialize it. It will assume
     *  that it is a volume, surface, or edge phase depending on its dimensionality. 
     *
     *  @param[in]     vp                      Previously initialized ThermoPhase object to be added to PhaseList
     *  @param[in]     vNode                   XML_Node pointer to the phase XML Node for the current object
     */
    void addPhase(ThermoPhase* const vp, XML_Node* vNode = 0);

    //! Add a volume phase to the volume PhaseList object given an XML file name
    /*!
     *  Add a volume phase to the PhaseList object. This routine doesn't add the kinetics information to the
     *  PhaseList object.
     *
     *  This routine will pick the first phase listed in the file, open it up and initialize it. It will assume
     *  that it is a volume phase. 
     *
     *  \todo Check that the phase added is unique and has a dimension of three
     *
     *   @param[in]    canteraFile         String containing the XML file description of a %Cantera volume phase. The phase
     *                                     may or may not have a kinetics object associated with it.
     *
     *   @param[in]    volID               String containing the ID() of the volume phase to be added. The default is "".
     *                                     If the default is used, the first volume phase in the file is used.
     *
     *   @param[in]          orderByDims         If true, the the volume phase will be placed at the end of the volume phases
     *                                           and before surface and edge phases. Defaults to true.
     *
     *  @return                            Returns a pointer to the ThermoPhase object
     */
    ThermoPhase* addVolPhase(const std::string& canteraFile, const std::string& volID = "", bool orderByDims = true);

    //! Add a surface phase to the list of surfaces within the object.
    /*!
     *  @param[in]   vp             Previously initialized ThermoPhase object to be added to PhaseList
     *  @param[in]   sNode          XML_Node pointer to the surface phase XML Node for the current object
     */
    void addSurPhase(ThermoPhase* const vp, XML_Node* sNode = 0);

    //! Add a surface phase to the PhaseList object given an XML file name
    /*!
     *  Add a surface phase to the PhaseList object. This routine doesn't add the kinetics information to the
     *  PhaseList object.
     *
     *  This routine will pick the first phase listed in the file, open it up and initialize it. It will assume
     *  that it is a surface phase. 
     *
     *   @param[in]    canteraFile         String containing the XML file description of a cantera surface phase. The phase
     *                                     may or may not have a kinetics object associated with it.
     *
     *   @param[in]    surID               String containing the ID() of the surface phase to be added. The default is "".
     *                                     If the default is used, the first surface phase in the file is used.
     *
     *  @return                                  Returns a pointer to the ThermoPhase object
     */
    ThermoPhase* addSurPhase(const std::string& canteraFile, const std::string& surID = "");

    //! Add an edge phase to the list of edges within the object.
    /*!
     *  @param[in]   ep             Previously initialized ThermoPhase object to be added to PhaseList
     *  @param[in]   eNode          XML_Node pointer to the edge phase XML Node for the current object
     */
    void addEdgePhase(ThermoPhase* const ep, XML_Node* eNode = 0);

    //! Add an edge phase to the PhaseList object given an XML file name
    /*!
     *  Add an edge phase to the PhaseList object. This routine doesn't add the kinetics information to the
     *  PhaseList object.
     *
     *  This routine will pick the first phase listed in the file, open it up and initialize it. It will assume
     *  that it is a edge phase. 
     *
     *
     *   @param[in]    canteraFile         String containing the XML file description of a cantera edge phase. The phase
     *                                     may or may not have a kinetics object associated with it.
     *
     *   @param[in]    edgeID               String containing the ID() of the surface phase to be added. The default is "".
     *                                     If the default is used, the first surface phase in the file is used.
     *
     *  @return                            Returns a pointer to the ThermoPhase object
     */
    ThermoPhase* addEdgePhase(const std::string& canteraFile, const std::string& edgeID = "");

    //! Add a phase to the PhaseList object given an XML file name
    /*!
     *  Add a volume, surface, or edge phase to the PhaseList object. This routine doesn't add the kinetics information to the
     *  PhaseList object.
     *
     *  This routine will pick the first phase listed in the file, open it up and initialize it. It will assume
     *  that it is a volume, surface, or edge phase depending on its dimensionality. 
     *
     *   @param[in]    canteraFile         String containing the XML file description of a %Cantera volume phase. The phase
     *                                     may or may not have a kinetics object associated with it.
     *
     *   @param[in]    volID               String containing the ID() of the volume phase to be added. The default is "".
     *                                     If the default is used, the first volume phase in the file is used.
     *
     *  @return                            Returns a pointer to the ThermoPhase object
     */
    ThermoPhase* addPhase(const std::string& canteraFile, const std::string& volID = "");

    //! Get the volume phase index of a volume phase given its %ThermoPhase pointer
    /*!
     *     This routine returns the phase index of a phase. This number is the index value of the phase in the VolPhaseList object.
     *
     *     @param[in]     ttp             Pointer to the volume phase ThermoPhase object
     *
     *     @return                        Returns the volume phase index number. Returns npos if the phase is not found.
     */
    size_t getVolPhaseIndex(const Zuzax::ThermoPhase* const ttp) const;

    //! Get the surface phase index given the pointer to the surface phase ThermoPhase object
    /*!
     *     @param[in]     sp              Pointer to a surface ThermoPhase
     *
     *     @return                        Returns the surface phase index number. Returns npos if the phase is not found.
     */
    size_t getSurPhaseIndex(const ThermoPhase* const sp) const;

    //! Get the edge phase index given the pointer to the edge phase ThermoPhase object
    /*!
     *     @param[in]     ep              Pointer to an edge ThermoPhase
     *
     *     @return                        Returns the edge phase index number. Returns npos if the phase is not found.
     */
    size_t getEdgePhaseIndex(const ThermoPhase* const ep) const;

    //! Given a pointer to the ThermoPhase object this returns the global phase index of the phase
    /*!
     *       @param[in]    tp          Pointer to the ThermoPhase object
     * 
     *       @return                   Returns the global phase index. If not found, this returns npos
     *
     */
    size_t getGlobalPhaseIndex(const ThermoPhase* const tp) const;

    //! Get a %ThermoPhase pointer to the named phase from the list.
    /*!
     *   @param[in]        phaseName       CString containing the name of the phase
     *
     *   @return                           Returns the pointer to the phase. On errors, a NULL pointer is returned.
     */
    ThermoPhase* getPhase(const char* phaseName) const;

    //! Get a %ThermoPhase pointer to the named phase from the list.
    /*!
     *   @param[in]        phaseName       string containing the name of the phase
     *
     *   @return                           Returns the pointer to the phase. On errors, a NULL pointer is returned.
     */
    ThermoPhase* getPhase(const std::string& phaseName) const;

    //! Get the name of the phase given its global phase index
    /*!
     * @param[in]            iphGlob             global phase Index of the volume or surface Phase.
     *
     * @return                                   Returns the global phase name as a string
     */
    std::string phaseName(size_t iphGlob) const;

    
    //! Get the ID of the phase given its global phase index
    /*!
     * @param[in]            iphGlob             global phase Index of the volume or surface Phase.
     *
     * @return                                   Returns the global phase ID() as a string
     */
    std::string phaseID(size_t iphGlob) const;

    //! Return the global index of a phase given its name
    /*!
     *  @param[in]           phaseName           Name of the phase
     *  @param[in]           phaseIDAfter
     *
     *  @return                                  Returns the global index
     */
    size_t globalPhaseIndex(const std::string& phaseName, bool phaseIDAfter = false) const;

    //! Get the global species index given the ThermoPhase and the local species index
    /*!
     *   @param[in]          ttp                 Pointer to the ThermoPhase index
     *   @param[in]          k                   Value of the local species index 
     *
     *   @return                                 Returns the value of the global species index.
     */
    size_t getGlobalSpeciesIndex(const ThermoPhase* const ttp, size_t k = 0) const;

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
     *                                       If the parameters don't identify an index negative 1 is returned, or the 
     *                                       an error is thrown.
     */
    size_t getGlobalSpeciesIndex(size_t globPhaseIndex, size_t k = 0) const;

    //! Return the global index of a species named 'name'
    /*!
     * @param[in]            speciesName         String name of the species
     * @param[in]            phaseName           String name of the phase. This defaults to "", in which case the index of the first match is returned
     *
     * @return  Returns the global index of the species.
     */
    size_t globalSpeciesIndex(const std::string& speciesName, const std::string phaseName="") const;

    //! Get the global species index for a volume phase
    /*!
     * PhaseList maintains a global species list, containing all
     * of the species in all of the phases within PhaseList.
     * The volume phases are first in this list, followed by
     * the surface phases. The global species index is the index in this list.
     *
     *
     * @param volPhaseIndex  volume Index of the volume Phase.
     * @param k              local species index within the phase.                 Default = 0.
     *
     * @return
     *   Returns the global species index within the PhaseList object.
     */
    size_t getGlobalSpeciesIndexVolPhaseIndex(size_t volPhaseIndex, size_t k = 0) const;

    //! Get the global species index for a surface phase
    /*!
     * PhaseList maintains a global species list, containing all
     * of the species in all of the phases within PhaseList.
     * The volume phases are first in this list, followed by
     * the surface phases. The global species index is the
     * index in this list.
     *
     * @param surPhaseIndex  surface  Index of the surface Phase.
     * @param k              local species index within the phase.   Default = 0.
     *
     * @return                                   Returns the global species index within the PhaseList object.
     */
    size_t getGlobalSpeciesIndexSurPhaseIndex(size_t surPhaseIndex, size_t k = 0) const;

    //! Get the global phase index from the global species index
    /*!
     *  @param[in]           kGlob               Global species index
     *  @return                                  return the phase index
     */
    size_t getPhaseIndexFromGlobalSpeciesIndex(size_t kGlob) const;

    //! Returns the dimensional phase index given the global phase index
    /*!
     *   @param[in]          iphGlob             Global phase index
     *   @param[out]         dimPtr              Pointer to an int, which on return contains the global dimension
     *
     *   @return                                 Returns the dimensional index of the phase 
     *                                           In other words if iphGlob is the second volume phase in the list
     *                                           then it returns 1 and *dimPtr= 3
     */
    size_t dimPhaseIndexFromGlobalPhaseIndex(size_t iphGlob, int *dimPtr = 0);

    //!  Get the local species and global phase index given the global species index
    /*!
     *   @param[in]          kGlob               Global species index
     *   @param[out]         iphGlob             On return this contains the global phase index
     *   @param[out]         k                   On return this contains the local species index
     */
    void getLocalIndecisesFromGlobalSpeciesIndex(size_t kGlob, size_t& iphGlob, size_t& k) const;

    //! Compare against other phase lists, coming up with a number for the matching of phases and species
    /*!
     *  This is a messy routine where we try to determine whether two Phaselist's are logically the same, or 
     *  whether one PhaseList is a superset of the other PhaseList. 
     *
     *  @param[in]           plGuest             Pointer to a const PhaseList which will be compared.
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
     */
    int compareOtherPL(const PhaseList* const plGuest) const;

    //! Return the reference to the %ThermoPhase of a single volume or surface phase
    /*!
     * 
     *  @param[in]        iphGlob           global phase index
     *
     *  @return                             Return a reference to the ThermoPhase object
     */
    ThermoPhase& thermo(int iphGlob) const;

    //! Return the reference to the %ThermoPhase of a single volume or surface phase
    /*!
     *  @param[in]        iphGlob           global phase index
     *
     *  @return                             Return a reference to the ThermoPhase object
     */
    ThermoPhase& thermo(size_t iphGlob) const;

    //! Return the reference to the %ThermoPhase of a single volume or surface phase
    /*!
     *  @param[in]        phaseNamePtr           Character string name of the phase
     *
     *  @return                                  Return a reference to the ThermoPhase object
     */
    ThermoPhase& thermo(const char* const phaseNamePtr) const;

    //! Return the reference to the %ThermoPhase of a single volume or surface phase
    /*!
     *  @param[in]        phaseName  String name of the phase
     *
     *  @return                      Return a reference to the ThermoPhase object
     */
    ThermoPhase& thermo(const std::string& phaseName) const;

    //! Return the reference to the %ThermoPhase of a single volume or surface phase
    /*!
     *  @param[in]        iphGlob           global phase index
     *
     *  @return                             Return a reference to the ThermoPhase object
     */
    ThermoPhase& phase(int iphGlob) const;

    //! Return the reference to the %ThermoPhase of a single volume or surface phase
    /*!
     *  @param[in]        iphGlob           global phase index
     *
     *  @return                             Return a reference to the ThermoPhase object
     */
    ThermoPhase& phase(size_t iphGlob) const;

    //! Return the common Elements object as a const pointer
    /*!
     *  @return                                  Return a pointer to a const Element object containing a description of the elements 
     *                                           in the problem.
     */
    const Elements* getGlobalElements() const;

    //! Return the element name of the eth global element
    /*!
     *  @param[in]           eGlobIndx           Element global index
     *  
     *  @return                                  String Name of the element
     */
    std::string elementName(size_t eGlobIndx) const;

    //! Index number of an element given its name
    /*!
     *    @param[in]           eName             String name of the element
     * 
     *    @return                                returns the element index of the element.
     */
    size_t elementIndex(const std::string& eName) const;

    //! Returns the number of global elements in all of the phases, combined
    /*!
     *  Note, we do not require all element objects to have the same number of elements nor have the same ordering of elements.
     *
     *  @return                                  Returns the number of distinct elements in all of the phases, combined.
     */
    size_t nElements() const;

    //! Returns the global element index given the local index
    /*!
     *  @param[in]           iphGlob             Global phase index within the %PhaseList
     *  @param[in]           eLocIndx            local element index within the phase
     *
     *  @return                                  Returns the global element index
     */
    size_t localToGlobalElementIndex(size_t iphGlob, size_t eLocIndx);

    //! Returns the local element index within the phase given the global phase index 
    /*!
     *  @param[in]           iphGlob             Global Phase Index within the %PhaseList object
     *  @param[in]           eGlobIndx           Global element index within the %PhaseList object
     *
     *  @return                                  Returns the local element index. If not present it returns npos
     */
    size_t globalToLocalElementIndex(size_t iphGlob, size_t eGlobIndx);

    //! Returns the Number of atoms of global element \a eGlob in global species \a kGlob.
    /*!
     *  @param[in]           kGlob               global species index
     *  @param[in]           eGlobIndx           global element index
     *
     *  @return                                  returns the number of atoms.
     */
    doublevalue nAtoms(const size_t kGlob, const size_t eGlobIndx) const;

    //! Return a pointer to the volume phase XML Node for a single volume phase
    /*!
     *    @param[in]         iVolIndx           Volume index of the volume phase
     *
     *    @return                                Returns a pointer to the volume phase XML Node
     */
    XML_Node* volPhaseXMLNode(size_t iVolIndx) const;

    //! Return a pointer to the surface phase XML Node for a single surface phase
    /*!
     * @param iSurIndx   Surface index of the surface phase
     *
     *   @return                  Return a pointer to the surface phase XML Node for a single surface phase.
     */
    XML_Node* surPhaseXMLNode(size_t iSurIndx) const;

    //! Return a pointer to the edge phase XML Node for a single edge phase
    /*!
     * @param iEdgeIndx   Edge index of the surface phase
     *
     *   @return                  Return a pointer to the edge phase XML Node for a single edge phase.
     */
    XML_Node* edgePhaseXMLNode(size_t iEdgeIndx) const;

    //! Returns a pointer to a single volume phase
    /*!
     *  @param iVolIndx   Volume index of the volume phase
     *
     *   @return                           Returns a reference to the %ThermoPhase object for the volume phase.
     */
    ThermoPhase& volPhase(size_t iVolIndx);

    //! Returns a reference to a single surface phase
    /*!
     *  @param[in]    iSurIndx            Surface index of the surface phase
     *
     *   @return                           Returns a reference to the %ThermoPhase object for the surface phase.
     */
    ThermoPhase& surPhase(size_t iSurIndx);

    //! Returns a reference to the ThermoPhase object of a single edge phase
    /*!
     *  @param[in]    iEdgeIndx            Surface index of the edge phase
     *
     *   @return                           Returns a reference to the %ThermoPhase object for the edge phase.
     */
    ThermoPhase& edgePhase(size_t iEdgeIndx);

    //! Return the total number of phases
    /*!
     *   @return                           returns the total number of phases of all kinds
     */
    size_t nPhases() const;

    //! Return the number of volume phases
    /*!
     *     @return                         Return the total number of volume phases
     */
    size_t nVolPhases() const;

    //! Return the number of surface phases
    /*!
     *     @return                         Return the total number of surface phases
     */
    size_t nSurPhases() const;

    //! Return the number of edge phases
    /*!
     *     @return                         Return the total number of edge phases
     */
    size_t nEdgePhases() const;

    //! Return the total number of volume species
    /*!
     *     @return                         Return the total number of volume species
     */
    size_t nVolSpecies() const;

    //! Return the total number of surface phase species
    /*!
     *     @return                         Return the total number of surface species
     */
    size_t nSurSpecies() const;

    //! Return the total number of edge phase species
    /*!
     *     @return                         Return the total number of edge species
     */
    size_t nEdgeSpecies() const;

    //! Return the total number of species in all volume, surface and edge phases
    /*!
     *     @return                          Return the total number of species in all phases
     */
    size_t nSpecies() const;

    //! Boolean indicating whether a volume phase has a kinetics object associated with it
    /*!
     *    @param[in]    iVolIndex          volume phase index 
     *
     *    @return                          Returns true if the selected volume phase has a kinetics object associated with it.
     */
    bool volPhaseHasKinetics(size_t iVolIndex) const ;

    //! Boolean indicating whether a surface phase has a surface kinetics object associated with it
    /*!
     *    @param[in]    iSurIndex          surface phase index 
     *
     *    @return                          Returns true if the selected surface phase has a kinetics object associated with it.
     */
    bool surPhaseHasKinetics(size_t iSurIndex) const ;
    
    //! Boolean indicating whether a edge phase has a edge kinetics object associated with it
    /*!
     *    @param[in]    iEdgeIndex         edge phase index 
     *
     *    @return                          Returns true if the selected edge phase has a kinetics object associated with it.
     */
    bool edgePhaseHasKinetics(size_t iEdgeIndex) const ;

    //! Return the species name given the global species index
    /*!
     *  @param[in]           iGlobSpeciesIndex   global species index within the PhaseList object
     *
     *  @return returns the species name
     */
    std::string speciesName(size_t iGlobSpeciesIndex) const;

    //! Set the state of all the phases within the PhaseList to a given temperature and pressure
    /*!
     *  Calls the underlying function \link Zuzax::ThermoPhase::setState_TP() setState_TP() \endlink
     *  for all phases in the group.
     *
     *   @param[in]          temperature         Temperature in Kelvin
     *   @param[in]          pressure            pressure in pascals
     */
    void setState_TP(doublevalue temperature, doublevalue pressure);

    //! Recalculate the maps from local to global elements
    /*!
     *  @param[in]           forceAll            Force the recalculation from scratch
     */
    void calcElementMaps(bool forceAll);

    /***********************************************************************/
    /*                BASIC INDEXING DATA                                  */
    /***********************************************************************/
protected:

    //! If true, then the phase are ordered by their dimensionality. This is the default
    /*!
     *  If this is true, then the phases are ordered by their dimensionality. All volume phases
     *  are first in PhaseList_. Then the surface phases come next. Then, the edge phases are last.
     *
     *  If false, then the PhaseList_[] order is mixed up. It's still true that the order of the volume
     *  phases within it's own list is ordered according to where they appear in the PhaseList_ []vector.
     */
    bool  m_OrderedByDimensionality;

    //! Total number of volume and surface phases defined in the object
    size_t m_NumTotPhases;

    //! Total number of species in all phases
    size_t m_NumTotSpecies;

    //!  Number of volume phases
    size_t NumVolPhases_;

    //!  Total number of volume phase species
    size_t m_totNumVolSpecies;

    //!  Vector of pointers to volume phases existing in the problem
    std::vector<ThermoPhase*> VolPhaseList;

    //! Vector of pointers to volume phase XML trees
    std::vector<XML_Node*> VolPhaseXMLNodes;

    //! Boolean vector indicating whether each volume phase has kinetics
    /*!
     *  Length is equal to the number of volume phases
     */
    std::vector<int> VolPhaseHasKinetics;

    //! Number of surface phases
    size_t m_NumSurPhases;

    //! Number of edge phases
    size_t m_NumEdgePhases;

    //! Total number of surface phase species.
    size_t m_totNumSurSpecies;

    //! Total number of edge phase species.
    size_t m_totNumEdgeSpecies;

    //! Vector of surface phases existing in the problem
    std::vector<ThermoPhase*> SurPhaseList;

    //! Pointer to the XML nodes corresponding to each surface phase
    std::vector<XML_Node*> SurPhaseXMLNodes;

    //! Boolean vector indicating whether surface phase has kinetics
    std::vector<int> SurPhaseHasKinetics;

    //! Vector of edge phases existing in the problem
    std::vector<ThermoPhase*> EdgePhaseList;

    //! Pointer to the XML nodes corresponding to each edge phase
    std::vector<XML_Node*> EdgePhaseXMLNodes;

    //! Boolean vector indicating whether edge phase has kinetics
    std::vector<int> EdgePhaseHasKinetics;

    //! Vector of phases in the problem
    std::vector<ThermoPhase*> PhaseList_;

    //! Vector of phase names in the problem
    /*!
     *  This is the string returned from the name() member function of %ThermoPhase.
     *  Note this is distinct from id(), which isn't unique when phases break apart into two.
     */
    std::vector<std::string> PhaseNames_;

    //! Number of unique elements in all of the phases.
    /*!
     * Note, we do not require all element objects have the same number of elements and the same ordering of elements.
     * 
     */
    size_t m_numElements;

    //! Index to the start of the species in a global species list
    /*!
     * length = m_NumTotPhases + 1
     */
    std::vector<size_t> m_PhaseSpeciesStartIndex;

private:

    /********************************************************** INTERNAL DATA *****************************************************/
    
    //! Boolean specifying whether this object owns its phase pointers.
    /*!
     *
     *   Note: normally, the volume and surface domains own the phase pointers. However, in thermodynamics programs,
     *         where this object is the top level organizing object,  this object will own phase pointers.
     *
     *         (ownership means that it is responsible for destroying  the object).
     */
    bool IOwnPhasePointers;

    //! Element object that represents all of the elements within all of the phases owned by this object
    Elements m_GlobalElementObj;

    //! Global Stoichiometric Coefficient array
    /*!
     *  This is a two dimensional array m_atoms(m, k). The first index is the
     *  global element index for the element in m_GlobalElementObj. The second index, k, is the global species
     *  index. The value is the number of atoms of type m in species k.
     */
    Array2D m_atoms;

    //! Global to local mapping of the Element index constraints
    /*!
     * Has a size of ( m_numElements, m_NumTotPhases).
     * The local element index in the global phase index corresponding to the global element index is given via:
     *
     *   m_globalToLocalEMap(mGlobalIndex, globPhaseIndex);
     *
     * If a phase doesn't have an element, the entry contains npos
     */ 
    Array2D m_globalToLocalEMap;

    //! Mapping between the local element index and the global element index
    /*!
     *  Outer vector is over m_NumTotPhases. The inner vector is over the number of local elements
     *
     *  The following returns the global element index:
     *
     *    m_localToGlobalEMap[globPhaseIndex][mLocalIndex];
     */
    std::vector< std::vector<size_t> > m_localToGlobalEMap;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
