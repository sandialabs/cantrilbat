/**
 * @file m1d_NodalVars.h
 *
 */

/*
 *  $Id: m1d_NodalVars.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
#ifndef M1D_NODALVARS_H
#define M1D_NODALVARS_H

#include "m1d_Eqn_Names.h"
#include "m1d_EqnVarTypes.h"

#include <vector>

#include <map>

namespace m1d
{

class DomainLayout;
class LocalNodeIndices;

//! Class that describes the environment at one node in the domain
/*!
 *
 *   The ordering of the equation at a single node is as follows.
 *
 *      domain_Left equations
 *
 *      domain_Right equations (if not mapped into domain_left equations)
 *
 *      domain_other equations*     (* not yet implemented)
 *
 *      surfaceDomain_0 equations
 *      surfaceDomain_1 equations
 *       . . .
 *      surfaceDomain_N-1 equations
 *
 *   The surfaceDomain_0 object determines what equations from the domain_Right
 *   bulk domain gets mapped into the domain_Left equations.
 *
 *   When a bulk equation gets mapped from one bulk domain to another, it means
 *   that there will be created a single conservation equation for the unknown
 *   striding both domains, and the variable is continuous across the interface
 *   as well.
 */
class NodalVars
{
public:

    //! Default constructor
    NodalVars(DomainLayout* dl_ptr);

    //! Constructor
    NodalVars(const int gbnode, DomainLayout* dl_ptr);

    //! Copy Constructor
    /*!
     *   @param r  object to be copied.
     */
    NodalVars(const NodalVars& r);

    //! Destructor
    ~NodalVars();

    //! Assignment Operator
    /*!
     *   @param r  object to be copied.
     */
    NodalVars& operator=(const NodalVars& r);

    //! This routine will search the global information to discover what
    //! domains are located at this node.
    /*!
     *
     */
    void DiscoverDomainsAtThisNode();

    //! This routine discovers and orders the equation unknowns at a node
    /*!
     * This is the ONE PLACE IN THE CODE that discovers and orders the
     * equations at a node.
     */
    void GenerateEqnOrder();

    //! Generate a name for the kth variable at the node
    /*!
     *
     * @param k  input the index of the variable on the local node
     * @return  Returns a string representing the variable.
     */
    std::string VariableName(int k);

    //! Change the node position
    /*!
     * @param xNodePos  new value of the node position
     */
    void changeNodePosition(double xNodePos);

    //! Set up the initial node positions
    /*!
     * @param xNodePos
     * @param x0NodePos
     * @param xFracNodePos
     */
    void setupInitialNodePosition(double x0NodePos, double xFracNodePos);

    //! Returns the bulk domain index given the bulk domain ID
    /*!
     *  @param[in]           myBBD_ID              Bulk domain ID
     *
     *  @return                                    Returns the index of the bulk domain
     *                                             If the bulk domain isn't on this node, an npos is returned.
     */
    size_t bindexBulkDomain_fromID(int myBDD_ID);

    //!  Return the starting index of a particular variable from the start of the nodal solution vector
    /*!
     *  These functions are important as they are used frequently to index into the
     *  solution vector.
     *
     *    @param variableType   VAR_TYPE to look up the index for
     *    @param subVarIndex    VAR_TYPE_SUBNUM subindex to look up the variable
     *
     *    @return  Returns npos if there isn't a variable of that type and subtype.
     *             Returns the index into the solution vector from the start of
     *             variables at that node.
     */
    size_t indexBulkDomainVar(VAR_TYPE variableType, VAR_TYPE_SUBNUM subVarIndex) const;

    size_t indexBulkDomainVar(const VarType& vt) const;

    //!  Find the offset index of a particular variable into the nodal solution vector
    /*!
     *   These functions are important as they are used frequently to index into the
     *   solution vector.
     *
     *    @param variableTypeS   VAR_TYPE size_t value to look up the index
     *
     *    @return  Returns npos=size_t(-1) if there isn't a variable of that type.
     *             Returns the index into the solution vector from the start of
     *             variables of that VAR_TYPE at that node.
     */
    inline size_t indexBulkDomainVar0(size_t variableTypeS) const;

    //!   Find the index of a particular bulk domain variable, given by its index, at this particular node
    /*!
     *   This function isn't important for speed
     *
     *    @param  indexBD   Index of the variable within the list of unknowns on the bulk domain
     *
     *    @return BDnum     Number of the bulk domain at the current node. Defaults to zero,
     *                      which is defined as the first bulk domain (MAY CHANGE)
     *
     *    @return           Returns the offset of the variable index at the current node.
     */
    size_t indexBulkDomainVar_BDIndex(size_t indexBDom, size_t BDnum = 0) const;


    size_t indexBulkDomainEqn(EQ_TYPE equationType, EQ_TYPE_SUBNUM subEqnIndex) const;

    //!  Find the offset index of a particular equation into the nodal solution vector
    /*!
     *   These functions are important as they are used frequently to index into the
     *   solution vector.
     *
     *    @param equationTypeS   EQN_TYPE size_t value to look up the index
     *
     *    @return  Returns npos=size_t(-1) if there isn't an equation of that type.
     *             Returns the index into the solution vector from the start of
     *             equations of that EQN_TYPE at that node.
     */
    inline size_t indexBulkDomainEqn0(size_t equationTypeS) const;

    inline size_t numberBulkDomainVar0(size_t variableTypeS) const;

    inline size_t numberBulkDomainEqn0(size_t equationTypeS) const;

    //! Returns the current node position
    /*!
     *  @return                                  Returns the current node position of the node
     */
    double xNodePos() const;

    //! Returns the initial node position
    /*!
     *  @return                                  Returns the original node position of the node
     */
    double x0NodePos() const;

    //! Return the fraction node position from the left boundary
    /*!
     *  @return                                  Returns the fractional node position of the node
     */
    double xFracNodePos() const;

    //  ------------------------------------------------------------------------------------------------------------

    //! Global node value
    int GbNode;

    //! What local node am I?
    /*!
     *  If I am not a local node on this processor, set this variable to -1.
     */
    int LcNode;

    //! Number of equations located at this node
    int NumEquations;

    //! Number of bulk domains located at the node
    int NumBulkDomains;

    //! Number of surface domains located at the node
    int NumSurfDomains;

    //! Starting global index for equations located at this node
    /*!
     * This is in terms of the global equation index, GbEqnIndex,
     * and is independent of processor number
     */
    int EqnStart_GbEqnIndex;

    //! Global index of the bulk domain at this node given its bulk domain number at this node
    /*!
     *  Length = number of bulk domains, NumBulkDomains, at the current node
     *
     *  Value = index of the Bulk Domain in the global vector of bulk domains
     */
    std::vector<int> BulkDomainIndex_BDN;

    //! Offset index for the equations corresponding
    //! to each domain located at this node
    /*!
     *  Length = number of domains
     *
     *
     *    @deprecated  This doesn't make any sense to have. The bulk equations
     *                 are intermixed between domains, and are not contiguous
     *                 within each bulk domain.
     */
    std::vector<int> OffsetIndex_BulkDomainEqnStart_BDN;

    //! Map from the bulk domain id to the order index at the current node
    /*!
     *   Given the bulk domain ID, this returns the local index of the bulk domain
     *   at the current node. If the bulk domain is not at the current node,
     *   this returns -1.
     */
    std::map<int, int> BulkDomainIndex_fromID;


    //! Listing of the surface domains at this node
    /*!
     *  Length = number of domains, NumDomains
     */
    std::vector<int> SurfDomainIndex_SDN;

    //! Offset index for the equations corresponding
    //! to each domain located at this node
    /*!
     *  Length = number of domains
     *
     *    This makes sense still to have because the surface unknowns
     *    if there are any, are contiguous.
     */
    std::vector<int> OffsetIndex_SurfDomainEqnStart_SDN;


    //! Map from the surface domain id to the order index at the current node
    /*!
     *  If the domain doesn't exist on this node, a -1 is returned.
     */
    std::map<int, int> SurfDomainIndex_fromID;


    //! Listing of the Variable Type for each dof at a node
    /*!
     *   This is a listing of the variables types, with the index
     *   being the equation number
     *
     *  Length = number of variables at the node
     */
    std::vector<VarType> VariableNameList_EqnNum;

    //! Listing of the variable subtype
    std::vector<VAR_TYPE_SUBNUM> VariableSubType_EqnNum;

    //! Listing of the Equation Types for each dof at a node
    /*!
     *   This is a listing of the variables types, with the index
     *   being the equation number
     *
     *  Length = number of variables at the node
     */
    std::vector<EqnType> EquationNameList_EqnNum;

    //! Listing of the variable subtype
    std::vector<VAR_TYPE_SUBNUM> EquationSubType_EqnNum;

    //! Map between the variable type and offset of the variable in the unknowns for the node
    /*!
     *       Offset_VarType[MoleFraction_Species] is the offset for the mole fraction variables.
     */
    std::map<VAR_TYPE, size_t> Offset_VarType;

    std::map<EQ_TYPE, size_t> Offset_EqnType;


    //! Vector between the variable type redefined as a size_t and the offset of the variable in the unknowns for the node
    /*!
     *       Offset_VarType[MoleFraction_Species] is the offset for the mole fraction variables.
     */
    std::vector<size_t> Offset_VarTypeVector;
    std::vector<size_t> Offset_EqnTypeVector;


    //! Map between the variable type and the number of variables of that VAR_TYPE in the unknowns for the node
    /*!
     *   Example:
     *      Number_VarType[MoleFraction_Species] is the number of mole fraction variables. Each will of course
     *      have different subtypes.
     */
    std::map<VAR_TYPE, size_t> Number_VarType;
    std::map<EQ_TYPE, size_t> Number_EqnType;

    //! Vector between the variable type redefined as a size_t and the number of variables of that
    //! VAR_TYPE in the unknowns for the node
    /*!
     *       Number_VarType[MoleFraction_Species] is the number of mole fraction variables.
     */
    std::vector<size_t> Number_VarTypeVector;
    std::vector<size_t> Number_EqnTypeVector;

protected:
    //! Current Spatial position of the node
    /*!
     *
     */
    double XNodePos;

    //! Initial Spatial position of the node
    /*!
     *
     */
    double X0NodePos;

    //! Fraction of the nodal position from the left to the right
    double XFracNodePos;

public:
    //! Owning Domain Layout for the node
    DomainLayout* DL_ptr_;

};

//====================================================================================================
inline size_t NodalVars::indexBulkDomainVar0(size_t variableTypeS) const
{
    return Offset_VarTypeVector[variableTypeS];
}
//====================================================================================================
inline size_t NodalVars::indexBulkDomainEqn0(size_t equationTypeS) const
{
    return Offset_EqnTypeVector[equationTypeS];
}
//====================================================================================================
inline size_t NodalVars::numberBulkDomainVar0(size_t variableTypeS) const
{
    return Number_VarTypeVector[variableTypeS];
}
//====================================================================================================
inline size_t NodalVars::numberBulkDomainEqn0(size_t equationTypeS) const
{
    return Number_EqnTypeVector[equationTypeS];
}
//====================================================================================================
}

#endif
