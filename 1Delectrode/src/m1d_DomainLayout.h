/**
 * @file m1d_DomainLayout.h  This is a top level object that describes
 *         the domain layout of the problem.
 */

/*
 *  $Id: m1d_DomainLayout.h 504 2013-01-07 22:32:48Z hkmoffa $
 */
#ifndef M1D_DOMAINLAYOUT_H
#define M1D_DOMAINLAYOUT_H
/*
 *
 *
 *     Formulation of the structure of the problem.
 *
 *
 *    We will create a system for solving one dimensional problems based
 *    on using Trilinos and Cantera.
 *
 *    The matrix stencil for the calculations will be based on  a set of
 *    of domains that are connected with each other by tie-regions.
 *    The tie-regions may have a set of equations associated with
 *    the tie region. This set of equations will have dependencies
 *    on the neighboring nodes of the bounding domains.
 *
 *
 *
 *    |  Tie        |  Domain  | Tie      |  Domain  |  Tie     |
 *    |  Region     |  Zero    | Region   |  One     |  Region  |
 *    |  0          |          | 1        |          |  2       |
 *    |             |          |          |          |          |
 *    |  Ne^T_0     |  Ne^D_0  | Ne^T_1   |  Ne^D_1  |  Ne_T_2  |
 *
 *
 *
 *
 *   The system will be able to do solved using a MP formulation. However,
 *   we will not build into the formulation any kind of extensibility. All
 *   nodes will know everything about the problem. All nodes will have all
 *   data. MP will be achieved by breaking the nodes of the problem up
 *   into sections. The residual for each section will belong to one processor
 *   The matrix row for that section will belong to one processor.
 *   Matrix solves will occur via a MP method.
 *
 *
 *
 *  For our first step we will assume one domain, with no tie regions !
 */

#include <cantera/equil/MultiPhase.h>

#include <string>
#include <vector>
#include <map>

namespace m1d
{


class GlobalIndices;
class ProblemResidEval;
class ProblemStatement;
class BulkDomainDescription;
class SurfDomainDescription;
class BulkDomain1D;
class SurDomain1D;

//! DomainLayout is a light class that describes the overall
//! disposition of the domains in the problem
/*!
 *
 */
class DomainLayout
{
public:

  //! Constructor
  DomainLayout(ProblemStatement* psInput_ptr = 0);

  //! Constructor
  DomainLayout(std::vector<std::string> domainList, ProblemStatement *psInput_ptr,
	       std::map<std::string, Cantera::MultiPhase*> bulkMap,
	       std::map<std::string, int> surfaceMap);

  //! Copy Constructor
  /*!
   *
   * @param r   Object to be copied
   */
  DomainLayout(const DomainLayout &r);

  //! destructor
  virtual
  ~DomainLayout();

  DomainLayout &
  operator=(const DomainLayout &r);

  //! Overarching structure for getting DomainLayout and
  //! DomainDescription Objects initialized and filled.
  /*!
   *
   */
  virtual void
  InitializeDomainPicture();

  //! Allocate the domain structure
  /*!
   *   Note since we don't know what will go here, this is realy
   *   a placeholder for child routines that will have
   *   the required information to execute the routine.
   */
  virtual void
  malloc_domains();

  //! Add a bulk domain with its nodes to the right end
  /*!
   *
   * @param bdd
   * @param numNodes
   * @param startZ
   * @param endZ
   */
  virtual void
  addBulkDomainToRightEnd(BulkDomainDescription *bdd,
                          int numNodes,
                          double startZ,
                          double endZ);

  //! Add a surface domain to the left end
  /*!
   *  @param  sddL  Left Surface domain description
   *   @param bdd    bulk domain description
   */
  void
  addSurfDomainToLeftEnd(SurfDomainDescription *sddL,
                         BulkDomainDescription *bdd);

  //! Add a surface domain to the right end
  /*!
   *  @param  sddR  Left Surface domain description
   *  @param  bdd   Bulk domain description
   */
  void
  addSurfDomainToRightEnd(SurfDomainDescription *sddR,
                          BulkDomainDescription *bdd);

  //! Tell the domain descriptions what they need to know
  void
  InfoToDomainDescriptions();

  //! Update the domain layout information to underlying domains
  /*!
   *   Note, that this routine is now basically an internal check routine.
   *   The functionality has been placed elsewhere
   */
  void
  updateAdjInfoInDomainDescriptions();

  //! Update the global node layout information to
  //! underlying domains
  void
  updateGbNodeInfoInDomainDescriptions();

  //! Tell the domain descriptions what their domain boundary positions are
  void
  updateXposInDomainDescriptions();

  //! Set the equation and variables lists
  void
  SetEquationsVariablesList();

  //! Set the equation descriptions for the domain layout
  void
  SetEqnDescriptions();

  void 
  setProblemResid(ProblemResidEval *problemResid_ptr);

  virtual void
  generateDomain1D(ProblemResidEval *problemResid_ptr);

  void
  InitializeXposNodes(GlobalIndices *gi_ptr);

  //! Set the temperature and pressure reference values for all domains in the problem
  //! statement
  /*!
   *  The temperature and pressure is passed down to all Domain1D objects within the
   *  DomainLayout
   *
   * @param temp_ref value of the temperature (Kelvin)
   * @param pres_ref value of the pressure  (Pascal)
   */
  void set_TP_Reference(const double temp_ref, const double pres_ref);

  //===============================================================================================================

  //! Total number of domains, both surface and bulk
  /*!
   * This number is equal to the sum of the number of bulk
   * domains and the number of surface domains.
   */
  int NumDomains;

  //! Number of Bulk domains.
  /*!
   *   These bulk domains are distributed across the mesh.
   *   There are one or more bulk domains in each problem
   */
  int NumBulkDomains;

  //! Number of surface domains
  /*!
   *   These surface domains are not distributed across the mesh.
   *   They are defined to exist at one mesh point. They are
   *   defined to exist at the top and bottom node, always, and
   *   there exists one surface domain between neighboring bulk
   *   domains.
   *   There are 2 or more surface domains in every problem.
   */
  int NumSurfDomains;

  //! Total number of global nodes
  int NumGbNodes;

  //! Vector of bulk domain Descriptions
  std::vector<BulkDomainDescription *> BulkDomainDesc_global;

  //! Vector of bulk domains in the problem
  /*!
   *   The index is from right to left within the domain
   */
  std::vector<BulkDomain1D *> BulkDomain1D_List;

  //! Vector of Surface Domain Descriptions
  std::vector<SurfDomainDescription *> SurfDomainDesc_global;

  //! Vector of surface domains in the problem
  /*!
   *   The index is from right to left within the domain
   */
  std::vector<SurDomain1D *> SurDomain1D_List;

  //! Starting global node for the bulk domain regions
  /*!
   *    Length = number of bulk domains
   */
  std::vector<int> StartGBNode_Domain;

  //! Ending node for the domain regions
  /*
   *  Length = number of bulk domains
   */
  std::vector<int> EndGBNode_Domain;

  //! Location node for the surface domain regions
  /*!
   * These are defined to occur at one global node.
   *
   *  Length = number of surface domain regions.
   */
  std::vector<int> locGBNode_SurfDomain_List;

  //! Starting left location of each of the bulk domains
  /*!
   *  Length is the number of bulk domains
   */
  std::vector<double> StartXLoc_Domain;

  //! Starting right location of each of the bulk domains
  /*!
   *  Length is the number of bulk domains
   */
  std::vector<double> EndXLoc_Domain;

  //! Left boundary of the first bulk domain
  doublereal XLoc_LeftBoundary;

  //! Right boundary of the last domain
  doublereal XLoc_RightBoundary;

  //! pointer to the Problem residual 
  ProblemResidEval *problemResid_;

  //! Pointer to the input file 
  ProblemStatement *psInput_ptr_;

 protected:

  //! Ordered list of names for domains
  std::vector<std::string> domainList_;

  //! Map between bulk domain names and 
  //!the corresponding Cantera MultiPhase objects
  /*!
   * HKM 1/3/2009
   * suggest moving this to its own layer. 
   * Also, MultiPhase should be replaced with another object
   * that more fully represents what's going on. Put it within Cantera
   * in a new porousMultiphase directory.
   */
  //std::map<std::string, Cantera::MultiPhase*> bulkMap_;

  //! Map between surface domain names and 
  //! the corresponding reacting electrode objects
  /*!
   * HKM 1/3/2009
   * suggest moving this to its own layer. 
   * Put electrode objects within Cantera in a new directory.
   */
  //std::map<std::string, int> surfaceMap_ ;

};
//=====================================================================================
//! Each processor contains a pointer to one global instance of this class.
//extern DomainLayout *DL_Global_ptr;
//=====================================================================================


// We will include several simple types of domain layouts within this file

//! Hard coded simple diffusion problem
/*!
 *
 */
class SimpleDiffusionLayout : public DomainLayout
{
public:
  //! Constructor
  SimpleDiffusionLayout(ProblemStatement *psInput_ptr = 0);

  //! Constructor - hard coded problem
  SimpleDiffusionLayout(int probNum, ProblemStatement *psInput_ptr);

  //! Copy constructor
  /*!
   * @param r  Object to be copied
   */
  SimpleDiffusionLayout(const SimpleDiffusionLayout &r);

  //! Destructor
  virtual
  ~SimpleDiffusionLayout();

  //! Assignment operator
  /*!
   * @param r  Object to be copied.
   * @return   returns a changeable reference to the current object
   */
  SimpleDiffusionLayout &
  operator=(const SimpleDiffusionLayout &r);

  //! Allocate the domain structure
  virtual void
  malloc_domains();


  //==================================================================================
  //! Hard coded problem number
  int ProbNum_;

};


//! Hard coded simple time dependent diffusion problem
/*!
 *
 */
class SimpleTimeDependentDiffusionLayout : public DomainLayout
{
public:
  //! Constructor
  SimpleTimeDependentDiffusionLayout(ProblemStatement *psInput_ptr = 0);

  //! Constructor - hard coded problem
  SimpleTimeDependentDiffusionLayout(int probNum, ProblemStatement *psInput_ptr);

  //! Copy constructor
  /*!
   * @param r  Object to be copied
   */
  SimpleTimeDependentDiffusionLayout(const SimpleTimeDependentDiffusionLayout &r);

  //! Destructor
  virtual
  ~SimpleTimeDependentDiffusionLayout();

  //! Assignment operator
  /*!
   * @param r  Object to be copied.
   * @return   returns a changeable reference to the current object
   */
  SimpleTimeDependentDiffusionLayout &
  operator=(const SimpleTimeDependentDiffusionLayout &r);

  //! Allocate the domain structure
  virtual void
  malloc_domains();


  //==================================================================================
  //! Hard coded problem number
  int ProbNum_;

};
//=====================================================================================
} // end namespace m1d
//=====================================================================================
#endif
