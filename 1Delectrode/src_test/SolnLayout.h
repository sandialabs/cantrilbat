/**
 * @file m1d_DomainLayout.h  This is a top level object that describes
 *         the domain layout of the problem.
 */

#ifndef M1D_SOLNLAYOUT_H
#define M1D_SOLNLAYOUT_H


#include <vector>
#include "cantera/base/global.h"

#ifdef useZuzaxNamespace
#ifndef ZZCantera
#define ZZCantera Zuzax
#endif
#else
#ifndef ZZCantera
#define ZZCantera Cantera
#endif
#endif


namespace m1d
{

  class SolnDomain;
  class SolnDomainBulk;
  class SolnDomainSurf;

  class SolnLayout {
  public:

    //! Constructor
    SolnLayout();

    //! Constructor with reading
    SolnLayout(ZZCantera::XML_Node *xmlSim);

    //! Copy Constructor
    /*!
     *
     * @param r   Object to be copied
     */
    SolnLayout(const SolnLayout &r);

    //! destructor
    virtual
    ~SolnLayout();

    SolnLayout &
    operator=(const SolnLayout &r);

    void readXML(ZZCantera::XML_Node *simXML_ptr);

    SolnDomainBulk* readBulkDomainXML(ZZCantera::XML_Node & domXML);

    SolnDomainSurf* readSurfDomainXML(ZZCantera::XML_Node & domXML);

    //! Number of  domains.
    /*!
     *   These are bulk and surface domains
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

    //! Vector of domain Descriptions
    std::vector<SolnDomain *> SolnDomain_List;

    //! Vector of bulk domain Descriptions
    std::vector<SolnDomainBulk *> SolnDomainBulk_List;

    //! Vector of surf domain Descriptions
    std::vector<SolnDomainSurf *> SolnDomainSurf_List;

  };

  //=====================================================================================
} // end namespace m1d
//=====================================================================================
#endif
