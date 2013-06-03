/**
 * @file m1d_DomainLayout.h  This is a top level object that describes
 *         the domain layout of the problem.
 */

#ifndef M1D_SOLNDOMAIN_H
#define M1D_SOLNDOMAIN_H

#include <vector>
#include "cantera/base/global.h"

namespace m1d
{

  class SolnDomain;
  class SolnDomainBulk;
  class SolnDomainSurf;
  //====================================================================================================================
  class SolnDomain {
  public:

    //! Constructor
    SolnDomain();

    //! Copy Constructor
    /*!
     *
     * @param r   Object to be copied
     */
    SolnDomain(const SolnDomain &r);

    //! destructor
    virtual  ~SolnDomain();

    SolnDomain & operator=(const SolnDomain &r);


    //! Order of the domain in the list
    /*!
     *  
     */
    int DomOrder;

  };

  //====================================================================================================================
  class SolnDomainBulk : public SolnDomain {
  public:

    SolnDomainBulk();
    SolnDomainBulk(Cantera::XML_Node & bulkXML);
    //! Copy Constructor
    /*!
     *
     * @param r   Object to be copied
     */
    SolnDomainBulk(const SolnDomainBulk &r);

    //! destructor
    virtual
    ~SolnDomainBulk();

    SolnDomainBulk &  operator=(const SolnDomainBulk &r);

    void readXML(Cantera::XML_Node & bulkXML);

    int NumVariables;
    int NumNodes;
    
    std::vector<double> XNodePos;
    std::vector<double> X0NodePos;

    std::vector<std::string> VarNames;

    std::vector<std::string> VarTypes;

    std::vector<std::vector<double> > DataArray;

  };

  //====================================================================================================================
  class SolnDomainSurf : public SolnDomain {
  public:

    SolnDomainSurf();
    SolnDomainSurf(Cantera::XML_Node & surfXML);

    //! Copy Constructor
    /*!
     *
     * @param r   Object to be copied
     */
    SolnDomainSurf(const SolnDomainSurf &r);

    //! destructor
    virtual
    ~SolnDomainSurf();

    SolnDomainSurf &  operator=(const SolnDomainSurf &r);

    void readXML(Cantera::XML_Node & surfXML);

    int NumVariables;

    double XNodePos;
    double X0NodePos;

    std::vector<std::string> VarNames;
    std::vector<std::string> VarTypes;
    
    std::vector<double> DataValues;

  };
  //====================================================================================================================


} // end namespace m1d
//======================================================================================================================
#endif
