/**
 * @file m1d_DomainLayout.cpp
 *
 */

/*
 *  $Id: SolnLayout.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */
#include "SolnLayout.h"
#include "SolnDomain.h"


#include <iostream>

#ifndef safeDelete
#define safeDelete(ptr) if (ptr) { delete ptr; ptr = 0; }
#endif

using namespace std;
using namespace Cantera;

namespace m1d
{
  


  SolnLayout::SolnLayout() :
    NumDomains(0), NumBulkDomains(0), NumSurfDomains(0), NumGbNodes(0)
  {
  }
  //====================================================================================================================
   //! Constructor with reading
  SolnLayout::SolnLayout(Cantera::XML_Node *xmlSim):
   NumDomains(0), NumBulkDomains(0), NumSurfDomains(0), NumGbNodes(0)
  {
    readXML(xmlSim);
    
  }
  //====================================================================================================================
  SolnLayout::~SolnLayout() {
    for (int ibd = 0; ibd < NumBulkDomains; ibd++) {
      safeDelete(SolnDomainBulk_List[ibd]);
    }

    for (int isd = 0; isd < NumSurfDomains; isd++) {
      safeDelete(SolnDomainSurf_List[isd]);
    }
  }
  //====================================================================================================================
  SolnLayout::SolnLayout(const SolnLayout &r) :
    NumDomains(0), NumBulkDomains(0), NumSurfDomains(0), NumGbNodes(0)
  {
    *this = r;
  }
  //====================================================================================================================
  SolnLayout &
  SolnLayout::operator=(const SolnLayout &r)
  {
    if (this == &r)
      return *this;

    for (int ibd = 0; ibd < NumBulkDomains; ibd++) {
      safeDelete(SolnDomainBulk_List[ibd]);
    }

    for (int isd = 0; isd < NumSurfDomains; isd++) {
      safeDelete(SolnDomainSurf_List[isd]);
    }

    NumDomains = r.NumDomains;
    NumBulkDomains = r.NumBulkDomains;
    NumSurfDomains = r.NumSurfDomains;
    NumGbNodes = r.NumGbNodes;

    SolnDomain_List = r.SolnDomain_List;
    SolnDomainBulk_List = r.SolnDomainBulk_List;
    SolnDomainSurf_List = r.SolnDomainSurf_List;

    for (int i = 0; i <   NumBulkDomains; i++) {
      SolnDomainBulk_List[i] = new SolnDomainBulk(*(r.SolnDomainBulk_List[i]));
      int iorder = (SolnDomainBulk_List[i])->DomOrder;
      SolnDomain_List[iorder] =  SolnDomainBulk_List[i];
    }
    for (int i = 0; i <   NumSurfDomains; i++) {
      SolnDomainSurf_List[i] = new SolnDomainSurf(*(r.SolnDomainSurf_List[i]));
      int iorder = (SolnDomainSurf_List[i])->DomOrder;
      SolnDomain_List[iorder] =  SolnDomainSurf_List[i];
    }

    return *this;
  }


  //====================================================================================================================
  void
  SolnLayout::readXML(Cantera::XML_Node *simXML) {
    std::vector<Cantera::XML_Node*> ccc;
    simXML->getChildren("domain", ccc);

    NumDomains = ccc.size();
    NumBulkDomains = 0;
    NumSurfDomains = 0;
    NumGbNodes = 0;

    for (int iDom = 0; iDom < NumDomains; iDom++) {
      XML_Node &domXML = *(ccc[iDom]);
      std::string ttt = domXML["type"];
      if (ttt == "bulk") {
	SolnDomainBulk * bulkD = readBulkDomainXML(domXML);
	SolnDomainBulk_List.push_back(bulkD);
	NumBulkDomains++;
	bulkD->DomOrder = iDom;
      } else if (ttt == "surface") {
	SolnDomainSurf * surfD = readSurfDomainXML(domXML);
	SolnDomainSurf_List.push_back(surfD);
	NumSurfDomains++;
	surfD->DomOrder = iDom;
      } else {
	throw Cantera::CanteraError("SolnLayout::readXML", "error");
      }
      
    }
  }
 //=====================================================================================================================
  SolnDomainBulk * SolnLayout::readBulkDomainXML(Cantera::XML_Node & domXML) {
    SolnDomainBulk *bulkD = new SolnDomainBulk(domXML);


    return bulkD;
  }
  //=====================================================================================================================
  SolnDomainSurf * SolnLayout::readSurfDomainXML(Cantera::XML_Node & domXML) {
   SolnDomainSurf *surfD = new SolnDomainSurf();

  

    return surfD;

  }
  //====================================================================================================================
}
//======================================================================================================================
