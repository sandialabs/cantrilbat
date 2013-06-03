/**
 * @file SolnDomain.cpp
 *
 */

/*
 *  $Id: SolnDomain.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */
#include "SolnDomain.h"
#include "SolnLayout.h"

#include "cantera/base/ctml.h"

#include <iostream>

#ifndef safeDelete
#define safeDelete(ptr) if (ptr) { delete ptr; ptr = 0; }
#endif

using namespace std;
using namespace Cantera;

namespace m1d {
  
  //====================================================================================================================
  SolnDomain::SolnDomain() :
    DomOrder(-1)
  {
  }

  //====================================================================================================================
  SolnDomain::~SolnDomain() {
 
  }
  //====================================================================================================================
  SolnDomain::SolnDomain(const SolnDomain &r) :
    DomOrder(-1)
  {
    *this = r;
  }
  //====================================================================================================================
  SolnDomain &
  SolnDomain::operator=(const SolnDomain &r)
  {
    if (this == &r)
      return *this;

    DomOrder = r.DomOrder;

    return *this;
  }

  //====================================================================================================================
  //====================================================================================================================
  //====================================================================================================================


  SolnDomainBulk::SolnDomainBulk() :
    NumVariables(0), NumNodes(0)
  {
  }
  //====================================================================================================================
  SolnDomainBulk::SolnDomainBulk(Cantera::XML_Node & bulkXML) :
    SolnDomain(),
    NumVariables(0), 
    NumNodes(0)
  {
    readXML(bulkXML);
  }
  //====================================================================================================================
  SolnDomainBulk::~SolnDomainBulk() {
 
  }
  //====================================================================================================================
  SolnDomainBulk::SolnDomainBulk(const SolnDomainBulk &r) 
  {
    *this = r;
  }
  //====================================================================================================================
  SolnDomainBulk &
  SolnDomainBulk::operator=(const SolnDomainBulk &r)
  {
    if (this == &r) {
      return *this;
    }
    SolnDomain::operator=(r);
    NumVariables = r.NumVariables;
    NumNodes = r.NumNodes;
    return *this;
  }
  //====================================================================================================================
  void SolnDomainBulk::readXML(Cantera::XML_Node & bulkXML) {
    string sss =  bulkXML["numVariables"];
    NumVariables = Cantera::atofCheck(sss.c_str());
    sss = bulkXML["points"];
    NumNodes= Cantera::atofCheck(sss.c_str());
    
    XML_Node * bulkgXML_ptr = bulkXML.findByName("grid_data");

    XML_Node *x0XML = ctml::getByTitle(*bulkgXML_ptr, "X0");
    if (!x0XML) {
      throw CanteraError("", "");
    }
    ctml::getFloatArray(*x0XML, X0NodePos);
    int sz = X0NodePos.size();
    if (sz != NumNodes) {
      throw CanteraError("SolnDomainBulk::readXML", "sz of X0Node not right: " + int2str(sz) + "  " + int2str(NumNodes));
    }

    XML_Node *xXML = ctml::getByTitle(*bulkgXML_ptr, "X");
    if (xXML) {
      ctml::getFloatArray(*xXML, XNodePos);
      sz = XNodePos.size();
      if (sz != NumNodes) {
	throw CanteraError("SolnDomainBulk::readXML", "sz of XNode not right: " + int2str(sz) + "  " + int2str(NumNodes));
      }
    }
    
    std::vector<XML_Node *> dataList; 
    bulkgXML_ptr->getChildren("floatArray", dataList);
    int numData = dataList.size();

    int numVars = 0;
    for (int i = 0; i < numData; i++) {
      XML_Node *dataXML = dataList[i];
      string title = (*dataXML)["title"];
      string type = (*dataXML)["type"];
      if (title == "X0") continue;
      if (title == "X") continue;
      numVars++;
      VarNames.push_back(title);
      VarTypes.push_back(type);
      std::vector<double> dataValues;
      ctml::getFloatArray(*dataXML, dataValues);
      int sz = dataValues.size();
      if (sz != NumNodes) {
	throw CanteraError("SolnDomainBulk::readXML", "sz of data not right: " + int2str(sz) + "  " + int2str(NumNodes));
      }
      DataArray.push_back(dataValues);
    } 



  }
 //=====================================================================================================================
 //=====================================================================================================================
 //=====================================================================================================================
  SolnDomainSurf::SolnDomainSurf()
  {
  }
  //====================================================================================================================
  SolnDomainSurf::SolnDomainSurf(Cantera::XML_Node & surfXML) :
    SolnDomain()
  {
    readXML(surfXML);
  }
  //====================================================================================================================
  SolnDomainSurf::~SolnDomainSurf() {
 
  }
  //====================================================================================================================
  SolnDomainSurf::SolnDomainSurf(const SolnDomainSurf &r) 
  {
    *this = r;
  }
  //====================================================================================================================
  SolnDomainSurf &
  SolnDomainSurf::operator=(const SolnDomainSurf &r)
  {
    if (this == &r) {
      return *this;
    }
    SolnDomain::operator=(r);

    return *this;
  }
  //====================================================================================================================
  void SolnDomainSurf::readXML(Cantera::XML_Node & surfXML) {
    string sss =  surfXML["numVariables"];
    NumVariables = Cantera::atofCheck(sss.c_str());
    sss = surfXML["points"];
    int numNodes= Cantera::atofCheck(sss.c_str());
    if (numNodes != 1) {
      throw CanteraError("SolnDomainSurf::readXML", "nodes not right: " + int2str(numNodes));
    }



    XNodePos = ctml::getFloat(surfXML, "X");
    X0NodePos = ctml::getFloat(surfXML, "X0");
 
    std::vector<XML_Node *> dataList; 
    surfXML.getChildren("", dataList);
    int numData = dataList.size();

    int numVars = 0;
    for (int i = 0; i < numData; i++) {
      XML_Node *dataXML = dataList[i];
      string name = dataXML->name();
      if (name == "X0") continue;
      if (name == "X") continue;
      string type = (*dataXML)["type"];
      numVars++;
      VarNames.push_back(name);
      VarTypes.push_back(type);
      std::vector<double> dataValues;
      double val = ctml::getFloatCurrent(*dataXML);
      DataValues.push_back(val);
    } 
  }
  //====================================================================================================================
}
