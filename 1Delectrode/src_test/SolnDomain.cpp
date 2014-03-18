/**
 * @file SolnDomain.cpp
 *
 */

/*
 *  $Id: SolnDomain.cpp 506 2013-01-07 22:43:59Z hkmoffa $
 */
#include "SolnDomain.h"
#include "SolnLayout.h"

#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"

#include <iostream>
#include <fstream>

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
    NumVariables = Cantera::fpValueCheck(sss);
    sss = bulkXML["points"];
    NumNodes= Cantera::fpValueCheck(sss);
    
    XML_Node * bulkgXML_ptr = bulkXML.findByName("grid_data");
    if (!bulkgXML_ptr) {
	throw CanteraError(" SolnDomainBulk::readXML", "Can't find grid_data");
    }

 
    ctml::getFloatArray(*bulkgXML_ptr, X0NodePos, true, "", "X0");
    int sz = X0NodePos.size();
    if (sz != NumNodes) {
      throw CanteraError("SolnDomainBulk::readXML()", "sz of X0Node not right: " + Cantera::int2str(sz) + "  " + Cantera::int2str(NumNodes));
    }

    XML_Node *xXML = bulkgXML_ptr->findByName("X");
    if (xXML) {
      ctml::getFloatArray(*xXML, XNodePos);
      sz = XNodePos.size();
      if (sz != NumNodes) {
	throw CanteraError("SolnDomainBulk::readXML", "sz of XNode not right: " + int2str(sz) + "  " + int2str(NumNodes));
      }
    }


    const std::vector<XML_Node *>& gridChildren = bulkgXML_ptr->children();
    std::vector<XML_Node *> dataList; 
    for (size_t iC = 0; iC < gridChildren.size(); iC++) {
       XML_Node* xx = gridChildren[iC]->findByAttr("vtype", "floatArray");
       if (xx) {
         dataList.push_back(xx);
       }
    }
    
    //bulkgXML_ptr->getChildren("floatArray", dataList);
    size_t numData = dataList.size();

    int numVars = 0;
    for (size_t i = 0; i < numData; i++) {
      XML_Node *dataXML = dataList[i];
      string vtitle = (*dataXML).name();
      string vtype = (*dataXML)["type"];
      if (vtitle == "X0") continue;
      if (vtitle == "X") continue;
      numVars++;
      VarNames.push_back(vtitle);
      VarTypes.push_back(vtype);
      std::vector<double> dataValues;
      size_t sz = ctml::getFloatArray(*dataXML, dataValues, false, "", vtitle);
      if (sz != (size_t) NumNodes) {
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
    NumVariables = Cantera::fpValueCheck(sss);
    sss = surfXML["points"];
    int numNodes= Cantera::fpValueCheck(sss);
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
