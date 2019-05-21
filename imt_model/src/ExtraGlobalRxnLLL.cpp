/**
 *  @file example2.cpp
 *
 */
/*
 * $Id: ExtraGlobalRxn.cpp 507 2013-01-07 22:48:29Z hkmoffa $
 * 
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

//  Example 2
//
//  Read a mechanism, and print to the standard output stream a
//  well-formatted Chemkin ELEMENT section.
//

#include "ExtraGlobalRxn.h"
//#include "IdealReactingGas.h"

//#include "TemperatureTable.h"

#include "LE_PickList.h"
#include "BE_MoleComp.h"

#include "mdp_allo.h"

//#include "IdealSolidSolnPhase.h"


#include "zuzax/numerics/DenseMatrix.h"

// Kinetics includes
#include "zuzax/kinetics.h"
#include "zuzax/kinetics/InterfaceKinetics.h"
#include "zuzax/thermo/SurfPhase.h"
//#include "SolidKinetics.h"
#include "zuzax/kinetics/KineticsFactory.h"

#include <iostream>
#include <new>
#include <string>
#include <vector>
#include <typeinfo>

//#include "electrodeCell_kin.h"


using namespace std;
using namespace Zuzax;

#define MIN(x,y)     (( (x) < (y) ) ? (x) : (y))

static void erase_vd( std::vector<doublereal> & m_vec, int index) {
  std::vector<double>::iterator ipos;
  ipos = m_vec.begin();
  ipos += index;
  m_vec.erase(ipos);
}

static void erase_vi( std::vector<int> & m_vec, int index) {
  std::vector<int>::iterator ipos;
  ipos = m_vec.begin();
  ipos += index;
  m_vec.erase(ipos);
}

static void addV(int kkinspec, double ps, std::vector<int> & m_Products, 
	    std::vector<doublereal> & m_ProductStoich) {
  int nsize = m_Products.size();
  for (int i = 0; i < nsize; i++) {
    if (m_Products[i] == kkinspec) {
      m_ProductStoich[i] += ps;
      return;
    }
  }
  m_Products.push_back(kkinspec);
  m_ProductStoich.push_back(ps);
}

ExtraGlobalRxn::ExtraGlobalRxn(Kinetics *k_ptr) :
  m_ThisIsASurfaceRxn(false),
  m_kinetics(k_ptr),
  m_InterfaceKinetics(0),
  m_nKinSpecies(0),
  m_nReactants(0),
  m_nProducts(0),
  m_nNetSpecies(0),
  m_nRxns(0),
  m_SpecialSpecies(-1),
  m_SpecialSpeciesProduct(true),
  iphaseKin(0),
  m_ok(false),
  m_reversible(true)
{

  m_InterfaceKinetics = dynamic_cast<InterfaceKinetics *>(k_ptr);
  if  (m_InterfaceKinetics) {
    m_ThisIsASurfaceRxn = true;
  } 
  m_nRxns = m_kinetics->nReactions();
  m_ElemRxnVector.resize(m_nRxns,0.0);
  m_nKinSpecies = m_kinetics->nTotalSpecies();
}

ExtraGlobalRxn::~ExtraGlobalRxn() {
}
void ExtraGlobalRxn::setupElemRxnVector(double *RxnVector, 
					int specialSpecies) {
  int i;
  int kkinspec;
  for (int i = 0; i < m_nRxns; i++) {
    m_ElemRxnVector[i] = RxnVector[i];
  }
  for (int i = 0; i < m_nRxns; i++) {
    if (RxnVector[i] > 0.0) {
      for (kkinspec = 0; kkinspec < m_nKinSpecies; kkinspec++) {
	double ps = m_kinetics->productStoichCoeff(kkinspec, i);
	if (ps > 0.0) {
	  addV(kkinspec, RxnVector[i]* ps, m_Products, m_ProductStoich);
	  addV(kkinspec, RxnVector[i]* ps, m_NetSpecies, m_netStoich);
	}
	double rs = m_kinetics->reactantStoichCoeff(kkinspec, i);
	if (rs > 0.0) {
	  addV(kkinspec,  RxnVector[i] * rs, m_Reactants, m_ReactantStoich);
	  addV(kkinspec, -RxnVector[i] * rs, m_NetSpecies, m_netStoich);
	}
      }
    } else if  (RxnVector[i] < 0.0) {
      for (kkinspec = 0; kkinspec < m_nKinSpecies; kkinspec++) {
	double ps = m_kinetics->productStoichCoeff(kkinspec, i);
	if (ps > 0.0) {
	  addV(kkinspec,- RxnVector[i]* ps, m_Reactants, m_ReactantStoich);
	  addV(kkinspec, RxnVector[i]* ps, m_NetSpecies, m_netStoich);
	}
	double rs = m_kinetics->reactantStoichCoeff(kkinspec, i);
	if (rs > 0.0) {
	  addV(kkinspec, -RxnVector[i] * rs, m_Products, m_ProductStoich);
	  addV(kkinspec, -RxnVector[i] * rs, m_NetSpecies, m_netStoich);
	}
      }
    }
  }
 Recheck:
  for (i = 0; i < static_cast<int>(m_Products.size()); i++) {
    if (m_ProductStoich[i] == 0.0) {
      erase_vi(m_Products, i);
      erase_vd(m_ProductStoich, i);
      goto Recheck ;
    }
  }
  for (i = 0; i < static_cast<int>(m_Reactants.size()); i++) {
    if (m_ReactantStoich[i] == 0.0) {
      erase_vi(m_Reactants, i);
      erase_vd(m_ReactantStoich, i);
      goto Recheck ;
    }
  }
 for (i = 0; i < static_cast<int>(m_NetSpecies.size()); i++) {
    if (m_netStoich[i] == 0.0) {
      erase_vi(m_NetSpecies, i);
      erase_vd(m_netStoich, i);
      goto Recheck ;
    }
  }

  for (i = 0; i < static_cast<int>(m_Products.size()); i++) {
    int ik = m_Products[i];
    for (int j = 0; j < static_cast<int>(m_Reactants.size()); j++) {
      int jk = m_Reactants[j];
      if (ik == jk) {
	if (m_ProductStoich[i] == m_ReactantStoich[j]) {
	  erase_vi(m_Products, i);
	  erase_vd(m_ProductStoich, i);
	  erase_vi(m_Reactants, j);
	  erase_vd(m_ReactantStoich, j);
	} else if (m_ProductStoich[i] > m_ReactantStoich[j]) {
	  m_ProductStoich[i] -= m_ReactantStoich[j];
	  erase_vi(m_Reactants, j);
	  erase_vd(m_ReactantStoich, j);
	} else {
	  m_ReactantStoich[j] -=   m_ProductStoich[i];
	  erase_vi(m_Products, i);
	  erase_vd(m_ProductStoich, i);
	}
	// We just screwed up the indexing -> restart.
	goto Recheck ;
      }

    }
  }
  m_nProducts = m_Products.size();
  m_nReactants = m_Reactants.size();
  m_nNetSpecies = m_NetSpecies.size();

  /*
   * Section to assign the special species
   */
  m_SpecialSpecies = specialSpecies;
  if (specialSpecies == -1) {
    m_SpecialSpecies = m_Products[0];
  } 
  bool ifound = false;
  for (i = 0; i < (int) m_NetSpecies.size(); i++) {
    int ik = m_NetSpecies[i];
    if (ik == m_SpecialSpecies) {
      if (m_netStoich[i] > 0.0) {
	m_SpecialSpeciesProduct = true;
      } else {
	m_SpecialSpeciesProduct = false;
      }
      m_SS_index = i;
      ifound = true;
      break;
    }
  }
  if (!ifound) {
    throw ZuzaxError(":setupElemRxnVector",
		       "Special species not a reactant or product: " 
		       + int2str( m_SpecialSpecies));
  }


  m_ok = true;
  
}

std::string ExtraGlobalRxn::reactionString() {
  string rs;
  int k, istoich;
  for (k = 0; k < m_nReactants; k++) {
    int kkinspecies = m_Reactants[k];
    double stoich =  m_ReactantStoich[k];
    if (stoich != 1.0) {
      istoich = (int) stoich;
      if (fabs((double)istoich - stoich) < 0.00001) {
	rs += int2str(istoich) + " ";
      } else {
	rs += fp2str(stoich) + " ";
      }
    }
    string sName = m_kinetics->kineticsSpeciesName(kkinspecies);
    rs += sName;
    if (k < (m_nReactants-1)) {
      rs += " + ";
    }
  }
  rs += " = ";
  for (k = 0; k < m_nProducts; k++) {
    int kkinspecies = m_Products[k];
    double stoich =  m_ProductStoich[k];
    if (stoich != 1.0) {
      istoich = (int) stoich;
      if (fabs((double)istoich - stoich) < 0.00001) {
	rs += int2str(istoich) + " ";
      } else {
	rs += fp2str(stoich) + " ";
      }
    }
    string sName = m_kinetics->kineticsSpeciesName(kkinspecies);
    rs += sName;
    if (k < (m_nProducts-1)) {
      rs += " + ";
    }
  }
  return rs;
}

std::vector<int> & ExtraGlobalRxn::reactants() {
  return m_Reactants;
}

std::vector<int> & ExtraGlobalRxn::products() {
  return m_Products;
}

bool ExtraGlobalRxn::isReversible() {
  return m_reversible;
}

double ExtraGlobalRxn::reactantStoichCoeff(int kKin) {
  for (int k = 0; k < m_nReactants; k++) {
    int kkinspec = m_Reactants[k];
    if (kkinspec == kKin) {
      return m_ReactantStoich[k];
    }
  }
  return 0.0;
}
double ExtraGlobalRxn::productStoichCoeff(int kKin) {
  for (int k = 0; k < m_nProducts; k++) {
    int kkinspec = m_Products[k];
    if (kkinspec == kKin) {
      return m_ProductStoich[k];
    }
  }
  return 0.0;
}
  
double ExtraGlobalRxn::deltaSpecValue(double *speciesVectorProperty) {
  int k;
  double val = 0;
  for (k = 0; k < m_nNetSpecies; k++) {
    int kkinspec = m_NetSpecies[k];
    val += speciesVectorProperty[kkinspec] * m_netStoich[k];
  }
  return val;
}

double ExtraGlobalRxn::deltaRxnVecValue(double *rxnVectorProperty) {
  double val = 0;
  for (int i = 0; i < m_nRxns; i++) {
    val += m_ElemRxnVector[i] * rxnVectorProperty[i];
  }
  return val;
}



double ExtraGlobalRxn::ROPValue(double *ROPElemKinVector) {
  double val = 0.0;
  for (int i = 0; i < m_nRxns; i++) {
     double kstoich = m_kinetics->productStoichCoeff(m_SpecialSpecies, i) -
       m_kinetics->reactantStoichCoeff(m_SpecialSpecies, i);
     //if (!m_SpecialSpeciesProduct) {
     // kstoich = -kstoich;
     //}
    if (m_ElemRxnVector[i] > 0.0) {
      val += ROPElemKinVector[i] * kstoich * m_ElemRxnVector[i];
    } else {
      val -= ROPElemKinVector[i] * kstoich * m_ElemRxnVector[i];
    }
  }
  if (!m_SpecialSpeciesProduct) {
    val = -val;
  }
  return val;
}

double ExtraGlobalRxn::FwdROPValue(double *FwdROPElemKinVector,
                                   double *RevROPElemKinVector) {
  double val = 0.0;
  for (int i = 0; i < m_nRxns; i++) {    
    double kstoich = m_kinetics->productStoichCoeff(m_SpecialSpecies, i) -  m_kinetics->reactantStoichCoeff(m_SpecialSpecies, i);
    if (m_ElemRxnVector[i] > 0.0) {
      val += FwdROPElemKinVector[i]  * kstoich *  m_ElemRxnVector[i];
    }
    if (m_ElemRxnVector[i] < 0.0) {
      val += RevROPElemKinVector[i]   * kstoich *  m_ElemRxnVector[i];
    }
  }
  if (!m_SpecialSpeciesProduct) {
    val = -val;
  }
  return val;
}

double ExtraGlobalRxn::RevROPValue(double *FwdROPElemKinVector,
                                   double *RevROPElemKinVector) {
  double val = 0.0;
  for (int i = 0; i < m_nRxns; i++) {
    double kstoich = m_kinetics->productStoichCoeff(m_SpecialSpecies, i)-  m_kinetics->reactantStoichCoeff(m_SpecialSpecies, i);
    if (m_ElemRxnVector[i] > 0.0) {
      val +=  RevROPElemKinVector[i]  * kstoich *  m_ElemRxnVector[i];
    }
    if (m_ElemRxnVector[i] < 0.0) {
      val +=  FwdROPElemKinVector[i]  *  kstoich * m_ElemRxnVector[i];
    }
  }
  if (!m_SpecialSpeciesProduct) {
    val = -val;
  }
  return val;
}

