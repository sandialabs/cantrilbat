/*
 * $Id: VoltageTable.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef VOLTAGE_TABLE_H
#define VOLTAGE_TABLE_H
//#include "sortAlgorithms.h"
#include "mdp_allo.h"
#include <vector>
#include <cstdio>
#include <algorithm>
using std::vector;

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/**
 *   This Class constructs a vector of voltages from which to make
 *   a table.
 */
class VoltageTable {

public:
  int    NPoints;
  bool   IncludeZero;
  bool   IncludeEzero;
  bool   IncludeEeq;
  double Vlow;                 //!<   Min temperature for thermo data fit
  double Vhigh;                //!<   Max temperature for thermo table
  double DeltaV;
  std::vector<double> V;
  int numAddedTs;
  std::vector<double> AddedVoltVector;
  double Ezero;
  double Eeq;
public:
  //====================================================================================================================
  /*
   * Default constructor for VoltageTable()
   */
  VoltageTable(const int nPts = 11, 
	       const double vlow = -1.0,
	       const double deltaV = 0.2,
	       const bool incZero = true,
	       const int numAdded = 0,
	       const double *addedVoltVector = 0) :
    NPoints(nPts), 
    IncludeZero(incZero),
    IncludeEzero(false),
    IncludeEeq(false),
    Vlow(vlow),
    Vhigh(0.0),
    DeltaV(deltaV),
    V(0),
    numAddedTs(numAdded),
    AddedVoltVector(0),
    Ezero(0.0),
    Eeq(0.0)
  {
    /****************************/
    int i;
    AddedVoltVector.resize(numAdded);
    for (i = 0; i < numAdded; i++) {
      AddedVoltVector[i] = addedVoltVector[i];
    }
 
    V.resize(NPoints + numAdded);

    double VCurrent = Vlow;
    for (i = 0; i < NPoints; i++) {
      if (fabs(VCurrent) < 1.0E-10) {
        VCurrent = 0.0;
      }
      V[i] = VCurrent;
      VCurrent += DeltaV;
    }
    for (i = 0; i < numAdded; i++) {
      V[NPoints + i] = addedVoltVector[i];
    }
    NPoints += numAdded;

    if (IncludeZero) {
      V.resize(NPoints + 1);
      V[NPoints] = 0.0;
      NPoints++;
    }
    //sort_dbl_1(&V[0], NPoints);
    std::sort(&V[0], &V[NPoints]);

  }
  //====================================================================================================================
  /*
   * Destructor
   */
  ~VoltageTable() {
  }
  //====================================================================================================================
  void AddEzero(double ezero) {
    int ipoint;
    if (IncludeEzero) {
      ipoint = findPoint(Ezero);
      if (ipoint == -1) {
        printf("AddEzero: Inconsistent logic Error exit\n");
	exit(-1);
      }
      removePoint(ipoint);
    }
    ipoint = findPoint(ezero);
    if (ipoint != -1) {
      return;
    }
    IncludeEzero = true;
    Ezero = ezero;
    V.push_back(Ezero);
    NPoints++;
    //sort_dbl_1(&V[0], NPoints);
    std::sort(&V[0], &V[NPoints]);
  }
  //====================================================================================================================
  void AddEeq(double eeq) {
    int ipoint;
    if (IncludeEeq) {
      ipoint = findPoint(Eeq);
      if (ipoint == -1) {
        printf("AddEeq: Inconsistent logic Error exit\n");
	exit(-1);
      }
      removePoint(ipoint);
    }
    ipoint =  findPoint(eeq);
    if (ipoint != -1) {
      return;
    }
    Eeq = eeq;
    IncludeEeq = true;
    V.push_back(Eeq);
    NPoints++;
    //sort_dbl_1(&V[0], NPoints);
    std::sort(&V[0], &V[NPoints]);
  }
  //====================================================================================================================
  int findPoint(const double Vval) {
    for (int i = 0; i < NPoints; i++) {
      if (fabs(Vval - V[i]) < 1.0E-5) {
	return i;
      }
    }
    return -1;  
  }
  //====================================================================================================================
  void removePoint(int ipoint) {
    if (ipoint >= 0 && ipoint < NPoints) {
      for (int jpoint = ipoint; jpoint < NPoints; jpoint++) {
	V[jpoint] = V[jpoint+1];
      }
      V.resize(NPoints -1);
      NPoints--;
    }
  }
  //====================================================================================================================
  void addPoint(double Vval) {
   
    V.resize(NPoints +1);
    V[NPoints] = Vval;
    NPoints++;
    //sort_dbl_1(&V[0], NPoints);
    std::sort(&V[0], &V[NPoints]);
  }
  //====================================================================================================================
  /*
   *  Overloaded operator[]
   * 
   *       return the array value in the vector
   */
  double operator[](const int i)  {
    return V[i];
  }
  //====================================================================================================================
  /*
   *  size()
   */
  int size() {
    return NPoints;
  }
  //====================================================================================================================
  /*
   * Block assignment and copy constructors: not needed.
   */
private:
  VoltageTable(const VoltageTable &);
  VoltageTable& operator=(const VoltageTable&);
};
  //====================================================================================================================
#endif
