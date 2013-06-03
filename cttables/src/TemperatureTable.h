/*
 * $Id: TemperatureTable.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef TEMPERATURE_TABLE_H
#define  TEMPERATURE_TABLE_H
#include "sortAlgorithms.h"
#include "mdp_allo.h"
#include <vector>
using std::vector;

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/**
 *   This Class constructs a vector of temperature from which to make
 *   a table.
 */
class TemperatureTable {

public:
    int    NPoints;
    bool   Include298;
    double Tlow;                 //!<   Min temperature for thermo data fit
    double Thigh;                //!<   Max temperature for thermo table
    double DeltaT;
    double *T;
    int numAddedTs;
    double *AddedTempVector;
public:
  /*
   * Default constructor for TemperatureTable()
   */
  TemperatureTable(const int nPts = 14, 
		   const bool inc298 = true,
		   const double tlow = 300.,
                   const double deltaT = 100.,
		   const int numAdded = 0,
                   const double *addedTempVector = 0) :
      NPoints(nPts), 
      Include298(inc298), 
      Tlow(tlow), 
      DeltaT(deltaT),
      T(0),
      numAddedTs(numAdded) {
    /****************************/
      int i;
      AddedTempVector = mdpUtil::mdp_alloc_dbl_1(numAdded, 0.0);
      mdpUtil::mdp_copy_dbl_1(AddedTempVector, addedTempVector, numAdded);
      T = mdpUtil::mdp_alloc_dbl_1(NPoints, 0.0);

      double TCurrent = Tlow;
      for (i = 0; i < NPoints; i++) {
	T[i] = TCurrent;
	TCurrent += DeltaT;
      }
      if (Include298) {
	mdpUtil::mdp_realloc_dbl_1(&T, NPoints+1, NPoints, 298.15);
	NPoints++;
      }
      if (numAdded > 0) {
	mdpUtil::mdp_realloc_dbl_1(&T, NPoints+numAdded, NPoints, 0.0);
	for (i = 0; i < numAdded; i++) {
	  T[i+NPoints] = addedTempVector[i];
	}
	NPoints += numAdded;
      }
 
      sort_dbl_1(T, NPoints);

    
  }
  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/
  /*
   * Destructor
   */
  ~TemperatureTable() {
      mdpUtil::mdp_safe_free((void **) &AddedTempVector);
      mdpUtil::mdp_safe_free((void **) &T);
  }

  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/
  /*
   *  Overloaded operator[]
   * 
   *       return the array value in the vector
   */
  double operator[](const int i)  {
    return T[i];
  }
  /***********************************************************************/
  /***********************************************************************/
  /***********************************************************************/
  /*
   *  size()
   */
  int size() {
      return NPoints;
  }
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
  /*
   * Block assignment and copy constructors: not needed.
   */
private:
  TemperatureTable(const TemperatureTable &);
  TemperatureTable& operator=(const TemperatureTable&);
};
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
#endif
