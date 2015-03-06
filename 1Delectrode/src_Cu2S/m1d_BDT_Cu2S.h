/**
 * @file m1d_BulkDomainTypes.h
 */
/*
 *   $Id: m1d_BDT_Cu2S.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef M1D_BULKDOMAINTYPES_CU2S_H_
#define M1D_BULKDOMAINTYPES_CU2S_H_

#include "m1d_BulkDomainTypes.h"

namespace m1d
{

//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================

//! This class consists of one species diffusing using time dependence
/*!
 *  This class is used to test the implementation
 */
class BDT_Cu2S : public BulkDomainDescription
{
public:

  //! Constructor
  BDT_Cu2S(DomainLayout *dl_ptr, int id);

  //! Destructor
  virtual ~BDT_Cu2S();

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  BDT_Cu2S(const BDT_Cu2S &r);

  //! Assignment operator
  /*!
   * @param r    Object to be copied
   * @return     Returns a changeable reference to the current object
   */
  BDT_Cu2S &
  operator=(const BDT_Cu2S &r);

  virtual void SetEquationsVariablesList();


  //! Malloc and Return the object that will calculate the residual efficiently
  /*!
   * @return  Returns a pointer to the object that will calculate the residual
   *          efficiently
   */
  virtual BulkDomain1D *
  mallocDomain1D();

  //! Equation type and var type to apply them
  std::vector<VarType> EquationID;

};

}

#endif /* M1D_BULKDOMAINTYPES_H_ */
