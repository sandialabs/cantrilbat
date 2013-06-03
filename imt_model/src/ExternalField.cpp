/*
 * $Id: InterfacialMassTransfer.cpp 208 2012-06-22 15:57:16Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"

#include "cantera/base/ctexceptions.h"

#include "ExternalField.h"


namespace Cantera {

  //======================================================================================================================
  ExternalFieldScalar::ExternalFieldScalar(double val) :
    value_init_init_(val), 
    value_init_(val),
    value_final_(val),
    value_final_final_(val),
    t_init_init_(0.0),
    t_init_(0.0),
    t_final_(0.0),
    t_final_final_(0.0),
    tb_(EF_TIMEBEHAVIOR_TFINALCONSTANT)
  {

  }
  //======================================================================================================================
  ExternalFieldScalar::ExternalFieldScalar(const ExternalFieldScalar &right) :
    value_init_init_(right.value_init_init_), 
    value_init_(right.value_init_),
    value_final_(right.value_final_),
    value_final_final_(right.value_final_final_),
    t_init_init_(right.t_init_init_),
    t_init_(right.t_init_),
    t_final_(right.t_final_),
    t_final_final_(right.t_final_final_),
    tb_(right.tb_)
  {

  }
  //======================================================================================================================
  ExternalFieldScalar & ExternalFieldScalar::operator=(const ExternalFieldScalar &right) 
  {
    if (this == &right) return *this;

    value_init_init_ = right.value_init_init_;
    value_init_ = right.value_init_;
    value_final_ = right.value_final_;
    value_final_final_ = right.value_final_final_;
    t_init_init_ = right.t_init_init_;
    t_init_ = right.t_init_;
    t_final_ = right.t_final_;
    t_final_final_ = right.t_final_final_;
    tb_ = right.tb_;

    return *this;
  }
  //======================================================================================================================
  void  ExternalFieldScalar::calcValue_final(double t_final)
  {
    t_final_ = t_final;
    if (t_final < t_init_init_) {
      throw CanteraError("ExternalFieldScalar::calcValue_final",
			 " Time periods are confused");
    }
    if (tb_ == EF_TIMEBEHAVIOR_TFINALCONSTANT) {
      value_final_ = value_final_final_;
    } else {
      double delTG =  t_final_final_ - t_init_init_;
      if (delTG <= 0.0) {
	value_final_ = value_final_final_;
      } else {
	double delT = t_final_ - t_init_init_;
	double incr = delT / delTG;
	value_final_ = incr * value_final_final_ + (1.0 - incr) * value_init_init_;
      }
    }
  }
  //======================================================================================================================
  void  ExternalFieldScalar::advanceValue_final(double t_final)
  {
    if (t_final > t_final_) {
      value_init_ = value_final_;
      t_init_ = t_final_;
    }
    calcValue_final(t_final);
  }
 //======================================================================================================================
  void  ExternalFieldScalar::advanceAndEndIntPeriod(double t_final)
  { 
    if (fabs(t_final - t_final_) > 1.0E-14) {
      throw CanteraError("ExternalFieldScalar::advanceAndEndIntPeriod()",
			 " Time periods are confused");
    }
    value_init_ = value_final_;
    t_init_ = t_final_;
  }
  //======================================================================================================================
  void  ExternalFieldScalar::calcValue_init(double t_init)
  {
    t_init_ = t_init;
    if (tb_ == EF_TIMEBEHAVIOR_TFINALCONSTANT) {
      value_init_ = value_final_final_;
    } else {
      double delTG =  t_final_final_ - t_init_init_;
      if (delTG <= 0.0) {
	value_init_ = value_final_final_;
      } else {
	double delT = t_final_ - t_init_init_;
	double incr = delT / delTG;
	value_init_ = incr * value_final_final_ + (1.0 - incr) * value_init_init_;
      }
    }
  }
  //======================================================================================================================
  void  ExternalFieldScalar::enterGlobalPeriod(double t_init_init, double val_init_init, double t_final_final,
					       double val_final_final)
  {
    value_init_init_ = val_init_init;
    value_final_final_ = val_final_final;
    t_init_init_ = t_init_init;
    t_final_final_ = t_final_final;
    value_init_ = val_final_final;
    value_final_ = val_final_final;
    t_init_ = t_init_init_;
    t_final_ = t_init_init_;
    if (tb_ == EF_TIMEBEHAVIOR_LINEARRAMP) {
      value_init_ = value_init_init_;
      value_final_ = value_init_init_;
    }
  }
 //======================================================================================================================
  void  ExternalFieldScalar::advanceAndEndPeriod(double t_final_final)
  {
    if (fabs(t_final_final - t_final_final_) > 1.0E-14) {
      throw CanteraError("ExternalFieldScalar::advanceAndEndPeriod()",
			 " Time periods are confused");
    }
    t_init_init_ = t_final_final_;
    t_init_ = t_init_init_;
    t_final_ = t_init_init_;
    value_init_init_ =  value_final_final_;
    value_init_ =  value_final_final_;
    value_final_ =  value_final_final_;
  }
  //======================================================================================================================
  void  ExternalFieldScalar::setNewGlobalTFinalValue(double t_final_final, double val_final_final)
  {
    t_final_final_ = t_final_final;
    value_final_final_ = val_final_final;
  } 
  //======================================================================================================================
  void  ExternalFieldScalar::check_t_globals(double t_init_init, double t_final_final)
  {
      if (fabs(t_final_final - t_final_final_) > 1.0E-14) {
      throw CanteraError("ExternalFieldScalar::check_t_globals()",
			 " Time periods are confused");
    } 
      if (fabs(t_init_init - t_init_init_) > 1.0E-14) {
	throw CanteraError("ExternalFieldScalar::check_t_globals()",
			   " Time periods are confused");
    }
  }

 //======================================================================================================================
 //======================================================================================================================
 //======================================================================================================================
  ExternalFieldVector::ExternalFieldVector(const std::vector<double> &val) :
    value_init_init_(val), 
    value_init_(val),
    value_final_(val),
    value_final_final_(val),
    t_init_init_(0.0),
    t_init_(0.0),
    t_final_(0.0),
    t_final_final_(0.0),
   
    tb_(EF_TIMEBEHAVIOR_TFINALCONSTANT)
  {
    len = value_init_init_.size();
  }
  //======================================================================================================================
  ExternalFieldVector::ExternalFieldVector(const ExternalFieldVector &right) :
    value_init_init_(right.value_init_init_), 
    value_init_(right.value_init_),
    value_final_(right.value_final_),
    value_final_final_(right.value_final_final_),
    t_init_init_(right.t_init_init_),
    t_init_(right.t_init_),
    t_final_(right.t_final_),
    t_final_final_(right.t_final_final_),
    len(right.len),
    tb_(right.tb_)
  {

  }
  //======================================================================================================================
  ExternalFieldVector & ExternalFieldVector::operator=(const ExternalFieldVector &right) 
  {
    if (this == &right) return *this;

    value_init_init_ = right.value_init_init_;
    value_init_ = right.value_init_;
    value_final_ = right.value_final_;
    value_final_final_ = right.value_final_final_;
    t_init_init_ = right.t_init_init_;
    t_init_ = right.t_init_;
    t_final_ = right.t_final_;
    t_final_final_ = right.t_final_final_;
    tb_ = right.tb_;
    len = right.len;

    return *this;
  }
  //======================================================================================================================
  void  ExternalFieldVector::calcValue_final(double t_final)
  {
    t_final_ = t_final;
    if (t_final < t_init_init_) {
      throw CanteraError("ExternalFieldVector::calcValue_final",
			 " Time periods are confused");
    }
    if (tb_ == EF_TIMEBEHAVIOR_TFINALCONSTANT) {
	for (size_t k = 0; k < len; k++) {
	  value_final_[k] = value_final_final_[k];
	}
    } else {
      double delTG =  t_final_final_ - t_init_init_;
      if (delTG <= 0.0) {
	for (size_t k = 0; k < len; k++) {
	  value_final_[k] = value_final_final_[k];
	}
      } else {
	double delT = t_final_ - t_init_init_;
	double incr = delT / delTG;
	for (size_t k = 0; k < len; k++) {
	  value_final_[k] = incr * value_final_final_[k] + (1.0 - incr) * value_init_init_[k];
	}
      }
    }
  }
  //======================================================================================================================
  void  ExternalFieldVector::advanceValue_final(double t_final)
  {
    if (t_final > t_final_) {
      value_init_ = value_final_;
      t_init_ = t_final_;
    }
    calcValue_final(t_final);
  }
 //======================================================================================================================
  void  ExternalFieldVector::advanceAndEndIntPeriod(double t_final)
  { 
    if (fabs(t_final - t_final_) > 1.0E-14) {
      throw CanteraError("ExternalFieldVector::advanceAndEndIntPeriod()",
			 " Time periods are confused");
    }
    value_init_ = value_final_;
    t_init_ = t_final_;
  }
  //======================================================================================================================
  void  ExternalFieldVector::calcValue_init(double t_init)
  {
    t_init_ = t_init;
    if (tb_ == EF_TIMEBEHAVIOR_TFINALCONSTANT) {
      for (size_t k = 0; k < len; k++) {
	value_init_[k] = value_final_final_[k];
      }
    } else {
      double delTG =  t_final_final_ - t_init_init_;
      if (delTG <= 0.0) {
	for (size_t k = 0; k < len; k++) {
	  value_init_[k] = value_final_final_[k];
	}
      } else {
	double delT = t_final_ - t_init_init_;
	double incr = delT / delTG;
	for (size_t k = 0; k < len; k++) {
	  value_init_[k] = incr * value_final_final_[k] + (1.0 - incr) * value_init_init_[k];
	}
      }
    }
  }
  //======================================================================================================================
  void  ExternalFieldVector::enterGlobalPeriod(double t_init_init, const std::vector<double> &val_init_init, double t_final_final,
					       const std::vector<double> &  val_final_final)
  {
    for (size_t k = 0; k < len; k++) {
      value_init_init_[k] = val_init_init[k];
      value_final_final_[k] = val_final_final[k];
      value_init_[k] = val_final_final[k];
      value_final_[k] = val_final_final[k];
    }
    t_init_init_ = t_init_init;
    t_final_final_ = t_final_final;
    t_init_ = t_init_init_;
    t_final_ = t_init_init_;
    if (tb_ == EF_TIMEBEHAVIOR_LINEARRAMP) {
      for (size_t k = 0; k < len; k++) {
	value_init_[k] = value_init_init_[k];
	value_final_[k] = value_init_init_[k];
      }
    }
  }
 //======================================================================================================================
  void  ExternalFieldVector::advanceAndEndPeriod(double t_final_final)
  {
    if (fabs(t_final_final - t_final_final_) > 1.0E-14) {
      throw CanteraError("ExternalFieldVector::advanceAndEndPeriod()",
			 " Time periods are confused");
    }
    t_init_init_ = t_final_final_;
    t_init_ = t_init_init_;
    t_final_ = t_init_init_;
    for (size_t k = 0; k < len; k++) {
      value_init_init_[k] =  value_final_final_[k];
      value_init_[k] =  value_final_final_[k];
      value_final_[k] =  value_final_final_[k];
    }
  }
  //======================================================================================================================
  void  ExternalFieldVector::setNewGlobalTFinalValue(double t_final_final, const std::vector<double> & val_final_final)
  {
    t_final_final_ = t_final_final; 
    for (size_t k = 0; k < len; k++) {
      value_final_final_[k] = val_final_final[k];
    }
  } 

  //======================================================================================================================
  void  ExternalFieldVector::check_t_globals(double t_init_init, double t_final_final)
  {
    if (fabs(t_final_final - t_final_final_) > 1.0E-14) {
      throw CanteraError("ExternalFieldVector::check_t_globals()",
			 " Time periods are confused");
    } 
    if (fabs(t_init_init - t_init_init_) > 1.0E-14) {
      throw CanteraError("ExternalFieldVector::check_t_globals()",
			 " Time periods are confused");
    }
  }
  
  //======================================================================================================================
}
//========================================================================================================================
