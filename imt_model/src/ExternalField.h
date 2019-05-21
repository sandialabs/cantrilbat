
 
#ifndef _EXTERNAL_FIELD_H
#define _EXTERNAL_FIELD_H

namespace Zuzax
{
  class IMT_KEY_INPUT;
}

namespace Zuzax
{

//! Enum type identifying the time behavior of the external fields
  enum EF_FieldTimeBehavior_Enum {
    //! backwards euler constant
    /*!
     *  Behavior within the global step is akin to backwards Euler. A step jump is 
     *  assumed to the global values at the end of the global time step even for
     *  intermediate times.
     */
    EF_TIMEBEHAVIOR_TFINALCONSTANT = 0,
    //! Linear Ramp
    /*!
     *  Behavior within the global step is treated as a linear function between the 
     *  beginning values and the end values
     */
    EF_TIMEBEHAVIOR_LINEARRAMP
  };




  class ExternalFieldScalar {
  public:
    //! constructor
    ExternalFieldScalar(double val = 0.0);

    ExternalFieldScalar(const ExternalFieldScalar &right);

    ExternalFieldScalar & operator=(const ExternalFieldScalar &right);

    void  calcValue_final(double t_final);

    void  advanceValue_final(double t_final);
      
    void  advanceAndEndIntPeriod(double t_final);

    void  calcValue_init(double t_init);

    void  enterGlobalPeriod(double t_init_init, double val_init_init, double t_final_final,
			    double val_final_final);
    
    void  advanceAndEndPeriod(double t_final_final);
    
    void  setNewGlobalTFinalValue(double t_final_final, double val_final_final);

    void check_t_globals(double t_init_init, double t_final_final);

    double value_init_init_;
    double value_init_;
    double value_final_;
    double value_final_final_;
    double t_init_init_;
    double t_init_;
    double t_final_;
    double t_final_final_;


    EF_FieldTimeBehavior_Enum tb_;
  };





 class ExternalFieldVector {
  public:
    //! constructor
    ExternalFieldVector(const std::vector<double> & val);

    ExternalFieldVector(const ExternalFieldVector &right);

    ExternalFieldVector & operator=(const ExternalFieldVector &right);

    void  calcValue_final(double t_final);

    void  advanceValue_final(double t_final);
      
    void  advanceAndEndIntPeriod(double t_final);

    void  calcValue_init(double t_init);

    void  enterGlobalPeriod(double t_init_init, const std::vector<double> & val_init_init, double t_final_final,
			    const std::vector<double> & val_final_final);
    
    void  advanceAndEndPeriod(double t_final_final);
    
    void  setNewGlobalTFinalValue(double t_final_final,   const std::vector<double> & vval_final_final);

    void check_t_globals(double t_init_init, double t_final_final);

    std::vector<double> value_init_init_;
    std::vector<double> value_init_;
    std::vector<double> value_final_;
    std::vector<double> value_final_final_;
    double t_init_init_;
    double t_init_;
    double t_final_;
    double t_final_final_;

    size_t len;

    EF_FieldTimeBehavior_Enum tb_;
  };

}

#endif
