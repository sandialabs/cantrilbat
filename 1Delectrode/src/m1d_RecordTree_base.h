/**
 * @file m1d_RecordTree_base.h
 *
 */

/*
 *  $Id: m1d_RecordTree_base.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
#include <string>
#include <vector>

#ifndef _M1D_RECORDTREE_BASE_H
#define _M1D_RECORDTREE_BASE_H

namespace m1d
{

  //! Record tree class
  /*!
   *   Searches are carried out by string. Searches can either be by a generic
   *   string or by a specific string. The generic string search refers to the
   *   subclass of the record. The specific string refers to the name of the object.
   */
  class RecordTree_base {
  public:

    //! Constructor
    RecordTree_base(RecordTree_base *parent = 0);


    RecordTree_base(const RecordTree_base &right, RecordTree_base *parent = 0);

    //! Destructor
    virtual ~RecordTree_base();



    RecordTree_base&operator=(const RecordTree_base &right);

    virtual RecordTree_base *duplicateMyselfAsRecordTree(RecordTree_base *parent = 0);

    //! Generic name of the record
    std::string GeneralName_;

    //! Specific name of the record
    std::string SpecificName_;

    //! Number of children in the current record
    int numSubBE_;

    //! Vector of children to the current record
    std::vector<RecordTree_base *> SubBE_;

    //! Parent to the current position
    RecordTree_base *parent_;

  };

}

#endif
