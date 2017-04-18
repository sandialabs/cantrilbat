/**
 * @file m1d_VBRIndices.h 
 *        Header for the class LocalRowNodeVBRIndices that manages the interaction with the VBR matrix and the storage
 *        for the blocks of the VBR matrix (see \ref matrixInteraction
 *        and class \link m1d::LocalRowNodeVBRIndices LocalRowNodeVBRIndices\endlink).
 */

#include "m1d_RecordTree_base.h"
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
RecordTree_base::RecordTree_base(RecordTree_base* parent) :
    GeneralName_(""),
    SpecificName_(""),
    numSubBE_(0),
    SubBE_(0),
    parent_(parent)
{
}
//==================================================================================================================================
RecordTree_base::~RecordTree_base()
{
    for (int i = 0; i < numSubBE_; i++) {
        delete SubBE_[i];
        SubBE_[i] = 0;
    }
}
//==================================================================================================================================
RecordTree_base::RecordTree_base(const RecordTree_base& right, RecordTree_base* parent) :
    GeneralName_(""),
    SpecificName_(""),
    numSubBE_(0),
    SubBE_(0),
    parent_(parent)
{
    RecordTree_base::operator=(right);
}
//==================================================================================================================================
RecordTree_base& RecordTree_base::
operator=(const RecordTree_base& right)
{
    if (this == &right) {
        return *this;
    }

    for (int i = 0; i < numSubBE_; i++) {
        delete SubBE_[i];
        SubBE_[i] = 0;
    }
    GeneralName_ = right.GeneralName_;
    SpecificName_ = right.SpecificName_;
    numSubBE_ = right.numSubBE_;
    SubBE_.resize(numSubBE_,0);
    for (int i = 0; i < numSubBE_; i++) {
        SubBE_[i] = (right.SubBE_[i])->duplicateMyselfAsRecordTree(this);
    }
    return *this;
}
//==================================================================================================================================
RecordTree_base* RecordTree_base::duplicateMyselfAsRecordTree(RecordTree_base* parent)
{
    RecordTree_base* rt = new RecordTree_base(*this, parent);
    return rt;
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
