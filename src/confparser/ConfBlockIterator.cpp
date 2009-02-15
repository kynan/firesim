//! \file ConfBlockIterator.cpp
//! Implementation of the ConfBlockIterator class

//! \date   Jan 19, 2009
//! \author Florian Rathgeber

#include <assert.h>

#include "ConfBlockIterator.h"
#include "ConfBlock.h"

namespace confparser {

  // ============================ //
  // Constructors and destructors //
  // ============================ //

  ConfBlockIterator::ConfBlockIterator( ConfBlock* startBlock,
                                        std::string name,
                                        int maxDepth,
                                        bool begin )
      : b_( 0 ),
        name_( name ),
        restricted_( name != "" ),
        recursive_( maxDepth != 0 ),
        level_( startBlock->level_ ),
        maxDepth_( maxDepth ) {

    if ( begin ) {
      assert( startBlock != NULL );
      b_ = startBlock;
      // Make sure that the iterator really points to a block with the given
      // name if restricted
      if ( restricted_ && b_->blockName_ != name ) this->operator++();
    }

  }

  // ==================== //
  // Overloaded operators //
  // ==================== //

  void ConfBlockIterator::operator++() {

    assert ( b_ != NULL );

    if ( recursive_ ) {
      // First try to descend to child if we have not yet reached maxDepth
      if ( b_->child_ != NULL && ( maxDepth_ == -1 || b_->level_ < level_ + maxDepth_ )  ) {
        b_ = b_->child_;
      // If there is no child try to advance to next sibling
      } else if ( b_->sibling_ != NULL ) {
        b_ = b_->sibling_;
      // If there is no more sibling and we are not at top level of iterator
      // yet, ascend to next sibling of parent
      } else if ( b_->level_ > level_ ) {
        b_ = b_->parent_->sibling_;
      // Else we cannot advance
      } else {
        b_ = NULL;
      }
    } else {
      b_ = b_->sibling_;
    }

    // If we are restricted and didn't find a block with given name, advance
    if ( restricted_ && b_ != NULL && b_->blockName_ != name_ ) operator++();
  }

  bool ConfBlockIterator::operator==( const ConfBlockIterator& a ) {
    // Equality means all members are equal
    return a.b_ == b_ && a.level_ == level_
                      && a.name_ == name_
                      && a.recursive_ == recursive_
                      && a.restricted_ == restricted_
                      && a.maxDepth_ == maxDepth_;
  }

  bool ConfBlockIterator::operator!=( const ConfBlockIterator& a ) {
    return !( a.b_ == b_ && a.level_ == level_
                         && a.name_ == name_
                         && a.recursive_ == recursive_
                         && a.restricted_ == restricted_
                         && a.maxDepth_ == maxDepth_ );
  }

  ConfBlock& ConfBlockIterator::operator*() {
    assert ( b_ != NULL );
    return *b_;
  }

}
