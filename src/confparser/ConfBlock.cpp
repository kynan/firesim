//! \file ConfBlock.cpp
//! Brief description

//! \date   Jan 17, 2009
//! \author fuan
//!
//! Detailed description

#include "ConfBlock.h"

using boost::lexical_cast;

namespace confparser {

  // ============================ //
  // Constructors and destructors //
  // ============================ //

  ConfBlock::ConfBlock( std::string blockName, int level, ConfBlock* parent )
      : blockName_( blockName ),
        level_( level ),
        props_(),
        child_( NULL ),
        sibling_( NULL ),
        parent_( parent ) {}

  ConfBlock::~ConfBlock() {

    if ( child_ != NULL ) delete child_;

    // Remove the block from the siblings unless it is the outermost block
    if ( parent_ != NULL ) {

      // If the block is the first child, then its first sibling becomes the first
      // child of the parent
      if ( parent_->child_ == this ) {
        parent_->child_ = sibling_;
      // Else find the previous sibling and unlink the block from the chain of
      // siblings
      } else {
        ConfBlock* prevSib = parent_->child_;
        while ( prevSib->sibling_ != this ) prevSib = prevSib->sibling_;
        prevSib->sibling_ = sibling_;
      }
    }
  }

  // ======= //
  // Getters //
  // ======= //

  ConfBlock* ConfBlock::findSibling( std::string name ) {

    ConfBlock* sib = this;

    while ( sib->sibling_ != NULL ) {
      sib = sib->sibling_;
      if ( sib->blockName_ == name ) return sib;
    }

    return NULL;
  }

  // ======= //
  // Setters //
  // ======= //

//  bool ConfBlock::advance() {
//    if ( sibling_ != NULL ) {
//      this = sibling_;
//      return true;
//    }
//    return false;
//  }
//
//  bool ConfBlock::stepIn() {
//    if ( child_ != NULL ) {
//      this = child_;
//      return true;
//    }
//    return false;
//  }
//
//  bool ConfBlock::stepOut() {
//    if ( parent_ != NULL ) {
//      this = parent_;
//      return true;
//    }
//    return false;
//  }

}
