//! \file ConfBlock.cpp
//! Brief description

//! \date   Jan 17, 2009
//! \author fuan
//!
//! Detailed description

#include <fstream>
#include <boost/algorithm/string_regex.hpp>

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

std::string ConfBlock::getQualifiedName() {

  std::string qName = blockName_;
  ConfBlock* searchBlock = this;

  while ( searchBlock->level_ > 0 ) {
    searchBlock = searchBlock->parent_;
    qName = searchBlock->blockName_ + "." + qName;
  }

  return qName;
}

ConfBlock* ConfBlock::findSibling( std::string name ) {

  ConfBlock* sib = this;

  while ( sib->sibling_ != NULL ) {
    sib = sib->sibling_;
    if ( sib->blockName_ == name ) return sib;
  }

  return NULL;
}

ConfBlock* ConfBlock::findRec( std::string name ) {

  ConfBlock* sib = this;
  int initLvl = level_;

  while ( sib->blockName_ != name ) {
    if ( sib->child_ != NULL ) sib = sib->child_;
    else if ( sib->sibling_ != NULL ) sib = sib->sibling_;
    else if ( sib->level_ > initLvl && sib->level_ > 0 ) sib = sib->parent_;
    else return NULL;
  }

  return sib;
}

// ======= //
// Setters //
// ======= //

void ConfBlock::writeConfigFile( std::string fileName ) {

  std::ofstream f( fileName.c_str() );
  int prevLvl = level_;

  for ( ConfBlockIterator it = begin(); it != end(); ++it ) {

    ConfBlock currBlock = *it;
    std::string indent( 2 * ( currBlock.level_ - level_ ), ' ' );

    // Close brackets for every block closed since the previous
    for ( int i = prevLvl; i > currBlock.level_; --i ) {
      f << std::string( 2 * ( i - level_ ), ' ' ) << '}' << '\n';
    }

    f << indent << currBlock.blockName_ << " {" << '\n';

    indent.append( 2, ' ' );
    f << '\n';
    std::map<std::string,std::string>::iterator mit;
    for ( mit = currBlock.props_.begin(); mit != currBlock.props_.end(); ++mit ) {
      f << indent << mit->first << " " << mit->second << '\n';
    }
    f << '\n';
  }

}

} // namespace confparser
