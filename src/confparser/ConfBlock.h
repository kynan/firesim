//! \file ConfBlock.h
//! Declaration of the ConfBlock class

//! \date   Jan 17, 2009
//! \author Florian Rathgeber

#ifndef CONFBLOCK_H_
#define CONFBLOCK_H_

#include <cstring>
#include <assert.h>
#include <string>
#include <map>
#include <boost/lexical_cast.hpp>

#include "ConfBlockIterator.h"

//! Common namespace for all classes related to parsing of configuration files

namespace confparser {

//! Exception class for non-existing parameters

//! Is capable of providing information about the parameter tried to retrieve
//! and the block tried to retrieve from

class ParameterNotFound : public std::exception {

public:

  // ============================ //
  // Constructors and destructors //
  // ============================ //

  //! Default constructor

  //! Sets no specific information

  ParameterNotFound() throw()
      : blockName_( "" ),
        paramName_( "" ),
        what_( "Tried to retrieve non-existing parameter" ) {}

  //! Constructor

  //! Sets information about parameter name tried to retrieve as well as block
  //! name tried to retrieve from

  ParameterNotFound( const char* paramName, const char* blockName ) throw()
      : blockName_( blockName ), paramName_( paramName ), what_( 0 ) {
    int buff_size = strlen( paramName_ ) + strlen( blockName_ ) + 55;
    what_ = new char[buff_size];
    std::snprintf( what_, buff_size,
        "Tried to retrieve non-existing parameter %s from block %s", paramName_, blockName_ );
  }

  //! Destructor

  virtual ~ParameterNotFound() throw() {
    delete [] blockName_;
    delete [] paramName_;
    delete [] what_;
  }

  // ======= //
  // Getters //
  // ======= //

  //! Returns information about the error

  //! \return Information about the error incorporating parameter name tried to
  //! retrieve and block tried to retrieve from

  virtual const char* what() const throw() { return what_; }

  //! Returns the block name

  //! \return Name of the block the parameter that was tried to retrieve from

  const char* getBlockName() const throw() { return blockName_; }

  //! Returns the parameter name

  //! \return Name of the parameter that was tried to retrieve

  const char* getParamName() const throw() { return paramName_; }

private:

  // ============ //
  // Data members //
  // ============ //

  //! Block where the parameter was tried to retrieve from

  const char* blockName_;

  //! Parameter name tried to retrieve

  const char* paramName_;

  //! Message returned by what()

  char* what_;

};

//! Data structure corresponding to one block in the hierarchical configuration
//! file

class ConfBlock {

  // Declare ConfParser and ConfBlockIterator as friend class

  friend class ConfParser;
  friend class ConfBlockIterator;

public:

  // ============================ //
  // Constructors and destructors //
  // ============================ //

  //! Constructor to create a new block with given name

  //! \param[in] blockName Name of the block to create
  //! \param[in] level     Nesting depth of the current block (optional,
  //!                      defaults to 0)
  //! \param[in] parent    Pointer to parent block (optional, defaults to NULL)

  ConfBlock( std::string blockName, int level = 0, ConfBlock* parent = NULL );

  //! Destructor

  //! Deletes all children removes the block from the siblings and finally
  //! clears the block itself

  virtual ~ConfBlock();

  //! Unchains the block from its siblings and deletes all its offspring

  void unlink();

  // ======= //
  // Getters //
  // ======= //

  //! Retrieve parameter as certain type

  //! \param[in] key Parameter name to retrieve
  //! \tparam    T   Type to retrieve parameter in
  //! \return Parameter specified by \em key as type \em T

  template < typename T >
  T getParam( std::string key ) throw ( ParameterNotFound ) {
    std::map< std::string, std::string >::iterator it = props_.find( key );
    if ( it == props_.end() )
      throw ParameterNotFound( key.c_str(), getQualifiedName().c_str() );
    return boost::lexical_cast<T>( it->second );
  }

  //! Get the name of the block

  //! \return Name of this block

  std::string getName() { return blockName_; }

  //! Get the full qualified name of the block

  //! \return Full qualified name of the block, separated by dots

  std::string getQualifiedName();

  //! Get the child of the current block, call only if there is one

  //! \return First child of this block

  ConfBlock* descend() { assert( child_ != NULL ); return child_; }

  //! Get the sibling of the current block, call only if there is one

  //! \return Next sibling of this block

  ConfBlock* proceed() { assert( sibling_ != NULL ); return sibling_; }

  //! Get the parent of the current block, call only if there is one

  //! \return Parent of this block

  ConfBlock* ascend() { assert( parent_ != NULL ); return parent_; }

  //! Find the next sibling block with given name

  //! \param[in] name Name of the block to search for
  //! \return ConfBlockIterator pointing to the first sibling block with given
  //!         name, to the end if none was found

  ConfBlockIterator findSibling( std::string name ) {
    return ConfBlockIterator( this, name, 0, true );
  }

  //! Find the next block in the subtree of the current one with given name

  //! \param[in] name Name of the block to search for
  //! \return ConfBlockIterator pointing to the first block with given name in
  //!         the subtree of the current block, to the end if none was found

  ConfBlockIterator findRec( std::string name ) {
    return ConfBlockIterator( this, name, -1, true );
  }

  // ======= //
  // Setters //
  // ======= //

  //! Add a key value pair

  //! \param[in] key   Parameter name to add
  //! \param[in] value Parameter value to add

  void addParam( std::string key, std::string value ) {
    props_.insert( std::make_pair( key, value ) );
  }

  //! Set parameter

  //! \param[in] key   Parameter name to set
  //! \param[in] value Parameter value to set
  //! \tparam    T     Type of the parameter to set (will be converted to string)

  template < typename T >
  void setParam( std::string key, T value ) throw ( ParameterNotFound ) {
    std::map< std::string, std::string >::iterator it = props_.find( key );
    if ( it == props_.end() )
      throw ParameterNotFound( key.c_str(), getQualifiedName().c_str() );
    it->second = boost::lexical_cast<std::string>( value );
  }

  //! Write out a configuration file of this block's subtree

  //! \param[in] fileName Name of the configuration file to write out

  void writeConfigFile( std::string fileName );

  // ========= //
  // Iterators //
  // ========= //

  //! Returns an iterator pointing to the current block

  //! \param[in] name      Restricts the iterator to only iterate over blocks
  //!                      with given name (optional, defaults to empty)
  //! \param[in] maxDepth  Maximum recursion depth, counted from level where
  //!                      the iterator is created at. A value of 0 makes the
  //!                      iterator non-recursive, -1 makes no restriction to
  //!                      the recursion depth (the default).

  ConfBlockIterator begin( std::string name = "", int maxDepth = -1 ) {
    return ConfBlockIterator( this, name, maxDepth, true );
  }

  //! Returns an iterator pointing to the end, i.e. an invalid block

  //! \param[in] name      Restricts the iterator to only iterate over blocks
  //!                      with given name (optional, defaults to empty)
  //! \param[in] maxDepth  Maximum recursion depth, counted from level where
  //!                      the iterator is created at. A value of 0 makes the
  //!                      iterator non-recursive, -1 makes no restriction to
  //!                      the recursion depth (the default).

  ConfBlockIterator end( std::string name = "", int maxDepth = -1 ) {
    return ConfBlockIterator( this, name, maxDepth, false );
  }

protected:

  // ========================== //
  // Protected helper functions //
  // ========================== //

  //! Recursive helper function to write the configuration of a subtree to a
  //! file

  //! \param[in] fileHandle Stream to write to
  //! \param[in] initLvl    Level of the block writeConfigFile was called from

  void writeConfigFileRec( std::ofstream &fileHandle, int initLvl );

  // ============ //
  // Data members //
  // ============ //

  //! Name of this configuration block

  std::string blockName_;

  //! Nesting level in the block hierarchy

  int level_;

  //! Pointer to to key / value pairs defined in this block (may be NULL)

  std::map< std::string, std::string > props_;

  //! Pointer to first subblock (may be NULL)

  ConfBlock* child_;

  //! Pointer to next sibling block (may be NULL)

  ConfBlock* sibling_;

  //! Pointer to parent block (NULL only for the default block)

  ConfBlock* parent_;

};

} // namespace confparser

#endif /* CONFBLOCK_H_ */
