//! \file ConfBlock.h
//! Declaration of the ConfBlock class

//! \date   Jan 17, 2009
//! \author Florian Rathgeber

#ifndef CONFBLOCK_H_
#define CONFBLOCK_H_

#include <cstring>
#include <string>
#include <map>
#include <stdexcept>
#include <boost/lexical_cast.hpp>

//#include "ConfBlockIterator.h"

//! Common namespace for all classes related to parsing of configuration files

namespace confparser {

//! Exception class for non-existing parameters

//! Is capable of providing information about the parameter tried to retrieve
//! and the block tried to retrieve from

class ParameterNotFound : public std::invalid_argument {

public:

  // ============================ //
  // Constructors and destructors //
  // ============================ //

  //! Default constructor

  //! Sets no specific information

  ParameterNotFound() throw()
      : std::invalid_argument( "Tried to retrieve non-existing parameter" ),
        paramName_( "" ) {}

  //! Constructor

  //! Sets information about parameter name tried to retrieve as well as block
  //! name tried to retrieve from

  ParameterNotFound( const std::string paramName ) throw()
      : std::invalid_argument( "Tried to retrieve non-existing parameter "
          + paramName ),
        paramName_( paramName ) {}

  //! Destructor

  virtual ~ParameterNotFound() throw() {}

  // ======= //
  // Getters //
  // ======= //

  //! Returns the parameter name

  //! \return Name of the parameter that was tried to retrieve

  const char* getParamName() const throw() { return paramName_.c_str(); }

private:

  // ============ //
  // Data members //
  // ============ //

  //! Parameter name tried to retrieve

  const std::string paramName_;

};

//! Data structure corresponding to one block in the hierarchical configuration
//! file

class ConfBlock {

  //! Declare ConfParser as friend class to give it access to internals

  friend class ConfParser;

public:

  //! Iterator type for the children of the block

  typedef std::multimap< std::string, ConfBlock >::iterator childIter;
  typedef std::map< std::string, std::string >::iterator propIter;

  // ============================ //
  // Constructors and destructors //
  // ============================ //

  //! Default constructor to create a new (empty) block

  //! Parameters and children will be empty

  ConfBlock() : props_(), children_() {}

  //! Destructor

  //! Does nothing

  virtual ~ConfBlock() {}

  // ======= //
  // Getters //
  // ======= //

  //! Check whether the block is empty

  //! \return \em true if block is empty, i.e. children and parameters are
  //!         empty, \em false otherwise

  bool empty() {
    return props_.empty() && children_.empty();
  }

  //! Find a block with given name among the children

  //! \param[in] blockName Name of the block to retrieve
  //! \return Pointer to searched for block, \em NULL if none was found

  ConfBlock* find( std::string blockName ) {
    childIter it = children_.find( blockName );
    if ( it == children_.end() ) return NULL;
    return &(it->second);
  }

  //! Find all blocks with given name among the children

  //! \param[in] blockName Name of the blocks to retrieve
  //! \return Pair of iterators into the children pointing to the first and
  //!         after the last block found, are the same if none was found. Blocks
  //!         can be retrieved by ->second property of the iterator.

  std::pair< childIter, childIter > findAll( std::string blockName ) {
    return children_.equal_range( blockName );
  }

  //! Retrieve parameter as certain type

  //! \param[in] key Parameter name to retrieve
  //! \tparam    T   Type to retrieve parameter in
  //! \return Parameter specified by \em key as type \em T

  template < typename T >
  T getParam( std::string key ) throw ( ParameterNotFound ) {
    std::map< std::string, std::string >::iterator it = props_.find( key );
    if ( it == props_.end() )
      throw ParameterNotFound( key.c_str() );
    return boost::lexical_cast<T>( it->second );
  }

  //! Retrieve parameter as certain type

  //! \param[in]  key   Parameter name to retrieve
  //! \param[out] value Variable to store value in
  //! \tparam    T   Type to retrieve parameter in
  //! \return \e true if parameter was found, \e false if not
  //! \note Will not throw an exception, instead return value must be checked

  template < typename T >
  bool getParam( std::string key, T& value ) {
    propIter it = props_.find( key );
    if ( it == props_.end() ) return false;
    value = boost::lexical_cast<T>( it->second );
    return true;
  }

  // ======= //
  // Setters //
  // ======= //

  //! Add a child with given name to children of current block

  //! \param[in] blockName Block name of child to add
  //! \return Reference to newly added child

  ConfBlock& addChild( std::string blockName ) {
    childIter bit = children_.insert( std::make_pair( blockName, ConfBlock() ) );
    return bit->second;
  }

  //! Remove all children with given name

  //! \param[in] blockName Block name of children to remove
  //! \return Number of blocks removed, 0 if none of the given name were found

  int removeChildren( std::string blockName ) {
    return children_.erase( blockName );
  }

  //! Add a key value pair

  //! \param[in] key   Parameter name to add
  //! \param[in] value Parameter value to add
  //! \return \e true if a new key value pair was actually created, \e false if
  //!         the key already existed

  bool addParam( std::string key, std::string value ) {
    std::pair<propIter,bool> p = props_.insert( std::make_pair( key, value ) );
    return p.second;
  }

  //! Set parameter

  //! \param[in] key   Parameter name to set
  //! \param[in] value Parameter value to set
  //! \tparam    T     Type of the parameter to set (will be converted to string)
  //! \return \e true if the parameter was actually set, \e false if it does
  //!         not exist

  template < typename T >
  bool setParam( std::string key, T value ) {
    propIter it = props_.find( key );
    if ( it == props_.end() ) return false;
    it->second = boost::lexical_cast<std::string>( value );
    return true;
  }

  //! Remove parameter with a given key

  //! \param[in] key Parameter key to remove
  //! \return \e true if the key / value pair was removed, \e false if it was
  //!         not found

  bool removeParam( std::string key ) {
    return props_.erase( key );
  }

  //! Write out a configuration file of this block's subtree

  //! \param[in] fileName Name of the configuration file to write out

  void writeConfigFile( std::string fileName );

protected:

  // ========================== //
  // Protected helper functions //
  // ========================== //

  //! Recursive helper function to write the configuration of a subtree to a
  //! file

  //! \param[in] fileHandle Stream to write to
  //! \param[in] indent     Level of indentation for the current block

  void writeConfigFileRec( std::ofstream &fileHandle, std::string indent );

  // ============ //
  // Data members //
  // ============ //

  //! Map containing parameter key / value pairs defined in this block

  std::map< std::string, std::string > props_;

  //! Multimap containing the children

  std::multimap< std::string, ConfBlock > children_;

};

} // namespace confparser

#endif /* CONFBLOCK_H_ */
