//! \file ConfBlock.h
//! Declaration of the ConfBlock class

//! \date   Jan 17, 2009
//! \author Florian Rathgeber

#ifndef CONFBLOCK_H_
#define CONFBLOCK_H_

#include <string>
#include <map>
#include <boost/lexical_cast.hpp>

#include "ConfBlockIterator.h"

//! Common namespace for all classes related to parsing of configuration files

namespace confparser {

  //! Data structure corresponding to one block in the hierarchical configuration
  //! file

  class ConfBlock {

    //! Declare ConfParser and ConfBlockIterator as friend class

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
    T getParam( std::string key ) {
      return boost::lexical_cast<T>( props_[ key ] );
    }

    //! Get the full qualified name of the block

    //! \return Full qualified name of the block, separated by dots

    std::string getQualifiedName();

    //! Find the next sibling block with given name

    //! \param[in] name Name of the block to search for
    //! \return Pointer to the first sibling block with given name, \em NULL if
    //!         none was found

    ConfBlock* findSibling( std::string name );

    //! Find the next block in the subtree of the current one with given name

    //! \param[in] name Name of the block to search for
    //! \return Pointer to the first block with given name in the subtree of
    //!         the current block, \em NULL if none was found

    ConfBlock* findRec( std::string name );

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
    void setParam( std::string key, T value ) {
      props_[ key ] = boost::lexical_cast<std::string>( value );
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

}

#endif /* CONFBLOCK_H_ */
