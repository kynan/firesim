//! \file ConfBlock.h
//! Declaration of the ConfBlock class

//! \date   Jan 17, 2009
//! \author Florian Rathgeber

#ifndef CONFBLOCK_H_
#define CONFBLOCK_H_

#include <string>
#include <map>
#include <boost/algorithm/string_regex.hpp>
#include <boost/lexical_cast.hpp>

//! Common namespace for all classes related to parsing of configuration files

namespace confparser {

  //! Data structure corresponding to one block in the hierarchical configuration
  //! file

  class ConfBlock {

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

    //! Find the next sibling block with given name

    //! \param[in] name Name of the block to search for
    //! \return Pointer to the first sibling block with given name, \em NULL if
    //!         none was found

    ConfBlock* findSibling( std::string name );

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

//    //! Move on to next sibling
//
//    //! After this operator \em this will point to the next sibling if any
//    //! \return \em true if succeeded, \em false if there is no further sibling
//
//    inline bool advance();
//
//    //! Move on to first child
//
//    //! After this operator \em this will point to the first child if any
//    //! \return \em true if succeeded, \em false if there is no child
//
//    inline bool stepIn();
//
//    //! Move out to parent
//
//    //! After this operator \em this will point to the parent if any
//    //! \return \em true if succeeded, \em false if there is no parent (i.e. it
//    //! is the outermost block)
//
//    inline bool stepOut();

  protected:

    //! Declare ConfParser as friend class

    friend class ConfParser;

    // ============ //
    // Data members //
    // ============ //

    // Name of this configuration block

    std::string blockName_;

    // Nesting level in the block hierarchy

    int level_;

    // Pointer to to key / value pairs defined in this block (may be NULL)

    std::map< std::string, std::string > props_;

    // Pointer to first subblock (may be NULL)

    ConfBlock* child_;

    // Pointer to next sibling block (may be NULL)

    ConfBlock* sibling_;

    //Pointer to parent block (NULL only for the default block)

    ConfBlock* parent_;

  };

}

#endif /* CONFBLOCK_H_ */
