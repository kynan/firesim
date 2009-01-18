//! \file ConfBlock.h
//! Declaration of the ConfBlock class

//! \date   Jan 17, 2009
//! \author Florian Rathgeber

#ifndef CONFBLOCK_H_
#define CONFBLOCK_H_

#include <string>
#include <map>
#include <boost/algorithm/string_regex.hpp>

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

    //! Get the parent

    //! \return Partent block of the current block

    ConfBlock* parent();

    //! Get the first child

    //! \return First child block of current block

    ConfBlock* child();

    //! Get the first sibling

    //! \return First sibling block of current block

    ConfBlock* sibling();

    // ======= //
    // Setters //
    // ======= //

    //! Add a key value pair

    //! \param[in] key   Parameter name to add
    //! \param[in] value Parameter value to add

    void addKey( std::string key, std::string value );

    //! Declare ConfParser as friend class
    friend class ConfParser;

  protected:

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
