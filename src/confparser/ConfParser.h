//! \file ConfParser.h
//! Declaration of the ConfParser class

//! \date   Jan 17, 2009
//! \author Florian Rathgeber

#ifndef CONFPARSER_H_
#define CONFPARSER_H_

#include <cstring>
#include <stdexcept>
#include <boost/lexical_cast.hpp>

#include "ConfBlock.h"

//! Common namespace for all classes related to parsing of configuration files

namespace confparser {

  //! Exception class for syntax errors

  //! Is capable of providing information about the file  and the line where
  //! the syntax error occured as well as a description of the error

  class BadSyntax : public std::logic_error {

  public:

    // ============================ //
    // Constructors and destructors //
    // ============================ //

    //! Default constructor

    //! Sets no specific information

    BadSyntax() throw()
        : std::logic_error( "Syntax error occured" ),
          file_( "" ),
          line_( -1 ),
          error_( "" ) {}

    //! Constructor

    //! Sets information about file and line of error occurrence as well as an
    //! error description

    BadSyntax( const std::string file, const int line, const std::string error ) throw()
        : std::logic_error( "Syntax error occured in file " + file + " on line "
            + boost::lexical_cast<std::string>( line ) + ": " + error ),
          file_( file ),
          line_( line ),
          error_( error ) {}

    //! Destructor

    virtual ~BadSyntax() throw() {}

    // ======= //
    // Getters //
    // ======= //

    //! Returns the file name

    //! \return Name of the file the syntax error occurred

    const char* getFileName() const throw() { return file_.c_str(); }

    //! Returns the line of the error

    //! \return Line where the syntax error occurred

    int getLine() const throw() { return line_; }

    //! Returns error description

    //! \return Description of the error

    const char* getErrorMsg() const throw() { return error_.c_str(); }

  private:

    // ============ //
    // Data members //
    // ============ //

    //! File where with erroneous syntax

    const std::string file_;

    //! Line of the syntax error

    const int line_;

    //! Description of the syntax error

    const std::string error_;

  };

  //! Parser for hierarchical configuration files

  //! Parses hierarchical configuration files of the following structure:
  //! \code
  //! key0 value0;
  //! ...
  //! block0 {
  //!   key1 value1;
  //!   ...
  //!   subblock0 {
  //!     key2 value2;
  //!     ...
  //!   }
  //!   subblock1 {
  //!     key3 value3;
  //!     ...
  //!   }
  //!   ...
  //! }
  //! block1 {
  //!   ...
  //! }
  //! ...
  //! \endcode
  //! The parsing result is a tree structure of ConfBlock objects and the actual
  //! key-value pairs are stored in a std::map

//  class ConfBlock;

  class ConfParser {

  public:

    // ============================ //
    // Constructors and destructors //
    // ============================ //

    //! Constructor

    //! Initializes root as empty block

    ConfParser() : outermost_() {}

    //! Destructor

    //! Does nothing

    virtual ~ConfParser() {}

    // ============= //
    // Functionality //
    // ============= //

    //! Parse a given configuration file

    //! \param[in] configFileName Path to the configuration file to parse
    //! \throw BadSyntax

    ConfBlock& parse( std::string configFileName ) throw( BadSyntax );

  protected:

    // ========================== //
    // Protected helper functions //
    // ========================== //

    //! Recursively parse a specific block

    //! \param[in] currBlock      Pointer to current block
    //! \param[in] level          Nesting depth of current block
    //! \param[in] nLine          Current line in the configuration file
    //! \param[in] configFile     Stream to write out to
    //! \param[in] configFileName Name of the configuration file parsed

    int parse_rec( ConfBlock* currBlock,
                   int level,
                   int nLine,
                   std::ifstream& configFile,
                   const std::string& configFileName );

    //! Root block of the parse tree that is built, will be returned by parse

    ConfBlock outermost_;

  };

} // namespace confparser

#endif /* CONFPARSER_H_ */
