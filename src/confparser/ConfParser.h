//! \file ConfParser.h
//! Declaration of the ConfParser class

//! \date   Jan 17, 2009
//! \author Florian Rathgeber

#ifndef CONFPARSER_H_
#define CONFPARSER_H_

//#include <cstdio>
#include <cstring>
#include <exception>

//! Common namespace for all classes related to parsing of configuration files

namespace confparser {

  //! Exception class for syntax errors

  //! Is capable of providing information about the file  and the line where
  //! the syntax error occured as well as a description of the error

  class BadSyntax : public std::exception {

  public:

    // ============================ //
    // Constructors and destructors //
    // ============================ //

    //! Default constructor

    //! Sets no specific information

    BadSyntax() throw() : file_( "" ), line_( -1 ), error_( "" ) {}

    //! Constructor

    //! Sets information about file and line of error occurrence as well as an
    //! error description

    BadSyntax( const char* file, const int line, const char* error ) throw()
        : file_( file ), line_( line ), error_( error ), what_( 0 ) {
      int buff_size = strlen( file_ ) + strlen( error_ ) + 50;
      what_ = new char[buff_size];
      std::snprintf( what_, buff_size,
          "Syntax error occurred in file %s on line %i: %s", file_, line_, error_ );
    }

    //! Destructor

    virtual ~BadSyntax() throw() {
      delete [] file_;
      delete [] error_;
      delete [] what_;
    }

    // ======= //
    // Getters //
    // ======= //

    //! Returns information about the syntax error

    //! \return Information about the syntax error incorporating file and line
    //! of the error as well as the description

    virtual const char* what() const throw() { return what_; }

    //! Returns the file name

    //! \return Name of the file the syntax error occurred

    const char* getFileName() const throw() { return file_; }

    //! Returns the line of the error

    //! \return Line where the syntax error occurred

    int getLine() const throw() { return line_; }

    //! Returns error description

    //! \return Description of the error

    const char* getErrorMsg() const throw() { return error_; }

  private:

    // ============ //
    // Data members //
    // ============ //

    //! File where with erroneous syntax

    const char* file_;

    //! Line of the syntax error

    const int line_;

    //! Description of the syntax error

    const char* error_;

    //! Message returned by what()

    char* what_;

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

  class ConfBlock;

  class ConfParser {

  public:

    // ============================ //
    // Constructors and destructors //
    // ============================ //

    //! Constructor

    //! Does nothing

    ConfParser() {}

    //! Destructor

    //! Does nothing

    virtual ~ConfParser() {}

    // ============= //
    // Functionality //
    // ============= //

    //! Parse a given configuration file

    //! \param[in] configFileName Path to the configuration file to parse
    //! \throw BadSyntax

    ConfBlock* parse( std::string configFileName ) throw( BadSyntax );

  };

} // namespace confparser

#endif /* CONFPARSER_H_ */
