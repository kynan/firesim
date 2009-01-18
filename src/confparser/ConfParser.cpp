//! \file ConfParser.cpp
//! Implementation of the ConfParser class

//! \date   Jan 17, 2009
//! \author Florian Rathgeber

#include <fstream>
#include <boost/algorithm/string_regex.hpp>

#include "ConfParser.h"

using namespace std;
using namespace boost;

namespace confparser {

  ConfBlock* ConfParser::parse( string configFileName ) throw( BadSyntax ) {

    // Open configuration file
    ifstream configFile( configFileName.c_str() );

    // Create default outermost block
    ConfBlock* outermost = new ConfBlock( "__default" );
    // ... which is of level 0
    int level = 0;
    // Pointer to currently processed block
    ConfBlock* currBlock = outermost;
    // Pointer to previous child of parent block
    ConfBlock* prevChild = NULL;

    // Read in config file linewise
    string l;
    int nLine = 0;
    while ( getline( configFile, l ) ) {

      ++nLine;

      // Remove leading and trailing whitespace
      trim( l );
      // If the line contained only whitspace we can skip it
      if ( l.empty() ) continue;

//      size_t n = l.find_first_not_of( " \r\n\t" );
//      if ( n != std::string::npos ) {
//        l.erase( 0, n );
//      // if we do not find any non-whitespace character we can continue anyway
//      } else {
//        continue;
//      }

      // Ignore lines beginning with # or // as comments
      if ( starts_with( l, "#") || starts_with( l, "//" ) ) continue;
//      if ( ( l.length() > 0 && l[0] == '#' )
//           || ( l.length() > 1 && l[0] == '/' && l[1] == '/' ) )
//        continue;

      // Check whether it is the end of a block
      if ( l[0] == '}' ) {

        // Syntax error if we are at level 0
        if ( level == 0 ) {
          throw BadSyntax( configFileName.c_str(), nLine,
              "Found closing block at outermost level" );
        }

        assert( currBlock->parent_ != NULL );
        currBlock = currBlock->parent_;
        level--;
        continue;
      }

      cmatch m;
      regex blockBegin( "(\\w+)\\s*\\{" );
      regex keyVal( "(\\w+)\\s+(.*)" );

      // Check whether it is the beginning of a new block
      if ( regex_match( l.c_str(), m, blockBegin ) ) {

        ++level;
        ConfBlock* newBlock = new ConfBlock( m[1], level, currBlock );

        // Check if the current block already has a child, then the new block
        // is the sibling of the previous child
        if ( currBlock->child_ != NULL ) {
          currBlock->child_ = newBlock;
        // Else the new block is the current block's first child
        } else {
          assert( prevChild != NULL );
          prevChild->sibling_ = newBlock;
        }
        currBlock = newBlock;
        prevChild = newBlock;

      // Check whether it is a key / value pair
      } else if ( regex_match( l.c_str(), m, keyVal ) ) {

        currBlock->addKey( m[1], m[2] );

      // Else we have a malformed expression and throw an exception
      } else {

        throw BadSyntax( configFileName.c_str(), nLine, "Malformed expression" );

      }

    }

    return outermost;

  }

}
