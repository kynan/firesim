//! \file ConfParser.cpp
//! Implementation of the ConfParser class

//! \date   Jan 17, 2009
//! \author Florian Rathgeber

#include <iostream>
#include <fstream>
#include <boost/algorithm/string_regex.hpp>

#include "ConfParser.h"
#include "ConfBlock.h"

using namespace std;
using namespace boost;

namespace confparser {

ConfBlock& ConfParser::parse( string configFileName ) throw( BadSyntax ) {

  // Open configuration file
  ifstream configFile( configFileName.c_str(), ios::out );
  if ( !configFile ) throw BadSyntax( configFileName, 0,
      "Specified configuration file does not exist" );

#ifdef DEBUG
    cout << "Parsing configuration file " << configFileName << endl;
#endif

  parse_rec( &outermost_, 0, 0, configFile, configFileName );

  return outermost_;

}

int ConfParser::parse_rec( ConfBlock* currBlock,
                           int level,
                           int nLine,
                           ifstream& configFile,
                           const string& configFileName ) {

  // Read in config file linewise
  string l;
  while ( getline( configFile, l ) ) {

    ++nLine;

    // Remove leading and trailing whitespace
    trim( l );
    // If the line contained only whitspace we can skip it
    if ( l.empty() ) continue;

    // Ignore lines beginning with # or // as comments
    if ( starts_with( l, "#") || starts_with( l, "//" ) ) continue;

    // Check whether it is the end of a block
    if ( l.length() == 1 && l[0] == '}' ) {

      // Syntax error if we are at level 0
      if ( level == 0 ) {
        throw BadSyntax( configFileName, nLine,
            "Found closing block at outermost level" );
      }

      return nLine;
    }

    cmatch m;
    regex blockBegin( "^(\\w+)\\s*\\{$" );
    regex keyVal( "^(\\w+)\\s+(.*);" );

    // Check whether it is the beginning of a new block
    if ( regex_match( l.c_str(), m, blockBegin ) ) {

#ifdef DEBUG
        cout << "Adding Block " << m[1] << " at level " << level << " (line " << nLine << ")" << endl;
#endif

      nLine = parse_rec( &currBlock->addChild( m[1] ), level + 1, nLine, configFile, configFileName );

    // Check whether it is a key / value pair
    } else if ( regex_match( l.c_str(), m, keyVal ) ) {

      currBlock->addParam( m[1], m[2] );

    // Else we have a malformed expression and throw an exception
    } else {

      throw BadSyntax( configFileName, nLine, "Malformed expression" );

    }

  }

  // check if we are at outermost level again at the end
  if ( level != 0 )
    throw BadSyntax( configFileName, nLine, "Unexpected end of configuration file" );

  return nLine;
}

} // namespace
