//! \file ConfBlock.cpp
//! Brief description

//! \date   Jan 17, 2009
//! \author fuan
//!
//! Detailed description

#include <fstream>

#include "ConfBlock.h"

namespace confparser {

// ======= //
// Setters //
// ======= //

void ConfBlock::writeConfigFile( std::string fileName ) {

  std::ofstream f( fileName.c_str() );
  f << "# Automatically generated config file\n\n";
  writeConfigFileRec( f, "" );
  f.close();

}

// ========================== //
// Protected helper functions //
// ========================== //

void ConfBlock::writeConfigFileRec( std::ofstream &fileHandle, std::string indent ) {

  // Write out all data of current block
  std::map<std::string,std::string>::iterator mit;
  for ( mit = props_.begin(); mit != props_.end(); ++mit ) {
    fileHandle << indent << mit->first << " " << mit->second << ";\n";
  }
  fileHandle << '\n';

  // Process children, if any
  for ( childIter bit = children_.begin(); bit != children_.end(); ++bit ) {
    // Start new block
    fileHandle << indent << bit->first << " {\n\n";
    // Recurse into block
    bit->second.writeConfigFileRec( fileHandle, indent + "  " );
    // End block
    fileHandle << indent << "}\n\n";
  }

}

} // namespace confparser
