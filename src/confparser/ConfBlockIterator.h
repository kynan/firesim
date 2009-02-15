//! \file ConfBlockIterator.h
//! Declaration of the ConfBlockIterator class

//! \date   Jan 19, 2009
//! \author Florian Rathgeber

#ifndef CONFBLOCKITERATOR_H_
#define CONFBLOCKITERATOR_H_

#include <string>

//! Common namespace for all classes related to parsing of configuration files

namespace confparser {

  class ConfBlock;

  //! Iterator for ConfBlock trees

  class ConfBlockIterator {

    //! Declare ConfBlock as friend class

    friend class ConfBlock;

  public:

    //! Destructor

    //! Does nothing

    virtual ~ConfBlockIterator() {}

    // ==================== //
    // Overloaded operators //
    // ==================== //

    //! Advance to next matching ConfBlock

    void operator++();

    //! Check if the given and this iterator match

    bool operator==( const ConfBlockIterator& a );

    //! Check if the given and this iterator do not match

    bool operator!=( const ConfBlockIterator& a );

    //! Deference the iterator

    //! \return Reference to ConfBlock the iterator currently points to

    ConfBlock& operator*();

  protected:

    // ============================ //
    // Constructors and destructors //
    // ============================ //

    //! Constructor to create iterator

    //! \param[in] startBlock Initial block the iterator points to
    //! \param[in] name       Restricts the iterator to only iterate over blocks
    //!                       with given name (optional, defaults to empty)
    //! \param[in] maxDepth   Maximum recursion depth, counted from level where
    //!                       the iterator is created at. A value of 0 makes the
    //!                       iterator non-recursive, -1 makes no restriction to
    //!                       the recursion depth.
    //! \param[in] begin      Determines whether the iterator points to the
    //!                       start, i.e. this block, or the end, i.e. an
    //!                       invalid block

    ConfBlockIterator( ConfBlock* startBlock,
                       std::string name,
                       int maxDepth,
                       bool begin );

    // ============ //
    // Data members //
    // ============ //

    //! ConfBlock the iterator currently points to

    ConfBlock* b_;

    //! Restriction for the iterator only to iterate over blocks with given name

    std::string name_;

    //! Determines whether iterator has a restriction

    bool restricted_;

    //! Determines whether iterator is recursive, i.e. descends to children

    bool recursive_;

    //! Level at which the iterator was initialized

    //! It will not advance to a parent of this level

    int level_;

    //! Maximum recursion depth, counted from level at which the iterator was
    //! initialized. Does not affect non-recursive iterators

    int maxDepth_;

  };

}

#endif /* CONFBLOCKITERATOR_H_ */
