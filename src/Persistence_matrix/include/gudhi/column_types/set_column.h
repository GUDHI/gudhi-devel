/*    This file is part of the MMA Library - https://gitlab.inria.fr/dloiseau/multipers - which is released under MIT.
 *    See file LICENSE for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SETCOLUMN_H
#define SETCOLUMN_H

#include <iostream>
#include <list>
#include <set>

#include "utilities.h"

namespace Vineyard {

class Set_column
{
public:
    Set_column();
    Set_column(boundary_type& boundary);
    Set_column(Set_column& column);
    Set_column(Set_column&& column) noexcept;

    void get_content(boundary_type& container);
    bool contains(unsigned int value) const;
    bool is_empty();
    dimension_type get_dimension() const;
    int get_pivot();
    void clear();
    void clear(unsigned int value);
    void reorder(std::vector<index>& valueMap);
    void add(Set_column& column);

    Set_column& operator=(Set_column other);

    friend void swap(Set_column& col1, Set_column& col2);

private:
    int dim_;
    std::set<unsigned int> column_;
};

inline Set_column::Set_column() : dim_(0)
{}

inline Set_column::Set_column(boundary_type &boundary)
    : dim_(boundary.size() == 0 ? 0 : boundary.size() - 1),
      column_(boundary.begin(), boundary.end())
{}

inline Set_column::Set_column(Set_column &column)
    : dim_(column.dim_),
      column_(column.column_)
{}

inline Set_column::Set_column(Set_column &&column) noexcept
    : dim_(std::exchange(column.dim_, 0)),
      column_(std::move(column.column_))
{}

inline void Set_column::get_content(boundary_type &container)
{
    std::copy(column_.begin(), column_.end(), std::back_inserter(container));
}

inline bool Set_column::contains(unsigned int value) const
{
    return column_.find(value) != column_.end();
}

inline bool Set_column::is_empty()
{
    return column_.empty();
}

inline dimension_type Set_column::get_dimension() const
{
    return dim_;
}

inline int Set_column::get_pivot()
{
    if (column_.empty()) return -1;
    return *(column_.rbegin());
}

inline void Set_column::clear()
{
    column_.clear();
}

inline void Set_column::clear(unsigned int value)
{
    column_.erase(value);
}

inline void Set_column::reorder(std::vector<index> &valueMap)
{
    std::set<unsigned int> newSet;
    for (const unsigned int& v : column_) newSet.insert(valueMap.at(v));
    column_.swap(newSet);
}

inline void Set_column::add(Set_column &column)
{
    for (const unsigned int& v : column.column_){
        if (column_.find(v) != column_.end())
            column_.erase(v);
        else
            column_.insert(v);
    }
}

inline Set_column &Set_column::operator=(Set_column other)
{
    std::swap(dim_, other.dim_);
    std::swap(column_, other.column_);
    return *this;
}

inline void swap(Set_column& col1, Set_column& col2)
{
    std::swap(col1.dim_, col2.dim_);
    col1.column_.swap(col2.column_);
}

}   //namespace Vineyard

#endif // SETCOLUMN_H
