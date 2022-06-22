/*    This file is part of the MMA Library - https://gitlab.inria.fr/dloiseau/multipers - which is released under MIT.
 *    See file LICENSE for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef HEAPCOLUMN_H
#define HEAPCOLUMN_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "utilities.h"

namespace Vineyard {

class Heap_column
{
public:
    Heap_column();
    Heap_column(boundary_type& boundary);
    Heap_column(Heap_column& column);
    Heap_column(Heap_column&& column) noexcept;

    void get_content(boundary_type& container);
    bool contains(unsigned int value) const;
    bool is_empty();
    dimension_type get_dimension() const;
    int get_pivot();
    void clear();
    void clear(unsigned int value);
    void reorder(std::vector<index>& valueMap);
    void add(Heap_column& column);

    Heap_column& operator=(Heap_column other);

    friend void swap(Heap_column& col1, Heap_column& col2);

private:
    int dim_;
    std::vector<unsigned int> column_;
    unsigned int insertsSinceLastPrune_;
    std::unordered_set<unsigned int> erasedValues_;

    void _prune();
    int _pop_pivot();
};

inline Heap_column::Heap_column() : dim_(0), insertsSinceLastPrune_(0)
{}

inline Heap_column::Heap_column(boundary_type& boundary)
    : dim_(boundary.size() == 0 ? 0 : boundary.size() - 1),
      column_(boundary),
      insertsSinceLastPrune_(0)
{
    std::make_heap(column_.begin(), column_.end());
}

inline Heap_column::Heap_column(Heap_column& column)
    : dim_(column.dim_),
      column_(column.column_),
      insertsSinceLastPrune_(column.insertsSinceLastPrune_),
      erasedValues_(column.erasedValues_)
{}

inline Heap_column::Heap_column(Heap_column&& column) noexcept
    : dim_(std::exchange(column.dim_, 0)),
      column_(std::move(column.column_)),
      insertsSinceLastPrune_(std::exchange(column.insertsSinceLastPrune_, 0)),
      erasedValues_(std::move(column.erasedValues_))
{}

inline void Heap_column::get_content(boundary_type &container)
{
    _prune();
    container = column_;
    std::sort_heap(container.begin(), container.end());
}

inline bool Heap_column::contains(unsigned int value) const
{
    if (erasedValues_.find(value) != erasedValues_.end()) return false;

    unsigned int c = 0;

    for (unsigned int v : column_){
        if (v == value) c++;
    }

    return c % 2 != 0;
}

inline bool Heap_column::is_empty()
{
    int pivot = _pop_pivot();
    if (pivot != -1){
        column_.push_back(pivot);
        std::push_heap(column_.begin(), column_.end());
        return false;
    }
    return true;
}

inline dimension_type Heap_column::get_dimension() const
{
    return dim_;
}

inline int Heap_column::get_pivot()
{
    int pivot = _pop_pivot();
    if (pivot != -1){
        column_.push_back(pivot);
        std::push_heap(column_.begin(), column_.end());
    }
    return pivot;
}

inline void Heap_column::clear()
{
    column_.clear();
    insertsSinceLastPrune_ = 0;
    erasedValues_.clear();
}

inline void Heap_column::clear(unsigned int value)
{
    erasedValues_.insert(value);
}

inline void Heap_column::reorder(std::vector<index> &valueMap)
{
    std::vector<unsigned int> tempCol;
    int pivot = _pop_pivot();
    while (pivot != -1) {
        tempCol.push_back(valueMap.at(pivot));
        pivot = _pop_pivot();
    }
    column_.swap(tempCol);
    std::make_heap(column_.begin(), column_.end());

    insertsSinceLastPrune_ = 0;
    erasedValues_.clear();
}

inline void Heap_column::add(Heap_column &column)
{
    std::vector<unsigned int>& colToAdd = column.column_;
    const unsigned int size = colToAdd.size();

    if (size == 0) return;

    for (unsigned int v : colToAdd) {
        if (column.erasedValues_.find(v) == column.erasedValues_.end()){
            column_.push_back(v);
            std::push_heap(column_.begin(), column_.end());
            erasedValues_.erase(v);
        }
    }
    insertsSinceLastPrune_ += size;

    if (2 * insertsSinceLastPrune_ > column_.size()) _prune();
}

inline Heap_column& Heap_column::operator=(Heap_column other)
{
    std::swap(dim_, other.dim_);
    std::swap(column_, other.column_);
    std::swap(insertsSinceLastPrune_, other.insertsSinceLastPrune_);
    std::swap(erasedValues_, other.erasedValues_);
    return *this;
}

inline void Heap_column::_prune()
{
    if (insertsSinceLastPrune_ == 0 && erasedValues_.empty()) return;

    std::vector<unsigned int> tempCol;
    int pivot = _pop_pivot();
    while (pivot != -1) {
        tempCol.push_back(pivot);
        pivot = _pop_pivot();
    }
    column_.swap(tempCol);
    std::make_heap(column_.begin(), column_.end());

    insertsSinceLastPrune_ = 0;
    erasedValues_.clear();
}

inline int Heap_column::_pop_pivot()
{
    if (column_.empty()) {
        return -1;
    }

    unsigned int pivot = column_.front();
    std::pop_heap(column_.begin(), column_.end());
    column_.pop_back();
    while (!column_.empty() && column_.front() == pivot)
    {
        std::pop_heap(column_.begin(), column_.end());
        column_.pop_back();

        if (column_.empty()) {
            return -1;
        }
        pivot = column_.front();
        std::pop_heap(column_.begin(), column_.end());
        column_.pop_back();
    }

    if (erasedValues_.find(pivot) != erasedValues_.end())
        pivot = _pop_pivot();

    return pivot;
}

inline void swap(Heap_column& col1, Heap_column& col2)
{
    std::swap(col1.dim_, col2.dim_);
    col1.column_.swap(col2.column_);
    std::swap(col1.insertsSinceLastPrune_, col2.insertsSinceLastPrune_);
    std::swap(col1.erasedValues_, col2.erasedValues_);
}

}   //namespace Vineyard

#endif // HEAPCOLUMN_H
