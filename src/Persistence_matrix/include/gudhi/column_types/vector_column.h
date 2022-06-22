/*    This file is part of the MMA Library - https://gitlab.inria.fr/dloiseau/multipers - which is released under MIT.
 *    See file LICENSE for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2022 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef VECTORCOLUMN_H
#define VECTORCOLUMN_H

#include <iostream>
#include <list>
#include <unordered_set>

#include "utilities.h"

namespace Vineyard {

class Vector_column
{
public:
    Vector_column();
    Vector_column(boundary_type& boundary);
    Vector_column(Vector_column& column);
    Vector_column(Vector_column&& column) noexcept;

    void get_content(boundary_type& container);
    bool contains(unsigned int value) const;
    bool is_empty();
    dimension_type get_dimension() const;
    int get_pivot();
    void clear();
    void clear(unsigned int value);
    void reorder(std::vector<index>& valueMap);
    void add(Vector_column& column);

    Vector_column& operator=(Vector_column other);

    friend void swap(Vector_column& col1, Vector_column& col2);

private:
    int dim_;
    std::vector<unsigned int> column_;
    std::unordered_set<unsigned int> erasedValues_;

    void _cleanValues();
};

inline Vector_column::Vector_column() : dim_(0)
{}

inline Vector_column::Vector_column(boundary_type &boundary)
    : dim_(boundary.size() == 0 ? 0 : boundary.size() - 1),
      column_(boundary)
{}

inline Vector_column::Vector_column(Vector_column &column)
    : dim_(column.dim_),
      column_(column.column_),
      erasedValues_(column.erasedValues_)
{}

inline Vector_column::Vector_column(Vector_column &&column) noexcept
    : dim_(std::exchange(column.dim_, 0)),
      column_(std::move(column.column_)),
      erasedValues_(std::move(column.erasedValues_))
{}

inline void Vector_column::get_content(boundary_type &container)
{
    _cleanValues();
    std::copy(column_.begin(), column_.end(), std::back_inserter(container));
}

inline bool Vector_column::contains(unsigned int value) const
{
    if (erasedValues_.find(value) != erasedValues_.end()) return false;

    for (unsigned int v : column_){
        if (v == value) return true;
    }
    return false;
}

inline bool Vector_column::is_empty()
{
    _cleanValues();
    return column_.empty();
}

inline dimension_type Vector_column::get_dimension() const
{
    return dim_;
}

inline int Vector_column::get_pivot()
{
    while (!column_.empty() &&
           erasedValues_.find(column_.back()) != erasedValues_.end()) {
        erasedValues_.erase(column_.back());
        column_.pop_back();
    }

    if (column_.empty()) return -1;

    return column_.back();
}

inline void Vector_column::clear()
{
    column_.clear();
    erasedValues_.clear();
}

inline void Vector_column::clear(unsigned int value)
{
    erasedValues_.insert(value);
}

inline void Vector_column::reorder(std::vector<index> &valueMap)
{
    std::vector<unsigned int> newColumn;
    for (unsigned int& v : column_) {
        if (erasedValues_.find(v) == erasedValues_.end())
            newColumn.push_back(valueMap.at(v));
    }
    std::sort(newColumn.begin(), newColumn.end());
    erasedValues_.clear();
    column_.swap(newColumn);
}

inline void Vector_column::add(Vector_column &column)
{
    if (column.is_empty()) return;
    if (column_.empty()){
        std::copy(column.column_.begin(), column.column_.end(), std::back_inserter(column_));
        return;
    }

    std::vector<unsigned int> newColumn;

    std::vector<unsigned int>::iterator itToAdd = column.column_.begin();
    std::vector<unsigned int>::iterator itTarget = column_.begin();
    unsigned int valToAdd = *itToAdd;
    unsigned int valTarget = *itTarget;

    while (itToAdd != column.column_.end() && itTarget != column_.end())
    {
        while (itToAdd != column.column_.end() &&
               column.erasedValues_.find(valToAdd) != column.erasedValues_.end()) {
            itToAdd++;
            valToAdd = *itToAdd;
        }

        while (itTarget != column_.end() &&
               erasedValues_.find(valTarget) != erasedValues_.end()) {
            itTarget++;
            valTarget = *itTarget;
        }

        if (itToAdd != column.column_.end() && itTarget != column_.end()){
            if (valToAdd == valTarget){
                itTarget++;
                itToAdd++;
            } else if (valToAdd < valTarget){
                newColumn.push_back(valToAdd);
                itToAdd++;
            } else {
                newColumn.push_back(valTarget);
                itTarget++;
            }
        }

        valToAdd = *itToAdd;
        valTarget = *itTarget;
    }

    while (itToAdd != column.column_.end()){
        newColumn.push_back(*itToAdd);
        itToAdd++;
    }

    while (itTarget != column_.end()){
        newColumn.push_back(*itTarget);
        itTarget++;
    }

    column_.swap(newColumn);
    erasedValues_.clear();
}

inline Vector_column &Vector_column::operator=(Vector_column other)
{
    std::swap(dim_, other.dim_);
    std::swap(column_, other.column_);
    std::swap(erasedValues_, other.erasedValues_);
    return *this;
}

inline void Vector_column::_cleanValues()
{
    std::vector<unsigned int> newColumn;
    for (unsigned int v : column_){
        if (erasedValues_.find(v) == erasedValues_.end())
            newColumn.push_back(v);
    }
    erasedValues_.clear();
    column_.swap(newColumn);
}

inline void swap(Vector_column& col1, Vector_column& col2)
{
    std::swap(col1.dim_, col2.dim_);
    col1.column_.swap(col2.column_);
    std::swap(col1.erasedValues_, col2.erasedValues_);
}

}   //namespace Vineyard

#endif // VECTORCOLUMN_H
