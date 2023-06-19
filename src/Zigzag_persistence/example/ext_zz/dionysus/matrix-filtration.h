#pragma once

#include <cassert>
#include <vector>
#include <iostream>

namespace dionysus
{


template<class MatrixFiltration_>
class MatrixFiltrationCell;


// adapt Matrix as a Filtration to make it possible to feed into reduction algorithms
template <class Matrix_, class CellValue_>
class MatrixFiltration
{
    public:
        using Matrix = Matrix_;
        using CellValue = CellValue_;
        using Dimensions = std::vector<short unsigned>;
        using Values = std::vector<CellValue>;
        using Cell = MatrixFiltrationCell<MatrixFiltration>;


    public:
                MatrixFiltration(Matrix m, Dimensions dimensions, Values values):
                    m_(std::move(m)),
                    dimensions_(dimensions),
                    values_(values)                     { assert(m_->size() == dimensions_.size()); assert(m_->size() == values_.size()); }

        Cell            operator[](size_t i) const      { return Cell(this, i); }
        size_t          size() const                    { return m_.size(); }

        size_t          index(const Cell& c) const;

        Cell            begin() const                   { return Cell(this, 0); }
        Cell            end() const                     { return Cell(this, size()); }

        const Dimensions&   dimensions() const          { return dimensions_; }
        const Values&       values() const              { return values_; }

    private:
        Matrix          m_;
        Dimensions      dimensions_;
        Values          values_;

        friend class MatrixFiltrationCell<MatrixFiltration>;
};


template<class MatrixFiltration_>
class MatrixFiltrationCell
{
    public:
        using MatrixFiltration = MatrixFiltration_;
        using Matrix = typename MatrixFiltration::Matrix;
        using Data = typename MatrixFiltration::CellValue;
        using Field = typename Matrix::Field;

        template<class Field_ = Field>
        using Entry = ChainEntry<Field_, MatrixFiltrationCell>;

        template<class Field_ = Field>
        using BoundaryChain = std::vector<Entry<Field_>>;

    public:
                MatrixFiltrationCell(const MatrixFiltration* mf, size_t i):
                    mf_(mf), i_(i)      {}

        short unsigned  dimension() const       { return mf_->dimensions_[i_]; }
        const Data&     data() const            { return mf_->values_[i_]; }

        bool            operator==(const MatrixFiltrationCell& other) const     { return i_ == other.i_; }
        bool            operator!=(const MatrixFiltrationCell& other) const     { return i_ != other.i_; }

        BoundaryChain<> boundary() const
        {
            BoundaryChain<> bdry;
            for (auto& entry : (mf_->m_)[i_])
                bdry.emplace_back(Entry<> { entry.e, MatrixFiltrationCell(mf_, entry.i) });
            return bdry;
        }

        template<class Field_>
        BoundaryChain<Field_>   boundary(const Field_& field) const
        {
            BoundaryChain<Field_> bdry;
            for (auto& entry : (mf_->m_)[i_])
                bdry.emplace_back(Entry<Field_> { field.init(entry.e), MatrixFiltrationCell(mf_, entry.i) });
            return bdry;
        }

        // iterator interface
        MatrixFiltrationCell    operator++(int)         { MatrixFiltrationCell copy = *this; i_++; return copy; }
        MatrixFiltrationCell&   operator++()            { ++i_; return *this; }

        const MatrixFiltrationCell& operator*() const   { return *this; }
        MatrixFiltrationCell&       operator*()         { return *this; }

        size_t          i() const       { return i_; }

        friend
        std::ostream&   operator<<(std::ostream& out, const MatrixFiltrationCell& c)
        { out << c.i_; return out; }

    private:
        const MatrixFiltration* mf_ = nullptr;
        size_t                  i_;
};

template<class Matrix, class CellValue>
size_t
MatrixFiltration<Matrix,CellValue>::index(const Cell& c) const
{
    return c.i();
}

}
