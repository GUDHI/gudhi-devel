template<class F, class I, class C, template<class E, class... A> class Col>
template<class ChainRange>
typename dionysus::SparseRowMatrix<F,I,C,Col>::Column
dionysus::SparseRowMatrix<F,I,C,Col>::
reduce(const ChainRange& chain_, IndexChain& trail)
{
    auto    row_cmp = [this](const Entry& e1, const Entry& e2)
                      { return this->cmp_(std::get<0>(e1.index()), std::get<0>(e2.index())); };

#define __DIONYSUS_USE_VECTOR_CHAINS    1

#if !(__DIONYSUS_USE_VECTOR_CHAINS)
    std::set<Entry,decltype(row_cmp)>   chain(row_cmp);
    for (auto x : chain_)
        chain.insert(Entry(x.element(), IndexPair(x.index(), 0)));
#else
    Column chain;
    for (auto x : chain_)
        chain.emplace_back(x.element(), IndexPair(x.index(), 0));
    std::sort(chain.begin(), chain.end(), row_cmp);
#endif

    typedef   Reduction<IndexPair>              ReductionIP;

    auto      chains   = [this](const IndexPair& rc) -> const Column&    { return this->col(std::get<1>(rc)); };
    auto      lows     = [this](const IndexPair& rc) -> IndexPair
                         {
                             Index r  = std::get<0>(rc);
                             auto  it = this->lows_.find(r);
                             if (it == this->lows_.end())
                                 return ReductionIP::unpaired;
                             else
                             {
                                 Index rr = std::get<0>(col(it->second).back().index());
                                 if (rr != r)
                                     std::cout << "Mismatch: " << rr << ' ' << r << std::endl;
                                 return IndexPair(r, it->second);
                             }
                         };

    auto      addto    = [&trail](FieldElement m, const IndexPair& rc)  { trail.emplace_back(m, std::get<1>(rc)); };

    ReductionIP::reduce(chain,
                        chains, lows,
                        field_, addto, row_cmp);

#if !(__DIONYSUS_USE_VECTOR_CHAINS)
    return Column(std::begin(chain), std::end(chain));
#else
    return chain;
#endif
}

template<class F, class I, class C, template<class E, class... A> class Col>
typename dionysus::SparseRowMatrix<F,I,C,Col>::Index
dionysus::SparseRowMatrix<F,I,C,Col>::
set(Index col, Column&& chain)
{
    Column& column = columns_.emplace(col, std::move(chain)).first->second;

    fix(col, column);

    Index r = std::get<0>(column.back().index());
    Index res;
    if (is_low(r))
        res = low(r);
    else
        res = col;
    lows_[r] = col;

    return res;
}

template<class F, class I, class C, template<class E, class... A> class Col>
void
dionysus::SparseRowMatrix<F,I,C,Col>::
fix(Index col, Column& column)
{
    for (auto& x : column)
    {
        std::get<1>(x.index()) = col;
        Index r = std::get<0>(x.index());
        row(r).push_back(x);
    }
}

template<class F, class I, class C, template<class E, class... A> class Col>
const typename dionysus::SparseRowMatrix<F,I,C,Col>::Row&
dionysus::SparseRowMatrix<F,I,C,Col>::
prepend_row(Index r, FieldElement m, const Row& chain)
{
    Row& new_row = row(r);

    for (auto& x : chain)
    {
        Index c = std::get<1>(x.index());
        Column& column = col(c);
        auto it = column.emplace(column.begin(), field().mul(x.element(), m), r, c);
        new_row.push_back(*it);
    }

    return new_row;
}
