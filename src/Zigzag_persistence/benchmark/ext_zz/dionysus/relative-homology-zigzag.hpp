template<class F, class I, class C>
template<class ChainRange>
void
dionysus::RelativeHomologyZigzag<F,I,C>::
add_both(const ChainRange& chain)
{
    zzp_.add(absolute_chain(chain));
    op_map_.insert( { zzp_op_++, op_ } );
    absolute_.left.insert( { cell_, zzp_cell_++ } );

    zzp_.add(relative_chain(cell_, chain));
    op_map_.insert( { zzp_op_++, op_ } );
    relative_.left.insert( { cell_, zzp_cell_++ } );

    cell_++;
    op_++;
}

template<class F, class I, class C>
void
dionysus::RelativeHomologyZigzag<F,I,C>::
remove_both(Index cell)
{
    Index   abs_cell = absolute_.left.find(cell)->second;
    Index   rel_cell = relative_.left.find(cell)->second;

    zzp_.remove(rel_cell);
    zzp_.remove(abs_cell);

    absolute_.left.erase(cell);
    relative_.left.erase(cell);

    op_map_.insert( { zzp_op_++, op_ } );
    op_map_.insert( { zzp_op_++, op_ } );

    op_++;
}

template<class F, class I, class C>
template<class ChainRange>
typename dionysus::RelativeHomologyZigzag<F,I,C>::Index
dionysus::RelativeHomologyZigzag<F,I,C>::
add(Index cell, const ChainRange& chain)
{
    Index pair  = zzp_.add(relative_chain(cell, chain));
    op_map_.insert( { zzp_op_++, op_++ } );
    relative_.left.insert( { cell, zzp_cell_++ } );

    return decode_pair(pair);
}


template<class F, class I, class C>
typename dionysus::RelativeHomologyZigzag<F,I,C>::Index
dionysus::RelativeHomologyZigzag<F,I,C>::
decode_pair(Index pair)
{
    if (pair == unpaired())
        return pair;

    Index decoded = op_map_.find(pair)->second;
    op_map_.erase(pair);
    return decoded;
}

template<class F, class I, class C>
template<class ChainRange>
typename dionysus::RelativeHomologyZigzag<F,I,C>::IndexChain
dionysus::RelativeHomologyZigzag<F,I,C>::
absolute_chain(const ChainRange& chain) const
{
    IndexChain  res;
    for (const auto& e : chain)
        res.push_back(ChainEntry(e.element(), abs_index(e.index())));
    return res;
}

template<class F, class I, class C>
template<class ChainRange>
typename dionysus::RelativeHomologyZigzag<F,I,C>::IndexChain
dionysus::RelativeHomologyZigzag<F,I,C>::
relative_chain(Index cell, const ChainRange& chain) const
{
    // NB: to compute the signs correctly,
    //     this assumes that the cone vertex w is the last vertex in some total order

    typedef     typename IndexChain::value_type     ChainEntry;

    IndexChain res;
    if (!chain.empty())
    {
        for (const auto& e : chain)
            res.push_back(ChainEntry(e.element(), rel_index(e.index())));

        FieldElement a = field().id();
        if (chain.size() % 2 == 0)      // TODO: double-check
            a = field().neg(a);
        res.push_back(ChainEntry(a, abs_index(cell)));       // add the base space cell
    } else
    {
        res.reserve(2);
        res.push_back(ChainEntry(field().id(), abs_index(cell)));
        res.push_back(ChainEntry(field().neg(field().id()), 0));
    }
    return res;
}


template<class F, class I, class C>
typename dionysus::RelativeHomologyZigzag<F,I,C>::Index
dionysus::RelativeHomologyZigzag<F,I,C>::
remove(Index cell)
{
    Index rel_cell  = rel_index(cell);
    Index pair      = zzp_.remove(rel_cell);
          pair      = decode_pair(pair);

    op_map_.insert( { zzp_op_++, op_++ } );
    relative_.left.erase(cell);

    return pair;
}
