template<class F, typename I, class C, template<class Self> class... V>
void
dionysus::ReducedMatrix<F,I,C,V...>::
resize(size_t s)
{
    reduced_.resize(s);
    pairs_.resize(s, unpaired());
    skip_.resize(s, false);
}

template<class F, typename I, class C, template<class Self> class... V>
typename dionysus::ReducedMatrix<F,I,C,V...>::Index
dionysus::ReducedMatrix<F,I,C,V...>::
add(Chain&& chain)
{
    // TODO: skip the computation entirely if we already know this is positive (in case of the clearing optimization)
    Index i = pairs_.size();
    pairs_.emplace_back(unpaired());
    reduced_.emplace_back();
    skip_.push_back(false);

    set(i, std::move(chain));

    return reduce(i);
}

template<class F, typename I, class C, template<class Self> class... V>
void
dionysus::ReducedMatrix<F,I,C,V...>::
add_skip()
{
    pairs_.emplace_back(unpaired());
    reduced_.emplace_back();
    skip_.push_back(true);
}

template<class F, typename I, class C, template<class Self> class... V>
void
dionysus::ReducedMatrix<F,I,C,V...>::
set(Index i, Chain&& c)
{
    sort(c);
    visitors_chain_initialized(c);
    reduced_[i] = std::move(c);
}

template<class F, typename I, class C, template<class Self> class... V>
typename dionysus::ReducedMatrix<F,I,C,V...>::Index
dionysus::ReducedMatrix<F,I,C,V...>::
reduce(Index i)
{
    Chain& c    = column(i);
    Index  pair = reduce(c);

    if (pair != unpaired())
        pairs_[pair] = i;

    pairs_[i]   = pair;
    visitors_reduction_finished<>();

    return pair;
}

template<class F, typename I, class C, template<class Self> class... V>
template<class ChainsLookup,
         class LowLookup>
typename dionysus::ReducedMatrix<F,I,C,V...>::Index
dionysus::ReducedMatrix<F,I,C,V...>::
reduce(      Chain&                c,
       const ChainsLookup&         chains,
       const LowLookup&            lows)
{
    auto entry_cmp = [this](const Entry& e1, const Entry& e2) { return this->cmp_(e1.index(), e2.index()); };
    return Reduction<Index>::reduce(c, chains, lows, field_,
                                    [this](FieldElement m, Index cl)
                                    { this->visitors_addto<>(m, cl); },
                                    entry_cmp);
}
