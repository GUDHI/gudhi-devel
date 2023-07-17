#ifndef DIONYSUS_PAIR_RECORDER_H
#define DIONYSUS_PAIR_RECORDER_H

namespace dionysus
{

template<class Persistence_>
struct PairRecorder: public Persistence_
{
    typedef             Persistence_                    Persistence;
    typedef             typename Persistence::Index     Index;


    using Persistence::Persistence;

    template<class ChainRange>
    Index               add(const ChainRange& chain)
    {
        Index p = Persistence::add(chain);
        pairs_.push_back(p);
        if (p != unpaired())
            pairs_[p] = pairs_.size() - 1;

        return p;
    }

    Index               pair(Index i) const             { return pairs_[i]; }

    void                resize(size_t s)                { Persistence::resize(s); pairs_.resize(s, unpaired()); }
    size_t              size() const                    { return pairs_.size(); }
    static const Index  unpaired()                      { return Reduction<Index>::unpaired; }

    std::vector<Index>  pairs_;
};

template<class Persistence_>
struct PairChainRecorder: public PairRecorder<Persistence_>
{
    using Persistence = Persistence_;
    using Parent      = PairRecorder<Persistence>;
    using Index       = typename Persistence_::Index;
    using Chain       = typename Persistence_::Chain;

    using Parent::Parent;

    template<class ChainRange>
    Index               add(const ChainRange& chain)
    {
        auto  p_chain = Persistence::add(chain, keep_cocycles);
        Index p       = std::get<0>(p_chain);

        pairs_.push_back(p);
        chains_.emplace_back();

        if (p != unpaired())
        {
            pairs_[p] = pairs_.size() - 1;
            chains_[p] = std::move(std::get<1>(p_chain));
        }

        return p;
    }

    using Parent::unpaired;

    Index               pair(Index i) const             { return pairs_[i]; }
    const Chain&        chain(Index i) const            { return chains_[i]; }      // chain that dies at i
    void                resize(size_t s)                { Parent::resize(s); chains_.resize(s); }

    std::vector<Chain>  chains_;
    using Parent::pairs_;

    bool                keep_cocycles = true;
};

}

#endif
