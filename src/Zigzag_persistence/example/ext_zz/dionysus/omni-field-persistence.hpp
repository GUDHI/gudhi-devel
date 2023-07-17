template<typename Index_, class Comparison_, class Q_, class Zp_>
void
dionysus::OmniFieldPersistence<Index_, Comparison_, Q_, Zp_>::
add(QChain&& chain)
{
    sort(chain);

    q_chains_.emplace_back(std::move(chain));
    q_pairs_.emplace_back(unpaired());
    Index i = q_chains_.size() - 1;

    QChain& c = q_chains_.back();

    auto reduce = [this,&c,i](BaseElement p)
    {
        auto zp_chain = convert(c, zp(p));

        this->reduce(zp_chain, p);

        if (!zp_chain.empty())
        {
            auto l = zp_chain.back().index();
            zp_lows_[l].emplace(p,i);
            set_pair(l,i,p);
        }

        zp_chains_[i].emplace(p, std::move(zp_chain));        // empty chain is still a valid indicator that we don't need to bother with this field
    };

    // reduce
    auto entry_cmp = [this](const QEntry& e1, const QEntry& e2) { return this->cmp_(e1.index(), e2.index()); };
    while (!c.empty())
    {
        auto& low = c.back();

        auto e = low.element();
        auto l = low.index();
        assert(!q_.is_zero(e));
        if (e != q_.id())
        {
            auto factors = factor(q_.numerator(e));
            for (auto p : factors)
            {
                if (!special(i, p))        // there is already a dedicated column over p
                    reduce(p);
            }
        }

        auto it_zp = zp_lows_.find(l);
        if (it_zp != zp_lows_.end())
            for (auto& x : it_zp->second)
            {
                auto p = x.first;
                if (!special(i,p))
                    reduce(p);
            }

        auto it_q = q_lows_.find(l);
        if (it_q != q_lows_.end())
        {
            Index j = it_q->second;

            // add the primes from j to i
            auto it_zp = zp_chains_.find(j);
            if (it_zp != zp_chains_.end())
                for (auto& x : it_zp->second)
                {
                    auto p = x.first;
                    if (!special(i,p))
                        reduce(p);
                }

            // reduce over Q
            auto  j_chain  = q_chains_[j];
            auto  j_e      = j_chain.back().element();

            auto  m = q_.neg(q_.div(e,j_e));
            Chain<QChain>::addto(c, m, j_chain, q_, entry_cmp);
            assert(c.empty() || !q_.is_zero(c.back().element()));
        } else
        {
            q_lows_.emplace(l,i);
            set_pair(l,i);
            break;
        }
    }
}

template<typename Index_, class Comparison_, class Q_, class Zp_>
void
dionysus::OmniFieldPersistence<Index_,Comparison_,Q_,Zp_>::
reduce(ZpChain& zp_chain, BaseElement p)
{
    auto& field = zp(p);

    auto entry_cmp = [this](const ZpEntry& e1, const ZpEntry& e2) { return this->cmp_(e1.index(), e2.index()); };

    while (!zp_chain.empty())
    {
        auto& low = zp_chain.back();
        auto  j   = low.index();

        auto it = zp_lows_.find(j);
        if (it != zp_lows_.end())
        {
            auto it2 = it->second.find(p);
            if (it2 != it->second.end())
            {
                const ZpChain& co = zp_chains_[it2->second][p];

                auto  m = field.neg(field.div(low.element(), co.back().element()));
                assert(m < p);
                Chain<ZpChain>::addto(zp_chain, m, co, field, entry_cmp);
                continue;
            }
        }

        auto qit = q_lows_.find(j);
        if (qit == q_lows_.end() || special(qit->second, p))       // no valid pivot over Q
            return;

        // TODO: this could be optimized (add and convert on the fly)
        auto& q_chain = q_chains_[qit->second];
        assert(q_chain.empty() || !q_.is_zero(q_chain.back().element()));

        auto  co      = convert(q_chain, field);
        auto  m       = field.neg(field.div(low.element(), co.back().element()));
        Chain<ZpChain>::addto(zp_chain, m, co, field, entry_cmp);

        assert(!zp_chain.empty() || zp_chain.back().index() != j);
    }
}

template<typename Index_, class Comparison_, class Q_, class Zp_>
typename dionysus::OmniFieldPersistence<Index_,Comparison_,Q_,Zp_>::ZpChain
dionysus::OmniFieldPersistence<Index_,Comparison_,Q_,Zp_>::
convert(const QChain& c, const Zp& field) const
{
    ZpChain result;
    result.reserve(c.size());
    auto p = field.prime();
    for (auto& x : c)
    {
        auto num = q_.numerator(x.element()) % p;
        if (num != 0)
        {
            while (num < 0) num += p;
            auto denom = q_.denominator(x.element()) % p;
            while (denom < 0) denom += p;
            assert(denom % p != 0);
            result.emplace_back(field.div(num, denom), x.index());
        }
    }
    return result;
}


template<typename Index_, class Comparison_, class Q_, class Zp_>
typename dionysus::OmniFieldPersistence<Index_,Comparison_,Q_,Zp_>::Factors
dionysus::OmniFieldPersistence<Index_, Comparison_,Q_,Zp_>::
factor(BaseElement x)
{
    if (x < 0)
        x = -x;
    Factors result;

    if (Q::is_prime(x))
    {
        result.push_back(x);
        return result;
    }

    BaseElement p { 2 };
    while (p*p <= x)
    {
        if (x % p == 0)
        {
            result.push_back(p);
            do { x /= p; } while (x % p == 0);
            if (Q::is_prime(x))
            {
                result.push_back(x);
                break;
            }
        }
        ++p;
    }
    if (x > 1)
        result.push_back(x);

    return result;
}

template<typename Index_, class Comparison_, class Q_, class Zp_>
typename dionysus::OmniFieldPersistence<Index_,Comparison_,Q_,Zp_>::Index
dionysus::OmniFieldPersistence<Index_, Comparison_, Q_, Zp_>::
pair(Index i, BaseElement p) const
{
    if (p == 1)
        return q_pairs_[i];
    else
    {
        auto it = zp_pairs_.find(p);
        if (it == zp_pairs_.end())
            return q_pairs_[i];
        else
        {
            auto pit = it->second.find(i);
            if (pit == it->second.end())
                return q_pairs_[i];
            else
                return pit->second;
        }
    }
}

template<typename Index_, class Comparison_, class Q_, class Zp_>
void
dionysus::OmniFieldPersistence<Index_, Comparison_, Q_, Zp_>::
set_pair(Index i, Index j, BaseElement p)
{
    auto& pairs = zp_pairs_[p];
    pairs[i] = j;
    pairs[j] = i;
}

template<typename Index_, class Comparison_, class Q_, class Zp_>
void
dionysus::OmniFieldPersistence<Index_, Comparison_, Q_, Zp_>::
set_pair(Index i, Index j)
{
    q_pairs_[i] = j;
    q_pairs_[j] = i;

    auto it = zp_chains_.find(j);
    if (it == zp_chains_.end())
        return;

    auto& chains = it->second;
    for (auto& x : chains)
    {
        auto  p     = x.first;
        auto& chain = x.second;
        if (chain.empty())
        {
            zp_pairs_[p][j] = unpaired();
            zp_pairs_[p][i] = unpaired();
        }
    }
}
