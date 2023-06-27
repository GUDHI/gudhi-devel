#include <map>

template<class F, class I, class C>
template<class ChainRange>
typename dionysus::ZigzagPersistence<F,I,C>::Index
dionysus::ZigzagPersistence<F,I,C>::
add_impl(const ChainRange& chain_)
{
    // std::cout << "add(" << cell_indices << ")" << std::endl;
    Index op = operations++;

    IndexChain cycles;      // chain_ -> Z*cycles
    Column     z_remainder = Z.reduce(chain_, cycles);
	// std::cout << "cycle: ";
	// for (auto& v : cycles){
	// 	std::cout << v.index() << " ";
	// }
	// std::cout << "\n";
    assert(z_remainder.empty());

    IndexChain  boundaries;
    DequeColumn b_remainder = B.reduce(cycles, boundaries);

    // add up columns of C indexed by boundaries
    typedef     typename Column::value_type             Entry;
    auto        row_cmp = [this](const Entry& e1, const Entry& e2)
                          { return this->cmp()(row(e1), row(e2)); };
    Column      chain;
    for (auto& x : boundaries)
        Chain<Column>::addto(chain, x.element(), C.col(x.index()), field(), row_cmp);
    chain.push_back(Entry(field().neg(field().id()), IndexPair(cell_indices++,0)));

    if (b_remainder.empty())        // birth
    {
        // std::cout << "  birth" << std::endl;
        Index z_col = z_indicies_last++;
        Z.set(z_col, std::move(chain));
        birth_index[z_col] = op;
        return unpaired();
    }
    else                            // death
    {
        // std::cout << "  death" << std::endl;
        Index b_col = b_indices++;
        Index pair  = row(b_remainder.back());
        B.set(b_col, std::move(b_remainder));
        C.set(b_col, std::move(chain));
        return birth_index[pair];
    }
}

template<class F, class I, class C>
typename dionysus::ZigzagPersistence<F,I,C>::Index
dionysus::ZigzagPersistence<F,I,C>::
remove_impl(Index cell)
{
    //std::cout << "remove(" << cell << ")" << std::endl;

    Index   op    = operations++;

    typedef     typename Column::value_type             Entry;
    auto        row_cmp = [this](const Entry& e1, const Entry& e2)
                          { return this->cmp()(row(e1), row(e2)); };
    typedef     typename DequeColumn::value_type        DequeEntry;
    auto        b_row_cmp = [this](const DequeEntry& e1, const DequeEntry& e2)
                            { return this->cmp()(row(e1), row(e2)); };

    IndexChain  z_row;
    for (auto& x : Z.row(cell))
        z_row.emplace_back(x.element(), col(x));

    if (z_row.empty())              // birth
    {
        //std::cout << "  birth" << std::endl;
        Row&    c_row   = C.row(cell);
        // c_row.front() may not be the first column in order, but that doesn't really matter, does it? (TODO)
        auto&   c_front = c_row.front();

        Index   j     = col(c_front);
        Index   l     = row(B.col(j).back());

        //std::cout << j << ' ' << l << std::endl;

        // cycle = ZB[j] = DC[j]
        Column  cycle;
        for (auto& x : B.col(j))
            Chain<Column>::addto(cycle, x.element(), Z.col(row(x)), field(), row_cmp);

        //std::cout << "Cycle:" << std::endl;
        //for (auto& x : cycle)
        //    std::cout << x.element() << ' ' << row(x) << std::endl;

        // 1: prepend the cycle
        Index   znew        = z_indicies_first--;
        Index   oth         = Z.set(znew, std::move(cycle));        // oth records our collision (used in step 6)
        birth_index[znew]   = op;

        //std::cout << "znew oth: " << znew << ' ' << oth << std::endl;
        //std::cout << "oth column:" << std::endl;
        //for (auto& x : Z.col(oth))
        //    std::cout << x.element() << ' ' << row(x) << std::endl;

        // 2: prepend the row to B
        FieldElement    m     = field().neg(field().inv(c_front.element()));        // m = -1/c
        const DequeRow& b_row = B.prepend_row(znew, m, c_row);
        //std::cout << "Prepended row with multiplier: " << m << " (" << b_row.size() << ")" << std::endl;

        // 3: subtract C[j] from every C[k]
        const Column&   Cj    = C.col(j);

        // use the copy of c_row in B, since c_row will be modified in the following loop
        for (auto it = std::next(b_row.begin()); it != b_row.end(); ++it)
        {
            Index c = col(*it);
            assert(c != j);
            //std::cout << "adding to " << c << " in C" << std::endl;
            Chain<Column>::addto(C.col(c), it->element(), Cj, field(), row_cmp);    // using it->element() since b_row = m*c_row
            C.fix(c);                                                               // old elements got removed via auto_unlink_hook
            // we don't need lows in C, so not updating them
        }
        //std::cout << "Done with step 3" << std::endl;

        // 4: subtract B[j] from every B[k] that has l
        //    (we don't need to update C because ZB[j] = 0 after step 2)
        DequeColumn&    Bj      = B.col(j);
        FieldElement    bm      = field().neg(field().inv(Bj.back().element()));    // bm = -1/B[l,j]
        IndexChain      Bl_row; // make a copy of Bl_row, since it will be changing
        for (auto& x : B.row(l))
        {
            if (col(x) == j)
                continue;
            Bl_row.emplace_back(x.element(), col(x));
        }
        for (auto& x : Bl_row)
        {
            Index c = x.index();
            assert(c != j);
            Chain<DequeColumn>::addto(B.col(c), field().mul(bm, x.element()), Bj, field(), b_row_cmp);
            B.fix(c);                                                               // old elements got removed via auto_unlink_hook
            // l cannot be the low in c, so no need to update lows
        }
        //std::cout << "Done with step 4" << std::endl;

        // 5: drop row l and column j from B; drop column l from Z; drop column j from C
        B.drop_col(j);
        assert(B.row(l).empty());
        B.drop_row(l);
        Index Zl_low = row(Z.col(l).back());
        Z.drop_col(l);
        birth_index.erase(l);
        C.drop_col(j);
        assert(Z.row(cell).empty());
        assert(C.row(cell).empty());
        C.drop_row(cell);
        Z.drop_row(cell);
        //std::cout << "Done with step 5" << std::endl;
        if (oth == l)       // we just dropped our collision in Z
            oth = znew;
        else
            Z.drop_low(Zl_low);

        // 6: reduce Z
        std::unordered_map<Index, DequeColumn>  b_changes;  // the columns to add in B to apply row changes
        Index cur = znew;
        while (oth != cur)
        {
            Column& cur_col = Z.col(cur);
            Column& oth_col = Z.col(oth);
            assert(row(cur_col.back()) == row(oth_col.back()));
            //std::cout << "--- " << cur << " (" << cur_col.size() << ") " << oth << " (" << oth_col.size() << ")" << std::endl;
            FieldElement m1 = cur_col.back().element();
            FieldElement m2 = oth_col.back().element();
            FieldElement m2_div_m1 = field().div(m2, m1);
            Chain<Column>::addto(oth_col, field().neg(m2_div_m1), cur_col, field(), row_cmp);
            Z.fix(oth, oth_col);

            // record the changes we need to make in B;
            //   because there is only one collision in the matrix during the reduction,
            //   once we use a row as the source, we never revisit it. This means once the row is updated in B,
            //   we never touch it again, so below record is fine.
            for (auto& x : this->B.row(oth))
                b_changes[col(x)].emplace_back(field().mul(x.element(), m2_div_m1), cur, col(x));

            cur = oth;
            Index low = row(oth_col.back());
            if (Z.is_low(low))
                oth = Z.low(low);
            //std::cout << "--- -- new low: " << low << ' ' << cur << ' ' << oth << std::endl;

            if (cmp()(oth, cur))
                std::swap(oth, cur);
            else
                Z.update_low(cur);
        }

        // apply changes in B (the complexity here could get ugly)
        for (auto& bx : b_changes)
        {
            std::sort(bx.second.begin(), bx.second.end(), b_row_cmp);
            Chain<DequeColumn>::addto(B.col(bx.first), field().id(), bx.second, field(), b_row_cmp);
            B.fix(bx.first);
            // no need to update low (additions from bottom up)
        }
        //std::cout << "Done with step 6" << std::endl;

        return unpaired();
    }
    else                            // death
    {
        //std::cout << "  death" << std::endl;

        auto index_chain_cmp = [this](const typename IndexChain::value_type& e1, const typename IndexChain::value_type& e2)
                               { return this->cmp()(e1.index(), e2.index()); };

        // 1: change basis to clear z_row
        std::sort(z_row.begin(), z_row.end(), index_chain_cmp);     // this adds a log factor, but it makes life easier
        Index        j = z_row.front().index();
        FieldElement e = z_row.front().element();

        if (z_row.size() > 1)
        {
            // figure out the columns we use for reduction
            typedef     typename IndexChain::const_iterator     RowIterator;
            std::vector<RowIterator>                            reducers;
            reducers.push_back(z_row.begin());
            for (RowIterator it = std::next(z_row.begin()); it != z_row.end(); ++it)
            {
                Index c = it->index();

                assert(Z.col_exists(c));
                assert(Z.col_exists(reducers.back()->index()));
                if (cmp()(row(Z.col(c).back()),
                          row(Z.col(reducers.back()->index()).back())))
                    reducers.push_back(it);
            }
            reducers.push_back(z_row.end());
            //std::cout << "reducers.size(): " << reducers.size() << std::endl;
            //std::cout << "z_row.size():    " << z_row.size() << std::endl;


            std::map<Index, IndexChain>  b_changes;  // the rows to add to B
            auto add_in_z = [this,&b_changes,&row_cmp,&index_chain_cmp](Index to, Index from, FieldElement m, FieldElement e)
                            {
                                //std::cout << "  add_in_z: " << from << ' ' << to << std::endl;

                                FieldElement    mult = this->field().mul(m, e);
                                assert(Z.col_exists(to));
                                assert(Z.col_exists(from));
                                Chain<Column>::addto(Z.col(to), mult, Z.col(from), this->field(), row_cmp);
                                assert(!Z.col(to).empty());
                                this->Z.fix(to);       // NB: rows will be linked in the back, so the iterators are Ok
                                this->Z.update_low(to);

                                // subtract B.row(to) from B.row(from)
                                IndexChain Bto_row;
                                for (auto& x : this->B.row(to))
                                    Bto_row.emplace_back(x.element(), col(x));
                                std::sort(Bto_row.begin(), Bto_row.end(), index_chain_cmp);

#if 0
                                for (auto& x : this->B.row(to))
                                    std::cout << x.element() << ' ' << row(x) << ' ' << col(x) << std::endl;

                                std::cout << "---\n";

                                for (auto& x : this->B.row(from))
                                    std::cout << x.element() << ' ' << row(x) << ' ' << col(x) << std::endl;
#endif

                                Chain<IndexChain>::addto(b_changes[from], this->field().neg(mult), Bto_row, this->field(), index_chain_cmp);

                                // if there is b_changes[to] add it, too
                                auto it = b_changes.find(to);
                                if (it != b_changes.end())
                                    Chain<IndexChain>::addto(b_changes[from], this->field().neg(mult), it->second, this->field(), index_chain_cmp);
                            };
            Index last_low = row(Z.col(reducers[reducers.size() - 2]->index()).back());
            for (int i = reducers.size() - 2; i >= 0; --i)
            {
                auto rit = reducers[i];
                FieldElement m = field().neg(field().inv(rit->element()));

                for (auto it  = std::next(rit); it != reducers[i+1]; ++it)
                    add_in_z(it->index(), rit->index(), m, it->element());

                if (static_cast<unsigned int>(i + 1) != reducers.size() - 1)
                {
                    auto it = reducers[i+1];
                    add_in_z(it->index(), rit->index(), m, it->element());
                }
            }
            if (reducers.size() > 2)
                Z.drop_low(last_low);

            // apply changes in b (the complexity here could get ugly)
            // Specifically, transpose b_changes and add it in
            std::unordered_map<Index, DequeColumn>  b_changes_transposed;
            for (auto& b_row : b_changes)
                for (auto& bx : b_row.second)
                    b_changes_transposed[bx.index()].emplace_back(bx.element(), b_row.first, bx.index());

            for (auto& b_col : b_changes_transposed)
            {
#if 0
                std::cout << "Adding:" << std::endl;
                for (auto& x : b_col.second)
                    std::cout << x.element() << ' ' << row(x) << ' ' << col(x) << std::endl;
#endif
                Chain<DequeColumn>::addto(B.col(b_col.first), field().id(), b_col.second, field(), b_row_cmp);
                assert(!B.col(b_col.first).empty());
                B.fix(b_col.first);
                // no need to update low (additions from bottom up)
            }
        }   // z_row.size() > 1

        // 2: subtract cycle from every chain in C
        const Column& Zj = Z.col(j);
        //std::cout << "Zj:" << std::endl;
        //for (auto& x : Zj)
        //    std::cout << x.element() << " * " << row(x) << std::endl;

        IndexChain Ccols;       // save the columns in C, we'll be modifying C.row(cell)
        for (auto& x : C.row(cell))
            Ccols.emplace_back(x.element(), col(x));

        for (auto& x : Ccols)
        {
            Index           c = x.index();
            FieldElement    m = field().neg(field().div(x.element(), e));      // m = -C[k][cell]/Z[j][cell]
            //std::cout << "Adding to C: " << c << std::endl;
            Chain<Column>::addto(C.col(c), m, Zj, field(), row_cmp);
            C.fix(c);
            // we don't care about lows in C, so don't update them
        }

        // 3: drop
        assert(Z.row(cell).size() == 1);
        Z.drop_col(j);
        assert(Z.row(cell).empty());
        assert(C.row(cell).empty());
        Z.drop_row(cell);
        C.drop_row(cell);
        assert(B.row(j).empty());
        B.drop_row(j);

        Index birth = birth_index[j];
        birth_index.erase(j);

        return birth;
    }
}


/* debug routines */
template<class F, class I, class C>
void
dionysus::ZigzagPersistence<F,I,C>::
check_b_cols() const
{
    // check that entries in B refer to existing Z columns
    bool stop = false;
    for (auto& b : B.columns())
        for (auto& x : b.second)
            if (!Z.col_exists(row(x)))
            {
                std::cout << "B refers to a non-existent column in Z: " << row(x) << std::endl;
                stop = true;
            }
    if (stop)
        assert(0);
}

template<class F, class I, class C>
template<class SimplexToIndex, class IndexToSimplex>
void
dionysus::ZigzagPersistence<F,I,C>::
check_cycles(const SimplexToIndex& s2i, const IndexToSimplex& i2s) const
{
    typedef     typename Column::value_type             Entry;
    auto        row_cmp = [this](const Entry& e1, const Entry& e2)
                          { return this->cmp()(row(e1), row(e2)); };

    for (auto& z : Z.columns())
    {
        Column res;
        for (auto& x : z.second)
        {
            Column bdry = boundary(row(x), s2i, i2s);
            Chain<Column>::addto(res, x.element(), bdry, field(), row_cmp);
        }
        assert(res.empty());
    }
}

template<class F, class I, class C>
template<class SimplexToIndex, class IndexToSimplex>
void
dionysus::ZigzagPersistence<F,I,C>::
check_boundaries(const SimplexToIndex& s2i, const IndexToSimplex& i2s) const
{
    check_cycles(s2i, i2s);

    for (auto& x : B.columns())
        if (!C.col_exists(x.first))
        {
            std::cout << x.first << " in B, but not in C" << std::endl;
            assert(0);
        }

    for (auto& x : C.columns())
        if (!B.col_exists(x.first))
        {
            std::cout << x.first << " in B, but not in C" << std::endl;
            assert(0);
        }

    for (auto& x : B.columns())
    {
        auto zb = zb_dot(x.first);
        auto dc = dc_dot(x.first, s2i, i2s);

        auto it_zb = zb.begin(),
             it_dc = dc.begin();
        for (; it_zb != zb.end(); ++it_zb, ++it_dc)
        {
            if (it_zb->element() != it_dc->element() || row(*it_zb) != row(*it_dc))
            {
                std::cout << "Boundary mismatch: " << x.first << std::endl;
                std::cout << "===" << std::endl;
                for (auto& x : zb)
                    std::cout << "   " << x.element() << ' ' << row(x) << std::endl;
                for (auto& y : B.col(x.first))
                {
                    std::cout << "   " << y.element() << " * " << row(y) << std::endl;
                    for (auto& z : Z.col(row(y)))
                        std::cout << "      " << z.element() << ' ' << row(z) << std::endl;
                    std::cout << "   ---" << std::endl;
                }
                std::cout << "===" << std::endl;
                for (auto& x : dc)
                    std::cout << "   " << x.element() << ' ' << row(x) << std::endl;
                for (auto& y : C.col(x.first))
                {
                    std::cout << "   " << y.element() << " * " << row(y) << std::endl;
                    for (auto& z : boundary(row(y), s2i, i2s))
                        std::cout << "     " << z.element() << ' ' << row(z) << std::endl;
                    std::cout << "   ---" << std::endl;
                }
                assert(0);
            }
        }
        if (it_zb != zb.end() || it_dc != dc.end())
        {
            std::cout << "zb.end() doesn't match dc.end()" << std::endl;
            assert(0);
        }
    }
}

template<class F, class I, class C>
typename dionysus::ZigzagPersistence<F,I,C>::Column
dionysus::ZigzagPersistence<F,I,C>::
zb_dot(Index c) const
{
    typedef     typename Column::value_type             Entry;
    auto        row_cmp = [this](const Entry& e1, const Entry& e2)
                          { return this->cmp()(row(e1), row(e2)); };
    Column res;
    for (auto& x : B.col(c))
        Chain<Column>::addto(res, x.element(), Z.col(row(x)), field(), row_cmp);

    return res;
}

template<class F, class I, class C>
template<class SimplexToIndex, class IndexToSimplex>
typename dionysus::ZigzagPersistence<F,I,C>::Column
dionysus::ZigzagPersistence<F,I,C>::
dc_dot(Index c, const SimplexToIndex& s2i, const IndexToSimplex& i2s) const
{
    typedef     typename Column::value_type             Entry;
    auto        row_cmp = [this](const Entry& e1, const Entry& e2)
                          { return this->cmp()(row(e1), row(e2)); };
    Column res;
    for (auto& x : C.col(c))
    {
        Column bdry = boundary(row(x), s2i, i2s);
        Chain<Column>::addto(res, x.element(), bdry, field(), row_cmp);
    }
    return res;
}

template<class F, class I, class C>
template<class SimplexToIndex, class IndexToSimplex>
typename dionysus::ZigzagPersistence<F,I,C>::Column
dionysus::ZigzagPersistence<F,I,C>::
boundary(Index i, const SimplexToIndex& s2i, const IndexToSimplex& i2s) const
{
    typedef     typename Column::value_type             Entry;
    auto        row_cmp = [this](const Entry& e1, const Entry& e2)
                          { return this->cmp()(row(e1), row(e2)); };
    Column bdry;
    auto s = i2s(i);
    for (auto y : s.boundary(field()))
        bdry.emplace_back(y.element(), s2i(y.index()), 0);
    std::sort(bdry.begin(), bdry.end(), row_cmp);
    return bdry;
}

template<class F, class I, class C>
void
dionysus::ZigzagPersistence<F,I,C>::
check_sorted() const
{
    typedef     typename Column::value_type             Entry;
    auto        row_cmp = [this](const Entry& e1, const Entry& e2)
                          { return this->cmp()(row(e1), row(e2)); };
    typedef     typename DequeColumn::value_type        DequeEntry;
    auto        b_row_cmp = [this](const DequeEntry& e1, const DequeEntry& e2)
                            { return this->cmp()(row(e1), row(e2)); };

    for (auto& x : Z.columns())
        if (!std::is_sorted(x.second.begin(), x.second.end(), row_cmp))
        {
            std::cout << "Z column not sorted: " << x.first << std::endl;
            assert(0);
        }
    for (auto& x : C.columns())
        if (!std::is_sorted(x.second.begin(), x.second.end(), row_cmp))
        {
            std::cout << "C column not sorted: " << x.first << std::endl;
            assert(0);
        }
    for (auto& x : B.columns())
        if (!std::is_sorted(x.second.begin(), x.second.end(), b_row_cmp))
        {
            std::cout << "B column not sorted: " << x.first << std::endl;
            assert(0);
        }
}

