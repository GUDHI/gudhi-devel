template<class T>
template<class C2, class Field, class Cmp, class Visitor_>
void
dionysus::Chain<std::list<T>>::
addto(std::list<T>& x, typename Field::Element a, const C2& y, const Field& field, const Cmp& cmp, const Visitor_& visitor)
{
    typedef typename Field::Element                     Element;

    auto cur_x = std::begin(x),
         end_x = std::end(x);
    auto cur_y = std::begin(y),
         end_y = std::end(y);

    while (cur_x != end_x && cur_y != end_y)
    {
        if (cmp(cur_x->index(), cur_y->index()))
        {
            visitor.first(cur_x++);
        } else if (cmp(cur_y->index(), cur_x->index()))
        {
            // multiply and add
            Element ay = field.mul(a, cur_y->element());
            auto nw_x = x.insert(cur_x, *cur_y);
            nw_x->set_element(ay);
            ++cur_y;
            visitor.second(nw_x);
        } else
        {
            Element ay = field.mul(a, cur_y->element());
            Element r  = field.add(cur_x->element(), ay);
            if (field.is_zero(r))
            {
                visitor.equal_drop(cur_x);
                x.erase(cur_x++);
            }
            else
            {
                cur_x->set_element(r);
                visitor.equal_keep(cur_x);
                ++cur_x;
            }
            ++cur_y;
        }
    }

    for (auto it = cur_y; it != end_y; ++it)
    {
        Element ay = field.mul(a, it->element());
        x.push_back(*it);
        x.back().set_element(ay);
        visitor.second(--x.end());
    }
}

template<class T, class TCmp>
template<class C2, class Field, class Cmp, class Visitor_>
void
dionysus::Chain<std::set<T,TCmp>>::
addto(std::set<T,TCmp>& x, typename Field::Element a, const C2& y, const Field& field, const Cmp&, const Visitor_& visitor)
{
    typedef typename Field::Element                     Element;

    auto cur_y = std::begin(y),
         end_y = std::end(y);

    while (cur_y != end_y)
    {
        auto cur_x = x.find(*cur_y);
        if (cur_x == x.end())
        {
            auto nw = x.insert(*cur_y).first;
            Element ay = field.mul(a, nw->element());
            const_cast<T&>(*nw).set_element(ay);
            visitor.second(nw);
        } else
        {
            Element ay = field.mul(a, cur_y->element());
            Element r  = field.add(cur_x->element(), ay);
            if (field.is_zero(r))
            {
                visitor.equal_drop(cur_x);
                x.erase(cur_x);
            }
            else
            {
                const_cast<T&>(*cur_x).set_element(r);
                visitor.equal_keep(cur_x);
            }
        }
        ++cur_y;
    }
}

template<class T, class TCmp>
template<class Field, class Cmp, class Visitor_>
void
dionysus::Chain<std::set<T,TCmp>>::
addto(std::set<T,TCmp>& x, typename Field::Element a, T&& y, const Field& field, const Cmp&, const Visitor_& visitor)
{
    typedef typename Field::Element                     Element;

    auto cur_x = x.find(y);
    if (cur_x == x.end())
    {
        auto nw = x.insert(std::move(y)).first;
        Element ay = field.mul(a, nw->element());
        const_cast<T&>(*nw).set_element(ay);
        visitor.second(nw);
    } else
    {
        Element ay = field.mul(a, y.element());
        Element r  = field.add(cur_x->element(), ay);
        if (field.is_zero(r))
        {
            visitor.equal_drop(cur_x);
            x.erase(cur_x);
        }
        else
        {
            const_cast<T&>(*cur_x).set_element(r);
            visitor.equal_keep(cur_x);
        }
    }
}

template<class C1>
template<class C2, class Field, class Cmp, class Visitor_>
void
dionysus::Chain<C1>::
addto(C1& x, typename Field::Element a, const C2& y, const Field& field, const Cmp& cmp, const Visitor_& visitor)
{
    typedef typename Field::Element                     Element;

    C1 res;

    auto cur_x = std::begin(x),
         end_x = std::end(x);
    auto cur_y = std::begin(y),
         end_y = std::end(y);

    while (cur_x != end_x && cur_y != end_y)
    {
        if (cmp(*cur_x, *cur_y))
        {
            res.emplace_back(std::move(*cur_x));
            visitor.first(--res.end());
            ++cur_x;
        } else if (cmp(*cur_y, *cur_x))
        {
            // multiply and add
            Element ay = field.mul(a, cur_y->element());
            res.emplace_back(ay, cur_y->index());
            visitor.second(--res.end());
            ++cur_y;
        } else
        {
            Element ay = field.mul(a, cur_y->element());
            Element r  = field.add(cur_x->element(), ay);
            if (field.is_zero(r))
                visitor.equal_drop(cur_x);
            else
            {
                res.emplace_back(std::move(*cur_x));
                res.back().set_element(r);
                visitor.equal_keep(--res.end());
            }
            ++cur_x;
            ++cur_y;
        }
    }

    while (cur_y != end_y)
    {
        Element ay = field.mul(a, cur_y->element());
        res.emplace_back(ay, cur_y->index());
        visitor.second(--res.end());
        ++cur_y;
    }

    while (cur_x != end_x)
    {
        res.emplace_back(std::move(*cur_x));
        visitor.first(--res.end());
        ++cur_x;
    }

    x.swap(res);
}
