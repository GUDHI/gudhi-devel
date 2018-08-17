#include <functional>
#include <vector>

// From https://codereview.stackexchange.com/a/164766
// C++ Concepts: template<ForwardIterator It>


template<typename It>
class SubsetChooser
{
    using Subset = std::vector<It>;
    
    It first;
    It last;
    size_t subset_size;

    Subset state;

public:
    SubsetChooser(It first, It last, size_t subset_size);
    SubsetChooser();
    const Subset& subset() const;
    bool advance();
};

//---
// factory methods

template<typename It>
SubsetChooser<It> make_chooser(It first, It last, size_t subset_size)
{
    return SubsetChooser<It>{first, last, subset_size};
}
template<typename Container>
SubsetChooser<typename Container::const_iterator> make_chooser(const Container& c, size_t subset_size)
{
    using std::begin;
    using std::end;
    return make_chooser(begin(c), end(c), subset_size);
}

// default predicate

template<typename Subset>
bool empty_predicate(const Subset& stack) {
  return true;
}


//---
// private helpers

// Calculate it+n==end, without requiring a BidirectionalIterator
template<typename Iter, typename Distance>
bool is_n_from(Iter it, Distance n, const Iter& end)
{
    std::advance(it, n);
    return it == end;
}

//---
// implementation

template<typename It>
SubsetChooser<It>::SubsetChooser(It first, It last, size_t subset_size)
    : first{first}, last{last},
      subset_size{subset_size},
      state{}
{
    state.reserve(subset_size);
}

template<typename It>
SubsetChooser<It>::SubsetChooser() : subset_size(0) { }


template<typename It>
const typename SubsetChooser<It>::Subset& SubsetChooser<It>::subset() const
{
    return state;
}

template<typename It>
bool SubsetChooser<It>::advance()
{
    do {
        if (state.empty()) {
            state.push_back(first);
        } else {
            if (state.size() < subset_size) {
                state.push_back(state.back());
            }

            // Roll over when the remaining elements wouldn't fill the subset.
            while (is_n_from(++state.back(), subset_size - state.size(), last)) {
                state.pop_back();
                if (state.empty())
                    // we have run out of possibilities
                    return false;
            }
        }
    } while (state.size() < subset_size);
    return true;
}
