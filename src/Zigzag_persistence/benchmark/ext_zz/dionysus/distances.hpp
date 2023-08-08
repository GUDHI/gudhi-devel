template<class Distances_>
dionysus::ExplicitDistances<Distances_>::
ExplicitDistances(const Distances& distances): 
    size_(distances.size()), distances_((distances.size() * (distances.size() + 1))/2)
{
    IndexType i = 0;
    for (typename Distances::IndexType a = distances.begin(); a != distances.end(); ++a)
        for (typename Distances::IndexType b = a; b != distances.end(); ++b)
        {
            distances_[i++] = distances(a,b);
        }
}

template<class Distances_>
typename dionysus::ExplicitDistances<Distances_>::DistanceType
dionysus::ExplicitDistances<Distances_>::
operator()(IndexType a, IndexType  b) const
{
    if (a > b) std::swap(a,b);
    return distances_[a*size_ - ((a*(a-1))/2) + (b-a)];
}

template<class Distances_>
typename dionysus::ExplicitDistances<Distances_>::DistanceType&
dionysus::ExplicitDistances<Distances_>::
operator()(IndexType a, IndexType  b)
{
    if (a > b) std::swap(a,b);
    return distances_[a*size_ - ((a*(a-1))/2) + (b-a)];
}
