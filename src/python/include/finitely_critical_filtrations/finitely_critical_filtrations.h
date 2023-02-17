#pragma once


namespace Gudhi{

template<typename T>
class Finitely_critical_multi_filtration : public std::vector<T> {
	// Class to prevent doing illegal stuff with the standard library, e.g., compare two vectors
public:
	// using std::vector<T> :: vector; // inherit constructor
	explicit Finitely_critical_multi_filtration(std::vector<T> * v) : std::vector<T>(v) {}; // I'm not sure if it does a copy ?
	explicit Finitely_critical_multi_filtration(int n) : std::vector<T>(n) {};
	explicit Finitely_critical_multi_filtration(int n, T value) : std::vector<T>(n,value) {};

	//TODO : multicritical -> iterator over filtrations
private:

};
template<typename T>
bool operator<(const Finitely_critical_multi_filtration<T>& v1, const Finitely_critical_multi_filtration<T>& v2)
{
	bool isSame = true;
	if (v1.size() != v2.size()){
		std::cerr << "Filtrations are not of the same size ! (" << v1.size() << " " << v2.size() << ").";
		throw;
	}
	int n = v1.size();
	for (int i = 0; i < n; ++i){
		if (v1[i] > v2[i]) return false;
		if (isSame && v1[i] != v2[i]) isSame = false;
	}
	if (isSame) return false;
	return true;
}



template<typename T>
bool operator<=(const Finitely_critical_multi_filtration<T>& v1, const Finitely_critical_multi_filtration<T>& v2)
{
	if (v1.size() != v2.size()){
		std::cerr << "Filtrations are not of the same size ! (" << v1.size() << " " << v2.size() << ").";
		throw;
	}
	int n = v1.size();
	for (int i = 0; i < n; ++i){
		if (v1[i] > v2[i]) return false;
	}
	return true;
}

template<typename T>
Finitely_critical_multi_filtration<T>& operator-=(Finitely_critical_multi_filtration<T> &result, const Finitely_critical_multi_filtration<T> &to_substract){
	std::transform(result.begin(), result.end(), to_substract.begin(),result.begin(), std::minus<T>());
	return result;
}


template<typename T>
Finitely_critical_multi_filtration<T>& operator+=(Finitely_critical_multi_filtration<T> &result, const Finitely_critical_multi_filtration<T> &to_add){
	std::transform(result.begin(), result.end(), to_add.begin(),result.begin(), std::plus<T>());
	return result;
}

} // namespace Gudhi
