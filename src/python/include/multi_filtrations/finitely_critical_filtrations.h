#ifndef FINITELY_CRITICAL_FILTRATIONS_H_
#define FINITELY_CRITICAL_FILTRATIONS_H_

#include <iostream>
#include <algorithm>
#include <limits>

namespace Gudhi::multiparameter::multi_filtrations{

template<typename T=float>
class Finitely_critical_multi_filtration : public std::vector<T> {
	// Class to prevent doing illegal stuff with the standard library, e.g., compare two vectors
public:
	// explicit Finitely_critical_multi_filtration(std::vector<T>& v) : ptr_(&v) {
	// }; // Conversion
	// using std::vector<T>::vector;
	Finitely_critical_multi_filtration() : std::vector<T>() {};
	Finitely_critical_multi_filtration(int n) : std::vector<T>(n, -std::numeric_limits<T>::infinity()) {}; // minus infinity by default
	Finitely_critical_multi_filtration(int n, T value) : std::vector<T>(n,value) {};
	Finitely_critical_multi_filtration(std::initializer_list<T> init) : std::vector<T>(init) {};
	Finitely_critical_multi_filtration(const std::vector<T>& v) : std::vector<T>(v) {};
	Finitely_critical_multi_filtration(typename std::vector<T>::iterator it_begin,typename std::vector<T>::iterator it_end) : std::vector<T>(it_begin, it_end) {};
	Finitely_critical_multi_filtration(typename std::vector<T>::const_iterator it_begin,typename std::vector<T>::const_iterator it_end) : std::vector<T>(it_begin, it_end) {};


	operator std::vector<T>&() const {
		return *this;
	}
	std::vector<T> get_vector() const{
		return static_cast<std::vector<T>>(*this);
	}



	//TODO : multicritical -> iterator over filtrations

	// LESS THAN OPERATORS
	friend bool operator<(const Finitely_critical_multi_filtration& a, const Finitely_critical_multi_filtration& b)
	{
		bool isSame = true;
		// if (a.size() != b.size()){
			
			
		// 	if (a.size()>1 && b.size() >1){
		// 		std::cerr << "Filtrations are not of the same size ! (" << a.size() << " " << b.size() << ").";
		// 		throw;
		// 	}
		// 	// {inf, inf, ...} can be stored as {inf}
		// 	if (b[0] == std::numeric_limits<T>::infinity()) //TODO FIXME, we have to check every coord of b instead of 1
		// 		return a[0] != std::numeric_limits<T>::infinity();
		// 	return false;
		// }
		int n = std::min(a.size(), b.size());
		for (int i = 0; i < n; ++i){
			if (a[i] > b[i]) return false;
			if (isSame && a[i] != b[i]) isSame = false;
		}
		if (isSame) return false;
		return true;
	}
	friend bool operator<=(const Finitely_critical_multi_filtration& a, const Finitely_critical_multi_filtration& b) 
	{
		// if (a.size() != b.size()){
		// 	if (a.size()>1 && b.size() >1){
		// 		std::cerr << "Filtrations are not of the same size ! (" << a.size() << " " << b.size() << ").";
		// 		throw;
		// 	}
		// 	// {inf, inf, ...} can be stored as {inf}
		// 	if (b[0] == std::numeric_limits<T>::infinity()) //TODO FIXME, we have to check every coord of b instead of 1
		// 		return a[0] != std::numeric_limits<T>::infinity();
		// 	return false;
		// }
		int n = std::min(a.size(), b.size());
		for (int i = 0; i < n; ++i){
			if (a[i] > b[i]) return false;
		}
		return true;
	}



	//GREATER THAN OPERATORS
	friend bool operator>(const Finitely_critical_multi_filtration& a, const Finitely_critical_multi_filtration& b)
	{
		return b<a; // C'est honteux.
	}
	friend bool operator>=(const Finitely_critical_multi_filtration& a, const Finitely_critical_multi_filtration& b) 
	{
		return b<=a; // C'est honteux.
	}

	Finitely_critical_multi_filtration& operator=(const Finitely_critical_multi_filtration& a){
		std::vector<T>::operator=(a);
		return *this;
	}

	std::vector<T>& _convert_back(){
		return *this;
	}




	friend Finitely_critical_multi_filtration& operator-=(Finitely_critical_multi_filtration &result, const Finitely_critical_multi_filtration &to_substract){
		std::transform(result.begin(), result.end(), to_substract.begin(),result.begin(), std::minus<T>());
		return result;
	}
	friend Finitely_critical_multi_filtration<T>& operator+=(Finitely_critical_multi_filtration<T> &result, const Finitely_critical_multi_filtration &to_add){
		std::transform(result.begin(), result.end(), to_add.begin(),result.begin(), std::plus<T>());
		return result;
	}

	friend Finitely_critical_multi_filtration& operator-=(Finitely_critical_multi_filtration &result, const T &to_substract){
		// std::transform(result.begin(), result.end(), to_substract.begin(),result.begin(), std::minus<T>());
		for (auto & truc : result){
			truc -= to_substract;
		}
		return result;
	}
	friend Finitely_critical_multi_filtration<T>& operator+=(Finitely_critical_multi_filtration<T> &result, const T &to_add){
		for (auto & truc : result){
			truc += to_add;
		}
		return result;	
	}

	static std::vector<std::vector<T>> to_python(const std::vector<Finitely_critical_multi_filtration<T>>& to_convert){
		return std::vector<std::vector<T>>(to_convert.begin(), to_convert.end());
	}
	

	static std::vector<Finitely_critical_multi_filtration<T>> from_python(const std::vector<std::vector<T>>& to_convert){
		return std::vector<Finitely_critical_multi_filtration<T>>(to_convert.begin(), to_convert.end());;
	}
	void push_to(const Finitely_critical_multi_filtration<T>& x){
		if (this->size() != x.size())
			{std::cerr << "Does only work with 1-critical filtrations ! Sizes " << this->size() << " and " << x.size() << "are different !" << std::endl; return;}
		for (unsigned int i = 0; i < x.size(); i++)
			this->at(i) = this->at(i) > x[i] ? this->at(i) : x[i];
	}
	// Warning, this function  assumes that the comparisons checks have already been made !
	void insert_new(Finitely_critical_multi_filtration to_concatenate){
		this->insert(
			this->end(), std::move_iterator(to_concatenate.begin()), std::move_iterator(to_concatenate.end())
		);
	}

	T linear_projection(const std::vector<T>& x){
		T projection=0;
		unsigned int size = std::min(x.size(), this->size());
		for (auto i =0u; i<size;i++)
			projection += x[i]*this->at(i);
		return projection;
	}

	// easy debug 
    friend std::ostream& operator<<(std::ostream& stream, const Finitely_critical_multi_filtration<T>& truc){
        stream << "[";
        for(unsigned int i = 0; i < truc.size()-1; i++){
            stream << truc[i] << ", ";
        }
        if(!truc.empty()) stream << truc.back();
        stream << "]";
        return stream;
    }




};


} // namespace Gudhi
#endif  // FINITELY_CRITICAL_FILTRATIONS_H_
